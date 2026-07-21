#include "cppdefs.h"
      MODULE extract_slice_mod
#ifdef SOLVE3D
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine extracts a slice of a 3D ROMS at requested constant    !
!  depths, like z=-100m, using vertical linear interpolation, It is    !
!  done for output puposes.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     tile       Domain partition (integer)                            !
!     LBi        I-dimension Lower bound (integer)                     !
!     UBi        I-dimension Upper bound (integer)                     !
!     LBj        J-dimension Lower bound (integer)                     !
!     UBj        J-dimension Upper bound (integer)                     !
!     LBk        K-dimension Lower bound (integer)                     !
!     UBk        K-dimension Upper bound (integer)                     !
!     Adata      3D field (real)                                       !
!     Adepth     3D field terrain following depths (real; negative, m) !
!     Amask      Land/sea masking (real)                               !
!     z_slice    Depths of requested slices (real; negative, m)        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aslice      Extracted field constant depth slices (real)         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
      implicit none
!
      PUBLIC  :: extract_slice
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE extract_slice (ng, model, tile, gtype,                 &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          Adata, Adepth,                          &
# ifdef MASKING
     &                          Amask,                                  &
# endif
     &                          z_slice, Aslice)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,  intent(in)  :: ng, model, tile, gtype
      integer,  intent(in)  :: LBi, UBi, LBj, UBj, LBk, UBk
      real(r8), intent(in)  :: z_slice(:)
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)  :: Adata (LBi:,LBj:,LBk:)
      real(r8), intent(in)  :: Adepth(LBi:,LBj:,LBk:)
#  ifdef MASKING
      real(r8), intent(in)  :: Amask (LBi:,LBj:)
#  endif
      real(r8), intent(out) :: Aslice(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in)  :: Adata (LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(in)  :: Adepth(LBi:UBi,LBj:UBj,LBk:UBk)
#  ifdef MASKING
      real(r8), intent(in)  :: Amask (LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: Aslice(LBi:UBi,LBj:UBj,:)
# endif
!
!  Local variable declarations.
!
      integer :: i, iz, j, k, k1, k2
      integer :: Imin, Imax, Jmin, Jmax
      integer :: mslice
!
      real(r8) :: my_mask
      real(r8) :: dz, r1, r2, zbot, ztop
!
!-----------------------------------------------------------------------
!  Extract field slice(s) from 3D field using linear interpolation.
!-----------------------------------------------------------------------
!
!  Select tile data bounds acordin to staggered C-grid type.
!
      SELECT CASE (gtype)
        CASE (p3dvar)
          Imin=BOUNDS(ng)%Istr (tile)
          Imax=BOUNDS(ng)%Iend (tile)
          Jmin=BOUNDS(ng)%Jstr (tile)
          Jmax=BOUNDS(ng)%Jend (tile)
        CASE (r3dvar, w3dvar)
          Imin=BOUNDS(ng)%IstrR(tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%JstrR(tile)
          Jmax=BOUNDS(ng)%JendR(tile)
        CASE (u3dvar)
          Imin=BOUNDS(ng)%Istr (tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%JstrR(tile)
          Jmax=BOUNDS(ng)%JendR(tile)
        CASE (v3dvar)
          Imin=BOUNDS(ng)%IstrR(tile)
          Imax=BOUNDS(ng)%IendR(tile)
          Jmin=BOUNDS(ng)%Jstr (tile)
          Jmax=BOUNDS(ng)%JendR(tile)
      END SELECT
!
!  Determine number of slice depths (negative, m) to extract.
!
      mslice = SIZE(z_slice)      
!
!  Extract horizontal slices via linear interpolation.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
# ifdef MASKING
          my_mask=Amask(i,j)
# else
          my_mask=1.0_r8
# endif
          DO iz=1,mslice
            IF (my_mask.gt.0.5_r8) THEN
              Ztop=Adepth(i,j,UBK)
              Zbot=Adepth(i,j,LBK)
              IF (z_slice(iz).ge.Ztop) THEN    ! shallower than top grid
                k1=UBk                         ! cell. The slice is in 
                k2=UBK                         ! the upper cell half
                r1=1.0_r8                      ! shallower than top grid
                r2=0.0_r8                      ! above its middle depth
              ELSE IF (Zbot.ge.z_slice(iz)) THEN
                r1=0.0_r8                      ! If deeper, ignore
                r2=0.0_r8
              ELSE
                DO k=UBk,LBk+1,-1              ! Otherwise, interpolate
                  Ztop=Adepth(i,j,k  )         ! field at slice depth
                  Zbot=Adepth(i,j,k-1)
                  IF ((Ztop.gt.z_slice(iz)).and.(z_slice(iz).ge.Zbot)) THEN
                    k1=k-1
                    k2=k
                  END IF
                END DO
                dz=Adepth(i,j,k2)-Adepth(i,j,k1)
                r2=(z_slice(iz)-Adepth(i,j,k1))/dz
                r1=1.0_r8-r2
              END IF
              IF ((r1+r2).gt.0.0_r8) THEN
                Aslice(i,j,iz)=r1*Adata(i,j,k1) + r2*Adata(i,j,k2)
              ELSE
                Aslice(i,j,iz)=spval           ! unbounded slice depth
              END IF
            ELSE
              Aslice(i,j,iz)=0.0_r8            ! land masked value
            END IF 
          END DO
        END DO
      END DO
!
      RETURN
      END SUBROUTINE extract_slice
#endif
      END MODULE extract_slice_mod

