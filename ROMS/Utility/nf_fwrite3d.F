#include "cppdefs.h"
      FUNCTION nf_fwrite3d (ng, model, ncid, ncvarid, tindex, gtype,    &
     &                      LBi, UBi, LBj, UBj, LBk, UBk, Ascl,         &
#ifdef MASKING
     &                      Amask,                                      &
#endif
     &                      A)
!
!svn $Id: nf_fwrite3d.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function writes out a generic floating point 3D array into an  !
!  output NetCDF file.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
!     ncvarid      NetCDF variable ID.                                 !
!     tindex       NetCDF time record index to write.                  !
!     gtype        Grid type. If negative, only write water points.    !
!     LBi          I-dimension Lower bound.                            !
!     UBi          I-dimension Upper bound.                            !
!     LBj          J-dimension Lower bound.                            !
!     UBj          J-dimension Upper bound.                            !
!     LBk          K-dimension Lower bound.                            !
!     UBk          K-dimension Upper bound.                            !
!     Amask        land/Sea mask, if any (real).                       !
!     Ascl         Factor to scale field before writing (real).        !
!     A            Field to write out (real).                          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nf_fwrite3d  Error flag (integer).                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_gather
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

      real(r8), intent(in) :: Ascl

#ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: i, j, k, ic, Npts
      integer :: Imin, Imax, Jmin, Jmax, Kmin, Kmax, Koff
      integer :: Ilen, Jlen, Klen, IJlen, MyType, status

      integer, dimension(4) :: start, total

      integer :: nf_fwrite3d

#if defined INLINE_2DIO && defined DISTRIBUTE
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Aout
#elif defined COMPOSED_GRID
      real(r8), dimension((Lm(ng)+7)*(Mm(ng)+7)*(UBk-LBk+1)) :: Aout
#elif defined REFINED_GRID
      real(r8), allocatable :: Aout(:)
#else
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: Aout
#endif
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Set first and last grid point according to staggered C-grid
!  classification. Set loops offsets.
!
#ifdef COMPOSED_GRID
      Imin=LOWER_BOUND_I
      Imax=Lm(ng)+GHOST_POINTS
      Jmin=LOWER_BOUND_J
      Jmax=Mm(ng)+GHOST_POINTS
      MyType=ABS(gtype)
      IF ((MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) THEN
        Imin=LOWER_BOUND_I+1
      END IF
      IF ((MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) THEN
        Jmin=LOWER_BOUND_J+1
      END IF
#elif defined REFINED_GRID
      IF (ng.eq.1) THEN
        Imin=0
        Imax=Lm(ng)+1
        Jmin=0
        Jmax=Mm(ng)+1
        MyType=ABS(gtype)
        IF ((MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) THEN
          Imin=1
        END IF
        IF ((MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) THEN
          Jmin=1
        END IF
      ELSE
        Imin=LOWER_BOUND_I
        Imax=Lm(ng)+GHOST_POINTS
        Jmin=LOWER_BOUND_J
        Jmax=Mm(ng)+GHOST_POINTS
        MyType=ABS(gtype)
        IF ((MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) THEN
          Imin=LOWER_BOUND_I+1
        END IF
        IF ((MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) THEN
          Jmin=LOWER_BOUND_J+1
        END IF
      END IF
#else
      Imin=0
      Imax=Lm(ng)+1
      Jmin=0
      Jmax=Mm(ng)+1
      MyType=ABS(gtype)
      IF ((MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) THEN
        Imin=1
      END IF
      IF ((MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) THEN
        Jmin=1
      END IF
#endif
      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF
      Ilen=Imax-Imin+1
      Jlen=Jmax-Jmin+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen

#ifdef REFINED_GRID
      Npts=IJlen*Klen
      IF (.not.allocated(Aout)) THEN
        allocate (Aout(IJlen*Klen))
      END IF
#endif

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  If distributed-memory set-up, collect tile data from all spawned
!  nodes and store it into a global scratch 1D array, packed in column-
!  major order. If applicable, remove land points.
!-----------------------------------------------------------------------
!
# ifdef INLINE_2DIO

!  If appropriate, process 3D data level by level to reduce memory
!  requirements.
!
      Kmin=1
      Kmax=1
      DO k=LBk,UBk
# else
        Kmin=LBk
        Kmax=UBk
# endif
        CALL mp_gather (ng, model, LBi, UBi, LBj, UBj, Kmin, Kmax,      &
     &                  gtype, Ascl,                                    &
# ifdef MASKING
     &                  Amask,                                          &
# endif
# ifdef INLINE_2DIO
     &                  A(LBi,LBj,k), Npts, Aout)
# else
     &                  A, Npts, Aout)
# endif
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
        nf_fwrite3d=nf90_noerr
        IF (OutThread) THEN
          IF (gtype.gt.0) THEN
            start(1)=1
            total(1)=Ilen
            start(2)=1
            total(2)=Jlen
# ifdef INLINE_2DIO
            start(3)=k-Koff+1
            total(3)=1
# else
            start(3)=1
            total(3)=Klen
# endif
            start(4)=tindex
            total(4)=1
# ifdef MASKING
          ELSE
            start(1)=1
            total(1)=Npts
            start(2)=tindex
            total(2)=1
# endif
          END IF
          status=nf90_put_var(ncid, ncvarid, Aout, start, total)
          nf_fwrite3d=status
        END IF
# ifdef INLINE_2DIO
      END DO
# endif
#else
!
!-----------------------------------------------------------------------
!  If serial or shared-memory applications and serial output, pack data
!  into a global 1D array in column-major order. If applicable, remove
!  land points.
!-----------------------------------------------------------------------
!
      IF (gtype.gt.0) THEN
        ic=0
        Npts=IJlen*Klen
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              Aout(ic)=A(i,j,k)*Ascl
            END DO
          END DO
        END DO
# ifdef MASKING
      ELSE
        Npts=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              IF (Amask(i,j).gt.0.0_r8) THEN
                Npts=Npts+1
                Aout(Npts)=A(i,j,k)*Ascl
              END IF
            END DO
          END DO
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
      nf_fwrite3d=nf90_noerr
      IF (OutThread) THEN
        IF (gtype.gt.0) THEN
          start(1)=1
          total(1)=Ilen
          start(2)=1
          total(2)=Jlen
          start(3)=1
          total(3)=Klen
          start(4)=tindex
          total(4)=1
# ifdef MASKING
        ELSE
          start(1)=1
          total(1)=Npts
          start(2)=tindex
          total(2)=1
# endif
        END IF
        status=nf90_put_var(ncid, ncvarid, Aout, start, total)
      END IF
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast IO error flag to all nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcasti (ng, model, status, 1)
#endif
      nf_fwrite3d=status
#ifdef REFINED_GRID
      IF (allocated(Aout)) THEN
        deallocate (Aout)
      END IF
#endif
      RETURN
      END FUNCTION nf_fwrite3d
