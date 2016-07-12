#include "cppdefs.h"
      MODULE uv_rotate_mod
!
!svn $Id: uv_rotate.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!=======================================================================
!                                                                      !
!  These routines average momentum component to RHO-points and then    !
!  rotates from (XI,ETA) coordinates to geographical Eastward and      !
!  Northward directions.                                               !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: uv_rotate2d
#ifdef SOLVE3D
      PUBLIC  :: uv_rotate3d
#endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE uv_rotate2d (ng, tile, add, Lboundary,                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CosAngler, SinAngler,                     &
#ifdef MASKING
     &                        rmask_full,                               &
#endif
     &                        Uinp, Vinp, Uout, Vout)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: add, Lboundary

      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
# ifdef MASKING
      real(r8), intent(in) :: rmask_full(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Uinp(LBi:,LBj:)
      real(r8), intent(in) :: Vinp(LBi:,LBj:)

      real(r8), intent(inout) :: Uout(LBi:,LBj:)
      real(r8), intent(inout) :: Vout(LBi:,LBj:)
#else
      real(r8), intent(in) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: SinAngler(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: rmask_full(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Uinp(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vinp(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: Uout(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Vout(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: Urho, Vrho

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Rotate 2D vector components to Eastward and Northward directions.
!-----------------------------------------------------------------------
!
      IF (add) THEN
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Urho=0.5_r8*(Uinp(i,j)+Uinp(i+1,j))
            Vrho=0.5_r8*(Vinp(i,j)+Vinp(i,j+1))
            Uout(i,j)=Uout(i,j)+                                        &
     &                Urho*CosAngler(i,j)-                              &
     &                Vrho*SinAngler(i,j)
            Vout(i,j)=Vout(i,j)+                                        &
     &                Vrho*CosAngler(i,j)+                              &
     &                Urho*SinAngler(i,j)
#ifdef MASKING
            Uout(i,j)=Uout(i,j)*rmask_full(i,j)
            Vout(i,j)=Vout(i,j)*rmask_full(i,j)
#endif
          END DO
        END DO
      ELSE
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Urho=0.5_r8*(Uinp(i,j)+Uinp(i+1,j))
            Vrho=0.5_r8*(Vinp(i,j)+Vinp(i,j+1))
            Uout(i,j)=Urho*CosAngler(i,j)-                              &
     &                Vrho*SinAngler(i,j)
            Vout(i,j)=Vrho*CosAngler(i,j)+                              &
     &                Urho*SinAngler(i,j)
#ifdef MASKING
            Uout(i,j)=Uout(i,j)*rmask_full(i,j)
            Vout(i,j)=Vout(i,j)*rmask_full(i,j)
#endif
          END DO
        END DO
      END IF
!
!  Exchange boundary data, if applicable.
!
      IF (Lboundary) THEN
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Uout)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Vout)
#ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, iNLM, 2,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        Uout, Vout)
#endif
        END IF
      END IF

      RETURN
      END SUBROUTINE uv_rotate2d

#ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE uv_rotate3d (ng, tile, add, Lboundary,                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        CosAngler, SinAngler,                     &
# ifdef MASKING
     &                        rmask_full,                                 &
# endif
     &                        Uinp, Vinp, Uout, Vout)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: add, Lboundary

      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask_full(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Uinp(LBi:,LBj:,LBk:)
      real(r8), intent(in) :: Vinp(LBi:,LBj:,LBk:)

      real(r8), intent(inout) :: Uout(LBi:,LBj:,LBk:)
      real(r8), intent(inout) :: Vout(LBi:,LBj:,LBk:)
# else
      real(r8), intent(in) :: CosAngler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: SinAngler(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask_full(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Uinp(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(in) :: Vinp(LBi:UBi,LBj:UBj,LBk:UBk)

      real(r8), intent(inout) :: Uout(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(inout) :: Vout(LBi:UBi,LBj:UBj,LBk,UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: Urho, Vrho

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Rotate 3D vector components to Eastward and Northward directions.
!-----------------------------------------------------------------------
!
      IF (add) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
            DO i=Istr,Iend
              Urho=0.5_r8*(Uinp(i,j,k)+Uinp(i+1,j,k))
              Vrho=0.5_r8*(Vinp(i,j,k)+Vinp(i,j+1,k))
              Uout(i,j,k)=Uout(i,j,k)+                                  &
     &                    Urho*CosAngler(i,j)-                          &
     &                    Vrho*SinAngler(i,j)
              Vout(i,j,k)=Vout(i,j,k)+                                  &
     &                    Vrho*CosAngler(i,j)+                          &
     &                    Urho*SinAngler(i,j)
# ifdef MASKING
              Uout(i,j,k)=Uout(i,j,k)*rmask_full(i,j)
              Vout(i,j,k)=Vout(i,j,k)*rmask_full(i,j)
# endif
            END DO
          END DO
        END DO
      ELSE
        DO k=LBk,UBk
          DO j=Jstr,Jend
            DO i=Istr,Iend
              Urho=0.5_r8*(Uinp(i,j,k)+Uinp(i+1,j,k))
              Vrho=0.5_r8*(Vinp(i,j,k)+Vinp(i,j+1,k))
              Uout(i,j,k)=Urho*CosAngler(i,j)-                          &
     &                    Vrho*SinAngler(i,j)
              Vout(i,j,k)=Vrho*CosAngler(i,j)+                          &
     &                    Urho*SinAngler(i,j)
# ifdef MASKING
              Uout(i,j,k)=Uout(i,j,k)*rmask_full(i,j)
              Vout(i,j,k)=Vout(i,j,k)*rmask_full(i,j)
# endif
            END DO
          END DO
        END DO
      END IF
!
!  Exchange boundary data.
!
      IF (Lboundary) THEN
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            Uout)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            Vout)
# ifdef DISTRIBUTE
          CALL mp_exchange3d (ng, tile, iNLM, 2,                        &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        Uout, Vout)
# endif
        END IF
      END IF

      RETURN
      END SUBROUTINE uv_rotate3d
#endif
      END MODULE uv_rotate_mod
