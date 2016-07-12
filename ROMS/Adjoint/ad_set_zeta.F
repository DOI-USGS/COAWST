#include "cppdefs.h"
      MODULE ad_set_zeta_mod

#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_set_zeta.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets adjoint free-surface to its fast-time averaged    !
!  value.                                                              !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ad_set_zeta
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_set_zeta (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 12)
# endif
      CALL ad_set_zeta_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       COUPLING(ng) % ad_Zt_avg1,                 &
     &                       OCEAN(ng) % ad_zeta,                       &
     &                       OCEAN(ng) % ad_zeta_sol)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 12)
# endif
      RETURN
      END SUBROUTINE ad_set_zeta
!
!***********************************************************************
      SUBROUTINE ad_set_zeta_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             ad_Zt_avg1, ad_zeta, ad_zeta_sol)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)

      real(r8), intent(out) :: ad_zeta_sol(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(out) :: ad_zeta_sol(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Prepare to time-step adjoint 2D equations:  set initial free-surface
!  to its fast-time averaged values (which corresponds to the time
!  step "n").
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_zeta(:,:,1),                               &
!>   &                    tl_zeta(:,:,2))
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,1),                            &
     &                       ad_zeta(:,:,2))
!
# endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_zeta(:,:,2))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,2))
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_zeta(:,:,1))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,1))
      END IF
!
!  Save the adjoint free-surface solution for the current time-step
!  (sum of time indices 1 and 2).
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta_sol(i,j)=ad_zeta(i,j,1)+ad_zeta(i,j,2)
        END DO
      END DO
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!>        tl_zeta(i,j,2)=tl_Zt_avg1(i,j)
!>        tl_zeta(i,j,1)=tl_Zt_avg1(i,j)
!>
          ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+                              &
     &                    ad_zeta(i,j,1)+ad_zeta(i,j,2)
          ad_zeta(i,j,2)=0.0_r8
          ad_zeta(i,j,1)=0.0_r8
        END DO
      END DO

      RETURN
      END SUBROUTINE ad_set_zeta_tile
#endif
      END MODULE ad_set_zeta_mod
