#include "cppdefs.h"
      MODULE tl_set_zeta_mod

#if defined TANGENT && defined SOLVE3D
!
!svn $Id: tl_set_zeta.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets tangent linear free-surface to its fast-time      !
!  averaged value.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_set_zeta
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_set_zeta (ng, tile)
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
      CALL wclock_on (ng, iTLM, 12)
# endif
      CALL tl_set_zeta_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       COUPLING(ng) % tl_Zt_avg1,                 &
     &                       OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 12)
# endif

      RETURN
      END SUBROUTINE tl_set_zeta
!
!***********************************************************************
      SUBROUTINE tl_set_zeta_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             tl_Zt_avg1, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: tl_Zt_avg1(LBi:,LBj:)

      real(r8), intent(out) :: tl_zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Prepare to time-step 2D equations:  set initial free-surface
!  to its fast-time averaged values (which corresponds to the time
!  step "n").
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!>        zeta(i,j,1)=Zt_avg1(i,j)
!>        zeta(i,j,2)=Zt_avg1(i,j)
!>
          tl_zeta(i,j,1)=tl_Zt_avg1(i,j)
          tl_zeta(i,j,2)=tl_Zt_avg1(i,j)
        END DO
      END DO

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          zeta(:,:,1))
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_zeta(:,:,1))
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          zeta(:,:,2))
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_zeta(:,:,2))
      END IF

# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    zeta(:,:,1),                                  &
!>   &                    zeta(:,:,2))
!>
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_zeta(:,:,1),                               &
     &                    tl_zeta(:,:,2))
# endif

      RETURN
      END SUBROUTINE tl_set_zeta_tile
#endif
      END MODULE tl_set_zeta_mod
