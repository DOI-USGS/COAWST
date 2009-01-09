#include "cppdefs.h"
      MODULE set_zeta_mod

#ifdef SOLVE3D
!
!svn $Id: set_zeta.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets free-surface to its fast-time averaged value.     !
!                                                                      !
!=======================================================================
!
!
      implicit none

      PRIVATE
      PUBLIC  :: set_zeta

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_zeta (ng, tile)
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
      CALL wclock_on (ng, iNLM, 12)
# endif
      CALL set_zeta_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    COUPLING(ng) % Zt_avg1,                       &
     &                    OCEAN(ng) % zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 12)
# endif
      RETURN
      END SUBROUTINE set_zeta
!
!***********************************************************************
      SUBROUTINE set_zeta_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Zt_avg1, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)

      real(r8), intent(out) :: zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
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
          zeta(i,j,1)=Zt_avg1(i,j)
          zeta(i,j,2)=Zt_avg1(i,j)
        END DO
      END DO
# if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,1))
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        zeta(:,:,2))
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    zeta(:,:,1),                                  &
     &                    zeta(:,:,2))
# endif

      RETURN
      END SUBROUTINE set_zeta_tile
#endif
      END MODULE set_zeta_mod
