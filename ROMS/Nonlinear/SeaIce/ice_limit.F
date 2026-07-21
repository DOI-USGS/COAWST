#include "cppdefs.h"
      MODULE ice_limit_mod

#if defined ICE_ADVECT || defined ICE_THERMO
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                                      !
!  It imposes positiveness on ice tracer variables. Thus, it does not  !
!  allow negative ice concentration or thickness.                      !
!                                                                      !
!======================================================================!
!
      USE mod_param
      USE mod_ice
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE ice_bc2d_mod,    ONLY : ice_bc2d_tile
      USE ice_tibc_mod,    ONLY : ice_tibc_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  ice_limit
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_limit (ng, tile, model)
!***********************************************************************
!
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL ice_limit_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     liold(ng), linew(ng),                        &
     &                     ICE(ng) % Si)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_limit
!
!***********************************************************************
      SUBROUTINE ice_limit_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           liold, linew,                          &
     &                           Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(inout) :: liold, linew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
# endif
!
! Local variable definitions
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Impose positive concentration and thicknesses.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Si(i,j,linew,isAice)=MIN(Si(i,j,linew,isAice), max_ai(ng))
          Si(i,j,linew,isAice)=MAX(Si(i,j,linew,isAice), 0.0_r8)
          Si(i,j,linew,isHice)=MAX(Si(i,j,linew,isHice), 0.0_r8)
          Si(i,j,linew,isHsno)=MAX(Si(i,j,linew,isHsno), 0.0_r8)
          Si(i,j,linew,isHmel)=MAX(Si(i,j,linew,isHmel), 0.0_r8)
          Si(i,j,linew,isTice)=MAX(Si(i,j,linew,isTice),-70.0_r8)
          Si(i,j,linew,isTice)=MIN(Si(i,j,linew,isTice), 0.0_r8)
          Si(i,j,linew,isIage)=MAX(Si(i,j,linew,isIage), 0.0_r8)
          IF ((Si(i,j,linew,isAice).le.0.0_r8).or.                      &
     &        (Si(i,j,linew,isHice).le.0.0_r8)) THEN
            Si(i,j,linew,isAice)=0.0_r8
            Si(i,j,linew,isHice)=0.0_r8
            Si(i,j,linew,isHmel)=0.0_r8
            Si(i,j,linew,isHsno)=0.0_r8
            Si(i,j,linew,isIage)=0.0_r8
          END IF
# ifdef ICE_BIO
          Si(i,j,linew,isIphy)=MAX(Si(i,j,linew,isIphy), 0.0_r8)
          Si(i,j,linew,isINO3)=MAX(Si(i,j,linew,isINO3), 0.0_r8)
          Si(i,j,linew,isINH4)=MAX(Si(i,j,linew,isINH4), 0.0_r8)
          IF ((Si(i,j,linew,isAice).le.0.5_r8).or.                      &
     &        (Si(i,j,linew,isHice).le.0.02_r8)) THEN
            Si(i,j,linew,isIlog)=-1.0_r8
          ELSE
            Si(i,j,linew,isIlog)=1.0_r8
          END IF
# endif
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL ice_bc2d_tile (ng, tile, model, isAice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isAice),                             &
     &                    LBC(:,ibICE(isAice),ng))

      CALL ice_bc2d_tile (ng, tile, model, isHice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    LBC(:,ibICE(isHice),ng))

      CALL ice_bc2d_tile (ng, tile, model, isHsno,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHsno),                             &
     &                    LBC(:,ibICE(isHsno),ng))

      CALL ice_bc2d_tile (ng, tile, model, isHmel,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHmel),                             &
     &                    LBC(:,ibICE(isHmel),ng))

      CALL ice_bc2d_tile (ng, tile, model, isIage,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isIage),                             &
     &                    LBC(:,ibICE(isIage),ng))

# ifdef ICE_BIO
      CALL ice_bc2d_tile (ng, tile, model, isIphy,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isIphy),                             &
     &                    LBC(:,ibICE(isIphy),ng))

      CALL ice_bc2d_tile (ng, tile, model, isINO3,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isINO3),                             &
     &                    LBC(:,ibICE(isINO3),ng))

      CALL ice_bc2d_tile (ng, tile, iNLM, isINH4,                       &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isINH4),                             &
     &                    LBC(:,ibICE(isINH4),ng))
# endif

      CALL ice_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold, linew,                                 &
     &                    Si(:,:,:,isUice),                             &
     &                    Si(:,:,:,isVice),                             &
     &                    Si(:,:,:,isHice),                             &
     &                    Si(:,:,:,isTice),                             &
     &                    Si(:,:,:,isEnth))
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isAice))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHice))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHsno))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isTice))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isHmel))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isIage))

# ifdef ICE_BIO
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isIphy))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isINO3))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isINH4))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isIlog))
# endif
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,linew,isAice),                         &
     &                    Si(:,:,linew,isHice),                         &
     &                    Si(:,:,linew,isHsno),                         &
     &                    Si(:,:,linew,isHmel))

      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,linew,isTice))

#  ifdef ICE_BIO
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,linew,isIphy),                         &
     &                    Si(:,:,linew,isINO3),                         &
     &                    Si(:,:,linew,isINH4))
#  endif
# endif
!
      RETURN
      END SUBROUTINE ice_limit_tile
#endif
      END MODULE ice_limit_mod
