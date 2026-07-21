      MODULE ice_advect_mod
!
!git $Id
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                            W. Paul Budgell    !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine performs advection of ice scalars using the         !
!  Smolarkiewicz second-order upwind scheme.                           !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Smolarkiewicz, P.K. and W.W. Grabowski, 1990: The multidimensional  !
!    Positive definite advection transport algorithm: Nonoscillatory   !
!    option, J. Comp. Phys., 86, 355-375.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_forces
      USE mod_grid
      USE mod_ice
      USE mod_ocean
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
      USE ice_bc2d_mod,    ONLY : ice_bc2d_tile
      USE ice_tibc_mod,    ONLY : ice_tibc_tile
!
      implicit none
!
      PUBLIC  :: ice_advect
      PRIVATE :: ice_advect_tile
      PRIVATE :: ice_mpdata_tile
!
      CONTAINS
!
!***********************************************************************
       SUBROUTINE ice_advect (ng, tile, model)
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
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
#endif
      CALL ice_advect_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew(ng), liold(ng), liunw(ng))
#ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
#endif
!
      RETURN
      END SUBROUTINE ice_advect
!
!***********************************************************************
      SUBROUTINE ice_advect_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            linew, liold, liunw)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: linew, liold, liunw

!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Advect the ice concentation, isAice.
!-----------------------------------------------------------------------
!
      CALL ice_mpdata_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew, liold, liunw,                        &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#endif
#ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
#endif
#ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % Si(:,:,:,isUice),                 &
     &                      ICE(ng) % Si(:,:,:,isVice),                 &
     &                      ICE(ng) % Si(:,:,:,isAice))
!
      CALL ice_bc2d_tile (ng, tile, model, isAice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isAice),                   &
     &                    LBC(:,ibICE(isAice),ng))
!
!-----------------------------------------------------------------------
!  Advect the ice thickness, isHice.
!-----------------------------------------------------------------------
!
      CALL ice_mpdata_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew, liold, liunw,                        &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#endif
#ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
#endif
#ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % Si(:,:,:,isUice),                 &
     &                      ICE(ng) % Si(:,:,:,isVice),                 &
     &                      ICE(ng) % Si(:,:,:,isHice))
!
      CALL ice_bc2d_tile (ng, tile, model, isHice,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isHice),                   &
     &                    LBC(:,ibICE(isHice),ng))
!
!  Compute rate of ice divergence (m3/s).
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ICE(ng)%Fi(i,j,icWdiv)=(ICE(ng)%Si(i,j,linew,isHice)-         &
       &                          ICE(ng)%Si(i,j,liold,isHice))/dt(ng)
        END DO
      END DO

#ifdef ICE_THERMO
!
!-----------------------------------------------------------------------
!  Advect the snow thickness, isHsno.
!-----------------------------------------------------------------------
!
      CALL ice_mpdata_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew, liold, liunw,                        &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
# ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % Si(:,:,:,isUice),                 &
     &                      ICE(ng) % Si(:,:,:,isVice),                 &
     &                      ICE(ng) % Si(:,:,:,isHsno))
!
      CALL ice_bc2d_tile (ng, tile, model, isHsno,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isHsno),                   &
     &                    LBC(:,ibICE(isHsno),ng))
!
!-----------------------------------------------------------------------
!  Advect the surface melt water thickness, isHmel.
!-----------------------------------------------------------------------
!
      CALL ice_mpdata_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew, liold, liunw,                        &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
# ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % Si(:,:,:,isUice),                 &
     &                      ICE(ng) % Si(:,:,:,isVice),                 &
     &                      ICE(ng) % Si(:,:,:,isHmel))
!
      CALL ice_bc2d_tile (ng, tile, model, isHmel,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isHmel),                   &
     &                    LBC(:,ibICE(isHmel),ng))
!
!-----------------------------------------------------------------------
!  Advect the interior ice temperature, isTice.
!-----------------------------------------------------------------------
!
      ICE(ng) % Si(:,:,:,isEnth) = -ICE(ng) % Si(:,:,:,isEnth)
      CALL ice_mpdata_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew, liold, liunw,                        &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
# ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % Si(:,:,:,isUice),                 &
     &                      ICE(ng) % Si(:,:,:,isVice),                 &
     &                      ICE(ng) % Si(:,:,:,isEnth))
!
      ICE(ng) % Si(:,:,:,isEnth) = -ICE(ng) % Si(:,:,:,isEnth)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ICE(ng)%Si(i,j,linew,isTice)=ICE(ng)%Si(i,j,linew,isEnth)/    &
     &                                 MAX(ICE(ng)%Si(i,j,linew,isHice),&
     &                                     1.0E-6_r8)
          IF (ICE(ng)%Si(i,j,linew,isAice).le.min_ai(ng)) THEN
            ICE(ng)%Si(i,j,linew,isEnth)=0.0_r8
            ICE(ng)%Si(i,j,linew,isTice)=-2.0_r8
          END IF
        END DO
      END DO
!
      CALL ice_tibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold, linew,                                 &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isHice),                   &
     &                    ICE(ng) % Si(:,:,:,isTice),                   &
     &                    ICE(ng) % Si(:,:,:,isEnth))
!
!-----------------------------------------------------------------------
!  Advect thickness associated with age of ice, isHage. Then, compute
!  ice age (isIage).
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ICE(ng)%Si(i,j,liold,isHage)=ICE(ng)%Si(i,j,liold,isHice)*    &
     &                                 ICE(ng)%Si(i,j,liold,isIage)
          ICE(ng)%Si(i,j,linew,isHage)=ICE(ng)%Si(i,j,linew,isHice)*    &
     &                                 ICE(ng)%Si(i,j,linew,isIage)
          IF (ICE(ng)%Si(i,j,liold,isHice).le.min_hi(ng)) THEN
            ICE(ng)%Si(i,j,liold,isHage)=0.0_r8
            ICE(ng)%Si(i,j,liold,isIage)=0.0_r8
          END IF
        END DO
      END DO
!
      CALL ice_mpdata_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew, liold, liunw,                        &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
# ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % Si(:,:,:,isUice),                 &
     &                      ICE(ng) % Si(:,:,:,isVice),                 &
     &                      ICE(ng) % Si(:,:,:,isHage))
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ICE(ng)%Si(i,j,linew,isIage)=ICE(ng)%Si(i,j,linew,isHage)/    &
     &                                 MAX(ICE(ng)%Si(i,j,linew,isHice),&
     &                                     1.0E-6_r8)
          IF (ICE(ng)%Si(i,j,linew,isHice).le.min_hi(ng)) THEN
            ICE(ng)%Si(i,j,linew,isHage)=0.0_r8
            ICE(ng)%Si(i,j,linew,isIage)=0.0_r8
          END IF
        END DO
      END DO
!
      CALL ice_bc2d_tile (ng, tile, model, isIage,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    liold, linew,                                 &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isIage),                   &
     &                    LBC(:,ibICE(isIage),ng))
#endif

#if defined ICE_THERMO && defined ICE_BIO
!
!-----------------------------------------------------------------------
!  Advect the ice algae concentration, IcePhL.
!-----------------------------------------------------------------------
!
      CALL ice_mpdata_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew, liold, liunw,                        &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
# ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % Si(:,:,:,isUice),                 &
     &                      ICE(ng) % Si(:,:,:,isVice),                 &
     &                      ICE(ng) % Si(:,:,:,isIphy))
!
! Need to change this to ice_bc2d_tile calls.
!
      CALL ice_bc2d_tile (ng, tile, model, isIphy,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold, linew,                                 &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isIphy),                   &
     &                    LBC(:,ibICE(isIphy),ng))
!
!-----------------------------------------------------------------------
!  Advect the ice nitrate concentration, isINO3.
!-----------------------------------------------------------------------
!
      CALL ice_mpdata_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew, liold, liunw,                        &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
# ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % Si(:,:,:,isUice),                 &
     &                      ICE(ng) % Si(:,:,:,isVice),                 &
     &                      ICE(ng) % Si(:,:,:,isINO3))
!
      CALL ice_bc2d_tile (ng, tile, model, isINO3,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold, linew,                                 &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isINO3),                   &
     &                    LBC(:,ibICE(isINO3),ng))
!
!-----------------------------------------------------------------------
!  Advect the ice ammonia concentration, isINH4.
!-----------------------------------------------------------------------
!
      CALL ice_mpdata_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew, liold, liunw,                        &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
# ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % Si(:,:,:,isUice),                 &
     &                      ICE(ng) % Si(:,:,:,isVice),                 &
     &                      ICE(ng) % Si(:,:,:,isINH4))
!
      CALL ice_bc2d_tile (ng, tile, model, isINH4,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold, linew,                                 &
     &                    ICE(ng) % Si(:,:,:,isUice),                   &
     &                    ICE(ng) % Si(:,:,:,isVice),                   &
     &                    ICE(ng) % Si(:,:,:,isINH4),                   &
     &                    LBC(:,ibICE(isINH4),ng))
#endif
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isAice))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isHice))

#ifdef ICE_THERMO
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isHsno))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isHmel))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isTice))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isEnth))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isIage))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isHage))

# ifdef ICE_BIO
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isIphy))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isINO3))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%Si(:,:,linew,isINH4))
# endif
#endif
      END IF

#ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%Si(:,:,linew,isAice),                 &
     &                    ICE(ng)%Si(:,:,linew,isHice))

# ifdef ICE_THERMO
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%Si(:,:,linew,isHsno),                 &
     &                    ICE(ng)%Si(:,:,linew,isHmel))

      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%Si(:,:,linew,isTice),                 &
     &                    ICE(ng)%Si(:,:,linew,isEnth),                 &
     &                    ICE(ng)%Si(:,:,linew,isIage),                 &
     &                    ICE(ng)%Si(:,:,linew,isHage))

#  if defined ICE_BIO
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%Si(:,:,linew,isIphy),                 &
     &                    ICE(ng)%Si(:,:,linew,isINO3),                 &
     &                    ICE(ng)%Si(:,:,linew,isINH4))
#  endif
# endif
#endif
!
      RETURN
      END SUBROUTINE ice_advect_tile
!
!***********************************************************************
      SUBROUTINE ice_mpdata_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            linew, liold, liunw,                  &
#ifdef MASKING
     &                            rmask,  umask, vmask,                 &
#endif
#ifdef WET_DRY
     &                            rmask_wet, umask_wet, vmask_wet,      &
#endif
#ifdef ICESHELF
     &                            zice,                                 &
#endif
#ifndef ICE_UPWIND
     &                            pm, pn,                               &
#endif
     &                            on_u, om_v, omn,                      &
     &                            ui, vi, field)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: linew, liold, liunw
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
# endif
# ifndef ICE_UPWIND
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
# endif
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: field(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
# endif
# ifndef ICE_UPWIND
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: field(LBi:UBi,LBj:UBj,2)
#endif
!
! Local variable definitions
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
!
      real(r8) :: Cu_crss, Cu
      real(r8) :: cff1, cff2, rateu, ratev, rateyiu, ratexiv
      real(r8) :: uspeed, vspeed
!
      real(r8), parameter :: epsil = 1.0E-15_r8
      real(r8), parameter :: add = 3.0E+3_r8
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aflxu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aflxv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Horizontal advection of generic ice field.
!-----------------------------------------------------------------------
!
#ifndef ICE_UPWIND
      IF (EWperiodic(ng)) THEN
        Imin=Istr-1
        Imax=Iend+1
      ELSE
        Imin=MAX(Istr-1,1)
        Imax=MIN(Iend+1,Lm(ng))
      END IF
      IF (NSperiodic(ng)) THEN
        Jmin=Jstr-1
        Jmax=Jend+1
      ELSE
        Jmin=MAX(Jstr-1,1)
        Jmax=MIN(Jend+1,Mm(ng))
      END IF
#else
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
#endif
!
!  Upstream fluxes. (HGA: the advection is not in flux form?)
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax+1
          cff1=MAX(0.0_r8,ui(i,j,liunw))
          cff2=MIN(0.0_r8,ui(i,j,liunw))
          aflxu(i,j)=on_u(i,j)*                                         &
     &               (cff1*field(i-1,j,liold)+                          &
     &                cff2*field(i  ,j,liold))
        END DO
      END DO
      DO j=Jmin,Jmax+1
        DO i=Imin,Imax
          cff1=MAX(0.0_r8,vi(i,j,liunw))
          cff2=MIN(0.0_r8,vi(i,j,liunw))
          aflxv(i,j)=om_v(i,j)*                                         &
     &               (cff1*field(i,j-1,liold)+                          &
     &                cff2*field(i,j  ,liold))
        END DO
      END DO
!
!  Step 1: first guess advection operator.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ar(i,j)=1.0_r8/omn(i,j)
          aif(i,j)=(field(i,j,liold)-                                   &
     &              dtice(ng)*(aflxu(i+1,j)-aflxu(i,j)+                 &
     &                         aflxv(i,j+1)-aflxv(i,j))*ar(i,j))
#ifdef MASKING
          aif(i,j)=aif(i,j)*rmask(i,j)
#endif
#ifdef WET_DRY
          aif(i,j)=aif(i,j)*rmask_wet(i,j)
#endif
#ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) aif(i,j)=0.0_r8
#endif
        END DO
      END DO
!
!  Set values at the open boundaries.
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jmin,Jmax
            aif(Istr-1,j)=aif(Istr,j)
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jmin,Jmax
            aif(Iend+1,j)=aif(Iend,j)
          END DO
        END IF
      END IF
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Imin,Imax
            aif(i,Jstr-1)=aif(i,Jstr)
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Imin,Imax
            aif(i,Jend+1)=aif(i,Jend)
          END DO
        END IF
        IF (.not.EWperiodic(ng)) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            aif(Istr-1,Jstr-1)=aif(Istr,Jstr)
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            aif(Istr-1,Jend+1)=aif(Istr,Jend)
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            aif(Iend+1,Jstr-1)=aif(Iend,Jstr)
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            aif(Iend+1,Jend+1)=aif(Iend,Jend)
          END IF
        END IF
      END IF

#ifdef MASKING
!
!  Apply land/sea mask.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          aif(i,j)=aif(i,j)*rmask(i,j)
# ifdef WET_DRY
          aif(i,j)=aif(i,j)*rmask_wet(i,j)
# endif
        END DO
      END DO
#endif
#ifdef ICESHELF
      DO j=Jmin,Jmax
        DO i=Imin,Jmax
          IF (zice(i,j).ne.0.0_r8) THEN
            aif(i,j)=0.0_r8
          END IF
        END DO
      END DO
#endif

#ifndef ICE_UPWIND
!
!  Step 2: Compute antidiffusive velocities .
!
!  This is needed to avoid touching "aif" under land mask.
!  Note that only aif(i,j) and aif(i-1,j) are allowed to appear
!  explicitly in the code segment below. This is okay because if
!  either of them is masked, then "ui" is zero at that point,
!  and therefore no additional masking is required.
!
      DO j=Jstr,Jend+1
        DO i=Istr-1,Iend+1
          FE(i,j)=0.5*                                                  &
# ifdef MASKING
     &            vmask(i,j)*                                           &
# endif
# ifdef WET_DRY
     &            vmask_wet(i,j)*                                       &
# endif
     &            (aif(i,j)-aif(i,j-1))
        END DO
      END DO
!
      DO j=Jstr-1,Jend+1
        DO i=Istr,Iend+1
          FX(i,j)=0.5*                                                  &
# ifdef MASKING
     &            umask(i,j)*                                           &
# endif
# ifdef WET_DRY
     &            umask_wet(i,j)*                                       &
# endif
     &            (aif(i,j)-aif(i-1,j))
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend+1
          rateu=(aif(i,j)-aif(i-1,j))/                                  &
     &          MAX(epsil, aif(i,j)+aif(i-1,j))

          rateyiu=(FE(i,j+1)+FE(i,j)+FE(i-1,j+1)+FE(i-1,j))/            &
     &            (MAX(epsil, aif(i  ,j)+FE(i  ,j+1)-FE(i  ,j)+         &
     &                        aif(i-1,j)+FE(i-1,j+1)-FE(i-1,j)))

          Cu=0.5*dtice(ng)*(pm(i,j)+pm(i-1,j))*ui(i,j,liunw)

          Cu_crss=0.5_r8*dtice(ng)*                                     &
     &            0.0625_r8*(pn(i-1,j+1)+pn(i,j+1)+                     &
     &                       pn(i-1,j-1)+pn(i,j-1))*                    &
     &            (vi(i-1,j+1,liunw)+vi(i,j+1,liunw)+                   &
     &             vi(i-1,j  ,liunw)+vi(i,j  ,liunw))

          uspeed=rateu*(ABS(ui(i,j,liunw))-Cu*ui(i,j,liunw))-           &
     &           rateyiu*Cu_crss*ui(i,j,liunw)

          cff1=MAX(0.0_r8,uspeed)
          cff2=MIN(0.0_r8,uspeed)
          aflxu(i,j)=on_u(i,j)*(cff1*aif(i-1,j)+                        &
     &                          cff2*aif(i  ,j))
        END DO
      END DO
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend
          ratev=(aif(i,j)-aif(i,j-1))/                                  &
     &          MAX(epsil, aif(i,j)+aif(i,j-1))

          ratexiv=(FX(i+1,j)+FX(i,j) +FX(i+1,j-1)+FX(i,j-1))/           &
     &            (MAX(epsil, aif(i,j  )+FX(i+1,j  )-FX(i,j  )+         &
     &                        aif(i,j-1)+FX(i+1,j-1)-FX(i,j-1)))

          Cu=0.5*dtice(ng)*(pn(i,j)+pn(i,j-1))*vi(i,j,liunw)

          Cu_crss=0.5_r8*dtice(ng)*                                     &
     &            0.0625_r8*(pm(i+1,j)+pm(i+1,j-1)+                     &
     &                       pm(i-1,j)+pm(i-1,j-1))*                    &
     &            (ui(i,j  ,liunw)+ui(i+1,j  ,liunw)+                   &
     &             ui(i,j-1,liunw)+ui(i+1,j-1,liunw))

          vspeed=ratev*(ABS(vi(i,j,liunw))-Cu*vi(i,j,liunw))-           &
     &           ratexiv*Cu_crss*vi(i,j,liunw)

          cff1=MAX(0.0_r8,vspeed)
          cff2=MIN(0.0_r8,vspeed)
          aflxv(i,j)=om_v(i,j)*(cff1*aif(i,j-1)+                        &
     &                          cff2*aif(i,j  ))
        END DO
      END DO
!
!  Correct advection operator by substraction the antidifusive fluxes.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          aif(i,j)=aif(i,j)-                                            &
     &             dtice(ng)*pm(i,j)*pn(i,j)*                           &
     &             (aflxu(i+1,j)-aflxu(i,j)+                            &
     &              aflxv(i,j+1)-aflxv(i,j))
# ifdef MASKING
          aif(i,j)=aif(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          aif(i,j)=aif(i,j)*rmask_wet(i,j)
# endif
# ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) aif(i,j)=0.
# endif
        END DO
      END DO
#endif
!
!  Load advected solution.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          field(i,j,linew)=aif(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ice_mpdata_tile

      END MODULE ice_advect_mod
