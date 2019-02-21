       SUBROUTINE ice_advect (ng, tile)
!
!*************************************************** W. Paul Budgell ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
!************************************************** Hernan G. Arango ***
!                                                                      !
!  This subroutine performs advection of ice scalars using the         !
!  Smolarkiewicz second-order upwind scheme.                           !
!  Reference:                                                          !
!  Smolarkiewicz and Grabowski (1990)                                  !
!***********************************************************************
!
      USE mod_param

      implicit none

      integer, intent(in) :: ng, tile

#include "tile.h"

#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 56)
#endif
! ---------------------------------------------------------------------
!  Advect the ice concentration.
! ---------------------------------------------------------------------
      CALL ice_advect_all_tile (ng, tile,                               &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS                  &
     &                      )
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 56)
#endif
      RETURN
      END SUBROUTINE ice_advect

      SUBROUTINE ice_advect_all_tile (ng, tile,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS                  &
     &                      )

      USE mod_param
      USE mod_ncparam
      USE mod_grid
      USE mod_ocean
      USE mod_ice
      USE mod_forces
      USE mod_scalars
      USE mod_stepping
      USE mod_boundary
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
      USE i2d_bc_mod
      USE tibc_mod, ONLY : tibc_tile
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variables
!
      integer :: i, j
      real(r8) :: wrk(LBi:UBi, LBj:UBj, 2)

#include "set_bounds.h"

      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
#ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
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
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % ai                                &
     &                      )
!
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
     &                  BOUNDARY(ng)%ai_west(LBj:UBj),                  &
     &                  BOUNDARY(ng)%ai_east(LBj:UBj),                  &
     &                  BOUNDARY(ng)%ai_north(LBi:UBi),                 &
     &                  BOUNDARY(ng)%ai_south(LBi:UBi),                 &
     &                  ICE(ng)%ui,                                     &
     &                  ICE(ng)%vi,                                     &
     &                  ICE(ng)%ai, LBC(:,isAice,ng))
!
! ---------------------------------------------------------------------
!  Advect the ice thickness.
! ---------------------------------------------------------------------
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
#ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
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
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % hi                                &
     &                      )
!
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
     &                  BOUNDARY(ng)%hi_west(LBj:UBj),                  &
     &                  BOUNDARY(ng)%hi_east(LBj:UBj),                  &
     &                  BOUNDARY(ng)%hi_north(LBi:UBi),                 &
     &                  BOUNDARY(ng)%hi_south(LBi:UBi),                 &
     &                  ICE(ng)%ui,                                     &
     &                  ICE(ng)%vi,                                     &
     &                  ICE(ng)%hi, LBC(:,isHice,ng))
!
! Store ice divergence rate
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ICE(ng)%wdiv(i,j) = (ICE(ng)%hi(i,j,linew(ng))-               &
       &                       ICE(ng)%hi(i,j,liold(ng)))/dt(ng)
        ENDDO
      ENDDO
!
! ---------------------------------------------------------------------
!  Advect the snow thickness.
! ---------------------------------------------------------------------
#ifdef ICE_THERMO
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % hsn                               &
     &                      )
!
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
     &                  BOUNDARY(ng)%hsn_west(LBj:UBj),                 &
     &                  BOUNDARY(ng)%hsn_east(LBj:UBj),                 &
     &                  BOUNDARY(ng)%hsn_north(LBi:UBi),                &
     &                  BOUNDARY(ng)%hsn_south(LBi:UBi),                &
     &                  ICE(ng)%ui,                                     &
     &                  ICE(ng)%vi,                                     &
     &                  ICE(ng)%hsn, LBC(:,isHsno,ng))
# if defined ICE_BIO && defined BERING_10K
!  ---------------------------------------------------------------------
!  Advect the icePhL
!  ---------------------------------------------------------------------
       
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % IcePhL                            &
     &                      )
!  FOOOO 
! Need to change this to i2d_bc_tile calls
      CALL IcePhLbc_tile (ng, tile,                                   &
     &                          LBi, UBi, LBj, UBj,                      &
     &                          liold(ng), linew(ng),                   &
     &                          ICE(ng)%ui,                             &
     &                          ICE(ng)%vi,                             &
     &                          ICE(ng)%IcePhL)

!  ---------------------------------------------------------------------
!  Advect the iceNH4
!  ---------------------------------------------------------------------

      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % IceNH4                            &
     &                      )
!
        CALL IceNH4bc_tile (ng, tile,                                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          liold(ng), linew(ng),                   &
     &                          ICE(ng)%ui,                             &
     &                          ICE(ng)%vi,                             &
     &                          ICE(ng)%IceNH4)
!

!
!  ---------------------------------------------------------------------
!  Advect the iceNO3.
!  ---------------------------------------------------------------------

      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % IceNO3                            &
     &                      )
!
        CALL IceNO3bc_tile (ng, tile,                                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          liold(ng), linew(ng),                   &
     &                          ICE(ng)%ui,                             &
     &                          ICE(ng)%vi,                             &
     &                          ICE(ng)%IceNO3)
!
# endif
!
# ifdef MELT_PONDS
! ---------------------------------------------------------------------
!  Advect the surface melt water.
! ---------------------------------------------------------------------
      DO j=LBj,UBj
        DO i=LBi,UBi
          wrk(i,j,liold(ng)) = ICE(ng)%ai(i,j,liold(ng))*                         &
       &                       ICE(ng)%apond(i,j,liold(ng))
        ENDDO
      ENDDO

      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
#  ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
#  endif
#  ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#  endif
#  ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#  endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      wrk                                         &
     &                      )

      DO j=LBj,UBj
        DO i=LBi,UBi
          IF (ICE(ng)%ai(i,j,linew(ng)) > 0.0_r8) THEN
            ICE(ng)%apond(i,j,linew) = wrk(i,j,linew(ng))/              &
       &             ICE(ng)%ai(i,j,linew(ng))
            ICE(ng)%apond(i,j,linew) =                                  &
       &             MAX(0.0_r8,ICE(ng)%apond(i,j,linew))
          ELSE
            ICE(ng)%apond(i,j,linew) = 0.0_r8
          END IF
        ENDDO
      ENDDO
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          wrk(i,j,liold(ng)) = ICE(ng)%ai(i,j,liold(ng))*               &
       &                       ICE(ng)%apond(i,j,liold(ng))*            &
       &                       ICE(ng)%hpond(i,j,liold(ng))
        ENDDO
      ENDDO

      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      wrk                                         &
     &                      )

      DO j=LBj,UBj
        DO i=LBi,UBi
          IF (ICE(ng)%ai(i,j,linew(ng))*ICE(ng)%apond(i,j,linew(ng))    &
     &                    > 0.0_r8) THEN
            ICE(ng)%hpond(i,j,linew) = wrk(i,j,linew(ng))/              &
     &          (ICE(ng)%ai(i,j,linew(ng))*ICE(ng)%apond(i,j,linew(ng)))
            ICE(ng)%hpond(i,j,linew) =                                  &
       &             MAX(0.0_r8,ICE(ng)%hpond(i,j,linew))
          ELSE
            ICE(ng)%hpond(i,j,linew) = 0.0_r8
          END IF
        ENDDO
      ENDDO
!
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
     &                  BOUNDARY(ng)%apond_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%apond_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%apond_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%apond_south(LBi:UBi),              &
     &                  ICE(ng)%ui,                                     &
     &                  ICE(ng)%vi,                                     &
     &                  ICE(ng)%apond, LBC(:,isApond,ng))
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
     &                  BOUNDARY(ng)%hpond_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%hpond_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%hpond_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%hpond_south(LBi:UBi),              &
     &                  ICE(ng)%ui,                                     &
     &                  ICE(ng)%vi,                                     &
     &                  ICE(ng)%hpond, LBC(:,isHpond,ng))
# endif
!
! ---------------------------------------------------------------------
!  Advect the interior ice temperature.
! ---------------------------------------------------------------------
!
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % enthalpi                          &
     &                      )
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ICE(ng)%ti(i,j,linew(ng)) = ICE(ng)%enthalpi(i,j,linew(ng))/  &
       &                  MAX(ICE(ng)%hi(i,j,linew(ng)),1.0E-6_r8)
          IF (ICE(ng)%hi(i,j,linew(ng)).LE.min_h(ng)) THEN
            ICE(ng)%enthalpi(i,j,linew(ng)) = 0.0_r8
            ICE(ng)%ti(i,j,linew(ng)) = 0.0_r8
          END IF
        ENDDO
      ENDDO
!
      CALL tibc_tile (ng, tile, iNLM,                                   &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold(ng), linew(ng),                             &
     &                ICE(ng)%ui,                                       &
     &                ICE(ng)%vi,                                       &
     &                ICE(ng)%hi,                                       &
     &                ICE(ng)%ti,                                       &
     &                ICE(ng)%enthalpi)
!
! ---------------------------------------------------------------------
!  Advect the ice age.
! ---------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ICE(ng)%hage(i,j,liold(ng)) = ICE(ng)%hi(i,j,liold(ng))*      &
     &                       ICE(ng)%ageice(i,j,liold(ng))
          ICE(ng)%hage(i,j,linew(ng)) = ICE(ng)%hi(i,j,linew(ng))*      &
     &                       ICE(ng)%ageice(i,j,linew(ng))
          IF (ICE(ng)%hi(i,j,liold(ng)).LE.min_h(ng)) THEN
            ICE(ng)%hage(i,j,liold(ng)) = 0.0_r8
            ICE(ng)%ageice(i,j,liold(ng)) = 0.0_r8
          END IF
        ENDDO
      ENDDO
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), linew(ng), liold(ng), liunw(ng),  &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % hage                              &
     &                      )
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ICE(ng)%ageice(i,j,linew(ng)) = ICE(ng)%hage(i,j,linew(ng))/  &
     &                  MAX(ICE(ng)%hi(i,j,linew(ng)),1.0E-6_r8)
          IF (ICE(ng)%hi(i,j,linew(ng)).LE.min_h(ng)) THEN
            ICE(ng)%hage(i,j,linew(ng)) = 0.0_r8
            ICE(ng)%ageice(i,j,linew(ng)) = 0.0_r8
          END IF
        ENDDO
      ENDDO
!
!        CALL i2d_bc_tile (ng, tile, iNLM,                               &
!     &                    LBi, UBi, LBj, UBj,                           &
!     &                    IminS, ImaxS, JminS, JmaxS,                   &
!     &                    liold(ng), linew(ng),                         &
!     &                    BOUNDARY(ng)%ageice_west(LBj:UBj),            &
!     &                    BOUNDARY(ng)%ageice_east(LBj:UBj),            &
!     &                    BOUNDARY(ng)%ageice_north(LBi:UBi),           &
!     &                    BOUNDARY(ng)%ageice_south(LBi:UBi),           &
!     &                    ICE(ng)%ui,                                   &
!     &                    ICE(ng)%vi,                                   &
!     &                    ICE(ng)%ageice, LBC(:,isHice,ng))
!      CALL ageicebc_tile (ng, tile,                                     &
!     &                          LBi, UBi, LBj, UBj,                     &
!     &                          liold(ng), linew(ng),                   &
!     &                          ICE(ng)%ui,                             &
!     &                          ICE(ng)%vi,                             &
!     &                          ICE(ng)%hi,                             &
!     &                          ICE(ng)%ageice,                         &
!     &                          ICE(ng)%hage)
!!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%ai(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%hi(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%hsn(:,:,linew(ng)))
# ifdef MELT_PONDS
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%apond(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%hpond(:,:,linew(ng)))
# endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%ti(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%enthalpi(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%ageice(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%hage(:,:,linew(ng)))
# if defined ICE_BIO && defined BERING_10K
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%IcePhL(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%IceNO3(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%IceNH4(:,:,linew(ng)))
# endif
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%ai(:,:,linew(ng)),                    &
     &                    ICE(ng)%hi(:,:,linew(ng)),                    &
     &                    ICE(ng)%hsn(:,:,linew(ng)))
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%ti(:,:,linew(ng)),                    &
     &                    ICE(ng)%enthalpi(:,:,linew(ng)),              &
     &                    ICE(ng)%ageice(:,:,linew(ng)),                &
     &                    ICE(ng)%hage(:,:,linew(ng)))
#  ifdef MELT_PONDS
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%apond(:,:,linew(ng)),                 &
     &                    ICE(ng)%hpond(:,:,linew(ng)))
#  endif
#  if defined ICE_BIO && defined BERING_10K
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ICE(ng)%IcePhL(:,:,linew(ng)),                &
     &                    ICE(ng)%IceNO3(:,:,linew(ng)),                &
     &                    ICE(ng)%IceNH4(:,:,linew(ng)))
#  endif
# endif
#endif
      RETURN
      END SUBROUTINE ice_advect_all_tile
!
!==========================================================================!
      SUBROUTINE ice_advect_tile (ng, tile,                             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs, linew, liold, liunw,                &
#ifdef MASKING
     &                        rmask,                                    &
#endif
#ifdef WET_DRY
     &                        rmask_wet,                                &
#endif
#ifdef ICESHELF
     &                        zice,                                     &
#endif
#ifndef ICE_UPWIND
     &                        pm, pn,                                   &
#endif
     &                        on_u, om_v, omn,                          &
     &                        ui, vi, scr)
!==========================================================================!

      USE mod_param
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, linew, liold, liunw

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
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
      real(r8), intent(inout) :: scr(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
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
      real(r8), intent(inout) :: scr(LBi:UBi,LBj:UBj,2)
#endif

!
! Local variable definitions
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aflxu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aflxv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE

      real(r8), parameter :: epsil = 1.0E-15_r8
      real(r8), parameter :: add = 3.0E+3_r8

      real(r8) :: Cu_crss, Cu
      real(r8) :: rateu
      real(r8) :: ratev
      real(r8) :: rateyiu
      real(r8) :: ratexiv
      real(r8) :: uspeed
      real(r8) :: vspeed

#include "set_bounds.h"

# ifndef ICE_UPWIND
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
# else
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
# endif
!
! upstream:
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax+1
          aflxu(i,j)=on_u(i,j)*                                         &
     &          (max(0.0_r8,ui(i,j,liunw))*scr(i-1,j,liold)             &
     &          +min(0.0_r8,ui(i,j,liunw))*scr(i,j,liold))
        END DO
      END DO
      DO j=Jmin,Jmax+1
        DO i=Imin,Imax
          aflxv(i,j)=om_v(i,j)*                                         &
     &          (max(0.0_r8,vi(i,j,liunw))*scr(i,j-1,liold)             &
     &          +min(0.0_r8,vi(i,j,liunw))*scr(i,j,liold))
!
        END DO
      END DO
!
! step number 1 in mpdata:
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
!
          ar(i,j)=1.0_r8/omn(i,j)
          aif(i,j)=(scr(i,j,liold)-dtice(ng)*(aflxu(i+1,j)-aflxu(i,j)   &
     &        +aflxv(i,j+1)-aflxv(i,j))*ar(i,j))
#ifdef MASKING
          aif(i,j) = aif(i,j)*rmask(i,j)
#endif
#ifdef WET_DRY
          aif(i,j) = aif(i,j)*rmask_wet(i,j)
#endif
#ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) aif(i,j) = 0.0_r8
#endif
        END DO
      END DO
!
! set values at the open boundaries
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jmin,Jmax
            aif(Istr-1,j)=aif(Istr,j)   !? scr(Istr-1,j,liold)
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jmin,Jmax
            aif(Iend+1,j)=aif(Iend,j)  !? scr(Iend+1,j,liold)
          END DO
        END IF
      END IF
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Imin,Imax
             aif(i,Jstr-1)=aif(i,Jstr)  !??? scr(i,Jstr-1,liold)
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Imin,Imax
             aif(i,Jend+1)=aif(i,Jend)  !??? scr(i,Jend+1,liold)
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
!
! mask ???
!
#ifdef MASKING
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          aif(i,j)=aif(i,j)*rmask(i,j)
        END DO
      END DO
#endif
#ifdef WET_DRY
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          aif(i,j)=aif(i,j)*rmask_wet(i,j)
        END DO
      END DO
#endif
#ifdef ICESHELF
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          IF (zice(i,j).ne.0.0_r8) THEN
            aif(i,j) = 0.0_r8
          END IF
        END DO
      END DO
#endif

#ifndef ICE_UPWIND
!
! Antidiffusive corrector step:
!-------------- --------- -----
! This is needed to avoid touching "aif" under land mask.
! Note that only aif(i,j) and aif(i-1,j) are allowed to appear 
! explicitly in the code segment below. This is OK 
! because if either of them masked, then "ui" is zero 
! at that point, and therefore no additional masking is required. 
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
          FE(i,j)=0.5*                                                  &
#  ifdef MASKING
     &                          vmask(i,j)*                             &
#  endif
#  ifdef WET_DRY
     &                          vmask_wet(i,j)*                         &
#  endif
     &                    (aif(i,j)-aif(i,j-1))
          FX(i,j)=0.5*                                                  &
#  ifdef MASKING
     &                          umask(i,j)*                             &
#  endif
#  ifdef WET_DRY
     &                          umask_wet(i,j)*                         &
#  endif
     &                    (aif(i,j)-aif(i-1,j))
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend+1
          rateu=(aif(i,j)-aif(i-1,j))/max(epsil, aif(i,j)+aif(i-1,j))

          rateyiu=(FE(i,j+1)+FE(i,j)  +FE(i-1,j+1)+FE(i-1,j))           &
     &        /( max( epsil,  aif(i  ,j)+FE(i  ,j+1)-FE(i  ,j)          &
     &                       +aif(i-1,j)+FE(i-1,j+1)-FE(i-1,j)          &
     &                                                      ))

          Cu=0.5*dtice(ng)*(pm(i,j)+pm(i-1,j))*ui(i,j,liunw)

          Cu_crss=0.5*dtice(ng) * 0.0625*( pn(i-1,j+1)+pn(i,j+1)        &
     &                                 +pn(i-1,j-1)+pn(i,j-1)           &
     &                  )*( vi(i-1,j+1,liunw)+vi(i,j+1,liunw)           &
     &                       +vi(i-1,j,liunw)  +vi(i,j,liunw)           &
     &                                                     )

          uspeed=rateu*(abs(ui(i,j,liunw)) -Cu*ui(i,j,liunw))           &
     &                       -rateyiu*Cu_crss * ui(i,j,liunw)

          aflxu(i,j)=on_u(i,j)*( max(0.,uspeed)*aif(i-1,j)              &
     &                        +min(0.,uspeed)*aif(i,j) )
        END DO
      END DO

      DO j=Jstr,Jend+1
        DO i=Istr,Iend
          ratev=(aif(i,j)-aif(i,j-1))/max(epsil, aif(i,j)+aif(i,j-1))

          ratexiv=(FX(i+1,j)+FX(i,j)  +FX(i+1,j-1)+FX(i,j-1))           &
     &        /( max( epsil,  aif(i,j  )+FX(i+1,j  )-FX(i,j  )          &
     &                       +aif(i,j-1)+FX(i+1,j-1)-FX(i,j-1)          &
     &                                                      ))

          Cu=0.5*dtice(ng)*(pn(i,j)+pn(i,j-1))*vi(i,j,liunw)

          Cu_crss=0.5*dtice(ng) * 0.0625*( pm(i+1,j)+pm(i+1,j-1)        &
     &                                 +pm(i-1,j)+pm(i-1,j-1)           &
     &                 )*(  ui(i,j,liunw)    +ui(i+1,j,liunw)           &
     &                     +ui(i,j-1,liunw)+ui(i+1,j-1,liunw)           &
     &                                                     )

          vspeed=ratev*(abs(vi(i,j,liunw)) -Cu*vi(i,j,liunw))           &
     &                       -ratexiv*Cu_crss * vi(i,j,liunw)

          aflxv(i,j)=om_v(i,j)*( max(0.,vspeed)*aif(i,j-1)              &
     &                        +min(0.,vspeed)*aif(i,j) )
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
          aif(i,j)=aif(i,j) -dtice(ng)*pm(i,j)*pn(i,j)*                 &
     &             (aflxu(i+1,j)-aflxu(i,j) + aflxv(i,j+1)-aflxv(i,j))
#  ifdef MASKING
          aif(i,j)=aif(i,j)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          aif(i,j)=aif(i,j)*rmask_wet(i,j)
#  endif
#  ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) aif(i,j)=0.
#  endif
        END DO
      END DO
# endif /* !ICE_UPWIND */

      DO j=Jstr,Jend
        DO i=Istr,Iend
            scr(i,j,linew) = aif(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ice_advect_tile
