       SUBROUTINE ice_advect (ng, tile)
!
!*************************************************** W. Paul Budgell ***
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
      CALL wclock_on (ng, iNLM, 49)
#endif
! ---------------------------------------------------------------------
!  Advect the ice concentration.
! ---------------------------------------------------------------------
      CALL ice_advect_all_tile (ng, tile,                               &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS                  &
     &                      )
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 49)
#endif
      RETURN
      END SUBROUTINE ice_advect

      SUBROUTINE ice_advect_all_tile (ng, tile,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS                  &
     &                      )

      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_ice
      USE mod_forces
      USE mod_scalars
      USE mod_stepping
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
      USE aibc_mod, ONLY : aibc_tile
      USE hibc_mod, ONLY : hibc_tile
      USE hsnbc_mod, ONLY : hsnbc_tile
      USE tibc_mod, ONLY : tibc_tile
      USE sfwatbc_mod, ONLY : sfwatbc_tile
      USE ageicebc_mod, ONLY : ageicebc_tile
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
        CALL aibc_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          liold(ng), linew(ng),                   &
     &                          ICE(ng)%ui,                             &
     &                          ICE(ng)%vi,                             &
     &                          ICE(ng)%ai)
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
        CALL hibc_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          liold(ng), linew(ng),                   &
     &                          ICE(ng)%ui,                             &
     &                          ICE(ng)%vi,                             &
     &                          ICE(ng)%hi)
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
        CALL hsnbc_tile (ng, tile,                                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          liold(ng), linew(ng),                   &
     &                          ICE(ng)%ui,                             &
     &                          ICE(ng)%vi,                             &
     &                          ICE(ng)%hsn)
!
! ---------------------------------------------------------------------
!  Advect the surface melt water.
! ---------------------------------------------------------------------
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
     &                      ICE(ng) % sfwat                             &
     &                      )
!
        CALL sfwatbc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     liold(ng), linew(ng),                        &
     &                     ICE(ng)%ui,                                  &
     &                     ICE(ng)%vi,                                  &
     &                     ICE(ng)%sfwat)
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
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ICE(ng)%ti(i,j,linew(ng)) = ICE(ng)%enthalpi(i,j,linew(ng))/  &
       &                  MAX(ICE(ng)%hi(i,j,linew(ng)),1.0E-6_r8)
          IF (ICE(ng)%hi(i,j,linew(ng)).LE.min_h(ng)) THEN
            ICE(ng)%enthalpi(i,j,linew(ng)) = 0.0_r8
            ICE(ng)%ti(i,j,linew(ng)) = 0.0_r8
          END IF
        ENDDO
      ENDDO
!
        CALL tibc_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          liold(ng), linew(ng), min_h(ng),        &
     &                          ICE(ng)%ui,                             &
     &                          ICE(ng)%vi,                             &
     &                          ICE(ng)%hi,                             &
     &                          ICE(ng)%ti,                             &
     &                          ICE(ng)%enthalpi)
!
! ---------------------------------------------------------------------
!  Advect the ice age.
! ---------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
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
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ICE(ng)%ageice(i,j,linew(ng)) = ICE(ng)%hage(i,j,linew(ng))/  &
     &                  MAX(ICE(ng)%hi(i,j,linew(ng)),1.0E-6_r8)
          IF (ICE(ng)%hi(i,j,linew(ng)).LE.min_h(ng)) THEN
            ICE(ng)%hage(i,j,linew(ng)) = 0.0_r8
            ICE(ng)%ageice(i,j,linew(ng)) = 0.0_r8
          END IF
        ENDDO
      ENDDO
!
        CALL ageicebc_tile (ng, tile,                                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          liold(ng), linew(ng), min_h(ng),        &
     &                          ICE(ng)%ui,                             &
     &                          ICE(ng)%vi,                             &
     &                          ICE(ng)%hi,                             &
     &                          ICE(ng)%ageice,                         &
     &                          ICE(ng)%hage)
!
# if defined EW_PERIODIC || defined NS_PERIODIC
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%ai(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%hi(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%hsn(:,:,linew(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ICE(ng)%sfwat(:,:,linew(ng)))
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
# endif
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      ICE(ng)%ai(:,:,linew(ng)),                  &
     &                      ICE(ng)%hi(:,:,linew(ng)),                  &
     &                      ICE(ng)%hsn(:,:,linew(ng)),                 &
     &                      ICE(ng)%sfwat(:,:,linew(ng)))
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      ICE(ng)%ti(:,:,linew(ng)),                  &
     &                      ICE(ng)%enthalpi(:,:,linew(ng)),            &
     &                      ICE(ng)%ageice(:,:,linew(ng)),              &
     &                      ICE(ng)%hage(:,:,linew(ng)))
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

      real(r8) :: fakt1
      real(r8) :: fakt2

      real(r8) :: aim1
      real(r8) :: ajm1
      real(r8) :: aim1jm1u
      real(r8) :: aim1jm1v
      real(r8) :: ajp1
      real(r8) :: aim1jp1
      real(r8) :: aip1
      real(r8) :: aip1jm1

      real(r8) :: Cu_crss, Cu
      real(r8) :: rateu
      real(r8) :: ratev
      real(r8) :: rateyiu
      real(r8) :: ratexiv
      real(r8) :: ratiou
      real(r8) :: ratiov
      real(r8) :: ratioyiu
      real(r8) :: ratioxiv
      real(r8) :: uiv
      real(r8) :: viu
      real(r8) :: uspeed
      real(r8) :: vspeed

#include "set_bounds.h"

# ifndef ICE_UPWIND
#  ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN
        Imin=Istr
      ELSE
        Imin=Istr-1
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=Iend
      ELSE
        Imax=Iend+1
      END IF
#  else
      Imin=Istr-1
      Imax=Iend+1
#  endif
#  ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
        Jmin=Jstr
      ELSE
        Jmin=Jstr-1
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=Jend
      ELSE
        Jmax=Jend+1
      END IF
#  else
      Jmin=Jstr-1
      Jmax=Jend+1
#  endif
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
# ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN
        DO j=Jmin,Jmax
          aif(Istr-1,j)=aif(Istr,j)   !? scr(Istr-1,j,liold)
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=Jmin,Jmax
          aif(Iend+1,j)=aif(Iend,j)  !? scr(Iend+1,j,liold)
        END DO
      END IF
# endif
# ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
        DO i=Imin,Imax
             aif(i,Jstr-1)=aif(i,Jstr)  !??? scr(i,Jstr-1,liold)
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=Imin,Imax
             aif(i,Jend+1)=aif(i,Jend)  !??? scr(i,Jend+1,liold)
        END DO
      END IF
#  ifndef EW_PERIODIC
      IF (WESTERN_EDGE .and. SOUTHERN_EDGE) THEN
        aif(Istr-1,Jstr-1)=aif(Istr,Jstr)
      END IF
      IF (WESTERN_EDGE .and. NORTHERN_EDGE) THEN
        aif(Istr-1,Jend+1)=aif(Istr,Jend)
      END IF
      IF (EASTERN_EDGE .and. SOUTHERN_EDGE) THEN
        aif(Iend+1,Jstr-1)=aif(Iend,Jstr)
      END IF
      IF (EASTERN_EDGE .and. NORTHERN_EDGE) THEN
        aif(Iend+1,Jend+1)=aif(Iend,Jend)
      END IF
#  endif
# endif
!
! mask ???
!
#ifdef FOO
#ifdef MASKING
      do j=J_RANGE
      do i=I_RANGE
        aif(i,j)=aif(i,j)*rmask(i,j)
      enddo
      enddo
#endif
#ifdef WET_DRY
      do j=J_RANGE
      do i=I_RANGE
        aif(i,j)=aif(i,j)*rmask_wet(i,j)
      enddo
      enddo
#endif
#ifdef ICESHELF
      do j=J_RANGE
      do i=I_RANGE
         IF (zice(i,j).ne.0.0_r8) THEN
            aif(i,j) = 0.0_r8
         ENDIF
      enddo
      enddo
#endif
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

          Cu=0.5*dtice*(pm(i,j)+pm(i-1,j))*ui(i,j,liunw)

          Cu_crss=0.5*dtice * 0.0625*( pn(i-1,j+1)+pn(i,j+1)            &
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

          Cu=0.5*dtice*(pn(i,j)+pn(i,j-1))*vi(i,j,liunw)

          Cu_crss=0.5*dtice * 0.0625*( pm(i+1,j)+pm(i+1,j-1)            &
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
          aif(i,j)=aif(i,j) -dtice*iAr(i,j)*( aflxu(i+1,j)-aflxu(i,j)   &
     &                                       +aflxv(i,j+1)-aflxv(i,j))
#  ifdef MASKING
          aif(i,j)=aif(i,j)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          aif(i,j)=aif(i,j)*rmask_wet(i,j)
#  endif
#  ifdef ICESHELF
          IF (zice(i,j).ne.0.) aif(i,j)=0.
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
