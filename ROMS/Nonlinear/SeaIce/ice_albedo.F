#include "cppdefs.h"
      MODULE ice_albedo_mod

#if defined ICE_MODEL && defined ICE_ALBEDO
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module computes the surface albedo over seawater, snow, or ice !
!  using selected formulation.                                         !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Briegleb, B.P., P. Minnis, V. Ramanathan, and E. Harrison,        !
!      1986: Comparison of Regional Clear-Sky Albedos Inferred         !
!      from Satellite Observations and Model Computations, J.          !
!      Climate and Appled Meteor., 25, 214-226.                        !
!                                                                      !
!    Ebert, E.E. and J.A. Curry, 1992: A parameterization of ice       !
!      cloud optical properties for climate models, J. Geophys.        !
!      Res., 97, 3831-3836, doi: doi.org/10.1029/91jD02472.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ice
      USE mod_scalars
!
      USE dateclock_mod,   ONLY : caldate
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PUBLIC  :: ice_albedo
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_albedo (ng, tile, model)
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
      CALL ice_albedo_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold(ng), linew(ng),                       &
# if defined SHORTWAVE && (defined ALBEDO_CURVE || defined ALBEDO_SZO)
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
# endif
     &                      ICE(ng) % Fi,                               &
     &                      ICE(ng) % Si,                               &
     &                      FORCES(ng) % albedo_ice,                    &
     &                      FORCES(ng) % albedo)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_albedo
!
!***********************************************************************
      SUBROUTINE ice_albedo_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            liold, linew,                         &
# if defined SHORTWAVE && (defined ALBEDO_CURVE || defined ALBEDO_SZO)
     &                            lonr, latr,                           &
# endif
     &                            Fi, Si,                               &
     &                            albedo_ice, albedo)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, linew
!
# ifdef ASSUMED_SHAPE
#  if defined SHORTWAVE && (defined ALBEDO_CURVE || defined ALBEDO_SZO)
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Fi(LBi:,LBj:,:)
      real(r8), intent(in) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(out) :: albedo_ice(LBi:,LBj:)
      real(r8), intent(out) :: albedo(LBi:,LBj:)

# else

#  if defined SHORTWAVE && (defined ALBEDO_CURVE || defined ALBEDO_SZO)
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(in) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
      real(r8), intent(out) :: albedo_ice(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: albedo(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, li_stp
# if defined ALBEDO_SZO
      integer :: iday, month, year
# endif
!
      real(r8), parameter :: alb_w = 0.06_r8

# ifdef ALBEDO_CSIM
      real(r8), parameter :: alb_i_thick = 0.54_r8
      real(r8), parameter :: alb_s_dry = 0.83_r8
      real(r8), parameter :: alb_s_wet = 0.70_r8
# else
      real(r8), parameter :: alb_i_dry = 0.65_r8
      real(r8), parameter :: alb_i_wet = 0.60_r8
      real(r8), parameter :: alb_s_dry = 0.85_r8
      real(r8), parameter :: alb_s_wet = 0.72_r8
# endif

# ifdef ICE_ALB_EC92
      real(r8) :: alb_ice, alb_snow, Hice, Hsnow
# endif
# if defined ALBEDO_SZO
      real(r8) :: Dangle, Hangle, LatRad, zenith
      real(dp) :: hour, yday
# endif
      real(r8) :: cff, cff1, cff2, sfc_temp
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick

# include "set_bounds.h"
!
!=======================================================================
!  Compute the surface albedo over seawater, snow or ice
!=======================================================================
!
!  Set ice model time level.
!
      IF (PerfectRST(ng).and.(iic(ng).eq.ntstart(ng))) THEN
        li_stp=liold
      ELSE
        li_stp=linew
      END IF

# if defined ALBEDO_SZO
!
!  Calculate the solar zenith angle and approximate the  ocean albedo
!  using Briegleb et al. (1986)  empirical formulation.
!
      CALL caldate (tdays(ng), h_dp=hour, yd_dp=yday)
!
!  Estimate solar declination angle (radians).
!
      Dangle=23.44_r8*COS((172.0_r8-yday)*2.0_r8*pi/365.25_r8)
      Dangle=Dangle*deg2rad
!
!  Compute hour angle (radians).
!
      Hangle=(12.0_r8-hour)*pi/12.0_r8
# endif
!
!-----------------------------------------------------------------------
!  Compute surface albedo.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
!
!  Compute the ice/snow albedo
!
          cff=1.0_r8/(Si(i,j,li_stp,isAice)+0.001_r8)
          ice_thick (i,j)=cff*Si(i,j,li_stp,isHice)
          snow_thick(i,j)=cff*Si(i,j,li_stp,isHsno)
          sfc_temp=Fi(i,j,icIsst)

# ifdef ICE_ALB_EC92
!
!  Ice and snow albedo is calculated from Ebert and Curry (1992).
!
          alb_ice=0.0_r8
          alb_snow=0.0_r8

          IF (Si(i,j,li_stp,isAice).gt.min_ai(ng)) THEN
            Hice =ice_thick(i,j)
            Hsnow=snow_thick(i,j)
!
!  The threshold test can still lead to a negative albedo, so we modify
!  the minimum ice thickness for this estimate to something more
!  reasonable (SMD).
!
            Hice=MAX(Hice, 0.01_r8)
            IF (Hice.ge.2.0_r8) THEN
              alb_ice=0.561632_r8
            ELSE IF (Hice.ge.1.0_r8) THEN
              alb_ice=0.07616_r8*Hice+0.414492_r8
            ELSE
              alb_ice=0.082409_r8*LOG(Hice)+0.485472_r8
            END IF
!
!  Approximated values for alb_snow depends on COSZ, but small
!  variation.
!
            IF (Si(i,j,li_stp,isHmel).gt.0.0_r8) THEN
              IF (Hsnow.ge.0.1_r8) THEN
                alb_snow=0.701009_r8
              ELSE
                alb_snow=alb_ice+(0.701009_r8-alb_ice)*Hsnow/0.1_r8
              END IF
            ELSE
              alb_snow=0.83_r8
            END IF
!
!  Ebert and Curry estimate the melt pond albedo to be between 0.14 and
!  0.26, depending on the thickness of the meltwater ponds.  But of
!  course, the meltwater ponds do not cover all the ice surface, except
!  maybe effectively when the ice consists of small pieces with little
!  freeboard. So let's try the idea that the effective albedo of the
!  ice/melt pond/snow is that of ice for bare ice, that of snow
!  (potentially modified by meltwater if snow is present), and the
!  average of the high-end melt pond (0.26) and the low-end ice bare
!  ice (0.56) when melt ponds of more than a few centimeters deep are
!  present (SMD).
!
            IF (Si(i,j,li_stp,isHsno).gt.0.0_r8) THEN
              albedo_ice(i,j)=alb_snow
            ELSE IF (Si(i,j,li_stp,isHmel).gt.0.02_r8) THEN
              albedo_ice(i,j)=0.42_r8
            ELSE
              albedo_ice(i,j)=alb_ice
            END IF
          ELSE
            albedo_ice(i,j)=alb_w             ! water albedo
          ENDIF

# elif defined ALBEDO_CSIM
!
!  Community Sea Ice Model (CSIM) formulation.
!
          fhi=MIN(ATAN(4.0_r8*ice_thick(i,j))/ATAN(2.0_r8), 1.0_r8)
          fsn=snow_thick(i,j)/(snow_thick(i,j)+0.02_r8)
          alb_i_dry=alb_w*(1-fh)+alb_i_thick*fh
          cff1=alb_s_wet-alb_s_dry
          cff2=-0.075                            ! alb_i_wet - alb_i_dry
          IF (Si(i,j,li_stp,isAice).gt.min_ai(ng)) THEN
            IF (sfc_temp-273.16_r8.gt. -1.0_r8) THEN
              alb_snow=cff1*(sfc_temp-272.16_r8)+alb_s_dry
            ELSE
              alb_snow=alb_s_dry
            ENDIF
            IF (sfc_temp-273.16_r8.gt. -1._r8) THEN
              alb_ice=cff2*(sfc_temp-272.16_r8)+alb_i_dry
            ELSE
              alb_ice=alb_i_dry
            END IF
            albedo_ice(i,j)=fsn*alb_snow+(1-fsn)*alb_ice
          ELSE
            albedo_ice(i,j)=alb_w             ! water albedo
          END IF
# else
!
!  Default formulation.
!
          cff1=alb_s_wet-alb_s_dry
          cff2=alb_i_wet-alb_i_dry
          IF (Si(i,j,li_stp,isAice).gt.min_ai(ng)) THEN
            IF (Si(i,j,li_stp,isHsno).gt.0.0_r8) THEN
              IF ((sfc_temp-273.16_r8).gt.-1.0_r8) THEN
                albedo_ice(i,j)=cff1*(sfc_temp-272.16_r8)+alb_s_dry
              ELSE
                albedo_ice(i,j)=alb_s_dry
              END IF
            ELSE
              IF ((sfc_temp-273.16_r8).gt.-1.0_r8) THEN
                albedo_ice(i,j)=cff2*(sfc_temp-272.16_r8)+alb_i_dry
              ELSE
                albedo_ice(i,j)=alb_i_dry
              END IF
            END IF
          ELSE
            albedo_ice(i,j)=alb_w             ! water albedo
          END IF
# endif

# ifdef ALBEDO_CURVE
!
!  Compute seawater albedo from curve.
!
          albedo(i,j)=0.069_r8-                                         &
     &                0.011_r8*COS(2.0_r8*deg2rad*latr(i,j))

# elif defined ALBEDO_SZO
!
!  Compute albedo over water based on the calendar day and time. It
!  ASSUMES that the forcing file contains the TOTAL INCIDENT SHORTWAVE
!  RADIATION (SMD).
!
!  Local daylight is a function of the declination (Dangle) and hour
!  angle adjusted for the local meridian (Hangle-lonr(i,j)/15.0).
!  The 15.0 factor is because the sun moves 15 degrees every hour.
!
          LatRad=latr(i,j)*deg2rad
          cff1=SIN(LatRad)*SIN(Dangle)
          cff2=COS(LatRad)*COS(Dangle)
          zenith=max(cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad             &
     &                             -pi/12.0_r8),0.0_r8)
!
!  Use albedo formula from Briegleb et al. (1986).
!
          albedo(i,j)=0.026_r8/(zenith**1.7_r8+0.065_r8)+               &
     &                0.15_r8*(zenith-0.1_r8)*                          &
     &                (zenith-0.5_r8)*(zenith-1.0_r8)
# else
!
!  Default constant seawater albedo.
!
          albedo(i,j)=alb_w
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo_ice)

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo)
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    albedo, albedo_ice)
# endif
!
      RETURN
      END SUBROUTINE ice_albedo_tile
!
#endif
      END module ice_albedo_mod
