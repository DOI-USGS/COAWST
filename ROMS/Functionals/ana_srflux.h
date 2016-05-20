      SUBROUTINE ana_srflux (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets kinematic surface solar shortwave radiation    !
!  flux "srflx" (degC m/s) using an analytical expression.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_srflux_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
#ifdef ALBEDO_CLOUD
     &                      FORCES(ng) % cloud,                         &
     &                      FORCES(ng) % Hair,                          &
     &                      FORCES(ng) % Tair,                          &
     &                      FORCES(ng) % Pair,                          &
#endif
     &                      FORCES(ng) % srflx)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(27)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_srflux
!
!***********************************************************************
      SUBROUTINE ana_srflux_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            lonr, latr,                           &
#ifdef ALBEDO_CLOUD
     &                            cloud, Hair, Tair, Pair,              &
#endif
     &                            srflx)
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
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# ifdef ALBEDO_CLOUD
      real(r8), intent(in) :: cloud(LBi:,LBj:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
# endif
      real(r8), intent(out) :: srflx(LBi:,LBj:)
#else
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# ifdef ALBEDO_CLOUD
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: srflx(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
#if defined ALBEDO_CLOUD || defined DIURNAL_SRFLUX
      integer :: iday, month, year
      real(r8) :: Dangle, Hangle, LatRad
      real(r8) :: cff1, cff2, hour, yday
# ifdef ALBEDO_CLOUD
      real(r8) :: Rsolar, e_sat, vap_p, zenith
# endif
#endif
      real(r8) :: cff

#include "set_bounds.h"

#if defined ALBEDO_CLOUD || defined DIURNAL_SRFLUX
!
!-----------------------------------------------------------------------
!  Compute shortwave radiation (degC m/s):
!
!  ALBEDO_CLOUD option: Compute shortwave radiation flux using the Laevastu
!                 cloud correction to the Zillman equation for cloudless
!  radiation (Parkinson and Washington 1979, JGR, 84, 311-337).  Notice
!  that flux is scaled from W/m2 to degC m/s by dividing by (rho0*Cp).
!
!  DIURNAL_SRFLUX option: Modulate shortwave radiation SRFLX (which
!                         read and interpolated elsewhere) by the local
!  diurnal cycle (a function of longitude, latitude and day-of-year).
!  This option is provided for cases where SRFLX computed by SET_DATA is
!  an average over >= 24 hours. For "diurnal_srflux" to work ana_srflux
!  must be undefined. If you want a strictly analytical diurnal cycle
!  enter it explicitly at the end of this subroutine or use the "albedo"
!  option.
!
!  For a review of shortwave radiation formulations check:
!
!    Niemela, S., P. Raisanen, and H. Savijarvi, 2001: Comparison of
!      surface radiative flux parameterizations, Part II, Shortwave
!      radiation, Atmos. Res., 58, 141-154.
!
!-----------------------------------------------------------------------
!
!  Assume time is in modified Julian day.  Get hour and year day.
!
      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
!
!  Estimate solar declination angle (radians).
!
      Dangle=23.44_r8*COS((172.0_r8-yday)*2.0_r8*pi/365.25_r8)
      Dangle=Dangle*deg2rad
!
!  Compute hour angle (radians).
!
      Hangle=(12.0_r8-hour)*pi/12.0_r8
!
# ifdef ALBEDO_CLOUD
      Rsolar=Csolar/(rho0*Cp)
# endif
      DO j=JstrT,JendT
        DO i=IstrT,IendT
!
!  Local daylight is a function of the declination (Dangle) and hour
!  angle adjusted for the local meridian (Hangle-lonr(i,j)).
!
          LatRad=latr(i,j)*deg2rad
          cff1=SIN(LatRad)*SIN(Dangle)
          cff2=COS(LatRad)*COS(Dangle)
# if defined ALBEDO_CLOUD
!
!  Estimate variation in optical thickness of the atmosphere over
!  the course of a day under cloudless skies (Zillman, 1972). To
!  obtain net incoming shortwave radiation multiply by (1.0-0.6*c**3),
!  where c is the fractional cloud cover.
!
!  The equation for saturation vapor pressure is from Gill (Atmosphere-
!  Ocean Dynamics, pp 606).
!
          srflx(i,j)=0.0_r8
          zenith=cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad)
          IF (zenith.gt.0.0_r8) THEN
            cff=(0.7859_r8+0.03477_r8*Tair(i,j))/                       &
     &          (1.0_r8+0.00412_r8*Tair(i,j))
            e_sat=10.0_r8**cff    ! saturation vapor pressure (hPa=mbar)
#  ifdef SPECIFIC_HUMIDITY
!  With this directive specific humidity is input as kg/kg
            vap_p=Pair(i,j)*Hair(i,j)/(0.62197_r8+0.378_r8*Hair(i,j))
#  else
            vap_p=e_sat*Hair(i,j) ! water vapor pressure (hPa=mbar)
#  endif
            srflx(i,j)=Rsolar*zenith*zenith*                            &
     &                 (1.0_r8-0.6_r8*cloud(i,j)**3)/                   &
     &                 ((zenith+2.7_r8)*vap_p*1.0E-3_r8+                &
     &                  1.085_r8*zenith+0.1_r8)
          END IF
# elif defined DIURNAL_SRFLUX
!
!  SRFLX is reset on each time step in subroutine SET_DATA which
!  interpolates values in the forcing file to the current date.
!  This DIURNAL_SRFLUX option is provided so that SRFLX values
!  corresponding to a greater or equal daily average can be modulated
!  by the local length of day to produce a diurnal cycle with the
!  same daily average as the original data.  This approach assumes
!  the net effect of clouds is incorporated into the SRFLX data.
!
!  Normalization = (1/2*pi)*INTEGRAL{ABS(a+b*COS(t)) dt}  from 0 to 2*pi
!                = (a*ARCCOS(-a/b)+SQRT(b**2-a**2))/pi    for |a| < |b|
!
          srflx(i,j)=MAX(0.0_r8, srflx(i,j))
          IF (ABS(cff1).gt.ABS(cff2)) THEN
            IF (cff1*cff2.gt.0.0_r8) THEN
              cff=cff1                                 ! All day case
              srflx(i,j)=MAX(0.0_r8,                                    &
     &                       srflx(i,j)/cff*                            &
     &                       (cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad)))
            ELSE
              srflx(i,j)=0.0_r8                        ! All night case
            END IF
          ELSE
            cff=(cff1*ACOS(-cff1/cff2)+SQRT((cff2+cff1)*(cff2-cff1)))/pi
	    IF (cff .lt. 10.e-10) THEN
              srflx(i,j)=0.0_r8
            ELSE
              srflx(i,j)=MAX(0.0_r8,                                      &
     &                     srflx(i,j)/cff*                              &
     &                     (cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad)))
            END IF
          END IF
# endif
        END DO
      END DO
#else
!
!-----------------------------------------------------------------------
!  Set incoming solar shortwave radiation (degC m/s).  Usually, the
!  shortwave radiation from input files is Watts/m2 and then converted
!  to degC m/s by multiplying by conversion factor 1/(rho0*Cp) during
!  reading (Fscale). However, we are already inside ROMS kernel here
!  and all the fluxes are kinematic so shortwave radiation units need
!  to be degC m/s.
!-----------------------------------------------------------------------
!
      cff=1.0_r8/(rho0*cp)
# if defined UPWELLING
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          srflx(i,j)=cff*150.0_r8
        END DO
      END DO
# else
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          srflx(i,j)=0.0_r8
        END DO
      END DO
# endif
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx)
#endif

      RETURN
      END SUBROUTINE ana_srflux_tile
