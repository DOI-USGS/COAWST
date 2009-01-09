      SUBROUTINE ana_srflux (ng, tile, model)
!
!! svn $Id: ana_srflux.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
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
      CALL ana_srflux_tile (ng, model, Istr, Iend, Jstr, Jend,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
#ifdef ALBEDO
     &                      FORCES(ng) % cloud,                         &
     &                      FORCES(ng) % Hair,                          &
     &                      FORCES(ng) % Tair,                          &
     &                      FORCES(ng) % Pair,                          &
#endif
     &                      FORCES(ng) % srflx)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(27),'(a,a)') TRIM(Adir), '/ana_srflux.h'
      END IF

      RETURN
      END SUBROUTINE ana_srflux
!
!***********************************************************************
      SUBROUTINE ana_srflux_tile (ng, model, Istr, Iend, Jstr, Jend,    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            lonr, latr,                           &
#ifdef ALBEDO 
     &                            cloud, Hair, Tair, Pair,              &
#endif
     &                            srflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# ifdef ALBEDO
      real(r8), intent(in) :: cloud(LBi:,LBj:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
# endif
      real(r8), intent(out) :: srflx(LBi:,LBj:)
#else
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# ifdef ALBEDO
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
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j
#if defined ALBEDO || defined DIURNAL_SRFLUX
      integer :: iday, month, year
      real(r8) :: Dangle, Hangle, LatRad
      real(r8) :: cff, cff1, cff2, hour, yday
# ifdef ALBEDO
      real(r8) :: Rsolar, e_sat, vap_p, zenith
# endif
#endif

#include "set_bounds.h"

#if defined ALBEDO || defined DIURNAL_SRFLUX
!
!-----------------------------------------------------------------------
!  Compute shortwave radiation (degC m/s):
!
!  ALBEDO option: Compute shortwave radiation flux using the Laevastu
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
# ifdef ALBEDO
      Rsolar=Csolar/(rho0*Cp)
# endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!
!  Local daylight is a function of the declination (Dangle) and hour 
!  angle adjusted for the local meridian (Hangle-lonr(i,j)/15.0). 
!  The 15.0 factor is because the sun moves 15 degrees every hour.
!
          LatRad=latr(i,j)*deg2rad
          cff1=SIN(LatRad)*SIN(Dangle)
          cff2=COS(LatRad)*COS(Dangle)
# if defined ALBEDO
!
!  Estimate variation in optical thickness of the atmosphere over
!  the course of a day under cloudless skies. To obtain net incoming 
!  shortwave radiation multiply by (1.0-0.6*c**3), where c is the 
!  fractional cloud cover.
!
          srflx(i,j)=0.0_r8
          zenith=cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad/15.0_r8)
          IF (zenith.gt.0.0_r8) THEN
            cff=(0.7859_r8+0.03477_r8*Tair(i,j))/                       &
     &          (1.0_r8+0.00412_r8*Tair(i,j))
            e_sat=10.0_r8**cff
!!
!! If relative humidity in kg/kg.
!!
!!          vap_p=Pair(i,j)*Hair(i,j)/(0.62197_r8+0.378_r8*Hair(i,j))
!!
            vap_p=e_sat*Hair(i,j)
            srflx(i,j)=Rsolar*zenith*zenith*                            &
     &                 (1.0_r8-0.6_r8*cloud(i,j)**3)/                   &
     &                 ((zenith+2.7_r8)*vap_p*1.0E-5_r8+                &
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
!  Normalization factor = INTEGRAL{ABS(a+b*COS(t)) dt} from 0 to 2*pi 
!                       = (a*ARCCOS(-a/b)+SQRT(b**2-a**2))/pi
!  
          IF ((ABS(cff1)+1.E-8_r8).gt.ABS(cff2)) THEN
            IF (cff1*cff2.gt.0.0_r8) THEN
              cff=cff1                                 ! All day case
              srflx(i,j)=MAX(0.0_r8,                                    &
     &                       srflx(i,j)/cff*                            &
     &                       (cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad)))
            ELSE
              srflx(i,j)=0.0_r8                        ! All night case
            END IF
          ELSE
            cff=(cff1*ACOS(-cff1/cff2)+SQRT(cff2*cff2-cff1*cff1))/pi
            srflx(i,j)=MAX(0.0_r8,                                      &
     &                     srflx(i,j)/cff*                              &
     &                     (cff1+cff2*COS(Hangle-lonr(i,j)*deg2rad)))
          END IF
# endif
        END DO
      END DO
#else
!
!-----------------------------------------------------------------------
!  Set incoming solar shortwave radiation (W/m2).
!-----------------------------------------------------------------------
!
# if defined MY_APPLICATION
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          srflx(i,j)=???
        END DO
      END DO
# else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          srflx(i,j)=0.0_r8
        END DO
      END DO
# endif
#endif

#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        srflx)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, model, 1, Istr, Iend, Jstr, Jend,         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    srflx)
#endif

      RETURN
      END SUBROUTINE ana_srflux_tile
