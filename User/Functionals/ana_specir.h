      SUBROUTINE ana_specir (ng, tile, model)
!
!! svn $Id: ana_specir.h 795 2016-05-11 01:42:43Z arango $
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets surface solar downwelling spectral irradiance  !
!  at just beneath the sea surface, Ed(lambda,0-), in micromol quanta  !
!  per meter squared per second.                                       !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Gregg, W.W. and K.L. Carder, 1990:  A simple spectral solar       !
!           irradiance model for cloudless maritime atmospheres,       !
!           Limmol. Oceanogr., 35(8), 1657-1675.                       !
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
      CALL ana_specir_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
     &                      FORCES(ng) % cloud,                         &
     &                      FORCES(ng) % Hair,                          &
     &                      FORCES(ng) % Tair,                          &
     &                      FORCES(ng) % Pair,                          &
     &                      FORCES(ng) % Uwind,                         &
     &                      FORCES(ng) % Vwind,                         &
     &                      FORCES(ng) % SpecIr,                        &
     &                      FORCES(ng) % avcos)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(25)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_specir
!
!***********************************************************************
      SUBROUTINE ana_specir_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            lonr, latr,                           &
     &                            cloud, Hair, Tair, Pair,              &
     &                            Uwind, Vwind,                         &
     &                            SpecIr, avcos)
!***********************************************************************
!
      USE mod_param
      USE mod_eclight
      USE mod_iounits
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
#endif
!
      implicit none
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
      real(r8), intent(in) :: cloud(LBi:,LBj:)
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
      real(r8), intent(out) :: SpecIr(LBi:,LBj:,:)
      real(r8), intent(out) :: avcos(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: UWind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: SpecIr(LBi:UBi,LBj:UBj,NBands)
      real(r8), intent(out) :: avcos(LBi:UBi,LBj:UBj,NBands)
#endif
!
!  Local constant declarations.
!
      real(r8) :: am = 1.0_r8        ! Aerosol type 1-10: ocean to land
      real(r8) :: betalam = 0.55_r8
      real(r8) :: p0 = 29.92_r8      ! Standard pressure (inches of Hg)
      real(r8) :: rex = -1.6364_r8
      real(r8) :: roair = 1200.0_r8  ! Density of air (g/m3)
      real(r8) :: rn = 1.341_r8      ! Index of refraction of pure seawater
      real(r8) :: vis = 15.0_r8      ! Visibility (km)
      real(r8) :: wv = 1.5_r8        ! Precipitable water (cm): water vapor
!
!  Local variable declarations.
!
      integer :: i, iband, ic, j, nc
      integer :: iday, month, year

      real(r8) :: Dangle, Hangle, LatRad, LonRad
      real(r8) :: cff, cff1, cff2, hour, yday
      real(r8) :: alpha, beta, gamma, theta, rtheta, rthetar
      real(r8) :: atra, gtra, otra, rtra, wtra
      real(r8) :: alg, arg, asymp, cosunz, Fa
      real(r8) :: frh, rh, rlam, rlogc
      real(r8) :: rm, rmin, rmo, rmp, rod, rof
      real(r8) :: ros, rospd, rosps, rpls
      real(r8) :: sumx, sumx2, sumxy, sumy
      real(r8) :: taa, tas, to3, wa, wspeed, zenith

      real(r8), dimension(NBands) :: Fo, Edir, Edif, Ed, qlam

      real(r8), dimension(3) :: a_arr, dndr
      real(r8), dimension(3) :: ro     = (/ 0.03_r8, 0.24_r8, 2.00_r8 /)
      real(r8), dimension(3) :: r_arr  = (/ 0.10_r8, 1.00_r8, 10.0_r8 /)
!
#include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  Compute spectral irradiance: Using RADTRAN formulations.
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
!  Conversion constant from E to micromol quanta.
!   1/(Plank*c*(avos#*1.0e6).
!
      cff=1.0E-9_r8/(6.6256E-34_r8*2.998E8_r8*6.023E17_r8)
      DO iband=1,NBands
        qlam(iband)=ec_wave_ab(iband)*cff
      END DO
!
!  Correct solar constant for Earth-Sun distance.
!
      cff=(1.0_r8+0.0167_r8*COS(2.0_r8*pi*(yday-3.0_r8)/365.0_r8))**2
      DO iband=1,NBands
        Fo(iband)=ec_Fobar(iband)*cff
      END DO
!
!  Compute spectral irradiance.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT

          LatRad=latr(i,j)*deg2rad
          LonRad=lonr(i,j)*deg2rad
!
!  Compute Climatological Ozone.
!
          to3=(235.0_r8+(150.0_r8+40.0_r8*                              &
     &         SIN(0.9865_r8*(yday-30.0_r8)*deg2rad)+                   &
     &         20.0_r8*SIN(3.0_r8*LonRad))*                             &
     &         SIN(1.28_r8*LatRad)*SIN(1.28_r8*LatRad))*                &
     &        0.001_r8                                 ! sco3 conversion
!
!  Local daylight is a function of the declination (Dangle) and hour
!  angle adjusted for the local meridian (Hangle-lonr(i,j)*deg2rad).
!
          cosunz=SIN(LatRad)*SIN(Dangle)+                               &
     &           COS(LatRad)*COS(Dangle)*COS(Hangle-lonr(i,j)*deg2rad)
          zenith=ACOS(cosunz)
          theta=zenith*rad2deg
!
!  Model for atmospheric transmittance of solar irradiance through
!  a maritime atmosphere.  Computes direct and diffuse separately.
!  Includes water vapor and oxygen absorption.
!
!  Compute atmospheric path lengths (air mass); pressure-corrected
!
          IF ((theta.ge.0.0_r8).and.(theta.le.90.0_r8)) THEN
!
!  Modified March, 1994 according to Kasten and Young 1989.
!
            rm=1.0_r8/(cosunz+0.50572_r8*(96.07995_r8-theta)**rex)
            rmp=rm*(Pair(i,j)*0.02952756_r8)/p0
            rmo=(1.0_r8+22.0_r8/6370.0_r8)/                             &
     &          SQRT(cosunz*cosunz+44.0_r8/6370.0_r8)
!
!  Computes aerosol parameters according to a simplified version
!  of the Navy marine aerosol model.
!
!  Compute wind speed (24 hour mean is equal to current wind).
!
            wspeed=SQRT(Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j))
!
!  Relative humidity factor, frh.
!
            rh=Hair(i,j)
            IF (rh.ge.100.0_r8) rh=99.9_r8
            frh=((2.0_r8-rh*0.01_r8)/                                   &
                 (6.0_r8*(1.0_r8-rh*0.01_r8)))**0.333_r8
!
!  Size distribution amplitude components.
!
            a_arr(1)=2000.0_r8*am*am
            a_arr(2)=5.866_r8*(wspeed-2.2_r8)
            a_arr(3)=0.01527_r8*(wspeed-2.2_r8)*0.05_r8   !from Hughes 1987
            IF (a_arr(2).lt.0.5_r8) a_arr(2)=0.5_r8
            IF (a_arr(3).lt.0.000014_r8) a_arr(3)=0.000014_r8
!
!  Compute size distribution at three selected radii according to
!  Navy method.
!
            cff=1.0_r8/frh
            DO nc=1,3
              dndr(nc)=0.0_r8
              DO ic=1,3
                arg=LOG(r_arr(nc)/(frh*ro(ic)))
                dndr(nc)=dndr(nc)+a_arr(ic)*EXP(-arg*arg)*cff
              END DO
            END DO
!
!  Least squares approximation
!
            sumx=0.0_r8
            sumy=0.0_r8
            sumxy=0.0_r8
            sumx2=0.0_r8
            DO ic=1,3
              cff1=LOG10(r_arr(ic))
              cff2=LOG10(dndr(ic))
              sumx=sumx+cff1
              sumy=sumy+cff2
              sumxy=sumxy+cff1*cff2
              sumx2=sumx2+cff1*cff1
            END DO
            gamma=sumxy/sumx2
            rlogc=sumy/3.0_r8-gamma*sumx/3.0_r8          ! no used
            alpha=-(gamma+3.0_r8)
!
!  Compute aerosol turbity coefficient, beta.
!
            beta=(3.91_r8/vis)*betalam**alpha
!
!  Compute asymmetry parameter -- a function of alpha.
!
            IF (alpha.gt.1.2_r8) THEN
              asymp=0.65_r8
            ELSE IF (alpha .lt. 0.0_r8) THEN
              asymp=0.82_r8
            ELSE
              asymp=-0.14167_r8*alpha+0.82_r8
            END IF
!
!  Single scattering albedo at 550; function of RH.
!
            wa=(-0.0032_r8*am+0.972_r8)*EXP(0.000306_r8*rh)
!
!  Forward scattering probability.
!
            alg=LOG(1.0_r8-asymp)
            Fa=1.0_r8-0.5_r8*                                           &
     &         EXP((alg*(1.459_r8+alg*(0.1595_r8+alg*0.4129_r8))+       &
     &              alg*(0.0783_r8+alg*(-0.3824_r8-alg*0.5874_r8))*     &
     &             cosunz)*cosunz)
!
!  Compute surface reflectance for direct (rod) and diffuse (ros)
!  components separately, as a function of theta, wind speed or
!  stress.
!
!  Foam and diffuse reflectance
!
            IF (wspeed.gt.4.0_r8) THEN
              IF (wspeed.le.7.0_r8) THEN
                rof=roair*(0.00062_r8+0.00156_r8/wspeed)*               &
     &              0.000022_r8*wspeed*wspeed-0.00040_r8
              ELSE
                rof=(roair*(0.00049_r8+0.000065_r8*wspeed)*             &
     &               0.000045_r8-0.000040_r8)*wspeed*wspeed
              END IF
              rosps=0.057_r8
            ELSE
              rof=0.0_r8
              rosps=0.066_r8
            END IF
!
!  Direct Fresnel reflectance for theta < 40, wspeed < 2 m/s.
!
            IF ((theta.lt.40.0_r8).or.(wspeed.lt.2.0_r8)) THEN
              IF (theta.eq.0.0_r8) THEN
                rospd=0.0211_r8
              ELSE
                rtheta=zenith
                rthetar=ASIN(SIN(rtheta)/rn)
                rmin=rtheta-rthetar
                rpls=rtheta+rthetar
                rospd=0.5_r8*((SIN(rmin)*SIN(rmin))/                    &
     &                        (SIN(rpls)*SIN(rpls))+                    &
     &                        (TAN(rmin)*TAN(rmin))/                    &
     &                        (TAN(rpls)*TAN(rpls)))
              END IF
!
!  Empirical fit otherwise.
!
            ELSE
              rospd=0.0253_r8*EXP((-0.000714_r8*wspeed+0.0618_r8)*      &
     &                            (theta-40.0_r8))
            END IF
!
!  Reflectance totals.
!
            rod=rospd+rof
            ros=rosps+rof
!
!  Compute spectral irradiance for selected spectral bands.
!
            DO iband=1,NBands
              rlam=ec_wave_ab(iband)*0.001_r8
!
!  Transmittance, Rayleigh, by the method of Bird.
!
              rtra=EXP(-rmp/(115.6406_r8*rlam**4-1.335_r8*rlam**2))
!
!  Ozone.
!
              otra=EXP(-ec_aoz(iband)*to3*rmo)
!
!  Aerosols.
!
              arg=beta*rm*rlam**(-alpha)
              atra=EXP(-arg)
              taa=EXP(-(1.0_r8-wa)*arg)
              tas=EXP(-wa*arg)
!
!  Oxygen/gases.
!
              gtra=EXP((-1.41_r8*ec_ag(iband)*rmp)/                     &
     &                 ((1.0_r8+118.3_r8*ec_ag(iband)*rmp)**0.45_r8))
!
!  Water Vapor.
!
              wtra=EXP((-0.2385_r8*ec_aw(iband)*wv*rm)/                 &
     &                 ((1.0_r8+20.07_r8*ec_aw(iband)*wv*rm)**0.45_r8))
!
!  Direct irradiance.
!
              Edir(iband)=Fo(iband)*cosunz*rtra*otra*atra*gtra*         &
     &                    wtra*(1.0_r8-rod)
!
!  Total diffuse irradiance.
!
              Edif(iband)=(1.0_r8-ros)*                                 &
     &                    Fo(iband)*cosunz*gtra*wtra*otra*              &
     &                    (taa*0.5_r8*(1.0_r8-rtra**0.95_r8)+           &
     &                     taa*Fa*(1.0_r8-tas)*rtra**1.5_r8)
!
!  Cloud effects approximations, Kasten and Czeplak (1980).
!  (See Hydrolight Technical Notes).
!
              IF (cloud(i,j).gt.0.25_r8) THEN
                Ed(iband)=(Edir(iband)+Edif(iband))*                    &
     &                    (1.0_r8-0.75_r8*cloud(i,j)**3.4_r8)
                Edif(iband)=Ed(iband)*                                  &
     &                      (0.3_r8+0.7_r8*cloud(i,j)**2.0_r8)
              ELSE
                Ed(iband)=Edir(iband)+Edif(iband)
              END IF
!
!  Convert from W/cm/um to micromole quanta/m2/s.
!
              SpecIr(i,j,iband)=Ed(iband)*10.0_r8*qlam(iband)
!
!  Correction of zenith angle after crossing air/sea interface.
!
              cff1=COS(ASIN((SIN(zenith))/rn))
              cff2=Edif(iband)/Ed(iband)
              avcos(i,j,iband)=cff1*(1.0_r8-cff2)+0.86_r8*cff2
            END DO
          ELSE
            DO iband=1,NBands
              SpecIr(i,j,iband)=0.0_r8
              avcos(i,j,iband)=0.66564_r8
            END DO
          END IF
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, NBands,          &
     &                          SpecIr)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, NBands,          &
     &                          avcos)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, NBands,                &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    SpecIr, avcos)
#endif

      RETURN
      END SUBROUTINE ana_specir_tile
