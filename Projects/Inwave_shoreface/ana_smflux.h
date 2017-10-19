      SUBROUTINE ana_smflux (ng, tile, model)
!
!! svn $Id: ana_smflux.h 429 2009-12-20 17:30:26Z arango $
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface momentum flux (wind stress)     !
!  "sustr" and "svstr" (m2/s2) using an analytical expression.         !
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
      CALL ana_smflux_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % angler,                          &
#ifdef SPHERICAL
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
#else
     &                      GRID(ng) % xr,                              &
     &                      GRID(ng) % yr,                              &
#endif
#ifdef TL_IOMS
     &                      FORCES(ng) % tl_sustr,                      &
     &                      FORCES(ng) % tl_svstr,                      &
#endif
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(24)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_smflux
!
!***********************************************************************
      SUBROUTINE ana_smflux_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            angler,                               &
#ifdef SPHERICAL
     &                            lonr, latr,                           &
#else
     &                            xr, yr,                               &
#endif
#ifdef TL_IOMS
     &                            tl_sustr, tl_svstr,                   &
#endif
     &                            sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
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
      real(r8), intent(in) :: angler(LBi:,LBj:)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# else
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# endif
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
# ifdef TL_IOMS
      real(r8), intent(out) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(out) :: tl_svstr(LBi:,LBj:)
# endif
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
# ifdef SPHERICAL
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
# ifdef TL_IOMS
      real(r8), intent(out) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_svstr(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: Ewind, Nwind, cff, val1, val2, windamp, winddir
#if defined LAKE_SIGNELL
      real(r8) :: cff1, mxst, ramp_u, ramp_time, ramp_d
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  XI-direction (m2/s2) at horizontal U-points.
!-----------------------------------------------------------------------
!
#ifdef BASIN
      val1=5.0E-05_r8*(1.0_r8+TANH((time(ng)-6.0_r8*86400.0_r8)/        &
     &                 (3.0_r8*86400.0_r8)))
      val2=2.0_r8*pi/el(ng)
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          sustr(i,j)=-val1*COS(val2*yr(i,j))
# ifdef TL_IOMS
          tl_sustr(i,j)=-val1*COS(val2*yr(i,j))
# endif
        END DO
      END DO
#elif defined BL_TEST
      Ewind=0.0_r8/rho0
      Nwind=0.3_r8/rho0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          sustr(i,j)=Ewind
# ifdef TL_IOMS
          tl_sustr(i,j)=Ewind
# endif
        END DO
      END DO
#elif defined CANYON
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          sustr(i,j)=5.0E-05_r8*SIN(2.0_r8*pi*tdays(ng)/10.0_r8)*       &
     &               (1.0_r8-TANH((yr(i,j)-0.5_r8*el(ng))/10000.0_r8))
# ifdef TL_IOMS
          tl_sustr(i,j)=5.0E-05_r8*SIN(2.0_r8*pi*tdays(ng)/10.0_r8)*    &
     &               (1.0_r8-TANH((yr(i,j)-0.5_r8*el(ng))/10000.0_r8))
# endif
        END DO
      END DO
#elif defined CHANNEL_NECK
!!    IF ((tdays(ng)-dstart).le.4.0_r8) THEN
!!      windamp=-0.01_r8*SIN(pi*(tdays(ng)-dstart)/8.0_r8)/rho0
!!    ELSE
        windamp=-0.01_r8/rho0
!!    END IF
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          sustr(i,j)=windamp
# ifdef TL_IOMS
          tl_sustr(i,j)=windamp
# endif
        END DO
      END DO
#elif defined MIXED_LAYER
      DO j=JstrT,JendT
         DO i=IstrP,IendT
           sustr(i,j)=0.0001_r8        ! m2/s2
# ifdef TL_IOMS
           tl_sustr(i,j)=0.0001_r8     ! m2/s2
# endif
         END DO
      END DO
#elif defined DOUBLE_GYRE
!!    windamp=user(1)/rho0
      windamp=-0.05_r8/rho0
      val1=2.0_r8*pi/el(ng)
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          sustr(i,j)=windamp*COS(val1*yr(i,j))
# ifdef TL_IOMS
          tl_sustr(i,j)=windamp*COS(val1*yr(i,j))
# endif
        END DO
      END DO
#elif defined FLT_TEST
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          sustr(i,j)=1.0E-03_r8
# ifdef TL_IOMS
          tl_sustr(i,j)=1.0E-03_r8
# endif
        END DO
      END DO
#elif defined LAKE_SIGNELL
      mxst=0.2500_r8          ! N/m2
      ramp_u=15.0_r8           ! start ramp UP at RAMP_UP hours
      ramp_time=10.0_r8       ! ramp from 0 to 1 over RAMP_TIME hours
      ramp_d=50.0_r8          ! start ramp DOWN at RAMP_DOWN hours
      DO j=JstrT,JendT
         DO i=IstrP,IendT
           cff1=MIN((0.5_r8*(TANH((time(ng)/3600.0_r8-ramp_u)/          &
     &                            (ramp_time/5.0_r8))+1.0_r8)),         &
     &              (1.0_r8-(0.5_r8*(TANH((time(ng)/3600.0_r8-ramp_d)/  &
     &                                    (ramp_time/5.0_r8))+1.0_r8))))
           sustr(i,j)=mxst/rho0*cff1
# ifdef TL_IOMS
           tl_sustr(i,j)=mxst/rho0*cff1
# endif
         END DO
      END DO
#elif defined LMD_TEST
      IF (time(ng).le.57600.0_r8) THEN
        windamp=-0.6_r8*SIN(pi*time(ng)/57600.0_r8)*                    &
     &                  SIN(2.0_r8*pi*time(ng)/57600.0_r8)/rho0
      ELSE
        windamp=0.0_r8
      END IF
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          sustr(i,j)=windamp
# ifdef TL_IOMS
          tl_sustr(i,j)=windamp
# endif
        END DO
      END DO
#elif defined NJ_BIGHT
!!    windamp=0.086824313_r8
!!    winddir=0.5714286_r8
!!    if ((tdays(ng)-dstart).le.0.5_r8) then
!!      Ewind=windamp*winddir*SIN(pi*(tdays(ng)-dstart))/rho0
!!      Nwind=windamp*SIN(pi*(tdays(ng)-dstart))/rho0
!!    else
!!      Ewind=windamp*winddir/rho0
!!      Nwind=windamp/rho0
!!    endif
      IF ((tdays(ng)-dstart).le.3.0_r8) THEN
         winddir=60.0_r8
         windamp=0.1_r8
      ELSE IF (((tdays(ng)-dstart).gt.3.0_r8).and.                      &
     &        ((tdays(ng)-dstart).le.4.0_r8)) THEN
         winddir= 60.0_r8*((tdays(ng)-dstart)-2.0_r8)-                  &
     &           120.0_r8*((tdays(ng)-dstart)-2.0_r8)
         windamp=0.0_r8
      ELSE
         winddir=-120.0_r8
         windamp=0.0_r8
      END IF
      Ewind=windamp*COS(pi*winddir/180.0_r8)/rho0
      Nwind=windamp*SIN(pi*winddir/180.0_r8)/rho0
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          val1=0.5_r8*(angler(i-1,j)+angler(i,j))
          sustr(i,j)=Ewind*COS(val1)+Nwind*SIN(val1)
# ifdef TL_IOMS
          tl_sustr(i,j)=Ewind*COS(val1)+Nwind*SIN(val1)
# endif
        END DO
      END DO
#elif defined SED_TOY
      DO j=JstrT,JendT
         DO i=IstrP,IendT
           cff=0.0001_r8
           IF (time(ng).gt.3000.0_r8) THEN
             cff=0.0_r8
           END IF
           sustr(i,j)=cff
# ifdef TL_IOMS
           tl_sustr(i,j)=cff
# endif
         END DO
      END DO
#elif defined INWAVE_SHOREFACE
      DO j=JstrT,JendT
         DO i=IstrP,IendT
          sustr(i,j)=0.0_r8
# ifdef TL_IOMS
          tl_sustr(i,j)=0.0_r8
# endif
         END DO
      END DO
#elif defined UPWELLING
      IF (NSperiodic(ng)) THEN
        DO j=JstrT,JendT
           DO i=IstrP,IendT
            sustr(i,j)=0.0_r8
# ifdef TL_IOMS
            tl_sustr(i,j)=0.0_r8
# endif
          END DO
        END DO
      ELSE IF (EWperiodic(ng)) THEN
        IF ((tdays(ng)-dstart).le.2.0_r8) THEN
          windamp=-0.1_r8*SIN(pi*(tdays(ng)-dstart)/4.0_r8)/rho0
        ELSE
          windamp=-0.1_r8/rho0
        END IF
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            sustr(i,j)=windamp
# ifdef TL_IOMS
            tl_sustr(i,j)=windamp
# endif
          END DO
        END DO
      END IF
#elif defined WINDBASIN
      IF ((tdays(ng)-dstart).le.2.0_r8) THEN
        windamp=-0.1_r8*SIN(pi*(tdays(ng)-dstart)/4.0_r8)/rho0
      ELSE
        windamp=-0.1_r8/rho0
      END IF
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          sustr(i,j)=windamp
# ifdef TL_IOMS
          tl_sustr(i,j)=windamp
# endif
        END DO
      END DO
#else
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          sustr(i,j)=0.0_r8
# ifdef TL_IOMS
          tl_sustr(i,j)=0.0_r8
# endif
        END DO
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  ETA-direction (m2/s2) at horizontal V-points.
!-----------------------------------------------------------------------
!
#if defined BL_TEST
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          svstr(i,j)=Nwind
# ifdef TL_IOMS
          tl_svstr(i,j)=Nwind
# endif
        END DO
      END DO
#elif defined LMD_TEST
      IF (time(ng).le.57600.0_r8) THEN
        windamp=-0.6_r8*SIN(pi*time(ng)/57600.0_r8)*                    &
     &                  COS(2.0_r8*pi*time(ng)/57600.0_r8)/rho0
      ELSE
        windamp=0.0_r8
      END IF
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          svstr(i,j)=windamp
# ifdef TL_IOMS
          tl_svstr(i,j)=windamp
# endif
        END DO
      END DO
#elif defined NJ_BIGHT
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          val1=0.5_r8*(angler(i,j)+angler(i,j-1))
          svstr(i,j)=-Ewind*SIN(val1)+Nwind*COS(val1)
# ifdef TL_IOMS
          tl_svstr(i,j)=-Ewind*SIN(val1)+Nwind*COS(val1)
# endif
        END DO
      END DO
#elif defined SED_TOY
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          svstr(i,j)=0.0_r8
# ifdef TL_IOMS
          tl_svstr(i,j)=0.0_r8
# endif
        END DO
      END DO
#elif defined INWAVE_SHOREFACE
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          svstr(i,j)=0.0_r8
# ifdef TL_IOMS
          tl_svstr(i,j)=0.0_r8
# endif
        END DO
      END DO
#elif defined UPWELLING
      IF (NSperiodic(ng)) THEN
        IF ((tdays(ng)-dstart).le.2.0_r8) THEN
          windamp=-0.1_r8*SIN(pi*(tdays(ng)-dstart)/4.0_r8)/rho0
        ELSE
          windamp=-0.1_r8/rho0
        END IF
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            svstr(i,j)=windamp
# ifdef TL_IOMS
            tl_svstr(i,j)=windamp
# endif
          END DO
        END DO
      ELSE IF (EWperiodic(ng)) THEN
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            svstr(i,j)=0.0_r8
# ifdef TL_IOMS
            tl_svstr(i,j)=0.0_r8
# endif
          END DO
        END DO
      END IF
#else
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          svstr(i,j)=0.0_r8
# ifdef TL_IOMS
          tl_svstr(i,j)=0.0_r8
# endif
        END DO
      END DO
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
#ifdef TL_IOMS
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_sustr)
#endif
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
#ifdef TL_IOMS
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_svstr)
#endif
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
# ifdef TL_IOMS
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_sustr, tl_svstr)
# endif
#endif

      RETURN
      END SUBROUTINE ana_smflux_tile
