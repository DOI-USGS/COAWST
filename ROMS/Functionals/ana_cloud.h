      SUBROUTINE ana_cloud (ng, tile, model)
!
!! svn $Id: ana_cloud.h 1054 2021-03-06 19:47:12Z arango $
!!======================================================================
!! Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets cloud fraction using an analytical expression.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
! Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#include "tile.h"
!
      CALL ana_cloud_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     FORCES(ng) % cloud)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 4)=MyFile
      END IF
!
      RETURN
      END SUBROUTINE ana_cloud
!
!***********************************************************************
      SUBROUTINE ana_cloud_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           cloud)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#ifdef PAPA_CLM
      USE dateclock_mod,   ONLY : caldate
#endif
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
      real(r8), intent(out) :: cloud(LBi:,LBj:)
#else
      real(r8), intent(out) :: cloud(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer  :: i, j
!
      real(r8) :: Cval
      real(dp) :: yday

#ifdef PAPA_CLM
!
      real(dp), dimension(14) :: Coktas =                               &
     &         (/ 6.29_r8, 6.26_r8, 6.31_r8, 6.31_r8, 6.32_r8,          &
     &            6.70_r8, 7.12_r8, 7.26_r8, 6.93_r8, 6.25_r8,          &
     &            6.19_r8, 6.23_r8, 6.31_r8, 6.29_r8          /)

      real(dp), dimension(14) :: Cyday =                                &
     &          (/  0.0_dp,  16.0_dp,  46.0_dp,  75.0_dp, 105.0_dp,     &
     &            136.0_dp, 166.0_dp, 197.0_dp, 228.0_dp, 258.0_dp,     &
     &            289.0_dp, 319.0_dp, 350.0_dp, 365.0_dp           /)
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical cloud fraction (%/100): 0=clear sky, 1:overcast sky.
!-----------------------------------------------------------------------
!
#if defined PAPA_CLM

!  OWS Papa cloud climatology.
!
      CALL caldate (tdays(ng), yd_dp=yday)
      DO i=1,13
        IF ((yday.ge.Cyday(i)).and.(yday.le.Cyday(i+1))) THEN
          Cval=0.125_r8*(Coktas(i  )*(Cyday(i+1)-yday)+                 &
     &                   Coktas(i+1)*(yday-Cyday(i)))/                  &
     &                  (Cyday(i+1)-Cyday(i))
        END IF
      END DO
#elif defined BENCHMARK
      Cval=0.6_r8
#elif defined NJ_BIGHT
      Cval=0.3_r8
#else
      Cval=0.0_r8
#endif

      DO j=JstrT,JendT
        DO i=IstrT,IendT
          cloud(i,j)=Cval
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          cloud)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    cloud)
#endif
!
      RETURN
      END SUBROUTINE ana_cloud_tile
