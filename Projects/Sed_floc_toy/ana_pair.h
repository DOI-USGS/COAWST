     SUBROUTINE ana_pair (ng, tile, model)
!
!! svn $Id: ana_pair.h 658 2013-04-18 22:17:05Z arango $
!!======================================================================
!! Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets surface air pressure (mb) using an analytical     !
!  expression.                                                         !
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

#include "tile.h"
!
     CALL ana_pair_tile (ng, tile, model,                              &
    &                    LBi, UBi, LBj, UBj,                           &
    &                    IminS, ImaxS, JminS, JmaxS,                   &
    &                    FORCES(ng) % Pair)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
     IF (Lanafile) THEN
#else
     IF (Lanafile.and.(tile.eq.0)) THEN
#endif
       ANANAME(17)=__FILE__
     END IF

     RETURN
     END SUBROUTINE ana_pair
!
!***********************************************************************
     SUBROUTINE ana_pair_tile (ng, tile, model,                        &
    &                          LBi, UBi, LBj, UBj,                     &
    &                          IminS, ImaxS, JminS, JmaxS,             &
    &                          Pair)
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
     real(r8), intent(out) :: Pair(LBi:,LBj:)
#else
     real(r8), intent(out) :: Pair(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
     integer :: i, j
       real(r8) :: dpdx, dx, omega, tstart

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical surface air pressure (mb).
!  (1 mb = 100 Pa = 1 hPa,  1 bar = 1.0e+5 N/m2 = 1.0e+5 dynes/cm2).
!-----------------------------------------------------------------------
!
#if defined BENCHMARK
     DO j=JstrT,JendT
       DO i=IstrT,IendT
         Pair(i,j)=1025.0_r8
       END DO
     END DO
#elif defined BL_TEST
     DO j=JstrT,JendT
       DO i=IstrT,IendT
         Pair(i,j)=1013.48_r8
       END DO
     END DO
#elif defined SED_FLOC_TOY
! 0.02 Pa/m corresponds to a slope of 1 m / 500 km
! dx and dy are 100 m, h = 10 m
!
     dx = 100.0_r8
! CASE 51q, 52q 53q
     dpdx = 0.04_r8

!     tstart = 3600.0_r8*12.0_r8
!     omega = 2.0_r8*3.14159_r8/(12.0_r8*3600.0_r8)
!     dpdx = 0.02_r8 !Pa (bottom stress will be h*dpdx)
     DO j=JstrR-2,JendR+2
       DO i=IstrR-2,IendR+2
          ! convert Pa to millibars with 0.01
          Pair(i,j)=1000.0_r8+0.01_r8*dpdx*FLOAT(i-1)*dx
!          if(j.eq.3)then
!            write(*,*) "Pair:", Pair(i,j)
!	  endif
       END DO
     END DO
#elif defined ONED_TOY
! 0.02 Pa/m corresponds to a slope of 1 m / 500 km
! dx and dy are 10 m, h = 10 m
!
     DO j=JstrT-2,JendT+2
       DO i=IstrT-2,IendT+2
Pair(i,j)=1000.0_r8 + 0.01_r8*user(6)*FLOAT(i)/Lm(ng)
       END DO
     END DO
#else
     ana_pair.h: no values provided for Pair.
#endif

#if !defined ONED_TOY && !defined SED_FLOC_TOY
!  Exchange boundary data. !! Remove this per JWilkin !! CRS
!
     IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
       CALL exchange_r2d_tile (ng, tile,                               &
    &                          LBi, UBi, LBj, UBj,                     &
    &                          Pair)
     END IF
#endif

#ifdef DISTRIBUTE
     CALL mp_exchange2d (ng, tile, model, 1,                           &
    &                    LBi, UBi, LBj, UBj,                           &
    &                    NghostPoints,                                 &
    &                    EWperiodic(ng), NSperiodic(ng),               &
    &                    Pair)
#endif

     RETURN
     END SUBROUTINE ana_pair_tile
