#include "cppdefs.h"
      MODULE obs_k2z_mod
#ifdef SOLVE3D
!
!svn $Id: obs_k2z.F 798 2016-05-21 17:04:58Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine converts observations vertical fractional coordinate   !
!  (Zobs => obs_grid) to depth in meters (obs_depths). The depths are  !
!  negative downwards. This needed by the Ensemble Kalman (EnKF) for   !
!  localization.                                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Imin       Global I-coordinate lower bound of RHO-points.        !
!     Imax       Global I-coordinate upper bound of RHO-points.        !
!     Jmin       Global J-coordinate lower bound of RHO-points.        !
!     Jmax       Global J-coordinate upper bound of RHO-points.        !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     Xmin       Global minimum fractional I-coordinate to consider.   !
!     Xmax       Global maximum fractional I-coordinate to consider.   !
!     Ymin       Global minimum fractional J-coordinate to consider.   !
!     Ymax       Global maximum fractional J-coordinate to consider.   !
!     Mobs       Number of observations.                               !
!     Xobs       Observations X-locations (fractional coordinates).    !
!     Yobs       Observations Y-locations (fractional coordinates).    !
!     Zobs       Observations Z-locations (fractional coordinates or   !
!                  or actual meters).                                  !
!     obs_scale  Observation screenning flag.                          !
!     z          Model grid depths of W-points (meters, 3D array).     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     obs_depths Observations depth (meters, negative downwards.       !
!                                                                      !
!  The interpolation weights matrix, Hmat(1:8), is as follows:         !
!                                                                      !
!                               8____________7                         !
!                               /.          /| (i2,j2,k2)              !
!                              / .         / |                         !
!                            5/___________/6 |                         !
!                             |  .        |  |                         !
!                             |  .        |  |         Grid Cell       !
!                             | 4.........|..|3                        !
!                             | .         |  /                         !
!                             |.          | /                          !
!                  (i1,j1,k1) |___________|/                           !
!                             1           2                            !
!                                                                      !
!  Notice that the indices i2 and j2 are reset when observations are   !
!  located exactly at the eastern and/or northern boundaries. This is  !
!  needed to avoid out-of-range array computations.                    !
!                                                                      !
!  All the observations are assumed to in fractional coordinates with  !
!  respect to RHO-points:                                              !
!                                                                      !
!                                                                      !
!  M      r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r      !
!         :                                                     :      !
!  Mm+.5  v  p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  Mm     r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  Mm-.5  v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!         r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!         v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!         r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  2.5    v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  2.0    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  1.5    v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  1.0    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r      !
!         :  +     |     |     |     |     |     |     |     +  :      !
!  0.5    v  p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p  v      !
!         :                                                     :      !
!  0.0    r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r      !
!                                                                      !
!           0.5   1.5   2.5                          Lm-.5 Lm+.5       !
!                                                                      !
!        0.0   1.0   2.0                                  Lm    L      !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE obs_k2z (ng, Imin, Imax, Jmin, Jmax,                   &
     &                    LBi, UBi, LBj, UBj, LBk, UBk,                 &
     &                    Xmin, Xmax, Ymin, Ymax,                       &
     &                    Mobs, Xobs, Yobs, Zobs, obs_scale,            &
     &                    z, obs_depths)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Mobs
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: obs_scale(:)
      real(r8), intent(in) :: Xobs(:)
      real(r8), intent(in) :: Yobs(:)
      real(r8), intent(in) :: Zobs(:)
      real(r8), intent(in) :: z(LBi:,LBj:,LBk:)

      real(r8), intent(out) :: obs_depths(:)
#  else
      real(r8), intent(in) :: obs_scale(Mobs)
      real(r8), intent(in) :: Xobs(Mobs)
      real(r8), intent(in) :: Yobs(Mobs)
      real(r8), intent(in) :: Zobs(Mobs)
      real(r8), intent(in) :: z(LBi:UBi,LBj:UBj,LBk:UBk)

      real(r8), intent(out) :: obs_depths(Mobs)
#  endif
!
!  Local variable declarations.
!
      logical :: Linterpolate

      integer :: i, ic, iobs, i1, i2, j1, j2, k, k1, k2

      real(r8) :: Zbot, Ztop, dz, p1, p2, q1, q2, r1, r2
      real(r8) :: w11, w12, w21, w22

      real(r8), dimension(8) :: Hmat
!
!-----------------------------------------------------------------------
!  Interpolate vertical fractional coordinate to depths.
!-----------------------------------------------------------------------
!
      DO iobs=1,Mobs
        IF (((Xmin.le.Xobs(iobs)).and.(Xobs(iobs).lt.Xmax)).and.        &
     &      ((Ymin.le.Yobs(iobs)).and.(Yobs(iobs).lt.Ymax))) THEN
          i1=INT(Xobs(iobs))
          j1=INT(Yobs(iobs))
          i2=i1+1
          j2=j1+1
          IF (i2.gt.Imax) THEN
            i2=i1                 ! Observation at the eastern boundary
          END IF
          IF (j2.gt.Jmax) THEN
            j2=j1                 ! Observation at the northern boundary
          END IF
          p2=REAL(i2-i1,r8)*(Xobs(iobs)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(Yobs(iobs)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          w11=p1*q1
          w21=p2*q1
          w22=p2*q2
          w12=p1*q2
          IF (Zobs(iobs).gt.0.0_r8) THEN
            IF (ABS(REAL(UBk,r8)-Zobs(iobs)).lt.1.0E-8_r8) THEN
              Linterpolate=.FALSE.            ! surface observation
              obs_depths(iobs)=0.0_r8
            ELSE
              Linterpolate=.TRUE.             ! fractional level
              k1=MAX(LBk,INT(Zobs(iobs)-0.5_r8))               ! W-point
              k2=MIN(k1+1,UBk)
              r2=REAL(k2-k1,r8)*(Zobs(iobs)-REAL(k1,r8))
              r1=1.0_r8-r2
            END IF
          ELSE
            Linterpolate=.FALSE.              ! already depths in meters
            obs_depths(iobs)=Zobs(iobs)
          END IF
          IF (Linterpolate) THEN
            IF ((r1+r2).gt.0.0_r8) THEN
              Hmat(1)=w11*r1
              Hmat(2)=w21*r1
              Hmat(3)=w22*r1
              Hmat(4)=w12*r1
              Hmat(5)=w11*r2
              Hmat(6)=w21*r2
              Hmat(7)=w22*r2
              Hmat(8)=w12*r2
              obs_depths(iobs)=Hmat(1)*z(i1,j1,k1)+                     &
     &                         Hmat(2)*z(i2,j1,k1)+                     &
     &                         Hmat(3)*z(i2,j2,k1)+                     &
     &                         Hmat(4)*z(i1,j2,k1)+                     &
     &                         Hmat(5)*z(i1,j1,k2)+                     &
     &                         Hmat(6)*z(i2,j1,k2)+                     &
     &                         Hmat(7)*z(i2,j2,k2)+                     &
     &                         Hmat(8)*z(i1,j2,k2)
            END IF
          END IF
        END IF
      END DO

      RETURN
      END SUBROUTINE obs_k2z
#endif
      END MODULE obs_k2z_mod
