#include "cppdefs.h"
      MODULE ad_extract_obs_mod
#if defined ADJOINT && defined FOUR_DVAR
!
!svn $Id: ad_extract_obs.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine computes the adjoint misfit forcing terms at the    !
!  observations locations (Xobs,Yobs,Zobs).  The extraction is done    !
!  via linear interpolation.  The (Xobs,Yobs) locations  must be in    !
!  grid coordinates. Zobs can be in grid coordinates (Zobs >= 0) or    !
!  actual depths (Zobs < 0), if applicable.                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Imin       Global I-coordinate lower bound threshold for         !
!                  requested state field type.                         !
!     Imax       Global I-coordinate upper bound threshold for         !
!                  requested state field type.                         !
!     Jmin       Global J-coordinate lower bound threshold for         !
!                  requested state field type.                         !
!     Jmax       Global J-coordinate upper bound threshold for         !
!                  requested state field type.                         !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     ifield     State field identification to process.                !
!     Mobs       Observation dimension in the calling program.         !
!     NobsSTR    Starting observation to process.                      !
!     NobsEND    Last observations to process.                         !
!     Xmin       Global minimum fractional I-coordinate to consider.   !
!     Xmax       Global maximum fractional I-coordinate to consider.   !
!     Ymin       Global minimum fractional J-coordinate to consider.   !
!     Ymax       Global maximum fractional J-coordinate to consider.   !
!     time       Current model time (secods).                          !
!     dt         Model baroclinic time-step (seconds).                 !
!     ObsType    Observations type.                                    !
!     Tobs       Observations time (days).                             !
!     Xobs       Observations X-locations (grid coordinates).          !
!     Yobs       Observations Y-locations (grid coordinates).          !
!     Zobs       Observations Z-locations (grid coordinates or meters).!
!     ad_Aobs    Adjoint forcing term at observations locations.       !
!     Adepth     Depths (meter; negative).                             !
!     Amask      Land-sea masking.                                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ObsScale   Observation screenning flag.                          !
!     ad_A       Adjoint equations forcing array (2D or 3D).           !
!     ad_Aobs    Adjoint forcing term at observations locations.       !
!     Zobs       Observations Z-locations (grid coordinates).          !
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

      implicit none

      PUBLIC :: ad_extract_obs2d
# ifdef SOLVE3D
      PUBLIC :: ad_extract_obs3d
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_extract_obs2d (ng, Imin, Imax, Jmin, Jmax,          &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ifield, Mobs, NobsSTR, NobsEND,      &
     &                             Xmin, Xmax, Ymin, Ymax,              &
     &                             time, dt,                            &
     &                             ObsType, ObsScale,                   &
     &                             Tobs, Xobs, Yobs,                    &
     &                             ad_A,                                &
# ifdef MASKING
     &                             Amask,                               &
# endif
     &                             ad_Aobs)
!***********************************************************************

# ifdef BGQC
!
      USE mod_scalars, ONLY: inner, outer
      USE mod_fourdvar, ONLY: wrtNLmod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ifield, Mobs, NobsSTR, NobsEND
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
      real(r8), intent(in) :: time, dt
!
# ifdef ASSUMED_SHAPE
      integer, intent(in) :: ObsType(:)

      real(r8), intent(in) :: Tobs(:)
      real(r8), intent(in) :: Xobs(:)
      real(r8), intent(in) :: Yobs(:)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ObsScale(:)
      real(r8), intent(inout) :: ad_Aobs(:)
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)

# else
      integer, intent(in) :: ObsType(Mobs)

      real(r8), intent(in) :: Tobs(Mobs)
      real(r8), intent(in) :: Xobs(Mobs)
      real(r8), intent(in) :: Yobs(Mobs)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ObsScale(Mobs)
      real(r8), intent(inout) :: ad_Aobs(Mobs)
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: ic, iobs, i1, i2, j1, j2
      real(r8) :: TimeLB, TimeUB, p1, p2, q1, q2, wsum

      real(r8), dimension(8) :: Hmat
!
!-----------------------------------------------------------------------
!  Interpolate from requested 2D state field when appropriate.
!-----------------------------------------------------------------------
!
      TimeLB=(time-0.5_r8*dt)/86400.0_r8
      TimeUB=(time+0.5_r8*dt)/86400.0_r8
!
      DO iobs=NobsSTR,NobsEND
        IF ((ObsType(iobs).eq.ifield).and.                              &
     &      ((TimeLB.le.Tobs(iobs)).and.(Tobs(iobs).lt.TimeUB)).and.    &
     &      ((Xmin.le.Xobs(iobs)).and.(Xobs(iobs).lt.Xmax)).and.        &
     &      ((Ymin.le.Yobs(iobs)).and.(Yobs(iobs).lt.Ymax))) THEN
          IF (ObsType(iobs).eq.2) THEN
            i1=INT(Xobs(iobs)+0.5_r8)        ! 2D U-grid type variable
            j1=INT(Yobs(iobs))
          ELSE IF (ObsType(iobs).eq.3) THEN
            i1=INT(Xobs(iobs))               ! 2D V-grid type variable
            j1=INT(Yobs(iobs)+0.5_r8)
          ELSE
            i1=INT(Xobs(iobs))               ! 2D RHO-grid type variable
            j1=INT(Yobs(iobs))
          END IF
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
          Hmat(1)=p1*q1
          Hmat(2)=p2*q1
          Hmat(3)=p2*q2
          Hmat(4)=p1*q2
# ifdef MASKING
          Hmat(1)=Hmat(1)*Amask(i1,j1)
          Hmat(2)=Hmat(2)*Amask(i2,j1)
          Hmat(3)=Hmat(3)*Amask(i2,j2)
          Hmat(4)=Hmat(4)*Amask(i1,j2)
          wsum=0.0_r8
          DO ic=1,4
            wsum=wsum+Hmat(ic)
          END DO
          IF (wsum.gt.0.0_r8) THEN
            wsum=1.0_r8/wsum
            DO ic=1,4
              Hmat(ic)=Hmat(ic)*wsum
            END DO
          END IF
# endif
!>        Aobs(iobs)=Hmat(1)*A(i1,j1)+                                  &
!>   &               Hmat(2)*A(i2,j1)+                                  &
!>   &               Hmat(3)*A(i2,j2)+                                  &
!>   &               Hmat(4)*A(i1,j2)
!>
          ad_A(i1,j1)=ad_A(i1,j1)+Hmat(1)*ad_Aobs(iobs)
          ad_A(i2,j1)=ad_A(i2,j1)+Hmat(2)*ad_Aobs(iobs)
          ad_A(i2,j2)=ad_A(i2,j2)+Hmat(3)*ad_Aobs(iobs)
          ad_A(i1,j2)=ad_A(i1,j2)+Hmat(4)*ad_Aobs(iobs)
          ad_Aobs(iobs)=0.0_r8
# ifdef BGQC
#  ifdef IS4DVAR
          IF (inner.eq.0.and.wrtNLmod(ng)) THEN
#  else
          IF (inner.eq.0.and.outer.eq.0) THEN
#  endif
# endif
# ifdef MASKING
             IF (wsum.gt.0.0_r8) ObsScale(iobs)=1.0_r8
# else
             ObsScale(iobs)=1.0_r8
# endif
# ifdef BGQC
          END IF
# endif
        END IF
      END DO

      RETURN
      END SUBROUTINE ad_extract_obs2d

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE ad_extract_obs3d (ng, Imin, Imax, Jmin, Jmax,          &
     &                             LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                             ifield, Mobs, NobsSTR, NobsEND,      &
     &                             Xmin, Xmax, Ymin, Ymax,              &
     &                             time, dt,                            &
     &                             ObsType, ObsScale,                   &
     &                             Tobs, Xobs, Yobs, Zobs,              &
     &                             ad_A, Adepth,                        &
#  ifdef MASKING
     &                             Amask,                               &
#  endif
     &                             ad_Aobs)
!***********************************************************************
!
      USE mod_param
#  ifdef BGQC
      USE mod_scalars, ONLY: inner, outer
      USE mod_fourdvar, ONLY: wrtNLmod
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, Mobs, NobsSTR, NobsEND
!
      real(r8), intent(in) :: Xmin, Xmax, Ymin, Ymax
      real(r8), intent(in) :: time, dt
!
#  ifdef ASSUMED_SHAPE
      integer, intent(in) :: ObsType(:)

      real(r8), intent(in) :: Tobs(:)
      real(r8), intent(in) :: Xobs(:)
      real(r8), intent(in) :: Yobs(:)
      real(r8), intent(in) :: Adepth(LBi:,LBj:,LBk:)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: ObsScale(:)
      real(r8), intent(inout) :: ad_Aobs(:)
      real(r8), intent(inout) :: Zobs(:)
      real(r8), intent(inout) :: ad_A(LBi:,LBj:,LBk:)
#  else
      integer, intent(in) :: ObsType(Mobs)

      real(r8), intent(in) :: Tobs(Mobs)
      real(r8), intent(in) :: Xobs(Mobs)
      real(r8), intent(in) :: Yobs(Mobs)
      real(r8), intent(in) :: Adepth(LBi:UBi,LBj:UBj,LBk:UBk)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: ObsScale(Mobs)
      real(r8), intent(inout) :: ad_Aobs(Mobs)
      real(r8), intent(inout) :: Zobs(Mobs)
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj,LBk:UBk)
#  endif
!
!  Local variable declarations.
!
      integer :: i, ic, iobs, i1, i2, j1, j2, k, k1, k2

      real(r8) :: TimeLB, TimeUB, Zbot, Ztop, dz, p1, p2, q1, q2, r1, r2
      real(r8) :: w11, w12, w21, w22, wsum

      real(r8), dimension(8) :: Hmat
!
!-----------------------------------------------------------------------
!  Interpolate from requested 3D state field.
!-----------------------------------------------------------------------
!
      TimeLB=(time-0.5_r8*dt)/86400.0_r8
      TimeUB=(time+0.5_r8*dt)/86400.0_r8
!
      DO iobs=NobsSTR,NobsEND
        IF ((ObsType(iobs).eq.ifield).and.                              &
     &      ((TimeLB.le.Tobs(iobs)).and.(Tobs(iobs).lt.TimeUB)).and.    &
     &      ((Xmin.le.Xobs(iobs)).and.(Xobs(iobs).lt.Xmax)).and.        &
     &      ((Ymin.le.Yobs(iobs)).and.(Yobs(iobs).lt.Ymax))) THEN
          IF (ObsType(iobs).eq.4) THEN
            i1=INT(Xobs(iobs)+0.5_r8)        ! 3D U-grid type variable
            j1=INT(Yobs(iobs))
          ELSE IF (ObsType(iobs).eq.5) THEN
            i1=INT(Xobs(iobs))               ! 3D V-grid type variable
            j1=INT(Yobs(iobs)+0.5_r8)
          ELSE
            i1=INT(Xobs(iobs))               ! 3D RHO-grid type variable
            j1=INT(Yobs(iobs))
          END IF
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
            k1=MAX(1,INT(Zobs(iobs)))          ! Positions in fractional
            k2=MIN(k1+1,N(ng))                 ! levels
            r2=REAL(k2-k1,r8)*(Zobs(iobs)-REAL(k1,r8))
            r1=1.0_r8-r2
          ELSE
            Ztop=Adepth(i1,j1,N(ng))
            Zbot=Adepth(i1,j1,1    )
            IF (Zobs(iobs).ge.Ztop) THEN
              r1=0.0_r8                        ! If shallower, ignore.
              r2=0.0_r8
# ifdef BGQC
#  ifdef IS4DVAR
              IF (inner.eq.0.and.wrtNLmod(ng)) ObsScale(iobs)=0.0_r8
#  else
              IF (inner.eq.0.and.outer.eq.0) ObsScale(iobs)=0.0_r8
#  endif
# else
              ObsScale(iobs)=0.0_r8
# endif
            ELSE IF (Zbot.ge.Zobs(iobs)) THEN
              r1=0.0_r8                        ! If deeper, ignore.
              r2=0.0_r8
# ifdef BGQC
#  ifdef IS4DVAR
              IF (inner.eq.0.and.wrtNLmod(ng)) ObsScale(iobs)=0.0_r8
#  else
              IF (inner.eq.0.and.outer.eq.0) ObsScale(iobs)=0.0_r8
#  endif
# else
              ObsScale(iobs)=0.0_r8
# endif
            ELSE
              DO k=N(ng),2,-1                  ! Otherwise, interpolate
                Ztop=Adepth(i1,j1,k  )         ! to fractional level
                Zbot=Adepth(i1,j1,k-1)
                IF ((Ztop.gt.Zobs(iobs)).and.(Zobs(iobs).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=Adepth(i1,j1,k2)-Adepth(i1,j1,k1)
              r2=(Zobs(iobs)-Adepth(i1,j1,k1))/dz
              r1=1.0_r8-r2
              Zobs(iobs)=REAL(k1,r8)+r2        ! overwrite
            END IF
          END IF
          IF ((r1+r2).gt.0.0_r8) THEN
            Hmat(1)=w11*r1
            Hmat(2)=w21*r1
            Hmat(3)=w22*r1
            Hmat(4)=w12*r1
            Hmat(5)=w11*r2
            Hmat(6)=w21*r2
            Hmat(7)=w22*r2
            Hmat(8)=w12*r2
#  ifdef MASKING
            Hmat(1)=Hmat(1)*Amask(i1,j1)
            Hmat(2)=Hmat(2)*Amask(i2,j1)
            Hmat(3)=Hmat(3)*Amask(i2,j2)
            Hmat(4)=Hmat(4)*Amask(i1,j2)
            Hmat(5)=Hmat(5)*Amask(i1,j1)
            Hmat(6)=Hmat(6)*Amask(i2,j1)
            Hmat(7)=Hmat(7)*Amask(i2,j2)
            Hmat(8)=Hmat(8)*Amask(i1,j2)
            wsum=0.0_r8
            DO ic=1,8
              wsum=wsum+Hmat(ic)
            END DO
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              DO ic=1,8
                Hmat(ic)=Hmat(ic)*wsum
              END DO
            END IF
#  endif
!>          Aobs(iobs)=Hmat(1)*A(i1,j1,k1)+                             &
!>   &                 Hmat(2)*A(i2,j1,k1)+                             &
!>   &                 Hmat(3)*A(i2,j2,k1)+                             &
!>   &                 Hmat(4)*A(i1,j2,k1)+                             &
!>   &                 Hmat(5)*A(i1,j1,k2)+                             &
!>   &                 Hmat(6)*A(i2,j1,k2)+                             &
!>   &                 Hmat(7)*A(i2,j2,k2)+                             &
!>   &                 Hmat(8)*A(i1,j2,k2)
!>
            ad_A(i1,j1,k1)=ad_A(i1,j1,k1)+Hmat(1)*ad_Aobs(iobs)
            ad_A(i2,j1,k1)=ad_A(i2,j1,k1)+Hmat(2)*ad_Aobs(iobs)
            ad_A(i2,j2,k1)=ad_A(i2,j2,k1)+Hmat(3)*ad_Aobs(iobs)
            ad_A(i1,j2,k1)=ad_A(i1,j2,k1)+Hmat(4)*ad_Aobs(iobs)
            ad_A(i1,j1,k2)=ad_A(i1,j1,k2)+Hmat(5)*ad_Aobs(iobs)
            ad_A(i2,j1,k2)=ad_A(i2,j1,k2)+Hmat(6)*ad_Aobs(iobs)
            ad_A(i2,j2,k2)=ad_A(i2,j2,k2)+Hmat(7)*ad_Aobs(iobs)
            ad_A(i1,j2,k2)=ad_A(i1,j2,k2)+Hmat(8)*ad_Aobs(iobs)
            ad_Aobs(iobs)=0.0_r8
# ifdef BGQC
#  ifdef IS4DVAR
            IF (inner.eq.0.and.wrtNLmod(ng)) THEN
#  else
            IF (inner.eq.0.and.outer.eq.0) THEN
#  endif
# endif
#  ifdef MASKING
               IF (wsum.gt.0.0_r8) ObsScale(iobs)=1.0_r8
#  else
               ObsScale(iobs)=1.0_r8
#  endif
# ifdef BGQC
            END IF
# endif
          END IF
        END IF
      END DO

      RETURN
      END SUBROUTINE ad_extract_obs3d
# endif
#endif
      END MODULE ad_extract_obs_mod
