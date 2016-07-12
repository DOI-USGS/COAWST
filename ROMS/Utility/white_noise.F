#include "cppdefs.h"
      MODULE white_noise_mod
!
!svn $Id: white_noise.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines generates "white noise" arrays with random numbers   !
!  over specified range.  These random numbers are scaled to insure,   !
!  approximately, zero mean expectation and unit variance.             !
!                                                                      !
!  Several random number generation schemes are allowed:               !
!                                                                      !
!  Rscheme = 0         Use F90 intrinsic random numbers generator,     !
!                        0 <= R < 1                                    !
!  Rscheme = 1         Random deviates with Gaussian distribution,     !
!                        -1 < R < 1                                    !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!  white_noise1d       Random numbers for 1D arrays                    !
!  white_noise2d       Random numbers for 2D arrays                    !
!  white_noise2d_bry   Random numbers for 2D boundary arrays           !
!  white_noise3d       Random numbers for 3D arrays                    !
!  white_noise3d_bry   Random numbers for 3D boundary arrays           !
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PRIVATE

      PUBLIC :: white_noise1d
      PUBLIC :: white_noise2d
      PUBLIC :: white_noise2d_bry
#ifdef SOLVE3D
      PUBLIC :: white_noise3d
      PUBLIC :: white_noise3d_bry
#endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE white_noise1d (ng, model, Rscheme,                     &
     &                          Imin, Imax, LBi, UBi,                   &
     &                          Rmin, Rmax, R)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
      USE nrutil, ONLY : gasdev
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Rscheme
      integer, intent(in) :: Imin, Imax, LBi, UBi

      real(r8), intent(out) :: Rmin, Rmax

#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: R(LBi:)
#else
      real(r8), intent(out) :: R(LBi:UBi)
#endif
!
!  Local variable declarations.
!
      integer :: i, ic

      real(r8), parameter :: fac = 2.0_r8 * 1.73205080756887720_r8

      real(r8), dimension((UBi-LBi+1)) :: random
!
!-----------------------------------------------------------------------
!  Generate an array of random numbers.
!-----------------------------------------------------------------------
!
!  Initialize output array.
!
      DO i=LBi,UBi
        R(i)=0.0_r8
      END DO
!
!  F90 intrinsic pseudorandom numbers with an uniform distribution
!  over the range 0 <= R < 1.
!
      IF (Rscheme.eq.0) THEN
        CALL random_number (R(LBi:))
!
!  Scale (fac=2*SQRT(3)) the random number to insure the expectation
!  mean to be approximately zero, E(R)=0, and the expectation variance
!  to be approximately unity, E(R^2)=1 (Bennett, 2002; book page 72).
!
!  Recall that,
!
!      E(R) = sum[R f(R)]
!
!  where F(R) is the random variable probability distribution.
!
        DO i=Imin,Imax
          R(i)=fac*(R(i)-0.5_r8)
        END DO
        Rmin=MINVAL(R)
        Rmax=MAXVAL(R)
!
!  Random deviates with a Gaussian (normal) distribuiton over the
!  range from -1 to 1.
!
      ELSE IF (Rscheme.eq.1) THEN
        IF (Master) THEN
          CALL gasdev (random)
          Rmin=MINVAL(random)
          Rmax=MAXVAL(random)
        END IF
        ic=0
        DO i=Imin,Imax
          ic=ic+1
          R(i)=random(ic)
        END DO
      END IF

      RETURN
      END SUBROUTINE white_noise1d

!
!***********************************************************************
      SUBROUTINE white_noise2d (ng, model, gtype, Rscheme,              &
     &                          Imin, Imax, Jmin, Jmax,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Rmin, Rmax, R)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_scatter2d
#endif
      USE nrutil, ONLY : gasdev
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, gtype, Rscheme
      integer, intent(in) :: Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj

      real(r8), intent(out) :: Rmin, Rmax

      real(r8), intent(out) :: R(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: Npts, i, ic, j

      real(r8), parameter :: fac = 2.0_r8 * 1.73205080756887720_r8

      real(r8), allocatable :: random(:)
!
!-----------------------------------------------------------------------
!  Generate an array of random numbers.
!-----------------------------------------------------------------------
!
!  Initialize output array.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          R(i,j)=0.0_r8
        END DO
      END DO
!
!  F90 intrinsic pseudorandom numbers with an uniform distribution
!  over the range 0 <= R < 1.
!
      IF (Rscheme.eq.0) THEN
        CALL random_number (R(LBi:,LBj:))
!
!  Scale (fac=2*SQRT(3)) the random number to insure the expectation
!  mean to be approximately zero, E(R)=0, and the expectation variance
!  to be approximately unity, E(R^2)=1 (Bennett, 2002; book page 72).
!
!  Recall that,
!
!      E(R) = sum[R f(R)]
!
!  where F(R) is the random variable probability distribution.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            R(i,j)=fac*(R(i,j)-0.5_r8)
          END DO
        END DO
        Rmin=MINVAL(R)
        Rmax=MAXVAL(R)
!
!  Random deviates with a Gaussian (normal) distribuiton with zero
!  mean and unit variance.
!
      ELSE IF (Rscheme.eq.1) THEN
        Npts=(Lm(ng)+2)*(Mm(ng)+2)
        IF (.not.allocated(random)) THEN
          allocate ( random(Npts+2) )
        END IF
        IF (Master) THEN
          CALL gasdev (random)
          Rmin=random(1)
          Rmax=random(1)
          DO ic=1,Npts
            Rmin=MIN(Rmin,random(ic))
            Rmax=MAX(Rmax,random(ic))
          END DO
        END IF
#ifdef DISTRIBUTE
        Npts=SIZE(random)
        CALL mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,               &
     &                     NghostPoints, gtype, Rmin, Rmax,             &
# if defined READ_WATER && defined MASKING
     &                     (Lm(ng)+2)*(Mm(ng)+2),                       &
     &                     SCALARS(ng)%IJwater(:,gtype),                &
# endif
     &                     Npts, random, R)
#else
        ic=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ic=ic+1
            R(i,j)=random(ic)
          END DO
        END DO
#endif
        IF (allocated(random)) THEN
          deallocate (random)
        END IF
      END IF

      RETURN
      END SUBROUTINE white_noise2d

!
!***********************************************************************
      SUBROUTINE white_noise2d_bry (ng, tile, model, boundary, Rscheme, &
     &                              Imin, Imax,                         &
     &                              LBij, UBij,                         &
     &                              Rmin, Rmax, R)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf
      USE mp_exchange_mod, ONLY : mp_exchange2d_bry
#endif
      USE nrutil, ONLY : gasdev
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, boundary, Rscheme
      integer, intent(in) :: Imin, Imax
      integer, intent(in) :: LBij, UBij

      real(r8), intent(out) :: Rmin, Rmax

#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: R(LBij:)
#else
      real(r8), intent(out) :: R(LBij:UBij)
#endif
!
!  Local variable declarations.
!
      integer :: Ioff, i, ic

      real(r8), parameter :: fac = 2.0_r8 * 1.73205080756887720_r8

      real(r8), dimension(UBij-LBij+1) :: random
!
!-----------------------------------------------------------------------
!  Generate an array of random numbers.
!-----------------------------------------------------------------------
!
!  Initialize output array.
!
      DO i=LBij,UBij
        R(i)=0.0_r8
      END DO
!
!  F90 intrinsic pseudorandom numbers with an uniform distribution
!  over the range 0 <= R < 1.
!
      IF (Rscheme.eq.0) THEN
        CALL random_number (R(LBij:))
!
!  Scale (fac=2*SQRT(3)) the random number to insure the expectation
!  mean to be approximately zero, E(R)=0, and the expectation variance
!  to be approximately unity, E(R^2)=1 (Bennett, 2002; book page 72).
!
!  Recall that,
!
!      E(R) = sum[R f(R)]
!
!  where F(R) is the random variable probability distribution.
!
        DO i=Imin,Imax
          R(i)=fac*(R(i)-0.5_r8)
        END DO
        Rmin=MINVAL(R)
        Rmax=MAXVAL(R)
!
!  Random deviates with a Gaussian (normal) distribuiton with zero
!  mean and unit variance.
!
      ELSE IF (Rscheme.eq.1) THEN
        IF (Master) THEN
          CALL gasdev (random)
        END IF
#ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, random)
#endif
        Rmin=MINVAL(random)
        Rmax=MAXVAL(random)

        Ioff=1-LBij
        DO i=Imin,Imax
          ic=i+Ioff
          R(i)=random(ic)
        END DO
#ifdef DISTRIBUTE
        CALL mp_exchange2d_bry (ng, tile, model, 1, boundary,           &
     &                          LBij, UBij,                             &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          R)
#endif
      END IF

      RETURN
      END SUBROUTINE white_noise2d_bry


#ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE white_noise3d (ng, model, gtype, Rscheme,              &
     &                          Imin, Imax, Jmin, Jmax,                 &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          Rmin, Rmax, R)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_scatter3d
# endif
      USE nrutil, ONLY : gasdev
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, gtype, Rscheme
      integer, intent(in) :: Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

      real(r8), intent(out) :: Rmin, Rmax

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: R(LBi:,LBj:,LBk:)
# else
      real(r8), intent(out) :: R(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Npts, i, ic, j, k

      real(r8), parameter :: fac = 2.0_r8 * 1.73205080756887720_r8

      real(r8), allocatable :: random(:)
!
!-----------------------------------------------------------------------
!  Generate an array of random numbers.
!-----------------------------------------------------------------------
!
!  Initialize output array.
!
      DO k=LBk,UBk
        DO j=LBj,UBj
          DO i=LBi,UBi
            R(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
!
!  F90 intrinsic pseudorandom numbers with an uniform distribution
!  over the range 0 <= R < 1.
!
      IF (Rscheme.eq.0) THEN
        CALL random_number (R(LBi:,LBj:,LBk:))
!
!  Scale (fac=2*SQRT(3)) the random number to insure the expectation
!  mean to be approximately zero, E(R)=0, and the expectation variance
!  to be approximately unity, E(R^2)=1.
!
!  Recall that,
!
!      E(R) = sum[R f(R)]
!
!  where F(R) is the random variable probability distribution
!
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              R(i,j,k)=fac*(R(i,j,k)-0.5_r8)
            END DO
          END DO
        END DO
        Rmin=MINVAL(R)
        Rmax=MAXVAL(R)
!
!  Random deviates with a Gaussian (normal) distribution with zero
!  mean and unit variance.
!
      ELSE IF (Rscheme.eq.1) THEN
        Npts=(Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)
        IF (.not.allocated(random)) THEN
          allocate ( random(Npts+2) )
        END IF
        IF (Master) THEN
          CALL gasdev (random)
          Rmin=random(1)
          Rmax=random(1)
          DO ic=1,Npts
            Rmin=MIN(Rmin,random(ic))
            Rmax=MAX(Rmax,random(ic))
          END DO
        END IF
# ifdef DISTRIBUTE
        Npts=SIZE(random)
        CALL mp_scatter3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk,     &
     &                     NghostPoints, gtype, Rmin, Rmax,             &
#  if defined READ_WATER && defined MASKING
     &                     (Lm(ng)+2)*(Mm(ng)+2),                       &
     &                     SCALARS(ng)%IJwater(:,gtype),                &
#  endif
     &                     Npts, random, R)
# else
        ic=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ic=ic+1
              R(i,j,k)=random(ic)
            END DO
          END DO
        END DO
# endif
        IF (allocated(random)) THEN
          deallocate (random)
        END IF
      END IF

      RETURN
      END SUBROUTINE white_noise3d

!
!***********************************************************************
      SUBROUTINE white_noise3d_bry (ng, tile, model, boundary, Rscheme, &
     &                              Imin, Imax,                         &
     &                              LBij, UBij, LBk, UBk,               &
     &                              Rmin, Rmax, R)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf
      USE mp_exchange_mod, ONLY : mp_exchange3d_bry
# endif
      USE nrutil, ONLY : gasdev
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, boundary, Rscheme
      integer, intent(in) :: Imin, Imax
      integer, intent(in) :: LBij, UBij, LBk, UBk

      real(r8), intent(out) :: Rmin, Rmax

# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: R(LBij:,LBk:)
# else
      real(r8), intent(out) :: R(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Ilen, Ioff, i, ic, k, kc

      real(r8), parameter :: fac = 2.0_r8 * 1.73205080756887720_r8

      real(r8), dimension((UBij-LBij+1)*(UBk-LBk+1)) :: random
!
!-----------------------------------------------------------------------
!  Generate an array of random numbers.
!-----------------------------------------------------------------------
!
!  Initialize output array.
!
      DO k=LBk,UBk
        DO i=LBij,UBij
          R(i,k)=0.0_r8
        END DO
      END DO
!
!  F90 intrinsic pseudorandom numbers with an uniform distribution
!  over the range 0 <= R < 1.
!
      IF (Rscheme.eq.0) THEN
        CALL random_number (R(LBij:,LBk:))
!
!  Scale (fac=2*SQRT(3)) the random number to insure the expectation
!  mean to be approximately zero, E(R)=0, and the expectation variance
!  to be approximately unity, E(R^2)=1 (Bennett, 2002; book page 72).
!
!  Recall that,
!
!      E(R) = sum[R f(R)]
!
!  where F(R) is the random variable probability distribution.
!
        DO k=LBk,UBk
          DO i=Imin,Imax
            R(i,k)=fac*(R(i,k)-0.5_r8)
          END DO
        END DO
        Rmin=MINVAL(R)
        Rmax=MAXVAL(R)
!
!  Random deviates with a Gaussian (normal) distribuiton with zero
!  mean and unit variance.
!
      ELSE IF (Rscheme.eq.1) THEN
        IF (Master) THEN
          CALL gasdev (random)
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, random)
# endif
        Rmin=MINVAL(random)
        Rmax=MAXVAL(random)

        Ilen=UBij-LBij+1
        Ioff=1-LBij

        DO k=LBk,UBk
          kc=(k-LBk)*Ilen
          DO i=Imin,Imax
            ic=i+Ioff+kc
            R(i,k)=random(ic)
          END DO
        END DO
# ifdef DISTRIBUTE
        CALL mp_exchange3d_bry (ng, tile, model, 1, boundary,           &
     &                          LBij, UBij, LBk, UBk,                   &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          R)
# endif
      END IF

      RETURN
      END SUBROUTINE white_noise3d_bry
#endif

      END MODULE white_noise_mod
