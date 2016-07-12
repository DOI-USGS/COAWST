#include "cppdefs.h"
      SUBROUTINE gasdev_s (harvest)
!
!svn $Id: gasdev.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Return in harvest a normally distributed deviate with zero mean     !
!  and unit variance, using RAN1 as the source of uniform deviates.    !
!                                                                      !
!  Scalar version adapted from Numerical Recipes.                      !
!                                                                      !
!  Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery,    !
!     1996:  Numerical Recipes in Fortran 90,  The Art of Parallel     !
!     Scientific Computing, 2nd Edition, Cambridge Univ. Press.        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE nrutil, ONLY : ran1
!
!  Imported variable declarations.
!
      real(r8), intent(out) :: harvest
!
!  Local variable declarations.
!
      logical, save :: gaus_stored = .FALSE.

      real(r8) :: rsq, v1, v2
      real(r8), save :: g
!
!-----------------------------------------------------------------------
!  Compute a normally distributed scalar deviate.
!-----------------------------------------------------------------------
!
!  We have an extra deviate handy, so return it, and unset the flag.
!
      IF (gaus_stored) THEN
        harvest=g
        gaus_stored=.FALSE.
!
!  We do not have an extra deviate handy, so pick two uniform numbers
!  in the square extending from -1 to +1 in each direction.
!
      ELSE
        DO
          CALL ran1 (v1)
          CALL ran1 (v2)
          v1=2.0_r8*v1-1.0_r8
          v2=2.0_r8*v2-1.0_r8
          rsq=v1*v1+v2*v2
!
!  See if they are in the unit circle, and if they are not, try again.
!
          IF ((rsq.gt.0.0_r8).and.(rsq.lt.1.0_r8)) EXIT
        END DO
!
!  Now make the Box-Muller transformation to get two normal deviates.
!  Return one and save the other for next time.
!
        rsq=SQRT(-2.0_r8*LOG(rsq)/rsq)
        harvest=v1*rsq
        g=v2*rsq
        gaus_stored=.TRUE.
      END IF

      RETURN
      END SUBROUTINE gasdev_s

      SUBROUTINE gasdev_v (harvest)
!
!=======================================================================
!                                                                      !
!  Return in harvest a normally distributed deviate with zero mean     !
!  and unit variance, using RAN1 as the source of uniform deviates.    !
!                                                                      !
!  Vector version adapted from Numerical Recipes.                      !
!                                                                      !
!  Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery,    !
!     1996:  Numerical Recipes in Fortran 90,  The Art of Parallel     !
!     Scientific Computing, 2nd Edition, Cambridge Univ. Press.        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE nrutil, ONLY : array_copy
      USE nrutil, ONLY : ran1
!
!  Imported variable declarations.
!
      real(r8), dimension(:), intent(out) :: harvest
!
!  Local variable declarations.
!
      logical, save :: gaus_stored = .FALSE.

      logical, dimension(SIZE(harvest)) :: mask

      integer(i4b), save :: last_allocated = 0
      integer(i4b) :: i, ii, m, mc, n, ng, nn

      real(r8), dimension(SIZE(harvest)) :: rsq, v1, v2, v3
      real(r8), allocatable, dimension(:), save :: g
!
!-----------------------------------------------------------------------
!  Compute a normally distributed vector deviate.
!-----------------------------------------------------------------------
!
!  We have an extra deviate handy, so return it, and unset the flag.
!
      n=SIZE(harvest)
      IF (n.ne.last_allocated) THEN
        IF (last_allocated.ne.0) DEALLOCATE (g)
        ALLOCATE ( g(n) )
        last_allocated=n
        gaus_stored=.FALSE.
      END IF
!
!  We do not have an extra deviate handy, so pick two uniform numbers
!  in the square extending from -1 to +1 in each direction.
!
      IF (gaus_stored) THEN
        harvest=g
        gaus_stored=.FALSE.
      ELSE
        ng=1
        DO
          IF (ng.gt.n) EXIT
          CALL ran1 (v1(ng:n))
          CALL ran1 (v2(ng:n))
          v1(ng:n)=2.0_r8*v1(ng:n)-1.0_r8
          v2(ng:n)=2.0_r8*v2(ng:n)-1.0_r8
!
!  See if they are in the unit circle, and if they are not, try again.
!  The original code was modified for portability with old F90 compilers
!  when using the PACK function (HGA/AMM).
!
          rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
          mask(ng:n)=((rsq(ng:n).gt.0.0_r8).and.(rsq(ng:n).lt.1.0_r8))
          mc=COUNT(mask(ng:n))

          ii=0
          DO i=ng,n
            IF (mask(i)) THEN
              ii=ii+1
              v3(ii)=v1(i)
            END IF
          END DO
          CALL array_copy (v3(1:mc), v1(ng:), nn, m)
          ii=ng-1
          DO i=ng,n
            IF (mask(i)) THEN
              ii=ii+1
              v2(ii)=v2(i)
              rsq(ii)=rsq(i)
            END IF
          END DO
          ng=ng+nn
        END DO
!
!  Make the Box-Muller transformation to get two normal deviates.
!  Return one and save the other for next time.
!
        rsq=SQRT(-2.0_r8*LOG(rsq)/rsq)
        harvest=v1*rsq
        g=v2*rsq
        gaus_stored=.TRUE.
      END IF

      RETURN
      END SUBROUTINE gasdev_v
