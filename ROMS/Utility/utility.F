#include "cppdefs.h"
      MODULE utility_mod
!
!svn $Id: utility.F 1336 2008-01-24 02:45:56Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several all purpuse generic routines:          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    nrng            NSWC Gaussian random number generator.            !
!    urng            NSWC uniform random number generator.             !
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PUBLIC

      CONTAINS

      SUBROUTINE nrng (ix, a, n, ierr)
!
!=======================================================================
!                                                                      !
!  Gaussian random-number generator from the NSWC Library. It calls    !
!  the NSWC uniform random-number generator, URNG.                     !
!                                                                      !
!  Modernised and included in ROMS by Mark Hadfield, NIWA.             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: n

      integer, intent(inout) :: ix

      integer, intent(out) :: ierr

#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: a(:)
#else
      real(r8), intent(out) :: a(n)
#endif
!
!  Local variable declarations.
!
      integer :: i, m

      real(r8), parameter ::pi2 = 6.2831853071796_r8

      real(r8) :: phi, r
      real(r8) :: temp(1)
!
!-----------------------------------------------------------------------
!  Generate Gaussian random numbers.
!-----------------------------------------------------------------------
!
      CALL urng (ix, a, n, ierr)
!
      IF (ierr.ne.0) RETURN
!
      IF (n.gt.1) THEN
        m=n/2
        m=m+m
        DO i=1,m,2
          r=SQRT(-2.0_r8*LOG(a(i)))
          phi=pi2*a(i+1)
          a(i  )=r*COS(phi)
          a(i+1)=r*SIN(phi)
        END DO
         IF (m.eq.n) RETURN
      END IF
!
      CALL urng (ix, temp, 1, ierr)
!
      r=SQRT(-2.0_r8*LOG(a(n)))
!
      a(n)=r*COS(pi2*temp(1))
!
      RETURN
      END SUBROUTINE nrng

      SUBROUTINE urng (ix, x, n, ierr)
!
!=======================================================================
!                                                                      !
!  Uniform random-number generator from the NSWC Library               !
!                                                                      !
!  Uses the recursion ix = ix*a mod p, where 0 < ix < p                !
!                                                                      !
!  Written by Linus Schrage, University of Chicago. Adapted for NSWC   !
!  Library by A. H. Morris. Modernised & included in ROMS by Mark      !
!  Hadfield, NIWA.                                                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: n

      integer, intent(inout) :: ix

      integer, intent(out) :: ierr

#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: x(:)
#else
      real(r8), intent(out) :: x(n)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: a = 16807          ! 7^5
      integer, parameter :: b15 = 32768        ! 2^15
      integer, parameter :: b16 = 65536        ! 2^16
      integer, parameter :: p = 2147483647     ! 2^31-1

      integer :: fhi, k, l, leftlo, xalo, xhi

      real(r8), parameter :: s = 0.465661E-09_r8
!
!-----------------------------------------------------------------------
!  Generate random numbers.
!-----------------------------------------------------------------------
!
      IF (n.le.0) THEN
        ierr=1
        RETURN
      END IF
      IF ((ix.le.0).or.(ix.ge.p)) THEN
        ierr=2
        RETURN
      END IF
!
      ierr=0
!
      DO l=1,n
!
! Get 15 high order bits of "ix".
!
        xhi=ix/b16
!
! Get 16 lower bits of ix and multiply with "a".
!
        xalo=(ix-xhi*b16)*a
!
! Get 15 high order bits of the product.
!
        leftlo=xalo/b16
!
! Form the 31 highest bits of "a*ix".
!
        fhi=xhi*a+leftlo
!
! Obtain the overflow past the 31st bit of "a*ix".
!
        k=fhi/b15
!
! Assemble all the parts and presubtract "p". The parentheses are
! essential.
!
        ix=(((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
!
! Add "p" if necessary.
!
        IF (ix.lt.0) ix=ix+p
!
! Rescale "ix", to interpret it as a value between 0 and 1.
! the scale factor "s" is selected to be as near "1/p" as is
! appropriate in order that the floating value for "ix = 1",
! namely "s", be roughly the same distance from 0 as "(p-1)*s"
! is from 1. The current value for "s" assures us that "x(l)"
! is less than 1 for any floating point arithmetic of 6
! or more digits.
!
         x(l)=REAL(ix,r8)*s
      END DO
      RETURN
      END SUBROUTINE urng

      END MODULE utility_mod
