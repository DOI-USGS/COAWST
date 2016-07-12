#include "cppdefs.h"
      MODULE erf_mod
!
!svn $Id: erf.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines compute the incomplete gamma function and error      !
!  function:                                                           !
!                                                                      !
!  ERF        Error function, ERF(x)                                   !
!  ERFC       Complementary error function, ERFC(x)                    !
!  ERFCC      Complementary error function, ERFCC(x): cheaper          !
!               Chebyshev fitting approximation.                       !
!  GAMMP      Incomplete gamma function, P(a,x)                        !
!  GAMMQ      Incomplete gamma function complement, Q(a,x)=1-P(a,x)    !
!                                                                      !
!  Adapted from Numerical Recipes:                                     !
!                                                                      !
!    Press, W.H, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery,   !
!     1986:  Numerical Recipes in Fortran 77,  The Art of Parallel     !
!     Scientific Computing, 1st Edition, Cambridge Univ. Press.        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PRIVATE
      PUBLIC :: ERF
      PUBLIC :: ERFC
      PUBLIC :: ERFCC
      PUBLIC :: GAMMP
      PUBLIC :: GAMMQ
!
      CONTAINS
!
      FUNCTION ERF (x) RESULT (value)
!
!=======================================================================
!                                                                      !
!  This routine computes the error function, ERF(x).                   !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: x
!
!  Local variable declarations.
!
      real(r8) :: value
!
!-----------------------------------------------------------------------
!  Compute error function in terms of the incomplete gamma function,
!  P(a,x).
!-----------------------------------------------------------------------
!
      IF (x.lt.0.0_r8) THEN
        value=-GAMMP(0.5_r8, x**2)
      ELSE
        value= GAMMP(0.5_r8, x**2)
      ENDIF

      RETURN
      END FUNCTION ERF
!
      FUNCTION ERFC (x) RESULT (value)
!
!=======================================================================
!                                                                      !
!  This routine computes the complementary error function, ERFC(x).    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: x
!
!  Local variable declarations.
!
      real(r8) :: value
!
!-----------------------------------------------------------------------
!  Compute complementary error function in terms of the incomplete
!  gamma functions, P(a,x) and Q(a,x).
!-----------------------------------------------------------------------
!
      IF (x.lt.0.0_r8) THEN
        value=1.0_r8+GAMMP(0.5_r8, x**2)
      ELSE
        value=GAMMQ(0.5_r8, x**2)
      endif

      RETURN
      END FUNCTION ERFC
!
      FUNCTION ERFCC (x) RESULT (value)
!
!=======================================================================
!                                                                      !
!  This routine computes the complementary error function, ERFCC(x),   !
!  with fractional error everywhere less than 1.2E-7.                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: x
!
!  Local variable declarations.
!
      real(r8) :: t, value, z
!
!-----------------------------------------------------------------------
!  Compute complementary error function based on Chebyshev fitting to
!  an inspired guess.
!-----------------------------------------------------------------------
!
      z=ABS(x)
      t=1.0_r8/(1.0_r8+0.5_r8*z)
      value=t*EXP(-z*z-1.26551223_r8+                                   &
     &            t*(1.00002368_r8+                                     &
     &               t*(0.37409196_r8+                                  &
     &                  t*(0.09678418_r8+                               &
     &                     t*(-0.18628806_r8+                           &
     &                        t*(0.27886807_r8+                         &
     &                           t*(-1.13520398_r8+                     &
     &                              t*(1.48851587_r8+                   &
     &                                 t*(-0.82215223_r8+               &
     &                                    t*.17087277)))))))))
      IF (x.lt.0.0_r8) value=2.0-value

      RETURN
      END FUNCTION ERFCC
!
      FUNCTION GAMMLN (xx) RESULT (value)
!
!=======================================================================
!                                                                      !
!  This routine computes the natural logarithm of the gamma function,  !
!  LN(Gamma(xx)). The natural logarithm is computed to avoid floating- !
!  point overflow.                                                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: xx
!
!  Local variable declarations.
!
      integer :: j

      real(r8) :: ser, tmp, x, y
      real(r8) :: value
      real(r8), save, dimension(6) :: cof =                             &
     &          (/ 76.18009172947146_r8,                                &
     &            -86.50532032941677_r8,                                &
     &             24.01409824083091_r8,                                &
     &             -1.231739572450155_r8,                               &
     &              0.1208650973866179E-2_r8,                           &
     &             -0.5395239384953E-5_r8 /)
      real(r8), save :: stp = 2.5066282746310005_r8
!
!-----------------------------------------------------------------------
!  Compute natural logarithm of the gamma function.
!-----------------------------------------------------------------------
!
      x=xx
      y=x
      tmp=x+5.5_r8
      tmp=(x+0.5_r8)*LOG(tmp)-tmp
      ser=1.000000000190015_r8
      DO j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      END DO
      value=tmp+LOG(stp*ser/x)

      RETURN
      END FUNCTION GAMMLN
!
      FUNCTION GAMMP (a,x) RESULT (value)
!
!=======================================================================
!                                                                      !
!  This routine computes the incomplete gamma function, P(a,x).        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_iounits
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: a, x
!
!  Local variable declarations.
!
      real(r8) :: gammcf, gamser, gln
      real(r8) :: value
!
!-----------------------------------------------------------------------
!  Compute incompleate gamma fucntion, P(a.x).
!-----------------------------------------------------------------------
!
      IF ((x.lt.0.0_r8).or.(a.le.0.0_r8)) THEN
        WRITE (stdout,10) a, x
  10    FORMAT (/,' GAMMAP - gamma function negative argument,',        &
     &          ' a = ',1pe13.6,'  x = ',1pe13.6)
        exit_flag=8
      END IF

      IF (x.lt.a+1.0_r8) THEN
        CALL gser (gamser, a, x, gln)
        value=gamser                 ! series representation
      ELSE
        CALL gcf (gammcf, a, x, gln)
        value=1.0_r8-gammcf          ! continued fraction representation
      ENDIF

      RETURN
      END FUNCTION GAMMP
!
      FUNCTION GAMMQ (a,x) RESULT (value)
!
!=======================================================================
!                                                                      !
!  This routine computes the incomplete gamma function complement,     !
!  Q(a,x)=1-P(a,x).                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_iounits
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: a, x
!
!  Local variable declarations.
!
      real(r8) :: gammcf, gamser, gln
      real(r8) :: value
!
!-----------------------------------------------------------------------
!  Compute incompleate gamma fucntion, Q(a,x)=1-P(a.x).
!-----------------------------------------------------------------------
!
      IF ((x.lt.0.0_r8).or.(a.le.0.0_r8)) THEN
        WRITE (stdout,10) a, x
  10    FORMAT (/,' GAMMAQ - gamma function negative argument,',        &
     &          ' a = ',1pe13.6,'  x = ',1pe13.6)
        exit_flag=8
      END IF

      IF (x.lt.a+1.0_r8) THEN
        CALL gser (gamser, a, x, gln)
        value=1.0_r8-gamser          ! series representation
      ELSE
        CALL gcf (gammcf, a, x, gln)
        value=gammcf                 ! continued fraction representation
      ENDIF

      RETURN
      END FUNCTION GAMMQ
!
      SUBROUTINE gser (gamser, a, x, gln)
!
!=======================================================================
!                                                                      !
!  This routine computes the incomplete gamma function P(a,x)          !
!  evaliuated by its series representation as gamser.                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_iounits
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: a, x
      real(r8), intent(out) :: gamser, gln
!
!  Local variable declarations.
!
      logical :: Converged
      integer, parameter :: ITMAX = 100
      integer :: i

      real(r8), parameter :: EPS =3.0E-7_r8
      real(r8) :: ap, del, my_sum
!
!-----------------------------------------------------------------------
!  Compue incomplete gamma function by its series representation.
!-----------------------------------------------------------------------
!
      gln=GAMMLN(a)
      IF (x.le.0.0_r8) THEN
        IF (x.lt.0.0) THEN
          WRITE (stdout,10) x
  10      FORMAT (/,' GSER - gamma function negative argument, x = ',   &
     &            1pe13.6)
          exit_flag=8
        END IF
        gamser=0.0_r8
        RETURN
      END IF
      ap=a
      my_sum=1.0_r8/a
      del=my_sum
      Converged=.FALSE.
      DO i=1,ITMAX
        ap=ap+1.0_r8
        del=del*x/ap
        my_sum=my_sum+del
        IF (ABS(del).lt.ABS(my_sum)*EPS) THEN
          Converged=.TRUE.
          EXIT
        END IF
      END DO
      IF (Converged) THEN
        gamser=my_sum*EXP(-x+a*LOG(x)-gln)
      ELSE
        WRITE (stdout,20) ITMAX
  20    FORMAT (/,' GSER - Gamma function not converged, ITMAX = ',     &
     &          i4.4,/,8x,'a is too large')
        exit_flag=8
      END IF

      RETURN
      END SUBROUTINE gser
!
      SUBROUTINE gcf (gammcf, a, x, gln)
!
!=======================================================================
!                                                                      !
!  This routine computes the incomplete gamma function P(a,x)          !
!  evaliuated by its series representation as gamser.                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_iounits
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: a, x
      real(r8), intent(out) :: gammcf, gln
!
!  Local variable declarations.
!
      logical :: Converged

      integer, parameter :: ITMAX = 100          ! number of iterations
      integer :: i

      real(r8), parameter :: EPS = 3.0E-7_r8     ! relative accuracy
      real(r8), parameter :: FPMIN = 1.0E-30_r8  ! smallest number
      real(r8) :: an, b, c, d, del, h
!
!-----------------------------------------------------------------------
!  Compue incomplete gamma function evaluated by its continued fraction
!  representation.
!-----------------------------------------------------------------------
!
      gln=GAMMLN(a)
      b=x+1.0_r8-a
      c=1.0_r8/FPMIN
      d=1.0_r8/b
      h=d
      Converged=.FALSE.
      DO i=1,ITMAX
        an=-REAL(i,r8)*(REAL(i,r8)-a)
        b=b+2.0_r8
        d=an*d+b
        IF (ABS(d).lt.FPMIN) d=FPMIN
        c=b+an/c
        IF (ABS(c).lt.FPMIN) c=FPMIN
        d=1.0_r8/d
        del=d*c
        h=h*del
        IF (ABS(del-1.0_r8).lt.EPS) THEN
          Converged=.TRUE.
          EXIT
        END IF
      END DO
      IF (Converged) THEN
        gammcf=EXP(-x+a*LOG(x)-gln)*h
      ELSE
        WRITE (stdout,10) ITMAX
  10    FORMAT (/,' GCF - Gamma function not converged, ITMAX = ',      &
     &          i4.4,/,8x,'a is too large')
        exit_flag=8
      END IF

      RETURN
      END SUBROUTINE gcf
!
      FUNCTION AERF (arg) RESULT (value)
!
!=======================================================================
!                                                                      !
!  The error function has been taken in approximate form from          !
!  Abramowitz and Stegun, page 299.                                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Abramowitz, M. and I.A. Stegun, 1972: Handbook of Mathematical      !
!       Functions with Formulas, Graphs, and Mathematical Tables,      !
!       National Bureau of Standards, Applied Mathematics Series       !
!       55, US Goverment Printing Office, Washington, D.C.             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: arg
!
!  Local variable declarations.
!
      real(r8), parameter :: a1 = +0.254829592_r8
      real(r8), parameter :: a2 = -0.284496736_r8
      real(r8), parameter :: a3 = -1.421413741_r8
      real(r8), parameter :: a4 = -1.453152027_r8
      real(r8), parameter :: a5 = +1.061405429_r8
      real(r8), parameter :: p  = +0.3275911_r8

      real(r8) :: big, c, my_sign, t, y
      real(r8) :: value
!
!-----------------------------------------------------------------------
!  Compute approximated error function.
!-----------------------------------------------------------------------
!
      big=20.0_r8
      my_sign=1.0_r8

      IF (arg.lt.0.0_r8) THEN
        my_sign=-1.0
      ELSE
        my_sign=1.0_r8
      END IF

      IF (arg.eq.0.0_r8) THEN
        value=0.0_r8
        RETURN
      END IF
!
!  Polynomial approximation.
!
      y=ABS(arg)
      t=1.0_r8/(1.0_r8+p*y)
      c=t*(a1+t*(a2+t*(a3+t*(a4+a5*t))))
      IF (y*y.lt.big) THEN
        value=my_sign*(1.0_r8-c*EXP(-y*y))
      ELSE
        value=my_sign
      END IF

      RETURN
      END FUNCTION AERF

      END MODULE erf_mod
