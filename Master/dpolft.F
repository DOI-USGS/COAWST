#include "cppdefs.h"
      MODULE dpolft_mod
#ifdef UV_BANIHASHEMI
!

!svn $Id: dpolft_mod.F 830 2024-04-24 21:21:11Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!    Submitted by John C. Warner                                       !
!=======================================================================
!                                                                      !
!  MODULE CONTAINS:                                                    !
!  Based on Supplement to the Paper 'Approximation of Wave Action Flux !
!  Velocity in Strongly Sheared Mean Flows',                           !
!  Banihashemi, Kirby, and Dong. (2017). Ocean Modelling 116, 33-47.   !
!                                                                      !
!  This function was obtained from:                                    !
!  https://netlib.org/slatec/src/                                      !
!  https://netlib.org/slatec/src/dpopt.f                               !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: DPOLFT, DPCOEF, DP1VLU
      CONTAINS

!     *DECK DPOLFT
!***********************************************************************
      SUBROUTINE DPOLFT (N, X, Y, W, MAXDEG, NDEG, EPS, R, IERR, A)
!***********************************************************************
!** 
!C** https://www.netlib.org/slatec/src/dpolft.f
!C** 
!C***BEGIN PROLOGUE  DPOLFT
!C***PURPOSE  Fit discrete data in a least squares sense by polynomials
!C            in one variable.
!C***LIBRARY   SLATEC
!C***CATEGORY  K1A1A2
!C***TYPE      DOUBLE PRECISION (POLFIT-S, DPOLFT-D)
!C***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
!C***AUTHOR  Shampine, L. F., (SNLA)
!C           Davenport, S. M., (SNLA)
!C           Huddleston, R. E., (SNLL)
!C***DESCRIPTION
!C
!C     Abstract
!C
!C     Given a collection of points X(I) and a set of values Y(I) which
!C     correspond to some function or measurement at each of the X(I),
!C     subroutine  DPOLFT  computes the weighted least-squares polynomial
!C     fits of all degrees up to some degree either specified by the user
!C     or determined by the routine.  The fits thus obtained are in
!C    orthogonal polynomial form.  Subroutine  DP1VLU  may then be
!C  called to evaluate the fitted polynomials and any of their
!C  derivatives at any point.  The subroutine  DPCOEF  may be used to
!C  express the polynomial fits as powers of (X-C) for any specified
!C  point C.
!C
!C  The parameters for  DPOLFT  are
!C
!C  Input -- All TYPE REAL variables are DOUBLE PRECISION
!C      N -      the number of data points.  The arrays X, Y and W
!C                  must be dimensioned at least  N  (N .GE. 1).
!C         X -      array of values of the independent variable.  These
!C                  values may appear in any order and need not all be
!C                  distinct.
!C         Y -      array of corresponding function values.
!C         W -      array of positive values to be used as weights.  If
!C                  W(1) is negative,  DPOLFT  will set all the weights
!C                  to 1.0, which means unweighted least squares error
!C                  will be minimized.  To minimize relative error, the
!C                  user should set the weights to:  W(I) = 1.0/Y(I)**2,
!C                  I = 1,...,N .
!C         MAXDEG - maximum degree to be allowed for polynomial fit.
!C                  MAXDEG  may be any non-negative integer less than  N.
!C                  Note -- MAXDEG  cannot be equal to  N-1  when a
!C                  statistical test is to be used for degree selection,
!C                  i.e., when input value of  EPS  is negative.
!C         EPS -    specifies the criterion to be used in determining
!C                  the degree of fit to be computed.
!C                  (1)  If  EPS  is input negative,  DPOLFT  chooses the
!C                       degree based on a statistical F test of
!C                       significance.  One of three possible
!C                       significance levels will be used:  .01, .05 or
!C                       .10.  If  EPS=-1.0 , the routine will
!C                       automatically select one of these levels based
!C                       on the number of data points and the maximum
!C                       degree to be considered.  If  EPS  is input as
!C                       -.01, -.05, or -.10, a significance level of
!C                       .01, .05, or .10, respectively, will be used.
!C                  (2)  If  EPS  is set to 0.,  DPOLFT  computes the
!C                       polynomials of degrees 0 through  MAXDEG .
!C                  (3)  If  EPS  is input positive,  EPS  is the RMS
!C                       error tolerance which must be satisfied by the
!C                       fitted polynomial.  DPOLFT  will increase the
!C                       degree of fit until this criterion is met or
!C                       until the maximum degree is reached.
!C
!C     Output -- All TYPE REAL variables are DOUBLE PRECISION
!C         NDEG -   degree of the highest degree fit computed.
!C         EPS -    RMS error of the polynomial of degree  NDEG .
!C         R -      vector of dimension at least NDEG containing values
!C                  of the fit of degree  NDEG  at each of the  X(I) .
!C                  Except when the statistical test is used, these
!C                  values are more accurate than results from subroutine
!C                  DP1VLU  normally are.
!C         IERR -   error flag with the following possible values.
!C             1 -- indicates normal execution, i.e., either
!C                  (1)  the input value of  EPS  was negative, and the
!C                       computed polynomial fit of degree  NDEG
!C                       satisfies the specified F test, or
!C                  (2)  the input value of  EPS  was 0., and the fits of
!C                       all degrees up to  MAXDEG  are complete, or
!C                  (3)  the input value of  EPS  was positive, and the
!C                       polynomial of degree  NDEG  satisfies the RMS
!C                       error requirement.
!C             2 -- invalid input parameter.  At least one of the input
!C                  parameters has an illegal value and must be corrected
!C                  before  DPOLFT  can proceed.  Valid input results
!C                  when the following restrictions are observed
!C                       N .GE. 1
!C                       0 .LE. MAXDEG .LE. N-1  for  EPS .GE. 0.
!C                       0 .LE. MAXDEG .LE. N-2  for  EPS .LT. 0.
!C                       W(1)=-1.0  or  W(I) .GT. 0., I=1,...,N .
!C             3 -- cannot satisfy the RMS error requirement with a
!C                  polynomial of degree no greater than  MAXDEG .  Best
!C                  fit found is of degree  MAXDEG .
!C             4 -- cannot satisfy the test for significance using
!C                  current value of  MAXDEG .  Statistically, the
!C                  best fit found is of order  NORD .  (In this case,
!C                  NDEG will have one of the values:  MAXDEG-2,
!C                  MAXDEG-1, or MAXDEG).  Using a higher value of
!C                  MAXDEG  may result in passing the test.
!C         A -      work and output array having at least 3N+3MAXDEG+3
!C                  locations
!C
!C     Note - DPOLFT  calculates all fits of degrees up to and including
!C            NDEG .  Any or all of these fits can be evaluated or
!C            expressed as powers of (X-C) using  DP1VLU  and  DPCOEF
!C            after just one call to  DPOLFT .
!C
!C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!C                 Curve fitting by polynomials in one variable, Report
!C                 SLA-74-0270, Sandia Laboratories, June 1974.
!C***ROUTINES CALLED  DP1VLU, XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   740601  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   891006  Cosmetic changes to prologue.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900911  Added variable YP to DOUBLE PRECISION declaration.  (WRB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C   920527  Corrected erroneous statements in DESCRIPTION.  (WRB)
!C***END PROLOGUE  DPOLFT
      INTEGER I,IDEGF,IERR,J,JP1,JPAS,K1,K1PJ,K2,K2PJ,K3,K3PI,K4,       &
     & K4PI,K5,K5PI,KSIG,M,MAXDEG,MOP1,NDEG,NDER,NFAIL
!jcw 
      INTEGER N,NDO
      DOUBLE PRECISION TEMD1,TEMD2
      DOUBLE PRECISION A(*),DEGF,DEN,EPS,ETST,F,FCRIT,R(*),SIG,SIGJ,    &
     & SIGJM1,SIGPAS,TEMP,X(*),XM,Y(*),W(*),W1,W11
! jcw made YP(:)
      DOUBLE PRECISION, ALLOCATABLE :: YP(:)
      DOUBLE PRECISION CO(4,3)
      SAVE CO
      DATA  CO(1,1), CO(2,1), CO(3,1), CO(4,1), CO(1,2), CO(2,2),       &
     &      CO(3,2), CO(4,2), CO(1,3), CO(2,3), CO(3,3),                &
     &  CO(4,3)/-13.086850D0,-2.4648165D0,-3.3846535D0,-1.2973162D0,    &
     &          -3.3381146D0,-1.7812271D0,-3.2578406D0,-1.6589279D0,    &
     &          -1.6282703D0,-1.3152745D0,-3.2640179D0,-1.9829776D0/
!C***FIRST EXECUTABLE STATEMENT  DPOLFT
      M = ABS(N)
      IF (M .EQ. 0) GO TO 30
      IF (MAXDEG .LT. 0) GO TO 30
      A(1) = MAXDEG
      MOP1 = MAXDEG + 1
      IF (M .LT. MOP1) GO TO 30
      IF (EPS .LT. 0.0D0 .AND.  M .EQ. MOP1) GO TO 30
      XM = M
      ETST = EPS*EPS*XM
      IF (W(1) .LT. 0.0D0) GO TO 2
      DO 1 I = 1,M
        IF (W(I) .LE. 0.0D0) GO TO 30
 1      CONTINUE
      GO TO 4
 2    DO I = 1,M
        W(I) = 1.0d0
      END DO
 4    IF (EPS .GE. 0.0D0) GO TO 8
!C
!C DETERMINE SIGNIFICANCE LEVEL INDEX TO BE USED IN STATISTICAL TEST FOR
!C CHOOSING DEGREE OF POLYNOMIAL FIT
!C
      IF (EPS .GT. (-.55D0)) GO TO 5
      IDEGF = M - MAXDEG - 1
      KSIG = 1
      IF (IDEGF .LT. 10) KSIG = 2
      IF (IDEGF .LT. 5) KSIG = 3
      GO TO 8
 5    KSIG = 1
      IF (EPS .LT. (-.03D0)) KSIG = 2
      IF (EPS .LT. (-.07D0)) KSIG = 3
!C
!C INITIALIZE INDEXES AND COEFFICIENTS FOR FITTING
!C
 8    K1 = MAXDEG + 1
      K2 = K1 + MAXDEG
      K3 = K2 + MAXDEG + 2
      K4 = K3 + M
      K5 = K4 + M
      DO I = 2,K4
        A(I) = 0.0D0
      END DO
      W11 = 0.0D0
      IF (N .LT. 0) GO TO 11
!C
!C UNCONSTRAINED CASE
!C
      DO I = 1,M
        K4PI = K4 + I
        A(K4PI) = 1.0D0
        W11 = W11 + W(I)
      END DO
      GO TO 13
!C
!C CONSTRAINED CASE
!C
 11   DO I = 1,M
        K4PI = K4 + I
        W11 = W11 + W(I)*A(K4PI)**2
      END DO
!C
!C COMPUTE FIT OF DEGREE ZERO
!C
 13   TEMD1 = 0.0D0
      DO I = 1,M
        K4PI = K4 + I
        TEMD1 = TEMD1 + W(I)*Y(I)*A(K4PI)
      END DO
      TEMD1 = TEMD1/W11
      A(K2+1) = TEMD1
      SIGJ = 0.0D0
      DO I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = TEMD1*A(K4PI)
        R(I) = TEMD2
        A(K5PI) = TEMD2 - R(I)
        SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
      END DO
      J = 0
!C
!C SEE IF POLYNOMIAL OF DEGREE 0 SATISFIES THE DEGREE SELECTION CRITERION
!C
!jcw  IF (EPS) 24,26,27
      IF (EPS .LT. 0.) THEN
        GOTO  24
      ELSE IF (EPS .EQ. 0.) THEN
        GOTO 26
      ELSE
        GOTO 27
      END IF
!C
!C INCREMENT DEGREE
!C
 16   J = J + 1
      JP1 = J + 1
      K1PJ = K1 + J
      K2PJ = K2 + J
      SIGJM1 = SIGJ
!C
!C COMPUTE NEW B COEFFICIENT EXCEPT WHEN J = 1
!C
      IF (J .GT. 1) A(K1PJ) = W11/W1
!C
!C COMPUTE NEW A COEFFICIENT
!C
      TEMD1 = 0.0D0
      DO I = 1,M
        K4PI = K4 + I
        TEMD2 = A(K4PI)
        TEMD1 = TEMD1 + X(I)*W(I)*TEMD2*TEMD2
      END DO
      A(JP1) = TEMD1/W11
!C
!C EVALUATE ORTHOGONAL POLYNOMIAL AT DATA POINTS
!C
      W1 = W11
      W11 = 0.0D0
      DO I = 1,M
        K3PI = K3 + I
        K4PI = K4 + I
        TEMP = A(K3PI)
        A(K3PI) = A(K4PI)
        A(K4PI) = (X(I)-A(JP1))*A(K3PI) - A(K1PJ)*TEMP
        W11 = W11 + W(I)*A(K4PI)**2
      END DO
!C
!C GET NEW ORTHOGONAL POLYNOMIAL COEFFICIENT USING PARTIAL DOUBLE
!C PRECISION
!C
      TEMD1 = 0.0D0
      DO I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = W(I)*((Y(I)-R(I))-A(K5PI))*A(K4PI)
        TEMD1 = TEMD1 + TEMD2
      END DO
      TEMD1 = TEMD1/W11
      A(K2PJ+1) = TEMD1
!C
!C UPDATE POLYNOMIAL EVALUATIONS AT EACH OF THE DATA POINTS, AND
!C ACCUMULATE SUM OF SQUARES OF ERRORS.  THE POLYNOMIAL EVALUATIONS ARE
!C COMPUTED AND STORED IN EXTENDED PRECISION.  FOR THE I-TH DATA POINT,
!C THE MOST SIGNIFICANT BITS ARE STORED IN  R(I) , AND THE LEAST
!C SIGNIFICANT BITS ARE IN  A(K5PI) .
!C
      SIGJ = 0.0D0
      DO I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = R(I) + A(K5PI) + TEMD1*A(K4PI)
        R(I) = TEMD2
        A(K5PI) = TEMD2 - R(I)
        SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
      END DO
!C
!C SEE IF DEGREE SELECTION CRITERION HAS BEEN SATISFIED OR IF DEGREE
!C MAXDEG  HAS BEEN REACHED
!C
!jcw  IF (EPS) 23,26,27
      IF (EPS .LT. 0.) THEN
        GOTO  23
      ELSE IF (EPS .EQ. 0.) THEN
        GOTO 26
      ELSE
        GOTO 27
      END IF
!C
!C COMPUTE F STATISTICS  (INPUT EPS .LT. 0.)
!C
 23   IF (SIGJ .EQ. 0.0D0) GO TO 29
      DEGF = M - J - 1
      DEN = (CO(4,KSIG)*DEGF + 1.0D0)*DEGF
      FCRIT = (((CO(3,KSIG)*DEGF) + CO(2,KSIG))*DEGF + CO(1,KSIG))/DEN
      FCRIT = FCRIT*FCRIT
      F = (SIGJM1 - SIGJ)*DEGF/SIGJ
      IF (F .LT. FCRIT) GO TO 25
!C
!C POLYNOMIAL OF DEGREE J SATISFIES F TEST
!C
 24   SIGPAS = SIGJ
      JPAS = J
      NFAIL = 0
      IF (MAXDEG .EQ. J) GO TO 32
      GO TO 16
!C
!C POLYNOMIAL OF DEGREE J FAILS F TEST.  IF THERE HAVE BEEN THREE
!C SUCCESSIVE FAILURES, A STATISTICALLY BEST DEGREE HAS BEEN FOUND.
!C
 25   NFAIL = NFAIL + 1
      IF (NFAIL .GE. 3) GO TO 29
      IF (MAXDEG .EQ. J) GO TO 32
      GO TO 16
!C
!C RAISE THE DEGREE IF DEGREE  MAXDEG  HAS NOT YET BEEN REACHED  (INPUT
!C EPS = 0.)
!C
 26   IF (MAXDEG .EQ. J) GO TO 28
      GO TO 16
!C
!C SEE IF RMS ERROR CRITERION IS SATISFIED  (INPUT EPS .GT. 0.)
!C
 27   IF (SIGJ .LE. ETST) GO TO 28
      IF (MAXDEG .EQ. J) GO TO 31
      GO TO 16
!C
!C RETURNS
!C
 28   IERR = 1
      NDEG = J
      SIG = SIGJ
      GO TO 33
 29   IERR = 1
      NDEG = JPAS
      SIG = SIGPAS
      GO TO 33
 30   IERR = 2
!     CALL XERMSG ('SLATEC', 'DPOLFT', 'INVALID INPUT PARAMETER.', 2,   &
!    &   1)
      WRITE(*,*) 'SLATEC DPOLFT INVALID INPUT PARAMETER.'
      GO TO 37
 31   IERR = 3
      NDEG = MAXDEG
      SIG = SIGJ
      GO TO 33
 32   IERR = 4
      NDEG = JPAS
      SIG = SIGPAS
!C
 33   A(K3) = NDEG
!C
!C WHEN STATISTICAL TEST HAS BEEN USED, EVALUATE THE BEST POLYNOMIAL AT
!C ALL THE DATA POINTS IF  R  DOES NOT ALREADY CONTAIN THESE VALUES
!C
      IF(EPS .GE. 0.0  .OR.  NDEG .EQ. MAXDEG) GO TO 36
      NDER = 0
!jcw
      NDO = MAX(NDER,0)
      NDO = MIN(NDO,NDEG)
      ALLOCATE (YP(NDO))
      DO I = 1,M
        CALL DP1VLU (NDEG,NDER,X(I),R(I),YP,A)
      END DO
      DEALLOCATE (YP)
!36   EPS = SQRT(SIG/XM)
!36   EPS = SQRT(SIG/XM+0.001D0)
 36   CONTINUE
 37   RETURN
      END
      
!*DECK DPCOEF
!***********************************************************************
      SUBROUTINE DPCOEF (L, C, TC, A)
!***********************************************************************
!C***BEGIN PROLOGUE  DPCOEF
!C***PURPOSE  Convert the DPOLFT coefficients to Taylor series form.
!C***LIBRARY   SLATEC
!C***CATEGORY  K1A1A2
!C***TYPE      DOUBLE PRECISION (PCOEF-S, DPCOEF-D)
!C***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
!C***AUTHOR  Shampine, L. F., (SNLA)
!C           Davenport, S. M., (SNLA)
!C***DESCRIPTION
!C
!C     Abstract
!C
!C     DPOLFT  computes the least squares polynomial fit of degree  L  as
!C     a sum of orthogonal polynomials.  DPCOEF  changes this fit to its
!C     Taylor expansion about any point  C , i.e. writes the polynomial
!C     as a sum of powers of (X-C).  Taking  C=0.  gives the polynomial
!C     in powers of X, but a suitable non-zero  C  often leads to
!C     polynomials which are better scaled and more accurately evaluated.
!C
!C     The parameters for  DPCOEF  are
!C
!C     INPUT -- All TYPE REAL variables are DOUBLE PRECISION
!C         L -      Indicates the degree of polynomial to be changed to
!C                  its Taylor expansion.  To obtain the Taylor
!C                  coefficients in reverse order, input  L  as the
!C                  negative of the degree desired.  The absolute value
!C                  of L  must be less than or equal to NDEG, the highest
!C                  degree polynomial fitted by  DPOLFT .
!C         C -      The point about which the Taylor expansion is to be
!C                  made.
!C         A -      Work and output array containing values from last
!C                  call to  DPOLFT .
!C
!C     OUTPUT -- All TYPE REAL variables are DOUBLE PRECISION
!C         TC -     Vector containing the first LL+1 Taylor coefficients
!C                  where LL=ABS(L).  If  L.GT.0 , the coefficients are
!C                  in the usual Taylor series order, i.e.
!C                    P(X) = TC(1) + TC(2)*(X-C) + ... + TC(N+1)*(X-C)**N
!C                  If L .LT. 0, the coefficients are in reverse order,
!C                  i.e.
!C                    P(X) = TC(1)*(X-C)**N + ... + TC(N)*(X-C) + TC(N+1)
!C
!C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!C                 Curve fitting by polynomials in one variable, Report
!C                 SLA-74-0270, Sandia Laboratories, June 1974.
!C***ROUTINES CALLED  DP1VLU
!C***REVISION HISTORY  (YYMMDD)
!C   740601  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   891006  Cosmetic changes to prologue.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  DPCOEF
!C
      INTEGER I,L,LL,LLP1,LLP2,NEW,NR
      DOUBLE PRECISION A(*),C,FAC,SAVE,TC(*)
!C***FIRST EXECUTABLE STATEMENT  DPCOEF
      LL = ABS(L)
      LLP1 = LL + 1
      CALL DP1VLU (LL,LL,C,TC(1),TC(2),A)
      IF (LL .LT. 2) GO TO 2
      FAC = 1.0d0
      DO I = 3,LLP1
        FAC = FAC*(I-1)
        TC(I) = TC(I)/FAC
      END DO
 2    IF (L .GE. 0) GO TO 4
      NR = LLP1/2
      LLP2 = LL + 2
      DO I = 1,NR
        SAVE = TC(I)
        NEW = LLP2 - I
        TC(I) = TC(NEW)
        TC(NEW) = SAVE
      END DO
 4    RETURN
      END
      
!*DECK DP1VLU
!***********************************************************************
      SUBROUTINE DP1VLU (L, NDER, X, YFIT, YP, A)
!***********************************************************************
!C***BEGIN PROLOGUE  DP1VLU
!C***PURPOSE  Use the coefficients generated by DPOLFT to evaluate the
!C            polynomial fit of degree L, along with the first NDER of
!C            its derivatives, at a specified point.
!C***LIBRARY   SLATEC
!C***CATEGORY  K6
!C***TYPE      DOUBLE PRECISION (PVALUE-S, DP1VLU-D)
!C***KEYWORDS  CURVE FITTING, LEAST SQUARES, POLYNOMIAL APPROXIMATION
!C***AUTHOR  Shampine, L. F., (SNLA)
!C           Davenport, S. M., (SNLA)
!C***DESCRIPTION
!C
!C     Abstract
!C
!C     The subroutine  DP1VLU  uses the coefficients generated by  DPOLFT
!C     to evaluate the polynomial fit of degree  L , along with the first
!C     NDER  of its derivatives, at a specified point.  Computationally
!C     stable recurrence relations are used to perform this task.
!C
!C     The parameters for  DP1VLU  are
!C
!C     Input -- ALL TYPE REAL variables are DOUBLE PRECISION
!C         L -      the degree of polynomial to be evaluated.  L  may be
!C                  any non-negative integer which is less than or equal
!C                  to  NDEG , the highest degree polynomial provided
!C                  by  DPOLFT .
!C         NDER -   the number of derivatives to be evaluated.  NDER
!C                  may be 0 or any positive value.  If NDER is less
!C                  than 0, it will be treated as 0.
!C         X -      the argument at which the polynomial and its
!C                  derivatives are to be evaluated.
!C         A -      work and output array containing values from last
!C                  call to  DPOLFT .
!C
!C     Output -- ALL TYPE REAL variables are DOUBLE PRECISION
!C         YFIT -   value of the fitting polynomial of degree  L  at  X
!C         YP -     array containing the first through  NDER  derivatives
!C                  of the polynomial of degree  L .  YP  must be
!C                  dimensioned at least  NDER  in the calling program.
!C
!C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!C                 Curve fitting by polynomials in one variable, Report
!C                 SLA-74-0270, Sandia Laboratories, June 1974.
!C***ROUTINES CALLED  XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   740601  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890911  Removed unnecessary intrinsics.  (WRB)
!C   891006  Cosmetic changes to prologue.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  DP1VLU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I,IC,ILO,IN,INP1,IUP,K1,K1I,K2,K3,K3P1,K3PN,K4,K4P1,K4PN,    &
     & KC,L,LM1,LP1,MAXORD,N,NDER,NDO,NDP1,NORD
      DOUBLE PRECISION A(*),CC,DIF,VAL,X,YFIT,YP(*)
      CHARACTER*8 XERN1, XERN2
!C***FIRST EXECUTABLE STATEMENT  DP1VLU
      IF (L .LT. 0) GO TO 12
      NDO = MAX(NDER,0)
      NDO = MIN(NDO,L)
      MAXORD = A(1) + 0.5D0
      K1 = MAXORD + 1
      K2 = K1 + MAXORD
      K3 = K2 + MAXORD + 2
      NORD = A(K3) + 0.5D0
      IF (L .GT. NORD) GO TO 11
      K4 = K3 + L + 1
      IF (NDER .LT. 1) GO TO 2
      DO I = 1,NDER
        YP(I) = 0.0D0
      END DO
 2    IF (L .GE. 2) GO TO 4
      IF (L .EQ. 1) GO TO 3
!C
!C L IS 0
!C
      VAL = A(K2+1)
      GO TO 10
!C
!C L IS 1
!C
 3    CC = A(K2+2)
      VAL = A(K2+1) + (X-A(2))*CC
      IF (NDER .GE. 1) YP(1) = CC
      GO TO 10
!C
!C L IS GREATER THAN 1
!C
 4    NDP1 = NDO + 1
      K3P1 = K3 + 1
      K4P1 = K4 + 1
      LP1 = L + 1
      LM1 = L - 1
      ILO = K3 + 3
      IUP = K4 + NDP1
      DO I = ILO,IUP
        A(I) = 0.0D0
      END DO
      DIF = X - A(LP1)
      KC = K2 + LP1
      A(K4P1) = A(KC)
      A(K3P1) = A(KC-1) + DIF*A(K4P1)
      A(K3+2) = A(K4P1)
!C
!C EVALUATE RECURRENCE RELATIONS FOR FUNCTION VALUE AND DERIVATIVES
!C
      DO I = 1,LM1
        IN = L - I
        INP1 = IN + 1
        K1I = K1 + INP1
        IC = K2 + IN
        DIF = X - A(INP1)
        VAL = A(IC) + DIF*A(K3P1) - A(K1I)*A(K4P1)
        IF (NDO .LE. 0) GO TO 8
        DO N = 1,NDO
          K3PN = K3P1 + N
          K4PN = K4P1 + N
          YP(N) = DIF*A(K3PN) + N*A(K3PN-1) - A(K1I)*A(K4PN)
        END DO
!C
!C SAVE VALUES NEEDED FOR NEXT EVALUATION OF RECURRENCE RELATIONS
!C
        DO N = 1,NDO
          K3PN = K3P1 + N
          K4PN = K4P1 + N
          A(K4PN) = A(K3PN)
          A(K3PN) = YP(N)
        END DO
 8      A(K4P1) = A(K3P1)
        A(K3P1) = VAL
      END DO
!C
!C NORMAL RETURN OR ABORT DUE TO ERROR
!C
 10   YFIT = VAL
      RETURN
!C
   11 WRITE (XERN1, '(I8)') L
      WRITE (XERN2, '(I8)') NORD
!     CALL XERMSG ('SLATEC', 'DP1VLU',                                  &
      WRITE(*,*)  'SLATEC  DP1VLU                                       &
     &    THE ORDER OF POLYNOMIAL EVALUATION,                           &
     &    REQUESTED EXCEEDS THE HIGHEST ORDER FIT, NORD =               &
     &    COMPUTED BY DPOLFT -- EXECUTION TERMINATED.'
      RETURN
!C
!  12 CALL XERMSG ('SLATEC', 'DP1VLU',                                  &
   12 CONTINUE
      WRITE(*,*)      'SLATEC DP1VLU                                    &
     &    INVALID INPUT PARAMETER.  ORDER OF POLYNOMIAL EVALUATION      &
     &    REQUESTED IS NEGATIVE.'
      RETURN
      END

#endif
      END MODULE dpolft_mod

      