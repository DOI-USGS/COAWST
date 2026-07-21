#include "cppdefs.h"
      MODULE lapack_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module includes modernized versions of the selected routines   !
!  from thw Linear Algebra Package (LAPACK) library, which are used    !
!  in ROMS 4D-Var algorithms.                                          !
!                                                                      !
!  Adapted from LAPACK library version 2.0                             !
!                                                                      !
!  NOTES:                                                              !
!                                                                      !
!   - The LAPACK library was written originally in Fortran-77 and has  !
!     been extensively tested and used in various compilers and        !
!     applications. We are modernizing a few of the functions used by  !
!     ROMS at the request of NOAA and NCEP to remove the undesirable   !
!     GOTOs statements with the modern capabilities of the Fortran     !
!     standard (1995 and 2003).                                        !
!                                                                      !
!   - We need to declare the non-scalar arguments to routines as       !
!     assumed-size array, A(*), instead of assumed-shape array, A(:).  !
!     It is relatively common in linear algebra to pass sub-matrices   !
!     as an argument. However, this feature was declared obsolescent   !
!     in Fortran95. The colon we need an explicit interface. Currently,!
!     if fails to compile if colon are used intead of asterisk, even   !
!     if it is inside this module.                                     !
!                                                                      !
!======================================================================!
!
      USE mod_kinds
!
!      implicit none
!
      PUBLIC  :: DSTEQR
!
      PRIVATE :: DLAE2
      PRIVATE :: DLAEV2
      PRIVATE :: DLAMC1
      PRIVATE :: DLAMC2
      PRIVATE :: DLAMC3
      PRIVATE :: DLAMC4
      PRIVATE :: DLAMC5
      PRIVATE :: DLAMCH
      PRIVATE :: DLANST
      PRIVATE :: DLAPY2
      PRIVATE :: DLARTG
      PRIVATE :: DLASCL
      PRIVATE :: DLASET
      PRIVATE :: DLASR
      PRIVATE :: DLASRT
      PRIVATE :: DLASSQ
      PRIVATE :: LSAME
      PRIVATE :: XERBLA
!
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!
      SUBROUTINE DSTEQR ( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!=======================================================================
!                                                                      !
!  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a  !
!  symmetric tridiagonal matrix using the implicit QL or QR method.    !
!  The eigenvectors of a full or band symmetric matrix can also be     !
!  found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this    !
!  matrix to tridiagonal form.                                         !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  COMPZ   (input) CHARACTER*1                                         !
!          = 'N':  Compute eigenvalues only.                           !
!          = 'V':  Compute eigenvalues and eigenvectors of the         !
!                  original symmetric matrix. On entry, Z must contain !
!                  the orthogonal matrix used to reduce the original   !
!                  matrix to tridiagonal form.                         !
!          = 'I':  Compute eigenvalues and eigenvectors of the         !
!                  tridiagonal matrix.  Z is initialized to the        !
!                  identity matrix.                                    !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The order of the matrix.  N >= 0.                           !
!                                                                      !
!  D       (input/output) DOUBLE PRECISION array, dimension (N)        !
!          On entry, the diagonal elements of the tridiagonal matrix.  !
!          On exit, if INFO = 0, the eigenvalues in ascending order.   !
!                                                                      !
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)      !
!          On entry, the (n-1) subdiagonal elements of the tridiagonal !
!          matrix.                                                     !
!          On exit, E has been destroyed.                              !
!                                                                      !
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)   !
!          On entry, if  COMPZ = 'V', then Z contains the orthogonal   !
!          matrix used in the reduction to tridiagonal form.           !
!          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the  !
!          orthonormal eigenvectors of the original symmetric matrix,  !
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors !
!          of the symmetric tridiagonal matrix.                        !
!          If COMPZ = 'N', then Z is not referenced.                   !
!                                                                      !
!  LDZ     (input) INTEGER                                             !
!          The leading dimension of the array Z.  LDZ >= 1, and if     !
!          eigenvectors are desired, then  LDZ >= max(1,N).            !
!                                                                      !
!  WORK    (workspace) DOUBLE PRECISION array, dimension               !
!          (max(1,2*N-2)).                                             !
!          If COMPZ = 'N', then WORK is not referenced.                !
!                                                                      !
!  INFO    (output) INTEGER                                            !
!          = 0:  successful exit                                       !
!          < 0:  if INFO = -i, the i-th argument had an illegal value  !
!          > 0:  the algorithm has failed to find all the eigenvalues  !
!                in a total of 30*N iterations; if INFO = i, then i    !
!                elements of E have not converged to zero; on exit, D  !
!                and E contain the elements of a symmetric tridiagonal !
!                matrix which is orthogonally similar to the original  !
!                matrix.                                               !
!                                                                      !
!  -- LAPACK routine (version 2.0) --                                  !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)  :: LDZ, N
      integer, intent(out) :: INFO
!
      real (dp), intent(inout) :: WORK(*)
      real (dp), intent(inout) :: D(*), E(*), Z(LDZ,*)
!
      character (len=1), intent(in) :: COMPZ
!
!  Local variable declarations.
!
      integer, parameter :: MAXIT = 30
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: THREE = 3.0_dp
!
      integer :: I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, LENDM1
      integer :: LENDP1, LENDSV, LM1, LSV, M, MM, MM1, NM1, NMAXIT
!
      real (dp) :: ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2
      real (dp) :: S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input parameters.
!
      INFO = 0
!
      IF (LSAME(COMPZ, 'N')) THEN
        ICOMPZ = 0
      ELSE IF (LSAME(COMPZ, 'V')) THEN
        ICOMPZ = 1
      ELSE IF (LSAME(COMPZ, 'I' )) THEN
        ICOMPZ = 2
      ELSE
        ICOMPZ = -1
      END IF
!
      IF (ICOMPZ.LT.0) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF ((LDZ.LT.1).OR.(ICOMPZ.GT.0 .AND. LDZ.LT.MAX(1,N))) THEN
        INFO = -6
      END IF
!
      IF (INFO.NE.0) THEN
        CALL XERBLA( 'DSTEQR', -INFO )
        RETURN
      END IF
!
!  Quick return if possible
!
      IF (N.EQ.0) RETURN
!
      IF ( N.EQ.1 ) THEN
        IF (ICOMPZ.EQ.2) Z(1,1) = ONE
        RETURN
      END IF
!
!  Determine the unit roundoff and over/underflow thresholds.
!
      EPS = DLAMCH('E')
      EPS2 = EPS**2
      SAFMIN = DLAMCH('S')
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT(SAFMAX) / THREE
      SSFMIN = SQRT(SAFMIN) / EPS2
!
!  Compute the eigenvalues and eigenvectors of the tridiagonal
!  matrix.
!
      IF (ICOMPZ .EQ. 2) THEN
        CALL DLASET ('Full', N, N, ZERO, ONE, Z, LDZ)
      END IF
!
      NMAXIT = N*MAXIT
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
      IF (L1 .GT. N) GO TO 160
      IF (L1 .GT. 1) E(L1-1) = ZERO
      IF (L1 .LE. NM1) THEN
        DO M = L1, NM1
          TST = ABS(E(M))
          IF (TST .EQ. ZERO) GO TO 30
          IF (TST .LE. (SQRT(ABS(D(M)))*SQRT(ABS(D(M+1))))*EPS) THEN
            E(M) = ZERO
            GO TO 30
          END IF
        END DO
      END IF
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF (LEND .EQ. L) GO TO 10
!
!  Scale submatrix in rows and columns L to LEND.
!
      ANORM = DLANST('I', LEND-L+1, D(L), E(L))
      ISCALE = 0
      IF (ANORM .EQ. ZERO) GO TO 10
      IF (ANORM .GT. SSFMAX) THEN
        ISCALE = 1
        CALL DLASCL ('G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D(L), N,    &
     &               INFO )
        CALL DLASCL ('G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E(L), N,      &
     &               INFO )
      ELSE IF (ANORM .LT. SSFMIN) THEN
        ISCALE = 2
        CALL DLASCL ('G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D(L), N,    &
     &               INFO )
        CALL DLASCL ('G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,    &
     &               INFO )
      END IF
!
!  Choose between QL and QR iteration
!
      IF (ABS(D(LEND)) .LT. ABS(D(L))) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF (LEND .GT. L) THEN
!
!  QL Iteration.
!
!    Look for small subdiagonal element.
!
   40   CONTINUE
        IF (L .NE. LEND) THEN
          LENDM1 = LEND - 1
          DO M = L, LENDM1
            TST = ABS(E(M))**2
            IF (TST .LE. (EPS2*ABS(D(M)))*ABS(D(M+1))+SAFMIN) GO TO 60
          END DO
        END IF
!
        M = LEND
!
   60   CONTINUE
        IF (M .LT. LEND) E(M) = ZERO
        P = D(L)
        IF (M .EQ. L) GO TO 80
!
!    If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!    to compute its eigensystem.
!
        IF (M .EQ. L+1) THEN
          IF (ICOMPZ .GT. 0) THEN
            CALL DLAEV2 ( D(L), E(L), D(L+1), RT1, RT2, C, S)
            WORK(L) = C
            WORK(N-1+L) = S
            CALL DLASR ('R', 'V', 'B', N, 2, WORK(L),                   &
     &                  WORK(N-1+L), Z(1, L), LDZ)
          ELSE
            CALL DLAE2 (D(L), E(L), D(L+1), RT1, RT2)
          END IF
          D(L) = RT1
          D(L+1) = RT2
          E(L) = ZERO
          L = L + 2
          IF (L .LE. LEND) GO TO 40
          GO TO 140
        END IF
!
        IF (JTOT .EQ. NMAXIT) GO TO 140
        JTOT = JTOT + 1
!
!    Form shift.
!
        G = (D(L+1)-P) / (TWO*E(L))
        R = DLAPY2(G, ONE)
        G = D(M) - P + (E(L) / (G + SIGN(R, G)))
!
        S = ONE
        C = ONE
        P = ZERO
!
!    Inner loop
!
        MM1 = M - 1
        DO I = MM1, L, -1
          F = S*E(I)
          B = C*E(I)
          CALL DLARTG (G, F, C, S, R)
          IF (I .NE. M-1) E(I+1) = R
          G = D(I+1) - P
          R = (D(I)-G)*S + TWO*C*B
          P = S*R
          D(I+1) = G+P
          G = C*R - B
!
!    If eigenvectors are desired, then save rotations.
!
          IF (ICOMPZ.GT.0) THEN
            WORK(I) = C
            WORK(N-1+I) = -S
          END IF
        END DO
!
!    If eigenvectors are desired, then apply saved rotations.
!
        IF (ICOMPZ.GT.0) THEN
          MM = M - L + 1
          CALL DLASR ('R', 'V', 'B', N, MM, WORK(L), WORK(N-1+L),       &
     &                Z(1,L), LDZ)
        END IF
!
        D(L) = D( L ) - P
        E(L) = G
        GO TO 40
!
!    Eigenvalue found.
!
   80   CONTINUE
        D(L) = P
 !
        L = L + 1
        IF (L .LE. LEND) GO TO 40
        GO TO 140
!
      ELSE
!
!  QR Iteration
!
!    Look for small superdiagonal element.
!
   90   CONTINUE
        IF (L .NE. LEND) THEN
          LENDP1 = LEND + 1
          DO M = L, LENDP1, -1
            TST = ABS(E( M-1))**2
            IF (TST .LE. (EPS2*ABS(D(M)))*ABS(D(M-1))+SAFMIN) GO TO 110
          END DO
        END IF
!
        M = LEND
!
  110   CONTINUE
        IF (M .GT. LEND) E(M-1) = ZERO
        P = D(L)
        IF (M .EQ. L) GO TO 130
!
!    If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!    to compute its eigensystem.
!
        IF (M .EQ. L-1) THEN
          IF (ICOMPZ .GT. 0) THEN
            CALL DLAEV2 (D(L-1), E(L-1), D(L), RT1, RT2, C, S)
            WORK(M) = C
            WORK(N-1+M) = S
            CALL DLASR ('R', 'V', 'F', N, 2, WORK(M),                   &
     &                  WORK(N-1+M), Z(1,L-1), LDZ)
          ELSE
            CALL DLAE2 (D(L-1), E(L-1), D(L), RT1, RT2)
          END IF
          D(L-1) = RT1
          D(L) = RT2
          E(L-1) = ZERO
          L = L - 2
          IF (L .GE. LEND) GO TO 90
          GO TO 140
        END IF
!
        IF (JTOT .EQ. NMAXIT) GO TO 140
        JTOT = JTOT + 1
!
!    Form shift.
!
        G = (D(L-1)-P) / (TWO*E(L-1))
        R = DLAPY2(G, ONE)
        G = D(M) - P + (E(L-1) / (G+SIGN(R, G)))
!
        S = ONE
        C = ONE
        P = ZERO
!
!    Inner loop
!
        LM1 = L - 1
        DO I = M, LM1
          F = S*E(I)
          B = C*E(I)
          CALL DLARTG (G, F, C, S, R)
          IF (I .NE. M) E(I-1) = R
          G = D(I) - P
          R = (D(I+1)-G)*S + TWO*C*B
          P = S*R
          D(I) = G + P
          G = C*R - B
!
!    If eigenvectors are desired, then save rotations.
!
          IF (ICOMPZ.GT.0) THEN
            WORK(I) = C
            WORK(N-1+I) = S
          END IF
        END DO
!
!    If eigenvectors are desired, then apply saved rotations.
!
        IF (ICOMPZ.GT.0) THEN
          MM = L - M + 1
          CALL DLASR ('R', 'V', 'F', N, MM, WORK(M), WORK(N-1+M),       &
     &                Z(1,M), LDZ)
        END IF
!
        D(L) = D(L) - P
        E(LM1) = G
        GO TO 90
!
!    Eigenvalue found.
!
  130   CONTINUE
        D(L) = P
!
        L = L - 1
        IF (L .GE. LEND) GO TO 90
        GO TO 140
!
      END IF
!
!  Undo scaling if necessary.
!
  140 CONTINUE
      IF (ISCALE .EQ. 1) THEN
        CALL DLASCL ('G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,         &
     &               D(LSV), N, INFO)
        CALL DLASCL ('G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E(LSV),   &
     &               N, INFO)
      ELSE IF (ISCALE .EQ. 2) THEN
        CALL DLASCL ('G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,         &
     &               D(LSV), N, INFO)
        CALL DLASCL ('G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E(LSV),   &
     &               N, INFO)
      END IF
!
!  Check for no convergence to an eigenvalue after a total
!  of N*MAXIT iterations.
!
      IF (JTOT .LT. NMAXIT) GO TO 10
      DO I = 1, N - 1
        IF (E(I) .NE. ZERO) INFO = INFO + 1
      END DO
      GO TO 190
!
!  Order eigenvalues and eigenvectors.
!
  160 CONTINUE
      IF (ICOMPZ .EQ. 0) THEN
!
!    Use Quick Sort.
!
        CALL DLASRT ('I', N, D, INFO)
!
      ELSE
!
!    Use Selection Sort to minimize swaps of eigenvectors.
!
        DO II = 2, N
          I = II - 1
          K = I
          P = D( I )
          DO  J = II, N
            IF (D(J) .LT. P) THEN
               K = J
               P = D(J)
            END IF
          END DO
          IF (K .NE. I) THEN
            D(K) = D(I)
            D(I) = P
            CALL DSWAP (N, Z(1,I), 1, Z(1,K), 1)
          END IF
        END DO
      END IF
!
  190 CONTINUE
      RETURN
!
      END SUBROUTINE DSTEQR
!
      SUBROUTINE DLAE2 (A, B, C, RT1, RT2)
!
!=======================================================================
!                                                                      !
!  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix:       !
!                                                                      !
!     [  A   B  ]                                                      !
!     [  B   C  ].                                                     !
!                                                                      !
!  On return, RT1 is the eigenvalue of larger absolute value, and RT2  !
!  is the eigenvalue of smaller absolute value.                        !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  A       (input) DOUBLE PRECISION                                    !
!          The (1,1) element of the 2-by-2 matrix.                     !
!                                                                      !
!  B       (input) DOUBLE PRECISION                                    !
!          The (1,2) and (2,1) elements of the 2-by-2 matrix.          !
!                                                                      !
!  C       (input) DOUBLE PRECISION                                    !
!          The (2,2) element of the 2-by-2 matrix.                     !
!                                                                      !
!  RT1     (output) DOUBLE PRECISION                                   !
!          The eigenvalue of larger absolute value.                    !
!                                                                      !
!  RT2     (output) DOUBLE PRECISION                                   !
!          The eigenvalue of smaller absolute value.                   !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  RT1 is accurate to a few ulps barring over/underflow.               !
!                                                                      !
!  RT2 may be inaccurate if there is massive cancellation in the       !
!  determinant A*C-B*B; higher precision or correctly rounded or       !
!  correctly truncated arithmetic would be needed to compute RT2       !
!  accurately in all cases.                                            !
!                                                                      !
!  Overflow is possible only if RT1 is within a factor of 5 of         !
!  overflow.                                                           !
!                                                                      !
!  Underflow is harmless if the input data is 0 or exceeds             !
!  underflow_threshold / macheps.                                      !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in ) :: A, B, C
      real (dp), intent(out) :: RT1, RT2
!
!  Local variable declarations.
!
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      real (dp) :: AB, ACMN, ACMX, ADF, DF, RT, SM, TB
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Compute the eigenvalues.
!
      SM = A + C
      DF = A - C
      ADF = ABS(DF)
      TB = B + B
      AB = ABS(TB)
!
      IF (ABS(A) .GT. ABS(C)) THEN
        ACMX = A
        ACMN = C
      ELSE
        ACMX = C
        ACMN = A
      END IF
!
      IF (ADF .GT. AB) THEN
        RT = ADF*SQRT(ONE+(AB/ADF )**2)
      ELSE IF (ADF .LT. AB) THEN
        RT = AB*SQRT(ONE+(ADF/AB )**2 )
      ELSE
        RT = AB*SQRT(TWO)                 ! Includes case AB=ADF=0
      END IF
!
      IF (SM .LT. ZERO) THEN
        RT1 = HALF*(SM-RT)
!
!  Order of execution important.
!  To get fully accurate smaller eigenvalue,
!  next line needs to be executed in higher precision.
!
        RT2 = (ACMX / RT1)*ACMN - (B / RT1)*B
      ELSE IF (SM .GT. ZERO ) THEN
        RT1 = HALF*(SM + RT)
!
!  Order of execution important.
!  To get fully accurate smaller eigenvalue,
!  next line needs to be executed in higher precision.
!
        RT2 = (ACMX / RT1)*ACMN - (B / RT1)*B
      ELSE
        RT1 = HALF*RT                     ! Includes case RT1 = RT2 = 0
        RT2 = -HALF*RT
      END IF
!
      RETURN
      END SUBROUTINE DLAE2
!
      SUBROUTINE DLAEV2 (A, B, C, RT1, RT2, CS1, SN1)
!
!=======================================================================
!                                                                      !
!  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric        !
!  matrix:                                                             !
!                                                                      !
!     [  A   B  ]                                                      !
!     [  B   C  ].                                                     !
!                                                                      !
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is   !
!  the eigenvalue of smaller absolute value, and (CS1,SN1) is the unit !
!  right eigenvector for RT1, giving the decomposition:                !
!                                                                      !
!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]            !
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].           !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  A       (input) DOUBLE PRECISION                                    !
!          The (1,1) element of the 2-by-2 matrix.                     !
!                                                                      !
!  B       (input) DOUBLE PRECISION                                    !
!          The (1,2) element and the conjugate of the (2,1) element of !
!          the 2-by-2 matrix.                                          !
!                                                                      !
!  C       (input) DOUBLE PRECISION                                    !
!          The (2,2) element of the 2-by-2 matrix.                     !
!                                                                      !
!  RT1     (output) DOUBLE PRECISION                                   !
!          The eigenvalue of larger absolute value.                    !
!                                                                      !
!  RT2     (output) DOUBLE PRECISION                                   !
!          The eigenvalue of smaller absolute value.                   !
!                                                                      !
!  CS1     (output) DOUBLE PRECISION                                   !
!  SN1     (output) DOUBLE PRECISION                                   !
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.  !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  RT1 is accurate to a few ulps barring over/underflow.               !
!                                                                      !
!  RT2 may be inaccurate if there is massive cancellation in the       !
!  determinant A*C-B*B; higher precision or correctly rounded or       !
!  correctly truncated arithmetic would be needed to compute RT2       !
!  accurately in all cases.                                            !
!                                                                      !
!  CS1 and SN1 are accurate to a few ulps barring over/underflow.      !
!                                                                      !
!  Overflow is possible only if RT1 is within a factor of 5 of         !
!  overflow.                                                           !
!                                                                      !
!  Underflow is harmless if the input data is 0 or exceeds             !
!  underflow_threshold / macheps.                                      !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in ) :: A, B, C
      real (dp), intent(out) :: CS1, RT1, RT2, SN1
!
!  Local variable declarations.
!
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      integer :: SGN1, SGN2
      real (dp) :: AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, TB, TN
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Compute the eigenvalues.
!
      SM = A + C
      DF = A - C
      ADF = ABS(DF)
      TB = B + B
      AB = ABS(TB)
!
      IF (ABS(A) .GT. ABS(C)) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
!
      IF (ADF .GT. AB) THEN
        RT = ADF*SQRT(ONE+(AB/ADF)**2)
      ELSE IF (ADF .LT. AB) THEN
        RT = AB*SQRT(ONE+(ADF/AB)**2 )
      ELSE
        RT = AB*SQRT(TWO)                 ! Includes case AB=ADF=0
      END IF
      IF (SM.LT.ZERO) THEN
        RT1 = HALF*(SM - RT)
        SGN1 = -1
!
!  Order of execution important.
!  To get fully accurate smaller eigenvalue,
!  next line needs to be executed in higher precision.
!
        RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE IF (SM .GT. ZERO) THEN
        RT1 = HALF*(SM + RT)
        SGN1 = 1
!
!  Order of execution important.
!  To get fully accurate smaller eigenvalue,
!  next line needs to be executed in higher precision.
!
        RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE
        RT1 = HALF*RT                     ! Includes case RT1 = RT2 = 0
        RT2 = -HALF*RT
        SGN1 = 1
      END IF
!
!  Compute the eigenvector.
!
      IF (DF .GE. ZERO) THEN
        CS = DF + RT
        SGN2 = 1
      ELSE
        CS = DF - RT
        SGN2 = -1
      END IF
      ACS = ABS(CS)
      IF (ACS .GT. AB) THEN
        CT = -TB / CS
        SN1 = ONE / SQRT(ONE+CT*CT)
        CS1 = CT*SN1
      ELSE
        IF (AB .EQ. ZERO) THEN
          CS1 = ONE
          SN1 = ZERO
        ELSE
          TN = -CS / TB
          CS1 = ONE / SQRT(ONE+TN*TN)
          SN1 = TN*CS1
        END IF
      END IF
      IF (SGN1 .EQ. SGN2) THEN
        TN = CS1
        CS1 = -SN1
        SN1 = TN
      END IF
!
      RETURN
      END SUBROUTINE DLAEV2
!
      SUBROUTINE DLAMC1 (BETA, T, RND, IEEE1)
!
!=======================================================================
!                                                                      !
!  DLAMC1 determines the machine parameters given by BETA, T, RND, and !
!  IEEE1.                                                              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  BETA    (output) INTEGER                                            !
!          The base of the machine.                                    !
!                                                                      !
!  T       (output) INTEGER                                            !
!          The number of ( BETA ) digits in the mantissa.              !
!                                                                      !
!  RND     (output) LOGICAL                                            !
!          Specifies whether proper rounding  ( RND = .TRUE. )  or     !
!          chopping  ( RND = .FALSE. )  occurs in addition. This may   !
!          not be a reliable guide to the way in which the machine     !
!          performs its arithmetic.                                    !
!                                                                      !
!  IEEE1   (output) LOGICAL                                            !
!          Specifies whether rounding appears to be done in the IEEE   !
!          'round to nearest' style.                                   !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  The routine is based on the routine  ENVRON  by Malcolm and         !
!  incorporates suggestions by Gentleman and Marovich. See             !
!                                                                      !
!     Malcolm M. A. (1972) Algorithms to reveal properties of          !
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.    !
!                                                                      !
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms     !
!        that reveal properties of floating point arithmetic units.    !
!        Comms. of the ACM, 17, 276-277.                               !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(out) :: IEEE1, RND
!
      integer, intent(out) :: BETA, T
!
!  Local variable declarations.
!
      logical, save :: FIRST = .TRUE.
      logical, save :: LIEEE1, LRND
!
      integer, save :: LBETA, LT
!
      real (dp) :: A, B, C, F, ONE, QTR, SAVEC, T1, T2
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      FIRST_PASS : IF (FIRST) THEN
        FIRST = .FALSE.
        ONE = 1
!
!  LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!  IEEE1, T and RND.
!
!  Throughout this routine  we use the function  DLAMC3  to ensure
!  that relevant values are  stored and not held in registers,  or
!  are not affected by optimizers.
!
!  Compute  a = 2.0**m  with the  smallest positive integer m such that
!
!    fl( a + 1.0 ) = a.
!
        A = 1
        C = 1
!
        DO WHILE (C .EQ. ONE)
          A = 2*A
          C = DLAMC3( A, ONE )
          C = DLAMC3( C, -A )
        END DO
!
!  Now compute  b = 2.0**m  with the smallest positive integer m
!  such that
!
!    fl( a + b ) .gt. a.
!
        B = 1
        C = DLAMC3(A, B)
!
        DO WHILE (C .EQ. A)
          B = 2*B
          C = DLAMC3(A, B)
        END DO
!
!  Now compute the base.  a and c  are neighbouring floating point
!  numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!  their difference is beta. Adding 0.25 to c is to ensure that it
!  is truncated to beta and not ( beta - 1 ).
!
        QTR = ONE / 4
        SAVEC = C
        C = DLAMC3(C, -A)
        LBETA = C + QTR
!
!  Now determine whether rounding or chopping occurs,  by adding a
!  bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
        B = LBETA
        F = DLAMC3(B/2, -B/100)
        C = DLAMC3(F, A)
        IF (C .EQ. A) THEN
          LRND = .TRUE.
        ELSE
          LRND = .FALSE.
        END IF
        F = DLAMC3(B/2, B/100)
        C = DLAMC3(F, A)
        IF ((LRND) .AND. (C.EQ.A)) LRND = .FALSE.
!
!  Try and decide whether rounding is done in the  IEEE  round to
!  nearest style. B/2 is half a unit in the last place of the two
!  numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!  zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!  A, but adding B/2 to SAVEC should change SAVEC.
!
        T1 = DLAMC3(B/2, A)
        T2 = DLAMC3(B/2, SAVEC)
        LIEEE1 = (T1.EQ.A) .AND. (T2.GT.SAVEC) .AND. LRND
!
!  Now find  the  mantissa, t.  It should  be the  integer part of
!  log to the base beta of a,  however it is safer to determine  t
!  by powering.  So we find t as the smallest positive integer for
!  which
!
!    fl( beta**t + 1.0 ) = 1.0.
!
        LT = 0
        A = 1
        C = 1
!
        DO WHILE (C .EQ. ONE)
          LT = LT + 1
          A = A*LBETA
          C = DLAMC3(A, ONE)
          C = DLAMC3(C, -A)
        END DO
!
      END IF FIRST_PASS
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
!
      RETURN
      END SUBROUTINE DLAMC1
!
      SUBROUTINE DLAMC2 (BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX)
!
!=======================================================================
!                                                                      !
!  DLAMC2 determines the machine parameters specified in its argument  !
!  list.                                                               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  BETA    (output) INTEGER                                            !
!          The base of the machine.                                    !
!                                                                      !
!  T       (output) INTEGER                                            !
!          The number of ( BETA ) digits in the mantissa.              !
!                                                                      !
!  RND     (output) LOGICAL                                            !
!          Specifies whether proper rounding  ( RND = .TRUE. )  or     !
!          chopping  ( RND = .FALSE. )  occurs in addition. This may   !
!          not be a reliable guide to the way in which the machine     !
!          performs its arithmetic.                                    !
!                                                                      !
!  EPS     (output) DOUBLE PRECISION                                   !
!          The smallest positive number such that                      !
!                                                                      !
!             fl( 1.0 - EPS ) .LT. 1.0,                                !
!                                                                      !
!          where fl denotes the computed value.                        !
!                                                                      !
!  EMIN    (output) INTEGER                                            !
!          The minimum exponent before (gradual) underflow occurs.     !
!                                                                      !
!  RMIN    (output) DOUBLE PRECISION                                   !
!          The smallest normalized number for the machine, given by    !
!          BASE**( EMIN - 1 ), where  BASE  is the floating point      !
!          value of BETA.                                              !
!                                                                      !
!  EMAX    (output) INTEGER                                            !
!          The maximum exponent before overflow occurs.                !
!                                                                      !
!  RMAX    (output) DOUBLE PRECISION                                   !
!          The largest positive number for the machine, given by       !
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating      !
!          point value of BETA.                                        !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  The computation of  EPS  is based on a routine PARANOIA by          !
!  W. Kahan of the University of California at Berkeley.               !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(out)   :: RND
!
      integer, intent(out)   :: BETA, EMAX, EMIN, T
!
      real (dp), intent(out) :: EPS, RMAX, RMIN
!
!  Local variable declarations.
!
      logical, save :: FIRST = .TRUE.
      logical, save :: IWARN = .FALSE.
      logical       :: IEEE, LIEEE1, LRND
!
      integer, save :: LBETA, LEMAX, LEMIN, LT
      integer       :: GNMIN, GPMIN, I, NGNMIN, NGPMIN
!
      real (dp), save :: LEPS, LRMAX, LRMIN
      real (dp)       :: A, B, C, HALF, ONE, RBASE,                     &
     &                   SIXTH, SMALL, THIRD, TWO, ZERO
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      FIRST_PASS : IF (FIRST) THEN
        FIRST = .FALSE.
        ZERO = 0
        ONE = 1
        TWO = 2
!
!  LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!  BETA, T, RND, EPS, EMIN and RMIN.
!
!  Throughout this routine  we use the function  DLAMC3  to ensure
!  that relevant values are stored  and not held in registers,  or
!  are not affected by optimizers.
!
!  DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
        CALL DLAMC1 (LBETA, LT, LRND, LIEEE1)
!
!  Start to find EPS.
!
        B = LBETA
        A = B**(-LT)
        LEPS = A
!
!  Try some tricks to see whether or not this is the correct EPS.
!
        B = TWO / 3
        HALF = ONE / 2
        SIXTH = DLAMC3(B, -HALF)
        THIRD = DLAMC3(SIXTH, SIXTH)
        B = DLAMC3(THIRD, -HALF)
        B = DLAMC3(B, SIXTH)
        B = ABS(B)
        IF (B .LT. LEPS) B = LEPS
!
        LEPS = 1

        DO WHILE ((LEPS.GT.B) .AND. (B.GT.ZERO))
          LEPS = B
          C = DLAMC3(HALF*LEPS, (TWO**5)*(LEPS**2))
          C = DLAMC3(HALF, -C)
          B = DLAMC3(HALF, C)
          C = DLAMC3(HALF, -B)
          B = DLAMC3(HALF, C)
        END DO
!
        IF (A .LT. LEPS) LEPS = A
!
!  Computation of EPS complete.
!
!  Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!  Keep dividing  A by BETA until (gradual) underflow occurs. This
!  is detected when we cannot recover the previous A.
!
        RBASE = ONE / LBETA
        SMALL = ONE
        DO I = 1, 3
          SMALL = DLAMC3(SMALL*RBASE, ZERO)
        END DO
        A = DLAMC3(ONE, SMALL)
        CALL DLAMC4(NGPMIN, ONE, LBETA)
        CALL DLAMC4(NGNMIN, -ONE, LBETA)
        CALL DLAMC4(GPMIN, A, LBETA)
        CALL DLAMC4(GNMIN, -A, LBETA)
        IEEE = .FALSE.
!
        IF ((NGPMIN.EQ.NGNMIN) .AND. (GPMIN.EQ.GNMIN)) THEN
          IF (NGPMIN .EQ. GPMIN) THEN
            LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
          ELSE IF ((GPMIN-NGPMIN) .EQ. 3) THEN
            LEMIN = NGPMIN - 1 + LT
            IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
          ELSE
            LEMIN = MIN(NGPMIN, GPMIN)
!            ( A guess; no known machine )
            IWARN = .TRUE.
          END IF
!
        ELSE IF ((NGPMIN.EQ.GPMIN) .AND. (NGNMIN.EQ.GNMIN)) THEN
          IF (ABS(NGPMIN-NGNMIN) .EQ. 1) THEN
            LEMIN = MAX(NGPMIN, NGNMIN)
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
          ELSE
            LEMIN = MIN(NGPMIN, NGNMIN)
!            ( A guess; no known machine )
            IWARN = .TRUE.
          END IF
!
        ELSE IF ((ABS(NGPMIN-NGNMIN).EQ.1) .AND. (GPMIN.EQ.GNMIN)) THEN
          IF ((GPMIN-MIN(NGPMIN, NGNMIN)) .EQ. 3) THEN
            LEMIN = MAX(NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
          ELSE
            LEMIN = MIN(NGPMIN, NGNMIN)
!            ( A guess; no known machine )
            IWARN = .TRUE.
          END IF
!
        ELSE
          LEMIN = MIN(NGPMIN, NGNMIN, GPMIN, GNMIN)
!         ( A guess; no known machine )
          IWARN = .TRUE.
        END IF
!
!  Comment out this if block if EMIN is okay.
!
        IF (IWARN) THEN
          FIRST = .TRUE.
          PRINT 10, LEMIN
        END IF
!
!  Assume IEEE arithmetic if we found denormalised  numbers above,
!  or if arithmetic seems to round in the  IEEE style,  determined
!  in routine DLAMC1. A true IEEE machine should have both  things
!  true; however, faulty machines may have one or the other.
!
        IEEE = IEEE .OR. LIEEE1
!
!  Compute  RMIN by successive division by  BETA. We could compute
!  RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!  this computation.
!
        LRMIN = 1
        DO I = 1, 1 - LEMIN
          LRMIN = DLAMC3(LRMIN*RBASE, ZERO)
        END DO
!
!  Finally, call DLAMC5 to compute EMAX and RMAX.
!
        CALL DLAMC5 (LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX)
      END IF FIRST_PASS
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
 10   FORMAT ( / / ' WARNING. The value EMIN may be incorrect:-',       &
     &        '  EMIN = ', i0, /,                                       &
     &        ' If, after inspection, the value EMIN looks',            &
     &        ' acceptable please comment out ',                        &
     &        /, ' the IF block as marked within the code of routine',  &
     &        ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
      RETURN
      END SUBROUTINE DLAMC2
!
      FUNCTION DLAMC3 (A, B) RESULT (AplusB)
!
!=======================================================================
!                                                                      !
!  DLAMC3  is intended to force A and B to be stored prior to doing    !
!  the addition of A and B , for use in situations where optimizers    !
!  might hold one of these in a register.                              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  A, B    (input) DOUBLE PRECISION                                    !
!          The values A and B.                                         !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in) ::  A, B
!
!  Local variable declarations.
!
      real (dp) :: AplusB
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      AplusB = A + B
!
      RETURN
      END FUNCTION DLAMC3
!
      SUBROUTINE DLAMC4 (EMIN, START, BASE)
!
!=======================================================================
!                                                                      !
!  DLAMC4 is a service routine for DLAMC2.                             !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  EMIN    (output) EMIN                                               !
!          The minimum exponent before (gradual) underflow, computed   !
!          by setting A = START and dividing by BASE until the         !
!          previous A can not be recovered.                            !
!                                                                      !
!  START   (input) DOUBLE PRECISION                                    !
!          The starting point for determining EMIN.                    !
!                                                                      !
!  BASE    (input) INTEGER                                             !
!          The base of the machine.                                    !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in ) :: BASE
      integer, intent(out) :: EMIN
!
      real (dp), intent(in) :: START
!
!  Local variable declarations.
!
      integer   :: I
!
      real (dp) :: A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3 (A*RBASE, ZERO)
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!
      DO WHILE ((C1.EQ.A) .AND. (C2.EQ.A) .AND.                         &
     &          (D1.EQ.A) .AND. (D2.EQ.A))
        EMIN = EMIN - 1
        A = B1
        B1 = DLAMC3 (A / BASE, ZERO)
        C1 = DLAMC3 (B1*BASE, ZERO)
        D1 = ZERO
        DO I = 1, BASE
          D1 = D1 + B1
        END DO
        B2 = DLAMC3 (A*RBASE, ZERO)
        C2 = DLAMC3 (B2 / RBASE, ZERO)
        D2 = ZERO
        DO I = 1, BASE
          D2 = D2 + B2
        END DO
      END DO
!
      RETURN
      END SUBROUTINE DLAMC4
!
      SUBROUTINE DLAMC5 (BETA, P, EMIN, IEEE, EMAX, RMAX)
!
!=======================================================================
!                                                                      !
!  DLAMC5 attempts to compute RMAX, the largest machine floating-point !
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum     !
!  approximately to a power of 2.  It will fail on machines where this !
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,!
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is !
!  too large (i.e. too close to zero), probably with overflow.         !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  BETA    (input) INTEGER                                             !
!          The base of floating-point arithmetic.                      !
!                                                                      !
!  P       (input) INTEGER                                             !
!          The number of base BETA digits in the mantissa of a         !
!          floating-point value.                                       !
!                                                                      !
!  EMIN    (input) INTEGER                                             !
!          The minimum exponent before (gradual) underflow.            !
!                                                                      !
!  IEEE    (input) LOGICAL                                             !
!          A logical flag specifying whether or not the arithmetic     !
!          system is thought to comply with the IEEE standard.         !
!                                                                      !
!  EMAX    (output) INTEGER                                            !
!          The largest exponent before overflow                        !
!                                                                      !
!  RMAX    (output) DOUBLE PRECISION                                   !
!          The largest machine floating-point number.                  !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in)  :: IEEE
!
      integer, intent(in)  :: BETA, EMIN, P
      integer, intent(out) :: EMAX
!
      real (dp), intent(out) :: RMAX
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      integer :: EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
!
      real (dp) :: OLDY, RECBAS, Y, Z
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  First compute LEXP and UEXP, two powers of 2 that bound
!  ABS(EMIN). We then assume that EMAX + ABS(EMIN) will sum
!  approximately to the bound that is closest to ABS(EMIN).
!  (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
#ifdef REMOVE_LAPACK_GOTOS
      TRY_LOOP : DO                              ! forever DO-loop
        TRY = LEXP*2
        IF (TRY .LE. (-EMIN)) THEN
          LEXP = TRY
          EXBITS = EXBITS + 1
          CYCLE TRY_LOOP                         ! iterate DO-loop again
        ELSE
          EXIT TRY_LOOP                          ! terminate DO-loop
        END IF
      END DO TRY_LOOP
#else
   10 CONTINUE
      TRY = LEXP*2
      IF (TRY .LE. (-EMIN)) THEN
        LEXP = TRY
        EXBITS = EXBITS + 1
        GO TO 10
      END IF
#endif
      IF (LEXP .EQ. -EMIN) THEN
        UEXP = LEXP
      ELSE
        UEXP = TRY
        EXBITS = EXBITS + 1
      END IF
!
!  Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!  than or equal to EMIN. EXBITS is the number of bits needed to
!  store the exponent.
!
      IF ((UEXP+EMIN) .GT. (-LEXP-EMIN)) THEN
        EXPSUM = 2*LEXP
      ELSE
        EXPSUM = 2*UEXP
      END IF
!
!  EXPSUM is the exponent range, approximately equal to
!  EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!  NBITS is the total number of bits needed to store a
!  floating-point number.
!
      IF ((MOD(NBITS, 2).EQ.1) .AND. (BETA.EQ.2)) THEN
!
!  Either there are an odd number of bits used to store a
!  floating-point number, which is unlikely, or some bits are
!  not used in the representation of numbers, which is possible,
!  (e.g. Cray machines) or the mantissa has an implicit bit,
!  (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!  most likely. We have to assume the last alternative.
!  If this is true, then we need to reduce EMAX by one because
!  there must be some way of representing zero in an implicit-bit
!  system. On machines like Cray, we are reducing EMAX by one
!  unnecessarily.
!
        EMAX = EMAX - 1
      END IF
!
      IF (IEEE) THEN
!
!  Assume we are on an IEEE machine which reserves one exponent
!  for infinity and NaN.
!
        EMAX = EMAX - 1
      END IF
!
!  Now create RMAX, the largest machine number, which should
!  be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!  First compute 1.0 - BETA**(-P), being careful that the
!  result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO I = 1, P
        Z = Z*RECBAS
        IF (Y .LT. ONE) OLDY = Y
        Y = DLAMC3(Y, Z)
      END DO
      IF (Y .GE. ONE) Y = OLDY
!
!  Now multiply by BETA**EMAX to get RMAX.
!
      DO I = 1, EMAX
        Y = DLAMC3(Y*BETA, ZERO)
      END DO
!
      RMAX = Y
!
      RETURN
      END SUBROUTINE DLAMC5
!
      FUNCTION DLAMCH (CMACH) RESULT (RMACH)
!
!=======================================================================
!                                                                      !
!  DLAMCH determines double precision machine parameters.              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  CMACH   (input) CHARACTER*1                                         !
!          Specifies the value to be returned by DLAMCH:               !
!          = 'E' or 'e',   DLAMCH := eps                               !
!          = 'S' or 's ,   DLAMCH := sfmin                             !
!          = 'B' or 'b',   DLAMCH := base                              !
!          = 'P' or 'p',   DLAMCH := eps*base                          !
!          = 'N' or 'n',   DLAMCH := t                                 !
!          = 'R' or 'r',   DLAMCH := rnd                               !
!          = 'M' or 'm',   DLAMCH := emin                              !
!          = 'U' or 'u',   DLAMCH := rmin                              !
!          = 'L' or 'l',   DLAMCH := emax                              !
!          = 'O' or 'o',   DLAMCH := rmax                              !
!                                                                      !
!          where                                                       !
!                                                                      !
!          eps   = relative machine precision                          !
!          sfmin = safe minimum, such that 1/sfmin does not overflow   !
!          base  = base of the machine                                 !
!          prec  = eps*base                                            !
!          t     = number of (base) digits in the mantissa             !
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise !
!          emin  = minimum exponent before (gradual) underflow         !
!          rmin  = underflow threshold - base**(emin-1)                !
!          emax  = largest exponent before overflow                    !
!          rmax  = overflow threshold  - (base**emax)*(1-eps)          !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: CMACH
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      logical, save   :: FIRST = .TRUE.
      logical         :: LRND
!
      integer         :: BETA, IMAX, IMIN, IT
!
      real (dp), save :: EPS, SFMIN, BASE, T, RND, EMIN, RMIN,           &
     &                   EMAX, RMAX, PREC
      real (dp)       :: RMACH, SMALL
!
!-----------------------------------------------------------------------
!  Executable Statements ..
!-----------------------------------------------------------------------
!
      IF (FIRST) THEN
        FIRST = .FALSE.
        CALL DLAMC2 (BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX)
        BASE = BETA
        T = IT
        IF (LRND) THEN
          RND = ONE
          EPS = (BASE**(1-IT)) / 2
        ELSE
          RND = ZERO
          EPS = BASE**(1-IT)
        END IF
        PREC = EPS*BASE
        EMIN = IMIN
        EMAX = IMAX
        SFMIN = RMIN
        SMALL = ONE / RMAX
        IF (SMALL.GE.SFMIN) THEN
!
!  Use SMALL plus a bit, to avoid the possibility of rounding
!  causing overflow when computing  1/sfmin.
!
          SFMIN = SMALL*(ONE+EPS)
        END IF
      END IF
!
      IF (LSAME(CMACH, 'E')) THEN
         RMACH = EPS
      ELSE IF (LSAME(CMACH, 'S')) THEN
         RMACH = SFMIN
      ELSE IF (LSAME(CMACH, 'B')) THEN
         RMACH = BASE
      ELSE IF (LSAME(CMACH, 'P')) THEN
         RMACH = PREC
      ELSE IF (LSAME(CMACH, 'N')) THEN
         RMACH = T
      ELSE IF (LSAME(CMACH, 'R')) THEN
         RMACH = RND
      ELSE IF (LSAME(CMACH, 'M')) THEN
         RMACH = EMIN
      ELSE IF (LSAME(CMACH, 'U')) THEN
         RMACH = RMIN
      ELSE IF (LSAME(CMACH, 'L')) THEN
         RMACH = EMAX
      ELSE IF (LSAME(CMACH, 'O')) THEN
         RMACH = RMAX
      END IF
!
      RETURN
      END FUNCTION DLAMCH
!
      FUNCTION DLANST (NORM, N, D, E)
!
!=======================================================================
!                                                                      !
!  DLANST returns the value of the one norm, or the Frobenius norm,    !
!  or the infinity norm, or the element of largest absolute value of   !
!  a real symmetric tridiagonal matrix A.                              !
!                                                                      !
!  Description:                                                        !
!                                                                      !
!  DLANST returns the value                                            !
!                                                                      !
!      ANORM = ( MAX(ABS(A(i,j))), NORM = 'M' or 'm'                   !
!              (                                                       !
!              ( norm1(A),         NORM = '1', 'O' or 'o'              !
!              (                                                       !
!              ( normI(A),         NORM = 'I' or 'i'                   !
!              (                                                       !
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'         !
!  where                                                               !
!                                                                      !
!  norm1 denotes the one norm of a matrix (maximum column sum),        !
!  normI denotes the infinity norm of a matrix (maximum row sum), and  !
!  normF denotes the Frobenius norm of a matrix (square root of sum of !
!        squares).                                                     !
!  Note that MAX(ABS(A(i,j))) is not a matrix norm.                    !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  NORM    (input) CHARACTER*1                                         !
!          Specifies the value to be returned in DLANST as described   !
!          above.                                                      !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The order of the matrix A.  N >= 0.  When N = 0, DLANST is  !
!          set to zero.                                                !
!                                                                      !
!  D       (input) DOUBLE PRECISION array, dimension (N)               !
!          The diagonal elements of A.                                 !
!                                                                      !
!  E       (input) DOUBLE PRECISION array, dimension (N-1)             !
!          The (n-1) sub-diagonal or super-diagonal elements of A.     !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     February 29, 1992                                                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: N
!
      real (dp), intent(in) :: D(*), E(*)
!
      character (len=1), intent(in) :: NORM
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      integer   :: I
      real (dp) :: ANORM, SCALE, SUM
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (N.LE.0) THEN
        ANORM = ZERO
      ELSE IF (LSAME(NORM, 'M')) THEN
!
!  Find MAX(ABS(A(i,j))).
!
        ANORM = ABS(D(N))
        DO I = 1, N - 1
          ANORM = MAX(ANORM, ABS(D(I)))
          ANORM = MAX(ANORM, ABS(E(I)))
        END DO
      ELSE IF (LSAME(NORM, 'O') .OR. NORM.EQ.'1' .OR.                   &
     &         LSAME(NORM, 'I')) THEN
!
!  Find norm1(A).
!
        IF (N.EQ.1) THEN
          ANORM = ABS(D(1))
        ELSE
          ANORM = MAX(ABS(D(1)) + ABS(E(1)),                            &
     &                ABS(E(N-1)) + ABS(D(N)))
          DO I = 2, N - 1
            ANORM = MAX(ANORM, ABS(D(I)) + ABS(E(I)) + ABS(E(I-1)))
          END DO
        END IF
      ELSE IF ((LSAME(NORM, 'F')) .OR. (LSAME(NORM, 'E'))) THEN
!
!  Find normF(A).
!
        SCALE = ZERO
        SUM = ONE
        IF (N .GT. 1) THEN
          CALL DLASSQ (N-1, E, 1, SCALE, SUM)
          SUM = 2*SUM
       END IF
       CALL DLASSQ (N, D, 1, SCALE, SUM)
       ANORM = SCALE*SQRT(SUM)
     END IF
!
     DLANST = ANORM
!
     RETURN
     END FUNCTION DLANST
!
     FUNCTION DLAPY2 (X, Y) RESULT (LAPY2)
!
!=======================================================================
!                                                                      !
!  DLAPY2 returns SQRT(x**2+y**2), taking care not to cause            !
!  unnecessary overflow.                                               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  X       (input) DOUBLE PRECISION                                    !
!  Y       (input) DOUBLE PRECISION                                    !
!          X and Y specify the values x and y.                         !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in) :: X, Y
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      real (dp) :: W, XABS, YABS, Z
      real (dp) :: LAPY2
!
!-----------------------------------------------------------------------
!  Executable Statements ..
!-----------------------------------------------------------------------
!
      XABS = ABS(X)
      YABS = ABS(Y)
      W = MAX(XABS, YABS)
      Z = MIN(XABS, YABS)
      IF (Z .EQ. ZERO) THEN
        LAPY2 = W
      ELSE
        LAPY2 = W*SQRT(ONE+(Z / W)**2)
      END IF
!
      RETURN
      END FUNCTION DLAPY2
!
      SUBROUTINE DLARTG (F, G, CS, SN, R)
!
!=======================================================================
!                                                                      !
!  DLARTG generate a plane rotation so that                            !
!                                                                      !
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.      !
!     [ -SN  CS  ]     [ G ]     [ 0 ]                                 !
!                                                                      !
!  This is a slower, more accurate version of the BLAS1 routine DROTG, !
!  with the following other differences:                               !
!                                                                      !
!     F and G are unchanged on return.                                 !
!     If G=0, then CS=1 and SN=0.                                      !
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any      !
!        floating point operations (saves work in DBDSQR when          !
!        there are zeros on the diagonal).                             !
!                                                                      !
!  If F exceeds G in magnitude, CS will be positive.                   !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  F       (input) DOUBLE PRECISION                                    !
!          The first component of vector to be rotated.                !
!                                                                      !
!  G       (input) DOUBLE PRECISION                                    !
!          The second component of vector to be rotated.               !
!                                                                      !
!  CS      (output) DOUBLE PRECISION                                   !
!          The cosine of the rotation.                                 !
!                                                                      !
!  SN      (output) DOUBLE PRECISION                                   !
!          The sine of the rotation.                                   !
!                                                                      !
!  R       (output) DOUBLE PRECISION                                   !
!          The nonzero component of the rotated vector.                !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in ) :: F, G
      real (dp), intent(out) :: CS, R, SN
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
!
      logical, save   :: FIRST = .TRUE.
!
      integer         :: COUNT, I
!
      real (dp), save :: SAFMX2, SAFMIN, SAFMN2
      real (dp)       :: EPS, F1, G1, SCALE
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (FIRST) THEN
        FIRST = .FALSE.
        SAFMIN = DLAMCH('S')
        EPS = DLAMCH('E')
        SAFMN2 = DLAMCH('B')**INT(LOG(SAFMIN / EPS) /                   &
     &           LOG(DLAMCH('B')) / TWO)
        SAFMX2 = ONE / SAFMN2
      END IF
      IF (G .EQ. ZERO) THEN
        CS = ONE
        SN = ZERO
        R = F
      ELSE IF (F .EQ. ZERO) THEN
        CS = ZERO
        SN = ONE
        R = G
      ELSE
        F1 = F
        G1 = G
        SCALE = MAX(ABS(F1), ABS(G1))
        IF (SCALE .GE. SAFMX2) THEN
          COUNT = 0
#ifdef REMOVE_LAPACK_GOTOS
          LOOP_1 : DO                            ! forever DO-loop
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX(ABS(F1), ABS(G1))
            IF (SCALE .GE. SAFMX2) CYCLE LOOP_1  ! iterate DO-loop again
            R = SQRT(F1**2 + G1**2)
            CS = F1 / R
            SN = G1 / R
            DO I = 1, COUNT
              R = R*SAFMX2
            END DO
            EXIT LOOP_1                          ! terminate DO-loop
          END DO LOOP_1
#else
   10     CONTINUE
          COUNT = COUNT + 1
          F1 = F1*SAFMN2
          G1 = G1*SAFMN2
          SCALE = MAX(ABS(F1), ABS(G1))
          IF (SCALE .GE. SAFMX2) GO TO 10
          R = SQRT(F1**2 + G1**2)
          CS = F1 / R
          SN = G1 / R
          DO I = 1, COUNT
            R = R*SAFMX2
          END DO
#endif
        ELSE IF (SCALE .LE. SAFMN2) THEN
          COUNT = 0
#ifdef REMOVE_LAPACK_GOTOS
          LOOP_2 : DO                            ! forever DO-loop
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX(ABS(F1), ABS(G1))
            IF (SCALE .LE. SAFMN2) CYCLE LOOP_2  ! iterate DO-loop again
            R = SQRT(F1**2 + G1**2)
            CS = F1 / R
            SN = G1 / R
            DO I = 1, COUNT
              R = R*SAFMN2
            END DO
            EXIT LOOP_2                          ! terminate DO-loop
          END DO LOOP_2
#else
   30     CONTINUE
          COUNT = COUNT + 1
          F1 = F1*SAFMX2
          G1 = G1*SAFMX2
          SCALE = MAX(ABS(F1), ABS(G1))
          IF (SCALE .LE. SAFMN2) GO TO 30
          R = SQRT(F1**2 + G1**2)
          CS = F1 / R
          SN = G1 / R
          DO I = 1, COUNT
            R = R*SAFMN2
          END DO
#endif
        ELSE
          R = SQRT(F1**2 + G1**2)
          CS = F1 / R
          SN = G1 / R
        END IF
        IF ((ABS(F).GT.ABS(G)) .AND. (CS.LT.ZERO)) THEN
          CS = -CS
          SN = -SN
          R = -R
        END IF
      END IF
!
      RETURN
      END SUBROUTINE DLARTG
!
      SUBROUTINE DLASCL (TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO)
!
!=======================================================================
!                                                                      !
!  DLASCL multiplies the M by N real matrix A by the real scalar       !
!  CTO/CFROM. This is done without over/underflow as long as the       !
!  final result CTO*A(I,J)/CFROM does not over/underflow. TYPE         !
!  specifies that A may be full, upper triangular, lower triangular,   !
!  upper Hessenberg, or banded.                                        !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  TYPE    (input) CHARACTER*1                                         !
!          TYPE indices the storage type of the input matrix.          !
!          = 'G':  A is a full matrix.                                 !
!          = 'L':  A is a lower triangular matrix.                     !
!          = 'U':  A is an upper triangular matrix.                    !
!          = 'H':  A is an upper Hessenberg matrix.                    !
!          = 'B':  A is a symmetric band matrix with lower bandwidth   !
!                  KL and upper bandwidth KU and with the only the     !
!                  lower half stored.                                  !
!          = 'Q':  A is a symmetric band matrix with lower bandwidth   !
!                  KL and upper bandwidth KU and with the only the     !
!                  upper half stored.                                  !
!          = 'Z':  A is a band matrix with lower bandwidth KL and      !
!                  upper bandwidth KU.                                 !
!                                                                      !
!  KL      (input) INTEGER                                             !
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',   !
!          'Q' or 'Z'.                                                 !
!                                                                      !
!  KU      (input) INTEGER                                             !
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',   !
!          'Q' or 'Z'.                                                 !
!                                                                      !
!  CFROM   (input) DOUBLE PRECISION                                    !
!  CTO     (input) DOUBLE PRECISION                                    !
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed !
!          without over/underflow if the final result CTO*A(I,J)/CFROM !
!          can be represented without over/underflow.  CFROM must be   !
!          nonzero.                                                    !
!                                                                      !
!  M       (input) INTEGER                                             !
!          The number of rows of the matrix A.  M >= 0.                !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The number of columns of the matrix A.  N >= 0.             !
!                                                                      !
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)    !
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the !
!          storage type.                                               !
!                                                                      !
!  LDA     (input) INTEGER                                             !
!          The leading dimension of the array A.  LDA >= max(1,M).     !
!                                                                      !
!  INFO    (output) INTEGER                                            !
!          0  - successful exit                                        !
!          <0 - if INFO = -i, the i-th argument had an illegal value.  !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     February 29, 1992                                                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
!
      integer, intent (in ) :: KL, KU, LDA, M, N
      integer, intent (out) :: INFO
!
      real (dp), intent(in   ) :: CFROM, CTO
      real (dp), intent(inout) :: A(LDA, *)
!
      character (len=1), intent(in) :: TYPE
!
!  Local variable declarations.
!
      real (r8), parameter :: ZERO = 0.0_dp
      real (r8), parameter :: ONE = 1.0_dp
!
      logical   :: DONE
!
      integer   :: I, ITYPE, J, K1, K2, K3, K4
!
      real (dp) :: BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input arguments.
!
      INFO = 0
!
      IF (LSAME(TYPE, 'G')) THEN
        ITYPE = 0
      ELSE IF (LSAME(TYPE, 'L')) THEN
        ITYPE = 1
      ELSE IF (LSAME(TYPE, 'U')) THEN
        ITYPE = 2
      ELSE IF (LSAME(TYPE, 'H')) THEN
        ITYPE = 3
      ELSE IF (LSAME(TYPE, 'B')) THEN
        ITYPE = 4
      ELSE IF (LSAME(TYPE, 'Q')) THEN
        ITYPE = 5
      ELSE IF (LSAME(TYPE, 'Z')) THEN
        ITYPE = 6
      ELSE
        ITYPE = -1
      END IF
!
      IF (ITYPE .EQ. -1) THEN
        INFO = -1
      ELSE IF (CFROM .EQ. ZERO) THEN
        INFO = -4
      ELSE IF (M .LT. 0) THEN
        INFO = -6
      ELSE IF (N.LT.0 .OR. (ITYPE.EQ.4 .AND. N.NE.M) .OR.               &
     &                     (ITYPE.EQ.5 .AND. N.NE.M)) THEN
        INFO = -7
      ELSE IF (ITYPE.LE.3 .AND. LDA.LT.MAX(1, M)) THEN
        INFO = -9
      ELSE IF (ITYPE .GE. 4) THEN
        IF (KL.LT.0 .OR. KL.GT.MAX(M-1, 0)) THEN
          INFO = -2
        ELSE IF (KU.LT.0 .OR. KU.GT.MAX(N-1, 0) .OR.                    &
     &           ((ITYPE.EQ.4 .OR. ITYPE.EQ.5) .AND. KL.NE.KU)) THEN
          INFO = -3
        ELSE IF ((ITYPE.EQ.4 .AND. LDA.LT.KL+1) .OR.                    &
     &           (ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.                   &
     &           (ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1)) THEN
          INFO = -9
        END IF
      END IF
!
      IF (INFO .NE. 0) THEN
         CALL XERBLA ('DLASCL', -INFO)
         RETURN
      END IF
!
!  Quick return if possible.
!
      IF ((N.EQ.0) .OR. (M.EQ.0))  RETURN
!
!  Get machine parameters
!
      SMLNUM = DLAMCH ('S')
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
      DONE = .FALSE.
!
      ITERATE : DO WHILE (.NOT. DONE)
        CFROM1 = CFROMC*SMLNUM
        CTO1 = CTOC / BIGNUM
        IF (ABS(CFROM1).GT.ABS(CTOC) .AND. CTOC.NE.ZERO) THEN
          MUL = SMLNUM
          DONE = .FALSE.
          CFROMC = CFROM1
        ELSE IF (ABS(CTO1) .GT. ABS(CFROMC)) THEN
          MUL = BIGNUM
          DONE = .FALSE.
          CTOC = CTO1
        ELSE
          MUL = CTOC / CFROMC
          DONE = .TRUE.
        END IF
!
        IF (ITYPE .EQ. 0) THEN
!
!  Full matrix.
!
          DO J = 1, N
            DO I = 1, M
               A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE .EQ. 1) THEN
!
!  Lower triangular matrix.
!
          DO J = 1, N
            DO I = J, M
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE.EQ.2) THEN
!
!  Upper triangular matrix.
!
          DO J = 1, N
            DO I = 1, MIN(J, M)
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE .EQ. 3) THEN
!
!  Upper Hessenberg matrix.
!
           DO J = 1, N
             DO I = 1, MIN(J+1, M)
               A(I,J) = A(I,J)*MUL
             END DO
           END DO
!
        ELSE IF (ITYPE .EQ. 4) THEN
!
!  Lower half of a symmetric band matrix.
!
          K3 = KL + 1
          K4 = N + 1
          DO J = 1, N
            DO I = 1, MIN(K3, K4-J)
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE .EQ. 5) THEN
!
!  Upper half of a symmetric band matrix.
!
          K1 = KU + 2
          K3 = KU + 1
          DO J = 1, N
            DO I = MAX(K1-J, 1), K3
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE .EQ. 6) THEN
!
!  Band matrix.
!
          K1 = KL + KU + 2
          K2 = KL + 1
          K3 = 2*KL + KU + 1
          K4 = KL + KU + 1 + M
          DO J = 1, N
            DO I = MAX(K1-J, K2), MIN(K3, K4-J)
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        END IF
      END DO ITERATE
!
      RETURN
      END SUBROUTINE DLASCL
!
      SUBROUTINE DLASET (UPLO, M, N, ALPHA, BETA, A, LDA)
!
!=======================================================================
!                                                                      !
!  DLASET initializes an m-by-n matrix A to BETA on the diagonal and   !
!  ALPHA on the offdiagonals.                                          !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  UPLO    (input) CHARACTER*1                                         !
!          Specifies the part of the matrix A to be set.               !
!          = 'U':      Upper triangular part is set; the strictly      !
!                      lower triangular part of A is not changed.      !
!          = 'L':      Lower triangular part is set; the strictly      !
!                      upper triangular part of A is not changed.      !
!          Otherwise:  All of the matrix A is set.                     !
!                                                                      !
!  M       (input) INTEGER                                             !
!          The number of rows of the matrix A.  M >= 0.                !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The number of columns of the matrix A.  N >= 0.             !
!                                                                      !
!  ALPHA   (input) DOUBLE PRECISION                                    !
!          The constant to which the offdiagonal elements are to be    !
!          set.                                                        !
!                                                                      !
!  BETA    (input) DOUBLE PRECISION                                    !
!          The constant to which the diagonal elements are to be set.  !
!                                                                      !
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)    !
!          On exit, the leading m-by-n submatrix of A is set as        !
!          follows:                                                    !
!                                                                      !
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,          !
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,          !
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,    !
!                                                                      !
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).           !
!                                                                      !
!  LDA     (input) INTEGER                                             !
!          The leading dimension of the array A.  LDA >= max(1,M).     !
!                                                                      !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)      :: LDA, M, N
!
      real (dp), intent(in)    :: ALPHA, BETA
      real (dp), intent(inout) :: A(LDA, *)
!
      character (len=1), intent(in) :: UPLO
!
!  Local variable declarations.
!
      integer :: I, J
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (LSAME(UPLO, 'U')) THEN
!
!  Set the strictly upper triangular or trapezoidal part of the
!  array to ALPHA.
!
        DO J = 2, N
          DO I = 1, MIN(J-1, M)
            A(I,J) = ALPHA
          END DO
        END DO
!
      ELSE IF (LSAME(UPLO, 'L')) THEN
!
!  Set the strictly lower triangular or trapezoidal part of the
!  array to ALPHA.
!
        DO J = 1, MIN(M, N)
          DO I = J + 1, M
            A(I,J) = ALPHA
          END DO
        END DO
!
      ELSE
!
!  Set the leading m-by-n submatrix to ALPHA.
!
        DO J = 1, N
          DO I = 1, M
            A(I,J) = ALPHA
          END DO
        END DO
      END IF
!
!  Set the first MIN(M,N) diagonal elements to BETA.
!
      DO I = 1, MIN(M, N)
        A(I,I) = BETA
      END DO
!
      RETURN
      END SUBROUTINE DLASET
!
      SUBROUTINE DLASR (SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA)
!
!=======================================================================
!                                                                      !
!  DLASR performs the transformation:                                  !
!                                                                      !
!     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )          !
!                                                                      !
!     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )          !
!                                                                      !
!  where A is an m by n real matrix and P is an orthogonal matrix,     !
!  consisting of a sequence of plane rotations determined by the       !
!  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L'      !
!  or 'l' and z = n when SIDE = 'R' or 'r' ):                          !
!                                                                      !
!  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then                !
!                                                                      !
!     P = P( z - 1 )*...*P( 2 )*P( 1 ),                                !
!                                                                      !
!  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then            !
!                                                                      !
!     P = P( 1 )*P( 2 )*...*P( z - 1 ),                                !
!                                                                      !
!  where  P( k ) is a plane rotation matrix for the following planes:  !
!                                                                      !
!     when  PIVOT = 'V' or 'v'  ( Variable pivot ),                    !
!        the plane ( k, k + 1 )                                        !
!                                                                      !
!     when  PIVOT = 'T' or 't'  ( Top pivot ),                         !
!        the plane ( 1, k + 1 )                                        !
!                                                                      !
!     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),                      !
!        the plane ( k, z )                                            !
!                                                                      !
!  c( k ) and s( k )  must contain the  cosine and sine that define    !
!  the matrix  P( k ).  The two by two plane rotation part of the      !
!  matrix P( k ), R( k ), is assumed to be of the form                 !
!                                                                      !
!     R( k ) = (  c( k )  s( k ) ).                                    !
!              ( -s( k )  c( k ) )                                     !
!                                                                      !
!  This version vectorises across rows of the array A when             !
!  SIDE = 'L'.                                                         !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  SIDE    (input) CHARACTER*1                                         !
!          Specifies whether the plane rotation matrix P is applied    !
!          to A on the left or the right.                              !
!          = 'L':  Left, compute A := P*A                              !
!          = 'R':  Right, compute A:= A*P'                             !
!                                                                      !
!  DIRECT  (input) CHARACTER*1                                         !
!          Specifies whether P is a forward or backward sequence of    !
!          plane rotations.                                            !
!          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )           !
!          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )          !
!                                                                      !
!  PIVOT   (input) CHARACTER*1                                         !
!          Specifies the plane for which P(k) is a plane rotation      !
!          matrix.                                                     !
!          = 'V':  Variable pivot, the plane (k,k+1)                   !
!          = 'T':  Top pivot, the plane (1,k+1)                        !
!          = 'B':  Bottom pivot, the plane (k,z)                       !
!                                                                      !
!  M       (input) INTEGER                                             !
!          The number of rows of the matrix A.  If m <= 1, an          !
!          immediate return is effected.                               !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The number of columns of the matrix A.  If n <= 1, an       !
!          immediate return is effected.                               !
!                                                                      !
!  C, S    (input) DOUBLE PRECISION arrays, dimension                  !
!                  (M-1) if SIDE = 'L'                                 !
!                  (N-1) if SIDE = 'R'                                 !
!          c(k) and s(k) contain the cosine and sine that define the   !
!          matrix P(k).  The two by two plane rotation part of the     !
!          matrix P(k), R(k), is assumed to be of the form             !
!          R( k ) = (  c( k )  s( k ) ).                               !
!                   ( -s( k )  c( k ) )                                !
!                                                                      !
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)    !
!          The m by n matrix A.  On exit, A is overwritten by P*A if   !
!          SIDE = 'R' or by A*P' if SIDE = 'L'.                        !
!                                                                      !
!  LDA     (input) INTEGER                                             !
!          The leading dimension of the array A.  LDA >= max(1,M).     !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)      :: LDA, M, N
!
      real (dp), intent(in   ) :: C(*), S(*)
      real (dp), intent(inout) :: A(LDA, *)
!
      character (len=1), intent(in) :: DIRECT, PIVOT, SIDE
!
!  Local variable declarations.
!
      real (r8), parameter :: ZERO = 0.0_dp
      real (r8), parameter :: ONE = 1.0_dp
!
      integer :: I, INFO, J
!
      real (dp) :: CTEMP, STEMP, TEMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input parameters.
!
      INFO = 0
      IF (.NOT. (LSAME(SIDE, 'L') .OR. LSAME(SIDE, 'R'))) THEN
        INFO = 1
      ELSE IF (.NOT.(LSAME(PIVOT, 'V') .OR. LSAME( PIVOT, 'T') .OR.     &
     &         LSAME(PIVOT, 'B'))) THEN
        INFO = 2
      ELSE IF (.NOT.(LSAME(DIRECT, 'F') .OR. LSAME(DIRECT, 'B'))) THEN
        INFO = 3
      ELSE IF (M .LT. 0) THEN
        INFO = 4
      ELSE IF (N .LT. 0) THEN
        INFO = 5
      ELSE IF (LDA .LT. MAX(1, M)) THEN
        INFO = 9
      END IF
      IF (INFO .NE. 0) THEN
        CALL XERBLA ('DLASR ', INFO)
        RETURN
      END IF
!
!  Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0)) RETURN
!
      IF (LSAME( SIDE, 'L')) THEN
!
!  Form P * A.
!
        IF (LSAME(PIVOT, 'V')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 1, M - 1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J+1,I)
                   A( J+1,I) = CTEMP*TEMP - STEMP*A(J,I)
                   A( J,  I) = STEMP*TEMP + CTEMP*A(J,I)
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = M - 1, 1, -1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. ( STEMP.NE.ZERO ) ) THEN
                DO I = 1, N
                  TEMP = A(J+1,I)
                  A(J+1,I) = CTEMP*TEMP - STEMP*A(J,I)
                  A(J,  I) = STEMP*TEMP + CTEMP*A(J,I)
                END DO
              END IF
            END DO
          END IF
        ELSE IF (LSAME(PIVOT, 'T')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 2, M
              CTEMP = C(J-1)
              STEMP = S(J-1)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J,I)
                  A(J,I) = CTEMP*TEMP - STEMP*A(1,I)
                  A(1,I) = STEMP*TEMP + CTEMP*A(1,I)
                END DO
              END IF
            END DO
          ELSE IF (LSAME( DIRECT, 'B')) THEN
            DO J = M, 2, -1
              CTEMP = C(J-1)
              STEMP = S(J-1)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J,I)
                  A(J,I) = CTEMP*TEMP - STEMP*A(1,I)
                  A(1,I) = STEMP*TEMP + CTEMP*A(1,I)
                END DO
              END IF
            END DO
          END IF
        ELSE IF (LSAME(PIVOT, 'B')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 1, M - 1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J,I)
                  A(J,I) = STEMP*A(M,I) + CTEMP*TEMP
                  A(M,I) = CTEMP*A(M,I) - STEMP*TEMP
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = M - 1, 1, -1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J,I)
                  A(J,I) = STEMP*A(M,I) + CTEMP*TEMP
                  A(M,I) = CTEMP*A(M,I) - STEMP*TEMP
                END DO
              END IF
            END DO
          END IF
        END IF
      ELSE IF (LSAME(SIDE, 'R')) THEN
!
!  Form A * P'.
!
        IF (LSAME(PIVOT, 'V')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 1, N - 1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J+1)
                  A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J)
                  A(I,J  ) = STEMP*TEMP + CTEMP*A(I,J)
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = N - 1, 1, -1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J+1)
                  A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J)
                  A(I,J  ) = STEMP*TEMP + CTEMP*A(I,J)
                END DO
              END IF
            END DO
          END IF
        ELSE IF (LSAME(PIVOT, 'T')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 2, N
              CTEMP = C(J-1)
              STEMP = S(J-1)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J)
                  A(I,J) = CTEMP*TEMP - STEMP*A(I,1)
                  A(I,1) = STEMP*TEMP + CTEMP*A(I,1)
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = N, 2, -1
              CTEMP = C(J-1)
              STEMP = S(J-1)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J)
                  A(I,J) = CTEMP*TEMP - STEMP*A(I,1)
                  A(I,1) = STEMP*TEMP + CTEMP*A(I,1)
                END DO
              END IF
            END DO
          END IF
        ELSE IF (LSAME(PIVOT, 'B')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 1, N - 1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE ) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J)
                  A(I,J) = STEMP*A(I,N) + CTEMP*TEMP
                  A(I,N) = CTEMP*A(I,N) - STEMP*TEMP
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = N - 1, 1, -1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J)
                  A(I,J) = STEMP*A(I,N) + CTEMP*TEMP
                  A(I,N) = CTEMP*A(I,N) - STEMP*TEMP
                END DO
              END IF
            END DO
          END IF
        END IF
      END IF
!
      RETURN
      END SUBROUTINE DLASR
!
      SUBROUTINE DLASRT (ID, N, D, INFO)
!
!=======================================================================
!                                                                      !
!  Sort the numbers in D in increasing order (if ID = 'I') or          !
!  in decreasing order (if ID = 'D' ).                                 !
!                                                                      !
!  Use Quick Sort, reverting to Insertion sort on arrays of            !
!  size <= 20. Dimension of STACK limits N to about 2**32.             !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  ID      (input) CHARACTER*1                                         !
!          = 'I': sort D in increasing order;                          !
!          = 'D': sort D in decreasing order.                          !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The length of the array D.                                  !
!                                                                      !
!  D       (input/output) DOUBLE PRECISION array, dimension (N)        !
!          On entry, the array to be sorted.                           !
!          On exit, D has been sorted into increasing order            !
!          (D(1) <= ... <= D(N) ) or into decreasing order             !
!          (D(1) >= ... >= D(N) ), depending on ID.                    !
!                                                                      !
!  INFO    (output) INTEGER                                            !
!          = 0:  successful exit                                       !
!          < 0:  if INFO = -i, the i-th argument had an illegal value  !
!                                                                      !
!  -- LAPACK routine (version 2.0) --                                  !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in ) :: N
      integer, intent(out) :: INFO
!
      real (dp), intent(inout) :: D(*)
!
      character (len=1), intent(in) :: ID
!
!  Local variable declarations.
!
      integer, parameter :: SELECT = 20
!
      integer :: DIR, ENDD, I, J, START, STKPNT
      integer :: STACK( 2, 32 )
!
      real (dp) :: D1, D2, D3, DMNMX, TMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input paramters.
!
      INFO = 0
      DIR = -1
      IF (LSAME(ID, 'D')) THEN
        DIR = 0
      ELSE IF (LSAME(ID, 'I')) THEN
        DIR = 1
      END IF
!
      IF (DIR .EQ. -1) THEN
        INFO = -1
      ELSE IF (N .LT. 0) THEN
        INFO = -2
      END IF
      IF (INFO .NE. 0) THEN
        CALL XERBLA ('DLASRT', -INFO)
        RETURN
      END IF
!
!  Quick return if possible
!
      IF (N .LE. 1) RETURN
!
      STKPNT = 1
      STACK(1,1) = 1
      STACK(2,1) = N
!
      SORT_LOOP : DO WHILE (STKPNT .GT. 0)
        START = STACK(1,STKPNT)
        ENDD  = STACK(2,STKPNT)
        STKPNT = STKPNT - 1
        IF ((ENDD-START.LE.SELECT) .AND. (ENDD-START.GT.0)) THEN
!
!  Do Insertion sort on D(START:ENDD).
!
          IF (DIR .EQ. 0) THEN
!
!    Sort into decreasing order.
!
            OUTER1 : DO I = START + 1, ENDD
              INNER1 : DO J = I, START + 1, -1
                IF (D(J) .GT. D(J-1)) THEN
                  DMNMX = D(J)
                  D(J  ) = D(J-1)
                  D(J-1) = DMNMX
                ELSE
                  EXIT INNER1
                END IF
              END DO INNER1
            END DO OUTER1
!
          ELSE
!
!    Sort into increasing order.
!
            OUTER2 : DO I = START + 1, ENDD
              INNER2 : DO J = I, START + 1, -1
                IF (D(J) .LT. D(J-1)) THEN
                  DMNMX = D(J)
                  D(J  ) = D(J-1)
                  D(J-1) = DMNMX
                ELSE
                  EXIT INNER2
                END IF
              END DO INNER2
            END DO OUTER2
!
          END IF
!
        ELSE IF (ENDD-START .GT. SELECT) THEN
!
!  Partition D( START:ENDD ) and stack parts, largest one first
!
!    Choose partition entry as median of 3
!
          D1 = D(START)
          D2 = D(ENDD)
          I = (START+ENDD) / 2
          D3 = D( I )
          IF (D1 .LT. D2) THEN
            IF (D3 .LT. D1) THEN
              DMNMX = D1
            ELSE IF (D3 .LT. D2) THEN
              DMNMX = D3
            ELSE
              DMNMX = D2
            END IF
          ELSE
            IF (D3 .LT. D2) THEN
              DMNMX = D2
            ELSE IF (D3 .LT. D1) THEN
              DMNMX = D3
            ELSE
              DMNMX = D1
            END IF
          END IF
!
          IF (DIR .EQ. 0) THEN
!
!    Sort into decreasing order.
!
            I = START - 1
            J = ENDD + 1
#ifdef REMOVE_LAPACK_GOTOS
            DecLOOP_3 : DO                           ! forever LOOP 3
              DecLOOP_1 : DO                         ! forever LOOP 1
                J = J - 1
                IF (D(J) .GE. DMNMX) EXIT DecLOOP_1  ! terminate LOOP 1
              END DO DecLOOP_1
!
              DecLOOP_2 : DO                         ! forever LOOP 2
                I = I + 1
                IF (D(I) .LE. DMNMX) EXIT DecLOOP_2  ! terminate LOOP 2
              END DO DecLOOP_2
!
              IF (I .LT. J) THEN
                TMP = D(I)
                D(I) = D(J)
                D(J) = TMP
                CYCLE DecLOOP_3                      ! iterate LOOP 3
              ELSE
                EXIT DecLOOP_3                       ! terminate LOOP 3
              END IF
            END DO DecLOOP_3
#else
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF (D(J) .LT. DMNMX) GO TO 70
   80       CONTINUE
            I = I + 1
            IF (D(I) .GT. DMNMX) GO TO 80
            IF (I .LT. J) THEN
              TMP = D(I)
              D(I) = D(J)
              D(J) = TMP
              GO TO 60
            END IF
#endif
            IF (J-START .GT. ENDD-J-1) THEN
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = START
              STACK(2, STKPNT) = J
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = J + 1
              STACK(2, STKPNT) = ENDD
            ELSE
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = J + 1
              STACK(2, STKPNT) = ENDD
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = START
              STACK(2, STKPNT) = J
            END IF
          ELSE
!
!    Sort into increasing order.
!
            I = START - 1
            J = ENDD + 1
#ifdef REMOVE_LAPACK_GOTOS
            IncLOOP_3 : DO                           ! forever LOOP 3
              IncLOOP_1 : DO                         ! forever LOOP 1
                J = J - 1
                IF (D(J) .LE. DMNMX) EXIT IncLOOP_1  ! terminate LOOP 1
              END DO IncLOOP_1
!
              IncLOOP_2 : DO                         ! forever LOOP 2
                I = I + 1
                IF (D(I) .GE. DMNMX) EXIT IncLOOP_2  ! terminate LOOP 2
              END DO IncLOOP_2
!
              IF (I .LT. J) THEN
                TMP = D(I)
                D(I) = D(J)
                D(J) = TMP
                CYCLE IncLOOP_3                      ! iterate LOOP 3
              ELSE
                EXIT IncLOOP_3                       ! terminate LOOP 3
              END IF
            END DO IncLOOP_3
#else
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF (D(J) .GT. DMNMX) GO TO 100
  110       CONTINUE
            I = I + 1
            IF (D(I) .LT. DMNMX) GO TO 110
            IF (I .LT. J) THEN
              TMP = D(I)
              D(I) = D(J)
              D(J) = TMP
              GO TO 90
            END IF
#endif
            IF (J-START .GT. ENDD-J-1) THEN
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = START
              STACK(2, STKPNT) = J
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = J + 1
              STACK(2, STKPNT) = ENDD
            ELSE
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = J + 1
              STACK(2, STKPNT) = ENDD
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = START
              STACK(2, STKPNT) = J
            END IF
          END IF
        END IF
      END DO SORT_LOOP
!
      RETURN
      END SUBROUTINE DLASRT
!
      SUBROUTINE DLASSQ (N, X, INCX, SCALE, SUMSQ)
!
!=======================================================================
!                                                                      !
!  DLASSQ  returns the values  scl  and  smsq  such that               !
!                                                                      !
!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq !
!                                                                      !
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is    !
!  assumed to be non-negative and  scl  returns the value              !
!                                                                      !
!     scl = max( scale, abs( x( i ) ) ).                               !
!                                                                      !
!  scale and sumsq must be supplied in SCALE and SUMSQ and             !
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.       !
!                                                                      !
!  The routine makes only one pass through the vector x.               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The number of elements to be used from the vector X.        !
!                                                                      !
!  X       (input) DOUBLE PRECISION                                    !
!          The vector for which a scaled sum of squares is computed.   !
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.          !
!                                                                      !
!  INCX    (input) INTEGER                                             !
!          The increment between successive values of the vector X.    !
!          INCX > 0.                                                   !
!                                                                      !
!  SCALE   (input/output) DOUBLE PRECISION                             !
!          On entry, the value  scale  in the equation above.          !
!          On exit, SCALE is overwritten with  scl, the scaling factor !
!          for the sum of squares.                                     !
!                                                                      !
!  SUMSQ   (input/output) DOUBLE PRECISION                             !
!          On entry, the value  sumsq  in the equation above.          !
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of !
!          squares from which  scl  has been factored out.             !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: INCX, N
!
      real (dp), intent(in   ) :: X(*)
      real (dp), intent(inout) :: SCALE, SUMSQ
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
!
      integer   :: IX
!
      real (dp) :: ABSXI
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (N .GT. 0) THEN
        DO IX = 1, 1 + ( N-1 )*INCX, INCX
          IF (X(IX) .NE. ZERO) THEN
            ABSXI = ABS(X(IX))
            IF (SCALE .LT. ABSXI) THEN
              SUMSQ = 1.0_dp + SUMSQ*(SCALE / ABSXI)**2
              SCALE = ABSXI
            ELSE
              SUMSQ = SUMSQ + (ABSXI / SCALE)**2
            END IF
          END IF
        END DO
      END IF
!
      RETURN
      END SUBROUTINE DLASSQ
!
      SUBROUTINE  DSWAP (N, DX, INCX, DY, INCY)
!
!=======================================================================
!                                                                      !
!  DSWAP interchanges two vectors. It uses unrolled loops for          !
!  increments equal one.                                               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N       (input) INTEGER                                             !
!          number of elements in input vector(s)                       !
!                                                                      !
!  DX      (input/output) DOUBLE PRECISION                             !
!          array, dimension ( 1 + ( N - 1 )*abs( INCX ) )              !
!                                                                      !
!  INCX    (input) INTEGER                                             !
!          storage spacing between elements of DX                      !
!                                                                      !
!  DY      (input/output) DOUBLE PRECISION                             !
!          array, dimension ( 1 + ( N - 1 )*abs( INCY ) )              !
!                                                                      !
!  INCY    (input) INTEGER                                             !
!          storage spacing between elements of DY                      !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!  Based on Jack Dongarra, LINPACK, 3/11/78.                           !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: INCX, INCY, N
!
      real (dp), intent(inout) :: DX(*), DY(*)
!
!  Local variable declarations.
!
      integer I, IX, IY, M, MP1
!
      real (dp) :: DTEMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (N .le. 0) RETURN
!
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!  Code for both increments equal to 1.
!
        M = MOD(N,3)
        IF (M .NE. 0 ) THEN
          DO I = 1, M
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
          END DO
          IF (N .LT. 3) RETURN
        END IF
        MP1 = M + 1
        DO I = MP1, N, 3
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
          DTEMP = DX(I+1)
          DX(I+1) = DY(I+1)
          DY(I+1) = DTEMP
          DTEMP = DX(I+2)
          DX(I+2) = DY(I+2)
          DY(I+2) = DTEMP
        END DO
      ELSE
!
!  Code for unequal increments or equal increments not equal to 1.
!
        IX = 1
        IY = 1
        IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
        IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
        DO I = 1, N
          DTEMP = DX(IX)
          DX(IX) = DY(IY)
          DY(IY) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
        END DO
      END IF
!
      RETURN
      END SUBROUTINE DSWAP
!
      FUNCTION LSAME (CA, CB) RESULT (IsSAME)
!
!=======================================================================
!                                                                      !
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of   !
!  case.                                                               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  CA      (input) CHARACTER*1                                         !
!  CB      (input) CHARACTER*1                                         !
!          CA and CB specify the single characters to be compared.     !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: CA, CB
!
!  Local variable declarations.
!
      logical :: IsSAME
!
      integer :: INTA, INTB, ZCODE
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test if the characters are equal.
!
      IsSAME = CA .EQ. CB
      IF( IsSAME ) RETURN
!
!  Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR('Z')
!
!  Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!  machines, on which ICHAR returns a value with bit 8 set.
!  ICHAR('A') on Prime machines returns 193 which is the same as
!  ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
!
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
!
!  ASCII is assumed - ZCODE is the ASCII code of either lower or
!  upper case 'Z'.
!
        IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
        IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
!
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
!
!  EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!  upper case 'Z'.
!
        IF (INTA.GE.129 .AND. INTA.LE.137 .OR.                          &
     &      INTA.GE.145 .AND. INTA.LE.153 .OR.                          &
     &      INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
        IF (INTB.GE.129 .AND. INTB.LE.137 .OR.                          &
     &      INTB.GE.145 .AND. INTB.LE.153 .OR.                          &
     &      INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
!
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
!
!  ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!  plus 128 of either lower or upper case 'Z'.
!
        IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
        IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
!
      IsSAME = INTA .EQ. INTB

      RETURN
      END FUNCTION LSAME
!
      SUBROUTINE XERBLA (SRNAME, INFO)
!
!=======================================================================
!                                                                      !
!  XERBLA  is an error handler for the LAPACK routines.                !
!  It is called by an LAPACK routine if an input parameter has an      !
!  invalid value.  A message is printed and execution stops.           !
!                                                                      !
!  Installers may consider modifying the STOP statement in order to    !
!  call system-specific exception-handling facilities.                 !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  SRNAME  (input) CHARACTER*6                                         !
!          The name of the routine which called XERBLA.                !
!                                                                      !
!  INFO    (input) INTEGER                                             !
!          The position of the invalid parameter in the parameter list !
!          of the calling routine.                                     !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: INFO
!
      character (len=*), intent(in) :: SRNAME
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      PRINT 10, SRNAME, INFO
  10  FORMAT (' LAPACK_MOD: ** On entry to ', a, ' parameter number ',  &
     &        i0,' had an illegal value.')
      STOP
!
      END SUBROUTINE XERBLA
!
      END MODULE lapack_mod
