      MODULE round_mod
!
!svn $Id: round.F 841 2017-04-19 21:42:22Z arango $
!====================================================== H. D. Knoble ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Floating point rounding function with a Fuzzy or Tolerant Floor     !
!  function.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     X         Double precision argument to be operated on. It is     !
!                 assumed that X is represented with m mantissa bits.  !
!                                                                      !
!     CT        Comparison Tolerance such that 0 < CT <= 3-SQRT(5)/2.  !
!                 If the relative difference between X and a whole     !
!                 number is less than CT, then TFLOOR is returned as   !
!                 this whole number. By treating the floating-point    !
!                 numbers as a finite ordered set, note that the       !
!                 heuristic EPS=2.**(-(m-1)) and CT=3*eps causes       !
!                 arguments of TFLOOR/TCEIL to be treated as whole     !
!                 numbers if they are exactly whole numbers or are     !
!                 immediately adjacent to whole number representations.!
!                 Since EPS, the  "distance"  between  floating-point  !
!                 numbers on the unit interval, and m, the number of   !
!                 bits in X mantissa, exist on every  floating-point   !
!                 computer, TFLOOR/TCEIL are consistently definable    !
!                 on every floating-point computer.                    !
!                                                                      !
!  Usage:                                                              !
!                                                                      !
!    CT = 3 * EPSILON(X)     That is, CT is about 1 bit on either      !
!                            side of X mantissa bits.                  !
!     Y = round(X, CT)                                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    P. E. Hagerty, 1978: More on Fuzzy Floor and Ceiling, APL QUOTE   !
!      QUAD 8(4):20-24. (The TFLOOR=FL5 took five years of refereed    !
!      evolution publication).                                         !
!                                                                      !
!    L. M. Breed, 1978: Definitions for Fuzzy Floor and Ceiling, APL   !
!        QUOTE QUAD 8(3):16-23.                                        !
!                                                                      !
!  Adapted from H.D. Knoble code (Penn State University).              !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PUBLIC  :: ROUND           ! Tolerant round function
      PUBLIC  :: TCEIL           ! Tolerant ceiling function
      PUBLIC  :: TFLOOR          ! Tolerant floor function
      PRIVATE :: UFLOOR          ! Unfuzzy floor function
!
      CONTAINS
!
!***********************************************************************
      DOUBLE PRECISION FUNCTION ROUND (X, CT) RESULT (Y)
!***********************************************************************
!
!  Imported variable declarations.
!
      real (r8), intent(in) :: X, CT
!
!------------------------------------------------------------------------
!  Compute tolerant round function.
!------------------------------------------------------------------------
!
      Y=TFLOOR(X+0.5_r8,CT)
!
      RETURN
      END FUNCTION ROUND
!
!***********************************************************************
      DOUBLE PRECISION FUNCTION TCEIL (X,CT) RESULT (Y)
!***********************************************************************
!
!  Imported variable declarations.
!
      real (r8), intent(in) :: X, CT
!
!------------------------------------------------------------------------
!  Compute tolerant ceiling function.
!------------------------------------------------------------------------
!
      Y=-TFLOOR(-X,CT)
!
      RETURN
      END FUNCTION TCEIL
!
!***********************************************************************
      DOUBLE PRECISION FUNCTION TFLOOR (X, CT) RESULT (Y)
!***********************************************************************
!
!  Imported variable declarations.
!
      real (r8), intent(in) :: X, CT
!
!  Local variable declarations.
!
      real (r8) :: Q, RMAX, EPS5
!
!------------------------------------------------------------------------
!  Compute tolerant floor function.
!------------------------------------------------------------------------
!
!  Hagerty FL5 function
!
      Q=1.0_r8
      IF (X.lt.0.0_r8) Q=1.0_r8-CT
      RMAX=Q/(2.0_r8-CT)
      EPS5=CT/Q
      Y=UFLOOR(X+MAX(CT,MIN(RMAX,EPS5*ABS(1.0_r8+UFLOOR(X)))))
      IF ((X.le.0.0_r8).or.(Y-X).lt.RMAX) RETURN
      Y=Y-1.0_r8
!
      RETURN
      END FUNCTION TFLOOR
!
!***********************************************************************
      DOUBLE PRECISION FUNCTION UFLOOR (X) RESULT (Y)
!***********************************************************************
!
!  Imported variable declarations.
!
      real (r8), intent(in) :: X
!
!-----------------------------------------------------------------------
!  Compute the largest integer algebraically less than or equal to X;
!  that is, the unfuzzy Floor Function.
!-----------------------------------------------------------------------
!
      Y=X-MOD(X,1.0_r8)-MOD(2.0_r8+SIGN(1.0_r8,X),3.0_r8)
!
      RETURN
      END FUNCTION UFLOOR
!
      END MODULE round_mod


