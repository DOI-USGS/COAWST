#include "cppdefs.h"
      SUBROUTINE ran1_s (harvest)
!
!svn $Id: ran1.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Lagged  Fibonacci  generator  combined with two  Marsaglia  shift   !
!  sequences. On output, returns as HARVEST a uniform random deviate   !
!  between  0.0  and  1.0  (exclusive of the endpoint values).  This   !
!  generator has the same calling and  initialization conventions as   !
!  Fortran 90 RANDOM_NUMBEER routine.  Use RAN_SEED to initialize or   !
!  reinitialize a particular sequence.  The period of this generator   !
!  is about 8.5E+37, and it fully vectorizes.  Validy of the integer   !
!  model assumend by this generator is tested at initialization.       !
!                                                                      !
!  Scalar version adapted from Numerical recipes.                      !
!                                                                      !
!  Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery,    !
!     1996:  Numerical Recipes in Fortran 90,  The Art of Parallel     !
!     Scientific Computing, 2nd Edition, Cambridge Univ. Press.        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE ran_state, ONLY : ran_init
      USE ran_state, ONLY : iran0, jran0, kran0, nran0, mran0
      USE ran_state, ONLY : amm, lenran, rans
!
!  Imported variable declarations.
!
      real(r8), intent(out) :: harvest
!
!-----------------------------------------------------------------------
!  Compute an uniform random deviate (scalar version).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      IF (lenran.lt.1) CALL ran_init (1_i4b)
!
!  Update Fibonacci generator, which has a period p*p+p+1 (p=2^(31)-69).
!
      rans=iran0-kran0
      IF (rans.lt.0) rans=rans+2147483579_i4b
      iran0=jran0
      jran0=kran0
      kran0=rans
!
!  Update Marsaglia shift sequence.
!
      nran0=IEOR(nran0,ISHFT(nran0,13))
      nran0=IEOR(nran0,ISHFT(nran0,-17))
      nran0=IEOR(nran0,ISHFT(nran0,5))
!
!  Once only per cycle, advance sequence by 1, shortening its period to
!  2^(32)-2.
!
      IF (nran0.eq.1) nran0=270369_i4b
!
!  Update Marsaglia shift sequence with perios 2^(32)-1.
!
      mran0=IEOR(mran0,ISHFT(mran0,5))
      mran0=IEOR(mran0,ISHFT(mran0,-13))
      mran0=IEOR(mran0,ISHFT(mran0,6))
!
!  Wrap=around addition.
!
      rans=IEOR(nran0,rans)+mran0
!
!  Make the results positive definite (note that AMM is negative).
!
      harvest=amm*MERGE(rans,NOT(rans), rans<0)

      RETURN
      END SUBROUTINE ran1_s

      SUBROUTINE ran1_v (harvest)
!
!=======================================================================
!                                                                      !
!  Lagged  Fibonacci  generator  combined with two  Marsaglia  shift   !
!  sequences. On output, returns as HARVEST a uniform random deviate   !
!  between  0.0  and  1.0  (exclusive of the endpoint values).  This   !
!  generator has the same calling and  initialization conventions as   !
!  Fortran 90 RANDOM_NUMBEER routine.  Use RAN_SEED to initialize or   !
!  reinitialize a particular sequence.  The period of this generator   !
!  is about 8.5E+37, and it fully vectorizes.  Validy of the integer   !
!  model assumend by this generator is tested at initialization.       !
!                                                                      !
!  Vector version adapted from Numerical recipes.                      !
!                                                                      !
!  Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery,    !
!     1996:  Numerical Recipes in Fortran 90,  The Art of Parallel     !
!     Scientific Computing, 2nd Edition, Cambridge Univ. Press.        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE ran_state, ONLY : ran_init
      USE ran_state, ONLY : iran, jran, kran, nran, mran
      USE ran_state, ONLY : amm, lenran, ranv
!
!  Imported variable declarations.
!
      real(r8), dimension(:), intent(out) :: harvest
!
!  Local variable declarations.
!
      integer(i4b) :: n
!
!-----------------------------------------------------------------------
!  Compute an uniform random deviate (scalar version).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      n=SIZE(harvest)
      IF (lenran.lt.n+1) CALL ran_init (n+1_i4b)
!
!  Update Fibonacci generator, which has a period p*p+p+1 (p=2^(31)-69).
!
      ranv(1:n)=iran(1:n)-kran(1:n)
      WHERE (ranv(1:n).lt.0)                                            &
     &  ranv(1:n)=ranv(1:n)+2147483579_i4b
      iran(1:n)=jran(1:n)
      jran(1:n)=kran(1:n)
      kran(1:n)=ranv(1:n)
!
!  Update Marsaglia shift sequence.
!
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),13))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),-17))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),5))
!
!  Once only per cycle, advance sequence by 1, shortening its period to
!  2^(32)-2.
!
      WHERE (nran(1:n).eq.1)                                            &
     &  nran(1:n)=270369_i4b
!
!  Update Marsaglia shift sequence with perios 2^(32)-1.
!
      mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
      mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
      mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
!
!  Wrap=around addition.
!
      ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
!
!  Make the results positive definite (note that AMM is negative).
!
      harvest=amm*MERGE(ranv(1:n),NOT(ranv(1:n)), ranv(1:n)<0 )

      RETURN
      END SUBROUTINE ran1_v
