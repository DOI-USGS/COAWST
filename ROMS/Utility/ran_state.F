        MODULE ran_state
!
!svn $Id: ran_state.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module supports the random number generation routines. It      !
!  provides each generator with five vectors integers, for use as      !
!  internal random state space.  The  first three integers (iran,      !
!  jran, kran) are maintained as nonnegative  values,  while  the      !
!  last two (mran, nran) have 32-bit nonzero values.  This module      !
!  also includes support for initializing  or  reinitializing the      !
!  state to a desired sequence number, hashing the initial values      !
!  to randon values, and allocating and deallocateing of internal      !
!  workspace.                                                          !
!                                                                      !
!  Adapted from Numerical Recepies.                                    !
!                                                                      !
!  Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery,    !
!     1996:  Numerical Recipes in Fortran 90,  The Art of Parallel     !
!     Scientific Computing, 2nd Edition, Cambridge Univ. Press.        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PUBLIC

      integer(i4b), parameter :: hg = HUGE(1_i4b)
      integer(i4b), parameter :: hgm = -hg
      integer(i4b), parameter :: hgng = hgm - 1

      integer(i4b), save :: lenran = 0
      integer(i4b), save :: seq = 0
      integer(i4b), save :: iran0
      integer(i4b), save :: jran0
      integer(i4b), save :: kran0
      integer(i4b), save :: nran0
      integer(i4b), save :: mran0
      integer(i4b), save :: rans

      integer(i4b), pointer, save :: iran(:)
      integer(i4b), pointer, save :: jran(:)
      integer(i4b), pointer, save :: kran(:)
      integer(i4b), pointer, save :: nran(:)
      integer(i4b), pointer, save :: mran(:)
      integer(i4b), pointer, save :: ranv(:)

      integer(i4b), pointer, save :: ranseeds(:,:)

      real(r8), save :: amm
!
      INTERFACE ran_hash
        MODULE PROCEDURE ran_hash_s, ran_hash_v
      END INTERFACE
!
      CONTAINS

      SUBROUTINE ran_init (length)
!
!=======================================================================
!                                                                      !
!  This routine initializes or reinitializes the random generator      !
!  state space to vectors of size LENGTH.  The saved variable SEQ      !
!  is hashed (via a call to RAN_HASH) to create  unique  starting      !
!  seeds, different for each vector component.                         !
!                                                                      !
!=======================================================================
!
      USE nrutil, ONLY : arth, nrerror, reallocate
!
!  Imported variable declarations.
!
      integer(i4b), intent(in) :: length
!
!  Local variable declarations.
!
      integer(i4b) :: hgt, j, new, sz
!
!-----------------------------------------------------------------------
!  Initialize randon number generator vectors.
!-----------------------------------------------------------------------
!
      IF (length.lt.lenran) RETURN
      hgt=hg
!
!  Check that kind value I4B is in fact a 32-bit integer with the usual
!  properties that we expect it to have (under negation and wrap-around
!  addition). If all these test are satisfied, then the routines that
!  use this module are portable, even though they go beyond F90 integer
!  model.
!
      IF (hg.ne.2147483647)                                             &
     &  CALL nrerror ('RAN_INIT: arith assump 1 fails')
      IF (hgng.ge.0)                                                    &
     &  CALL nrerror ('RAN_INIT: arith assump 2 fails')
      IF ((hgt+1).ne.hgng)                                              &
     &  CALL nrerror ('RAN_INIT: arith assump 3 fails')
      IF (NOT(hg).ge.0)                                                 &
     &  CALL nrerror ('RAN_INIT: arith assump 4 fails')
      IF (NOT(hgng).lt.0)                                               &
     &  CALL nrerror ('RAN_INIT: arith assump 5 fails')
      IF ((hg+hgng).ge.0)                                               &
     &  CALL nrerror ('RAN_INIT: arith assump 6 fails')
      IF (NOT(-1_i4b).lt.0)                                             &
     &  CALL nrerror ('RAN_INIT: arith assump 7 fails')
      IF (NOT(0_i4b).ge.0)                                              &
     &  CALL nrerror ('RAN_INIT: arith assump 8 fails')
      IF (NOT(1_i4b).ge.0)                                              &
     &  CALL nrerror ('RAN_INIT: arith assump 9 fails')
!
!  Reallocate or allocate state space.
!
      IF (lenran.gt.0) THEN
        ranseeds => reallocate (ranseeds, length, 5_i4b)
        ranv => reallocate (ranv, length-1_i4b)
        new=lenran+1
      ELSE
        ALLOCATE (ranseeds(length,5))
        ALLOCATE (ranv(length-1))
        new=1
        amm=NEAREST(1.0_r8,-1.0_r8)/hgng
        IF ((amm*hgng.ge.1.0_r8).or.(amm*hgng.le.0.0_r8))               &
     &    CALL nrerror ('RAN_INIT: arth assump 10 fails')
      END IF
!
!  Set starting values, unique by SEQ and vector component.
!
      ranseeds(new:,1)=seq
      sz=SIZE(ranseeds(new:,1))
      ranseeds(new:,2:5)=SPREAD(arth(new,1_i4b,sz),2,4)
!
!  Hash them.
!
      DO j=1,4
        CALL ran_hash (ranseeds(new:,j), ranseeds(new:,j+1))
      END DO
!
!  Enforce nonnegativity.
!
      WHERE (ranseeds(new:,1:3).lt.0)                                   &
     &  ranseeds(new:,1:3)=NOT(ranseeds(new:,1:3))
!
!  Enforce nonzero.
!
      WHERE (ranseeds(new:,4:5).eq.0)                                   &
     &   ranseeds(new:,4:5)=1
!
!  Set scalar seeds.
!
      IF (new.eq.1) THEN
        iran0=ranseeds(1,1)
        jran0=ranseeds(1,2)
        kran0=ranseeds(1,3)
        mran0=ranseeds(1,4)
        nran0=ranseeds(1,5)
        rans=nran0
      END IF
!
!  Point to vector seeds.
!
      IF (length.gt.1) THEN
        iran => ranseeds(2:,1)
        jran => ranseeds(2:,2)
        kran => ranseeds(2:,3)
        mran => ranseeds(2:,4)
        nran => ranseeds(2:,5)
        ranv = nran
      END IF
      lenran=length

      END SUBROUTINE ran_init

      SUBROUTINE ran_deallocate
!
!=======================================================================
!                                                                      !
!  User interface to release the workspace used by random number       !
!  routines.                                                           !
!                                                                      !
!=======================================================================
!
      IF (lenran.gt.0) THEN
        DEALLOCATE (ranseeds, ranv)
        NULLIFY (ranseeds, ranv, iran, jran, kran, mran, nran)
        lenran=0
      END IF
      END SUBROUTINE ran_deallocate

      SUBROUTINE ran_seed (sequence, size, put, get)
!
!=======================================================================
!                                                                      !
!  User interface for seeding the random number routines.  Syntax is   !
!  exactly like Fortran 90 RANDOM_SEED, with one additional argument   !
!  keyword: SEQUENCE, set to any integer value,  causes an immediate   !
!  new initialization, seeded by that integer.                         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, optional, intent(in) :: sequence
      integer, optional, intent(out) :: size

      integer, optional, intent(in) :: put(:)
      integer, optional, intent(out) :: get(:)
!
!-----------------------------------------------------------------------
!  Set random number seeds.
!-----------------------------------------------------------------------
!
      IF (PRESENT(size)) THEN
        size=5*lenran
      ELSE IF (PRESENT(put)) THEN
        IF (lenran.eq.0) RETURN
        ranseeds=RESHAPE(put,SHAPE(ranseeds))
        WHERE (ranseeds(:,1:3).lt.0)                                    &
     &    ranseeds(:,1:3)=NOT(ranseeds(:,1:3))
        WHERE (ranseeds(:,4:5).eq.0)                                    &
     &    ranseeds(:,4:5)=1
          iran0=ranseeds(1,1)
          jran0=ranseeds(1,2)
          kran0=ranseeds(1,3)
          mran0=ranseeds(1,4)
          nran0=ranseeds(1,5)
      ELSE IF (present(get)) THEN
        IF (lenran.eq.0) RETURN
        ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
        get=RESHAPE(ranseeds,SHAPE(get))
      ELSE IF (PRESENT(sequence)) THEN
        CALL ran_deallocate
        seq=sequence
      END IF

      RETURN
      END SUBROUTINE ran_seed

      SUBROUTINE ran_hash_s (il, ir)
!
!=======================================================================
!                                                                      !
!  DES-like hashing of 32-bit integer, using shifts, xor, and adds to  !
!  make the interval nonlinear function. Scalar version.               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer(i4b), intent(inout) :: il, ir
!
!  Local variable declarations.
!
      integer(i4b) :: is, j
!
!-----------------------------------------------------------------------
!  Bit mixing. The various constants should not be changed.
!-----------------------------------------------------------------------
!
      DO j=1,4
        is=ir
        ir=IEOR(ir,ISHFT(ir,5))+1422217823
        ir=IEOR(ir,ISHFT(ir,-16))+1842055030
        ir=IEOR(ir,ISHFT(ir,9))+80567781
        ir=IEOR(il,ir)
        il=is
      END DO

      RETURN
      END SUBROUTINE ran_hash_s

      SUBROUTINE ran_hash_v (il, ir)
!
!=======================================================================
!                                                                      !
!  DES-like hashing of 32-bit integer, using shifts, xor, and adds to  !
!  make the interval nonlinear function. Vector version.               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer(i4b), intent(inout) :: il(:)
      integer(i4b), intent(inout) :: ir(:)
!
!  Local variable declarations.
!
      integer(i4b) :: j
      integer(i4b), dimension(SIZE(il)) :: is
!
!-----------------------------------------------------------------------
!  Bit mixing. The various constants should not be changed.
!-----------------------------------------------------------------------
!
      DO j=1,4
        is=ir
        ir=IEOR(ir,ISHFT(ir,5))+1422217823
        ir=IEOR(ir,ISHFT(ir,-16))+1842055030
        ir=IEOR(ir,ISHFT(ir,9))+80567781
        ir=IEOR(il,ir)
        il=is
      END DO

      RETURN
      END SUBROUTINE ran_hash_v

      END MODULE ran_state

