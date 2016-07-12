        MODULE nrutil
!
!svn $Id: nrutil.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Numerical Recipies Utility.                                         !
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

      integer(i4b), parameter :: NPAR_ARTH = 16
      integer(i4b), parameter :: NPAR2_ARTH = 8

      INTERFACE array_copy
        MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
      END INTERFACE

      INTERFACE arth
        MODULE PROCEDURE arth_r, arth_d, arth_i
      END INTERFACE

      INTERFACE reallocate
        MODULE PROCEDURE reallocate_rv, reallocate_rm,                  &
     &                   reallocate_iv, reallocate_im,                  &
     &                   reallocate_hv
      END INTERFACE

      INTERFACE gasdev
        SUBROUTINE gasdev_s (harvest)
        USE mod_kinds
        real(r8), intent(out) :: harvest
        END SUBROUTINE gasdev_s
!
        SUBROUTINE gasdev_v (harvest)
        USE mod_kinds
        real(r8), dimension(:), intent(out) :: harvest
        END SUBROUTINE gasdev_v
      END INTERFACE

      INTERFACE ran1
        SUBROUTINE ran1_s (harvest)
        USE mod_kinds
        real(r8), intent(out) :: harvest
        END SUBROUTINE ran1_s
!
        SUBROUTINE ran1_v (harvest)
        USE mod_kinds
        real(r8), dimension(:), intent(out) :: harvest
        END SUBROUTINE ran1_v
      END INTERFACE

      CONTAINS

      SUBROUTINE array_copy_r (src, dest, n_copied, n_not_copied)
!
!=======================================================================
!                                                                      !
!  Copy single precision array where size of source not known in       !
!  advance.                                                            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r4), intent(in) :: src(:)
      real(r4), intent(out) :: dest(:)

      integer(i4b), intent(out) :: n_copied, n_not_copied
!
!-----------------------------------------------------------------------
!  Copy single precision array.
!-----------------------------------------------------------------------
!
      n_copied=MIN(SIZE(src), SIZE(dest))
      n_not_copied=SIZE(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)

      RETURN
      END SUBROUTINE array_copy_r

      SUBROUTINE array_copy_d (src, dest, n_copied, n_not_copied)
!
!=======================================================================
!                                                                      !
!  Copy double precision array where size of source not known in       !
!  advance.                                                            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: src(:)
      real(r8), intent(out) :: dest(:)

      integer(i4b), intent(out) :: n_copied, n_not_copied
!
!-----------------------------------------------------------------------
!  Copy double precision array.
!-----------------------------------------------------------------------
!
      n_copied=MIN(SIZE(src), SIZE(dest))
      n_not_copied=size(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)

      RETURN
      END SUBROUTINE array_copy_d

      SUBROUTINE array_copy_i (src, dest, n_copied, n_not_copied)
!
!=======================================================================
!                                                                      !
!  Copy integer array where size of source not known in advance.       !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer(i4b), intent(in) :: src(:)
      integer(i4b), intent(out) :: dest(:)

      integer(i4b), intent(out) :: n_copied, n_not_copied
!
!-----------------------------------------------------------------------
!  Copy integer array.
!-----------------------------------------------------------------------
!
      n_copied=min(size(src),size(dest))
      n_not_copied=size(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)

      RETURN
      END SUBROUTINE array_copy_i

      FUNCTION arth_r (first, increment, n)
!
!=======================================================================
!                                                                      !
!  Array function returning an arithmetic progression, single          !
!  precision.                                                          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer(i4b), intent(in) :: n

      real(r4), intent(in) :: first, increment
      real(r4), dimension(n) :: arth_r
!
!  Local variable declarations.
!
      integer(i4b) :: k, k2

      real(r4) :: temp
!
!----------------------------------------------------------------------
!  Set arithmetic progression.
!----------------------------------------------------------------------
!
      IF (n.gt.0) arth_r(1)=first
      IF (n.le.NPAR_ARTH) THEN
        DO k=2,n
          arth_r(k)=arth_r(k-1)+increment
        END DO
      ELSE
        DO k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
        END DO
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        DO
          IF (k.ge.n) EXIT
          k2=k+k
          arth_r(k+1:MIN(k2,n))=temp+arth_r(1:MIN(k,n-k))
          temp=temp+temp
          k=k2
        END DO
      END IF

      RETURN
      END FUNCTION arth_r

      FUNCTION arth_d (first, increment, n)
!
!=======================================================================
!                                                                      !
!  Array function returning an arithmetic progression, double          !
!  precision.                                                          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer(i4b), intent(in) :: n

      real(r8), intent(in) :: first, increment
      real(r8), dimension(n) :: arth_d
!
!  Local variable declarations.
!
      integer(i4b) :: k, k2
      real(r8) :: temp
!
!----------------------------------------------------------------------
!  Set arithmetic progression.
!----------------------------------------------------------------------
!
      IF (n.gt.0) arth_d(1)=first
      IF (n.le.NPAR_ARTH) THEN
        DO k=2,n
          arth_d(k)=arth_d(k-1)+increment
        END DO
      ELSE
        DO k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
        END DO
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        DO
          IF (k.ge.n) EXIT
          k2=k+k
          arth_d(k+1:MIN(k2,n))=temp+arth_d(1:MIN(k,n-k))
          temp=temp+temp
          k=k2
        END DO
      END IF

      RETURN
      END FUNCTION arth_d

      FUNCTION arth_i (first, increment, n)
!
!=======================================================================
!                                                                      !
!  Integer array function returning an arithmetic progression.         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer(i4b), intent(in) :: first, increment, n
      integer(i4b), dimension(n) :: arth_i
!
!  Local variable declarations.
!
      integer(i4b) :: k, k2, temp
!
!----------------------------------------------------------------------
!  Set arithmetic progression.
!----------------------------------------------------------------------
!
      IF (n.gt.0) arth_i(1)=first
      IF (n.le.NPAR_ARTH) THEN
        DO k=2,n
          arth_i(k)=arth_i(k-1)+increment
        END DO
      ELSE
        DO k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
        END DO
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        DO
          IF (k.ge.n) EXIT
          k2=k+k
          arth_i(k+1:MIN(k2,n))=temp+arth_i(1:MIN(k,n-k))
          temp=temp+temp
          k=k2
        END DO
      END IF

      RETURN
      END FUNCTION arth_i

      SUBROUTINE nrerror (string)
!
!=======================================================================
!                                                                      !
!  Report an error message and the die.                                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character(len=*), intent(in) :: string
!
!-----------------------------------------------------------------------
!  Report error message to standard output and terminate execution.
!-----------------------------------------------------------------------
!
      PRINT 10, string, 'program terminated by NRERROR'
 10   FORMAT (/,1x,a,/20x,a)
      STOP

      END SUBROUTINE nrerror

      FUNCTION reallocate_rv (p, n)
!
!=======================================================================
!                                                                      !
!  Reallocate a pointer of a single precision vector to a new size,    !
!  preserving its previous content.                                    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r4), pointer :: p(:)
      real(r4), pointer :: reallocate_rv(:)

      integer(i4b), intent(in) :: n
!
!  Local variable declarations.
!
      integer(i4b) :: nold, ierr
!
!-----------------------------------------------------------------------
!  Reallocate pointer for a single precision vector.
!-----------------------------------------------------------------------
!
      ALLOCATE (reallocate_rv(n), STAT=ierr)
      IF (ierr.ne.0)                                                    &
     &  CALL nrerror ('REALLOCATE_RV: error while allocating memory')
      IF (.not.ASSOCIATED(p)) RETURN
      nold=SIZE(p)
      reallocate_rv(1:MIN(nold,n))=p(1:MIN(nold,n))
      DEALLOCATE (p)

      RETURN
      END FUNCTION reallocate_rv

      FUNCTION reallocate_iv (p,n)
!
!=======================================================================
!                                                                      !
!  Reallocate a pointer of a integer vector to a new size, preserving  !
!  its previous content.                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer(i4b), pointer :: p(:)
      integer(i4b), pointer :: reallocate_iv(:)

      integer(i4b), intent(in) :: n
!
!  Local variable declarations.
!
      integer(i4b) :: nold, ierr
!
!-----------------------------------------------------------------------
!  Reallocate pointer for a integer vector.
!-----------------------------------------------------------------------
!
      ALLOCATE (reallocate_iv(n), STAT=ierr)
      IF (ierr.ne.0)                                                    &
     &  CALL nrerror ('REALLOCATE_IV: error while allocating memory')
      IF (.not.ASSOCIATED(p)) RETURN
      nold=SIZE(p)
      reallocate_iv(1:MIN(nold,n))=p(1:MIN(nold,n))
      DEALLOCATE (p)

      RETURN
      END FUNCTION reallocate_iv

      FUNCTION reallocate_hv (p, n)
!
!=======================================================================
!                                                                      !
!  Reallocate a pointer of a character vector to a new size,           !
!  preserving its previous content.                                    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), pointer :: p(:)
      character (len=1), pointer :: reallocate_hv(:)

      integer(i4b), intent(in) :: n
!
!  Local variable declarations.
!
      integer(i4b) :: nold, ierr
!
!-----------------------------------------------------------------------
!  Reallocate pointer for a integer vector.
!-----------------------------------------------------------------------
!
      ALLOCATE (reallocate_hv(n),stat=ierr)
      IF (ierr.ne.0)                                                    &
     &  CALL nrerror ('REALLOCATE_HV: error while allocating memory')
      IF (.not.ASSOCIATED(p)) RETURN
      nold=SIZE(p)
      reallocate_hv(1:MIN(nold,n))=p(1:MIN(nold,n))
      DEALLOCATE (p)

      RETURN
      END FUNCTION reallocate_hv

      FUNCTION reallocate_rm (p, n, m)
!
!=======================================================================
!                                                                      !
!  Reallocate a pointer of a single precision matrix to a new size,    !
!  preserving its previous content.                                    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r4), pointer :: p(:,:)
      real(r4), pointer :: reallocate_rm(:,:)

      integer(i4b), intent(in) :: n, m
!
!  Local variable declarations.
!
      integer(i4b) :: nold, mold, ierr
!
!-----------------------------------------------------------------------
!  Reallocate pointer for a single precision matrix.
!-----------------------------------------------------------------------
!
      ALLOCATE (reallocate_rm(n,m), STAT=ierr)
      IF (ierr.ne.0)                                                    &
     &  CALL nrerror ('REALLOCATE_RM: error while allocating memory')
      IF (.not.ASSOCIATED(p)) RETURN
      nold=SIZE(p,1)
      mold=SIZE(p,2)
      reallocate_rm(1:MIN(nold,n),1:MIN(mold,m))=                       &
     &                                   p(1:MIN(nold,n),1:MIN(mold,m))
      DEALLOCATE (p)

      RETURN
      END FUNCTION reallocate_rm

      FUNCTION reallocate_im (p, n, m)
!
!=======================================================================
!                                                                      !
!  Reallocate a pointer of a integer matrix to a new size, preserving  !
!  its previous content.                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer(i4b), pointer :: p(:,:)
      integer(i4b), pointer :: reallocate_im(:,:)

      integer(i4b), intent(in) :: n, m
!
!  Local variable declarations.
!
      integer(i4b) :: nold, mold, ierr
!
!-----------------------------------------------------------------------
!  Reallocate pointer for a integer matrix.
!-----------------------------------------------------------------------
!
      ALLOCATE (reallocate_im(n,m), STAT=ierr)
      IF (ierr.ne.0)                                                    &
     &  CALL nrerror ('REALLOCATE_IM: error while allocating memory')
      IF (.not.ASSOCIATED(p)) RETURN
      nold=SIZE(p,1)
      mold=SIZE(p,2)
      reallocate_im(1:MIN(nold,n),1:MIN(mold,m))=                       &
     &                                   p(1:MIN(nold,n),1:MIN(mold,m))
      DEALLOCATE (p)

      RETURN
      END FUNCTION reallocate_im

      END MODULE nrutil

