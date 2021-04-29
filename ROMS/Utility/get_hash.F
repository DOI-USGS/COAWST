#include "cppdefs.h"
      MODULE get_hash_mod
!
!svn $Id: get_hash.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module includes several routines to compute the "checksum" of  !
!  a floating-point array using one of the following methods:          !
!                                                                      !
!    adler32         Fortran, 32-bit Adler algorithm                   !
!    crc32           Fortran, 32-bit Cyclic Redundancy Check algorithm !
!                                                                      !
!  The available methods compute the "checksum" from characters and    !
!  integers. For floating-point data, its values are interpreted as    !
!  unsigned bytes. Here, we have the problem that Fortran does not     !
!  support unsigned integers.  Therefore, the intrinsic function       !
!  TRANSFER is used to convert for 32-bit reals to 32-bit integers.    !
!  We could use the C-binding capabilities of Fortran-2003 to call     !
!  a C-based "checksum" function.  Perhaps, we will do than in the     !
!  future.                                                             !
!                                                                      !
!  The "checksum" value can be used during debugging to compare        !
!  input data solutions from different versions of ROMS when           !
!  implementing new algorithms.                                        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PUBLIC  :: get_hash
      PRIVATE :: adler32
      PRIVATE :: crc32
      PRIVATE :: crc32_ini
!
!  Declare module internal parameters.
!
      integer(i8b) :: crc32_table(0:255)
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_hash (A, Asize, hash)
!***********************************************************************
!
      USE mod_scalars, ONLY : HashMethod, exit_flag
      USE mod_iounits, ONLY : stdout
      USE strings_mod, ONLY : uppercase
!
!  Imported variable declarations.
!
      integer, intent(in) :: Asize
      integer(i8b), intent(out) :: hash
!
      real(r8), intent(in) :: A(:)
!
!  Local variable declarations.
!
      logical, save :: first = .TRUE.
!
!-----------------------------------------------------------------------
!  Compute checksum for the requested floating point vector.
!-----------------------------------------------------------------------
!
      hash=0_i8b
!
      SELECT CASE (uppercase(TRIM(HashMethod)))
        CASE ('ADLER32')
          CALL adler32 (A, Asize, hash)
        CASE ('CRC32')
          IF (first) THEN
            first=.FALSE.
            CALL crc32_ini                  ! compute CRC-32 look table
          END IF
          CALL crc32 (A, Asize, hash)
        CASE DEFAULT
          WRITE (stdout,10) TRIM(HashMethod)
          exit_flag=5
      END SELECT
!
  10  FORMAT (/,' GET_HASH - Illegal checksum method: ',a)
!
      RETURN
      END SUBROUTINE get_hash
!
!***********************************************************************
      SUBROUTINE adler32 (A, Asize, hash)
!***********************************************************************
!                                                                      !
!  Computes the checksum of a 1D array using the 32-bit algorithm from !
!  Mark Adler (Adler-32).                                              !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: Asize

      integer(i8b) :: hash
!
      real(r8), intent(in) :: A(:)
!
!  Local variable declarations.
!
      integer :: alpha, beta, i, j
!
      integer(i8b), parameter :: mod_adler = 65521_i8b

      integer(i8b), allocatable :: Awrk(:)
!
!-----------------------------------------------------------------------
!  Compute ADLER-32 checksum.
!-----------------------------------------------------------------------
!
      alpha=1_i8b
      beta=0_i8b
!
!  Awrk will be an integer array sufficient to hold A(i).
!
      DO i=1,Asize
        Awrk=TRANSFER(A(i), [0])
        DO j=1,SIZE(Awrk)
          alpha=MOD(alpha+Awrk(j), mod_adler)
          beta=MOD(beta+alpha, mod_adler)
        END DO
      END DO
      hash=IOR(beta*65536_i8b, alpha)
!
      RETURN
      END SUBROUTINE adler32
!
!***********************************************************************
      SUBROUTINE crc32 (A, Asize, hash)
!***********************************************************************
!                                                                      !
!  Computes the checksum of a 1D array using the 32-bits (8 bytes)     !
!  cyclic redundancy check (CRC-32) algorithm.                         !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: Asize
      integer(i8b), intent(inout) :: hash
!
      real(r8), intent(in) :: A(:)
!
!  Local variable declarations.
!
      integer :: i
      integer(i8b) :: Ai
!
!-----------------------------------------------------------------------
!  Compute CRC-32 checksum.
!-----------------------------------------------------------------------
!
      hash=NOT(hash)
      DO i=1,Asize
        Ai=TRANSFER(A(i), 1_i8b)      ! 32-bit reals to 32-bit integers
        hash=IEOR(SHIFTR(hash, 8_i8b),                                  &
     &            crc32_table(IAND(IEOR(hash, Ai), 255_i8b)))
      END DO
      hash=NOT(hash)
!
      RETURN
      END SUBROUTINE crc32
!
!***********************************************************************
      SUBROUTINE crc32_ini
!***********************************************************************
!
!  Local variable declarations.
!
      integer :: i, j
      integer(i8b) :: k
!
!-----------------------------------------------------------------------
!  Compute CRC-32 look table.
!-----------------------------------------------------------------------
!
      DO i=0,255
        k=i
        DO j=1,8
          IF (BTEST(k, 0)) THEN
            k=IEOR(SHIFTR(k, 1), -306674912_i8b)
          ELSE
            k=SHIFTR(k, 1_i8b)
          END IF
        END DO
        crc32_table(i)=k
      END DO
!
      RETURN
      END SUBROUTINE crc32_ini
!
      END MODULE get_hash_mod
