#include "cppdefs.h"
#undef ADLER32_CHECKSUM
#undef CRC32_CHECKSUM

      MODULE get_hash_mod
!
!git $Id$
!svn $Id: get_hash.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module includes several routines to compute the "checksum" of  !
!  a floating-point array using one of the following methods:          !
!                                                                      !
#ifdef ADLER32_CHECKSUM
!    adler32       Fortran, 32-bit Adler algorithm                     !
#endif
!    bitsum        Simple bit-by-bit order-invariant sum algorithm     !
#ifdef CRC32_CHECKSUM
!    crc32         Fortran, 32-bit Cyclic Redundancy Check algorithm   !
#endif
!                                                                      !
!  The available methods compute the "checksum" from characters and    !
!  integers. For floating-point data, its values are interpreted as    !
!  unsigned bytes. Here, we have the problem that Fortran does not     !
!  support unsigned integers.  Therefore, the intrinsic function       !
!  TRANSFER is used to convert for 32-bit reals to 32-bit integers.    !
!                                                                      !
!  The "checksum" value can be used during debugging to compare        !
!  input data solutions from different versions of ROMS when           !
!  implementing new algorithms. It is only available for reading       !
!  and writting data in input/output NetCDF files.                     !
!                                                                      !
!  The function "bitsum" is the default method in ROMS since it        !
!  allows tiled I/O data when the PIO library is used. Notice that     !
!  reduction communications are not required with the standard NetCDF  !
!  library since all the data is processed by the master.              !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
#if defined DISTRIBUTE && (defined JEDI || defined PIO_LIB)
      USE distribute_mod, ONLY : mp_reduce
#endif
      USE mod_scalars,    ONLY : HashMethod, exit_flag
      USE mod_iounits,    ONLY : stdout
      USE strings_mod,    ONLY : uppercase
!
      implicit none
!
      PUBLIC  :: get_hash
#ifdef ADLER32_CHECKSUM
      PRIVATE :: adler32
#endif
      PRIVATE :: bitsum
#ifdef CRC32_CHECKSUM
      PRIVATE :: crc32
      PRIVATE :: crc32_ini
!
!  Declare module internal parameters.
!
      integer(i8b) :: crc32_table(0:255)
#endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_hash (A, Asize, hash, Lreduce)
!***********************************************************************
!
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: Lreduce
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
#ifdef ADLER32_CHECKSUM
        CASE ('ADLER32')
          CALL adler32 (A, Asize, hash)
#endif
        CASE ('BITSUM')
          CALL bitsum (A, Asize, hash, Lreduce)
#ifdef CRC32_CHECKSUM
        CASE ('CRC32')
          IF (first) THEN
            first=.FALSE.
            CALL crc32_ini                  ! compute CRC-32 look table
          END IF
          CALL crc32 (A, Asize, hash)
#endif
        CASE DEFAULT
          WRITE (stdout,10) TRIM(HashMethod)
          exit_flag=5
      END SELECT
!
  10  FORMAT (/,' GET_HASH - Illegal checksum method: ',a)
!
      RETURN
      END SUBROUTINE get_hash

#ifdef ADLER32_CHECKSUM
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
#endif
!
!***********************************************************************
      SUBROUTINE bitsum (A, Asize, hash, Lreduce)
!***********************************************************************
!                                                                      !
!  Computes the checksum of a 1D floating-point by casting each value  !
!  to an integer to faciliate the invariant order of the sum in tiled  !
!  parallel applications. A real number can be represented with a set  !
!  64-bit integers (Hallberg and Adcroft, 2014).                       !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Hallberg, R. and A. Adcroft, 2014: An order-invariant real-to-    !
!      integer conversion sum, Parallel Computing, 40, 140-143,        !
!      doi:10.1016/j.parco.2014.04.007.                                !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: Lreduce
!
      integer, intent(in) :: Asize

      integer(i8b) :: hash
!
      real(r8), intent(in) :: A(:)
!
!  Local variable declarations.
!
      integer, parameter :: Ak = KIND(A)
      integer :: i, j
!
      integer(i8b) :: Ac, Asum
#if defined DISTRIBUTE && (defined JEDI || defined PIO_LIB)
      integer(i8b) :: ibuffer(1)
!
      character (len=3) :: op_handle(1)
#endif
!
!-----------------------------------------------------------------------
!  Compute checksum by counting bit-by-bit and summing.
!-----------------------------------------------------------------------
!
!  Here, the "POPCNT" function counts the number of set bits in a
!  machine instruction.  For example, for two 8-bit words operated
!  with XOR, we get
!
!    00100110
!    01100000
!   ----------
!    01000110
!
!    POPCNT(01000110) = 3       'counts the number of bits set to 1'
!
!  The POPCNT is available in all modern Fortran compilers and CPU
!  architectures.
!
      Asum=0_i8b
      DO i=1,Asize
        Ac=POPCNT(TRANSFER(ABS(A(i)), 1_Ak))
        Asum=Asum+Ac
      END DO
      hash=Asum

#if defined DISTRIBUTE && (defined JEDI || defined PIO_LIB)
!
!  If PIO library data processing, sum across all processes. Recall that
!  one or several PETs are reading and or writing data, so we need a
!  global reduction of the checksum, which is order invariant for parallel
!  tiles. The integer arithmetic has no truncation errors. Notice that
!  such reduction is not required in I/0 processed by the master node
!  and Lreduce=.FALSE in such cases.
!
      IF (PRESENT(Lreduce)) THEN
        IF (Lreduce) THEN
          ibuffer(1)=Asum
          op_handle(1)='SUM'
          CALL mp_reduce (1, 1, 1, ibuffer, op_handle)
          hash=ibuffer(1)
        END IF
      END IF
#endif
!
      RETURN
      END SUBROUTINE bitsum

#ifdef CRC32_CHECKSUM
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
#endif
!
      END MODULE get_hash_mod
