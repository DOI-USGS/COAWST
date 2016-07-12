      PROGRAM types
!
!svn $Id: types.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This program can be used to check the precision and range of real   !
!  type variables in a particular computer architecture.               !
!                                                                      !
!  There are several compilers that do not support REAL*16 data types. !
!  For example PGF90, G95, and GFORTRAN. However, other compilers like !
!  IFORT do have REAL*16 support.                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      integer, parameter :: i = 1

      integer, parameter :: r4  = selected_real_kind(6,30)
      integer, parameter :: r8  = selected_real_kind(12,300)
#ifdef REAL16
# if defined SUN || defined AIX
      integer, parameter :: r16 = selected_real_kind(16,3000)
# else
      integer, parameter :: r16 = selected_real_kind(15,3000)
# endif
#else
# if defined SUN || defined AIX
      integer, parameter :: r16 = selected_real_kind(16,300)
# else
      integer, parameter :: r16 = selected_real_kind(15,300)
# endif
#endif
      real(r4)  :: a04
      real(r8)  :: a08
      real(r16) :: a16

      real (kind=4)  :: b04
      real (kind=8)  :: b08
#ifdef CRAY
      real (kind=8)  :: b16
#else
# ifdef REAL16
      real (kind=16) :: b16
# else
      real(r16) :: b16
# endif
#endif

      WRITE(*,*)
      WRITE(*,*) ' This computer has the following precision and range'
      WRITE(*,*) ' ==================================================='
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) ' TINY Intrinsic Function: smallest positive number'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', TINY(b04)
      WRITE(*,*) ' (kind=8):    ', TINY(b08)
      WRITE(*,*) ' (kind=16):   ', TINY(b16)
      WRITE(*,*) ' K( 6,38 ):   ', TINY(a04)
      WRITE(*,*) ' K(12,300):   ', TINY(a08)
#ifdef REAL16
      WRITE(*,*) ' K(16,3000):  ', TINY(a16)
#else
      WRITE(*,*) ' K(16,300):   ', TINY(a16)
#endif

      WRITE(*,*)
      WRITE(*,*) ' HUGE Intrinsic Function: largest number'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', HUGE(b04)
      WRITE(*,*) ' (kind=8):    ', HUGE(b08)
      WRITE(*,*) ' (kind=16):   ', HUGE(b16)
      WRITE(*,*) ' K( 6,38 ):   ', HUGE(a04)
      WRITE(*,*) ' K(12,300):   ', HUGE(a08)
#ifdef REAL16
      WRITE(*,*) ' K(16,3000):  ', HUGE(a16)
#else
      WRITE(*,*) ' K(16,300):   ', HUGE(a16)
#endif

      WRITE(*,*)
      WRITE(*,*) ' RANGE Intrinsic Function: decimal exponent range'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', RANGE(b04)
      WRITE(*,*) ' (kind=8):    ', RANGE(b08)
      WRITE(*,*) ' (kind=16):   ', RANGE(b16)
      WRITE(*,*) ' K( 6,38 ):   ', RANGE(a04)
      WRITE(*,*) ' K(12,300):   ', RANGE(a08)
#ifdef REAL16
      WRITE(*,*) ' K(16,3000):  ', RANGE(a16)
#else
      WRITE(*,*) ' K(16,300):   ', RANGE(a16)
#endif

      WRITE(*,*)
      WRITE(*,*) ' PRECISION Intrinsic Function: decimal precision'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', PRECISION(b04)
      WRITE(*,*) ' (kind=8):    ', PRECISION(b08)
      WRITE(*,*) ' (kind=16):   ', PRECISION(b16)
      WRITE(*,*) ' K( 6,38 ):   ', PRECISION(a04)
      WRITE(*,*) ' K(12,300):   ', PRECISION(a08)
#ifdef REAL16
      WRITE(*,*) ' K(16,3000):  ', PRECISION(a16)
#else
      WRITE(*,*) ' K(16,300):   ', PRECISION(a16)
#endif

      WRITE(*,*)
      WRITE(*,*) ' EPSILON Intrinsic Function: negligible number, zero'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', EPSILON(b04)
      WRITE(*,*) ' (kind=8):    ', EPSILON(b08)
      WRITE(*,*) ' (kind=16):   ', EPSILON(b16)
      WRITE(*,*) ' K( 6,38 ):   ', EPSILON(a04)
      WRITE(*,*) ' K(12,300):   ', EPSILON(a08)
#ifdef REAL16
      WRITE(*,*) ' K(16,3000):  ', EPSILON(a16)
#else
      WRITE(*,*) ' K(16,300):   ', EPSILON(a16)
#endif

      WRITE(*,*)
      WRITE(*,*) ' Integer (i=1) to Real Conversion:'
      WRITE(*,*)
      WRITE(*,*) ' FLOAT(i):    ', FLOAT(i)
      WRITE(*,*) ' REAL(i):     ', REAL(i)
      WRITE(*,*) ' REAL(i,r8):  ', REAL(i,r8)
      WRITE(*,*) ' REAL(i,r16): ', REAL(i,r16)

      WRITE(*,*)
      WRITE(*,*) ' Values for selected_real_kind parameter:'
      WRITE(*,*)
      WRITE(*,*) ' r4:          ', r4
      WRITE(*,*) ' r8:          ', r8
      WRITE(*,*) ' r16:         ', r16

      STOP

      END PROGRAM types

