      PROGRAM types
!
!svn $Id: types.F 503 2008-01-10 00:11:51Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This program can be used to check the precision and range of real   !
!  type variables in a particular computer architecture.               !
!                                                                      !
!=======================================================================
!
      implicit none

      integer, parameter :: i = 1

      integer, parameter :: r4  = selected_real_kind(6,300)
      integer, parameter :: r8  = selected_real_kind(12,300)
#if defined SUN || defined AIX
      integer, parameter :: r16 = selected_real_kind(16,300)
#else
      integer, parameter :: r16 = selected_real_kind(15,300)
#endif

      real(r4)  :: a04
      real(r8)  :: a08
      real(r16) :: a16

      real (kind=4)  :: b04
      real (kind=8)  :: b08
#ifdef CRAY
      real (kind=8)  :: b16
#else
      real (kind=16) :: b16
#endif

      WRITE(*,*)
      WRITE(*,*) ' This computer has the following precision and range'
      WRITE(*,*) ' ==================================================='
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) ' TINY Intrinsic Function:'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', TINY(b04)
      WRITE(*,*) ' (kind=8):    ', TINY(b08)
      WRITE(*,*) ' (kind=16):   ', TINY(b16)
      WRITE(*,*) ' K( 6,38 ):   ', TINY(a04)
      WRITE(*,*) ' K(12,300):   ', TINY(a08)
      WRITE(*,*) ' K(16,300):   ', TINY(a16)

      WRITE(*,*)
      WRITE(*,*) ' HUGE Intrinsic Function:'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', HUGE(b04)
      WRITE(*,*) ' (kind=8):    ', HUGE(b08)
      WRITE(*,*) ' (kind=16):   ', HUGE(b16)
      WRITE(*,*) ' K( 6,38 ):   ', HUGE(a04)
      WRITE(*,*) ' K(12,300):   ', HUGE(a08)
      WRITE(*,*) ' K(16,300):   ', HUGE(a16)

      WRITE(*,*)
      WRITE(*,*) ' RANGE Intrinsic Function:'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', RANGE(b04)
      WRITE(*,*) ' (kind=8):    ', RANGE(b08)
      WRITE(*,*) ' (kind=16):   ', RANGE(b16)
      WRITE(*,*) ' K( 6,38 ):   ', RANGE(a04)
      WRITE(*,*) ' K(12,300):   ', RANGE(a08)
      WRITE(*,*) ' K(16,300):   ', RANGE(a16)

      WRITE(*,*)
      WRITE(*,*) ' PRECISION Intrinsic Function:'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', PRECISION(b04)
      WRITE(*,*) ' (kind=8):    ', PRECISION(b08)
      WRITE(*,*) ' (kind=16):   ', PRECISION(b16)
      WRITE(*,*) ' K( 6,38 ):   ', PRECISION(a04)
      WRITE(*,*) ' K(12,300):   ', PRECISION(a08)
      WRITE(*,*) ' K(16,300):   ', PRECISION(a16)

      WRITE(*,*)
      WRITE(*,*) ' EPSILON Intrinsic Function:'
      WRITE(*,*)
      WRITE(*,*) ' (kind=4):    ', EPSILON(b04)
      WRITE(*,*) ' (kind=8):    ', EPSILON(b08)
      WRITE(*,*) ' (kind=16):   ', EPSILON(b16)
      WRITE(*,*) ' K( 6,38 ):   ', EPSILON(a04)
      WRITE(*,*) ' K(12,300):   ', EPSILON(a08)
      WRITE(*,*) ' K(16,300):   ', EPSILON(a16)

      WRITE(*,*)
      WRITE(*,*) ' Integer to Real Conversion:'
      WRITE(*,*)
      WRITE(*,*) ' FLOAT(i):    ', FLOAT(i)
      WRITE(*,*) ' REAL(i):     ', REAL(i)
      WRITE(*,*) ' REAL(i,r8):  ', REAL(i,r8), r4, r8, r16
      WRITE(*,*) ' REAL(i,r16): ', REAL(i,r16)

      STOP

      END PROGRAM types

