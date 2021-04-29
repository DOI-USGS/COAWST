#include "cppdefs.h"
      MODULE mod_kinds
!
!svn $Id: mod_kinds.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
#ifdef ESMF_LIB
        USE ESMF
!
#endif
        implicit none
!
#ifdef ESMF_LIB
!!      integer, parameter :: i1b= ESMF_KIND_I1                !  8-bit
        integer, parameter :: i1b= SELECTED_INT_KIND(1)        !  8-bit
!!      integer, parameter :: i2b= ESMF_KIND_I2                !  8-bit
        integer, parameter :: i2b= SELECTED_INT_KIND(2)        !  8-bit
        integer, parameter :: i4b= ESMF_KIND_I4                ! 16-bit
        integer, parameter :: i8b= ESMF_KIND_I8                ! 32-bit
        integer, parameter :: c8 = SELECTED_REAL_KIND(6,30)    ! 32-bit
        integer, parameter :: dp = ESMF_KIND_R8                ! 64-bit
        integer, parameter :: r4 = ESMF_KIND_R4                ! 32-bit
# ifdef SINGLE_PRECISION
        integer, parameter :: r8 = ESMF_KIND_R4                ! 32-bit
# else
        integer, parameter :: r8 = ESMF_KIND_R8                ! 64-bit
# endif
#else
        integer, parameter :: i1b= SELECTED_INT_KIND(1)        !  8-bit
        integer, parameter :: i2b= SELECTED_INT_KIND(2)        !  8-bit
        integer, parameter :: i4b= SELECTED_INT_KIND(4)        ! 16-bit
        integer, parameter :: i8b= SELECTED_INT_KIND(8)        ! 32-bit
        integer, parameter :: c8 = SELECTED_REAL_KIND(6,30)    ! 32-bit
        integer, parameter :: dp = SELECTED_REAL_KIND(12,300)  ! 64-bit
        integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)    ! 32-bit
# ifdef SINGLE_PRECISION
        integer, parameter :: r8 = SELECTED_REAL_KIND(6,30)    ! 32-bit
# else
        integer, parameter :: r8 = SELECTED_REAL_KIND(12,300)  ! 64-bit
# endif
#endif
#if defined SUN || defined AIX || defined NEC || defined SGI || \
    defined CRAYX1 || defined DEC
        integer, parameter :: r16 = SELECTED_REAL_KIND(24,270) !128-bit
#else
        integer, parameter :: r16 = SELECTED_REAL_KIND(15,300) !128-bit
#endif

      END MODULE mod_kinds
