      MODULE mod_kinds
!
!svn $Id: mod_kinds.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
#ifdef ESMF_LIB
        USE ESMF_mod
!
#endif
        implicit none
!
#ifdef ESMF_LIB
        integer, parameter :: i1b= ESMF_KIND_I1                !  1-byte
        integer, parameter :: i2b= ESMF_KIND_I2                !  2-byte
        integer, parameter :: i4b= ESMF_KIND_I4                !  4-byte
        integer, parameter :: i8b= ESMF_KIND_I8                !  8-byte
        integer, parameter :: c8 = selected_real_kind(6,30)    ! 32-bit
        integer, parameter :: r4 = ESMF_KIND_R4                ! 32-bit
        integer, parameter :: r8 = ESMF_KIND_R8                ! 64-bit
#else
        integer, parameter :: i1b= selected_int_kind(2)        !  1-byte
        integer, parameter :: i2b= selected_int_kind(4)        !  2-byte
        integer, parameter :: i4b= selected_int_kind(9)        !  4-byte
        integer, parameter :: c8 = selected_real_kind(6,30)    ! 32-bit
        integer, parameter :: r4 = selected_real_kind(6,30)    ! 32-bit
        integer, parameter :: r8 = selected_real_kind(12,300)  ! 64-bit
#endif
#if defined SUN || defined AIX || defined NEC || defined SGI || \
    defined CRAYX1 || defined DEC
        integer, parameter :: r16 = selected_real_kind(24,270) !128-bit
#else
        integer, parameter :: r16 = selected_real_kind(15,300) !128-bit
#endif
        integer, parameter :: char_len = 256

      END MODULE mod_kinds
