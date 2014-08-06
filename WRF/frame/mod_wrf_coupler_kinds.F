      MODULE mod_wrf_coupler_kinds
!
!svn $Id: mod_coupler_kinds.F 1589 2008-07-24 04:02:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
        implicit none
!
        integer, parameter :: m4 = selected_real_kind(6,30)    ! 32-bit
        integer, parameter :: m8 = selected_real_kind(12,300)  ! 64-bit
      END MODULE mod_wrf_coupler_kinds
