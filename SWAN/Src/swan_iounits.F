#include "swancpp.h"
      MODULE swan_iounits
!
!svn $Id: swan_iounits.F 755 2008-09-14 19:07:08Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Wname       Wave model stadard input file name.                     !
!                                                                      !
!=======================================================================
!
        implicit none

        integer, parameter :: IOnamesize = 80

#if !defined MODEL_COUPLING
        integer :: NoError = 0
        integer :: exit_flag = 0
        integer, parameter :: stdout = 6
#endif

        integer, parameter :: Ngridss = NestedGrids

#ifdef REFINED_GRID
        integer, dimension(Ngridss) :: mysparent
        integer, dimension(Ngridss) :: myschild
#endif

        character (len=IOnamesize), dimension(Ngridss) :: Wname

      END MODULE swan_iounits
