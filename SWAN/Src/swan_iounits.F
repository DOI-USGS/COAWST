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


        integer, parameter :: Ngridss = NestedGrids

#ifdef REFINED_GRID
        integer, dimension(Ngridss) :: mysparent
        integer, dimension(Ngridss) :: myschild
#endif

        character (len=IOnamesize), dimension(Ngridss) :: Wname

      END MODULE swan_iounits
