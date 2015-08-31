#include "cppdefs.h"
      MODULE inwave_iounits
!
!svn $Id: inwave_iounits.F 755 2008-09-14 19:07:08Z jcwarner $
!================================================== John C. Warner   ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  name       Wave model stadard input file name.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param

      implicit none

      integer, parameter :: IOnamesize = 160
      integer  :: ncSPSid

      integer, allocatable :: ncIWBRYid(:)     ! input boundary conditions for Ac
      integer, allocatable :: ncIWINIid(:)     ! input initial conditions for Ac,cx,cy,ct
      integer, allocatable :: ncIWSWNid(:)     ! input swan 2d spec file

      character (len=IOnamesize), allocatable :: IWBNDname(:)
      character (len=IOnamesize), allocatable :: IWINIname(:)
      character (len=IOnamesize), allocatable :: IWSWNname(:)

      CONTAINS

      SUBROUTINE allocate_inwave_iounits
!
!-----------------------------------------------------------------------
!  Allocate I/O files.
!-----------------------------------------------------------------------
!
      integer :: i,j
      character (len=1), parameter :: blank = ' '

      allocate ( ncIWBRYid(Ngrids) )
      allocate ( ncIWINIid(Ngrids) )
      allocate ( ncIWSWNid(Ngrids) )

      allocate ( IWBNDname(Ngrids) )
      allocate ( IWINIname(Ngrids) )
      allocate ( IWSWNname(Ngrids) )

      DO j=1,Ngrids
        DO i=1,IOnamesize
          IWBNDname(j)(i:i)=blank
          IWINIname(j)(i:i)=blank
          IWSWNname(j)(i:i)=blank
        END DO
      END DO

      RETURN
      END SUBROUTINE allocate_inwave_iounits

      END MODULE inwave_iounits
