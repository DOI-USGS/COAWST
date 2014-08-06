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
      USE M_COUPLING

      implicit none

      integer, parameter :: IOnamesize = 160

!     integer, parameter :: Ngridss = NUM_SGRIDS

#ifdef NESTING
      integer, dimension(:), pointer :: mysparent
      integer, dimension(:), pointer :: myschild
#endif
      character (len=IOnamesize), allocatable :: Wname(:)

      CONTAINS

      SUBROUTINE allocate_swan_iounits
!
!-----------------------------------------------------------------------
!  Allocate I/O files.
!-----------------------------------------------------------------------
!
      integer :: i,j
      character (len=1), parameter :: blank = ' '

      allocate ( Wname(NUM_SGRIDS) )

      DO j=1,NUM_SGRIDS
        DO i=1,IOnamesize
          Wname(j)(i:i)=blank
        END DO
      END DO

#ifdef NESTING
      allocate (mysparent(NUM_SGRIDS))
      allocate (myschild(NUM_SGRIDS))
#endif

      RETURN
      END SUBROUTINE allocate_swan_iounits

      END MODULE swan_iounits
