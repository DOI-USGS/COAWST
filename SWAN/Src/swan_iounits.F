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

      real, dimension(:), pointer :: dtswanrst

      integer, parameter :: IOnamesize = 160
#ifdef NESTING
      integer, dimension(:), pointer :: mysparent
      integer, dimension(:), pointer :: myschild
#endif
      integer, dimension(:), pointer :: SwanRstFnum

      character (len=IOnamesize), allocatable :: Wname(:)
      character (len=IOnamesize), allocatable :: SwanRstName(:)

      CONTAINS

      SUBROUTINE allocate_swan_iounits
!
!-----------------------------------------------------------------------
!  Allocate I/O files.
!-----------------------------------------------------------------------
!
      integer :: i,j
      character (len=1), parameter :: blank = ' '

#ifdef NESTING
      allocate (mysparent(NUM_SGRIDS))
      allocate (myschild(NUM_SGRIDS))
#endif

      allocate (dtswanrst(NUM_SGRIDS))
      allocate (SwanRstFnum(NUM_SGRIDS))

      allocate ( Wname(NUM_SGRIDS) )
      allocate ( SwanRstName(NUM_SGRIDS) )

      DO j=1,NUM_SGRIDS
        DO i=1,IOnamesize
          Wname(j)(i:i)=blank
          SwanRstName(j)(i:i)=blank
        END DO
      END DO

      RETURN
      END SUBROUTINE allocate_swan_iounits

      END MODULE swan_iounits
