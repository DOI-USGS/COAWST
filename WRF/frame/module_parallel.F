      MODULE module_parallel
!
!svn $Id:$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains all variables used for parallelization         !
!                                                                      !
!=======================================================================
!
      implicit none
!
!  Ocean model MPI group communicator handle.
!
      integer :: MyComm
      integer :: ATM_COMM_WORLD

!     real, pointer :: hwave(:,:)
!     real, pointer :: lwave(:,:)

      CONTAINS

      SUBROUTINE allocate_mctvars (grid)
      USE module_domain

      implicit none
      TYPE(domain) , INTENT (IN) :: grid 
      integer :: i, j, is, ie, js, je, ij, ierr

      is = grid%sm31
      ie = grid%em31
      js = grid%sm33
      je = grid%em33

!     ALLOCATE(hwave(is:ie,js:je),STAT=ierr)
!     ALLOCATE(lwave(is:ie,js:je),STAT=ierr)

!      if (ierr.ne.0) then
!        CALL wrf_error_fatal (                                          &
!    &               'frame/module_parallel.F: Failed to allocate HWAVE')
!      endif

      DO j=js,je
        DO i=is,ie
!         hwave(i,j)=0.001
!         lwave(i,j)=0.001
        END DO
      END DO


      END SUBROUTINE allocate_mctvars



      END MODULE module_parallel

