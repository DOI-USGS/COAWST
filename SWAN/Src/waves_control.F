#include "swancpp.h"
      MODULE waves_control_mod
!
!svn $Id: waves_control.F 814 2008-10-29 01:42:17Z jcwarner $
!================================================== John C. Warner  ====
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  SWAN model:                                                         !
!                                                                      !
!  This driver executes SWAN by controlling initialization,            !
!  time-stepping, and finalization for nested grids.                   !
!                                                                      !
!     SWAN_initialize                                                  !
!     SWAN_run                                                         !
!     SWAN_finalize                                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: SWAN_driver_init
      PUBLIC  :: SWAN_driver_run
      PUBLIC  :: SWAN_driver_finalize

      CONTAINS

      SUBROUTINE SWAN_driver_init (MyCOMM, CouplingTime)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes SWAN variables               !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE swan_iounits
      USE mod_coupler_kinds
      USE M_PARALL
      USE M_PARALL2
      USE SWCOMM3
      USE TIMECOMM
      USE M_MPI
#ifdef REFINED_GRID
      USE INTERP_SWAN_MOD, ONLY: swan_ref_init
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: MyCOMM
      real(m4), intent(in)  :: CouplingTime             ! single precision
!
!  Local variable declarations.
!
      integer :: ng, i, MyRank, MyError
      integer :: ngp, ngc
#ifdef REFINED_GRID
!
!  Set refined grid parameters.
!
      mysparent(1)=1
      myschild(1)=2
      DO ng=2,Ngridss
        mysparent(ng)=ng-1
        myschild(ng)=ng+1
      END DO
      myschild(Ngridss)=Ngridss
!
! Allocate some arrays for refinement.
!

      IF (.not.allocated(Numspec)) allocate (Numspec(Ngridss))
      IF (.not.allocated(ac2size)) allocate (ac2size(Ngridss))
      IF (.not.allocated(sstp)) allocate (sstp(Ngridss))
      IF (.not.allocated(snew)) allocate (snew(Ngridss))
#endif
#if !defined MODEL_COUPLING
!
!  Establish MPI rank.
!
      CALL mpi_comm_rank (MyCOMM, MyRank, MyError)
!
!     here we need to get name of wname.
!
      DO ng=1,Ngridss
        IF (MyRank.eq.0) CALL getarg (ng,Wname(ng))
      END DO
      i=IOnamesize*ng
      CALL mpi_bcast (Wname, i, MPI_BYTE, 0, MyCOMM, MyError)
#endif
!
!  Initialize the grids.
!
      DO ng=1,Ngridss
#ifdef REFINED_GRID
        ngp=mysparent(ng)
        ngc=myschild(ng)
#else
        ngp=1
        ngc=1
#endif
        CALL SWAN_INITIALIZE (ng, ngc, ngp, Ngridss, MyCOMM, Wname(ng))
#ifdef REFINED_GRID
!
! Call routine to initialize the child refined grid bc arrays.
!
        IF (ng.gt.1) THEN
          CALL swan_ref_init (ng, ngp)
        END IF
#endif
        CALL SWSYNC
      END DO
!
!  The call to run here does not do a time step, it fills the bc arrays,
!  fill AC2 array of bound spec data for child grids, and enters into MCT.
!
      DO ng=1,Ngridss
        CALL SWSYNC
#ifdef REFINED_GRID
        ngp=mysparent(ng)
        ngc=myschild(ng)
        sstp(ng)=1
        snew(ng)=3-sstp(ng)
#else
        ngp=1
        ngc=1
#endif
!
        CALL SWAN_RUN (CouplingTime, 0, ng, ngc, ngp, 0, Ngridss)
!
      END DO

      RETURN
      END SUBROUTINE SWAN_driver_init

      SUBROUTINE SWAN_driver_run (CouplingTime)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE swan_iounits
      USE mod_coupler_kinds
      USE M_PARALL
      USE M_PARALL2
      USE SWCOMM3
      USE TIMECOMM
      USE M_MPI
!
!  Imported variable declarations.
!
      real(m4), intent(in)  :: CouplingTime             ! single precision
!
!  Local variable declarations.
!
      integer :: i, ng, MyRank, MyError
      integer, allocatable :: run_grid(:)
      integer, allocatable :: count(:)
      real :: rtime, rtime_start, cff
      integer :: ngp, ngc
      CALL mpi_comm_rank (WAV_COMM_WORLD, MyRank, MyError)
!
!  Set some initial run time parameters here.
!
      IF (.not.ALLOCATED(run_grid)) ALLOCATE (run_grid(Ngridss))
      IF (.not.ALLOCATED(count)) ALLOCATE (count(Ngridss))
      DO ng=1,Ngridss
        run_grid(ng)=1
        count(ng)=0
      END DO
      rtime_start=0.
      rtime=rtime_start
!
!  Main job control loop here.
!
      DO WHILE (count(Ngridss).lt.MTC_G(Ngridss))
!
!  Advance grids in time that have run_grid flag == 1
!  For the first entry, all grids step individual dts.
!
        DO ng=1,Ngridss
          IF (run_grid(ng).eq.1) THEN
            count(ng)=count(ng)+1
#ifdef REFINED_GRID
            ngp=mysparent(ng)
            ngc=myschild(ng)
            sstp(ng)=1+MOD(count(ng),2)
            snew(ng)=3-sstp(ng)
#else
            ngp=1
            ngc=1
#endif
            CALL SWAN_RUN (CouplingTime, count(ng), ng, ngc, ngp,       &
     &                     1, Ngridss)
            run_grid(ng)=0
          END IF
        END DO
!
!  Advance the time counter by the smallest dt.
!
        rtime=rtime+DT_G(Ngridss)
!
!  Determine what grids can be time stepped. This is determined
!  by comparing dt(each grid) to global time rtime.
!
        DO ng=1,Ngridss
          cff=rtime-rtime_start
          IF (MOD(cff,REAL(DT_G(ng))).eq.0) THEN
            run_grid(ng)=1
          END IF
        END DO
      END DO
      IF (ALLOCATED(run_grid)) DEALLOCATE (run_grid)
      IF (ALLOCATED(count)) DEALLOCATE (count)

      RETURN
      END SUBROUTINE SWAN_driver_run

      SUBROUTINE SWAN_driver_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates SWAN.                                       !
!                                                                      !
!=======================================================================
!
      USE swan_iounits
!
!  Local variable declarations.
!
      integer :: ng

      DO ng=1,Ngridss
        CALL SWAN_FINALIZE (ng, Ngridss)
      END DO

      RETURN
      END SUBROUTINE SWAN_driver_finalize

      END MODULE waves_control_mod
