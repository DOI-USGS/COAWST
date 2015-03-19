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
# ifdef COAWST_COUPLING
      USE mct_coupler_params
# endif
      USE M_COUPLING

      implicit none


      PRIVATE
      PUBLIC  :: SWAN_driver_init
      PUBLIC  :: SWAN_driver_run
      PUBLIC  :: SWAN_driver_finalize

      CONTAINS

      SUBROUTINE SWAN_driver_init (MyCOMM)
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
#ifdef COAWST_COUPLING
      USE WAVES_COUPLER_MOD, ONLY: INITIALIZE_WAV_ROUTERS
#endif
#ifdef NESTING
      USE INTERP_SWAN_MOD, ONLY: swan_ref_init
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: MyCOMM
!     real(m4), intent(in)  :: CouplingTime             ! single precision
!
!  Local variable declarations.
!
      integer :: ng, i, iw, MyRank, MyError
      integer :: ngp, ngc
      character (len=160) :: tempname
#if !defined COAWST_COUPLING
!
!  Establish MPI rank.
!
      CALL mpi_comm_rank (MyCOMM, MyRank, MyError)
!
!     here we need to get name of wname.
!     Read the first file to see how many grids there are.
!
      IF (MyRank.eq.0) CALL getarg (1,tempname)
      i=IOnamesize
      CALL mpi_bcast (tempname, i, MPI_BYTE, 0, MyCOMM, MyError)
!
!     go read the first swan grid to get Num_swan grids
!
      CALL get_numswan_grids(tempname)
      call allocate_swan_iounits
!
!     Now we have the total number of grids so get any more
!     input files.
!
      DO ng=1,NUM_SGRIDS
        IF (MyRank.eq.0) CALL getarg (ng,Wname(ng))
      END DO
      i=IOnamesize*NUM_SGRIDS
      CALL mpi_bcast (Wname, i, MPI_BYTE, 0, MyCOMM, MyError)
#endif
#ifdef NESTING
!
!  Set refined grid parameters.
!
      mysparent(1)=1
      myschild(1)=2
      DO ng=2,NUM_SGRIDS
        mysparent(ng)=ng-1
        myschild(ng)=ng+1
      END DO
      myschild(NUM_SGRIDS)=NUM_SGRIDS
!
! Allocate some arrays for refinement.
!
      IF (.not.allocated(Numspec)) allocate (Numspec(NUM_SGRIDS))
      IF (.not.allocated(ac2size)) allocate (ac2size(NUM_SGRIDS))
      IF (.not.allocated(sstp)) allocate (sstp(NUM_SGRIDS))
      IF (.not.allocated(snew)) allocate (snew(NUM_SGRIDS))
#endif
!
!  Initialize the time counter.
!
      IF (.not.ALLOCATED(iics)) ALLOCATE (iics(NUM_SGRIDS))
      DO iw=1,NUM_SGRIDS
        iics(iw)=0
      END DO
!
!  Initialize restart counters.
!
      DO iw=1,NUM_SGRIDS
        dtswanrst(iw)=0.
        SwanRstFnum(iw)=0
      END DO
!
!  Initialize the grids.
!
      DO ng=1,NUM_SGRIDS
#ifdef NESTING
        ngp=mysparent(ng)
        ngc=myschild(ng)
#else
        ngp=1
        ngc=1
#endif
        CALL SWAN_INITIALIZE (ng, ngc, ngp, NUM_SGRIDS, MyCOMM, Wname(ng))
#ifdef NESTING
!
! Call routine to initialize the child refined grid bc arrays.
!
        IF (ng.gt.1) THEN
          CALL swan_ref_init (ng, ngp)
        END IF
#endif
        CALL SWSYNC
      END DO

#ifdef COAWST_COUPLING
!
!  Initialize the MCT routers to ROMS. This has to be
!  done outside an 'ng' loop so that it is synchronous 
!  with ROMS.
!
      CALL INITIALIZE_WAV_ROUTERS
#endif
!
!  The call to run here does not do a time step, it fills the bc arrays,
!  fill AC2 array of bound spec data for child grids, and enters into MCT.
!
      DO ng=1,NUM_SGRIDS
        CALL SWSYNC
#ifdef NESTING
        ngp=mysparent(ng)
        ngc=myschild(ng)
        sstp(ng)=1
        snew(ng)=3-sstp(ng)
#else
        ngp=1
        ngc=1
#endif
        CALL SWAN_RUN (0, ng, ngc, ngp, 0, NUM_SGRIDS)
      END DO
#ifdef COAWST_COUPLING
!
!     --- couple models during output computations
!
        CALL SWSYNC
        CALL SWAN_CPL (0)
#endif

      RETURN
      END SUBROUTINE SWAN_driver_init

      SUBROUTINE SWAN_driver_run
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
!     real(m4), intent(in)  :: CouplingTime             ! single precision
!
!  Local variable declarations.
!
      integer :: i, iw, ng, MyRank, MyError
      integer, allocatable :: run_grid(:)
      integer, allocatable :: run_cplr(:)
      real :: rtime, rtime_start, cff
      integer :: ngp, ngc
      CALL mpi_comm_rank (WAV_COMM_WORLD, MyRank, MyError)
!
!  Set some initial run time parameters here.
!
      IF (.not.ALLOCATED(run_grid)) ALLOCATE (run_grid(NUM_SGRIDS))
      IF (.not.ALLOCATED(run_cplr)) ALLOCATE (run_cplr(NUM_SGRIDS))
      DO iw=1,NUM_SGRIDS
        run_grid(iw)=1
        run_cplr(iw)=1
      END DO
      rtime_start=0.
      rtime=rtime_start
!
!  Main job control loop here.
!
      DO WHILE (iics(NUM_SGRIDS).lt.MTC_G(NUM_SGRIDS))
!
!  Advance grids in time that have run_grid flag == 1
!  For the first entry, all grids step individual dts.
!
        DO ng=1,NUM_SGRIDS
          IF (run_grid(ng).eq.1) THEN
            iics(ng)=iics(ng)+1
#ifdef NESTING
            ngp=mysparent(ng)
            ngc=myschild(ng)
            sstp(ng)=1+MOD(iics(ng),2)
            snew(ng)=3-sstp(ng)
#else
            ngp=1
            ngc=1
#endif
            CALL SWAN_RUN (iics(ng), ng, ngc, ngp, 1, NUM_SGRIDS)
            run_grid(ng)=0
          END IF
        END DO
!
!  Advance the time counter by the smallest dt.
!
        rtime=rtime+DT_G(NUM_SGRIDS)
!
!  Determine what grids can be time stepped. This is determined
!  by comparing dt(each grid) to global time rtime.
!
        DO ng=1,NUM_SGRIDS
          cff=rtime-rtime_start
          IF (MOD(cff,REAL(DT_G(ng))).eq.0) THEN
            run_grid(ng)=1
          END IF
        END DO
#ifdef COAWST_COUPLING
        IF (run_grid(1).eq.1) THEN
!
!     --- receive data from other coupled models.
          CALL SWSYNC
          CALL SWAN_CPL (1)
         END IF
#endif
      END DO
      IF (ALLOCATED(run_grid)) DEALLOCATE (run_grid)
      IF (ALLOCATED(run_cplr)) DEALLOCATE (run_cplr)
      IF (ALLOCATED(iics)) DEALLOCATE (iics)

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

      DO ng=1,NUM_SGRIDS
        CALL SWAN_FINALIZE (ng, NUM_SGRIDS)
      END DO

      RETURN
      END SUBROUTINE SWAN_driver_finalize

      END MODULE waves_control_mod
