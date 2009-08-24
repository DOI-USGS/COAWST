#include "cppdefs.h"
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
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_iounits
      USE mod_scalars
      USE M_PARALL
      USE M_PARALL2
      USE SWCOMM3
      USE TIMECOMM
      USE M_MPI
!
!  Imported variable declarations.
!
      integer, intent(in) :: MyCOMM
      real(r4), intent(in)  :: CouplingTime             ! single precision
!
!  Local variable declarations.
!
      integer :: i, ng, MyRank, MyError
      integer :: ngp, ngc
#ifdef REFINED_GRID
      integer :: ID
      REAL, ALLOCATABLE :: AC_EXCH(:)
#endif
#ifdef MODEL_COUPLING
!
! Allocate some arrays
!
      IF (.not.allocated(Numspec)) allocate (Numspec(0:Ngrids))
      IF (.not.allocated(ac2size)) allocate (ac2size(Ngrids))
      IF (.not.allocated(BGRIDPG)) allocate (BGRIDPG(Ngrids,100000))
      Numspec(0)=1
      DO ng=1,Ngrids
        BGRIDPG(ng,:)=9999
        Numspec(ng)=46
      END DO
!     Numspec(1)=112
!     Numspec(2)=270
!     Numspec(3)=1  !294
!     Numspec(4)=1

      NBGRPT_last=0
#endif
      ngp=1
      ngc=1
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communicator.
!-----------------------------------------------------------------------
!
!  Initialize the grids.
!
      CALL mpi_comm_rank (MyCOMM, MyRank, MyError)
      DO ng=1,Ngrids
#ifdef REFINED_GRID
        ngp=myparent(ng)
        ngc=mychild(ng)
#endif
        CALL SWAN_INITIALIZE (ng, ngp, Ngrids, MyCOMM, Wname(ng))
!
!  Compute number of boundary arrays to be passed.
!
        ac2size(ng)=MSC*MDC*Numspec(ng)
        CALL SWSYNC
      END DO
!
!  The call to run here does not do a time step, it fills the bc arrays,
!  fill AC2 array of bound spec data for child grids, and enters into MCT.
!
!  AC2_child_G(ng) needs to be the same size as AC2_parent2_G(ngp).
!
      DO ng=1,Ngrids
        CALL SWSYNC
#ifdef REFINED_GRID
        ngp=myparent(ng)
        ngc=mychild(ng)
        IF (ng.gt.1) THEN
          DO i=1,ac2size(ngp)
            PARALL2_MOD(ng)%AC2_child_G(i)=                         &
     &      PARALL2_MOD(ngp)%AC2_parent2_G(i)
           END DO
        END IF
#endif

        CALL SWAN_RUN (CouplingTime, 0, ng, ngp, 0)

#ifdef REFINED_GRID
        IF (ng.lt.Ngrids) THEN
          ALLOCATE ( AC_EXCH(MSC*MDC*Numspec(ng)))
          IF (MyRank.eq.0) THEN
            DO ID=1,ac2size(ng)
              AC_EXCH(ID)=PARALL2_MOD(ng)%AC2_parent2_G(ID)
            END DO
          ENDIF
          CALL MPI_BCAST(AC_EXCH,ac2size(ng),SWREAL,0,MyCOMM,MyError)
          DO ID=1,ac2size(ng)
            PARALL2_MOD(ng)%AC2_parent2_G(ID)=AC_EXCH(ID)
          END DO
          DEALLOCATE (AC_EXCH)
        END IF
#endif
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
      USE mod_param
      USE mod_iounits
      USE mod_scalars
      USE M_PARALL
      USE M_PARALL2
      USE SWCOMM3
      USE TIMECOMM
      USE M_MPI
!
!  Imported variable declarations.
!
      real(r4), intent(in)  :: CouplingTime             ! single precision
!
!  Local variable declarations.
!
      integer :: i, ng, MyRank, MyError
      logical, allocatable :: run_grid(:)
      integer, allocatable :: count(:)
      real :: rtime, rtime_start, rtime_end, cff
      integer :: ngp, ngc
#ifdef REFINED_GRID
      integer :: ID
      REAL, ALLOCATABLE :: AC_EXCH(:)
#endif
      ngp=1
      ngc=1
      CALL mpi_comm_rank (WAV_COMM_WORLD, MyRank, MyError)
!
!  Set some initial run time parameters here.
!
      IF (.not.ALLOCATED(run_grid)) ALLOCATE (run_grid(Ngrids))
      IF (.not.ALLOCATED(count)) ALLOCATE (count(Ngrids))
      DO ng=1,Ngrids
        run_grid(ng)=.TRUE.
        count(ng)=0
      END DO
      rtime_start=0.
      rtime_end=9999999.00
      rtime=rtime_start
!
!  Main job control loop here.
!
!     DO WHILE (rtime.lt.rtime_end)
      DO WHILE (count(Ngrids).lt.MTC_G(Ngrids))
!
!  Advance grids in time that have run_grid flag == True.
!  For the first entry, all grids step individual dts.
!
        DO ng=1,Ngrids
          IF (run_grid(ng).eq..TRUE.) THEN
#ifdef REFINED_GRID
            ngp=myparent(ng)
            ngc=mychild(ng)
            IF (ng.gt.1) THEN
              DO i=1,ac2size(ngp)
                PARALL2_MOD(ng)%AC2_child_G(i)=                         &
     &          PARALL2_MOD(ngp)%AC2_parent2_G(i)
              END DO
            END IF
#endif
            count(ng)=count(ng)+1
            CALL SWAN_RUN (CouplingTime, count(ng), ng, ngp, count(ng))
#ifdef REFINED_GRID
!
!  Get AC data from parent.
!
            IF (ng.lt.Ngrids) THEN
              ALLOCATE ( AC_EXCH(MSC*MDC*Numspec(ng)))
              IF (MyRank.eq.0) THEN
                DO ID=1,ac2size(ng)
                  AC_EXCH(ID)=PARALL2_MOD(ng)%AC2_parent2_G(ID)
                END DO
              ENDIF
              CALL MPI_BCAST(AC_EXCH,ac2size(ng),SWREAL,0,              &
     &                       WAV_COMM_WORLD,MyError)
              DO ID=1,ac2size(ng)
                PARALL2_MOD(ng)%AC2_parent2_G(ID)=AC_EXCH(ID)
              END DO
              DEALLOCATE (AC_EXCH)
            END IF
#endif
          END IF
          run_grid(ng)=.FALSE.
        END DO
!
!  Advance the time counter by the smallest dt.
!
        rtime=rtime+DT_G(Ngrids)
!
!  Determine what grids can be time stepped. This is determined
!  by comparing dt(each grid) to global time rtime.
!
        DO ng=1,Ngrids
          cff=rtime-rtime_start
          IF (MOD(cff,REAL(DT_G(ng))).eq.0) THEN
            run_grid(ng)=.TRUE.
          END IF
        END DO
      END DO

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
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: ng

      DO ng=1,Ngrids
        CALL SWAN_FINALIZE (ng, Ngrids)
      END DO

      RETURN
      END SUBROUTINE SWAN_driver_finalize

      END MODULE waves_control_mod
