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
      PUBLIC  :: SWAN_driver
      PUBLIC  :: SWAN_driver_run
      PUBLIC  :: SWAN_driver_finalize

      CONTAINS

      SUBROUTINE SWAN_driver (MyCOMM, CouplingTime)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!     USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_coupler
      USE M_GENARR
      USE M_PARALL
      USE M_PARALL2
      USE SWCOMM3
      USE TIMECOMM
      USE M_MPI
!     USE M_COUPLING
!
!  Imported variable declarations.
!
      integer, intent(in) :: MyCOMM
      real(r4), intent(in)  :: CouplingTime             ! single precision
!
!  Local variable declarations.
!
      integer :: i, inp, ModelId, ng
      integer :: MyColor, MyKey, MyRank
      INTEGER :: MySize, MyError, MyStatus, MyHandle
      logical, allocatable :: run_grid(:)
      integer, allocatable :: count(:)
      real :: rtime, rtime_start, rtime_end, cff
#ifdef REFINED_GRID
      integer :: ngp, ngc
#endif
      MyKey=0
      inp=22
#ifdef REFINED_GRID
!
! Allocate some arrays
!
      IF (.not.allocated(Numspec)) allocate (Numspec(Ngrids))
      IF (.not.allocated(ac2size)) allocate (ac2size(Ngrids))
!     IF (.not.allocated(BGRIDPL)) allocate (BGRIDPL(Ngrids,100000))
      IF (.not.allocated(BGRIDPG)) allocate (BGRIDPG(Ngrids,100000))
      DO ng=1,Ngrids
        BGRIDPG(ng,:)=9999
      END DO
      Numspec(1)=46
      Numspec(2)=46
      NBGRPT_last=0
#endif
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communicator.
!-----------------------------------------------------------------------
!
      CALL mpi_comm_rank (MyCOMM, MyRank, MyError)
      CALL mpi_comm_size (MyCOMM, MySize, MyError)
!
!  Initialize the grids.
!
      DO ng=1,Ngrids
        CALL SWAN_INITIALIZE (ng, Ngrids, MyCOMM, Wname(ng))
!
!  Compute number of boundary arrays to be passed.
!
#ifdef REFINED_GRID
        ac2size(ng)=MSC*MDC*Numspec(ng)
#endif
      END DO
!
!  The call to run here does not do a time step, it fills the bc arrays,
!  fill AC2 array of bound spec data for child grids, and enters into MCT.
!
      DO ng=1,Ngrids

#ifdef REFINED_GRID
            IF (ng.gt.1) THEN
              ngp=myparent(ng)
              ngc=mychild(ng)
              DO i=1,ac2size(ng)
                PARALL2_MOD(ng)%AC2_parent2_G(i)=                       &
     &          PARALL2_MOD(ngp)%AC2_parent2_G(i)
              END DO
            END IF
#endif

        CALL SWAN_RUN (CouplingTime, 0, ng, 0)

#ifdef REFINED_GRID
        IF (ng.lt.Ngrids) THEN
          CALL MPI_BCAST(PARALL2_MOD(ng)%AC2_parent2_G,                 &
     &      ac2size(ng),SWREAL,0,MyCOMM,MyError)
        END IF
#endif
      END DO
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
      DO WHILE (rtime.lt.rtime_end)
!
!  Advance grids in time that have run_grid flag == True.
!  For the first entry, all grids step individual dts.
!
        DO ng=1,Ngrids
          IF (run_grid(ng).eq..TRUE.) THEN
#ifdef REFINED_GRID
            IF (ng.gt.1) THEN
              ngp=myparent(ng)
              ngc=mychild(ng)
              DO i=1,ac2size(ng)
                PARALL2_MOD(ng)%AC2_parent2_G(i)=                          &
     &          PARALL2_MOD(ngp)%AC2_parent2_G(i)
              END DO
            END IF
#endif
            count(ng)=count(ng)+1
            CALL SWAN_RUN (CouplingTime, count(ng), ng, count(ng))
!
!  Get AC data from parent.
!
#ifdef REFINED_GRID
            IF (ng.lt.Ngrids) THEN
              CALL MPI_BCAST(PARALL2_MOD(ng)%AC2_parent2_G,             &
     &                         ac2size(ng),SWREAL,0,MyCOMM,MyError)
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
      END SUBROUTINE SWAN_driver

      SUBROUTINE SWAN_driver_run (CouplingTime)
!
!=======================================================================
!                                                                      !
!  This routine runs SWAN.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef VERIFICATION
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_scalars
      USE mod_coupler
      USE M_PARALL
!
!  Imported variable declarations.
!
      real(r4), intent(in)  :: CouplingTime             ! single precision
!
!  Local variable declarations.
!
      integer :: ng, my_iic, MyError
      integer :: count1, MyCOMM_wave
#if defined REFINED_GRID
      integer :: count2
#endif
!
!-----------------------------------------------------------------------
!  Run model for all nested grids, if any.
!-----------------------------------------------------------------------

      NL_LOOP : DO count1=0,1


        IF (exit_flag.ne.NoError) THEN
          IF (MyRank.eq.0) THEN
            WRITE (stdout,'(/,a,i3,/)') Rerror(exit_flag), exit_flag
          END IF
          RETURN
        END IF

      END DO NL_LOOP


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
        CALL SWAN_FINALIZE
      END DO

      RETURN
      END SUBROUTINE SWAN_driver_finalize

      END MODULE waves_control_mod
