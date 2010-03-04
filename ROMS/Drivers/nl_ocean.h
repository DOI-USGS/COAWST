      MODULE ocean_control_mod
!
!svn $Id: nl_ocean.h 814 2008-10-29 01:42:17Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Nonlinear model:                                          !
!                                                                      !
!  This driver executes ROMS/TOMS standard nonlinear model.  It        !
!  controls the initialization, time-stepping, and finalization        !
!  of the nonlinear model execution following ESMF conventions:        !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize

      CONTAINS

      SUBROUTINE ROMS_initialize (first, MyCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
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
!
#ifdef AIR_OCEAN 
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_coupling
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_coupling
#endif
#ifdef INWAVE_COUPLING
      USE driver_inwave_mod, ONLY : inwave_init
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: MyCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: ng, thread

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(MyCOMM)) THEN
        OCN_COMM_WORLD=MyCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, MySize, MyError)
#endif
!
!-----------------------------------------------------------------------
!  On first pass, initialize model parameters a variables for all
!  nested/composed grids.  Notice that the logical switch "first"
!  is used to allow multiple calls to this routine during ensemble
!  configurations.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
!
!  Initialize parallel parameters.
!
        CALL initialize_parallel
!
!  Initialize wall clocks.
!
        IF (Master) THEN
          WRITE (stdout,10)
        END IF
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            CALL wclock_on (ng, iNLM, 0)
          END DO
!$OMP END PARALLEL DO
        END DO
!
!  Read in model tunable parameters from standard input. Initialize
!  "mod_param", "mod_ncparam" and "mod_scalar" modules.
!
        CALL inp_par (iNLM)
        IF (exit_flag.ne.NoError) RETURN
!
!  Allocate and initialize modules variables.
!
        CALL mod_arrays (allocate_vars)
!
#if defined MCT_LIB && (defined AIR_OCEAN || defined WAVES_OCEAN)
!
!  Initialize coupling streams between model(s).
!
        DO ng=1,Ngrids
# ifdef WAVES_OCEAN
          CALL initialize_ocn2wav_coupling (ng, MyRank)
# endif
# ifdef AIR_OCEAN
          CALL initialize_ocn2atm_coupling (ng, MyRank)
# endif
        END DO
#endif
#ifdef VERIFICATION
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Initialize model state variables for all nested/composed grids.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL initial (ng)
        IF (exit_flag.ne.NoError) RETURN
#ifdef DISTRIBUTE
        CALL mpi_barrier(OCN_COMM_WORLD, MyError)
#endif
      END DO
#ifdef INWAVE_COUPLING
        DO ng=1,Ngrids
	CALL inwave_init (ng, TILE)
        END DO
#endif
!
!  Initialize run or ensemble counter.
!
      Nrun=1

#ifdef VERIFICATION
!
!  Create out NetCDF file containing model solution at observation
!  locations.
!
      IF (Nrun.eq.1) THEN
        DO ng=1,Ngrids
          LdefMOD(ng)=.TRUE.
          wrtNLmod(ng)=.TRUE.
          CALL def_mod (ng)
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END IF
#endif

#ifdef IOM
!
!  Set the current outer loop iteration.
!
      outer=Nouter
      IF (Master) THEN
        WRITE (stdout,'(/,a,i3,/)')                                     &
     &        'NL ROMS/TOMS: Outer Loop Iteration = ', outer
      END IF
#endif
!
!  Substract a time-step to model time after initialization because the
!  main time-stepping driver always add a single time-step.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,20) ng, ntstart(ng), ntend(ng)
        END IF
        time(ng)=time(ng)-dt(ng)
      END DO

 10   FORMAT (' Process Information:',/)
 20   FORMAT ('NL ROMS/TOMS: started time-stepping:',                   &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')')

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Tend)
!
!=======================================================================
!                                                                      !
!  This routine runs ROMS/TOMS nonlinear model from specified starting !
!  (Tstr) to ending (Tend) time-steps.                                 !
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
#ifdef INWAVE_COUPLING
      USE driver_inwave_mod, ONLY : inwave_run
#endif
!
!  Imported variable declarations.
!
      integer, dimension(Ngrids) :: Tstr   ! starting time-step
      integer, dimension(Ngrids) :: Tend   ! ending   time-step
!
!  Local variable declarations.
!
	integer :: ng, my_iic, MyError
#if defined REFINED_GRID
      logical, allocatable :: run_grid(:)
      integer, allocatable :: count(:)
      real(r8) :: rtime, rtime_start, rtime_end, cff
#endif
!
!-----------------------------------------------------------------------
!  Run model for all nested grids, if any.
!-----------------------------------------------------------------------
!

#if defined REFINED_GRID

      IF (.not.ALLOCATED(run_grid)) ALLOCATE (run_grid(Ngrids))
      IF (.not.ALLOCATED(count)) ALLOCATE (count(Ngrids))
      DO ng=1,Ngrids
        run_grid(ng)=.TRUE.
!       count(ng)=my_iic-Tstr(1)
        count(ng)=0
      END DO
      rtime_start=0.
      rtime_end=9999999.00
      rtime=rtime_start
!
!  Main job control loop here.
!
        DO WHILE (count(Ngrids).lt.Tend(Ngrids))
!
!  Advance grids in time that have run_grid flag == True.
!  For the first entry, all grids step individual dts.
!
        DO ng=1,Ngrids
          IF (run_grid(ng).eq..TRUE.) THEN
            count(ng)=count(ng)+1
            iic(ng)=count(ng)
# ifdef SOLVE3D
            CALL main3d (ng)
# else
            CALL main2d (ng)
# endif
          END IF
          run_grid(ng)=.FALSE.
        END DO
!
!  Advance the time counter by the smallest dt.
!
        rtime=rtime+dt(Ngrids)
!
!  Determine what grids can be time stepped. This is determined
!  by comparing dt(each grid) to global time rtime.
!
        DO ng=1,Ngrids
          cff=rtime-rtime_start
          IF (MOD(cff,dt(ng)).eq.0.0_r8) THEN
            run_grid(ng)=.TRUE.
          END IF
        END DO
      END DO
      IF (exit_flag.ne.NoError) THEN
        IF (Master) THEN
          WRITE (stdout,'(/,a,i3,/)') Rerror(exit_flag), exit_flag
        END IF
        RETURN
      END IF

#elif defined COMPOSED_GRID

      NL_LOOP : DO my_iic=Tstr(1),Tend(1)

        NEST_LOOP : DO ng=1,Ngrids
          iic(ng)=my_iic
        END DO NEST_LOOP

# ifdef SOLVE3D
            CALL main3d
# else
            CALL main2d
# endif

          IF (exit_flag.ne.NoError) THEN
            IF (Master) THEN
              WRITE (stdout,'(/,a,i3,/)') Rerror(exit_flag), exit_flag
            END IF
            RETURN
          END IF


      END DO NL_LOOP

#else

      NEST_LOOP : DO ng=1,Ngrids

        NL_LOOP : DO my_iic=Tstr(ng),Tend(ng)

          iic(ng)=my_iic
# ifdef SOLVE3D
!         CALL main3d (ng)
          CALL main3d
# else
!         CALL main2d (ng)
          CALL main2d
# endif
#ifdef INWAVE_COUPLING
          CALL inwave_run (ng)
#endif
          IF (exit_flag.ne.NoError) THEN
            IF (Master) THEN
              WRITE (stdout,'(/,a,i3,/)') Rerror(exit_flag), exit_flag
            END IF
            RETURN
          END IF

        END DO NL_LOOP

      END DO NEST_LOOP
#endif

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear model execution.        !
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
      integer :: ng, thread

#ifdef VERIFICATION
!
!-----------------------------------------------------------------------
!  Compute and report model-observation comparison statistics.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL stats_modobs (ng)
      END DO
#endif
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
!
      DO ng=1,Ngrids
        IF (LwrtRST(ng).and.(exit_flag.eq.1)) THEN
          IF (Master) WRITE (stdout,10)
 10       FORMAT (/,' Blowing-up: Saving latest model state into ',     & 
     &              ' RESTART file',/)
          IF (LcycleRST(ng).and.(NrecRST(ng).ge.2)) THEN
            tRSTindx(ng)=2
            LcycleRST(ng)=.FALSE.
          END IF
          blowup=exit_flag
          exit_flag=NoError
          CALL wrt_rst (ng)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks.  Close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,' Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL DO
      END DO
!
!  Close IO files.
!
      CALL close_io

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
