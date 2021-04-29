      MODULE ocean_control_mod
!
!svn $Id: nl_ocean.h 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Nonlinear Model Driver:                                   !
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
!
      PRIVATE
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize
!
      CONTAINS
!
      SUBROUTINE ROMS_initialize (first, mpiCOMM)
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
      USE inp_par_mod,       ONLY : inp_par
#ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_coupling
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_routers
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_coupling
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_routers
#endif
#ifdef HYDRO_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2hyd_coupling
      USE ocean_coupler_mod, ONLY : initialize_ocn2hyd_routers
#endif
#ifdef INWAVE_MODEL
      USE driver_inwave_mod, ONLY : inwave_init
#endif
      USE strings_mod,       ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first
!
      integer, intent(in), optional :: mpiCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.
!
#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: chunk_size, ng, thread
#ifdef _OPENMP
      integer :: my_threadnum
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_initialize"

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (mpi) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(mpiCOMM)) THEN
        OCN_COMM_WORLD=mpiCOMM
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
!  Initialize parallel control switches. These scalars switches are
!  independent from standard input parameters.
!
        CALL initialize_parallel
!
!  Read in model tunable parameters from standard input. Allocate and
!  initialize variables in several modules after the number of nested
!  grids and dimension parameters are known.
!
        CALL inp_par (iNLM)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set domain decomposition tile partition range.  This range is
!  computed only once since the "first_tile" and "last_tile" values
!  are private for each parallel thread/node.
!
!$OMP PARALLEL
#if defined _OPENMP
      MyThread=my_threadnum()
#elif defined DISTRIBUTE
      MyThread=MyRank
#else
      MyThread=0
#endif
      DO ng=1,Ngrids
        chunk_size=(NtileX(ng)*NtileE(ng)+numthreads-1)/numthreads
        first_tile(ng)=MyThread*chunk_size
        last_tile (ng)=first_tile(ng)+chunk_size-1
      END DO
!$OMP END PARALLEL
!
!  Initialize internal wall clocks. Notice that the timings does not
!  includes processing standard input because several parameters are
!  needed to allocate clock variables.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (/,' Process Information:',/)
        END IF
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO thread=THREAD_RANGE
            CALL wclock_on (ng, iNLM, 0, __LINE__, MyFile)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Allocate and initialize all model state arrays.
!
!$OMP PARALLEL
        CALL mod_arrays (allocate_vars)
!$OMP END PARALLEL

#ifdef VERIFICATION
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar
#endif
      END IF

#if defined MCT_LIB && (defined AIR_OCEAN || defined WAVES_OCEAN || \
    defined HYDRO_OCEAN)
!
!-----------------------------------------------------------------------
!  Initialize coupling streams between model(s).
!-----------------------------------------------------------------------
!
# ifdef WAVES_OCEAN
        DO ng=1,Ngrids
          CALL initialize_ocn2wav_coupling (ng, MyRank)
        END DO
        CALL initialize_ocn2wav_routers (MyRank)
# endif
# ifdef AIR_OCEAN
        DO ng=1,Ngrids
          CALL initialize_ocn2atm_coupling (ng, MyRank)
        END DO
        CALL initialize_ocn2atm_routers (MyRank)
# endif
# ifdef HYDRO_OCEAN
        DO ng=1,Ngrids
          CALL initialize_ocn2hyd_coupling (ng, MyRank)
        END DO
        CALL initialize_ocn2hyd_routers (MyRank)
# endif
#endif
#ifdef INWAVE_MODEL
        DO ng=1,Ngrids
!         CALL inwave_init (ng, MyRank)
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Initialize nonlinear model state variables over all nested grids,
!  if applicable.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL
      CALL initial
!$OMP END PARALLEL
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize run or ensemble counter.
!
      Nrun=1

#ifdef VERIFICATION
!
!  Create NetCDF file for model solution at observation locations.
!
      IF (Nrun.eq.1) THEN
        DO ng=1,Ngrids
          LdefMOD(ng)=.TRUE.
          wrtNLmod(ng)=.TRUE.
          wrtObsScale(ng)=.TRUE.
          CALL def_mod (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
#endif
#ifdef ENKF_RESTART
!
!  Create Ensenble Kalman Filter (EnKF) reastart NetCDF file.
!
      IF (Nrun.eq.1) THEN
        DO ng=1,Ngrids
          LdefDAI(ng)=.TRUE.
          CALL def_dai (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
#endif
!
      RETURN
      END SUBROUTINE ROMS_initialize
!
      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine runs ROMS/TOMS nonlinear model for the specified time  !
!  interval (seconds), RunInterval.  It RunInterval=0, ROMS advances   !
!  one single time-step.                                               !
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
      USE strings_mod, ONLY : FoundError
#ifdef INWAVE_MODEL
      USE driver_inwave_mod, ONLY : inwave_run
#endif
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      integer :: ng
#if defined MODEL_COUPLING && !defined MCT_LIB
      integer :: NstrStep, NendStep, extra
!
      real(dp) :: ENDtime, NEXTtime
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_run"
!
!-----------------------------------------------------------------------
!  Time-step nonlinear model over nested grids, if applicable.
#if defined MODEL_COUPLING && !defined MCT_LIB
!  Since the ROMS kernel has a delayed output and line diagnostics by
!  one timestep, subtact an extra value to the report of starting and
!  ending timestep for clarity. Usually, the model coupling interval
!  is of the same size as ROMS timestep.
#endif
!-----------------------------------------------------------------------
!
      MyRunInterval=RunInterval
      IF (Master) WRITE (stdout,'(1x)')
      DO ng=1,Ngrids
#if defined MODEL_COUPLING && !defined MCT_LIB
        NEXTtime=time(ng)+RunInterval
        ENDtime=INItime(ng)+(ntimes(ng)-1)*dt(ng)
        IF ((NEXTtime.eq.ENDtime).and.(ng.eq.1)) THEN
          extra=0                                   ! last time interval
        ELSE
          extra=1
        END IF
        step_counter(ng)=0
        NstrStep=iic(ng)
        NendStep=NstrStep+INT((MyRunInterval)/dt(ng))-extra
        IF (Master) WRITE (stdout,10) 'NL', ng, NstrStep, NendStep
#else
        IF (Master) WRITE (stdout,10) 'NL', ng, ntstart(ng), ntend(ng)
#endif
      END DO
      IF (Master) WRITE (stdout,'(1x)')
!
!$OMP PARALLEL
#ifdef SOLVE3D
      CALL main3d (MyRunInterval)
#else
      CALL main2d (MyRunInterval)
#endif
#ifdef INWAVE_MODEL
!     CALL inwave_run
#endif
!$OMP END PARALLEL

      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
 10   FORMAT (1x,a,1x,'ROMS/TOMS: started time-stepping:',              &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i12.12,' - ',i12.12,')')
!
      RETURN
      END SUBROUTINE ROMS_run
!
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
      integer :: Fcount, ng, thread
#ifdef ENKF_RESTART
      integer :: tile
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_finalize"

#ifdef ENKF_RESTART
!
!-----------------------------------------------------------------------
!  Write out initial conditions for the next time window of the Ensemble
!  Kalman (EnKF) filter.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      tile=MyRank
# else
      tile=-1
# endif
!
      IF (exit_flag.eq.NoError) THEN
        DO ng=1,Ngrids
          CALL wrt_dai (ng, tile)
        END DO
      END IF
#endif
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
      IF (exit_flag.eq.1) THEN
        DO ng=1,Ngrids
          IF (LwrtRST(ng)) THEN
            IF (Master) WRITE (stdout,10) TRIM(blowup_string)
 10         FORMAT (/,' Blowing-up: Saving latest model state into ',   &
     &                ' RESTART file',/,'     REASON: ',a,/)
            Fcount=RST(ng)%load
            IF (LcycleRST(ng).and.(RST(ng)%Nrec(Fcount).ge.2)) THEN
              RST(ng)%Rindex=2
              LcycleRST(ng)=.FALSE.
            END IF
            blowup=exit_flag
            exit_flag=NoError
            CALL wrt_rst (ng)
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks, report memory requirements, and
!  close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,'Elapsed wall CPU time for each process (seconds):',/)
      END IF
!
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 0, __LINE__, MyFile)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Report dynamic memory and automatic memory requirements.
!
!$OMP PARALLEL
      CALL memory
!$OMP END PARALLEL
!
!  Close IO files.
!
      DO ng=1,Ngrids
        CALL close_inp (ng, iNLM)
      END DO
      CALL close_out
!
      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
