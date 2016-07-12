      MODULE ocean_control_mod
!
!svn $Id: picard_ocean.h 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Picard Iterations Driver:                                 !
!                                                                      !
!  This driver is used to perform the Picard iterations test for the   !
!  representers tangent linear model used in IOMs weak constraint 4D   !
!  variational data assimilation (W4DVAR).  Recall that all  tangent   !
!  linear variables are in term of the full fields and the model can   !
!  expressed symbolically as:                                          !
!                                                                      !
!         d(S')/d(t) = N(So) + A(S' - So)                              !
!                                                                      !
!  where S' is the tangent linear state and So is the "basic state".   !
!  The "basic state" here is the solution of previous tangent linear   !
!  model iteration.                                                    !
!                                                                      !
!  This driver uses ESMF conventions for the  initialization,  time-   !
!  stepping,  and  finalization  of the representer  tangent  linear   !
!  model via:                                                          !
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
      USE mod_iounits
      USE mod_scalars

#ifdef MCT_LIB
!
# ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_coupling
# endif
# ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_coupling
# endif
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: mpiCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: chunk_size, ng, thread
#ifdef _OPENMP
      integer :: my_threadnum
#endif

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
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
        IF (exit_flag.ne.NoError) RETURN
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
            CALL wclock_on (ng, iNLM, 0)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Allocate and initialize modules variables.
!
!$OMP PARALLEL
        CALL mod_arrays (allocate_vars)
!$OMP END PARALLEL

      END IF

#if defined MCT_LIB && (defined AIR_OCEAN || defined WAVES_OCEAN)
!
!-----------------------------------------------------------------------
!  Initialize coupling streams between model(s).
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
# ifdef AIR_OCEAN
        CALL initialize_ocn2atm_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
        CALL initialize_ocn2wav_coupling (ng, MyRank)
# endif
      END DO
#endif

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS representer tangent linear model  !
!  for the specified time interval (seconds), RunInterval.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
!  Imported variable declarations
!
      real(r8), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      integer :: ng
!
!-----------------------------------------------------------------------
!  Run Picard iteratons.
!-----------------------------------------------------------------------
!
!  Use ensemble parameters for Picard itereations.
!
      ITER_LOOP : DO Nrun=ERstr,ERend
!
!  Cycle history and forward file names in such a way that the history
!  from the previous iteration becomes the basic state for the next.
!
        DO ng=1,Ngrids
          WRITE (TLM(ng)%name,10) TRIM(TLM(ng)%base), Nrun
          WRITE (FWD(ng)%name,10) TRIM(TLM(ng)%base), Nrun-1

          IF (Master) THEN
            WRITE (stdout,20) 'ROMS/TOMS Picard Iteration: ', Nrun, ng, &
     &                        TRIM(TLM(ng)%name),                       &
     &                        TRIM(FWD(ng)%name)
          END IF
        END DO
!
!  Activate defining history an restart files on each iteration. The
!  restart file is used to the store the solution of each iteration.
!
        DO ng=1,Ngrids
          iic(ng)=0
          LdefTLM(ng)=.TRUE.
          LwrtTLM(ng)=.TRUE.
          LdefRST(ng)=.FALSE.
        END DO
!
!  Initialize representer tangent linear model.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          CALL rp_initial (ng)
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Time-step representers tangent linear model
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,30) 'RP', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#ifdef SOLVE3D
        CALL rp_main3d (RunInterval)
#else
        CALL rp_main2d (RunInterval)
#endif
!$OMP END PARALLEL
        IF (exit_flag.ne.NoError) RETURN
!
!  Close IO and re-initialize NetCDF switches.
!
        DO ng=1,Ngrids
          IF (TLM(ng)%ncid.ne.-1) THEN
            CALL netcdf_close (ng, iRPM, TLM(ng)%ncid)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF (FWD(ng)%ncid.ne.-1) THEN
            CALL netcdf_close (ng, iRPM, FWD(ng)%ncid)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

      END DO ITER_LOOP
!
 10   FORMAT (a,'_',i2.2,'.nc')
 20   FORMAT (/,a,i3,2x,'(Grid: ',i2.2,')',/,                           &
     &        /,5x,'  History file: ',a,                                &
     &        /,5x,'  Forward file: ',a,/)
 30   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')',/)

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
      integer :: Fcount, ng, thread
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
            IF (Master) WRITE (stdout,10)
 10         FORMAT (/,' Blowing-up: Saving latest model state into ',   &
     &                ' RESTART file',/)
            Fcount=RST(ng)%Fcount
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
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL
      END DO

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
