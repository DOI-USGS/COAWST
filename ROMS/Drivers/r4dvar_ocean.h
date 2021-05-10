      MODULE ocean_control_mod
!
!svn $Id: r4dvar_ocean.h 1054 2021-03-06 19:47:12Z arango $
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS Strong/Weak Constraint 4-Dimensional Variational Data          !
!            Assimilation Driver: Indirect Representer Approach        !
!                                 (R4D-Var).                           !
!                                                                      !
!  This driver is used for the dual formulation (observation space),   !
!  strong or weak constraint 4D-Var where errors may be considered     !
!  in both model and observations.                                     !
!                                                                      !
!  The routines in this driver control the initialization,  time-      !
!  stepping, and finalization of ROMS model following ESMF/NUOPC       !
!  conventions:                                                        !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Moore, A.M., H.G. Arango, G. Broquet, B.S. Powell, A.T. Weaver,   !
!      and J. Zavala-Garay, 2011: The Regional Ocean Modeling System   !
!      (ROMS)  4-dimensional variational data assimilations systems,   !
!      Part I - System overview and formulation, Prog. Oceanogr., 91,  !
!      34-49, doi:10.1016/j.pocean.2011.05.004.                        !
!                                                                      !
!    Moore, A.M., H.G. Arango, G. Broquet, C. Edward, M. Veneziani,    !
!      B. Powell, D. Foley, J.D. Doyle, D. Costa, and P. Robinson,     !
!      2011: The Regional Ocean Modeling System (ROMS) 4-dimensional   !
!      variational data assimilations systems, Part II - Performance   !
!      and application to the California Current System, Prog.         !
!      Oceanogr., 91, 50-73, doi:10.1016/j.pocean.2011.05.003.         !
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
!  This routine allocates and initializes ROMS state variables and     !
!  internal parameters. It reads standard input parameters.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
      USE inp_par_mod,       ONLY : inp_par
#ifdef MCT_LIB
# ifdef ATM_COUPLING
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_coupling
# endif
# ifdef WAV_COUPLING
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_coupling
# endif
#endif
      USE r4dvar_mod,        ONLY : prior_error
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
!  Initialize counters. The 'Nrun' counter will be recomputed in the
!  RBL4D-Var phases to process the obervation operator correctly.
!
        Nrun=1                ! run counter
        ERstr=1               ! ensemble start counter
        ERend=Nouter          ! ensemble end counter
!
!  Set domain decomposition tile partition range.  This range is
!  computed only once since the "first_tile" and "last_tile" values
!  are private for each parallel thread/node.
!
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
          DO thread=THREAD_RANGE
            CALL wclock_on (ng, iNLM, 0, __LINE__, MyFile)
          END DO
        END DO
!
!  Allocate and initialize modules variables.
!
        CALL mod_arrays (allocate_vars)
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar

      END IF

#if defined MCT_LIB && (defined ATM_COUPLING || defined WAV_COUPLING)
!
!-----------------------------------------------------------------------
!  Initialize coupling streams between model(s).
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
# ifdef ATM_COUPLING
        CALL initialize_ocn2atm_coupling (ng, MyRank)
# endif
# ifdef WAV_COUPLING
        CALL initialize_ocn2wav_coupling (ng, MyRank)
# endif
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Set application grid, metrics, and associated variables. Then,
!  Proccess background and model prior error covariance standard
!  deviations and normalization coefficients.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL prior_error (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        SetGridConfig(ng)=.FALSE.
      END DO
!
      RETURN
      END SUBROUTINE ROMS_initialize
!
      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This subroutine runs the Strong or Weak constraint, Indirect        !
!  Representers 4D-Var data assimilation (R4D-Var) algorithm. It       !
!  time-steps ROMS nonlinear, representer, tangent linear, and         !
!  adjoint kernels.                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     RunInterval     Execution time stepping window (seconds)         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
      USE r4dvar_mod,  ONLY : background, increment, analysis

#if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F || \
    defined POSTERIOR_EOFS
      USE r4dvar_mod,  ONLY : posterior_error
#endif
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: my_outer, ng
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_run"
!
!=======================================================================
!  Run R4D-Var Data Assimilation algorithm.
!=======================================================================
!
!  Initialize several global parameters.
!
      DO ng=1,Ngrids
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Lfinp(ng)=1         ! forcing index for input
        Lfout(ng)=1         ! forcing index for output history files
#endif
#ifdef ADJUST_BOUNDARY
        Lbinp(ng)=1         ! boundary index for input
        Lbout(ng)=1         ! boundary index for output history files
#endif
        Lold(ng)=1          ! old minimization time index
        Lnew(ng)=2          ! new minimization time index
      END DO
!
      outer=0               ! outer-loop counter
      inner=0               ! inner-loop counter
!
!  Compute nonlinear background state trajectory, Xb(t)|n-1. Interpolate
!  the background at the observation locations, and compute the quality
!  control accept/reject flag, ObsScale. The background state is used
!  to linearize the tangent linear and adjoint models during the
!  minimization.
!
      CALL background (outer, RunInterval)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Start outer loop iterations.
!
      OUTER_LOOP : DO my_outer=1,Nouter
        outer=my_outer
        inner=0
!
!  Compute 4D-Var data assimilation increment, dXa, by iterating over
!  the inner loops, and minimizing the cost function.
!
        CALL increment (my_outer, RunInterval)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Compute 4D-Var data assimilation analysis, Xa = Xb + dXa.  Set
!  nonlinear model initial conditions for next outer loop.
!
        CALL analysis (my_outer, RunInterval)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END DO OUTER_LOOP

#if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F || \
    defined POSTERIOR_EOFS
!
!  Compute full (diagonal) posterior analysis error covariance matrix.
!  (NOTE: Currently, this code only works for a single outer-loop).
!
      CALL posterior_error (RunInterval)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
      RETURN
      END SUBROUTINE ROMS_run
!
      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS R4D-Var execution.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      USE strings_mod, ONLY : FoundError
!
!  Local variable declarations.
!
      integer :: Fcount, InpRec, Nfiles, Tindex
      integer :: ifile, lstr, ng, tile, thread
!
      character (len=10) :: suffix

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_finalize"
!
!-----------------------------------------------------------------------
!  Write out 4D-Var analysis fields that can be used as the initial
!  conditions for the next data assimilation cycle. Here, use the
!  last record of the RPM for the final outer loop.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
      tile=MyRank
#else
      tile=-1
#endif
!
      IF (exit_flag.eq.NoError) THEN
        DO ng=1,Ngrids
          LdefDAI(ng)=.TRUE.
          CALL def_dai (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          WRITE (TLM(ng)%name,10) TRIM(FWD(ng)%head), Nouter
 10       FORMAT (a,'_outer',i0,'.nc')
          lstr=LEN_TRIM(TLM(ng)%name)
          TLM(ng)%base=TLM(ng)%name(1:lstr-3)
          IF (TLM(ng)%Nfiles.gt.1) THEN
            Nfiles=TLM(ng)%Nfiles
            DO ifile=1,Nfiles
              WRITE (suffix,"('_',i4.4,'.nc')") ifile
              TLM(ng)%files(ifile)=TRIM(TLM(ng)%base)//TRIM(suffix)
            END DO
            TLM(ng)%name=TRIM(TLM(ng)%files(Nfiles))
          ELSE
            TLM(ng)%files(1)=TRIM(TLM(ng)%name)
          END IF
!
          CALL netcdf_get_dim (ng, iRPM, TLM(ng)%name,                  &
     &                         DimName = 'ocean_time',                  &
     &                         DimSize = InpRec)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          Tindex=1
          CALL get_state (ng, iRPM, 1, TLM(ng)%name, InpRec, Tindex)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          KOUT=Tindex
          NOUT=Tindex
          CALL wrt_dai (ng, tile)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Compute and report model-observation comparison statistics.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL stats_modobs (ng)
      END DO
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into record 3.
!
      IF (exit_flag.eq.1) THEN
        DO ng=1,Ngrids
          IF (LwrtRST(ng)) THEN
            IF (Master) WRITE (stdout,20)
 20         FORMAT (/,' Blowing-up: Saving latest model state into ',   &
     &                ' RESTART file',/)
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
!  Stop model and time profiling clocks, report memory requirements,
!  and close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,30)
 30     FORMAT (/,'Elapsed wall CPU time for each process (seconds):',/)
      END IF
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 0, __LINE__, MyFile)
        END DO
      END DO
!
!  Report dynamic memory and automatic memory requirements.
!
      CALL memory
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
