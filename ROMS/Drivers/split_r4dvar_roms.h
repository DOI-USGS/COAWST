      MODULE roms_kernel_mod
!
!git $Id$
!svn $Id: split_r4dvar_roms.h 1151 2023-02-09 03:08:53Z arango $
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS Strong/Weak Constraint Split 4-Dimensional Variational Data    !
!            Assimilation Driver: Indirect Representer Approach        !
!                                 (R4D-Var).                           !
!                                                                      !
!  This driver is used for the dual formulation (observation space),   !
!  strong or weak constraint 4D-Var where errors may be considered     !
!  in both model and observations.                                     !
!                                                                      !
!  The R4D-Var algorithm is split into multiple executables to         !
!  facilitate various configurations:                                  !
!                                                                      !
!    (1) Executable A computes ROMS nonlinear trajectory used to       !
!        linearize the tangent linear and adjoint models used in       !
!        the iterations of the inner loop for the minimization of      !
!        the cost function. It allows the nonlinear trajectory to      !
!        be part of a coupling system and or include nested grids.     !
!        It calls either the RBL4D-Var "background" or "analysis"      !
!        routines.                                                     !
!                                                                      !
!    (2) Executable B calls either R4D-Var "increment" or              !
!        "posterior_error". The R4D-Var increment is obtained by       !
!        minimizing the cost function over Ninner loops. It is         !
!        possible to use a coarser grid resolution in the inner        !
!        loop.  If so, the finer background trajectory needs to        !
!        be interpolated into the coarser grid. Then, at the end       !
!        of inner loops, the coarse grid increment needs to be         !
!        interpolated to the finer grid.  The increment phase          !
!        may be run at a lower precision.                              !
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
      USE mod_param
      USE mod_parallel
      USE mod_arrays
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
#if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
#endif
      USE mod_scalars
      USE mod_stepping
!
      USE r4dvar_mod
!
      USE close_io_mod,      ONLY : close_inp, close_out
      USE def_dai_mod,       ONLY : def_dai
      USE get_state_mod,     ONLY : get_state
      USE inp_par_mod,       ONLY : inp_par
#ifdef MCT_LIB
# ifdef ATM_COUPLING
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_coupling
# endif
# ifdef WAV_COUPLING
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_coupling
# endif
#endif
      USE stats_modobs_mod,  ONLY : stats_modobs
      USE stdinp_mod,        ONLY : getpar_s
      USE strings_mod,       ONLY : FoundError, uppercase
      USE wrt_dai_mod,       ONLY : wrt_dai
      USE wrt_rst_mod,       ONLY : wrt_rst
!
      implicit none
!
      PUBLIC :: ROMS_initialize
      PUBLIC :: ROMS_run
      PUBLIC :: ROMS_finalize
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
!  Get 4D-Var phase from APARNAM input script file.
!
        CALL getpar_s (MyRank, aparnam, 'APARNAM')
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL getpar_s (MyRank, Phase4DVAR, 'Phase4DVAR', aparnam)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Determine ROMS standard output append switch. It is only relevant if
!  "ROMS_STDINP" is activated. The standard output is created in the
!  "background" phase and open to append in the other phases. Set
!  switch so the "stiffness" routine is only called in the "background"
!  phase.
!
        IF (INDEX(TRIM(uppercase(Phase4DVAR)),'BACKG').ne.0) THEN
          Lappend=.FALSE.
          Lstiffness=.TRUE.
        ELSE
          Lappend=.TRUE.
          Lstiffness=.FALSE.
        END IF
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
        CALL ROMS_allocate_arrays (allocate_vars)
        CALL ROMS_initialize_arrays
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
!  proccess background and mode prior error covariance standard
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
!  Run Split R4D-Var Data Assimilation algorithm.
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
      Ldone=.FALSE.         ! 4D-Var cycle finish switch
!
!  Select R4D-Var phase to execute.
!
      SELECT CASE (uppercase(Phase4DVAR(1:6)))

!
!  Compute nonlinear background state trajectory, Xb(t)|n-1. Interpolate
!  the background at the observation locations, and compute the quality
!  control accept/reject flag, ObsScale. The background state is used
!  to linearize the tangent linear and adjoint models during the
!  minimization.
!
        CASE ('BACKGR')

          my_outer=0
          outer=0
          inner=0

          CALL background (outer, RunInterval)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Compute 4D-Var data assimilation increment, dXa, by iterating over
!  the inner loops, and minimizing the cost function.
!
        CASE ('INCREM')

          my_outer=OuterLoop
          outer=OuterLoop
          inner=0

          CALL increment (my_outer, RunInterval)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Compute 4D-Var data assimilation analysis, Xa = Xb + dXa.  Set
!  nonlinear model initial conditions for next outer loop.
!
        CASE ('ANALYS')

          my_outer=OuterLoop
          outer=OuterLoop
          inner=Ninner

          CALL analysis (my_outer, RunInterval)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN


#if defined POSTERIOR_ERROR_I || \
    defined POSTERIOR_ERROR_F || \
    defined POSTERIOR_EOFS
!
!  Compute full (diagonal) posterior analysis error covariance matrix.
!  (NOTE: Currently, this code only works for a single outer-loop).
!
        CASE ('POST_E')

          CALL posterior_error (RunInterval)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!  Issue an error if incorrect 4D-Var phase.
!
        CASE DEFAULT

          IF (Master) THEN
            WRITE (stdout,10) TRIM(Phase4DVAR)
 10         FORMAT (' ROMS_run - illegal 4D-Var phase: ''',a,'''')
          END IF
          exit_flag=5
          RETURN

      END SELECT
!
!  Set finish R4D-Var cycle switch.
!
      IF ((OuterLoop.eq.Nouter).and.                                    &
#if defined POSTERIOR_ERROR_I || \
    defined POSTERIOR_ERROR_F || \
    defined POSTERIOR_EOFS
     &    (INDEX(TRIM(uppercase(Phase4DVAR)),'POST_E').ne.0)) THEN
        Ldone=.TRUE.
#else
     &    (INDEX(TRIM(uppercase(Phase4DVAR)),'ANALYS').ne.0)) THEN
        Ldone=.TRUE.
#endif
      END IF
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
!  Local variable declarations.
!
      integer :: Fcount, InpRec, Nfiles, Tindex
      integer :: ifile, lstr, ng, thread
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
      IF (Ldone.and.(exit_flag.eq.NoError)) THEN
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
          SELECT CASE (TLM(ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_get_dim (ng, iRPM, TLM(ng)%name,              &
     &                             DimName = 'ocean_time',              &
     &                             DimSize = InpRec)

#if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_get_dim (ng, iRPM, TLM(ng)%name,          &
     &                                 DimName = 'ocean_time',          &
     &                                 DimSize = InpRec)
#endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          Tindex=1
          CALL get_state (ng, iRPM, 1, TLM(ng), InpRec, Tindex)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          KOUT=Tindex
          NOUT=Tindex
#ifdef DISTRIBUTE
          CALL wrt_dai (ng, MyRank)
#else
          CALL wrt_dai (ng, -1)
#endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Compute and report model-observation comparison statistics.
!-----------------------------------------------------------------------
!
      IF (Ldone.or.(exit_flag.eq.1)) THEN
        DO ng=1,Ngrids
#ifdef DISTRIBUTE
          CALL stats_modobs (ng, MyRank)
#else
          CALL stats_modobs (ng, -1)
#endif
        END DO
      END IF
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
#ifdef DISTRIBUTE
            CALL wrt_rst (ng, MyRank)
#else
            CALL wrt_rst (ng, -1)
#endif
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

      END MODULE roms_kernel_mod
