      MODULE roms_kernel_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  ROMS 4DVAR Background-error Correlation Model                       !
!                                                                      !
!  This driver is used to build and test the 4DVAR background-error    !
!  correlation model using a generalized difussion operator:           !
!                                                                      !
!            B = S C S                                                 !
!                                                                      !
!            C = C^(1/2) C^(T/2)                                       !
!                                                                      !
!      C^(1/2) = G L^(1/2) W^(-1/2)                                    !
!                                                                      !
!      C^(T/2) = W^(T/2) L^(T/2) G                                     !
!                                                                      !
!  where                                                               !
!                                                                      !
!         B : background-error covariance                              !
!         C : background-error correlation                             !
!         G : normalization coefficient matrix                         !
!         L : self-adjoint diffusion filter                            !
!         S : background-error standard deviation                      !
!         W : Grid cell area or volume diagonal matrix                 !
!                                                                      !
!  The routines in this driver control the initialization,  time-      !
!  stepping, and finalization of  ROMS  model following ESMF           !
!  conventions:                                                        !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!  Reference                                                           !
!                                                                      !
!     Weaver, A. and P. Courtier, 2001: Correlation modelling on       !
!       the sphere using a generalized diffusion equation, Q. J.       !
!       R. Meteorol. Soc., 127, 1815-1845.                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_arrays
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
!
      USE analytical_mod,       ONLY : ana_perturb
      USE close_io_mod,         ONLY : close_inp,                       &
     &                                 close_out
      USE  convolve_mod,        ONLY : convolve
      USE def_norm_mod,         ONLY : def_norm
      USE get_state_mod,        ONLY : get_state
      USE inp_par_mod,          ONLY : inp_par
#ifdef MULTI_SCALE_B
      USE multiscale_eigen_mod, ONLY : multiscale_eigen,                &
     &                                 multiscale_eigen_read,           &
     &                                 multiscale_eigen_write
#endif
      USE normalization_mod,    ONLY : normalization
#ifdef MULTI_SCALE_B
      USE roms_multiscale_mod,  ONLY : MSB
# ifdef NONUNIFORM_SCALES
      USE roms_multiscale_mod,  ONLY : multiscale_get_scales
# endif
#endif
      USE stdout_mod,           ONLY : Set_StdOutUnit,                  &
     &                                 stdout_unit
      USE strings_mod,          ONLY : FoundError
      USE tl_def_his_mod,       ONLY : tl_def_his
      USE tl_wrt_his_mod,       ONLY : tl_wrt_his
      USE wrt_rst_mod,          ONLY : wrt_rst
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
!  This routine allocates and initializes ROMS state variables         !
!  and internal and external parameters.                               !
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
      integer :: NRMrec, STDrec, Tindex, ifac
      integer :: chunk_size, ng, thread, tile
#ifdef _OPENMP
      integer :: my_threadnum
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_initialize"

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (mpi) world communicator.
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
!  Set the ROMS standard output unit to write verbose execution info.
!  Notice that the default standard out unit in Fortran is 6.
!
!  In some applications like coupling or disjointed mpi-communications,
!  it is advantageous to write standard output to a specific filename
!  instead of the default Fortran standard output unit 6. If that is
!  the case, it opens such formatted file for writing.
!
        IF (Set_StdOutUnit) THEN
          stdout=stdout_unit(Master)
          Set_StdOutUnit=.FALSE.
        END IF
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
!
!-----------------------------------------------------------------------
!  Initialize metrics over all nested grids, if applicable.
!-----------------------------------------------------------------------
!
      CALL initial
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Adjust "time" variable since we are not time-stepping.
!
      DO ng=1,Ngrids
        time(ng)=time(ng)+dt(ng)
      END DO
!
!  Initialize run or ensemble counter.
!
      Nrun=1
!
!-----------------------------------------------------------------------
!  Read in standard deviation factors for error covariance.
!-----------------------------------------------------------------------
!
#ifdef WEAK_CONSTRAINT
      NSA=2               ! include weak constraint error hypothesis
#else
      NSA=1               ! only strong constraint error hypothesis
#endif
!
!  Initial conditions standard deviation. They are loaded in Tindex=1
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 10, 10, STD(1,ng), STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  Model error standard deviation. They are loaded in Tindex=2
!  of the e_var(...,Tindex) state variables.
!
      IF (NSA.eq.2) THEN
        STDrec=1
        Tindex=2
        DO ng=1,Ngrids
          CALL get_state (ng, 11, 11, STD(2,ng), STDrec, Tindex)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF

#ifdef ADJUST_BOUNDARY
!
!  Open boundary conditions standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 12, 12, STD(3,ng), STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Surface forcing standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 13, 13, STD(4,ng), STDrec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Compute or read in error covariance normalization factors.
!-----------------------------------------------------------------------
!
      NESTED_LOOP : DO ng=1,Ngrids

#if defined MULTI_SCALE_B && defined NONUNIFORM_SCALES
!
!  Read in horizontal, spatially-varying correlation length scales.
!
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL multiscale_get_scales (MSB(ng), ng, tile, iTLM)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
#endif
!
!  Process normalization coefficients.
!
        GET_NORMALIZATION : IF (ANY(LwrtNRM(:,ng))) THEN
!
!  If computing, define output normalization NetCDF file(s).
!
          IF (LdefNRM(1,ng).or.LwrtNRM(1,ng)) THEN
            CALL def_norm (ng, iNLM, 1)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF ((LdefNRM(2,ng).or.LwrtNRM(2,ng)).and.(NSA.eq.2)) THEN
            CALL def_norm (ng, iNLM, 2)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

#ifdef ADJUST_BOUNDARY
!
          IF (LdefNRM(3,ng).or.LwrtNRM(3,ng)) THEN
            CALL def_norm (ng, iNLM, 3)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
          IF (LdefNRM(4,ng).or.LwrtNRM(4,ng)) THEN
            CALL def_norm (ng, iNLM, 4)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
#endif

#ifdef MULTI_SCALE_B
!
!  Compute the extrema eigenvalues of the K-Laplacian operator required
!  by the Chebyshev Iterations (CI) solver, which is applied to implicit
!  diffusion operators in the modeling of the multiscale background-
!  error covariance for variables in the control vector.
!
!  The eigenvalue spectrum of the K-Laplacian operator remains invariant
!  for a fixed application grid and a given value of K. Consequently,
!  estimates can be precomputed via the Lanczos formulation of the
!  Conjugate Gradient (CG) method, initialized with random vectors.
!
          ifac=1
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL multiscale_eigen (ng, tile, Lnew(ng), ifac)
          END DO
!
!  Write out extrema eigenvalues into output normalizations NetCDF
!  file(s).
!
          CALL multiscale_eigen_write (MSB(ng), ng, iTLM)
#endif
!
!  Compute normalization factors. The ifac=2 indicates the squared-root
!  operator, so the spatial convolution is applied for only half of the
!  pseudo-diffusion steps.
!
          ifac=2
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL normalization (ng, tile, ifac)
          END DO
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          LdefNRM(1:4,ng)=.FALSE.
          LwrtNRM(1:4,ng)=.FALSE.
!
!  Otherwise, read in normalization coefficients.
!
        ELSE
        
          NRMrec=1
          CALL get_state (ng, 14, 14, NRM(1,ng), NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          IF (NSA.eq.2) THEN
            CALL get_state (ng, 15, 15, NRM(2,ng), NRMrec, 2)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

#ifdef ADJUST_BOUNDARY
!
          CALL get_state (ng, 16, 16, NRM(3,ng), NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
          CALL get_state (ng, 17, 17, NRM(4,ng), NRMrec, 1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef MULTI_SCALE_B
!
!  Read in extrema Ritz eigenvalues frot the normalization file, which
!  are required by the Implicit Chebyshev Iterations (CI) solver.
!
          CALL multiscale_eigen_read (MSB(ng), ng, iTLM)
#endif
        END IF GET_NORMALIZATION 
      END DO NESTED_LOOP
!
      RETURN
      END SUBROUTINE ROMS_initialize
!
      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes background-error correlations.                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      integer :: i, ng, tile
      integer :: Lini = 1
!
      character (len=8) :: driver = 'rbl4dvar'
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_run"
!
!-----------------------------------------------------------------------
!  Test correlation model: Dirac Delta Functions.
!-----------------------------------------------------------------------
!
!  Initialize adjoint model state with a delta function at specified
!  point. Use USER parameters from standard input to perturb solution
!  in routine "ana_perturb". Then, convolve solution with the adjoint
!  diffusion operator.
!
      ADmodel=.FALSE.
      TLmodel=.TRUE.
!
      DO ng=1,Ngrids
        Lold(ng)=1
        Lnew(ng)=2
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_perturb (ng, tile, iTLM)
        END DO
      END DO
!
!  Apply background-error covariance convolutions.
!
      CALL convolve (driver, Lini, Lold, Lnew)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out background error correlation in adjoint history NetCDF
!  file.
!
      DO ng=1,Ngrids
        kstp(ng)=Lnew(ng)
#ifdef SOLVE3D
        nrhs(ng)=Lnew(ng)
#endif
        Lfout(ng)=Lnew(ng)
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        LwrtState2d(ng)=.TRUE.
        CALL tl_def_his (ng, LdefTLM(ng))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Ladjusted(ng)=.TRUE.
#endif
#ifdef DISTRIBUTE
        CALL tl_wrt_his (ng, MyRank)
#else
        CALL tl_wrt_his (ng, -1)
#endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Ladjusted(ng)=.FALSE.
#endif
      END DO
!
      RETURN
      END SUBROUTINE ROMS_run
!
      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS driver execution.                      !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: Fcount, ng, thread
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_finalize"
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
            IF (Master) WRITE (stdout,10)
 10         FORMAT (/,' Blowing-up: Saving latest model state into ',   &
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
