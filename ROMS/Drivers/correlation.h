      MODULE ocean_control_mod
!
!svn $Id: correlation.h 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS 4DVAR Background-error Correlation Model                  !
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
!  stepping, and finalization of  ROMS/TOMS  model following ESMF      !
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
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE normalization_mod, ONLY : normalization
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
      integer :: STDrec, Tindex
      integer :: chunk_size, ng, thread, tile
#ifdef _OPENMP
      integer :: my_threadnum
#endif

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communicator.
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
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar

      END IF
!
!-----------------------------------------------------------------------
!  Initialize metrics over all nested grids, if applicable.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL
      CALL initial
!$OMP END PARALLEL
      IF (exit_flag.ne.NoError) RETURN
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
!  Initial conditions standard deviation. They are loaded in Tindex=1
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        IF (LdefNRM(1,ng).or.LwrtNRM(1,ng)) THEN
          CALL get_state (ng, 6, 6, STD(1,ng)%name, STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
!
!  Model error standard deviation. They are loaded in Tindex=2
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=2
      DO ng=1,Ngrids
        IF ((LdefNRM(2,ng).or.LwrtNRM(2,ng)).and.(NSA.eq.2)) THEN
          CALL get_state (ng, 6, 6, STD(2,ng)%name, STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO

#ifdef ADJUST_BOUNDARY
!
!  Open boundary conditions standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        IF (LdefNRM(3,ng).or.LwrtNRM(3,ng)) THEN
          CALL get_state (ng, 8, 8, STD(3,ng)%name, STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Surface forcing standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        IF (LdefNRM(4,ng).or.LwrtNRM(4,ng)) THEN
          CALL get_state (ng, 9, 9, STD(4,ng)%name, STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Compute or read in error covariance normalization factors.
!-----------------------------------------------------------------------
!
!  If computing, write out factors to NetCDF. This is an expensive
!  computation and needs to be computed once for a particular
!  application grid.
!
      DO ng=1,Ngrids
        IF (ANY(LwrtNRM(:,ng))) THEN
          IF (LdefNRM(1,ng).or.LwrtNRM(1,ng)) THEN
            CALL def_norm (ng, iNLM, 1)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          IF ((LdefNRM(2,ng).or.LwrtNRM(2,ng)).and.(NSA.eq.2)) THEN
            CALL def_norm (ng, iNLM, 2)
          IF (exit_flag.ne.NoError) RETURN
          END IF
#ifdef ADJUST_BOUNDARY
          IF (LdefNRM(3,ng).or.LwrtNRM(3,ng)) THEN
            CALL def_norm (ng, iNLM, 3)
            IF (exit_flag.ne.NoError) RETURN
          END IF
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          IF (LdefNRM(4,ng).or.LwrtNRM(4,ng)) THEN
            CALL def_norm (ng, iNLM, 4)
            IF (exit_flag.ne.NoError) RETURN
          END IF
#endif
          IF (exit_flag.ne.NoError) RETURN
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL normalization (ng, tile, 2)
          END DO
!$OMP END PARALLEL
          LdefNRM(1:4,ng)=.FALSE.
          LwrtNRM(1:4,ng)=.FALSE.
        END IF
      END DO

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes background-error correlations.                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
#ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE analytical_mod, ONLY : ana_perturb
      USE ini_adjust_mod, ONLY : load_ADtoTL
      USE ini_adjust_mod, ONLY : load_TLtoAD
#ifdef BALANCE_OPERATOR
      USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod, ONLY: balance_ref, biconj
#endif
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      logical :: Lweak, add
      integer :: i, ng, tile
#ifdef BALANCE_OPERATOR
      integer :: Lbck = 1
#endif
!
!-----------------------------------------------------------------------
!  Test correlation model.
!-----------------------------------------------------------------------

#ifdef BALANCE_OPERATOR
!
!  Read background state, use initial conditions.
!
      DO ng=1,Ngrids
        CALL get_state (ng, iNLM, 9, INI(ng)%name, Lbck, Lbck)
        IF (exit_flag.ne.NoError) RETURN
      END DO

# ifdef ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
      IF (balance(isFsur)) THEN
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL balance_ref (ng, tile, Lbck)
            CALL biconj (ng, tile, iNLM, Lbck)
          END DO
!$OMP END PARALLEL
          wrtZetaRef(ng)=.TRUE.
        END DO
      END IF
# endif
#endif
!
!  Initialize adjoint model state with a delta function at specified
!  point. Use USER parameters from standard input to perturb solution
!  in routine "ana_perturb". Then, convolve solution with the adjoint
!  diffusion operator.
!
      ADmodel=.TRUE.
      Lweak=.FALSE.

      DO ng=1,Ngrids
        Lnew(ng)=1
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_perturb (ng, tile, iADM)
#ifdef BALANCE_OPERATOR
          CALL ad_balance (ng, tile, Lbck, Lnew(ng))
          CALL ad_variability (ng, tile, Lnew(ng), Lweak)
#endif
          CALL ad_convolution (ng, tile, Lnew(ng), Lweak, 2)
        END DO
!$OMP END PARALLEL

        ADmodel=.FALSE.
!
!  Initialize tangent linear model with convolved adjoint solution.
!  Then, apply tangent linear convolution.
!
        add=.FALSE.
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL load_ADtoTL (ng, tile, Lnew(ng), Lnew(ng), add)
          CALL tl_convolution (ng, tile, Lnew(ng), Lweak, 2)
#ifdef BALANCE_OPERATOR
          CALL tl_variability (ng, tile, Lnew(ng), Lweak)
          CALL tl_balance (ng, tile, Lbck, Lnew(ng))
#endif
          CALL load_TLtoAD (ng, tile, Lnew(ng), Lnew(ng), add)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Write out background error correlation in adjoint history NetCDF
!  file.
!
      DO ng=1,Ngrids
        kstp(ng)=Lnew(ng)
#ifdef SOLVE3D
        nstp(ng)=Lnew(ng)
#endif
        LdefADJ(ng)=.TRUE.
        LwrtADJ(ng)=.TRUE.
        LwrtState2d(ng)=.TRUE.
        CALL ad_def_his (ng, LdefADJ(ng))
        IF (exit_flag.ne.NoError) RETURN
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Ladjusted(ng)=.TRUE.
#endif
        CALL ad_wrt_his (ng)
        IF (exit_flag.ne.NoError) RETURN
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Ladjusted(ng)=.FALSE.
#endif
      END DO

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS driver execution.                 !
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
!  If cycling restart records, write solution into record 3.
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
 20     FORMAT (/,'Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Close IO files.
!
      CALL close_out

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
