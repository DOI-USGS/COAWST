      MODULE ocean_control_mod
!
!svn $Id: correlation.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
#ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_atmos_coupling
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_waves_coupling
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

      integer :: STDrec, Tindex, ng, thread

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
 10       FORMAT (' Process Information:',/)
        END IF
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            CALL wclock_on (ng, iNLM, 0)
          END DO
!$OMP END PARALLEL DO
        END DO

#if defined AIR_OCEAN || defined WAVES_OCEAN
!
!  Initialize coupling streams between model(s).
!
        DO ng=1,Ngrids
# ifdef AIR_OCEAN
          CALL initialize_atmos_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
          CALL initialize_waves_coupling (ng, MyRank)
# endif
        END DO
#endif
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
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar
!
!  Read in standard deviation factors for initial conditions
!  error covariance.  They are loaded in Tindex=1 of the
!  e_var(...,Tindex) state variables.
!
        STDrec=1
        Tindex=1
        DO ng=1,Ngrids
          CALL get_state (ng, 6, 6, STDname(1,ng), STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Read in standard deviation factors for model error covariance.
!  They are loaded in Tindex=2 of the e_var(...,Tindex) state
!  variables.
!
        STDrec=1
        Tindex=2
        DO ng=1,Ngrids
          IF (NSA.eq.2) THEN
            CALL get_state (ng, 6, 6, STDname(2,ng), STDrec, Tindex)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

#ifdef ADJUST_BOUNDARY
!
!  Read in standard deviation factors for boundary conditions
!  error covariance.
!
        STDrec=1
        Tindex=1
        DO ng=1,Ngrids
          CALL get_state (ng, 8, 8, STDname(3,ng), STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END DO
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Read in standard deviation factors for surface forcing
!  error covariance.
!
        STDrec=1
        Tindex=1
        DO ng=1,Ngrids
          CALL get_state (ng, 9, 9, STDname(4,ng), STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END DO
#endif
      END IF

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Tend)
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
      USE mod_ncparam
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
      USE normalization_mod, ONLY : normalization
#ifdef BALANCE_OPERATOR
      USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
!
!  Imported variable declarations
!
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      logical :: Lweak, add
      integer :: i, ng, subs, tile, thread
#ifdef BALANCE_OPERATOR
      integer :: Lbck = 1
#endif
      integer :: NRMrec
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!-----------------------------------------------------------------------
!  Initialize metrics.
!-----------------------------------------------------------------------
!
        CALL initial (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Get background-error covariance normalization matrix.
!-----------------------------------------------------------------------
!
!  Compute or read in background-error covariance normalization factors.
!  If computing, write out factors to NetCDF. This is an expensive
!  computation and needs to be computed once for an application grid.
!
        IF (ANY(LwrtNRM(:,ng))) THEN
          CALL def_norm (ng, iNLM, 1)
          IF (exit_flag.ne.NoError) RETURN

          IF (NSA.eq.2) THEN
            CALL def_norm (ng, iNLM, 2)
          IF (exit_flag.ne.NoError) RETURN
          END IF
#ifdef ADJUST_BOUNDARY
          CALL def_norm (ng, iNLM, 3)
          IF (exit_flag.ne.NoError) RETURN
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL def_norm (ng, iNLM, 4)
          IF (exit_flag.ne.NoError) RETURN
#endif
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL normalization (ng, TILE, 2)
            END DO
          END DO
!$OMP END PARALLEL DO
          LdefNRM(1:4,ng)=.FALSE.
          LwrtNRM(1:4,ng)=.FALSE.
        ELSE
          NRMrec=1
          CALL get_state (ng, 5, 5, NRMname(1,ng), NRMrec, 1)
          IF (exit_flag.ne.NoError) RETURN

          IF (NSA.eq.2) THEN
            CALL get_state (ng, 5, 5, NRMname(2,ng), NRMrec, 2)
            IF (exit_flag.ne.NoError) RETURN
          END IF
#ifdef ADJUST_BOUNDARY
          CALL get_state (ng, 10, 10, NRMname(3,ng), NRMrec, 1)
          IF (exit_flag.ne.NoError) RETURN
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL get_state (ng, 11, 11, NRMname(4,ng), NRMrec, 1)
          IF (exit_flag.ne.NoError) RETURN
#endif
        END IF

#ifdef BALANCE_OPERATOR
!
!-----------------------------------------------------------------------
!  Read background state.
!-----------------------------------------------------------------------
!
        CALL get_state (ng, iNLM, 9, FWDname(ng), Lbck, Lbck)
        IF (exit_flag.ne.NoError) RETURN
#endif
!
!-----------------------------------------------------------------------
!  Test correlation model.
!-----------------------------------------------------------------------
!
!  Initialize adjoint model state with a delta function at specified
!  point. Use USER parameters from standard input to perturb solution
!  in routine "ana_perturb". Then, convolve solution with the adjoint
!  diffusion operator.
!
        ADmodel=.TRUE.
        Lweak=.FALSE.
        Lnew(ng)=1

!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lbck) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL ana_perturb (ng, TILE, iADM)
#ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, TILE, Lbck, Lnew(ng))
            CALL ad_variability (ng, TILE, Lnew(ng), Lweak)
#endif
            CALL ad_convolution (ng, TILE, Lnew(ng), Lweak, 2)
          END DO
        END DO
!$OMP END PARALLEL DO
        ADmodel=.FALSE.
!
!  Initialize tangent linear model with convolved adjoint solution.
!  Then, apply tangent linear convolution.
!
        add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,add,thread,subs,tile,Lbck)
!$OMP&            SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL load_ADtoTL (ng, TILE, Lnew(ng), Lnew(ng), add)
            CALL tl_convolution (ng, TILE, Lnew(ng), Lweak, 2)
#ifdef BALANCE_OPERATOR
            CALL tl_variability (ng, TILE, Lnew(ng), Lweak)
            CALL tl_balance (ng, TILE, Lbck, Lnew(ng))
#endif
            CALL load_TLtoAD (ng, TILE, Lnew(ng), Lnew(ng), add)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Write out background error correlation in adjoint history NetCDF
!  file.
!
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

      END DO NEST_LOOP

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear, tangent linear, and    !
!  adjoint models execution.                                           !
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
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into record 3.
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
 20     FORMAT (/,'Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(ng,thread) SHARED(numthreads)
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
