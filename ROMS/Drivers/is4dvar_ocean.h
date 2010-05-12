#include "cppdefs.h"
      MODULE ocean_control_mod
!
!svn $Id: is4dvar_ocean.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Strong Constraint 4-Dimensional Variational (4DVar)       !
!            Data Assimilation Driver, Incremental Approach:           !
!                                                                      !
!  This driver is used for strong constraint 4DVar where the only      !
!  errors considered are those for the observations. The model is      !
!  assumed to be perfect.  This is the incremental method and the      !
!  nonlinear, tangent linear and adjoint models are needed.            !
!                                                                      !
!  The misfit (squared difference) between model and observations      !
!  is defined as:                                                      !
!                                                                      !
!         J  = Jb + Jo                                                 !
!                                                                      !
!         Jb = 1/2 transpose(deltaX) * B^(-1) * (deltaX)               !
!                                                                      !
!         Jo = 1/2 transpose(H deltaX - d) * O^(-1) * (H deltaX - d)   !
!                                                                      !
!          d = Xo - Xb                                                 !
!                                                                      !
!     deltaX = X - Xb                                                  !
!                                                                      !
!  where                                                               !
!                                                                      !
!         B : background error covariance                              !
!         d : innovation vector                                        !
!    deltaX : increment vector propagated in time by the TLM.          !
!         H : linearized observation operator                          !
!         O : observation error covariance                             !
!         X : model at observation or background points                !
!        Xb : background or reference vector                           !
!        Xo : observations vector                                      !
!                                                                      !
!  The routines in this driver control the initialization,  time-      !
!  stepping, and finalization of  ROMS/TOMS  model following ESMF      !
!  conventions:                                                        !
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
!  This routine time-steps ROMS/TOMS nonlinear, tangent linear and     !
!  adjoint models.                                                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
#ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE back_cost_mod, ONLY : back_cost
      USE cgradient_mod, ONLY : cgradient
      USE cost_grad_mod, ONLY : cost_grad
      USE ini_adjust_mod, ONLY : ini_adjust
      USE ini_fields_mod, ONLY : ini_fields
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces, ONLY : initialize_forces
#endif
      USE mod_ocean, ONLY : initialize_ocean
      USE normalization_mod, ONLY : normalization
      USE sum_grad_mod, ONLY : sum_grad
#ifdef BALANCE_OPERATOR
      USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod, ONLY: balance_ref, biconj
#endif
!
!  Imported variable declarations
!
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      logical :: converged
      logical :: Lweak = .FALSE.

      integer :: my_inner, my_outer
      integer :: AdjRec, Lbck, Lini, Lsav, Rec1, Rec2, Rec3, Rec4
      integer :: i, my_iic, ng, subs, tile, thread
      integer :: Lcon, LTLM1, LTLM2, LTLM3, LADJ1, LADJ2
      integer :: NRMrec
      integer :: lstr, status

      real(r8) :: rate
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!-----------------------------------------------------------------------
!  OUTER LOOP: time-step nonlinear model.
!-----------------------------------------------------------------------
!
!  Initialize relevant parameters.
!
#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
        Lfinp(ng)=1         ! forcing index for input
        Lfout(ng)=1         ! forcing index for output history files
#endif
#ifdef ADJUST_BOUNDARY
        Lbinp(ng)=1         ! boundary index for input
        Lbout(ng)=1         ! boundary index for output history files
#endif
        Lold(ng)=1          ! old minimization time index
        Lnew(ng)=2          ! new minimization time index
        LTLM1=1             ! trial x-space TLM IC record in ITLname
        LTLM2=2             ! previous v-space TLM IC record in ITLname
        LTLM3=3             ! trial v-space TLM IC record in ITLname
        LADJ1=1             ! initial cost gradient
        LADJ2=2             ! new cost gradient (not normalized)
        Lini=1              ! NLM initial conditions record in INIname
        Lbck=2              ! background record in INIname
        Rec1=1
        Rec2=2
        Rec3=3
        Rec4=4
        Nrun=1
        ERstr=1
        ERend=Nouter

        OUTER_LOOP : DO my_outer=1,Nouter
          outer=my_outer
          inner=0
!
!  Set nonlinear output history file name. Create a basic state file
!  for each outher loop.
!
          LdefHIS(ng)=.TRUE.
          LwrtHIS(ng)=.TRUE.
          lstr=LEN_TRIM(FWDbase(ng))
          WRITE (HISname(ng),10) FWDbase(ng)(1:lstr-3), outer-1

#if defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!  Set file name containing the nonlinear model bulk fluxes to be read
!  and processed by other algorithms.
!
          IF (outer.eq.1) THEN
            BLKname(ng)=HISname(ng)
          END IF
#endif
!
!  Initialize nonlinear model. If outer=1, the model is initialized
!  with the background or reference state. Otherwise, the model is
!  initialized with the estimated initial conditions from previous
!  iteration, X(0) = X(0) + deltaX(0).
!
          wrtNLmod(ng)=.TRUE.
          wrtTLmod(ng)=.FALSE.
          tRSTindx(ng)=0
          NrecRST(ng)=0
          CALL initial (ng)
          IF (exit_flag.ne.NoError) RETURN
!
!  If first pass, save nonlinear initial conditions (currently in time
!  index 1, background) into next record (Lbck) of INIname NetCDF file.
!  The record "Lbck" becomes the background state record and the record
!  "Lini" becomes current nonlinear initial conditions.  Both record
!  are used in the algorithm below.
!
          IF (Nrun.eq.1) THEN
            tINIindx(ng)=1
            NrecINI(ng)=1
            CALL wrt_ini (ng, 1)
            IF (exit_flag.ne.NoError) RETURN
          END IF

#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
          IF (balance(isFsur)) THEN
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1
                CALL balance_ref (ng, TILE, Lini)
                CALL biconj (ng, TILE, iNLM, Lini)
              END DO
            END DO
!$OMP END PARALLEL DO
            wrtZetaRef(ng)=.TRUE.
          END IF
#endif
!
!  If first pass, compute or read in background-error covariance
!  normalization factors. If computing, write out factors to
!  NetCDF. This is an expensive computation that needs to be
!  computed only once for a particular application grid.
!
          IF (Nrun.eq.1) THEN
            IF (ANY(LwrtNRM(:,ng))) THEN
              CALL def_norm (ng, iNLM, 1)
              IF (exit_flag.ne.NoError) RETURN

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

#ifdef ADJUST_BOUNDARY
              CALL get_state (ng, 10, 10, NRMname(3,ng), NRMrec, 1)
              IF (exit_flag.ne.NoError) RETURN
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
              CALL get_state (ng, 11, 11, NRMname(4,ng), NRMrec, 1)
              IF (exit_flag.ne.NoError) RETURN
#endif
            END IF
          END IF
!
!  If first pass, define output 4DVAR NetCDF file containing all
!  processed data at observation locations.
!
          IF (Nrun.eq.1) THEN
            LdefMOD(ng)=.TRUE.
            CALL def_mod (ng)
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  Run nonlinear model. Save nonlinear tracjectory needed by the
!  adjoint and tangent linear models. Interpolate nonlinear model
!  to boservation locations (compute and save H x).
!
          IF (Master) THEN
            WRITE (stdout,20) 'NL', ntstart(ng), ntend(ng)
          END IF

          time(ng)=time(ng)-dt(ng)

          NL_LOOP1 : DO my_iic=ntstart(ng),ntend(ng)+1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL main3d (ng)
#else
            CALL main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO NL_LOOP1
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.TRUE.

#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
!
!  Write out initial and background surface forcing into initial
!  INIname NetCDF file for latter use.
!
          CALL wrt_frc (ng, Lfout(ng), Lini)
          IF (exit_flag.ne.NoError) RETURN
          IF (Nrun.eq.1) THEN
            CALL wrt_frc (ng, Lfout(ng), Lbck)
            IF (exit_flag.ne.NoError) RETURN
          END IF
#endif
!
!  Write out nonlinear model misfit cost function into MODname NetCDF
!  file.
!
          SourceFile='is4dvar_ocean.h, ROMS_run'

          CALL netcdf_put_fvar (ng, iNLM, MODname(ng),                  &
     &                          'NLcost_function',                      &
     &                          FOURDVAR(ng)%NLobsCost(0:),             &
     &                          (/1,outer/), (/NstateVar(ng)+1,1/),     &
     &                          ncid = ncMODid(ng))
          IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  INNER LOOP: iterate using tangent linear model increments.
!-----------------------------------------------------------------------
!
!  The minimization algorithm requires to save all the gradient
!  solutions for each inner loop iteration.  They are used for
!  orthogonalization in the conjugate gradient algorithm.  Thus,
!  we need to reset adjoint file record indices.
!
          tADJindx(ng)=0
          NrecADJ(ng)=0
!
!  An adjoint NetCDF is created for each outer loop.
!
          LdefADJ(ng)=.TRUE.
          lstr=LEN_TRIM(ADJbase(ng))
          WRITE (ADJname(ng),10) ADJbase(ng)(1:lstr-3), outer
!
!  Define output Hessian NetCDF file containing the eigenvectors
!  approximation to the Hessian matrix computed from the Lanczos
!  algorithm. Notice that the file name is a function of the
!  outer loop. That is, a file is created for each outer loop.
!
          lstr=LEN_TRIM(HSSbase(ng))
          WRITE (HSSname(ng),10) HSSbase(ng)(1:lstr-3), outer
          LdefHSS(ng)=.TRUE.
          CALL def_hessian (ng)
          IF (exit_flag.ne.NoError) RETURN
!
!  Notice that inner loop iteration start from zero. This is needed to
!  compute the minimization initial increment deltaX(0), its associated
!  gradient G(0), and descent direction d(0) used in the conjugate
!  gradient algorithm.
!
          INNER_LOOP : DO my_inner=0,Ninner
            inner=my_inner
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time-step tangent linear model: compute cost function.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  If first pass inner=0, initialize tangent linear state (increments,
!  deltaX) from rest. Otherwise, use trial initial conditions estimated
!  by the conjugate gradient algorithm in previous inner loop. The TLM
!  initial conditions are read from ITLname, record 1.
!
            tITLindx(ng)=1
            CALL tl_initial (ng)
            IF (exit_flag.ne.NoError) RETURN
!
!  On first pass, initialize records 2, 3 and 4 of the ITL file to zero.
!
            IF (inner.eq.0.and.outer.eq.1) THEN
              CALL tl_wrt_ini (ng, LTLM1, Rec2)
              IF (exit_flag.ne.NoError) RETURN
              CALL tl_wrt_ini (ng, LTLM1, Rec3)
              IF (exit_flag.ne.NoError) RETURN
              CALL tl_wrt_ini (ng, LTLM1, Rec4)
              IF (exit_flag.ne.NoError) RETURN
            END IF

#ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, activate writing and determine
!  output file name. The multiple file option is use to perturb initial
!  state and create ensembles.  The TLM final  trajectory is written for
!  each inner loop on separated NetCDF files.
!
            LdefTLM(ng)=.TRUE.
            LwrtTLM(ng)=.TRUE.
            lstr=LEN_TRIM(TLMbase(ng))
            WRITE (TLMname(ng),10) TLMbase(ng)(1:lstr-3), Nrun
#endif
!
!  Activate switch to write out initial and final misfit between
!  model and observations.
!
            wrtMisfit(ng)=.FALSE.
            IF (((outer.eq.1).and.(inner.eq.0)).or.                     &
     &          ((outer.eq.Nouter).and.(inner.eq.Ninner))) THEN
              wrtMisfit(ng)=.TRUE.
            END IF
!
!  Run tangent linear model. Compute misfit observation cost function,
!  Jo.
!
            IF (Master) THEN
              WRITE (stdout,20) 'TL', ntstart(ng), ntend(ng)
            END IF

            time(ng)=time(ng)-dt(ng)

            TL_LOOP : DO my_iic=ntstart(ng),ntend(ng)+1

              iic(ng)=my_iic
#ifdef SOLVE3D
              CALL tl_main3d (ng)
#else
              CALL tl_main2d (ng)
#endif
              IF (exit_flag.ne.NoError) RETURN

            END DO TL_LOOP

#ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, close current NetCDF file.
!
            IF (ncTLMid(ng).ne.-1) THEN
              SourceFile='is4dvar_ocean.h, ROMS_run'

              CALL netcdf_close (ng, iTLM, ncTLMid(ng))
              IF (exit_flag.ne.NoError) RETURN
            END IF
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time step adjoint model backwards: compute cost function gradient.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize the adjoint model always from rest.
!
            CALL ad_initial (ng)
            IF (exit_flag.ne.NoError) RETURN
!
!  Time-step adjoint model backwards. The adjoint model is forced with
!  the adjoint of the observation misfit (Jo) term.
!
            IF (Master) THEN
              WRITE (stdout,20) 'AD', ntstart(ng), ntend(ng)
            END IF

            time(ng)=time(ng)+dt(ng)

            AD_LOOP : DO my_iic=ntstart(ng),ntend(ng),-1

              iic(ng)=my_iic
#ifdef SOLVE3D
              CALL ad_main3d (ng)
#else
              CALL ad_main2d (ng)
#endif
              IF (exit_flag.ne.NoError) RETURN

            END DO AD_LOOP
!
!  Clear adjoint arrays.  Is it needed?
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
            DO thread=0,numthreads-1
#if defined _OPENMP || defined DISTRIBUTE
              subs=NtileX(ng)*NtileE(ng)/numthreads
#else
              subs=1
#endif
              DO tile=subs*thread,subs*(thread+1)-1
                CALL initialize_ocean (ng, TILE, iADM)
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
                CALL initialize_forces (ng, TILE, iADM)
#endif
              END DO
            END DO
!$OMP END PARALLEL DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Descent algorithm.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read TLM v-space initial conditions, record 3 in ITLname, and
!  load it into time index LTLM1. This is needed to compute background
!  cost function. Also read in new (x-space) gradient vector, GRADx(Jo),
!  from adjoint history file ADJname.  Read in the sum of all the
!  previous outer-loop increments which are always in record 4 of
!  the ITL file.
!
            IF (inner.eq.0) THEN
              CALL get_state (ng, iTLM, 8, ITLname(ng), Rec1, LTLM1)
              IF (exit_flag.ne.NoError) RETURN
            ELSE
              CALL get_state (ng, iTLM, 8, ITLname(ng), Rec3, LTLM1)
              IF (exit_flag.ne.NoError) RETURN
            END IF
            CALL get_state (ng, iTLM, 8, ITLname(ng), Rec4, LTLM2)
            IF (exit_flag.ne.NoError) RETURN
            CALL get_state (ng, iADM, 4, ADJname(ng), tADJindx(ng),     &
     &                      LADJ2)
            IF (exit_flag.ne.NoError) RETURN
#ifdef BALANCE_OPERATOR
            CALL get_state (ng, iNLM, 2, INIname(ng), Lini, Lini)
            IF (exit_flag.ne.NoError) RETURN
            nrhs(ng)=Lini
#endif
!
!  Convert observation cost function gradient, GRADx(Jo), from model
!  space (x-space) to minimization space (v-space):
!
!     GRADv(Jo) = B^(T/2) GRADx(Jo),  operator: S G L^(T/2) W^(-1/2)
!
!  First, multiply the adjoint solution, GRADx(Jo), by the background-
!  error standard deviations, S.  Second, convolve result with the
!  adjoint diffusion operator, G L^(T/2) W^(-1/2). Then, backgound
!  cost function contribution (BackCost) and cost function gradient
!  (v-space) by adding background and observation contributions:
!
!     GRADv(J) = GRADv(Jb) + GRADv(Jo) = deltaV + GRADv(Jo)
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1
#ifdef BALANCE_OPERATOR
                CALL ad_balance (ng, TILE, Lini, LADJ2)
#endif
                CALL ad_variability (ng, TILE, LADJ2, Lweak)
                CALL ad_convolution (ng, TILE, LADJ2, Lweak, 2)
                CALL cost_grad (ng, TILE, LTLM1, LTLM2, LADJ2)
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  Compute current total cost function.
!
            IF (Nrun.eq.1) THEN
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostNorm(i)
                FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostNorm(i)
              END DO
            ELSE
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
              END DO
            END IF
!
!  Prepare for background cost function (Jb) calculation:
!
!  Read the convolved gradient from inner=0 (which is permanently
!  saved in record 1 of the adjoint file)  ALWAYS into record 1.
!
            IF (inner.gt.0) THEN
              CALL get_state (ng, iADM, 3, ADJname(ng), LADJ1,          &
     &                        LADJ1)
              IF (exit_flag.ne.NoError) RETURN
            END IF
!
!  Compute background cost function (Jb) for inner=0:
!
!  If first pass of inner loop, read in the sum of previous v-space
!  gradients from record 4 of ITL file using the TLM model variables
!  as temporary storage. Also add background cost function to Cost0.
!
            IF (inner.eq.0) THEN
              CALL get_state (ng, iTLM, 2, ITLname(ng), Rec4, LTLM2)
              IF (exit_flag.ne.NoError) RETURN
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1
                  CALL back_cost (ng, TILE, LTLM2)
                END DO
              END DO
!$OMP END PARALLEL DO
!
              FOURDVAR(ng)%Cost0(outer)=FOURDVAR(ng)%Cost0(outer)+      &
     &                                  FOURDVAR(ng)%BackCost(0)
            END IF
!
!  Compute current total cost function.
!
            IF (Nrun.eq.1) THEN
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostNorm(i)=FOURDVAR(ng)%CostNorm(i)+      &
     &                                   FOURDVAR(ng)%BackCost(i)
                FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostNorm(i)
                FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostNorm(i)
              END DO
            ELSE
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
              END DO
            END IF
!
!  Determine the descent direction in which the quadractic total cost
!  function decreases. Then, determine the TLM initial conditions,
!  deltaV(LTLM1), and its gradient, GRADv{J(Lnew)} at the new
!  direction.  Also, Compute TLM v-space trial initial conditions for
!  next inner loop, deltaV(LTLM2). The new gradient minimize the
!  quadratic cost function spanned by current and previous inner loop
!  iterations.  This is achieved by orthogonalizing (Gramm-Schmidt
!  algorithm) against all previous inner loop gradients.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,outer,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1
                CALL cgradient (ng, TILE, iTLM, inner, outer)
              END DO
            END DO
!$OMP END PARALLEL DO
            IF (exit_flag.ne.NoError) RETURN
!
!  Report background (Jb) and observations (Jo) cost function values
!  normalized by their first minimization value. It also reports the
!  percentage change on total cost function value with respect to
!  previous iteration. Compute the optimality of the minimization to
!  check the statistical hypotheses between the background and
!  observations errors: the cost function value at the minimum, Jmin,
!  is idealy equal to half the number of observations assimilated
!  (Optimality=1=2*Jmin/Nobs), for a linear system.
!
            IF (Master) THEN
              IF (Nrun.gt.1) THEN
                rate=100.0_r8*ABS(FOURDVAR(ng)%CostFun(0)-              &
     &                            FOURDVAR(ng)%CostFunOld(0))/          &
     &                        FOURDVAR(ng)%CostFunOld(0)
              ELSE
                rate=0.0_r8
              END IF
              Optimality(ng)=2.0_r8*FOURDVAR(ng)%CostFun(0)/            &
     &                       (FOURDVAR(ng)%ObsCount(0)-                 &
     &                        FOURDVAR(ng)%ObsReject(0))
              WRITE (stdout,30) outer, inner,                           &
     &                          FOURDVAR(ng)%BackCost(0)/               &
     &                          FOURDVAR(ng)%CostNorm(0),               &
     &                          FOURDVAR(ng)%ObsCost(0)/                &
     &                          FOURDVAR(ng)%CostNorm(0),               &
     &                          rate
              IF (inner.eq.0) THEN
                DO i=0,NstateVar(ng)
                  IF (FOURDVAR(ng)%NLobsCost(i).ne.0.0_r8) THEN
                    IF (i.eq.0) THEN
                      WRITE (stdout,40) outer, inner,                   &
     &                                  FOURDVAR(ng)%NLobsCost(i)/      &
     &                                  FOURDVAR(ng)%CostNorm(i)
                    ELSE
                      WRITE (stdout,50) outer, inner,                   &
     &                                  FOURDVAR(ng)%NLobsCost(i)/      &
     &                                  FOURDVAR(ng)%CostNorm(i),       &
     &                                  TRIM(Vname(1,idSvar(i)))
                    END IF
                  END IF
                  FOURDVAR(ng)%NLobsCost(i)=0.0
                END DO
              END IF
              WRITE (stdout,60) outer, inner, Optimality(ng)
            END IF
!
!  Save total v-space cost function gradient, GRADv{J(Lnew)}, into
!  ADJname history NetCDF file. Noticed that the lastest adjoint
!  solution record is over-written in the NetCDF file for future use.
!  The switch "LwrtState2d" is activated to write out state arrays
!  instead ad_*_sol arrays.
!
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
            Lfout(ng)=LADJ2
#endif
#ifdef ADJUST_BOUNDARY
            Lbout(ng)=LADJ2
#endif
            kstp(ng)=LADJ2
#ifdef SOLVE3D
            nstp(ng)=LADJ2
#endif
            tADJindx(ng)=tADJindx(ng)-1
            LwrtState2d(ng)=.TRUE.
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
            LwrtState2d(ng)=.FALSE.
!
!  Write out trial v-space TLM initial conditions, currently in time
!  index LTM2, into record 3 of ITLname NetCDF file.
!
            CALL tl_wrt_ini (ng, LTLM2, Rec3)
            IF (exit_flag.ne.NoError) RETURN
!
!  Read current outer loop nonlinear model initial conditions and
!  background state vectors.
!
            CALL get_state (ng, iNLM, 2, INIname(ng), Lini, Lini)
            IF (exit_flag.ne.NoError) RETURN
            CALL get_state (ng, iNLM, 9, INIname(ng), Lbck, Lbck)
            IF (exit_flag.ne.NoError) RETURN
!
!  Convert increment vector, deltaV, from minimization space (v-space)
!  to model space (x-space):
!
!     deltaX = B^(1/2) deltaV
!  or
!     deltaX = W^(1/2) L^(1/2) G S
!
!  First, convolve estimated increment vector (v-space) by with the
!  tangent linear diffusion operator, W^(1/2) L^(1/2) G.  Second,
!  multiply result by the background-error standard deviation, S.
!
            Lcon=LTLM2
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL tl_convolution (ng, TILE, Lcon, Lweak, 2)
                CALL tl_variability (ng, TILE, Lcon, Lweak)
#ifdef BALANCE_OPERATOR
                CALL tl_balance (ng, TILE, Lini, Lcon)
#endif
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  Write out trial x-space (convolved) TLM initial conditions, currently
!  in time index Lcon, into record 1 of ITLname NetCDF file.
!
            CALL tl_wrt_ini (ng, Lcon, Rec1)
            IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Update counters.
!-----------------------------------------------------------------------
!
            Lsav=Lnew(ng)
            Lnew(ng)=Lold(ng)
            Lold(ng)=Lsav
            Nrun=Nrun+1

          END DO INNER_LOOP
!
!  Close adjoint NetCDF file.
!
          IF (ncADJid(ng).ne.-1) THEN
            SourceFile='is4dvar_ocean.h, ROMS_run'

            CALL netcdf_close (ng, iADM, ncADJid(ng))
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  Close Hessian NetCDF file.
!
          IF (ncHSSid(ng).ne.-1) THEN
            SourceFile='is4dvar_ocean.h, ROMS_run'

            CALL netcdf_close (ng, iADM, ncHSSid(ng))
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!-----------------------------------------------------------------------
!  Clear nonlinear state variables.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
#if defined _OPENMP || defined DISTRIBUTE
            subs=NtileX(ng)*NtileE(ng)/numthreads
#else
            subs=1
#endif
            DO tile=subs*thread,subs*(thread+1)-1
              CALL initialize_ocean (ng, TILE, iNLM)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Compute new nonlinear initial conditions by adding minimization
!  increment to previous outer loop initial conditions:
!
!         Xi(outer+1) = Xi(outer) + deltaX(Lcon)
!
!-----------------------------------------------------------------------
!
!  Notice that "ini_fields" is called here for output purposes only.
!  It computes the vertically integrated momentum in 3D applications.
!  In order to use the correct fields, the model time indices are set
!  to Lini.
!  The appropriate tl correction for the NL model resides in record 1
!  of the ITL file.
!
          kstp(ng)=Lini
#ifdef SOLVE3D
          nstp(ng)=Lini
#endif
          CALL get_state (ng, iNLM, 1, INIname(ng), Lini, Lini)
          IF (exit_flag.ne.NoError) RETURN
          CALL get_state (ng, iTLM, 1, ITLname(ng), LTLM1, LTLM1)
          IF (exit_flag.ne.NoError) RETURN

!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL ini_adjust (ng, TILE, LTLM1, Lini)
              CALL ini_fields (ng, TILE, iNLM)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  Write out new nonlinear model initial conditions into record Lini
!  of INIname.
!
          tINIindx(ng)=0
          NrecINI(ng)=1
          CALL wrt_ini (ng, Lini)
          IF (exit_flag.ne.NoError) RETURN
!
! Gather the v-space increments from the final inner-loop and
! save in record 4 of the ITL file. The current v-space increment
! is in record 3 and the sum so far is in record 4.
!
          CALL get_state (ng, iTLM, 8, ITLname(ng), Rec3, LTLM1)
          IF (exit_flag.ne.NoError) RETURN
          CALL get_state (ng, iTLM, 8, ITLname(ng), Rec4, LTLM2)
          IF (exit_flag.ne.NoError) RETURN
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL sum_grad (ng, TILE, LTLM1, LTLM2)
            END DO
          END DO
!$OMP END PARALLEL DO
!
! Write the current sum into record 4 of the ITL file.
!
          CALL tl_wrt_ini (ng, LTLM2, Rec4)
          IF (exit_flag.ne.NoError) RETURN

#if defined ADJUST_STFLUX   || defined ADJUST_WSTRESS || \
    defined ADJUST_BOUNDARY
!
!  Set index containing the surface forcing increments used the run
!  the nonlinear model in the outer loop and read the forcing
!  increments. For bulk fluxes, we read Rec1 because the stress
!  fluxes change by virtue of the changing initial conditions.
!  When not using bulk fluxes, we read Rec4 because the background
!  stress and flux is prescribed by input files which are not
!  overwritten so we need to correct the background using the
!  sum of the increments from all previous outer-loops.
!  If using Rec4 we need to convert from v-space to x-space
!  by applying the convolution.
!  Note that Lfinp=Lbinp so the the forcing and boundary
!  adjustments are both processsed correctly.
# ifdef BALANCE_OPERATOR
!  Currently, We don't need the call to tl_balance below, but we
!  might later if we impose a balance constraint on the wind stress
!  corrections.
# endif
!
!  AMM: CHECK WHAT HAPPENS WITH SECONDARY PRECONDITIONING.
!
          Lfinp(ng)=LTLM1
# ifdef BULK_FLUXES
          CALL get_state (ng, iTLM, 1, ITLname(ng), Rec1, Lfinp(ng))
# else
          CALL get_state (ng, iTLM, 1, ITLname(ng), Rec4, Lfinp(ng))
          Lcon=Lfinp(ng)
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL tl_convolution (ng, TILE, Lcon, Lweak, 2)
              CALL tl_variability (ng, TILE, Lcon, Lweak)
# ifdef BALANCE_OPERATOR
!!            CALL tl_balance (ng, TILE, Lini, Lcon)
# endif
            END DO
          END DO
!$OMP END PARALLEL DO
# endif
          IF (exit_flag.ne.NoError) RETURN
#endif
!
!-----------------------------------------------------------------------
!  Clear tangent linear state variables.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
#if defined _OPENMP || defined DISTRIBUTE
            subs=NtileX(ng)*NtileE(ng)/numthreads
#else
            subs=1
#endif
            DO tile=subs*thread,subs*(thread+1)-1
              CALL initialize_ocean (ng, TILE, iTLM)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  Close current forward NetCDF file.
!
          SourceFile='is4dvar_ocean.h, ROMS_run'

          CALL netcdf_close (ng, iNLM, ncFWDid(ng))
          IF (exit_flag.ne.NoError) RETURN

        END DO OUTER_LOOP
!
!-----------------------------------------------------------------------
!  Done with data assimilation. Initialize the nonlinear model with
!  estimated initial conditions. Save nonlinear solution at observation
!  points for posterior analysis.
!-----------------------------------------------------------------------
!
!  Set nonlinear output history file name. Create a basic state file
!  for each outher loop.
!
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        tHISindx(ng)=0
        NrecHIS(ng)=0
        lstr=LEN_TRIM(FWDbase(ng))
        WRITE (HISname(ng),10) FWDbase(ng)(1:lstr-3), Nouter
!
!  Initialize nonlinear model with estimated initial conditions.
!
        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.
        wrtMisfit(ng)=.FALSE.
        tRSTindx(ng)=0
        NrecRST(ng)=0
        CALL initial (ng)
        IF (exit_flag.ne.NoError) RETURN
!
! Clear NLobsCost.
!
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%NLobsCost(i)=0.0
        END DO
!
!  Run nonlinear model. Interpolate nonlinear model to observation
!  locations.
!
        IF (Master) THEN
          WRITE (stdout,20) 'NL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        NL_LOOP2 : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL main3d (ng)
#else
          CALL main2d (ng)
#endif
          IF (exit_flag.ne.NoError) RETURN

        END DO NL_LOOP2
!
!  Write out nonlinear model final misfit cost function into MODname
!  NetCDF file. Notice that it is written in the Nouter+1 record.
!
        SourceFile='is4dvar_ocean.h, ROMS_run'

        CALL netcdf_put_fvar (ng, iNLM, MODname(ng), 'NLcost_function', &
     &                        FOURDVAR(ng)%NLobsCost(0:),               &
     &                        (/1,Nouter+1/), (/NstateVar(ng)+1,1/),    &
     &                        ncid = ncMODid(ng))
        IF (exit_flag.ne.NoError) RETURN
!
!  Report the final value of the nonlinear model misfit cost function.
!
        IF (Master) THEN
          DO i=0,NstateVar(ng)
            IF (FOURDVAR(ng)%NLobsCost(i).ne.0.0_r8) THEN
              IF (i.eq.0) THEN
                WRITE (stdout,40) outer, inner,                         &
     &                            FOURDVAR(ng)%NLobsCost(i)/            &
     &                            FOURDVAR(ng)%CostNorm(i)
              ELSE
                WRITE (stdout,50) outer, inner,                         &
     &                            FOURDVAR(ng)%NLobsCost(i)/            &
     &                            FOURDVAR(ng)%CostNorm(i),             &
     &                            TRIM(Vname(1,idSvar(i)))
              END IF
            END IF
          END DO
        END IF
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
        ncHISid(ng)=-1
!
!  Compute and report model-observation comparison statistics.
!
        CALL stats_modobs (ng)

      END DO NEST_LOOP
!
 10   FORMAT (a,'_',i3.3,'.nc')
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 30   FORMAT (/,' (',i3.3,',',i3.3,'): TLM Cost Jb, J  = ',             &
     &        1p,e16.10,0p,1x,1p,e16.10,0p,t68,1p,e10.4,' %')
 40   FORMAT (/,'>(',i3.3,',',i3.3,'): NLM Cost     J  = ',             &
     &        17x,1p,e16.10,0p)
 50   FORMAT (' (',i3.3,',',i3.3,'): NLM Cost     J  = ',               &
     &        17x,1p,e16.10,0p,t68,a)
 60   FORMAT (/,1x,'(',i3.3,',',i3.3,'): Optimality (2*J/Nobs) = ',     &
     &        1p,e16.10,/)

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
