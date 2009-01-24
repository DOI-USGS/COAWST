      MODULE ocean_control_mod
!
!svn $Id: is4dvar_ocean_old.h 652 2008-07-24 23:20:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
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

      integer :: STDrec, ng, thread

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
!  Read in background-error standard deviation factors and spatial
!  convolution diffusion coefficients.
!  
        STDrec=1
        DO ng=1,Ngrids
          CALL get_state (ng, 6, 6, STDname(ng), STDrec, 1)
          IF (exit_flag.ne.NoError) RETURN
        END DO

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
#ifdef MULTIPLE_TLM
      USE mod_netcdf
#endif
      USE mod_scalars
      USE mod_stepping
!
#ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE back_cost_mod, ONLY : back_cost
      USE back_step_mod, ONLY : back_step
      USE cost_grad_mod, ONLY : cost_grad
      USE cost_norm_mod, ONLY : cost_norm
      USE descent_mod, ONLY : descent
      USE ini_adjust_mod, ONLY : ini_adjust
      USE ini_fields_mod, ONLY : ini_fields
      USE mod_ocean, ONLY : initialize_ocean
      USE normalization_mod, ONLY : normalization
#ifdef BALANCE_OPERATOR
      USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_ini_adjust_mod, ONLY : tl_ini_adjust
      USE tl_variability_mod, ONLY : tl_variability
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

      integer :: AdjRec, Iold, Inew, Isav, Lbck, Lini, Litl, Lsav
      integer :: i, lstr, my_iic, ng, status, subs, tile, thread

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
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Lfout(ng)=1         ! forcing index for output history files
#endif
        Lold(ng)=1          ! old minimization time index
        Lnew(ng)=2          ! new minimization time index
        Lini=1              ! NLM initial conditions record in INIname
        Lbck=2              ! background record in INIname
        Litl=1              ! TLM initial conditions record in ITLname
        Iold=2              ! old increment (v-space) record in ITLname
        Inew=3              ! new increment (v-space) record in ITLname
        Nrun=1
        ERstr=1
        ERend=Nouter
        IterSD=0
        Ipass=0

        OUTER_LOOP : DO outer=1,Nouter
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
!  index 1, background) into record "Lini" of INIname NetCDF file. The
!  record "Lbck" becomes the background state record and the record
!  "Lini" becomes current nonlinear initial conditions.  Both record
!  are used in the algorithm below.
!
          IF (Nrun.eq.1) THEN
            IF (LcycleINI(ng)) THEN
              tINIindx(ng)=1
              NrecINI(ng)=1
            END IF
            CALL wrt_ini (ng, 1)
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  If first pass, compute or read in background-error covariance
!  normalization factors. If computing, write out factors to
!  NetCDF. This is an expensive computation that needs to be
!  computed only once for a particular application grid.
!  
          IF (Nrun.eq.1) THEN
            IF (LwrtNRM(ng)) THEN
              CALL def_norm (ng)
              IF (exit_flag.ne.NoError) RETURN
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1
                  CALL normalization (ng, TILE, 2)
                END DO
              END DO
!$OMP END PARALLEL DO
              LdefNRM(ng)=.FALSE.
              LwrtNRM(ng)=.FALSE.
            ELSE
              tNRMindx(ng)=1
              CALL get_state (ng, 5, 5, NRMname(ng), tNRMindx(ng), 1)
              IF (exit_flag.ne.NoError) RETURN
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
!
!-----------------------------------------------------------------------
!  INNER LOOP: iterate using tangent linear model increments.
!-----------------------------------------------------------------------
!
          IterSD=0
          CGstepF=CGstepI
          converged=.FALSE.

          INNER_LOOP : DO inner=1,Ninner

            CONVERGE : IF (.not.converged) THEN
!
!  Set cost function convergence switch.
!
              IF (inner.gt.1) THEN
                converged=ABS(FOURDVAR(ng)%CostFunOld(0)-               &
     &                        FOURDVAR(ng)%CostFun(0)).le.              &
     &                    CostFunFac*FOURDVAR(ng)%CostFun(0)
              END IF
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time-step tangent linear model: compute cost function.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  If first pass (inner=1, Ipass=1), initialize tangent linear state
!  (increments, deltaX) from rest. Otherwise, initialize with previous 
!  estimate.
! 
              Ipass=1
              CALL tl_initial (ng)
              IF (exit_flag.ne.NoError) RETURN

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
              WRITE (TLMname(ng),110) TLMbase(ng)(1:lstr-3), Nrun
#endif
!
!  Activate switch to write out initial and final misfit between
!  model and observations.
!
              wrtMisfit(ng)=.FALSE.
              IF (((outer.eq.1).and.(inner.eq.1)).or.                   &
     &            ((outer.eq.Nouter).and.(inner.eq.Ninner)).or.         &
     &            ((outer.eq.Nouter).and.converged)) THEN
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

              TL_LOOP1 : DO my_iic=ntstart(ng),ntend(ng)+1

                iic(ng)=my_iic
#ifdef SOLVE3D
                CALL tl_main3d (ng)
#else
                CALL tl_main2d (ng)
#endif
                IF (exit_flag.ne.NoError) RETURN

              END DO TL_LOOP1
              wrtMisfit(ng)=.FALSE.
!
!  Report background (Jb, previous value) and observations (Jo) cost
!  function values normalized by their first minimizationvalue. It
!  also reports the percentage change on total cost function value
!  with respect to previous iteration. Compute the optimality of the
!  minimization to check the statistical hypotheses between the
!  background and observations errors: the cost function value at
!  the minimum, Jmin, is idealy equal to half the number of
!  observations assimilated (Optimality=1=2*Jmin/Nobs), for a
!  linear system.
!
              IF (Master) THEN        
                IF (Nrun.gt.1) THEN
                  rate=100.0_r8*ABS(FOURDVAR(ng)%BackCost(0)+           &
     &                              FOURDVAR(ng)%ObsCost(0)-            &
     &                              FOURDVAR(ng)%CostFun(0))/           &
     &                              FOURDVAR(ng)%CostFun(0)
                ELSE
                  rate=0.0_r8
                END IF
                Optimality(ng)=2.0_r8*(FOURDVAR(ng)%BackCost(0)+        &
     &                                 FOURDVAR(ng)%ObsCost(0))/        &
     &                                (FOURDVAR(ng)%ObsCount(0)-        &
     &                                 FOURDVAR(ng)%ObsReject(0))
                WRITE (stdout,30) outer, inner,                         &
     &                            FOURDVAR(ng)%BackCost(0)/             &
     &                            FOURDVAR(ng)%CostNorm(0),             &
     &                            FOURDVAR(ng)%ObsCost(0)/              &
     &                            FOURDVAR(ng)%CostNorm(0),             &
     &                            rate
                WRITE (stdout,40) outer, inner, Optimality(ng)
!
                DO i=1,NstateVar(ng)
                  IF (FOURDVAR(ng)%ObsCost(i).ne.0.0_r8) THEN
                    IF (Nrun.gt.1) THEN
                      rate=100.0_r8*ABS(FOURDVAR(ng)%BackCost(i)+       &
     &                                  FOURDVAR(ng)%ObsCost(i)-        &
     &                                  FOURDVAR(ng)%CostFun(i))/       &
     &                                  FOURDVAR(ng)%CostFun(i)
                    ELSE
                      rate=0.0_r8
                    END IF
                    IF (i.eq.1) THEN
                      WRITE (stdout,50) outer, inner,                   &
     &                                  FOURDVAR(ng)%BackCost(i)/       &
     &                                  FOURDVAR(ng)%CostNorm(i),       &
     &                                  FOURDVAR(ng)%ObsCost(i)/        &
     &                                  FOURDVAR(ng)%CostNorm(i),       &
     &                                  TRIM(Vname(1,idSvar(i))),       &
     &                                  rate
                    ELSE
                      WRITE (stdout,60) outer, inner,                   &
     &                                  FOURDVAR(ng)%BackCost(i)/       &
     &                                  FOURDVAR(ng)%CostNorm(i),       &
     &                                  FOURDVAR(ng)%ObsCost(i)/        &
     &                                  FOURDVAR(ng)%CostNorm(i),       &
     &                                  TRIM(Vname(1,idSvar(i))),       &
     &                                  rate
                    END IF
                  END IF
                END DO
              END IF
!
!  Compute current total cost function.
!
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
                FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%BackCost(i)+       &
     &                                  FOURDVAR(ng)%ObsCost(i)
              END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time step adjoint model backwards: compute cost function gradients.
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
                END DO
              END DO
!$OMP END PARALLEL DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Descent algorithm: first pass.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  For efficiency, reinitialize decent algorithm every NiterSD
!  iterations to avoid searching in directions which have been
!  seached previously. A steep descent is performed (BetaK=0).
!
              IF (MOD(inner-IterSD,NiterSD).eq.0) CGstepF=CGstepI
              IF (Master) THEN
                WRITE (stdout,70) Nrun, Ipass, CGstepF
              END IF
!
!  Read in previous increment vector deltaV(Lold), Iold record in 
!  ITLname. Also read in old and new total gradient vector
!  GRADv{J(Lold)} and GRADv{J(Lnew)} from adjoint history file
!  ADJname.
!
              IF (inner.gt.1) THEN
                IF (LcycleADJ(ng)) THEN
                  AdjRec=3-tADJindx(ng)
                ELSE
                  AdjRec=tADJindx(ng)-1
                END IF
                CALL get_state (ng, iADM, 3, ADJname(ng), AdjRec,       &
     &                          Lold(ng))
                IF (exit_flag.ne.NoError) RETURN
              END IF
              CALL get_state (ng, iTLM, 8, ITLname(ng), Iold, Lold(ng))
              IF (exit_flag.ne.NoError) RETURN
              CALL get_state (ng, iADM, 4, ADJname(ng), tADJindx(ng),   &
     &                        Lnew(ng))
              IF (exit_flag.ne.NoError) RETURN

#ifdef BALANCE_OPERATOR
              CALL get_state (ng, iNLM, 9, INIname(ng), Lini, Lini)
              IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Convert observation cost function gradient, GRADx(Jo), from model
!  space (x-space) to minimization space (v-space):
!
!     GRADv(Jo) = B^(T/2) GRADx(Jo) = S G L^(T/2) W^(-1/2) 
!
!  First, multiply the adjoint solution, GRADx(Jo), by the background-
!  error standard deviations, S.  Second, convolve result with the
!  adjoint diffusion operator, G L^(T/2) W^(-1/2). Then, compute total
!  cost function gradient (v-space) by adding background and observation
!  contributions:
!
!     GRADv(J) = GRADv(Jb) + GRADv(Jo) = deltaV + GRADv(Jo)
!
!  Compute total cost function gradient norm as:
!
!     CostGrad = SQRT ( transpose{GRADv(J)} * GRADv(J) )
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1
#ifdef BALANCE_OPERATOR
                  CALL ad_balance (ng, TILE, Lini, Lnew(ng))
#endif
                  CALL ad_variability (ng, TILE, Lnew(ng), Lweak)
                  CALL ad_convolution (ng, TILE, Lnew(ng), 2)
                  CALL cost_grad (ng, TILE, Lold(ng), Lnew(ng))
                  CALL cost_norm (ng, TILE, Lnew(ng))
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  Report total cost function gradient norm.
!
              IF (Master) THEN        
                IF (Nrun.gt.1) THEN
                  rate=100.0_r8*ABS(FOURDVAR(ng)%CostGrad(0)-           &
     &                              FOURDVAR(ng)%CostGradOld(0))/       &
     &                          FOURDVAR(ng)%CostGradOld(0)
                ELSE
                  rate=0.0_r8
                END IF
                WRITE (stdout,80) outer, inner,                         &
     &                            FOURDVAR(ng)%CostGrad(0), rate
                DO i=1,NstateVar(ng)
                  IF (FOURDVAR(ng)%CostGrad(i).gt.0.0_r8) THEN
                    IF (Nrun.gt.1) THEN
                      rate=100.0_r8*ABS(FOURDVAR(ng)%CostGrad(i)-       &
     &                                  FOURDVAR(ng)%CostGradOld(i))/   &
     &                              FOURDVAR(ng)%CostGradOld(i)
                    ELSE
                      rate=0.0_r8
                    END IF
                    IF (i.eq.1) THEN
                      WRITE (stdout,90) outer, inner,                   &
     &                                  FOURDVAR(ng)%CostGrad(i),       &
     &                                  TRIM(Vname(1,idSvar(i))), rate
                    ELSE
                      WRITE (stdout,100) outer, inner,                  &
     &                                   FOURDVAR(ng)%CostGrad(i),      &
     &                                   TRIM(Vname(1,idSvar(i))), rate
                    END IF
                  END IF
                END DO
              END IF
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostGradOld(i)=FOURDVAR(ng)%CostGrad(i)
              END DO
!
!  Save total cost function gradient (v-space) into ADJname history
!  NetCDF file. Noticed that the lastest adjoint solution record is
!  over-written in the NetCDF file for future use. The switch
!  "LwrtState2d" is activated to write out state arrays instead
!  ad_*_sol arrays.
!
              kstp(ng)=Lnew(ng)
#ifdef SOLVE3D
              nstp(ng)=Lnew(ng)
#endif
              tADJindx(ng)=tADJindx(ng)-1
              LwrtState2d(ng)=.TRUE.
              CALL ad_wrt_his (ng)
              IF (exit_flag.ne.NoError) RETURN
              LwrtState2d(ng)=.FALSE.
!
!  Estimate new increment vector, deltaV(Lnew), using a CONJUGATE
!  GRADIENT algorithm with first guess step size (CGstepF).
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,CGstepF,numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1
                  CALL descent (ng, TILE, iTLM, inner, CGstepF)
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  Write out trial increment vector, deltaV(Lnew), into record "Inew" of
!  ITLname NetCDF file. It is needed to compute the backgound step size
!  contribution in Ipass=2.
!
              CALL tl_wrt_ini (ng, Lnew(ng), Inew)
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
!     deltaX = B^(1/2) deltaV - (Xi(outer) - Xb)
!  or
!     deltaX = W^(1/2) L^(1/2) G S  - (Xi(outer) - Xb)
!
!  First, convolve estimated increment vector (v-space) by with the
!  tangent linear diffusion operator, W^(1/2) L^(1/2) G.  Second,
!  multiply result by the background-error standard deviation, S.
!  Then, substract current nonlinear initial conditions from the
!  background state.
! 
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1,+1
                  CALL tl_convolution (ng, TILE, Lnew(ng), 2)
                  CALL tl_variability (ng, TILE, Lnew(ng), Lweak)
#ifdef BALANCE_OPERATOR
                  CALL tl_balance (ng, TILE, Lini, Lnew(ng))
#endif
                  CALL tl_ini_adjust (ng, TILE, Lbck, Lini, Lnew(ng))
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  Write out estimated tangent linear initial conditions, deltaX(Lnew), 
!  into record "Litl" of ITLname NetCDF file. It is needed to initialize
!  the tangent linear model.
!
              CALL tl_wrt_ini (ng, Lnew(ng), Litl)
              IF (exit_flag.ne.NoError) RETURN
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time step tangent linear model to determine optimum step size due
!  to observations.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize tangent linear model with just estimated initial
!  conditions, deltaX(Lnew).
!             
              Ipass=2
              CALL tl_initial (ng)
              IF (exit_flag.ne.NoError) RETURN

#ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, turn off history output during
!  second pass.
!
              LdefTLM(ng)=.FALSE.
              LwrtTLM(ng)=.FALSE.
#endif
!
!  Time-step tangent linear model:  Compute optimum step size norms
!  due to observations.
!
              IF (Master) THEN
                WRITE (stdout,20) 'TL', ntstart(ng), ntend(ng)
              END IF

              time(ng)=time(ng)-dt(ng)

              TL_LOOP2 : DO my_iic=ntstart(ng),ntend(ng)+1

                iic(ng)=my_iic
#ifdef SOLVE3D
                CALL tl_main3d (ng)
#else
                CALL tl_main2d (ng)
#endif
                IF (exit_flag.ne.NoError) RETURN

              END DO TL_LOOP2

#ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, close current NetCDF file.
!
              IF (ncTLMid(ng).ne.-1) THEN
                status=nf90_close(ncTLMid(ng))
                ncTLMid(ng)=-1
              END IF
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Descent algorithm: second pass.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read in previous (Iold) and first guess (Inew) increment vectors, 
!  deltaV(Lold) and deltaV(Lnew) respectively.
!
              CALL get_state (ng, iTLM, 8, ITLname(ng), Iold, Lold(ng))
              IF (exit_flag.ne.NoError) RETURN
              CALL get_state (ng, iTLM, 8, ITLname(ng), Inew, Lnew(ng))
              IF (exit_flag.ne.NoError) RETURN
!
!  Compute background contribution (v-space) to the conjugate gradient
!  step size using old and new increments (deltaV) estimates.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1
                  CALL back_step (ng, TILE, Lold(ng), Lnew(ng))
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  Use Derber method (1985) to compute the optimum step size (CGstepR)
!  for the conjugate gradient algorithm.
!
              CGstepR=-CGstepF*((StepTopBck+StepTopObs)/                &
     &                          (StepBotBck+StepBotObs))
              IF (Master) THEN
                WRITE (stdout,70) Nrun, Ipass, CGstepR
              END IF
!
!  Compute new increment vector, deltaV(Lnew), using optimun step size,
!  CGstepR. Then, compute backgound cost function for current inner
!  loop.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,CGstepR,numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1
                  CALL descent (ng, TILE, iTLM, inner, CGstepR)
                  CALL back_cost (ng, TILE, Lnew(ng))
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  Write out new increment vector, deltaV(Lnew), into record "Inew" of
!  ITLname NetCDF file. It is needed by descent in the next minimization
!  iteration. Notice "Inew" record is over-written.
!
              CALL tl_wrt_ini (ng, Lnew(ng), Inew)
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
!  Convert optimum increment vector, deltaV, from minimization space
!  (v-space) to model space (x-space):
!
!     deltaX = B^(1/2) deltaV - (Xi(outer) - Xb)
!  or
!     deltaX = W^(1/2) L^(1/2) G S  - (Xi(outer) - Xb)
!
!  First, convolve estimated increment vector (v-space) by with the
!  tangent linear diffusion operator, W^(1/2) L^(1/2) G.  Second,
!  multiply result by the background-error standard deviation, S.
!  Then, substract current nonlinear initial conditions from the
!  background state.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,CGstepF,numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1,+1
                  CALL tl_convolution (ng, TILE, Lnew(ng), 2)
                  CALL tl_variability (ng, TILE, Lnew(ng), Lweak)
#ifdef BALANCE_OPERATOR
                  CALL tl_balance (ng, TILE, Lini, Lnew(ng))
#endif
                  CALL tl_ini_adjust (ng, TILE, Lbck, Lini, Lnew(ng))
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  Write out optimum tangent linear initial conditions, deltaX(Lnew), 
!  into record "Litl" of ITLname NetCDF file. It is needed to initialize
!  the tangent linear model in the next iteration.
!
              CALL tl_wrt_ini (ng, Lnew(ng), Litl)
              IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Update counters.
!-----------------------------------------------------------------------
!
              Lsav=Lnew(ng)
              Lnew(ng)=Lold(ng)
              Lold(ng)=Lsav
              Isav=Inew
              Inew=Iold
              Iold=Isav
              Nrun=Nrun+1
              CGstepF=CGstepR

            END IF CONVERGE

          END DO INNER_LOOP
!
!-----------------------------------------------------------------------
!  Clear nonlinear state variables.
!-----------------------------------------------------------------------
!
          Ipass=0
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
!         Xi(outer+1) = Xi(outer) + deltaX(Lold)
!
!-----------------------------------------------------------------------
!
!  Notice that "ini_fields" is called here for output purposes only. 
!  It computes the vertically integrated momentum in 3D applications.
!  In order to use the correct fields, the model time indices are set
!  to Lini.
!
          kstp(ng)=Lini
#ifdef SOLVE3D
          nstp(ng)=Lini
#endif
          CALL get_state (ng, iNLM, 1, INIname(ng), Lini, Lini)
          IF (exit_flag.ne.NoError) RETURN

!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL ini_adjust (ng, TILE, Lold(ng), Lini)
              CALL ini_fields (ng, TILE, iNLM)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  Write out new nonlinear model initial conditions.
!
          IF (LcycleINI(ng)) THEN
            tINIindx(ng)=0
            NrecINI(ng)=1
          END IF
          CALL wrt_ini (ng, Lini)
          IF (exit_flag.ne.NoError) RETURN

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Set index containing the surface forcing increments used the run
!  the nonlinear model in the outer loop.
!
          Lfinp(ng)=Lold(ng)
#endif
        END DO OUTER_LOOP
!
!-----------------------------------------------------------------------
!  Done with data assimilation. Initialize the nonlinear model with 
!  estimated initial conditions. Save nonlinear solution at observation
!  points for posterior analysis.
!-----------------------------------------------------------------------
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
!  Compute and report model-observation comparison statistics.
!
        CALL stats_modobs (ng)

      END DO NEST_LOOP
!
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 30   FORMAT (/,'>(',i3.3,',',i3.3,'): Cost Jb, Jo = ',                 &
     &        1p,e16.10,0p,1x,1p,e16.10,0p,t68,1p,e10.4,' %')
 40   FORMAT (1x,'(',i3.3,',',i3.3,'): Optimality  = ',                 &
     &        1p,e16.10)
 50   FORMAT ('<(',i3.3,',',i3.3,'): cost Jb, Jo = ',                   &
     &        1p,e16.10,0p,1x,1p,e16.10,0p,1x,a,t68,1p,e10.4,' %')
 60   FORMAT (1x,'(',i3.3,',',i3.3,'): cost Jb, Jo = ',                 &
     &        1p,e16.10,0p,1x,1p,e16.10,0p,1x,a,t68,1p,e10.4,' %')
 70   FORMAT (/, ' <<<< Descent Algorithm >>>>',/,                      &
     &        6x,'Iteration = ',i5.5,', Ipass = ',i1,                   &
     &        ', Step Size = ',1p,e15.8,/)
 80   FORMAT (/,'>(',i3.3,',',i3.3,'): Gradient Norm = ',               &
     &        1p,e16.10,0p,1x,t68,1p,e10.4,' %')
 90   FORMAT ('<(',i3.3,',',i3.3,'): gradient norm = ',                 &
     &        1p,e16.10,0p,1x,a,t68,1p,e10.4,' %')
100   FORMAT (1x,'(',i3.3,',',i3.3,'): gradient norm = ',               &
     &        1p,e16.10,0p,1x,a,t68,1p,e10.4,' %')
110   FORMAT (a,'_',i3.3,'.nc')

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
