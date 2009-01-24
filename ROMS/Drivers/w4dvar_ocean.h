      MODULE ocean_control_mod
!
!svn $Id: w4dvar_ocean.h 652 2008-07-24 23:20:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group   Emanuele Di Lorenzo   !
!    Licensed under a MIT/X style license            Andrew M. Moore   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Weak Constraint 4-Dimensional Variational (4DVar)         !
!            Data Assimilation Driver: Indirect Representer Approach   !
!                                                                      !
!  This driver is used for weak constraint 4DVar where errors are      !
!  considered in both model and observations.                          !
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
!  Read in model-error standard deviation factors and spatial
!  convolution diffusion convolution.
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
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
#ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE impulse_mod, ONLY : impulse
      USE ini_adjust_mod, ONLY : rp_ini_adjust
      USE ini_adjust_mod, ONLY : load_ADtoTL
      USE ini_adjust_mod, ONLY : load_TLtoAD
      USE mod_ocean, ONLY : initialize_ocean
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
      logical :: add, converged, outer_impulse
      logical :: Lweak

      integer :: ADrec, Lbck, Lini, Nrec, Rec1, Rec2
      integer :: i, lstr, my_iic, ng, rec, status, subs, tile, thread

      real(r8) :: MyTime, LB_time, UB_time
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!  Initialize relevant parameters.
!
        Lold(ng)=1          ! old minimization time index
        Lnew(ng)=2          ! new minimization time index
        Lini=1              ! NLM initial conditions record in INIname
        Lbck=2              ! background record in INIname
        Rec1=1
        Rec2=2 
        Nrun=1
        Ipass=1
        outer=0
        inner=0
        ERstr=1
        ERend=Nouter
!
!-----------------------------------------------------------------------
!  Configure weak constraint 4DVAR algorithm: Indirect Representer
!  Approach.
!-----------------------------------------------------------------------
!
!  Initialize and set nonlinear model initial conditions.
!
        wrtNLmod(ng)=.TRUE.
        wrtRPmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
        CALL initial (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!  Save nonlinear initial conditions (currently in time index 1,
!  background) into record "Lini" of INIname NetCDF file. The record
!  "Lbck" becomes the background state record and the record "Lini"
!  becomes current nonlinear initial conditions.
!
        IF (LcycleINI(ng)) THEN
          tINIindx(ng)=1
          NrecINI(ng)=1
        END IF
        CALL wrt_ini (ng, 1)
        IF (exit_flag.ne.NoError) RETURN
!
!  Set nonlinear output history file as the initial basic state
!  trajectory.
!
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        lstr=LEN_TRIM(FWDbase(ng))
        WRITE (HISname(ng),10) FWDbase(ng)(1:lstr-3), outer
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Model-error covariance normalization and stardard deviation factors.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Compute or read in the model-error correlation normalization factors.
!  If computing, write out factors to NetCDF. This is an expensive
!  computation that needs to be computed only once for a particular
!  application grid and decorrelation scales.
!  
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
!
!  Define TLM/RPM impulse forcing NetCDF file
!
        LdefTLF(ng)=.TRUE.
        CALL def_impulse (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define output 4DVAR NetCDF file containing all processed data
!  at observation locations.
!
        LdefMOD(ng)=.TRUE.
        CALL def_mod (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire IDs of tangent linear and representer model initial
!  conditions NetCDF files.
!
        LdefITL(ng)=.FALSE.
        LdefIRP(ng)=.FALSE.
        CALL tl_def_ini (ng)
        IF (exit_flag.ne.NoError) RETURN
        CALL rp_def_ini (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Run nonlinear model and compute basic state trajectory.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
        IF (Master) THEN
          WRITE (stdout,20) 'NL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        NL_LOOP : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL main3d (ng)
#else
          CALL main2d (ng)
#endif
          IF (exit_flag.ne.NoError) RETURN

        END DO NL_LOOP
        wrtNLmod(ng)=.FALSE.
!
!  Set forward basic state NetCDF ID to nonlinear model trajectory to
!  avoid the inquiring stage. 
!
        ncFWDid(ng)=ncHISid(ng)
!
!-----------------------------------------------------------------------
!  Solve the system:
!
!              (R_n + Cobs) * Beta_n = h_n
!
!              h_n = Xo - H * X_n
!
!  where R_n is the representer matrix, Cobs is the observation-error
!  covariance, Beta_n are the representer coefficients, h_n is the
!  misfit between observations (Xo) and model (H * X_n), and H is
!  the linearized observation operator. Here, _n denotes a sequence
!  of estimates.
!
!  The system does not need to be solved explicitly by inverting the
!  symmetric stabilized representer matrix, P_n:
!
!              P_n = R_n + Cobs
!
!  but by computing the action of P_n on any vector PSI, such that
!
!              P_n * PSI = R_n * PSI + Cobs * PSI
!
!  The representer matrix is not explicitly computed but evaluated by
!  one integration backward of the adjoint model and one integration
!  forward of the tangent linear model for any forcing vector PSI.
!
!  A preconditioned conjugate gradient algorithm is used to compute
!  an approximation PSI for Beta_n.
!
!-----------------------------------------------------------------------
!
        OUTER_LOOP : DO outer=1,Nouter
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Run representer model and compute a "prior estimate" state
!  trajectory, X_n(t). Use linearized state trajectory (X_n-1) as
!  basic state.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set representer model basic state trajectory file to previous outer
!  loop file (outer-1). If outer=1, the basic state trajectory is the
!  nonlinear model. 
!
          lstr=LEN_TRIM(FWDbase(ng))
          WRITE (FWDname(ng),10) FWDbase(ng)(1:lstr-3), outer-1
!
!  Activate switch to write the representer model at observation points.
!  Turn off writing into history file and turn off impulse forcing.
!
          wrtRPmod(ng)=.TRUE.
          LdefTLM(ng)=.FALSE.
          LwrtTLM(ng)=.FALSE.
          SporadicImpulse=.FALSE.
          FrequentImpulse=.FALSE.
!
!  As in the nonlinear model, initialize always the representer model
!  here with the background or reference state (IRPname, record Rec1).
!
          tIRPindx(ng)=Rec1
          CALL rp_initial (ng)
          IF (exit_flag.ne.NoError) RETURN
!
!  Run representer model using the nonlinear trajectory as a basic
!  state.  Compute model solution at observation points, H * X_n.
!
          IF (Master) THEN
            WRITE (stdout,20) 'RP', ntstart(ng), ntend(ng)
          END IF

          time(ng)=time(ng)-dt(ng)

          RP_LOOP1 : DO my_iic=ntstart(ng),ntend(ng)+1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL rp_main3d (ng)
#else
            CALL rp_main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO RP_LOOP1
          wrtRPmod(ng)=.FALSE.
!
!  Set approximation vector PSI to representer coefficients Beta_n.
!  Here, PSI is set to misfit between observations and model, H_n.
!
          inner=0
          CALL congrad (ng, outer, inner, Ninner, converged)

          INNER_LOOP : DO inner=1,Ninner
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate adjoint model forced with any vector PSI at the observation
!  locations and generate adjoint trajectory, Lambda_n(t).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize the adjoint model from rest.
!
            CALL ad_initial (ng)
            IF (exit_flag.ne.NoError) RETURN
            wrtMisfit(ng)=.FALSE.
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file only once to avoid opening too many files.
!
            IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
            NrecADJ(ng)=0
            tADJindx(ng)=0
!
!  Time-step adjoint model backwards forced with current PSI vector.
!
            IF (Master) THEN
              WRITE (stdout,20) 'AD', ntstart(ng), ntend(ng)
            END IF

            time(ng)=time(ng)+dt(ng)

            AD_LOOP1 : DO my_iic=ntstart(ng),ntend(ng),-1

              iic(ng)=my_iic
#ifdef SOLVE3D
              CALL ad_main3d (ng)
#else
              CALL ad_main2d (ng)
#endif
              IF (exit_flag.ne.NoError) RETURN

            END DO AD_LOOP1

#ifdef CONVOLVE
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Convolve adjoint trajectory with error covariances and convert
!  to impulse forcing.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
            Nrec=NrecADJ(ng)
            NrecADJ(ng)=0
            tADJindx(ng)=0
            LwrtState2d(ng)=.TRUE.
            LwrtTime(ng)=.FALSE.
            IF (Master) THEN
              WRITE (stdout,30) outer, inner
            END IF
!
!  Clear adjoint state arrays.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1
                CALL initialize_ocean (ng, TILE, iADM)
              END DO
            END DO
!$OMF END PARALLEL DO
!
!  Convolve initial conditions record (ADJname, record Nrec) with
!  initial conditions background error covariance. Since routine
!  "get_state" loads data into the ghost points, the adjoint
!  solution is read into the tangent linear state arrays by using
!  iTLM instead of iADM in the calling arguments.
!
            ADrec=Nrec
            Lweak=.FALSE.
            CALL get_state (ng, iTLM, 4, ADJname(ng), ADrec, Lold(ng))
            IF (exit_flag.ne.NoError) RETURN
!
!  Load interior solution, read above, into adjoint state arrays. 
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the 
!  squared-root adjoint diffusion operator which impose initial
!  conditions background error covaraince. Notice that the spatial
!  convolution is only done for half of the diffusion steps
!  (squared-root filter). Clear tangent linear state arrays when
!  done.
!
            add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1
                CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
# ifdef BALANCE_OPERATOR_NOT_YET
                CALL ad_balance (ng, TILE, Lini, Lold(ng))
# endif
                CALL ad_variability (ng, TILE, Lold(ng), Lweak)
                CALL ad_convolution (ng, TILE, Lold(ng), 2)
                CALL initialize_ocean (ng, TILE, iTLM)
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  background-error standard deviations. Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points.
!
            add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL load_ADtoTL (ng, TILE, Lold(ng), Lold(ng), add)
                CALL tl_convolution (ng, TILE, Lold(ng), 2)
                CALL tl_variability (ng, TILE, Lold(ng), Lweak)
# ifdef BALANCE_OPERATOR_NOT_YET
                CALL tl_balance (ng, TILE, Lini, Lold(ng))
# endif
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  Write out tangent linear model initial conditions for next inner
!  loop into ITLname (record Rec1). The tangent model initial
!  conditions are set to the convolved adjoint solution.
!
            CALL tl_wrt_ini (ng, Lold(ng), Rec1) 
            IF (exit_flag.ne.NoError) RETURN
!
!  If weak constraint, convolve all adjoint records in ADJname and
!  impose model error covariance.
!
            IF (Nrec.gt.1) THEN
              DO rec=1,Nrec
                Lweak=.TRUE.
!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
                ADrec=rec
                CALL get_state (ng, iTLM, 4, ADJname(ng), ADrec,        &
     &                          Lold(ng))
                IF (exit_flag.ne.NoError) RETURN
!
!  Load interior solution, read above, into adjoint state arrays. 
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the 
!  squared-root adjoint diffusion operator which impose the model-error
!  spatial correlations. Notice that the spatial convolution is only
!  done for half of the diffusion steps (squared-root filter). Clear
!  tangent linear state arrays when done.
!
                add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
                DO thread=0,numthreads-1
                  subs=NtileX(ng)*NtileE(ng)/numthreads
                  DO tile=subs*thread,subs*(thread+1)-1
                    CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
# ifdef BALANCE_OPERATOR_NOT_YET
                    CALL ad_balance (ng, TILE, Lini, Lold(ng))
# endif
                    CALL ad_variability (ng, TILE, Lold(ng), Lweak)
                    CALL ad_convolution (ng, TILE, Lold(ng), 2)
                    CALL initialize_ocean (ng, TILE, iTLM)
                  END DO
                END DO
!$OMP END PARALLEL DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  background-error standard deviations.  Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points. Copy back to adjoint state
!  arrays when done with the convolution for output purposes.
!
                add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
                DO thread=0,numthreads-1
                  subs=NtileX(ng)*NtileE(ng)/numthreads
                  DO tile=subs*thread,subs*(thread+1)-1,+1
                    CALL load_ADtoTL (ng, TILE, Lold(ng), Lold(ng), add)
                    CALL tl_convolution (ng, TILE, Lold(ng), 2)
                    CALL tl_variability (ng, TILE, Lold(ng), Lweak)
# ifdef BALANCE_OPERATOR_NOT_YET
                    CALL tl_balance (ng, TILE, Lini, Lold(ng))
# endif
                    CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
                  END DO
                END DO
!$OMP END PARALLEL DO
!
!  Overwrite ADJname history NetCDF file with convolved adjoint
!  solution.
!
                kstp(ng)=Lold(ng)
# ifdef SOLVE3D
                nstp(ng)=Lold(ng)
# endif
                CALL ad_wrt_his (ng)
                IF (exit_flag.ne.NoError) RETURN
              END DO
              LwrtState2d(ng)=.FALSE.
              LwrtTime(ng)=.TRUE.
            END IF
#endif
!
!  Convert the current adjoint solution in ADJname to impulse forcing.
!  Write out impulse forcing into TLFname NetCDF file. To facilitate
!  the forcing to the TLM and RPM, the forcing is processed and written
!  in increasing time coordinates (recall that the adjoint solution
!  in ADJname is backwards in time).
!
            IF (Master) THEN
              WRITE (stdout,40) outer, inner
            END IF
            tTLFindx(ng)=0
            outer_impulse=.FALSE.
#ifdef DISTRIBUTE
            tile=MyRank
#else
            tile=-1
#endif
            CALL impulse (ng, tile, iADM, outer_impulse, ADJname(ng))
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory (impulse forcing) to compute R_n * PSI at observation
!  points.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
            TLMname(ng)=TLMbase(ng)
            wrtNLmod(ng)=.FALSE.
            wrtTLmod(ng)=.TRUE.
!
!  Initialize tangent linear model from ITLname, record Rec1. 
!
            tITLindx(ng)=Rec1
            CALL tl_initial (ng)
            IF (exit_flag.ne.NoError) RETURN
!
!  Activate switch to write out initial misfit between model and
!  observations.
!
            IF ((outer.eq.1).and.(inner.eq.1)) THEN
              wrtMisfit(ng)=.TRUE.
            END IF
!
!  Set tangent linear history NetCDF parameters.  Define tangent linear
!  history file at the beggining of each inner loop  to avoid opening
!  too many NetCDF files.
!
            IF (inner.gt.1) LdefTLM(ng)=.FALSE.
            NrecTLM(ng)=0
            tTLMindx(ng)=0
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute R_n * PSI at observation points which
!  are used in the conjugate gradient algorithm.
!
            IF (Master) THEN
              WRITE (stdout,20) 'TL', ntstart(ng), ntend(ng)
            END IF

            MyTime=time(ng)
            time(ng)=time(ng)-dt(ng)

            TL_LOOP : DO my_iic=ntstart(ng),ntend(ng)+1

              iic(ng)=my_iic
!
!  Set impulse forcing switches.  Notice that at this point, time(ng)
!  has not been incremented, so for the following check we need to
!  add dt(ng).
!
              LB_time=time(ng)+0.5_r8*dt(ng)
              UB_time=time(ng)+1.5_r8*dt(ng)
              SporadicImpulse=(LB_time.le.FrcTime(ng)).and.             &
     &                        (FrcTime(ng).lt.UB_time).and.             &
     &                        (NrecFrc(ng).gt.1)
              FrequentImpulse=.FALSE.
#ifdef SOLVE3D
              CALL tl_main3d (ng)
#else
              CALL tl_main2d (ng)
#endif
              MyTime=time(ng)

              IF (exit_flag.ne.NoError) RETURN

            END DO TL_LOOP
            wrtNLmod(ng)=.FALSE.
            wrtTLmod(ng)=.FALSE.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Use conjugate gradient algorithm to find a better approximation
!  PSI to representer coefficients Beta_n. Exit inner loop if
!  convergence is achieved.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
            Nrun=Nrun+1
            CALL congrad (ng, outer, inner, Ninner, converged)
            IF (converged) EXIT INNER_LOOP

          END DO INNER_LOOP
!
!  Close tangent linear NetCDF file.
!
          status=nf90_close(ncTLMid(ng))
          ncTLMid(ng)=-1
!
!-----------------------------------------------------------------------
!  Once that the representer coefficients, Beta_n, have been
!  approximated with sufficient accuracy, compute estimates of
!  Lambda_n and Xhat_n by carrying out one backward intergration
!  of the adjoint model and one forward itegration of the representer
!  model.
!-----------------------------------------------------------------------
!
!  Initialize the adjoint model always from rest.
!
          CALL ad_initial (ng)
          IF (exit_flag.ne.NoError) RETURN
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
          IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
          NrecADJ(ng)=0
          tADJindx(ng)=0
!
!  Time-step adjoint model backwards forced with estimated representer
!  coefficients, Beta_n.
!
          IF (Master) THEN
            WRITE (stdout,20) 'AD', ntstart(ng), ntend(ng)
          END IF

          time(ng)=time(ng)+dt(ng)

          AD_LOOP2 : DO my_iic=ntstart(ng),ntend(ng),-1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL ad_main3d (ng)
#else
            CALL ad_main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO AD_LOOP2

#ifdef CONVOLVE
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Convolve adjoint trajectory with model-error covariance and convert
!  to impulse forcing.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
          Nrec=NrecADJ(ng)
          NrecADJ(ng)=0
          tADJindx(ng)=0
          LwrtState2d(ng)=.TRUE.
          LwrtTime(ng)=.FALSE.
          IF (Master) THEN
            WRITE (stdout,30) outer, inner
          END IF
!
!  Clear adjoint state arrays.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL initialize_ocean (ng, TILE, iADM)
            END DO
          END DO
!$OMF END PARALLEL DO
!
!  Convolve initial conditions record (ADJname, record Nrec) with
!  initial conditions background error covariance. Since routine
!  "get_state" loads data into the ghost points, the adjoint
!  solution is read into the tangent linear state arrays by using
!  iTLM instead of iADM in the calling arguments.
!
          ADrec=Nrec
          Lweak=.FALSE.
          CALL get_state (ng, iTLM, 4, ADJname(ng), ADrec, Lold(ng))
          IF (exit_flag.ne.NoError) RETURN
!
!  Load interior solution, read above, into adjoint state arrays. 
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the 
!  squared-root adjoint diffusion operator which impose initial
!  conditions background error covaraince. Notice that the spatial
!  convolution is only done for half of the diffusion steps
!  (squared-root filter). Clear tangent linear state arrays when
!  done.
!
          add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
# ifdef BALANCE_OPERATOR_NOT_YET
              CALL ad_balance (ng, TILE, Lini, Lold(ng))
# endif
              CALL ad_variability (ng, TILE, Lold(ng), Lweak)
              CALL ad_convolution (ng, TILE, Lold(ng), 2)
              CALL initialize_ocean (ng, TILE, iTLM)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  background-error standard deviations.  Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points. Copy back to adjoint state
!  arrays when done with the convolution. Compute representer model
!  initial conditions by adding convolved adjoint solution to the
!  reference nonlinear state (INIname, record Lbck).
!
          CALL get_state (ng, iNLM, 9, INIname(ng), Lbck, Lnew(ng))
          IF (exit_flag.ne.NoError) RETURN

          add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL load_ADtoTL (ng, TILE, Lold(ng), Lold(ng), add)
              CALL tl_convolution (ng, TILE, Lold(ng), 2)
              CALL tl_variability (ng, TILE, Lold(ng), Lweak)
# ifdef BALANCE_OPERATOR_NOT_YET
              CALL tl_balance (ng, TILE, Lini, Lold(ng))
# endif
              CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
              CALL rp_ini_adjust (ng, TILE, Lold(ng), Lbck)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  Write out representer model initial conditions into IRPname, record
!  Rec2.
!
          CALL rp_wrt_ini (ng, Lold(ng), Rec2) 
          IF (exit_flag.ne.NoError) RETURN
!
!  If weak constraint, convolve adjoint records in ADJname and impose
!  model error covariance.
!
          IF (Nrec.gt.1) THEN
            DO rec=1,Nrec
              Lweak=.TRUE.
!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
              ADrec=rec
              CALL get_state (ng, iTLM, 4, ADJname(ng), ADrec, Lold(ng))
              IF (exit_flag.ne.NoError) RETURN
!
!  Load interior solution, read above, into adjoint state arrays. 
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the 
!  squared-root adjoint diffusion operator which impose the model-error
!  spatial correlations. Notice that the spatial convolution is only
!  done for half of the diffusion steps (squared-root filter). Clear
!  tangent linear state arrays when done.
!
!  WARNING: We need to add logic here for new decorrelation scales
!           and normalization coefficients.
!  
              add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1
                  CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
# ifdef BALANCE_OPERATOR_NOT_YET
                  CALL ad_balance (ng, TILE, Lini, Lold(ng))
# endif
                  CALL ad_variability (ng, TILE, Lold(ng), Lweak)
                  CALL ad_convolution (ng, TILE, Lold(ng), 2)
                  CALL initialize_ocean (ng, TILE, iTLM)
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  background-error standard deviations.  Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points. Copy back to adjoint state
!  arrays when done with the convolution for output purposes.
!
              add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
              DO thread=0,numthreads-1
                subs=NtileX(ng)*NtileE(ng)/numthreads
                DO tile=subs*thread,subs*(thread+1)-1,+1
                  CALL load_ADtoTL (ng, TILE, Lold(ng), Lold(ng), add)
                  CALL tl_convolution (ng, TILE, Lold(ng), 2)
                  CALL tl_variability (ng, TILE, Lold(ng), Lweak)
# ifdef BALANCE_OPERATOR_NOT_YET
                  CALL tl_balance (ng, TILE, Lini, Lold(ng))
# endif
                  CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
                END DO
              END DO
!$OMP END PARALLEL DO
!
!  Overwrite ADJname history NetCDF file with convolved adjoint
!  solution.
!
              kstp(ng)=Lold(ng)
# ifdef SOLVE3D
              nstp(ng)=Lold(ng)
# endif
              CALL ad_wrt_his (ng)
              IF (exit_flag.ne.NoError) RETURN
            END DO
            LwrtState2d(ng)=.FALSE.
            LwrtTime(ng)=.TRUE.
          END IF
#endif
!
!  Convert the current adjoint solution in ADJname to impulse forcing.
!  Write out impulse forcing into TLFname NetCDF file. To facilitate
!  the forcing to the TLM and RPM, the forcing is processed and written
!  in increasing time coordinates (recall that the adjoint solution
!  in ADJname is backwards in time).
!
          IF (Master) THEN
            WRITE (stdout,40) outer, inner
          END IF
          tTLFindx(ng)=0
          outer_impulse=.TRUE.
#ifdef DISTRIBUTE
          tile=MyRank
#else
          tile=-1
#endif
          CALL impulse (ng, tile, iADM, outer_impulse, ADJname(ng))
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Run representer model and compute a "new estimate" of the state
!  trajectory, X_n(t).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set new basic state trajectory for next outer loop.
!
          LdefTLM(ng)=.TRUE.
          LwrtTLM(ng)=.TRUE.
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.TRUE.
          lstr=LEN_TRIM(FWDbase(ng))
          WRITE (TLMname(ng),10) FWDbase(ng)(1:lstr-3), outer
!
!  If weak constraint, the impulses are time-interpolated at each
!  time-steps.
!
          SporadicImpulse=.FALSE.
          IF (FrcRec(ng).gt.1) THEN
            FrequentImpulse=.TRUE. 
          END IF
!
!  Initialize representer model IRPname file, record Rec2.
!
          tIRPindx(ng)=Rec2
          CALL rp_initial (ng)
          IF (exit_flag.ne.NoError) RETURN
!
!  Activate switch to write out final misfit between model and
!  observations.
!
          IF (outer.eq.Nouter) THEN
            wrtMisfit(ng)=.TRUE.
          END IF
!
!  Run representer model using previous linearized trajectory, X_n-1, as
!  basic state and forced with convolved adjoint trajectory impulses.
!
          IF (Master) THEN
            WRITE (stdout,20) 'RP', ntstart(ng), ntend(ng)
          END IF

          time(ng)=time(ng)-dt(ng)

          RP_LOOP2 : DO my_iic=ntstart(ng),ntend(ng)+1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL rp_main3d (ng)
#else
            CALL rp_main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO RP_LOOP2
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.FALSE.
!
!  Close current forward NetCDF file.
!
          status=nf90_close(ncFWDid(ng))
          ncFWDid(ng)=-1

        END DO OUTER_LOOP
!
!  Compute and report model-observation comparison statistics.
!
        CALL stats_modobs (ng)

      END DO NEST_LOOP
!
 10   FORMAT (a,'_',i3.3,'.nc')
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 30   FORMAT (/,' Convolving Adjoint Trajectory: Outer = ',i3.3,        &
     &          ' Inner = ',i3.3)
 40   FORMAT (/,' Converting Convolved Adjoint Trajectory to',          &
     &          ' Impulses: Outer = ',i3.3,' Inner = ',i3.3,/)

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
