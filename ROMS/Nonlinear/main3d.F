#include "cppdefs.h"
#if defined NONLINEAR && defined SOLVE3D
      SUBROUTINE main3d (RunInterval)
!
!svn $Id: main3d.F 1054 2021-03-06 19:47:12Z arango $
!=======================================================================
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as a full 3D baroclinic ocean model.  It  advances     !
!  forward the primitive equations for all  nested  grids, if any,     !
!  for the specified time interval (seconds), RunInterval.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined MODEL_COUPLING && defined MCT_LIB && !defined COAWST_MODEL
      USE mod_coupler
# endif
# ifdef COAWST_COUPLING
      USE mct_coupler_params
# endif
      USE mod_iounits
# ifdef NESTING
      USE mod_nesting
# endif
      USE mod_scalars
      USE mod_stepping
!
# ifdef ANA_VMIX
      USE analytical_mod,       ONLY : ana_vmix
# endif
# ifdef BIOLOGY
      USE biology_mod,          ONLY : biology
# endif
# ifdef BBL_MODEL
      USE bbl_mod,              ONLY : bblm
# endif
# ifdef BULK_FLUXES
      USE bulk_flux_mod,        ONLY : bulk_flux
# endif
# ifdef ATM2OCN_FLUXES
      USE atm2ocn_flux_mod,     ONLY : atm2ocn_flux
# endif
# ifdef BVF_MIXING
      USE bvf_mix_mod,          ONLY : bvf_mix
# endif
      USE dateclock_mod,        ONLY : time_string
      USE diag_mod,             ONLY : diag
# ifdef TLM_CHECK
      USE dotproduct_mod,       ONLY : nl_dotproduct
# endif
# if defined NLM_OUTER                || \
     defined RBL4DVAR                 || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined SP4DVAR
      USE forcing_mod,          ONLY : forcing
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE frc_adjust_mod,       ONLY : frc_adjust, load_frc
# endif
# ifdef GLS_MIXING
      USE gls_corstep_mod,      ONLY : gls_corstep
      USE gls_prestep_mod,      ONLY : gls_prestep
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
      USE hmixing_mod,          ONLY : hmixing
# endif
      USE ini_fields_mod,       ONLY : ini_fields, ini_zeta
#ifdef INWAVE_MODEL
      USE driver_inwave_mod,    ONLY : inwave_run
      USE mod_inwave_swan
#endif
# ifdef LMD_MIXING
      USE lmd_vmix_mod,         ONLY : lmd_vmix
# endif
# ifdef MY25_MIXING
      USE my25_corstep_mod,     ONLY : my25_corstep
      USE my25_prestep_mod,     ONLY : my25_prestep
# endif
# ifdef NESTING
      USE nesting_mod,          ONLY : nesting
#  ifndef ONE_WAY
      USE nesting_mod,          ONLY : do_twoway
#  endif
# endif
# if defined ADJUST_BOUNDARY
      USE obc_adjust_mod,       ONLY : obc_adjust, load_obc
# endif
# if defined ROMS_COUPLING && defined MCT_LIB
      USE ocean_coupler_mod,    ONLY : ocean_coupling
# endif
      USE omega_mod,            ONLY : omega
# ifdef WEC_MELLOR
      USE radiation_stress_mod, ONLY : radiation_stress
# endif
# ifdef WEC_VF
      USE wec_vf_mod, ONLY : wec_vf
#  if defined WDISS_THORGUZA || defined WDISS_CHURTHOR
      USE wec_dissip_mod,       ONLY : wec_dissip
#  endif
#  ifdef WEC_ROLLER
      USE wec_roller_mod,       ONLY : wec_roller
#  endif
      USE wec_stokes_mod,       ONLY : wec_stokes
#  ifdef WAVE_MIXING
      USE wec_wave_mix_mod,     ONLY : wec_wave_mix
#  endif
# endif
# ifdef SEDIMENT
      USE sediment_mod,         ONLY : sediment
# endif
# if defined VEGETATION && defined MARSH_DYNAMICS
      USE mod_vegetation 
      USE mod_vegarr
      USE marsh_dynamics_mod, ONLY : marsh_dynamics
# endif
# ifndef TS_FIXED
      USE rho_eos_mod,          ONLY : rho_eos
# endif
      USE rhs3d_mod,            ONLY : rhs3d
# ifdef AVERAGES
      USE set_avg_mod,          ONLY : set_avg
# endif
      USE set_depth_mod,        ONLY : set_depth
      USE set_massflux_mod,     ONLY : set_massflux
# if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
      USE set_tides_mod,        ONLY : set_tides
# endif
      USE set_vbc_mod,          ONLY : set_vbc
      USE set_zeta_mod,         ONLY : set_zeta
      USE step2d_mod,           ONLY : step2d
# ifndef TS_FIXED
      USE step3d_t_mod,         ONLY : step3d_t
# endif
      USE step3d_uv_mod,        ONLY : step3d_uv
# ifdef FLOATS
      USE step_floats_mod,      ONLY : step_floats
# endif
      USE strings_mod,          ONLY : FoundError
      USE wvelocity_mod,        ONLY : wvelocity
# ifdef WEC
      USE wstvelocity_mod,      ONLY : wstvelocity
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: DoNestLayer, Time_Step
!
      integer :: Nsteps, Rsteps
      integer :: ig, iw, ia, il, istep, ng, nl, tile
      integer :: my_iif, next_indx1
# ifdef FLOATS
      integer :: Lend, Lstr, chunk_size
# endif
# ifdef INWAVE_MODEL
      integer :: instep
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Time-step nonlinear 3D primitive equations by the specified time.
!=======================================================================
!
!  Time-step the 3D kernel for the specified time interval (seconds),
!  RunInterval.
!
      Time_Step=.TRUE.
      DoNestLayer=.TRUE.
!
      KERNEL_LOOP : DO WHILE (Time_Step)
!
!  In nesting applications, the number of nesting layers (NestLayers) is
!  used to facilitate refinement grids and composite/refinament grids
!  combinations. Otherwise, the solution it is looped once for a single
!  grid application (NestLayers = 1).
!
        nl=0
#ifdef NESTING
        TwoWayInterval(1:Ngrids)=0.0_r8
#endif
!
        NEST_LAYER : DO WHILE (DoNestLayer)
!
!  Determine number of time steps to compute in each nested grid layer
!  based on the specified time interval (seconds), RunInterval. Non
!  nesting applications have NestLayers=1. Notice that RunInterval is
!  set in the calling driver. Its value may span the full period of the
!  simulation, a multi-model coupling interval (RunInterval > ifac*dt),
!  or just a single step (RunInterval=0).
!
          CALL ntimesteps (iNLM, RunInterval, nl, Nsteps, Rsteps)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF ((nl.le.0).or.(nl.gt.NestLayers)) EXIT
!
!  Time-step governing equations for Nsteps.
!
          STEP_LOOP : DO istep=1,Nsteps
!
!  Set time indices and time clock.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iic(ng)=iic(ng)+1
              nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
              nnew(ng)=3-nstp(ng)
              nrhs(ng)=nstp(ng)
              time(ng)=time(ng)+dt(ng)
              tdays(ng)=time(ng)*sec2day
              CALL time_string (time(ng), time_code(ng))
              IF (step_counter(ng).eq.Rsteps) Time_Step=.FALSE.
            END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
!$OMP MASTER
              CALL get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL set_data (ng, tile)
              END DO
!$OMP BARRIER
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined NLM_OUTER                || \
     defined RBL4DVAR                 || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution impulse forcing to
!  the nonlinear model solution. Notice that the forcing is only needed
!  after finishing all the inner loops. The forcing is continuous.
!  That is, it is time interpolated at every time-step from available
!  snapshots (FrequentImpulse=TRUE).
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
#  if defined WEAK_NOINTERP
              IF ((iic(ng).gt.1).and.(iic(ng).ne.ntend(ng)+1).and.      &
     &            (MOD(iic(ng)-1,nadj(ng)).eq.0)) THEN
                IF (Master) THEN
                  WRITE (stdout,*) ' FORCING NLM at iic = ',iic(ng)
                END IF
#  endif
                IF (FrequentImpulse(ng)) THEN
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL forcing (ng, tile, kstp(ng), nstp(ng))
                    CALL set_depth (ng, tile, iNLM)
                  END DO
!$OMP BARRIER
                END IF
#  if defined WEAK_NOINTERP
              END IF
#  endif
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface and compute initial level thicknesses and
!  depths.
!
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ini_zeta (ng, tile, iNLM)
                  CALL set_depth (ng, tile, iNLM)
                END DO
!$OMP BARRIER
!
!  Initialize other state variables.
!
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ini_fields (ng, tile, iNLM)
                END DO
!$OMP BARRIER

# ifdef NESTING
!
!  Extract donor grid initial data at contact points and store it in
!  REFINED structure so it can be used for the space-time interpolation.
!
                IF (RefinedGrid(ng)) THEN
                  CALL nesting (ng, iNLM, ngetD)
                END IF
# endif
              END IF
            END DO
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL set_massflux (ng, tile, iNLM)
# ifndef TS_FIXED
                CALL rho_eos (ng, tile, iNLM)
# endif
                CALL diag (ng, tile)
# ifdef TLM_CHECK
                CALL nl_dotproduct (ng, tile, Lnew(ng))
# endif
              END DO
!$OMP BARRIER
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined ROMS_COUPLING && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple ocean to atmosphere and wave model grids.
!-----------------------------------------------------------------------
            IF (iic(ng).ne.ntstart(ng)) THEN
              CALL ocean_coupling (nl)
            END IF
# endif
# if defined INWAVE_MODEL
!
!-----------------------------------------------------------------------
!  Couple ocean to infragravity wave model grids.
!-----------------------------------------------------------------------
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO instep=1,INT(1.0_r8/WAVES(ng)%dtifac)
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL inwave_run (ng, tile)
                END DO
              END DO
            END DO
# endif

# ifdef WEC_VF
!
!-----------------------------------------------------------------------
!  Compute wave dissipation, roller, and stokes terms.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
#  if defined WDISS_THORGUZA || defined WDISS_CHURTHOR
                CALL wec_dissip (ng, tile)
#  endif
#  ifdef WEC_ROLLER
                CALL wec_roller (ng, tile)
#  endif
                CALL wec_stokes (ng, tile)
              END DO
!$OMP BARRIER
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef WEC_MELLOR
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL radiation_stress (ng, tile)
              END DO
!$OMP BARRIER
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
# ifdef BULK_FLUXES
#  if defined FOUR_DVAR && defined PRIOR_BULK_FLUXES
                IF (Nrun.eq.1) CALL bulk_flux (ng, tile)
#  else
                CALL bulk_flux (ng, tile)
#  endif
# endif
# ifdef ATM2OCN_FLUXES
                CALL atm2ocn_flux (ng, tile)
# endif
# ifdef BBL_MODEL
                CALL bblm (ng, tile)
# endif
                CALL set_vbc (ng, tile)
# if defined SSH_TIDES || defined UV_TIDES
                CALL set_tides (ng, tile)
# endif
              END DO
!$OMP BARRIER
            END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for bottom stress variables.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, nbstr)
              END IF
            END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundary.
!  Load open boundary into storage arrays. Skip the last output
!  timestep.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).lt.(ntend(ng)+1)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL obc_adjust (ng, tile, Lbinp(ng))
                  CALL load_obc (ng, tile, Lbout(ng))
                END DO
!$OMP BARRIER
              END IF
            END DO
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Load surface forcing into storage arrays. Skip the last output
!  timestep.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).lt.(ntend(ng)+1)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL frc_adjust (ng, tile, Lfinp(ng))
                  CALL load_frc (ng, tile, Lfout(ng))
                END DO
!$OMP BARRIER
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute time-dependent vertical/horizontal mixing coefficients for
!  momentum and tracers. Compute S-coordinate vertical velocity,
!  diagnostically from horizontal mass divergence.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
# if defined ANA_VMIX
                CALL ana_vmix (ng, tile, iNLM)
# elif defined LMD_MIXING
                CALL lmd_vmix (ng, tile)
# elif defined BVF_MIXING
                CALL bvf_mix (ng, tile)
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
                CALL hmixing (ng, tile)
# endif
                CALL omega (ng, tile, iNLM)
                CALL wvelocity (ng, tile, nstp(ng))
# ifdef WEC
                CALL wstvelocity (ng, TILE, nstp(ng))
# endif
              END DO
!$OMP BARRIER
            END DO
!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1     ! irreversible
                CALL set_zeta (ng, tile)
# ifdef DIAGNOSTICS
                CALL set_diags (ng, tile)
# endif
# ifdef AVERAGES
#   if defined FILTERED
                CALL set_filter (ng, tile)
#   else
                CALL set_avg (ng, tile)
#   endif
# endif
              END DO
!$OMP BARRIER
            END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for 3D kernel free-surface.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, nzeta)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
!$OMP MASTER
              CALL output (ng)
!$OMP END MASTER
!$OMP BARRIER
              IF ((FoundError(exit_flag, NoError, __LINE__, MyFile)).or.&
     &            ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) THEN
                RETURN
              END IF
            END DO

# ifdef NESTING
!
!-----------------------------------------------------------------------
!  If refinement grid, interpolate (space, time) state variables
!  contact points from donor grid extracted data.
#  ifdef NESTING_DEBUG
!
!  Also, fill BRY_CONTACT(:,:)%Mflux to check for mass conservation
!  between coarse and fine grids.  This is only done for diagnostic
!  purposes. Also, debugging is possible with very verbose output
!  to fort.300 is allowed by activating uppercase(nesting_debug).
#  endif
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                CALL nesting (ng, iNLM, nputD)
#  ifdef NESTING_DEBUG
                CALL nesting (ng, iNLM, nmflx)
#  endif
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL rhs3d (ng, tile)
# ifdef MY25_MIXING
                CALL my25_prestep (ng, tile)
# elif defined GLS_MIXING
                CALL gls_prestep (ng, tile)
# endif
              END DO
!$OMP BARRIER
            END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for right-hand-side terms
!  (tracers).
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, nrhst)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and barotropic momentum components.
!-----------------------------------------------------------------------
!
            LOOP_2D : DO my_iif=1,MAXVAL(nfast)+1
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                next_indx1=3-indx1(ng)
                IF (.not.PREDICTOR_2D_STEP(ng).and.                     &
     &              my_iif.le.(nfast(ng)+1)) THEN
                  PREDICTOR_2D_STEP(ng)=.TRUE.
                  iif(ng)=my_iif
                  IF (FIRST_2D_STEP) THEN
                    kstp(ng)=indx1(ng)
                  ELSE
                    kstp(ng)=3-indx1(ng)
                  END IF
                  knew(ng)=3
                  krhs(ng)=indx1(ng)
                END IF
!
!  Predictor step - Advance barotropic equations using 2D time-step
!  ==============   predictor scheme.  No actual time-stepping is
!  performed during the auxiliary (nfast+1) time-step. It is needed
!  to finalize the fast-time averaging of 2D fields, if any, and
!  compute the new time-evolving depths.
!
                IF (my_iif.le.(nfast(ng)+1)) THEN
                  DO tile=last_tile(ng),first_tile(ng),-1
                    CALL step2d (ng, tile)
                  END DO
!$OMP BARRIER
                END IF
              END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine Predictor Step section.
!  If refinement, check mass flux conservation between coarse and
!  fine grids during debugging.
#  ifdef NESTING_DEBUG
!  Warning: very verbose output to fort.300 ascii file to check
!           mass flux conservation.
#  endif
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (ANY(CompositeGrid(:,ng))) THEN
                  CALL nesting (ng, iNLM, n2dPS)
                END IF
#  ifdef NESTING_DEBUG
                IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                  CALL nesting (ng, iNLM, nmflx)
                END IF
#  endif
              END DO
# endif
!
!  Set time indices for corrector step.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (PREDICTOR_2D_STEP(ng)) THEN
                  PREDICTOR_2D_STEP(ng)=.FALSE.
                  knew(ng)=next_indx1
                  kstp(ng)=3-knew(ng)
                  krhs(ng)=3
                  IF (iif(ng).lt.(nfast(ng)+1)) indx1(ng)=next_indx1
                END IF
!
!  Corrector step - Apply 2D time-step corrector scheme.  Notice that
!  ==============   there is not need for a corrector step during the
!  auxiliary (nfast+1) time-step.
!
                IF (iif(ng).lt.(nfast(ng)+1)) THEN
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL step2d (ng, tile)
                  END DO
!$OMP BARRIER
                END IF
              END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine Corrector step section.
!  If refinement, check mass flux conservation between coarse and
!  fine grids during debugging.
#  ifdef NESTING_DEBUG
!  Warning: very verbose output to fort.300 ascii file to check
!           mass flux conservation.
#  endif
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (ANY(CompositeGrid(:,ng))) THEN
                  CALL nesting (ng, iNLM, n2dCS)
                END IF
#  ifdef NESTING_DEBUG
                IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                  CALL nesting (ng, iNLM, nmflx)
                END IF
#  endif
              END DO
# endif
            END DO LOOP_2D

# ifdef NESTING
#  if defined MASKING && defined WET_DRY
!
!  If nesting and wetting and drying, scale horizontal interpolation
!  weights to account for land/sea masking in contact areas. This needs
!  to be done at very time-step since the Land/Sea masking is time
!  dependent.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL nesting (ng, iNLM, nmask)
            END DO
#  endif
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the time-averaged
!  momentum fluxes (DU_avg1, DV_avg1) and free-surface (Zt_avg).
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, n2dfx)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Recompute depths and thicknesses using the new time filtered
!  free-surface.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL set_depth (ng, tile, iNLM)
              END DO
!$OMP BARRIER
            END DO

# ifdef NESTING
!
!  If nesting, determine vertical indices and vertical interpolation
!  weights in the contact zone using new depth arrays.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL nesting (ng, iNLM, nzwgt)
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Time-step 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Time-step 3D momentum equations and couple with vertically
!  integrated equations.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL step3d_uv (ng, tile)
              END DO
!$OMP BARRIER
            END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for 3D momentum (u,v),
!  adjusted 2D momentum (ubar,vbar), and fluxes (Huon, Hvom).
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, n3duv)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL omega (ng, tile, iNLM)
# ifdef MY25_MIXING
                CALL my25_corstep (ng, tile)
# elif defined GLS_MIXING
                CALL gls_corstep (ng, tile)
# endif
# if defined WEC_VF && defined WAVE_MIXING
            CALL wec_wave_mix (ng, tile)
# endif
# ifdef BIOLOGY
                CALL biology (ng, tile)
# endif
# ifdef SEDIMENT
                CALL sediment (ng, tile)
# endif
              END DO
!$OMP BARRIER
            END DO

# if defined VEGETATION && defined MARSH_DYNAMICS
!
!-----------------------------------------------------------------------
!  Calculate marsh erosion due to wave thrust. 
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL marsh_dynamics (ng, tile)
              END DO
            END DO
# endif 
# ifndef TS_FIXED
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL step3d_t (ng, tile)
              END DO
!$OMP BARRIER
            END DO

#  ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for Tracer Variables.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL nesting (ng, iNLM, n3dTV)
              END IF
            END DO
#  endif
# endif

# ifdef NESTING
#  ifndef ONE_WAY
!
!-----------------------------------------------------------------------
!  If refinement grids, perform two-way coupling between fine and
!  coarse grids. Correct coarse grid tracers values at the refinement
!  grid with refined accumulated fluxes.  Then, replace coarse grid
!  state variable with averaged refined grid values (two-way nesting).
!  Update coarse grid depth variables.
!
!  The two-way exchange of infomation between nested grids needs to be
!  done at the correct time-step and in the right sequence.
!-----------------------------------------------------------------------
!
            DO il=NestLayers,1,-1
              DO ig=1,GridsInLayer(il)
                ng=GridNumber(ig,il)
                IF (do_twoway(iNLM, nl, il, ng, istep)) THEN
                  CALL nesting (ng, iNLM, n2way)
                END IF
              END DO
            END DO
#  endif
!
!-----------------------------------------------------------------------
!  If donor to a finer grid, extract data for the external contact
!  points. This is the latest solution for the coarser grid.
!
!  It is stored in the REFINED structure so it can be used for the
!  space-time interpolation when "nputD" argument is used above.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (DonorToFiner(ng)) THEN
                CALL nesting (ng, iNLM, ngetD)
              END IF
            END DO
# endif

# ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories: Split all the drifters
!  between all the computational threads, except in distributed-memory
!  and serial configurations. In distributed-memory, the parallel node
!  containing the drifter is selected internally since the state
!  variables do not have a global scope.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (Lfloats(ng)) THEN
#  ifdef _OPENMP
                chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
                Lstr=1+MyThread*chunk_size
                Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
#  else
                Lstr=1
                Lend=Nfloats(ng)
#  endif
                CALL step_floats (ng, Lstr, Lend)
!$OMP BARRIER
!
!  Shift floats time indices.
!
                nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
                nf  (ng)=MOD(nf  (ng)+1,NFT+1)
                nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
                nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
                nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
              END IF
            END DO
# endif

          END DO STEP_LOOP

        END DO NEST_LAYER

      END DO KERNEL_LOOP

      RETURN
      END SUBROUTINE main3d
#else
      SUBROUTINE main3d
      RETURN
      END SUBROUTINE main3d
#endif
