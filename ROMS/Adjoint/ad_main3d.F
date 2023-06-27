#include "cppdefs.h"
#if defined ADJOINT && defined SOLVE3D
      SUBROUTINE ad_main3d (RunInterval)
!
!git $Id$
!svn $Id: ad_main3d.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine is the main driver for adjoint  ROMS/TOMS when      !
!  configurated as a full 3D baroclinic  ocean model. It advances      !
!  backwards the adjoint model equations for all nested grids, if      !
!  any, by the specified time interval (seconds), RunInterval.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined MODEL_COUPLING && defined MCT_LIB
      USE mod_coupler
# endif
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
# ifdef SENSITIVITY_4DVAR
      USE mod_ncparam
# endif
# ifdef NESTING
      USE mod_nesting
# endif
      USE mod_scalars
      USE mod_stepping
# ifdef SO_SEMI
      USE mod_storage
# endif
!
# if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR
      USE adsen_force_mod,         ONLY : adsen_force
# endif
# ifdef ANA_VMIX
      USE analytical_mod,          ONLY : ana_vmix
# endif
# ifdef BIOLOGY
      USE ad_biology_mod,          ONLY : ad_biology
# endif
# ifdef BBL_MODEL_NOT_YET
!!    USE ad_bbl_mod,              ONLY : ad_bblm
# endif
# if defined BULK_FLUXES_NOT_YET && !defined PRIOR_BULK_FLUXES
!!    USE ad_bulk_flux_mod,        ONLY : ad_bulk_flux
# endif
# ifdef BVF_MIXING_NOT_YET
!!    USE ad_bvf_mix_mod,          ONLY : ad_bvf_mix
# endif
      USE ad_diag_mod, ONLY : ad_diag
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE ad_frc_adjust_mod,       ONLY : ad_frc_adjust
# endif
      USE ad_ini_fields_mod,       ONLY : ad_ini_fields, ad_ini_zeta
# ifdef AD_OUTPUT_STATE
      USE ad_ini_fields_mod,       ONLY : ad_out_fields, ad_out_zeta
# endif
# ifdef WEAK_CONSTRAINT
      USE ad_force_dual_mod,       ONLY : ad_force_dual
# endif
# ifdef FORCING_SV
      USE ad_forcing_mod,          ONLY : ad_forcing
# endif
# ifdef GLS_MIXING_NOT_YET
!!    USE ad_gls_corstep_mod,      ONLY : ad_gls_corstep
!!    USE ad_gls_prestep_mod,      ONLY : ad_gls_prestep
# endif
# ifdef LMD_MIXING_NOT_YET
!!    USE ad_lmd_vmix_mod,         ONLY : ad_lmd_vmix
# endif
# if defined FOUR_DVAR && defined OBSERVATIONS
#  ifdef WEAK_CONSTRAINT
      USE ad_htobs_mod,            ONLY : ad_htobs
#  else
      USE ad_misfit_mod,           ONLY : ad_misfit
#  endif
# endif
# ifdef MY25_MIXING
!!    USE ad_my25_corstep_mod,     ONLY : ad_my25_corstep
!!    USE ad_my25_prestep_mod,     ONLY : ad_my25_prestep
# endif
# ifdef NESTING
      USE nesting_mod,             ONLY : nesting
      USE ad_nesting_mod,          ONLY : ad_nesting
#  ifndef ONE_WAY
      USE nesting_mod,             ONLY : do_twoway
#  endif
# endif
# ifdef ADJUST_BOUNDARY
      USE ad_obc_adjust_mod,       ONLY : ad_obc_adjust
      USE ad_obc_adjust_mod,       ONLY : ad_obc2d_adjust
      USE ad_set_depth_mod,        ONLY : ad_set_depth_bry
# endif
      USE ad_omega_mod,            ONLY : ad_omega
# ifdef NEARSHORE_MELLOR_NOT_YET
!!    USE ad_radiation_stress_mod, ONLY : ad_radiation_stress
# endif
# ifndef TS_FIXED
      USE ad_rho_eos_mod,          ONLY : ad_rho_eos
# endif
      USE ad_rhs3d_mod,            ONLY : ad_rhs3d
# ifdef SEDIMENT_NOT_YET
!!    USE ad_sediment_mod,         ONLY : ad_sediment
# endif
# ifdef AD_AVERAGES
      USE ad_set_avg_mod,          ONLY : ad_set_avg
# endif
      USE ad_set_depth_mod,        ONLY : ad_set_depth
      USE ad_set_massflux_mod,     ONLY : ad_set_massflux
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE ad_set_tides_mod,        ONLY : ad_set_tides
# endif
      USE ad_set_vbc_mod,          ONLY : ad_set_vbc
      USE ad_set_zeta_mod,         ONLY : ad_set_zeta
      USE ad_step2d_mod,           ONLY : ad_step2d
# ifndef TS_FIXED
      USE ad_step3d_t_mod,         ONLY : ad_step3d_t
# endif
      USE ad_step3d_uv_mod,        ONLY : ad_step3d_uv
# ifdef FLOATS_NOT_YET
!!    USE ad_step_floats_mod,      ONLY : ad_step_floats
# endif
# if defined BULK_FLUXES && !defined PRIOR_BULK_FLUXES
      USE bulk_flux_mod,           ONLY : bulk_flux
# endif
      USE dateclock_mod,           ONLY : time_string
# ifdef TIDE_GENERATING_FORCES
      USE equilibrium_tide_mod,    ONLY : equilibrium_tide
# endif
# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod,       ONLY : atmos_coupling
# endif
# ifdef WEAK_CONSTRAINT
      USE frc_weak_mod,            ONLY : frc_ADgather, frc_clear
# endif
# if (defined FOUR_DVAR    && !defined I4DVAR_ANA_SENSITIVITY) && \
      defined OBSERVATIONS
      USE obs_read_mod,            ONLY : obs_read
# endif
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod,       ONLY : waves_coupling
# endif
      USE omega_mod,               ONLY : omega
# ifdef SO_SEMI
      USE packing_mod,             ONLY : ad_pack
# endif
      USE rho_eos_mod,             ONLY : rho_eos
      USE set_depth_mod,           ONLY : set_depth
      USE set_massflux_mod
      USE strings_mod,             ONLY : FoundError
# ifdef SENSITIVITY_4DVAR
      USE ad_wvelocity_mod,        ONLY : ad_wvelocity
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
      logical :: backward = .TRUE.
      logical :: ad_advance
      logical :: DoNestLayer, Time_Step
      integer :: Nsteps, Rsteps
      integer :: ig, il, is, istep, ng, nl, tile
      integer :: my_iif, next_indx1
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
      integer :: ks, kt
# ifdef NESTING
      integer :: icount, itcount
# endif
!
      real(r8) :: HalfDT, MaxDT, my_StepTime
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Time-step adjoint 3D primitive equations backwards.
!=======================================================================

# ifdef NESTING
!
!  Compute nonlinear model timestepping sequence and load it into
!  StepInfo array.  The sequence of values are reversed to timestep
!  the adjoint model backwards.
!
      CALL nlm_step_sequence (RunInterval, icount)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
# endif
!
!  Initialize.
!
      Time_Step=.TRUE.
      DoNestLayer=.TRUE.
# ifdef NESTING
      itcount=icount+1
# endif
!
      KERNEL_LOOP : DO WHILE (Time_Step)
!
!  In nesting applications, the number of nesting layers (NestLayers) is
!  used to facilitate refinement grids and composite/refinament grids
!  combinations. Otherwise, the solution it is looped once for a single
!  grid application (NestLayers = 1).
!
        nl=0
# ifdef NESTING
        DO ng=1,Ngrids
          TwoWayInterval(ng)=dt(1)
        END DO
# endif
!
        NEST_LAYER : DO WHILE (DoNestLayer)

# ifdef NESTING
!
!  Determine current nested layer (nl), number of timesteps (Nsteps)
!  to execute for grids in current nested layer, number of timesteps
!  (Rsteps) to complete the RunInterval time window, and timestep
!  counter (step_counter) for each grid. Their values are reversed
!  from the saved nonlinear model time stepping sequence.
!
          itcount=itcount-1             ! count backwards for adjoint
          IF (itcount.gt.0) THEN
            nl=StepInfo(itcount,1)
            Nsteps=StepInfo(itcount,2)
            IF (itcount.eq.icount) Nsteps=1
            Rsteps=StepInfo(itcount,3)
            DO ng=1,Ngrids
              step_counter(ng)=StepInfo(itcount,ng+3)
            END DO
          ELSE
            nl=0
          END IF
# else
!
!  Determine number of time steps to compute in each nested grid layer
!  based on the specified time interval (seconds), RunInterval. Non
!  nesting applications have NestLayers=1. Notice that RunInterval is
!  set in the calling driver. Its value may span the full period of the
!  simulation, a multi-model coupling interval (RunInterval > ifac*dt),
!  or just a single step (RunInterval=0).
!
          CALL ntimesteps (iADM, RunInterval, nl, Nsteps, Rsteps)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
          IF ((nl.le.0).or.(nl.gt.NestLayers)) EXIT
!
!  Time-step governing equations for Nsteps.
!
          STEP_LOOP : DO istep=Nsteps,1,-1
!
!  Set time indices and time clock.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
# ifdef JEDI
              jic(ng)=jic(ng)-1
              time4jedi(ng)=time4jedi(ng)-dt(ng)
# endif
              tdays(ng)=time(ng)*sec2day
              IF (step_counter(ng).eq.1) Time_Step=.FALSE.
              IF ((iic(ng).eq.ntstart(ng)).and.(ng.eq.Ngrids)) THEN
                ad_advance=.FALSE.
# ifdef NESTING
                DO ng=1,Ngrids
                  TwoWayInterval(ng)=0.0_r8
                END DO
# endif
              ELSE
                ad_advance=.TRUE.
              END IF
            END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL ad_get_data (ng)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!-----------------------------------------------------------------------
!  Process input data, if any: time interpolate between snapshots.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL ad_set_data (ng, tile)
# ifdef FORWARD_READ
                CALL set_depth (ng, tile, iADM)
# endif
# ifdef TIDE_GENERATING_FORCES
                CALL equilibrium_tide (ng, tile, iADM)
# endif
              END DO
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FORWARD_READ
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE horizontal mass fluxes (Hz*u/n and Hz*v/m).
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL set_massflux (ng, tile, iADM)
                CALL rho_eos (ng, tile, iADM)
#  if defined BULK_FLUXES && !defined PRIOR_BULK_FLUXES
!!              CALL bulk_flux (ng, tile)
#  endif
              END DO
            END DO
# endif
!
!  Avoid time-stepping if additional delayed IO time-step.
!
            ADVANCE : IF (ad_advance) THEN

# ifdef FLOATS_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories: Split all the drifters
!  between all the computational threads, except in distributed-memory
!  and serial configurations. In distributed-memory, the parallel node
!  containing the drifter is selected internally since the state
!  variables do not have a global scope.
!-----------------------------------------------------------------------
!
!  Shift floats time indices.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (Lfloats(Ng)) THEN
                  nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
                  nf  (ng)=MOD(nf  (ng)+1,NFT+1)
                  nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
                  nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
                  nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
!
#  ifdef _OPENMP
                  chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
                  Lstr=1+MyThread*chunk_size
                  Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
#  else
                  Lstr=1
                  Lend=Nfloats(ng)
#  endif
                  CALL ad_step_floats (ng, Lstr, Lend)
                END IF
              END DO
# endif

# ifdef NESTING
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
                  CALL ad_nesting (ng, iADM, ngetD)
                END IF
              END DO
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
              DO il=1,NestLayers
                DO ig=1,GridsInLayer(il)
                  ng=GridNumber(ig,il)
                  IF (do_twoway(iADM,nl,il,ng,istep)) THEN
                    CALL ad_nesting (ng, iADM, n2way)
                  END IF
                END DO
              END DO
#  endif
# endif

# ifndef TS_FIXED
!
!-----------------------------------------------------------------------
!  Time-step adjoint tracer equations.
!-----------------------------------------------------------------------
!
#  ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for Tracer Variables.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (ANY(CompositeGrid(:,ng))) THEN
                  CALL ad_nesting (ng, iADM, n3dTV)
                END IF
              END DO
#  endif
!  Compute intermediate BASIC STATE mass fluxes (Huon,Hvom) for use in
!  the adjoint horizontal advection (ad_step3d_t) and adjoint vertical
!  velocity (ad_omega).
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL reset_massflux (ng, tile, iADM)    ! intermediate
                END DO                                    ! fluxes
              END DO
!
!  Compute basic STATE omega vertical velocity with intermediate mass
!  fluxes. Time-step adjoint tracer equations.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL omega (ng, tile, iADM)   ! BASIC STATE w-velocity
                  CALL ad_step3d_t (ng, tile)
                END DO
              END DO
# endif
!
!-----------------------------------------------------------------------
!  Time-step adjoint vertical mixing turbulent equations and passive
!  tracer source and sink terms, if applicable. Reinstate original
!  BASIC STATE (Huon,Hvom).
!-----------------------------------------------------------------------
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                DO tile=first_tile(ng),last_tile(ng),+1
# ifdef SEDIMENT_NOT_YET
                  CALL ad_sediment (ng, tile)
# endif
# ifdef BIOLOGY
                  CALL ad_biology (ng, tile)
# endif
# ifdef MY25_MIXING_NOT_YET
                  CALL ad_my25_corstep (ng, tile)
# elif defined GLS_MIXING_NOT_YET
                  CALL ad_gls_corstep (ng, tile)
# endif
                  CALL ad_omega (ng, tile, iADM)
                  CALL set_massflux (ng, tile, iADM)     ! BASIC STATE
                END DO                                   ! mass fluxes
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
                  CALL ad_nesting (ng, iADM, n3duv)
                END IF
              END DO
# endif
!
!-----------------------------------------------------------------------
!  Time-step adjoint 3D equations.
!-----------------------------------------------------------------------
!
!  Reinstate original BASIC STATE omega vertical velocity. Time-step
!  3D adjoint momentum equations and couple with vertically integrated
!  equations.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL omega (ng, tile, iADM)   ! BASIC STATE w-velocity
                  CALL ad_step3d_uv (ng, tile)
                END DO
              END DO
!
!-----------------------------------------------------------------------
!  Adjoint of recompute depths and thicknesses using the new time
!  filtered free-surface. This call was moved from "ad_step2d" to here.
!-----------------------------------------------------------------------
!
# ifdef NESTING
!
!  If nesting, determine vertical indices and vertical interpolation
!  weights in the contact zone using new depth arrays.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                CALL ad_nesting (ng, iADM, nzwgt)
              END DO
# endif
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ad_set_depth (ng, tile, iADM)
                END DO
              END DO

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the time-averaged
!  momentum fluxes (DU_avg1, DV_avg1) and free-surface (Zt_avg).
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                IF (ANY(CompositeGrid(:,ng))) THEN
                  CALL ad_nesting (ng, iADM, n2dfx)
                END IF
              END DO

#  if defined MASKING && defined WET_DRY
!
!  If nesting and wetting and drying, scale horizontal interpolation
!  weights to account for land/sea masking in contact areas. This needs
!  to be done at very time-step since the Land/Sea masking is time
!  dependent.
!
              DO ig=1,GridsInLayer(nl)
                ng=GridNumber(ig,nl)
                CALL ad_nesting (ng, iADM, nmask)
              END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Solve adjoint vertically integrated primitive equations for the
!  free-surface and barotropic momentum components.
!-----------------------------------------------------------------------
!
              LOOP_2D : DO my_iif=MAXVAL(nfast)+1,1,-1

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine Corrector step section.
!  If refinement, check mass flux conservation between coarse and
!  fine grids during debugging.
!
                DO ig=1,GridsInLayer(nl)
                  ng=GridNumber(ig,nl)
                  IF (ANY(CompositeGrid(:,ng))) THEN
                    CALL ad_nesting (ng, iADM, n2dCS)
                  END IF
                  IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                    CALL ad_nesting (ng, iADM, nmflx)
                  END IF
                END DO
# endif
!
!  Corrector step - Apply 2D adjoint time-step corrector scheme.  Notice
!  ==============    that there is not need for a corrector step during
!  the auxiliary (nfast+1) time-step.
!
                DO ig=1,GridsInLayer(nl)
                  ng=GridNumber(ig,nl)
                  iif(ng)=my_iif
                  IF (iif(ng).lt.(nfast(ng)+1)) THEN
                    DO tile=first_tile(ng),last_tile(ng),+1
                      CALL ad_step2d (ng, tile)
                    END DO
                  END IF
!
!  Set time indices for adjoint predictor step.
!
                  next_indx1=3-indx1(ng)
                  IF (.not.PREDICTOR_2D_STEP(ng)) THEN
                    PREDICTOR_2D_STEP(ng)=.TRUE.
!^                  knew(ng)=next_indx1
!^                  kstp(ng)=3-knew(ng)
!^                  krhs(ng)=3
!^
                    kt=knew(ng)
                    ks=kstp(ng)
                    knew(ng)=krhs(ng)
                    kstp(ng)=kt
                    krhs(ng)=ks
!^                  IF (my_iif.lt.(nfast(ng)+1)) indx1(ng)=next_indx1
                  END IF
                END DO

# ifdef NESTING

!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for the state variables
!  associated with the 2D engine Predictor Step section.
!  If refinement, check mass flux conservation between coarse and
!  fine grids during debugging.
!
                DO ig=1,GridsInLayer(nl)
                  ng=GridNumber(ig,nl)
                  IF (ANY(CompositeGrid(:,ng))) THEN
                    CALL ad_nesting (ng, iADM, n2dPS)
                  END IF
                  IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                    CALL ad_nesting (ng, iADM, nmflx)
                  END IF
                END DO
# endif
!
!  Predictor step - Advance adjoint barotropic equations using 2D
!  ==============   time-step predictor scheme.  No actual time-
!  stepping is performed during the auxiliary (nfast+1) time-step.
!  It is needed to finalize the fast-time averaging of 2D fields,
!  if any, and compute the new time-evolving depths.
!
                DO ig=1,GridsInLayer(nl)
                  ng=GridNumber(ig,nl)
                  IF (my_iif.le.(nfast(ng)+1)) THEN
                    DO tile=last_tile(ng),first_tile(ng),-1
                      CALL ad_step2d (ng, tile)
                    END DO
                  END IF
!
!  Set time indices for next adjoint corrector step. The
!  PREDICTOR_2D_STEP switch it is assumed to be false before the
!  first time-step.
!
                  IF (PREDICTOR_2D_STEP(ng).and.                        &
     &                my_iif.le.(nfast(ng)+1)) THEN
                    PREDICTOR_2D_STEP(ng)=.FALSE.
!^                  IF (FIRST_2D_STEP) THEN
!^                    kstp(ng)=indx1(ng)
!^                  ELSE
!^                    kstp(ng)=3-indx1(ng)
!^                  END IF
!^                  knew(ng)=3
!^                  krhs(ng)=indx1(ng)
!^
                    ks=knew(ng)
                    knew(ng)=krhs(ng)
                    krhs(ng)=ks
                  END IF
                END DO

              END DO LOOP_2D

            END IF ADVANCE

# if (defined FOUR_DVAR    && !defined I4DVAR_ANA_SENSITIVITY) && \
      defined OBSERVATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, read observation and model state at observation
!  locations.  Then, compute adjoint forcing terms due to observations.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
#  ifdef SENSITIVITY_4DVAR
#   ifdef RBL4DVAR_FCT_SENSITIVITY
              IF (.not.Lsen4DVAR(ng).and.LsenFCT(ng)) THEN
#   else
              IF (.not.Lsen4DVAR(ng)) THEN
#   endif
#  endif
              HalfDT=0.5_r8*dt(ng)
              IF (((time(ng)-HalfDT).le.ObsTime(ng)).and.               &
     &             (ObsTime(ng).lt.(time(ng)+HalfDT))) THEN
                ProcessObs(ng)=.TRUE.
                CALL obs_read (ng, iADM, backward)
              ELSE
                ProcessObs(ng)=.FALSE.
              END IF
!
#  ifdef SP4DVAR
!
!  Skip assimilation of obs on first timestep unless inter=Nsaddle.
!
              IF ((iic(ng).ne.ntstart(ng)).or.Lsadd(ng)) THEN
#  endif

              DO tile=first_tile(ng),last_tile(ng),+1
#  ifdef WEAK_CONSTRAINT
                CALL ad_htobs (ng, tile, iADM)
#  else
                CALL ad_misfit (ng, tile, iADM)
#  endif
              END DO
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
#  ifdef SENSITIVITY_4DVAR
              END IF
#  endif
#  ifdef SP4DVAR
              END IF
#  endif
            END DO
# endif

# ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  If appropriate, add representer coefficients (Beta hat) impulse
!  forcing to adjoint solution. Read next impulse record, if available.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ProcessObs(ng)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_force_dual (ng, tile, kstp(ng), nstp(ng))
                END DO
              END IF
            END DO
# endif

# ifdef NESTING

!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for right-hand-side terms
!  (tracers).
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL ad_nesting (ng, iADM, nrhst)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute adjoint right-hand-side terms for 3D equations. If
!  applicable, time-step turbulent mixing schemes.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
# ifdef MY25_MIXING_NOT_YET
                CALL ad_my25_prestep (ng, tile)
# elif defined GLS_MIXING_NOT_YET
                CALL ad_gls_prestep (ng, tile)
# endif
                CALL ad_rhs3d (ng, tile)
              END DO
            END DO

# ifdef NESTING
!
!-----------------------------------------------------------------------
!  If refinement grid, interpolate (space, time) state variables
!  contact points from donor grid extracted data.
!
!  Also, fill BRY_CONTACT(:,:)%Mflux to check for mass conservation
!  between coarse and fine grids.  This is only done for diagnostic
!  purposes.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                CALL ad_nesting (ng, iADM, nmflx)
                CALL ad_nesting (ng, iADM, nputD)
              END IF
            END DO
# endif
!!          END IF ADVANCE

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for 3D kernel free-surface.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL ad_nesting (ng, iADM, nzeta)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Set adjoint free-surface to it time-averaged value.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
# ifdef DIAGNOSTICS
!!              CALL set_diags (ng, tile)
# endif
                CALL ad_set_zeta (ng, tile)
              END DO
            END DO
!
!-----------------------------------------------------------------------
!  Compute adjoint vertical mixing coefficients for momentum and
!  tracers. Compute adjoint S-coordinate vertical velocity,
!  diagnostically from horizontal mass divergence.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
# ifdef SENSITIVITY_4DVAR
                IF (SCALARS(ng)%Lstate(isWvel)) THEN
                  CALL ad_wvelocity (ng, tile, nstp(ng))
                END IF
# endif
                CALL ad_omega (ng, tile, iADM)
# if defined ANA_VMIX_NOT_YET
                CALL ad_ana_vmix (ng, tile, iADM)
# elif defined LMD_MIXING_NOT_YET
                CALL ad_lmd_vmix (ng, tile)
# elif defined BVF_MIXING
                CALL ad_bvf_mix (ng, tile)
# endif
              END DO
            END DO

# ifdef SO_SEMI
!
!-----------------------------------------------------------------------
!  If stochastic optimals with respect the seminorm of chosen
!  functional, pack adjoint state surface forcing needed by the
!  dynamical propagator.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (MOD(iic(ng)-1,nADJ(ng)).eq.0) THEN
                SOrec(ng)=SOrec(ng)+1
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_pack (ng, tile, Nstr(ng), Nend(ng),           &
     &                          STORAGE(ng)%so_state(:,SOrec(ng)))
                END DO
              END IF
            END DO
# endif

# ifdef NESTING
!
!  If composite or mosaic grids, process additional points in the
!  contact zone between connected grids for bottom stress variables.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (ANY(CompositeGrid(:,ng))) THEN
                CALL ad_nesting (ng, iADM, nbstr)
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Set adjoint fields for vertical boundary conditions. Process tidal
!  forcing, if any.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
                CALL ad_set_tides (ng, tile)
# endif
                CALL ad_set_vbc (ng, tile)
# ifdef BBL_MODEL_NOT_YET
                CALL ad_bblm (ng, tile)
# endif
# if defined BULK_FLUXES_NOT_YET && !defined PRIOR_BULK_FLUXES
                CALL ad_bulk_flux (ng, tile)
# endif
              END DO
            END DO

# ifdef NEARSHORE_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL ad_radiation_stress (ng, tile)
              END DO
            END DO
# endif

# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple to waves model every CoupleSteps(Iwaves) timesteps: get
!  waves/sea fluxes.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF ((iic(ng).ne.ntstart(ng)).and.                         &
     &            MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL waves_coupling (ng, tile)
                END DO
              END IF
            END DO
# endif

# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple to atmospheric model every CoupleSteps(Iatmos) timesteps: get
!  air/sea fluxes.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF ((iic(ng).ne.ntstart(ng)).and.                         &
     &            MOD(iic(ng)-1,CoupleSteps(Iatmos,ng)).eq.0) THEN
                DO tile=last_tile(ng),first_tile(ng),-1
                      CALL atmos_coupling (ng, tile)
                END DO
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute adjoint density related fields and horizontal mass fluxes
!  (Hz*u/n and Hz*v/m). If applicable, compute and report diagnostics
!  and accumulate time-averaged adjoint fields which needs a
!  irreversible loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1     ! irreversible
# ifndef TS_FIXED
                CALL ad_rho_eos (ng, tile, iADM)
# endif
                CALL ad_set_massflux (ng, tile, iADM)
                CALL ad_diag (ng, tile)
# ifdef AD_AVERAGES
                CALL ad_set_avg (ng, tile)
# endif
              END DO
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FORCING_SV
!
!-----------------------------------------------------------------------
!  Compute the adjoint forcing for the forcing singular vectors.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL ad_set_depth (ng, tile, iADM)
                CALL ad_forcing (ng, tile, kstp(ng), nstp(ng))
              END DO
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Update 3D time-level indices.
!-----------------------------------------------------------------------
!
!  The original forward time-stepping indices are advanced as follows:
!
!     nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
!     nnew(ng)=3-nstp(ng)
!     nrhs(ng)=nstp(ng)
!
!  This yields the only 2 permutations:
!
!     nstp  nnew  nrhs
!      1     2     1
!      2     1     2
!
!   With nstp=1, nnew=1 and nrhs=2 at time zero, this is equivalent to
!   the following:
!
!     nstp(ng)=nnew(ng)
!     nnew(ng)=nrhs(ng)
!     nrhs(ng)=nstp(ng)
!
!   The adjoint of this is as follows:
!
!     nstp(ng)=nrhs(ng)
!     nrhs(ng)=nnew(ng)
!     nnew(ng)=nstp(ng)
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).ne.ntend(ng)) THEN
                nrhs(ng)=nnew(ng)
                nnew(ng)=nstp(ng)
                nstp(ng)=nrhs(ng)
              END IF
            END DO

# ifdef AD_OUTPUT_STATE
!
!-----------------------------------------------------------------------
!  Set full adjoint output arrays. Due to the exact discrete adjoint,
!  the predictor/corrector time-stepping scheme with multiple time
!  levels, pieces of the adjoint solution are in two-time levels and
!  need to be added in the "_sol" arrays for output purposes.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).ne.ntend(ng)) THEN
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ad_out_fields (ng, tile, iADM)
                END DO
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_set_depth (ng, tile, iADM)
                  CALL ad_out_zeta (ng, tile, iADM)
                END DO
              END IF
            END DO
# endif
# ifndef FORCING_SV
!
!-----------------------------------------------------------------------
!  If not a restart, initialize all time levels and compute other
!  initial fields.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).eq.ntend(ng)) THEN
#  ifdef NESTING
!
!  Extract donor grid initial data at contact points and store it in
!  REFINED structure so it can be used for the space-time interpolation.
!
                IF (RefinedGrid(ng)) THEN
                  CALL ad_nesting (ng, iADM, ngetD)
                END IF
#  endif
!
!  Adjoint of initialize other state variables.
!
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ad_ini_fields (ng, tile, iADM)
                END DO
!
!  Adjoint of initialize free-surface and compute initial level
!  thicknesses and depths.
!
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_set_depth (ng, tile, iADM)
                  CALL ad_ini_zeta (ng, tile, iADM)
                END DO
              END IF
            END DO
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Skip first timestep.
!-----------------------------------------------------------------------
!
#  ifdef RBL4DVAR_FCT_SENSITIVITY
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (.not.Lsen4DVAR(ng)) THEN          ! ignore in forecast
                IF (iic(ng).ne.ntstart(ng)) THEN    !           interval
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL ad_frc_adjust (ng, tile, Lfout(ng))
                  END DO
                END IF
              END IF
            END DO
#  else
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).ne.ntstart(ng)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_frc_adjust (ng, tile, Lfout(ng))
                END DO
              END IF
            END DO
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundaries.
!  Skip first timestep.
!-----------------------------------------------------------------------
!
#  ifdef RBL4DVAR_FCT_SENSITIVITY
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (.not.Lsen4DVAR(ng)) THEN          ! ignore in forecast
                IF (iic(ng).ne.ntstart(ng)) THEN    !           interval
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL ad_obc2d_adjust (ng, tile, Lbout(ng))
                    CALL ad_set_depth_bry (ng, tile, iADM)
                    CALL ad_obc_adjust (ng, tile, Lbout(ng))
                  END DO
                END IF
              END IF
            END DO
#  else
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).ne.ntstart(ng)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_obc2d_adjust (ng, tile, Lbout(ng))
                  CALL ad_set_depth_bry (ng, tile, iADM)
                  CALL ad_obc_adjust (ng, tile, Lbout(ng))
                END DO
              END IF
            END DO
#  endif
# endif

# if defined WEAK_CONSTRAINT && !defined SP4DVAR
!
!-----------------------------------------------------------------------
!  Gather weak constraint forcing to storage arrays.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).ne.ntstart(ng)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ad_set_depth (ng, tile, iADM)
                  CALL frc_ADgather (ng, tile)
                END DO
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              CALL ad_output (ng)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO

# if defined WEAK_CONSTRAINT && !defined SP4DVAR
!
!-----------------------------------------------------------------------
!  Copy storage arrays index 1 into index 2, and clear index 1.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (MOD(iic(ng)-1,nADJ(ng)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL frc_clear (ng, tile)
                END DO
              END IF
            END DO
# endif

# if (defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR)     && \
     !defined OBS_SPACE
!
!-----------------------------------------------------------------------
!  Add appropriate forcing terms to the adjoint model. The form of the
!  forcing depends on the functional whose sensitivity is required.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
#  ifdef SENSITIVITY_4DVAR
              IF (Lsen4DVAR(ng)) THEN
#  endif
#  if !defined AD_IMPULSE
                IF ((DendS(ng).ge.tdays(ng)).and.                       &
     &              (tdays(ng).ge.DstrS(ng))) THEN
#  endif
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL adsen_force (ng, tile)
                  END DO
#  if !defined AD_IMPULSE
                END IF
#  endif
#  ifdef SENSITIVITY_4DVAR
              END IF
#  endif
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Set adjoint time indices and time clock.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iic(ng)=iic(ng)-1
# ifndef NESTING
              step_counter(ng)=step_counter(ng)-1
# endif
              time(ng)=time(ng)-dt(ng)
              CALL time_string (time(ng), time_code(ng))
            END DO

          END DO STEP_LOOP

        END DO NEST_LAYER

      END DO KERNEL_LOOP
!
      RETURN
      END SUBROUTINE ad_main3d
#else
      SUBROUTINE ad_main3d
      RETURN
      END SUBROUTINE ad_main3d
#endif
