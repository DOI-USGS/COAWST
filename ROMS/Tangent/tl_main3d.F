#include "cppdefs.h"
#if defined TANGENT && defined SOLVE3D
      SUBROUTINE tl_main3d (RunInterval)
!
!git $Id$
!svn $Id: tl_main3d.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the main driver for tangent linear ROMS/TOMS when   !
!  configurated as a full  3D baroclinic  ocean model.  It  advances   !
!  forward the tangent linear model equations  for all nested grids,   !
!  if any, by the specified time interval (seconds), RunInterval.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined MODEL_COUPLING && defined MCT_LIB
      USE mod_coupler
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
      USE dateclock_mod,        ONLY : time_string
# ifdef TLM_CHECK
      USE dotproduct_mod,       ONLY : tl_dotproduct
# endif
# ifdef TIDE_GENERATING_FORCES
      USE equilibrium_tide_mod, ONLY : equilibrium_tide
# endif
# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod,    ONLY : atmos_coupling
# endif
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod,    ONLY : waves_coupling
# endif
# ifdef FORWARD_READ
      USE omega_mod,            ONLY : omega
      USE set_depth_mod,        ONLY : set_depth
      USE set_massflux_mod,     ONLY : set_massflux
# endif
      USE strings_mod,          ONLY : FoundError
# ifdef BIOLOGY
      USE tl_biology_mod,       ONLY : tl_biology
# endif
# ifdef BBL_MODEL_NOT_YET
!!    USE tl_bbl_mod,           ONLY : tl_bblm
# endif
# if defined BULK_FLUXES_NOT_YET && !defined PRIOR_BULK_FLUXES
!!    USE tl_bulk_flux_mod,     ONLY : tl_bulk_flux
# endif
# ifdef BVF_MIXING_NOT_YET
!!    USE tl_bvf_mix_mod,       ONLY : tl_bvf_mix
# endif
      USE tl_diag_mod,          ONLY : tl_diag
# if defined WEAK_CONSTRAINT || defined FORCING_SV
      USE tl_forcing_mod,       ONLY : tl_forcing
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE tl_frc_adjust_mod,    ONLY : tl_frc_adjust
# endif
# ifdef GLS_MIXING_NOT_YET
!!    USE tl_gls_corstep_mod,   ONLY : tl_gls_corstep
!!    USE tl_gls_prestep_mod,   ONLY : tl_gls_prestep
# endif
# ifndef FORCING_SV
      USE tl_ini_fields_mod,    ONLY : tl_ini_fields, tl_ini_zeta
# endif
# ifdef LMD_MIXING_NOT_YET
!!    USE tl_lmd_vmix_mod,      ONLY : tl_lmd_vmix
# endif
# ifdef MY25_MIXING_NOT_YET
!!    USE tl_my25_corstep_mod,  ONLY : tl_my25_corstep
!!    USE tl_my25_prestep_mod,  ONLY : tl_my25_prestep
# endif
# ifdef NESTING
      USE nesting_mod,          ONLY : nesting
      USE tl_nesting_mod,       ONLY : tl_nesting
#  ifndef ONE_WAY
      USE nesting_mod,          ONLY : do_twoway
#  endif
# endif
# ifdef ADJUST_BOUNDARY
      USE tl_obc_adjust_mod,    ONLY : tl_obc_adjust
      USE tl_obc_adjust_mod,    ONLY : tl_obc2d_adjust
      USE tl_set_depth_mod,     ONLY : tl_set_depth_bry
# endif
      USE tl_omega_mod,         ONLY : tl_omega
# ifdef NEARSHORE_MELLOR_NOT_YET
!!    USE tl_radiation_stress_mod, ONLY : tl_radiation_stress
# endif
# ifndef TS_FIXED
      USE tl_rho_eos_mod,       ONLY : tl_rho_eos
# endif
      USE tl_rhs3d_mod,         ONLY : tl_rhs3d
# ifdef SEDIMENT_NOT_YET
!!    USE tl_sediment_mod,      ONLY : tl_sediment
# endif
# ifdef TL_AVERAGES
      USE tl_set_avg_mod,       ONLY : tl_set_avg
# endif
      USE tl_set_depth_mod,     ONLY : tl_set_depth
      USE tl_set_massflux_mod,  ONLY : tl_set_massflux
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE tl_set_tides_mod,     ONLY : tl_set_tides
# endif
      USE tl_set_vbc_mod,       ONLY : tl_set_vbc
      USE tl_set_zeta_mod,      ONLY : tl_set_zeta
      USE tl_step2d_mod,        ONLY : tl_step2d
# ifndef TS_FIXED
      USE tl_step3d_t_mod,      ONLY : tl_step3d_t
# endif
      USE tl_step3d_uv_mod,     ONLY : tl_step3d_uv
# ifdef FLOATS_NOT_YET
!!    USE tl_step_floats_mod,   ONLY : tl_step_floats
# endif
!!    USE wvelocity_mod,        ONLY : wvelocity
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

      integer :: Nsteps, Rsteps
      integer :: ig, il, istep, ng, nl, tile
      integer :: my_iif, next_indx1
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
!
      real(r8) :: MaxDT, my_StepTime
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Time-step tangent linear 3D primitive equations.
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
          CALL ntimesteps (iTLM, RunInterval, nl, Nsteps, Rsteps)
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
              nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
              nnew(ng)=3-nstp(ng)
              nrhs(ng)=nstp(ng)
# ifdef JEDI
              jic(ng)=jic(ng)+1
              time4jedi(ng)=time4jedi(ng)+dt(ng)
# endif
              tdays(ng)=time(ng)*sec2day
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
              CALL tl_get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots. Compute BASIC STATE depths and thickness.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL tl_set_data (ng, tile)
# ifdef FORWARD_READ
                CALL set_depth (ng, tile, iTLM)
# endif
              END DO
!$OMP BARRIER
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
                CALL set_massflux (ng, tile, iTLM)
              END DO
!$OMP BARRIER
            END DO
# endif

# if (defined WEAK_CONSTRAINT || defined FORCING_SV) && \
     !defined SP4DVAR
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution impulse forcing to
!  the representer model solution. Notice that the forcing is only
!  needed after finishing all inner loops. The forcing is continuous.
!  That is, it is time interpolated at every time-step from available
!  snapshots (FrequentImpulse=TRUE).
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
#  ifdef WEAK_CONSTRAINT
#   ifdef WEAK_NOINTERP
              IF ((iic(ng).gt.1).and.(iic(ng).ne.ntend(ng)+1).and.      &
     &            (MOD(iic(ng)-1,nadj(ng)).eq.0)) THEN
                IF (Master) THEN
                  WRITE (stdout,*) ' FORCING TLM at iic = ',iic(ng)
                END IF
#   endif
                IF (FrequentImpulse(ng)) THEN
                  DO tile=first_tile(ng),last_tile(ng),+1
                    CALL tl_forcing (ng, tile, kstp(ng), nstp(ng))
                    CALL tl_set_depth (ng, tile, iTLM)
                  END DO
!$OMP BARRIER
                END IF
#   ifdef WEAK_NOINTERP
              END IF
#   endif
#  else
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL tl_forcing (ng, tile, kstp(ng), nstp(ng))
                CALL tl_set_depth (ng, tile, iTLM)
              END DO
!$OMP BARRIER
#  endif
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
              IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface and compute initial level thicknesses and
!  depths.
!
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL tl_ini_zeta (ng, tile, iTLM)
                  CALL tl_set_depth (ng, tile, iTLM)
                END DO
!$OMP BARRIER
!
!  Initialize other state variables.
!
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL tl_ini_fields (ng, tile, iTLM)
                END DO
!$OMP BARRIER

#  ifdef NESTING
!
!  Extract donor grid initial data at contact points and store it in
!  REFINED structure so it can be used for the space-time interpolation.
!
                IF (RefinedGrid(ng)) THEN
                  CALL tl_nesting (ng, iTLM, ngetD)
                END IF
#  endif
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics. Compute BASIC STATE omega
!  vertical velocity.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL tl_set_massflux (ng, tile, iTLM)
# ifndef TS_FIXED
                CALL tl_rho_eos (ng, tile, iTLM)
# endif
# ifdef TIDE_GENERATING_FORCES
                CALL equilibrium_tide (ng, tile, iTLM)
# endif
                CALL tl_diag (ng, tile)
# ifdef FORWARD_READ
                CALL omega (ng, tile, iTLM)
# endif
# ifdef TLM_CHECK
                CALL tl_dotproduct (ng, tile, Lnew(ng))
# endif
              END DO
!$OMP BARRIER
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
!$OMP BARRIER
              END IF
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
!$OMP BARRIER
              END IF
            END DO
# endif

# ifdef NEARSHORE_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL tl_radiation_stress (ng, tile)
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
# if defined BULK_FLUXES_NOT_YET && !defined PRIOR_BULK_FLUXES
                CALL tl_bulk_flux (ng, tile)
# endif
# ifdef BBL_MODEL_NOT_YET
                CALL tl_bblm (ng, tile)
# endif
                CALL tl_set_vbc (ng, tile)
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
                CALL tl_set_tides (ng, tile)
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
                CALL tl_nesting (ng, iTLM, nbstr)
              END IF
            END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundaries.
!  Skip the last output timestep.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).lt.(ntend(ng)+1)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL tl_obc_adjust (ng, tile, Lbinp(ng))
                  CALL tl_set_depth_bry (ng, tile, iTLM)
                  CALL tl_obc2d_adjust (ng, tile, Lbinp(ng))
                END DO
!$OMP BARRIER
              END IF
            END DO
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Skip the last output timestep.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).lt.(ntend(ng)+1)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL tl_frc_adjust (ng, tile, Lfinp(ng))
                END DO
!$OMP BARRIER
              END IF
            END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute tangent linear vertical mixing coefficients for momentum and
!  tracers. Compute S-coordinate vertical velocity, diagnostically from
!  horizontal mass divergence.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
# if defined ANA_VMIX_NOT_YET
                CALL tl_ana_vmix (ng, tile)
# elif defined LMD_MIXING_NOT_YET
                CALL tl_lmd_vmix (ng, tile)
# elif defined BVF_MIXING_NOT_YET
                CALL tl_bvf_mix (ng, tile)
# endif
                CALL tl_omega (ng, tile, iTLM)
!!              CALL wvelocity (ng, tile, nstp(ng))
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
                CALL tl_set_zeta (ng, tile)
#  ifdef DIAGNOSTICS
!!              CALL set_diags (ng, tile)
#  endif
#  ifdef TL_AVERAGES
                CALL tl_set_avg (ng, tile)
#  endif
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
                CALL tl_nesting (ng, iTLM, nzeta)
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
              CALL tl_output (ng)
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
                CALL tl_nesting (ng, iTLM, nputD)
#  ifdef NESTING_DEBUG
                CALL tl_nesting (ng, iTLM, nmflx)
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
                CALL tl_rhs3d (ng, tile)
#  ifdef MY25_MIXING_NOT_YET
                CALL tl_my25_prestep (ng, tile)
#  elif defined GLS_MIXING_NOT_YET
                CALL tl_gls_prestep (ng, tile)
#  endif
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
                CALL tl_nesting (ng, iTLM, nrhst)
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
                    CALL tl_step2d (ng, tile)
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
                  CALL tl_nesting (ng, iTLM, n2dPS)
                END IF
                IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                  CALL tl_nesting (ng, iTLM, nmflx)
                END IF
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
                    CALL tl_step2d (ng, tile)
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
                  CALL tl_nesting (ng, iTLM, n2dCS)
                END IF
                IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
                  CALL tl_nesting (ng, iTLM, nmflx)
                END IF
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
              CALL tl_nesting (ng, iTLM, nmask)
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
                CALL tl_nesting (ng, iTLM, n2dfx)
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
                CALL tl_set_depth (ng, tile, iTLM)
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
              CALL tl_nesting (ng, iTLM, nzwgt)
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
                CALL tl_step3d_uv (ng, tile)
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
                CALL tl_nesting (ng, iTLM, n3duv)
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
                CALL tl_omega (ng, tile, iTLM)
# ifdef MY25_MIXING_NOT_YET
                CALL tl_my25_corstep (ng, tile)
# elif defined GLS_MIXING_NOT_YET
                CALL tl_gls_corstep (ng, tile)
# endif
# ifdef BIOLOGY
                CALL tl_biology (ng, tile)
# endif
# ifdef SEDIMENT_NOT_YET
                CALL tl_sediment (ng, tile)
# endif
              END DO
!$OMP BARRIER
            END DO

# ifndef TS_FIXED
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL tl_step3d_t (ng, tile)
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
                CALL tl_nesting (ng, iTLM, n3dTV)
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
                IF (do_twoway(iTLM, nl, il, ng, istep)) THEN
                  CALL tl_nesting (ng, iTLM, n2way)
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
                CALL tl_nesting (ng, iTLM, ngetD)
              END IF
            END DO
# endif

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
                CALL tl_step_floats (ng, Lstr, Lend)
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
!
!-----------------------------------------------------------------------
!  Advance time index and time clock.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iic(ng)=iic(ng)+1
              time(ng)=time(ng)+dt(ng)
              step_counter(ng)=step_counter(ng)-1
              CALL time_string (time(ng), time_code(ng))
            END DO

          END DO STEP_LOOP

        END DO NEST_LAYER

      END DO KERNEL_LOOP

      RETURN
      END SUBROUTINE tl_main3d
#else
      SUBROUTINE tl_main3d
      RETURN
      END SUBROUTINE tl_main3d
#endif
