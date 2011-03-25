#include "cppdefs.h"
#if defined NONLINEAR && defined SOLVE3D
# if !defined REFINED_GRID
      SUBROUTINE main3d
# else
      SUBROUTINE main3d(ng)
# endif
!
!svn $Id: main3d.F 838 2008-11-17 04:22:18Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as a full 3D baroclinic ocean model.  It  advances     !
!  forward the primitive equations for a single time step.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef MODEL_COUPLING
      USE mod_coupler
# endif
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
# ifdef ANA_VMIX
      USE analytical_mod, ONLY : ana_vmix
# endif
# ifdef BIOLOGY
      USE biology_mod, ONLY : biology
# endif
# ifdef BBL_MODEL
      USE bbl_mod, ONLY : bblm
# endif
# ifdef BULK_FLUXES
      USE bulk_flux_mod, ONLY : bulk_flux
# endif
# ifdef ATM2OCN_FLUXES
      USE atm2ocn_flux_mod, ONLY : atm2ocn_flux
# endif
# ifdef BVF_MIXING
      USE bvf_mix_mod, ONLY : bvf_mix
# endif
      USE diag_mod, ONLY : diag
# ifdef TLM_CHECK
      USE dotproduct_mod, ONLY : nl_dotproduct
# endif
# if defined W4DPSAS || defined NLM_OUTER || \
     defined W4DPSAS_SENSITIVITY
      USE forcing_mod, ONLY : forcing
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE frc_adjust_mod, ONLY : frc_adjust
# endif
# ifdef GLS_MIXING
      USE gls_corstep_mod, ONLY : gls_corstep
      USE gls_prestep_mod, ONLY : gls_prestep
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
      USE hmixing_mod, ONLY : hmixing
# endif
      USE ini_fields_mod, ONLY : ini_fields, ini_zeta
# ifdef LMD_MIXING
      USE lmd_vmix_mod, ONLY : lmd_vmix
# endif
# ifdef MY25_MIXING
      USE my25_corstep_mod, ONLY : my25_corstep
      USE my25_prestep_mod, ONLY : my25_prestep
# endif
# if defined ADJUST_BOUNDARY
      USE obc_adjust_mod, ONLY : obc_adjust, load_obc
# endif
# ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2atm_coupling
# endif
# ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2wav_coupling
# endif
# ifdef ASSIMILATION
      USE oi_update_mod, ONLY : oi_update
# endif
      USE omega_mod, ONLY : omega
# ifdef WEC_MELLOR
      USE radiation_stress_mod, ONLY : radiation_stress
# endif
# ifdef WEC_VF
      USE wec_vf_mod, ONLY : wec_vf
#  if defined WDISS_THORGUZA || defined WDISS_CHURTHOR
      USE wec_dissip_mod, ONLY : wec_dissip
#  endif
#  ifdef WEC_ROLLER
      USE wec_roller_mod, ONLY : wec_roller
#  endif
      USE wec_stokes_mod, ONLY : wec_stokes
#  ifdef BOTTOM_STREAMING
      USE wec_streaming_mod, ONLY : wec_streaming
#  endif
#  ifdef WAVE_MIXING
      USE wec_wave_mix_mod, ONLY : wec_wave_mix
#  endif
# endif
# ifndef TS_FIXED
      USE rho_eos_mod, ONLY : rho_eos
# endif
      USE rhs3d_mod, ONLY : rhs3d
# ifdef SEDIMENT
      USE sediment_mod, ONLY : sediment
# endif
# if defined AVERAGES && !defined ADJOINT
      USE set_avg_mod, ONLY : set_avg
# endif
# ifdef MOVE_SET_DEPTH
      USE set_depth_mod, ONLY : set_depth
# endif
      USE set_massflux_mod, ONLY : set_massflux
# if defined SSH_TIDES || defined UV_TIDES
      USE set_tides_mod, ONLY : set_tides
# endif
      USE set_vbc_mod, ONLY : set_vbc
      USE set_zeta_mod, ONLY : set_zeta
      USE step2d_mod, ONLY : step2d
# ifndef TS_FIXED
      USE step3d_t_mod, ONLY : step3d_t
# endif
      USE step3d_uv_mod, ONLY : step3d_uv
# ifdef FLOATS
      USE step_floats_mod, ONLY : step_floats
# endif
      USE wvelocity_mod, ONLY : wvelocity
# ifdef WEC
      USE wstvelocity_mod, ONLY : wstvelocity
# endif
# ifdef COMPOSED_GRID
      USE nesting_mod, ONLY : nesting
# endif
# ifdef REFINED_GRID
      USE get_2dparent_data_mod, ONLY : get_2dparent_data
      USE get_3dparent_data_mod, ONLY : get_3dparent_data
      USE set_2dparent_data_mod, ONLY : set_2dparent_data
      USE set_3dparent_data_mod, ONLY : set_3dparent_data
      USE set_2dchild_data_mod, ONLY : set_2dchild_data
      USE set_3dchild_data_mod, ONLY : set_3dchild_data
# endif
!
      implicit none
!
!  Imported variable declarations.
!
# if !defined REFINED_GRID
      integer :: ng
# else
      integer, intent(in) :: ng
# endif
!  Local variable declarations.
!
      integer :: my_iif, next_indx1, subs, tile, thread
# ifdef FLOATS
      integer :: Lend, Lstr, chunk_size
# endif
# if defined REFINED_GRID
      integer :: ngp
      ngp=myparent(ng)
# endif
!
!=======================================================================
!  Time-step nonlinear 3D primitive equations.
!=======================================================================
!
!  Set time indices and time clock.
!
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
        nnew(ng)=3-nstp(ng)
        nrhs(ng)=nstp(ng)
        time(ng)=time(ng)+dt(ng)
        tdays(ng)=time(ng)*sec2day
        CALL time_string (time(ng), time_code(ng))
# if !defined REFINED_GRID
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        CALL get_data (ng)
        IF (exit_flag.ne.NoError) RETURN
# if !defined REFINED_GRID
      END DO
# endif
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL set_data (ng, TILE)
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Get data from parent grids.
!-----------------------------------------------------------------------
!
# ifdef REFINED_GRID
! get data from parent grid at time of parent. Only need to 
! do this once per child loop.
        IF (ng.gt.1) THEN
          IF (time(ng).eq.time(ngp)) THEN
            CALL get_2dparent_data (ng, TILE)
            CALL get_3dparent_data (ng, TILE)
          END IF
        END IF
! put child data back to parent 
        IF ((ng.gt.Ngrids).and.(iic(ng).gt.1)) THEN
!          CALL set_2dchild_data (ng, TILE)
!          CALL set_3dchild_data (ng, TILE)
        END IF
# endif
!
# if defined W4DPSAS || defined NLM_OUTER || \
     defined W4DPSAS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution impulse forcing to
!  the nonlinear model solution. Notice that the forcing is only needed
!  after finishing all the inner loops. The forcing is continuous.
!  That is, it is time interpolated at every time-step from available
!  snapshots (FrequentImpulse=TRUE).
!-----------------------------------------------------------------------
!
      IF (FrequentImpulse) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
        DO ng=1,Ngrids
# endif
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL forcing (ng, TILE, kstp(ng), nstp(ng))
            END DO
          END DO
# if !defined REFINED_GRID
        END DO
# endif
!$OMP END PARALLEL DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
!
!  Initialize free-surface.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        IF (iic(ng).eq.ntstart(ng)) THEN
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL ini_zeta (ng, TILE, iNLM)
            END DO
          END DO
        END IF
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
!
!  Initialize other state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        IF (iic(ng).eq.ntstart(ng)) THEN
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*(thread+1)-1,subs*thread,-1
              CALL ini_fields (ng, TILE, iNLM)
            END DO
          END DO
        END IF
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL set_massflux (ng, TILE)
# ifndef TS_FIXED
            CALL rho_eos (ng, TILE)
# endif
            CALL diag (ng, TILE)
# ifdef TLM_CHECK
            CALL nl_dotproduct (ng, TILE, Lnew(ng))
# endif
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
      IF (exit_flag.ne.NoError) RETURN
#  ifdef COMPOSED_GRID
      CALL nesting (5)
#  endif

# ifdef AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Couple ocean to atmosphere model every "CoupleSteps(Iatmos)"
!  timesteps: get air/sea fluxes.
!-----------------------------------------------------------------------
!
#  if !defined REFINED_GRID
      DO ng=1,Ngrids
#  endif
      IF ((iic(ng).ne.ntstart(ng)).and.                                 &
     &    MOD(iic(ng),nOCN_ATM(ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*(thread+1)-1,subs*thread,-1
              CALL ocn2atm_coupling (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
#  if !defined REFINED_GRID
      END DO
#  endif
# endif

# ifdef WAVES_OCEAN
!
!-----------------------------------------------------------------------
!  Couple to ocean to waves model every "CoupleSteps(Iwaves)"
!  timesteps: get waves/ocean fluxes.
!-----------------------------------------------------------------------
!
#  if !defined REFINED_GRID
      DO ng=1,Ngrids
#  endif
      IF ((iic(ng).ne.ntstart(ng)).and.                                 &
#ifdef REFINED_GRID
     &    MOD(iic(ng)-1,nOCN_WAV(ng)).eq.0) THEN
#else
     &    MOD(iic(ng),nOCN_WAV(ng)).eq.0) THEN
#endif
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL ocn2wav_coupling (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
#  if !defined REFINED_GRID
      END DO
#  endif
# endif
# ifdef WEC_VF
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
#  if !defined REFINED_GRID
      DO ng=1,Ngrids
#  endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
#  if defined WDISS_THORGUZA || defined WDISS_CHURTHOR
            CALL wec_dissip (ng, TILE)
#  endif
#  ifdef WEC_ROLLER
            CALL wec_roller (ng, TILE)
#  endif
            CALL wec_stokes (ng, TILE)
            CALL wec_vf (ng, TILE)
#  ifdef BOTTOM_STREAMING
            CALL wec_streaming (ng, TILE)
#  endif
          END DO
        END DO
#  if !defined REFINED_GRID
      END DO
#  endif
!$OMP END PARALLEL DO
#  ifdef COMPOSED_GRID
      CALL nesting (8)
#  endif
# endif
# ifdef WEC_MELLOR
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
#  if !defined REFINED_GRID
      DO ng=1,Ngrids
#  endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL radiation_stress (ng, TILE)
          END DO
        END DO
#  if !defined REFINED_GRID
      END DO
#  endif
!$OMP END PARALLEL DO
#  ifdef COMPOSED_GRID
      CALL nesting (8)
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
# ifdef BULK_FLUXES
#  if defined FOUR_DVAR && defined NL_BULK_FLUXES
            IF (Nrun.eq.1) CALL bulk_flux (ng, TILE)
#  else
            CALL bulk_flux (ng, TILE)
#  endif

# elif defined ATM2OCN_FLUXES
            CALL atm2ocn_flux (ng, TILE)
# endif
# ifdef BBL_MODEL
            CALL bblm (ng, TILE)
# endif
            CALL set_vbc (ng, TILE)
# if defined SSH_TIDES || defined UV_TIDES
            CALL set_tides (ng, TILE)
# endif
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
#  ifdef COMPOSED_GRID
      CALL nesting (9)
#  endif
# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundary.
!  Load open boundary into storage arrays. Skip the last output
!  timestep.
!-----------------------------------------------------------------------
!
      IF (iic(ng).lt.(ntend(ng)+1)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL obc_adjust (ng, TILE, Lbinp(ng))
            CALL load_obc (ng, TILE, Lbout(ng))
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Load surface forcing into storage arrays. Skip the last output
!  timestep.
!-----------------------------------------------------------------------
!
      IF (iic(ng).lt.(ntend(ng)+1)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL frc_adjust (ng, TILE, Lfinp(ng))
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute time-dependent vertical/horizontal mixing coefficients for
!  momentum and tracers. Compute S-coordinate vertical velocity,
!  diagnostically from horizontal mass divergence.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,nstp,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
# if defined ANA_VMIX
            CALL ana_vmix (ng, TILE, iNLM)
# elif defined LMD_MIXING
            CALL lmd_vmix (ng, TILE)
# elif defined BVF_MIXING
            CALL bvf_mix (ng, TILE)
# endif
# if defined DIFF_3DCOEF || defined VISC_3DCOEF
          CALL hmixing (ng, TILE)
# endif
            CALL omega (ng, TILE)
            CALL wvelocity (ng, TILE, nstp(ng))
# ifdef WEC
            CALL wstvelocity (ng, TILE, nstp(ng))
# endif
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
#  ifdef COMPOSED_GRID
      CALL nesting (10)
#  endif

!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL set_zeta (ng, TILE)
#  ifdef DIAGNOSTICS
            CALL set_diags (ng, TILE)
#  endif
#  if defined AVERAGES && !defined ADJOINT
            CALL set_avg (ng, TILE)
#  endif
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
#  ifdef COMPOSED_GRID
      CALL nesting (11)
#  endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        CALL output (ng)
        IF ((exit_flag.ne.NoError).or.                                    &
     &      ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) RETURN
# if !defined REFINED_GRID
      END DO
# endif
# ifdef REFINED_GRID
! interpolate the parent data to child time.
        IF (ng.gt.1) THEN
          CALL set_2dparent_data (ng, TILE)
          CALL set_3dparent_data (ng, TILE)
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL rhs3d (ng, TILE)
#  ifdef MY25_MIXING
            CALL my25_prestep (ng, TILE)
#  elif defined GLS_MIXING
            CALL gls_prestep (ng, TILE)
#  endif
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
#  ifdef COMPOSED_GRID
      CALL nesting (12)
#  endif
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and barotropic momentum components.
!-----------------------------------------------------------------------
!
# if !defined REFINED_GRID
      DO my_iif=1,MAXVAL(nfast)+1
        DO ng=1,Ngrids
# else
      DO my_iif=1,nfast(ng)+1
# endif
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
          next_indx1=3-indx1(ng)
          IF (.not.PREDICTOR_2D_STEP(ng).and.                           &
     &        my_iif.le.(nfast(ng)+1)) THEN
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*(thread+1)-1,subs*thread,-1
                CALL step2d (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
        END IF
# if !defined REFINED_GRID
        END DO
# endif

# ifdef COMPOSED_GRID
        CALL nesting (14)
# endif
!
!  Set time indices for corrector step.
!
# if !defined REFINED_GRID
        DO ng=1,Ngrids
# endif
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL step2d (ng, TILE)
              END DO
            END DO
!$OMP END PARALLEL DO
          END IF
# if !defined REFINED_GRID
        END DO
# endif
# ifdef COMPOSED_GRID
        CALL nesting (15)
# endif
      END DO
# ifdef COMPOSED_GRID
      CALL nesting (21)
# endif
!
!-----------------------------------------------------------------------
!  Recompute depths and thicknesses using the new time filtered
!  free-surface.  This call was moved from "step2d" to here.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL set_depth (ng, TILE)
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Time-step 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Time-step 3D momentum equations and couple with vertically
!  integrated equations.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL step3d_uv (ng, TILE)
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
# ifdef COMPOSED_GRID
      CALL nesting (16)
# endif
!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,nnew,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL omega (ng, TILE)
# ifdef MY25_MIXING
            CALL my25_corstep (ng, TILE)
# elif defined GLS_MIXING
            CALL gls_corstep (ng, TILE)
# endif
# if defined WEC_VF && defined WAVE_MIXING
            CALL wec_wave_mix (ng, TILE)
# endif
# ifdef BIOLOGY
            CALL biology (ng, TILE)
# endif
# ifdef SEDIMENT
            CALL sediment (ng, TILE)
# endif
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
# ifdef COMPOSED_GRID
      CALL nesting (17)
# endif

# ifndef TS_FIXED
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL step3d_t (ng, TILE)
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
#  ifdef COMPOSED_GRID
      CALL nesting (18)
#  endif
# endif

# ifdef ASSIMILATION
!
!-----------------------------------------------------------------------
!  Assimilate observations via Optimal Interpolation.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL oi_update (ng, TILE)
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO
# endif

# ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories.
!-----------------------------------------------------------------------
!
#  if !defined REFINED_GRID
      DO ng=1,Ngrids
#  endif
        IF (Lfloats(Ng)) THEN
!$OMP PARALLEL DO PRIVATE(thread,chunk_size,Lstr,Lend)                  &
!$OMP&            SHARED(ng,numthreads,Nfloats)
          DO thread=0,numthreads-1
            chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
            Lstr=1+thread*chunk_size
            Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
            CALL step_floats (ng, Lstr, Lend)
          END DO
!$OMP END PARALLEL DO
!
!  Shift floats time indices.
!
          nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
          nf(ng)  =MOD(nf(ng)  +1,NFT+1)
          nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
          nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
          nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
        END IF
#  if !defined REFINED_GRID
      END DO
#  endif
# endif
      RETURN
      END SUBROUTINE main3d
#else
      SUBROUTINE main3d
      RETURN
      END SUBROUTINE main3d
#endif
