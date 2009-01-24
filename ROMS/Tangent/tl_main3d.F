#include "cppdefs.h"
#if defined TANGENT && defined SOLVE3D
      SUBROUTINE tl_main3d (ng)
!
!svn $Id: tl_main3d.F 704 2008-08-13 06:20:10Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine is the main driver for tangent linear ROMS/TOMS when   !
!  configurated as a full  3D baroclinic  ocean model.  It  advances   !
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
# ifdef TLM_CHECK
      USE dotproduct_mod, ONLY : tl_dotproduct
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE frc_adjust_mod, ONLY : frc_TLadjust, load_TLforcing
# endif
# ifdef AIR_OCEAN_NOT_YET
      USE ocean_coupler_mod, ONLY : atmos_coupling
# endif
# ifdef WAVES_OCEAN_NOT_YET
      USE ocean_coupler_mod, ONLY : waves_coupling
# endif
# ifdef ASSIMILATION
!!    USE oi_update_mod, ONLY : oi_update
# endif
# ifdef FORWARD_READ
      USE omega_mod, ONLY : omega
      USE set_depth_mod, ONLY : set_depth
      USE set_massflux_mod, ONLY : set_massflux
# endif
# ifdef BIOLOGY
      USE tl_biology_mod, ONLY : tl_biology
# endif
# ifdef BBL_MODEL_NOT_YET
!!    USE tl_bbl_mod, ONLY : tl_bblm
# endif
# ifdef BULK_FLUXES_NOT_YET
!!    USE tl_bulk_flux_mod, ONLY : tl_bulk_flux
# endif
# ifdef BVF_MIXING_NOT_YET
!!    USE tl_bvf_mix_mod, ONLY : tl_bvf_mix
# endif
      USE tl_diag_mod, ONLY : tl_diag
# if defined WEAK_CONSTRAINT || defined IOM
      USE tl_forcing_mod, ONLY : tl_forcing
# endif
# ifdef GLS_MIXING_NOT_YET
!!    USE tl_gls_corstep_mod, ONLY : tl_gls_corstep
!!    USE tl_gls_prestep_mod, ONLY : tl_gls_prestep
# endif
      USE tl_ini_fields_mod, ONLY : tl_ini_fields, tl_ini_zeta
# ifdef LMD_MIXING_NOT_YET
!!    USE tl_lmd_vmix_mod, ONLY : tl_lmd_vmix
# endif
# ifdef MY25_MIXING_NOT_YET
!!    USE tl_my25_corstep_mod, ONLY : tl_my25_corstep
!!    USE tl_my25_prestep_mod, ONLY : tl_my25_prestep
# endif
      USE tl_omega_mod, ONLY : tl_omega
# ifdef NEARSHORE_MELLOR_NOT_YET
!!    USE tl_radiation_stress_mod, ONLY : tl_radiation_stress
# endif
# ifndef TS_FIXED
      USE tl_rho_eos_mod, ONLY : tl_rho_eos
# endif
      USE tl_rhs3d_mod, ONLY : tl_rhs3d
# ifdef SEDIMENT_NOT_YET
!!    USE tl_sediment_mod, ONLY : tl_sediment
# endif
# ifdef AVERAGES
!!    USE tl_set_avg_mod, ONLY : tl_set_avg
# endif
      USE tl_set_massflux_mod, ONLY : tl_set_massflux
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE tl_set_tides_mod, ONLY : tl_set_tides
# endif
      USE tl_set_vbc_mod, ONLY : tl_set_vbc
# ifdef SET_ZETA
      USE tl_set_zeta_mod, ONLY : tl_set_zeta
# endif
      USE tl_step2d_mod, ONLY : tl_step2d
# ifndef TS_FIXED
      USE tl_step3d_t_mod, ONLY : tl_step3d_t
# endif
      USE tl_step3d_uv_mod, ONLY : tl_step3d_uv
# ifdef FLOATS_NOT_YET
!!    USE tl_step_floats_mod, ONLY : tl_step_floats
# endif
!!    USE wvelocity_mod, ONLY : wvelocity
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: my_iif, next_indx1, subs, tile, thread
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
# if defined WEAK_CONSTRAINT || defined IOM
      real(r8) :: HalfDT
# endif
!
!=======================================================================
!  Time-step tangent linear 3D primitive equations.
!=======================================================================
!
!  Set time indices and time clock.
!
      nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
      nnew(ng)=3-nstp(ng)
      nrhs(ng)=nstp(ng)
      time(ng)=time(ng)+dt(ng)
      tdays(ng)=time(ng)*sec2day
      CALL time_string (time(ng), time_code(ng))

#if defined WEAK_CONSTRAINT || defined IOM
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution as impulse forcing to
!  the nonlinear solution. Notice that, in weak constraint, the forcing
!  is added intermittently (SporadicImpulse=TRUE) at available times
!  in the inner loop.
!-----------------------------------------------------------------------
!
      IF (SporadicImpulse) THEN
!
!  Add impulse forcing.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_forcing (ng, TILE, kstp(ng), nstp(ng))
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Read in next impulse forcing record, if any.
!
        FrcRec(ng)=FrcRec(ng)+1
        IF (FrcRec(ng).le.NrecFrc(ng)) THEN
          CALL get_state (ng, 7, 7, TLFname(ng), FrcRec(ng), 1)
        END IF
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Read in required data, if any, data from input NetCDF files.
!-----------------------------------------------------------------------
!
      CALL tl_get_data (ng)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots. Compute BASIC STATE depths and thickness.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL tl_set_data (ng, TILE)
# ifdef FORWARD_READ
          CALL set_depth (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
      IF (exit_flag.ne.NoError) RETURN

# ifdef FORWARD_READ
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE horizontal mass fluxes (Hz*u/n and Hz*v/m).
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL set_massflux (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  If not a restart, initialize all time levels and compute other
!  initial fields.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_ini_zeta (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Initialize other state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL tl_ini_fields (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics. Compute BASIC STATE omega
!  vertical velocity.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL tl_set_massflux (ng, TILE)
# ifndef TS_FIXED
          CALL tl_rho_eos (ng, TILE)
# endif
          CALL tl_diag (ng, TILE)
# ifdef FORWARD_READ
          CALL omega (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
      IF (exit_flag.ne.NoError) RETURN

# ifdef AIR_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Couple to atmospheric model every CoupleSteps(Iatmos) timesteps: get
!  air/sea fluxes.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).ne.ntstart(ng)).and.                                 &
     &    MOD(iic(ng)-1,CoupleSteps(Iatmos,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL atmos_coupling (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# ifdef WAVES_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Couple to waves model every CoupleSteps(Iwaves) timesteps: get
!  waves/sea fluxes.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).ne.ntstart(ng)).and.                                 &
     &    MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL waves_coupling (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# ifdef NEARSHORE_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL tl_radiation_stress (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!  Interpolate surface forcing increments.
# endif
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
# ifdef BULK_FLUXES_NOT_YET
          CALL tl_bulk_flux (ng, TILE)
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          CALL frc_TLadjust (ng, TILE, Lfinp(ng))
# endif
# ifdef BBL_MODEL_NOT_YET
          CALL tl_bblm (ng, TILE)
# endif
          CALL tl_set_vbc (ng, TILE)
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
          CALL tl_set_tides (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Load surface forcing into storage arrays.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,Lfout,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL load_TLforcing (ng, TILE, Lfout(ng))
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  Compute tangent linear vertical mixing coefficients for momentum and
!  tracers. Compute S-coordinate vertical velocity, diagnostically from
!  horizontal mass divergence.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
# if defined ANA_VMIX_NOT_YET
          CALL tl_ana_vmix (ng, TILE)
# elif defined LMD_MIXING_NOT_YET
          CALL tl_lmd_vmix (ng, TILE)
# elif defined BVF_MIXING_NOT_YET
          CALL tl_bvf_mix (ng, TILE)
# endif
          CALL tl_omega (ng, TILE)
!!        CALL wvelocity (ng, TILE, nstp(ng))
        END DO
      END DO
!$OMP END PARALLEL DO

# ifdef SET_ZETA
!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1     ! irreversible loop
          CALL tl_set_zeta (ng, TILE)
#  if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_UV
!!        CALL set_diags (ng, TILE)
#  endif
#  if defined AVERAGES && !defined ADJOINT
!!        CALL tl_set_avg (ng, TILE)
#  endif
        END DO
      END DO
!$OMP END PARALLEL DO

#  if !(defined WEAK_CONSTRAINT || defined IOM)
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
      CALL tl_output (ng)
      IF ((exit_flag.ne.NoError).or.(iic(ng).eq.(ntend(ng)+1))) RETURN
#  endif
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL tl_rhs3d (ng, TILE)
#  ifdef MY25_MIXING_NOT_YET
          CALL tl_my25_prestep (ng, TILE)
#  elif defined GLS_MIXING_NOT_YET
          CALL tl_gls_prestep (ng, TILE)
#  endif
        END DO
      END DO
!$OMP END PARALLEL DO

# else
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.  If applicable,
!  accumulate time-averaged output data.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&         SHARED(ng,Lnew,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
#  if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_UV
!!        CALL set_diags (ng, TILE)
#  endif
          CALL tl_rhs3d (ng, TILE)
#  ifdef MY25_MIXING_NOT_YET
          CALL tl_my25_prestep (ng, TILE)
#  elif defined GLS_MIXING_NOT_YET
          CALL tl_gls_prestep (ng, TILE)
#  endif
#  ifdef AVERAGES
!!        CALL tl_set_avg (ng, TILE)
#  endif
#  ifdef TLM_CHECK
          CALL tl_dotproduct (ng, TILE, Lnew(ng))
#  endif
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
      CALL tl_output (ng)
      IF ((exit_flag.ne.NoError).or.(iic(ng).eq.(ntend(ng)+1))) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and barotropic momentum components.
!-----------------------------------------------------------------------
!
      DO my_iif=1,nfast(ng)+1
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
        next_indx1=3-indx1(ng)
        IF (.not.PREDICTOR_2D_STEP(ng)) THEN
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL tl_step2d (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Set time indices for corrector step.
!
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
              CALL tl_step2d (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Time-step 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Time-step 3D momentum equations and couple with vertically
!  integrated equations.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL tl_step3d_uv (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL tl_omega (ng, TILE)
# ifdef MY25_MIXING_NOT_YET
          CALL tl_my25_corstep (ng, TILE)
# elif defined GLS_MIXING_NOT_YET
          CALL tl_gls_corstep (ng, TILE)
# endif
# ifdef BIOLOGY
          CALL tl_biology (ng, TILE)
# endif
# ifdef SEDIMENT_NOT_YET
          CALL tl_sediment (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO

# ifndef TS_FIXED
!
!-----------------------------------------------------------------------
!  Time-step tracer equations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL tl_step3d_t (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# ifdef FLOATS_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,chunk_size,Lstr,Lend)                  &
!$OMP&            SHARED(ng,numthreads,Nfloats)
      DO thread=0,numthreads-1
        chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
        Lstr=1+thread*chunk_size
        Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
        CALL tl_step_floats (ng, Lstr, Lend)
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
# endif

      RETURN
      END SUBROUTINE tl_main3d
#else
      SUBROUTINE tl_main3d
      RETURN
      END SUBROUTINE tl_main3d
#endif
