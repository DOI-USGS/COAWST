#include "cppdefs.h"
#if defined TL_IOMS && defined SOLVE3D
      SUBROUTINE rp_main3d (ng)
!
!svn $Id: rp_main3d.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the main driver for representers tangent linear     !
!  ROMS/TOMS when  configure as a  full 3D baroclinic ocean model.     !
!  It advances RPM primitive equations for a single time step.         !
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
      USE rp_biology_mod, ONLY : rp_biology
# endif
# ifdef BBL_MODEL_NOT_YET
!!    USE rp_bbl_mod, ONLY : rp_bblm
# endif
# if defined BULK_FLUXES && !defined NL_BULK_FLUXES
      USE rp_bulk_flux_mod, ONLY : rp_bulk_flux
# endif
# ifdef BVF_MIXING_NOT_YET
!!    USE rp_bvf_mix_mod, ONLY : rp_bvf_mix
# endif
      USE rp_diag_mod, ONLY : rp_diag
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE rp_frc_adjust_mod, ONLY : rp_frc_adjust
# endif
# ifdef GLS_MIXING_NOT_YET
!!    USE rp_gls_corstep_mod, ONLY : rp_gls_corstep
!!    USE rp_gls_prestep_mod, ONLY : rp_gls_prestep
# endif
      USE rp_ini_fields_mod, ONLY : rp_ini_fields, rp_ini_zeta
# ifdef LMD_MIXING_NOT_YET
!!    USE rp_lmd_vmix_mod, ONLY : rp_lmd_vmix
# endif
# ifdef MY25_MIXING
!!    USE rp_my25_corstep_mod, ONLY : rp_my25_corstep
!!    USE rp_my25_prestep_mod, ONLY : rp_my25_prestep
# endif
# ifdef ADJUST_BOUNDARY
      USE rp_obc_adjust_mod, ONLY : rp_obc_adjust
      USE rp_obc_adjust_mod, ONLY : rp_obc2d_adjust
      USE rp_set_depth_mod, ONLY : rp_set_depth_bry
# endif
      USE rp_omega_mod, ONLY : rp_omega
# ifdef NEARSHORE_MELLOR_NOT_YET
!!    USE rp_radiation_stress_mod, ONLY : rp_radiation_stress
# endif
# ifndef TS_FIXED
      USE rp_rho_eos_mod, ONLY : rp_rho_eos
# endif
      USE rp_rhs3d_mod, ONLY : rp_rhs3d
# ifdef SEDIMENT_NOT_YET
!!    USE rp_sediment_mod, ONLY : rp_sediment
# endif
# ifdef AVERAGES
!!    USE rp_set_avg_mod, ONLY : rp_set_avg
# endif
# ifdef MOVE_SET_DEPTH
      USE rp_set_depth_mod, ONLY : rp_set_depth
# endif
      USE rp_set_massflux_mod, ONLY : rp_set_massflux
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE rp_set_tides_mod, ONLY : rp_set_tides
# endif
      USE rp_set_vbc_mod, ONLY : rp_set_vbc
      USE rp_set_zeta_mod, ONLY : rp_set_zeta
      USE rp_step2d_mod, ONLY : rp_step2d
# ifndef TS_FIXED
      USE rp_step3d_t_mod, ONLY : rp_step3d_t
# endif
      USE rp_step3d_uv_mod, ONLY : rp_step3d_uv
# ifdef FLOATS_NOT_YET
!!    USE rp_step_floats_mod, ONLY : rp_step_floats
# endif
# ifdef WEAK_CONSTRAINT
      USE tl_forcing_mod, ONLY : tl_forcing
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
!
!-----------------------------------------------------------------------
!  Read in required data, if any, data from input NetCDF files.
!-----------------------------------------------------------------------
!
      CALL rp_get_data (ng)
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
          CALL rp_set_data (ng, TILE)
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
# ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution impulse forcing to
!  the representer model solution. Notice that the forcing is only
!  needed after finishing all inner loops. The forcing is continuous.
!  That is, it is time interpolated at every time-step from available
!  snapshots (FrequentImpulse=TRUE).
!-----------------------------------------------------------------------
!
      IF (FrequentImpulse) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_forcing (ng, TILE, kstp(ng), nstp(ng))
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
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
            CALL rp_ini_zeta (ng, TILE, iRPM)
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
            CALL rp_ini_fields (ng, TILE, iRPM)
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
          CALL rp_set_massflux (ng, TILE)
# ifndef TS_FIXED
          CALL rp_rho_eos (ng, TILE)
# endif
          CALL rp_diag (ng, TILE)
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
     &    MOD(iic(ng),CoupleSteps(Iwaves,ng)).eq.0) THEN
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
          CALL rp_radiation_stress (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  Set fields for vertical boundary conditions. Process tidal forcing,
!  if any.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
# if defined BULK_FLUXES && !defined NL_BULK_FLUXES
          CALL rp_bulk_flux (ng, TILE)
# endif
# ifdef BBL_MODEL_NOT_YET
          CALL rp_bblm (ng, TILE)
# endif
          CALL rp_set_vbc (ng, TILE)
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
          CALL rp_set_tides (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundaries.
!  Skip the last output timestep.
!-----------------------------------------------------------------------
!
      IF ((Nrun.ne.1).and.(iic(ng).lt.(ntend(ng)+1))) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL rp_obc_adjust (ng, TILE, Lbinp(ng))
            CALL rp_set_depth_bry (ng, TILE)
            CALL rp_obc2d_adjust (ng, TILE, Lbinp(ng))
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Skip the last output timestep.
!-----------------------------------------------------------------------
!
      IF (iic(ng).lt.(ntend(ng)+1)) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL rp_frc_adjust (ng, TILE, Lfinp(ng))
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
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
          CALL rp_ana_vmix (ng, TILE)
# elif defined LMD_MIXING_NOT_YET
          CALL rp_lmd_vmix (ng, TILE)
# elif defined BVF_MIXING_NOT_YET
          CALL rp_bvf_mix (ng, TILE)
# endif
          CALL rp_omega (ng, TILE)
!!        CALL wvelocity (ng, TILE, nstp(ng))
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Set free-surface to it time-averaged value.  If applicable,
!  accumulate time-averaged output data which needs a irreversible
!  loop in shared-memory jobs.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&         SHARED(ng,Lnew,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1     ! irreversible loop
          CALL rp_set_zeta (ng, TILE)
#  ifdef DIAGNOSTICS
!!        CALL rp_set_diags (ng, TILE)
#  endif
#  ifdef AVERAGES
!!        CALL rp_set_avg (ng, TILE)
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
      CALL rp_output (ng)
      IF ((exit_flag.ne.NoError).or.(iic(ng).eq.(ntend(ng)+1))) RETURN
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&         SHARED(ng,Lnew,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL rp_rhs3d (ng, TILE)
#  ifdef MY25_MIXING_NOT_YET
          CALL rp_my25_prestep (ng, TILE)
#  elif defined GLS_MIXING_NOT_YET
          CALL rp_gls_prestep (ng, TILE)
#  endif
        END DO
      END DO
!$OMP END PARALLEL DO
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
            CALL rp_step2d (ng, TILE)
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
              CALL rp_step2d (ng, TILE)
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
      END DO

# ifdef MOVE_SET_DEPTH
!
!-----------------------------------------------------------------------
!  Recompute depths and thicknesses using the new time filtered
!  free-surface.  This call was moved from "rp_step2d" to here.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL rp_set_depth (ng, TILE)
        END DO
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL rp_step3d_uv (ng, TILE)
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
          CALL rp_omega (ng, TILE)
# ifdef MY25_MIXING_NOT_YET
          CALL rp_my25_corstep (ng, TILE)
# elif defined GLS_MIXING_NOT_YET
          CALL rp_gls_corstep (ng, TILE)
# endif
# ifdef BIOLOGY
          CALL rp_biology (ng, TILE)
# endif
# ifdef SEDIMENT_NOT_YET
          CALL rp_sediment (ng, TILE)
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
          CALL rp_step3d_t (ng, TILE)
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
      IF (Lfloats(Ng)) THEN
!$OMP PARALLEL DO PRIVATE(thread,chunk_size,Lstr,Lend)                  &
!$OMP&            SHARED(ng,numthreads,Nfloats)
        DO thread=0,numthreads-1
          chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
          Lstr=1+thread*chunk_size
          Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
          CALL rp_step_floats (ng, Lstr, Lend)
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
# endif

      RETURN
      END SUBROUTINE rp_main3d
#else
      SUBROUTINE rp_main3d
      RETURN
      END SUBROUTINE rp_main3d
#endif
