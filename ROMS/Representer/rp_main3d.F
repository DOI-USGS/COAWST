#include "cppdefs.h"
#if defined TL_IOMS && defined SOLVE3D
      SUBROUTINE rp_main3d (RunInterval)
!
!svn $Id: rp_main3d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the main driver for representers tangent linear     !
!  ROMS/TOMS when  configure as a  full 3D baroclinic ocean model.     !
!  It advances forward the  representer model equations  for  all      !
!  nested grids, if any, for the specified time interval (seconds),    !
!  RunInterval.                                                        !
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
      USE rp_set_depth_mod, ONLY : rp_set_depth
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
# ifdef RP_AVERAGES
      USE tl_set_avg_mod, ONLY : tl_set_avg
# endif
!!    USE wvelocity_mod, ONLY : wvelocity
!
      implicit none
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: ng, tile
      integer :: my_iif, next_indx1
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
      real(r8) :: MaxDT, my_StepTime
!
!=======================================================================
!  Time-step representers tangent linear 3D primitive equations.
!=======================================================================
!
      my_StepTime=0.0_r8
      MaxDT=MAXVAL(dt)

      STEP_LOOP : DO WHILE (my_StepTime.le.(RunInterval+0.5_r8*MaxDT))

        my_StepTime=my_StepTime+MaxDT
!
!  Set time indices and time clock.
!
        DO ng=1,Ngrids
          iic(ng)=iic(ng)+1
          nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
          nnew(ng)=3-nstp(ng)
          nrhs(ng)=nstp(ng)
!$OMP MASTER
          time(ng)=time(ng)+dt(ng)
          tdays(ng)=time(ng)*sec2day
          CALL time_string (time(ng), time_code(ng))
!$OMP END MASTER
        END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Read in required data, if any, data from input NetCDF files.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP MASTER
          CALL rp_get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots. Compute BASIC STATE depths and thickness.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL rp_set_data (ng, tile)
# ifdef FORWARD_READ
            CALL set_depth (ng, tile)
# endif
          END DO
!$OMP BARRIER
        END DO
        IF (exit_flag.ne.NoError) RETURN

# ifdef FORWARD_READ
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE horizontal mass fluxes (Hz*u/n and Hz*v/m).
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL set_massflux (ng, tile)
          END DO
!$OMP BARRIER
        END DO
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
      DO ng=1,Ngrids
        IF (FrequentImpulse(ng)) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tl_forcing (ng, tile, kstp(ng), nstp(ng))
          END DO
!$OMP BARRIER
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  If not a restart, initialize all time levels and compute other
!  initial fields.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface and compute initial level thicknesses and
!  depths.
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL rp_ini_zeta (ng, tile, iRPM)
              CALL rp_set_depth (ng, tile)
            END DO
!$OMP BARRIER
!
!  Initialize other state variables.
!
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL rp_ini_fields (ng, tile, iRPM)
            END DO
!$OMP BARRIER
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Compute horizontal mass fluxes (Hz*u/n and Hz*v/m), density related
!  quatities and report global diagnostics. Compute BASIC STATE omega
!  vertical velocity.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL rp_set_massflux (ng, tile)
# ifndef TS_FIXED
            CALL rp_rho_eos (ng, tile)
# endif
            CALL rp_diag (ng, tile)
# ifdef FORWARD_READ
            CALL omega (ng, tile)
# endif
          END DO
!$OMP BARRIER
        END DO
        IF (exit_flag.ne.NoError) RETURN

# ifdef AIR_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Couple to atmospheric model every CoupleSteps(Iatmos) timesteps: get
!  air/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng)-1,CoupleSteps(Iatmos,ng)).eq.0) THEN
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL atmos_coupling (ng, tile)
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif

# ifdef WAVES_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Couple to waves model every CoupleSteps(Iwaves) timesteps: get
!  waves/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng),CoupleSteps(Iwaves,ng)).eq.0) THEN
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
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL rp_radiation_stress (ng, tile)
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
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
# if defined BULK_FLUXES && !defined NL_BULK_FLUXES
            CALL rp_bulk_flux (ng, tile)
# endif
# ifdef BBL_MODEL_NOT_YET
            CALL rp_bblm (ng, tile)
# endif
            CALL rp_set_vbc (ng, tile)
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
            CALL rp_set_tides (ng, tile)
# endif
          END DO
!$OMP BARRIER
        END DO

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundaries.
!  Skip the last output timestep.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((Nrun.ne.1).and.(iic(ng).lt.(ntend(ng)+1))) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL rp_obc_adjust (ng, tile, Lbinp(ng))
              CALL rp_set_depth_bry (ng, tile)
              CALL rp_obc2d_adjust (ng, tile, Lbinp(ng))
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
        DO ng=1,Ngrids
          IF (iic(ng).lt.(ntend(ng)+1)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL rp_frc_adjust (ng, tile, Lfinp(ng))
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
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
# if defined ANA_VMIX_NOT_YET
            CALL rp_ana_vmix (ng, tile)
# elif defined LMD_MIXING_NOT_YET
            CALL rp_lmd_vmix (ng, tile)
# elif defined BVF_MIXING_NOT_YET
            CALL rp_bvf_mix (ng, tile)
# endif
            CALL rp_omega (ng, tile)
!!          CALL wvelocity (ng, tile, nstp(ng))
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
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1         ! irreversible
            CALL rp_set_zeta (ng, tile)
#  ifdef DIAGNOSTICS
!!          CALL rp_set_diags (ng, tile)
#  endif
#  ifdef RP_AVERAGES
            CALL tl_set_avg (ng, tile)
#  endif
          END DO
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP MASTER
          CALL rp_output (ng)
!$OMP END MASTER
!$OMP BARRIER
          IF ((exit_flag.ne.NoError).or.                                &
     &        ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for 3D equations.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL rp_rhs3d (ng, tile)
#  ifdef MY25_MIXING_NOT_YET
            CALL rp_my25_prestep (ng, tile)
#  elif defined GLS_MIXING_NOT_YET
            CALL rp_gls_prestep (ng, tile)
#  endif
          END DO
!$OMP BARRIER
        END DO
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
          DO ng=1,Ngrids
            next_indx1=3-indx1(ng)
            IF (.not.PREDICTOR_2D_STEP(ng).and.                         &
     &          my_iif.le.(nfast(ng)+1)) THEN
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
                CALL rp_step2d (ng, tile)
              END DO
!$OMP BARRIER
            END IF
          END DO
!
!  Set time indices for corrector step.
!
          DO ng=1,Ngrids
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
                CALL rp_step2d (ng, tile)
              END DO
!$OMP BARRIER
            END IF
          END DO

        END DO LOOP_2D
!
!-----------------------------------------------------------------------
!  Recompute depths and thicknesses using the new time filtered
!  free-surface.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL rp_set_depth (ng, tile)
          END DO
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  Time-step 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Time-step 3D momentum equations and couple with vertically
!  integrated equations.
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL rp_step3d_uv (ng, tile)
          END DO
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  Time-step vertical mixing turbulent equations and passive tracer
!  source and sink terms, if applicable.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL rp_omega (ng, tile)
# ifdef MY25_MIXING_NOT_YET
            CALL rp_my25_corstep (ng, tile)
# elif defined GLS_MIXING_NOT_YET
            CALL rp_gls_corstep (ng, tile)
# endif
# ifdef BIOLOGY
            CALL rp_biology (ng, tile)
# endif
# ifdef SEDIMENT_NOT_YET
            CALL rp_sediment (ng, tile)
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
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL rp_step3d_t (ng, tile)
          END DO
        END DO
!$OMP BARRIER
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
        DO ng=1,Ngrids
          IF (Lfloats(Ng)) THEN
#  ifdef _OPENMP
            chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
            Lstr=1+my_thread*chunk_size
            Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
#  else
            Lstr=1
            Lend=Nfloats(ng)
#  endif
            CALL rp_step_floats (ng, Lstr, Lend)
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

      RETURN
      END SUBROUTINE rp_main3d
#else
      SUBROUTINE rp_main3d
      RETURN
      END SUBROUTINE rp_main3d
#endif
