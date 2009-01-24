#include "cppdefs.h"
#ifdef TL_IOMS
      SUBROUTINE rp_initial (ng)
!
!svn $Id: rp_initial.F 594 2008-04-01 18:11:31Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes representers tangent linear model.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BBL_MODEL_NOT_YET
      USE mod_bbl
# endif
# ifdef SOLVE3D
      USE mod_coupling      
# endif
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_grid
      USE mod_forces
      USE mod_iounits
      USE mod_ncparam
# ifdef SOLVE3D
      USE mod_mixing
# endif
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
# ifdef ANALYTICAL
      USE analytical_mod
# endif
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
# if defined GRADIENT_CHECK || defined TLM_CHECK
      USE ini_adjust_mod, ONLY : tl_ini_perturb
# endif
      USE metrics_mod, ONLY : metrics
# if !(defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET)
      USE rp_set_depth_mod, ONLY : rp_bath
# endif
# ifdef SOLVE3D
      USE rp_set_depth_mod, ONLY : rp_set_depth
      USE rp_omega_mod, ONLY : rp_omega
      USE rp_rho_eos_mod, ONLY : rp_rho_eos
      USE rp_set_massflux_mod, ONLY : rp_set_massflux
      USE set_depth_mod, ONLY : set_depth
      USE omega_mod, ONLY : omega
      USE rho_eos_mod, ONLY : rho_eos
      USE set_massflux_mod, ONLY : set_massflux
# endif

      USE stiffness_mod, ONLY : stiffness
# if defined PROPAGATOR || \
    (defined MASKING && (defined READ_WATER || defined WRITE_WATER ))
      USE wpoints_mod, ONLY : wpoints
# endif
# ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : waves_coupling
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: update = .FALSE.

      integer :: LBi, UBi, LBj, UBj
      integer :: IniRec, Tindex, subs, tile, thread, wrtRec

      integer :: my_numthreads
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
      IF (Master) THEN
# if defined PERTURBATION
        WRITE (stdout,10) Nrun
 10     FORMAT (/,' <<<< Ensemble/Perturbation Run: ',i5.5,' >>>>',/)
# elif defined W4DVAR
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
# endif
        WRITE (stdout,20) 'RP_INITIAL: Configuring and ',               &
     &                    'initializing representer model ...'
  20    FORMAT (/,1x,a,a,/)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize time stepping indices and counters.
!-----------------------------------------------------------------------
!
      iif(ng)=1
      indx1(ng)=1
      kstp(ng)=1
      krhs(ng)=1
      knew(ng)=1
      PREDICTOR_2D_STEP(ng)=.FALSE.
      synchro_flag(ng)=.TRUE.
      first_time(ng)=0
!
      iic(ng)=0
      nstp(ng)=1
      nrhs(ng)=1
      nnew(ng)=1
# ifdef FLOATS_NOT_YET
      nf(ng)=0
      nfp1(ng)=1
      nfm1(ng)=4
      nfm2(ng)=3
      nfm3(ng)=2
# endif
      tdays(ng)=dstart
      time(ng)=tdays(ng)*day2sec
      ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
      ntend(ng)=ntimes(ng)
      ntfirst(ng)=ntstart(ng)
      CALL time_string (time(ng), time_code(ng))

      IniRec=nrrec(ng)
      Tindex=1

      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        CALL wclock_on (ng, iRPM, 1)
      END DO
!$OMP END PARALLEL DO
# endif

# if defined WEAK_CONSTRAINT || defined IOM
!
!-----------------------------------------------------------------------
!  If weak constraint variational data assimilation, reset several IO
!  switches and variables.
!-----------------------------------------------------------------------
!
!  Set switch to create (TRUE) representer model initial conditions
!  NetCDF file or append (FALSE) to existing NetCDF files.
!
      IF (Nrun.eq.ERstr) THEN
#  ifdef ANA_INITIAL
        LdefIRP(ng)=.TRUE.
#  endif
        CALL rp_def_ini (ng)
        IF (exit_flag.ne.NoError) RETURN
      ELSE
      END IF
      IniRec=tIRPindx(ng)

#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Initialize surface forcing counter for storage arrays.
!
      SFcount(ng)=0
#  endif
!
!  Reset representer model history time record counters. These
!  counters are reset in every iteration pass but the NetCDF is
!  created on the first iteration pass.
!
      LcycleTLM(ng)=.FALSE.
      tTLMindx(ng)=0
      NrecTLM(ng)=0
!
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations.
!
      CALL obs_initial (ng, iRPM, .FALSE.)
# endif
!
!=======================================================================
!  On first pass of ensemble run loop, initialize model configuration.
!=======================================================================
!
      IF (Nrun.eq.ERstr) THEN
!
!-----------------------------------------------------------------------
!  Set horizontal grid, bathymetry, and Land/Sea masking (if any).
!  Use analytical functions or read in from a grid NetCDF.
!-----------------------------------------------------------------------
!
# ifdef ANA_GRID
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_grid (ng, TILE, iRPM)
#  ifdef MASKING
            CALL ana_mask (ng, TILE, iRPM)
#  endif
          END DO
        END DO
!$OMP END PARALLEL DO
# else
        CALL get_grid (ng, iRPM)
        if (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set vertical S-coordinate transformation function.
!-----------------------------------------------------------------------
!
        CALL set_scoord (ng)
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
        CALL set_weights (ng)
# endif

# if defined VISC_GRID || defined DIFF_GRID || defined SPONGE
!
!-----------------------------------------------------------------------
!  Set horizontal mixing coefficients. Rescale according to the local
!  grid size. If applicable, increases horizontal mixing in sponge
!  areas.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_hmixcoef (ng, TILE, iRPM)
          END DO
        END DO
!$OMP END PARALLEL DO
# endif

# ifdef NUDGING_COFF
!
!-----------------------------------------------------------------------
!  If appropriate, set nudging coefficiests time scales.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_nudgcoef (ng, TILE, iRPM)
          END DO
        END DO
!$OMP END PARALLEL DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL metrics (ng, TILE, iRPM)
# if defined PROPAGATOR || \
    (defined MASKING && (defined READ_WATER || defined WRITE_WATER ))
          CALL wpoints (ng, TILE, iRPM)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration pass.
!=======================================================================

# if defined PICARD_TEST || defined WEAK_CONSTRAINT || defined IOM
!
!  Clear nonlinear (background) and tangent linear state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL initialize_ocean (ng, TILE, iNLM)
          CALL initialize_ocean (ng, TILE, iRPM)
#  ifdef SOLVE3D
          CALL initialize_coupling (ng, TILE, 0)
          CALL initialize_mixing (ng, TILE, iRPM)
#  endif
          CALL initialize_forces (ng, TILE, 0)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# if defined SOLVE3D && !defined INI_FILE
!
!-----------------------------------------------------------------------
!  If analytical initial conditions, compute initial time-evolving
!  depths with zero free-surface.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL set_depth (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# if !(defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET)
!
!-----------------------------------------------------------------------
!  Initialize tangent linear bathymetry tl_h(i,j) to h(i,j) so some of
!  the terms are cancelled in the barotropic pressure gradient.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs)                                  &
!$OMP&            SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL rp_bath (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.  Use analytical
!  functions or read from an initial or restart NetCDF file.
!-----------------------------------------------------------------------
!
# ifdef ANA_INITIAL
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_initial (ng, TILE, iRPM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_passive (ng, TILE, iRPM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology tracers.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_biology (ng, TILE, iRPM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# if defined ANA_SEDIMENT_NOT_YET && defined SOLVE3D
!
!  Analytical initial conditions for sediment tracers.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_sediment (ng, TILE, iRPM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif
!
!  Read in representer model initial conditions.
!
# ifdef INI_FILE
      CALL get_state (ng, iRPM, 1, IRPname(ng), IniRec, Tindex)
# else
      IF (nrrec(ng).ne.0) THEN
        CALL get_state (ng, iRPM, 1, IRPname(ng), IniRec, Tindex)
      END IF
# endif
      IF (exit_flag.ne.NoError) RETURN

# if defined ANA_PERTURB && defined SANITY_CHECK
!
!-----------------------------------------------------------------------
!  Perturb tangent linear initial conditions with analitical
!  expressions.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL ana_perturb (ng, TILE, iRPM)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute initial time-evolving depths.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL set_depth (ng, TILE)
          CALL rp_set_depth (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL rp_set_massflux (ng, TILE)
          CALL set_massflux (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Compute initial representer tangent linear and basic state
!  S-coordinates vertical velocity.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL rp_omega (ng, TILE)
          CALL omega (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
# ifdef TIMELESS_DATA
      CALL rp_get_idata (ng)
# endif
      CALL rp_get_data (ng)
      IF (exit_flag.ne.NoError) RETURN

# if defined WEAK_CONSTRAINT_NOT || defined IOM_NOT
!
!-----------------------------------------------------------------------
!  Read in convolved adjoint impulse forcing (first record) and its
!  application time.
!-----------------------------------------------------------------------
!
      IF (outer.ge.1) THEN
        FrcRec(ng)=1
        CALL get_state (ng, 7, 7, TLFname(ng), FrcRec(ng), 1)
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iTLM, exit_flag, 1)
#  endif
        IF (exit_flag.ne.NoError) RETURN
      END IF
# endif

# ifdef SOLVE3D
!-----------------------------------------------------------------------
!  Compute initial representer tangent linear and basic state equation
!  of state related quantities.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL rp_rho_eos (ng, TILE)
          CALL rho_eos (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      IF (Lstiffness) THEN
        Lstiffness=.FALSE.
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL stiffness (ng, TILE, iRPM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF

# if defined FLOATS_NOT_YET || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
      CALL grid_coords (ng, iRPM)
# endif

# ifdef WAVES_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled wave model.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL waves_coupling (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initiialization time wall clock.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        CALL wclock_off (ng, iRPM, 1)
      END DO
!$OMP END PARALLEL DO
# endif
      RETURN
      END SUBROUTINE rp_initial
#else
      SUBROUTINE rp_initial
      RETURN
      END SUBROUTINE rp_initial
#endif
