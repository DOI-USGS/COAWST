#include "cppdefs.h"
#ifdef TL_IOMS
      SUBROUTINE rp_initial (ng)
!
!svn $Id: rp_initial.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      USE analytical_mod
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
# ifdef TLM_CHECK
      USE ini_adjust_mod, ONLY : tl_ini_perturb
# endif
      USE ini_hmixcoef_mod, ONLY : ini_hmixcoef
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
# ifdef MASKING
      USE set_masks_mod, ONLY : set_masks
# endif
      USE stiffness_mod, ONLY : stiffness
# ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : waves_coupling
# endif
# ifdef WET_DRY
      USE wetdry_mod, ONLY : wetdry
# endif
# if defined PROPAGATOR || \
    (defined MASKING    && (defined READ_WATER || defined WRITE_WATER))
      USE wpoints_mod, ONLY : wpoints
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
      integer :: Fcount, IniRec, Tindex
      integer :: thread, tile
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
!$OMP MASTER
      IF (Master) THEN
# if defined PERTURBATION
        WRITE (stdout,10) Nrun
 10     FORMAT (/,' <<<< Ensemble/Perturbation Run: ',i5.5,' >>>>',/)
# elif defined TL_W4DVAR || defined W4DVAR || defined W4DVAR_SENSITIVITY
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
# endif
        WRITE (stdout,20) 'RP_INITIAL: Configuring and ',               &
     &                    'initializing representer model ...'
 20     FORMAT (/,1x,a,a,/)
      END IF
!$OMP END MASTER
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
!
      synchro_flag(ng)=.TRUE.
      first_time(ng)=0
      tdays(ng)=dstart
      time(ng)=tdays(ng)*day2sec
!$OMP MASTER
      ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
      ntend(ng)=ntimes(ng)
      ntfirst(ng)=ntstart(ng)
!$OMP END MASTER
!$OMP BARRIER

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
      DO thread=THREAD_RANGE
        CALL wclock_on (ng, iRPM, 2)
      END DO
!$OMP BARRIER
# endif

# ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  If weak constraint variational data assimilation, reset several IO
!  switches and variables.
!-----------------------------------------------------------------------
!
!  Set switch to create (TRUE) representer model initial conditions
!  NetCDF file or append (FALSE) to existing NetCDF files.
!
!$OMP MASTER
      IF (Nrun.eq.ERstr) THEN
#  ifdef ANA_INITIAL
        LdefIRP(ng)=.TRUE.
#  endif
        CALL rp_def_ini (ng)
        IF (exit_flag.ne.NoError) RETURN
      ELSE
      END IF
      IniRec=IRP(ng)%Rindex

# ifdef ADJUST_BOUNDARY
!
!  Initialize open boundary counter for storage arrays.
!
      OBCcount(ng)=0
# endif
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
      TLM(ng)%Rindex=0
      Fcount=TLM(ng)%Fcount
      TLM(ng)%Nrec(Fcount)=0
!
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations.
!
      CALL obs_initial (ng, iRPM, .FALSE.)
!$OMP END MASTER
!$OMP BARRIER
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
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_grid (ng, tile, iRPM)
#  ifdef MASKING
          CALL ana_mask (ng, tile, iRPM)
#  endif
        END DO
!$OMP BARRIER
# else
!$OMP MASTER
        CALL get_grid (ng, iRPM)
!$OMP END MASTER
!$OMP BARRIER
        if (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set vertical S-coordinate transformation function.
!-----------------------------------------------------------------------
!
!$OMP MASTER
        CALL set_scoord (ng)
!$OMP END MASTER
!$OMP BARRIER
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
!$OMP MASTER
        CALL set_weights (ng)
!$OMP END MASTER
!$OMP BARRIER
# endif

# if defined WTYPE_GRID && defined ANA_WTYPE     && \
    (defined LMD_SKPP   || defined SOLAR_SOURCE)
!
!-----------------------------------------------------------------------
!  Set spatially varying Jerlov water type.
!-----------------------------------------------------------------------
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_wtype (ng, tile, iRPM)
        END DO
!$OMP BARRIER
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, set spatially varying nudging coefficients time
!  scales.
!-----------------------------------------------------------------------
!
# ifdef ANA_NUDGCOEF
        IF (Lnudging(ng)) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_nudgcoef (ng, tile, iRPM)
          END DO
!$OMP BARRIER
        END IF
# else
        IF (Lnudging(ng)) THEN
!$OMP MASTER
          CALL get_nudgcoef (ng, iRPM)
!$OMP END MASTER
#  ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iRPM, exit_flag)
#  endif
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients. If applicable, scale
!  mixing coefficients according to the grid size (smallest area).
# ifndef ANA_SPONGE
!  Also increase their values in sponge areas using the "visc_factor"
!  and/or "diff_factor" read from input Grid NetCDF file.
# endif
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ini_hmixcoef (ng, tile, iRPM)
      END DO
!$OMP BARRIER

# ifdef ANA_SPONGE
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing coefficients in sponge areas using
!  analytical functions.
!-----------------------------------------------------------------------
!
      IF (Lsponge(ng)) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_sponge (ng, tile, iRPM)
        END DO
!$OMP BARRIER
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL metrics (ng, tile, iRPM)
      END DO
!$OMP BARRIER
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration pass.
!=======================================================================

# if defined PICARD_TEST || defined WEAK_CONSTRAINT
!
!  Clear nonlinear (background) and tangent linear state variables.
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL initialize_ocean (ng, tile, iNLM)
        CALL initialize_ocean (ng, tile, iRPM)
#  ifdef SOLVE3D
        CALL initialize_coupling (ng, tile, 0)
        CALL initialize_mixing (ng, tile, iRPM)
#  endif
        CALL initialize_forces (ng, tile, iRPM)
        CALL initialize_forces (ng, tile, iTLM)
        CALL initialize_forces (ng, tile, iNLM)
        CALL initialize_forces (ng, tile, iADM)
      END DO
!$OMP BARRIER
# endif

# if defined SOLVE3D && !defined INI_FILE
!
!-----------------------------------------------------------------------
!  If analytical initial conditions, compute initial time-evolving
!  depths with zero free-surface.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_depth (ng, tile)
      END DO
!$OMP BARRIER
# endif

# if !(defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET)
!
!-----------------------------------------------------------------------
!  Initialize tangent linear bathymetry tl_h(i,j) to h(i,j) so some of
!  the terms are cancelled in the barotropic pressure gradient.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL rp_bath (ng, tile)
      END DO
!$OMP BARRIER
# endif
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.  Use analytical
!  functions or read from an initial or restart NetCDF file.
!-----------------------------------------------------------------------
!
# ifdef ANA_INITIAL
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_initial (ng, tile, iRPM)
        END DO
!$OMP BARRIER
      END IF
# endif

# if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_passive (ng, tile, iRPM)
        END DO
!$OMP BARRIER
      END IF
# endif

# if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_biology (ng, tile, iRPM)
        END DO
!$OMP BARRIER
      END IF
# endif

# if defined ANA_SEDIMENT_NOT_YET && defined SOLVE3D
!
!  Analytical initial conditions for sediment tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_sediment (ng, tile, iRPM)
        END DO
!$OMP BARRIER
      END IF
# endif
!
!  Read in representer model initial conditions.
!
# ifdef INI_FILE
!$OMP MASTER
      CALL get_state (ng, iRPM, 1, IRP(ng)%name, IniRec, Tindex)
!$OMP END MASTER
!$OMP BARRIER
      time(ng)=io_time                       ! needed for shared-memory
# else
      IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
        CALL get_state (ng, iRPM, 1, IRP(ng)name, IniRec, Tindex)
!$OMP END MASTER
!$OMP BARRIER
        time(ng)=io_time                     ! needed for shared-memory
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
      END IF
# endif
      IF (exit_flag.ne.NoError) RETURN

# ifdef WET_DRY
!
!-----------------------------------------------------------------------
!  Process initial wet/dry masks.
!-----------------------------------------------------------------------
!
!  If restart, read in wet/dry masks.
!
      IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
        CALL get_wetdry (ng, iRPM, IniRec(ng))
!$OMP END MASTER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iRPM, exit_flag)
#  endif
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
      ELSE
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
        END DO
!$OMP BARRIER
      END IF
# endif

# if defined ANA_PERTURB && defined SANITY_CHECK
!
!-----------------------------------------------------------------------
!  Perturb tangent linear initial conditions with analitical
!  expressions.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ana_perturb (ng, tile, iRPM)
      END DO
!$OMP BARRIER
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute initial time-evolving depths.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_depth (ng, tile)
        CALL rp_set_depth (ng, tile)
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL rp_set_massflux (ng, tile)
        CALL set_massflux (ng, tile)
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Compute initial representer tangent linear and basic state
!  S-coordinates vertical velocity.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL rp_omega (ng, tile)
        CALL omega (ng, tile)
      END DO
!$OMP BARRIER
# endif

#ifdef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks position, direction, special flag, and mass
!  transport nondimensional shape profile with analytcal expressions.
!  Point sources are at U- and V-points. We need to get their positions
!  to process internal Land/Sea masking arrays during initialization.
!-----------------------------------------------------------------------
!
      IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_psource (ng, tile, iRPM)
        END DO
      END IF
!$OMP BARRIER
#endif
!
!-----------------------------------------------------------------------
!  If applicable, check input multi-files where the time records for
!  a particular field are split into several NetCDF files.  Initialize
!  several parameters in the file structure so the appropriate input
!  file is selected during initialization/restart.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL check_multifile (ng, iRPM)
!$OMP END MASTER
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iRPM, exit_flag)
# endif
!$OMP BARRIER
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL close_inp (ng, iRPM)
      CALL rp_get_idata (ng)
      CALL rp_get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
      IF (exit_flag.ne.NoError) RETURN

# ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_masks (ng, tile, iRPM)
      END DO
!$OMP BARRIER
# endif

# if defined PROPAGATOR || \
    (defined MASKING    && (defined READ_WATER || defined WRITE_WATER ))
!
!-----------------------------------------------------------------------
!  Set variables associated with the processing water points and/or
!  size of packed state arrays.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL wpoints (ng, tile, iRPM)
      END DO
!$OMP BARRIER
# endif

# ifdef SOLVE3D
!-----------------------------------------------------------------------
!  Compute initial representer tangent linear and basic state equation
!  of state related quantities.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL rp_rho_eos (ng, tile)
        CALL rho_eos (ng, tile)
      END DO
!$OMP BARRIER
# endif

# if defined ANA_DRAG && defined UV_DRAG_GRID
!
!-----------------------------------------------------------------------
!  Set analytical spatially varying bottom friction parameter.
!-----------------------------------------------------------------------
!
      IF (Nrun.eq.ERstr) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_drag (ng, tile, iRPM)
        END DO
!$OMP BARRIER
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      IF (Lstiffness) THEN
        Lstiffness=.FALSE.
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL stiffness (ng, tile, iRPM)
        END DO
!$OMP BARRIER
      END IF

# if defined FLOATS_NOT_YET || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL grid_coords (ng, iRPM)
!$OMP END MASTER
!$OMP BARRIER
# endif

# ifdef WAVES_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled wave model.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL waves_coupling (ng, tile)
      END DO
!$OMP BARRIER
# endif
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter and clock.
!-----------------------------------------------------------------------
!
!  Subsract one time unit to avoid special case due to initialization
!  in the main time-stepping routine.
!
      iic(ng)=ntstart(ng)-1
!$OMP MASTER
      time(ng)=time(ng)-dt(ng)
!$OMP END MASTER
!$OMP BARRIER

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initiialization time wall clock.
!-----------------------------------------------------------------------
!
      DO thread=THREAD_RANGE
        CALL wclock_off (ng, iRPM, 2)
      END DO
!$OMP BARRIER
# endif
      RETURN
      END SUBROUTINE rp_initial
#else
      SUBROUTINE rp_initial
      RETURN
      END SUBROUTINE rp_initial
#endif
