#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_initial (ng)
!
!svn $Id: ad_initial.F 995 2020-01-10 04:01:28Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes all adjoint model variables.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BBL_MODEL_NOT_YET
      USE mod_bbl
# endif
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_grid
      USE mod_iounits
# ifdef SOLVE3D
      USE mod_mixing
# endif
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE analytical_mod
      USE dateclock_mod,     ONLY : time_string
# ifdef DISTRIBUTE
      USE distribute_mod,    ONLY : mp_bcasti
# endif
# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
      USE adsen_initial_mod, ONLY : adsen_initial
# endif
      USE ini_hmixcoef_mod,  ONLY : ini_hmixcoef
      USE metrics_mod,       ONLY : metrics
# ifdef SOLVE3D
      USE set_depth_mod,     ONLY : set_depth
      USE omega_mod,         ONLY : omega
      USE rho_eos_mod,       ONLY : rho_eos
      USE set_massflux_mod,  ONLY : set_massflux
# endif
# ifdef MASKING
      USE set_masks_mod,     ONLY : set_masks
# endif
      USE stiffness_mod,     ONLY : stiffness
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod, ONLY : ocn2wav_coupling
# endif
      USE strings_mod,       ONLY : FoundError
# ifdef WET_DRY
      USE wetdry_mod,        ONLY : wetdry
# endif
# if defined PROPAGATOR || \
    (defined MASKING    && (defined READ_WATER || defined WRITE_WATER))
      USE wpoints_mod,       ONLY : wpoints
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

      real(dp) :: my_dstart
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
# elif defined IS4DVAR || defined SENSITIVITY_4DVAR || \
       defined W4DPSAS || defined W4DVAR
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
# endif
        WRITE (stdout,20) 'AD_INITIAL: Configuring and ',               &
     &                    'initializing adjoint model ...'
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
      krhs(ng)=3
      knew(ng)=2
      PREDICTOR_2D_STEP(ng)=.FALSE.
!
      iic(ng)=0
# ifdef SOLVE3D
      nstp(ng)=1
      nnew(ng)=2
      nrhs(ng)=nstp(ng)
# endif
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
      ad_ubar_xs=0.0_r8

# ifdef GENERIC_DSTART
!
!  Rarely, the adjoint model is initialized from a NetCDF file, so we do
!  not know its actual initialization time. Usually, it is computed from
!  DSTART, implying that its value is correct in the ROMS input script.
!  Therefore, the user needs to check and update its value to every time
!  that ROMS is executed. Alternatively, if available, we can use the
!  initialization time from the nonlinear model, INItime. This variable
!  is assigned when computing or processing the basic state trajectory
!  needed to linearize the adjoint model.
!
      IF (INItime(ng).lt.0.0_dp) THEN
        my_dstart=dstart                    ! ROMS input script
      ELSE
        my_dstart=INItime(ng)/86400.0_dp    ! NLM IC time is known
      END IF
      tdays(ng)=my_dstart+dt(ng)*REAL(ntimes(ng),r8)*sec2day
      time(ng)=tdays(ng)*day2sec
!$OMP MASTER
      ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
      ntend(ng)=ntstart(ng)-ntimes(ng)-1
      ntfirst(ng)=ntend(ng)
!$OMP END MASTER
!$OMP BARRIER
#else
      tdays(ng)=dstart+dt(ng)*REAL(ntimes(ng)-ntfirst(ng)+1,r8)*sec2day
      time(ng)=tdays(ng)*day2sec
!$OMP MASTER
      ntstart(ng)=ntimes(ng)+1
      ntend(ng)=ntfirst(ng)
      ntfirst(ng)=ntend(ng)
!$OMP END MASTER
!$OMP BARRIER
#endif

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
        CALL wclock_on (ng, iADM, 2, __LINE__, __FILE__)
      END DO
!$OMP BARRIER
# endif

# if defined FOUR_DVAR            && \
    !(defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS)
!
!-----------------------------------------------------------------------
!  If variational data assimilation, reset several IO switches and
!  variables.
!-----------------------------------------------------------------------
!
!  Set switch to create adjoint NetCDF file or append to an existing
!  adjoint NetCDF file.
!
!$OMP MASTER
      IF (Nrun.eq.ERstr) THEN
        LdefADJ(ng)=.TRUE.
      END IF
!
!  Activate switch to write adjoint NetCDF file.
!
      LwrtADJ(ng)=.TRUE.

#  ifdef ADJUST_BOUNDARY
!
!  Initialize open boundary counter for storage arrays.
!
      OBCcount(ng)=0
#  endif
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Initialize surface forcing counter for storage arrays.
!
      SFcount(ng)=Nfrec(ng)+1
#  endif
!$OMP END MASTER
!$OMP BARRIER
# endif
!
!=======================================================================
!  On first pass of ensemble/perturbation/iteration loop, initialize
!  model configuration.
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
          CALL ana_grid (ng, tile, iADM)
#  ifdef MASKING
          CALL ana_mask (ng, tile, iADM)
#  endif
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI             || \
      defined SENSITIVITY_4DVAR
          CALL ana_scope (ng, tile, iADM)
#  endif
        END DO
!$OMP BARRIER
# else
!$OMP MASTER
        CALL get_grid (ng, iADM)
!$OMP END MASTER
!$OMP BARRIER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, exit_flag)
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
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
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL metrics (ng, tile, iADM)
        END DO
!$OMP BARRIER

# if defined WTYPE_GRID && defined ANA_WTYPE     && \
    (defined LMD_SKPP   || defined SOLAR_SOURCE)
!
!-----------------------------------------------------------------------
!  Set spatially varying Jerlov water type.
!-----------------------------------------------------------------------
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_wtype (ng, tile, iADM)
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
            CALL ana_nudgcoef (ng, tile, iADM)
          END DO
!$OMP BARRIER
        END IF
# else
        IF (Lnudging(ng)) THEN
!$OMP MASTER
          CALL get_nudgcoef (ng, iADM)
!$OMP END MASTER
#  ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iADM, exit_flag)
#  endif
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
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
        CALL ini_hmixcoef (ng, tile, iADM)
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
          CALL ana_sponge (ng, tile, iADM)
        END DO
!$OMP BARRIER
      END IF
# endif
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration pass.
!=======================================================================

# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
!
!-----------------------------------------------------------------------
!  Clear all adjoint and nonlinear variables.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL initialize_ocean (ng, tile, iNLM)
        CALL initialize_ocean (ng, tile, iADM)
#  if defined SOLVE3D
        CALL initialize_coupling (ng, tile, iADM)
        CALL initialize_mixing (ng, tile, iADM)
#  endif
        CALL initialize_forces (ng, tile, iADM)
        CALL initialize_grid (ng, tile, iADM)
#  ifdef ADJUST_BOUNDARY
        CALL initialize_boundary (ng, tile, iADM)
#  endif
      END DO
!$OMF BARRIER

# elif defined FOUR_DVAR && !defined IS4DVAR_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Clear all adjoint variables.  In variational data assimilation the
!  initial condition are always zero and the forcing is only via the
!  (model-observations) misfit terms.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL initialize_ocean (ng, tile, iNLM)
        CALL initialize_ocean (ng, tile, iADM)
#  if defined SOLVE3D
        CALL initialize_coupling (ng, tile, iADM)
        CALL initialize_mixing (ng, tile, iADM)
#  endif
        CALL initialize_forces (ng, tile, iADM)
        CALL initialize_grid (ng, tile, iADM)
#  ifdef ADJUST_BOUNDARY
        CALL initialize_boundary (ng, tile, iADM)
#  endif
      END DO
!$OMF BARRIER

# else

#  if defined SOLVE3D && !defined INI_FILE
!
!-----------------------------------------------------------------------
!  If analytical initial conditions, compute initial time-evolving
!  depths with zero free-surface.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_depth (ng, tile, iADM)
      END DO
!$OMP BARRIER
#  endif
!
!-----------------------------------------------------------------------
!  Set adjoint primitive variables initial conditions.
!-----------------------------------------------------------------------
!
#  ifdef ANA_INITIAL
!
!  Analytical initial conditions for momentum and active tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_initial (ng, tile, iADM)
        END DO
!$OMP BARRIER
      END IF
#  endif

#  if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_passive (ng, tile, iADM)
        END DO
!$OMP BARRIER
      END IF
#  endif

#  if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_biology (ng, tile, iADM)
        END DO
!$OMP BARRIER
      END IF
#  endif

#  if defined ANA_SEDIMENT_NOT_YET && defined SOLVE3D
!
!  Analytical initial conditions for sediment tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_sediment (ng, tile, iADM)
        END DO
!$OMP BARRIER
      END IF
#  endif
!
!  Read in initial conditions for initial or restart NetCDF file.
!
#  ifdef INI_FILE
!$OMP MASTER
      CALL get_state (ng, iADM, 1, IAD(ng)%name, IniRec, Tindex)
!$OMP END MASTER
!$OMP BARRIER
      time(ng)=io_time                       ! needed for shared-memory
#   ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iADM, exit_flag)
#   endif
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
#  else
      IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
        CALL get_state (ng, iADM, 1, IAD(ng)%name, IniRec, Tindex)
!$OMP END MASTER
!$OMP BARRIER
        time(ng)=io_time                     ! needed for shared-memory
#   ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, exit_flag)
#   endif
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END IF
#  endif
# endif

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
        CALL get_wetdry (ng, iADM, IniRec(ng))
!$OMP END MASTER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, exit_flag)
#  endif
!$OMP BARRIER
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      ELSE
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
        END DO
!$OMP BARRIER
      END IF
# endif

# ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  Initialize various variables needed for processing observations
!  backwards in time. Need to be done after processing initial
!  conditions since the correct initial time is needed to determine
!  the first "ObsTime" to process.
!-----------------------------------------------------------------------
!
!$OMP MASTER
# ifdef W4DPSAS_FCT_SENSITIVITY
      IF (.not.LsenPSAS(ng).and.LsenFCT(ng)) THEN
#  ifdef OBS_SPACE
        IF (Lobspace(ng)) THEN
#  endif
# endif
          CALL obs_initial (ng, iADM, .TRUE.)
# ifdef OBS_SPACE
        END IF
# endif
# ifdef W4DPSAS_FCT_SENSITIVITY
      END IF
# endif
!$OMP END MASTER
!$OMP BARRIER
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
# endif

# if defined ANA_PERTURB  && \
    (defined SANITY_CHECK || defined R_SYMMETRY)
!
!-----------------------------------------------------------------------
!  Perturb adjoint initial conditions with analitical expressions.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ana_perturb (ng, tile, iADM)
      END DO
!$OMP BARRIER
# endif

# ifdef SOLVE3D
!!
!!----------------------------------------------------------------------
!!  Compute initial time-evolving depths.
!!----------------------------------------------------------------------
!!
!!    DO tile=first_tile(ng),last_tile(ng),+1
!!      CALL ad_set_depth (ng, tile, iADM)
!!    END DO
!!$OMP BARRIER
!!
!!----------------------------------------------------------------------
!!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!!----------------------------------------------------------------------
!!
!!    DO tile=first_tile(ng),last_tile(ng),+1
!!      CALL ad_set_massflux (ng, tile, iADM)
!!    END DO
!!$OMP BARRIER
!!
!!----------------------------------------------------------------------
!!  Compute initial S-coordinates vertical velocity. Compute initial
!!  density anomaly from potential temperature and salinity via equation
!!  of state for seawater.  Also compute other equation of state related
!!  quatities.
!!----------------------------------------------------------------------
!!
!!    DO tile=first_tile(ng),last_tile(ng),+1
!!      CALL ad_omega (ng, tile, iADM)
!!      CALL ad_rho_eos (ng, tile, iADM)
!!    END DO
!!!$OMP BARRIER
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
          CALL ana_psource (ng, tile, iADM)
        END DO
      END IF
!$OMP BARRIER
#endif
!
!-----------------------------------------------------------------------
!  If applicable, close all input boundary, climatology, and forcing
!  NetCDF files and set associated parameters to the closed state. This
!  step is essential in iterative algorithms that run the full TLM
!  repetitively. Then, Initialize several parameters in their file
!  structure, so the appropriate input single or multi-file is selected
!  during initialization/restart.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL close_inp (ng, iADM)
      CALL check_multifile (ng, iADM)
!$OMP END MASTER
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iADM, exit_flag)
# endif
!$OMP BARRIER
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
!$OMP MASTER
      CALL ad_get_idata (ng)
      CALL ad_get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iADM, exit_flag)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

# ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_masks (ng, tile, iADM)
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
        CALL wpoints (ng, tile, iADM)
      END DO
!$OMP BARRIER
# endif

# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
!
!-----------------------------------------------------------------------
!  Initialize adjoint state with the functional whose sensitivity is
!  is required.
!-----------------------------------------------------------------------
!
#  ifdef SENSITIVITY_4DVAR
      IF (LsenPSAS(ng)) THEN
#  endif
#  if !defined AD_IMPULSE
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL adsen_initial (ng, tile)
        END DO
!$OMP BARRIER
#  endif
#  ifdef SENSITIVITY_4DVAR
      END IF
#  endif
# endif

# if defined ANA_DRAG && defined UV_DRAG_GRID
!
!-----------------------------------------------------------------------
!  Set analytical spatially varying bottom friction parameter.
!-----------------------------------------------------------------------
!
      IF (Nrun.eq.ERstr) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_drag (ng, tile, iADM)
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
          CALL stiffness (ng, tile, iADM)
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
      CALL grid_coords (ng, iADM)
!$OMP END MASTER
!$OMP BARRIER
# endif

# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled wave model.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ocn2wav_coupling (ng, tile)
      END DO
!$OMP BARRIER
# endif
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter and clock.
!-----------------------------------------------------------------------
!
!  Add one time unit to avoid special case due to initialization
!  in the main time-stepping routine.
!
      iic(ng)=ntstart(ng)+1
!$OMP MASTER
      time(ng)=time(ng)+dt(ng)
!$OMP END MASTER
!$OMP BARRIER

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initiialization time wall clock.
!-----------------------------------------------------------------------
!
      DO thread=THREAD_RANGE
        CALL wclock_off (ng, iADM, 2, __LINE__, __FILE__)
      END DO
!$OMP BARRIER
# endif
      RETURN
      END SUBROUTINE ad_initial
#else
      SUBROUTINE ad_initial
      RETURN
      END SUBROUTINE ad_initial
#endif
