#include "cppdefs.h"
#ifdef TL_IOMS
      SUBROUTINE rp_initial (ng)
!
!svn $Id: rp_initial.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
      Use mod_boundary
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
      USE dateclock_mod,       ONLY : time_string
# ifdef DISTRIBUTE
      USE distribute_mod,      ONLY : mp_bcasti
# endif
# ifdef TLM_CHECK
      USE ini_adjust_mod,      ONLY : tl_ini_perturb
# endif
      USE ini_hmixcoef_mod,    ONLY : ini_hmixcoef
# if !(defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET)
      USE rp_set_depth_mod,    ONLY : rp_bath
# endif
# ifdef SOLVE3D
      USE rp_set_depth_mod,    ONLY : rp_set_depth
      USE rp_omega_mod,        ONLY : rp_omega
      USE rp_rho_eos_mod,      ONLY : rp_rho_eos
      USE rp_set_massflux_mod, ONLY : rp_set_massflux
      USE set_depth_mod,       ONLY : set_depth
      USE omega_mod,           ONLY : omega
      USE rho_eos_mod,         ONLY : rho_eos
      USE set_massflux_mod,    ONLY : set_massflux
# endif
# ifdef MASKING
      USE set_masks_mod,       ONLY : set_masks
# endif
      USE stiffness_mod,       ONLY : stiffness
      USE strings_mod,         ONLY : FoundError
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod,   ONLY : ocn2wav_coupling
# endif
# ifdef WET_DRY
      USE wetdry_mod,          ONLY : wetdry
# endif
# if defined PROPAGATOR || \
    (defined MASKING    && (defined READ_WATER || defined WRITE_WATER))
      USE wpoints_mod,         ONLY : wpoints
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
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, IniRec, Tindex
      integer :: thread, tile
!
      real(dp) :: my_dstart
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
      IF (Master) THEN
# if defined PERTURBATION
        WRITE (stdout,10) Nrun
 10     FORMAT (/,' <<<< Ensemble/Perturbation Run: ',i5.5,' >>>>',/)
# elif defined R4DVAR    || defined R4DVAR_ANA_SENSITIVITY || \
       defined TL_R4DVAR
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
# endif
        WRITE (stdout,20) 'RP_INITIAL: Configuring and ',               &
     &                    'initializing representer model ...'
 20     FORMAT (/,1x,a,a,/)
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

# ifdef GENERIC_DSTART
!
!  Rarely, the tangent linear model is initialized from a NetCDF file,
!  so we do not know its actual initialization time. Usually, it is
!  computed from DSTART, implying that its value is correct in the ROMS
!  input script. Therefore, the user needs to check and update its value
!  to every time that ROMS is executed. Alternatively, if available, we
!  can use the initialization time from the nonlinear model, INItime.
!  This variable is assigned when computing or processing the basic
!  state trajectory needed to linearize the adjoint model.
!
      IF (INItime(ng).lt.0.0_dp) THEN
        my_dstart=dstart                    ! ROMS input script
      ELSE
        my_dstart=INItime(ng)/86400.0_dp    ! NLM IC time is known
      END IF
      tdays(ng)=my_dstart
# else
      tdays(ng)=dstart
# endif
      time(ng)=tdays(ng)*day2sec
      ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
      ntend(ng)=ntstart(ng)+ntimes(ng)-1
      ntfirst(ng)=ntstart(ng)

      CALL time_string (time(ng), time_code(ng))

      IniRec=nrrec(ng)
      Tindex=1

      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Initialize global diagnostics variables.
!
      avgke=0.0_dp
      avgpe=0.0_dp
      avgkp=0.0_dp
      volume=0.0_dp

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
      DO thread=THREAD_RANGE
        CALL wclock_on (ng, iRPM, 2, __LINE__, MyFile)
      END DO
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
      IF (Nrun.eq.ERstr) THEN
#  ifdef ANA_INITIAL
        LdefIRP(ng)=.TRUE.
#  endif
        CALL rp_def_ini (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
# endif
!
!-----------------------------------------------------------------------
!  Set application grid, metrics, and associated variables and
!  parameters.
!-----------------------------------------------------------------------
!
      IF (SetGridConfig(ng)) THEN
        CALL set_grid (ng, iTLM)
        SetGridConfig(ng)=.FALSE.
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
      END IF
# endif
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
!!      CALL initialize_boundary (ng, tile, iRPM)
        CALL initialize_forces (ng, tile, iADM)
        CALL initialize_forces (ng, tile, iNLM)
        CALL initialize_forces (ng, tile, iRPM)
      END DO
# endif

# if defined SOLVE3D && !defined INI_FILE
!
!-----------------------------------------------------------------------
!  If analytical initial conditions, compute initial time-evolving
!  depths with zero free-surface.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_depth (ng, tile, iRPM)
      END DO
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
# endif
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.  Use analytical
!  functions or read from an initial or restart NetCDF file.
!-----------------------------------------------------------------------

# ifdef ANA_INITIAL
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_initial (ng, tile, iRPM)
        END DO
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
      END IF
# endif
!
!  Read in representer model initial conditions.
!
# ifdef INI_FILE
      CALL get_state (ng, iRPM, 1, IRP(ng)%name, IniRec, Tindex)
      time(ng)=io_time                       ! needed for shared-memory
# else
      IF (nrrec(ng).ne.0) THEN
        CALL get_state (ng, iRPM, 1, IRP(ng)name, IniRec, Tindex)
        time(ng)=io_time                     ! needed for shared-memory
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
      END IF
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef WET_DRY
!
!-----------------------------------------------------------------------
!  Process initial wet/dry masks.
!-----------------------------------------------------------------------
!
!  If restart, read in wet/dry masks.
!
      IF (nrrec(ng).ne.0) THEN
        CALL get_wetdry (ng, iRPM, IniRec(ng))
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iRPM, exit_flag)
#  endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
        END DO
      END IF
# endif

# ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  Open observations NetCDF file and initialize various variables
!  needed for processing the tangent linear state solution at
!  observation locations. Need to be done after processing initial
!  conditions since the correct initial time is needed to determine
!  the first "ObsTime" to process.
!-----------------------------------------------------------------------
!
      CALL obs_initial (ng, iRPM, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute initial time-evolving depths.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_depth (ng, tile, iRPM)
        CALL rp_set_depth (ng, tile, iRPM)
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL rp_set_massflux (ng, tile, iRPM)
        CALL set_massflux (ng, tile, iRPM)
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial representer tangent linear and basic state
!  S-coordinates vertical velocity.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL rp_omega (ng, tile, iRPM)
        CALL omega (ng, tile, iRPM)
      END DO
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
      CALL close_inp (ng, iRPM)
      CALL check_multifile (ng, iRPM)
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iRPM, exit_flag)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
      CALL rp_get_idata (ng)
      CALL rp_get_data (ng)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_masks (ng, tile, iRPM)
      END DO
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
# endif

# ifdef SOLVE3D
!-----------------------------------------------------------------------
!  Compute initial representer tangent linear and basic state equation
!  of state related quantities.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL rp_rho_eos (ng, tile, iRPM)
        CALL rho_eos (ng, tile, iRPM)
      END DO
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

# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled wave model.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ocn2wav_coupling (ng, tile)
      END DO
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
      time(ng)=time(ng)-dt(ng)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initiialization time wall clock.
!-----------------------------------------------------------------------
!
      DO thread=THREAD_RANGE
        CALL wclock_off (ng, iRPM, 2, __LINE__, MyFile)
      END DO
# endif
!
      RETURN
      END SUBROUTINE rp_initial
#else
      SUBROUTINE rp_initial
      RETURN
      END SUBROUTINE rp_initial
#endif
