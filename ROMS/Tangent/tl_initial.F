#include "cppdefs.h"
#ifdef TANGENT
      SUBROUTINE tl_initial (ng)
!
!git $Id$
!svn $Id: tl_initial.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes all tangent linear model variables.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BBL_MODEL_NOT_YET
      USE mod_bbl
# endif
# ifdef FOUR_DVAR
#  ifdef SOLVE3D
      USE mod_coupling
#  endif
#  ifdef SP4DVAR
      USE mod_forces
#  endif
      USE mod_fourdvar
# endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
# ifdef ADJUST_BOUNDARY
      USE mod_boundary,        ONLY : initialize_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces,          ONLY : initialize_forces
# endif
!
      USE analytical_mod
      USE close_io_mod,        ONLY : close_inp
      USE dateclock_mod,       ONLY : time_string
# ifdef DISTRIBUTE
      USE distribute_mod,      ONLY : mp_bcasti
# endif
      USE get_state_mod,       ONLY : get_state
# ifdef WET_DRY
      USE get_wetdry_mod,      ONLY : get_wetdry
# endif
# ifdef TLM_CHECK
      USE ini_adjust_mod,      ONLY : tl_ini_perturb
# endif
      USE ini_hmixcoef_mod,    ONLY : ini_hmixcoef
# ifdef I4DVAR_ANA_SENSITIVITY
      USE ini_lanczos_mod,     ONLY : ini_lanczos
# endif
# ifdef OBSERVATIONS
      USE obs_initial_mod,     ONLY : obs_initial
# endif
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod,   ONLY : ocn2wav_coupling
# endif
# ifdef SOLVE3D
      USE set_depth_mod,       ONLY : set_depth
# endif
# ifdef MASKING
      USE set_masks_mod,       ONLY : set_masks
# endif
      USE stiffness_mod,       ONLY : stiffness
      USE strings_mod,         ONLY : FoundError
# if defined RBL4DVAR           || defined R4DVAR      || \
     defined SENSITIVITY_4DVAR  || defined TL_RBL4DVAR || \
     defined TL_R4DVAR
      USE tl_set_depth_mod,    ONLY : tl_bath
# endif
# ifdef FOUR_DVAR
      USE tl_def_ini_mod,      ONLY : tl_def_ini
# endif
# ifdef SOLVE3D
      USE tl_omega_mod,        ONLY : tl_omega
      USE tl_rho_eos_mod,      ONLY : tl_rho_eos
      USE tl_set_depth_mod,    ONLY : tl_set_depth
      USE tl_set_massflux_mod, ONLY : tl_set_massflux
# endif
# ifdef FOUR_DVAR
      USE tl_wrt_ini_mod,      ONLY : tl_wrt_ini
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
      integer :: Fcount, IniRec, Tindex, wrtRec
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
# elif defined I4DVAR      || defined RBL4DVAR          || \
       defined R4DVAR      || defined SENSITIVITY_4DVAR || \
       defined TL_RBL4DVAR || defined TL_R4DVAR
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
# endif
        WRITE (stdout,20) 'TL_INITIAL: Configuring and ',               &
     &                    'initializing tangent linear model ...'
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
      IF (ANY(tl_VolCons(:,ng))) THEN
        tl_ubar_xs=0.0_r8
      END IF

# if defined GENERIC_DSTART
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
      step_counter(ng)=0
!
      IniRec=nrrec(ng)
      Tindex=1
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
        CALL wclock_on (ng, iTLM, 2, __LINE__, MyFile)
      END DO
# endif

# if defined OPT_OBSERVATIONS
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
!  Set initial conditions time record to read.
!
       IniRec=1

# elif defined FOUR_DVAR    && \
     !(defined HESSIAN_SV   || defined HESSIAN_SO || \
       defined HESSIAN_FSV) || defined TLM_CHECK
!
!-----------------------------------------------------------------------
!  If variational data assimilation, reset several IO switches and
!  variables.
!-----------------------------------------------------------------------

#  ifndef I4DVAR_ANA_SENSITIVITY
#   ifdef I4DVAR
!
!  Set switch to create (TRUE) tangent linear initial conditions and
!  history NetCDF files or append (FALSE) to existing files. Then,
!  create tangent linear model initialization file and write zero
!  initial conditions for records 1 and 2.
!
      IF ((Nrun.eq.ERstr).and.(inner.eq.0)) THEN
        LdefITL(ng)=.TRUE.
        CALL tl_def_ini (ng)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    ifdef DISTRIBUTE
        CALL tl_wrt_ini (ng, MyRank, Tindex, 1)
#    else
        CALL tl_wrt_ini (ng, -1, Tindex, 1)
#    endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    ifdef DISTRIBUTE
        CALL tl_wrt_ini (ng, MyRank, Tindex, 2)
#    else
        CALL tl_wrt_ini (ng, -1, Tindex, 2)
#    endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   ifndef WEAK_CONSTRAINT
!
!  Set switch to create tangent linear history file.
!
      IF (Nrun.eq.ERstr) THEN
        LdefTLM(ng)=.TRUE.
      END IF
#   endif
!
!  Set record to read from initial tangent linear NetCDF file.
!
      IniRec=ITL(ng)%Rindex

#   ifdef ADJUST_BOUNDARY
!
!  Initialize open boundary counter for storage arrays.
!
      OBCcount(ng)=0
#   endif
#   if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Initialize surface forcing counter for storage arrays.
!
      SFcount(ng)=0
#   endif
#   if !defined WEAK_CONSTRAINT
!
!  Reset tangent linear model history time record counters. These
!  counters are reset in every iteration pass. This file is created
!  on the first iteration pass.
!
      TLM(ng)%Rindex=0
      Fcount=TLM(ng)%Fcount
      TLM(ng)%Nrec(Fcount)=0
      LwrtTLM(ng)=.TRUE.
#   endif
#  endif
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
        CALL ini_hmixcoef (ng, tile, iTLM)
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
          CALL ana_sponge (ng, tile, iTLM)
        END DO
      END IF
# endif
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration pass.
!=======================================================================

# if defined FOUR_DVAR              && \
    !defined I4DVAR_ANA_SENSITIVITY && \
    !(defined HESSIAN_SV            || \
      defined HESSIAN_SO            || \
      defined HESSIAN_FSV)
#  if defined OPT_OBSERVATIONS || defined TLM_CHECK || \
      defined WEAK_CONSTRAINT
!
!  Clear tangent linear state variables.
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL initialize_ocean (ng, tile, iTLM)
#   ifdef SOLVE3D
        CALL initialize_coupling (ng, tile, 0)
#   endif
#   ifdef SP4DVAR
        CALL initialize_ocean (ng, tile, iNLM)
        CALL initialize_forces (ng, tile, iNLM)
#   endif
      END DO

#  else
#   ifndef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  If first interation of the inner loop, clear all tangent linear
!  variables. In incrementatal 4DVAR, the tangent linear model is
!  started from rest on the first pass of the inner loop for each
!  outer loop iteration.
!-----------------------------------------------------------------------
!
      IF (inner.eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
#    ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, tile, iTLM)
#    endif
#    if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          CALL initialize_forces (ng, tile, iTLM)
#    endif
          CALL initialize_ocean (ng, tile, iTLM)
        END DO
!
!  Rewrite tangent linear initial NetCDF (record 1) with zero initial
!  conditions since the model needs to be started from at the first
!  pass of the inner loop.
!
        IF (Nrun.gt.1) THEN
          wrtRec=1
#    ifdef DISTRIBUTE
          CALL tl_wrt_ini (ng, MyRank, Tindex, wrtRec)
#    else
          CALL tl_wrt_ini (ng, -1, Tindex, wrtRec)
#    endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
#   endif
#  endif
# endif

# if defined SOLVE3D && !defined INI_FILE
!
!-----------------------------------------------------------------------
!  If analytical initial conditions, compute initial time-evolving
!  depths with zero free-surface.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_depth (ng, tile, iTLM)
      END DO
# endif

# if defined RBL4DVAR          || defined R4DVAR      || \
     defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
     defined TL_R4DVAR
!
!-----------------------------------------------------------------------
!  Initialize tangent linear bathymetry to zero.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL tl_bath (ng, tile)
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Set tangent linear model state variables initial conditions.
!-----------------------------------------------------------------------

# ifdef ANA_INITIAL
!
!  Analytical initial conditions for momentum and active tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_initial (ng, tile, iTLM)
        END DO
      END IF
# endif

# if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_passive (ng, tile, iTLM)
        END DO
      END IF
# endif

# if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_biology (ng, tile, iTLM)
        END DO
      END IF
# endif

# if defined ANA_SEDIMENT_NOT_YET && defined SOLVE3D
!
!  Analytical initial conditions for sediment.
!
      IF (nrrec(ng).eq.0) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_sediment (ng, tile, iTLM)
        END DO
      END IF
# endif

# ifdef I4DVAR_ANA_SENSITIVITY
!
!  Initialize with the weighted sum of all Lanczos vectors computed
!  from the first outer loop of the I4D-VAR Lanczos algorithm.
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ini_lanczos (ng, tile, Lnew(ng), Tindex)
      END DO

# else
!
!  Read in initial conditions for initial or restart NetCDF file.
!
#  ifdef INI_FILE
      CALL get_state (ng, iTLM, 1, ITL(ng), IniRec, Tindex)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      time(ng)=io_time                       ! needed for shared-memory
#  else
      IF (nrrec(ng).ne.0) THEN
        CALL get_state (ng, iTLM, 1, ITL(ng), IniRec, Tindex)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        time(ng)=io_time                     ! needed for shared-memory
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
#  ifdef DISTRIBUTE
        CALL get_wetdry (ng, MyRank, iTLM, IniRec(ng))
        CALL mp_bcasti (ng, iTLM, exit_flag)
#  else
        CALL get_wetdry (ng, -1, iTLM, IniRec(ng))
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
      CALL obs_initial (ng, iTLM, .FALSE.)
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
        CALL ana_perturb (ng, tile, iTLM)
      END DO
# endif

# ifdef TLM_CHECK
!
!-----------------------------------------------------------------------
!  Perturb tangent linear state variable according to the outer loop
!  iteration with the steepest descent direction of the gradient
!  (adjoint state).
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL tl_ini_perturb (ng, tile, Lnew(ng), Tindex)
      END DO
#  endif

# ifdef SOLVE3D
!!
!!----------------------------------------------------------------------
!!  Compute initial time-evolving depths.
!!----------------------------------------------------------------------
!!
!!    DO tile=first_tile(ng),last_tile(ng),+1
!!      CALL tl_set_depth (ng, tile, iTLM)
!!    END DO
!!
!!----------------------------------------------------------------------
!!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!!----------------------------------------------------------------------
!!
!!    DO tile=first_tile(ng),last_tile(ng),+1
!!      CALL tl_set_massflux (ng, tile, iTLM)
!!    END DO
!!
!!----------------------------------------------------------------------
!!  Compute initial S-coordinates vertical velocity. Compute initial
!!  density anomaly from potential temperature and salinity via equation
!!  of state for seawater.  Also compute other equation of state related
!!  quatities.
!!----------------------------------------------------------------------
!!
!!    DO tile=first_tile(ng),last_tile(ng),+1
!!      CALL tl_omega (ng, tile, iTLM)
!!      CALL tl_rho_eos (ng, tile, iTLM)
!!    END DO
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
          CALL ana_psource (ng, tile, iTLM)
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
      CALL close_inp (ng, iTLM)
      CALL check_multifile (ng, iTLM)
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
      CALL tl_get_idata (ng)
      CALL tl_get_data (ng)
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_masks (ng, tile, iTLM)
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
        CALL wpoints (ng, tile, iTLM)
      END DO
# endif

# ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  If available, read in first TLM impulse forcing and its application
!  time. In true weak constraint applications, the impulse records
!  after the initial are associated with the model error and are
!  processed with different statistics. If there is only one (initial)
!  impulse forcing available, the assimilation tis similar to strong
!  constraint but in observation space.
!-----------------------------------------------------------------------
!
      IF (nADJ(ng).lt.ntimes(ng)) THEN
        IniRec=1
        CALL get_state (ng, 7, 7, TLF(ng), IniRec, 1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif

# if defined ANA_DRAG && defined UV_DRAG_GRID
!
!-----------------------------------------------------------------------
!  Set analytical spatially varying bottom friction parameter.
!-----------------------------------------------------------------------
!
      IF (Nrun.eq.ERstr) THEN
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ana_drag (ng, tile, iTLM)
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
          CALL stiffness (ng, tile, iTLM)
        END DO
      END IF

# if defined FLOATS_NOT_YET || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
      CALL grid_coords (ng, iTLM)
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
!  Initialize time-stepping counter.
!-----------------------------------------------------------------------
!
      iic(ng)=ntstart(ng)
      CALL time_string (time(ng), time_code(ng))

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initiialization time wall clock.
!-----------------------------------------------------------------------
!
      DO thread=THREAD_RANGE
        CALL wclock_off (ng, iTLM, 2, __LINE__, MyFile)
      END DO
# endif
!
      RETURN
      END SUBROUTINE tl_initial
#else
      SUBROUTINE tl_initial
      RETURN
      END SUBROUTINE tl_initial
#endif
