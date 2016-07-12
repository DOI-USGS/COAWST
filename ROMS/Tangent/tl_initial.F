#include "cppdefs.h"
#ifdef TANGENT
      SUBROUTINE tl_initial (ng)
!
!svn $Id: tl_initial.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      USE mod_fourdvar
# endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
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
# ifdef IS4DVAR_SENSITIVITY
      USE ini_lanczos_mod, ONLY : ini_lanczos
# endif
      USE metrics_mod, ONLY : metrics
# ifdef ADJUST_BOUNDARY
      USE mod_boundary, ONLY : initialize_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces, ONLY : initialize_forces
# endif
# if defined SENSITIVITY_4DVAR || \
     defined TL_W4DPSAS        || defined TL_W4DVAR  || \
     defined W4DPSAS           || defined W4DVAR
      USE tl_set_depth_mod, ONLY : tl_bath
# endif
# ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth
      USE tl_set_depth_mod, ONLY : tl_set_depth
      USE tl_omega_mod, ONLY : tl_omega
      USE tl_rho_eos_mod, ONLY : tl_rho_eos
      USE tl_set_massflux_mod, ONLY : tl_set_massflux
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
      integer :: Fcount, IniRec, Tindex, wrtRec
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
# elif defined IS4DVAR    || defined SENSITIVITY_4DVAR || \
       defined TL_W4DPSAS || defined TL_W4DVAR         || \
       defined W4DPSAS    || defined W4DVAR
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
# endif
        WRITE (stdout,20) 'TL_INITIAL: Configuring and ',               &
     &                    'initializing tangent linear model ...'
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
      IF (ANY(tl_VolCons(:,ng))) THEN
        tl_ubar_xs=0.0_r8
      END IF
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
        CALL wclock_on (ng, iTLM, 2)
      END DO
!$OMP BARRIER
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
!
!$OMP MASTER

#  ifndef IS4DVAR_SENSITIVITY
#   ifdef IS4DVAR
!
!  Set switch to create (TRUE) tangent linear initial conditions and
!  history NetCDF files or append (FALSE) to existing files. Then,
!  create tangent linear model initialization file and write zero
!  initial conditions for records 1 and 2.
!
      IF ((Nrun.eq.ERstr).and.(inner.eq.0)) THEN
        LdefITL(ng)=.TRUE.
        CALL tl_def_ini (ng)
        IF (exit_flag.ne.NoError) RETURN
        CALL tl_wrt_ini (ng, Tindex, 1)
        IF (exit_flag.ne.NoError) RETURN
        CALL tl_wrt_ini (ng, Tindex, 2)
        IF (exit_flag.ne.NoError) RETURN
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
!
!  Insure that forward and history file names are the same. In 4DVar,
!  the forward solution is computed by the nonlinear model and stored
!  on HISNAME NetCDF file.
!
#    ifdef TLM_CHECK
      FWD(ng)%name=TRIM(HIS(ng)%base)//'.nc'
      FWD(ng)%ncid=-1
#    else
      FWD(ng)%name=HIS(ng)%name
      FWD(ng)%ncid=HIS(ng)%ncid
#    endif
#   endif
#  endif
#  ifdef OBSERVATIONS
!
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations.
!
      CALL obs_initial (ng, iTLM, .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif
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
          CALL ana_grid (ng, tile, iTLM)
#  ifdef MASKING
          CALL ana_mask (ng, tile, iTLM)
#  endif
        END DO
!$OMP BARRIER
# else
!$OMP MASTER
        CALL get_grid (ng, iTLM)
!$OMP END MASTER
!$OMP BARRIER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iTLM, exit_flag)
#  endif
        IF (exit_flag.ne.NoError) RETURN
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
          CALL metrics (ng, tile, iTLM)
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
          CALL ana_wtype (ng, tile, iTLM)
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
            CALL ana_nudgcoef (ng, tile, iTLM)
          END DO
!$OMP BARRIER
        END IF
# else
        IF (Lnudging(ng)) THEN
!$OMP MASTER
          CALL get_nudgcoef (ng, iTLM)
!$OMP END MASTER
#  ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iTLM, exit_flag)
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
        CALL ini_hmixcoef (ng, tile, iTLM)
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
          CALL ana_sponge (ng, tile, iTLM)
        END DO
!$OMP BARRIER
      END IF
# endif
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration pass.
!=======================================================================

# if defined FOUR_DVAR  && !defined IS4DVAR_SENSITIVITY && \
   !(defined HESSIAN_SV ||  defined HESSIAN_SO || defined HESSIAN_FSV)
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
      END DO
!$OMP BARRIER

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
!$OMP BARRIER
!
!  Rewrite tangent linear initial NetCDF (record 1) with zero initial
!  conditions since the model needs to be started from at the first
!  pass of the inner loop.
!
        IF (Nrun.gt.1) THEN
          wrtRec=1
!$OMP MASTER
          CALL tl_wrt_ini (ng, Tindex, wrtRec)
!$OMP END MASTER
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
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
        CALL set_depth (ng, tile)
      END DO
!$OMP BARRIER
# endif

# if defined SENSITIVITY_4DVAR || \
     defined TL_W4DPSAS        || defined TL_W4DVAR || \
     defined W4DPSAS           || defined W4DVAR
!
!-----------------------------------------------------------------------
!  Initialize tangent linear bathymetry to zero.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL tl_bath (ng, tile)
      END DO
!$OMP BARRIER
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
!$OMP BARRIER
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
!$OMP BARRIER
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
!$OMP BARRIER
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
!$OMP BARRIER
      END IF
# endif

# ifdef IS4DVAR_SENSITIVITY
!
!  Initialize with the weighted sum of all Lanczos vectors computed
!  from the first outer loop of the IS4DVAR Lanczos algorithm.
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ini_lanczos (ng, tile, Lnew(ng), Tindex)
      END DO
!$OMP BARRIER

# else
!
!  Read in initial conditions for initial or restart NetCDF file.
!
#  ifdef INI_FILE
!$OMP MASTER
      CALL get_state (ng, iTLM, 1, ITL(ng)%name, IniRec, Tindex)
!$OMP END MASTER
!$OMP BARRIER
      IF (exit_flag.ne.NoError) RETURN
      time(ng)=io_time                       ! needed for shared-memory
#  else
      IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
        CALL get_state (ng, iTLM, 1, ITL(ng)%name, IniRec, Tindex)
!$OMP END MASTER
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
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
!$OMP MASTER
        CALL get_wetdry (ng, iTLM, IniRec(ng))
!$OMP END MASTER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iTLM, exit_flag)
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
        CALL ana_perturb (ng, tile, iTLM)
      END DO
!$OMP BARRIER
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
!$OMP BARRIER
#  endif

# ifdef SOLVE3D
!!
!!----------------------------------------------------------------------
!!  Compute initial time-evolving depths.
!!----------------------------------------------------------------------
!!
!!    DO tile=first_tile(ng),last_tile(ng),+1
!!      CALL tl_set_depth (ng, tile)
!!    END DO
!!$OMP BARRIER
!!
!!----------------------------------------------------------------------
!!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!!----------------------------------------------------------------------
!!
!!    DO tile=first_tile(ng),last_tile(ng),+1
!!      CALL tl_set_massflux (ng, tile)
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
!!      CALL tl_omega (ng, tile)
!!      CALL tl_rho_eos (ng, tile)
!!    END DO
!!$OMP BARRIER
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
      CALL check_multifile (ng, iTLM)
!$OMP END MASTER
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag)
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
      CALL close_inp (ng, iTLM)
      CALL tl_get_idata (ng)
      CALL tl_get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag)
# endif
      IF (exit_flag.ne.NoError) RETURN

# ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_masks (ng, tile, iTLM)
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
        CALL wpoints (ng, tile, iTLM)
      END DO
!$OMP BARRIER
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
!$OMP MASTER
        CALL get_state (ng, 7, 7, TLF(ng)%name, IniRec, 1)
!$OMP END MASTER
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
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
          CALL stiffness (ng, tile, iTLM)
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
      CALL grid_coords (ng, iTLM)
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
        CALL wclock_off (ng, iTLM, 2)
      END DO
!$OMP BARRIER
# endif
      RETURN
      END SUBROUTINE tl_initial
#else
      SUBROUTINE tl_initial
      RETURN
      END SUBROUTINE tl_initial
#endif
