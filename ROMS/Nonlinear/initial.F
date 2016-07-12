#include "cppdefs.h"
      SUBROUTINE initial
!
!svn $Id: initial.F 790 2008-10-14 19:32:07Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes all model variables.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BBL_MODEL
      USE mod_bbl
#endif
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_grid
#ifdef CICE_MODEL
      USE CICE_InitMod
#endif
      USE mod_iounits
      USE mod_ncparam
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE analytical_mod
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
#ifdef TLM_CHECK
      USE ini_adjust_mod, ONLY : ini_perturb
#endif
      USE ini_hmixcoef_mod, ONLY : ini_hmixcoef
      USE metrics_mod, ONLY : metrics
#ifdef NESTING
      USE nesting_mod, ONLY : nesting
#endif
#ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth0, set_depth
      USE omega_mod, ONLY : omega
      USE rho_eos_mod, ONLY : rho_eos
      USE set_massflux_mod, ONLY : set_massflux
#endif
#ifdef ICE_STRENGTH_QUAD
      USE ini_strengthcoef_mod, ONLY : ini_strengthcoef
#endif
#ifdef MASKING
      USE set_masks_mod, ONLY : set_masks
#endif
      USE stiffness_mod, ONLY : stiffness
# ifdef COAWST_COUPLING
      USE mct_coupler_params
# endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2wav_coupling
      USE ocean_coupler_mod, ONLY : ocnfwav_coupling
#endif
#ifdef WET_DRY
      USE wetdry_mod, ONLY : wetdry
#endif
#if defined PROPAGATOR || \
    (defined MASKING    && (defined READ_WATER || defined WRITE_WATER))
      USE wpoints_mod, ONLY : wpoints
#endif
#ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2atm_coupling
      USE ocean_coupler_mod, ONLY : ocnfatm_coupling
#endif
#ifdef INWAVE_MODEL
      USE driver_inwave_mod, ONLY : inwave_init
#endif
!
      implicit none
!
!  Local variable declarations.
!
      logical, save :: First = .TRUE.
      logical :: update = .FALSE.

      integer :: Fcount, iw, ia, io
      integer :: ng, thread, tile
#ifdef NESTING
      integer :: ig, nl
      integer :: cr, i, m
#endif

      integer, dimension(Ngrids) :: IniRec, Tindex

#if defined ADJUST_BOUNDARY || \
    defined ADJUST_STFLUX   || defined ADJUST_WSTRESS
      integer :: irec
#endif
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
!$OMP MASTER
      IF (Master) THEN
#if defined PERTURBATION
        WRITE (stdout,10) Nrun
 10     FORMAT (/,' <<<< Ensemble/Perturbation Run: ',i5.5,' >>>>',/)
#elif defined IS4DVAR    || defined SENSITIVITY_4DVAR || \
      defined TL_W4DPSAS || defined TL_W4DVAR         || \
      defined W4DPSAS    || defined W4DVAR
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
#endif
        WRITE (stdout,20) 'INITIAL: Configuring and initializing ',     &
     &                    'forward nonlinear model ...'
 20     FORMAT (/,1x,a,a,/,1x,'*******')
      END IF
!$OMP END MASTER
!
!-----------------------------------------------------------------------
!  Initialize time stepping indices and counters.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
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
#ifdef FLOATS
        nf(ng)=0
        nfp1(ng)=1
        nfm1(ng)=4
        nfm2(ng)=3
        nfm3(ng)=2
#endif
!
        IniRec(ng)=nrrec(ng)
        Tindex(ng)=1
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
        step_counter(ng)=0
        CALL time_string (time(ng), time_code(ng))
#ifdef CICE_MODEL
! kca - calculate timesteps for restart output
        tspy(ng) = 60*60*24*365/dt(ng)
        tspd(ng) = 60*60*24/dt(ng)
#endif
      END DO

#ifdef CICE_MODEL
! kca - calculate timesteps for restart output
      dleftinSep = 273-dstart
      Oct = dleftinSep
      Nov = Oct + 31
      Dec = Nov + 30
      Jan = Dec + 31
      Feb = Jan + 31
      Mar = Feb + 28
      Apr = Mar + 31
      May = Apr + 30
      Jun = May + 31
      Jul = Jun + 30
      Aug = Jul + 31
      Sep = Aug + 31
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 2)
        END DO
      END DO
!$OMP BARRIER
#endif

#ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  If variational data assimilation, reset several IO switches and
!  variables.
!-----------------------------------------------------------------------
!
!  Set switch to create (true) nonlinear model initial conditions and
!  histroy NetCDF files or append (false) to an existing file.  Set
!  record to read from initial NetCDF file.
!
!$OMP MASTER
      IF (First) THEN
        First=.FALSE.
        DO ng=1,Ngrids
# ifdef ANA_INITIAL
          LdefINI(ng)=.TRUE.
# endif
          LdefHIS(ng)=.TRUE.
          CALL def_ini (ng)
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
          IF (exit_flag.ne.NoError) RETURN
          IniRec(ng)=nrrec(ng)
          INI(ng)%Rindex=IniRec(ng)
        END DO
      ELSE
        DO ng=1,Ngrids
          IniRec=INI(ng)%Rindex
        END DO
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Initialize open boundary counter for storage arrays.
!
      DO ng=1,Ngrids
        OBCcount(ng)=0
      END DO
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Initialize surface forcing counter for storage arrays.
!
      DO ng=1,Ngrids
        SFcount(ng)=0
      END DO
# endif
!
!  Reset nonlinear history time record counters. These counters are
!  reset on every iteration pass. This file is created on the first
!  iteration pass.
!
      DO ng=1,Ngrids
        HIS(ng)%Rindex=0
        Fcount=HIS(ng)%Fcount
        HIS(ng)%Nrec(Fcount)=0
      END DO

# ifdef IS4DVAR
!
!  Activate switches to writting data into average, history and
!  restart files.
!
      DO ng=1,Ngrids
        LwrtAVG(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        LwrtRST(ng)=.TRUE.
      END DO
# endif
!$OMP END MASTER
!$OMP BARRIER
#endif
#ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL obs_initial (ng, iNLM, .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!$OMP END MASTER
      END DO
!$OMP BARRIER
#endif
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
#ifdef ANA_GRID
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_grid (ng, tile, iNLM)
# ifdef MASKING
            CALL ana_mask (ng, tile, iNLM)
# endif
          END DO
        END DO
!$OMP BARRIER
#else
        DO ng=1,Ngrids
!$OMP MASTER
          CALL get_grid (ng, iNLM)
!$OMP END MASTER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
        END DO
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set vertical S-coordinate transformation function.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP MASTER
          CALL set_scoord (ng)
!$OMP END MASTER
        END DO
!$OMP BARRIER
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP MASTER
          CALL set_weights (ng)
!$OMP END MASTER
        END DO
!$OMP BARRIER
#endif
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL metrics (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END DO
#ifdef ICE_STRENGTH_QUAD
!
!-----------------------------------------------------------------------
!  Compute ice strength as a function of grid size.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ini_strengthcoef (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END DO
#endif
#ifdef NESTING
!
!-----------------------------------------------------------------------
!  If nesting, initialize grid spacing (on_u and om_v) in REFINED(:)
!  structure.  They are used to impose mass flux at the finer grid
!  physical boundaries.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL nesting (ng, iNLM, ndxdy)
        END DO
#endif

#if defined WTYPE_GRID && defined ANA_WTYPE     && \
   (defined LMD_SKPP   || defined SOLAR_SOURCE)
!
!-----------------------------------------------------------------------
!  Set spatially varying Jerlov water type.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_wtype (ng, tile, iNLM)
          END DO
        END DO
!$OMP BARRIER
#endif
!
!-----------------------------------------------------------------------
!  If appropriate, set spatially varying nudging coefficients time
!  scales.
!-----------------------------------------------------------------------
!
#ifdef ANA_NUDGCOEF
        DO ng=1,Ngrids
          IF (Lnudging(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ana_nudgcoef (ng, tile, iNLM)
            END DO
!$OMP BARRIER
          END IF
        END DO
#else
        DO ng=1,Ngrids
          IF (Lnudging(ng)) THEN
!$OMP MASTER
            CALL get_nudgcoef (ng, iNLM)
!$OMP END MASTER
# ifdef DISTRIBUTE
            CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients. If applicable, scale
!  mixing coefficients according to the grid size (smallest area).
#ifndef ANA_SPONGE
!  Also increase their values in sponge areas using the "visc_factor"
!  and/or "diff_factor" read from input Grid NetCDF file.
#endif
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ini_hmixcoef (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO

#ifndef OFFLINE_FLOATS
#ifdef ANA_SPONGE
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing coefficients in sponge areas using
!  analytical functions.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (Lsponge(ng)) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_sponge (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif
#endif /* OFFLINE_FLOATS */
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration run.
!=======================================================================

#ifdef TLM_CHECK
!
!  Clear state variables.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
#endif

#if defined SOLVE3D && !defined INI_FILE
!
!-----------------------------------------------------------------------
!  If analytical initial conditions, compute initial time-evolving
!  depths with zero free-surface.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_depth (ng, tile)
        END DO
!$OMP BARRIER
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.
!-----------------------------------------------------------------------

#ifdef ANA_INITIAL
!
!  Analytical initial conditions for momentum and active tracers.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_initial (ng, tile, iNLM)
          END DO
          io_time = time(ng)
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_passive (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology tracers.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_biology (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined ANA_SEDIMENT && defined SOLVE3D
!
!  Analytical initial conditions for sediment tracers.
!

      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_sediment (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined ANA_VEGETATION && defined SOLVE3D
!
!  Analytical initial conditions for vegetation types.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_vegetation (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#ifndef OFFLINE_FLOATS
#ifdef INI_FILE
!
!  Read in initial conditions from initial NetCDF file.
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL get_state (ng, iNLM, 1, INI(ng)%name,                      &
     &                  IniRec(ng), Tindex(ng))
!$OMP END MASTER
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
        time(ng)=io_time                     ! needed for shared-memory
      END DO
#else
!
!  If restart, read in initial conditions restart NetCDF file.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
          CALL get_state (ng, 0, 1, INI(ng)%name,                      &
     &                    IniRec(ng), Tindex(ng))
!$OMP END MASTER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
          time(ng)=io_time                   ! needed for shared-memory
        END IF
      END DO
#endif

#ifdef ANA_ICE
!
!  Analytical initial conditions for sea ice.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_ice (ng, TILE, iNLM)
          END DO
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
          time(ng)=io_time                   ! needed for shared-memory
        END IF
      END DO
#endif
#ifdef WET_DRY
!
!-----------------------------------------------------------------------
!  Process initial wet/dry masks.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!
!  If restart, read in wet/dry masks.
!
        IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
          CALL get_wetdry (ng, iNLM, IniRec(ng))
!$OMP END MASTER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
        ELSE
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if (defined W4DPSAS || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY)          && \
    (defined ADJUST_BOUNDARY || defined ADJUST_WSTRESS ||\
     defined ADJUST_STFLUX)
!
!-----------------------------------------------------------------------
!  Read in the surface forcing and/or open boundary conditions
!  increments for PSAS from record IniRec of the NLM initial NetCDF
!  file.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        DO ng=1,Ngrids
!$OMP MASTER
          CALL get_state (ng, 12, 12, INI(ng)%name,                     &
     &                    IniRec(ng), Tindex(ng))
!$OMP END MASTER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END IF
#endif

#ifdef TLM_CHECK
!
!-----------------------------------------------------------------------
!  Add a perturbation to nonlinear state variable according to the outer
!  loop iteration with the steepest descent direction of the gradient
!  (adjoint state).
!-----------------------------------------------------------------------
!
      IF (outer.ge.1) THEN
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ini_perturb (ng, tile, Lnew(ng), Tindex)
          END DO
!$OMP BARRIER
        END DO
      END IF
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute initial unperturbed depths. Notice that during the
!  initialization "Zt_avg1" is always zero.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_depth0 (ng, tile)
          CALL set_depth  (ng, tile)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_massflux (ng, tile)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial S-coordinates vertical velocity. Compute initial
!  density anomaly from potential temperature and salinity via equation
!  of state for seawater.  Also compute other equation of state related
!  quatities.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL omega (ng, tile)
          CALL rho_eos (ng, tile)
        END DO
!$OMP BARRIER
      END DO
#endif

#ifdef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks position, direction, special flag, and mass
!  transport nondimensional shape profile with analytcal expressions.
!  Point sources are at U- and V-points. We need to get their positions
!  to process internal Land/Sea masking arrays during initialization.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_psource (ng, tile, iNLM)
          END DO
        END IF
!$OMP BARRIER
      END DO
#endif
#endif /* !OFFLINE_FLOATS */

#if defined FOUR_DVAR || !defined TANGENT || !defined ADJOINT
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------

# ifdef ADJUST_BOUNDARY
!
!  If first pass of iteration loop, set time of open boundary
!  adjustment.
!
!$OMP MASTER
      IF (Nrun.eq.ERstr) THEN
        DO ng=1,Ngrids
          OBC_time(1,ng)=time(ng)
          DO irec=2,Nbrec(ng)
            OBC_time(irec,ng)=OBC_time(irec-1,ng)+nOBC(ng)*dt(ng)
          END DO
        END DO
      END IF
!$OMP END MASTER
!$OMP BARRIER
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  If first pass of iteration loop, set time of surface forcing
!  adjustment.
!
!$OMP MASTER
      IF (Nrun.eq.ERstr) THEN
        DO ng=1,Ngrids
          SF_time(1,ng)=time(ng)
          DO irec=2,Nfrec(ng)
            SF_time(irec,ng)=SF_time(irec-1,ng)+nSFF(ng)*dt(ng)
          END DO
        END DO
      END IF
!$OMP END MASTER
!$OMP BARRIER
# endif
!
!  If applicable, check input multi-files where the time records for
!  a particular field are split into several NetCDF files.  Initialize
!  several parameters in the file structure so the appropriate input
!  file is selected during initialization/restart.
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL check_multifile (ng, iNLM)
!$OMP END MASTER
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  If applicable, read in input data.
!
      DO ng=1,Ngrids
        CALL close_inp (ng, iNLM)
!$OMP MASTER
        CALL get_idata (ng)
        CALL get_data (ng)
!$OMP END MASTER
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif

#ifdef MASKING
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_masks (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
#endif

#ifdef NESTING
# if defined MASKING || defined WET_DRY
!
!-----------------------------------------------------------------------
!  If nesting and Land/Sea masking, scale horizontal interpolation
!  weights to account for land contact points.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL nesting (ng, iNLM, nmask)
        END DO
!$OMP BARRIER
      END DO
# endif
!
!-----------------------------------------------------------------------
!  If nesting, process state fields initial conditons in the contact
!  regions.
!-----------------------------------------------------------------------
!
!  Free-surface and 2D-momentum.
!
      DO nl=1,NestLayers
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          IF (ANY(CompositeGrid(:,ng))) THEN
            CALL nesting (ng, iNLM, nFSIC)        ! free-surface
# ifndef SOLVE3D
            CALL nesting (ng, iNLM, n2dIC)        ! 2d momentum
# endif
          END IF
        END DO
      END DO

# ifdef SOLVE3D
!
!  Determine vertical indices and vertical interpolation weights in
!  the contact zone using initial unperturbed depth arrays.
!
      DO ng=1,Ngrids
        CALL nesting (ng, iNLM, nzwgt)
      END DO
!
!  3D-momentum and tracers.
!
      DO nl=1,NestLayers
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          IF (ANY(CompositeGrid(:,ng))) THEN
            CALL nesting (ng, iNLM, n3dIC)        ! 3D momentum
            CALL nesting (ng, iNLM, nTVIC)        ! Tracer variables
          END IF
        END DO
      END DO
# endif
#endif

#if defined PROPAGATOR || \
   (defined MASKING    && (defined READ_WATER || defined WRITE_WATER ))
!
!-----------------------------------------------------------------------
!  Set variables associated with the processing water points and/or
!  size of packed state arrays.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL wpoints (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
#endif

#if defined NLM_OUTER || defined TL_W4DPSAS          || \
    defined W4DPSAS   || defined W4DPSAS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Read in convolved adjoint impulse forcing (first record) and its
!  application time.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (SporadicImpulse(ng)) THEN
          FrcRec(ng)=1
!$OMP MASTER
          CALL get_state (ng, 7, 7, TLF(ng)%name, FrcRec(ng), 1)
!$OMP END MASTER
!$OMP BARRIER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iTLM, exit_flag)
# endif
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
#endif

#if defined ANA_DRAG && defined UV_DRAG_GRID
!
!-----------------------------------------------------------------------
!  Set analytical spatially varying bottom friction parameter.
!-----------------------------------------------------------------------
!
      IF (Nrun.eq.ERstr) THEN
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_drag (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      IF (Lstiffness) THEN
        Lstiffness=.FALSE.
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL stiffness (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END DO
      END IF

#if defined FLOATS || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL grid_coords (ng, iNLM)
!$OMP END MASTER
!$OMP BARRIER
      END DO
#endif

#ifdef WAVES_OCEAN
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled wave model.
!-----------------------------------------------------------------------
!
      allocate(roms_fwcoup(Nocn_grids))
      allocate(roms_2wcoup(Nocn_grids))
      DO io=1,Nocn_grids
        roms_fwcoup(io)=0
        roms_2wcoup(io)=0
      END DO
      DO iw=1,Nwav_grids
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ocnfwav_coupling (ng, iw, tile)
          END DO
!$OMP BARRIER
          IF (Master) WRITE (stdout,'(/)')
        END DO
      END DO
      DO iw=1,Nwav_grids
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ocn2wav_coupling (ng, iw, tile)
          END DO
!$OMP BARRIER
          IF (Master) WRITE (stdout,'(/)')
        END DO
      END DO
#endif
#ifdef AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled atm model.
!-----------------------------------------------------------------------
!
      allocate(roms_facoup(Nocn_grids))
      allocate(roms_2acoup(Nocn_grids))
      DO io=1,Nocn_grids
        roms_facoup(io)=0
        roms_2acoup(io)=0
      END DO
      DO ng=1,Ngrids
        DO ia=1,Natm_grids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ocnfatm_coupling (ng, ia, tile)
          END DO
!$OMP BARRIER
          IF (Master) WRITE (stdout,'(/)')
        END DO
      END DO
      DO ng=1,Ngrids
        DO ia=1,Natm_grids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ocn2atm_coupling (ng, ia, tile)
          END DO
!$OMP BARRIER
          IF (Master) WRITE (stdout,'(/)')
        END DO
      END DO
#endif
#ifdef INWAVE_MODEL
        DO ng=1,Ngrids
          CALL inwave_init (ng, MyRank)
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter and clock.
!-----------------------------------------------------------------------
!
!  Subtract one time unit to avoid special case due to initialization
!  in the main time-stepping routine.
!
      DO ng=1,Ngrids
        iic(ng)=ntstart(ng)-1
        time(ng)=time(ng)-dt(ng)
      END DO

#ifdef CICE_MODEL
!
!-----------------------------------------------------------------------
! Initialize the ice model
!-----------------------------------------------------------------------
!
      CALL cice_init
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initialization time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 2)
        END DO
!$OMP BARRIER
      END DO
#endif

      RETURN
      END SUBROUTINE initial
