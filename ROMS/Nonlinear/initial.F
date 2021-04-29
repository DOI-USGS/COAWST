#include "cppdefs.h"
      SUBROUTINE initial
!
!svn $Id: initial.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
      USE dateclock_mod,     ONLY : time_string
# ifdef DISTRIBUTE
      USE distribute_mod,    ONLY : mp_bcasti
# endif
#ifdef TLM_CHECK
      USE ini_adjust_mod,    ONLY : ini_perturb
#endif
      USE ini_hmixcoef_mod,  ONLY : ini_hmixcoef
      USE metrics_mod,       ONLY : metrics
#ifdef NESTING
      USE nesting_mod,       ONLY : nesting
#endif
#ifdef SOLVE3D
      USE set_depth_mod,     ONLY : set_depth0, set_depth
      USE omega_mod,         ONLY : omega
      USE rho_eos_mod,       ONLY : rho_eos
      USE set_massflux_mod,  ONLY : set_massflux
#endif
#ifdef MASKING
      USE set_masks_mod,     ONLY : set_masks
#endif
      USE stiffness_mod,     ONLY : stiffness
# ifdef COAWST_COUPLING
      USE mct_coupler_params
# endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2wav_coupling
      USE ocean_coupler_mod, ONLY : ocnfwav_coupling
#endif
      USE strings_mod,       ONLY : FoundError
#ifdef WET_DRY
      USE wetdry_mod,        ONLY : wetdry
#endif
#if defined PROPAGATOR || \
    (defined MASKING    && (defined READ_WATER || defined WRITE_WATER))
      USE wpoints_mod,       ONLY : wpoints
#endif
#ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2atm_coupling
      USE ocean_coupler_mod, ONLY : ocnfatm_coupling
#endif
#ifdef HYDRO_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2hyd_coupling
      USE ocean_coupler_mod, ONLY : ocnfhyd_coupling
#endif
#ifdef INWAVE_MODEL
      USE driver_inwave_mod, ONLY : inwave_init
#endif
!
      implicit none
!
!  Local variable declarations.
!
      logical :: update = .FALSE.

      integer :: Fcount, iw, ia, io, ih
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
#ifdef SP4DVAR
      real(dp) :: my_dstart
!
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
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
#elif defined I4DVAR      || defined RBL4DVAR          || \
      defined R4DVAR      || defined SENSITIVITY_4DVAR || \
      defined TL_RBL4DVAR || defined TL_R4DVAR
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
#ifdef SP4DVAR
        IF (INItime(ng).lt.0.0_dp) THEN
          my_dstart=dstart
        ELSE
          my_dstart=INItime(ng)/86400.0_dp
        END IF
#else
        tdays(ng)=dstart
#endif
        time(ng)=tdays(ng)*day2sec
!$OMP MASTER
        ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
        ntend(ng)=ntstart(ng)+ntimes(ng)-1
        ntfirst(ng)=ntstart(ng)
!$OMP END MASTER
!$OMP BARRIER
        step_counter(ng)=0
        CALL time_string (time(ng), time_code(ng))
      END DO
!
!  Initialize global diagnostics variables.
!
      avgke=0.0_dp
      avgpe=0.0_dp
      avgkp=0.0_dp
      volume=0.0_dp

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 2, __LINE__, MyFile)
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
!  Set initial conditions record to process. If applicable open existing
!  nonlinear model initial conditions NetCDF file and, if needed, define
!  new variables. Then, inquire about available variables.
!
      DO ng=1,Ngrids
        IF (LdefINI(ng)) THEN
          LdefINI(ng)=.FALSE.          ! needed to inquire variables IDs
# ifndef CORRELATION
          CALL def_ini (ng)
# endif
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IniRec(ng)=nrrec(ng)
          INI(ng)%Rindex=IniRec(ng)
        ELSE
          IniRec(ng)=INI(ng)%Rindex
        END IF
      END DO

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

# ifdef I4DVAR
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
!$OMP BARRIER
#endif
!
!-----------------------------------------------------------------------
!  Set application grid, metrics, and associated variables and
!  parameters.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (SetGridConfig(ng)) THEN
          CALL set_grid (ng, iNLM)
          SetGridConfig(ng)=.FALSE.
        END IF
      END DO
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
          CALL set_depth (ng, tile, iNLM)
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
#if defined INI_FILE && !defined RBL4DVAR_FCT_SENSITIVITY
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        time(ng)=io_time                     ! needed for shared-memory
      END DO
#else
!
!  If restart, read in initial conditions restart NetCDF file.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
          CALL get_state (ng, 0, 1, INI(ng)%name,                       &
     &                    IniRec(ng), Tindex(ng))
!$OMP END MASTER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
          END DO
!$OMP BARRIER
        END IF
      END DO
#endif

#if defined OBSERVATIONS && !defined RBL4DVAR_FCT_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations. Need to be done after processing initial conditions since
!  the correct initial time is needed to determine the first "ObsTime"
!  to process.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL obs_initial (ng, iNLM, .FALSE.)
!$OMP END MASTER
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!$OMP BARRIER
#endif

#if (defined ADJUST_BOUNDARY          || \
     defined ADJUST_STFLUX            || \
     defined ADJUST_WSTRESS)          && \
    (defined RBL4DVAR                 || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined TL_RBL4DVAR)
!
!-----------------------------------------------------------------------
!  Read in the surface forcing and or open boundary conditions
!  increments for RBL4D-Var from record IniRec of the NLM initial
!  NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        DO ng=1,Ngrids
!$OMP MASTER
          CALL get_state (ng, 5, 5, INI(ng)%name,                       &
     &                    IniRec(ng), Tindex(ng))
!$OMP END MASTER
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!  Compute time independent (Zt_avg1=0) anf initial time dependent
!  depths and level thicknesses.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_depth0 (ng, tile, iNLM)
          CALL set_depth  (ng, tile, iNLM)
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
          CALL set_massflux (ng, tile, iNLM)
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
          CALL omega (ng, tile, iNLM)
          CALL rho_eos (ng, tile, iNLM)
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
# if !defined CORRELATION
!
!  If applicable, close all input boundary, climatology, and forcing
!  NetCDF files and set associated parameters to the closed state. This
!  step is essential in iterative algorithms that run the full TLM
!  repetitively. Then, Initialize several parameters in their file
!  structure, so the appropriate input single or multi-file is selected
!  during initialization/restart.
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL close_inp (ng, iNLM)
        CALL check_multifile (ng, iNLM)
!$OMP END MASTER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
#  endif
!$OMP BARRIER
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!  If applicable, read in input data.
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL get_idata (ng)
        CALL get_data (ng)
!$OMP END MASTER
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
#  endif
!$OMP BARRIER
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
# endif
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

#if !defined CORRELATION
# ifdef NESTING
#  if defined MASKING || defined WET_DRY
!
!-----------------------------------------------------------------------
!  If nesting and Land/Sea masking, scale horizontal interpolation
!  weights to account for land contact points.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL nesting (ng, iNLM, nmask)
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  If nesting, process state fields initial conditions in the contact
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
#  ifndef SOLVE3D
            CALL nesting (ng, iNLM, n2dIC)        ! 2d momentum
#  endif
          END IF
        END DO
      END DO

#  ifdef SOLVE3D
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
#  endif
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

#if defined NLM_OUTER                || \
    defined RBL4DVAR                 || \
    defined RBL4DVAR_ANA_SENSITIVITY || \
    defined RBL4DVAR_FCT_SENSITIVITY || \
    defined TL_RBL4DVAR
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
#ifdef HYDRO_OCEAN
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled atm model.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO ih=1,Nhyd_grids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ocnfhyd_coupling (ng, ih, tile)
          END DO
!$OMP BARRIER
          IF (Master) WRITE (stdout,'(/)')
        END DO
      END DO
      DO ng=1,Ngrids
        DO ih=1,Nhyd_grids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ocn2hyd_coupling (ng, ih, tile)
          END DO
!$OMP BARRIER
          IF (Master) WRITE (stdout,'(/)')
        END DO
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
        INItime(ng)=time(ng)         ! Save NLM initial conditions time
        time(ng)=time(ng)-dt(ng)
      END DO

#ifdef INWAVE_MODEL
        DO ng=1,Ngrids
          CALL inwave_init (ng, MyRank, IniRec(ng))
        END DO
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initialization time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 2, __LINE__, MyFile)
        END DO
!$OMP BARRIER
      END DO
#endif
!
      RETURN
      END SUBROUTINE initial
