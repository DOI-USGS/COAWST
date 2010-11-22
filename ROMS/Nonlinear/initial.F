#include "cppdefs.h"
      SUBROUTINE initial (ng)
!
!svn $Id: initial.F 790 2008-10-14 19:32:07Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
#ifdef TLM_CHECK
      USE ini_adjust_mod, ONLY : ini_perturb
#endif
      USE ini_hmixcoef_mod, ONLY : ini_hmixcoef
      USE metrics_mod, ONLY : metrics
#ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth
      USE omega_mod, ONLY : omega
      USE rho_eos_mod, ONLY : rho_eos
      USE set_massflux_mod, ONLY : set_massflux
#endif
      USE stiffness_mod, ONLY : stiffness
#if defined PROPAGATOR || \
   (defined MASKING && (defined READ_WATER || defined WRITE_WATER ))
      USE wpoints_mod, ONLY : wpoints
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2wav_coupling
#endif
#ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : ocn2atm_coupling
#endif
#ifdef REFINED_GRID
      USE init_child_hindices_mod, ONLY : init_child_hindices
      USE init_parent_hindices_mod, ONLY : init_parent_hindices
      USE get_2dparent_data_mod, ONLY : get_2dparent_data
# ifdef SOLVE3D
      USE get_3dparent_data_mod, ONLY : get_3dparent_data
# endif
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, save :: First = .TRUE.
      logical :: update = .FALSE.

      integer :: LBi, UBi, LBj, UBj
      integer :: IniRec, Tindex, subs, tile, thread
#if defined ADJUST_BOUNDARY || \
    defined ADJUST_STFLUX   || defined ADJUST_WSTRESS
      integer :: irec
#endif
      integer :: my_numthreads
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
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
#ifdef FLOATS
      nf(ng)=0
      nfp1(ng)=1
      nfm1(ng)=4
      nfm2(ng)=3
      nfm3(ng)=2
#endif
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

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        CALL wclock_on (ng, iNLM, 1)
      END DO
!$OMP END PARALLEL DO
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
      IF (First) THEN
        First=.FALSE.
# ifdef ANA_INITIAL
        LdefINI(ng)=.TRUE.
# endif
        LdefHIS(ng)=.TRUE.
        CALL def_ini (ng)
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
        IniRec=nrrec(ng)
        tINIindx(ng)=IniRec
      ELSE
        IniRec=tINIindx(ng)
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Initialize open boundary counter for storage arrays.
!
      OBCcount(ng)=0
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Initialize surface forcing counter for storage arrays.
!
      SFcount(ng)=0
# endif
!
!  Reset nonlinear history time record counters. These counters are
!  reset on every iteration pass. This file is created on the first
!  iteration pass.
!
      tHISindx(ng)=0
      NrecHIS(ng)=0

# ifdef IS4DVAR
!
!  Activate switches to writting data into average, history and
!  restart files.
!
      LwrtAVG(ng)=.TRUE.
      LwrtHIS(ng)=.TRUE.
      LwrtRST(ng)=.TRUE.
# endif
#endif
#ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations.
!-----------------------------------------------------------------------
!
      CALL obs_initial (ng, iNLM, .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_grid (ng, TILE, iNLM)
# ifdef MASKING
            CALL ana_mask (ng, TILE, iNLM)
# endif
          END DO
        END DO
!$OMP END PARALLEL DO
#else
        CALL get_grid (ng, iNLM)
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        if (exit_flag.ne.NoError) RETURN
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set vertical S-coordinate transformation function.
!-----------------------------------------------------------------------
!
        CALL set_scoord (ng)
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set barotropic time-steps average weighting function.
!-----------------------------------------------------------------------
!
        CALL set_weights (ng)
#endif
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL metrics (ng, TILE, iNLM)
#if defined PROPAGATOR || \
   (defined MASKING && (defined READ_WATER || defined WRITE_WATER ))
            CALL wpoints (ng, TILE, iNLM)
#endif
          END DO
        END DO
!$OMP END PARALLEL DO

#ifdef NUDGING_COFF
!
!-----------------------------------------------------------------------
!  If appropriate, set nudging coefficiests time scales.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_nudgcoef (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL ini_hmixcoef (ng, TILE, iNLM)
        END DO
      END DO
!$OMP END PARALLEL DO

#if defined VISC_GRID || defined DIFF_GRID || defined SPONGE
!
!-----------------------------------------------------------------------
!  Rescale horizontal mixing coefficients according to the local
!  grid size. If applicable, increases horizontal mixing in sponge
!  areas.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_hmixcoef (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
#endif
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration run.
!=======================================================================

#ifdef TLM_CHECK
!
!  Clear state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL initialize_ocean (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
#endif

#if defined SOLVE3D && !defined INI_FILE
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
#endif
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.
!-----------------------------------------------------------------------

#ifdef ANA_INITIAL
!
!  Analytical initial conditions for momentum and active tracers.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_initial (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif

#if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_passive (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif

#if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology tracers.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_biology (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif

#if defined ANA_SEDIMENT && defined SOLVE3D
!
!  Analytical initial conditions for sediment tracers.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_sediment (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#endif

#ifdef INI_FILE
!
!  Read in initial conditions from initial NetCDF file.
!
      CALL get_state (ng, iNLM, 1, INIname(ng), IniRec, Tindex)
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
      IF (exit_flag.ne.NoError) RETURN
#else
!
!  If restart, read in initial conditions restart NetCDF file.
!
      IF (nrrec(ng).ne.0) THEN
        CALL get_state (ng, 0, 1, INIname(ng), IniRec, Tindex)
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif

#if (defined W4DPSAS || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY)          && \
    (defined ADJUST_BOUNDARY || defined ADJUST_WSTRESS ||\
     defined ADJUST_STFLUX)
!
!-----------------------------------------------------------------------
!  Read in the surface forcing and/or obc increments for PSAS from
!  record 1 of the TLM initial NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Nrun.gt.1) THEN
        CALL get_state (ng, 12, 1, ITLname(ng), 1, Tindex)
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile,Tindex)                      &
!$OMP             SHARED(ng,numthreads,Lnew)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ini_perturb (ng, TILE, Lnew(ng), Tindex)
          END DO
        END DO
      END IF
!$OMP END PARALLEL DO
#endif

#ifdef SOLVE3D
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
          CALL set_massflux (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Compute initial S-coordinates vertical velocity. Compute initial
!  density anomaly from potential temperature and salinity via equation
!  of state for seawater.  Also compute other equation of state related
!  quatities.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL omega (ng, TILE)
          CALL rho_eos (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
#endif

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
      IF (Nrun.eq.ERstr) THEN
        OBC_time(1,ng)=time(ng)
        DO irec=2,Nbrec(ng)
          OBC_time(irec,ng)=OBC_time(irec-1,ng)+nOBC(ng)*dt(ng)
        END DO
      END IF
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  If first pass of iteration loop, set time of surface forcing
!  adjustment.
!
      IF (Nrun.eq.ERstr) THEN
        SF_time(1,ng)=time(ng)
        DO irec=2,Nfrec(ng)
          SF_time(irec,ng)=SF_time(irec-1,ng)+nSFF(ng)*dt(ng)
        END DO
      END IF
# endif
!
!  If applicable, read in input data.
!
# ifdef TIMELESS_DATA
      CALL get_idata (ng)
# endif
      CALL get_data (ng)
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
      IF (exit_flag.ne.NoError) RETURN
#endif

#if defined NLM_OUTER || defined TL_W4DPSAS          || \
    defined W4DPSAS   || defined W4DPSAS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Read in convolved adjoint impulse forcing (first record) and its
!  application time.
!-----------------------------------------------------------------------
!
      IF (SporadicImpulse) THEN
        FrcRec(ng)=1
        CALL get_state (ng, 7, 7, TLFname(ng), FrcRec(ng), 1)
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iTLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif
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
            CALL stiffness (ng, TILE, iNLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF

#if defined FLOATS || defined STATIONS
!
!-----------------------------------------------------------------------
!  If applicable, convert initial locations to fractional grid
!  coordinates.
!-----------------------------------------------------------------------
!
      CALL grid_coords (ng, iNLM)
#endif

#ifdef WAVES_OCEAN
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled wave model.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL ocn2wav_coupling (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
      IF (Master) WRITE (stdout,'(/)')
#endif

#ifdef AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Read in initial forcing from coupled wave model.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL ocn2atm_coupling (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
      IF (Master) WRITE (stdout,'(/)')
#endif
#ifdef REFINED_GRID
!-----------------------------------------------------------------------
!  Compute indices of children grid locations in the parent grid.
!  For 2-way, this needs to be done for all grids, not just ng>1.
!-----------------------------------------------------------------------
      IF (ng.gt.1) THEN
        CALL init_parent_hindices (ng, TILE)
      END IF
      IF (ng.lt.NestedGrids) THEN
        CALL init_child_hindices (ng, TILE)
      ENd IF
!-----------------------------------------------------------------------
!  Obtain initial boudnary conditions from the parent data.
!-----------------------------------------------------------------------
      IF (ng.gt.1) THEN
        CALL get_2dparent_data (ng, TILE)
# ifdef SOLVE3D
        CALL get_3dparent_data (ng, TILE)
# endif
      END IF
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off initiialization time wall clock.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        CALL wclock_off (ng, iNLM, 1)
      END DO
!$OMP END PARALLEL DO
#endif
      RETURN
      END SUBROUTINE initial
