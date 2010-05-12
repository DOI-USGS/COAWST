#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_initial (ng)
!
!svn $Id: ad_initial.F 431 2009-12-26 20:36:20Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
# ifdef ANALYTICAL
      USE analytical_mod
# endif
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
# if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
     defined SO_SEMI
      USE adsen_initial_mod, ONLY : adsen_initial
# endif
      USE ini_hmixcoef_mod, ONLY : ini_hmixcoef
      USE metrics_mod, ONLY : metrics
# ifdef SOLVE3D
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
      integer :: IniRec, Tindex, subs, tile, thread

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
      synchro_flag(ng)=.TRUE.
      first_time(ng)=0
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

      tdays(ng)=dstart+dt(ng)*FLOAT(ntimes(ng)-ntfirst(ng)+1)*sec2day
      time(ng)=tdays(ng)*day2sec
      ntstart(ng)=ntimes(ng)+1
      ntend(ng)=ntfirst(ng)
      ntfirst(ng)=ntend(ng)
      CALL time_string (time(ng), time_code(ng))

      IniRec=nrrec(ng)
      Tindex=1

      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
! Initialize few parameters.
!
      ad_ubar_xs=0.0_r8

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Start time wall clocks.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        CALL wclock_on (ng, iADM, 1)
      END DO
!$OMP END PARALLEL DO
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
      IF (Nrun.eq.ERstr) THEN
        LdefADJ(ng)=.TRUE.
      END IF
!
!  Activate switch to write adjoint NetCDF file.
!
      LwrtADJ(ng)=.TRUE.

#  ifndef WEAK_CONSTRAINT
!
!  Insure that forward and history file names are the same. In 4DVar,
!  the forward solution is computed by the nonlinear model and stored
!  HISNAME NetCDF file.
!
#   ifdef TLM_CHECK
      FWDname(ng)=HISbase(ng)
      ncFWDid(ng)=-1
#   else
      FWDname(ng)=HISname(ng)
      ncFWDid(ng)=ncHISid(ng)
#   endif
#  endif
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
#  ifdef OBSERVATIONS
!
!  Initialize various variables needed for processing observations
!  backwards in time.
!
      CALL obs_initial (ng, iADM, .TRUE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_grid (ng, TILE, iADM)
#  ifdef MASKING
            CALL ana_mask (ng, TILE, iADM)
#  endif
#  if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS || defined SO_SEMI             || \
      defined SENSITIVITY_4DVAR
            CALL ana_scope (ng, TILE, iADM)
#  endif
          END DO
        END DO
!$OMP END PARALLEL DO
# else
        CALL get_grid (ng, iADM)
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, exit_flag)
#  endif
        IF (exit_flag.ne.NoError) RETURN
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
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL metrics (ng, TILE, iADM)
# if defined PROPAGATOR || \
    (defined MASKING && (defined READ_WATER || defined WRITE_WATER ))
            CALL wpoints (ng, TILE, iADM)
# endif
          END DO
        END DO
!$OMP END PARALLEL DO

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
            CALL ana_nudgcoef (ng, TILE, iADM)
          END DO
        END DO
!$OMP END PARALLEL DO
# endif
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
          CALL ini_hmixcoef (ng, TILE, iADM)
        END DO
      END DO
!$OMP END PARALLEL DO

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
          CALL ana_hmixcoef (ng, TILE, iADM)
        END DO
      END DO
!$OMP END PARALLEL DO
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
!  Clear all adjoint variables.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL initialize_ocean (ng, TILE, iADM)
#  if defined SOLVE3D
          CALL initialize_coupling (ng, TILE, iADM)
          CALL initialize_mixing (ng, TILE, iADM)
#  endif
          CALL initialize_forces (ng, TILE, iADM)
          CALL initialize_grid (ng, TILE, iADM)
#  ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, TILE, iADM)
#  endif
        END DO
      END DO
!$OMF END PARALLEL DO

# elif defined FOUR_DVAR && !defined IS4DVAR_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Clear all adjoint variables.  In variational data assimilation the
!  initial condition are always zero and the forcing is only via the
!  (model-observations) misfit terms.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL initialize_ocean (ng, TILE, iNLM)
          CALL initialize_ocean (ng, TILE, iADM)
#  if defined SOLVE3D
          CALL initialize_coupling (ng, TILE, iADM)
          CALL initialize_mixing (ng, TILE, iADM)
#  endif
          CALL initialize_forces (ng, TILE, iADM)
          CALL initialize_grid (ng, TILE, iADM)
#  ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, TILE, iADM)
#  endif
        END DO
      END DO
!$OMF END PARALLEL DO

# else

#  if defined SOLVE3D && !defined INI_FILE
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_initial (ng, TILE, iADM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#  endif

#  if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_passive (ng, TILE, iADM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#  endif

#  if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology tracers.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_biology (ng, TILE, iADM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#  endif

#  if defined ANA_SEDIMENT_NOT_YET && defined SOLVE3D
!
!  Analytical initial conditions for sediment tracers.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_sediment (ng, TILE, iADM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
#  endif
!
!  Read in initial conditions for initial or restart NetCDF file.
!
#  ifdef INI_FILE
      CALL get_state (ng, iADM, 1, IADname(ng), IniRec, Tindex)
#   ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iADM, exit_flag)
#   endif
      IF (exit_flag.ne.NoError) RETURN
#  else
      IF (nrrec(ng).ne.0) THEN
        CALL get_state (ng, iADM, 1, IADname(ng), IniRec, Tindex)
#   ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, exit_flag)
#   endif
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
# endif

# if defined ANA_PERTURB  && \
    (defined SANITY_CHECK || defined R_SYMMETRY)
!
!-----------------------------------------------------------------------
!  Perturb adjoint initial conditions with analitical expressions.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL ana_perturb (ng, TILE, iADM)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# ifdef SOLVE3D
!!
!!----------------------------------------------------------------------
!!  Compute initial time-evolving depths.
!!----------------------------------------------------------------------
!!
!!$OMP PARALLEL DO PRIVATE(thread,subs)                                 &
!!$OMP&            SHARED(ng,numthreads)
!!    DO thread=0,numthreads-1
!!      subs=NtileX(ng)*NtileE(ng)/numthreads
!!      DO tile=subs*thread,subs*(thread+1)-1
!!        CALL ad_set_depth (ng, TILE)
!!      END DO
!!    END DO
!!$OMP END PARALLEL DO
!!
!!----------------------------------------------------------------------
!!  Compute initial horizontal mass fluxes, Hz*u/n and Hz*v/m.
!!----------------------------------------------------------------------
!!
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
!!    DO thread=0,numthreads-1
!!      subs=NtileX(ng)*NtileE(ng)/numthreads
!!      DO tile=subs*thread,subs*(thread+1)-1
!!        CALL ad_set_massflux (ng, TILE)
!!      END DO
!!    END DO
!!$OMP END PARALLEL DO
!!
!!----------------------------------------------------------------------
!!  Compute initial S-coordinates vertical velocity. Compute initial
!!  density anomaly from potential temperature and salinity via equation
!!  of state for seawater.  Also compute other equation of state related
!!  quatities.
!!----------------------------------------------------------------------
!!
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
!!    DO thread=0,numthreads-1
!!      subs=NtileX(ng)*NtileE(ng)/numthreads
!!      DO tile=subs*thread,subs*(thread+1)-1
!!        CALL ad_omega (ng, TILE)
!!        CALL ad_rho_eos (ng, TILE)
!!      END DO
!!    END DO
!!!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
# ifdef TIMELESS_DATA
      CALL ad_get_idata (ng)
# endif
      CALL ad_get_data (ng)
#  ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iADM, exit_flag)
#  endif
      IF (exit_flag.ne.NoError) RETURN

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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL adsen_initial (ng, TILE)
          END DO
        END DO
!$OMF END PARALLEL DO
#  endif
#  ifdef SENSITIVITY_4DVAR
      END IF
#  endif
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
            CALL stiffness (ng, TILE, iADM)
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
      CALL grid_coords (ng, iADM)
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
        CALL wclock_off (ng, iADM, 1)
      END DO
!$OMP END PARALLEL DO
# endif
      RETURN
      END SUBROUTINE ad_initial
#else
      SUBROUTINE ad_initial
      RETURN
      END SUBROUTINE ad_initial
#endif
