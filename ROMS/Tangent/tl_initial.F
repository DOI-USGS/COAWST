#include "cppdefs.h"
#ifdef TANGENT
      SUBROUTINE tl_initial (ng)
!
!svn $Id: tl_initial.F 652 2008-07-24 23:20:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
# ifdef ANALYTICAL
      USE analytical_mod
# endif
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
# if defined GRADIENT_CHECK || defined TLM_CHECK
      USE ini_adjust_mod, ONLY : tl_ini_perturb
# endif
# ifdef OBS_SENSITIVITY
      USE ini_lanczos_mod, ONLY : ini_lanczos
# endif
      USE metrics_mod, ONLY : metrics
# if defined W4DPSAS || defined W4DVAR
      USE tl_set_depth_mod, ONLY : tl_bath
# endif
# ifdef SOLVE3D
      USE set_depth_mod, ONLY : set_depth
      USE tl_set_depth_mod, ONLY : tl_set_depth
      USE tl_omega_mod, ONLY : tl_omega
      USE tl_rho_eos_mod, ONLY : tl_rho_eos
      USE tl_set_massflux_mod, ONLY : tl_set_massflux
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
      logical, save :: FirstPass = .FALSE.

      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: update = .FALSE.

      integer :: LBi, UBi, LBj, UBj
      integer :: IniRec, InpRec, Tindex, subs, tile, thread, wrtRec

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
# elif defined IS4DVAR || defined IS4DVAR_OLD || defined W4DPSAS || \
       defined W4DVAR
        WRITE (stdout,10) outer, inner
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Outer = ',i3.3, ', Inner = ',i3.3,' >>>>',/)
# elif defined S4DVAR
        WRITE (stdout,10) Nrun, Ipass
 10     FORMAT (/,' <<<< 4D Variational Data Assimilation, ',           &
     &          'Iteration = ',i5.5, ', Ipass = ',i1,' >>>>',/)
# endif
        WRITE (stdout,20) 'TL_INITIAL: Configuring and ',               &
     &                    'initializing tangent linear model ...'
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
        CALL wclock_on (ng, iTLM, 1)
      END DO
!$OMP END PARALLEL DO
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

# elif defined FOUR_DVAR || defined GRADIENT_CHECK || defined TLM_CHECK
!
!-----------------------------------------------------------------------
!  If variational data assimilation, reset several IO switches and
!  variables.
!-----------------------------------------------------------------------
#  ifndef OBS_SENSITIVITY
!
!  Set switches to create (TRUE) tangent linear initial conditions and
!  history NetCDF files or append (FALSE) to existing NetCDF files.
!
#   if defined WEAK_CONSTRAINT
      IF (Nrun.eq.ERstr) THEN
#   elif defined IS4DVAR
      IF ((Nrun.eq.ERstr).and.(inner.eq.0)) THEN
#   else    
      IF ((Nrun.eq.ERstr).and.((inner.eq.1).and.(Ipass.le.1))) THEN
#   endif
        LdefITL(ng)=.TRUE.
#   ifndef WEAK_CONSTRAINT
        LdefTLM(ng)=.TRUE.
#   endif
!
!  It first pass, create tangent linear model initialization file and
!  write zero initial conditions for records 1 and 2.
!
        CALL tl_def_ini (ng)
        IF (exit_flag.ne.NoError) RETURN
        CALL tl_wrt_ini (ng, Tindex, 1)
        IF (exit_flag.ne.NoError) RETURN
        CALL tl_wrt_ini (ng, Tindex, 2)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Set record to read from initial tangent linear NetCDF file.
!
      IniRec=tITLindx(ng)

#   if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Initialize surface forcing counter for storage arrays.
!
      SFcount(ng)=0
#   endif
#   ifndef WEAK_CONSTRAINT
!
!  Reset tangent linear model history time record counters. These
!  counters are reset in every iteration pass. This file is created
!  on the first iteration pass.
!
      tTLMindx(ng)=0
      NrecTLM(ng)=0

#    ifdef IS4DVAR_OLD
!  
!  Avoid writing tangent linear history file on the second pass
!  (Ipass=2).
!
      IF (Ipass.eq.1) THEN
        LwrtTLM(ng)=.TRUE.
      ELSEif (Ipass.eq.2) THEN
        LwrtTLM(ng)=.FALSE.
      END IF
#    else
      LwrtTLM(ng)=.TRUE.
#    endif
!
!  Insure that forward and history file names are the same. In 4DVar,
!  the forward solution is computed by the nonlinear model and stored
!  on HISNAME NetCDF file.
!
#    ifdef TLM_CHECK
      FWDname(ng)=HISbase(ng)
      ncFWDid(ng)=-1
#    else
      FWDname(ng)=HISname(ng)
      ncFWDid(ng)=ncHISid(ng)
#    endif
#   else
!
!  If weak constraint 4DVAR, activate writting of history NetCDF
!  file.
!
!!    LdefTLM(ng)=.TRUE.
!!    LwrtTLM(ng)=.TRUE.
#   endif
#  endif
!
!  Open observations NetCDF file and initialize various variables
!  needed for processing the nonlinear state solution at observation
!  locations.
!
      CALL obs_initial (ng, iTLM, .FALSE.)
#  ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag, 1)
#  endif
      IF (exit_flag.ne.NoError) RETURN
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
            CALL ana_grid (ng, TILE, iTLM)
#  ifdef MASKING
            CALL ana_mask (ng, TILE, iTLM)
#  endif
          END DO
        END DO
!$OMP END PARALLEL DO
# else
        CALL get_grid (ng, iTLM)
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iTLM, exit_flag, 1)
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
            CALL metrics (ng, TILE, iTLM)
# if defined PROPAGATOR || \
    (defined MASKING && (defined READ_WATER || defined WRITE_WATER ))
            CALL wpoints (ng, TILE, iTLM)
# endif
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
            CALL ana_hmixcoef (ng, TILE, iTLM)
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
            CALL ana_nudgcoef (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
# endif
      END IF
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration pass.
!=======================================================================

# if defined FOUR_DVAR && !defined OBS_SENSITIVITY
#  if defined GRADIENT_CHECK || defined OPT_OBSERVATIONS || \
      defined TLM_CHECK      || defined WEAK_CONSTRAINT
!
!  Clear tangent linear state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL initialize_ocean (ng, TILE, iTLM)
#   ifdef SOLVE3D
          CALL initialize_coupling (ng, TILE, 0)
#   endif
        END DO
      END DO
!$OMP END PARALLEL DO

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
#    ifdef IS4DVAR_OLD
      IF ((inner.eq.1).and.(Ipass.le.1)) THEN
#    else
      IF (inner.eq.0) THEN
#    endif
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL initialize_ocean (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Rewrite tangent linear initial NetCDF (record 1) with zero initial
!  conditions since the model needs to be started from at the first
!  pass of the inner loop.
!
        IF (Nrun.gt.1) THEN
          wrtRec=1
          CALL tl_wrt_ini (ng, Tindex, wrtRec)
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
!$OMP PARALLEL DO PRIVATE(thread,subs)                                  &
!$OMP&            SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL set_depth (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# if defined W4DPSAS || defined W4DVAR
!
!-----------------------------------------------------------------------
!  Initialize tangent linear bathymetry to zero.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs)                                  &
!$OMP&            SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL tl_bath (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_initial (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# if defined ANA_PASSIVE && defined SOLVE3D
!
!  Analytical initial conditions for inert passive tracers.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_passive (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# if defined ANA_BIOLOGY && defined SOLVE3D
!
!  Analytical initial conditions for biology.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_biology (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# if defined ANA_SEDIMENT_NOT_YET && defined SOLVE3D
!
!  Analytical initial conditions for sediment.
!
      IF (nrrec(ng).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ana_sediment (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# ifdef OBS_SENSITIVITY
!
!  Initialize with the weighted sum of all Lanczos vectors computed
!  from the first outer loop of the IS4DVAR Lanczos algorithm.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile,Tindex)                      &
!$OMP&            SHARED(ng,numthreads,Lnew)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL ini_lanczos (ng, TILE, Lnew(ng), Tindex)
        END DO
      END DO
!$OMP END PARALLEL DO

# else
!
!  Read in initial conditions for initial or restart NetCDF file.
!
#  ifdef INI_FILE
      CALL get_state (ng, iTLM, 1, ITLname(ng), IniRec, Tindex)
      IF (exit_flag.ne.NoError) RETURN
#  else
      IF (nrrec(ng).ne.0) THEN
        CALL get_state (ng, iTLM, 1, ITLname(ng), IniRec, Tindex)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
# endif

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
          CALL ana_perturb (ng, TILE, iTLM)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# if defined GRADIENT_CHECK || defined TLM_CHECK
!
!-----------------------------------------------------------------------
!  Perturb tangent linear state variable according to the outer loop
!  iteration with the steepest descent direction of the gradient
!  (adjoint state).
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile,Tindex)                      &
!$OMP             SHARED(ng,numthreads,Lnew)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1
          CALL tl_ini_perturb (ng, TILE, Lnew(ng), Tindex)
        END DO
      END DO
!$OMP END PARALLEL DO
#  endif

# ifdef SOLVE3D
!!
!!----------------------------------------------------------------------
!!  Compute initial time-evolving depths.
!!----------------------------------------------------------------------
!!
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                            &
!!$OMP&            SHARED(ng,numthreads)
!!    DO thread=0,numthreads-1
!!      subs=NtileX(ng)*NtileE(ng)/numthreads
!!      DO tile=subs*thread,subs*(thread+1)-1
!!        CALL tl_set_depth (ng, TILE)
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
!!        CALL tl_set_massflux (ng, TILE)
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
!!        CALL tl_omega (ng, TILE)
!!        CALL tl_rho_eos (ng, TILE)
!!      END DO
!!    END DO
!!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
# ifdef TIMELESS_DATA
      CALL tl_get_idata (ng)
# endif
      CALL tl_get_data (ng)
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag, 1)
# endif
      IF (exit_flag.ne.NoError) RETURN

# if defined WEAK_CONSTRAINT || defined IOM
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
        FrcRec(ng)=1
        CALL get_state (ng, 7, 7, TLFname(ng), FrcRec(ng), 1)
        IF (exit_flag.ne.NoError) RETURN
      END IF
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
            CALL stiffness (ng, TILE, iTLM)
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
      CALL grid_coords (ng, iTLM)
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
        CALL wclock_off (ng, iTLM, 1)
      END DO
!$OMP END PARALLEL DO
# endif
      RETURN
      END SUBROUTINE tl_initial
#else
      SUBROUTINE tl_initial
      RETURN
      END SUBROUTINE tl_initial
#endif
