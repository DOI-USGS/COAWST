#include "cppdefs.h"
#if defined ADJOINT && !defined SOLVE3D
      SUBROUTINE ad_main2d (RunInterval)
!
!git $Id$
!svn $Id: ad_main2d.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine is the main driver for adjoint  ROMS/TOMS  when     !
!  configurated as shallow water (barotropic) ocean model only. It     !
!  advances backward  the  adjoint model  for all nested grids, if     !
!  any, by the specified time interval (seconds), RunInterval.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined MODEL_COUPLING && defined MCT_LIB
      USE mod_coupler
# endif
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
# ifdef SO_SEMI
      USE mod_storage
# endif
!
# if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR
      USE adsen_force_mod,         ONLY : adsen_force
# endif
      USE ad_diag_mod,             ONLY : ad_diag
# if defined WEAK_CONSTRAINT || defined FORCING_SV
      USE ad_forcing_mod,          ONLY : ad_forcing
# endif
# ifdef ADJUST_WSTRESS
      USE ad_frc_adjust_mod,       ONLY : ad_frc_adjust
# endif
      USE ad_ini_fields_mod,       ONLY : ad_ini_fields, ad_ini_zeta
# ifdef AD_OUTPUT_STATE
      USE ad_ini_fields_mod,       ONLY : ad_out_fields, ad_out_zeta
# endif
# if defined FOUR_DVAR && defined OBSERVATIONS
#  ifdef WEAK_CONSTRAINT
      USE ad_htobs_mod,            ONLY : ad_htobs
#  else
      USE ad_misfit_mod,           ONLY : ad_misfit
#  endif
# endif
# ifdef ADJUST_BOUNDARY
      USE ad_obc_adjust_mod,       ONLY : ad_obc_adjust
# endif
# ifdef NEARSHORE_MELLOR_NOT_YET
!!    USE ad_radiation_stress_mod, ONLY : ad_radiation_stress
# endif
# ifdef AD_AVERAGES
      USE ad_set_avg_mod,          ONLY : ad_set_avg
# endif
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE ad_set_tides_mod,        ONLY : ad_set_tides
# endif
      USE ad_set_vbc_mod,          ONLY : ad_set_vbc
      USE ad_step2d_mod,           ONLY : ad_step2d
# ifdef FLOATS_NOT_YET
!!    USE ad_step_floats_mod,      ONLY : tl_step_floats
# endif
      USE dateclock_mod,           ONLY : time_string
# ifdef TIDE_GENERATING_FORCES
      USE equilibrium_tide_mod,    ONLY : equilibrium_tide
# endif
# ifdef WEAK_CONSTRAINT
      USE frc_weak_mod,            ONLY : frc_ADgather, frc_clear
# endif
# if (defined FOUR_DVAR    && !defined I4DVAR_ANA_SENSITIVITY) && \
      defined OBSERVATIONS
      USE obs_read_mod,            ONLY : obs_read
# endif
# ifdef AIR_OCEAN_NOT_YET
      USE ocean_coupler_mod,       ONLY : atmos_coupling
# endif
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod,       ONLY : waves_coupling
# endif
# ifdef SO_SEMI
      USE packing_mod,             ONLY : ad_pack
# endif
      USE strings_mod,             ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical:: backward = .TRUE.
!
      integer :: ng, tile
      integer :: next_indx1
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
      integer :: ksav, ktmp
!
      real(r8) :: HalfDT, MaxDT, my_StepTime
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Time-step adjoint vertically integrated equations.
!=======================================================================
!
      my_StepTime=0.0_r8
      MaxDT=MAXVAL(dt)

      STEP_LOOP : DO WHILE (my_StepTime.le.(RunInterval+0.5_r8*MaxDT))

        my_StepTime=my_StepTime+MaxDT
!
!  Set time clock.
!
        DO ng=1,Ngrids
          iic(ng)=iic(ng)-1
          time(ng)=time(ng)-dt(ng)
          tdays(ng)=time(ng)*sec2day
          CALL time_string (time(ng), time_code(ng))
        END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL ad_get_data (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Process input data, if any: time interpolate between snapshots.
!  If appropriate, compute and report diagnostics.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1         ! irreversible
            CALL ad_set_data (ng, tile)
# ifdef AD_AVERAGES
            CALL ad_set_avg (ng, tile)
# endif
# ifdef TIDE_GENERATING_FORCES
            CALL equilibrium_tide (ng, tile, iADM)
# endif
# ifdef DIAGNOSTICS
!!          CALL ad_set_diags (ng, tile)
# endif
            CALL ad_diag (ng, tile)
          END DO
        END DO
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if (defined FOUR_DVAR    && !defined I4DVAR_ANA_SENSITIVITY) && \
      defined OBSERVATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, read observation and model state at observation
!  locations.  Then, compute adjoint misfit forcing terms.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
#  ifdef SENSITIVITY_4DVAR
          IF (.not.Lsen4DVAR(ng)) THEN
#  endif
            HalfDT=0.5_r8*dt(ng)
            IF (((time(ng)-HalfDT).le.ObsTime(ng)).and.                 &
     &          (ObsTime(ng).lt.(time(ng)+HalfDT))) THEN
              ProcessObs(ng)=.TRUE.
              CALL obs_read (ng, iADM, backward)
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            ELSE
              ProcessObs(ng)=.FALSE.
            END IF
!
#  ifdef SP4DVAR
!
!  Skip assimilation of obs on first timestep unless inter=Nsaddle.
!
          IF ((iic(ng).ne.ntstart(ng)).or.Lsadd(ng)) THEN
#  endif
            DO tile=first_tile(ng),last_tile(ng),+1
#  ifdef WEAK_CONSTRAINT
              CALL ad_htobs (ng, tile, iADM)
#  else
              CALL ad_misfit (ng, tile, iADM)
#  endif
            END DO
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  ifdef SENSITIVITY_4DVAR
          END IF
#  endif
#  ifdef SP4DVAR
          END IF
#  endif
        END DO
# endif

# ifdef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  If appropriate, add representer coefficients (Beta hat) impulse
!  forcing to adjoint solution. Read next impulse record, if available.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (ProcessObs(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ad_forcing (ng, tile, knew(ng), nnew(ng))
            END DO
          END IF
        END DO
# endif
!
!  Avoid time-stepping if additional delayed IO time-step.
!
        TIME_STEP: IF (MINVAL(iic).ne.MINVAL(ntstart)) THEN

# ifdef FLOATS_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories.
!-----------------------------------------------------------------------
!
!  Shift floats time indices.
!
          DO ng=1,Ngrids
            IF (Lfloats(Ng)) THEN
              nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
              nf  (ng)=MOD(nf  (ng)+1,NFT+1)
              nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
              nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
              nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
!
#  ifdef _OPENMP
              chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
              Lstr=1+MyThread*chunk_size
              Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
#  else
              Lstr=1
              Lend=Nfloats(ng)
#  endif
              CALL ad_step_floats (ng, Lstr, Lend)
            END IF
          END DO
# endif
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and momentum components.
!-----------------------------------------------------------------------
!
!  Corrector step - Apply 2D time-step corrector scheme.  Notice that
!  ==============   there is not need for a corrector step during the
!  auxiliary (nfast+1) time-step.
!
          DO ng=1,Ngrids
            iif(ng)=1
            nfast(ng)=1
            IF (iif(ng).lt.(nfast(ng)+1)) THEN
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL ad_step2d (ng, tile)
              END DO
            END IF
!
!  Set time indices for corrector step.
!
            next_indx1=3-indx1(ng)
            IF (.not.PREDICTOR_2D_STEP(ng)) THEN
              PREDICTOR_2D_STEP(ng)=.TRUE.
              ktmp=knew(ng)
              ksav=kstp(ng)
              knew(ng)=krhs(ng)
              kstp(ng)=ktmp
              krhs(ng)=ksav
            END IF
          END DO
!
!  Predictor step - Advance barotropic equations using 2D time-step
!  ==============   predictor scheme.
!
          DO ng=1,Ngrids
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL ad_step2d (ng, tile)
            END DO
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
            IF (PREDICTOR_2D_STEP(ng).and.                              &
     &          (iic(ng).ne.ntend(ng))) THEN
              PREDICTOR_2D_STEP(ng)=.FALSE.
              ksav=knew(ng)
              knew(ng)=krhs(ng)
              krhs(ng)=ksav
            END IF
          END DO

        END IF TIME_STEP

# ifdef SO_SEMI
!
!-----------------------------------------------------------------------
!  If stochastic optimals with respect the seminorm of chosen
!  functional, pack adjoint state surface forcing needed by the
!  dynamical propagator.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (MOD(iic(ng)-1,nADJ(ng)).eq.0) THEN
            SOrec(ng)=SOrec(ng)+1
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ad_pack (ng, tile, Nstr(ng), Nend(ng),               &
     &                      STORAGE(ng)%so_state(:,SOrec(ng)))
            END DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Set vertical boundary conditions. Process tidal forcing.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
            CALL ad_set_tides (ng, tile)
# endif
            CALL ad_set_vbc (ng, tile)
          END DO
        END DO

# ifdef NEARSHORE_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL ad_radiation_stress (ng, tile)
          END DO
        END DO
# endif

# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple to waves model every CoupleSteps(Iwaves,ng) timesteps: get
!  waves/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL waves_coupling (ng, tile)
            END DO
          END IF
        END DO
# endif

# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple to atmospheric model every CoupleSteps(Iatmos) timesteps: get
!  air/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng)-1,CoupleSteps(Iatmos,ng)).eq.0) THEN
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL atmos_coupling (ng, tile)
            END DO
          END IF
        END DO
# endif
# ifdef AD_OUTPUT_STATE
!
!-----------------------------------------------------------------------
!  Set full adjoint output arrays. Due to the predictor/corrector and
!  multiple time level schemes, pieces of the adjoint solution are in
!  two-time levels and need to be added in the "_sol" arrays for output
!  purposes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).ne.ntend(ng)) THEN
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL ad_out_fields (ng, tile, iADM)
            END DO
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ad_out_zeta (ng, tile, iADM)
            END DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  If not a restart, initialize all time levels and compute other
!  initial fields.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).eq.ntend(ng)) THEN
!
!  Initialize other state variables.
!
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL ad_ini_fields (ng, tile, iADM)
            END DO
!
!  Initialize free-surface.
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ad_ini_zeta (ng, tile, iADM)
            END DO
          END IF
        END DO

# ifdef FORCING_SV
!
!-----------------------------------------------------------------------
!  Compute adjoint forcing for forcing singular vectors.
!-----------------------------------------------------------------------
!
        IF (iic(ng).ne.ntend(ng)) THEN
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ad_forcing (ng, tile, kstp(ng), knew(ng), nstp(ng))
            END DO
          END DO
        ELSE
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ad_forcing (ng, tile, kstp(ng), krhs(ng), nstp(ng))
            END DO
          END DO
        END IF
# endif

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Skip first timestep.
!-----------------------------------------------------------------------
!
#  ifdef RBL4DVAR_FCT_SENSITIVITY
        DO ng=1,Ngrids
          IF (.not.Lsen4DVAR(ng)) THEN     ! ignore in forecast interval
            IF (iic(ng).ne.ntstart(ng)) THEN
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL ad_frc_adjust (ng, tile, Lfout(ng))
              END DO
            END IF
          END IF
        END DO
#  else
        DO ng=1,Ngrids
          IF (iic(ng).ne.ntstart(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ad_frc_adjust (ng, tile, Lfout(ng))
            END DO
          END IF
        END DO
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundaries.
!  Skip first timestep.
!-----------------------------------------------------------------------
!
#  ifdef RBL4DVAR_FCT_SENSITIVITY
        DO ng=1,Ngrids
          IF (.not.Lsen4DVAR(ng)) THEN     ! ignore in forecast interval
            IF (iic(ng).ne.ntstart(ng)) THEN
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL ad_obc_adjust (ng, tile, Lbout(ng))
              END DO
            END IF
          END IF
        END DO
#  else
        DO ng=1,Ngrids
          IF (iic(ng).ne.ntstart(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ad_obc_adjust (ng, tile, Lbout(ng))
            END DO
          END IF
        END DO
#  endif
# endif

# if defined WEAK_CONSTRAINT && !defined SP4DVAR
!
!-----------------------------------------------------------------------
!  Gather weak constraint forcing to storage arrays.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).ne.ntstart(ng)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL frc_ADgather (ng, tile)
            END DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          CALL ad_output (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO

# if defined WEAK_CONSTRAINT && !defined SP4DVAR
!
!-----------------------------------------------------------------------
!  Copy storage arrays index 1 into index 2, and clear index 1.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (MOD(iic(ng)-1,nADJ(ng)).eq.0) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL frc_clear (ng, tile)
            END DO
          END IF
        END DO
# endif

# if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR
!
!-----------------------------------------------------------------------
!  Add appropriate forcing terms to the adjoint model. The form of the
!  forcing depends on the functional whose sensitivity is required.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
#  ifdef SENSITIVITY_4DVAR
          IF (Lsen4DVAR(ng)) THEN
#  endif
#  if !defined AD_IMPULSE
            IF ((DendS(ng).ge.tdays(ng)).and.                           &
     &          (tdays(ng).ge.DstrS(ng))) THEN
#  endif
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL adsen_force (ng, tile)
              END DO
#  if !defined AD_IMPULSE
            END IF
#  endif
#  ifdef SENSITIVITY_4DVAR
          END IF
#  endif
        END DO
# endif
      END DO STEP_LOOP
!
      RETURN
      END SUBROUTINE ad_main2d
#else
      SUBROUTINE ad_main2d
      RETURN
      END SUBROUTINE ad_main2d
#endif
