#include "cppdefs.h"
#if defined TANGENT && !defined SOLVE3D
      SUBROUTINE tl_main2d (RunInterval)
!
!git $Id$
!svn $Id: tl_main2d.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the main driver for tangent linear  ROMS/TOMS when  !
!  configurated as shallow water  (barotropic ) ocean model only.  It  !
!  advances forward  the  tangent linear model  for all nested grids,  !
!  if any, by the specified time interval (seconds), RunInterval.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined MODEL_COUPLING && defined MCT_LIB
      USE mod_coupler
# endif
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
      USE dateclock_mod,        ONLY : time_string
# ifdef TLM_CHECK
      USE dotproduct_mod,       ONLY : tl_dotproduct
# endif
# ifdef TIDE_GENERATING_FORCES
      USE equilibrium_tide_mod, ONLY : equilibrium_tide
# endif
# if defined ATM_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod,    ONLY : atmos_coupling
# endif
# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
      USE ocean_coupler_mod,    ONLY : waves_coupling
# endif
      USE strings_mod,          ONLY : FoundError
      USE tl_diag_mod,          ONLY : tl_diag
# if defined WEAK_CONSTRAINT || defined FORCING_SV
      USE tl_forcing_mod,       ONLY : tl_forcing
# endif
# ifdef ADJUST_WSTRESS
      USE tl_frc_adjust_mod,    ONLY : tl_frc_adjust
# endif
      USE tl_ini_fields_mod,    ONLY : tl_ini_fields, tl_ini_zeta
# ifdef ADJUST_BOUNDARY
      USE tl_obc_adjust_mod,    ONLY : tl_obc_adjust
# endif
# ifdef NEARSHORE_MELLOR_NOT_YET
!!    USE tl_radiation_stress_mod, ONLY : tl_radiation_stress
# endif
# ifdef TL_AVERAGES
      USE tl_set_avg_mod,       ONLY : tl_set_avg
# endif
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE tl_set_tides_mod,     ONLY : tl_set_tides
# endif
      USE tl_set_vbc_mod,       ONLY : tl_set_vbc
      USE tl_step2d_mod,        ONLY : tl_step2d
# ifdef FLOATS_NOT_YET
!!    USE tl_step_floats_mod,   ONLY : tl_step_floats
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: ng, tile
      integer :: next_indx1
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
!
      real(r8) :: MaxDT, my_StepTime
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Time-step tangent linear vertically integrated equations.
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
          iic(ng)=iic(ng)+1
!$OMP MASTER
          time(ng)=time(ng)+dt(ng)
          tdays(ng)=time(ng)*sec2day
          CALL time_string (time(ng), time_code(ng))
!$OMP END MASTER
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP MASTER
          CALL tl_get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tl_set_data (ng, tile)
          END DO
!$OMP BARRIER
        END DO
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if (defined WEAK_CONSTRAINT || defined FORCING_SV) && \
     !defined SP4DVAR
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution impulse forcing to
!  the representer model solution. Notice that the forcing is only
!  needed after finishing all inner loops. The forcing is continuous.
!  That is, it is time interpolated at every time-step from available
!  snapshots (FrequentImpulse=TRUE).
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
#  ifdef WEAK_CONSTRAINT
#   ifdef WEAK_NOINTERP
          IF ((iic(ng).gt.1).and.(iic(ng).ne.ntend(ng)+1).and.          &
     &        (MOD(iic(ng)-1,nadj(ng)).eq.0)) THEN
            IF (Master) THEN
              WRITE (stdout,*) ' FORCING TLM at iic = ',iic(ng)
            END IF
#   endif
            IF (FrequentImpulse(ng)) THEN
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL tl_forcing (ng, tile, kstp(ng), nstp(ng))
              END DO
!$OMP BARRIER
            END IF
#   ifdef WEAK_NOINTERP
          END IF
#   endif
#  else
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tl_forcing (ng, tile, kstp(ng), nstp(ng))
          END DO
!$OMP BARRIER
#  endif
        END DO
# endif
!
!-----------------------------------------------------------------------
!  If not a restart, initialize all time levels and compute other
!  initial fields.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface.
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL tl_ini_zeta (ng, tile, iTLM)
            END DO
!$OMP BARRIER
!
!  Initialize other state variables.
!
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL tl_ini_fields (ng, tile, iTLM)
            END DO
!$OMP BARRIER
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Compute and report diagnostics. If appropriate, accumulate time-
!  averaged output data which needs a irreversible loop in shared-memory
!  jobs.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1         ! irreversible
# ifdef TL_AVERAGES
            CALL tl_set_avg (ng, tile)
# endif
# ifdef DIAGNOSTICS
!!          CALL set_diags (ng, tile)
# endif
# ifdef TIDE_GENERATING_FORCES
            CALL equilibrium_tide (ng, tile, iTLM)
# endif
            CALL tl_diag (ng, tile)
# ifdef TLM_CHECK
            CALL tl_dotproduct (ng, tile, Lnew(ng))
# endif
          END DO
!$OMP BARRIER
        END DO

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
!$OMP BARRIER
          END IF
        END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Interpolate open boundary increments and adjust open boundaries.
!  Skip the last output timestep.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).lt.(ntend(ng)+1)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL tl_obc_adjust (ng, tile, Lbinp(ng))
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Skip the last output timestep.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (iic(ng).lt.(ntend(ng)+1)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL tl_frc_adjust (ng, tile, Lfinp(ng))
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif

# if defined WAV_COUPLING_NOT_YET && defined MCT_LIB
!
!-----------------------------------------------------------------------
!  Couple to waves model every CoupleSteps(Iwaves) timesteps: get
!  waves/sea fluxes.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL waves_coupling (ng, tile)
            END DO
!$OMP BARRIER
          END IF
        END DO
# endif

# ifdef NEARSHORE_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL tl_radiation_stress (ng, tile)
          END DO
!$OMP BARRIER
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Set vertical boundary conditions. Process tidal forcing.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tl_set_vbc (ng, tile)
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
            CALL tl_set_tides (ng, tile)
# endif
          END DO
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP MASTER
          CALL tl_output (ng)
!$OMP END MASTER
!$OMP BARRIER
          IF ((FoundError(exit_flag, NoError, __LINE__, MyFile)).or.    &
     &        ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and momentum components.
!-----------------------------------------------------------------------
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
        DO ng=1,Ngrids
          iif(ng)=1
          nfast(ng)=1
          next_indx1=3-indx1(ng)
          IF (.not.PREDICTOR_2D_STEP(ng)) THEN
            PREDICTOR_2D_STEP(ng)=.TRUE.
            IF (FIRST_2D_STEP) THEN
              kstp(ng)=indx1(ng)
            ELSE
              kstp(ng)=3-indx1(ng)
            END IF
            knew(ng)=3
            krhs(ng)=indx1(ng)
          END IF
!
!  Predictor step - Advance barotropic equations using 2D time-step
!  ==============   predictor scheme.
!
          DO tile=last_tile(ng),first_tile(ng),-1
            CALL tl_step2d (ng, tile)
          END DO
!$OMP BARRIER
        END DO
!
!  Set time indices for corrector step.
!
        DO ng=1,Ngrids
          IF (PREDICTOR_2D_STEP(ng)) THEN
            PREDICTOR_2D_STEP(ng)=.FALSE.
            knew(ng)=next_indx1
            kstp(ng)=3-knew(ng)
            krhs(ng)=3
            IF (iif(ng).lt.(nfast(ng)+1)) indx1(ng)=next_indx1
          END IF
!
!  Corrector step - Apply 2D time-step corrector scheme.  Notice that
!  ==============   there is not need for a corrector step during the
!  auxiliary (nfast+1) time-step.
!
          IF (iif(ng).lt.(nfast(ng)+1)) THEN
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL tl_step2d (ng, tile)
            END DO
!$OMP BARRIER
          END IF
        END DO

# ifdef FLOATS_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories: Split all the drifters
!  between all the computational threads, except in distributed-memory
!  and serial configurations. In distributed-memory, the parallel node
!  containing the drifter is selected internally since the state
!  variables do not have a global scope.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (Lfloats(Ng)) THEN
#  ifdef _OPENMP
            chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
            Lstr=1+my_thread*chunk_size
            Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
#  else
            Lstr=1
            Lend=Nfloats(ng)
#  endif
            CALL tl_step_floats (ng, Lstr, Lend)
!$OMP BARRIER
!
!  Shift floats time indices.
!
            nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
            nf  (ng)=MOD(nf  (ng)+1,NFT+1)
            nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
            nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
            nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
          END IF
        END DO
# endif
      END DO STEP_LOOP
!
      RETURN
      END SUBROUTINE tl_main2d
#else
      SUBROUTINE tl_main2d
      RETURN
      END SUBROUTINE tl_main2d
#endif
