#include "cppdefs.h"
#if defined TANGENT && !defined SOLVE3D
      SUBROUTINE tl_main2d (ng)
!
!svn $Id: tl_main2d.F 704 2008-08-13 06:20:10Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine is the main driver for tangent linear  ROMS/TOMS when  !
!  configurated as shallow water  (barotropic ) ocean model only.  It  !
!  advances forward the vertically integrated primitive equations for  !
!  a single time step.                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef MODEL_COUPLING
      USE mod_coupler
# endif
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
# ifdef TLM_CHECK
      USE dotproduct_mod, ONLY : tl_dotproduct
# endif
# ifdef ADJUST_WSTRESS
      USE frc_adjust_mod, ONLY : frc_TLadjust, load_TLforcing
# endif
# ifdef AIR_OCEAN_NOT_YET
      USE ocean_coupler_mod, ONLY : atmos_coupling
# endif
# ifdef WAVES_OCEAN_NOT_YET
      USE ocean_coupler_mod, ONLY : waves_coupling
# endif
# ifdef ASSIMILATION
!!    USE oi_update_mod, ONLY : oi_update
# endif
      USE tl_diag_mod, ONLY : tl_diag
# if defined WEAK_CONSTRAINT || defined IOM
      USE tl_forcing_mod, ONLY : tl_forcing
# endif
      USE tl_ini_fields_mod, ONLY : tl_ini_fields, tl_ini_zeta
# ifdef NEARSHORE_MELLOR_NOT_YET
!!    USE tl_radiation_stress_mod, ONLY : tl_radiation_stress
# endif
# ifdef AVERAGES
!!    USE tl_set_avg_mod, ONLY : tl_set_avg
# endif
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
!!    USE tl_set_tides_mod, ONLY : tl_set_tides
# endif
      USE tl_set_vbc_mod, ONLY: tl_set_vbc
      USE tl_step2d_mod, ONLY : tl_step2d
# ifdef FLOATS_NOT_YET
!!    USE tl_step_floats_mod, ONLY : tl_step_floats
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
      integer :: next_indx1, subs, tile, thread
# ifdef FLOATS_NOT_YET
      integer :: Lend, Lstr, chunk_size
# endif
# if defined WEAK_CONSTRAINT || defined IOM
      real(r8) :: HalfDT
# endif
!
!=======================================================================
!  Time-step tangent linear vertically integrated equations.
!=======================================================================
!
!  Set time clock.
!
      time(ng)=time(ng)+dt(ng)
      tdays(ng)=time(ng)*sec2day
      CALL time_string (time(ng), time_code(ng))

#if defined WEAK_CONSTRAINT || defined IOM
!
!-----------------------------------------------------------------------
!  If appropriate, add convolved adjoint solution as impulse forcing to
!  the tangent linear solution. Read next impulse record, if available.
!-----------------------------------------------------------------------
!
      HalfDT=0.5_r8*dt(ng)
      IF (((time(ng)-HalfDT).le.FrcTime(ng)).and.                       &
     &    (FrcTime(ng).lt.(time(ng)+HalfDT))) THEN
!
!  Add impulse forcing.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_forcing (ng, TILE, kstp(ng), nstp(ng))
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Read in impulse forcing record, if any.
!
        FrcRec(ng)=FrcRec(ng)+1
        IF (FrcRec(ng).le.NrecFrc(ng)) THEN
          CALL get_state (ng, 7, 7, TLFname(ng), FrcRec(ng), 1)
        END IF
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
      CALL tl_get_data (ng)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL tl_set_data (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  If not a restart, initialize all time levels and compute other
!  initial fields.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_ini_zeta (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Initialize other state variables.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL tl_ini_fields (ng, TILE, iTLM)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
!
!-----------------------------------------------------------------------
!  Compute and report diagnostics. If appropriate, accumulate time-
!  averaged output data which needs a irreversible loop in shared-memory
!  jobs.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,Lnew,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1     ! irreversible loop
# ifdef AVERAGES
!!        CALL tl_set_avg (ng, TILE)
# endif
# ifdef DIAGNOSTICS_UV
!!        CALL set_diags (ng, TILE)
# endif
          CALL tl_diag (ng, TILE)
# ifdef TLM_CHECK
          CALL tl_dotproduct (ng, TILE, Lnew(ng))
# endif
        END DO
      END DO
!$OMP END PARALLEL DO

# ifdef AIR_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Couple to atmospheric model every CoupleSteps(Iatmos) timesteps: get
!  air/sea fluxes.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).ne.ntstart(ng)).and.                                 &
     &    MOD(iic(ng)-1,CoupleSteps(Iatmos,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL atmos_coupling (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  Interpolate surface forcing increments and adjust surface forcing.
!  Load surface forcing to storage arrays.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,Lfout,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL frc_TLadjust (ng, TILE, Lfinp(ng))
          CALL load_TLforcing (ng, TILE, Lfout(ng))
        END DO
      END DO
!$OMP END PARALLEL DO
# endif

# ifdef WAVES_OCEAN_NOT_YET
!
!-----------------------------------------------------------------------
!  Couple to waves model every CoupleSteps(Iwaves) timesteps: get
!  waves/sea fluxes.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).ne.ntstart(ng)).and.                                 &
     &    MOD(iic(ng)-1,CoupleSteps(Iwaves,ng)).eq.0) THEN
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL waves_coupling (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
# endif

# ifdef NEARSHORE_MELLOR_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute radiation stress terms.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL tl_radiation_stress (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
!  Set vertical boundary conditions. Process tidal forcing.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL tl_set_vbc (ng, TILE)
# if defined SSH_TIDES_NOT_YET || defined UV_TIDES_NOT_YET
          CALL tl_set_tides (ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
      CALL tl_output (ng)
      IF ((exit_flag.ne.NoError).or.(iic(ng).eq.(ntend(ng)+1))) RETURN
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and momentum components.
!-----------------------------------------------------------------------
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*(thread+1)-1,subs*thread,-1
          CALL tl_step2d (ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
!
!  Set time indices for corrector step.
!
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_step2d (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF

# ifdef FLOATS_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute Lagrangian drifters trajectories.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,chunk_size,Lstr,Lend)                  &
!$OMP&            SHARED(ng,numthreads,Nfloats)
      DO thread=0,numthreads-1
        chunk_size=(Nfloats(ng)+numthreads-1)/numthreads
        Lstr=1+thread*chunk_size
        Lend=MIN(Nfloats(ng),Lstr+chunk_size-1)
        CALL tl_step_floats (ng, Lstr, Lend)
      END DO
!$OMP END PARALLEL DO
!
!  Shift floats time indices.
!
      nfp1(ng)=MOD(nfp1(ng)+1,NFT+1)
      nf(ng)  =MOD(nf(ng)  +1,NFT+1)
      nfm1(ng)=MOD(nfm1(ng)+1,NFT+1)
      nfm2(ng)=MOD(nfm2(ng)+1,NFT+1)
      nfm3(ng)=MOD(nfm3(ng)+1,NFT+1)
# endif

      RETURN
      END SUBROUTINE tl_main2d
#else
      SUBROUTINE tl_main2d
      RETURN
      END SUBROUTINE tl_main2d
#endif
