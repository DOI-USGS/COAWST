      SUBROUTINE propagator (RunInterval, state, ad_state)
!
!svn $Id: propagator_hop.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Optimal Perturbation (Hessian singular vectors) Propagator:         !
!                                                                      !
!  This routine solves the singular vectors of the propagator R(0,t)   !
!  which measure the  fastest  growing of all possible perturbations   !
!  over a given time interval. It involves a single integration of a   !
!  perturbation forward in time with the  tangent linear model  over   !
!  [0,t],  followed by an integration of the result backward in time   !
!  with the adjoint model over [t,0].                                  !
!                                                                      !
!   Reference:                                                         !
!                                                                      !
!     Moore, A.M. et al., 2004: A comprehensive ocean prediction and   !
!       analysis system based on the tangent linear and adjoint of a   !
!       regional ocean model, Ocean Modelling, 7, 227-258.             !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
#ifdef SOLVE3D
      USE mod_coupling
#endif
      USE mod_iounits
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE dotproduct_mod, ONLY : tl_statenorm
      USE ini_adjust_mod, ONLY : ad_ini_perturb
      USE inner2state_mod, ONLY : ad_inner2state, tl_inner2state
      USE inner2state_mod, ONLY : ini_C_norm
#ifdef SOLVE3D
      USE set_depth_mod, ONLY: set_depth
#endif
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: RunInterval

      TYPE (T_GST), intent(in) :: state(Ngrids)
      TYPE (T_GST), intent(inout) :: ad_state(Ngrids)
!
!  Local variable declarations.
!
      integer :: ng, tile
      integer :: ktmp, ntmp, Lini

      real(r8) :: StateNorm(Ngrids)
!
!=======================================================================
!  Forward integration of the tangent linear model.
!=======================================================================
!
!$OMP MASTER
      Nrun=Nrun+1
      IF (Master) THEN
        DO ng=1,Ngrids
          WRITE (stdout,10) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Iteration: ', Nrun,                     &
     &                      ',  number converged RITZ values: ',        &
     &                      Nconv(ng)
        END DO
      END IF
!$OMP END MASTER
!
!  Initialize time stepping indices and counters.
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
!
        synchro_flag(ng)=.TRUE.
        tdays(ng)=dstart
        time(ng)=tdays(ng)*day2sec
!$OMP MASTER
        ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
        ntend(ng)=ntimes(ng)
        ntfirst(ng)=ntstart(ng)
!$OMP END MASTER
      END DO
!$OMP BARRIER
!
      Lini=1
!
!-----------------------------------------------------------------------
!  Clear tangent linear state variables. There is not need to clean
!  the basic state arrays since they were zeroth out at initialization
!  and bottom of previous iteration.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iTLM)
        END DO
!$OMP BARRIER
      END DO

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute basic state initial level thicknesses used for state norm
!  scaling. It uses zero time-averaged free-surface (rest state).
!  Therefore, the norm scaling is time invariant.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL set_depth (ng, tile)
        END DO
!$OMP BARRIER
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Compute tangent linear initial conditions from state vector.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL tl_inner2state (ng, tile, Lini, state(ng)%vector)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial tangent linear state analysis error norm.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL ini_C_norm (ng, tile, kstp(ng), nstp(ng),                &
     &                     StateNorm(ng))
        END DO
!$OMP BARRIER

!$OMP MASTER
        IF (Master) THEN
          WRITE (stdout,20) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Tangent Initial Norm: ', StateNorm(ng)
        END IF
!$OMP END MASTER
      END DO
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL close_inp (ng, iTLM)
        IF (exit_flag.ne.NoError) RETURN
        CALL tl_get_idata (ng)
        IF (exit_flag.ne.NoError) RETURN
        CALL tl_get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Time-step the tangent linear model.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        IF (Master) THEN
          WRITE (stdout,30) 'TL', ng, ntstart(ng), ntend(ng)
        END IF
        time(ng)=time(ng)-dt(ng)
!$OMP END MASTER
        iic(ng)=ntstart(ng)-1
      END DO
!$OMP BARRIER

#ifdef SOLVE3D
      CALL tl_main3d (RunInterval)
#else
      CALL tl_main2d (RunInterval)
#endif
!$OMP BARRIER
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Clear nonlinear (basic state) and adjoint state variables.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iNLM)
          CALL initialize_ocean (ng, tile, iADM)
#ifdef SOLVE3D
          CALL initialize_coupling (ng, tile, 0)
#endif
        END DO
!$OMP BARRIER
      END DO

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute basic state final level thicknesses used for state norm
!  scaling. It uses zero time-averaged free-surface (rest state).
!  Therefore, the norm scaling is time invariant.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL set_depth (ng, tile)
        END DO
!$OMP BARRIER
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Compute final tangent linear energy norm.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL tl_statenorm (ng, tile, knew(ng), nstp(ng),              &
     &                       StateNorm(ng))
        END DO
!$OMP BARRIER

!$OMP MASTER
        IF (Master) THEN
          WRITE (stdout,20) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Tangent   Final Norm: ', StateNorm(ng)
        END IF
!$OMP END MASTER
      END DO
!
!=======================================================================
!  Backward integration with the adjoint model.
!=======================================================================
!
!  Initialize time stepping indices and counters.
!
      DO ng=1,Ngrids
        iif(ng)=1
        indx1(ng)=1
        ktmp=knew(ng)
        kstp(ng)=1
        krhs(ng)=3
        knew(ng)=2
        PREDICTOR_2D_STEP(ng)=.FALSE.
!
        iic(ng)=0
        ntmp=nstp(ng)
        nstp(ng)=1
        nrhs(ng)=1
        nnew(ng)=2
!
        synchro_flag(ng)=.TRUE.
        tdays(ng)=dstart+dt(ng)*REAL(ntimes(ng),r8)*sec2day
        time(ng)=tdays(ng)*day2sec
!$OMP MASTER
        ntstart(ng)=ntimes(ng)+1
        ntend(ng)=1
        ntfirst(ng)=ntend(ng)
!$OMP END MASTER
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Initialize adjoint model with the final tangent linear solution
!  scaled by the energy norm.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL ad_ini_perturb (ng, tile,                                &
     &                         ktmp, knew(ng), ntmp, nstp(ng))
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL close_inp (ng, iADM)
        IF (exit_flag.ne.NoError) RETURN
        CALL ad_get_idata (ng)
        IF (exit_flag.ne.NoError) RETURN
        CALL ad_get_data (ng)
!$OMP END MASTER
        IF (exit_flag.ne.NoError) RETURN
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Time-step the adjoint model backwards.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        IF (Master) THEN
          WRITE (stdout,30) 'AD', ng, ntstart(ng), ntend(ng)
        END IF
        time(ng)=time(ng)+dt(ng)
!$OMP END MASTER
        iic(ng)=ntstart(ng)+1
      END DO
!$OMP BARRIER

#ifdef SOLVE3D
      CALL ad_main3d (RunInterval)
#else
      CALL ad_main2d (RunInterval)
#endif
!$OMP BARRIER
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Clear nonlinear state (basic state) variables for next iteration
!  and to insure a rest state time averaged free-surface before adjoint
!  state norm scaling.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iNLM)
#ifdef SOLVE3D
          CALL initialize_coupling (ng, tile, 0)
#endif
        END DO
!$OMP BARRIER
      END DO

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute basic state initial level thicknesses used for state norm
!  scaling. It uses zero free-surface (rest state).  Therefore, the
!  norm scaling is time invariant.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL set_depth (ng, tile)
        END DO
!$OMP BARRIER
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Compute adjoint state vector.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ad_inner2state (ng, tile, Lini, ad_state(ng)%vector)
        END DO
!$OMP BARRIER
      END DO
!
 10   FORMAT (/,a,i2.2,a,i3.3,a,i3.3/)
 20   FORMAT (/,a,i2.2,a,1p,e15.6,/)
 30   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')')

      RETURN
      END SUBROUTINE propagator
