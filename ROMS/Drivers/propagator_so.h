      SUBROUTINE propagator (RunInterval, Iter, state, ad_state)
!
!svn $Id: propagator_so.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Stochastic Optimals Propagator for white noise forcing.             !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_netcdf
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
      USE mod_forces, ONLY : initialize_forces
#ifdef STOCH_OPT_WHITE
      USE packing_mod, ONLY : ad_so_pack, tl_unpack
#else
      USE packing_mod, ONLY : ad_so_pack_red, tl_unpack
#endif
#ifdef SOLVE3D
      USE set_depth_mod, ONLY: set_depth
#endif
!
!  Imported variable declarations.
!
      integer :: Iter

      real(r8), intent(in) :: RunInterval

      TYPE (T_GST), intent(in) :: state(Ngrids)
      TYPE (T_GST), intent(inout) :: ad_state(Ngrids)
!
!  Local variable declarations.
!
      integer :: ng, tile
      integer :: ktmp, ntmp
      integer :: kout, nout
      integer :: Fcount, Interval, IntTrap

      real(r8) :: StateNorm(Ngrids)
      real(r8) :: so_run_time

      logical :: SOrunTL
#ifdef STOCH_OPT_WHITE
      logical :: SOrunAD
#endif
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
!  Loop over the required numger if trapezoidal intervals in time.
!
      INTERVAL_LOOP : DO Interval=1,Nintervals+1

        SOrunTL=.TRUE.
#ifdef STOCH_OPT_WHITE
        SOrunAD=.TRUE.
#endif
        IF (Interval.eq.Nintervals+1) THEN
          SOrunTL=.FALSE.
#ifdef STOCH_OPT_WHITE
          SOrunAD=.FALSE.
#endif
        END IF
        IntTrap=Interval
!
        IF (Master) THEN
          WRITE (stdout,20) ' Stochastic Optimals Time Interval = ',    &
     &                      Interval
        END IF
!
!  Initialize time stepping indices and counters.
!
        DO ng=1,Ngrids
#ifndef STOCH_OPT_WHITE
          IF (Interval.eq.1) THEN
            SOinitial(ng)=.TRUE.
          ELSE
            SOinitial(ng)=.FALSE.
          END IF
#endif
          IF (Interval.eq.1) THEN
            LwrtTLM(ng)=.TRUE.
          ELSE
            LwrtTLM(ng)=.FALSE.
          END IF
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
          tdays(ng)=dstart+REAL(ntimes(ng),r8)*REAL(Interval-1,r8)*     &
     &                     dt(ng)*sec2day/REAL(Nintervals,r8)
          time(ng)=tdays(ng)*day2sec
!$OMP MASTER
          ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
          ntend(ng)=ntimes(ng)
          ntfirst(ng)=ntstart(ng)
          so_run_time=dt(ng)*REAL(ntend(ng)-ntstart(ng)+1,r8)
!$OMP END MASTER
        END DO
!$OMP BARRIER
!
!  Set switches and counters to manage output adjoint and tangent linear
!  history NetCDF files.
!
        DO ng=1,Ngrids
          IF (Iter.gt.0) THEN                ! Arnoldi iteration loop
            IF ((Iter.eq.1).and.(Interval.eq.1)) THEN
              LdefADJ(ng)=.TRUE.
              LdefTLM(ng)=.TRUE.             ! NetCDF files are created
            ELSE                             ! on the first pass
              LdefADJ(ng)=.FALSE.
              LdefTLM(ng)=.FALSE.
            END IF
            Fcount=ADM(ng)%Fcount
            ADM(ng)%Nrec(Fcount)=0
            ADM(ng)%Rindex=0
            Fcount=TLM(ng)%Fcount
            TLM(ng)%Nrec(Fcount)=0
            TLM(ng)%Rindex=0
          ELSE                               ! Computing eigenvectors
            IF (LmultiGST.and.(Interval.eq.1)) THEN
              LdefTLM(ng)=.TRUE.
            ELSE
              LdefTLM(ng)=.FALSE.
            END IF
#ifdef STOCH_OPT_WHITE
            IF (Interval.le.Nintervals) THEN
              Fcount=ADM(ng)%Fcount
              ADM(ng)%Nrec(Fcount)=0
              ADM(ng)%Rindex=0
            END IF
#else
            Fcount=ADM(ng)%Fcount
            ADM(ng)%Nrec(Fcount)=0
            ADM(ng)%Rindex=0
#endif
            IF ((LmultiGST.or.(ABS(Iter).eq.1)).and.                    &
     &          (Interval.eq.1)) THEN
              Fcount=TLM(ng)%Fcount
              TLM(ng)%Nrec(Fcount)=0
              TLM(ng)%Rindex=0
            END IF
          END IF
        END DO
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
!  Unpack tangent linear initial conditions from state vector.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tl_unpack (ng, tile, Nstr(ng), Nend(ng),               &
     &                      state(ng)%vector)
          END DO
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  Compute initial tangent linear state dot product norm.
!-----------------------------------------------------------------------
!
        IF (Interval.eq.Nintervals+1) THEN
          DO ng=1,Ngrids
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL tl_statenorm (ng, tile, kstp(ng), nstp(ng),          &
     &                           StateNorm(ng))
            END DO
!$OMP BARRIER

!$OMP MASTER
            IF (Master) THEN
              WRITE (stdout,30) ' PROPAGATOR - Grid: ', ng,             &
     &                          ',  Tangent Initial Norm: ',            &
     &                          StateNorm(ng)
            END IF
!$OMP END MASTER
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
        IF (SOrunTL) THEN              ! do not run TLM on last interval
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
              WRITE (stdout,40) 'TL', ng, ntstart(ng), ntend(ng)
            END IF
            time(ng)=time(ng)-dt(ng)
!$OMP END MASTER
            iic(ng)=ntstart(ng)-1
          END DO
!$OMP BARRIER

#ifdef SOLVE3D
          CALL tl_main3d (so_run_time)
#else
          CALL tl_main2d (so_run_time)
#endif
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
        END IF
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
!  Compute final tangent linear state dot product norm.
!-----------------------------------------------------------------------
!
        IF (Interval.eq.Nintervals+1) THEN
          DO ng=1,Ngrids
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL tl_statenorm (ng, tile, kstp(ng), nstp(ng),          &
     &                           StateNorm(ng))
            END DO
!$OMP BARRIER

!$OMP MASTER
            IF (Master) THEN
              WRITE (stdout,30) ' PROPAGATOR - Grid: ', ng,             &
     &                          ',  Tangent   Final Norm: ',            &
     &                          StateNorm(ng)
            END IF
!$OMP END MASTER
          END DO
        END IF
!
!=======================================================================
!  Backward integration with the adjoint model.
!=======================================================================
!
!  Initialize time stepping indices and counters.
!
        DO ng=1,Ngrids
#ifndef STOCH_OPT_WHITE
          LwrtState2d(ng)=.FALSE.
          LwrtState3d(ng)=.FALSE.
#endif
          iif(ng)=1
          indx1(ng)=1
          ktmp=knew(ng)
          IF (IntTrap.eq.Nintervals+1) THEN
            ktmp=kstp(ng)
          END IF
          kstp(ng)=1
          krhs(ng)=3
          knew(ng)=2
          kout=knew(ng)
#ifdef STOCH_OPT_WHITE
          IF (IntTrap.eq.Nintervals+1) THEN
            kout=kstp(ng)
          END IF
#endif
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
# ifdef STOCH_OPT_WHITE
          ntend(ng)=1+(Interval-1)*ntimes(ng)/Nintervals
# else
          ntend(ng)=1
# endif
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
            CALL ad_ini_perturb (ng, tile,                              &
     &                           ktmp, kout, ntmp, nstp(ng))
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
#ifdef STOCH_OPT_WHITE
        IF (SOrunAD) THEN              ! do not run ADM on last interval
#endif
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
              WRITE (stdout,40) 'AD', ng, ntstart(ng), ntend(ng)
            END IF
            time(ng)=time(ng)+dt(ng)
!$OMP END MASTER
            iic(ng)=ntstart(ng)+1
          END DO
!$OMP BARRIER

#ifdef SOLVE3D
# ifdef STOCH_OPT_WHITE
          CALL ad_main3d (so_run_time)
# else
          CALL ad_main3d (RunInterval)
# endif
#else
# ifdef STOCH_OPT_WHITE
          CALL ad_main2d (so_run_time)
# else
          CALL ad_main2d (RunInterval)
# endif
#endif
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
#ifdef STOCH_OPT_WHITE
        END IF
#endif
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
!  Pack final adjoint solution into adjoint state vector.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
# ifdef STOCH_OPT_WHITE
            CALL ad_so_pack (ng, tile, Nstr(ng), Nend(ng), IntTrap,     &
     &                       ad_state(ng)%vector)
# else
            CALL ad_so_pack_red (ng, tile, Nstr(ng), Nend(ng),          &
     &                           IntTrap, ad_state(ng)%vector)
# endif
          END DO
!$OMP BARRIER
        END DO
!
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Clear forcing variables for next iteration.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
            CALL initialize_forces (ng, tile, iADM)
          END DO
!$OMP BARRIER
        END DO

      END DO INTERVAL_LOOP
!
 10   FORMAT (/,a,i2.2,a,i3.3,a,i3.3/)
 20   FORMAT (/,a,i2.2)
 30   FORMAT (/,a,i2.2,a,1p,e15.6,/)
 40   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')')

      RETURN
      END SUBROUTINE propagator
