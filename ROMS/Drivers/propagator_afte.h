      SUBROUTINE propagator (RunInterval, state, ad_state)
!
!svn $Id: propagator_afte.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Adjoint Finite Time Eigenvalues Propagator:                         !
!                                                                      !
!  This routine is used during the computation of the eigenvectors of  !
!  the adjoint propagator, transpose[R(t,0)]. They are computed in an  !
!  analogous way to those of  R(t,0).  A  single  integration  of  an  !
!  arbitrary perturbation state vector  "u" backward in time over the  !
!  interval [t,0] by the adjoint model: transpose[R(t,0)]*u.           !
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
      USE dotproduct_mod, ONLY : ad_statenorm
      USE packing_mod, ONLY : ad_unpack, ad_pack
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
#ifdef SOLVE3D
      logical :: FirstPass = .TRUE.
#endif
      integer :: ng, tile

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
        krhs(ng)=3
        knew(ng)=2
        PREDICTOR_2D_STEP(ng)=.FALSE.
!
        iic(ng)=0
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
!  Clear adjoint state variables.  There is not need to clean the basic
!  state arrays since they were zeroth out at initialization and bottom
!  of previous iteration.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iADM)
        END DO
!$OMP BARRIER
      END DO

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute basic state initial level thicknesses used for state norm
!  scaling. It uses zero time averaged free-surface (rest state).
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
!  Unpack adjoint initial conditions from state vector.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ad_unpack (ng, tile, Nstr(ng), Nend(ng),                 &
     &                    state(ng)%vector)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Compute initial adjoint state dot product norm.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL ad_statenorm (ng, tile, knew(ng), nstp(ng),              &
     &                       StateNorm(ng))
        END DO
!$OMP BARRIER

!$OMP MASTER
        IF (Master) THEN
          WRITE (stdout,20) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Adjoint Initial Norm: ', StateNorm(ng)
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
        CALL close_inp (ng, iADM)
        IF (exit_flag.ne.NoError) RETURN
        CALL ad_get_idata (ng)
        IF (exit_flag.ne.NoError) RETURN
        CALL ad_get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
      END DO
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
!  Clear nonlinear state (basic state) variables and insure that the
!  time averaged free-surface is zero for scaling below and next
!  iteration.
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
!  Compute basic state final level thicknesses used for state norm
!  scaling. It uses zero time averaged free-surface (rest state).
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
!  Compute final adjoint state dot product norm.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ad_statenorm (ng, tile, knew(ng), nstp(ng),              &
     &                       StateNorm(ng))
        END DO
!$OMP BARRIER

!$OMP MASTER
        IF (Master) THEN
          WRITE (stdout,20) ' PROPAGATOR - Grid: ', ng,                 &
     &                      ',  Adjoint   Final Norm: ', StateNorm(ng)
        END IF
!$OMP END MASTER
      END DO
!
!-----------------------------------------------------------------------
!  Pack final adjoint solution into adjoint state vector.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL ad_pack (ng, tile, Nstr(ng), Nend(ng),                   &
     &                  ad_state(ng)%vector)
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
