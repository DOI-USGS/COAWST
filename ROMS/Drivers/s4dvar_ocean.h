      MODULE ocean_control_mod
!
!svn $Id: s4dvar_ocean.h 652 2008-07-24 23:20:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Strong Constraint 4-Dimensional Variational (4DVar)       !
!            Data Assimilation Driver:                                 !
!                                                                      !
!  This driver is used for strong constraint 4DVar where the only      !
!  errors considered are those for the observations. The model is      !
!  assumed to be perfect.  This is the "full" method and only the      !
!  nonlinear and adjoint models are needed.                            !
!                                                                      !
!  The misfit (squared difference) between model and observations      !
!  is defined as:                                                      !
!                                                                      !
!         J  = Jb + Jo                                                 !
!                                                                      !
!  where                                                               !
!                                                                      !
!         Jb = transpose(X - Xb) * B^(-1) * (X -Xb)                    !
!                                                                      !
!         Jo = transpose(X - Xo) * O^(-1) * (X -Xo)                    !
!                                                                      !
!         Xb : background state (first guess)                          !
!         Xo : observations                                            !
!         X  : model at observation or background points               !
!         B  : background error covariance                             !
!         O  : observations error covariance (assigned weight)         !
!                                                                      !
!  The routines in this driver control the initialization,  time-      !
!  stepping, and finalization of  ROMS/TOMS  model following ESMF      !
!  conventions:                                                        !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize

      CONTAINS

      SUBROUTINE ROMS_initialize (first, MyCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars
!
#ifdef AIR_OCEAN 
      USE ocean_coupler_mod, ONLY : initialize_atmos_coupling
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_waves_coupling
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: MyCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

      integer :: ng, thread

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(MyCOMM)) THEN
        OCN_COMM_WORLD=MyCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
#endif
!
!-----------------------------------------------------------------------
!  On first pass, initialize model parameters a variables for all
!  nested/composed grids.  Notice that the logical switch "first"
!  is used to allow multiple calls to this routine during ensemble
!  configurations.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
!
!  Initialize parallel parameters.
!
        CALL initialize_parallel
!
!  Initialize wall clocks.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (' Process Information:',/)
        END IF
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            CALL wclock_on (ng, iNLM, 0)
          END DO
!$OMP END PARALLEL DO
        END DO

#if defined AIR_OCEAN || defined WAVES_OCEAN
!
!  Initialize coupling streams between model(s).
!
        DO ng=1,Ngrids
# ifdef AIR_OCEAN
          CALL initialize_atmos_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
          CALL initialize_waves_coupling (ng, MyRank)
# endif
        END DO
#endif
!
!  Read in model tunable parameters from standard input. Initialize
!  "mod_param", "mod_ncparam" and "mod_scalar" modules.
!
        CALL inp_par (iNLM)
        IF (exit_flag.ne.NoError) RETURN
!
!  Allocate and initialize module variables.
!
        CALL mod_arrays (allocate_vars)
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar
      END IF

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Tend)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS nonlinear and adjoint models.     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
!
      USE back_cov_mod, ONLY : back_cov
      USE downhill_mod, ONLY : downhill
      USE ini_fields_mod, ONLY : ini_fields
!
!  Imported variable declarations
!
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      integer :: AdjRec, IniRec, Lsav
      integer :: i, my_iic, ng, subs, tile, thread
#ifdef AVOID_ADJOINT
      integer :: ntend_sav, ntstr_sav
#endif
      real(r8) :: rate
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids

        Lold(ng)=1
        Lnew(ng)=2
        IterSD=0
        CGstepF=CGstepI
        wrtNLmod(ng)=.TRUE.

        ITER_LOOP : DO Nrun=ERstr,ERend
!
!-----------------------------------------------------------------------
!  Time step nonlinear model:  Compute misfit cost function.
!-----------------------------------------------------------------------
!
!  Initialize nonlinear model with first guess initial conditions.
!
          Ipass=1
          CALL initial (ng)
          IF (exit_flag.ne.NoError) RETURN
!
!  If first pass, define output 4DVAR NetCDF file containing all
!  processed data at observation locations.
!
          IF (Nrun.eq.ERstr) THEN
            LdefMOD(ng)=.TRUE.
            CALL def_mod (ng)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          wrtMisfit(ng)=.FALSE.
!
!  Time-step nonlinear model: Compute misfit cost function.
!
          IF (Master) THEN
            WRITE (stdout,10) 'NL', ntstart(ng), ntend(ng)
          END IF

          time(ng)=time(ng)-dt(ng)

          NL_LOOP1 : DO my_iic=ntstart(ng),ntend(ng)+1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL main3d (ng)
#else
            CALL main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO NL_LOOP1
!
!  Report cost function between nonlinear model and observations.
!
          IF (Master) THEN        
            IF (Nrun.gt.ERstr) THEN
              rate=100.0_r8*ABS(FOURDVAR(ng)%CostFun(0)-                &
     &                          FOURDVAR(ng)%CostFunOld(0))/            &
     &                      FOURDVAR(ng)%CostFunOld(0)
            ELSE
              rate=0.0_r8
            END IF
            WRITE (stdout,20) Nrun, FOURDVAR(ng)%CostFun(0), rate
            DO i=1,NstateVar(ng)
              IF (FOURDVAR(ng)%CostFun(i).gt.0.0_r8) THEN
                IF (Nrun.gt.ERstr) THEN
                  rate=100.0_r8*ABS(FOURDVAR(ng)%CostFun(0)-            &
     &                              FOURDVAR(ng)%CostFunOld(0))/        &
     &                          FOURDVAR(ng)%CostFunOld(0)
                ELSE
                  rate=0.0_r8
                END IF
                IF (i.eq.1) THEN
                  WRITE (stdout,30) Nrun, FOURDVAR(ng)%CostFun(i),      &
     &                              TRIM(Vname(1,idSvar(i))), rate
                ELSE
                  WRITE (stdout,40) Nrun, FOURDVAR(ng)%CostFun(i),      &
     &                              TRIM(Vname(1,idSvar(i))), rate
                END IF
              END IF
            END DO
          END IF
          DO i=0,NstateVar(ng)
            FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
          END DO
!
!-----------------------------------------------------------------------
!  Time step adjoint model backwards.
!-----------------------------------------------------------------------
!
!  Initialize adjoint model from rest.
!
          CALL ad_initial (ng)
          IF (exit_flag.ne.NoError) RETURN

#ifdef AVOID_ADJOINT
!
!  Use model-observation misfit as a good approximation of the gradient.
!
          ntstr_sav=ntstart(ng)
          ntend_sav=ntend(ng)
          ntend(ng)=ntimes(ng)
#endif
!
!  Time-step adjoint model: Compute model state gradient, GRAD(J).
!  Force the adjoint model with the adjoint misfit between nonlinear
!  model and observations.
!
          IF (Master) THEN
            WRITE (stdout,10) 'AD', ntstart(ng), ntend(ng)
          END IF

          time(ng)=time(ng)+dt(ng)

          AD_LOOP : DO my_iic=ntstart(ng),ntend(ng),-1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL ad_main3d (ng)
#else
            CALL ad_main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO AD_LOOP
#ifdef AVOID_ADJOINT
          ntstart(ng)=ntstr_sav
          ntend(ng)=ntend_sav
#endif
!
!-----------------------------------------------------------------------
!  Descent algorithm: first pass.
!-----------------------------------------------------------------------
!
!  For efficiency, reinitialize decent algorithm every NSK iterations to
!  avoid searching in directions which have been seached previously. 
!
          IF (MOD(Nrun-IterSD,NiterSD).eq.0) CGstepF=CGstepI
          IF (Master) THEN
            WRITE (stdout,50) Nrun, Ipass, CGstepF
          END IF
!
!  Read in previous state initial conditions and adjoint solution and
!  load it into the appropriate state arrays at descent index "Lold".
!  Also read in latest adjoint solution and load it to descent index
!  "Lnew".
!
          CALL get_state (ng, iNLM, 2, INIname(ng), tINIindx(ng),       &
     &                    Lold(ng))
          IF (exit_flag.ne.NoError) RETURN
          IF (Nrun.gt.ERstr) THEN
            IF (LcycleADJ(ng)) THEN
              AdjRec=3-tADJindx(ng)
            ELSE
              AdjRec=tADJindx(ng)-1
            END IF
            CALL get_state (ng, iADM, 3, ADJname(ng), AdjRec,           &
     &                      Lold(ng))
            IF (exit_flag.ne.NoError) RETURN
          END IF
          CALL get_state (ng, iADM, 4, ADJname(ng), tADJindx(ng),       &
     &                    Lnew(ng))
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute new intial conditions using first guess step size, CGstepF.
!
!$OMP PARALLEL DO PRIVATE(ng,CGstepF,thread,subs,tile)                  &
!$OMP&            SHARED(iNLM,Nrun,numthreads,Lnew)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL back_cov (ng, TILE)
              CALL downhill (ng, TILE, iNLM, Nrun, CGstepF)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  Process new initial conditions: apply boundary conditions, exchange
!  boundary information, and compute (ubar,vbar) from the vertical
!  integral of (u,v), if applicable.
!
          kstp(ng)=Lnew(ng)
# ifdef SOLVE3D
          nstp(ng)=Lnew(ng)
# endif
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(iNLM,,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL ini_fields (ng, TILE, iNLM)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  Write out first guess initial conditions.  This step is needed
!  because of the IO flow. The strategy here is read data instead of
!  increasing memory with additional state arrays.
!
          CALL wrt_ini (ng, Lnew(ng))
          IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Time-step nonlinear model: Compute change in misfit cost function.
!-----------------------------------------------------------------------
!
!  Initialize nonlinear model with estimated initial conditions.
!
          Ipass=2
          CALL initial (ng)
          IF (exit_flag.ne.NoError) RETURN
!
!  Activate switch to write out initial and final misfit between
!  model and observations.
!
          IF ((Nrun.eq.ERstr).or.(Nrun.eq.ERend)) THEN
            wrtMisfit(ng)=.TRUE.
          END IF
!
!  Time-step nonlinear model: Compute change in misfit cost function.
!
          IF (Master) THEN
            WRITE (stdout,10) 'NL', ntstart(ng), ntend(ng)
          END IF

          time(ng)=time(ng)-dt(ng)

          NL_LOOP2 : DO my_iic=ntstart(ng),ntend(ng)+1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL main3d (ng)
#else
            CALL main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO NL_LOOP2
!
!-----------------------------------------------------------------------
!  Descent algorithm: second pass.
!-----------------------------------------------------------------------
!
!  Use Derber method (1985) to compute the refined optimum step size
!  for the conjugate gradient algorithm.
!
          CGstepR=-CGstepF*((StepTopBck+StepTopObs)/                    &
     &                      (StepBotBck+StepBotObs))
          IF (Master) THEN
            WRITE (stdout,50) Nrun, Ipass, CGstepR
          END IF
!
!  Read in previous initial conditions and load it into the appropriate
!  state arrays at descent index "Lold". Reset initial condition time
!  index to write since initial conditions will be over-written with
!  new adjusted values.
!
          IF (LcycleINI(ng)) THEN
            IniRec=3-tINIindx(ng)
          ELSE
            IniRec=tINIindx(ng)-1
          END IF
!!        tINIindx(ng)=IniRec
          CALL get_state (ng, iNLM, 2, INIname(ng), IniRec, Lold(ng))
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute new initial conditions using refined step size, CGstepR.
!
!$OMP PARALLEL DO PRIVATE(ng,CGstepF,thread,subs,tile)                  &
!$OMP&            SHARED(iNLM,Nrun,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL back_cov (ng, TILE)
              CALL downhill (ng, TILE, iNLM, Nrun, CGstepR)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  Process new initial conditions: apply boundary conditions, exchange
!  boundary information, and compute (ubar,vbar) from the vertical
!  integral of (u,v), if applicable.
!
          kstp(ng)=Lnew(ng)
# ifdef SOLVE3D
          nstp(ng)=Lnew(ng)
# endif
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(iNLM,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1,+1
              CALL ini_fields (ng, TILE, iNLM)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!  Write out new adjusted initial conditions.  Notice that the first
!  guess initial conditions are over written.
!
          CALL wrt_ini (ng, Lnew(ng))
          IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Update iteration counters and conjugate gradient step size.
!-----------------------------------------------------------------------
!
          Lsav=Lnew(ng)
          Lnew(ng)=Lold(ng)
          Lold(ng)=Lsav
          CGstepF=CGstepR

        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Done with data assimilation. Initialize the nonlinear model with
!  estimated initial conditions. Save nonlinear solution at observation
!  points for posterior analysis.
!-----------------------------------------------------------------------
!
!  Initialize nonlinear model with estimated initial conditions.
!
        wrtNLmod(ng)=.TRUE.
        wrtMisfit(ng)=.FALSE.
        tRSTindx(ng)=0
        NrecRST(ng)=0
        CALL initial (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!  Run nonlinear model. Interpolate nonlinear model to observation
!  locations.
!
        IF (Master) THEN
          WRITE (stdout,10) 'NL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        NL_LOOP3 : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL main3d (ng)
#else
          CALL main2d (ng)
#endif
          IF (exit_flag.ne.NoError) RETURN

        END DO NL_LOOP3
!
!  Compute and report model-observation comparison statistics.
!
        CALL stats_modobs (ng)

      END DO NEST_LOOP
!
 10   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 20   FORMAT (/,' Iteration = ',i5.5,1x,'Cost Function = ',1p,e17.10,   &
     &        t62,0p,f15.10,' %')
 30   FORMAT (' ----------- ',i5.5,1x,'cost function = ',1p,e17.10,     &
     &        1x,a,t62,0p,f15.10,' %')
 40   FORMAT (13x,i5.5,1x,'cost function = ',1pe17.10,1x,a,             &
     &        t62,0p,f15.10,' %')
 50   FORMAT (/, ' <<<< Descent Algorithm, Iteration = ',i5.5,          &
     &        ', Ipass = ',i1,' >>>>',/,25x,'Step Size = ',1p,e15.8,/)

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear and adjoint models      !
!  execution.                                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: ng, thread
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
!
      DO ng=1,Ngrids
        IF (LwrtRST(ng).and.(exit_flag.eq.1)) THEN
          IF (Master) WRITE (stdout,10)
 10       FORMAT (/,' Blowing-up: Saving latest model state into ',     &
     &              ' RESTART file',/)
          IF (LcycleRST(ng).and.(NrecRST(ng).ge.2)) THEN
            tRSTindx(ng)=2
            LcycleRST(ng)=.FALSE.
          END IF
          blowup=exit_flag
          exit_flag=NoError
          CALL wrt_rst (ng)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks.  Close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,' Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL DO
      END DO
!
!  Close IO files.
!
      CALL close_io

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
