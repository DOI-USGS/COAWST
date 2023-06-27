#include "cppdefs.h"
      SUBROUTINE ntimesteps (model, RunInterval, nl, Nsteps, Rsteps)
!
!git $Id$
!svn $Id: ntimestep.F 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine set the number of time-steps to compute. In nesting    !
!  applications,  the number of time-steps depends on nesting layer    !
!  configuration.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     model        Calling model identifier (integer)                  !
!     RunInterval  Time interval (seconds) to advance forward or       !
!                    backwards the primitive equations (scalar)        !
!     nl           Nesting layer number (integer)                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     nl           Updated nesting layer number (integer)              !
!     Nsteps       Number of time-steps to solve for all grids in      !
!                    current nesting layer "nl" (integer)              !
!     Rsteps       Number of time-steps to complete RunInterval        !
!                    time window (integer)                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in)    :: model
      integer, intent(inout) :: nl
      integer, intent(out)   :: Nsteps, Rsteps
!
      real(dp), intent(in)   :: RunInterval
!
!  Local variable declarations.
!
      integer :: extra, gn, ig, il, ng, ngm1
      integer, dimension(Ngrids) :: WindowSteps, my_Nsteps

#if defined MODEL_COUPLING && !defined MCT_LIB
!
      real(dp) :: ENDtime, NEXTtime
#endif
!
!=======================================================================
!  Set number of time steps to execute for grid "ng".
!=======================================================================
!
!  Initialize.
!
      WindowSteps=0
      my_Nsteps=0
      Nsteps=0
      Rsteps=0
!
!  If appropriate, advance the nested layer counter.
!
      nl=nl+1

#ifdef NESTING
!
!  If refinement and telescoping refined grids, reset the nested layer
!  counter to advance all the grid to the current time of the coarser
!  of all grids (ng=1).
!
     IF (ANY(Telescoping).and.(nl.gt.NestLayers)) THEN
       DO il=1,NestLayers
         DO ig=1,GridsInLayer(il)
           ng=GridNumber(ig,il)
           IF (Telescoping(ng).and.                                     &
     &         (RefineStepsCounter(ng).lt.RefineSteps(ng))) THEN
             nl=il
           END IF
         END DO
       END DO
     END IF
#endif
!
     IF ((0.lt.nl).and.(nl.le.NestLayers)) THEN
!
!  In composite and mosaic grids or all grids in the same nesting
!  layer, it is assumed that the donor and receiver grids have the
!  same time-step size. This is done to avoid the time interpolation
!  between donor and receiver grids. Only spatial interpolations
!  are possible in the current nesting design.
!
!  In grid refinement, it is assumed that the donor and receiver grids
!  are an interger factor of the grid size and time-step size.
!
        WindowSteps=0
!
!  Loop over all grids in current layer nesting layer.
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)

#if defined MODEL_COUPLING && !defined MCT_LIB
!
!  Set extra step parameter needed to finish the simulation due to ROMS
!  delayed output until the next half step. If RunInterval (seconds) is
!  less than full simulation interval because of model coupling, extra=1
!  for last coupling interval and extra=0 otherwise.
!
          NEXTtime=time(ng)+RunInterval
          ENDtime=INItime(ng)+(ntimes(ng)-1)*dt(ng)
          IF (NEXTtime.eq.ENDtime) THEN
            extra=1
          ELSE
            extra=0
          END IF

#elif defined JEDI
!
! OOPS will advance ROMS kernels by smaller intervals, usually a single
! timestep. The strategy here is different from that used for coupling
! since OOPS controls every NLM, TLM, and ADM timestep. The exchange of
! ROMS-to-JEDI state fields is also done one at a time and, in a delay
! mode, to account for the modification to initial conditions and an
! extra half step to complete the solution.
!
          IF (model.eq.iNLM) THEN
            NEXTtime=time(ng)+RunInterval
            ENDtime=INItime(ng)+ntimes(ng)*dt(ng)
            extra=0
          ELSE IF (model.eq.iTLM) THEN
            NEXTtime=time(ng)+RunInterval
            ENDtime=INItime(ng)+ntimes(ng)*dt(ng)
            IF (iic(ng).eq.ntend(ng)) THEN
              extra=1
            ELSE
              extra=0
            END IF
          ELSE IF (model.eq.iADM) THEN
            NEXTtime=time(ng)-RunInterval
            ENDtime=INItime(ng)+ntimes(ng)*dt(ng)
            IF (iic(ng).eq.ntstart(ng)) THEN
              extra=1
            ELSE
              extra=0
            END IF
          END IF
#else
!
!  Here, extra=1, indicates that the RunInterval is the same as
!  simulation interval.
!
          extra=1
#endif
!
!  Determine number of steps in time-interval window.
!
          WindowSteps(ig)=INT((RunInterval+0.5_r8*dt(ng))/dt(ng))
!
!  Advancing model forward: Nonlinear, tangent linear, and representer
!  models.
!
          IF ((model.eq.iNLM).or.                                       &
     &        (model.eq.iTLM).or.                                       &
     &        (model.eq.iRPM)) THEN

#ifdef NESTING
            IF (ANY(CompositeGrid(:,ng))) THEN
              IF (step_counter(ng).le.(WindowSteps(ig)+extra)) THEN
                my_Nsteps(ig)=1
                step_counter(ng)=step_counter(ng)+1
              ELSE
                my_Nsteps(ig)=0
              END IF
            ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).eq.0)) THEN
              RefineStepsCounter=0    ! The coarser grid, reset counters
              IF (step_counter(ng).le.(WindowSteps(ig)+extra)) THEN
                my_Nsteps(ig)=1
                step_counter(ng)=step_counter(ng)+1
                RefineStepsCounter(ng)=RefineStepsCounter(ng)+1
              ELSE
                my_Nsteps(ig)=0
              END IF
            ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
              IF (step_counter(ng).le.(WindowSteps(ig)+extra)) THEN
                IF (Telescoping(ng)) THEN
                  my_Nsteps(ig)=1
                  step_counter(ng)=step_counter(ng)+1
                  RefineStepsCounter(ng)=RefineStepsCounter(ng)+1
                  DO il=nl+1,NestLayers         ! When a parent steps,
                    gn=GridNumber(ig,il)        ! set all its telescoped
                    IF (Telescoping(gn)) THEN   ! children counters to
                      RefineStepsCounter(gn)=0  ! zero
                    END IF
                  END DO
                ELSE
                  my_Nsteps(ig)=RefineSteps(ng)
                  step_counter(ng)=step_counter(ng)+RefineSteps(ng)
                  RefineStepsCounter(ng)=RefineStepsCounter(ng)+        &
     &                                   RefineSteps(ng)
                END IF
              ELSE
                my_Nsteps(ig)=0
              END IF
            END IF
#else
            my_Nsteps(ig)=MAX(my_Nsteps(ig), WindowSteps(ig)+extra)
            step_counter(ng)=step_counter(ng)+WindowSteps(ig)+extra
#endif

#if defined ADJOINT && !defined NESTING
!
!  Advancing model backwards: Adjoint model.
!
          ELSE IF (model.eq.iADM) THEN
            my_Nsteps(ig)=MAX(my_Nsteps(ig), WindowSteps(ig)+extra)
            step_counter(ng)=step_counter(ng)+WindowSteps(ig)+extra
#endif
          END IF
        END DO
!
!  Insure that the steps per time-window are the same.
!
        IF (GridsInLayer(nl).gt.1) THEN
          DO ig=2,GridsInLayer(nl)
            IF (WindowSteps(ig).ne.WindowSteps(ig-1)) THEN
              ngm1=GridNumber(ig-1,nl)
              ng  =GridNumber(ig  ,nl)
              IF (Master) THEN
                WRITE (stdout,10) nl, ngm1, dt(ngm1), ng, dt(ng)
  10            FORMAT (/,' NTIMESTEPS - timestep size are not the ',   &
     &                    ' same in nesting layer: ',i2,                &
     &                  2(/,14x,'Grid ',i2.2,3x,'dt = ',f11.3))
              END IF
              exit_flag=5
            END IF
          END DO
        END IF
!
!  Set number of time-steps to execute. Choose minimum values.
!
        Nsteps=my_Nsteps(1)
        Rsteps=WindowSteps(1)+extra
        DO ig=2,GridsInLayer(nl)
          Nsteps=MIN(Nsteps, my_Nsteps(ig))
          Rsteps=MIN(Rsteps, WindowSteps(ig)+extra)
        END DO
      END IF

      RETURN
      END SUBROUTINE ntimesteps

#if defined NESTING && defined ADJOINT
!
      SUBROUTINE nlm_step_sequence (RunInterval, icount)
!
!git $Id$
!svn $Id: ntimestep.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Given a time interval window (seconds), this routine stores all the !
!  parameters needed to timestep the adjoint model backwards in nested !
!  applications.  It is an elegant way to have the correct parameters  !
!  in complex setups.                                                  !
!                                                                      !
!  StepInfo(:,1)     NL, current nested layer number                   !
!  StepInfo(:,2)     Nsteps, number of timesteps for grids in NL       !
!  StepInfo(:,3)     Rsteps, number of timesteps for simulation window !
!  StepInfo(:,ng+3)  step_counter(ng), timestep counter for each grid  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     RunInterval    Time interval window (seconds) to advance forward !
!                      or backwards the primitive equations (scalar)   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     icount         Number of stored values in "StepInfo" (integer)   !
!                      StepInfo(icount,Ngrid+3)                        !
!                                                                      !
!  It calls repetitively routine "ntimestep".                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_nesting
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(out) :: icount
!
      real(dp), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: DoNestLayer, Time_Step
      integer :: NL, Nsteps, Rsteps
      integer :: dg, ig, istep, ng
!
!=======================================================================
!  Compute the NLM sequence of nested timesteps parameters.
!=======================================================================
!
!  If applicable, allocate and intialize parameter sequence array.  It
!  can only be allocated in "ROMS_run" phase when the RunInterval value
!  is known, and niot before. That value is specified by the Master or
!  Coupler subroutine. Notice that we use the (rs-1) factor to discount
!  the coarser grid step that other grids need to execute to arrive to
!  the same simulation time.
!
      IF (.not.allocated(StepInfo)) THEN
        dg=1
        rs=1
        DO ng=1,Ngrids
          IF (RefineScale(ng).gt.0) THEN
            rs=rs*RefineScale(ng)
          END IF
        END DO
        nStepInfo=INT((RunInterval+0.5_r8*dt(dg))/dt(dg))*(rs-1)+       &
     &            Ngrids
        allocate ( StepInfo(nStepInfo,Ngrids+3) )
        StepInfo=0
      END IF
!
!  Initialize.
!
      Time_Step=.TRUE.
      DoNestLayer=.TRUE.
      icount=0
!
!  Compute and store nonlinear model sequence of NL, Nstep, Rstep, and
!  step_counter parameters. The sequence of values will be reversed
!  when timestepping the adjoint model backwards.
!
      DO WHILE (Time_Step)
        NL=0
        DO WHILE (DoNestLayer)
          CALL ntimesteps (iNLM, RunInterval, NL, Nsteps, Rsteps)
          icount=icount+1
          IF (icount.le.nStepInfo) THEN
            StepInfo(icount,1)=NL
            StepInfo(icount,2)=Nsteps
            StepInfo(icount,3)=Rsteps
            DO ng=1,Ngrids
              StepInfo(icount,ng+3)=step_counter(ng)
            END DO
!
            IF ((NL.le.0).or.(NL.gt.NestLayers)) EXIT
!
            DO istep=1,Nsteps
              DO ig=1,GridsInLayer(NL)
                ng=GridNumber(ig,NL)
                IF (step_counter(ng).eq.Rsteps) Time_Step=.FALSE.
              END DO
            END DO
            IF (.not.Time_Step.and.NL.eq.NestLayers) EXIT
          ELSE
            IF (Master) WRITE(stdout,10) ', nStepInfo = ', nStepInfo,   &
     &                                   icount
            exit_flag=5
            RETURN
          END IF
        END DO
      END DO
!
  10  FORMAT (/,' NLM_STEP_SEQUENCE - too small dimension parameter ',  &
     &        a,i0,3x,i0,/,21x,'Recompute dimensions for ''StepInfo''')
!
      RETURN
      END SUBROUTINE nlm_step_sequence
#endif
