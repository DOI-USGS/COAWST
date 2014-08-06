#include "cppdefs.h"
      SUBROUTINE ntimesteps (nl, ngl, model, RunInterval,               &
     &                       Nsteps, Rsteps)
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
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
!     nl           Nesting layer number (integer)                      !
!     ngl          Number of grids in nested layer (integer)           !
!     model        Calling model identifier (integer)                  !
!     RunInterval  Time interval (seconds) to advance forward or       !
!                    backwards the primitive equations (scalar)        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Nsteps       Number of time-steps to solve for all grids in      !
!                    current nesting layer "nl" (integer)              !
!     Rsteps       Number of time-steps to complete RunInterval        !
!                    time window (integer)                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_iounits
!
!  Imported variable declarations.
!
      integer, intent(in) :: nl, ngl, model
      integer, intent(out) :: Nsteps, Rsteps
!
      real(r8), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      integer :: ig, ng, ngm1
      integer, dimension(ngl) :: WindowSteps, my_Nsteps
!
!=======================================================================
!  Set number of time steps to execute for grid "ng".
!=======================================================================
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
!
!  Determine number of steps in time-interval window.
!
        WindowSteps(ig)=INT((RunInterval+0.5_r8*dt(ng))/dt(ng))
!
!  Advancing model forward: Nonlinear, tangent linear, and representer
!  models.
!
        IF ((model.eq.iNLM).or.                                         &
     &      (model.eq.iTLM).or.                                         &
     &      (model.eq.iRPM)) THEN
          IF (ANY(CompositeGrid(:,ng))) THEN
            IF (step_counter(ng).le.(WindowSteps(ig)+1)) THEN
              my_Nsteps(ig)=1
              step_counter(ng)=step_counter(ng)+1
            ELSE
              my_Nsteps(ig)=0
            END IF
          ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).eq.0)) THEN
            IF (step_counter(ng).le.(WindowSteps(ig)+1)) THEN
              my_Nsteps(ig)=1
              step_counter(ng)=step_counter(ng)+1
            ELSE
              my_Nsteps(ig)=0
            END IF
          ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
            IF (step_counter(ng).le.(WindowSteps(ig)+1)) THEN
!jcw          my_Nsteps(ig)=RefineScale(ng)
              my_Nsteps(ig)=INT(dt(ng-1)/dt(ng))
              step_counter(ng)=step_counter(ng)+RefineScale(ng)
            ELSE
              my_Nsteps(ig)=0
            END IF
          ELSE
            my_Nsteps(ig)=MAX(my_Nsteps(ig), WindowSteps(ig)+1)
            step_counter(ng)=step_counter(ng)+WindowSteps(ig)+1
          END IF

# ifdef ADJOINT
!
!  Advancing model backwards: Adjoint model.
!
        ELSE IF (model.eq.iADM) THEN
          IF (ANY(CompositeGrid(:,ng))) THEN
          ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).eq.0)) THEN
          ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          ELSE
            my_Nsteps(ig)=MAX(my_Nsteps(ig), WindowSteps(ig)+1)
          END IF
# endif
        END IF
      END DO
!
!  Insure that the steps per time-window are the same.
!
      IF (ngl.gt.1) THEN
        DO ig=2,ngl
          IF (WindowSteps(ig).ne.WindowSteps(ig-1)) THEN
            ngm1=GridNumber(ig-1,nl)
            ng  =GridNumber(ig  ,nl)
            IF (Master) THEN
              WRITE (stdout,10) nl, ngm1, dt(ngm1), ng, dt(ng)
  10          FORMAT (/,' NTIMESTEPS - timestep size are not the ',     &
     &                  ' same in nesting layer: ',i2,                  &
     &                2(/,14x,'Grid ',i2.2,3x,'dt = ',f11.3))
            END IF
            exit_flag=5
          END IF
        END DO
      END IF
!
!  Set number of time-steps to execute.
!
      Nsteps=MINVAL(my_Nsteps)
      Rsteps=MINVAL(WindowSteps)+1

      RETURN
      END SUBROUTINE ntimesteps
