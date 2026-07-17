#include "cppdefs.h"
      MODULE esmf_roms_mod

#if defined CICE_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_ice_cice.h 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module sets CICE as the sea-ice model gridded component        !
!  using the generic ESMF/NUOPC layer:                                 !
!                                                                      !
!    ICE_SetServices         Sets CICE component shared-object entry   !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                                                                      !
!    CICE_SetInitializeP1    CICE component phase 1 initialization:    !
!                            sets import and export fields long and    !
!                            short names into its respective state.    !
!                                                                      !
!    CICE_SetInitializeP2    CICE component phase 2 initialization:    !
!                            Initializes component (CICE_Initialize),  !
!                            sets component grid (CICE_SetGridArrays), !
!                            and adds fields into import and export    !
!                            into respective states (CICE_SetStates).  !
!                                                                      !
!    CICE_DataInit           Exports CICE component fields during      !
!                            initialization or restart.                !
!                                                                      !
!    CICE_SetClock           Sets CICE component date calendar, start  !
!                            and stop times, and coupling interval.    !
!                                                                      !
!    CICE_CheckImport        Checks if CICE component import field is  !
!                            at the correct time.                      !
!                                                                      !
!    CICE_SetGridArrays      Sets CICE component horizontal grid       !
!                            arrays, grid area, and land/sea mask.     !
!                                                                      !
!    CICE_SetStates          Adds CICE component export and import     !
!                            fields into its respective state.         !
!                                                                      !
!    CICE_ModelAdvance       Advances CICE component for a coupling    !
!                            interval. It calls import and export      !
!                            routines.                                 !
!                                                                      !
!    CICE_SetFinalize        Finalizes CICE component execution.       !
!                                                                      !
!    CICE_Import             Imports fields into CICE from other       !
!                            gridded components.                       !
!                                                                      !
!    CICE_Export             Exports CICE fields to other gridded      !
!                            components.                               !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!  CICE:   Los Alamos Sea Ice Model                                    !
!            http://oceans11.lanl.gov/trac/CICE                        !
!            https://esgf.esrl.noaa.gov/projects/couplednems/cice_cap  !
!                                                                      !
!=======================================================================
!
      USE ESMF
      USE NUOPC
      USE NUOPC_Model,                                                  &
     &    NUOPC_SetServices          => SetServices,                    &
     &    NUOPC_Label_Advance        => label_Advance,                  &
     &    NUOPC_Label_DataInitialize => label_DataInitialize,           &
     &    NUOPC_Label_SetClock       => label_SetClock
!
      USE mod_esmf_esm          ! ESM coupling structures and variables
!
      USE CICE_InitMod,  ONLY : CICE_Initialize
      USE CICE_RunMod,   ONLY : CICE_Run
      USE CICE_FinalMod, ONLY : CICE_Finalize
!
      implicit none
!
      PUBLIC  :: ICE_SetServices

      PRIVATE :: CICE_SetInitializeP1
      PRIVATE :: CICE_SetInitializeP2
      PRIVATE :: CICE_DataInit
      PRIVATE :: CICE_SetClock
      PRIVATE :: CICE_SetGridArrays
      PRIVATE :: CICE_SetStates
      PRIVATE :: CICE_ModelAdvance
      PRIVATE :: CICE_SetFinalize
      PRIVATE :: CICE_Import
      PRIVATE :: CICE_Export
!
      CONTAINS
!
      SUBROUTINE ICE_SetServices (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets CICE component shared-object entry points for "initialize",    !
!  "run", and "finalize" by using NUOPC generic methods.               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ICE_SetServices"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ICE_SetServices',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS

!-----------------------------------------------------------------------
!  Register NUOPC generic routines.
!-----------------------------------------------------------------------
!
      CALL NUOPC_CompDerive (model,                                     &
     &                       NUOPC_SetServices,                         &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Register initialize routines.
!-----------------------------------------------------------------------
!
!  Set routine for Phase 1 initialization (import and export fields).
!
      CALL NUOPC_CompSetEntryPoint (model,                              &
     &                              methodflag=ESMF_METHOD_INITIALIZE,  &
     &                              phaseLabelList=(/"IPDv00p1"/),      &
     &                              userRoutine=CICE_InitializeP1,      &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for Phase 2 initialization (exchange arrays).
!
      CALL NUOPC_CompSetEntryPoint (model,                              &
     &                              methodflag=ESMF_METHOD_INITIALIZE,  &
     &                              phaseLabelList=(/"IPDv00p2"/),      &
     &                              userRoutine=CICE_InitializeP2,      &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Attach CICE component phase independent specializing methods.
!-----------------------------------------------------------------------
!
!  Set routine for export initial/restart fields.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_DataInitialize,  &
     &                           specRoutine=CICE_DataInit,             &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for setting CICE clock.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_SetClock,        &
     &                           specRoutine=CICE_SetClock,             &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for time-stepping CICE component.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_Advance,         &
     &                           specRoutine=CICE_ModelAdvance,         &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Register CICE finalize routine.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompSetEntryPoint (model,                           &
     &                                 methodflag=ESMF_METHOD_FINALIZE, &
     &                                 userRoutine=CICE_SetFinalize,    &
     &                                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ICE_SetServices',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE ICE_SetServices
!
      SUBROUTINE CICE_SetInitializeP1 (model,                           &
     &                                 ImportState, ExportState,        &
     &                                 clock, rc)
!
!=======================================================================
!                                                                      !
!  CICE component Phase 1 initialization: sets import and export       !
!  fields long and short names into its respective state.              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
      TYPE (ESMF_State)    :: ImportState
      TYPE (ESMF_State)    :: ExportState
      TYPE (ESMF_Clock)    :: clock
!
!  Local variable declarations.
!
      integer :: i, ng
!
      character (len=100) :: CoupledSet, StateLabel
      character (len=240) :: StandardName, ShortName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_SetInitializeP1"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_SetInitializeP1',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set CICE import state and fields.
!-----------------------------------------------------------------------
!
!  Add CICE import state(s). If nesting, each grid has its own import
!  state.
!
      IMPORTING : IF (Nimport(Iseaice).gt.0) THEN
        DO ng=1,MODELS(Iseaice)%Ngrids
          IF (ANY(COUPLED(Iseaice)%LinkedGrid(ng,:))) THEN
            CoupledSet=TRIM(COUPLED(Iseaice)%SetLabel(ng))
            StateLabel=TRIM(COUPLED(Iseaice)%ImpLabel(ng))
            CALL NUOPC_AddNestedState (ImportState,                     &
     &                                 CplSet=TRIM(CoupledSet),         &
     &                                 nestedStateName=TRIM(StateLabel),&
     &                                 nestedState=MODELS(Iseaice)%     &
     &                                                 ImportState(ng), &
                                       rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Add fields import state.
!
            DO i=1,Nimport(Iseaice)
              StandardName=MODELS(Iseaice)%ImportField(i)%standard_name
              ShortName   =MODELS(Iseaice)%ImportField(i)%short_name
              CALL NUOPC_Advertise (MODELS(Iseaice)%ImportState(ng),    &
     &                              StandardName=TRIM(StandardName),    &
     &                              name=TRIM(ShortName),               &
     &                              rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END DO
          END IF
        END DO
      END IF IMPORTING
!
!-----------------------------------------------------------------------
!  Set CICE export state and fields.
!-----------------------------------------------------------------------
!
!  Add CICE import state. If nesting, each grid has its own import
!  state.
!
      EXPORTING : IF (Nexport(Iseaice).gt.0) THEN
        DO ng=1,MODELS(Iseaice)%Ngrids
          IF (ANY(COUPLED(Iseaice)%LinkedGrid(ng,:))) THEN
            CoupledSet=TRIM(COUPLED(Iseaice)%SetLabel(ng))
            StateLabel=TRIM(COUPLED(Iseaice)%ExpLabel(ng))
            CALL NUOPC_AddNestedState (ExportState,                     &
     &                                 CplSet=TRIM(CoupledSet),         &
     &                                 nestedStateName=TRIM(StateLabel),&
     &                                 nestedState=MODELS(Iseaice)%     &
     &                                                 ExportState(ng), &
                                       rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Add fields to export state.
!
            DO i=1,Nexport(Iseaice)
              StandardName=MODELS(Iseaice)%ExportField(i)%standard_name
              ShortName   =MODELS(Iseaice)%ExportField(i)%short_name
              CALL NUOPC_Advertise (MODELS(Iseaice)%ExportState(ng),    &
     &                              StandardName=TRIM(StandardName),    &
     &                              name=TRIM(ShortName),               &
     &                              rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END DO
          END IF
        END DO
      END IF EXPORTING
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_SetInitializeP1',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE CICE_SetInitializeP1
!
      SUBROUTINE CICE_SetInitializeP2 (model,                           &
     &                                 ImportState, ExportState,        &
     &                                 clock, rc)
!
!=======================================================================
!                                                                      !
!  CICE component Phase 2 initialization: Initializes CICE, sets       !
!  component grid, and adds import and export fields to respective     !
!  states.                                                             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
      TYPE (ESMF_State)    :: ImportState
      TYPE (ESMF_State)    :: ExportState
      TYPE (ESMF_Clock)    :: clock
!
!  Local variable declarations.
!
      integer :: MyComm, localPET, ng, PETcount
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_SetInitializeP2"
!
      TYPE (ESMF_Time) :: CurrentTime, StartTime
      TYPE (ESMF_VM)   :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_SetInitializeP2',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the MPI
!  communicator handle and current node rank.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       vm=vm,                                     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_VMGet (vm,                                              &
     &                 localPet=localPET,                               &
     &                 petCount=PETcount,                               &
     &                 mpiCommunicator=MyComm,                          &
     &                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Initialize CICE component.
!-----------------------------------------------------------------------
!
      CALL CICE_Initialize (MyComm)
!
!-----------------------------------------------------------------------
!  Set-up grid and load coordinate data.
!-----------------------------------------------------------------------
!
      DO ng=1,MODELS(Iseaice)%Ngrids
        IF (ANY(COUPLED(Iseaice)%LinkedGrid(ng,:))) THEN
          CALL CICE_SetGridArrays (ng, model, rc)
          IF (ESQMF_LogFoundError(rcToCheck=rc,                         &
     &                            msg=ESMF_LOGERR_PASSTHRU,             &
     &                            line=__LINE__,                        &
     &                            file=MyFile)) THEN
            RETURN
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Set-up fields and register to import/export states.
!-----------------------------------------------------------------------
!
      DO ng=1,MODELS(Iseaice)%Ngrids
        IF (ANY(COUPLED(Iseaice)%LinkedGrid(ng,:))) THEN
          CALL CICE_SetStates (ng, model, rc)
          IF (ESQMF_LogFoundError(rcToCheck=rc,                         &
     &                            msg=ESMF_LOGERR_PASSTHRU,             &
     &                            line=__LINE__,                        &
     &                            file=MyFile)) THEN
            RETURN
          END IF
        END IF
      END DO
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_SetInitializeP2',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE CICE_SetInitializeP2
!
      SUBROUTINE CICE_DataInit (model, rc)
!
!=======================================================================
!                                                                      !
!  Exports CICE component fields during initialization or restart.     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: ng
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_DataInit"
!
      TYPE (ESMF_Time)  :: CurrentTime
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_DataInit',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get gridded component clock current time.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (ClockInfo(Iseaice)%Clock,                     &
     &                    currTime=CurrentTime,                         &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Export initialization or restart fields.
!-----------------------------------------------------------------------
!
!  Run CICE component only for one time-step to fill variables.
!
      CALL CICE_Run ()
!
!  Put export fields.
!
      IF (Nexport(Iseaice).gt.0) THEN
        DO ng=1,MODELS(Iseaice)%Ngrids
          IF (ANY(COUPLED(Iseaice)%LinkedGrid(ng,:))) THEN
            CALL CICE_Export (ng, model, rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
        END DO
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_DataInit',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE CICE_DataInit
!
      SUBROUTINE CICE_SetClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets CICE component date calendar, start and stop time, and         !
!  coupling interval.                                                  !
!                                                                      !
!=======================================================================
!
!!    USE coamnl_mod, ONLY : ktaust   ! starting time (hour, min, sec)
!!    USE coamnl_mod, ONLY : ktauf    ! ending   time (hour, min, sec)
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: PETcount, localPET
      integer :: TimeFrac, ig
# ifdef REGRESS_STARTCLOCK
      integer :: RegressStartDate(7)
# endif
!
# ifdef REGRESS_STARTCLOCK
      character (len= 20) :: RegressStartString
# endif
      character (len= 20) :: Calendar
      character (len=160) :: message

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_SetClock"
!
      TYPE (ESMF_CalKind_Flag) :: CalType
      TYPE (ESMF_Time)         :: StartTime
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_SetClock',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the MPI
!  communicator handle and current node rank.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       localPet=localPET,                         &
     &                       petCount=PETcount,                         &
     &                       vm=vm,                                     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Create CICE component clock.
!-----------------------------------------------------------------------
!
      Calendar=TRIM(ClockInfo(Iseaice)%CalendarString)
      IF (TRIM(Calendar).eq.'gregorian') THEN
        CalType=ESMF_CALKIND_GREGORIAN
      ELSE
        CalType=ESMF_CALKIND_GREGORIAN
      END IF
!
      ClockInfo(Iseaice)%Calendar=ESMF_CalendarCreate(CalType,          &
     &                                             name=TRIM(Calendar), &
     &                                                rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set reference time. Use driver configuration values.
!
      CALL ESMF_TimeSet (ClockInfo(Iseaice)%ReferenceTime,              &
     &                   yy=ReferenceDate(1),                           &
     &                   mm=ReferenceDate(2),                           &
     &                   dd=ReferenceDate(3),                           &
     &                   h =ReferenceDate(4),                           &
     &                   m =ReferenceDate(5),                           &
     &                   s =ReferenceDate(6),                           &
     &                   calendar=ClockInfo(Iseaice)%Calendar,          &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF

# ifdef REGRESS_STARTCLOCK
!
!  Use the same as driver. A coupling interval is substracted to the
!  driver clock to properly initialize all the ESM components.
!
      ClockInfo(Iseaice)%StartTime=ClockInfo(Idriver)%StartTime
!
      CALL ESMF_TimeGet (ClockInfo(Iseaice)%StartTime,                  &
     &                   yy=RegressStartDate(1),                        &
     &                   mm=RegressStartDate(2),                        &
     &                   dd=RegressStartDate(3),                        &
     &                   h= RegressStartDate(4),                        &
     &                   m= RegressStartDate(5),                        &
     &                   s= RegressStartDate(6),                        &
     &                   ms=RegressStartDate(7),                        &
     &                   timeString=RegressStartString,                 &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# else
!
!  Set start time. Use driver configuration values.
!
      CALL ESMF_TimeSet (ClockInfo(Iseaice)%StartTime,                  &
                         yy=StartDate(1),                               &
                         mm=StartDate(2),                               &
                         dd=StartDate(3),                               &
                         h =StartDate(4),                               &
                         m =StartDate(5),                               &
                         s =StartDate(6),                               &
                         calendar=ClockInfo(Iseaice)%Calendar,          &
                         rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# endif
!
!  Set stop time. Use driver configuration values.
!
      CALL ESMF_TimeSet (ClockInfo(Iseaice)%StopTime,                   &
     &                   yy=StopDate(1),                                &
     &                   mm=StopDate(2),                                &
     &                   dd=StopDate(3),                                &
     &                   h =StopDate(4),                                &
     &                   m =StopDate(5),                                &
     &                   s =StopDate(6),                                &
     &                   calendar=ClockInfo(Iseaice)%Calendar,          &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get component clock.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       clock=ClockInfo(Iseaice)%Clock,            &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_ClockGet (ClockInfo(Iseaice)%Clock,                     &
     &                    timeStep=ClockInfo(Iseaice)%TimeStep,         &
     &                    currTime=ClockInfo(Iseaice)%CurrentTime,      &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Compare driver time against CICE component time.
!-----------------------------------------------------------------------
!
      IF (ClockInfo(Idriver)%Restarted) THEN
        StartTime=ClockInfo(Idriver)%RestartTime
      ELSE
        StartTime=ClockInfo(Idriver)%StartTime
      END IF
!
      IF (ClockInfo(Iseaice)%StartTime.ne.StartTime) THEN
        CALL ESMF_TimePrint (ClockInfo(Iseaice)%StartTime,              &
     &                       options="string",                          &
     &                       rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        CALL ESMF_TimePrint (StartTime,                                 &
     &                       options="string",                          &
     &                       rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        message='Driver and CICE start times do not match: '//          &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ClockInfo(Iseaice)%StopTime.ne.                               &
     &    ClockInfo(Idriver)%StopTime) THEN
        CALL ESMF_TimePrint (ClockInfo(Iseaice)%StopTime,               &
     &                       options="string",                          &
     &                       rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        CALL ESMF_TimePrint (ClockInfo(Idriver)%StopTime,               &
     &                       options="string",                          &
     &                       rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        message='Driver and CICE stop times do not match: '//           &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ClockInfo(Iseaice)%Calendar.ne.                               &
     &    ClockInfo(Idriver)%Calendar) THEN
        CALL ESMF_CalendarPrint (ClockInfo(Iseaice)%Calendar,           &
     &                           options="calkindflag",                 &
     &                           rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        CALL ESMF_CalendarPrint (ClockInfo(Idriver)%Calendar,           &
     &                           options="calkindflag",                 &
     &                           rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        message='Driver and CICE calendars do not match: '//            &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Modify component clock time step.
!-----------------------------------------------------------------------
!
      TimeFrac=0
      DO ig=1,MODELS(Iseaice)%Ngrids
        TimeFrac=MAX(TimeFrac,                                          &
     &               MAXVAL(MODELS(Iseaice)%TimeFrac(ig,:),             &
     &                      mask=MODELS(:)%IsActive))
      END DO
      IF (TimeFrac.lt.1) THEN              ! needs to be 1 or greater
        rc=ESMF_RC_NOT_SET                 ! cannot be 0
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
      ClockInfo(Iseaice)%TimeStep=ClockInfo(Idriver)%TimeStep/TimeFrac
!
      ClockInfo(Iseaice)%Name='CICE_clock'
      CALL ESMF_ClockSet (ClockInfo(Iseaice)%Clock,                     &
     &                    name=TRIM(ClockInfo(Iseaice)%Name),           &
     &                    refTime  =ClockInfo(Iseaice)%ReferenceTime,   &
     &                    timeStep =ClockInfo(Iseaice)%TimeStep,        &
     &                    startTime=ClockInfo(Iseaice)%StartTime,       &
     &                    stopTime =ClockInfo(Iseaice)%StopTime,        &
                          rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_SetClock',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE CICE_SetClock
!
      SUBROUTINE CICE_SetGridArrays (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Sets CICE component staggered, horizontal grids arrays, grid area,  !
!  and land/sea mask.                                                  !
!                                                                      !
!=======================================================================
!
      USE ice_blocks,       ONLY : nblocks_tot
      USE ice_blocks,       ONLY : block
      USE ice_blocks,       ONLY : get_block, get_block_parameter
      USE ice_constants,    ONLY : rad_to_deg
      USE ice_distribution, ONLY : ice_distributionGetBlockLoc
      USE ice_domain,       ONLY : nblocks, blocks_ice, distrb_info
      USE ice_domain_size,  ONLY : nx_global, ny_global
      USE ice_grid,         ONLY : TLAT, TLON, ULAT, ULON, tarea,       &
     &                             hm, uvm
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: blk, i, ii, ilo, ihi, j, jj, jlo, jhi
      integer :: gtype, ivar, localDE, n
      integer :: locID, peID
      integer :: lbnd(2),ubnd(2)
!
      integer, pointer :: deLabelList(:) => NULL()
      integer, pointer :: deBlockList(:,:,:) => NULL()
      integer, pointer :: i_glob(:) => NULL()
      integer, pointer :: j_glob(:) => NULL()
      integer, pointer :: PETmap(:) => NULL()
!
      integer (i4b), pointer :: ptrM(:,:) => NULL()
!
      real (dp), pointer :: ptrA(:,:) => NULL()
      real (dp), pointer :: ptrX(:,:) => NULL()
      real (dp), pointer :: ptrY(:,:) => NULL()
!
      character (len=40) :: name

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_SetGridArrays"
!
      TYPE (block)           :: my_block
      TYPE (ESMF_DELayout)   :: delayout
      TYPE (ESMF_DistGrid)   :: distGrid
      TYPE (ESMF_StaggerLoc) :: staggerLoc
!
      TYPE (ESMF_DistGridConnection), allocatable :: connectionList(:)
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_SetGridArrays',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set upper and lower bounds for each decomposition element (DE),
!  create layout, and boundary conditions.
!-----------------------------------------------------------------------
!
      allocate ( deBlockList(2,2,nblocks_tot) )
      allocate ( deLabelList(nblocks_tot) )
      allocate ( PETmap(nblocks_tot) )
!
      DO n=1,nblocks_tot
        deLabelList(i)=n
        CALL get_block_parameter (n, ilo=ilo, ihi=ihi,                  &
     &                               jlo=jlo, jhi=jhi,                  &
     &                               i_glob=i_glob, j_glob=j_glob)
        deBlockList(1,1,n)=i_glob(ilo)
        deBlockList(1,2,n)=i_glob(ihi)
        deBlockList(2,1,n)=j_glob(jlo)
        deBlockList(2,2,n)=j_glob(jhi)
        CALL ice_distributionGetBlockLoc (distrb_info, n, peID, locID)
        PETmap(n)=peID-1
      END DO
!
!  Create decomposition elements layout.
!
      DElayout=ESMF_DELayoutCreate(PETmap,                              &
     &                             rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Connection between tiles: bipolar boundary condition at top row
!                            (nyg).
!
      allocate (connectionList(2))
!
      CALL ESMF_DistGridConnectionSet (connectionList(1),
     &                                 tileIndexA=1,                    &
     &                                 tileIndexB=1,                    &
     &                                 positionVector=(/ nx_global+1,   &
     &                                                 2*ny_global+1/), &
     &                                 orientationVector=(/-1, -2/),    &
     &                                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Connectivity between tiles: periodic boundary condition along first
!                              dimension.
!
      CALL ESMF_DistGridConnectionSet (connectionList(2),               &
     &                                 tileIndexA=1,                    &
     &                                 tileIndexB=1,                    &
     &                                 positionVector=(/nx_global, 0/), &
     &                                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      deallocate (connectionList)
!
!-----------------------------------------------------------------------
!  Create DistGrid based on model domain decomposition.
!-----------------------------------------------------------------------
!
      distgrid=ESMF_DistGridCreate(minIndex=(/ 1, 1 /),                 &
     &                             maxIndex=(/ nx_global, ny_global /), &
     &                             deBlockList=deBlockList,             &
     &                             delayout=delayout,                   &
     &                             connectionList=connectionList,       &
     &                             rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      deallocate (deLabelList)
      deallocate (deBlockList)
      deallocate (PETmap)
!
!-----------------------------------------------------------------------
!  Set component grid coordinates.
!-----------------------------------------------------------------------
!
!  Define component grid location type: Arakawa B-grid.
!
      IF (.not.allocated(MODELS(Iseaice)%mesh)) THEN
        allocate ( MODELS(Iseaice)%mesh(2) )
        MODELS(Iseaice)%mesh(1)%gtype=Icenter            ! T-cell
        MODELS(Iseaice)%mesh(2)%gtype=Icorner            ! UV-cell
      END IF
!
!  Create ESMF Grid.
!
      MODELS(Iseaice)%grid(ng)=ESMF_GridCreate(distgrid=distGrid,       &
     &                                  coordSys=ESMF_COORDSYS_SPH_DEG, &
     &                                  gridEdgeLWidth=(/ 0, 0 /),      &
     &                                  gridEdgeUWidth=(/ 0, 1 /),      &
     &                                         name="cice_grid",        &
     &                                         rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get number of local decomposition elements (DEs). Usually, a single
!  DE is associated with each Persistent Execution Thread (PETs). Thus,
!  localDEcount=1.
!
      CALL ESMF_GridGet (MODELS(Iseaice)%grid(ng),                      &
     &                   localDECount=localDEcount,                     &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Mesh coordinates for each variable type.
!
      MESH_LOOP : DO ivar=1,UBOUND(MODELS(Iseaice)%mesh, DIM=1)
!
!  Set staggering type, Arakawa B-grid.
!
        SELECT CASE (MODELS(Iseaice)%mesh(ivar)%gtype)
          CASE (Icorner)
            staggerLoc=ESMF_STAGGERLOC_CORNER
          CASE (Icenter)
            staggerLoc=ESMF_STAGGERLOC_CENTER
!
!  Allocate storage for masking.
!
            CALL ESMF_GridAddItem (MODELS(Iseaice)%grid(ng),            &
     &                             staggerLoc=staggerLoc,               &
     &                             itemflag=ESMF_GRIDITEM_MASK,         &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
            MODELS(Iseaice)%LandValue=0
            MODELS(Iseaice)%SeaValue=1
!
!  Allocate storage for grid area.
!
            CALL ESMF_GridAddItem (MODELS(Iseaice)%grid(ng),            &
     &                             staggerLoc=staggerLoc,               &
     &                             itemflag=ESMF_GRIDITEM_AREA,         &
     &                              rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
        END SELECT
!
!  Allocate coordinate storage associated with staggered grid type.
!  No coordinate values are set yet.
!
        CALL ESMF_GridAddCoord (MODELS(Iseaice)%grid(ng),               &
     &                          staggerLoc=staggerLoc,                  &
     &                          rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get pointers and set coordinates for the grid using the block
!  decomposition.
!
        BLOCK_LOOP : DO blk=1,nblocks
          localDE=blk-1
          my_block=get_block(blocks_ice(blk), blk)
          ilo=my_block%ilo
          ihi=my_block%ihi
          jlo=my_block%jlo
          jhi=my_block%jhi
!
          CALL ESMF_GridGetCoord (MODELS(Iseaice)%grid(ng),             &
     &                            coordDim=1,                           &
     &                            localDE=localDE,                      &
     &                            staggerLoc=staggerLoc,                &
     &                            computationalLBound=lbnd,             &
     &                            computationalUBound=ubnd,             &
     &                            farrayPtr=ptrX,                       &
     &                            rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_GridGetCoord (MODELS(Iseaice)%grid(ng),             &
     &                            coordDim=2,                           &
     &                            localDE=localDE,                      &
     &                            staggerLoc=staggerLoc,                &
     &                            farrayPtr=ptrY,                       &
     &                            rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_GridGetItem (MODELS(Iseaice)%grid(ng),              &
     &                           localDE=localDE,                       &
     &                           staggerLoc=staggerLoc,                 &
     &                           itemflag=ESMF_GRIDITEM_MASK,           &
     &                           farrayPtr=ptrM,                        &
     &                           rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_GridGetItem (MODELS(Iseaice)%grid(ng),              &
     &                           localDE=localDE,                       &
     &                           staggerLoc=staggerLoc,                 &
     &                           itemflag=ESMF_GRIDITEM_AREA,           &
     &                           farrayPtr=ptrA,                        &
     &                           rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
!  Fill grid pointers.
!
          SELECT CASE (MODELS(Iseaice)%mesh(ivar)%gtype)
            CASE (Icorner)
              DO jj=lbnd(2),ubnd(2)
                j=jj+jlo-lbnd(2)
                DO ii=lbnd(1),ubnd(1)
                  i=ii+ilo-lbnd(1)
                  ptrX(ii,jj)=ULON(i-1,j-1,blk)*rad_to_deg
                  ptrY(ii,jj)=ULAT(i-1,j-1,blk)*rad_to_deg
                  ptrM(ii,jj)=NINT(uvm(i,j,blk))
                  ptrA(ii,jj)=tarea(i,j,blk)
                END DO
              END DO
            CASE (Icenter)
              DO jj=lbnd(2),ubnd(2)
                j=jj+jlo-lbnd(2)
                DO ii=lbnd(1),ubnd(1)
                  i=ii+ilo-lbnd(1)
                  ptrX(ii,jj)=TLON(i,j,blk)*rad_to_deg
                  ptrY(ii,jj)=TLAT(i,j,blk)*rad_to_deg
                  ptrM(ii,jj)=NINT(hm(i,j,blk))
                  ptrA(ii,jj)=tarea(i,j,blk)
                END DO
              END DO
          END SELECT
!
!  Nullify pointers.
!
          IF ( associated(ptrX) ) nullify (ptrX)
          IF ( associated(ptrY) ) nullify (ptrY)
          IF ( associated(ptrM) ) nullify (ptrM)
          IF ( associated(ptrA) ) nullify (ptrA)
        END DO BLOCK_LOOP
!
!  Debugging: write out component grid in VTK format.
!
        IF (DebugLevel.ge.4) THEN
          gtype=MODELS(Iseaice)%mesh(ivar)%gtype
          CALL ESMF_GridWriteVTK (MODELS(Iseaice)%grid(ng),             &
     &                            filename="cice_"//                    &
     &                                     TRIM(GridType(gtype))//      &
     &                                     "_point",                    &
     &                            staggerLoc=staggerLoc,                &
     &                            rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
        END IF
      END DO MESH_LOOP
!
!  Assign grid to gridded component.
!
      CALL ESMF_GridCompSet (model,                                     &
     &                       grid=MODELS(Iseaice)%grid(ng),             &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_SetGridArrays',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE CICE_SetGridArrays
!
      SUBROUTINE CICE_SetStates (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Adds CICE component export and import fields into its respective    !
!  state.                                                              !
!                                                                      !
!=======================================================================
!
      USE ice_domain,       ONLY : nblocks
      USE ice_domain_size,  ONLY : max_blocks
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: i, id
      integer :: blk, localDE
      integer :: localPET
      integer :: ExportCount, ImportCount
!
      real (dp), pointer :: ptr3d(:,:,:) => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_SetStates"

      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)
!
      TYPE (ESMF_ArraySpec)  :: arraySpec3d
      TYPE (ESMF_Field)      :: field
      TYPE (ESMF_StaggerLoc) :: staggerLoc
      TYPE (ESMF_VM)         :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_SetStates',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Query gridded component.
!-----------------------------------------------------------------------
!
!  Get import and export states.
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       localPet=localPET,                         &
     &                       vm=vm,                                     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Set a 3D floating-point array descriptor.  CICE import and export
!  fields are dimensioned (nx_global, ny_global, max_blocks).
!-----------------------------------------------------------------------
!
      CALL ESMF_ArraySpecSet (arraySpec3d,                              &
     &                        typekind=ESMF_TYPEKIND_R8,                &
     &                        rank=3,                                   &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Add export fields into export state.
!-----------------------------------------------------------------------
!
      EXPORTING : IF (Nexport(Iseaice).gt.0) THEN
!
!  Get number of fields to export.
!
        CALL ESMF_StateGet (MODELS(Iseaice)%ExportState(ng),            &
     &                      itemCount=ExportCount,                      &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get a list of export fields names.
!
        IF (.not.allocated(ExportNameList)) THEN
          allocate ( ExportNameList(ExportCount) )
        END IF
        CALL ESMF_StateGet (MODELS(Iseaice)%ExportState(ng),            &
     &                      itemNameList=ExportNameList,                &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Set export field(s).
!
        DO i=1,ExportCount
          id=field_index(MODELS(Iseaice)%ExportField, ExportNameList(i))
!
          IF (NUOPC_IsConnected(MODELS(Iseaice)%ExportState(ng),        &
     &                          fieldName=TRIM(ExportNameList(i)),      &
     &                          rc=rc)) THEN
!
!  Set staggering type.
!
            SELECT CASE (MODELS(Iseaice)%ExportField(id)%gtype)
              CASE (Icenter)
                staggerLoc=ESMF_STAGGERLOC_CENTER
              CASE (Icorner)
                staggerLoc=ESMF_STAGGERLOC_CORNER
            END SELECT
!
!  Create 2D field from the Grid and arraySpec.
!
            field=ESMF_FieldCreate(MODELS(Iseaice)%grid(ng),            &
     &                             arraySpec3d,                         &
     &                             indexflag=ESMF_INDEX_DELOCAL,        &
     &                             staggerloc=staggerLoc,               &
     &                             ungriddedLBound=(/1/),               &
     &                             ungriddedUBound=(/max_blocks/),      &
     &                             name=TRIM(ExportNameList(i)),        &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Put data into state. Use CICE block decomposition.
!
            DO blk=1,nblocks
              localDE=blk-1
!
!  Get pointer to DE-local memory allocation within field.
!
              CALL ESMF_FieldGet (field,                                &
     &                            localDe=localDE,                      &
     &                            farrayPtr=ptr3d,                      &
     &                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  Initialize pointer.
!
              ptr3d=MISSING_dp
!
!  Nullify pointer to make sure that it does not point on a random part
!  in the memory.
!
              IF ( associated(ptr3d) ) nullify (ptr3d)
            END DO
!
!  Add field export state.
!
            CALL NUOPC_Realize (ExportState,                            &
     &                          field=field,                            &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Remove field from export state because it is not connected.
!
          ELSE
            IF (localPET.eq.0) THEN
              WRITE (cplout,10) TRIM(ExportNameList(i)),                &
     &                          'Export State: ',                       &
     &                          TRIM(COUPLED(Iseaice)%ExpLabel(ng))
            END IF
            CALL ESMF_StateRemove (MODELS(Iseaice)%ExportState(ng),     &
     &                             (/ TRIM(ExportNameList(i)) /),       &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
        END DO
!
!  Deallocate arrays.
!
        IF ( allocated(ExportNameList) ) deallocate (ExportNameList)
!
      END IF EXPORTING
!
!-----------------------------------------------------------------------
!  Add import fields into import state.
!-----------------------------------------------------------------------
!
      IMPORTING : IF (Nimport(Iseaice).gt.0) THEN
!
!  Get number of fields to import.
!
        CALL ESMF_StateGet (MODELS(Iseaice)%ImportState(ng),            &
     &                      itemCount=ImportCount,                      &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get a list of import fields names.
!
        IF (.not.allocated(ImportNameList)) THEN
          allocate (ImportNameList(ImportCount))
        END IF
        CALL ESMF_StateGet (MODELS(Iseaice)%ImportState(ng),            &
     &                      itemNameList=ImportNameList,                &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Set import field(s).
!
        DO i=1,ImportCount
          id=field_index(MODELS(Iseaice)%ImportField, ImportNameList(i))
!
          IF (NUOPC_IsConnected(MODELS(Iseaice)%ImportState(ng),        &
     &                          fieldName=TRIM(ImportNameList(i)),      &
     &                          rc=rc)) THEN

!
!  Set staggering type.
!
            SELECT CASE (MODELS(Iseaice)%ImportField(id)%gtype)
              CASE (Icenter)
                staggerLoc=ESMF_STAGGERLOC_CENTER
              CASE (Icorner)
                staggerLoc=ESMF_STAGGERLOC_CORNER
            END SELECT
!
!  Create field from the Grid, arraySpec.
!
            field=ESMF_FieldCreate(MODELS(Iseaice)%grid(ng),            &
     &                             arraySpec3d,                         &
     &                             indexflag=ESMF_INDEX_DELOCAL,        &
     &                             staggerloc=staggerLoc,               &
     &                             ungriddedLBound=(/1/),               &
     &                             ungriddedUBound=(/max_blocks/),      &
     &                             name=TRIM(ImportNameList(i)),        &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Put data into state. Use CICE block decomposition.
!
            DO blk=1,nblocks
              localDE=blk-1
!
!  Get pointer to DE-local memory allocation within field.
!
              CALL ESMF_FieldGet (field,                                &
     &                            localDe=localDE,                      &
     &                            farrayPtr=ptr3d,                      &
     &                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  Initialize pointer.
!
              ptr3d=MISSING_dp
!
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!
              IF (associated(ptr3d)) nullify (ptr3d)
            END DO
!
!  Add field import state.
!
            CALL NUOPC_Realize (MODELS(Iseaice)%ImportState(ng),        &
     &                          field=field,                            &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Remove field from import state because it is not connected.
!
          ELSE
            IF (localPET.eq.0) THEN
              WRITE (cplout,10) TRIM(ImportNameList(i)),                &
     &                          'Import State: ',                       &
     &                          TRIM(COUPLED(Iseaice)%ImpLabel(ng))
            END IF
            CALL ESMF_StateRemove (MODELS(Iseaice)%ImportState(ng),     &
     &                             TRIM(ImportNameList(i)),             &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
        END DO
!
!  Deallocate arrays.
!
        IF (allocated(ImportNameList)) deallocate (ImportNameList)
!
      END IF IMPORTING
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_SetStates',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE CICE_SetStates
!
      SUBROUTINE CICE_ModelAdvance (model, rc)
!
!=======================================================================
!                                                                      !
!  Advance CICE component for a coupling interval (seconds) using      !
!  "CICE_run". It also calls "CICE_Import" and "CICE_Export" to import !
!  and export coupling fields, respectively.                           !
!                                                                      !
!=======================================================================
!
!!    USE mod_runparams, ONLY : ifrest, ktau, dtsrf
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      logical :: Ladvance
      integer :: is, ng
      integer :: localPET, phase
!
      real (dp) :: CouplingInterval, RunInterval
      real (dp) :: TcurrentInSeconds, TstopInSeconds
!
      character (len=22) :: Cinterval
      character (len=22) :: CurrTimeString, StopTimeString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_ModelAdvance"
!
      TYPE (ESMF_Clock)        :: clock
      TYPE (ESMF_State)        :: ExportState, ImportState
      TYPE (ESMF_Time)         :: ReferenceTime
      TYPE (ESMF_Time)         :: CurrentTime, StopTime
      TYPE (ESMF_TimeInterval) :: TimeStep
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_ModelAdvance',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get information about the gridded component.
!-----------------------------------------------------------------------
!
!  Inquire about CICE component.
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       importState=ImportState,                   &
     &                       exportState=ExportState,                   &
     &                       clock=clock,                               &
     &                       localPet=localPET,                         &
     &                       currentPhase=phase,                        &
     &                       vm=vm,                                     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get time step interval, stopping time, reference time, and current
!  time.
!
      CALL ESMF_ClockGet (clock,                                        &
     &                    timeStep=TimeStep,                            &
     &                    stopTime=StopTime,                            &
     &                    refTime=ReferenceTime,                        &
     &                    currTime=ClockInfo(Iseaice)%CurrentTime,      &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Current CICE time (seconds).
!
      CALL ESMF_TimeGet (ClockInfo(Iseaice)%CurrentTime,                &
     &                   s_r8=TcurrentInSeconds,                        &
     &                   timeStringISOFrac=CurrTimeString,              &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(CurrTimeString, 'T')                 ! remove 'T' in
      IF (is.gt.0) CurrTimeString(is:is)=' '        ! ISO 8601 format
!
!  CICE stop time (seconds) for this coupling window.
!
      CALL ESMF_TimeGet (ClockInfo(Iseaice)%CurrentTime+TimeStep,       &
     &                   s_r8=TstopInSeconds,                           &
     &                   timeStringISOFrac=StopTimeString,              &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(StopTimeString, 'T')                 ! remove 'T' in
      IF (is.gt.0) StopTimeString(is:is)=' '        ! ISO 8601 form
!
!  Get coupling time interval (seconds, double precision).
!
      CALL ESMF_TimeIntervalGet (TimeStep,                              &
     &                           s_r8=CouplingInterval,                 &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set CICE running interval (seconds) for the current coupling window.
!
      Ladvance=.TRUE.
      RunInterval=CouplingInterval
!
!-----------------------------------------------------------------------
!  Report time information strings (YYYY-MM-DD hh:mm:ss).
!-----------------------------------------------------------------------
!
      IF (localPET.eq.0) THEN
        WRITE (Cinterval,'(f15.2)') CouplingInterval
        WRITE (cplout,10) TRIM(CurrTimeString), TRIM(StopTimeString),   &
     &                    phase, TRIM(ADJUSTL(Cinterval))
      END IF
!
!-----------------------------------------------------------------------
!  Get import fields from other ESM components.
!-----------------------------------------------------------------------
!
      IF ((Nimport(Iseaice).gt.0).and.                                  &
     &    (TcurrentInSeconds.gt.ClockInfo(Idriver)%Time_Start)) THEN
        DO ng=1,MODELS(Iseaice)%Ngrids
          IF (ANY(COUPLED(Iseaice)%LinkedGrid(ng,:))) THEN
            CALL CICE_Import (ng, model, rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
        END DO
      ELSE
        Ladvance=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  Run CICE component.
!-----------------------------------------------------------------------
!
      IF (Ladvance)) THEN
        CALL CICE_Run
      END IF
!
!-----------------------------------------------------------------------
!  Put export fields.
!-----------------------------------------------------------------------
!
      IF (Nexport(Iseaice).gt.0) THEN
        DO ng=1,MODELS(Iseaice)%Ngrids
          IF (ANY(COUPLED(Iseaice)%LinkedGrid(ng,:))) THEN
            CALL CICE_Export (ng, model, rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
        END DO
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_ModelAdvance',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (3x,'ModelAdvance - ESMF, Running CICE:',t42,a,            &
     &        ' => ',a,', Phase: ',i1,' [',a,' s]')
!
      RETURN
      END SUBROUTINE CICE_ModelAdvance
!
      SUBROUTINE CICE_SetFinalize (model,                               &
     &                             ImportState, ExportState,            &
     &                             clock, rc)
!
!=======================================================================
!                                                                      !
!  Finalize CICE component execution. It calls CICE_finalize.          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_Clock)    :: clock
      TYPE (ESMF_GridComp) :: model
      TYPE (ESMF_State)    :: ExportState
      TYPE (ESMF_State)    :: ImportState
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_SetFinalize"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_SetFinalize',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Finalize CICE component.
!-----------------------------------------------------------------------
!
      CALL CICE_Finalize
      CALL my_flush (6)                   ! flush standard output buffer
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_SetFinalize',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE CICE_SetFinalize
!
      SUBROUTINE CICE_Import (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Imports fields into CICE array structure from other coupled         !
!  gridded components.                                                 !
!                                                                      !
!=======================================================================
!
      USE ice_blocks,       ONLY : block
      USE ice_blocks,       ONLY : get_block
      USE ice_domain,       ONLY : nblocks, blocks_ice
      USE ice_domain_size,  ONLY : max_blocks, nx_global, ny_global
      USE ice_grid,         ONLY : ANGLET
      USE ice_flux
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      logical :: got_pair, got_tair
      logical :: got_current(2), got_swfx(4), got_wind(2), got_wstr(2)
!
      integer :: id, ifld
      integer :: blk, i, ii, j, jj
      integer :: iyear, iday, imonth, ihour
      integer :: ImportCount
      integer :: localPET
      integer :: year, month, day, hour, minutes, seconds, sN, SD
!
      real (dp) :: ciceScale, scale, add_offset
      real (dp) :: TimeInDays, Time_Current, Tmin, Tmax, Tstr, Tend
      real (dp) :: sigma_c, sigma_l, sigma_r, slopex, slopey
!
      real (dp), dimension(nx_global,ny_global,max_blocks) :: Pair
!
      real (dp), pointer :: ptr3d(:,:,:) => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_Import"

      character (ESMF_MAXSTR) :: ofile
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)
!
      TYPE (block)      :: my_block
      TYPE (ESMF_Clock) :: clock
      TYPE (ESMF_Field) :: field
      TYPE (ESMF_Time)  :: CurrentTime
      TYPE (ESMF_VM)    :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_Import',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get information about the gridded component.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       clock=clock,                               &
     &                       localPet=localPET,                         &
     &                       vm=vm,                                     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get current time.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (clock,                                        &
     &                    currTime=CurrentTime,                         &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeGet (CurrentTime,                                   &
     &                   yy=year,                                       &
     &                   mm=month,                                      &
     &                   dd=day,                                        &
     &                   h =hour,                                       &
     &                   m =minutes,                                    &
     &                   s =seconds,                                    &
     &                   sN=sN,                                         &
     &                   sD=sD,                                         &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeGet (CurrentTime,                                   &
     &                   s_r8=Time_Current,                             &
     &                   timeString=Time_CurrentString,                 &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      TimeInDays=Time_Current/86400.0_dp
      is=INDEX(Time_CurrentString, 'T')              ! remove 'T' in
      IF (is.gt.0) Time_CurrentString(is:is)=' '     ! ISO 8601 format
!
!-----------------------------------------------------------------------
!  Get list of import fields.
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet (MODELS(Iseaice)%ImportState(ng),              &
     &                    itemCount=ImportCount,                        &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (.not.allocated(ImportNameList)) THEN
        allocate ( ImportNameList(ImportCount) )
      END IF
      CALL ESMF_StateGet (MODELS(Iseaice)%ImportState(ng),              &
     &                    itemNameList=ImportNameList,                  &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get import fields.
!-----------------------------------------------------------------------
!
      got_pair=.FALSE.
      got_tair=.FALSE.
      got_current(1:2)=.FALSE.
      got_swfx(1:4)=.FALSE.
      got_wind(1:2)=.FALSE.
      got_wstr(1:2)=.FALSE.
!
      FLD_LOOP : DO ifld=1,ImportCount
        id=field_index(MODELS(Iseaice)%ImportField,ImportNameList(ifld))
!
!  Get field from import state.
!
        CALL ESMF_StateGet (MODELS(Iseaice)%ImportState(ng),            &
     &                      TRIM(ImportNameList(ifld)),                 &
     &                      field,                                      &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get field pointer.
!
        CALL ESMF_FieldGet (field,                                      &
     &                      farrayPtr=ptr3d,                            &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Load import data into CICE component variable.
!
        scale=MODELS(Iseaice)%ImportField(id)%scale_factor
        add_offset=MODELS(Iseaice)%ImportField(id)%add_offset
!
          MyFmin(1)= MISSING_dp
          MyFmax(1)=-MISSING_dp
          MyFmin(2)= MISSING_dp
          MyFmax(2)=-MISSING_dp
!
        SELECT CASE (TRIM(ADJUSTL(itemNameList(ifld))))
!
!  Atmospheric height of the lowest level (m), from ATM.
!
          CASE ('zlvl', 'inst_height_lowest')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  zlvl(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Air density (kg m-3) at surface defined by inst_height_lowest (near
!  surface; maybe lowest level), from ATM.
!
          CASE ('rhoa', 'air_density_height_lowest')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  rhoa(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Air pressure (N m-2) at surface defined by inst_height_lowest (near
!  surface; maybe lowest level), from ATM.
!
          CASE ('Pair', 'ips', 'inst_pres_height_lowest')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  Pair(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_pair=.TRUE.
!
!  Air temperature (K) at surface defined by inst_height_lowest (near
!  surface; maybe lowest level), from ATM.
!
          CASE ('Tair', 'its', 'inst_temp_height_lowest')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  Tair(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_tair=.TRUE.
!
!  Air humidity (kg kg-1), at surface defined by inst_height_lowest
! (near surface; maybe lowest level), from ATM.
!
          CASE ('Qair', 'Qa', 'ishh', 'inst_spec_humid_height_lowest')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  Qa(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Downwelling longwave flux (W m-2), averaged over coupling interval,
!  from ATM.
!
          CASE ('flw', 'mdlwfx', 'mean_down_lw_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  flw(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Visible direct band of downward shortwave flux (W m-2), averaged
!  over the coupling interval, from ATM.
!
          CASE ('swvdr', 'sw_flux_vis_dir', 'mean_down_sw_vis_dir_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  swvdr(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_swfx(1)=.TRUE.
!
!  Visible diffusive band of downward shortwave flux (W m-2), averaged
!  over the coupling interval, from ATM.
!
          CASE ('swvdf', 'sw_flux_vis_dif', 'mean_down_sw_vis_dif_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  swvdf(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_swfx(2)=.TRUE.
!
!  Infrared direct band of downward shortwave flux (W m-2), averaged
!  over the coupling interval, from ATM.
!
          CASE ('swidr', 'sw_flux_nir_dir', 'mean_down_sw_ir_dir_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  swidr(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_swfx(3)=.TRUE.
!
!  Infrared diffusive band of downward shortwave flux (W m-2), averaged
!  over the coupling interval, from ATM.
!
          CASE ('swidf', 'sw_flux_nir_dif', 'mean_down_sw_ir_dif_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  swidf(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_swfx(4)=.TRUE.
!
!  Near surface (maybe lowest level) U-wind component (m s-1), from ATM.
!  Needs to be rotated from east/north to i,j coordinates after it is
!  loaded.
!
          CASE ('Uwind', 'uatm', 'inst_zonal_wind_height_lowest')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  uatm(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_wind(1)=.TRUE.
!
!  Near surface (maybe lowest level) V-wind component (m s-1), from ATM.
!  Needs to be rotated from east/north to i,j coordinates after it is
!  loaded.
!
          CASE ('Vwind', 'vatm', 'inst_merid_wind_height_lowest')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  vatm(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_wind(2)=.TRUE.
!
!  Near surface (maybe lowest level) U-wind stress component (N m-2),
!  averaged over the coupling interval, from ATM.  Needs to be rotated
!  from east/north to i,j coordinates after it is loaded.
!
          CASE ('Ustr', 'strax', 'mzmfx', 'mean_zonal_moment_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  strax(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_wstr(1)=.TRUE.
!
!  Near surface (maybe lowest level) V-wind stress component (N m-2)),
!  averaged over the coupling interval, from ATM.  Needs to be rotated
!  from east/north to i,j coordinates after it is loaded.
!
          CASE ('Vstr', 'stray', 'mmmfx', 'mean_merid_momentum_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  stray(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_wstr(2)=.TRUE.
!
!  Liquid precipitation rate (kg m-2 s-1), averaged over the coupling
!  interval, from ATM.
!
          CASE ('frain', 'lprec', 'mean_prec_rate')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  frain(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Frozen/snow precipitation rate (kg m-2 s-1), averaged over the coupling
!  interval, from ATM.
!
          CASE ('fsnow', 'fprec', 'mean_fprec_rate')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  fsnow(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Sea surface heigh (m), from OCN.
!
          CASE ('ssh', 'sea_lev')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
!                                       zonal sea surface slope
!
                  sigma_r=0.5_dp*(ptr3d(ii+1,jj+1,blk)-                 &
     &                            ptr3d(ii  ,jj+1,blk)+                 &
     &                            ptr3d(ii+1,jj  ,blk)-                 &
     &                            ptr3d(ii  ,jj  ,blk))/dxt(i,j,blk)
                  sigma_l=0.5_dp*(ptr3d(ii  ,jj+1,blk)-                 &
     &                            ptr3d(ii-1,jj+1,blk)+                 &
     &                            ptr3d(ii  ,jj  ,blk)-                 &
     &                            ptr3d(ii-1,jj  ,blk))/dxt(i,j,blk)
                  sigma_c=0.5_dp*(sigma_r+sigma_l)
                  IF ((sigma_r*sigma_l).gt.0.0_dp) THEN
                    slopex=SIGN(MIN(2.0_dp*MIN(ABS(sigma_l),            &
     &                                         ABS(sigma_r)),           &
     &                              ABS(sigma_c)),                      &
     &                          sigma_c)
                  ELSE
                    slopex=0.0_dp
                  ENDIF
!                                       meridional sea surface slope
!
                  sigma_r=0.5_dp*(ptr3d(ii+1,jj+1,blk)-                 &
     &                            ptr3d(ii+1,jj  ,blk)+                 &
     &                            ptr3d(ii  ,jj+1,blk)-                 &
     &                            ptr3d(ii  ,jj  ,blk))/dyt(i,j,blk)
                  sigma_l=0.5_dp*(ptr3d(ii+1,jj  ,blk)-                 &
     &                            ptr3d(ii+1,jj-1,blk)+                 &
     &                            ptr3d(ii  ,jj  ,blk)-                 &
     &                            ptr3d(ii  ,jj-1,blk))/dyt(i,j,blk)
                  sigma_c=0.5_dp*(sigma_r+sigma_l)
                  IF ((sigma_r*sigma_l).gt.0.0_dp) THEN
                    slopey=SIGN(MIN(2.0_dp*MIN(ABS(sigma_l),            &
     &                                         ABS(sigma_r)),           &
     &                              ABS(sigma_c)),                      &
     &                          sigma_c)
                  ELSE
                    slopey(i,j,blk)=0.0_dp
                  ENDIF
!                                       rotate onto local basis vectors
!
                  ss_tltx(i,j,blk)= slopex*COS(ANGLET(i,j,blk))+        &
     &                              slopey*SIN(ANGLET(i,j,blk))
                  ss_tlty(i,j,blk)=-slopex*SIN(ANGLET(i,j,blk))+        &
     &                              slopey*COS(ANGLET(i,j,blk))
!
                  CALL t2ugrid_vector (ss_tltx)
                  CALL t2ugrid_vector (ss_tlty)
                END DO
              END DO
            END DO
!
!  Ocean mixed layer depth (m), from OCN.
!
          CASE ('hmix', 'mixed_layer_depth')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  hmix(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Freezing/Melting potential (W m-2), from OCN.
!
          CASE ('frzmlt', 'freezing_melting_potential')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  frzmlt(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Sea surface temperature (Celsius), maybe not needed, from OCN.
!
          CASE ('sst', 'sea_surface_temperature')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  sst(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Sea surface salinity (maybe for mushy layer), from OCN.
!
          CASE ('sss', 's_surf', 's_surf_ppt')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  sss(i,j,blk)=Fval
                END DO
              END DO
            END DO
!
!  Zonal surface ocean current (m s-1), from OCN.  Needs to be
!  rotated from east/north to i,j after it is loaded.
!
          CASE ('Usur', 'uocn', 'ocn_current_zonal')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  uocn(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_current(1)=.TRUE.
!
!  Meridional surface ocean current (m s-1), from OCN.  Needs to be
!  rotated from east/north to i,j after it is loaded.
!
          CASE ('Vsur', 'vocn', 'ocn_current_merid')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                  Fval=scale*ptr3d(ii,jj,blk)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  uocn(i,j,blk)=Fval
                END DO
              END DO
            END DO
            got_current(2)=.TRUE.
!
!  Import field not found.
!
          CASE DEFAULT
            IF (localPET.eq.0) THEN
              WRITE (cplout,10) TRIM(ImportNameList(ifld)),             &
     &                          TRIM(Time_CurrentString),               &
     &                          TRIM(CinpName)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     MyFile)) THEN
              rc=ESMF_RC_NOT_FOUND
              RETURN
            END IF
        END SELECT
!
!  Print pointer information.
!
        IF (DebugLevel.eq.4) THEN
          WRITE (cplout,20) localPET                                    &
     &                      LBOUND(ptr3d, DIM=1), UBOUND(ptr3d, DIM=1), &
     &                      LBOUND(ptr3d, DIM=2), UBOUND(ptr3d, DIM=2), &
     &                      LBOUND(ptr3d, DIM=3), UBOUND(ptr3d, DIM=3)
        END IF
!
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!
        IF (associated(ptr3d)) nullify (ptr3d)
!
!  Get import field minimun and maximum values.
!
        CALL ESMF_VMAllReduce (vm,                                      &
     &                         sendData=MyFmin,                         &
     &                         recvData=Fmin,                           &
     &                         count=2,                                 &
     &                         reduceflag=ESMF_REDUCE_MIN,              &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        CALL ESMF_VMAllReduce (vm,                                      &
     &                         sendData=MyFmax,                         &
     &                         recvData=Fmax,                           &
     &                         count=2,                                 &
     &                         reduceflag=ESMF_REDUCE_MAX,              &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Write out import field information.
!
        IF ((DebugLevel.ge.0).and.(localPET.eq.0)) THEN
          WRITE (cplout,30) TRIM(ImportNameList(ifld)),                 &
     &                      TRIM(Time_CurrentString), ng,               &
     &                      Fmin(1), Fmax(1)
          IF (ciceScale.ne.1.0_dp) THEN
            WRITE (cplout,40) Fmin(2), Fmax(2),                         &
     &                        ' ciceScale = ', ciceScale
          ELSE IF (add_offset.ne.0.0_dp) THEN
            WRITE (cplout,40) Fmin(2), Fmax(2),                         &
     &                        ' AddOffset = ', add_offset
          END IF
        END IF
!
!  Debugging: write out import field into NetCDF file.
!
        IF ((DebugLevel.ge.3).and.                                      &
     &      MODELS(Iseaice)%ImportField(ifld)%debug_write) THEN
          WRITE (ofile,50) ng, TRIM(ImportNameList(ifld)),              &
     &                     year, month, day, hour, minutes, seconds
          CALL ESMF_FieldWrite (field,                                  &
     &                          TRIM(ofile),                            &
     &                          overwrite=.TRUE.,                       &
     &                          rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
        END IF

      END DO FLD_LOOP
!
!  Deallocate local arrays.
!
      IF (allocated(ImportNameList)) deallocate (ImportNameList)
!
!  Update ROMS import calls counter.
!
      IF (ImportCount.gt.0) THEN
        MODELS(Iseaice)%ImportCalls=MODELS(Iseaice)%ImportCalls+1
      END IF
!
!  Rotate wind components from east/north to i,j coordinates and compute
!  wind magnitude.
!
      IF (ALL(got_wind)) THEN
        DO blk=1,nblocks
          my_block=get_block(blocks_ice(blk), blk)
          DO j=my_block%jlo,my_block%jhi
            DO i=my_block%ilo,my_block%ihi
              uvel=uatm(i,j,blk)
              vvel=vatm(i,j,blk)
              uatm(i,j,blk)= uvel*COS(ANGLET(i,j,blk))+                 &
     &                       vvel*SIN(ANGLET(i,j,blk))
              vatm(i,j,blk)=-uvel*SIN(ANGLET(i,j,blk))+                 &
     &                       vvel*COS(ANGLET(i,j,blk))
              wind(i,j,blk)=SQRT(uvel*uvel+vvel*vvel)
            END DO
          END DO
        END DO
      END IF
!
!  Rotate wind stress components from east/north to i,j coordinates.
!
      IF (ALL(got_wind)) THEN
        DO blk=1,nblocks
          my_block=get_block(blocks_ice(blk), blk)
          DO j=my_block%jlo,my_block%jhi
            DO i=my_block%ilo,my_block%ihi
              uvel=strax(i,j,blk)
              vvel=stray(i,j,blk)
              strax(i,j,blk)= uvel*COS(ANGLET(i,j,blk))+                &
     &                        vvel*SIN(ANGLET(i,j,blk))
              stray(i,j,blk)=-uvel*SIN(ANGLET(i,j,blk))+                &
     &                        vvel*COS(ANGLET(i,j,blk))
            END DO
          END DO
        END DO
      END IF
!
!  Rotate ocean current components from east/north to i,j coordinates.
!
      IF (ALL(got_current)) THEN
        DO blk=1,nblocks
          my_block=get_block(blocks_ice(blk), blk)
          DO j=my_block%jlo,my_block%jhi
            DO i=my_block%ilo,my_block%ihi
              uvel=uocn(i,j,blk)
              vvel=vocn(i,j,blk)
              uocn(i,j,blk)= uvel*COS(ANGLET(i,j,blk))+                 &
     &                       vvel*SIN(ANGLET(i,j,blk))
              vocn(i,j,blk)=-uvel*SIN(ANGLET(i,j,blk))+                 &
     &                       vvel*COS(ANGLET(i,j,blk))
            END DO
          END DO
        END DO
      END IF
!
!  Compute potential air temperature (K).
!
      IF (got_pair.and.got_tair) THEN
        DO blk=1,nblocks
          my_block=get_block(blocks_ice(blk), blk)
          DO j=my_block%jlo,my_block%jhi
            DO i=my_block%ilo,my_block%ihi
               potT(i,j,blk)=Tair(i,j,blk)*                             &
     &                       (100000.0_dp/Pair(i,j,blk))**0.286_dp
            END DO
          END DO
        END DO
      END IF
!
!  Compute net incomming shortwave radiation (W m-2).

      IF (ALL(got_swfx)) THEN
        DO blk=1,nblocks
          my_block=get_block(blocks_ice(blk), blk)
          DO j=my_block%jlo,my_block%jhi
            DO i=my_block%ilo,my_block%ihi
              fsw(i,j,blk)=swvdr(i,j,blk)+                              &
     &                     swvdf(i,j,blk)+                              &
     &                     swidr(i,j,blk)+                              &
     &                     swidf(i,j,blk)
            END DO
          END DO
        END DO
      END IF
!
!  Deallocate arrays.
!
      IF (allocated(ImportNameList)) deallocate (ImportNameList)
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_Import',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (/,3x,' CICE_Import - unable to find option to import: ',  &
     &        a,t68,a,/,18x,'check ''Import(roms)'' in input script: ', &
     &        a)
  20  FORMAT (18x,'PET [',i3.3,'],  Pointer Size: ',6i8)
  30  FORMAT (3x,' CICE_Import - ESMF: importing field ''',a,'''',      &
     &        t72,a,2x,'Grid ',i2.2,                                    &
     &        /,19x,'(Dmin = ', 1p,e15.8,0p,' Dmax = ',1p,e15.8,0p,')')
  40  FORMAT (19x,'(Cmin = ', 1p,e15.8,0p,' Cmax = ',1p,e15.8,0p,       &
     &        a,1p,e15.8,0p,')')
  50  FORMAT ('cice_',i2.2,'_import_',a,'_',i4.4,2('-',i2.2),'_',       &
     &        i2.2,2('.',i2.2),'.nc')

      RETURN
      END SUBROUTINE CICE_Import
!
      SUBROUTINE CICE_Export (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Exports CICE fields to other coupled gridded components.            !
!                                                                      !
!=======================================================================
!
      USE ice_blocks,       ONLY : block
      USE ice_blocks,       ONLY : get_block
      USE ice_constants,    ONLY : Tffresh
      USE ice_domain,       ONLY : nblocks, blocks_ice
      USE ice_grid,         ONLY : hm, ANGLET
      USE ice_flux
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: id, ifld
      integer :: blk, i, ii, j, jj
      integer :: ExportCount
      integer :: localPET
      integer :: year, month, day, hour, minutes, seconds, sN, SD
!
      real (dp) :: Fmin(1), Fmax(1), Fval, MyFmin(1), MyFmax(1)
      real (dp) :: wdir
!
      real (dp), pointer :: ptr3d(:,:,:) => NULL()
!
      character (len=22)      :: Time_CurrentString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CICE_Export"

      character (ESMF_MAXSTR) :: cname, ofile
      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
!
      TYPE (block)      :: my_block
      TYPE (ESMF_Field) :: field
      TYPE (ESMF_Time)  :: CurrentTime
      TYPE (ESMF_VM)    :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering CICE_Export',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get information about the gridded component.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       localPet=localPET,                         &
     &                       vm=vm,                                     &
     &                       name=cname,                                &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get current time.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (ClockInfo(Iseaice)%Clock,                     &
     &                    currTime=CurrentTime,                         &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeGet (CurrentTime,                                   &
     &                   yy=year,                                       &
     &                   mm=month,                                      &
     &                   dd=day,                                        &
     &                   h =hour,                                       &
     &                   m =minutes,                                    &
     &                   s =seconds,                                    &
     &                   sN=sN,                                         &
     &                   sD=sD,                                         &
     &                   timeString=Time_CurrentString,                 &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(Time_CurrentString, 'T')              ! remove 'T' in
      IF (is.gt.0) Time_CurrentString(is:is)=' '     ! ISO 8601 format
!
!-----------------------------------------------------------------------
!  Get list of export fields.
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet (MODELS(Iseaice)%ExportState(ng),              &
     &                    itemCount=ExportCount,                        &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (.not. allocated(ExportNameList)) THEN
        allocate ( ExportNameList(ExportCount) )
      END IF
!
      CALL ESMF_StateGet (MODELS(Iseaice)%ExportState(ng),              &
     &                    itemNameList=ExportNameList,                  &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Load export fields.
!-----------------------------------------------------------------------
!
      FLD_LOOP : DO ifld=1,ExportCount
!
!   Get field from export state.
!
        CALL ESMF_StateGet (MODELS(Iseaice)%ExportState(ng),            &
     &                      TRIM(ExportNameList(ifld)),                 &
     &                      field,                                      &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get field pointer.
!
        CALL ESMF_FieldGet (field,                                      &
     &                      farrayPtr=ptr3d,                            &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Initialize pointer to missing value.
!
        ptr3d=MISSING_dp
        Fmin(1)= MISSING_dp
        Tmax(1)=-MISSING_dp
!
!  Load field data into export state.
!
        SELECT CASE (TRIM(ADJUSTL(ExportNameList(ifld))))
!
!  Ice mask at cell center (T-cell), computed from land-boundary mask.
!
          CASE ('mask', 'hm', 'ice_mask')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  IF (hm(i,j,blk).gt.0.5_dp) THEN
                    ptr3d(ii,jj,blk)=1.0_dp
                  ELSE
                    ptr3d(ii,jj,blk)=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Fractional ice area (nondimensional; 0.0 - 1.0).
!
          CASE ('ifrac', 'ice_fraction')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=aice(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Surface temperature of ice/snow covered portion (K), to ATM.
!
          CASE ('sit', 'sea_ice_temperature')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  IF (aice(i,j,blk).gt.0.0_dp) THEN
                    ptr3d(ii,jj,blk)=Tffresh+trcr(i,j,1,blk)
                  ELSE
                    ptr3d(ii,jj,blk)=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Fraction of visible band, direct albedo aggregated over ice
!  categories (nondimesional), to ATM.
!
          CASE ('alvdr', 'inst_ice_vis_dir_albedo')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=alvdr(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Fraction of visible band, diffusive albedo aggregated over ice
!  categories (nondimesional), to ATM.
!
          CASE ('alvdf', 'inst_ice_vis_dif_albedo')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=alvdf(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Fraction of near-infrared band, direct albedo aggregated over
!  ice categories (nondimesional), to ATM.
!
          CASE ('alidr', 'inst_ice_ir_dir_albedo')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=alidr(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Fraction of near-infrared band, diffusive albedo aggregated over
!  ice categories (nondimesional), to ATM.
!
          CASE ('alidf', 'inst_ice_ir_dif_albedo')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=alidf(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Shortwave flux penetrating through ice to ocean (W m-2), to OCN.
!
          CASE ('fswthru', 'sw_pen_to_ocean')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=fswthru(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Visible direct band of net shortwave flux penetrating through
!  ice to ocean (W m-2), to OCN.
!
          CASE ('fswthruvdr', 'net_sw_vis_dir_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=fswthruvdr(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Visible diffusive band of net shortwave flux penetrating through
!  ice to ocean (W m-2), to OCN.
!
          CASE ('fswthruvdf', 'net_sw_vis_dif_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=fswthruvdf(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Infrared direct band of net shortwave flux penetrating through
!  ice to ocean (W m-2), to OCN.
!
          CASE ('fswthruidr', 'net_sw_ir_dir_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=fswthruidr(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Infrared diffusive band of net shortwave flux penetrating through
!  ice to ocean (W m-2), to OCN.
!
          CASE ('fswthruidf', 'net_sw_ir_dif_flx')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=fswthruidf(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Outgoing upward longwave ratiation (W m-2), averaged over ice
!  fraction only, to ATM.
!
          CASE ('flwout', 'mean_up_lw_flx_ice')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=flwout(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Ice sensible heat flux (W m-2), to ATM.
!
          CASE ('fsens', 'mean_sensi_heat_flx_atm_into_ice')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=fsens(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Ice latent heat flux (W m-2), to ATM.
!
          CASE ('flat', 'mean_laten_heat_flx_atm_into_ice')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=flat(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Evaporative water flux (kg m-2 s-1), to ATM.
!
          CASE ('evap', 'mean_evap_rate_atm_into_ice')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=evap(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Net heat flux to ocean (W m-2), to OCN.
!
          CASE ('fhocn', 'net_heat_flx_to_ocn')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=fhocn(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Fresh water flux to ocean (kg m-2 s-1), to OCN.
!
          CASE ('fresh', 'fresh_water_flx_to_ocean')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=fresh(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Salt flux to ocean (kg m-2 s-1), to OCN.
!
          CASE ('fsalt', 'salt_flx_to_ocean')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=fsalt(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Ice volume per unit area (m).
!
          CASE ('vice', 'mean_ice_volume')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=vice(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Snow volume per unit area (m).
!
          CASE ('vsno', 'mean_snow_volume')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ptr3d(ii,jj,blk)=vsno(i,j,blk)
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Zonal stress on ice by air (N m-2), to ATM.
!
          CASE ('strairxT', 'stress_on_air_ice_zonal')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ui=strairxT(i,j,blk)
                  vj=strairyT(i,j,blk)
                  ptr3d(ii,jj,blk)=ui*COS(ANGLET(i,j,blk))-             &
     &                             vj*SIN(ANGLET(i,j,blk))
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Meridional stress on ice by air (N m-2), to ATM.
!
          CASE ('strairyT', 'stress_on_air_ice_merid')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ui=strairxT(i,j,blk)
                  vj=strairyT(i,j,blk)
                  ptr3d(ii,jj,blk)=ui*SIN(ANGLET(i,j,blk))+             &
     &                             vj*COS(ANGLET(i,j,blk))
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Zonal stress on ice by ocean (N m-2), to OCN.
!
          CASE ('strocnxT', 'stress_on_ocn_ice_zonal')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ui=-strocnxT(i,j,blk)
                  vj=-strocnyT(i,j,blk)
                  ptr3d(ii,jj,blk)=ui*COS(ANGLET(i,j,blk))-             &
     &                             vj*SIN(ANGLET(i,j,blk))
                  MyFmin(1)=MIN(MyFmin(1),ptr3d(ii,jj,blk))
                  MyFmax(1)=MAX(MyFmax(1),ptr3d(ii,jj,blk))
                END DO
              END DO
            END DO
!
!  Meridional stress on ice by ocean (N m-2), to OCN.
!
          CASE ('strocnyT', 'stress_on_ocn_ice_merid')
            DO blk=1,nblocks
              my_block=get_block(blocks_ice(blk), blk)
              DO j=my_block%jlo,my_block%jhi
                jj=j-my_block%jlo+1
                DO i=my_block%ilo,my_block%ihi
                  ii=i-my_block%ilo+1
                  ui=-strocnxT(i,j,blk)
                  vj=-strocnyT(i,j,blk)
                  ptr(ii,jj,blk)=ui*SIN(ANGLET(i,j,blk))+               &
     &                           vj*COS(ANGLET(i,j,blk))
                END DO
              END DO
            END DO
!
!  Export field not found.
!
          CASE DEFAULT
            IF (localPET.eq.0) THEN
              WRITE (cplout,10) TRIM(ADJUSTL(ExportNameList(ifld))),    &
     &                          TRIM(CinpName)
            END IF
            rc=ESMF_RC_NOT_FOUND
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
        END SELECT
!
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!
          IF (associated(ptr3d)) nullify (ptr3d)
        END DO DE_LOOP
!
!  Get export field minimun and maximum values.
!
        CALL ESMF_VMAllReduce (vm,                                      &
     &                         sendData=MyFmin,                         &
     &                         recvData=Fmin,                           &
     &                         count=1,                                 &
     &                         reduceflag=ESMF_REDUCE_MIN,              &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        CALL ESMF_VMAllReduce (vm,                                      &
     &                         sendData=MyFmax,                         &
     &                         recvData=Fmax,                           &
     &                         count=1,                                 &
     &                         reduceflag=ESMF_REDUCE_MAX,              &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) TRIM(ExportNameList(ifld)),                 &
     &                      TRIM(Time_CurrentString), ng,               &
     &                      Fmin(1), Fmax(1)
        END IF
!
!  Debugging: write out field into a NetCDF file.
!
        IF ((DebugLevel.ge.3).and.                                      &
     &      MODELS(Iseaice)%ExportField(ifld)%debug_write) THEN
          WRITE (ofile,10) ng, TRIM(ExportNameList(ifld)),              &
     &                     year, month, day, hour, minutes, seconds
          CALL ESMF_FieldWrite (field,                                  &
     &                          TRIM(ofile),                            &
     &                          overwrite=.TRUE.,                       &
     &                          rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
        END IF
      END DO FLD_LOOP
!
!  Deallocate local arrays.
!
      IF (allocated(ExportNameList)) deallocate (ExportNameList)
!
!  Update CICE export calls counter.
!
      IF (ExportCount.gt.0) THEN
        MODELS(Iseaice)%ExportCalls=MODELS(Iseaice)%ExportCalls+1
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  CICE_Export',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      CALL my_flush (cplout)
!
  10  FORMAT (/,3x,' CICE_Export - unable to find option to export: ',  &
     &        a,/,18x,'check ''Export(cice)'' in input script: ',a)
  20  FORMAT (3x,' CICE_Export - ESMF: exporting field ''',a,'''',      &
     &        t72,a,2x,'Grid ',i2.2,/,                                  &
     &        18x,'(Cmin = ', 1p,e15.8,0p,' Cmax = ',1p,e15.8,0p,')')
  30  FORMAT ('cice_',i2.2,'_export_',a,'_',i4.4,2('-',i2.2),'_',       &
     &        i2.2,2('.',i2.2),'.nc')

      RETURN
      END SUBROUTINE CICE_Export
!
#endif
      END MODULE esmf_cice_mod
