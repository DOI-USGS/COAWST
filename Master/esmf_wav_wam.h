      MODULE esmf_wam_mod

#if defined WAM_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_wav_wam.h 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module sets WAM as the wave model gridded component using      !
!  the generic ESMF/NUOPC layer:                                       !
!                                                                      !
!    WAV_SetServices         Sets WAV component shared-object entry    !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                                                                      !
!    WAM_SetInitializeP1     WAM component phase 1 initialization:     !
!                            sets import and export fields long and    !
!                            short names into its respective state.    !
!                                                                      !
!    WAM_SetInitializeP2     WAM component phase 2 initialization:     !
!                            Initializes component (WAM_Initialize),   !
!                            sets component grid (WAM_SetGridArrays),  !
!                            and adds fields into import and export    !
!                            into respective states (WAM_SetStates).   !
!                                                                      !
!    WAM_DataInit            Exports WAM component fields during       !
!                            initialization or restart.                !
!                                                                      !
!    WAM_SetClock            Sets WAM component date calendar, start   !
!                            and stop times, and coupling interval.    !
!                                                                      !
!    WAM_CheckImport         Checks if WAM component import field is   !
!                            at the correct time.                      !
!                                                                      !
!    WAM_SetGridArrays       Sets WAM component horizontal grid        !
!                            arrays, grid area, and land/sea mask.     !
!                                                                      !
!    WAM_SetStates           Adds WAM component export and import      !
!                            fields into its respective state.         !
!                                                                      !
!    WAM_ModelAdvance        Advances WAM component for a coupling     !
!                            interval. It calls import and export      !
!                            routines.                                 !
!                                                                      !
!    WAM_SetFinalize         Finalizes WAM component execution.        !
!                                                                      !
!    WAM_import              Imports fields into WAM from other        !
!                            gridded components.                       !
!                                                                      !
!    WAM_Export              Exports WAM fields to other gridded       !
!                            components.                               !
!                                                                      !
!    WAM_Unpack              Unpacks WAM component export field by     !
!                            collecting data from each MPI node.       !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!  WAM:    ECMWF Wave Model (WAM), Cycle_4.5.3_MPI modified by RegESM  !
!            https://github.com/uturuncoglu/RegESM                     !
!                                                                      !
!=======================================================================
!
      USE ESMF
      USE NUOPC
      USE NUOPC_Model,                                                  &
     &    NUOPC_SetServices          => SetServices,                    &
     &    NUOPC_Label_Advance        => label_Advance,                  &
     &    NUOPC_Label_DataInitialize => label_DataInitialize,           &
     &    NUOPC_Label_SetClock       => label_SetClock,                 &
     &    NUOPC_Label_CheckImport    => label_CheckImport
!
      USE mod_esmf_esm          ! ESM coupling structures and variables
!
      USE wam_user_interface, ONLY : WAM_Initialize => WAM_init,        &
     &                               WAM_Run,                           &
     &                               WAM_Finalize
!
      implicit none
!
      PUBLIC  :: WAV_SetServices

      PRIVATE :: WAM_SetInitializeP1
      PRIVATE :: WAM_SetInitializeP2
      PRIVATE :: WAM_DataInit
      PRIVATE :: WAM_SetClock
      PRIVATE :: WAM_CheckImport
      PRIVATE :: WAM_SetGridArrays
      PRIVATE :: WAM_SetStates
      PRIVATE :: WAM_ModelAdvance
      PRIVATE :: WAM_SetFinalize
      PRIVATE :: WAM_Import
      PRIVATE :: WAM_Export
      PRIVATE :: WAM_Unpack
!
      CONTAINS
!
      SUBROUTINE WAV_SetServices (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets WAM component shared-object entry points for "initialize",     !
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
     &  __FILE__//", WAV_SetServices"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAV_SetServices',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
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
     &                              userRoutine=WAM_SetInitializeP1,    &
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
     &                              userRoutine=WAM_SetInitializeP2,    &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Attach WAM component phase independent specializing methods.
!-----------------------------------------------------------------------
!
!  Set routine for export initial/restart fields.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_DataInitialize,  &
     &                           specRoutine=WAM_DataInit,              &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for setting WAM clock.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_SetClock,        &
     &                           specRoutine=WAM_SetClock,              &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for checking import state.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_CheckImport,     &
     &                           specPhaseLabel="RunPhase1",            &
     &                           specRoutine=WAM_CheckImport,           &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for time-stepping WAM component.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_Advance,         &
     &                           specRoutine=WAM_ModelAdvance,          &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Register WAM finalize routine.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompSetEntryPoint (model,                           &
     &                                 methodflag=ESMF_METHOD_FINALIZE, &
     &                                 userRoutine=WAM_SetFinalize,     &
     &                                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAV_SetServices',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_SetServices
!
      SUBROUTINE WAM_SetInitializeP1 (model,                            &
     &                                ImportState, ExportState,         &
     &                                clock, rc)
!
!=======================================================================
!                                                                      !
!  WAM component Phase 1 initialization: sets import and export        !
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
     &  __FILE__//", WAM_SetInitializeP1"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_SetInitializeP1',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set WAM import state and fields.
!-----------------------------------------------------------------------
!
!  Add WAM import state(s). If nesting, each grid has its own import
!  state.
!
      IMPORTING : IF (Nimport(Iwave).gt.0) THEN
        DO ng=1,MODELS(Iwave)%Ngrids
          IF (ANY(COUPLED(Iwave)%LinkedGrid(ng,:))) THEN
            CoupledSet=TRIM(COUPLED(Iwave)%SetLabel(ng))
            StateLabel=TRIM(COUPLED(Iwave)%ImpLabel(ng))
            CALL NUOPC_AddNestedState (ImportState,                     &
     &                                 CplSet=TRIM(CoupledSet),         &
     &                                 nestedStateName=TRIM(StateLabel),&
     &                                 nestedState=MODELS(Iwave)%       &
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
            DO i=1,Nimport(Iwave)
              StandardName=MODELS(Iwave)%ImportField(i)%standard_name
              ShortName   =MODELS(Iwave)%ImportField(i)%short_name
              CALL NUOPC_Advertise (MODELS(Iwave)%ImportState(ng),      &
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
!  Set WAM export state and fields.
!-----------------------------------------------------------------------
!
!  Add WAM import state. If nesting, each grid has its own import
!  state.
!
      EXPORTING : IF (Nexport(Iwave).gt.0) THEN
        DO ng=1,MODELS(Iwave)%Ngrids
          IF (ANY(COUPLED(Iwave)%LinkedGrid(ng,:))) THEN
            CoupledSet=TRIM(COUPLED(Iwave)%SetLabel(ng))
            StateLabel=TRIM(COUPLED(Iwave)%ExpLabel(ng))
            CALL NUOPC_AddNestedState (ExportState,                     &
     &                                 CplSet=TRIM(CoupledSet),         &
     &                                 nestedStateName=TRIM(StateLabel),&
     &                                 nestedState=MODELS(Iwave)%       &
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
            DO i=1,Nexport(Iwave)
              StandardName=MODELS(Iwave)%ExportField(i)%standard_name
              ShortName   =MODELS(Iwave)%ExportField(i)%short_name
              CALL NUOPC_Advertise (MODELS(Iwave)%ExportState(ng),      &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_SetInitializeP1',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_SetInitializeP1
!
      SUBROUTINE WAM_SetInitializeP2 (model,                            &
     &                                ImportState, ExportState,         &
     &                                clock, rc)
!
!=======================================================================
!                                                                      !
!  WAM component Phase 2 initialization: Initializes WAM, sets         !
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
     &  __FILE__//", WAM_SetInitializeP2"
!
      TYPE (ESMF_VM) :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_SetInitializeP2',     &
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
!  Initialize WAM component.
!-----------------------------------------------------------------------
!
      CALL WAM_Initialize (MyComm)
!
!-----------------------------------------------------------------------
!  Set-up grid and load coordinate data.
!-----------------------------------------------------------------------
!
      DO ng=1,MODELS(Iwave)%Ngrids
        IF (ANY(COUPLED(Iwave)%LinkedGrid(ng,:))) THEN
          CALL WAM_SetGridArrays (ng, model, localPET, rc)
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
      DO ng=1,MODELS(Iwave)%Ngrids
        IF (ANY(COUPLED(Iwave)%LinkedGrid(ng,:))) THEN
          CALL WAM_SetStates (ng, model, rc)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_SetInitializeP2',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_SetInitializeP2
!
      SUBROUTINE WAM_DataInit (model, rc)
!
!=======================================================================
!                                                                      !
!  Exports WAM component fields during initialization or restart.      !
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
     &  __FILE__//", WAM_DataInit"
!
      TYPE (ESMF_Time)  :: CurrentTime
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_DataInit',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get gridded component clock current time.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (ClockInfo(Iwave)%Clock,                       &
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
      IF (Nexport(Iwave).gt.0) THEN
        DO ng=1,MODELS(Iwave)%Ngrids
          IF (ANY(COUPLED(Iwave)%LinkedGrid(ng,:))) THEN
            CALL WAM_Export (ng, model, rc=rc)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_DataInit',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_DataInit
!
      SUBROUTINE WAM_SetClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets WAM component date calendar, start and stop time, and coupling !
!  interval.                                                           !
!                                                                      !
!=======================================================================
!
      USE wam_timopt_module, ONLY : cdatea, cdatee, coldstart
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: ref_year,   start_year,   stop_year
      integer :: ref_month,  start_month,  stop_month
      integer :: ref_day,    start_day,    stop_day
      integer :: ref_hour,   start_hour,   stop_hour
      integer :: ref_minute, start_minute, stop_minute
      integer :: ref_second, start_second, stop_second
      integer :: localPET, PETcount
      integer :: TimeFrac, ig
!
      real(r8) :: hour, minute, yday
!
      character (len= 80) :: Calendar
      character (len=160) :: message

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WAM_SetClock"
!
      TYPE (ESMF_CalKind_Flag) :: CalType
      TYPE (ESMF_Time)         :: CurrentTime, StartTime
      TYPE (ESMF_TimeInterval) :: TimeStep
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_SetClock',
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
!  Create WAM component clock.
!-----------------------------------------------------------------------
!
!  Set calendar.
!
      Calendar=TRIM(ClockInfo(Iwave)%CalendarString)
      IF (TRIM(Calendar).eq.'gregorian') THEN
        CalType=ESMF_CALKIND_GREGORIAN
      ELSE
        CalType=ESMF_CALKIND_GREGORIAN
      END IF
      ClockInfo(Iwave)%Calendar=ESMF_CalendarCreate(CalType,            &
     &                                              name=TRIM(calendar),&
     &                                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set reference time.
!
      IF (coldstart) THEN
        READ (cdatea,'(i4,5i2)') ref_year, ref_month, ref_day,          &
     &                           ref_hour, ref_minute, ref_second
!
        CALL ESMF_TimeSet(ClockInfo(Iwave)%ReferenceTime,               &
     &                    yy=ref_year,                                  &
     &                    mm=ref_month,                                 &
     &                    dd=ref_day,                                   &
     &                    h=ref_hour,                                   &
     &                    m=ref_minute,                                 &
     &                    s=ref_second,                                 &
     &                    calendar=ClockInfo(Iwave)%Calendar,           &
     &                    rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
!
!  Set start time.
!
      READ (cdatea,'(i4,5i2)') start_year, start_month, start_day,      &
     &                         start_hour, start_minute, start_second
!
      CALL ESMF_TimeSet (ClockInfo(Iwave)%StartTime,                    &
                         yy=start_year,                                 &
                         mm=start_month,                                &
                         dd=start_day,                                  &
                         h=start_hour,                                  &
                         m=start_minute,                                &
                         s=start_second,                                &
                         calendar=ClockInfo(Iwave)%Calendar,            &
                         rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set stop time.
!
      READ (cdatee,'(i4,5i2)') stop_year, stop_month, stop_day,         &
     &                         stop_hour, stop_minute, stop_second
!
      CALL ESMF_TimeSet (ClockInfo(Iwave)%StopTime,                     &
     &                   yy=stop_year,                                  &
     &                   mm=stop_month,                                 &
     &                   dd=stop_day,                                   &
     &                   h=stop_hour,                                   &
     &                   m=stop_minute,                                 &
     &                   s=stop_second,                                 &
     &                   calendar=ClockInfo(Iwave)%Calendar,            &
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
     &                       clock=ClockInfo(Iwave)%Clock,              &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_ClockGet (ClockInfo(Iwave)%Clock,                       &
     &                    timeStep=ClockInfo(Iwave)%TimeStep,           &
     &                    currTime=ClockInfo(Iwave)%CurrentTime,        &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Compare driver time against WAM component time.
!-----------------------------------------------------------------------
!
      IF (ClockInfo(Idriver)%Restarted) THEN
        StartTime=ClockInfo(Idriver)%RestartTime
      ELSE
        StartTime=ClockInfo(Idriver)%StartTime
      END IF
!
      IF (ClockInfo(Iwave)%StartTime.ne.StartTime) THEN
        CALL ESMF_TimePrint (ClockInfo(Iwave)%StartTime,                &
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
        message='Driver and WAM start times do not match: '//           &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ClockInfo(Iwave  )%StopTime.ne.                               &
     &    ClockInfo(Idriver)%StopTime) THEN
        CALL ESMF_TimePrint (ClockInfo(Iwave)%StopTime,                 &
     &                       options="string",                          &
     &                       rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,
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
        message='Driver and WAM stop times do not match: '//            &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ClockInfo(Iwave  )%Calendar.ne.                               &
     &    ClockInfo(Idriver)%Calendar) THEN
        CALL ESMF_CalendarPrint (ClockInfo(Iwave)%Calendar,             &
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
        message='Driver and WAM calendars do not match: '//             &
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
      DO ig=1,MODELS(Iwave)%Ngrids
        TimeFrac=MAX(TimeFrac,                                          &
     &               MAXVAL(MODELS(Iwave)%TimeFrac(ig,:),               &
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
      ClockInfo(Iwave)%TimeStep=ClockInfo(Idriver)%TimeStep/TimeFrac
!
      IF (coldstart) THEN
        CALL ESMF_ClockSet (ClockInfo(Iwave)%Clock,                     &
     &                      name=TRIM(ClockInfo(Iwave)%Name),           &
     &                      refTime  =ClockInfo(Iwave)%ReferenceTime,   &
     &                      timeStep =ClockInfo(Iwave)%TimeStep,        &
     &                      startTime=ClockInfo(Iwave)%StartTime,       &
     &                      stopTime =ClockInfo(Iwave)%StopTime,        &
                            rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      ELSE
        CALL ESMF_ClockSet (ClockInfo(Iwave)%Clock,                     &
     &                      name=TRIM(ClockInfo(Iwave)%Name),           &
     &                      timeStep =ClockInfo(Iwave)%TimeStep,        &
     &                      startTime=ClockInfo(Iwave)%StartTime,       &
     &                      stopTime =ClockInfo(Iwave)%StopTime,        &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_SetClock',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_SetClock
!
      SUBROUTINE WAM_CheckImport (model, rc)
!
!=======================================================================
!                                                                      !
!  Checks if WAM component import field is at the correct time.        !
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
      logical :: atCorrectTime
!
      integer :: ImportCount, localPET
!
      character(ESMF_MAXSTR), allocatable :: ImportNameList(:)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WAM_CheckImport"
!
      TYPE (ESMF_Clock)        :: driverClock
      TYPE (ESMF_Field)        :: field
      TYPE (ESMF_State)        :: ImportState
      TYPE (ESMF_Time)         :: StartTime, CurrentTime
      TYPE (ESMF_TimeInterval) :: TimeStep
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_CheckImport',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Query component for the driver clock.
!-----------------------------------------------------------------------
!
      CALL NUOPC_ModelGet (model,                                       &
     &                     driverClock=driverClock,                     &
     &                     rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
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
     &                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get the start time and current time from driver clock.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (driverClock,                                  &
     &                    startTime=StartTime,                          &
     &                    currTime=CurrentTime,                         &
     &                    timeStep=TimeStep,                            &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Query WAM component for its import state.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       importState=ImportState,                   &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get list of import fields.
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet (MODELS(Iwave)%ImportState(ng),                &
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
!
      CALL ESMF_StateGet (MODELS(Iwave)%ImportState(ng),                &
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
!  Check fields in the ImportState object.
!-----------------------------------------------------------------------
!
      IF (ImportCount.gt.0) THEN
        CALL ESMF_StateGet (MODELS(Iwave)%ImportState(ng),              &
     &                      itemName=TRIM(ImportNameList(1)),           &
     &                      field=field,                                &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Check if import field is at the correct time.
!
        atCorrectTime=NUOPC_IsAtTime(field,                             &
     &                               CurrentTime,                       &
     &                               rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF

        CALL report_timestamp (field, CurrentTime,                      &
     &                         localPET, "WAM", rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        IF (.not.atCorrectTime) THEN
          CALL ESMF_LogSetError(ESMF_RC_ARG_BAD,                        &
     &                          msg="NUOPC INCOMPATIBILITY DETECTED:"// &
     &                          " Import Fields not at correct time",   &
     &                          line=__LINE__,                          &
     &                          file=MyFile,                    &
     &                          rcToReturn=rc)
          RETURN
        END IF
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_CheckImport',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_CheckImport
!
      SUBROUTINE WAM_SetGridArrays (ng, model, localPET, rc)
!
!=======================================================================
!                                                                      !
!  Sets WAM component staggered, horizontal grids arrays, grid area,   !
!  and land/sea mask, if any.                                          !
!                                                                      !
!=======================================================================
!
      USE wam_grid_module, ONLY : nx, amowep, amoeap, xdello
      USE wam_grid_module, ONLY : ny, amosop, amonop, xdella, l_s_mask
      USE wam_mpi_module , ONLY : petotal, irank, nstart, nend
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng, localPET
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp), intent(inout) :: model
!
!  Local variable declarations.
!
      integer :: i, ivar, j, localDECount, tile
      integer :: imin, imax, jmin, jmax
!
      integer, allocatable :: deBlockList(:,:,:)
!
      integer (i4b), pointer :: ptrM(:,:) => NULL()
!
      real (r8), pointer :: ptrX(:,:) => NULL()
      real (r8), pointer :: ptrY(:,:) => NULL()
!
      character (ESMF_MAXSTR) :: name

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WAM_SetGridArrays"
!
      TYPE (ESMF_DistGrid)   :: distGrid1, distGrid2
      TYPE (ESMF_StaggerLoc) :: staggerLoc
      TYPE (ESMF_VM)         :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_SetGridArrays',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get limits of the grid arrays (based on PET and nest level)
!-----------------------------------------------------------------------
!
      IF (.not. allocated(deBlockList)) THEN
        allocate ( deBlockList(1,2,petotal) )
      END IF
!
      DO tile=1,petotal
        deBlockList(1,1,tile)=nstart(tile)
        deBlockList(1,2,tile)=nend(tile)
      END DO
!
!-----------------------------------------------------------------------
!  Create ESMF DistGrid based on model domain decomposition.
!-----------------------------------------------------------------------
!
      distGrid1=ESMF_DistGridCreate(minIndex=(/ 1 /),                   &
     &                              maxIndex=(/ nend(petotal) /),       &
     &                              deBlockList=deBlockList,            &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      distGrid2=ESMF_DistGridCreate(minIndex=(/ 1, 1 /),                &
     &                              maxIndex=(/ nx, ny /),              &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Set component grid coordinates.
!-----------------------------------------------------------------------
!
!  Define component grid location type.
!
      IF (.not.allocated(MODELS(Iwave)%mesh)) THEN
        allocate ( MODELS(Iwave)%mesh(1) )
        MODELS(Iwave)%mesh(1)%gtype = Icenter
      END IF
!
!  Create ESMF Grid. The array indices are global.
!
      MODELS(Iwave)%grid(ng)=ESMF_GridCreate(distgrid=distGrid2,        &
     &                                     indexflag=ESMF_INDEX_GLOBAL, &
     &                                       name="wam_grid",           &
     &                                       rc=rc)
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
      CALL ESMF_GridGet (MODELS(Iwave)%grid(ng),                        &
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
      MESH_LOOP : DO ivar=1,UBOUND(MODELS(Iwave)%mesh, DIM=1)
!
        SELECT CASE (MODELS(Iwave)%mesh(ivar)%gtype)
          CASE (Icenter)
            staggerLoc=ESMF_STAGGERLOC_CENTER
        END SELECT
!
!  Allocate coordinate storage associated with staggered grid type.
!  No coordinate values are set yet.
!
        CALL ESMF_GridAddCoord (MODELS(Iwave)%grid(ng),                 &
     &                          staggerLoc=staggerLoc                   &
     &                          rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Allocate storage for land/sea masking.
!
        CALL ESMF_GridAddItem (MODELS(Iwave)%grid(ng),                  &
     &                         staggerLoc=staggerLoc,                   &
     &                         itemflag=ESMF_GRIDITEM_MASK,             &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
        MODELS(Iwave)%LandValue=0
        MODELS(Iwave)%SeaValue=1
!
!  Get pointers and set coordinates for the grid.  Usually, the DO-loop
!  is executed once since localDEcount=1.
!
        DE_LOOP : DO localDE=0,localDEcount-1
          CALL ESMF_GridGetCoord (MODELS(Iwave)%grid(ng),               &
     &                            localDE=localDE,                      &
     &                            staggerLoc=staggerLoc,                &
     &                            coordDim=1,                           &
     &                            farrayPtr=ptrX,                       &
     &                            rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_GridGetCoord (MODELS(Iwave)%grid(ng),               &
     &                            localDE=localDE,                      &
     &                            staggerLoc=staggerLoc,                &
     &                            coordDim=2,                           &
     &                            farrayPtr=ptrY,                       &
     &                            rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_GridGetItem (MODELS(Iwave)%grid(ng),                &
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
!  Fill the pointers.
!
          imin=LBOUND(ptrX, DIM=1)
          imax=UBOUND(ptrX, DIM=1)
          jmin=LBOUND(ptrX, DIM=2)
          jmax=UBOUND(ptrX, DIM=2)
!
          SELECT CASE (MODELS(Iwave)%mesh(ivar)%gtype)
            CASE (Icenter)
              DO i=imin,imax
                ptrX(i,jmin:jmax)=REAL(i-1,r8)*xdello+amowep
              END DO
              DO j=jmin,jmax
                ptrY(imin:imax,j)=REAL(j-1,r8)*xdella+amosop
              END DO
              DO i=imin,imax
                DO j=jmin,jmax
                  IF (l_s_mask(i,j)) THEN
                    ptrM(i,j)=MODELS(Iwave)%SeaValue
                  ELSE
                    ptrM(i,j)=MODELS(Iwave)%LandValue
                  END IF
                END DO
              END DO
          END SELECT
!
!  Nullify pointers.
!
          IF (associated(ptrX)) nullify (ptrX)
          IF (associated(ptrY)) nullify (ptrY)
          IF (associated(ptrM)) nullify (ptrM)
        END DO DE_LOOP
!
!  Debugging: write out component grid in VTK format.
!
        IF (DebugLevel.ge.4) THEN
          gtype=MODELS(Iwave)%mesh(ivar)%gtype)
          CALL ESMF_GridWriteVTK (MODELS(Iwave)%grid(ng),               &
     &                            filename="wam_"//                     &
     &                                     TRIM(GridType(gtype)//       &
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
     &                       grid=MODELS(Iwave)%grid(ng),               &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_SetGridArrays',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_SetGridArrays
!
      SUBROUTINE WAM_SetStates (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Adds WAM component export and import fields into its respective     !
!  state.                                                              !
!                                                                      !
!=======================================================================
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
      integer :: localDE, localDEcount
      integer :: ExportCount, ImportCount
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WAM_SetStates"
!
      TYPE (ESMF_ArraySpec)  :: arraySpec
      TYPE (ESMF_Field)      :: field
      TYPE (ESMF_StaggerLoc) :: staggerLoc
      TYPE (ESMF_VM)         :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_SetStates',           &
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
!  Get number of local decomposition elements (DEs). Usually, a single
!  Decomposition Element (DE) is associated with each Persistent
!  Execution Thread (PETs). Thus, localDEcount=1.
!
      CALL ESMF_GridGet (MODELS(Iwave)%grid(ng),                        &
     &                   localDECount=localDEcount,                     &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Set a 2D floating-point array descriptor.
!-----------------------------------------------------------------------
!
      CALL ESMF_ArraySpecSet (arraySpec,                                &
     &                        typekind=ESMF_TYPEKIND_R8,                &
     &                        rank=2,                                   &
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
      EXPORTING : IF (Nexport(Iwave).gt.0) THEN
!
!  Get number of fields to export.
!
        CALL ESMF_StateGet (MODELS(Iwave)%ExportState(ng),              &
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
        CALL ESMF_StateGet (MODELS(Iwave)%ExportState(ng),              &
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
          id=field_index(MODELS(Iwave)%ExportField, ExportNameList(i))
!
          IF (NUOPC_IsConnected(MODELS(Iwave)%ExportState(ng),          &
     &                          fieldName=TRIM(ExportNameList(i)),      &
     &                          rc=rc)) THEN
!
!  Set staggering type.
!
            SELECT CASE (MODELS(Iwave)%ExportField(id)%gtype)
              CASE (Icenter)
                staggerLoc=ESMF_STAGGERLOC_CENTER
            END SELECT
!
!  Create 2D field from the Grid and arraySpec.
!
            field=ESMF_FieldCreate(MODELS(Iwave)%grid(ng),              &
     &                             arraySpec,                           &
     &                             staggerloc=staggerLoc,               &
     &                             name=TRIM(ExportNameList(i)),        &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Put data into state. Usually, the DO-loop is executed once since
!  localDEcount=1.
!
            DO localDE=0,localDEcount-1
!
!  Get pointer to DE-local memory allocation within field.
!
              CALL ESMF_FieldGet (field,                                &
     &                            localDe=localDE,                      &
     &                            farrayPtr=ptr2d,                      &
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
              ptr2d=MISSING_dp
!
!  Nullify pointer to make sure that it does not point on a random part
!  in the memory.
!
              IF ( associated(ptr2d) ) nullify (ptr2d)
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
     &                          TRIM(COUPLED(Iwave)%ExpLabel(ng))
            END IF
            CALL ESMF_StateRemove (MODELS(Iwave)%ExportState(ng),       &
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
      IMPORTING : IF (Nimport(Iwave).gt.0) THEN
!
!  Get number of fields to import.
!
        CALL ESMF_StateGet (MODELS(Iwave)%ImportState(ng),              &
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
        CALL ESMF_StateGet (MODELS(Iwave)%ImportState(ng),              &
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
          id=field_index(MODELS(Iwave)%ImportField, ImportNameList(i))
!
          IF (NUOPC_IsConnected(MODELS(Iwave)%ImportState(ng),          &
     &                          fieldName=TRIM(ImportNameList(i)),      &
     &                          rc=rc)) THEN

!
!  Set staggering type.
!
            SELECT CASE (MODELS(Iwave)%ImportField(id)%gtype)
              CASE (Icenter)
                staggerLoc=ESMF_STAGGERLOC_CENTER
            END SELECT
!
!  Create 2D field from the Grid, arraySpec. The array indices are
!  global.
!
            field=ESMF_FieldCreate(MODELS(Iwave)%grid(ng),              &
     &                             arraySpec,                           &
     &                             staggerloc=staggerLoc,               &
     &                             indexflag=ESMF_INDEX_GLOBAL,         &
     &                             name=TRIM(ImportNameList(i)),        &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Put data into state. Usually, the DO-loop is executed once since
!  localDEcount=1.
!
            DO localDE=0,localDEcount-1
!
!  Get pointer to DE-local memory allocation within field.
!
              CALL ESMF_FieldGet (field,                                &
     &                            localDe=localDE,                      &
     &                            farrayPtr=ptr2d,                      &
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
              ptr2d=MISSING_dp
!
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!
              IF (associated(ptr2d)) nullify (ptr2d)
            END DO
!
!  Add field import state.
!
            CALL NUOPC_Realize (MODELS(Iwave)%ImportState(ng),          &
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
     &                          TRIM(COUPLED(Iwave)%ImpLabel(ng))
            END IF
            CALL ESMF_StateRemove (MODELS(Iwave)%ImportState(ng),       &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_SetStates',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_SetStates
!
      SUBROUTINE WAM_ModelAdvance (model, rc)
!
!=======================================================================
!                                                                      !
!  Advance WAM component for a coupling interval (seconds) using       !
!  "WAM_Run". It also calls "WAM_Import" and "WAM_Export" to           !
!  import and export coupling fields, respectively.                    !
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
      integer :: localPET, PETcount, phase
!
      real (r8) :: trun
!
      character (ESMF_MAXSTR) :: str1, str2

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WAM_ModelAdvance"
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
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_ModelAdvance',        &
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
     &                       importState=ImportState,                   &
     &                       exportState=ExportState,                   &
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
!  Get ID for local PET and number of PETs.
!
      CALL ESMF_VMGet (vm,                                              &
     &                 localPet=localPET,                               &
     &                 petCount=PETcount,                               &
     &                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get driver time step interval, stopping time, reference time, and
!  current time.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (clock,                                        &
     &                    timeStep=TimeStep,                            &
     &                    stopTime=StopTime,                            &
     &                    refTime=RefTime,                              &
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
!  Get time interval (seconds, double precision).
!-----------------------------------------------------------------------
!
      CALL ESMF_TimeIntervalGet (TimeStep,                              &
     &                           s_r8=trun,                             &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Debugging: report time information string (YYYY-MM-DD hh:mm:ss).
!-----------------------------------------------------------------------
!
      IF ((DebugLevel.ge.0).and.(localPET.eq.0)) THEN
!
!  Current driver time.
!
        CALL ESMF_TimeGet (CurrentTime,                                 &
     &                     timeStringISOFrac=str1,                      &
     &                     rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Next driver coupling time.
!
        CALL ESMF_TimeGet (CurrentTime+TimeStep,                        &
     &                     timeStringISOFrac=str2,                      &
     &                     rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        IF (DebugLevel.eq.0) THEN
          WRITE (cplout,10) TRIM(str1), TRIM(str2), phase
        ELSE
          WRITE (cplout,20) TRIM(str1), TRIM(str2), phase, trun
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Get import fields.
!-----------------------------------------------------------------------
!
      IF ((Nimport(Iwave).gt.0).and.                                    &
     &    (CurrentTime.ne.RefTime).or.Restarted)) THEN
        DO ng=1,MODELS(Iwave)%Ngrids
          IF (ANY(COUPLED(Iwave)%LinkedGrid(ng,:))) THEN
            CALL WAM_Import (ng, model, rc=rc)
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
!-----------------------------------------------------------------------
!  Run WAM component.
!-----------------------------------------------------------------------
!
      CALL WAM_Run ()
!
!-----------------------------------------------------------------------
!  Put export fields.
!-----------------------------------------------------------------------
!
      IF (Nexport(Iwave).gt.0) THEN
        DO ng=1,MODELS(Iwave)%Ngrids
          IF (ANY(COUPLED(Iwave)%LinkedGrid(ng,:))) THEN
            CALL WAM_Export (ng, model, rc)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_ModelAdvance',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      CALL my_flush (cplout)
!
  10  FORMAT (/,' ESMF, Running WAM: ',a,' --> ',a,' Phase: ',i1)
  20  FORMAT (/,' ESMF, Running WAV: ',a,' --> ',a,' Phase: ',i1,       &
     &        ' [', f15.2, ' s]')

      RETURN
      END SUBROUTINE WAM_ModelAdvance
!
      SUBROUTINE WAM_SetFinalize (model,                                &
     &                            ImportState, ExportState,             &
     &                            clock, rc)
!
!=======================================================================
!                                                                      !
!  Finalize WAM component execution. It calls WAM_finalize.            !
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
     &  __FILE__//", WAM_SetFinalize"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_SetFinalize',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Finalize WAM component.
!-----------------------------------------------------------------------
!
      CALL WAM_Finalize ()
      CALL my_flush (6)                   ! flush standard output buffer
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_SetFinalize',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_SetFinalize
!
      SUBROUTINE WAM_Import (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Imports fields into WAM array structures. The fields aew loaded     !
!  into the snapshot storage arrays to allow time interpolation in     !
!  WAM kernel.                                                         !
!                                                                      !
!=======================================================================
!
      USE wam_grid_module,    ONLY : nx, ny, nsea, l_s_mask
      USE wam_user_interface, ONLY : us_esmf, vs_esmf
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
      integer :: id, ifld, tile
      integer :: iyear, iday, imonth, ihour
      integer :: localPET
!
      real (dp) :: add_offset, scale
!
      real (dp), allocatable :: arr2d(:,:)
!
      character (ESMF_MAXSTR) :: ofile
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WAM_Import"
!
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
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_Import',              &
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
     &                   yy=iyear,                                      &
     &                   mm=imonth,                                     &
     &                   dd=iday,                                       &
     &                   h=ihour,                                       &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get list of import fields.
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet (MODELS(Iwave)%ImportState(ng),                &
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
      CALL ESMF_StateGet (MODELS(Iwave)%ImportState(ng),                &
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
      FLD_LOOP : DO ifld=1,ImportCount
        id=field_index(MODELS(Iwave)%ImportField, ImportNameList(ifld))
!
!  Get field from import state.
!
        CALL ESMF_StateGet (MODELS(Iwave)%ImportState(ng),              &
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
!  Collect field from all PETs.
!
        IF (.not. allocated(arr2d)) THEN
          allocate ( arr2d(nx,ny) )
        END IF
!
        DO tile=0,PETcount-1
          CALL ESMF_FieldGather (field,                                 &
     &                           arr2d,                                 &
     &                           rootPet=tile,                          &
     &                           vm=vm,                                 &
     &                           rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
        END DO
!
!  Debugging: write size of pointer.
!
        IF (DebugLevel.gt.1) THEN
          WRITE (cplout,10) localPET, tile,                             &
     &                 ADJUSTL("IND/WAV/IMP/"//ImportNameList(ifld)),   &
     &                 1, nx, 1, ny
        END IF
!
!  Load import data into WAM component variable.
!  (HGA: It is kind of weird that everything is loaded into us_esmf
!   and vs_esmf)
!
        scale=MODELS(Iwave)%ImportField(id)%scale_factor
        add_offset=MODELS(Iwave)%ImportField(id)%add_offset
!
        SELECT CASE (TRIM(ADJUSTL(ImportNameList(ifld))))
!
!  Surface eastward wind component (m s-1).
!
          CASE ('wndu')
            us_esmf(1:nsea)=PACK(arr2d, l_s_mask)
            us_esmf(1:nsea)=(us_esmf(1:nsea)*scale)+add_offset
!
!  Surface northward wind component (m s-1).
!
          CASE ('wndv')
            vs_esmf(1:nsea)=PACK(arr2d, l_s_mask)
            vs_esmf(1:nsea)=(vs_esmf(1:nsea)*scale)+add_offset
!
!  Friction velocity (m s-1).
!
          CASE ('ustr')
            us_esmf(1:nsea)=PACK(arr2d, l_s_mask)
            us_esmf(1:nsea)=(us_esmf(1:nsea)*scale)+add_offset
!
!  Direction. (HGA: ???)
!
          CASE ('wdir')
            vs_esmf(1:nsea)=PACK(arr2d, l_s_mask)
            vs_esmf(1:nsea)=(vs_esmf(1:nsea)*scale)+add_offset
        END SELECT
!
!  Debugging: write field into a NetCDF file.
!
        IF ((DebugLevel.ge.3).and.                                      &
     &      MODELS(Iwave)%ImportField(ifld)%debug_write) THEN
          WRITE (ofile,20) 'wam_import', TRIM(ImportNameList(ifld)),    &
     &                     iyear, imonth, iday, ihour
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
!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------
!
      IF (allocated(ImportNameList)) deallocate (ImportNameList)
      IF (allocated(arr2d)) deallocate (arr2d)
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_Import',              &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (' PET(',I3,') - tile(',i2,') - ', a20, ' : ', 4i8)
  20  FORMAT (a,'_',a,'_',i4.4,3('-',i2.2),'.nc')

      RETURN
      END SUBROUTINE WAM_Import
!
      SUBROUTINE WAM_Export (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Exports WAM fields to other coupled gridded components.             !
!                                                                      !
!=======================================================================
!
      USE wam_model_module, ONLY : z0, ustar, tauw
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: ifld, imin, imax, jmin, jmax
      integer :: iyear, iday, imonth, ihour
      integer :: ExportCount
      integer :: localDE, localDEcount, localPET, PETcount
!
      real (r8), pointer :: ptr(:,:) => NULL()
!
      character (ESMF_MAXSTR) :: ofile
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WAM_Export"
!
      TYPE (ESMF_Clock) :: clock
      TYPE (ESMF_Field) :: field
      TYPE (ESMF_State) :: ExportState
      TYPE (ESMF_Time)  :: CurrentTime
      TYPE (ESMF_VM)    :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_Export',              &
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
     &                       vm=vm,                                     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get ID for local PET and number of PETs.
!
      CALL ESMF_VMGet (vm,                                              &
     &                 localPet=localPET,                               &
     &                 petCount=PETcount,                               &
     &                 rc=rc)
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
      CALL ESMF_GridGet (MODELS(Iwave)%grid(ng),                        &
     &                   localDECount=localDEcount,                     &
     &                   rc=rc)
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
     &                   yy=iyear,                                      &
     &                   mm=imonth,                                     &
     &                   dd=iday,                                       &
     &                   h=ihour,                                       &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get list of export fields.
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet (MODELS(Iwave)%ExportState(ng),                &
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
      CALL ESMF_StateGet (MODELS(Iwave)%ExportState(ng),                &
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
!   Get field from export field.
!
        CALL ESMF_StateGet (MODELS(Iwave)%ExportState(ng),              &
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
!  Get field pointer.  Usually, the DO-loop is executed once since
!  localDEcount=1.
!
        DE_LOOP : DO localDE=0,localDEcount-1
          CALL ESMF_FieldGet (field,                                    &
     &                        localDE=localDE,                          &
     &                        farrayPtr=ptr,                            &
     &                        rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
!  Initialize pointer to missing value.
!
          ptr=MISSING_r8
!
!  Load field data into export state.
!
          imin=LBOUND(ptr, DIM=1)
          imax=UBOUND(ptr, DIM=1)
          jmin=LBOUND(ptr, DIM=2)
          jmax=UBOUND(ptr, DIM=2)
!
          SELECT CASE (TRIM(ADJUSTL(ExportNameList(ifld))))
!
!  Surface roughness length scale (m).
!
            CASE ('zo', 'Zo')
              CALL WAM_Unpack (vm, ptr, imin, imax, jmin, jmax,         &
     &                         z0, rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  Friction velocity (m s-1).
!
            CASE ('ustar')
              CALL WAM_Unpack (vm, ptr, imin, imax, jmin, jmax,         &
     &                         ustar, rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!
            CASE ('tauw')
              CALL WAM_Unpack (vm, ptr, imin, imax, jmin, jmax,         &
     &                         tauw, rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
          END SELECT
!
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!
          IF (associated(ptr)) nullify (ptr)
        END DO DE_LOOP
!
!  Debugging: write out field into a netCDF format
!
        IF ((DebugLevel.ge.3).and.                                      &
     &      MODELS(Iwave)%ExportField(ifld)%debug_write) THEN
          WRITE (ofile,10) ng, TRIM(ExportNameList(ifld)),              &
     &                      iyear, imonth, iday, ihour
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
!  Deallocate array.
!
      IF (allocated(ExportNameList)) deallocate (ExportNameList)
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exitinh  WAM_Export',              &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT ('wam_',i2.2,'_export_',a,'_',i4.4,3('-',i2.2),'.nc')

      RETURN
      END SUBROUTINE WAM_Export
!
      SUBROUTINE WAM_Unpack (vm, ptr, imin, imax, jmin, jmax, var, rc)
!
!=======================================================================
!                                                                      !
!  Unpacks WAM component export variable and load it into pointer      !
!  after collecting data from all MPI nodes.                           !
!                                                                      !
!=======================================================================
!
      USE wam_mpi_module,  ONLY : nstart, nend, pelocal, irank
      USE wam_grid_module, ONLY : nx, ny, nsea, l_s_mask
!
!  Imported variable declarations.
!
      integer, intent(in) :: imin, imax, jmin, jmax
      integer, intent(inout) :: rc
!
      real (r4), intent(in) :: var(1:nsea)
      real (r8), intent(inout) :: ptr(imin:imax,jmin:jmax)
!
      TYPE (ESMF_VM), intent(in) :: vm
!
!  Local variable declarations.
!
      integer :: localPET, PETcount
      integer :: i, j, k, ii, jj, ijs, ijl
!
      integer (i4b), allocatable :: offsets_recv(:)
      integer (i4b), allocatable :: blocksize(:)
!
      real (r4), allocatable :: work1(:), work2(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WAM_Unpack"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WAM_Unpack',              &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the MPI
!  current node rank and number of nodes.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMGet (vm,                                              &
     &                 localPet=localPET,                               &
     &                 petCount=PETcount,                               &
     &                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Allocate local work arrays.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(work1)) THEN
        allocate ( work1(1:nsea) )
        work1=0.0_r4
      end if
      IF (.not. allocated(work2)) THEN
        allocate ( work2(nx,ny) )
        work2=0.0_r4
      END IF
!
      IF (.not. allocated(blocksize)) THEN
        allocate ( blocksize(PETcount) )
        blocksize=0_i4b
      END IF
      IF (.not.allocated(offsets_recv)) THEN
        allocate ( offsets_recv(PETcount) )
        offsets_recv=0_i4b
      END IF
!
!-----------------------------------------------------------------------
!  Collect block size from each PET.
!-----------------------------------------------------------------------
!
      ijs=nstart(irank)
      ijl=nend(irank)
!
      CALL ESMF_VMAllGatherV (vm,
     &                        sendData=(/ ijl-ijs+1 /),                 &
     &                        sendCount=1,                              &
     &                        recvData=blocksize,                       &
     &                        recvCounts=(/ (1, k=0,petCount-1) /),     &
     &                        recvOffsets=(/ (k, k=0,petCount-1) /),    &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_VMAllGatherV (vm,
     &                        sendData= (/ ijs-1 /),                    &
     &                        sendCount=1,                              &
     &                        recvData=offsets_recv,                    &
     &                        recvCounts=(/ (1, k=0,petCount-1) /),     &
     &                        recvOffsets=(/ (k, k=0,petCount-1) /),    &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from each PET.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMAllGatherV (vm,                                       &
     &                        sendData=var(1:blocksize(irank)),         &
     &                        sendCount=blocksize(irank),               &
     &                        recvData=work1,                           &
     &                        recvCounts=blocksize,                     &
     &                        recvOffsets=offsets_recv,                 &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Unpack data and fill pointer.
!-----------------------------------------------------------------------
!
      work2=UNPACK(work1, l_s_mask, MISSING_r4)
!
      DO j =jmin,jmax
        DO i=imin,imax
          IF (work2(i,j).lt.TOL_r4) THEN
            ptr(i,j)=work2(i,j)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Deallocate local work arrays.
!-----------------------------------------------------------------------
!
      IF (allocated(work1)) deallocate (work1)
      IF (allocated(work2)) deallocate (work2)
      IF (allocated(blocksize)) deallocate (blocksize)
      IF (allocated(offsets_recv)) deallocate (offsets_recv)
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WAM_Unpack',              &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WAM_Unpack
!
#endif
      END MODULE esmf_wam_mod
