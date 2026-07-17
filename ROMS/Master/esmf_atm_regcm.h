      MODULE esmf_regcm_mod

#if defined REGCM_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_atm_regcm.h 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module sets RegCM as the atmospheric model gridded component   !
!  using the generic ESMF/NUOPC layer:                                 !
!                                                                      !
!    RegCM_SetServices       Sets ATM component shared-object entry    !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                                                                      !
!    RegCM_SetInitializeP1   RegCM component phase 1 initialization:   !
!                            sets import and export fields long and    !
!                            short names into its respective state.    !
!                                                                      !
!    RegCM_SetInitializeP2   RegCM component phase 2 initialization:   !
!                            Initializes component (RegCM_Initialize), !
!                            sets component grid (RegCM_SetGridArrays),!
!                            and adds fields into import and export    !
!                            into respective states (RegCM_SetStates). !
!                                                                      !
!    RegCM_DataInit          Exports RegCM component fields during     !
!                            initialization or restart.                !
!                                                                      !
!    RegCM_SetClock          Sets RegCM component date calendar, start !
!                            and stop times, and coupling interval.    !
# ifdef ESM_SETRUNCLOCK
!                                                                      !
!    RegCM_SetRunClock       Sets RegCM run clock manually.            !
# endif
!                                                                      !
!    RegCM_CheckImport       Checks if RegCM component import field is !
!                            at the correct time.                      !
!                                                                      !
!    RegCM_SetGridArrays     Sets RegCM component horizontal grid      !
!                            arrays, grid area, and land/sea mask.     !
!                                                                      !
!    RegCM_SetStates         Adds RegCM component export and import    !
!                            fields into its respective state.         !
!                                                                      !
!    RegCM_ModelAdvance      Advances RegCM component for a coupling   !
!                            interval. It calls import and export      !
!                            routines.                                 !
!                                                                      !
!    RegCM_SetFinalize       Finalizes RegCM component execution.      !
!                                                                      !
!    RegCM_Import            Imports fields into RegCM from other      !
!                            gridded components.                       !
!                                                                      !
!    RegCM_Export            Exports RegCM fields to other gridded     !
!                            components.                               !
!                                                                      !
!    RegCM_uvrot             Rotates RegCM vector components to true   !
!                            EAST and true NORTH.                      !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!  RegCM:  ICTP Regional Climate Model (Version 4.5)                   !
!            http://gforge.ictp.it/gf/project/regcm                    !
!                                                                      !
!=======================================================================
!
      USE ESMF
      USE NUOPC
      USE NUOPC_Model,                                                  &
     &    NUOPC_SetServices          => SetServices,                    &
     &    NUOPC_Label_Advance        => label_Advance,                  &
     &    NUOPC_Label_DataInitialize => label_DataInitialize,           &
# ifdef ESM_SETRUNCLOCK
     &    NUOPC_Label_SetRunClock    => label_SetRunClock,              &
# endif
     &    NUOPC_Label_SetClock       => label_SetClock,                 &
     &    NUOPC_Label_CheckImport    => label_CheckImport
!
      USE mod_esmf_esm          ! ESM coupling structures and variables
!
      USE mod_regcm_interface, ONLY :                                   &
     &    RegCM_Initialize => RCM_initialize,                           &
     &    RegCM_Run        => RCM_run,                                  &
     &    RegCM_Finalize   => RCM_finalize
!
      implicit none
!
      PUBLIC  :: ATM_SetServices

      PRIVATE :: RegCM_SetInitializeP1
      PRIVATE :: RegCM_SetInitializeP2
      PRIVATE :: RegCM_DataInit
      PRIVATE :: RegCM_SetClock
# ifdef ESM_SETRUNCLOCK
      PRIVATE :: RegCM_SetRunClock
# endif
      PRIVATE :: RegCM_CheckImport
      PRIVATE :: RegCM_SetGridArrays
      PRIVATE :: RegCM_SetStates
      PRIVATE :: RegCM_ModelAdvance
      PRIVATE :: RegCM_SetFinalize
      PRIVATE :: RegCM_Import
      PRIVATE :: RegCM_Export
      PRIVATE :: RegCM_uvrot
!
      CONTAINS
!
      SUBROUTINE ATM_SetServices (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets RegCM component shared-object entry points for "initialize",   !
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
     &  __FILE__//", ATM_SetServices"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ATM_SetServices',         &
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
     &                              userRoutine=RegCM_InitializeP1,     &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for Phase 2 initialization.
!
      CALL NUOPC_CompSetEntryPoint (model,                              &
     &                              methodflag=ESMF_METHOD_INITIALIZE,  &
     &                              phaseLabelList=(/"IPDv00p2"/),      &
     &                              userRoutine=RegCM_InitializeP2,     &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Attach RegCM component phase independent specializing methods.
!-----------------------------------------------------------------------
!
!  Set routine for export initial/restart fields.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_DataInitialize,  &
     &                           specRoutine=RegCM_DataInit,            &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for setting RegCM clock.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_SetClock,        &
     &                           specRoutine=RegCM_SetClock,            &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF

# ifdef ESM_SETRUNCLOCK
!
!  Set routine for setting COAMPS run clock manually. First, remove the
!  default.
!
      CALL ESMF_MethodRemove (model,                                    &
     &                        NUOPC_label_SetRunClock,                  &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_SetRunClock,     &
     &                           specRoutine=RegCM_SetRunClock,         &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# endif
!
!  Set routine for checking import state.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_CheckImport,     &
     &                           specPhaseLabel="RunPhase1",            &
     &                           specRoutine=RegCM_CheckImport,         &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for time-stepping RegCM component.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_Advance,         &
     &                           specRoutine=RegCM_ModelAdvance,        &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Register RegCM finalize routine.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompSetEntryPoint (model,                           &
     &                                 methodflag=ESMF_METHOD_FINALIZE, &
     &                                 userRoutine=RegCM_SetFinalize,   &
     &                                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_SetServices',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE ATM_SetServices
!
      SUBROUTINE RegCM_SetInitializeP1 (model,                          &
     &                                  ImportState, ExportState,       &
     &                                  clock, rc)
!
!=======================================================================
!                                                                      !
!  RegCM component Phase 1 initialization: sets import and export      !
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
     &  __FILE__//", RegCM_SetInitializeP1"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_SetInitializeP1',   &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set RegCM import state and fields.
!-----------------------------------------------------------------------
!
      IMPORTING : IF (Nimport(Iatmos).gt.0) THEN
        DO ng=1,MODELS(Iatmos)%Ngrids
          IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
            CoupledSet=TRIM(COUPLED(Iatmos)%SetLabel(ng))
            StateLabel=TRIM(COUPLED(Iatmos)%ImpLabel(ng))
            CALL NUOPC_AddNestedState (ImportState,                     &
     &                                 CplSet=TRIM(CoupledSet),         &
     &                                 nestedStateName=TRIM(StateLabel),&
     &                                 nestedState=MODELS(Iatmos)%      &
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
            DO i=1,Nimport(Iatmos)
              StandardName=MODELS(Iatmos)%ImportField(i)%standard_name
              ShortName   =MODELS(Iatmos)%ImportField(i)%short_name
              CALL NUOPC_Advertise (MODELS(Iatmos)%ImportState(ng),     &
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
!  Set RegCM export state and fields.
!-----------------------------------------------------------------------
!
      EXPORTING : IF (Nexport(Iatmos).gt.0) THEN
        DO ng=1,MODELS(Iatmos)%Ngrids
          IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
            CoupledSet=TRIM(COUPLED(Iatmos)%SetLabel(ng))
            StateLabel=TRIM(COUPLED(Iatmos)%ExpLabel(ng))
            CALL NUOPC_AddNestedState (ExportState,                     &
     &                                 CplSet=TRIM(CoupledSet),         &
     &                                 nestedStateName=TRIM(StateLabel),&
     &                                 nestedState=MODELS(Iatmos)%      &
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
            DO i=1,Nexport(Iatmos)
              StandardName=MODELS(Iatmos)%ExportField(i)%standard_name
              ShortName   =MODELS(Iatmos)%ExportField(i)%short_name
              CALL NUOPC_Advertise (MODELS(Iatmos)%ExportState(ng),     &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_SetInitializeP1',   &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE RegCM_SetInitializeP1
!
      SUBROUTINE RegCM_SetInitializeP2 (model,                          &
     &                                  ImportState, ExportState,       &
     &                                  clock, rc)
!
!=======================================================================
!                                                                      !
!  RegCM component Phase 2 initialization: Initializes RegCM, sets     !
!  component grid, and adds import and export fields to respective     !
!  states.                                                             !
!                                                                      !
!=======================================================================
!
      USE mod_runparams, ONLY : dtsrf
      USE mod_update,    ONLY : RegCM_Allocate => RCM_Allocate
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
      integer :: ng
      integer :: MyComm, localPET, PETcount
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_SetInitializeP2"
!
      TYPE (ESMF_Time) :: CurrentTime, StartTime
      TYPE (ESMF_VM)   :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_SetInitializeP2',   &
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
!  Initialize RegCM component.
!-----------------------------------------------------------------------
!
      CALL RegCM_Initialize (mpiCommunicator=MyComm)
      CALL RegCM_Allocate ()
!
!-----------------------------------------------------------------------
!  Set-up grid and load coordinate data.
!-----------------------------------------------------------------------
!
      DO ng=1,MODELS(Iatmos)%Ngrids
        IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
          CALL RegCM_SetGridArrays (ng, model, rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Set-up fields and register to import/export states.
!-----------------------------------------------------------------------
!
      DO ng=1,MODELS(Iatmos)%Ngrids
        IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
          CALL RegCM_SetStates (ng, model, rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
        END IF
      END DO
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_SetInitializeP2',   &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE RegCM_SetInitializeP2
!
      SUBROUTINE RegCM_DataInit (model, rc)
!
!=======================================================================
!                                                                      !
!  Exports RegCM component fields during initialization or restart.    !
!                                                                      !
!=======================================================================
!
      USE mod_runparams, ONLY : dtsrf
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: is, ng
      integer :: localPET, PETcount, phase
!
      real (dp) :: tstr
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_DataInit"

      character (ESMF_MAXSTR) :: str1, str2
!
      TYPE (ESMF_Clock)        :: clock
      TYPE (ESMF_Time)         :: CurrentTime
      TYPE (ESMF_TimeInterval) :: TimeStep
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_DataInit',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get gridded component clock.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       clock=clock,                               &
     &                       localPet=localPET,                         &
     &                       petCount=PETcount,                         &
     &                       currentPhase=phase,                        &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_ClockGet (clock,                                        &
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
!  If initialization or restart, put export fields.
!-----------------------------------------------------------------------
!
      IF (Restarted.and.(CurrentTime.eq.ESM_RestartTime)) THEN
!
!  Debugging: report time information.
!
        IF ((DebugLevel.ge.0).and.(localPET.eq.0)) THEN
          CALL ESMF_TimeGet (CurrentTime,                               &
     &                       timeStringISOFrac=str1,                    &
     &                       rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          is=INDEX(str1, 'T')                        ! remove 'T' in
          IF (is.gt.0) str1(is:is)=' '               ! ISO 8601 format
!
          CALL ESMF_TimeGet (CurrentTime+TimeStep,                      &
     &                       timeStringISOFrac=str2,                    &
     &                       rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          is=INDEX(str2, 'T')                        ! remove 'T' in
          IF (is.gt.0) str2(is:is)=' '               ! ISO 8601 format
!
          IF (DebugLevel.eq.0) THEN
            WRITE (cplout,10) TRIM(str1), TRIM(str2), phase
          ELSE
            WRITE (cplout,20) TRIM(str1), TRIM(str2), phase, 0.0_r8,    &
     &                        dtsrf
          END IF
        END IF
!
!  Run RegCM component only for one time-step to fill variables.
!
        CALL RegCM_Run (0.0_r8, dtsrf)
!
!  Put export fields.
!
        IF (Nexport(Iatmos).gt.0) THEN
          DO ng=1,MODELS(Iatmos)%Ngrids
            CALL RegCM_Export (ng, model, rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END DO
        END IF
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_DataInit',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (/,' ESMF, Running RegCM: ',a,' --> ',a,' Phase: ',i1)
  20  FORMAT (/,' ESMF, Running RegCM: ',a,' --> ',a,' Phase: ',i1,     &
     &        ' [',f15.2, '-', f15.2, ' s]')

      RETURN
      END SUBROUTINE RegCM_DataInit
!
      SUBROUTINE RegCM_SetClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets RegCM component date calendar, start and stop time, and        !
!  coupling interval.                                                  !
!                                                                      !
!=======================================================================
!
      USE mod_dynparam , ONLY : calendar
      USE mod_runparams, ONLY : idate0, idate1, idate2, dtsec
      USE mod_date,      ONLY : split_idate
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
      integer :: ref_year,   start_year,   stop_year
      integer :: ref_month,  start_month,  stop_month
      integer :: ref_day,    start_day,    stop_day
      integer :: ref_hour,   start_hour,   stop_hour
      integer :: ref_minute, start_minute, stop_minute
      integer :: ref_second, start_second, stop_second
      integer :: TimeFrac
!
      character (len= 20) :: Calendar
      character (len=160) :: message

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_SetClock"
!
      TYPE (ESMF_CalKind_Flag) :: CalType
      TYPE (ESMF_Time)         :: StartTime
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_SetClock',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Create RegCM component clock.
!-----------------------------------------------------------------------
!
      Calendar=TRIM(ClockInfo(Iatmos)%CalendarString)
      IF (TRIM(Calendar).eq.'gregorian') THEN
        CalType=ESMF_CALKIND_GREGORIAN
      ELSE IF ((Calendar.eq.'noleap').or.(Calendar.eq.'365_day')) THEN
        CalType=ESMF_CALKIND_NOLEAP
      ELSE IF (Calendar.eq.'360_day') THEN
        CalType=ESMF_CALKIND_360DAY
      ELSE
        CalType=ESMF_CALKIND_GREGORIAN
      END IF
!
      ClockInfo(Iatmos)%Calendar=ESMF_CalendarCreate(CalType,           &
     &                                             name=TRIM(Calendar), &
     &                                               rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set reference time.
!
      CALL split_idate (idate0, ref_year, ref_month, ref_day, ref_hour)
      ref_minute=0
      ref_second=0
!
      CALL ESMF_TimeSet (ClockInfo(Iatmos)%ReferenceTime,               &
     &                   yy=ref_year,                                   &
     &                   mm=ref_month,                                  &
     &                   dd=ref_day,                                    &
     &                   h =ref_hour,                                   &
     &                   m =ref_minute,                                 &
     &                   s =ref_second,                                 &
     &                   calendar=ClockInfo(Iatmos)%Calendar,           &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set start time.
!
      CALL split_idate (idate1, start_year, start_month, start_day,     &
     &                  start_hour)
      start_minute=0
      start_second=0
!
      CALL ESMF_TimeSet (ClockInfo(Iatmos)%StartTime,                   &
     &                   yy=start_year,                                 &
     &                   mm=start_month,                                &
     &                   dd=start_day,                                  &
     &                   h =start_hour,                                 &
     &                   m =start_minute,                               &
     &                   s =start_second,                               &
     &                   calendar=ClockInfo(Iatmos)%Calendar,           &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set stop time.
!
      CALL split_idate (idate2, stop_year, stop_month, stop_day,        &
     &                  stop_hour)
      stop_minute=0
      stop_second=0
!
      CALL ESMF_TimeSet (ClockInfo(Iatmos)%StopTime,                    &
     &                   yy=stop_year,                                  &
     &                   mm=stop_month,                                 &
     &                   dd=stop_day,                                   &
     &                   h =stop_hour,                                  &
     &                   m =stop_minute,                                &
     &                   s =stop_second,                                &
     &                   calendar=ClockInfo(Iatmos)%Calendar,           &
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
     &                       clock=ClockInfo(Iatmos)%Clock,             &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_ClockGet (ClockInfo(Iatmos)%Clock,                      &
     &                    timeStep=ClockInfo(Iatmos)%TimeStep,          &
     &                    currTime=ClockInfo(Iatmos)%CurrentTime,       &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Compare driver time against RegCM component time.
!-----------------------------------------------------------------------
!
      IF (ClockInfo(Idriver)%Restarted) THEN
        StartTime=ClockInfo(Idriver)%RestartTime
      ELSE
        StartTime=ClockInfo(Idriver)%StartTime
      END IF
!
      IF (ClockInfo(Iatmos)%StartTime.ne.StartTime) THEN
        CALL ESMF_TimePrint (ClockInfo(Iatmos)%StartTime,               &
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
        message='Driver and RegCM start times do not match: '//         &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ClockInfo(Iatmos )%StopTime.ne.                               &
     &    ClockInfo(Idriver)%StopTime) THEN
        CALL ESMF_TimePrint (ClockInfo(Iatmos)%StopTime,                &
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
        message='Driver and RegCM stop times do not match: '//          &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ClockInfo(Iatmos )%Calendar.ne.                               &
     &    ClockInfo(Idriver)%Calendar) THEN
        CALL ESMF_CalendarPrint (ClockInfo(Iatmos)%Calendar,            &
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
        message='Driver and RegCM calendars do not match: '//           &
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
      DO ng=1,MODELS(Iatmos)%Ngrids
        IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
          TimeFrac=MAX(TimeFrac,                                        &
     &                 MAXVAL(MODELS(Iatmos)%TimeFrac(ng,:),            &
     &                        mask=MODELS(:)%IsActive))
        END IF
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
      ClockInfo(Iatmos)%TimeStep=ClockInfo(Idriver)%TimeStep/TimeFrac
!
      ClockInfo(Iatmos)%Name='RegCM_clock'
      CALL ESMF_ClockSet (ClockInfo(Iatmos)%Clock,                      &
     &                    name=TRIM(ClockInfo(Iatmos)%Name),            &
     &                    refTime  =ClockInfo(Iatmos)%ReferenceTime,    &
     &                    timeStep =ClockInfo(Iatmos)%TimeStep,         &
     &                    startTime=ClockInfo(Iatmos)%StartTime,        &
     &                    stopTime =ClockInfo(Iatmos)%StopTime,         &
                          rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_SetClock',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE RegCM_SetClock

# ifdef ESM_SETRUNCLOCK
!
      SUBROUTINE RegCM_SetRunClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets RegCM run clock manually to avoid getting zero time stamps at  !
!  the first regridding call.                                          !
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
     &  __FILE__//", COAMPS_SetRunClock"
!
      TYPE (ESMF_Clock) :: driverClock, modelClock
      TYPE (ESMF_Time)  :: currTime
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_SetRunClock',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set ROMS run clock manually.
!-----------------------------------------------------------------------
!
!  Inquire driver and model clock.
!
      CALL NUOPC_ModelGet (model,                                       &
     &                     driverClock=driverClock,                     &
     &                     modelClock=modelClock,                       &
     &                     rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set model clock to have the current start time as the driver clock.
!
      CALL ESMF_ClockGet (driverClock,                                  &
     &                    currTime=currTime,                            &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_ClockSet (modelClock,                                   &
     &                    currTime=currTime,                            &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Check and set the component clock against the driver clock.
!
      CALL NUOPC_CompCheckSetClock (model,                              &
     &                              driverClock,                        &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_SetRunClock',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE RegCM_SetRunClock
# endif
!
      SUBROUTINE RegCM_CheckImport (model, rc)
!
!=======================================================================
!                                                                      !
!  Checks if RegCM component import field is at the correct time.      !
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
      logical :: IsValid, atCorrectTime
!
      integer :: ImportCount, i, is, localPET, ng
!
      real (dp) :: TcurrentInSeconds
!
      character (len=22) :: DriverTimeString, FieldTimeString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_CheckImport"
!
      character (ESMF_MAXSTR) :: string, FieldName
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)
!
      TYPE (ESMF_Clock)        :: DriverClock
      TYPE (ESMF_Field)        :: field
      TYPE (ESMF_Time)         :: StartTime, CurrentTime
      TYPE (ESMF_Time)         :: DriverTime, FieldTime
      TYPE (ESMF_TimeInterval) :: TimeStep
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_CheckImport',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Query component.
!-----------------------------------------------------------------------
!
      CALL NUOPC_ModelGet (model,                                       &
     &                     driverClock=DriverClock,                     &
     &                     rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
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
!  Get the start time and current time from driver clock.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (DriverClock,                                  &
     &                    timeStep=TimeStep,                            &
     &                    startTime=StartTime,                          &
     &                    currTime=DriverTime,                          &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeGet (DriverTime,                                    &
     &                   s_r8=TcurrentInSeconds,                        &
     &                   timeStringISOFrac=DriverTimeString,            &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(DriverTimeString, 'T')                 ! remove 'T' in
      IF (is.gt.0) DriverTimeString(is:is)=' '        ! ISO 8601 format
!
!-----------------------------------------------------------------------
!  Get list of import fields.
!-----------------------------------------------------------------------
!
      IF (Nimport(Iatmos).gt.0) THEN
        NESTED_LOOP : DO ng=1,MODELS(Iatmos)%Ngrids
          IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
            CALL ESMF_StateGet (MODELS(Iatmos)%ImportState(ng),         &
     &                          itemCount=ImportCount,                  &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            IF (.not.allocated(ImportNameList)) THEN
              allocate ( ImportNameList(ImportCount) )
            END IF
!
            CALL ESMF_StateGet (MODELS(Iatmos)%ImportState(ng),         &
     &                          itemNameList=ImportNameList,            &
     &                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!-----------------------------------------------------------------------
!  Only check fields in the ImportState object.
!-----------------------------------------------------------------------
!
            FIELD_LOOP : DO i=1,ImportCount
              FieldName=TRIM(ImportNameList(i))
              CALL ESMF_StateGet (MODELS(Iatmos)%ImportState(ng),       &
     &                            itemName=TRIM(FieldName),             &
     &                            field=field,                          &
     &                            rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  If debugging, report field timestamp.
!
              IF (DebugLevel.gt.1) THEN
                CALL NUOPC_GetTimeStamp (field,                         &
     &                                   isValid = IsValid,             &
     &                                   time = FieldTime,              &
     &                                   rc = rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
                IF (IsValid) THEN
                  CALL ESMF_TimeGet (FieldTime,                         &
     &                             timeStringISOFrac = FieldTimeString, &
     &                               rc=rc)
                  IF (ESMF_LogFoundError(rcToCheck=rc,                  &
     &                                   msg=ESMF_LOGERR_PASSTHRU,      &
     &                                   line=__LINE__,                 &
     &                                   file=MyFile)) THEN
                    RETURN
                  END IF
                  is=INDEX(FieldTimeString, 'T')            ! remove 'T'
                  IF (is.gt.0) FieldTimeString(is:is)=' '
!
                  IF (localPET.eq.0) THEN
                    WRITE (cplout,10) TRIM(FieldName),                  &
     &                                TRIM(FieldTimeString),            &
     &                                TRIM(DriverTimeString)
                  END IF
                END IF
              END IF
!
!  Check if import field is at the correct time.
!
              string='COAMPS_CheckImport - '//TRIM(FieldName)//' field'
              CurrentTime=DriverTime
!
              atCorrectTime=NUOPC_IsAtTime(field,                       &
     &                                     CurrentTime,                 &
     &                                     rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
              IF (.not.atCorrectTime) THEN
                CALL report_timestamp (field, CurrentTime,              &
     &                                 localPET, TRIM(string), rc)
!
                string='NUOPC INCOMPATIBILITY DETECTED: Import '//      &
     &                 'Fields not at correct time'
                CALL ESMF_LogSetError(ESMF_RC_NOT_VALID,                &
     &                                msg=TRIM(string),                 &
     &                                line=__LINE__,                    &
     &                                file=MyFile,                      &
     &                                rcToReturn=rc)
                RETURN
              END IF
            END DO FIELD_LOOP
            IF (allocated(ImportNameList)) deallocate (ImportNameList)
          END IF
        END DO NESTED_LOOP
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_CheckImport',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (1x,'RegCM_CheckImport - ',a,':',t32,'TimeStamp = ',a,     &
     &        ',  DriverTime = ',a)
!
      RETURN
      END SUBROUTINE RegCM_CheckImport
!
      SUBROUTINE RegCM_SetGridArrays (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Sets RegCM component staggered, horizontal grids arrays, grid area, !
!  and land/sea mask.                                                  !
!                                                                      !
!=======================================================================
!
      USE mod_mppparam,      ONLY : ma
      USE mod_runparams,     ONLY : dxsq
      USE mod_atm_interface, ONLY : mddom
      USE mod_dynparam,      ONLY : iy, jx, nproc,                      &
     &                              ide1, ide2, jde1, jde2,             &
     &                              idi1, idi2, jdi1, jdi2,             &
     &                              ice1, ice2, jce1, jce2
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp), intent(inout) :: model
!
!  Local variable declarations.
!
      integer :: gtype, i, ivar, j
      integer :: localDE, localDEcount, localPET, PETcount
      integer :: cpus_per_dim(2)
!
      integer (i4b), pointer :: ptrM(:,:) => NULL()
!
      real (dp), pointer :: ptrA(:,:) => NULL()
      real (dp), pointer :: ptrX(:,:) => NULL()
      real (dp), pointer :: ptrY(:,:) => NULL()
!
      character (len=40) :: name

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_SetGridArrays"
!
      TYPE (ESMF_DistGrid)   :: distGrid
      TYPE (ESMF_StaggerLoc) :: staggerLoc
      TYPE (ESMF_VM)         :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_SetGridArrays',     &
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
!  Calculate number of CPUs in each direction.
!-----------------------------------------------------------------------
!
      IF (nproc.lt.4) THEN
        cpus_per_dim(2)=1
        cpus_per_dim(1)=nproc
      ELSE IF (nproc.ge.4) THEN
        cpus_per_dim(2)=(NINT(SQRT(DBLE(nproc)))/2)*2
        IF (iy.gt.INT(1.5*DBLE(jx))) THEN
          cpus_per_dim(2)=cpus_per_dim(2)-1
          DO WHILE (MOD(nproc,cpus_per_dim(2)).ne.0)
            cpus_per_dim(2)=cpus_per_dim(2)-1
          end do
        ELSE IF (jx.gt.INT(1.5*dble(iy))) THEN
          cpus_per_dim(2)=cpus_per_dim(2)+1
          DO WHILE (MOD(nproc,cpus_per_dim(2)).ne.0)
            cpus_per_dim(2)=cpus_per_dim(2)+1
          END DO
        ELSE
          DO WHILE (MOD(nproc,cpus_per_dim(2)).ne.0)
            cpus_per_dim(2)=cpus_per_dim(2)+1
          END DO
        END IF
        cpus_per_dim(1)=nproc/cpus_per_dim(2)
      END IF
!
!-----------------------------------------------------------------------
!  Create DistGrid based on model domain decomposition.
!
!  ESMF is basically using a right handed coordinate system, and
!  using the Fortran way of using the smallest stride to the first
!  dimension but RegCM not. The order of dimension is reversed
!  because of this limitation.
!-----------------------------------------------------------------------
!
      distGrid=ESMF_DistGridCreate(minIndex=(/ 1, 1 /),                 &
     &                             maxIndex=(/ iy, jx /),               &
     &                             regDecomp=cpus_per_dim,              &
     &                             rc=rc)
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
!  Define component grid location type: Arakawa B-grid.
!
      IF (.not.allocated(MODELS(Iatmos)%mesh)) THEN
        allocate ( MODELS(Iatmos)%mesh(2) )
        MODELS(Iatmos)%mesh(1)%gtype=Icenter
        MODELS(Iatmos)%mesh(2)%gtype=Icorner
      END IF
!
!  Create ESMF Grid.
!
      MODELS(Iatmos)%grid(ng)=ESMF_GridCreate(distgrid=distGrid,        &
     &                                  indexflag=ESMF_INDEX_GLOBAL,    &
     &                                  name=TRIM(MODELS(Iatmos)%name), &
     &                                        rc=rc)
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
      CALL ESMF_GridGet (MODELS(Iatmos)%grid(ng),                       &
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
      MESH_LOOP : DO ivar=1,UBOUND(MODELS(Iatmos)%mesh, DIM=1)
!
!  Set staggering type, Arakawa B-grid.
!
        SELECT CASE (MODELS(Iatmos)%mesh(ivar)%gtype)
          CASE (Icenter)
            staggerLoc=ESMF_STAGGERLOC_CENTER
          CASE (Icorner)
            staggerLoc=ESMF_STAGGERLOC_CORNER
        END SELECT
!
!  Allocate coordinate storage associated with staggered grid type.
!  No coordinate values are set yet.
!
        CALL ESMF_GridAddCoord (MODELS(Iatmos)%grid(ng),                &
     &                          staggerLoc=staggerLoc,                  &
     &                          rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Allocate storage for masking.
!
        CALL ESMF_GridAddItem (MODELS(Iatmos)%grid(ng),                 &
     &                         staggerLoc=staggerLoc,                   &
     &                         itemflag=ESMF_GRIDITEM_MASK,             &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  The RegCM masking is as follows,  0: ocean
!                                    2: land
!
        MODELS(Iatmos)%LandValue=2
        MODELS(Iatmos)%SeaValue=0
!
!  Allocate storage for grid area.
!
        CALL ESMF_GridAddItem (MODELS(Iatmos)%grid(ng),                 &
     &                         staggerLoc=staggerLoc,                   &
     &                         itemflag=ESMF_GRIDITEM_AREA,             &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Get pointers and set coordinates for the grid.  Usually, the DO-loop
!  is executed once since localDEcount=1.
!
        DE_LOOP : DO localDE=0,localDEcount-1
          CALL ESMF_GridGetCoord (MODELS(Iatmos)%grid(ng),              &
     &                            coordDim=1,                           &
     &                            staggerLoc=staggerLoc,                &
     &                            localDE=localDE,                      &
     &                            farrayPtr=ptrX,                       &
     &                            rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_GridGetCoord (MODELS(Iatmos)%grid(ng),              &
     &                            coordDim=2,                           &
     &                            staggerLoc=staggerLoc,                &
     &                            localDE=localDE,                      &
     &                            farrayPtr=ptrY,                       &
     &                            rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_GridGetItem (MODELS(Iatmos)%grid(ng),               &
     &                           itemflag=ESMF_GRIDITEM_MASK,           &
     &                           staggerLoc=staggerLoc,                 &
     &                           localDE=localDE,                       &
     &                           farrayPtr=ptrM,                        &
     &                           rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_GridGetItem (MODELS(Iatmos)%grid(ng),               &
     &                           itemflag=ESMF_GRIDITEM_AREA,           &
     &                           staggerLoc=staggerLoc,                 &
     &                           localDE=localDE,                       &
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
          SELECT CASE (MODELS(Iatmos)%mesh(ivar)%gtype)
            CASE (Icorner)
              DO i=ide1,ide2
                DO j=jde1,jde2
                  ptrX(i,j)=mddom%dlon(j,i)
                  ptrY(i,j)=mddom%dlat(j,i)
                END DO
              END DO
              ptrA=dxsq
!
              IF (ma%has_bdyright) THEN
                ptrX(:,jde2+1)=ptrX(:,jde2)+                            &
     &                         (ptrX(:,jde2)-ptrX(:,jde2-1))
                ptrY(:,jde2+1)=ptrY(:,jde2)+                            &
     &                         (ptrY(:,jde2)-ptrY(:,jde2-1))
              END IF
!
              IF (ma%has_bdytop) THEN
                ptrX(ide2+1,:)=ptrX(ide2,:)+                            &
     &                         (ptrX(ide2,:)-ptrX(ide2-1,:))
                ptrY(ide2+1,:)=ptrY(ide2,:)+                            &
     &                         (ptrY(ide2,:)-ptrY(ide2-1,:))
              END IF
            CASE (Icenter)
              DO i=ice1,ice2
                DO j=jce1,jce2
                  ptrX(i,j)=mddom%xlon(j,i)
                  ptrY(i,j)=mddom%xlat(j,i)
                  ptrM(i,j)=INT(mddom%mask(j,i))
                END DO
              END DO
              ptrA=dxsq
!
              IF (ma%has_bdyright) THEN
                ptrX(:,jce2+1)=ptrX(:,jce2)+                            &
     &                         (ptrX(:,jce2)-ptrX(:,jce2-1))
                ptrY(:,jce2+1)=ptrY(:,jce2)+                            &
     &                         (ptrY(:,jce2)-ptrY(:,jce2-1))
              END IF
!
              IF (ma%has_bdytop) THEN
                ptrX(ice2+1,:)=ptrX(ice2,:)+                            &
     &                         (ptrX(ice2,:)-ptrX(ice2-1,:))
                ptrY(ice2+1,:)=ptrY(ice2,:)+                            &
     &                         (ptrY(ice2,:)-ptrY(ice2-1,:))
              END IF
          END SELECT
!
!  Nullify pointers.
!
          IF ( associated(ptrX) ) nullify (ptrX)
          IF ( associated(ptrY) ) nullify (ptrY)
          IF ( associated(ptrM) ) nullify (ptrM)
          IF ( associated(ptrA) ) nullify (ptrA)
        END DO DE_LOOP
!
!  Debugging: write out component grid in VTK format.
!
        IF (DebugLevel.ge.4) THEN
          gtype=MODELS(Iatmos)%mesh(ivar)%gtype
          CALL ESMF_GridWriteVTK (MODELS(Iatmos)%grid(ng),              &
     &                            filename="regcm_"//                   &
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
     &                       grid=MODELS(Iatmos)%grid(ng),              &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_SetGridArrays',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE RegCM_SetGridArrays
!
      SUBROUTINE RegCM_SetStates (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Adds RegCM component export and import fields into its respective   !
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
      integer :: localPET, PETcount
      integer :: ExportCount, ImportCount
!
      real (dp), dimension(:,:), pointer :: ptr2d => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_SetStates"
!
      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)
!
      TYPE (ESMF_ArraySpec)  :: arraySpec2d
      TYPE (ESMF_Field)      :: field
      TYPE (ESMF_StaggerLoc) :: staggerLoc
      TYPE (ESMF_VM)         :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_SetStates',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get gridded component information.
!-----------------------------------------------------------------------
!
!  Get import and export states.
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
!  Get number of local decomposition elements (DEs). Usually, a single
!  Decomposition Element (DE) is associated with each Persistent
!  Execution Thread (PETs). Thus, localDEcount=1.
!
      CALL ESMF_GridGet (MODELS(Iatmos)%grid(ng),                       &
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
      CALL ESMF_ArraySpecSet (arraySpec2d,                              &
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
      EXPORTING : IF (Nexport(Iatmos).gt.0) THEN
!
!  Get number of fields to export.
!
        CALL ESMF_StateGet (MODELS(Iatmos)%ExportState(ng),             &
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
        CALL ESMF_StateGet (MODELS(Iatmos)%ExportState(ng),             &
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
          id=field_index(MODELS(Iatmos)%ExportField, ExportNameList(i))
!
          IF (NUOPC_IsConnected(MODELS(Iatmos)%ExportState(ng),         &
     &                          fieldName=TRIM(ExportNameList(i)),      &
     &                          rc=rc)) THEN
!
!  Set staggering type.
!
            SELECT CASE (MODELS(Iatmos)%ExportField(id)%gtype)
              CASE (Icenter)
                staggerLoc=ESMF_STAGGERLOC_CENTER
              CASE (Icorner)
                staggerLoc=ESMF_STAGGERLOC_CORNER
            END SELECT
!
!  Create 2D field from the Grid and arraySpec.
!
            field=ESMF_FieldCreate(MODELS(Iatmos)%grid(ng),             &
     &                             arraySpec2d,                         &
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
            CALL NUOPC_Realize (MODELS(Iatmos)%ExportState(ng),         &
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
     &                          TRIM(COUPLED(Iatmos)%ExpLabel(ng))
            END IF
            CALL ESMF_StateRemove (MODELS(Iatmos)%ExportState(ng),      &
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
      IMPORTING : IF (Nimport(Iatmos).gt.0) THEN
!
!  Get number of fields to import.
!
        CALL ESMF_StateGet (MODELS(Iatmos)%ImportState(ng),             &
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
        CALL ESMF_StateGet (MODELS(Iatmos)%ImportState(ng),             &
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
          id=field_index(MODELS(Iatmos)%ImportField, ImportNameList(i))
!
          IF (NUOPC_IsConnected(MODELS(Iatmos)%ImportState(ng),         &
     &                          fieldName=TRIM(ImportNameList(i)),      &
     &                          rc=rc)) THEN
!
!  Set staggering type.
!
            SELECT CASE (MODELS(Iatmos)%ImportField(id)%gtype)
              CASE (Icenter)
                staggerLoc=ESMF_STAGGERLOC_CENTER
              CASE (Icorner)
                staggerLoc=ESMF_STAGGERLOC_CORNER
            END SELECT
!
!  Create 2D field from the Grid, arraySpec.
!
            field=ESMF_FieldCreate(MODELS(Iatmos)%grid(ng),             &
     &                             arraySpec2d,                         &
     &                             staggerloc=staggerLoc,               &
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
            CALL NUOPC_Realize (MODELS(Iatmos)%ImportState(ng),         &
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
     &                          TRIM(COUPLED(Iatmos)%ImpLabel(ng))
            END IF
            CALL ESMF_StateRemove (MODELS(Iatmos)%ImportState(ng),      &
     &                             (/ TRIM(ImportNameList(i)) /),       &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_SetStates',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
!
 10   FORMAT ('RegCM_SetStates - Removing field ''',a,''' from ',a,     &
     &        '''',a,'''',/,18x,'because it is not connected.')
!
      RETURN
      END SUBROUTINE RegCM_SetStates
!
      SUBROUTINE RegCM_ModelAdvance (model, rc)
!
!=======================================================================
!                                                                      !
!  Advance RegCM component for a coupling interval (seconds) using     !
!  "RegCM_run". It also calls "RegCM_Import" and "RegCM_Export" to     !
!  import and export coupling fields, respectively.                    !
!                                                                      !
!=======================================================================
!
      USE mod_runparams, ONLY : ifrest, ktau, dtsrf
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc

      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: is, ng
      integer :: localPET, PETcount, phase
!
      real (dp) :: tstr, tend
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_SetModelAdvance"
!
      character (ESMF_MAXSTR) :: str1, str2
!
      TYPE (ESMF_Clock)        :: clock
      TYPE (ESMF_State)        :: ExportState, ImportState
      TYPE (ESMF_TimeInterval) :: TimeFrom, TimeTo, TimeStep
      TYPE (ESMF_Time)         :: ReferenceTime
      TYPE (ESMF_Time)         :: CurrentTime, StartTime, StopTime
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_ModelAdvance',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get information about the gridded component.
!-----------------------------------------------------------------------
!
!  Inquire about RegCM component.
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       importState=ImportState,                   &
     &                       exportState=ExportState,                   &
     &                       clock=clock,                               &
     &                       localPet=localPET,                         &
     &                       petCount=PETcount,                         &
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
!  Get driver time step interval, start time, stopping time, reference
!  time, and current time.
!
      CALL ESMF_ClockGet (clock,                                        &
     &                    timeStep=TimeStep,                            &
     &                    startTime=StartTime,                          &
     &                    stopTime=StopTime,                            &
     &                    refTime=ReferenceTime,                        &
     &                    currTime=ClockInfo(Iatmos)%CurrentTime,       &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Calculate run time.
!-----------------------------------------------------------------------
!
      IF (ClockInfo(Idriver)%Restarted) THEN
        TimeFrom=CurrentTime-ClockInfo(Iatmos)%RestartTime
      ELSE
        TimeFrom=CurrentTime-ClockInfo(Iatmos)%StartTime
      END IF
!
      CALL ESMF_TimeIntervalGet (TimeFrom,                              &
     &                           s_r8=tstr,                             &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      TimeTo=TimeFrom+TimeStep
      CALL ESMF_TimeIntervalGet (TimeTo,                                &
     &                           s_r8=tend,                             &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ClockInfo(Idriver)%Restarted.and.                             &
     &    (StartTime.eq.CurrentTime)) THEN
        tstr=tstr+dtsrf
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
          WRITE (cplout,20) TRIM(str1), TRIM(str2), phase, tstr, tend
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Get import fields.
!-----------------------------------------------------------------------
!
      IF ((CurrentTime.ne.RefTime).or.Restarted) THEN
        IF (Nimport(Iatmos).gt.0) THEN
          DO ng=1,MODELS(Iatmos)%Ngrids
            IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
              CALL RegCM_Import (ng, model, rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Run RegCM component. Notice that atmosphere component is advanced
!  when ng=1.  In nested application, its numerical kernel will advance
!  all the nested grids in their logical order.
!-----------------------------------------------------------------------
!
      CALL RegCM_Run (tstr, tend)
!
!-----------------------------------------------------------------------
!  Put export fields.
!-----------------------------------------------------------------------
!
      IF (Nexport(Iatmos).gt.0) THEN
        DO ng=1,MODELS(Iatmos)%Ngrids
          IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
            CALL RegCM_Export (ng, model, rc)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_ModelAdvance',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (' Running RegCM Component: ',a,' --> ',a,' Phase: ',i1)
  20  FORMAT (' Running RegCM Component: ',A,' --> ',A,' Phase: ',i1,   &
     &        ' [',f15.2, '-', f15.2, ']')

      RETURN
      END SUBROUTINE RegCM_ModelAdvance
!
      SUBROUTINE RegCM_SetFinalize (model,                              &
     &                              ImportState, ExportState,           &
     &                              clock, rc)
!
!=======================================================================
!                                                                      !
!  Finalize RegCM component execution. It calls RegCM_finalize.        !
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
     &  __FILE__//", RegCM_SetFinalize"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_SetFinalize',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Finalize RegCM component.
!-----------------------------------------------------------------------
!
      CALL RegCM_Finalize ()
      CALL my_flush (6)                   ! flush standard output buffer
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_SetFinalize',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE RegCM_SetFinalize
!
      SUBROUTINE RegCM_Import (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Imports fields into RegCM array structure from other coupled        !
!  gridded components.                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_update,   ONLY : importFields
      USE mod_dynparam, ONLY : ici1, ici2, jci1, jci2
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
      integer :: id, ifld, i, is, j
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: ImportCount
      integer :: localDE, localDEcount, localPET, PETcount
      integer :: LBi, UBi, LBj, UBj
      integer :: IminP, ImaxP, JminP, JmaxP
!
      real (dp) :: Fseconds, TimeInDays, Time_Current

      real (dp) :: MyFmax(2), MyFmin(2), Fmin(2), Fmax(2), Fval
      real (dp) :: scale, add_offset
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      character (len=22 )     :: Time_CurrentString
      character (len=100)     :: FieldName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_Import"

      character (ESMF_MAXSTR) :: cname, ofile
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)
!
      TYPE (ESMF_Clock) :: clock
      TYPE (ESMF_Field) :: field
      TYPE (ESMF_State) :: ImportState
      TYPE (ESMF_Time)  :: CurrentTime
      TYPE (ESMF_VM)    :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_Import',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Compute RegCM lower and upper bounds (non-overlapping) for physical
!  area per nested grid and tile.
!-----------------------------------------------------------------------
!
      IminP=ice1
      ImaxP=ice2
      JminP=jce1
      JmaxP=jce2
!
!-----------------------------------------------------------------------
!  Get information about the gridded component.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       importState=ImportState,                   &
     &                       clock=clock,                               &
     &                       localPet=localPET,                         &
     &                       petCount=PETcount,                         &
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
!  Get number of local decomposition elements (DEs). Usually, a single
!  DE is associated with each Persistent Execution Thread (PETs). Thus,
!  localDEcount=1.
!
      CALL ESMF_GridGet (MODELS(Iatmos)%grid(ng),                       &
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
      Fseconds=REAL(seconds,dp)+REAL(sN,dp)/REAL(sD,dp)
      TimeInDays=Time_Current/86400.0_dp
      is=INDEX(Time_CurrentString, 'T')              ! remove 'T' in
      IF (is.gt.0) Time_CurrentString(is:is)=' '     ! ISO 8601 format
!
!-----------------------------------------------------------------------
!  Get list of import fields.
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet (MODELS(Iatmos)%ImportState(ng),               &
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
      CALL ESMF_StateGet (MODELS(Iatmos)%ImportState(ng),               &
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
        id=field_index(MODELS(Iatmos)%ImportField, ImportNameList(ifld))
!
!  Get field from import state.
!
        CALL ESMF_StateGet (MODELS(Iatmos)%ImportState(ng),             &
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
!  Get field pointer.  Usually, the DO-loop is executed once since
!  localDEcount=1.
!
        DE_LOOP : DO localDE=0,localDEcount-1
          CALL ESMF_FieldGet (field,                                    &
     &                        localDE=localDE,                          &
     &                        farrayPtr=ptr2d,                          &
     &                        rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          LBi=LBOUND(ptr2d,1)
          UBi=UBOUND(ptr2d,1)
          LBj=LBOUND(ptr2d,2)
          UBj=UBOUND(ptr2d,2)
!
!  Initialize import field parameters.  Set "scale" and "add_offset"
!  values need to convert imported fields to RegCM requirements.
!
          scale     =MODELS(Iatmos)%ImportField(id)%scale_factor
          add_offset=MODELS(Iatmos)%ImportField(id)%add_offset
!
          Fval=ptr2d(IminP,JminP)
          MyFmin(1)=Fval
          MyFmax(1)=Fval
          MyFmin(2)=Fval
          MyFmax(2)=Fval
!
!  Load import data into RegCM component variable.
!
          SELECT CASE (TRIM(ADJUSTL(ImportNameList(ifld))))
!
!  Sea surface temperature from OCN component (C).
!
            CASE ('sst', 'SST')
              FieldName=TRIM(ImportNameList(ifld))
              DO i=IminP,ImaxP
                DO j=JminP,JmaxP
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  Fval=scale*ptr2d(i,j)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  importFields%sst(j,i)=Fval
                END DO
              END DO

# if defined SEA_ICE || defined OCNICE
!
!  Sea ice thickness (m), from either OCN or ICE.
!
            CASE ('sit')
              DO i=IminP,ImaxP
                DO j=JminP,JmaxP
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  Fval=scale*ptr2d(i,j)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  importFields%sit(j,i)=Fval
                END DO
              END DO
# endif
!
!  Land-sea or wet-dry mask (nondimensional), from OCN.
!
            CASE ('msk')
              DO i=IminP,ImaxP
                DO j=JminP,JmaxP
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  Fval=scale*ptr2d(i,j)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  importFields%msk(j,i)=Fval
                END DO
              END DO
!
!  Surface roughness length scale (m), from WAV.
!
            CASE ('zo', 'Zo')
              DO i=IminP,ImaxP
                DO j=JminP,JmaxP
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  Fval=scale*ptr2d(i,j)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  importFields%zo(j,i)=Fval
                END DO
              END DO
!
!  Friction velocity (m s-1), from WAV.
!
            CASE ('ustar')
              DO i=IminP,ImaxP
                DO j=JminP,JmaxP
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  Fval=scale*ptr2d(i,j)+add_offset
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  importFields%ustar(j,i)=Fval
                END DO
              END DO
!
!  Import field not found.
!
            CASE DEFAULT
              IF (localPET.eq.0) THEN
                WRITE (cplout,10) TRIM(ImportNameList(ifld)),           &
     &                            TRIM(CinpName)
              END IF
              rc=ESMF_RC_NOT_FOUND
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
          IF (associated(ptr2d)) nullify (ptr2d)
        END DO DE_LOOP
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
!  Report import field information.
!
        IF ((DebugLevel.ge.0).and.(localPET.eq.0)) THEN
          WRITE (cplout,20) TRIM(ImportNameList(ifld)),                 &
     &                      TRIM(Time_CurrentString), ng,               &
     &                      Fmin(1), Fmax(1)
          IF (scale.ne.1.0_dp) THEN
            WRITE (cplout,30) Fmin(2), Fmax(2),                         &
     &                        ' regcmScale = ', scale
          ELSE IF (add_offset.ne.0.0_dp) THEN
            WRITE (cplout,30) Fmin(2), Fmax(2),                         &
     &                        ' AddOffset   = ', add_offset
          END IF
        END IF
!
!  Debugging: write out import field into a NetCDF file.
!
        IF ((DebugLevel.ge.3).and.                                      &
     &      MODELS(Iatmos)%ImportField(ifld)%debug_write) THEN
          WRITE (ofile,40) ng, TRIM(ImportNameList(ifld)),              &
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
!  Update RegCM import calls counter.
!
      IF (ImportCount.gt.0) THEN
        MODELS(Iatmos)%ImportCalls=MODELS(Iatmos)%ImportCalls+1
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_Import',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (/,2x,'RegCM_Import - unable to find option to import: ',  &
     &        a,/,18x,'check ''Import(atmos)'' in input script: ', a)
  20  FORMAT (2x,'RegCM_Import - ESMF: importing field ''',a,'''',      &
     &        t72,a,2x,'Grid ',i2.2,/,                                  &
     &        19x,'(Dmin = ', 1p,e15.8,0p,' Dmax = ',1p,e15.8,0p,')')
  30  FORMAT (19x,'(Cmin = ', 1p,e15.8,0p,' Cmax = ',1p,e15.8,0p,       &
     &        a,1p,e15.8,0p,')')
  40  FORMAT ('regcm_',i2.2,'_import_',a,'_',i4.4,2('-',i2.2),'_',      &
     &        i2.2,2('.',i2.2),'.nc')

      RETURN
      END SUBROUTINE RegCM_Import
!
      SUBROUTINE RegCM_Export (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Exports RegCM fields to other coupled gridded components.           !
!                                                                      !
!=======================================================================
!
      USE mod_update,   ONLY : exportFields
      USE mod_dynparam, ONLY : ici1, ici2, jci1, jci2
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng
      integer, intent(out) :: rc

      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: ifld, i, is, j
      integer :: Istr, Iend, Jstr, Jend
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: ExportCount
      integer :: localDE, localDEcount, localPET, PETcount
!
      real (dp), parameter :: pi = 3.14159265358979323846_dp
      real (dp) :: Fseconds, TimeInDays, Time_Current

      real (dp) :: MyFmax(1), MyFmin(1), Fmin(1), Fmax(1), Fval
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      character (len=22)      :: Time_CurrentString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_Export"
!
      character (ESMF_MAXSTR) :: cname, ofile
      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
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
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_Export',            &
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
     &                       petCount=PETcount,                         &
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
!  Get number of local decomposition elements (DEs). Usually, a single
!  DE is associated with each Persistent Execution Thread (PETs). Thus,
!  localDEcount=1.
!
      CALL ESMF_GridGet (MODELS(Iatmos)%grid(ng),                       &
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
      Fseconds=REAL(seconds,dp)+REAL(sN,dp)/REAL(sD,dp)
      TimeInDays=Time_Current/86400.0_dp
      is=INDEX(Time_CurrentString, 'T')              ! remove 'T' in
      IF (is.gt.0) Time_CurrentString(is:is)=' '     ! ISO 8601 format
!
!-----------------------------------------------------------------------
!  Get list of export fields.
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet (MODELS(Iatmos)%ExportState(ng),               &
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
      CALL ESMF_StateGet (MODELS(Iatmos)%ExportState(ng),               &
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
!  Rotate vector components (wind, wind stress) to a rectangular
!  north-south and east-west oriented grid.
!-----------------------------------------------------------------------
!
      CALL RegCM_uvrot (exportFields%wndu, exportFields%wndv)
      CALL RegCM_uvrot (exportFields%taux, exportFields%tauy)
!
!-----------------------------------------------------------------------
!  Load export fields.
!-----------------------------------------------------------------------
!
      FLD_LOOP : DO ifld=1,ExportCount
!
!   Get field from export state.
!
        CALL ESMF_StateGet (MODELS(Iatmos)%ExportState(ng),             &
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
     &                        farrayPtr=ptr2d,                          &
     &                        rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          Istr=ici1
          Iend=ici2
          Jstr=jci1
          Jend=jci2
!
!  Initialize pointer to missing value.
!
          ptr2d=MISSING_dp
!
!  Load field data into export state.  Notice that all export fields
!  are kept as computed by RegCM. The imported component does the
!  proper scaling, physical units conversion, and other manipulations.
!  It is done to avoid applying such transformations twice.
!
          SELECT CASE (TRIM(ADJUSTL(ExportNameList(ifld))))
!
!  Surface air pressure (hPa or mb).
!
            CASE ('psfc', 'Pair')
              Fval=exportFields%psfc(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%psfc(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  2 meter surface air temperature (K).
!
            CASE ('tsfc', 'Tair')
              Fval=exportFields%tsfc(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%tsfc(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  2 meter surface specific humidity (kg/kg)
!
            CASE ('qsfc', 'Hair')
              Fval=exportFields%qsfc(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%qsfc(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface net longwave radiation (W m-2)
!
            CASE ('lwrd', 'LWrad')
              Fval=exportFields%lwrd(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%lwrd(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface downward longwave radiation (W m-2).
!
            CASE ('dlwrd', 'dLWrad', 'lwrad_down')
              Fval=exportFields%dlwr(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%dlwr(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface latent heat flux (W m-2)
!
            CASE ('lhfx', 'LHfx')
              Fval=exportFields%lhfx(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%lhfx(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface sensible heat flux (W m-2)
!
            CASE ('shfx')
              Fval=exportFields%shfx(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%shfx(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Total precipitation rate (m s-1).
!
            CASE ('prec')
              Fval=exportFields%prec(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%prec(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface eastward wind component (m s-1).
!
            CASE ('Uwind', 'u10', 'wndu')
              Fval=exportFields%wndu(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%wndu(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface northward wind component (m s-1).
!
            CASE ('wndv')
              Fval=exportFields%wndv(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%wndv(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface net solar shortwave radiation (W m-2).
!
            CASE ('swrd')
              Fval=exportFields%swrd(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%swrd(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Downward shortwave radiation (W m-2).
!
            CASE ('dswr')
              Fval=exportFields%dswr(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%dswr(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface runoff over land (m s-1).
!
            CASE ('rnof')
              Fval=exportFields%rnof(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%rnof(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Sub-surface runoff over land (m s-1).
!
            CASE ('snof')
              Fval=exportFields%snof(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%snof(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface eastward wind stress component (N m-2 or Pa).
!
            CASE ('taux')
              Fval=exportFields%taux(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%taux(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface eastward wind stress component (N m-2 or Pa).
!
            CASE ('tauy')
              Fval=exportFields%tauy(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%tauy(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface wind speed ( m s-1).
!
            CASE ('wspd')
              Fval=exportFields%wspd(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%wspd(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface wind direction (radians).
!
            CASE ('wdir')
              Fval=ATAN2(exportFields%wndu(Jstr,Istr),                  &
     &                   exportFields%wndv(Jstr,Istr))
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=ATAN2(exportFields%wndu(j,i),                    &
     &                       exportFields%wndv(j,i))
                  IF (dd.lt.0.0_r8) Fval=Fval+2.0_r8*pi
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface net heat flux (W m-2).
!
            CASE ('nflx')
              Fval=exportFields%nflx(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%nflx(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface net freshwater flux: E-P (m s-1).
!
            CASE ('sflx')
              Fval=exportFields%sflx(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%sflx(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface snow concentration.
!
            CASE ('snow')
              Fval=exportFields%snow(Jstr,Istr)
              MyFmin(1)=Fval
              MyFmax(1)=Fval
              DO i=Istr,Iend
                DO j=Jstr,Jend
                  Fval=exportFields%snow(j,i)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Export field not found.
!
            CASE DEFAULT
              IF (localPET.eq.0) THEN
                WRITE (cplout,10) TRIM(ADJUSTL(ExportNameList(ifld))),  &
     &                            TRIM(CinpName)
              END IF
              rc=ESMF_RC_NOT_FOUND
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
          IF (associated(ptr2d)) nullify (ptr2d)
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
!  Report export field information.
!
        IF ((DebugLevel.ge.0).and.(localPET.eq.0)) THEN
          WRITE (cplout,20) TRIM(ExportNameList(ifld)),                 &
     &                      TRIM(Time_CurrentString), ng,               &
     &                      Fmin(1), Fmax(1)
        END IF
!
!  Debugging: write out export field into a NetCDF file.
!
        IF ((DebugLevel.ge.3).and.                                      &
     &      MODELS(Iatmos)%ExportField(ifld)%debug_write) THEN
          WRITE (ofile,30) ng, TRIM(ExportNameList(ifld)),              &
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
      IF (allocated(ExportNameList)) deallocate(ExportNameList)
!
!  Update RegCM export calls counter.
!
      IF (ExportCount.gt.0) THEN
        MODELS(Iatmos)%ExportCalls=MODELS(Iatmos)%ExportCalls+1
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_Export',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (/,2x,'RegCM_Export - unable to find option to export: ',  &
     &        a,/,18x,'check ''Export(atmos)'' in input script: ',a)
  20  FORMAT (2x,'RegCM_Export - ESMF: exporting field ''',a,'''',      &
     &        t72,a,2x,'Grid ',i2.2,/,                                  &
     &        19x,'(Cmin = ', 1p,e15.8,0p,' Cmax = ',1p,e15.8,0p,')')
  30  FORMAT ('regcm',i2.2,'_export_',a,'_',i4.4,2('-',i2.2),'_',       &
     &        i2.2,2('.',i2.2),'.nc')

      RETURN
      END SUBROUTINE RegCM_Export
!
      SUBROUTINE RegCM_uvrot (u, v)
!
!=======================================================================
!                                                                      !
!  Rotates RegCM vector components to a rectangular north-south and    !
!  east-west oriented grid.                                            !
!                                                                      !
!=======================================================================
!
      USE mod_constants,     ONLY : degrad
      USE mod_atm_interface, ONLY : mddom
      USE mod_dynparam,      ONLY : iproj
      USE mod_dynparam,      ONLY : clon, clat, plon, plat, xcone
      USE mod_dynparam,      ONLY : ici1, ici2, jci1, jci2
!
!  Imported variable declarations.
!
      real (r8), intent(inout) :: u(jci1:jci2,ici1:ici2)
      real (r8), intent(inout) :: v(jci1:jci2,ici1:ici2)
!
!  Local variable declarations.
!
      integer :: i, j
!
      real (r8) :: x, xs, xc, d, us, vs, sindel, cosdel
      real (r8) :: pollam, polphi, polcphi, polsphi
      real (r8) :: zarg1, zarg2, znorm, zphi, zrla, zrlap
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", RegCM_uvrot"
!
!-----------------------------------------------------------------------
!  Rotate vector components to true EAST and true NORTH.
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering RegCM_uvrot',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
!  Rotated Mercator (ROTMER) or Normal Mercator (NORMER) projections.
!
      IF ((iproj.eq.'ROTMER').or.(iproj.eq.'NORMER')) THEN
        IF (plat.gt.0.0_r8) THEN
          pollam=plon+180.0_r8
          polphi=90.0_r8-plat
        ELSE
          polphi=90.0_r8+plat
          pollam=plon
        END if
        IF (pollam.gt.180.0_r8) pollam=pollam-360.0_r8
!
        polcphi=DCOS(degrad*polphi)
        polsphi=DSIN(degrad*polphi)
!
        DO j=jci1,jci2
          DO i=ici1,ici2
            zphi=mddom%dlat(j,i)*degrad
            zrla=mddom%dlon(j,i)*degrad
            IF (mddom%dlat(j,i).gt.89.999999_r8) zrla=0.0_r8
            zrlap=pollam*degrad-zrla
            zarg1=polcphi*dsin(zrlap)
            zarg2=polsphi*COS(zphi)-polcphi*SIN(zphi)*dcos(zrlap)
            znorm=1.0_r8/SQRT(zarg1**2+zarg2**2)
            sindel=zarg1*znorm
            cosdel=zarg2*znorm
!
            us= u(j,i)*cosdel+v(j,i)*sindel
            vs=-u(j,i)*sindel+v(j,i)*cosdel
            u(j,i)=us
            v(j,i)=vs
          END DO
        END DO
!
!  Otherwise, rotate from Lambert conformal (LAMCON) projection.
!
      ELSE
        DO i=ici1,ici2
          DO j=jci1,jci2
            IF (((clon.ge.0.0_r8).and.(mddom%xlon(j,i).ge.0.0_r8)).or.  &
     &          ((clon.lt.0.0_r8).and.(mddom%xlon(j,i).lt.0.0_r8))) THEN
              x=(clon-mddom%xlon(j,i))*degrad*xcone
            ELSE
              IF (clon.ge.0.0_r8) THEN
                IF (ABS(clon-(mddom%xlon(j,i)+360.0_r8)).lt.            &
     &              ABS(clon-mddom%xlon(j,i))) THEN
                  x=(clon-(mddom%xlon(j,i)+360.0d0))*degrad*xcone
                ELSE
                  x=(clon-mddom%xlon(j,i))*degrad*xcone
                END IF
              ELSE
                IF (ABS(clon-(mddom%xlon(j,i)-360.0_r8)).lt.            &
     &              ABS(clon-mddom%xlon(j,i))) THEN
                  x=(clon-(mddom%xlon(j,i)-360.0_r8))*degrad*xcone
                ELSE
                  x=(clon-mddom%xlon(j,i))*degrad*xcone
                END IF
              END IF
            END IF
!
            xs=SIN(x)
            xc=COS(x)
!
            IF (clat.ge.0.0_r8) THEN
              d=u(j,i)*xc-v(j,i)*xs
              v(j,i)=u(j,i)*xs+v(j,i)*xc
              u(j,i)=d
            ELSE
              d=u(j,i)*xc+v(j,i)*xs
              v(j,i)=v(j,i)*xc-u(j,i)*xs
              u(j,i)=d
            END IF
          END DO
        END DO
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  RegCM_uvrot',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE RegCM_uvrot
!
#endif
      END MODULE esmf_regcm_mod
