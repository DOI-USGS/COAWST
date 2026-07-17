      MODULE esmf_wrf_mod

#if defined WRF_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_atm_wrf.h 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module sets WRF as the atmospheric model gridded component     !
!  generic ESMF/NUOPC layer:                                           !
!                                                                      !
!    ATM_SetServices         Sets ATM component shared-object entry    !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                                                                      !
!    WRF_SetInitializeP1     WRF component phase 1 initialization:     !
!                            sets import and export fields long and    !
!                            short names into its respective state.    !
!                                                                      !
!    WRF_SetInitializeP2     WRF component phase 2 initialization:     !
!                            Initializes component (WRF_Initialize),   !
!                            sets component grid (WRF_SetGridArrays),  !
!                            and adds fields into import and export    !
!                            into respective states (WRF_SetStates).   !
!                                                                      !
!    WRF_DataInit            Exports WRF component fields during       !
!                            initialization or restart.                !
!                                                                      !
!    WRF_SetClock            Sets WRF component date calendar,         !
!                            start and stop times, and coupling        !
!                            interval.                                 !
# ifdef ESM_SETRUNCLOCK
!                                                                      !
!    WRF_SetRunClock         Sets WRF run clock manually.              !
# endif
!                                                                      !
!    WRF_CheckImport         Checks if ROMS component import field is  !
!                            at the correct time.                      !
!                                                                      !
!    WRF_SetGridArrays       Sets WRF component horizontal grid        !
!                            arrays, grid area, and land/sea mask.     !
!                                                                      !
!    WRF_SetStates           Adds WRF component export and import      !
!                            fields into its respective state.         !
!                                                                      !
!    WRF_ModelAdvance        Advances WRF  component for a coupling    !
!                            interval. It calls import and export      !
!                            routines.                                 !
!                                                                      !
!    WRF_SetFinalize         Finalizes WRF component execution.        !
!                                                                      !
!    WRF_Import              Imports fields into WRF from other        !
!                            gridded components.                       !
!                                                                      !
!    WRF_ProcessImport       Loads or merges import WRF fields.        !
!                                                                      !
!    WRF_Export              Exports WRF fields to other gridded       !
!                            components.                               !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!  WRF:    Weather Research and Forecasting model:                     !
!            http://www.wrf-model.org                                  !
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
      USE module_wrf_top, ONLY : WRF_Initialize => wrf_init,            &
     &                           WRF_Run,                               &
     &                           WRF_Finalize
!
      implicit none
!
      PUBLIC  :: ATM_SetServices

      PRIVATE :: WRF_SetInitializeP1
      PRIVATE :: WRF_SetInitializeP2
      PRIVATE :: WRF_DataInit
      PRIVATE :: WRF_SetClock
# ifdef ESM_SETRUNCLOCK
      PRIVATE :: WRF_SetRunClock
# endif
      PRIVATE :: WRF_CheckImport
      PRIVATE :: WRF_SetGridArrays
      PRIVATE :: WRF_SetStates
      PRIVATE :: WRF_ModelAdvance
      PRIVATE :: WRF_SetFinalize
      PRIVATE :: WRF_Import
      PRIVATE :: WRF_ProcessImport
      PRIVATE :: WRF_Export
!
      CONTAINS
!
      SUBROUTINE ATM_SetServices (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets WRF component shared-object entry points for "initialize",     !
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
     &                              userRoutine=WRF_SetInitializeP1,    &
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
     &                              userRoutine=WRF_SetInitializeP2,    &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Attach WRF component phase independent specializing methods.
!-----------------------------------------------------------------------
!
!  Set routine for export initial/restart fields.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_DataInitialize,  &
     &                           specRoutine=WRF_DataInit,              &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for setting WRF clock.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_SetClock,        &
     &                           specRoutine=WRF_SetClock,              &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF

# ifdef ESM_SETRUNCLOCK
!
!  Set routine for setting WRF run clock manually. First, remove the
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
     &                           specRoutine=WRF_SetRunClock,           &
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
     &                           specRoutine=WRF_CheckImport,           &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for time-stepping WRF component.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_Advance,         &
     &                           specRoutine=WRF_ModelAdvance,          &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Register WRF finalize routine.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompSetEntryPoint (model,                           &
     &                                 methodflag=ESMF_METHOD_FINALIZE, &
     &                                 userRoutine=WRF_SetFinalize,     &
     &                                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ATM_SetServices',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE ATM_SetServices
!
      SUBROUTINE WRF_SetInitializeP1 (model,                            &
     &                                ImportState, ExportState,         &
     &                                clock, rc)
!
!=======================================================================
!                                                                      !
!  WRF component Phase 1 initialization: sets import and export fields !
!  long and short names into its respective state.                     !
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
     &  __FILE__//", WRF_SetInitializeP1"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_SetInitializeP1',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set WRF import state and fields.
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
!  Set WRF export state and fields.
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_SetInitializeP1',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WRF_SetInitializeP1
!
      SUBROUTINE WRF_SetInitializeP2 (model,                            &
     &                                ImportState, ExportState,         &
     &                                clock, rc)
!
!=======================================================================
!                                                                      !
!  WRF component Phase 2 initialization: Initializes WRF, sets         !
!  component grid, and adds import and export fields to respective     !
!  states.                                                             !
!                                                                      !
!=======================================================================
!
      USE module_domain
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
      integer :: MyComm, ng
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_SetInitializeP2"
!
      TYPE (domain), pointer :: nest               ! pointer to WRF nest
!
      TYPE (ESMF_VM)   :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_SetInitializeP2',     &
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
!  Initialize WRF component. In nested applications, WRF kernel
!  will allocate and initialize all grids with a single call to
!  "WRF_initialize".
!-----------------------------------------------------------------------
!
      CALL wrf_set_dm_communicator (MyComm)
      CALL WRF_Initialize ()
!
!-----------------------------------------------------------------------
!  Set-up grid and load coordinate data.
!-----------------------------------------------------------------------
!
      DO ng=1,MODELS(Iatmos)%Ngrids
        IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
          IF (ng.eq.1) THEN
            CALL WRF_SetGridArrays (head_grid, model, rc)
          ELSE
            CALL find_grid_by_id (ng, head_grid, nest)
            CALL WRF_SetGridArrays (nest, model, rc)
          END IF
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
          IF (ng.eq.1) THEN
            CALL WRF_SetStates (head_grid, model, rc)
          ELSE
            CALL find_grid_by_id (ng, head_grid, nest)
            CALL WRF_SetStates (nest, model, rc)
          END IF
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_SetInitializeP2',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WRF_SetInitializeP2
!
      SUBROUTINE WRF_DataInit (model, rc)
!
!=======================================================================
!                                                                      !
!  Exports WRF component fields during initialization or restart.      !
!                                                                      !
!=======================================================================
!
      USE module_domain
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
      integer :: localPET, phase
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_DataInit"
!
      TYPE (domain), pointer :: nest               ! pointer to WRF nest
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_DataInit',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Inquire gridded component.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       localPet=localPET,                         &
     &                       currentPhase=phase,                        &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  If explicit coupling from atmosphere to ocean, export initialization
!  or restart fields.
!-----------------------------------------------------------------------
!
      IF ((CouplingType.eq.0).and.(Nexport(Iatmos).gt.0)) THEN
        DO ng=1,MODELS(Iatmos)%Ngrids
          IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
            IF (ng.eq.1) THEN
              CALL WRF_Export (head_grid, model, rc)
            ELSE
              CALL find_grid_by_id (ng, head_grid, nest)
              CALL WRF_Export (nest, model, rc)
            END IF
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_DataInit',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WRF_DataInit
!
      SUBROUTINE WRF_SetClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets WRF component date calendar, start and stop time, and coupling !
!  interval.                                                           !
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
      integer :: ig, is, ng
      integer :: ref_year,   start_year,   stop_year
      integer :: ref_month,  start_month,  stop_month
      integer :: ref_day,    start_day,    stop_day
      integer :: ref_hour,   start_hour,   stop_hour
      integer :: ref_minute, start_minute, stop_minute
      integer :: ref_second, start_second, stop_second
      integer :: localPET
      integer :: TimeFrac
!
      character (len= 22) :: Calendar
      character (len= 22) :: StartTimeString, StopTimeString
      character (len=160) :: message

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_SetClock"
!
      TYPE (ESMF_CalKind_Flag) :: CalType
      TYPE (ESMF_Clock)        :: clock
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_SetClock',            &
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
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Create WRF component clock.
!-----------------------------------------------------------------------
!
      Calendar=TRIM(ClockInfo(Iatmos)%CalendarString)
      IF (TRIM(Calendar).eq.'gregorian') THEN
        CalType=ESMF_CALKIND_GREGORIAN
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
      ng=1
      CALL nl_get_simulation_start_year   (ng, ref_year)
      CALL nl_get_simulation_start_month  (ng, ref_month)
      CALL nl_get_simulation_start_day    (ng, ref_day)
      CALL nl_get_simulation_start_hour   (ng, ref_hour)
      CALL nl_get_simulation_start_minute (ng, ref_minute)
      CALL nl_get_simulation_start_second (ng, ref_second)
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
      ng=1
      CALL nl_get_start_year   (ng, start_year)
      CALL nl_get_start_month  (ng, start_month)
      CALL nl_get_start_day    (ng, start_day)
      CALL nl_get_start_hour   (ng, start_hour)
      CALL nl_get_start_minute (ng, start_minute)
      CALL nl_get_start_second (ng, start_second)
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

# ifdef REGRESS_STARTCLOCK
!
!  Substract a coupling interval because the driver clock was regressed
!  by that amount to properly initialize all ESM components. Notice that
!  above routines returns the values specified in its "namelist" input
!  file.
!
      ClockInfo(Iatmos)%StartTime=ClockInfo(Iatmos)%StartTime-          &
     &                            ClockInfo(Iatmos)%TimeStep
      CALL ESMF_TimeGet (ClockInfo(Iatmos)%StartTime,                   &
     &                   timeStringISOFrac=StartTimeString,             &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# else
!
      CALL ESMF_TimeGet (ClockInfo(Iatmos)%StartTime,                   &
     &                   timeStringISOFrac=StartTimeString,             &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# endif
      is=INDEX(StartTimeString, 'T')                ! remove 'T' in
      IF (is.gt.0) StartTimeString(is:is)=' '       ! ISO 8601 format
      ClockInfo(Iatmos)%Time_StartString=StartTimeString
!
!  Set stop time.
!
      ng=1
      CALL nl_get_end_year   (ng, stop_year)
      CALL nl_get_end_month  (ng, stop_month)
      CALL nl_get_end_day    (ng, stop_day)
      CALL nl_get_end_hour   (ng, stop_hour)
      CALL nl_get_end_minute (ng, stop_minute)
      CALL nl_get_end_second (ng, stop_second)
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
      CALL ESMF_TimeGet (ClockInfo(Iatmos)%StopTime,                    &
     &                   timeStringISOFrac=StopTimeString,              &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(StopTimeString, 'T')                 ! remove 'T' in
      IF (is.gt.0) StopTimeString(is:is)=' '        ! ISO 8601 format
      ClockInfo(Iatmos)%Time_StopString=StopTimeString
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
!-----------------------------------------------------------------------
!  Create WRF component clock.
!-----------------------------------------------------------------------
!
      ClockInfo(Iatmos)%Name='WRF_clock'
      clock=ESMF_ClockCreate(ClockInfo(Iatmos)%TimeStep,                &
     &                       ClockInfo(Iatmos)%StartTime,               &
     &                       stopTime =ClockInfo(Iatmos)%StopTime,      &
     &                       refTime  =ClockInfo(Iatmos)%ReferenceTime, &
     &                       name     =TRIM(ClockInfo(Iatmos)%Name),    &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      ClockInfo(Iatmos)%Clock=clock
!
!  Set ROMS component clock.
!
      CALL ESMF_GridCompSet (model,                                     &
     &                       clock=ClockInfo(Iatmos)%Clock,             &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get current time.
!
      CALL ESMF_ClockGet (ClockInfo(Iatmos)%Clock,                      &
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
!  Compare driver time against WRF component time.
!-----------------------------------------------------------------------
!
      IF (ClockInfo(Idriver)%Restarted) THEN
        StartTimeString=ClockInfo(Idriver)%Time_RestartString
      ELSE
        StartTimeString=ClockInfo(Idriver)%Time_StartString
      END IF
!
      IF (ClockInfo(Iatmos)%Time_StartString.ne.                        &
     &    StartTimeString) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,10) 'WRF    Start Time: ',                      &
     &                      TRIM(ClockInfo(Iatmos)%Time_StartString),   &
     &                      'Driver Start Time: ',                      &
     &                      TRIM(StartTimeString),                      &
     &                      '                   are not equal!'
        END IF
        message='Driver and WRF start times do not match: '//           &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ClockInfo(Iatmos )%Time_StopString(1:19).ne.                  &
     &    ClockInfo(Idriver)%Time_StopString(1:19)) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,10) 'WRF    Stop Time: ',                       &
     &                      TRIM(ClockInfo(Iatmos )%Time_StopString),   &
     &                      'Driver Stop Time: ',                       &
     &                      TRIM(ClockInfo(Idriver)%Time_StopString),   &
     &                      '                   are not equal!'
        END IF
        message='Driver and WRF stop times do not match: '//            &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (TRIM(ClockInfo(Iatmos )%CalendarString).ne.                   &
     &    TRIM(ClockInfo(Idriver)%CalendarString)) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,10) 'WRF    Calendar: ',                        &
     &                      TRIM(ClockInfo(Iatmos )%CalendarString),    &
     &                      'Driver Calendar: ',                        &
     &                      TRIM(ClockInfo(Idriver)%CalendarString),    &
     &                      '                  are not equal!'
        END IF
        message='Driver and WRF calendars do not match: '//             &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_SetClock',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
 10   FORMAT (/,1x,a,a,/,1x,a,a,/,1x,a)
!
      RETURN
      END SUBROUTINE WRF_SetClock

# ifdef ESM_SETRUNCLOCK
!
      SUBROUTINE WRF_SetRunClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets WRF run clock manually to avoid getting zero time stamps at    !
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
     &  __FILE__//", WRF_SetRunClock"
!
      TYPE (ESMF_Clock) :: driverClock, modelClock
      TYPE (ESMF_Time)  :: currTime
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_SetRunClock',         &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_SetRunClock',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WRF_SetRunClock
# endif
!
      SUBROUTINE WRF_CheckImport (model, rc)
!
!=======================================================================
!                                                                      !
!  Checks if WRF component import field is at the correct time.        !
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
     &  __FILE__//", WRF_CheckImport"
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
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_CheckImport',         &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_CheckImport',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (1x,'WRF_CheckImport - ',a,':',t32,'TimeStamp = ',a,       &
     &        ',  DriverTime = ',a)
!
      RETURN
      END SUBROUTINE WRF_CheckImport
!
      SUBROUTINE WRF_SetGridArrays (grid, model, rc)
!
!=======================================================================
!                                                                      !
!  Sets WRF component staggered, horizontal grids arrays, grid area,   !
!  and land/sea mask.                                                  !
!                                                                      !
!=======================================================================
!
      USE module_domain, ONLY : domain,                                 &
     &                          get_ijk_from_grid
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (domain), intent(in)           :: grid
      TYPE (ESMF_GridComp), intent(inout) :: model
!
!  Local variable declarations.
!
      integer :: i, j, k, ng
      integer :: gtype, ivar
      integer :: ids, ide, jds, jde, kds, kde
      integer :: ims, ime, jms, jme, kms, kme
      integer :: ips, ipe, jps, jpe, kps, kpe
      integer :: Im, Jm
      integer :: IstrD, IendD, JstrD, JendD
      integer :: IstrM, IendM, JstrM, JendM
      integer :: IstrP, IendP, JstrP, JendP
!
      integer :: localPET, PETcount, node
      integer :: localDE, localDEcount
      integer :: NumProcsX, NumProcsY
      integer :: LBi, UBi, LBj, UBj
      integer :: cLB(2), cUB(2), eLB(2), eUB(2), tLB(2), tUB(2)
      integer :: minIndex(2), maxIndex(2)
!
      integer, allocatable :: IpatchStarts(:), JpatchStarts(:)
      integer, allocatable :: IpatchEnds(:),   JpatchEnds(:)
      integer, allocatable :: deBlockList(:,:,:)
!
      integer (i4b), pointer :: ptrM(:,:) => NULL()
!
      real (dp), pointer :: ptrA(:,:) => NULL()
      real (dp), pointer :: ptrX(:,:) => NULL()
      real (dp), pointer :: ptrY(:,:) => NULL()
!
      character (len=40) :: name

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_SetGridArrays"
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
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_SetGridArrays',       &
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
      ng=grid%grid_id
!
!-----------------------------------------------------------------------
!  Get variables related with grid partitioning: global indecing
!
!    ids, ide, jds, jde, kds, kde  =>  domain extent
!    ims, ime, jms, jme, kms, kme  =>  memory extent
!    ips, ipe, jps, jpe, kps, kpe  =>  patch extent
!-----------------------------------------------------------------------
!
      CALL get_ijk_from_grid (grid, ids, ide, jds, jde, kds, kde,       &
     &                              ims, ime, jms, jme, kms, kme,       &
     &                              ips, ipe, jps, jpe, kps, kpe)
!
!  For reference, get the indices from WRF state structure.
!
      IstrD=grid%sd31        ! Domain extent
      IendD=grid%ed31
      JstrD=grid%sd33
      JendD=grid%ed33
!
      IstrM=grid%sm31        ! Memory extent
      IendM=grid%em31
      JstrM=grid%sm33
      JendM=grid%em33
!
      IstrP=grid%sp31        ! Patch extent
      IendP=grid%ep31
      JstrP=grid%sp33
      JendP=grid%ep33
!
!-----------------------------------------------------------------------
!  Calculate patch and number of CPUs in each direction.
!-----------------------------------------------------------------------
!
!  Starting patch.
!
      IF (.not.allocated(IpatchStarts)) THEN
        allocate (IpatchStarts(0:PETcount-1))
      END IF
      IF (.not.allocated(JpatchStarts)) THEN
        allocate (JpatchStarts(0:PETcount-1))
      END IF
!
      CALL ESMF_VMAllGatherV (vm,                                       &
     &                        sendData=(/ips/),                         &
     &                        sendCount=1,                              &
     &                        recvData=IpatchStarts,                    &
     &                        recvCounts =(/(1, k=0, PETcount-1)/),     &
     &                        recvOffsets=(/(k, k=0, PETcount-1)/),     &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_VMAllGatherV (vm,                                       &
     &                        sendData=(/jps/),                         &
     &                        sendCount=1,                              &
     &                        recvData=JpatchStarts,                    &
     &                        recvCounts =(/(1, k=0, PETcount-1)/),     &
     &                        recvOffsets=(/(k, k=0, PETcount-1)/),     &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      NumProcsX=0
      NumProcsY=0
      DO node=0,PETcount-1
        IF (ips.eq.IpatchStarts(node)) THEN
          NumProcsY=NumProcsY+1
        END IF
        IF (jps.eq.JpatchStarts(node)) THEN
          NumProcsX=NumProcsX+1
        END IF
      END DO
!
!  Ending patch.
!
      IF (.not.allocated(IpatchEnds)) THEN
        allocate (IpatchEnds(0:PETcount-1))
      END IF
      IF (.not.allocated(JpatchEnds)) THEN
        allocate (JpatchEnds(0:PETcount-1))
      END IF
!
      CALL ESMF_VMAllGatherV (vm,                                       &
     &                        sendData=(/MIN(ide-1,ipe)/),              &
     &                        sendCount=1,                              &
     &                        recvData=IpatchEnds,                      &
     &                        recvCounts =(/(1, k=0, PETcount-1)/),     &
     &                        recvOffsets=(/(k, k=0, PETcount-1)/),     &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_VMAllGatherV (vm,                                       &
     &                        sendData=(/MIN(jde-1,jpe)/),              &
     &                        sendCount=1,                              &
     &                        recvData=JpatchEnds,                      &
     &                        recvCounts =(/(1, k=0, PETcount-1)/),     &
     &                        recvOffsets=(/(k, k=0, PETcount-1)/),     &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Create ESMF DistGrid based on model domain decomposition.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(deBlockList)) THEN
        allocate (deBlockList(2,2,PETcount))
      END IF
!
      DO node=1,PETCount
        deBlockList(1,1,node)=IpatchStarts(node-1)
        deBlockList(2,1,node)=JpatchStarts(node-1)
        deBlockList(1,2,node)=IpatchEnds(node-1)
        deBlockList(2,2,node)=JpatchEnds(node-1)
      END DO
!
      Im=MAXVAL(IpatchEnds)
      Jm=MAXVAL(JpatchEnds)
      minIndex=(/1, 1/)
      maxIndex=(/Im, Jm/)
!
      distGrid=ESMF_DistGridCreate(minIndex=minIndex,                   &
     &                             maxIndex=maxIndex,                   &
     &                             deBlockList=deBlockList,             &
     &                             rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Report WRF DistGrid based on model domain decomposition.
!
      IF ((localPET.eq.0).and.(DebugLevel.gt.0)) THEN
        WRITE (cplout,10) ng, TRIM(GridType(Icenter))//" Point",        &
     &                    NumProcsX, NumProcsY
        DO node=1,PETCount
          WRITE (cplout,20) node-1, deBlockList(1,1,node),              &
     &                              deBlockList(1,2,node),              &
     &                              deBlockList(2,1,node),              &
     &                              deBlockList(2,2,node)
        END DO
      END IF
      IF (allocated(deBlockList)) deallocate (deBlockList)

# ifdef DATA_COUPLING
!
!  Read in melding weights coefficients needed by WRF to merge imported
!  fields from DATA and other ESM components at the specified nested
!  grid because of incongruent grids.
!
      IF ((MODELS(Idata)%IsActive).and.                                 &
     &    (ng.eq.WEIGHTS(Iatmos)%NestedGrid)) THEN
        CALL get_weights (Iatmos, Im, Jm, vm, rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Set component grid coordinates.
!-----------------------------------------------------------------------
!
!  Define component grid location type: center cell.
!
      IF (.not.allocated(MODELS(Iatmos)%mesh)) THEN
        allocate ( MODELS(Iatmos)%mesh(1) )
        MODELS(Iatmos)%mesh(1)%gtype=Icenter
      END IF
!
!  Create ESMF Grid.
!
      MODELS(Iatmos)%grid(ng)=ESMF_GridCreate(distgrid=distGrid,        &
     &                                  coordSys=ESMF_COORDSYS_SPH_DEG, &
     &                                  coordTypeKind=ESMF_TYPEKIND_R8, &
     &                                  indexflag=ESMF_INDEX_GLOBAL,    &
     &                                  name=TRIM(MODELS(Iatmos)%name), &
     &                                  rc=rc)
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
!  Set staggering type, Arakawa C-grid.
!
        SELECT CASE (MODELS(Iatmos)%mesh(ivar)%gtype)
          CASE (Icenter)
            staggerLoc=ESMF_STAGGERLOC_CENTER
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
!  The WRF masking is as follows,  1: land
!                                  0: ocean
!
        MODELS(Iatmos)%LandValue=1
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
     &                            exclusiveLBound=eLB,                  &
     &                            exclusiveUBound=eUB,                  &
     &                            computationalLBound=cLB,              &
     &                            computationalUBound=cUB,              &
     &                            totalLBound=tLB,                      &
     &                            totalUBound=tUB,                      &
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
     &                            exclusiveLBound=eLB,                  &
     &                            exclusiveUBound=eUB,                  &
     &                            computationalLBound=cLB,              &
     &                            computationalUBound=cUB,              &
     &                            totalLBound=tLB,                      &
     &                            totalUBound=tUB,                      &
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
            CASE (Icenter)
              LBi=LBOUND(ptrX,1)
              UBi=UBOUND(ptrX,1)
              LBj=LBOUND(ptrX,2)
              UBj=UBOUND(ptrX,2)
              DO j=LBj,UBj
                DO i=LBi,UBi
                  ptrX(i,j)=REAL(grid%xlong(i,j),dp)
                  ptrY(i,j)=REAL(grid%xlat(i,j),dp)
                  ptrM(i,j)=INT(grid%landmask(i,j))
                  ptrA(i,j)=REAL(grid%dx*grid%dy,dp)
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
        END DO DE_LOOP
!
!  Debugging: write out component grid in VTK format.
!
        IF (DebugLevel.ge.4) THEN
          gtype=MODELS(Iatmos)%mesh(ivar)%gtype
          CALL ESMF_GridWriteVTK (MODELS(Iatmos)%grid(ng),              &
     &                            filename="wrf_"//                     &
     &                                      TRIM(GridType(gtype))//     &
     &                                      "_point",                   &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_SetGridArrays',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (3x,'WRF_DistGrid - Grid = ',i2.2,',',3x,'Mesh = ',a,',',  &
     &        3x,'Partition = ',i0,' x ',i0)
  20  FORMAT (18x,'node = ',i0,t32,'Istr = ',i0,t45,'Iend = ',i0,       &
     &                         t58,'Jstr = ',i0,t71,'Jend = ',i0)
!
      RETURN
      END SUBROUTINE WRF_SetGridArrays
!
      SUBROUTINE WRF_SetStates (grid, model, rc)
!
!=======================================================================
!                                                                      !
!  Adds WRF component export and import fields into its respective     !
!  state.                                                              !
!                                                                      !
!=======================================================================
!
      USE module_domain, ONLY : domain,                                 &
     &                          get_ijk_from_grid
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (domain), intent(in) :: grid
      TYPE (ESMF_GridComp)      :: model
!
!  Local variable declarations.
!
      integer :: i, id, ng
      integer :: localDE, localDEcount
      integer :: localPET, PETcount
      integer :: ExportCount, ImportCount
!
      real (dp), dimension(:,:), pointer :: ptr2d => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_SetStates"
!
      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)
!
      TYPE (ESMF_ArraySpec)  :: arraySpec2d
      TYPE (ESMF_Field)      :: field
      TYPE (ESMF_StaggerLoc) :: staggerLoc
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_SetStates',           &
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
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      ng=grid%grid_id
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
              CASE (Iupoint)
                staggerLoc=ESMF_STAGGERLOC_EDGE1
              CASE (Ivpoint)
                staggerLoc=ESMF_STAGGERLOC_EDGE2
            END SELECT
!
!  Create 2D field from the Grid and arraySpec.
!
            field=ESMF_FieldCreate(MODELS(Iatmos)%grid(ng),             &
     &                             arraySpec2d,                         &
     &                             indexflag=ESMF_INDEX_GLOBAL,         &
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
              CASE (Iupoint)
                staggerLoc=ESMF_STAGGERLOC_EDGE1
              CASE (Ivpoint)
                staggerLoc=ESMF_STAGGERLOC_EDGE2
            END SELECT
!
!  Create 2D field from the Grid, arraySpec.
!
            field=ESMF_FieldCreate(MODELS(Iatmos)%grid(ng),             &
     &                             arraySpec2d,                         &
     &                             indexflag=ESMF_INDEX_GLOBAL,         &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_SetStates',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
!
 10   FORMAT ('WRF_SetStates - Removing field ''',a,''' from ',a,       &
     &        '''',a,'''',/,16x,'because it is not connected.')
!
      RETURN
      END SUBROUTINE WRF_SetStates
!
      SUBROUTINE WRF_ModelAdvance (model, rc)
!
!=======================================================================
!                                                                      !
!  Advance WRF component for a coupling interval (seconds) using       !
!  "WRF_run". It also calls "WRF_Import" and "WRF_Export" to import    !
!  and export coupling fields, respectively.                           !
!                                                                      !
!=======================================================================
!
      USE module_domain
      USE module_symbols_util, ONLY : WRFU_TimeSet
      USE WRF_ESMF_TimeMod,    ONLY : WRF_ESMF_Time => ESMF_Time
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
!
      integer :: ig, is, ng
      integer :: localPET, PETcount, phase
      integer :: TimeStart(6), TimeStop(6)
!
      real (dp) :: CouplingInterval, RunInterval
      real (dp) :: TcurrentInSeconds, TstopInSeconds
!
      character (len=22) :: Cinterval
      character (len=22) :: CurrTimeString, StopTimeString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_ModelAdvance"
!
      TYPE (domain), pointer :: nest               ! pointer to WRF nest
!
      TYPE (ESMF_Clock)        :: clock
      TYPE (ESMF_State)        :: ExportState, ImportState
      TYPE (ESMF_TimeInterval) :: TimeStep
      TYPE (ESMF_Time)         :: ReferenceTime
      TYPE (ESMF_Time)         :: CurrentTime, StopTime
      TYPE (ESMF_Time)         :: TimeFrom, TimeTo
      TYPE (ESMF_VM)           :: vm
!
      TYPE (WRF_ESMF_Time)     :: MyTimeFrom, MyTimeTo
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_ModelAdvance',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get information about the gridded component.
!-----------------------------------------------------------------------
!
!  Inquire about WRF component.
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
!  Get time step interval, stopping time, reference time, and current
!  time.
!
      CALL ESMF_ClockGet (clock,                                        &
     &                    timeStep=TimeStep,                            &
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
!  Current WRF time (seconds).
!
      CALL ESMF_TimeGet (ClockInfo(Iatmos)%CurrentTime,                 &
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
!  WRF stop time (seconds) for this coupling window.
!
      CALL ESMF_TimeGet (ClockInfo(Iatmos)%CurrentTime+TimeStep,        &
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

# ifdef REGRESS_STARTCLOCK
!
!  If regressed driver starting clock, avoid timestepping COAMPS during
!  the regressed coupling interval.
!
      IF (TcurrentInSeconds.gt.ClockInfo(Idriver)%Time_Start) THEN
        Ladvance=.TRUE.
      ELSE
        Ladvance=.FALSE.
      END IF
# else
!
!  Set model advance switch.
!
      Ladvance=.TRUE.
# endif
!
!-----------------------------------------------------------------------
!  Calculate run time for the current coupling window. Convert ESMF_Time
!  to WRF_ESMF_Time. Here, WRF_ESMF_Time is the old version of ESMF_Time
!  implemented in WRF.
!-----------------------------------------------------------------------
!
      IF (ClockInfo(Idriver)%Restarted) THEN
        TimeFrom=ClockInfo(Iatmos)%RestartTime
      ELSE
        TimeFrom=ClockInfo(Iatmos)%CurrentTime
      END IF
      TimeTo=TimeFrom+TimeStep
!
!  Coupling window starting time.
!
      CALL ESMF_TimeGet (TimeFrom,                                      &
     &                   yy=TimeStart(1),                               &
     &                   mm=TimeStart(2),                               &
     &                   dd=TimeStart(3),                               &
     &                   h =TimeStart(4),                               &
     &                   m =TimeStart(5),                               &
     &                   s =TimeStart(6),                               &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL WRFU_TimeSet (MyTimeFrom,                                    &
     &                   yy=TimeStart(1),                               &
     &                   mm=TimeStart(2),                               &
     &                   dd=TimeStart(3),                               &
     &                   h =TimeStart(4),                               &
     &                   m =TimeStart(5),                               &
     &                   s =TimeStart(6),                               &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Coupling window ending time.
!
      CALL ESMF_TimeGet (TimeTo,                                        &
     &                   yy=TimeStop(1),                                &
     &                   mm=TimeStop(2),                                &
     &                   dd=TimeStop(3),                                &
     &                   h =TimeStop(4),                                &
     &                   m =TimeStop(5),                                &
     &                   s =TimeStop(6),                                &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL WRFU_TimeSet (MyTimeTo,                                      &
     &                   yy=TimeStop(1),                                &
     &                   mm=TimeStop(2),                                &
     &                   dd=TimeStop(3),                                &
     &                   h =TimeStop(4),                                &
     &                   m =TimeStop(5),                                &
     &                   s =TimeStop(6),                                &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Report time information strings (YYYY-MM-DD hh:mm:ss).
!-----------------------------------------------------------------------
!
      IF (localPET.eq.0) THEN
        WRITE (Cinterval,'(f15.2)') CouplingInterval
        WRITE (cplout,10) TRIM(CurrTimeString), TRIM(StopTimeString),   &
     &                    TRIM(ADJUSTL(Cinterval)), Ladvance
      END IF
!
!-----------------------------------------------------------------------
!  Get import fields from other ESM components.
!-----------------------------------------------------------------------
!
      IF (Nimport(Iatmos).gt.0) THEN
        DO ng=1,MODELS(Iatmos)%Ngrids
          IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
            IF (ng.eq.1) THEN
              CALL WRF_Import (head_grid, model, rc)
            ELSE
              CALL find_grid_by_id (ng, head_grid, nest)
              CALL WRF_Import (nest, model, rc)
            END IF
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
!  Run WRF component. Notice that atmosphere component is advanced
!  when ng=1.  In nested application, its numerical kernel will advance
!  all the nested grids in their logical order.
!-----------------------------------------------------------------------
!
      IF (Ladvance) THEN
        head_grid%start_subtime=MyTimeFrom
        head_grid%stop_subtime =MyTimeTo
        IF (ESM_track) THEN
          WRITE (trac,'(a,a,i0)') '==> Entering WRF_Run',               &
     &                            ', PET', PETrank
          CALL my_flush (trac)
        END IF
        CALL WRF_Run ()
        IF (ESM_track) THEN
          WRITE (trac,'(a,a,i0)') '==> Exiting  WRF_Run',               &
     &                            ', PET', PETrank
          CALL my_flush (trac)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Put export fields.
!-----------------------------------------------------------------------
!
      IF (Nexport(Iatmos).gt.0) THEN
        DO ng=1,MODELS(Iatmos)%Ngrids
          IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
            IF (ng.eq.1) THEN
              CALL WRF_Export (head_grid, model, rc)
            ELSE
              CALL find_grid_by_id (ng, head_grid, nest)
              CALL WRF_Export (nest, model, rc)
            END IF
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_ModelAdvance',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (3x,'ModelAdvance - ESMF, Running WRF:',t42,a,             &
     &        ' => ',a,', [',a,' s], Advance: ',l1)
!
      RETURN
      END SUBROUTINE WRF_ModelAdvance
!
      SUBROUTINE WRF_SetFinalize (model,                                &
     &                            ImportState, ExportState,             &
     &                            clock, rc)
!
!=======================================================================
!                                                                      !
!  Finalize WRF component execution. It calls WRF_Finalize.            !
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
      logical :: no_shutdown
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_SetFinalize"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_SetFinalize',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Finalize WRF component.
!-----------------------------------------------------------------------
!
!  Set switch to avoid quiting the MPI communication servers.  It will
!  be done by the ESMF driver for all components.
!
      no_shutdown=.TRUE.
      CALL WRF_Finalize (no_shutdown)
      CALL my_flush (6)                   ! flush standard output buffer
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_SetFinalize',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE WRF_SetFinalize
!
      SUBROUTINE WRF_Import (grid, model, rc)
!
!=======================================================================
!                                                                      !
!  Imports fields into WRF array structure from other coupled gridded  !
!  components.                                                         !
!                                                                      !
!=======================================================================
!
      USE module_domain, ONLY : domain,                                 &
     &                          get_ijk_from_grid
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (domain), pointer :: grid
!
      TYPE (ESMF_GridComp)   :: model
!
!  Local variable declarations.
!
      logical :: got_sst(2)
!
      integer :: id, ifld, i, is, j, ng
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: ImportCount
      integer :: localDE, localDEcount, localPET, PETcount
      integer :: LBi, UBi, LBj, UBj
      integer :: IminP, ImaxP, JminP, JmaxP
      integer :: ids, ide, jds, jde, kds, kde
      integer :: ims, ime, jms, jme, kms, kme
      integer :: ips, ipe, jps, jpe, kps, kpe
      integer :: ifield(2)
!
      real (dp) :: Fseconds, TimeInDays, Time_Current

      real (dp) :: MyFmax(2), MyFmin(2), Fmin(2), Fmax(2), Fval
      real (dp) :: scale, add_offset
!
      real (dp), pointer     :: ptr2d(:,:)
!
      real (dp), allocatable :: dat_sst(:,:), ocn_sst(:,:)
!
      character (len=22 )     :: Time_CurrentString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_Import"
!
      character (ESMF_MAXSTR) :: FieldName, fld_name(2)
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
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_Import',              &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Compute WRF lower and upper bounds (non-overlapping) for physical
!  area per nested grid and patch.
!-----------------------------------------------------------------------
!
      IminP=grid%sp31
      ImaxP=grid%ep31
      JminP=grid%sp33
      JmaxP=grid%ep33
      IF (grid%ed31.eq.ImaxP) THEN
        ImaxP=ImaxP-1
      END IF
      IF (grid%ed33.eq.JmaxP) THEN
        JmaxP=JmaxP-1
      END IF
!
!  Get variables related with grid partitioning: head grid.
!
!    ids, ide, jds, jde, kds, kde  =>  domain extent
!    ims, ime, jms, jme, kms, kme  =>  memory extent
!    ips, ipe, jps, jpe, kps, kpe  =>  patch extent
!
      CALL get_ijk_from_grid (grid, ids, ide, jds, jde, kds, kde,       &
     &                              ims, ime, jms, jme, kms, kme,       &
     &                              ips, ipe, jps, jpe, kps, kpe)
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
      ng=grid%grid_id
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
      got_sst(1)=.FALSE.               ! SST from OCN  component
      got_sst(2)=.FALSE.               ! SST from DATA component
      ifield(1)=0                      ! SST from OCN  index
      ifield(2)=0                      ! SST from DATA index
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
!  values need to convert imported fields to WRF requirements.
!
          scale     =MODELS(Iatmos)%ImportField(id)%scale_factor
          add_offset=MODELS(Iatmos)%ImportField(id)%add_offset
!
          MyFmin= MISSING_dp
          MyFmax=-MISSING_dp
!
!  Load import data into WRF component variable.
!
          FieldName=ADJUSTL(ImportNameList(ifld))
!
          SELECT CASE (TRIM(FieldName))
!
!  Sea surface temperature from OCN component (C).
!
            CASE ('sst', 'SST')
              IF (.not.allocated(ocn_sst)) THEN
                allocate ( ocn_sst(LBi:UBi,LBj:UBj) )
                ocn_sst=MISSING_dp
              END IF
              IF (.not.allocated(dat_sst)) THEN
                allocate ( dat_sst(LBi:UBi,LBj:UBj) )
                dat_sst=MISSING_dp
              END IF
              got_sst(1)=.TRUE.
              ifield(1)=ifld
              fld_name(1)=TRIM(FieldName)
              DO j=JminP,JmaxP
                DO i=IminP,ImaxP
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                    MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                    Fval=scale*ptr2d(i,j)+add_offset
                    MyFmin(2)=MIN(MyFmin(2),Fval)
                    MyFmax(2)=MAX(MyFmax(2),Fval)
                    ocn_sst(i,j)=Fval
                  END IF
                END DO
              END DO
!
!  Sea surface temperature from DATA component (C).  It is used to
!  fill values in cells not covered by the OCN component.
!
            CASE ('dsst', 'dSST')
              IF (.not.allocated(ocn_sst)) THEN
                allocate ( ocn_sst(LBi:UBi,LBj:UBj) )
                ocn_sst=MISSING_dp
              END IF
              IF (.not.allocated(dat_sst)) THEN
                allocate ( dat_sst(LBi:UBi,LBj:UBj) )
                dat_sst=MISSING_dp
              END IF
              got_sst(2)=.TRUE.
              ifield(2)=ifld
              fld_name(2)=TRIM(FieldName)
              DO j=JminP,JmaxP
                DO i=IminP,ImaxP
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                    MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                    Fval=scale*ptr2d(i,j)+add_offset
                    MyFmin(2)=MIN(MyFmin(2),Fval)
                    MyFmax(2)=MAX(MyFmax(2),Fval)
                    dat_sst(i,j)=Fval
                  END IF
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
     &                        ' wrfScale = ', scale
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
!  Load or merge sea surface temperature into WRF structure variable:
!  adom(ng)%tsea
!
      IF (ANY(got_sst)) THEN
        CALL WRF_ProcessImport (grid, model,                            &
     &                          got_sst, ifield, fld_name,              &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ocn_sst, dat_sst,                       &
     &                          rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
!
!  Deallocate local arrays.
!
      IF (allocated(ImportNameList)) deallocate (ImportNameList)
      IF (allocated(ocn_sst)) deallocate (ocn_sst)
      IF (allocated(dat_sst)) deallocate (dat_sst)
!
!  Update WRF import calls counter.
!
      IF (ImportCount.gt.0) THEN
        MODELS(Iatmos)%ImportCalls=MODELS(Iatmos)%ImportCalls+1
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_Import',              &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (/,5x,'WRF_Import - unable to find option to import: ',    &
     &        a,/,18x,'check ''Import(atmos)'' in input script: ', a)
  20  FORMAT (5x,'WRF_Import - ESMF:  importing field ''',a,'''',       &
     &        t72,a,2x,'Grid ',i2.2,/,                                  &
     &        19x,'(InpMin = ', 1p,e15.8,0p,' InpMax = ',1p,e15.8,0p,   &
     &        ')')
  30  FORMAT (19x,'(OutMin = ', 1p,e15.8,0p,' OutMax = ',1p,e15.8,0p,   &
     &        1x,a,1p,e15.8,0p,')')
  40  FORMAT ('wrf_',i2.2,'_import_',a,'_',i4.4,2('-',i2.2),'_',        &
     &        i2.2,2('.',i2.2),'.nc')

      RETURN
      END SUBROUTINE WRF_Import
!
      SUBROUTINE WRF_ProcessImport (grid, model,                        &
     &                              got, ifield, FieldName,             &
     &                              LBi, UBi, LBj, UBj, Focn, Fdat,     &
     &                              rc)
!
!=======================================================================
!                                                                      !
!  If both import fields Focn and Fdat are avaliable, it merges        !
!  its values. Otherwise, it loads available data into ouput field,    !
!  Fout. Only sea-water or sea-ice points are processed. It is         !
!  used when atmosphere and ocean grids are incongruent. The DATA      !
!  component provides values on those grid points not covered by       !
!  the OCEAN component.                                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     grid       WRF grid state structure (TYPE domain)                !
!     model      Gridded component object (TYPE ESMF_GridComp)         !
!     got        Switches indicating source and availability of        !
!                import data (logical vector):                         !
!                    got(1)      OCEAN component switch (T/F)          !
!                    got(2)      DATA  component switch (T/F)          !
!     ifield     Import field index (integer vector)                   !
!                    ifield(1)   OCEAN component field index           !
!                    ifield(2)   DATA  component field index           !
!     FieldName  Field short name (string array)                       !
!     ifld       Import field index (integer)                          !
!     LBi        I-dimension lower bound  (integer)                    !
!     UBi        I-dimension upper bound  (integer)                    !
!     LBj        J-dimension lower bound  (integer)                    !
!     UBj        J-dimension upper bound  (integer)                    !
!     Focn       Import field from ocean component (2D real array)     !
!     Fdat       Import field from data  component (2D real array)     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     rc         Return code (integer)                                 !
!                                                                      !
!=======================================================================
!
      USE module_domain, ONLY : domain
      USE strings_mod,   ONLY : lowercase
!
!  Imported variable declarations.
!
      logical, intent(in)  :: got(2)
!
      integer, intent(in)  :: ifield(2)
      integer, intent(in)  :: LBi, UBi, LBj, UBj
      integer, intent(out) :: rc
!
      real (dp), intent(in) :: Focn(LBi:UBi,LBj:UBj)
      real (dp), intent(in) :: Fdat(LBi:UBi,LBj:UBj)
!
      character (len=*), intent(in) :: FieldName(:)
!
      TYPE (domain), pointer :: grid
      TYPE (ESMF_GridComp)   :: model
!
!  Local variable declarations.
!
      logical :: DebugWrite(2) = (/ .FALSE., .FALSE. /)
!
      integer :: i, ic, is, j, ng
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: LakeValue, LandValue
      integer :: localDE, localDEcount, localPET, PETcount
      integer :: IminP, ImaxP, JminP, JmaxP
!
      real (dp) :: Fseconds, TimeInDays, Time_Current

      real (dp) :: Fval, MyFmax(3), MyFmin(3), Fmin(3), Fmax(3)
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      real (KIND(grid%sst)), pointer :: Fout(:,:) => NULL()
!
      character (len=22 )     :: Time_CurrentString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_ProcessImport"
!
      character (ESMF_MAXSTR) :: cname, fld_string, ofile
!
      TYPE (ESMF_ArraySpec)  :: arraySpec2d
      TYPE (ESMF_Clock)      :: clock
      TYPE (ESMF_Field)      :: Fmerge
      TYPE (ESMF_StaggerLoc) :: staggerLoc
      TYPE (ESMF_Time)       :: CurrentTime
      TYPE (ESMF_VM)         :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_ProcessImport',       &
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
      ng=grid%grid_id
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
!  Get current time.
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
      is=INDEX(Time_CurrentString, 'T')            ! remove 'T' in
      IF (is.gt.0) Time_CurrentString(is:is)=' '   ! ISO 8601 format
!
!-----------------------------------------------------------------------
!  Create merged field.
!-----------------------------------------------------------------------
!
!  Set a 2D floating-point array descriptor.
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
!  Create 2D merge field from the Grid and arraySpec.
!
      IF (.not.got(2).and.got(1)) THEN
        DebugWrite(1)=MODELS(Iatmos)%ImportField(ifield(1))%debug_write
        fld_string=TRIM(FieldName(1))
      ELSE IF (.not.got(1).and.got(2)) THEN
        DebugWrite(2)=MODELS(Iatmos)%ImportField(ifield(2))%debug_write
        fld_string=TRIM(FieldName(2))
      ELSE IF (got(1).and.got(2)) THEN
        DebugWrite(1)=MODELS(Iatmos)%ImportField(ifield(1))%debug_write
        DebugWrite(2)=MODELS(Iatmos)%ImportField(ifield(2))%debug_write
        fld_string=TRIM(FieldName(1))//'-'//TRIM(FieldName(2))
      END IF
      staggerLoc=ESMF_STAGGERLOC_CENTER
!
      Fmerge=ESMF_FieldCreate(MODELS(Iatmos)%grid(ng),                  &
     &                        arraySpec2d,                              &
     &                        staggerloc=staggerLoc,                    &
     &                        name=TRIM(fld_string),                    &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Get merge field pointer.
!
      CALL ESMF_FieldGet (Fmerge,                                       &
     &                    farrayPtr=ptr2d,                              &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      ptr2d=MISSING_dp
!
!-----------------------------------------------------------------------
!  Create pointer to WRF export field target.
!-----------------------------------------------------------------------
!
      SELECT CASE (lowercase(TRIM(fld_string)))
        CASE ('sst', 'dsst', 'sst-dsst', 'dsst-sst')
          Fout => grid%sst
        CASE DEFAULT
          IF (localPET.eq.0) THEN
            WRITE (cplout,10) TRIM(fld_string), TRIM(CinpName)
          END IF
          rc=ESMF_RC_NOT_FOUND
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
      END SELECT
!
!  Set patch indices range.
!
      IminP=grid%sp31
      ImaxP=grid%ep31
      JminP=grid%sp33
      JmaxP=grid%ep33
      IF (grid%ed31.eq.ImaxP) THEN
        ImaxP=ImaxP-1
      END IF
      IF (grid%ed33.eq.JmaxP) THEN
        JmaxP=JmaxP-1
      END IF
!
!-----------------------------------------------------------------------
!  Set WRF mask values:
!
!  lakemask >   0: elsewhere (land, ocean)     1: lake
!  landmask >   0: elsewhere (ocean, lakes)    1: land
!
!-----------------------------------------------------------------------
!
      LakeValue=1
      LandValue=1
!
!-----------------------------------------------------------------------
!  If only one field is available, load field into output array at
!  seawater points. Notice that Fout has the same precision as the
!  WRF variable. It can be single or double precision.
!-----------------------------------------------------------------------
!
      IF (.not.got(2).and.got(1)) THEN
        MyFmin= MISSING_dp
        MyFmax=-MISSING_dp
        DO j=JminP,JmaxP
          DO i=IminP,ImaxP
            IF ((INT(grid%landmask(i,j)).ne.LandValue).and.             &
     &          (INT(grid%lakemask(i,j)).ne.LakeValue)) THEN
              Fout(i,j)=REAL(Focn(i,j), KIND(grid%sst))
            END IF
            ptr2d(i,j)=REAL(Fout(i,j), dp)
            MyFmin(1)=MIN(MyFmin(1),Fout(i,j))
            MyFmax(1)=MAX(MyFmax(1),Fout(i,j))
          END DO
        END DO
      ELSE IF (.not.got(1).and.got(2)) THEN
        MyFmin= MISSING_dp
        MyFmax=-MISSING_dp
        DO j=JminP,JmaxP
          DO i=IminP,ImaxP
            IF ((INT(grid%landmask(i,j)).ne.LandValue).and.             &
     &          (INT(grid%lakemask(i,j)).ne.LakeValue)) THEN
              Fout(i,j)=REAL(Fdat(i,j), KIND(grid%sst))
            END IF
            ptr2d(i,j)=REAL(Fout(i,j), dp)
            MyFmin(1)=MIN(MyFmin(1),Fout(i,j))
            MyFmax(1)=MAX(MyFmax(1),Fout(i,j))
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Otherwise, merge imported fields.
!-----------------------------------------------------------------------
!
      IF (got(1).and.got(2)) THEN
!
!  Merge Focn and Fdat at sea-water and sea-ice points.  Notice that
!  the ESMF regridding will not fill unbounded interpolation points.
!  Such grid cells still have the pointer initialized value MISSING_dp.
!  The TOL_dp is used to identify such values. The user has full
!  control of how the merging is done from the weights coefficients
!  provided from input NetCDF file specified in "WeightsFile(atmos)".
!
        MyFmin= MISSING_dp
        MyFmax=-MISSING_dp
        DO j=JminP,JmaxP
          DO i=IminP,ImaxP
            IF ((INT(grid%landmask(i,j)).ne.LandValue).and.             &
     &          (INT(grid%lakemask(i,j)).ne.LakeValue)) THEN
              IF (ABS(Fdat(i,j)).lt.TOL_dp) THEN
                MyFmin(2)=MIN(MyFmin(2),Fdat(i,j))
                MyFmax(2)=MAX(MyFmax(2),Fdat(i,j))
                Fval=Fdat(i,j)                   ! initialize with DATA
                IF (ABS(Focn(i,j)).lt.TOL_dp) THEN
                  MyFmin(1)=MIN(MyFmin(1),Focn(i,j))
                  MyFmax(1)=MAX(MyFmax(1),Focn(i,j))
                  Fval=WEIGHTS(Iatmos)%Cdat(i,j)*Fval+                  &
     &                 WEIGHTS(Iatmos)%Cesm(i,j)*Focn(i,j)
                END IF
                Fout(i,j)=REAL(Fval, KIND(grid%sst))
                ptr2d(i,j)=REAL(Fval, dp)
                MyFmin(3)=MIN(MyFmin(3),Fval)
                MyFmax(3)=MAX(MyFmax(3),Fval)
              END IF
            END IF
          END DO
        END DO
      END IF
!
!  Get merged fields minimun and maximum values.
!
      IF (got(1).and.got(2)) THEN
        ic=3
      ELSE
        ic=1
      END IF
      CALL ESMF_VMAllReduce (vm,                                        &
     &                       sendData=MyFmin,                           &
     &                       recvData=Fmin,                             &
     &                       count=ic,                                  &
     &                       reduceflag=ESMF_REDUCE_MIN,                &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_VMAllReduce (vm,                                        &
     &                       sendData=MyFmax,                           &
     &                       recvData=Fmax,                             &
     &                       count=ic,                                  &
     &                       reduceflag=ESMF_REDUCE_MAX,                &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Report merged import field information.
!
      IF (got(1).and.got(2)) THEN
        IF ((DebugLevel.ge.0).and.(localPET.eq.0)) THEN
          WRITE (cplout,20) TRIM(fld_string),                           &
     &                      TRIM(Time_CurrentString), ng,               &
     &                      Fmin(1), Fmax(1),                           &
     &                      Fmin(2), Fmax(2),                           &
     &                      Fmin(3), Fmax(3)
        END IF
      ELSE
        IF ((DebugLevel.ge.0).and.(localPET.eq.0)) THEN
          WRITE (cplout,30) Fmin(1), Fmax(1)
        END IF
      END IF
!
!  Debugging: write out export field into a NetCDF file.
!
      IF ((DebugLevel.ge.3).and.ANY(DebugWrite)) THEN
        WRITE (ofile,40) ng, TRIM(fld_string),                          &
     &                   year, month, day, hour, minutes, seconds
        CALL ESMF_FieldWrite (Fmerge,                                   &
     &                        TRIM(ofile),                              &
     &                        overwrite=.TRUE.,                         &
     &                        rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
!
!  Nullify pointer to make sure that it does not point on a random
!  part in the memory.
!
      IF (associated(ptr2d)) nullify (ptr2d)
      IF (associated(Fout )) nullify (Fout)
!
!  Destroy merged field.
!
      CALL ESMF_FieldDestroy (Fmerge,                                   &
     &                        noGarbage=.FALSE.,                        &
     &                        rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_ProcessImport',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (/,5x,'WRF_ProcessImport - ',                              &
     &         'unable to find option to import: ',a,                   &
     &        /,25x,'check ''Import(atmos)'' in input script: ',a)
  20  FORMAT (1x,' WRF_ProcessImport - ESMF: merging field ''',a,'''',  &
     &        t72,a,2x,'Grid ',i2.2,                                    &
     &        /,19x,'(OcnMin = ', 1p,e15.8,0p,                          &
     &              ' OcnMax = ', 1p,e15.8,0p,')',                      &
     &        /,19x,'(DatMin = ', 1p,e15.8,0p,                          &
     &              ' DatMax = ', 1p,e15.8,0p,')',                      &
     &        /,19x,'(OutMin = ', 1p,e15.8,0p,                          &
     &              ' OutMax = ', 1p,e15.8,0p,')')
  30  FORMAT (19x,  '(OutMin = ', 1p,e15.8,0p,                          &
     &              ' OutMax = ', 1p,e15.8,0p,') WRF_ProcessImport')
  40  FORMAT ('wrf_',i2.2,'_merged_',a,'_',i4.4,2('-',i2.2),'_',        &
     &         i2.2,2('.',i2.2),'.nc')
!
      RETURN
      END SUBROUTINE WRF_ProcessImport
!
      SUBROUTINE WRF_Export (grid, model, rc)
!
!=======================================================================
!                                                                      !
!  Exports WRF fields to other coupled gridded components.             !
!                                                                      !
!=======================================================================
!
      USE module_domain, ONLY : domain
# ifdef WRF_TIMEAVG
      USE strings_mod,   ONLY : uppercase
# endif
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (domain), intent(in) :: grid
      TYPE (ESMF_GridComp)      :: model
!
!  Local variable declarations.
!
      integer :: ifld, i, is, j, ng
      integer :: Istr, Iend, Jstr, Jend
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: ExportCount
      integer :: localDE, localDEcount, localPET, PETcount
# ifdef WRF_TIMEAVG
      integer :: mean_interval
# endif
!
      real (dp), parameter :: eps = 1.0E-10_dp
      real (dp), parameter :: StBolt = 5.67051E-8_dp  ! Stefan-Boltzmann
      real (dp), parameter :: z1 = 3.0_dp             ! layer thickness
!
      real (dp) :: Fseconds, TimeInDays, Time_Current
      real (dp) :: cff1, cff2, cff3, f1, scale
      real (dp) :: MyFmax(1), MyFmin(1), Fmin(1), Fmax(1), Fval
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      character (len=22)      :: Time_CurrentString
# ifdef WRF_TIMEAVG
      character (len=35)      :: Istring
# endif
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WRF_Export"
!
      character (ESMF_MAXSTR) :: cname, ofile
      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
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
        WRITE (trac,'(a,a,i0)') '==> Entering WRF_Export',              &
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
     &                       exportState=ExportState,                   &
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
      ng=grid%grid_id
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

# ifdef WRF_TIMEAVG
!
!-----------------------------------------------------------------------
!  If time-averaged WRF fluxes, check if the fields are time-averaged
!  over the same coupling window. For this option, the user needs the
!  RASM Climate Diagnostics option in input namelist (time_control):
!
!    mean_diag               flag to turn on the mean diagnostic output
!                            (1 = on, 0 = off)
!
!  then, pick only one of the following time interval parameters for
!  the desired units:
!
!    mean_diag_interval      minutes
!    mean_diag_interval_s    seconds
!    mean_diag_interval_m    minutes
!    mean_diag_interval_h    hours
!    mean_diag_interval_d    days
!    mean_diag_interval_mo   months
!
!  See WRF run/README.rams_diag for details and Registry/registry.diags
!  for full listing of variables.
!-----------------------------------------------------------------------
!
!  Check if time-averaged fluxes were activate in WRF.
!
      IF (grid%mean_diag.ne.1) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,10) 'namelist &time_control, mean_diag = ',     &
     &                      grid%mean_diag, uppercase('wrf_timeavg')
        END IF
        rc=ESMF_RC_NOT_VALID
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
!  Check time average interval (seconds) against coupling interval
!  (seconds)
!
      ELSE
        IF (grid%mean_diag_interval.gt.0) THEN
          Istring='namelist: mean_diag_interval = '
          mean_interval=grid%mean_diag_interval*60
        ELSE IF (grid%mean_diag_interval_s.gt.0) THEN
          Istring='namelist: mean_diag_interval_s = '
          mean_interval=grid%mean_diag_interval_s
        ELSE IF (grid%mean_diag_interval_m.gt.0) THEN
          Istring='namelist: mean_diag_interval_m = '
          mean_interval=grid%mean_diag_interval_m*60
        ELSE IF (grid%mean_diag_interval_h.gt.0) THEN
          Istring='namelist: mean_diag_interval_h = '
          mean_interval=grid%mean_diag_interval_h*3600
        ELSE IF (grid%mean_diag_interval_d.gt.0) THEN
          Istring='namelist: mean_diag_interval_d = '
          mean_interval=grid%mean_diag_interval_d*86400
        ELSE IF (grid%mean_diag_interval_mo.gt.0) THEN
          Istring='namelist: mean_diag_interval_mo = '
          mean_interval=grid%mean_diag_interval_mo*30*86400
        END IF
!
        IF (mean_interval.ne.INT(ClockInfo(Iatmos)%Time_Step)) THEN
          IF (localPET.eq.0) THEN
            WRITE (cplout,20) TRIM(Istring),                            &
     &                        mean_interval,                            &
     &                        TRIM(CinpName),                           &
     &                        INT(ClockInfo(Iatmos)%Time_Step)
          END IF
          rc=ESMF_RC_VAL_WRONG
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
        END IF
      END IF
# endif
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
          Istr=LBOUND(ptr2d,1)        ! iminf(ng)
          Iend=UBOUND(ptr2d,1)        ! imaxf(ng)
          Jstr=LBOUND(ptr2d,2)        ! jminf(ng)
          Jend=UBOUND(ptr2d,2)        ! jmaxf(ng)
!
!  Initialize pointer.
!
          ptr2d=MISSING_dp
!
!  Load field data into export state.  Notice that all export fields
!  are kept as computed by WRF. The imported component does the
!  proper scaling, physical units conversion, and other manipulations.
!  It is done to avoid applying such transformations twice.
!
!  Some of the computation for export fields are adapted from COAWST.
!
          SELECT CASE (TRIM(ADJUSTL(ExportNameList(ifld))))
!
!  Surface atmospheric pressure (Pa).  Use the hyposometric equation
!  to reduce surface pressure to mean sea level pressure.
!
            CASE ('psfc', 'Pair')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=REAL(grid%psfc_mean(i,j),dp)*                    &
     &                 EXP((9.81_dp*REAL(grid%ht(i,j),dp))/             &
     &                     (287.0_dp*REAL(grid%t2_mean(i,j),dp)*        &
     &                     (1.0_dp+0.61_dp*REAL(grid%q2_mean(i,j),dp))))
# else
                  Fval=REAL(grid%psfc(i,j),dp)*                         &
     &                 EXP((9.81_dp*REAL(grid%ht(i,j),dp))/             &
     &                     (287.0_dp*REAL(grid%t2(i,j),dp)*             &
     &                     (1.0_dp+0.61_dp*REAL(grid%q2(i,j),dp))))
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (2m) air temperature (K).
!
            CASE ('tsfc', 'Tair')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=REAL(grid%t2_mean(i,j),dp)
# else
                  Fval=REAL(grid%t2(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (2m) specific humidity (kg/kg). Compute the specific humidity
!  from the water vapor mixing ratio (q2) and atmospheric pressure at
!  2m. Assume the surface temperature to be equal to that at 2m (t2)
!  when computing the atmospheric pressuare at 2m using the hypsometry
!  equation.
!
            CASE ('Hair')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  cff1=REAL(grid%psfc_mean(i,j),dp)/                    &
     &                 (EXP((9.81_dp*2.0_dp)/                           &
     &                      (287.0_dp*REAL(grid%t2_mean(i,j),dp))))
                  Fval=REAL(grid%q2_mean(i,j),dp)*cff1/                 &
     &                 (REAL(grid%q2_mean(i,j),dp)*                     &
     &                  (1.0_dp-0.622_dp)+0.622_dp)
# else
                  cff1=REAL(grid%psfc(i,j),dp)/                         &
     &                 (EXP((9.81_dp*2.0_dp)/                           &
     &                      (287.0_dp*REAL(grid%t2(i,j),dp))))
                  Fval=REAL(grid%q2(i,j),dp)*cff1/                      &
     &                 (REAL(grid%q2(i,j),dp)*                          &
     &                  (1.0_dp-0.622_dp)+0.622_dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (2m) relative humidity (percentage). Compute the specific
!  humidity from the water vapor mixing ratio (q2) at 2m, as documented
!  above. Then, compute the saturation specific humidity using Bolton
!  equation.
!
            CASE ('qsfc', 'Qair')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  cff1=REAL(grid%psfc_mean(i,j),dp)/                    &
     &                 (EXP((9.81_dp*2.0_dp)/                           &
     &                      (287.0_dp*REAL(grid%t2_mean(i,j),dp))))
                  cff2=REAL(grid%q2_mean(i,j),dp)*cff1/                 &
     &                 (REAL(grid%q2_mean(i,j),dp)*                     &
     &                  (1.0_dp-0.622_dp)+0.622_dp)
                  cff3=6.112_dp*                                        &
     &                 EXP((17.67_dp*(REAL(grid%t2_mean(i,j),dp)-       &
     &                                273.15_dp))/                      &
     &                     ((REAL(grid%t2_mean(i,j),dp)-273.15_dp)+     &
     &                      243.5_dp))
# else
                  cff1=REAL(grid%psfc(i,j),dp)/                         &
     &                 (EXP((9.81_dp*2.0_dp)/                           &
     &                      (287.0_dp*REAL(grid%t2(i,j),dp))))
                  cff2=REAL(grid%q2(i,j),dp)*cff1/                      &
     &                 (REAL(grid%q2(i,j),dp)*                          &
     &                  (1.0_dp-0.622_dp)+0.622_dp)
                  cff3=6.112_dp*                                        &
     &                 EXP((17.67_dp*(REAL(grid%t2(i,j),dp)-            &
     &                                273.15_dp))/                      &
     &                     ((REAL(grid%t2(i,j),dp)-273.15_dp)+          &
     &                      243.5_dp))
# endif
                  Fval=cff2/cff3
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Net heat flux (W/m2) at the surface. Use net shortwave (downward
!  minus upward), net longwave (downward minus upward), latent, and
!  sensible fluxes to compute the surface net heat flux for the
!  ocean.
!
!  The documentation of the heat flux components is very confusing in
!  WRF. We need a consistent an unambiguous metadata model in WRF.
!
!  In WRF routine 'sst_skin_update', there is an f1=0.6352 to represent
!  the shortwave flux mean absorption in the cool-skin layer (a kludge),
!  but this value is king of high.  A formal approach is presented in
!  Zeng and Beljaars (2005; GRL).  Also, ROMS 'bulk_flux' routine
!  shows a formal cool skin correction.
!
!  In WRF, we have downwelling shortwave flux (swdnb, downward) and
!  upwelling shortwave flux (swupb, upward) at the atmosphere bottom.
!  Both are positive. The net shortwave flux (gsw = swdnb - swupb) is
!  positive and downward since 'swupb' is around 8 percent the magnitude
!  of 'swdnb'. Notice that swdown = swdnb, what a mess. Here, we use
!  the 'swdnb' and 'swupd' variables since they are available in the
!  RAMS time-averaged diagnotics.
!
!  Similarly, in WRF we have downwelling longwave flux (lwdnb, downward)
!  and upwelling longwave flux (lwupb, upward). The 'lwupb' variable is
!  almost identical to StefBo * emiss * sst^4. The 'sst' the is botton
!  land/ocean temperature (in Kelvin). Also, glw = lwdnb. Therefore,
!  'glw' is used here because it time-average RAMS diagonostic variable
!  is 'glw_mean' (downward, positive, W/m2). The metadata nomenclature
!  in WRF is such a mess and confusing.
!
!  In WRF, 'hfx' is the sensible heat flux (W/m2) and 'lh' is the
!  latent heat flux (W/m2). Both fluxes are positive if upward
!  from the surface. The same applies for the time-averaged values,
!  'hfx_mean' and 'lh_mean'. Therefore, we need to flip the sign and
!  use MINUS when adding such components to the ocean net heat flux
!  below. WARNING: The documentation for the direction of these fluxes
!  is very confusing in WRF (like negative upward).
!
            CASE ('nflx', 'shflux')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
!!            f1=1.0_dp-0.27_dp*EXP(-2.80_dp*z1)-                       &
!!   &                  0.45_dp*EXP(-0.07_dp*z1)
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=(REAL(grid%swdnb_mean(i,j),dp)-                  &
                        REAL(grid%swupb_mean(i,j),dp))+                 &
     &                 (REAL(grid%glw_mean(i,j),dp)-                    &
     &                  REAL(grid%lwupb_mean(i,j),dp))-                 &
     &                 REAL(grid%lh_mean (i,j),dp)-                     &
     &                 REAL(grid%hfx_mean(i,j),dp)
# else
                  Fval=(REAL(grid%swdnb(i,j),dp)-                       &
     &                  REAL(grid%swupb(i,j),dp))+                      &
     &                 (REAL(grid%glw(i,j),dp)-                         &
     &                  REAL(grid%lwupb(i,j),dp))-                      &
     &                 REAL(grid%lh (i,j),dp)-                          &
     &                 REAL(grid%hfx(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Net longwave radiation (W m-2) at the surface: downweling minus
!  upwelling fluxes at the bottom of the atmosphere.
!
            CASE ('lwrd', 'LWrad')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=REAL(grid%glw_mean(i,j),dp)-                     &
     &                 REAL(grid%lwupb_mean(i,j),dp)
# else
                  Fval=REAL(grid%glw(i,j),dp)-                          &
     &                 REAL(grid%lwupb(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Downward longwave radiation flux (W m-2) at surface.
!
            CASE ('dlwrd', 'dLWrad', 'lwrad_down')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=REAL(grid%glw_mean(i,j),dp)
# else
                  Fval=REAL(grid%glw(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Net shortwave radiation flux (W m-2) at the surface: downweling
!  minus  upwelling fluxes at the bottom of the atmosphere.
!
            CASE ('swrd', 'SWrad')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=REAL(grid%swdnb_mean(i,j),dp)-                   &
     &                 REAL(grid%swupb_mean(i,j),dp)
# else
                  Fval=REAL(grid%swdnb(i,j),dp)-                        &
     &                 REAL(grid%swupb(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Downward shortwave radiation flux (W m-2) at the surface.
!
            CASE ('dswrd', 'dSWrad')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=REAL(grid%swdnb_mean(i,j),dp)
# else
                  Fval=REAL(grid%swdnb(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Latent heat flux (W m-2) at surface.
!
            CASE ('lhfx', 'LHfx')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
# ifndef BULK_FLUXES
              scale=-1.0_dp           ! upward positive flux, flip sign
# else
              scale=1.0_dp
# endif
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=scale*REAL(grid%lh_mean(i,j),dp)
# else
                  Fval=scale*REAL(grid%lh(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface sensible heat flux (W m-2).
!
            CASE ('shfx', 'SHfx')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
# ifndef BULK_FLUXES
              scale=-1.0_dp           ! upward positive flux, flip sign
# else
              scale=1.0_dp
# endif
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=scale*REAL(grid%hfx_mean(i,j),dp)
# else
                  Fval=scale*REAL(grid%hfx(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface moisture (E-P) flux (kg m-2 s-1; downward).
!
            CASE ('swflux')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=REAL(grid%qfx(i,j),dp)-                          &
     &                 (REAL(grid%raincv(i,j),dp)+                      &
     &                  REAL(grid%rainncv(i,j),dp))/REAL(grid%dt,dp)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Total precipitation rate (kg m-2 s-1).
!
            CASE ('rain')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=(REAL(grid%raincv(i,j),dp)+                      &
     &                  REAL(grid%rainncv(i,j),dp))/REAL(grid%dt,dp)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Upward moisture flux at surface (kg m-2 s-1) at the surface,
!  evaporation rate.
!
            CASE ('evap')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=REAL(grid%qfx(i,j),dp)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Cloud fraction (unitless, range 0.0 to 1.0).
!
            CASE ('cloud')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=REAL(grid%cldfra(i,1,j),dp)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface eastward wind stress component (Pa). Here grid%alt is inverse
!  density (m3 kg-1), grid%u_2 and grid%v_2 are the wind components
!  (m s-1) at time level 2 rotated to geographical EAST and NORTH
!  components, and grid%ust is frictional wind magnitude
!  (m s-1) in similarity theory (Ustar).
!
            CASE ('taux', 'taux10', 'sustr')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  cff1=1.0_dp/(REAL(grid%alt(i,1,j),dp)+eps)
                  cff2=1.0_dp/                                          &
     &                 (SQRT((0.5_dp*                                   &
     &                        (REAL(grid%u_2(i  ,1,j  ),dp)+            &
     &                         REAL(grid%u_2(i+1,1,j  ),dp)))**2+       &
     &                       (0.5_dp*                                   &
     &                        (REAL(grid%v_2(i  ,1,j  ),dp)+            &
     &                         REAL(grid%v_2(i  ,1,j+1),dp)))**2)+      &
     &                  eps)
                  cff3=0.5_dp*((REAL(grid%u_2(i  ,1,j  ),dp)+           &
     &                          REAL(grid%u_2(i+1,1,j  ),dp))*          &
     &                         REAL(grid%cosa(i,j),dp)-                 &
     &                         (REAL(grid%v_2(i  ,1,j  ),dp)+           &
     &                          REAL(grid%v_2(i  ,1,j+1),dp))*          &
     &                         REAL(grid%sina(i,j),dp))
                  Fval=cff1*cff2*(REAL(grid%ust(i,j),dp)**2)*cff3
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface northward wind stress component (Pa).
!
            CASE ('tauy', 'tauy10', 'svstr')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  cff1=1.0_dp/(REAL(grid%alt(i,1,j),dp)+eps)
                  cff2=1.0_dp/                                          &
     &                 (SQRT((0.5_dp*                                   &
     &                        (REAL(grid%u_2(i  ,1,j),dp)+              &
     &                         REAL(grid%u_2(i+1,1,j),dp)))**2+         &
     &                       (0.5_dp*                                   &
     &                        (REAL(grid%v_2(i,1,j  ),dp)+              &
     &                         REAL(grid%v_2(i,1,j+1),dp)))**2)+        &
     &                  eps)
                  cff3=0.5_dp*((REAL(grid%v_2(i,1,j  ),dp)+             &
     &                          REAL(grid%v_2(i,1,j+1),dp))*            &
     &                         REAL(grid%cosa(i,j),dp)+                 &
     &                         (REAL(grid%u_2(i  ,1,j),dp)+             &
     &                          REAL(grid%u_2(i+1,1,j),dp))*            &
     &                         REAL(grid%sina(i,j),dp))
                  Fval=cff1*cff2*(REAL(grid%ust(i,j),dp)**2)*cff3
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface air density (km m-3) at grid center.
!
            CASE ('RhoAir')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=1.0_dp/(REAL(grid%alt(i,1,j),dp)+eps)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Eastward wind component (m s-1) at surface boundary layer (model
!  level 2) rotated to geographical EAST at grid center (RHO-center).
!
            CASE ('Uwind_sbl', 'u_2')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=0.5_dp*((REAL(grid%u_2(i  ,1,j  ),dp)+           &
     &                          REAL(grid%u_2(i+1,1,j  ),dp))*          &
     &                         REAL(grid%cosa(i,j),dp)-                 &
     &                         (REAL(grid%v_2(i  ,1,j  ),dp)+           &
     &                          REAL(grid%v_2(i  ,1,j+1),dp))*          &
     &                         REAL(grid%sina(i,j),dp))
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Northward wind component (m s-1) at surface boundary layer (model
!  level 2) rotated to geographical EAST at grid center (RHO-center).
!
            CASE ('Vwind_sbl', 'v_2')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=0.5_dp*((REAL(grid%v_2(i,1,j  ),dp)+             &
     &                          REAL(grid%v_2(i,1,j+1),dp))*            &
     &                         REAL(grid%cosa(i,j),dp)+                 &
     &                         (REAL(grid%u_2(i  ,1,j),dp)+             &
     &                          REAL(grid%u_2(i+1,1,j),dp))*            &
     &                         REAL(grid%sina(i,j),dp))
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (10m) eastward wind component (m s-1) rotated to
!  geographical EAST.
!
            CASE ('Uwind', 'u10', 'wndu')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=REAL(grid%u10_mean(i,j),dp)*                     &
     &                 REAL(grid%cosa(i,j),dp)-                         &
     &                 REAL(grid%v10_mean(i,j),dp)*                     &
     &                 REAL(grid%sina(i,j),dp)
# else
                  Fval=REAL(grid%u10(i,j),dp)*                          &
     &                 REAL(grid%cosa(i,j),dp)-                         &
     &                 REAL(grid%v10(i,j),dp)*                          &
     &                 REAL(grid%sina(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (10m) northward wind component (m s-1) rotated to
!  geographical NORTH.
!
            CASE ('Vwind', 'v10', 'wndv')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
# ifdef WRF_TIMEAVG
                  Fval=REAL(grid%v10_mean(i,j),dp)*                     &
     &                 REAL(grid%cosa(i,j),dp)+                         &
     &                 REAL(grid%u10_mean(i,j),dp)*                     &
     &                 REAL(grid%sina(i,j),dp)
# else
                  Fval=REAL(grid%v10(i,j),dp)*                          &
     &                 REAL(grid%cosa(i,j),dp)+                         &
     &                 REAL(grid%u10(i,j),dp)*                          &
     &                 REAL(grid%sina(i,j),dp)
# endif
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface frictional wind magnitude (m s-1) from similarity theory
!  at grid center.
!
            CASE ('Wstar')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=REAL(grid%ust(i,j),dp)
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
                WRITE (cplout,30) TRIM(ADJUSTL(ExportNameList(ifld))),  &
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
          WRITE (cplout,40) TRIM(ExportNameList(ifld)),                 &
     &                      TRIM(Time_CurrentString), ng,               &
     &                      Fmin(1), Fmax(1)
        END IF
!
!  Debugging: write out export field into a NetCDF file.
!
        IF ((DebugLevel.ge.3).and.                                      &
     &      MODELS(Iatmos)%ExportField(ifld)%debug_write) THEN
          WRITE (ofile,50) ng, TRIM(ExportNameList(ifld)),              &
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
!  Update WRF export calls counter.
!
      IF (ExportCount.gt.0) THEN
        MODELS(Iatmos)%ExportCalls=MODELS(Iatmos)%ExportCalls+1
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  WRF_Export',              &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
# ifdef WRF_TIMEAVG
  10  FORMAT (/,5x,'WRF_Export - illegal configuration: ',a,            &
     &        /,18x,a,' CPP option requires ''mean_diag = 1'' in ',     &
     &        'input ''namelist''',/,18x,'for time-averaged fluxes.')
  20  FORMAT (/,5x,'WRF_Export - inconsistent input parameters:',       &
     &        /,18x,a,1x,i0,/,18x,a,': TimeStep = ',i0)
# endif
  30  FORMAT (/,5x,'WRF_Export - unable to find option to export: ',    &
     &        a,/,18x,'check ''Export(atmos)'' in input script: ',a)
  40  FORMAT (5x,'WRF_Export - ESMF: exporting field ''',a,'''',        &
     &        t72,a,2x,'Grid ',i2.2,/,                                  &
     &        19x,'(OutMin = ', 1p,e15.8,0p,' OutMax = ',1p,e15.8,0p,   &
     &        ')')
  50  FORMAT ('wrf_',i2.2,'_export_',a,'_',i4.4,2('-',i2.2),'_',        &
     &         i2.2,2('.',i2.2),'.nc')
!
      RETURN
      END SUBROUTINE WRF_Export
!
#endif
      END MODULE esmf_wrf_mod
