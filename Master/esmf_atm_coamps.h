      MODULE esmf_coamps_mod

#if defined COAMPS_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_atm_coamps.h 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets COAMPS as the atmospheric model gridded component  !
!  using the generic ESMF/NUOPC layer:                                 !
!                                                                      !
!    ATM_SetServices         Sets ATM component shared-object entry    !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                                                                      !
!    COAMPS_SetInitializeP1  COAMPS component phase 1 initialization:  !
!                            sets import and export fields long and    !
!                            short names into its respective state.    !
!                                                                      !
!    COAMPS_SetInitializeP2  COAMPS component phase 2 initialization:  !
!                            Initializes component (COAMPS_Initialize),!
!                            sets component grid (COAMPS_SetGridArrays)!
!                            and adds fields into import and export    !
!                            into respective states (COAMPS_SetStates).!
!                                                                      !
!    COAMPS_DataInit         Exports COAMPS component fields during    !
!                            initialization or restart.                !
!                                                                      !
!    COAMPS_SetClock         Sets COAMPS component date calendar,      !
!                            start and stop times, and coupling        !
!                            interval.                                 !
# ifdef ESM_SETRUNCLOCK
!                                                                      !
!    COAMPS_SetRunClock      Sets COAMPS run clock manually.           !
# endif
!                                                                      !
!    COAMPS_CheckImport      Checks if ROMS component import field is  !
!                            at the correct time.                      !
!                                                                      !
!    COAMPS_SetGridArrays    Sets COAMPS component horizontal grid     !
!                            arrays, grid area, and land/sea mask.     !
!                                                                      !
!    COAMPS_SetStates        Adds COAMPS component export and import   !
!                            fields into its respective state.         !
!                                                                      !
!    COAMPS_ModelAdvance     Advances COAMPS component for a coupling  !
!                            interval. It calls import and export      !
!                            routines.                                 !
!                                                                      !
!    COAMPS_SetFinalize      Finalizes COAMPS component execution.     !
!                                                                      !
!    COAMPS_Import           Imports fields into COAMPS from other     !
!                            gridded components.                       !
!                                                                      !
!    COAMPS_ProcessImport    Loads or merges import COAMPS fields.     !
!                                                                      !
!    COAMPS_Export           Exports COAMPS fields to other gridded    !
!                            components.                               !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!  COAMPS: Coupled Ocean-Atmosphere Mesoscale Prediction System        !
!            https://www.nrlmry.navy.mil/coamps-web/web/home           !
!                                                                      !
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
      USE atmos_forecast, ONLY : COAMPS_Initialize => atmos_init,       &
     &                           COAMPS_Run        => atmos_run,        &
     &                           COAMPS_Finalize   => atmos_finalize
!
      implicit none
!
      PUBLIC  :: ATM_SetServices

      PRIVATE :: COAMPS_SetInitializeP1
      PRIVATE :: COAMPS_SetInitializeP2
      PRIVATE :: COAMPS_DataInit
      PRIVATE :: COAMPS_SetClock
# ifdef ESM_SETRUNCLOCK
      PRIVATE :: COAMPS_SetRunClock
# endif
      PRIVATE :: COAMPS_CheckImport
      PRIVATE :: COAMPS_SetGridArrays
      PRIVATE :: COAMPS_SetStates
      PRIVATE :: COAMPS_ModelAdvance
      PRIVATE :: COAMPS_SetFinalize
      PRIVATE :: COAMPS_Import
      PRIVATE :: COAMPS_ProcessImport
      PRIVATE :: COAMPS_Export
!
      CONTAINS
!
      SUBROUTINE ATM_SetServices (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets COAMPS component shared-object entry points for "initialize",  !
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
     &                              userRoutine=COAMPS_SetInitializeP1, &
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
     &                              userRoutine=COAMPS_SetInitializeP2, &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Attach COAMPS component phase independent specializing methods.
!-----------------------------------------------------------------------
!
!  Set routine for export initial/restart fields.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_DataInitialize,  &
     &                           specRoutine=COAMPS_DataInit,           &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for setting COAMPS clock.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_SetClock,        &
     &                           specRoutine=COAMPS_SetClock,           &
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
     &                           specRoutine=COAMPS_SetRunClock,        &
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
     &                           specRoutine=COAMPS_CheckImport,        &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for time-stepping COAMPS component.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_Advance,         &
     &                           specRoutine=COAMPS_ModelAdvance,       &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Register COAMPS finalize routine.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompSetEntryPoint (model,                           &
     &                                 methodflag=ESMF_METHOD_FINALIZE, &
     &                                 userRoutine=COAMPS_SetFinalize,  &
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
      SUBROUTINE COAMPS_SetInitializeP1 (model,                         &
     &                                   ImportState, ExportState,      &
     &                                   clock, rc)
!
!=======================================================================
!                                                                      !
!  COAMPS component Phase 1 initialization: sets import and export     !
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
     &  __FILE__//", COAMPS_SetInitializeP1"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_SetInitializeP1',  &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set COAMPS import state and fields.
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
!  Set COAMPS export state and fields.
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_SetInitializeP1',  &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE COAMPS_SetInitializeP1
!
      SUBROUTINE COAMPS_SetInitializeP2 (model,                         &
     &                                   ImportState, ExportState,      &
     &                                   clock, rc)
!
!=======================================================================
!                                                                      !
!  COAMPS component Phase 2 initialization: Initializes COAMPS, sets   !
!  component grid, and adds import and export fields to respective     !
!  states.                                                             !
!                                                                      !
!=======================================================================
!
      USE avg_mod,      ONLY : avg_init, avg_init_fld, avg_set_ptr
      USE avg_mod,      ONLY : fld_name, navg_fields
      USE avg_mod,      ONLY : ifld_airrhm,                             &
     &                         ifld_airshm,                             &
     &                         ifld_airtmp,                             &
     &                         ifld_heaflx,                             &
     &                         ifld_lahflx,                             &
     &                         ifld_lonflx,                             &
     &                         ifld_lwdown,                             &
     &                         ifld_mstflx,                             &
     &                         ifld_sehflx,                             &
     &                         ifld_slpres,                             &
     &                         ifld_solflx,                             &
     &                         ifld_stress_u_true,                      &
     &                         ifld_stress_v_true,                      &
     &                         ifld_swdown,                             &
     &                         ifld_ttlprr,                             &
     &                         ifld_u10_true,                           &
     &                         ifld_v10_true
      USE coamm_memm,   ONLY : t_nest2d_ptr
      USE coamnl_mod,   ONLY : locean
      USE coamps_parms, ONLY : max_grids
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
      logical :: got_heaflx, got_lwdown
      logical :: ltau_0
!
      integer :: StepCount, ng
      integer :: MyComm, localPET, PETcount
      integer :: ExportCount, Findex, ifld
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", COAMPS_SetInitializeP@"
!
      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
!
      TYPE (ESMF_Time)    :: CurrentTime, StartTime
      TYPE (ESMF_VM)      :: vm
!
      TYPE (t_nest2d_ptr) :: ExportPointer(NgridsA,navg_fields)
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_SetInitializeP2',  &
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
!  Initialize COAMPS component. In nested applications, COAMPS kernel
!  will allocate and initialize all grids with a single call to
!  "COAMPS_Initialize".
!-----------------------------------------------------------------------
!
      CALL COAMPS_Initialize (MyComm, .FALSE., rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Allocate COAMPS time-averaged export fields "avg" structure in terms
!  of the number of nested grids.
!-----------------------------------------------------------------------
!
      CALL avg_init (MODELS(Iatmos)%Ngrids, .TRUE.)
!
!  Get list of export fields.
!
      NESTED_LOOP : DO ng=1,MODELS(Iatmos)%Ngrids
        IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
          CALL ESMF_StateGet (MODELS(Iatmos)%ExportState(ng),           &
     &                        itemCount=ExportCount,                    &
     &                        rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          IF (.not. allocated(ExportNameList)) THEN
            allocate ( ExportNameList(ExportCount) )
          END IF
          CALL ESMF_StateGet (MODELS(Iatmos)%ExportState(ng),           &
     &                        itemNameList=ExportNameList,              &
     &                        rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
!  Allocate time-averaged export fields pointers in "avg" structure.
!  (See coamps/src/atmos/libsrc/amlib/avg_mod.F)
!
          got_heaflx=.FALSE.
          got_lwdown=.FALSE.
          DO ifld=1,ExportCount
            SELECT CASE (TRIM(ADJUSTL(ExportNameList(ifld))))
              CASE ('psfc', 'Pair')
                Findex=ifld_slpres        ! sea level pressure
              CASE ('tsfc', 'Tair')
                Findex=ifld_airtmp        ! air temperature
              CASE ('Hair')
                Findex=ifld_airshm        ! specific humidity
              CASE ('qsfc', 'Qair')
                Findex=ifld_airrhm        ! relative humidity
              CASE ('nflx', 'shflux')
                Findex=ifld_heaflx        ! net heat flux
                got_heaflx=.TRUE.
              CASE ('lwrd', 'LWrad')
                Findex=ifld_lonflx        ! longwave flux
              CASE ('dlwrd', 'dLWrad', 'lwrad_down')
                Findex=ifld_lwdown        ! downward longwave flux
                got_lwdown=.TRUE.
              CASE ('swrd', 'SWrad')
                Findex=ifld_solflx        ! shortwave flux
              CASE ('dswrd', 'dSWrad')
                Findex=ifld_swdown        ! downward shortwave flux
              CASE ('lhfx', 'LHfx')
                Findex=ifld_lahflx        ! latent heat flux
              CASE ('shfx', 'SHfx')
                Findex=ifld_sehflx        ! sensible heat flux
              CASE ('swflx', 'swflux')
                Findex=ifld_mstflx        ! moisture (E-P) flux
              CASE ('rain')
                Findex=ifld_ttlprr        ! total precipitation rate
              CASE ('taux', 'taux10', 'sustr')
                Findex=ifld_stress_u_true ! eastward wind stress
              CASE ('tauy', 'tauy10', 'svstr')
                Findex=ifld_stress_v_true ! northward wind stress
              CASE ('Uwind', 'u10', 'wndu')
                Findex=ifld_u10_true      ! eastward wind
              CASE ('Vwind', 'v10', 'wndv')
                Findex=ifld_v10_true      ! northward wind
              CASE DEFAULT
                IF (localPET.eq.0) THEN
                  WRITE (cplout,10) TRIM(ExportNameList(ifld))
                END IF
                rc=ESMF_RC_NOT_FOUND
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
            END SELECT
            CALL avg_init_fld (ng, Findex)
            CALL avg_set_ptr (ng, Findex, ExportPointer(ng,Findex)%p)
          END DO
        END IF
        IF (allocated(ExportNameList)) deallocate (ExportNameList)
!
! If computing net heat flux, allocate time-averaged downward longwave
! radiation for export.
!
        IF (.not.got_lwdown.and.got_heaflx) THEN
          CALL avg_init_fld (ng, ifld_lwdown)
          CALL avg_set_ptr (ng, ifld_lwdown, ExportPointer(ng,Findex)%p)
        END IF
      END DO NESTED_LOOP
!
!-----------------------------------------------------------------------
!  Run COAMPS with no time-stepping to finalize the initialization.
!-----------------------------------------------------------------------
!
      ltau_0=.TRUE.
      StepCount=0
      CALL COAMPS_Run (ltau_0, StepCount)
!
!  Activate "locean" to indicate that COAMPS is part of a coupled
!  system. It implies that COAMPS is invoked from the ESMF/NUOPC
!  driver as a coupled component. It is used to compute time-averaged
!  export fields in subroutine "coamm".
!
      locean=.TRUE.
!
!-----------------------------------------------------------------------
!  Set-up grid and load coordinate data.
!-----------------------------------------------------------------------
!
      DO ng=1,MODELS(Iatmos)%Ngrids
        IF (ANY(COUPLED(Iatmos)%LinkedGrid(ng,:))) THEN
          CALL COAMPS_SetGridArrays (ng, model, localPET, rc)
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
          CALL COAMPS_SetStates (ng, model, rc)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_SetInitializeP2',  &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (/,' COAMPS_SetInitializeP2 - unable to find time-',       &
     &          'averaged index for Export Field: ',a)
!
      RETURN
      END SUBROUTINE COAMPS_SetInitializeP2
!
      SUBROUTINE COAMPS_DataInit (model, rc)
!
!=======================================================================
!                                                                      !
!  Exports COAMPS component fields during initialization or restart.   !
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
      integer :: is, ng
      integer :: localPET, PETcount, phase
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", COAMPS_DataInit"
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
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_DataInit',         &
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
!  If explicit coupling from atmosphere to ocean, export initialization
!  or restart fields.
!-----------------------------------------------------------------------
!
      IF ((CouplingType.eq.0).and.(Nexport(Iatmos).gt.0)) THEN
        DO ng=1,MODELS(Iatmos)%Ngrids
          CALL COAMPS_Export (ng, model, rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
        END DO
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_DataInit',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE COAMPS_DataInit
!
      SUBROUTINE COAMPS_SetClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets COAMPS component date calendar, start and stop time, and       !
!  coupling interval.                                                  !
!                                                                      !
!=======================================================================
!
      USE coamnl_mod, ONLY : ktaust   ! starting time (hour, min, sec)
      USE coamnl_mod, ONLY : ktauf    ! ending   time (hour, min, sec)
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
      integer :: localPET, PETcount
      integer :: TimeFrac
# ifdef REGRESS_STARTCLOCK
      integer :: RegressStartDate(7)
# endif
!
      character (len= 22) :: Calendar
# ifdef REGRESS_STARTCLOCK
      character (len= 22) :: RegressStartString
# endif
      character (len= 22) :: StartTimeString, StopTimeString
      character (len=160) :: message

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", COAMPS_SetClock"
!
      TYPE (ESMF_CalKind_Flag) :: CalType
      TYPE (ESMF_Clock)        :: clock
      TYPE (ESMF_Time)         :: StartTime
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_SetClock',         &
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
!  Create COAMPS component clock.
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
!  Set reference time. Use driver configuration values.
!
      CALL ESMF_TimeSet (ClockInfo(Iatmos)%ReferenceTime,               &
     &                   yy=ReferenceDate(1),                           &
     &                   mm=ReferenceDate(2),                           &
     &                   dd=ReferenceDate(3),                           &
     &                   h =ReferenceDate(4),                           &
     &                   m =ReferenceDate(5),                           &
     &                   s =ReferenceDate(6),                           &
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
!  Use the same as driver. A coupling interval is substracted to the
!  driver clock to properly initialize all the ESM components.
!
      ClockInfo(Iatmos)%StartTime=ClockInfo(Idriver)%StartTime
!
      CALL ESMF_TimeGet (ClockInfo(Iatmos)%StartTime,                   &
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
      CALL ESMF_TimeSet (ClockInfo(Iatmos)%StartTime,                   &
                         yy=StartDate(1),                               &
                         mm=StartDate(2),                               &
                         dd=StartDate(3),                               &
                         h =StartDate(4),                               &
                         m =StartDate(5),                               &
                         s =StartDate(6),                               &
                         calendar=ClockInfo(Iatmos)%Calendar,           &
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
      CALL ESMF_TimeSet (ClockInfo(Iatmos)%StopTime,                    &
     &                   yy=StopDate(1),                                &
     &                   mm=StopDate(2),                                &
     &                   dd=StopDate(3),                                &
     &                   h =StopDate(4),                                &
     &                   m =StopDate(5),                                &
     &                   s =StopDate(6),                                &
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
!  Create COAMPS component clock.
!-----------------------------------------------------------------------
!
      ClockInfo(Iatmos)%Name='COAMPS_clock'
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
!  Compare driver time against COAMPS component time.
!-----------------------------------------------------------------------
!
      IF (ClockInfo(Idriver)%Restarted) THEN
        StartTimeString=ClockInfo(Idriver)%Time_RestartString
      ELSE
        StartTimeString=ClockInfo(Idriver)%Time_StartString
      END IF
!
!  Report start and stop time clocks.
!
      IF (localPET.eq.0) THEN
        WRITE (cplout,10) 'COAMPS Calendar:    ',                       &
     &                    TRIM(ClockInfo(Iatmos)%CalendarString),       &
     &                    'COAMPS Start Clock: ',                       &
     &                    TRIM(ClockInfo(Iatmos)%Time_StartString),     &
     &                    'COAMPS Stop  Clock: ',                       &
     &                    TRIM(ClockInfo(Iatmos )%Time_StopString)
      END IF
!
!  Compare Driver and COAMPS clocks.
!
      IF (ClockInfo(Iatmos)%Time_StartString.ne.                        &
     &    StartTimeString) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) 'COAMPS Start Time: ',                      &
     &                      TRIM(ClockInfo(Iatmos)%Time_StartString),   &
     &                      'Driver Start Time: ',                      &
     &                      TRIM(StartTimeString),                      &
     &                      '                   are not equal!'
        END IF
        message='Driver and COAMPS start times do not match: '//        &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ClockInfo(Iatmos )%Time_StopString(1:19).ne.                  &
     &    ClockInfo(Idriver)%Time_StopString(1:19)) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) 'COAMPS Stop Time: ',                       &
     &                      TRIM(ClockInfo(Iatmos )%Time_StopString),   &
     &                      'Driver Stop Time: ',                       &
     &                      TRIM(ClockInfo(Idriver)%Time_StopString),   &
     &                      '                   are not equal!'
        END IF
        message='Driver and COAMPS stop times do not match: '//         &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (TRIM(ClockInfo(Iatmos )%CalendarString).ne.                   &
     &    TRIM(ClockInfo(Idriver)%CalendarString)) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) 'COAMPS Calendar: ',                        &
     &                      TRIM(ClockInfo(Iatmos )%CalendarString),    &
     &                      'Driver Calendar: ',                        &
     &                      TRIM(ClockInfo(Idriver)%CalendarString),    &
     &                      '                  are not equal!'
        END IF
        message='Driver and COAMPS calendars do not match: '//          &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_SetClock',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
 10   FORMAT (2x,a,2x,a/,2x,a,2x,a,/,2x,a,2x,a,/)
 20   FORMAT (/,2x,a,a,/,2x,a,a,/,2x,a)
!
      RETURN
      END SUBROUTINE COAMPS_SetClock

# ifdef ESM_SETRUNCLOCK
!
      SUBROUTINE COAMPS_SetRunClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets COAMPS run clock manually to avoid getting zero time stamps at !
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
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_SetRunClock',      &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_SetRunClock',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE COAMPS_SetRunClock
# endif
!
      SUBROUTINE COAMPS_CheckImport (model, rc)
!
!=======================================================================
!                                                                      !
!  Checks if COAMPS component import field is at the correct time.     !
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
     &  __FILE__//", COAMPS_CheckImport"
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
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_CheckImport',      &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_CheckImport',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (1x,'COAMPS_CheckImport - ',a,':',t32,'TimeStamp = ',a,    &
     &        ',  DriverTime = ',a)
!
      RETURN
      END SUBROUTINE COAMPS_CheckImport
!
      SUBROUTINE COAMPS_SetGridArrays (ng, model, localPET, rc)
!
!=======================================================================
!                                                                      !
!  Sets COAMPS component staggered, horizontal grids arrays, grid      !
!  area, and land/sea mask.                                            !
!                                                                      !
!  COAMPS Grid Decomposition:                                          !
!  ==========================                                          !
!                                                                      !
!  COAMPS global horizontal domain for nest "ng":                      !
!                                                                      !
!      full-extent area: [0:m(ng)+1, 0:n(ng)+1]                        !
!         physical area: [1:m(ng)  , 1:n(ng)  ]                        !
!    computational area: [2:m(ng)-1, 2:n(ng)-1]                        !
!                                                                      !
!  COAMPS number of horizontal subdomains is "nprdom".                 !
!                                                                      !
!  COAMPS number of horizontal subdomains along each dimension for     !
!  nest "ng" is:                                                       !
!                                                                      !
!      domdec%ndx(ng)                                                  !
!      domdec%ndy(ng)                                                  !
!                                                                      !
!  COAMPS supports only one horizontal subdomain (tile) per process    !
!  (that is, one DE per PET).                                          !
!                                                                      !
!  COAMPS physical area bounds for horizontal subdomain (tile) in      !
!  nest "ng" are:                                                      !
!                                                                      !
!      [nlimx(ng)%bp(tile) : nlimx(ng)%ep(tile),                       !
!       nlimy(ng)%bp(tile) : nlimy(ng)%ep(tile)]                       !
!                                                                      !
!  COAMPS computational area bounds for horizontal subdomain (tile) in !
!  nest "ng" are:                                                      !
!                                                                      !
!      [nlimx(ng)%b (tile) : nlimx(ng)%e (tile),                       !
!       nlimy(nn)%b (tile) : nlimy(ng)%e (tile)]                       !
!                                                                      !
!  COAMPS local horizontal subdomain (tile) area bounds for nest       !
!  "ng" are:                                                           !
!                                                                      !
!   Full-extent grid:   [iminf(ng) : imaxf(ng),                        !
!                        jminf(ng) : jmaxf(ng)]                        !
!                                                                      !
!   Physical grid:      [iminp_nest(ng) : imaxp_nest(ng),              !
!                        jminp_nest(ng) : jmaxp_nest(ng)]              !
!                                                                      !
!   Interior grid:      [imini(ng) : imaxi(ng),                        !
!                        jmini(ng) : jmaxi(ng)]                        !
!                                                                      !
!  Relationship between COAMPS array and ESMF_Array subdomain regions. !
!  The second index in the ESMF LBound/UBound arrays is the DE index,  !
!  which is always one since COAMPS only supports one DE per PET.      !
!  COAMPS array indexing is based on global grid index.  ESMF array    !
!  indexing mirrors the COAMPS array indexing (ESMF_INDEX_GLOBAL).     !
!                                                                      !
!  * ESMF Exclusive Region  <=> COAMPS Physical Area                   !
!      Partial halos for subdomains that contain the physical boundary !
!      No halos otherwise                                              !
!                                                                      !
!      Array bounds:                                                   !
!                                                                      !
!      ESMF     [exclusiveLBound(1,1) : exclusiveUBound(1,1),          !
!                exclusiveLBound(2,1) : exclusiveUBound(2,1)]          !
!                                                                      !
!      COAMPS   [iminp : imaxp,                                        !
!                jminp : jmaxp]                                        !
!                                                                      !
!  * ESMF Computational Region  <=> COAMPS Physical Area               !
!      Partial halos for subdomains that contain the physical boundary !
!      No halos otherwise                                              !
!                                                                      !
!      Array bounds:                                                   !
!                                                                      !
!      ESMF     [computationalLBound(1,1) : computationalUBound(1,1),  !
!                computationalLBound(2,1) : computationalUBound(2,1)]  !
!                                                                      !
!      COAMPS   [iminp : imaxp,                                        !
!                jminp : jmaxp]                                        !
!                                                                      !
!  * ESMF Total Region  <=> COAMPS Full Extent Area                    !
!      Full halos                                                      !
!                                                                      !
!      Array bounds:                                                   !
!                                                                      !
!      ESMF     [totalLBound(1,1) : totalUBound(1,1),                  !
!                totalLBound(2,1) : totalUBound(2,1)]                  !
!                                                                      !
!      COAMPS   [iminf : imaxf,                                        !
!                jminf : jmaxf]                                        !
!                                                                      !
!=======================================================================
!
      USE coamm_memm, ONLY : adom
      USE domdec,     ONLY : iminf, imaxf, jminf, jmaxf,                &
     &                       nlimx, nlimy, nprdom, ndx, ndy
      USE gridnl_mod, ONLY : delx, dely, m, n
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
      integer :: gtype, i, ivar, j, node, tile
      integer :: localDE, localDEcount
      integer :: LBi, UBi, LBj, UBj
      integer :: cLB(2), cUB(2), eLB(2), eUB(2), tLB(2), tUB(2)
!
      integer, allocatable :: deBlockList(:,:,:)
      integer (i4b), pointer :: ptrM(:,:) => NULL()
!
      real (dp), pointer :: ptrA(:,:) => NULL()
      real (dp), pointer :: ptrX(:,:) => NULL()
      real (dp), pointer :: ptrY(:,:) => NULL()
!
      character (len=40) :: name

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", COAMPS_SetGridArrays"
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
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_SetGridArrays',    &
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
!-----------------------------------------------------------------------
!  Set tiles lower and upper bounds for each decomposition element
!  (DE=1) in terms of global indices using the physical grid area
!  (iminp:imaxp, jminp:jmaxp; no overlap between tiles).
!-----------------------------------------------------------------------
!
      IF (.not.allocated(deBlockList)) THEN
        allocate ( deBlockList(2,2,nprdom) )
      END IF
      DO tile=1,nprdom
        deBlockList(1,1,tile)=nlimx(ng)%bp(tile)       ! iminp
        deBlockList(1,2,tile)=nlimx(ng)%ep(tile)       ! imaxp
        deBlockList(2,1,tile)=nlimy(ng)%bp(tile)       ! jminp
        deBlockList(2,2,tile)=nlimy(ng)%ep(tile)       ! jmaxp
      END DO
!
!-----------------------------------------------------------------------
!  Create ESMF DistGrid object based on model domain decomposition.
!-----------------------------------------------------------------------
!
!  A single Decomposition Element (DE) per Persistent Execution Thread
!  (PET).
!
      distGrid=ESMF_DistGridCreate(minIndex=(/ 1, 1 /),                 &
     &                             maxIndex=(/ m(ng), n(ng) /),         &
     &                             deBlockList=deBlockList,             &
     &                             rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Report COAMPS DistGrid based on model domain decomposition.
!
      IF ((localPET.eq.0).and.(DebugLevel.gt.0)) THEN
        WRITE (cplout,10) ng, TRIM(GridType(Icenter))//" Point",        &
     &                    ndx(ng), ndy(ng)
        DO node=1,nprdom
          WRITE (cplout,20) node-1, deBlockList(1,1,node),              &
     &                              deBlockList(1,2,node),              &
     &                              deBlockList(2,1,node),              &
     &                              deBlockList(2,2,node)
        END DO
      END IF
      IF (allocated(deBlockList)) deallocate (deBlockList)

# ifdef DATA_COUPLING
!
!  Read in melding weights coefficients needed by COAMPS to merge
!  imported fields from DATA and other ESM components at the specified
!  nested grid because of incongruent grids.
!
      IF ((MODELS(Idata)%IsActive).and.                                 &
     &    (ng.eq.WEIGHTS(Iatmos)%NestedGrid)) THEN
        CALL get_weights (Iatmos, m(ng), n(ng), vm, rc)
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
!  Define component grid location type: Although COAMPS is discritased
!  on an Arakawa C-grid, it exports and imports fields at the grid cell
!  center.
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
     &                                  gridEdgeLWidth=(/0,0/),         &
     &                                  gridEdgeUWidth=(/0,0/),         &
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
!  The COAMPS masking is as follows, -1: inland lake
!                                     0: sea water
!                                     1: land
!                                     2: sea ice
!                                     3: land ice
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
!  is executed once since localDEcount=1. Notice that the indices for
!  the exclusive, computational, and total regions
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
                  ptrX(i,j)=adom(ng)%aln(i,j)
                  ptrY(i,j)=adom(ng)%phi(i,j)
                  ptrM(i,j)=adom(ng)%xland(i,j)
                  ptrA(i,j)=delx(ng)*dely(ng)
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
     &                            filename="coamps_"//                  &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_SetGridArrays',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT ('COAMPS_DistGrid - Grid = ',i2.2,',',3x,'Mesh = ',a,',',  &
     &        3x,'Partition = ',i0,' x ',i0)
  20  FORMAT (18x,'node = ',i0,t32,'Istr = ',i0,t45,'Iend = ',i0,       &
     &                         t58,'Jstr = ',i0,t71,'Jend = ',i0)
!
      RETURN
      END SUBROUTINE COAMPS_SetGridArrays
!
      SUBROUTINE COAMPS_SetStates (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Adds COAMPS component export and import fields into its respective  !
!  state.                                                              !
!                                                                      !
!=======================================================================
!
      USE domdec, ONLY : iminf, imaxf, jminf, jmaxf,                    &
     &                   ndom,  nlimx, nlimy
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
      integer :: IminP, ImaxP, JminP, JmaxP
      integer :: haloLW(2), haloUW(2)
!
      real (dp), dimension(:,:), pointer :: ptr2d => NULL()
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", COAMPS_SetStates"
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
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_SetStates',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Compute lower and upper bound tile halo widths for ESMF fields.
!-----------------------------------------------------------------------
!
      IminP=nlimx(ng)%bp(ndom)
      ImaxP=nlimx(ng)%ep(ndom)
      JminP=nlimy(ng)%bp(ndom)
      JmaxP=nlimy(ng)%ep(ndom)
!
      haloLW(1)=IminP-iminf(ng)
      haloLW(2)=JminP-jminf(ng)
      haloUW(1)=imaxf(ng)-ImaxP
      haloUW(2)=jmaxf(ng)-JmaxP
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
     &                             totalLWidth=haloLW,                  &
     &                             totalUWidth=haloUW,                  &
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
     &                             totalLWidth=haloLW,                  &
     &                             totalUWidth=haloUW,                  &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_SetStates',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
!
 10   FORMAT ('COAMPS_SetStates - Removing field ''',a,''' from ',a,    &
     &        '''',a,'''',/,19x,'because it is not connected.')
!
      RETURN
      END SUBROUTINE COAMPS_SetStates
!
      SUBROUTINE COAMPS_ModelAdvance (model, rc)
!
!=======================================================================
!                                                                      !
!  Advance COAMPS component for a coupling interval (seconds) using    !
!  "COAMPS_Run". It also calls "COAMPS_Import" and "COAMPS_Export" to  !
!  import and export coupling fields, respectively.                    !
!                                                                      !
!=======================================================================
!
      USE coamnl_mod, ONLY : delta             ! timestep in seconds
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      logical :: Ladvance, ltau_0
!
      integer :: is, ng
      integer :: localPET, PETcount, phase
      integer :: NstrStep, NendStep, StepCount
!
      real (dp) :: CouplingInterval, SecondsSinceStart
      real (dp) :: TcurrentInSeconds, TstopInSeconds
!
      character (len=22) :: Cinterval
      character (len=22) :: CurrTimeString, StopTimeString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", COAMPS_SetModelAdvance"
!
      TYPE (ESMF_Clock)        :: clock
      TYPE (ESMF_State)        :: ExportState, ImportState
      TYPE (ESMF_TimeInterval) :: TimeStep
      TYPE (ESMF_Time)         :: ReferenceTime
      TYPE (ESMF_Time)         :: CurrentTime, StartTime, StopTime
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_ModelAdvance',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get information about the gridded component.
!-----------------------------------------------------------------------
!
!  Inquire about COAMPS component.
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
!  Current COAMPS time (seconds).
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
!  COAMPS stop time (seconds) for this coupling window.
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
!
!-----------------------------------------------------------------------
!  Calculate run time for the current coupling window.
!-----------------------------------------------------------------------
!
!  Get elapsed time in seconds since start.
!
      IF (ClockInfo(Idriver)%Restarted) THEN
        SecondsSinceStart=TcurrentInSeconds-                            &
     &                    ClockInfo(Iatmos)%Time_Restart
      ELSE
        SecondsSinceStart=TcurrentInSeconds-                            &
     &                    ClockInfo(Iatmos)%Time_Start
      END IF
!
!  Set number of COAMPS timesteps to run.
!
      NstrStep=INT((SecondsSinceStart+0.001_dp)/delta)+1
      NendStep=INT(SecondsSinceStart+CouplingInterval+0.001_dp)/delta
      StepCount=NendStep-NstrStep+1

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
            CALL COAMPS_Import (ng, model, rc=rc)
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
!  Run COAMPS component. Notice that atmosphere component is advanced
!  when ng=1.  In nested application, its numerical kernel will advance
!  all the nested grids in their logical order.
!-----------------------------------------------------------------------
!
      IF (Ladvance) THEN
        ltau_0=.FALSE.
        IF (ESM_track) THEN
          WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_Run',            &
     &                            ', PET', PETrank
          CALL my_flush (trac)
        END IF
        CALL COAMPS_Run (ltau_0, StepCount)
        IF (ESM_track) THEN
          WRITE (trac,'(a,a,i0)') '==> Exiting  COAMPS_Run',            &
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
            CALL COAMPS_Export (ng, model, rc=rc)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_ModelAdvance',     &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (3x,'ModelAdvance - ESMF, Running COAMPS:',t42,a,          &
     &        ' => ',a,', [',a,' s], Advance: ',l1)
!
      RETURN
      END SUBROUTINE COAMPS_ModelAdvance
!
      SUBROUTINE COAMPS_SetFinalize (model,                             &
     &                               ImportState, ExportState,          &
     &                               clock, rc)
!
!=======================================================================
!                                                                      !
!  Finalize COAMPS component execution. It calls COAMPS_Finalize.      !
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
     &  __FILE__//", COAMPS_SetFinalize"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_SetFinalize',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Finalize COAMPS component.
!-----------------------------------------------------------------------
!
      CALL COAMPS_Finalize ()
      CALL my_flush (6)                   ! flush standard output buffer
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_SetFinalize',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE COAMPS_SetFinalize
!
      SUBROUTINE COAMPS_Import (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Imports fields into COAMPS array structure from other coupled       !
!  gridded components.                                                 !
!                                                                      !
!=======================================================================
!
      USE coamm_memm, ONLY : adom
      USE domdec,     ONLY : iminf, imaxf, jminf, jmaxf,                &
     &                       ndom,  nlimx, nlimy
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
      logical :: got_sst(2)
!
      integer :: id, ifld, i, is, j
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: SeaIce, SeaWater
      integer :: ImportCount
      integer :: localDE, localDEcount, localPET, PETcount
      integer :: LBi, UBi, LBj, UBj
      integer :: IminP, ImaxP, JminP, JmaxP
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
     &  __FILE__//", COAMPS_Import"
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
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_Import',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Compute COAMPS lower and upper bounds (non-overlapping) for physical
!  area per nested grid and tile.
!-----------------------------------------------------------------------
!
      IminP=nlimx(ng)%bp(ndom)
      ImaxP=nlimx(ng)%ep(ndom)
      JminP=nlimy(ng)%bp(ndom)
      JmaxP=nlimy(ng)%ep(ndom)
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
!  If the regridding includes an extrapolation option to fill unmapped
!  grid cells due to incongruents ESM grids, the land/sea mask arras
!  are used to load only the needed data.
!
      SeaWater=0                       ! COAMPS sea water mask value
      SeaIce=2                         ! COAMPS sea ice   mask value
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
!  values need to convert imported fields to COAMPS requirements.
!
          scale     =MODELS(Iatmos)%ImportField(id)%scale_factor
          add_offset=MODELS(Iatmos)%ImportField(id)%add_offset
!
          MyFmin= MISSING_dp
          MyFmax=-MISSING_dp
!
!  Load import data into COAMPS component variable.
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
                  IF (((NINT(adom(ng)%xland(i,j)).eq.SeaWater).or.      &
     &                 (NINT(adom(ng)%xland(i,j)).eq.SeaIce)).and.      &
     &                 (ABS(ptr2d(i,j)).lt.TOL_dp)) THEN
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
                  IF (((NINT(adom(ng)%xland(i,j)).eq.SeaWater).or.      &
     &                 (NINT(adom(ng)%xland(i,j)).eq.SeaIce)).and.      &
     &                 (ABS(ptr2d(i,j)).lt.TOL_dp)) THEN
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
!  Wave-induced Charnock parameter.
!
            CASE ('charno', 'Charnock')
              DO j=JminP,JmaxP
                DO i=IminP,ImaxP
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  adom(ng)%charnock(i,j)=Fval
                END DO
              END DO
!
!  Surface, wave-induced eastward stress.
!
            CASE ('Wustr')
              DO j=JminP,JmaxP
                DO i=IminP,ImaxP
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  adom(ng)%wvsu(i,j)=Fval
                END DO
              END DO
!
!  Surface, wave-induced eastward stress.
!
            CASE ('Wvstr')
              DO j=JminP,JmaxP
                DO i=IminP,ImaxP
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  adom(ng)%wvsv(i,j)=Fval
                END DO
              END DO
!
!  Surface, wave-induced stress magnitude.
!
            CASE ('Wstr')
              DO j=JminP,JmaxP
                DO i=IminP,ImaxP
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  adom(ng)%wvst(i,j)=Fval
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
     &                        ' coampsScale = ', scale
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
!  Load or merge sea surface temperature into COAMPS structure variable:
!  adom(ng)%tsea
!
      IF (ANY(got_sst)) THEN
        CALL COAMPS_ProcessImport (ng, model,                           &
     &                             got_sst, ifield, fld_name,           &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ocn_sst, dat_sst,                    &
     &                             rc)
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
!  Update COAMPS import calls counter.
!
      IF (ImportCount.gt.0) THEN
        MODELS(Iatmos)%ImportCalls=MODELS(Iatmos)%ImportCalls+1
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_Import',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (/,2x,'COAMPS_Import - unable to find option to import: ', &
     &        a,/,18x,'check ''Import(atmos)'' in input script: ', a)
  20  FORMAT (2x,'COAMPS_Import - ESMF: importing field ''',a,'''',     &
     &        t72,a,2x,'Grid ',i2.2,/,                                  &
     &        19x,'(InpMin = ', 1p,e15.8,0p,' InpMax = ',1p,e15.8,0p,   &
     &        ')')
  30  FORMAT (19x,'(OutMin = ', 1p,e15.8,0p,' OutMax = ',1p,e15.8,0p,   &
     &        1x,a,1p,e15.8,0p,')')
  40  FORMAT ('coamps_',i2.2,'_import_',a,'_',i4.4,2('-',i2.2),'_',     &
     &        i2.2,2('.',i2.2),'.nc')

      RETURN
      END SUBROUTINE COAMPS_Import
!
      SUBROUTINE COAMPS_ProcessImport (ng, model,                       &
     &                                 got, ifield, FieldName,          &
     &                                 LBi, UBi, LBj, UBj, Focn, Fdat,  &
     &                                 rc)
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
!     ng         Nested grid number (integer)                          !
!     model      Gridded component object (TYPE ESMF_GridComp)         !
!     got        Switches indicating source and availability of        !
!                import data (logical vector):                         !
!                    got(1)      OCEAN component switch (T/F)          !
!                    got(2)      DATA  component switch (T/F)          !
!     ifield     Import field index (integer vector)                   !
!                    ifield(1)   OCEAN component field index           !
!                    ifield(2)   DATA  component field index           !
!     FieldName  Field short name (string array)                       !
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
      USE coamm_memm,  ONLY : adom
      USE domdec,      ONLY : iminf, imaxf, jminf, jmaxf,               &
     &                        ndom,  nlimx, nlimy
      USE strings_mod, ONLY : lowercase
!
!  Imported variable declarations.
!
      logical, intent(in)  :: got(2)
!
      integer, intent(in)  :: ng, ifield(2)
      integer, intent(in)  :: LBi, UBi, LBj, UBj
      integer, intent(out) :: rc
!
      real (dp), intent(in) :: Focn(LBi:UBi,LBj:UBj)
      real (dp), intent(in) :: Fdat(LBi:UBi,LBj:UBj)
!
      character (len=*), intent(in) :: FieldName(:)
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      logical :: DebugWrite(2) = (/ .FALSE., .FALSE. /)
!
      integer :: i, ic, is, j
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: SeaIce, SeaWater
      integer :: localDE, localDEcount, localPET, PETcount
      integer :: IminP, ImaxP, JminP, JmaxP
!
      real (dp) :: Fseconds, TimeInDays, Time_Current

      real (dp) :: Fval, MyFmax(3), MyFmin(3), Fmin(3), Fmax(3)
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      real (KIND(adom(1)%tsea)), pointer :: Fout(:,:) => NULL()
!
      character (len=22 )     :: Time_CurrentString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", COAMPS_ProcessImport"
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
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_ProcessImport',    &
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
!  Create pointer to COAMPS export field target. Here, adom(ng)%tsea
!  has surface surface temperature values in land and ocean points.
!  Only the ocean points are updated.
!-----------------------------------------------------------------------
!
      SELECT CASE (lowercase(TRIM(fld_string)))
        CASE ('sst', 'dsst', 'sst-dsst', 'dsst-sst')
          Fout => adom(ng)%tsea
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
!  Set COAMPS lower and upper bounds (non-overlapping) for physical
!  area per nested grid and tile.
!
      IminP=nlimx(ng)%bp(ndom)
      ImaxP=nlimx(ng)%ep(ndom)
      JminP=nlimy(ng)%bp(ndom)
      JmaxP=nlimy(ng)%ep(ndom)
!
!-----------------------------------------------------------------------
!  Set COAMPS mask value at seawater and seaice points:
!
!    -1: inland lake   0: sea water   1: land   2: sea ice   3: land ice
!-----------------------------------------------------------------------
!
      SeaWater=0
      SeaIce=2
!
!-----------------------------------------------------------------------
!  If only one field is available, load field into output array at
!  seawater points. Notice that Fout has the same precision as the
!  COAMPS variable. It can be single or double precision.
!-----------------------------------------------------------------------
!
      IF (.not.got(2).and.got(1)) THEN
        MyFmin= MISSING_dp
        MyFmax=-MISSING_dp
        DO j=JminP,JmaxP
          DO i=IminP,ImaxP
            IF (((NINT(adom(ng)%xland(i,j)).eq.SeaWater).or.            &
     &           (NINT(adom(ng)%xland(i,j)).eq.SeaIce)).and.            &
     &           (ABS(Focn(i,j)).lt.TOL_dp)) THEN
              Fout(i,j)=REAL(Focn(i,j), KIND(adom(ng)%tsea))
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
            IF (((NINT(adom(ng)%xland(i,j)).eq.SeaWater).or.            &
     &           (NINT(adom(ng)%xland(i,j)).eq.SeaIce)).and.            &
     &           (ABS(Fdat(i,j)).lt.TOL_dp)) THEN
              Fout(i,j)=REAL(Fdat(i,j), KIND(adom(1)%tsea))
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
            IF ((NINT(adom(ng)%xland(i,j)).eq.SeaWater).or.             &
     &          (NINT(adom(ng)%xland(i,j)).eq.SeaIce)) THEN
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
                Fout(i,j)=REAL(Fval, KIND(adom(ng)%tsea))
                ptr2d(i,j)=REAL(Fval, dp)
                MyFmin(3)=MIN(MyFmin(3),Fval)
                MyFmax(3)=MAX(MyFmax(3),Fval)
              END IF
            ELSE
              ptr2d(i,j)=REAL(Fout(i,j), dp)     ! include land values
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_ProcessImport',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (/,5x,'COAMPS_ProcessImport - ',                           &
     &         'unable to find option to import: ',a,                   &
     &        /,25x,'check ''Import(atmos)'' in input script: ',a)
  20  FORMAT (1x,' COAMPS_ProcessImport - ESMF merging field ''',       &
     &        a,'''',t72,a,2x,'Grid ',i2.2,                             &
     &        /,19x,'(OcnMin = ', 1p,e15.8,0p,                          &
     &              ' OcnMax = ', 1p,e15.8,0p,')',                      &
     &        /,19x,'(DatMin = ', 1p,e15.8,0p,                          &
     &              ' DatMax = ', 1p,e15.8,0p,')',                      &
     &        /,19x,'(OutMin = ', 1p,e15.8,0p,                          &
     &              ' OutMax = ', 1p,e15.8,0p,')')
  30  FORMAT (19x,  '(OutMin = ', 1p,e15.8,0p,                          &
     &              ' OutMax = ', 1p,e15.8,0p,') COAMPS_ProcessImport')
  40  FORMAT ('coamps_',i2.2,'_merged_',a,'_',i4.4,2('-',i2.2),'_',     &
     &        i2.2,2('.',i2.2),'.nc')
!
      RETURN
      END SUBROUTINE COAMPS_ProcessImport
!
      SUBROUTINE COAMPS_Export (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Exports COAMPS fields to other coupled gridded components. The      !
!  fields in COAMPS are time-averaged over the coupling interval.      !
!                                                                      !
!  The time-averaging of exported surface fields is done in COAMPS     !
!  file: ROOT_DIR/coamps/src/atmos/libsrc/amlib/avg_mod.F              !
!                                                                      !
!=======================================================================
!
      USE avg_mod,    ONLY : avg
      USE avg_mod,    ONLY : ifld_airrhm,                               &
     &                       ifld_airshm,                               &
     &                       ifld_airtmp,                               &
     &                       ifld_heaflx,                               &
     &                       ifld_lahflx,                               &
     &                       ifld_lonflx,                               &
     &                       ifld_lwdown,                               &
     &                       ifld_mstflx,                               &
     &                       ifld_sehflx,                               &
     &                       ifld_slpres,                               &
     &                       ifld_solflx,                               &
     &                       ifld_stress_u_true,                        &
     &                       ifld_stress_v_true,                        &
     &                       ifld_swdown,                               &
     &                       ifld_ttlprr,                               &
     &                       ifld_u10_true,                             &
     &                       ifld_v10_true
      USE coamm_memm, ONLY : adom
      USE domdec,     ONLY : iminf, imaxf, jminf, jmaxf
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
      integer :: ifld, i, is, j
      integer :: Istr, Iend, Jstr, Jend
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: ExportCount
      integer :: localDE, localDEcount, localPET, PETcount
!
      real (dp), parameter :: Emiss = 0.97_dp         ! IR emissivity
      real (dp), parameter :: StBolt = 5.67051E-8_dp  ! Stefan-Boltzmann
      real (dp), parameter :: z1 = 3.0_dp             ! layer thickness
!
      real (dp) :: Fseconds, TimeInDays, Time_Current
      real (dp) :: cff1, cff2, f1, scale

      real (dp) :: MyFmax(1), MyFmin(1), Fmin(1), Fmax(1), Fval
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      character (len=22)      :: Time_CurrentString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", COAMPS_Export"
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
        WRITE (trac,'(a,a,i0)') '==> Entering COAMPS_Export',           &
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
!  are kept as computed by COAMPS. The imported component does the
!  proper scaling, physical units conversion, and other manipulations.
!  It is done to avoid applying such transformations twice.
!
          SELECT CASE (TRIM(ADJUSTL(ExportNameList(ifld))))
!
!  Sea level pressure (Pa).
!
            CASE ('psfc', 'Pair')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_slpres)%p(i,j)
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
                  Fval=avg(ng)%fld(ifld_airtmp)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (2m) specific humidity (kg/kg).
!
            CASE ('Hair')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_airshm)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (2m) relative humidity (percentage).
!
            CASE ('qsfc', 'Qair')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_airrhm)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Net heat flux (W m-2) at the surface. Use shortwave, longwave,
!  latent, sensible fluxes to compute net heat flux.  Remove outgoing
!  IR from ocean sea surface temperature (K) using infrared emissivity
!  (unitless) and Stefan-Boltzmann constant (W m-2 K-4).  As in COAMPS
!  routine 'sst_skin_update', the f1 represents the shortwave flux
!  mean absorption in the cool-skin layer (an approximation kludge).
!  A formal approach is presented in Zeng and Beljaars (2005; GRL).
!  Also, ROMS 'bulk_flux' routine shows a formal cool skin correction.
!
!  The latent and sensible flux computed in COAMPS need to have the
!  sign reversed because of COAMPS convention of positive to for
!  upward flux and negative for downward flux. In the ocean is the
!  opposite.
!
            CASE ('nflx', 'shflux')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              f1=1.0_dp-0.27_dp*EXP(-2.80_dp*z1)-                       &
     &                  0.45_dp*EXP(-0.07_dp*z1)
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  cff1=adom(ng)%tsea(i,j)*adom(ng)%tsea(i,j)*           &
     &                 adom(ng)%tsea(i,j)*adom(ng)%tsea(i,j)
                  cff2=Emiss*StBolt*cff1
                  Fval=avg(ng)%fld(ifld_solflx)%p(i,j)*f1+              &
     &                 avg(ng)%fld(ifld_lwdown)%p(i,j)-cff2-            &
     &                 avg(ng)%fld(ifld_lahflx)%p(i,j)-                 &
     &                 avg(ng)%fld(ifld_sehflx)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface net longwave radiation flux (W m-2; positive upward).
!
            CASE ('lwrd', 'LWrad')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_lonflx)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface downward longwave radiation flux (W m-2).
!
            CASE ('dlwrd', 'dLWrad', 'lwrad_down')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_lwdown)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface net shortwave radiation (W m-2; positive into ocean).
!
            CASE ('swrd', 'SWrad')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_solflx)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface downward shortwave radiation flux (W m-2).
!
            CASE ('dswrd', 'dSWrad')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_swdown)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface latent heat flux (W m-2). In COAMPS, the latent heat flux
!  is a positive upward flux. For the ocean, it is the reverse and
!  needs to be switched to negative.
!
!
            CASE ('lhfx', 'LHfx')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=-avg(ng)%fld(ifld_lahflx)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface sensible heat flux (W m-2).  In COAMPS, the sensible heat
!  flux is a positive upward flux. For the ocean, it is the reverse and
!  needs to be switched to negative.
!
            CASE ('shfx', 'SHfx')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=-avg(ng)%fld(ifld_sehflx)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface moisture (E-P) flux (kg m-2 s-1). In COAMPS, the evaporation
!  is a positive upward flux. For the ocean, it is the reverse so the
!  moisture flux needs to be switched to negative.
!
            CASE ('swflux')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=-avg(ng)%fld(ifld_mstflx)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Precipitation tendency rate (kg m-2 s-1). In COAMPS, precipitation
!  is averaged with cm/s units.
!
            CASE ('rain')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              scale=10.0_dp       ! cm/s to kg m-2 s-1 (rhow=1000 km/m3)
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_ttlprr)%p(i,j)*scale
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (10m) eastward wind stress component (millibar, mb).
!
            CASE ('taux', 'taux10', 'sustr')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_stress_u_true)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (10m) northward wind stress component (Pa).
!
            CASE ('tauy', 'tauy10', 'svstr')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_stress_v_true)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (10m) eastward wind component (m s-1).
!
            CASE ('Uwind', 'u10', 'wndu')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_u10_true)%p(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Surface (10m) northward wind component (m s-1).
!
            CASE ('Vwind', 'v10', 'wndv')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Fval=avg(ng)%fld(ifld_v10_true)%p(i,j)
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
!  Update COAMPS export calls counter.
!
      IF (ExportCount.gt.0) THEN
        MODELS(Iatmos)%ExportCalls=MODELS(Iatmos)%ExportCalls+1
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  COAMPS_Export',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
  10  FORMAT (/,2x,'COAMPS_Export - unable to find option to export: ', &
     &        a,/,18x,'check ''Export(atmos)'' in input script: ',a)
  20  FORMAT (2x,'COAMPS_Export - ESMF: exporting field ''',a,'''',     &
     &        t72,a,2x,'Grid ',i2.2,/,                                  &
     &        19x,'(OutMin = ', 1p,e15.8,0p,' OutMax = ',1p,e15.8,0p,   &
     &        ')')
  30  FORMAT ('coamps_',i2.2,'_export_',a,'_',i4.4,2('-',i2.2),'_',     &
     &        i2.2,2('.',i2.2),'.nc')

      RETURN
      END SUBROUTINE COAMPS_Export
!
#endif
      END MODULE esmf_coamps_mod
