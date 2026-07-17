#include "cppdefs.h"
      MODULE esmf_roms_mod

#if defined MODEL_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_roms.F 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module sets ROMS as the ocean gridded component using generic  !
!  ESMF/NUOPC layer:                                                   !
!                                                                      !
!    ROMS_SetServices        Sets ROMS component shared-object entry   !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                                                                      !
!    ROMS_SetInitializeP1    ROMS component phase 1 initialization:    !
!                            sets import and export fields long and    !
!                            short names into its respective state.    !
!                                                                      !
!    ROMS_SetInitializeP2    ROMS component phase 2 initialization:    !
!                            Initializes component (ROMS_initialize),  !
!                            sets component grid (ROMS_SetGridArrays), !
!                            and adds fields into import and export    !
!                            into respective states.                   !
!                                                                      !
!    ROMS_DataInit           Exports ROMS component fields during      !
!                            initialization or restart.                !
!                                                                      !
!    ROMS_SetClock           Sets ROMS component date calendar, start  !
!                            and stop times, and coupling interval.    !
# ifdef ESM_SETRUNCLOCK
!                                                                      !
!    ROMS_SetRunClock        Sets ROMS run clock manually.             !
# endif
!                                                                      !
!    ROMS_CheckImport        Checks if ROMS component import field is  !
!                            at the correct time.                      !
!                                                                      !
!    ROMS_SetGridArrays      Sets ROMS component staggered, horizontal !
!                            grid arrays, grid area, and land/sea mask !
!                            if any.                                   !
!                                                                      !
!    ROMS_SetStates          Adds ROMS component export and import     !
!                            fields into its respective state.         !
!                                                                      !
!    ROMS_ModelAdvance       Advances ROMS component for a coupling    !
!                            interval. It calls import and export      !
!                            routines.                                 !
!                                                                      !
!    ROMS_SetFinalize        Finalizes ROMS component execution.       !
!                                                                      !
!    ROMS_Import             Imports fields into ROMS. The fields are  !
!                            loaded into the snapshot storage arrays   !
!                            to allow time interpolation elsewhere.    !
!                                                                      !
!    ROMS_Export             Exports ROMS fields to other gridded      !
!                            components.                               !
!                                                                      !
!    ROMS_Rotate             Rotates exchanged vector components from  !
!                            computational grid to geographical EAST   !
!                            and NORTH directions or vice versa.       !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!  ROMS:   Regional Ocean Modeling System                              !
!            https://www.myroms.org                                    !
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
!-----------------------------------------------------------------------
!  ROMS module association: parameters, variables, derived-type objects.
!-----------------------------------------------------------------------
!
      USE roms_kernel_mod,  ONLY : ROMS_initialize,                     &
     &                             ROMS_run,                            &
     &                             ROMS_finalize
!
      USE bc_2d_mod,        ONLY : bc_r2d_tile
      USE dateclock_mod,    ONLY : ROMS_clock, caldate, time_string
      USE exchange_2d_mod,  ONLY : exchange_r2d_tile,                   &
     &                             exchange_u2d_tile,                   &
     &                             exchange_v2d_tile
      USE get_metadata_mod, ONLY : CouplingField,                       &
     &                             cmeps_metadata
      USE mod_kinds,        ONLY : dp, i4b, i8b, r4, r8
      USE mod_forces,       ONLY : FORCES
      USE mod_grid,         ONLY : GRID
      USE mod_iounits,      ONLY : Iname, SourceFile, stdout
      USE mod_mixing,       ONLY : MIXING
      USE mod_ncparam,      ONLY : Iinfo,  idLdwn, idLrad, idPair,      &
     &                             idQair, idrain, idSrad, idTair,      &
     &                             idTsur, idUair, idUsms, idVair,      &
     &                             idVsms
# ifdef TIME_INTERP
      USE mod_netcdf,       ONLY : netcdf_get_ivar,                     &
     &                             netcdf_get_svar,                     &
     &                             netcdf_get_time
# endif
      USE mod_ocean,        ONLY : OCEAN
      USE mod_param,        ONLY : BOUNDS, Lm, Mm, N, NghostPoints,     &
     &                             Ngrids, NtileI, NtileJ, iNLM,        &
     &                             r2dvar, u2dvar, v2dvar
      USE mod_scalars,      ONLY : Cp, EWperiodic, NSperiodic, NoError, &
     &                             Rclock, dt, exit_flag, itemp, isalt, &
     &                             ntfirst, ntend, ntimes,              &
     &                             rho0, sec2day, tdays, time_ref
      USE mod_stepping,     ONLY : nstp, knew
      USE mp_exchange_mod,  ONLY : mp_exchange2d
      USE stdinp_mod,       ONLY : getpar_i
      USE strings_mod,      ONLY : FoundError, assign_string
!
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
      PUBLIC  :: ROMS_SetServices

      PRIVATE :: ROMS_SetInitializeP1
      PRIVATE :: ROMS_SetInitializeP2
      PRIVATE :: ROMS_DataInit
      PRIVATE :: ROMS_SetClock
# ifdef ESM_SETRUNCLOCK
      PRIVATE :: ROMS_SetRunClock
# endif
      PRIVATE :: ROMS_CheckImport
      PRIVATE :: ROMS_SetGridArrays
      PRIVATE :: ROMS_SetStates
      PRIVATE :: ROMS_ModelAdvance
      PRIVATE :: ROMS_SetFinalize
      PRIVATE :: ROMS_Import
      PRIVATE :: ROMS_Export
      PRIVATE :: ROMS_Rotate
!
      PRIVATE
!
!  Define parameters to rotate exchanged fields from geographical (EAST,
!  NORTH) to computational directions or vice versa. The resulting
!  vector components can be staggered (U- and V-points) or at cell
!  center (RHO-points: full or interior grid).
!
      integer, parameter :: geo2grid     = 0        ! U- and V-points
      integer, parameter :: geo2grid_rho = 0        ! RHO-points
      integer, parameter :: grid2geo_rho = 1        ! export vector
!
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!
      SUBROUTINE ROMS_SetServices (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets ROMS component shared-object entry points for "initialize",    !
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
     &  __FILE__//", ROMS_SetServices"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_SetServices',        &
     &                          ', PET', PETrank
        FLUSH (trac)
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
     &                              userRoutine=ROMS_SetInitializeP1,   &
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
     &                              userRoutine=ROMS_SetInitializeP2,   &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Attach ROMS component phase independent specializing methods.
!-----------------------------------------------------------------------
!
!  Set routine for export initial/restart fields.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_DataInitialize,  &
     &                           specRoutine=ROMS_DataInit,             &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for setting ROMS clock.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_SetClock,        &
     &                           specRoutine=ROMS_SetClock,             &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF

# ifdef ESM_SETRUNCLOCK
!
!  Set routine for setting ROMS run clock manually. First, remove the
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
     &                           specRoutine=ROMS_SetRunClock,          &
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
     &                           specRoutine=ROMS_CheckImport,          &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for time-stepping ROMS component.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_Advance,         &
     &                           specRoutine=ROMS_ModelAdvance,         &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Register ROMS finalize routine.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompSetEntryPoint (model,                           &
     &                                 methodflag=ESMF_METHOD_FINALIZE, &
     &                                 userRoutine=ROMS_SetFinalize,    &
     &                                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_SetServices',        &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
      RETURN
      END SUBROUTINE ROMS_SetServices
!
      SUBROUTINE ROMS_SetInitializeP1 (model,                           &
     &                                 ImportState, ExportState,        &
     &                                 clock, rc)
!
!=======================================================================
!                                                                      !
!  ROMS component Phase 1 initialization: sets import and export       !
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
      integer :: i, ng, localPET
!
      character (len=100) :: CoupledSet, StateLabel
      character (len=240) :: StandardName, ShortName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_SetInitializeP1"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_SetInitializeP1',    &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the MPI
!  current node rank.
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
!  Set ROMS Import State metadata.
!-----------------------------------------------------------------------
!
!  Add ROMS import state(s). If nesting, each grid has its own import
!  state.
!
      IMPORTING : IF (Nimport(Iroms).gt.0) THEN
        DO ng=1,MODELS(Iroms)%Ngrids
          IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
            CoupledSet=TRIM(COUPLED(Iroms)%SetLabel(ng))
            StateLabel=TRIM(COUPLED(Iroms)%ImpLabel(ng))
            CALL NUOPC_AddNestedState (ImportState,                     &
     &                                 CplSet=TRIM(CoupledSet),         &
     &                                 nestedStateName=TRIM(StateLabel),&
     &                                 nestedState=MODELS(Iroms)%       &
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
            DO i=1,Nimport(Iroms)
              StandardName=MODELS(Iroms)%ImportField(i)%standard_name
              ShortName   =MODELS(Iroms)%ImportField(i)%short_name
              CALL NUOPC_Advertise (MODELS(Iroms)%ImportState(ng),      &
     &                              StandardName=TRIM(StandardName),    &
     &                              name=TRIM(ShortName),               &
     &                              rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF

# ifdef LONGWAVE_OUT
!
              IF (TRIM(ShortName).eq.'LWrad') THEN
                rc=ESMF_RC_NOT_VALID
                IF (localPET.eq.0) THEN
                  WRITE (cplout,10) TRIM(ShortName), 'LONGWAVE_OUT',    &
     &                           'downward longwave radiation: dLWrad', &
     &                           'LONGWAVE_OUT'
                END IF
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
              END IF
# endif
            END DO
          END IF
        END DO
      END IF IMPORTING
!
!-----------------------------------------------------------------------
!  Set ROMS Export State metadata.
!-----------------------------------------------------------------------
!
!  Add ROMS import state. If nesting, each grid has its own import
!  state.
!
      EXPORTING : IF (Nexport(Iroms).gt.0) THEN
        DO ng=1,MODELS(Iroms)%Ngrids
          IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
            CoupledSet=TRIM(COUPLED(Iroms)%SetLabel(ng))
            StateLabel=TRIM(COUPLED(Iroms)%ExpLabel(ng))
            CALL NUOPC_AddNestedState (ExportState,                     &
     &                                 CplSet=TRIM(CoupledSet),         &
     &                                 nestedStateName=TRIM(StateLabel),&
     &                                 nestedState=MODELS(Iroms)%       &
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
            DO i=1,Nexport(Iroms)
              StandardName=MODELS(Iroms)%ExportField(i)%standard_name
              ShortName   =MODELS(Iroms)%ExportField(i)%short_name
              CALL NUOPC_Advertise (MODELS(Iroms)%ExportState(ng),      &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_SetInitializeP1',    &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
# ifdef LONGWAVE_OUT
  10  FORMAT (/,' ROMS_SetInitializeP1 - incorrect field to process: ', &
     &        a,/,24x,'when activating option: ',a,/,24x,               &
     &        'use instead ',a,/,24x,'or deactivate option: ',a,/)
!
# endif
!
      RETURN
      END SUBROUTINE ROMS_SetInitializeP1
!
      SUBROUTINE ROMS_SetInitializeP2 (model,                           &
     &                                 ImportState, ExportState,        &
     &                                 clock, rc)
!
!=======================================================================
!                                                                      !
!  ROMS component Phase 2 initialization: Initializes ROMS, sets       !
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
      logical, save :: first
!
      integer :: LBi, UBi, LBj, UBj
      integer :: MyComm
      integer :: ng, localPET, PETcount, tile
!
      real (dp) :: driverDuration, romsDuration
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_SetInitializeP2"
!
      TYPE (ESMF_TimeInterval) :: RunDuration, TimeStep
      TYPE (ESMF_Time)         :: CurrTime, startTime
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_SetInitializeP2',    &
     &                          ', PET', PETrank
        FLUSH (trac)
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
      tile=localPET
      ESMcomm(Iroms)=MyComm
!
!-----------------------------------------------------------------------
!  Initialize ROMS component.  In nested applications, ROMS kernel will
!  allocate and initialize all grids with a single call to
!  "ROMS_initialize".
!-----------------------------------------------------------------------
!
      first=.TRUE.
      CALL ROMS_initialize (first, mpiCOMM=MyComm)
      IF (exit_flag.ne.NoError) THEN
        rc=ESMF_RC_OBJ_INIT
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF

# ifdef TIME_INTERP
!
!-----------------------------------------------------------------------
!  Create field time interpolation variable attributes NetCDF file. It
!  needs to be done after ROMS initialization since the NetCDF and
!  mpi interface use several variables from ROMS profiling that need
!  to be allocated.
!-----------------------------------------------------------------------
!
      IF (PETlayoutOption.eq.'CONCURRENT') THEN
        CALL def_FieldAtt (vm, rc)
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
!  Check ROMS simulation length and compare with that of the coupling
!  driver.  We need to use the driver clock here since the ROMS
!  component clock has been not created before this intialization
!  phase.
!-----------------------------------------------------------------------
!
      IF (MODELS(Iroms)%IsActive) THEN
        CALL ESMF_ClockGet (ClockInfo(Idriver)%Clock,                   &
     &                      currTime=CurrTime,                          &
     &                      timeStep=TimeStep,                          &
     &                      runDuration=RunDuration,                    &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
# ifdef REGRESS_STARTCLOCK
        CALL ESMF_TimeIntervalGet (RunDuration-TimeStep,                &
     &                             s_r8=driverDuration,                 &
     &                             rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
# else
        CALL ESMF_TimeIntervalGet (RunDuration,                         &
     &                             s_r8=driverDuration,                 &
     &                             rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
# endif
!
        DO ng=1,MODELS(Iroms)%Ngrids
          IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
            romsDuration=(ntend(ng)-ntfirst(ng)+1)*dt(ng)
            IF (romsDuration.ne.driverDuration) THEN
              IF (localPET.eq.0) THEN
                WRITE (cplout,10) romsDuration, driverDuration,         &
     &                          TRIM(INPname(Iroms))
              END IF
              rc=ESMF_RC_NOT_VALID
              RETURN
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Set-up grid and load coordinate data.
!-----------------------------------------------------------------------
!
      DO ng=1,MODELS(Iroms)%Ngrids
        IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
          CALL ROMS_SetGridArrays (ng, tile, model, rc)
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
      DO ng=1,MODELS(Iroms)%Ngrids
        IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
          CALL ROMS_SetStates (ng, tile, model, rc)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_SetInitializeP2',    &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
  10  FORMAT (/,' ROMS_SetInitializeP2 - inconsitent configuration ',   &
     &        'run duration',/,24x,                                     &
     &        'ROMS Duration     = ',f20.2,' seconds',/,24x,            &
     &        'Coupling Duration = ',f20.2,' seconds',/,24x,            &
     &        'Check paramenter NTIMES in ''',a,'''',a)
!
      RETURN
      END SUBROUTINE ROMS_SetInitializeP2
!
      SUBROUTINE ROMS_DataInit (model, rc)
!
!=======================================================================
!                                                                      !
!  Exports ROMS component fields during initialization or restart.     !
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
     &  __FILE__//", ROMS_DataInit"
!
      TYPE (ESMF_Time)  :: CurrentTime
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_DataInit',           &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get gridded component clock current time.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (ClockInfo(Iroms)%Clock,                       &
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
      IF (Nexport(Iroms).gt.0) THEN
        DO ng=1,MODELS(Iroms)%Ngrids
          IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
            CALL ROMS_Export (ng, model, rc)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_DataInit',           &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
      RETURN
      END SUBROUTINE ROMS_DataInit
!
      SUBROUTINE ROMS_SetClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets ROMS component date calendar, start and stop time, and         !
!  coupling interval.  At initilization, the variable "tdays" is       !
!  the initial time meassured in fractional days since the reference   !
!  time.                                                               !
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
      integer :: ref_year,   start_year,   stop_year
      integer :: ref_month,  start_month,  stop_month
      integer :: ref_day,    start_day,    stop_day
      integer :: ref_hour,   start_hour,   stop_hour
      integer :: ref_minute, start_minute, stop_minute
      integer :: ref_second, start_second, stop_second
      integer :: PETcount, localPET
      integer :: TimeFrac
!
      real(dp) :: MyStartTime, MyStopTime
!
      character (len= 22) :: Calendar
      character (len= 22) :: StartTimeString, StopTimeString
      character (len=160) :: message

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_SetClock"
!
      TYPE (ESMF_CalKind_Flag) :: CalType
      TYPE (ESMF_Clock)        :: clock
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_SetClock',           &
     &                          ', PET', PETrank
        FLUSH (trac)
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
!  Create ROMS component clock.
!-----------------------------------------------------------------------
!
!  Set ROMS time reference: model time is meassured as seconds since
!  reference time.  ESMF does not support the Proleptic Gregorian
!  Calendar that extends backward the dates preceeding 15 October 1582
!  which always have a year length of 365.2425 days.
!
      ref_year  =Rclock%year
      ref_month =Rclock%month
      ref_day   =Rclock%day
      ref_hour  =Rclock%hour
      ref_minute=Rclock%minutes
      ref_second=Rclock%seconds
      Calendar  =TRIM(Rclock%calendar)
!
      IF (INT(time_ref).eq.-1) THEN
        CalType=ESMF_CALKIND_360DAY
      ELSE
        CalType=ESMF_CALKIND_GREGORIAN
      END IF
!
      ClockInfo(Iroms)%Calendar=ESMF_CalendarCreate(CalType,            &
     &                                              name=TRIM(Calendar),&
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
      CALL ESMF_TimeSet (ClockInfo(Iroms)%ReferenceTime,                &
     &                   yy=ref_year,                                   &
     &                   mm=ref_month,                                  &
     &                   dd=ref_day,                                    &
     &                   h =ref_hour,                                   &
     &                   m =ref_minute,                                 &
     &                   s =ref_second,                                 &
     &                   calendar=ClockInfo(Iroms)%Calendar,            &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF

# ifdef REGRESS_STARTCLOCK
!
!  Set start time, use the minimum value of all nested grids. Notice
!  that a coupling interval is substracted since the driver clock was
!  regressed by that amount to properly initialize all ESM components.
!
      MyStartTime=MINVAL(tdays)-ClockInfo(Iroms)%Time_Step/86400.0_dp
# else
!
!  Set start time, use the minimum value of all nested grids.
!
      MyStartTime=MINVAL(tdays)
# endif
!
      ClockInfo(Iroms)%Time_Start=MyStartTime*86400.0_dp
      CALL caldate (MyStartTime,                                        &
     &              yy_i=start_year,                                    &
     &              mm_i=start_month,                                   &
     &              dd_i=start_day,                                     &
     &              h_i =start_hour,                                    &
     &              m_i =start_minute,                                  &
     &              s_i =start_second)
      CALL time_string (ClockInfo(Iroms)%Time_Start,                    &
     &                  ClockInfo(Iroms)%Time_StartString)
!
      CALL ESMF_TimeSet (ClockInfo(Iroms)%StartTime,                    &
     &                   yy=start_year,                                 &
     &                   mm=start_month,                                &
     &                   dd=start_day,                                  &
     &                   h =start_hour,                                 &
     &                   m =start_minute,                               &
     &                   s =start_second,                               &
     &                   ms=0,                                          &
     &                   calendar=ClockInfo(Iroms)%Calendar,            &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set stop time, use the maximum value of all nested grids.
!
      MyStopTime=0.0_dp
      DO ng=1,MODELS(Iroms)%Ngrids
        IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
          MyStopTime=MAX(MyStopTime,                                    &
     &                   tdays(ng)+(REAL(ntimes(ng),dp)*dt(ng))*sec2day)
        END IF
      END DO
      ClockInfo(Iroms)%Time_Stop=MyStopTime*86400.0_dp
      CALL caldate (MyStopTime,                                         &
     &              yy_i=stop_year,                                     &
     &              mm_i=stop_month,                                    &
     &              dd_i=stop_day,                                      &
     &              h_i =stop_hour,                                     &
     &              m_i =stop_minute,                                   &
     &              s_i =stop_second)
      CALL time_string (ClockInfo(Iroms)%Time_Stop,                     &
     &                  ClockInfo(Iroms)%Time_StopString)
!
      CALL ESMF_TimeSet (ClockInfo(Iroms)%StopTime,                     &
     &                   yy=stop_year,                                  &
     &                   mm=stop_month,                                 &
     &                   dd=stop_day,                                   &
     &                   h =stop_hour,                                  &
     &                   m =stop_minute,                                &
     &                   s =stop_second,                                &
     &                   calendar=ClockInfo(Iroms)%Calendar,            &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Modify component clock time step.
!-----------------------------------------------------------------------
!
      TimeFrac=0
      DO ng=1,MODELS(Iroms)%Ngrids
        IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
          TimeFrac=MAX(TimeFrac,                                        &
     &                 MAXVAL(MODELS(Iroms)%TimeFrac(ng,:),             &
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
      ClockInfo(Iroms)%TimeStep=ClockInfo(Idriver)%TimeStep/TimeFrac
!
!-----------------------------------------------------------------------
!  Create ROMS component clock.
!-----------------------------------------------------------------------
!
      ClockInfo(Iroms)%Name='ROMS_clock'
      clock=ESMF_ClockCreate(ClockInfo(Iroms)%TimeStep,                 &
     &                       ClockInfo(Iroms)%StartTime,                &
     &                       stopTime =ClockInfo(Iroms)%StopTime,       &
     &                       refTime  =ClockInfo(Iroms)%ReferenceTime,  &
     &                       name     =TRIM(ClockInfo(Iroms)%Name),     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      ClockInfo(Iroms)%Clock=clock
!
!  Set ROMS component clock.
!
      CALL ESMF_GridCompSet (model,                                     &
     &                       clock=ClockInfo(Iroms)%Clock,              &
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
      CALL ESMF_ClockGet (ClockInfo(Iroms)%Clock,                       &
     &                    currTime=ClockInfo(Iroms)%CurrentTime,        &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Compare driver time against ROMS component time.
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
        WRITE (cplout,'(/)')
        WRITE (cplout,10) 'DRIVER Calendar:    ',                       &
     &                    TRIM(ClockInfo(Idriver)%CalendarString),      &
     &                    'DRIVER Start Clock: ',                       &
     &                    TRIM(ClockInfo(Idriver)%Time_StartString),    &
     &                    'DRIVER Stop  Clock: ',                       &
     &                    TRIM(ClockInfo(Idriver)%Time_StopString)
!
        WRITE (cplout,10) 'ROMS Calendar:      ',                       &
     &                    TRIM(ClockInfo(Iroms)%CalendarString),        &
     &                    'ROMS Start Clock:   ',                       &
     &                    TRIM(ClockInfo(Iroms)%Time_StartString),      &
     &                    'ROMS Stop  Clock:   ',                       &
     &                    TRIM(ClockInfo(Iroms)%Time_StopString)
      END IF
!
!  Compare Driver and ROMS clocks.
!
      IF (ClockInfo(Iroms)%Time_StartString(1:19).ne.                   &
     &    StartTimeString(1:19)) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) 'ROMS   Start Time: ',                      &
     &                      ClockInfo(Iroms)%Time_StartString(1:19),    &
     &                      'Driver Start Time: ',                      &
     &                      TRIM(StartTimeString),                      &
     &                      '                   are not equal!'
        END IF
        message='Driver and ROMS start times do not match: '//          &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ClockInfo(Iroms  )%Time_StopString(1:19).ne.                  &
     &    ClockInfo(Idriver)%Time_StopString(1:19)) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) 'ROMS   Stop Time: ',                       &
     &                      ClockInfo(Iroms  )%Time_StopString(1:19),   &
     &                      'Driver Stop Time: ',                       &
     &                      TRIM(ClockInfo(Idriver)%Time_StopString),   &
     &                      '                   are not equal!'
        END IF
        message='Driver and ROMS stop times do not match: '//           &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (TRIM(ClockInfo(Iroms  )%CalendarString).ne.                   &
     &    TRIM(ClockInfo(Idriver)%CalendarString)) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) 'ROMS   Calendar: ',                        &
     &                      TRIM(ClockInfo(Iroms  )%CalendarString),    &
     &                      'Driver Calendar: ',                        &
     &                      TRIM(ClockInfo(Idriver)%CalendarString),    &
     &                      '                  are not equal!'
        END IF
        message='Driver and ROMS calendars do not match: '//            &
     &          'please check the config files.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_SetClock',           &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
 10   FORMAT (2x,a,2x,a/,2x,a,2x,a,/,2x,a,2x,a,/)
 20   FORMAT (/,2x,a,a,/,2x,a,a,/,2x,a)
!
      RETURN
      END SUBROUTINE ROMS_SetClock

# ifdef ESM_SETRUNCLOCK
!
      SUBROUTINE ROMS_SetRunClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets ROMS run clock manually to avoid getting zero time stamps at   !
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
     &  __FILE__//", ROMS_SetRunClock"
!
      TYPE (ESMF_Clock) :: driverClock, modelClock
      TYPE (ESMF_Time)  :: currTime
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_SetRunClock',        &
     &                          ', PET', PETrank
        FLUSH (trac)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_SetRunClock',        &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
      RETURN
      END SUBROUTINE ROMS_SetRunClock
# endif
!
      SUBROUTINE ROMS_CheckImport (model, rc)
!
!=======================================================================
!                                                                      !
!  Checks if ROMS component import field is at the correct time.       !
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
     &  __FILE__//", ROMS_CheckImport"
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
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_CheckImport',        &
     &                          ', PET', PETrank
        FLUSH (trac)
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
!  Adjust driver clock for semi-implicit coupling.

      IF (CouplingType.eq.1) THEN
        CurrentTime=DriverTime                  ! explicit coupling
      ELSE
        CurrentTime=DRiverTime+TimeStep         ! semi-implicit coupling
      END IF
!
      CALL ESMF_TimeGet (CurrentTime,                                   &
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
      IF (Nimport(Iroms).gt.0) THEN
        NESTED_LOOP : DO ng=1,MODELS(Iroms)%Ngrids
          IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
            CALL ESMF_StateGet (MODELS(Iroms)%ImportState(ng),          &
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
            CALL ESMF_StateGet (MODELS(Iroms)%ImportState(ng),          &
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
              CALL ESMF_StateGet (MODELS(Iroms)%ImportState(ng),        &
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
              string='ROMS_CheckImport - '//TRIM(FieldName)//' field'
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_CheckImport',        &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
  10  FORMAT (1x,'ROMS_CheckImport - ',a,':',t32,'TimeStamp = ',a,      &
     &        ',  DriverTime = ',a)
!
      RETURN
      END SUBROUTINE ROMS_CheckImport
!
      SUBROUTINE ROMS_SetGridArrays (ng, tile, model, rc)
!
!=======================================================================
!                                                                      !
!  Sets ROMS component staggered, horizontal grids arrays, grid area,  !
!  and land/sea mask, if any.                                          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng, tile
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp), intent(inout) :: model
!
!  Local variable declarations.
!
      integer :: MyTile, gtype, i, ivar, j, node
      integer :: Istr,  Iend,  Jstr,  Jend
      integer :: IstrR, IendR, JstrR, JendR
      integer :: localDE, localDEcount
      integer :: staggerEdgeLWidth(2)
      integer :: staggerEdgeUWidth(2)
!
      integer, allocatable :: deBlockList(:,:,:)
      integer (i4b), pointer :: ptrM(:,:) => NULL()     ! land/sea mask
!
      real (dp), pointer :: ptrA(:,:) => NULL()         ! area
      real (dp), pointer :: ptrX(:,:) => NULL()         ! longitude
      real (dp), pointer :: ptrY(:,:) => NULL()         ! latitude
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_SetGridArrays"
!
      TYPE (ESMF_DistGrid)   :: distGrid
      TYPE (ESMF_StaggerLoc) :: staggerLoc
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_SetGridArrays',      &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Set limits of the grid arrays based on tile decomposition (MPI rank)
!  and nested grid number.
!-----------------------------------------------------------------------
!
      IstrR=BOUNDS(ng)%IstrR(tile)   ! Full range I-starting (RHO)
      IendR=BOUNDS(ng)%IendR(tile)   ! Full range I-ending   (RHO)
      JstrR=BOUNDS(ng)%JstrR(tile)   ! Full range J-starting (RHO)
      JendR=BOUNDS(ng)%JendR(tile)   ! Full range J-ending   (RHO)
!
      Istr=BOUNDS(ng)%Istr(tile)     ! Full range I-starting (PSI, U)
      Iend=BOUNDS(ng)%Iend(tile)     ! Full range I-ending   (PSI)
      Jstr=BOUNDS(ng)%Jstr(tile)     ! Full range J-starting (PSI, V)
      Jend=BOUNDS(ng)%Jend(tile)     ! Full range J-ending   (PSI)
!
!  Set tiles lower and upper bounds for each decomposition element.
!  In ROMS, the "exclusive region" for each decomposition element or
!  horizontal tile ranges is bounded by (Istr:Iend, Jstr:Jend). Each
!  tiled array is dimensioned as (LBi:UBi, LBj:UBj) which includes
!  halo regions (usually 2 ghost points) and padding when appropriate
!  (total/memory region). All ROMS arrays are horizontally dimensioned
!  with the same bounds regardless if they are variables located at
!  RHO-, PSI-, U-, or V-points. There is no halos at the boundary edges.
!  The physical boundary is a U-points (east/west edges) and V-points
!  (south/north edges). The boundary for RHO-points variables are
!  located at half grid (dx,dy) distance away from the physical boundary
!  at array indices(i=0; i=Lm+1) and (j=0; j=Mm+1).
!
!               --------------------- UBj      ESMF uses a very
!              |                     |         complicated array
!              | Jend __________     |         regions:
!              |     |          |    |
!              |     |          |    |         * interior region
!              |     |          |    |         * exclusive region
!              | Jstr|__________|    |         * computational region
!              |     Istr    Iend    |         * total (memory) region
!              |                     |
!               --------------------- LBj
!               LBi               UBi
!
      IF (.not.allocated(deBlockList)) THEN
        allocate ( deBlockList(2,2,NtileI(ng)*NtileJ(ng)) )
      END IF
      DO MyTile=0,NtileI(ng)*NtileJ(ng)-1
        deBlockList(1,1,MyTile+1)=BOUNDS(ng)%Istr(MyTile)
        deBlockList(1,2,MyTile+1)=BOUNDS(ng)%Iend(MyTile)
        deBlockList(2,1,MyTile+1)=BOUNDS(ng)%Jstr(MyTile)
        deBlockList(2,2,MyTile+1)=BOUNDS(ng)%Jend(MyTile)
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
     &                             maxIndex=(/ Lm(ng), Mm(ng) /),       &
     &                             deBlockList=deBlockList,             &
     &                             rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Report ROMS DistGrid based on model domain decomposition.
!
      IF ((tile.eq.0).and.(DebugLevel.gt.0)) THEN
        WRITE (cplout,10) ng, TRIM(GridType(Icenter))//" Point",        &
     &                    NtileI(ng), NtileJ(ng)
        DO MyTile=1,NtileI(ng)*NtileJ(ng)
          WRITE (cplout,20) MyTile-1, deBlockList(1,1,MyTile),          &
     &                                deBlockList(1,2,MyTile),          &
     &                                deBlockList(2,1,MyTile),          &
     &                                deBlockList(2,2,MyTile)
        END DO
      END IF
      IF (allocated(deBlockList)) deallocate (deBlockList)
!
!-----------------------------------------------------------------------
!  Set component grid coordinates.
!-----------------------------------------------------------------------
!
!  Define component grid location type: Arakawa C-grid.
!
!    Icenter:  RHO-point, cell center
!    Icorner:  PSI-point, cell corners
!    Iupoint:  U-point,   cell west/east sides
!    Ivpoint:  V-point,   cell south/north sides
!
      IF (.not.allocated(MODELS(Iroms)%mesh)) THEN
        allocate ( MODELS(Iroms)%mesh(4) )
        MODELS(Iroms)%mesh(1)%gtype=Icenter
        MODELS(Iroms)%mesh(2)%gtype=Icorner
        MODELS(Iroms)%mesh(3)%gtype=Iupoint
        MODELS(Iroms)%mesh(4)%gtype=Ivpoint
      END IF
!
!  Create ESMF Grid. The array indices are global following ROMS
!  design.
!
      MODELS(Iroms)%grid(ng)=ESMF_GridCreate(distgrid=distGrid,         &
     &                                   gridEdgeLWidth=(/1,1/),        &
     &                                   gridEdgeUWidth=(/1,1/),        &
     &                                   indexflag=ESMF_INDEX_GLOBAL,   &
     &                                   name=TRIM(MODELS(Iroms)%name), &
     &                                   rc=rc)
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
      CALL ESMF_GridGet (MODELS(Iroms)%grid(ng),                        &
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
      MESH_LOOP : DO ivar=1,UBOUND(MODELS(Iroms)%mesh, DIM=1)
!
!  Set staggering type, Arakawa C-grid.
!
        SELECT CASE (MODELS(Iroms)%mesh(ivar)%gtype)
          CASE (Icenter)
            staggerLoc=ESMF_STAGGERLOC_CENTER
            staggerEdgeLWidth=(/1,1/)
            staggerEdgeUWidth=(/1,1/)
          CASE (Icorner)
            staggerLoc=ESMF_STAGGERLOC_CORNER
            staggerEdgeLWidth=(/0,0/)
            staggerEdgeUWidth=(/1,1/)
          CASE (Iupoint)
            staggerLoc=ESMF_STAGGERLOC_EDGE1
            staggerEdgeLWidth=(/0,1/)
            staggerEdgeUWidth=(/1,1/)
          CASE (Ivpoint)
            staggerLoc=ESMF_STAGGERLOC_EDGE2
            staggerEdgeLWidth=(/1,0/)
            staggerEdgeUWidth=(/1,1/)
        END SELECT
!
!  Allocate coordinate storage associated with staggered grid type.
!  No coordinate values are set yet.
!
        CALL ESMF_GridAddCoord (MODELS(Iroms)%grid(ng),                 &
     &                          staggerLoc=staggerLoc,                  &
     &                          staggerEdgeLWidth=staggerEdgeLWidth,    &
     &                          staggerEdgeUWidth=staggerEdgeUWidth,    &
     &                          rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF

# ifdef MASKING
!
!  Allocate storage for land/sea masking.
!
        CALL ESMF_GridAddItem (MODELS(Iroms)%grid(ng),                  &
     &                         staggerLoc=staggerLoc,                   &
     &                         itemflag=ESMF_GRIDITEM_MASK,             &
     &                         rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
        MODELS(Iroms)%LandValue=0
        MODELS(Iroms)%SeaValue=1
# endif
!
!  Allocate storage for grid area.
!
        CALL ESMF_GridAddItem (MODELS(Iroms)%grid(ng),                  &
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
          CALL ESMF_GridGetCoord (MODELS(Iroms)%grid(ng),               &
     &                            coordDim=1,                           &
     &                            localDE=localDE,                      &
     &                            staggerLoc=staggerLoc,                &
     &                            farrayPtr=ptrX,                       &
     &                            rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_GridGetCoord (MODELS(Iroms)%grid(ng),               &
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
          CALL ESMF_GridGetItem (MODELS(Iroms)%grid(ng),                &
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
          CALL ESMF_GridGetItem (MODELS(Iroms)%grid(ng),                &
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
          SELECT CASE (MODELS(Iroms)%mesh(ivar)%gtype)
            CASE (Icenter)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  ptrX(i,j)=GRID(ng)%lonr(i,j)
                  ptrY(i,j)=GRID(ng)%latr(i,j)
# ifdef MASKING
                  ptrM(i,j)=INT(GRID(ng)%rmask(i,j))
# else
                  ptrM(i,j)=1
# endif
                  ptrA(i,j)=GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j)
                END DO
              END DO
            CASE (Icorner)
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  ptrX(i,j)=GRID(ng)%lonp(i,j)
                  ptrY(i,j)=GRID(ng)%latp(i,j)
# ifdef MASKING
                  ptrM(i,j)=INT(GRID(ng)%pmask(i,j))
# else
                  ptrM(i,j)=1
# endif
                  ptrA(i,j)=GRID(ng)%om_p(i,j)*GRID(ng)%on_p(i,j)
                END DO
              END DO
            CASE (Iupoint)
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  ptrX(i,j)=GRID(ng)%lonu(i,j)
                  ptrY(i,j)=GRID(ng)%latu(i,j)
# ifdef MASKING
                  ptrM(i,j)=INT(GRID(ng)%umask(i,j))
# else
                  ptrM(i,j)=1
# endif
                  ptrA(i,j)=GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)
                END DO
              END DO
            CASE (Ivpoint)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  ptrX(i,j)=GRID(ng)%lonv(i,j)
                  ptrY(i,j)=GRID(ng)%latv(i,j)
# ifdef MASKING
                  ptrM(i,j)=INT(GRID(ng)%vmask(i,j))
# else
                  ptrM(i,j)=1
# endif
                  ptrA(i,j)=GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)
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
          gtype=MODELS(Iroms)%mesh(ivar)%gtype
          CALL ESMF_GridWriteVTK (MODELS(Iroms)%grid(ng),               &
     &                            filename="roms_"//                    &
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
     &                       grid=MODELS(Iroms)%grid(ng),               &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_SetGridArrays',      &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      IF (DebugLevel.gt.0) FLUSH (cplout)
!
  10  FORMAT (2x,'ROMS_DistGrid - Grid = ',i2.2,',',3x,'Mesh = ',a,',', &
     &        3x,'Partition = ',i0,' x ',i0)
  20  FORMAT (18x,'node = ',i0,t32,'Istr = ',i0,t45,'Iend = ',i0,       &
     &                         t58,'Jstr = ',i0,t71,'Jend = ',i0)
!
      RETURN
      END SUBROUTINE ROMS_SetGridArrays
!
      SUBROUTINE ROMS_SetStates (ng, tile, model, rc)
!
!=======================================================================
!                                                                      !
!  Adds ROMS component export and import fields into its respective    !
!  state.                                                              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng, tile
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: id, ifld
      integer :: localDE, localDEcount, localPET
      integer :: ExportCount, ImportCount
      integer :: staggerEdgeLWidth(2)
      integer :: staggerEdgeUWidth(2)
      integer :: haloLW(2), haloUW(2)
!
      real (dp), dimension(:,:), pointer :: ptr2d => NULL()
!
      character (len=10) :: AttList(1)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_SetStates"
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
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_SetStates',          &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds tile halo widths for ESMF fields.
!-----------------------------------------------------------------------
!
      haloLW(1)=BOUNDS(ng)%Istr(tile)-BOUNDS(ng)%LBi (tile)
      haloLW(2)=BOUNDS(ng)%Jstr(tile)-BOUNDS(ng)%LBj (tile)
      haloUW(1)=BOUNDS(ng)%UBi (tile)-BOUNDS(ng)%Iend(tile)
      haloUW(2)=BOUNDS(ng)%UBj (tile)-BOUNDS(ng)%Jend(tile)
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
      CALL ESMF_GridGet (MODELS(Iroms)%grid(ng),                        &
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
      EXPORTING : IF (Nexport(Iroms).gt.0) THEN
!
!  Get number of fields to export.
!
        CALL ESMF_StateGet (MODELS(Iroms)%ExportState(ng),              &
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
        CALL ESMF_StateGet (MODELS(Iroms)%ExportState(ng),              &
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
        DO ifld=1,ExportCount
          id=field_index(MODELS(Iroms)%ExportField,ExportNameList(ifld))
!
          IF (NUOPC_IsConnected(MODELS(Iroms)%ExportState(ng),          &
     &                          fieldName=TRIM(ExportNameList(ifld)),   &
     &                          rc=rc)) THEN
!
!  Set staggering type.
!
            SELECT CASE (MODELS(Iroms)%ExportField(id)%gtype)
              CASE (Icenter)                                ! RHO-points
                staggerLoc=ESMF_STAGGERLOC_CENTER
                staggerEdgeLWidth=(/1,1/)
                staggerEdgeUWidth=(/1,1/)
              CASE (Icorner)                                ! PSI-points
                staggerLoc=ESMF_STAGGERLOC_CORNER
                staggerEdgeLWidth=(/0,0/)
                staggerEdgeUWidth=(/1,1/)
              CASE (Iupoint)                                ! U-points
                staggerLoc=ESMF_STAGGERLOC_EDGE1
                staggerEdgeLWidth=(/0,1/)
                staggerEdgeUWidth=(/1,1/)
              CASE (Ivpoint)                                ! V-points
                staggerLoc=ESMF_STAGGERLOC_EDGE2
                staggerEdgeLWidth=(/1,0/)
                staggerEdgeUWidth=(/1,1/)
            END SELECT
!
!  Create 2D field from the Grid and arraySpec.
!
            field=ESMF_FieldCreate(MODELS(Iroms)%grid(ng),              &
     &                             arraySpec2d,                         &
     &                             indexflag=ESMF_INDEX_GLOBAL,         &
     &                             staggerloc=staggerLoc,               &
     &                             totalLWidth=haloLW,                  &
     &                             totalUWidth=haloUW,                  &
     &                             name=TRIM(ExportNameList(ifld)),     &
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
            CALL NUOPC_Realize (MODELS(Iroms)%ExportState(ng),          &
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
              WRITE (cplout,10) TRIM(ExportNameList(ifld)),             &
     &                          'Export State: ',                       &
     &                          TRIM(COUPLED(Iroms)%ExpLabel(ng))
            END IF
            CALL ESMF_StateRemove (MODELS(Iroms)%ExportState(ng),       &
     &                             (/ TRIM(ExportNameList(ifld)) /),    &
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
      IMPORTING : IF (Nimport(Iroms).gt.0) THEN
!
!  Get number of fields to import.
!
        CALL ESMF_StateGet (MODELS(Iroms)%ImportState(ng),              &
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
        CALL ESMF_StateGet (MODELS(Iroms)%ImportState(ng),              &
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
        DO ifld=1,ImportCount
          id=field_index(MODELS(Iroms)%ImportField,ImportNameList(ifld))
!
          IF (NUOPC_IsConnected(MODELS(Iroms)%ImportState(ng),          &
     &                          fieldName=TRIM(ImportNameList(ifld)),   &
     &                          rc=rc)) THEN
!
!  Set staggering type.
!
            SELECT CASE (MODELS(Iroms)%ImportField(id)%gtype)
              CASE (Icenter)                                ! RHO-points
                staggerLoc=ESMF_STAGGERLOC_CENTER
                staggerEdgeLWidth=(/1,1/)
                staggerEdgeUWidth=(/1,1/)
              CASE (Icorner)                                ! PSI-points
                staggerLoc=ESMF_STAGGERLOC_CORNER
                staggerEdgeLWidth=(/0,0/)
                staggerEdgeUWidth=(/1,1/)
              CASE (Iupoint)                                ! U-points
                staggerLoc=ESMF_STAGGERLOC_EDGE1
                staggerEdgeLWidth=(/0,1/)
                staggerEdgeUWidth=(/1,1/)
              CASE (Ivpoint)                                ! V-points
                staggerLoc=ESMF_STAGGERLOC_EDGE2
                staggerEdgeLWidth=(/1,0/)
                staggerEdgeUWidth=(/1,1/)
            END SELECT
!
!  Create 2D field from the Grid, arraySpec, total tile size including
!  halos.  The array indices are global following ROMS design.
!
            field=ESMF_FieldCreate(MODELS(Iroms)%grid(ng),              &
     &                             arraySpec2d,                         &
     &                             indexflag=ESMF_INDEX_GLOBAL,         &
     &                             staggerloc=staggerLoc,               &
     &                             totalLWidth=haloLW,                  &
     &                             totalUWidth=haloUW,                  &
     &                             name=TRIM(ImportNameList(ifld)),     &
     &                             rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF

# ifdef TIME_INTERP_NOT
!
!  Create standard Attribute Package for each export field. Then, nest
!  custom Attribute Package around it.
!
            CALL ESMF_AttributeAdd (field,                              &
     &                              convention='ESMF',                  &
     &                              purpose='General',                  &
     &                              rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
            AttList(1)='TimeInterp'
            CALL ESMF_AttributeAdd (field,                              &
     &                              convention='CustomConvention',      &
     &                              purpose='General',                  &
!!   &                              purpose='Instance',                 &
     &                              attrList=AttList,                   &
     &                              nestConvention='ESMF',              &
     &                              nestPurpose='General',              &
     &                              rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
# endif
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
            CALL NUOPC_Realize (MODELS(Iroms)%ImportState(ng),          &
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
              WRITE (cplout,10) TRIM(ImportNameList(ifld)),             &
     &                          'Import State: ',                       &
     &                          TRIM(COUPLED(Iroms)%ImpLabel(ng))
            END IF
            CALL ESMF_StateRemove (MODELS(Iroms)%ImportState(ng),       &
     &                             (/ TRIM(ImportNameList(ifld)) /),    &
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_SetStates',          &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
 10   FORMAT (1x,'ROMS_SetStates - Removing field ''',a,''' from ',a,   &
     &        '''',a,'''',/,18x,'because it is not connected.')
!
      RETURN
      END SUBROUTINE ROMS_SetStates
!
      SUBROUTINE ROMS_ModelAdvance (model, rc)
!
!=======================================================================
!                                                                      !
!  Advance ROMS component for a coupling interval (seconds) using      !
!  "ROMS_run". It also calls "ROMS_Import" and "ROMS_Export" to        !
!  import and export coupling fields, respectively.                    !
!                                                                      !
!  During configuration, the driver clock was decreased by a single    !
!  coupling interval (TimeStep) to allow the proper initialization     !
!  of the import and export fields pointers.  ROMS is not advanced     !
!  on the first call to this routine, so the time stepping is over     !
!  the specified application start and ending dates.                   !
!                                                                      !
# if defined TIME_INTERP
!  On the first pass, it imports the LOWER time snapshot fields,       !
!  but cannot time-step ROMS until the next call after importing       !
!  the UPPER snapshot.  Therefore, it starts time-stepping when        !
!  both LOWER and UPPER time snapshot fields are exchanged so that     !
!  ROMS can perform time interpolation.                                !
# else
!  ROMS is actually advanced on the second call to this routine.       !
# endif
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
      logical :: Ladvance
      integer :: is, ng
      integer :: MyTask, PETcount, localPET, phase
!
      real (dp) :: CouplingInterval, RunInterval
      real (dp) :: TcurrentInSeconds, TstopInSeconds
!
      character (len=22) :: Cinterval
      character (len=22) :: CurrTimeString, StopTimeString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_SetModelAdvance"
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
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_ModelAdvance',       &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get information about the gridded component.
!-----------------------------------------------------------------------
!
!  Inquire about ROMS component.
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
     &                    currTime=ClockInfo(Iroms)%CurrentTime,        &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Current ROMS time (seconds).
!
      CALL ESMF_TimeGet (ClockInfo(Iroms)%CurrentTime,                  &
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
!  ROMS stop time (seconds) for this coupling window.
!
      CALL ESMF_TimeGet (ClockInfo(Iroms)%CurrentTime+TimeStep,         &
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
!  Set ROMS running interval (seconds) for the current coupling window.
!
      RunInterval=CouplingInterval
!
!  Set local model advance time stepping switch.
!
      Ladvance=.TRUE.
# ifdef TIME_INTERP
      IF ((MODELS(Iroms)%ImportCalls.eq.0).and.                         &
     &    (Nimport(Iroms).gt.0)) THEN
        Ladvance=.FALSE.
      END IF
# else
#  ifdef REGRESS_STARTCLOCK
      IF (TcurrentInSeconds.eq.ClockInfo(Idriver)%Time_Start) THEN
        Ladvance=.FALSE.
      END IF
#  endif
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
      IF (Nimport(Iroms).gt.0) THEN
        DO ng=1,MODELS(Iroms)%Ngrids
          IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
            CALL ROMS_Import (ng, model, rc)
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
!  Run ROMS component. Notice that ROMS component is advanced when
!  ng=1.  In nested application, ROMS kernel (main2d or main3d) will
!  advance all the nested grid in their logical order.  In nesting,
!  the execution order of the grids is critical since nesting is
!  two-way by default.
!-----------------------------------------------------------------------
!
      IF (Ladvance) THEN
        IF (ESM_track) THEN
          WRITE (trac,'(a,a,i0)') '==> Entering ROMS_Run',              &
     &                            ', PET', PETrank
          FLUSH (trac)
        END IF
        CALL ROMS_run (RunInterval)
        IF (ESM_track) THEN
          WRITE (trac,'(a,a,i0)') '==> Exiting  ROMS_Run',              &
     &                            ', PET', PETrank
          FLUSH (trac)
        END IF
      END IF
!
      IF (exit_flag.ne.NoError) then
        IF (localPET.eq.0) then
          WRITE (cplout,'(a,i1)') 'ROMS component exit with flag = ',   &
     &                            exit_flag
        END IF
        CALL ROMS_finalize
        CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
      END IF
!
!-----------------------------------------------------------------------
!  Put export fields.
!-----------------------------------------------------------------------
!
      IF (Nexport(Iroms).gt.0) THEN
        DO ng=1,MODELS(Iroms)%Ngrids
          IF (ANY(COUPLED(Iroms)%LinkedGrid(ng,:))) THEN
            CALL ROMS_Export (ng, model, rc)
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
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_ModelAdvance',       &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
  10  FORMAT (3x,'ModelAdvance - ESMF, Running ROMS:',t42,a,            &
     &        ' => ',a,', [',a,' s], Advance: ',l1)
!
      RETURN
      END SUBROUTINE ROMS_ModelAdvance
!
      SUBROUTINE ROMS_SetFinalize (model,                               &
     &                             ImportState, ExportState,            &
     &                             clock, rc)
!
!=======================================================================
!                                                                      !
!  Finalize ROMS component execution. It calls ROMS_finalize.          !
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
     &  __FILE__//", ROMS_SetFinalize"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_SetFinalize',        &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  If ng=1, finalize ROMS component. In nesting applications this step
!  needs to be done only once.
!-----------------------------------------------------------------------
!
      CALL ROMS_finalize
      FLUSH (stdout)                      ! flush standard output buffer
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_SetFinalize',        &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
      RETURN
      END SUBROUTINE ROMS_SetFinalize
!
      SUBROUTINE ROMS_Import (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Imports fields into ROMS array structures. The fields aew loaded    !
!  into the snapshot storage arrays to allow time interpolation in     !
!  ROMS kernel.                                                        !
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
      logical :: LoadIt, isPresent
      logical :: got_stress(2), got_wind(2)
# if defined WIND_MINUS_CURRENT && !defined BULK_FLUXES
      logical :: got_RhoAir, got_Wstar, got_wind_sbl(2)
# endif
!
      integer :: Istr,  Iend,  Jstr,  Jend
      integer :: IstrR, IendR, JstrR, JendR
      integer :: LBi, UBi, LBj, UBj
      integer :: ImportCount, Tindex
      integer :: localDE, localDEcount, localPET, tile
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: gtype, id, ifield, ifld, i, is, j
!
# ifdef TIME_INTERP
      integer, save :: record = 0
!
# endif
      real (dp), parameter :: eps = 1.0E-10_dp
!
      real (dp) :: TimeInDays, Time_Current, Tmin, Tmax, Tstr, Tend
# ifdef TIME_INTERP
      real (dp) :: MyTimeInDays
# endif
      real (dp) :: Fseconds, ROMSclockTime
      real (dp) :: MyTintrp(2), MyVtime(2)

      real (dp) :: MyFmax(2), MyFmin(2), Fmin(2), Fmax(2), Fval
      real (dp) :: add_offset, romsScale, scale, cff1, cff2, cff3
      real (dp) :: FreshWaterScale, StressScale, TracerFluxScale
# if defined WIND_MINUS_CURRENT && !defined BULK_FLUXES
      real (dp) :: Urel, Vrel, Wmag, Wrel
# endif
      real (dp) :: AttValues(14)
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
# if defined WIND_MINUS_CURRENT && !defined BULK_FLUXES
      real (dp), allocatable  :: RhoAir(:,:),  Wstar(:,:)
      real (dp), allocatable  :: Uwrk(:,:),    Vwrk(:,:)
      real (dp), allocatable  :: Xwind(:,:),   Ywind(:,:)
# endif
      real (dp), allocatable  :: Ustress(:,:), Vstress(:,:)
      real (dp), allocatable  :: Uwind(:,:),   Vwind(:,:)
!
      character (len=22)      :: MyDate(2)
# ifdef TIME_INTERP
      character (len=22)      :: MyDateString(1,1,1)
# endif
      character (len=22)      :: Time_CurrentString
      character (len=40)      :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_Import"

      character (ESMF_MAXSTR) :: cname, ofile
      character (ESMF_MAXSTR), allocatable :: ImportNameList(:)
!
      TYPE (ESMF_AttPack) :: AttPack
      TYPE (ESMF_Clock)   :: clock
      TYPE (ESMF_Field)   :: field
      TYPE (ESMF_Time)    :: CurrentTime
      TYPE (ESMF_VM)      :: vm

# ifdef TIME_INTERP
!
      SourceFile=MyFile
# endif
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_Import',             &
     &                          ', PET', PETrank
        FLUSH (trac)
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
      CALL ESMF_GridGet (MODELS(Iroms)%grid(ng),                        &
     &                   localDECount=localDEcount,                     &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set size of imported tiled-arrays.
!
      tile=localPET
!
      LBi=BOUNDS(ng)%LBi(tile)        ! lower bound I-direction
      UBi=BOUNDS(ng)%UBi(tile)        ! upper bound I-direction
      LBj=BOUNDS(ng)%LBj(tile)        ! lower bound J-direction
      UBj=BOUNDS(ng)%UBj(tile)        ! upper bound J-direction
!
      IstrR=BOUNDS(ng)%IstrR(tile)    ! Full range I-starting (RHO)
      IendR=BOUNDS(ng)%IendR(tile)    ! Full range I-ending   (RHO)
      JstrR=BOUNDS(ng)%JstrR(tile)    ! Full range J-starting (RHO)
      JendR=BOUNDS(ng)%JendR(tile)    ! Full range J-ending   (RHO)
!
      Istr=BOUNDS(ng)%Istr(tile)      ! Full range I-starting (PSI, U)
      Iend=BOUNDS(ng)%Iend(tile)      ! Full range I-ending   (PSI)
      Jstr=BOUNDS(ng)%Jstr(tile)      ! Full range J-starting (PSI, V)
      Jend=BOUNDS(ng)%Jend(tile)      ! Full range J-ending   (PSI)
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
!  Convert CurrentTime into ROMS clock ellapsed time since
!  initialization in seconds from reference time.
!  (The routine "ROMS_clock" is located in ROMS/Utility/dateclock.F)
!-----------------------------------------------------------------------
!
      Fseconds=REAL(seconds,dp)+REAL(sN,dp)/REAL(sD,dp)
      CALL ROMS_clock (year, month, day, hour, minutes, Fseconds,       &
     &                 ROMSclockTime)
!
!-----------------------------------------------------------------------
!  Get list of import fields.
!-----------------------------------------------------------------------
!
      CALL ESMF_StateGet (MODELS(Iroms)%ImportState(ng),                &
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
      CALL ESMF_StateGet (MODELS(Iroms)%ImportState(ng),                &
     &                    itemNameList=ImportNameList,                  &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF

# ifdef TIME_INTERP
!
!-----------------------------------------------------------------------
!  Advance unlimited dimension counter.
!-----------------------------------------------------------------------
!
      IF (PETlayoutOption.eq.'CONCURRENT') THEN
        record=record+1
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Get import fields.
!-----------------------------------------------------------------------
!
!  Set switches to rotate wind stress and wind component for curvilinear
!  ROMS grid applications.
!
      got_stress(1:2)=.FALSE.
      got_wind(1:2)=.FALSE.
# if defined WIND_MINUS_CURRENT && !defined BULK_FLUXES
      got_RhoAir=.FALSE.
      got_Wstar=.FALSE.
      got_wind_sbl(1:2)=.FALSE.
# endif
!
!  Loop over all import fields to process.
!
      FLD_LOOP : DO ifld=1,ImportCount
        id=field_index(MODELS(Iroms)%ImportField, ImportNameList(ifld))
!
!  Get field from import state.
!
        CALL ESMF_StateGet (MODELS(Iroms)%ImportState(ng),              &
     &                      TRIM(ImportNameList(ifld)),                 &
     &                      field,                                      &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF

# ifdef TIME_INTERP
!
!  If cuncurrent coupling and importing time snapshots, update values
!  in the MODELS(Iroms)%ImportField structure by reading import field
!  interpolation attributes from source NetCDF file. It is very tricky
!  to perform inter VM communications. It is easier to read them from
!  a NetCDF file.  ROMS needs these attributes to perform the time
!  interpolation between snapshots in its kernel.
!  (HGA: need to figure out how to do inter VM communications)
!
        IF (PETlayoutOption.eq.'CONCURRENT') THEN
          CALL netcdf_get_ivar (ng, iNLM, AttFileName, 'Tindex',        &
     &                          MODELS(Iroms)%ImportField(id)%Tindex,   &
     &                          start=(/Iroms,id,record/),              &
     &                          total=(/1,1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_READ
            RETURN
          END IF
!
          is=MODELS(Iroms)%ImportField(id)%Tindex
          CALL netcdf_get_svar (ng, iNLM, AttFileName, 'Date',          &
     &                          MyDateString,                           &
     &                          start=(/1,Iroms,id,record/),            &
     &                          total=(/22,1,1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_READ
            RETURN
          END IF
          MODELS(Iroms)%ImportField(id)%DateString(is)=                 &
     &                                               MyDateString(1,1,1)
!
          CALL netcdf_get_time (ng, iNLM, AttFileName, 'Tcurrent',      &
     &                          Rclock%DateNumber, MyTimeInDays,        &
     &                          start=(/Iroms,id,record/),              &
     &                          total=(/1,1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_READ
            RETURN
          END IF
!
          CALL netcdf_get_time (ng, iNLM, AttFileName, 'Tstr',          &
     &                          Rclock%DateNumber,                      &
     &                          MODELS(Iroms)%ImportField(id)%Tstr,     &
     &                          start=(/Iroms,id,record/),              &
     &                          total=(/1,1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_READ
            RETURN
          END IF
!
          CALL netcdf_get_time (ng, iNLM, AttFileName, 'Tend',          &
     &                          Rclock%DateNumber,                      &
     &                          MODELS(Iroms)%ImportField(id)%Tend,     &
     &                          start=(/Iroms,id,record/),              &
     &                          total=(/1,1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_READ
            RETURN
          END IF
!
          CALL netcdf_get_time (ng, iNLM, AttFileName, 'Tintrp',        &
     &                          Rclock%DateNumber,                      &
     &                        MODELS(Iroms)%ImportField(id)%Tintrp(is), &
     &                          start=(/Iroms,id,record/),              &
     &                          total=(/1,1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_READ
            RETURN
          END IF
!
          CALL netcdf_get_time (ng, iNLM, AttFileName, 'Vtime',         &
     &                          Rclock%DateNumber,                      &
     &                        MODELS(Iroms)%ImportField(id)%Vtime(is),  &
     &                          start=(/Iroms,id,record/),              &
     &                          total=(/1,1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_READ
            RETURN
          END IF
          CALL netcdf_get_time (ng, iNLM, AttFileName, 'Tmin',          &
     &                          Rclock%DateNumber,                      &
     &                          MODELS(Iroms)%ImportField(id)%Tmin,     &
     &                          start=(/Iroms,id,record/),              &
     &                          total=(/1,1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_READ
            RETURN
          END IF
!
          CALL netcdf_get_time (ng, iNLM, AttFileName, 'Tmax',          &
     &                          Rclock%DateNumber,                      &
     &                          MODELS(Iroms)%ImportField(id)%Tmax,     &
     &                          start=(/Iroms,id,record/),              &
     &                          total=(/1,1,1/))
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_READ
            RETURN
          END IF
        END IF
# endif
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

# ifdef TIME_INTERP_NOT_WORKING
!
!  Retrieve custom Attribute Package.
!
          CALL ESMF_AttributeGetAttPack (field,                         &
     &                                   'CustomConvention',            &
     &                                   'General',                     &
!!   &                                   'Instance',                    &
     &                                   attpack=AttPack,               &
     &                                   isPresent=IsPresent,           &
     &                                   rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
!  Get field custom attribute for field for time interpolation.
!
          CALL ESMF_AttributeGet (field,                                &
     &                            name='TimeInterp',                    &
     &                            valueList=AttValues,                  &
     &                            attpack=AttPack,                      &
     &                            isPresent=IsPresent,                  &
     &                            rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
# endif
!
!  Load import data into ROMS component variable.
# ifdef TIME_INTERP
!  If time interpolating in ROMS kernel, loaded import data into
!  snapshot storage arrays so time interpolating is carry out.
!  It is a generic strategy for the case that coupling interval
!  is greater than ROMS time-step size. Usually, time persisting
!  of coupling data may alter ocean solution. For example, it may
!  affect the ocean circulation/energetics if atmospheric forcing
!  is persisted during infrequent coupling (like every 3, 6, or
!  24 hours and so on).
# endif
!
          LoadIt=.TRUE.
          scale      =MODELS(Iroms)%ImportField(id)%scale_factor
          add_offset =MODELS(Iroms)%ImportField(id)%add_offset
          Tindex     =MODELS(Iroms)%ImportField(id)%Tindex
# ifdef TIME_INTERP
          Tmin       =MODELS(Iroms)%ImportField(id)%Tmin
          Tmax       =MODELS(Iroms)%ImportField(id)%Tmax
          Tstr       =MODELS(Iroms)%ImportField(id)%Tstr
          Tend       =MODELS(Iroms)%ImportField(id)%Tend
          MyTintrp(1)=MODELS(Iroms)%ImportField(id)%Tintrp(1)
          MyTintrp(2)=MODELS(Iroms)%ImportField(id)%Tintrp(2)
          MyVtime(1) =MODELS(Iroms)%ImportField(id)%Vtime(1)
          MyVtime(2) =MODELS(Iroms)%ImportField(id)%Vtime(2)
          MyDate(1)  =MODELS(Iroms)%ImportField(id)%DateString(1)
          MyDate(2)  =MODELS(Iroms)%ImportField(id)%DateString(2)
# endif
!
!  Set ROMS momentum fluxes and tracer flux scales to kinematic values.
!  Recall, that all the fluxes are kinematic.
!
          FreshWaterScale=1.0_dp/rho0           ! Kg m-2 s-1 to m/s
          StressScale=1.0_dp/rho0               ! Pa=N m-2 to m2/s2
          TracerFluxScale=1.0_dp/(rho0*Cp)      ! Watts m-2 to C m/s
!
          Fval=ptr2d(IstrR,JstrR)
          MyFmin(1)= MISSING_dp
          MyFmax(1)=-MISSING_dp
          MyFmin(2)= MISSING_dp
          MyFmax(2)=-MISSING_dp
!
          SELECT CASE (TRIM(ADJUSTL(ImportNameList(ifld))))

# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Surface air pressure or mean sea level pressure (mb).
!
            CASE ('psfc', 'Pair', 'Pmsl')
              romsScale=scale
              ifield=idPair
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%PairG(i,j,Tindex)=Fval
#  else
                  FORCES(ng)%Pair(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%Pair)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%Pair)
              END IF
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM        || \
    (defined SHORTWAVE   && defined ANA_SRFLUX    && defined ALBEDO)
!
!  Surface air temperature (Celsius).
!
            CASE ('tsfc', 'Tair')
              romsScale=scale
              ifield=idTair
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%TairG(i,j,Tindex)=Fval
#  else
                  FORCES(ng)%Tair(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%Tair)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%Tair)
              END IF
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM
!
!  Surface air relative humidity (percentage).
!
            CASE ('Qair')
              romsScale=scale
              ifield=idQair
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%HairG(i,j,Tindex)=Fval
#  else
                  FORCES(ng)%Hair(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%Hair)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%Hair)
              END IF
#  endif
# endif
# if defined BULK_FLUXES
!
!  Surface air specific humidity (kg kg-1).
!
            CASE ('Hair', 'qsfc')
              romsScale=scale
              ifield=idQair
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%HairG(i,j,Tindex)=Fval
#  else
                  FORCES(ng)%Hair(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%Hair)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%Hair)
              END IF
#  endif
# endif
# if defined BULK_FLUXES
!
!  Surface net longwave radiation (Celcius m s-1).
!
            CASE ('lwrd', 'LWrad')
              romsScale=TracerFluxScale
              ifield=idLrad
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%lrflxG(i,j,Tindex)=Fval
#  else
                  FORCES(ng)%lrflx(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%lrflx)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%lrflx)
              END IF
#  endif
# endif
# if defined BULK_FLUXES && defined LONGWAVE_OUT
!
!  Surface downward longwave radiation (Celcius m s-1).  ROMS will
!  substract the outgoing IR from model sea surface temperature.
!
            CASE ('dlwr', 'dLWrad', 'lwrad_down')
              romsScale=TracerFluxScale
              ifield=idLdwn
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%lrflxG(i,j,Tindex)=Fval
#  else
                  FORCES(ng)%lrflx(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%lrflx)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%lrflx)
              END IF
#  endif
# endif
# if defined BULK_FLUXES
!
!  Rain fall rate (kg m-2 s-1).
!
            CASE ('prec', 'rain')
              romsScale=scale
              ifield=idrain
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%rainG(i,j,Tindex)=Fval
#  else
                  FORCES(ng)%rain(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%rain)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%rain)
              END IF
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM
!
!  Surface eastward wind component (m s-1). Imported wind component
!  is at RHO-points.
!
            CASE ('wndu', 'Uwind')
              IF (.not.allocated(Uwind)) THEN
                allocate ( Uwind(LBi:UBi,LBj:UBj) )
                Uwind=MISSING_dp
              END IF
              got_wind(1)=.TRUE.
              romsScale=scale
              ifield=idUair
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%UwindG(i,j,Tindex)=Fval
#  else
                  Uwind(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    Uwind)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              Uwind)
              END IF
#  endif
# endif
# if defined BULK_FLUXES || defined ECOSIM
!
!  Surface northward wind component (m s-1). Imported wind component
!  is at RHO-points.
!
            CASE ('wndv', 'Vwind')
              IF (.not.allocated(Vwind)) THEN
                allocate ( Vwind(LBi:UBi,LBj:UBj) )
                Vwind=MISSING_dp
              END IF
              got_wind(2)=.TRUE.
              romsScale=scale
              ifield=idVair
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%VwindG(i,j,Tindex)=Fval
#  else
                  Vwind(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    Vwind)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              Vwind)
              END IF
#  endif
# endif
# if defined SHORTWAVE
!
!  Surface solar shortwave radiation (Celsius m s-1).
!
            CASE ('swrd', 'swrad', 'SWrad', 'SWrad_daily')
              romsScale=TracerFluxScale
              ifield=idSrad
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%srflxG(i,j,Tindex)=Fval
#  else
                  FORCES(ng)%srflx(i,j)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%srflx)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%srflx)
              END IF
#  endif
# endif
# if !defined BULK_FLUXES
!
!  Net longwave radiation flux(W m-2). Used for debugging and plotting
!  purposes to check the fluxes used for the computation of the surface
!  net heat flux in NUOPC cap file "esmf_atm.F".
!
            CASE ('lwr', 'LWrad')
              romsScale=TracerFluxScale
              ifield=idLrad
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  FORCES(ng)%lrflx(i,j)=Fval
                END DO
              END DO
!
!  Surface downward longwave radiation flux(W m-2). Used for debugging
!  and plotting purposes to check the fluxes used for the computation
!  of the surface net heat flux in NUOPC cap file "esmf_atm.F".
!
            CASE ('dlwr', 'dLWrad', 'lwrad_down')
              romsScale=TracerFluxScale
              ifield=idLdwn
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  FORCES(ng)%lrflx(i,j)=Fval
                END DO
              END DO
!
!  Surface latent heat flux (W m-2). Used for plotting and debugging
!  purposes (DebugLevel=3) to check the components of the net surface
!  net heat flux computation.
!
            CASE ('latent', 'LHfx')
              romsScale=TracerFluxScale
              gtype=r2dvar
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  FORCES(ng)%lhflx(i,j)=Fval
                END DO
              END DO
!
!  Surface sensible heat flux (W m-2). Used for plotting and debugging
!  purposes (DebugLevel=3) to check the components of the net surface
!  net heat flux computation.
!
            CASE ('sensible', 'SHfx')
              romsScale=TracerFluxScale
              gtype=r2dvar
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  FORCES(ng)%shflx(i,j)=Fval
                END DO
              END DO
!
!  Surface net heat flux (Celsius m s-1).
!
            CASE ('nflx', 'shflux')
              romsScale=TracerFluxScale
              ifield=idTsur(itemp)
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%stfluxG(i,j,Tindex,itemp)=Fval
#  else
                  FORCES(ng)%stflux(i,j,itemp)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%stflux(:,:,itemp))
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%stflux(:,:,itemp))
              END IF
#  endif
# endif
# if !defined BULK_FLUXES && defined SALINITY
!
!  Surface net freshwater flux: E-P (m s-1).
!
            CASE ('sflx', 'swflux')
              romsScale=FreshWaterScale
              ifield=idTsur(isalt)
              gtype=r2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%stfluxG(i,j,Tindex,isalt)=Fval
#  else
                  FORCES(ng)%stflux(i,j,isalt)=Fval
#  endif
                END DO
              END DO
#  ifndef TIME_INTERP
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%stflux(:,:,isalt))
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%stflux(:,:,isalt))
              END IF
#  endif
# endif
# if !defined BULK_FLUXES
!
!  Surface eastward wind stress component (m2 s-2). Imported stress
!  component is at RHO-points.
!
            CASE ('taux', 'sustr')
              IF (.not.allocated(Ustress)) THEN
                allocate ( Ustress(LBi:UBi,LBj:UBj) )
                Ustress=MISSING_dp
              END IF
              got_stress(1)=.TRUE.
              romsScale=StressScale
              ifield=idUsms
              gtype=u2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%sustrG(i,j,Tindex)=Fval
#  else
                  Ustress(i,j)=Fval
#  endif
                END DO
              END DO
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    Ustress)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              Ustress)
              END IF
# endif
# if !defined BULK_FLUXES
!
!  Surface northward wind stress component (m2 s-2). Imported stress
!  component is at RHO-points.
!
            CASE ('tauy', 'svstr')
              IF (.not.allocated(Vstress)) THEN
                allocate ( Vstress(LBi:UBi,LBj:UBj) )
                Vstress=MISSING_dp
              END IF
              got_stress(2)=.TRUE.
              romsScale=StressScale
              ifield=idVsms
              gtype=v2dvar
              Tindex=3-Iinfo(8,ifield,ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
#  ifdef TIME_INTERP
                  FORCES(ng)%svstrG(i,j,Tindex)=Fval
#  else
                  Vstress(i,j)=Fval
#  endif
                END DO
              END DO
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    Vstress)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              Vstress)
              END IF
# endif
# if defined WIND_MINUS_CURRENT && !defined BULK_FLUXES
!
!  Surface air density (kg/m3).
!
            CASE ('RhoAir')
              IF (.not.allocated(RhoAir)) THEN
                allocate ( RhoAir(LBi:UBi,LBj:UBj) )
                RhoAir=MISSING_dp
              END IF
              got_RhoAir=.TRUE.
              romsScale=scale
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  RhoAir(i,j)=Fval
                END DO
              END DO
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    RhoAir)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              RhoAir)
              END IF
!
!  Eastward wind component (m s-1) at surface boundary layer. Imported
!  wind component is at RHO-points.
!
            CASE ('Uwind_sbl')
              IF (.not.allocated(Xwind)) THEN
                allocate ( Xwind(LBi:UBi,LBj:UBj) )
                Xwind=MISSING_dp
              END IF
              got_wind_sbl(1)=.TRUE.
              romsScale=scale
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  Xwind(i,j)=Fval
                END DO
              END DO
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    Xwind)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              Xwind)
              END IF
!
!  Northward wind component (m s-1) at surface boundary layer. Imported
!  wind component is at RHO-points.
!
            CASE ('Vwind_sbl')
              IF (.not.allocated(Ywind)) THEN
                allocate ( Ywind(LBi:UBi,LBj:UBj) )
                Ywind=MISSING_dp
              END IF
              got_wind_sbl(2)=.TRUE.
              romsScale=scale
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  Ywind(i,j)=Fval
                END DO
              END DO
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    Ywind)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              Ywind)
              END IF
!
!  Surface frictional wind magnitude (m s-1) from similarity theory.
!  Imported wind magnitude is at RHO-points.
!
            CASE ('Wstar')
              IF (.not.allocated(Wstar)) THEN
                allocate ( Wstar(LBi:UBi,LBj:UBj) )
                Wstar=MISSING_dp
              END IF
              got_Wstar=.TRUE.
              romsScale=scale
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (ABS(ptr2d(i,j)).lt.TOL_dp) THEN
                    Fval=scale*ptr2d(i,j)+add_offset
                  ELSE
                    Fval=0.0_dp
                  END IF
                  MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                  MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  Fval=Fval*romsScale
                  MyFmin(2)=MIN(MyFmin(2),Fval)
                  MyFmax(2)=MAX(MyFmax(2),Fval)
                  Wstar(i,j)=Fval
                END DO
              END DO
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    Wstar)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              Wstar)
              END IF
#  endif
!
!  Import field not found.
!
            CASE DEFAULT
              IF (localPET.eq.0) THEN
                WRITE (cplout,10) TRIM(ImportNameList(ifld)),           &
     &                            TRIM(Time_CurrentString),             &
     &                            TRIM(CinpName)
              END IF
              exit_flag=9
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                     MyFile)) THEN
                rc=ESMF_RC_NOT_FOUND
                RETURN
              END IF
          END SELECT
!
!  Print pointer information.
!
          IF (DebugLevel.eq.4) THEN
            WRITE (cplout,20) localPET, localDE,                        &
     &                        LBOUND(ptr2d,DIM=1), UBOUND(ptr2d,DIM=1), &
     &                        LBOUND(ptr2d,DIM=2), UBOUND(ptr2d,DIM=2), &
     &                        IstrR, IendR, JstrR, JendR
          END IF
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
!  Write out import field information.
!
        IF ((DebugLevel.ge.0).and.(localPET.eq.0)) THEN
          WRITE (cplout,30) TRIM(ImportNameList(ifld)),                 &
# ifdef TIME_INTERP
     &                      TRIM(MyDate(Tindex)), ng,                   &
     &                      Fmin(1), Fmax(1), Tindex
# else
     &                      TRIM(Time_CurrentString), ng,               &
     &                      Fmin(1), Fmax(1)
# endif
          IF (romsScale.ne.1.0_dp) THEN
            WRITE (cplout,40) Fmin(2), Fmax(2),                         &
     &                        ' romsScale = ', romsScale
          ELSE IF (add_offset.ne.0.0_dp) THEN
            WRITE (cplout,40) Fmin(2), Fmax(2),                         &
     &                        ' AddOffset = ', add_offset
          END IF
        END IF

# ifdef TIME_INTERP
!
!  Load ROMS metadata information needed for time interpolation and
!  reporting.
!
        IF (Loadit) THEN
          Linfo(1,ifield,ng)=.TRUE.                          ! Lgrided
          Linfo(3,ifield,ng)=.FALSE.                         ! Lonerec
          Iinfo(1,ifield,ng)=gtype
          Iinfo(8,ifield,ng)=Tindex
          Finfo(1,ifield,ng)=Tmin
          Finfo(2,ifield,ng)=Tmax
          Finfo(3,ifield,ng)=Tstr
          Finfo(4,ifield,ng)=Tend
          Finfo(8,ifield,ng)=Fmin(1)
          Finfo(9,ifield,ng)=Fmax(1)
          Vtime(Tindex,ifield,ng)=MyVtime(Tindex)
          Tintrp(Tindex,ifield,ng)=MyTintrp(Tindex)*86400.0_dp
        END IF
# endif
!
!  Debugging: write out import field into NetCDF file.
!
        IF ((DebugLevel.ge.3).and.                                      &
     &      MODELS(Iroms)%ImportField(ifld)%debug_write) THEN
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

# if defined BULK_FLUXES || defined ECOSIM
!
!  If applicable, rotate wind components to ROMS curvilinear grid.
!
      IF (got_wind(1).and.got_wind(2)) THEN
        CALL ROMS_Rotate (ng, tile, geo2grid_rho,                       &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Uwind, Vwind,                                 &
     &                    FORCES(ng)%Uwind, FORCES(ng)%Vwind)
        deallocate (Uwind)
        deallocate (Vwind)
      END IF
#  endif
#  if !defined BULK_FLUXES
!
!  If applicable, rotate wind stress components to ROMS curvilinear
!  grid.
!
      IF (got_stress(1).and.got_stress(2)) THEN
        CALL ROMS_Rotate (ng, tile, geo2grid,                            &
     &                    LBi, UBi, LBj, UBj,                            &
     &                    Ustress, Vstress,                              &
     &                    FORCES(ng)%sustr, FORCES(ng)%svstr)
        deallocate (Ustress)
        deallocate (Vstress)
      END IF
# endif
# if defined WIND_MINUS_CURRENT && !defined BULK_FLUXES
!
!  If applicable, compute surface wind stress components. The surface
!  ocean currents are substracted to the wind.
!
!  The wind stress component are computed as:
!
!          taux/rho0 = RhoAir * Cd * Wrel * Urel
!          tauy/rho0 = RhoAir * Cd * Wrel * Vrel
!  where
!          Cd = Wstr**2 / Wmag**2
!
!  so the magnitude is diminished by the weaker relative (wind minus
!  current) components. The coupling is incompleate becasue there is
!  not feeback to the atmosphere (wind is not modified by currents).
!
      MyFmin= MISSING_dp
      MyFmax=-MISSING_dp
!
      IF (got_RhoAir.and.got_Wstar.and.                                 &
     &    got_wind_sbl(1).and.got_wind_sbl(2)) THEN
        IF (.not.allocated(Uwrk)) THEN
          allocate ( Uwrk(LBi:UBi,LBj:UBj) )
          Uwrk=MISSING_dp
        END IF
        IF (.not.allocated(Vwrk)) THEN
          allocate ( Vwrk(LBi:UBi,LBj:UBj) )
          Vwrk=MISSING_dp
        END IF
!
        CALL ROMS_Rotate (ng, tile, grid2geo_rho,                       &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    OCEAN(ng)%u(:,:,N(ng),nstp(ng)),              &
     &                    OCEAN(ng)%v(:,:,N(ng),nstp(ng)),              &
     &                    Uwrk, Vwrk)         ! rotated currents to E-N
!
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            romsScale=StressScale                   ! m3/kg
            Urel=Xwind(i,j)-Uwrk(i,j)               ! relative wind:
            Vrel=Ywind(i,j)-Vwrk(i,j)               ! wind minus current
            Wmag=SQRT(Xwind(i,j)*Xwind(i,j)+                            &
     &                Ywind(i,j)*Ywind(i,j))        ! ATM wind magnitude
            Wrel=SQRT(Urel*Urel+Vrel*Vrel)          ! relative magmitude
            cff1=romsScale*RhoAir(i,j)
            cff2=Wstar(i,j)*Wstar(i,j)/(Wmag*Wmag+eps)
            cff3=cff1*cff2*Wrel                     ! m/s
            Uwrk(i,j)=cff3*Urel                     ! m2/s2
            Vwrk(i,j)=cff3*Vrel                     ! m2/s2
            MyFmin(1)=MIN(MyFmin(1),Uwrk(i,j))
            MyFmin(2)=MIN(MyFmin(2),Vwrk(i,j))
            MyFmax(1)=MAX(MyFmax(1),Uwrk(i,j))
            MyFmax(2)=MAX(MyFmax(2),Vwrk(i,j))
          END DO
        END DO
        deallocate (RhoAir)
        deallocate (Wstar)
        deallocate (Xwind)
        deallocate (Ywind)
!                                            ! rotate stress to grid
        CALL ROMS_Rotate (ng, tile, geo2grid,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Uwrk, Vwrk,                                   &
     &                    FORCES(ng)%sustr,                             &
     &                    FORCES(ng)%svstr)
        deallocate (Uwrk)
        deallocate (Vwrk)
!
!  Report computed wind stress minimum and maximum values.
!
        IF (DebugLevel.ge.0) THEN
          CALL ESMF_VMAllReduce (vm,                                    &
     &                           sendData=MyFmin,                       &
     &                           recvData=Fmin,                         &
     &                           count=2,                               &
     &                           reduceflag=ESMF_REDUCE_MIN,            &
     &                           rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_VMAllReduce (vm,                                    &
     &                           sendData=MyFmax,                       &
     &                           recvData=Fmax,                         &
     &                           count=2,                               &
     &                           reduceflag=ESMF_REDUCE_MAX,            &
     &                           rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          IF (localPET.eq.0) THEN
            WRITE (cplout,60) 'sustr',                                  &
     &                        TRIM(Time_CurrentString), ng,             &
     &                        Fmin(1)/StressScale,                      &
     &                        Fmax(1)/StressScale
            WRITE (cplout,40) Fmin(1), Fmax(1),                         &
     &                        ' romsScale = ', StressScale
!
            WRITE (cplout,60) 'svstr',                                  &
     &                        TRIM(Time_CurrentString), ng,             &
     &                        Fmin(2)/StressScale,                      &
     &                        Fmax(2)/StressScale
            WRITE (cplout,40) Fmin(2), Fmax(2),                         &
     &                        ' romsScale = ', StressScale
          END IF
        END IF
      END IF
# endif
!
!  Deallocate local arrays.
!
      IF (allocated(ImportNameList)) deallocate (ImportNameList)
!
!  Update ROMS import calls counter.
!
      IF (ImportCount.gt.0) THEN
        MODELS(Iroms)%ImportCalls=MODELS(Iroms)%ImportCalls+1
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_Import',             &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      IF (DebugLevel.gt.0) FLUSH (cplout)
!
  10  FORMAT (/,3x,' ROMS_Import - unable to find option to import: ',  &
     &        a,t72,a,/,18x,'check ''Import(roms)'' in input script: ', &
     &        a)
  20  FORMAT (18x,'PET/DE [',i3.3,'/',i2.2,'],  Pointer Size: ',4i8,    &
     &        /,36x,'Tiling Range: ',4i8)
  30  FORMAT (3x,' ROMS_Import - ESMF: importing field ''',a,'''',      &
     &        t72,a,2x,'Grid ',i2.2,                                    &
# ifdef TIME_INTERP
     &        /,19x,'(InpMin = ', 1p,e15.8,0p,' InpMax = ',1p,e15.8,0p, &
     &        '  SnapshotIndex = ',i1,')')
# else
     &        /,19x,'(InpMin = ', 1p,e15.8,0p,' InpMax = ',1p,e15.8,0p, &
     &        ')')
# endif
  40  FORMAT (19x,'(OutMin = ', 1p,e15.8,0p,' OutMax = ',1p,e15.8,0p,   &
     &        1x,a,1p,e15.8,0p,')')
  50  FORMAT ('roms_',i2.2,'_import_',a,'_',i4.4,2('-',i2.2),'_',       &
     &        i2.2,2('.',i2.2),'.nc')
  60  FORMAT (3x,' ROMS_Import - ESMF: computing field ''',a,'''',      &
     &        t72,a,2x,'Grid ',i2.2,                                    &
     &        /,19x,'(InpMin = ', 1p,e15.8,0p,' InpMax = ',1p,e15.8,0p, &
     &        ')')
!
      RETURN
      END SUBROUTINE ROMS_Import
!
      SUBROUTINE ROMS_Export (ng, model, rc)
!
!=======================================================================
!                                                                      !
!  Exports ROMS fields to other coupled gridded components.            !
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
      logical :: get_barotropic
      logical :: get_SurfaceCurrent
!
      integer :: Istr,  Iend,  Jstr,  Jend
      integer :: IstrR, IendR, JstrR, JendR
      integer :: LBi, UBi, LBj, UBj
      integer :: ExportCount
      integer :: localDE, localDEcount, localPET, tile
      integer :: year, month, day, hour, minutes, seconds, sN, SD
      integer :: ifld, i, is, j
!
      real (dp) :: Fmin(1), Fmax(1), Fval, MyFmin(1), MyFmax(1)
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
!
      real (dp), allocatable  :: Ubar(:,:), Vbar(:,:)
      real (dp), allocatable  :: Usur(:,:), Vsur(:,:)
!
      character (len=22)      :: Time_CurrentString

      character (len=:), allocatable :: fldname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_Export"

      character (ESMF_MAXSTR) :: cname, ofile
      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
!
      TYPE (ESMF_Field) :: field
      TYPE (ESMF_Time)  :: CurrentTime
      TYPE (ESMF_VM)    :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_Export',             &
     &                          ', PET', PETrank
        FLUSH (trac)
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
!  Get number of local decomposition elements (DEs). Usually, a single
!  DE is associated with each Persistent Execution Thread (PETs). Thus,
!  localDEcount=1.
!
      CALL ESMF_GridGet (MODELS(Iroms)%grid(ng),                        &
     &                   localDECount=localDEcount,                     &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set horizontal tile bounds.
!
      tile=localPET
!
      LBi=BOUNDS(ng)%LBi(tile)       ! lower bound I-direction
      UBi=BOUNDS(ng)%UBi(tile)       ! upper bound I-direction
      LBj=BOUNDS(ng)%LBj(tile)       ! lower bound J-direction
      UBj=BOUNDS(ng)%UBj(tile)       ! upper bound J-direction
!
      IstrR=BOUNDS(ng)%IstrR(tile)   ! Full range I-starting (RHO)
      IendR=BOUNDS(ng)%IendR(tile)   ! Full range I-ending   (RHO)
      JstrR=BOUNDS(ng)%JstrR(tile)   ! Full range J-starting (RHO)
      JendR=BOUNDS(ng)%JendR(tile)   ! Full range J-ending   (RHO)
!
      Istr=BOUNDS(ng)%Istr(tile)     ! Full range I-starting (PSI, U)
      Iend=BOUNDS(ng)%Iend(tile)     ! Full range I-ending   (PSI)
      Jstr=BOUNDS(ng)%Jstr(tile)     ! Full range J-starting (PSI, V)
      Jend=BOUNDS(ng)%Jend(tile)     ! Full range J-ending   (PSI)
!
!-----------------------------------------------------------------------
!  Get current time.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (ClockInfo(Iroms)%Clock,                       &
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
      CALL ESMF_StateGet (MODELS(Iroms)%ExportState(ng),                &
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
      CALL ESMF_StateGet (MODELS(Iroms)%ExportState(ng),                &
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
      get_barotropic=.TRUE.
      get_SurfaceCurrent=.TRUE.
!
      FLD_LOOP : DO ifld=1,ExportCount
!
!   Get field from export state.
!
        CALL ESMF_StateGet (MODELS(Iroms)%ExportState(ng),              &
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
!
!  Initialize pointer to missing value.
!
          ptr2d=MISSING_dp
!
!  Load field data into export state. Notice that all export fields
!  are kept as computed by ROMS. The imported component does the
!  proper scaling, physical units conversion, and other manipulations.
!  It is done to avoid applying such transformations twice.
!
          SELECT CASE (TRIM(ADJUSTL(ExportNameList(ifld))))
!
!  Sea surface temperature (C).
# if defined EXCLUDE_SPONGE && \
    (defined DATA_COUPLING  && !defined ANA_SPONGE)
!  If using a diffusion sponge, remove the SST points in the sponge area
!  to supress the spurious influence of open boundary conditions in the
!  computation of the net heat flux. The SST values in the sponge are
!  from the large scale DATA component in the merged ocean/data field
!  imported by the atmosphere model.
# endif
!
            CASE ('sst', 'SST')
              MyFmin(1)= MISSING_dp
              MyFmax(1)=-MISSING_dp
              DO j=JstrR,JendR
                DO i=IstrR,IendR
# if defined EXCLUDE_SPONGE && \
    (defined DATA_COUPLING  && !defined ANA_SPONGE)
                  IF (LtracerSponge(itemp,ng).and.                      &
     &                MIXING(ng)%diff_factor(i,j).gt.1.0_dp) THEN
                    Fval=MISSING_dp
                  ELSE
                    Fval=OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)
#  ifdef MASKING
                    IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                      MyFmin(1)=MIN(MyFmin(1),Fval)
                      MyFmax(1)=MAX(MyFmax(1),Fval)
                    END IF
#  else
                    MyFmin(1)=MIN(MyFmin(1),Fval)
                    MyFmax(1)=MAX(MyFmax(1),Fval)
#  endif
                  END IF
# else
                  Fval=OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)
#  ifdef MASKING
                  IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                    MyFmin(1)=MIN(MyFmin(1),Fval)
                    MyFmax(1)=MAX(MyFmax(1),Fval)
                  END IF
#  else
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
#  endif
# endif
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Sea surface height (m).
!
            CASE ('ssh', 'SSH')
              MyFmin(1)=1.0_dp
              MyFmax(1)=0.0_dp
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  Fval=OCEAN(ng)%zeta(i,j,knew(ng))
# ifdef MASKING
                  IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                    MyFmin(1)=MIN(MyFmin(1),Fval)
                    MyFmax(1)=MAX(MyFmax(1),Fval)
                  END IF
# else
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
# endif
                  ptr2d(i,j)=Fval
                END DO
              END DO
!
!  Depth-integrated (barotropic) currents (m/s) at interior RHO-points
!  (East/North direction).
!
            CASE ('Ubar', 'Vbar')
              IF (FoundError(assign_string(fldname,                     &
     &                                     ExportNameList(ifld)),       &
     &            NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_NOT_FOUND
                RETURN
              END IF
!
              IF (get_barotropic) THEN
                get_barotropic=.FALSE.
                IF (.not.allocated(Ubar)) THEN
                  allocate ( Ubar(LBi:UBi,LBj:UBj) )
                  Ubar=MISSING_dp
                END IF
                IF (.not.allocated(Vbar)) THEN
                  allocate ( Vbar(LBi:UBi,LBj:UBj) )
                  Vbar=MISSING_dp
                END IF
                CALL ROMS_Rotate (ng, tile, grid2geo_rho,               &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OCEAN(ng)%ubar(:,:,knew(ng)),         &
     &                            OCEAN(ng)%vbar(:,:,knew(ng)),         &
     &                            Ubar, Vbar)
              END IF
!
              IF (fldname.eq.'Ubar') THEN
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    Fval=Ubar(i,j)
                    MyFmin(1)=MIN(MyFmin(1),Fval)
                    MyFmax(1)=MAX(MyFmax(1),Fval)
                    ptr2d(i,j)=Fval
                  END DO
                END DO
                deallocate (Ubar)
              ELSE
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    Fval=Vbar(i,j)
                    MyFmin(1)=MIN(MyFmin(1),Fval)
                    MyFmax(1)=MAX(MyFmax(1),Fval)
                    ptr2d(i,j)=Fval
                  END DO
                END DO
                deallocate (Vbar)
              END IF
!
!  Surface currents (m/s) at interior RHO-points (East/North direction).
!
            CASE ('Usur', 'Vsur')
              IF (FoundError(assign_string(fldname,                     &
     &                                     ExportNameList(ifld)),       &
     &            NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_NOT_FOUND
                RETURN
              END IF
!
              IF (get_SurfaceCurrent) THEN
                get_SurfaceCurrent=.FALSE.
                IF (.not.allocated(Ubar)) THEN
                  allocate ( Usur(LBi:UBi,LBj:UBj) )
                  Usur=MISSING_dp
                END IF
                IF (.not.allocated(Vbar)) THEN
                  allocate ( Vsur(LBi:UBi,LBj:UBj) )
                  Vsur=MISSING_dp
                END IF
                CALL ROMS_Rotate (ng, tile, grid2geo_rho,              &
     &                            LBi, UBi, LBj, UBj,                  &
     &                            OCEAN(ng)%u(:,:,N(ng),nstp(ng)),     &
     &                            OCEAN(ng)%v(:,:,N(ng),nstp(ng)),     &
     &                            Usur, Vsur)
              END IF
!
              IF (fldname.eq.'Usur') THEN
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    Fval=Usur(i,j)
                    MyFmin(1)=MIN(MyFmin(1),Fval)
                    MyFmax(1)=MAX(MyFmax(1),Fval)
                    ptr2d(i,j)=Fval
                  END DO
                END DO
                deallocate (Usur)
              ELSE
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    Fval=Vsur(i,j)
                    MyFmin(1)=MIN(MyFmin(1),Fval)
                    MyFmax(1)=MAX(MyFmax(1),Fval)
                    ptr2d(i,j)=Fval
                  END DO
                END DO
                deallocate (Vsur)
              END IF
!
!  Bathymetry (m). It can be time dependent due sediment morphology.
!
            CASE ('bath')
              MyFmin(1)=1.0_dp
              MyFmax(1)=0.0_dp
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  Fval=GRID(ng)%h(i,j)
                  MyFmin(1)=MIN(MyFmin(1),Fval)
                  MyFmax(1)=MAX(MyFmax(1),Fval)
                  ptr2d(i,j)=Fval
                END DO
              END DO

# if defined MASKING && defined WET_DRY
!
!  Update wet point land/sea mask, if differs from static mask.
!
            CASE ('msk')
              MyFmin(1)=1.0_dp
              MyFmax(1)=0.0_dp
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                    IF (GRID(ng)%rmask(i,j).ne.                         &
     &                  GRID(ng)%rmask_wet(i,j)) THEN
                      ptr2d(i,j)=GRID(ng)%rmask_wet(i,j)
                    ELSE
                      ptr2d(i,j)=GRID(ng)%rmask(i,j)
                    END IF
                    MyFmin(1)=MIN(MyFmin(1),ptr2d(i,j))
                    MyFmax(1)=MAX(MyFmax(1),ptr2d(i,j))
                  END IF
                END DO
              END DO
# endif
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
        IF (localPET.eq.0) THEN
          WRITE (cplout,20) TRIM(ExportNameList(ifld)),                 &
     &                      TRIM(Time_CurrentString), ng,               &
     &                      Fmin(1), Fmax(1)
        END IF
!
!  Debugging: write out field into a NetCDF file.
!
        IF ((DebugLevel.ge.3).and.                                      &
     &      MODELS(Iroms)%ExportField(ifld)%debug_write) THEN
          WRITE (ofile,30) ng, TRIM(ExportNameList(ifld)),              &
                           year, month, day, hour, minutes, seconds
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
!  Update ROMS export calls counter.
!
      IF (ExportCount.gt.0) THEN
        MODELS(Iroms)%ExportCalls=MODELS(Iroms)%ExportCalls+1
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_Export',             &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      FLUSH (cplout)
!
  10  FORMAT (/,3x,' ROMS_Export - unable to find option to export: ',  &
     &        a,/,18x,'check ''Export(roms)'' in input script: ',a)
  20  FORMAT (3x,' ROMS_Export - ESMF:  exporting field ''',a,'''',     &
     &        t72,a,2x,'Grid ',i2.2,/,                                  &
     &        18x,'(OutMin = ', 1p,e15.8,0p,' OutMax = ',1p,e15.8,0p,   &
     &        ')')
  30  FORMAT ('roms_',i2.2,'_export_',a,'_',i4.4,2('-',i2.2),'_',       &
     &        i2.2,2('.',i2.2),'.nc')

      RETURN
      END SUBROUTINE ROMS_Export
!
      SUBROUTINE ROMS_Rotate (ng, tile, Lrotate,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Uinp, Vinp,                               &
     &                        Uout, Vout)
!
!=======================================================================
!                                                                      !
!  It rotates exchanged vector components from computational grid to   !
!  geographical EAST and NORTH directions or vice versa acccording to  !
!  Lrotate flag:                                                       !
!                                                                      !
!    Lrotate =  geo2grid_rho       RHO-points rotation                 !
!    Lrotate =  grid2geo_rho       Exporting interior RHO-points       !
!    Lrotate =  geo2grid           U- and V-points staggered rotation  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ng, tile, Lrotate
      integer, intent(in)  :: LBi, UBi, LBj, UBj
!
      real (dp), intent(in)  :: Uinp(LBi:UBi,LBj:UBj)
      real (dp), intent(in)  :: Vinp(LBi:UBi,LBj:UBj)
      real (r8), intent(out) :: Uout(LBi:UBi,LBj:UBj)
      real (r8), intent(out) :: Vout(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: IstrR, IendR, JstrR, JendR
      integer :: Istr, Iend, Jstr, Jend
!
      real :: Urho, Vrho
!
      real (r8) :: Urot(LBi:UBi,LBj:UBj)
      real (r8) :: Vrot(LBi:UBi,LBj:UBj)
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_Rotate',             &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
!  Set horizontal tile bounds.
!
      IstrR=BOUNDS(ng)%IstrR(tile)   ! Full range I-starting (RHO)
      IendR=BOUNDS(ng)%IendR(tile)   ! Full range I-ending   (RHO)
      JstrR=BOUNDS(ng)%JstrR(tile)   ! Full range J-starting (RHO)
      JendR=BOUNDS(ng)%JendR(tile)   ! Full range J-ending   (RHO)
!
      Istr=BOUNDS(ng)%Istr(tile)     ! Full range I-starting (PSI, U)
      Iend=BOUNDS(ng)%Iend(tile)     ! Full range I-ending   (PSI)
      Jstr=BOUNDS(ng)%Jstr(tile)     ! Full range J-starting (PSI, V)
      Jend=BOUNDS(ng)%Jend(tile)     ! Full range J-ending   (PSI)

# ifdef CURVGRID
!
!-----------------------------------------------------------------------
!  Rotate from geographical (EAST, NORTH) to computational grid
!  directions (ROMS import case).
!-----------------------------------------------------------------------
!
      IF ((Lrotate.eq.geo2grid).or.(lrotate.eq.geo2grid_rho)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Urot(i,j)=Uinp(i,j)*GRID(ng)%CosAngler(i,j)+                &
     &                Vinp(i,j)*GRID(ng)%SinAngler(i,j)
            Vrot(i,j)=Vinp(i,j)*GRID(ng)%CosAngler(i,j)-                &
     &                Uinp(i,j)*GRID(ng)%SinAngler(i,j)
          END DO
        END DO
!
!  There is an option to import the rotated vector to staggered U- and
!  V-locations (arithmetic avererage) or import vector at its native
!  cell center (RHO-points).
!
        IF (Lrotate.eq.geo2grid_rho) THEN             ! RHO-points
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              Uout(i,j)=Urot(i,j)
              Vout(i,j)=Vrot(i,j)
            END DO
          END DO
!
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Uout)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Vout)
          END IF

        ELSE IF (Lrotate.eq.geo2grid) THEN            ! U- and V-points
          DO j=JstrR,JendR
            DO i=Istr,IendR
              Uout(i,j)=0.5_r8*(Urot(i-1,j)+Urot(i,j))
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              Vout(i,j)=0.5_r8*(Vrot(i,j-1)+Vrot(i,j))
            END DO
          END DO

          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Uout)
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Vout)
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Rotate from computational grid to geographical (EAST, NORTH)
!  directions (ROMS Export case: vector at RHO-points).
!-----------------------------------------------------------------------
!
      ELSE IF (Lrotate.eq.grid2geo_rho) THEN
        Uout=0.0_r8
        Vout=0.0_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Urho=0.5_r8*(Uinp(i,j)+Uinp(i+1,j))
            Vrho=0.5_r8*(Vinp(i,j)+Vinp(i,j+1))
            Uout(i,j)=Urho*GRID(ng)%CosAngler(i,j)-                     &
     &                Vrho*GRID(ng)%SinAngler(i,j)
            Vout(i,j)=Vrho*GRID(ng)%CosAngler(i,j)+                     &
     &                Urho*GRID(ng)%SinAngler(i,j)
#  ifdef MASKING
            Uout(i,j)=Uout(i,j)*GRID(ng)%rmask(i,j)
            Vout(i,j)=Vout(i,j)*GRID(ng)%rmask(i,j)
#  endif
          END DO
        END DO
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Uout)
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Vout)
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Uout)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Vout)
        END IF
      END IF
# else
!
!-----------------------------------------------------------------------
!  Otherwise, load unrotated components to staggered location. ROMS grid
!  is not curvilinear (ROMS import case). It is very unlikely to have
!  realistic applications that are not curvilinear and rotated).
!-----------------------------------------------------------------------
!
      IF (Lrotate.eq.geo2grid_rho) THEN               ! RHO-points
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Uout(i,j)=Uinp(i,j)
            Vout(i,j)=Vinp(i,j)
          END DO
        END DO
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Uout)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Vout)
        END IF

      ELSE IF (Lrotate.eq.geo2grid) THEN              ! U- and V-points
        DO j=JstrR,JendR
          DO i=Istr,IendR
            Uout(i,j)=0.5_r8*(Uinp(i-1,j)+Uinp(i,j))
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            Vout(i,j)=0.5_r8*(Vinp(i,j-1)+Vinp(i,j))
          END DO
        END DO
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Uout)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Vout)
        END IF
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Distributed-memory tile (halo) exchange.
!-----------------------------------------------------------------------
!
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Uout, Vout)
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_Rotate',             &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
      END SUBROUTINE ROMS_Rotate
!
#endif
      END MODULE esmf_roms_mod
