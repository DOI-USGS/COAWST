#include "cppdefs.h"
      MODULE esmf_data_mod

#if defined DATA_COUPLING && defined MODEL_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_data.F 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets a DATA component using generic ESMF/NUOPC layer.   !
!  It is used to provide data to ESM components uncoupled fields or to !
!  provide values at locations not imported because of noncoincident   !
!  grids. The data comes from specified input files, usually NetCDF    !
!  files.                                                              !
!                                                                      !
!    DATA_SetServices        Sets DATA component shared-object entry   !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                                                                      !
!    DATA_SetInitializeP1    DATA component phase 1 initialization:    !
!                            sets import and export fields long and    !
!                            short names into its respective state.    !
!                                                                      !
!    DATA_SetInitializeP2    DATA component phase 2 initialization:    !
!                            Initializes component (DATA_initialize),  !
!                            sets component grid (DATA_SetGridArrays), !
!                            and adds fields into import and export    !
!                            into respective states.                   !
!                                                                      !
!    DATA_Initialize         Initializes DATA component IFS structure, !
!                            inquire the contents of source NetCDF     !
!                            files, reads (lon,lat) coordinates, and   !
!                            land/sea mask for each export field, and  !
!                            reads LOWER time-snapshots data.          !
# ifdef TIME_INTERP
!                                                                      !
!    DATA_DataInit           Exports DATA component LOWER time snapshot!
!                            (if exporting the two time levels to ESM  !
!                            destination) or reads UPPER time-snapshot !
!                            (if time interpolating in DATA_Export).   !
# endif
!                                                                      !
!    DATA_SetClock           Sets DATA component date calendar, start  !
!                            and stop times, and coupling interval.    !
!                                                                      !
!    DATA_SetGridArrays      Sets DATA component staggered, horizontal !
!                            grid arrays, grid area, and land/sea mask !
!                            if any.                                   !
!                                                                      !
!    DATA_SetStates          Adds DATA component export and import     !
!                            fields into its respective state.         !
!                                                                      !
!    DATA_ModelAdvance       Advances DATA component for a coupling    !
!                            interval. It calls import and export      !
!                            routines.                                 !
!                                                                      !
!    DATA_SetFinalize        Finalizes DATA component execution.       !
!                                                                      !
!    DATA_Export             Exports DATA fields to other gridded      !
!                            components.                               !
!                                                                      !
!    DATA_TimeInterp         It time interpolates from the snapshots   !
!                            the fields that the DATA component exports!
!                            other ESM components.                     !
!                                                                      !
!    DATA_ncread             Reads DATA component fields to export     !
!                            from source NetCDF files.                 !
!                                                                      !
!    DATA_multifile          Initializes Input File Structure (IFS)    !
!                            and sets various parameters for single    !
!                            and multi-file source data.               !
!                                                                      !
!    DATA_checkfile          Scans the variables of a NetCDF file for  !
!                            the time record variable and gets its     !
!                            range values.                             !
!                                                                      !
!    DATA_inquiry            Inquires source input NetCDF files for    !
!                            export field to process.                  !
!                                                                      !
!    DATA_ncvarcoords        Reads export field longitude, latitude,   !
!                            and land/sea mask, if available.          !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!                                                                      !
!=======================================================================
!
      USE ESMF
      USE NUOPC
      USE NUOPC_Model,                                                  &
     &    NUOPC_SetServices           => SetServices,                   &
     &    NUOPC_Label_Advance         => label_Advance,                 &
# ifdef TIME_INTERP
     &    NUOPC_Label_DataInitialize  => label_DataInitialize,          &
# endif
     &    NUOPC_Label_SetClock        => label_SetClock
!
      USE mod_esmf_esm          ! ESM coupling structures and variables
!
      implicit none
!
      PUBLIC  :: DATA_SetServices
!
      PRIVATE :: DATA_SetInitializeP1
      PRIVATE :: DATA_SetInitializeP2
      PRIVATE :: DATA_Initialize
# ifdef TIME_INTERP
      PRIVATE :: DATA_DataInit
# endif
      PRIVATE :: DATA_SetClock
      PRIVATE :: DATA_SetGridArrays
      PRIVATE :: DATA_SetStates
      PRIVATE :: DATA_ModelAdvance
      PRIVATE :: DATA_SetFinalize
      PRIVATE :: DATA_Export
      PRIVATE :: DATA_TimeInterp
      PRIVATE :: DATA_ncread
      PRIVATE :: DATA_multifile
      PRIVATE :: DATA_checkfile
      PRIVATE :: DATA_inquiry
      PRIVATE :: DATA_ncvarcoords
!
      CONTAINS
!
      SUBROUTINE DATA_SetServices (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets DATA component shared-object entry points for "initialize",    !
!  "run", and "finalize" by using NUOPC generic methods.               !
!                                                                      !
!=======================================================================
!
      implicit none
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
     &  __FILE__//", DATA_SetServices"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_SetServices',        &
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
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Register initialize routines.
!-----------------------------------------------------------------------
!
!  Set routine for Phase 1 initialization (advertise export fields).
!
      CALL NUOPC_CompSetEntryPoint (model,                              &
     &                              methodflag=ESMF_METHOD_INITIALIZE,  &
     &                              phaseLabelList=(/"IPDv00p1"/),      &
     &                              userRoutine=DATA_SetInitializeP1,   &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
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
     &                              userRoutine=DATA_SetInitializeP2,   &
     &                              rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Attach DATA component phase independent specializing methods.
!-----------------------------------------------------------------------

# ifdef TIME_INTERP
!
!  Set routine for export initial/restart fields.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_DataInitialize,  &
     &                           specRoutine=DATA_DataInit,             &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# endif
!
!  Set routine for setting DATA component clock.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_SetClock,        &
     &                           specRoutine=DATA_SetClock,             &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set routine for advancing DATA component.
!
      CALL NUOPC_CompSpecialize (model,                                 &
     &                           specLabel=NUOPC_Label_Advance,         &
     &                           specRoutine=DATA_ModelAdvance,         &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Register DATA component finalize routine.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompSetEntryPoint (model,                           &
     &                                 methodflag=ESMF_METHOD_FINALIZE, &
     &                                 userRoutine=DATA_SetFinalize,    &
     &                                 rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_SetServices',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE DATA_SetServices
!
      SUBROUTINE DATA_SetInitializeP1 (model,                           &
     &                                 ImportState, ExportState,        &
     &                                 clock, rc)
!
!=======================================================================
!                                                                      !
!  DATA component Phase 1 initialization: sets export fields long and  !
!  short names into its respective state.  Currently, the DATA model   !
!  does not need to import fields.                                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc

      TYPE (ESMF_GridComp) :: model
      TYPE (ESMF_State)    :: ImportState
      TYPE (ESMF_State)    :: ExportState
      TYPE (ESMF_Clock)    :: clock
!
!  Local variable declarations.
!
      integer :: id, ifld, localPET, nd, ng
      integer :: Icomp, Nfields, Nfiles
!
      character (len=100) :: CoupledSet, StateLabel
      character (len=240) :: StandardName, ShortName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_SetInitializeP1"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_SetInitializeP1',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry about current node rank.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       localPet=localPET,                         &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Set DATA export state and fields.
!-----------------------------------------------------------------------
!
!  Add DATA export state for connected components.
!
       DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
            DO ng=1,MODELS(Icomp)%Ngrids
              IF (COUPLED(Icomp)%LinkedGrid(ng,Idata)) THEN
                nd=COUPLED(Idata)%DataCoupledSets(ng,Icomp)
                CoupledSet=TRIM(COUPLED(Icomp)%SetLabel(ng))
                StateLabel=TRIM(COUPLED(Idata)%ExpLabel(nd))
                CALL NUOPC_AddNestedState (ExportState,                 &
     &                                     CplSet=TRIM(CoupledSet),     &
     &                               nestedStateName=TRIM(StateLabel),  &
     &                               nestedState=COUPLED(Idata)%        &
     &                                           ExportState(nd,Icomp), &
     &                                     rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=__FILE__)) THEN
                  RETURN
                END IF
!
!  Add fields to export state.
!
                DO ifld=1,Nfields
                  ShortName=DataSet(Icomp)%Field(ifld)
                  id=field_index(Models(Idata)%ExportField,             &
     &                           TRIM(ShortName))
                  IF (id.gt.0) THEN
                    StandardName=MODELS(Idata)%ExportField(id)%         &
     &                                                    standard_name
                    CALL NUOPC_Advertise (COUPLED(Idata)%               &
     &                                           ExportState(nd,Icomp), &
     &                                 StandardName=TRIM(StandardName), &
     &                                    name=TRIM(ShortName),         &
     &                                    rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
                  ELSE
                    IF (localPET.eq.0) THEN
                      WRITE (dataout,10) TRIM(ShortName)
                    END IF
                    rc=ESMF_RC_NOT_FOUND
                    RETURN
                  END IF
                END DO
              END IF
            END DO
          END IF
        END IF
      END DO
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_SetInitializeP1',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
 10   FORMAT (1x,'DATA_SetInitializeP1 - unable to find field ''',a,    &
     &        ''' in ''Models(Idata)%ExportField'' list')
!
      RETURN
      END SUBROUTINE DATA_SetInitializeP1
!
      SUBROUTINE DATA_SetInitializeP2 (model,                           &
     &                                 ImportState, ExportState,        &
     &                                 clock, rc)
!
!=======================================================================
!                                                                      !
!  DATA component Phase 2 initialization: Initializes DATA structure,  !
!  sets export fields grid, and adds fields to export states.          !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mod_scalars, ONLY : NoError, exit_flag
      USE mod_strings, ONLY : Nregion, my_cpu, my_fc, my_fflags,        &
     &                        my_fort, my_os, Rdir
!
      USE strings_mod, ONLY : FoundError
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
      integer :: is, localPET, lstr, PETcount, MyComm
      integer :: ExportCount
!
      real(dp) :: TimeInDays, Time_Current
!
      character (len=20) :: Time_CurrentString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_SetInitializeP2"
!
      TYPE (ESMF_TimeInterval) :: TimeStep
      TYPE (ESMF_Time)         :: CurrentTime
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_SetInitializeP2',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the
!  mpi communicator handle and current node rank.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       vm=vm,                                     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
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
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      ESMcomm(Idata)=MyComm
!
!-----------------------------------------------------------------------
!  If concurrent PET layout, call ROMS 'allocate parallel' routine over
!  all DATA component PETs.  It is needed because the DATA component
!  uses ROMS NetCDF- and mpi-framework.  We need to allocate all the
!  ROMS profiling arrays in such PETs. Set Master, InpThread, and
!  OutThread switches for this communicator needed for processing
!  NetCDF files.  Notice that OCN_COMM_WORLD is set to the DATA
!  component communicator.
!-----------------------------------------------------------------------
!
      IF (PETlayoutOption.eq.'CONCURRENT') THEN
        IF (.not.allocated(proc)) THEN
          allocate ( proc(0:1,4,NgridsR) )
        END IF
        proc(0:1,1:4,1:NgridsR)=0
!
        IF (.not.allocated(Cstr)) THEN
          allocate ( Cstr(0:Nregion,4,NgridsR) )
        END IF
        Cstr(0:Nregion,1:4,1:NgridsR)=0.0_r8
!
        IF (.not.allocated(Cend)) THEN
          allocate ( Cend(0:Nregion,4,NgridsR) )
        END IF
        Cend(0:Nregion,1:4,1:NgridsR)=0.0_r8
!
        IF (.not.allocated(Csum)) THEN
          allocate ( Csum(0:Nregion,4,NgridsR) )
        END IF
        Csum(0:Nregion,1:4,1:NgridsR)=0.0_r8
!
        OCN_COMM_WORLD=MyComm
        MyRank=localPET
        Ctotal=0.0_r8
        total_cpu=0.0_r8
        total_model=0.0_r8
        Lwclock=.TRUE.
        CALL initialize_parallel
!
!  The standard output is redirected to an specific file for clarity and
!  it unit is redifined.
!
        dataout=101                    ! overwite Fortran default unit 6
!
        IF (localPET.eq.0) THEN
          OPEN (dataout, FILE='log.data', FORM='formatted',             &
     &          STATUS='replace')
          lstr=INDEX(my_fflags, 'free')-2
          IF (lstr.le.0) lstr=LEN_TRIM(my_fflags)
          WRITE (dataout,10) TRIM(ESMF_VERSION_STRING),                 &
     &                       TRIM(TodayDateString),                     &
     &                       TRIM(Rdir),                                &
     &                       TRIM(my_os),                               &
     &                       TRIM(my_cpu),                              &
     &                       TRIM(my_fort),                             &
     &                       TRIM(my_fc),                               &
     &                       my_fflags(1:lstr),                         &
     &                       MyComm, PETcount
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Get driver current time. Notice that the DATA component clock has
!  not been created before this initialization phase.
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet (ClockInfo(Idriver)%Clock,                     &
     &                    timeStep=TimeStep,                            &
     &                    currTime=CurrentTime,                         &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
         RETURN
      END IF

# ifdef REGRESS_STARTCLOCK
!
!  The starting time was regressed during configuration to allow the
!  proper initialization of import and export states.  We need to add
!  the coupling interval here to have the correct values for the
!  internal monotonic time coordinate (Tmono) in "DATA_inquiry" and
!  "DATA_ncread".
!
      CALL ESMF_TimeGet (CurrentTime+TimeStep,                          &
     &                   s_r8=Time_Current,                             &
     &                   timeString=Time_CurrentString)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# else
!
!  Get current time in seconds.
!
      CALL ESMF_TimeGet (CurrentTime,                                   &
     &                   s_r8=Time_Current,                             &
     &                   timeString=Time_CurrentString)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# endif
      TimeInDays=(Time_Current-                                         &
     &            ClockInfo(Idata)%Time_Reference)/86400.0_dp
      is=INDEX(Time_CurrentString, 'T')              ! remove 'T' in
      IF (is.gt.0) Time_CurrentString(is:is)=' '     ! ISO 8601 format
!
!-----------------------------------------------------------------------
!  Initilize DATA component.
!-----------------------------------------------------------------------
!
      CALL DATA_Initialize (model, TimeInDays, localPET, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Set-up grid and load coordinate data.
!-----------------------------------------------------------------------
!
      ExportCount=UBOUND(MODELS(Idata)%ExportField, DIM=1)
!
      CALL DATA_SetGridArrays (model, ExportCount, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Set DATA component land/sea mask as follows:  0: land
!                                                1: ocean
!
      MODELS(Idata)%LandValue=0
      MODELS(Idata)%SeaValue=1
!
!-----------------------------------------------------------------------
!  Set-up fields and register export state.
!-----------------------------------------------------------------------
!
      CALL DATA_SetStates (model, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_SetInitializeP2',    &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (80('-'),/,                                                &
     &        ' Earth System Models Coupling: ESMF/NUOPC Library,',     &
     &        ' Version ',a,/,31x,a,/,                                  &
     &        80('-'),                                                  &
     &        /,1x,'Repository Root  : ',a,                             &
     &        /,1x,'Operating System : ',a,                             &
     &        /,1x,'CPU Hardware     : ',a,                             &
     &        /,1x,'Compiler System  : ',a,                             &
     &        /,1x,'Compiler Command : ',a,                             &
     &        /,1x,'Compiler Flags   : ',a,                             &
     &        /,1x,'MPI communicator : ',i0,2x,'PET size = ',i0,        &
     &        /,80('-'),/)
!
      RETURN
      END SUBROUTINE DATA_SetInitializeP2
!
      SUBROUTINE DATA_Initialize (model, Tcurrent, localPET, rc)
!
!=======================================================================
!                                                                      !
!  Initializes DATA component upper level structure "DataSet", which   !
!  includes complete information about the fields to export.           !
!                                                                      !
!  The infomation is gathered as follows:                              !
!                                                                      !
!    (1) Initializes the Input Files Structure, IFS(1:Nmodels), so the !
!        associated single or multi file dataset is selected during    !
!        initiliazation or restart.                                    !
!    (2) Inquire the contents of NetCDF files associated with each     !
!        export field.                                                 !
!    (3) Reads source data longitude and latitude for each export      !
!        field and, if avalable, read land/sea mask.                   !
!    (4) Reads LOWER time-snapshot data for each export field such     !
!        that: LowerSnapshot < Tcurrent < UpperSnapshot. This is       !
!        needed for the time interpolation elsewhere.                  !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     model      DATA component object (TYPE ESMF_GridComp)            !
!     Tcurrent   Current time in days since reference date (real)      !
!     localPET   Local Persistent Execution Thread (integer)           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     DataSet    Updated DATA component structure                      !
!     rc         Return code flag (integer)                            !
!                                                                      !
!  WARNING:                                                            !
!                                                                      !
!  This routine uses ROMS NetCDF managing framework.                   !
!                                                                      !
!=======================================================================
!
      USE mod_scalars, ONLY : NoError, exit_flag
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: localPET
      integer, intent(out) :: rc
!
      real(dp) :: Tcurrent
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      logical :: FirstPass, IsUpdated, Lmulti
!
      integer :: ExportCount, NfieldsTotal
      integer :: Icomp, Nfields, Nfiles, ifld
      integer :: ifile, iMulti, nMultiFiles
      integer :: id, is
!
      character (len=20 ) :: FieldName, nc_vname, nc_tname
      character (len=100) :: vunits
      character (len=256) :: mfile, ncfile, longname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_Initialize"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_Initialize',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Initialize DATA model input files structures (IFS), so the
!  appropriate single or multi-file is selected during initialization
!  or restart.
!-----------------------------------------------------------------------
!
      DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            CALL DATA_multifile (Tcurrent, DataSet(Icomp)%IFS, Nfiles,  &
     &                           localPET, rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Inquire DATA model NetCDF files associated with each export field.
!-----------------------------------------------------------------------
!
      ExportCount=UBOUND(MODELS(Idata)%ExportField, DIM=1)
      NfieldsTotal=0                ! total number of fields to export.
!                                     It should be equal to ExportCount
      Lmulti=.FALSE.                ! multi-file switch: needs to be
!                                     false during initialization
      DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
            NfieldsTotal=NfieldsTotal+Nfields
            DO ifld=1,Nfields
              DataSet(Icomp)%Export(ifld)%Lmulti=.FALSE.    ! initilize
              FieldName=DataSet(Icomp)%Field(ifld)
              id=field_index(Models(Idata)%ExportField,                 &
     &                       TRIM(FieldName))
              IF (id.gt.0) THEN
                nc_vname=TRIM(Models(Idata)%ExportField(id)%nc_vname)
                nc_tname=TRIM(Models(Idata)%ExportField(id)%nc_tname)

                longname=TRIM(Models(Idata)%ExportField(id)%long_name)
                vunits  =TRIM(Models(Idata)%ExportField(id)%src_units)
                DataSet(Icomp)%Export(ifld)%Vdescriptor=TRIM(longname)
                DataSet(Icomp)%Export(ifld)%Vunits=TRIM(vunits)
                CALL DATA_inquiry (ifld, nc_vname, nc_tname,            &
     &                             Tcurrent,                            &
     &                             DataSet(Icomp)%Export, Nfields,      &
     &                             DataSet(Icomp)%IFS, Nfiles,          &
     &                             Lmulti, localPET, rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
              ELSE
                IF (localPET.eq.0) WRITE (dataout,10) TRIM(FieldName)
                exit_flag=5
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         MyFile)) THEN
                  rc=ESMF_RC_NOT_FOUND
                  RETURN
                END IF
              END IF
            END DO
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Read in export field longitude and latitude. If available, read
!  land/sea mask or compute it if possible.
!-----------------------------------------------------------------------
!
!  Several variables in structure "DataSet(Icomp)%Export(ifld)" will
!  be updated.
!
      DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
            DO ifld=1,Nfields
              FieldName=DataSet(Icomp)%Field(ifld)
              CALL DATA_ncvarcoords (ifld, FieldName, Nfields,          &
     &                               DataSet(Icomp)%Export,             &
     &                               localPET, rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
              DataSet(Icomp)%Export(ifld)%Icomp=Icomp
            END DO
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Read in DATA component fields from source NetCDF files and load
!  into first record of snapshot arrays. This is the first pass to
!  the processing of data needed for the time interpolation of fields
!  to export.
!-----------------------------------------------------------------------
!
      IF ((localPET.eq.0).and.(PETlayoutOption.eq.'CONCURRENT')) THEN
        WRITE (dataout,20)
      END IF
      FirstPass=.TRUE.
      CALL DATA_ncread (Tcurrent, FirstPass, localPET, IsUpdated, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Initilize the multifile switch indicating that the DATA time records
!  are split into several files, like monthly, annual. etc.  It needs
!  to be done after the first call to "DATA_ncread".  The logic for it
!  is complicated because of the various possibilities.
!-----------------------------------------------------------------------
!
      DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
            FLD_LOOP : DO ifld=1,Nfields
              ncfile=DataSet(Icomp)%Export(ifld)%ncfile
              DataSet(Icomp)%Export(ifld)%Lmulti=.FALSE.    ! initilize
              IFS_LOOP : DO ifile=1,Nfiles
                nMultiFiles=DataSet(Icomp)%IFS(ifile)%Nfiles
                IF (nMultiFiles.gt.1) THEN
                  DO iMulti=1,nMultiFiles
                    mfile=DataSet(Icomp)%IFS(ifile)%files(iMulti)
                    IF (TRIM(ncfile).eq.TRIM(mfile)) THEN
                      DataSet(Icomp)%Export(ifld)%Lmulti=.TRUE.
                      EXIT IFS_LOOP
                    END IF
                  END DO
                END IF
              END DO IFS_LOOP
            END DO FLD_LOOP
          END IF
        END IF
      END DO
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_Initialize',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (/,' DATA_Initialize - cannot find export field: ',a,      &
     &        /,19x,"in structure 'Models(Idata)%ExportField'")
  20  FORMAT (/,' DATA Component Processing Log:',/,1x,29('='),/)

      RETURN
      END SUBROUTINE DATA_Initialize

# ifdef TIME_INTERP
!
      SUBROUTINE DATA_DataInit (model, rc)
!
!=======================================================================
!                                                                      !
!  Exports DATA component fields during initialization or restart.     !
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
      logical :: FirstPass, IsUpdated, Lreport
!
      integer :: is, localPET, PETcount
!
      real(dp) :: TimeInDays, Time_Current
!
      character (len=20) :: Time_CurrentString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_DataInit"
!
      TYPE (ESMF_Clock) :: clock
      TYPE (ESMF_Time)  :: CurrentTime
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_DataInit',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get DATA component clock.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       clock=clock,                               &
     &                       localPet=localPET,                         &
     &                       petCount=PETcount,                         &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_ClockGet (clock,                                        &
     &                    currTime=CurrentTime,                         &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeGet (CurrentTime,                                   &
     &                   s_r8=Time_Current,                             &
     &                   timeString=Time_CurrentString)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      TimeInDays=(Time_Current-                                         &
     &            ClockInfo(Idata)%Time_Reference)/86400.0_dp
      is=INDEX(Time_CurrentString, 'T')              ! remove 'T' in
      IF (is.gt.0) Time_CurrentString(is:is)=' '     ! ISO 8601 format
!
!-----------------------------------------------------------------------
!  Export LOWER time-snapshot. The target ESM component will time
!  interpolates the needed field using the exported two time-level
!  data snapshots, internally.
!-----------------------------------------------------------------------
!
      Lreport=.TRUE.
      CALL DATA_Export (model, Lreport, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_DataInit',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE DATA_DataInit
# endif
!
      SUBROUTINE DATA_SetClock (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets DATA component date calendar, start and stop time, and         !
!  coupling interval.  At initilization, the variable "tdays" is       !
!  the initial time meassured in fractional days since the reference   !
!  time.                                                               !
!                                                                      !
!=======================================================================
!
      USE dateclock_mod, ONLY : caldate, time_string
      USE strings_mod,   ONLY : lowercase
!
!  Imported variable declarations.
!
      integer, intent(out):: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      integer :: Icomp, Nfields, Nfiles, ifld, ig
      integer :: localPET, PETcount
      integer :: TimeFrac
      integer :: MyStartTime(6), MyStopTime(6)
!
      real(dp) :: Tmin, Tmax, Tstr, Tend
      real(dp) :: Time_Stop
!
      character (len=22) :: Calendar

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_SetClock"
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
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_SetClock',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the
!  mpi communicator handle and current node rank.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       localPet=localPET,                         &
     &                       petCount=PETcount,                         &
     &                       vm=vm,                                     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Create DATA component clock.
!-----------------------------------------------------------------------
!
!  Create model calendar.
!
      SELECT CASE (TRIM(lowercase(ClockInfo(Idata)%CalendarString)))
        CASE ('gregorian')
          CalType=ESMF_CALKIND_GREGORIAN
          Calendar=ClockInfo(Idata)%CalendarString
        CASE ('year_360_day', '360_day')
          CalType=ESMF_CALKIND_360DAY
          Calendar=ClockInfo(Idata)%CalendarString
      END SELECT
!
      ClockInfo(Idata)%Calendar=ESMF_CalendarCreate(CalType,            &
     &                                              name=TRIM(Calendar),&
     &                                              rc=rc)
!
!  Inquire DATA component high-level structure for the minimum and
!  maximum value of available times.  It assumes that the reference
!  is the same as Time_Reference set during configuration.
!
      Tmin= MISSING_dp
      Tmax=-MISSING_dp
      Tstr= MISSING_dp
      Tend=-MISSING_dp
      DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
            DO ifld=1,Nfields
              Tmin=MIN(Tmin, DataSet(Icomp)%Export(ifld)%Tmin)
              Tmax=MAX(Tmax, DataSet(Icomp)%Export(ifld)%Tmax)
              Tstr=MIN(Tstr, DataSet(Icomp)%Export(ifld)%Tstr)
              Tend=MAX(Tend, DataSet(Icomp)%Export(ifld)%Tend)
            END DO
          END IF
        END IF
      END DO
!
!  Set DATA component starting time. Notice that is value can be less
!  than driver Time_Start since at initialization it represent the LOWER
!  time-snapshot used for time interpolation. ROMS routine "caldate"
!  adds the reference time internally.
!
      IF ((Tstr+ClockInfo(Idata)%Time_Reference/86400.0_dp).le.         &
     &    (ClockInfo(Idriver)%Time_Start/86400.0_dp)) THEN
        ClockInfo(Idata)%Time_Start=Tstr*86400.0_dp
        CALL caldate (Tstr,                                             &
     &                yy_i=MyStartTime(1),                              &
     &                mm_i=MyStartTime(2),                              &
     &                dd_i=MyStartTime(3),                              &
     &                h_i =MyStartTime(4),                              &
     &                m_i =MyStartTime(5),                              &
     &                s_i =MyStartTime(6))
        CALL time_string (ClockInfo(Idata)%Time_Start,                  &
     &                    ClockInfo(Idata)%Time_StartString)
!
        CALL ESMF_TimeSet (ClockInfo(Idata)%StartTime,                  &
     &                     yy=MyStartTime(1),                           &
     &                     mm=MyStartTime(2),                           &
     &                     dd=MyStartTime(3),                           &
     &                     h =MyStartTime(4),                           &
     &                     m =MyStartTime(5),                           &
     &                     s =MyStartTime(6),                           &
     &                     calendar=ClockInfo(Idata)%Calendar,          &
     &                     rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
!
!  Set DATA component stopping time.  Use coupling simulation stopping
!  time. ROMS routine "caldate" adds the reference time internally.
!
      ClockInfo(Idata)%Time_Stop=ClockInfo(Idriver)%Time_Stop
      Time_Stop=(ClockInfo(Idata)%Time_Stop-                            &
     &           ClockInfo(Idata)%Time_Reference)/86400.0_dp
      CALL caldate (Time_Stop,                                          &
     &              yy_i=MyStopTime(1),                                 &
     &              mm_i=MyStopTime(2),                                 &
     &              dd_i=MyStopTime(3),                                 &
     &              h_i =MyStopTime(4),                                 &
     &              m_i =MyStopTime(5),                                 &
     &              s_i =MyStopTime(6))
      CALL time_string (ClockInfo(Idata)%Time_Stop-                     &
     &                  ClockInfo(Idata)%Time_Reference,                &
     &                  ClockInfo(Idata)%Time_StopString)
!
      CALL ESMF_TimeSet (ClockInfo(Idata)%StopTime,                     &
     &                   yy=MyStopTime(1),                              &
     &                   mm=MyStopTime(2),                              &
     &                   dd=MyStopTime(3),                              &
     &                   h =MyStopTime(4),                              &
     &                   m =MyStopTime(5),                              &
     &                   s =MyStopTime(6),                              &
     &                   calendar=ClockInfo(Idata)%Calendar,            &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
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
      TimeFrac=MAX(TimeFrac,                                            &
     &             MAXVAL(MODELS(Idata)%TimeFrac(1,:),                  &
     &                    mask=MODELS(:)%IsActive))
      IF (TimeFrac.lt.1) THEN              ! needs to be 1 or greater
        rc=ESMF_RC_NOT_SET                 ! cannot be 0
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
      ClockInfo(Idata)%TimeStep=ClockInfo(Idriver)%TimeStep/TimeFrac
!
!-----------------------------------------------------------------------
!  Create ROMS component clock.
!-----------------------------------------------------------------------
!
      ClockInfo(Idata)%Name='DATA_clock'
      clock=ESMF_ClockCreate(ClockInfo(Idata)%TimeStep,                 &
     &                       ClockInfo(Idata)%StartTime,                &
     &                       stopTime =ClockInfo(Idata)%StopTime,       &
     &                       refTime  =ClockInfo(Idata)%ReferenceTime,  &
     &                       name     =TRIM(ClockInfo(Idata)%Name),     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      ClockInfo(Idata)%Clock=clock
!
!  Get current time.
!
      CALL ESMF_ClockGet (ClockInfo(Idata)%Clock,                       &
     &                    currTime=ClockInfo(Idata)%CurrentTime,        &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_SetClock',           &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE DATA_SetClock
!
      SUBROUTINE DATA_SetGridArrays (model, ExportCount, rc)
!
!=======================================================================
!                                                                      !
!  Sets DATA component horizontal grids arrays for each export fields  !
!  and land/sea mask, if any.                                          !
!                                                                      !
!=======================================================================
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ExportCount
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp), intent(inout) :: model
!
!  Local variable declarations.
!
      integer :: Icomp, Nfields, Nfiles, ifld
      integer :: Im, Istr, Iend, Jm, Jstr, Jend, i, j
      integer :: localDE, localDEcount
!
      integer (i4b), pointer :: ptrM(:,:) => NULL()     ! land/sea mask
!
      real (dp), pointer :: ptrX(:,:) => NULL()         ! longitude
      real (dp), pointer :: ptrY(:,:) => NULL()         ! latitude
!
      character (len=40) :: GridName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_SetGridArrays"
!
      TYPE (ESMF_Decomp_Flag) :: decompflag(3)
      TYPE (ESMF_DistGrid)    :: distGrid
      TYPE (ESMF_Grid)        :: grid
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_SetGridArrays',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Create ESMF DistGrid object for each DATA component export field.
!  Recall that source data may come from different grids.
!-----------------------------------------------------------------------
!
!  Set decomposition flag: divide the elements of DEs and assign the
!  rest of the division to the last DE.
!
      decompflag=(/ ESMF_DECOMP_RESTLAST,                               &
     &              ESMF_DECOMP_RESTLAST,                               &
                    ESMF_DECOMP_RESTLAST /)
!
!  Create grid object for each export field. Currently, all grids have
!  two spatial dimensions (lon,lat).
!
      FIELD_LOOP : DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
!
!  Create grid decomposition object.
!
            DO ifld=1,Nfields
              Im=DataSet(Icomp)%Export(ifld)%Vsize(1)
              Jm=DataSet(Icomp)%Export(ifld)%Vsize(2)
              distGrid=ESMF_DistGridCreate(minIndex=(/1,1/),            &
     &                                     maxIndex=(/Im,Jm/),          &
     &                                     regDecomp=(/ItileD,JtileD/), &
     &                                     decompflag=decompflag(1:2),  &
     &                                     rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  Create export field associated grid object. The array indices are
!  global.
!
              GridName=TRIM(DataSet(Icomp)%Field(ifld))//'_'//          &
     &                 TRIM(MODELS (Icomp)%name)
              grid=ESMF_GridCreate(distgrid=distgrid,                   &
     &                             indexflag=ESMF_INDEX_GLOBAL,         &
     &                             name=TRIM(GridName),                 &
     &                             rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  Get number of local decomposition elements (DEs). Usually, a single
!  DE is associated with each Persistent Execution Thread (PETs). Thus,
!  localDEcount=1.
!
              CALL ESMF_GridGet (grid,                                  &
     &                           localDECount=localDEcount,             &
     &                           rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  Allocate coordinate storage associated with staggered grid type.
!  No coordinate values are set yet.
!
              CALL ESMF_GridAddCoord (grid,                             &
     &                                staggerEdgeLWidth=(/0,0/),        &
     &                                staggerEdgeUWidth=(/0,0/),        &
     &                                rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  Allocate storage for land/sea masking.
!
              CALL ESMF_GridAddItem (grid,                              &
     &                               itemflag=ESMF_GRIDITEM_MASK,       &
     &                               rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
              DataSet(Icomp)%export(ifld)%LandValue=0
              DataSet(Icomp)%export(ifld)%SeaValue=1
!
!  Get pointers and set coordinates for the grid.  Usually, the DO-loop
!  is executed once since localDEcount=1.
!
              DE_LOOP : DO localDE=0,localDEcount-1
                CALL ESMF_GridGetCoord (grid,                           &
     &                                  localDE=localDE,                &
     &                                  coordDim=1,                     &
     &                                  farrayPtr=ptrX,                 &
     &                                  rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
                CALL ESMF_GridGetCoord (grid,                           &
     &                                  localDE=localDE,                &
     &                                  coordDim=2,                     &
     &                                  farrayPtr=ptrY,                 &
     &                                  rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
                CALL ESMF_GridGetItem (grid,                            &
     &                                 localDE=localDE,                 &
     &                                 itemflag=ESMF_GRIDITEM_MASK,     &
     &                                 farrayPtr=ptrM,                  &
     &                                 rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Fill grid pointers. In the DATA model the longitude, latitude, and
!  mask are all of the same size with identical parallel decomposition.
!
                Istr=LBOUND(ptrX,1)
                Iend=UBOUND(ptrX,1)
                Jstr=LBOUND(ptrX,2)
                Jend=UBOUND(ptrX,2)
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    ptrX(i,j)=DataSet(Icomp)%Export(ifld)%lon(i,j)
                    ptrY(i,j)=DataSet(Icomp)%Export(ifld)%lat(i,j)
                    ptrM(i,j)=INT(DataSet(Icomp)%Export(ifld)%mask(i,j))
                  END DO
                END DO
!
!  Save grid object in data strcuture.
!
                DataSet(Icomp)%export(ifld)%grid=grid
!
!  Nullify pointers.
!
                IF ( associated(ptrX) ) nullify (ptrX)
                IF ( associated(ptrY) ) nullify (ptrY)
                IF ( associated(ptrM) ) nullify (ptrM)
              END DO DE_LOOP
!!
!!  Assign grid to gridded component.
!!  (HGA: how this is done for this particular case)
!!
!!            CALL ESMF_GridCompSet (model,                             &
!!   &                               grid=grid,                         &
!!   &                               rc=rc)
!!            IF (ESMF_LogFoundError(rcToCheck=rc,                      &
!!   &                               msg=ESMF_LOGERR_PASSTHRU,          &
!!   &                               line=__LINE__,                     &
!!   &                               file=MyFile)) THEN
!!              RETURN
!!            END IF
            END DO
          END IF
        END IF
      END DO FIELD_LOOP
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_SetGridArrays',      &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE DATA_SetGridArrays
!
      SUBROUTINE DATA_SetStates (model, rc)
!
!=======================================================================
!                                                                      !
!  Adds DATA component fields into export state.                       !
!                                                                      !
!=======================================================================
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp), intent(inout) :: model
!
!  Local variable declarations.
!
      integer :: Icomp, Nfields, Nfiles, ifld, nd, ng
      integer :: localDE, localDEcount, localPET
      integer :: ExportCount
!
      real (dp), dimension(:,:), pointer :: ptr2d => NULL()
!
      character (len=10) :: AttList(1)
      character (len=20) :: FieldName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_SetStates"
!
      character (ESMF_MAXSTR), allocatable :: ExportNameList(:)
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
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_SetStates',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Query gridded component.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (model,                                     &
     &                       localPet=localPET,                         &
     &                       vm=vm,                                     &
     &                       rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
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
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      staggerLoc=ESMF_STAGGERLOC_CENTER
!
!-----------------------------------------------------------------------
!  Add export fields into export state.
!-----------------------------------------------------------------------
!
!  Set export field(s).
!
      FIELD_LOOP : DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
            DO ng=1,MODELS(Icomp)%Ngrids
              IF (COUPLED(Icomp)%LinkedGrid(ng,Idata)) THEN
                nd=COUPLED(Idata)%DataCoupledSets(ng,Icomp)
!
!  For debugging, inquire state about the number of fields to export.
!  It should be the same as Nfields.
!
                CALL ESMF_StateGet (COUPLED(Idata)%                     &
     &                                     ExportState(nd,Icomp),       &
     &                             itemCount=ExportCount,               &
     &                             rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  For debugging, inquire state about the list of export field names. It
!  should be the same as the ones advertised in "DATA_SetInitializeP1".
!
                IF (.not.allocated(ExportNameList)) THEN
                  allocate ( ExportNameList(ExportCount) )
                END IF
                CALL ESMF_StateGet (COUPLED(Idata)%                     &
     &                                     ExportState(nd,Icomp),       &
     &                              itemNameList=ExportNameList,        &
     &                              rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
!
!  Add export fields to the state set.
!
                DO ifld=1,Nfields
                  FieldName=DataSet(Icomp)%Field(ifld)
!
                  IF (NUOPC_IsConnected(COUPLED(Idata)%                 &
     &                                         ExportState(nd,Icomp),   &
     &                                  fieldName=TRIM(FieldName),      &
     &                                  rc=rc)) THEN
!
!  Create 2D field from the Grid and arraySpec.
!
                    field=ESMF_FieldCreate(DataSet(Icomp)%              &
     &                                            export(ifld)%grid,    &
     &                                     arraySpec,                   &
     &                                     indexflag=ESMF_INDEX_GLOBAL, &
     &                                     staggerloc=staggerLoc,       &
     &                                     name=TRIM(FieldName),        &
     &                                     rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
!
!  Get number of local decomposition elements (DEs). Usually, a single
!  Decomposition Element (DE) is associated with each Persistent
!  Execution Thread (PETs). Thus, localDEcount=1.
!
                    CALL ESMF_GridGet (DataSet(Icomp)%                  &
     &                                        export(ifld)%grid,        &
     &                                 localDECount=localDEcount,       &
     &                                 rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF

# ifdef TIME_INTERP_NOT_WORKING
!
!  Create standard Attribute Package for each export field. Then, nest
!  custom Attribute Package around it.
!
                    CALL ESMF_AttributeAdd (field,                      &
     &                                      convention='ESMF',          &
     &                                      purpose='General',          &
     &                                      rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
!
                    AttList(1)='TimeInterp'
                    CALL ESMF_AttributeAdd (field,                      &
     &                                  convention='CustomConvention',  &
     &                                      purpose='General',          &
!!   &                                      purpose='Instance',         &
     &                                      attrList=AttList,           &
     &                                      nestConvention='ESMF',      &
     &                                      nestPurpose='General',      &
     &                                      rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
!
                    CALL ESMF_AttributeLink (ExportState, field, rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
# endif
!
!  Get pointer to DE-local memory allocation within field. Usually, the
!  DO-loop is executed once since localDEcount=1.
!
                    DO localDE=0,localDEcount-1
                      CALL ESMF_FieldGet (field,                        &
     &                                    localDe=localDE,              &
     &                                    farrayPtr=ptr2d,              &
     &                                    rc=rc)
                      IF (ESMF_LogFoundError(rcToCheck=rc,              &
     &                                       msg=ESMF_LOGERR_PASSTHRU,  &
     &                                       line=__LINE__,             &
     &                                       file=MyFile)) THEN
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
                    CALL NUOPC_Realize (COUPLED(Idata)%                 &
     &                                         ExportState(nd,Icomp),   &
     &                                  field=field,                    &
     &                                  rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
!
!  Save field pointer. It will be used by the coupler for regridding. It
!  needs to be saved only for the coarser grid since the other nested
!  grids uses the same dataset. There is a problem with the FieldBundles
!  in the DATA component that yields empty field objects. Perhaps it is
!  because the FieldBundle stores similar fields discretized on the same
!  Grid. The DATA component is generic, and a different grid object is
!  defined for each exported field.
!
                    IF (ng.eq.1) THEN
                       DataSet(Icomp)%export(ifld)%field=field
                    END IF
!
!  Remove field from export state because it is not connected.
!
                  ELSE
                    IF (localPET.eq.0) THEN
                      WRITE (dataout,10) TRIM(FieldName),               &
     &                                   'Export State: ',              &
     &                                   TRIM(COUPLED(Idata)%           &
     &                                               ExpLabel(nd))
                    END IF
                    CALL ESMF_StateRemove (COUPLED(Idata)%              &
     &                                            ExportState(nd,Icomp),&
     &                                     (/ TRIM(FieldName) /),       &
     &                                     rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
                  END IF
                END DO
              END IF
            END DO
!
!  Deallocate temporary variable.
!
            IF ( allocated(ExportNameList) ) THEN
              deallocate (ExportNameList)
            END IF
          END IF
        END IF
      END DO FIELD_LOOP
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_SetStates',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
 10   FORMAT (1x,'DATA_SetStates - Removing field ''',a,''' from ',a,   &
     &        '''',a,'''',/,18x,'because it is not connected.')
!
      RETURN
      END SUBROUTINE DATA_SetStates
!
      SUBROUTINE DATA_ModelAdvance (model, rc)
!
!=======================================================================
!                                                                      !
!  Advance DATA component for a coupling interval (seconds) using      !
!  "DATA_run". It also calls "DATA_Export" to export coupling fields.  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc

      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      logical :: FirstPass, IsUpdated, Lreport
!
      integer :: MyTask, PETcount, is, localPET, phase
!
      real (dp) :: CouplingInterval, RunInterval
      real (dp) :: TcurrentInSeconds, TstopInSeconds
      real (dp) :: TcurrentInDays
!
      character (len=22) :: Cinterval
      character (len=22) :: CurrTimeString, StopTimeString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_ModelAdvance"
!
      TYPE (ESMF_Clock)        :: clock
      TYPE (ESMF_State)        :: ExportState, ImportState
      TYPE (ESMF_Time)         :: ReferenceTime
      TYPE (ESMF_Time)         :: StartTime, StopTime
      TYPE (ESMF_TimeInterval) :: TimeFrom, TimeTo, TimeStep
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_ModelAdvance',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get information about the gridded component.
!-----------------------------------------------------------------------
!
!  Inquire about DATA component.
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
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
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
     &                    currTime=ClockInfo(Idata)%CurrentTime,        &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!  Current DATA component time (seconds).
!
      CALL ESMF_TimeGet (ClockInfo(Idata)%CurrentTime,                  &
     &                   s_r8=TcurrentInSeconds,                        &
     &                   timeStringISOFrac=CurrTimeString,              &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(CurrTimeString, 'T')                 ! remove 'T' in
      IF (is.gt.0) CurrTimeString(is:is)=' '        ! ISO 8601 format
!
!  DATA component stop time (seconds) for this coupling window.
!
      CALL ESMF_TimeGet (ClockInfo(Idata)%CurrentTime+TimeStep,         &
     &                   s_r8=TstopInSeconds,                           &
     &                   timeStringISOFrac=StopTimeString,              &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(StopTimeString, 'T')                 ! remove 'T' in
      IF (is.gt.0) StopTimeString(is:is)=' '        ! ISO 8601 form
!
!  Get coupling interval (seconds, double precision).
!
      CALL ESMF_TimeIntervalGet (TimeStep,                              &
     &                           s_r8=CouplingInterval,                 &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
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
        WRITE (dataout,10) TRIM(CurrTimeString), TRIM(StopTimeString),  &
     &                    phase, TRIM(ADJUSTL(Cinterval))
      END IF
!
!-----------------------------------------------------------------------
!  Run DATA component.
!-----------------------------------------------------------------------

# ifdef REGRESS_STARTCLOCK
!
!  If applicable, read in the next data time-snapshot. Recall that the
!  current time was adjusted during configuration to allow the exporting
!  of both time snapshots. We need to add here the coupling interval to
!  have the correct values for the internal monotonic time coordinate
!  (Tmono) in "DATA_ncread" with day units since reference date.
!
      FirstPass=.FALSE.
      TcurrentInDays=(TcurrentInSeconds+CouplingInterval-               &
     &                ClockInfo(Idata)%Time_Reference)/86400.0_dp
# else
!
!  If applicable, read in the next data time-snapshot.
!
      FirstPass=.FALSE.
      TcurrentInDays=(TcurrentInSeconds-                                &
     &                ClockInfo(Idata)%Time_Reference)/86400.0_dp
# endif
      CALL DATA_ncread (TcurrentInDays, FirstPass, localPET, IsUpdated, &
     &                  rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Send export fields to destination ESM component.
!-----------------------------------------------------------------------
# ifdef TIME_INTERP
!
!  Since the DATA component is exporting data-snaphots of the fields
!  for time interpolation in the destination ESM component, it needs
!  to export when reading new data (IsUpdate=.TRUE.).
!
      IF (IsUpdated) THEN
        Lreport=.TRUE.
        CALL DATA_Export (model, Lreport, rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
      END IF
# else
!
!  The Data component always export fields at every coupling timestep
!  since "DATA_Export" time interpolates the fields from available
!  LOWER and UPPER snapshots.
!
      IF (TcurrentInSeconds.eq.ClockInfo(Idriver)%Time_Start) THEN
        Lreport=.FALSE.           ! export state pointers will give
      ELSE                        ! MISSING_r8 values, avoid report
        Lreport=.TRUE.
      END IF
      CALL DATA_Export (model, Lreport, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# endif
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_ModelAdvance',       &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (3x,'ModelAdvance - ESMF, Running DATA:',t42,a,            &
     &        ' => ',a,', Phase: ',i1,' [',a,' s]')

      RETURN
      END SUBROUTINE DATA_ModelAdvance
!
      SUBROUTINE DATA_SetFinalize (model,                               &
     &                             ImportState, ExportState,            &
     &                             clock, rc)
!
!=======================================================================
!                                                                      !
!  Finalize DATA component execution. It calls DATA_finalize.          !
!                                                                      !
!=======================================================================
!
      USE mod_scalars,  ONLY : NoError, exit_flag
!
      USE mod_netcdf,   ONLY : netcdf_close
      USE strings_mod,  ONLY : FoundError
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
      integer :: Icomp, Nfiles, ifile, ncid
      integer :: ROMScomm
!
      integer, parameter :: imodel = 1    ! for compatibility with ROMS
      integer, parameter :: ng = 1        ! used routines
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_SetFinalize"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_SetFinalize',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Finalize DATA component. Close all input NetCDF files.
!-----------------------------------------------------------------------
!
      DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            DO ifile=1,Nfiles
              ncid=DataSet(Icomp)%IFS(ifile)%ncid
              IF (ncid.ne.-1) THEN
                CALL netcdf_close (ng, imodel, ncid)
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         MyFile)) THEN
                  rc=ESMF_RC_FILE_CLOSE
                  RETURN
                END IF
              END IF
            END DO
          END IF
        END IF
      END DO
!
      CALL my_flush (dataout)             ! flush standard output buffer
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_SetFinalize',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
      RETURN
      END SUBROUTINE DATA_SetFinalize
!
      SUBROUTINE DATA_Export (model, Lreport, rc)
!
!=======================================================================
!                                                                      !
!  Exports ROMS fields to other coupled gridded components.            !
!                                                                      !
!  Reporting is supressed for the first pass (Lreport=.FALSE) since    !
!  the clocks were regressed by one coupling interval to facilitate    !
!  initilazation phases.                                               !
!                                                                      !
!=======================================================================
!
# ifdef TIME_INTERP
      USE mod_iounits,  ONLY : SourceFile
      USE mod_netcdf,   ONLY : netcdf_put_fvar, netcdf_put_ivar
      USE mod_netcdf,   ONLY : netcdf_put_svar
      USE mod_scalars,  ONLY : NoError, exit_flag
      USE strings_mod,  ONLY : FoundError
# endif
!
!  Imported variable declarations.
!
      logical, intent(in)  :: Lreport
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      logical :: IsPresent
!
      integer :: ExportCount, localPET, PETcount
      integer :: localDE, localDEcount
# ifdef TIME_INTERP
      integer :: ROMScomm
# endif
      integer :: Icomp, Nfields, Nfiles, Nvdim, ifld, is, nd, ng
      integer :: Istr, Iend, Jstr, Jend, Kstr, Kend, i, j, k
      integer :: Tindex, id
      integer :: CurrDate(9)
      integer :: MyDateVec(9)
# ifdef TIME_INTERP
!
      integer, save :: record = 0
!
      integer, parameter :: iNLM = 1            ! ROMS framework usage
# endif
!
      real (dp) :: TimeInDays, Time_Current
      real (dp) :: Tintrp, Tmin, Tmax, Tstr, Tend, Vtime
      real (dp) :: Fmin(1), Fmax(1), MyFmin(1), MyFmax(1), Fval
      real (dp) :: MyAttValues(14)
!
      real (dp), pointer :: ptr2d(:,:) => NULL()
      real (dp), pointer :: ptr3d(:,:,:) => NULL()
!
      character (len=20)      :: ShortName
# ifdef TIME_INTERP
      character (len=20)      :: MyShortName(1,1)
      character (len=22)      :: MyDateString(1,1,1)
# endif
      character (len=22)      :: MyDate, Time_CurrentString
      character (len=40)      :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_Export"

      character (ESMF_MAXSTR) :: cname, ofile
!
      TYPE (ESMF_AttPack)      :: AttPack
      TYPE (ESMF_Clock)        :: clock
      TYPE (ESMF_Field)        :: field
      TYPE (ESMF_Time)         :: CurrentTime
      TYPE (ESMF_TimeInterval) :: TimeStep
      TYPE (ESMF_VM)           :: vm
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_Export',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
# ifdef TIME_INTERP
      SourceFile=MyFile
# endif
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
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
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
     &                    timeStep=TimeStep,                            &
     &                    currTime=CurrentTime,                         &
     &                    rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CurrDate(1:9)=0                                  ! initialize
!
      CALL ESMF_TimeGet (CurrentTime,                                   &
     &                   yy=CurrDate(1),                                &
     &                   mm=CurrDate(2),                                &
     &                   dd=CurrDate(3),                                &
     &                   h =CurrDate(4),                                &
     &                   m =CurrDate(5),                                &
     &                   s =CurrDate(6),                                &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeGet (CurrentTime,                                   &
     &                   s_r8=Time_Current,                             &
     &                   timeString=Time_CurrentString,                 &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF

# ifdef REGRESS_STARTCLOCK
!
!  Compute current time in days needed for time interpolation. Add one
!  coupling timestep since it was substracted at configuration. It is
!  needed to keep the interpolation bounded between LOWER and UPPER
!  snapshots.
!
      TimeInDays=(Time_Current+                                         &
     &            ClockInfo(Idata)%Time_Step-                           &
     &            ClockInfo(Idata)%Time_Reference)/86400.0_dp
# else
!
!  Compute current time in days needed for time interpolation.
!
      TimeInDays=(Time_Current-                                         &
     &            ClockInfo(Idata)%Time_Reference)/86400.0_dp
# endif
      is=INDEX(Time_CurrentString, 'T')              ! remove 'T' in
      IF (is.gt.0) Time_CurrentString(is:is)=' '     ! ISO 8601 format

# ifndef TIME_INTERP
!
!-----------------------------------------------------------------------
!  If no exporting the snapshots source data, time interpolate field to
!  export from snapshots records read from source NetCDF files. The
!  receiving ESM component will time interpolate the field internally.
!-----------------------------------------------------------------------
!
      CALL DATA_TimeInterp (TimeInDays, localPET, rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# endif
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
!  Load export fields.
!-----------------------------------------------------------------------
!
      FIELD_LOOP : DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
            DO ng=1,MODELS(Icomp)%Ngrids
              IF (COUPLED(Icomp)%LinkedGrid(ng,Idata)) THEN
                nd=COUPLED(Idata)%DataCoupledSets(ng,Icomp)
                DO ifld=1,Nfields
                  Nvdim=DataSet(Icomp)%Export(ifld)%Nvdim
                  ShortName=DataSet(Icomp)%Field(ifld)
!
!  Get target componend import field ID.
!
                  id=field_index(MODELS(Icomp)%ImportField, ShortName)
!
!  Get field from export state.  Use field name from DataSet structure
!  to follow the order of export data in storage.
!
                  CALL ESMF_StateGet (COUPLED(Idata)%                   &
     &                                       ExportState(nd,Icomp),     &
     &                                TRIM(ShortName),                  &
     &                                field,                            &
     &                                rc=rc)
                  IF (ESMF_LogFoundError(rcToCheck=rc,                  &
     &                                   msg=ESMF_LOGERR_PASSTHRU,      &
     &                                   line=__LINE__,                 &
     &                                   file=MyFile)) THEN
                    RETURN
                  END IF
!
!  Get number of local decomposition elements (DEs). Usually, a single
!  DE is associated with each Persistent Execution Thread (PETs). Thus,
!  localDEcount=1.
!
                  CALL ESMF_GridGet (DataSet(Icomp)%Export(ifld)%grid,  &
     &                               localDECount=localDEcount,         &
     &                               rc=rc)
                  IF (ESMF_LogFoundError(rcToCheck=rc,                  &
     &                                   msg=ESMF_LOGERR_PASSTHRU,      &
     &                                   line=__LINE__,                 &
     &                                   file=MyFile)) THEN
                    RETURN
                  END IF
!
!  Get field pointer.  Usually, the DO-loop is executed once since
!  localDEcount=1.
!
                  DE_LOOP : DO localDE=0,localDEcount-1
                    IF (Nvdim.eq.2) THEN                  ! 2D field
                      CALL ESMF_FieldGet (field,                        &
     &                                    localDE=localDE,              &
     &                                    farrayPtr=ptr2d,              &
     &                                    rc=rc)
                      IF (ESMF_LogFoundError(rcToCheck=rc,              &
     &                                       msg=ESMF_LOGERR_PASSTHRU,  &
     &                                       line=__LINE__,             &
     &                                       file=MyFile)) THEN
                        RETURN
                      END IF
                      ptr2d=MISSING_dp
                    ELSE IF (Nvdim.eq.3) THEN             ! 3D field
                      CALL ESMF_FieldGet (field,                        &
     &                                    localDE=localDE,              &
     &                                    farrayPtr=ptr3d,              &
     &                                    rc=rc)
                      IF (ESMF_LogFoundError(rcToCheck=rc,              &
     &                                       msg=ESMF_LOGERR_PASSTHRU,  &
     &                                       line=__LINE__,             &
     &                                      file=MyFile)) THEN
                        RETURN
                      END IF
                      ptr3d=MISSING_dp
                    END IF
!
!  Load field data into export state.  Notice that all export fields
!  are kept as computed by ROMS. The imported component does the
!  proper scaling, physical units conversion, and other manipulations.
!  It is done to avoid applying such transformations twice.
!
                    IF (Nvdim.eq.2) THEN                  ! 2D field
                      Istr=LBOUND(ptr2d,1)
                      Iend=UBOUND(ptr2d,1)
                      Jstr=LBOUND(ptr2d,2)
                      Jend=UBOUND(ptr2d,2)
                      Tindex=DataSet(Icomp)%Export(ifld)%Tindex
# ifdef TIME_INTERP
                      Fval=DataSet(Icomp)%Export(ifld)%                 &
     &                                          A2dG(Istr,Jstr,Tindex)
# else
                      Fval=DataSet(Icomp)%Export(ifld)%                 &
     &                                          A2d(Istr,Jstr)
# endif
                      MyFmin(1)=Fval
                      MyFmax(1)=Fval
                      DO j=Jstr,Jend
                        DO i=Istr,Iend
# ifdef TIME_INTERP
                          Fval=DataSet(Icomp)%Export(ifld)%             &
     &                                              A2dG(i,j,Tindex)
# else
                          Fval=DataSet(Icomp)%Export(ifld)%             &
     &                                              A2d(i,j)
# endif
                          MyFmin(1)=MIN(MyFmin(1),Fval)
                          MyFmax(1)=MAX(MyFmax(1),Fval)
                          ptr2d(i,j)=Fval
                        END DO
                      END DO
                      IF (associated(ptr2d)) nullify (ptr2d)
                    ELSE IF (Nvdim.eq.3) THEN             ! 3D field
                      Istr=LBOUND(ptr3d,1)
                      Iend=UBOUND(ptr3d,1)
                      Jstr=LBOUND(ptr3d,2)
                      Jend=UBOUND(ptr3d,2)
                      Kstr=LBOUND(ptr3d,3)
                      Kend=UBOUND(ptr3d,3)
                      Tindex=DataSet(Icomp)%Export(ifld)%Tindex
# ifdef TIME_INTERP
                      Fval=DataSet(Icomp)%Export(ifld)%                 &
     &                                      A3dG(Istr,Jstr,Kstr,Tindex)
# else
                      Fval=DataSet(Icomp)%Export(ifld)%                 &
     &                                      A3d(Istr,Jstr,Kstr)
# endif
                      MyFmin(1)=Fval
                      MyFmax(1)=Fval
                      DO k=Kstr,Kend
                        DO j=Jstr,Jend
                          DO i=Istr,Iend
# ifdef TIME_INTERP
                            Fval=DataSet(Icomp)%Export(ifld)%           &
     &                                            A3dG(i,j,k,Tindex)
# else
                            Fval=DataSet(Icomp)%Export(ifld)%           &
     &                                            A3d(i,j,k)
# endif
                            MyFmin(1)=MIN(MyFmin(1),Fval)
                            MyFmax(1)=MAX(MyFmax(1),Fval)
                            ptr3d(i,j,k)=Fval
                          END DO
                        END DO
                      END DO
                      IF (associated(ptr3d)) nullify (ptr3d)
                    END IF
                  END DO DE_LOOP
!
!  Get export field minimun and maximum values.
!
                  CALL ESMF_VMAllReduce (vm,                            &
     &                                   sendData=MyFmin,               &
     &                                   recvData=Fmin,                 &
     &                                   count=1,                       &
     &                                   reduceflag=ESMF_REDUCE_MIN,    &
     &                                   rc=rc)
                  IF (ESMF_LogFoundError(rcToCheck=rc,                  &
     &                                   msg=ESMF_LOGERR_PASSTHRU,      &
     &                                   line=__LINE__,                 &
     &                                   file=MyFile)) THEN
                    RETURN
                  END IF
!
                  CALL ESMF_VMAllReduce (vm,                            &
     &                                   sendData=MyFmax,               &
     &                                   recvData=Fmax,                 &
     &                                   count=1,                       &
     &                                   reduceflag=ESMF_REDUCE_MAX,    &
     &                                   rc=rc)
                  IF (ESMF_LogFoundError(rcToCheck=rc,                  &
     &                                   msg=ESMF_LOGERR_PASSTHRU,      &
     &                                   line=__LINE__,                 &
     &                                   file=MyFile)) THEN
                    RETURN
                  END IF

# ifdef TIME_INTERP
!
!  Set parameters in the destination component that need for the time
!  interpolation between snapshots and associated information.
!
                  Tmin  =DataSet(Icomp)%Export(ifld)%Tmin
                  Tmax  =DataSet(Icomp)%Export(ifld)%Tmax
                  Tstr  =DataSet(Icomp)%Export(ifld)%Tstr
                  Tend  =DataSet(Icomp)%Export(ifld)%Tend
                  Tintrp=DataSet(Icomp)%Export(ifld)%Tintrp(Tindex)
                  Vtime =DataSet(Icomp)%Export(ifld)%Vtime(Tindex)
                  MyDate=DataSet(Icomp)%Export(ifld)%DateString(Tindex)
!
                  MODELS(Icomp)%ImportField(id)%Tmin=Tmin
                  MODELS(Icomp)%ImportField(id)%Tmax=Tmax
                  MODELS(Icomp)%ImportField(id)%Tstr=Tstr
                  MODELS(Icomp)%ImportField(id)%Tend=Tend
                  MODELS(Icomp)%ImportField(id)%Tindex=Tindex
                  MODELS(Icomp)%ImportField(id)%Tintrp(Tindex)=Tintrp
                  MODELS(Icomp)%ImportField(id)%Vtime(Tindex)=Vtime
                  MODELS(Icomp)%ImportField(id)%DateString(Tindex)=MyDate
!
!  If concurrent coupling and exporting time snapshots, write time
!  interpolation metadata into a NetCDF file. It is very tricky
!  to perform inter VM communications. It is easier to write them into
!  a NetCDF file.  The importing ESM component needs these vaiables to
!  perform the time interpolation between snapshots in its kernel.
!
                  IF (PETlayoutOption.eq.'CONCURRENT') THEN
                    MyShortName(1,1)=TRIM(ShortName)
                    CALL netcdf_put_svar (ng, iNLM, AttFileName,        &
     &                                    'field',                      &
     &                                    MyShortName,                  &
     &                                    (/1,Icomp,id/),               &
     &                                    (/20,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
!
                    CALL netcdf_put_ivar (ng, iNLM, AttFileName,        &
     &                                    'Tindex',                     &
     &                                    Tindex,                       &
     &                                    (/Icomp,id,record/),          &
     &                                    (/1,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
!
                    MyDateString(1,1,1)=MyDate
                    CALL netcdf_put_svar (ng, iNLM, AttFileName,        &
     &                                    'Date',                       &
     &                                    MyDateString,                 &
     &                                    (/1,Icomp,id,record/),        &
     &                                    (/22,1,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
!
                    CALL netcdf_put_fvar (ng, iNLM, AttFileName,        &
     &                                    'Tcurrent',                   &
     &                                    TimeInDays,                   &
     &                                    (/Icomp,id,record/),          &
     &                                    (/1,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
!
                    CALL netcdf_put_fvar (ng, iNLM, AttFileName,        &
     &                                    'Tstr',                       &
     &                                    Tstr,                         &
     &                                    (/Icomp,id,record/),          &
     &                                    (/1,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
!
                    CALL netcdf_put_fvar (ng, iNLM, AttFileName,        &
     &                                    'Tend',                       &
     &                                    Tend,                         &
     &                                    (/Icomp,id,record/),          &
     &                                    (/1,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
!
                    CALL netcdf_put_fvar (ng, iNLM, AttFileName,        &
     &                                    'Tintrp',                     &
     &                                    Tintrp,                       &
     &                                    (/Icomp,id,record/),          &
     &                                    (/1,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
!
                    CALL netcdf_put_fvar (ng, iNLM, AttFileName,        &
     &                                    'Vtime',                      &
     &                                    Vtime,                        &
     &                                    (/Icomp,id,record/),          &
     &                                    (/1,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
!
                    CALL netcdf_put_fvar (ng, iNLM, AttFileName,        &
     &                                    'Tmin',                       &
     &                                    Tmin,                         &
     &                                    (/Icomp,id,record/),          &
     &                                    (/1,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
!
                    CALL netcdf_put_fvar (ng, iNLM, AttFileName,        &
     &                                    'Tmax',                       &
     &                                    Tmax,                         &
     &                                    (/Icomp,id,record/),          &
     &                                    (/1,1,1/))
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      rc=ESMF_RC_FILE_WRITE
                      RETURN
                    END IF
                  END IF
# endif
!
!  Overwrite the "TimeStamp" attribute with the correct date time
!  value.
!
# ifdef TIME_INTERP
                  MyDateVec(1:9)=0
                  DO i=1,6
                    MyDateVec(i)=INT(DataSet(Icomp)%Export(ifld)%       &
     &                                              Date(i,Tindex))
                  END DO
# endif
                  CALL ESMF_AttributeSet (field,                        &
     &                                    name='TimeStamp',             &
# ifdef TIME_INTERP
     &                                    valueList=MyDateVec,          &
# else
     &                                    valueList=CurrDate,           &
# endif
     &                                    convention='NUOPC',           &
     &                                    purpose='Instance',           &
     &                                    rc=rc)
                  IF (ESMF_LogFoundError(rcToCheck=rc,                  &
     &                                   msg=ESMF_LOGERR_PASSTHRU,      &
     &                                   line=__LINE__,                 &
     &                                   file=MyFile)) THEN
                    RETURN
                  END IF

# ifdef TIME_INTERP_NOT_WORKING
!
!  Set field attributes for time interpolation
!
                  MyAttValues( 1)=REAL(MyDateVec(1), dp)
                  MyAttValues( 2)=REAL(MyDateVec(2), dp)
                  MyAttValues( 3)=REAL(MyDateVec(3), dp)
                  MyAttValues( 4)=REAL(MyDateVec(4), dp)
                  MyAttValues( 5)=REAL(MyDateVec(5), dp)
                  MyAttValues( 6)=REAL(MyDateVec(6), dp)
                  MyAttValues( 7)=REAL(Tindex, dp)
                  MyAttValues( 8)=Tstr
                  MyAttValues( 9)=TimeInDays
                  MyAttValues(10)=Tend
                  MyAttValues(11)=Tintrp
                  MyAttValues(12)=Vtime
                  MyAttValues(13)=Tmin
                  MyAttValues(14)=Tmax
!
!  Retrieve custom Attribute Package.
!
                  CALL ESMF_AttributeGetAttPack (field,                 &
     &                                           'CustomConvention',    &
     &                                           'General',             &
!!   &                                           'Instance',            &
     &                                           attpack=AttPack,       &
     &                                           isPresent=IsPresent,   &
     &                                           rc=rc)
                  IF (ESMF_LogFoundError(rcToCheck=rc,                  &
     &                                   msg=ESMF_LOGERR_PASSTHRU,      &
     &                                   line=__LINE__,                 &
     &                                   file=MyFile)) THEN
                    RETURN
                  END IF
!
!  Set "TimeInterp" attribute for export field.
!
                  CALL ESMF_AttributeSet (field,                        &
     &                                    name='TimeInterp',            &
     &                                    valueList=MyAttValues,        &
     &                                    attpack=AttPack,              &
     &                                    rc=rc)
                  IF (ESMF_LogFoundError(rcToCheck=rc,                  &
     &                                   msg=ESMF_LOGERR_PASSTHRU,      &
     &                                   line=__LINE__,                 &
     &                                   file=MyFile)) THEN
                    RETURN
                  END IF
# endif
!
                  IF (Lreport.and.(localPET.eq.0)) THEN
                    WRITE (dataout,20) TRIM(ShortName),                 &
     &                                 TRIM(DataSet(Icomp)%             &
     &                                             Ctarget(ifld)), ng,  &
# ifdef TIME_INTERP
     &                                 TRIM(MyDate),                    &
     &                                 Fmin(1), Fmax(1), Tindex
# else
     &                                 TRIM(Time_CurrentString),        &
     &                                 Fmin(1), Fmax(1)
# endif
                  END IF
!
!  Debugging: write out field information.
!
                  IF (Lreport.and.(DebugLevel.ge.4)) THEN
                    CALL ESMF_FieldPrint(field,                         &
     &                                   rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
                  END IF
!
!  Debugging: write out field data into a NetCDF file.
!
                  IF ((DebugLevel.ge.3).and.                            &
     &                MODELS(Idata)%ExportField(ifld)%debug_write) THEN
                    WRITE (ofile,30) 'data_export', TRIM(ShortName),    &
# ifdef TIME_INTERP
     &                                MyDateVec(1:6)
# else
     &                                CurrDate(1:6)
# endif
                    CALL ESMF_FieldWrite (field,                        &
     &                                    TRIM(ofile),                  &
     &                                    variableName=TRIM(ShortName), &
     &                                    overwrite = .TRUE.,           &
     &                                    rc=rc)
                    IF (ESMF_LogFoundError(rcToCheck=rc,                &
     &                                     msg=ESMF_LOGERR_PASSTHRU,    &
     &                                     line=__LINE__,               &
     &                                     file=MyFile)) THEN
                      RETURN
                    END IF
                  END IF
                END DO
!
!  Update export sets counter.
!
                IF (Nfields.gt.0) THEN
                  MODELS(Idata)%ExportCalls=MODELS(Idata)%ExportCalls+1
                END IF
              END IF
            END DO
          END IF
        END IF
      END DO FIELD_LOOP
!
!  Flux DATA component standard out unit.
!
      CALL my_flush (dataout)
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_Export',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (i4.4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',a)
  20  FORMAT (3x,' DATA_Export - ESMF: exporting ''',a,'''',            &
     &        t50,'-> ''',a,''' Grid ',i2.2,',',t72,a,/ 19x,            &
# ifdef TIME_INTERP
     &        '(Dmin= ', 1p,e15.8,0p,' Dmax = ',1p,e15.8,0p,            &
     &        ' SnapshotIndex = ',i1,')')
# else
     &        '(DMin= ', 1p,e15.8,0p,' Dmax= ',1p,e15.8,0p,')')
# endif
  30  FORMAT (a,'_',a,'_',i4.4,2('-',i2.2),'_',i2.2,2('.',i2.2),'.nc')

      RETURN
      END SUBROUTINE DATA_Export
!
      SUBROUTINE Data_TimeInterp (Tcurrent, localPET, rc)
!
!=======================================================================
!                                                                      !
!  This routine time interpolates the fields that the DATA component   !
!  exports to other ESM components. The data is loaded to the arrays   !
!  in the structure:                                                   !
!                                                                      !
!    DataSet(Icomp)%Export(ifld)%A2d(:,:)                   2D field   !
!    DataSet(Icomp)%Export(ifld)%A3d(:,:,:)                 3D field   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Tcurrent   Current time in days (real)                           !
!     localPET   Local Persistent Execution Thread (integer)           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     DataSet    Updated DATA component structure.                     !
!     rc         Return code flag (integer)                            !
!                                                                      !
!=======================================================================
!
      USE mod_scalars, ONLY : NoError, exit_flag
!
      USE strings_mod,   ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: localPET
      integer, intent(out) :: rc
!
      real(dp), intent(in) :: Tcurrent
!
!  Local variable declarations.
!
      integer :: Icomp, Nfields, Nfiles, Nvdim, ifld
      integer :: Tindex, it1, it2
      integer :: Im, Jm, Km, i, j, k
!
      real(dp) :: Tstr, Tend, Tmin, Tmax
      real(dp) :: Tintrp(2)
      real(dp) :: DayScale, fac, fac1, fac2, w1, w2
      real(dp) :: Fval
!
      character (len=100) :: Vname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_TimeInterp"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_TimeInterp',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Time interpolate fields to export from data snapshots.
!-----------------------------------------------------------------------
!
      FIELD_LOOP : DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
            DO ifld=1,Nfields
              Nvdim=DataSet(Icomp)%Export(ifld)%Nvdim
!
!  Load properties of export field to read from source NetCDF file.
!
              Vname    =DataSet(Icomp)%Export(ifld)%Vname
              Tmin     =DataSet(Icomp)%Export(ifld)%Tmin
              Tmax     =DataSet(Icomp)%Export(ifld)%Tmax
              Tstr     =DataSet(Icomp)%Export(ifld)%Tstr
              Tend     =DataSet(Icomp)%Export(ifld)%Tend
              Tindex   =DataSet(Icomp)%Export(ifld)%Tindex
              Tintrp(1)=DataSet(Icomp)%Export(ifld)%Tintrp(1)
              Tintrp(2)=DataSet(Icomp)%Export(ifld)%Tintrp(2)
!
!  Set linear-interpolation factors.  To avoid roundoff, the fractional
!  days interval are rounded to the nearest millisecond interger toward
!  zero in the time interpolation weights.
!
              DayScale=86400.0_dp*1000.0_dp   ! days to milliseconds
              it1=3-Tindex
              it2=Tindex
              fac1=ANINT((Tintrp(it2)-Tcurrent)*DayScale)
              fac2=ANINT((Tcurrent-Tintrp(it1))*DayScale)
!
!  Time-interpolate from gridded or point data.
!
              IF (((fac1*fac2).ge.0.0_dp).and.                          &
     &            ((fac1+fac2).gt.0.0_dp)) THEN
                fac=1.0_dp/(fac1+fac2)
                w1=fac*fac1
                w2=fac*fac2
                IF (Nvdim.eq.2) THEN                      ! 2D variable
                  Im=DataSet(Icomp)%Export(ifld)%Vsize(1)
                  Jm=DataSet(Icomp)%Export(ifld)%Vsize(2)
                  IF (.not.allocated(DataSet(Icomp)%                    &
     &                                      Export(ifld)%A2d)) THEN
                    allocate ( DataSet(Icomp)%                          &
     &                                Export(ifld)%A2d(Im,Jm) )
                  END IF
                  DO j=1,Jm
                    DO i=1,Im
                      Fval=w1*DataSet(Icomp)%Export(ifld)%              &
     &                                             A2dG(i,j,it1)+       &
     &                     w2*DataSet(Icomp)%Export(ifld)%              &
     &                                             A2dG(i,j,it2)
                      DataSet(Icomp)%Export(ifld)%A2d(i,j)=Fval
                    END DO
                  END DO
                ELSE IF (Nvdim.eq.3) THEN                 ! 3D variable
                  Im=DataSet(Icomp)%Export(ifld)%Vsize(1)
                  Jm=DataSet(Icomp)%Export(ifld)%Vsize(2)
                  Km=DataSet(Icomp)%Export(ifld)%Vsize(3)
                  IF (.not.allocated(DataSet(Icomp)%                    &
     &                                       Export(ifld)%A3d)) THEN
                    allocate ( DataSet(Icomp)%                          &
     &                                Export(ifld)%A3d(Im,Jm,Km) )
                  END IF
                  DO k=1,Km
                    DO j=1,Jm
                      DO i=1,Im
                        Fval=w1*DataSet(Icomp)%Export(ifld)%            &
     &                                               A3dG(i,j,k,it1)+   &
     &                       w2*DataSet(Icomp)%Export(ifld)%            &
     &                                               A3dG(i,j,k,it2)
                        DataSet(Icomp)%Export(ifld)%A3d(i,j,k)=Fval
                      END DO
                    END DO
                  END DO
                END IF
!
!  Unable to set-up requested field.  Activate error flag to quit.
!
              ELSE
                IF (localPET.eq.0) THEN
                  WRITE (dataout,10) TRIM(Vname), Tcurrent,             &
     &                               Tmin, Tmax,                        &
     &                               Tstr, Tend,                        &
     &                               Tintrp(it1), Tintrp(it2),          &
     &                               fac1, fac2
                END IF
                exit_flag=2
                IF (FoundError(exit_flag, NoError, __LINE__,            &
     &                         MyFile)) THEN
                  rc=ESMF_RC_VAL_OUTOFRANGE
                  RETURN
                END IF
              END IF
            END DO
          END IF
        END IF
      END DO FIELD_LOOP
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_TimeInterp',         &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (/,' DATA_TimeInterp - current coupling time',             &
     &          ' exceeds ending value for variable: ',a,               &
     &        /,14x,'Tcurrent  = ',f15.4,                               &
     &        /,14x,'Data Tmin = ',f15.4,2x,'Data Tmax = ',f15.4,       &
     &        /,14x,'Data Tstr = ',f15.4,2x,'Data Tend = ',f15.4,       &
     &        /,14x,'TINTRP1   = ',f15.4,2x,'TINTRP2   = ',f15.4,       &
     &        /,14x,'FAC1      = ',f15.4,2x,'FAC2      = ',f15.4)

      RETURN
      END SUBROUTINE Data_TimeInterp
!
      SUBROUTINE DATA_ncread (Tcurrent, FirstPass, localPET,            &
     &                        IsUpdated, rc)
!
!=======================================================================
!                                                                      !
!  This routine read DATA component field to export from NetCDF source !
!  file at the appropriate time.  The data is loaded to the snapshots  !
!  arrays for time interpolation elsewhere:                            !
!                                                                      !
!    DataSet(Icomp)%Export(ifld)%A2dG(:,:,Tindex)           2D field   !
!    DataSet(Icomp)%Export(ifld)%A3dG(:,:,:,Tindex)         3D field   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Tcurrent   Current time in days since reference date (real)      !
!     FirstPass  Switch indicating initialization or restart phase     !
!                  (logical)                                           !
!     localPET   Local Persistent Execution Thread (integer)           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     DataSet    Updated DATA component structure                      !
!     IsUpdated  Set to TRUE if new fields have been read (logical)    !
!     rc         Return code flag (integer)                            !
!                                                                      !
!  WARNING:                                                            !
!                                                                      !
!  This routine uses ROMS NetCDF processing framework.                 !
!                                                                      !
!=======================================================================
!
      USE mod_netcdf
!
      USE mod_iounits,   ONLY : SourceFile
      USE mod_scalars,   ONLY : NoError, exit_flag
!
      USE dateclock_mod, ONLY : caldate, time_string
      USE strings_mod,   ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(in)  :: FirstPass
      logical, intent(out) :: IsUpdated
!
      integer, intent(in)  :: localPET
      integer, intent(out) :: rc
!
      real(dp), intent(in) :: Tcurrent
!
!  Local variable declarations.
!
      logical :: Lcycle, Linquire, Lmulti, LZL, RecLast
!
      integer :: Icomp, Ndims, Nfields, Nfiles, Nvdim, ZL, ifld
      integer :: ncid, Tid, Vid
      integer :: Nrec, Tindex, Trec
      integer :: Nx, Ny, Nz, i
      integer :: lend, lstr, lvar
      integer :: ROMScomm
      integer :: MyDateVec(5)

      integer, parameter :: imodel = 1    ! for compatibility with ROMS
      integer, parameter :: ng = 1        ! used routines
!
      real(r8) :: Vmax, Vmin

      real(dp) :: Clength, Tdelta, Tmax, Tmin, Tmono, Tscale, Tsec, Tval
      real(dp) :: MySeconds
      real(dp) :: V_time(2)
!
      character (len=15 ) :: Zlabel
      character (len=20 ) :: Ctarget, nc_vname, nc_tname, shortname
      character (len=22 ) :: Tcode
      character (len=100) :: T_name, V_name, Vunits
      character (len=256) :: ncname, longname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_ncread"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_ncread',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
      SourceFile=MyFile
!
!=======================================================================
!  If appropriate, read in new data.
!=======================================================================
!
!  Read in DATA component fields from source NetCDF files. The DATA
!  model can be connected to any of the activated ESM components.
!
      IsUpdated=.FALSE.
!
      FIELD_LOOP : DO Icomp=1,Nmodels
        IF (Icomp.ne.Idata) THEN
          Nfiles=DataSet(Icomp)%Nfiles
          IF (MODELS(Icomp)%IsActive.and.(Nfiles.gt.0)) THEN
            Nfields=DataSet(Icomp)%Nfields
            DO ifld=1,Nfields
!
!  If appropriate, inquire about the contents of input NetCDF file and
!  fill Information File Structure (IFS).
!
!  If the switch "RecLast" is true, we need to inquire information about
!  the next multifile for the UPPER time snapshot data. It implies that
!  the last record in the file was processed for LOWER time snapshot,
!  previously.  Notice that sfter the inquiry, "LastRec" is deactivated.
!
              Lcycle =DataSet(Icomp)%Export(ifld)%Lcycle
              Lmulti =DataSet(Icomp)%Export(ifld)%Lmulti
              RecLast=DataSet(Icomp)%Export(ifld)%LastRec
              Tmax   =DataSet(Icomp)%Export(ifld)%Tmax
              Tmono  =DataSet(Icomp)%Export(ifld)%Tmono
!
              Linquire=Lmulti.and.                                      &
     &                 (RecLast.or.(.not.Lcycle.and.(Tmax.lt.Tcurrent)))
!

              IF (Linquire) THEN
                nc_vname=Models(Idata)%ExportField(ifld)%nc_vname
                nc_tname=Models(Idata)%ExportField(ifld)%nc_tname
                CALL DATA_inquiry (ifld, nc_vname, nc_tname,            &
     &                             Tcurrent,                            &
     &                             DataSet(Icomp)%Export, Nfields,      &
     &                             DataSet(Icomp)%IFS, Nfiles,          &
     &                             Lmulti, localPET, rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
                  RETURN
                END IF
                IF (RecLast) THEN
                  DataSet(Icomp)%Export(ifld)%LastRec=.FALSE.
                END IF
              END IF
!
!-----------------------------------------------------------------------
!  If needed, read in the field data time-snapshot.
!-----------------------------------------------------------------------
!
              IF ((Tmono.le.Tcurrent).or.FirstPass) THEN
!
!  Load properties of export field to read from source NetCDF file.
!
                Ctarget  =DataSet(Icomp)%Ctarget(ifld)
                ncname   =DataSet(Icomp)%Export(ifld)%ncfile
                T_name   =DataSet(Icomp)%Export(ifld)%Tname
                V_name   =DataSet(Icomp)%Export(ifld)%Vname
                shortname=DataSet(Icomp)%Field(ifld)
                longname =DataSet(Icomp)%Export(ifld)%Vdescriptor
                Vunits   =DataSet(Icomp)%Export(ifld)%Vunits
                ncid     =DataSet(Icomp)%Export(ifld)%ncid
                Vid      =DataSet(Icomp)%Export(ifld)%Vid
                Tid      =DataSet(Icomp)%Export(ifld)%Tid
                Tindex   =DataSet(Icomp)%Export(ifld)%Tindex
                Nrec     =DataSet(Icomp)%Export(ifld)%Nrec
                Trec     =DataSet(Icomp)%Export(ifld)%Trec
                Tscale   =DataSet(Icomp)%Export(ifld)%Tscale
                Clength  =DataSet(Icomp)%Export(ifld)%Clength
                Tmin     =DataSet(Icomp)%Export(ifld)%Tmin
                V_time(1)=DataSet(Icomp)%Export(ifld)%Vtime(1)
                V_time(2)=DataSet(Icomp)%Export(ifld)%Vtime(2)
!
!  Advance time record to process.
!
                IF (Lcycle) THEN
                  Trec=MOD(Trec,Nrec)+1
                ELSE
                  Trec=Trec+1
                END IF
                DataSet(Icomp)%Export(ifld)%Trec=Trec
                LZL=.FALSE.
                Tval=0.0_dp
!
!  Process if time record is available in NetCDF file.
!
                IF (Trec.le.Nrec) THEN
!
!  Set rolling index for two-time record storage of input data for the
!  time interpolation elsewhere.
!
                  Tindex=3-Tindex
                  DataSet(Icomp)%Export(ifld)%Tindex=Tindex
!
!  Read in time coordinate and scale it to day units.
!
                  CALL netcdf_get_time (ng, imodel, ncname,             &
     &                                  TRIM(T_name),                   &
     &                                  ReferenceDateNumber,            &
     &                                  Tval,                           &
     &                                  ncid = ncid,                    &
     &                                  start = (/Trec/),               &
     &                                  total = (/1/))
                  IF (FoundError(exit_flag, NoError, __LINE__,          &
     &                           MyFile)) THEN
                    IF (localPET.eq.0) WRITE (dataout,10) TRIM(T_name), &
     &                                                    Trec,         &
     &                                                    TRIM(ncname)
                    rc=ESMF_RC_FILE_READ
                    RETURN
                  END IF
                  Tval=Tval*Tscale             ! scaled to day units
                  V_time(Tindex)=Tval
!
!  Activate "LastRec" switch if processing the LAST record of the file
!  for the LOWER time snapshot. We need to get the UPPER time snapshot
!  from NEXT multifile.
!
                  IF ((Trec.eq.Nrec).and.(Tval.le.Tcurrent)) THEN
                    DataSet(Icomp)%Export(ifld)%LastRec=.TRUE.
                  END IF
!
!  Read in field. Allocate snapshot array, if necessary. Notice that it
!  is possible to read a particular depth level (ZL) from a 3D field.
!  For example, reading sea surface temperature for a 3D temperature
!  variable where ZL is the depth level index for the surface.
!
                  Nvdim=DataSet(Icomp)%Export(ifld)%Nvdim
                  Ndims=SIZE(DataSet(Icomp)%Export(ifld)%Vsize)-1
!
                  IF ((Nvdim.eq.2).and.(Ndims.eq.2)) THEN       ! 2D var
                    Nx=DataSet(Icomp)%Export(ifld)%Vsize(1)
                    Ny=DataSet(Icomp)%Export(ifld)%Vsize(2)
                    IF (.not.allocated(DataSet(Icomp)%                  &
     &                                        Export(ifld)%A2dG)) THEN
                      allocate ( DataSet(Icomp)%                        &
     &                                  Export(ifld)%A2dG(Nx,Ny,2) )
                    END IF
                    Vmin=0.0_r8
                    Vmax=0.0_r8
                    CALL netcdf_get_fvar (ng, imodel, ncname,           &
     &                                    TRIM(V_name),                 &
     &                  DataSet(Icomp)%Export(ifld)%A2dG(:,:,Tindex),   &
     &                                    ncid = ncid,                  &
     &                                    start = (/1,1,Trec/),         &
     &                                    total = (/Nx,Ny,1/),          &
     &                                    min_val = Vmin,               &
     &                                    max_val = Vmax)
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      IF (localPET.eq.0)                                &
     &                  WRITE (dataout,10) TRIM(V_name), Trec,          &
     &                                     TRIM(ncname)
                      rc=ESMF_RC_FILE_READ
                      RETURN
                    END IF
                    DataSet(Icomp)%Export(ifld)%Vmin=Vmin
                    DataSet(Icomp)%Export(ifld)%Vmax=Vmax
                  ELSE IF ((Nvdim.eq.2).and.(Ndims.eq.3)) THEN  ! 2D var
                    Nx=DataSet(Icomp)%Export(ifld)%Vsize(1)
                    Ny=DataSet(Icomp)%Export(ifld)%Vsize(2)
                    ZL=DataSet(Icomp)%Export(ifld)%Zlevel
                    IF (.not.allocated(DataSet(Icomp)%                  &
     &                                        Export(ifld)%A2dG)) THEN
                      allocate ( DataSet(Icomp)%                        &
     &                                  Export(ifld)%A2dG(Nx,Ny,2) )
                    END IF
                    Vmin=0.0_r8
                    Vmax=0.0_r8
                    CALL netcdf_get_fvar (ng, imodel, ncname,           &
     &                                    TRIM(V_name),                 &
     &                    DataSet(Icomp)%Export(ifld)%A2dG(:,:,Tindex), &
     &                                    ncid = ncid,                  &
     &                                    start = (/1,1,ZL,Trec/),      &
     &                                    total = (/Nx,Ny,1,1/),        &
     &                                    min_val = Vmin,               &
     &                                    max_val = Vmax)
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      IF (localPET.eq.0)                                &
                        WRITE (dataout,10) TRIM(V_name), Trec,          &
     &                                     TRIM(ncname)
                      rc=ESMF_RC_FILE_READ
                      RETURN
                    END IF
                    DataSet(Icomp)%Export(ifld)%Vmin=Vmin
                    DataSet(Icomp)%Export(ifld)%Vmax=Vmax
                    WRITE (Zlabel,'(a,i2.2)') 'Level = ', ZL
                    LZL=.TRUE.
                  ELSE IF ((Nvdim.eq.3).and.(Ndims.eq.3)) THEN  ! 3D var
                    Nx=DataSet(Icomp)%Export(ifld)%Vsize(1)
                    Ny=DataSet(Icomp)%Export(ifld)%Vsize(2)
                    Nz=DataSet(Icomp)%Export(ifld)%Vsize(3)
                    IF (.not.allocated(DataSet(Icomp)%                  &
     &                                        Export(ifld)%A3dG)) THEN
                      allocate ( DataSet(Icomp)%                        &
     &                                  Export(ifld)%A3dG(Nx,Ny,Nz,2) )
                    END IF
                    Vmin=0.0_r8
                    Vmax=0.0_r8
                    CALL netcdf_get_fvar (ng, imodel, ncname,           &
     &                                    TRIM(V_name),                 &
     &                  DataSet(Icomp)%Export(ifld)%A3dG(:,:,:,Tindex), &
     &                                    ncid = ncid,                  &
     &                                    start = (/1,1,1,Trec/),       &
     &                                    total = (/Nx,Ny,Nz,1/),       &
     &                                    min_val = Vmin,               &
     &                                    max_val = Vmax)
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) THEN
                      IF (localPET.eq.0)                                &
     &                  WRITE (dataout,10) TRIM(V_name), Trec,          &
     &                                     TRIM(ncname)
                      rc=ESMF_RC_FILE_READ
                      RETURN
                    END IF
                    DataSet(Icomp)%Export(ifld)%Vmin=Vmin
                    DataSet(Icomp)%Export(ifld)%Vmax=Vmax
                  END IF
                  CALL caldate (Tval,                                   &
     &                          yy_i = MyDateVec(1),                    &
     &                          mm_i = MyDateVec(2),                    &
     &                          dd_i = MyDateVec(3),                    &
     &                          h_i  = MyDateVec(4),                    &
     &                          m_i  = MyDateVec(5),                    &
     &                          s_dp = MySeconds)
                  DO i=1,5
                    DataSet(Icomp)%Export(ifld)%Date(i,Tindex)=         &
     &                                            REAL(MyDateVec(i),dp)
                  END DO
                  DataSet(Icomp)%Export(ifld)%Date(6,Tindex)=MySeconds
                  IsUpdated=.TRUE.
                  lstr=SCAN(ncname,'/',BACK=.TRUE.)+1
                  lend=LEN_TRIM(ncname)
                  lvar=MIN(43,LEN_TRIM(longname))
                  Tsec=Tval*86400.0_dp              ! scaled to seconds
                  CALL time_string (Tsec, Tcode)
                  DataSet(Icomp)%Export(ifld)%DateString(Tindex)=Tcode
                  IF (localPET.eq.0) THEN
                    IF (LZL) THEN
                      WRITE (dataout,20) TRIM(V_name), Tcode,           &
     &                                   TRIM(shortname),               &
     &                                   TRIM(longname),                &
     &                                   TRIM(Vunits), TRIM(Ctarget),   &
     &                                   Trec, Tindex,                  &
     &                                   ncname(lstr:lend),             &
     &                                   Tmin, Tmax, Tval, Vmin, Vmax,  &
     &                                   TRIM(Zlabel)
                    ELSE
                      WRITE (dataout,30) TRIM(V_name), Tcode,           &
     &                                   TRIM(shortname),               &
     &                                   TRIM(longname),                &
     &                                   TRIM(Vunits), TRIM(Ctarget),   &
     &                                   Trec, Tindex,                  &
     &                                   ncname(lstr:lend),             &
     &                                   Tmin, Tmax, Tval, Vmin, Vmax
                    END IF
                  END IF
                END IF
!
!  Increment the local time variable "Tmono" by the interval between
!  snapshots. If the interval is negative, indicating cycling, add in
!  a cycle length.  Load time value (sec) into "Tintrp" which used
!  during interpolation between snapshots.
!
                Tdelta=V_time(Tindex)-V_time(3-Tindex)
                IF (Lcycle.and.(Tdelta.lt.0.0_dp)) THEN
                  Tdelta=Tdelta+Clength
                END IF
                Tmono=Tmono+Tdelta
                DataSet(Icomp)%Export(ifld)%Tmono=Tmono
                DataSet(Icomp)%Export(ifld)%Tintrp(Tindex)=Tmono
                DataSet(Icomp)%Export(ifld)%Vtime(Tindex)=Tval
              END IF
            END DO
          END IF
        END IF
      END DO FIELD_LOOP
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_ncread',             &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
  10  FORMAT (/,' DATA_ncread - error while reading variable: ',a,2x,   &
     &        ' at TIME record = ',i7,/,15x,'in file: ',a)
  20  FORMAT (3x,' DATA_ncread - ESMF: reading ''',a,''',',t68,a,/,     &
     &        7x,17x,'''',a,''': ',a,2x,'(',a,')',/,2x,17x,             &
     &        '(Target: ',a,', Rec=',i7.7,', SnapshotIndex=',i1,        &
     &        ', File: ',a,')',/,19x,                                   &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',                   &
     &        t71, 't = ', f15.4 ,/, 19x,                               &
     &        '(Dmin= ', 1p,e15.8,0p,' Dmax= ',1p,e15.8,0p,')',         &
     &        t71,a)
  30  FORMAT (3x,' DATA_ncread - ESMF: reading ''',a,''',',t68,a,/,     &
     &        7x,17x,'''',a,''': ',a,2x,'(',a,')',/,2x,17x,             &
     &        '(Target: ',a,', Rec=',i7.7,', SnapshotIndex=',i1,        &
     &        ', File: ',a,')',/,19x,                                   &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',                   &
     &        t71, 't = ', f15.4 ,/, 19x,                               &
     &        '(Dmin= ', 1p,e15.8,0p,' Dmax= ',1p,e15.8,0p,')')

      RETURN

      END SUBROUTINE DATA_ncread
!
      SUBROUTINE DATA_multifile (Tcurrent, IFS, Nfiles, localPET, rc)
!
!=======================================================================
!                                                                      !
!  This routine checks DATA model input NetCDF multi-files and         !
!  sets several parameters in the file information structure so        !
!  the appropriate file is selected during initialization or restart   !
!                                                                      !
!  Here, multi-file implies that the time records for a particular     !
!  field can be split into several NetCDF files.                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Tcurrent    Current time in days since reference date (real)     !
!     IFS         Input Files Structure, TYPE(T_IO)                    !
!     Nfiles      Number of files in structure (vector)                !
!     localPET    Local Persistent Execution Thread (integer)          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     IFS         Updated Input Files Structure, TYPE(T_IO)            !
!     rc          Return code flag (integer)                           !
!                                                                      !
!=======================================================================
!
      USE mod_iounits, ONLY : SourceFile
      USE mod_scalars, ONLY : NoError, exit_flag, spval
!
      USE dateclock_mod, ONLY : time_string
      USE strings_mod,   ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: Nfiles, localPET
      integer, intent(out) :: rc
!
      real (dp) :: Tcurrent
!
      TYPE(T_IO), intent(inout) :: IFS(Nfiles)
!
!  Local variable declarations.
!
      logical :: Lcheck, foundit
!
      integer :: Fcount, Mfiles, i, ifile, lstr
!
      real(dp) :: TimeStrDay, TimeEndDay
      real(dp) :: TimeStrSec, TimeEndSec
      real(dp) :: Tmax, Tmin, Tscale
!
      character (len=1), parameter :: blank = ' '
      character (len= 22) :: F_code, I_code, Tmin_code, Tmax_code
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_multifile"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_multifile',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Process DATA model file structure.
!-----------------------------------------------------------------------
!
!  Get simulation starting and ending times: seconds since reference
!  time.  Then, compute days since reference time.
!
      TimeStrSec=Tcurrent*86400.0_dp
      TimeEndSec=ClockInfo(Idata)%Time_Stop-                            &
     &           ClockInfo(Idata)%Time_Reference
!
      TimeStrDay=TimeStrSec/86400.0_dp
      TimeEndDay=TimeEndSec/86400.0_dp
!
!  Get simulation start and ending time string.
!
      CALL time_string (TimeStrSec, I_code)
      CALL time_string (TimeEndSec, F_code)
!
!  Set available minimum and maximum time coordinates.
!
      DO i=1,Nfiles
        Mfiles=IFS(i)%Nfiles  ! number of multi-files within source file
        DO ifile=1,Mfiles
          ncname=IFS(i)%files(ifile)
          foundit=DATA_checkfile(ncname, Tmin, Tmax, Tscale, localPET,  &
     &                           Lcheck)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_CANNOT_SET
            RETURN
          END IF
          IFS(i)%time_min(ifile)=Tmin
          IFS(i)%time_max(ifile)=Tmax
        END DO
!
!  Set the appropriate file counter to use during initialization or
!  restart.
!
        Fcount=0
        IF (Lcheck) THEN
          DO ifile=1,Mfiles
            Tmin=Tscale*IFS(i)%time_min(ifile)
            IF (TimeStrDay.ge.Tmin) THEN
              Fcount=ifile
            END IF
          END DO
        ELSE
          Fcount=1
        END IF
!
!  Initialize other structure parameters or issue an error if data does
!  not include initalization time.
!
        IF (Fcount.gt.0) THEN
          IFS(i)%Fcount=Fcount
          ncname=IFS(i)%files(Fcount)
          lstr=LEN_TRIM(ncname)
          IFS(i)%name=TRIM(ncname)
          IFS(i)%base=ncname(1:lstr-3)
        ELSE
          IF ((localPET.eq.0).and.Lcheck) THEN
            WRITE (dataout,10) 'Data Model', I_code
            DO ifile=1,Mfiles
              Tmin=Tscale*IFS(i)%time_min(ifile)
              Tmax=Tscale*IFS(i)%time_max(ifile)
              CALL time_string (Tmin*86400.0_dp, Tmin_code)
              CALL time_string (Tmax*86400.0_dp, Tmax_code)
              WRITE (dataout,20) Tmin_code, Tmax_code,                  &
     &                           TRIM(IFS(i)%files(ifile))
            END DO
          END IF
          exit_flag=4
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_CANNOT_SET
            RETURN
          END IF
        END IF
!
!  Check if there is forcing data up to the end of the simulation.
!
        IF (Lcheck) THEN
          Tmax=Tscale*IFS(i)%time_max(Mfiles)
          IF (TimeEndDay.gt.Tmax) THEN
            CALL time_string (Tmax*86400.0_dp, Tmax_code)
            IF (localPET.eq.0) THEN
              WRITE (dataout,30) 'Data Model',                          &
     &                           TRIM(IFS(i)%files(Mfiles)),            &
     &                           'last  ', Tmax_code, F_code
            END IF
            exit_flag=4
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     MyFile)) THEN
              rc=ESMF_RC_CANNOT_SET
              RETURN
            END IF
          END IF
        END IF
      END DO
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_multifile',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
 10   FORMAT (/,' DATA_MULTIFILE - Error while processing ', a,         &
     &        ' multi-files: ',/,18x,'data does not include',           &
     &        ' initialization time = ', a,/)
 20   FORMAT (3x,a,2x,a,5x,a)
 30   FORMAT (/,' DATA_MULTIFILE - Error while checking input ', a,     &
     &        ' file:',/,18x,a,/,18x,                                   &
     &        a,'data time record available is for day: ',a,/,18x,      &
     &        'but  data is needed to finish run until day: ',a)

      RETURN
      END SUBROUTINE DATA_multifile
!
      FUNCTION DATA_checkfile (ncname, Tmin, Tmax, Tscale, localPET,    &
     &                         Lcheck) RESULT (foundit)
!
!=======================================================================
!                                                                      !
!  This logical function scans the variables of the provided input     !
!  NetCDF for the time record variable and gets its range of values.   !
!  It used elsewhere to determine which input NetCDF multi-file is     !
!  needed for initialization or restart.                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ncname       NetCDF file name to process (string)                !
!     localPET     Local Persistent Execution Thread (integer)         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Tmin         Available minimum time variable value               !
!     Tmax         Available maximum time variable value               !
!     Tscale       Scale to convert time variable units to days        !
!     Lcheck       Switch to indicate that the time range needs to be  !
!                    checked by the calling routine                    !
!     foundit      The value of the result is TRUE/FALSE if the        !
!                    time variable is found or not                     !
!                                                                      !
!  WARNING:                                                            !
!                                                                      !
!  This routine uses ROMS NetCDF processing framework.                 !
!                                                                      !
!=======================================================================
!
      USE mod_netcdf
!
      USE mod_iounits,  ONLY : SourceFile
      USE mod_scalars,  ONLY : NoError, exit_flag
      USE strings_mod,  ONLY : FoundError, lowercase
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: Lcheck
!
      integer, intent(in ) :: localPET
!
      character (*), intent(in) :: ncname
!
      real(dp), intent(out) :: Tmin, Tmax, Tscale
!
!  Local variable declarations.
!
      logical :: Lcycle, Lperpetual, foundit
!
      integer :: Nrec, TvarID, i, j, ncid, ncvid, nvdim, nvatt
      integer :: ROMScomm
!
      integer, parameter :: imodel = 1    ! for compatibility with ROMS
      integer, parameter :: ng = 1        ! used routines
!
      character (len= 40) :: Tunits, TvarName
      character (len=100) :: blank, long_name, units

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_checkfile"
!
!-----------------------------------------------------------------------
!  Check if requested time is within the NetCDF file dataset.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_checkfile',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      SourceFile=MyFile
!
      foundit=.FALSE.
      Lcheck=.TRUE.
      Lcycle=.FALSE.
      Lperpetual=.FALSE.
      Tscale=1.0_dp                        ! days
      Tmin= MISSING_dp
      Tmax=-MISSING_dp
      DO i=1,LEN(blank)
        blank(i:i)=' '
      END DO
      IF (localPET.eq.0) THEN
        WRITE (dataout,10) TRIM(ncname)
      END IF
!
!  Open NetCDF file for reading.
!
      CALL netcdf_open (ng, imodel, ncname, 0, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               MyFile)) THEN
        IF (localPET.eq.0) THEN
          WRITE (dataout,20) TRIM(ncname)
        END IF
        RETURN
      END IF
!
!  Inquire about all the variables
!
      CALL netcdf_inq_var (ng, imodel, ncname,                          &
     &                     ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               MyFile)) RETURN
!
!  Search for the time variable: any 1D array variable with the string
!  'time' in the variable name.
!
      DO i=1,n_var
        IF ((INDEX(TRIM(lowercase(var_name(i))),'time').ne.0).and.      &
     &            (var_ndim(i).eq.1)) THEN
          TvarName=TRIM(var_name(i))
          foundit=.TRUE.
          EXIT
        END IF
      END DO
!
!  If not found, scan all the 1D array variables and inquire the
!  'long_name' attribute for the string time and the 'units' attribute
!  for the strin 'day' or 'second'.
!
      IF (.not.foundit) THEN
        DO i=1,n_var
          IF (var_ndim(i).eq.1) THEN
            CALL netcdf_inq_var (ng, imodel, ncname,                    &
     &                           ncid = ncid,                           &
     &                           MyVarName = TRIM(var_name(i)),         &
     &                           VarID = ncvid,                         &
     &                           nVarDim = nvdim,                       &
     &                           nVarAtt = nvatt)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     MyFile)) RETURN
            long_name=blank
            units=blank
            DO j=1,nvatt
              IF (TRIM(var_Aname(j)).eq.'long_name') THEN
                long_name=TRIM(var_Achar(j))
              ELSE IF (TRIM(var_Aname(j)).eq.'units') THEN
                units=TRIM(var_Achar(j))
              END IF
            END DO
            IF ((INDEX(TRIM(lowercase(long_name)),'time').ne.0).and.    &
     &          ((INDEX(TRIM(lowercase(units)),'day').ne.0).or.         &
     &           (INDEX(TRIM(lowercase(units)),'second').ne.0))) THEN
              TvarName=TRIM(var_name(i))
              foundit=.TRUE.
              EXIT
            END IF
          END IF
        END DO
      END IF
!
!  Issue and error if time variable not found.
!
      IF (.not.foundit) THEN
        IF (localPET.eq.0) THEN
          WRITE (dataout,30) TRIM(ncname)
        END IF
        exit_flag=4
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) RETURN
      END IF
!
!  Inquire about time variable.
!
      CALL netcdf_inq_var (ng, imodel, ncname,                          &
     &                     ncid = ncid,                                 &
     &                     MyVarName = TRIM(TvarName),                  &
     &                     VarID = TvarID,                              &
     &                     nVarDim = nvdim,                             &
     &                     nVarAtt = nvatt)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               MyFile)) RETURN
!
!  Set number of records available and check the 'units' attribute.
!  Also, set output logical switch 'Lcheck' for the calling to check
!  the available data time range. For example, we need to check it
!  there is enough data to finish the simulation.  Notice that for
!  data with 'cycle_length', Lcheck = FALSE.  Also,  Lcheck = FALSE
!  for perpetual time axis: the 'calendar' attribute is 'none' or
!  the number of records in the time dimension is one (Nrec=1).
!
      Nrec=var_Dsize(1)              ! time is a 1D array
      DO i=1,nvatt
        IF (TRIM(var_Aname(i)).eq.'units') THEN
          Tunits=TRIM(var_Achar(i))
          IF (INDEX(TRIM(var_Achar(i)),'day').ne.0) THEN
            Tscale=1.0_dp
          ELSE IF (INDEX(TRIM(var_Achar(i)),'hour').ne.0) THEN
            Tscale=24.0_dp
          ELSE IF (INDEX(TRIM(var_Achar(i)),'second').ne.0) THEN
            Tscale=86400.0_dp
          END IF
        ELSE IF (TRIM(var_Aname(i)).eq.'calendar') THEN
          IF ((Nrec.eq.1).or.                                           &
     &        (INDEX(TRIM(var_Achar(i)),'none').ne.0)) THEN
            Lperpetual=.TRUE.
          END IF
        ELSE IF (TRIM(var_Aname(i)).eq.'cycle_length') THEN
          Lcycle=.TRUE.
        END IF
      END DO
!
!  Turn off the checking of time range if cycling, perpectual, or
!  spectral time axis.
!
      IF (Lcycle.or.Lperpetual.or.(Nrec.eq.1)) THEN
        Lcheck=.FALSE.
      END IF
!
!  Read in time variable minimum and maximun values (input time units).
!
      CALL netcdf_get_time (ng, imodel, ncname, TvarName,               &
     &                      ReferenceDateNumber,                        &
     &                      Tmin,                                       &
     &                      ncid = ncid,                                &
     &                      start = (/1/),                              &
     &                      total = (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               MyFile)) RETURN
!
      CALL netcdf_get_time (ng, imodel, ncname, TvarName,               &
     &                      ReferenceDateNumber,                        &
     &                      Tmax,                                       &
     &                      ncid = ncid,                                &
     &                      start = (/Nrec/),                           &
     &                      total = (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               MyFile)) RETURN
!
!  Close NetCDF file.
!
      CALL netcdf_close (ng, imodel, ncid, ncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               MyFile)) RETURN
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_checkfile',          &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
 10   FORMAT (' DATA_checkfile   - inquiring range of time records in', &
     &        ' input NetCDF file:',/,20x,a)
 20   FORMAT (/, ' DATA_CHECKFILE - unable to open grid NetCDF file: ', &
     &        a)
 30   FORMAT (/, ' DATA_CHECKFILE - unable to find time variable in ',  &
     &        ' input NetCDF file:',/,18x,a,/,18x,                      &
     &        'variable name does not contains the "time" string.')

      RETURN
      END FUNCTION DATA_checkfile
!
      SUBROUTINE DATA_inquiry (ifield, ncvname, nctname, Tcurrent,      &
     &                         Export, Nfields, IFS, Nfiles,            &
     &                         Lmulti, localPET, rc)
!
!=======================================================================
!                                                                      !
!  This routine inquires DATA model variable to process from input     !
!  NetCDF files.                                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ifield     Field index in Export structure to process (integer)  !
!     ncvname    DATA model NetCDF field variable to process (string)  !
!     nctname    DATA model NetCDF time variable to process (string)   !
!     Tcurrent   Current time in days since reference date (real)      !
!     Export     DATA component structure, TYPE(ESM_Data)              !
!     Nfields    Number of fields in Export structure (interger)       !
!     IFS        DATA model Input Files Structure, TYPE(T_IO)          !
!     Nfiles     Number of files in IFS structure (integer)            !
!     Lmuti      Switch indicating multifiles for field (logical)      !
!     localPET   Local Persistent Execution Thread (integer)           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Export     Updated DATA component structure                      !
!     rc         Return code flag (integer)                            !
!                                                                      !
!  WARNING:                                                            !
!                                                                      !
!  This routine uses ROMS NetCDF managing framework.                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_netcdf
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,  ONLY : FoundError, lowercase, uppercase
!
!  Imported variable declarations.
!
      logical, intent(in)  :: Lmulti
!
      integer, intent(in)  :: ifield, Nfields, Nfiles
      integer, intent(in)  :: localPET
      integer, intent(out) :: rc
!
      real(dp) :: Tcurrent
!
      TYPE(ESM_Data), intent(inout) :: Export(Nfields)
      TYPE(T_IO),     intent(inout) :: IFS(Nfiles)
!
      character (len=*), intent(in) :: ncvname, nctname
!
!  Local variable declarations.
!
      logical :: CloseFile, Lcycle, Linside, LowerBound, Upperbound
      logical :: foundAtt(1), foundit, got_coord, got_var, got_time
!
      integer :: Fcount, Nrec, Tid, Tindex, Trec, Vid, Zlevel
      integer :: i, ifile, j, lstr
      integer :: ncid, ntatt, ntdim, nvatt, nvdim
      integer :: clen, iblank, ie, is
      integer :: ROMScomm
      integer :: Vsize(4)
!
      integer, parameter :: imodel = 1    ! for compatibility with ROMS
      integer, parameter :: ng = 1        ! used routines
!
      real(dp) :: Clength, Tday, Tend, Tmax, Tmin, Tmono, Tscale, Tstr
      real(dp) :: scale, tstart
!
      real(dp), allocatable :: TimeValue(:)
!
      character (len=1  ), parameter :: blank = ' '
      character (len=3  )  :: label
      character (len=20 )  :: coordinates(5)
      character (len=40 )  :: AttName(1), T_name
      character (len=100)  :: Cstring, Tunits
      character (len=256)  :: Fname
      character (len=2048) :: AttValue(1)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_inquiry"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_inquiry',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about DATA model variable to process
!-----------------------------------------------------------------------
!
!  Initialize local variables.
!
      Lcycle=.FALSE.
      Linside=.FALSE.
      LowerBound=.FALSE.
      Upperbound=.FALSE.
      got_coord=.FALSE.
      got_time=.FALSE.
      got_var=.FALSE.
      Fcount=0                                        ! initialize?
      Nrec=0
      Trec=0
      Zlevel=-1
      DO i=1,LEN(ncfile)
        ncfile(i:i)=blank
      END DO
      Clength=MISSING_dp
      IF (localPET.eq.0) THEN
        WRITE (dataout,5) TRIM(ncvname)
      END IF
!
!  If multi-files, increase (decrease if backward logic) file counter
!  and set new file names.
!
      IF (Lmulti) THEN
        DO ifile=1,Nfiles
          IF (TRIM(Export(ifield)%ncfile).eq.TRIM(IFS(ifile)%name)) THEN
            Fcount=IFS(ifile)%Fcount+1
            IF ((1.gt.Fcount).and.(Fcount.gt.IFS(ifile)%Nfiles)) THEN
              IF (localPET.eq.0) THEN
                WRITE (dataout,10) TRIM(ncvname),                       &
     &                             Fcount, IFS(ifile)%Nfiles
              END IF
              exit_flag=4
              rc=ESMF_RC_NOT_VALID
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       MyFile)) RETURN
            END IF
            IFS(ifile)%Fcount=Fcount
            IFS(ifile)%name=TRIM(IFS(ifile)%files(Fcount))
            lstr=LEN_TRIM(IFS(ifile)%name)
            IFS(ifile)%base=IFS(ifile)%name(1:lstr-3)
            CALL netcdf_close (ng, imodel, IFS(ifile)%ncid)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     MyFile)) THEN
              rc=ESMF_RC_FILE_CLOSE
              RETURN
            END IF
            EXIT
          ELSE                         ! IFS(ifile)%name and Fcount
            Fcount=IFS(ifile)%Fcount   ! already updated in first field
          END IF                       ! processed for current new file
        END DO
      ELSE
        Fcount=IFS(1)%Fcount
      END IF
      IF (Fcount.eq.0) THEN
        IF (localPET.eq.0) THEN
          WRITE (dataout,20) Fcount, label, TRIM(ncvname)
        END IF
        rc=ESMF_RC_NOT_SET
        exit_flag=4
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  If several input NetCDF files (Nfiles>1), scan files until the
!  requested variable is found.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
      QUERY: DO ifile=1,Nfiles
        Fname=IFS(ifile)%name
!
!  Open NetCDF file for reading.
!
        IF (IFS(ifile)%ncid.eq.-1) THEN
          CALL netcdf_open (ng, imodel, Fname, 0, ncid)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_OPEN
            IF (localPET.eq.0) WRITE (dataout,30) TRIM(Fname)
            RETURN
          END IF
          CloseFile=.TRUE.
        ELSE
          ncid=IFS(ifile)%ncid
          CloseFile=.FALSE.
        END IF
!
!  Inquire about requested variable.
!
        CALL netcdf_inq_var (ng, imodel, Fname,                         &
     &                       ncid = ncid,                               &
     &                       MyVarName = TRIM(ncvname),                 &
     &                       SearchVar = foundit,                       &
     &                       VarID = Vid,                               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)!
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) THEN
          rc=ESMF_RC_FILE_READ
          RETURN
        END IF
!
!  Set variable information. So far, we are exporting only 2D fields
!  from the DATA Model, Export(ifield)%Nvdim=2. Only spatial dimensions
!  are counted.  All dataset have an additional time dimension. Notice
!  that 2D fields can be extracted from the full 3D field at the desired
!  depth level index.
!
!  HGA:  We need logic here when exporting 3D full fields in the future.
!        I hate hardwired code but this is postponed. Perhaps, it needs
!        to be specified in input metadata.
!
        PROBE: IF (foundit) THEN
          IF (localPET.eq.0) THEN
            WRITE (dataout,35) TRIM(Fname)
          END IF
          got_var=.TRUE.
          ncfile=Fname
          IFS(ifile)%ncid=ncid
          Export(ifield)%ncid=ncid
          Export(ifield)%Vid=Vid
          Export(ifield)%Vname=TRIM(ncvname)
          Export(ifield)%Nvdim=2                ! only 2D fields for now

          IF (.not.allocated(Export(ifield)%Vsize)) THEN
            allocate ( Export(ifield)%Vsize(nvdim) )
            Export(ifield)%Vsize(1:nvdim)=0
          END IF
          IF (.not.allocated(Export(ifield)%Dname)) THEN
            allocate ( Export(ifield)%Dname(nvdim) )
            Export(ifield)%Dname(1:nvdim)=' '
          END IF
          IF (.not.allocated(Export(ifield)%Vcoord)) THEN
            allocate ( Export(ifield)%Vcoord(nvdim) )
            Export(ifield)%Vcoord(1:nvdim)=' '
          END IF

          DO i=1,nvdim
            Export(ifield)%Dname(i)=var_Dname(i)
            Export(ifield)%Vsize(i)=var_Dsize(i)
          END DO
!
!  If singleton depth dimension, set level index to process.
!
          IF ((nvdim.eq.4).and.(var_Dsize(3).eq.1)) THEN
            Zlevel=1
          END IF
!
!  If the NetCDF file is CF compliant, the variable dimensions and the
!  space and time coordinates have the same names. Therefo, check if
!  a time varible name  with same time dimension name exist.
!
          DO i=1,n_var
            IF (TRIM(var_name(i)).eq.TRIM(var_Dname(nvdim))) THEN
              T_name=TRIM(var_name(i))
              Nrec=var_Dsize(nvdim)
              got_time=.TRUE.
              EXIT
            END IF
          END DO
!
!  Check NetCDF file global attributes to identify data source. If
!  If HyCOM Tripolar grid, reset horizontal dimensions by removing last
!  row at j=Jmax. The global HyCOM grid is rectilinear from j=1 to
!  j=2172, and we can check this by confirming that the MINVAL and
!  MAXVAL of lat(1:4500,j) are identical there. The grid is curvilinear
!  (bi-polar patch) from j=2173 to j=3298. Note that the j=3298 row is
!  a permuted copy of j=3297 (because of the way the tripole grid is
!  implemented in HYCOM), so we can discard j=3298.
!
          Export(ifield)%SpecialAction='NONE'
          DO i=1,n_gatt
            AttName(1)=TRIM(att_name(i))
            IF (att_kind(i).eq.NF90_CHAR) THEN
              CALL netcdf_get_satt (ng, imodel, Fname, nf90_global,     &
     &                              AttName, AttValue,                  &
     &                              foundAtt,                           &
     &                              ncid = ncid)
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       MyFile)) THEN
                rc=ESMF_RC_FILE_READ
                RETURN
              END IF
              IF ((INDEX(TRIM(uppercase(AttValue(1))),                  &
     &                   'HYCOM').ne.0).and.                            &
     &            ((var_Dsize(1).eq.4500).and.                          &
     &             (var_Dsize(2).eq.3298))) THEN
                Export(ifield)%SpecialAction='HYCOM TRIPOLAR GRID'
                EXport(ifield)%Vsize(2)=var_Dsize(2)-1  ! discard j=3298
                IF ((nvdim.eq.4).and.(var_Dsize(3).gt.1)) THEN
                  Zlevel=1                  ! Depth = 0 m is at level 1
                END IF
                EXIT
              END IF
            END IF
          END DO
!
!  Load depth level to process, if any.
!
          Export(ifield)%Zlevel=Zlevel
!
!  Check variable for several attributes:
!
!     "add_offset"      if present, value to add the data
!     "scale_factor"    if present, factor to mutliply the data
!     "_FillValue"      fill value for missing data
!     "coordinates"     variable spatial and temporal coordinates
!     "time"            associated time variable
!     "units"           variable units (overwrite metadata value)
!
          DO i=1,nvatt
            SELECT CASE (TRIM(var_Aname(i)))
              CASE ('add_offset')
                Export(ifield)%add_offset=var_Afloat(i)
              CASE ('scale_factor')
                Export(ifield)%scale_factor=var_Afloat(i)
              CASE ('_FillValue', 'missing_value')
                Export(ifield)%FillValue=var_Afloat(i)
              CASE ('coordinates')
                Cstring=TRIM(ADJUSTL(var_Achar(i)))
                Export(ifield)%Vcoord=TRIM(Cstring)
                got_coord=.TRUE.
              CASE ('time')
                IF (.not.got_time) THEN
                  T_name=TRIM(var_Achar(i))
                  lstr=LEN_TRIM(T_name)
                  DO j=1,n_vdim
                    IF (TRIM(var_Dname(j)).eq.T_name(1:lstr)) THEN
                      Nrec=var_Dsize(j)
                      EXIT
                    END IF
                  END DO
                  got_time=.TRUE.
                END IF
              CASE ('long_name')
                Export(ifield)%Vlongname=TRIM(var_Achar(i))
              CASE ('units')
                Export(ifield)%Vunits=TRIM(var_Achar(i))
            END SELECT
          END DO
!
!  If the "coordinates" is present, extract variables strings.
!
          Export(ifield)%Lcoord=got_coord
          IF (got_coord) THEN
            clen=LEN_TRIM(Cstring)
            is=1
            DO i=1,nvdim
              iblank=INDEX(Cstring(is:clen),' ')
              IF (iblank.eq.0) THEN
                ie=clen+1                 ! last value, add 1
              ELSE
                ie=iblank+is-1            ! includes blank index
              END IF
              coordinates(i)=Cstring(is:ie-1)
              Export(ifield)%Vcoord(i)=TRIM(coordinates(i))
              is=ie+1
            END DO
!
!  If found associated time variable, overwrite time coordinate to
!  insure that we have the correct variable.  For example, HyCOM data
!  use 'Date' in the coordinate attribute instead of the actul time
!  variable.
!
            IF (got_time) THEN
              coordinates(nvdim)=TRIM(T_name)
              Export(ifield)%Vcoord(nvdim)=TRIM(T_name)
            END IF
!
!  If not found time variable, inquire the coordinates for a variable
!  that contains the "time" string. If unsucessful, use the last
!  string in the coordinates attribute which is usually the time
!  variable.
!
            IF (.not.got_time) THEN
              DO i=1,nvdim
                IF (INDEX(TRIM(lowercase(coordinates(i))),'time')       &
     &              .ne.0) THEN
                  T_name=TRIM(coordinates(i))
                  Nrec=var_Dsize(i)
                  got_time=.TRUE.
                END IF
              END DO
              IF (.not.got_time) THEN
                T_name=TRIM(coordinates(nvdim))
                Nrec=var_Dsize(nvdim)
                got_time=.TRUE.
              END IF
            END IF
          END IF
!
!  If Nrec=0, input file is not CF compliant, check variable dimension
!  to see if the dimension contains the "time" string.
!
          IF (.not.(got_time.or.got_coord)) THEN
            DO i=1,nvdim
              IF (INDEX(TRIM(lowercase(var_Dname(i))),'time').ne.0) THEN
                T_name=TRIM(var_Dname(i))
                Nrec=var_Dsize(i)
                got_time=.TRUE.
              END IF
            END DO
          END IF
          IF (.not.got_time.and.(Nrec.eq.0)) THEN
            IF (localPET.eq.0) WRITE (dataout,40) TRIM(T_name),         &
     &                                            TRIM(ncvname),        &
     &                                            TRIM(Fname)
            rc=ESMF_RC_NOT_FOUND
            exit_flag=4
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     MyFile)) RETURN
          END IF
!
!  Inquire about associated time variable.
!
          IF (got_time.and.(Nrec.ge.1)) THEN
            CALL netcdf_inq_var (ng, imodel, Fname,                     &
     &                           ncid = ncid,                           &
     &                           MyVarName = TRIM(T_name),              &
     &                           VarID = Tid,                           &
     &                           nVarDim = ntdim,                       &
     &                           nVarAtt = ntatt)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     MyFile)) THEN
              rc=ESMF_RC_FILE_READ
              RETURN
            END IF
            Export(ifield)%Tname=TRIM(T_name)
            Export(ifield)%Tid=Tid
            Export(ifield)%Nrec=Nrec
!
!  Check associated time variable attributes. The internal processing of
!  all fields requires time in day units
!
            DO i=1,ntatt
              SELECT CASE (TRIM(var_Aname(i)))
                CASE ('cycle_length')
                  Lcycle=.TRUE.                      ! time cycling data
                  IF (var_Afloat(i).gt.0.0_r8) THEN
                    Clength=var_Afloat(i)
                  ELSE IF (var_Aint(i).gt.0) THEN    ! no CF compliance
                    Clength=REAL(var_Aint(i),r8)     ! attribute is an
                  ELSE                               ! integer
                    IF (localPET.eq.0)                                  &
                      WRITE (dataout,50) TRIM(var_Aname(i)),            &
     &                                   TRIM(T_name)
                    rc=ESMF_RC_VAL_WRONG
                    exit_flag=2
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) RETURN
                  END IF
                  Export(ifield)%Lcycle=Lcycle
                  Export(ifield)%Clength=Clength
                CASE ('units')
                  Tunits=TRIM(var_Achar(i))
                  IF (Tunits(1:6).eq.'second') THEN
                    Tscale=1.0_dp/86400.0_dp         ! seconds to days
                  ELSE
                    Tscale=1.0_dp                    ! day units
                  END IF
                  Export(ifield)%Tunits=TRIM(Tunits)
                  Export(ifield)%Tscale=Tscale
              END SELECT
            END DO
!
!  Read associated time variable.
!
            IF (.not.allocated(TimeValue)) THEN
              allocate ( TimeValue(Nrec) )
            END IF
            CALL netcdf_get_time (ng, imodel, Fname, T_name,            &
     &                            ReferenceDateNumber,                  &
     &                            TimeValue,                            &
     &                            ncid = ncid,                          &
     &                            start = (/1/),                        &
     &                            total = (/Nrec/),                     &
     &                            min_val = Tmin,                       &
     &                            max_val = Tmax)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     MyFile)) THEN
              rc=ESMF_RC_FILE_READ
              RETURN
            END IF
!
!  Scale time variable to days.  Determine the minimum and maximum time
!  values available.
!
            DO i=1,Nrec
              TimeValue(i)=TimeValue(i)*Tscale
            END DO
            Tmin=Tmin*Tscale
            Tmax=Tmax*Tscale
            Tstr=Tmin                   ! initialize for the case of
            Tend=Tmax                   ! Nrec=1
            IF (Lcycle) THEN
              Tday=MOD(Tcurrent,Clength)
            ELSE
              Tday=Tcurrent
            END IF
!
!  Is the model time inside the data time range? If not, check if the
!  data just has the LOWER- or the UPPER-snapshot interpolant.
!
            IF ((Tmin.le.Tday).and.(Tday.le.Tmax)) THEN
              Linside=.TRUE.
            ELSE IF (Tday.ge.Tmax) THEN
              LowerBound=.TRUE.
            ELSE IF (Tday.le.Tmin) THEN
              UpperBound=.TRUE.
            END IF
!
!  If processing field split in several files, find UPPER time-snapshot
!  and its associated time record (Trec).
!
            IF (Lmulti) THEN
              DO i=1,Nrec
                IF (TimeValue(i).gt.Tday) THEN
                  Trec=i-1               ! one is added when processing
                  Tend=TimeValue(i)
                  EXIT
                END IF
              END DO
!
!  If not processing a multi-file field or initialization, find LOWER
!  time-snapshot and its associated time record (Trec).  Notice that the
!  conditional below uses (Tstr .le. Tday .le. Tend) when bracketing the
!  LOWER time-snapshot instead the usual (Tstr .le. Tday .lt. Tend). It
!  is to transition smoothly from the end of a multifile to the next
!  during initialization when Lmulti is still false.  The logic is
!  tricky.  We need to check if Tday is equal to the time of the last
!  record in the file to compute Trec correctly if "TIME_INTERP" is
!  activated.
!
            ELSE
              IF (Linside) THEN
                tstart=Tmin
                DO i=2,Nrec
                  IF ((tstart.le.Tday).and.(Tday.le.TimeValue(i))) THEN
                    IF ((Tday.eq.TimeValue(i)).and.(i.ne.Nrec)) THEN
                      Tstr=TimeValue(i)
                      Trec=i             ! one is added when processing
                    ELSE
                      Tstr=tstart
                      Trec=i-1           ! one is added when processing
                    END IF
                    EXIT
                  END IF
                  tstart=TimeValue(i)
                END DO
              ELSE
                Tstr=Tmax        ! LowerBound for next multifile or
                Trec=Nrec        ! time cycling
              END IF
            END IF
!
!  If processing a multi-file field, set LOWER time-snapshot. It
!  is the last value from previous file. Otherwise, set UPPER
!  time-snapshot.
!
            IF (Lmulti) THEN
              Tstr=Export(ifield)%Tmax    ! Tmax from previous file
            ELSE
              IF (Lcycle.and.(Trec.eq.Nrec)) THEN
                Tend=Tmin
              ELSE
                i=MIN(Nrec,Trec+1)
                Tend=TimeValue(i)
              END IF
            END IF
            IF (allocated(TimeValue)) THEN
              deallocate (TimeValue)
            END IF
            Export(ifield)%Tmin=Tmin
            Export(ifield)%Tmax=Tmax
            Export(ifield)%Tstr=Tstr
            Export(ifield)%Tend=Tend
!
!  If not cycling, stop execution if there is not field data
!  available for current model time. Avoid check on tidal data
!  since time is in terms of frequencies.
!
            IF (.not.Lcycle.and.(Nrec.gt.1)) THEN
              IF (Lmulti) THEN
                IF (Tcurrent.gt.Tmax) THEN
                  IF (localPET.eq.0) WRITE (dataout,60) TRIM(T_name),   &
     &                                                  Tmax, Tcurrent
                  rc=ESMF_RC_VAL_WRONG
                  exit_flag=2
                  IF (FoundError(exit_flag, NoError, __LINE__,          &
     &                           MyFile)) RETURN
                ELSE IF (Export(ifield)%LastRec) THEN
                  IF (Tmin.lt.Tcurrent) THEN
                    IF (localPET.eq.0) THEN
                      WRITE (dataout,70)                                &
     &                  'Upper snapshot time for multi-file variable:', &
     &                                  TRIM(T_name),                   &
     &                                  TRIM(ncvname),                  &
     &                  'is less than current model time.',             &
     &                                  'Tmin = ', Tmin, Tcurrent
                    END IF
                    rc=ESMF_RC_VAL_WRONG
                    exit_flag=2
                    IF (FoundError(exit_flag, NoError, __LINE__,        &
     &                             MyFile)) RETURN
                  END IF
                END IF
              ELSE
                IF (.not.UpperBound.and.(Tcurrent.lt.Tmin)) THEN
                  IF (localPET.eq.0) THEN
                     WRITE (dataout,70)                                 &
     &                 'starting time for variable:',                   &
     &                                  TRIM(T_name),                   &
     &                                  TRIM(ncvname),                  &
     &                 'is greater than current model time.',           &
     &                                  'Tmin = ', Tmin, Tcurrent
                  END IF
                  rc=ESMF_RC_VAL_WRONG
                  exit_flag=2
                  IF (FoundError(exit_flag, NoError, __LINE__,          &
     &                           MyFile)) RETURN
                END IF
              END IF
            END IF
          END IF
          Export(ifield)%ncfile=TRIM(Fname)    ! sucess query, need to
          EXIT QUERY                           ! exit and keep last
        ELSE                                   ! value for "ifile"
          IF (.not.Lmulti) THEN
            ncfile=Fname                       ! need for error report
          END IF
        END IF PROBE
!
!  Close input NetCDF file if opened during the query.  Files opened
!  outside the query loop remain open.  This is done to avoid opening
!  too many files.
!
        IF (CloseFile) THEN
          CALL netcdf_close (ng, imodel, ncid, Fname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   MyFile)) THEN
            rc=ESMF_RC_FILE_CLOSE
            RETURN
          END IF
        END IF
      END DO QUERY
!
!-----------------------------------------------------------------------
!  Terminate execution requested variables are not found.
!-----------------------------------------------------------------------
!
      IF (.not.got_var) THEN
        lstr=LEN_TRIM(ncfile)
        IF (localPET.eq.0) THEN
          WRITE (dataout,80) TRIM(ncvname), 'file:'
          IF (lstr.gt.0) THEN
            WRITE (dataout,'(15x,a)') TRIM(ncfile)
          ELSE
            WRITE (dataout,'(15x,a,a)') 'file name is blank, ',         &
     &                                  'cannot be determined.'
          END IF
        END IF
        rc=ESMF_RC_NOT_FOUND
        exit_flag=2
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) RETURN
      END IF
!
      IF (.not.got_time) THEN
        lstr=LEN_TRIM(ncfile)
        IF (localPET.eq.0) THEN
          WRITE (dataout,80) TRIM(T_name), 'file:'
          IF (lstr.gt.0) THEN
             WRITE (dataout,'(15x,a)') TRIM(ncfile)
          ELSE
            WRITE (dataout,'(15x,a,a)') 'file name is blank, ',         &
     &                                  'cannot be determined.'
          END IF
        END IF
        rc=ESMF_RC_NOT_FOUND
        exit_flag=2
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  If appropriate, open input NetCDF file for reading.  Notice that the
!  "ifile" is correct here because of EXIT QUERY command.
!  HGA: Why is there a need to open the appropiate file?  It is still
!       open when the variable is found and processed.
!-----------------------------------------------------------------------
!
      IF (IFS(ifile)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, imodel, ncfile, 0, ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) THEN
          rc=ESMF_RC_FILE_OPEN
          IF (localPET.eq.0) WRITE (dataout,60) TRIM(ncfile)
          RETURN
        END IF
        IFS(ifile)%ncid=ncid
      END IF
!
!-----------------------------------------------------------------------
!  The strategy here is to create a local, monotonically increasing
!  time variable so the interpolation between snapshots is trivial
!  when time cycling data. Notice that a one is substracted to time
!  record counter "Trec" to avoid doing special case at initialization.
!-----------------------------------------------------------------------
!
      IF (.not.Lmulti) THEN
        IF (Lcycle) THEN                 ! time cycling of data, like
          IF (Trec.eq.Nrec) THEN         ! for perpectual annual forcing
            IF (Tcurrent.lt.Tmax) THEN
              Tmono=Tstr-Clength
            ELSE
              Tmono=Tcurrent+(Tstr-Clength)
              IF (Tstr.eq.Tmax) THEN
                Tmono=Tmono+(Tmin-MOD(Tcurrent+Tmin,Clength))
              ELSE
                Tmono=Tmono+(Tstr-MOD(Tcurrent+Tstr,Clength))
              END IF
            END IF
          ELSE
            IF (Tcurrent.gt.Clength) THEN
              Tmono=Tcurrent-MOD(tdays(ng)-Tstr,Clength)
            ELSE
              Tmono=Tstr
            END IF
          END IF
        ELSE
          Tmono=Tstr
        END IF
        Tindex=2
        Trec=Trec-1
        Export(ifield)%Tindex=Tindex
        Export(ifield)%Trec=Trec
        Export(ifield)%Tmono=Tmono
        Export(ifield)%Vtime(Tindex)=Tstr
      ELSE
        Tindex=2
        Export(ifield)%Tindex=Tindex
        Export(ifield)%Trec=Trec
        Export(ifield)%Vtime(Tindex)=Tstr
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_inquiry',            &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
   5  FORMAT (' DATA_inquiry     - inquiring NetCDF variable ''',a,     &
     &        ''' in input file(s) ...')
  35  FORMAT (20x,'found in:  ',a)
  10  FORMAT (/,' DATA_INQUIRY - out of range multi-files counter ',    &
     &        'for variable: ',a,/,16x,'Fcount = ',i2.2,                &
     &        ',   Expected range: 1 - ',i2.2)
  20  FORMAT (/,' DATA_INQUIRY - unable to assign file counter, ',      &
     &        'Fcount = ',i4,/,15x,'while processing structure: ',a,    &
     &        /,16x,'and variable; ',a)
  30  FORMAT (/,' DATA_INQUIRY - unable to open input NetCDF file: ',a)
  40  FORMAT (/,' DATA_INQUIRY - unable to find dimension ',a,          &
     &        /,16x,'for variable: ',a,/,16x,'in file: ',a,             &
     &        /,16x,'file is not CF compliant...')
  50  FORMAT (/,' DATA_INQUIRY - unable to get value for attribute: ',  &
     &        a,/,16x,'in variable: ',a,                                &
     &        /,16x,'This attribute value is expected to be of',        &
     &        /,16x,'the same external type as the variable.')
  60  FORMAT (/,' DATA_INQUIRY - ending time for multi-file variable: ',&
     &        a,/,16x,'is less than current model time. ',              &
     &        /,16x,'Tmax = ',f15.4,2x,'Tcurrent = ',f15.4)
  70  FORMAT (/,' DATA_INQUIRY - ',a,1x,a,2x,'(',a,')',/,16x,a,         &
     &        /,16x,a,f15.4,2x,'Tcurrent = ',f15.4)
  80  FORMAT (/,' DATA_INQUIRY - unable to find requested variable: ',  &
     &        a,/,16x,'in ',a)

      RETURN
      END SUBROUTINE DATA_inquiry
!
      SUBROUTINE DATA_ncvarcoords (ifield, FieldName,                   &
     &                             Nfields, Export,                     &
     &                             localPET, rc)
!
!=======================================================================
!                                                                      !
!  This routine reads the spatial locations of DATA model export       !
!  variable from associated input NetCDF file. If available, it        !
!  also reads the land/se mask.                                        !
!                                                                      !
!  It assumes that the NetCDF variable has the attribute "coordinates",!
!  as specified by CF compatability standard.                          !
!                                                                      !
!  For example, in CDL syntax:                                         !
!                                                                      !
!       float my_var(time, lat, lon) ;                                 !
!               my_var:long_name = "my variable long name" ;           !
!               my_var:units = "my variable units" ;                   !
!               my_var:coordinates = "lon lat time" ;                  !
!               my_var:time = "time" ;                                 !
!                                                                      !
!  The following "coordinates" attribute is also allowed:              !
!                                                                      !
!               my_var:coordinates = "lon lat" ;                       !
!                                                                      !
!  That is, the time variable "time" is missing in the "coordinates"   !
!  attribute.                                                          !
!                                                                      !
!  Notice that the associated coordinate names "lon" and "lat" are     !
!  separated by a single blank space.  Both "lon" and "lat" can be     !
!  1D or 2D arrays. If 1D array, the positions are rectangular and     !
!  and full 2D arrays are filled with the same values.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ifield     Field index in Export structure to process (integer)  !
!     FieldName  DATA model field short name to process (string)       !
!     Export     DATA model export structure, TYPE(ESM_Data)           !
!     Nfields    Number of fields in Export structure (interger)       !
!     localPET   Local Persistent Execution Thread (integer)           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Export     DATA model export structure:                          !
!                  Export(ifield)%lon       longitude 2D locations     !
!                  Export(ifield)%lat       latitude  2D locations     !
!                  Export(ifield)%mask      land/sea mask              !
!                  Export(ifield)%LonMin    minimum longitude          !
!                  Export(ifield)%LonMax    maximum longitude          !
!                  Export(ifield)%LatMin    minimum latitude           !
!                  Export(ifield)%LatMax    maximum latitude           !
!     rc         Return code flag (integer)                            !
!                                                                      !
!  WARNING:                                                            !
!                                                                      !
!  This routine uses ROMS NetCDF managing framework.                   !
!                                                                      !
!=======================================================================
!
      USE mod_netcdf
!
      USE mod_iounits,  ONLY : SourceFile
      USE mod_scalars,  ONLY : NoError, exit_flag, spval
!
      USE strings_mod, ONLY : FoundError, lowercase
!
!  Imported variable declarations.
!
      integer, intent(in)  :: ifield, Nfields
      integer, intent(in)  :: localPET
      integer, intent(out) :: rc
!
      TYPE(ESM_Data), intent(inout) :: Export(Nfields)
!
      character (len=*), intent(in) :: FieldName
!
!  Local variable declarations
!
      logical :: Lcoord, got_lon, got_lat, got_mask
!
      integer :: ng, model
      integer :: Imax, Jmax, i, j
      integer :: Nvdim, ncid, ncvid, nlatatt, nlatdim, nlonatt, nlondim
      integer :: nmaskdim, nmaskatt
      integer :: ROMScomm
!
      real(r8) :: Lon_Min, Lon_Max, Lat_Min, Lat_Max
!
      real(r8), allocatable :: LonWrk(:)
      real(r8), allocatable :: LatWrk(:)
!
      character (len=5  ) :: lstr
      character (len=20 ) :: Dname(2), LonName, LatName, MaskName
      character (len=100) :: ncvname
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", DATA_ncvarcoords"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering DATA_ncvarcoords',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
      rc=ESMF_SUCCESS
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in variable spatial coordinates.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      ng=1                                  ! Needed for ROMS interface
      model=1                               ! Needed for ROMS interface
      got_lon=.FALSE.
      got_lat=.FALSE.
      got_mask=.FALSE.
!
      IF (localPET.eq.0) THEN
        WRITE (dataout,10) TRIM(Export(ifield)%Vname),                  &
     &                     TRIM(Export(ifield)%ncfile)
      END IF
!
!  Load information variables that were set when calling "DATA_inquiry".
!
      ncname=TRIM(Export(ifield)%ncfile)
      ncvname=TRIM(Export(ifield)%Vname)
      Dname(1)=TRIM(Export(ifield)%Dname(1))
      Dname(2)=TRIM(Export(ifield)%Dname(2))
      Lcoord=Export(ifield)%Lcoord
      Nvdim=Export(ifield)%Nvdim

      IF (Lcoord.and.(Nvdim.ge.2)) THEN
        got_lon=.TRUE.
        got_lat=.TRUE.
        LonName=TRIM(Export(ifield)%Vcoord(1))
        LatName=TRIM(Export(ifield)%Vcoord(2))
        Imax=Export(ifield)%Vsize(1)
        Jmax=Export(ifield)%Vsize(2)
      END IF
!
!  If applicable, open input NetCDF for reading.
!
      IF (Export(ifield)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) THEN
          WRITE (dataout,20) TRIM(ncname)
          RETURN
        END IF
        Export(ifield)%ncid=ncid
      ELSE
        ncid=Export(ifield)%ncid            ! already open
      END IF
!
!  Inquire NetCDF file variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               MyFile)) THEN
        rc=ESMF_RC_FILE_READ
        RETURN
      END IF
!
!  Check NetCDF variable and look for longitude and latitude variables
!  guess names.
!
      IF (.not.Lcoord) THEN
        DO i=1,n_var
          SELECT CASE (TRIM(lowercase(var_name(i))))
             CASE ('longitude', 'lon', 'lon_rho', 'lon_u', 'lon_v')
               IF (.not.got_lon) THEN
                 CALL netcdf_inq_var (ng, model, ncname,                &
     &                                ncid = ncid,                      &
     &                                MyVarName = TRIM(var_name(i)),    &
     &                                nVarDim = nlondim,                &
     &                                nVarAtt = nlonatt)
                 IF (FoundError(exit_flag, NoError, __LINE__,           &
     &                          MyFile)) THEN
                   rc=ESMF_RC_FILE_READ
                   RETURN
                 END IF
                 IF ((TRIM(var_Dname(1)).eq.Dname(1)).and.              &
     &               (TRIM(var_Dname(2)).eq.Dname(2))) THEN
                   LonName=TRIM(var_name(i))
                   got_lon=.TRUE.
                 END IF
               END IF
             CASE ('latitude',  'lat', 'lat_rho', 'lat_u', 'lat_v')
               IF (.not.got_lat) THEN
                 CALL netcdf_inq_var (ng, model, ncname,                &
     &                                ncid = ncid,                      &
     &                                MyVarName = TRIM(var_name(i)),    &
     &                                nVarDim = nlatdim,                &
     &                                nVarAtt = nlatatt)
                 IF (FoundError(exit_flag, NoError, __LINE__,           &
     &                          MyFile)) THEN
                   rc=ESMF_RC_FILE_READ
                   RETURN
                 END IF
                 IF ((TRIM(var_Dname(1)).eq.Dname(1)).and.              &
     &               (TRIM(var_Dname(2)).eq.Dname(2))) THEN
                   LatName=TRIM(var_name(i))
                   got_lat=.TRUE.
                 END IF
               END IF
          END SELECT
        END DO
!
        IF (.not.(got_lon.or.got_lat)) THEN
          rc=ESMF_RC_CANNOT_GET
          exit_flag=2
          WRITE (lstr,'(i5)') __LINE__
          IF (localPET.eq.0) THEN
            WRITE (dataout,10) TRIM(ncvname), TRIM(ncname),             &
     &                         exit_flag, ADJUSTL(TRIM(lstr)),          &
     &                         MyFile
          END IF
          RETURN
        END IF
      END IF
!
!  Check for land/sea mask using various guess names.
!
      DO i=1,n_var
        SELECT CASE (TRIM(lowercase(var_name(i))))
          CASE ('mask', 'mask_rho', 'mask_u', 'mask_v')
            IF (.not.got_mask) THEN
              CALL netcdf_inq_var (ng, model, ncname,                   &
     &                             ncid = ncid,                         &
     &                             MyVarName = TRIM(var_name(i)),       &
     &                             nVarDim = nmaskdim,                  &
     &                             nVarAtt = nmaskatt)
              IF (FoundError(exit_flag, NoError, __LINE__,              &
     &                       MyFile)) THEN
                rc=ESMF_RC_FILE_READ
                RETURN
              END IF
              IF ((TRIM(var_Dname(1)).eq.Dname(1)).and.                 &
     &            (TRIM(var_Dname(2)).eq.Dname(2))) THEN
                MaskName=TRIM(var_name(i))
                got_mask=.TRUE.
              END IF
            END IF
        END SELECT
      END DO
!
!  Read in longitude coordinate.
!
      IF (.not.allocated(Export(ifield)%lon)) THEN
        allocate ( Export(ifield)%lon(Imax,Jmax) )
      END IF
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = ncid,                                 &
     &                     MyVarName = TRIM(LonName),                   &
     &                     VarID = ncvid,                               &
     &                     nVarDim = nlondim,                           &
     &                     nVarAtt = nlonatt)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               MyFile)) THEN
        rc=ESMF_RC_FILE_READ
        RETURN
      END IF
!
      IF (nlondim.eq.1) THEN
        IF (.not.allocated(LonWrk)) THEN
          allocate ( LonWrk(Imax) )
        END IF
        CALL netcdf_get_fvar (ng, model, ncname, LonName,               &
     &                        LonWrk,                                   &
     &                        ncid = ncid,                              &
     &                        start = (/1/),                            &
     &                        total = (/Imax/))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) RETURN

        DO j=1,Jmax
          DO i=1,Imax
            Export(ifield)%lon(i,j)=LonWrk(i)
          END DO
        END DO
        deallocate (LonWrk)
      ELSE
        CALL netcdf_get_fvar (ng, model, ncname, LonName,               &
     &                        Export(ifield)%lon,                       &
     &                        ncid = ncid,                              &
     &                        start = (/1,1/),                          &
     &                        total = (/Imax,Jmax/))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) RETURN
      END IF
!
!  Read in latitute coordinate.
!
      IF (.not.allocated(Export(ifield)%lat)) THEN
        allocate ( Export(ifield)%lat(Imax,Jmax) )
      END IF
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = ncid,                                 &
     &                     MyVarName = TRIM(LatName),                   &
     &                     VarID = ncvid,                               &
     &                     nVarDim = nlatdim,                           &
     &                     nVarAtt = nlatatt)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               MyFile)) THEN
        rc=ESMF_RC_FILE_READ
        RETURN
      END IF
!
      IF (nlatdim.eq.1) THEN
        IF (.not.allocated(LatWrk)) THEN
          allocate ( LatWrk(Jmax) )
        END IF
        CALL netcdf_get_fvar (ng, model, ncname, LatName,               &
     &                        LatWrk,                                   &
     &                        ncid = ncid,                              &
     &                        start = (/1/),                            &
     &                        total = (/Jmax/))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) RETURN

        DO j=1,Jmax
          DO i=1,Imax
            Export(ifield)%lat(i,j)=LatWrk(j)
          END DO
        END DO
        deallocate (LatWrk)
      ELSE
        CALL netcdf_get_fvar (ng, model, ncname, LatName,               &
     &                        Export(ifield)%lat,                       &
     &                        ncid = ncid,                              &
     &                        start = (/1,1/),                          &
     &                        total = (/Imax,Jmax/))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) RETURN
      END IF
!
!  Determine data minimum and maximum longitude/latitude.
!
      Lon_Min= spval
      Lon_Max=-spval
      Lat_Min= spval
      Lat_Max=-spval
      DO j=1,Jmax
        DO i=1,Imax
          Lon_Min=MIN(Lon_Min, Export(ifield)%lon(i,j))
          Lon_Max=MAX(Lon_Max, Export(ifield)%lon(i,j))
          Lat_Min=MIN(Lat_Min, Export(ifield)%lat(i,j))
          Lat_Max=MAX(Lat_Max, Export(ifield)%lat(i,j))
        END DO
      END DO
      Export(ifield)%LonMin=Lon_Min
      Export(ifield)%LonMax=Lon_Max
      Export(ifield)%LatMin=Lat_Min
      Export(ifield)%LatMax=Lat_Max
!
!  Read or set land/sea mask.
!
      IF (.not.allocated(Export(ifield)%mask)) THEN
        allocate ( Export(ifield)%mask(Imax,Jmax) )
      END IF
      Export(ifield)%Lmask=got_mask
!
      IF (got_mask) THEN
        CALL netcdf_get_fvar (ng, model, ncname, MaskName,              &
     &                        Export(ifield)%mask,                      &
     &                        ncid = ncid,                              &
     &                        start = (/1,1/),                          &
     &                        total = (/Imax,Jmax/))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 MyFile)) RETURN
      ELSE
        DO j=1,Jmax
          DO i=1,Imax
            Export(ifield)%mask(i,j)=1.0_r8
          END DO
        END DO
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  DATA_ncvarcoords',        &
     &                          ', PET', PETrank
        CALL my_flush (trac)
      END IF
!
 10   FORMAT (' DATA_ncvarcoords - setting spatial coordinates for',    &
     &        ' NetCDF variable ''',a,'''',/,20x,'from file: ',a)
 20   FORMAT (/,' DATA_ncvarcoords - Cannot find "coordinates" ',       &
     &        'attribute for variable:',2x,a,/,20x,'in file:',2x,a,/,   &
     &     /,20x,'This attribute is needed to interpolate input data',  &
     &     /,20x,'to model grid. Following CF compliance, we need:',/,  &
     &     /,20x,'float my_var(time, lat, lon) ;',                      &
     &     /,20x,'      my_var:long_name = "my variable long name" ;',  &
     &     /,20x,'      my_var:units = "my variable units" ;',          &
     &     /,20x,'      my_var:coordinates = "lon lat my_var_time" ;',  &
     &     /,20x,'      my_var:time = "my_var_time" ;',/,               &
     &        /, ' Found Error: ', i2.2, t20, 'Line: ',a,               &
     &        t35, 'Source: ', a)

      RETURN
      END SUBROUTINE DATA_ncvarcoords

#endif
      END MODULE esmf_data_mod
