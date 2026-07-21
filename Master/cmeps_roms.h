#include "cppdefs.h"
      MODULE cmeps_roms_mod

#if defined MODEL_COUPLING && defined ESMF_LIB && defined CMEPS
!
!git $Id$
!>
!! \brief   **ROMS** ESMF/NUOPC Cap file for **CMEPS**
!!
!! \details A NUOPC-compliant software layer on top of ROMS to couple to
!!          any Earth System Model (ESM) component through the Community
!!          Mediator for Earth Prediction Systems (CMEPS).
!!
!! \author  Hernan G. Arango (Rutgers University)
!! \date    January 2022
!!
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.md                          Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This ESMF/NUOPC module sets ROMS as an ocean gridded component in   !
!  the cdeps/cmeps system:                                             !
!                                                                      !
!                                                                      !
!    ROMS_SetServices        Sets ROMS component shared-object entry   !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                                                                      !
!    ROMS_Create             Allocates module structures and process   !
!                            configuration from input YAML file.       !
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
!                            grid arrays, and land/sea mask if any.    !
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
! 'CMEPS': Community Mediator for Earth Prediction Systems             !
!            'https://escomp.github.io/CMEPS'                          !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            'https://www.earthsystemcog.org/projects/esmf'            !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            'https://www.earthsystemcog.org/projects/nuopc'           !
!                                                                      !
!  ROMS:   Regional Ocean Modeling System                              !
!            'https://www.myroms.org'                                  !
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
!-----------------------------------------------------------------------
!  ROMS module association: parameters, variables, derived-type objects.
!-----------------------------------------------------------------------
!
      USE roms_kernel_mod,  ONLY : ROMS_initialize,                     &
     &                             ROMS_run,                            &
     &                             ROMS_finalize
!
      USE bc_2d_mod,        ONLY : bc_r2d_tile
      USE dateclock_mod,    ONLY : ROMS_clock, caldate, get_date,       &
     &                             time_string
      USE exchange_2d_mod,  ONLY : exchange_r2d_tile,                   &
     &                             exchange_u2d_tile,                   &
     &                             exchange_v2d_tile
      USE get_metadata_mod, ONLY : CouplingField,                       &
     &                             cmeps_metadata,                      &
     &                             metadata_has
      USE mod_kinds,        ONLY : dp, i4b, i8b, r4, r8
      USE mod_forces,       ONLY : FORCES
      USE mod_grid,         ONLY : GRID
      USE mod_iounits,      ONLY : Iname, SourceFile, stdout
      USE mod_mixing,       ONLY : MIXING
      USE mod_ncparam,      ONLY : Iinfo,  idLdwn, idLrad, idPair,      &
     &                             idPBLh, idQair, idrain, idSrad,      &
     &                             idTair, idTsur, idUair, idUsms,      &
     &                             idVair, idVsms
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
      USE stdout_mod,       ONLY : Set_StdOutUnit, stdout_unit
      USE strings_mod,      ONLY : FoundError, assign_string, lowercase
      USE yaml_parser_mod,  ONLY : yaml_AssignString,                   &
     &                             yaml_Error,                          &
     &                             yaml_Svec,                           &
     &                             yaml_get,                            &
     &                             yaml_initialize,                     &
     &                             yaml_tree
!
!-----------------------------------------------------------------------
!  Derived-type structures.
!-----------------------------------------------------------------------
!
      implicit none
!
!  ESM coupling time managing variables and ESMF objects.
!
      TYPE, PRIVATE :: ESM_Clock
        logical :: Restarted

        integer (i8b) :: AdvanceCount               ! advance counter

        real (dp) :: Current_Time                   ! seconds
        real (dp) :: Time_Reference                 ! seconds
        real (dp) :: Time_Restart                   ! seconds
        real (dp) :: Time_Start                     ! seconds
        real (dp) :: Time_Stop                      ! seconds
        real (dp) :: Time_Step                      ! seconds
!
        character (len=22) :: Name
        character (len=22) :: CalendarString        ! 360_day, gregorian
        character (len=22) :: Time_ReferenceString
        character (len=22) :: Time_RestartString
        character (len=22) :: Time_StartString
        character (len=22) :: Time_StopString
!
        TYPE (ESMF_Calendar)       :: Calendar
        TYPE (ESMF_Clock)          :: Clock
        TYPE (ESMF_Direction_flag) :: Direction
        TYPE (ESMF_Time)           :: CurrentTime
        TYPE (ESMF_Time)           :: ReferenceTime
        TYPE (ESMF_Time)           :: RestartTime
        TYPE (ESMF_Time)           :: StartTime
        TYPE (ESMF_Time)           :: StopTime
        TYPE (ESMF_TimeInterval)   :: TimeStep
      END TYPE ESM_Clock
!
      TYPE (ESM_Clock), allocatable, target :: ClockInfo(:)
!
! ESM coupled state sets. If appropriate, it includes the logic for
! connecting nested grids.
!
! In other ESM components, the import and export states are defined and
! allocated in the MODELS structure as vectors (see below).
!
      TYPE, PRIVATE :: ESM_CplSet
        logical, allocatable :: LinkedGrid(:,:)       ! connected grid
!
        character (len=100), allocatable :: SetLabel(:) ! set label
        character (len=100), allocatable :: ExpLabel(:) ! export label
        character (len=100), allocatable :: ImpLabel(:) ! import label
      END TYPE ESM_CplSet
!
      TYPE (ESM_CplSet), allocatable, target :: COUPLED(:)
!
!  Import and export fields metadata information.
!
      TYPE, PRIVATE :: ESM_Field
        logical :: connected                  ! connected to coupler
        logical :: debug_write                ! write exchanged field
!
        integer :: gtype                      ! field grid mesh type
        integer :: itype                      ! field interpolation flag
        integer :: Tindex                     ! rolling two-time indices
!
        character (len=:), allocatable :: short_name    ! short name
        character (len=:), allocatable :: standard_name ! standard name
        character (len=:), allocatable :: long_name     ! long name
        character (len=:), allocatable :: dst_gtype     ! DST grid type
        character (len=:), allocatable :: dst_units     ! DST units
        character (len=:), allocatable :: src_gtype     ! SRC grid type
        character (len=:), allocatable :: src_units     ! SRC units
        character (len=:), allocatable :: nc_vname      ! DATA Vname
        character (len=:), allocatable :: nc_tname      ! DATA Tname
        character (len=:), allocatable :: map_norm      ! mapping norm
        character (len=:), allocatable :: map_type      ! regrid method
!
        character (len=22)  :: DateString(2)  ! snapshots date
!
        real (dp) :: scale_factor             ! field scale factor
        real (dp) :: add_offset               ! field add offset value
        real (dp) :: Tmin                     ! DATA time minimum value
        real (dp) :: Tmax                     ! DATA time maximum value
        real (dp) :: Tstr                     ! DATA lower time-snapshot
        real (dp) :: Tend                     ! DATA upper time-snapshot
        real (dp) :: Tintrp(2)                ! interpolation time (day)
        real (dp) :: Vtime(2)                 ! latest two-time values
!
        TYPE (ESMF_RouteHandle) :: rhandle    ! field RouteHandle
      END TYPE ESM_Field
!
!  Import and export fields mesh data.
!
      TYPE, PRIVATE :: ESM_Mesh
        integer :: gid                           ! grid ID
        integer :: gtype                         ! grid mesh type

        integer (i4b), allocatable :: mask(:,:)  ! grid land/sea mask

        real (r8), allocatable :: lon(:,:)       ! grid longitude
        real (r8), allocatable :: lat(:,:)       ! grid latitude
        real (r8), allocatable :: area(:,:)      ! grid area
      END TYPE ESM_Mesh
!
!  Coupled models high-level object, [Nmodels].
!
      TYPE, PRIVATE :: ESM_Model
        logical :: IsActive                       ! active for coupling
!
        integer (i4b) :: LandValue                ! land mask value
        integer (i4b) :: SeaValue                 ! sea  mask value

        integer :: Ngrids                         ! number nested grids

        integer :: ExportCalls                    ! export CALL counter
        integer :: ImportCalls                    ! import CALL counter

        integer :: nPETs                          ! number model PETs
        integer, allocatable :: PETlist(:)        ! component PETs list

        integer, allocatable :: TimeFrac(:,:)     ! driver time fraction
!
        character (len=:), allocatable :: name    ! component name
!
        TYPE (ESMF_Grid),  allocatable :: grid(:)        ! grid object
        TYPE (ESM_Mesh),   allocatable :: mesh(:)        ! mesh
        TYPE (ESM_Field),  allocatable :: ImportField(:) ! import fields
        TYPE (ESM_Field),  allocatable :: ExportField(:) ! export fields
        TYPE (ESMF_State), allocatable :: ImportState(:) ! import state
        TYPE (ESMF_State), allocatable :: ExportState(:) ! export state
      END TYPE ESM_Model
!
      TYPE (ESM_Model), allocatable, target :: MODELS(:)
!
!-----------------------------------------------------------------------
!  Internal module paramters and variables.
!-----------------------------------------------------------------------
!
!  Switch to trace/track run sequence during debugging.  All information
!  is written to Fortan unit trac. For now, use standard output unit.
!
      logical :: ESM_track = .TRUE.  ! trace/track CALL sequence switch
!
!  Number of coupled ESM gridded components and identification index
!
      integer, parameter :: Nmodels = 1
      integer, parameter :: Idriver = 0
      integer, parameter :: Iroms   = 1
!
!  Number of ROMS export and import fields per component.
!
      integer, allocatable :: Nexport(:)
      integer, allocatable :: Nimport(:)
!
!  Model coupling type: [1] Explicit, [otherwise] Semi-Implicit.
!
!  In explicit coupling, exchange fields at the next time-step are
!  defined using known values from the time-step before it. Explicit
!  methods require less computational effort and are accurate for
!  small coupling time-steps.
!
!  In implicit coupling, exchange fields at the next time-step are
!  defined by including values at the next time-step. Implicit methods
!  are stable and allow longer coupling time-steps but are expensier.
!
!  In semi-implicit coupling, ROMS -> ATM is explicit, ATM -> ROMS is
!  implicit.
!
      integer :: CouplingType = 1
!
!  Linked/coupled ROMS nested grid number.
!
      integer :: linked_grid
!
!  Distributed-memory communicator handle for each component, rank of
!  each PET, and PET layout (sequential or concurrent).
!
      integer, allocatable :: ESMcomm(:)
      integer :: PETrank
      character (len=10), allocatable :: PETlayoutOption(:)
!
!  Coupling debugging flag:
!
!    [0] no debugging
!    [1] reports informative messages
!    [2] '1' plus writes grid information in VTK format
!    [3] '2' plus writes exchage fields into NetCDF files
!
      integer :: DebugLevel = 0
!
!  Execution tracing level flag:
!
!    [0] no tracing
!    [1] reports sequence of coupling subroutine calls
!    [2] <1> plus writes voluminous ESMF library tracing
!            information which slowdown performace, and
!            creates large log file
!
      integer :: TraceLevel = 0
!
!  Standard output units.
!
      integer :: cplout = 77          ! coupling driver
      integer :: trac   = 6           ! trace/track CALL sequence unit
!
      character (len=11), parameter :: CouplerLog = 'log.coupler'
!
!  Coupled model staggered grid-cell type indices:
!
!     Arakawa B-grid    Arakawa C-grid
!
!     q --------- q     q --- v --- q
!     |           |     |           |
!     |     c     |     u     c     u
!     |           |     |           |
!     q --------- q     q --- v --- q
!
!  COAMPS, C-grid
!  RegCM,  B-grid
!  ROMS,   C-grid        (c = RHO-point,  q = PSI-point)
!  WRF,    C-grid
!
      integer, parameter :: Inan    = 0     ! unstaggered, cell center
      integer, parameter :: Icenter = 1     ! cell center
      integer, parameter :: Icorner = 2     ! cell corners
      integer, parameter :: Iupoint = 3     ! right and left cell faces
      integer, parameter :: Ivpoint = 4     ! upper and lower cell faces
!
      character (len=6), dimension(0:4) :: GridType =                   &
     &                                     (/ 'N/A   ',                 &
     &                                        'Center',                 &
     &                                        'Corner',                 &
     &                                        'U     ',                 &
     &                                        'V     ' /)
!
!  REGRID interpolation method between source and destination fields.
!
      integer, parameter :: Inone   = 0     ! none
      integer, parameter :: Ibilin  = 1     ! bilinear
      integer, parameter :: Ipatch  = 2     ! high-order patch recovery
      integer, parameter :: Iconsvd = 3     ! 1st-order conservative, D
      integer, parameter :: Iconsvf = 4     ! 1st-order conservative, F
      integer, parameter :: Ifcopy  = 5     ! redist
      integer, parameter :: InStoD  = 6     ! nearest Src 2 Dst
      integer, parameter :: InStoDd = 7     ! nearest Src 2 Dst, consv D
      integer, parameter :: InStoDf = 8     ! nearest Src 2 Dst, consv F
!
!  Standard input filename for each coupled model, [Nmodels].
!
      character (len=256), allocatable :: INPname(:)
!
!  ROMS coupling YAML configuration filename.
!
      character (len=:), allocatable :: CPLname
!
!  Today date string.
!
      character (len=44) :: TodayDateString
!
!  Git URL and Revision.
!
      character (len=80)  :: git_url
      character (len=256) :: git_rev
!
!  ESM strings.
!
      character (len=:), allocatable :: CoupledSet
      character (len=:), allocatable :: ExportStateName
      character (len=:), allocatable :: ImportStateName
!
!  ESM constants.
!
      integer (i4b), parameter :: MAPPED_MASK = 99_i4b
      integer (i4b), parameter :: UNMAPPED_MASK = 98_i4b

      real (dp), parameter :: MISSING_dp = 1.0E20_dp
      real (r4), parameter :: MISSING_r4 = 1.0E20_r4
      real (r8), parameter :: MISSING_r8 = 1.0E20_r8

      real (dp), parameter :: TOL_dp = 0.001E20_dp
      real (r4), parameter :: TOL_r4 = 0.001E20_r4
      real (r8), parameter :: TOL_r8 = 0.001E20_r8
!
!  Specific parameters for cmeps.
!
      integer :: ScalarFieldCount
      integer :: ScalarFieldIdxGridNX
      integer :: ScalarFieldIdxGridNY

      character (len=256) :: ScalarFieldName
!
!-----------------------------------------------------------------------
      PUBLIC  :: ROMS_SetServices
!-----------------------------------------------------------------------
!
      PRIVATE :: ROMS_Create
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
      PRIVATE :: CreateScalarField
      PRIVATE :: SetScalarFieldValues
      PRIVATE :: field_index
      PRIVATE :: report_timestamp
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
      SUBROUTINE ROMS_Create (localPET, PETcount, MyComm, rc)
!
!=======================================================================
!                                                                      !
!  It allocates module structures and process configuration from input !
!  YAML file.                                                          !
!                                                                      !
!=======================================================================
!
      USE mod_strings
!
!  Imported variable declarations.
!
      integer, intent(in ) :: localPET, PETcount, MyComm
      integer, intent(out) :: rc
!
!  Local variable declarations.
!
# ifdef METADATA_REPORT
      logical :: Lreport = .TRUE.                ! dumps YAML dictionary
# else
      logical :: Lreport = .FALSE.
# endif
      logical :: Lexist, MasterPET
!
      integer :: Findex, i, layout, ng
!
      TYPE (CouplingField), allocatable :: Export(:), Import(:)
      TYPE (yaml_Svec), allocatable     :: Estring(:), Istring(:)
      TYPE (yaml_tree)                  :: YML
!
      character (len=240) :: StandardName, ShortName
      character (len=256) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_Create"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '==> Entering ROMS_Create',             &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Create YAML dictionary object (TYPE yaml_tree).
!-----------------------------------------------------------------------
!
      IF (yaml_Error(yaml_initialize(YML, TRIM(CPLname), Lreport),      &
     &               NoError, __LINE__, MyFile)) THEN
        IF (localPET.eq.0) WRITE (cplout,10) TRIM(CPLname)
        rc=ESMF_RC_FILE_READ
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Get ROMS standard input filename.
!-----------------------------------------------------------------------
!
     IF (YML%has('standard_input.OCN_component')) THEN
         IF (FoundError(yaml_get(YML, 'standard_input.OCN_component',    &
     &                           string),                                &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
        Iname=TRIM(string)
      ELSE
        rc=ESMF_RC_NOT_FOUND
        IF (localPET.eq.0) WRITE (cplout,20)                            &
     &                                  'standard_input.OCN_component', &
     &                                       TRIM(CPLname)
        RETURN
      END IF
!
!  Get size of number of nested grids, 'Ngrids' parameter from ROMS
!  standard input file.
!
      MasterPET=localPET.eq.0
      CALL getpar_i (MasterPET, Ngrids, 'Ngrids', TRIM(Iname))
!
!-----------------------------------------------------------------------
!  Get ROMS linked/coupled nested grid number for current application.
!-----------------------------------------------------------------------
!
      IF (YML%has('linked_grid')) THEN
        IF (FoundError(yaml_get(YML, 'linked_grid',                     &
     &                          linked_grid),                           &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
      ELSE
        linked_grid=1
      END IF
!
!-----------------------------------------------------------------------
!  Allocate 'cap' module structure and variables.
!-----------------------------------------------------------------------
!
!  Coupled model(s) high-level objects.
!
      IF (.not.allocated(ClockInfo)) THEN
        allocate ( ClockInfo(0:Nmodels) )         ! TYPE ESM_Clock
      END IF
!
      IF (.not.allocated(COUPLED)) THEN
        allocate ( COUPLED(Nmodels) )             ! TYPE ESM_CplSet
      END IF
!
      IF (.not.allocated(MODELS)) THEN
        allocate ( MODELS(Nmodels) )              ! TYPE ESM_Model
      END IF
!
!  Number of export and import fields per component.
!
      IF (.not.allocated(Nexport)) THEN
        allocate ( Nexport(Nmodels) )
      END IF
!
      IF (.not.allocated(Nimport)) THEN
        allocate ( Nimport(Nmodels) )
      END IF
!
!  Allocate and set several variables that depend on the number of ROMS
!  nested grids.
!
      IF (.not.allocated(ESMcomm)) THEN
        allocate ( ESMcomm(Nmodels) )             ! mpi-communicator
      END IF
      ESMcomm(Iroms)=MyComm
!
      IF (.not.allocated(PETlayoutOption)) THEN
        allocate ( PETlayoutOption(Nmodels) )     ! PET layout
      END IF
!
      IF (.not.allocated(INPname)) THEN
        allocate ( INPname(Nmodels) )             ! standard input files
      END IF
      INPname(Iroms)=TRIM(Iname)
!
      IF (.not.allocated(MODELS(Iroms)%grid)) THEN
        allocate ( MODELS(Iroms)%grid(Ngrids) )
      END IF
!
      IF (.not.allocated(COUPLED(Iroms)%LinkedGrid)) THEN
        allocate ( COUPLED(Iroms)%LinkedGrid(Ngrids,Nmodels) )
        COUPLED(Iroms)%LinkedGrid=.FALSE.
      END IF
      COUPLED(Iroms)%LinkedGrid(linked_grid,:)=.TRUE.
!
      IF (.not.allocated(COUPLED(Iroms)%SetLabel)) THEN
        allocate ( COUPLED(Iroms)%SetLabel(Ngrids) )
      END IF
!
      IF (.not.allocated(COUPLED(Iroms)%ExpLabel)) THEN
        allocate ( COUPLED(Iroms)%ExpLabel(Ngrids) )
      END IF
!
      IF (.not.allocated(COUPLED(Iroms)%ImpLabel)) THEN
        allocate ( COUPLED(Iroms)%ImpLabel(Ngrids) )
      END IF
!
      IF (.not.allocated(MODELS(Iroms)%TimeFrac)) THEN
        allocate ( MODELS(Iroms)%TimeFrac(Ngrids,Nmodels) )
        MODELS(Iroms)%TimeFrac=0
      END IF
      MODELS(Iroms)%TimeFrac(linked_grid,:)=1
!
!  Initialize module various variables.
!
      MODELS(Iroms)%Ngrids=Ngrids
      MODELS(Iroms)%IsActive=.TRUE.
      ClockInfo(Idriver)%Restarted=.FALSE.
      ClockInfo(Iroms  )%Restarted=.FALSE.
!
      IF (FoundError(assign_string(MODELS(Iroms)%name,                  &
     &                             'ROMS'),                             &
     &               NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Process export field(s) metadata from YAML object.
!-----------------------------------------------------------------------
!
!  Allocate ExportState structure. It is still needed if there are no
!  fields to export. It is needed for the 'scalar field' named
!  'cpl_scalars' (see nems.configure) used by cmeps.
!
      IF (.not.allocated(MODELS(Iroms)%ExportState)) THEN
        allocate ( MODELS(Iroms)%ExportState(Ngrids) )
      END IF
!
!  Get export variables short name to process.
!
      IF (YML%has('export_variables')) THEN
        IF (FoundError(yaml_get(YML, 'export_variables',                 &
     &                          Estring),                                &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
        Nexport(Iroms)=SIZE(Estring)
      ELSE
        Nexport(Iroms)=0                           ! no fields to export
      END IF
!
!  Extract export field metadata from dictionary.
!
      IF (Nexport(Iroms).gt.0) THEN
        IF (YML%has('export')) THEN
          CALL cmeps_metadata (YML, TRIM(CPLname), 'export',            &
     &                         Export)
          IF (yaml_Error(exit_flag, NoError, __LINE__, MyFile)) THEN
            rc=ESMF_RC_VAL_WRONG
            RETURN
          END IF
        ELSE
          rc=ESMF_RC_NOT_FOUND
          IF (localPET.eq.0) WRITE (cplout,20) 'export',                &
     &                                         TRIM(CPLname)
          RETURN
        END IF
!
!  Allocate export fields structure (TYPE ESM_Fields).
!
        IF (.not.allocated(MODELS(Iroms)%ExportField)) THEN
          allocate ( MODELS(Iroms)%ExportField(Nexport(Iroms)) )
        END IF
!
!  Load export field(s) metadata.
!
        DO i=1,Nexport(Iroms)
          ShortName=Estring(i)%value
          Findex=metadata_has(Export, TRIM(ShortName))
          IF (Findex.gt.0) THEN
            MODELS(Iroms)%ExportField(i)%connected=                     &
     &                                   Export(Findex)%connected
            MODELS(Iroms)%ExportField(i)%debug_write=                   &
     &                                   Export(Findex)%debug_write
!
            MODELS(Iroms)%ExportField(i)%add_offset=                    &
     &                                   Export(Findex)%add_offset
            MODELS(Iroms)%ExportField(i)%scale_factor=                  &
     &                                   Export(Findex)%scale
!
!                                      field short name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%short_name,     &
     &                     Export(Findex)%short_name),                  &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      field standard name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%standard_name,  &
     &                     Export(Findex)%standard_name),               &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      field descriptive long name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%long_name,      &
     &                     Export(Findex)%long_name),                   &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      field mapping normalization type
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%map_norm,       &
     &                     Export(Findex)%map_norm),                    &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      field reggriding method
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%map_type,       &
     &                     Export(Findex)%map_type),                    &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      destination field grid-cell type
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%dst_gtype,      &
     &                     Export(Findex)%destination_grid),            &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      destination field units
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%dst_units,      &
     &                     Export(Findex)%destination_units),           &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      source field grid-cell type
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%src_gtype,      &
     &                     Export(Findex)%source_grid),                 &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      source field units
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%src_units,      &
     &                     Export(Findex)%source_units),                &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      DATA NetCDF variable name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%nc_vname,       &
     &                     Export(Findex)%data_netcdf_vname),           &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      DATA NetCDF time variable name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ExportField(i)%nc_tname,       &
     &                     Export(Findex)%data_netcdf_tname),           &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!
!  Set grid type flag.
!
            SELECT CASE (lowercase(Export(Findex)%source_grid))
              CASE ('center_cell', 'cell_center', 'center')
                MODELS(Iroms)%ExportField(i)%gtype=Icenter
              CASE ('corner_cell', 'cell_corner', 'corner')
                MODELS(Iroms)%ExportField(i)%gtype=Icorner
              CASE ('left_right_edge', 'right_left_edge', 'edge1')
                MODELS(Iroms)%ExportField(i)%gtype=Iupoint
              CASE ('lower_upper_edge', 'upper_lower_edge', 'edge2')
                MODELS(Iroms)%ExportField(i)%gtype=Ivpoint
              CASE DEFAULT
                MODELS(Iroms)%ExportField(i)%gtype=Icenter
            END SELECT
!
!  Set map type flag.
!
            SELECT CASE (lowercase(Export(Findex)%map_type))
              CASE ('mapbilnr')
                MODELS(Iroms)%ExportField(i)%itype=Ibilin
              CASE ('mappatch')
                MODELS(Iroms)%ExportField(i)%itype=Ipatch
              CASE ('mapconsd')
                MODELS(Iroms)%ExportField(i)%itype=Iconsvd
              CASE ('mapconsf')
                MODELS(Iroms)%ExportField(i)%itype=Iconsvf
              CASE ('mapfcopy')
                MODELS(Iroms)%ExportField(i)%itype=Ifcopy
              CASE ('mapnstod')
                MODELS(Iroms)%ExportField(i)%itype=InStoD
              CASE ('mapnstod_consd')
                MODELS(Iroms)%ExportField(i)%itype=InStoDd
              CASE ('mapnstod_consf')
                MODELS(Iroms)%ExportField(i)%itype=InStoDf
              CASE DEFAULT
                MODELS(Iroms)%ExportField(i)%itype=Inone
            END SELECT
!
!  Check if field exits in NUOPC dictionary.
!
            Lexist=NUOPC_FieldDictionaryHasEntry(                       &
     &                     MODELS(Iroms)%ExportField(i)%standard_name,  &
     &                                        rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  If appropriate, add field to NUOPC dictionary.
!
            IF (.not.Lexist) THEN
              CALL NUOPC_FieldDictionaryAddEntry(                       &
     &                     MODELS(Iroms)%ExportField(i)%standard_name,  &
     &                                         canonicalUnits =         &
     &                     MODELS(Iroms)%ExportField(i)%src_units,      &
     &                                         rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
          ELSE
            IF (localPET.eq.0) THEN
              WRITE (cplout,30) 'export field short_name: ',            &
     &                          TRIM(ShortName), TRIM(CPLname)
            END IF
            rc=ESMF_RC_NOT_FOUND
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
!  Process import field(s) metadata from YAML object.
!-----------------------------------------------------------------------
!
!  Get import variables short name to process.
!
      IF (YML%has('import_variables')) THEN
        IF (FoundError(yaml_get(YML, 'import_variables',                 &
     &                          Istring),                                &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
        Nimport(Iroms)=SIZE(Istring)
      ELSE
        Nimport(Iroms)=0                           ! no fields to import
      END IF

!  Extract ROMS import field(s) metadata for YML dictionary.
!
      IF (Nimport(Iroms).gt.0) THEN
        IF (YML%has('import')) THEN
          CALL cmeps_metadata (YML, TRIM(CPLname), 'import',            &
     &                         Import)
          IF (yaml_Error(exit_flag, NoError, __LINE__, MyFile)) THEN
            rc=ESMF_RC_VAL_WRONG
            RETURN
          END IF
        ELSE
          rc=ESMF_RC_NOT_FOUND
          IF (localPET.eq.0) WRITE (cplout,20) 'import',                &
     &                                         TRIM(CPLname)
          RETURN
        END IF
!
!  Allocate import fields structure (TYPE ESM_Fields).
!
        IF (.not.allocated(MODELS(Iroms)%ImportField)) THEN
          allocate ( MODELS(Iroms)%ImportField(Nimport(Iroms)) )
        END IF
!
        IF (.not.allocated(MODELS(Iroms)%ImportState)) THEN
          allocate ( MODELS(Iroms)%ImportState(Ngrids) )
        END IF
!
!  Load import field(s) metadata.
!
        DO i=1,Nimport(Iroms)
          ShortName=Istring(i)%value
          Findex=metadata_has(Import, TRIM(ShortName))
          IF (Findex.gt.0) THEN
            MODELS(Iroms)%ImportField(i)%connected=                     &
     &                                   Import(Findex)%connected
            MODELS(Iroms)%ImportField(i)%debug_write=                   &
     &                                   Import(Findex)%debug_write
!
            MODELS(Iroms)%ImportField(i)%add_offset=                    &
     &                                   Import(Findex)%add_offset
            MODELS(Iroms)%ImportField(i)%scale_factor=                  &
     &                                   Import(Findex)%scale
!
!                                      field short name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%short_name,     &
     &                     Import(Findex)%short_name),                  &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      field standard name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%standard_name,  &
     &                     Import(Findex)%standard_name),               &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      field descriptive long name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%long_name,      &
     &                     Import(Findex)%long_name),                   &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      field mapping normalization type
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%map_norm,       &
     &                     Import(Findex)%map_norm),                    &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      field reggriding method
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%map_type,       &
     &                     Import(Findex)%map_type),                    &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      destination field grid-cell type
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%dst_gtype,      &
     &                     Import(Findex)%destination_grid),            &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      destination field units
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%dst_units,      &
     &                     Import(Findex)%destination_units),           &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      source field grid-cell type
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%src_gtype,      &
     &                     Import(Findex)%source_grid),                 &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      source field units
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%src_units,      &
     &                     Import(Findex)%source_units),                &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      DATA NetCDF variable name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%nc_vname,       &
     &                     Import(Findex)%data_netcdf_vname),           &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!                                      DATA NetCDF time variable name
            IF (FoundError(assign_string(                               &
     &                     MODELS(Iroms)%ImportField(i)%nc_tname,       &
     &                     Import(Findex)%data_netcdf_tname),           &
     &                   NoError, __LINE__, MyFile)) THEN
              rc=ESMF_RC_COPY_FAIL
              RETURN
            END IF
!
!  Set grid type flag.
!
            SELECT CASE (lowercase(Import(Findex)%destination_grid))
              CASE ('center_cell', 'cell_center', 'center')
                MODELS(Iroms)%ImportField(i)%gtype=Icenter
              CASE ('corner_cell', 'cell_corner', 'corner')
                MODELS(Iroms)%ImportField(i)%gtype=Icorner
              CASE ('left_right_edge', 'right_left_edge', 'edge1')
                MODELS(Iroms)%ImportField(i)%gtype=Iupoint
              CASE ('lower_upper_edge', 'upper_lower_edge', 'edge2')
                MODELS(Iroms)%ImportField(i)%gtype=Ivpoint
              CASE DEFAULT
                MODELS(Iroms)%ImportField(i)%gtype=Icenter
            END SELECT
!
!  Set map type flag.
!
            SELECT CASE (lowercase(Import(Findex)%map_type))
              CASE ('mapbilnr')
                MODELS(Iroms)%ImportField(i)%itype=Ibilin
              CASE ('mappatch')
                 MODELS(Iroms)%ImportField(i)%itype=Ipatch
              CASE ('mapconsd')
                MODELS(Iroms)%ImportField(i)%itype=Iconsvd
              CASE ('mapconsf')
                MODELS(Iroms)%ImportField(i)%itype=Iconsvf
              CASE ('mapfcopy')
                MODELS(Iroms)%ImportField(i)%itype=Ifcopy
              CASE ('mapnstod')
                MODELS(Iroms)%ImportField(i)%itype=InStoD
              CASE ('mapnstod_consd')
                MODELS(Iroms)%ImportField(i)%itype=InStoDd
              CASE ('mapnstod_consf')
                MODELS(Iroms)%ImportField(i)%itype=InStoDf
              CASE DEFAULT
                MODELS(Iroms)%ImportField(i)%itype=Inone
            END SELECT
!
!  Check if field exits in NUOPC dictionary.
!
            Lexist=NUOPC_FieldDictionaryHasEntry(                       &
     &                     MODELS(Iroms)%ImportField(i)%standard_name,  &
     &                                        rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  If appropriate, add field to NUOPC dictionary.
!
            IF (.not.Lexist) THEN
              CALL NUOPC_FieldDictionaryAddEntry(                       &
     &                     MODELS(Iroms)%ImportField(i)%standard_name,  &
     &                                         canonicalUnits =         &
     &                     MODELS(Iroms)%ImportField(i)%src_units,      &
     &                                         rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
          ELSE
            IF (localPET.eq.0) THEN
              WRITE (cplout,30) 'import field short_name: ',            &
     &                          TRIM(ShortName), TRIM(CPLname)
            END IF
            rc=ESMF_RC_NOT_FOUND
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
!  Extract other configuration parameters from YAML object.
!-----------------------------------------------------------------------
!
!  Coupling type: explicit or semi-implicit.
!
      IF (YML%has('CouplingType')) THEN
        IF (FoundError(yaml_get(YML, 'CouplingType',                    &
     &                          CouplingType),                          &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
      ELSE
        rc=ESMF_RC_NOT_FOUND
        IF (localPET.eq.0) WRITE (cplout,20) 'CouplingType',            &
     &                                       TRIM(CPLname)
        RETURN
      END IF
!
!  PET layout: sequential or concurrent.
!
      IF (YML%has('PETlayoutOption')) THEN
        IF (FoundError(yaml_get(YML, 'PETlayoutOption',                 &
     &                          layout),                                &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
      END IF
!
      SELECT CASE (layout)
        CASE (0)
          PETlayoutOption(Iroms)='SEQUENTIAL'
        CASE (1)
          PETlayoutOption(Iroms)='CONCURRENT'
      END SELECT
!
!  Coupled state labels.
!
      IF (YML%has('CoupledSet')) THEN
        IF (FoundError(yaml_get(YML, 'CoupledSet',                      &
     &                          string),  &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
        COUPLED(Iroms)%SetLabel(linked_grid)=TRIM(string)
      END IF
!
      IF (YML%has('ExportState')) THEN
        IF (FoundError(yaml_get(YML, 'ExportState',                     &
     &                          string),                                &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
        COUPLED(Iroms)%ExpLabel(linked_grid)=TRIM(string)
      END IF
!
      IF (YML%has('ImportState')) THEN
        IF (FoundError(yaml_get(YML, 'ImportState',                     &
     &                          string),  &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
        COUPLED(Iroms)%ImpLabel(linked_grid)=TRIM(string)
      END IF
!
!  Coupling debugging flag.
!
      IF (YML%has('DebugLevel')) THEN
        IF (FoundError(yaml_get(YML, 'DebugLevel',                      &
     &                          DebugLevel),                            &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
      END IF
!
!  Execution tracing flag.
!
      IF (YML%has('TraceLevel')) THEN
        IF (FoundError(yaml_get(YML, 'TraceLevel',                      &
     &                          TraceLevel),                            &
     &                 NoError, __LINE__, MyFile)) THEN
          rc=ESMF_RC_COPY_FAIL
          RETURN
        END IF
      END IF
!
!  Destroy YAML obeject and deallocate local structures.
!
      CALL YML%destroy ()
      IF (allocated(Export))  deallocate (Export)
      IF (allocated(Import))  deallocate (Import)
      IF (allocated(Estring)) deallocate (Estring)
      IF (allocated(Istring)) deallocate (Istring)
!
!-----------------------------------------------------------------------
!  Report specified import and export states.
!-----------------------------------------------------------------------
!
      git_url=GIT_URL
      git_rev=GIT_REV
!
      IF (LocalPET.eq.0) THEN
        CALL get_date (TodayDateString)
        WRITE (cplout,40) TRIM(ESMF_VERSION_STRING),                    &
     &                    TRIM(TodayDateString),                        &
     &                    TRIM(git_url),                                &
     &                    TRIM(git_rev),                                &
     &                    TRIM(my_os),                                  &
     &                    TRIM(my_cpu),                                 &
     &                    TRIM(my_fort),                                &
     &                    TRIM(my_fc),                                  &
     &                    TRIM(my_fflags),                              &
     &                    MyComm, PETcount
        WRITE (cplout,50) 'Coupling Input Metadata Filename = ',        &
     &                    TRIM(CPLname)
        WRITE (cplout,50) '  ROMS Input Parameters Filename = ',        &
     &                    TRIM(Iname)
        WRITE (cplout,'(a)') CHAR(32)              ! blank space
        WRITE (cplout,60) '          Number of nested grids = ',        &
                          Ngrids
        WRITE (cplout,60) '             Coupled nested grid = ',        &
     &                    linked_grid
        IF (CouplingType.eq.1) THEN
          WRITE (cplout,70) '                   Coupling flag = ',      &
     &                      CouplingType,                               &
     &                      ',  Explicit coupling method'
        ELSE IF (CouplingType.eq.2) THEN
          WRITE (cplout,70) '                   Coupling flag = ',      &
     &                      CouplingType,                               &
     &                      ',  Semi-Implicit coupling method'
        END IF
        IF (layout.eq.0) THEN
          WRITE (cplout,70) '              PETs layout option = ',      &
     &                      layout,                                     &
     &                      ',  Sequential'
        ELSE IF (layout.eq.1) THEN
          WRITE (cplout,70) '              PETs layout option = ',      &
     &                      layout,                                     &
     &                      ',  Concurrent'
        END IF
        WRITE (cplout,'(a)') CHAR(32)
        WRITE (cplout,60) '   Coupling debugging level flag = ',        &
     &                    DebugLevel
        WRITE (cplout,60) '    Execution tracing level flag = ',        &
     &                    TraceLevel
!
        IF (Nimport(Iroms).gt.0) THEN
          WRITE (cplout,80) 'ROMS IMPORT Fields Metadata:'
          DO i=1,Nimport(Iroms)
            WRITE (cplout,90)                                           &
     &              TRIM(MODELS(Iroms)%ImportField(i)%short_name),      &
     &              TRIM(MODELS(Iroms)%ImportField(i)%standard_name),   &
     &              MODELS(Iroms)%ImportField(i)%gtype,                 &
     &              MODELS(Iroms)%ImportField(i)%itype,                 &
     &              MODELS(Iroms)%ImportField(i)%connected,             &
     &              MODELS(Iroms)%ImportField(i)%debug_write,           &
     &              MODELS(Iroms)%ImportField(i)%add_offset,            &
     &              MODELS(Iroms)%ImportField(i)%scale_factor
          END DO
        END IF
!
        IF (Nexport(Iroms).gt.0) THEN
          WRITE (cplout,80) 'ROMS EXPORT Fields Metadata:'
          DO i=1,Nexport(Iroms)
            WRITE (cplout,90)                                           &
     &              TRIM(MODELS(Iroms)%ExportField(i)%short_name),      &
     &              TRIM(MODELS(Iroms)%ExportField(i)%standard_name),   &
     &              MODELS(Iroms)%ExportField(i)%gtype,                 &
     &              MODELS(Iroms)%ExportField(i)%itype,                 &
     &              MODELS(Iroms)%ExportField(i)%connected,             &
     &              MODELS(Iroms)%ExportField(i)%debug_write,           &
     &              MODELS(Iroms)%ExportField(i)%add_offset,            &
     &              MODELS(Iroms)%ExportField(i)%scale_factor
          END DO
        END IF
        WRITE (cplout,100)
      END IF
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_Create',             &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
  10  FORMAT (/,' ROMS_CREATE - Unable to create YAML object for',      &
     &        ' ROMS/CMEPS configuration metadata file: ',/,15x,a,/,    &
     &        15x,'Default file is located in source directory.')
  20  FORMAT (/," ROMS_CREATE - Unable to find key: '",a,"'",           &
     &        ' ROMS/CMEPS configuration metadata file: ',/,15x,a,/,    &
     &        /,15x,'YAML file: ',a)
  30  FORMAT (/,' ROMS_CREATE - cannot find metadata for',              &
     &        1x,a,'''',a,'''.',/,15x,                                  &
     &        'Add entry to metadata file: ',a)
  40  FORMAT (80('-'),/,                                                &
     &        ' Earth System Models Coupling: ESMF/NUOPC Library,',     &
     &        ' Version ',a,/,31x,a,/,                                  &
     &        80('-'),                                                  &
     &        /,1x,'GIT Root URL     : ',a,                             &
     &        /,1x,'GIT Revision     : ',a,                             &
     &        /,1x,'Operating System : ',a,                             &
     &        /,1x,'CPU Hardware     : ',a,                             &
     &        /,1x,'Compiler System  : ',a,                             &
     &        /,1x,'Compiler Command : ',a,                             &
     &        /,1x,'Compiler Flags   : ',a,                             &
     &        /,1x,'MPI Communicator : ',i0,2x,'PET size = ',i0,        &
     &        /,80('-'),/)
  50  FORMAT (1x,a,a)
  60  FORMAT (1x,a,i0)
  70  FORMAT (1x,a,i0,a)
  80  FORMAT (/,a,/, 27('='),/,/, 'Short Name',                         &
     &        t15,'Standard Name', t74,'G', t77,'R', t80,'C', t83,'W',  &
     &        t87,'add_offset', t99,'scale_factor',/, 111('-'))
  90  FORMAT (a, t15,a, t74,i1, t77,i1, t80,l1, t83,l1,                 &
     &        t86,1p,e12.5, t100,1p,e12.5)
 100  FORMAT (/,' G: Grid cell location,     1=Center,',                &
     &                                     ' 2=Corner,',                &
     &                                     ' 3=U-point,',               &
     &                                     ' 4=V-point',                &
     &        /,' R: Regridding method,      1=bilin,',                 &
     &                                     ' 2=patch,',                 &
     &                                     ' 3=consD,',                 &
     &                                     ' 4=consF,',                 &
     &                                     ' 5=redist,',                &
     &                                     ' 6=nStoD,',                 &
     &                                     ' 7=nStoD-consD,',           &
     &                                     ' 8=nStoD-consF',            &
     &        /,' C: Connected to coupler,   F=derived from other,',    &
     &                                     ' T=exchanged/regridded',    &
     &        /,' W: Field write to NetCDF,  F=false, T=true',          &
     &                                     ' (used if DebugLevel > 2)'/)
!
      RETURN
      END SUBROUTINE ROMS_Create
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
      USE mod_parallel, ONLY : Master
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
      logical :: MasterPET, isPresent, isSet
!
      integer :: i
      integer :: ng = 1
      integer :: MyComm, PETcount, localPET
!
      TYPE (CouplingField), allocatable :: ROMSexport(:), ROMSimport(:)
      TYPE (ESMF_VM)                    :: vm
      TYPE (yaml_tree)                  :: YML
!
# ifdef ADD_NESTED_STATE
      character (len=100) :: CoupledSet, StateLabel
# endif
      character (len=240) :: StandardName, ShortName
      character (len=240) :: cfgValue
      character (len=160) :: message
!
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
      MasterPET=localPET.eq.0
      PETrank=localPET
!
!-----------------------------------------------------------------------
!  Set ROMS standard ouput unit and file
!-----------------------------------------------------------------------
!
!  Sets the ROMS standard output unit to write verbose execution info.
!  Notice that the default standard out unit in Fortran is 6.
!
!  In some applications like coupling or disjointed mpi-communications,
!  it is advantageous to write standard output to a specific filename
!  instead of the default Fortran standard output unit 6. If that is
!  the case, it opens such formatted file for writing.
!
      IF (Set_StdOutUnit) THEN
        Master=MasterPET
        stdout=stdout_unit(MasterPET)
        Set_StdOutUnit=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  Open standard output unit for ROMS cap information and messages.
!-----------------------------------------------------------------------
!
      OPEN (cplout, FILE=TRIM(CouplerLog), FORM='formatted',            &
     &      STATUS='replace')
!
!-----------------------------------------------------------------------
!  Query driver for attributes
!-----------------------------------------------------------------------
!
!  Get ROMS coupling configuration YAML filename.
!
      CALL NUOPC_CompAttributeGet (model,                               &
     &                             name='CouplingConfig',               &
     &                             value=cfgValue,                      &
     &                             isPresent=isPresent,                 &
     &                             isSet=isSet,                         &
     &                             rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (isPresent.and.isSet) THEN
        CPLname= TRIM(cfgValue)
        message='CouplingConfig = '//TRIM(CPLname)
        CALL ESMF_LogWrite (TRIM(message), ESMF_LOGMSG_INFO)
      ELSE
        message='CouplingConfig needs to be provided: '//               &
     &          'please check the top level ESMF config file.'
        CALL ESMF_LogSetError (ESMF_FAILURE, rcToReturn=rc,             &
     &                         msg=TRIM(message))
        RETURN
      END IF
!
!  Get 'scalar_field' parameters used by cmeps to check exchanged
!  fields.
!
      CALL NUOPC_CompAttributeGet (model,                               &
     &                             name='ScalarFieldName',              &
     &                             value=cfgValue,                      &
     &                             isPresent=isPresent,                 &
     &                             isSet=isSet,                         &
     &                             rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (isPresent.and.isSet) THEN
        READ (cfgValue,*) ScalarFieldName
      ELSE
        ScalarFieldName='cpl_scalars'
      END IF
      WRITE (message, '(a)') 'ScalarFieldName = '//                     &
     &                       TRIM(ScalarFieldName)
      CALL ESMF_LogWrite(TRIM(message), ESMF_LOGMSG_INFO)
!
      CALL NUOPC_CompAttributeGet(model,                                &
     &                            name='ScalarFieldCount',              &
     &                            value=cfgValue,                       &
     &                            isPresent=isPresent,                  &
     &                            isSet=isSet,                          &
     &                            rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (isPresent.and.isSet) THEN
        READ (cfgValue,*) ScalarFieldCount
      ELSE
        ScalarFieldCount=0
      END IF
      WRITE(message, '(a,i0)') 'ScalarFieldCount = ',                   &
     &                         ScalarFieldCount
      CALL ESMF_LogWrite(TRIM(message), ESMF_LOGMSG_INFO)
!
      CALL NUOPC_CompAttributeGet(model,                                &
     &                            name='ScalarFieldIdxGridNX',          &
     &                            value=cfgValue,                       &
     &                            isPresent=isPresent,                  &
     &                            isSet=isSet,                          &
     &                            rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (isPresent.and.isSet) THEN
        READ(cfgValue,*) ScalarFieldIdxGridNX
      ELSE
        ScalarFieldIdxGridNX = 0
      END IF
      WRITE(message, '(a,i0)') 'ScalarFieldIdxGridNX = ',               &
     &                         ScalarFieldIdxGridNX
      CALL ESMF_LogWrite (TRIM(message), ESMF_LOGMSG_INFO)
!
      CALL NUOPC_CompAttributeGet(model,                                &
     &                            name='ScalarFieldIdxGridNY',          &
     &                            value=cfgValue,                       &
     &                            isPresent=isPresent,                  &
     &                            isSet=isSet,                          &
     &                            rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (isPresent.and.isSet) THEN
        READ (cfgValue,*) ScalarFieldIdxGridNY
      ELSE
        ScalarFieldIdxGridNY=0
      END IF
      WRITE (message, '(a,i0)') 'ScalarFieldIdxGridNY = ',              &
     &                          ScalarFieldIdxGridNY
      CALL ESMF_LogWrite(TRIM(message), ESMF_LOGMSG_INFO)
!
!-----------------------------------------------------------------------
!  Create, allocate, and populate module structures. The export and
!  import fields metadata are read from input YAML configuration file.
!-----------------------------------------------------------------------
!
      CALL ROMS_Create (localPET, PETcount, MyComm, rc)
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

# ifdef ADD_NESTED_STATE
!
!  Add ROMS Import state. If nesting, each grid has its own import
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
            IF (LocalPET.eq.0) THEN
              WRITE (cplout,10) 'ROMS adding Import Nested State: ',    &
     &                          TRIM(StateLabel), ng
            END IF
!
!  Add fields import state.
!
            DO i=1,Nimport(Iroms)
              StandardName=MODELS(Iroms)%ImportField(i)%standard_name
              ShortName   =MODELS(Iroms)%ImportField(i)%short_name
              IF (LocalPET.eq.0) THEN
                WRITE (cplout,20) 'Advertising Import Field: ',         &
     &                            TRIM(ShortName), TRIM(StandardName)
              END IF
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

#  ifdef LONGWAVE_OUT
!
              IF (TRIM(ShortName).eq.'LWrad') THEN
                rc=ESMF_RC_NOT_VALID
                IF (localPET.eq.0) THEN
                  WRITE (cplout,30) TRIM(ShortName), 'LONGWAVE_OUT',    &
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
#  endif
            END DO
          END IF
        END DO
      END IF IMPORTING
# else
!
!  Add fields to ROMS Import state.  Coupled NestedStates are not
!  supported in cdeps/cmeps.
!
      IMPORTING : IF (Nimport(Iroms).gt.0) THEN
        ng=linked_grid
        MODELS(Iroms)%ImportState(ng)=ImportState
!
        IF (LocalPET.eq.0) THEN
          WRITE (cplout,10) 'ROMS Import STATE: ', ng
        END IF
!
        DO i=1,Nimport(Iroms)
          StandardName=MODELS(Iroms)%ImportField(i)%standard_name
          ShortName   =MODELS(Iroms)%ImportField(i)%short_name
          IF (LocalPET.eq.0) THEN
            WRITE (cplout,20) 'Advertising Import Field: ',             &
     &                        TRIM(ShortName), TRIM(StandardName)
          END IF
          CALL NUOPC_Advertise (MODELS(Iroms)%ImportState(ng),          &
     &                          StandardName=TRIM(StandardName),        &
     &                          name=TRIM(ShortName),                   &
     &                          rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF

#  ifdef LONGWAVE_OUT
!
          IF (TRIM(ShortName).eq.'LWrad') THEN
            rc=ESMF_RC_NOT_VALID
            IF (localPET.eq.0) THEN
              WRITE (cplout,30) TRIM(ShortName), 'LONGWAVE_OUT',        &
     &                          'downward longwave radiation: dLWrad',  &
     &                          'LONGWAVE_OUT'
            END IF
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
          END IF
#  endif
        END DO
      END IF IMPORTING
# endif
!
!-----------------------------------------------------------------------
!  Set ROMS Export State metadata.
!-----------------------------------------------------------------------

# ifdef ADD_NESTED_STATE
!
!  Add ROMS Export state. If nesting, each grid has its own export
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
            IF (LocalPET.eq.0) THEN
              WRITE (cplout,10) 'ROMS adding Export Nested State: ',    &
     &                          TRIM(StateLabel), ng
            END IF
!
!  Add fields to export state.
!
            DO i=1,Nexport(Iroms)
              StandardName=MODELS(Iroms)%ExportField(i)%standard_name
              ShortName   =MODELS(Iroms)%ExportField(i)%short_name
              IF (LocalPET.eq.0) THEN
                WRITE (cplout,20) 'Advertising Export Field: ',         &
     &                            TRIM(ShortName), TRIM(StandardName)
              END IF
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
# else
!
!  Add fields to ROMS Export state.  Coupled NestedStates are not
!  supported in cdeps/cmeps.
!
      EXPORTING : IF (Nexport(Iroms).gt.0) THEN
        ng=linked_grid
        MODELS(Iroms)%ExportState(ng)=ExportState
!
        IF (LocalPET.eq.0) THEN
          WRITE (cplout,10) 'ROMS Export STATE: ', ng
        END IF
!
        DO i=1,Nexport(Iroms)
          StandardName=MODELS(Iroms)%ExportField(i)%standard_name
          ShortName   =MODELS(Iroms)%ExportField(i)%short_name
          IF (LocalPET.eq.0) THEN
            WRITE (cplout,20) 'Advertising Export Field: ',             &
     &                        TRIM(ShortName), TRIM(StandardName)
          END IF
          CALL NUOPC_Advertise (MODELS(Iroms)%ExportState(ng),          &
     &                          StandardName=TRIM(StandardName),        &
     &                          name=TRIM(ShortName),                   &
     &                          rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
        END DO
      END IF EXPORTING
# endif
!
      IF (ESM_track) THEN
        WRITE (trac,'(a,a,i0)') '<== Exiting  ROMS_SetInitializeP1',    &
     &                          ', PET', PETrank
        FLUSH (trac)
      END IF
!
# ifdef ADD_NESTED_STATE
  10  FORMAT (/,a,a,', ng = ',i0,/,31('='),/)
# else
  10  FORMAT (/,a,'ng = ',i0,/,17('='),/)
# endif
  20  FORMAT (2x,a,"'",a,"'",t45,a)
# ifdef LONGWAVE_OUT
  30  FORMAT (/,' ROMS_SetInitializeP1 - incorrect field to process: ', &
     &        a,/,24x,'when activating option: ',a,/,24x,               &
     &        'use instead ',a,/,24x,'or deactivate option: ',a,/)
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
      integer :: ng, is, localPET, PETcount, tile
!
      real (dp) :: driverDuration, romsDuration
!
      character (len=20) :: TimeStartString
      character (len=20) :: TimeStopString

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_SetInitializeP2"
!
      TYPE (ESMF_TimeInterval) :: RunDuration, TimeStep
      TYPE (ESMF_Time)         :: CurrTime, startTime, stopTime
      TYPE (ESMF_CalKind_Flag) :: calkindflag
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
!  Query the Virtual Machine (VM) parallel environmemt for the MPI
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
!  Query driver clock.
!-----------------------------------------------------------------------
!
      IF (MODELS(Iroms)%IsActive) THEN
        CALL NUOPC_ModelGet (model,                                     &
     &                       driverClock=ClockInfo(Idriver)%Clock,      &
                             rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        CALL ESMF_ClockGet (ClockInfo(Idriver)%Clock,                   &
     &                      startTime=startTime,                        &
     &                      stopTime=stopTime,                          &
     &                      currTime=CurrTime,                          &
     &                      timeStep=TimeStep,                          &
     &                      calkindflag=calkindflag,                    &
     &                      runDuration=RunDuration,                    &
     &                      rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        ClockInfo(Idriver)%TimeStep = TimeStep
        ClockInfo(Idriver)%CurrentTime = CurrTime
        ClockInfo(Idriver)%StartTime = StartTime
        ClockInfo(Idriver)%StopTime = StopTime
!
        IF (calkindflag == ESMF_CALKIND_GREGORIAN) THEN
          ClockInfo(Idriver)%CalendarString = 'gregorian'
        ELSE IF (calkindflag == ESMF_CALKIND_360DAY) THEN
          ClockInfo(Idriver)%CalendarString = '360_day'
        END IF
!
        CALL ESMF_TimeGet (ClockInfo(Idriver)%StartTime,                &
     &                     s_r8=ClockInfo(Idriver)%Time_Start,          &
     &                     timeString=TimeStartString,                  &
     &                     rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        is=INDEX(TimeStartString, 'T')
        IF (is.gt.0) TimeStartString(is:is)=' '
        ClockInfo(Idriver)%Time_StartString=TRIM(TimeStartString)
!
        CALL ESMF_TimeGet (ClockInfo(Idriver)%StopTime,                 &
     &                     s_r8=ClockInfo(Idriver)%Time_Stop,           &
     &                     timeString=TimeStopString,                   &
     &                     rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        is=INDEX(TimeStopString, 'T')
        IF (is.gt.0) TimeStopString(is:is)=' '
        ClockInfo(Idriver)%Time_StopString=TRIM(TimeStopString)
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
     &                            TRIM(INPname(Iroms))
              END IF
              rc=ESMF_RC_NOT_VALID
              RETURN
            END IF
          END IF
        END DO
      END IF
!
!  Report Clock:
!
      IF (localPET.eq.0) THEN
        WRITE (cplout,20) TimeStartString, TimeStopString,              &
     &                    INT(driverDuration), INT(romsDuration)
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
  20  FORMAT (/,'Coupling Clock: ROMS_SetInitializeP2',/,15('='),/,     &
     &        /,2x,'DRIVER Starting Date = ',a,                         &
     &        /,2x,'DRIVER Ending   Date = ',a,                         &
     &        /,2x,'DRIVER Duration (s)  = ',i0,                        &
     &        /,2x,'ROMS   Duration (s)  = ',i0)
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
!  Using driver provided calendar since UFS Weather Model does not
!  support Proleptic Gregorian Calendar.
!
      ref_year  =Rclock%year
      ref_month =Rclock%month
      ref_day   =Rclock%day
      ref_hour  =Rclock%hour
      ref_minute=Rclock%minutes
      ref_second=Rclock%seconds
      Calendar  =TRIM(ClockInfo(Idriver)%CalendarString)
      ClockInfo(Iroms)%CalendarString = TRIM(Calendar)
!
      IF (ClockInfo(Iroms)%CalendarString == '360_day') THEN
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
!
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
        StartTimeString=TRIM(ClockInfo(Idriver)%Time_RestartString)
      ELSE
        StartTimeString=TRIM(ClockInfo(Idriver)%Time_StartString)
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
!  Sets ROMS component staggered, horizontal grids arrays and          !
!  land/sea mask, if any.                                              !
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
     &                                   gridEdgeLWidth=(/2,2/),        &
     &                                   gridEdgeUWidth=(/2,2/),        &
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
            staggerEdgeLWidth=(/1,1/)
            staggerEdgeUWidth=(/2,2/)
          CASE (Iupoint)
            staggerLoc=ESMF_STAGGERLOC_EDGE1
            staggerEdgeLWidth=(/1,1/)
            staggerEdgeUWidth=(/2,1/)
          CASE (Ivpoint)
            staggerLoc=ESMF_STAGGERLOC_EDGE2
            staggerEdgeLWidth=(/1,1/)
            staggerEdgeUWidth=(/1,2/)
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
!                                                        U-points
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
!                                                        PSI-points
            CASE (Icorner)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
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
!                              Extrapolate PSI-points at bottom edge
!
              IF (tile.lt.NtileI(ng)) THEN
                ptrX(:,Jstr-1)=2.0_dp*ptrX(:,Jstr)-ptrX(:,Jstr+1)
                ptrY(:,Jstr-1)=2.0_dp*ptrY(:,Jstr)-ptrY(:,Jstr+1)
                ptrM(:,Jstr-1)=ptrM(:,Jstr)
                ptrA(:,Jstr-1)=ptrA(:,Jstr)
              END IF
!                              Extrapolate PSI-points at left edge
!
              IF (MOD(tile,NtileI(ng)).eq.0) THEN
                ptrX(Istr-1,:)=2.0_dp*ptrX(Istr,:)-ptrX(Istr+1,:)
                ptrY(Istr-1,:)=2.0_dp*ptrY(Istr,:)-ptrY(Istr+1,:)
                ptrM(Istr-1,:)=ptrM(Istr,:)
                ptrA(Istr-1,:)=ptrA(Istr,:)
              END IF
!                              Extrapolate PSI-points at top edge
!
              IF (tile.ge.(NtileI(ng)*(NtileJ(ng)-1))) THEN
                ptrX(:,Jend+2)=2.0_dp*ptrX(:,Jend+1)-ptrX(:,Jend)
                ptrY(:,Jend+2)=2.0_dp*ptrY(:,Jend+1)-ptrY(:,Jend)
                ptrM(:,Jend+2)=ptrM(:,Jend+1)
                ptrA(:,Jend+2)=ptrA(:,Jend+1)
              END IF
!                             Extrapolate PSI-points at right edge
!
              IF (MOD(tile+1,NtileI(ng)).eq.0) THEN
                ptrX(Iend+2,:)=2.0_dp*ptrX(Iend+1,:)-ptrX(Iend,:)
                ptrY(Iend+2,:)=2.0_dp*ptrY(Iend+1,:)-ptrY(Iend,:)
                ptrM(Iend+2,:)=ptrM(Iend+1,:)
                ptrA(Iend+2,:)=ptrA(Iend+1,:)
              END IF
!                                                        U-points
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
!                             Extrapolate U-points at left edge
!
              IF (MOD(tile,NtileI(ng)).eq.0) THEN
                ptrX(Istr-1,:)=2.0_dp*ptrX(Istr,:)-ptrX(Istr+1,:)
                ptrY(Istr-1,:)=2.0_dp*ptrY(Istr,:)-ptrY(Istr+1,:)
                ptrM(Istr-1,:)=ptrM(Istr,:)
                ptrA(Istr-1,:)=ptrA(Istr,:)
              END IF
!                             Extrapolate U-points at right edge
!
              IF (MOD(tile+1,NtileI(ng)).eq.0) THEN
                ptrX(Iend+2,:)=2.0_dp*ptrX(Iend+1,:)-ptrX(Iend,:)
                ptrY(Iend+2,:)=2.0_dp*ptrY(Iend+1,:)-ptrY(Iend,:)
                ptrM(Iend+2,:)=ptrM(Iend+1,:)
                ptrA(Iend+2,:)=ptrA(Iend+1,:)
              END IF
!                                                        V-points
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
!                             Extrapolate V-points at bottom edge
!
              IF (tile.lt.NtileI(ng)) THEN
                ptrX(:,Jstr-1)=2.0_dp*ptrX(:,Jstr)-ptrX(:,Jstr+1)
                ptrY(:,Jstr-1)=2.0_dp*ptrY(:,Jstr)-ptrY(:,Jstr+1)
                ptrM(:,Jstr-1)=ptrM(:,Jstr)
                ptrA(:,Jstr-1)=ptrA(:,Jstr)
              END IF
!                             Extrapolate V-points at top edge
!
              IF (tile.ge.(NtileI(ng)*(NtileJ(ng)-1))) THEN
                ptrX(:,Jend+2)=2.0_dp*ptrX(:,Jend+1)-ptrX(:,Jend)
                ptrY(:,Jend+2)=2.0_dp*ptrY(:,Jend+1)-ptrY(:,Jend)
                ptrM(:,Jend+2)=ptrM(:,Jend+1)
                ptrA(:,Jend+2)=ptrA(:,Jend+1)
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
  10  FORMAT (/,'ROMS Domain Decomposition:',/,25('='),/,               &
              /,2x,'ROMS_DistGrid - Grid = ',i2.2,',',3x,'Mesh = ',a,   &
     &        ',',3x,'Partition = ',i0,' x ',i0)
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
      TYPE (ESMF_Field)      :: field, field_scalar
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
!  If cmeps scalar field, 'cpl_scalars', add it to Export State.
!
            IF (TRIM(ExportNameList(ifld)).eq.                          &
     &          TRIM(ScalarFieldName)) THEN
!
!  Create scalar field.
!
              CALL CreateScalarField(field,                             &
     &                               ScalarFieldName,                   &
     &                               ScalarFieldCount,                  &
     &                               rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  Set values and fill scalar field.
!
              CALL SetScalarFieldValues(field,                          &
     &                                  (/ Lm(ng)+2, Mm(ng)+2 /),       &
     &                                  (/ ScalarFieldIdxGridNX,        &
     &                                     ScalarFieldIdxGridNY /),     &
     &                                  rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
!
!  Otherwise, proccess regular ROMS export field(s).
!
            ELSE
!
!  Set staggering type.
!
              SELECT CASE (MODELS(Iroms)%ExportField(id)%gtype)
                CASE (Icenter)                                ! RHO-points
                  staggerLoc=ESMF_STAGGERLOC_CENTER
                CASE (Icorner)                                ! PSI-points
                  staggerLoc=ESMF_STAGGERLOC_CORNER
                CASE (Iupoint)                                ! U-points
                  staggerLoc=ESMF_STAGGERLOC_EDGE1
                CASE (Ivpoint)                                ! V-points
                  staggerLoc=ESMF_STAGGERLOC_EDGE2
              END SELECT
!
!  Create 2D field from the Grid and arraySpec.
!
              field=ESMF_FieldCreate(MODELS(Iroms)%grid(ng),            &
     &                               arraySpec2d,                       &
     &                               indexflag=ESMF_INDEX_GLOBAL,       &
     &                               staggerloc=staggerLoc,             &
     &                               name=TRIM(ExportNameList(ifld)),   &
     &                               rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
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
                CALL ESMF_FieldGet (field,                              &
     &                              localDe=localDE,                    &
     &                              farrayPtr=ptr2d,                    &
     &                              rc=rc)
                IF (ESMF_LogFoundError(rcToCheck=rc,                    &
     &                                 msg=ESMF_LOGERR_PASSTHRU,        &
     &                                 line=__LINE__,                   &
     &                                 file=MyFile)) THEN
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
            END IF
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
            IF (TRIM(ExportNameList(id)).ne.TRIM(ScalarFieldName)) THEN
              IF (localPET.eq.0) THEN
                WRITE (cplout,10) TRIM(ExportNameList(ifld)),           &
     &                            'Export State: ',                     &
     &                            TRIM(COUPLED(Iroms)%ExpLabel(ng))
              END IF
              CALL ESMF_StateRemove (MODELS(Iroms)%ExportState(ng),     &
     &                               (/ TRIM(ExportNameList(ifld)) /),  &
     &                               rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
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
              CASE (Icorner)                                ! PSI-points
                staggerLoc=ESMF_STAGGERLOC_CORNER
              CASE (Iupoint)                                ! U-points
                staggerLoc=ESMF_STAGGERLOC_EDGE1
              CASE (Ivpoint)                                ! V-points
                staggerLoc=ESMF_STAGGERLOC_EDGE2
            END SELECT
!
!  Create 2D field from the Grid, arraySpec, total tile size.
!  The array indices are global following ROMS design.
!
            field=ESMF_FieldCreate(MODELS(Iroms)%grid(ng),              &
     &                             arraySpec2d,                         &
     &                             indexflag=ESMF_INDEX_GLOBAL,         &
     &                             staggerloc=staggerLoc,               &
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
      FLUSH (cplout)                      ! flush coupling output buffer
      CLOSE (cplout)                      ! close coupling log file
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
!  Surface air relative humidity (percentage). Notice that as the
!  specific humidity, it is loaded to FORCES(ng)%Hair and "bulk_flux.F"
!  will compute the specific humidity (kg/kg).
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
# if defined BULK_FLUXES & defined PBLH
!
!  Atmosphere boundary layer thickness (m).
!
            CASE ('pblh', 'PBLh')
              romsScale=scale
              ifield=idPBLh
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
                  FORCES(ng)%PBLh(i,j)=Fval
                END DO
              END DO
              IF (localDE.eq.localDEcount-1) THEN
                IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                  CALL exchange_r2d_tile (ng, tile,                     &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    FORCES(ng)%PBLh)
                END IF
                CALL mp_exchange2d (ng, tile, iNLM, 1,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              FORCES(ng)%PBLh)
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
     &                            TRIM(CPLname)
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
     &      MODELS(Iroms)%ImportField(id)%debug_write) THEN
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
!  Skip if it is scalar field.
!
        IF (TRIM(ExportNameList(ifld)) == TRIM(ScalarFieldName)) THEN
          CYCLE
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

# if defined MASKING
!
!  Update wet point land/sea mask, if differs from static mask.
!
            CASE ('mask_rho', 'rmask', 'msk')
              MyFmin(1)=1.0_dp
              MyFmax(1)=0.0_dp
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
# ifdef WET_DRY
                    IF (GRID(ng)%rmask(i,j).ne.                         &
     &                  GRID(ng)%rmask_wet(i,j)) THEN
                      ptr2d(i,j)=GRID(ng)%rmask_wet(i,j)
                    ELSE
                      ptr2d(i,j)=GRID(ng)%rmask(i,j)
                    END IF
# else
                    ptr2d(i,j)=GRID(ng)%rmask(i,j)
# endif
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
     &                            TRIM(CPLname)
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
     &        a,/,18x,'check ''Export(roms)'' in input YAML: ',a)
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
#  ifdef MASKING
              Uout(i,j)=Uout(i,j)*GRID(ng)%umask(i,j)
#  endif
#  ifdef WET_DRY
              Uout(i,j)=Uout(i,j)*GRID(ng)%umask_wet(i,j)
#  endif
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              Vout(i,j)=0.5_r8*(Vrot(i,j-1)+Vrot(i,j))
#  ifdef MASKING
              Vout(i,j)=Vout(i,j)*GRID(ng)%vmask(i,j)
#  endif
#  ifdef WET_DRY
              Vout(i,j)=Vout(i,j)*GRID(ng)%vmask_wet(i,j)
#  endif
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
#  ifdef WET_DRY
            Uout(i,j)=Uout(i,j)*GRID(ng)%rmask_wet(i,j)
            Vout(i,j)=Vout(i,j)*GRID(ng)%rmask_wet(i,j)
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
#  ifdef MASKING
            Uout(i,j)=Uout(i,j)*GRID(ng)%umask(i,j)
#  endif
#  ifdef WET_DRY
            Uout(i,j)=Uout(i,j)*GRID(ng)%umask_wet(i,j)
#  endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            Vout(i,j)=0.5_r8*(Vinp(i,j-1)+Vinp(i,j))
#  ifdef MASKING
              Vout(i,j)=Vout(i,j)*GRID(ng)%vmask(i,j)
#  endif
#  ifdef WET_DRY
              Vout(i,j)=Vout(i,j)*GRID(ng)%vmask_wet(i,j)
#  endif
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
      INTEGER FUNCTION field_index (Fnames, Fvalue) RESULT (Findex)
!
!=======================================================================
!                                                                      !
!  This integer function scans an array structure of type ESM_Field    !
!  containing fields short_name list for specific field value and      !
!  returns its location index in the list.                             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: Fvalue

      TYPE (ESM_Field), intent(in)  :: Fnames(:)
!
!  Local variable declarations.
!
      integer :: Mfields
      integer :: i
!
!-----------------------------------------------------------------------
!  Find index of specified field from names list.
!-----------------------------------------------------------------------
!
      Mfields=SIZE(Fnames, DIM=1)
      Findex=-1
!
      DO i=1,Mfields
        IF (TRIM(Fnames(i)%short_name).eq.TRIM(Fvalue)) THEN
          Findex=i
          EXIT
        END IF
      END DO
!
      RETURN
      END FUNCTION field_index
!
      SUBROUTINE report_timestamp (field, CurrTime, localPET, string,   &
     &                             rc)
!
!=======================================================================
!                                                                      !
!  Reports coupling time-stamp.                                        !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)  :: localPET
      integer, intent(out) :: rc
!
      character (len=*), intent(in) :: string
!
      TYPE (ESMF_Field), intent(in) :: field
      TYPE (ESMF_Time),  intent(in) :: CurrTime
!
!  Local variable declarations.
!
      logical :: IsValid
      integer :: vtime1(10), vtime2(10)
!
      TYPE (ESMF_Time) :: FieldTime
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", report_timestamp"

      character (len=22) :: str1, str2
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get driver current time.
!-----------------------------------------------------------------------
!
      CALL ESMF_TimeGet (CurrTime,                                      &
     &                   yy=vtime1(1),                                  &
     &                   mm=vtime1(2),                                  &
     &                   dd=vtime1(3),                                  &
     &                   h =vtime1(4),                                  &
     &                   m =vtime1(5),                                  &
     &                   s =vtime1(6),                                  &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      WRITE (str1,10) vtime1(1), vtime1(2), vtime1(3),                  &
     &                vtime1(4), vtime1(5), vtime1(6)
!
!-----------------------------------------------------------------------
!  Get field TimeStamp.
!-----------------------------------------------------------------------
!
      CALL NUOPC_GetTimeStamp (field,                                   &
     &                         isValid = IsValid,                       &
     &                         time = FieldTime,                        &
     &                         rc = rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (IsValid) THEN
        CALL ESMF_TimeGet (FieldTime,                                   &
     &                     yy=vtime2(1),                                &
     &                     mm=vtime2(2),                                &
     &                     dd=vtime2(3),                                &
     &                     h =vtime2(4),                                &
     &                     m =vtime2(5),                                &
     &                     s =vtime2(6),                                &
     &                     rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc,                            &
     &                         msg=ESMF_LOGERR_PASSTHRU,                &
     &                         line=__LINE__,                           &
     &                         file=MyFile)) THEN
          RETURN
        END IF
!
        WRITE (str2,10) vtime2(1), vtime2(2), vtime2(3),                &
     &                  vtime2(4), vtime2(5), vtime2(6)
      END IF
!
!-----------------------------------------------------------------------
!  Report TimeStamp.
!-----------------------------------------------------------------------
!
      IF (IsValid) THEN
        IF (TRIM(str1).ne.TRIM(str2)) THEN
          IF (localPET.eq.0) THEN
            WRITE (cplout,20) TRIM(string), ': TimeStamp = ',           &
     &                        TRIM(str2), ' not equal ' ,               &
     &                        TRIM(str1)
          END IF
          rc=ESMF_RC_VAL_WRONG
          RETURN
        ELSE
          IF (localPET.eq.0) THEN
!!          WRITE (cplout,30) TRIM(string), ': TimeStamp = ', TRIM(str2)
          END IF
        END IF
      ELSE
        IF (localPET.eq.0) THEN
          WRITE (cplout,30) TRIM(string), ': TimeStamp is not valid',   &
     &                      ', DriverTime = '//TRIM(str1)
        END IF
      END IF
!
  10  FORMAT (i4.4,2('-',i2.2),1x,i2.2,':',i2.2,':',i2.2)
  20  FORMAT (/,1x,a,a,a,a,a)
  30  FORMAT (1x,a,a,a)
!
      RETURN
      END SUBROUTINE report_timestamp
!
      SUBROUTINE CreateScalarField (Field, Fname, Fcount, rc)
!
!=======================================================================
!                                                                      !
!  Creeates ESMF field to store set scalars                            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)  :: Fcount
      integer, intent(out) :: rc
!
      character (len=*), intent(in) :: Fname
!
      TYPE (ESMF_Field) , intent(inout) :: Field
!
!  Local variable declarations.
!
      TYPE (ESMF_Grid)     :: grid
      TYPE (ESMF_DistGrid) :: distgrid
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", CreateScalarField"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Create ESMF field
!-----------------------------------------------------------------------
!
      distgrid = ESMF_DistGridCreate(minIndex=(/1/),                    &
     &                               maxIndex=(/1/),                    &
     &                               rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      grid = ESMF_GridCreate(distgrid, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      Field = ESMF_FieldCreate(name=trim(Fname),                        &
     &                         grid=grid,                               &
     &                         typekind=ESMF_TYPEKIND_R8,               &
     &                         gridToFieldMap=(/2/),                    &
     &                         ungriddedLBound=(/1/),                   &
     &                         ungriddedUBound=(/Fcount/),              &
     &                         rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      RETURN
      END SUBROUTINE CreateScalarField
!
      SUBROUTINE SetScalarFieldValues (field, vals, idxs, rc)
!
!=======================================================================
!                                                                      !
!  Sets scalar field.                                                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)    :: vals(:)
      integer, intent(in)    :: idxs(:)
      integer, intent(inout) :: rc
!
      TYPE (ESMF_Field), intent(inout) :: field
!
!  Local variable declarations.
!
      integer                     :: I
      integer                     :: ungriddedLBound(2)
      integer                     :: ungriddedUBound(2)
!
      real (r8), pointer          :: fPtr(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", SetScalarFieldValues"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Get pointer and fill the field.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldGet(field,                                         &
     &                   ungriddedLBound=ungriddedLBound,               &
     &                   ungriddedUBound=ungriddedUBound,               &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      IF (ANY(idxs(:).lt.ungriddedLBound(1)) .or.                       &
          ANY(idxs(:).gt.ungriddedUBound(1))) THEN
        CALL ESMF_LogSetError(ESMF_RC_ARG_BAD,                          &
     &    msg="ScalarFieldIdxGrid[NX|NY] outside ScalarFieldCount",     &
     &    line=__LINE__, file=__FILE__, rcToReturn=rc)
      END IF
!
      CALL ESMF_FieldGet(field, farrayPtr=fPtr, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      DO I=LBOUND(idxs, dim=1), UBOUND(idxs, dim=1)
        fPtr(idxs(I),1) = DBLE(vals(I))
      END DO
!
      RETURN
      END SUBROUTINE SetScalarFieldValues
!
#endif
      END MODULE cmeps_roms_mod
