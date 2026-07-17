#include "cppdefs.h"
      MODULE mod_esmf_esm
#if defined MODEL_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: mod_esmf_esm.F 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module defines the structures and variables for Earth System   !
!  Model (ESM) coupling using the ESMF/NUOPC library. It also includes !
!  several support routines:                                           !
!                                                                      !
!  allocate_esmf_esm     Allocatess and initializes module structures  !
!                        variables.                                    !
!                                                                      !
!  clock_report          Inquire ESM driver and component clocks and   !
!                        reports current information.                  !
!                                                                      !
!  def_FieldAtt          Creates output NetCDF file containing the     !
!                        field attribute varaibles needed to perform   !
!                        the time interpolation in concurrent coupling !
!                        and when ESM components import time snapshots.!
!                                                                      !
!  field_index           Scans an array structure of type ESM_Field    !
!                        containing a list of export or import fields  !
!                        short  names (short_name) for requested field !
!                        and returns its location index in the list.   !
!                                                                      !
!  get_atm_Ngrids        Reads atmosphere model number of nested grids !
!                        (NgridsA) from input script or namelist.      !
!                                                                      !
!  get_roms_Ngrids       Reads ROMS number of nested grids (NgridsR)   !
!                        from input script (like ocean.in).            !
!                                                                      !
!  get_weights           Sets or reads in melding weights coefficients !
!                        needed by the atmosphere model to merge fields!
!                        from DATA and other ESM components because of !
!                        incongruent grids.                            !
!                                                                      !
!  load_IFS              Loads DATA Model source files information     !
!                        into the DataSet(:)%IFS structure.            !
!                                                                      !
!  read_ESMconfig        Reads in and reports ESM coupling coupling    !
!                        configuration parameter from standard input   !
!                        script (like coupling_esmf.in).               !
!                                                                      !
!  report_timestamp      Reports coupling time-stamp for requested     !
!                        field.                                        !
!                                                                      !
!  set_metadata          Process ESM coupling import and export fields !
!                        metadata and adds fields to NUOPC dictionary. !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!=======================================================================
!
      USE ESMF
      USE NUOPC
!
      USE mod_kinds
      USE mod_iounits, ONLY : T_IO
!
      implicit none
!
!-----------------------------------------------------------------------
!  ESM support routines.
!-----------------------------------------------------------------------
!
      PRIVATE :: load_IFS
!
      PUBLIC  :: allocate_esmf_esm
      PUBLIC  :: clock_report
      PUBLIC  :: def_FieldAtt
      PUBLIC  :: field_index
      PUBLIC  :: get_atm_Ngrids
      PUBLIC  :: get_roms_Ngrids
      PUBLIC  :: get_weights
      PUBLIC  :: read_ESMconfig
      PUBLIC  :: report_timestamp
      PUBLIC  :: set_metadata
!
!-----------------------------------------------------------------------
!  ESM generic data types structures.
!-----------------------------------------------------------------------
!
!  ESM coupling time managing variables and ESMF objects.

      TYPE ESM_Clock
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
! ESM coupled state sets. If appropriate, it includes the logic for
! connecting nested grids. The export state COUPLED(:)%ExportState(:,:)
! will be only allocated and used here by the DATA component, for
! example exporting to a specified atmosphere nested grid 'ng':
!
!    COUPLED(Idata)%ExportState(ng,Iatmos)
!
! Notice that it is matrix because it is possible that the same dataset
! will be used multiple times for different coupled sets (during
! nesting), and likely, also needed by other ESM component(s).
!
! In other ESM components, the import and export states are defined and
! allocated in the MODELS structure as vectors (see below).
!
      TYPE ESM_CplSet
        logical, allocatable :: LinkedGrid(:,:)       ! connected grid
!
        integer, allocatable :: DataCoupledSets(:,:)  ! DATA linked sets
!
        character (len=100), allocatable :: SetLabel(:) ! set label
        character (len=100), allocatable :: ExpLabel(:) ! export label
        character (len=100), allocatable :: ImpLabel(:) ! import label
!
        TYPE(ESMF_State), allocatable :: ExportState(:,:) ! export set
      END TYPE ESM_CplSet
!
      TYPE (ESM_CplSet), allocatable, target :: COUPLED(:)
!
!  DATA model field processing information.
!
      TYPE ESM_Data
        logical :: Lcycle                     ! cycling time coordinate
        logical :: Lcoord                     ! coordinates attribute
        logical :: Lmask                      ! land/sea mask
        logical :: Lmulti                     ! field across multi-files
        logical :: LastRec                    ! processed last record
!
        integer :: Icomp                      ! target component index
        integer :: ncid                       ! NetCDF file ID
        integer :: Tid                        ! time variable ID
        integer :: Vid                        ! field variable ID
        integer :: Nvdim                      ! number spatial dimensions
        integer :: Zlevel                     ! level index to process
        integer :: Nrec                       ! number of time records
        integer :: Trec                       ! latest read time record
        integer :: Tindex                     ! rolling two-time indices
        integer :: LandValue                  ! Masking land value
        integer :: SeaValue                   ! Masking sea value
!
        integer, allocatable :: Vsize(:)      ! dimensions size
!
        real(dp) :: add_offset                ! add_offset attribute
        real(dp) :: FillValue                 ! _FillValue attribute
        real(dp) :: scale_factor              ! scale_factor attribute
        real(dp) :: Clength                   ! time cycling length
        real(r8) :: LonMin                    ! grid minimum longitude
        real(r8) :: LonMax                    ! grid maximum longitude
        real(r8) :: LatMin                    ! grid minimum latitude
        real(r8) :: LatMax                    ! grid maximum latitude
        real(dp) :: Tscale                    ! time scale to day
        real(dp) :: Tmono                     ! monotonic time (days)
        real(dp) :: Tmin                      ! time minimum value
        real(dp) :: Tmax                      ! time maximum value
        real(dp) :: Tstr                      ! lower time-snapshot
        real(dp) :: Tend                      ! upper time-snapshot
        real(r8) :: Vmin                      ! variable minimum value
        real(r8) :: Vmax                      ! variable maximum value
        real(dp) :: Tintrp(2)                 ! interpolation time (days)
        real(dp) :: Vtime(2)                  ! latest two-time values
        real(dp) :: Date(6,2)                 ! time-snapshots dates
!                                               YYYY,MM,DD hh:mm:ss.ss
        character(len=30), allocatable :: Dname(:)  ! variable
!                                                     dimensions names
        character(len=20), allocatable :: Vcoord(:) ! variable
!                                                     coordinates names
        character(len=22 ) :: DateString(2)   ! date-snapshots string
        character(len=30 ) :: SpecialAction   ! special processing
        character(len=100) :: Tname           ! time variable name
        character(len=100) :: Tunits          ! time variable units
        character(len=100) :: Vname           ! variable name
        character(len=100) :: Vunits          ! variable units
        character(len=256) :: Vdescriptor     ! Variable descriptive name
        character(len=256) :: Vlongname       ! long_name attribute
        character(len=256) :: ncfile          ! NetCDF filename
!
        real(r8), allocatable :: lon(:,:)     ! field longitude
        real(r8), allocatable :: lat(:,:)     ! field latitude
        real(r8), allocatable :: mask(:,:)    ! field land/sea mask
!
        real(r8), allocatable :: A2d(:,:)     ! time interpolated 2D data
        real(r8), allocatable :: A3d(:,:,:)   ! time interpolated 3d data
        real(r8), allocatable :: A2dG(:,:,:)  ! latest 2D data snapsnots
        real(r8), allocatable :: A3dG(:,:,:,:)! latest 3D data snapsnots
!
        TYPE (ESMF_Field)       :: field      ! field object
        TYPE (ESMF_Grid)        :: grid       ! field grid object
        TYPE (ESMF_RouteHandle) :: rhandle    ! field RouteHandle
      END TYPE ESM_Data
!
!  DATA model high-level structure.
!
      TYPE ESM_DataSet
        integer :: Nfields                    ! number of fields
        integer :: Nfiles                     ! number of input files
!
        character(len=20), allocatable :: Ctarget(:)  ! component target
        character(len=20), allocatable :: Field(:)    ! field short-name
!
        TYPE (ESM_Data),   allocatable :: Export(:)   ! Export field
!                                                       structure
        TYPE (T_IO),       allocatable :: IFS(:)      ! Input Files
!                                                       Structure
      END TYPE ESM_DataSet
!
!  Import and export fields metadata information.
!
      TYPE ESM_Field
        logical :: connected                  ! connected to coupler
        logical :: debug_write                ! write exchanged field
        logical :: enable_integral_adj        ! area integral adjusted

        integer :: fid                        ! internal field ID
        integer :: gtype                      ! field grid mesh type
        integer :: etype                      ! field extrapolation flag
        integer :: itype                      ! field interpolation flag
        integer :: Tindex                     ! rolling two-time indices

        character (len=20)  :: Ctarget        ! component destination
        character (len=22)  :: DateString(2)  ! date-snapshots string

        character (len=:), allocatable :: short_name    ! short name
        character (len=:), allocatable :: standard_name ! standard name
        character (len=:), allocatable :: long_name     ! long name
        character (len=:), allocatable :: dst_gtype     ! DST grid type
        character (len=:), allocatable :: dst_units     ! DST units
        character (len=:), allocatable :: src_gtype     ! SRC grid type
        character (len=:), allocatable :: src_units     ! SRC units
        character (len=:), allocatable :: nc_vname      ! DATA Vname
        character (len=:), allocatable :: nc_tname      ! DATA Tname
        character (len=:), allocatable :: RegridMethod  ! regrid method
        character (len=:), allocatable :: ExtrapMethod  ! extrapolate

        real (r8) :: scale_factor             ! field scale factor
        real (r8) :: add_offset               ! field add offset value
        real (r8) :: Tmin                     ! DATA time minimum value
        real (r8) :: Tmax                     ! DATA time maximum value
        real (r8) :: Tstr                     ! DATA lower time-snapshot
        real (r8) :: Tend                     ! DATA upper time-snapshot
        real (r8) :: Tintrp(2)                ! interpolation time (days)
        real (r8) :: Vtime(2)                 ! latest two-time values

        TYPE (ESMF_RouteHandle) :: rhandle    ! field RouteHandle
      END TYPE ESM_Field
!
!  Import and export fields mesh data.
!
      TYPE ESM_Mesh
        integer :: gid                           ! grid ID
        integer :: gtype                         ! grid mesh type

        integer (i4b), allocatable :: mask(:,:)  ! grid land/sea mask

        real (r8), allocatable :: lon(:,:)       ! grid longitude
        real (r8), allocatable :: lat(:,:)       ! grid latitude
        real (r8), allocatable :: area(:,:)      ! grid area
      END TYPE ESM_Mesh
!
!  Melding coefficients used to combine fields from DATA and ESM
!  components. The weight factors are read from the input NetCDF
!  specified in the "WeightsFile(atmos)"  keyword.  The user has
!  full control of how the merging is done. It is recommended to
!  provide a gradual transition between the two components.
!
!  Recall that the DATA component supplies needed data to a particular
!  ESM component.  For example, it may export data to the atmosphere
!  model at locations not covered by the other ESM components because
!  of smaller grid coverage.  If the atmosphere and ocean model grids
!  are incongruent, the atmosphere component needs to import sea surface
!  temperature (SST) on those grid points not covered by the ocean
!  component. Thus, the weighting coefficients are used to merge the
!  SST data:
!
!  SST_atm(:,:) = Cesm(:,:) * SST_esm(;,;) + Cdat(:,:) * SST_dat(:,:)
!
!  where Cesm(:,:) + Cdat(:,:) = 1.
!
      TYPE ESM_Meld
        integer :: NestedGrid               ! grid needing merged field

        character(len=100) :: VnameDATA     ! DATA weights variable name
        character(len=100) :: VnameESM      ! ESM  weights variable name
        character(len=256) :: ncfile        ! Weights NetCDF filename

        real (r8), allocatable :: Cdat(:,:) ! coefficients for DATA
        real (r8), allocatable :: Cesm(:,:) ! coefficients for ESM
      END TYPE ESM_Meld
!
      TYPE (ESM_Meld), allocatable, target :: WEIGHTS(:)
!
      real(dp) :: WeightDAT = 0.0_dp      ! DATA component weight
      real(dp) :: WeightESM = 1.0_dp      ! ESM  component weight
!
!  Coupled models high-level data structure, [Nmodels].
!
!  The coupling field (Import/Export) exchange is determined by the
!  TimeStep/TimeFrac. Where TimeStep is the coupling driver interval
!  and TimeFrac is the fraction from Time/Step between the connector
!  between two ESM components.
!
      TYPE ESM_Model
        logical :: IsActive                       ! active for coupling

        integer (i4b) :: LandValue                ! land mask value
        integer (i4b) :: SeaValue                 ! sea  mask value

        integer :: Ngrids                         ! number nested grids

        integer :: ExportCalls                    ! export CALL counter
        integer :: ImportCalls                    ! import CALL counter

        integer :: nPETs                          ! number model PETs
        integer, allocatable :: PETlist(:)        ! model PETs list

        integer, allocatable :: TimeFrac(:,:)     ! driver time fraction

        character (len=100) :: name               ! model name

        TYPE (ESMF_Grid), allocatable :: grid(:)         ! grid object
        TYPE (ESM_Mesh),  allocatable :: mesh(:)         ! mesh
        TYPE (ESM_Field), allocatable :: ImportField(:)  ! import fields
        TYPE (ESM_Field), allocatable :: ExportField(:)  ! export fields
        TYPE (ESMF_State), allocatable :: ImportState(:) ! import state
        TYPE (ESMF_State), allocatable :: ExportState(:) ! export state
      END TYPE ESM_Model
!
      TYPE (ESM_Model), allocatable, target :: MODELS(:)
!
!  Coupling models connector used for the interpolation/extrapolaton
!  between source and destination fields, [Nmodels, Nmodels].
!
      TYPE ESM_Conn
        logical :: IsActive                 ! active connector

        integer :: divDT
        integer :: MaskInteraction          ! connector mask interaction

        integer :: nPETs                    ! number of connector PETs
        integer, allocatable :: PETlist(:)  ! connector PETs list

        character (len=100) :: name         ! connector name
      END TYPE ESM_Conn
!
      TYPE (ESM_Conn), allocatable, target :: CONNECTORS(:,:)
!
!  ESM import and export fields dictionary.  The fields are read from
!  metadata file (CPLname).
!
      integer, parameter :: MaxNumberFields = 200
      integer :: Nfields                          ! processed and loaded
!
      TYPE (ESM_Field) :: Fields(MaxNumberFields)
!
!  Define DATA Model field processing information from input data files,
!  [1:Nmodels]. Currently, the DATA Model only export fields and it
!  only supports input NetCDF files.
!
      TYPE (ESM_DataSet), allocatable, target :: DataSet(:)
!
!  ESM clock for driver (zeroth element) and coupled components,
!  [0:Nmodels].
!
      TYPE (ESM_Clock), allocatable, target :: ClockInfo(:)
!
!-----------------------------------------------------------------------
!  ESM coupling parameters.
!-----------------------------------------------------------------------
!
!  Number of coupled ESM gridded components. Currently, five
!  ESM components are supported (ROMS, DATA, Atmosphere, Sea-ice, and
!  wave model).
!
!  All supported components are accounted here even if we are running
!  an application with less number of models. The IsActive switches
!  are use to operate only on the desired coupled components.  This
!  is done to have complete infornation in the above structures.  The
!  gridded arrays are never allocated if a particular component is not
!  active.
!
      integer :: Nmodels = 5
!
!  Number of nested grids.  Initialize to just one grid; its values
!  are overwritten during processing.
!
!  (An additional variable NgridsR is created to avoid using ROMS
!   module "mod_param" in the generic interface.  Both Ngrids and
!   NgridsR have the same value)
!
      integer :: NgridsA = 1                  ! Atmosphere Model
      integer :: NgridsD = 1                  ! DATA Model
      integer :: NgridsI = 1                  ! Sea-ice Model
      integer :: NgridsR = 1                  ! ROMS
      integer :: NgridsW = 1                  ! Wave Model
!
!  Coupled models identification indices.  ROMS needs to be the
!  first index since we are using several of its modules to
!  initialize the coupled system.
!
      integer :: Idriver = 0
      integer :: Iroms   = 1
      integer :: Iatmos  = 2
      integer :: Idata   = 3
      integer :: Iseaice = 4
      integer :: Iwave   = 5
!
!  Generic ESM component labels used in the CASE constructs.  We cannot
!  use the identification indices because the vector Iroms(:) cannot
!  be defined as a parameter and a non constant expression is illegal:
!
!       CASE ( Iroms(1) : Iroms(NgridsR) )
!
      character (len= 3), allocatable :: Clabel(:)
      character (len=10), allocatable :: Cmodel(:)
!
!  Number of ESM import and export fields, [Nmodels].
!
      integer, allocatable :: Nimport(:)
      integer, allocatable :: Nexport(:)
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
!  Driver virtual Machine (VM) parallel enviroment object.
!
      TYPE (ESMF_VM) :: VMdriver
!
!  PET layout: sequential or concurrent.
!
      character (len=10) :: PETlayoutOption
!
!  Total number of PETs needed in concurrent PET layout and rank for
!  each PET.
!
      integer :: sumPETs
      integer :: PETrank
!
!  MPI Communicator handle for each ESM component.
!
      integer, allocatable :: ESMcomm(:)
!
!  Driver clock parameters specified in configuration script. A integer
!  vector with six elements:
!
!     (1) year including century, like 2017
!     (2) month of the year, 1 to 12
!     (3) day of the month
!     (4) hour of the day, 0 to 23
!     (5) minutes of the hour, 0 to 59
!     (6) seconds of the minute, 0 to 59
!
      integer :: ReferenceDate(6)          ! reference date
      integer :: RestartDate(6)            ! restarting date
      integer :: StartDate(6)              ! starting date
      integer :: StopDate(6)               ! stopping date
      integer :: TimeStep(6)               ! coupling interval
!
!  Today date string.
!
      character (len=44) :: TodayDateString
!
!  ESM coupling simulation reference date number:
!
!      (1) seconds
!      (2) fractional days
!
       real(dp) :: ReferenceDateNumber(2)
!
!  DATA component parallel distributed-memory domain partions in the
!  I- and J-directions (lon,lat).
!
      integer :: ItileD
      integer :: JtileD
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
!  Execution tracing flag:
!
!    [0] no tracing
!    [1] reports sequence of coupling subroutine calls
!    [2] <1> plus writes voluminous ESMF library tracing
!            information which slowdown performace, and
!            creates large log file
!
      integer :: TraceLevel = 0
!
!  Switch to trace/track run sequence during debugging.  All information
!  is written to Fortan unit trac. For now, use standard output unit.
!
      logical :: ESM_track = .FALSE.  ! trace/track CALL sequence switch
      integer :: trac = 6             ! trace/track CALL sequence unit
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
      integer, parameter :: Iconsv1 = 3     ! first-order conservative
      integer, parameter :: Iconsv2 = 4     ! second-order conservative
      integer, parameter :: InStoD  = 5     ! nearest neighbor Src 2 Dst
      integer, parameter :: InDtoS  = 6     ! nearest neighbor Dst 2 Src
!
      character (len=4), dimension(0:6) :: IntrpType =                  &
     &                                     (/ 'NONE',                   &
     &                                        'BLIN',                   &
     &                                        'PTCH',                   &
     &                                        'CNS1',                   &
     &                                        'CNS2',                   &
     &                                        'NS2D',                   &
     &                                        'ND2S' /)
!
!  Extrapolation method for unmapped destination points.
!
      integer, parameter :: Enone   = 0     ! none
      integer, parameter :: ExStoD  = 1     ! nearear neighbor Src 2 Dst
      integer, parameter :: Eidavg  = 2     ! inverse distance average
      integer, parameter :: Ecreep  = 3     ! creep fill
      integer, parameter :: E2steps = 4     ! Turuncoglu two steps
!
      character (len=4), dimension(0:4) :: ExtrpType =                  &
     &                                     (/ 'NONE',                   &
     &                                        'NS2D',                   &
     &                                        'IDAV',                   &
     &                                        'CREE',                   &
     &                                        '2STP' /)
!
!  The number of levels to output for the extrapolation methods that
!  fill levels, like creep fill (ESMF_EXTRAPMETHOD_CREEP). Unmapped
!  destination points are supplied by repeatedly moving data from
!  mapped locations to neighboring unmapped locations for a user-
!  specified number of levels.  For each creeped point, its value is
!  the average of the values of the immediate neighbors from the
!  mapped points from regridding (ESMF Reference Manual, v 8.0.0).
!
      integer :: extrapNumLevels = 1
!
!  Interpolation connectors mask interaction flags.
!
      integer, parameter :: OverLand  = 1
      integer, parameter :: OverOcean = 2
      integer, parameter :: OverAll   = 3
!
      character (len=3), dimension(3) :: MaskType =                     &
     &                                   (/ 'LND',                      &
     &                                      'OCN',                      &
     &                                      'ALL'/)
!
!  Coupling run mode: sequential or concurrent.
!
      integer, parameter :: Iseq = 1
      integer, parameter :: Ipar = 2
!
      character (len=10), dimension(2) :: RunMode =                     &
     &                                   (/ 'SEQUENTIAL',               &
     &                                      'CONCURRENT' /)
!
!  Compling standard input parameters filename.
!
      character (len=256) :: CinpName
!
!  ESM free-format run sequence configuration filename
!
      character (len=256) :: CONFname
!
!  Coupling Import/Export variable metadata filename.
!
      character (len=256) :: CPLname
!
!  Standard input filename for each coupled model, [Nmodels].
!
      character (len=256), allocatable :: INPname(:)
!
!  Standard output units coupler and log messages filename
!  for coupler and ESMF library.
!
      integer :: cplout  = 77             ! coupling driver
      integer :: dataout = 77             ! data component
!
      character (len= 8), parameter :: ESMnameLog = 'log.esmf'
      character (len=11), parameter :: CouplerLog = 'log.coupler'
!
!  Output NetCDF file used to store field snapshot attributes needed for
!  time interpolation by the ESM component kernel during concurrent
!  coupling.
!
      character (len=17), parameter :: AttFileName = 'time_intrp_att.nc'
!
!-----------------------------------------------------------------------
!  ESM constants
!-----------------------------------------------------------------------
!
      integer (i4b), parameter :: MAPPED_MASK = 99_i4b
      integer (i4b), parameter :: UNMAPPED_MASK = 98_i4b

      real (dp), parameter :: MISSING_dp = 1.0E20_dp
      real (r4), parameter :: MISSING_r4 = 1.0E20_r4
      real (r8), parameter :: MISSING_r8 = 1.0E20_r8

      real (dp), parameter :: TOL_dp = 0.5E20_dp
      real (r4), parameter :: TOL_r4 = 0.5E20_r4
      real (r8), parameter :: TOL_r8 = 0.5E20_r8
!
      CONTAINS
!
      SUBROUTINE allocate_esmf_esm
!
!=======================================================================
!                                                                      !
!  This routine allocates module coupling structures.                  !
!                                                                      !
!=======================================================================
!
!  Local variable definitions.
!
      integer :: i, j, ng

      character (len= 1), parameter :: blank = ' '
      character (len=50) :: MyLabel
!
!-----------------------------------------------------------------------
!  Allocate coupling structures.
!-----------------------------------------------------------------------
!
!  Allocate coupled models high-level data structure.
!
      IF (.not.allocated(MODELS)) THEN
        allocate ( MODELS(Nmodels) )
      END IF
      DO i=1,Nmodels
        MODELS(i) % IsActive = .FALSE.
        MODELS(i) % nPETs = 0
        MODELS(i) % Ngrids = 0
      END DO
!
!  Allocate coupled state set high-level structure.
!
      IF (.not.allocated(COUPLED)) THEN
        allocate ( COUPLED(Nmodels) )
      END IF

# ifdef DATA_COUPLING
!
!  Allocate melding weights structure.
!
      IF (.not.allocated(WEIGHTS)) THEN
        allocate ( WEIGHTS(Nmodels) )
      END IF
      DO i=1,Nmodels
        WEIGHTS(i) % VnameDATA = blank
        WEIGHTS(i) % VnameESM = blank
        WEIGHTS(i) % ncfile = blank
        WEIGHTS(i) % NestedGrid = 0
      END DO
# endif
!
!  Set counter for ESM component calls to export and import routines.
!  It is used to check if the ESM component has the two time-snapshots
!  to performs time interpolation of export and import fields.
!
      DO i=1,Nmodels
        MODELS(i) % ExportCalls = 0
        MODELS(i) % ImportCalls = 0
      END DO
!
!  Allocate coupled models connector structure.
!
      IF (.not.allocated(CONNECTORS)) THEN
        allocate ( CONNECTORS(Nmodels,Nmodels) )
      END IF
      DO j=1,Nmodels
        DO i=1,Nmodels
          CONNECTORS(i,j) % IsActive = .FALSE.
          CONNECTORS(i,j) % MaskInteraction = OverOcean
          CONNECTORS(i,j) % nPETs = 0
        END DO
      END DO
!
!  Allocate drivers and coupled models clock information structure.
!
      IF (.not.allocated(ClockInfo)) THEN
        allocate ( ClockInfo(0:Nmodels) )
      END IF
!
!  Allocate MPI communicator handle for each ESM component.
!
      IF (.not.allocated(ESMcomm)) THEN
        allocate ( ESMcomm(0:Nmodels) )
      END IF
      ESMcomm=0
!
!  Allocate number of import and export ESM fields.
!
      IF (.not.allocated(Nimport)) THEN
        allocate ( Nimport(Nmodels) )
      END IF
      Nimport=0
!
      IF (.not.allocated(Nexport)) THEN
        allocate ( Nexport(Nmodels) )
      END IF
      Nexport=0
!
!  Allocate DATA Model high-level structure.
!
      IF (.not.allocated(DataSet)) THEN
        allocate ( DataSet(Nmodels) )
      END IF
      DO i=1,Nmodels
        DataSet(i)%Nfields=0
        DataSet(i)%Nfiles=0
      END DO
!
!  Standard input filename for each ESM component.
!
      IF (.not.allocated(INPname)) THEN
        allocate ( INPname(Nmodels) )
      END IF
      DO i=1,Nmodels
        INPname(i)=blank
      END DO
!
!  Generic ESM component labels used in the CASE constructs.
!
      IF (.not.allocated(Cmodel)) THEN
        allocate ( Clabel(Nmodels) )
      END IF
      Clabel(Iroms  )='OCN'
      Clabel(Iatmos )='ATM'
      Clabel(Idata  )='DAT'
      Clabel(Iseaice)='ICE'
      Clabel(Iwave  )='WAV'
!
      IF (.not.allocated(Cmodel)) THEN
        allocate ( Cmodel(Nmodels) )
      END IF
      Cmodel(Iroms  )='ROMS'
# if defined COAMPS_COUPLING
      Cmodel(Iatmos )='COAMPS'
# elif defined REGCM_COUPLING
      Cmodel(Iatmos )='RegCM'
# elif defined WRF_COUPLING
      Cmodel(Iatmos )='WRF'
# else
      Cmodel(Iatmos )='ATMOS'
# endif
      Cmodel(Idata  )='DATA'
# if defined CICE_COUPLING
      Cmodel(Iseaice)='CICE'
# else
      Cmodel(Iseaice)='SEAICE'
# endif
# if defined WAM_COUPLING
      Cmodel(Iwave  )='WAM'
# else
      Cmodel(Iwave  )='WAVES'
# endif
!
!  Set ESM component name.
!
      DO i=1,Nmodels
        MODELS(i)%name=TRIM(Cmodel(i))
      END DO
!
      RETURN
      END SUBROUTINE allocate_esmf_esm
!
      SUBROUTINE clock_report (model, Icomp, localPET, source, rc)
!
!=======================================================================
!                                                                      !
!  Inquire ESM driver and component clock and report information.      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     model     ESM component ESMF object (ESMF_GridComp)              !
!     Icomp     Component index in ESM structures (integer)            !
!     localPET  Local Persistent Execution Thread (integer)            !
!     source    calling routine (string)                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     DataSet   Updata DATA Model high-level structure in module       !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)  :: Icomp, localPET
      integer, intent(out) :: rc
!
      character (len=*) :: source
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      logical :: IsActive
!
      integer :: i, is
      integer :: AlarmCount(0:Nmodels)

      integer (i8b) :: AdvanceCount(0:Nmodels)
!
      real (dp) :: scale
      real (dp) :: RunTimeStepCount(0:Nmodels)
      real (dp) :: TimeCurr(0:Nmodels)
      real (dp) :: TimeNext(0:Nmodels)
      real (dp) :: TimePrev(0:Nmodels)
      real (dp) :: TimeStart(0:Nmodels)
      real (dp) :: TimeStop(0:Nmodels)
      real (dp) :: Time_Dura(0:Nmodels)
      real (dp) :: Time_Step(0:Nmodels)
!
      character (len=22) :: ClockName(0:Nmodels)
      character (len=22) :: DirectionString(0:Nmodels)
      character (len=22) :: TimeCurrString(0:Nmodels)
      character (len=22) :: TimeNextString(0:Nmodels)
      character (len=22) :: TimePrevString(0:Nmodels)
      character (len=22) :: TimeStartString(0:Nmodels)
      character (len=22) :: TimeStopString(0:Nmodels)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", clock_report"
!
      TYPE (ESMF_Clock)          :: Clock(0:Nmodels)
      TYPE (ESMF_Direction_flag) :: Direction(0:Nmodels)
      TYPE (ESMF_TimeInterval)   :: RunDuration(0:Nmodels)
      TYPE (ESMF_TimeInterval)   :: TimeStep(0:Nmodels)
      TYPE (ESMF_Time)           :: CurrTime(0:Nmodels)
      TYPE (ESMF_Time)           :: PrevTime(0:Nmodels)
      TYPE (ESMF_Time)           :: StartTime(0:Nmodels)
      TYPE (ESMF_Time)           :: StopTime(0:Nmodels)
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Inquire driver and compones clocks.
!-----------------------------------------------------------------------
!
      scale=1.0_dp/86400.0_dp
!
      DO i=0,Nmodels
        IF (i.eq.0) THEN
          IsActive=.TRUE.                           ! ESM driver
        ELSE
          IsActive=MODELS(i)%IsActive               ! ESM components
        END IF
        IF (IsActive) THEN
          Clock(i)=ClockInfo(i)%Clock
          CALL ESMF_ClockGet (Clock(i),                                 &
     &                        timeStep = TimeStep(i),                   &
     &                        startTime = StartTime(i),                 &
     &                        stopTime = StopTime(i),                   &
     &                        runDuration = RunDuration(i),             &
     &                        runTimeStepCount = RunTimeStepCount(i),   &
     &                        currTime = CurrTime(i),                   &
     &                        prevTime = PrevTime(i),                   &
     &                        advanceCount = AdvanceCount(i),           &
     &                        alarmCount = AlarmCount(i),               &
     &                        direction = Direction(i),                 &
     &                        name = ClockName(i),                      &
     &                        rc = rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          IF (Direction(i).eq.ESMF_DIRECTION_FORWARD) THEN
            DirectionString(i)='FORWARD'
          ELSE
            DirectionString(i)='REVERSE'
          END IF
          ClockInfo(i)%AdvanceCount=AdvanceCount(i)
!
          CALL ESMF_TimeGet (StartTime(i),                              &
     &                       s_r8=TimeStart(i),                         &
     &                       timeString=TimeStartString(i),             &
     &                       rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          is=INDEX(TimeStartString(i), 'T')           ! remove 'T' in
          IF (is.gt.0) TimeStartString(i)(is:is)=' '  ! ISO 8601 format
!
          CALL ESMF_TimeGet (StopTime(i),                               &
     &                       s_r8=TimeStop(i),                          &
     &                       timeString=TimeStopString(i),              &
     &                       rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          is=INDEX(TimeStopString(i), 'T')            ! remove 'T' in
          IF (is.gt.0) TimeStopString(i)(is:is)=' '   ! ISO 8601 format
!
          CALL ESMF_TimeGet (PrevTime(i),                               &
     &                       s_r8=TimePrev(i),                          &
     &                       timeString=TimePrevString(i),              &
     &                       rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          is=INDEX(TimePrevString(i), 'T')            ! remove 'T' in
          IF (is.gt.0) TimePrevString(i)(is:is)=' '   ! ISO 8601 format
!
          CALL ESMF_TimeGet (CurrTime(i),                               &
     &                       s_r8=TimeCurr(i),                          &
     &                       timeString=TimeCurrString(i),              &
     &                       rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          is=INDEX(TimeCurrString(i), 'T')            ! remove 'T' in
          IF (is.gt.0) TimeCurrString(i)(is:is)=' '   ! ISO 8601 format
!
          CALL ESMF_TimeGet (CurrTime(i)+TimeStep(i),                   &
     &                       s_r8=TimeNext(i),                          &
     &                       timeString=TimeNextString(i),              &
     &                       rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
          is=INDEX(TimeNextString(i), 'T')           ! remove 'T' in
          IF (is.gt.0) TimeNextString(i)(is:is)=' '  ! ISO 8601 format
!
          CALL ESMF_TimeIntervalGet (TimeStep(i),                       &
     &                               s_r8=Time_Step(i),                 &
     &                               rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
          CALL ESMF_TimeIntervalGet (RunDuration(i),                    &
     &                               s_r8=Time_Dura(i),                 &
     &                               rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc,                          &
     &                           msg=ESMF_LOGERR_PASSTHRU,              &
     &                           line=__LINE__,                         &
     &                           file=MyFile)) THEN
            RETURN
          END IF
!
!  Report clock information
!
          IF ((DebugLevel.gt.0).and.(localPET.eq.0)) THEN
            WRITE (cplout,10)                                           &
     &        '       Clock Name: ', TRIM(clockName(i))//' in '//source,&
     &        '       Start Time: ', TRIM(TimeStartString(i)),          &
     &                             TimeStart(i)*scale,                  &
     &        '        Stop Time: ', TRIM(TimeStopString(i)),           &
     &                             TimeStop(i)*scale,                   &
     &        '        Prev Time: ', TRIM(TimePrevString(i)),           &
     &                             TimePrev(i)*scale,                   &
     &        '     Current Time: ', TRIM(TimeCurrString(i)),           &
     &                             TimeCurr(i)*scale,                   &
     &        '        Next Time: ', TRIM(TimeNextString(i)),           &
     &                             TimeNext(i)*scale,                   &
     &        '    Time Interval: ', Time_Step(i), Time_Step(i)*scale,  &
     &        '     Run Duration: ', Time_Dura(i), Time_Dura(i)*scale,  &
     &        '   Stepping Count: ', RunTimeStepCount(i),               &
     &        '    Advance Count: ', AdvanceCount(i),                   &
     &        'Current Direction: ', TRIM(DirectionString(i)),          &
     &        '      Alarm Count: ', AlarmCount(i)
          END IF
        END IF
      END DO
      IF (DebugLevel.gt.0) CALL my_flush (cplout)
!
!  Store the clock advance counter.
!
  10  FORMAT(/,4x,a,a,/,5(4x,a,a,2x,f15.8,/),2(4x,a,f19.8,2x,f15.8,/),  &
     &       4x,a,f19.8,/,4x,a,2x,i8,/,4x,a,a,/,4x,a,2x,i8,/)
!
      RETURN
      END SUBROUTINE clock_report
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
      INTEGER FUNCTION get_atm_Ngrids (Sname,localPET) RESULT (MyValue)
!
!=======================================================================
!                                                                      !
!  This function reads number of nested grids from atmosphere model    !
!  input script or namelist.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_scalars,  ONLY : NoError, exit_flag
!
      USE strings_mod,  ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: localPET
!
      character (len=*), intent(in) :: Sname
!
!  Local variable declarations.
!
      integer :: equal, ie, is, inp, io_err, out
!
      character (len= 40) :: keyword
      character (len=256) :: io_errmsg, line

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_atm_Ngrids"

#if defined COAMPS_COUPLING || defined WRF_COUPLING
!
!-----------------------------------------------------------------------
!  Read atmosphere model number of nested grids from input script or
!  namelist.
!-----------------------------------------------------------------------
!
      inp=2
      out=cplout
!
# if defined COAMPS_COUPLING
      keyword='nnest'             ! COAMPS "gridnl" namelist parameter
# elif defined WRF_COUPLING
      keyword='max_dom'           ! WRF "domains" namelist parameter
# endif
      OPEN (inp, FILE=TRIM(Sname), FORM='formatted', STATUS='old',      &
     &      IOSTAT=io_err, IOMSG=io_errmsg)
      IF (io_err.ne.0) THEN
        IF (localPET.eq.0) WRITE (out,30) TRIM(Sname), TRIM(io_errmsg)
        exit_flag=5
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Read script/namelist keyword parameter.
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        is=INDEX(line, TRIM(keyword))
        equal=INDEX(line, CHAR(61), BACK=.FALSE.)        ! equal sign
        IF ((is.ne.0).and.(equal.ne.0)) THEN
          is=equal+1
          ie=LEN_TRIM(line)
          READ (line(is:ie),*) MyValue
          IF (MyValue.le.0) THEN
            IF (localPET.eq.0) WRITE (out,40) TRIM(keyword), MyValue,   &
     &                         'must be greater than zero.'
            exit_flag=5
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
      END DO
  10  IF (localPET.eq.0) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)
!
  30  FORMAT (/,' get_atm_Ngrids - Unable to open ROMS/TOMS input ',    &
     &            'script file.',/,18x,a,/,18x,'ERROR: ',a)
  40  FORMAT (/,' get_atm_Ngrids - Invalid input parameter, ',a,i4,/,   &
     &        18x,a)
  50  FORMAT (/,' get_atm_Ngrids - Error while processing line: ',/,a)

#else
!
!-----------------------------------------------------------------------
!  Set atmosphere model nested grid parameter to unity.
!-----------------------------------------------------------------------
!
      MyValue=1
#endif
!
      RETURN
      END FUNCTION get_atm_Ngrids
!
      INTEGER FUNCTION get_roms_Ngrids (Sname,localPET) RESULT (MyValue)
!
!=======================================================================
!                                                                      !
!  This function reads number of nested grids (Ngrids) from ROMS input !
!  script.                                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param,     ONLY : Ngrids
      USE mod_scalars,   ONLY : NoError, exit_flag, time_ref
!
      USE inp_decode_mod
!
      USE dateclock_mod, ONLY : ref_clock
      USE strings_mod,   ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: localPET
!
      character (len=*), intent(in) :: Sname
!
!  Local variable declarations.
!
      integer :: Npts, Nval, inp, io_err, out, status
      integer :: Ivalue(1)
!
      real(r8) :: Rvalue(1)

      real(dp), dimension(nRval) :: Rval
!
      character (len= 40) :: KeyWord
      character (len=256) :: io_errmsg, line
      character (len=256), dimension(nCval) :: Cval

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_roms_Ngrids"
!
!-----------------------------------------------------------------------
!  Read ROMS application standard input file to determine the value
!  of Ngrids. Also, read time reference to check for compactability
!  with specified coupling value.
!-----------------------------------------------------------------------
!
      inp=2
      out=cplout
!
      OPEN (inp, FILE=TRIM(Sname), FORM='formatted', STATUS='old',      &
     &      IOSTAT=io_err, IOMSG=io_errmsg)
      IF (io_err.ne.0) THEN
        IF (localPET.eq.0) WRITE (out,30) TRIM(Sname), TRIM(io_errmsg)
        exit_flag=5
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Ngrids')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              MyValue=Ivalue(1)
              IF (MyValue.le.0) THEN
                IF (localPET.eq.0) WRITE (out,40) 'Ngrids', MyValue,    &
     &                             'must be greater than zero.'
                exit_flag=5
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
              END IF
              Ngrids=MyValue
            CASE ('TIME_REF')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              time_ref=Rvalue(1)
              CALL ref_clock (time_ref)
          END SELECT
        END IF
      END DO
  10  IF (localPET.eq.0) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)
!
  30  FORMAT (/,' get_roms_Ngrids - Unable to open ROMS/TOMS input ',   &
     &          'script file.',/,19x,a,/,19x,a)
  40  FORMAT (/,' get_roms_Ngrids - Invalid input parameter, ',a,i4,/,  &
     &        14x,a)
  50  FORMAT (/,' get_roms_Ngrids - Error while processing line: ',/,a)
!
      RETURN
      END FUNCTION get_roms_Ngrids
!
      SUBROUTINE load_IFS (S, Ifiles, Ngrids, Nfiles, Icomp)
!
!=======================================================================
!                                                                      !
!  This routine loads the DATA Model source files into Input File      !
!  Structure (IFS).                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     S         Input DATA Model source files as read from coupling    !
!                 script (TYPE T_IO)                                   !
!     Ifiles    First dimension of S-structure (integer)               !
!     Ngrids    Number of nested grids in ESM component (integer)      !
!     Nfiles    Number of files per nested grid (integer vector)       !
!     Icomp     ESM component index per neste grid (integer)           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     DataSet   Updata DATA Model high-level structure in module       !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)  :: Ifiles, Ngrids
      integer, intent(in)  :: Icomp, Nfiles(Ngrids)
!
      TYPE (T_IO), intent(in) :: S(Ifiles,Ngrids)
!
!  Local variable declarations.
!
      integer :: ic, if, j, ng, mf
!
!-----------------------------------------------------------------------
!  Load DATA Model source files into IFS.
!-----------------------------------------------------------------------
!
      ic=Icomp                                  ! ESM component index
      DO ng=1,Ngrids
        IF (MODELS(Icomp)%IsActive) THEN
          DO if=1,Nfiles(ng)
            mf=S(if,ng)%Nfiles                  ! number of multi-files
!
!  Allocate various variables dimensioned by the number of multi-files.
!
            allocate ( DataSet(ic)%IFS(if)%Nrec(mf) )
            allocate ( DataSet(ic)%IFS(if)%time_min(mf) )
            allocate ( DataSet(ic)%IFS(if)%time_max(mf) )
            allocate ( DataSet(ic)%IFS(if)%Vid(mf) )
            allocate ( DataSet(ic)%IFS(if)%Tid(mf) )
            allocate ( DataSet(ic)%IFS(if)%files(mf) )
!
!  Initialize and load fields into structure.
!
            DataSet(ic)%IFS(if)%Nfiles=mf       ! number of multi-files
            DataSet(ic)%IFS(if)%Fcount=1        ! multi-file counter
            DataSet(ic)%IFS(if)%Rindex=0        ! time index
            DataSet(ic)%IFS(if)%ncid=-1         ! closed NetCDF state
            DataSet(ic)%IFS(if)%Vid=-1          ! NetCDF variable ID
            DataSet(ic)%IFS(if)%Tid=-1          ! NetCDF time ID
            DataSet(ic)%IFS(if)%Nrec=0          ! number of time records
            DataSet(ic)%IFS(if)%time_min=0.0_r8 ! starting time
            DataSet(ic)%IFS(if)%time_max=0.0_r8 ! ending time
!
!  Enter multi-filenames.
!
            DO j=1,mf
              DataSet(ic)%IFS(if)%files(j)=TRIM(S(if,ng)%files(j))
            END DO
            DataSet(ic)%IFS(if)%label=TRIM(S(if,ng)%label)
            DataSet(ic)%IFS(if)%name=TRIM(S(if,ng)%name)
            DataSet(ic)%IFS(if)%base=TRIM(S(if,ng)%base)
          END DO
        END IF
      END DO
!
      RETURN
      END SUBROUTINE load_IFS
!
      SUBROUTINE read_ESMconfig (vm, rc)
!
!=======================================================================
!                                                                      !
!  This routine reads in the Earth System Models (ESM) coupling        !
!  configuration parameters from standard input file.                  !
!                                                                      !
!=======================================================================
!
      USE mod_parallel,   ONLY : OCN_COMM_WORLD, MyRank
      USE mod_parallel,   ONLY : allocate_parallel, initialize_parallel
      USE mod_iounits,    ONLY : Iname, SourceFile
      USE mod_scalars,    ONLY : NoError, Rclock, exit_flag
      USE mod_strings,    ONLY : my_cpu, my_fc, my_fflags, my_fort,     &
     &                           my_os, Rdir
!
      USE inp_decode_mod
!
      USE dateclock_mod,  ONLY : datenum
      USE dateclock_mod,  ONLY : get_date
      USE dateclock_mod,  ONLY : time_string
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasts
# endif
      USE strings_mod,    ONLY : lowercase
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_VM) :: vm
!
!  Local variable declarations.
!
      logical :: doit, first
      logical :: Lwrite
      logical :: Lvalue(1)
      logical, allocatable :: LvalueA(:), LvalueI(:)
      logical, allocatable :: LvalueR(:), LvalueW(:)
!
      integer :: Icomp, Jcomp, MyError, MySize, Nfields, Npts, Nval
      integer :: ESMcount, Igrids, Jgrids, Ncplsets, Ngrd, Nstates
      integer :: ExpLstr, ImpLstr, SetLstr, lstr
      integer :: i, ic, ie, ig, is, j, jg, k, inp, ng, out, status
      integer :: ifile, igrid, io_err
      integer :: maxD2A_files, maxD2I_files, maxD2R_files, maxD2W_files
      integer :: localPET, PETcount, MyComm, nPETs
      integer :: Cdim, Rdim, TimeFrac
      integer :: Ivalue(1), nD2A(1), nD2I(1), nD2R(1), nD2W(1)
      integer :: NewDate(7)
      integer :: ESMorder(Nmodels)
      integer :: location(1)

      integer, allocatable :: IvalueA(:)
      integer, allocatable :: IvalueI(:)
      integer, allocatable :: IvalueR(:)
      integer, allocatable :: IvalueW(:)

      integer, allocatable :: Ncount(:,:)
!
      real(r8) :: Rvalue(1)

      real(dp), dimension(nRval) :: Rval
!
      character (len=  1), parameter :: blank = ' '
      character (len=  5) :: Pstr, Pend
      character (len= 20) :: DateCalendar
      character (len= 20) :: TimeReferenceString
      character (len= 20) :: TimeRestartString
      character (len= 20) :: TimeStartString
      character (len= 20) :: TimeStopString
      character (len= 40) :: frmt
      character (len= 40) :: KeyWord
      character (len= 50) :: MyLabel, label
      character (len= 80) :: string
      character (len=256) :: fname, io_errmsg, line
      character (len=256) :: ATMname, ICEname, OCNname, WAVname
      character (len=256), dimension(nCval) :: Cval

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", read_ESMconfig"
!
      character (len=  6), allocatable  :: SetLabel(:)
      character (len= 13), allocatable  :: ExpLabel(:), ImpLabel(:)
      character (len= 22), dimension(2) :: Fcode
!
      TYPE (ESMF_CalKind_Flag) :: CalType
      TYPE (ESMF_TIME)         :: MyStartTime, MyRestartTime
!
      TYPE(T_IO), allocatable  :: D2R(:,:)  ! DATA -> ROMS source files
      TYPE(T_IO), allocatable  :: D2A(:,:)  ! DATA -> ATM  source files
      TYPE(T_IO), allocatable  :: D2I(:,:)  ! DATA -> ICE  source files
      TYPE(T_IO), allocatable  :: D2W(:,:)  ! DATA -> WAV  source files
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Query gridded component.
!-----------------------------------------------------------------------
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
!  Assign driver communicatior handle to ROMS temporarily. It is needed
!  since here are using ROMS generic NetCDF and distributed-memory
!  interface.
!
      OCN_COMM_WORLD=MyComm
!
!  Set the rank for each driver PET. In sequential and concurrent setup,
!  each couputational node has it unique PET associated with the driver.
!
      PETrank=localPET
!
!  Open standard output unit for ESM coupler information and messages.
!  It is advisable to have such infomation separated from other standard
!  output from the coupled models.
!
      IF (localPET.eq.0) THEN
        OPEN (cplout, FILE=TRIM(CouplerLog), FORM='formatted',          &
     &        STATUS='replace')
      END IF
!
!-----------------------------------------------------------------------
!  Determine coupling standard input filename.  In distributed-memory,
!  this name is assigned at the executtion command line and processed
!  with the Unix routine GETARG.  The ROMS input parameter script name
!  is specified in this coupling script.
!-----------------------------------------------------------------------
!
      Lwrite=localPET .eq. 0
      inp=100
      out=cplout
!
      IF (localPET.eq.0) THEN
        CALL my_getarg (1, CinpName)
        CALL get_date (TodayDateString)
      END IF
      CALL mp_bcasts (1, 1, CinpName)
      CALL mp_bcasts (1, 1, TodayDateString)
      OPEN (inp, FILE=TRIM(CinpName), FORM='formatted', STATUS='old',   &
     &      IOSTAT=io_err, IOMSG=io_errmsg)
      IF (io_err.ne.0) THEN
        IF (localPET.eq.0) WRITE (out,10) TRIM(io_errmsg)
        exit_flag=5
        rc=ESMF_RC_FILE_OPEN
        RETURN
 10     FORMAT (/,' read_ESMconfig - Unable to open coupling input',    &
     &          ' script.',/,19x,'ERROR: ',a,/,                         &
     &          /,19x,'In distributed-memory applications, the input',  &
     &          /,19x,'script file is processed in parallel. The Unix', &
     &          /,19x,'routine GETARG is used to get script filename.', &
     &          /,19x,'For example, in MPI applications make sure that',&
     &          /,19x,'command line is something like:',/,              &
     &          /,19x,'mpirun -np 4 romsM coupling_esmf.in',/,          &
     &          /,19x,'and not',/,                                      &
     &          /,19x,'mpirun -np 4 romsM < coupling_esmf.in',/)
      END IF
!
!-----------------------------------------------------------------------
!  Read in multiple models coupling parameters. Then, load input
!  data into module. Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      ifile=1                            ! multiple file counter
      igrid=1                            ! nested grid counter
      DO i=1,LEN(label)
        label(i:i)=blank
      END DO
      Cdim=SIZE(Cval,1)
      Rdim=SIZE(Rval,1)
!
!  Read in coupling parameters.
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('CPLname')
              CPLname=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('CONFname')
              CONFname=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('INPname(roms)')
              OCNname=TRIM(ADJUSTL(Cval(Nval)))
              NgridsR=get_roms_Ngrids(OCNname,localPET)
              IF (exit_flag.ne.NoError) THEN
                rc=ESMF_RC_FILE_READ
                RETURN
              END IF
              IF (.not.allocated(LvalueR)) THEN
                allocate ( LvalueR(NgridsR) )
                LvalueR(1:NgridsR)=.FALSE.
              END IF
              IF (.not.allocated(IvalueR)) THEN
                allocate ( IvalueR(NgridsR) )
                IvalueR(1:NgridsR)=0
              END IF
            CASE ('INPname(atmos)')
              ATMname=TRIM(ADJUSTL(Cval(Nval)))
              NgridsA=get_atm_Ngrids(ATMname,localPET)
              IF (exit_flag.ne.NoError) THEN
                rc=ESMF_RC_FILE_READ
                RETURN
              END IF
              IF (.not.allocated(LvalueA)) THEN
                allocate ( LvalueA(NgridsA) )
                LvalueA(1:NgridsA)=.FALSE.
              END IF
              IF (.not.allocated(IvalueA)) THEN
                allocate ( IvalueA(NgridsA) )
                IvalueA(1:NgridsA)=0
              END IF
            CASE ('INPname(seaice)')
              ICEname=TRIM(ADJUSTL(Cval(Nval)))
              IF (.not.allocated(LvalueI)) THEN
                allocate ( LvalueI(NgridsI) )
                LvalueI(1:NgridsI)=.FALSE.
              END IF
              IF (.not.allocated(IvalueI)) THEN
                allocate ( IvalueI(NgridsI) )
                IvalueI(1:NgridsI)=0
              END IF
            CASE ('INPname(waves)')
              WAVname=TRIM(ADJUSTL(Cval(Nval)))
              IF (.not.allocated(LvalueW)) THEN
                allocate ( LvalueW(NgridsW) )
                LvalueW(1:NgridsW)=.FALSE.
              END IF
              IF (.not.allocated(IvalueW)) THEN
                allocate ( IvalueW(NgridsW) )
                IvalueW(1:NgridsW)=0
              END IF
            CASE ('IsActive(roms)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              CALL allocate_esmf_esm
              ESMcomm(Idriver)=MyComm
              Iname=TRIM(ADJUSTL(OCNname))              ! needed by ROMS
              INPname(Iroms)=TRIM(ADJUSTL(OCNname))
              MODELS(Iroms)%IsActive=Lvalue(1)
              IF (MODELS(Iroms)%IsActive) THEN
                MODELS(Iroms)%Ngrids=NgridsR
                IF (.not.allocated(MODELS(Iroms)%ExportState)) THEN
                  allocate ( MODELS(Iroms)%ExportState(NgridsR) )
                END IF
                IF (.not.allocated(MODELS(Iroms)%ImportState)) THEN
                  allocate ( MODELS(Iroms)%ImportState(NgridsR) )
                END IF
                IF (.not.allocated(COUPLED(Iroms)%LinkedGrid)) THEN
                  allocate ( COUPLED(Iroms)%LinkedGrid(NgridsR,         &
     &                                                 Nmodels) )
                  COUPLED(Iroms)%LinkedGrid=.FALSE.
                END IF
                IF (.not.allocated(MODELS(Iroms)%TimeFrac)) THEN
                  allocate ( MODELS(Iroms)%TimeFrac(NgridsR,Nmodels) )
                  MODELS(Iroms)%TimeFrac=1
                END IF
              END IF
!
!             In some setups, the information below is needed before
!             ROMS initialization (if active or not) because of calls
!             from other components to the distribute and netcdf
!             modules.
!
              CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
              CALL mpi_comm_size (OCN_COMM_WORLD, MySize, MyError)
              CALL allocate_parallel (NgridsR)
              CALL initialize_parallel
            CASE ('IsActive(atmos)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              INPname(Iatmos)=TRIM(ADJUSTL(ATMname))
              MODELS(Iatmos)%IsActive=Lvalue(1)
              IF (MODELS(Iatmos)%IsActive) THEN
                MODELS(Iatmos)%Ngrids=NgridsA
                IF (.not.allocated(MODELS(Iatmos)%ExportState)) THEN
                  allocate ( MODELS(Iatmos)%ExportState(NgridsA) )
                END IF
                IF (.not.allocated(MODELS(Iatmos)%ImportState)) THEN
                  allocate ( MODELS(Iatmos)%ImportState(NgridsA) )
                END IF
                IF (.not.allocated(COUPLED(Iatmos)%LinkedGrid)) THEN
                  allocate ( COUPLED(Iatmos)%LinkedGrid(NgridsA,        &
     &                                                  Nmodels) )
                  COUPLED(Iatmos)%LinkedGrid=.FALSE.
                END IF
                IF (.not.allocated(MODELS(Iatmos)%TimeFrac)) THEN
                  allocate ( MODELS(Iatmos)%TimeFrac(NgridsA,Nmodels) )
                  MODELS(Iatmos)%TimeFrac=1
                END IF
              END IF
            CASE ('IsActive(seaice)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              InpName(Iseaice)=TRIM(ADJUSTL(ICEname))
              MODELS(Iseaice)%IsActive=Lvalue(1)
              IF (MODELS(Iseaice)%IsActive) THEN
                MODELS(Iseaice)%Ngrids=NgridsI
                IF (.not.allocated(MODELS(Iseaice)%ExportState)) THEN
                  allocate ( MODELS(Iseaice)%ExportState(NgridsI) )
                END IF
                IF (.not.allocated(MODELS(Iseaice)%ImportState)) THEN
                  allocate ( MODELS(Iseaice)%ImportState(NgridsI) )
                END IF
                IF (.not.allocated(COUPLED(Iseaice)%LinkedGrid)) THEN
                  allocate( COUPLED(Iseaice)%LinkedGrid(NgridsI,        &
     &                                                  Nmodels) )
                  COUPLED(Iseaice)%LinkedGrid=.FALSE.
                END IF
                IF (.not.allocated(MODELS(Iseaice)%TimeFrac)) THEN
                  allocate ( MODELS(Iseaice)%TimeFrac(NgridsI,Nmodels) )
                  MODELS(Iseaice)%TimeFrac=1
                END IF
              END IF
            CASE ('IsActive(waves)')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              InpName(Iwave)=TRIM(ADJUSTL(WAVname))
              MODELS(Iwave)%IsActive=Lvalue(1)
              IF (MODELS(Iwave)%IsActive) THEN
                MODELS(Iwave)%Ngrids=NgridsW
                IF (.not.allocated(MODELS(Iwave)%ExportState)) THEN
                  allocate ( MODELS(Iwave)%ExportState(NgridsW) )
                END IF
                IF (.not.allocated(MODELS(Iwave)%ImportState)) THEN
                  allocate ( MODELS(Iwave)%ImportState(NgridsW) )
                END IF
                IF (.not.allocated(COUPLED(Iwave)%LinkedGrid)) THEN
                  allocate ( COUPLED(Iwave)%LinkedGrid(NgridsW,         &
     &                                                 Nmodels) )
                  COUPLED(Iwave)%LinkedGrid=.FALSE.
                END IF
                IF (.not.allocated(MODELS(Iwave)%TimeFrac)) THEN
                  allocate ( MODELS(Iwave)%TimeFrac(NgridsW,Nmodels) )
                  MODELS(Iwave)%TimeFrac=1
                END IF
              END IF
            CASE ('IsActive(data)')
# ifdef DATA_COUPLING
              Npts=load_l(Nval, Cval, 1, Lvalue)
              MODELS(Idata)%IsActive=Lvalue(1)
              IF (MODELS(Idata)%IsActive) THEN
                MODELS(Idata)%Ngrids=NgridsD
                IF (.not.allocated(MODELS(Idata)%TimeFrac)) THEN
                  allocate ( MODELS(Idata)%TimeFrac(NgridsD,Nmodels) )
                  MODELS(Idata)%TimeFrac=1
                END IF
              END IF
# endif
            CASE ('Coupled(ATM2OCN)')                   ! ESM to OCN
              Npts=load_l(Nval, Cval, NgridsR, LvalueR)
              IF (MODELS(Iatmos)%IsActive.and.                          &
     &            MODELS(Iroms )%IsActive) THEN
                DO ng=1,NgridsR
                  COUPLED(Iroms)%LinkedGrid(ng,Iatmos)=LvalueR(ng)
                END DO
                IF (ANY(LvalueR)) THEN
                  CONNECTORS(Iatmos,Iroms)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(ICE2OCN)')
              Npts=load_l(Nval, Cval, NgridsR, LvalueR)
              IF (MODELS(Iseaice)%IsActive.and.                         &
     &            MODELS(Iroms  )%IsActive) THEN
                DO ng=1,NgridsR
                  COUPLED(Iroms)%LinkedGrid(ng,Iseaice)=LvalueR(ng)
                END DO
                IF (ANY(LvalueR)) THEN
                  CONNECTORS(Iseaice,Iroms)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(WAV2OCN)')
              Npts=load_l(Nval, Cval, NgridsR, LvalueR)
              IF (MODELS(Iwave)%IsActive.and.                           &
     &            MODELS(Iroms)%IsActive) THEN
                DO ng=1,NgridsR
                  COUPLED(Iroms)%LinkedGrid(ng,Iwave)=LvalueR(ng)
                END DO
                IF (ANY(LvalueR)) THEN
                  CONNECTORS(Iwave,Iroms)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(DAT2OCN)')
              Npts=load_l(Nval, Cval, NgridsR, LvalueR)
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Iroms)%IsActive) THEN
                DO ng=1,NgridsR
                  COUPLED(Iroms)%LinkedGrid(ng,Idata)=LvalueR(ng)
                END DO
                IF (ANY(LvalueR)) THEN
                  CONNECTORS(Idata,Iroms)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(OCN2ATM)')                   ! ESM to ATM
              Npts=load_l(Nval, Cval, NgridsA, LvalueA)
              IF (MODELS(Iroms )%IsActive.and.                          &
     &            MODELS(Iatmos)%IsActive) THEN
                DO ng=1,NgridsA
                  COUPLED(Iatmos)%LinkedGrid(ng,Iroms)=LvalueA(ng)
                END DO
                IF (ANY(LvalueA)) THEN
                  CONNECTORS(Iroms,Iatmos)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(ICE2ATM)')
              Npts=load_l(Nval, Cval, NgridsA, LvalueA)
              IF (MODELS(Iseaice)%IsActive.and.                         &
     &            MODELS(Iatmos )%IsActive) THEN
                DO ng=1,NgridsA
                  COUPLED(Iatmos)%LinkedGrid(ng,Iseaice)=LvalueA(ng)
                END DO
                IF (ANY(LvalueA)) THEN
                  CONNECTORS(Iseaice,Iatmos)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(WAV2ATM)')
              Npts=load_l(Nval, Cval, NgridsA, LvalueA)
              IF (MODELS(Iwave )%IsActive.and.                          &
     &            MODELS(Iatmos)%IsActive) THEN
                DO ng=1,NgridsA
                  COUPLED(Iatmos)%LinkedGrid(ng,Iwave)=LvalueA(ng)
                END DO
                IF (ANY(LvalueA)) THEN
                  CONNECTORS(Iwave,Iatmos)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(DAT2ATM)')
              Npts=load_l(Nval, Cval, NgridsA, LvalueA)
              IF (MODELS(Idata )%IsActive.and.                          &
     &            MODELS(Iatmos)%IsActive) THEN
                DO ng=1,NgridsA
                  COUPLED(Iatmos)%LinkedGrid(ng,Idata)=LvalueA(ng)
                END DO
                IF (ANY(LvalueA)) THEN
                  CONNECTORS(Idata,Iatmos)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(ATM2ICE)')                   ! ESM to ICE
              Npts=load_l(Nval, Cval, NgridsI, LvalueI)
              IF (MODELS(Iatmos )%IsActive.and.                         &
     &            MODELS(Iseaice)%IsActive) THEN
                DO ng=1,NgridsI
                  COUPLED(Iseaice)%LinkedGrid(ng,Iatmos)=LvalueI(ng)
                END DO
                IF (ANY(LvalueI)) THEN
                  CONNECTORS(Iatmos,Iseaice)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(OCN2ICE)')
              Npts=load_l(Nval, Cval, NgridsI, LvalueI)
              IF (MODELS(Iroms  )%IsActive.and.                         &
     &            MODELS(Iseaice)%IsActive) THEN
                DO ng=1,NgridsI
                  COUPLED(Iseaice)%LinkedGrid(ng,Iroms)=LvalueI(ng)
                END DO
                IF (ANY(LvalueI)) THEN
                  CONNECTORS(Iroms,Iseaice)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(WAV2ICE)')
              Npts=load_l(Nval, Cval, NgridsI, LvalueI)
              IF (MODELS(Iwave  )%IsActive.and.                         &
     &            MODELS(Iseaice)%IsActive) THEN
                DO ng=1,NgridsI
                  COUPLED(Iseaice)%LinkedGrid(ng,Iwave)=LvalueI(ng)
                END DO
                IF (ANY(LvalueI)) THEN
                  CONNECTORS(Iwave,Iseaice)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(DAT2ICE)')
              Npts=load_l(Nval, Cval, NgridsI, LvalueI)
              IF (MODELS(Idata  )%IsActive.and.                         &
     &            MODELS(Iseaice)%IsActive) THEN
                DO ng=1,NgridsI
                  COUPLED(Iseaice)%LinkedGrid(ng,Idata)=LvalueI(ng)
                END DO
                IF (ANY(LvalueI)) THEN
                  CONNECTORS(Idata,Iseaice)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(ATM2WAV)')                   ! ESM to WAV
              Npts=load_l(Nval, Cval, NgridsW, LvalueW)
              IF (MODELS(Iatmos)%IsActive.and.                          &
     &            MODELS(Iwave )%IsActive) THEN
                DO ng=1,NgridsW
                  COUPLED(Iwave)%LinkedGrid(ng,Iatmos)=LvalueW(ng)
                END DO
                IF (ANY(LvalueW)) THEN
                  CONNECTORS(Iatmos,Iwave)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(ICE2WAV)')
              Npts=load_l(Nval, Cval, NgridsW, LvalueW)
              IF (MODELS(Iseaice)%IsActive.and.                         &
     &            MODELS(Iwave  )%IsActive) THEN
                DO ng=1,NgridsW
                  COUPLED(Iwave)%LinkedGrid(ng,Iseaice)=LvalueW(ng)
                END DO
                IF (ANY(LvalueW)) THEN
                  CONNECTORS(Iseaice,Iwave)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(OCN2WAV)')
              Npts=load_l(Nval, Cval, NgridsW, LvalueW)
              IF (MODELS(Iroms)%IsActive.and.                           &
     &            MODELS(Iwave)%IsActive) THEN
                DO ng=1,NgridsW
                  COUPLED(Iwave)%LinkedGrid(ng,Iroms)=LvalueW(ng)
                END DO
                IF (ANY(LvalueW)) THEN
                  CONNECTORS(Iroms,Iwave)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('Coupled(DAT2WAV)')
              Npts=load_l(Nval, Cval, NgridsW, LvalueW)
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Iwave)%IsActive) THEN
                DO ng=1,NgridsW
                  COUPLED(Iwave)%LinkedGrid(ng,Idata)=LvalueW(ng)
                END DO
                IF (ANY(LvalueW)) THEN
                  CONNECTORS(Idata,Iwave)%IsActive=.TRUE.
                END IF
              END IF
            CASE ('CouplingType')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              CouplingType=Ivalue(1)
            CASE ('PETlayoutOption')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (Ivalue(1).eq.0) THEN
                PETlayoutOption='SEQUENTIAL'
              ELSE IF (Ivalue(1).eq.1) THEN
                PETlayoutOption='CONCURRENT'
              END IF
            CASE ('ItileD')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              ItileD=Ivalue(1)
            CASE ('JtileD')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              JtileD=Ivalue(1)
              IF (MODELS(Idata)%IsActive) THEN
                MODELS(Idata)%nPETs=ItileD*JtileD
                IF (.not.allocated(MODELS(Idata)%PETlist)) THEN
                  allocate ( MODELS(Idata)%PETlist(ItileD*JtileD) )
                END IF
                MODELS(Idata)%PETlist=0
              ELSE
                MODELS(Idata)%nPETs=0
              END IF
            CASE ('Nthreads(roms)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (MODELS(Iroms)%IsActive) THEN
                MODELS(Iroms)%nPETs=Ivalue(1)
                IF (.not.allocated(MODELS(Iroms)%PETlist)) THEN
                  allocate ( MODELS(Iroms)%PETlist(Ivalue(1)) )
                END IF
                MODELS(Iroms)%PETlist=0
              ELSE
                MODELS(Iroms)%nPETs=0
              END IF
            CASE ('Nthreads(atmos)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (MODELS(Iatmos)%IsActive) THEN
                MODELS(Iatmos)%nPETs=Ivalue(1)
                IF (.not.allocated(MODELS(Iatmos)%PETlist)) THEN
                  allocate ( MODELS(Iatmos)%PETlist(Ivalue(1)) )
                END IF
                MODELS(Iatmos)%PETlist=0
              ELSE
                MODELS(Iatmos)%nPETs=0
              END IF
            CASE ('Nthreads(seaice)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (MODELS(Iseaice)%IsActive) THEN
                MODELS(Iseaice)%nPETs=Ivalue(1)
                IF (.not.allocated(MODELS(Iseaice)%PETlist)) THEN
                  allocate ( MODELS(Iseaice)%PETlist(Ivalue(1)) )
                END IF
                MODELS(Iseaice)%PETlist=0
              ELSE
                MODELS(Iseaice)%nPETs=0
              END IF
            CASE ('Nthreads(waves)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (MODELS(Iwave)%IsActive) THEN
                MODELS(Iwave)%nPETs=Ivalue(1)
                IF (.not.allocated(MODELS(Iwave)%PETlist)) THEN
                  allocate ( MODELS(Iwave)%PETlist(Ivalue(1)) )
                END IF
                MODELS(Iwave)%PETlist=0
              ELSE
                MODELS(Iwave)%nPETs=0
              END IF
            CASE ('ReferenceTime')
              Npts=load_i(Nval, Rval, 6, ReferenceDate)
            CASE ('StartTime')
              Npts=load_i(Nval, Rval, 6, StartDate)
            CASE ('RestartTime')
              Npts=load_i(Nval, Rval, 6, RestartDate)
            CASE ('StopTime')
              Npts=load_i(Nval, Rval, 6, StopDate)
            CASE ('TimeStep')
              Npts=load_i(Nval, Rval, 6, TimeStep)
            CASE ('Calendar')
              DateCalendar=TRIM(ADJUSTL(Cval(Nval)))
              DO i=0,Nmodels
                ClockInfo(i)%CalendarString=TRIM(DateCalendar)
              END DO
            CASE ('TimeFrac(ATM2OCN)')                  ! ESM to OCN
              Npts=load_i(Nval, Rval, NgridsR, IvalueR)
              IF (MODELS(Iatmos)%IsActive.and.                          &
     &            MODELS(Iroms )%IsActive) THEN
                DO ng=1,NgridsR
                  MODELS(Iroms)%TimeFrac(ng,Iatmos)=IvalueR(ng)
                END DO
              END IF
            CASE ('TimeFrac(ICE2OCN)')
              Npts=load_i(Nval, Rval, NgridsR, IvalueR)
              IF (MODELS(Iseaice)%IsActive.and.                         &
     &            MODELS(Iroms  )%IsActive) THEN
                DO ng=1,NgridsR
                  MODELS(Iroms)%TimeFrac(ng,Iseaice)=IvalueR(ng)
                END DO
              END IF
            CASE ('TimeFrac(WAV2OCN)')
              Npts=load_i(Nval, Rval, NgridsR, IvalueR)
              IF (MODELS(Iwave)%IsActive.and.                           &
     &            MODELS(Iroms)%IsActive) THEN
                DO ng=1,NgridsR
                  MODELS(Iroms)%TimeFrac(ng,Iwave)=IvalueR(ng)
                END DO
              END IF
            CASE ('TimeFrac(DAT2OCN)')
              Npts=load_i(Nval, Rval, NgridsR, IvalueR)
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Iroms)%IsActive) THEN
                DO ng=1,NgridsR
                  MODELS(Iroms)%TimeFrac(ng,Idata)=IvalueR(ng)
                END DO
              END IF
            CASE ('TimeFrac(OCN2ATM)')                  ! ESM to ATM
              Npts=load_i(Nval, Rval, NgridsA, IvalueA)
              IF (MODELS(Iroms )%IsActive.and.                          &
     &            MODELS(Iatmos)%IsActive) THEN
                DO ng=1,NgridsA
                  MODELS(Iatmos)%TimeFrac(ng,Iroms)=IvalueA(ng)
                END DO
              END IF
            CASE ('TimeFrac(ICE2ATM)')
              Npts=load_i(Nval, Rval, NgridsA, IvalueA)
              IF (MODELS(Iseaice)%IsActive.and.                         &
     &            MODELS(Iatmos )%IsActive) THEN
                DO ng=1,NgridsA
                  MODELS(Iatmos)%TimeFrac(ng,Iseaice)=IvalueA(ng)
                END DO
              END IF
            CASE ('TimeFrac(WAV2ATM)')
              Npts=load_i(Nval, Rval, NgridsA, IvalueA)
              IF (MODELS(Iwave )%IsActive.and.                          &
     &            MODELS(Iatmos)%IsActive) THEN
                DO ng=1,NgridsA
                  MODELS(Iatmos)%TimeFrac(ng,Iwave)=IvalueA(ng)
                END DO
              END IF
            CASE ('TimeFrac(DAT2ATM)')
              Npts=load_i(Nval, Rval, NgridsA, IvalueA)
              IF (MODELS(Idata )%IsActive.and.                          &
     &            MODELS(Iatmos)%IsActive) THEN
                DO ng=1,NgridsA
                  MODELS(Iatmos)%TimeFrac(ng,Idata)=IvalueA(ng)
                END DO
              END IF
            CASE ('TimeFrac(ATM2ICE)')                  ! ESM to ICE
              Npts=load_i(Nval, Rval, NgridsI, IvalueI)
              IF (MODELS(Iatmos )%IsActive.and.                         &
     &            MODELS(Iseaice)%IsActive) THEN
                DO ng=1,NgridsI
                  MODELS(Iseaice)%TimeFrac(ng,Iatmos)=IvalueI(ng)
                END DO
              END IF
            CASE ('TimeFrac(OCN2ICE)')
              Npts=load_i(Nval, Rval, NgridsI, IvalueI)
              IF (MODELS(Iroms  )%IsActive.and.                         &
     &            MODELS(Iseaice)%IsActive) THEN
                DO ng=1,NgridsI
                  MODELS(Iseaice)%TimeFrac(ng,Iroms)=IvalueI(ng)
                END DO
              END IF
            CASE ('TimeFrac(WAV2ICE)')
              Npts=load_i(Nval, Rval, NgridsI, IvalueI)
              IF (MODELS(Iwave  )%IsActive.and.                         &
     &            MODELS(Iseaice)%IsActive) THEN
                DO ng=1,NgridsI
                  MODELS(Iseaice)%TimeFrac(ng,Iwave)=IvalueI(ng)
                END DO
              END IF
            CASE ('TimeFrac(DAT2ICE)')
              Npts=load_i(Nval, Rval, NgridsI, IvalueI)
              IF (MODELS(Idata  )%IsActive.and.                         &
     &            MODELS(Iseaice)%IsActive) THEN
                DO ng=1,NgridsI
                  MODELS(Iseaice)%TimeFrac(ng,Idata)=IvalueI(ng)
                END DO
              END IF
            CASE ('TimeFrac(ATM2WAV)')                  ! ESM to WAV
              Npts=load_i(Nval, Rval, NgridsW, IvalueW)
              IF (MODELS(Iatmos)%IsActive.and.                          &
     &            MODELS(Iwave )%IsActive) THEN
                DO ng=1,NgridsW
                  MODELS(Iwave)%TimeFrac(ng,Iatmos)=IvalueW(ng)
                END DO
              END IF
            CASE ('TimeFrac(ICE2WAV)')
              Npts=load_i(Nval, Rval, NgridsW, IvalueW)
              IF (MODELS(Iseaice)%IsActive.and.                         &
     &            MODELS(Iwave  )%IsActive) THEN
                DO ng=1,NgridsW
                  MODELS(Iwave)%TimeFrac(ng,Iseaice)=IvalueW(ng)
                END DO
              END IF
            CASE ('TimeFrac(OCN2WAV)')
              Npts=load_i(Nval, Rval, NgridsW, IvalueW)
              IF (MODELS(Iroms)%IsActive.and.                           &
     &            MODELS(Iwave)%IsActive) THEN
                DO ng=1,NgridsW
                  MODELS(Iwave)%TimeFrac(ng,Iroms)=IvalueW(ng)
                END DO
              END IF
            CASE ('TimeFrac(DAT2WAV)')
              Npts=load_i(Nval, Rval, NgridsW, IvalueW)
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Iwave)%IsActive) THEN
                DO ng=1,NgridsW
                  MODELS(Iwave)%TimeFrac(ng,Idata)=IvalueW(ng)
                END DO
              END IF
            CASE ('extrapNumLevels')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              extrapNumLevels=Ivalue(1)
# ifdef DATA_COUPLING
            CASE ('WeightsFile(atmos)')
              IF ((MODELS(Iatmos)%IsActive).and.                        &
     &            (MODELS(Idata )%IsActive)) THEN
                WEIGHTS(Iatmos)%ncfile=TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('VnameDATA(atmos)')
              IF ((MODELS(Iatmos)%IsActive).and.                        &
     &            (MODELS(Idata )%IsActive)) THEN
                WEIGHTS(Iatmos)%VnameDATA=TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('VnameESM(atmos)')
              IF ((MODELS(Iatmos)%IsActive).and.                        &
     &            (MODELS(Idata )%IsActive)) THEN
                WEIGHTS(Iatmos)%VnameESM=TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('NestedGrid(atmos)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF ((MODELS(Iatmos)%IsActive).and.                        &
     &            (MODELS(Idata )%IsActive)) THEN
                WEIGHTS(Iatmos)%NestedGrid=Ivalue(1)
              END IF
# endif
            CASE ('DebugLevel')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              DebugLevel=Ivalue(1)
            CASE ('TraceLevel')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              TraceLevel=Ivalue(1)
              IF (TraceLevel.gt.0) THEN
                ESM_track=.TRUE.
              END IF
            CASE ('Nimport(roms)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (Ivalue(1).gt.0) THEN
                IF (.not.allocated(MODELS(Iroms)%ImportField)) THEN
                  allocate ( MODELS(Iroms)%ImportField(Ivalue(1)) )
                END IF
                Nimport(Iroms)=Ivalue(1)
              END IF
            CASE ('Nexport(roms)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (Ivalue(1).gt.0) THEN
                IF (.not.allocated(MODELS(Iroms)%ExportField)) THEN
                  allocate ( MODELS(Iroms)%ExportField(Ivalue(1)) )
                END IF
                Nexport(Iroms)=Ivalue(1)
              END IF
            CASE ('Import(roms)')
              IF ((Nimport(Iroms).gt.0).and.                            &
     &            (Nval.le.Nimport(Iroms))) THEN
                MODELS(Iroms)%ImportField(Nval)%short_name=             &
     &                                   TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('Export(roms)')
              IF ((Nexport(Iroms).gt.0).and.                            &
     &            (Nval.le.Nexport(Iroms))) THEN
                MODELS(Iroms)%ExportField(Nval)%short_name=             &
     &                                   TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('Nimport(atmos)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (Ivalue(1).gt.0) THEN
                IF (MODELS(Iatmos)%IsActive) THEN
                  IF (.not.allocated(MODELS(Iatmos)%ImportField)) THEN
                    allocate ( MODELS(Iatmos)%ImportField(Ivalue(1)) )
                  END IF
                  Nimport(Iatmos)=Ivalue(1)
                END IF
              END IF
            CASE ('Nexport(atmos)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (Ivalue(1).gt.0) THEN
                IF (MODELS(Iatmos)%IsActive) THEN
                  IF (.not.allocated(MODELS(Iatmos)%ExportField)) THEN
                    allocate ( MODELS(Iatmos)%ExportField(Ivalue(1)) )
                  END IF
                  Nexport(Iatmos)=Ivalue(1)
                END IF
              END IF
            CASE ('Import(atmos)')
              IF (MODELS(Iatmos)%IsActive) THEN
                IF ((Nimport(Iatmos).gt.0).and.                         &
     &              (Nval.le.Nimport(Iatmos))) THEN
                  MODELS(Iatmos)%ImportField(Nval)%short_name=          &
     &                                      TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('Export(atmos)')
              IF (MODELS(Iatmos)%IsActive) THEN
                IF ((Nexport(Iatmos).gt.0).and.                         &
     &              (Nval.le.Nexport(Iatmos))) THEN
                  MODELS(Iatmos)%ExportField(Nval)%short_name=          &
     &                                      TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('Nimport(seaice)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (Ivalue(1).gt.0) THEN
                IF (MODELS(Iseaice)%IsActive) THEN
                  IF (.not.allocated(MODELS(Iseaice)%ImportField)) THEN
                    allocate ( MODELS(Iseaice)%ImportField(Ivalue(1)) )
                  END IF
                  Nimport(Iseaice)=Ivalue(1)
                END IF
              END IF
            CASE ('Nexport(seaice)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (Ivalue(1).gt.0) THEN
                IF (MODELS(Iseaice)%IsActive) THEN
                  IF (.not.allocated(MODELS(Iseaice)%ExportField)) THEN
                    allocate ( MODELS(Iseaice)%ExportField(Ivalue(1)) )
                  END IF
                  Nexport(Iseaice)=Ivalue(1)
                END IF
              END IF
            CASE ('Import(seaice)')
              IF (MODELS(Iseaice)%IsActive) THEN
                IF ((Nimport(Iseaice).gt.0).and.                        &
     &              (Nval.le.Nimport(Iseaice))) THEN
                  MODELS(Iseaice)%ImportField(Nval)%short_name=         &
     &                                       TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('Export(seaice)')
              IF (MODELS(Iseaice)%IsActive) THEN
                IF ((Nexport(Iseaice).gt.0).and.                        &
     &              (Nval.le.Nexport(Iseaice))) THEN
                  MODELS(Iseaice)%ExportField(Nval)%short_name=         &
     &                                       TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('Nimport(waves)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (Ivalue(1).gt.0) THEN
                IF (MODELS(Iwave)%IsActive) THEN
                  IF (.not.allocated(MODELS(Iwave)%ImportField)) THEN
                    allocate ( MODELS(Iwave)%ImportField(Ivalue(1)) )
                  END IF
                  Nimport(Iwave)=Ivalue(1)
                END IF
              END IF
            CASE ('Nexport(waves)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              IF (Ivalue(1).gt.0) THEN
                IF (MODELS(Iwave)%IsActive) THEN
                  IF (.not.allocated(MODELS(Iwave)%ExportField)) THEN
                    allocate ( MODELS(Iwave)%ExportField(Ivalue(1)) )
                  END IF
                  Nexport(Iwave)=Ivalue(1)
                END IF
              END IF
            CASE ('Import(waves)')
              IF (MODELS(Iwave)%IsActive) THEN
                IF ((Nimport(Iwave).gt.0).and.                          &
     &              (Nval.le.Nimport(Iwave))) THEN
                  MODELS(Iwave)%ImportField(Nval)%short_name=           &
     &                                     TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('Export(waves)')
              IF (MODELS(Iwave)%IsActive) THEN
                IF ((Nexport(Iwave).gt.0).and.                          &
     &              (Nval.le.Nexport(Iwave))) THEN
                  MODELS(Iwave)%ExportField(Nval)%short_name=           &
     &                                     TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF

# ifdef DATA_COUPLING
            CASE ('nDataExport(roms)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Icomp=Iroms
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Icomp)%IsActive) THEN
                DataSet(Icomp)%Nfields=Ivalue(1)
                IF (DataSet(Icomp)%Nfields.gt.0) THEN
                  Nexport(Idata)=Nexport(Idata)+Ivalue(1)
                  IF (.not.allocated(DataSet(Icomp)%Export)) THEN
                    allocate ( DataSet(Icomp)%Export(Ivalue(1)) )
                  END IF
                  IF (.not.allocated(DataSet(Icomp)%Field)) THEN
                    allocate ( DataSet(Icomp)%Field(Ivalue(1)) )
                  END IF
                  IF (.not.allocated(DataSet(Icomp)%Ctarget)) THEN
                    allocate ( DataSet(Icomp)%Ctarget(Ivalue(1)) )
                  END IF
                  DO i=1,DataSet(Icomp)%Nfields
                    DataSet(Icomp)%Ctarget(i)=TRIM(MODELS(Icomp)%name)
                  END DO
                END IF
              END IF
            CASE ('nDataFiles(roms)')
              Npts=load_i(Nval, Rval, 1, nD2R)
              Icomp=Iroms
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Icomp)%IsActive) THEN
                DataSet(Icomp)%Nfiles=nD2R(1)
                IF (DataSet(Icomp)%Nfiles.gt.0) THEN
                  IF (.not.allocated(DataSet(Icomp)%IFS)) THEN
                    allocate( DataSet(Icomp)%IFS(nD2R(1)) )
                  END IF
                END IF
              END IF
              IF (DataSet(Icomp)%Nfiles.gt.0) THEN
                IF (allocated(Ncount)) THEN
                  deallocate ( Ncount )
                END IF
                allocate ( Ncount(nD2R(1),1) )
                Ncount=0
                IF (.not.allocated(D2R)) THEN
                  allocate ( D2R(nD2R(1),1) )    ! 2D structure so we
                END IF                           ! use "load_s2d"
              END IF
            CASE ('DataExport(roms)')
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Iroms)%IsActive) THEN
                Nfields=DataSet(Icomp)%Nfields
                IF ((Nfields.gt.0).and.(Nval.le.Nfields)) THEN
                  DataSet(Iroms)%Field(Nval)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('DataFiles(roms)')
              IF (MODELS(Idata)%IsActive.and.                           &
     &            (DataSet(Iroms)%Nfiles.gt.0)) THEN
                label='D2R - Data Model export fields to ROMS'
                Npts=load_s2d(Nval, Cval, Cdim, line, label, ifile,     &
     &                        igrid, 1, nD2R, Ncount, nD2R(1), inp_lib, &
     &                        D2R)
              END IF
            CASE ('nDataExport(atmos)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Icomp=Iatmos
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Icomp)%IsActive) THEN
                DataSet(Icomp)%Nfields=Ivalue(1)
                IF (DataSet(Icomp)%Nfields.gt.0) THEN
                  Nexport(Idata)=Nexport(Idata)+Ivalue(1)
                  IF (.not.allocated(DataSet(Icomp)%Field)) THEN
                    allocate ( DataSet(Icomp)%Field(Ivalue(1)) )
                  END IF
                  IF (.not.allocated(DataSet(Icomp)%Export)) THEN
                    allocate ( DataSet(Icomp)%Export(Ivalue(1)) )
                  END IF
                  IF (.not.allocated(DataSet(Icomp)%Ctarget)) THEN
                    allocate ( DataSet(Icomp)%Ctarget(Ivalue(1)) )
                  END IF
                  DO i=1,DataSet(Icomp)%Nfields
                    DataSet(Icomp)%Ctarget(i)=TRIM(MODELS(Icomp)%name)
                  END DO
                END IF
              END IF
            CASE ('nDataFiles(atmos)')
              Npts=load_i(Nval, Rval, 1, nD2A)
              Icomp=Iatmos
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Icomp)%IsActive) THEN
                DataSet(Icomp)%Nfiles=nD2A(1)
                IF (DataSet(Icomp)%Nfiles.gt.0) THEN
                  IF (.not.allocated(DataSet(Icomp)%IFS)) THEN
                    allocate ( DataSet(Icomp)%IFS(nD2A(1)) )
                  END IF
                END IF
              END IF
              IF (DataSet(Icomp)%Nfiles.gt.0) THEN
                IF (allocated(Ncount)) THEN
                  deallocate ( Ncount )
                END IF
                allocate ( Ncount(nD2A(1),1) )
                Ncount=0
                IF (.not.allocated(D2A)) THEN
                  allocate ( D2A(nD2A(1),1) )    ! 2D structure so we
                END IF                           ! use "load_s2d"
              END IF
            CASE ('DataExport(atmos)')
              IF (MODELS(Idata )%IsActive.and.                          &
     &            MODELS(Iatmos)%IsActive) THEN
                Nfields=DataSet(Iatmos)%Nfields
                IF ((Nfields.gt.0).and.(Nval.le.Nfields)) THEN
                  DataSet(Iatmos)%Field(Nval)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('DataFiles(atmos)')
              IF (MODELS(Idata)%IsActive.and.                           &
     &            (DataSet(Iatmos)%Nfiles.gt.0)) THEN
                label='D2A - Data Model export fields to ATM model'
                Npts=load_s2d(Nval, Cval, Cdim, line, label, ifile,     &
     &                        igrid, 1, nD2A, Ncount, nD2A(1), inp_lib, &
     &                        D2A)
              END IF
            CASE ('nDataExport(seaice)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Icomp=Iseaice
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Icomp)%IsActive) THEN
                DataSet(Icomp)%Nfields=Ivalue(1)
                IF (DataSet(Icomp)%Nfields.gt.0) THEN
                  Nexport(Idata)=Nexport(Idata)+Ivalue(1)
                  IF (.not.allocated(DataSet(Icomp)%Field)) THEN
                    allocate ( DataSet(Icomp)%Field(Ivalue(1)) )
                  END IF
                  IF (.not.allocated(DataSet(Icomp)%Export)) THEN
                    allocate ( DataSet(Icomp)%Export(Ivalue(1)) )
                  END IF
                  IF (.not.allocated(DataSet(Icomp)%Ctarget)) THEN
                    allocate ( DataSet(Icomp)%Ctarget(Ivalue(1)) )
                  END IF
                  DO i=1,DataSet(Icomp)%Nfields
                    DataSet(Icomp)%Ctarget(i)=TRIM(MODELS(Icomp)%name)
                  END DO
                END IF
              END IF
            CASE ('nDataFiles(seaice)')
              Npts=load_i(Nval, Rval, 1, nD2I)
              Icomp=Iseaice
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Icomp)%IsActive) THEN
                DataSet(Icomp)%Nfiles=nD2I(1)
                IF (DataSet(Icomp)%Nfiles.gt.0) THEN
                  IF (.not.allocated(DataSet(Icomp)%IFS)) THEN
                    allocate( DataSet(Icomp)%IFS(nD2I(1)) )
                  END IF
                END IF
              END IF
              IF (DataSet(Icomp)%Nfiles.gt.0) THEN
                IF (allocated(Ncount)) THEN
                  deallocate ( Ncount )
                END IF
                allocate ( Ncount(nD2I(1),1) )
                Ncount=0
                IF (.not.allocated(D2I)) THEN
                  allocate ( D2I(nD2I(1),1) )    ! 2D structure so we
                END IF                           ! use "load_s2d"
              END IF
            CASE ('DataExport(seaice)')
              IF (MODELS(Idata  )%IsActive.and.                         &
     &            MODELS(Iseaice)%IsActive) THEN
                Nfields=DataSet(Iseaice)%Nfields
                IF ((Nfields.gt.0).and.(Nval.le.Nfields)) THEN
                  DataSet(IseaIce)%Field(Nval)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('DataFiles(seaice)')
              IF (MODELS(Idata)%IsActive.and.                           &
     &            (DataSet(Iseaice)%Nfiles.gt.0)) THEN
                label='D2I - Data Model export fields to Sea ICE model'
                Npts=load_s2d(Nval, Cval, Cdim, line, label, ifile,     &
     &                        igrid, 1, nD2I, Ncount, nD2I(1), inp_lib, &
     &                        D2I)
              END IF
            CASE ('nDataExport(waves)')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Icomp=Iwave
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Icomp)%IsActive) THEN
                DataSet(Icomp)%Nfields=Ivalue(1)
                IF (DataSet(Icomp)%Nfields.gt.0) THEN
                  Nexport(Idata)=Nexport(Idata)+Ivalue(1)
                  IF (.not.allocated(DataSet(Icomp)%Field)) THEN
                    allocate ( DataSet(Icomp)%Field(Ivalue(1)) )
                  END IF
                  IF (.not.allocated(DataSet(Icomp)%Export)) THEN
                    allocate ( DataSet(Icomp)%Export(Ivalue(1)) )
                  END IF
                  IF (.not.allocated(DataSet(Icomp)%Ctarget)) THEN
                    allocate ( DataSet(Icomp)%Ctarget(Ivalue(1)) )
                  END IF
                  DO i=1,DataSet(Icomp)%Nfields
                    DataSet(Icomp)%Ctarget(i)=TRIM(MODELS(Icomp)%name)
                  END DO
                END IF
              END IF
            CASE ('nDataFiles(waves)')
              Npts=load_i(Nval, Rval, 1, nD2W)
              Icomp=Iwave
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Icomp)%IsActive) THEN
                DataSet(Icomp)%Nfiles=nD2W(1)
                IF (DataSet(Icomp)%Nfiles.gt.0) THEN
                  IF (.not.allocated(DataSet(Icomp)%IFS)) THEN
                    allocate( DataSet(Icomp)%IFS(nD2W(1)) )
                  END IF
                END IF
              END IF
              IF (DataSet(Icomp)%Nfiles.gt.0) THEN
                IF (allocated(Ncount)) THEN
                  deallocate ( Ncount )
                END IF
                allocate ( Ncount(nD2W(1),1) )
                Ncount=0
                IF (.not.allocated(D2W)) THEN
                  allocate ( D2W(nD2W(1),1) )    ! 2D structure so we
                END IF                           ! use "load_s2d"
              END IF
            CASE ('DataExport(waves)')
              IF (MODELS(Idata)%IsActive.and.                           &
     &            MODELS(Iwave)%IsActive) THEN
                Nfields=DataSet(Iwave)%Nfields
                IF ((Nfields.gt.0).and.(Nval.le.Nfields)) THEN
                  DataSet(Iwave)%Field(Nval)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END IF
            CASE ('DataFiles(waves)')
              IF (MODELS(Idata)%IsActive.and.                           &
     &            (DataSet(Iwave)%Nfiles.gt.0)) THEN
                label='D2W - Data Model export fields to Wave model'
                Npts=load_s2d(Nval, Cval, Cdim, line, label, ifile,     &
     &                        igrid, 1, nD2W, Ncount, nD2W(1), inp_lib, &
     &                        D2W)
              END IF
# endif
          END SELECT
        END IF
      END DO
  20  IF (localPET.eq.0) WRITE (out,40) line
      exit_flag=4
      rc=ESMF_RC_FILE_READ
      RETURN
  30  CLOSE (inp)
  40  FORMAT (/,' read_ESMconfig - Error while processing line: ',/,a)
!
!-----------------------------------------------------------------------
!  Set ESM, import state, and export state unique labels. They are
!  used for adding and advertising import/export state including
!  nested grids.
!-----------------------------------------------------------------------
!
!  Determine the number of DATA component states when connected to other
!  ESM components nested grids.  Recall that the DATA component has
!  a single grid with no nesting.
!
      IF (MODELS(Idata)%IsActive) THEN
        Nstates=1
        DO i=1,Nmodels
          ic=0
          IF (i.ne.Idata) THEN
            DO ng=1,MODELS(i)%Ngrids
              IF (COUPLED(i)%LinkedGrid(ng,Idata)) THEN
                ic=ic+1
              END IF
            END DO
            Nstates=MAX(Nstates,ic)
          END IF
        END DO
        IF (.not.allocated(COUPLED(Idata)%LinkedGrid)) THEN
          allocate ( COUPLED(Idata)%LinkedGrid(Nstates,Nmodels) )
          COUPLED(Idata)%LinkedGrid=.FALSE.
          DO i=1,Nmodels
            ic=0
            IF (i.ne.Idata) THEN
              DO ng=1,MODELS(i)%Ngrids
                IF (COUPLED(i)%LinkedGrid(ng,Idata)) THEN
                  ic=ic+1
                  COUPLED(Idata)%LinkedGrid(ic,i)=.TRUE.
                END IF
              END DO
            END IF
          END DO
        END IF
        IF (.not.allocated(COUPLED(Idata)%DataCoupledSets)) THEN
          allocate ( COUPLED(Idata)%DataCoupledSets(Nstates,Nmodels) )
          COUPLED(Idata)%DataCoupledSets=0
        END IF
        IF (.not.allocated(COUPLED(Idata)%ExportState)) THEN
          allocate ( COUPLED(Idata)%ExportState(Nstates,Nmodels) )
        END IF
      END IF
!
!  Allocate coupled sets and states variables, which are needed when
!  calling "NUOPC_AddNestedState"
!
      DO i=1,Nmodels
        IF (MODELS(i)%IsActive) THEN
          IF (i.eq.Idata) THEN
             Ngrd=Nstates
          ELSE
             Ngrd=MODELS(i)%Ngrids
          END IF
          IF (.not.allocated(MODELS(i)%grid)) THEN
            allocate ( MODELS(i)%grid(Ngrd) )
          END IF
          IF (.not.allocated(COUPLED(i)%SetLabel)) THEN
            allocate ( COUPLED(i)%SetLabel(Ngrd) )
            DO ng=1,Ngrd
              COUPLED(i)%SetLabel(ng)=blank
            END DO
          END IF
          IF (.not.allocated(COUPLED(i)%ExpLabel)) THEN
            allocate ( COUPLED(i)%ExpLabel(Ngrd) )
            DO ng=1,Ngrd
              COUPLED(i)%ExpLabel(ng)=blank
            END DO
          END IF
          IF (.not.allocated(COUPLED(i)%ImpLabel)) THEN
            allocate ( COUPLED(i)%ImpLabel(Ngrd) )
            DO ng=1,Ngrd
              COUPLED(i)%ImpLabel(ng)=blank
            END DO
          END IF
        END IF
      END DO
!
!  Determine the number of coupled state sets as the maximum number of
!  connected nested grids in a ESM component and track its location.
!  If Ncplsets=1, choose as location the component with the maximum
!  connections to other ESM component.
!
      Ncplsets=MAXVAL(MODELS(:)%Ngrids, mask=MODELS(:)%IsActive)
      IF ((Ncplsets.eq.1).and.(Nstates.eq.1)) THEN
        ic=0
        location(1)=1
        DO i=1,Nmodels
          IF (MODELS(i)%IsActive) THEN
            j=COUNT(COUPLED(i)%LinkedGrid)
            IF (j.gt.ic) THEN
              ic=j
              location(1)=i
            END IF
          END IF
        END DO
      ELSE
        location=MAXLOC(MODELS(:)%Ngrids, mask=MODELS(:)%IsActive)
      END IF
      ESMorder=0
!
!  The ESM component with the maximum number of connected nested grids
!  determines the different coupled sets. Order the active and connected
!  components accordingly to facilitate the setting of import and
!  export state labels.
!
      ic=1
      ESMorder(ic)=location(1)
      DO i=1,Nmodels
        IF (MODELS(i)%IsActive.and.(i.ne.location(1))) THEN
          ic=ic+1
          ESMorder(ic)=i
        END IF
      END DO
      ESMcount=ic
!
!  Define coupled sets labels.
!
      IF (.not.allocated(SetLabel)) THEN
        allocate ( SetLabel(Ncplsets) )
      END IF
      IF (.not.allocated(ExpLabel)) THEN
        allocate ( ExpLabel(Ncplsets) )
      END IF
      IF (.not.allocated(ImpLabel)) THEN
        allocate ( ImpLabel(Ncplsets) )
      END IF
      DO i=1,Ncplsets
        WRITE (SetLabel(i), '(a,i2.2)') 'ESM_', i
        WRITE (ExpLabel(i), '(a,i2.2)') 'Export_ESM_', i
        WRITE (ImpLabel(i), '(a,i2.2)') 'Import_ESM_', i
      END DO
!
!  Set coupled sets and import/export state labels.
!
      DO j=1,ESMcount
        Jcomp=ESMorder(j)
        IF (Jcomp.eq.Idata) THEN
          Jgrids=Nstates
        ELSE
          Jgrids=MODELS(Jcomp)%Ngrids
        END IF
        DO i=1,ESMcount
          Icomp=ESMorder(i)
          IF (Icomp.eq.Idata) THEN
            Igrids=Nstates
          ELSE
            Igrids=MODELS(Icomp)%Ngrids
          END IF
          IF (Icomp.ne.Jcomp) THEN
            DO jg=1,Jgrids
              DO ig=1,Igrids
                IF (COUPLED(Icomp)%LinkedGrid(ig,Jcomp).or.             &
     &              COUPLED(Jcomp)%LinkedGrid(jg,Icomp)) THEN
                  ng=MAX(ig,jg)
                  SetLstr=LEN_TRIM(COUPLED(Jcomp)%SetLabel(jg))
                  IF (SetLstr.eq.0) THEN
                    COUPLED(Jcomp)%SetLabel(jg)=TRIM(SetLabel(ng))
                  END IF
!
                  ExpLstr=LEN_TRIM(COUPLED(Jcomp)%ExpLabel(jg))
                  IF (Nexport(Jcomp).gt.0) THEN
                    IF (ExpLstr.eq.0) THEN
                      COUPLED(Jcomp)%ExpLabel(jg)=TRIM(ExpLabel(ng))
                    END IF
                  ELSE
                    IF (ExpLstr.eq.0) THEN
                      COUPLED(Jcomp)%ExpLabel(jg)='NONE'
                    END IF
                  END IF
!
                  ImpLstr=LEN_TRIM(COUPLED(Jcomp)%ImpLabel(jg))
                  IF (Nimport(Jcomp).gt.0) THEN
                    IF (ImpLstr.eq.0) THEN
                      COUPLED(Jcomp)%ImpLabel(jg)=TRIM(ImpLabel(ng))
                    END IF
                  ELSE
                    IF (ImpLstr.eq.0) THEN
                      COUPLED(Jcomp)%ImpLabel(jg)='NONE'
                    END IF
                  END IF
!
                  doit=COUPLED(Idata)%DataCoupledSets(jg,Jcomp).eq.0
                  IF (COUPLED(Jcomp)%LinkedGrid(jg,Idata).and.doit) THEN
                    COUPLED(Idata)%DataCoupledSets(jg,Jcomp)=ng
                  END IF
                END IF
              END DO
            END DO
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set ESM time managing variables.
!-----------------------------------------------------------------------
!
!  Create ESM calendar for driver and components.
!
      SELECT CASE (TRIM(lowercase(DateCalendar)))
        CASE ('gregorian')
          CalType=ESMF_CALKIND_GREGORIAN
        CASE ('year_360_day', '360_day')
          CalType=ESMF_CALKIND_360DAY
      END SELECT
!
      ClockInfo(Idriver)%Calendar=ESMF_CalendarCreate(CalType, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      ClockInfo(Idriver)%Name='Driver_clock'
!
!  Set ESM coupling simulation reference date number.
!
      CALL datenum (ReferenceDateNumber,                                &
     &              ReferenceDate(1),                                   &
     &              ReferenceDate(2),                                   &
     &              ReferenceDate(3),                                   &
     &              ReferenceDate(4),                                   &
     &              ReferenceDate(5),                                   &
     &              REAL(ReferenceDate(6),dp))
!
!  Set ESM coupling simulation reference time.
!
      CALL ESMF_TimeSet (ClockInfo(Idriver)%ReferenceTime,              &
     &                   yy=ReferenceDate(1),                           &
     &                   mm=ReferenceDate(2),                           &
     &                   dd=ReferenceDate(3),                           &
     &                   h= ReferenceDate(4),                           &
     &                   m= ReferenceDate(5),                           &
     &                   s= ReferenceDate(6),                           &
     &                   calkindflag=CalType,                           &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeGet (ClockInfo(Idriver)%ReferenceTime,              &
     &                   s_r8=ClockInfo(Idriver)%Time_Reference,        &
     &                   timeString=TimeReferenceString)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(TimeReferenceString, 'T')             ! remove 'T' in
      IF (is.gt.0) TimeReferenceString(is:is)=' '    ! ISO 8601 format
!
      ClockInfo(Idriver)%Time_ReferenceString=TRIM(TimeReferenceString)
      DO i=1,Nmodels
        ClockInfo(i)%ReferenceTime=ClockInfo(Idriver)%ReferenceTime
        ClockInfo(i)%Time_Reference=ClockInfo(Idriver)%Time_Reference
        ClockInfo(i)%Time_ReferenceString=TRIM(TimeReferenceString)
      END DO
!
!  Set ESM coupling driver interval.
!
      CALL ESMF_TimeIntervalSet (ClockInfo(Idriver)%TimeStep,           &
     &                           calendar=ClockInfo(Idriver)%Calendar,  &
     &                           yy=TimeStep(1),                        &
     &                           mm=TimeStep(2),                        &
     &                           d_r8=REAL(TimeStep(3),dp),             &
     &                           h= TimeStep(4),                        &
     &                           m= TimeStep(5),                        &
     &                           s= TimeStep(6),                        &
     &                           rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeIntervalGet (ClockInfo(Idriver)%TimeStep,           &
     &                           s_r8=ClockInfo(Idriver)%Time_Step)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      DO i=1,Nmodels
        ClockInfo(i)%TimeStep=ClockInfo(Idriver)%TimeStep
        ClockInfo(i)%Time_Step=ClockInfo(Idriver)%Time_Step
      END DO
!
!  Set ESM coupling driver Start Time.
!
      CALL ESMF_TimeSet (MyStartTime,                                   &
     &                   yy=StartDate(1),                               &
     &                   mm=StartDate(2),                               &
     &                   dd=StartDate(3),                               &
     &                   h= StartDate(4),                               &
     &                   m= StartDate(5),                               &
     &                   s= StartDate(6),                               &
     &                   calkindflag=CalType,                           &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF

# ifdef REGRESS_STARTCLOCK
!
!  If regressing start time clock, SUBSTRACT a coupling interval to
!  the provided start date to allow the proper initialization of the
!  import and export states. Also, it facilitates the processing and
!  exchange of the two time-level (LOWER and UPPER) snapshots between
!  ESM components before time stepping. The LOWER snapshop is exchanged
!  during initialization. The UPPER snapshot is exchanged on the first
!  call to the "XXXX_ModelAdvance" routine in the component NUOPC cap
!  module, but no time-stepping occurs.
!
      CALL ESMF_TimeGet (MyStartTime-ClockInfo(Idriver)%TimeStep,       &
     &                   yy=NewDate(1),                                 &
     &                   mm=NewDate(2),                                 &
     &                   dd=NewDate(3),                                 &
     &                   h= NewDate(4),                                 &
     &                   m= NewDate(5),                                 &
     &                   s= NewDate(6),                                 &
     &                   ms=NewDate(7),                                 &
     &                   timeString=TimeStartString,                    &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!                                                Regressed Start Time
      CALL ESMF_TimeSet (ClockInfo(Idriver)%StartTime,                  &
     &                   yy=NewDate(1),                                 &
     &                   mm=NewDate(2),                                 &
     &                   dd=NewDate(3),                                 &
     &                   h= NewDate(4),                                 &
     &                   m= NewDate(5),                                 &
     &                   s= NewDate(6),                                 &
     &                   ms=NewDate(7),                                 &
     &                   calkindflag=CalType,                           &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# else
      ClockInfo(Idriver)%StartTime=MyStartTime
# endif
!
      CALL ESMF_TimeGet (ClockInfo(Idriver)%StartTime,                  &
     &                   s_r8=ClockInfo(Idriver)%Time_Start,            &
     &                   timeString=TimeStartString)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(TimeStartString, 'T')                 ! remove 'T' in
      IF (is.gt.0) TimeStartString(is:is)=' '        ! ISO 8601 format
!
      ClockInfo(Idriver)%Time_StartString=TRIM(TimeStartString)
      DO i=1,Nmodels
        ClockInfo(i)%StartTime=ClockInfo(Idriver)%StartTime
        ClockInfo(i)%Time_Start=ClockInfo(Idriver)%Time_Start
        ClockInfo(i)%Time_StartString=TRIM(TimeStartString)
      END DO
!
!  Set ESM coupling driver Restart Time.
!
      CALL ESMF_TimeSet (MyRestartTime,                                 &
     &                   yy=RestartDate(1),                             &
     &                   mm=RestartDate(2),                             &
     &                   dd=RestartDate(3),                             &
     &                   h= RestartDate(4),                             &
     &                   m= RestartDate(5),                             &
     &                   s= RestartDate(6),                             &
     &                   calkindflag=CalType,                           &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF

# ifdef REGRESS_STARTCLOCK
!
!  If regressing start time clock, SUBSTRACT coupling interval to the
!  Restart Time for the same reasons mentioned above.
!
      CALL ESMF_TimeGet (MyRestartTime-ClockInfo(Idriver)%TimeStep,     &
     &                   yy=NewDate(1),                                 &
     &                   mm=NewDate(2),                                 &
     &                   dd=NewDate(3),                                 &
     &                   h= NewDate(4),                                 &
     &                   m= NewDate(5),                                 &
     &                   s= NewDate(6),                                 &
     &                   ms=NewDate(7),                                 &
     &                   timeString=TimeRestartString,                  &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!                                                Regressed Restart Time
      CALL ESMF_TimeSet (ClockInfo(Idriver)%RestartTime,                &
     &                   yy=NewDate(1),                                 &
     &                   mm=NewDate(2),                                 &
     &                   dd=NewDate(3),                                 &
     &                   h= NewDate(4),                                 &
     &                   m= NewDate(5),                                 &
     &                   s= NewDate(6),                                 &
     &                   ms=NewDate(7),                                 &
     &                   calkindflag=CalType,                           &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
# else
      ClockInfo(Idriver)%RestartTime=MyRestartTime
# endif
!
      CALL ESMF_TimeGet (ClockInfo(Idriver)%RestartTime,                &
     &                   s_r8=ClockInfo(Idriver)%Time_Restart,          &
     &                   timeString=TimeRestartString)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(TimeRestartString, 'T')               ! remove 'T' in
      IF (is.gt.0) TimeRestartString(is:is)=' '      ! ISO 8601 format
!
      ClockInfo(Idriver)%Time_RestartString=TRIM(TimeRestartString)
      DO i=1,Nmodels
        ClockInfo(i)%RestartTime=ClockInfo(Idriver)%RestartTime
        ClockInfo(i)%Time_Restart=ClockInfo(Idriver)%Time_Restart
        ClockInfo(i)%Time_RestartString=TRIM(TimeRestartString)
      END DO
!
!  ESM coupling driver stop time.
!
      CALL ESMF_TimeSet (ClockInfo(Idriver)%StopTime,                   &
     &                   yy=StopDate(1),                                &
     &                   mm=StopDate(2),                                &
     &                   dd=StopDate(3),                                &
     &                   h= StopDate(4),                                &
     &                   m= StopDate(5),                                &
     &                   s= StopDate(6),                                &
     &                   calkindflag=CalType,                           &
     &                   rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
!
      CALL ESMF_TimeGet (ClockInfo(Idriver)%StopTime,                   &
     &                   s_r8=ClockInfo(Idriver)%Time_Stop,             &
     &                   timeString=TimeStopString)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) THEN
        RETURN
      END IF
      is=INDEX(TimeStopString, 'T')                  ! remove 'T' in
      IF (is.gt.0) TimeStopString(is:is)=' '         ! ISO 8601 format
!
      ClockInfo(Idriver)%Time_StopString=TRIM(TimeStopString)
      DO i=1,Nmodels
        ClockInfo(i)%StopTime=ClockInfo(Idriver)%StopTime
        ClockInfo(i)%Time_Stop=ClockInfo(Idriver)%Time_Stop
        ClockInfo(i)%Time_StopString=TRIM(TimeStopString)
      END DO
!
!  Get time string from ROMS routine for debugging purposes.

      CALL time_string (ClockInfo(Idriver)%Time_Start-                  &
     &                  ClockInfo(Idriver)%Time_Reference, Fcode(1))
      CALL time_string (ClockInfo(Idriver)%Time_Stop-                   &
     &                  ClockInfo(Idriver)%Time_Reference, Fcode(2))

# ifdef DATA_COUPLING
!
!-----------------------------------------------------------------------
!  Load DATA model source files into Input File Structure (IFS).
!-----------------------------------------------------------------------
!
      IF (MODELS(Idata)%IsActive.and.                                   &
     &    MODELS(Iroms)%IsActive.and.                                   &
     &    (DataSet(Iroms)%Nfiles.gt.0)) THEN
        CALL load_IFS (D2R, nD2R(1), 1, nD2R, Iroms)
      END IF

      IF (MODELS(Idata )%IsActive.and.                                  &
     &    MODELS(Iatmos)%IsActive.and.                                  &
     &    (DataSet(Iatmos)%Nfiles.gt.0)) THEN
        CALL load_IFS (D2A, nD2A(1), 1, nD2A, Iatmos)
      END IF

      IF (MODELS(Idata  )%IsActive.and.                                 &
     &    MODELS(Iseaice)%IsActive.and.                                 &
     &    (DataSet(Iseaice)%Nfiles.gt.0)) THEN
        CALL load_IFS (D2I, nD2I(1), 1, nD2I, Iseaice)
      END IF

      IF (MODELS(Idata)%IsActive.and.                                   &
          MODELS(Iwave)%IsActive.and.                                   &
     &    (DataSet(Iwave)%Nfiles.gt.0)) THEN
        CALL load_IFS (D2W, nD2W(1), 1, nD2W, Iwave)
      END IF
!
!  Set total of DATA model export fields.
!
      Nimport(Idata)=0            ! DATA model does not to import fields
      Nexport(Idata)=0
      IF (MODELS(Idata)%IsActive) THEN
        DO i=1,Nmodels
          IF (MODELS(i)%IsActive.and.(i.ne.Idata)) THEN
            Nfields=DataSet(i)%Nfields
            Nexport(Idata)=Nexport(Idata)+Nfields
          END IF
        END DO
        IF (.not.allocated(MODELS(Idata)%ExportField)) THEN
          allocate ( MODELS(Idata)%ExportField(Nexport(Idata)) )
        END IF
        ic=0
        DO i=1,Nmodels
          IF (MODELS(i)%IsActive.and.(i.ne.Idata)) THEN
            Nfields=DataSet(i)%Nfields
            DO j=1,Nfields
              ic=ic+1
              MODELS(Idata)%ExportField(ic)%short_name=                 &
     &                                      TRIM(DataSet(i)%Field(j))
            END DO
          END IF
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Set several parameters.
!-----------------------------------------------------------------------
!
!  Assign PET list for each active ESM component.
!
      SELECT CASE (TRIM(PETlayoutOption))
        CASE ('SEQUENTIAL')
          DO j=1,Nmodels
            IF (MODELS(j)%IsActive) THEN
              DO i=1,MODELS(j)%nPETs
                MODELS(j)%PETlist(i)=i-1
              END DO
            END IF
          END DO
        CASE ('CONCURRENT')
          ic=-1
          sumPETs=0
          DO j=1,Nmodels
            IF (MODELS(j)%IsActive) THEN
              DO i=1,MODELS(j)%nPETs
                ic=ic+1
                MODELS(j)%PETlist(i)=ic
              END DO
              sumPETs=sumPETs+MODELS(j)%nPETs
            END IF
          END DO
      END SELECT
!
!  Assign PET list to connectors.  Notice that the DATA component only
!  exports fields and the connection is one way (DATA-TO-XXXX) and the
!  importing of fields is very unlikely (XXXX-TO-DATA). Thefore, the
!  import connector to the DATA component is never active.
!
      DO i=1,Nmodels
        DO j=1,Nmodels
          IF (CONNECTORS(i,j)%IsActive) THEN
            CONNECTORS(i,j)%name=TRIM(MODELS(i)%name)//'-TO-'//         &
     &                           TRIM(MODELS(j)%name)
            SELECT CASE (TRIM(PETlayoutOption))
              CASE ('SEQUENTIAL')
                nPETs=MODELS(i)%nPETs
                CONNECTORS(i,j)%nPETs=nPETs
                IF (.not.allocated(CONNECTORS(i,j)%PETlist)) THEN
                  allocate ( CONNECTORS(i,j)%PETlist(nPETs) )
                END IF
                DO k=1,nPETs
                  CONNECTORS(i,j)%PETlist(k)=MODELS(i)%PETlist(k)
                END DO
              CASE ('CONCURRENT')
                nPETs=MODELS(i)%nPETs+MODELS(j)%nPETs
                CONNECTORS(i,j)%nPETs=nPETs
                IF (.not.allocated(CONNECTORS(i,j)%PETlist)) THEN
                  allocate ( CONNECTORS(i,j)%PETlist(nPETs) )
                END IF
                DO k=1,MODELS(i)%nPETs
                  CONNECTORS(i,j)%PETlist(k)=MODELS(i)%PETlist(k)
                END DO
                ic=MODELS(i)%nPETs
                DO k=1,MODELS(i)%nPETs
                  ic=ic+1
                  CONNECTORS(i,j)%PETlist(ic)=MODELS(j)%PETlist(k)
                END DO
            END SELECT
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Report coupling input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        lstr=INDEX(my_fflags, 'free')-2
        IF (lstr.le.0) lstr=LEN_TRIM(my_fflags)
        WRITE (out,70) TRIM(ESMF_VERSION_STRING),                       &
     &                 TRIM(TodayDateString),                           &
     &                 TRIM(Rdir),                                      &
     &                 TRIM(my_os),                                     &
     &                 TRIM(my_cpu),                                    &
     &                 TRIM(my_fort),                                   &
     &                 TRIM(my_fc),                                     &
     &                 my_fflags(1:lstr),                               &
     &                 MyComm, PETcount
        WRITE (out,80) '     Coupling Input Parameters Filename = ',    &
     &                 TRIM(CinpName)
        WRITE (out,80) '         Coupling Run Sequence Filename = ',    &
     &                 TRIM(CONFname)
        WRITE (out,80) '       Coupling Input Metadata Filename = ',    &
     &                 TRIM(CPLname)
        WRITE (out,80) '         ROMS Input Parameters Filename = ',    &
     &                 TRIM(INPname(Iroms))

        IF (MODELS(Iatmos)%IsActive) THEN
          WRITE (out,80) '    ATM Model Input Parameters Filename = ',  &
     &                   TRIM(INPname(Iatmos))
        END IF
        IF (MODELS(Iseaice)%IsActive) THEN
          WRITE (out,80) ' SEAICE Model Input Parameters Filename = ',  &
     &                   TRIM(INPname(Iseaice))
        END IF
        IF (MODELS(Iwave)%IsActive) THEN
          WRITE (out,80) '   WAVE Model Input Parameters Filename = ',  &
     &                   TRIM(INPname(Iwave))
        END IF
!
        WRITE (out,'(a)') CHAR(10)                 ! new line

        IF (MODELS(Iroms)%IsActive) THEN
          IF (NgridsR.gt.1) THEN
            string=TRIM(MODELS(Iroms)%name)//                           &
     &             ' component with nested grids is coupled.'
          ELSE
            string=TRIM(MODELS(Iroms)%name)//' component is coupled.'
          END IF
        ELSE
          string=TRIM(MODELS(Iroms)%name)//' component is not coupled.'
        END IF
        WRITE (out,100) MODELS(Iroms)%IsActive, 'IsActive(roms)',       &
     &                  TRIM(string)

        IF (MODELS(Iatmos)%IsActive) THEN
          IF (NgridsA.gt.1) THEN
             string=TRIM(MODELS(Iatmos)%name)//                         &
     &              ' components with nested grids is coupled.'
          ELSE
             string=TRIM(MODELS(Iatmos)%name)//' component is coupled.'
          END IF
        ELSE
          string=TRIM(MODELS(Iatmos)%name)//' component is not coupled.'
        END IF
        WRITE (out,100) MODELS(Iatmos)%IsActive, 'IsActive(atmos)',     &
     &                  TRIM(string)

        IF (MODELS(Idata)%IsActive) THEN
           string=TRIM(MODELS(Idata)%name)//' component is coupled.'
        ELSE
           string=TRIM(MODELS(Idata)%name)//' component is not coupled.'
        END IF
        WRITE (out,100) MODELS(Idata)%IsActive, 'IsActive(data)',       &
     &                  TRIM(string)

        IF (MODELS(Iseaice)%IsActive) THEN
          IF (NgridsI.gt.1) THEN
            string=TRIM(MODELS(Iseaice)%name)//                         &
     &             ' components with nested grids is coupled.'
          ELSE
            string=TRIM(MODELS(Iseaice)%name)//' component is coupled.'
          END IF
        ELSE
         string=TRIM(MODELS(Iseaice)%name)//' component is not coupled.'
        END IF
        WRITE (out,100) MODELS(Iseaice)%IsActive, 'IsActive(seaice)',   &
     &                  TRIM(string)

        IF (MODELS(Iwave)%IsActive) THEN
          IF (NgridsW.gt.1) THEN
            string=TRIM(MODELS(Iwave)%name)//                           &
     &             ' components with nested grids is coupled.'
          ELSE
            string=TRIM(MODELS(Iwave)%name)//' component is coupled.'
          END IF
        ELSE
          string=TRIM(MODELS(Iwave)%name)//' component is not coupled.'
        END IF
        WRITE (out,100) MODELS(Iseaice)%IsActive, 'IsActive(waves)',    &
     &                  TRIM(string)
!
!  Notice that if the DATA component is activated, the reporting is
!  suppressed for COUPLED(Idata)%LinkedGrid(:,:) since the DATA
!  component does not import fields from others (like ATM -> DATA).
!  It only exports fields (like DATA -> ATM). That is, its connector
!  is only possible in the export direction. Its "LinkedGrid" switches
!  are only used for setting the coupled set and export state labels
!  above.
!
        DO j=1,Nmodels
          IF (MODELS(j)%IsActive.and.(j.ne.Idata)) THEN
            Ngrd=MODELS(j)%Ngrids
            DO i=1,Nmodels
              DO ng=1,Ngrd
                IF (COUPLED(j)%LinkedGrid(ng,i).and.(i.ne.j)) THEN
                  WRITE (out,110) COUPLED(j)%LinkedGrid(ng,i),          &
     &                            'Coupled('//                          &
     &                            TRIM(Clabel(i))//'2'//                &
     &                            TRIM(Clabel(j))//')',                 &
     &                            TRIM(MODELS(i)%name),                 &
     &                            TRIM(MODELS(j)%name), ng
                END IF
              END DO
            END DO
          END IF
        END DO
!
        IF (MODELS(Idata  )%IsActive.or.                                &
     &      MODELS(Iatmos )%IsActive.or.                                &
     &      MODELS(Iseaice)%IsActive.or.                                &
     &      MODELS(Iwave  )%IsActive) THEN
          IF (CouplingType.eq.1) THEN
            string='Explicit coupling method.'
          ELSE IF (CouplingType.eq.2) THEN
            string='Semi-Implicit coupling method.'
          ELSE
            WRITE (out,210) 'CouplingType', CouplingType
            rc=ESMF_RC_VAL_WRONG
            exit_flag=7
            RETURN
          END IF
          WRITE (out,120) CouplingType, 'CouplingType', TRIM(string)
        END IF
!
        IF (MODELS(Iatmos )%IsActive.or.                                &
     &      MODELS(Idata  )%IsActive.or.                                &
     &      MODELS(Iseaice)%IsActive.or.                                &
     &      MODELS(Iwave  )%IsActive) THEN
          SELECT CASE(TRIM(PETlayoutOption))
            CASE ('SEQUENTIAL')
              string='Sequential, models run on all PETs.'
            CASE ('CONCURRENT')
              string='Concurrent, each model runs on a subset of PETs.'
            CASE DEFAULT
              WRITE (out,210) 'PETlayoutOption', PETlayoutOption
              rc=ESMF_RC_VAL_WRONG
              exit_flag=7
              RETURN
          END SELECT
          WRITE (out,130) TRIM(PETlayoutOption), 'PETlayoutOption',     &
     &                    TRIM(string)
        END IF
!
        IF (MODELS(Idata )%IsActive) THEN
          WRITE (out,120) ItileD, 'ItileD',                             &
     &          'DATA model tile partition in the I-direction.'
          WRITE (out,120) JtileD, 'JtileD',                             &
     &          'DATA model tile partition in the J-direction.'
        END IF
!
        IF (MODELS(Iatmos )%IsActive.or.                                &
     &      MODELS(Idata  )%IsActive.or.                                &
     &      MODELS(Iseaice)%IsActive.or.                                &
     &      MODELS(Iwave  )%IsActive) THEN
          SELECT CASE(TRIM(PETlayoutOption))
            CASE ('SEQUENTIAL')
              IF (MODELS(Iroms)%IsActive) THEN
                Icomp=Iroms
                nPETs=MODELS(Icomp)%nPETs
                IF (nPETs.eq.PETcount) THEN
                  WRITE (out,150) nPETs, 'Nthreads(roms)',              &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//'.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                ELSE
                  WRITE (out,220) 'Nthreads(roms)', ng, nPETs,          &
     &                            PETcount
                  rc=ESMF_RC_VAL_WRONG
                  exit_flag=7
                  RETURN
                END IF
              END IF
              IF (MODELS(Iatmos)%IsActive) THEN
                Icomp=Iatmos
                nPETs=MODELS(Icomp)%nPETs
                IF (nPETs.eq.PETcount) THEN
                  WRITE (out,150) nPETs, 'Nthreads(atmos)',             &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//'.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                ELSE
                  WRITE (out,230) 'Nthreads(atmos)', nPETs, PETcount
                  rc=ESMF_RC_VAL_WRONG
                  exit_flag=7
                  RETURN
                END IF
              END IF
              IF (MODELS(Iseaice)%IsActive) THEN
                Icomp=Iseaice
                nPETs=MODELS(Icomp)%nPETs
                IF (nPETs.eq.PETcount) THEN
                  WRITE (out,150) nPETs, 'Nthreads(seaice)',            &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//'.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                ELSE
                  WRITE (out,230) 'Nthreads(seaice)', nPETs, PETcount
                  rc=ESMF_RC_VAL_WRONG
                  exit_flag=7
                  RETURN
                END IF
              END IF
              IF (MODELS(Iwave)%IsActive) THEN
                Icomp=Iwave
                nPETs=MODELS(Iwave)%nPETs
                IF (nPETs.eq.PETcount) THEN
                  WRITE (out,150) nPETs, 'Nthreads(waves)',             &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//'.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                ELSE
                  WRITE (out,230) 'Nthreads(waves)', nPETs, PETcount
                  rc=ESMF_RC_VAL_WRONG
                  exit_flag=7
                  RETURN
                END IF
              END IF
              IF (MODELS(Idata)%IsActive) THEN
                Icomp=Idata
                nPETs=MODELS(Icomp)%nPETs
                IF (nPETs.eq.PETcount) THEN
                  WRITE (out,150) nPETs, 'Nthreads(data)',              &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//            &
     &                            ', ItileD * JtileD.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                ELSE
                  WRITE (out,230) 'ItileD * JtileD', nPETs, PETcount
                  rc=ESMF_RC_VAL_WRONG
                  exit_flag=7
                  RETURN
                END IF
              END IF
            CASE ('CONCURRENT')
              IF (sumPETs.eq.PETcount) THEN
                IF (MODELS(Iroms)%IsActive) THEN
                  Icomp=Iroms
                  nPETs=MODELS(Icomp)%nPETs
                  WRITE (out,150) nPETs, 'Nthreads(roms)',              &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//'.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                END IF
                IF (MODELS(Idata)%IsActive) THEN
                  Icomp=Idata
                  nPETs=MODELS(Icomp)%nPETs
                  WRITE (out,150) nPETs, 'Nthreads(data)',              &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//            &
     &                            ', ItileD * JtileD.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                END IF
                IF (MODELS(Iatmos)%IsActive) THEN
                  Icomp=Iatmos
                  nPETs=MODELS(Icomp)%nPETs
                  WRITE (out,150) nPETs, 'Nthreads(atmos)',             &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//'.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                END IF
                IF (MODELS(Iseaice)%IsActive) THEN
                  Icomp=Iseaice
                  nPETs=MODELS(Icomp)%nPETs
                  WRITE (out,150) nPETs, 'Nthreads(seaice)',            &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//'.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                END IF
                IF (MODELS(Iwave)%IsActive) THEN
                  Icomp=Iwave
                  nPETs=MODELS(Icomp)%nPETs
                  WRITE (out,150) nPETs, 'Nthreads(waves)',             &
     &                            'Assigned number of PETs for '//      &
     &                            TRIM(MODELS(Icomp)%name)//'.'
                  WRITE (Pstr,'(i5)') MODELS(Icomp)%PETlist(1)
                  WRITE (Pend,'(i5)') MODELS(Icomp)%PETlist(nPETS)
                  IF (nPETS.gt.1) THEN
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))//' to '//       &
     &                              TRIM(ADJUSTL(Pend))
                  ELSE
                    WRITE (out,290) TRIM(ADJUSTL(Pstr))
                  END IF
                END IF
              ELSE
                WRITE (out,240) sumPETs, PETcount
                rc=ESMF_RC_VAL_WRONG
                exit_flag=7
                RETURN
              END IF
          END SELECT
        ELSE                                          ! ROMS standalone
          IF (MODELS(Iroms)%IsActive) THEN
            nPETs=MODELS(Iroms)%nPETs
            IF (nPETs.eq.PETcount) THEN
              WRITE (out,150) nPETs, 'Nthreads(roms)',                  &
     &                        'Assigned number of PETs for ROMS.'
            ELSE
              WRITE (out,220) 'Nthreads(roms)', ng, nPETs, PETcount
              rc=ESMF_RC_VAL_WRONG
              exit_flag=7
              RETURN
            END IF
          END IF
        END IF
!
        SELECT CASE (TRIM(DateCalendar))
          CASE ('gregorian')
            WRITE (out,170) '   gregorian', 'Calendar',                 &
     &                      'ESM components date calendar.'
          CASE ('year_360_day')
            WRITE (out,170) 'year_360_day', 'Calendar',                 &
     &                      'ESM components date calendar.'
          CASE DEFAULT
            WRITE (out,250) 'Calendar', TRIM(DateCalendar)
            rc=ESMF_RC_VAL_WRONG
            exit_flag=7
            RETURN
        END SELECT
!
        WRITE (out,160) (ReferenceDate(i), i=1,6), 'ReferenceTime',     &
     &                  'Coupling driver reference time.'
        WRITE (out,160) (StartDate(i), i=1,6), 'StartTime',             &
     &                  'Simulation start time.'
        WRITE (out,160) (RestartDate(i), i=1,6), 'RestartTime',         &
     &                  'Simulation re-start time.'
        WRITE (out,160) (StopDate(i), i=1,6), 'StopTime',               &
     &                  'Simulation stop time.'
        IF (MODELS(Iatmos )%IsActive.or.                                &
     &      MODELS(Idata  )%IsActive.or.                                &
     &      MODELS(Iseaice)%IsActive.or.                                &
     &      MODELS(Iwave  )%IsActive) THEN
          WRITE (out,160) (TimeStep(i), i=1,6), 'TimeStep',             &
     &                    'Driver coupling time interval.'
        ELSE
          WRITE (out,160) (TimeStep(i), i=1,6), 'TimeStep',             &
     &                    'Simulation elapsed time interval.'
        END IF
!
        DO j=1,Nmodels
          IF (MODELS(j)%IsActive) THEN
            IF (j.eq.Idata) THEN
              Ngrd=1                          ! DATA has no nested grids
            ELSE
              Ngrd=SIZE(COUPLED(j)%LinkedGrid,DIM=1)
            END IF
            DO i=1,Nmodels
              DO ng=1,Ngrd
                IF (COUPLED(j)%LinkedGrid(ng,i).and.(i.ne.j)) THEN
                  TimeFrac=MODELS(j)%TimeFrac(ng,i)
                  string='TimeFrac('//                                  &
     &                   TRIM(Clabel(i))//'2'//                         &
     &                   TRIM(Clabel(j))//')'
                  IF (TimeFrac.gt.0) THEN
                    WRITE (out,180) TimeFrac, TRIM(string),             &
     &                              'Coupling TimeStep fraction for '// &
     &                              TRIM(MODELS(i)%name)//' -> '//      &
     &                              TRIM(MODELS(j)%name)//', grid: ', ng
                  ELSE
                    WRITE (out,210) TRIM(string), TimeFrac
                    rc=ESMF_RC_VAL_WRONG
                    exit_flag=7
                    RETURN
                  END IF
                END IF
              END DO
            END DO
          END IF
        END DO
!
        WRITE (out,120) extrapNumLevels, 'extrapNumLevels',             &
     &        'Number of levels for creep fill extrapolation.'

# ifdef DATA_COUPLING
!
        IF (MODELS(Iatmos)%IsActive.and.                                &
     &      MODELS(Idata )%IsActive) THEN
          WRITE (out,310) 1, 'WeightsFile(atmos)',                      &
     &                    TRIM(MODELS(Iatmos)%name)//                   &
     &          ' component melding weights coefficients NetCDF file:', &
     &                    TRIM(WEIGHTS(Iatmos)%ncfile)
          WRITE (out,320) TRIM(WEIGHTS(Iatmos)%VnameDATA),              &
     &                    'VnameDATA(atmos)',                           &
     &          'DATA component weights NetCDF variable name.'
          WRITE (out,320) TRIM(WEIGHTS(Iatmos)%VnameESM),               &
     &                    'VnameESM(atmos)',                            &
     &          'ESM  component weights NetCDF variable name.'
          WRITE (out,190) WEIGHTS(Iatmos)%NestedGrid,                   &
     &                    'NestedGrid(atmos)',                          &
     &          'Grid needing merged fields from DATA-ESM components.'
        END IF
# endif
!
        IF ((0.le.DebugLevel).and.(DebugLevel.le.4)) THEN
          WRITE (out,190) DebugLevel, 'DebugLevel',                     &
     &          'Coupling debugging level flag.'
        ELSE
          WRITE (out,210) 'DebugLevel', DebugLevel
          rc=ESMF_RC_VAL_WRONG
          exit_flag=7
          RETURN
        END IF
!
        WRITE (out,190) TraceLevel, 'traceLevel',                       &
     &        'Execution tracing level flag.'
!
        DO j=1,Nmodels
          IF (MODELS(j)%IsActive.and.(j.ne.Idata)) THEN
            string=lowercase(Cmodel(j))
            WRITE (out,190) Nimport(j),                                 &
     &                      'Nimport('//TRIM(string)//')',              &
     &                      'Number of '//TRIM(MODELS(j)%name)//        &
     &                      ' component import fields:'
            DO i=1,Nimport(j)
              WRITE (out,200) i,                                        &
     &                        TRIM(MODELS(j)%ImportField(i)%short_name)
            END DO
            WRITE (out,190) Nexport(j),                                 &
     &                      'Nexport('//TRIM(string)//')',              &
     &                      'Number of '//TRIM(MODELS(j)%name)//        &
     &                      ' component export fields:'
            DO i=1,Nexport(j)
              WRITE (out,200) i,                                        &
     &                        TRIM(MODELS(j)%ExportField(i)%short_name)
            END DO
!
            IF ((Nimport(j)+Nexport(j)).lt.1) THEN
              WRITE (out,260) 'ERROR: The '//TRIM(MODELS(j)%name)//     &
     &              ' component is not importing or exporting fields.', &
     &              'Revise '''//TRIM(CinpName)//''' script for '//     &
     &              'Nimport('//TRIM(string)//') or '//                 &
     &              'Nexport('//TRIM(string)//').'
              rc=ESMF_RC_VAL_WRONG
              exit_flag=7
              RETURN
            END IF
          END IF
        END DO
!
        DO j=1,Nmodels
          string=lowercase(Cmodel(j))
          IF (j.ne.Idata) THEN
            IF (MODELS(Idata)%IsActive.and.                             &
     &          MODELS(j)%IsActive.and.(DataSet(j)%Nfields.gt.0)) THEN
              WRITE (out,270) DataSet(j)%Nfields,                       &
     &                        'nDataExport('//TRIM(string)//')',        &
     &                        'Number of export DATA model fields to '//&
     &                        TRIM(MODELS(j)%name)//' component:'
              DO i=1,DataSet(j)%Nfields
                WRITE (out,200) i, TRIM(DataSet(j)%Field(i))
              END DO
              WRITE (out,270) DataSet(j)%Nfiles,                        &
     &                        'nDataFiles('//TRIM(string)//')',         &
     &                        'Number of source DATA model files for '//&
     &                        TRIM(MODELS(j)%name)//' component:'
              KeyWord='DataFiles('//TRIM(string)//')'
              DO i=1,DataSet(j)%Nfiles
                DO ifile=1,DataSet(j)%IFS(i)%Nfiles
                  fname=DataSet(j)%IFS(i)%files(ifile)
                  IF (find_file(ng, out, fname, TRIM(KeyWord))) THEN
                    IF (ifile.eq.1) THEN
                      WRITE (out,280) i, ': ', TRIM(fname)
                    ELSE
                      WRITE (out,'(t42,6x,a)') TRIM(fname)
                    END IF
                  ELSE
                    rc=ESMF_RC_NOT_FOUND
                    exit_flag=4
                    WRITE (out,300) TRIM(MODELS(j)%name), ng,           &
     &                              TRIM(fname)
                  END IF
                END DO
              END DO
            END IF
          END IF
        END DO
!
        WRITE (out,340) 'Coupled Import and Export States Name Sets:'
        DO i=1,Nmodels
          IF (MODELS(i)%IsActive) THEN
            IF (i.eq.Idata) THEN
              Ngrd=Nstates
            ELSE
              Ngrd=MODELS(i)%Ngrids
            END IF
            DO ng=1,Ngrd
              is=1
              ie=1
              string=' '
              DO j=1,Nmodels
                IF (MODELS(j)%IsActive.and.(j.ne.i)) THEN
                  IF (i.eq.Idata) THEN
                    ig=1                       ! single grid, no nesting
                    IF (COUPLED(j)%LinkedGrid(ig,i)) THEN
                      ie=is+LEN_TRIM(MODELS(j)%name)
                      string(is:ie)=TRIM(MODELS(j)%name)
                      is=ie+2
                    END IF
                  ELSE
                    ig=ng
                    IF (COUPLED(i)%LinkedGrid(ng,j)) THEN
                      ie=is+LEN_TRIM(MODELS(j)%name)
                      string(is:ie)=TRIM(MODELS(j)%name)
                      is=ie+2
                    END IF
                  END IF
                END IF
              END DO
              WRITE (out,350) TRIM(MODELS(i)%name), ig,                 &
     &                        TRIM(COUPLED(i)%SetLabel(ng)),            &
     &                        TRIM(COUPLED(i)%ImpLabel(ng)),            &
     &                        TRIM(COUPLED(i)%ExpLabel(ng)),            &
     &                        TRIM(string)
            END DO
          END IF
        END DO
!
      END IF
!
!  Flush standard output buffer.
!
      CALL my_flush (out)
!
  70  FORMAT (80('-'),/,                                                &
     &        ' Earth System Models Coupling: ESMF/NUOPC Library,',     &
     &        ' Version ',a,/,31x,a,/,                                  &
     &        80('-'),                                                  &
     &        /,1x,'Repository Root  : ',a,                             &
     &        /,1x,'Operating System : ',a,                             &
     &        /,1x,'CPU Hardware     : ',a,                             &
     &        /,1x,'Compiler System  : ',a,                             &
     &        /,1x,'Compiler Command : ',a,                             &
     &        /,1x,'Compiler Flags   : ',a,                             &
     &        /,1x,'MPI Communicator : ',i0,2x,'PET size = ',i0,        &
     &        /,80('-'),/)
  80  FORMAT (1x,a,a)
  90  FORMAT (18x,l1,2x,a,t42,a,i2.2,'.')
 100  FORMAT (18x,l1,2x,a,t42,a)
 110  FORMAT (18x,l1,2x,a,t42,'Connector ',a,' -> ',a,                  &
     &        ' is activated, grid: ',i2.2,'.')
 120  FORMAT (9x,i10,2x,a,t42,a)
 130  FORMAT (9x,a,2x,a,t42,a)
 150  FORMAT (15x,i4,2x,a,t42,a)
 160  FORMAT (i4.4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2,2x,      &
     &        a,t42,a)
 170  FORMAT (7x,a,2x,a,t42,a)
 180  FORMAT (18x,i1,2x,a,t42,a,i2.2,'.')
 190  FORMAT (18x,i1,2x,a,t42,a,i2.2)
 200  FORMAT (t42,2x,i2.2,':',1x,a)
 210  FORMAT (/,' read_ESMconfig - Invalid input parameter, ',a,        &
     &          ': ',i1)
 220  FORMAT (/,' read_ESMconfig - Illegal input parameter, ',a,        &
     &          ' for grid ',i2.2,': ',2i4)
 230  FORMAT (/,' read_ESMconfig - Illegal input parameter, ',a,        &
     &          ': ',2i4)
 240  FORMAT (/,' read_ESMconfig - Inconsistent number of PETs for',    &
     &          ' concurrent execution:',/,18x,'sumPETs = ',i4,2x,      &
     &          ' PETcount = ',i4)
 250  FORMAT (/,' read_ESMconfig - Invalid input parameter, ',a,        &
     &          ': ',a)
 260  FORMAT (/,21x,a,/,21x,a)
 270  FORMAT (18x,i1,2x,a,t42,a,i2.2,'.')
 280  FORMAT (t42,2x,i2.2,a,a)
 290  FORMAT (t44,'Coupling Driver PETs: ',a)
 300  FORMAT (/,' read_ESMconfig - ',a,' Grid ',i2.2,                   &
     &        ', could not find input file:',/,18x,a)
 310  FORMAT (17x,i2,2x,a,t42,a,/,t44,a)
 320  FORMAT (a19,t22,a,t42,a)
 330  FORMAT (6x,1p,e13.6,2x,a,t42,a)
 340  FORMAT (/,a,/, 42('='),/,/, 'Component', t13,'Grid',              &
     &        t20,'CoupledSet', t33,'ImportState', t50,'ExportState',   &
     &        t67,'ConnectedTo',/, 114('-'))
 350  FORMAT (a,t14,i0,t20,a,t33,a,t50,a,t67,a)
!
      RETURN
      END SUBROUTINE read_ESMconfig
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
      SUBROUTINE set_metadata (vm, rc)
!
!=======================================================================
!                                                                      !
!  Process ESM coupling import and export fields metadata. It Adds     !
!  fields to the NUOPC dictionary.                                     !
!                                                                      !
!=======================================================================
!
      USE mod_scalars,      ONLY : NoError, exit_flag
      USE get_metadata_mod, ONLY : CouplingField,                       &
     &                             coupling_metadata,                   &
     &                             metadata_has
      USE strings_mod,      ONLY : FoundError, assign_string, uppercase
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_VM) :: vm
!
!  Local variable declarations.
!
      TYPE (CouplingField), allocatable :: S(:)
!
      logical :: Exist, Lreport, connected
!
      integer :: i, id, io_err, j, lvar, varid
      integer :: localPET, PETcount
      integer :: etype, gtype, itype
!
      integer, parameter :: inp = 10
!
      real(r8) :: add_offset, scale
!
      character (len=:), allocatable :: Smodel

      character (len=100) :: Sname, units, val
      character (len=100) :: ShortName, line
      character (len=256) :: io_errmsg
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", set_metadata"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Query gridded component.
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
!  Get ESM coupling export/import variable metadata.
!-----------------------------------------------------------------------
!
!  The metadata structure can be read either a YAML file (extention
!  .yaml) or deprecated ASCII file (extension .dat).
!
      CALL coupling_metadata (TRIM(CPLname), S)
!
!-----------------------------------------------------------------------
!  Search IMPORT fields in metadata dictionary and load information.
!-----------------------------------------------------------------------
!
      DO i=1,Nmodels
        IF (MODELS(i)%IsActive) THEN
          DO j=1,Nimport(i)
            ShortName=MODELS(i)%ImportField(j)%short_name
            id=metadata_has(S, TRIM(ShortName))
            IF (id.gt.0) THEN
              MODELS(i)%ImportField(j)%connected    = S(id)%connected
              MODELS(i)%ImportField(j)%debug_write  = S(id)%debug_write
              MODELS(i)%ImportField(j)%add_offset   = S(id)%add_offset
              MODELS(i)%ImportField(j)%scale_factor = S(id)%scale
!
!                                      field short name keyword
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%short_name,     &
     &                                     S(id)%short_name),           &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      field standard name
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%standard_name,  &
     &                                     S(id)%standard_name),        &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      field descriptive long name
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%long_name,      &
     &                                     S(id)%long_name),            &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      destination field grid-cell type
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%dst_gtype,      &
     &                                     S(id)%destination_grid),     &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      destination field units
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%dst_units,      &
     &                                     S(id)%destination_units),    &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      source field grid-cell type
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%src_gtype,      &
     &                                     S(id)%source_grid),          &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      source field units
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%src_units,      &
     &                                     S(id)%source_units),         &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      DATA NetCDF variable name
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%nc_vname,       &
     &                                     S(id)%data_netcdf_vname),    &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      DATA NetCDF time variable name
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%nc_tname,       &
     &                                     S(id)%data_netcdf_tname),    &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      field reggriding method
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%RegridMethod,   &
     &                                     S(id)%regrid_method),        &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      unmapped extrapolation method
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ImportField(j)%ExtrapMethod,   &
     &                                     S(id)%extrapolate_method),   &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      source model (used for reporting)
              IF (FoundError(assign_string(Smodel,                      &
     &                                     S(id)%connected_to),         &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!
              val=uppercase(TRIM(MODELS(i)%ImportField(j)%src_gtype))
              SELECT CASE (TRIM(val))
                CASE ('CENTER_CELL', 'CELL_CENTER', 'CENTER')
                  gtype=Icenter
                CASE ('CORNER_CELL', 'CELL_CORNER', 'CORNER')
                  gtype=Icorner
                CASE ('U','LEFT_RIGHT_EDGE', 'RIGHT_LEFT_EDGE')
                  gtype=Iupoint
                CASE ('V','LOWER_UPPER_EDGE', 'UPPER_LOWER_EDGE')
                  gtype=Ivpoint
                CASE DEFAULT
                  gtype=Icenter
              END SELECT
              MODELS(i)%ImportField(j)%gtype=gtype
!
              val=uppercase(TRIM(MODELS(i)%ImportField(j)%RegridMethod))
              SELECT CASE (TRIM(val))
                CASE ('BILINEAR')
                  itype=Ibilin
                CASE ('PATCH')
                  itype=Ipatch
                CASE ('CONSERVATIVE1')
                  itype=Iconsv1
                CASE ('CONSERVATIVE2')
                  itype=Iconsv2
                CASE ('NEAREST')
                  itype=InStoD
                CASE DEFAULT
                  itype=Ibilin
              END SELECT
              MODELS(i)%ImportField(j)%itype=itype
!
              val=uppercase(TRIM(MODELS(i)%ImportField(j)%ExtrapMethod))
              SELECT CASE (TRIM(val))
                CASE ('NONE')
                  etype=Enone
                CASE ('NEAREST')
                  etype=ExStoD
                CASE ('IDAVG')
                  etype=Eidavg
                CASE ('CREEP')
                  etype=Ecreep
                CASE ('2STEPS')
                  etype=E2steps
                CASE DEFAULT
                  etype=Enone
              END SELECT
              MODELS(i)%ImportField(j)%etype=etype
            ELSE
              IF (localPET.eq.0) THEN
                WRITE (cplout,30) 'import field short_name: ',          &
     &                             TRIM(ShortName), TRIM(CPLname)
              END IF
              rc=ESMF_RC_NOT_FOUND
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) RETURN
            END IF
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Search EXPORT fields in metadata dictionary and load information.
!-----------------------------------------------------------------------
!
!  Notice that all export fields are kept in their original form.
!  The imported component does the proper scaling, physical units
!  conversion, and other manipulations. It is done to avoid applying
!  such transformations twice.
!
      DO i=1,Nmodels
        IF (MODELS(i)%IsActive) THEN
          DO j=1,Nexport(i)
            ShortName=MODELS(i)%ExportField(j)%short_name
            id=metadata_has(S, TRIM(ShortName))
            IF (id.gt.0) THEN
              MODELS(i)%ExportField(j)%connected    = S(id)%connected
              MODELS(i)%ExportField(j)%debug_write  = S(id)%debug_write
              MODELS(i)%ExportField(j)%add_offset   = 0.0_r8
              MODELS(i)%ExportField(j)%scale_factor = 1.0_r8 !no scaling
!
!                                      field short name keyword
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%short_name,     &
     &                                     S(id)%short_name),           &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      field standard name
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%standard_name,  &
     &                                     S(id)%standard_name),        &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      field descriptive long name
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%long_name,      &
     &                                     S(id)%long_name),            &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      destination field grid-cell type
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%dst_gtype,      &
     &                                     S(id)%destination_grid),     &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      destination field units
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%dst_units,      &
     &                                     S(id)%destination_units),    &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      source field grid-cell type
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%src_gtype,      &
     &                                     S(id)%source_grid),          &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      source field units
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%src_units,      &
     &                                     S(id)%source_units),         &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      DATA NetCDF variable name
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%nc_vname,       &
     &                                     S(id)%data_netcdf_vname),    &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      DATA NetCDF time variable name
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%nc_tname,       &
     &                                     S(id)%data_netcdf_tname),    &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      field reggriding method
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%RegridMethod,   &
     &                                     S(id)%regrid_method),        &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      unmapped extrapolation method
              IF (FoundError(assign_string(                             &
     &                         MODELS(i)%ExportField(j)%ExtrapMethod,   &
     &                                     S(id)%extrapolate_method),   &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!                                      source model (used for reporting)
              IF (FoundError(assign_string(Smodel,                      &
     &                                     S(id)%connected_to),         &
     &                       NoError, __LINE__, MyFile)) THEN
                rc=ESMF_RC_COPY_FAIL
                RETURN
              END IF
!
              val=uppercase(TRIM(MODELS(i)%ExportField(j)%src_gtype))
              SELECT CASE (TRIM(val))
                CASE ('CENTER_CELL', 'CELL_CENTER', 'CENTER')
                  gtype=Icenter
                CASE ('CORNER_CELL', 'CELL_CORNER', 'CORNER')
                  gtype=Icorner
                CASE ('U','LEFT_RIGHT_EDGE', 'RIGHT_LEFT_EDGE')
                  gtype=Iupoint
                CASE ('V','LOWER_UPPER_EDGE', 'UPPER_LOWER_EDGE')
                  gtype=Ivpoint
                CASE DEFAULT
                  gtype=Icenter
              END SELECT
              MODELS(i)%ExportField(j)%gtype=gtype
!
              val=uppercase(TRIM(MODELS(i)%ExportField(j)%RegridMethod))
              SELECT CASE (TRIM(val))
                CASE ('BILINEAR')
                  itype=Ibilin
                CASE ('PATCH')
                  itype=Ipatch
                CASE ('CONSERVATIVE1')
                  itype=Iconsv1
                CASE ('CONSERVATIVE2')
                  itype=Iconsv2
                CASE ('NEAREST')
                  itype=InStoD
                CASE DEFAULT
                  itype=Ibilin
              END SELECT
              MODELS(i)%ExportField(j)%itype=itype
!
              val=uppercase(TRIM(MODELS(i)%ExportField(j)%ExtrapMethod))
              SELECT CASE (TRIM(val))
                CASE ('NONE')
                  etype=Enone
                CASE ('NEAREST')
                  etype=ExStoD
                CASE ('IDAVG')
                  etype=Eidavg
                CASE ('CREEP')
                  etype=Ecreep
                CASE ('2STEPS')
                  etype=E2steps
                CASE DEFAULT
                  etype=Enone
              END SELECT
              MODELS(i)%ExportField(j)%etype=etype
            ELSE
              IF (localPET.eq.0) THEN
                WRITE (cplout,30) 'export field short_name: ',          &
     &                             TRIM(ShortName), TRIM(CPLname)
              END IF
              rc=ESMF_RC_NOT_FOUND
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) RETURN
            END IF
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Add IMPORT fields to NUOPC field dictionary.
!-----------------------------------------------------------------------
!
      DO i=1,Nmodels
        IF (MODELS(i)%IsActive) THEN
          DO j=1,Nimport(i)
            Sname=MODELS(i)%ImportField(j)%standard_name
            units=MODELS(i)%ImportField(j)%src_units
!
!  Check if field already exists.
!
            Exist=NUOPC_FieldDictionaryHasEntry(TRIM(Sname),            &
     &                                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Add field to dictionary.
!
            IF (.not.Exist) THEN
              CALL NUOPC_FieldDictionaryAddEntry(TRIM(Sname),           &
     &                                     canonicalUnits=TRIM(units),  &
     &                                           rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
          END DO
        END IF
      END DO
!
!  Report import fields metadata.
!
      IF ((localPET.eq.0).and.(MODELS(Iatmos )%IsActive.or.             &
     &                         MODELS(Idata  )%IsActive.or.             &
     &                         MODELS(Iseaice)%IsActive.or.             &
     &                         MODELS(Iwave  )%IsActive)) THEN

        WRITE (cplout,40) 'ESM IMPORT Fields Metadata Dictionary:'
!
        DO i=1,Nmodels
          Lreport=.FALSE.
          IF (MODELS(i)%IsActive) THEN
            Lreport=.TRUE.
            Smodel=TRIM(MODELS(i)%name)
          END IF
          IF (Lreport) THEN
            DO j=1,Nimport(i)
              WRITE (cplout,50)                                         &
     &              TRIM(Smodel),                                       &
     &              TRIM(MODELS(i)%ImportField(j)%short_name),          &
     &              TRIM(MODELS(i)%ImportField(j)%standard_name),       &
     &              MODELS(i)%ImportField(j)%gtype,                     &
     &              MODELS(i)%ImportField(j)%itype,                     &
     &              MODELS(i)%ImportField(j)%etype,                     &
     &              MODELS(i)%ImportField(j)%connected,                 &
     &              MODELS(i)%ImportField(j)%debug_write,               &
     &              MODELS(i)%ImportField(j)%add_offset,                &
     &              MODELS(i)%ImportField(j)%scale_factor
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Add EXPORT fields to NUOPC field dictionary.
!-----------------------------------------------------------------------
!
      DO i=1,Nmodels
        IF (MODELS(i)%IsActive) THEN
          DO j=1,Nexport(i)
            Sname=MODELS(i)%ExportField(j)%standard_name
            units=MODELS(i)%ExportField(j)%src_units
!
!  Check if field already exists.
!
            Exist=NUOPC_FieldDictionaryHasEntry(TRIM(Sname),            &
     &                                          rc=rc)
            IF (ESMF_LogFoundError(rcToCheck=rc,                        &
     &                             msg=ESMF_LOGERR_PASSTHRU,            &
     &                             line=__LINE__,                       &
     &                             file=MyFile)) THEN
              RETURN
            END IF
!
!  Add field to dictionary.
!
            IF (.not.Exist) THEN
              CALL NUOPC_FieldDictionaryAddEntry(TRIM(Sname),           &
     &                                     canonicalUnits=TRIM(units),  &
     &                                           rc=rc)
              IF (ESMF_LogFoundError(rcToCheck=rc,                      &
     &                               msg=ESMF_LOGERR_PASSTHRU,          &
     &                               line=__LINE__,                     &
     &                               file=MyFile)) THEN
                RETURN
              END IF
            END IF
          END DO
        END IF
      END DO
!
!  Report export fields metadata.
!
      IF ((localPET.eq.0).and.(MODELS(Iatmos )%IsActive.or.             &
     &                         MODELS(Idata  )%IsActive.or.             &
     &                         MODELS(Iseaice)%IsActive.or.             &
     &                         MODELS(Iwave  )%IsActive)) THEN

        WRITE (cplout,40) 'ESM EXPORT Fields Metadata Dictionary:'
!
        DO i=1,Nmodels
          Lreport=.FALSE.
          IF (MODELS(i)%IsActive) THEN
            Lreport=.TRUE.
            Smodel=TRIM(MODELS(i)%name)
          END IF
          IF (Lreport) THEN
            DO j=1,Nexport(i)
              WRITE (cplout,50)                                         &
     &              TRIM(Smodel),                                       &
     &              TRIM(MODELS(i)%ExportField(j)%short_name),          &
     &              TRIM(MODELS(i)%ExportField(j)%standard_name),       &
     &              MODELS(i)%ExportField(j)%gtype,                     &
     &              MODELS(i)%ExportField(j)%itype,                     &
     &              MODELS(i)%ExportField(j)%etype,                     &
     &              MODELS(i)%ExportField(j)%connected,                 &
     &              MODELS(i)%ExportField(j)%debug_write,               &
     &              MODELS(i)%ExportField(j)%add_offset,                &
     &              MODELS(i)%ExportField(j)%scale_factor
            END DO
          END IF
        END DO
      END IF
!
      IF (MODELS(Iatmos )%IsActive.or.                                  &
     &    MODELS(Idata  )%IsActive.or.                                  &
     &    MODELS(Iseaice)%IsActive.or.                                  &
     &    MODELS(Iwave  )%IsActive) THEN
        IF (localPET.eq.0) THEN
          WRITE (cplout,60)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Read in RunSequence input file and report.
!-----------------------------------------------------------------------
!
!  Open input coupling variable information file.
!
      OPEN (inp, FILE=TRIM(CONFname), FORM='formatted', STATUS='old',   &
     &      IOSTAT=io_err, IOMSG=io_errmsg)
      IF (io_err.ne.0) THEN
        IF (localPET.eq.0) WRITE(cplout,70) TRIM(CONFname),             &
     &                                      TRIM(io_errmsg)
        exit_flag=5
        rc=ESMF_RC_FILE_OPEN
        RETURN
      END IF
!
!  Read in and report RunSequence file. Ignore blank and comments
!  [char(35)=#] input lines.
!
      IF (localPET.eq.0) THEN
        WRITE (cplout,80) 'ESM Coupling RunSequence:', TRIM(CONFname)
      END IF
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        Lvar=LEN_TRIM(line)
        IF ((Lvar.gt.0).and.(line (1:1).ne.CHAR(35))) THEN
          IF (localPET.eq.0) WRITE (cplout,'(a)') TRIM(line)
        END IF
      END DO
  10  IF (localPET.eq.0) WRITE (cplout,90) TRIM(CONFname)
      rc=ESMF_RC_FILE_READ
      RETURN
  20  CLOSE (inp)
      IF (localPET.eq.0) WRITE (cplout,'(a)') CHAR(10)    ! new line
!
!  Flush standard output buffer.
!
      FLUSH (cplout)
!
!  Deallocate.
!
      IF (allocated(S))      deallocate (S)
      IF (allocated(Smodel)) deallocate (Smodel)
!
  30  FORMAT (/,' SET_METADATA - cannot find metadata for',             &
     &        1x,a,'''',a,'''.',/,16x,                                  &
     &        'Add entry to metadata file: ',a)
  40  FORMAT (/,a,/, 37('='),/,/, 'Model', t11,'Short Name',            &
     &        t25,'Standard Name', t74,'G', t77,'I', t80,'E' t83,'C',   &
     &        t86,'W',t90,'add_offset', t102,'scale_factor',/, 114('-'))
  50  FORMAT (a, t11,a, t25,a, t74,i1, t77,i1, t80,i1, t83,l1, t86,l1,  &
     &        t89,1p,e12.5, t103,1p,e12.5)
  60  FORMAT (/,' G: Grid cell location,     1=Center,',                &
     &                                     ' 2=Corner,',                &
     &                                     ' 3=U-point,',               &
     &                                     ' 4=V-point',                &
     &        /,' I: Regridding method,      1=bilinear,',              &
     &                                     ' 2=patch,',                 &
     &                                     ' 3=conservative1',          &
     &                                     ' 4=conservative2',          &
     &                                     ' 5=nearest',                &
     &        /,' E: Extrapolation method,   0=none,',                  &
     &                                     ' 1=nearest,',               &
     &                                     ' 2=inv distance avg,',      &
     &                                     ' 3=creep fill,',            &
     &                                     ' 4=two-steps (old way)',    &
     &        /,' C: Connected to coupler,   F=derived from other,',    &
     &                                     ' T=exchanged/regridded',    &
     &        /,' W: Field write to NetCDF,  F=false, T=true',          &
     &                                     ' (used if DebugLevel > 2)'/)
  70  FORMAT (/,' SET_METADATA - Unable to open RunSequence ',          &
     &          ' configuration file: ',/,16x,a,/,16x,'ERROR: ',a,      &
     &          /,16x,'Prototypes are located in the ESM directory.')
  80  FORMAT (/,a,2x,a,/, 25('='),/)
  90  FORMAT (/,' SET_METADATA - error while reading RunSequence',      &
     &        ' configuration file: ',a)
!
      RETURN
      END SUBROUTINE set_metadata
!
      SUBROUTINE def_FieldAtt (vm, rc)
!
!=======================================================================
!                                                                      !
!  If concurrent coupling and ESM components importing time snapshots, !
!  create output NetCDF file containing the field attributes needed to !
!  perform the time interpolation.                                     !
!                                                                      !
!=======================================================================
!
      USE mod_netcdf
!
      USE mod_parallel,   ONLY : OCN_COMM_WORLD
      USE mod_iounits,    ONLY : SourceFile, ioerror
      USE mod_ncparam,    ONLY : MaxLen
      USE mod_scalars,    ONLY : NoError, exit_flag
!
      USE def_dim_mod,    ONLY : def_dim
      USE def_var_mod,    ONLY : def_var
      USE distribute_mod, ONLY : mp_bcasti
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc

      TYPE (ESMF_VM) :: vm
!
!  Local variable declarations.
!
      integer :: CdimID, DdimID, FdimID, SdimID, SNdimID
      integer :: MyComm, ROMScomm, PETcount, localPET
      integer :: Ncomp, Nfields, i, j, ncid, status, varid
      integer :: ibuffer(3)
!
      integer, parameter :: Natt = 25           ! ROMS framework usage
      integer, parameter :: ng = 1              ! ROMS framework usage
      integer, parameter :: iKernel = 1         ! ROMS framework usage
!
      real(r8) :: Aval(6)
!
      character (len= 80)    :: text
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_FieldAtt"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the
!  mpi communicator handle and current node rank.
!-----------------------------------------------------------------------
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
!  Temporarily replace ROMS mpi communicator handle with that for DATA
!  component. Recall that the DATA component uses both ROMS NetCDF- and
!  mpi-frameworks.
!
      ROMScomm=OCN_COMM_WORLD
      OCN_COMM_WORLD=MyComm
!
!-----------------------------------------------------------------------
!  Create field time interpolation variable attributes NetCDF file.
!-----------------------------------------------------------------------
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
      Nfields=SUM(Nexport)
!
!  Create NetCDF file.
!
      CALL netcdf_create (ng, iKernel, AttFileName, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (localPET.eq.0) WRITE (cplout,10) TRIM(AttFileName)
        rc=ESMF_RC_FILE_CREATE
        RETURN
      END IF
!
!  Define dimensions.
!
      status=def_dim(ng, iKernel, ncid, AttFileName, 'Nmodels',         &
     &               Nmodels, CdimID)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      status=def_dim(ng, iKernel, ncid, AttFileName, 'field',           &
     &               Nfields, FdimID)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!                                                  length of short-name
      status=def_dim(ng, iKernel, ncid, AttFileName, 'char20',          &
     &               20, SNdimID)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!                                                  length snapshot date
      status=def_dim(ng, iKernel, ncid, AttFileName, 'char22',          &
     &               22, DdimID)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      status=def_dim(ng, iKernel, ncid, AttFileName, 'snapshot',        &
     &               nf90_unlimited, SdimID)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Define global attributes.
!
      text='ROMS ESMF/NUOPC coupling, '//                               &
     &     'field metadata for time interpolation between snapshots'
      IF (localPET.eq.0) THEN
        status=nf90_put_att(ncid, nf90_global, 'type',                  &
     &                      TRIM(text))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          WRITE (cplout,20) 'type', TRIM(AttFileName)
          rc=ESMF_RC_FILE_CREATE
          exit_flag=3
          ioerror=status
        END IF
      END IF
      ibuffer(1)=rc
      ibuffer(2)=exit_flag
      ibuffer(3)=ioerror
      CALL mp_bcasti (ng, iKernel, ibuffer)
      rc=ibuffer(1)
      exit_flag=ibuffer(2)
      ioerror=ibuffer(3)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) RETURN
!
      IF (localPET.eq.0) THEN
        status=nf90_put_att(ncid, nf90_global, 'file',                  &
     &                      TRIM(AttFileName))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          WRITE (cplout,20) 'file', TRIM(AttFileName)
          rc=ESMF_RC_FILE_CREATE
          exit_flag=3
          ioerror=status
        END IF
      END IF
      ibuffer(1)=rc
      ibuffer(2)=exit_flag
      ibuffer(3)=ioerror
      CALL mp_bcasti (ng, iKernel, ibuffer)
      rc=ibuffer(1)
      exit_flag=ibuffer(2)
      ioerror=ibuffer(3)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) RETURN
!
      text='CF-1.4'
      IF (localPET.eq.0) THEN
        status=nf90_put_att(ncid, nf90_global, 'Conventions',           &
     &                    TRIM(text))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          WRITE (cplout,20) 'Conventions', TRIM(AttFileName)
          rc=ESMF_RC_FILE_CREATE
          exit_flag=3
          ioerror=status
        END IF
      END IF
      ibuffer(1)=rc
      ibuffer(2)=exit_flag
      ibuffer(3)=ioerror
      CALL mp_bcasti (ng, iKernel, ibuffer)
      rc=ibuffer(1)
      exit_flag=ibuffer(2)
      ioerror=ibuffer(3)
      IF (ESMF_LogFoundError(rcToCheck=rc,                              &
     &                       msg=ESMF_LOGERR_PASSTHRU,                  &
     &                       line=__LINE__,                             &
     &                       file=MyFile)) RETURN
!
!  Define variables.
!
      Vinfo(1)='field'
      Vinfo(2)='field short name'
      status=def_var(ng, iKernel, ncid, varid, nf90_char,               &
     &               3, (/SNdimID,CdimID,FdimID/), Aval, Vinfo,         &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      Vinfo(1)='Tindex'
      Vinfo(2)='current time snapshot rolling index'
      status=def_var(ng, iKernel, ncid, varid, nf90_int,                &
     &               3, (/CdimID,FdimID,SdimID/), Aval, Vinfo,          &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      Vinfo(1)='Date'
      Vinfo(2)='field snapshot date'
      status=def_var(ng, iKernel, ncid, varid, nf90_char,               &
     &               4, (/DdimID,CdimID,FdimID,SdimID/), Aval, Vinfo,   &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      Vinfo(1)='Tcurrent'
      Vinfo(2)='current coupling time'
      Vinfo(3)='days since '//                                          &
     &         TRIM(ClockInfo(Idriver)%Time_ReferenceString)
      status=def_var(ng, iKernel, ncid, varid, NF_TOUT,                 &
     &               3, (/CdimID,FdimID,SdimID/), Aval, Vinfo,          &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      Vinfo(1)='Tstr'
      Vinfo(2)='field lower time snapshot'
      Vinfo(3)='days since '//                                          &
     &         TRIM(ClockInfo(Idriver)%Time_ReferenceString)
      status=def_var(ng, iKernel, ncid, varid, NF_TOUT,                 &
     &               3, (/CdimID,FdimID,SdimID/), Aval, Vinfo,          &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      Vinfo(1)='Tend'
      Vinfo(2)='field upper time snapshot'
      status=def_var(ng, iKernel, ncid, varid, NF_TOUT,                 &
     &               3, (/CdimID,FdimID,SdimID/), Aval, Vinfo,          &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      Vinfo(1)='Tintrp'
      Vinfo(2)='field monotonic time for interpolation'
      Vinfo(3)='days since '//                                          &
     &         TRIM(ClockInfo(Idriver)%Time_ReferenceString)
      status=def_var(ng, iKernel, ncid, varid, NF_TOUT,                 &
     &               3, (/CdimID,FdimID,SdimID/), Aval, Vinfo,          &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      Vinfo(1)='Vtime'
      Vinfo(2)='field variable current time snapshots'
      Vinfo(3)='days since '//                                          &
     &         TRIM(ClockInfo(Idriver)%Time_ReferenceString)
      status=def_var(ng, iKernel, ncid, varid, NF_TOUT,                 &
     &               3, (/CdimID,FdimID,SdimID/), Aval, Vinfo,          &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      Vinfo(1)='Tmin'
      Vinfo(2)='field time minimum in data'
      Vinfo(3)='days since '//                                          &
     &         TRIM(ClockInfo(Idriver)%Time_ReferenceString)
      status=def_var(ng, iKernel, ncid, varid, NF_TOUT,                 &
     &               3, (/CdimID,FdimID,SdimID/), Aval, Vinfo,          &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      Vinfo(1)='Tmax'
      Vinfo(2)='field time maximum in data'
      Vinfo(3)='days since '//                                          &
     &         TRIM(ClockInfo(Idriver)%Time_ReferenceString)
      status=def_var(ng, iKernel, ncid, varid, NF_TYPE,                 &
     &               3, (/CdimID,FdimID,SdimID/), Aval, Vinfo,          &
     &               AttFileName, SetParAccess = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Leave definition mode.
!
      CALL netcdf_enddef (ng, iKernel, AttFileName, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        RETURN
      END IF
!
!  Close file to allow ESM components to manipulate the data across
!  different VM in concurrent coupling. The file will be opened and
!  closed by during data input and output processing allowing
!  synchronization and access to data immidiately after is written.
!
      CALL netcdf_close (ng, iKernel, ncid, AttFileName)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        rc=ESMF_RC_FILE_CREATE
        RETURN
      END IF
!
!  Restore ROMS mpi communicator handle.
!
      OCN_COMM_WORLD=ROMScomm
!
  10  FORMAT (/,' def_FieldAtt - unable to create NetCDF file: ',a)
  20  FORMAT (/,' def_FieldAtt - error while creating global',          &
     &          ' attribute: ', a,/,16x,a)
!
      RETURN
      END SUBROUTINE def_FieldAtt
!
      SUBROUTINE get_weights (Icomp, Ix, Iy, vm, rc)
!
!=======================================================================
!                                                                      !
!  This routine Sets or reads in melding weights coefficients needed   !
!  by the atmosphere model to merge fields from DATA and other ESM     !
!  components because of incongruent grids.                            !
!                                                                      !
!=======================================================================
!
      USE mod_netcdf
!
      USE mod_parallel,   ONLY : OCN_COMM_WORLD
      USE mod_iounits,    ONLY : SourceFile
      USE mod_scalars,    ONLY : NoError, exit_flag
!
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: Icomp, Ix, Iy
      integer, intent(out) :: rc
!
      TYPE (ESMF_VM) :: vm
!
!  Local variable declarations.
!
      integer :: MyComm, ROMScomm, PETcount, localPET
      integer :: Nx, Ny
!
      integer, parameter :: ng = 1              ! ROMS framework usage
      integer, parameter :: iKernel = 1         ! ROMS framework usage
!
      real(r8) :: Vmin, Vmax
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_weights"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the
!  mpi communicator handle and current node rank.
!-----------------------------------------------------------------------
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
!  Temporarily replace ROMS mpi communicator handle with that for DATA
!  component. Recall that the DATA component uses both ROMS NetCDF- and
!  mpi-frameworks.
!
      ROMScomm=OCN_COMM_WORLD
      OCN_COMM_WORLD=MyComm
!
!-----------------------------------------------------------------------
!  Read in melding weighting coefficients.
!-----------------------------------------------------------------------
!
      IF (MODELS(Icomp)%IsActive) THEN
        ncname=TRIM(WEIGHTS(Icomp)%ncfile)
!
!  Inquire about the weights variables.
!
        CALL netcdf_inq_var (ng, iKernel, ncname,                       &
     &                       MyVarName = TRIM(WEIGHTS(Icomp)%VnameDATA))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (localPET.eq.0) THEN
            WRITE (cplout,10) TRIM(WEIGHTS(Icomp)%VnameDATA),           &
     &                        TRIM(ncname)
          END IF
          rc=ESMF_RC_FILE_READ
          RETURN
        ELSE
          Nx=var_Dsize(1)
          Ny=var_Dsize(2)
          IF (Ix.ne.Nx) THEN
            IF (localPET.eq.0) THEN
              WRITE (cplout,20) TRIM(var_Dname(1)), Nx, Ix, TRIM(ncname)
            END IF
            rc=ESMF_RC_NOT_VALID
            RETURN
          END IF
          IF (Iy.ne.Ny) THEN
            IF (localPET.eq.0) THEN
              WRITE (cplout,20) TRIM(var_Dname(2)), Ny, Iy, TRIM(ncname)
            END IF
            rc=ESMF_RC_NOT_VALID
            RETURN
          END IF
        END IF
!
!  Allocate weights arrays.
!
        IF (.not.allocated(WEIGHTS(Icomp)%Cdat)) THEN
          allocate ( WEIGHTS(Icomp)%Cdat(Ix,Iy) )
          WEIGHTS(Icomp)%Cdat = 0.0_r8
        END IF
        IF (.not.allocated(WEIGHTS(Icomp)%Cesm)) THEN
          allocate ( WEIGHTS(Icomp)%Cesm(Ix,Iy) )
          WEIGHTS(Icomp)%Cesm = 0.0_r8
        END IF
!
!  Read in weights and broadcast full arrays to all the members of the
!  communicator.
!
        CALL netcdf_get_fvar (ng, iKernel, ncname,                      &
     &                        TRIM(WEIGHTS(Icomp)%VnameDATA),           &
     &                        WEIGHTS(Icomp)%Cdat,                      &
     &                        start = (/1,1/),                          &
     &                        total = (/Ix,Iy/),                        &
     &                        min_val = Vmin,                           &
     &                        max_val = Vmax)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (localPET.eq.0) THEN
            WRITE (cplout,30) TRIM(WEIGHTS(Icomp)%VnameDATA),           &
     &                        TRIM(ncname)
          END IF
          rc=ESMF_RC_FILE_READ
          RETURN
        ELSE
          IF (localPET.eq.0) THEN
            WRITE (cplout,40) TRIM(WEIGHTS(Icomp)%VnameDATA),           &
     &                        Vmin, Vmax
          END IF
        END IF
!
        CALL netcdf_get_fvar (ng, iKernel, ncname,                      &
     &                        TRIM(WEIGHTS(Icomp)%VnameESM),            &
     &                        WEIGHTS(Icomp)%Cesm,                      &
     &                        start = (/1,1/),                          &
     &                        total = (/Ix,Iy/),                        &
     &                        min_val = Vmin,                           &
     &                        max_val = Vmax)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (localPET.eq.0) THEN
            WRITE (cplout,30) TRIM(WEIGHTS(Icomp)%VnameESM),            &
     &                        TRIM(ncname)
          END IF
          rc=ESMF_RC_FILE_READ
          RETURN
        ELSE
          IF (localPET.eq.0) THEN
            WRITE (cplout,40) TRIM(WEIGHTS(Icomp)%VnameESM),            &
     &                        Vmin, Vmax
          END IF
        END IF
      END IF
!
!  Restore ROMS mpi communicator handle.
!
      OCN_COMM_WORLD=ROMScomm
!
  10  FORMAT (/,' GET_WEIGHTS - error while inquiring variable: ',a,    &
     &        /,15x,'in file: ''',a,'''')
  20  FORMAT (/,' GET_WEIGHTS - incorrect value for dimension: ',a,     &
     &        ' = ',i0,', need ',i0,' instead'/,15x,'in file: ',a)
  30  FORMAT (/,' GET_WEIGHTS - error while reading variable: ''',a,    &
     &        '''',/,15x,'in file: ',a)
  40  FORMAT (3x,' GET_WEIGHTS - ESMF: reading ''',a,'''',/,19x,        &
     &        '(Wmin= ', 1p,e15.8,0p,' Wmax= ',1p,e15.8,0p,')')
!
      RETURN
      END SUBROUTINE get_weights
#endif
      END MODULE mod_esmf_esm
