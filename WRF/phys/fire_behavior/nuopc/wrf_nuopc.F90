
#include "Fire_behavior_NUOPC_Macros.h"

module wrf_nuopc

  !-----------------------------------------------------------------------------
  ! offline WRF component
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS    => SetServices
  use wrfdata_mod, only : wrfdata_t
  use namelist_mod, only : namelist_t
  use datetime_mod, only : datetime_t
  use initialize_mod, only: Init_atm_state
  use stderrout_mod, only : Print_message, Stop_simulation

  implicit none

  private

  public SetVM, SetServices

  type (wrfdata_t) :: state
  real(ESMF_KIND_R8), pointer     :: ptr_z0(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_q2(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_psfc(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_rain(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_t2(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_u10(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_v10(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_u3d(:,:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_v3d(:,:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_phl(:,:,:)
!  real(ESMF_KIND_R8), pointer     :: ptr_pres(:,:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_hflx_fire(:,:)
  integer                         :: clb(2), cub(2), clb3(3), cub3(3)
  integer ::  igs, jgs, ige, jge
  type (namelist_t) :: config_flags

  logical, parameter :: DEBUG_ALL = .false.

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(model, rc)

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    logical, parameter :: DEBUG_LOCAL = .false.


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering SetServices wrf-data...')

    rc = ESMF_SUCCESS

    ! derive from NUOPC_Model
    call NUOPC_CompDerive(model, modelSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! specialize model
    call NUOPC_CompSpecialize(model, specLabel=label_Advertise, &
      specRoutine=Advertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_RealizeProvided, &
      specRoutine=Realize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize (model, specLabel = label_SetClock, specRoutine = SetClock, rc = rc)
    if (ESMF_LogFoundError (rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line = __LINE__, file = __FILE__)) &
        return

    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=Advance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Leaving SetServices wrf-data...')

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Advertise(model, rc)

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)        :: importState, exportState
    logical, parameter :: DEBUG_LOCAL = .false.


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering Advertise wrf-data...')

    rc = ESMF_SUCCESS

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Eventually, you will advertise your model's import and
    ! export fields in this phase.  For now, however, call
    ! your model's initialization routine(s).

    call config_flags%Init_time_block ('namelist.fire')
    call config_flags%Init_atm_block ('namelist.fire')

    call Init_atm_state(state, config_flags)

    if (.not. allocated (state%q2)) &
        allocate (state%q2(size(state%lats, dim=1), size(state%lats, dim=2)))

    if (.not. allocated (state%t2)) &
        allocate (state%t2(size(state%lats, dim=1), size(state%lats, dim=2)))

    if (.not. allocated (state%u10)) &
        allocate (state%u10(size(state%lats, dim=1), size(state%lats, dim=2)))

    if (.not. allocated (state%v10)) &
        allocate (state%v10(size(state%lats, dim=1), size(state%lats, dim=2)))

    if (.not. allocated (state%z0)) &
        allocate (state%z0(size(state%lats, dim=1), size(state%lats, dim=2)))

    if (.not. allocated (state%psfc)) &
        allocate (state%psfc(size(state%lats, dim=1), size(state%lats, dim=2)))

    if (.not. allocated (state%rain)) &
        allocate (state%rain(size(state%lats, dim=1), size(state%lats, dim=2)))

    allocate (state%u3d(size(state%lats, dim=1), size(state%lats, dim=2), state%kde - 1))
    allocate (state%v3d(size(state%lats, dim=1), size(state%lats, dim=2), state%kde - 1))
    allocate (state%phl(size(state%lats, dim=1), size(state%lats, dim=2), state%kde - 1))
!    allocate (state%pres(size(state%lats, dim=1), size(state%lats, dim=2), state%kde - 1))

    ! Import/ Export Variables -----------------------------------------------------

    ! Disabling the following macro, e.g. renaming to WITHIMPORTFIELDS_disable,
    ! will result in a model component that does not advertise any importable
    ! Fields. Use this if you want to drive the model independently.

#define WITHEXPORTFIELDS
#ifdef WITHEXPORTFIELDS
    ! 3D
    ! exportable field: inst_zonal_wind_levels
    call NUOPC_Advertise(exportState, &
      StandardName="inst_zonal_wind_levels", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_merid_wind_levels
    call NUOPC_Advertise(exportState, &
      StandardName="inst_merid_wind_levels", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_geop_levels
    call NUOPC_Advertise(exportState, &
      StandardName="inst_geop_levels", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    ! exportable field: inst_pres_levels
!    call NUOPC_Advertise(exportState, &
!      StandardName="inst_pres_levels", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    ! 2D
    ! exportable field: inst_zonal_wind_levels
    call NUOPC_Advertise(exportState, &
      StandardName="inst_zonal_wind_height10m", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_merid_wind_levels
    call NUOPC_Advertise(exportState, &
      StandardName="inst_merid_wind_height10m", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_surface_roughness
    call NUOPC_Advertise(exportState, &
      StandardName="inst_surface_roughness", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_spec_humid_height2m
    call NUOPC_Advertise(exportState, &
      StandardName="inst_spec_humid_height2m", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_pres_height_surface
    call NUOPC_Advertise(exportState, &
      StandardName="inst_pres_height_surface", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: accumulated_lwe_thickness_of_precipitation_amount
    ! rainfall is accumulated rainfall not instantaneous
    call NUOPC_Advertise(exportState, &
      StandardName="accumulated_lwe_thickness_of_precipitation_amount", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_temp_height2m
    call NUOPC_Advertise(exportState, &
      StandardName="inst_temp_height2m", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_pres_height_lowest_from_phys
    call NUOPC_Advertise(exportState, &
      StandardName="inst_pres_height_lowest_from_phys", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_spec_humid_height_lowest_from_phys
    call NUOPC_Advertise(exportState, &
      StandardName="inst_spec_humid_height_lowest_from_phys", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: inst_temp_height_lowest_from_phys
    call NUOPC_Advertise(exportState, &
      StandardName="inst_temp_height_lowest_from_phys", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Advertise(importState, &
      StandardName="hflx_fire", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Advertise(importState, &
      StandardName="evap_fire", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Advertise(importState, &
      StandardName="smoke_fire", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Leaving Advertise wrf-data...')

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Realize(model, rc)

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)        :: importState, exportState
    ! type(ESMF_TimeInterval) :: stabilityTimeStep
    type(ESMF_Field)        :: field
    type(ESMF_DistGrid)     :: distgrid
    type(ESMF_Grid)         :: grid

    type(ESMF_VM)                   :: vm
    integer :: localDE, petCount, localPet, localDECount

    ! working local variables
    integer                        :: lbnd(2),ubnd(2)
    integer, allocatable :: minIndexPDe(:,:)
    integer, allocatable :: maxIndexPDe(:,:)
    integer, allocatable :: localDeToDeMap(:)
    integer :: dimCount, deCount
    real(ESMF_KIND_COORD), pointer :: coordXcenter(:,:)
    real(ESMF_KIND_COORD), pointer :: coordYcenter(:,:)
    real(ESMF_KIND_COORD), pointer :: coordXcorner(:,:)
    real(ESMF_KIND_COORD), pointer :: coordYcorner(:,:)
    integer                        :: i, j, nx, ny
    type (datetime_t) :: datetime_now
    logical, parameter :: DEBUG_LOCAL = .false.
    integer :: gDE, iglobal, jglobal


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering Realize wrf-data...')

    rc = ESMF_SUCCESS

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)

    nx = size(state%lats, dim=1)
    ny = size(state%lons, dim=2)

    ! Create distgrid based on the state grid
    distgrid = ESMF_DistGridCreate( &
      minIndex=(/1,1/), maxIndex=(/nx,ny/), &
      rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    grid = ESMF_GridCreate(name='ATM', &
      distgrid=distgrid, coordSys = ESMF_COORDSYS_SPH_DEG, &
      rc = rc)
    if(ESMF_STDERRORCHECK(rc)) return

    call ESMF_GridGet(grid, localDECount=localDECount, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return
    allocate (localDeToDeMap(localDECount))

    call ESMF_DistGridGet(distgrid, dimCount=dimCount, deCount=deCount, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return
    allocate (minIndexPDe(dimCount, deCount))
    allocate (maxIndexPDe(dimCount, deCount))

    call ESMF_DistGridGet(distgrid, localDeToDeMap=localDeToDeMap, &
        minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    ! Add Center Coordinates to Grid
    call ESMF_GridAddCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    localDE = 0
!do localDE = 0, localDECount - 1
    call ESMF_GridGetCoord(grid, coordDim=1, localDE=localDE, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        computationalLBound=lbnd, computationalUBound=ubnd, &
        farrayPtr=coordXcenter, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_GridGetCoord(grid, coordDim=2, localDE=localDE, &
        staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=coordYcenter, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    gDE = localDeToDeMap(localDE + 1)
    igs = minIndexPDe(1, gDE + 1)
    jgs = minIndexPde(2, gDE + 1)
    ige = maxIndexPDe(1, gDE + 1)
    jge = maxIndexPde(2, gDE + 1)
    do j = lbnd(2), ubnd(2)
      jglobal = jgs + (j - lbnd(2))
      do i = lbnd(1), ubnd(1)
        iglobal = igs + (i - lbnd(1))
        coordXcenter(i,j) = state%lons(iglobal, jglobal)
        coordYcenter(i,j) = state%lats(iglobal, jglobal)
      end do
    end do
!end do

    ! CORNERS

    ! Add Corner Coordinates to Grid
    call ESMF_GridAddCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

!do localDE = 0, localDECount - 1
    call ESMF_GridGetCoord(grid, coordDim=1, localDE=localDE, &
      staggerloc=ESMF_STAGGERLOC_CORNER, &
      computationalLBound=lbnd, computationalUBound=ubnd, &
      farrayPtr=coordXcorner, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    call ESMF_GridGetCoord(grid, coordDim=2, localDE=localDE, &
      staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=coordYcorner, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    gDE = localDeToDeMap(localDE + 1)

    do j = lbnd(2), ubnd(2)
      jglobal = minIndexPDe(2, gDE + 1) + (j - lbnd(2))
      do i = lbnd(1), ubnd(1)
        iglobal = minIndexPDe(1, gDE + 1) + (i - lbnd(1))
        coordXcorner(i, j) = state%lons_c(iglobal, jglobal)
        coordYcorner(i, j) = state%lats_c(iglobal, jglobal)
      end do
    end do
!end do

    ! importable field on Grid: hflx_fire
    field = ESMF_FieldCreate(name="hflx_fire", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!do localDE = 0, localDECount - 1
    call ESMF_FieldGet(field, localDe=localDE, farrayPtr=ptr_hflx_fire, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

    ! importable field on Grid: evap_fire
    field = ESMF_FieldCreate(name="evap_fire", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_evap_fire, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    ! importable field on Grid: smoke_fire
    field = ESMF_FieldCreate(name="smoke_fire", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_evap_fire, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out


#ifdef WITHEXPORTFIELDS

     ! exportable field on Grid: inst_zonal_wind_levels
     field = ESMF_FieldCreate(name="inst_zonal_wind_levels", grid=grid, &
       gridToFieldMap=(/1,2/), ungriddedLBound=(/1/), &
       ungriddedUBound=(/state%kde - 1/), &
       typekind=ESMF_TYPEKIND_R8, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
     call NUOPC_Realize(exportState, field=field, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
     ! Get Field memory
!do localDE = 0, localDECount - 1
     call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_u3d, &
       computationalLBound=clb3, computationalUBound=cub3, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

     ! exportable field on Grid: inst_merid_wind_levels
     field = ESMF_FieldCreate(name="inst_merid_wind_levels", grid=grid, &
       gridToFieldMap=(/1,2/), ungriddedLBound=(/1/), &
       ungriddedUBound=(/state%kde - 1/), &
       typekind=ESMF_TYPEKIND_R8, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
     call NUOPC_Realize(exportState, field=field, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
     ! Get Field memory
!do localDE = 0, localDECount - 1
     call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_v3d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

     ! exportable field on Grid: inst_geop_levels
     field = ESMF_FieldCreate(name="inst_geop_levels", grid=grid, &
       gridToFieldMap=(/1,2/), ungriddedLBound=(/1/), &
       ungriddedUBound=(/state%kde - 1/), &
       typekind=ESMF_TYPEKIND_R8, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
     call NUOPC_Realize(exportState, field=field, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
     ! Get Field memory
!do localDE = 0, localDECount - 1
     call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_phl, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

     ! exportable field on Grid: inst_pres_levels
!     field = ESMF_FieldCreate(name="inst_pres_levels", grid=grid, &
!       gridToFieldMap=(/1,2/), ungriddedLBound=(/1/), &
!       ungriddedUBound=(/state%kde - 1/), &
!       typekind=ESMF_TYPEKIND_R8, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!       line=__LINE__, &
!       file=__FILE__)) &
!       return  ! bail out
!     call NUOPC_Realize(exportState, field=field, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!       line=__LINE__, &
!       file=__FILE__)) &
!       return  ! bail out
!     ! Get Field memory
!     call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_pres, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    ! 2D
    ! exportable field on Grid: inst_surface_roughness
    field = ESMF_FieldCreate(name="inst_surface_roughness", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!do localDE = 0, localDECount - 1
    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_z0, &
      computationalLBound=clb, computationalUBound=cub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

    ! exportable field on Grid: inst_spec_humid_height2m
    field = ESMF_FieldCreate(name="inst_spec_humid_height2m", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!do localDE = 0, localDECount - 1
    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_q2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

    ! exportable field on Grid: inst_pres_height_surface
    field = ESMF_FieldCreate(name="inst_pres_height_surface", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!do localDE = 0, localDECount - 1
    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_psfc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

    ! exportable field on Grid: accumulated_lwe_thickness_of_precipitation_amount
    field = ESMF_FieldCreate(name="accumulated_lwe_thickness_of_precipitation_amount", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!do localDE = 0, localDECount - 1
    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_rain, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

    ! exportable field on Grid: inst_temp_height2m
    field = ESMF_FieldCreate(name="inst_temp_height2m", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!do localDE = 0, localDECount - 1
    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_t2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

    ! exportable field on Grid: inst_zonal_wind_height10m
    field = ESMF_FieldCreate(name="inst_zonal_wind_height10m", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!do localDE = 0, localDECount - 1
    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_u10, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

    ! exportable field on Grid: inst_merid_wind_height10m
    field = ESMF_FieldCreate(name="inst_merid_wind_height10m", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Get Field memory
!do localDE = 0, localDECount - 1
    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=ptr_v10, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!end do

    ! exportable field on Grid: inst_pres_height_lowest_from_phys
    field = ESMF_FieldCreate(name="inst_pres_height_lowest_from_phys", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field on Grid: inst_spec_humid_height_lowest_from_phys
    field = ESMF_FieldCreate(name="inst_spec_humid_height_lowest_from_phys", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field on Grid: inst_temp_height_lowest_from_phys
    field = ESMF_FieldCreate(name="inst_temp_height_lowest_from_phys", grid=grid, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    datetime_now = datetime_t (config_flags%start_year, config_flags%start_month, &
      config_flags%start_day, config_flags%start_hour, config_flags%start_minute, &
      config_flags%start_second)

      ! "Initialize" atmospheric model
     call Update_atm_state(datetime_now)
#endif

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Check_grid_cells(grid)

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Leaving Realize wrf-data...')

  end subroutine

  subroutine SetClock (model, rc)

    implicit none

    type(ESMF_GridComp) :: model
    integer, intent(out) :: rc

    type (ESMF_Clock) :: modelClock
    type (ESMF_Time) :: startTime
    type (ESMF_Time) :: stopTime
    type (ESMF_TimeInterval) :: timeStep
    logical, parameter :: DEBUG_LOCAL = .false.


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering SetClock wrf-data...')

    rc = ESMF_SUCCESS

    call ESMF_TimeIntervalSet (timeStep, s = config_flags%interval_atm, rc = rc)
    if (ESMF_LogFoundError (rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, line = __LINE__, file = __FILE__)) &
        return

    call ESMF_TimeSet (startTime, yy = config_flags%start_year, mm = config_flags%start_month, &
        dd = config_flags%start_day, h = config_flags%start_hour, m = config_flags%start_minute, &
        s = config_flags%start_second, calkindflag = ESMF_CALKIND_GREGORIAN, rc = rc)
    if (ESMF_LogFoundError (rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, line = __LINE__, file = __FILE__)) &
        return

    call ESMF_TimeSet (stopTime, yy = config_flags%end_year, mm = config_flags%end_month, &
        dd = config_flags%end_day, h = config_flags%end_hour, m = config_flags%end_minute, &
        s = config_flags%end_second, calkindflag = ESMF_CALKIND_GREGORIAN, rc = rc)
    if (ESMF_LogFoundError (rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, line = __LINE__, file = __FILE__)) &
        return

    modelClock = ESMF_ClockCreate (name = "WRFdata Clock", timeStep = timeStep, startTime = startTime, stopTime = stopTime, rc = rc)
    if (ESMF_LogFoundError (rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, line = __LINE__, file = __FILE__)) &
        return

    call ESMF_GridCompSet (model, clock = modelClock, rc = rc)
    if (ESMF_LogFoundError (rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, line = __LINE__, file = __FILE__)) &
        return

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Leaving SetClock wrf-data...')

  end subroutine

  subroutine Advance(model, rc)

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: nextTime
!    type(ESMF_Time)             :: currTime
    ! type(ESMF_TimeInterval)     :: timeStep
    ! type(ESMF_VM)               :: vm
    ! integer                     :: currentSsiPe
     character(len=160)          :: msgString
    integer :: yy_now, mm_now, dd_now, h_now, m_now, s_now, ms
    real (kind = ESMF_KIND_R8) :: s_r8
    type (datetime_t) :: datetime_now
    logical, parameter :: DEBUG_LOCAL = .false.


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering Advance wrf-data...')

    rc = ESMF_SUCCESS

    ! query for clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

      ! Get current time stamp
!    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet(time=nextTime, yy  = yy_now, mm = mm_now, dd = dd_now, h = h_now, m = m_now, s = s_now, s_r8 = s_r8, ms = ms, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    datetime_now = datetime_t (yy_now, mm_now, dd_now, h_now, m_now, s_now)
    call datetime_now%Print_datetime ()

    ! ! Query for VM
    ! call ESMF_GridCompGet(model, vm=vm, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !   line=__LINE__, &
    !   file=__FILE__)) &
    !   return  ! bail out

    ! call ESMF_VMLog(vm, "LUMO Advance(): ", ESMF_LOGMSG_INFO, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !   line=__LINE__, &
    !   file=__FILE__)) &
    !   return  ! bail out

!      ! "Run" atmospheric model
     call Update_atm_state(datetime_now)
!    call state%Destroy_t2 ()

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing WRFdata model from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockPrint(clock, options="stopTime", &
      preString="---------------------> to: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Leaving Advance wrf-data...')

  end subroutine

  subroutine Update_atm_state (datetime)

    type (datetime_t), intent (in) :: datetime

    logical, parameter :: DEBUG_LOCAL = .false.
    character (len = 256) :: message


    call state%Get_z0(datetime)
    call state%Get_t2(datetime)
    call state%Get_q2(datetime)
    call state%Get_psfc(datetime)
    call state%Get_rain(datetime)
    call state%Get_u10(datetime)
    call state%Get_v10(datetime)
    call state%Get_u3d(datetime)
    call state%Get_v3d(datetime)
    call state%Get_phl(datetime)
!    call state%Get_pres(datetime)

    if (DEBUG_LOCAL .or. DEBUG_ALL) then
      write (message, *) 'clb(1), cub(1) = ', clb(1), cub(1)
      call Print_message (message)

      write (message, *) 'clb(2), cub(2) = ', clb(2), cub(2)
      call Print_message (message)

      write (message, *) 'igs, ige = ', igs, ige
      call Print_message (message)

      write (message, *) 'jgs, jge = ', jgs, jge
      call Print_message (message)
    end if

      ! 2D vars
      ! convert z0 from [m] to [cm]
    ptr_z0(clb(1):cub(1),clb(2):cub(2)) = state%z0(igs:ige,jgs:jge) * 100.0
    ptr_q2(clb(1):cub(1),clb(2):cub(2)) = state%q2(igs:ige, jgs:jge)
    ptr_psfc(clb(1):cub(1),clb(2):cub(2)) = state%psfc(igs:ige, jgs:jge)
    ptr_rain(clb(1):cub(1),clb(2):cub(2)) = state%rain(igs:ige, jgs:jge)
    ptr_t2(clb(1):cub(1),clb(2):cub(2)) = state%t2(igs:ige, jgs:jge)
    ptr_u10(clb(1):cub(1),clb(2):cub(2)) = state%u10(igs:ige, jgs:jge)
    ptr_v10(clb(1):cub(1),clb(2):cub(2)) = state%v10(igs:ige, jgs:jge)

      ! 3D vars
    ptr_u3d(clb3(1):cub3(1),clb3(2):cub3(2),clb3(3):cub3(3))= &
        state%u3d(igs:ige, jgs:jge, 1:state%kde - 1)
    ptr_v3d(clb3(1):cub3(1),clb3(2):cub3(2),clb3(3):cub3(3))= &
        state%v3d(igs:ige, jgs:jge, 1:state%kde - 1)
    ptr_phl(clb3(1):cub3(1),clb3(2):cub3(2),clb3(3):cub3(3))= &
        state%phl(igs:ige, jgs:jge, 1:state%kde - 1)

!    ptr_pres(clb3(1):cub3(1),clb3(2):cub3(2),clb3(3):cub3(3))= &
!        state%pres(igs:ige, jgs:jge, 1:state%kde - 1)

  end subroutine

  subroutine Check_grid_cells(grid)

    implicit none

    type(ESMF_Grid), intent(in) :: grid
    integer :: rc
    integer :: localDE, localDECount
    integer :: lbnd(2), ubnd(2)
    real(ESMF_KIND_COORD), pointer :: Xc(:,:), Yc(:,:)        ! centers
    real(ESMF_KIND_COORD), pointer :: Xcor(:,:), Ycor(:,:)    ! corners
    integer :: i, j
    real :: xmin, xmax, ymin, ymax
    logical :: bad
    character (len = 256) :: message


    rc = ESMF_SUCCESS

    call ESMF_GridGet(grid, localDECount=localDECount, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    do localDE = 0, localDECount - 1
      call ESMF_GridGetCoord(grid, coordDim=1, localDE=localDE, &
          staggerloc=ESMF_STAGGERLOC_CENTER, &
          computationalLBound=lbnd, computationalUBound=ubnd, &
          farrayPtr=Xc, rc=rc)

      call ESMF_GridGetCoord(grid, coordDim=2, localDE=localDE, &
          staggerloc=ESMF_STAGGERLOC_CENTER, &
          farrayPtr=Yc, rc=rc)

      call ESMF_GridGetCoord(grid, coordDim=1, localDE=localDE, &
          staggerloc=ESMF_STAGGERLOC_CORNER, &
          farrayPtr=Xcor, rc=rc)

      call ESMF_GridGetCoord(grid, coordDim=2, localDE=localDE, &
          staggerloc=ESMF_STAGGERLOC_CORNER, &
          farrayPtr=Ycor, rc=rc)

      write (message, *) 'Checking localDE =', localDE
      call Print_message (trim (message))

      do j = lbnd(2), ubnd(2) - 1
        do i = lbnd(1), ubnd(1) - 1
            ! Bounding box from 4 corners
          xmin = min (Xcor(i, j), Xcor(i + 1,j), Xcor(i, j + 1), Xcor(i + 1, j + 1))
          xmax = max (Xcor(i, j), Xcor(i + 1,j), Xcor(i, j + 1), Xcor(i + 1, j + 1))
          ymin = min (Ycor(i, j), Ycor(i + 1,j), Ycor(i, j + 1), Ycor(i + 1, j + 1))
          ymax = max (Ycor(i, j), Ycor(i + 1,j), Ycor(i, j + 1), Ycor(i + 1, j + 1))

          bad = .false.

          if (Xc(i, j) < xmin .or. Xc(i, j) > xmax) bad = .true.
          if (Yc(i, j) < ymin .or. Yc(i, j) > ymax) bad = .true.

            ! Also catch zeros (your main issue)
          if (abs (Ycor(i, j)) < 1.0e-10) bad = .true.

          if (bad) then
            write (message, *) 'BAD CELL at (i, j)=', i, j
            call Print_message (trim (message))

            write (message, *) ' center =', Xc(i, j), Yc(i, j)
            call Print_message (trim (message))

            write (message, *) ' corners:'
            call Print_message (trim (message))
            write (message, *) '  (i,j)     =', Xcor(i, j), Ycor(i, j)
            call Print_message (trim (message))
            write (message, *) '  (i+1,j)   =', Xcor(i + 1, j), Ycor(i + 1, j)
            call Print_message (trim (message))
            write (message, *) '  (i,j+1)   =', Xcor(i, j + 1), Ycor(i, j + 1)
            call Print_message (trim (message))
            write (message, *) '  (i+1,j+1) =', Xcor(i + 1, j + 1), Ycor(i + 1, j + 1)
            call Print_message (trim (message))

            call Stop_simulation ('ERROR: grid centers not within grid corners')
          end if
        end do
      end do
    end do

  end subroutine Check_grid_cells

end module wrf_nuopc
