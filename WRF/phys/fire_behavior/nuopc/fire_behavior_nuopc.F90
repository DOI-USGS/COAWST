
#include "Fire_behavior_NUOPC_Macros.h"

module fire_behavior_nuopc

#ifdef DM_PARALLEL
  use mpi
#endif

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS    => SetServices

  use state_mod, only : state_fire_t
  use namelist_mod, only : namelist_t
  use initialize_mod, only : Init_fire_state
  use advance_mod, only : Advance_state
  use constants_mod, only : G, XLV, CP, FVIRT, R_D
  use stderrout_mod, only : Stop_simulation, Print_message
  use coupling_mod, only : Calc_fire_wind
  use interp_mod, only: VINTERP_WINDS_FROM_3D_WINDS, VINTERP_WINDS_FROM_10M_WINDS

  implicit none

  private

  public SetVM, SetServices

  type (state_fire_t) :: grid
  type (namelist_t) :: config_flags
  real(ESMF_KIND_R8), pointer     :: ptr_z0(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_t2(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_psfc(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_rainrte(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_rainacc(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_q2(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_lowest_q(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_lowest_t(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_lowest_pres(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_u3d(:,:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_v3d(:,:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_ph(:,:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_hflx_fire(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_evap_fire(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_smoke_fire(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_u10(:,:)
  real(ESMF_KIND_R8), pointer     :: ptr_v10(:,:)

  integer :: clb(2), cub(2), clb3(3), cub3(3)
  logical :: imp_rainrte = .FALSE.
  logical :: imp_rainacc = .FALSE.

  logical, parameter :: DEBUG_ALL = .false.
  integer :: mpi_comm_cfbm

  contains

  subroutine SetServices(model, rc)

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    logical, parameter :: DEBUG_LOCAL = .false.


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering SetServices fire...')

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

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Leaving SetServices fire...')

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Advertise(model, rc)

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)        :: importState, exportState
    integer :: rank, ierr
    logical, parameter :: DEBUG_LOCAL = .false.

    type(ESMF_VM) :: vm


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering Advertise fire...')

    rc = ESMF_SUCCESS

      ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
        exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out

    call ESMF_GridCompGet (model, vm = vm, rc = rc)
    if (ESMF_LogFoundError (rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
        line = __LINE__, file = __FILE__)) return

    call ESMF_VMGet (vm, mpiCommunicator = mpi_comm_cfbm, rc = rc)
    if (ESMF_LogFoundError (rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, &
        line = __LINE__, file = __FILE__)) return

      ! Read namelist
#ifdef DM_PARALLEL
    call Mpi_comm_rank (mpi_comm_cfbm, rank, ierr)
    if (ierr /= MPI_SUCCESS) call Stop_simulation ('ERROR: mpi_comm_rank failed')
#else
    rank = 0
#endif

    if (rank == 0) call config_flags%Initialization (file_name = 'namelist.fire')

#ifdef DM_PARALLEL
    call config_flags%Broadcast_nml (mpi_comm_cfbm)
#endif

#ifdef DM_PARALLEL
    call grid%Set_mpi_comm_cfbm (mpi_comm_cfbm)
#endif

    call Init_fire_state (grid, config_flags)

    ! Import/ Export Variables -----------------------------------------------------

    ! Disabling the following macro, e.g. renaming to WITHIMPORTFIELDS_disable,
    ! will result in a model component that does not advertise any importable
    ! Fields. Use this if you want to drive the model independently.

#define WITHIMPORTFIELDS
#ifdef WITHIMPORTFIELDS
    ! 3D fields

    ! importable field: inst_zonal_wind_levels
    call NUOPC_Advertise(importState, &
      StandardName="inst_zonal_wind_levels", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_merid_wind_levels
    call NUOPC_Advertise(importState, &
      StandardName="inst_merid_wind_levels", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_geop_levels
    call NUOPC_Advertise(importState, &
      StandardName="inst_geop_levels", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !  2D fields

    ! importable field: inst_surface_roughness
    call NUOPC_Advertise(importState, &
      StandardName="inst_surface_roughness", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: accumulated_lwe_thickness_of_precipitation_amount
    call NUOPC_Advertise(importState, &
      StandardName="accumulated_lwe_thickness_of_precipitation_amount", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: mean_prec_rate
    call NUOPC_Advertise(importState, &
      StandardName="mean_prec_rate", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_spec_humid_height2m
    call NUOPC_Advertise(importState, &
      StandardName="inst_spec_humid_height2m", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_pres_height_surface
    call NUOPC_Advertise(importState, &
      StandardName="inst_pres_height_surface", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_temp_height2m
    call NUOPC_Advertise(importState, &
      StandardName="inst_temp_height2m", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_pres_height_lowest_from_phys
    call NUOPC_Advertise(importState, &
      StandardName="inst_pres_height_lowest_from_phys", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_spec_humid_height_lowest_from_phys
    call NUOPC_Advertise(importState, &
      StandardName="inst_spec_humid_height_lowest_from_phys", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_temp_height_lowest_from_phys
    call NUOPC_Advertise(importState, &
      StandardName="inst_temp_height_lowest_from_phys", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_zonal_wind_height10m
    call NUOPC_Advertise(importState, &
      StandardName="inst_zonal_wind_height10m", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable field: inst_merid_wind_height10m
    call NUOPC_Advertise(importState, &
      StandardName="inst_merid_wind_height10m", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


#endif

!#define WITHEXPORTFIELDS_disable
!#ifdef WITHEXPORTFIELDS
    ! exportable field: hflx_fire
    call NUOPC_Advertise(exportState, &
      StandardName="hflx_fire", name="hflx_fire", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: evap_fire
    call NUOPC_Advertise(exportState, &
      StandardName="evap_fire", name="evap_fire", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: smoke_fire
    call NUOPC_Advertise(exportState, &
      StandardName="smoke_fire", name="smoke_fire", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!#endif

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Leaving Advertise fire...')

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Realize(model, rc)

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)        :: importState, exportState
    type(ESMF_Field)        :: field
    type(ESMF_DistGrid)     :: fire_distgrid
    type(ESMF_Grid)         :: fire_grid

    ! working local variables
    integer                        :: lbnd(2),ubnd(2)
    real(ESMF_KIND_COORD), pointer :: coordXcenter(:,:)
    real(ESMF_KIND_COORD), pointer :: coordYcenter(:,:)
    real(ESMF_KIND_COORD), pointer :: coordXcorner(:,:)
    real(ESMF_KIND_COORD), pointer :: coordYcorner(:,:)
    integer                        :: i, j, iglobal, jglobal
    integer, dimension(:, :, :), allocatable :: deBlockList
    logical, parameter :: DEBUG_LOCAL = .false.
    type(ESMF_VM)                   :: vm
    integer :: ierr, petCount, de, localDE
    integer :: sendbuf(4)
    integer, dimension(:), allocatable :: recvbuf
    character (len = 256) :: message


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering Realize fire...')

    rc = ESMF_SUCCESS

      ! Get VM from model
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out

      ! Get number of MPI tasks
    call ESMF_VMGet(vm, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out

    if (DEBUG_LOCAL .or. DEBUG_ALL) then
      write (message, *) 'petCount in fire = ', petCount
      call Print_message (trim (message))
    end if

#ifdef DM_PARALLEL
    allocate (recvbuf(4 * petCount))
    sendbuf = (/ grid%ifps, grid%ifpe, grid%jfps, grid%jfpe /)
    call MPI_Allgather(sendbuf, 4, MPI_INTEGER, &
        recvbuf, 4, MPI_INTEGER, mpi_comm_cfbm, ierr)

    allocate (deBlockList(2, 2, petCount))

      ! deBlockList, 1st index is the coordenate
      ! deBlockList, 2nd index is lower, upper index
      ! deBlockList, 3rd index is patch
    do de = 1, petCount
      deBlockList(1, 1, de) = recvbuf(4 * (de - 1) + 1)
      deBlockList(1, 2, de) = recvbuf(4 * (de - 1) + 2)
      deBlockList(2, 1, de) = recvbuf(4 * (de - 1) + 3)
      deBlockList(2, 2, de) = recvbuf(4 * (de - 1) + 4)
    end do
#else
    allocate (deBlockList(2, 2, 1))
    deBlockList (1, 1, 1) = 1
    deBlockList (1, 2, 1) = grid%nx
    deBlockList (2, 1, 1) = 1
    deBlockList (2, 2, 1) = grid%ny
#endif

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Create distgrid based on the state grid
    fire_distgrid = ESMF_DistGridCreate( &
      minIndex=(/1,1/), maxIndex=(/grid%nx,grid%ny/), &
      deBlockList = deBlockList, &
      rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    fire_grid = ESMF_GridCreate(name='FIRE_BEHAVIOR', &
      distgrid=fire_distgrid, coordSys = ESMF_COORDSYS_SPH_DEG, &
!      coordTypeKind=ESMF_TYPEKIND_COORD, & ?
!      gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
      rc = rc)
    if(ESMF_STDERRORCHECK(rc)) return

    if (allocated(grid%lats) .and. allocated (grid%lons)) then

      ! CENTERS
      localDE = 0
!do localDE = 0, localDECount - 1
      ! Add Center Coordinates to Grid
      call ESMF_GridAddCoord(fire_grid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return
      !
      call ESMF_GridGetCoord(fire_grid, coordDim=1, localDE=localDE, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        computationalLBound=lbnd, computationalUBound=ubnd, &
        farrayPtr=coordXcenter, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_GridGetCoord(fire_grid, coordDim=2, localDE=localDE, &
        staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=coordYcenter, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      do j = lbnd(2), ubnd(2)
        jglobal = grid%jfps + j - 1
        do i = lbnd(1), ubnd(1)
          iglobal = grid%ifps + i - 1
          coordXcenter(i,j) = grid%lons(iglobal, jglobal)
          coordYcenter(i,j) = grid%lats(iglobal, jglobal)
        end do
      end do
!end do

      ! CORNERS
      localDE = 0
      ! Add Corner Coordinates to Grid
      call ESMF_GridAddCoord(fire_grid, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return
      !
!do localDE = 0, localDECount - 1
      call ESMF_GridGetCoord(fire_grid, coordDim=1, localDE=localDE, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        computationalLBound=lbnd, computationalUBound=ubnd, &
        farrayPtr=coordXcorner, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      call ESMF_GridGetCoord(fire_grid, coordDim=2, localDE=localDE, &
        staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=coordYcorner, rc=rc)
      if (ESMF_STDERRORCHECK(rc)) return
      do j = lbnd(2), ubnd(2)
        jglobal = grid%jfps + j - 1
        do i = lbnd(1), ubnd(1)
          iglobal = grid%ifps + i - 1
          coordXcorner(i, j) = grid%lons_c(iglobal, jglobal)
          coordYcorner(i, j) = grid%lats_c(iglobal, jglobal)
        end do
      end do
!end do
    end if

#ifdef WITHIMPORTFIELDS
     !  3D fields

     ! importable field on Grid: inst_zonal_wind_levels
     field = ESMF_FieldCreate(name="inst_zonal_wind_levels", grid=fire_grid, &
       gridToFieldMap=(/1,2/), ungriddedLBound=(/1/), &
       ungriddedUBound=(/grid%kfde - 1/), &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_u3d, &
       computationalLBound=clb3, computationalUBound=cub3, rc=rc)

     ! importable field on Grid: inst_zonal_wind_levels
     field = ESMF_FieldCreate(name="inst_merid_wind_levels", grid=fire_grid, &
       gridToFieldMap=(/1,2/), ungriddedLBound=(/1/), &
       ungriddedUBound=(/grid%kfde - 1/), &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_v3d, rc=rc)

     ! importable field on Grid: inst_geop_levels
     field = ESMF_FieldCreate(name="inst_geop_levels", grid=fire_grid, &
       gridToFieldMap=(/1,2/), ungriddedLBound=(/1/), &
       ungriddedUBound=(/grid%kfde - 1/), &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_ph, rc=rc)

     !  2D fields

     ! importable field on Grid: inst_surface_roughness
     field = ESMF_FieldCreate(name="inst_surface_roughness", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_z0, &
       computationalLBound=clb, computationalUBound=cub, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
!    print *, 'clb(1), cub(1), clb(2), cub(2)', clb(1), cub(1), clb(2), cub(2)

     if (NUOPC_IsConnected(importState, fieldName="mean_prec_rate")) then
       imp_rainrte = .TRUE.
       ! importable field on Grid: mean_prec_rate
       field = ESMF_FieldCreate(name="mean_prec_rate", grid=fire_grid, &
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
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_rainrte, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
     else
       imp_rainrte = .FALSE.
       call ESMF_StateRemove(importState, (/"mean_prec_rate"/), rc=rc)
     endif

     if (NUOPC_IsConnected(importState, fieldName="accumulated_lwe_thickness_of_precipitation_amount")) then
       imp_rainacc = .TRUE.
       ! importable field on Grid: inst_rainfall_amount
       field = ESMF_FieldCreate(name="accumulated_lwe_thickness_of_precipitation_amount", grid=fire_grid, &
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
       call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_rainacc, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
     else
       imp_rainacc = .FALSE.
       call ESMF_StateRemove(importState, (/"accumulated_lwe_thickness_of_precipitation_amount"/), rc=rc)
     endif

     ! importable field on Grid: inst_spec_humid_height2m
     field = ESMF_FieldCreate(name="inst_spec_humid_height2m", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_q2, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

     ! importable field on Grid: inst_pres_height_surface
     field = ESMF_FieldCreate(name="inst_pres_height_surface", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_psfc, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

     ! importable field on Grid: inst_temp_height2m
     field = ESMF_FieldCreate(name="inst_temp_height2m", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_t2, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

     ! importable field on Grid: inst_pres_height_lowest_from_phys
     field = ESMF_FieldCreate(name="inst_pres_height_lowest_from_phys", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_lowest_pres, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

     ! importable field on Grid: inst_spec_humid_height_lowest_from_phys
     field = ESMF_FieldCreate(name="inst_spec_humid_height_lowest_from_phys", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_lowest_q, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

     ! importable field on Grid: inst_temp_height_lowest_from_phys
     field = ESMF_FieldCreate(name="inst_temp_height_lowest_from_phys", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_lowest_t, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

     ! importable field on Grid: inst_zonal_wind_height10m
     field = ESMF_FieldCreate(name="inst_zonal_wind_height10m", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_u10, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

     ! importable field on Grid: inst_merid_wind_height10m
     field = ESMF_FieldCreate(name="inst_merid_wind_height10m", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_v10, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

#endif

! #ifdef WITHEXPORTFIELDS
     ! exportable field on Grid: hflx_fire
     field = ESMF_FieldCreate(name="hflx_fire", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_hflx_fire, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

     ! exportable field on Grid: evap_fire
     field = ESMF_FieldCreate(name="evap_fire", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_evap_fire, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

     ! exportable field on Grid: smoke_fire
     field = ESMF_FieldCreate(name="smoke_fire", grid=fire_grid, &
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
     call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr_smoke_fire, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out
! #endif

    ! Initialize fire behavior exports
    ptr_hflx_fire = 0.
    ptr_evap_fire = 0.
    ptr_smoke_fire = 0.

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Check_fire_grid_cells (fire_grid)

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Leaving Realize fire...')

  end subroutine

  subroutine SetClock(model, rc)

    implicit none

    type(ESMF_GridComp) :: model
    integer, intent(out) :: rc

    type (ESMF_Clock) :: modelClock
    type (ESMF_Time) :: startTime
    type (ESMF_Time) :: stopTime
    type (ESMF_TimeInterval) :: timeStep
    logical, parameter :: DEBUG_LOCAL = .false.


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering SetClock fire...')

    rc = ESMF_SUCCESS

    call ESMF_TimeIntervalSet (timeStep, s_r8 = real (config_flags%dt, kind = ESMF_KIND_R8), rc = rc)
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

    modelClock = ESMF_ClockCreate (name = "Fire Clock", timeStep = timeStep, startTime = startTime, stopTime = stopTime, rc = rc)
    if (ESMF_LogFoundError (rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, line = __LINE__, file = __FILE__)) &
        return

    call ESMF_GridCompSet (model, clock = modelClock, rc = rc)
    if (ESMF_LogFoundError (rcToCheck = rc, msg = ESMF_LOGERR_PASSTHRU, line = __LINE__, file = __FILE__)) &
        return

    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Leaving SetClock fire...')

  end subroutine

  subroutine Advance(model, rc)

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_TimeInterval)     :: timeStep
    real(ESMF_KIND_R8)          :: ts
    type(ESMF_State)            :: importState, exportState
    integer                     :: i, j
    real                        :: q0, rho
    character(len=160)          :: msgString
    real, dimension(:, :, :), allocatable :: atm_u3d, atm_v3d, atm_ph
    real, dimension(:, :), allocatable :: atm_lowest_t, atm_lowest_q, atm_lowest_pres
    real, dimension(:, :), allocatable :: grnhfx_kinematic, grnqfx_kinematic, smoke
    real :: dtratio
    integer :: iims, iime, jims, jime, kims, kime, ioms, iome, joms, jome, iops, iope, jops, jope
    logical, parameter :: DEBUG_LOCAL = .false.


    if (DEBUG_LOCAL .or. DEBUG_ALL) call Print_message ('Entering Advance fire...')

    rc = ESMF_SUCCESS

    ! ratio of fire to atmosphere time step
    dtratio = config_flags%dt / config_flags%interval_atm

    ! query for clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeIntervalGet(timeStep, s_r8=ts, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#ifdef WITHIMPORTFIELDS
    ! Update atmospheric fields
    ! convert cm to m
    grid%fz0(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe) = ptr_z0(clb(1):cub(1),clb(2):cub(2)) * 0.01
    grid%fire_q2(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe) = ptr_q2(clb(1):cub(1),clb(2):cub(2))
    grid%fire_t2(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe) = ptr_t2(clb(1):cub(1),clb(2):cub(2))
    grid%fire_psfc(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe) = ptr_psfc(clb(1):cub(1),clb(2):cub(2))
    if (imp_rainrte) then
      ! convert m s-1 to m and accumulate
      grid%fire_rain(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe) = grid%fire_rain(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe) + &
          ( ptr_rainrte(clb(1):cub(1),clb(2):cub(2)) * ts )
    elseif (imp_rainacc) then
      grid%fire_rain(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe) = ptr_rainacc(clb(1):cub(1),clb(2):cub(2))
    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="missing rainfall import", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    do j = grid%jfps, grid%jfpe
      do i = grid%ifps, grid%ifpe
        grid%fire_q2(i,j) = max (grid%fire_q2(i,j), .001)
        grid%fire_t2(i,j) = max (grid%fire_t2(i,j), 123.4) ! avoid arithmatic error
        grid%fire_psfc(i,j) = max (grid%fire_psfc(i,j), .001)
      end do
    end do

    allocate (atm_u3d(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe, 1:grid%kfde - 1))
    allocate (atm_v3d(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe, 1:grid%kfde - 1))
    allocate (atm_ph(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe, 1:grid%kfde - 1))

    atm_u3d(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe, 1:grid%kfde - 1)  = ptr_u3d(clb3(1):cub3(1),clb3(2):cub3(2),clb3(3):cub3(3))
    atm_v3d(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe, 1:grid%kfde - 1)  = ptr_v3d(clb3(1):cub3(1),clb3(2):cub3(2),clb3(3):cub3(3))
    atm_ph(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe, 1:grid%kfde - 1)   = ptr_ph(clb3(1):cub3(1),clb3(2):cub3(2),clb3(3):cub3(3))

#endif

    select case (config_flags%wind_vinterp_opt)
      case (VINTERP_WINDS_FROM_3D_WINDS)

      iims = grid%ifps
      iime = grid%ifpe
      jims = grid%jfps
      jime = grid%jfpe
      kims = grid%kfds
      kime = grid%kfde - 1

      ioms = grid%ifms
      iome = grid%ifme
      joms = grid%jfms
      jome = grid%jfme

      iops = grid%ifps
      iope = grid%ifpe
      jops = grid%jfps
      jope = grid%jfpe
                                                           ! pass the z0 array without halos
                                                           ! for compatibility with offline sims
      call Calc_fire_wind (atm_u3d, atm_v3d, atm_ph / 9.81, grid%fz0(iims:iime, jims:jime), iims, iime, jims, jime, kims, kime, &
          config_flags%fire_lsm_zcoupling,  config_flags%fire_lsm_zcoupling_ref, config_flags%fire_wind_height, &
          ioms, iome, joms, jome, iops, iope, jops, jope, grid%uf, grid%vf, cap_winds = .true.)

      case (VINTERP_WINDS_FROM_10M_WINDS)
        do j = grid%jfps, grid%jfpe
          do i = grid%ifps, grid%ifpe
            grid%uf(i,j) = grid%fuels%waf(int(grid%nfuel_cat(i,j))) * ptr_u10(i,j) 
            grid%vf(i,j) = grid%fuels%waf(int(grid%nfuel_cat(i,j))) * ptr_v10(i,j)
          end do
        end do
      case default
        call Stop_simulation ('Error: wrong wind_vinterp_opt')

    end select

    if (grid%datetime_now == grid%datetime_start) call grid%Save_state ()

    If_reset_fluxes: if (grid%datetime_now == grid%datetime_next_atm_update) then

      call grid%datetime_now%Print_datetime ()
      call grid%datetime_next_atm_update%Add_seconds (config_flags%interval_atm)

      ptr_hflx_fire = 0.
      ptr_evap_fire = 0.
      ptr_smoke_fire = 0.

    end if If_reset_fluxes

    call Advance_state (grid, config_flags)

    allocate (grnhfx_kinematic(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe))
    allocate (grnqfx_kinematic(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe))
    allocate (smoke(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe))
    allocate (atm_lowest_t(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe))
    allocate (atm_lowest_q(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe))
    allocate (atm_lowest_pres(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe))

    atm_lowest_t(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe)    = ptr_lowest_t(clb(1):cub(1),clb(2):cub(2))
    atm_lowest_q(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe)    = ptr_lowest_q(clb(1):cub(1),clb(2):cub(2))
    atm_lowest_pres(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe) = ptr_lowest_pres(clb(1):cub(1),clb(2):cub(2))

    do j = grid%jfps, grid%jfpe
      do i = grid%ifps, grid%ifpe
        q0   = max(atm_lowest_q(i,j)/(1.-atm_lowest_q(i,j)), 1.e-8)
        rho = atm_lowest_pres(i,j) / (R_D * atm_lowest_t(i,j) * &
            (1.0 + FVIRT * q0))
        if (rho > 0.) then ! avoid unpredictable behavior on the edges
           ! convert [W m-2] to [K m s-1]
          grnhfx_kinematic(i,j)  = grid%fgrnhfx(i,j) / (CP * rho)
           ! convert [W m-2] to [kg kg-1 m s-1]
          grnqfx_kinematic(i,j) = grid%fgrnqfx(i,j) / (XLV * rho)
           ! convert [kg smoke m-2] to [kg smoke kg-1 air]
          smoke(i,j) = grid%emis_smoke(i,j) / ((atm_ph(i,j,2) - atm_ph(i,j,1)) / G * rho)
        end if
      enddo
    enddo

    deallocate (atm_u3d, atm_v3d, atm_ph) !, atm_pres)

    ptr_hflx_fire(clb(1):cub(1),clb(2):cub(2)) = ptr_hflx_fire(clb(1):cub(1),clb(2):cub(2)) \
                 + grnhfx_kinematic(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe) * config_flags%fire_atm_feedback * dtratio
    ptr_evap_fire(clb(1):cub(1),clb(2):cub(2)) = ptr_evap_fire(clb(1):cub(1),clb(2):cub(2)) \
                 + grnqfx_kinematic(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe) * config_flags%fire_atm_feedback * dtratio
    ptr_smoke_fire(clb(1):cub(1),clb(2):cub(2)) = ptr_smoke_fire(clb(1):cub(1),clb(2):cub(2)) \
                 + smoke(grid%ifps:grid%ifpe,grid%jfps:grid%jfpe)

    deallocate(grnhfx_kinematic, grnqfx_kinematic, smoke)

    call grid%Handle_output (config_flags)

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing Fire model from: ", unit=msgString, rc=rc)
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

    if (DEBUG_LOCAL .or. DEBUG_ALL)  call Print_message ('Leaving Advance fire...')

  end subroutine

  subroutine Check_fire_grid_cells (fire_grid)

    implicit none

    type(ESMF_Grid), intent(in) :: fire_grid
    integer :: rc
    integer :: localDE, localDECount
    integer :: lbnd(2), ubnd(2)
    real(ESMF_KIND_COORD), pointer :: Xc(:, :), Yc(:, :) ! centers
    real(ESMF_KIND_COORD), pointer :: Xcor(:, :), Ycor(:, :) ! corners
    integer :: i, j
    real :: xmin, xmax, ymin, ymax
    logical :: bad
    character (len = 256) :: message


    rc = ESMF_SUCCESS

    call ESMF_GridGet(fire_grid, localDECount=localDECount, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    Loop_localde: do localDE = 0, localDECount - 1
      call ESMF_GridGetCoord(fire_grid, coordDim=1, localDE=localDE, &
          staggerloc=ESMF_STAGGERLOC_CENTER, &
          computationalLBound=lbnd, computationalUBound=ubnd, &
          farrayPtr=Xc, rc=rc)

      call ESMF_GridGetCoord(fire_grid, coordDim=2, localDE=localDE, &
          staggerloc=ESMF_STAGGERLOC_CENTER, &
          farrayPtr=Yc, rc=rc)

      call ESMF_GridGetCoord(fire_grid, coordDim=1, localDE=localDE, &
          staggerloc=ESMF_STAGGERLOC_CORNER, &
          farrayPtr=Xcor, rc=rc)

      call ESMF_GridGetCoord(fire_grid, coordDim=2, localDE=localDE, &
          staggerloc=ESMF_STAGGERLOC_CORNER, &
          farrayPtr=Ycor, rc=rc)

      write (message, *) 'Checking localDE =', localDE
      call Print_message (trim (message))

      Loop_jbnd: do j = lbnd(2), ubnd(2) - 1
        Loop_ibnd: do i = lbnd(1), ubnd(1) - 1

          ! Bounding box from 4 corners
          xmin = min (Xcor(i, j), Xcor(i + 1, j), Xcor(i, j + 1), Xcor(i + 1, j + 1))
          xmax = max (Xcor(i, j), Xcor(i + 1, j), Xcor(i, j + 1), Xcor(i + 1, j + 1))
          ymin = min (Ycor(i, j), Ycor(i + 1, j), Ycor(i, j + 1), Ycor(i + 1, j + 1))
          ymax = max (Ycor(i, j), Ycor(i + 1, j), Ycor(i, j + 1), Ycor(i + 1, j + 1))

          bad = .false.

          if (Xc(i, j) < xmin .or. Xc(i, j) > xmax) bad = .true.
          if (Yc(i, j) < ymin .or. Yc(i, j) > ymax) bad = .true.

          ! Also catch zeros (your main issue)
          if (abs (Ycor(i,j)) < 1.0e-10) bad = .true.

          if (bad) then
            write (message, *) 'BAD CELL at (i,j)=', i, j
            call Print_message (trim (message))

            write (message, *) ' center =', Xc(i,j), Yc(i,j)
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

        end do Loop_ibnd
      end do Loop_jbnd
    end do Loop_localde

  end subroutine Check_fire_grid_cells

end module

