#include "cppdefs.h"
#define ESMF_FILENAME "esmf_roms_comp.F"
#define ESMF_CONTEXT __LINE__,ESMF_FILENAME,ESMF_METHOD

      MODULE ESMF_ROMS

#if defined ESMF_COUPLING || defined ESMF_LIB
!
!svn $Id: esmf_roms.F 607 2008-04-23 04:02:27Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS, ESMF Gridded Component:                                  !
!                                                                      !
!  Framework registration of ROMS component  "initialize", "run", and  !
!  "finalize" routines. Notice that only the ROMS_SetServices routine  !
!  needs to be public.                                                 !
!                                                                      !
!=======================================================================
!
!  Required data modules.
!
      USE ESMF_Mod
      USE mod_coupler
!
!  ROMS Component routines.
!
      USE ocean_control_mod, ONLY : ROMS_initialize
      USE ocean_control_mod, ONLY : ROMS_run
      USE ocean_control_mod, ONLY : ROMS_finalize
!
      implicit none
!
!  Public/Private routines.
!
      PUBLIC  :: ROMS_SetServices
      PUBLIC  :: ROMS_SetGridArrays
      PUBLIC  :: ROMS_SetStates
      PUBLIC  :: ROMS_SetClock

      PUBLIC  :: ROMS_GetImportData
      PUBLIC  :: ROMS_PutExportData
      PUBLIC  :: ROMS_PutGridData

      PRIVATE :: ROMS_SetInitialize
      PRIVATE :: ROMS_SetRun
      PRIVATE :: ROMS_SetFinalize
!
!  DistGrid object that describes how arrays are decomposed and
!  distributed over all Decomposition Elements (DEs).
!
      TYPE (ESMF_DistGrid) :: romsDistGrid(Ngrids)
!
!  DELayout object that keeps track of the relationship between
!  its DEs and the resources of the associated VM object.
!
      TYPE (ESMF_DELayout) :: romsDELayout(Ngrids)
!
!  Array descriptor.
!
      TYPE (ESMF_ArraySpec) :: romsArray2dSpec(Ngrids)
!
!  Number of exchange grid arrays.
!
# ifdef MASKING
      integer, parameter :: NromsGrdArrays = 10
# else
      integer, parameter :: NromsGrdArrays = 7
# endif
!
!  Exchange grids information structure.
!
      TYPE (T_MESH) :: romsMesh(NromsGrdArrays)
!
!  Exchange grid array object.
!
      TYPE (ESMF_Array) :: romsGrdArray(NromsGrdArrays,Ngrids)
!
!  Exchange grid array data pointer.
!
      TYPE (T_DATA2D) :: romsGrdData(NromsGrdArrays,Ngrids)
!
!  Export/Import array objects.
!
      TYPE (ESMF_Array), allocatable :: romsImpArray(:,:)
      TYPE (ESMF_Array), allocatable :: romsExpArray(:,:)
!
!  Export/Import fields data pointers.
!
      TYPE (T_DATA2D), allocatable :: romsImpData(:,:)
      TYPE (T_DATA2D), allocatable :: romsExpData(:,:)

      CONTAINS

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_SetServices"

      SUBROUTINE ROMS_SetServices (comp, status)
!
!=======================================================================
!                                                                      !
!  Set routines to be called by ESMF as the "initialize", "run", and   !
!  "finalize" for the ROMS gridded component.                          !
!                                                                      !
!=======================================================================
!
!  Imported variable definitions.
!
      integer, intent(out) :: status

      TYPE (ESMF_GridComp), intent(inout) :: comp
!
!  Local variable definitions.
!
      logical :: CheckError
!
!-----------------------------------------------------------------------
!  Register ROMS gridded components routines.
!-----------------------------------------------------------------------
!
!  Register "initialize" routine.
!
      CALL ESMF_GridCompSetEntryPoint (comp,                            &
     &                                 ESMF_SETINIT,                    &
     &                                 ROMS_SetInitialize,              &
     &                                 ESMF_SINGLEPHASE,                &
     &                                 status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetServices',                &
     &               'registering ROMS initialize routine')) RETURN
!
!  Register "run" routine.
!
      CALL ESMF_GridCompSetEntryPoint (comp,                            &
     &                                 ESMF_SETRUN,                     &
     &                                 ROMS_SetRun,                     &
     &                                 ESMF_SINGLEPHASE,                &
     &                                 status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetServices',                &
     &               'registering ROMS run routine')) RETURN
!
!  Register "finalize" routine.
!
      CALL ESMF_GridCompSetEntryPoint (comp,                            &
     &                                 ESMF_SETFINAL,                   &
     &                                 ROMS_SetFinalize,                &
     &                                 ESMF_SINGLEPHASE,                &
     &                                 status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetServices',                &
     &               'registering ROMS finalize routine')) RETURN
!
!  Set return flag to success.
!
      status=ESMF_SUCCESS

      END SUBROUTINE ROMS_SetServices

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_SetGridArrays"

      SUBROUTINE ROMS_SetGridArrays (comp, status)
!
!=======================================================================
!                                                                      !
!  This routines set-up ROMS grid arrays needed to exchange data       !
!  with other coupled model using ESMF protocols.  This needs to       !
!  be done only once during initialization.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable definitions.
!
      integer, intent(out) :: status

      TYPE (ESMF_GridComp), intent(inout) :: comp
!
!  Local variable declarations.
!
      logical :: CheckError, found

      integer :: MyMPIcomm, Nnodes
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: LBi, UBi, LBj, UBj
      integer :: i, igrid, gtype, j, ng, tile

      integer, dimension(2) :: deCountList
      integer, dimension(2) :: minCorner, maxCorner
      integer, dimension(2) :: CLW, CUW, TLW, TUW

      integer, allocatable :: deBlockList(:,:,:)
      integer, allocatable :: TLWidth(:,:), TUWidth(:,:)
      integer, allocatable :: CLW_r(:,:), CUW_r(:,:)
      integer, allocatable :: CLW_u(:,:), CUW_u(:,:)
      integer, allocatable :: CLW_v(:,:), CUW_v(:,:)

      character (len=40) :: code
      character (len=80) :: name

      TYPE (ESMF_ARRAY) :: GrdArray
!
!-----------------------------------------------------------------------
!  Querry the Virtual Machine (VM) parallel environmemt for the MPI
!  communicator handle and current node rank.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet (comp,                                      &
     &                       vm=VM(Iocean),                             &
     &                       rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetGridArrays',              &
     &               'inquiring VM parallel environment')) RETURN

      CALL ESMF_VMGet (VM(Iocean),                                      &
     &                 localPet=MyRank,                                 &
     &                 petCount=Nnodes,                                 &
     &                 mpiCommunicator=MyMPIcomm,                       &
     &                 rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetGridArrays',              &
     &               'getting mpiCommunicator handle')) RETURN
!
!-----------------------------------------------------------------------
!  Set ROMS exchange grid information structure.
!-----------------------------------------------------------------------
!
      IF (spherical) THEN
        romsMesh( 1) % GridID   = 1
        romsMesh( 1) % GridType = r2dvar
        romsMesh( 1) % code     = 'XposR'
        romsMesh( 1) % Variable = 'lonr'
        romsMesh( 1) % name     = 'longitude, RHO-points'
        romsMesh( 1) % units    = 'degrees_east'

        romsMesh( 2) % GridID   = 2
        romsMesh( 2) % GridType = r2dvar
        romsMesh( 2) % code     = 'YposR'
        romsMesh( 2) % Variable = 'latr'
        romsMesh( 2) % name     = 'latitude, RHO-points'
        romsMesh( 2) % units    = 'degrees_north'

        romsMesh( 3) % GridID   = 3
        romsMesh( 3) % GridType = u2dvar
        romsMesh( 3) % code     = 'XposU'
        romsMesh( 3) % Variable = 'lonu'
        romsMesh( 3) % name     = 'longitude, U-points'
        romsMesh( 3) % units    = 'degrees_east'

        romsMesh( 4) % GridID   = 4
        romsMesh( 4) % GridType = u2dvar
        romsMesh( 4) % code     = 'YposU'
        romsMesh( 4) % Variable = 'latu'
        romsMesh( 4) % name     = 'latitude, U-points'
        romsMesh( 4) % units    = 'degrees_north'

        romsMesh( 5) % GridID   = 5
        romsMesh( 5) % GridType = v2dvar
        romsMesh( 5) % code     = 'XposV'
        romsMesh( 5) % Variable = 'lonv'
        romsMesh( 5) % name     = 'longitude, V-points'
        romsMesh( 5) % units    = 'degrees_east'

        romsMesh( 6) % GridID   = 6
        romsMesh( 6) % GridType = v2dvar
        romsMesh( 6) % code     = 'YposV'
        romsMesh( 6) % Variable = 'latv'
        romsMesh( 6) % name     = 'latitude, V-points'
        romsMesh( 6) % units    = 'degrees_north'
      ELSE
        romsMesh( 1) % GridID   = 1
        romsMesh( 1) % GridType = r2dvar
        romsMesh( 1) % code     = 'XposR'
        romsMesh( 1) % Variable = 'xr'
        romsMesh( 1) % name     = 'Cartesian X-positions, RHO-points'
        romsMesh( 1) % units    = 'meter'

        romsMesh( 2) % GridID   = 2
        romsMesh( 2) % GridType = r2dvar
        romsMesh( 2) % code     = 'YposR'
        romsMesh( 2) % Variable = 'yr'
        romsMesh( 2) % name     = 'Cartesian Y-positions, RHO-points'
        romsMesh( 2) % units    = 'meter'

        romsMesh( 3) % GridID   = 3
        romsMesh( 3) % GridType = u2dvar
        romsMesh( 3) % code     = 'XposU'
        romsMesh( 3) % Variable = 'xu'
        romsMesh( 3) % name     = 'Cartesian X-positions, U-points'
        romsMesh( 3) % units    = 'meter'

        romsMesh( 4) % GridID   = 4
        romsMesh( 4) % GridType = u2dvar
        romsMesh( 4) % code     = 'YposU'
        romsMesh( 4) % Variable = 'yu'
        romsMesh( 4) % name     = 'Cartesian Y-positions, U-points'
        romsMesh( 4) % units    = 'meter'

        romsMesh( 5) % GridID   = 5
        romsMesh( 5) % GridType = v2dvar
        romsMesh( 5) % code     = 'XposV'
        romsMesh( 5) % Variable = 'xv'
        romsMesh( 5) % name     = 'Cartesian X-positions, V-points'
        romsMesh( 5) % units    = 'meter'

        romsMesh( 6) % GridID   = 6
        romsMesh( 6) % GridType = v2dvar
        romsMesh( 6) % code     = 'YposV'
        romsMesh( 6) % Variable = 'yv'
        romsMesh( 6) % name     = 'Cartesian Y-positions, V-points'
        romsMesh( 6) % units    = 'meter'
      END IF
      romsMesh( 7) % GridID   = 7
      romsMesh( 7) % GridType = r2dvar
      romsMesh( 7) % code     = 'angleR'
      romsMesh( 7) % Variable = 'angler'
      romsMesh( 7) % name     = 'curvilinear rotation angle, RHO-points'
      romsMesh( 7) % units    = 'radians'
# ifdef MASKING
      romsMesh( 8) % GridID   = 8
      romsMesh( 8) % GridType = r2dvar
      romsMesh( 8) % code     = 'maskR'
      romsMesh( 8) % Variable = 'rmask'
      romsMesh( 8) % name     = 'land/sea mask, RHO-points'
      romsMesh( 8) % units    = 'nondimensional'

      romsMesh( 9) % GridID   = 9
      romsMesh( 9) % GridType = u2dvar
      romsMesh( 9) % code     = 'maskU'
      romsMesh( 9) % Variable = 'umask'
      romsMesh( 9) % name     = 'land/sea mask, U-points'
      romsMesh( 9) % units    = 'nondimensional'

      romsMesh(10) % GridID   = 10
      romsMesh(10) % GridType = v2dvar
      romsMesh(10) % code     = 'maskV'
      romsMesh(10) % Variable = 'vmask'
      romsMesh(10) % name     = 'land/sea mask, V-points'
      romsMesh(10) % units    = 'nondimensional'
# endif
!
!-----------------------------------------------------------------------
!  Set ROMS domain decomposition variables
!-----------------------------------------------------------------------
!
!  Loop over number of nested/composed grids. The strategy here is
!  to set all exchange arrays for all grids in case of composed or
!  mosaic applications.
!
      DO ng=1,Ngrids
!
!  Set tiles lower and upper bounds for each decomposition element.
!  In ROMS, the "exclusive region" for each decomposition element or
!  horizontal tile ranges is bounded by (Istr:Iend, Jstr:Jend). Each
!  tiled array is dimensioned as (LBi:UBi, LBj:UBj) which includes
!  halo regions (usually 2 ghost points) and padding when appropriate
!  (total/memory region). All ROMS arrays are horizontally dimensioned
!  with the same bounds regardless if they are variables located at
!  RHO-, PSI-, U-, or V-points. There is no halos at the boundary edges.
!  The physical boundary is a U-points (east/west edged) and V-points
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
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          deBlockList(1,1,tile+1)=BOUNDS(ng)%Istr(tile)
          deBlockList(1,2,tile+1)=BOUNDS(ng)%Iend(tile)
          deBlockList(2,1,tile+1)=BOUNDS(ng)%Jstr(tile)
          deBlockList(2,2,tile+1)=BOUNDS(ng)%Jend(tile)
        END DO
!
!  Set computational and total regions lower and upper corner widths
!  for each staggered C-grid ROMS arrays.  The total region width is
!  the same for variables at RHO-, U-, and V-points since the exclusive
!  region bounded between Istr:Iend and Jstr:Jend. The starting indices
!  IstrU and JstrV are not used.  This implies that the computational
!  widths are set to zero at those edges for U- and V-type fields.
!
        IF (.not.allocated(TLWidth)) THEN
          allocate ( TLWidth(2,0:NtileI(ng)*NtileJ(ng)-1) )
          allocate ( TUWidth(2,0:NtileI(ng)*NtileJ(ng)-1) )
          allocate ( CLW_r(2,0:NtileI(ng)*NtileJ(ng)-1) )
          allocate ( CUW_r(2,0:NtileI(ng)*NtileJ(ng)-1) )
          allocate ( CLW_u(2,0:NtileI(ng)*NtileJ(ng)-1) )
          allocate ( CUW_u(2,0:NtileI(ng)*NtileJ(ng)-1) )
          allocate ( CLW_v(2,0:NtileI(ng)*NtileJ(ng)-1) )
          allocate ( CUW_v(2,0:NtileI(ng)*NtileJ(ng)-1) )
        END IF
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          TLWidth(1,tile)=BOUNDS(ng)%Istr(tile)-BOUNDS(ng)%LBi(tile)
          TLWidth(2,tile)=BOUNDS(ng)%Jstr(tile)-BOUNDS(ng)%LBj(tile)
          TUWidth(1,tile)=BOUNDS(ng)%UBi(tile)-BOUNDS(ng)%Iend(tile)
          TUWidth(2,tile)=BOUNDS(ng)%UBj(tile)-BOUNDS(ng)%Jend(tile)

          CLW_r(1,tile)=BOUNDS(ng)%Istr(tile)-BOUNDS(ng)%IstrR(tile)
          CLW_r(2,tile)=BOUNDS(ng)%Jstr(tile)-BOUNDS(ng)%JstrR(tile)
          CUW_r(1,tile)=BOUNDS(ng)%IendR(tile)-BOUNDS(ng)%Iend(tile)
          CUW_r(2,tile)=BOUNDS(ng)%JendR(tile)-BOUNDS(ng)%Jend(tile)

          CLW_u(1,tile)=0
          CLW_u(2,tile)=CLW_r(2,tile)
          CUW_u(1,tile)=CUW_r(1,tile)
          CUW_u(2,tile)=CUW_r(2,tile)

          CLW_v(1,tile)=CLW_r(1,tile)
          CLW_v(2,tile)=0
          CUW_v(1,tile)=CUW_r(1,tile)
          CUW_v(2,tile)=CUW_r(2,tile)
        END DO
        TLW=(/TLWidth(1,MyRank), TLWidth(2,MyRank)/)
        TUW=(/TUWidth(1,MyRank), TUWidth(2,MyRank)/)
!
        IF (Master.eq.0) THEN
          PRINT *, ' '
          PRINT *, 'Horizontal decomposition indices per tile:'
          PRINT *, ' '
          PRINT 10, 'Istr   = ',(BOUNDS(ng)%Istr(tile),tile=0,Nnodes-1)
          PRINT 10, 'IstrU  = ',(BOUNDS(ng)%IstrU(tile),tile=0,Nnodes-1)
          PRINT 10, 'Iend   = ',(BOUNDS(ng)%Iend(tile),tile=0,Nnodes-1)
          PRINT 10, 'Jstr   = ',(BOUNDS(ng)%Jstr(tile),tile=0,Nnodes-1)
          PRINT 10, 'JstrV  = ',(BOUNDS(ng)%JstrV(tile),tile=0,Nnodes-1)
          PRINT 10, 'Jend   = ',(BOUNDS(ng)%Jend(tile),tile=0,Nnodes-1)
          PRINT *, ' '
          PRINT 10, 'LBi    = ',(BOUNDS(ng)%LBi(tile),tile=0,Nnodes-1)
          PRINT 10, 'UBi    = ',(BOUNDS(ng)%UBi(tile),tile=0,Nnodes-1)
          PRINT 10, 'LBj    = ',(BOUNDS(ng)%LBj(tile),tile=0,Nnodes-1)
          PRINT 10, 'UBj    = ',(BOUNDS(ng)%UBj(tile),tile=0,Nnodes-1)
          PRINT *, ' '
          PRINT 10, 'TLWi   = ',(TLWidth(1,tile),tile=0,Nnodes-1)
          PRINT 10, 'TLWj   = ',(TLWidth(2,tile),tile=0,Nnodes-1)
          PRINT 10, 'TUWi   = ',(TUWidth(1,tile),tile=0,Nnodes-1)
          PRINT 10, 'TUWj   = ',(TUWidth(2,tile),tile=0,Nnodes-1)
          PRINT *, ' '
          PRINT 10, 'CLWi_r = ',(CLW_r(1,tile),tile=0,Nnodes-1)
          PRINT 10, 'CLWj_r = ',(CLW_r(2,tile),tile=0,Nnodes-1)
          PRINT 10, 'CUWi_r = ',(CUW_r(1,tile),tile=0,Nnodes-1)
          PRINT 10, 'CUWj_r = ',(CUW_r(2,tile),tile=0,Nnodes-1)
          PRINT *, ' '
          PRINT 10, 'CLWi_u = ',(CLW_u(1,tile),tile=0,Nnodes-1)
          PRINT 10, 'CLWj_u = ',(CLW_u(2,tile),tile=0,Nnodes-1)
          PRINT 10, 'CUWi_u = ',(CUW_u(1,tile),tile=0,Nnodes-1)
          PRINT 10, 'CUWj_u = ',(CUW_u(2,tile),tile=0,Nnodes-1)
          PRINT *, ' '
          PRINT 10, 'CLWi_v = ',(CLW_u(1,tile),tile=0,Nnodes-1)
          PRINT 10, 'CLWj_v = ',(CLW_u(2,tile),tile=0,Nnodes-1)
          PRINT 10, 'CUWi_v = ',(CUW_u(1,tile),tile=0,Nnodes-1)
          PRINT 10, 'CUWj_v = ',(CUW_u(2,tile),tile=0,Nnodes-1)
 10       FORMAT (1x,a,64i5)
        END IF
!
!-----------------------------------------------------------------------
!  Set ESMF Layout and distribution objects.
!-----------------------------------------------------------------------
!
!  Create ESMF DELayout which describes how the data elements (DE) are
!  laid out over the VM computational resources.
!
        deCountList=(/NtileI(ng), NtileJ(ng)/)
        romsDELayout(ng)=ESMF_DELayoutCreate(VM(Iocean),                &
     &                                       deCountList=deCountList,   &
     &                                       petList=pets(Iocean)%val,  &
     &                                       rc=status)
        IF (CheckError(status, 'ESMF', 'ROMS_SetGridArrays',            &
     &                 'creating DELayout')) RETURN
!
!  Create ESMF DistGrid based on ROMS domain decomposition.
!
        minCorner=(/1, 1/)
        maxCorner=(/Lm(ng), Mm(ng)/)
        romsDistGrid(ng)=ESMF_DistGridCreate(minCorner,                 &
     &                                       maxCorner,                 &
     &                                       deBlockList=deBlockList,   &
     &                                       deLayout=romsDELayout(ng), &
     &                                       vm=VM(Iocean),             &
     &                                       rc=status)
        IF (CheckError(status, 'ESMF', 'ROMS_SetGridArrays',            &
     &                 'creating DistGrid')) RETURN
        deallocate (deBlockList)
!
!-----------------------------------------------------------------------
!  Set array descriptor.
!-----------------------------------------------------------------------
!
!  Currently, all ROMS arrays used during coupling are floating-point,
!  double-precision and two-dimensional.
!
        CALL ESMF_ArraySpecSet (romsArray2dSpec(ng),                    &
     &                          rank=2,                                 &
     &                          typekind=ESMF_TYPEKIND_R8,              &
     &                          rc=status)
        IF (CheckError(status, 'ESMF', 'ROMS_SetGridArrays',            &
     &                 'setting 2D array descriptor')) RETURN
!
!-----------------------------------------------------------------------
!  Create ROMS exchange arrays. Notice that after initialization all
!  ROMS grid arrays are allocated and initialized.
!-----------------------------------------------------------------------
!
        DO igrid=1,NromsGrdArrays
          code=romsMesh(igrid)%code
          name=romsMesh(igrid)%name
          gtype=romsMesh(igrid)%GridType
!
!  Set computational region widths.
!
          IF (gtype.eq.u2dvar) THEN
            CLW=(/CLW_u(1,MyRank), CLW_u(2,MyRank)/)
            CUW=(/CUW_u(1,MyRank), CUW_u(2,MyRank)/)
          ELSE IF (gtype.eq.v2dvar) THEN
            CLW=(/CLW_v(1,MyRank), CLW_v(2,MyRank)/)
            CUW=(/CUW_v(1,MyRank), CUW_v(2,MyRank)/)
          ELSE
            CLW=(/CLW_r(1,MyRank), CLW_r(2,MyRank)/)
            CUW=(/CUW_r(1,MyRank), CUW_r(2,MyRank)/)
          END IF
!
!  Create array.
!
          GrdArray=ESMF_ArrayCreate(romsArray2dSpec(ng),                &
     &                              distgrid=romsDistGrid(ng),          &
     &                              computationalLWidth=CLW,            &
     &                              computationalUWidth=CUW,            &
     &                              totalLWidth=TLW,                    &
     &                              totalUWidth=TUW,                    &
     &                              indexflag=ESMF_INDEX_GLOBAL,        &
     &                              rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetGridArrays',          &
     &                   'creating grid array '//TRIM(name))) RETURN
          romsGrdArray(igrid,ng)=GrdArray
!
!  Get data pointer from array.
!
          CALL ESMF_ArrayGet (romsGrdArray(igrid,ng),                   &
     &                        romsGrdData(igrid,ng)%field,              &
     &                        rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetGridArrays',          &
     &                   'getting array pointer '//TRIM(name))) RETURN
!
!  Set array name.
!
          CALL ESMF_ArraySet (romsGrdArray(igrid,ng),                   &
     &                        TRIM(name),                               &
     &                        rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetGridArrays',          &
     &                   'setting grid array '//TRIM(name))) RETURN
!
!  Add array to state.
!
          CALL ESMF_StateAddArray (StateExport(Iocean),                 &
     &                             romsGrdArray(igrid,ng),              &
     &                             rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetGridArrays',          &
     &                   'adding grid array '//TRIM(name))) RETURN
!
!  Initialize exchange grid data to zero to avoid infinities or NaNs.
!
          romsGrdData(igrid,ng)%field=0.0_r8
        END DO
!
!  Deallocate grid widths.
!
        deallocate (TLWidth)
        deallocate (TUWidth)
        deallocate (CLW_r)
        deallocate (CUW_r)
        deallocate (CLW_u)
        deallocate (CUW_u)
        deallocate (CLW_v)
        deallocate (CUW_v)
      END DO
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS

      END SUBROUTINE ROMS_SetGridArrays

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_SetStates"

      SUBROUTINE ROMS_SetStates (ImportState, ExportState,              &
     &                           MyRank, status)
!
!=======================================================================
!                                                                      !
!  This routines set-up ROMS grid arrays needed to exchange data       !
!  with other coupled model using ESMF protocols.  This needs to       !
!  be done only once during initialization.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
!  Imported variable definitions.
!
      integer, intent(in)  :: MyRank
      integer, intent(out) :: status

      TYPE (ESMF_State), intent(inout) :: ImportState
      TYPE (ESMF_State), intent(inout) :: ExportState
!
!  Local variable declarations.
!
      logical :: CheckError

      integer :: gtype, i, id, ng

      integer, dimension(2) :: CLW, CUW, TLW, TUW
      integer, dimension(2) :: CLW_r, CUW_r
      integer, dimension(2) :: CLW_u, CUW_u
      integer, dimension(2) :: CLW_v, CUW_v

      character (len=40) :: Aname

      TYPE (ESMF_ARRAY) :: ImpArray, ExpArray
!
!-----------------------------------------------------------------------
!  Allocate.
!-----------------------------------------------------------------------
!
!  Export/Import array objects.
!
      IF (.not.allocated(romsExpArray)) THEN
        allocate ( romsExpArray(Nexport(Iocean),Ngrids) )
      END IF
      IF (.not.allocated(romsImpArray)) THEN
        allocate ( romsImpArray(Nimport(Iocean),Ngrids) )
      END IF
!
!  Export/Import fields data pointers.
!
      IF (.not.allocated(romsExpData)) THEN
        allocate ( romsExpData(Nexport(Iocean),Ngrids) )
      END IF
      IF (.not.allocated(romsImpData)) THEN
        allocate ( romsImpData(Nimport(Iocean),Ngrids) )
      END IF
!
!-----------------------------------------------------------------------
!  Set ROMS domain decomposition variables
!-----------------------------------------------------------------------
!
!  Loop over number of nested/composed grids. The strategy here is
!  to set all import and export states for all grids in case of
!  composed or mosaic applications.
!
      DO ng=1,Ngrids
!
!  Set computational and total regions lower and upper corner widths
!  for each staggered C-grid ROMS arrays.
!
        TLW(1)=BOUNDS(ng)%Istr(MyRank)-BOUNDS(ng)%LBi(MyRank)
        TLW(2)=BOUNDS(ng)%Jstr(MyRank)-BOUNDS(ng)%LBj(MyRank)
        TUW(1)=BOUNDS(ng)%UBi(MyRank)-BOUNDS(ng)%Iend(MyRank)
        TUW(2)=BOUNDS(ng)%UBj(MyRank)-BOUNDS(ng)%Jend(MyRank)
!
        CLW_r(1)=BOUNDS(ng)%Istr(MyRank)-BOUNDS(ng)%IstrR(MyRank)
        CLW_r(2)=BOUNDS(ng)%Jstr(MyRank)-BOUNDS(ng)%JstrR(MyRank)
        CUW_r(1)=BOUNDS(ng)%IendR(MyRank)-BOUNDS(ng)%Iend(MyRank)
        CUW_r(2)=BOUNDS(ng)%JendR(MyRank)-BOUNDS(ng)%Jend(MyRank)
!
        CLW_u(1)=0
        CLW_u(2)=CLW_r(2)
        CUW_u(1)=CUW_r(1)
        CUW_u(2)=CUW_r(2)
!
        CLW_v(1)=CLW_r(1)
        CLW_v(2)=0
        CUW_v(1)=CUW_r(1)
        CUW_v(2)=CUW_r(2)
!
!-----------------------------------------------------------------------
!  Create export state arrays.
!-----------------------------------------------------------------------
!
        DO i=1,Nexport(Iocean)
          id=ExportID(Iocean)%val(i)
          Aname=Fields(id)%name
!
!  Set computational region width.
!
          gtype=Fields(id)%GridType
          IF ((gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
            CLW=(/CLW_u(1), CLW_u(2)/)
            CUW=(/CUW_u(1), CUW_u(2)/)
          ELSE IF ((gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
            CLW=(/CLW_v(1), CLW_v(2)/)
            CUW=(/CUW_v(1), CUW_v(2)/)
          ELSE
            CLW=(/CLW_r(1), CLW_r(2)/)
            CUW=(/CUW_r(1), CUW_r(2)/)
          END IF
!
!  Create array.
!
          ExpArray=ESMF_ArrayCreate(romsArray2dSpec(ng),                  &
     &                              distgrid=romsDistGrid(ng),          &
     &                              computationalLWidth=CLW,            &
     &                              computationalUWidth=CUW,            &
     &                              totalLWidth=TLW,                    &
     &                              totalUWidth=TUW,                    &
     &                              indexflag=ESMF_INDEX_GLOBAL,        &
     &                              rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetStates',              &
     &                   'creating export array '//TRIM(Aname))) RETURN
          romsExpArray(i,ng)=ExpArray
!
!  Get data pointer from array.
!
          CALL ESMF_ArrayGet (romsExpArray(i,ng),                       &
     &                        romsExpData(i,ng)%field,                  &
     &                        rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetStates',              &
     &                   'getting array pointer '//TRIM(Aname))) RETURN
!
!  Set array name.
!
          CALL ESMF_ArraySet (romsExpArray(i,ng),                       &
     &                        TRIM(Fields(id)%variable),                &
     &                        rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetStates',              &
     &                   'setting export array '//TRIM(Aname))) RETURN
!
!  Add array to state.
!
          CALL ESMF_StateAddArray (ExportState,                         &
     &                             romsExpArray(i,ng),                  &
     &                             rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetStates',              &
     &                   'adding export array '//TRIM(Aname))) RETURN
!
!  Initialize import field to zero to avoid infinities or NaNs.
!
          romsExpData(i,ng)%field=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Create import state arrays.
!-----------------------------------------------------------------------
!
        DO i=1,Nimport(Iocean)
          id=ImportID(Iocean)%val(i)
          Aname=Fields(id)%name
!
!  Set computational region width.
!
          gtype=Fields(id)%GridType
          IF ((gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
            CLW=(/CLW_u(1), CLW_u(2)/)
            CUW=(/CUW_u(1), CUW_u(2)/)
          ELSE IF ((gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
            CLW=(/CLW_v(1), CLW_v(2)/)
            CUW=(/CUW_v(1), CUW_v(2)/)
          ELSE
            CLW=(/CLW_r(1), CLW_r(2)/)
            CUW=(/CUW_r(1), CUW_r(2)/)
          END IF
!
!  Create array.
!
          ImpArray=ESMF_ArrayCreate(romsArray2dSpec(ng),                &
     &                              distgrid=romsDistGrid(ng),          &
     &                              computationalLWidth=CLW,            &
     &                              computationalUWidth=CUW,            &
     &                              totalLWidth=TLW,                    &
     &                              totalUWidth=TUW,                    &
     &                              indexflag=ESMF_INDEX_GLOBAL,        &
     &                              rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetStates',              &
     &                   'creating import array '//TRIM(Aname))) RETURN
          romsImpArray(i,ng)=ImpArray
!
!  Get data pointer from array.
!
          CALL ESMF_ArrayGet (romsImpArray(i,ng),                       &
     &                        romsImpData(i,ng)%field,                  &
     &                        rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetStates',              &
     &                   'getting array pointer '//TRIM(Aname))) RETURN
!
!  Set array name.
!
          CALL ESMF_ArraySet (romsImpArray(i,ng),                       &
     &                        TRIM(Fields(id)%variable),                &
     &                        rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetStates',              &
     &                   'setting import array '//TRIM(Aname))) RETURN
!
!  Add array to state.
!
          CALL ESMF_StateAddArray (ImportState,                         &
     &                             romsImpArray(i,ng),                  &
     &                             rc=status)
          IF (CheckError(status, 'ESMF', 'ROMS_SetStates',              &
     &                   'adding import array '//TRIM(Aname))) RETURN
!
!  Initialize import field to zero to avoid infinities or NaNs.
!
          romsImpData(i,ng)%field=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS

      END SUBROUTINE ROMS_SetStates

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_SetClock"

      SUBROUTINE ROMS_SetClock (clock, status)
!
!=======================================================================
!                                                                      !
!  Set routine sets ROMS date calendar, start and stop time, and time  !
!  interval.  At initilization,  the variable "tdays" is the initial   !
!  time meassured in fractional days since the reference time.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable definitions.
!
      integer, intent(out) :: status

      TYPE(ESMF_Clock), intent(inout) :: clock
!
!  Local variable declarations.
!
      logical :: CheckError

      integer :: ref_year,   start_year,   stop_year
      integer :: ref_month,  start_month,  stop_month
      integer :: ref_day,    start_day,    stop_day
      integer :: ref_hour,   start_hour,   stop_hour
      integer :: ref_minute, start_minute, stop_minute
      integer :: ref_second, start_second, stop_second
      integer :: MyTimeStep, ng

      real(r8) :: MyStartTime, MyStopTime
      real(r8) :: hour, minute, yday

      character (len=80) :: name
!
!-----------------------------------------------------------------------
!  Set coupled models time clocks in seconds. For now, use ROMS time
!  managing variables to set coupled models time clocks.
!-----------------------------------------------------------------------
!
!  Set time reference: model time is meassured as seconds since
!  reference time.
!
      IF (INT(time_ref).eq.-2) THEN
        ref_year=1968
        ref_month=5
        ref_day=23
        ref_hour=0
        ref_minute=0
        ref_second=0
        name='Modified Julian day number, Gregorian Calendar'
        TimeCalendar(Iocean)=ESMF_CalendarCreate(TRIM(name),            &
     &                                           ESMF_CAL_GREGORIAN,    &
     &                                           status)
      ELSE IF (INT(time_ref).eq.-1) THEN
        ref_year=1
        ref_month=1
        ref_day=1
        ref_hour=0
        ref_minute=0
        ref_second=0
        name='360-day, 30 days per month'
        TimeCalendar(Iocean)=ESMF_CalendarCreate(TRIM(name),            &
     &                                           ESMF_CAL_360DAY,       &
     &                                           status)
      ELSE IF (INT(time_ref).eq.0) THEN
        ref_year=1
        ref_month=1
        ref_day=1
        ref_hour=0
        ref_minute=0
        ref_second=0
        name='Julian Calendar, leap year if divisible by 4'
        TimeCalendar(Iocean)=ESMF_CalendarCreate(TRIM(name),            &
     &                                           ESMF_CAL_JULIAN,       &
     &                                           status)
      ELSE IF (time_ref.gt.0.0_r8) THEN
        ref_year=INT(r_date(2))
        ref_month=INT(r_date(4))
        ref_day=INT(r_date(5))
        ref_hour=INT(r_date(6))
        ref_minute=INT(r_date(7))
        ref_second=INT(r_date(8))
        name='Gregorian Calendar'
        TimeCalendar(Iocean)=ESMF_CalendarCreate(TRIM(name),            &
     &                                           ESMF_CAL_GREGORIAN,    &
     &                                           status)
      END IF
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'creating '//TRIM(name))) RETURN
!
!  Set Reference time.
!
      CALL ESMF_TimeSet (ReferenceTime(Iocean),                         &
     &                   yy=ref_year,                                   &
     &                   mm=ref_month,                                  &
     &                   dd=ref_day,                                    &
     &                   h=ref_hour,                                    &
     &                   m=ref_minute,                                  &
     &                   s=ref_second,                                  &
     &                   calendar=TimeCalendar(Iocean),                 &
     &                   rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'setting reference time')) RETURN
!
!  Set start time, use the minimum value of all nested grids.
!
      MyStartTime=MINVAL(tdays)
      CALL caldate (r_date, MyStartTime, start_year, yday, start_month, &
     &              start_day, hour)
      minute=(hour-AINT(hour))*60.0_r8
      start_hour=INT(hour)
      start_minute=INT(minute)
      start_second=INT((minute-AINT(minute))*60.0_r8)
!
      CALL ESMF_TimeSet (StartTime(Iocean),                             &
     &                   yy=start_year,                                 &
     &                   mm=start_month,                                &
     &                   dd=start_day,                                  &
     &                   h=start_hour,                                  &
     &                   m=start_minute,                                &
     &                   s=start_second,                                &
     &                   calendar=TimeCalendar(Iocean),                 &
     &                   rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'setting start time')) RETURN
!
!  Set stop time, use the maximum value of all nested grids.
!
      MyStopTime=0.0_r8
      DO ng=1,Ngrids
        MyStopTime=MAX(MyStopTime,                                      &
     &                 tdays(ng)+(REAL(ntimes(ng),r8)*dt(ng))*sec2day)
      END DO
      CALL caldate (r_date, MyStopTime, stop_year, yday, stop_month,    &
     &              stop_day, hour)
      minute=(hour-AINT(hour))*60.0_r8
      stop_hour=INT(hour)
      stop_minute=INT(minute)
      stop_second=INT((minute-AINT(minute))*60.0_r8)
!
      CALL ESMF_TimeSet (StopTime(Iocean),                              &
     &                   yy=stop_year,                                  &
     &                   mm=stop_month,                                 &
     &                   dd=stop_day,                                   &
     &                   h=stop_hour,                                   &
     &                   m=stop_minute,                                 &
     &                   s=stop_second,                                 &
     &                   calendar=TimeCalendar(Iocean),                 &
     &                   rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'setting stop time')) RETURN
!
!  Set time interval, use the maximum value of all nested grids.
!  Currently, only integer type variables are implemented in ESMF.
!
      MyTimeStep=INT(MAXVAL(dt))
      CALL ESMF_TimeIntervalSet (TimeStep(Iocean),                      &
     &                           s=MyTimeStep,                          &
     &                           rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'setting time interval')) RETURN
!
!  Create time clock.
!
      name='ROMS Time Clock'
      TimeClock(Iocean)=ESMF_ClockCreate(TRIM(name),                    &
     &                                   refTime=ReferenceTime(Iocean), &
     &                                   timeStep=TimeStep(Iocean),     &
     &                                   startTime=StartTime(Iocean),   &
     &                                   stopTime=StopTime(Iocean),     &
     &                                   rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'creating time clock')) RETURN
!
!  Validate time clock.
!
      CALL ESMF_ClockValidate (TimeClock(Iocean),                       &
     &                         rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'validating time clock')) RETURN
!
!  Get ROMS internal clock current time.
!
      CALL ESMF_ClockGet (TimeClock(Iocean),                            &
     &                    currTime=CurrTime(Iocean),                    &
     &                    rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'getting ROMS current time object')) RETURN
      CALL ESMF_TimeGet (CurrTime(Iocean),                              &
     &                   yy=CurrentTime(Iocean)%year,                   &
     &                   mm=CurrentTime(Iocean)%month,                  &
     &                   dd=CurrentTime(Iocean)%day,                    &
     &                   h=CurrentTime(Iocean)%hour,                    &
     &                   m=CurrentTime(Iocean)%minute,                  &
     &                   s=CurrentTime(Iocean)%second,                  &
     &                   timeZone=CurrentTime(Iocean)%TimeZone,         &
     &                   timeStringISOFrac=CurrentTime(Iocean)%string,  &
     &                   dayOfYear=CurrentTime(Iocean)%YearDay,         &
     &                   rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetClock',                   &
     &               'getting ROMS current time')) RETURN
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS

      END SUBROUTINE ROMS_SetClock

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_GetImportData"

      SUBROUTINE ROMS_GetImportData (MyRank, Kimp, Nimp, status)
!
!=======================================================================
!                                                                      !
!  This routines process import data from "ROMS_ImpData" pointer       !
!  structure and interpolates/loads to respective state arrays.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     MyRank     Local PET                                             !
!     Kimp       Import time index for 2D ROMS state arrays,           !
!                  [1:Ngrids].                                         !
!     Nimp       Import time index for 3D ROMS state arrays,           !
!                  [1:Ngrids].                                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     status     Success/failure error flag.                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ocean
      USE mod_forces
      USE mod_scalars
!
      USE ROMS_import_mod, ONLY : ROMS_import2d
!
!  Imported variable definitions.
!
      integer, intent(in)  :: MyRank
      integer, intent(in)  :: Kimp(1:Ngrids)
      integer, intent(in)  :: Nimp(1:Ngrids)
      integer, intent(out) :: status
!
!  Local variable declarations.
!
      logical :: CheckError

      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: LBi, UBi, LBj, UBj
      integer :: gtype, i, id, ifield, j, ng

      real(r8) :: add_offset, scale
      character (len=40) :: code
!
!-----------------------------------------------------------------------
!  Get import fields.
!-----------------------------------------------------------------------
!
!  Loop over all nested, composed, or mosaic grids.
!
      DO ng=1,Ngrids
        Istr =BOUNDS(ng)%Istr (MyRank)
        Iend =BOUNDS(ng)%Iend (MyRank)
        Jstr =BOUNDS(ng)%Jstr (MyRank)
        Jend =BOUNDS(ng)%Jend (MyRank)
        IstrR=BOUNDS(ng)%IstrR(MyRank)
        IendR=BOUNDS(ng)%IendR(MyRank)
        IstrU=BOUNDS(ng)%IstrU(MyRank)
        JstrR=BOUNDS(ng)%JstrR(MyRank)
        JendR=BOUNDS(ng)%JendR(MyRank)
        JstrV=BOUNDS(ng)%JstrV(MyRank)
        LBi=BOUNDS(ng)%LBi(MyRank)
        UBi=BOUNDS(ng)%UBi(MyRank)
        LBj=BOUNDS(ng)%LBj(MyRank)
        UBj=BOUNDS(ng)%UBj(MyRank)

        DO ifield=1,Nimport(Iocean)
          id=ImportID(Iocean)%val(ifield)
          code=Fields(id)%code
          gtype=Fields(id)%GridType
          scale=Fields(id)%scale
          add_offset=Fields(id)%AddOffset

          SELECT CASE (TRIM(ADJUSTL(code)))
# ifdef CLOUDS
            CASE ('cloud')                    ! cloud fraction
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%cloud,                     &
     &                            status)
# endif
# ifdef SHORTWAVE
            CASE ('SWrad')                    ! shortwave radiation
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%swrad,                     &
     &                            status)
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
            CASE ('Pair')                     ! surface air pressure
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Pair,                      &
     &                            status)
# endif
# if defined BULK_FLUXES || defined ECOSIM     || \
    (defined SHORTWAVE   && defined ANA_SRFLUX)
            CASE ('Tair')                     ! surface air temperature
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Tair,                      &
     &                            status)
            CASE ('Hair')                     ! surface air humidity
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Hair,                      &
     &                            status)
# endif
# if defined BULK_FLUXES || defined ECOSIM
            CASE ('Uwind')                    ! U-wind component
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            Istr, IendR, JstrR, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Uwind,                     &
     &                            status)
            CASE ('Vwind')                    ! V-wind component
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Vwind,                     &
     &                            status)
# endif
# ifdef BULK_FLUXES
            CASE ('LWrad')                    ! longwave radiation
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%lwrad,                     &
     &                            status)
            CASE ('rain')                     ! precipitation
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%rain,                      &
     &                            status)
# else
            CASE ('heat')                     ! surface net heat flux
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, JstrR, JendR,           &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%stflx(:,:,itemp),          &
     &                            status)
            CASE ('Ustr')                     ! surface U-wind stress
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            Istr, IendR, JstrR, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%sustr,                     &
     &                            status)
            CASE ('Vstr')                     ! surface V-wind stress
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%svstr,                     &
     &                            status)
# endif
# ifdef SWAN_COUPLING
#  ifdef WAVES_DIR
            CASE ('Wdir')                      ! wave direction
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Dwave,                     &
     &                            status)
#  endif
#  ifdef WAVES_HEIGHT
            CASE ('Wamp')                      ! wave height
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Hwave,                     &
     &                            status)
#  endif
#  ifdef WAVES_LENGTH
            CASE ('Wlen')                      ! wave length
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Lwave,                     &
     &                            status)
#  endif
#  ifdef WAVES_TOP_PERIOD
            CASE ('Wptop')                     ! surface wave period
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Pwave_top,                 &
     &                            status)
#  endif
#  ifdef WAVES_BOT_PERIOD
            CASE ('Wpbot')                     ! bottom wave period
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%Pwave_bot,                 &
     &                            status)
#  endif
#  if defined TKE_WAVEDISS || defined WAVES_OCEAN
            CASE ('Wdiss')                     ! wave dissipation
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%wave_dissip,               &
     &                            status)
#  endif
#  if defined SVENDSEN_ROLLER
            CASE ('Wbrk')                      ! wave breaking
              CALL ROMS_import2d (ng, MyRank,                           &
     &                            id, gtype, scale, add_offset,         &
     &                            romsImpData(ifield,ng)%field,         &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fields(id)%ImpMin, Fields(id)%ImpMax, &
     &                            FORCES(ng)%wave_break,                &
     &                            status)
#  endif
# endif
          END SELECT
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS

      END SUBROUTINE ROMS_GetImportData

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_PutExportData"

      SUBROUTINE ROMS_PutExportData (MyRank, Kexp, Nexp, status)
!
!=======================================================================
!                                                                      !
!  This routines loads ROMS export data into "ROMS_ExpData" pointer    !
!  structure.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     MyRank     Local PET                                             !
!     Kexp       Export time index for 2D ROMS state arrays,           !
!                  [1:Ngrids].                                         !
!     Nexp       Export time index for 3D ROMS state arrays,           !
!                  [1:Ngrids].                                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     status     Success/failure error flag.                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ocean
      USE mod_scalars
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable definitions.
!
      integer, intent(in)  :: MyRank
      integer, intent(in)  :: Kexp(Ngrids)
      integer, intent(in)  :: Nexp(Ngrids)
      integer, intent(out) :: status
!
!  Local variable declarations.
!
      logical :: found

      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, id, ifield, j, ng

      real(r8), parameter :: Large = 1.0E+20_r8

      real(r8) :: Fmax, Fmin, add_offset, scale
# ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len= 3), dimension(2) :: op_handle
# endif
      character (len=40) :: code
      character (len=80) :: name
!
!-----------------------------------------------------------------------
!  Load export fields.
!-----------------------------------------------------------------------
!
!  Loop over all nested, composed, or mosaic grids.
!
      DO ng=1,Ngrids
        Istr =BOUNDS(ng)%Istr (MyRank)
        Iend =BOUNDS(ng)%Iend (MyRank)
        Jstr =BOUNDS(ng)%Jstr (MyRank)
        Jend =BOUNDS(ng)%Jend (MyRank)
        IstrR=BOUNDS(ng)%IstrR(MyRank)
        IendR=BOUNDS(ng)%IendR(MyRank)
        IstrU=BOUNDS(ng)%IstrU(MyRank)
        JstrR=BOUNDS(ng)%JstrR(MyRank)
        JendR=BOUNDS(ng)%JendR(MyRank)
        JstrV=BOUNDS(ng)%JstrV(MyRank)
        DO ifield=1,Nexport(Iocean)
          id=ExportID(Iocean)%val(ifield)
          code=Fields(id)%code
          name=Fields(id)%name
          scale=Fields(id)%scale
          add_offset=Fields(id)%AddOffset
          Fmin= Large
          Fmax=-Large
          found=.FALSE.

          SELECT CASE (TRIM(ADJUSTL(code)))
            CASE ('SST')                     ! sea surface temperature
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  romsExpData(ifield,ng)%field(i,j)=                    &
     &                OCEAN(ng)%t(i,j,N(ng),Nexp(ng),itemp)+add_offset
                  Fmin=MIN(Fmin,romsExpData(ifield,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsExpData(ifield,ng)%field(i,j))
                END DO
              END DO
            CASE ('bath')                    ! bathymetry
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  romsExpData(ifield,ng)%field(i,j)=                    &
     &                GRID(ng)%h(i,j)*scale
                  Fmin=MIN(Fmin,romsExpData(ifield,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsExpData(ifield,ng)%field(i,j))
                END DO
              END DO
            CASE ('SSH')                     ! free-surface
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  romsExpData(ifield,ng)%field(i,j)=                    &
     &                OCEAN(ng)%zeta(i,j,Kexp(ng))*scale
                  Fmin=MIN(Fmin,romsExpData(ifield,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsExpData(ifield,ng)%field(i,j))
                END DO
              END DO
            CASE ('Ubar')                    ! 2D u-momentum
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  romsExpData(ifield,ng)%field(i,j)=                    &
     &                OCEAN(ng)%ubar(i,j,Kexp(ng))*scale
                  Fmin=MIN(Fmin,romsExpData(ifield,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsExpData(ifield,ng)%field(i,j))
                END DO
              END DO
            CASE ('Vbar')                    ! 2D v-momentum
              found=.TRUE.
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  romsExpData(ifield,ng)%field(i,j)=                    &
     &                OCEAN(ng)%vbar(i,j,Kexp(ng))*scale
                  Fmin=MIN(Fmin,romsExpData(ifield,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsExpData(ifield,ng)%field(i,j))
                END DO
              END DO
          END SELECT
# ifdef DISTRIBUTE
          buffer(1)=Fmin
          buffer(2)=Fmax
          op_handle(1)='MIN'
          op_handle(2)='MAX'
          CALL mp_reduce (ng, iNLM, 2, buffer, op_handle)
          Fmin=buffer(1)
          Fmax=buffer(2)
# endif
          IF (Master.and.found) THEN
            IF ((ng.eq.1).and.(ifield.eq.1)) WRITE (stdout,'(/)')
            WRITE (stdout,10) TRIM(name), tdays(ng), Fmin, Fmax
 10         FORMAT (' PutExportData - ',a,',',t64,'t = ',f12.4,/,19x,   &
     &              '(Min = ',1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS

      END SUBROUTINE ROMS_PutExportData

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_PutExportData"

      SUBROUTINE ROMS_PutGridData (MyRank, status)
!
!=======================================================================
!                                                                      !
!  This routines loads ROMS exchange grid data into "romsGrdData"      !
!  pointer structure.                                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_scalars
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable definitions.
!
      integer, intent(in)  :: MyRank
      integer, intent(out) :: status
!
!  Local variable declarations.
!
      logical :: found

      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, igrid, j, ng

      real(r8), parameter :: Large = 1.0E+20_r8

      real(r8) :: Fmax, Fmin
# ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len= 3), dimension(2) :: op_handle
# endif
      character (len=40) :: variable
      character (len=80) :: name
!
!-----------------------------------------------------------------------
!  Load export fields.
!-----------------------------------------------------------------------
!
!  Loop over all nested, composed, or mosaic grids.
!
      DO ng=1,Ngrids
        Istr =BOUNDS(ng)%Istr (MyRank)
        Iend =BOUNDS(ng)%Iend (MyRank)
        Jstr =BOUNDS(ng)%Jstr (MyRank)
        Jend =BOUNDS(ng)%Jend (MyRank)
        IstrR=BOUNDS(ng)%IstrR(MyRank)
        IendR=BOUNDS(ng)%IendR(MyRank)
        IstrU=BOUNDS(ng)%IstrU(MyRank)
        JstrR=BOUNDS(ng)%JstrR(MyRank)
        JendR=BOUNDS(ng)%JendR(MyRank)
        JstrV=BOUNDS(ng)%JstrV(MyRank)
        DO igrid=1,NromsGrdArrays
          name=romsMesh(igrid)%name
          variable=romsMesh(igrid)%variable
          Fmin= Large
          Fmax=-Large
          found=.FALSE.

          SELECT CASE (TRIM(ADJUSTL(variable)))
            CASE ('lonr')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%lonr(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('latr')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%latr(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('lonu')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%lonu(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('latu')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%latu(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('lonv')
              found=.TRUE.
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%lonv(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('latv')
              found=.TRUE.
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%latv(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('xr')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%xr(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('yr')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%yr(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('xu')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%xu(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('yu')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%yu(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('xv')
              found=.TRUE.
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%xv(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('yv')
              found=.TRUE.
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%yv(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('angler')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%angler(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
# ifdef MASKING
            CASE ('rmask')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%rmask(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('umask')
              found=.TRUE.
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%umask(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
            CASE ('vmask')
              found=.TRUE.
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  romsGrdData(igrid,ng)%field(i,j)=GRID(ng)%vmask(i,j)
                  Fmin=MIN(Fmin,romsGrdData(igrid,ng)%field(i,j))
                  Fmax=MAX(Fmax,romsGrdData(igrid,ng)%field(i,j))
                END DO
              END DO
# endif
          END SELECT
# ifdef DISTRIBUTE
          buffer(1)=Fmin
          buffer(2)=Fmax
          op_handle(1)='MIN'
          op_handle(2)='MAX'
          CALL mp_reduce (ng, iNLM, 2, buffer, op_handle)
          Fmin=buffer(1)
          Fmax=buffer(2)
# endif
          IF (Master.and.found) THEN
            IF ((ng.eq.1).and.(igrid.eq.1)) WRITE (stdout,'(/)')
            WRITE (stdout,10)  TRIM(name), tdays(ng), Fmin, Fmax
 10         FORMAT ('   PutGridData - ',a,',',t64,'t = ',f12.4,/,19x,   &
     &              '(Min = ',1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS

      END SUBROUTINE ROMS_PutGridData

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_SetInitialize"

      SUBROUTINE ROMS_SetInitialize (comp, ImportState, ExportState,    &
     &                               clock, status)
!
!=======================================================================
!                                                                      !
!  Private routine interface to ROMS "initialize" routine. It is used  !
!  by ESMF_GridComp to initialize ROMS gridded component.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
!  Imported variable definitions.
!
      integer, intent(out) :: status

      TYPE (ESMF_GridComp), intent(inout) :: comp
      TYPE (ESMF_State), intent(inout) :: ImportState
      TYPE (ESMF_State), intent(inout) :: ExportState
      TYPE (ESMF_Clock), intent(inout) :: clock
!
!  Local variable declarations.
!
      logical :: CheckError

      integer :: MyMPIcomm, MyRank, Nnodes
!
!-----------------------------------------------------------------------
!  Call ROMS initialization routines.
!-----------------------------------------------------------------------
!
!  Querry the Virtual Machine (VM) parallel environmemt for the MPI
!  communicator handle.
!
      CALL ESMF_GridCompGet (comp,                                      &
     &                       config=config,                             &
     &                       vm=VM(Iocean),                             &
     &                       rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetInitialize',              &
     &               'inquiring VM parallel environment')) RETURN

      CALL ESMF_VMGet (VM(Iocean),                                      &
     &                 localPet=MyRank,                                 &
     &                 petCount=Nnodes,                                 &
     &                 mpiCommunicator=MyMPIcomm,                       &
     &                 rc=status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetInitialize',              &
     &               'getting mpiCommunicator handle')) RETURN
!
!  Initialize ROMS.
!
      CALL ROMS_initialize (.TRUE., mpiCOMM=MyMPIcomm)

      IF (CheckError(exit_flag, 'ROMS', 'ROMS_SetInitialize',           &
                     'initializing ROMS')) THEN
        status=ESMF_FAILURE
        RETURN
      END IF
!
!  Set-up ESMF internal clock for ROMS component.
!
      CALL ROMS_SetClock (clock, status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetInitialize',              &
     &               'setting-up clock')) RETURN
!
!  Set-up ROMS exchange grid.
!
      CALL ROMS_SetGridArrays (comp, status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetInitialize',              &
     &               'setting-up exchange grid')) RETURN
!
!  Load ROMS exchange grid arrays.
!
      CALL ROMS_PutGridData (Myrank, status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetInitialize',              &
     &               'loading exchange grid arrays')) RETURN
!
!  Set-up import/export states.
!
      CALL ROMS_SetStates (ImportState, ExportState, MyRank, status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetInitialize',              &
     &               'setting-up import/export states')) RETURN
!
!  Load export initial conditions data.
!
      CALL ROMS_PutExportData (Myrank, kstp, nstp, status)
      IF (CheckError(status, 'ESMF', 'ROMS_SetInitialize',              &
     &               'loading initial conditions export data')) RETURN
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS

      RETURN
      END SUBROUTINE ROMS_SetInitialize

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_SetRun"

      SUBROUTINE ROMS_SetRun (comp, ImportState, ExportState,           &
     &                        clock, status)
!
!=======================================================================
!                                                                      !
!  Private routine interface to ROMS "run" routine. It is used  by     !
!  ESMF_GridComp to run ROMS gridded component.                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
!
!  Imported variable definitions.
!
      integer, intent(out) :: status

      TYPE (ESMF_GridComp), intent(inout) :: comp
      TYPE (ESMF_State), intent(inout) :: ImportState
      TYPE (ESMF_State), intent(inout) :: ExportState
      TYPE (ESMF_Clock), intent(inout) :: clock
!
!  Local variable definitions.
!
      logical :: CheckError
      logical :: GetImport = .FALSE.
      logical :: PutExport = .FALSE.

      TYPE (ESMF_TIME) :: CheckTime
!
!-----------------------------------------------------------------------
!  Call ROMS time-stepping driver.
!-----------------------------------------------------------------------
!
!  Time-step only if external time clock is within ROMS start and stop
!  times. Notice that the external time clock, CurrTime(0), is ahead
!  by an external time step, TimeStep(0). Also, ROMS start time and
!  stop time are, repectively, the minimum and maximum values of all
!  nested grids (see routine ROMS_SetClock above).  Similarly, ROMS
!  time step is the maximum value of all nested grids.
!
      CheckTime=CurrTime(0)-TimeStep(Iocean)
      IF ((StartTime(Iocean) <= CheckTime).and.                         &
     &    (Checktime <= StopTime(Iocean))) THEN
!
!  Get ROMS internal clock current time.
!
        CALL ESMF_ClockGet (TimeClock(Iocean),                          &
     &                      currTime=CurrTime(Iocean),                  &
     &                      rc=status)
        IF (CheckError(status, 'ESMF', 'ROMS_SetRun',                   &
     &                 'getting ROMS current time object')) RETURN
        CALL ESMF_TimeGet (CurrTime(Iocean),                            &
     &                     yy=CurrentTime(Iocean)%year,                 &
     &                     mm=CurrentTime(Iocean)%month,                &
     &                     dd=CurrentTime(Iocean)%day,                  &
     &                     h=CurrentTime(Iocean)%hour,                  &
     &                     m=CurrentTime(Iocean)%minute,                &
     &                     s=CurrentTime(Iocean)%second,                &
     &                     timeZone=CurrentTime(Iocean)%TimeZone,       &
     &                     timeStringISOFrac=CurrentTime(Iocean)%string,&
     &                     dayOfYear=CurrentTime(Iocean)%YearDay,       &
     &                     rc=status)
        IF (CheckError(status, 'ESMF', 'ROMS_SetRun',                   &
     &                 'getting ROMS current time')) RETURN
!
!  Advance ROMS internal time clock.
!
        CALL ESMF_ClockAdvance (TimeClock(Iocean),                      &
     &                          rc=status)
        IF (CheckError(status, 'ESMF', 'ROMS_SetRun',                   &
     &                 'advancing ROMS time clock')) RETURN
!
!  Get import data. (HGA:  need to decide about time indices latter.
!  Currently, none of the imported fields require time index.)
!
        IF (GetImport) THEN
          CALL ROMS_GetImportData (MyRank, kstp, nstp, status)
          IF (CheckError(exit_flag, 'ESMF', 'ROMS_SetRun',              &
     &                   'importing data from coupled models')) RETURN
        END IF
!
!  Run ROMS for specified time interval (seconds).
!
        RunInterval=CurrTime(0)-CurrTime(Iocean)
        CALL ROMS_run (RunInterval)
        IF (CheckError(exit_flag, 'ROMS', 'ROMS_SetRun',                &
     &                 'running ROMS')) THEN
          status=ESMF_FAILURE
          RETURN
        END IF
!
!  Get import data. (HGA:  need to decide about time indices latter).
!
        IF (PutExport) THEN
          CALL ROMS_PutExportData (MyRank, kstp, nstp, status)
          IF (CheckError(exit_flag, 'ESMF', 'ROMS_SetRun',              &
     &                   'exporting data to coupled models')) RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS

      END SUBROUTINE ROMS_SetRun

# undef  ESMF_METHOD
# define ESMF_METHOD "ROMS_SetFinalize"

      SUBROUTINE ROMS_SetFinalize (comp, ImportState, ExportState,      &
     &                             MyClock, status)
!
!=======================================================================
!                                                                      !
!  Private routine interface to ROMS "finalize" routine. It is used    !
!  by ESMF_GridComp to terminate ROMS gridded component.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable definitions.
!
      integer, intent(out) :: status

      TYPE (ESMF_GridComp), intent(inout) :: comp
      TYPE (ESMF_State), intent(inout) :: ImportState
      TYPE (ESMF_State), intent(inout) :: ExportState
      TYPE (ESMF_Clock), intent(inout) :: MyClock
!
!  Local variable definitions.
!
      logical :: CheckError
!
!-----------------------------------------------------------------------
!  Initialize return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Terminate ROMS execution.  Close all NetCDF files.
!-----------------------------------------------------------------------
!
      CALL ROMS_finalize

      IF (CheckError(exit_flag, 'ROMS', 'ROMS_SetFinalize',             &
                     'finalizing ROMS')) THEN
        status=ESMF_FAILURE
        RETURN
      END IF

      END SUBROUTINE ROMS_SetFinalize
#endif
      END MODULE ESMF_ROMS
