#include "wrfhydro_cpp.h"
      MODULE HYDRO_COUPLER

#if defined COAWST_COUPLING && defined MCT_LIB
!
!svn $Id: hydro_coupler.F 756 2008-09-14 20:18:28Z jcwarner $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group      Hernan G. Arango   !
!   Licensed under a MIT/X style license                               !
!   See License_ROMS.txt                                               !
!=======================================================================
!                                                                      !
!  This module is used to communicate n exchange data between WRF_Hydro!
!  other coupled model(s) using the Model Coupling Toolkit (MCT).      !
!                                                                      !
!=======================================================================
!
!  Componenet model registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain decompositin descriptor datatype and assocoiated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field storage data types and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_indxR => indexRA
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communitcations scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY : MCT_isend => isend
      USE m_Transfer, ONLY : MCT_irecv => irecv
      USE m_Transfer, ONLY : MCT_waitr => waitrecv
      USE m_Transfer, ONLY : MCT_waits => waitsend
!
# if defined MCT_INTERP_WV2AT
!
!  Sparse Matrix DataType and associated methods.
!
      USE m_SparseMatrix, ONLY : SparseMatrix
      USE m_SparseMatrix, ONLY : SparseMatrix_init => init
      USE m_SparseMatrix, ONLY : SparseMatrix_importGRowInd =>          &
     &                           importGlobalRowIndices
      USE m_SparseMatrix, ONLY : SparseMatrix_importGColInd =>          &
     &                           importGlobalColumnIndices
      USE m_SparseMatrix, ONLY : SparseMatrix_importMatrixElts =>       &
     &                           importMatrixElements
      USE m_SparseMatrix, only : SparseMatrix_lsize => lsize
      USE m_SparseMatrix, only : SparseMatrix_clean => clean
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_init => init
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_clean => clean
!
!  Decompose matrix by row.
!
      USE m_SparseMatrixPlus, ONLY : Xonly
!
!  Matrix-Vector multiply methods.
!
      USE m_MatAttrVectMul, ONLY : MCT_MatVecMul => sMatAvMult
# endif

      implicit none
!
      PRIVATE

      PUBLIC :: initialize_hyd_coupling
      PUBLIC :: initialize_hyd_routers
# ifdef ROMS_COUPLING
      PUBLIC :: hyd2ocn_coupling
      PUBLIC :: hydfocn_coupling
# endif
      PUBLIC :: finalize_hyd_coupling
!
!  Declarations.
!
      TYPE T_GlobalSegMap_G
        TYPE(GlobalSegMap) :: GSMapHYDRO         ! GloabalSegMap variables
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:)

      TYPE T_AttrVect_G
# ifdef ROMS_COUPLING
        TYPE(AttrVect) :: hyd2ocn_AV            ! AttrVect variables
        TYPE(AttrVect) :: ocn2hyd_AV
# endif
      END TYPE T_AttrVect_G
      TYPE (T_AttrVect_G), ALLOCATABLE :: AttrVect_G(:)

# ifdef ROMS_COUPLING
      TYPE T_Router_O
        type(Router)   :: HYDROtoROMS            ! Router variables
      END TYPE T_Router_O
      TYPE (T_Router_O), ALLOCATABLE :: Router_O(:,:)
# endif

      CONTAINS

      SUBROUTINE INITIALIZE_HYD_COUPLING (ih)
!
!=======================================================================
!                                                                      !
!  Initialize hydro and ocean models coupling stream.  This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE module_NoahMP_hrldas_driver
      USE module_mpp_land, only: global_nx, global_ny
!
!      include 'mpif.h'
      integer, intent(in) :: ih
!
!  Local variable declarations.
!
      integer :: MyError, MyRank
      integer :: gsmsize, nprocs
      integer :: is, ie, js, je, jc
      integer :: Isize, Jsize
      integer :: cid, cad
      character (len=70)  :: nc_name
      character (len=20)  :: to_add
      character (len=120) :: hostring
      character (len=120) :: ohstring

      real :: cff

      integer, allocatable :: start(:), length(:)
      Nhyd_grids=1
!
!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (HYDRO_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (HYDRO_COMM_WORLD, nprocs, MyError)
!
!  Initialize MCT coupled model registry.
!
      IF (ih.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Nhyd_grids))
        ALLOCATE(AttrVect_G(Nhyd_grids))
      END IF
!
      HYDid=hydids(ih)
      IF (Nhyd_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      HYDRO_COMM_WORLD,myids=hydids)
      ELSE
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      HYDRO_COMM_WORLD,HYDid)
      END IF
!
!  Initialize a Global Segment Map for non-haloed transfer of data for
!  HYDRO. Determine non-haloed start and length arrays for this
!  processor.
!
!  global_nx is total number of x cells, global_ny is total number y cells
!  xstart xend and ystart yend are local indices on this tile
      is = xstart
      ie = xend
      js = ystart
      je = yend
!
!  Determine tile size
!
      Isize=ie-is+1
      Jsize=je-js+1
      allocate( start(Jsize) )
      allocate( length(Jsize) )
      jc=0
      DO j=js,je
        jc=jc+1
        start(jc)=(j-1)*(global_nx)+is
        length(jc)=Isize
      END DO
      gsmsize=Isize*Jsize
      CALL GlobalSegMap_init (GlobalSegMap_G(ih)%GSMapHYDRO, start,     &
     &                        length, 0, HYDRO_COMM_WORLD, HYDid)
      deallocate (start)
      deallocate (length)

!
!  Initialize attribute vector holding the export data code strings of
!  the hydro model.
!
      cad=LEN(hostring)
      DO i=1,cad
        hostring(i:i)=''
      END DO
      cid=1
!
# ifdef ROMS_COUPLING
      to_add='QRIVER'
      cad=LEN_TRIM(to_add)
      write(hostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WATERLEVEL'
      cad=LEN_TRIM(to_add)
      write(hostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
!  Finalize and remove trailing spaces from the hostring
!  for the rlist.
!
      cad=LEN_TRIM(hostring)
      hostring=hostring(1:cad)
!
      CALL AttrVect_init(AttrVect_G(ih)%hyd2ocn_AV,                     &
     &                   rList=TRIM(hostring),lsize=gsmsize)
      CALL AttrVect_zero(AttrVect_G(ih)%hyd2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      cad=LEN(ohstring)
      DO i=1,cad
        ohstring(i:i)=''
      END DO
      cid=1
!
      to_add='BATH'
      cad=LEN_TRIM(to_add)
      write(ohstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':ZETA'
      cad=LEN_TRIM(to_add)
      write(ohstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
!  Finalize and remove trailing spaces from the ohstring
!  for the rlist.
!
      cad=LEN_TRIM(ohstring)
      ohstring=ohstring(1:cad)
!
      CALL AttrVect_init (AttrVect_G(ih)%ocn2hyd_AV,                    &
     &                    rList=TRIM(ohstring),lsize=gsmsize)
      CALL AttrVect_zero (AttrVect_G(ih)%ocn2hyd_AV)
!
# endif

      RETURN
      END SUBROUTINE INITIALIZE_HYD_COUPLING

      SUBROUTINE INITIALIZE_HYD_ROUTERS
!
!=======================================================================
!                                                                      !
!  Initialize waves routers for wave model.                            !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
!
!      include 'mpif.h'
!
!  Local variable declarations.
!
      integer :: MyError, MyRank
      integer :: ng, ih
!
!  Initialize MCT Routers.
!
# ifdef ROMS_COUPLING
      ALLOCATE(Router_O(Nhyd_grids,Nocn_grids))
!
!  Initialize a router to the ocean model component.
!
      DO ng=1,Nocn_grids
        DO ih=1,Nhyd_grids
          OCNid=ocnids(ng)
          CALL Router_init (OCNid, GlobalSegMap_G(ih)%GSMapHYDRO,       &
     &                    HYDRO_COMM_WORLD, Router_O(ih,ng)%HYDROtoROMS)
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE INITIALIZE_HYD_ROUTERS

# ifdef ROMS_COUPLING
      SUBROUTINE HYD2OCN_COUPLING (ih, io)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!  Fields exported to the OCEAN model:                                 !
!                                                                      !
!     * Wave direction mean (degrees)                                  !
!     * Wave direction peak (degrees)                                  !
!                                                                      !
!  Fields imported from the OCEAN Model:                               !
!                                                                      !
!     * Bathymetry, bottom elevation (m)                               !
!     * Free-surface, water surface elevation (m)                      !
!     * Depth integrated u-momentum (m/s)                              !
!     * Depth integrated v-momentum (m/s)                              !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE module_NoahMP_hrldas_driver
!
      implicit none
!
!  Imported variable declarations.
!
      integer :: ih, io
!
!  Local variable declarations.
!
      integer :: MyStatus, MyError, MySize, MyRank
      integer :: is, ie, js, je, ij, gsmsize, Tag
      integer :: Isize, Jsize, NPROCS
!      integer, pointer :: points(:)

      real, pointer :: avdata(:)
!
!-----------------------------------------------------------------------
!  Send hydro fields to ROMS.
!-----------------------------------------------------------------------
!
!  ix is number of global x cells, iy is number global y cells
!  xstart xend and ystart yend are local indices on this tile
      is = xstart
      ie = xend
      js = ystart
      je = yend
!
!  Determine tile size
!
      Isize=ie-is+1
      Jsize=je-js+1
!
      CALL MPI_COMM_RANK (HYDRO_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (HYDRO_COMM_WORLD, NPROCS, MyError)
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(ih)%GSMapHYDRO,         &
     &                           HYDRO_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      avdata=0.0
!
!  Ask for points in this tile.
!
!      CALL GlobalSegMap_Ordpnts (GlobalSegMap_G(ih)%GSMapHYDRO,         &
!     &                           MyRank, points)
!
!-----------------------------------------------------------------------
!  Send Qriver to ROMS.
!-----------------------------------------------------------------------
!     QRIVER        River discharge (m3 s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          avdata(ij)=0.0
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ih)%hyd2ocn_AV,             &
     &                             "QRIVER",avdata,gsmsize)
!
!-----------------------------------------------------------------------
!  Send Water levels to ROMS.
!-----------------------------------------------------------------------
!     WATERLEVEL        Water level elevation (m).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          avdata(ij)=0.0
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ih)%hyd2ocn_AV,             &
     &                             "WATERLEVEL",avdata,gsmsize)
!-----------------------------------------------------------------------
!
      Tag=io*100+ih*10+0
      CALL MCT_isend (AttrVect_G(ih)%hyd2ocn_AV,                        &
     &                Router_O(ih,io)%HYDROtoROMS, Tag)
      CALL MCT_waits (Router_O(ih,io)%HYDROtoROMS)
      IF (MyRank.EQ.0) THEN
        WRITE (*,36)' == HYDRO grid ',ih,                               &
     &                   ' sent hydro data to ROMS grid ', io
 36     FORMAT (a14,i2,a29,i2)
      END IF
      IF (MyError.ne.0) THEN
        WRITE (*,*)'coupling send fail hydrocplr, Error= ', MyError
        CALL FINALIZE_HYD_COUPLING(ih)
      END IF

      deallocate (avdata)
!     deallocate (avdata, points)
!
      RETURN
      END SUBROUTINE HYD2OCN_COUPLING
# endif
# ifdef ROMS_COUPLING
      SUBROUTINE HYDFOCN_COUPLING (ng, ih)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!  Fields exported to the OCEAN model:                                 !
!                                                                      !
!     * Qriver                                                         !
!     * Waterlevels                                                    !
!                                                                      !
!  Fields imported from the OCEAN Model:                               !
!                                                                      !
!     * Bathymetry, bottom elevation (m)                               !
!     * Free-surface, water surface elevation (m)                      !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE module_NoahMP_hrldas_driver
!
      implicit none
!
!  Imported variable declarations.
!
      integer :: ng, ih
!
!  Local variable declarations.
!
      integer :: MyStatus, MyError, MySize, MyRank
      integer :: gsmsize, indx, Tag
      integer :: is, ie, js, je, IP
      integer :: Isize, Jsize, INDXG, NPROCS
      integer :: iddep, idwlv
      integer, dimension(MPI_STATUS_SIZE,4) :: status

      real :: cff, retval
      real, parameter ::  Large = 1.0E+20
      real, dimension(2) :: range
      real(m8), pointer :: avdata(:)
!
!-----------------------------------------------------------------------
!  Send wave fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (HYDRO_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (HYDRO_COMM_WORLD, NPROCS, MyError)
!
!  Get the number of grid point on this processor.
!
      gsmsize=GlobalSegMap_lsize(GlobalSegMap_G(ih)%GSMapHYDRO,         &
     &                           HYDRO_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( avdata(gsmsize),stat=ierr )
      avdata=0.0_m8
!
!-----------------------------------------------------------------------
!  Receive from ROMS: Depth, Water Level, VELX, and VELY.
!-----------------------------------------------------------------------
!
!  Schedule receiving field from ocean model.
!
      Tag=ng*100+ih*10+1
      CALL MCT_irecv (AttrVect_G(ih)%ocn2hyd_AV,                        &
     &                Router_O(ih,ng)%HYDROtoROMS, Tag)
!
!     Wait to make sure the OCN data has arrived.
!
      CALL MCT_waitr (AttrVect_G(ih)%ocn2hyd_AV,                        &
     &                 Router_O(ih,ng)%HYDROtoROMS)
!
      IF (MyRank.EQ.0) THEN
        WRITE (*,35) ' == HYDRO grid ',ih,                              &
     &                    ' recv data from ROMS grid ', ng
      END IF
      IF (MyError.ne.0) THEN
       WRITE (*,*) 'coupling fail hydrocplr, MyStatus= ', MyError
        CALL FINALIZE_HYD_COUPLING(ih)
      END IF
 35   FORMAT (a14,i2,a26,i2)
!
!  Extract data and load into model arrays.
!
!  Bottom elevation.
!
        CALL AttrVect_exportRAttr (AttrVect_G(ih)%ocn2hyd_AV,           &
     &                             "BATH",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        IP=0
        DO js=ystart,yend
          DO is=xstart,xend
            IP=IP+1
!           some hydr var =avdata(IP)
            range(1)=MIN(range(1),REAL(avdata(IP)))
            range(2)=MAX(range(2),REAL(avdata(IP)))
          END DO
        END DO
        retval=0.0
        CALL mpi_allreduce ( range(1), retval , 1, MPI_REAL,            &
     &                       MPI_MIN, HYDRO_COMM_WORLD, ierr )
        range(1)=retval
        CALL mpi_allreduce ( range(2), retval , 1, MPI_REAL,            &
     &                       MPI_MAX, HYDRO_COMM_WORLD, ierr )
        range(2)=retval
        IF (MyRank.eq.0) THEN
          write(*,40) 'ROMStoHYDRO Min/Max BATH   (m):     ',           &
     &                      range(1),range(2)
        END IF
!
!  Water surface elevation.
!
        CALL AttrVect_exportRAttr (AttrVect_G(ih)%ocn2hyd_AV,           &
     &                             "ZETA",avdata,gsmsize)
        range(1)= Large
        range(2)=-Large
        IP=0
        DO js=ystart,yend
          DO is=xstart,xend
            IP=IP+1
!  Here is where we load the zeta from roms into a Hydro variable
!
!           SFCHEADRT(is,js) = REAL(avdata(IP))
!
            range(1)=MIN(range(1),REAL(avdata(IP)))
            range(2)=MAX(range(2),REAL(avdata(IP)))
          END DO
        END DO
        CALL mpi_allreduce ( range(1), retval , 1, MPI_REAL,            &
     &                       MPI_MIN, HYDRO_COMM_WORLD, ierr )
        range(1)=retval
        CALL mpi_allreduce ( range(2), retval , 1, MPI_REAL,            &
     &                       MPI_MAX, HYDRO_COMM_WORLD, ierr )
        range(2)=retval
        IF (MyRank.eq.0) THEN
          write(*,40) 'ROMStoHYDRO Min/Max WLEV    (m):     ',          &
     &                      range(1),range(2)
        END IF
!
      deallocate (avdata)
!
 40   FORMAT (a36,1x,2(1pe14.6))

      RETURN
      END SUBROUTINE HYDFOCN_COUPLING
# endif

      SUBROUTINE FINALIZE_HYD_COUPLING(ih)
!
!=======================================================================
!                                                                    ===
!  This routines terminates execution during coupling error.         ===
!                                                                    ===
!=======================================================================
      USE mct_coupler_params
!     USE module_NoahMP_hrldas_driver
!
!  Local variable declarations.
!
      integer :: ih, io, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
# ifdef ROMS_COUPLING
      DO io=1,Nhyd_grids
        DO ih=1,Nhyd_grids
          CALL Router_clean (Router_O(ih,io)%HYDROtoROMS, MyError)
        END DO
      END DO
      CALL AttrVect_clean (AttrVect_G(ih)%hyd2ocn_AV, MyError)
# endif
      CALL GlobalSegMap_clean (GlobalSegMap_G(ih)%GSMapHYDRO, MyError)

      END SUBROUTINE FINALIZE_HYD_COUPLING
#endif
      END MODULE HYDRO_COUPLER
