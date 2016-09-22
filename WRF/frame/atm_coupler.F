#include "wrfcpp.h"
      MODULE atm_coupler_mod
#ifdef WRF_COUPLING
!svn $Id$
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group      Hernan G. Arango   !
!   Licensed under a MIT/X style license                               !
!   See License_ROMS.txt                                               !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between WRF    !
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
      USE m_Transfer, ONLY: MCT_send => send
      USE m_Transfer, ONLY: MCT_recv => recv
      USE m_Transfer, ONLY: MCT_isend => isend
      USE m_Transfer, ONLY: MCT_irecv => irecv
      USE m_Transfer, ONLY: MCT_waitr => waitrecv
      USE m_Transfer, ONLY: MCT_waits => waitsend

# if defined MCT_INTERP_WV2AT || defined MCT_INTERP_OC2AT
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
      USE m_SparseMatrix, only : SparseMatrix_GNumElem =>               &
     &                           GlobalNumElements
      USE m_SparseMatrix, only : SparseMatrix_clean => clean
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_init => init
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_clean => clean
!
!  Decompose matrix by row.
!
      USE m_SparseMatrixPlus, ONLY : Xonly
!     USE m_SparseMatrixPlus, ONLY : Yonly
!
!  Matrix-Vector multiply methods.
!
      USE m_MatAttrVectMul, ONLY : MCT_MatVecMul => sMatAvMult
# endif
!
      USE mct_wrf_coupler_params
      USE module_parallel
!
      implicit none
!
      PRIVATE

      PUBLIC :: INITIALIZE_ATM_ROUTERS
      PUBLIC :: initialize_atm_coupling
      PUBLIC :: atm_coupling
# ifdef ROMS_COUPLING
      PUBLIC :: atm2ocn_coupling
      PUBLIC :: atmfocn_coupling
      PUBLIC :: atm_coupling_aux4
# endif
# ifdef SWAN_COUPLING
      PUBLIC :: atm2wav_coupling
      PUBLIC :: atmfwav_coupling
# endif
      PUBLIC :: finalize_atm_coupling

      include 'mpif.h'
!
!  Declarations.
!

      TYPE T_GlobalSegMap_G
        TYPE(GlobalSegMap) :: GSMapWRF         ! GloabalSegMap variables
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:)

!# ifdef MCT_INTERP_WV2AT
!     TYPE T_GSMapInterp_W
!       TYPE(GlobalSegMap) :: GSMapSWAN       ! GloabalSegMap variables
!     END TYPE T_GSMapInterp_W
!     TYPE (T_GSMapInterp_W), ALLOCATABLE :: GSMapInterp_W(:,:)
!# endif
!# ifdef MCT_INTERP_OC2AT
!     TYPE T_GSMapInterp_O
!       TYPE(GlobalSegMap) :: GSMapROMS       ! GloabalSegMap variables
!     END TYPE T_GSMapInterp_O
!     TYPE (T_GSMapInterp_O), ALLOCATABLE :: GSMapInterp_O(:,:)
!# endif
!# if defined MCT_INTERP_OC2AT || defined MCT_INTERP_WV2AT
# ifdef MOVE_NESTS
      TYPE T_GSMapInterp_R
        TYPE(GlobalSegMap) :: GSMapStatic     ! GloabalSegMap variables
      END TYPE T_GSMapInterp_R
      TYPE (T_GSMapInterp_R), ALLOCATABLE :: GSMapInterp_R(:)
# endif

# if defined ROMS_COUPLING
      TYPE T_AttrVect_O
        TYPE(AttrVect) :: atm2ocn_AV            ! AttrVect variables
        TYPE(AttrVect) :: ocn2atm_AV            ! AttrVect variables
      END TYPE T_AttrVect_O
      TYPE (T_AttrVect_O), ALLOCATABLE :: AttrVect_O(:,:)
# endif
# if defined SWAN_COUPLING
      TYPE T_AttrVect_W
       TYPE(AttrVect) :: atm2wav_AV            ! AttrVect variables
       TYPE(AttrVect) :: wav2atm_AV            ! AttrVect variables
      END TYPE T_AttrVect_W
      TYPE (T_AttrVect_W), ALLOCATABLE :: AttrVect_W(:,:)
# endif

# if defined ROMS_COUPLING
      TYPE T_Router_O
        type(Router)   :: WRFtoROMS           ! Router variables
      END TYPE T_Router_O
      TYPE (T_Router_O), ALLOCATABLE :: Router_O(:,:)
# endif
# if defined SWAN_COUPLING
      TYPE T_Router_W
        type(Router)   :: WRFtoSWAN           ! Router variables
      END TYPE T_Router_W
      TYPE (T_Router_W), ALLOCATABLE :: Router_W(:,:)
# endif

!# ifdef MCT_INTERP_WV2AT
!      TYPE T_AV2_W
!        TYPE(AttrVect) :: wav2atm_AV2         ! AttrVect variables
!        TYPE(AttrVect) :: atm2wav_AV2 
!      END TYPE T_AV2_W
!      TYPE (T_AV2_W), ALLOCATABLE :: AV2_W(:,:)

!      TYPE(SparseMatrix) :: sMatA             ! Sparse matrix elements
!      TYPE(SparseMatrix) :: sMatW             ! Sparse matrix elements

!      TYPE T_SMPlus_W
!        TYPE(SparseMatrixPlus) :: W2AMatPlus  ! Sparse matrix plus elements
!        TYPE(SparseMatrixPlus) :: A2WMatPlus  ! Sparse matrix plus elements
!      END TYPE T_SMPlus_W
!      TYPE (T_SMPlus_W), ALLOCATABLE :: SMPlus_W(:,:)
!# endif

!# ifdef MCT_INTERP_OC2AT
!      TYPE T_AV2_O
!        TYPE(AttrVect) :: atm2ocn_AV2           ! AttrVect variables
!        TYPE(AttrVect) :: ocn2atm_AV2
!      END TYPE T_AV2_O
!      TYPE (T_AV2_O), ALLOCATABLE :: AV2_O(:,:)

!#  if !defined MCT_INTERP_WV2AT
!      TYPE(SparseMatrix) :: sMatA             ! Sparse matrix elements
!#  endif
!      TYPE(SparseMatrix) :: sMatO             ! Sparse matrix elements

!      TYPE T_SMPlus_O
!        TYPE(SparseMatrixPlus) :: A2OMatPlus  ! Sparse matrix plus elements
!        TYPE(SparseMatrixPlus) :: O2AMatPlus  ! Sparse matrix plus elements
!      END TYPE T_SMPlus_O
!      TYPE (T_SMPlus_O), ALLOCATABLE :: SMPlus_O(:,:)
!# endif
!# if defined MCT_INTERP_OC2AT || defined MCT_INTERP_WV2AT
# ifdef MOVE_NESTS
      TYPE T_AV2_R
        TYPE(AttrVect) :: crs2fin_AV2         ! AttrVect variables
        TYPE(AttrVect) :: fin2crs_AV2 
      END TYPE T_AV2_R
      TYPE (T_AV2_R), ALLOCATABLE :: AV2_RO(:,:)
      TYPE (T_AV2_R), ALLOCATABLE :: AV2_RW(:,:)

      TYPE(SparseMatrix) :: sMatF             ! Sparse matrix elements
      TYPE(SparseMatrix) :: sMatC             ! Sparse matrix elements

      TYPE T_SMPlus_R
        TYPE(SparseMatrixPlus) :: C2FMatPlus  ! Sparse matrix plus elements
        TYPE(SparseMatrixPlus) :: F2CMatPlus  ! Sparse matrix plus elements
      END TYPE T_SMPlus_R
      TYPE (T_SMPlus_R), ALLOCATABLE :: SMPlus_R(:)
!
!      integer :: par_grid_ratio,full_Isize,full_Jsize,parent_id

# endif

      CONTAINS

      SUBROUTINE initialize_atm_coupling(grid, ia)
!
!=======================================================================
!                                                                      !
!  Initialize waves and ocean models coupling stream.  This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
      USE module_domain

      implicit none
!
!  Imported variable definitions.
!
      TYPE(domain) , INTENT (IN) :: grid 
      integer, intent(in) :: ia

      include 'mpif.h'
!
!  Local variable declarations.  
!
      integer :: MyError, MyRank
      integer :: npoints, gsmsize, nprocs, localsize
      integer :: j, jc, Isize, Jsize, count, Asize
      integer :: i, is, ie, js, je, cid, cad, io, iw 
      integer :: nRows, nCols, num_sparse_elems

      integer, pointer :: start(:), length(:)
      character (len=120)  :: to_add, avstring

!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (ATM_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (ATM_COMM_WORLD, nprocs, MyError)
!
      IF (ia.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Natm_grids))
# if defined ROMS_COUPLING
        ALLOCATE(AttrVect_O(Natm_grids,Nocn_grids))
!#  ifdef MCT_INTERP_OC2AT
!       ALLOCATE(SMPlus_O(Natm_grids,Nocn_grids))
!       ALLOCATE(AV2_O(Natm_grids,Nocn_grids))
!       ALLOCATE(GSMapInterp_O(Natm_grids,Nocn_grids))
!#  endif
# endif
# if defined SWAN_COUPLING
        ALLOCATE(AttrVect_W(Natm_grids,Nwav_grids))
!#  ifdef MCT_INTERP_WV2AT
!       ALLOCATE(SMPlus_W(Natm_grids,Nwav_grids))
!       ALLOCATE(AV2_W(Natm_grids,Nwav_grids))
!       ALLOCATE(GSMapInterp_W(Natm_grids,Nwav_grids))
!#  endif
# endif
# ifdef MOVE_NESTS
        ALLOCATE(GSMapInterp_R(Natm_grids))
        ALLOCATE(SMPlus_R(Natm_grids))
        ALLOCATE(AV2_RO(Natm_grids,Nocn_grids))
        ALLOCATE(AV2_RW(Natm_grids,Nwav_grids))
# endif
      END IF
!
!  Initialize MCT coupled model registry.
!
      ATMid=atmids(ia)
      IF (Natm_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels,MPI_COMM_WORLD,                 &
     &                      ATM_COMM_WORLD,myids=atmids)
      ELSE
        CALL MCTWorld_init (N_mctmodels,MPI_COMM_WORLD,                 &
     &                      ATM_COMM_WORLD,ATMid)
      END IF
!
!  Initialize a Global Segment Map for non-haloed transfer of data out
!  of WRF. Determine non-haloed start and length arrays for this
!  processor.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
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
        start(jc)=(j-1)*(grid%ed31-1)+is
        length(jc)=Isize
      END DO
      gsmsize=Isize*Jsize
!
      CALL GlobalSegMap_init (GlobalSegMap_G(ia)%GSMapWRF,              &
     &                        start, length, 0, ATM_COMM_WORLD, ATMid)
      deallocate(start)
      deallocate(length)
!
# ifdef MOVE_NESTS
!  Save some grid information.
      par_grid_ratio(ia)=grid%parent_grid_ratio
      full_Isize(ia)=grid%e_we-1
      full_Jsize(ia)=grid%e_sn-1
      parent_id(ia)=grid%parent_id
!
      IF ((Myrank.eq.0).and.(moving_nest(ia).eq.1)) THEN
!  For staitc grids, that interpolation occurs in the mct_roms_wrf.h.
!  For moving WRF nests,for now, this will take effect here and
!  we need to interpolate the small moving nest to a defined
!  static larger grid that is the size of
!  parent_grid_ratio . num_cellsx_par . num_cellsy_par
!
!!!!!!!!!!!!!!!!!!!!!!
! First work on atm_moving nest to static grid.
! This really means to interpolate the moving nest to a larger
! static grid that is a factor par_grid_ratio larger.
!!!!!!!!!!!!!!!!!!!!!!
!

          src_grid_dims(1)=full_Isize(ia)
          src_grid_dims(2)=full_Jsize(ia)
          dst_grid_dims(1)=full_Isize(parent_id(ia))*par_grid_ratio(ia)
          dst_grid_dims(2)=full_Jsize(parent_id(ia))*par_grid_ratio(ia)
!
! Init the sparse matrix. Rows = dst grid cell indices (these change).
!                         Cols = src grid cell indices.
!
          num_sparse_elems=src_grid_dims(1)*src_grid_dims(2)
          nRows=dst_grid_dims(1)*dst_grid_dims(2)
          nCols=src_grid_dims(1)*src_grid_dims(2)
!
          allocate ( sparse_rows(num_sparse_elems) )
          allocate ( sparse_cols(num_sparse_elems) )
          allocate ( sparse_weights(num_sparse_elems) )
!
! Create sparse matrix.
!
          DO i=1,num_sparse_elems
            sparse_cols(i)=i
!           sparse_weights(i)=1.0
            sparse_weights(i)=0.0
          END DO

!  Do not use first 3 (or 5) rows around child perimeter.
          jc=0
          DO j=par_grid_ratio(ia)+1,src_grid_dims(2)-par_grid_ratio(ia)
           DO i=par_grid_ratio(ia)+1,src_grid_dims(1)-par_grid_ratio(ia)
              jc=(j-1)*(src_grid_dims(1))+i
              sparse_weights(jc)=1.0
            END DO
          END DO

          is=(grid%i_parent_start-1)*par_grid_ratio(ia)+1
          js=(grid%j_parent_start-1)*par_grid_ratio(ia)+1
          jc=0
          DO j=1,src_grid_dims(2)
            count=(js-1+j-1)*(dst_grid_dims(1))+is-1
            DO i=1,src_grid_dims(1)
              jc=jc+1
              sparse_rows(jc)=count+i
            END DO
          END DO

          call SparseMatrix_init(sMatF,nRows,nCols,num_sparse_elems)
          call SparseMatrix_importGRowInd(sMatF, sparse_rows,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importGColInd(sMatF, sparse_cols,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importMatrixElts(sMatF, sparse_weights,     &
     &                                    num_sparse_elems)
!
! Deallocate arrays.
!
          deallocate ( sparse_rows )
          deallocate ( sparse_cols )
          deallocate ( sparse_weights )
!
!!!!!!!!!!!!!!!!!!!!!!
! Second work on full static grid to moving nest.
!!!!!!!!!!!!!!!!!!!!!!
!
          dst_grid_dims(1)=full_Isize(ia)
          dst_grid_dims(2)=full_Jsize(ia)
          src_grid_dims(1)=full_Isize(parent_id(ia))*par_grid_ratio(ia)
          src_grid_dims(2)=full_Jsize(parent_id(ia))*par_grid_ratio(ia)
!
! Init the sparse matrix. Rows = dst grid cell indices.
!                         Cols = src grid cell indices (these change).
!
          num_sparse_elems=dst_grid_dims(1)*dst_grid_dims(2)
          nRows=dst_grid_dims(1)*dst_grid_dims(2)
          nCols=src_grid_dims(1)*src_grid_dims(2)
!
          allocate ( sparse_rows(num_sparse_elems) )
          allocate ( sparse_cols(num_sparse_elems) )
          allocate ( sparse_weights(num_sparse_elems) )
!
! Create sparse matrix.
!
          DO i=1,num_sparse_elems
            sparse_rows(i)=i
            sparse_weights(i)=1.0
          END DO
          is=(grid%i_parent_start-1)*par_grid_ratio(ia)+1
          js=(grid%j_parent_start-1)*par_grid_ratio(ia)+1
          jc=0
          DO j=1,dst_grid_dims(2)
            count=(js-1+j-1)*(src_grid_dims(1))+is-1
            DO i=1,dst_grid_dims(1)
!              count=(js-1)*(src_grid_dims(1))+is+1
              jc=jc+1
              sparse_cols(jc)=count+i
            END DO
          END DO
!
          call SparseMatrix_init(sMatC,nRows,nCols,num_sparse_elems)
          call SparseMatrix_importGRowInd(sMatC, sparse_rows,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importGColInd(sMatC, sparse_cols,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importMatrixElts(sMatC, sparse_weights,     &
     &                                    num_sparse_elems)
!
! Deallocate arrays.
!
          deallocate ( sparse_rows )
          deallocate ( sparse_cols )
          deallocate ( sparse_weights )
        END IF
        IF (moving_nest(ia).eq.1) THEN
!
!  Create Global Seg Map for the static larger grid.
!  Determine start and lengths for domain decomposition.
!
          dst_grid_dims(1)=full_Isize(parent_id(ia))*par_grid_ratio(ia)
          dst_grid_dims(2)=full_Jsize(parent_id(ia))*par_grid_ratio(ia)
          Isize=INT(dst_grid_dims(1)/nprocs)
          IF (MyRank.eq.nprocs-1) THEN
            Isize=dst_grid_dims(1)-Isize*(nprocs-1)
          ENDIF
          allocate( start(1) )
          allocate( length(1) )
          start=(MyRank*INT(dst_grid_dims(1)/nprocs))*dst_grid_dims(2)+1
          length=Isize*dst_grid_dims(2)
!
          CALL GlobalSegMap_init (GSMapInterp_R(ia)%GSMapStatic,        &
     &                            start,length,0,ATM_COMM_WORLD,ATMid)
          deallocate (start)
          deallocate (length)
!
! Create ATM sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
          call SparseMatrixPlus_init(SMPlus_R(ia)%F2CMatPlus, sMatF,    &
     &                               GlobalSegMap_G(ia)%GSMapWRF,       &
     &                               GSMapInterp_R(ia)%GSMapStatic,     &
     &                               Xonly, 0, ATM_COMM_WORLD, ATMid)
!         call SparseMatrix_clean(sMatA)
!
! Create Ocean sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
          call SparseMatrixPlus_init(SMPlus_R(ia)%C2FMatPlus, sMatC,    &
     &                               GSMapInterp_R(ia)%GSMapStatic,     &
     &                               GlobalSegMap_G(ia)%GSMapWRF,       &
     &                               Xonly, 0, ATM_COMM_WORLD, ATMid)
!         call SparseMatrix_clean(sMatO)
        END IF
# endif
# ifdef ROMS_COUPLING
!
!  Initialize attribute vector holding the export data code strings of
!  the atmosphere model.
      cad=LEN(avstring)
      DO j=1,cad
        avstring(j:j)=''
      END DO
      cid=1
!
      to_add='GSW'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':GLW'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
#  ifdef ATM2OCN_FLUXES
      to_add=':LH'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':HFX'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':USTRESS'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VSTRESS'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      to_add=':MSLP'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if defined BULK_FLUXES || defined ECOSIM || \
   (defined SHORTWAVE && defined ANA_SRFLUX)
      to_add=':RELH'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':T2'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if defined BULK_FLUXES || defined ECOSIM
      to_add=':U10'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':V10'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  ifdef CLOUDS
      to_add=':CLDFRA'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if !defined ANA_RAIN && defined EMINUSP
      to_add=':RAIN'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if defined EMINUSP
      to_add=':EVAP'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
      cad=LEN_TRIM(avstring)
      avstring=avstring(1:cad)
!
      DO io=1,Nocn_grids
        OCNid=ocnids(io)
        CALL AttrVect_init (AttrVect_O(ia,io)%atm2ocn_AV,               &
     &                      rlist=TRIM(avstring),lsize=gsmsize)
        CALL AttrVect_zero (AttrVect_O(ia,io)%atm2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
#  ifdef MCT_INTERP_OC2AT
        CALL AttrVect_init (AttrVect_O(ia,io)%ocn2atm_AV,               &
     &                      rList="SST:CPL_MASK",lsize=gsmsize)
#  else
        CALL AttrVect_init (AttrVect_O(ia,io)%ocn2atm_AV,               &
     &                      rList="SST",lsize=gsmsize)
#  endif
        CALL AttrVect_zero (AttrVect_O(ia,io)%ocn2atm_AV)

#  ifdef MOVE_NESTS
        IF (moving_nest(ia).eq.1) THEN
!
!  Initialize attribute vector holding the export data of
!  the atm model. The Asize is the number of grid point on this
!  processor.
!
          Asize=GlobalSegMap_lsize(GSMapInterp_R(ia)%GSMapStatic,       &
     &                             ATM_COMM_WORLD)
          CALL AttrVect_init (AV2_RO(ia,io)%crs2fin_AV2,                &
     &                        rList="SST:CPL_MASK",lsize=Asize)
          CALL AttrVect_zero (AV2_RO(ia,io)%crs2fin_AV2)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
          CALL AttrVect_init (AV2_RO(ia,io)%fin2crs_AV2,                &
     &                        rlist=TRIM(avstring),lsize=Asize)
          CALL AttrVect_zero (AV2_RO(ia,io)%fin2crs_AV2)
        END IF
#  endif
      END DO
# endif
# ifdef SWAN_COUPLING
      DO iw=1,Nwav_grids
        WAVid=wavids(iw)
!
!  Initialize attribute vector holding the export data code string of
!  the wave model.
!
        CALL AttrVect_init (AttrVect_W(ia,iw)%atm2wav_AV,               &
     &                      rlist="U10:V10",lsize=gsmsize)
        CALL AttrVect_zero (AttrVect_W(ia,iw)%atm2wav_AV)
!
!  Initialize attribute vector for data from SWAN.
!
#  ifdef MCT_INTERP_WV2AT
        CALL AttrVect_init (AttrVect_W(ia,iw)%wav2atm_AV,               &
     &                      rList="HSIGN:WLENP:RTP:CPL_MASK",           &
     &                      lsize=gsmsize)
#  else
        CALL AttrVect_init (AttrVect_W(ia,iw)%wav2atm_AV,               &
     &                      rList="HSIGN:WLENP:RTP",lsize=gsmsize)
#  endif
        CALL AttrVect_zero (AttrVect_W(ia,iw)%wav2atm_AV)
#  ifdef MOVE_NESTS
        IF (moving_nest(ia).eq.1) THEN
!
!  Initialize attribute vector holding the export data of
!  the atm model. The Asize is the number of grid point on this
!  processor.
!
          Asize=GlobalSegMap_lsize(GSMapInterp_R(ia)%GSMapStatic,       &
     &                             ATM_COMM_WORLD)
          CALL AttrVect_init (AV2_RW(ia,iw)%crs2fin_AV2,                &
     &                        rList="HSIGN:WLENP:RTP:CPL_MASK",         &
     &                        lsize=Asize)
          CALL AttrVect_zero (AV2_RW(ia,iw)%crs2fin_AV2)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
          CALL AttrVect_init (AV2_RW(ia,iw)%fin2crs_AV2,                &
     &                        rList="U10:V10",lsize=Asize)
          CALL AttrVect_zero (AV2_RW(ia,iw)%fin2crs_AV2)
        END IF
#  endif
      END DO
# endif

      RETURN
      END SUBROUTINE initialize_atm_coupling

      SUBROUTINE INITIALIZE_ATM_ROUTERS
!
!=======================================================================
!                                                                      !
!  Initialize atm routers.                                             !
!                                                                      !
!=======================================================================

      implicit none
!
!  Local variable declarations.
!
      integer :: io, iw, ia
!
!  Initialize MCT Routers.
!
# ifdef ROMS_COUPLING
      ALLOCATE(Router_O(Natm_grids,Nocn_grids))
!
!  Initialize a router to the ocean model component.
!
      DO io=1,Nocn_grids
        DO ia=1,Natm_grids
          OCNid=ocnids(io)
#  ifdef MOVE_NESTS
          IF (moving_nest(ia).eq.1) THEN
            CALL Router_init (OCNid, GSMapInterp_R(ia)%GSMapStatic,     &
     &                        ATM_COMM_WORLD, Router_O(ia,io)%WRFtoROMS)
          ELSE
#  endif
            CALL Router_init (OCNid, GlobalSegMap_G(ia)%GSMapWRF,       &
     &                        ATM_COMM_WORLD, Router_O(ia,io)%WRFtoROMS)
#  ifdef MOVE_NESTS
          END IF 
#  endif
        END DO
      END DO
# endif
# ifdef SWAN_COUPLING
      ALLOCATE(Router_W(Natm_grids,Nwav_grids))
!
!  Initialize a router to the wav model component.
!
      DO ia=1,Natm_grids
        DO iw=1,Nwav_grids
          WAVid=wavids(iw)
#  ifdef MOVE_NESTS
          IF (moving_nest(ia).eq.1) THEN
            CALL Router_init (WAVid, GSMapInterp_R(ia)%GSMapStatic,     &
     &                        ATM_COMM_WORLD, Router_W(ia,iw)%WRFtoSWAN)
          ELSE
#  endif
            CALL Router_init (WAVid, GlobalSegMap_G(ia)%GSMapWRF,       &
     &                        ATM_COMM_WORLD, Router_W(ia,iw)%WRFtoSWAN)
#  ifdef MOVE_NESTS
          END IF
#  endif
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE INITIALIZE_ATM_ROUTERS

      SUBROUTINE atm_coupling(grid, num_steps)
!
!=======================================================================
!                                                                      !
!  Determine which other model to do the coupling at this time step.   !
!                                                                      !
!=======================================================================
!
      USE module_domain

      implicit none

!
!  Imported variable definitions.
!
      integer, intent(in) :: num_steps
!     TYPE(domain) , INTENT (IN) :: grid
      TYPE(domain) , POINTER     :: grid
      TYPE(domain) , POINTER     :: grid_ptr
!
!  Local variable declarations.  
!
      integer :: io, iw, ia, offset
      integer :: kid, num_ksteps
!
      IF (num_steps.eq.0) THEN
        offset=0
      ELSE
!       offset=-1
        offset=0
      END IF
# ifdef ROMS_COUPLING
      IF (MOD(num_steps+offset, nATM2OCN(1,1)).eq.0) THEN
        DO io=1,Nocn_grids
          ia=1
          CALL atm2ocn_coupling(grid,ia,io)
          DO ia=2,Natm_grids
            CALL find_grid_by_id(ia, grid, grid_ptr)
            CALL atm2ocn_coupling(grid_ptr,ia,io)
          END DO
        END DO
      END IF
!
      IF (MOD(num_steps+offset, nATMFOCN(1,1)).eq.0) THEN
        DO io=1,Nocn_grids
          ia=1
          CALL atmfocn_coupling(grid,ia,io,0)
          DO ia=2,Natm_grids
            CALL find_grid_by_id(ia, grid, grid_ptr)
            CALL atmfocn_coupling(grid_ptr,ia,io,0)
          END DO
        END DO
      END IF
# endif
# ifdef SWAN_COUPLING
      IF (MOD(num_steps+offset, nATM2WAV(1,1)).eq.0) THEN
        DO iw=1,Nwav_grids
          ia=1
          CALL atm2wav_coupling(grid,ia,iw)
          DO ia=2,Natm_grids
            CALL find_grid_by_id(ia, grid, grid_ptr)
            CALL atm2wav_coupling(grid_ptr,ia,iw)
          END DO
        END DO
      END IF
!
      IF (MOD(num_steps+offset, nATMFWAV(1,1)).eq.0) THEN
        DO iw=1,Nwav_grids
          ia=1
          CALL atmfwav_coupling(grid,ia,iw)
          DO ia=2,Natm_grids
            CALL find_grid_by_id(ia, grid, grid_ptr)
            CALL atmfwav_coupling(grid_ptr,ia,iw)
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE atm_coupling

# ifdef ROMS_COUPLING
      SUBROUTINE atm_coupling_aux4(grid, num_steps)
!
!=======================================================================
!                                                                      !
!  Determine which other model to do the coupling at this time step.   !
!                                                                      !
!=======================================================================
!
      USE module_domain

      implicit none
!
!  Imported variable definitions.
!
      integer, intent(in) :: num_steps
      TYPE(domain) , POINTER     :: grid 
      TYPE(domain) , POINTER     :: grid_ptr
!
!  Local variable declarations.  
!
      integer :: io, ia
      integer :: kid, num_ksteps
!
      DO io=1,Nocn_grids
        ia=1
        CALL atmfocn_coupling (grid,ia,io,1)
        DO ia=2,Natm_grids
          CALL find_grid_by_id(ia, grid, grid_ptr)
          CALL atmfocn_coupling(grid_ptr,ia,io,1)
        END DO
      END DO

      RETURN
      END SUBROUTINE atm_coupling_aux4
# endif

# ifdef ROMS_COUPLING
      SUBROUTINE atm2ocn_coupling (grid, ia, io)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupled data streams.          !
!  Currently, the following data streams are processed:                !
!                                                                      !
!  Possible fields exported to the OCEAN Model:                        !
!                                                                      !
!     * GSW        Short wave raditaion  (Watts/m2)                    !
!     * GLW        Long wave raditaion  (Watts/m2)                     !
!     * LH         Latent heat flux     (Watts/m2)                     !
!     * HFX        Sensible heat flux   (Watts/m2)                     !
!     * USTRESS    Surface U-wind stress (Pa)                          !
!     * VSTRESS    Surface v-stress      (Pa)                          !
!     * MSLP       Mean Sea Level Pressure (Pa)                        !
!     * RELH       Surface air relative humidity (percent)             !
!     * T2         Surface 2m air temperature (Celsius)                !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!     * CLDFRA     Cloud fraction       (percent/100)                  !
!     * RAIN       Precipitation        (m/s)                          !
!     * EVAP       Evaporation          (m/s)                          !
!                                                                      !
!  Fields exported to the WAVE Model:                                  !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!                                                                      !
!  Fields acquired from the WAVE Model:                                !
!                                                                      !
!     * HISGN      Significant wave heigth (m)                         !
!     * WLENP      Peak wave length (m)                                !
!     * RTP        Peak wave period (s)                                !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * SST        Sea surface temperature                             !
!=======================================================================
!
      USE module_domain
!
      implicit none
!     TYPE(domain) , INTENT (IN) :: grid
      TYPE(domain) , POINTER     :: grid            ! jcw check this
#  ifdef MOVE_NESTS
      TYPE(domain) , POINTER     :: grid_ptr
#  endif
      integer, intent(in) :: ia, io
!
!  Local variable declarations.
!
      integer :: is, ie, js, je, ij, Tag
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, nprocs, offset
#  ifdef MOVE_NESTS
      integer :: count, jc, num_sparse_elems
#  endif

      real, parameter :: eps=1.0e-10
      real :: cff, cff1, cff2, cff3, rnum, rden, c04, c05
      real, pointer :: AA(:)
#  ifdef MOVE_NESTS
      real, pointer :: MOV_MASK(:,:)
#  endif
!
!  Set grid range.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
      Isize=ie-is+1
      Jsize=je-js+1
!
!-----------------------------------------------------------------------
!  Send atmosphere fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (ATM_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (ATM_COMM_WORLD, nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ia)%GSMapWRF,             &
     &                         ATM_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( AA(Asize), stat=ierr )
!
#  ifdef MOVE_NESTS
!  If we have a moving nest, then need to mask out parent location of
!  the child nest.
!
      allocate ( MOV_MASK(is:ie,js:je), stat=ierr )
      DO i=is,ie
        DO j=js,je
          MOV_MASK(i,j)=1.0
        END DO
      END DO
!  Get Istr and Jstr of moving child in parent grid.
!  Allow parent to provide one more row around child perimeter.


!  want to know:
!  1 - do you have a child that is moving?
!  2  - if so, then where is it now?

        IF (ia.eq.parent_id(Natm_grids)) THEN
          CALL find_grid_by_id(Natm_grids, grid, grid_ptr)
          offset=1
          Istr=grid_ptr%i_parent_start+offset
          Jstr=grid_ptr%j_parent_start+offset
          ij=full_Isize(Natm_grids)/par_grid_ratio(Natm_grids)
          Iend=grid_ptr%i_parent_start+ij-1-offset
          ij=full_Jsize(Natm_grids)/par_grid_ratio(Natm_grids)
          Jend=grid_ptr%j_parent_start+ij-1-offset
!
          IF ((Istr.ge.is).or.(Istr.le.ie)) THEN
            IF ((Jstr.ge.js).or.(Jstr.le.je)) THEN
              DO i=MAX(is,Istr),MIN(Iend,ie)
                DO j=MAX(js,Jstr),MIN(Jend,je)
                  MOV_MASK(i,j)=0.0
                END DO
              END DO
            END IF
          END IF
        END IF
!
!  Adjust the sparse matrices here.
!
!       IF (ia.eq.Natm_grids) THEN
        IF ((Myrank.eq.0).and.(moving_nest(ia).eq.1)) THEN

!  First the fine to coarse
          src_grid_dims(1)=full_Isize(ia)
          src_grid_dims(2)=full_Jsize(ia)
          dst_grid_dims(1)=full_Isize(parent_id(ia))*par_grid_ratio(ia)
          dst_grid_dims(2)=full_Jsize(parent_id(ia))*par_grid_ratio(ia)
          num_sparse_elems=src_grid_dims(1)*src_grid_dims(2)
          allocate ( sparse_rows(num_sparse_elems) )
!
          Istr=(grid%i_parent_start-1)*par_grid_ratio(ia)+1
          Jstr=(grid%j_parent_start-1)*par_grid_ratio(ia)+1
          jc=0
          DO j=1,src_grid_dims(2)
            count=(Jstr-1+j-1)*(dst_grid_dims(1))+Istr-1
            DO i=1,src_grid_dims(1)
              jc=jc+1
              sparse_rows(jc)=count+i
            END DO
          END DO
          call SparseMatrix_importGRowInd(sMatF, sparse_rows,           &
     &                                    num_sparse_elems)
          deallocate ( sparse_rows )
!
! Second the coarse to fine.
!
          dst_grid_dims(1)=full_Isize(ia)
          dst_grid_dims(2)=full_Jsize(ia)
          src_grid_dims(1)=full_Isize(parent_id(ia))*par_grid_ratio(ia)
          src_grid_dims(2)=full_Jsize(parent_id(ia))*par_grid_ratio(ia)
          num_sparse_elems=dst_grid_dims(1)*dst_grid_dims(2)
          allocate ( sparse_cols(num_sparse_elems) )
          Istr=(grid%i_parent_start-1)*par_grid_ratio(ia)+1
          Jstr=(grid%j_parent_start-1)*par_grid_ratio(ia)+1
          jc=0
          DO j=1,dst_grid_dims(2)
            count=(Jstr-1+j-1)*(src_grid_dims(1))+Istr-1
            DO i=1,dst_grid_dims(1)
              jc=jc+1
              sparse_cols(jc)=count+i
            END DO
          END DO
          call SparseMatrix_importGColInd(sMatC, sparse_cols,           &
     &                                    num_sparse_elems)
          deallocate ( sparse_cols )
        END IF
        IF (moving_nest(ia).eq.1) THEN
!  reinit the matrices
          call SparseMatrixPlus_clean(SMPlus_R(ia)%F2CMatPlus)
          call SparseMatrixPlus_init(SMPlus_R(ia)%F2CMatPlus, sMatF,    &
     &                               GlobalSegMap_G(ia)%GSMapWRF,       &
     &                               GSMapInterp_R(ia)%GSMapStatic,     &
     &                               Xonly, 0, ATM_COMM_WORLD, ATMid)
!
          call SparseMatrixPlus_clean(SMPlus_R(ia)%C2FMatPlus)
          call SparseMatrixPlus_init(SMPlus_R(ia)%C2FMatPlus, sMatC,    &
     &                               GSMapInterp_R(ia)%GSMapStatic,     &
     &                               GlobalSegMap_G(ia)%GSMapWRF,       &
     &                               Xonly, 0, ATM_COMM_WORLD, ATMid)
        END IF
#  endif
!
!-----------------------------------------------------------------------
!  Export fields from atmosphere (WRF) to ocean (ROMS) model.
!-----------------------------------------------------------------------
!     GSW        Short wave raditaion (W m-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%GSW(i,j)
#  ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#  endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "GSW",   &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     GLW        Long wave raditaion (W m-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%GLW(i,j)
#  ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#  endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "GLW",   &
     &                           AA, Asize)
#  ifdef ATM2OCN_FLUXES
!-----------------------------------------------------------------------
!     LH     Latent heat flux (W m-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%LH(i,j)
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "LH",    &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     HFX     Sensible heat flux (W m-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%HFX(i,j)
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "HFX",   &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     USTRESS    Surface u-stress (m2 s-2).
!     Use Uearth and Vearth rotated velocities to partition UST.
!
      ij=0
      DO j=js,je
        DO i=is,ie
          cff1=1.0/(grid%alt(i,1,j)+eps)
          cff2=2.0/(((grid%u_2(i,1,j)+grid%u_2(i+1,1,j))**2+            &
     &               (grid%v_2(i,1,j)+grid%v_2(i,1,j+1))**2)**0.5+eps)
          ij=ij+1
          cff3=0.5*(grid%u_2(i,1,j)+grid%u_2(i+1,1,j))*grid%cosa(i,j)-  &
     &         0.5*(grid%v_2(i,1,j)+grid%v_2(i,1,j+1))*grid%sina(i,j)
          cff=cff1*cff2*(grid%UST(i,j)**2)*cff3
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV,"USTRESS",&
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     VSTRESS    Surface v-stress (m2 s-2).
!     Use Uearth and Vearth rotated velocities to partition UST.
!
      ij=0
      DO j=js,je
        DO i=is,ie
          cff1=1.0/(grid%alt(i,1,j)+eps)
          cff2=2.0/(((grid%u_2(i,1,j)+grid%u_2(i+1,1,j))**2+            &
     &               (grid%v_2(i,1,j)+grid%v_2(i,1,j+1))**2)**0.5+eps)
          ij=ij+1
          cff3=0.5*(grid%v_2(i,1,j)+grid%v_2(i,1,j+1))*grid%cosa(i,j)+  &
     &         0.5*(grid%u_2(i,1,j)+grid%u_2(i+1,1,j))*grid%sina(i,j)
          cff=cff1*cff2*(grid%UST(i,j)**2)*cff3
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV,"VSTRESS",&
     &                           AA, Asize)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!-----------------------------------------------------------------------
!     MSLP       Surface atmospheric pressure (Pa).
!     Use the hypsometric equation to reduce 
!     surface pressure to mean sea level pressure.
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%PSFC(i,j)*                                           &
     &           exp((9.81*grid%ht(i,j))/                               &
     &           (287.0*grid%T2(i,j)*(1.0+0.61*grid%Q2(i,j))))
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "MSLP",  &
     &                           AA, Asize)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || \
     (defined SHORTWAVE && defined ANA_SRFLUX)
!-----------------------------------------------------------------------
!     RELH       Surface air relative humidity (-).
!
      ij=0
      DO j=js,je
        DO i=is,ie
!
!         Calculate 2-m pressure using hypsometric equation. 
!         Assume temp at 2m = temp at 0m.
!
          cff1 = grid%PSFC(i,j)/(exp((9.81*2.0)/(287.0*grid%T2(i,j))))
!
!         Compute specific humidity using the 2-m mixing ratio and 
!         2-m pressure.
!
          rnum = grid%Q2(i,j)*cff1
          rden  = (grid%Q2(i,j)*(1.-0.622)+0.622)
          cff2 = rnum/rden
!
!         Compute saturation specific humidity using Bolton equation 10.
!
          c04 = 17.67*(grid%T2(i,j)-273.15)
          c05 = (grid%T2(i,j)-273.15) + 243.5
          cff3  = 6.112*exp(c04/c05)
!
          ij=ij+1
          cff=cff2/cff3
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "RELH",  &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     T2         Surface 2m air temperature (Convert to C).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%T2(i,j)-273.15
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "T2",    &
     &                           AA, Asize)
#  endif
!-----------------------------------------------------------------------
#  if defined BULK_FLUXES || defined ECOSIM
!     U10        U-Wind speed at 10 m (m s-1).
!     Rotate to earth lon lat.
!     
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%U10(i,j)*grid%cosa(i,j)-grid%V10(i,j)*grid%sina(i,j)
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "U10",   &
     &                           AA, Asize)
!
!-----------------------------------------------------------------------
!     V10        V-Wind speed at 10 m (m s-1).
!     Rotate to earth lon lat.
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%V10(i,j)*grid%cosa(i,j)+grid%U10(i,j)*grid%sina(i,j)
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "V10",   &
     &                           AA, Asize) 
#  endif
#  ifdef CLOUDS
!-----------------------------------------------------------------------
!     CLDFRA     Cloud fraction (--, 0-1.0).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%CLDFRA(i,1,j)
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "CLDFRA",&
     &                           AA, Asize)
#  endif
#  if !defined ANA_RAIN && defined EMINUSP
!-----------------------------------------------------------------------
!     RAIN       Precipitation (Convert to m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=0.001*(grid%RAINCV(i,j)+grid%RAINNCV(i,j))/grid%dt
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "RAIN",  &
     &                           AA, Asize)
#  endif
#  if defined EMINUSP
!-----------------------------------------------------------------------
!     EVAP      Evaporation (kg/m2 to Convert to m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
!         AA(ij)=0.001*(grid%SFCEVP(i,j))/grid%dt
          cff=0.001*(grid%QFX(i,j))
#   ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#   endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "EVAP",  &
     &                           AA, Asize)
#  endif
!-----------------------------------------------------------------------
!  Send fields to ocean model.
!
      Tag=io*100+ia*10+0
#  ifdef MOVE_NESTS
      IF (moving_nest(ia).eq.1) THEN
        CALL MCT_MatVecMul(AttrVect_O(ia,io)%atm2ocn_AV,                &
     &                     SMPlus_R(ia)%F2CMatPlus,                     &
     &                     AV2_RO(ia,io)%fin2crs_AV2)
!
        CALL MCT_isend (AV2_RO(ia,io)%fin2crs_AV2,                      &
     &                  Router_O(ia,io)%WRFtoROMS, Tag)
      ELSE
#  endif
        CALL MCT_isend (AttrVect_O(ia,io)%atm2ocn_AV,                   &
     &                  Router_O(ia,io)%WRFtoROMS, Tag)
#  ifdef MOVE_NESTS
      END IF
#  endif
!
      CALL MCT_waits (Router_O(ia,io)%WRFtoROMS)
      IF (MYRANK.EQ.0) THEN
        WRITE (*,36) ' ## WRF grid ',ia,                                &
     &                    ' sent data to ROMS grid ',io
 36     FORMAT (a14,i2,a24,i2)
      ENDIF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupling send fail atm_coupler, error= ', MyError
        CALL finalize_atm_coupling
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (AA)
#  ifdef MOVE_NESTS
      deallocate (MOV_MASK)
#  endif
      RETURN
      END SUBROUTINE atm2ocn_coupling
# endif

# ifdef SWAN_COUPLING
      SUBROUTINE atm2wav_coupling (grid, ia, iw)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupled data streams.          !
!  Currently, the following data streams are processed:                !
!                                                                      !
!  Possible fields exported to the OCEAN Model:                        !
!                                                                      !
!     * GSW        Short wave raditaion  (Watts/m2)                    !
!     * GLW        Long wave raditaion  (Watts/m2)                     !
!     * LH         Latent heat flux     (Watts/m2)                     !
!     * HFX        Sensible heat flux   (Watts/m2)                     !
!     * USTRESS    Surface U-wind stress (Pa)                          !
!     * VSTRESS    Surface v-stress      (Pa)                          !
!     * MSLP       Mean Sea Level Pressure (Pa)                        !
!     * RELH       Surface air relative humidity (percent)             !
!     * T2         Surface 2m air temperature (Celsius)                !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!     * CLDFRA     Cloud fraction       (percent/100)                  !
!     * RAIN       Precipitation        (m/s)                          !
!     * EVAP       Evaporation          (m/s)                          !
!                                                                      !
!  Fields exported to the WAVE Model:                                  !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!                                                                      !
!  Fields acquired from the WAVE Model:                                !
!                                                                      !
!     * HISGN      Significant wave heigth (m)                         !
!     * WLENP      Peak wave length (m)                                !
!     * RTP        Peak wave period (s)                                !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * SST        Sea surface temperature                             !
!=======================================================================
!
      USE module_domain
!
      implicit none
!     TYPE(domain) , INTENT (IN) :: grid 
      TYPE(domain) , POINTER     :: grid            ! jcw check this
#  ifdef MOVE_NESTS
      TYPE(domain) , POINTER     :: grid_ptr
#  endif
      integer, intent(in) :: ia, iw
!
!  Local variable declarations.
!
      integer :: is, ie, js, je, ij, Tag
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, nprocs, offset
#  ifdef MOVE_NESTS
      integer :: count, jc, num_sparse_elems
#  endif

      real, parameter :: eps=1.0e-10
      real :: cff
      real, pointer :: AA(:)
#  ifdef MOVE_NESTS
      real, pointer :: MOV_MASK(:,:)
#  endif
!
!  Set grid range.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
!
!-----------------------------------------------------------------------
!  Send atmosphere fields to SWAN.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (ATM_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (ATM_COMM_WORLD, nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ia)%GSMapWRF,             &
     &                         ATM_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( AA(Asize), stat=ierr )
!
#  ifdef MOVE_NESTS
!  If we have a moving nest, then need to mask out parent location of
!  the child nest.
!
      allocate ( MOV_MASK(is:ie,js:je), stat=ierr )
      DO i=is,ie
        DO j=js,je
          MOV_MASK(i,j)=1.0
        END DO
      END DO
!  Get Istr and Jstr of moving child in parent grid.
!  Allow parent to provide one more row around child perimeter.


!  want to know:
!  1 - do you have a child that is moving?
!  2  - if so, then where is it now?

        IF (ia.eq.parent_id(Natm_grids)) THEN
          CALL find_grid_by_id(Natm_grids, grid, grid_ptr)
          offset=1
          Istr=grid_ptr%i_parent_start+offset
          Jstr=grid_ptr%j_parent_start+offset
          ij=full_Isize(Natm_grids)/par_grid_ratio(Natm_grids)
          Iend=grid_ptr%i_parent_start+ij-1-offset
          ij=full_Jsize(Natm_grids)/par_grid_ratio(Natm_grids)
          Jend=grid_ptr%j_parent_start+ij-1-offset
!
          IF ((Istr.ge.is).or.(Istr.le.ie)) THEN
            IF ((Jstr.ge.js).or.(Jstr.le.je)) THEN
              DO i=MAX(is,Istr),MIN(Iend,ie)
                DO j=MAX(js,Jstr),MIN(Jend,je)
                  MOV_MASK(i,j)=0.0
                END DO
              END DO
            END IF
          END IF
        END IF
!
!  Adjust the sparse matrices here.
!
!       IF (ia.eq.Natm_grids) THEN
        IF ((Myrank.eq.0).and.(moving_nest(ia).eq.1)) THEN

!  First the fine to coarse
          src_grid_dims(1)=full_Isize(ia)
          src_grid_dims(2)=full_Jsize(ia)
          dst_grid_dims(1)=full_Isize(parent_id(ia))*par_grid_ratio(ia)
          dst_grid_dims(2)=full_Jsize(parent_id(ia))*par_grid_ratio(ia)
          num_sparse_elems=src_grid_dims(1)*src_grid_dims(2)
          allocate ( sparse_rows(num_sparse_elems) )
!
          Istr=(grid%i_parent_start-1)*par_grid_ratio(ia)+1
          Jstr=(grid%j_parent_start-1)*par_grid_ratio(ia)+1
          jc=0
          DO j=1,src_grid_dims(2)
            count=(Jstr-1+j-1)*(dst_grid_dims(1))+Istr-1
            DO i=1,src_grid_dims(1)
              jc=jc+1
              sparse_rows(jc)=count+i
            END DO
          END DO
          call SparseMatrix_importGRowInd(sMatF, sparse_rows,           &
     &                                    num_sparse_elems)
          deallocate ( sparse_rows )
!
! Second the coarse to fine.
!
          dst_grid_dims(1)=full_Isize(ia)
          dst_grid_dims(2)=full_Jsize(ia)
          src_grid_dims(1)=full_Isize(parent_id(ia))*par_grid_ratio(ia)
          src_grid_dims(2)=full_Jsize(parent_id(ia))*par_grid_ratio(ia)
          num_sparse_elems=dst_grid_dims(1)*dst_grid_dims(2)
          allocate ( sparse_cols(num_sparse_elems) )
          Istr=(grid%i_parent_start-1)*par_grid_ratio(ia)+1
          Jstr=(grid%j_parent_start-1)*par_grid_ratio(ia)+1
          jc=0
          DO j=1,dst_grid_dims(2)
            count=(Jstr-1+j-1)*(src_grid_dims(1))+Istr-1
            DO i=1,dst_grid_dims(1)
              jc=jc+1
              sparse_cols(jc)=count+i
            END DO
          END DO
          call SparseMatrix_importGColInd(sMatC, sparse_cols,           &
     &                                    num_sparse_elems)
          deallocate ( sparse_cols )
        END IF
        IF (moving_nest(ia).eq.1) THEN
!  reinit the matrices
          call SparseMatrixPlus_clean(SMPlus_R(ia)%F2CMatPlus)
          call SparseMatrixPlus_init(SMPlus_R(ia)%F2CMatPlus, sMatF,    &
     &                               GlobalSegMap_G(ia)%GSMapWRF,       &
     &                               GSMapInterp_R(ia)%GSMapStatic,     &
     &                               Xonly, 0, ATM_COMM_WORLD, ATMid)
!
          call SparseMatrixPlus_clean(SMPlus_R(ia)%C2FMatPlus)
          call SparseMatrixPlus_init(SMPlus_R(ia)%C2FMatPlus, sMatC,    &
     &                               GSMapInterp_R(ia)%GSMapStatic,     &
     &                               GlobalSegMap_G(ia)%GSMapWRF,       &
     &                               Xonly, 0, ATM_COMM_WORLD, ATMid)
        END IF
#  endif
!-----------------------------------------------------------------------
!  Send U10 and V10 to SWAN.
!-----------------------------------------------------------------------
!     U10        U-Wind speed at 10 m (m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%U10(i,j)*grid%cosa(i,j)-grid%V10(i,j)*grid%sina(i,j)
#  ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#  endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_W(ia,iw)%atm2wav_AV, "U10",   &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     V10        V-Wind speed at 10 m (m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=grid%V10(i,j)*grid%cosa(i,j)+grid%U10(i,j)*grid%sina(i,j)
#  ifdef MOVE_NESTS
          cff=cff*MOV_MASK(i,j)
#  endif
          AA(ij)=cff
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_W(ia,iw)%atm2wav_AV, "V10",   &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!
      Tag=0+ia*10+iw
#  ifdef MOVE_NESTS
      IF (moving_nest(ia).eq.1) THEN
        CALL MCT_MatVecMul(AttrVect_W(ia,iw)%atm2wav_AV,                &
     &                     SMPlus_R(ia)%F2CMatPlus,                     &
     &                     AV2_RW(ia,iw)%fin2crs_AV2)
!
        CALL MCT_isend (AV2_RW(ia,iw)%fin2crs_AV2,                      &
     &                  Router_W(ia,iw)%WRFtoSWAN, Tag)
      ELSE
#  endif
        CALL MCT_isend (AttrVect_W(ia,iw)%atm2wav_AV,                   &
     &                  Router_W(ia,iw)%WRFtoSWAN, Tag)
#  ifdef MOVE_NESTS
      END IF
#  endif
!
      CALL MCT_waits (Router_W(ia,iw)%WRFtoSWAN)
      IF (MYRANK.EQ.0) THEN
        WRITE (*,37) ' ## WRF grid ',ia,                                &
     &                    ' sent data to SWAN grid ',iw
 37     FORMAT (a14,i2,a24,i2)
      ENDIF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupling send fail atm_coupler, error= ', MyError
        CALL finalize_atm_coupling
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (AA)
      RETURN
      END SUBROUTINE atm2wav_coupling
# endif

# ifdef ROMS_COUPLING
      SUBROUTINE atmfocn_coupling (grid, ia, io, aux4flag)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupled data streams.          !
!  Currently, the following data streams are processed:                !
!                                                                      !
!  Possible fields exported to the OCEAN Model:                        !
!                                                                      !
!     * GSW        Short wave raditaion  (Watts/m2)                    !
!     * GLW        Long wave raditaion  (Watts/m2)                     !
!     * LH         Latent heat flux     (Watts/m2)                     !
!     * HFX        Sensible heat flux   (Watts/m2)                     !
!     * USTRESS    Surface U-wind stress (Pa)                          !
!     * VSTRESS    Surface v-stress      (Pa)                          !
!     * MSLP       Mean Sea Level Pressure (Pa)                        !
!     * RELH       Surface air relative humidity (percent)             !
!     * T2         Surface 2m air temperature (Celsius)                !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!     * CLDFRA     Cloud fraction       (percent/100)                  !
!     * RAIN       Precipitation        (m/s)                          !
!     * EVAP       Evaporation          (m/s)                          !
!                                                                      !
!  Fields exported to the WAVE Model:                                  !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!                                                                      !
!  Fields acquired from the WAVE Model:                                !
!                                                                      !
!     * HISGN      Significant wave heigth (m)                         !
!     * WLENP      Peak wave length (m)                                !
!     * RTP        Peak wave period (s)                                !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * SST        Sea surface temperature                             !
!=======================================================================
!
      USE module_domain
!
      implicit none
      TYPE(domain) , INTENT (IN) :: grid 
      integer, intent(in) :: ia, io, aux4flag
!
!  Local variable declarations.
!
      integer :: is, ie, js, je, ij, Tag
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, nprocs

      real, parameter :: eps=1.0e-10
      real, parameter ::  Large = 1.0E+20
      real :: cff, inval, retval
      real, pointer :: AA(:)
      real, pointer :: Amask(:)
      real, dimension(2) :: range
!  Set grid range.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
!
!-----------------------------------------------------------------------
!  Send atmosphere fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (ATM_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (ATM_COMM_WORLD, nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ia)%GSMapWRF,             &
     &                         ATM_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( AA(Asize), stat=ierr )
      allocate ( Amask(Asize), stat=ierr )
!
!  Schedule receiving fields from ocean model.
!  For now, use this aux4flag to do update if sst_update is used.
!
      IF (aux4flag.eq.0) THEN
        Tag=io*100+ia*10+0
#  ifdef MOVE_NESTS
        IF (moving_nest(ia).eq.1) THEN
          CALL MCT_irecv (AV2_RO(ia,io)%crs2fin_AV2,                    &
     &                    Router_O(ia,io)%WRFtoROMS, Tag)
! Wait to make sure the WRF data has arrived.
          CALL MCT_waitr (AV2_RO(ia,io)%crs2fin_AV2,                    &
     &                    Router_O(ia,io)%WRFtoROMS)
          CALL MCT_MatVecMul(AV2_RO(ia,io)%crs2fin_AV2,                 &
     &                       SMPlus_R(ia)%C2FMatPlus,                   &
     &                       AttrVect_O(ia,io)%ocn2atm_AV)
        ELSE
#  endif
          CALL MCT_irecv (AttrVect_O(ia,io)%ocn2atm_AV,                 &
     &                    Router_O(ia,io)%WRFtoROMS,Tag)
! Wait to make sure the OCN data has arrived.
          CALL MCT_waitr (AttrVect_O(ia,io)%ocn2atm_AV,                 &
     &                    Router_O(ia,io)%WRFtoROMS)
#  ifdef MOVE_NESTS
        END IF
#  endif
!
        IF (MYRANK.EQ.0) THEN
          WRITE (*,38) ' ## WRF grid ',ia,                              &
     &                      ' recvd data from ROMS grid ',io
 38       FORMAT (a13,i2,a27,i2)
        END IF
        IF (MyError.ne.0) THEN
          WRITE (*,*) 'coupling fail wrfcplr, MyStatus= ', MyError
          CALL finalize_atm_coupling
        END IF
      END IF
!
 40         FORMAT (a36,1x,2(1pe14.6))
!  SST (Convert to K).
!  USE CPLMASK:  coupling mask (0 for data read in wrflowinput,
!                               1 data received from the coupler)
#  ifdef MCT_INTERP_OC2AT
      CALL AttrVect_exportRAttr (AttrVect_O(ia,io)%ocn2atm_AV,          &
     &                           "CPL_MASK", AA, Asize)
#  endif
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
#  ifdef SST_CONST
          Amask(ij)=0.
#  else
#   ifdef MCT_INTERP_OC2AT
          Amask(ij)=AA(ij)
#   else
          Amask(ij)=1.
#   endif
#  endif
        END DO
      END DO
!
      CALL AttrVect_exportRAttr (AttrVect_O(ia,io)%ocn2atm_AV,          &
     &                           "SST", AA, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=(AA(ij)+273.15)*Amask(ij)
          grid%sst(i,j)=cff+grid%sst(i,j)*(1.0-Amask(ij))
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      IF (aux4flag.eq.0) THEN
        CALL mpi_allreduce ( range(1), retval , 1, MPI_REAL,            &
     &                       MPI_MIN, ATM_COMM_WORLD, ierr )
        range(1)=retval
        CALL mpi_allreduce ( range(2), retval , 1, MPI_REAL,            &
     &                       MPI_MAX, ATM_COMM_WORLD, ierr )
        range(2)=retval
        IF (Myrank.eq.0) THEN
          write(*,40) 'ROMStoWRF  Min/Max SST     (K):     ',           &
     &                      range(1),range(2)
        END IF
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (AA)
      deallocate (Amask)
      RETURN
      END SUBROUTINE atmfocn_coupling
# endif

# ifdef SWAN_COUPLING
      SUBROUTINE atmfwav_coupling (grid, ia, iw)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupled data streams.          !
!  Currently, the following data streams are processed:                !
!                                                                      !
!  Possible fields exported to the OCEAN Model:                        !
!                                                                      !
!     * GSW        Short wave raditaion  (Watts/m2)                    !
!     * GLW        Long wave raditaion  (Watts/m2)                     !
!     * LH         Latent heat flux     (Watts/m2)                     !
!     * HFX        Sensible heat flux   (Watts/m2)                     !
!     * USTRESS    Surface U-wind stress (Pa)                          !
!     * VSTRESS    Surface v-stress      (Pa)                          !
!     * MSLP       Mean Sea Level Pressure (Pa)                        !
!     * RELH       Surface air relative humidity (percent)             !
!     * T2         Surface 2m air temperature (Celsius)                !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!     * CLDFRA     Cloud fraction       (percent/100)                  !
!     * RAIN       Precipitation        (m/s)                          !
!     * EVAP       Evaporation          (m/s)                          !
!                                                                      !
!  Fields exported to the WAVE Model:                                  !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!                                                                      !
!  Fields acquired from the WAVE Model:                                !
!                                                                      !
!     * HISGN      Significant wave heigth (m)                         !
!     * WLENP      Peak wave length (m)                                !
!     * RTP        Peak wave period (s)                                !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * SST        Sea surface temperature                             !
!=======================================================================
!
      USE module_domain
!
      implicit none
      TYPE(domain) , INTENT (IN) :: grid 
      integer, intent(in) :: ia, iw
!
!  Local variable declarations.
!
      integer :: is, ie, js, je, ij, Tag
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, nprocs

      real, parameter :: eps=1.0e-10
      real, pointer :: AA(:)
      real, parameter ::  Large = 1.0E+20
      real :: cff, inval, retval
      real, dimension(2) :: range
#  ifdef MCT_INTERP_WV2AT
      real, pointer :: Amask(:)
#  endif
!
!  Set grid range.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
!
!-----------------------------------------------------------------------
!  Send atmosphere fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (ATM_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (ATM_COMM_WORLD, nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ia)%GSMapWRF,             &
     &                         ATM_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( AA(Asize), stat=ierr )
#  ifdef MCT_INTERP_WV2AT
      allocate ( Amask(Asize), stat=ierr )
#  endif
!
!  Schedule receiving fields from wave model.
!
      Tag=0+ia*10+iw
#  ifdef MOVE_NESTS
      IF (moving_nest(ia).eq.1) THEN
        CALL MCT_irecv (AV2_RW(ia,iw)%crs2fin_AV2,                      &
     &                  Router_W(ia,iw)%WRFtoSWAN, Tag)
!  Wait to make sure the WRF data has arrived.
        CALL MCT_waitr (AV2_RW(ia,iw)%crs2fin_AV2,                      &
     &                  Router_W(ia,iw)%WRFtoSWAN)
        CALL MCT_MatVecMul(AV2_RW(ia,iw)%crs2fin_AV2,                   &
     &                     SMPlus_R(ia)%C2FMatPlus,                     &
     &                     AttrVect_W(ia,iw)%wav2atm_AV)
      ELSE
#  endif
        CALL MCT_irecv (AttrVect_W(ia,iw)%wav2atm_AV,                   &
     &                  Router_W(ia,iw)%WRFtoSWAN, Tag)
!  Wait to make sure the WAV data has arrived.
        CALL MCT_waitr (AttrVect_W(ia,iw)%wav2atm_AV,                   &
     &                   Router_W(ia,iw)%WRFtoSWAN)
#  ifdef MOVE_NESTS
      END IF
#  endif
!
      IF (MYRANK.EQ.0) THEN
        WRITE (*,30) ' ## WRF grid ',ia,                                &
     &                    ' recv data from SWAN grid ',iw
 30     FORMAT (a14,i2,a26,i2)
      END IF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupl fail wrfcplr, MyStatus= ', MyError
        CALL finalize_atm_coupling
      END IF
 40         FORMAT (a36,1x,2(1pe14.6))
!
!  CPL_MASK
!
#  ifdef MCT_INTERP_WV2AT
      CALL AttrVect_exportRAttr (AttrVect_W(ia,iw)%wav2atm_AV,          &
     &                           "CPL_MASK", AA, Asize)
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          Amask(ij)=AA(ij)
        END DO
      END DO
#  endif
!
!  Hwave.
!
      CALL AttrVect_exportRAttr (AttrVect_W(ia,iw)%wav2atm_AV, "HSIGN", &
     &                           AA, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=MAX(0.0,AA(ij))
#  ifdef MCT_INTERP_WV2AT
          grid%hwave(i,j)=cff*Amask(ij)+                                &
     &                    grid%hwave(i,j)*(1.0-Amask(ij))
#  else
          grid%hwave(i,j)=cff
#  endif
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mpi_allreduce ( range(1), retval , 1, MPI_REAL,              &
     &                     MPI_MIN, ATM_COMM_WORLD, ierr )
      range(1)=retval
      CALL mpi_allreduce ( range(2), retval , 1, MPI_REAL,              &
     &                     MPI_MAX, ATM_COMM_WORLD, ierr )
      range(2)=retval
      IF (Myrank.eq.0) THEN
        write(*,40) 'SWANtoWRF  Min/Max HSIGN   (m):     ',             &
     &                    range(1),range(2)
      END IF
!
!  Lwave.
!
      CALL AttrVect_exportRAttr (AttrVect_W(ia,iw)%wav2atm_AV, "WLENP", &
     &                           AA, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=MAX(0.0,AA(ij))
#  ifdef MCT_INTERP_WV2AT
          grid%lwavep(i,j)=cff*Amask(ij)+                   &
     &                     grid%lwavep(i,j)*(1.0-Amask(ij))
#  else
          grid%lwavep(i,j)=cff
#  endif
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mpi_allreduce ( range(1), retval , 1, MPI_REAL,              &
     &                     MPI_MIN, ATM_COMM_WORLD, ierr )
      range(1)=retval
      CALL mpi_allreduce ( range(2), retval , 1, MPI_REAL,              &
     &                     MPI_MAX, ATM_COMM_WORLD, ierr )
      range(2)=retval
      IF (Myrank.eq.0) THEN
        write(*,40) 'SWANtoWRF  Min/Max WLENP   (m):     ',             &
     &                    range(1),range(2)
      END IF
!
!  Pwave_top.
!
      CALL AttrVect_exportRAttr (AttrVect_W(ia,iw)%wav2atm_AV, "RTP",   &
     &                           AA, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          cff=MAX(0.0,AA(ij))
#  ifdef MCT_INTERP_WV2AT
          grid%pwave(i,j)=cff*Amask(ij)+                    &
     &                     grid%pwave(i,j)*(1.0-Amask(ij))
#  else
          grid%pwave(i,j)=cff
#  endif
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
      CALL mpi_allreduce ( range(1), retval , 1, MPI_REAL,              &
     &                     MPI_MIN, ATM_COMM_WORLD, ierr )
      range(1)=retval
      CALL mpi_allreduce ( range(2), retval , 1, MPI_REAL,              &
     &                     MPI_MAX, ATM_COMM_WORLD, ierr )
      range(2)=retval
      IF (Myrank.eq.0) THEN
        write(*,40) 'SWANtoWRF  Min/Max RTP     (m):     ',             &
     &                    range(1),range(2)
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (AA)
#  ifdef MCT_INTERP_WV2AT
      deallocate (Amask)
#  endif
      RETURN
      END SUBROUTINE atmfwav_coupling
#endif

      SUBROUTINE finalize_atm_coupling
!
!=======================================================================
!                                                                    ===
!  This routines terminates execution during coupling error.         ===
!                                                                    ===
!=======================================================================
!
      implicit none
!
!  Local variable declarations.
!
      integer :: ia, io, iw, MyStatus
!
!-----------------------------------------------------------------------
!  Terminate MPI execution environment.
!-----------------------------------------------------------------------
!
!     CALL Router_clean (RoutWRFtoROMS)
      DO ia=1,Natm_grids
        CALL GlobalSegMap_clean (GlobalSegMap_G(ia)%GSMapWRF)
# ifdef ROMS_COUPLING
        DO io=1,Nocn_grids
          CALL AttrVect_clean (AttrVect_O(ia,io)%atm2ocn_AV)
          CALL AttrVect_clean (AttrVect_O(ia,io)%ocn2atm_AV)
        END DO
# endif
# ifdef SWAN_COUPLING
        DO iw=1,Nwav_grids
          CALL AttrVect_clean (AttrVect_W(ia,iw)%atm2wav_AV)
          CALL AttrVect_clean (AttrVect_W(ia,iw)%wav2atm_AV)
        END DO
# endif
      END DO
!     CALL MCTWorld_clean ()

      END SUBROUTINE finalize_atm_coupling
#endif
      END MODULE atm_coupler_mod
