/*
** svn $Id: mct_roms_hydro.h 830 2017-01-24 21:21:11Z jwarner $
***************************************************** John C. Warner ***
** Copyright (c) 2002-2020 The ROMS/TOMS Group      Hernan G. Arango  **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** These routines are use couple WRF Hydro to ROMS ocean model        **
** using the Model Coupling Toolkit (MCT).                            **
**                                                                    **
************************************************************************
*/

      SUBROUTINE initialize_ocn2hyd_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  Initialize ocean and atmosphere models coupling stream. This is     !
!  the training phase used to constuct  MCT parallel interpolators     !
!  and stablish communication patterns.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mct_coupler_params
      USE mod_kinds
      USE mod_scalars
      USE mod_iounits
#ifdef MCT_INTERP_OC2HY
      USE mod_coupler_iounits
#endif
!
!  Imported variable definitions.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrT, IendT, JstrT, JendT
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, Isize, Jsize, MyError
      integer :: i, ic, ih, j, jc, nprocs, cid, cad
      integer :: nRows, nCols, num_sparse_elems

      real(r8) :: cff

      integer, allocatable  :: length(:)
      integer, allocatable  :: start(:)
      character (len=70) :: nc_name
      character (len=20) :: to_add
      character (len=120) :: hostring
      character (len=120) :: ohstring
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IstrR=BOUNDS(ng)%Istr(tile)-1
      ELSE
        IstrR=BOUNDS(ng)%Istr(tile)
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IendR=BOUNDS(ng)%Iend(tile)+1
      ELSE
        IendR=BOUNDS(ng)%Iend(tile)
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        JstrR=BOUNDS(ng)%Jstr(tile)-1
      ELSE
        JstrR=BOUNDS(ng)%Jstr(tile)
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        JendR=BOUNDS(ng)%Jend(tile)+1
      ELSE
        JendR=BOUNDS(ng)%Jend(tile)
      END IF
!
!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, nprocs, MyError)
!
#ifdef MCT_INTERP_OC2HY
      IF (ng.eq.1) THEN
        ALLOCATE(SMPlus_H(Nocn_grids,Nhyd_grids))
        ALLOCATE(AV2_H(Nocn_grids,Nhyd_grids))
        ALLOCATE(GSMapInterp_H(Nocn_grids,Nhyd_grids))
      END IF
#endif
      OCNid=ocnids(ng)
#if !defined WAVES_OCEAN
      IF (ng.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Nocn_grids))
        ALLOCATE(AttrVect_G(Nocn_grids))
      END IF
!
!  Initialize MCT coupled model registry.
!
      IF (Nocn_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      OCN_COMM_WORLD,myids=ocnids)
      ELSE
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      OCN_COMM_WORLD,OCNid)
      END IF
!
!  Determine start and lengths for domain decomposition.
!
      Jsize=JendR-JstrR+1
      IF (.not.allocated(start)) THEN
        allocate ( start(Jsize) )
      END IF
      IF (.not.allocated(length)) THEN
        allocate ( length(Jsize) )
      END IF
      jc=0
      DO j=JstrR,JendR
        jc=jc+1
        start (jc)=j*(Lm(ng)+2)+IstrR+1
        length(jc)=(IendR-IstrR+1)
      END DO
      CALL GlobalSegMap_init (GlobalSegMap_G(ng)%GSMapROMS,             &
     &                        start, length, 0, OCN_COMM_WORLD, OCNid)
!
!  Deallocate working arrays.
!
      IF (allocated(start)) THEN
        deallocate (start)
      END IF
      IF (allocated(length)) THEN
        deallocate (length)
      END IF
#endif
#ifdef MCT_INTERP_OC2HY
!
!  If ocean grid and atm grids are different sizes, then
!  develop sparse matrices for interpolation.
!
  35  FORMAT(a3,i1,a7,i1,a11)
      DO ih=1,Nhyd_grids
!
!!!!!!!!!!!!!!!!!!!!!!
! First work on hyd to ocean.
!!!!!!!!!!!!!!!!!!!!!!
!
        IF (Myrank.eq.MyMaster) THEN
          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'hyd',ih,'_to_ocn',ng,'_weights.nc'
          ELSE
            nc_name=H2Oname(ih,ng)
          END IF
          call get_sparse_matrix (ng, nc_name, num_sparse_elems,        &
     &                            src_grid_dims, dst_grid_dims)
!
! Init the sparse matrix.
!
          nRows=dst_grid_dims(1)*dst_grid_dims(2)
          nCols=src_grid_dims(1)*src_grid_dims(2)
!
! Create sparse matrix.
!
!         Sparse rows is the dst address. Multiply the interp weights
!         by the dst masking.
!
          DO i=1,num_sparse_elems
            j=sparse_rows(i)
            cff=REAL(dst_grid_imask(j),r8)
            sparse_weights(i)=sparse_weights(i)*cff
          END DO
!
! Load the dst grid as a coupling mask.
!
!          allocate(H2O_CPLMASK(ih,ng)%dst_mask(nRows))
!          DO i=1,nRows
!            H2O_CPLMASK(ih,ng)%dst_mask(i)=dst_grid_imask(i)
!          END DO


          call SparseMatrix_init(sMatH,nRows,nCols,num_sparse_elems)
          call SparseMatrix_importGRowInd(sMatH, sparse_rows,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importGColInd(sMatH, sparse_cols,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importMatrixElts(sMatH, sparse_weights,     &
     &                                    num_sparse_elems)
!
! Deallocate arrays.
!
          deallocate ( sparse_rows )
          deallocate ( sparse_cols )
          deallocate ( sparse_weights )
          deallocate ( dst_grid_imask )
        END IF
!
!
        CALL mpi_bcast(dst_grid_dims, 2, MPI_INTEGER, MyMaster,         &
     &                 OCN_COMM_WORLD, MyError)
!
! scatter dst_grid_imask to be used as cpl_mask
!
!        IF (Myrank.ne.MyMaster) THEN
!          nRows=dst_grid_dims(1)*dst_grid_dims(2)
!          allocate(H2O_CPLMASK(ih,ng)%dst_mask(nRows))
!        END IF
!        CALL mpi_bcast(H2O_CPLMASK(ih,ng)%dst_mask,nRows,               &
!     &                 MPI_INTEGER, MyMaster,                           &
!     &                 OCN_COMM_WORLD, MyError)


!!!!!!!!!!!!!!!!!!!!!!
! Second work on ocean to hyd.
!!!!!!!!!!!!!!!!!!!!!!
!
        IF (Myrank.eq.MyMaster) THEN
          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'ocn',ng,'_to_hyd',ih,'_weights.nc'
          ELSE
            nc_name=O2Hname(ng,ih)
          END IF
          call get_sparse_matrix (ng, nc_name, num_sparse_elems,        &
     &                            src_grid_dims, dst_grid_dims)
!
! Init the sparse matrix.
!
          nRows=dst_grid_dims(1)*dst_grid_dims(2)
          nCols=src_grid_dims(1)*src_grid_dims(2)
!
! Create sparse matrix.
!
          DO i=1,num_sparse_elems
            j=sparse_rows(i)
            cff=REAL(dst_grid_imask(j),r8)
            sparse_weights(i)=sparse_weights(i)*cff
          END DO
!
! Load the dst grid as a coupling mask.
!
!          allocate(O2H_CPLMASK(ng,ih)%dst_mask(nRows))
!          DO i=1,nRows
!            O2H_CPLMASK(ng,ih)%dst_mask(i)=dst_grid_imask(i)
!          END DO
!
          call SparseMatrix_init(sMatO,nRows,nCols,num_sparse_elems)
          call SparseMatrix_importGRowInd(sMatO, sparse_rows,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importGColInd(sMatO, sparse_cols,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importMatrixElts(sMatO, sparse_weights,     &
     &                                    num_sparse_elems)
!
! Deallocate arrays.
!
          deallocate ( sparse_rows )
          deallocate ( sparse_cols )
          deallocate ( sparse_weights )
          deallocate ( dst_grid_imask )
        END IF
!
        CALL mpi_bcast(dst_grid_dims, 2, MPI_INTEGER, MyMaster,         &
     &                 OCN_COMM_WORLD, MyError)
!
! scatter dst_grid_imask to be used as cpl_mask
!
!        IF (Myrank.ne.MyMaster) THEN
!          nRows=dst_grid_dims(1)*dst_grid_dims(2)
!          allocate(O2H_CPLMASK(ng,ih)%dst_mask(nRows))
!        END IF
!        CALL mpi_bcast(O2H_CPLMASK(ng,ih)%dst_mask,nRows,               &
!     &                 MPI_INTEGER, MyMaster,                           &
!     &                 OCN_COMM_WORLD, MyError)
!
!  Create Global Seg Map for hyd model.
!  Determine start and lengths for domain decomposition
!  of the atm model.
!
        Isize=INT(dst_grid_dims(1)/nprocs)
        IF (MyRank.eq.nprocs-1) THEN
          Isize=dst_grid_dims(1)-Isize*(nprocs-1)
        ENDIF
        IF (.not.allocated(start)) THEN
          allocate ( start(1) )
        END IF
        IF (.not.allocated(length)) THEN
          allocate ( length(1) )
        END IF
        start=(MyRank*INT(dst_grid_dims(1)/nprocs))*dst_grid_dims(2)+1
        length=Isize*dst_grid_dims(2)
!
        CALL GlobalSegMap_init (GSMapInterp_H(ng,ih)%GSMapHYD,          &
     &                          start, length, 0, OCN_COMM_WORLD, OCNid)
!
!  Deallocate working arrays.
!
        IF (allocated(start)) THEN
          deallocate (start)
        END IF
        IF (allocated(length)) THEN
          deallocate (length)
        END IF
!
! Create ATM sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
        call SparseMatrixPlus_init(SMPlus_H(ng,ih)%H2OMatPlus, sMatH,   &
     &                             GSMapInterp_H(ng,ih)%GSMapHYD,       &
     &                             GlobalSegMap_G(ng)%GSMapROMS,        &
     &                             Xonly,MyMaster,OCN_COMM_WORLD, OCNid)
        call SparseMatrix_clean(sMatH)
!
! Create Ocean sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
        call SparseMatrixPlus_init(SMPlus_H(ng,ih)%O2HMatPlus, sMatO,   &
     &                             GlobalSegMap_G(ng)%GSMapROMS,        &
     &                             GSMapInterp_H(ng,ih)%GSMapHYD,       &
     &                             Xonly,MyMaster,OCN_COMM_WORLD, OCNid)
        call SparseMatrix_clean(sMatO)
      END DO
#endif
!
!  Initialize attribute vector holding the data code strings from
!  the hydro model.
!
      cad=LEN(hostring)
      DO i=1,cad
        hostring(i:i)=''
      END DO
      cid=1
!
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
!  Initialize attribute vector holding the export data code strings of
!  the hydrology model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &                         OCN_COMM_WORLD)
      CALL AttrVect_init(AttrVect_G(ng)%hyd2ocn_AV,                     &
     &                   rList=TRIM(hostring),lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%hyd2ocn_AV)
!
!  Initialize attribute vector holding the export data of
!  the hyd model.
!
!  Initialize attribute vector that contain the data strings from
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
!#ifdef MCT_INTERP_OC2AT
!      to_add=':CPL_MASK'
!      cad=LEN_TRIM(to_add)
!      write(ohstring(cid:cid+cad-1),'(a)') to_add(1:cad)
!      cid=cid+cad
!#endif
!
!  Finalize and remove trailing spaces from the ohstring
!  for the rlist.
!
      cad=LEN_TRIM(ohstring)
      ohstring=ohstring(1:cad)
!
      CALL AttrVect_init (AttrVect_G(ng)%ocn2hyd_AV,                    &
     &                    rList=TRIM(ohstring),lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%ocn2hyd_AV)
!
#ifdef MCT_INTERP_OC2HY
      DO ih=1,Nhyd_grids
        HYDid=hydids(ih)
!
!  Initialize attribute vector holding the export data of
!  the hyd model. The Asize is the number of grid point on this
!  processor.
!
        Asize=GlobalSegMap_lsize(GSMapInterp_H(ng,ih)%GSMapHYD,        &
     &                           OCN_COMM_WORLD)
        CALL AttrVect_init (AV2_H(ng,ih)%hyd2ocn_AV2,                   &
     &                      rList=TRIM(hostring),lsize=Asize)
        CALL AttrVect_zero (AV2_H(ng,ih)%hyd2ocn_AV2)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
        CALL AttrVect_init (AV2_H(ng,ih)%ocn2hyd_AV2,                   &
     &                      rList=TRIM(ohstring),lsize=Asize)
        CALL AttrVect_zero (AV2_H(ng,ih)%ocn2hyd_AV2)
      END DO
#endif
      RETURN
      END SUBROUTINE initialize_ocn2hyd_coupling


      SUBROUTINE initialize_ocn2hyd_routers (tile)
!
!=======================================================================
!                                                                      !
!  Initialize ocean and hyd models coupling stream.  This is the      !
!  training phase used to constuct MCT parallel interpolators and      !
!  and stablish communication patterns.                                !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mct_coupler_params
!
!  Imported variable definitions.
!
      integer, intent(in) :: tile
!
!  Local variable declarations.
!
      integer :: MyError, nprocs
      integer :: ng, ih
!
!-----------------------------------------------------------------------
!  Establish MCT router.
!-----------------------------------------------------------------------
!
      ALLOCATE(Router_H(Nocn_grids,Nhyd_grids))
!
!  Initialize routers to the hyd model component.
!
      DO ng=1,Nocn_grids
        DO ih=1,Nhyd_grids
          HYDid=hydids(ih)
#ifdef MCT_INTERP_OC2HY
          CALL Router_init (HYDid, GSMapInterp_H(ng,ih)%GSMapHYD,       &
     &                      OCN_COMM_WORLD, Router_H(ng,ih)%ROMStoHYD)
#else
          CALL Router_init (HYDid, GlobalSegMap_G(ng)%GSMapROMS,        &
     &                      OCN_COMM_WORLD, Router_H(ng,ih)%ROMStoHYD)
#endif
        END DO
      END DO
      RETURN
      END SUBROUTINE initialize_ocn2hyd_routers

      SUBROUTINE ocn2hyd_coupling (ng, ih, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine acquires the coupling data streams from the ocean   !
!  to the hydrology model. Currently, the following data streams are   !
!  coded:                                                              !
!                                                                      !
!     (...) WRFHydro  units                                            !
!     [...] ROMS units                                                 !
!                                                                      !
!  Fields imported from hydro model:                                   !
!                                                                      !
!     * QRIVER  River flow [m3/s]                                      !
!     * WATERLEVEL Water surface elevation [m]                         !
!                                                                      !
!  Fields exported to Hydro model:                                     !
!                                                                      !
!     * BATH  bathymetry [m]                                           !
!     * ZETA  water surface elevation [m]                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ih, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocn2hyd_coupling_tile (ng, ih, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
#ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 48)
#endif

      RETURN
      END SUBROUTINE ocn2hyd_coupling
!
!***********************************************************************
      SUBROUTINE ocn2hyd_coupling_tile (ng, ih, tile,                   &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mct_coupler_params
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
#if defined CURVGRID || defined MASKING
      USE mod_grid
#endif
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ih, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Asize, Iimport, Iexport, MyError, Tag
      integer :: gtype, i, id, ifield, j, ij,  status

      real(r8) :: add_offset, cff, scale
      real(r8) :: RecvTime, SendTime, buffer(2), wtime(2)

      real(r8), pointer :: A(:)
#ifdef MCT_INTERP_OC2HY
      integer, pointer :: points(:)
      integer, pointer :: indices(:)
      real(r8), pointer :: Amask(:)
#endif
      real(r8) :: cff1, cff2
      character (len=40) :: code
!
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,         &
     &                          OCN_COMM_WORLD)
      allocate ( A(Asize) )
      A=0.0_r8
!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to hydrology (HYD) model.
!-----------------------------------------------------------------------
!
!  Bathymetry  (m)
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=GRID(ng)%h(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2hyd_AV, "BATH", A,  &
     &                           Asize)
!
!  Water level (m)
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=OCEAN(ng)%zeta(i,j,1)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2hyd_AV, "ZETA", A,  &
     &                           Asize)
!
!  Send ocean fields to atmosphere model.
!
      MyError=0
      Tag=ng*100+ih*10+1
!
#ifdef MCT_INTERP_OC2HY
      CALL MCT_MatVecMul(AttrVect_G(ng)%ocn2hyd_AV,                     &
     &                   SMPlus_H(ng,ih)%O2HMatPlus,                    &
     &                   AV2_H(ng,ih)%ocn2hyd_AV2)
!
!  Now add in the CPL_MASK before we send it over to wrf.
!  Get the number of grid points on this processor.
!
      Asize=GlobalSegMap_lsize (GSMapInterp_H(ng,ih)%GSMapHYD,          &
     &                          OCN_COMM_WORLD)
      allocate (Amask(Asize))
      Amask=0.0_r8
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GSMapInterp_H(ng,ih)%GSMapHYD,          &
     &                           MyRank, points)
!
!  Load the dst grid cpl mask into the attr vect.
!
!      DO i=1,Asize
!        Amask(i)=REAL(O2H_CPLMASK(ng,ih)%dst_mask(points(i)))
!      END DO
!      CALL AttrVect_importRAttr (AV2_H(ng,ih)%ocn2hyd_AV2, "CPL_MASK",  &
!     &                           Amask, Asize)
!
      CALL MCT_isend (AV2_H(ng,ih)%ocn2hyd_AV2,                         &
     &                Router_H(ng,ih)%ROMStoHYD, Tag)
#else
      CALL MCT_isend (AttrVect_G(ng)%ocn2hyd_AV,                        &
     &                Router_H(ng,ih)%ROMStoHYD, Tag)
#endif
      CALL MCT_waits (Router_H(ng,ih)%ROMStoHYD)
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) 'hydrology model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,37) ' ## ROMS grid ',ng,                        &
     &                      ' sent data to HYD grid ',ih
 37       FORMAT (a14,i2,a24,i2)
        END IF
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (A)
#ifdef MCT_INTERP_OC2HY
      deallocate (points, Amask)
      if (associated (indices)) then
        deallocate (indices)
      endif
#endif
!
 10   FORMAT (' OCN2HYD_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2HYD_COUPLING - error while sending fields to ',     &
     &        a, i4)
      RETURN
      END SUBROUTINE ocn2hyd_coupling_tile

      SUBROUTINE ocnfhyd_coupling (ng, ih, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine acquires the coupling data streams between ocean    !
!  and hydrology models. Currently, the following data streams are     !
!  coded:                                                              !
!                                                                      !
!     (...) HYD  units                                                 !
!     [...] ROMS units                                                 !
!                                                                      !
!                                                                      !
!  Fields imported from hydro model:                                   !
!                                                                      !
!     * QRIVER  River flow [m3/s]                                      !
!     * WATERLEVEL Water surface elevation [m]                         !
!                                                                      !
!  Fields exported to Hydro model:                                     !
!                                                                      !
!     * BATH  bathymetry [m]                                           !
!     * ZETA  water surface elevation [m]                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ih, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocnfhyd_coupling_tile (ng, ih, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
#ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 48)
#endif

      RETURN
      END SUBROUTINE ocnfhyd_coupling
!
!***********************************************************************
      SUBROUTINE ocnfhyd_coupling_tile (ng, ih, tile,                   &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
      USE mct_coupler_params
#if defined CURVGRID || defined MASKING
      USE mod_grid
#endif
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
#ifdef DISTRIBUTE
      USE distribute_mod,  ONLY : mp_reduce
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ih, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Asize, Iimport, Iexport, MyError, Tag
      integer :: gtype, i, id, ifield, j, ij,  status
#ifdef MCT_INTERP_OC2HY
      integer, pointer :: points(:)
      integer, pointer :: indices(:)
#endif

      real(r8) :: add_offset, cff, fac, scale
      real(r8) :: RecvTime, SendTime, buffer(2), wtime(2)
      real(r8) :: BBR, cff1, cff2
!     real(r8), parameter ::  Large = 1.0E+20_r8
      real(r8), pointer :: A(:)
      real(r8), dimension(2) :: range

      character (len=40) :: code
#ifdef DISTRIBUTE
      character (len=3), dimension(2) :: op_handle
#endif
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
#include "set_bounds.h"
#ifdef DISTRIBUTE
      op_handle(1)='MIN'
      op_handle(2)='MAX'
#endif
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,         &
     &                          OCN_COMM_WORLD)
      allocate ( A(Asize) )
      A=0.0_r8

#ifdef MCT_INTERP_OC2HY
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GlobalSegMap_G(ng)%GSMapROMS,          &
     &                           MyRank, points)
#endif

!
!-----------------------------------------------------------------------
!  Import fields from atmosphere model (HYD) to ocean model (ROMS).
!-----------------------------------------------------------------------
!
!  Receive fields from hydrology model.
!
      MyError=0
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      Tag=ng*100+ih*10+0
#ifdef MCT_INTERP_OC2HY
      CALL MCT_irecv (AV2_H(ng,ih)%hyd2ocn_AV2,                         &
     &                Router_H(ng,ih)%ROMStoHYD, Tag)
!     Wait to make sure the HYD data has arrived.
      CALL MCT_waitr (AV2_H(ng,ih)%hyd2ocn_AV2,                         &
     &                Router_H(ng,ih)%ROMStoHYD)
      CALL MCT_MatVecMul(AV2_H(ng,ih)%hyd2ocn_AV2,                      &
     &                   SMPlus_H(ng,ih)%H2OMatPlus,                    &
     &                   AttrVect_G(ng)%hyd2ocn_AV)
#else
      CALL MCT_irecv (AttrVect_G(ng)%hyd2ocn_AV,                        &
     &                Router_H(ng,ih)%ROMStoHYD, Tag)
!     Wait to make sure the HYD data has arrived.
      CALL MCT_waitr (AttrVect_G(ng)%hyd2ocn_AV,                        &
     &                Router_H(ng,ih)%ROMStoHYD)
#endif
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'hydrology model, MyError = ', MyError
        END IF
        CALL finalize_ocn2hyd_coupling
      ELSE
        IF (Master) THEN
          WRITE (stdout,38) ' ## ROMS grid ',ng,                        &
     &                      ' recv data from HYD grid ',ih
 38       FORMAT (a14,i2,a25,i2)
        END IF
      END IF
!
!  Receive fields from hydrology model.
 40         FORMAT (a36,1x,2(1pe14.6))
!
	!  River Flow          (m3/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%hyd2ocn_AV, "QRIVER",   &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=1.0_r8
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ih.eq.1) THEN
!            FORCES(ng)%Qriver(i,j)=cff*fac
          ELSE
!            FORCES(ng)%Qriver(i,j)=FORCES(ng)%Qriver(i,j)+cff*fac
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'HYDtoROMS Min/Max  Qriver     (m3/s):  ',     &
     &                    range(1),range(2)
      END IF
!
!  Waterlevel          (m)
!
      CALL AttrVect_exportRAttr(AttrVect_G(ng)%hyd2ocn_AV, "WATERLEVEL",&
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=1.0_r8
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
!#ifdef MCT_INTERP_OC2AT
!          cff=A(ij)*REAL(A2O_CPLMASK(ih,ng)%dst_mask(points(ij)))
!#else
          cff=A(ij)
!#endif
          IF (ih.eq.1) THEN
!            FORCES(ng)%Waterlevel(i,j)=cff*fac
          ELSE
!            FORCES(ng)%Waterlevel(i,j)=FORCES(ng)%Waterlevel(i,j)+      &
!     &                                 cff*fac
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'HYDtoROMS Min/Max  Waterlevle    (m):  ',     &
     &                    range(1),range(2)
      END IF
!
!  Apply boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
!        CALL exchange_r2d_tile (ng, tile,                               &
!     &                          LBi, UBi, LBj, UBj,                     &
!     &                          FORCES(ng)%Qriver)
!        CALL exchange_r2d_tile (ng, tile,                               &
!     &                          LBi, UBi, LBj, UBj,                     &
!     &                          FORCES(ng)%Waterlevel)
      END IF
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile boundaries.
!-----------------------------------------------------------------------
!
!      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!     &                    LBi, UBi, LBj, UBj,                           &
!     &                    NghostPoints,                                 &
!     &                    EWperiodic(ng), NSperiodic(ng),               &
!     &                    FORCES(ng)%Qriver, FORCES(ng)%Waterlevel)
#endif
!
!  Deallocate communication arrays.
!
      deallocate (A)
#ifdef MCT_INTERP_OC2HY
      deallocate (points)
      if (associated (indices)) then
        deallocate (indices)
      endif
#endif
!
 10   FORMAT (' OCNFHYD_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCNFHYD_COUPLING - error while sending fields to ',     &
     &        a, i4)
      RETURN
      END SUBROUTINE ocnfhyd_coupling_tile

      SUBROUTINE finalize_ocn2hyd_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and atmosphere models coupling data     !
!  streams.                                                             !
!                                                                       !
!========================================================================
      USE mod_scalars
      USE mct_coupler_params
!
!  Local variable declarations.
!
      integer :: ng, ih, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
#if !defined WAVES_OCEAN
      deallocate ( hydids )
      deallocate ( ocnids )
#endif
      DO ng=1,Nocn_grids
        CALL AttrVect_clean (AttrVect_G(ng)%ocn2hyd_AV, MyError)
        CALL GlobalSegMap_clean (GlobalSegMap_G(ng)%GSMapROMS, MyError)
        DO ih=1,Nhyd_grids
          CALL Router_clean (Router_H(ng,ih)%ROMStoHYD, MyError)
        END DO
      END DO
      RETURN

      END SUBROUTINE finalize_ocn2hyd_coupling
