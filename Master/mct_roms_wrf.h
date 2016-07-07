/*
** svn $Id: mct_roms_wrf.h 756 2008-09-14 20:18:28Z jcwarner $
***************************************************** John C. Warner ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group      Hernan G. Arango  **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** These routines are use couple ROMS/TOMS to WRF atmosphere model    **
** using the Model Coupling Toolkit (MCT).                            **
**                                                                    **
************************************************************************
*/

      SUBROUTINE initialize_ocn2atm_coupling (ng, tile)
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
#ifdef MCT_INTERP_OC2AT
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
      integer :: i, ic, ia, j, jc, nprocs, cid, cad
      integer :: nRows, nCols, num_sparse_elems

      real(r8) :: cff

      integer, allocatable  :: length(:)
      integer, allocatable  :: start(:)
      character (len=70) :: nc_name
      character (len=20) :: to_add
      character (len=120) :: aostring
      character (len=120) :: oastring
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
#ifdef MCT_INTERP_OC2AT
      IF (ng.eq.1) THEN
        ALLOCATE(SMPlus_A(Nocn_grids,Natm_grids))
        ALLOCATE(AV2_A(Nocn_grids,Natm_grids))
        ALLOCATE(GSMapInterp_A(Nocn_grids,Natm_grids))
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
#ifdef MCT_INTERP_OC2AT
!
!  If ocean grid and atm grids are different sizes, then
!  develop sparse matrices for interpolation.
!
  35  FORMAT(a3,i1,a7,i1,a11)
      DO ia=1,Natm_grids
!
!!!!!!!!!!!!!!!!!!!!!!
! First work on atm to ocean.
!!!!!!!!!!!!!!!!!!!!!!
!
        IF (Myrank.eq.MyMaster) THEN
          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'atm',ia,'_to_ocn',ng,'_weights.nc'
          ELSE
            nc_name=A2Oname(ia,ng)
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
          allocate(A2O_CPLMASK(ia,ng)%dst_mask(nRows))
          DO i=1,nRows
            A2O_CPLMASK(ia,ng)%dst_mask(i)=dst_grid_imask(i)
          END DO


          call SparseMatrix_init(sMatA,nRows,nCols,num_sparse_elems)
          call SparseMatrix_importGRowInd(sMatA, sparse_rows,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importGColInd(sMatA, sparse_cols,           &
     &                                    num_sparse_elems)
          call SparseMatrix_importMatrixElts(sMatA, sparse_weights,     &
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
        IF (Myrank.ne.MyMaster) THEN
          nRows=dst_grid_dims(1)*dst_grid_dims(2)
          allocate(A2O_CPLMASK(ia,ng)%dst_mask(nRows))
        END IF
        CALL mpi_bcast(A2O_CPLMASK(ia,ng)%dst_mask,nRows,               &
     &                 MPI_INTEGER, MyMaster,                           &
     &                 OCN_COMM_WORLD, MyError)


!!!!!!!!!!!!!!!!!!!!!!
! Second work on ocean to atm.
!!!!!!!!!!!!!!!!!!!!!!
!
        IF (Myrank.eq.MyMaster) THEN
          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'ocn',ng,'_to_atm',ia,'_weights.nc'
          ELSE
            nc_name=O2Aname(ng,ia)
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
          allocate(O2A_CPLMASK(ng,ia)%dst_mask(nRows))
          DO i=1,nRows
            O2A_CPLMASK(ng,ia)%dst_mask(i)=dst_grid_imask(i)
          END DO
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
        IF (Myrank.ne.MyMaster) THEN
          nRows=dst_grid_dims(1)*dst_grid_dims(2)
          allocate(O2A_CPLMASK(ng,ia)%dst_mask(nRows))
        END IF
        CALL mpi_bcast(O2A_CPLMASK(ng,ia)%dst_mask,nRows,               &
     &                 MPI_INTEGER, MyMaster,                           &
     &                 OCN_COMM_WORLD, MyError)
!
!  Create Global Seg Map for atm model.
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
        CALL GlobalSegMap_init (GSMapInterp_A(ng,ia)%GSMapWRF,          &
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
        call SparseMatrixPlus_init(SMPlus_A(ng,ia)%A2OMatPlus, sMatA,   &
     &                             GSMapInterp_A(ng,ia)%GSMapWRF,       &
     &                             GlobalSegMap_G(ng)%GSMapROMS,        &
     &                             Xonly,MyMaster,OCN_COMM_WORLD, OCNid)
        call SparseMatrix_clean(sMatA)
!
! Create Ocean sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
        call SparseMatrixPlus_init(SMPlus_A(ng,ia)%O2AMatPlus, sMatO,   &
     &                             GlobalSegMap_G(ng)%GSMapROMS,        &
     &                             GSMapInterp_A(ng,ia)%GSMapWRF,       &
     &                             Xonly,MyMaster,OCN_COMM_WORLD, OCNid)
        call SparseMatrix_clean(sMatO)
      END DO
#endif
!
!  Initialize attribute vector holding the data code strings from
!  the atmosphere model.
!
      cad=LEN(aostring)
      DO i=1,cad
        aostring(i:i)=''
      END DO
      cid=1
!
      to_add='GSW'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':GLW'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
#ifdef ATM2OCN_FLUXES
      to_add=':LH'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':HFX'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':USTRESS'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VSTRESS'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
!
#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      to_add=':MSLP'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
!
#if defined BULK_FLUXES || defined ECOSIM || \
   (defined SHORTWAVE && defined ANA_SRFLUX)
      to_add=':RELH'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':T2'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
!
#if defined BULK_FLUXES || defined ECOSIM
      to_add=':U10'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':V10'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
!
#ifdef CLOUDS
      to_add=':CLDFRA'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
!
#if !defined ANA_RAIN && defined EMINUSP
      to_add=':RAIN'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
!
#if defined EMINUSP
      to_add=':EVAP'
      cad=LEN_TRIM(to_add)
      write(aostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
!
!  Finalize and remove trailing spaces from the aostring
!  for the rlist.
!
      cad=LEN_TRIM(aostring)
      aostring=aostring(1:cad)
!
!  Initialize attribute vector holding the export data code strings of
!  the atmosphere model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &                         OCN_COMM_WORLD)
      CALL AttrVect_init(AttrVect_G(ng)%atm2ocn_AV,                     &
     &                   rList=TRIM(aostring),lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%atm2ocn_AV)
!
!  Initialize attribute vector holding the export data of
!  the atm model.
!
!
!  Initialize attribute vector that contain the data strings from
!  the ocean model.
!
      cad=LEN(oastring)
      DO i=1,cad
        oastring(i:i)=''
      END DO
      cid=1
!
      to_add='SST'
      cad=LEN_TRIM(to_add)
      write(oastring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
#ifdef MCT_INTERP_OC2AT
      to_add=':CPL_MASK'
      cad=LEN_TRIM(to_add)
      write(oastring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
!
!  Finalize and remove trailing spaces from the oastring
!  for the rlist.
!
      cad=LEN_TRIM(oastring)
      oastring=oastring(1:cad)
!
      CALL AttrVect_init (AttrVect_G(ng)%ocn2atm_AV,                    &
     &                    rList=TRIM(oastring),lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%ocn2atm_AV)
!
#ifdef MCT_INTERP_OC2AT
      DO ia=1,Natm_grids
        ATMid=atmids(ia)
!
!  Initialize attribute vector holding the export data of
!  the atm model. The Asize is the number of grid point on this
!  processor.
!
        Asize=GlobalSegMap_lsize(GSMapInterp_A(ng,ia)%GSMapWRF,         &
     &                           OCN_COMM_WORLD)
        CALL AttrVect_init (AV2_A(ng,ia)%atm2ocn_AV2,                   &
     &                      rList=TRIM(aostring),lsize=Asize)
        CALL AttrVect_zero (AV2_A(ng,ia)%atm2ocn_AV2)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
        CALL AttrVect_init (AV2_A(ng,ia)%ocn2atm_AV2,                   &
     &                      rList=TRIM(oastring),lsize=Asize)
        CALL AttrVect_zero (AV2_A(ng,ia)%ocn2atm_AV2)
      END DO
#endif
      RETURN
      END SUBROUTINE initialize_ocn2atm_coupling


      SUBROUTINE initialize_ocn2atm_routers (tile)
!
!=======================================================================
!                                                                      !
!  Initialize ocean and atm models coupling stream.  This is the      !
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
      integer :: ng, ia
!
!-----------------------------------------------------------------------
!  Establish MCT router.
!-----------------------------------------------------------------------
!
      ALLOCATE(Router_A(Nocn_grids,Natm_grids))
!
!  Initialize routers to the wave model component.
!
      DO ng=1,Nocn_grids
        DO ia=1,Natm_grids
          ATMid=atmids(ia)
#ifdef MCT_INTERP_OC2AT
          CALL Router_init (ATMid, GSMapInterp_A(ng,ia)%GSMapWRF,       &
     &                      OCN_COMM_WORLD, Router_A(ng,ia)%ROMStoWRF)
#else
          CALL Router_init (ATMid, GlobalSegMap_G(ng)%GSMapROMS,        &
     &                      OCN_COMM_WORLD, Router_A(ng,ia)%ROMStoWRF)
#endif
        END DO
      END DO
      RETURN
      END SUBROUTINE initialize_ocn2atm_routers

      SUBROUTINE ocn2atm_coupling (ng, ia, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine acquires the coupling data streams from the ocean   !
!  to the atmosphere model. Currently, the following data streams are  !
!  coded:                                                              !
!                                                                      !
!     (...) WRF  units                                                 !
!     [...] ROMS units                                                 !
!                                                                      !
!  Fields imported WRF model:                                          !
!                                                                      !
!     * GSW     Net shortwave radiation (Watts/m2), [Celsius m/s]      !
!     * GLW     Long wave raditaion (Watts/m2), [Celsius m/s]          !
!     * LH      Latent heat flux (Watts/m2), [Celsius m/s]             !
!     * HFX     Sensible heat flux (Watts/m2), [Celsius m/s]           !
!     * USTRESS Surface U-wind stress (Pa), [m2/s2]                    !
!     * VSTRESS Surface V-wind stress (Pa), [m2/s2]                    !
!     * MSLP    Mean Sea Level Pressure (Pa), [mb]                     !
!     * RELH    Surface air relative humidity (percent), [fraction]    !
!     * T2      Surface (2 m) air temperature (Celsius), [Celsius]     !
!     * U10     Surface (10 m) U-wind speed (m/s), [m/s]               !
!     * V10     Surface (10 m) V-wind speed (m/s), [m/s]               !
!     * CLDFRA  Cloud fraction (percent/100), [percent/100]            !
!     * RAIN    Precipitation (m/s), [kg/m2/s]                         !
!     * EVAP    Evaporation (m/s), [kg/m2/s]                           !
!                                                                      !
!  Fields exported to WRF model:                                       !
!                                                                      !
!     * SST     Sea surface potential temperature (Kelvin), [Celsius]  !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ia, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocn2atm_coupling_tile (ng, ia, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif

      RETURN
      END SUBROUTINE ocn2atm_coupling
!
!***********************************************************************
      SUBROUTINE ocn2atm_coupling_tile (ng, ia, tile,                   &
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
      integer, intent(in) :: ng, ia, tile
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
#ifdef MCT_INTERP_OC2AT
      integer, pointer :: points(:)
      integer, pointer :: indices(:)
      real(r8), pointer :: Amask(:)
#endif
      real(r8) :: BBR, cff1, cff2
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
!  Export fields from ocean (ROMS) to atmosphere (WRF) model.
!-----------------------------------------------------------------------
!
!  Sea surface temperature       (degC)
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2atm_AV, "SST", A,   &
     &                           Asize)
!
!  Send ocean fields to atmosphere model.
!
      MyError=0
      Tag=ng*100+ia*10+0
!
#ifdef MCT_INTERP_OC2AT
      CALL MCT_MatVecMul(AttrVect_G(ng)%ocn2atm_AV,                     &
     &                   SMPlus_A(ng,ia)%O2AMatPlus,                    &
     &                   AV2_A(ng,ia)%ocn2atm_AV2)
!
!  Now add in the CPL_MASK before we send it over to wrf.
!  Get the number of grid points on this processor.
!
      Asize=GlobalSegMap_lsize (GSMapInterp_A(ng,ia)%GSMapWRF,          &
     &                          OCN_COMM_WORLD)
      allocate (Amask(Asize))
      Amask=0.0_r8
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GSMapInterp_A(ng,ia)%GSMapWRF,          &
     &                           MyRank, points)
!
!  Load the dst grid cpl mask into the attr vect.
!
      DO i=1,Asize
        Amask(i)=REAL(O2A_CPLMASK(ng,ia)%dst_mask(points(i)))
      END DO
      CALL AttrVect_importRAttr (AV2_A(ng,ia)%ocn2atm_AV2, "CPL_MASK",  &
     &                           Amask, Asize)
!
      CALL MCT_isend (AV2_A(ng,ia)%ocn2atm_AV2,                         &
     &                Router_A(ng,ia)%ROMStoWRF, Tag)
#else
      CALL MCT_isend (AttrVect_G(ng)%ocn2atm_AV,                        &
     &                Router_A(ng,ia)%ROMStoWRF, Tag)
#endif
      CALL MCT_waits (Router_A(ng,ia)%ROMStoWRF)
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) 'atmosphere model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,37) ' ## ROMS grid ',ng,                        &
     &                      ' sent data to WRF grid ',ia
 37       FORMAT (a14,i2,a24,i2)
        END IF
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (A)
#ifdef MCT_INTERP_OC2AT
      deallocate (points, Amask)
      if (associated (indices)) then
        deallocate (indices)
      endif
#endif
!
 10   FORMAT (' OCN2ATM_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2ATM_COUPLING - error while sending fields to ',     &
     &        a, i4)
      RETURN
      END SUBROUTINE ocn2atm_coupling_tile

      SUBROUTINE ocnfatm_coupling (ng, ia, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine acquires the coupling data streams between ocean    !
!  and atmosphere models. Currently, the following data streams are    !
!  coded:                                                              !
!                                                                      !
!     (...) WRF  units                                                 !
!     [...] ROMS units                                                 !
!                                                                      !
!  Fields imported WRF model:                                          !
!                                                                      !
!     * GSW     Net shortwave radiation (Watts/m2), [Celsius m/s]      !
!     * GLW     Long wave raditaion (Watts/m2), [Celsius m/s]          !
!     * LH      Latent heat flux (Watts/m2), [Celsius m/s]             !
!     * HFX     Sensible heat flux (Watts/m2), [Celsius m/s]           !
!     * USTRESS Surface U-wind stress (Pa), [m2/s2]                    !
!     * VSTRESS Surface V-wind stress (Pa), [m2/s2]                    !
!     * MSLP    Mean Sea Level Pressure (Pa), [mb]                     !
!     * RELH    Surface air relative humidity (percent), [fraction]    !
!     * T2      Surface (2 m) air temperature (Celsius), [Celsius]     !
!     * U10     Surface (10 m) U-wind speed (m/s), [m/s]               !
!     * V10     Surface (10 m) V-wind speed (m/s), [m/s]               !
!     * CLDFRA  Cloud fraction (percent/100), [percent/100]            !
!     * RAIN    Precipitation (m/s), [kg/m2/s]                         !
!     * EVAP    Evaporation (m/s), [kg/m2/s]                           !
!                                                                      !
!  Fields exported to WRF model:                                       !
!                                                                      !
!     * SST     Sea surface potential temperature (Kelvin), [Celsius]  !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ia, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocnfatm_coupling_tile (ng, ia, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif

      RETURN
      END SUBROUTINE ocnfatm_coupling
!
!***********************************************************************
      SUBROUTINE ocnfatm_coupling_tile (ng, ia, tile,                   &
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
      integer, intent(in) :: ng, ia, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Asize, Iimport, Iexport, MyError, Tag
      integer :: gtype, i, id, ifield, j, ij,  status
#ifdef MCT_INTERP_OC2AT
      integer, pointer :: points(:)
      integer, pointer :: indices(:)
#endif

      real(r8) :: add_offset, cff, fac, scale
      real(r8) :: RecvTime, SendTime, buffer(2), wtime(2)
      real(r8) :: BBR, cff1, cff2
      real(r8), parameter ::  Large = 1.0E+20_r8
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

#ifdef MCT_INTERP_OC2AT
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GlobalSegMap_G(ng)%GSMapROMS,          &
     &                           MyRank, points)
#endif

!
!-----------------------------------------------------------------------
!  Import fields from atmosphere model (WRF) to ocean model (ROMS).
!-----------------------------------------------------------------------
!
!  Receive fields from atmosphere model.
!
      MyError=0
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      Tag=ng*100+ia*10+0
#ifdef MCT_INTERP_OC2AT
      CALL MCT_irecv (AV2_A(ng,ia)%atm2ocn_AV2,                         &
     &                Router_A(ng,ia)%ROMStoWRF, Tag)
!     Wait to make sure the WRF data has arrived.
      CALL MCT_waitr (AV2_A(ng,ia)%atm2ocn_AV2,                         &
     &                Router_A(ng,ia)%ROMStoWRF)
      CALL MCT_MatVecMul(AV2_A(ng,ia)%atm2ocn_AV2,                      &
     &                   SMPlus_A(ng,ia)%A2OMatPlus,                    &
     &                   AttrVect_G(ng)%atm2ocn_AV)
#else
      CALL MCT_irecv (AttrVect_G(ng)%atm2ocn_AV,                        &
     &                Router_A(ng,ia)%ROMStoWRF, Tag)
!     Wait to make sure the WRF data has arrived.
      CALL MCT_waitr (AttrVect_G(ng)%atm2ocn_AV,                        &
     &                Router_A(ng,ia)%ROMStoWRF)
#endif
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'atmosphere model, MyError = ', MyError
        END IF
        CALL finalize_ocn2atm_coupling
      ELSE
        IF (Master) THEN
          WRITE (stdout,38) ' ## ROMS grid ',ng,                        &
     &                      ' recv data from WRF grid ',ia
 38       FORMAT (a14,i2,a25,i2)
        END IF
      END IF
!
!  Receive fields from atmosphere model.
 40         FORMAT (a36,1x,2(1pe14.6))
!
!  Short wave radiation          (from W/m^2 to Celsius m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "GSW",      &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=1.0_r8/(rho0*Cp)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ia.eq.1) THEN
            FORCES(ng)%srflx(i,j)=cff*fac
          ELSE
            FORCES(ng)%srflx(i,j)=FORCES(ng)%srflx(i,j)+cff*fac
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  GSW     (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
!
!  Long wave radiation          (from W/m^2 to Celsius m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "GLW",      &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=1.0_r8/(rho0*Cp)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
!         BBR=OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)+273.16_r8
!         BBR=BBR*BBR*BBR*BBR
!         BBR=0.97_r8*5.67E-8_r8*BBR
!         A(ij)=A(ij)-BBR
#ifdef MCT_INTERP_OC2AT
          cff=A(ij)*REAL(A2O_CPLMASK(ia,ng)%dst_mask(points(ij)))
#else
          cff=A(ij)
#endif
          IF (ia.eq.1) THEN
            FORCES(ng)%lrflx(i,j)=cff*fac
          ELSE
            FORCES(ng)%lrflx(i,j)=FORCES(ng)%lrflx(i,j)+cff*fac
          END IF
          IF (ia.eq.Natm_grids) THEN
            BBR=OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)+273.16_r8
            BBR=BBR*BBR*BBR*BBR
            BBR=0.97_r8*5.67E-8_r8*BBR
            FORCES(ng)%lrflx(i,j)=FORCES(ng)%lrflx(i,j)-BBR*fac
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  GLW     (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
#ifdef ATM2OCN_FLUXES
!
!  Latent heat flux            (from W/m^2 to Celsius m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "LH",       &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=-1.0_r8/(rho0*Cp)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ia.eq.1) THEN
            FORCES(ng)%lhflx(i,j)=cff*fac
          ELSE
            FORCES(ng)%lhflx(i,j)=FORCES(ng)%lhflx(i,j)+cff*fac
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  LH      (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
!
!  Sensible heat flux            (from W/m^2 to Celsius m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "HFX",      &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=-1.0_r8/(rho0*Cp)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ia.eq.1) THEN
            FORCES(ng)%shflx(i,j)=cff*fac
          ELSE
            FORCES(ng)%shflx(i,j)=FORCES(ng)%shflx(i,j)+cff*fac
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  HFX     (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
!
!  Surface u-stress              (m2/s2)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "USTRESS",  &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=1.0_r8/rho0
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ia.eq.1) THEN
            FORCES(ng)%Taux(i,j)=cff*fac
          ELSE
            FORCES(ng)%Taux(i,j)=FORCES(ng)%Taux(i,j)+cff*fac
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  USTRESS (Nm-2):  ',        &
     &                    range(1),range(2)
      END IF
!
!  Surface v-stress              (m2/s2)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "VSTRESS",  &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=1.0_r8/rho0
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ia.eq.1) THEN
            FORCES(ng)%Tauy(i,j)=cff*fac
          ELSE
            FORCES(ng)%Tauy(i,j)=FORCES(ng)%Tauy(i,j)+cff*fac
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  VSTRESS (Nm-2):  ',        &
     &                    range(1),range(2)
      END IF
# ifdef CURVGRID
!
!  Rotate to curvilinear grid.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff1=FORCES(ng)%Taux(i,j)*GRID(ng)%CosAngler(i,j)+            &
     &         FORCES(ng)%Tauy(i,j)*GRID(ng)%SinAngler(i,j)
          cff2=FORCES(ng)%Tauy(i,j)*GRID(ng)%CosAngler(i,j)-            &
     &         FORCES(ng)%Taux(i,j)*GRID(ng)%SinAngler(i,j)
          FORCES(ng)%Taux(i,j)=cff1
          FORCES(ng)%Tauy(i,j)=cff2
        END DO
      END DO
# endif
#endif
#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Mean seal level pressure, convert from Pa to mb.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "MSLP",     &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      fac=0.01_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)*fac
          IF (ia.eq.1) THEN
            FORCES(ng)%Pair(i,j)=cff
          ELSE
            FORCES(ng)%Pair(i,j)=FORCES(ng)%Pair(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  MSLP    (mb):    ',        &
     &                    range(1),range(2)
      END IF
#endif
#if defined BULK_FLUXES || defined ECOSIM || \
   (defined SHORTWAVE && defined ANA_SRFLUX)
!
!  Surface air relative humidity (-)
!  Convert RELH from percent to fraction.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "RELH",     &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      fac=0.01_r8
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)*fac
          IF (ia.eq.1) THEN
            FORCES(ng)%Hair(i,j)=cff
          ELSE
            FORCES(ng)%Hair(i,j)=FORCES(ng)%Hair(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  RELH    (-):     ',        &
     &                    range(1),range(2)
      END IF
!
!  Surface 2m air temperature    (degC)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "T2",       &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ia.eq.1) THEN
            FORCES(ng)%Tair(i,j)=cff
          ELSE
            FORCES(ng)%Tair(i,j)=FORCES(ng)%Tair(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  T2      (C):     ',        &
     &                    range(1),range(2)
      END IF
#endif
#if defined BULK_FLUXES || defined ECOSIM
!
!  U-Wind speed at 10 m          (m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "U10",      &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ia.eq.1) THEN
            FORCES(ng)%Uwind(i,j)=cff
          ELSE
            FORCES(ng)%Uwind(i,j)=FORCES(ng)%Uwind(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  U10     (ms-1):  ',        &
     &                    range(1),range(2)
      END IF
!
!  V-Wind speed at 10 m          (m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "V10",      &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ia.eq.1) THEN
            FORCES(ng)%Vwind(i,j)=cff
          ELSE
            FORCES(ng)%Vwind(i,j)=FORCES(ng)%Vwind(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  V10     (ms-1):  ',        &
     &                    range(1),range(2)
      END IF
# ifdef CURVGRID
!
!  Rotate to curvilinear grid.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=FORCES(ng)%Uwind(i,j)*GRID(ng)%CosAngler(i,j)+         &
     &           FORCES(ng)%Vwind(i,j)*GRID(ng)%SinAngler(i,j)
            cff2=FORCES(ng)%Vwind(i,j)*GRID(ng)%CosAngler(i,j)-         &
     &           FORCES(ng)%Uwind(i,j)*GRID(ng)%SinAngler(i,j)
            FORCES(ng)%Uwind(i,j)=cff1
            FORCES(ng)%Vwind(i,j)=cff2
          END DO
        END DO
# endif
#endif
#ifdef CLOUDS
!
!  Cloud fraction                (percent/100, so 0-1)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "CLDFRA",   &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)
          IF (ia.eq.1) THEN
            FORCES(ng)%cloud(i,j)=cff
          ELSE
            FORCES(ng)%cloud(i,j)=FORCES(ng)%cloud(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  CLDFRA  (-):     ',        &
     &                    range(1),range(2)
      END IF
#endif
#if !defined ANA_RAIN && defined EMINUSP
!
!  Precipitation                 (convert to kg/m2/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "RAIN",     &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=rho0
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)*fac
          IF (ia.eq.1) THEN
            FORCES(ng)%rain(i,j)=cff
          ELSE
            FORCES(ng)%rain(i,j)=FORCES(ng)%rain(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  RAIN  (kgm-2s-1):',        &
     &                    range(1),range(2)
      END IF
#endif
#if defined EMINUSP
!
!  Evaporation                 (convert to kg/m2/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "EVAP",     &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      fac=rho0
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=A(ij)*fac
          IF (ia.eq.1) THEN
            FORCES(ng)%evap(i,j)=cff
          ELSE
            FORCES(ng)%evap(i,j)=FORCES(ng)%evap(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'WRFtoROMS Min/Max  EVAP  (kgm-2s-1):',        &
     &                    range(1),range(2)
      END IF
#endif
!
!  Apply boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%srflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%lrflx)
# ifdef ATM2OCN_FLUXES
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%lhflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%shflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Taux)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Tauy)
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Pair)
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Hair)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Tair)
# endif
# if defined BULK_FLUXES || defined ECOSIM
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Uwind)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Vwind)
# endif
# ifdef CLOUDS
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%cloud)
# endif
# if !defined ANA_RAIN && defined EMINUSP
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%rain)
# endif
# if defined EMINUSP
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%evap)
# endif
      END IF
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile boundaries.
!-----------------------------------------------------------------------
!
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%srflx, FORCES(ng)%lrflx)
# if defined ATM2OCN_FLUXES
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%lhflx, FORCES(ng)%shflx)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                           &
     &                    LBi, UBi, LBj, UBj,                          &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Taux, FORCES(ng)%Tauy)
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Pair)
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Hair, FORCES(ng)%Tair)
# endif
# if defined BULK_FLUXES || defined ECOSIM
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Uwind, FORCES(ng)%Vwind)
# endif
# ifdef CLOUDS
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%cloud)
# endif
# if !defined ANA_RAIN && defined EMINUSP
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%rain)
# endif
# if defined EMINUSP
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%evap)
# endif
#endif
!
!  Deallocate communication arrays.
!
      deallocate (A)
#ifdef MCT_INTERP_OC2AT
      deallocate (points)
      if (associated (indices)) then
        deallocate (indices)
      endif
#endif
!
 10   FORMAT (' OCNFATM_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCNFATM_COUPLING - error while sending fields to ',     &
     &        a, i4)
      RETURN
      END SUBROUTINE ocnfatm_coupling_tile

      SUBROUTINE finalize_ocn2atm_coupling
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
      integer :: ng, ia, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
#if !defined WAVES_OCEAN
      deallocate ( atmids )
      deallocate ( ocnids )
#endif
      DO ng=1,Nocn_grids
        CALL AttrVect_clean (AttrVect_G(ng)%ocn2atm_AV, MyError)
        CALL GlobalSegMap_clean (GlobalSegMap_G(ng)%GSMapROMS, MyError)
        DO ia=1,Natm_grids
          CALL Router_clean (Router_A(ng,ia)%ROMStoWRF, MyError)
        END DO
      END DO
      RETURN

      END SUBROUTINE finalize_ocn2atm_coupling
