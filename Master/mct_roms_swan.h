/*
** svn $Id: mct_roms_swan.h 756 2008-09-14 20:18:28Z jcwarner $
***************************************************** John C. Warner ***
** Copyright (c) 2002-2008 The ROMS/TOMS Group      Hernan G. Arango  **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** These routines are use couple ROMS/TOMS to SWAN wave model using   **
** the Model Coupling Toolkit (MCT).                                  **
**                                                                    **
************************************************************************
*/

      SUBROUTINE initialize_ocn2wav_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  Initialize ocean and wave models coupling stream.  This is the      !
!  training phase used to constuct MCT parallel interpolators and      !
!  and stablish communication patterns.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_kinds
      USE mod_scalars
      USE mod_iounits
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
      integer :: i, ic, j, jc, nprocs
      integer :: nRows, nCols, num_sparse_elems
!     integer, dimension(:), pointer :: points
      integer :: ioff, joff, ieff
#ifdef REFINED_GRID
      integer, dimension(:), pointer :: ocnids
      integer, dimension(:), pointer :: wavids
#endif
      integer, allocatable :: length(:)
      integer, allocatable :: start(:)
      integer, dimension(2) :: src_grid_dims, dst_grid_dims
      character (len=70)    :: nc_name
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
      IF (WESTERN_EDGE) THEN
        IstrR=BOUNDS(ng)%Istr(tile)-1
      ELSE
        IstrR=BOUNDS(ng)%Istr(tile)
      END IF
      IF (EASTERN_EDGE) THEN
        IendR=BOUNDS(ng)%Iend(tile)+1
      ELSE
        IendR=BOUNDS(ng)%Iend(tile)
      END IF
      IF (SOUTHERN_EDGE) THEN
        JstrR=BOUNDS(ng)%Jstr(tile)-1
      ELSE
        JstrR=BOUNDS(ng)%Jstr(tile)
      END IF
      IF (NORTHERN_EDGE) THEN
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
!  Initialize MCT coupled model registry.
!
#ifdef REFINED_GRID
      allocate ( ocnids(Ngrids) )
      allocate ( wavids(Ngrids) )
      DO i=1,Ngrids
        ocnids(i)=i
        wavids(i)=Ngrids+i
      END DO
      WAVid=wavids(ng)
      OCNid=ocnids(ng) 
      CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                    OCN_COMM_WORLD, OCNid, ocnids)
#else
      CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD, OCN_COMM_WORLD,  &
     &                    OCNid)
#endif
#ifdef MCT_INTERP_OC2WV
!
!  If ocean grid and wave grids are different sizes, then
!  develop sparse matrices for interpolation.
!
!!!!!!!!!!!!!!!!!!!!!!
! First work on waves to ocean.
!!!!!!!!!!!!!!!!!!!!!!
!
      IF (Myrank.eq.MyMaster) THEN
        nc_name=SP1name(ng)
        call get_sparse_matrix (ng, nc_name, num_sparse_elems,          &
     &                          src_grid_dims, dst_grid_dims)
!
! Init the sparse matrix.
!
        nRows=dst_grid_dims(1)*dst_grid_dims(2)
        nCols=src_grid_dims(1)*src_grid_dims(2)
!
! Create sparse matrix.
!
        call SparseMatrix_init(sMatW,nRows,nCols,num_sparse_elems)
        call SparseMatrix_importGRowInd(sMatW, sparse_rows,             &
     &                                  size(sparse_rows))
        call SparseMatrix_importGColInd(sMatW, sparse_cols,             &
     &                                  size(sparse_cols))
        call SparseMatrix_importMatrixElts(sMatW, sparse_weights,       &
     &                                     size(sparse_weights))
!
! Deallocate arrays.
!
!        IF (allocated(sparse_rows)) THEN
          deallocate ( sparse_rows )
!        END IF
!        IF (allocated(sparse_cols)) THEN
          deallocate ( sparse_cols )
!        END IF
!        IF (allocated(sparse_weights)) THEN
          deallocate ( sparse_weights )
!        END IF
!        IF (allocated(sparse_weights)) THEN
          deallocate ( dst_grid_imask )
!        END IF
      END IF

!!!!!!!!!!!!!!!!!!!!!!
! Second work on ocean to waves.
!!!!!!!!!!!!!!!!!!!!!!
!
      IF (Myrank.eq.MyMaster) THEN
        nc_name=SP2name(ng)
        call get_sparse_matrix (ng, nc_name, num_sparse_elems,          &
     &                          src_grid_dims, dst_grid_dims)
!
! Init the sparse matrix.
!
        nRows=dst_grid_dims(1)*dst_grid_dims(2)
        nCols=src_grid_dims(1)*src_grid_dims(2)
!
! Create sparse matrix.
!
        call SparseMatrix_init(sMatO,nRows,nCols,num_sparse_elems)
        call SparseMatrix_importGRowInd(sMatO, sparse_rows,              &
     &                                  size(sparse_rows))
        call SparseMatrix_importGColInd(sMatO, sparse_cols,              &
     &                                  size(sparse_cols))
        call SparseMatrix_importMatrixElts(sMatO, sparse_weights,        &
     &                                     size(sparse_weights))
!
! Deallocate arrays.
!
!        IF (allocated(sparse_rows)) THEN
          deallocate ( sparse_rows )
!        END IF
!        IF (allocated(sparse_cols)) THEN
          deallocate ( sparse_cols )
!        END IF
!        IF (allocated(sparse_weights)) THEN
          deallocate ( sparse_weights )
!        END IF
!        IF (allocated(sparse_weights)) THEN
          deallocate ( dst_grid_imask )
!        END IF
      END IF
#endif
#if !defined MCT_INTERP_OC2WV
!
!  Determine start and lengths for roms domain decomposition.
!
      Jsize=JendT-JstrT+1
      IF (.not.allocated(start)) THEN
        allocate ( start(Jsize) )
      END IF
      IF (.not.allocated(length)) THEN
        allocate ( length(Jsize) )
      END IF
      jc=0
      ioff=0
      joff=0
      ieff=0
#ifdef REFINED_GRID
      IF (ng.gt.1) THEN
        ioff=5
        joff=3
        ieff=3
      END IF
#endif
      DO j=JstrT,JendT
        jc=jc+1
        start (jc)=(j+joff)*(Lm(ng)+2+ioff)+(IstrT+ieff)+1
        length(jc)=(IendT-IstrT+1)
      END DO
#ifdef REFINED_GRID
      IF (ng.eq.1) THEN
        CALL GlobalSegMap_init (GSMapROMS1, start, length, 0,           &
     &                        OCN_COMM_WORLD, OCNid)
      ELSE IF (ng.eq.2) THEN
        CALL GlobalSegMap_init (GSMapROMS2, start, length, 0,           &
     &                        OCN_COMM_WORLD, OCNid)
      END IF
#else
      CALL GlobalSegMap_init (GSMapROMS, start, length, 0,              &
     &                        OCN_COMM_WORLD, OCNid)
#endif
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
#ifdef MCT_INTERP_OC2WV
!
!  Determine start and lengths for roms domain decomposition.
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
      CALL GlobalSegMap_init (GSMapROMS, start, length, 0,              &
     &                        OCN_COMM_WORLD, OCNid)
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
!  Determine start and lengths for domain decomposition
!  of the wave model.
!
      Jsize=dst_grid_dims(1)
      IF (.not.allocated(start)) THEN
        allocate ( start(Jsize) )
      END IF
      IF (.not.allocated(length)) THEN
        allocate ( length(Jsize) )
      END IF
      jc=0
      DO j=0,Jsize-1
        jc=jc+1
        start (jc)=j*dst_grid_dims(2)+1
        length(jc)=dst_grid_dims(2)
      END DO

      CALL GlobalSegMap_init (GSMapSWAN, start, length, 0,              &
     &                        OCN_COMM_WORLD, OCNid)
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
! Create Waves sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
        call SparseMatrixPlus_init(W2OMatPlus, sMatW,                   &
     &                             GSMapSWAN, GSMapROMS, Xonly,         &
     &                             MyMaster, OCN_COMM_WORLD, OCNid)
        call SparseMatrix_clean(sMatW)
!
! Create Ocean sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
         call SparseMatrixPlus_init(O2WMatPlus, sMatO,                  &
     &                              GSMapROMS, GSMapSWAN, Xonly,        &
     &                              MyMaster, OCN_COMM_WORLD, OCNid)
        call SparseMatrix_clean(sMatO)
#endif
#ifdef MCT_INTERP_OC2WV
!
!  Initialize attribute vector holding the export data code strings of
!  the wave model. The Asize is the number of grid point on this
!  processor.
!
!     call GlobalSegMap_Ordpnts(xPrimeGSMap,MyRank,points)
      Asize=GlobalSegMap_lsize(GSMapSWAN, OCN_COMM_WORLD)
!     CALL AttrVect_init (wav2ocn_AV2, rList=TRIM(ExportList(Iwaves)),  &
!    &                    lsize=Asize)
      CALL AttrVect_init(wav2ocn_AV2,                                   &
     &     rList="DISSIP:HSIGN:RTP:SETUP:TMBOT:UBOT:DIR:WLEN:TM01:QB",  &
     &     lsize=Asize)
      CALL AttrVect_zero (wav2ocn_AV2)
!
      Asize=GlobalSegMap_lsize(GSMapROMS, OCN_COMM_WORLD)
!     CALL AttrVect_init (wav2ocn_AV, rList=TRIM(ExportList(Iwaves)),   &
!    &                    lsize=Asize)
      CALL AttrVect_init(wav2ocn_AV2,                                   &
     &     rList="DISSIP:HSIGN:RTP:SETUP:TMBOT:UBOT:DIR:WLEN:TM01:QB",  &
     &     lsize=Asize)
      CALL AttrVect_zero (wav2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      Asize=GlobalSegMap_lsize(GSmapSWAN, OCN_COMM_WORLD)
!     CALL AttrVect_init (ocn2wav_AV2, rList=TRIM(ExportList(Iocean)),  &
!    &                    lsize=Asize)
      CALL AttrVect_init (ocn2wav_AV2,                                  &
     &                    rList="DEPTH:WLEV:VELX:VELY:ZO",              &
     &                    lsize=Asize)
      CALL AttrVect_zero (ocn2wav_AV2)
!
      Asize=GlobalSegMap_lsize(GSmapROMS, OCN_COMM_WORLD)
!     CALL AttrVect_init (ocn2wav_AV, rList=TRIM(ExportList(Iocean)),   &
!    &                    lsize=Asize)
      CALL AttrVect_init (ocn2wav_AV,                                   &
     &                    rList="DEPTH:WLEV:VELX:VELY:ZO",              &
     &                    lsize=Asize)
      CALL AttrVect_zero (ocn2wav_AV)
!
!  Initialize a router to the wave model component.
!
      CALL Router_init (WAVid, GSMapSWAN, OCN_COMM_WORLD, ROMStoSWAN)
#else
# ifdef REFINED_GRID
!
!  Initialize attribute vector holding the export data code strings of
!  the wave model. The Asize is the number of grid point on this
!  processor.
!
      IF (ng.eq.1) THEN
        Asize=GlobalSegMap_lsize(GSMapROMS1, OCN_COMM_WORLD)
!       CALL AttrVect_init (wav2ocn_AV1, rList=TRIM(ExportList(Iwaves)),  &
!    &                      lsize=Asize)
        CALL AttrVect_init(wav2ocn_AV1,                                   &
     &  rList="DISSIP:HSIGN:RTP:SETUP:TMBOT:UBOT:DIR:WLEN:TM01:QB",       &
     &  lsize=Asize)
        CALL AttrVect_zero (wav2ocn_AV1)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
!       CALL AttrVect_init (ocn2wav_AV1, rList=TRIM(ExportList(Iocean)),  &
!    &                    lsize=Asize)
        CALL AttrVect_init (ocn2wav_AV1,                                  &
     &                    rList="DEPTH:WLEV:VELX:VELY:ZO",                &
     &                    lsize=Asize)
        CALL AttrVect_zero (ocn2wav_AV1)
!
!  Initialize a router to the wave model component.
!
        CALL Router_init (WAVid, GSMapROMS1, OCN_COMM_WORLD, ROMStoSWAN1)
      ELSE IF (ng.eq.2) THEN
        Asize=GlobalSegMap_lsize(GSMapROMS2, OCN_COMM_WORLD)
!       CALL AttrVect_init (wav2ocn_AV2, rList=TRIM(ExportList(Iwaves)),  &
!    &                      lsize=Asize)
        CALL AttrVect_init(wav2ocn_AV2,                                   &
     &  rList="DISSIP:HSIGN:RTP:SETUP:TMBOT:UBOT:DIR:WLEN:TM01:QB",       &
     &  lsize=Asize)
        CALL AttrVect_zero (wav2ocn_AV2)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
!       CALL AttrVect_init (ocn2wav_AV2, rList=TRIM(ExportList(Iocean)),  &
!    &                    lsize=Asize)
        CALL AttrVect_init (ocn2wav_AV2,                                  &
     &                    rList="DEPTH:WLEV:VELX:VELY:ZO",                &
     &                    lsize=Asize)
        CALL AttrVect_zero (ocn2wav_AV2)
!
!  Initialize a router to the wave model component.
!
        CALL Router_init (WAVid, GSMapROMS2, OCN_COMM_WORLD, ROMStoSWAN2)
      END IF
# else
!
!  Initialize attribute vector holding the export data code strings of
!  the wave model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GSMapROMS, OCN_COMM_WORLD)
!     CALL AttrVect_init (wav2ocn_AV, rList=TRIM(ExportList(Iwaves)),     &
!    &                    lsize=Asize)
      CALL AttrVect_init(wav2ocn_AV,                                      &
     &  rList="DISSIP:HSIGN:RTP:SETUP:TMBOT:UBOT:DIR:WLEN:TM01:QB",       &
     &  lsize=Asize)
      CALL AttrVect_zero (wav2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
!     CALL AttrVect_init (ocn2wav_AV, rList=TRIM(ExportList(Iocean)),     &
!    &                    lsize=Asize)
      CALL AttrVect_init (ocn2wav_AV,                                     &
     &                    rList="DEPTH:WLEV:VELX:VELY:ZO",                &
     &                    lsize=Asize)
      CALL AttrVect_zero (ocn2wav_AV)
!
!  Initialize a router to the wave model component.
!
      CALL Router_init (WAVid, GSMapROMS, OCN_COMM_WORLD, ROMStoSWAN)
# endif
#endif
#ifdef REFINED_GRID
      deallocate ( wavids )
      deallocate ( ocnids )
#endif
      RETURN
      END SUBROUTINE initialize_ocn2wav_coupling

      SUBROUTINE ocn2wav_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine acquires the coupling data streams between waves       !
!  and ocean models.  Currently,  the following data streams are       !
!  coded:                                                              !
!                                                                      !
!     (...) SWAN units                                                 !
!     [...] ROMS units                                                 !
!                                                                      !
!  Fields imported from SWAN model:                                    !
!                                                                      !
!     * Wave direction (degrees), [radians]                            !
!     * Significant wave height (m), [m]                               !
!     * Average wave length (m), [m]                                   !
!     * Surface wave relative peak period (s), [s]                     !
!     * Bottom wave period (s), [s]                                    !
!     * Percent of breakig waves (nondimensional), [nondimensional]    !
!     * Wave energy dissipation (W/m2), [m3/s3]                        !
!     * Wave bottom orbital velocity (m/s), [m/s]                      !
!                                                                      !
!  Fields exported to SWAN model:                                      !
!                                                                      !
!     * Bathymetry, bottom elevation (m), [m]                          !
!     * Free-surface, water surface elevation (m), [m]                 !
!     * Depth integrated u-momentum (m/s), [m/s]                       !
!     * Depth integrated v-momentum (m/s), [m/s]                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocn2wav_coupling_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif

      RETURN
      END SUBROUTINE ocn2wav_coupling
!
!***********************************************************************
      SUBROUTINE ocn2wav_coupling_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
      USE mod_sediment
!
      USE distribute_mod, ONLY : mp_reduce
      USE ROMS_import_mod, ONLY : ROMS_import2d
      USE ROMS_export_mod, ONLY : ROMS_export2d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrT, IendT, JstrT, JendT
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, Iimport, Iexport, MyError
      integer :: gtype, i, id, ifield, ij, j, k, status

      real(r8), parameter ::  Lwave_min = 1.0_r8
      real(r8), parameter ::  Lwave_max = 500.0_r8

      real(r8) :: add_offset, scale
      real(r8) :: RecvTime, SendTime, buffer(2), wtime(2)

      real(r8) :: my_wtime, cff, waven

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: AA
!     real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: uavg
!     real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: vavg

      real(r8), pointer :: A(:)

      character (len=3 ), dimension(2) :: op_handle
      character (len=40) :: code
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
      IF (WESTERN_EDGE) THEN
        IstrR=BOUNDS(ng)%Istr(tile)-1
      ELSE
        IstrR=BOUNDS(ng)%Istr(tile)
      END IF
      IF (EASTERN_EDGE) THEN
        IendR=BOUNDS(ng)%Iend(tile)+1
      ELSE
        IendR=BOUNDS(ng)%Iend(tile)
      END IF
      IF (SOUTHERN_EDGE) THEN
        JstrR=BOUNDS(ng)%Jstr(tile)-1
      ELSE
        JstrR=BOUNDS(ng)%Jstr(tile)
      END IF
      IF (NORTHERN_EDGE) THEN
        JendR=BOUNDS(ng)%Jend(tile)+1
      ELSE
        JendR=BOUNDS(ng)%Jend(tile)
      END IF
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
#ifdef REFINED_GRID
      IF (ng.eq.1) THEN
        Asize=GlobalSegMap_lsize (GSMapROMS1, OCN_COMM_WORLD)
      ELSE IF (ng.eq.2) THEN
        Asize=GlobalSegMap_lsize (GSMapROMS2, OCN_COMM_WORLD)
      END IF
#else
      Asize=GlobalSegMap_lsize (GSMapROMS, OCN_COMM_WORLD)
#endif
      allocate ( A(Asize) )
      A=0.0_r8
!
!  Initialize coupling wait time clocks.
!
      RecvTime=0.0_r8
      SendTime=0.0_r8
!
!-----------------------------------------------------------------------
!  Import fields from wave model (SWAN) to ocean model (ROMS).
!  Currently, both waves and ocean model grids are the same.
!  We need to revisit this logic to allow interpolation.
!-----------------------------------------------------------------------
!
!  Schedule receiving fields from wave model.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      buffer(1)=my_wtime(wtime)
#ifdef MCT_INTERP_OC2WV
      CALL MCT_Recv (wav2ocn_AV2, ROMStoSWAN, MyError)
      CALL MCT_MatVecMul(wav2ocn_AV2, W2OMatPlus, wav2ocn_AV)
#else
# ifdef REFINED_GRID
      CALL AttrVect_init (wav2ocn_AV, rList=TRIM(ExportList(Iwaves)),   &
     &                    lsize=Asize)
      CALL AttrVect_zero (wav2ocn_AV)
      IF (ng.eq.1) THEN
        CALL MCT_Recv (wav2ocn_AV1, ROMStoSWAN1, MyError)
        CALL AttrVect_copy(wav2ocn_AV1,wav2ocn_AV)
      ELSE IF (ng.eq.2) THEN
        CALL MCT_Recv (wav2ocn_AV2, ROMStoSWAN2, MyError)
        CALL AttrVect_copy(wav2ocn_AV2,wav2ocn_AV)
      END IF
# else
      CALL MCT_Recv (wav2ocn_AV, ROMStoSWAN, MyError)
# endif
#endif
      RecvTime=RecvTime+my_wtime(wtime)-buffer(1)
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'wave model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      END IF
!
!  Receive fields from wave model.
!
!
!  Set ramp coefficient.
!
!!    ramp=MIN((tdays(ng)-dstart)*4.0_r8,1.0_r8)
      ramp=1.0_r8
!
!  Wave dissipation.
!
      CALL AttrVect_exportRAttr (FrWAVToOCNAV, "DISSIP", A, Asize)
      ij=0
      cff=1.0_r8/rho0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          Wave_dissip(i,j)=MAX(0.0_r8,A(ij)*ramp)*cff
        END DO
      END DO
!
!  Wave height.
!
      CALL AttrVect_exportRAttr (FrWAVToOCNAV, "HSIGN", A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          Hwave(i,j)=MAX(0.0_r8,A(ij)*ramp)
        END DO
      END DO
!
!  Surface wave period.
!
      CALL AttrVect_exportRAttr(FrWAVToOCNAV, "RTP", A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          Pwave_top(i,j)=MAX(0.0_r8,A(ij))
        END DO
      END DO
!
!  Bottom wave period.
!
      CALL AttrVect_exportRAttr (FrWAVToOCNAV, "TMBOT", A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          Pwave_bot(i,j)=MAX(0.0_r8,A(ij))
        END DO
      END DO
!
!  Bottom orbital velocity (m/s).
!
      CALL AttrVect_exportRAttr(FrWAVToOCNAV, "UBOT", A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          Ub_swan(i,j)=MAX(0.0_r8,A(ij)*ramp)
        END DO
      END DO
!
!  Wave direction (radians).
!
      CALL AttrVect_exportRAttr (FrWAVToOCNAV, "DIR", A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          Dwave(i,j)=MAX(0.0_r8,A(ij))*deg2rad
        END DO
      END DO
!
!  Wave length (m).
!
      CALL AttrVect_exportRAttr (FrWAVToOCNAV, "WLEN", A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
            Lwave(i,j)=MAX(1.0_r8,A(ij))
!            IF (Lwave(i,j).eq.1.0_r8/0.0_r8) THEN
!              Lwave(i,j)=Lwave_max
!            END IF
            LWave(i,j)=MIN(Lwave_max,A(ij))
        END DO
      END DO

#ifdef SVENDSEN_ROLLER
!
!  Percent wave breaking.
!  
      CALL AttrVect_exportRAttr (FrWAVToOCNAV, "QB", A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          Wave_break(i,j)=MAX(0.0_r8,A(ij))
        END DO
      END DO
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Wave_dissip)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Hwave)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Pwave_top)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Pwave_bot)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Ub_swan)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Dwave)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Lwave)
# ifdef SVENDSEN_ROLLER
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Wave_break)
# endif
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile boundaries.
!-----------------------------------------------------------------------
!
      CALL mp_exchange2d (ng, iNLM, 4, tile,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Wave_dissip, Hwave, Dwave, Lwave)
      CALL mp_exchange2d (ng, iNLM, 3, tile,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Pwave_top, Pwave_bot, Ub_swan)
# ifdef SVENDSEN_ROLLER
      CALL mp_exchange2d (ng, iNLM, 1, tile,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    wave_break)
# endif
#endif
!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to wave (SWAN) model.
!-----------------------------------------------------------------------
# ifdef REFINED_GRID
      CALL AttrVect_init (ocn2wav_AV, rList=TRIM(ExportList(Iocean)),   &
     &                    lsize=Asize)
      CALL AttrVect_zero (ocn2wav_AV)
# endif
!
!  Schedule sending fields to the wave model.
!
!
!  Depth (bathymetry).
!
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          A(ij)=h(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ocn2wav_AV, "DEPTH", A, Asize)
!
!  Water level (free-surface).
!
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          A(ij)=zeta(i,j,knew)
        END DO
      END DO
      CALL AttrVect_importRAttr (ocn2wav_AV, "WLEV", A, Asize)
!
!  Vertically-integrated U-velocity at RHO-points.
!
#  ifdef REFINED_GRID
        DO j=JstrT,JendT
          DO i=IstrTU+1,IendTU
            ubar_rho(i,j)=0.5_r8*                                       &
     &                  (u(i,j,N(ng),nstp)+u(i+1,j,N(ng),nstp))
          END DO
        END DO
        IF (WESTERN_EDGE) THEN
          DO j=JstrT,JendT
            ubar_rho(IstrT,j)=ubar_rho(IstrT+1,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=JstrT,JendT
            ubar_rho(IendT,j)=ubar_rho(IendT-1,j)
          END DO
        END IF
#  else
        DO j=JstrR,JendR
          DO i=Istr,Iend
            ubar_rho(i,j)=0.5_r8*                                       &
     &                  (u(i,j,N(ng),nstp)+u(i+1,j,N(ng),nstp))
          END DO
        END DO
        IF (WESTERN_EDGE) THEN
          DO j=Jstr,Jend
            ubar_rho(Istr-1,j)=ubar_rho(Istr,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=Jstr,Jend
            ubar_rho(Iend+1,j)=ubar_rho(Iend,j)
          END DO
        END IF
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          ubar_rho(Istr-1,Jstr-1)=0.5_r8*(ubar_rho(Istr  ,Jstr-1)+      &
     &                                 ubar_rho(Istr-1,Jstr  ))
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          ubar_rho(Iend+1,Jstr-1)=0.5_r8*(ubar_rho(Iend  ,Jstr-1)+      &
     &                                 ubar_rho(Iend+1,Jstr  ))
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          ubar_rho(Istr-1,Jend+1)=0.5_r8*(ubar_rho(Istr-1,Jend  )+      &
     &                                 ubar_rho(Istr  ,Jend+1))
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          ubar_rho(Iend+1,Jend+1)=0.5_r8*(ubar_rho(Iend+1,Jend  )+      &
     &                                 ubar_rho(Iend  ,Jend+1))
        END IF
#  endif
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
#ifdef UV_CONST
          A(ij)=0.0_r8
#else
          A(ij)=ubar_rho(i,j)
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (ocn2wav_AV, "VELX", A, Asize)
!
!  Vertically-integrated V-velocity at RHO-points.
!
#  ifdef REFINED_GRID
        DO j=JstrTV+1,JendTV
          DO i=IstrT,IendT
            vbar_rho(i,j)=0.5_r8*                                       &
     &                  (v(i,j,N(ng),nstp)+v(i,j+1,N(ng),nstp))
          END DO
        END DO
        IF (NORTHERN_EDGE) THEN
          DO i=IstrT,IendT
            vbar_rho(i,JendT)=vbar_rho(i,JendT-1)
          END DO
        END IF
        IF (SOUTHERN_EDGE) THEN
          DO i=IstrT,IendT
            vbar_rho(i,JstrT)=vbar_rho(i,JstrT+1)
          END DO
        END IF
#  else
        DO j=Jstr,Jend
          DO i=IstrR,IendR
            vbar_rho(i,j)=0.5_r8*                                       &
     &                  (v(i,j,N(ng),nstp)+v(i,j+1,N(ng),nstp))
          END DO
        END DO
        IF (NORTHERN_EDGE) THEN
          DO i=Istr,Iend
            vbar_rho(i,Jend+1)=vbar_rho(i,Jend)
          END DO
        END IF
        IF (SOUTHERN_EDGE) THEN
          DO i=Istr,Iend
            vbar_rho(i,Jstr-1)=vbar_rho(i,Jstr)
          END DO
        END IF
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          vbar_rho(Istr-1,Jstr-1)=0.5_r8*(vbar_rho(Istr  ,Jstr-1)+            &
     &                                 vbar_rho(Istr-1,Jstr  ))
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          vbar_rho(Iend+1,Jstr-1)=0.5_r8*(vbar_rho(Iend  ,Jstr-1)+            &
     &                                 vbar_rho(Iend+1,Jstr  ))
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          vbar_rho(Istr-1,Jend+1)=0.5_r8*(vbar_rho(Istr-1,Jend  )+            &
     &                                 vbar_rho(Istr  ,Jend+1))
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          vbar_rho(Iend+1,Jend+1)=0.5_r8*(vbar_rho(Iend+1,Jend  )+            &
     &                                 vbar_rho(Iend  ,Jend+1))
        END IF
#  endif
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
#ifdef UV_CONST
          A(ij)=0.0_r8
#else
          A(ij)=vbar_rho(i,j)
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (ocn2wav_AV, "VELY", A, Asize)
!
!  bottom roughness.
!
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          A(ij)=h(i,j)
#ifdef BBL_MODEL
                AA(ij)=MAX(0.0001_r8,                                  &
     &                      OCEAN(ng)%bottom(i,j,izNik)*30.0_r8)
#else
                AA(ij)=MAX(0.0001_r8,rdrg2(ng))
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (ocn2wav_AV, "ZO", A, Asize)
!
!  Send ocean fields to wave model.
!
#ifdef MCT_INTERP_OC2WV
        call MCT_MatVecMul(ocn2wav_AV,O2WMatPlus,ocn2wav_AV2)
        CALL MCT_Send (ocn2wav_AV2, ROMStoSWAN, MyError)
#else
# ifdef REFINED_GRID
      IF (ng.eq.1) THEN
        CALL AttrVect_copy(ocn2wav_AV,ocn2wav_AV1)
        CALL MCT_Send (ocn2wav_AV1, ROMStoSWAN1, MyError)
      ELSE IF (ng.eq.2) THEN
        CALL AttrVect_copy(ocn2wav_AV,ocn2wav_AV2)
        CALL MCT_Send (ocn2wav_AV2, ROMStoSWAN2, MyError)
      END IF
!
      CALL AttrVect_clean (ocn2wav_AV, MyError)
      CALL AttrVect_clean (wav2ocn_AV, MyError)
# else
        CALL MCT_Send (ocn2wav_AV, ROMStoSWAN, MyError)
# endif
#endif
        IF (MyError.ne.0) THEN
          IF (Master) THEN
            WRITE (stdout,20) 'wave model, MyError = ', MyError
          END IF
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (A)
!
 10   FORMAT (' OCN2WAV_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2WAV_COUPLING - error while sending fields to: ',    &
     &        a, i4)
 30   FORMAT (6x,'OCN2WAV   - (', i2.2, ') imported and (', i2.2,       &
     &        ') exported fields,', t62, 't = ', a,/, 16x,              &
     &        '- ROMS coupling exchanges wait clock (s):',/, 19x,       &
     &        '(Recv= ', 1p,e14.8,0p, ' Send= ', 1p,e14.8,0p,')')
 40   FORMAT (16x,'- ',a,a,                                             &
     &        /,19x,'(Min= ',1p,e15.8,0p,' Max= ',1p,e15.8,0p,')')

      RETURN
      END SUBROUTINE ocn2wav_coupling_tile

      SUBROUTINE finalize_ocn2wav_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and wave models coupling data streams.  !
!                                                                       !
!========================================================================
!
!  Local variable declarations.
!
      integer :: MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
#ifdef REFINED_GRID
!      IF (ng.eq.1) THEN
        CALL Router_clean (ROMStoSWAN1, MyError)
        CALL AttrVect_clean (ocn2wav_AV1, MyError)
        CALL GlobalSegMap_clean (GSMapROMS1, MyError)
!      ELSE IF (ng.eq.2) THEN
        CALL Router_clean (ROMStoSWAN2, MyError)
        CALL AttrVect_clean (ocn2wav_AV2, MyError)
        CALL GlobalSegMap_clean (GSMapROMS2, MyError)
!      END IF
#else
      CALL Router_clean (ROMStoSWAN, MyError)
      CALL AttrVect_clean (ocn2wav_AV, MyError)
      CALL GlobalSegMap_clean (GSMapROMS, MyError)
#endif
      RETURN

      END SUBROUTINE finalize_ocn2wav_coupling
