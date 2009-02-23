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
      CALL AttrVect_init (wav2ocn_AV2, rList=TRIM(ExportList(Iwaves)),   &
     &                    lsize=Asize)
      CALL AttrVect_zero (wav2ocn_AV2)
!
      Asize=GlobalSegMap_lsize(GSMapROMS, OCN_COMM_WORLD)
      CALL AttrVect_init (wav2ocn_AV, rList=TRIM(ExportList(Iwaves)),   &
     &                    lsize=Asize)
      CALL AttrVect_zero (wav2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      Asize=GlobalSegMap_lsize(GSmapSWAN, OCN_COMM_WORLD)
      CALL AttrVect_init (ocn2wav_AV2, rList=TRIM(ExportList(Iocean)),   &
     &                    lsize=Asize)
      CALL AttrVect_zero (ocn2wav_AV2)
!
      Asize=GlobalSegMap_lsize(GSmapROMS, OCN_COMM_WORLD)
      CALL AttrVect_init (ocn2wav_AV, rList=TRIM(ExportList(Iocean)),   &
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
        CALL AttrVect_init (wav2ocn_AV1, rList=TRIM(ExportList(Iwaves)),  &
     &                      lsize=Asize)
        CALL AttrVect_zero (wav2ocn_AV1)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
        CALL AttrVect_init (ocn2wav_AV1, rList=TRIM(ExportList(Iocean)),  &
     &                    lsize=Asize)
        CALL AttrVect_zero (ocn2wav_AV1)
!
!  Initialize a router to the wave model component.
!
        CALL Router_init (WAVid, GSMapROMS1, OCN_COMM_WORLD, ROMStoSWAN1)
      ELSE IF (ng.eq.2) THEN
        Asize=GlobalSegMap_lsize(GSMapROMS2, OCN_COMM_WORLD)
        CALL AttrVect_init (wav2ocn_AV2, rList=TRIM(ExportList(Iwaves)),  &
     &                      lsize=Asize)
        CALL AttrVect_zero (wav2ocn_AV2)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
        CALL AttrVect_init (ocn2wav_AV2, rList=TRIM(ExportList(Iocean)),  &
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
      CALL AttrVect_init (wav2ocn_AV, rList=TRIM(ExportList(Iwaves)),   &
     &                    lsize=Asize)
      CALL AttrVect_zero (wav2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      CALL AttrVect_init (ocn2wav_AV, rList=TRIM(ExportList(Iocean)),   &
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
      Iimport=0
      DO ifield=1,Nimport(Iocean)
        id=ImportID(Iocean)%val(ifield)
        code=ADJUSTL(Fields(id)%code)
        gtype=Fields(id)%GridType
        scale=Fields(id)%scale
        add_offset=Fields(id)%AddOffset

        SELECT CASE (TRIM(code))

          CASE ('Wdir')                   ! wave direction

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=deg2rad                 ! degress to radians
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrT, IendT, JstrT, JendT,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Dwave,                       &
     &                          status)

          CASE ('Wamp')                   ! significant wave hight

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8                  ! m
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrT, IendT, JstrT, JendT,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Hwave,                       &
     &                          status)

          CASE ('Wlen')                   ! wave length

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(Lwave_min,A(i))
              IF (A(i).eq.1.0_r8/0.0_r8) THEN
                A(i)=Lwave_max
              END IF
            END DO
            scale=1.0_r8                  ! m
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrT, IendT, JstrT, JendT,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Lwave,                       &
     &                          status)

          CASE ('Wptop')                  ! peak surface wave period

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrT, IendT, JstrT, JendT,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Pwave_top,                   &
     &                          status)

          CASE ('Wpbot')                  ! mean bottom wave period

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrT, IendT, JstrT, JendT,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Pwave_bot,                   &
     &                          status)

          CASE ('Wdiss')                  ! wave dissipation

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8/rho0
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrT, IendT, JstrT, JendT,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Wave_dissip,                 &
     &                          status)

          CASE ('Wubot')                  ! bottom orbital velocity

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8                  ! m/s
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrT, IendT, JstrT, JendT,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Ub_swan,                     &
     &                          status)

#ifdef SVENDSEN_ROLLER

          CASE ('Wbrk')                   ! percent wave breaking

            CALL AttrVect_exportRAttr (wav2ocn_AV, TRIM(code), A, Asize)
            Iimport=Iimport+1
            DO i=1,Asize
              A(i)=MAX(0.0_r8,A(i))
            END DO
            scale=1.0_r8
            add_offset=0.0_r8
            CALL ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          Asize, A,                               &
     &                          IstrT, IendT, JstrT, JendT,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fields(id)%ImpMin, Fields(id)%ImpMax,   &
     &                          FORCES(ng)%Wave_break,                  &
     &                          status)
#endif
        END SELECT
      END DO
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
      Iexport=0
      DO ifield=1,Nexport(Iocean)
        id=ExportID(Iocean)%val(ifield)
        code=ADJUSTL(Fields(id)%code)
        gtype=Fields(id)%GridType
        scale=Fields(id)%scale
        add_offset=Fields(id)%AddOffset

        SELECT CASE (TRIM(code))

          CASE ('bath')                   ! bathymetry (depth)

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          GRID(ng)%h,                             &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

          CASE ('SSH')                    ! free-surface (water level)

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng)%zeta(:,:,KOUT),               &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

          CASE ('Ubar')                   ! 2D U-momentum
!
!           Vertically integrate the velicties according to             !
!           Dingemans (1997) Water wave propagation over uneven bottoms,!
!           Advanced series on Ocean Engineering, p.86 eq 2.144.        !
!
!            DO j=JstrR,JendR
!              DO i=Istr,Iend
!                waven=2.0_r8*pi/MAX(Lwave_min,FORCES(ng)%Lwave(i,j))
!                cff=0.0_r8
!                DO k=1,N(ng)
!                  cff=cff+OCEAN(ng)%u(i,j,k,NOUT)*cosh(2.0_r8*waven*    &
!     &                (GRID(ng)%z_r(i,j,k)+GRID(ng)%h(i,j)))*           &
!     &                GRID(ng)%Hz(i,j,k)
!                END DO
!                uavg(i,j)=2.0_r8*waven*cff/                             &
!     &               sinh(2.0_r8*waven*GRID(ng)%h(i,j))
!                uavg(i,j)=0.5_r8*(OCEAN(ng)%u(i,j,N(ng),NOUT)+          &
!     &                            OCEAN(ng)%u(i+1,j,N(ng),NOUT))
!              END DO
!            END DO
!     &                          uavg(LBi:UBi,LBj:UBj),                  &

!    &                          uavg(LBi:UBi,LBj:UBj),                  &
            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
#ifdef SOLVE3D
     &                          OCEAN(ng)%u(:,:,N(ng),NOUT),            &
#else
     &                          OCEAN(ng)%ubar(:,:,KOUT),               &
#endif
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

          CASE ('Vbar')                   ! 2D V-momentum

!            DO j=Jstr,Jend+1
!              DO i=IstrR,IendR
!                waven=2.0_r8*pi/FORCES(ng)%Lwave(i,j)
!                cff=0.0_r8
!                DO k=1,N(ng)
!                  cff=cff+OCEAN(ng)%v(i,j,k,NOUT)*cosh(2.0_r8*waven*    &
!     &                (GRID(ng)%z_r(i,j,k)+GRID(ng)%h(i,j)))*           &
!     &                GRID(ng)%Hz(i,j,k)
!                END DO
!                vavg(i,j)=2.0_r8*waven*cff/                             &
!     &                    sinh(2.0_r8*waven*GRID(ng)%h(i,j))
!              END DO
!            END DO
!     &                          vavg(LBi:UBi,LBj:UBj),                  &

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
#ifdef SOLVE3D
     &                          OCEAN(ng)%v(:,:,N(ng),NOUT),            &
#else
     &                          OCEAN(ng)%vbar(:,:,KOUT),               &
#endif
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

          CASE ('ZO')                   ! bottom roughness

            DO j=JstrT,JendT
              DO i=IstrT,IendT
#ifdef BBL_MODEL
                AA(i,j)=MAX(0.0001_r8,                                  &
     &                      OCEAN(ng)%bottom(i,j,izNik)*30.0_r8)
#else
                AA(i,j)=MAX(0.0001_r8,rdrg2(ng))
#endif
              END DO
            END DO
            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          AA(LBi:UBi,LBj:UBj),                    &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
            CALL AttrVect_importRAttr (ocn2wav_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

        END SELECT
      END DO
!
!  Send ocean fields to wave model.
!
      IF (Iexport.gt.0) THEN
        buffer(2)=my_wtime(wtime)
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
        SendTime=SendTime+my_wtime(wtime)-buffer(2)
        IF (MyError.ne.0) THEN
          IF (Master) THEN
            WRITE (stdout,20) 'wave model, MyError = ', MyError
          END IF
          exit_flag=2
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Report.
!-----------------------------------------------------------------------
!
      IF (Nthreads(Iocean).gt.1) THEN
        buffer(1)=RecvTime
        buffer(2)=SendTime
        op_handle(1)='SUM'
        op_handle(2)='SUM'
        CALL mp_reduce (ng, iNLM, 2, buffer, op_handle)
        RecvTime=buffer(1)
        SendTime=buffer(2)
      END IF
      IF (Master.and.((Iimport.gt.0).or.(Iexport.gt.0))) THEN
        WRITE (stdout,30) Iimport, Iexport, time_code(ng),              &
     &                    RecvTime, SendTime
        IF (Lreport) THEN
          DO ifield=1,Nimport(Iocean)
            id=ImportID(Iocean)%val(ifield)
            WRITE (stdout,40) 'ROMS Import: ',TRIM(fields(id)%name),    &
     &                        Fields(id)%ImpMin, Fields(id)%ImpMax
          END DO
          DO ifield=1,Nexport(Iocean)
            id=ExportID(Iocean)%val(ifield)
            WRITE (stdout,40) 'ROMS Export: ',TRIM(fields(id)%name),    &
     &                        Fields(id)%ExpMin, Fields(id)%ExpMax
          END DO
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
