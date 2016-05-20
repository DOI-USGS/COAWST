/*
** svn $Id: mct_roms_swan.h 756 2008-09-14 20:18:28Z jcwarner $
***************************************************** John C. Warner ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group      Hernan G. Arango  **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** These routines couple ROMS/TOMS to the SWAN wave model using       **
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
      USE mct_coupler_params
      USE mod_kinds
      USE mod_scalars
      USE mod_iounits
#ifdef MCT_INTERP_OC2WV
      USE mod_coupler_iounits
#endif
#if defined VEGETATION && defined VEG_SWAN_COUPLING
      USE mod_vegetation
      USE mod_vegarr 
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
      integer :: i, ic, iw, j, jc, nprocs
      integer :: nRows, nCols, num_sparse_elems
      integer :: cid, cad
      integer, allocatable :: length(:)
      integer, allocatable :: start(:)
!      integer, dimension(2) :: src_grid_dims, dst_grid_dims
      character (len=70)    :: nc_name
      character (len=20)    :: to_add
      character (len=120)   :: wostring
      character (len=120)   :: owstring
      real(r8) :: cff
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
!  Establish MCT communicator.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, nprocs, MyError)
!
      IF (ng.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Nocn_grids))
        ALLOCATE(AttrVect_G(Nocn_grids))
#ifdef MCT_INTERP_OC2WV
        ALLOCATE(SMPlus_W(Nocn_grids,Nwav_grids))
        ALLOCATE(AV2_W(Nocn_grids,Nwav_grids))
        ALLOCATE(GSMapInterp_W(Nocn_grids,Nwav_grids))
#endif
      END IF
!
!  Initialize MCT coupled model registry.
!
      OCNid=ocnids(ng)
      IF (Nocn_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      OCN_COMM_WORLD,myids=ocnids)
      ELSE
        CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD,                &
     &                      OCN_COMM_WORLD,OCNid)
      END IF
!
!  Determine the part of the grid we are working on and develop
!  local segment of the global map.
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
        start (jc)=(j)*(Lm(ng)+2)+IstrR+1
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
!
#ifdef MCT_INTERP_OC2WV
!  If ocean grid and wave grids are different sizes, then
!  develop sparse matrices for interpolation.
!
  35  FORMAT(a3,i1,a7,i1,a11)
      DO iw=1,Nwav_grids
!
!  Prepare sparse matrices. First work on waves to ocean.
!
        IF (Myrank.eq.MyMaster) THEN
          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'wav',iw,'_to_ocn',ng,'_weights.nc'
          ELSE
            nc_name=W2Oname(iw,ng)
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

          call SparseMatrix_init(sMatW,nRows,nCols,num_sparse_elems)
          call SparseMatrix_importGRowInd(sMatW, sparse_rows,           &
     &                                    size(sparse_rows))
          call SparseMatrix_importGColInd(sMatW, sparse_cols,           &
     &                                    size(sparse_cols))
          call SparseMatrix_importMatrixElts(sMatW, sparse_weights,     &
     &                                       size(sparse_weights))
!
! Deallocate arrays.
!
          deallocate ( sparse_rows )
          deallocate ( sparse_cols )
          deallocate ( sparse_weights )
          deallocate ( dst_grid_imask )
!
!  Prepare sparse matrices. Second work on ocean to waves.
!
          IF (scrip_opt.eq.1) THEN
            write(nc_name,35) 'ocn',ng,'_to_wav',iw,'_weights.nc'
          ELSE
            nc_name=O2Wname(ng,iw)
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

          call SparseMatrix_init(sMatO,nRows,nCols,num_sparse_elems)
          call SparseMatrix_importGRowInd(sMatO, sparse_rows,           &
     &                                    size(sparse_rows))
          call SparseMatrix_importGColInd(sMatO, sparse_cols,           &
     &                                    size(sparse_cols))
          call SparseMatrix_importMatrixElts(sMatO, sparse_weights,     &
     &                                       size(sparse_weights))
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
!  Determine start and lengths for domain decomposition
!  of the wave model.
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
        CALL GlobalSegMap_init (GSMapInterp_W(ng,iw)%GSMapSWAN,         &
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
! Create Waves sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
        call SparseMatrixPlus_init(SMPlus_W(ng,iw)%W2OMatPlus, sMatW,   &
     &                             GSMapInterp_W(ng,iw)%GSMapSWAN,      &
     &                             GlobalSegMap_G(ng)%GSMapROMS,        &
     &                             Xonly,MyMaster,OCN_COMM_WORLD,OCNid)
        call SparseMatrix_clean(sMatW)
!
! Create Ocean sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
         call SparseMatrixPlus_init(SMPlus_W(ng,iw)%O2WMatPlus, sMatO,  &
     &                              GlobalSegMap_G(ng)%GSMapROMS,       &
     &                              GSMapInterp_W(ng,iw)%GSMapSWAN,     &
     &                              Xonly,MyMaster,OCN_COMM_WORLD,OCNid)
        call SparseMatrix_clean(sMatO)
      END DO
#endif
!
!  Initialize the list of fields from the wave model.
!
      cad=LEN(wostring)
      DO i=1,cad
        wostring(i:i)=''
      END DO
      cid=1
!
      to_add='DISBOT'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DISSURF'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DISWCAP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':HSIGN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':RTP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':TMBOT'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':UBOT'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':DIRE'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad

      to_add=':DIRN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WLEN'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WLENP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':QB'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WDSPR'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WQP'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
#if defined WAVES_OCEAN && defined WEC_VF && \
    defined BOTTOM_STREAMING && defined VEGETATION &&  \
    defined VEG_SWAN_COUPLING && defined VEG_STREAMING
      to_add=':DISVEG'
      cad=LEN_TRIM(to_add)
      write(wostring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif

!
!  Finalize and remove trailing spaces from the wostring
!  for the rlist.
!
      cad=LEN_TRIM(wostring)
      wostring=wostring(1:cad)
!
!  Initialize attribute vector holding the export data of
!  the wav model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &                         OCN_COMM_WORLD)
      CALL AttrVect_init(AttrVect_G(ng)%wav2ocn_AV,                     &
     &                   rList=TRIM(wostring),lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%wav2ocn_AV)
!
!  Initialize attribute vector that contain the data strings from
!  the ocean model.
!
      cad=LEN(owstring)
      DO i=1,cad
        owstring(i:i)=''
      END DO
      cid=1
!
      to_add='DEPTH'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':WLEV'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VELX'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VELY'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':ZO'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#if defined VEGETATION && defined VEG_SWAN_COUPLING 
!
      to_add=':VEGDENS'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
#ifdef ICE_MODEL
!
      to_add=':SEAICE'
      cad=LEN_TRIM(to_add)
      write(owstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#endif
!
!  Finalize and remove trailing spaces from the owstring
!  for the rlist.
!
      cad=LEN_TRIM(owstring)
      owstring=owstring(1:cad)
!
      CALL AttrVect_init(AttrVect_G(ng)%ocn2wav_AV,                     &
     &                   rList=TRIM(owstring),lsize=Asize)
      CALL AttrVect_zero(AttrVect_G(ng)%ocn2wav_AV)
!
#ifdef MCT_INTERP_OC2WV
      DO iw=1,Nwav_grids
        WAVid=wavids(iw)
!
!  Initialize attribute vectors that are used for the sparse
!  interpolation. These vectors are _AV2 and are the same vales
!  as _AV but different size grids.
!
        Asize=GlobalSegMap_lsize(GSMapInterp_W(ng,iw)%GSMapSWAN,        &
     &                           OCN_COMM_WORLD)
        CALL AttrVect_init(AV2_W(ng,iw)%wav2ocn_AV2,                    &
     &                     rList=TRIM(wostring),lsize=Asize)
        CALL AttrVect_zero(AV2_W(ng,iw)%wav2ocn_AV2)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
        CALL AttrVect_init(AV2_W(ng,iw)%ocn2wav_AV2,                    &
     &                     rList=TRIM(owstring),lsize=Asize)
        CALL AttrVect_zero(AV2_W(ng,iw)%ocn2wav_AV2)
      END DO
#endif
      RETURN
      END SUBROUTINE initialize_ocn2wav_coupling

      SUBROUTINE initialize_ocn2wav_routers (tile)
!
!=======================================================================
!                                                                      !
!  Initialize ocean and wave models coupling stream.  This is the      !
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
      integer :: ng, iw
!
!-----------------------------------------------------------------------
!  Establish MCT router.
!-----------------------------------------------------------------------
!
      ALLOCATE(Router_W(Nocn_grids,Nwav_grids))
!
!  Initialize routers to the wave model component.
!
      DO ng=1,Nocn_grids
        DO iw=1,Nwav_grids
          WAVid=wavids(iw)
#ifdef MCT_INTERP_OC2WV
          CALL Router_init (WAVid, GSMapInterp_W(ng,iw)%GSMapSWAN,      &
     &                      OCN_COMM_WORLD, Router_W(ng,iw)%ROMStoSWAN)
#else
          CALL Router_init (WAVid, GlobalSegMap_G(ng)%GSMapROMS,        &
     &                      OCN_COMM_WORLD, Router_W(ng,iw)%ROMStoSWAN)
#endif
        END DO
      END DO
      RETURN
      END SUBROUTINE initialize_ocn2wav_routers

      SUBROUTINE ocn2wav_coupling (ng, iw, tile)
!
!=======================================================================
!                                                                      !
!  This routine acquires the coupling data streams between waves       !
!  and ocean models.                                                   !
!  coded:                                                              !
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
      integer, intent(in) :: ng, iw, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocn2wav_coupling_tile (ng, iw, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif
      RETURN
      END SUBROUTINE ocn2wav_coupling
!
!***********************************************************************
      SUBROUTINE ocn2wav_coupling_tile (ng, iw, tile,                   &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS)
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
      USE mod_sedbed
      USE mod_sediment
#ifdef UV_KIRBY
      USE mod_coupling
#endif
#if defined VEGETATION && defined VEG_SWAN_COUPLING
      USE mod_vegetation
      USE mod_vegarr 
#endif
!
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
      integer, intent(in) :: ng, iw, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Asize, MyError, Tag
      integer :: gtype, i, id, ifield, ij, j, k, status
#if defined VEGETATION && defined VEG_SWAN_COUPLING
      integer :: iveg
#endif	

      real(r8), parameter ::  Lwave_min = 1.0_r8
      real(r8), parameter ::  Lwave_max = 500.0_r8

      real(r8) :: add_offset, scale
      real(r8) :: cff, ramp
      real(r8) :: cff1, cff2, cff3, cff4, kwn, prof, u_cff, v_cff

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ubar_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vbar_rho

      real(r8), pointer :: A(:)
#ifdef MCT_INTERP_OC2WV
      integer, pointer :: indices(:)
#endif
!
#include "set_bounds.h"
!
!  Modify ranges to allow full exchange of fields for periodic applications.
!
      IF (EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IstrR=Istr-1
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IendR=Iend+1
        END IF
      END IF
      IF (NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          JstrR=Jstr-1
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          JendR=Jend+1
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,           &
     &      OCN_COMM_WORLD)
      allocate ( A(Asize) )
      A=0.0_r8
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to wave (SWAN) model.
!-----------------------------------------------------------------------
!
!  Depth (bathymetry).
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=GRID(ng)%h(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "DEPTH",    &
     &                           A, Asize)
!
!  Water level (free-surface).
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
#ifdef ZETA_CONST
          A(ij)=0.0_r8
#else
          A(ij)=OCEAN(ng)%zeta(i,j,knew(ng))
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "WLEV",     &
     &                           A, Asize)
!
!  U-velocity at RHO-points.
!
      DO j=JstrR,JendR
        DO i=Istr,Iend
#ifdef SOLVE3D
# ifdef UV_KIRBY
!
! Compute the coupling current according to Kirby and Chen (1989).
!
          kwn=2.0_r8*pi/FORCES(ng)%Lwave(i,j)
          prof=GRID(ng)%h(i,j)+COUPLING(ng)%Zt_avg1(i,j)
          cff1=0.0_r8
          cff2=2.0_r8*kwn*prof
          IF (cff2.lt.700.0_r8) THEN
            cff2=2.0_r8*kwn
          ELSE
            cff2=700.0_r8/prof
          ENDIF
          cff3=0.0_r8
          DO k=1,N(ng)
            u_cff=0.5_r8*(OCEAN(ng)%u(i,  j,k,NOUT)+                    &
     &                    OCEAN(ng)%u(i+1,j,k,NOUT))
            cff4=cosh(cff2*(GRID(ng)%h(i,j)+GRID(ng)%z_r(i,j,k)))*      &
     &           GRID(ng)%Hz(i,j,k)
            cff1=cff1+cff4*u_cff
            cff3=cff3+cff4
          END DO
          ubar_rho(i,j)=cff1/cff3
# else
          ubar_rho(i,j)=0.5_r8*(OCEAN(ng)%u(i,  j,N(ng),NOUT)+          &
     &                          OCEAN(ng)%u(i+1,j,N(ng),NOUT))
# endif
#else
          ubar_rho(i,j)=0.5_r8*(OCEAN(ng)%ubar(i,  j,KOUT)+             &
     &                          OCEAN(ng)%ubar(i+1,j,KOUT))
#endif
        END DO
      END DO
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=JstrR,JendR
          ubar_rho(IstrR,j)=ubar_rho(IstrR+1,j)
        END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=JstrR,JendR
          ubar_rho(IendR,j)=ubar_rho(IendR-1,j)
        END DO
      END IF
#ifdef SOLVE3D
# ifdef UV_KIRBY
        DO j=JstrR,JendR
          DO i=IstrR,IendR
             OCEAN(ng)%uwave(i,j)=ubar_rho(i,j)
          ENDDO
        ENDDO
# endif
#endif
!
!  V-velocity at RHO-points.
!
      DO j=Jstr,Jend
        DO i=IstrR,IendR
#ifdef SOLVE3D
# ifdef UV_KIRBY
!
! Compute the coupling current according to Kirby and Chen (1989).
!
          kwn=2.0_r8*pi/FORCES(ng)%Lwave(i,j)
          prof=GRID(ng)%h(i,j)+COUPLING(ng)%Zt_avg1(i,j)
          cff1=0.0_r8
          cff2=2.0_r8*kwn*prof
          IF (cff2.lt.700.0_r8) THEN
            cff2=2.0_r8*kwn
          ELSE
            cff2=700.0_r8/prof
          ENDIF
          cff3=0.0_r8
          DO k=1,N(ng)
             v_cff=0.5_r8*(OCEAN(ng)%v(i,  j,k,NOUT)+                   &
     &                     OCEAN(ng)%v(i,j+1,k,NOUT))
             cff4=cosh(cff2*(GRID(ng)%h(i,j)+GRID(ng)%z_r(i,j,k)))*     &
     &            GRID(ng)%Hz(i,j,k)
             cff1=cff1+cff4*v_cff
             cff3=cff3+cff4
          END DO
          vbar_rho(i,j)=cff1/cff3
# else
          vbar_rho(i,j)=0.5_r8*(OCEAN(ng)%v(i,j  ,N(ng),NOUT)+          &
     &                          OCEAN(ng)%v(i,j+1,N(ng),NOUT))
# endif
#else
          vbar_rho(i,j)=0.5_r8*(OCEAN(ng)%vbar(i,j  ,KOUT)+             &
     &                          OCEAN(ng)%vbar(i,j+1,KOUT))
#endif
        END DO
      END DO
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=IstrR,IendR
          vbar_rho(i,JendR)=vbar_rho(i,JendR-1)
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=IstrR,IendR
          vbar_rho(i,JstrR)=vbar_rho(i,JstrR+1)
        END DO
      END IF
#ifdef SOLVE3D
# ifdef UV_KIRBY
      DO j=JstrR,JendR
        DO i=IstrR,Iend
          OCEAN(ng)%vwave(i,j)=vbar_rho(i,j)
        ENDDO
      ENDDO
# endif
#endif
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
#ifdef UV_CONST
          A(ij)=0.0_r8
#else
# ifdef CURVGRID
!
! Rotate velocity to be East positive.
!
          cff1=ubar_rho(i,j)*GRID(ng)%CosAngler(i,j)-                   &
     &         vbar_rho(i,j)*GRID(ng)%SinAngler(i,j)
# else
          cff1=ubar_rho(i,j)
# endif
          A(ij)=cff1
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "VELX",     &
     &                           A, Asize)
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
#ifdef UV_CONST
          A(ij)=0.0_r8
#else
# ifdef CURVGRID
!
! Rotate velocity to be North positive.
!
          cff1=ubar_rho(i,j)*GRID(ng)%SinAngler(i,j)+                   &
     &         vbar_rho(i,j)*GRID(ng)%CosAngler(i,j)
# else
          cff1=vbar_rho(i,j)
# endif
          A(ij)=cff1
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "VELY",     &
     &                           A, Asize)
!
!  bottom roughness.
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
#ifdef BBL_MODEL
# if defined VEGETATION && defined VEG_SWAN_COUPLING
!         Specify 0.0015 to be consistent with Z0_min in ROMS.
          A(ij)=MAX(0.0015_r8, SEDBED(ng)%bottom(i,j,izNik)*30.0_r8)
# else
!         Specify this to be Madsen 0.05 minimum.
          A(ij)=MAX(0.05_r8, SEDBED(ng)%bottom(i,j,izNik)*30.0_r8)
# endif 
#else
!               This value will be replaced by the value entered in the
!               SWAN INPUT file. See SWAN/Src/waves_coupler.F.
                A(ij)=0.05_r8
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "ZO",       &
     &                           A, Asize)
#if defined VEGETATION && defined VEG_SWAN_COUPLING 
!
!  Equivalent Plant density.
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=0.0
          DO iveg=1,NVEG
            cff=VEG(ng)%plant(i,j,iveg,pdens)+cff
          END DO
          A(ij)=cff/NVEG
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "VEGDENS",  &
     &                           A, Asize)
#endif
#ifdef ICE_MODEL
!
!  sea ice.
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=0.0_r8
        ! A(ij)=MAX(ICE(ng)%sfw(i,j,1) ! NEED TO PROVIDE CORRECT VAR
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "SEAICE",   &
     &                           A, Asize)
#endif
!
!  Send ocean fields to wave model.
!
      Tag=ng*100+0*10+iw
#ifdef MCT_INTERP_OC2WV
      CALL MCT_MatVecMul(AttrVect_G(ng)%ocn2wav_AV,                     &
     &                   SMPlus_W(ng,iw)%O2WMatPlus,                    &
     &                   AV2_W(ng,iw)%ocn2wav_AV2)
      CALL MCT_isend (AV2_W(ng,iw)%ocn2wav_AV2,                         &
     &                Router_W(ng,iw)%ROMStoSWAN, Tag)
#else
      CALL MCT_isend (AttrVect_G(ng)%ocn2wav_AV,                        &
     &                Router_W(ng,iw)%ROMStoSWAN, Tag)
#endif
      CALL MCT_waits (Router_W(ng,iw)%ROMStoSWAN)
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) 'wave model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
        WRITE (stdout,36) ' ** ROMS grid ',ng,                          &
     &                    ' sent data to SWAN grid ',iw
 36     FORMAT (a14,i2,a24,i2)
        END IF
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (A)
#ifdef MCT_INTERP_OC2WV
      if (associated (indices)) then
        deallocate (indices)
      endif
#endif
!
 10   FORMAT (' OCN2WAV_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2WAV_COUPLING - error while sending fields to: ',    &
     &        a, i4)
      RETURN
      END SUBROUTINE ocn2wav_coupling_tile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ocnfwav_coupling (ng, iw, tile)
!
!=======================================================================
!                                                                      !
!  This routine acquires the coupling data streams between waves       !
!  and ocean models.                                                   !
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
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, iw, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocnfwav_coupling_tile (ng, iw, tile,                         &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif
      RETURN
      END SUBROUTINE ocnfwav_coupling
!
!***********************************************************************
      SUBROUTINE ocnfwav_coupling_tile (ng, iw, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mct_coupler_params
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits
      USE mod_sedbed
      USE mod_sediment
#if defined VEGETATION && defined VEG_SWAN_COUPLING
      USE mod_vegetation
      USE mod_vegarr
#endif

#ifdef UV_KIRBY
      USE mod_coupling
#endif
!
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
      integer, intent(in) :: ng, iw, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Asize, MyError, Tag
      integer :: gtype, i, id, ifield, ij, j, k, status

      real(r8), parameter ::  Lwave_min = 1.0_r8
      real(r8), parameter ::  Lwave_max = 500.0_r8
      real(r8), parameter ::  Large = 1.0E+20_r8

      real(r8) :: add_offset, scale
      real(r8) :: cff, fac, ramp
      real(r8) :: cff1, cff2, cff3, cff4, kwn, prof, u_cff, v_cff

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ubar_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vbar_rho
      real(r8), dimension(2) :: range

      real(r8), pointer :: A(:)
      real(r8), pointer :: A1(:)
#ifdef MCT_INTERP_OC2WV
      integer, pointer :: indices(:)
#endif
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
!  Modify ranges to allow full exchange of fields for periodic applications.
!
      IF (EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IstrR=Istr-1
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IendR=Iend+1
        END IF
      END IF
      IF (NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          JstrR=Jstr-1
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          JendR=Jend+1
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,           &
     &      OCN_COMM_WORLD)
      allocate ( A(Asize) )
      allocate ( A1(Asize) )
      A=0.0_r8
      A1=0.0_r8
!
!-----------------------------------------------------------------------
!  Import fields from wave model (SWAN) to ocean model (ROMS).
!-----------------------------------------------------------------------
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      Tag=ng*100+0*10+iw
#ifdef MCT_INTERP_OC2WV
      CALL MCT_irecv (AV2_W(ng,iw)%wav2ocn_AV2,                         &
     &                Router_W(ng,iw)%ROMStoSWAN, Tag)
!     Wait to make sure the SWAN data has arrived.
      CALL MCT_waitr (AV2_W(ng,iw)%wav2ocn_AV2,                         &
     &                Router_W(ng,iw)%ROMStoSWAN)
      CALL MCT_MatVecMul(AV2_W(ng,iw)%wav2ocn_AV2,                      &
     &                   SMPlus_W(ng,iw)%W2OMatPlus,                    &
     &                   AttrVect_G(ng)%wav2ocn_AV)
#else
      CALL MCT_irecv (AttrVect_G(ng)%wav2ocn_AV,                        &
     &                Router_W(ng,iw)%ROMStoSWAN, Tag)
!     Wait to make sure the SWAN data has arrived.
      CALL MCT_waitr (AttrVect_G(ng)%wav2ocn_AV,                        &
     &                Router_W(ng,iw)%ROMStoSWAN)
#endif
!
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'wave model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
        WRITE (stdout,36) ' ** ROMS grid ',ng,                          &
     &                    ' recv data from SWAN grid ',iw
 36     FORMAT (a14,i2,a26,i2)
        END IF
      END IF
!
!  Set ramp coefficient.
!
      ramp=1.0_r8
!
!  Receive fields from wave model.
 40         FORMAT (a36,1x,2(1pe14.6))
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
!
!  Wave dissipation due to bottom friction.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DISBOT",   &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      fac=1.0_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)*fac
          IF (iw.eq.1) THEN
            FORCES(ng)%Dissip_fric(i,j)=cff
          ELSE
            FORCES(ng)%Dissip_fric(i,j)=FORCES(ng)%Dissip_fric(i,j)+    &
     &                                  cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DISBOT  (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
#endif
#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
!
!  Wave dissipation due to surface breaking.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DISSURF",  &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      fac=1.0_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)*fac
          IF (iw.eq.1) THEN
            FORCES(ng)%Dissip_break(i,j)=cff
          ELSE
            FORCES(ng)%Dissip_break(i,j)=FORCES(ng)%Dissip_break(i,j)+  &
     &                                   cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DISSURF (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
!
!  Wave dissipation due to white capping.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DISWCAP",  &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      fac=1.0_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)*fac
          IF (iw.eq.1) THEN
            FORCES(ng)%Dissip_wcap(i,j)=cff
          ELSE
            FORCES(ng)%Dissip_wcap(i,j)=FORCES(ng)%Dissip_wcap(i,j)+    &
     &                                  cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DISWCAP (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
#endif
!
!  Wave height.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "HSIGN",    &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)
          IF (iw.eq.1) THEN
            FORCES(ng)%Hwave(i,j)=cff
          ELSE
            FORCES(ng)%Hwave(i,j)=FORCES(ng)%Hwave(i,j)+                &
     &                            cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max HSIGN   (m):     ',        &
     &                    range(1),range(2)
      END IF
!
!  Surface peak wave period.
!
      CALL AttrVect_exportRAttr(AttrVect_G(ng)%wav2ocn_AV, "RTP",       &
     &                          A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij))
          IF (iw.eq.1) THEN
            FORCES(ng)%Pwave_top(i,j)=cff
          ELSE
            FORCES(ng)%Pwave_top(i,j)=FORCES(ng)%Pwave_top(i,j)+        &
     &                                cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max RTP     (s):     ',        &
     &                    range(1),range(2)
      END IF
!
!  Bottom mean wave period.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "TMBOT",    &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij))
          IF (iw.eq.1) THEN
            FORCES(ng)%Pwave_bot(i,j)=cff
          ELSE
            FORCES(ng)%Pwave_bot(i,j)=FORCES(ng)%Pwave_bot(i,j)+        &
     &                                cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max TMBOT   (s):     ',        &
     &                    range(1),range(2)
      END IF
#if defined BBL_MODEL
!
!  Bottom orbital velocity (m/s).
!
      CALL AttrVect_exportRAttr(AttrVect_G(ng)%wav2ocn_AV, "UBOT",      &
     &                          A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)
          IF (iw.eq.1) THEN
            FORCES(ng)%Uwave_rms(i,j)=cff
          ELSE
            FORCES(ng)%Uwave_rms(i,j)=FORCES(ng)%Uwave_rms(i,j)+        &
     &                                cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max UBOT    (ms-1):  ',        &
     &                    range(1),range(2)
      END IF
#endif
!
!  Wave direction (radians).
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DIRE",     &
     &                           A, Asize)
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DIRN",     &
     &                           A1, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=ATAN2(A(ij),A1(ij))
          IF (cff.lt.0.0_r8) cff=cff+2.0_r8*pi
          IF (iw.eq.1) THEN
            FORCES(ng)%Dwave(i,j)=cff
          ELSE
            FORCES(ng)%Dwave(i,j)=FORCES(ng)%Dwave(i,j)+cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DIR     (deg):   ',        &
     &                    range(1),range(2)
      END IF
!
!  Wave length (m).
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "WLEN",     &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MIN(Lwave_max,MAX(1.0_r8,A(ij)))
          IF (iw.eq.1) THEN
            FORCES(ng)%Lwave(i,j)=cff
          ELSE
            FORCES(ng)%Lwave(i,j)=FORCES(ng)%Lwave(i,j)+                &
     &                            cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max WLEN    (m):     ',        &
     &                    range(1),range(2)
      END IF
#ifdef WAVES_LENGTHP
!
!  Wave length (m).
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "WLENP",    &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MIN(Lwave_max,MAX(1.0_r8,A(ij)))
          IF (iw.eq.1) THEN
            FORCES(ng)%Lwavep(i,j)=cff
          ELSE
            FORCES(ng)%Lwavep(i,j)=FORCES(ng)%Lwavep(i,j)+              &
     &                             cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max WLENP   (m):     ',        &
     &                    range(1),range(2)
      END IF
#endif
!
#ifdef ROLLER_SVENDSEN
!
!  Percent wave breaking.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "QB",       &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij))
          IF (iw.eq.1) THEN
            FORCES(ng)%Wave_break(i,j)=cff
          ELSE
            FORCES(ng)%Wave_break(i,j)=FORCES(ng)%Wave_break(i,j)+      &
     &                                 cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max QB      (%):     ',        &
     &                    range(1),range(2)
      END IF
#endif
#ifdef WAVES_DSPR
!
!  wave directional spreading
!  
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "WDSPR",    &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij))
          IF (iw.eq.1) THEN
            FORCES(ng)%Wave_ds(i,j)=cff
          ELSE
            FORCES(ng)%Wave_ds(i,j)=FORCES(ng)%Wave_ds(i,j)+            &
     &                              cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max WDSPR   (deg):   ',        &
     &                    range(1),range(2)
      END IF
!
!  wave spectrum peakedness
!  
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "WQP",      &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij))
          IF (iw.eq.1) THEN
            FORCES(ng)%Wave_qp(i,j)=cff
          ELSE
            FORCES(ng)%Wave_qp(i,j)=FORCES(ng)%Wave_qp(i,j)+            &
     &                              cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max WQP     (-):     ',        &
     &                    range(1),range(2)
      END IF
#endif 
!
#if defined WAVES_OCEAN && defined WEC_VF && \
    defined BOTTOM_STREAMING && defined VEGETATION &&  \
    defined VEG_SWAN_COUPLING && defined VEG_STREAMING
!
!  Wave dissipation due to vegetation.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DISVEG",   &
     &                           A, Asize)
      range(1)= Large
      range(2)=-Large
      ij=0
      fac=1.0_r8/rho0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          cff=MAX(0.0_r8,A(ij)*ramp)*fac
          IF (iw.eq.1) THEN 
            VEG(ng)%Dissip_veg(i,j)=cff
          ELSE
            VEG(ng)%Dissip_veg(i,j)=VEG(ng)%Dissip_veg(i,j)+            &    
     &                              cff
          END IF
          range(1)=MIN(range(1),cff)
          range(2)=MAX(range(2),cff)
        END DO
      END DO
# ifdef DISTRIBUTE
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
# endif
      IF (Myrank.eq.MyMaster) THEN
        write(stdout,40) 'SWANtoROMS Min/Max DISVEG  (Wm-2):  ',        &
     &                    range(1),range(2)
      END IF
#endif
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dissip_break)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dissip_fric)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dissip_wcap)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Hwave)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Pwave_top)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Pwave_bot)
# if defined BBL_MODEL
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Uwave_rms)
# endif
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dwave)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Lwave)
# ifdef WAVES_LENGTHP
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Lwavep)
# endif
# ifdef ROLLER_SVENDSEN
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Wave_break)
# endif
# ifdef WAVES_DSPR
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Wave_ds)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Wave_qp)
# endif
      END IF
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile boundaries.
!-----------------------------------------------------------------------
!
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Dissip_break,                      &
     &                    FORCES(ng)%Dissip_fric,                       &
     &                    FORCES(ng)%Dissip_wcap)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Hwave, FORCES(ng)%Pwave_top,       &
     &                    FORCES(ng)%Pwave_bot)
# if defined BBL_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Uwave_rms)
# endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Dwave, FORCES(ng)%Lwave)
# ifdef WAVES_LENGTHP
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Lwavep)
# endif
# ifdef ROLLER_SVENDSEN
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Wave_break)
# endif
# ifdef WAVES_DSPR
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    FORCES(ng)%Wave_ds, FORCES(ng)%Wave_qp)
# endif
#endif
!
!  Deallocate communication arrays.
!
      deallocate (A)
      deallocate (A1)
#ifdef MCT_INTERP_OC2WV
      if (associated (indices)) then
        deallocate (indices)
      endif
#endif
!
 10   FORMAT (' OCN2WAV_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2WAV_COUPLING - error while sending fields to: ',    &
     &        a, i4)
      RETURN
      END SUBROUTINE ocnfwav_coupling_tile

      SUBROUTINE finalize_ocn2wav_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and wave models coupling data streams.  !
!                                                                       !
!========================================================================
      USE mod_scalars
      USE mct_coupler_params
!
!  Local variable declarations.
!
      integer :: ng, iw, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      deallocate ( wavids )
      deallocate ( ocnids )
      DO ng=1,Nocn_grids
!       CALL GlobalSegMap_clean (GlobalSegMap_G(ng)%GSMapROMS,          &
!    &                           MyError)
        DO iw=1,Nwav_grids
          CALL Router_clean (Router_W(ng,iw)%ROMStoSWAN, MyError)
          CALL AttrVect_clean (AttrVect_G(ng)%ocn2wav_AV, MyError)
        END DO
      END DO
      RETURN

      END SUBROUTINE finalize_ocn2wav_coupling
