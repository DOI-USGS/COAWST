/*
** svn $Id: mct_roms_swan.h 756 2008-09-14 20:18:28Z jcwarner $
***************************************************** John C. Warner ***
** Copyright (c) 2002-2010 The ROMS/TOMS Group      Hernan G. Arango  **
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
      IF (ng.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Ngrids))
        ALLOCATE(AttrVect_G(Ngrids))
        ALLOCATE(Router_G(Ngrids))
#ifdef MCT_INTERP_OC2AT
        ALLOCATE(SMPlus_G(Ngrids))
#endif
      END IF
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
     &                    OCN_COMM_WORLD,OCNid, myids=ocnids)
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
        nc_name=W2Oname(ng)
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
        nc_name=O2Wname(ng)
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
#ifdef MCT_INTERP_OC2WV
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

      CALL GlobalSegMap_init (GlobalSegMap_G(ng)%GSMapSWAN,             &
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
! Create Waves sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
        call SparseMatrixPlus_init(W2OMatPlus, sMatW,                   &
     &                             GlobalSegMap_G(ng)%GSMapSWAN,        &
     &                             GlobalSegMap_G(ng)%GSMapROMS,        &
     &                             Xonly,MyMaster,OCN_COMM_WORLD,OCNid)
        call SparseMatrix_clean(sMatW)
!
! Create Ocean sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
         call SparseMatrixPlus_init(O2WMatPlus, sMatO,                  &
     &                              GlobalSegMap_G(ng)%GSMapROMS,       &
     &                              GlobalSegMap_G(ng)%GSMapSWAN,       &
     &                              Xonly,MyMaster,OCN_COMM_WORLD,OCNid)
        call SparseMatrix_clean(sMatO)
#endif
#ifdef MCT_INTERP_OC2WV
!
!  Initialize attribute vector holding the export data code strings of
!  the wave model. The Asize is the number of grid point on this
!  processor.
!
!     call GlobalSegMap_Ordpnts(xPrimeGSMap,MyRank,points)
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapSWAN,            &
     &      OCN_COMM_WORLD)
      CALL AttrVect_init(wav2ocn_AV2,                                   &
     &     rList="DISSIP:HSIGN:RTP:TMBOT:UBOT:DIR:WLEN:QB",             &
     &     lsize=Asize)
      CALL AttrVect_zero (wav2ocn_AV2)
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &      OCN_COMM_WORLD)
      CALL AttrVect_init(AttrVect_G(ng)%wav2ocn_AV,                     &
     &     rList="DISSIP:HSIGN:RTP:TMBOT:UBOT:DIR:WLEN:QB",             &
     &     lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%wav2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapSWAN,            &
     &      OCN_COMM_WORLD)
      CALL AttrVect_init (ocn2wav_AV2,                                  &
     &                    rList="DEPTH:WLEV:VELX:VELY:ZO",              &
     &                    lsize=Asize)
      CALL AttrVect_zero (ocn2wav_AV2)
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &      OCN_COMM_WORLD)
      CALL AttrVect_init (AttrVect_G(ng)%ocn2wav_AV,                    &
     &                    rList="DEPTH:WLEV:VELX:VELY:ZO",              &
     &                    lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%ocn2wav_AV)
!
!  Initialize a router to the wave model component.
!
      CALL Router_init (WAVid, GlobalSegMap_G(ng)%GSMapSWAN,            &
     &                  OCN_COMM_WORLD, Router_G(ng)%ROMStoSWAN)
#else
!
!  Initialize attribute vector holding the export data code strings of
!  the wave model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &                         OCN_COMM_WORLD)
      CALL AttrVect_init(AttrVect_G(ng)%wav2ocn_AV,                     &
     &  rList="DISSIP:HSIGN:RTP:TMBOT:UBOT:DIR:WLEN:QB",                &
     &  lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%wav2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      CALL AttrVect_init (AttrVect_G(ng)%ocn2wav_AV,                    &
     &                    rList="DEPTH:WLEV:VELX:VELY:ZO",              &
     &                    lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%ocn2wav_AV)
!
!  Initialize a router to the wave model component.
!
      CALL Router_init (WAVid, GlobalSegMap_G(ng)%GSMapROMS,  &
     &                  OCN_COMM_WORLD, Router_G(ng)%ROMStoSWAN)
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
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif
      RETURN
      END SUBROUTINE ocn2wav_coupling
!
!***********************************************************************
      SUBROUTINE ocn2wav_coupling_tile (ng, tile,                       &
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
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrT, IendT, JstrT, JendT
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: IstrTU, IendTU, JstrTV, JendTV
      integer :: Asize, Iimport, Iexport, MyError
      integer :: gtype, i, id, ifield, ij, j, k, status

      real(r8), parameter ::  Lwave_min = 1.0_r8
      real(r8), parameter ::  Lwave_max = 500.0_r8

      real(r8) :: add_offset, scale
      real(r8) :: cff, ramp
      real(r8) :: cff1, cff2, cff3, cff4, kwn, prof, u_cff, v_cff

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: ubar_rho
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: vbar_rho

      real(r8), pointer :: A(:)
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
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
      IstrTU=BOUNDS(ng)%IstrTU(tile)
      IendTU=BOUNDS(ng)%IendTU(tile)
      JstrTV=BOUNDS(ng)%JstrTV(tile)
      JendTV=BOUNDS(ng)%JendTV(tile)
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
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,           &
     &      OCN_COMM_WORLD)
      allocate ( A(Asize) )
      A=0.0_r8
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
#ifdef MCT_INTERP_OC2WV
      CALL MCT_Recv (wav2ocn_AV2, Router_G(ng)%ROMStoSWAN, MyError)
      CALL MCT_MatVecMul(wav2ocn_AV2, W2OMatPlus,                       &
     &                   AttrVect_G(ng)%wav2ocn_AV)
#else
      CALL MCT_Recv (AttrVect_G(ng)%wav2ocn_AV, Router_G(ng)%ROMStoSWAN,&
     &               MyError)
#endif
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'wave model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
        WRITE (stdout,36) ' ** ROMS grid ',ng,' recv data from SWAN'
 36     FORMAT (a14,i2,a20)
        END IF
      END IF
!
!  Set ramp coefficient.
!
!!    ramp=MIN((tdays(ng)-dstart)*4.0_r8,1.0_r8)
      ramp=1.0_r8
!
!  Receive fields from wave model.
!
!  Wave dissipation.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DISSIP",   &
     &                           A, Asize)
      ij=0
      cff=1.0_r8/rho0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Wave_dissip(i,j)=MAX(0.0_r8,A(ij)*ramp)*cff
        END DO
      END DO
!
!  Wave height.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "HSIGN",    &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Hwave(i,j)=MAX(0.0_r8,A(ij)*ramp)
        END DO
      END DO
!
!  Surface peak wave period.
!
      CALL AttrVect_exportRAttr(AttrVect_G(ng)%wav2ocn_AV, "RTP",       &
     &                          A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Pwave_top(i,j)=MAX(0.0_r8,A(ij))
        END DO
      END DO
!
!  Bottom mean wave period.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "TMBOT",    &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Pwave_bot(i,j)=MAX(0.0_r8,A(ij))
        END DO
      END DO
!
!  Bottom orbital velocity (m/s).
!
      CALL AttrVect_exportRAttr(AttrVect_G(ng)%wav2ocn_AV, "UBOT",      &
     &                          A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Ub_swan(i,j)=MAX(0.0_r8,A(ij)*ramp)
        END DO
      END DO
!
!  Wave direction (radians).
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "DIR",      &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Dwave(i,j)=MAX(0.0_r8,A(ij))*deg2rad
        END DO
      END DO
!
!  Wave length (m).
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "WLEN",     &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Lwave(i,j)=MIN(Lwave_max,MAX(1.0_r8,A(ij)))
        END DO
      END DO
#ifdef SVENDSEN_ROLLER
!
!  Percent wave breaking.
!  
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%wav2ocn_AV, "QB",       &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Wave_break(i,j)=MAX(0.0_r8,A(ij))
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
     &                        FORCES(ng)%Wave_dissip)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Hwave)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Pwave_top)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Pwave_bot)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Ub_swan)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Dwave)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Lwave)
# ifdef SVENDSEN_ROLLER
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Wave_break)
# endif
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile boundaries.
!-----------------------------------------------------------------------
!
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%Wave_dissip, FORCES(ng)%Hwave,     &
     &                    FORCES(ng)%Dwave, FORCES(ng)%Lwave)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%Pwave_top, FORCES(ng)%Pwave_bot,   &
     &                    FORCES(ng)%Ub_swan)
# ifdef SVENDSEN_ROLLER
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%wave_break)
# endif
#endif
!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to wave (SWAN) model.
!-----------------------------------------------------------------------
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
          A(ij)=GRID(ng)%h(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "DEPTH",    &
     &                           A, Asize)
!
!  Water level (free-surface).
!
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
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
      DO j=JstrT,JendT
        DO i=IstrTU+1,IendTU
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
#ifdef SOLVE3D
# ifdef UV_KIRBY
        DO j=JstrT,JendT
          DO i=IstrT,IendT
             OCEAN(ng)%uwave(i,j)=ubar_rho(i,j)
          ENDDO
        ENDDO
# endif  
#endif 
!
!  V-velocity at RHO-points.
!
      DO j=JstrTV+1,JendTV
        DO i=IstrT,IendT
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
          vbar_rho(i,j)=0.5_r8*(OCEAN(ng)%v(i,j  ,KOUT)+                &
     &                          OCEAN(ng)%v(i,j+1,KOUT))
#endif
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
#ifdef SOLVE3D
# ifdef UV_KIRBY
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          OCEAN(ng)%vwave(i,j)=vbar_rho(i,j)
        ENDDO
      ENDDO
# endif
#endif
!
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
#ifdef UV_CONST
          A(ij)=0.0_r8
#else
!
! Rotate velocity to be East positive.
!
# ifdef CURVGRID
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
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
#ifdef UV_CONST
          A(ij)=0.0_r8
#else
!
! Rotate velocity to be North positive.
!
# ifdef CURVGRID
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
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
#ifdef BBL_MODEL
                A(ij)=MAX(0.0001_r8,                                    &
     &                    SEDBED(ng)%bottom(i,j,izNik)*30.0_r8)
#else
                A(ij)=MAX(0.0001_r8,rdrg2(ng))
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2wav_AV, "ZO",       &
     &                           A, Asize)
!
!  Send ocean fields to wave model.
!
#ifdef MCT_INTERP_OC2WV
      CALL MCT_MatVecMul(AttrVect_G(ng)%ocn2wav_AV,O2WMatPlus,          &
     &                   ocn2wav_AV2)
      CALL MCT_Send (ocn2wav_AV2, Router_G(ng)%ROMStoSWAN, MyError)
#else
      CALL MCT_Send (AttrVect_G(ng)%ocn2wav_AV, Router_G(ng)%ROMStoSWAN,&
     &               MyError)
#endif
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) 'wave model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
        WRITE (stdout,35) ' ** ROMS grid ',ng,' sent data to SWAN'
 35     FORMAT (a14,i2,a18)
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
      RETURN
      END SUBROUTINE ocn2wav_coupling_tile

      SUBROUTINE finalize_ocn2wav_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and wave models coupling data streams.  !
!                                                                       !
!========================================================================
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: ng, MyError
!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL Router_clean (Router_G(ng)%ROMStoSWAN, MyError)
        CALL AttrVect_clean (AttrVect_G(ng)%ocn2wav_AV, MyError)
        CALL GlobalSegMap_clean (GlobalSegMap_G(ng)%GSMapROMS,          &
     &                           MyError)
      END DO
      RETURN

      END SUBROUTINE finalize_ocn2wav_coupling
