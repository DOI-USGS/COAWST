/*
** svn $Id: mct_roms_wrf.h 756 2008-09-14 20:18:28Z jcwarner $
***************************************************** John C. Warner ***
** Copyright (c) 2002-2008 The ROMS/TOMS Group      Hernan G. Arango  **
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
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, Isize, Jsize, MyError
      integer :: i, ic, j, jc, nprocs
      integer :: nRows, nCols, num_sparse_elems

      integer, allocatable  :: length(:)
      integer, allocatable  :: start(:)
      integer, dimension(2) :: src_grid_dims, dst_grid_dims
      character (len=70)    :: nc_name
      character (len=70)   :: avstring
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
#if !defined WAVES_OCEAN
      IF (ng.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Ngrids))
        ALLOCATE(AttrVect_G(Ngrids))
        ALLOCATE(Router_G(Ngrids))
      END IF
!
!  Initialize MCT coupled model registry.
!
      CALL MCTWorld_init (N_mctmodels, MPI_COMM_WORLD, OCN_COMM_WORLD,  &
     &                    OCNid)
#endif
#ifdef MCT_INTERP_OC2AT
!
!  If ocean grid and atm grids are different sizes, then
!  develop sparse matrices for interpolation.
!
!!!!!!!!!!!!!!!!!!!!!!
! First work on atm to ocean.
!!!!!!!!!!!!!!!!!!!!!!
!
      IF (Myrank.eq.MyMaster) THEN
       nc_name=A2Oname(ng)
       call get_sparse_matrix (ng, nc_name, num_sparse_elems,           &
     &                          src_grid_dims, dst_grid_dims)
!
! Init the sparse matrix.
!
        nRows=dst_grid_dims(1)*dst_grid_dims(2)
        nCols=src_grid_dims(1)*src_grid_dims(2)
!
! Create sparse matrix.
!
        call SparseMatrix_init(sMatA,nRows,nCols,num_sparse_elems)
        call SparseMatrix_importGRowInd(sMatA, sparse_rows,             &
     &                                  size(sparse_rows))
        call SparseMatrix_importGColInd(sMatA, sparse_cols,             &
     &                                  size(sparse_cols))
        call SparseMatrix_importMatrixElts(sMatA, sparse_weights,       &
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
! Second work on ocean to atm.
!!!!!!!!!!!!!!!!!!!!!!
!
      IF (Myrank.eq.MyMaster) THEN
        nc_name=O2Aname(ng)
        call get_sparse_matrix (ng, nc_name, num_sparse_elems,          &
     &                          src_grid_dims, dst_grid_dims)
!
! Init the sparse matrix.
!
        nRows=dst_grid_dims(1)*dst_grid_dims(2)
        nCols=src_grid_dims(1)*src_grid_dims(2)
!
! Zero out the destination cells with masking.
!
        ic=1
        DO i=1,nRows
          DO ic=i*4-3,i*4
            sparse_weights(ic)=sparse_weights(ic)*                      &
     &                         REAL(dst_grid_imask(i),r8)
          END DO
        END DO
!
! Create sparse matrix.
!
        call SparseMatrix_init(sMatO,nRows,nCols,num_sparse_elems)
        call SparseMatrix_importGRowInd(sMatO, sparse_rows,             &
     &                                  size(sparse_rows))
        call SparseMatrix_importGColInd(sMatO, sparse_cols,             &
     &                                  size(sparse_cols))
        call SparseMatrix_importMatrixElts(sMatO, sparse_weights,       &
     &                                  size(sparse_weights))
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
!
      CALL mpi_bcast(dst_grid_dims, 2, MPI_INTEGER, MyMaster,           &
     &               OCN_COMM_WORLD, MyError)

#endif
#if !defined WAVES_OCEAN && !defined MCT_INTERP_OC2AT
!
!  Determine start and lengths for domain decomposition.
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
#endif
#ifdef MCT_INTERP_OC2AT
# if !defined WAVES_OCEAN 
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
# endif
!
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
      CALL GlobalSegMap_init (GSMapWRF, start, length, 0,              &
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
! Create ATM sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
        call SparseMatrixPlus_init(A2OMatPlus, sMatA,                  &
     &                             GSMapWRF,                           &
     &                             GlobalSegMap_G(ng)%GSMapROMS, Xonly,&
     &                             MyMaster, OCN_COMM_WORLD, OCNid)
        call SparseMatrix_clean(sMatA)
!
! Create Ocean sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
!
         call SparseMatrixPlus_init(O2AMatPlus, sMatO,                 &
     &                              GlobalSegMap_G(ng)%GSMapROMS,      &
     &                              GSMapWRF, Xonly,                   &
     &                              MyMaster, OCN_COMM_WORLD, OCNid)
        call SparseMatrix_clean(sMatO)
#endif
!
!  Initialize attribute vector holding the export data code strings of
!  the atmosphere model.
!
      avstring(1:4)='PSFC'
      avstring(5:9)=':RELH'
      avstring(10:12)=':T2'
      avstring(13:16)=':U10'
      avstring(17:20)=':V10'
      avstring(21:27)=':CLDFRA'
      avstring(28:32)=':RAIN'
      avstring(33:39)=':SWDOWN'
      avstring(40:43)=':GLW'
      avstring(44:51)=':USTRESS'
      avstring(52:59)=':VSTRESS'
      avstring(60:62)=':LH'
      avstring(63:66)=':HFX'
      avstring(67:70)=':GSW'
!
#ifdef MCT_INTERP_OC2AT
!
!  Initialize attribute vector holding the export data code strings of
!  the atm model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GSMapWRF, OCN_COMM_WORLD)
      CALL AttrVect_init (atm2ocn_AV2, rList=TRIM(avstring),lsize=Asize)
      CALL AttrVect_zero (atm2ocn_AV2)
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &                         OCN_COMM_WORLD)
      CALL AttrVect_init (AttrVect_G(ng)%atm2ocn_AV,                    &
     &                    rList=TRIM(avstring),lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%atm2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      Asize=GlobalSegMap_lsize(GSmapWRF, OCN_COMM_WORLD)
      CALL AttrVect_init (ocn2atm_AV2,rList="SST",lsize=Asize)
      CALL AttrVect_zero (ocn2atm_AV2)
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSmapROMS,            &
     &                         OCN_COMM_WORLD)
      CALL AttrVect_init (AttrVect_G(ng)%ocn2atm_AV,rList="SST",        &
     &                    lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%ocn2atm_AV)
!
!  Initialize a router to the wave model component.
!
      CALL Router_init (ATMid, GSMapWRF, OCN_COMM_WORLD,                &
     &                  Router_G(ng)%ROMStoWRF)
#else
!
!  Initialize attribute vector holding the export data code strings of
!  the atmosphere model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ng)%GSMapROMS,            &
     &                         OCN_COMM_WORLD)
      CALL AttrVect_init (AttrVect_G(ng)%atm2ocn_AV,                    &
     &                    rList=TRIM(avstring),lsize=Asize)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      CALL AttrVect_init (AttrVect_G(ng)%ocn2atm_AV, rList="SST",       &
     &                    lsize=Asize)
      CALL AttrVect_zero (AttrVect_G(ng)%ocn2atm_AV)
!
!  Initialize a router to the atmosphere model component.
!
      CALL Router_init (ATMid, GlobalSegMap_G(ng)%GSMapROMS,            &
     &                  OCN_COMM_WORLD, Router_G(ng)%ROMStoWRF)
#endif

      RETURN
      END SUBROUTINE initialize_ocn2atm_coupling

      SUBROUTINE ocn2atm_coupling (ng, tile)
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
!     * PSFC    Surface atmospheric pressure (Pa), [mb]                !
!     * RELH    Surface air relative humidity (percent), [fraction]    !
!     * T2      Surface (2 m) air temperature (Celsius), [Celsius]     !
!     * U10     Surface (10 m) U-wind speed (m/s), [m/s]               !
!     * V10     Surface (10 m) V-wind speed (m/s), [m/s]               !
!     * CLDFRA  Cloud fraction (percent/100), [percent/100]            !
!     * RAIN    Precipitation (m/s), [kg/m2/s]                         !
!     * SWDOWN  Shortwave radiation (Watts/m2), [Celsius m/s]          !
!     * GSW     Net shortwave radiation (Watts/m2), [Celsius m/s]      !
!     * GLW     Long wave raditaion (Watts/m2), [Celsius m/s]          !
!     * USTRESS Surface U-wind stress (Pa), [m2/s2]                    !
!     * VSTRESS Surface V-wind stress (Pa), [m2/s2]                    !
!     * LH      Latent heat flux (Watts/m2), [Celsius m/s]             !
!     * HFX     Sensible heat flux (Watts/m2), [Celsius m/s]           !
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
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 48)
#endif
      CALL ocn2atm_coupling_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 48)
#endif

      RETURN
      END SUBROUTINE ocn2atm_coupling
!
!***********************************************************************
      SUBROUTINE ocn2atm_coupling_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj)
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
#ifdef CURVGRID
      USE mod_grid
#endif
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
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrT, IendT, JstrT, JendT
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, Iimport, Iexport, MyError
      integer :: gtype, i, id, ifield, j, ij,  status

      real(r8) :: add_offset, cff, scale, ramp
      real(r8) :: RecvTime, SendTime, buffer(2), wtime(2)

      real(r8), pointer :: A(:)
      real(r8) :: BBR, cff1, cff2

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
      Asize=GlobalSegMap_lsize (GlobalSegMap_G(ng)%GSMapROMS,         &
     &                          OCN_COMM_WORLD)
      allocate ( A(Asize) )
      A=0.0_r8
!
!  Initialize coupling wait time clocks.
!
!-----------------------------------------------------------------------
!  Import fields from atmosphere model (WRF) to ocean model (ROMS).
!-----------------------------------------------------------------------
!
!  Receive fields from atmosphere model.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
#ifdef MCT_INTERP_OC2AT
      CALL MCT_Recv (atm2ocn_AV2, Router_G(ng)%ROMStoWRF, MyError)
      CALL MCT_MatVecMul(atm2ocn_AV2, A2OMatPlus,                       &
     &                   AttrVect_G(ng)%atm2ocn_AV)
#else
      CALL MCT_Recv (AttrVect_G(ng)%atm2ocn_AV, Router_G(ng)%ROMStoWRF, &
     &               MyError)
#endif
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'atmosphere model, MyError = ', MyError
        END IF
        CALL finalize_ocn2atm_coupling
      ELSE
        IF (Master) THEN
          WRITE (stdout,'a') 'ROMS recv Atm fields'
        END IF
      END IF
!
!  Set ramp coefficient.
!
!!    ramp=MIN((tdays(ng)-dstart)*4.0_r8,1.0_r8)
      ramp=1.0_r8
!
!  Receive fields from atmosphere model.
!
#if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Surface atmospheric pressure  (mb).
!  Need to scale PSFC in Pa to mb.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "PSFC",     &
     &                           A, Asize)
      ij=0
      cff=0.01_r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Pair(i,j)=A(ij)*cff
        END DO
      END DO
#endif
#if defined BULK_FLUXES || defined ECOSIM || \
   (defined SHORTWAVE && defined ANA_SRFLUX)
!
!  Surface air relative humidity (-)
!  Convert RELH from percent to fraction.
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "RELH",     &
     &                           A, Asize)
      ij=0
      cff=0.01_r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Hair(i,j)=A(ij)*cff
        END DO
      END DO
!
!  Surface 2m air temperature    (degC)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "T2",       &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Tair(i,j)=A(ij)
        END DO
      END DO
#endif
#if defined BULK_FLUXES || defined ECOSIM
!
!  U-Wind speed at 10 m          (m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "U10",      &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Uwind(i,j)=A(ij)
        END DO
      END DO
!
!  V-Wind speed at 10 m          (m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "V10",      &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%Vwind(i,j)=A(ij)
        END DO
      END DO
# ifdef CURVGRID
!
!  If input point surface winds or interpolated from coarse data, rotate
!  to curvilinear grid.
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
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
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%cloud(i,j)=A(ij)
        END DO
      END DO
#endif
#ifdef BULK_FLUXES
!
!  Precipitation                 (kg/m2/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "RAIN",     &
     &                           A, Asize)
      cff=rho0
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%rain(i,j)=A(ij)*cff
        END DO
      END DO
!
!  Short wave radiation          (Celsius m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "GSW",      &
     &                           A, Asize)
      cff=1.0_r8/(rho0*Cp)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%srflx(i,j)=A(ij)*cff
        END DO
      END DO
!
!  Long wave radiation          (Celsius m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "GLW",      &
     &                           A, Asize)
      cff=1.0_r8/(rho0*Cp)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          BBR=OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)+273.16_r8
          BBR=BBR*BBR*BBR*BBR
          BBR=0.97_r8*5.67E-8_r8*BBR
          A(ij)=A(ij)-BBR
          FORCES(ng)%lrflx(i,j)=A(ij)*cff
        END DO
      END DO
# ifdef RUOYING_CASE1
!
!  Latent heat flux            (W/m^2)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "LH",       &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%lhflx(i,j)=A(ij)
        END DO
      END DO
!
!  Sensible heat flux            (W/m^2)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "HFX",      &
     &                           A, Asize)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%shflx(i,j)=A(ij)
        END DO
      END DO
# endif
#endif
#ifdef SHORTWAVE
!
!  Short wave radiation          (Celsius m/s)
!
      CALL AttrVect_exportRAttr (AttrVect_G(ng)%atm2ocn_AV, "SWDOWN",   &
     &                           A, Asize)
      cff=1.0_r8/(rho0*Cp)
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
          FORCES(ng)%srflx(i,j)=A(ij)*cff
        END DO
      END DO
#endif
!
!  Surface u-stress              (m2/s2)
!
! THIS NEEDS TO BE AVERAGED TO U POINTS

!      CALL AttrVect_exportRAttr (atm2ocn_AV, "USTRESS", A, Asize)
!      cff=1.0_r8/rho0
!      ij=0
!      DO j=JstrT,JendT
!        DO i=IstrU,IendR
!          ij=ij+1
!          sustr(i,j)=A(ij)*cff
!        END DO
!      END DO
!
!  Surface v-stress              (m2/s2)
!
! THIS NEEDS TO BE AVERAGED TO V POINTS

!      CALL AttrVect_exportRAttr (atm2ocn_AV, "VSTRESS", A, Asize)
!      cff=1.0_r8/rho0
!      ij=0
!      DO j=JstrV,JendR
!        DO i=IstrT,IendT
!          ij=ij+1
!          svstr(i,j)=A(ij)*cff
!        END DO
!      END DO
!
!  Apply boundary conditions.
!
#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Pair)
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Hair)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Tair)
# endif
# if defined BULK_FLUXES || defined ECOSIM
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Uwind)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%Vwind)
# endif
# ifdef CLOUDS
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%cloud)
# endif
# ifdef BULK_FLUXES
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%rain)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%lrflx)
#  ifdef RUOYING_CASE1
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%lhflx)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%shflx)
#  endif
# endif
# ifdef SHORTWAVE
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        FORCES(ng)%srflx)
# endif
!      CALL exchange_u2d_tile (ng, tile,                                &
!     &                        LBi, UBi, LBj, UBj,                      &
!     &                        sustr)
!      CALL exchange_v2d_tile (ng, tile,                                &
!     &                        LBi, UBi, LBj, UBj,                      &
!     &                        svstr)
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile boundaries.
!-----------------------------------------------------------------------
!
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%Pair)
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
    (defined SHORTWAVE && defined ANA_SRFLUX)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%Hair, FORCES(ng)%Tair)
# endif
# if defined BULK_FLUXES || defined ECOSIM
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%Uwind, FORCES(ng)%Vwind)
# endif
# ifdef CLOUDS
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%cloud)
# endif
# ifdef BULK_FLUXES
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%rain, FORCES(ng)%lrflx)
#  ifdef RUOYING_CASE1
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%lhflx, FORCES(ng)%shflx)
#  endif
# endif
# ifdef SHORTWAVE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    FORCES(ng)%srflx)
# endif
!      CALL mp_exchange2d (ng, tile, iNLM, 2,                           &
!     &                    LBi, UBi, LBj, UBj,                          &
!     &                    NghostPoints, EWperiodic, NSperiodic,        &
!     &                    sustr, svstr)
#endif
!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to atmosphere (WRF) model.
!-----------------------------------------------------------------------
!
!
!  Sea surface temperature       (degC)
!
      ij=0
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ij=ij+1
#ifdef SST_CONST
          A(ij)=29.0_r8  ! exp A. rhe 03/13/08
#else       
          A(ij)=OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)
#endif
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_G(ng)%ocn2atm_AV, "SST", A, &
     &                             Asize)
!
!  Send ocean fields to atmosphere model.
!
#ifdef MCT_INTERP_OC2AT
      CALL MCT_MatVecMul(AttrVect_G(ng)%ocn2atm_AV,O2AMatPlus,        &
     &                     ocn2atm_AV2)
      CALL MCT_Send (ocn2atm_AV2, Router_G(ng)%ROMStoWRF, MyError)
#else
      CALL MCT_Send (AttrVect_G(ng)%ocn2atm_AV, Router_G(ng)%ROMStoWRF, &
     &               MyError)
#endif
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) 'atmosphere model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,'a') 'ROMS sent data to WRF '
        END IF
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (A)
!
 10   FORMAT (' OCN2ATM_COUPLING - error while receiving fields from ', &
     &        a, i4)
 20   FORMAT (' OCN2ATM_COUPLING - error while sending fields to ',     &
     &        a, i4)
      RETURN
      END SUBROUTINE ocn2atm_coupling_tile

      SUBROUTINE finalize_ocn2atm_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and atmosphere models coupling data     !
!  streams.                                                             !
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
        CALL Router_clean (Router_G(ng)%ROMStoWRF, MyError)
        CALL AttrVect_clean (AttrVect_G(ng)%ocn2atm_AV, MyError)
        CALL GlobalSegMap_clean (GlobalSegMap_G(ng)%GSMapROMS, MyError)
      END DO
      RETURN

      END SUBROUTINE finalize_ocn2atm_coupling
