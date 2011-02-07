!
!svn $Id: mct_roms_refdif.h 1009 2007-08-19 16:51:29Z jcwarner $
!==================================================== John C. Warner ===
!=================================================== Kevin Haas     ====
!  Copyright (c) 2002-2007 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module is for coupling ROMS/TOMS to RefDif wave model using    !
!  the Model Coupling Toolkit (MCT), developed at the Argonne National !
!  Laboratory.                                                         !
!                                                                      !
!=======================================================================
!
!
!  Component Model Registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain Decomposition Descriptor DataType and associated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field Storage DataType and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_lsize => lsize
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communications scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY: MCT_Send => send
      USE m_Transfer, ONLY: MCT_Recv => recv
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
!
      implicit none
!
      PRIVATE

      PUBLIC :: initialize_waves_coupling
      PUBLIC :: waves_coupling
      PUBLIC :: finalize_waves_coupling
!
!  Declarations.
!
      TYPE(GlobalSegMap) :: GSMapROMS         ! GloabalSegMap variables

      TYPE(AttrVect) :: FromWavesAV           ! AttrVect variables
      TYPE(AttrVect) :: ToWavesAV 

      TYPE(Router) :: RoutROMS                ! Router variables

      CONTAINS

      SUBROUTINE initialize_waves_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  Initialize waves and ocean modeles coupling stream. This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_kinds
      USE mod_scalars
!
!  Imported variable definitions.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.  
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, Jsize, MyError
      integer :: j, jc, nprocs

      integer, allocatable :: length(:)
      integer, allocatable :: start(:)

# include "tile.h"
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
        IstrR=Istr-1
      ELSE
        IstrR=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        IendR=Iend+1
      ELSE
        IendR=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        JstrR=Jstr-1
      ELSE
        JstrR=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        JendR=Jend+1
      ELSE
        JendR=Jend
      END IF
!
!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.

      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, nprocs, MyError)
!
!  Initialize MCT coupled model registry.
!
      CALL MCTWorld_init (2, MPI_COMM_WORLD, OCN_COMM_WORLD, OcnId)
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
      CALL GlobalSegMap_init (GSMapROMS, start, length, 0,              &
     &                        OCN_COMM_WORLD, OcnId)
!
!  Initialize Attribute vectors.  Their size is the number of grid point
!  on this processor.
!
      Asize=GlobalSegMap_lsize(GSMapROMS, OCN_COMM_WORLD)
      CALL AttrVect_init(FromWavesAV,                                   &
     &     rList="DWAVE:HWAVE:LWAVE:WAVE_BREAK:WAVE_DISSIP",            &
     &     lsize=Asize)
!
!  Initialize oceanAv with one real attribute for now.
!
      CALL AttrVect_init (ToWavesAV,                                    &
     &                    rList="XR:YR:DEPTH:UBAR:VBAR:ZETA",           &
     &                    lsize=Asize)
      CALL AttrVect_zero (ToWavesAV)
!
!  Initialize a router to the Waves component.
!
      CALL Router_init (WavId, GSMapROMS, OCN_COMM_WORLD, RoutROMS)
!
!  Deallocate working arrays.
!
      IF (allocated(start)) THEN
        deallocate (start)
      END IF
      IF (allocated(length)) THEN
        deallocate (length)
      END IF

      RETURN
      END SUBROUTINE initialize_waves_coupling

      SUBROUTINE waves_coupling (ng, tile)
!
!=======================================================================
!                                                                      !
!  This subroutine acquires the coupling data streams between waves    !
!  and ocean models. Currently, the following data streams are         !
!  processed:                                                          !
!                                                                      !
!  Fields acquired from the WAVE Model:                                !
!                                                                      !
!     * Dwave      Wave direction.                                     !
!     * Hwave      Wave height.                                        !
!     * Lwave      Wave length.                                        !
!     * Wave_break Percent of breakig waves.                           !
!     * Wave_dissip Wave energy dissipation.                           !
!                                                                      !
!  Fields sent to the WAVE Model:                                      !
!                                                                      !
!     * xr         x-rho coordinates                                   !
!     * yr         y-rho coordinates                                   !
!     * h          Bottom elevation.                                   !
!     * ubar       Depth integrated xi-direction velocity.             !
!     * vbar       Depth integrated eta-direction velocity.            !
!     * zeta       Water surface elevation.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 36)
#endif
      CALL waves_coupling_tile (ng, Istr, Iend, Jstr, Jend,             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          knew(ng),                               &
     &                          GRID(ng) % h,                           &
     &                          GRID(ng) % xr,                          &
     &                          GRID(ng) % yr,                          &
     &                          OCEAN(ng) % ubar,                       &
     &                          OCEAN(ng) % vbar,                       &
     &                          OCEAN(ng) % zeta,                       &
     &                          FORCES(ng) % Dwave,                     &
     &                          FORCES(ng) % Hwave,                     &
     &                          FORCES(ng) % Lwave,                     &
#ifdef ROLLER_SVENDSEN
     &                          FORCES(ng) % Wave_break,                &
#endif
     &                          FORCES(ng) % Wave_dissip)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 36)
#endif

      RETURN
      END SUBROUTINE waves_coupling
!
!***********************************************************************
      SUBROUTINE waves_coupling_tile (ng, Istr, Iend, Jstr, Jend,       &
     &                                LBi, UBi, LBj, UBj,               &
     &                                knew, h, xr, yr,                  &
     &                                ubar, vbar, zeta,                 &
     &                                Dwave, Hwave, Lwave,              &
#ifdef ROLLER_SVENDSEN
     &                                Wave_break,                       &
#endif
     &                                Wave_dissip)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_iounits
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
      USE bc_2d_mod
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: knew

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)

      real(r8), intent(inout) :: Dwave(LBi:,LBj:)
      real(r8), intent(inout) :: Hwave(LBi:,LBj:)
      real(r8), intent(inout) :: Lwave(LBi:,LBj:)
# ifdef ROLLER_SVENDSEN
      real(r8), intent(inout) :: Wave_break(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: Wave_dissip(LBi:,LBj:)

#else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(inout) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Lwave(LBi:UBi,LBj:UBj)
# ifdef ROLLER_SVENDSEN
      real(r8), intent(inout) :: Wave_break(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: Wave_dissip(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
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

      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: Asize, MyError
      integer :: i, ij, j

      integer, pointer :: points(:)

      real(r8) :: cff, ramp, sqrt2
      real(r8), parameter ::  Lwave_max=500.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: ubar_rho
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: vbar_rho
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: A2

      real(r8), pointer :: A(:)
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
        IstrR=Istr-1
      ELSE
        IstrR=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        IendR=Iend+1
      ELSE
        IendR=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        JstrR=Jstr-1
      ELSE
        JstrR=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        JendR=Jend+1
      ELSE
        JendR=Jend
      END IF
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize (GSMapROMS, OCN_COMM_WORLD)
      allocate ( A(Asize) )
!
!-----------------------------------------------------------------------
!  Receive needed fields from wave model.
!-----------------------------------------------------------------------
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL MCT_Recv (FromWavesAV, RoutROMS, MyError)
      IF (Master) THEN
        WRITE (stdout,10) tdays(ng)
      END IF
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,30) MyError
        END IF
        CALL finalize_waves_coupling
      END IF
!
!  Set ramp coefficient.
!
!!    ramp=MIN((tdays(ng)-dstart)*4.0_r8,1.0_r8)
      ramp=1.0_r8
      sqrt2=sqrt(2.0_r8)
      DO j=JstrR-1,JendR+1
        DO i=IstrR-1,IendR+1
          A2(i,j)=0.0_r8
        END DO
      END DO
!
!  Wave direction (radians).
!
      CALL AttrVect_exportRAttr (FromWavesAV, "DWAVE", A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          Dwave(i,j)=270.0_r8*3.14159/180.0_r8-A(ij)*deg2rad
        END DO
      END DO
!
!  Wave height.
!
      CALL AttrVect_exportRAttr (FromWavesAV, "HWAVE", A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A2(i,j)=MAX(0.0_r8,A(ij)*ramp)*sqrt2
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
         Hwave(i,j)=0.2_r8*(A2(i  ,j  )+A2(i-1,j  )+                    &
    &                       A2(i+1,j  )+A2(i  ,j-1)+                    &
    &                       A2(i  ,j+1))
        END DO
      END DO
      CALL bc_r2d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Hwave)
!
!  Wave length (m).
!
      CALL AttrVect_exportRAttr (FromWavesAV, "LWAVE", A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
            A2(i,j)=2.0_r8*pi/MAX(eps,A(ij))
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Lwave(i,j)=0.2_r8*(A2(i  ,j  )+A2(i-1,j  )+                    &
     &                       A2(i+1,j  )+A2(i  ,j-1)+                    &
     &                       A2(i  ,j+1))
        END DO
      END DO
      CALL bc_r2d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Lwave)
#ifdef ROLLER_SVENDSEN
!
!  Percent wave breaking.
!  
      CALL AttrVect_exportRAttr (FromWavesAV, "WAVE_BREAK", A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A2(i,j)=MAX(0.0_r8,A(ij))
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Wave_break(i,j)=0.2_r8*(A2(i  ,j  )+A2(i-1,j  )+               &
     &                            A2(i+1,j  )+A2(i  ,j-1)+               &
     &                            A2(i  ,j+1))
        END DO
      END DO
      CALL bc_r2d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Wave_break)
#endif
!
!  Wave dissipation.
!
      CALL AttrVect_exportRAttr (FromWavesAV, "WAVE_DISSIP", A, Asize)
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A2(i,j)=MAX(0.0_r8,A(ij)*ramp)
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Wave_dissip(i,j)=0.2_r8*(A2(i  ,j  )+A2(i-1,j  )+              &
     &                             A2(i+1,j  )+A2(i  ,j-1)+              &
     &                             A2(i  ,j+1))
        END DO
      END DO
      CALL bc_r2d_tile (ng, Istr, Iend, Jstr, Jend,                     &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Wave_dissip)
#if defined EW_PERIODIC || defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Dwave)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Hwave)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Lwave)
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Wave_dissip)
# ifdef ROLLER_SVENDSEN
      CALL exchange_r2d_tile (ng, Istr, Iend, Jstr, Jend,               &
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
      CALL mp_exchange2d (ng, iNLM, 3, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Dwave, Hwave, Lwave)
# ifdef ROLLER_SVENDSEN
      CALL mp_exchange2d (ng, iNLM, 1, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Wave_break)
# endif
      CALL mp_exchange2d (ng, iNLM, 1, Istr, Iend, Jstr, Jend,          &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Wave_dissip)
#endif
!
!-----------------------------------------------------------------------
!  Schedule and send required fields to wave model.
!-----------------------------------------------------------------------
!
!  xr rho point.
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=xr(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToWavesAV, "XR", A, Asize)
!
!  yr rho point.
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=yr(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToWavesAV, "YR", A, Asize)
!
!  Depth (bathymetry).
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=h(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToWavesAV, "DEPTH", A, Asize)
!
!  Water level (free-surface).
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=0.2_r8*(zeta(i  ,j  ,2)+zeta(i-1,j  ,2)+                &
     &                  zeta(i+1,j  ,2)+zeta(i  ,j-1,2)+                &
     &                  zeta(i  ,j+1,2))
        END DO
      END DO
      CALL AttrVect_importRAttr (ToWavesAV, "ZETA", A, Asize)
!
!  Vertically-integrated U-velocity at RHO-points.
!
      DO j=JstrR,JendR
        DO i=Istr,Iend
          ubar_rho(i,j)=0.5_r8*(ubar(i,j,knew)+ubar(i+1,j,knew))
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
        ubar_rho(Istr-1,Jstr-1)=0.5_r8*(ubar_rho(Istr  ,Jstr-1)+        &
     &                                  ubar_rho(Istr-1,Jstr  ))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar_rho(Iend+1,Jstr-1)=0.5_r8*(ubar_rho(Iend  ,Jstr-1)+        &
     &                                  ubar_rho(Iend+1,Jstr  ))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        ubar_rho(Istr-1,Jend+1)=0.5_r8*(ubar_rho(Istr-1,Jend  )+        &
     &                                  ubar_rho(Istr  ,Jend+1))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar_rho(Iend+1,Jend+1)=0.5_r8*(ubar_rho(Iend+1,Jend  )+        &
     &                                  ubar_rho(Iend  ,Jend+1))
      END IF
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=ubar_rho(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToWavesAV, "UBAR", A, Asize)
!
!  Vertically-integrated V-velocity at RHO-points.
!
      DO j=Jstr,Jend
        DO i=IstrR,IendR
          vbar_rho(i,j)=0.5_r8*(vbar(i,j,knew)+vbar(i,j+1,knew))
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
        vbar_rho(Istr-1,Jstr-1)=0.5_r8*(vbar_rho(Istr  ,Jstr-1)+        &
     &                                  vbar_rho(Istr-1,Jstr  ))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        vbar_rho(Iend+1,Jstr-1)=0.5_r8*(vbar_rho(Iend  ,Jstr-1)+        &
     &                                  vbar_rho(Iend+1,Jstr  ))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        vbar_rho(Istr-1,Jend+1)=0.5_r8*(vbar_rho(Istr-1,Jend  )+        &
     &                                  vbar_rho(Istr  ,Jend+1))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        vbar_rho(Iend+1,Jend+1)=0.5_r8*(vbar_rho(Iend+1,Jend  )+        &
     &                                  vbar_rho(Iend  ,Jend+1))
      END IF
!
      ij=0
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ij=ij+1
          A(ij)=vbar_rho(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToWavesAV, "VBAR", A, Asize)
!
!  Send the data.
!
      CALL MCT_Send (ToWavesAV, RoutROMS, MyError)
      IF (Master) THEN
        WRITE (stdout,20) tdays(ng)
      END IF
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,40) MyError
        END IF
        CALL finalize_waves_coupling
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (A)

 10   FORMAT (1x,'Wave-Ocean models coupling, receive',t64,'t= ',f12.4)
 20   FORMAT (1x,'Wave-Ocean models coupling, send,',t64,'t= ',f12.4)
 30   FORMAT (/,' WAVES_COUPLING - error while receiving data,',         &
     &          ' MyError = ',i4)   
 40   FORMAT (/,' WAVES_COUPLING - error while sending data,',           &
     &          ' MyError = ',i4)   

      RETURN
      END SUBROUTINE waves_coupling_tile

      SUBROUTINE finalize_waves_coupling
!
!========================================================================
!                                                                     ===
!  This routine terminates execution during coupling error.           ===
!                                                                     ===
!========================================================================
!
!  Local variable declarations.
!
      integer :: MyStatus
!
!-----------------------------------------------------------------------
!  Terminate MPI execution environment.
!-----------------------------------------------------------------------
!
      CALL Router_clean (RoutROMS)
      CALL AttrVect_clean (ToWavesAV)
      CALL AttrVect_clean (FromWavesAV)
      CALL GlobalSegMap_clean (GSMapROMS)
      CALL MCTWorld_clean ()
      CALL mpi_finalize (MyStatus)

      STOP
      END SUBROUTINE finalize_waves_coupling

