#include "cppdefs.h"
      MODULE ocean_coupler_mod

#if defined ROMS_COUPLING && defined MCT_LIB
!
!svn $Id: ocean_coupler.F 755 2008-09-14 19:07:08Z jcwarner $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between        !
!  ROMS/TOMS and other coupled model(s)  via the Model Coupling        !
!  Toolkit (MCT), developed at the Argonne National Laboratory.        !
!                                                                      !
!=======================================================================
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
      USE m_AttrVect, ONLY : AttrVect_copy => copy
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
      USE m_Transfer, ONLY: MCT_send => send
      USE m_Transfer, ONLY: MCT_recv => recv
      USE m_Transfer, ONLY: MCT_isend => isend
      USE m_Transfer, ONLY: MCT_irecv => irecv
      USE m_Transfer, ONLY: MCT_waitr => waitrecv
      USE m_Transfer, ONLY: MCT_waits => waitsend

# if defined MCT_INTERP_OC2WV || defined MCT_INTERP_OC2AT
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
      implicit none
!
      PRIVATE

      PUBLIC :: ocean_coupling
# ifdef SWAN_COUPLING
      PUBLIC :: initialize_ocn2wav_coupling
      PUBLIC :: initialize_ocn2wav_routers
      PUBLIC :: ocn2wav_coupling
      PUBLIC :: ocnfwav_coupling
      PUBLIC :: finalize_ocn2wav_coupling
# endif
# ifdef WRF_COUPLING
      PUBLIC :: initialize_ocn2atm_coupling
      PUBLIC :: initialize_ocn2atm_routers
      PUBLIC :: ocn2atm_coupling
      PUBLIC :: ocnfatm_coupling
      PUBLIC :: finalize_ocn2atm_coupling
# endif
!
!  Declarations.
!
      TYPE T_GlobalSegMap_G
        TYPE(GlobalSegMap) :: GSMapROMS       ! GloabalSegMap variables
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:)

# ifdef MCT_INTERP_OC2WV
      TYPE T_GSMapInterp_W
        TYPE(GlobalSegMap) :: GSMapSWAN       ! GloabalSegMap variables
      END TYPE T_GSMapInterp_W
      TYPE (T_GSMapInterp_W), ALLOCATABLE :: GSMapInterp_W(:,:)
# endif
# ifdef MCT_INTERP_OC2AT
      TYPE T_GSMapInterp_A
        TYPE(GlobalSegMap) :: GSMapWRF        ! GloabalSegMap variables
      END TYPE T_GSMapInterp_A
      TYPE (T_GSMapInterp_A), ALLOCATABLE :: GSMapInterp_A(:,:)
# endif

      TYPE T_AttrVect_G
# ifdef SWAN_COUPLING
        TYPE(AttrVect) :: wav2ocn_AV          ! AttrVect variables
        TYPE(AttrVect) :: ocn2wav_AV
# endif
# ifdef WRF_COUPLING
      TYPE(AttrVect) :: atm2ocn_AV            ! AttrVect variables
      TYPE(AttrVect) :: ocn2atm_AV 
# endif
      END TYPE T_AttrVect_G
      TYPE (T_AttrVect_G), ALLOCATABLE :: AttrVect_G(:)

# ifdef SWAN_COUPLING
      TYPE T_Router_W
        TYPE(Router)   :: ROMStoSWAN          ! Router variables
      END TYPE T_Router_W
      TYPE (T_Router_W), ALLOCATABLE :: Router_W(:,:)
# endif
# ifdef WRF_COUPLING
      TYPE T_Router_A
        TYPE(Router)   :: ROMStoWRF           ! Router variables
      END TYPE T_Router_A
      TYPE (T_Router_A), ALLOCATABLE :: Router_A(:,:)
# endif

# ifdef MCT_INTERP_OC2WV
      TYPE T_AV2_W
        TYPE(AttrVect) :: wav2ocn_AV2         ! AttrVect variables
        TYPE(AttrVect) :: ocn2wav_AV2 
      END TYPE T_AV2_W
      TYPE (T_AV2_W), ALLOCATABLE :: AV2_W(:,:)

      TYPE(SparseMatrix) :: sMatO             ! Sparse matrix elements
      TYPE(SparseMatrix) :: sMatW             ! Sparse matrix elements

      TYPE T_SMPlus_W
        TYPE(SparseMatrixPlus) :: W2OMatPlus  ! Sparse matrix plus elements
        TYPE(SparseMatrixPlus) :: O2WMatPlus  ! Sparse matrix plus elements
      END TYPE T_SMPlus_W
      TYPE (T_SMPlus_W), ALLOCATABLE :: SMPlus_W(:,:)
# endif

# ifdef MCT_INTERP_OC2AT
      TYPE T_AV2_A
        TYPE(AttrVect) :: atm2ocn_AV2           ! AttrVect variables
        TYPE(AttrVect) :: ocn2atm_AV2
      END TYPE T_AV2_A
      TYPE (T_AV2_A), ALLOCATABLE :: AV2_A(:,:)

#  if !defined MCT_INTERP_OC2WV
      TYPE(SparseMatrix) :: sMatO             ! Sparse matrix elements
#  endif
      TYPE(SparseMatrix) :: sMatA             ! Sparse matrix elements

      TYPE T_SMPlus_A
        TYPE(SparseMatrixPlus) :: A2OMatPlus  ! Sparse matrix plus elements
        TYPE(SparseMatrixPlus) :: O2AMatPlus  ! Sparse matrix plus elements
      END TYPE T_SMPlus_A
      TYPE (T_SMPlus_A), ALLOCATABLE :: SMPlus_A(:,:)
# endif

      CONTAINS
/*
************************************************************************
*  Include model specific communication routines.
************************************************************************
*/

# ifdef SWAN_COUPLING
#  include "mct_roms_swan.h"
# endif
# ifdef REFDIF_COUPLING
#  include "mct_roms_refdif.h"
# endif
# ifdef WRF_COUPLING
#  include "mct_roms_wrf.h"
# endif


/*
************************************************************************
*  Include main driver to determin which grids are performing excahgnes.
************************************************************************
*/

      SUBROUTINE ocean_coupling (nl)
!
!=======================================================================
!                                                                      !
!  Determine which roms grids are going to exchange data to otehr      !
!  model grids and call those exchange.                                !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mct_coupler_params
      USE mod_scalars
!
!  Imported variable definitions.
!
      integer, intent(in) :: nl
!
!  Local variable declarations.
!
      integer :: MyError, nprocs, tile
      integer :: ng, iw, ia, ig, nlay, offset
!
!-----------------------------------------------------------------------
!
# ifdef AIR_OCEAN
!
!-----------------------------------------------------------------------
!  Couple ocean to atmosphere model every nOCN_ATM timesteps.
!-----------------------------------------------------------------------
!
      IF (nl.eq.1) THEN
        DO nlay=1,NestLayers
          DO ig=1,GridsInLayer(nlay)
            ng=GridNumber(ig,nlay)
            DO ia=1,Natm_grids
              offset=-1 !nlay-NestLayers
              IF (MOD(iic(1)+offset,nOCNFATM(1,1)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ocnfatm_coupling (ng, ia, tile)
                END DO
              END IF
            END DO
          END DO
        END DO
      END IF
!
      IF (nl.eq.1) THEN
        DO nlay=1,NestLayers
          DO ig=1,GridsInLayer(nlay)
            ng=GridNumber(ig,nlay)
            DO ia=1,Natm_grids
              offset=-1 !nlay-NestLayers
              IF (MOD(iic(1)+offset,nOCN2ATM(1,1)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ocn2atm_coupling (ng, ia, tile)
                END DO
              END IF
            END DO
          END DO
        END DO
      END IF
# endif
# ifdef WAVES_OCEAN
!
!-----------------------------------------------------------------------
!  Exchange fields from ocn to wav every TI_OCN2WAV steps and
!                  from wav to ocn every TI_WAV2OCN steps.
!-----------------------------------------------------------------------
!
      IF (nl.eq.1) THEN
        DO iw=1,Nwav_grids
          DO nlay=1,NestLayers
            DO ig=1,GridsInLayer(nlay)
              ng=GridNumber(ig,nlay)
              offset=-1 !nlay-NestLayers
              IF (MOD(iic(1)+offset,nOCNFWAV(1,1)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ocnfwav_coupling (ng, iw, tile)
                END DO
              END IF
            END DO
          END DO
        END DO
      END IF
!
      IF (nl.eq.1) THEN
        DO iw=1,Nwav_grids
          DO nlay=1,NestLayers
            DO ig=1,GridsInLayer(nlay)
              ng=GridNumber(ig,nlay)
              offset=-1 !nlay-NestLayers
              IF (MOD(iic(1)+offset,nOCN2WAV(1,1)).eq.0) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ocn2wav_coupling (ng, iw, tile)
                END DO
              END IF
            END DO
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE ocean_coupling

#endif
      END MODULE ocean_coupler_mod
