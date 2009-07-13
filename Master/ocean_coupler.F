#include "cppdefs.h"
      MODULE ocean_coupler_mod

#if defined MODEL_COUPLING && defined MCT_LIB
!
!svn $Id: ocean_coupler.F 755 2008-09-14 19:07:08Z jcwarner $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group      Hernan G. Arango   !
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
      USE m_Transfer, ONLY: MCT_Send => send
      USE m_Transfer, ONLY: MCT_Recv => recv

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

# ifdef SWAN_COUPLING
      PUBLIC :: initialize_ocn2wav_coupling
      PUBLIC :: ocn2wav_coupling
      PUBLIC :: finalize_ocn2wav_coupling
# endif
# ifdef WRF_COUPLING
      PUBLIC :: initialize_ocn2atm_coupling
      PUBLIC :: ocn2atm_coupling
      PUBLIC :: finalize_ocn2atm_coupling
# endif
!
!  Declarations.
!
      TYPE T_GlobalSegMap_G
        TYPE(GlobalSegMap) :: GSMapROMS       ! GloabalSegMap variables
# ifdef MCT_INTERP_OC2WV
        TYPE(GlobalSegMap) :: GSMapSWAN       ! GloabalSegMap variables
#endif
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:)

# ifdef SWAN_COUPLING
      TYPE(AttrVect) :: wav2ocn_AV            ! AttrVect variables
      TYPE(AttrVect) :: ocn2wav_AV 
      TYPE(Router)   :: ROMStoSWAN            ! Router variables
#  ifdef REFINED_GRID
      TYPE(AttrVect) :: wav2ocn_AV1           ! AttrVect variables
      TYPE(AttrVect) :: ocn2wav_AV1 
      TYPE(Router)   :: ROMStoSWAN1           ! Router variables
      TYPE(AttrVect) :: wav2ocn_AV2           ! AttrVect variables
      TYPE(AttrVect) :: ocn2wav_AV2 
      TYPE(Router)   :: ROMStoSWAN2           ! Router variables
#  endif
# endif
# ifdef WRF_COUPLING
      TYPE(AttrVect) :: atm2ocn_AV            ! AttrVect variables
      TYPE(AttrVect) :: ocn2atm_AV 
      TYPE(Router)   :: ROMStoWRF             ! Router variables
# endif
# ifdef MCT_INTERP_OC2WV
      TYPE(AttrVect) :: wav2ocn_AV2           ! AttrVect variables
      TYPE(AttrVect) :: ocn2wav_AV2 

      TYPE(SparseMatrix) :: sMatO             ! Sparse matrix elements
      TYPE(SparseMatrix) :: sMatW             ! Sparse matrix elements
      TYPE(SparseMatrixPlus) :: W2OMatPlus    ! Sparse matrix plus elements
      TYPE(SparseMatrixPlus) :: O2WMatPlus    ! Sparse matrix plus elements
# endif

# ifdef MCT_INTERP_OC2AT
      TYPE(GlobalSegMap) :: GSMapROMS2        ! GloabalSegMap variables
      TYPE(GlobalSegMap) :: GSMapWRF          ! GloabalSegMap variables
      TYPE(AttrVect) :: atm2ocn_AV2           ! AttrVect variables
      TYPE(AttrVect) :: ocn2atm_AV2 

#  if !defined MCT_INTERP_OC2WV
      TYPE(SparseMatrix) :: sMatO             ! Sparse matrix elements
#  endif
      TYPE(SparseMatrix) :: sMatA             ! Sparse matrix elements
      TYPE(SparseMatrixPlus) :: A2OMatPlus    ! Sparse matrix plus elements
      TYPE(SparseMatrixPlus) :: O2AMatPlus    ! Sparse matrix plus elements
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

#endif
      END MODULE ocean_coupler_mod
