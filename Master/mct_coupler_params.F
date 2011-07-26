#include "cppdefs.h"

      MODULE mct_coupler_params

#if defined MODEL_COUPLING && defined MCT_LIB
      USE mod_kinds
      implicit none
!
!  Number of coupling models.
!
      integer :: N_mctmodels

# if defined MCT_INTERP_WV2AT
!
!  Sparse matrix weights
!
      integer, dimension(:), pointer :: sparse_rows2
      integer, dimension(:), pointer :: sparse_cols2
      real(r8), dimension(:), pointer :: sparse_weights2
# endif
# if defined MCT_INTERP_OC2WV || defined MCT_INTERP_OC2AT
!
!  Sparse matrix weights
!
      integer, dimension(:), pointer :: sparse_rows
      integer, dimension(:), pointer :: sparse_cols
      integer, dimension(:), pointer :: dst_grid_imask
      real(r8), dimension(:), pointer :: sparse_weights
# endif
!
!  Number of parallel nodes assigned to each model in the coupled
!  system.
!
      integer :: NnodesATM
      integer :: NnodesWAV
      integer :: NnodesOCN
!
!  Time interval (seconds) between coupling of models.
!
      real(r8) :: TI_ATM_WAV        ! atmosphere-wave coupling interval
      real(r8) :: TI_ATM_OCN        ! atmosphere-ocean coupling interval
      real(r8) :: TI_WAV_OCN        ! wave-ocean coupling interval
!
!  Number of atmosphere model time-steps and atmosphere model ID.
!
      integer, dimension(NestedGrids) :: nOCN_ATM
      integer, dimension(NestedGrids) :: nOCN_WAV
!
!  Coupled model components IDs.
!
# ifdef REFINED_GRID
      integer, dimension(:), pointer :: ocnids
      integer, dimension(:), pointer :: wavids
# endif
      integer :: OCNid
      integer :: WAVid
      integer :: ATMid
!
!  WRF grid number to couple to roms and or swan grids. 
!  This logic is currently being updated to allow multiple 
!  grids from multiple models to interact.
!
      integer :: WRF_CPL_GRID = 1

!     CONTAINS
#endif
      END MODULE mct_coupler_params
