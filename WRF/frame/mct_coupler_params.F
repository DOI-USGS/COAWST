#include "wrfcpp.h"

      MODULE mct_coupler_params

#if defined MODEL_COUPLING && defined MCT_LIB

      implicit none
!
!  Number of coupling models.
!
      integer :: N_mctmodels
!
!  Time interval (seconds) between coupling of models.
!
      real :: TI_ATM_WAV        ! atmosphere-wave coupling interval
      real :: TI_ATM_OCN        ! atmosphere-ocean coupling interval
!
!  Coupled model components IDs.
!
# ifdef REFINED_GRID
      integer, dimension(NestedGrids) :: ocnids
      integer, dimension(NestedGrids) :: wavids
# endif
      integer :: OCNid
      integer :: WAVid
      integer :: ATMid
!
!  WRF grid number to couple to roms and or swan grids. 
!  This logic is currently being updated to allow multiple 
!  grids from multiple models to interact.
!
      integer :: WRF_CPL_GRID

#endif

      END MODULE mct_coupler_params
