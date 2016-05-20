#include "wrfcpp.h"

      MODULE mct_wrf_coupler_params

#if defined COAWST_COUPLING && defined MCT_LIB


      USE mod_wrf_coupler_kinds

      implicit none
!
!  Number of coupling models.
!
      integer :: N_mctmodels

# if defined MCT_INTERP_OC2WV || defined MCT_INTERP_OC2AT || \
     defined MCT_INTERP_WV2AT
!
!  Sparse matrix weights
!
      integer, dimension(:), pointer :: sparse_rows
      integer, dimension(:), pointer :: sparse_cols
      integer, dimension(:), pointer :: dst_grid_imask
      integer, dimension(2) :: src_grid_dims, dst_grid_dims
      real(m8), dimension(:), pointer :: sparse_weights

      TYPE T_DST_GRID
        integer, pointer :: dst_mask(:)
      END TYPE T_DST_GRID
      TYPE (T_DST_GRID), allocatable :: O2A_CPLMASK(:,:)
      TYPE (T_DST_GRID), allocatable :: A2O_CPLMASK(:,:)
      TYPE (T_DST_GRID), allocatable :: W2A_CPLMASK(:,:)

# endif
!
!  Number of parallel nodes assigned to each model in the coupled
!  system.
!
      integer :: NnodesATM
      integer :: NnodesWAV
      integer :: NnodesOCN

# ifdef WRF_COUPLING
!
!  Parallel nodes assined to the atmosphere model.
!
        integer :: peATM_frst          ! first atmosphere parallel node
        integer :: peATM_last          ! last  atmosphere parallel node
# endif
# ifdef SWAN_COUPLING
!
!  Parallel nodes assined to the wave model.
!
        integer :: peWAV_frst          ! first atmosphere parallel node
        integer :: peWAV_last          ! last  atmosphere parallel node
# endif
# ifdef ROMS_COUPLING
!
!  Parallel nodes assined to the ocean model.
!
        integer :: peOCN_frst          ! first ocean parallel node
        integer :: peOCN_last          ! last  ocean parallel node
        integer, dimension(:), pointer :: roms_fwcoup
        integer, dimension(:), pointer :: roms_2wcoup
        integer, dimension(:), pointer :: roms_facoup
        integer, dimension(:), pointer :: roms_2acoup
# endif
!
!  Time interval (seconds) between coupling of models.
!
      real(m8) :: TI_ATM2WAV           ! atmosphere to wave coupling interval
      real(m8) :: TI_ATM2OCN           ! atmosphere to ocean coupling interval
      real(m8) :: TI_WAV2ATM           ! wave to atmosphere coupling interval
      real(m8) :: TI_WAV2OCN           ! wave to ocean coupling interval
      real(m8) :: TI_OCN2WAV           ! ocean to wave coupling interval
      real(m8) :: TI_OCN2ATM           ! ocean to atmosphere coupling interval
!
!  Number of atmosphere model time-steps and atmosphere model ID.
!
      integer :: Natm_grids
      integer :: Nocn_grids
      integer :: Nwav_grids

      real(m8), dimension(:), pointer :: dtocn
      real(m8), dimension(:), pointer :: dtwav
      real(m8), dimension(:), pointer :: dtatm
!
# ifdef WAVES_OCEAN
      integer, dimension(:,:), pointer :: nOCN2WAV
      integer, dimension(:,:), pointer :: nWAV2OCN
      integer, dimension(:,:), pointer :: nOCNFWAV
      integer, dimension(:,:), pointer :: nWAVFOCN
# endif
# ifdef AIR_OCEAN
      integer, dimension(:,:), pointer :: nOCN2ATM
      integer, dimension(:,:), pointer :: nATM2OCN
      integer, dimension(:,:), pointer :: nOCNFATM
      integer, dimension(:,:), pointer :: nATMFOCN
# endif
# ifdef AIR_WAVES
      integer, dimension(:,:), pointer :: nATM2WAV
      integer, dimension(:,:), pointer :: nWAV2ATM
      integer, dimension(:,:), pointer :: nATMFWAV
      integer, dimension(:,:), pointer :: nWAVFATM
# endif
!
!  Coupled model components IDs.
!
      integer, dimension(:), pointer :: ocnids
      integer, dimension(:), pointer :: wavids
      integer, dimension(:), pointer :: atmids
      integer :: OCNid
      integer :: WAVid
      integer :: ATMid

# ifdef MOVE_NESTS
!  If WRF has a moving nest then we need to keep some information.
      integer, dimension(:), pointer :: par_grid_ratio
      integer, dimension(:), pointer :: full_Isize
      integer, dimension(:), pointer :: full_Jsize
      integer, dimension(:), pointer :: parent_id
      integer, dimension(:), pointer :: moving_nest
# endif

      CONTAINS

      SUBROUTINE allocate_coupler_params
!=======================================================================
!                                                                      !
!  This routine initialize all the coupler vars.                       !
!                                                                      !
!=======================================================================

      integer :: i

# ifdef WAVES_OCEAN
      allocate (nOCN2WAV(Nocn_grids,Nwav_grids))
      allocate (nWAV2OCN(Nwav_grids,Nocn_grids))
      allocate (nOCNFWAV(Nocn_grids,Nwav_grids))
      allocate (nWAVFOCN(Nwav_grids,Nocn_grids))
# endif
# ifdef AIR_OCEAN
      allocate (nOCN2ATM(Nocn_grids,Natm_grids))
      allocate (nATM2OCN(Natm_grids,Nocn_grids))
      allocate (nOCNFATM(Nocn_grids,Natm_grids))
      allocate (nATMFOCN(Natm_grids,Nocn_grids))
# endif
# ifdef AIR_WAVES
      allocate (nATM2WAV(Natm_grids,Nwav_grids))
      allocate (nWAV2ATM(Nwav_grids,Natm_grids))
      allocate (nATMFWAV(Natm_grids,Nwav_grids))
      allocate (nWAVFATM(Nwav_grids,Natm_grids))
# endif

# if defined MCT_INTERP_OC2WV || defined MCT_INTERP_OC2AT || \
     defined MCT_INTERP_WV2AT
#  ifdef AIR_OCEAN
      allocate(O2A_CPLMASK(Nocn_grids,Natm_grids))
      allocate(A2O_CPLMASK(Natm_grids,Nocn_grids))
#  endif
#  ifdef AIR_WAVES
      allocate(W2A_CPLMASK(Nwav_grids,Natm_grids))
#  endif
# endif
# ifdef MOVE_NESTS
! If WRF has a moving nest, for now, it has to be the last grid.
      allocate (par_grid_ratio(Natm_grids))
      allocate (full_Isize(Natm_grids))
      allocate (full_Jsize(Natm_grids))
      allocate (parent_id(Natm_grids))
      allocate (moving_nest(Natm_grids))
      DO i=1,Natm_grids-1
       moving_nest(i)=0
      END DO
      moving_nest(Natm_grids)=1

# endif

      RETURN
      END SUBROUTINE allocate_coupler_params

#endif
      END MODULE mct_wrf_coupler_params
