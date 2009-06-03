#include "cppdefs.h"
      MODULE mod_iparam
!
!svn $Id: mod_iparam.F 001 2008-11-17 04:22:18Z jcwarner $
!================================================== Hernan G. Arango ===
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Grid parameters:                                                    !
!                                                                      !
!  Im         Number of global grid points in the XI-direction         !
!               for each nested grid.                                  !
!  Jm         Number of global grid points in the ETA-direction        !
!               for each nested grid.                                  !
!  Lm         Number of interior grid points in the XI-direction       !
!               for each nested grid.                                  !
!  Mm         Number of internal grid points in the ETA-direction.     !
!               for each nested grid.                                  !
!  N          Number of vertical levels for each nested grid.          !
!  Ngrids     Number of nested and/or connected grids to solve.        !
!  NtileI     Number of XI-direction tiles or domain partitions for    !
!               each nested grid. Values used to compute tile ranges.  !
!  NtileJ     Number of ETA-direction tiles or domain partitions for   !
!               each nested grid. Values used to compute tile ranges.  !
!  NtileX     Number of XI-direction tiles or domain partitions for    !
!               each nested grid. Values used in parallel loops.       !
!  NtileE     Number of ETA-direction tiles or domain partitions for   !
!               each nested grid. Values used in parallel loops.       !
!  HaloSizeI  Maximum halo size, in grid points, in XI-direction.      ! 
!  HaloSizeJ  Maximum halo size, in grid points, in ETA-direction.     ! 
!  TileSide   Maximun tile side length in XI- or ETA-directions.       ! 
!  TileSize   Maximum tile size.                                       ! 
!                                                                      !
!=======================================================================
!
        USE mod_kinds
!
        implicit none
!
!-----------------------------------------------------------------------
!  Lower and upper bounds indices per domain partition for all grids.
!-----------------------------------------------------------------------
!
        TYPE T_BOUNDS
          integer, pointer :: tile(:)  ! tile partition

          integer, pointer :: LBi(:)   ! lower bound I-dimension
          integer, pointer :: UBi(:)   ! upper bound I-dimension
          integer, pointer :: LBj(:)   ! lower bound J-dimension
          integer, pointer :: UBj(:)   ! upper bound J-dimension

          integer, pointer :: Istr(:)  ! starting tile I-direction
          integer, pointer :: Iend(:)  ! ending   tile I-direction
          integer, pointer :: Jstr(:)  ! starting tile J-direction
          integer, pointer :: Jend(:)  ! ending   tile J-direction

          integer, pointer :: IstrR(:) ! starting tile I-direction (RHO)
          integer, pointer :: IstrT(:) ! starting nest I-direction (RHO)
          integer, pointer :: IstrU(:) ! starting tile I-direction (U)
          integer, pointer :: IendR(:) ! ending   tile I-direction (RHO)
          integer, pointer :: IendT(:) ! ending   nest I-direction (RHO)

          integer, pointer :: JstrR(:) ! starting tile J-direction (RHO)
          integer, pointer :: JstrT(:) ! starting nest J-direction (RHO)
          integer, pointer :: JstrV(:) ! starting tile J-direction (V)
          integer, pointer :: JendR(:) ! ending   tile J-direction (RHO)
          integer, pointer :: JendT(:) ! ending   nest J-direction (RHO)

          integer, pointer :: IstrTU(:) ! 
          integer, pointer :: JstrTV(:) ! 
        END TYPE T_BOUNDS

        TYPE (T_BOUNDS), allocatable :: BOUNDS(:)
!
!-----------------------------------------------------------------------
!  Model grid(s) parameters.
!-----------------------------------------------------------------------
!
!  Number of nested and/or connected grids to solve.
!
        integer, parameter :: Ngrids = NestedGrids
!
!  Number of interior RHO-points in the XI- and ETA-directions. The
!  size of models state variables (C-grid) at input and output are:
!
!    RH0-type variables:  [0:Lm+1, 0:Mm+1]        ----v(i,j+1)----
!    PSI-type variables:  [1:Lm+1, 1:Mm+1]        |              |
!      U-type variables:  [1:Lm+1, 0:Mm+1]     u(i,j)  r(i,j)  u(i+1,j)
!      V-type variables:  [0:Lm+1, 1:Mm+1]        |              |
!                                                 -----v(i,j)-----
        integer, dimension(Ngrids) :: Lm
        integer, dimension(Ngrids) :: Mm
!
!  Global horizontal size of model arrays including padding.  All the
!  model state arrays are of same size to facilitate parallelization.
!
        integer, dimension(Ngrids) :: Im
        integer, dimension(Ngrids) :: Jm
!
!  Number of vertical levels. The vertical ranges of model state
!  variables are:
!                                                 -----W(i,j,k)-----
!    RHO-, U-, V-type variables: [1:N]            |                |
!              W-type variables: [0:N]            |    r(i,j,k)    |  
!                                                 |                |
!                                                 ----W(i,j,k-1)----
        integer, dimension(Ngrids) :: N
!
# ifdef REFINED_GRID
        integer, dimension(Ngrids) :: nrefined
        integer, dimension(Ngrids) :: myparent
        integer, dimension(Ngrids) :: mychild
# endif
!
!-----------------------------------------------------------------------
!  Domain partition parameters.
!-----------------------------------------------------------------------
!
!  Number of tiles or domain partitions in the XI- and ETA-directions.
!  These values are used to compute tile ranges [Istr:Iend, Jstr:Jend].
!
        integer, dimension(Ngrids) :: NtileI
        integer, dimension(Ngrids) :: NtileJ
!
!  Number of tiles or domain partitions in the XI- and ETA-directions.
!  These values are used to parallel loops to differentiate between
!  shared-memory and distributed-memory.  Notice that in distributed
!  memory both values are set to one.
!
        integer, dimension(Ngrids) :: NtileX
        integer, dimension(Ngrids) :: NtileE
!
!  Maximum number of points in the halo region in the XI- and
!  ETA-directions.
!
        integer, dimension(Ngrids) :: HaloSizeI
        integer, dimension(Ngrids) :: HaloSizeJ
!
!  Maximum tile side length in XI- or ETA-directions.
!
        integer, dimension(Ngrids) :: TileSide
!
!  Maximum number of points in a tile partition.
!
        integer, dimension(Ngrids) :: TileSize
!
!  Set number of ghost-points in the halo region.  It is only used
!  in distributed-memory applications.
!
        integer :: NghostPoints = GHOST_POINTS

      CONTAINS

        SUBROUTINE initialize_param
!
!=======================================================================
!                                                                      !
!  This routine initializes several parameters in module "mod_param"   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
        integer :: ng
!
!-----------------------------------------------------------------------
!  Derived dimension parameters.
!-----------------------------------------------------------------------
!
#ifdef REFINED_GRID
!
!  Set refined grid parameters.
!
      nrefined(1)=1
      myparent(1)=0
      mychild(1)=2
# if defined DOGBONE
      nrefined(2)=3
      myparent(2)=1
# elif defined ISABEL
      nrefined(2)=5
      nrefined(3)=5
!      nrefined(4)=5
      myparent(2)=1
      myparent(3)=2
!      myparent(4)=3
      mychild(2)=3
      mychild(3)=4
!      mychild(4)=5
# elif defined SOLITON_REFINED
      nrefined(2)=5
      myparent(2)=1
      mychild(2)=3
# elif defined INLET_TEST
      nrefined(2)=5
      myparent(2)=1
      mychild(2)=3
# endif
#endif

        RETURN
        END SUBROUTINE initialize_param

      END MODULE mod_param
