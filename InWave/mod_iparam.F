#include "cppdefs.h"
      MODULE mod_param
!
!svn $Id: mod_param.F 838 2008-11-17 04:22:18Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
!  Configuration parameters:                                           !
!                                                                      !
!  Nfloats    Number of floats tracjectories.                          !
!  Nstation   Number of output stations.                               !
!  MTC        Maximum number of tidal components.                      !
!  NSV        Number of model state variables.                         !
!                                                                      !
!  Tracer parameters:                                                  !
!                                                                      !
!  NAT        Number of active tracer type variables (usually,         !
!               NAT=2 for potential temperature and salinity).         !
!  NBT        Number of biological tracer type variables.              !
!  NST        Number of sediment tracer type variables (NCS+NNS).      !
!  NPT        Number of extra passive tracer type variables to         !
!               advect and diffuse only (dyes, etc).                   !
!  NT         Total number of tracer type variables.                   !
!  MT         Maximum number of tracer type variables.                 !
!                                                                      !
!  Nbed       Number of sediment bed layers.                           !
!  NCS        Number of cohesive (mud) sediment tracers.               !
!  NNS        Number of non-cohesive (sand) sediment tracers.          !
!                                                                      !
!  NBands     Number of spectral irradiance bands in EcoSim.           !
!  Nbac       Number of bacteria constituents in EcoSim.               !
!  Ndom       Number of dissolved matter constituents in EcoSim.       !
!  Nfec       Number of fecal matter constituents in EcoSim.           !
!  Nphy       Number of phytoplankton constituents in EcoSim.          !
!  Npig       Number of pigment constituents in EcoSim.                !
!  PHY        EcoSim indices of phytoplankton species considered.      !
!  PIG        EcoSim phytoplankton-pigment matrix.                     !
!                                                                      !
!  Diagnostic fields parameters:                                       !
!                                                                      !
!  NDbio2d    Number of diagnostic 2D biology fields.                  !
!  NDbio3d    Number of diagnostic 3D biology fields.                  !
!  NDT        Number of diagnostic tracer fields.                      !
!  NDM2d      Number of diagnostic 2D momentum fields.                 !
!  NDM3d      Number of diagnostic 3D momentum fields.                 !
!  NDrhs      Number of diagnostic 3D right-hand-side fields.          !
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
!  Tracers parameters.
!-----------------------------------------------------------------------
!
!  Total number of tracer type variables, NT(:) = NAT + NBT + NPT + NST.
!  The MT corresponds to the maximum number of tracers between all
!  nested grids.
!
        integer, dimension(Ngrids) :: NT
        integer :: MT
!
!  Number of active tracers. Usually, NAT=2 for potential temperature
!  and salinity.
!
        integer :: NAT = 0
!
!  Total number of inert passive tracers to advect and diffuse only
!  (like dyes, etc). This parameter is independent of the number of
!  biological and/or sediment tracers.
!
        integer :: NPT = 0
!
!-----------------------------------------------------------------------
!  Sediment tracers parameters.
!-----------------------------------------------------------------------
!
!  Number of sediment bed layes.
!
        integer :: Nbed = 0
!
!  Total number of sediment tracers, NST = NCS + NNS.
!
        integer :: NST = 0
!
!  Number of cohesive (mud) sediments.
!
        integer :: NCS                
!
!  Number of non-cohesive (sand) sediments.
!
        integer :: NNS

#ifdef BIOLOGY
!
!-----------------------------------------------------------------------
!  Biological tracers parameters.
!-----------------------------------------------------------------------

# if defined BIO_FASHAM
!
!  Number of tracers for Fasham-type, Nitrogen-based biological model.
!
#  ifdef CARBON
#    ifdef OXYGEN
        integer, parameter :: NBT = 12
#    else
        integer, parameter :: NBT = 11
#    endif
#  else
#    ifdef OXYGEN
        integer, parameter :: NBT = 8
#    else
        integer, parameter :: NBT = 7
#    endif
#  endif
# elif defined BIO_LIMADONEY
!
!  Number of tracers for Lima and Doney (2004) ecological model.
!
        integer, parameter :: NBT = 14
# elif defined NPZD_FRANKS || defined NPZD_POWELL
!
!  Number of tracers for NPZD model.
!
        integer, parameter :: NBT = 4
# elif defined ECOSIM
!
!  Bio-optical EcoSim parameters.
!
        integer, parameter :: NBands = 60    ! spectral bands
        integer, parameter :: Nbac = 1       ! bacteria constituents
        integer, parameter :: Ndom = 2       ! DOM constituents
        integer, parameter :: Nfec = 2       ! Fecal constituents
        integer, parameter :: Nphy = 4       ! Phytoplankton groups
        integer, parameter :: Npig = 7       ! Pigments
!
!  Determine number of EcoSim biological tracer. Currently, there is a
!  maximum of seven phytoplankton species and seven different pigments:
!
! [1] small diatom           [1] chlorophyll-a
! [2] large diatom           [2] chlorophyll-b
! [3] small dinoflagellate   [3] chlorophyll-c
! [4] large dinoflagellate   [4] photosythetic carotenoids
! [5] synechococcus          [5] photoprotective carotenoids
! [6] small prochlorococcus  [6] low  urobilin phycoeurythin carotenoids
! [7] large prochlorococcus  [7] high urobilin phycoeurythin carotenoids
!
!  The phytoplankton/pigment matrix is as follows:
!
!               P h y t o p l a n k t o n
!              [1]   [2]   [3]   [4]   [5]   [6]   [7]
!
!       t [7]   0     0     0     0     1     0     0
!       n [6]   0     0     0     0     0     0     0
!       e [5]   1     1     1     1     1     1     1
!       m [4]   1     1     1     1     0     0     0
!       g [3]   1     1     1     1     0     0     0
!       i [2]   0     0     0     0     0     1     1
!       P [1]   1     1     1     1     1     1     1
!
        integer, parameter, dimension(7,7) :: PIG = reshape (           &
     &                                      (/ 1, 1, 1, 1, 1, 1, 1,     &
     &                                         0, 0, 0, 0, 0, 1, 1,     &
     &                                         1, 1, 1, 1, 0, 0, 0,     &
     &                                         1, 1, 1, 1, 0, 0, 0,     &
     &                                         1, 1, 1, 1, 1, 1, 1,     &
     &                                         0, 0, 0, 0, 0, 0, 0,     &
     &                                         0, 0, 0, 0, 1, 0, 0 /),  &
     &                                      (/ 7, 7 /) )
!
!  Set phytoplankton species to consider (see above classification):
!
        integer, parameter, dimension(Nphy) :: PHY = (/ 1, 2, 4, 5 /)
!
!  Number of biological tracer type variables is set at run time
!
        integer :: NBT
# elif defined NEMURO
!
!  Number of tracers for the Nemuro lower trophic level ecosystem model.
!
        integer, parameter :: NBT = 11
# endif
#else
        integer, parameter :: NBT = 0
#endif
#ifdef FLOATS
! 
!-----------------------------------------------------------------------
!  Floats tracjectories parameters.
!-----------------------------------------------------------------------
!
!  Number of trajectory time-stepping levels, [0:NFT].
!
        integer, parameter :: NFT = 4
!
!  Total number of floats to track.
!
        integer, dimension(Ngrids) :: Nfloats
!
!  Total number of float variables to process and output.
!
        integer, dimension(Ngrids) :: NFV
#endif
#ifdef STATIONS
!
!-----------------------------------------------------------------------
!  Stations parameters.
!-----------------------------------------------------------------------
!
!  Number of output stations.
!
        integer, dimension(Ngrids) :: Nstation
#endif
!
!-----------------------------------------------------------------------
!  Maximum number of tidal constituents to process.
!-----------------------------------------------------------------------
!
        integer :: MTC

#ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  Diagnostic fields parameters.
!-----------------------------------------------------------------------
!
!  Number of diagnostic tracer fields.
!
        integer :: NDT
!
!  Number of diagnostic momentum fields.
!
        integer :: NDM2d                  ! 2D momentum
        integer :: NDM3d                  ! 3D momentum
!
!  Number of diagnostic biology fields.  Currenly, only available for
!  the Fasham model.
!
        integer :: NDbio2d                ! 2D fields
        integer :: NDbio3d                ! 3D fields
!
!  Number of diagnostic 3D right-hand-side fields.
!
        integer :: NDrhs
#endif
!
!-----------------------------------------------------------------------
!  Model state parameters.
!-----------------------------------------------------------------------
!
!  Number of model state variables.
!
        integer, dimension(Ngrids) :: NSV
!
!  Set nonlinear, tangent linear, and adjoint models identifiers.
!
        integer :: iNLM = 1
        integer :: iTLM = 2
        integer :: iRPM = 3
        integer :: iADM = 4
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
        integer :: I_padd, J_padd, ng
#ifdef ECOSIM
        integer :: i, j
!
!-----------------------------------------------------------------------
!  Determine number of EcoSim total bio-optical constituents:
!-----------------------------------------------------------------------
!
!       Nutrients: NO3, NO4, PO4, FeO, SiO, DIC  (6)
!        Bacteria: C, Fe, N, P                   (Nbac*4)
!             DOM: CDM, C, N, P                  (Ndom*4)
!           Fecal: C, Fe, N, P, Si               (Nfec*5)
!    Phytoplakton: C, Fe, N, P                   (Nfec*4 + Si)
!        Pigments: look table
!
        NBT=6+(Nbac*4)+(Ndom*4)+(Nfec*5)+(Nphy*4)
!
!  Add phytoplankton silica constituents.
!
        DO i=1,Nphy
          IF (PHY(i).le.2) NBT=NBT+1
        END DO
!
!  Add pigments.  Check phytoplankton-pigment table for values greater
!  than zero.
!
        DO j=1,Npig
          DO i=1,Nphy
            IF (PIG(PHY(i),j).eq.1) NBT=NBT+1
          END DO
        END DO
#endif

#ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  Determine number of diagnostic variables.
!-----------------------------------------------------------------------

# ifdef DIAGNOSTICS_TS
!
!  Tracer diagnostics.
!
        NDT=4          ! Acceleration, advection, vertical diffusion
#  if defined TS_DIF2 || defined TS_DIF4
        NDT=NDT+1      ! Horizontal diffusion
#  endif
# else
        NDT=0          ! No tracer diagnostics
# endif
# ifdef DIAGNOSTICS_UV
!
!  2D Momentum diagnostics.
!
        NDM2d=4        ! Acceleration, 2D P-Gradient, stresses
#   ifdef UV_ADV
        NDM2d=NDM2d+1  ! Horizontal advection
#   endif
#   ifdef NEARSHORE_MELLOR
        NDM2d=NDM2d+1  ! Horizontal radiation stresses
#   endif
#   ifdef UV_COR
        NDM2d=NDM2d+1  ! Coriolis
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
        NDM2d=NDM2d+1  ! Horizontal viscosity
#   endif
#  ifdef SOLVE3D
!
!  3D Momentum diagnostics and right-hand-side terms.
!
        NDM3d=3        ! Acceleration, 3D P-Gradient, vertical viscosity
        NDrhs=1        ! 3D P-Gradient
#   ifdef UV_ADV
        NDM3d=NDM3d+2  ! Horizontal and vertical advection
        NDrhs=NDrhs+2  
#   endif
#   ifdef NEARSHORE_MELLOR
        NDM3d=NDM3d+2  ! Horizontal and vertical radiation stresses
        NDrhs=NDrhs+2  
#   endif
#   ifdef UV_COR
        NDM3d=NDM3d+1  ! Coriolis
        NDrhs=NDrhs+1
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
        NDM3d=NDM3d+1  ! Horizontal viscosity
#   endif
#   ifdef BODYFORCE
        NDrhs=NDrhs+1  ! Vertical viscosity
#   endif
#  else
        NDM3d=0        ! No 3D momentum diagnostics
        NDrhs=0
#  endif
# endif
# if defined BIO_FASHAM && defined DIAGNOSTICS_BIO
!
!  Source and sink biology diagnostic terms.
!
        NDbio3d=2
        NDbio2d=0
#  ifdef DENITRIFICATION
        NDbio2d=NDbio2d+1
#  endif
#  ifdef CARBON
        NDbio2d=NDbio2d+2
#  endif
#  ifdef OXYGEN
        NDbio2d=NDbio2d+1
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Derived dimension parameters.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          I_padd=(Lm(ng)+2)/2-(Lm(ng)+1)/2
          J_padd=(Mm(ng)+2)/2-(Mm(ng)+1)/2
          Im(ng)=Lm(ng)+I_padd
          Jm(ng)=Mm(ng)+J_padd
          NT(ng)=NAT+NBT+NST+NPT
#ifdef FLOATS
# ifdef FLOAT_VWALK
          NFV(ng)=NT(ng)+12
# else
          NFV(ng)=NT(ng)+10
# endif
#endif
          NSV(ng)=NT(ng)+5
#ifdef ADJUST_WSTRESS
          NSV(ng)=NSV(ng)+2
#endif
#ifdef ADJUST_STFLUX
          NSV(ng)=NSV(ng)+NT(ng)
#endif
        END DO
!
!  Set maximum number of tracer between all nested grids.
!
        MT=MAX(2,MAXVAL(NT))
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
