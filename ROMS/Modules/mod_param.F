#include "cppdefs.h"
      MODULE mod_param
!
!svn $Id: mod_param.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
!  HaloBry    Buffers halo size for exchanging boundary arrays.        !
!  HaloSizeI  Maximum halo size, in grid points, in XI-direction.      !
!  HaloSizeJ  Maximum halo size, in grid points, in ETA-direction.     !
!  TileSide   Maximun tile side length in XI- or ETA-directions.       !
!  TileSize   Maximum tile size.                                       !
!                                                                      !
!  Configuration parameters:                                           !
!                                                                      !
!  Nbico      Number of balanced SSH elliptic equation iterations.     !
!  Nfloats    Number of floats trajectories.                           !
!  Nstation   Number of output stations.                               !
!  MTC        Maximum number of tidal components.                      !
!  MTY        Maximum number of tidal nodal correction years           !
!                                                                      !
!  State variables parameters:                                         !
!                                                                      !
#ifdef PROPAGATOR
!  Mstate     Size of FULL state vector (water points only).           !
!  Nstate     Size of NODE partition state vector: Nstate=Mstate in    !
!               serial applications.                                   !
!  Nstr       State vector node partition starting index.              !
!  Nend       State vector node partition ending   index.              !
!  Nsize      Size of the eigenvalue problem: Nend-Nstr+1.             !
# ifdef SO_SEMI
!  Nsemi      Number of time record saved in seminorm strochastic      !
!               optimals adjoint state vector.                         !
# endif
#endif
!  NSA        Number of state array for error covariance.              !
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
!  MT         Maximum number of tracer type variables.                 !
!  NT         Total number of tracer type variables.                   !
!  NTCLM      Total number of climatology tracer type variables to     !
!               process.                                               !
!                                                                      !
!  NBAND      Number of bands in the spectral irradiance model.        !
!  Nbed       Number of sediment bed layers.                           !
!  NCS        Number of cohesive (mud) sediment tracers.               !
!  NNS        Number of non-cohesive (sand) sediment tracers.          !
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
!  Grid nesting parameters.
!-----------------------------------------------------------------------
!
!  Number of nested and/or connected grids to solve.
!
      integer :: Ngrids
!
!  Number of grid nesting layers. This parameter allows refinement and
!  composite grid combinations.
!
      integer :: NestLayers = 1
!
!  Number of grids in each nested layer.
!
      integer, allocatable :: GridsInLayer(:)   ! [NestLayers]
!
!  Application grid number as a function of the maximum number of grids
!  in each nested layer and number of nested layers.
!
      integer, allocatable :: GridNumber(:,:)   ! [Ngrids,NestLayers]
!
!-----------------------------------------------------------------------
!  Define derived type structure, T_IO, used to store information about
!  the input and output files.
!-----------------------------------------------------------------------
!
!  This structure is used to declare the variables associated with
!  input and output files, like TYPE(IO) :: HIS(Ngrids). It is a
!  compact way to store a lot of information.
!
      TYPE T_IO
        integer :: Nfiles                        ! number of multi-files
        integer :: Fcount                        ! multi-file counter
        integer :: Rindex                        ! NetCDF record index
        integer :: ncid                          ! NetCDF file ID
        integer,  pointer :: Nrec(:)             ! NetCDF record size
        integer,  pointer :: Vid(:)              ! NetCDF variables IDs
        integer,  pointer :: Tid(:)              ! NetCDF tracers IDs
        real(r8), pointer :: time_min(:)         ! starting time
        real(r8), pointer :: time_max(:)         ! ending time
        character (len=50 ) :: label             ! structure label
        character (len=256) :: base              ! base file name
        character (len=256) :: name              ! current name
        character (len=256), pointer :: files(:) ! multi-file names
      END TYPE T_IO
!
!-----------------------------------------------------------------------
!  Lower and upper bounds indices per domain partition for all grids.
!-----------------------------------------------------------------------
!
!  Notice that these indices have different values in periodic and
!  nesting applications, and on tiles next to the boundaries. Special
!  indices are required to process overlap regions (suffices P and T)
!  lateral boundary conditions (suffices B and M) in nested grid
!  applications. The halo indices are used in private computations
!  which include ghost-points and are limited by MAX/MIN functions
!  on the tiles next to the  model boundaries. For more information
!  about all these indices, see routine "var_bounds" in file
!  "Utility/get_bounds.F".
!
!  All the 1D array indices are of size -1:NtileI(ng)*NtileJ(ng)-1. The
!  -1 index include the values for the full (no partitions) grid.
!
!  Notice that the starting (Imin, Jmin) and ending (Imax, Jmax) indices
!  for I/O processing are 3D arrays. The first dimension (1:4) is for
!  1=PSI, 2=RHO, 3=u, 4=v points; the second dimension (0:1) is number
!  of ghost points (0: no ghost points, 1: Nghost points), and the
!  the third dimension is for 0:NtileI(ng)*NtileJ(ng)-1.
!
      TYPE T_BOUNDS
        integer, pointer :: tile(:)  ! tile partition

        integer, pointer :: LBi(:)   ! lower bound I-dimension
        integer, pointer :: UBi(:)   ! upper bound I-dimension
        integer, pointer :: LBj(:)   ! lower bound J-dimension
        integer, pointer :: UBj(:)   ! upper bound J-dimension

        integer :: LBij              ! lower bound MIN(I,J)-dimension
        integer :: UBij              ! upper bound MAX(I,J)-dimension

        integer :: edge(4,4)         ! boundary edges I- or J-indices

        integer, pointer :: Istr(:)  ! starting tile I-direction
        integer, pointer :: Iend(:)  ! ending   tile I-direction
        integer, pointer :: Jstr(:)  ! starting tile J-direction
        integer, pointer :: Jend(:)  ! ending   tile J-direction

        integer, pointer :: IstrR(:) ! starting tile I-direction (RHO)
        integer, pointer :: IendR(:) ! ending   tile I-direction (RHO)
        integer, pointer :: IstrU(:) ! starting tile I-direction (U)

        integer, pointer :: JstrR(:) ! starting tile J-direction (RHO)
        integer, pointer :: JendR(:) ! ending   tile J-direction (RHO)
        integer, pointer :: JstrV(:) ! starting tile J-direction (V)

        integer, pointer :: IstrB(:) ! starting obc I-direction (RHO,V)
        integer, pointer :: IendB(:) ! ending   obc I-direction (RHO,V)
        integer, pointer :: IstrM(:) ! starting obc I-direction (PSI,U)

        integer, pointer :: JstrB(:) ! starting obc J-direction (RHO,U)
        integer, pointer :: JendB(:) ! ending   obc J-direction (RHO,U)
        integer, pointer :: JstrM(:) ! starting obc J-direction (PSI,V)

        integer, pointer :: IstrP(:) ! starting nest I-direction (PSI,U)
        integer, pointer :: IendP(:) ! ending   nest I-direction (PSI)
        integer, pointer :: JstrP(:) ! starting nest J-direction (PSI,V)
        integer, pointer :: JendP(:) ! ending   nest J-direction (PSI)

        integer, pointer :: IstrT(:) ! starting nest I-direction (RHO)
        integer, pointer :: IendT(:) ! ending   nest I-direction (RHO)
        integer, pointer :: JstrT(:) ! starting nest J-direction (RHO)
        integer, pointer :: JendT(:) ! ending   nest J-direction (RHO)

        integer, pointer :: Istrm3(:)    ! starting I-halo, Istr-3
        integer, pointer :: Istrm2(:)    ! starting I-halo, Istr-2
        integer, pointer :: Istrm1(:)    ! starting I-halo, Istr-1
        integer, pointer :: IstrUm2(:)   ! starting I-halo, IstrU-2
        integer, pointer :: IstrUm1(:)   ! starting I-halo, IstrU-1
        integer, pointer :: Iendp1(:)    ! ending   I-halo, Iend+1
        integer, pointer :: Iendp2(:)    ! ending   I-halo, Iend+2
        integer, pointer :: Iendp2i(:)   ! ending   I-halo, Iend+2
        integer, pointer :: Iendp3(:)    ! ending   I-halo, Iend+3

        integer, pointer :: Jstrm3(:)    ! starting J-halo, Jstr-3
        integer, pointer :: Jstrm2(:)    ! starting J-halo, Jstr-2
        integer, pointer :: Jstrm1(:)    ! starting J-halo, Jstr-1
        integer, pointer :: JstrVm2(:)   ! starting J-halo, JstrV-2
        integer, pointer :: JstrVm1(:)   ! starting J-halo, JstrV-1
        integer, pointer :: Jendp1(:)    ! ending   J-halo, Jend+1
        integer, pointer :: Jendp2(:)    ! ending   J-halo, Jend+2
        integer, pointer :: Jendp2i(:)   ! ending   J-halo, Jend+2
        integer, pointer :: Jendp3(:)    ! ending   J-halo, Jend+3

        integer, pointer :: Imin(:,:,:)  ! starting ghost I-direction
        integer, pointer :: Imax(:,:,:)  ! ending   ghost I-direction
        integer, pointer :: Jmin(:,:,:)  ! starting ghost J-direction
        integer, pointer :: Jmax(:,:,:)  ! ending   ghost J-direction
      END TYPE T_BOUNDS

      TYPE (T_BOUNDS), allocatable :: BOUNDS(:)
!
!-----------------------------------------------------------------------
!  Lower and upper bounds in NetCDF files.
!-----------------------------------------------------------------------
!
      TYPE T_IOBOUNDS
        integer :: ILB_psi       ! I-direction lower bound (PSI)
        integer :: IUB_psi       ! I-direction upper bound (PSI)
        integer :: JLB_psi       ! J-direction lower bound (PSI)
        integer :: JUB_psi       ! J-direction upper bound (PSI)

        integer :: ILB_rho       ! I-direction lower bound (RHO)
        integer :: IUB_rho       ! I-direction upper bound (RHO)
        integer :: JLB_rho       ! J-direction lower bound (RHO)
        integer :: JUB_rho       ! J-direction upper bound (RHO)

        integer :: ILB_u         ! I-direction lower bound (U)
        integer :: IUB_u         ! I-direction upper bound (U)
        integer :: JLB_u         ! J-direction lower bound (U)
        integer :: JUB_u         ! J-direction upper bound (U)

        integer :: ILB_v         ! I-direction lower bound (V)
        integer :: IUB_v         ! I-direction upper bound (V)
        integer :: JLB_v         ! J-direction lower bound (V)
        integer :: JUB_v         ! J-direction upper bound (V)

        integer :: IorJ          ! number of MAX(I,J)-direction points

        integer :: xi_psi        ! number of I-direction points (PSI)
        integer :: xi_rho        ! number of I-direction points (RHO)
        integer :: xi_u          ! number of I-direction points (U)
        integer :: xi_v          ! number of I-direction points (V)

        integer :: eta_psi       ! number of J-direction points (PSI)
        integer :: eta_rho       ! number of J-direction points (RHO)
        integer :: eta_u         ! number of I-direction points (U)
        integer :: eta_v         ! number of I-direction points (V)

#if (defined READ_WATER || defined WRITE_WATER) && defined MASKING
        integer :: xy_psi        ! number of IJ-water points (PSI)
        integer :: xy_rho        ! number of IJ-water points (RHO)
        integer :: xy_u          ! number of IJ-water points (U)
        integer :: xy_v          ! number of IJ-water points (V)
#endif
      END TYPE T_IOBOUNDS

      TYPE (T_IOBOUNDS), allocatable :: IOBOUNDS(:)
!
!-----------------------------------------------------------------------
!  Domain boundary edges switches and tiles minimum and maximum
!  fractional grid coordinates.
!-----------------------------------------------------------------------
!
      TYPE T_DOMAIN
        logical, pointer :: Eastern_Edge(:)
        logical, pointer :: Western_Edge(:)
        logical, pointer :: Northern_Edge(:)
        logical, pointer :: Southern_Edge(:)

        logical, pointer :: NorthEast_Corner(:)
        logical, pointer :: NorthWest_Corner(:)
        logical, pointer :: SouthEast_Corner(:)
        logical, pointer :: SouthWest_Corner(:)

        logical, pointer :: NorthEast_Test(:)
        logical, pointer :: NorthWest_Test(:)
        logical, pointer :: SouthEast_Test(:)
        logical, pointer :: SouthWest_Test(:)

        real(r8), pointer :: Xmin_psi(:)
        real(r8), pointer :: Xmax_psi(:)
        real(r8), pointer :: Ymin_psi(:)
        real(r8), pointer :: Ymax_psi(:)

        real(r8), pointer :: Xmin_rho(:)
        real(r8), pointer :: Xmax_rho(:)
        real(r8), pointer :: Ymin_rho(:)
        real(r8), pointer :: Ymax_rho(:)

        real(r8), pointer :: Xmin_u(:)
        real(r8), pointer :: Xmax_u(:)
        real(r8), pointer :: Ymin_u(:)
        real(r8), pointer :: Ymax_u(:)

        real(r8), pointer :: Xmin_v(:)
        real(r8), pointer :: Xmax_v(:)
        real(r8), pointer :: Ymin_v(:)
        real(r8), pointer :: Ymax_v(:)
      END TYPE T_DOMAIN

      TYPE (T_DOMAIN), allocatable :: DOMAIN(:)
!
!-----------------------------------------------------------------------
!  Lateral Boundary Conditions (LBC) switches structure.
!-----------------------------------------------------------------------
!
!  The lateral boundary conditions are specified by logical switches
!  to facilitate applications with nested grids. It also allows to
!  set different boundary conditions between the nonlinear model and
!  the adjoint/tangent models. The LBC structure is allocated as:
!
!       LBC(1:4, nLBCvar, Ngrids)
!
!  where 1:4 are the number boundary edges, nLBCvar are the number
!  LBC state variables, and Ngrids is the number of nested grids.
!  For example, for free-surface gradient boundary conditions we
!  have:
!
!       LBC(iwest,  isFsur, ng) % gradient
!       LBC(ieast,  isFsur, ng) % gradient
!       LBC(isouth, isFsur, ng) % gradient
!       LBC(inorth, isFsur, ng) % gradient
!
      integer :: nLBCvar
!
      TYPE T_LBC
        logical :: acquire        ! process lateral boundary data

        logical :: Chapman_explicit
        logical :: Chapman_implicit
        logical :: clamped
        logical :: closed
        logical :: Flather
        logical :: gradient
        logical :: mixed
        logical :: nested
        logical :: nudging
        logical :: periodic
        logical :: radiation
        logical :: reduced
        logical :: Shchepetkin
      END TYPE T_LBC

      TYPE (T_LBC), allocatable :: LBC(:,:,:)

#if defined ADJOINT || defined TANGENT || defined TL_IOMS
      TYPE (T_LBC), allocatable :: ad_LBC(:,:,:)
      TYPE (T_LBC), allocatable :: tl_LBC(:,:,:)
#endif

#ifdef PROPAGATOR
!
!-----------------------------------------------------------------------
!  Generalized Stability Theory (GST) nested grid pointers for the
!  state vector. Mostly used for memory management with the ARPACK
!  library.
!-----------------------------------------------------------------------
!
      TYPE T_GST
        real(r8), pointer :: vector(:)
      END TYPE
#endif
!
!-----------------------------------------------------------------------
!  Model grid(s) parameters.
!-----------------------------------------------------------------------
!
!  Number of interior RHO-points in the XI- and ETA-directions. The
!  size of models state variables (C-grid) at input and output are:
!
!    RH0-type variables:  [0:Lm+1, 0:Mm+1]        ----v(i,j+1)----
!    PSI-type variables:  [1:Lm+1, 1:Mm+1]        |              |
!      U-type variables:  [1:Lm+1, 0:Mm+1]     u(i,j)  r(i,j)  u(i+1,j)
!      V-type variables:  [0:Lm+1, 1:Mm+1]        |              |
!                                                 -----v(i,j)-----
      integer, allocatable :: Lm(:)
      integer, allocatable :: Mm(:)
!
!  Global horizontal size of model arrays including padding.  All the
!  model state arrays are of same size to facilitate parallelization.
!
      integer, allocatable :: Im(:)
      integer, allocatable :: Jm(:)
!
!  Number of vertical levels. The vertical ranges of model state
!  variables are:
!                                                 -----W(i,j,k)-----
!    RHO-, U-, V-type variables: [1:N]            |                |
!              W-type variables: [0:N]            |    r(i,j,k)    |
!                                                 |                |
!                                                 ----W(i,j,k-1)----
      integer, allocatable :: N(:)
!
!-----------------------------------------------------------------------
!  Tracers parameters.
!-----------------------------------------------------------------------
!
!  Total number of tracer type variables, NT(:) = NAT + NBT + NPT + NST.
!  The MT corresponds to the maximum number of tracers between all
!  nested grids.
!
      integer, allocatable :: NT(:)
      integer :: MT
!
!  Total number of climatology tracer type variables to process.
!
      integer, allocatable :: NTCLM(:)
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

#ifdef FILTERED_RST
# ifdef ICE_MODEL
        integer, parameter :: nfile = 10
# else
        integer, parameter :: nfile = 9
# endif
#endif
!
!  Number of biological tracers.
!
      integer :: NBT = 0
#if defined ESTUARYBGC && defined SPECTRAL_LIGHT
      integer :: NBAND = 0
#endif
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
      integer :: NCS = 0
!
!  Number of non-cohesive (sand) sediments.
!
      integer :: NNS = 0

#ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Floats trajectories parameters.
!-----------------------------------------------------------------------
!
!  Number of trajectory time-stepping levels, [0:NFT].
!
      integer, parameter :: NFT = 4
!
!  Total number of floats to track.
!
      integer, allocatable :: Nfloats(:)
!
!  Total number of float variables to process and output.
!
      integer, allocatable :: NFV(:)
#endif
#ifdef STATIONS
!
!-----------------------------------------------------------------------
!  Stations parameters.
!-----------------------------------------------------------------------
!
!  Number of output stations.
!
      integer, allocatable :: Nstation(:)
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
!  the Fennel model.
!
      integer :: NDbio2d                ! 2D fields
      integer :: NDbio3d                ! 3D fields
!
!  Number of diagnostic 3D right-hand-side fields.
!
      integer :: NDrhs
#endif

#if defined FOUR_DVAR || defined IMPULSE
!
!-----------------------------------------------------------------------
!  4DVAR parameters.
!-----------------------------------------------------------------------
!
!  Number of state arrays (size of additional dimension) for error
!  covariance normalization and standard deviation factors. This
!  paramenter is NSA=1 for strong constraint and NSA=2 for weak
!  constraint. In weak constraint, it contains fields for initial
!  conditions and model error covariance states.
!
      integer :: NSA
#endif

#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Number of biconjugate gradient algorithm iteractions.
!
      integer, allocatable :: Nbico(:)
#endif

#if defined TIME_CONV && defined WEAK_CONSTRAINT
!
!  Number of full state records needed for time convolutions. It is
!  assign to (ntimes(ng)/nADJ(ng))+1 in "read_phypar.F"
!
      integer, allocatable :: NrecTC(:)
#endif
!
!-----------------------------------------------------------------------
!  Model state parameters.
!-----------------------------------------------------------------------
!
!  Number of model state variables.
!
      integer, allocatable :: NSV(:)

#ifdef PROPAGATOR
!
!  Size of the full packed state vector (water points only), its
!  parallel node partition, and associated node starting and ending
!  indices.
!
      integer, allocatable :: Mstate(:)        ! full state
      integer, allocatable :: Nstate(:)        ! partioned state
      integer, allocatable :: Nstr(:)          ! starting index
      integer, allocatable :: Nend(:)          ! ending index
      integer, allocatable :: Nsize(:)         ! Nend-Nstr+1

# ifdef SO_SEMI
!
!  Number of time record saved in seminorm strochastic optimals adjoint
!  state vector.
!
      integer, allocatable :: Nsemi(:)
# endif
#endif
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
      integer, allocatable :: NtileI(:)
      integer, allocatable :: NtileJ(:)
!
!  Number of tiles or domain partitions in the XI- and ETA-directions.
!  These values are used to parallel loops to differentiate between
!  shared-memory and distributed-memory.  Notice that in distributed
!  memory both values are set to one.
!
      integer, allocatable :: NtileX(:)
      integer, allocatable :: NtileE(:)
!
!  Maximun number of points in the halo region for exchanging state
!  boundary arrays during convolutions.
!
      integer, allocatable :: HaloBry(:)
!
!  Maximum number of points in the halo region in the XI- and
!  ETA-directions.
!
      integer, allocatable :: HaloSizeI(:)
      integer, allocatable :: HaloSizeJ(:)
!
!  Maximum tile side length in XI- or ETA-directions.
!
      integer, allocatable :: TileSide(:)
!
!  Maximum number of points in a tile partition.
!
      integer, allocatable :: TileSize(:)
!
!  Set number of ghost-points in the halo region.  It is used in
!  periodic, nested, and distributed-memory applications.
!
      integer :: NghostPoints
!
      CONTAINS
!
      SUBROUTINE allocate_param
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Allocate dimension parameters.
!-----------------------------------------------------------------------
!
      allocate ( Lm(Ngrids) )
      allocate ( Mm(Ngrids) )

      allocate ( Im(Ngrids) )
      allocate ( Jm(Ngrids) )

      allocate ( N(Ngrids) )

      allocate ( NT(Ngrids) )

#if defined TIME_CONV && defined WEAK_CONSTRAINT
      allocate ( NrecTC(Ngrids) )
#endif

      allocate ( NTCLM(Ngrids) )

      allocate ( NSV(Ngrids) )

      allocate ( NtileI(Ngrids) )
      allocate ( NtileJ(Ngrids) )

      allocate ( NtileX(Ngrids) )
      allocate ( NtileE(Ngrids) )

      allocate ( HaloBry(Ngrids) )

      allocate ( HaloSizeI(Ngrids) )
      allocate ( HaloSizeJ(Ngrids) )

      allocate ( TileSide(Ngrids) )
      allocate ( TileSize(Ngrids) )

#ifdef PROPAGATOR
      allocate ( Mstate(Ngrids) )
      allocate ( Nstate(Ngrids) )
      allocate ( Nstr(Ngrids) )
      allocate ( Nend(Ngrids) )
      allocate ( Nsize(Ngrids) )
# ifdef SO_SEMI
      allocate ( Nsemi(Ngrids) )
# endif
#endif

#ifdef FLOATS
      allocate ( Nfloats(Ngrids) )
      allocate ( NFV(Ngrids) )
#endif

#ifdef STATIONS
      allocate ( Nstation(Ngrids) )
#endif

#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      allocate ( Nbico(Ngrids) )
#endif

      RETURN
      END SUBROUTINE allocate_param

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
      integer :: I_padd, J_padd, Ntiles
      integer :: ibry, ivar, ng
!
!-----------------------------------------------------------------------
!  Now that we know the values for the tile partitions (NtileI,NtileJ),
!  allocate module structures.
!-----------------------------------------------------------------------
!
!  Allocate lower and upper bounds indices structure.
!
      IF (.not.allocated(BOUNDS)) THEN
        allocate ( BOUNDS(Ngrids) )
        DO ng=1,Ngrids

          Ntiles=NtileI(ng)*NtileJ(ng)-1
          allocate ( BOUNDS(ng) % tile(-1:Ntiles) )

          allocate ( BOUNDS(ng) % LBi(-1:Ntiles) )
          allocate ( BOUNDS(ng) % UBi(-1:Ntiles) )
          allocate ( BOUNDS(ng) % LBj(-1:Ntiles) )
          allocate ( BOUNDS(ng) % UBj(-1:Ntiles) )

          allocate ( BOUNDS(ng) % Istr(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Iend(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jstr(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jend(-1:Ntiles) )

          allocate ( BOUNDS(ng) % IstrR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrU(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrV(-1:Ntiles) )

          allocate ( BOUNDS(ng) % IstrB(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendB(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrM(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrB(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendB(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrM(-1:Ntiles) )

          allocate ( BOUNDS(ng) % IstrP(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendP(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrP(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendP(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendT(-1:Ntiles) )

          allocate ( BOUNDS(ng) % Istrm3(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Istrm2(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Istrm1(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrUm2(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrUm1(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Iendp1(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Iendp2(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Iendp2i(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Iendp3(-1:Ntiles) )

          allocate ( BOUNDS(ng) % Jstrm3(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jstrm2(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jstrm1(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrVm2(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrVm1(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jendp1(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jendp2(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jendp2i(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jendp3(-1:Ntiles) )

          allocate ( BOUNDS(ng) % Imin(4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % Imax(4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % Jmin(4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % Jmax(4,0:1,0:Ntiles) )
        END DO
      END IF
!
!  DOMAIN structure containing boundary edges switches and fractional
!  grid lower/upper bounds for each tile.
!
      IF (.not.allocated(DOMAIN)) THEN
        allocate ( DOMAIN(Ngrids) )
        DO ng=1,Ngrids
          Ntiles=NtileI(ng)*NtileJ(ng)-1

          allocate ( DOMAIN(ng) % Eastern_Edge(-1:Ntiles) )
          allocate ( DOMAIN(ng) % Western_Edge(-1:Ntiles) )
          allocate ( DOMAIN(ng) % Northern_Edge(-1:Ntiles) )
          allocate ( DOMAIN(ng) % Southern_Edge(-1:Ntiles) )

          allocate ( DOMAIN(ng) % NorthEast_Corner(-1:Ntiles) )
          allocate ( DOMAIN(ng) % NorthWest_Corner(-1:Ntiles) )
          allocate ( DOMAIN(ng) % SouthEast_Corner(-1:Ntiles) )
          allocate ( DOMAIN(ng) % SouthWest_Corner(-1:Ntiles) )

          allocate ( DOMAIN(ng) % NorthEast_Test(-1:Ntiles) )
          allocate ( DOMAIN(ng) % NorthWest_Test(-1:Ntiles) )
          allocate ( DOMAIN(ng) % SouthEast_Test(-1:Ntiles) )
          allocate ( DOMAIN(ng) % SouthWest_Test(-1:Ntiles) )

          allocate ( DOMAIN(ng) % Xmin_psi(0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmax_psi(0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymin_psi(0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymax_psi(0:Ntiles) )

          allocate ( DOMAIN(ng) % Xmin_rho(0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmax_rho(0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymin_rho(0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymax_rho(0:Ntiles) )

          allocate ( DOMAIN(ng) % Xmin_u(0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmax_u(0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymin_u(0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymax_u(0:Ntiles) )

          allocate ( DOMAIN(ng) % Xmin_v(0:Ntiles) )
          allocate ( DOMAIN(ng) % Xmax_v(0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymin_v(0:Ntiles) )
          allocate ( DOMAIN(ng) % Ymax_v(0:Ntiles) )
        END DO
      END IF
!
!  Allocate lower and upper bounds structure for I/O NetCDF files.
!
      allocate ( IOBOUNDS(Ngrids) )

#ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  Determine number of diagnostic variables.
!-----------------------------------------------------------------------

# ifdef DIAGNOSTICS_TS
!
!  Tracer diagnostics.
!
      NDT=6          ! Acceleration, advection, vertical diffusion
#  if defined TS_DIF2 || defined TS_DIF4
      NDT=NDT+3      ! Horizontal (total, X-, Y-) diffusion
#   if defined MIX_GEO_TS || defined MIX_ISO_TS
      NDT=NDT+1      ! Horizontal S-diffusion due to rotated tensor
#   endif
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
      NDM2d=NDM2d+3  ! Horizontal total-, X-, and Y-advection
#   endif
#   ifdef WEC_MELLOR
      NDM2d=NDM2d+1  ! Horizontal radiation stresses
#   endif
#   ifdef WEC_VF
#    ifdef UV_COR
      NDM2d=NDM2d+1  ! Coriolis
#    endif
#    ifdef BOTTOM_STREAMING
      NDM2d=NDM2d+1  ! Bottom streaming
#    endif
#    ifdef SURFACE_STREAMING
      NDM2d=NDM2d+1  ! Surface streaming
#    endif
      NDM2d=NDM2d+8  ! zeta, zetaw, zqsp, zbeh (all 4 = prsgrd),
                       ! K- and Horiz- VF, roller, break
#   endif
#   if defined VEGETATION && defined VEG_DRAG 
      NDM2d=NDM2d+1  ! Resistance imposed on the flow by vegetation
#   endif
#   ifdef UV_COR
      NDM2d=NDM2d+1  ! Coriolis
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
      NDM2d=NDM2d+3  ! Horizontal total-, X-, and Y-viscosity
#   endif
#  ifdef SOLVE3D
!
!  3D Momentum diagnostics and right-hand-side terms.
!
      NDM3d=3        ! Acceleration, 3D P-Gradient, vertical viscosity
      NDrhs=1        ! 3D P-Gradient
#   ifdef UV_ADV
      NDM3d=NDM3d+4  ! Horizontal (total, X, Y) and vertical advection
      NDrhs=NDrhs+4
#   endif
#   ifdef WEC_MELLOR
      NDM3d=NDM3d+2  ! Horizontal and vertical radiation stresses
      NDrhs=NDrhs+2
#   endif
#   ifdef WEC_VF
#    ifdef UV_COR
      NDM3d=NDM3d+1  ! Coriolis
      NDrhs=NDrhs+1
#    endif
#    ifdef BOTTOM_STREAMING
      NDM3d=NDM3d+1  ! Bottom streaming
      NDrhs=NDrhs+1
#    endif
#    ifdef SURFACE_STREAMING
      NDM3d=NDM3d+1  ! Surface streaming
      NDrhs=NDrhs+1
#    endif
      NDM3d=NDM3d+5  ! wave break, roller, K- Horiz and Jvert- VF
      NDrhs=NDrhs+5
#   endif
#   if defined VEGETATION && defined VEG_DRAG
      NDM3d=NDM3d+1  ! Resistance imposed on the flow by vegetation
      NDrhs=NDrhs+1
#   endif
#   ifdef UV_COR
      NDM3d=NDM3d+1  ! Coriolis
      NDrhs=NDrhs+1
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
      NDM3d=NDM3d+3  ! Horizontal (total, X, Y) viscosity
#   endif
#   ifdef BODYFORCE
      NDrhs=NDrhs+1  ! Vertical viscosity
#   endif
#  else
      NDM3d=0        ! No 3D momentum diagnostics
      NDrhs=0
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
# ifdef FLOAT_OYSTER
#  ifdef FLOAT_VWALK
        NFV(ng)=NT(ng)+18
#  else
        NFV(ng)=NT(ng)+16
#  endif
# else
#  ifdef FLOAT_VWALK
        NFV(ng)=NT(ng)+12
#  else
        NFV(ng)=NT(ng)+10
#  endif
# endif
#ifdef OFFLINE_FLOATS_LATLON
! add 4 variables: x-space, y-space, water depth and Hz 
          NFV(ng)=NFV(ng)+4
#endif
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
#if defined IS4DVAR || defined W4DVAR || defined W4DPSAS
      MT=MAX(NAT,MAXVAL(NT))
#else
      MT=MAX(2,MAXVAL(NT))
#endif
!
!-----------------------------------------------------------------------
!  Allocate Lateral Boundary Conditions switches structure.
!-----------------------------------------------------------------------
!
#ifdef SOLVE3D
      nLBCvar=5+MT
# if defined GLS_MIXING || defined MY25_MIXING
      nLBCvar=nLBCvar+1
# endif
# ifdef ICE_MODEL
      nLBCvar=nLBCvar+11
# endif
# ifdef WEC
      nLBCvar=nLBCvar+4
# endif
# ifdef INWAVE_MODEL
      nLBCvar=nLBCvar+4
# endif
#else
      nLBCvar=3
#endif
!
      IF (.not.allocated(LBC)) THEN
        allocate ( LBC(4,nLBCvar,Ngrids) )
        DO ng=1,Ngrids
          DO ivar=1,nLBCvar
            DO ibry=1,4
              LBC(ibry,ivar,ng)%acquire = .FALSE.
              LBC(ibry,ivar,ng)%Chapman_explicit = .FALSE.
              LBC(ibry,ivar,ng)%Chapman_implicit = .FALSE.
              LBC(ibry,ivar,ng)%clamped = .FALSE.
              LBC(ibry,ivar,ng)%closed = .FALSE.
              LBC(ibry,ivar,ng)%Flather = .FALSE.
              LBC(ibry,ivar,ng)%gradient = .FALSE.
              LBC(ibry,ivar,ng)%mixed    =.FALSE.
              LBC(ibry,ivar,ng)%nested = .FALSE.
              LBC(ibry,ivar,ng)%nudging = .FALSE.
              LBC(ibry,ivar,ng)%periodic = .FALSE.
              LBC(ibry,ivar,ng)%radiation = .FALSE.
              LBC(ibry,ivar,ng)%Shchepetkin = .FALSE.
            END DO
          END DO
        END DO
      END IF
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
      IF (.not.allocated(ad_LBC)) THEN
        allocate ( ad_LBC(4,nLBCvar,Ngrids) )
        DO ng=1,Ngrids
          DO ivar=1,nLBCvar
            DO ibry=1,4
              ad_LBC(ibry,ivar,ng)%acquire = .FALSE.
              ad_LBC(ibry,ivar,ng)%Chapman_implicit = .FALSE.
              ad_LBC(ibry,ivar,ng)%Chapman_explicit = .FALSE.
              ad_LBC(ibry,ivar,ng)%clamped = .FALSE.
              ad_LBC(ibry,ivar,ng)%closed = .FALSE.
              ad_LBC(ibry,ivar,ng)%Flather = .FALSE.
              ad_LBC(ibry,ivar,ng)%gradient = .FALSE.
              ad_LBC(ibry,ivar,ng)%nested = .FALSE.
              ad_LBC(ibry,ivar,ng)%nudging = .FALSE.
              ad_LBC(ibry,ivar,ng)%periodic = .FALSE.
              ad_LBC(ibry,ivar,ng)%radiation = .FALSE.
              ad_LBC(ibry,ivar,ng)%reduced = .FALSE.
              ad_LBC(ibry,ivar,ng)%Shchepetkin = .FALSE.
            END DO
          END DO
        END DO
      END IF
#endif
#if defined TANGENT || defined TL_IOMS
      IF (.not.allocated(tl_LBC)) THEN
        allocate ( tl_LBC(4,nLBCvar,Ngrids) )
        DO ng=1,Ngrids
          DO ivar=1,nLBCvar
            DO ibry=1,4
              tl_LBC(ibry,ivar,ng)%acquire = .FALSE.
              tl_LBC(ibry,ivar,ng)%Chapman_implicit = .FALSE.
              tl_LBC(ibry,ivar,ng)%Chapman_explicit = .FALSE.
              tl_LBC(ibry,ivar,ng)%clamped = .FALSE.
              tl_LBC(ibry,ivar,ng)%closed = .FALSE.
              tl_LBC(ibry,ivar,ng)%Flather = .FALSE.
              tl_LBC(ibry,ivar,ng)%gradient = .FALSE.
              tl_LBC(ibry,ivar,ng)%nested = .FALSE.
              tl_LBC(ibry,ivar,ng)%nudging = .FALSE.
              tl_LBC(ibry,ivar,ng)%periodic =.FALSE.
              tl_LBC(ibry,ivar,ng)%radiation = .FALSE.
              tl_LBC(ibry,ivar,ng)%reduced = .FALSE.
              tl_LBC(ibry,ivar,ng)%Shchepetkin = .FALSE.
            END DO
          END DO
        END DO
      END IF
#endif

      RETURN
      END SUBROUTINE initialize_param

      END MODULE mod_param
