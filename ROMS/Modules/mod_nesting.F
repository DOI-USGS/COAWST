#include "cppdefs.h"
      MODULE mod_nesting

#ifdef NESTING
!
!git $Id$
!svn $Id: mod_nesting.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module defines structures for composite and refinement grids.  !
!                                                                      !
!  Composite Grids Structure: Donor grid data at contact points        !
!  =========================                                           !
!                                                                      !
!  bustr     Kinematic bottom momentum flux (bottom stress) in the     !
!              XI-direction (m2/s2)  at U-points.                      !
!  bvstr     Kinematic bottom momentum flux (bottom stress) in the     !
!              ETA-direction (m2/s2) at V-points.                      !
!  rzeta     Right-hand-side of free surface equation (m3/s).          !
!  ubar      Vertically integrated U-momentum component (m/s).         !
!  vbar      Vertically integrated V-momentum component (m/s).         !
!  zeta      Free surface (m).                                         !
# ifdef SOLVE3D
!                                                                      !
!  DU_avg1   Time averaged U-flux for 2D equations (m3/s).             !
!  DV_avg1   Time averaged V-flux for 2D equations (m3/s).             !
!  Huon      Total U-momentum flux term, Hz*u/pn.                      !
!  Hvom      Total V-momentum flux term, Hz*v/pm.                      !
!  Zt_avg1   Free-surface averaged over all short time-steps (m).      !
!  t         Tracer type variables (active and passive).               !
!  u         3D U-momentum component (m/s).                            !
!  v         3D U-momentum component (m/s).                            !
# endif
!                                                                      !
!  REFINED Grids Structure: Donor grid data at contact points          !
!  =======================  (two-time rolling snapshots)               !
!                                                                      !
!  ubar      Vertically integrated U-momentum component (m/s).         !
!  vbar      Vertically integrated V-momentum component (m/s).         !
!  zeta      Free surface (m).                                         !
# ifdef SOLVE3D
!                                                                      !
!  DU_avg2   Time averaged U-flux for 3D equations coupling (m3/s).    !
!  DV_avg2   Time averaged V-flux for 3D equations coupling (m3/s).    !
!  t         Tracer type variables (active and passive).               !
!  u         3D U-momentum component (m/s).                            !
!  v         3D U-momentum component (m/s).                            !
# endif
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PUBLIC :: allocate_nesting
      PUBLIC :: deallocate_nesting
      PUBLIC :: initialize_nesting
!
!-----------------------------------------------------------------------
!  Nesting identification index of variables to process.
!-----------------------------------------------------------------------
!
!  The following identification indices are used in "initial" or
!  "main2d/main3d" to specify the variables that are processed in
!  each sub-timestep section.  Negative indices are used in grid
!  refinement whereas positive indices are used in composite grids.
!
      integer, parameter :: nmflx = -6    ! check mass flux conservation
      integer, parameter :: ndxdy = -5    ! extract on_u and om_v
      integer, parameter :: ngetD = -4    ! extract donor grid data
      integer, parameter :: nmask = -3    ! scale interpolation weights
      integer, parameter :: nputD = -2    ! fill contact points
      integer, parameter :: n2way = -1    ! fine to course coupling
!
      integer, parameter :: nFSIC =  1    ! free surface initialization
      integer, parameter :: n2dIC =  2    ! 2D momentum initialization
      integer, parameter :: n3dIC =  3    ! 3D momentum initialization
      integer, parameter :: nTVIC =  4    ! tracers initialization
      integer, parameter :: nbstr =  5    ! bottom stress (bustr,bvstr)
      integer, parameter :: nrhst =  6    ! RHS terms (tracers)
      integer, parameter :: nzeta =  7    ! 3D kernel free-surface
      integer, parameter :: nzwgt =  8    ! 3D vertical weights
      integer, parameter :: n2dPS =  9    ! 2D engine Predictor Step
      integer, parameter :: n2dCS = 10    ! 2D engine Corrector Step
      integer, parameter :: n2dfx = 11    ! time-averaged 2D fluxes
      integer, parameter :: n3duv = 12    ! 3D momentum and fluxes
      integer, parameter :: n3dTV = 13    ! 3D tracer variables
!
!-----------------------------------------------------------------------
!  Nesting parameters.
!-----------------------------------------------------------------------
!
!  Nested grid connectivity switches. It is used to determine the
!  dimensions of the numerical kernel allocatable arrays. The arrays
!  have extra points due to the contact regions in any of sides
!  of the physical grid (1=iwest, 2=isouth, 3=ieast, 4=inorth).
!
      logical, allocatable :: ContactRegion(:,:)    ! [4,Ngrids]
!
!  Logical switch indicating which coarser grid is a donor to a
!  finer receiver grid (RefineScale(rg) > 0) external contact points.
!  This switch is in terms of the donor coarser grid.
!
      logical, allocatable :: DonorToFiner(:)       ! {Ngrids]
!$OMP THREADPRIVATE (DonorToFiner)
!
!  Switch indicating which refined grid(s), with RefineScale(ng) > 0,
!  include finer refined grids inside: telescoping refinement.
!
      logical, allocatable :: Telescoping(:)        ! [Ngrids]
!$OMP THREADPRIVATE (Telescoping)
!
!  Switch to compute depth-dependent, vertical interpolation weights.
!  Currently, vertical weights are used in composite grids because
!  their grids are not coincident. They are not needed in refinement
!  grids because the donor and receiver grids have the same number of
!  vertical levels and have matching bathymetry.  However, in the
!  future, it is possible to have configurations that require vertical
!  weights in refinement. The switch "get_Vweights" controls if such
!  weights are computed or not. If false, it will accelerate the
!  computations because of less distributed-memory communications.
!
      logical :: get_Vweights
!$OMP THREADPRIVATE (get_Vweights)
!
!  If refinement, it contains the coarser donor grid number to finer
!  receiver grid external contact points.  The donor grid is always
!  coarser that receiver grid.  This variable is in terms of the
!  finer receiver grid.
!
      integer, allocatable :: CoarserDonor(:)       ! [Ngrids]
!$OMP THREADPRIVATE (CoarserDonor)
!
!  If refinement and two-way exchange, it contains the donor finer
!  grid number to a coarser receiver grid.  This variable is in
!  terms of the coarser donor grid.
!
      integer, allocatable :: FinerDonor(:)        ! [Ngrids]
!$OMP THREADPRIVATE (FinerDonor)
!
!  Number of refined time-steps. In most cases, the number of refined
!  time-step is the same as the refinement scale ratio for numerical
!  stability. However, the user is allowed to take larger divisible
!  time-step with respect to the donor grid. The variable below is
!  computed donor and receiver time-step ratio from standard input.
!  It is up to the user to determine the appropiate time-step for
!  stability.
!
      integer, allocatable :: RefineSteps(:)        ! [Ngrids]
!$OMP THREADPRIVATE (RefineSteps)
!
!  Refined time-steps counter with respect the coarse grid (ng=1)
!  single time-step.
!
      integer, allocatable :: RefineStepsCounter(:) ! [Ngrids]
!
!  Interval used in the two-way exchange between fine and coarse
!  grids.
!
      real(r8), allocatable :: TwoWayInterval(:)    ! [Ngrids]
!
!  Donor and reciver grids for each contact region. These paremeters
!  are also duplicated in the T_NGC structure.
!
      integer, allocatable :: donor_grid(:)         ! [Ncontact]
      integer, allocatable :: receiver_grid(:)      ! [Ncontact]
!
!  Rolling index and time (seconds) used in the temporal interpolation
!  of contact point data.
!
      integer,  allocatable :: RollingIndex(:)      ! [Ncontact]
      real(dp), allocatable :: RollingTime(:,:)     ! [Ncontact]
!$OMP THREADPRIVATE (RollingIndex, RollingTime)
!
!  If refinement, donor grid (I,J) indices at PSI points used to extract
!  refined grid. Values are set to -999 if not applicable.
!
!          +---------+    J_top
!          |         |
!          | Refined |
!          |    grid |
!          |         |
!          +---------+    J_bottom
!       I_left   I_right
!
      integer, allocatable :: I_left(:)             ! [Ngrids]
      integer, allocatable :: I_right(:)            ! [Ngrids]
      integer, allocatable :: J_bottom(:)           ! [Ngrids]
      integer, allocatable :: J_top(:)              ! [Ngrids]
!
!  Compact arrays used to unpack data from nested grids contact points
!  NetCDF file. They are allocated to the size "datum" dimension in
!  routine "set_contact".  The start and end indices for each C-type
!  variable are used to unpack from compact vector.
!
      integer :: NCdatum
      integer, allocatable :: NCpoints(:)           ! [Ncontact]
      integer, allocatable :: NstrR(:), NendR(:)    ! [Ncontact]
      integer, allocatable :: NstrU(:), NendU(:)    ! [Ncontact]
      integer, allocatable :: NstrV(:), NendV(:)    ! [Ncontact]
!
      integer, allocatable :: contact_region(:)     ! [NCdatum]
      integer, allocatable :: on_boundary(:)        ! [NCdatum]
      integer, allocatable :: Idg_cp(:)             ! [NCdatum]
      integer, allocatable :: Jdg_cp(:)             ! [NCdatum]
      integer, allocatable :: Irg_cp(:)             ! [NCdatum]
      integer, allocatable :: Jrg_cp(:)             ! [NCdatum]
!
!-----------------------------------------------------------------------
!  Nested grid connectivity (NGC) structure.
!-----------------------------------------------------------------------
!
!  This structure is used to store all the connectivity information
!  between nested grids.  It will be used extensively when processing
!  contact region points between data donor and data receiver grids.
!  The nested grid contact region information is processed outside of
!  ROMS and read from a NetCDF file for functionality and efficiency.
!
!  In nested grids, the value in the contact region are interpolated
!  from the data donor grid cell using the following conventions at
!  the horizontal location in receiver grid (Irg,Jrg) and donor grid
!  cell (Idg,Jdg):
!
!  suffix 'dg' = donor grid
!         'rg' = receiver grid
!
!            4---------------3 (Idg+1,Jdg+1)   weight(1) = (1-p) * (1-q)
!            |        .      |                 weight(2) =    p  * (1-q)
!            |    1-q .      |                 weight(3) =    p  *  q
!            |        .      |                 weight(4) = (1-p) *  q
!            |        .   p  |
!        Jrg |....... x .....|  Linear interpolation:
!            |  1-p   .      |
!            |        . q    |  V(Irg,Jrg) = weight(1) * F(Idg  ,Jdg  )+
!            |        .      |               weight(2) * F(Idg+1,Jdg  )+
!  (Idg,Jdg) 1---------------2               weight(3) * F(Idg+1,Jdg+1)+
!                    Irg                     weight(4) * F(Idg  ,Jdg+1)
!
!  Notice that if p=0 and q=0 at all contact points, the donor and
!  receiver grids are coincident since weight(1)=1.0 and weight(2:3)=0.
!  Therefore, the above formula is generic for any nested grid
!  configuration.
!
!  If Land/Sea masking, the interpolation weights are rescaled in
!  "mask_weights" during initialization to account masked points in
!  the contact regions. If wetting and drying, the rescaling is done
!  at every time step since the land/sea masking is time dependent.
!
      integer :: Ncontact          ! total number of contact regions
!
      TYPE T_NGC
        logical :: coincident      ! coincident donor/receiver, p=q=0
        logical :: interpolate     ! perform vertical interpolation
        integer :: donor_grid      ! data donor grid number
        integer :: receiver_grid   ! data receiver grid number
        integer :: Npoints         ! number of points in contact region
        integer, pointer :: Irg(:) ! receiver grid, I-contact point
        integer, pointer :: Jrg(:) ! receiver grid, J-contact point
        integer, pointer :: Idg(:) ! donor grid, cell I-left   index
        integer, pointer :: Jdg(:) ! donor grid, cell J-bottom index
# ifdef SOLVE3D
        integer, pointer :: Kdg(:,:)         ! donor grid, cell K-index
# endif
        real(r8), pointer :: Lweight(:,:)          ! linear weights
# ifdef WET_DRY
        real(r8), pointer :: LweightUnmasked(:,:)  ! Unmasked Lweight
# endif
# ifdef QUADRATIC_WEIGHTS
        real(r8), pointer :: Qweight(:,:)          ! quadratic weights
#  ifdef WET_DRY
        real(r8), pointer :: QweightUnmasked(:,:)  ! Unmasked Qweight
#  endif
# endif
# ifdef SOLVE3D
        real(r8), pointer :: Vweight(:,:,:)        ! vertical weights
# endif
      END TYPE T_NGC
!
      TYPE (T_NGC), allocatable :: Rcontact(:)  ! RHO-points, [Ncontact]
      TYPE (T_NGC), allocatable :: Ucontact(:)  ! U-points,   [Ncontact]
      TYPE (T_NGC), allocatable :: Vcontact(:)  ! V-points,   [Ncontact]
!
!-----------------------------------------------------------------------
!  Boundary Contact Points (BCP) structure, allocated as (4,Ncontact).
!  The first dimension is for domain edge (1=iwest,2=isouth, 3=ieast,
!  4=inorth).
!-----------------------------------------------------------------------
!
!  Currently, this structure is only used in refinement grids where the
!  coarser (donor) and finer (receiver) grids have coincident boundaries
!  but with different I- and J-indices.  However, it can be used in the
!  future for composite grids with coincient boundaries.
!
!  The variable "C2Bindex" is used to tell us which contact points in
!  the "Ucontact" and "Vcontact" structure are located at the physical
!  boundary of the relevant nested grid. For example at the boundary
!  edge of a grid with contact region "cr", we can get the mapping
!  between contact point "m" and grid physical boundary edge index
!  "i" or "j" as:
!
!     m = BRY_CONTACT(iwest, cr) % C2Bindex(j)
!     m = BRY_CONTACT(isouth,cr) % C2Bindex(i)
!     m = BRY_CONTACT(ieast, cr) % C2Bindex(j)
!     m = BRY_CONTACT(inorth,cr) % C2Bindex(i)
!
!  This mapping is set during intialization and facilitates efficient
!  processing of nesting contact data.
!
      TYPE T_BCP
        integer :: spv                     ! fill value, unwanted index
        integer :: Ibmin                   ! viable minimum Ib
        integer :: Ibmax                   ! viable maximum Ib
        integer :: Jbmin                   ! viable minimum Jb
        integer :: Jbmax                   ! viable maximum Jb
        integer, pointer :: Ib(:)          ! I-boundary index
        integer, pointer :: Jb(:)          ! J-boundary index

        integer, pointer :: C2Bindex(:)    ! contact to boundary index

# ifdef NESTING_DEBUG
        real(r8), pointer :: Mflux(:)      ! perimeter mass flux
# endif
# ifdef SOLVE3D
        real(r8), pointer :: Tflux(:,:,:)  ! perimeter tracer flux
# endif
      END TYPE T_BCP
!
      TYPE (T_BCP), allocatable :: BRY_CONTACT(:,:)    ! [4,Ncontact]
!
!-----------------------------------------------------------------------
!  Nested Grid Metrics (NGM) structure for contact regions.  Usually,
!  there are contact points outside of the regular (physical) nested
!  grid domain. That is, such contact points are located in the
!  extended (numerical) regions.  These metrics values are computed
!  when designing and generating the application grids.
!
!  It is recommended to build an intermediary fine resolution grid
!  encompassing the study area first and extract/sample all the ROMS
!  application nested grids from it.  This would give a better handle
!  on volume conservation, bathymetry, land/sea masking and other
!  issues.
!-----------------------------------------------------------------------
!
!  These metrics are written the contact points NetCDF file and save
!  separated here. It is very tricky to load these values directly
!  to global grid metrics because of parallelization.
!
       TYPE T_NGM
         real(r8), pointer :: angler(:)    ! angle between XI and EAST
         real(r8), pointer :: dndx(:)      ! d(1/pn)/d(XI)
         real(r8), pointer :: dmde(:)      ! d(1/pm)/d(ETA)
         real(r8), pointer :: f(:)         ! Coriolis parameter
         real(r8), pointer :: h(:)         ! bathymetry
         real(r8), pointer :: rmask(:)     ! land/sea RHO-mask
         real(r8), pointer :: umask(:)     ! land/sea U-mask
         real(r8), pointer :: vmask(:)     ! land/sea V-mask
         real(r8), pointer :: pm(:)        ! XI-coordinate metric
         real(r8), pointer :: pn(:)        ! ETA-coordinate metric
         real(r8), pointer :: Xr(:)        ! X RHO-coordinate (m or deg)
         real(r8), pointer :: Yr(:)        ! Y RHO-coordinate (m or deg)
         real(r8), pointer :: Xu(:)        ! X U-coordinate (m or deg)
         real(r8), pointer :: Yu(:)        ! Y U-coordinate (m or deg)
         real(r8), pointer :: Xv(:)        ! X V-coordinate (m or deg)
         real(r8), pointer :: Yv(:)        ! Y V-coordinate (m or deg)
       END TYPE T_NGM
!
       TYPE (T_NGM), allocatable :: CONTACT_METRIC(:)    ! [Ncontact]
!
!-----------------------------------------------------------------------
!  Composite grids structure.  It contains the donor grid data at the
!  receiver grid contact points. The donor grid data is extracted for
!  the cell containing the contact point: 4 horizontal values to
!  facilitate spatial interpolation.
!-----------------------------------------------------------------------
!
      TYPE T_COMPOSITE
        real(r8), pointer :: bustr(:,:)       ! [4,Npoints]
        real(r8), pointer :: bvstr(:,:)       ! [4,Npoints)

        real(r8), pointer :: ubar(:,:,:)      ! [4,Npoints,2]
        real(r8), pointer :: vbar(:,:,:)      ! [4,Npoints,2]
        real(r8), pointer :: zeta(:,:,:)      ! [4,Npoints,2]

        real(r8), pointer :: rzeta(:,:)       ! [4,Npoints]

# ifdef SOLVE3D
        real(r8), pointer :: DU_avg1(:,:)     ! [4,Npoints]
        real(r8), pointer :: DV_avg1(:,:)     ! [4,Npoints]
        real(r8), pointer :: Zt_avg1(:,:)     ! [4,Npoints]

        real(r8), pointer :: u(:,:,:)         ! [4,k,Npoints]
        real(r8), pointer :: v(:,:,:)         ! [4,k,Npoints]

        real(r8), pointer :: Huon(:,:,:)      ! [4,k,Npoints]
        real(r8), pointer :: Hvom(:,:,:)      ! [4,k,Npoints]

        real(r8), pointer :: t(:,:,:,:)       ! [4,k,Npoints,itrc]
# endif

      END TYPE T_COMPOSITE
!
      TYPE (T_COMPOSITE), allocatable :: COMPOSITE(:)  ! [Ncontact]
!
!-----------------------------------------------------------------------
!  Refinement grids structure: It contains the coarser grid data at the
!  finer grid contact points. The finer grid data is extracted for the
!  cell containing the contact point: 4 horizontal values and 2 time
!  records (t1:t2) to facilitate the space-time interpolation.
!-----------------------------------------------------------------------
!
      TYPE T_REFINED
        real(r8), pointer :: ubar(:,:,:)      ! [4,Npoints,t1:t2]
        real(r8), pointer :: vbar(:,:,:)      ! [4,Npoints,t1:t2]
        real(r8), pointer :: zeta(:,:,:)      ! [4,Npoints,t1:t2]

        real(r8), pointer :: DU_avg2(:,:,:)   ! [4,Npoints,t1:t2]
        real(r8), pointer :: DV_avg2(:,:,:)   ! [4,Npoints,t1:t2]

        real(r8), pointer :: on_u(:)          ! [Npoints]
        real(r8), pointer :: om_v(:)          ! [Npoints]

# ifdef SOLVE3D
        real(r8), pointer :: u(:,:,:,:)       ! [4,k,Npoints,t1:t2]
        real(r8), pointer :: v(:,:,:,:)       ! [4,k,Npoints,t1:t2]

        real(r8), pointer :: t(:,:,:,:,:)     ! [4,k,Npoints,t1:t2,itrc]
# endif
      END TYPE T_REFINED
!
      TYPE (T_REFINED), allocatable :: REFINED(:)      ! [Ncontact]
!
      CONTAINS
!
      SUBROUTINE allocate_nesting
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 2D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Imin, Imax, Jmin, Jmax
      integer :: CCR, cr, dg, ng, rg
      integer :: i, ibry, ic, id, ir, j, jd, jr, m, my_tile
      integer :: ispval

      integer, allocatable :: Ibmin(:,:), Ibmax(:,:)
      integer, allocatable :: Jbmin(:,:), Jbmax(:,:)
!
!-----------------------------------------------------------------------
!  Unpack Boundary Contact Points structure (type T_BCP).
!-----------------------------------------------------------------------
!
!  Allocate boundary connectivity (type T_BCP) structure.
!
      allocate ( BRY_CONTACT(4,Ncontact) )
!
!  Allocate arrays in boundary connectivity structure.
!
      my_tile=-1                           ! for global values
      DO cr=1,Ncontact
        rg=receiver_grid(cr)
        LBi=BOUNDS(rg)%LBi(my_tile)
        UBi=BOUNDS(rg)%UBi(my_tile)
        LBj=BOUNDS(rg)%LBj(my_tile)
        UBj=BOUNDS(rg)%UBj(my_tile)
        DO ibry=1,4
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              allocate ( BRY_CONTACT(ibry,cr) % Ib(LBj:UBj) )
              Dmem(rg)=Dmem(rg)+REAL(UBj-LBj,r8)

              allocate ( BRY_CONTACT(ibry,cr) % Jb(LBj:UBj) )
              Dmem(rg)=Dmem(rg)+REAL(UBj-LBj,r8)

              allocate ( BRY_CONTACT(ibry,cr) % C2Bindex(LBj:UBj) )
              Dmem(rg)=Dmem(rg)+REAL(UBj-LBj,r8)

# ifdef NESTING_DEBUG
              allocate ( BRY_CONTACT(ibry,cr) % Mflux(LBj:UBj) )
              Dmem(rg)=Dmem(rg)+REAL(UBj-LBj,r8)
# endif

# ifdef SOLVE3D
              allocate ( BRY_CONTACT(ibry,cr) % Tflux(LBj:UBj,          &
     &                                                N(rg),NT(rg)) )
              Dmem(rg)=Dmem(rg)+REAL((UBj-LBj)*N(rg)*NT(rg),r8)
# endif
            CASE (isouth, inorth)
              allocate ( BRY_CONTACT(ibry,cr) % Ib(LBi:UBi) )
              Dmem(rg)=Dmem(rg)+REAL(UBi-LBi,r8)

              allocate ( BRY_CONTACT(ibry,cr) % Jb(LBi:UBi) )
              Dmem(rg)=Dmem(rg)+REAL(UBi-LBi,r8)

              allocate ( BRY_CONTACT(ibry,cr) % C2Bindex(LBi:UBi) )
              Dmem(rg)=Dmem(rg)+REAL(UBi-LBi,r8)

# ifdef NESTING_DEBUG
              allocate ( BRY_CONTACT(ibry,cr) % Mflux(LBi:UBi) )
              Dmem(rg)=Dmem(rg)+REAL(UBi-LBi,r8)
# endif

# ifdef SOLVE3D
              allocate ( BRY_CONTACT(ibry,cr) % Tflux(LBi:UBi,          &
     &                                                N(rg),NT(rg)) )
              Dmem(rg)=Dmem(rg)+REAL((UBi-LBi)*N(rg)*NT(rg),r8)
# endif
          END SELECT
        END DO
      END DO
!
!  Initialize boundary connectivity structure: Boundary indices array
!  are initialized to its special value.
!
      ispval=-999
!
      IF (.not.allocated(Ibmin)) THEN
        allocate ( Ibmin(4,Ncontact) )
        DO ng=1,Ngrids
          Dmem(ng)=Dmem(ng)+4.0_r8*REAL(Ncontact,r8)
        END DO
        Ibmin = -ispval
      END IF
      IF (.not.allocated(Ibmax)) THEN
        allocate ( Ibmax(4,Ncontact) )
        DO ng=1,Ngrids
          Dmem(ng)=Dmem(ng)+4.0_r8*REAL(Ncontact,r8)
        END DO
        Ibmax = ispval
      END IF
      IF (.not.allocated(Jbmin)) THEN
        allocate ( Jbmin(4,Ncontact) )
        DO ng=1,Ngrids
          Dmem(ng)=Dmem(ng)+4.0_r8*REAL(Ncontact,r8)
        END DO
        Jbmin = -ispval
      END IF
      IF (.not.allocated(Jbmax)) THEN
        allocate ( Jbmax(4,Ncontact) )
        DO ng=1,Ngrids
          Dmem(ng)=Dmem(ng)+4.0_r8*REAL(Ncontact,r8)
        END DO
        Jbmax = ispval
      END IF
!
      DO cr=1,Ncontact
        DO ibry=1,4
          BRY_CONTACT(ibry,cr) % spv = ispval
          BRY_CONTACT(ibry,cr) % Ib = ispval
          BRY_CONTACT(ibry,cr) % Jb = ispval
          BRY_CONTACT(ibry,cr) % C2Bindex = ispval
# ifdef NESTING_DEBUG
          BRY_CONTACT(ibry,cr) % Mflux = 0.0_r8
# endif
# ifdef SOLVE3D
          BRY_CONTACT(ibry,cr) % Tflux = 0.0_r8
# endif
        END DO
      END DO
!
!  Identify contact points located on the grid boundary.  Notice that
!  the conjugate contact region (CCR) is also processed but it is not
!  yet used. Also, the CCR indices (Ib,Jb) are in refinement overwriten
!  in the m-loop below because several finer grid contact points
!  (RefineScale) are contained in the coarser grid cell. The C2Bindex
!  in this case has the value for the last processed contact point
!  with contact region "cr".
!
      DO m=1,NCdatum
        cr=contact_region(m)
        dg=donor_grid(cr)
        rg=receiver_grid(cr)
        ibry=on_boundary(m)
        DO ic=1,Ncontact
          IF ((dg.eq.receiver_grid(ic)).and.                            &
     &        (rg.eq.donor_grid(ic))) THEN
            CCR=ic                            ! conjugate contact region
            EXIT
          END IF
        END DO
        IF ((ibry.eq.iwest ).or.(ibry.eq.ieast )) THEN
          ir=Irg_cp(m)
          jr=Jrg_cp(m)
          Ibmin(ibry,cr )=MIN(ir,Ibmin(ibry,cr ))
          Ibmax(ibry,cr )=MAX(ir,Ibmax(ibry,cr ))
          Jbmin(ibry,cr )=MIN(jr,Jbmin(ibry,cr ))
          Jbmax(ibry,cr )=MAX(jr,Jbmax(ibry,cr ))
          BRY_CONTACT(ibry,cr ) % Ib(jr) = ir
          BRY_CONTACT(ibry,cr ) % Jb(jr) = jr
          BRY_CONTACT(ibry,cr ) % C2Bindex(jr) = m-NstrU(cr)+1
!
          id=Idg_cp(m)
          jd=Jdg_cp(m)
          Ibmin(ibry,CCR)=MIN(id,Ibmin(ibry,CCR))
          Ibmax(ibry,CCR)=MAX(id,Ibmax(ibry,CCR))
          Jbmin(ibry,CCR)=MIN(jd,Jbmin(ibry,CCR))
          Jbmax(ibry,CCR)=MAX(jd,Jbmax(ibry,CCR))
          BRY_CONTACT(ibry,CCR) % Ib(jd) = id
          BRY_CONTACT(ibry,CCR) % Jb(jd) = jd
          BRY_CONTACT(ibry,CCR) % C2Bindex(jd) = m-NstrU(cr)+1  ! same
        ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
          ir=Irg_cp(m)
          jr=Jrg_cp(m)
          Ibmin(ibry,cr )=MIN(ir,Ibmin(ibry,cr))
          Ibmax(ibry,cr )=MAX(ir,Ibmax(ibry,cr))
          Jbmin(ibry,cr )=MIN(jr,Jbmin(ibry,cr))
          Jbmax(ibry,cr )=MAX(jr,Jbmax(ibry,cr))
          BRY_CONTACT(ibry,cr ) % Ib(ir) = ir
          BRY_CONTACT(ibry,cr ) % Jb(ir) = jr
          BRY_CONTACT(ibry,cr ) % C2Bindex(ir) = m-NstrV(cr)+1
!
          id=Idg_cp(m)
          jd=Jdg_cp(m)
          Ibmin(ibry,CCR)=MIN(id,Ibmin(ibry,CCR))
          Ibmax(ibry,CCR)=MAX(id,Ibmax(ibry,CCR))
          Jbmin(ibry,CCR)=MIN(jd,Jbmin(ibry,CCR))
          Jbmax(ibry,CCR)=MAX(jd,Jbmax(ibry,CCR))
          BRY_CONTACT(ibry,CCR) % Ib(id) = id
          BRY_CONTACT(ibry,CCR) % Jb(id) = jd
          BRY_CONTACT(ibry,CCR) % C2Bindex(id) = m-NstrV(cr)+1  ! same
        END IF
      END DO
!
!  Set minimum and maximum indices to process at each boundary.
!
      DO cr=1,Ncontact
        DO ibry=1,4
          IF (ABS(Ibmin(ibry,cr)).eq.ABS(ispval)) THEN
            BRY_CONTACT(ibry,cr) % Ibmin = ispval
          ELSE
            BRY_CONTACT(ibry,cr) % Ibmin = Ibmin(ibry,cr)
          END IF

          IF (ABS(Ibmax(ibry,cr)).eq.ABS(ispval)) THEN
            BRY_CONTACT(ibry,cr) % Ibmax = ispval
          ELSE
            BRY_CONTACT(ibry,cr) % Ibmax = Ibmax(ibry,cr)
          END IF

          IF (ABS(Jbmin(ibry,cr)).eq.ABS(ispval)) THEN
            BRY_CONTACT(ibry,cr) % Jbmin = ispval
          ELSE
            BRY_CONTACT(ibry,cr) % Jbmin = Jbmin(ibry,cr)
          END IF

          IF (ABS(Jbmax(ibry,cr)).eq.ABS(ispval)) THEN
            BRY_CONTACT(ibry,cr) % Jbmax = ispval
          ELSE
            BRY_CONTACT(ibry,cr) % Jbmax = Jbmax(ibry,cr)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Deactivate boundary condition switches if contact point lay on the
!  physical nested grid boundary.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
        rg=receiver_grid(cr)
        IF (RefinedGrid(rg)) THEN
          IF (RefineScale(rg).gt.0) THEN
            LBC_apply(rg) % west  = .FALSE.    ! This is a refinement
            LBC_apply(rg) % south = .FALSE.    ! grid, so we need to
            LBC_apply(rg) % east  = .FALSE.    ! avoid applying lateral
            LBC_apply(rg) % north = .FALSE.    ! boundary conditions
          END IF
        ELSE
          DO ibry=1,4
            Imin=BRY_CONTACT(ibry,cr) % Ibmin  ! Deactivate full or
            Imax=BRY_CONTACT(ibry,cr) % Ibmax  ! partial lateral
            Jmin=BRY_CONTACT(ibry,cr) % Jbmin  ! boundary conditions
            Jmax=BRY_CONTACT(ibry,cr) % Jbmax
            SELECT CASE (ibry)
              CASE (iwest)
                IF ((Jmin.ne.ispval).and.(Jmax.ne.ispval)) THEN
                  DO j=Jmin,Jmax
                    LBC_apply(rg) % west (j) = .FALSE.
                  END DO
                END IF
              CASE (isouth)
                IF ((Imin.ne.ispval).and.(Imax.ne.ispval)) THEN
                  DO i=Imin,Imax
                    LBC_apply(rg) % south(i) = .FALSE.
                  END DO
                END IF
              CASE (ieast)
                IF ((Jmin.ne.ispval).and.(Jmax.ne.ispval)) THEN
                  DO j=Jmin,Jmax
                    LBC_apply(rg) % east (j) = .FALSE.
                  END DO
                END IF
              CASE (inorth)
                IF ((Imin.ne.ispval).and.(Imax.ne.ispval)) THEN
                  DO i=Imin,Imax
                    LBC_apply(rg) % north(i) = .FALSE.
                  END DO
                END IF
            END SELECT
          END DO
        END IF
      END DO
!
      RETURN
      END SUBROUTINE allocate_nesting
!
      SUBROUTINE deallocate_nesting
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 2D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================

# ifdef SUBOBJECT_DEALLOCATION
!
      USE destroy_mod, ONLY : destroy
# endif
!
!  Local variable declarations.
!
      integer ::  cr, ibry
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", deallocate_nesting"

# ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_BCP boundary
!  connectivity structure separately.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
        DO ibry=1,4
          IF (.not.destroy(ng, BRY_CONTACT(ibry,cr)%Ib, MyFile,         &
     &                     __LINE__, 'BRY_CONTACT(ibry,cr)%Ib')) RETURN

          IF (.not.destroy(ng, BRY_CONTACT(ibry,cr)%Jb, MyFile,         &
     &                     __LINE__, 'BRY_CONTACT(ibry,cr)%Jb')) RETURN

          IF (.not.destroy(ng, BRY_CONTACT(ibry,cr)%C2Bindex, MyFile,   &
     &                __LINE__, 'BRY_CONTACT(ibry,cr)%C2Bindex')) RETURN

#  ifdef NESTING_DEBUG
          IF (.not.destroy(ng, BRY_CONTACT(ibry,cr)%Mflux, MyFile,      &
     &                __LINE__, 'BRY_CONTACT(ibry,cr)%Mflux')) RETURN
#  endif

#  ifdef SOLVE3D
          IF (.not.destroy(ng, BRY_CONTACT(ibry,cr)%Tflux, MyFile,      &
     &                __LINE__, 'BRY_CONTACT(ibry,cr)%Tflux')) RETURN
#  endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_NGC grid connectivity
!  structure separately.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
        IF (.not.destroy(ng, Rcontact(cr)%Irg, MyFile,                  &
     &                   __LINE__, 'Rcontact(cr)%Irg')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%Irg, MyFile,                  &
     &                   __LINE__, 'Ucontact(cr)%Irg')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%Irg, MyFile,                  &
     &                   __LINE__, 'Vcontact(cr)%Irg')) RETURN
!
        IF (.not.destroy(ng, Rcontact(cr)%Jrg, MyFile,                  &
     &                   __LINE__, 'Rcontact(cr)%Jrg')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%Jrg, MyFile,                  &
     &                   __LINE__, 'Ucontact(cr)%Jrg')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%Jrg, MyFile,                  &
     &                   __LINE__, 'Vcontact(cr)%Jrg')) RETURN
!
        IF (.not.destroy(ng, Rcontact(cr)%Idg, MyFile,                  &
     &                   __LINE__, 'Rcontact(cr)%Idg')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%Idg, MyFile,                  &
     &                   __LINE__, 'Ucontact(cr)%Idg')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%Idg, MyFile,                  &
     &                   __LINE__, 'Vcontact(cr)%Idg')) RETURN
!
        IF (.not.destroy(ng, Rcontact(cr)%Jdg, MyFile,                  &
     &                   __LINE__, 'Rcontact(cr)%Jdg')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%Jdg, MyFile,                  &
     &                   __LINE__, 'Ucontact(cr)%Jdg')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%Jdg, MyFile,                  &
     &                   __LINE__, 'Vcontact(cr)%Jdg')) RETURN

#  ifdef SOLVE3D
!
        IF (.not.destroy(ng, Rcontact(cr)%Kdg, MyFile,                  &
     &                   __LINE__, 'Rcontact(cr)%Kdg')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%Kdg, MyFile,                  &
     &                   __LINE__, 'Ucontact(cr)%Kdg')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%Kdg, MyFile,                  &
     &                   __LINE__, 'Vcontact(cr)%Kdg')) RETURN
#  endif
!
        IF (.not.destroy(ng, Rcontact(cr)%Lweight, MyFile,              &
     &                   __LINE__, 'Rcontact(cr)%Lweight')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%Lweight, MyFile,              &
     &                   __LINE__, 'Ucontact(cr)%Lweight')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%Lweight, MyFile,              &
     &                   __LINE__, 'Vcontact(cr)%Lweight')) RETURN
!
#  ifdef WET_DRY
        IF (.not.destroy(ng, Rcontact(cr)%LweightUnmasked, MyFile,      &
     &              __LINE__, 'Rcontact(cr)%LweightUnmasked')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%LweightUnmasked, MyFile,      &
     &              __LINE__, 'Ucontact(cr)%LweightUnmasked')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%LweightUnmasked, MyFile,      &
     &              __LINE__, 'Vcontact(cr)%LweightUnmasked')) RETURN
#  endif

#  ifdef QUADRATIC_WEIGHTS
!
        IF (.not.destroy(ng, Rcontact(cr)%Qweight, MyFile,              &
     &                   __LINE__, 'Rcontact(cr)%Qweight')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%Qweight, MyFile,              &
     &                   __LINE__, 'Ucontact(cr)%Qweight')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%Qweight, MyFile,              &
     &                   __LINE__, 'Vcontact(cr)%Qweight')) RETURN

#   ifdef WET_DRY
!
        IF (.not.destroy(ng, Rcontact(cr)%QweightUnmasked, MyFile,      &
     &              __LINE__, 'Rcontact(cr)%QweightUnmasked')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%QweightUnmasked, MyFile,      &
     &              __LINE__, 'Ucontact(cr)%QweightUnmasked')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%QweightUnmasked, MyFile,      &
     &              __LINE__, 'Vcontact(cr)%QweightUnmasked')) RETURN
#   endif
#  endif

#  ifdef SOLVE3D
!
        IF (.not.destroy(ng, Rcontact(cr)%Vweight, MyFile,              &
     &                   __LINE__, 'Rcontact(cr)%Vweight')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%Vweight, MyFile,              &
     &                   __LINE__, 'Ucontact(cr)%Vweight')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%Vweight, MyFile,              &
     &                   __LINE__, 'Vcontact(cr)%Vweight')) RETURN

#   if defined TANGENT || defined TL_IOMS
!
        IF (.not.destroy(ng, Rcontact(cr)%tl_Vweight, MyFile,           &
     &                   __LINE__, 'Rcontact(cr)%tl_Vweight')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%tl_Vweight, MyFile,           &
     &                   __LINE__, 'Ucontact(cr)%tl_Vweight')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%tl_Vweight, MyFile,           &
     &                   __LINE__, 'Vcontact(cr)%tl_Vweight')) RETURN
#   endif

#   ifdef ADJOINT
!
        IF (.not.destroy(ng, Rcontact(cr)%ad_Vweight, MyFile,           &
     &                   __LINE__, 'Rcontact(cr)%ad_Vweight')) RETURN
        IF (.not.destroy(ng, Ucontact(cr)%ad_Vweight, MyFile,           &
     &                   __LINE__, 'Ucontact(cr)%ad_Vweight')) RETURN
        IF (.not.destroy(ng, Vcontact(cr)%ad_Vweight, MyFile,           &
     &                   __LINE__, 'Vcontact(cr)%ad_Vweight')) RETURN
#   endif
#  endif
      END DO
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_NGC contact region
!  metrics structure separately.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
        IF (.not.destroy(ng, CONTACT_METRIC(cr)%angler, MyFile,         &
     &                   __LINE__, 'CONTACT_METRIC(cr)%angler')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%dndx, MyFile,           &
     &                   __LINE__, 'CONTACT_METRIC(cr)%dndx')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%dmde, MyFile,           &
     &                   __LINE__, 'CONTACT_METRIC(cr)%dmde')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%f, MyFile,              &
     &                   __LINE__, 'CONTACT_METRIC(cr)%f')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%h, MyFile,              &
     &                   __LINE__, 'CONTACT_METRIC(cr)%h')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%rmask, MyFile,          &
     &                   __LINE__, 'CONTACT_METRIC(cr)%rmask')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%umask, MyFile,          &
     &                   __LINE__, 'CONTACT_METRIC(cr)%umask')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%vmask, MyFile,          &
     &                   __LINE__, 'CONTACT_METRIC(cr)%vmask')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%pm, MyFile,             &
     &                   __LINE__, 'CONTACT_METRIC(cr)%pm')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%pn, MyFile,             &
     &                   __LINE__, 'CONTACT_METRIC(cr)%pn')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%Xr, MyFile,             &
     &                   __LINE__, 'CONTACT_METRIC(cr)%Xr')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%Yr, MyFile,             &
     &                   __LINE__, 'CONTACT_METRIC(cr)%Yr')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%Xu, MyFile,             &
     &                   __LINE__, 'CONTACT_METRIC(cr)%Xu')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%Yu, MyFile,             &
     &                   __LINE__, 'CONTACT_METRIC(cr)%Yu')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%Xv, MyFile,             &
     &                   __LINE__, 'CONTACT_METRIC(cr)%Xv')) RETURN

        IF (.not.destroy(ng, CONTACT_METRIC(cr)%Yv, MyFile,             &
     &                   __LINE__, 'CONTACT_METRIC(cr)%Yv')) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_COMPOSITE for
!  composite grids contact region structure separately.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
        IF (.not.destroy(ng, COMPOSITE(cr)%bustr, MyFile,               &
     &                   __LINE__, 'COMPOSITE(cr)%bustr')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%bvstr, MyFile,               &
     &                   __LINE__, 'COMPOSITE(cr)%bvstr')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%ubar, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%ubar')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%vbar, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%vbar')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%zeta, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%zeta')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%rzeta, MyFile,               &
     &                   __LINE__, 'COMPOSITE(cr)%rzeta')) RETURN

#   if defined TANGENT || defined TL_IOMS
        IF (.not.destroy(ng, COMPOSITE(cr)%tl_bustr, MyFile,            &
     &                   __LINE__, 'COMPOSITE(cr)%tl_bustr')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%tl_bvstr, MyFile,            &
     &                   __LINE__, 'COMPOSITE(cr)%tl_bvstr')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%tl_ubar, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%tl_ubar')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%tl_vbar, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%tl_vbar')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%tl_zeta, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%tl_zeta')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%tl_rzeta, MyFile,            &
     &                   __LINE__, 'COMPOSITE(cr)%tl_rzeta')) RETURN
#   endif

#   ifdef ADJOINT
        IF (.not.destroy(ng, COMPOSITE(cr)%ad_bustr, MyFile,            &
     &                   __LINE__, 'COMPOSITE(cr)%ad_bustr')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%ad_bvstr, MyFile,            &
     &                   __LINE__, 'COMPOSITE(cr)%ad_bvstr')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%ad_ubar, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%ad_ubar')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%ad_vbar, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%ad_vbar')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%ad_zeta, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%ad_zeta')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%ad_rzeta, MyFile,            &
     &                   __LINE__, 'COMPOSITE(cr)%ad_rzeta')) RETURN
#   endif

#  ifdef SOLVE3D
        IF (.not.destroy(ng, COMPOSITE(cr)%DU_avg1, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%DU_avg1')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%DV_avg1, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%DV_avg1')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%Zt_avg1, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%Zt_avg1')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%u, MyFile,                   &
     &                   __LINE__, 'COMPOSITE(cr)%u')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%v, MyFile,                   &
     &                   __LINE__, 'COMPOSITE(cr)%v')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%Huon, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%Huon')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%Hvom, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%Hvom')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%t, MyFile,                   &
     &                   __LINE__, 'COMPOSITE(cr)%t')) RETURN

#   if defined TANGENT || defined TL_IOMS
        IF (.not.destroy(ng, COMPOSITE(cr)%tl_DU_avg1, MyFile,          &
     &                   __LINE__, 'COMPOSITE(cr)%tl_DU_avg1')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%tl_DV_avg1, MyFile,          &
     &                   __LINE__, 'COMPOSITE(cr)%tl_DV_avg1')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%tl_Zt_avg1, MyFile,          &
     &                   __LINE__, 'COMPOSITE(cr)%tl_Zt_avg1')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%tl_u, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%tl_u')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%tl_v, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%tl_v')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%tl_Huon, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%tl_Huon')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%tl_Hvom, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%tl_Hvom')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%tl_t, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%tl_t')) RETURN
#   endif

#   ifdef ADJOINT
        IF (.not.destroy(ng, COMPOSITE(cr)%ad_DU_avg1, MyFile,          &
     &                   __LINE__, 'COMPOSITE(cr)%ad_DU_avg1')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%ad_DV_avg1, MyFile,          &
     &                   __LINE__, 'COMPOSITE(cr)%ad_DV_avg1')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%ad_Zt_avg1, MyFile,          &
     &                   __LINE__, 'COMPOSITE(cr)%ad_Zt_avg1')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%ad_u, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%ad_u')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%ad_v, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%ad_v')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%ad_Huon, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%ad_Huon')) RETURN
        IF (.not.destroy(ng, COMPOSITE(cr)%ad_Hvom, MyFile,             &
     &                   __LINE__, 'COMPOSITE(cr)%ad_Hvom')) RETURN

        IF (.not.destroy(ng, COMPOSITE(cr)%ad_t, MyFile,                &
     &                   __LINE__, 'COMPOSITE(cr)%ad_t')) RETURN
#   endif
#  endif
      END DO
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type T_REFINED for
!  refinement grids contact region structure separately.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
        IF (.not.destroy(ng, REFINED(cr)%ubar, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%ubar')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%vbar, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%vbar')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%zeta, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%zeta')) RETURN

        IF (.not.destroy(ng, REFINED(cr)%DU_avg2, MyFile,               &
     &                   __LINE__, 'REFINED(cr)%DU_avg2')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%DV_avg2, MyFile,               &
     &                   __LINE__, 'REFINED(cr)%DV_avg2')) RETURN

        IF (.not.destroy(ng, REFINED(cr)%on_u, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%on_u')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%om_v, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%om_v')) RETURN

#  if defined TANGENT || defined TL_IOMS
        IF (.not.destroy(ng, REFINED(cr)%tl_ubar, MyFile,               &
     &                   __LINE__, 'REFINED(cr)%tl_ubar')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%tl_vbar, MyFile,               &
     &                   __LINE__, 'REFINED(cr)%tl_vbar')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%tl_zeta, MyFile,               &
     &                   __LINE__, 'REFINED(cr)%tl_zeta')) RETURN

        IF (.not.destroy(ng, REFINED(cr)%tl_DU_avg2, MyFile,            &
     &                   __LINE__, 'REFINED(cr)%tl_DU_avg2')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%tl_DV_avg2, MyFile,            &
     &                   __LINE__, 'REFINED(cr)%tl_DV_avg2')) RETURN
#  endif

#  ifdef ADJOINT
        IF (.not.destroy(ng, REFINED(cr)%ad_ubar, MyFile,               &
     &                   __LINE__, 'REFINED(cr)%ad_ubar')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%ad_vbar, MyFile,               &
     &                   __LINE__, 'REFINED(cr)%ad_vbar')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%ad_zeta, MyFile,               &
     &                   __LINE__, 'REFINED(cr)%ad_zeta')) RETURN

        IF (.not.destroy(ng, REFINED(cr)%ad_DU_avg2, MyFile,            &
     &                   __LINE__, 'REFINED(cr)%ad_DU_avg2')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%ad_DV_avg2, MyFile,            &
     &                   __LINE__, 'REFINED(cr)%ad_DV_avg2')) RETURN
#  endif

#  ifdef SOLVE3D
        IF (.not.destroy(ng, REFINED(cr)%u, MyFile,                     &
     &                   __LINE__, 'REFINED(cr)%u')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%v, MyFile,                     &
     &                   __LINE__, 'REFINED(cr)%v')) RETURN

        IF (.not.destroy(ng, REFINED(cr)%t, MyFile,                     &
     &                   __LINE__, 'REFINED(cr)%t')) RETURN

#   if defined TANGENT || defined TL_IOMS
        IF (.not.destroy(ng, REFINED(cr)%tl_u, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%tl_u')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%tl_v, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%tl_v')) RETURN

        IF (.not.destroy(ng, REFINED(cr)%tl_t, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%tl_t')) RETURN
#   endif

#   ifdef ADJOINT
        IF (.not.destroy(ng, REFINED(cr)%ad_u, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%ad_u')) RETURN
        IF (.not.destroy(ng, REFINED(cr)%ad_v, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%ad_v')) RETURN

        IF (.not.destroy(ng, REFINED(cr)%ad_t, MyFile,                  &
     &                   __LINE__, 'REFINED(cr)%ad_t')) RETURN
#   endif
#  endif
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type structures:
!-----------------------------------------------------------------------
!
!  Boundary connectivity.

      IF (allocated(BRY_CONTACT)) deallocate ( BRY_CONTACT )
!
!  Grid connectivity.
!
      IF (allocated(Rcontact)) deallocate ( Rcontact )
      IF (allocated(Ucontact)) deallocate ( Ucontact )
      IF (allocated(Vcontact)) deallocate ( Vcontact )
!
!  Contact region metrics.
!
      IF (allocated(CONTACT_METRIC)) deallocate ( CONTACT_METRIC )
!
!  Composite grid contact regions.
!
      IF (allocated(COMPOSITE)) deallocate ( COMPOSITE )
!
!  Deallocate refinement grids contact regions structure.
!
      IF (allocated(REFINED)) deallocate ( REFINED )
!
!-----------------------------------------------------------------------
!  Deallocate other variables in module.
!-----------------------------------------------------------------------
!
      IF (allocated(ContactRegion)) THEN
        deallocate ( CoarserDonor )
      END IF

      IF (allocated(FinerDonor)) THEN
        deallocate ( FinerDonor )
      END IF

      IF (allocated(DonorToFiner)) THEN
        deallocate ( DonorToFiner )
      END IF

      IF (allocated(RefineSteps)) THEN
        deallocate ( RefineSteps )
      END IF

      IF (allocated(RefineStepsCounter)) THEN
        deallocate ( RefineStepsCounter )
      END IF

      IF (allocated(TwoWayInterval)) THEN
        deallocate ( TwoWayInterval )
      END IF

      IF (allocated(Telescoping)) THEN
        deallocate ( Telescoping )
      END IF

      IF (allocated(RollingIndex)) THEN
        deallocate ( RollingIndex )
      END IF

      IF (allocated(RollingTime)) THEN
        deallocate ( RollingTime )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_nesting
!
      SUBROUTINE initialize_nesting
!
!=======================================================================
!                                                                      !
!  This routine initializes time varying nesting structures.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: cr

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Initialize time-varying contact regions structures.  They are used
!  to process values from contact regions to global kernel arrays and
!  vice versa.
!-----------------------------------------------------------------------
!
!  Composite grids contact region structure.
!
      IF (ANY(CompositeGrid)) THEN
        DO cr=1,Ncontact
          COMPOSITE(cr) % bustr = IniVal
          COMPOSITE(cr) % bvstr = IniVal

          COMPOSITE(cr) % ubar = IniVal
          COMPOSITE(cr) % vbar = IniVal
          COMPOSITE(cr) % zeta = IniVal

          COMPOSITE(cr) % rzeta = IniVal

# ifdef SOLVE3D
          COMPOSITE(cr) % DU_avg1 = IniVal
          COMPOSITE(cr) % DV_avg1 = IniVal
          COMPOSITE(cr) % Zt_avg1 = IniVal

          COMPOSITE(cr) % u = IniVal
          COMPOSITE(cr) % v = IniVal

          COMPOSITE(cr) % Huon = IniVal
          COMPOSITE(cr) % Hvom = IniVal

          COMPOSITE(cr) % t = IniVal
# endif
        END DO
      END IF
!
!  Refinement grids contact region structure.
!
      IF (ANY(RefinedGrid(:))) THEN
        DO cr=1,Ncontact
          REFINED(cr) % ubar = IniVal
          REFINED(cr) % vbar = IniVal
          REFINED(cr) % zeta = IniVal

          REFINED(cr) % DU_avg2 = IniVal
          REFINED(cr) % DV_avg2 = IniVal

          REFINED(cr) % on_u = IniVal
          REFINED(cr) % om_v = IniVal

# ifdef SOLVE3D
          REFINED(cr) % u = IniVal
          REFINED(cr) % v = IniVal

          REFINED(cr) % t = IniVal
# endif
        END DO
      END IF

      RETURN
      END SUBROUTINE initialize_nesting
#endif
      END MODULE mod_nesting
