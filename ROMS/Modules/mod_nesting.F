#include "cppdefs.h"
      MODULE mod_nesting

#ifdef NESTING
!
!svn $Id: mod_nesting.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
!  Counter for two-way exchange between fine and coarse grids.
!
      integer, allocatable :: TwoWayCounter(:)      ! [Ngrids]
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
      real(r8), allocatable :: RollingTime(:,:)     ! [Ncontact]
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

        real(r8), pointer :: Mflux(:)      ! perimeter mass flux
# ifdef SOLVE3D
        real(r8), pointer :: Tflux(:,:,:)  ! perimeter tracer flux
# endif
      END TYPE T_BCP

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
              allocate ( BRY_CONTACT(ibry,cr) % Jb(LBj:UBj) )
              allocate ( BRY_CONTACT(ibry,cr) % C2Bindex(LBj:UBj) )
              allocate ( BRY_CONTACT(ibry,cr) % Mflux(LBj:UBj) )
# ifdef SOLVE3D
              allocate ( BRY_CONTACT(ibry,cr) % Tflux(LBj:UBj,          &
     &                                                N(rg),NT(rg)) )
# endif
            CASE (isouth, inorth)
              allocate ( BRY_CONTACT(ibry,cr) % Ib(LBi:UBi) )
              allocate ( BRY_CONTACT(ibry,cr) % Jb(LBi:UBi) )
              allocate ( BRY_CONTACT(ibry,cr) % C2Bindex(LBi:UBi) )
              allocate ( BRY_CONTACT(ibry,cr) % Mflux(LBi:UBi) )
# ifdef SOLVE3D
              allocate ( BRY_CONTACT(ibry,cr) % Tflux(LBi:UBi,          &
     &                                                N(rg),NT(rg)) )
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
        Ibmin = -ispval
      END IF
      IF (.not.allocated(Ibmax)) THEN
        allocate ( Ibmax(4,Ncontact) )
        Ibmax = ispval
      END IF
      IF (.not.allocated(Jbmin)) THEN
        allocate ( Jbmin(4,Ncontact) )
        Jbmin = -ispval
      END IF
      IF (.not.allocated(Jbmax)) THEN
        allocate ( Jbmax(4,Ncontact) )
        Jbmax = ispval
      END IF
!
      DO cr=1,Ncontact
        DO ibry=1,4
          BRY_CONTACT(ibry,cr) % spv = ispval
          BRY_CONTACT(ibry,cr) % Ib = ispval
          BRY_CONTACT(ibry,cr) % Jb = ispval
          BRY_CONTACT(ibry,cr) % C2Bindex = ispval
          BRY_CONTACT(ibry,cr) % Mflux = 0.0_r8
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

      RETURN
      END SUBROUTINE allocate_nesting
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
