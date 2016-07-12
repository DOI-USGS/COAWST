#include "cppdefs.h"
      MODULE nesting_mod

#ifdef NESTING
!
!svn $Id: nesting.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                             John C. Warner   !
!=======================================================================
!                                                                      !
!  This module contains several routines  to process the connectivity  !
!  between nested grids. It process the contact region points between  !
!  data donor and data receiver grids.                                 !
!                                                                      !
!  The locations of the linear interpolation weights in the donor      !
!  grid with respect the receiver grid contact region at contact       !
!  point x(Irg,Jrg,Krg) are:                                           !
!                                                                      !
!                       8___________7   (Idg+1,Jdg+1,Kdg)              !
!                      /.          /|                                  !
!                     / .         / |                                  !
!  (Idg,Jdg+1,Kdg)  5/___________/6 |                                  !
!                    |  .        |  |                                  !
!                    |  .   x    |  |                                  !
!                    | 4.........|..|3  (Idg+1,Jdg+1,Kdg-1)            !
!                    | .         |  /                                  !
!                    |.          | /                                   !
!                    |___________|/                                    !
!  (Idg,Jdg,Kdg-1)   1           2                                     !
!                                                                      !
!                                        Suffix:   dg = donor grid     !
!                                                  rg = receiver grid  !
!                                                                      !
!  Routines:                                                           !
!  ========                                                            !
!                                                                      !
!  nesting          Public interface to time-stepping kernel           !
!                                                                      !
!  get_composite    Composite  grid, extract contact points donor data !
!  get_metrics      Extract grid spacing metrics, on_u and om_v.       !
!  get_refine       Refinement grid, extract contact points donor data !
!  put_composite    Composite  grid, fill contact points (interpolate) !
!  put_refine       Refinement grid, fill contact points (interpolate) !
!                                                                      !
!  bry_fluxes       Extracts horizontat advective fluxes the contact   !
!                     boundary of donor and receiver grids.            !
!  check_massflux   If refinement, check mass fluxes between coarse    !
!                     and fine grids for volume conservation. It is    !
!                     use only for debugging and diagnostics.          !
!  correct_tracer   Correct coarse grid tracer at the refinement grid  !
!                     boundary with the refined accumulated fluxes     !
!  do_twoway        Logical function to determine at which time-step   !
!                     the two-way exchange takes place.                !
!  fill_contact     Used to Fill grid metrics at contact points.       !
!  fine2coarse      Replace coarse grid state variables with the       !
!                     averaged fine grid values (two-way nesting)      !
!                                                                      !
!  get_contact2d    Get 2D field donor grid cell holding contact point !
!  get_contact3d    Get 3D field donor grid cell holding contact point !
!  get_persisted2d  Get 2D field persisted values on contact points    !
!  mask_hweights    Scale horizontal interpolation weights with masking!
!  put_contact2d    Set 2D field contact points, spatial interpolation !
!  put_contact3d    Set 3D field contact points, spatial interpolation !
!                                                                      !
!  put_refine2d     Interpolate (space-time) 2D state variables        !
!  put_refine3d     Interpolate (space-time) 3D state variables        !
!                                                                      !
!  z_weights        Set donor grid vertical indices (cell holding      !
!                     contact point) and vertical interpolation        !
!                     weights                                          !
!                                                                      !
!  WARNINGS:                                                           !
!  ========                                                            !
!                                                                      !
!  All the routines contained in this module are inside of a parallel  !
!  region, except the main driver routine "nesting",  which is called  !
!  serially several times from main2d or main3d to perform  different  !
!  tasks.  Notice that the calls to private "get_***"  and  "put_***"  !
!  routines need to be in separated  parallel loops because of serial  !
!  with partitions and  shared-memory rules.  Furthermore,  the donor  !
!  and receiver grids may have different tile partitions. There is no  !
!  I/O management inside the nesting routines.                         !
!                                                                      !
!  The connectivity between donor and receiver grids can be  complex.  !
!  The horizontal mapping between grids is static and done outside of  !
!  ROMS.  Only the time-dependent  vertical interpolation weights are  !
!  computed here.  The contact region points  I- and  J-cell  indices  !
!  between donor and receiver grids, and the horizontal interpolation  !
!  weights are read from  the input nesting connectivity NetCDF file.  !
!  It makes the nesting efficient and greatly simplifies parallelism.  !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PUBLIC  :: nesting
      PUBLIC  :: bry_fluxes
# ifndef ONE_WAY
      PUBLIC  :: do_twoway
# endif
      PUBLIC  :: fill_contact
      PRIVATE :: check_massflux
# ifdef SOLVE3D
      PRIVATE :: correct_tracer
      PRIVATE :: correct_tracer_tile
# endif
      PRIVATE :: fine2coarse
      PRIVATE :: fine2coarse2d
# ifdef SOLVE3D
      PRIVATE :: fine2coarse3d
# endif
      PRIVATE :: get_contact2d
# ifdef SOLVE3D
      PRIVATE :: get_contact3d
# endif
      PRIVATE :: get_composite
      PRIVATE :: get_metrics
      PRIVATE :: get_persisted2d
      PRIVATE :: get_refine
# if defined MASKING || defined WET_DRY
      PRIVATE :: mask_hweights
# endif
      PRIVATE :: put_composite
      PRIVATE :: put_refine
      PRIVATE :: put_refine2d
# ifdef SOLVE3D
      PRIVATE :: put_refine3d
      PRIVATE :: z_weights
# endif
!
      CONTAINS
!
      SUBROUTINE nesting (ng, model, isection)
!
!=======================================================================
!                                                                      !
!  This routine process the contact region points between composite    !
!  grids.  In composite grids, it is possible to have more than one    !
!  contact region.                                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Data receiver grid number (integer)                   !
!     model      Calling model identifier (integer)                    !
!     isection   Governing equations time-stepping section in          !
!                  main2d or main3d indicating which state             !
!                  variables to process (integer)                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_nesting
      USE mod_scalars

# ifdef SOLVE3D
!
      USE set_depth_mod, ONLY : set_depth
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, isection
!
!  Local variable declarations.
!
      logical :: LputFsur
      integer :: subs, tile, thread
# ifdef SOLVE3D
      integer :: ngc
!
!-----------------------------------------------------------------------
!  Process vertical indices and interpolation weigths associated with
!  depth.
!-----------------------------------------------------------------------
!
      IF (isection.eq.nzwgt) THEN
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL z_weights (ng, model, tile)
        END DO
!$OMP BARRIER
        RETURN
      END IF
# endif

# if defined MASKING || defined WET_DRY
!
!-----------------------------------------------------------------------
!  If Land/Sea masking, scale horizontal interpolation weights to
!  account for land contact points. If wetting and drying, the scaling
!  is done at every time-step because masking is time dependent.
!-----------------------------------------------------------------------
!
      IF (isection.eq.nmask) THEN
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL mask_hweights (ng, model, tile)
        END DO
!$OMP BARRIER
        RETURN
      END IF
# endif
!
!-----------------------------------------------------------------------
!  If refinement grid, process contact points.
!-----------------------------------------------------------------------
!
      IF (RefinedGrid(ng)) THEN
!
!  Extract grid spacing metrics (on_u and om_v) and load then to
!  REFINE(:) structure.  These metrics are needed to impose mass
!  flux at the finer grid physical boundaries.  It need to be done
!  separately because parallelism partions between all nested grid.
!
        IF (isection.eq.ndxdy) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL get_metrics (ng, model, tile)
          END DO
!$OMP BARRIER
!
!  Extract and store donor grid data at contact points.
!
        ELSE IF (isection.eq.ngetD) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL get_refine (ng, model, tile)
          END DO
!$OMP BARRIER
!
!  Fill refinement grid contact points variables by interpolating
!  (space, time) from extracted donor grid data.  The free-surface
!  needs to be processed first and in a separate parallel region
!  because of shared-memory applications.
!
        ELSE IF (isection.eq.nputD) THEN
          LputFsur=.TRUE.
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL put_refine (ng, model, tile, LputFsur)
          END DO
!$OMP BARRIER
!
          LputFsur=.FALSE.
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL put_refine (ng, model, tile, LputFsur)
          END DO
!$OMP BARRIER
!
!  If refinement, check mass flux conservation between coarser and
!  finer grids.
!
        ELSE IF (isection.eq.nmflx) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL check_massflux (ng, model, tile)
          END DO

# ifndef ONE_WAY
!
!  Fine to coarse coupling: two-way nesting.
!
        ELSE IF (isection.eq.n2way) THEN

          ngc=CoarserDonor(ng)               ! coarse grid number

#  if defined SOLVE3D && !defined NO_CORRECT_TRACER
!
!  Correct coarse grid tracer values at the refinement grid, ng,
!  boundary with the refined accumulated fluxes (Hz*u*T/n, Hz*v*T/m).
!
          DO tile=first_tile(ngc),last_tile(ngc),+1
            CALL correct_tracer (ngc, ng, model, tile)
          END DO
!$OMP BARRIER
#  endif
!
!  Replace coarse grid 2D state variables with the averaged fine grid
!  values (two-way coupling).
!
          DO tile=last_tile(ngc),first_tile(ngc),-1
            CALL fine2coarse (ng, model, r2dvar, tile)
          END DO
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN

#  ifdef SOLVE3D
!
!  Update coarse grid depth variables. We have a new coarse grid
!  adjusted free-surface, Zt_avg1.
!
          DO tile=first_tile(ngc),last_tile(ngc),+1
            CALL set_depth (ngc, tile)
          END DO
!$OMP BARRIER
!
!  Replace coarse grid 3D state variables with the averaged fine grid
!  values (two-way coupling).
!
          DO tile=last_tile(ngc),first_tile(ngc),-1
            CALL fine2coarse (ng, model, r3dvar, tile)
          END DO
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
#  endif
# else
!
!  Fine to coarse coupling (two-way nesting) is not activated!
!
        ELSE IF (isection.eq.n2way) THEN
# endif
        END IF
!
!-----------------------------------------------------------------------
!  Otherwise, process contact points in composite grid.
!-----------------------------------------------------------------------
!
      ELSE
!
!  Get composite grid contact points data from donor grid. It extracts
!  the donor grid cell data necessary to interpolate state variables
!  at each contact point.
!
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL get_composite (ng, model, isection, tile)
        END DO
!$OMP BARRIER
!
!  Fill composite grid contact points variables by interpolating from
!  extracted donor grid data.
!
        DO tile=last_tile(ng),first_tile(ng),-1
          CALL put_composite (ng, model, isection, tile)
        END DO
!$OMP BARRIER

      END IF

      RETURN
      END SUBROUTINE nesting
!
      SUBROUTINE bry_fluxes (dg, rg, cr, model, tile,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       ILB, IUB, JLB, JUB,                        &
     &                       scale, FX, FE,                             &
     &                       F_west, F_east, F_south, F_north)
!
!=======================================================================
!                                                                      !
!  This routine extracts tracer horizontal advective fluxes (Hz*u*T/n, !
!  Hz*v*T/m) at the grid contact boundary (physical domain perimeter). !
!  The data source is either the coarse or finer grid.  These fluxes   !
!  are used for in two-way nesting.         b                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     dg         Donor grid number (integer)                           !
!     rg         Receiver grid number (integer)                        !
!     cr         Contact region number to process (integer)            !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     scale      Advective flux scale (floating-point)                 !
!     IminS      Advective flux, I-dimension Lower bound (integer)     !
!     ImaxS      Advective flux, I-dimension Upper bound (integer)     !
!     JminS      Advective flux, J-dimension Lower bound (integer)     !
!     JmaxS      Advective flux, J-dimension Upper bound (integer)     !
!     ILB        Western/Eastern   boundary flux Lower bound (integer) !
!     IUB        Western/Eastern   boundary flux Upper bound (integer) !
!     JLB        Southern/Northern boundary flux Lower bound (integer) !
!     JUB        Southern/Northern boundary flux Lower bound (integer) !
!     FX         Horizontal advetive flux in the XI-direction (array)  !
!     FE         Horizontal advetive flux in the ETA-direction (array) !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     F_west     Western  boundary advective flux (1D array)           !
!     F_east     Eastern  boundary advective flux (1D array)           !
!     F_south    Southern boundary advective flux (1D array)           !
!     F_north    Northerb boundary advective flux (1D array)           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_nesting
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: dg, rg, cr, model, tile
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: ILB, IUB, JLB, JUB

      real(r8), intent(in) :: scale
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: FX(IminS:,JminS:)
      real(r8), intent(in) :: FE(IminS:,JminS:)

      real(r8), intent(inout) :: F_west (JLB:)
      real(r8), intent(inout) :: F_east (JLB:)
      real(r8), intent(inout) :: F_south(ILB:)
      real(r8), intent(inout) :: F_north(ILB:)
# else
      real(r8), intent(in) :: FX(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: FE(IminS:ImaxS,JminS:JmaxS)

      real(r8), intent(inout) :: F_west (JLB:JUB)
      real(r8), intent(inout) :: F_east (JLB:JUB)
      real(r8), intent(inout) :: F_south(ILB:IUB)
      real(r8), intent(inout) :: F_north(ILB:IUB)
# endif
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: Ib_east, Ib_west, Jb_north, Jb_south
      integer :: i, j, m

# ifdef DISTRIBUTE
      integer :: NptsWE, NptsSN

      real(r8), parameter :: Fspv = 0.0_r8
# endif
!
!-----------------------------------------------------------------------
!  Initialize local variables.
!-----------------------------------------------------------------------
!
!  Set tile starting and ending indices.
!
      Istr=BOUNDS(rg)%Istr(tile)
      Iend=BOUNDS(rg)%Iend(tile)
      Jstr=BOUNDS(rg)%Jstr(tile)
      Jend=BOUNDS(rg)%Jend(tile)

# ifdef DISTRIBUTE
!
!  Initialize arrays to facilitate collective communications.
!
      NptsWE=JUB-JLB+1
      NptsSN=IUB-ILB+1
!
      F_west =Fspv
      F_east =Fspv
      F_south=Fspv
      F_north=Fspv
# endif
!
!-----------------------------------------------------------------------
!  If "rg" is the finer grid, extract advective tracer flux at its
!  physical domain boundaries (grid perimeter).
!-----------------------------------------------------------------------
!
!  Receiver finer grid number is greater than donor coaser grid number
!  because of refinement nesting layers.
!
      IF (rg.gt.dg) THEN
!
!  Western boundary.
!
        IF (DOMAIN(dg)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            F_west(j)=FX(Istr,j)*scale
          END DO
        END IF
!
!  Eastern boundary.
!
        IF (DOMAIN(dg)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            F_east(j)=FX(Iend+1,j)*scale
          END DO
        END IF
!
!  Southern boundary.
!
        IF (DOMAIN(dg)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            F_south(i)=FE(i,Jstr)*scale
          END DO
        END IF
!
!  Northern boundary.
!
        IF (DOMAIN(dg)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            F_north(i)=FE(i,Jend+1)*scale
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  If "rg" is the coarser grid, extract coarser grid advective tracer
!  flux at the location of the finer grid physical domain boundaries
!  (grid perimeter).
!-----------------------------------------------------------------------
!
!  Receiver coarser grid number is smaller than donor finer grid number
!  because of refinement nesting layers.
!
      ELSE IF (rg.lt.dg) THEN
!
!  Western/Eastern boundaries.
!
        Ib_west=I_left(dg)
        Ib_east=I_right(dg)
        DO j=Jstr,Jend
          IF ((Istr.le.Ib_west).and.(Ib_west.le.Iend)) THEN
            F_west(j)=FX(Ib_west,j)*scale
          END IF
!
          IF ((Istr.le.Ib_east).and.(Ib_east.le.Iend)) THEN
            F_east(j)=FX(Ib_east,j)*scale
          END IF
        END DO
!
!  Southern/Northern boundaries.
!
        Jb_south=J_bottom(dg)
        Jb_north=J_top(dg)
        DO i=Istr,Iend
          IF ((Jstr.le.Jb_south).and.(Jb_south.le.Jend)) THEN
            F_south(i)=FE(i,Jb_south)*scale
          END IF
!
          IF ((Jstr.le.Jb_north).and.(Jb_north.le.Jend)) THEN
            F_north(i)=FE(i,Jb_north)*scale
          END IF
        END DO
      END IF

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Gather and broadcast data from all nodes.
!-----------------------------------------------------------------------
!
      CALL mp_assemble (dg, model, NptsWE, Fspv, F_west (JLB:))
      IF (exit_flag.ne.NoError) RETURN

      CALL mp_assemble (dg, model, NptsWE, Fspv, F_east (JLB:))
      IF (exit_flag.ne.NoError) RETURN

      CALL mp_assemble (dg, model, NptsSN, Fspv, F_south(ILB:))
      IF (exit_flag.ne.NoError) RETURN

      CALL mp_assemble (dg, model, NptsSN, Fspv, F_north(ILB:))
      IF (exit_flag.ne.NoError) RETURN
# endif

      RETURN
      END SUBROUTINE bry_fluxes
!
      SUBROUTINE check_massflux (ngf, model, tile)
!
!=======================================================================
!                                                                      !
!  If refinement, this routine check mass fluxes between coarse and    !
!  fine grids for mass and volume conservation. It is only used for    !
!  diagnostic purposes.                                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ngf          Finer grid number (integer)                         !
!     model        Calling model identifier (integer)                  !
!     tile         Domain tile partition (integer)                     !
!                                                                      !
!  On Output:    (mod_nesting)                                         !
!                                                                      !
!     BRY_CONTACT  Updated Mflux in structure.                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_nesting
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ngf, model, tile
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
      integer :: ILB, IUB, JLB, JUB, NptsSN, NptsWE, my_tile
# endif
      integer :: Iedge, Ibc, Ibc_min, Ibc_max, Ibf, Io
      integer :: Jedge, Jbc, Jbc_min, Jbc_max, Jbf, Jo
      integer :: Istr, Iend, Jstr, Jend
      integer :: cjcr, cr, dg, half, icr, isum, jsum, m, rg
      integer :: tnew, told

# ifdef DISTRIBUTE
      real(r8), parameter :: spv = 0.0_r8
# endif
      real(r8) :: EastSum, NorthSum, SouthSum, WestSum
# ifdef NESTING_DEBUG
      real(r8) :: MFratio
# endif
!
!-----------------------------------------------------------------------
!  Check mass and volume conservation during refinement between coarse
!  and fine grids.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
!
!  Get data donor and data receiver grid numbers.
!
        dg=Rcontact(cr)%donor_grid
        rg=Rcontact(cr)%receiver_grid
!
!  Process only contact region data for requested nested finer grid
!  "ngf". Notice that the donor grid is coarser than receiver grid.
!
        IF ((rg.eq.ngf).and.(DXmax(dg).gt.DXmax(rg))) THEN
!
!  Set tile starting and ending indices for donor coarser grid.
!
          Istr=BOUNDS(dg)%Istr(tile)
          Iend=BOUNDS(dg)%Iend(tile)
          Jstr=BOUNDS(dg)%Jstr(tile)
          Jend=BOUNDS(dg)%Jend(tile)
!
!  Set time rolling indices and conjugate region where the coarser
!  donor grid becomes the receiver grid.
!
          told=3-RollingIndex(cr)
          tnew=RollingIndex(cr)
          DO icr=1,Ncontact
            IF ((rg.eq.Rcontact(icr)%donor_grid).and.                   &
     &          (dg.eq.Rcontact(icr)%receiver_grid)) THEN
              cjcr=icr
              EXIT
            END IF
          END DO

# ifdef DISTRIBUTE
!
!  Set global size of boundary edges for coarse grid (donor index).
!
          my_tile=-1
          ILB=BOUNDS(dg)%LBi(my_tile)
          IUB=BOUNDS(dg)%UBi(my_tile)
          JLB=BOUNDS(dg)%LBj(my_tile)
          JUB=BOUNDS(dg)%UBj(my_tile)
          NptsWE=JUB-JLB+1
          NptsSN=IUB-ILB+1
!
!  If distributed-memory, initialize arrays used to check mass flux
!  conservation with special value (zero) to facilitate the global
!  reduction when collecting data between all nodes.
!
          BRY_CONTACT(iwest ,cjcr)%Mflux=spv
          BRY_CONTACT(ieast ,cjcr)%Mflux=spv
          BRY_CONTACT(isouth,cjcr)%Mflux=spv
          BRY_CONTACT(inorth,cjcr)%Mflux=spv
# endif
!
!  Set finer grid center (half) and offset indices (Io and Jo) for
!  coarser grid (I,J) coordinates.
!
          half=(RefineScale(ngf)-1)/2
          Io=half+1
          Jo=half+1
!
!-----------------------------------------------------------------------
!  Average finer grid western boundary mass fluxes and load them to the
!  BRY_CONTACT structure.
!-----------------------------------------------------------------------
!
          Ibc=I_left(ngf)
          Jbc_min=J_bottom(ngf)
          Jbc_max=J_top(ngf)-1            ! interior points, no top
!                                           left corner
# ifdef NESTING_DEBUG
          IF (DOMAIN(ngf)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (300,10) 'Western Boundary Mass Fluxes: ',          &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
            END IF
          END IF
!
# endif
          DO Jbc=Jstr,Jend
            IF (((Istr.le.Ibc).and.(Ibc.le.Iend)).and.                  &
     &          ((Jbc_min.le.Jbc).and.(Jbc.le.Jbc_max))) THEN
!
!  Sum finer grid western boundary mass fluxes within coarser grid cell.
!
              WestSum=0.0_r8
              Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
              DO jsum=-half,half
                Jbf=Jedge+jsum
                WestSum=WestSum+BRY_CONTACT(iwest,cr)%Mflux(Jbf)
              END DO
              m=BRY_CONTACT(iwest,cr)%C2Bindex(Jbf)     ! pick last one
!
!  Load coarser grid western boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
              BRY_CONTACT(iwest,cjcr)%Mflux(Jbc)=WestSum

# ifdef NESTING_DEBUG
              IF (WestSum.ne.0) THEN
                MFratio=REFINED(cr)%DU_avg2(1,m,tnew)/WestSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (300,30) Jbc, REFINED(cr)%DU_avg2(1,m,tnew),        &
     &                       WestSum, MFratio
# endif
            END IF
          END DO
!
!-----------------------------------------------------------------------
!  Average finer grid eastern boundary mass fluxes and load them to the
!  BRY_CONTACT structure.
!-----------------------------------------------------------------------
!
          Ibc=I_right(ngf)
          Jbc_min=J_bottom(ngf)
          Jbc_max=J_top(ngf)-1            ! interior points, no top
!                                           right corner
# ifdef NESTING_DEBUG
          IF (DOMAIN(ngf)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (300,10) 'Eastern Boundary Mass Fluxes: ',          &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
            END IF
          END IF
!
# endif
          DO Jbc=Jstr,Jend
            IF (((Istr.le.Ibc).and.(Ibc.le.Iend)).and.                  &
     &          ((Jbc_min.le.Jbc).and.(Jbc.le.Jbc_max))) THEN
!
!  Sum finer grid eastern boundary mass fluxes within coarser grid cell.
!
              EastSum=0.0_r8
              Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
              DO jsum=-half,half
                Jbf=Jedge+jsum
                EastSum=EastSum+BRY_CONTACT(ieast,cr)%Mflux(Jbf)
              END DO
              m=BRY_CONTACT(ieast,cr)%C2Bindex(Jbf)     ! pick last one
!
!  Load coarser grid eastern boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
              BRY_CONTACT(ieast,cjcr)%Mflux(Jbc)=EastSum

# ifdef NESTING_DEBUG
              IF (EastSum.ne.0) THEN
                MFratio=REFINED(cr)%DU_avg2(1,m,tnew)/EastSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (300,30) Jbc, REFINED(cr)%DU_avg2(1,m,tnew),        &
     &                       EastSum, MFratio
# endif
            END IF
          END DO
!
!-----------------------------------------------------------------------
!  Average finer grid southern boundary mass fluxes and load them to the
!  BRY_CONTACT structure.
!-----------------------------------------------------------------------
!
          Jbc=J_bottom(ngf)
          Ibc_min=I_left(ngf)
          Ibc_max=I_right(ngf)-1          ! interior points, no bottom
!                                           right corner
# ifdef NESTING_DEBUG
          IF (DOMAIN(ngf)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (300,20) 'Southern Boundary Mass Fluxes: ',         &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
            END IF
          END IF
!
# endif
          DO Ibc=Istr,Iend
            IF (((Ibc_min.le.Ibc).and.(Ibc.le.Ibc_max)).and.            &
     &          ((Jstr.le.Jbc).and.(Jbc.le.Jend))) THEN
!
!  Sum finer grid southern boundary mass fluxes within coarser grid
!  cell.
!
              SouthSum=0.0_r8
              Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
              DO isum=-half,half
                Ibf=Iedge+isum
                SouthSum=SouthSum+BRY_CONTACT(isouth,cr)%Mflux(Ibf)
              END DO
              m=BRY_CONTACT(isouth,cr)%C2Bindex(Ibf)    ! pick last one
!
!  Load coarser grid southern boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
              BRY_CONTACT(isouth,cjcr)%Mflux(Ibc)=SouthSum

# ifdef NESTING_DEBUG
              IF (SouthSum.ne.0) THEN
                MFratio=REFINED(cr)%DV_avg2(1,m,tnew)/SouthSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (300,30) Ibc, REFINED(cr)%DV_avg2(1,m,tnew),        &
     &                       SouthSum, MFratio
# endif
            END IF
          END DO
!
!-----------------------------------------------------------------------
!  Average finer grid northern boundary mass fluxes and load them to the
!  BRY_CONTACT structure.
!-----------------------------------------------------------------------
!
          Jbc=J_top(ngf)
          Ibc_min=I_left(ngf)
          Ibc_max=I_right(ngf)-1          ! interior points, no top
!                                           right corner
# ifdef NESTING_DEBUG
          IF (DOMAIN(ngf)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (300,20) 'Northern Boundary Mass Fluxes: ',         &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
            END IF
          END IF
!
# endif
          DO Ibc=Istr,Iend
            IF (((Ibc_min.le.Ibc).and.(Ibc.le.Ibc_max)).and.            &
     &          ((Jstr.le.Jbc).and.(Jbc.le.Jend))) THEN
!
!  Sum finer grid northern boundary mass fluxes within coarser grid
!  cell.
!
              NorthSum=0.0_r8
              Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
              DO isum=-half,half
                Ibf=Iedge+isum
                NorthSum=NorthSum+BRY_CONTACT(inorth,cr)%Mflux(Ibf)
              END DO
              m=BRY_CONTACT(inorth,cr)%C2Bindex(Ibf)    ! pick last one
!
!  Load coarser grid northern boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
              BRY_CONTACT(inorth,cjcr)%Mflux(Ibc)=NorthSum

# ifdef NESTING_DEBUG
              IF (NorthSum.ne.0) THEN
                MFratio=REFINED(cr)%DV_avg2(1,m,tnew)/NorthSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (300,30) Ibc, REFINED(cr)%DV_avg2(1,m,tnew),        &
     &                       NorthSum, MFratio
# endif
            END IF
          END DO

# ifdef DISTRIBUTE
!
!  Collect data from all nodes.
!
          CALL mp_assemble (dg, model, NptsWE, spv,                      &
     &                      BRY_CONTACT(iwest ,cjcr)%Mflux(JLB:))
          CALL mp_assemble (dg, model, NptsWE, spv,                      &
     &                      BRY_CONTACT(ieast ,cjcr)%Mflux(JLB:))
          CALL mp_assemble (dg, model, NptsSN, spv,                      &
     &                      BRY_CONTACT(isouth,cjcr)%Mflux(ILB:))
          CALL mp_assemble (dg, model, NptsSN, spv,                      &
     &                      BRY_CONTACT(inorth,cjcr)%Mflux(ILB:))
# endif
        END IF
      END DO

# ifdef NESTING_DEBUG
!
      CALL my_flush (300)
!
  10  FORMAT (/,1x,a,/,/,4x,'cr = ',i2.2,4x,'dg = ',i2.2,4x,'rg = ',    &
     &        i2.2,4x,'iif(rg) = ',i3.3,4x,'iic(rg) = ',i7.7,4x,        &
     &        'time(rg) = ',i7.7,/,/,2x,'Coarse',6x,'Coarse Grid',8x,   &
     &        'Fine Grid',11x,'Ratio',/,4x,'Jb',9x,'DU_avg2',9x,        &
     &        'SUM(DU_avg2)',/)
  20  FORMAT (/,1x,a,/,/,4x,'cr = ',i2.2,4x,'dg = ',i2.2,4x,'rg = ',    &
     &        i2.2,4x,'iif(rg) = ',i3.3,4x,'iic(rg) = ',i7.7,4x,        &
     &        'time(rg) = ',i7.7,/,/,2x,'Coarse',6x,'Coarse Grid',8x,   &
     &        'Fine Grid',11x,'Ratio',/,4x,'Ib',9x,'DV_avg2',9x,        &
     &        'SUM(DV_avg2)',/)
  30  FORMAT (4x,i4.4,3(3x,1p,e15.8))
# endif

      RETURN
      END SUBROUTINE check_massflux

# ifndef ONE_WAY
!
      FUNCTION do_twoway (nl, il, ng, istep) RESULT (doit)
!
!=======================================================================
!                                                                      !
!  This function determines if the two-way exchange between finer to   !
!  coarser grid is appropriate. In complex nesting applications with   !
!  telescoping grids, grid with RefineScale > 0 including finer grids  !
!  inside, the two-way feedback between finer and coarse grids needs   !
!  to be in the correct sequence.                                      !
!                                                                      !
!  This function is called from either main2d or main3d.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     nl         Latest time-stepped nested layer index (integer)      !
!     il         Current nested layer index in the two-way DO loop     !
!                  call from main2d or main3d (integer)                !
!     ng         Finer grid number to process in the fine to coarse    !
!                  feedback (integer)                                  !
!     istep      Current rime-step counter from STEP_LOOP in main2d    !
!                  or main3d (integer)                                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     doit       The value of the result is TRUE/FALSE if the exchage  !
!                  is required or not (logical)                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_nesting
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: nl, il, ng, istep
!
!  Local variable declarations.
!
      logical :: doit

      integer :: dgc, dgf
!
!-----------------------------------------------------------------------
!  Determine if two-way exchage is required at the corrent time-step.
!-----------------------------------------------------------------------
!
      doit=.FALSE.
!
!  Chech if ng is a refined grid, RefineScale(ng) > 0.
!
      IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
!
!  If telescoping grid, the sequence of exchanges when all the nested
!  grid reach the coarser grid (ng=1) time needs to be done in the
!  correct sequence: finest grid exchage first.
!
        IF (Telescoping(ng)) THEN
          dgc=coarserDonor(ng)                        ! for debbuging
          dgf=FinerDonor(ng)
          IF (TwoWayCounter(dgf).eq.RefineSteps(dgf)) THEN
            doit=.TRUE.
          END IF
!
!  If not telescoping, the exchange is in terms of the "istep" counter
!  but also at the right sequence.
!
        ELSE
          IF ((istep.eq.RefineSteps(ng)).and.(il.eq.nl)) THEN
            doit=.TRUE.
          END IF
        END IF
      END IF
!
!  Update two-way exchange counter between fine and coarse grids.
!
      IF (doit) THEN
        TwoWayCounter(ng)=TwoWayCounter(ng)+1
      END IF

      RETURN
      END FUNCTION do_twoway
# endif
!
      SUBROUTINE fill_contact (rg, model, tile,                         &
     &                         cr, Npoints, contact,                    &
     &                         gtype, mvname, SpValCheck,               &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Ac, Ar)
!
!=======================================================================
!                                                                      !
!  This routine is used during initialization to fill the contact      !
!  points of a specified grid metric array. We need to have metric     !
!  data in all the extended computational points of the grid. No       !
!  attempt is done here to interpolate such values since the are       !
!  read in "set_contact" from input contact points NetCDF file.        !
!  This routine just unpack data into global arrays and check if       !
!  all needed values are filled.                                       !
!                                                                      !
!  During allocation these special metric grid arrays are initialized  !
!  to "spval" to avoid resetting those values processed already from   !
!  the regular Grid NetCDF file.  That is, only those contact points   !
!  outside the physical grid are processed here.  This is a good way   !
!  to check if all the extra numerical points have been processed.     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     rg         Receiver grid number (integer)                        !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     cr         Contact region number to process (integer)            !
!     Npoints    Number of points in the contact region (integer)      !
!     contact    Contact region information variables (T_NGC structure)!
!     gtype      C-grid variable type (integer)                        !
!     mvname     Metrics variable name (string)                        !
!     SpValCheck Special value used to check if the contact point      !
!                  needs to be processed.                              !
!     LBi        Receiver grid, I-dimension Lower bound (integer)      !
!     UBi        Receiver grid, I-dimension Upper bound (integer)      !
!     LBj        Receiver grid, J-dimension Lower bound (integer)      !
!     UBj        Receiver grid, J-dimension Upper bound (integer)      !
!     Ac         Metric data at Contact points.                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ar         Updated metric grid array                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_nesting
!
!  Imported variable declarations.
!
      integer, intent(in) :: rg, model, tile
      integer, intent(in) :: cr, gtype, Npoints
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(in) :: SpValCheck
!
      character(len=*), intent(in) :: mvname
!
      TYPE (T_NGC), intent(in) :: contact(:)
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Ac(:)
      real(r8), intent(inout) :: Ar(LBi:,LBj:)
# else
      real(r8), intent(in) :: Ac(Npoints)
      real(r8), intent(inout) :: Ar(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, m
      integer :: Istr, Iend, Jstr, Jend
!
!-----------------------------------------------------------------------
!  Interpolate 2D data from donor grid to receiver grid contact points.
!-----------------------------------------------------------------------
!
!
!  Set starting and ending tile indices for the receiver grid.
!
      SELECT CASE (gtype)
        CASE (p2dvar)
          Istr=BOUNDS(rg) % IstrP(tile)
          Iend=BOUNDS(rg) % IendP(tile)
          Jstr=BOUNDS(rg) % JstrP(tile)
          Jend=BOUNDS(rg) % JendP(tile)
        CASE (r2dvar)
          Istr=BOUNDS(rg) % IstrT(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrT(tile)
          Jend=BOUNDS(rg) % JendT(tile)
        CASE (u2dvar)
          Istr=BOUNDS(rg) % IstrP(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrT(tile)
          Jend=BOUNDS(rg) % JendT(tile)
        CASE (v2dvar)
          Istr=BOUNDS(rg) % IstrT(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrP(tile)
          Jend=BOUNDS(rg) % JendT(tile)
      END SELECT
!
!  Interpolate.
!
      DO m=1,Npoints
        i=contact(cr)%Irg(m)
        j=contact(cr)%Jrg(m)
        IF (((Istr.le.i).and.(i.le.Iend)).and.                          &
     &      ((Jstr.le.j).and.(j.le.Jend))) THEN
!!        IF (Ar(i,j).gt.SpValCheck) THEN        ! Only process contact
            Ar(i,j)=Ac(m)                        ! points outside in the
!!        END IF                                 ! regular grid
        END IF
      END DO

      RETURN
      END SUBROUTINE fill_contact

# if defined MASKING || defined WET_DRY
!
      SUBROUTINE mask_hweights (ng, model, tile)
!
!=======================================================================
!                                                                      !
!  This routine scales the horizontal interpolation weights to account !
!  for Land/Sea masking in  the nested contact region.  If wetting and !
!  drying,  the scaling is done at every time step since the  Land/Sea !
!  masking is time dependent.                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     rg         Receiver grid number (integer)                        !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Lweight    Updated linear interpolation weights                  !
#  ifdef QUADRATIC_WEIGHTS
!     Qweight    Updated quadratic interpolation weights               !
#  endif
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_nesting
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
!
!  Local variable declarations.
!
      integer :: cr, dg,  m, rg
      integer :: Istr, Iend, Jstr, Jend
      integer :: Idg, Idgp1, Jdg, Jdgp1
      integer :: NpointsR, NpointsU, NpointsV
#  ifdef DISTRIBUTE
      integer :: Lpoints, Qpoints
#  endif

      real(r8) :: cff
      real(r8) :: LWsum, LWsum_check, MaskLWsum
#  ifdef QUADRATIC_WEIGHTS
      real(r8) :: QWsum, QWsum_check, MaskQWsum
#  endif
#  ifdef DISTRIBUTE
      real(r8), parameter :: spv = 0.0_r8
#  endif
      real(r8), dimension(4) :: Lweight
#  ifdef QUADRATIC_WEIGHTS
      real(r8), dimension(9) :: Qweight
#  endif
#  ifdef DISTRIBUTE
      real(r8), allocatable :: LW(:,:)
#   ifdef QUADRATIC_WEIGHTS
      real(r8), allocatable :: QW(:,:)
#   endif
#  endif
!
!=======================================================================
!  If appropriate, scale horizontal interpolation weights to account
!  for Land/Sea masking.
!=======================================================================
!
      DO cr=1,Ncontact
!
!  Get data donor and data receiver grid numbers.
!
        dg=Rcontact(cr)%donor_grid
        rg=Rcontact(cr)%receiver_grid
!
        SCALE_WEIGHTS : IF (RefinedGrid(rg).and.                        &
     &                      (dg.eq.ng).and.(dg.lt.rg)) THEN
!
!  Scale interpolation weigths for RHO-contact points.
!
          Istr=BOUNDS(dg) % Istr(tile)
          Iend=BOUNDS(dg) % Iend(tile)
          Jstr=BOUNDS(dg) % Jstr(tile)
          Jend=BOUNDS(dg) % Jend(tile)
!
!-----------------------------------------------------------------------
!  Scale horizontal interpolation weigths for RHO-contact points.
!-----------------------------------------------------------------------
!
          NpointsR=Rcontact(cr)%Npoints

#  ifdef DISTRIBUTE
!
!  If distributed-memory, allocate and initialize working arrays
!  with special value (zero) to facilitate the global reduction
!  when collecting data between all nodes.
!
          Lpoints=4*NpointsR
          IF (.not.allocated(LW)) THEN
            allocate ( LW(4,NpointsR) )
          END IF
          LW=spv

#   ifdef QUADRATIC_WEIGHTS
          Qpoints=9*NpointsR
          IF (.not.allocated(QW)) THEN
            allocate ( QW(9,NpointsR) )
          END IF
          QW=spv
#   endif
#  endif
!
!  Scale interpolation weights.
!
          DO m=1,NpointsR
            Idg  =Rcontact(cr)%Idg(m)
            Idgp1=MIN(Idg+1, BOUNDS(dg)%UBi(-1))
            Jdg  =Rcontact(cr)%Jdg(m)
            Jdgp1=MIN(Jdg+1, BOUNDS(dg)%UBj(-1))
            IF (((Istr.le.Idg).and.(Idg.le.Iend)).and.                  &
     &          ((Jstr.le.Jdg).and.(Jdg.le.Jend))) THEN
!
!  Linear interpolation weights.
!
              MaskLWsum=GRID(dg)%rmask_full(Idg  ,Jdg  )+               &
     &                  GRID(dg)%rmask_full(Idgp1,Jdg  )+               &
     &                  GRID(dg)%rmask_full(Idgp1,Jdgp1)+               &
     &                  GRID(dg)%rmask_full(Idg  ,Jdgp1)
              IF (MaskLWsum.lt.4) THEN
#  ifdef WET_DRY
                Lweight(1)=Rcontact(cr)%LweightUnmasked(1,m)*           &
     &                     GRID(dg)%rmask_full(Idg  ,Jdg  )
                Lweight(2)=Rcontact(cr)%LweightUnmasked(2,m)*           &
     &                     GRID(dg)%rmask_full(Idgp1,Jdg  )
                Lweight(3)=Rcontact(cr)%LweightUnmasked(3,m)*           &
     &                     GRID(dg)%rmask_full(Idgp1,Jdgp1)
                Lweight(4)=Rcontact(cr)%LweightUnmasked(4,m)*           &
     &                     GRID(dg)%rmask_full(Idg  ,Jdgp1)
#  else
                Lweight(1)=Rcontact(cr)%Lweight(1,m)*                   &
     &                     GRID(dg)%rmask_full(Idg  ,Jdg  )
                Lweight(2)=Rcontact(cr)%Lweight(2,m)*                   &
     &                     GRID(dg)%rmask_full(Idgp1,Jdg  )
                Lweight(3)=Rcontact(cr)%Lweight(3,m)*                   &
     &                     GRID(dg)%rmask_full(Idgp1,Jdgp1)
                Lweight(4)=Rcontact(cr)%Lweight(4,m)*                   &
     &                     GRID(dg)%rmask_full(Idg  ,Jdgp1)
#  endif
                LWsum=SUM(Lweight)
                IF (LWsum.gt.0) THEN
                  cff=1.0_r8/LWsum
                  Lweight(1)=cff*Lweight(1)    ! using only water points
                  Lweight(2)=cff*Lweight(2)    ! in the interpolation of
                  Lweight(3)=cff*Lweight(3)    ! of contact points
                  Lweight(4)=cff*Lweight(4)
                ELSE
                  Lweight=0.0_r8               ! all the donor grid
                END IF                         ! corners are on land
                LWsum_check=SUM(Lweight)
#  ifdef DISTRIBUTE
                LW(1,m)=Lweight(1)
                LW(2,m)=Lweight(2)
                LW(3,m)=Lweight(3)
                LW(4,m)=Lweight(4)
              ELSE
                LW(1,m)=Rcontact(cr)%Lweight(1,m)   ! we need to load
                LW(2,m)=Rcontact(cr)%Lweight(2,m)   ! unscaled values
                LW(3,m)=Rcontact(cr)%Lweight(3,m)   ! to facilitate
                LW(4,m)=Rcontact(cr)%Lweight(4,m)   ! parallel reduction
#  else
                Rcontact(cr)%Lweight(1,m)=Lweight(1)
                Rcontact(cr)%Lweight(2,m)=Lweight(2)
                Rcontact(cr)%Lweight(3,m)=Lweight(3)
                Rcontact(cr)%Lweight(4,m)=Lweight(4)
#  endif
              END IF

#  ifdef QUADRATIC_WEIGHTS
!
!  Quadratic interpolation weights.
!
              MaskQWsum=GRID(dg)%rmask_full(Idg-1, Jdg-1)+              &
     &                  GRID(dg)%rmask_full(Idg  , Jdg-1)+              &
     &                  GRID(dg)%rmask_full(Idgp1, Jdg-1)+              &
     &                  GRID(dg)%rmask_full(Idg-1, Jdg  )+              &
     &                  GRID(dg)%rmask_full(Idg  , Jdg  )+              &
     &                  GRID(dg)%rmask_full(Idgp1, Jdg  )+              &
     &                  GRID(dg)%rmask_full(Idg-1, Jdgp1)+              &
     &                  GRID(dg)%rmask_full(Idg  , Jdgp1)+              &
     &                  GRID(dg)%rmask_full(Idgp1, Jdgp1)
              IF (MaskQWsum.lt.9) THEN
#   ifdef WET_DRY
                Qweight(1)=Rcontact(cr)%QweightUnmasked(1,m)*           &
     &                     GRID(dg)%rmask_full(Idg-1,Jdg-1)
                Qweight(2)=Rcontact(cr)%QweightUnmasked(2,m)*           &
     &                     GRID(dg)%rmask_full(Idg  ,Jdg-1)
                Qweight(3)=Rcontact(cr)%QweightUnmasked(3,m)*           &
     &                     GRID(dg)%rmask_full(Idgp1,Jdg-1)
                Qweight(4)=Rcontact(cr)%QweightUnmasked(4,m)*           &
     &                     GRID(dg)%rmask_full(Idg-1,Jdg  )
                Qweight(5)=Rcontact(cr)%QweightUnmasked(5,m)*           &
     &                     GRID(dg)%rmask_full(Idg  ,Jdg  )
                Qweight(6)=Rcontact(cr)%QweightUnmasked(6,m)*           &
     &                     GRID(dg)%rmask_full(Idgp1,Jdg  )
                Qweight(7)=Rcontact(cr)%QweightUnmasked(7,m)*           &
     &                     GRID(dg)%rmask_full(Idg-1,Jdgp1)
                Qweight(8)=Rcontact(cr)%QweightUnmasked(8,m)*           &
     &                     GRID(dg)%rmask_full(Idg  ,Jdgp1)
                Qweight(9)=Rcontact(cr)%QweightUnmasked(9,m)*           &
     &                     GRID(dg)%rmask_full(Idgp1,Jdgp1)
#   else
                Qweight(1)=Rcontact(cr)%Qweight(1,m)*                   &
     &                     GRID(dg)%rmask_full(Idg-1,Jdg-1)
                Qweight(2)=Rcontact(cr)%Qweight(2,m)*                   &
     &                     GRID(dg)%rmask_full(Idg  ,Jdg-1)
                Qweight(3)=Rcontact(cr)%Qweight(3,m)*                   &
     &                     GRID(dg)%rmask_full(Idgp1,Jdg-1)
                Qweight(4)=Rcontact(cr)%Qweight(4,m)*                   &
     &                     GRID(dg)%rmask_full(Idg-1,Jdg  )
                Qweight(5)=Rcontact(cr)%Qweight(5,m)*                   &
     &                     GRID(dg)%rmask_full(Idg  ,Jdg  )
                Qweight(6)=Rcontact(cr)%Qweight(6,m)*                   &
     &                     GRID(dg)%rmask_full(Idgp1,Jdg  )
                Qweight(7)=Rcontact(cr)%Qweight(7,m)*                   &
     &                     GRID(dg)%rmask_full(Idg-1,Jdgp1)
                Qweight(8)=Rcontact(cr)%Qweight(8,m)*                   &
     &                     GRID(dg)%rmask_full(Idg  ,Jdgp1)
                Qweight(9)=Rcontact(cr)%Qweight(9,m)*                   &
     &                     GRID(dg)%rmask_full(Idgp1,Jdgp1)
#   endif
                QWsum=SUM(Qweight)
                IF (QWsum.gt.0) THEN
                  cff=1.0_r8/QWsum
                  Qweight(1)=cff*Qweight(1)    ! using only water points
                  Qweight(2)=cff*Qweight(2)    ! in the interpolation of
                  Qweight(3)=cff*Qweight(3)    ! of contact points
                  Qweight(4)=cff*Qweight(4)
                  Qweight(5)=cff*Qweight(5)
                  Qweight(6)=cff*Qweight(6)
                  Qweight(7)=cff*Qweight(7)
                  Qweight(8)=cff*Qweight(8)
                  Qweight(9)=cff*Qweight(9)
                ELSE
                  Qweight=0.0_r8               ! all the donor grid
                END IF                         ! corners are on land
                QWsum_check=SUM(Qweight)
#   ifdef DISTRIBUTE
                QW(1,m)=Qweight(1)
                QW(2,m)=Qweight(2)
                QW(3,m)=Qweight(3)
                QW(4,m)=Qweight(4)
                QW(5,m)=Qweight(5)
                QW(6,m)=Qweight(6)
                QW(7,m)=Qweight(7)
                QW(8,m)=Qweight(8)
                QW(9,m)=Qweight(9)
              ELSE
                QW(1,m)=Rcontact(cr)%Qweight(1,m)   ! we need to load
                QW(2,m)=Rcontact(cr)%Qweight(2,m)   ! unscaled values
                QW(3,m)=Rcontact(cr)%Qweight(3,m)   ! to facilitate
                QW(4,m)=Rcontact(cr)%Qweight(4,m)   ! parallel reduction
                QW(5,m)=Rcontact(cr)%Qweight(5,m)
                QW(6,m)=Rcontact(cr)%Qweight(6,m)
                QW(7,m)=Rcontact(cr)%Qweight(7,m)
                QW(8,m)=Rcontact(cr)%Qweight(8,m)
                QW(9,m)=Rcontact(cr)%Qweight(9,m)
#   else
                Rcontact(cr)%Qweight(1,m)=Qweight(1)
                Rcontact(cr)%Qweight(2,m)=Qweight(2)
                Rcontact(cr)%Qweight(3,m)=Qweight(3)
                Rcontact(cr)%Qweight(4,m)=Qweight(4)
                Rcontact(cr)%Qweight(5,m)=Qweight(5)
                Rcontact(cr)%Qweight(6,m)=Qweight(6)
                Rcontact(cr)%Qweight(7,m)=Qweight(7)
                Rcontact(cr)%Qweight(8,m)=Qweight(8)
                Rcontact(cr)%Qweight(9,m)=Qweight(9)
#   endif
              END IF
#  endif
            END IF
          END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
          CALL mp_assemble (dg, model, Lpoints, spv, LW)
          IF (exit_flag.ne.NoError) RETURN

#   ifdef QUADRATIC_WEIGHTS
          CALL mp_assemble (dg, model, Qpoints, spv, QW)
          IF (exit_flag.ne.NoError) RETURN
#   endif
!
!  Load exchanged weights.
!
          Rcontact(cr)%Lweight(1:4,1:NpointsR)=LW(1:4,1:NpointsR)
#   ifdef QUADRATIC_WEIGHTS
          Rcontact(cr)%Qweight(1:9,1:NpointsR)=QW(1:9,1:NpointsR)
#   endif
!
!  Deallocate local work arrays.
!
          IF (allocated(LW)) THEN
            deallocate (LW)
          END IF

#   ifdef QUADRATIC_WEIGHTS
          IF (allocated(QW)) THEN
            deallocate (QW)
          END IF
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Scale interpolation weigths for U-contact points.
!-----------------------------------------------------------------------
!
          NpointsU=Ucontact(cr)%Npoints

#  ifdef DISTRIBUTE
!
!  If distributed-memory, allocate and initialize working arrays
!  with special value (zero) to facilitate the global reduction
!  when collecting data between all nodes.
!
          Lpoints=4*NpointsU
          IF (.not.allocated(LW)) THEN
            allocate ( LW(4,NpointsU) )
          END IF
          LW=spv

#   ifdef QUADRATIC_WEIGHTS
          Qpoints=9*NpointsU
          IF (.not.allocated(QW)) THEN
            allocate ( QW(9,NpointsU) )
          END IF
          QW=spv
#   endif
#  endif
!
!  Scale interpolation weights.
!
          DO m=1,NpointsU
            Idg  =Ucontact(cr)%Idg(m)
            Idgp1=MIN(Idg+1, BOUNDS(dg)%UBi(-1))
            Jdg  =Ucontact(cr)%Jdg(m)
            Jdgp1=MIN(Jdg+1, BOUNDS(dg)%UBj(-1))
            IF (((Istr.le.Idg).and.(Idg.le.Iend)).and.                  &
     &          ((Jstr.le.Jdg).and.(Jdg.le.Jend))) THEN
!
!  Linear interpolation weights.
!
              MaskLWsum=GRID(dg)%umask_full(Idg  ,Jdg  )+               &
     &                  GRID(dg)%umask_full(Idgp1,Jdg  )+               &
     &                  GRID(dg)%umask_full(Idgp1,Jdgp1)+               &
     &                  GRID(dg)%umask_full(Idg  ,Jdgp1)
              IF (MaskLWsum.lt.4) THEN
#  ifdef WET_DRY
                Lweight(1)=Ucontact(cr)%LweightUnmasked(1,m)*           &
     &                     GRID(dg)%umask_full(Idg  ,Jdg  )
                Lweight(2)=Ucontact(cr)%LweightUnmasked(2,m)*           &
     &                     GRID(dg)%umask_full(Idgp1,Jdg  )
                Lweight(3)=Ucontact(cr)%LweightUnmasked(3,m)*           &
     &                     GRID(dg)%umask_full(Idgp1,Jdgp1)
                Lweight(4)=Ucontact(cr)%LweightUnmasked(4,m)*           &
     &                     GRID(dg)%umask_full(Idg  ,Jdgp1)
#  else
                Lweight(1)=Ucontact(cr)%Lweight(1,m)*                   &
     &                     GRID(dg)%umask_full(Idg  ,Jdg  )
                Lweight(2)=Ucontact(cr)%Lweight(2,m)*                   &
     &                     GRID(dg)%umask_full(Idgp1,Jdg  )
                Lweight(3)=Ucontact(cr)%Lweight(3,m)*                   &
     &                     GRID(dg)%umask_full(Idgp1,Jdgp1)
                Lweight(4)=Ucontact(cr)%Lweight(4,m)*                   &
     &                     GRID(dg)%umask_full(Idg  ,Jdgp1)
#  endif
                LWsum=SUM(Lweight)
                IF (LWsum.gt.0) THEN
                  cff=1.0_r8/LWsum
                  Lweight(1)=cff*Lweight(1)    ! using only water points
                  Lweight(2)=cff*Lweight(2)    ! in the interpolation of
                  Lweight(3)=cff*Lweight(3)    ! of contact points
                  Lweight(4)=cff*Lweight(4)
                ELSE
                  Lweight=0.0_r8               ! All the donor grid
                END IF                         ! corners are on land
                LWsum_check=SUM(Lweight)
#  ifdef DISTRIBUTE
                LW(1,m)=Lweight(1)
                LW(2,m)=Lweight(2)
                LW(3,m)=Lweight(3)
                LW(4,m)=Lweight(4)
              ELSE
                LW(1,m)=Ucontact(cr)%Lweight(1,m)   ! we need to load
                LW(2,m)=Ucontact(cr)%Lweight(2,m)   ! unscaled values
                LW(3,m)=Ucontact(cr)%Lweight(3,m)   ! to facilitate
                LW(4,m)=Ucontact(cr)%Lweight(4,m)   ! parallel reduction
#  else
                Ucontact(cr)%Lweight(1,m)=Lweight(1)
                Ucontact(cr)%Lweight(2,m)=Lweight(2)
                Ucontact(cr)%Lweight(3,m)=Lweight(3)
                Ucontact(cr)%Lweight(4,m)=Lweight(4)
#  endif
              END IF

#  ifdef QUADRATIC_WEIGHTS
!
!  Quadratic interpolation weights.
!
              MaskQWsum=GRID(dg)%umask_full(Idg-1, Jdg-1)+              &
     &                  GRID(dg)%umask_full(Idg  , Jdg-1)+              &
     &                  GRID(dg)%umask_full(Idgp1, Jdg-1)+              &
     &                  GRID(dg)%umask_full(Idg-1, Jdg  )+              &
     &                  GRID(dg)%umask_full(Idg  , Jdg  )+              &
     &                  GRID(dg)%umask_full(Idgp1, Jdg  )+              &
     &                  GRID(dg)%umask_full(Idg-1, Jdgp1)+              &
     &                  GRID(dg)%umask_full(Idg  , Jdgp1)+              &
     &                  GRID(dg)%umask_full(Idgp1, Jdgp1)
              IF (MaskQWsum.lt.9) THEN
#   ifdef WET_DRY
                Qweight(1)=Ucontact(cr)%QweightUnmasked(1,m)*           &
     &                     GRID(dg)%umask_full(Idg-1,Jdg-1)
                Qweight(2)=Ucontact(cr)%QweightUnmasked(2,m)*           &
     &                     GRID(dg)%umask_full(Idg  ,Jdg-1)
                Qweight(3)=Ucontact(cr)%QweightUnmasked(3,m)*           &
     &                     GRID(dg)%umask_full(Idgp1,Jdg-1)
                Qweight(4)=Ucontact(cr)%QweightUnmasked(4,m)*           &
     &                     GRID(dg)%umask_full(Idg-1,Jdg  )
                Qweight(5)=Ucontact(cr)%QweightUnmasked(5,m)*           &
     &                     GRID(dg)%umask_full(Idg  ,Jdg  )
                Qweight(6)=Ucontact(cr)%QweightUnmasked(6,m)*           &
     &                     GRID(dg)%umask_full(Idgp1,Jdg  )
                Qweight(7)=Ucontact(cr)%QweightUnmasked(7,m)*           &
     &                     GRID(dg)%umask_full(Idg-1,Jdgp1)
                Qweight(8)=Ucontact(cr)%QweightUnmasked(8,m)*           &
     &                     GRID(dg)%umask_full(Idg  ,Jdgp1)
                Qweight(9)=Ucontact(cr)%QweightUnmasked(9,m)*           &
     &                     GRID(dg)%umask_full(Idgp1,Jdgp1)
#   else
                Qweight(1)=Ucontact(cr)%Qweight(1,m)*                   &
     &                     GRID(dg)%umask_full(Idg-1,Jdg-1)
                Qweight(2)=Ucontact(cr)%Qweight(2,m)*                   &
     &                     GRID(dg)%umask_full(Idg  ,Jdg-1)
                Qweight(3)=Ucontact(cr)%Qweight(3,m)*                   &
     &                     GRID(dg)%umask_full(Idgp1,Jdg-1)
                Qweight(4)=Ucontact(cr)%Qweight(4,m)*                   &
     &                     GRID(dg)%umask_full(Idg-1,Jdg  )
                Qweight(5)=Ucontact(cr)%Qweight(5,m)*                   &
     &                     GRID(dg)%umask_full(Idg  ,Jdg  )
                Qweight(6)=Ucontact(cr)%Qweight(6,m)*                   &
     &                     GRID(dg)%umask_full(Idgp1,Jdg  )
                Qweight(7)=Ucontact(cr)%Qweight(7,m)*                   &
     &                     GRID(dg)%umask_full(Idg-1,Jdgp1)
                Qweight(8)=Ucontact(cr)%Qweight(8,m)*                   &
     &                     GRID(dg)%umask_full(Idg  ,Jdgp1)
                Qweight(9)=Ucontact(cr)%Qweight(9,m)*                   &
     &                     GRID(dg)%umask_full(Idgp1,Jdgp1)
#   endif
                QWsum=SUM(Qweight)
                IF (QWsum.gt.0) THEN
                  cff=1.0_r8/QWsum
                  Qweight(1)=cff*Qweight(1)    ! using only water points
                  Qweight(2)=cff*Qweight(2)    ! in the interpolation of
                  Qweight(3)=cff*Qweight(3)    ! of contact points
                  Qweight(4)=cff*Qweight(4)
                  Qweight(5)=cff*Qweight(5)
                  Qweight(6)=cff*Qweight(6)
                  Qweight(7)=cff*Qweight(7)
                  Qweight(8)=cff*Qweight(8)
                  Qweight(9)=cff*Qweight(9)
                ELSE
                  Qweight=0.0_r8               ! All the donor grid
                END IF                         ! corners are on land
                QWsum_check=SUM(Qweight)
#   ifdef DISTRIBUTE
                QW(1,m)=Qweight(1)
                QW(2,m)=Qweight(2)
                QW(3,m)=Qweight(3)
                QW(4,m)=Qweight(4)
                QW(5,m)=Qweight(5)
                QW(6,m)=Qweight(6)
                QW(7,m)=Qweight(7)
                QW(8,m)=Qweight(8)
                QW(9,m)=Qweight(9)
              ELSE
                QW(1,m)=Ucontact(cr)%Qweight(1,m)   ! we need to load
                QW(2,m)=Ucontact(cr)%Qweight(2,m)   ! unscaled values
                QW(3,m)=Ucontact(cr)%Qweight(3,m)   ! to facilitate
                QW(4,m)=Ucontact(cr)%Qweight(4,m)   ! parallel reduction
                QW(5,m)=Ucontact(cr)%Qweight(5,m)
                QW(6,m)=Ucontact(cr)%Qweight(6,m)
                QW(7,m)=Ucontact(cr)%Qweight(7,m)
                QW(8,m)=Ucontact(cr)%Qweight(8,m)
                QW(9,m)=Ucontact(cr)%Qweight(9,m)
#   else
                Ucontact(cr)%Qweight(1,m)=Qweight(1)
                Ucontact(cr)%Qweight(2,m)=Qweight(2)
                Ucontact(cr)%Qweight(3,m)=Qweight(3)
                Ucontact(cr)%Qweight(4,m)=Qweight(4)
                Ucontact(cr)%Qweight(5,m)=Qweight(5)
                Ucontact(cr)%Qweight(6,m)=Qweight(6)
                Ucontact(cr)%Qweight(7,m)=Qweight(7)
                Ucontact(cr)%Qweight(8,m)=Qweight(8)
                Ucontact(cr)%Qweight(9,m)=Qweight(9)
#   endif
              END IF
#  endif
            END IF
          END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
          CALL mp_assemble (dg, model, Lpoints, spv, LW)
          IF (exit_flag.ne.NoError) RETURN

#   ifdef QUADRATIC_WEIGHTS
          CALL mp_assemble (dg, model, Qpoints, spv, QW)
          IF (exit_flag.ne.NoError) RETURN
#   endif
!
!  Load exchanged weights.
!
          Ucontact(cr)%Lweight(1:4,1:NpointsU)=LW(1:4,1:NpointsU)
#   ifdef QUADRATIC_WEIGHTS
          Ucontact(cr)%Qweight(1:9,1:NpointsU)=QW(1:9,1:NpointsU)
#   endif
!
!  Deallocate local work arrays.
!
          IF (allocated(LW)) THEN
            deallocate (LW)
          END IF

#   ifdef QUADRATIC_WEIGHTS
          IF (allocated(QW)) THEN
            deallocate (QW)
          END IF
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Scale interpolation weigths for V-contact points.
!-----------------------------------------------------------------------
!
          NpointsV=Vcontact(cr)%Npoints

#  ifdef DISTRIBUTE
!
!  If distributed-memory, allocate and initialize working arrays
!  with special value (zero) to facilitate the global reduction
!  when collecting data between all nodes.
!
          Lpoints=4*NpointsV
          IF (.not.allocated(LW)) THEN
            allocate ( LW(4,NpointsV) )
          END IF
          LW=spv

#   ifdef QUADRATIC_WEIGHTS
          Lpoints=9*NpointsV
          IF (.not.allocated(QW)) THEN
            allocate ( QW(9,NpointsV) )
          END IF
          QW=spv
#   endif
#  endif
!
!  Scale interpolation weights.
!
          DO m=1,NpointsV
            Idg  =Vcontact(cr)%Idg(m)
            Idgp1=MIN(Idg+1, BOUNDS(dg)%UBi(-1))
            Jdg  =Vcontact(cr)%Jdg(m)
            Jdgp1=MIN(Jdg+1, BOUNDS(dg)%UBj(-1))
            IF (((Istr.le.Idg).and.(Idg.le.Iend)).and.                  &
     &          ((Jstr.le.Jdg).and.(Jdg.le.Jend))) THEN
!
!  Linear interpolation weights.
!
              MaskLWsum=GRID(dg)%vmask_full(Idg  ,Jdg  )+               &
     &                  GRID(dg)%vmask_full(Idgp1,Jdg  )+               &
     &                  GRID(dg)%vmask_full(Idgp1,Jdgp1)+               &
     &                  GRID(dg)%vmask_full(Idg  ,Jdgp1)
              IF (MaskLWsum.lt.4) THEN
#  ifdef WET_DRY
                Lweight(1)=Vcontact(cr)%LweightUnmasked(1,m)*           &
     &                     GRID(dg)%vmask_full(Idg  ,Jdg  )
                Lweight(2)=Vcontact(cr)%LweightUnmasked(2,m)*           &
     &                     GRID(dg)%vmask_full(Idgp1,Jdg  )
                Lweight(3)=Vcontact(cr)%LweightUnmasked(3,m)*           &
     &                     GRID(dg)%vmask_full(Idgp1,Jdgp1)
                Lweight(4)=Vcontact(cr)%LweightUnmasked(4,m)*           &
     &                     GRID(dg)%vmask_full(Idg  ,Jdgp1)
#  else
                Lweight(1)=Vcontact(cr)%Lweight(1,m)*                   &
     &                     GRID(dg)%vmask_full(Idg  ,Jdg  )
                Lweight(2)=Vcontact(cr)%Lweight(2,m)*                   &
     &                     GRID(dg)%vmask_full(Idgp1,Jdg  )
                Lweight(3)=Vcontact(cr)%Lweight(3,m)*                   &
     &                     GRID(dg)%vmask_full(Idgp1,Jdgp1)
                Lweight(4)=Vcontact(cr)%Lweight(4,m)*                   &
     &                     GRID(dg)%vmask_full(Idg  ,Jdgp1)
#  endif
                LWsum=SUM(Lweight)
                IF (LWsum.gt.0) THEN
                  cff=1.0_r8/LWsum
                  Lweight(1)=cff*Lweight(1)    ! using only water points
                  Lweight(2)=cff*Lweight(2)    ! in the interpolation of
                  Lweight(3)=cff*Lweight(3)    ! of contact points
                  Lweight(4)=cff*Lweight(4)
                ELSE
                  Lweight=0.0_r8               ! All the donor grid
                END IF                         ! corners are on land
                LWsum_check=SUM(Lweight)
#  ifdef DISTRIBUTE
                LW(1,m)=Lweight(1)
                LW(2,m)=Lweight(2)
                LW(3,m)=Lweight(3)
                LW(4,m)=Lweight(4)
              ELSE
                LW(1,m)=Vcontact(cr)%Lweight(1,m)   ! we need to load
                LW(2,m)=Vcontact(cr)%Lweight(2,m)   ! unscaled values
                LW(3,m)=Vcontact(cr)%Lweight(3,m)   ! to facilitate
                LW(4,m)=Vcontact(cr)%Lweight(4,m)   ! parallel reduction
#  else
                Vcontact(cr)%Lweight(1,m)=Lweight(1)
                Vcontact(cr)%Lweight(2,m)=Lweight(2)
                Vcontact(cr)%Lweight(3,m)=Lweight(3)
                Vcontact(cr)%Lweight(4,m)=Lweight(4)
#  endif
              END IF

#  ifdef QUADRATIC_WEIGHTS
!
!  Quadratic interpolation weights.
!
              MaskQWsum=GRID(dg)%vmask_full(Idg-1, Jdg-1)+              &
     &                  GRID(dg)%vmask_full(Idg  , Jdg-1)+              &
     &                  GRID(dg)%vmask_full(Idgp1, Jdg-1)+              &
     &                  GRID(dg)%vmask_full(Idg-1, Jdg  )+              &
     &                  GRID(dg)%vmask_full(Idg  , Jdg  )+              &
     &                  GRID(dg)%vmask_full(Idgp1, Jdg  )+              &
     &                  GRID(dg)%vmask_full(Idg-1, Jdgp1)+              &
     &                  GRID(dg)%vmask_full(Idg  , Jdgp1)+              &
     &                  GRID(dg)%vmask_full(Idgp1, Jdgp1)
              IF (MaskQWsum.lt.9) THEN
#   ifdef WET_DRY
                Qweight(1)=Vcontact(cr)%QweightUnmasked(1,m)*           &
     &                     GRID(dg)%vmask_full(Idg-1,Jdg-1)
                Qweight(2)=Vcontact(cr)%QweightUnmasked(2,m)*           &
     &                     GRID(dg)%vmask_full(Idg  ,Jdg-1)
                Qweight(3)=Vcontact(cr)%QweightUnmasked(3,m)*           &
     &                     GRID(dg)%vmask_full(Idgp1,Jdg-1)
                Qweight(4)=Vcontact(cr)%QweightUnmasked(4,m)*           &
     &                     GRID(dg)%vmask_full(Idg-1,Jdg  )
                Qweight(5)=Vcontact(cr)%QweightUnmasked(5,m)*           &
     &                     GRID(dg)%vmask_full(Idg  ,Jdg  )
                Qweight(6)=Vcontact(cr)%QweightUnmasked(6,m)*           &
     &                     GRID(dg)%vmask_full(Idgp1,Jdg  )
                Qweight(7)=Vcontact(cr)%QweightUnmasked(7,m)*           &
     &                     GRID(dg)%vmask_full(Idg-1,Jdgp1)
                Qweight(8)=Vcontact(cr)%QweightUnmasked(8,m)*           &
     &                     GRID(dg)%vmask_full(Idg  ,Jdgp1)
                Qweight(9)=Vcontact(cr)%QweightUnmasked(9,m)*           &
     &                     GRID(dg)%vmask_full(Idgp1,Jdgp1)
#   else
                Qweight(1)=Vcontact(cr)%Qweight(1,m)*                   &
     &                     GRID(dg)%vmask_full(Idg-1,Jdg-1)
                Qweight(2)=Vcontact(cr)%Qweight(2,m)*                   &
     &                     GRID(dg)%vmask_full(Idg  ,Jdg-1)
                Qweight(3)=Vcontact(cr)%Qweight(3,m)*                   &
     &                     GRID(dg)%vmask_full(Idgp1,Jdg-1)
                Qweight(4)=Vcontact(cr)%Qweight(4,m)*                   &
     &                     GRID(dg)%vmask_full(Idg-1,Jdg  )
                Qweight(5)=Vcontact(cr)%Qweight(5,m)*                   &
     &                     GRID(dg)%vmask_full(Idg  ,Jdg  )
                Qweight(6)=Vcontact(cr)%Qweight(6,m)*                   &
     &                     GRID(dg)%vmask_full(Idgp1,Jdg  )
                Qweight(7)=Vcontact(cr)%Qweight(7,m)*                   &
     &                     GRID(dg)%vmask_full(Idg-1,Jdgp1)
                Qweight(8)=Vcontact(cr)%Qweight(8,m)*                   &
     &                     GRID(dg)%vmask_full(Idg  ,Jdgp1)
                Qweight(9)=Vcontact(cr)%Qweight(9,m)*                   &
     &                     GRID(dg)%vmask_full(Idgp1,Jdgp1)
#   endif
                QWsum=SUM(Qweight)
                IF (QWsum.gt.0) THEN
                  cff=1.0_r8/QWsum
                  Qweight(1)=cff*Qweight(1)    ! using only water points
                  Qweight(2)=cff*Qweight(2)    ! in the interpolation of
                  Qweight(3)=cff*Qweight(3)    ! of contact points
                  Qweight(4)=cff*Qweight(4)
                  Qweight(5)=cff*Qweight(5)
                  Qweight(6)=cff*Qweight(6)
                  Qweight(7)=cff*Qweight(7)
                  Qweight(8)=cff*Qweight(8)
                  Qweight(9)=cff*Qweight(9)
                ELSE
                  Qweight=0.0_r8               ! All the donor grid
                END IF                         ! corners are on land
                END IF
                QWsum_check=SUM(Qweight)
#   ifdef DISTRIBUTE
                QW(1,m)=Qweight(1)
                QW(2,m)=Qweight(2)
                QW(3,m)=Qweight(3)
                QW(4,m)=Qweight(4)
                QW(5,m)=Qweight(5)
                QW(6,m)=Qweight(6)
                QW(7,m)=Qweight(7)
                QW(8,m)=Qweight(8)
                QW(9,m)=Qweight(9)
              ELSE
                QW(1,m)=Vcontact(cr)%Qweight(1,m)   ! we need to load
                QW(2,m)=Vcontact(cr)%Qweight(2,m)   ! unscaled values
                QW(3,m)=Vcontact(cr)%Qweight(3,m)   ! to facilitate
                QW(4,m)=Vcontact(cr)%Qweight(4,m)   ! parallel reduction
                QW(5,m)=Vcontact(cr)%Qweight(5,m)
                QW(6,m)=Vcontact(cr)%Qweight(6,m)
                QW(7,m)=Vcontact(cr)%Qweight(7,m)
                QW(8,m)=Vcontact(cr)%Qweight(8,m)
                QW(9,m)=Vcontact(cr)%Qweight(9,m)
#   else
                Vcontact(cr)%Qweight(1,m)=Qweight(1)
                Vcontact(cr)%Qweight(2,m)=Qweight(2)
                Vcontact(cr)%Qweight(3,m)=Qweight(3)
                Vcontact(cr)%Qweight(4,m)=Qweight(4)
                Vcontact(cr)%Qweight(5,m)=Qweight(5)
                Vcontact(cr)%Qweight(6,m)=Qweight(6)
                Vcontact(cr)%Qweight(7,m)=Qweight(7)
                Vcontact(cr)%Qweight(8,m)=Qweight(8)
                Vcontact(cr)%Qweight(9,m)=Qweight(9)
#   endif
              END IF
#  endif
            END IF
          END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
          CALL mp_assemble (dg, model, Lpoints, spv, LW)
          IF (exit_flag.ne.NoError) RETURN

#   ifdef QUADRATIC_WEIGHTS
          CALL mp_assemble (dg, model, Qpoints, spv, QW)
          IF (exit_flag.ne.NoError) RETURN
#   endif
!
!  Load exchanged weights.
!
          Vcontact(cr)%Lweight(1:4,1:NpointsV)=LW(1:4,1:NpointsV)
#   ifdef QUADRATIC_WEIGHTS
          Vcontact(cr)%Qweight(1:9,1:NpointsV)=QW(1:9,1:NpointsV)
#   endif
!
!  Deallocate local work arrays.
!
          IF (allocated(LW)) THEN
            deallocate (LW)
          END IF

#   ifdef QUADRATIC_WEIGHTS
          IF (allocated(QW)) THEN
            deallocate (QW)
          END IF
#   endif
#  endif

        END IF SCALE_WEIGHTS
      END DO

      RETURN
      END SUBROUTINE mask_hweights
# endif
!
      SUBROUTINE get_composite (ng, model, isection, tile)
!
!=======================================================================
!                                                                      !
!  This routine gets the donor grid data required to process the       !
!  contact points of the current composite grid. It extracts the       !
!  donor cell points containing each contact point. In composite       !
!  grids, it is possible to have more than one contact region.         !
!                                                                      !
!  The interpolation of composite grid contact points from donor       !
!  grid data is carried out in a different parallel region using       !
!  'put_composite'.                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Composite grid number (integer)                       !
!     model      Calling model identifier (integer)                    !
!     isection   Governing equations time-stepping section in          !
!                  main2d or main3d indicating which state             !
!                  variables to process (integer)                      !
!     tile       Domain tile partition (integer)                       !
!                                                                      !
!  On Output:    (mod_nesting)                                         !
!                                                                      !
!     COMPOSITE  Updated contact points structure.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_nesting
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, isection, tile
!
!  Local variable declarations.
!
      integer :: cr, dg, rg, nrec, rec
# ifdef SOLVE3D
      integer :: itrc
# endif
      integer :: LBi, UBi, LBj, UBj
      integer :: Tindex
!
!-----------------------------------------------------------------------
!  Get donor grid data needed to process composite grid contact points.
!  Only process those variables associated with the governing equation
!  time-stepping section.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
!
!  Get data donor and data receiver grid numbers.
!
        dg=Rcontact(cr)%donor_grid
        rg=Rcontact(cr)%receiver_grid
!
!  Process only contact region data for requested nested grid "ng".
!
        IF (rg.eq.ng) THEN
!
!  Set donor grid lower and upper array indices.
!
          LBi=BOUNDS(dg)%LBi(tile)
          UBi=BOUNDS(dg)%UBi(tile)
          LBj=BOUNDS(dg)%LBj(tile)
          UBj=BOUNDS(dg)%UBj(tile)
!
!  Process bottom stress (bustr, bvstr).
!
          IF (isection.eq.nbstr) THEN
            CALL get_contact2d (dg, model, tile,                        &
     &                          u2dvar, Vname(1,idUbms),                &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(dg) % bustr,                     &
     &                          COMPOSITE(cr) % bustr)
            CALL get_contact2d (dg, model, tile,                        &
     &                          v2dvar, Vname(1,idVbms),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(dg) % bvstr,                     &
     &                          COMPOSITE(cr) % bvstr)
          END IF
!
!  Process free-surface (zeta) at the appropriate time index.
!
          IF ((isection.eq.nFSIC).or.                                   &
     &        (isection.eq.nzeta).or.                                   &
     &        (isection.eq.n2dPS).or.                                   &
     &        (isection.eq.n2dCS)) THEN
            IF (isection.eq.nzeta) THEN
              nrec=2                   ! process time records 1 and 2
            ELSE
              nrec=1                   ! process knew record
            END IF
            DO rec=1,nrec
              IF (isection.eq.nzeta) THEN
                Tindex=rec
              ELSE
                Tindex=knew(dg)
              END IF
              CALL get_contact2d (dg, model, tile,                      &
     &                            r2dvar, Vname(1,idFsur),              &
     &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OCEAN(dg) % zeta(:,:,Tindex),         &
     &                            COMPOSITE(cr) % zeta(:,:,rec))
            END DO
          END IF
!
!  Process free-surface equation rigth-hand-side (rzeta) term.
!
          IF (isection.eq.n2dPS) THEN
            Tindex=1
            CALL get_contact2d (dg, model, tile,                        &
     &                          r2dvar, Vname(1,idRzet),                &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(dg) % rzeta(:,:,Tindex),          &
     &                          COMPOSITE(cr) % rzeta)
          END IF
!
!  Process 2D momentum components (ubar,vbar) at the appropriate time
!  index.
!
          IF ((isection.eq.n2dIC).or.                                   &
     &        (isection.eq.n2dPS).or.                                   &
     &        (isection.eq.n2dCS).or.                                   &
     &        (isection.eq.n3duv)) THEN
            IF (isection.eq.n3duv) THEN
              nrec=2                   ! process time records 1 and 2
            ELSE
              nrec=1                   ! process knew record
            END IF
            DO rec=1,nrec
              IF (isection.eq.n3duv) THEN
                Tindex=rec
              ELSE
                Tindex=knew(dg)
              END IF
              CALL get_contact2d (dg, model, tile,                      &
     &                            u2dvar, Vname(1,idUbar),              &
     &                            cr, Ucontact(cr)%Npoints, Ucontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OCEAN(dg) % ubar(:,:,Tindex),         &
     &                            COMPOSITE(cr) % ubar(:,:,rec))
              CALL get_contact2d (dg, model, tile,                      &
     &                            v2dvar, Vname(1,idVbar),              &
     &                            cr, Vcontact(cr)%Npoints, Vcontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OCEAN(dg) % vbar(:,:,Tindex),         &
     &                            COMPOSITE(cr) % vbar(:,:,rec))
            END DO
          END IF

# ifdef SOLVE3D
!
!  Process time averaged free-surface (Zt_avg1) and 2D momentum fluxes
!  (DU_avg1, DV_avg1).
!
          IF (isection.eq.n2dfx) THEN
            CALL get_contact2d (dg, model, tile,                        &
     &                          r2dvar, 'Zt_avg1',                      &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % Zt_avg1,                 &
     &                          COMPOSITE(cr) % Zt_avg1)
            CALL get_contact2d (dg, model, tile,                        &
     &                          u2dvar, 'DU_avg1',                      &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % DU_avg1,                 &
     &                          COMPOSITE(cr) % DU_avg1)
            CALL get_contact2d (dg, model, tile,                        &
     &                          v2dvar, 'DV_avg1',                      &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % DV_avg1,                 &
     &                          COMPOSITE(cr) % DV_avg1)
          END IF

#  if !defined TS_FIXED
!
!  Process tracer variables (t) at the appropriate time index.
!
          IF ((isection.eq.nTVIC).or.                                   &
     &        (isection.eq.nrhst).or.                                   &
     &        (isection.eq.n3dTV)) THEN
            DO itrc=1,NT(ng)
              IF (isection.eq.nrhst) THEN
                Tindex=3
              ELSE
                Tindex=nnew(dg)
              END IF
              CALL get_contact3d (dg, model, tile,                      &
     &                            r3dvar, Vname(1,idTvar(itrc)),        &
     &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
     &                            LBi, UBi, LBj, UBj, 1, N(dg),         &
     &                            OCEAN(dg) % t(:,:,:,Tindex,itrc),     &
     &                            COMPOSITE(cr) % t(:,:,:,itrc))
            END DO
          END IF
#  endif
!
!  Process 3D momentum (u, v) at the appropriate time-index.
!
          IF ((isection.eq.n3dIC).or.                                   &
     &        (isection.eq.n3duv)) THEN
            Tindex=nnew(dg)
            CALL get_contact3d (dg, model, tile,                        &
     &                          u3dvar, Vname(1,idUvel),                &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          OCEAN(dg) % u(:,:,:,Tindex),            &
     &                          COMPOSITE(cr) % u)
            CALL get_contact3d (dg, model, tile,                        &
     &                          v3dvar, Vname(1,idVvel),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          OCEAN(dg) % v(:,:,:,Tindex),            &
     &                          COMPOSITE(cr) % v)
          END IF
!
!  Process 3D momentum fluxes (Huon, Hvom).
!
          IF (isection.eq.n3duv) THEN
            CALL get_contact3d (dg, model, tile,                        &
     &                          u3dvar, 'Huon',                         &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          GRID(dg) % Huon,                        &
     &                          COMPOSITE(cr) % Huon)
            CALL get_contact3d (dg, model, tile,                        &
     &                          v3dvar, 'Hvom',                         &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          GRID(dg) % Hvom,                        &
     &                          COMPOSITE(cr) % Hvom)
          END IF
# endif

        END IF
      END DO

      RETURN
      END SUBROUTINE get_composite
!
      SUBROUTINE get_metrics (ng, model, tile)
!
!=======================================================================
!                                                                      !
!  This routine extracts grid spacing metrics "on_u" and "om_v",       !
!  which are used to impose mass flux at the finer grid physical       !
!  boundaries in refinement applications.                              !
!                                                                      !
!  The extracted metrics is stored in REFINED structure and are        !
!  needed in 'put_refine2d'.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Refinement grid number (integer)                      !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!                                                                      !
!  On Output:    (mod_nesting)                                         !
!                                                                      !
!     REFINED    Updated contact points structure.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_nesting
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
!
!  Local variable declarations.
!
      integer :: cr, dg, i, j, m
      integer :: Istr, Iend, Jstr, Jend
      integer :: NpointsU, NpointsV

      real(r8), parameter :: spv = 0.0_r8
!
!-----------------------------------------------------------------------
!  Extract grid spacing metrics.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
!
!  Get data donor grid number.
!
        dg=Rcontact(cr)%donor_grid
!
!  Extract grid spacing at U-points.
!
        Istr=BOUNDS(dg) % IstrP(tile)
        Iend=BOUNDS(dg) % IendT(tile)
        Jstr=BOUNDS(dg) % JstrT(tile)
        Jend=BOUNDS(dg) % JendT(tile)
!
        NpointsU=Ucontact(cr) % Npoints
        REFINED(cr) % on_u(1:NpointsU) = spv
!
        DO m=1,NpointsU
          i=Ucontact(cr) % Idg(m)
          j=Ucontact(cr) % Jdg(m)
          IF (((Istr.le.i).and.(i.le.Iend)).and.                        &
     &        ((Jstr.le.j).and.(j.le.Jend))) THEN
            REFINED(cr) % on_u(m) = GRID(dg) % on_u(i,j)
          END IF
        END DO
!
!  Extract grid spacing at V-points.
!
        Istr=BOUNDS(dg) % IstrT(tile)
        Iend=BOUNDS(dg) % IendT(tile)
        Jstr=BOUNDS(dg) % JstrP(tile)
        Jend=BOUNDS(dg) % JendT(tile)
!
        NpointsV=Vcontact(cr) % Npoints
        REFINED(cr) % om_v(1:NpointsV) = spv
!
        DO m=1,NpointsV
          i=Vcontact(cr) % Idg(m)
          j=Vcontact(cr) % Jdg(m)
          IF (((Istr.le.i).and.(i.le.Iend)).and.                        &
     &        ((Jstr.le.j).and.(j.le.Jend))) THEN
            REFINED(cr) % om_v(m) = GRID(dg) % om_v(i,j)
          END IF
        END DO

# ifdef DISTRIBUTE
!
!  Exchange data between all nodes.
!
        CALL mp_assemble (dg, model, NpointsU, spv, REFINED(cr) % on_u)
        IF (exit_flag.ne.NoError) RETURN

        CALL mp_assemble (dg, model, NpointsV, spv, REFINED(cr) % om_v)
        IF (exit_flag.ne.NoError) RETURN
#  endif

      END DO

      RETURN
      END SUBROUTINE get_metrics
!
      SUBROUTINE get_refine (ng, model, tile)
!
!=======================================================================
!                                                                      !
!  This routine gets the donor grid data required to process the       !
!  contact points of the  current  refinement  grid. It extracts       !
!  the donor cell points containing each contact point.                !
!                                                                      !
!  The extracted data is stored in two-time rolling records which      !
!  are needed for the space and time interpolation in 'put_refine'.    !
!                                                                      !
!  Except for initialization, this routine is called at the bottom     !
!  of the donor grid time step so all the values are updated for the   !
!  time(dg)+dt(dg). That is, in 2D applications it is called after     !
!  "step2d" corrector step and in 3D applications it is called after   !
!  "step3d_t". This is done to have the coarser grid snapshots at      !
!  time(dg) and time(dg)+dt(dg) to bound the interpolation of the      !
!  finer grid contact points.                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Refinement grid number (integer)                      !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!                                                                      !
!  On Output:    (mod_nesting)                                         !
!                                                                      !
!     REFINED    Updated contact points structure.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
      USE mod_ncparam
      USE mod_nesting
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
!
!  Local variable declarations.
!
# ifdef NESTING_DEBUG
      logical, save :: first = .TRUE.
# endif
      integer :: Tindex2d, cr, dg, ir, rg, tnew
# ifdef SOLVE3D
      integer :: Tindex3d, itrc
# endif
      integer :: LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Get donor grid data needed to process refinement grid contact points.
!  The extracted contact point data is stored in two time records to
!  facilitate the space-time interpolation elsewhere.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
!
!  Get data donor and data receiver grid numbers.
!
        dg=Rcontact(cr)%donor_grid
        rg=Rcontact(cr)%receiver_grid
!
!  Process only contact region data for requested nested grid "ng".
!
        IF ((dg.eq.CoarserDonor(rg)).and.(dg.eq.ng)) THEN
!
!  Set donor grid lower and upper array indices.
!
          LBi=BOUNDS(dg)%LBi(tile)
          UBi=BOUNDS(dg)%UBi(tile)
          LBj=BOUNDS(dg)%LBj(tile)
          UBj=BOUNDS(dg)%UBj(tile)
!
!  Update rolling time indices. The contact data is stored in two time
!  levels.  We need a special case for ROMS initialization in "main2d"
!  or "main3d" after the processing "ini_fields".  Notice that a dt(dg)
!  is added because this routine is called after the end of the time
!  step.
!
          IF (RollingIndex(cr).eq.0) THEN
            tnew=1                                  ! ROMS
            RollingIndex(cr)=tnew                   ! initialization
            RollingTime(tnew,cr)=time(dg)           ! step
          ELSE
            tnew=3-RollingIndex(cr)
            RollingIndex(cr)=tnew
            RollingTime(tnew,cr)=time(dg)+dt(dg)
          END IF
!
!  Set donor grid time index to process. In 3D applications, the 2D
!  record index to use can be either 1 or 2 since both ubar(:,:,1:2)
!  and vbar(:,:,1:2) are set to its time-averaged values in "step3d_uv".
!  That is, we can use Tindex2d=kstp(dg) or Tindex2d=knew(dg). However,
!  in 2D applications we need to use Tindex2d=knew(dg).
!
          Tindex2d=knew(dg)
# ifdef SOLVE3D
          Tindex3d=nnew(dg)
# endif

# ifdef NESTING_DEBUG
!
!  If debugging, write information into Fortran unit 100 to check the
!  logic of processing donor grid data.
!
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              IF (first) THEN
                first=.FALSE.
                WRITE (100,10)
              END IF
              WRITE (100,20) ng, cr, dg, rg, iic(dg), iic(rg),          &
     &                       3-tnew, tnew, Tindex2d, Tindex3d,          &
     &                       INT(RollingTime(3-tnew,cr)),               &
     &                       INT(RollingTime(tnew,cr)),                 &
     &                       INT(time(ng))
 10           FORMAT (2x,'ng',2x,'cr',2x,'dg',2x,'rg',5x,'iic',4x,'iic',&
     &                3x,'told',3x,'tnew',2x,'Tindex',1x,'Tindex',      &
     &                2x,'time',3x,'time',3x,'time',/,20x,'(dg)',       &
     &                3x,'(rg)',18x,'2D',5x,'3D',4x,'told',3x,'tnew',   &
     &                3x,'(ng)',/)
 20           FORMAT (4i4,9i7)
              CALL my_flush (100)
            END IF
          END IF
# endif
!
!  Extract free-surface.
!
          CALL get_contact2d (dg, model, tile,                          &
     &                        r2dvar, 'Zt_avg1',                        &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        COUPLING(dg) % Zt_avg1,                   &
     &                        REFINED(cr) % zeta(:,:,tnew))
!
!  Extract 2D momentum components (ubar, vbar).
!
          CALL get_contact2d (dg, model, tile,                          &
     &                        u2dvar, Vname(1,idUbar),                  &
     &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(dg) % ubar(:,:,Tindex2d),           &
     &                        REFINED(cr) % ubar(:,:,tnew))

          CALL get_contact2d (dg, model, tile,                          &
     &                        v2dvar, Vname(1,idVbar),                  &
     &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(dg) % vbar(:,:,Tindex2d),           &
     &                        REFINED(cr) % vbar(:,:,tnew))

# ifdef SOLVE3D
!
!  Extract time-veraged fluxes (DU_avg2, DV_avg2).  We will use latter
!  only the values at the finer grid physical boundary to impose mass
!  flux conservation in routine "put_refine2d".
!
          CALL get_persisted2d (dg, rg, model, tile,                    &
     &                          u2dvar, 'DU_avg2',                      &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % DU_avg2,                 &
     &                          REFINED(cr) % DU_avg2(:,:,tnew))

          CALL get_persisted2d (dg, rg, model, tile,                    &
     &                          v2dvar, 'DV_avg2',                      &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % DV_avg2,                 &
     &                          REFINED(cr) % DV_avg2(:,:,tnew))
!
!  Tracer-type variables.
!
          DO itrc=1,NT(dg)
            CALL get_contact3d (dg, model, tile,                        &
     &                          r3dvar, Vname(1,idTvar(itrc)),          &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          OCEAN(dg) % t(:,:,:,Tindex3d,itrc),     &
     &                          REFINED(cr) % t(:,:,:,tnew,itrc))
          END DO
!
!  Extract 3D momentum components (u, v).
!
          CALL get_contact3d (dg, model, tile,                          &
     &                        u3dvar, Vname(1,idUvel),                  &
     &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
     &                        LBi, UBi, LBj, UBj, 1, N(dg),             &
     &                        OCEAN(dg) % u(:,:,:,Tindex3d),            &
     &                        REFINED(cr) % u(:,:,:,tnew))

          CALL get_contact3d (dg, model, tile,                          &
     &                        v3dvar, Vname(1,idVvel),                  &
     &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
     &                        LBi, UBi, LBj, UBj, 1, N(dg),             &
     &                        OCEAN(dg) % v(:,:,:,Tindex3d),            &
     &                        REFINED(cr) % v(:,:,:,tnew))
# endif
        END IF
      END DO

      RETURN
      END SUBROUTINE get_refine
!
      SUBROUTINE put_composite (ng, model, isection, tile)
!
!=======================================================================
!                                                                      !
!  This routine interpolates composite grid contact points from donor  !
!  grid data extracted in routine 'get_composite'.                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Composite grid number (integer)                       !
!     model      Calling model identifier (integer)                    !
!     isection   Governing equations time-stepping section in          !
!                  main2d or main3d indicating which state             !
!                  variables to process (integer)                      !
!     tile       Domain tile partition (integer)                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_nesting
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping

# ifdef DISTRIBUTE
!
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#  endif
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, isection, tile
!
!  Local variable declarations.
!
      integer :: dg, rg, cr, nrec, rec
# ifdef SOLVE3D
      integer :: itrc
# endif
      integer :: LBi, UBi, LBj, UBj
      integer :: Tindex
!
!-----------------------------------------------------------------------
!  Interpolate composite grid contact points from donor grid data.
!  Only process those variables associated with the governing equation
!  time-stepping section.
!-----------------------------------------------------------------------
!
      CR_LOOP : DO cr=1,Ncontact
!
!  Get data donor and data receiver grid numbers.
!
        dg=Rcontact(cr)%donor_grid
        rg=Rcontact(cr)%receiver_grid
!
!  Process only contact region data for requested nested grid "ng".
!
        IF (rg.eq.ng) THEN
!
!  Set receiver grid lower and upper array indices.
!
          LBi=BOUNDS(rg)%LBi(tile)
          UBi=BOUNDS(rg)%UBi(tile)
          LBj=BOUNDS(rg)%LBj(tile)
          UBj=BOUNDS(rg)%UBj(tile)
!
!  Process bottom stress (bustr, bvstr).
!
          IF (isection.eq.nbstr) THEN
            CALL put_contact2d (rg, model, tile,                        &
     &                          u2dvar, Vname(1,idUbms),                &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
# ifdef MASKING
     &                          GRID(rg) % umask,                       &
# endif
     &                          COMPOSITE(cr) % bustr,                  &
     &                          FORCES(rg) % bustr)
            CALL put_contact2d (rg, model, tile,                        &
     &                          v2dvar, Vname(1,idVbms),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
# ifdef MASKING
     &                          GRID(rg) % vmask,                       &
# endif
     &                          COMPOSITE(cr) % bvstr,                  &
     &                          FORCES(rg) % bvstr)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (rg, tile, model, 2,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          FORCES(rg) % bustr,                     &
     &                          FORCES(rg) % bvstr)
# endif
          END IF
!
!  Process free-surface (zeta) at the appropriate time index.
!
          IF ((isection.eq.nFSIC).or.                                   &
     &        (isection.eq.nzeta).or.                                   &
     &        (isection.eq.n2dPS).or.                                   &
     &        (isection.eq.n2dCS)) THEN
            IF (isection.eq.nzeta) THEN
              nrec=2                   ! process time records 1 and 2
            ELSE
              nrec=1                   ! process knew record
            END IF
            DO rec=1,nrec
              IF (isection.eq.nzeta) THEN
                Tindex=rec
              ELSE
                Tindex=knew(rg)
              END IF
              CALL put_contact2d (rg, model, tile,                      &
     &                            r2dvar, Vname(1,idFsur),              &
     &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
# ifdef MASKING
     &                            GRID(rg) % rmask,                     &
# endif
     &                            COMPOSITE(cr) % zeta(:,:,rec),        &
     &                            OCEAN(rg) % zeta(:,:,Tindex))
# ifdef DISTRIBUTE
              CALL mp_exchange2d (rg, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(rg), NSperiodic(rg),       &
     &                            OCEAN(rg) % zeta(:,:,Tindex))
# endif
            END DO
          END IF
!
!  Process free-surface equation rigth-hand-side (rzeta) term.
!
          IF (isection.eq.n2dPS) THEN
            Tindex=1
            CALL put_contact2d (rg, model, tile,                        &
     &                          r2dvar, Vname(1,idRzet),                &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
# ifdef MASKING
     &                          GRID(rg) % rmask,                       &
# endif
     &                          COMPOSITE(cr) % rzeta,                  &
     &                          OCEAN(rg) % rzeta(:,:,Tindex))
# ifdef DISTRIBUTE
            CALL mp_exchange2d (rg, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg) % rzeta(:,:,Tindex))
# endif
          END IF
!
!  Process 2D momentum components (ubar,vbar) at the appropriate time
!  index.
!
          IF ((isection.eq.n2dIC).or.                                   &
     &        (isection.eq.n2dPS).or.                                   &
     &        (isection.eq.n2dCS).or.                                   &
     &        (isection.eq.n3duv)) THEN
            IF (isection.eq.n3duv) THEN
              nrec=2                   ! process time records 1 and 2
            ELSE
              nrec=1                   ! process knew record
            END IF
            DO rec=1,nrec
              IF (isection.eq.n3duv) THEN
                Tindex=rec
              ELSE
                Tindex=knew(rg)
              END IF
              CALL put_contact2d (rg, model, tile,                      &
     &                            u2dvar, Vname(1,idUbar),              &
     &                            cr, Ucontact(cr)%Npoints, Ucontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
# ifdef MASKING
     &                            GRID(rg) % umask,                     &
# endif
     &                            COMPOSITE(cr) % ubar(:,:,rec),        &
     &                            OCEAN(rg) % ubar(:,:,Tindex))
              CALL put_contact2d (rg, model, tile,                      &
     &                            v2dvar, Vname(1,idVbar),              &
     &                            cr, Vcontact(cr)%Npoints, Vcontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
# ifdef MASKING
     &                            GRID(rg) % vmask,                     &
# endif
     &                            COMPOSITE(cr) % vbar(:,:,rec),        &
     &                            OCEAN(rg) % vbar(:,:,Tindex))
# ifdef DISTRIBUTE
              CALL mp_exchange2d (rg, tile, model, 2,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(rg), NSperiodic(rg),       &
     &                            OCEAN(rg) % ubar(:,:,Tindex),         &
     &                            OCEAN(rg) % vbar(:,:,Tindex))
# endif
            END DO
          END IF

# ifdef SOLVE3D
!
!  Process time averaged free-surface (Zt_avg1) and 2D momentum fluxes
!  (DU_avg1, DV_avg1).
!
          IF (isection.eq.n2dfx) THEN
            CALL put_contact2d (rg, model, tile,                        &
     &                          r2dvar, 'Zt_avg1',                      &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
#  ifdef MASKING
     &                          GRID(rg) % rmask,                       &
#  endif
     &                          COMPOSITE(cr) % Zt_avg1,                &
     &                          COUPLING(rg) % Zt_avg1)
            CALL put_contact2d (rg, model, tile,                        &
     &                          u2dvar, Vname(1,idUfx1),                &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
#  ifdef MASKING
     &                          GRID(rg) % umask,                       &
#  endif
     &                          COMPOSITE(cr) % DU_avg1,                &
     &                          COUPLING(rg) % DU_avg1)
            CALL put_contact2d (rg, model, tile,                        &
     &                          v2dvar, Vname(1,idVfx1),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
#  ifdef MASKING
     &                          GRID(rg) % vmask,                       &
#  endif
     &                          COMPOSITE(cr) % DV_avg1,                &
     &                          COUPLING(rg) % DV_avg1)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (rg, tile, model, 3,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          COUPLING(rg) % Zt_avg1,                 &
     &                          COUPLING(rg) % DU_avg1,                 &
     &                          COUPLING(rg) % DV_avg1)
# endif
          END IF

#  if !defined TS_FIXED
!
!  Process tracer variables (t) at the appropriate time index.
!
          IF ((isection.eq.nTVIC).or.                                   &
     &        (isection.eq.nrhst).or.                                   &
     &        (isection.eq.n3dTV)) THEN
            DO itrc=1,NT(ng)
              IF (isection.eq.nrhst) THEN
                Tindex=3
              ELSE
                Tindex=nnew(rg)
              END IF
              CALL put_contact3d (rg, model, tile,                      &
     &                            r3dvar, Vname(1,idTvar(itrc)),        &
     &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
     &                            LBi, UBi, LBj, UBj, 1, N(rg),         &
#   ifdef MASKING
     &                            GRID(rg) % rmask,                     &
#   endif
     &                            COMPOSITE(cr) % t(:,:,:,itrc),        &
     &                            OCEAN(rg) % t(:,:,:,Tindex,itrc))
            END DO
#   ifdef DISTRIBUTE
            CALL mp_exchange4d (rg, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg), 1, NT(rg),&
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg) % t(:,:,:,Tindex,:))
#   endif
          END IF
#  endif
!
!  Process 3D momentum (u, v) at the appropriate time-index.
!
          IF ((isection.eq.n3dIC).or.                                   &
     &        (isection.eq.n3duv)) THEN
            Tindex=nnew(rg)
            CALL put_contact3d (rg, model, tile,                        &
     &                          u3dvar, Vname(1,idUvel),                &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
#  ifdef MASKING
     &                          GRID(rg) % umask,                       &
#  endif
     &                          COMPOSITE(cr) % u,                      &
     &                          OCEAN(rg) % u(:,:,:,Tindex))
            CALL put_contact3d (rg, model, tile,                        &
     &                          v3dvar, Vname(1,idVvel),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
#  ifdef MASKING
     &                          GRID(rg) % vmask,                       &
#  endif
     &                          COMPOSITE(cr) % v,                      &
     &                          OCEAN(rg) % v(:,:,:,Tindex))
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (rg, tile, model, 2,                     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg) % u(:,:,:,Tindex),            &
     &                          OCEAN(rg) % v(:,:,:,Tindex))
#  endif
          END IF
!
!  Process 3D momentum fluxes (Huon, Hvom).
!
          IF (isection.eq.n3duv) THEN
            CALL put_contact3d (rg, model, tile,                        &
     &                          u3dvar, 'Huon',                         &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
#  ifdef MASKING
     &                          GRID(rg) % umask,                       &
#  endif
     &                          COMPOSITE(cr) % Huon,                   &
     &                          GRID(rg) % Huon)
            CALL put_contact3d (rg, model, tile,                        &
     &                          v3dvar, 'Hvom',                         &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
#  ifdef MASKING
     &                          GRID(rg) % vmask,                       &
#  endif
     &                          COMPOSITE(cr) % Hvom,                   &
     &                          GRID(rg) % Hvom)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (rg, tile, model, 2,                     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          GRID(rg) % Huon,                        &
     &                          GRID(rg) % Hvom)
#  endif
          END IF
# endif

        END IF
      END DO CR_LOOP

      RETURN
      END SUBROUTINE put_composite
!
      SUBROUTINE put_refine (ng, model, tile, LputFsur)
!
!=======================================================================
!                                                                      !
!  This routine interpolates refinement grid contact points from donor !
!  grid data extracted in routine 'get_refine'. Notice that because of !
!  shared-memory parallelism,  the free-surface is processed first and !
!  in a different parallel region.
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Refinement grid number (integer)                      !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     LputFsur   Switch to process or not free-surface (logical)       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_nesting
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
!  Imported variable declarations.
!
      logical, intent(in) :: LputFsur
      integer, intent(in) :: ng, model, tile
!
!  Local variable declarations.
!
      integer :: dg, rg, cr, nrec, rec
# ifdef SOLVE3D
      integer :: itrc
# endif
      integer :: LBi, UBi, LBj, UBj
      integer :: Tindex
!
!-----------------------------------------------------------------------
!  Interpolate refinement grid contact points from donor grid data
!  (space-time interpolation)
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
!
!  Get data donor and data receiver grid numbers.
!
        dg=Rcontact(cr)%donor_grid
        rg=Rcontact(cr)%receiver_grid
!
!  Process only contact region data for requested nested grid "ng", if
!  donor grid is coarser than receiver grid.  That is, we are only
!  processing external contact points areas.
!
        IF ((rg.eq.ng).and.(DXmax(dg).gt.DXmax(rg))) THEN
!
!  Set receiver grid lower and upper array indices.
!
          LBi=BOUNDS(rg)%LBi(tile)
          UBi=BOUNDS(rg)%UBi(tile)
          LBj=BOUNDS(rg)%LBj(tile)
          UBj=BOUNDS(rg)%UBj(tile)
!
!  Fill free-surface separatelly.
!
          IF (LputFsur) THEN
            CALL put_refine2d (ng, dg, cr, model, tile, LputFsur,       &
     &                         LBi, UBi, LBj, UBj)
          ELSE
!
!  Fill other 2D state variables (like momentum) contact points.
!
            CALL put_refine2d (ng, dg, cr, model, tile, LputFsur,       &
     &                         LBi, UBi, LBj, UBj)

# ifdef SOLVE3D
!
!  Fill 3D state variables contact points.
!
            CALL put_refine3d (ng, dg, cr, model, tile,                 &
     &                         LBi, UBi, LBj, UBj)
# endif
          END IF
        END IF
      END DO

      RETURN
      END SUBROUTINE put_refine

# ifdef SOLVE3D
!
      SUBROUTINE correct_tracer (ng, ngf, model, tile)
!
!=======================================================================
!                                                                      !
!  This routine corrects the tracer values in the coarser grid at the  !
!  location of the finer grid physical domain perimeter by comparing   !
!  vertically accumulated horizontal tracer flux (Hz*u*T/n, Hz*v*T/m)  !
!  in two-way nesting refinement:                                      !
!                                                                      !
!  coarse grid,  t(:,jb,:,nstp,:) = t(:,jb,:,nstp,:) - FacJ    (west,  !
!                                                               east)  !
!                t(ib,:,:,nstp,:) = t(ib,:,:,nstp,:) - FacI    (south, !
!                                                               north) !
!  where                                                               !
!                                                                      !
!                FacJ = (TFF(jb,itrc) - TFC(jb,itrc)) *                !
!                       pm(:,jb) * pn(:,jb) / D(:,jb)                  !
!                                                                      !
!                TFF(ib,itrc) = SUM[SUM[Tflux(ib,k,itrc)]]     finer   !
!                                                              grid    !
!                               for  k=1:N, 1:RefineScale      flux    !
!                                                                      !
!                TFC(ib,itrc) = SUM[Tflux(ib,k,itrc)]          coarser !
!                                                              grid    !
!                               for  k=1:N                     flux    !
!                                                                      !
!  Similarly, for the southern and northern tracer fluxes.             !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ngc        Coarser grid number (integer)                         !
!     ngf        Finer grid number (integer)                           !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!                                                                      !
!  On Output:    (mod_ocean)                                           !
!                                                                      !
!     t          Updated coarse grid tracer values at finer grid       !
!                perimeter                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ngf, model, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
      CALL correct_tracer_tile (ng, ngf, model, tile,                   &
                                LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS)

      RETURN
      END SUBROUTINE correct_tracer
!
!***********************************************************************
      SUBROUTINE correct_tracer_tile (ngc, ngf, model, tile,            &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_grid
      USE mod_ocean
      USE mod_nesting
      USE mod_scalars
      USE mod_stepping

#  ifdef DISTRIBUTE
!
      USE mp_exchange_mod, ONLY : mp_exchange4d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ngc, ngf, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Iedge, Ibc, Ibc_min, Ibc_max, Ibf, Io
      integer :: Jedge, Jbc, Jbc_min, Jbc_max, Jbf, Jo
      integer :: Istr, Iend, Jstr, Jend
      integer :: Istrm2, Iendp2, Jstrm2, Jendp2
      integer :: Tindex, i, ic, isum, itrc, j, jsum, k, half
      integer :: cr, dg, dgcr, rg, rgcr

      real(r8) :: TFC, TFF, Tvalue, cff

      real(r8) :: Dinv(IminS:ImaxS,JminS:JmaxS)
!
!-----------------------------------------------------------------------
!  Correct coarser grid tracer values at finer grid perimeter.
!-----------------------------------------------------------------------
!
!  Determine contact regions where coarse grid is the donor and coarse
!  grid is the receiver..
!
      DO cr=1,Ncontact
        dg=donor_grid(cr)
        rg=receiver_grid(cr)
        IF ((ngc.eq.dg).and.(ngf.eq.rg)) THEN
          dgcr=cr                                   ! coarse is donor
        ELSE IF ((ngc.eq.rg).and.(ngf.eq.dg)) THEN
          rgcr=cr                                   ! coarse is receiver
        END IF
      END DO
!
!  Set tile starting and ending indices for coarser grid.
!
      Istr  =BOUNDS(ngc)%Istr  (tile)
      Iend  =BOUNDS(ngc)%Iend  (tile)
      Jstr  =BOUNDS(ngc)%Jstr  (tile)
      Jend  =BOUNDS(ngc)%Jend  (tile)
!
      Istrm2=BOUNDS(ngc)%Istrm2(tile)
      Iendp2=BOUNDS(ngc)%Iendp2(tile)
      Jstrm2=BOUNDS(ngc)%Jstrm2(tile)
      Jendp2=BOUNDS(ngc)%Jendp2(tile)
!
!  Compute coarser grid inverse water colunm thickness.
!
      DO j=Jstrm2,Jendp2
        DO i=Istrm2,Iendp2
          cff=GRID(ngc)%Hz(i,j,1)
          DO k=2,N(rg)
            cff=cff+GRID(ngc)%Hz(i,j,k)
          END DO
          Dinv(i,j)=1.0_r8/cff
        END DO
      END DO
!
!  Set finer grid center (half) and offset indices (Io and Jo) for
!  coarser grid (I,J) coordinates.
!
      half=(RefineScale(ngf)-1)/2
      Io=half+1
      Jo=half+1
!
!  Set coarse grid tracer index to correct. Since the exchange of data
!  is done at the bottom of main3d, we need to use the newest time
!  index, I think.
!
      Tindex=nstp(ngc)                ! HGA: Why this index is stable?
!!    Tindex=nnew(ngc)                ! Gets a lot of noise at boundary
!
!=======================================================================
!  Compute vertically integrated horizontal advective tracer flux for
!  coarser at the finer grid physical boundary.  Then, correct coarser
!  grid tracer values at that boundary.
!=======================================================================
!
!  Initialize tracer counter index. The "tclm" array is only allocated
!  to the NTCLM fields that need to be processed. This is done to
!  reduce memory.
!
      ic=0
!
      T_LOOP : DO itrc=1,NT(ngc)
        ic=ic+1
!
!-----------------------------------------------------------------------
!  Finer grid western boundary.
!-----------------------------------------------------------------------
!
        Ibc=I_left(ngf)
        Jbc_min=J_bottom(ngf)
        Jbc_max=J_top(ngf)-1              ! interior points, no top
!                                           left corner
        DO Jbc=Jstr,Jend
          IF (((Istr.le.Ibc-1).and.(Ibc-1.le.Iend)).and.                &
     &        ((Jbc_min.le.Jbc).and.(Jbc.le.Jbc_max))) THEN
!
!  Sum vertically coarse grid horizontal advective tracer flux,
!  Hz*u*T/n, from last time-step.
!
            TFC=0.0_r8
            DO k=1,N(ngc)
              TFC=TFC+BRY_CONTACT(iwest,rgcr)%Tflux(Jbc,k,itrc)
            END DO
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal J-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the J-edge.
!
            TFF=0.0_r8
            Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
            DO jsum=-half,half
              Jbf=Jedge+jsum
              DO k=1,N(ngf)
                TFF=TFF+BRY_CONTACT(iwest,dgcr)%Tflux(Jbf,k,itrc)
              END DO
            END DO
!
!  Correct coarse grid tracer at the finer grid western boundary.
!
            cff=GRID(ngc)%pm(Ibc-1,Jbc)*                                &
     &          GRID(ngc)%pn(Ibc-1,Jbc)*                                &
     &          Dinv(Ibc-1,Jbc)
            DO k=1,N(ngc)
              Tvalue=MAX(0.0_r8,                                        &
     &                   OCEAN(ngc)%t(Ibc-1,Jbc,k,Tindex,itrc)-         &
     &                   cff*(TFF-TFC))
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
                Tvalue=Tvalue+                                          &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc-1,Jbc,k,itrc)*   &
     &                 (CLIMA(ngc)%tclm(Ibc-1,Jbc,k,ic)-Tvalue)
              END IF
#  ifdef MASKING
              Tvalue=Tvalue*GRID(ngc)%rmask(Ibc-1,Jbc)
#  endif
              OCEAN(ngc)%t(Ibc-1,Jbc,k,Tindex,itrc)=Tvalue
            END DO
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Finer grid eastern boundary.
!-----------------------------------------------------------------------
!
        Ibc=I_right(ngf)
        Jbc_min=J_bottom(ngf)
        Jbc_max=J_top(ngf)-1              ! interior points, no top
!                                           right corner
        DO Jbc=Jstr,Jend
          IF (((Istr.le.Ibc).and.(Ibc.le.Iend)).and.                    &
     &        ((Jbc_min.le.Jbc).and.(Jbc.le.Jbc_max))) THEN
!
!  Sum vertically coarse grid horizontal advective tracer flux,
!  Hz*u*T/n, from last time-step.
!
            TFC=0.0_r8
            DO k=1,N(ngc)
              TFC=TFC+BRY_CONTACT(ieast,rgcr)%Tflux(Jbc,k,itrc)
            END DO
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal J-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the J-edge.
!
            TFF=0.0_r8
            Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
            DO jsum=-half,half
              Jbf=Jedge+jsum
              DO k=1,N(ngf)
                TFF=TFF+BRY_CONTACT(ieast,dgcr)%Tflux(Jbf,k,itrc)
              END DO
            END DO
!
!  Correct coarse grid tracer at the finer grid eastern boundary.
!
            cff=GRID(ngc)%pm(Ibc,Jbc)*                                  &
     &          GRID(ngc)%pn(Ibc,Jbc)*                                  &
     &          Dinv(Ibc,Jbc)
            DO k=1,N(ngc)
              Tvalue=MAX(0.0_r8,                                        &
     &                   OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-           &
     &                   cff*(TFF-TFC))
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
                Tvalue=Tvalue+                                          &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc)*     &
     &                 (CLIMA(ngc)%tclm(Ibc,Jbc,k,ic)-Tvalue)
              END IF
#  ifdef MASKING
              Tvalue=Tvalue*GRID(ngc)%rmask(Ibc,Jbc)
#  endif
              OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)=Tvalue
            END DO
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Finer grid southern boundary.
!-----------------------------------------------------------------------
!
        Jbc=J_bottom(ngf)
        Ibc_min=I_left(ngf)
        Ibc_max=I_right(ngf)-1            ! interior points, no bottom
!                                           right corner
        DO Ibc=Istr,Iend
          IF (((Ibc_min.le.Ibc).and.(Ibc.le.Ibc_max)).and.              &
     &        ((Jstr.le.Jbc-1).and.(Jbc-1.le.Jend))) THEN
!
!  Sum vertically coarse grid horizontal advective tracer flux,
!  Hz*v*T/m, from last time-step.
!
            TFC=0.0_r8
            DO k=1,N(ngc)
              TFC=TFC+BRY_CONTACT(isouth,rgcr)%Tflux(Ibc,k,itrc)
            END DO
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal I-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the I-edge.
!
            TFF=0.0_r8
            Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
            DO isum=-half,half
              Ibf=Iedge+isum
              DO k=1,N(ngf)
                TFF=TFF+BRY_CONTACT(isouth,dgcr)%Tflux(Ibf,k,itrc)
              END DO
            END DO
!
!  Correct coarse grid tracer at the finer grid southern boundary.
!
            cff=GRID(ngc)%pm(Ibc,Jbc-1)*                                &
     &          GRID(ngc)%pn(Ibc,Jbc-1)*                                &
     &          Dinv(Ibc,Jbc-1)
            DO k=1,N(ngc)
              Tvalue=MAX(0.0_r8,                                        &
     &                   OCEAN(ngc)%t(Ibc,Jbc-1,k,Tindex,itrc)-         &
     &                   cff*(TFF-TFC))
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
                Tvalue=Tvalue+                                          &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc-1,k,itrc)*   &
     &                 (CLIMA(ngc)%tclm(Ibc,Jbc-1,k,ic)-Tvalue)
              END IF
#  ifdef MASKING
              Tvalue=Tvalue*GRID(ngc)%rmask(Ibc,Jbc-1)
#  endif
              OCEAN(ngc)%t(Ibc,Jbc-1,k,Tindex,itrc)=Tvalue
            END DO
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Finer grid northern boundary.
!-----------------------------------------------------------------------
!
        Jbc=J_top(ngf)
        Ibc_min=I_left(ngf)
        Ibc_max=I_right(ngf)-1            ! interior points, no top
!                                           right corner
        DO Ibc=Istr,Iend
          IF (((Ibc_min.le.Ibc).and.(Ibc.le.Ibc_max)).and.              &
     &        ((Jstr.le.Jbc).and.(Jbc.le.Jend))) THEN
!
!  Sum vertically coarse grid horizontal advective tracer flux,
!  Hz*v*T/m, from last time-step.
!
            TFC=0.0_r8
            DO k=1,N(ngc)
              TFC=TFC+BRY_CONTACT(inorth,rgcr)%Tflux(Ibc,k,itrc)
            END DO
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal I-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the I-edge.
!
            TFF=0.0_r8
            Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
            DO isum=-half,half
              Ibf=Iedge+isum
              DO k=1,N(ngf)
                TFF=TFF+BRY_CONTACT(inorth,dgcr)%Tflux(Ibf,k,itrc)
              END DO
            END DO
!
!  Correct coarse grid tracer at the finer grid northern boundary.
!
            cff=GRID(ngc)%pm(Ibc,Jbc)*                                  &
     &          GRID(ngc)%pn(Ibc,Jbc)*                                  &
     &          Dinv(Ibc,Jbc)
            DO k=1,N(ngc)
              Tvalue=MAX(0.0_r8,                                        &
     &                   OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-           &
     &                   cff*(TFF-TFC))
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
                Tvalue=Tvalue+                                          &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc)*     &
     &                 (CLIMA(ngc)%tclm(Ibc,Jbc,k,ic)-Tvalue)
              END IF
#  ifdef MASKING
              Tvalue=Tvalue*GRID(ngc)%rmask(Ibc,Jbc)
#  endif
              OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)=Tvalue
            END DO
          END IF
        END DO
      END DO T_LOOP

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL mp_exchange4d (ngc, tile, model, 1,                          &
     &                    LBi, UBi, LBj, UBj, 1, N(ngc),                &
     &                    1, NT(ngc),                                   &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ngc), NSperiodic(ngc),             &
     &                    OCEAN(ngc)%t(:,:,:,Tindex,:))
#  endif

      RETURN
      END SUBROUTINE correct_tracer_tile
# endif
!
      SUBROUTINE fine2coarse (ng, model, vtype, tile)
!
!=======================================================================
!                                                                      !
!  This routine replaces interior coarse grid data with the refined    !
!  averaged values: two-way nesting.                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Refinement grid number (integer)                      !
!     model      Calling model identifier (integer)                    !
!     vtype      State variables to process (integer):                 !
!                  vtype = r2dvar      2D state variables              !
!                  vtype = r3dvar      3D state variables              !
!     tile       Domain tile partition (integer)                       !
!                                                                      !
!  On Output:    (mod_coupling, mod_ocean)                             !
!                                                                      !
!                Updated state variable with average refined grid      !
!                  solution                                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_nesting
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE exchange_2d_mod
# ifdef SOLVE3D
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#  endif
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, vtype, tile
!
!  Local variable declarations.
!
      logical :: AreaAvg
      integer :: LBiD, UBiD, LBjD, UBjD
      integer :: LBiR, UBiR, LBjR, UBjR
      integer :: Dindex2d, Rindex2d
# ifdef SOLVE3D
      integer :: Dindex3d, Rindex3d
# endif
      integer :: cr, dg, k, rg, nrec, rec
# ifdef SOLVE3D
      integer :: itrc
# endif
!
!-----------------------------------------------------------------------
!  Average interior fine grid state variable data to the coarse grid
!  location. Then, replace coarse grid values with averaged data.
!-----------------------------------------------------------------------
!
      DO cr=1,Ncontact
!
!  Get data donor and data receiver grid numbers.
!
        dg=Rcontact(cr)%donor_grid
        rg=Rcontact(cr)%receiver_grid
!
!  Process contact region if the current refinement grid "ng" is the
!  donor grid.  The coarse grid "rg" is the receiver grid and the
!  contact structure has all the information necessary for fine to
!  coarse coupling. The donor grid size is always smaller than the
!  receiver coarser grid.
!
        IF ((ng.eq.dg).and.(DXmax(dg).lt.DXmax(rg))) THEN
!
!  Set donor and receiver grids lower and upper array indices.
!
          LBiD=BOUNDS(dg)%LBi(tile)
          UBiD=BOUNDS(dg)%UBi(tile)
          LBjD=BOUNDS(dg)%LBj(tile)
          UBjD=BOUNDS(dg)%UBj(tile)
!
          LBiR=BOUNDS(rg)%LBi(tile)
          UBiR=BOUNDS(rg)%UBi(tile)
          LBjR=BOUNDS(rg)%LBj(tile)
          UBjR=BOUNDS(rg)%UBj(tile)
!
!  Report.
!
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            IF (Master.and.(vtype.eq.r2dvar)) THEN
              WRITE (stdout,10) dg, rg, cr
  10          FORMAT (6x,'FINE2COARSE - exchanging data between grids:',&
     &                ' dg = ',i2.2,' and rg = ',i2.2,'  at cr = ',i2.2)
            END IF
          END IF
!
!  Set state variable indices to process for donor and receiver grids.
!  Since the exchange of data is done at the bottom of main2d/main3d,
!  we need to use the newest time indices.
!
          Dindex2d=knew(dg)         ! Donor    2D variables index
          Rindex2d=knew(rg)         ! Receiver 3D variables index
# ifdef SOLVE3D
          Dindex3d=nnew(dg)         ! Donor    3D variables index
          Rindex3d=nnew(rg)         ! Receiver 3D variables index
# endif
!
!-----------------------------------------------------------------------
!  Process 2D state variables.
!-----------------------------------------------------------------------
!
          IF (vtype.eq.r2dvar) THEN
!
!  Free-surface.
!
            AreaAvg=.FALSE.
# ifdef SOLVE3D
            CALL fine2coarse2d (rg, dg, model, tile,                    &
     &                          r2dvar, 'Zt_avg1',                      &
     &                          AreaAvg, RefineScale(dg),               &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBiD, UBiD, LBjD, UBjD,                 &
     &                          LBiR, UBiR, LBjR, UBjR,                 &
     &                          GRID(dg)%om_r,                          &
     &                          GRID(dg)%on_r,                          &
     &                          GRID(rg)%pm,                            &
     &                          GRID(rg)%pn,                            &
#  ifdef MASKING
     &                          GRID(dg)%rmask_full,                    &
     &                          GRID(rg)%rmask_full,                    &
#  endif
     &                          COUPLING(dg)%Zt_avg1,                   &
     &                          COUPLING(rg)%Zt_avg1)
# else
            CALL fine2coarse2d (rg, dg, model, tile,                    &
     &                          r2dvar, Vname(1,idFsur),                &
     &                          AreaAvg, RefineScale(dg),               &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBiD, UBiD, LBjD, UBjD,                 &
     &                          LBiR, UBiR, LBjR, UBjR,                 &
     &                          GRID(dg)%om_r,                          &
     &                          GRID(dg)%on_r,                          &
     &                          GRID(rg)%pm,                            &
     &                          GRID(rg)%pn,                            &
#  ifdef MASKING
     &                          GRID(dg)%rmask_full,                    &
     &                          GRID(rg)%rmask_full,                    &
#  endif
     &                          OCEAN(dg)%zeta(:,:,Dindex2d),           &
     &                          OCEAN(rg)%zeta(:,:,Rindex2d))
# endif
            IF (exit_flag.ne.NoError) RETURN
!
!  Process 2D momentum components (ubar,vbar).
!
            AreaAvg=.FALSE.
            CALL fine2coarse2d (rg, dg, model, tile,                    &
     &                          u2dvar, Vname(1,idUbar),                &
     &                          AreaAvg, RefineScale(dg),               &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBiD, UBiD, LBjD, UBjD,                 &
     &                          LBiR, UBiR, LBjR, UBjR,                 &
     &                          GRID(dg)%om_u,                          &
     &                          GRID(dg)%on_u,                          &
     &                          GRID(rg)%pm,                            &
     &                          GRID(rg)%pn,                            &
# ifdef MASKING
     &                          GRID(dg)%umask_full,                    &
     &                          GRID(rg)%umask_full,                    &
# endif
     &                          OCEAN(dg)%ubar(:,:,Dindex2d),           &
# ifdef SOLVE3D
     &                          OCEAN(rg)%ubar(:,:,1),                  &
     &                          OCEAN(rg)%ubar(:,:,2))
# else
     &                          OCEAN(rg)%ubar(:,:,Rindex2d))
# endif
            IF (exit_flag.ne.NoError) RETURN
!
            CALL fine2coarse2d (rg, dg, model, tile,                    &
     &                          v2dvar, Vname(1,idVbar),                &
     &                          AreaAvg, RefineScale(dg),               &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBiD, UBiD, LBjD, UBjD,                 &
     &                          LBiR, UBiR, LBjR, UBjR,                 &
     &                          GRID(dg)%om_v,                          &
     &                          GRID(dg)%on_v,                          &
     &                          GRID(rg)%pm,                            &
     &                          GRID(rg)%pn,                            &
# ifdef MASKING
     &                          GRID(dg)%vmask_full,                    &
     &                          GRID(rg)%vmask_full,                    &
# endif
     &                          OCEAN(dg)%vbar(:,:,Dindex2d),           &
# ifdef SOLVE3D
     &                          OCEAN(rg)%vbar(:,:,1),                  &
     &                          OCEAN(rg)%vbar(:,:,2))
# else
     &                          OCEAN(rg)%vbar(:,:,Rindex2d))
# endif
            IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Process 3D state variables.
!-----------------------------------------------------------------------
!
          ELSE IF (vtype.eq.r3dvar) THEN
!
!  Tracer type-variables.
!
            AreaAvg=.FALSE.
            DO itrc=1,NT(rg)
              CALL fine2coarse3d (rg, dg, model, tile,                  &
     &                            r3dvar, Vname(1,idTvar(itrc)),        &
     &                            AreaAvg, RefineScale(dg),             &
     &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
     &                            LBiD, UBiD, LBjD, UBjD, 1, N(dg),     &
     &                            LBiR, UBiR, LBjR, UBjR, 1, N(rg),     &
     &                            GRID(dg)%om_r,                        &
     &                            GRID(dg)%on_r,                        &
     &                            GRID(rg)%pm,                          &
     &                            GRID(rg)%pn,                          &
#  ifdef MASKING
     &                            GRID(dg)%rmask_full,                  &
     &                            GRID(rg)%rmask_full,                  &
#  endif
     &                            OCEAN(dg)%t(:,:,:,Dindex3d,itrc),     &
     &                            OCEAN(rg)%t(:,:,:,Rindex3d,itrc))
              IF (exit_flag.ne.NoError) RETURN
            END DO
!
!  Process 3D momentum components (u, v).
!
            AreaAvg=.FALSE.
            CALL fine2coarse3d (rg, dg, model, tile,                    &
     &                          u3dvar, Vname(1,idUvel),                &
     &                          AreaAvg, RefineScale(dg),               &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBiD, UBiD, LBjD, UBjD, 1, N(dg),       &
     &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
     &                          GRID(dg)%om_u,                          &
     &                          GRID(dg)%on_u,                          &
     &                          GRID(rg)%pm,                            &
     &                          GRID(rg)%pn,                            &
#  ifdef MASKING
     &                          GRID(dg)%umask_full,                    &
     &                          GRID(rg)%umask_full,                    &
#  endif
     &                          OCEAN(dg)%u(:,:,:,Dindex3d),            &
     &                          OCEAN(rg)%u(:,:,:,Rindex3d))
            IF (exit_flag.ne.NoError) RETURN
!
            CALL fine2coarse3d (rg, dg, model, tile,                    &
     &                          v3dvar, Vname(1,idVvel),                &
     &                          AreaAvg, RefineScale(dg),               &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBiD, UBiD, LBjD, UBjD, 1, N(dg),       &
     &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
     &                          GRID(dg)%om_v,                          &
     &                          GRID(dg)%on_v,                          &
     &                          GRID(rg)%pm,                            &
     &                          GRID(rg)%pn,                            &
#  ifdef MASKING
     &                          GRID(dg)%vmask_full,                    &
     &                          GRID(rg)%vmask_full,                    &
#  endif
     &                          OCEAN(dg)%v(:,:,:,Dindex3d),            &
     &                          OCEAN(rg)%v(:,:,:,Rindex3d))
            IF (exit_flag.ne.NoError) RETURN
# endif
          END IF
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
          IF (EWperiodic(rg).or.NSperiodic(rg)) THEN
# ifdef SOLVE3D
            CALL exchange_r2d_tile (rg, tile,                           &
     &                              LBiR, UBiR, LBjR, UBjR,             &
     &                              COUPLING(rg)%Zt_avg1)
            DO k=1,2
              CALL exchange_u2d_tile (rg, tile,                         &
     &                                LBiR, UBiR, LBjR, UBjR,           &
     &                                OCEAN(rg)%ubar(:,:,k))
              CALL exchange_v2d_tile (rg, tile,                         &
     &                                LBiR, UBiR, LBjR, UBjR,           &
     &                                OCEAN(rg)%vbar(:,:,k))
            END DO
# else
            CALL exchange_r2d_tile (rg, tile,                           &
     &                              LBiR, UBiR, LBjR, UBjR,             &
     &                              OCEAN(rg)%zeta(:,:,Rindex2d))
            CALL exchange_u2d_tile (rg, tile,                           &
     &                              LBiR, UBiR, LBjR, UBjR,             &
     &                              OCEAN(rg)%ubar(:,:,Rindex2d))
            CALL exchange_v2d_tile (rg, tile,                           &
     &                              LBiR, UBiR, LBjR, UBjR,             &
     &                              OCEAN(rg)%vbar(:,:,Rindex2d))

# endif
# ifdef SOLVE3D
            CALL exchange_u3d_tile (rg, tile,                           &
     &                              LBiR, UBiR, LBjR, UBjR, 1, N(rg),   &
     &                              OCEAN(rg)%u(:,:,:,Rindex3d))
            CALL exchange_v3d_tile (rg, tile,                           &
     &                              LBiR, UBiR, LBjR, UBjR, 1, N(rg),   &
     &                              OCEAN(rg)%v(:,:,:,Rindex3d))
            DO itrc=1,NT(rg)
              CALL exchange_r3d_tile (rg, tile,                         &
     &                                LBiR, UBiR, LBjR, UBjR, 1, N(rg), &
     &                                OCEAN(rg)%t(:,:,:,Rindex3d,itrc))
            END DO
# endif
          END IF

# ifdef DISTRIBUTE
!
#  ifdef SOLVE3D
          CALL mp_exchange2d (rg, tile, model, 1,                       &
     &                        LBiR, UBiR, LBjR, UBjR,                   &
     &                        NghostPoints,                             &
     &                        EWperiodic(rg), NSperiodic(rg),           &
     &                        COUPLING(rg)%Zt_avg1)
          CALL mp_exchange2d (rg, tile, model, 4,                       &
     &                        LBiR, UBiR, LBjR, UBjR,                   &
     &                        NghostPoints,                             &
     &                        EWperiodic(rg), NSperiodic(rg),           &
     &                        OCEAN(rg)%ubar(:,:,1),                    &
     &                        OCEAN(rg)%vbar(:,:,1),                    &
     &                        OCEAN(rg)%ubar(:,:,2),                    &
     &                        OCEAN(rg)%vbar(:,:,2))
#  else
          CALL mp_exchange2d (rg, tile, model, 3,                       &
     &                        LBiR, UBiR, LBjR, UBjR,                   &
     &                        NghostPoints,                             &
     &                        EWperiodic(rg), NSperiodic(rg),           &
     &                        OCEAN(rg)%zeta(:,:,Rindex2d),             &
     &                        OCEAN(rg)%ubar(:,:,Rindex2d),             &
     &                        OCEAN(rg)%vbar(:,:,Rindex2d))
#  endif
#  ifdef SOLVE3D
          CALL mp_exchange3d (rg, tile, model, 2,                       &
     &                        LBiR, UBiR, LBjR, UBjR, 1, N(rg),         &
     &                        NghostPoints,                             &
     &                        EWperiodic(rg), NSperiodic(rg),           &
     &                        OCEAN(rg)%u(:,:,:,Rindex3d),              &
     &                        OCEAN(rg)%v(:,:,:,Rindex3d))
          CALL mp_exchange4d (rg, tile, model, 1,                       &
     &                        LBiR, UBiR, LBjR, UBjR, 1, N(rg),         &
     &                        1, NT(rg),                                &
     &                        NghostPoints,                             &
     &                        EWperiodic(rg), NSperiodic(rg),           &
     &                        OCEAN(rg)%t(:,:,:,Rindex3d,:))
#  endif
# endif
        END IF
      END DO

      RETURN
      END SUBROUTINE fine2coarse
!
      SUBROUTINE fine2coarse2d (ng, dg, model, tile,                    &
     &                          gtype, svname,                          &
     &                          AreaAvg, Rscale,                        &
     &                          cr, Npoints, contact,                   &
     &                          LBiF, UBiF, LBjF, UBjF,                 &
     &                          LBiC, UBiC, LBjC, UBjC,                 &
# ifdef DISTRIBUTE
     &                          Adx, Ady,                               &
# else
     &                          dxF, dyF,                               &
# endif
     &                          pmC, pnC,                               &
# ifdef MASKING
#  ifdef DISTRIBUTE
     &                          Amsk,                                   &
#  else
     &                          Fmsk,                                   &
#  endif
     &                          Cmsk,                                   &
# endif
# ifdef DISTRIBUTE
     &                          A,                                      &
# else
     &                          F,                                      &
# endif
     &                          C1, C2)
!
!=======================================================================
!                                                                      !
!  This routine replaces the coarse grid data inside the refinement    !
!  grid interior for a 2D state variable with its refined averaged     !
!  values: two-way nesting.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Coarser grid number (integer)                         !
!     dg         Finer grid number (integer)                           !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     gtype      C-grid variable type (integer)                        !
!     svname     State variable name (string)                          !
!     AreaAvg    Switch for area averaging (logical)                   !
!     Rscale     Refinement grid scale (integer)                       !
!     cr         Contact region number to process (integer)            !
!     Npoints    Number of points in the contact zone (integer)        !
!     contact    Contact zone information variables (T_NGC structure)  !
!     LBiF       Finer grid, I-dimension Lower bound (integer)         !
!     UBiF       Finer grid, I-dimension Upper bound (integer)         !
!     LBjF       Finer grid, J-dimension Lower bound (integer)         !
!     UBjF       Finer grid, J-dimension Upper bound (integer)         !
!     LBiC       Coarser grid, I-dimension Lower bound (integer)       !
!     UBiC       Coarser grid, I-dimension Upper bound (integer)       !
!     LBjC       Coarser grid, J-dimension Lower bound (integer)       !
!     UBjC       Coarser grid, J-dimension Upper bound (integer)       !
# ifdef DISTRIBUTE
!     Adx        Finer grid, X-grid spacing (1/pm: om_r, om_u, om_v)   !
!     Ady        Finer grid, Y-grid spacing (1/pn: on_r, on_u, on_v)   !
# else
!     dxF        Finer grid, X-grid spacing (1/pm: om_r, om_u, om_v)   !
!     dyF        Finer grid, Y-grid spacing (1/pn: on_r, on_u, on_v)   !
# endif
!     pmC        Coarser grid, inverse X-grid spacing (1/dx) at RHO    !
!     pnC        Coarser grid, inverse Y-grid spacing (1/dy) at RHO    !
# ifdef MASKING
#  ifdef DISTRIBUTE
!     Amsk       Finer grid land/sea masking (2D array)                !
#  else
!     Fmsk       Finer grid land/sea masking (2D array)                !
#  endif
!     Cmsk       Coarser grid land/sea masking (2D array)              !
# endif
# ifdef DISTRIBUTE
!     A          Finer grid 2D data                                    !
# else
!     F          Finer grid 2D data                                    !
# endif
!     C1         Coarser grid 2D data, record 1                        !
!     C2         Coarser grid 2D data, record 2 (OPTIONAL)             !
!                                                                      !
!  On Output:    (mod_nesting)                                         !
!                                                                      !
!     C1         Updated Coarser grid 2D data, record 1                !
!     C2         Uodated Coarser grid 2D data, record 2 (OPTIONAL)     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_nesting
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_aggregate2d
# endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: AreaAvg
      integer, intent(in) :: ng, dg, model, tile
      integer, intent(in) :: gtype, cr, Npoints, Rscale
      integer, intent(in) :: LBiF, UBiF, LBjF, UBjF
      integer, intent(in) :: LBiC, UBiC, LBjC, UBjC
!
      character(len=*), intent(in) :: svname
!
      TYPE (T_NGC), intent(in) :: contact(:)
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pmC(LBiC:,LBjC:)
      real(r8), intent(in) :: pnC(LBiC:,LBjC:)
#  ifdef MASKING
      real(r8), intent(in) :: Cmsk(LBiC:,LBjC:)
#   ifdef DISTRIBUTE
      real(r8), intent(in) :: Amsk(LBiF:,LBjF:)
#   else
      real(r8), intent(in) :: Fmsk(LBiF:,LBjF:)
#   endif
#  endif
#  ifdef DISTRIBUTE
      real(r8), intent(in) :: A(LBiF:,LBjF:)
      real(r8), intent(in) :: Adx(LBiF:,LBjF:)
      real(r8), intent(in) :: Ady(LBiF:,LBjF:)
#  else
      real(r8), intent(in) :: F(LBiF:,LBjF:)
      real(r8), intent(in) :: dxF(LBiF:,LBjF:)
      real(r8), intent(in) :: dyF(LBiF:,LBjF:)
#  endif
      real(r8), intent(inout) :: C1(LBiC:,LBjC:)
      real(r8), intent(inout), optional :: C2(LBiC:,LBjC:)
# else
      real(r8), intent(in) :: pmC(LBiC:UBiC,LBjC:UBjC)
      real(r8), intent(in) :: pnC(LBiC:UBiC,LBjC:UBjC)
#  ifdef MASKING
      real(r8), intent(in) :: Cmsk(LBiC:UBiC,LBjC:UBjC)
#   ifdef DISTRIBUTE
      real(r8), intent(in) :: Amsk(LBiF:UBiF,LBjF:UBjF)
#   else
      real(r8), intent(in) :: Fmsk(LBiF:UBiF,LBjF:UBjF)
#   endif
#  endif
#  ifdef DISTRIBUTE
      real(r8), intent(in) :: A(LBiF:UBiF,LBjF:UBjF)
      real(r8), intent(in) :: Adx(LBiF:UBiF,LBjF:UBjF)
      real(r8), intent(in) :: Adx(LBiF:UBiF,LBjF:UBjF)
#  else
      real(r8), intent(in) :: F(LBiF:UBiF,LBjF:UBjF)
      real(r8), intent(in) :: dxF(LBiF:UBiF,LBjF:UBjF)
      real(r8), intent(in) :: dyF(LBiF:UBiF,LBjF:UBjF)
#  endif
      real(r8), intent(inout) :: C1(LBiC:UBiC,LBjC:UBjC)
      real(r8), intent(inout), optional :: C2(LBiC:UBiC,LBjC:UBjC)
# endif
!
!  Local variable declarations.
!
      integer :: Iadd, Ic, Jadd, Jc, half, i, j, m
# ifdef DISTRIBUTE
      integer :: LBi, UBi, LBj, UBj
# endif

      real(r8) :: areaC_inv, my_area, my_areasum, ratio
      real(r8) :: my_avg, my_count, my_sum

# ifdef DISTRIBUTE
      real(r8), allocatable :: F(:,:)
      real(r8), allocatable :: dxF(:,:)
      real(r8), allocatable :: dyF(:,:)
#  ifdef MASKING
      real(r8), allocatable :: Fmsk(:,:)
#  endif
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Average interior fine grid state variable data to the coarse grid
!  location. Then, replace coarse grid values with averaged data.
!-----------------------------------------------------------------------

# ifdef DISTRIBUTE
!
!  Allocate global work array(s).
!
      LBi=BOUNDS(dg)%LBi(-1)
      UBi=BOUNDS(dg)%UBi(-1)
      LBj=BOUNDS(dg)%LBj(-1)
      UBj=BOUNDS(dg)%UBj(-1)
      IF (.not.allocated(F)) THEN
        allocate ( F(LBi:UBi,LBj:UBj) )
      END IF
      IF (AreaAvg) THEN
        IF (.not.allocated(dxF)) THEN
          allocate ( dxF(LBi:UBi,LBj:UBj) )
        END IF
        IF (.not.allocated(dyF)) THEN
          allocate ( dyF(LBi:UBi,LBj:UBj) )
        END IF
      END IF
#  ifdef MASKING
      IF (.not.allocated(Fmsk)) THEN
        allocate ( Fmsk(LBi:UBi,LBj:UBj) )
      END IF
#  endif
!
!  Gather finer grid data from all nodes in the group to build a global
!  array.
!
      CALL mp_aggregate2d (dg, model, gtype,                            &
     &                     LBiF, UBiF, LBjF, UBjF,                      &
     &                     LBi,  UBi,  LBj,  UBj,                       &
     &                     A, F)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (AreaAvg) THEN
        CALL mp_aggregate2d (dg, model, gtype,                          &
     &                       LBiF, UBiF, LBjF, UBjF,                    &
     &                       LBi,  UBi,  LBj,  UBj,                     &
     &                       Adx, dxF)
        IF (exit_flag.ne.NoError) RETURN

        CALL mp_aggregate2d (dg, model, gtype,                          &
     &                       LBiF, UBiF, LBjF, UBjF,                    &
     &                       LBi,  UBi,  LBj,  UBj,                     &
     &                       Ady, dyF)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  ifdef MASKING
!
      CALL mp_aggregate2d (dg, model, gtype,                            &
     &                     LBiF, UBiF, LBjF, UBjF,                      &
     &                     LBi,  UBi,  LBj,  UBj,                       &
     &                     Amsk, Fmsk)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
!
!  Average finer grid data to coarse grid according to the refinement
!  ratio.
!
      half=(Rscale-1)/2
      IF (AreaAvg) THEN               ! area averaging
        DO m=1,Npoints
          i=contact(cr)%Idg(m)
          j=contact(cr)%Jdg(m)
          Ic=contact(cr)%Irg(m)
          Jc=contact(cr)%Jrg(m)
          IF (((Istr.le.Ic).and.(Ic.le.Iend)).and.                      &
     &        ((Jstr.le.Jc).and.(Jc.le.Jend))) THEN
            my_count=0.0_r8
            my_sum=0.0_r8
            my_areasum=0.0_r8
            DO Jadd=-half,half
              DO Iadd=-half,half
                my_area=dxF(i+Iadd,j+Jadd)*dyF(i+Iadd,j+Jadd)
                my_areasum=my_areasum+my_area
# ifdef MASKING
                my_sum=my_sum+                                          &
     &                 F(i+Iadd,j+Jadd)*my_area*                        &
     &                 MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
                my_count=my_count+MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
# else
                my_sum=my_sum+                                          &
     &                 F(i+Iadd,j+Jadd)*my_area
# endif
              END DO
            END DO
            SELECT CASE (gtype)              ! coarse grid inverse area
              CASE (r2dvar)
                areaC_inv=pmC(Ic,Jc)*pnC(Ic,Jc)
              CASE (u2dvar)
                areaC_inv=0.25_r8*(pmC(Ic-1,Jc)+pmC(Ic,Jc))*            &
     &                            (pnC(Ic-1,Jc)+pnC(Ic,Jc))
              CASE (v2dvar)
                areaC_inv=0.25_r8*(pmC(Ic,Jc-1)+pmC(Ic,Jc))*            &
     &                            (pnC(Ic,Jc-1)+pnC(Ic,Jc))
              CASE DEFAULT
                areaC_inv=pmC(Ic,Jc)*pnC(Ic,Jc)
            END SELECT
            ratio=my_areasum*areaC_inv       ! for debugging purposes
            my_avg=my_sum*areaC_inv
# ifdef MASKING
            IF (my_count.gt.0.0_r8) THEN
              my_avg=my_avg*Rscale*Rscale/my_count
            END IF
            my_avg=my_avg*Cmsk(Ic,Jc)
# endif
            C1(Ic,Jc)=my_avg
            IF (PRESENT(C2)) THEN
              C2(Ic,Jc)=my_avg
            END IF
          END IF
        END DO
      ELSE                            ! simple averaging
        DO m=1,Npoints
          i=contact(cr)%Idg(m)
          j=contact(cr)%Jdg(m)
          Ic=contact(cr)%Irg(m)
          Jc=contact(cr)%Jrg(m)
          IF (((Istr.le.Ic).and.(Ic.le.Iend)).and.                      &
     &        ((Jstr.le.Jc).and.(Jc.le.Jend))) THEN
            my_count=0.0_r8
            my_avg=0.0_r8
            my_sum=0.0_r8
            DO Jadd=-half,half
              DO Iadd=-half,half
# ifdef MASKING
                my_sum=my_sum+                                          &
     &                 F(i+Iadd,j+Jadd)*Fmsk(i+Iadd,j+Jadd)
                my_count=my_count+MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
# else
                my_sum=my_sum+                                          &
     &                 F(i+Iadd,j+Jadd)
                my_count=my_count+1.0_r8
# endif
              END DO
            END DO
            IF (my_count.gt.0.0_r8) my_avg=my_sum/my_count
# ifdef MASKING
            my_avg=my_avg*Cmsk(Ic,Jc)
# endif
            C1(Ic,Jc)=my_avg
            IF (PRESENT(C2)) THEN
              C2(Ic,Jc)=my_avg
            END IF
          END IF
        END DO
      END IF

# ifdef DISTRIBUTE
!
!  Deallocate work array.
!
      IF (allocated(F)) THEN
        deallocate (F)
      END IF
      IF (AreaAvg) THEN
        IF (allocated(dxF)) THEN
          deallocate (dxF)
        END IF
        IF (allocated(dyF)) THEN
          deallocate (dyF)
        END IF
      END IF
#  ifdef MASKING
      IF (allocated(Fmsk)) THEN
        deallocate (Fmsk)
      END IF
#  endif
# endif

      RETURN
      END SUBROUTINE fine2coarse2d
!
# ifdef SOLVE3D
      SUBROUTINE fine2coarse3d (ng, dg, model, tile,                    &
     &                          gtype, svname,                          &
     &                          AreaAvg, Rscale,                        &
     &                          cr, Npoints, contact,                   &
     &                          LBiF, UBiF, LBjF, UBjF, LBkF, UBkF,     &
     &                          LBiC, UBiC, LBjC, UBjC, LBkC, UBkC,     &
#  ifdef DISTRIBUTE
     &                          Adx, Ady,                               &
#  else
     &                          dxF, dyF,                               &
#  endif
     &                          pmC, pnC,                               &
#  ifdef MASKING
#   ifdef DISTRIBUTE
     &                          Amsk,                                   &
#   else
     &                          Fmsk,                                   &
#   endif
     &                          Cmsk,                                   &
#  endif
#  ifdef DISTRIBUTE
     &                          A,                                      &
#  else
     &                          F,                                      &
#  endif
     &                          C)
!
!=======================================================================
!                                                                      !
!  This routine replaces the coarse grid data inside the refinement    !
!  grid interior for a 3D state variable with its refined averaged     !
!  values: two-way nesting.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Coarser grid number (integer)                         !
!     dg         Finer grid number (integer)                           !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     gtype      C-grid variable type (integer)                        !
!     svname     State variable name (string)                          !
!     AreaAvg    Switch for area averaging (logical)                   !
!     Rscale     Refinement grid scale (integer)                       !
!     cr         Contact region number to process (integer)            !
!     Npoints    Number of points in the contact zone (integer)        !
!     contact    Contact zone information variables (T_NGC structure)  !
!     LBiF       Finer grid, I-dimension Lower bound (integer)         !
!     UBiF       Finer grid, I-dimension Upper bound (integer)         !
!     LBjF       Finer grid, J-dimension Lower bound (integer)         !
!     UBjF       Finer grid, J-dimension Upper bound (integer)         !
!     LBkF       Finer grid, K-dimension Lower bound (integer)         !
!     UBkF       Finer grid, K-dimension Upper bound (integer)         !
!     LBiC       Coarser grid, I-dimension Lower bound (integer)       !
!     UBiC       Coarser grid, I-dimension Upper bound (integer)       !
!     LBjC       Coarser grid, J-dimension Lower bound (integer)       !
!     UBjC       Coarser grid, J-dimension Upper bound (integer)       !
!     LBkC       Coarser grid, K-dimension Lower bound (integer)       !
!     UBkC       Coarser grid, K-dimension Upper bound (integer)       !
#  ifdef DISTRIBUTE
!     Adx        Finer grid, X-grid spacing (1/pm: om_r, om_u, om_v)   !
!     Ady        Finer grid, Y-grid spacing (1/pn: on_r, on_u, on_v)   !
#  else
!     dxF        Finer grid, X-grid spacing (1/pm: om_r, om_u, om_v)   !
!     dyF        Finer grid, Y-grid spacing (1/pn: on_r, on_u, on_v)   !
#  endif
!     pmC        Coarser grid, inverse X-grid spacing (1/dx) at RHO    !
!     pnC        Coarser grid, inverse Y-grid spacing (1/dy) at RHO    !
#  ifdef MASKING
#   ifdef DISTRIBUTE
!     Amsk       Finer grid land/sea masking (2D array)                !
#   else
!     Fmsk       Finer grid land/sea masking (2D array)                !
#   endif
!     Cmsk       Coarser grid land/sea masking (2D array)              !
#  endif
#  ifdef DISTRIBUTE
!     A          Finer grid 2D data                                    !
#  else
!     F          Finer grid 2D data                                    !
#  endif
!     C          Coarser grid 3D data                                  !
!                                                                      !
!  On Output:    (mod_nesting)                                         !
!                                                                      !
!     C          Updated Coarser grid 3D data                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_nesting
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_aggregate2d
      USE distribute_mod, ONLY : mp_aggregate3d
#  endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: AreaAvg
      integer, intent(in) :: ng, dg, model, tile
      integer, intent(in) :: gtype, cr, Npoints, Rscale
      integer, intent(in) :: LBiF, UBiF, LBjF, UBjF, LBkF, UBkF
      integer, intent(in) :: LBiC, UBiC, LBjC, UBjC, LBkC, UBkC
!
      character(len=*), intent(in) :: svname
!
      TYPE (T_NGC), intent(in) :: contact(:)
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pmC(LBiC:,LBjC:)
      real(r8), intent(in) :: pnC(LBiC:,LBjC:)
#   ifdef MASKING
      real(r8), intent(in) :: Cmsk(LBiC:,LBjC:)
#    ifdef DISTRIBUTE
      real(r8), intent(in) :: Amsk(LBiF:,LBjF:)
#    else
      real(r8), intent(in) :: Fmsk(LBiF:,LBjF:)
#    endif
#   endif
#   ifdef DISTRIBUTE
      real(r8), intent(in) :: A(LBiF:,LBjF:,LBkF:)
      real(r8), intent(in) :: Adx(LBiF:,LBjF:)
      real(r8), intent(in) :: Ady(LBiF:,LBjF:)
#   else
      real(r8), intent(in) :: F(LBiF:,LBjF:,LBkF:)
      real(r8), intent(in) :: dxF(LBiF:,LBjF:)
      real(r8), intent(in) :: dyF(LBiF:,LBjF:)
#   endif
      real(r8), intent(inout) :: C(LBiC:,LBjC:,LBkC:)
#  else
      real(r8), intent(in) :: pmC(LBiC:UBiC,LBjC:UBjC)
      real(r8), intent(in) :: pnC(LBiC:UBiC,LBjC:UBjC)
#   ifdef MASKING
      real(r8), intent(in) :: Cmsk(LBiC:UBiC,LBjC:UBjC)
#    ifdef DISTRIBUTE
      real(r8), intent(in) :: Amsk(LBiF:UBiF,LBjF:UBjF)
#    else
      real(r8), intent(in) :: Fmsk(LBiF:UBiF,LBjF:UBjF)
#    endif
#   endif
#   ifdef DISTRIBUTE
      real(r8), intent(in) :: A(LBiF:UBiF,LBjF:UBjF,LBkF:UBkF)
      real(r8), intent(in) :: Adx(LBiF:UBiF,LBjF:UBjF)
      real(r8), intent(in) :: Adx(LBiF:UBiF,LBjF:UBjF)
#   else
      real(r8), intent(in) :: F(LBiF:UBiF,LBjF:UBjF,LBkF:UBkF)
      real(r8), intent(in) :: dxF(LBiF:UBiF,LBjF:UBjF)
      real(r8), intent(in) :: dyF(LBiF:UBiF,LBjF:UBjF)
#   endif
      real(r8), intent(inout) :: C(LBiC:UBiC,LBjC:UBjC,LBkC:UBkC)
#  endif
!
!  Local variable declarations.
!
      integer :: Iadd, Ic, Jadd, Jc, half, i, j, k, m
#  ifdef DISTRIBUTE
      integer :: LBi, UBi, LBj, UBj
#  endif

      real(r8) :: areaC_inv, my_area, my_areasum, ratio
      real(r8) :: my_avg, my_count, my_sum

#  ifdef DISTRIBUTE
      real(r8), allocatable :: F(:,:,:)
      real(r8), allocatable :: dxF(:,:)
      real(r8), allocatable :: dyF(:,:)
#   ifdef MASKING
      real(r8), allocatable :: Fmsk(:,:)
#   endif
#  endif

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Average interior fine grid state variable data to the coarse grid
!  location. Then, replace coarse grid values with averaged data.
!-----------------------------------------------------------------------

#  ifdef DISTRIBUTE
!
!  Allocate global work array(s).
!
      LBi=BOUNDS(dg)%LBi(-1)
      UBi=BOUNDS(dg)%UBi(-1)
      LBj=BOUNDS(dg)%LBj(-1)
      UBj=BOUNDS(dg)%UBj(-1)
      IF (.not.allocated(F)) THEN
        allocate ( F(LBi:UBi,LBj:UBj,LBkF:UBkF) )
      END IF
      IF (AreaAvg) THEN
        IF (.not.allocated(dxF)) THEN
          allocate ( dxF(LBi:UBi,LBj:UBj) )
        END IF
        IF (.not.allocated(dyF)) THEN
          allocate ( dyF(LBi:UBi,LBj:UBj) )
        END IF
      END IF
#   ifdef MASKING
      IF (.not.allocated(Fmsk)) THEN
        allocate ( Fmsk(LBi:UBi,LBj:UBj) )
      END IF
#   endif
!
!  Gather finer grid data from all nodes in the group to build a global
!  array.
!
      CALL mp_aggregate3d (dg, model, gtype,                            &
     &                     LBiF, UBiF, LBjF, UBjF,                      &
     &                     LBi,  UBi,  LBj,  UBj,                       &
     &                     LBkF, UBkF,                                  &
     &                     A, F)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (AreaAvg) THEN
        CALL mp_aggregate2d (dg, model, gtype,                          &
     &                       LBiF, UBiF, LBjF, UBjF,                    &
     &                       LBi,  UBi,  LBj,  UBj,                     &
     &                       Adx, dxF)
        IF (exit_flag.ne.NoError) RETURN

        CALL mp_aggregate2d (dg, model, gtype,                          &
     &                       LBiF, UBiF, LBjF, UBjF,                    &
     &                       LBi,  UBi,  LBj,  UBj,                     &
     &                       Ady, dyF)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#   ifdef MASKING
!
      CALL mp_aggregate2d (dg, model, gtype,                            &
     &                     LBiF, UBiF, LBjF, UBjF,                      &
     &                     LBi,  UBi,  LBj,  UBj,                       &
     &                     Amsk, Fmsk)
      IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
!
!  Average finer grid data to coarse grid according to the refinement
!  ratio.
!
      half=(Rscale-1)/2
      IF (AreaAvg) THEN               ! area averaging
        DO k=LBkC,UBkC
          DO m=1,Npoints
            i=contact(cr)%Idg(m)
            j=contact(cr)%Jdg(m)
            Ic=contact(cr)%Irg(m)
            Jc=contact(cr)%Jrg(m)
            IF (((Istr.le.Ic).and.(Ic.le.Iend)).and.                    &
     &          ((Jstr.le.Jc).and.(Jc.le.Jend))) THEN
              my_count=0.0_r8
              my_sum=0.0_r8
              my_areasum=0.0_r8
              DO Jadd=-half,half
                DO Iadd=-half,half
                  my_area=dxF(i+Iadd,j+Jadd)*dyF(i+Iadd,j+Jadd)
                  my_areasum=my_areasum+my_area
#  ifdef MASKING
                  my_sum=my_sum+                                        &
     &                   F(i+Iadd,j+Jadd,k)*my_area*                    &
     &                   MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
                  my_count=my_count+MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
#  else
                  my_sum=my_sum+                                        &
     &                   F(i+Iadd,j+Jadd,k)*my_area
#  endif
                END DO
              END DO
              SELECT CASE (gtype)            ! coarse grid inverse area
                CASE (r3dvar)
                  areaC_inv=pmC(Ic,Jc)*pnC(Ic,Jc)
                CASE (u3dvar)
                  areaC_inv=0.25_r8*(pmC(Ic-1,Jc)+pmC(Ic,Jc))*          &
     &                              (pnC(Ic-1,Jc)+pnC(Ic,Jc))
                CASE (v3dvar)
                  areaC_inv=0.25_r8*(pmC(Ic,Jc-1)+pmC(Ic,Jc))*          &
     &                              (pnC(Ic,Jc-1)+pnC(Ic,Jc))
                CASE DEFAULT
                  areaC_inv=pmC(Ic,Jc)*pnC(Ic,Jc)
              END SELECT
              ratio=my_areasum*areaC_inv     ! for debugging purposes
              my_avg=my_sum*areaC_inv
#  ifdef MASKING
              IF (my_count.gt.0.0_r8) THEN
                my_avg=my_avg*Rscale*Rscale/my_count
              END IF
              my_avg=my_avg*Cmsk(Ic,Jc)
#  endif
              C(Ic,Jc,k)=my_avg
            END IF
          END DO
        END DO
      ELSE                            ! simple averaging
        DO k=LBkC,UBkC
          DO m=1,Npoints
            i=contact(cr)%Idg(m)
            j=contact(cr)%Jdg(m)
            Ic=contact(cr)%Irg(m)
            Jc=contact(cr)%Jrg(m)
            IF (((Istr.le.Ic).and.(Ic.le.Iend)).and.                    &
     &          ((Jstr.le.Jc).and.(Jc.le.Jend))) THEN
              my_count=0.0_r8
              my_avg=0.0_r8
              my_sum=0.0_r8
              DO Jadd=-half,half
                DO Iadd=-half,half
#  ifdef MASKING
                  my_sum=my_sum+                                        &
     &                   F(i+Iadd,j+Jadd,k)*Fmsk(i+Iadd,j+Jadd)
                  my_count=my_count+MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
#  else
                  my_sum=my_sum+                                        &
     &                   F(i+Iadd,j+Jadd,k)
                  my_count=my_count+1.0_r8
#  endif
                END DO
              END DO
              IF (my_count.gt.0.0_r8) my_avg=my_sum/my_count
#  ifdef MASKING
              my_avg=my_avg*Cmsk(Ic,Jc)
#  endif
              C(Ic,Jc,k)=my_avg
            END IF
          END DO
        END DO
      END IF

#  ifdef DISTRIBUTE
!
!  Deallocate work array.
!
      IF (allocated(F)) THEN
        deallocate (F)
      END IF
      IF (AreaAvg) THEN
        IF (allocated(dxF)) THEN
          deallocate (dxF)
        END IF
        IF (allocated(dyF)) THEN
          deallocate (dyF)
        END IF
      END IF
#   ifdef MASKING
      IF (allocated(Fmsk)) THEN
        deallocate (Fmsk)
      END IF
#   endif
#  endif

      RETURN
      END SUBROUTINE fine2coarse3d
# endif
!
      SUBROUTINE get_contact2d (dg, model, tile,                        &
     &                          gtype, svname,                          &
     &                          cr, Npoints, contact,                   &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Ad, Ac)
!
!=======================================================================
!                                                                      !
!  This routine gets the donor grid data (Ac) necessary  to process    !
!  the contact points for a 2D state variable (Ad). It extracts the    !
!  donor cell points containing each contact point, Ac(1:4,:).         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     dg         Donor grid number (integer)                           !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     gtype      C-grid variable type (integer)                        !
!     svname     State variable name (string)                          !
!     cr         Contact region number to process (integer)            !
!     Npoints    Number of points in the contact region (integer)      !
!     contact    Contact region information variables (T_NGC structure)!
!     LBi        Donor grid, I-dimension Lower bound (integer)         !
!     UBi        Donor grid, I-dimension Upper bound (integer)         !
!     LBj        Donor grid, J-dimension Lower bound (integer)         !
!     UBj        Donor grid, J-dimension Upper bound (integer)         !
!     Ad         Donor grid data (2D array)                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ac         2D state variable contact point data                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_nesting

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: dg, model, tile
      integer, intent(in) :: gtype, cr, Npoints
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      character(len=*), intent(in) :: svname
!
      TYPE (T_NGC), intent(in) :: contact(:)
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Ad(LBi:,LBj:)
      real(r8), intent(inout) :: Ac(:,:)
# else
      real(r8), intent(in) :: Ad(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Ac(Npoints,4)
# endif
!
!  Local variable declarations.
!
      integer :: i, ip1, j, jp1, m
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Istr, Iend, Jstr, Jend
# ifdef DISTRIBUTE
      integer :: Npts
# endif

      real(r8), parameter :: Aspv = 0.0_r8
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------

# ifdef DISTRIBUTE
!
!  Initialize contact points array to special value to facilite
!  distribute-memory data collection from all nodes.
!
      Npts=4*Npoints
      DO m=1,Npoints
        Ac(1,m)=Aspv
        Ac(2,m)=Aspv
        Ac(3,m)=Aspv
        Ac(4,m)=Aspv
      END DO
# endif
!
!  Set starting and ending tile indices for the donor grids.
!
      SELECT CASE (gtype)
        CASE (r2dvar)
          Imin=BOUNDS(dg) % IstrT(-1)    ! full RHO-grid range
          Imax=BOUNDS(dg) % IendT(-1)
          Jmin=BOUNDS(dg) % JstrT(-1)
          Jmax=BOUNDS(dg) % JendT(-1)
!
          Istr=BOUNDS(dg) % IstrT(tile)  ! domain partition range
          Iend=BOUNDS(dg) % IendT(tile)
          Jstr=BOUNDS(dg) % JstrT(tile)
          Jend=BOUNDS(dg) % JendT(tile)
        CASE (u2dvar)
          Imin=BOUNDS(dg) % IstrP(-1)    ! full U-grid range
          Imax=BOUNDS(dg) % IendT(-1)
          Jmin=BOUNDS(dg) % JstrT(-1)
          Jmax=BOUNDS(dg) % JendT(-1)
!
          Istr=BOUNDS(dg) % IstrP(tile)  ! domain partition range
          Iend=BOUNDS(dg) % IendT(tile)
          Jstr=BOUNDS(dg) % JstrT(tile)
          Jend=BOUNDS(dg) % JendT(tile)
        CASE (v2dvar)
          Imin=BOUNDS(dg) % IstrT(-1)    ! full V-grid range
          Imax=BOUNDS(dg) % IendT(-1)
          Jmin=BOUNDS(dg) % JstrP(-1)
          Jmax=BOUNDS(dg) % JendT(-1)
!
          Istr=BOUNDS(dg) % IstrT(tile)  ! domain partition range
          Iend=BOUNDS(dg) % IendT(tile)
          Jstr=BOUNDS(dg) % JstrP(tile)
          Jend=BOUNDS(dg) % JendT(tile)
      END SELECT
!
!-----------------------------------------------------------------------
!  Extract donor grid data at contact points.
!-----------------------------------------------------------------------
!
!  Notice that the indices i+1 and j+1 are bounded the maximum values
!  of the grid. This implies that contact point lies on the grid
!  boundary.
!
      DO m=1,Npoints
        i=contact(cr)%Idg(m)
        j=contact(cr)%Jdg(m)
        ip1=MIN(i+1,Imax)
        jp1=MIN(j+1,Jmax)
        IF (((Istr.le.i).and.(i.le.Iend)).and.                          &
     &      ((Jstr.le.j).and.(j.le.Jend))) THEN
          Ac(1,m)=Ad(i  ,j  )
          Ac(2,m)=Ad(ip1,j  )
          Ac(3,m)=Ad(ip1,jp1)
          Ac(4,m)=Ad(i  ,jp1)
        END IF
      END DO

# ifdef DISTRIBUTE
!
!  Gather and broadcast data from all nodes.
!
      CALL mp_assemble (dg, model, Npts, Aspv, Ac)
# endif

      RETURN
      END SUBROUTINE get_contact2d

# ifdef SOLVE3D
!
      SUBROUTINE get_contact3d (dg, model, tile,                        &
     &                          gtype, svname,                          &
     &                          cr, Npoints, contact,                   &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          Ad, Ac)
!
!=======================================================================
!                                                                      !
!  This routine gets the donor grid data (Ac) necessary  to process    !
!  the contact points for a 3D state variable (Ad). It extracts the    !
!  donor cell points containing each contact point, Ac(1:4,k,:).       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     dg         Donor grid number (integer)                           !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     gtype      C-grid variable type (integer)                        !
!     svname     State variable name (string)                          !
!     cr         Contact region number to process (integer)            !
!     Npoints    Number of points in the contact region (integer)      !
!     contact    Contact region information variables (T_NGC structure)!
!     LBi        Donor grid, I-dimension Lower bound (integer)         !
!     UBi        Donor grid, I-dimension Upper bound (integer)         !
!     LBj        Donor grid, J-dimension Lower bound (integer)         !
!     UBj        Donor grid, J-dimension Upper bound (integer)         !
!     Ad         Donor grid data (3D array)                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ac         3D state variable contact point data                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_nesting

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: dg, model, tile
      integer, intent(in) :: gtype, cr, Npoints
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      character(len=*), intent(in) :: svname
!
      TYPE (T_NGC), intent(in) :: contact(:)
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Ad(LBi:,LBj:,LBk:)
      real(r8), intent(inout) :: Ac(:,LBk:,:)
#  else
      real(r8), intent(in) :: Ad(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(inout) :: Ac(4,LBk:UBk,Npoints)
#  endif
!
!  Local variable declarations.
!
      integer :: i, ip1, j, jp1, k, m
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Istr, Iend, Jstr, Jend
#  ifdef DISTRIBUTE
      integer :: Npts
#  endif

      real(r8), parameter :: Aspv = 0.0_r8
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------

#  ifdef DISTRIBUTE
!
!  Initialize contact points array to special value to facilite
!  distribute-memory data collection from all nodes.
!
      Npts=4*(UBk-LBk+1)*Npoints
      DO k=LBk,UBk
        DO m=1,Npoints
          Ac(1,k,m)=Aspv
          Ac(2,k,m)=Aspv
          Ac(3,k,m)=Aspv
          Ac(4,k,m)=Aspv
        END DO
      END DO
#  endif
!
!  Set starting and ending tile indices for the donor grid.
!
      SELECT CASE (gtype)
        CASE (r3dvar)
          Imin=BOUNDS(dg) % IstrT(-1)    ! full RHO-grid range
          Imax=BOUNDS(dg) % IendT(-1)
          Jmin=BOUNDS(dg) % JstrT(-1)
          Jmax=BOUNDS(dg) % JendT(-1)
!
          Istr=BOUNDS(dg) % IstrT(tile)  ! domain partition range
          Iend=BOUNDS(dg) % IendT(tile)
          Jstr=BOUNDS(dg) % JstrT(tile)
          Jend=BOUNDS(dg) % JendT(tile)
        CASE (u3dvar)
          Imin=BOUNDS(dg) % IstrP(-1)    ! full U-grid range
          Imax=BOUNDS(dg) % IendT(-1)
          Jmin=BOUNDS(dg) % JstrT(-1)
          Jmax=BOUNDS(dg) % JendT(-1)
!
          Istr=BOUNDS(dg) % IstrP(tile)  ! domain partition range
          Iend=BOUNDS(dg) % IendT(tile)
          Jstr=BOUNDS(dg) % JstrT(tile)
          Jend=BOUNDS(dg) % JendT(tile)
        CASE (v3dvar)
          Imin=BOUNDS(dg) % IstrT(-1)    ! full V-grid range
          Imax=BOUNDS(dg) % IendT(-1)
          Jmin=BOUNDS(dg) % JstrP(-1)
          Jmax=BOUNDS(dg) % JendT(-1)
!
          Istr=BOUNDS(dg) % IstrT(tile)  ! domain partition range
          Iend=BOUNDS(dg) % IendT(tile)
          Jstr=BOUNDS(dg) % JstrP(tile)
          Jend=BOUNDS(dg) % JendT(tile)
      END SELECT
!
!-----------------------------------------------------------------------
!  Extract donor grid data at contact points.
!-----------------------------------------------------------------------
!
!  Notice that the indices i+1 and j+1 are bounded the maximum values
!  of the grid. This implies that contact point lies on the grid
!  boundary.
!
      DO k=LBk,UBk
        DO m=1,Npoints
          i=contact(cr)%Idg(m)
          j=contact(cr)%Jdg(m)
          ip1=MIN(i+1,Imax)
          jp1=MIN(j+1,Jmax)
          IF (((Istr.le.i).and.(i.le.Iend)).and.                        &
     &        ((Jstr.le.j).and.(j.le.Jend))) THEN
            Ac(1,k,m)=Ad(i  ,j  ,k)
            Ac(2,k,m)=Ad(ip1,j  ,k)
            Ac(3,k,m)=Ad(ip1,jp1,k)
            Ac(4,k,m)=Ad(i  ,jp1,k)
          END IF
        END DO
      END DO

#  ifdef DISTRIBUTE
!
!  Gather and broadcast data from all nodes.
!
      CALL mp_assemble (dg, model, Npts, Aspv, Ac(:,LBk:,:))
#  endif

      RETURN
      END SUBROUTINE get_contact3d
# endif
!
      SUBROUTINE get_persisted2d (dg, rg, model, tile,                  &
     &                            gtype, svname,                        &
     &                            cr, Npoints, contact,                 &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Ad, Ac)
!
!=======================================================================
!                                                                      !
!  This routine gets the donor grid data (Ac) necessary  to process    !
!  the contact points for a 2D flux variable (Ad). It extracts the     !
!  donor cell points containing each contact point, Ac(1:4,:).         !
!                                                                      !
!  This routine is different that 'get_contact2d'.  It is used in      !
!  refinement to impose the appropriate coarser grid flux to insure    !
!  volume and mass conservation.  The value of the coarse grid cell    !
!  is presisted over the refined grid points along its physical        !
!  boundary.  This will facilitate that the sum of all the refined     !
!  grid point is the same as that of the coarse grid containing such   !
!  points.  The spatial interpolation as set in 'get_contact2d' will   !
!  not conserve volume and mass.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     dg         Donor grid number (integer)                           !
!     rg         Receiver grid number (integer)                        !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     gtype      C-grid variable type (integer)                        !
!     svname     State variable name (string)                          !
!     cr         Contact region number to process (integer)            !
!     Npoints    Number of points in the contact region (integer)      !
!     contact    Contact region information variables (T_NGC structure)!
!     LBi        Donor grid, I-dimension Lower bound (integer)         !
!     UBi        Donor grid, I-dimension Upper bound (integer)         !
!     LBj        Donor grid, J-dimension Lower bound (integer)         !
!     UBj        Donor grid, J-dimension Upper bound (integer)         !
!     Ad         Donor grid data (2D array)                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ac         2D flux variable contact point data                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_nesting
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: dg, rg, model, tile
      integer, intent(in) :: gtype, cr, Npoints
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      character(len=*), intent(in) :: svname
!
      TYPE (T_NGC), intent(in) :: contact(:)
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Ad(LBi:,LBj:)
      real(r8), intent(inout) :: Ac(:,:)
# else
      real(r8), intent(in) :: Ad(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Ac(Npoints,4)
# endif
!
!  Local variable declarations.
!
      integer :: Idg, Ip1, Irg, Jdg, Jp1, Jrg
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Istr, Iend, Jstr, Jend
      integer :: i, i_add, j, j_add, m, m_add
# ifdef DISTRIBUTE
      integer :: Npts
# endif

      real(r8), parameter :: Aspv = 0.0_r8
      real(r8):: Rscale
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------

# ifdef DISTRIBUTE
!
!  Initialize contact points array to special value to facilite
!  distribute-memory data collection from all nodes.
!
      Npts=4*Npoints
      DO m=1,Npoints
        Ac(1,m)=Aspv
        Ac(2,m)=Aspv
        Ac(3,m)=Aspv
        Ac(4,m)=Aspv
      END DO
# endif
!
!  Set starting and ending tile indices for the donor grids.
!
      SELECT CASE (gtype)
        CASE (r2dvar)
          Imin=BOUNDS(dg) % IstrT(-1)    ! full RHO-grid range
          Imax=BOUNDS(dg) % IendT(-1)
          Jmin=BOUNDS(dg) % JstrT(-1)
          Jmax=BOUNDS(dg) % JendT(-1)
!
          Istr=BOUNDS(dg) % IstrT(tile)  ! domain partition range
          Iend=BOUNDS(dg) % IendT(tile)
          Jstr=BOUNDS(dg) % JstrT(tile)
          Jend=BOUNDS(dg) % JendT(tile)
!
          m_add=NstrR(cr)-1
        CASE (u2dvar)
          Imin=BOUNDS(dg) % IstrP(-1)    ! full U-grid range
          Imax=BOUNDS(dg) % IendT(-1)
          Jmin=BOUNDS(dg) % JstrT(-1)
          Jmax=BOUNDS(dg) % JendT(-1)
!
          Istr=BOUNDS(dg) % IstrP(tile)  ! domain partition range
          Iend=BOUNDS(dg) % IendT(tile)
          Jstr=BOUNDS(dg) % JstrT(tile)
          Jend=BOUNDS(dg) % JendT(tile)
!
          m_add=NstrU(cr)-1
        CASE (v2dvar)
          Imin=BOUNDS(dg) % IstrT(-1)    ! full V-grid range
          Imax=BOUNDS(dg) % IendT(-1)
          Jmin=BOUNDS(dg) % JstrP(-1)
          Jmax=BOUNDS(dg) % JendT(-1)
!
          Istr=BOUNDS(dg) % IstrT(tile)  ! domain partition range
          Iend=BOUNDS(dg) % IendT(tile)
          Jstr=BOUNDS(dg) % JstrP(tile)
          Jend=BOUNDS(dg) % JendT(tile)
!
          m_add=NstrV(cr)-1
      END SELECT
!
!-----------------------------------------------------------------------
!  Extract donor grid data at contact points.
!-----------------------------------------------------------------------
!
!  Notice that the indices i+1 and j+1 are bounded the maximum values
!  of the grid. This implies that contact point lies on the grid
!  boundary.
!
      Rscale=1.0_r8/REAL(RefineScale(rg))
      DO m=1,Npoints
        Idg=contact(cr)%Idg(m)
        Jdg=contact(cr)%Jdg(m)
        Irg=contact(cr)%Irg(m)
        Jrg=contact(cr)%Jrg(m)
        Ip1=MIN(Idg+1,Imax)
        Jp1=MIN(Jdg+1,Jmax)
        IF (((Istr.le.Idg).and.(Idg.le.Iend)).and.                      &
     &      ((Jstr.le.Jdg).and.(Jdg.le.Jend))) THEN
          IF (on_boundary(m+m_add).gt.0) THEN
            IF ((on_boundary(m+m_add).eq.1).or.                         &
     &          (on_boundary(m+m_add).eq.3)) THEN       ! western and
              j_add=INT(REAL(Jrg-1,r8)*Rscale)          ! eastern edges
              j=J_bottom(rg)+j_add
              Ac(1,m)=Ad(Idg,j)
              Ac(2,m)=Ad(Idg,j)
              Ac(3,m)=Ad(Idg,j)
              Ac(4,m)=Ad(Idg,j)
            ELSE IF ((on_boundary(m+m_add).eq.2).or.                    &
     &               (on_boundary(m+m_add).eq.4)) THEN  ! southern and
              i_add=INT(REAL(Irg-1,r8)*Rscale)          ! northern edges
              i=I_left(rg)+i_add
              Ac(1,m)=Ad(i,Jdg)
              Ac(2,m)=Ad(i,Jdg)
              Ac(3,m)=Ad(i,Jdg)
              Ac(4,m)=Ad(i,Jdg)
            END IF
          ELSE
            Ac(1,m)=Ad(Idg,Jdg)               ! contact point is not at
            Ac(2,m)=Ad(Ip1,Jdg)               ! physical boundary, just
            Ac(3,m)=Ad(Ip1,Jp1)               ! set values for spatial
            Ac(4,m)=Ad(Idg,Jp1)               ! interpolation (not used)
          END IF
        END IF
      END DO

# ifdef DISTRIBUTE
!
!  Gather and broadcast data from all nodes.
!
      CALL mp_assemble (dg, model, Npts, Aspv, Ac)
      IF (exit_flag.ne.NoError) RETURN
# endif

      RETURN
      END SUBROUTINE get_persisted2d
!
      SUBROUTINE put_contact2d (rg, model, tile,                        &
     &                          gtype, svname,                          &
     &                          cr, Npoints, contact,                   &
     &                          LBi, UBi, LBj, UBj,                     &
# ifdef MASKING
     &                          Amask,                                  &
# endif
     &                          Ac, Ar)
!
!=======================================================================
!                                                                      !
!  This routine uses extracted donor grid data (Ac) to spatially       !
!  interpolate a 2D state variable  at the receiver grid contact       !
!  points.  If the donor and receiver grids are coincident,  the       !
!  Lweight(1,:) is unity and Lweight(2:4,:) are zero.                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     rg         Receiver grid number (integer)                        !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     gtype      C-grid variable type (integer)                        !
!     svname     State variable name (string)                          !
!     cr         Contact region number to process (integer)            !
!     Npoints    Number of points in the contact region (integer)      !
!     contact    Contact region information variables (T_NGC structure)!
!     LBi        Receiver grid, I-dimension Lower bound (integer)      !
!     UBi        Receiver grid, I-dimension Upper bound (integer)      !
!     LBj        Receiver grid, J-dimension Lower bound (integer)      !
!     UBj        Receiver grid, J-dimension Upper bound (integer)      !
# ifdef MASKING
!     Amask      Receiver grid land/sea masking                        !
# endif
!     Ac         Contact point data extracted from donor grid          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ar         Updated receiver grid 2D state array                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_nesting
!
!  Imported variable declarations.
!
      integer, intent(in) :: rg, model, tile
      integer, intent(in) :: gtype, cr, Npoints
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      character(len=*), intent(in) :: svname
!
      TYPE (T_NGC), intent(in) :: contact(:)
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Ac(:,:)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: Ar(LBi:,LBj:)
# else
      real(r8), intent(in) :: Ac(4,Npoints)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: Ar(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, m
      integer :: Istr, Iend, Jstr, Jend
!
!-----------------------------------------------------------------------
!  Interpolate 2D data from donor grid to receiver grid contact points.
!-----------------------------------------------------------------------
!
!  Set starting and ending tile indices for the receiver grid.
!
      SELECT CASE (gtype)
        CASE (r2dvar)
          Istr=BOUNDS(rg) % IstrT(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrT(tile)
          Jend=BOUNDS(rg) % JendT(tile)
        CASE (u2dvar)
          Istr=BOUNDS(rg) % IstrP(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrT(tile)
          Jend=BOUNDS(rg) % JendT(tile)
        CASE (v2dvar)
          Istr=BOUNDS(rg) % IstrT(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrP(tile)
          Jend=BOUNDS(rg) % JendT(tile)
      END SELECT
!
!  Interpolate.
!
      DO m=1,Npoints
        i=contact(cr)%Irg(m)
        j=contact(cr)%Jrg(m)
        IF (((Istr.le.i).and.(i.le.Iend)).and.                          &
     &      ((Jstr.le.j).and.(j.le.Jend))) THEN
          Ar(i,j)=contact(cr)%Lweight(1,m)*Ac(1,m)+                     &
     &            contact(cr)%Lweight(2,m)*Ac(2,m)+                     &
     &            contact(cr)%Lweight(3,m)*Ac(3,m)+                     &
     &            contact(cr)%Lweight(4,m)*Ac(4,m)
# ifdef MASKING
          Ar(i,j)=Ar(i,j)*Amask(i,j)
# endif
        END IF
      END DO

      RETURN
      END SUBROUTINE put_contact2d

# ifdef SOLVE3D
!
      SUBROUTINE put_contact3d (rg, model, tile,                        &
     &                          gtype, svname,                          &
     &                          cr, Npoints, contact,                   &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
#  ifdef MASKING
     &                          Amask,                                  &
#  endif
     &                          Ac, Ar)
!
!=======================================================================
!                                                                      !
!  This routine uses extracted donor grid data (Ac) to spatially       !
!  interpolate a 3D state variable  at the receiver grid contact       !
!  points.  If the donor and receiver grids  are concident,  the       !
!  Lweight(1,:) is unity and Lweight(2:4,:) are zero.                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     rg         Receiver grid number (integer)                        !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     gtype      C-grid variable type (integer)                        !
!     svname     State variable name (string)                          !
!     cr         Contact region number to process (integer)            !
!     Npoints    Number of points in the contact region (integer)      !
!     contact    Contact region information variables (T_NGC structure)!
!     LBi        Receiver grid, I-dimension Lower bound (integer)      !
!     UBi        Receiver grid, I-dimension Upper bound (integer)      !
!     LBj        Receiver grid, J-dimension Lower bound (integer)      !
!     UBj        Receiver grid, J-dimension Upper bound (integer)      !
!     LBk        Receiver grid, K-dimension Lower bound (integer)      !
!     UBk        Receiver grid, K-dimension Upper bound (integer)      !
#  ifdef MASKING
!     Amask      Receiver grid land/sea masking                        !
#  endif
!     Ac         Contact point data extracted from donor grid          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ar         Updated receiver grid 3D state array                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_nesting
!
!  Imported variable declarations.
!
      integer, intent(in) :: rg, model, tile
      integer, intent(in) :: gtype, cr, Npoints
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      character(len=*), intent(in) :: svname
!
      TYPE (T_NGC), intent(in) :: contact(:)
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Ac(:,:,:)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: Ar(LBi:,LBj:,LBk:)
#  else
      real(r8), intent(in) :: Ac(Npoints,LBk:UBk,4)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: Ar(LBi:UBi,LBj:UBj,LBk:UBk)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k, kdg, kdgm1, m
      integer :: Istr, Iend, Jstr, Jend, Kmin

      real(r8), dimension(8) :: cff
!
!-----------------------------------------------------------------------
!  Interpolate 3D data from donor grid to receiver grid contact points.
!-----------------------------------------------------------------------
!
!  Set starting and ending tile indices for the receiver grid.
!
      SELECT CASE (gtype)
        CASE (r3dvar)
          Istr=BOUNDS(rg) % IstrT(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrT(tile)
          Jend=BOUNDS(rg) % JendT(tile)
          Kmin=1
        CASE (u3dvar)
          Istr=BOUNDS(rg) % IstrP(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrT(tile)
          Jend=BOUNDS(rg) % JendT(tile)
          Kmin=1
        CASE (v3dvar)
          Istr=BOUNDS(rg) % IstrT(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrP(tile)
          Jend=BOUNDS(rg) % JendT(tile)
          Kmin=1
        CASE (w3dvar)
          Istr=BOUNDS(rg) % IstrT(tile)
          Iend=BOUNDS(rg) % IendT(tile)
          Jstr=BOUNDS(rg) % JstrT(tile)
          Jend=BOUNDS(rg) % JendT(tile)
          Kmin=0
      END SELECT
!
!  Interpolate.
!
      DO k=LBk,UBk
        DO m=1,Npoints
          i=contact(cr)%Irg(m)
          j=contact(cr)%Jrg(m)
          kdg=contact(cr)%Kdg(k,m)
          kdgm1=MAX(kdg-1,Kmin)
          IF (((Istr.le.i).and.(i.le.Iend)).and.                        &
     &        ((Jstr.le.j).and.(j.le.Jend))) THEN
            cff(1)=contact(cr)%Lweight(1,m)*contact(cr)%Vweight(1,k,m)
            cff(2)=contact(cr)%Lweight(2,m)*contact(cr)%Vweight(1,k,m)
            cff(3)=contact(cr)%Lweight(3,m)*contact(cr)%Vweight(1,k,m)
            cff(4)=contact(cr)%Lweight(4,m)*contact(cr)%Vweight(1,k,m)
            cff(5)=contact(cr)%Lweight(1,m)*contact(cr)%Vweight(2,k,m)
            cff(6)=contact(cr)%Lweight(2,m)*contact(cr)%Vweight(2,k,m)
            cff(7)=contact(cr)%Lweight(3,m)*contact(cr)%Vweight(2,k,m)
            cff(8)=contact(cr)%Lweight(4,m)*contact(cr)%Vweight(2,k,m)
            Ar(i,j,k)=cff(1)*Ac(1,kdgm1,m)+                             &
     &                cff(2)*Ac(2,kdgm1,m)+                             &
     &                cff(3)*Ac(3,kdgm1,m)+                             &
     &                cff(4)*Ac(4,kdgm1,m)+                             &
     &                cff(5)*Ac(1,kdg  ,m)+                             &
     &                cff(6)*Ac(2,kdg  ,m)+                             &
     &                cff(7)*Ac(3,kdg  ,m)+                             &
     &                cff(8)*Ac(4,kdg  ,m)
#  ifdef MASKING
            Ar(i,j,k)=Ar(i,j,k)*Amask(i,j)
#  endif
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE put_contact3d
# endif
!
      SUBROUTINE put_refine2d (ng, dg, cr, model, tile, LputFsur,       &
     &                         LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine interpolates (space, time) refinement grid 2D state    !
!  variables contact points using data from the donor grid.  Notice    !
!  that because of shared-memory parallelism,  the  free-surface is    !
!  processed first and in a different parallel region.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Refinement (receiver) grid number (integer)           !
!     dg         Donor grid number (integer)                           !
!     cr         Contact region number to process (integer)            !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     LputFsur   Switch to process or not free-surface (logical)       !
!     LBi        Receiver grid, I-dimension Lower bound (integer)      !
!     UBi        Receiver grid, I-dimension Upper bound (integer)      !
!     LBj        Receiver grid, J-dimension Lower bound (integer)      !
!     UBj        Receiver grid, J-dimension Upper bound (integer)      !
!                                                                      !
!  On Output:    OCEAN(ng) structure                                   !
!                                                                      !
!     zeta       Updated free-surface                                  !
!     ubar       Updated 2D momentum in the XI-direction               !
!     vbar       Updated 2D momentum in the ETA-direction              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
      USE mod_grid
      USE mod_nesting
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: LputFsur
      integer, intent(in) :: ng, dg, cr, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      logical :: Uboundary, Vboundary

# ifdef DISTRIBUTE
      integer :: ILB, IUB, JLB, JUB, NptsSN, NptsWE, my_tile
# endif
      integer :: NSUB, i, irec, j, m, tnew, told
      integer :: Idg, Jdg

# ifdef DISTRIBUTE
      real(r8), parameter :: spv = 0.0_r8
# endif
      real(r8) :: Wnew, Wold
      real(r8) :: cff, cff1, fac
      real(r8) :: my_value

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Interpolate (space, time) refinement grid contact points for 2D state
!  variables from donor grid.
!-----------------------------------------------------------------------

# ifdef DISTRIBUTE
!
!  Set global size of boundary edges.
!
      IF (.not.LputFsur) THEN
        my_tile=-1
        ILB=BOUNDS(ng)%LBi(my_tile)
        IUB=BOUNDS(ng)%UBi(my_tile)
        JLB=BOUNDS(ng)%LBj(my_tile)
        JUB=BOUNDS(ng)%UBj(my_tile)
        NptsWE=JUB-JLB+1
        NptsSN=IUB-ILB+1
!
!  If distributed-memory, initialize arrays used to check mass flux
!  conservation with special value (zero) to facilitate the global
!  reduction when collecting data between all nodes.
!
        BRY_CONTACT(iwest ,cr)%Mflux=spv
        BRY_CONTACT(ieast ,cr)%Mflux=spv
        BRY_CONTACT(isouth,cr)%Mflux=spv
        BRY_CONTACT(inorth,cr)%Mflux=spv
      END IF
# endif
!
!  Set time snapshot indices for the donor grid data.
!
      told=3-RollingIndex(cr)
      tnew=RollingIndex(cr)
!
!  Set time interpolation weights.
!
      Wold=RollingTime(tnew,cr)-time(ng)
      Wnew=time(ng)-RollingTime(told,cr)
      fac=1.0_r8/(Wold+Wnew)
      Wold=fac*Wold
      Wnew=fac*Wnew
!
      IF (((Wold*Wnew).lt.0.0_r8).or.((Wold+Wnew).le.0.0_r8)) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) cr, dg, ng,                               &
     &                        iic(dg), told, tnew,                      &
     &                        iic(ng), Wold, Wnew,                      &
     &                        INT(time(ng)),                            &
     &                        INT(RollingTime(told,cr)),                &
     &                        INT(RollingTime(tnew,cr))
          END IF
          exit_flag=8
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Process free-surface.
!-----------------------------------------------------------------------
!
      FREE_SURFACE : IF (LputFsur) THEN
        DO m=1,Rcontact(cr)%Npoints
          i=Rcontact(cr)%Irg(m)
          j=Rcontact(cr)%Jrg(m)
          IF (((IstrT.le.i).and.(i.le.IendT)).and.                      &
     &        ((JstrT.le.j).and.(j.le.JendT))) THEN
            my_value=Wold*(Rcontact(cr)%Lweight(1,m)*                   &
     &                                  REFINED(cr)%zeta(1,m,told)+     &
     &                     Rcontact(cr)%Lweight(2,m)*                   &
     &                                  REFINED(cr)%zeta(2,m,told)+     &
     &                     Rcontact(cr)%Lweight(3,m)*                   &
     &                                  REFINED(cr)%zeta(3,m,told)+     &
     &                     Rcontact(cr)%Lweight(4,m)*                   &
     &                                  REFINED(cr)%zeta(4,m,told))+    &
     &               Wnew*(Rcontact(cr)%Lweight(1,m)*                   &
     &                                  REFINED(cr)%zeta(1,m,tnew)+     &
     &                     Rcontact(cr)%Lweight(2,m)*                   &
     &                                  REFINED(cr)%zeta(2,m,tnew)+     &
     &                     Rcontact(cr)%Lweight(3,m)*                   &
     &                                  REFINED(cr)%zeta(3,m,tnew)+     &
     &                     Rcontact(cr)%Lweight(4,m)*                   &
     &                                  REFINED(cr)%zeta(4,m,tnew))
# ifdef MASKING
            my_value=my_value*GRID(ng)%rmask(i,j)
# endif
# ifdef WET_DRY
            IF (my_value.le.(Dcrit(ng)-GRID(ng)%h(i,j))) THEN
              my_value=Dcrit(ng)-GRID(ng)%h(i,j)
            END IF
# endif
            OCEAN(ng)%zeta(i,j,1)=my_value
            OCEAN(ng)%zeta(i,j,2)=my_value
            OCEAN(ng)%zeta(i,j,3)=my_value
# ifdef SOLVE3D
            COUPLING(ng)%Zt_avg1(i,j)=my_value
# endif
          END IF
        END DO

      ELSE
!
!-----------------------------------------------------------------------
!  Process 2D momentum.
!-----------------------------------------------------------------------
!
!  2D momentum in the XI-direction.
# ifdef SOLVE3D
!
!  Notice that contact points at the domain western and eastern
!  boundaries are avoided for indx1(ng) time record. They are be
!  assigned in the mass flux computations below. This exception is
!  done for adjoint correctness.
# endif
!
        DO m=1,Ucontact(cr)%Npoints
          i=Ucontact(cr)%Irg(m)
          j=Ucontact(cr)%Jrg(m)
          IF (((IstrP.le.i).and.(i.le.IendT)).and.                      &
     &        ((JstrT.le.j).and.(j.le.JendT))) THEN
            my_value=Wold*(Ucontact(cr)%Lweight(1,m)*                   &
     &                                  REFINED(cr)%ubar(1,m,told)+     &
     &                     Ucontact(cr)%Lweight(2,m)*                   &
     &                                  REFINED(cr)%ubar(2,m,told)+     &
     &                     Ucontact(cr)%Lweight(3,m)*                   &
     &                                  REFINED(cr)%ubar(3,m,told)+     &
     &                     Ucontact(cr)%Lweight(4,m)*                   &
     &                                  REFINED(cr)%ubar(4,m,told))+    &
     &               Wnew*(Ucontact(cr)%Lweight(1,m)*                   &
     &                                  REFINED(cr)%ubar(1,m,tnew)+     &
     &                     Ucontact(cr)%Lweight(2,m)*                   &
     &                                  REFINED(cr)%ubar(2,m,tnew)+     &
     &                     Ucontact(cr)%Lweight(3,m)*                   &
     &                                  REFINED(cr)%ubar(3,m,tnew)+     &
     &                     Ucontact(cr)%Lweight(4,m)*                   &
     &                                  REFINED(cr)%ubar(4,m,tnew))
# ifdef MASKING
            my_value=my_value*GRID(ng)%umask(i,j)
# endif
# ifdef WET_DRY
            my_value=my_value*GRID(ng)%umask_wet(i,j)
# endif
            DO irec=1,3
# ifdef SOLVE3D
              Uboundary=(m.eq.BRY_CONTACT(iwest,cr)%C2Bindex(j)).or.    &
     &                  (m.eq.BRY_CONTACT(ieast,cr)%C2Bindex(j))
              IF(.not.(Uboundary.and.(irec.eq.indx1(ng)))) THEN
                OCEAN(ng)%ubar(i,j,irec)=my_value
!!            ELSE                                   ! for debugging
!!              OCEAN(ng)%ubar(i,j,irec)=0.0_r8      ! purposes
              END IF
# else
              OCEAN(ng)%ubar(i,j,irec)=my_value
# endif
            END DO
          END IF
        END DO
!
!  2D momentum in the ETA-direction.
# ifdef SOLVE3D
!
!  Notice that contact points at the domain southern and northern
!  boundaries are avoided for indx1(ng) time record. They are be
!  assigned in the mass flux computations below. This exception is
!  done for adjoint correctness.
# endif
!
        DO m=1,Vcontact(cr)%Npoints
          i=Vcontact(cr)%Irg(m)
          j=Vcontact(cr)%Jrg(m)
          IF (((IstrT.le.i).and.(i.le.IendT)).and.                      &
     &        ((JstrP.le.j).and.(j.le.JendT))) THEN
            my_value=Wold*(Vcontact(cr)%Lweight(1,m)*                   &
     &                                  REFINED(cr)%vbar(1,m,told)+     &
     &                     Vcontact(cr)%Lweight(2,m)*                   &
     &                                  REFINED(cr)%vbar(2,m,told)+     &
     &                     Vcontact(cr)%Lweight(3,m)*                   &
     &                                  REFINED(cr)%vbar(3,m,told)+     &
     &                     Vcontact(cr)%Lweight(4,m)*                   &
     &                                  REFINED(cr)%vbar(4,m,told))+    &
     &               Wnew*(Vcontact(cr)%Lweight(1,m)*                   &
     &                                  REFINED(cr)%vbar(1,m,tnew)+     &
     &                     Vcontact(cr)%Lweight(2,m)*                   &
     &                                  REFINED(cr)%vbar(2,m,tnew)+     &
     &                     Vcontact(cr)%Lweight(3,m)*                   &
     &                                  REFINED(cr)%vbar(3,m,tnew)+     &
     &                     Vcontact(cr)%Lweight(4,m)*                   &
     &                                  REFINED(cr)%vbar(4,m,tnew))
# ifdef MASKING
            my_value=my_value*GRID(ng)%vmask(i,j)
# endif
# ifdef WET_DRY
            my_value=my_value*GRID(ng)%vmask_wet(i,j)
# endif
            DO irec=1,3
# ifdef SOLVE3D
              Vboundary=(m.eq.BRY_CONTACT(isouth,cr)%C2Bindex(i)).or.   &
     &                  (m.eq.BRY_CONTACT(inorth,cr)%C2Bindex(i))
              IF(.not.(Vboundary.and.(irec.eq.indx1(ng)))) THEN
                OCEAN(ng)%vbar(i,j,irec)=my_value
!!            ELSE                                   ! for debugging
!!              OCEAN(ng)%vbar(i,j,irec)=0.0_r8      ! purposes
              END IF
# else
              OCEAN(ng)%vbar(i,j,irec)=my_value
# endif
            END DO
          END IF
        END DO

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Impose mass flux at the finer grid physical boundaries. This is only
!  done for indx1(ng) time record.
!
!  Western/Eastern boundary:
!
!    ubar(Ibry,:,indx1) = DU_avg2(Ibry,:) * pn(Ibry,:) / D(Ibry,:)
!
!  Southern/Northern boundary:
!
!    vbar(:,Jbry,indx1) = DV_avg2(:,Jbry) * pm(:,Jbry) / D(:,Jbry)
!
!  We use the latest coarse grid mass flux REFINED(cr)%DU_avg(1,:,tnew)
!  with a linear variation (cff1) to ensure that the sum of the refined
!  grid fluxes equals the coarse grid flux.
!-----------------------------------------------------------------------
!
!  Western edge.
!
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            m=BRY_CONTACT(iwest,cr)%C2Bindex(j)
            Idg=Ucontact(cr)%Idg(m)                 ! for debugging
            Jdg=Ucontact(cr)%Jdg(m)                 ! purposes
            cff=0.5_r8*GRID(ng)%on_u(Istr,j)*                           &
                (GRID(ng)%h(Istr-1,j)+                                  &
     &           OCEAN(ng)%zeta(Istr-1,j,indx1(ng))+                    &
     &           GRID(ng)%h(Istr  ,j)+                                  &
     &           OCEAN(ng)%zeta(Istr  ,j,indx1(ng)))
            cff1=GRID(ng)%on_u(Istr,j)/REFINED(cr)%on_u(m)
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DU_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DU_avg2(1,m,tnew))/cff
#  else
            my_value=cff1*REFINED(cr)%DU_avg2(1,m,tnew)/cff
#  endif
#  ifdef WEC 
            my_value=my_value-OCEAN(ng)%ubar_stokes(Istr,j) 
#  endif 
#  ifdef MASKING
            my_value=my_value*GRID(ng)%umask(Istr,j)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%umask_wet(Istr,j)
#  endif
            BRY_CONTACT(iwest,cr)%Mflux(j)=cff*my_value
            OCEAN(ng)%ubar(Istr,j,indx1(ng))=my_value
          END DO
        END IF
!
!  Eastern edge.
!
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            m=BRY_CONTACT(ieast,cr)%C2Bindex(j)
            Idg=Ucontact(cr)%Idg(m)                 ! for debugging
            Jdg=Ucontact(cr)%Jdg(m)                 ! purposes
            cff=0.5_r8*GRID(ng)%on_u(Iend+1,j)*                         &
     &          (GRID(ng)%h(Iend+1,j)+                                  &
     &           OCEAN(ng)%zeta(Iend+1,j,indx1(ng))+                    &
     &           GRID(ng)%h(Iend  ,j)+                                  &
     &           OCEAN(ng)%zeta(Iend  ,j,indx1(ng)))
            cff1=GRID(ng)%on_u(Iend+1,j)/REFINED(cr)%on_u(m)
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DU_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DU_avg2(1,m,tnew))/cff
#  else
            my_value=cff1*REFINED(cr)%DU_avg2(1,m,tnew)/cff
#  endif
#  ifdef WEC 
            my_value=my_value-OCEAN(ng)%ubar_stokes(Iend+1,j) 
#  endif 
#  ifdef MASKING
            my_value=my_value*GRID(ng)%umask(Iend+1,j)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%umask_wet(Iend+1,j)
#  endif
            BRY_CONTACT(ieast,cr)%Mflux(j)=cff*my_value
            OCEAN(ng)%ubar(Iend+1,j,indx1(ng))=my_value
          END DO
        END IF
!
!  Southern edge.
!
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            m=BRY_CONTACT(isouth,cr)%C2Bindex(i)
            Idg=Vcontact(cr)%Idg(m)                 ! for debugging
            Jdg=Vcontact(cr)%Jdg(m)                 ! purposes
            cff=0.5_r8*GRID(ng)%om_v(i,Jstr)*                           &
     &          (GRID(ng)%h(i,Jstr-1)+                                  &
     &           OCEAN(ng)%zeta(i,Jstr-1,indx1(ng))+                    &
     &           GRID(ng)%h(i,Jstr  )+                                  &
     &           OCEAN(ng)%zeta(i,Jstr  ,indx1(ng)))
            cff1=GRID(ng)%om_v(i,Jstr)/REFINED(cr)%om_v(m)
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DV_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DV_avg2(1,m,tnew))/cff
#  else
            my_value=cff1*REFINED(cr)%DV_avg2(1,m,tnew)/cff
#  endif
#  ifdef WEC 
            my_value=my_value-OCEAN(ng)%vbar_stokes(i,Jstr) 
#  endif 
#  ifdef MASKING
            my_value=my_value*GRID(ng)%vmask(i,Jstr)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%vmask_wet(i,Jstr)
#  endif
            BRY_CONTACT(isouth,cr)%Mflux(i)=cff*my_value
            OCEAN(ng)%vbar(i,Jstr,indx1(ng))=my_value
          END DO
        END IF
!
!  Northern edge.
!
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            m=BRY_CONTACT(inorth,cr)%C2Bindex(i)
            Idg=Vcontact(cr)%Idg(m)                 ! for debugging
            Jdg=Vcontact(cr)%Jdg(m)                 ! purposes
            cff=0.5_r8*GRID(ng)%om_v(i,Jend+1)*                         &
     &          (GRID(ng)%h(i,Jend+1)+                                  &
     &           OCEAN(ng)%zeta(i,Jend+1,indx1(ng))+                    &
     &           GRID(ng)%h(i,Jend  )+                                  &
     &           OCEAN(ng)%zeta(i,Jend  ,indx1(ng)))
            cff1=GRID(ng)%om_v(i,Jend+1)/REFINED(cr)%om_v(m)
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DV_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DV_avg2(1,m,tnew))/cff
#  else
            my_value=cff1*REFINED(cr)%DV_avg2(1,m,tnew)/cff
#  endif
#  ifdef WEC 
            my_value=my_value-OCEAN(ng)%vbar_stokes(i,Jend+1) 
#  endif 
#  ifdef MASKING
            my_value=my_value*GRID(ng)%vmask(i,Jend+1)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%vmask_wet(i,Jend+1)
#  endif
            BRY_CONTACT(inorth,cr)%Mflux(i)=cff*my_value
            OCEAN(ng)%vbar(i,Jend+1,indx1(ng))=my_value
          END DO
        END IF
# endif
      END IF FREE_SURFACE

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile information.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      IF (LputFsur) THEN
        CALL mp_exchange2d (ng, tile, model,                            &
#  ifdef SOLVE3D
     &                      4,                                          &
#  else
     &                      3,                                          &
#  endif
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
#  ifdef SOLVE3D
     &                      COUPLING(ng)%Zt_avg1,                       &
#  endif
     &                      OCEAN(ng)%zeta(:,:,1),                      &
     &                      OCEAN(ng)%zeta(:,:,2),                      &
     &                      OCEAN(ng)%zeta(:,:,3))
!
!  2D momentum.
!
      ELSE
        CALL mp_exchange2d (ng, tile, model, 3,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      OCEAN(ng)%ubar(:,:,1),                      &
     &                      OCEAN(ng)%ubar(:,:,2),                      &
     &                      OCEAN(ng)%ubar(:,:,3))

        CALL mp_exchange2d (ng, tile, model, 3,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      OCEAN(ng)%vbar(:,:,1),                      &
     &                      OCEAN(ng)%vbar(:,:,2),                      &
     &                      OCEAN(ng)%vbar(:,:,3))
!
        CALL mp_assemble (ng, model, NptsWE, spv,                        &
     &                    BRY_CONTACT(iwest ,cr)%Mflux(JLB:))
        CALL mp_assemble (ng, model, NptsWE, spv,                        &
     &                    BRY_CONTACT(ieast ,cr)%Mflux(JLB:))
        CALL mp_assemble (ng, model, NptsSN, spv,                        &
     &                    BRY_CONTACT(isouth,cr)%Mflux(ILB:))
        CALL mp_assemble (ng, model, NptsSN, spv,                        &
     &                    BRY_CONTACT(inorth,cr)%Mflux(ILB:))
      END IF
# endif
!
  10  FORMAT (/,' PUT_REFINE2D - unbounded contact points temporal: ',  &
     &          ' interpolation:',                                      &
     &        /,2x,  'cr       = ',i2.2,                                &
     &          8x,'dg         = ',i2.2,                                &
     &          8x,'ng         = ',i2.2,                                &
     &        /,2x,  'iic(dg)  = ',i7.7,                                &
     &          3x,'told       = ',i1,                                  &
     &          9x,'tnew       = ',i1,                                  &
     &        /,2x,  'iic(ng)  = ',i7.7,                                &
     &          3x,'Wold       = ',f8.5,                                &
     &          2x,'Wnew       = ',f8.5,                                &
     &        /,2x,  'time(ng) = ',i7.7,                                &
     &          3x,'time(told) = ',i7.7,                                &
     &          3x,'time(tnew) = ',i7.7)

      RETURN
      END SUBROUTINE put_refine2d

# ifdef SOLVE3D
!
      SUBROUTINE put_refine3d (ng, dg, cr, model, tile,                 &
     &                         LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine interpolates (space, time) refinement grid 3D state    !
!  variables contact points using data from the donor grid.            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Refinement (receiver) grid number (integer)           !
!     dg         Donor grid number (integer)                           !
!     cr         Contact region number to process (integer)            !
!     model      Calling model identifier (integer)                    !
!     tile       Domain tile partition (integer)                       !
!     LBi        Receiver grid, I-dimension Lower bound (integer)      !
!     UBi        Receiver grid, I-dimension Upper bound (integer)      !
!     LBj        Receiver grid, J-dimension Lower bound (integer)      !
!     UBj        Receiver grid, J-dimension Upper bound (integer)      !
!                                                                      !
!  On Output:    OCEAN(ng) structure                                   !
!                                                                      !
!     t          Updated tracer-type variables                         !
!     u          Updated 3D momentum in the XI-direction               !
!     v          Updated 3D momentum in the ETA-direction              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_nesting
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits

#  ifdef DISTRIBUTE
!
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, dg, cr, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
# ifdef NESTING_DEBUG
      logical, save :: first = .TRUE.
# endif
      integer :: i, itrc, j, k, m, tnew, told

      real(r8) :: Wnew, Wold, fac
      real(r8) :: my_value

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Interpolate (space, time) refinement grid contact points for 2D state
!  variables from donor grid.
!-----------------------------------------------------------------------
!
!  Set time snapshot indices for the donor grid data.
!
      told=3-RollingIndex(cr)
      tnew=RollingIndex(cr)
!
!  Set time interpolation weights.
!
      Wold=RollingTime(tnew,cr)-time(ng)
      Wnew=time(ng)-RollingTime(told,cr)
      fac=1.0_r8/(Wold+Wnew)
      Wold=fac*Wold
      Wnew=fac*Wnew
!
      IF (((Wold*Wnew).lt.0.0_r8).or.((Wold+Wnew).le.0.0_r8)) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) cr, dg, ng,                               &
     &                        iic(dg), told, tnew,                      &
     &                        iic(ng), Wold, Wnew,                      &
     &                        INT(time(ng)),                            &
     &                        INT(RollingTime(told,cr)),                &
     &                        INT(RollingTime(tnew,cr))
          END IF
          exit_flag=8
        END IF
      END IF

# ifdef NESTING_DEBUG
!
!  If debugging, write information into Fortran unit 200 to check the
!  logic of interpolating from donor grid data.
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        IF (Master) THEN
          IF (first) THEN
            first=.FALSE.
            WRITE (200,20)
          END IF
          WRITE (200,30) cr, dg, ng, iic(dg), iic(ng), told,  tnew,     &
     &                   INT(time(ng)),                                 &
     &                   INT(RollingTime(told,cr)),                     &
     &                   INT(RollingTime(tnew,cr)),                     &
     &                   Wold, Wnew
 20       FORMAT (3x,'cr',3x,'dg',3x,'ng',3x,'iic',2x,'iic',1x,'told',  &
     &            1x,'tnew',3x,'time',3x,'time',3x,'time',5x,'Wold',    &
     &            7x,'Wnew',/,17x,'(dg)',1x,'(ng)',13x,'(ng)',3x,       &
     &            'told',3x,'tnew',/)
 30       FORMAT (7i5,3i7,2f11.4)
          CALL my_flush (200)
        END IF
      END IF
# endif
!
!  Tracer-type variables.
!
      DO m=1,Rcontact(cr)%Npoints
        i=Rcontact(cr)%Irg(m)
        j=Rcontact(cr)%Jrg(m)
        IF (((IstrT.le.i).and.(i.le.IendT)).and.                        &
     &      ((JstrT.le.j).and.(j.le.JendT))) THEN
          DO itrc=1,NT(ng)
            DO k=1,N(ng)
              my_value=Wold*                                            &
     &                 (Rcontact(cr)%Lweight(1,m)*                      &
     &                               REFINED(cr)%t(1,k,m,told,itrc)+    &
     &                  Rcontact(cr)%Lweight(2,m)*                      &
     &                               REFINED(cr)%t(2,k,m,told,itrc)+    &
     &                  Rcontact(cr)%Lweight(3,m)*                      &
     &                               REFINED(cr)%t(3,k,m,told,itrc)+    &
     &                  Rcontact(cr)%Lweight(4,m)*                      &
     &                               REFINED(cr)%t(4,k,m,told,itrc))+   &
     &                 Wnew*                                            &
     &                 (Rcontact(cr)%Lweight(1,m)*                      &
     &                               REFINED(cr)%t(1,k,m,tnew,itrc)+    &
     &                  Rcontact(cr)%Lweight(2,m)*                      &
     &                               REFINED(cr)%t(2,k,m,tnew,itrc)+    &
     &                  Rcontact(cr)%Lweight(3,m)*                      &
     &                               REFINED(cr)%t(3,k,m,tnew,itrc)+    &
     &                  Rcontact(cr)%Lweight(4,m)*                      &
     &                               REFINED(cr)%t(4,k,m,tnew,itrc))
#  ifdef MASKING
              my_value=my_value*GRID(ng)%rmask(i,j)
#  endif
              OCEAN(ng)%t(i,j,k,1,itrc)=my_value
              OCEAN(ng)%t(i,j,k,2,itrc)=my_value
              OCEAN(ng)%t(i,j,k,3,itrc)=my_value
            END DO
          END DO
        END IF
      END DO
!
!  3D momentum in the XI-direction.
!
      DO m=1,Ucontact(cr)%Npoints
        i=Ucontact(cr)%Irg(m)
        j=Ucontact(cr)%Jrg(m)
        IF (((IstrP.le.i).and.(i.le.IendT)).and.                        &
     &      ((JstrT.le.j).and.(j.le.JendT))) THEN
          DO k=1,N(ng)
            my_value=Wold*                                              &
     &               (Ucontact(cr)%Lweight(1,m)*                        &
     &                             REFINED(cr)%u(1,k,m,told)+           &
     &                Ucontact(cr)%Lweight(2,m)*                        &
     &                             REFINED(cr)%u(2,k,m,told)+           &
     &                Ucontact(cr)%Lweight(3,m)*                        &
     &                             REFINED(cr)%u(3,k,m,told)+           &
     &                Ucontact(cr)%Lweight(4,m)*                        &
     &                             REFINED(cr)%u(4,k,m,told))+          &
     &               Wnew*                                              &
     &               (Ucontact(cr)%Lweight(1,m)*                        &
     &                             REFINED(cr)%u(1,k,m,tnew)+           &
     &                Ucontact(cr)%Lweight(2,m)*                        &
     &                             REFINED(cr)%u(2,k,m,tnew)+           &
     &                Ucontact(cr)%Lweight(3,m)*                        &
     &                             REFINED(cr)%u(3,k,m,tnew)+           &
     &                Ucontact(cr)%Lweight(4,m)*                        &
     &                             REFINED(cr)%u(4,k,m,tnew))
#  ifdef MASKING
            my_value=my_value*GRID(ng)%umask(i,j)
#  endif
            OCEAN(ng)%u(i,j,k,1)=my_value
            OCEAN(ng)%u(i,j,k,2)=my_value
          END DO
        END IF
      END DO
!
!  3D momentum in the ETA-direction.
!
      DO m=1,Vcontact(cr)%Npoints
        i=Vcontact(cr)%Irg(m)
        j=Vcontact(cr)%Jrg(m)
        IF (((IstrT.le.i).and.(i.le.IendT)).and.                        &
     &      ((JstrP.le.j).and.(j.le.JendT))) THEN
          DO k=1,N(ng)
            my_value=Wold*                                              &
     &               (Vcontact(cr)%Lweight(1,m)*                        &
     &                             REFINED(cr)%v(1,k,m,told)+           &
     &                Vcontact(cr)%Lweight(2,m)*                        &
     &                             REFINED(cr)%v(2,k,m,told)+           &
     &                Vcontact(cr)%Lweight(3,m)*                        &
     &                             REFINED(cr)%v(3,k,m,told)+           &
     &                Vcontact(cr)%Lweight(4,m)*                        &
     &                             REFINED(cr)%v(4,k,m,told))+          &
     &               Wnew*                                              &
     &               (Vcontact(cr)%Lweight(1,m)*                        &
     &                             REFINED(cr)%v(1,k,m,tnew)+           &
     &                Vcontact(cr)%Lweight(2,m)*                        &
     &                             REFINED(cr)%v(2,k,m,tnew)+           &
     &                Vcontact(cr)%Lweight(3,m)*                        &
     &                             REFINED(cr)%v(3,k,m,tnew)+           &
     &                Vcontact(cr)%Lweight(4,m)*                        &
     &                             REFINED(cr)%v(4,k,m,tnew))
#  ifdef MASKING
            my_value=my_value*GRID(ng)%vmask(i,j)
#  endif
            OCEAN(ng)%v(i,j,k,1)=my_value
            OCEAN(ng)%v(i,j,k,2)=my_value
          END DO
        END IF
      END DO

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile information.
!-----------------------------------------------------------------------
!
      CALL mp_exchange4d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    OCEAN(ng)%t(:,:,:,1,:),                       &
     &                    OCEAN(ng)%t(:,:,:,2,:),                       &
     &                    OCEAN(ng)%t(:,:,:,3,:))

      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    OCEAN(ng)%u(:,:,:,1),                         &
     &                    OCEAN(ng)%u(:,:,:,2),                         &
     &                    OCEAN(ng)%v(:,:,:,1),                         &
     &                    OCEAN(ng)%v(:,:,:,2))
#  endif
!
  10  FORMAT (/,' PUT_REFINE3D - unbounded contact points temporal: ',  &
     &          ' interpolation:',                                      &
     &        /,2x,  'cr       = ',i2.2,                                &
     &          8x,'dg         = ',i2.2,                                &
     &          8x,'ng         = ',i2.2,                                &
     &        /,2x,  'iic(dg)  = ',i7.7,                                &
     &          3x,'told       = ',i1,                                  &
     &          9x,'tnew       = ',i1,                                  &
     &        /,2x,  'iic(ng)  = ',i7.7,                                &
     &          3x,'Wold       = ',f8.5,                                &
     &          2x,'Wnew       = ',f8.5,                                &
     &        /,2x,  'time(ng) = ',i7.7,                                &
     &          3x,'time(told) = ',i7.7,                                &
     &          3x,'time(tnew) = ',i7.7)

      RETURN
      END SUBROUTINE put_refine3d
# endif

# ifdef SOLVE3D
!
      SUBROUTINE z_weights (ng, model, tile)
!
!=======================================================================
!                                                                      !
!  This routine determines the vertical indices and interpolation      !
!  weights associated with depth,  which are needed to process 3D      !
!  fields in the contact region.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     model      Calling model identifier (integer)                    !
!     tile       Domain partition for composite grid ng (integer)      !
!                                                                      !
!  On Output:    Updated T_NGC type structures in mod_param:           !
!                                                                      !
!     Rcontact   Updated values for Kdg(:,:) and Vweigths (:,:,:)      !
!     Ucontact   Updated values for Kdg(:,:) and Vweigths (:,:,:)      !
!     Vcontact   Updated values for Kdg(:,:) and Vweigths (:,:,:)      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_nesting
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_assemble
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
!
!  Local variable declarations.
!
      integer :: cr, dg, rg, i, j, k, m
      integer :: Idg, Jdg, Kdg, IminD, ImaxD, JminD, JmaxD
      integer :: Irg, Jrg, Krg, IminR, ImaxR, JminR, JmaxR
      integer :: Idgm1, Idgp1, Jdgm1, Jdgp1
      integer :: Npoints
#  ifdef DISTRIBUTE
      integer :: Nkpts, Nwpts, Nzpts

      integer, parameter :: ispv = 0
#  endif
      real(r8), parameter :: spv = 0.0_r8

      real(r8) :: Zbot, Zr, Ztop, dz, r1, r2

      real(r8), allocatable :: Zd(:,:,:)
!
!=======================================================================
!  Compute vertical indices and weights for each contact region.
!=======================================================================
!
      DO cr=1,Ncontact
!
!  Get donor and receiver grid numbers.
!
        dg=Rcontact(cr)%donor_grid
        rg=Rcontact(cr)%receiver_grid
!
!  Process only contact region data for requested nested grid "ng".
!
        IF (rg.eq.ng) THEN
!
!-----------------------------------------------------------------------
!  Process variables in structure Rcontact(cr).
!-----------------------------------------------------------------------
!
!  Get number of contact points to process.
!
          Npoints=Rcontact(cr)%Npoints
!
!  Set starting and ending tile indices for the donor and receiver
!  grids.
!
          IminD=BOUNDS(dg) % IstrT(tile)
          ImaxD=BOUNDS(dg) % IendT(tile)
          JminD=BOUNDS(dg) % JstrT(tile)
          JmaxD=BOUNDS(dg) % JendT(tile)
!
          IminR=BOUNDS(rg) % IstrT(tile)
          ImaxR=BOUNDS(rg) % IendT(tile)
          JminR=BOUNDS(rg) % JstrT(tile)
          JmaxR=BOUNDS(rg) % JendT(tile)

#  ifdef DISTRIBUTE
!
!  If distributed-memory, initialize with special value (zero) to
!  facilitate the global reduction when collecting data between all
!  nodes.
!
          Nkpts=N(rg)*Npoints
          Nwpts=2*Nkpts
          Nzpts=4*Nkpts

          Rcontact(cr)%Kdg(1:N(rg),1:Npoints)=ispv
          Rcontact(cr)%Vweight(1:2,1:N(rg),1:Npoints)=spv
#  endif
!
!  If coincident grids and requested, avoid vertical interpolation.
!
          R_CONTACT : IF (.not.Rcontact(cr)%interpolate.and.            &
     &                    Rcontact(cr)%coincident) THEN
            DO Krg=1,N(rg)
              DO m=1,Npoints
                Irg=Rcontact(cr)%Irg(m)
                Jrg=Rcontact(cr)%Jrg(m)
                IF (((IminR.le.Irg).and.(Irg.le.ImaxR)).and.            &
     &              ((JminR.le.Jrg).and.(Jrg.le.JmaxR))) THEN
                  Rcontact(cr)%Kdg(Krg,m)=Krg
                  Rcontact(cr)%Vweight(1,Krg,m)=1.0_r8
                  Rcontact(cr)%Vweight(2,Krg,m)=0.0_r8
                END IF
              END DO
            END DO
!
!  Otherwise, vertically interpolate because donor and receiver grids
!  are not coincident.
!
          ELSE
!
!  Allocate and initialize local working arrays.
!
            IF (.not.allocated(Zd)) THEN
              allocate ( Zd(4,N(dg),Npoints) )
            END IF
            Zd=spv
!
!  Extract donor grid depths for each cell containing the receiver grid
!  contact point.  Notice that indices i+1 and j+1 are bounded to the
!  maximum possible values in contact points at the edge of the contact
!  region.  In such cases, Lweight(1,m)=1 and Lweight(2:3,m)=0.  This is
!  done to avoid out of range errors. We need to take care of this in
!  the adjoint code.
!
            DO Kdg=1,N(dg)
              DO m=1,Npoints
                Idg  =Rcontact(cr)%Idg(m)
                Idgp1=MIN(Idg+1, BOUNDS(dg)%UBi(-1))
                Jdg  =Rcontact(cr)%Jdg(m)
                Jdgp1=MIN(Jdg+1, BOUNDS(dg)%UBj(-1))
                IF (((IminD.le.Idg).and.(Idg.le.ImaxD)).and.            &
     &              ((JminD.le.Jdg).and.(Jdg.le.JmaxD))) THEN
                  Zd(1,Kdg,m)=GRID(dg)%z_r(Idg  ,Jdg  ,Kdg)
                  Zd(2,Kdg,m)=GRID(dg)%z_r(Idgp1,Jdg  ,Kdg)
                  Zd(3,Kdg,m)=GRID(dg)%z_r(Idgp1,Jdgp1,Kdg)
                  Zd(4,Kdg,m)=GRID(dg)%z_r(Idg  ,Jdgp1,Kdg)
                END IF
              END DO
            END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
            CALL mp_assemble (dg, model, Nzpts, spv, Zd)
            IF (exit_flag.ne.NoError) RETURN
#  endif
!
!  Determine donor grid vertical indices (Kdg) and weights (Vweight)
!  needed for the interpolation of data at the receiver grid contact
!  points.
!
            DO Krg=1,N(rg)
              DO m=1,Npoints
                Irg=Rcontact(cr)%Irg(m)
                Jrg=Rcontact(cr)%Jrg(m)
                IF (((IminR.le.Irg).and.(Irg.le.ImaxR)).and.            &
     &              ((JminR.le.Jrg).and.(Jrg.le.JmaxR))) THEN
                  Ztop=Rcontact(cr)%Lweight(1,m)*Zd(1,N(dg),m)+         &
     &                 Rcontact(cr)%Lweight(2,m)*Zd(2,N(dg),m)+         &
     &                 Rcontact(cr)%Lweight(3,m)*Zd(3,N(dg),m)+         &
     &                 Rcontact(cr)%Lweight(4,m)*Zd(4,N(dg),m)
                  Zbot=Rcontact(cr)%Lweight(1,m)*Zd(1,1    ,m)+         &
     &                 Rcontact(cr)%Lweight(2,m)*Zd(2,1    ,m)+         &
     &                 Rcontact(cr)%Lweight(3,m)*Zd(3,1    ,m)+         &
     &                 Rcontact(cr)%Lweight(4,m)*Zd(4,1    ,m)
                  Zr=GRID(rg)%z_r(Irg,Jrg,Krg)
                  IF (Zr.ge.Ztop) THEN           ! If shallower, use top
                    Rcontact(cr)%Kdg(Krg,m)=N(dg)! donor grid cell value
                    Rcontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Rcontact(cr)%Vweight(2,Krg,m)=1.0_r8
                  ELSE IF (Zbot.ge.Zr) THEN      ! If deeper, use bottom
                    Rcontact(cr)%Kdg(Krg,m)=1    ! donor grid cell value
                    Rcontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Rcontact(cr)%Vweight(2,Krg,m)=1.0_r8
                  ELSE                           ! bounded, interpolate
                    DO Kdg=N(dg),2,-1
                      Ztop=Rcontact(cr)%Lweight(1,m)*Zd(1,Kdg  ,m)+     &
     &                     Rcontact(cr)%Lweight(2,m)*Zd(2,Kdg  ,m)+     &
     &                     Rcontact(cr)%Lweight(3,m)*Zd(3,Kdg  ,m)+     &
     &                     Rcontact(cr)%Lweight(4,m)*Zd(4,Kdg  ,m)
                      Zbot=Rcontact(cr)%Lweight(1,m)*Zd(1,Kdg-1,m)+     &
     &                     Rcontact(cr)%Lweight(2,m)*Zd(2,Kdg-1,m)+     &
     &                     Rcontact(cr)%Lweight(3,m)*Zd(3,Kdg-1,m)+     &
     &                     Rcontact(cr)%Lweight(4,m)*Zd(4,Kdg-1,m)
                      IF ((Ztop.gt.Zr).and.(Zr.ge.Zbot)) THEN
                        dz=Ztop-Zbot
                        r2=(Zr-Zbot)/dz
                        r1=1.0_r8-r2
                        Rcontact(cr)%Kdg(Krg,m)=Kdg
                        Rcontact(cr)%Vweight(1,Krg,m)=r1
                        Rcontact(cr)%Vweight(2,Krg,m)=r2
                      END IF
                    END DO
                  END IF
                END IF
              END DO
            END DO
          END IF R_CONTACT

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
          CALL mp_assemble (rg, model, Nkpts, ispv, Rcontact(cr)%Kdg)
          IF (exit_flag.ne.NoError) RETURN

          CALL mp_assemble (rg, model, Nwpts, spv, Rcontact(cr)%Vweight)
          IF (exit_flag.ne.NoError) RETURN
#  endif
!
!  Deallocate local work arrays.
!
          IF (allocated(Zd)) THEN
            deallocate (Zd)
          END IF
!
!-----------------------------------------------------------------------
!  Process variables in structure Ucontact(cr).
!-----------------------------------------------------------------------
!
!  Get number of contact points to process.
!
          Npoints=Ucontact(cr)%Npoints
!
!  Set starting and ending tile indices for the donor and receiver
!  grids.
!
          IminD=BOUNDS(dg) % IstrP(tile)
          ImaxD=BOUNDS(dg) % IendT(tile)
          JminD=BOUNDS(dg) % JstrT(tile)
          JmaxD=BOUNDS(dg) % JendT(tile)
!
          IminR=BOUNDS(rg) % IstrP(tile)
          ImaxR=BOUNDS(rg) % IendT(tile)
          JminR=BOUNDS(rg) % JstrT(tile)
          JmaxR=BOUNDS(rg) % JendT(tile)

#  ifdef DISTRIBUTE
!
!  If distributed-memory, initialize with special value (zero) to
!  facilitate the global reduction when collecting data between all
!  nodes.
!
          Nkpts=N(rg)*Npoints
          Nwpts=2*Nkpts
          Nzpts=4*Nkpts

          Ucontact(cr)%Kdg(1:N(rg),1:Npoints)=ispv
          Ucontact(cr)%Vweight(1:2,1:N(rg),1:Npoints)=spv
#  endif
!
!  If coincident grids and requested, avoid vertical interpolation.
!
          U_CONTACT : IF (.not.Ucontact(cr)%interpolate.and.            &
     &                    Ucontact(cr)%coincident) THEN
            DO Krg=1,N(rg)
              DO m=1,Npoints
                Irg=Ucontact(cr)%Irg(m)
                Jrg=Ucontact(cr)%Jrg(m)
                IF (((IminR.le.Irg).and.(Irg.le.ImaxR)).and.            &
     &              ((JminR.le.Jrg).and.(Jrg.le.JmaxR))) THEN
                  Ucontact(cr)%Kdg(Krg,m)=Krg
                  Ucontact(cr)%Vweight(1,Krg,m)=1.0_r8
                  Ucontact(cr)%Vweight(2,Krg,m)=0.0_r8
                END IF
              END DO
            END DO
!
!  Otherwise, vertically interpolate because donor and receiver grids
!  are not coincident.
!
          ELSE
!
!  Allocate and initialize local working arrays.
!
            IF (.not.allocated(Zd)) THEN
              allocate (Zd(4,N(dg),Npoints))
            END IF
            Zd=spv
!
!  Extract donor grid depths for each cell containing the receiver grid
!  contact point.  Notice that indices i-1, i+1 and j-1, j+1 are bounded
!  the minimum/maximum possible values in contact points at the edge of
!  the contact region.  In such cases, the interpolation weights
!  Lweight(1,m)=1 and Lweight(2:3,m)=0.  This is done to avoid out of
!  range errors. We need to take care of this in the adjoint code.
!
            DO Kdg=1,N(dg)
              DO m=1,Npoints
                Idg  =Ucontact(cr)%Idg(m)
                Idgm1=MAX(Idg-1, BOUNDS(dg)%LBi(-1))
                Idgp1=MIN(Idg+1, BOUNDS(dg)%UBi(-1))
                Jdg  =Ucontact(cr)%Jdg(m)
                Jdgp1=MIN(Jdg+1, BOUNDS(dg)%UBj(-1))
                IF (((IminD.le.Idg).and.(Idg.le.ImaxD)).and.            &
     &              ((JminD.le.Jdg).and.(Jdg.le.JmaxD))) THEN
                  Zd(1,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idgm1,Jdg  ,Kdg)+    &
     &                                GRID(dg)%z_r(Idg  ,Jdg  ,Kdg))
                  Zd(2,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idg  ,Jdg  ,Kdg)+    &
     &                                GRID(dg)%z_r(Idgp1,Jdg  ,Kdg))
                  Zd(3,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idg  ,Jdgp1,Kdg)+    &
     &                                GRID(dg)%z_r(Idgp1,Jdgp1,Kdg))
                  Zd(4,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idgm1,Jdgp1,Kdg)+    &
     &                                GRID(dg)%z_r(Idg  ,Jdgp1,Kdg))
                END IF
              END DO
            END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
            CALL mp_assemble (dg, model, Nzpts, spv, Zd)
            IF (exit_flag.ne.NoError) RETURN
#  endif
!
!  Determine donor grid vertical indices (Kdg) and weights (Vweight)
!  needed for the interpolation of data at the receiver grid contact
!  points.
!
            DO Krg=1,N(rg)
              DO m=1,Npoints
                Irg=Ucontact(cr)%Irg(m)
                Jrg=Ucontact(cr)%Jrg(m)
                IF (((IminR.le.Irg).and.(Irg.le.ImaxR)).and.            &
     &              ((JminR.le.Jrg).and.(Jrg.le.JmaxR))) THEN
                  Ztop=Ucontact(cr)%Lweight(1,m)*Zd(1,N(dg),m)+         &
     &                 Ucontact(cr)%Lweight(2,m)*Zd(2,N(dg),m)+         &
     &                 Ucontact(cr)%Lweight(3,m)*Zd(3,N(dg),m)+         &
     &                 Ucontact(cr)%Lweight(4,m)*Zd(4,N(dg),m)
                  Zbot=Ucontact(cr)%Lweight(1,m)*Zd(1,1    ,m)+         &
     &                 Ucontact(cr)%Lweight(2,m)*Zd(2,1    ,m)+         &
     &                 Ucontact(cr)%Lweight(3,m)*Zd(3,1    ,m)+         &
     &                 Ucontact(cr)%Lweight(4,m)*Zd(4,1    ,m)
                  Zr=0.5_r8*(GRID(rg)%z_r(Irg  ,Jrg,Krg)+               &
     &                       GRID(rg)%z_r(Irg-1,Jrg,Krg))
                  IF (Zr.ge.Ztop) THEN           ! If shallower, use top
                    Ucontact(cr)%Kdg(Krg,m)=N(dg)! donor grid cell value
                    Ucontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Ucontact(cr)%Vweight(2,Krg,m)=1.0_r8
                  ELSE IF (Zbot.ge.Zr) THEN      ! If deeper, use bottom
                    Ucontact(cr)%Kdg(Krg,m)=1    ! donor grid cell value
                    Ucontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Ucontact(cr)%Vweight(2,Krg,m)=1.0_r8
                  ELSE                           ! bounded, interpolate
                    DO Kdg=N(dg),2,-1
                      Ztop=Ucontact(cr)%Lweight(1,m)*Zd(1,Kdg  ,m)+     &
     &                     Ucontact(cr)%Lweight(2,m)*Zd(2,Kdg  ,m)+     &
     &                     Ucontact(cr)%Lweight(3,m)*Zd(3,Kdg  ,m)+     &
     &                     Ucontact(cr)%Lweight(4,m)*Zd(4,Kdg  ,m)
                      Zbot=Ucontact(cr)%Lweight(1,m)*Zd(1,Kdg-1,m)+     &
     &                     Ucontact(cr)%Lweight(2,m)*Zd(2,Kdg-1,m)+     &
     &                     Ucontact(cr)%Lweight(3,m)*Zd(3,Kdg-1,m)+     &
     &                     Ucontact(cr)%Lweight(4,m)*Zd(4,Kdg-1,m)
                      IF ((Ztop.gt.Zr).and.(Zr.ge.Zbot)) THEN
                        dz=Ztop-Zbot
                        r2=(Zr-Zbot)/dz
                        r1=1.0_r8-r2
                        Ucontact(cr)%Kdg(Krg,m)=Kdg
                        Ucontact(cr)%Vweight(1,Krg,m)=r1
                        Ucontact(cr)%Vweight(2,Krg,m)=r2
                      END IF
                    END DO
                  END IF
                END IF
              END DO
            END DO
          END IF U_CONTACT

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
          CALL mp_assemble (rg, model, Nkpts, ispv, Ucontact(cr)%Kdg)
          IF (exit_flag.ne.NoError) RETURN

          CALL mp_assemble (rg, model, Nwpts, spv, Ucontact(cr)%Vweight)
          IF (exit_flag.ne.NoError) RETURN
#  endif
!
!  Deallocate local work arrays.
!
          IF (allocated(Zd)) THEN
            deallocate (Zd)
          END IF
!
!-----------------------------------------------------------------------
!  Process variables in structure Vcontact(cr).
!-----------------------------------------------------------------------
!
!  Get number of contact points to process.
!
          Npoints=Vcontact(cr)%Npoints
!
!  Set starting and ending tile indices for the donor and receiver
!  grids.
!
          IminD=BOUNDS(dg) % IstrT(tile)
          ImaxD=BOUNDS(dg) % IendT(tile)
          JminD=BOUNDS(dg) % JstrP(tile)
          JmaxD=BOUNDS(dg) % JendT(tile)
!
          IminR=BOUNDS(rg) % IstrT(tile)
          ImaxR=BOUNDS(rg) % IendT(tile)
          JminR=BOUNDS(rg) % JstrP(tile)
          JmaxR=BOUNDS(rg) % JendT(tile)

#  ifdef DISTRIBUTE
!
!  If distributed-memory, initialize with special value (zero) to
!  facilitate the global reduction when collecting data between all
!  nodes.
!
          Nkpts=N(rg)*Npoints
          Nwpts=2*Nkpts
          Nzpts=4*Nkpts

          Vcontact(cr)%Kdg(1:N(rg),1:Npoints)=ispv
          Vcontact(cr)%Vweight(1:2,1:N(rg),1:Npoints)=spv
#  endif
!
!  If coincident grids and requested, avoid vertical interpolation.
!
          V_CONTACT : IF (.not.Vcontact(cr)%interpolate.and.            &
     &                    Vcontact(cr)%coincident) THEN
            DO Krg=1,N(rg)
              DO m=1,Npoints
                Irg=Vcontact(cr)%Irg(m)
                Jrg=Vcontact(cr)%Jrg(m)
                IF (((IminR.le.Irg).and.(Irg.le.ImaxR)).and.            &
     &              ((JminR.le.Jrg).and.(Jrg.le.JmaxR))) THEN
                  Vcontact(cr)%Kdg(Krg,m)=Krg
                  Vcontact(cr)%Vweight(1,Krg,m)=1.0_r8
                  Vcontact(cr)%Vweight(2,Krg,m)=0.0_r8
                END IF
              END DO
            END DO
!
!  Otherwise, vertically interpolate because donor and receiver grids
!  are not coincident.
!
          ELSE
!
!  Allocate and initialize local working arrays.
!
            IF (.not.allocated(Zd)) THEN
              allocate (Zd(4,N(dg),Npoints))
            END IF
            Zd=spv
!
!  Extract donor grid depths for each cell containing the receiver grid
!  contact point.
!
            DO Kdg=1,N(dg)
              DO m=1,Npoints
                Idg=Vcontact(cr)%Idg(m)
                Idgp1=MIN(Idg+1, BOUNDS(dg)%UBi(-1))
                Jdg=Vcontact(cr)%Jdg(m)
                Jdgm1=MAX(Jdg-1, BOUNDS(dg)%LBj(-1))
                Jdgp1=MIN(Jdg+1, BOUNDS(dg)%UBj(-1))
                IF (((IminD.le.Idg).and.(Idg.le.ImaxD)).and.            &
     &              ((JminD.le.Jdg).and.(Jdg.le.JmaxD))) THEN
                  Zd(1,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idg  ,Jdgm1,Kdg)+    &
     &                                GRID(dg)%z_r(Idg  ,Jdg  ,Kdg))
                  Zd(2,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idgp1,Jdgm1,Kdg)+    &
     &                                GRID(dg)%z_r(Idgp1,Jdg  ,Kdg))
                  Zd(3,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idgp1,Jdg  ,Kdg)+    &
     &                                GRID(dg)%z_r(Idgp1,Jdgp1,Kdg))
                  Zd(4,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idg  ,Jdg  ,Kdg)+    &
     &                                GRID(dg)%z_r(Idg  ,Jdgp1,Kdg))
                END IF
              END DO
            END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
            CALL mp_assemble (dg, model, Nzpts, spv, Zd)
            IF (exit_flag.ne.NoError) RETURN
#  endif
!
!  Determine donor grid vertical indices (Kdg) and weights (Vweight)
!  needed for the interpolation of data at the receiver grid contact
!  points.
!
            DO Krg=1,N(rg)
              DO m=1,Npoints
                Irg=Vcontact(cr)%Irg(m)
                Jrg=Vcontact(cr)%Jrg(m)
                IF (((IminR.le.Irg).and.(Irg.le.ImaxR)).and.            &
     &              ((JminR.le.Jrg).and.(Jrg.le.JmaxR))) THEN
                  Ztop=Vcontact(cr)%Lweight(1,m)*Zd(1,N(dg),m)+         &
     &                 Vcontact(cr)%Lweight(2,m)*Zd(2,N(dg),m)+         &
     &                 Vcontact(cr)%Lweight(3,m)*Zd(3,N(dg),m)+         &
     &                 Vcontact(cr)%Lweight(4,m)*Zd(4,N(dg),m)
                  Zbot=Vcontact(cr)%Lweight(1,m)*Zd(1,1    ,m)+         &
     &                 Vcontact(cr)%Lweight(2,m)*Zd(2,1    ,m)+         &
     &                 Vcontact(cr)%Lweight(3,m)*Zd(3,1    ,m)+         &
     &                 Vcontact(cr)%Lweight(4,m)*Zd(4,1    ,m)
                  Zr=0.5_r8*(GRID(rg)%z_r(Irg,Jrg  ,Krg)+               &
     &                       GRID(rg)%z_r(Irg,Jrg-1,Krg))
                  IF (Zr.ge.Ztop) THEN           ! If shallower, use top
                    Vcontact(cr)%Kdg(Krg,m)=N(dg)! donor grid cell value
                    Vcontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Vcontact(cr)%Vweight(2,Krg,m)=1.0_r8
                  ELSE IF (Zbot.ge.Zr) THEN      ! If deeper, use bottom
                    Vcontact(cr)%Kdg(Krg,m)=1    ! donor grid cell value
                    Vcontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Vcontact(cr)%Vweight(2,Krg,m)=1.0_r8
                  ELSE                           ! bounded, interpolate
                    DO Kdg=N(dg),2,-1
                      Ztop=Vcontact(cr)%Lweight(1,m)*Zd(1,Kdg  ,m)+     &
     &                     Vcontact(cr)%Lweight(2,m)*Zd(2,Kdg  ,m)+     &
     &                     Vcontact(cr)%Lweight(3,m)*Zd(3,Kdg  ,m)+     &
     &                     Vcontact(cr)%Lweight(4,m)*Zd(4,Kdg  ,m)
                      Zbot=Vcontact(cr)%Lweight(1,m)*Zd(1,Kdg-1,m)+     &
     &                     Vcontact(cr)%Lweight(2,m)*Zd(2,Kdg-1,m)+     &
     &                     Vcontact(cr)%Lweight(3,m)*Zd(3,Kdg-1,m)+     &
     &                     Vcontact(cr)%Lweight(4,m)*Zd(4,Kdg-1,m)
                      IF ((Ztop.gt.Zr).and.(Zr.ge.Zbot)) THEN
                        dz=Ztop-Zbot
                        r2=(Zr-Zbot)/dz
                        r1=1.0_r8-r2
                        Vcontact(cr)%Kdg(Krg,m)=Kdg
                        Vcontact(cr)%Vweight(1,Krg,m)=r1
                        Vcontact(cr)%Vweight(2,Krg,m)=r2
                      END IF
                    END DO
                  END IF
                END IF
              END DO
            END DO
          END IF V_CONTACT

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
          CALL mp_assemble (rg, model, Nkpts, ispv, Vcontact(cr)%Kdg)
          IF (exit_flag.ne.NoError) RETURN

          CALL mp_assemble (rg, model, Nwpts, spv, Vcontact(cr)%Vweight)
          IF (exit_flag.ne.NoError) RETURN
#  endif
!
!  Deallocate local work arrays.
!
          IF (allocated(Zd)) THEN
            deallocate (Zd)
          END IF

        END IF
      END DO

      RETURN

      END SUBROUTINE z_weights
# endif
#endif
      END MODULE nesting_mod
