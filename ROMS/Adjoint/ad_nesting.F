#include "cppdefs.h"
      MODULE ad_nesting_mod

#if defined NESTING && defined ADJOINT
!
!svn $Id: ad_nesting.F 1058 2021-03-25 19:04:23Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
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
!  ad_nesting         Public interface to time-stepping kernel         !
!                                                                      !
!  ad_get_composite   Composite  grid, extracts contact points donor   !
!                       data                                           !
!  ad_get_refine      Refinement grid, extracts contact points donor   !
!                       data                                           !
!  ad_put_composite   Composite  grid, fills contact points (spatial   !
!                       interpolation)                                 !
!  ad_put_refine      Refinement grid, fills contact points (spatial   !
!                       and temporal interpolation)                    !
!                                                                      !
!  ad_bry_fluxes      Extracts horizontat advective fluxes the contact !
!                       boundary of donor and receiver grids           !
# ifdef NESTING_DEBUG
!  ad_check_massflux  If refinement, checks mass fluxes between coarse !
!                       and fine grids for volume conservation. It is  !
!                       use only for debugging and diagnostics.        !
# endif
!  ad_correct_tracer  Corrects coarse grid tracer at the refinement    !
!                       grid boundary with the refined accumulated     !
!                       fluxes                                         !
!  ad_fine2coarse     Replace coarse grid state variables with the     !
!                       averaged fine grid values (two-way nesting)    !
!                                                                      !
!  ad_get_contact2d   Gets 2D field donor grid cell holding contact    !
!                       point                                          !
!  ad_get_contact3d   Gets 3D field donor grid cell holding contact    !
!                       point                                          !
!  ad_get_persisted2d Gets 2D field persisted values on contact points !
!  ad_put_contact2d   Sets 2D field contact points, spatial            !
!                       interpolation                                  !
!  ad_put_contact3d   Sets 3D field contact points, spatial            !
!                       interpolation                                  !
!                                                                      !
!  ad_put_refine2d    Interpolates (space-time) 2D state variables     !
!  ad_put_refine3d    Interpolates (space-time) 3D state variables     !
!                                                                      !
!  ad_z_weights       Sets donor grid vertical indices (cell holding   !
!                       contact point) and vertical interpolation      !
!                       weights                                        !
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
      PUBLIC  :: ad_nesting
      PUBLIC  :: ad_bry_fluxes
# ifdef NESTING_DEBUG
      PRIVATE :: ad_check_massflux
# endif
# ifdef SOLVE3D
      PRIVATE :: ad_correct_tracer
      PRIVATE :: ad_correct_tracer_tile
# endif
      PRIVATE :: ad_fine2coarse
      PRIVATE :: ad_fine2coarse2d
# ifdef SOLVE3D
      PRIVATE :: ad_fine2coarse3d
# endif
      PRIVATE :: ad_get_contact2d
# ifdef SOLVE3D
      PRIVATE :: ad_get_contact3d
# endif
      PRIVATE :: ad_get_composite
      PRIVATE :: ad_get_persisted2d
      PRIVATE :: ad_get_refine
      PRIVATE :: ad_put_composite
      PRIVATE :: ad_put_refine
      PRIVATE :: ad_put_refine2d
# ifdef SOLVE3D
      PRIVATE :: ad_put_refine3d
      PRIVATE :: ad_z_weights
# endif
!
      CONTAINS
!
      SUBROUTINE ad_nesting (ng, model, isection)
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
!
# ifdef SOLVE3D
      USE set_depth_mod,    ONLY : set_depth
# endif
      USE ad_set_depth_mod, ONLY : ad_set_depth
      USE nesting_mod,      ONLY : get_metrics
      USE nesting_mod,      ONLY : mask_hweights
      USE nesting_mod,      ONLY : z_weights
      USE strings_mod,      ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, isection
!
!  Local variable declarations.
!
      logical :: LputFsur
      integer :: subs, tile, thread
      integer :: ngc
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Process vertical indices and interpolation weigths associated with
!  depth. Currently, vertical weights are used in composite grids
!  because their grids are not coincident. They are not needed in
!  refinement grids because the donor and receiver grids have the same
!  number of vertical levels and matching bathymetry.  However, in
!  the future, it is possible to have configurations that require
!  vertical weights in refinement. The switch "get_Vweights" controls
!  if such weights are computed or not.  If false, it will accelerate
!  computations because of less distributed-memory communications.
!-----------------------------------------------------------------------
!
      IF ((isection.eq.nzwgt).and.get_Vweights) THEN
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
!>          CALL tl_get_refine (ng, model, tile)
!>
            CALL ad_get_refine (ng, model, tile)
          END DO
!$OMP BARRIER
!
!  Fill refinement grid contact points variables by interpolating
!  (space, time) from extracted donor grid data.  The free-surface
!  needs to be processed first and in a separate parallel region
!  because of shared-memory applications.
!
        ELSE IF (isection.eq.nputD) THEN
!$OMP BARRIER
!
          LputFsur=.FALSE.
          DO tile=first_tile(ng),last_tile(ng),+1
!>          CALL tl_put_refine (ng, model, tile, LputFsur)
!>
            CALL ad_put_refine (ng, model, tile, LputFsur)
          END DO
!$OMP BARRIER
          LputFsur=.TRUE.
          DO tile=first_tile(ng),last_tile(ng),+1
!>          CALL tl_put_refine (ng, model, tile, LputFsur)
!>
            CALL ad_put_refine (ng, model, tile, LputFsur)
          END DO

# ifdef NESTING_DEBUG
!
!  If refinement, check mass flux conservation between coarser and
!  finer grids. DIAGNOSTIC ONLY.
!
        ELSE IF (isection.eq.nmflx) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
!>          CALL tl_check_massflux (ng, model, tile)
!>
            CALL ad_check_massflux (ng, model, tile)
          END DO
# endif

# ifndef ONE_WAY
!
!  Fine to coarse coupling: two-way nesting.
!
        ELSE IF (isection.eq.n2way) THEN

          ngc=CoarserDonor(ng)               ! coarse grid number
#  ifdef SOLVE3D
!$OMP BARRIER
!
!  Replace coarse grid 3D state variables with the averaged fine grid
!  values (two-way coupling).
!
          DO tile=last_tile(ngc),first_tile(ngc),-1
!>          CALL tl_fine2coarse (ng, model, r3dvar, tile)
!>
            CALL ad_fine2coarse (ng, model, r3dvar, tile)
          END DO
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Update coarse grid depth variables. We have a new coarse grid
!  adjusted free-surface, Zt_avg1.
!
          DO tile=first_tile(ngc),last_tile(ngc),+1
!>          CALL tl_set_depth (ngc, tile, model)
!>
            CALL ad_set_depth (ngc, tile, model)
          END DO
#  endif
!
!  Replace coarse grid 2D state variables with the averaged fine grid
!  values (two-way coupling).
!
          DO tile=last_tile(ngc),first_tile(ngc),-1
!>          CALL fine2coarse (ng, model, r2dvar, tile)
!>
            CALL ad_fine2coarse (ng, model, r2dvar, tile)
          END DO
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined SOLVE3D && !defined NO_CORRECT_TRACER
!
!  Correct coarse grid tracer values at the refinement grid, ng,
!  boundary with the refined accumulated fluxes (Hz*u*T/n, Hz*v*T/m).
!
          DO tile=first_tile(ngc),last_tile(ngc),+1
!>          CALL tl_correct_tracer (ngc, ng, model, tile)
!>
            CALL ad_correct_tracer (ngc, ng, model, tile)
          END DO
!$OMP BARRIER
#  endif
# else
!
!  Fine to coarse coupling (two-way nesting) is not activated!
!
        ELSE IF (isection.eq.n2way) THEN
# endif
        END IF

      ELSE
!
!-----------------------------------------------------------------------
!  Otherwise, process contact points in composite grid.
!-----------------------------------------------------------------------
!

!$OMP BARRIER
!
!  Fill composite grid contact points variables by interpolating from
!  extracted donor grid data.
!
        DO tile=last_tile(ng),first_tile(ng),-1
!>        CALL tl_put_composite (ng, model, isection, tile)
!>
          CALL ad_put_composite (ng, model, isection, tile)
        END DO
!$OMP BARRIER

!
!  Get composite grid contact points data from donor grid. It extracts
!  the donor grid cell data necessary to interpolate state variables
!  at each contact point.
!
        DO tile=first_tile(ng),last_tile(ng),+1
!>        CALL tl_get_composite (ng, model, isection, tile)
!>
          CALL ad_get_composite (ng, model, isection, tile)
        END DO

      END IF

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Process vertical indices and interpolation weigths associated with
!  depth.
!-----------------------------------------------------------------------
!
      IF (isection.eq.nzwgt) THEN
        DO tile=last_tile(ng),first_tile(ng),-1
!>        CALL tl_z_weights (ng, model, tile)
!>
          CALL ad_z_weights (ng, model, tile)
        END DO
!$OMP BARRIER
        RETURN
      END IF
# endif

      RETURN
      END SUBROUTINE ad_nesting
!
      SUBROUTINE ad_get_composite (ng, model, isection, tile)
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
      USE nesting_mod, ONLY : get_contact2d
# ifdef SOLVE3D
      USE nesting_mod, ONLY : get_contact3d
# endif

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
!>          CALL get_contact2d (dg, model, tile,                        &
!>   &                          u2dvar, Vname(1,idUbms),                &
!>   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          FORCES(dg) % tl_bustr,                  &
!>   &                          COMPOSITE(cr) % tl_bustr)
!>
            CALL ad_get_contact2d (dg, model, tile,                     &
     &                             u2dvar, Vname(1,idUbms),             &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             FORCES(dg) % ad_bustr,               &
     &                             COMPOSITE(cr) % ad_bustr)

!>          CALL get_contact2d (dg, model, tile,                        &
!>   &                          v2dvar, Vname(1,idVbms),                &
!>   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          FORCES(dg) % tl_bvstr,                  &
!>   &                          COMPOSITE(cr) % tl_bvstr)
!>
            CALL ad_get_contact2d (dg, model, tile,                     &
     &                             v2dvar, Vname(1,idVbms),             &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             FORCES(dg) % ad_bvstr,               &
     &                             COMPOSITE(cr) % ad_bvstr)
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
!>            CALL get_contact2d (dg, model, tile,                      &
!>   &                            r2dvar, Vname(1,idFsur),              &
!>   &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            OCEAN(dg) % tl_zeta(:,:,Tindex),      &
!>   &                            COMPOSITE(cr) % tl_zeta(:,:,rec))
!>
              CALL ad_get_contact2d (dg, model, tile,                   &
     &                               r2dvar, Vname(1,idFsur),           &
     &                               cr, Rcontact(cr)%Npoints, Rcontact,&
     &                               LBi, UBi, LBj, UBj,                &
     &                               OCEAN(dg) % ad_zeta(:,:,Tindex),   &
     &                               COMPOSITE(cr) % ad_zeta(:,:,rec))
            END DO
          END IF
!
!  Process free-surface equation rigth-hand-side (rzeta) term.
!
          IF (isection.eq.n2dPS) THEN
            Tindex=1
!>          CALL get_contact2d (dg, model, tile,                        &
!>   &                          r2dvar, Vname(1,idRzet),                &
!>   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          OCEAN(dg) % tl_rzeta(:,:,Tindex),       &
!>   &                          COMPOSITE(cr) % tl_rzeta)
!>
            CALL ad_get_contact2d (dg, model, tile,                     &
     &                             r2dvar, Vname(1,idRzet),             &
     &                             cr, Rcontact(cr)%Npoints, Rcontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             OCEAN(dg) % ad_rzeta(:,:,Tindex),    &
     &                             COMPOSITE(cr) % ad_rzeta)
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
!>            CALL get_contact2d (dg, model, tile,                      &
!>   &                            u2dvar, Vname(1,idUbar),              &
!>   &                            cr, Ucontact(cr)%Npoints, Ucontact,   &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            OCEAN(dg) % tl_ubar(:,:,Tindex),      &
!>   &                            COMPOSITE(cr) % tl_ubar(:,:,rec))
!>
              CALL ad_get_contact2d (dg, model, tile,                   &
     &                               u2dvar, Vname(1,idUbar),           &
     &                               cr, Ucontact(cr)%Npoints, Ucontact,&
     &                               LBi, UBi, LBj, UBj,                &
     &                               OCEAN(dg) % ad_ubar(:,:,Tindex),   &
     &                               COMPOSITE(cr) % ad_ubar(:,:,rec))

!>            CALL get_contact2d (dg, model, tile,                      &
!>   &                            v2dvar, Vname(1,idVbar),              &
!>   &                            cr, Vcontact(cr)%Npoints, Vcontact,   &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            OCEAN(dg) % tl_vbar(:,:,Tindex),      &
!>   &                            COMPOSITE(cr) % tl_vbar(:,:,rec))
!>
              CALL ad_get_contact2d (dg, model, tile,                   &
     &                               v2dvar, Vname(1,idVbar),           &
     &                               cr, Vcontact(cr)%Npoints, Vcontact,&
     &                               LBi, UBi, LBj, UBj,                &
     &                               OCEAN(dg) % ad_vbar(:,:,Tindex),   &
     &                               COMPOSITE(cr) % ad_vbar(:,:,rec))
            END DO
          END IF

# ifdef SOLVE3D
!
!  Process time averaged free-surface (Zt_avg1) and 2D momentum fluxes
!  (DU_avg1, DV_avg1).
!
          IF (isection.eq.n2dfx) THEN
!>          CALL get_contact2d (dg, model, tile,                        &
!>   &                          r2dvar, 'Zt_avg1',                      &
!>   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          COUPLING(dg) % tl_Zt_avg1,              &
!>   &                          COMPOSITE(cr) % tl_Zt_avg1)
!>
            CALL ad_get_contact2d (dg, model, tile,                     &
     &                             r2dvar, 'Zt_avg1',                   &
     &                             cr, Rcontact(cr)%Npoints, Rcontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             COUPLING(dg) % ad_Zt_avg1,           &
     &                             COMPOSITE(cr) % ad_Zt_avg1)
!
!  Do we need to get DU_avg1 and DV_avg1 here? YES.
!
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

!>          CALL get_contact2d (dg, model, tile,                        &
!>   &                          u2dvar, 'DU_avg1',                      &
!>   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          COUPLING(dg) % tl_DU_avg1,              &
!>   &                          COMPOSITE(cr) % tl_DU_avg1)
!>
            CALL ad_get_contact2d (dg, model, tile,                     &
     &                             u2dvar, 'DU_avg1',                   &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             COUPLING(dg) % ad_DU_avg1,           &
     &                             COMPOSITE(cr) % ad_DU_avg1)

!>          CALL get_contact2d (dg, model, tile,                        &
!>   &                          v2dvar, 'DV_avg1',                      &
!>   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          COUPLING(dg) % tl_DV_avg1,              &
!>   &                          COMPOSITE(cr) % tl_DV_avg1)
!>
            CALL ad_get_contact2d (dg, model, tile,                     &
     &                             v2dvar, 'DV_avg1',                   &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             COUPLING(dg) % ad_DV_avg1,           &
     &                             COMPOSITE(cr) % ad_DV_avg1)
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
!>            CALL get_contact3d (dg, model, tile,                      &
!>   &                            r3dvar, Vname(1,idTvar(itrc)),        &
!>   &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
!>   &                            LBi, UBi, LBj, UBj, 1, N(dg),         &
!>   &                            OCEAN(dg) % tl_t(:,:,:,Tindex,itrc),  &
!>   &                            COMPOSITE(cr) % tl_t(:,:,:,itrc))
!>
              CALL ad_get_contact3d (dg, model, tile,                   &
     &                               r3dvar, Vname(1,idTvar(itrc)),     &
     &                               cr, Rcontact(cr)%Npoints, Rcontact,&
     &                               LBi, UBi, LBj, UBj, 1, N(dg),      &
     &                               OCEAN(dg) % ad_t(:,:,:,Tindex,     &
     &                                                itrc),            &
     &                               COMPOSITE(cr) % ad_t(:,:,:,itrc))
            END DO
          END IF
#  endif
!
!  Process 3D momentum (u, v) at the appropriate time-index.
!
          IF ((isection.eq.n3dIC).or.                                   &
     &        (isection.eq.n3duv)) THEN
            Tindex=nnew(dg)
!>          CALL get_contact3d (dg, model, tile,                        &
!>   &                          u3dvar, Vname(1,idUvel),                &
!>   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!>   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!>   &                          OCEAN(dg) % tl_u(:,:,:,Tindex),         &
!>   &                          COMPOSITE(cr) % tl_u)
!>
            CALL ad_get_contact3d (dg, model, tile,                     &
     &                             u3dvar, Vname(1,idUvel),             &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(dg),        &
     &                             OCEAN(dg) % ad_u(:,:,:,Tindex),      &
     &                             COMPOSITE(cr) % ad_u)
!>          CALL get_contact3d (dg, model, tile,                        &
!>   &                          v3dvar, Vname(1,idVvel),                &
!>   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!>   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!>   &                          OCEAN(dg) % tl_v(:,:,:,Tindex),         &
!>   &                          COMPOSITE(cr) % tl_v)
!>
            CALL ad_get_contact3d (dg, model, tile,                     &
     &                             v3dvar, Vname(1,idVvel),             &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(dg),        &
     &                             OCEAN(dg) % ad_v(:,:,:,Tindex),      &
     &                             COMPOSITE(cr) % ad_v)
          END IF
!
!  Process 3D momentum fluxes (Huon, Hvom).
!
          IF (isection.eq.n3duv) THEN
!>          CALL get_contact3d (dg, model, tile,                        &
!>   &                          u3dvar, 'Huon',                         &
!>   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!>   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!>   &                          GRID(dg) % tl_Huon,                     &
!>   &                          COMPOSITE(cr) % tl_Huon)
!>
            CALL ad_get_contact3d (dg, model, tile,                     &
     &                             u3dvar, 'Huon',                      &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(dg),        &
     &                             GRID(dg) % ad_Huon,                  &
     &                             COMPOSITE(cr) % ad_Huon)
!>          CALL get_contact3d (dg, model, tile,                        &
!>   &                          v3dvar, 'Hvom',                         &
!>   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!>   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!>   &                          GRID(dg) % tl_Hvom,                     &
!>   &                          COMPOSITE(cr) % tl_Hvom)
!>
            CALL ad_get_contact3d (dg, model, tile,                     &
     &                             v3dvar, 'Hvom',                      &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(dg),        &
     &                             GRID(dg) % ad_Hvom,                  &
     &                             COMPOSITE(cr) % ad_Hvom)
          END IF
# endif

        END IF
      END DO

      RETURN
      END SUBROUTINE ad_get_composite
!
      SUBROUTINE ad_get_refine (ng, model, tile)
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

      USE nesting_mod, ONLY : get_persisted2d
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
      integer :: Tindex2d, cr, dg, ir, rg, told, tnew
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
!         tnew=3-RollingIndex(cr)
          tnew=RollingIndex(cr)
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
!  If debugging, write information into Fortran unit 102 to check the
!  logic of processing donor grid data.
!
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              IF (first) THEN
                first=.FALSE.
                WRITE (102,10)
              END IF
              WRITE (102,20) ng, cr, dg, rg, iic(dg), iic(rg),          &
     &                       3-tnew, tnew, Tindex2d, Tindex3d,          &
     &                       INT(time(rg)),                             &
     &                       INT(RollingTime(3-tnew,cr)),               &
     &                       INT(time(ng)),                             &
     &                       INT(RollingTime(tnew,cr))
 10           FORMAT (2x,'ng',2x,'cr',2x,'dg',2x,'rg',5x,'iic',5x,'iic',&
     &                2x,'told',2x,'tnew',2x,'Tindex',2x,'Tindex',      &
     &                9x,'time',8x,'time',8x,'time',8x,'time',/,        &
     &                20x,'(dg)',4x,'(rg)',18x,'2D',6x,'3D',9x,'(rg)',  &
     &                8x,'told',8x,'(ng)',8x,'tnew',/)
 20           FORMAT (4(1x,i3),2(1x,i7),2(2x,i4),2(4x,i4),1x,4(2x,i10))
              CALL my_flush (102)
            END IF
          END IF
# endif
!
!  Extract free-surface.
!
# ifdef SOLVE3D
!>        CALL get_contact2d (dg, model, tile,                          &
!>   &                        r2dvar, 'Zt_avg1',                        &
!>   &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        COUPLING(dg) % tl_Zt_avg1,                &
!>   &                        REFINED(cr) % tl_zeta(:,:,tnew))
!>
          CALL ad_get_contact2d (dg, model, tile,                       &
     &                           r2dvar, 'Zt_avg1',                     &
     &                           cr, Rcontact(cr)%Npoints, Rcontact,    &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           COUPLING(dg) % ad_Zt_avg1,             &
     &                           REFINED(cr) % ad_zeta(:,:,tnew))
# else
!>        CALL get_contact2d (dg, model, tile,                          &
!>   &                        r2dvar, 'zeta',                           &
!>   &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        OCEAN(dg) % tl_zeta(:,:,Tindex2d),        &
!>   &                        REFINED(cr) % tl_zeta(:,:,tnew))
!>
          CALL ad_get_contact2d (dg, model, tile,                       &
     &                           r2dvar, 'zeta',                        &
     &                           cr, Rcontact(cr)%Npoints, Rcontact,    &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           OCEAN(dg) % ad_zeta(:,:,Tindex2d),     &
     &                           REFINED(cr) % ad_zeta(:,:,tnew))
# endif
!
!  Extract 2D momentum components (ubar, vbar).
!
!>        CALL get_contact2d (dg, model, tile,                          &
!>   &                        u2dvar, Vname(1,idUbar),                  &
!>   &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        OCEAN(dg) % tl_ubar(:,:,Tindex2d),        &
!>   &                        REFINED(cr) % tl_ubar(:,:,tnew))
!>
          CALL ad_get_contact2d (dg, model, tile,                       &
     &                           u2dvar, Vname(1,idUbar),               &
     &                           cr, Ucontact(cr)%Npoints, Ucontact,    &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           OCEAN(dg) % ad_ubar(:,:,Tindex2d),     &
     &                           REFINED(cr) % ad_ubar(:,:,tnew))


!>        CALL get_contact2d (dg, model, tile,                          &
!>   &                        v2dvar, Vname(1,idVbar),                  &
!>   &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        OCEAN(dg) % tl_vbar(:,:,Tindex2d),        &
!>   &                        REFINED(cr) % tl_vbar(:,:,tnew))
!>
          CALL ad_get_contact2d (dg, model, tile,                       &
     &                           v2dvar, Vname(1,idVbar),               &
     &                           cr, Vcontact(cr)%Npoints, Vcontact,    &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           OCEAN(dg) % ad_vbar(:,:,Tindex2d),     &
     &                           REFINED(cr) % ad_vbar(:,:,tnew))

# ifdef SOLVE3D
!
!  Extract time-averaged fluxes (DU_avg2, DV_avg2).  We will use latter
!  only the values at the finer grid physical boundary to impose mass
!  flux conservation in routine "put_refine2d".
!
          CALL get_persisted2d (dg, rg, model, tile,                    &
     &                          u2dvar, 'DU_avg2',                      &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % DU_avg2,                 &
     &                          REFINED(cr) % DU_avg2(:,:,tnew))
!>        CALL get_persisted2d (dg, rg, model, tile,                    &
!>   &                          u2dvar, 'DU_avg2',                      &
!>   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          COUPLING(dg) % tl_DU_avg2,              &
!>   &                          REFINED(cr) % tl_DU_avg2(:,:,tnew))
!>
          CALL ad_get_persisted2d (dg, rg, model, tile,                 &
     &                             u2dvar, 'DU_avg2',                   &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             COUPLING(dg) % ad_DU_avg2,           &
     &                             REFINED(cr) % ad_DU_avg2(:,:,tnew))

          CALL get_persisted2d (dg, rg, model, tile,                    &
     &                          v2dvar, 'DV_avg2',                      &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % DV_avg2,                 &
     &                          REFINED(cr) % DV_avg2(:,:,tnew))

!>        CALL get_persisted2d (dg, rg, model, tile,                    &
!>   &                          v2dvar, 'DV_avg2',                      &
!>   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          COUPLING(dg) % tl_DV_avg2,              &
!>   &                          REFINED(cr) % tl_DV_avg2(:,:,tnew))
!>
          CALL ad_get_persisted2d (dg, rg, model, tile,                 &
     &                             v2dvar, 'DV_avg2',                   &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             COUPLING(dg) % ad_DV_avg2,           &
     &                             REFINED(cr) % ad_DV_avg2(:,:,tnew))
!
!  Tracer-type variables.
!
          DO itrc=1,NT(dg)
!>          CALL get_contact3d (dg, model, tile,                        &
!>   &                          r3dvar, Vname(1,idTvar(itrc)),          &
!>   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!>   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!>   &                          OCEAN(dg) % tl_t(:,:,:,Tindex3d,itrc),  &
!>   &                          REFINED(cr) % tl_t(:,:,:,tnew,itrc))
!>
            CALL ad_get_contact3d (dg, model, tile,                     &
     &                             r3dvar, Vname(1,idTvar(itrc)),       &
     &                             cr, Rcontact(cr)%Npoints, Rcontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(dg),        &
     &                             OCEAN(dg) % ad_t(:,:,:,Tindex3d,     &
     &                                              itrc),              &
     &                             REFINED(cr) % ad_t(:,:,:,tnew,itrc))
          END DO
!
!  Extract 3D momentum components (u, v).
!
!>        CALL get_contact3d (dg, model, tile,                          &
!>   &                        u3dvar, Vname(1,idUvel),                  &
!>   &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
!>   &                        LBi, UBi, LBj, UBj, 1, N(dg),             &
!>   &                        OCEAN(dg) % tl_u(:,:,:,Tindex3d),         &
!>   &                        REFINED(cr) % tl_u(:,:,:,tnew))
!>
          CALL ad_get_contact3d (dg, model, tile,                       &
     &                           u3dvar, Vname(1,idUvel),               &
     &                           cr, Ucontact(cr)%Npoints, Ucontact,    &
     &                           LBi, UBi, LBj, UBj, 1, N(dg),          &
     &                           OCEAN(dg) % ad_u(:,:,:,Tindex3d),      &
     &                           REFINED(cr) % ad_u(:,:,:,tnew))

!>        CALL get_contact3d (dg, model, tile,                          &
!>   &                        v3dvar, Vname(1,idVvel),                  &
!>   &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
!>   &                        LBi, UBi, LBj, UBj, 1, N(dg),             &
!>   &                        OCEAN(dg) % tl_v(:,:,:,Tindex3d),         &
!>   &                        REFINED(cr) % tl_v(:,:,:,tnew))
!>
          CALL ad_get_contact3d (dg, model, tile,                       &
     &                           v3dvar, Vname(1,idVvel),               &
     &                           cr, Vcontact(cr)%Npoints, Vcontact,    &
     &                           LBi, UBi, LBj, UBj, 1, N(dg),          &
     &                           OCEAN(dg) % ad_v(:,:,:,Tindex3d),      &
     &                           REFINED(cr) % ad_v(:,:,:,tnew))
# endif
        END IF
      END DO

      RETURN
      END SUBROUTINE ad_get_refine
!
      SUBROUTINE ad_put_composite (ng, model, isection, tile)
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
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d, ad_mp_exchange4d
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
# ifdef DISTRIBUTE
!>          CALL mp_exchange2d (rg, tile, model, 2,                     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          FORCES(rg) % tl_bustr,                  &
!>   &                          FORCES(rg) % tl_bvstr)
!>
            CALL ad_mp_exchange2d (rg, tile, model, 2,                  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             FORCES(rg) % ad_bustr,               &
     &                             FORCES(rg) % ad_bvstr)
# endif
!>          CALL put_contact2d (rg, model, tile,                        &
!>   &                          u2dvar, Vname(1,idUbms),                &
!>   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
# ifdef MASKING
!>   &                          GRID(rg) % umask,                       &
# endif
!>   &                          COMPOSITE(cr) % tl_bustr,               &
!>   &                          FORCES(rg) % tl_bustr)
!>
            CALL ad_put_contact2d (rg, model, tile,                     &
     &                             u2dvar, Vname(1,idUbms),             &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
# ifdef MASKING
     &                             GRID(rg) % umask,                    &
# endif
     &                             COMPOSITE(cr) % ad_bustr,            &
     &                             FORCES(rg) % ad_bustr)
!>          CALL put_contact2d (rg, model, tile,                        &
!>   &                          v2dvar, Vname(1,idVbms),                &
!>   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
# ifdef MASKING
!>   &                          GRID(rg) % vmask,                       &
# endif
!>   &                          COMPOSITE(cr) % tl_bvstr,               &
!>   &                          FORCES(rg) % tl_bvstr)
!>
            CALL ad_put_contact2d (rg, model, tile,                     &
     &                             v2dvar, Vname(1,idVbms),             &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
# ifdef MASKING
     &                             GRID(rg) % vmask,                    &
# endif
     &                             COMPOSITE(cr) % ad_bvstr,            &
     &                             FORCES(rg) % ad_bvstr)
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
# ifdef DISTRIBUTE
!>            CALL mp_exchange2d (rg, tile, model, 1,                   &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            NghostPoints,                         &
!>   &                            EWperiodic(rg), NSperiodic(rg),       &
!>   &                            OCEAN(rg) % tl_zeta(:,:,Tindex))
!>
              CALL ad_mp_exchange2d (rg, tile, model, 1,                &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               EWperiodic(rg), NSperiodic(rg),    &
     &                               OCEAN(rg) % ad_zeta(:,:,Tindex))
# endif
!>            CALL put_contact2d (rg, model, tile,                      &
!>   &                            r2dvar, Vname(1,idFsur),              &
!>   &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
!>   &                            LBi, UBi, LBj, UBj,                   &
# ifdef MASKING
!>   &                            GRID(rg) % rmask,                     &
# endif
!>   &                            COMPOSITE(cr) % tl_zeta(:,:,rec),     &
!>   &                            OCEAN(rg) % tl_zeta(:,:,Tindex))
!>
              CALL ad_put_contact2d (rg, model, tile,                   &
     &                               r2dvar, Vname(1,idFsur),           &
     &                               cr, Rcontact(cr)%Npoints, Rcontact,&
     &                               LBi, UBi, LBj, UBj,                &
# ifdef MASKING
     &                               GRID(rg) % rmask,                  &
# endif
     &                               COMPOSITE(cr) % ad_zeta(:,:,rec),  &
     &                               OCEAN(rg) % ad_zeta(:,:,Tindex))
            END DO
          END IF
!
!  Process free-surface equation rigth-hand-side (rzeta) term.
!
          IF (isection.eq.n2dPS) THEN
            Tindex=1
# ifdef DISTRIBUTE
!>          CALL mp_exchange2d (rg, tile, model, 1,                     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          OCEAN(rg) % tl_rzeta(:,:,Tindex))
!>
            CALL ad_mp_exchange2d (rg, tile, model, 1,                  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             OCEAN(rg) % ad_rzeta(:,:,Tindex))
# endif
!>          CALL put_contact2d (rg, model, tile,                        &
!>   &                          r2dvar, Vname(1,idRzet),                &
!>   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
# ifdef MASKING
!>   &                          GRID(rg) % rmask,                       &
# endif
!>   &                          COMPOSITE(cr) % tl_rzeta,               &
!>   &                          OCEAN(rg) % tl_rzeta(:,:,Tindex))
!>
            CALL ad_put_contact2d (rg, model, tile,                     &
     &                             r2dvar, Vname(1,idRzet),             &
     &                             cr, Rcontact(cr)%Npoints, Rcontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
# ifdef MASKING
     &                             GRID(rg) % rmask,                    &
# endif
     &                             COMPOSITE(cr) % ad_rzeta,            &
     &                             OCEAN(rg) % ad_rzeta(:,:,Tindex))
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
# ifdef DISTRIBUTE
!>            CALL mp_exchange2d (rg, tile, model, 2,                   &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            NghostPoints,                         &
!>   &                            EWperiodic(rg), NSperiodic(rg),       &
!>   &                            OCEAN(rg) % tl_ubar(:,:,Tindex),      &
!>   &                            OCEAN(rg) % tl_vbar(:,:,Tindex))
!>
              CALL ad_mp_exchange2d (rg, tile, model, 2,                &
     &                               LBi, UBi, LBj, UBj,                &
     &                               NghostPoints,                      &
     &                               EWperiodic(rg), NSperiodic(rg),    &
     &                               OCEAN(rg) % ad_ubar(:,:,Tindex),   &
     &                               OCEAN(rg) % ad_vbar(:,:,Tindex))
# endif
!>            CALL put_contact2d (rg, model, tile,                      &
!>   &                            u2dvar, Vname(1,idUbar),              &
!>   &                            cr, Ucontact(cr)%Npoints, Ucontact,   &
!>   &                            LBi, UBi, LBj, UBj,                   &
# ifdef MASKING
!>   &                            GRID(rg) % umask,                     &
# endif
!>   &                            COMPOSITE(cr) % tl_ubar(:,:,rec),     &
!>   &                            OCEAN(rg) % tl_ubar(:,:,Tindex))
!>
              CALL ad_put_contact2d (rg, model, tile,                   &
     &                               u2dvar, Vname(1,idUbar),           &
     &                               cr, Ucontact(cr)%Npoints, Ucontact,&
     &                               LBi, UBi, LBj, UBj,                &
# ifdef MASKING
     &                               GRID(rg) % umask,                  &
# endif
     &                               COMPOSITE(cr) % ad_ubar(:,:,rec),  &
     &                               OCEAN(rg) % ad_ubar(:,:,Tindex))
!>            CALL put_contact2d (rg, model, tile,                      &
!>   &                            v2dvar, Vname(1,idVbar),              &
!>   &                            cr, Vcontact(cr)%Npoints, Vcontact,   &
!>   &                            LBi, UBi, LBj, UBj,                   &
# ifdef MASKING
!>   &                            GRID(rg) % vmask,                     &
# endif
!>   &                            COMPOSITE(cr) % tl_vbar(:,:,rec),     &
!>   &                            OCEAN(rg) % tl_vbar(:,:,Tindex))
!>
              CALL ad_put_contact2d (rg, model, tile,                   &
     &                               v2dvar, Vname(1,idVbar),           &
     &                               cr, Vcontact(cr)%Npoints, Vcontact,&
     &                               LBi, UBi, LBj, UBj,                &
# ifdef MASKING
     &                               GRID(rg) % vmask,                  &
# endif
     &                               COMPOSITE(cr) % ad_vbar(:,:,rec),  &
     &                               OCEAN(rg) % ad_vbar(:,:,Tindex))
            END DO
          END IF

# ifdef SOLVE3D
!
!  Process time averaged free-surface (Zt_avg1) and 2D momentum fluxes
!  (DU_avg1, DV_avg1).
!
          IF (isection.eq.n2dfx) THEN
#  ifdef DISTRIBUTE
!>          CALL mp_exchange2d (rg, tile, model, 3,                     &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          COUPLING(rg) % tl_Zt_avg1,              &
!>   &                          COUPLING(rg) % tl_DU_avg1,              &
!>   &                          COUPLING(rg) % tl_DV_avg1)
!>
            CALL ad_mp_exchange2d (rg, tile, model, 3,                  &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             COUPLING(rg) % ad_Zt_avg1,           &
     &                             COUPLING(rg) % ad_DU_avg1,           &
     &                             COUPLING(rg) % ad_DV_avg1)
#  endif
!>          CALL put_contact2d (rg, model, tile,                        &
!>   &                          r2dvar, 'Zt_avg1',                      &
!>   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
#  ifdef MASKING
!>   &                          GRID(rg) % rmask,                       &
#  endif
!>   &                          COMPOSITE(cr) % tl_Zt_avg1,             &
!>   &                          COUPLING(rg) % tl_Zt_avg1)a
!>
            CALL ad_put_contact2d (rg, model, tile,                     &
     &                             r2dvar, 'Zt_avg1',                   &
     &                             cr, Rcontact(cr)%Npoints, Rcontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
#  ifdef MASKING
     &                             GRID(rg) % rmask,                    &
#  endif
     &                             COMPOSITE(cr) % ad_Zt_avg1,          &
     &                             COUPLING(rg) % ad_Zt_avg1)
!>          CALL put_contact2d (rg, model, tile,                        &
!>   &                          u2dvar, Vname(1,idUfx1),                &
!>   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
#  ifdef MASKING
!>   &                          GRID(rg) % umask,                       &
#  endif
!>   &                          COMPOSITE(cr) % tl_DU_avg1,             &
!>   &                          COUPLING(rg) % tl_DU_avg1)
!>
            CALL ad_put_contact2d (rg, model, tile,                     &
     &                             u2dvar, Vname(1,idUfx1),             &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
#  ifdef MASKING
     &                             GRID(rg) % umask,                    &
#  endif
     &                             COMPOSITE(cr) % ad_DU_avg1,          &
     &                             COUPLING(rg) % ad_DU_avg1)
!>          CALL put_contact2d (rg, model, tile,                        &
!>   &                          v2dvar, Vname(1,idVfx1),                &
!>   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!>   &                          LBi, UBi, LBj, UBj,                     &
#  ifdef MASKING
!>   &                          GRID(rg) % vmask,                       &
#  endif
!>   &                          COMPOSITE(cr) % tl_DV_avg1,             &
!>   &                          COUPLING(rg) % tl_DV_avg1)
!>
            CALL ad_put_contact2d (rg, model, tile,                     &
     &                             v2dvar, Vname(1,idVfx1),             &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj,                  &
#  ifdef MASKING
     &                             GRID(rg) % vmask,                    &
#  endif
     &                             COMPOSITE(cr) % ad_DV_avg1,          &
     &                             COUPLING(rg) % ad_DV_avg1)
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
#   ifdef DISTRIBUTE
!>          CALL mp_exchange4d (rg, tile, model, 1,                     &
!>   &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
!>   &                          1, NT(rg),                              &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          OCEAN(rg) % tl_t(:,:,:,Tindex,:))
!>
              CALL ad_mp_exchange4d (rg, tile, model, 1,                &
     &                               LBi, UBi, LBj, UBj, 1, N(rg),      &
     &                               1, NT(rg),                         &
     &                               NghostPoints,                      &
     &                               EWperiodic(rg), NSperiodic(rg),    &
     &                               OCEAN(rg) % ad_t(:,:,:,Tindex,:))
#   endif
!>            CALL tl_put_contact3d (rg, model, tile,                   &
!>   &                               r3dvar, Vname(1,idTvar(itrc)),     &
!>   &                               cr, Rcontact(cr)%Npoints, Rcontact,&
!>   &                               LBi, UBi, LBj, UBj, 1, N(rg),      &
#   ifdef MASKING
!>   &                               GRID(rg) % rmask,                  &
#   endif
!>   &                               COMPOSITE(cr) % t(:,:,:,itrc),     &
!>   &                               COMPOSITE(cr) % tl_t(:,:,:,itrc),  &
!>   &                               OCEAN(rg) % tl_t(:,:,:,Tindex,     &
!>   &                                                itrc))
!>
              CALL ad_put_contact3d (rg, model, tile,                   &
     &                               r3dvar, Vname(1,idTvar(itrc)),     &
     &                               cr, Rcontact(cr)%Npoints, Rcontact,&
     &                               LBi, UBi, LBj, UBj, 1, N(rg),      &
#   ifdef MASKING
     &                               GRID(rg) % rmask,                  &
#   endif
     &                               COMPOSITE(cr) % t(:,:,:,itrc),     &
     &                               COMPOSITE(cr) % ad_t(:,:,:,itrc),  &
     &                               OCEAN(rg) % ad_t(:,:,:,Tindex,     &
     &                                                itrc))
            END DO
          END IF
#  endif
!
!  Process 3D momentum (u, v) at the appropriate time-index.
!
          IF ((isection.eq.n3dIC).or.                                   &
     &        (isection.eq.n3duv)) THEN
            Tindex=nnew(rg)
#  ifdef DISTRIBUTE
!>         CALL mp_exchange3d (rg, tile, model, 2,                      &
!>   &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          OCEAN(rg) % tl_u(:,:,:,Tindex),         &
!>   &                          OCEAN(rg) % tl_v(:,:,:,Tindex))
!>
            CALL ad_mp_exchange3d (rg, tile, model, 2,                  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             OCEAN(rg) % ad_u(:,:,:,Tindex),      &
     &                             OCEAN(rg) % ad_v(:,:,:,Tindex))
#  endif
!>          CALL tl_put_contact3d (rg, model, tile,                     &
!>   &                             u3dvar, Vname(1,idUvel),             &
!>   &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
!>   &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
!>   &                             GRID(rg) % umask,                    &
#  endif
!>   &                             COMPOSITE(cr) % u,                   &
!>   &                             COMPOSITE(cr) % tl_u,                &
!>   &                             OCEAN(rg) % tl_u(:,:,:,Tindex))
!>
            CALL ad_put_contact3d (rg, model, tile,                     &
     &                             u3dvar, Vname(1,idUvel),             &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
     &                             GRID(rg) % umask,                    &
#  endif
     &                             COMPOSITE(cr) % u,                   &
     &                             COMPOSITE(cr) % ad_u,                &
     &                             OCEAN(rg) % ad_u(:,:,:,Tindex))
!>          CALL tl_put_contact3d (rg, model, tile,                     &
!>   &                             v3dvar, Vname(1,idVvel),             &
!>   &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
!>   &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
!>   &                             GRID(rg) % vmask,                    &
#  endif
!>   &                             COMPOSITE(cr) % v,                   &
!>   &                             COMPOSITE(cr) % tl_v,                &
!>   &                             OCEAN(rg) % tl_v(:,:,:,Tindex))
!>
            CALL ad_put_contact3d (rg, model, tile,                     &
     &                             v3dvar, Vname(1,idVvel),             &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
     &                             GRID(rg) % vmask,                    &
#  endif
     &                             COMPOSITE(cr) % v,                   &
     &                             COMPOSITE(cr) % ad_v,                &
     &                             OCEAN(rg) % ad_v(:,:,:,Tindex))
          END IF
!
!  Process 3D momentum fluxes (Huon, Hvom).
!
          IF (isection.eq.n3duv) THEN
#  ifdef DISTRIBUTE
!>          CALL mp_exchange3d (rg, tile, model, 2,                     &
!>   &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          GRID(rg) % tl_Huon,                     &
!>   &                          GRID(rg) % tl_Hvom)
!>
            CALL ad_mp_exchange3d (rg, tile, model, 2,                  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             GRID(rg) % ad_Huon,                  &
     &                             GRID(rg) % ad_Hvom)
#  endif
!>          CALL tl_put_contact3d (rg, model, tile,                     &
!>   &                             u3dvar, 'Huon',                      &
!>   &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
!>   &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
!>   &                             GRID(rg) % umask,                    &
#  endif
!>   &                             COMPOSITE(cr) % Huon,                &
!>   &                             COMPOSITE(cr) % tl_Huon,             &
!>   &                             GRID(rg) % tl_Huon)
!>
            CALL ad_put_contact3d (rg, model, tile,                     &
     &                             u3dvar, 'Huon',                      &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
     &                             GRID(rg) % umask,                    &
#  endif
     &                             COMPOSITE(cr) % Huon,                &
     &                             COMPOSITE(cr) % ad_Huon,             &
     &                             GRID(rg) % ad_Huon)
!>          CALL tl_put_contact3d (rg, model, tile,                     &
!>   &                             v3dvar, 'Hvom',                      &
!>   &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
!>   &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
!>   &                             GRID(rg) % vmask,                    &
#  endif
!>   &                             COMPOSITE(cr) % Hvom,                &
!>   &                             COMPOSITE(cr) % tl_Hvom,             &
!>   &                             GRID(rg) % tl_Hvom)
!>
            CALL ad_put_contact3d (rg, model, tile,                     &
     &                             v3dvar, 'Hvom',                      &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
     &                             GRID(rg) % vmask,                    &
#  endif
     &                             COMPOSITE(cr) % Hvom,                &
     &                             COMPOSITE(cr) % ad_Hvom,             &
     &                             GRID(rg) % ad_Hvom)
          END IF
# endif

        END IF
      END DO CR_LOOP

      RETURN
      END SUBROUTINE ad_put_composite
!
      SUBROUTINE ad_put_refine (ng, model, tile, LputFsur)
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
      integer :: dg, rg, cr, nrec, rec, tnew, told
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
!  Update adjoint rolling time indices. The contact data is stored in
!  two time levels.
!
          IF (.not.LputFsur) THEN
            IF (time(dg).eq.time(rg)) THEN
              RollingIndex(cr)=3-RollingIndex(cr)
            END IF
            tnew=RollingIndex(cr)
            told=3-tnew
!           told=RollingIndex(cr)
!           tnew=3-told
!           IF (time(dg).eq.time(rg)) THEN
!             RollingIndex(cr)=3-RollingIndex(cr)
!           END IF
            RollingTime(tnew,cr)=time(dg)+dt(dg)
            RollingTime(told,cr)=time(dg)
          END IF
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
!>          CALL tl_put_refine2d (ng, dg, cr, model, tile, LputFsur,    &
!>   &                            LBi, UBi, LBj, UBj)
!>
            CALL ad_put_refine2d (ng, dg, cr, model, tile, LputFsur,    &
     &                            LBi, UBi, LBj, UBj)
          ELSE
# ifdef SOLVE3D
!
!  Fill 3D state variables contact points.
!
!>          CALL tl_put_refine3d (ng, dg, cr, model, tile,              &
!>   &                            LBi, UBi, LBj, UBj)
!>
            CALL ad_put_refine3d (ng, dg, cr, model, tile,              &
     &                            LBi, UBi, LBj, UBj)
# endif
!
!  Fill other 2D state variables (like momentum) contact points.
!
!>          CALL tl_put_refine2d (ng, dg, cr, model, tile, LputFsur,    &
!>   &                            LBi, UBi, LBj, UBj)
!>
            CALL ad_put_refine2d (ng, dg, cr, model, tile, LputFsur,    &
     &                            LBi, UBi, LBj, UBj)
          END IF
        END IF
      END DO

      RETURN
      END SUBROUTINE ad_put_refine

# ifdef SOLVE3D
!
      SUBROUTINE ad_correct_tracer (ng, ngf, model, tile)
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
      CALL ad_correct_tracer_tile (ng, ngf, model, tile,                &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS)
      RETURN
!
      END SUBROUTINE ad_correct_tracer
!
!***********************************************************************
      SUBROUTINE ad_correct_tracer_tile (ngc, ngf, model, tile,         &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS)
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
      USE mp_exchange_mod, ONLY : ad_mp_exchange4d
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
      real(r8) :: ad_TFC, ad_TFF, ad_Tvalue, ad_cff, adfac

      real(r8) :: Dinv(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ad_Dinv(IminS:ImaxS,JminS:JmaxS)

!
!  Clear adjoint constants.
!
      ad_TFC=0.0_r8
      ad_TFF=0.0_r8
      ad_Tvalue=0.0_r8
      ad_cff=0.0_r8
      ad_Dinv=0.0_r8

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
!>    CALL mp_exchange4d (ngc, tile, model, 1,                          &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ngc),                &
!>   &                    1, NT(ngc),                                   &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ngc), NSperiodic(ngc),             &
!>   &                    OCEAN(ngc)%tl_t(:,:,:,Tindex,:))
!>
      CALL ad_mp_exchange4d (ngc, tile, model, 1,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ngc),             &
     &                       1, NT(ngc),                                &
     &                       NghostPoints,                              &
     &                       EWperiodic(ngc), NSperiodic(ngc),          &
     &                       OCEAN(ngc)%ad_t(:,:,:,Tindex,:))
#  endif
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
!  Adjoint Finer grid northern boundary.
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
!  Zeroth order correction to fine grid time integral.
!
            TFF=TFF*dt(ngc)/dt(ngf)
!
            cff=GRID(ngc)%pm(Ibc,Jbc)*                                  &
     &          GRID(ngc)%pn(Ibc,Jbc)*                                  &
     &          Dinv(Ibc,Jbc)
            DO k=1,N(ngc)
!>            OCEAN(ngc)%tl_t(Ibc,Jbc,k,Tindex,itrc)=tl_Tvalue
!>
              ad_Tvalue=ad_Tvalue+                                      &
     &                  OCEAN(ngc)%ad_t(Ibc,Jbc,k,Tindex,itrc)
              OCEAN(ngc)%ad_t(Ibc,Jbc,k,Tindex,itrc)=0.0_r8
#  ifdef MASKING
!>            tl_Tvalue=tl_Tvalue*GRID(ngc)%rmask(Ibc,Jbc)
!>
              ad_Tvalue=ad_Tvalue*GRID(ngc)%rmask(Ibc,Jbc)

#  endif
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
!>              tl_Tvalue=tl_Tvalue-                                    &
!>   &                    dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc)*  &
!>   &                    tl_Tvalue
!>
                ad_Tvalue=ad_Tvalue*                                    &
     &                    (1.0_r8-dt(ngc)*                              &
     &                            CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc))
              END IF
!>            tl_Tvalue=(0.5_r8-                                        &
!>   &                   SIGN(0.5_r8,                                   &
!>   &                        -(OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-    &
!>   &                        cff*(TFF-TFC))))*                         &
!>   &                   (OCEAN(ngc)%tl_t(Ibc,Jbc,k,Tindex,itrc)-       &
!>   &                   tl_cff*(TFF-TFC)-cff*(tl_TFF-tl_TFC))
!>
              adfac=(0.5_r8-                                            &
     &               SIGN(0.5_r8,                                       &
     &                    -(OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-        &
     &                    cff*(TFF-TFC))))*ad_Tvalue
              OCEAN(ngc)%ad_t(Ibc,Jbc,k,Tindex,itrc)=                   &
     &                   OCEAN(ngc)%ad_t(Ibc,Jbc,k,Tindex,itrc)+adfac
              ad_cff=ad_cff-(TFF-TFC)*adfac
              ad_TFF=ad_TFF-cff*adfac
              ad_TFC=ad_TFC+cff*adfac
              ad_Tvalue=0.0_r8
            END DO
!
!  Correct coarse grid tracer at the finer grid northern boundary.
!
!>          tl_cff=GRID(ngc)%pm(Ibc,Jbc)*                               &
!>   &             GRID(ngc)%pn(Ibc,Jbc)*                               &
!>   &             tl_Dinv(Ibc,Jbc)
!>
            ad_Dinv(Ibc,Jbc)=ad_Dinv(Ibc,Jbc)+                          &
     &                       GRID(ngc)%pm(Ibc,Jbc)*                     &
     &                       GRID(ngc)%pn(Ibc,Jbc)*ad_cff
            ad_cff=0.0_r8
!
!  Zeroth order correction to fine grid time integral.
!
!>         tl_TFF=tl_TFF*dt(ngc)/dt(ngf)
!>
           ad_TFF=ad_TFF*dt(ngc)/dt(ngf)
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal I-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the I-edge.
!
            Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
            DO isum=-half,half
              Ibf=Iedge+isum
              DO k=1,N(ngf)
!>              tl_TFF=tl_TFF+                                          &
!>   &                 BRY_CONTACT(inorth,dgcr)%tl_Tflux(Ibf,k,itrc)
!>
                BRY_CONTACT(inorth,dgcr)%ad_Tflux(Ibf,k,itrc)=          &
     &              BRY_CONTACT(inorth,dgcr)%ad_Tflux(Ibf,k,itrc)+ad_TFF
              END DO
            END DO
!>          tl_TFF=0.0_r8
!>
            ad_TFF=0.0_r8
!
!  Sum vertically coarse grid horizontal advective tracer flux,
!  Hz*v*T/m, from last time-step.
!
            DO k=1,N(ngc)
!>            tl_TFC=tl_TFC+                                            &
!>   &                  BRY_CONTACT(inorth,rgcr)%tl_Tflux(Ibc,k,itrc)
!>
              BRY_CONTACT(inorth,rgcr)%ad_Tflux(Ibc,k,itrc)=            &
     &            BRY_CONTACT(inorth,rgcr)%ad_Tflux(Ibc,k,itrc)+ad_TFC
            END DO
!>          tl_TFC=0.0_r8
!>
            ad_TFC=0.0_r8
          END IF
        END DO

!
!-----------------------------------------------------------------------
!  Adjoint Finer grid southern boundary.
!-----------------------------------------------------------------------
!
        Jbc=J_bottom(ngf)
        Ibc_min=I_left(ngf)
        Ibc_max=I_right(ngf)-1            ! interior points, no bottom
!                                           right corner
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
!  Zeroth order correction to fine grid time integral (RIL, 2016).
!
            TFF=TFF*dt(ngc)/dt(ngf)

            cff=GRID(ngc)%pm(Ibc,Jbc-1)*                                &
     &          GRID(ngc)%pn(Ibc,Jbc-1)*                                &
     &          Dinv(Ibc,Jbc-1)

            DO k=1,N(ngc)
!>            OCEAN(ngc)%tl_t(Ibc,Jbc-1,k,Tindex,itrc)=tl_Tvalue
!>
              ad_Tvalue=ad_Tvalue+                                      &
     &                  OCEAN(ngc)%ad_t(Ibc,Jbc-1,k,Tindex,itrc)
              OCEAN(ngc)%ad_t(Ibc,Jbc-1,k,Tindex,itrc)=0.0_r8
#  ifdef MASKING
!>            tl_Tvalue=tl_Tvalue*GRID(ngc)%rmask(Ibc,Jbc-1)
!>
              ad_Tvalue=ad_Tvalue*GRID(ngc)%rmask(Ibc,Jbc-1)
#  endif
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
!>              tl_Tvalue=tl_Tvalue-                                    &
!>   &                    dt(ngc)*                                      &
!>   &                    CLIMA(ngc)%Tnudgcof(Ibc,Jbc-1,k,itrc)*        &
!>   &                    tl_Tvalue
!>
                ad_Tvalue=ad_Tvalue*                                    &
     &                    (1.0_r8-dt(ngc)*                              &
     &                            CLIMA(ngc)%Tnudgcof(Ibc,Jbc-1,k,itrc))
              END IF
!>            tl_Tvalue=(0.5_r8-                                        &
!>   &                   SIGN(0.5_r8,                                   &
!>   &                        -(OCEAN(ngc)%t(Ibc,Jbc-1,k,Tindex,itrc)-  &
!>   &                        cff*(TFF-TFC))))*                         &
!>   &                   (OCEAN(ngc)%tl_t(Ibc,Jbc-1,k,Tindex,itrc)-     &
!>   &                   tl_cff*(TFF-TFC)-cff*(tl_TFF-tl_TFC))
!>
              adfac=(0.5_r8-                                            &
     &               SIGN(0.5_r8,                                       &
     &                    -(OCEAN(ngc)%t(Ibc,Jbc-1,k,Tindex,itrc)-      &
     &                    cff*(TFF-TFC))))*ad_Tvalue
              OCEAN(ngc)%ad_t(Ibc,Jbc-1,k,Tindex,itrc)=                 &
     &                   OCEAN(ngc)%ad_t(Ibc,Jbc-1,k,Tindex,itrc)+adfac
              ad_cff=ad_cff-(TFF-TFC)*adfac
              ad_TFF=ad_TFF-cff*adfac
              ad_TFC=ad_TFC+cff*adfac
              ad_Tvalue=0.0_r8
            END DO
!
!  Correct coarse grid tracer at the finer grid southern boundary.
!
!>          tl_cff=GRID(ngc)%pm(Ibc,Jbc-1)*                             &
!>   &             GRID(ngc)%pn(Ibc,Jbc-1)*                             &
!>   &             tl_Dinv(Ibc,Jbc-1)
!>
            ad_Dinv(Ibc,Jbc-1)=ad_Dinv(Ibc,Jbc-1)+                      &
     &                         GRID(ngc)%pm(Ibc,Jbc-1)*                 &
     &                         GRID(ngc)%pn(Ibc,Jbc-1)*ad_cff
            ad_cff=0.0_r8
!
!  Zeroth order correction to fine grid time integral (RIL, 2016).
!
!>          tl_TFF=tl_TFF*dt(ngc)/dt(ngf)
!>
            ad_TFF=ad_TFF*dt(ngc)/dt(ngf)
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal I-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the I-edge.
!
            Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
            DO isum=-half,half
              Ibf=Iedge+isum
              DO k=1,N(ngf)
!>              tl_TFF=tl_TFF+                                          &
!>   &                 BRY_CONTACT(isouth,dgcr)%tl_Tflux(Ibf,k,itrc)
!>
                BRY_CONTACT(isouth,dgcr)%ad_Tflux(Ibf,k,itrc)=          &
     &              BRY_CONTACT(isouth,dgcr)%ad_Tflux(Ibf,k,itrc)+ad_TFF
              END DO
            END DO
!>          tl_TFF=0.0_r8
!>
            ad_TFF=0.0_r8
!
!  Sum vertically coarse grid horizontal advective tracer flux,
!  Hz*v*T/m, from last time-step.
!
            DO k=1,N(ngc)
!>            tl_TFC=tl_TFC+                                            &
!>   &               BRY_CONTACT(isouth,rgcr)%tl_Tflux(Ibc,k,itrc)
!>
              BRY_CONTACT(isouth,rgcr)%ad_Tflux(Ibc,k,itrc)=            &
     &            BRY_CONTACT(isouth,rgcr)%ad_Tflux(Ibc,k,itrc)+ad_TFC
            END DO
!>          tl_TFC=0.0_r8
!>
            ad_TFC=0.0_r8
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
!  Zeroth order correction to fine grid time integral (RIL, 2016).
!
            TFF=TFF*dt(ngc)/dt(ngf)
!
!  Correct coarse grid tracer at the finer grid eastern boundary.
!
            cff=GRID(ngc)%pm(Ibc,Jbc)*                                  &
     &          GRID(ngc)%pn(Ibc,Jbc)*                                  &
     &          Dinv(Ibc,Jbc)
            DO k=1,N(ngc)
!>            OCEAN(ngc)%tl_t(Ibc,Jbc,k,Tindex,itrc)=tl_Tvalue
!>
              ad_Tvalue=ad_Tvalue+                                      &
     &                  OCEAN(ngc)%ad_t(Ibc,Jbc,k,Tindex,itrc)
              OCEAN(ngc)%ad_t(Ibc,Jbc,k,Tindex,itrc)=0.0_r8
#  ifdef MASKING
!>            tl_Tvalue=tl_Tvalue*GRID(ngc)%rmask(Ibc,Jbc)
!>
              ad_Tvalue=ad_Tvalue*GRID(ngc)%rmask(Ibc,Jbc)
#  endif
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
!>              tl_Tvalue=tl_Tvalue-                                    &
!>   &                    dt(ngc)*                                      &
!>   &                    CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc)*          &
!>   &                    tl_Tvalue
!>
                ad_Tvalue=ad_Tvalue*                                    &
     &                    (1.0_r8-dt(ngc)*                              &
     &                            CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc))
              END IF
!>            tl_Tvalue=(0.5_r8-                                        &
!>   &                   SIGN(0.5_r8,                                   &
!>   &                        -(OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-    &
!>   &                        cff*(TFF-TFC))))*                         &
!>   &                   (OCEAN(ngc)%tl_t(Ibc,Jbc,k,Tindex,itrc)-       &
!>   &                    tl_cff*(TFF-TFC)-cff*(tl_TFF-tl_TFC))
!>
              adfac=(0.5_r8-                                            &
     &               SIGN(0.5_r8,                                       &
     &                    -(OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-        &
     &                    cff*(TFF-TFC))))*ad_Tvalue
              OCEAN(ngc)%ad_t(Ibc,Jbc,k,Tindex,itrc)=                   &
     &                   OCEAN(ngc)%ad_t(Ibc,Jbc,k,Tindex,itrc)+adfac
              ad_cff=ad_cff-(TFF-TFC)*adfac
              ad_TFF=ad_TFF-cff*adfac
              ad_TFC=ad_TFC+cff*adfac
              ad_Tvalue=0.0_r8
            END DO
!
!  Correct coarse grid tracer at the finer grid eastern boundary.
!
!>          tl_cff=GRID(ngc)%pm(Ibc,Jbc)*                               &
!>   &             GRID(ngc)%pn(Ibc,Jbc)*                               &
!>   &             tl_Dinv(Ibc,Jbc)
!>
            ad_Dinv(Ibc,Jbc)=ad_Dinv(Ibc,Jbc)+                          &
     &                       GRID(ngc)%pm(Ibc,Jbc)*                     &
     &                       GRID(ngc)%pn(Ibc,Jbc)*ad_cff
            ad_cff=0.0_r8
!
!  Zeroth order correction to fine grid time integral (RIL, 2016).
!
!>          tl_TFF=tl_TFF*dt(ngc)/dt(ngf)
!>
            ad_TFF=ad_TFF*dt(ngc)/dt(ngf)
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal J-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the J-edge.
!
            Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
            DO jsum=-half,half
              Jbf=Jedge+jsum
              DO k=1,N(ngf)
!>              tl_TFF=tl_TFF+                                          &
!>   &                 BRY_CONTACT(ieast,dgcr)%tl_Tflux(Jbf,k,itrc)
!>
                BRY_CONTACT(ieast,dgcr)%ad_Tflux(Jbf,k,itrc)=           &
     &              BRY_CONTACT(ieast,dgcr)%ad_Tflux(Jbf,k,itrc)+ad_TFF
              END DO
            END DO
!>          tl_TFF=0.0_r8
!>
            ad_TFF=0.0_r8
!
!  Sum vertically coarse grid horizontal advective tracer flux,
!  Hz*u*T/n, from last time-step.
!
            DO k=1,N(ngc)
!>            tl_TFC=tl_TFC+                                            &
!>   &               BRY_CONTACT(ieast,rgcr)%tl_Tflux(Jbc,k,itrc)
!>
              BRY_CONTACT(ieast,rgcr)%ad_Tflux(Jbc,k,itrc)=             &
     &            BRY_CONTACT(ieast,rgcr)%ad_Tflux(Jbc,k,itrc)+ad_TFC
            END DO
!>          tl_TFC=0.0_r8
!>
            ad_TFC=0.0_r8
          END IF
        END DO
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
!  Zeroth order correction to fine grid time integral (RIL, 2016).
!
            TFF=TFF*dt(ngc)/dt(ngf)
!
!  Correct coarse grid tracer at the finer grid western boundary.
!
            cff=GRID(ngc)%pm(Ibc-1,Jbc)*                                &
     &          GRID(ngc)%pn(Ibc-1,Jbc)*                                &
     &          Dinv(Ibc-1,Jbc)
            DO k=1,N(ngc)
!>            OCEAN(ngc)%tl_t(Ibc-1,Jbc,k,Tindex,itrc)=tl_Tvalue
!>
              ad_Tvalue=ad_Tvalue+                                      &
     &                  OCEAN(ngc)%ad_t(Ibc-1,Jbc,k,Tindex,itrc)
              OCEAN(ngc)%ad_t(Ibc-1,Jbc,k,Tindex,itrc)=0.0_r8
#  ifdef MASKING
!>            tl_Tvalue=tl_Tvalue*GRID(ngc)%rmask(Ibc-1,Jbc)
!>
              ad_Tvalue=ad_Tvalue*GRID(ngc)%rmask(Ibc-1,Jbc)
#  endif
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
!>              tl_Tvalue=tl_Tvalue-                                    &
!>   &                    dt(ngc)*                                      &
!>   &                    CLIMA(ngc)%Tnudgcof(Ibc-1,Jbc,k,itrc)*        &
!>   &                    tl_Tvalue
!>
                ad_Tvalue=ad_Tvalue*                                    &
     &                    (1.0_r8-dt(ngc)*                              &
     &                            CLIMA(ngc)%Tnudgcof(Ibc-1,Jbc,k,itrc))
              END IF
!>            tl_Tvalue=(0.5_r8-                                        &
!>   &                   SIGN(0.5_r8,                                   &
!>   &                        -(OCEAN(ngc)%t(Ibc-1,Jbc,k,Tindex,itrc)-  &
!>   &                        cff*(TFF-TFC))))*                         &
!>   &                   (OCEAN(ngc)%tl_t(Ibc-1,Jbc,k,Tindex,itrc)-     &
!>   &                   tl_cff*(TFF-TFC)-cff*(tl_TFF-tl_TFC))
!>
              adfac=(0.5_r8-                                            &
     &               SIGN(0.5_r8,                                       &
     &                    -(OCEAN(ngc)%t(Ibc-1,Jbc,k,Tindex,itrc)-      &
     &                    cff*(TFF-TFC))))*ad_Tvalue
              OCEAN(ngc)%ad_t(Ibc-1,Jbc,k,Tindex,itrc)=                 &
     &                   OCEAN(ngc)%ad_t(Ibc-1,Jbc,k,Tindex,itrc)+adfac
              ad_cff=ad_cff-(TFF-TFC)*adfac
              ad_TFF=ad_TFF-cff*adfac
              ad_TFC=ad_TFC+cff*adfac
              ad_Tvalue=0.0_r8
            END DO
!
!  Correct coarse grid tracer at the finer grid western boundary.
!
!>          tl_cff=GRID(ngc)%pm(Ibc-1,Jbc)*                             &
!>   &             GRID(ngc)%pn(Ibc-1,Jbc)*                             &
!>   &             tl_Dinv(Ibc-1,Jbc)
!>
            ad_Dinv(Ibc-1,Jbc)=ad_Dinv(Ibc-1,Jbc)+                      &
     &                         GRID(ngc)%pm(Ibc-1,Jbc)*                 &
     &                         GRID(ngc)%pn(Ibc-1,Jbc)*ad_cff
            ad_cff=0.0_r8
!
!  Zeroth order correction to fine grid time integral (RIL, 2016).
!
!>          tl_TFF=tl_TFF*dt(ngc)/dt(ngf)
!>
            ad_TFF=ad_TFF*dt(ngc)/dt(ngf)
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal J-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the J-edge.
!
            Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
            DO jsum=-half,half
              Jbf=Jedge+jsum
              DO k=1,N(ngf)
!>              tl_TFF=tl_TFF+                                          &
!>   &                 BRY_CONTACT(iwest,dgcr)%tl_Tflux(Jbf,k,itrc)
!>
                BRY_CONTACT(iwest,dgcr)%ad_Tflux(Jbf,k,itrc)=           &
     &              BRY_CONTACT(iwest,dgcr)%ad_Tflux(Jbf,k,itrc)+ad_TFF
              END DO
            END DO
!>          tl_TFF=0.0_r8
!>
            ad_TFF=0.0_r8
!
!  Sum vertically coarse grid horizontal advective tracer flux,
!  Hz*u*T/n, from last time-step.
!
            DO k=1,N(ngc)
!>            tl_TFC=tl_TFC+                                            &
!>   &               BRY_CONTACT(iwest,rgcr)%tl_Tflux(Jbc,k,itrc)
!>
              BRY_CONTACT(iwest,rgcr)%ad_Tflux(Jbc,k,itrc)=             &
     &            BRY_CONTACT(iwest,rgcr)%ad_Tflux(Jbc,k,itrc)+ad_TFC
            END DO
!>          tl_TFC=0.0_r8
!>
            ad_TFC=0.0_r8
          END IF
        END DO

      END DO T_LOOP
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
!>        tl_Dinv(i,j)=-tl_cff*Dinv(i,j)/cff
!>
          ad_cff=ad_cff-ad_Dinv(i,j)*Dinv(i,j)/cff
          ad_Dinv(i,j)=0.0_r8
          DO k=2,N(rg)
!>          tl_cff=tl_cff+GRID(ngc)%tl_Hz(i,j,k)
!>
            GRID(ngc)%ad_Hz(i,j,k)=GRID(ngc)%ad_Hz(i,j,k)+ad_cff
          END DO
!>        tl_cff=GRID(ngc)%tl_Hz(i,j,1)
!>
          GRID(ngc)%ad_Hz(i,j,1)=GRID(ngc)%ad_Hz(i,j,1)+ad_cff
          ad_cff=0.0_r8
        END DO
      END DO

      RETURN
      END SUBROUTINE ad_correct_tracer_tile
# endif
!
      SUBROUTINE ad_fine2coarse (ng, model, vtype, tile)
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
      USE ad_exchange_2d_mod
# ifdef SOLVE3D
      USE ad_exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d, ad_mp_exchange4d
#  endif
# endif
      USE strings_mod,     ONLY : FoundError
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
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_fine2coarse"
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
  10          FORMAT (6x,'AD_FINE2COARSE - exchanging data between ',   &
     &                'grids: dg = ',i2.2,' and rg = ',i2.2,            &
     &                '  at cr = ',i2.2)
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
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
          IF (EWperiodic(rg).or.NSperiodic(rg)) THEN
            IF (vtype.eq.r2dvar) THEN
# ifdef SOLVE3D
!>            CALL exchange_r2d_tile (rg, tile,                         &
!>   &                                LBiR, UBiR, LBjR, UBjR,           &
!>   &                                COUPLING(rg)%tl_Zt_avg1)
!>
              CALL ad_exchange_r2d_tile (rg, tile,                      &
     &                                   LBiR, UBiR, LBjR, UBjR,        &
     &                                   COUPLING(rg)%ad_Zt_avg1)
              DO k=1,2
!>              CALL exchange_u2d_tile (rg, tile,                       &
!>   &                                  LBiR, UBiR, LBjR, UBjR,         &
!>   &                                  OCEAN(rg)%tl_ubar(:,:,k))
!>
                CALL ad_exchange_u2d_tile (rg, tile,                    &
     &                                     LBiR, UBiR, LBjR, UBjR,      &
     &                                     OCEAN(rg)%ad_ubar(:,:,k))
!>              CALL exchange_v2d_tile (rg, tile,                       &
!>   &                                  LBiR, UBiR, LBjR, UBjR,         &
!>   &                                  OCEAN(rg)%tl_vbar(:,:,k))
!>
                CALL ad_exchange_v2d_tile (rg, tile,                    &
     &                                     LBiR, UBiR, LBjR, UBjR,      &
     &                                     OCEAN(rg)%ad_vbar(:,:,k))
              END DO
# else
!>            CALL exchange_r2d_tile (rg, tile,                         &
!>   &                                LBiR, UBiR, LBjR, UBjR,           &
!>   &                                OCEAN(rg)%tl_zeta(:,:,Rindex2d))
!>
              CALL ad_exchange_r2d_tile (rg, tile,                      &
     &                                   LBiR, UBiR, LBjR, UBjR,        &
     &                                   OCEAN(rg)%ad_zeta(:,:,         &
     &                                                     Rindex2d))
!>            CALL exchange_u2d_tile (rg, tile,                         &
!>   &                                LBiR, UBiR, LBjR, UBjR,           &
!>   &                                OCEAN(rg)%tl_ubar(:,:,Rindex2d))
!>
              CALL ad_exchange_u2d_tile (rg, tile,                      &
     &                                   LBiR, UBiR, LBjR, UBjR,        &
     &                                   OCEAN(rg)%ad_ubar(:,:,         &
     &                                                     Rindex2d))
!>            CALL exchange_v2d_tile (rg, tile,                         &
!>   &                                LBiR, UBiR, LBjR, UBjR,           &
!>   &                                OCEAN(rg)%tl_vbar(:,:,Rindex2d))
!>
              CALL ad_exchange_v2d_tile (rg, tile,                      &
     &                                   LBiR, UBiR, LBjR, UBjR,        &
     &                                   OCEAN(rg)%ad_vbar(:,:,         &
     &                                                     Rindex2d))
# endif
# ifdef SOLVE3D
            ELSE IF (vtype.eq.r3dvar) THEN
!>            CALL exchange_u3d_tile (rg, tile,                         &
!>   &                                LBiR, UBiR, LBjR, UBjR, 1, N(rg), &
!>   &                                OCEAN(rg)%tl_u(:,:,:,Rindex3d))
!>
              CALL ad_exchange_u3d_tile (rg, tile,                      &
     &                                   LBiR, UBiR, LBjR, UBjR,        &
     &                                   1, N(rg),                      &
     &                                   OCEAN(rg)%ad_u(:,:,:,Rindex3d))
!>            CALL exchange_v3d_tile (rg, tile,                         &
!>   &                                LBiR, UBiR, LBjR, UBjR, 1, N(rg), &
!>   &                                OCEAN(rg)%tl_v(:,:,:,Rindex3d))
!>
              CALL ad_exchange_v3d_tile (rg, tile,                      &
     &                                   LBiR, UBiR, LBjR, UBjR,        &
     &                                   1, N(rg),                      &
     &                                   OCEAN(rg)%ad_v(:,:,:,Rindex3d))
              DO itrc=1,NT(rg)
!>              CALL exchange_r3d_tile (rg, tile,                       &
!>   &                                  LBiR, UBiR, LBjR, UBjR,         &
!>   &                                  1, N(rg),                       &
!>   &                                  OCEAN(rg)%tl_t(:,:,:,Rindex3d,  &
!>   &                                                 itrc))
!>
                CALL ad_exchange_r3d_tile (rg, tile,                    &
     &                                     LBiR, UBiR, LBjR, UBjR,      &
     &                                     1, N(rg),                    &
     &                                     OCEAN(rg)%ad_t(:,:,:,        &
     &                                                    Rindex3d,     &
     &                                                    itrc))
              END DO
# endif
            END IF
          END IF

# ifdef DISTRIBUTE
!
          IF (vtype.eq.r2dvar) THEN
#  ifdef SOLVE3D
!>          CALL mp_exchange2d (rg, tile, model, 1,                     &
!>   &                          LBiR, UBiR, LBjR, UBjR,                 &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          COUPLING(rg)%tl_Zt_avg1)
!>
            CALL ad_mp_exchange2d (rg, tile, model, 1,                  &
     &                             LBiR, UBiR, LBjR, UBjR,              &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             COUPLING(rg)%ad_Zt_avg1)
!>          CALL mp_exchange2d (rg, tile, model, 4,                     &
!>   &                          LBiR, UBiR, LBjR, UBjR,                 &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          OCEAN(rg)%tl_ubar(:,:,1),               &
!>   &                          OCEAN(rg)%tl_vbar(:,:,1),               &
!>   &                          OCEAN(rg)%tl_ubar(:,:,2),               &
!>   &                          OCEAN(rg)%tl_vbar(:,:,2))
!>
            CALL ad_mp_exchange2d (rg, tile, model, 4,                  &
     &                             LBiR, UBiR, LBjR, UBjR,              &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             OCEAN(rg)%ad_ubar(:,:,1),            &
     &                             OCEAN(rg)%ad_vbar(:,:,1),            &
     &                             OCEAN(rg)%ad_ubar(:,:,2),            &
     &                             OCEAN(rg)%ad_vbar(:,:,2))
#  else
!>          CALL mp_exchange2d (rg, tile, model, 3,                     &
!>   &                          LBiR, UBiR, LBjR, UBjR,                 &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          OCEAN(rg)%tl_zeta(:,:,Rindex2d),        &
!>   &                          OCEAN(rg)%tl_ubar(:,:,Rindex2d),        &
!>   &                          OCEAN(rg)%tl_vbar(:,:,Rindex2d))
!>
            CALL ad_mp_exchange2d (rg, tile, model, 3,                  &
     &                             LBiR, UBiR, LBjR, UBjR,              &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             OCEAN(rg)%ad_zeta(:,:,Rindex2d),     &
     &                             OCEAN(rg)%ad_ubar(:,:,Rindex2d),     &
     &                             OCEAN(rg)%ad_vbar(:,:,Rindex2d))
#  endif
#  ifdef SOLVE3D
          ELSE IF (vtype.eq.r3dvar) THEN
!>          CALL mp_exchange3d (rg, tile, model, 2,                     &
!>   &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          OCEAN(rg)%tl_u(:,:,:,Rindex3d),         &
!>   &                          OCEAN(rg)%tl_v(:,:,:,Rindex3d))
!>
            CALL ad_mp_exchange3d (rg, tile, model, 2,                  &
     &                             LBiR, UBiR, LBjR, UBjR, 1, N(rg),    &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             OCEAN(rg)%ad_u(:,:,:,Rindex3d),      &
     &                             OCEAN(rg)%ad_v(:,:,:,Rindex3d))
!>          CALL mp_exchange4d (rg, tile, model, 1,                     &
!>   &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
!>   &                          1, NT(rg),                              &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(rg), NSperiodic(rg),         &
!>   &                          OCEAN(rg)%tl_t(:,:,:,Rindex3d,:))
!>
            CALL ad_mp_exchange4d (rg, tile, model, 1,                  &
     &                             LBiR, UBiR, LBjR, UBjR, 1, N(rg),    &
     &                             1, NT(rg),                           &
     &                             NghostPoints,                        &
     &                             EWperiodic(rg), NSperiodic(rg),      &
     &                             OCEAN(rg)%ad_t(:,:,:,Rindex3d,:))
#  endif
          END IF
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
!>          CALL fine2coarse2d (rg, dg, model, tile,                    &
!>   &                          r2dvar, 'Zt_avg1',                      &
!>   &                          AreaAvg, RefineScale(dg),               &
!>   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!>   &                          LBiD, UBiD, LBjD, UBjD,                 &
!>   &                          LBiR, UBiR, LBjR, UBjR,                 &
!>   &                          GRID(dg)%om_r,                          &
!>   &                          GRID(dg)%on_r,                          &
!>   &                          GRID(rg)%pm,                            &
!>   &                          GRID(rg)%pn,                            &
#  ifdef MASKING
!>   &                          GRID(dg)%rmask_full,                    &
!>   &                          GRID(rg)%rmask,                         &
#  endif
!>   &                          COUPLING(dg)%tl_Zt_avg1,                &
!>   &                          COUPLING(rg)%tl_Zt_avg1)
!>
            CALL ad_fine2coarse2d (rg, dg, model, tile,                 &
     &                             r2dvar, 'Zt_avg1',                   &
     &                             AreaAvg, RefineScale(dg),            &
     &                             cr, Rcontact(cr)%Npoints, Rcontact,  &
     &                             LBiD, UBiD, LBjD, UBjD,              &
     &                             LBiR, UBiR, LBjR, UBjR,              &
     &                             GRID(dg)%om_r,                       &
     &                             GRID(dg)%on_r,                       &
     &                             GRID(rg)%pm,                         &
     &                             GRID(rg)%pn,                         &
#  ifdef MASKING
     &                             GRID(dg)%rmask_full,                 &
     &                             GRID(rg)%rmask,                      &
#  endif
     &                             COUPLING(dg)%ad_Zt_avg1,             &
     &                             COUPLING(rg)%ad_Zt_avg1)
# else
!>          CALL fine2coarse2d (rg, dg, model, tile,                    &
!>   &                          r2dvar, Vname(1,idFsur),                &
!>   &                          AreaAvg, RefineScale(dg),               &
!>   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!>   &                          LBiD, UBiD, LBjD, UBjD,                 &
!>   &                          LBiR, UBiR, LBjR, UBjR,                 &
!>   &                          GRID(dg)%om_r,                          &
!>   &                          GRID(dg)%on_r,                          &
!>   &                          GRID(rg)%pm,                            &
!>   &                          GRID(rg)%pn,                            &
#  ifdef MASKING
!>   &                          GRID(dg)%rmask,                         &
!>   &                          GRID(rg)%rmask,                         &
#  endif
!>   &                          OCEAN(dg)%tl_zeta(:,:,Dindex2d),        &
!>   &                          OCEAN(rg)%tl_zeta(:,:,Rindex2d))
!>
            CALL ad_fine2coarse2d (rg, dg, model, tile,                 &
     &                             r2dvar, Vname(1,idFsur),             &
     &                             AreaAvg, RefineScale(dg),            &
     &                             cr, Rcontact(cr)%Npoints, Rcontact,  &
     &                             LBiD, UBiD, LBjD, UBjD,              &
     &                             LBiR, UBiR, LBjR, UBjR,              &
     &                             GRID(dg)%om_r,                       &
     &                             GRID(dg)%on_r,                       &
     &                             GRID(rg)%pm,                         &
     &                             GRID(rg)%pn,                         &
#  ifdef MASKING
     &                             GRID(dg)%rmask,                      &
     &                             GRID(rg)%rmask,                      &
#  endif
     &                             OCEAN(dg)%ad_zeta(:,:,Dindex2d),     &
     &                             OCEAN(rg)%ad_zeta(:,:,Rindex2d))
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Process 2D momentum components (ubar,vbar).
!
            AreaAvg=.FALSE.
!>          CALL fine2coarse2d (rg, dg, model, tile,                    &
!>   &                          u2dvar, Vname(1,idUbar),                &
!>   &                          AreaAvg, RefineScale(dg),               &
!>   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!>   &                          LBiD, UBiD, LBjD, UBjD,                 &
!>   &                          LBiR, UBiR, LBjR, UBjR,                 &
!>   &                          GRID(dg)%om_u,                          &
!>   &                          GRID(dg)%on_u,                          &
!>   &                          GRID(rg)%pm,                            &
!>   &                          GRID(rg)%pn,                            &
# ifdef MASKING
!>   &                          GRID(dg)%umask_full,                    &
!>   &                          GRID(rg)%umask_full,                    &
# endif
!>   &                          OCEAN(dg)%tl_ubar(:,:,Dindex2d),        &
# ifdef SOLVE3D
!>   &                          OCEAN(rg)%tl_ubar(:,:,1),               &
!>   &                          OCEAN(rg)%tl_ubar(:,:,2))
# else
!>   &                          OCEAN(rg)%tl_ubar(:,:,Rindex2d))
# endif
!>
            CALL ad_fine2coarse2d (rg, dg, model, tile,                 &
     &                             u2dvar, Vname(1,idUbar),             &
     &                             AreaAvg, RefineScale(dg),            &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBiD, UBiD, LBjD, UBjD,              &
     &                             LBiR, UBiR, LBjR, UBjR,              &
     &                             GRID(dg)%om_u,                       &
     &                             GRID(dg)%on_u,                       &
     &                             GRID(rg)%pm,                         &
     &                             GRID(rg)%pn,                         &
# ifdef MASKING
     &                             GRID(dg)%umask_full,                 &
     &                             GRID(rg)%umask_full,                 &
# endif
     &                             OCEAN(dg)%ad_ubar(:,:,Dindex2d),     &
# ifdef SOLVE3D
     &                             OCEAN(rg)%ad_ubar(:,:,1),            &
     &                             OCEAN(rg)%ad_ubar(:,:,2))
# else
     &                             OCEAN(rg)%ad_ubar(:,:,Rindex2d))
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!>          CALL fine2coarse2d (rg, dg, model, tile,                    &
!>   &                          v2dvar, Vname(1,idVbar),                &
!>   &                          AreaAvg, RefineScale(dg),               &
!>   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!>   &                          LBiD, UBiD, LBjD, UBjD,                 &
!>   &                          LBiR, UBiR, LBjR, UBjR,                 &
!>   &                          GRID(dg)%om_v,                          &
!>   &                          GRID(dg)%on_v,                          &
!>   &                          GRID(rg)%pm,                            &
!>   &                          GRID(rg)%pn,                            &
# ifdef MASKING
!>   &                          GRID(dg)%vmask_full,                    &
!>   &                          GRID(rg)%vmask_full,                    &
# endif
!>   &                          OCEAN(dg)%tl_vbar(:,:,Dindex2d),        &
# ifdef SOLVE3D
!>   &                          OCEAN(rg)%tl_vbar(:,:,1),               &
!>   &                          OCEAN(rg)%tl_vbar(:,:,2))
# else
!>   &                          OCEAN(rg)%tl_vbar(:,:,Rindex2d))
# endif
!>
            CALL ad_fine2coarse2d (rg, dg, model, tile,                 &
     &                             v2dvar, Vname(1,idVbar),             &
     &                             AreaAvg, RefineScale(dg),            &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBiD, UBiD, LBjD, UBjD,              &
     &                             LBiR, UBiR, LBjR, UBjR,              &
     &                             GRID(dg)%om_v,                       &
     &                             GRID(dg)%on_v,                       &
     &                             GRID(rg)%pm,                         &
     &                             GRID(rg)%pn,                         &
# ifdef MASKING
     &                             GRID(dg)%vmask_full,                 &
     &                             GRID(rg)%vmask_full,                 &
# endif
     &                             OCEAN(dg)%ad_vbar(:,:,Dindex2d),     &
# ifdef SOLVE3D
     &                             OCEAN(rg)%ad_vbar(:,:,1),            &
     &                             OCEAN(rg)%ad_vbar(:,:,2))
# else
     &                             OCEAN(rg)%ad_vbar(:,:,Rindex2d))
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
!>            CALL fine2coarse3d (rg, dg, model, tile,                  &
!>   &                            r3dvar, Vname(1,idTvar(itrc)),        &
!>   &                            AreaAvg, RefineScale(dg),             &
!>   &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
!>   &                            LBiD, UBiD, LBjD, UBjD, 1, N(dg),     &
!>   &                            LBiR, UBiR, LBjR, UBjR, 1, N(rg),     &
!>   &                            GRID(dg)%om_r,                        &
!>   &                            GRID(dg)%on_r,                        &
!>   &                            GRID(rg)%pm,                          &
!>   &                            GRID(rg)%pn,                          &
#  ifdef MASKING
!>   &                            GRID(dg)%rmask,                       &
!>   &                            GRID(rg)%rmask,                       &
#  endif
!>   &                            OCEAN(dg)%tl_t(:,:,:,Dindex3d,itrc),  &
!>   &                            OCEAN(rg)%tl_t(:,:,:,Rindex3d,itrc))
!>
              CALL ad_fine2coarse3d (rg, dg, model, tile,               &
     &                               r3dvar, Vname(1,idTvar(itrc)),     &
     &                               AreaAvg, RefineScale(dg),          &
     &                               cr, Rcontact(cr)%Npoints, Rcontact,&
     &                               LBiD, UBiD, LBjD, UBjD, 1, N(dg),  &
     &                               LBiR, UBiR, LBjR, UBjR, 1, N(rg),  &
     &                               GRID(dg)%om_r,                     &
     &                               GRID(dg)%on_r,                     &
     &                               GRID(rg)%pm,                       &
     &                               GRID(rg)%pn,                       &
#  ifdef MASKING
     &                               GRID(dg)%rmask,                    &
     &                               GRID(rg)%rmask,                    &
#  endif
     &                               OCEAN(dg)%ad_t(:,:,:,Dindex3d,     &
     &                                              itrc),              &
     &                               OCEAN(rg)%ad_t(:,:,:,Rindex3d,     &
     &                                              itrc))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END DO
!
!  Process 3D momentum components (u, v).
!
            AreaAvg=.FALSE.
!>          CALL fine2coarse3d (rg, dg, model, tile,                    &
!>   &                          u3dvar, Vname(1,idUvel),                &
!>   &                          AreaAvg, RefineScale(dg),               &
!>   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!>   &                          LBiD, UBiD, LBjD, UBjD, 1, N(dg),       &
!>   &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
!>   &                          GRID(dg)%om_u,                          &
!>   &                          GRID(dg)%on_u,                          &
!>   &                          GRID(rg)%pm,                            &
!>   &                          GRID(rg)%pn,                            &
#  ifdef MASKING
!>   &                          GRID(dg)%umask_full,                    &
!>   &                          GRID(rg)%umask_full,                    &
#  endif
!>   &                          OCEAN(dg)%tl_u(:,:,:,Dindex3d),         &
!>   &                          OCEAN(rg)%tl_u(:,:,:,Rindex3d))
!>
            CALL ad_fine2coarse3d (rg, dg, model, tile,                 &
     &                             u3dvar, Vname(1,idUvel),             &
     &                             AreaAvg, RefineScale(dg),            &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBiD, UBiD, LBjD, UBjD, 1, N(dg),    &
     &                             LBiR, UBiR, LBjR, UBjR, 1, N(rg),    &
     &                             GRID(dg)%om_u,                       &
     &                             GRID(dg)%on_u,                       &
     &                             GRID(rg)%pm,                         &
     &                             GRID(rg)%pn,                         &
#  ifdef MASKING
     &                             GRID(dg)%umask_full,                 &
     &                             GRID(rg)%umask_full,                 &
#  endif
     &                             OCEAN(dg)%ad_u(:,:,:,Dindex3d),      &
     &                             OCEAN(rg)%ad_u(:,:,:,Rindex3d))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!>          CALL fine2coarse3d (rg, dg, model, tile,                    &
!>   &                          v3dvar, Vname(1,idVvel),                &
!>   &                          AreaAvg, RefineScale(dg),               &
!>   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!>   &                          LBiD, UBiD, LBjD, UBjD, 1, N(dg),       &
!>   &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
!>   &                          GRID(dg)%om_v,                          &
!>   &                          GRID(dg)%on_v,                          &
!>   &                          GRID(rg)%pm,                            &
!>   &                          GRID(rg)%pn,                            &
#  ifdef MASKING
!>   &                          GRID(dg)%vmask_full,                    &
!>   &                          GRID(rg)%vmask_full,                    &
#  endif
!>   &                          OCEAN(dg)%tl_v(:,:,:,Dindex3d),         &
!>   &                          OCEAN(rg)%tl_v(:,:,:,Rindex3d))
!>
            CALL ad_fine2coarse3d (rg, dg, model, tile,                 &
     &                             v3dvar, Vname(1,idVvel),             &
     &                             AreaAvg, RefineScale(dg),            &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBiD, UBiD, LBjD, UBjD, 1, N(dg),    &
     &                             LBiR, UBiR, LBjR, UBjR, 1, N(rg),    &
     &                             GRID(dg)%om_v,                       &
     &                             GRID(dg)%on_v,                       &
     &                             GRID(rg)%pm,                         &
     &                             GRID(rg)%pn,                         &
#  ifdef MASKING
     &                             GRID(dg)%vmask_full,                 &
     &                             GRID(rg)%vmask_full,                 &
#  endif
     &                             OCEAN(dg)%ad_v(:,:,:,Dindex3d),      &
     &                             OCEAN(rg)%ad_v(:,:,:,Rindex3d))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
          END IF
!
        END IF
      END DO

      RETURN
      END SUBROUTINE ad_fine2coarse
!
      SUBROUTINE ad_put_refine2d (ng, dg, cr, model, tile, LputFsur,    &
     &                            LBi, UBi, LBj, UBj)
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
      USE distribute_mod,  ONLY : mp_assemble
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
      USE strings_mod,     ONLY : FoundError
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
      integer :: NSUB, i, irec, j, m, tnew, told, ii
      integer :: Idg, Jdg

# ifdef DISTRIBUTE
      real(r8), parameter :: spv = 0.0_r8
# endif
      real(dp) :: Wnew, Wold, SecScale, fac
      real(r8) :: cff, cff1, my_value
      real(r8) :: ad_cff, adfac, adfac1, adfac2, ad_my_value
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_put_refine2d"

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Interpolate (space, time) refinement grid contact points for 2D state
!  variables from donor grid.
!-----------------------------------------------------------------------
!
!  Clear adjoint constants.
!
      ad_cff=0.0_r8
      adfac=0.0_r8
      adfac1=0.0_r8
      adfac2=0.0_r8
      ad_my_value=0.0_r8

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

#  ifdef NESTING_DEBUG
!
!  If distributed-memory, initialize arrays used to check mass flux
!  conservation with special value (zero) to facilitate the global
!  reduction when collecting data between all nodes.
!
        BRY_CONTACT(iwest ,cr)%Mflux=spv
        BRY_CONTACT(ieast ,cr)%Mflux=spv
        BRY_CONTACT(isouth,cr)%Mflux=spv
        BRY_CONTACT(inorth,cr)%Mflux=spv
#  endif
      END IF
# endif
!
!  Set time snapshot indices for the donor grid data.
!
      told=3-RollingIndex(cr)
      tnew=RollingIndex(cr)
!
!  Set linear time interpolation weights. Fractional seconds are
!  rounded to the nearest milliseconds integer towards zero in the
!  time interpolation weights.
!
      SecScale=1000.0_dp              ! seconds to milliseconds
!
      Wold=ANINT((RollingTime(tnew,cr)-time(ng))*SecScale,dp)
      Wnew=ANINT((time(ng)-RollingTime(told,cr))*SecScale,dp)
      fac=1.0_dp/(Wold+Wnew)
      Wold=fac*Wold
      Wnew=fac*Wnew
!
!     IF (((Wold*Wnew).lt.0.0_dp).or.((Wold+Wnew).le.0.0_dp)) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) cr, dg, ng,                               &
     &                        iic(dg), told, tnew,                      &
     &                        iic(ng), Wold, Wnew,                      &
     &                        INT(time(ng)),                            &
     &                        INT(RollingTime(told,cr)),                &
     &                        INT(RollingTime(tnew,cr))
          END IF
!         exit_flag=8
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!     END IF

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile information.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      IF (LputFsur) THEN
!>      CALL mp_exchange2d (ng, tile, model,                            &
#  ifdef SOLVE3D
!>   &                      4,                                          &
#  else
!>   &                      3,                                          &
#  endif
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
#  ifdef SOLVE3D
!>   &                      COUPLING(ng)%tl_Zt_avg1,                    &
#  endif
!>   &                      OCEAN(ng)%tl_zeta(:,:,1),                   &
!>   &                      OCEAN(ng)%tl_zeta(:,:,2),                   &
!>   &                      OCEAN(ng)%tl_zeta(:,:,3))
!>
        CALL ad_mp_exchange2d (ng, tile, model,                         &
#  ifdef SOLVE3D
     &                         4,                                       &
#  else
     &                         3,                                       &
#  endif
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
#  ifdef SOLVE3D
     &                         COUPLING(ng)%ad_Zt_avg1,                 &
#  endif
     &                         OCEAN(ng)%ad_zeta(:,:,1),                &
     &                         OCEAN(ng)%ad_zeta(:,:,2),                &
     &                         OCEAN(ng)%ad_zeta(:,:,3))
!
!  2D momentum.
!
      ELSE

#  ifdef NESTING_DEBUG
!
!  No action required for the adjoint of mp_assemble (AMM).
!
!>      CALL mp_assemble (ng, model, NptsSN, spv,                       &
!>   &                    BRY_CONTACT(inorth,cr)%tl_Mflux(ILB:))
!>
!!      CALL ad_mp_assemble (ng, model, NptsSN, spv,                    &
!!   &                       BRY_CONTACT(inorth,cr)%ad_Mflux(ILB:))
!>      CALL mp_assemble (ng, model, NptsSN, spv,                       &
!>   &                    BRY_CONTACT(isouth,cr)%tl_Mflux(ILB:))
!>
!!      CALL ad_mp_assemble (ng, model, NptsSN, spv,                    &
!!   &                       BRY_CONTACT(isouth,cr)%ad_Mflux(ILB:))
!>      CALL mp_assemble (ng, model, NptsWE, spv,                       &
!>   &                    BRY_CONTACT(ieast ,cr)%tl_Mflux(JLB:))
!!      CALL ad_mp_assemble (ng, model, NptsWE, spv,                    &
!!   &                       BRY_CONTACT(ieast ,cr)%ad_Mflux(JLB:))
!>      CALL mp_assemble (ng, model, NptsWE, spv,                       &
!>   &                    BRY_CONTACT(iwest ,cr)%tl_Mflux(JLB:))
!>
!!      CALL ad_mp_assemble (ng, model, NptsWE, spv,                    &
!!   &                       BRY_CONTACT(iwest ,cr)%ad_Mflux(JLB:))
#  endif

!>      CALL mp_exchange2d (ng, tile, model, 3,                         &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      OCEAN(ng)%tl_vbar(:,:,1),                   &
!>   &                      OCEAN(ng)%tl_vbar(:,:,2),                   &
!>   &                      OCEAN(ng)%tl_vbar(:,:,3))
!>
        CALL ad_mp_exchange2d (ng, tile, model, 3,                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         OCEAN(ng)%ad_vbar(:,:,1),                &
     &                         OCEAN(ng)%ad_vbar(:,:,2),                &
     &                         OCEAN(ng)%ad_vbar(:,:,3))

!>      CALL mp_exchange2d (ng, tile, model, 3,                         &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      OCEAN(ng)%tl_ubar(:,:,1),                   &
!>   &                      OCEAN(ng)%tl_ubar(:,:,2),                   &
!>   &                      OCEAN(ng)%tl_ubar(:,:,3))
!>
        CALL ad_mp_exchange2d (ng, tile, model, 3,                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         OCEAN(ng)%ad_ubar(:,:,1),                &
     &                         OCEAN(ng)%ad_ubar(:,:,2),                &
     &                         OCEAN(ng)%ad_ubar(:,:,3))

      END IF
# endif

      FREE_SURFACE : IF (LputFsur) THEN
        DO m=1,Rcontact(cr)%Npoints
          i=Rcontact(cr)%Irg(m)
          j=Rcontact(cr)%Jrg(m)
          IF (((IstrT.le.i).and.(i.le.IendT)).and.                      &
     &        ((JstrT.le.j).and.(j.le.JendT))) THEN
# ifdef SOLVE3D
!>          COUPLING(ng)%tl_Zt_avg1(i,j)=tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  COUPLING(ng)%ad_Zt_avg1(i,j)
            COUPLING(ng)%ad_Zt_avg1(i,j)=0.0_r8
# endif
!>          OCEAN(ng)%tl_zeta(i,j,1)=tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  OCEAN(ng)%ad_zeta(i,j,1)
            OCEAN(ng)%ad_zeta(i,j,1)=0.0_r8
!>          OCEAN(ng)%tl_zeta(i,j,2)=tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  OCEAN(ng)%ad_zeta(i,j,2)
            OCEAN(ng)%ad_zeta(i,j,2)=0.0_r8
!>          OCEAN(ng)%tl_zeta(i,j,3)=tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  OCEAN(ng)%ad_zeta(i,j,3)
            OCEAN(ng)%ad_zeta(i,j,3)=0.0_r8
# ifdef WET_DRY
            IF (my_value.le.(Dcrit(ng)-GRID(ng)%h(i,j))) THEN
!>            tl_my_value=-GRID(ng)%tl_h(i,j)
!>
              GRID(ng)%ad_h(i,j)=GRID(ng)%ad_h(i,j)-ad_my_value
              ad_my_value=0.0_r8
            END IF
# endif
# ifdef MASKING
!>          tl_my_value=tl_my_value*GRID(ng)%rmask(i,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%rmask(i,j)
# endif
!>          tl_my_value=Wold*                                           &
!>   &                  (Rcontact(cr)%Lweight(1,m)*                     &
!>   &                            REFINED(cr)%tl_zeta(1,m,told)+        &
!>   &                   Rcontact(cr)%Lweight(2,m)*                     &
!>   &                            REFINED(cr)%tl_zeta(2,m,told)+        &
!>   &                   Rcontact(cr)%Lweight(3,m)*                     &
!>   &                            REFINED(cr)%tl_zeta(3,m,told)+        &
!>   &                   Rcontact(cr)%Lweight(4,m)*                     &
!>   &                            REFINED(cr)%tl_zeta(4,m,told))+       &
!>   &                  Wnew*                                           &
!>                       (Rcontact(cr)%Lweight(1,m)*                    &
!>   &                             REFINED(cr)%tl_zeta(1,m,tnew)+       &
!>   &                    Rcontact(cr)%Lweight(2,m)*                    &
!>   &                             REFINED(cr)%tl_zeta(2,m,tnew)+       &
!>   &                    Rcontact(cr)%Lweight(3,m)*                    &
!>   &                             REFINED(cr)%tl_zeta(3,m,tnew)+       &
!>   &                    Rcontact(cr)%Lweight(4,m)*                    &
!>   &                             REFINED(cr)%tl_zeta(4,m,tnew))
!>
            DO ii=1,4
              adfac1=Wold*Rcontact(cr)%Lweight(ii,m)*ad_my_value
              adfac2=Wnew*Rcontact(cr)%Lweight(ii,m)*ad_my_value
              REFINED(cr)%ad_zeta(ii,m,told)=                           &
     &                REFINED(cr)%ad_zeta(ii,m,told)+adfac1
              REFINED(cr)%ad_zeta(ii,m,tnew)=                           &
     &                REFINED(cr)%ad_zeta(ii,m,tnew)+adfac2
            END DO
            ad_my_value=0.0_r8
          END IF
        END DO
      ELSE

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
#  ifdef MASKING
            my_value=my_value*GRID(ng)%vmask(i,Jend+1)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%vmask_wet(i,Jend+1)
#  endif
!>          OCEAN(ng)%tl_vbar(i,Jend+1,indx1(ng))=tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  OCEAN(ng)%ad_vbar(i,Jend+1,indx1(ng))
            OCEAN(ng)%ad_vbar(i,Jend+1,indx1(ng))=0.0_r8
#  ifdef NESTING_DEBUG
!>          BRY_CONTACT(inorth,cr)%tl_Mflux(i)=tl_cff*my_value+         &
!>   &                                         cff*tl_my_value
!>
            ad_cff=ad_cff+                                              &
     &             my_value*BRY_CONTACT(inorth,cr)%ad_Mflux(i)
            ad_my_value=ad_my_value+                                    &
     &                  cff*BRY_CONTACT(inorth,cr)%ad_Mflux(i)
            BRY_CONTACT(inorth,cr)%ad_Mflux(i)=0.0_r8
#  endif
#  ifdef WET_DRY
!>          tl_my_value=tl_my_value*GRID(ng)%vmask_wet(i,Jend+1)
!>
            ad_my_value=ad_my_value*GRID(ng)%vmask_wet(i,Jend+1)
#  endif
#  ifdef MASKING
!>          tl_my_value=tl_my_value*GRID(ng)%vmask(i,Jend+1)
!>
            ad_my_value=ad_my_value*GRID(ng)%vmask(i,Jend+1)
#  endif
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DV_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DV_avg2(1,m,tnew))/cff
!>          tl_my_value=cff1*                                           &
!>   &                  (Wold*REFINED(cr)%tl_DV_avg2(1,m,told)+         &
!>   &                   Wnew*REFINED(cr)%tl_DV_avg2(1,m,tnew))/cff-    &
!>   &                  tl_cff*my_value/cff
!>
            adfac=ad_my_value/cff
            adfac1=cff1*adfac
            REFINED(cr)%ad_DV_avg2(1,m,told)=                           &
     &             REFINED(cr)%ad_DV_avg2(1,m,told)+Wold*adfac1
            REFINED(cr)%ad_DV_avg2(1,m,tnew)=                           &
     &             REFINED(cr)%ad_DV_avg2(1,m,tnew)+Wnew*adfac1
            ad_cff=ad_cff-                                              &
     &             my_value*adfac
            ad_my_value=0.0_r8
#  else
            my_value=cff1*REFINED(cr)%DV_avg2(1,m,tnew)/cff
!>          tl_my_value=cff1*REFINED(cr)%tl_DV_avg2(1,m,tnew)/cff-      &
!>   &                  tl_cff*my_value/cff
!>
            adfac=ad_my_value/cff
            REFINED(cr)%ad_DV_avg2(1,m,tnew)=                           &
     &             REFINED(cr)%ad_DV_avg2(1,m,tnew)+cff1*adfac
            ad_cff=ad_cff-                                              &
     &             my_value*adfac
            ad_my_value=0.0_r8
#  endif
!>          tl_cff=0.5_r8*GRID(ng)%om_v(i,Jend+1)*                      &
!>   &             (GRID(ng)%tl_h(i,Jend+1)+                            &
!>   &              OCEAN(ng)%tl_zeta(i,Jend+1,indx1(ng))+              &
!>   &              GRID(ng)%tl_h(i,Jend  )+                            &
!>   &              OCEAN(ng)%tl_zeta(i,Jend  ,indx1(ng)))
!>
            adfac=0.5_r8*GRID(ng)%om_v(i,Jend+1)*ad_cff
            GRID(ng)%ad_h(i,Jend  )=GRID(ng)%ad_h(i,Jend  )+adfac
            GRID(ng)%ad_h(i,Jend+1)=GRID(ng)%ad_h(i,Jend+1)+adfac
            OCEAN(ng)%ad_zeta(i,Jend  ,indx1(ng))=                      &
     &                OCEAN(ng)%ad_zeta(i,Jend  ,indx1(ng))+adfac
            OCEAN(ng)%ad_zeta(i,Jend+1,indx1(ng))=                      &
     &                OCEAN(ng)%ad_zeta(i,Jend+1,indx1(ng))+adfac
            ad_cff=0.0_r8
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
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DV_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DV_avg2(1,m,tnew))/cff
#  else
            my_value=cff1*REFINED(cr)%DV_avg2(1,m,tnew)/cff
#  endif
#  ifdef MASKING
            my_value=my_value*GRID(ng)%vmask(i,Jstr)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%vmask_wet(i,Jstr)
#  endif
!>          OCEAN(ng)%tl_vbar(i,Jstr,indx1(ng))=tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  OCEAN(ng)%ad_vbar(i,Jstr,indx1(ng))
            OCEAN(ng)%ad_vbar(i,Jstr,indx1(ng))=0.0_r8
#  ifdef NESTING_DEBUG
!>          BRY_CONTACT(isouth,cr)%tl_Mflux(i)=tl_cff*my_value+         &
!>   &                                         cff*tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  cff*BRY_CONTACT(isouth,cr)%ad_Mflux(i)
            ad_cff=ad_cff+                                              &
     &             my_value*BRY_CONTACT(isouth,cr)%ad_Mflux(i)
            BRY_CONTACT(isouth,cr)%ad_Mflux(i)=0.0_r8
#  endif
#  ifdef WET_DRY
!>          tl_my_value=tl_my_value*GRID(ng)%vmask_wet(i,Jstr)
!>
            ad_my_value=ad_my_value*GRID(ng)%vmask_wet(i,Jstr)
#  endif
#  ifdef MASKING
!>          tl_my_value=tl_my_value*GRID(ng)%vmask(i,Jstr)
!>
            ad_my_value=ad_my_value*GRID(ng)%vmask(i,Jstr)
#  endif
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DV_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DV_avg2(1,m,tnew))/cff
!>          tl_my_value=cff1*                                           &
!>   &                  (Wold*REFINED(cr)%tl_DV_avg2(1,m,told)+         &
!>   &                   Wnew*REFINED(cr)%tl_DV_avg2(1,m,tnew))/cff-    &
!>   &                  tl_cff*my_value/cff
!>
            adfac=ad_my_value/cff
            adfac1=cff1*adfac
            REFINED(cr)%ad_DV_avg2(1,m,told)=                           &
     &             REFINED(cr)%ad_DV_avg2(1,m,told)+Wold*adfac1
            REFINED(cr)%tl_DV_avg2(1,m,tnew)=                           &
     &             REFINED(cr)%ad_DV_avg2(1,m,tnew)+Wnew*adfac1
            ad_cff=ad_cff-                                              &
     &             my_value*adfac
            ad_my_value=0.0_r8
#  else
            my_value=cff1*REFINED(cr)%DV_avg2(1,m,tnew)/cff
!>          tl_my_value=cff1*REFINED(cr)%tl_DV_avg2(1,m,tnew)/cff-      &
!>   &                  tl_cff*my_value/cff
!>
            adfac=ad_my_value/cff
            REFINED(cr)%ad_DV_avg2(1,m,tnew)=                           &
     &             REFINED(cr)%ad_DV_avg2(1,m,tnew)+cff1*adfac
            ad_cff=ad_cff-                                              &
     &             my_value*adfac
            ad_my_value=0.0_r8
#  endif
!>          tl_cff=0.5_r8*GRID(ng)%om_v(i,Jstr)*                        &
!>   &             (GRID(ng)%tl_h(i,Jstr-1)+                            &
!>   &              OCEAN(ng)%tl_zeta(i,Jstr-1,indx1(ng))+              &
!>   &              GRID(ng)%tl_h(i,Jstr  )+                            &
!>   &              OCEAN(ng)%tl_zeta(i,Jstr  ,indx1(ng)))
!>
            adfac=0.5_r8*GRID(ng)%om_v(i,Jstr)*ad_cff
            GRID(ng)%ad_h(i,Jstr-1)=GRID(ng)%ad_h(i,Jstr-1)+adfac
            GRID(ng)%ad_h(i,Jstr  )=GRID(ng)%ad_h(i,Jstr  )+adfac
            OCEAN(ng)%ad_zeta(i,Jstr-1,indx1(ng))=                      &
     &                OCEAN(ng)%ad_zeta(i,Jstr-1,indx1(ng))+adfac
            OCEAN(ng)%ad_zeta(i,Jstr  ,indx1(ng))=                      &
     &                OCEAN(ng)%ad_zeta(i,Jstr  ,indx1(ng))+adfac
            ad_cff=0.0_r8
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
#  ifdef MASKING
            my_value=my_value*GRID(ng)%umask(Iend+1,j)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%umask_wet(Iend+1,j)
#  endif
!>          OCEAN(ng)%tl_ubar(Iend+1,j,indx1(ng))=tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  OCEAN(ng)%ad_ubar(Iend+1,j,indx1(ng))
            OCEAN(ng)%ad_ubar(Iend+1,j,indx1(ng))=0.0_r8
#  ifdef NESTING_DEBUG
!>          BRY_CONTACT(ieast,cr)%tl_Mflux(j)=tl_cff*my_value+          &
!>   &                                        cff*tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  cff*BRY_CONTACT(ieast,cr)%ad_Mflux(j)
            ad_cff=ad_cff+                                              &
     &             my_value*BRY_CONTACT(ieast,cr)%ad_Mflux(j)
            BRY_CONTACT(ieast,cr)%ad_Mflux(j)=0.0_r8
#  endif
#  ifdef WET_DRY
!>          tl_my_value=tl_my_value*GRID(ng)%umask_wet(Iend+1,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%umask_wet(Iend+1,j)
#  endif
#  ifdef MASKING
!>          tl_my_value=tl_my_value*GRID(ng)%umask(Iend+1,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%umask(Iend+1,j)
#  endif
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DU_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DU_avg2(1,m,tnew))/cff
!>          tl_my_value=cff1*                                           &
!>   &                  (Wold*REFINED(cr)%tl_DU_avg2(1,m,told)+         &
!>   &                   Wnew*REFINED(cr)%tl_DU_avg2(1,m,tnew))/cff-    &
!>   &                  tl_cff*my_value/cff
!>
            adfac=ad_my_value/cff
            adfac1=cff1*adfac
            REFINED(cr)%ad_DU_avg2(1,m,told)=                           &
     &              REFINED(cr)%ad_DU_avg2(1,m,told)+Wold*adfac1
            REFINED(cr)%ad_DU_avg2(1,m,tnew)=                           &
     &              REFINED(cr)%ad_DU_avg2(1,m,tnew)+Wnew*adfac1
            ad_cff=ad_cff-                                              &
     &             my_value*adfac
            ad_my_value=0.0_r8
#  else
            my_value=cff1*REFINED(cr)%DU_avg2(1,m,tnew)/cff
!>          tl_my_value=cff1*REFINED(cr)%tl_DU_avg2(1,m,tnew)/cff-      &
!>   &                  tl_cff*my_value/cff
!>
            adfac=ad_my_value/cff
            REFINED(cr)%ad_DU_avg2(1,m,tnew)=                           &
     &              REFINED(cr)%ad_DU_avg2(1,m,tnew)+cff1*adfac
            ad_cff=ad_cff-                                              &
     &             my_value*adfac
            ad_my_value=0.0_r8
#  endif
!>          tl_cff=0.5_r8*GRID(ng)%on_u(Iend+1,j)*                      &
!>   &             (GRID(ng)%tl_h(Iend+1,j)+                            &
!>   &              OCEAN(ng)%tl_zeta(Iend+1,j,indx1(ng))+              &
!>   &              GRID(ng)%tl_h(Iend  ,j)+                            &
!>   &              OCEAN(ng)%tl_zeta(Iend  ,j,indx1(ng)))
!>
            adfac=0.5_r8*GRID(ng)%on_u(Iend+1,j)*ad_cff
            GRID(ng)%ad_h(Iend  ,j)=GRID(ng)%ad_h(Iend  ,j)+adfac
            GRID(ng)%ad_h(Iend+1,j)=GRID(ng)%ad_h(Iend+1,j)+adfac
            OCEAN(ng)%ad_zeta(Iend  ,j,indx1(ng))=                      &
     &                OCEAN(ng)%ad_zeta(Iend  ,j,indx1(ng))+adfac
            OCEAN(ng)%ad_zeta(Iend+1,j,indx1(ng))=                      &
     &                OCEAN(ng)%ad_zeta(Iend+1,j,indx1(ng))+adfac
            ad_cff=0.0_r8
          END DO
        END IF
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
#  ifdef MASKING
            my_value=my_value*GRID(ng)%umask(Istr,j)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%umask_wet(Istr,j)
#  endif
!>          OCEAN(ng)%tl_ubar(Istr,j,indx1(ng))=tl_my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  OCEAN(ng)%ad_ubar(Istr,j,indx1(ng))
            OCEAN(ng)%ad_ubar(Istr,j,indx1(ng))=0.0_r8
#  ifdef NESTING_DEBUG
!>          BRY_CONTACT(iwest,cr)%tl_Mflux(j)=cff*tl_my_value+          &
!>   &                                        tl_cff*my_value
!>
            ad_my_value=ad_my_value+                                    &
     &                  cff*BRY_CONTACT(iwest,cr)%ad_Mflux(j)
            ad_cff=ad_cff+                                              &
     &             my_value*BRY_CONTACT(iwest,cr)%ad_Mflux(j)
            BRY_CONTACT(iwest,cr)%ad_Mflux(j)=0.0_r8
#  endif
#  ifdef WET_DRY
!>          tl_my_value=tl_my_value*GRID(ng)%umask_wet(Istr,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%umask_wet(Istr,j)
#  endif
#  ifdef MASKING
!>          tl_my_value=tl_my_value*GRID(ng)%umask(Istr,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%umask(Istr,j)
#  endif
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DU_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DU_avg2(1,m,tnew))/cff
!>          tl_my_value=cff1*                                           &
!>   &                  (Wold*REFINED(cr)%tl_DU_avg2(1,m,told)+         &
!>   &                   Wnew*REFINED(cr)%tl_DU_avg2(1,m,tnew))/cff-    &
!>   &                  tl_cff*my_value/cff
!>
            adfac=ad_my_value/cff
            adfac1=cff1*adfac
            REFINED(cr)%ad_DU_avg2(1,m,told)=                           &
     &              REFINED(cr)%ad_DU_avg2(1,m,told)+Wold*adfac1
            REFINED(cr)%ad_DU_avg2(1,m,tnew)=                           &
     &              REFINED(cr)%ad_DU_avg2(1,m,tnew)+Wnew*adfac1
            ad_cff=ad_cff-my_value*adfac
            ad_my_value=0.0_r8
#  else
            my_value=cff1*REFINED(cr)%DU_avg2(1,m,tnew)/cff
!>          tl_my_value=cff1*REFINED(cr)%tl_DU_avg2(1,m,tnew)/cff-      &
!>   &                  tl_cff*my_value/cff
!>
            adfac=ad_my_value/cff
            REFINED(cr)%ad_DU_avg2(1,m,tnew)=                           &
     &              REFINED(cr)%ad_DU_avg2(1,m,tnew)+cff1*adfac
            ad_cff=ad_cff-                                              &
     &             my_value*adfac
            ad_my_value=0.0_r8
#  endif
!>          tl_cff=0.5_r8*GRID(ng)%on_u(Istr,j)*                        &
!>   &             (GRID(ng)%tl_h(Istr-1,j)+                            &
!>   &              OCEAN(ng)%tl_zeta(Istr-1,j,indx1(ng))+              &
!>   &              GRID(ng)%tl_h(Istr  ,j)+                            &
!>   &              OCEAN(ng)%tl_zeta(Istr  ,j,indx1(ng)))
!>
            adfac=0.5_r8*GRID(ng)%on_u(Istr,j)*ad_cff
            GRID(ng)%ad_h(Istr-1,j)=GRID(ng)%ad_h(Istr-1,j)+adfac
            GRID(ng)%ad_h(Istr  ,j)=GRID(ng)%ad_h(Istr  ,j)+adfac
            OCEAN(ng)%ad_zeta(Istr-1,j,indx1(ng))=                      &
     &                OCEAN(ng)%ad_zeta(Istr-1,j,indx1(ng))+adfac
            OCEAN(ng)%ad_zeta(Istr  ,j,indx1(ng))=                      &
     &                OCEAN(ng)%ad_zeta(Istr  ,j,indx1(ng))+adfac
            ad_cff=0.0_r8
          END DO
        END IF
# endif
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
            DO irec=1,3
# ifdef SOLVE3D
              Vboundary=(m.eq.BRY_CONTACT(isouth,cr)%C2Bindex(i)).or.   &
     &                  (m.eq.BRY_CONTACT(inorth,cr)%C2Bindex(i))
              IF(.not.(Vboundary.and.(irec.eq.indx1(ng)))) THEN
!>              OCEAN(ng)%tl_vbar(i,j,irec)=tl_my_value
!>
                ad_my_value=ad_my_value+                                &
     &                      OCEAN(ng)%ad_vbar(i,j,irec)
                OCEAN(ng)%ad_vbar(i,j,irec)=0.0_r8
!!            ELSE                                   ! for debugging
!!              OCEAN(ng)%vbar(i,j,irec)=0.0_r8      ! purposes
              END IF
# else
!>            OCEAN(ng)%tl_vbar(i,j,irec)=tl_my_value
!>
              ad_my_value=ad_my_value+                                  &
     &                    OCEAN(ng)%ad_vbar(i,j,irec)
              OCEAN(ng)%ad_vbar(i,j,irec)=0.0_r8
# endif
            END DO
# ifdef WET_DRY
!>          tl_my_value=tl_my_value*GRID(ng)%vmask_wet(i,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%vmask_wet(i,j)
# endif
# ifdef MASKING
!>          tl_my_value=tl_my_value*GRID(ng)%vmask(i,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%vmask(i,j)
# endif
!>          tl_my_value=Wold*                                           &
!>   &                  (Vcontact(cr)%Lweight(1,m)*                     &
!>   &                            REFINED(cr)%tl_vbar(1,m,told)+        &
!>   &                   Vcontact(cr)%Lweight(2,m)*                     &
!>   &                            REFINED(cr)%tl_vbar(2,m,told)+        &
!>   &                   Vcontact(cr)%Lweight(3,m)*                     &
!>   &                            REFINED(cr)%tl_vbar(3,m,told)+        &
!>   &                   Vcontact(cr)%Lweight(4,m)*                     &
!>   &                            REFINED(cr)%tl_vbar(4,m,told))+       &
!>   &                  Wnew*                                           &
!>   &                  (Vcontact(cr)%Lweight(1,m)*                     &
!>   &                            REFINED(cr)%tl_vbar(1,m,tnew)+        &
!>   &                   Vcontact(cr)%Lweight(2,m)*                     &
!>   &                            REFINED(cr)%tl_vbar(2,m,tnew)+        &
!>   &                   Vcontact(cr)%Lweight(3,m)*                     &
!>   &                            REFINED(cr)%tl_vbar(3,m,tnew)+        &
!>   &                   Vcontact(cr)%Lweight(4,m)*                     &
!>   &                            REFINED(cr)%tl_vbar(4,m,tnew))
!>
            DO ii=1,4
              adfac1=Wold*Vcontact(cr)%Lweight(ii,m)*ad_my_value
              adfac2=Wnew*Vcontact(cr)%Lweight(ii,m)*ad_my_value
              REFINED(cr)%ad_vbar(ii,m,told)=                           &
     &                REFINED(cr)%ad_vbar(ii,m,told)+adfac1
              REFINED(cr)%ad_vbar(ii,m,tnew)=                           &
     &                REFINED(cr)%ad_vbar(ii,m,tnew)+adfac2
            END DO
            ad_my_value=0.0_r8
          END IF
        END DO
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
            DO irec=1,3
# ifdef SOLVE3D
              Uboundary=(m.eq.BRY_CONTACT(iwest,cr)%C2Bindex(j)).or.    &
     &                  (m.eq.BRY_CONTACT(ieast,cr)%C2Bindex(j))
              IF(.not.(Uboundary.and.(irec.eq.indx1(ng)))) THEN
!>              OCEAN(ng)%tl_ubar(i,j,irec)=tl_my_value
!>
                ad_my_value=ad_my_value+                                &
     &                      OCEAN(ng)%ad_ubar(i,j,irec)
                OCEAN(ng)%ad_ubar(i,j,irec)=0.0_r8
!!            ELSE                                   ! for debugging
!!              OCEAN(ng)%ubar(i,j,irec)=0.0_r8      ! purposes
              END IF
# else
!>            OCEAN(ng)%tl_ubar(i,j,irec)=tl_my_value
!>
              ad_my_value=ad_my_value+                                  &
     &                    OCEAN(ng)%ad_ubar(i,j,irec)
              OCEAN(ng)%ad_ubar(i,j,irec)=0.0_r8
# endif
            END DO
# ifdef WET_DRY
!>          tl_my_value=tl_my_value*GRID(ng)%umask_wet(i,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%umask_wet(i,j)
# endif
# ifdef MASKING
!>          tl_my_value=tl_my_value*GRID(ng)%umask(i,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%umask(i,j)
# endif
!>          tl_my_value=Wold*                                           &
!>   &                  (Ucontact(cr)%Lweight(1,m)*                     &
!>   &                            REFINED(cr)%tl_ubar(1,m,told)+        &
!>   &                   Ucontact(cr)%Lweight(2,m)*                     &
!>   &                            REFINED(cr)%tl_ubar(2,m,told)+        &
!>   &                   Ucontact(cr)%Lweight(3,m)*                     &
!>   &                            REFINED(cr)%tl_ubar(3,m,told)+        &
!>   &                   Ucontact(cr)%Lweight(4,m)*                     &
!>   &                            REFINED(cr)%tl_ubar(4,m,told))        &
!>   &                  Wnew*                                           &
!>   &                  (Ucontact(cr)%Lweight(1,m)*                     &
!>   &                            REFINED(cr)%tl_ubar(1,m,tnew)+        &
!>   &                   Ucontact(cr)%Lweight(2,m)*                     &
!>   &                            REFINED(cr)%tl_ubar(2,m,tnew)+        &
!>   &                   Ucontact(cr)%Lweight(3,m)*                     &
!>   &                            REFINED(cr)%tl_ubar(3,m,tnew)+        &
!>   &                   Ucontact(cr)%Lweight(4,m)*                     &
!>   &                            REFINED(cr)%tl_ubar(4,m,tnew))
!>
            DO ii=1,4
              adfac1=Wold*Ucontact(cr)%Lweight(ii,m)*ad_my_value
              adfac2=Wnew*Ucontact(cr)%Lweight(ii,m)*ad_my_value
              REFINED(cr)%ad_ubar(ii,m,told)=                           &
     &                REFINED(cr)%ad_ubar(ii,m,told)+adfac1
              REFINED(cr)%ad_ubar(ii,m,tnew)=                           &
     &                REFINED(cr)%ad_ubar(ii,m,tnew)+adfac2
            END DO
            ad_my_value=0.0_r8
          END IF
        END DO

      END IF FREE_SURFACE
!
  10  FORMAT (/,'AD_PUT_REFINE2D - unbounded contact points temporal: ',&
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
     &        /,2x,  'time(ng) = ',i10,                                 &
     &          3x,'time(told) = ',i10,                                 &
     &          3x,'time(tnew) = ',i10)

      RETURN
      END SUBROUTINE ad_put_refine2d

# ifdef SOLVE3D
!
      SUBROUTINE ad_put_refine3d (ng, dg, cr, model, tile,              &
     &                            LBi, UBi, LBj, UBj)
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
!
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d, ad_mp_exchange4d
#  endif
      USE strings_mod,     ONLY : FoundError
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
      integer :: i, itrc, j, k, m, tnew, told, ii

      real(dp) :: Wnew, Wold, SecScale, fac
      real(r8) :: my_value, ad_my_value, adfac1, adfac2
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_put_refine3d"

#  include "set_bounds.h"
!
!  Clear adjoint constants.
!
      ad_my_value=0.0_r8
      adfac1=0.0_r8
      adfac2=0.0_r8
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
!  Set linear time interpolation weights. Fractional seconds are
!  rounded to the nearest milliseconds integer towards zero in the
!  time interpolation weights.
!
      SecScale=1000.0_dp              ! seconds to milliseconds
!
      Wold=ANINT((RollingTime(tnew,cr)-time(ng))*SecScale,dp)
      Wnew=ANINT((time(ng)-RollingTime(told,cr))*SecScale,dp)
      fac=1.0_dp/(Wold+Wnew)
      Wold=fac*Wold
      Wnew=fac*Wnew
!
!     IF (((Wold*Wnew).lt.0.0_dp).or.((Wold+Wnew).le.0.0_dp)) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) cr, dg, ng,                               &
     &                        iic(dg), told, tnew,                      &
     &                        iic(ng), Wold, Wnew,                      &
     &                        INT(time(ng)),                            &
     &                        INT(RollingTime(told,cr)),                &
     &                        INT(RollingTime(tnew,cr))
          END IF
!         exit_flag=8
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!     END IF

# ifdef NESTING_DEBUG
!
!  If debugging, write information into Fortran unit 202 to check the
!  logic of interpolating from donor grid data.
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        IF (Master) THEN
          IF (first) THEN
            first=.FALSE.
            WRITE (202,20)
          END IF
          WRITE (202,30) cr, dg, ng, iic(dg), iic(ng), told,  tnew,     &
     &                   INT(time(dg)),                                 &
     &                   INT(RollingTime(told,cr)),                     &
     &                   INT(time(ng)),                                 &
     &                   INT(RollingTime(tnew,cr)),                     &
     &                   Wold, Wnew
 20       FORMAT (3x,'cr',3x,'dg',3x,'ng',4x,'iic',4x,'iic',2x,'told',  &
     &            2x,'tnew',7x,'time',7x,'time',7x,'time',7x,'time',    &
     &            7x,'Wold',7x,'Wnew',/,18x,'(dg)',3x,'(ng)',           &
     &            19x,'(dg)',7x,'told',7x,'(ng)',7x,'tnew',/)
 30       FORMAT (3i5,2i7,2i6,4(2x,i9),2f11.4)
          CALL my_flush (202)
        END IF
      END IF
# endif

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile information.
!-----------------------------------------------------------------------
!
!>    CALL mp_exchange3d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    OCEAN(ng)%tl_u(:,:,:,1),                      &
!>   &                    OCEAN(ng)%tl_u(:,:,:,2),                      &
!>   &                    OCEAN(ng)%tl_v(:,:,:,1),                      &
!>   &                    OCEAN(ng)%tl_v(:,:,:,2))
!>
      CALL ad_mp_exchange3d (ng, tile, model, 4,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       OCEAN(ng)%ad_u(:,:,:,1),                   &
     &                       OCEAN(ng)%ad_u(:,:,:,2),                   &
     &                       OCEAN(ng)%ad_v(:,:,:,1),                   &
     &                       OCEAN(ng)%ad_v(:,:,:,2))
!>    CALL mp_exchange3d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    OCEAN(ng)%tl_u(:,:,:,1),                      &
!>   &                    OCEAN(ng)%tl_u(:,:,:,2),                      &
!>   &                    OCEAN(ng)%tl_v(:,:,:,1),                      &
!>   &                    OCEAN(ng)%tl_v(:,:,:,2))
!>
      CALL ad_mp_exchange4d (ng, tile, model, 3,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       OCEAN(ng)%ad_t(:,:,:,1,:),                 &
     &                       OCEAN(ng)%ad_t(:,:,:,2,:),                 &
     &                       OCEAN(ng)%ad_t(:,:,:,3,:))
#  endif
!
!  3D momentum in the XI-direction.
!
      DO m=1,Ucontact(cr)%Npoints
        i=Ucontact(cr)%Irg(m)
        j=Ucontact(cr)%Jrg(m)
        IF (((IstrP.le.i).and.(i.le.IendT)).and.                        &
     &      ((JstrT.le.j).and.(j.le.JendT))) THEN
          DO k=1,N(ng)
!>          OCEAN(ng)%tl_u(i,j,k,1)=tl_my_value
!>
            ad_my_value=ad_my_value+OCEAN(ng)%ad_u(i,j,k,1)
            OCEAN(ng)%ad_u(i,j,k,1)=0.0_r8
!>          OCEAN(ng)%tl_u(i,j,k,2)=tl_my_value
!>
            ad_my_value=ad_my_value+OCEAN(ng)%ad_u(i,j,k,2)
            OCEAN(ng)%ad_u(i,j,k,2)=0.0_r8
#  ifdef MASKING
!>          tl_my_value=tl_my_value*GRID(ng)%umask(i,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%umask(i,j)
#  endif
!>          tl_my_value=Wold*                                           &
!>   &               (Ucontact(cr)%Lweight(1,m)*                        &
!>   &                         REFINED(cr)%tl_u(1,k,m,told)+            &
!>   &                Ucontact(cr)%Lweight(2,m)*                        &
!>   &                         REFINED(cr)%tl_u(2,k,m,told)+            &
!>   &                Ucontact(cr)%Lweight(3,m)*                        &
!>   &                         REFINED(cr)%tl_u(3,k,m,told)+            &
!>   &                Ucontact(cr)%Lweight(4,m)*                        &
!>   &                         REFINED(cr)%tl_u(4,k,m,told))+           &
!>   &               Wnew*                                              &
!>   &               (Ucontact(cr)%Lweight(1,m)*                        &
!>   &                         REFINED(cr)%tl_u(1,k,m,tnew)+            &
!>   &                Ucontact(cr)%Lweight(2,m)*                        &
!>   &                         REFINED(cr)%tl_u(2,k,m,tnew)+            &
!>   &                Ucontact(cr)%Lweight(3,m)*                        &
!>   &                         REFINED(cr)%tl_u(3,k,m,tnew)+            &
!>   &                Ucontact(cr)%Lweight(4,m)*                        &
!>   &                         REFINED(cr)%tl_u(4,k,m,tnew))
            DO ii=1,4
              adfac1=Wold*Ucontact(cr)%Lweight(ii,m)*ad_my_value
              adfac2=Wnew*Ucontact(cr)%Lweight(ii,m)*ad_my_value
              REFINED(cr)%ad_u(ii,k,m,told)=                            &
     &                REFINED(cr)%ad_u(ii,k,m,told)+adfac1
              REFINED(cr)%ad_u(ii,k,m,tnew)=                            &
     &                REFINED(cr)%ad_u(ii,k,m,tnew)+adfac2
            END DO
            ad_my_value=0.0_r8
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
!>          OCEAN(ng)%tl_v(i,j,k,1)=tl_my_value
!>
            ad_my_value=ad_my_value+OCEAN(ng)%ad_v(i,j,k,1)
            OCEAN(ng)%ad_v(i,j,k,1)=0.0_r8
!>          OCEAN(ng)%tl_v(i,j,k,2)=tl_my_value
!>
            ad_my_value=ad_my_value+OCEAN(ng)%ad_v(i,j,k,2)
            OCEAN(ng)%ad_v(i,j,k,2)=0.0_r8
#  ifdef MASKING
!>          tl_my_value=tl_my_value*GRID(ng)%vmask(i,j)
!>
            ad_my_value=ad_my_value*GRID(ng)%vmask(i,j)
#  endif
!>          tl_my_value=Wold*                                           &
!>   &                  (Vcontact(cr)%Lweight(1,m)*                     &
!>   &                            REFINED(cr)%tl_v(1,k,m,told)+         &
!>   &                   Vcontact(cr)%Lweight(2,m)*                     &
!>   &                            REFINED(cr)%tl_v(2,k,m,told)+         &
!>   &                   Vcontact(cr)%Lweight(3,m)*                     &
!>   &                            REFINED(cr)%tl_v(3,k,m,told)+         &
!>   &                   Vcontact(cr)%Lweight(4,m)*                     &
!>   &                            REFINED(cr)%tl_v(4,k,m,told))+        &
!>   &                  Wnew*                                           &
!>   &                  (Vcontact(cr)%Lweight(1,m)*                     &
!>   &                            REFINED(cr)%tl_v(1,k,m,tnew)+         &
!>   &                   Vcontact(cr)%Lweight(2,m)*                     &
!>   &                            REFINED(cr)%tl_v(2,k,m,tnew)+         &
!>   &                   Vcontact(cr)%Lweight(3,m)*                     &
!>   &                            REFINED(cr)%tl_v(3,k,m,tnew)+         &
!>   &                   Vcontact(cr)%Lweight(4,m)*                     &
!>   &                            REFINED(cr)%tl_v(4,k,m,tnew))
!>
            DO ii=1,4
              adfac1=Wold*Vcontact(cr)%Lweight(ii,m)*ad_my_value
              adfac2=Wnew*Vcontact(cr)%Lweight(ii,m)*ad_my_value
              REFINED(cr)%ad_v(ii,k,m,told)=                            &
     &                REFINED(cr)%ad_v(ii,k,m,told)+adfac1
              REFINED(cr)%ad_v(ii,k,m,tnew)=                            &
     &                REFINED(cr)%ad_v(ii,k,m,tnew)+adfac2
            END DO
            ad_my_value=0.0_r8
          END DO
        END IF
      END DO
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
!>            OCEAN(ng)%tl_t(i,j,k,1,itrc)=tl_my_value
!>
              ad_my_value=ad_my_value+OCEAN(ng)%ad_t(i,j,k,1,itrc)
              OCEAN(ng)%ad_t(i,j,k,1,itrc)=0.0_r8
!>            OCEAN(ng)%tl_t(i,j,k,2,itrc)=tl_my_value
!>
              ad_my_value=ad_my_value+OCEAN(ng)%ad_t(i,j,k,2,itrc)
              OCEAN(ng)%ad_t(i,j,k,2,itrc)=0.0_r8
!>            OCEAN(ng)%tl_t(i,j,k,3,itrc)=tl_my_value
!>a
              ad_my_value=ad_my_value+OCEAN(ng)%ad_t(i,j,k,3,itrc)
              OCEAN(ng)%ad_t(i,j,k,3,itrc)=0.0_r8
#  ifdef MASKING
!>            tl_my_value=tl_my_value*GRID(ng)%rmask(i,j)
              ad_my_value=ad_my_value*GRID(ng)%rmask(i,j)
#  endif
!>            tl_my_value=Wold*                                         &
!>   &                    (Rcontact(cr)%Lweight(1,m)*                   &
!>   &                              REFINED(cr)%tl_t(1,k,m,told,itrc)+  &
!>   &                     Rcontact(cr)%Lweight(2,m)*                   &
!>   &                              REFINED(cr)%tl_t(2,k,m,told,itrc)+  &
!>   &                     Rcontact(cr)%Lweight(3,m)*                   &
!>   &                              REFINED(cr)%tl_t(3,k,m,told,itrc)+  &
!>   &                     Rcontact(cr)%Lweight(4,m)*                   &
!>   &                              REFINED(cr)%tl_t(4,k,m,told,itrc))+ &
!>   &                    Wnew*                                         &
!>   &                    (Rcontact(cr)%Lweight(1,m)*                   &
!>   &                              REFINED(cr)%tl_t(1,k,m,tnew,itrc)+  &
!>   &                     Rcontact(cr)%Lweight(2,m)*                   &
!>   &                              REFINED(cr)%tl_t(2,k,m,tnew,itrc)+  &
!>   &                     Rcontact(cr)%Lweight(3,m)*                   &
!>   &                              REFINED(cr)%tl_t(3,k,m,tnew,itrc)+  &
!>   &                     Rcontact(cr)%Lweight(4,m)*                   &
!>   &                              REFINED(cr)%tl_t(4,k,m,tnew,itrc))
!>
              DO ii=1,4
                adfac1=Wold*Rcontact(cr)%Lweight(ii,m)*ad_my_value
                adfac2=Wnew*Rcontact(cr)%Lweight(ii,m)*ad_my_value
                REFINED(cr)%ad_t(ii,k,m,told,itrc)=                     &
     &                  REFINED(cr)%ad_t(ii,k,m,told,itrc)+adfac1
                REFINED(cr)%ad_t(ii,k,m,tnew,itrc)=                     &
     &                  REFINED(cr)%ad_t(ii,k,m,tnew,itrc)+adfac2
              END DO
              ad_my_value=0.0_r8
            END DO
          END DO
        END IF
      END DO
!
  10  FORMAT (/,'AD_PUT_REFINE3D - unbounded contact points temporal: ',&
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
     &        /,2x,  'time(ng) = ',i10,                                 &
     &          3x,'time(told) = ',i10,                                 &
     &          3x,'time(tnew) = ',i10)

      RETURN
      END SUBROUTINE ad_put_refine3d
# endif

# ifdef SOLVE3D
!
      SUBROUTINE ad_z_weights (ng, model, tile)
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
!
#  ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_assemble
#  endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
!
!  Local variable declarations.
!
      integer :: cr, dg, rg, i, j, k, m, ii
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
      real(r8) :: ad_Zbot, ad_Zr, ad_Ztop, ad_dz, ad_r1, ad_r2
      real(r8) :: adfac, adfac1

      real(r8), allocatable :: Zd(:,:,:)
      real(r8), allocatable :: ad_Zd(:,:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_z_weights"
!
!=======================================================================
!  Adjoint compute vertical indices and weights for each contact region.
!=======================================================================
!
!
!  Clear adjoint constants.
!
      ad_Zbot=0.0_r8
      ad_Zr=0.0_r8
      ad_Ztop=0.0_r8
      ad_dz=0.0_r8
      ad_r1=0.0_r8
      ad_r2=0.0_r8
      adfac=0.0_r8
!
!  Compute vertical indices and weights.
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

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes. No action required for the
!  adjoint of mp_assemble (AMM).
!
          Nkpts=N(rg)*Npoints
          Nwpts=2*Nkpts
          Nzpts=4*Nkpts
!
!>        CALL mp_assemble (rg, model, Nkpts, ispv, Vcontact(cr)%Kdg)
!>
!!        CALL ad_mp_assemble (rg, model, Nkpts, ispv,                  &
!!                             Vcontact(cr)%Kdg)
!!        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

!!        CALL mp_assemble (rg, model, Nwpts, spv,                      &
!!   &                      Vcontact(cr)%Vweight)
!!        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

!!        CALL ad_mp_assemble (rg, model, Nwpts, spv,                   &
!!   &                         Vcontact(cr)%ad_Vweight)
!!        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
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

          V_CONTACT : IF (.not.Vcontact(cr)%interpolate.and.            &
     &                    Vcontact(cr)%coincident) THEN
            DO Krg=1,N(rg)
              DO m=1,Npoints
                Irg=Vcontact(cr)%Irg(m)
                Jrg=Vcontact(cr)%Jrg(m)
                IF (((IminR.le.Irg).and.(Irg.le.ImaxR)).and.            &
     &              ((JminR.le.Jrg).and.(Jrg.le.JmaxR))) THEN
!>                Vcontact(cr)%Kdg(Krg,m)=Krg
!>                Vcontact(cr)%Vweight(1,Krg,m)=1.0_r8
!>                Vcontact(cr)%Vweight(2,Krg,m)=0.0_r8

!>                Vcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
!>
                  Vcontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                Vcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
!>
                  Vcontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
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
            IF (.not.allocated(ad_Zd)) THEN
              allocate (ad_Zd(4,N(dg),Npoints))
            END IF
            ad_Zd=0.0_r8
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
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
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
!>                  Vcontact(cr)%Kdg(Krg,m)=N(dg)! donor grid cell value
!>                  Vcontact(cr)%Vweight(1,Krg,m)=0.0_r8
!>                  Vcontact(cr)%Vweight(2,Krg,m)=1.0_r8

!>                  Vcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
!>
                    Vcontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                  Vcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
!>
                    Vcontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
                  ELSE IF (Zbot.ge.Zr) THEN      ! If deeper, use bottom
!>                  Vcontact(cr)%Kdg(Krg,m)=1    ! donor grid cell value
!>                  Vcontact(cr)%Vweight(1,Krg,m)=0.0_r8
!>                  Vcontact(cr)%Vweight(2,Krg,m)=1.0_r8

!>                  Vcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
!>
                    Vcontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                  Vcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
!>
                    Vcontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
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
!>                      Vcontact(cr)%Kdg(Krg,m)=Kdg
!>                      Vcontact(cr)%Vweight(1,Krg,m)=r1
!>                      Vcontact(cr)%Vweight(2,Krg,m)=r2

!>                      Vcontact(cr)%tl_Vweight(1,Krg,m)=tl_r1
!>
                        ad_r1=ad_r1+Vcontact(cr)%ad_Vweight(1,Krg,m)
                        Vcontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                      Vcontact(cr)%tl_Vweight(2,Krg,m)=tl_r2
!>
                        ad_r2=ad_r2+Vcontact(cr)%ad_Vweight(2,Krg,m)
                        Vcontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
!>                      tl_r1=-tl_r2
!>
                        ad_r2=ad_r2-ad_r1
                        ad_r1=0.0_r8
!>                      tl_r2=(tl_Zr-tl_Zbot)/dz-tl_dz*r2/dz
!>
                        adfac=ad_r2/dz
                        ad_Zr=ad_Zr+adfac
                        ad_Zbot=ad_Zbot-adfac
                        ad_dz=ad_dz-r2*adfac
                        ad_r2=0.0_r8

!>                      tl_dz=tl_Ztop-tl_Zbot
!>
                        ad_Ztop=ad_Ztop+ad_dz
                        ad_Zbot=ad_Zbot-ad_dz
                        ad_dz=0.0_r8
                      END IF
!>                    tl_Zbot=Vcontact(cr)%Lweight(1,m)*                &
!>   &                                     tl_Zd(1,Kdg-1,m)+            &
!>   &                        Vcontact(cr)%Lweight(2,m)*                &
!>   &                                     tl_Zd(2,Kdg-1,m)+            &
!>   &                        Vcontact(cr)%Lweight(3,m)*                &
!>   &                                     tl_Zd(3,Kdg-1,m)+            &
!>   &                        Vcontact(cr)%Lweight(4,m)*                &
!>   &                                     tl_Zd(4,Kdg-1,m)
!>
                      DO ii=1,4
                        adfac=Vcontact(cr)%Lweight(ii,m)*ad_Zbot
                        ad_Zd(ii,Kdg-1,m)=ad_Zd(ii,Kdg-1,m)+adfac
                      END DO
                      ad_Zbot=0.0_r8
!>                    tl_Ztop=Vcontact(cr)%Lweight(1,m)*                &
!>   &                                     tl_Zd(1,Kdg  ,m)+            &
!>   &                        Vcontact(cr)%Lweight(2,m)*                &
!>   &                                     tl_Zd(2,Kdg  ,m)+            &
!>   &                        Vcontact(cr)%Lweight(3,m)*                &
!>   &                                     tl_Zd(3,Kdg  ,m)+            &
!>   &                        Vcontact(cr)%Lweight(4,m)*
!>   &                                     tl_Zd(4,Kdg  ,m)
!>
                      DO ii=1,4
                        adfac=Vcontact(cr)%Lweight(ii,m)*ad_Ztop
                        ad_Zd(ii,Kdg  ,m)=ad_Zd(ii,Kdg  ,m)+adfac
                      END DO
                      ad_Ztop=0.0_r8
                    END DO
                  END IF
!>                tl_Zr=0.5_r8*                                         &
!>                      (GRID(rg)%tl_z_r(Irg,Jrg  ,Krg)+                &
!>   &                   GRID(rg)%tl_z_r(Irg,Jrg-1,Krg))
!>
                  GRID(rg)%ad_z_r(Irg,Jrg  ,Krg)=                       &
     &                     GRID(rg)%ad_z_r(Irg,Jrg  ,Krg)+0.5_r8*ad_Zr
                  GRID(rg)%ad_z_r(Irg,Jrg-1,Krg)=                       &
     &                     GRID(rg)%ad_z_r(Irg,Jrg-1,Krg)+0.5_r8*ad_Zr
                  ad_Zr=0.0_r8
!>                tl_Zbot=Vcontact(cr)%Lweight(1,m)*tl_Zd(1,1    ,m)+   &
!>   &                    Vcontact(cr)%Lweight(2,m)*tl_Zd(2,1    ,m)+   &
!>   &                    Vcontact(cr)%Lweight(3,m)*tl_Zd(3,1    ,m)+   &
!>   &                    Vcontact(cr)%Lweight(4,m)*tl_Zd(4,1    ,m)
!>
                  DO ii=1,4
                    adfac=Vcontact(cr)%Lweight(ii,m)*ad_Zbot
                    ad_Zd(ii,1    ,m)=ad_Zd(ii,1    ,m)+adfac
                  END DO
                  ad_Zbot=0.0_r8
!>                tl_Ztop=Vcontact(cr)%Lweight(1,m)*tl_Zd(1,N(dg),m)+   &
!>   &                    Vcontact(cr)%Lweight(2,m)*tl_Zd(2,N(dg),m)+   &
!>   &                    Vcontact(cr)%Lweight(3,m)*tl_Zd(3,N(dg),m)+   &
!>   &                    Vcontact(cr)%Lweight(4,m)*tl_Zd(4,N(dg),m)
!>
                  DO ii=1,4
                    adfac=Vcontact(cr)%Lweight(ii,m)*ad_Ztop
                    ad_Zd(ii,N(dg),m)=ad_Zd(ii,N(dg),m)+adfac
                  END DO
                  ad_Ztop=0.0_r8
                END IF
              END DO
            END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes. No action required for the
!  adjoint of mp_assemble (AMM).
!
!>          CALL mp_assemble (dg, model, Nzpts, spv, Zd)
!>
!!          CALL ad_mp_assemble (dg, model, Nzpts, spv, ad_Zd)
!!          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
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
                  adfac=0.5_r8*ad_Zd(1,Kdg,m)
!>                tl_Zd(1,Kdg,m)=0.5_r8*                                &
!>   &                           (GRID(dg)%tl_z_r(Idg  ,Jdgm1,Kdg)+     &
!>   &                            GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg))
!>
                  GRID(dg)%ad_z_r(Idg  ,Jdgm1,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idg  ,Jdgm1,Kdg)+adfac
                  GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)+adfac
                  ad_Zd(1,Kdg,m)=0.0
!>                tl_Zd(2,Kdg,m)=0.5_r8*                                &
!>   &                           (GRID(dg)%tl_z_r(Idgp1,Jdgm1,Kdg)+     &
!>   &                            GRID(dg)%tl_z_r(Idgp1,Jdg  ,Kdg))
!>
                  adfac=0.5_r8*ad_Zd(2,Kdg,m)
                  GRID(dg)%ad_z_r(Idgp1,Jdgm1,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idgp1,Jdgm1,Kdg)+adfac
                  GRID(dg)%ad_z_r(Idgp1,Jdg  ,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idgp1,Jdg  ,Kdg)+adfac
                  ad_Zd(2,Kdg,m)=0.0_r8
!>                tl_Zd(3,Kdg,m)=0.5_r8*                                &
!>   &                           (GRID(dg)%tl_z_r(Idgp1,Jdg  ,Kdg)+     &
!>   &                            GRID(dg)%tl_z_r(Idgp1,Jdgp1,Kdg))
!>
                  adfac=0.5_r8*ad_Zd(3,Kdg,m)
                  GRID(dg)%ad_z_r(Idgp1,Jdg  ,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idgp1,Jdg  ,Kdg)+adfac
                  GRID(dg)%ad_z_r(Idgp1,Jdgp1,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idgp1,Jdgp1,Kdg)+adfac
                  ad_Zd(3,Kdg,m)=0.0_r8
!>                tl_Zd(4,Kdg,m)=0.5_r8*                                &
!>   &                           (GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg)+     &
!>   &                            GRID(dg)%tl_z_r(Idg  ,Jdgp1,Kdg))
!>
                  adfac=0.5_r8*ad_Zd(4,Kdg,m)
                  GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)+adfac
                  GRID(dg)%ad_z_r(Idg  ,Jdgp1,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idg  ,Jdgp1,Kdg)+adfac
                  ad_Zd(4,Kdg,m)=0.0_r8
                END IF
              END DO
            END DO
          END IF V_CONTACT

#  ifdef DISTRIBUTE
!
!  If distributed-memory, initialize with special value (zero) to
!  facilitate the global reduction when collecting data between all
!  nodes.
!
!>        Vcontact(cr)%tl_Vweight(1:2,1:N(rg),1:Npoints)=0.0_r8
!>
          Vcontact(cr)%ad_Vweight(1:2,1:N(rg),1:Npoints)=0.0_r8
#  endif
!
!  Deallocate local work arrays.
!
          IF (allocated(Zd)) THEN
            deallocate (Zd)
          END IF
          IF (allocated(ad_Zd)) THEN
            deallocate (ad_Zd)
          END IF

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes. No action required for the
!  adjoint of mp_assemble (AMM).
!
!>        CALL mp_assemble (rg, model, Nwpts, spv,                      &
!>   &                      Ucontact(cr)%Vweight)
!>
!!        CALL ad_mp_assemble (rg, model, Nwpts, spv,                   &
!!   &                         Ucontact(cr)%ad_Vweight)
!!        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
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
!>                Ucontact(cr)%Kdg(Krg,m)=Krg
!>                Ucontact(cr)%Vweight(1,Krg,m)=1.0_r8
!>                Ucontact(cr)%Vweight(2,Krg,m)=0.0_r8

!>                Ucontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
!>
                  Ucontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                Ucontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
!>
                  Ucontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
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
            IF (.not.allocated(ad_Zd)) THEN
              allocate (ad_Zd(4,N(dg),Npoints))
            END IF
            ad_Zd=0.0_r8
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
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!>                  Ucontact(cr)%Kdg(Krg,m)=N(dg)! donor grid cell value
!>                  Ucontact(cr)%Vweight(1,Krg,m)=0.0_r8
!>                  Ucontact(cr)%Vweight(2,Krg,m)=1.0_r8

!>                  Ucontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
!>
                    Ucontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                  Ucontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
!>
                    Ucontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
                  ELSE IF (Zbot.ge.Zr) THEN      ! If deeper, use bottom
!>                  Ucontact(cr)%Kdg(Krg,m)=1    ! donor grid cell value
!>                  Ucontact(cr)%Vweight(1,Krg,m)=0.0_r8
!>                  Ucontact(cr)%Vweight(2,Krg,m)=1.0_r8

!>                  Ucontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
!>
                    Ucontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                  Ucontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
!>
                    Ucontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
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
!>                      Ucontact(cr)%Kdg(Krg,m)=Kdg
!>                      Ucontact(cr)%Vweight(1,Krg,m)=r1
!>                      Ucontact(cr)%Vweight(2,Krg,m)=r2

!>                      Ucontact(cr)%tl_Vweight(1,Krg,m)=tl_r1
!>
                        ad_r1=ad_r1+Ucontact(cr)%ad_Vweight(1,Krg,m)
                        Ucontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                      Ucontact(cr)%tl_Vweight(2,Krg,m)=tl_r2
!>
                        ad_r2=ad_r2+Ucontact(cr)%ad_Vweight(2,Krg,m)
                        Ucontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
!>                      tl_r1=-tl_r2
!>
                        ad_r2=ad_r2-ad_r1
                        ad_r1=0.0_r8
!>                      tl_r2=(tl_Zr-tl_Zbot)/dz-tl_dz*r2/dz
!>
                        adfac=ad_r2/dz
                        ad_Zr=ad_Zr+adfac
                        ad_Zbot=ad_Zbot-adfac
                        ad_dz=ad_dz-r2*adfac
                        ad_r2=0.0_r8
!>                      tl_dz=tl_Ztop-tl_Zbot
!>
                        ad_Ztop=ad_Ztop+ad_dz
                        ad_Zbot=ad_Zbot-ad_dz
                        ad_dz=0.0_r8
                      END IF
!>                    tl_Zbot=Ucontact(cr)%Lweight(1,m)*                &
!>   &                                     tl_Zd(1,Kdg-1,m)+            &
!>   &                        Ucontact(cr)%Lweight(2,m)*                &
!>   &                                     tl_Zd(2,Kdg-1,m)+            &
!>   &                        Ucontact(cr)%Lweight(3,m)*                &
!>   &                                     tl_Zd(3,Kdg-1,m)+            &
!>   &                        Ucontact(cr)%Lweight(4,m)*
!>   &                                     tl_Zd(4,Kdg-1,m)
!>
                      DO ii=1,4
                        adfac=Ucontact(cr)%Lweight(ii,m)*ad_Zbot
                        ad_Zd(ii,Kdg-1,m)=ad_Zd(ii,Kdg-1,m)+adfac
                      END DO
                      ad_Zbot=0.0_r8
!>                    tl_Ztop=Ucontact(cr)%Lweight(1,m)*                &
!>   &                                     tl_Zd(1,Kdg  ,m)+            &
!>   &                        Ucontact(cr)%Lweight(2,m)*                &
!>   &                                     tl_Zd(2,Kdg  ,m)+            &
!>   &                        Ucontact(cr)%Lweight(3,m)*                &
!>   &                                     tl_Zd(3,Kdg  ,m)+            &
!>   &                        Ucontact(cr)%Lweight(4,m)*
!>   &                                     tl_Zd(4,Kdg  ,m)
                      DO ii=1,4
                        adfac=Ucontact(cr)%Lweight(ii,m)*ad_Ztop
                        ad_Zd(ii,Kdg  ,m)=ad_Zd(ii,Kdg  ,m)+adfac
                      END DO
                      ad_Ztop=0.0_r8
                    END DO
                  END IF
!>                tl_Zr=0.5_r8*(GRID(rg)%tl_z_r(Irg  ,Jrg,Krg)+         &
!>   &                          GRID(rg)%tl_z_r(Irg-1,Jrg,Krg))
!>
                  adfac=0.5_r8*ad_Zr
                  GRID(rg)%ad_z_r(Irg  ,Jrg,Krg)=                       &
     &                     GRID(rg)%ad_z_r(Irg  ,Jrg,Krg)+adfac
                  GRID(rg)%ad_z_r(Irg-1,Jrg,Krg)=                       &
     &                     GRID(rg)%ad_z_r(Irg-1,Jrg,Krg)+adfac
                  ad_Zr=0.0_r8
!>                tl_Zbot=Ucontact(cr)%Lweight(1,m)*tl_Zd(1,1    ,m)+   &
!>   &                    Ucontact(cr)%Lweight(2,m)*tl_Zd(2,1    ,m)+   &
!>   &                    Ucontact(cr)%Lweight(3,m)*tl_Zd(3,1    ,m)+   &
!>   &                    Ucontact(cr)%Lweight(4,m)*tl_Zd(4,1    ,m)
!>
                  DO ii=1,4
                    adfac=Ucontact(cr)%Lweight(ii,m)*ad_Zbot
                    ad_Zd(ii,1    ,m)=ad_Zd(ii,1    ,m)+adfac
                  END DO
                  ad_Zbot=0.0_r8
!>                tl_Ztop=Ucontact(cr)%Lweight(1,m)*tl_Zd(1,N(dg),m)+   &
!>   &                    Ucontact(cr)%Lweight(2,m)*tl_Zd(2,N(dg),m)+   &
!>   &                    Ucontact(cr)%Lweight(3,m)*tl_Zd(3,N(dg),m)+   &
!>   &                    Ucontact(cr)%Lweight(4,m)*tl_Zd(4,N(dg),m)
!>
                  DO ii=1,4
                    adfac=Ucontact(cr)%Lweight(ii,m)*ad_Ztop
                    ad_Zd(ii,N(dg),m)=ad_Zd(ii,N(dg),m)+adfac
                  END DO
                  ad_Ztop=0.0_r8
                END IF
              END DO
            END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.  No action required for the
!  adjoint of mp_assemble (AMM)
!
!>          CALL mp_assemble (dg, model, Nzpts, spv, ad_Zd)
!>
!!          CALL ad_mp_assemble (dg, model, Nzpts, spv, ad_Zd)
!!          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
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
!>                tl_Zd(1,Kdg,m)=0.5_r8*                                &
!>   &                           (GRID(dg)%tl_z_r(Idgm1,Jdg  ,Kdg)+     &
!>   &                            GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg))
!>
                  adfac=0.5_r8*ad_Zd(1,Kdg,m)
                  GRID(dg)%ad_z_r(Idgm1,Jdg  ,Kdg)=                     &
     &                  GRID(dg)%ad_z_r(Idgm1,Jdg  ,Kdg)+adfac
                  GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)=                     &
     &                  GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)+adfac
                  ad_Zd(1,Kdg,m)=0.0_r8
!>                tl_Zd(2,Kdg,m)=0.5_r8*                                &
!>   &                           (GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg)+     &
!>   &                            GRID(dg)%tl_z_r(Idgp1,Jdg  ,Kdg))
!>
                  adfac=0.5_r8*ad_Zd(2,Kdg,m)
                  GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)=                     &
     &                  GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)+adfac
                  GRID(dg)%ad_z_r(Idgp1,Jdg  ,Kdg)=                     &
     &                  GRID(dg)%ad_z_r(Idgp1,Jdg  ,Kdg)+adfac
                  ad_Zd(2,Kdg,m)=0.0_r8
!>                tl_Zd(3,Kdg,m)=0.5_r8*                                &
!>   &                           (GRID(dg)%tl_z_r(Idg  ,Jdgp1,Kdg)+     &
!>   &                            GRID(dg)%tl_z_r(Idgp1,Jdgp1,Kdg))
!>
                  adfac=0.5_r8*ad_Zd(3,Kdg,m)
                  GRID(dg)%ad_z_r(Idg  ,Jdgp1,Kdg)=                     &
     &                  GRID(dg)%ad_z_r(Idg  ,Jdgp1,Kdg)+adfac
                  GRID(dg)%ad_z_r(Idgp1,Jdgp1,Kdg)=                     &
     &                  GRID(dg)%ad_z_r(Idgp1,Jdgp1,Kdg)+adfac
                  ad_Zd(3,Kdg,m)=0.0_r8
!>                tl_Zd(4,Kdg,m)=0.5_r8*                                &
!>   &                           (GRID(dg)%tl_z_r(Idgm1,Jdgp1,Kdg)+     &
!>   &                            GRID(dg)%tl_z_r(Idg  ,Jdgp1,Kdg))
!>
                  adfac=0.5_r8*ad_Zd(4,Kdg,m)
                  GRID(dg)%ad_z_r(Idgm1,Jdgp1,Kdg)=                     &
     &                  GRID(dg)%ad_z_r(Idgm1,Jdgp1,Kdg)+adfac
                  GRID(dg)%ad_z_r(Idg  ,Jdgp1,Kdg)=                     &
     &                  GRID(dg)%ad_z_r(Idg  ,Jdgp1,Kdg)+adfac
                  ad_Zd(4,Kdg,m)=0.0_r8
                END IF
              END DO
            END DO
         END IF U_CONTACT

#  ifdef DISTRIBUTE
!
!  If distributed-memory, initialize with special value (zero) to
!  facilitate the global reduction when collecting data between all
!  nodes.
!
!>        Ucontact(cr)%tl_Vweight(1:2,1:N(rg),1:Npoints)=0.0_r8
!>
          Ucontact(cr)%ad_Vweight(1:2,1:N(rg),1:Npoints)=0.0_r8
#  endif
!
!  Deallocate local work arrays.
!
          IF (allocated(Zd)) THEN
            deallocate (Zd)
          END IF
          IF (allocated(ad_Zd)) THEN
            deallocate (ad_Zd)
          END IF

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes. No action required for
!  the adjoint of mp_assemble.
!
          Npoints=Rcontact(cr)%Npoints
          Nkpts=N(rg)*Npoints
          Nwpts=2*Nkpts
          Nzpts=4*Nkpts

!>        CALL _mp_assemble (rg, model, Nwpts, spv,                     &
!>   &                       Rcontact(cr)Vweight)
!>
!!        CALL ad_mp_assemble (rg, model, Nwpts, spv,                   &
!!   &                         Rcontact(cr)%ad_Vweight)
!!        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
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
!>                Rcontact(cr)%Kdg(Krg,m)=Krg
!>                Rcontact(cr)%Vweight(1,Krg,m)=1.0_r8
!>                Rcontact(cr)%Vweight(2,Krg,m)=0.0_r8

!>                Rcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
!>
                  Rcontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                Rcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
!>
                  Rcontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
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
            IF (.not.allocated(ad_Zd)) THEN
              allocate ( ad_Zd(4,N(dg),Npoints) )
            END IF
            ad_Zd=0.0_r8
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
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!>                  Rcontact(cr)%Kdg(Krg,m)=N(dg)! donor grid cell value
!>                  Rcontact(cr)%Vweight(1,Krg,m)=0.0_r8
!>                  Rcontact(cr)%Vweight(2,Krg,m)=1.0_r8

!>                  Rcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
!>
                    Rcontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                  Rcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
!>
                    Rcontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
                  ELSE IF (Zbot.ge.Zr) THEN      ! If deeper, use bottom
!>                  Rcontact(cr)%Kdg(Krg,m)=1    ! donor grid cell value
!>                  Rcontact(cr)%Vweight(1,Krg,m)=0.0_r8
!>                  Rcontact(cr)%Vweight(2,Krg,m)=1.0_r8

!>                  Rcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
!>
                    Rcontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                  Rcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
!>
                    Rcontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
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
!>                      Rcontact(cr)%Kdg(Krg,m)=Kdg
!>                      Rcontact(cr)%Vweight(1,Krg,m)=r1
!>                      Rcontact(cr)%Vweight(2,Krg,m)=r2

!>                      Rcontact(cr)%tl_Vweight(1,Krg,m)=tl_r1
!>
                        ad_r1=ad_r1+Rcontact(cr)%ad_Vweight(1,Krg,m)
                        Rcontact(cr)%ad_Vweight(1,Krg,m)=0.0_r8
!>                      Rcontact(cr)%tl_Vweight(2,Krg,m)=tl_r2
!>
                        ad_r2=ad_r2+Rcontact(cr)%ad_Vweight(2,Krg,m)
                        Rcontact(cr)%ad_Vweight(2,Krg,m)=0.0_r8
!>                      tl_r1=-tl_r2
!>
                        ad_r2=ad_r2-ad_r1
                        ad_r1=0.0_r8
!>                      tl_r2=(tl_Zr-tl_Zbot)/dz-tl_dz*r2/dz
!>
                        adfac=ad_r1/dz
                        ad_Zr=ad_Zr+adfac
                        ad_Zbot=ad_Zbot-adfac
                        ad_dz=ad_dz-r2*adfac
                        ad_r2=0.0_r8
!>                      tl_dz=tl_Ztop-tl_Zbot
!>
                        ad_Ztop=ad_Ztop+ad_dz
                        ad_Zbot=ad_Zbot-ad_dz
                        ad_dz=0.0_r8
                      END IF
!>                    tl_Zbot=Rcontact(cr)%Lweight(1,m)*                &
!>   &                                     tl_Zd(1,Kdg-1,m)+            &
!>   &                        Rcontact(cr)%Lweight(2,m)*                &
!>   &                                     tl_Zd(2,Kdg-1,m)+            &
!>   &                        Rcontact(cr)%Lweight(3,m)*                &
!>   &                                     tl_Zd(3,Kdg-1,m)+            &
!>   &                        Rcontact(cr)%Lweight(4,m)*                &
!>   &                                     tl_Zd(4,Kdg-1,m)
!>
                      DO ii=1,4
                        adfac=Rcontact(cr)%Lweight(ii,m)*ad_Zbot
                        ad_Zd(ii,Kdg-1,m)=ad_Zd(ii,Kdg-1,m)+adfac
                      END DO
                      ad_Zbot=0.0_r8
!>                    tl_Ztop=Rcontact(cr)%Lweight(1,m)*                &
!>   &                                     tl_Zd(1,Kdg  ,m)+            &
!>   &                        Rcontact(cr)%Lweight(2,m)*                &
!>   &                                     tl_Zd(2,Kdg  ,m)+            &
!>   &                        Rcontact(cr)%Lweight(3,m)*                &
!>   &                                     tl_Zd(3,Kdg  ,m)+            &
!>   &                        Rcontact(cr)%Lweight(4,m)*
!>   &                                     tl_Zd(4,Kdg  ,m)
!>
                      DO ii=1,4
                        adfac=Rcontact(cr)%Lweight(ii,m)*ad_Ztop
                        ad_Zd(ii,Kdg  ,m)=ad_Zd(ii,Kdg  ,m)+adfac
                      END DO
                      ad_Ztop=0.0_r8
                    END DO
                  END IF
!>                tl_Zr=GRID(rg)%tl_z_r(Irg,Jrg,Krg)
!>
                  GRID(rg)%ad_z_r(Irg,Jrg,Krg)=                         &
     &                     GRID(rg)%ad_z_r(Irg,Jrg,Krg)+ad_Zr
                  ad_Zr=0.0_r8
!>                tl_Ztop=Rcontact(cr)%Lweight(1,m)*tl_Zd(1,N(dg),m)+   &
!>   &                    Rcontact(cr)%Lweight(2,m)*tl_Zd(2,N(dg),m)+   &
!>   &                    Rcontact(cr)%Lweight(3,m)*tl_Zd(3,N(dg),m)+   &
!>   &                    Rcontact(cr)%Lweight(4,m)*tl_Zd(4,N(dg),m)
!>
                  DO ii=1,4
                    adfac=Rcontact(cr)%Lweight(ii,m)*ad_Ztop
                    ad_Zd(ii,N(dg),m)=ad_Zd(ii,N(dg),m)+adfac
                  END DO
                  ad_Ztop=0.0_r8
!>                tl_Zbot=Rcontact(cr)%Lweight(1,m)*tl_Zd(1,1    ,m)+   &
!>   &                    Rcontact(cr)%Lweight(2,m)*tl_Zd(2,1    ,m)+   &
!>   &                    Rcontact(cr)%Lweight(3,m)*tl_Zd(3,1    ,m)+   &
!>   &                    Rcontact(cr)%Lweight(4,m)*tl_Zd(4,1    ,m)
!>
                  DO ii=1,4
                    adfac=Rcontact(cr)%Lweight(ii,m)*ad_Zbot
                    ad_Zd(ii,1    ,m)=ad_Zd(ii,1    ,m)+adfac
                  END DO
                  ad_Zbot=0.0_r8
                END IF
              END DO
            END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.  No action required for the
!  adjoint of mp_assemble (AMM).
!
!>          CALL ad_mp_assemble (dg, model, Nzpts, spv, Zd)
!>
!!          CALL ad_mp_assemble (dg, model, Nzpts, spv, ad_Zd)
!!          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
            DO Kdg=1,N(dg)
              DO m=1,Npoints
                Idg  =Rcontact(cr)%Idg(m)
                Idgp1=MIN(Idg+1, BOUNDS(dg)%UBi(-1))
                Jdg  =Rcontact(cr)%Jdg(m)
                Jdgp1=MIN(Jdg+1, BOUNDS(dg)%UBj(-1))
                IF (((IminD.le.Idg).and.(Idg.le.ImaxD)).and.            &
     &              ((JminD.le.Jdg).and.(Jdg.le.JmaxD))) THEN
!>                tl_Zd(1,Kdg,m)=GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg)
!>
                  GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idg  ,Jdg  ,Kdg)+            &
     &                     ad_Zd(1,Kdg,m)
                  ad_Zd(1,Kdg,m)=0.0_r8
!>                tl_Zd(2,Kdg,m)=GRID(dg)%tl_z_r(Idgp1,Jdg  ,Kdg)
!>
                  GRID(dg)%ad_z_r(Idgp1,Jdg  ,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idgp1,Jdg  ,Kdg)+            &
     &                     ad_Zd(2,Kdg,m)
                  ad_Zd(2,Kdg,m)=0.0_r8
!>                tl_Zd(3,Kdg,m)=GRID(dg)%tl_z_r(Idgp1,Jdgp1,Kdg)
!>
                  GRID(dg)%ad_z_r(Idgp1,Jdgp1,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idgp1,Jdgp1,Kdg)+            &
     &                     ad_Zd(3,Kdg,m)
                  ad_Zd(3,Kdg,m)=0.0_r8
!>                tl_Zd(4,Kdg,m)=GRID(dg)%tl_z_r(Idg  ,Jdgp1,Kdg)
!>
                  GRID(dg)%ad_z_r(Idg  ,Jdgp1,Kdg)=                     &
     &                     GRID(dg)%ad_z_r(Idg  ,Jdgp1,Kdg)+            &
     &                     ad_Zd(4,Kdg,m)
                  ad_Zd(4,Kdg,m)=0.0_r8
                END IF
              END DO
            END DO
          END IF R_CONTACT

#  ifdef DISTRIBUTE
!
!  If distributed-memory, initialize with special value (zero) to
!  facilitate the global reduction when collecting data between all
!  nodes.
!
!>        Rcontact(cr)%tl_Vweight(1:2,1:N(rg),1:Npoints)=0.0_r8
!>
          Rcontact(cr)%ad_Vweight(1:2,1:N(rg),1:Npoints)=0.0_r8
#  endif
!
!  Deallocate local work arrays.
!
          IF (allocated(Zd)) THEN
            deallocate (Zd)
          END IF
          IF (allocated(ad_Zd)) THEN
            deallocate (ad_Zd)
          END IF

        END IF

      END DO

      RETURN

      END SUBROUTINE ad_z_weights
# endif

# ifdef SOLVE3D
!
      SUBROUTINE ad_put_contact3d (rg, model, tile,                     &
     &                             gtype, svname,                       &
     &                             cr, Npoints, contact,                &
     &                             LBi, UBi, LBj, UBj, LBk, UBk,        &
#  ifdef MASKING
     &                             Amask,                               &
#  endif
     &                             Ac, ad_Ac, ad_Ar)
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
      TYPE (T_NGC), intent(inout) :: contact(:)
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Ac(:,:,:)
      real(r8), intent(inout) :: ad_Ac(:,:,:)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: ad_Ar(LBi:,LBj:,LBk:)
#  else
      real(r8), intent(in) :: Ac(Npoints,LBk:UBk,4)
      real(r8), intent(inout) :: ad_Ac(Npoints,LBk:UBk,4)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: ad_Ar(LBi:UBi,LBj:UBj,LBk:UBk)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k, kdg, kdgm1, m, ii
      integer :: Istr, Iend, Jstr, Jend, Kmin

      real(r8), dimension(8) :: cff
      real(r8), dimension(8) :: ad_cff
!
!     Clear adjoint constants.
!
      DO ii=1,8
        ad_cff(ii)=0.0_r8
      END DO
!
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
!
#  ifdef MASKING
!>          tl_Ar(i,j,k)=tl_Ar(i,j,k)*Amask(i,j)
!>
            ad_Ar(i,j,k)=ad_Ar(i,j,k)*Amask(i,j)
#  endif
!>          tl_Ar(i,j,k)=tl_cff(1)*Ac(1,kdgm1,m)+                       &
!>   &                   cff(1)*tl_Ac(1,kdgm1,m)+                       &
!>   &                   tl_cff(2)*Ac(2,kdgm1,m)+                       &
!>   &                   cff(2)*tl_Ac(2,kdgm1,m)+                       &
!>   &                   tl_cff(3)*Ac(3,kdgm1,m)+                       &
!>   &                   cff(3)*tl_Ac(3,kdgm1,m)+                       &
!>   &                   tl_cff(4)*Ac(4,kdgm1,m)+                       &
!>   &                   cff(4)*tl_Ac(4,kdgm1,m)+                       &
!>   &                   tl_cff(5)*Ac(1,kdg  ,m)+                       &
!>   &                   cff(5)*tl_Ac(1,kdg  ,m)+                       &
!>   &                   tl_cff(6)*Ac(2,kdg  ,m)+                       &
!>   &                   cff(6)*tl_Ac(2,kdg  ,m)+                       &
!>   &                   tl_cff(7)*Ac(3,kdg  ,m)+                       &
!>   &                   cff(7)*tl_Ac(3,kdg  ,m)+                       &
!>   &                   tl_cff(8)*Ac(4,kdg  ,m)+                       &
!>   &                   cff(8)*tl_Ac(4,kdg  ,m)
!>
            ad_cff(1)=ad_cff(1)+Ac(1,kdgm1,m)*ad_Ar(i,j,k)
            ad_cff(2)=ad_cff(2)+Ac(2,kdgm1,m)*ad_Ar(i,j,k)
            ad_cff(3)=ad_cff(3)+Ac(3,kdgm1,m)*ad_Ar(i,j,k)
            ad_cff(4)=ad_cff(4)+Ac(4,kdgm1,m)*ad_Ar(i,j,k)
            ad_cff(5)=ad_cff(5)+Ac(1,kdg  ,m)*ad_Ar(i,j,k)
            ad_cff(6)=ad_cff(6)+Ac(2,kdg  ,m)*ad_Ar(i,j,k)
            ad_cff(7)=ad_cff(7)+Ac(3,kdg  ,m)*ad_Ar(i,j,k)
            ad_cff(8)=ad_cff(8)+Ac(4,kdg  ,m)*ad_Ar(i,j,k)

            ad_Ac(1,kdgm1,m)=ad_Ac(1,kdgm1,m)+cff(1)*ad_Ar(i,j,k)
            ad_Ac(2,kdgm1,m)=ad_Ac(2,kdgm1,m)+cff(2)*ad_Ar(i,j,k)
            ad_Ac(3,kdgm1,m)=ad_Ac(3,kdgm1,m)+cff(3)*ad_Ar(i,j,k)
            ad_Ac(4,kdgm1,m)=ad_Ac(4,kdgm1,m)+cff(4)*ad_Ar(i,j,k)
            ad_Ac(1,kdg  ,m)=ad_Ac(1,kdg  ,m)+cff(5)*ad_Ar(i,j,k)
            ad_Ac(2,kdg  ,m)=ad_Ac(2,kdg  ,m)+cff(6)*ad_Ar(i,j,k)
            ad_Ac(3,kdg  ,m)=ad_Ac(3,kdg  ,m)+cff(7)*ad_Ar(i,j,k)
            ad_Ac(4,kdg  ,m)=ad_Ac(4,kdg  ,m)+cff(8)*ad_Ar(i,j,k)

            ad_Ar(i,j,k)=0.0_r8

!>          tl_cff(1)=contact(cr)%Lweight(1,m)*                         &
!>   &                contact(cr)%tl_Vweight(1,k,m)
!>          tl_cff(2)=contact(cr)%Lweight(2,m)*                         &
!>   &                contact(cr)%tl_Vweight(1,k,m)
!>          tl_cff(3)=contact(cr)%Lweight(3,m)*                         &
!>   &                contact(cr)%tl_Vweight(1,k,m)
!>          tl_cff(4)=contact(cr)%Lweight(4,m)*                         &
!>   &                contact(cr)%tl_Vweight(1,k,m)
!>          tl_cff(5)=contact(cr)%Lweight(1,m)*                         &
!>   &                contact(cr)%tl_Vweight(2,k,m)
!>          tl_cff(6)=contact(cr)%Lweight(2,m)*                         &
!>   &                contact(cr)%tl_Vweight(2,k,m)
!>          tl_cff(7)=contact(cr)%Lweight(3,m)*                         &
!>   &                contact(cr)%tl_Vweight(2,k,m)
!>          tl_cff(8)=contact(cr)%Lweight(4,m)*                         &
!>   &                contact(cr)%tl_Vweight(2,k,m)
!>
            contact(cr)%ad_Vweight(1,k,m)=                              &
     &                  contact(cr)%ad_Vweight(1,k,m)+                  &
     &                  contact(cr)%Lweight(1,m)*ad_cff(1)+             &
     &                  contact(cr)%Lweight(2,m)*ad_cff(2)+             &
     &                  contact(cr)%Lweight(3,m)*ad_cff(3)+             &
     &                  contact(cr)%Lweight(4,m)*ad_cff(4)
            ad_cff(1)=0.0_r8
            ad_cff(2)=0.0_r8
            ad_cff(3)=0.0_r8
            ad_cff(4)=0.0_r8

            contact(cr)%ad_Vweight(2,k,m)=                              &
                        contact(cr)%ad_Vweight(2,k,m)+                  &
     &                  contact(cr)%Lweight(1,m)*ad_cff(5)+             &
     &                  contact(cr)%Lweight(2,m)*ad_cff(6)+             &
     &                  contact(cr)%Lweight(3,m)*ad_cff(7)+             &
     &                  contact(cr)%Lweight(4,m)*ad_cff(8)
            ad_cff(5)=0.0_r8
            ad_cff(6)=0.0_r8
            ad_cff(7)=0.0_r8
            ad_cff(8)=0.0_r8
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE ad_put_contact3d
# endif
!
      SUBROUTINE ad_put_contact2d (rg, model, tile,                     &
     &                             gtype, svname,                       &
     &                             cr, Npoints, contact,                &
     &                             LBi, UBi, LBj, UBj,                  &
# ifdef MASKING
     &                             Amask,                               &
# endif
     &                             Ac, Ar)
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
      real(r8), intent(inout) :: Ac(:,:)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: Ar(LBi:,LBj:)
# else
      real(r8), intent(inout) :: Ac(4,Npoints)
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: Ar(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, m, ii
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
# ifdef MASKING
          Ar(i,j)=Ar(i,j)*Amask(i,j)
# endif
!>        Ar(i,j)=contact(cr)%Lweight(1,m)*Ac(1,m)+                     &
!>   &            contact(cr)%Lweight(2,m)*Ac(2,m)+                     &
!>   &            contact(cr)%Lweight(3,m)*Ac(3,m)+                     &
!>   &            contact(cr)%Lweight(4,m)*Ac(4,m)
!>
          DO ii=1,4
            Ac(ii,m)=Ac(ii,m)+contact(cr)%Lweight(ii,m)*Ar(i,j)
          END DO
          Ar(i,j)=0.0_r8
        END IF
      END DO

      RETURN
      END SUBROUTINE ad_put_contact2d
!
      SUBROUTINE ad_fine2coarse2d (ng, dg, model, tile,                 &
     &                             gtype, svname,                       &
     &                             AreaAvg, Rscale,                     &
     &                             cr, Npoints, contact,                &
     &                             LBiF, UBiF, LBjF, UBjF,              &
     &                             LBiC, UBiC, LBjC, UBjC,              &
# ifdef DISTRIBUTE
     &                             Adx, Ady,                            &
# else
     &                             dxF, dyF,                            &
# endif
     &                             pmC, pnC,                            &
# ifdef MASKING
#  ifdef DISTRIBUTE
     &                             Amsk,                                &
#  else
     &                             Fmsk,                                &
#  endif
     &                             Cmsk,                                &
# endif
# ifdef DISTRIBUTE
     &                             A,                                   &
# else
     &                             F,                                   &
# endif
     &                             C1, C2)
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
      USE strings_mod,    ONLY : FoundError
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
      real(r8), intent(inout) :: A(LBiF:,LBjF:)
      real(r8), intent(in) :: Adx(LBiF:,LBjF:)
      real(r8), intent(in) :: Ady(LBiF:,LBjF:)
#  else
      real(r8), intent(inout) :: F(LBiF:,LBjF:)
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
      real(r8), intent(inout) :: A(LBiF:UBiF,LBjF:UBjF)
      real(r8), intent(in) :: Adx(LBiF:UBiF,LBjF:UBjF)
      real(r8), intent(in) :: Adx(LBiF:UBiF,LBjF:UBjF)
#  else
      real(r8), intent(inout) :: F(LBiF:UBiF,LBjF:UBjF)
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
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_fine2coarse2d"

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Average interior fine grid state variable data to the coarse grid
!  location. Then, replace coarse grid values with averaged data.
!-----------------------------------------------------------------------
!
!  Clear constants.
!
      my_avg=0.0_r8
      my_sum=0.0_r8
      my_count=0.0_r8
      my_area=0.0_r8
      my_areasum=0.0_r8

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
      IF (AreaAvg) THEN
        CALL mp_aggregate2d (dg, model, gtype,                          &
     &                       LBiF, UBiF, LBjF, UBjF,                    &
     &                       LBi,  UBi,  LBj,  UBj,                     &
     &                       Adx, dxF)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL mp_aggregate2d (dg, model, gtype,                          &
     &                       LBiF, UBiF, LBjF, UBjF,                    &
     &                       LBi,  UBi,  LBj,  UBj,                     &
     &                       Ady, dyF)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  ifdef MASKING
!
      CALL mp_aggregate2d (dg, model, gtype,                            &
     &                     LBiF, UBiF, LBjF, UBjF,                      &
     &                     LBi,  UBi,  LBj,  UBj,                       &
     &                     Amsk, Fmsk)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
# ifdef MASKING
            DO Jadd=-half,half
              DO Iadd=-half,half
                my_count=my_count+MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
              END DO
            END DO
# endif
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
            IF (PRESENT(C2)) THEN
!>            C2(Ic,Jc)=my_avg
              my_avg=my_avg+C2(Ic,Jc)
              C2(Ic,Jc)=0.0_r8
            END IF
!>          C1(Ic,Jc)=my_avg
            my_avg=my_avg+C1(Ic,Jc)
            C1(Ic,Jc)=0.0_r8
# ifdef MASKING
            my_avg=my_avg*Cmsk(Ic,Jc)
            IF (my_count.gt.0.0_r8) THEN
              my_avg=my_avg*Rscale*Rscale/my_count
            END IF
# endif
!>          my_avg=my_sum*areaC_inv

            my_sum=my_sum+areaC_inv*my_avg
            my_avg=0.0_r8

!!          ratio=my_areasum*areaC_inv       ! for debugging purposes


            DO Jadd=-half,half
              DO Iadd=-half,half
                my_area=dxF(i+Iadd,j+Jadd)*dyF(i+Iadd,j+Jadd)
                my_areasum=my_areasum+my_area
# ifdef MASKING
!>              my_sum=my_sum+                                          &
!>   &                 F(i+Iadd,j+Jadd)*my_area*                        &
!>   &                 MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
                F(i+Iadd,j+Jadd)=F(i+Iadd,j+Jadd)+my_area*              &
     &                 MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))*my_sum
# else
!>              my_sum=my_sum+                                          &
!>   &                 F(i+Iadd,j+Jadd)*my_area
                F(i+Iadd,j+Jadd)=F(i+Iadd,j+Jadd)+my_area*my_sum
# endif
              END DO
            END DO
            my_sum=0.0_r8
            my_areasum=0.0_r8
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
!
!  Compute my_count first.
!
            my_count=0.0_r8
            DO Jadd=-half,half
              DO Iadd=-half,half
# ifdef MASKING
                my_count=my_count+MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
# else
                my_count=my_count+1.0_r8
# endif
              END DO
            END DO
            IF (PRESENT(C2)) THEN
!>            C2(Ic,Jc)=my_avg
              my_avg=my_avg+C2(Ic,Jc)
              C2(Ic,Jc)=0.0_r8
            END IF
            my_avg=my_avg+C1(Ic,Jc)
            C1(Ic,Jc)=0.0_r8
# ifdef MASKING
            my_avg=my_avg*Cmsk(Ic,Jc)
# endif
            IF (my_count.gt.0.0_r8) THEN
!>             my_avg=my_sum/my_count
               my_sum=my_sum+my_avg/my_count
               my_avg=0.0_r8
            END IF
            DO Jadd=-half,half
              DO Iadd=-half,half
# ifdef MASKING
!>              my_sum=my_sum+                                          &
!>   &                 F(i+Iadd,j+Jadd)*Fmsk(i+Iadd,j+Jadd)
                F(i+Iadd,j+Jadd)=F(i+Iadd,j+Jadd)+Fmsk(i+Iadd,j+Jadd)*  &
     &                           my_sum
# else
!>              my_sum=my_sum+                                          &
!>   &                 F(i+Iadd,j+Jadd)
                F(i+Iadd,j+Jadd)=F(i+Iadd,j+Jadd)+my_sum
# endif
              END DO
            END DO
            my_sum=0.0_r8
          END IF
        END DO
      END IF

# ifdef DISTRIBUTE

!AMM
!
!  This next loop represents the adjoint of mp_aggregate2d.
!
      DO j=LBjF,UBjF
        DO i=LBiF,UBiF
          A(i,j)=A(i,j)+F(i,j)
          F(i,j)=0.0_r8
        END DO
      END DO
!
!AMM
!

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
      END SUBROUTINE ad_fine2coarse2d
!
# ifdef SOLVE3D
      SUBROUTINE ad_fine2coarse3d (ng, dg, model, tile,                 &
     &                             gtype, svname,                       &
     &                             AreaAvg, Rscale,                     &
     &                             cr, Npoints, contact,                &
     &                             LBiF, UBiF, LBjF, UBjF, LBkF, UBkF,  &
     &                             LBiC, UBiC, LBjC, UBjC, LBkC, UBkC,  &
#  ifdef DISTRIBUTE
     &                             Adx, Ady,                            &
#  else
     &                             dxF, dyF,                            &
#  endif
     &                             pmC, pnC,                            &
#  ifdef MASKING
#   ifdef DISTRIBUTE
     &                             Amsk,                                &
#   else
     &                             Fmsk,                                &
#   endif
     &                             Cmsk,                                &
#  endif
#  ifdef DISTRIBUTE
     &                             A,                                   &
#  else
     &                             F,                                   &
#  endif
     &                             C)
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
!
#  ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_aggregate2d
      USE distribute_mod, ONLY : mp_aggregate3d
#  endif
      USE strings_mod,    ONLY : FoundError
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
      real(r8), intent(inout) :: A(LBiF:,LBjF:,LBkF:)
      real(r8), intent(in) :: Adx(LBiF:,LBjF:)
      real(r8), intent(in) :: Ady(LBiF:,LBjF:)
#   else
      real(r8), intent(inout) :: F(LBiF:,LBjF:,LBkF:)
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
      real(r8), intent(inout) :: A(LBiF:UBiF,LBjF:UBjF,LBkF:UBkF)
      real(r8), intent(in) :: Adx(LBiF:UBiF,LBjF:UBjF)
      real(r8), intent(in) :: Adx(LBiF:UBiF,LBjF:UBjF)
#   else
      real(r8), intent(inout) :: F(LBiF:UBiF,LBjF:UBjF,LBkF:UBkF)
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
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_fine2coarse3d"

#  include "set_bounds.h"
!
!  Clear constants.
!
      my_area=0.0_r8
      my_areasum=0.0_r8
      my_avg=0.0_r8
      my_count=0.0_r8
      my_sum=0.0_r8
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
      IF (AreaAvg) THEN
        CALL mp_aggregate2d (dg, model, gtype,                          &
     &                       LBiF, UBiF, LBjF, UBjF,                    &
     &                       LBi,  UBi,  LBj,  UBj,                     &
     &                       Adx, dxF)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL mp_aggregate2d (dg, model, gtype,                          &
     &                       LBiF, UBiF, LBjF, UBjF,                    &
     &                       LBi,  UBi,  LBj,  UBj,                     &
     &                       Ady, dyF)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   ifdef MASKING
!
      CALL mp_aggregate2d (dg, model, gtype,                            &
     &                     LBiF, UBiF, LBjF, UBjF,                      &
     &                     LBi,  UBi,  LBj,  UBj,                       &
     &                     Amsk, Fmsk)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!
!  Compute my_count first.
!
              my_count=0.0_r8
#  ifdef MASKING
              DO Jadd=-half,half
                DO Iadd=-half,half
                  my_count=my_count+MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
                END DO
              END DO
#  endif
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
!
!>            C(Ic,Jc,k)=my_avg
!>
              my_avg=my_avg+C(Ic,Jc,k)
              C(Ic,Jc,k)=0.0_r8
#  ifdef MASKING
              my_avg=my_avg*Cmsk(Ic,Jc)
              IF (my_count.gt.0.0_r8) THEN
                my_avg=my_avg*Rscale*Rscale/my_count
              END IF
#  endif
!>            my_avg=my_sum*areaC_inv

              my_sum=my_sum+areaC_inv*my_avg
              my_avg=0.0_r8

!>            ratio=my_areasum*areaC_inv     ! for debugging purposes
!>
              DO Jadd=-half,half
                DO Iadd=-half,half
                  my_area=dxF(i+Iadd,j+Jadd)*dyF(i+Iadd,j+Jadd)
                  my_areasum=my_areasum+my_area
#  ifdef MASKING
!>                my_sum=my_sum+                                        &
!>   &                   F(i+Iadd,j+Jadd,k)*my_area*                    &
!>   &                   MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
!>
                  F(i+Iadd,j+Jadd,k)=F(i+Iadd,j+Jadd,k)+                &
     &                               my_area*                           &
     &                               MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))*   &
     &                               my_sum
#  else
!>                my_sum=my_sum+                                        &
!>   &                   F(i+Iadd,j+Jadd,k)*my_area
!>
                  F(i+Iadd,j+Jadd,k)=F(i+Iadd,j+Jadd,k)+                &
     &                               my_area*my_sum
#  endif
                END DO
              END DO
              my_count=0.0_r8
              my_sum=0.0_r8
              my_areasum=0.0_r8
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
!
!  Compute my_count first.
!
              my_count=0.0_r8
              DO Jadd=-half,half
                DO Iadd=-half,half
#  ifdef MASKING
                  my_count=my_count+MIN(1.0_r8,Fmsk(i+Iadd,j+Jadd))
#  else
                  my_count=my_count+1.0_r8
#  endif
                END DO
              END DO
!
!>            C(Ic,Jc,k)=my_avg
!>
              my_avg=my_avg+C(Ic,Jc,k)
              C(Ic,Jc,k)=0.0_r8
#  ifdef MASKING
              my_avg=my_avg*Cmsk(Ic,Jc)
#  endif
              IF (my_count.gt.0.0_r8) THEN
!>               my_avg=my_sum/my_count
!>
                 my_sum=my_sum+my_avg/my_count
                 my_avg=0.0_r8
              END IF

              DO Jadd=-half,half
                DO Iadd=-half,half
#  ifdef MASKING
!>                my_sum=my_sum+                                        &
!>   &                   F(i+Iadd,j+Jadd,k)*Fmsk(i+Iadd,j+Jadd)
!>
                  F(i+Iadd,j+Jadd,k)=F(i+Iadd,j+Jadd,k)+                &
     &                               Fmsk(i+Iadd,j+Jadd)*my_sum
#  else
!>                my_sum=my_sum+                                        &
!>   &                   F(i+Iadd,j+Jadd,k)
!>
                  F(i+Iadd,j+Jadd,k)=F(i+Iadd,j+Jadd,k)+my_sum
#  endif
                END DO
              END DO
              my_count=0.0_r8
              my_avg=0.0_r8
              my_sum=0.0_r8
            END IF
          END DO
        END DO
      END IF

#  ifdef DISTRIBUTE
!
!  The following loop represents the adjoint of mp_aggregate3d (AMM).
!
      DO k=LBkF,UBkF
        DO j=LBjF,UBjF
          DO i=LBiF,UBiF
             A(i,j,k)=A(i,j,k)+F(i,j,k)
             F(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
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
      END SUBROUTINE ad_fine2coarse3d
# endif
!
      SUBROUTINE ad_get_contact2d (dg, model, tile,                     &
     &                             gtype, svname,                       &
     &                             cr, Npoints, contact,                &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             Ad, Ac)
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
!      USE distribute_mod, ONLY : ad_mp_assemble
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
      real(r8), intent(inout) :: Ad(LBi:,LBj:)
      real(r8), intent(inout) :: Ac(:,:)
# else
      real(r8), intent(inout) :: Ad(LBi:UBi,LBj:UBj)
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
!  Adjoint of extract donor grid data at contact points.
!-----------------------------------------------------------------------

# ifdef DISTRIBUTE
!
!  Gather and broadcast data from all nodes.  No action required for the
!  adjoint of mp_assemble (AMM).
!
!!    Npts=4*Npoints
!>    CALL mp_assemble (dg, model, Npts, Aspv, Ac)
!>
!!    CALL ad_mp_assemble (dg, model, Npts, Aspv, Ac)
# endif
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
!>        Ac(1,m)=Ad(i  ,j  )
!>
          Ad(i  ,j  )=Ad(i  ,j  )+Ac(1,m)
          Ac(1,m)=0.0_r8
!>        Ac(2,m)=Ad(ip1,j  )
!>
          Ad(ip1,j  )=Ad(ip1,j  )+Ac(2,m)
          Ac(2,m)=0.0_r8
!>        Ac(3,m)=Ad(ip1,jp1)
!>
          Ad(ip1,jp1)=Ad(ip1,jp1)+Ac(3,m)
          Ac(3,m)=0.0_r8
!>        Ac(4,m)=Ad(i  ,jp1)
!>
          Ad(i  ,jp1)=Ad(i  ,jp1)+Ac(4,m)
          Ac(4,m)=0.0_r8
        END IF
      END DO

# ifdef DISTRIBUTE
!
!  Adjoint of initialize contact points array to special value to
!  facilitate distribute-memory data collection from all nodes.
!
      DO m=1,Npoints
        Ac(1,m)=0.0_r8
        Ac(2,m)=0.0_r8
        Ac(3,m)=0.0_r8
        Ac(4,m)=0.0_r8
      END DO
# endif

      RETURN
      END SUBROUTINE ad_get_contact2d

# ifdef SOLVE3D
!
      SUBROUTINE ad_get_contact3d (dg, model, tile,                     &
     &                             gtype, svname,                       &
     &                             cr, Npoints, contact,                &
     &                             LBi, UBi, LBj, UBj, LBk, UBk,        &
     &                             Ad, Ac)
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
!
#  ifdef DISTRIBUTE
!!    USE distribute_mod, ONLY : ad_mp_assemble
#  endif
      USE strings_mod,    ONLY : FoundError
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
      real(r8), intent(inout) :: Ad(LBi:,LBj:,LBk:)
      real(r8), intent(inout) :: Ac(:,LBk:,:)
#  else
      real(r8), intent(inout) :: Ad(LBi:UBi,LBj:UBj,LBk:UBk)
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
!
      real(r8), parameter :: Aspv = 0.0_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_get_contact3d"
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
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
!  Adjoint of extract donor grid data at contact points.
!-----------------------------------------------------------------------

#  ifdef DISTRIBUTE
!
!  Gather and broadcast data from all nodes.  No action required for
!  the adjoint of mp_assemble.
!
!!    Npts=4*(UBk-LBk+1)*Npoints
!>    CALL mp_assemble (dg, model, Npts, Aspv, Ac(:,LBk:,:))
!>
!!    CALL ad_mp_assemble (dg, model, Npts, Aspv, Ac(:,LBk:,:))
!!    IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
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
!>          Ac(1,k,m)=Ad(i  ,j  ,k)
!>
            Ad(i  ,j  ,k)=Ad(i  ,j  ,k)+Ac(1,k,m)
            Ac(1,k,m)=0.0_r8
!>          Ac(2,k,m)=Ad(ip1,j  ,k)
!>
            Ad(ip1,j  ,k)=Ad(ip1,j  ,k)+Ac(2,k,m)
            Ac(2,k,m)=0.0_r8
!>          Ac(3,k,m)=Ad(ip1,jp1,k)
!>
            Ad(ip1,jp1,k)=Ad(ip1,jp1,k)+Ac(3,k,m)
            Ac(3,k,m)=0.0_r8
!>          Ac(4,k,m)=Ad(i  ,jp1,k)
!>
            Ad(i  ,jp1,k)=Ad(i  ,jp1,k)+Ac(4,k,m)
            Ac(4,k,m)=0.0_r8
          END IF
        END DO
      END DO

#  ifdef DISTRIBUTE
!
!  Adjoint of initialize contact points array to special value to
!  facilitate distribute-memory data collection from all nodes.
!
      DO k=LBk,UBk
        DO m=1,Npoints
          Ac(1,k,m)=0.0_r8
          Ac(2,k,m)=0.0_r8
          Ac(3,k,m)=0.0_r8
          Ac(4,k,m)=0.0_r8
        END DO
      END DO
#  endif
!
      RETURN
      END SUBROUTINE ad_get_contact3d
# endif
!
      SUBROUTINE ad_get_persisted2d (dg, rg, model, tile,               &
     &                               gtype, svname,                     &
     &                               cr, Npoints, contact,              &
     &                               LBi, UBi, LBj, UBj,                &
     &                               Ad, Ac)
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
!
# ifdef DISTRIBUTE
!!    USE distribute_mod, ONLY : ad_mp_assemble
# endif
      USE strings_mod,    ONLY : FoundError
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
      real(r8), intent(inout) :: Ad(LBi:,LBj:)
      real(r8), intent(inout) :: Ac(:,:)
# else
      real(r8), intent(inout) :: Ad(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Ac(Npoints,4)
# endif
!
!  Local variable declarations.
!
      integer :: Idg, Ip1, Irg, Jdg, Jp1, Jrg, ii
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Istr, Iend, Jstr, Jend
      integer :: i, i_add, j, j_add, m, m_add
# ifdef DISTRIBUTE
      integer :: Npts
# endif
!
      real(r8), parameter :: Aspv = 0.0_r8
      real(r8):: Rscale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_get_persisted2d"
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
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
!  Adjoint of extract donor grid data at contact points.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!
!  Gather and broadcast data from all nodes. No action required for the
!  adjoint of mp_assemble (AMM).
!
!>    CALL mp_assemble (dg, model, Npts, Aspv, Ac)
!>
!!    CALL ad_mp_assemble (dg, model, Npts, Aspv, Ac)
!!    IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
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
!>            Ac(1,m)=Ad(Idg,j)
!>            Ac(2,m)=Ad(Idg,j)
!>            Ac(3,m)=Ad(Idg,j)
!>            Ac(4,m)=Ad(Idg,j)
!>
              DO ii=1,4
                Ad(Idg,j)=Ad(Idg,j)+Ac(ii,m)
                Ac(ii,m)=0.0_r8
              END DO
            ELSE IF ((on_boundary(m+m_add).eq.2).or.                    &
     &               (on_boundary(m+m_add).eq.4)) THEN  ! southern and
              i_add=INT(REAL(Irg-1,r8)*Rscale)          ! northern edges
              i=I_left(rg)+i_add
!>            Ac(1,m)=Ad(i,Jdg)
!>            Ac(2,m)=Ad(i,Jdg)
!>            Ac(3,m)=Ad(i,Jdg)
!>            Ac(4,m)=Ad(i,Jdg)
!>
              DO ii=1,4
                Ad(i,Jdg)=Ad(i,Jdg)+Ac(ii,m)
                Ac(ii,m)=0.0_r8
              END DO
            END IF
!
!  Contact point is not at physical boundary, just set values for spatial
!  interpolation (not used).
!
          ELSE
!>          Ac(1,m)=Ad(Idg,Jdg)
!>          Ac(2,m)=Ad(Ip1,Jdg)
!>          Ac(3,m)=Ad(Ip1,Jp1)
!>          Ac(4,m)=Ad(Idg,Jp1)
!>
            Ad(Idg,Jdg)=Ad(Idg,Jdg)+Ac(1,m)
            Ac(1,m)=0.0_r8
            Ad(Ip1,Jdg)=Ad(Ip1,Jdg)+Ac(2,m)
            Ac(2,m)=0.0_r8
            Ad(Ip1,Jp1)=Ad(Ip1,Jp1)+Ac(3,m)
            Ac(3,m)=0.0_r8
            Ad(Idg,Jp1)=Ad(Idg,Jp1)+Ac(4,m)
            Ac(4,m)=0.0_r8
          END IF
        END IF
      END DO

# ifdef DISTRIBUTE
!
!  Initialize contact points array to special value to facilite
!  distribute-memory data collection from all nodes.
!
      DO m=1,Npoints
        Ac(1,m)=0.0_r8
        Ac(2,m)=0.0_r8
        Ac(3,m)=0.0_r8
        Ac(4,m)=0.0_r8
      END DO
# endif
!
      RETURN
      END SUBROUTINE ad_get_persisted2d
!
      SUBROUTINE ad_bry_fluxes (dg, rg, cr, model, tile,                &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          ILB, IUB, JLB, JUB,                     &
     &                          scale, FX, FE,                          &
     &                          F_west, F_east, F_south, F_north)
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
!
# ifdef DISTRIBUTE
!!    USE distribute_mod, ONLY : ad_mp_assemble
# endif
      USE strings_mod,    ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: dg, rg, cr, model, tile
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: ILB, IUB, JLB, JUB

      real(r8), intent(in) :: scale
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: FX(IminS:,JminS:)
      real(r8), intent(inout) :: FE(IminS:,JminS:)

      real(r8), intent(inout) :: F_west (JLB:)
      real(r8), intent(inout) :: F_east (JLB:)
      real(r8), intent(inout) :: F_south(ILB:)
      real(r8), intent(inout) :: F_north(ILB:)
# else
      real(r8), intent(inout) :: FX(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(inout) :: FE(IminS:ImaxS,JminS:JmaxS)

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
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_bry_fluxes"

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Gather and broadcast data from all nodes.
!-----------------------------------------------------------------------
!
!  No action required for the adjoint of mp_assemble (AMM).
!
!>    CALL mp_assemble (dg, model, NptsWE, Fspv, F_west (JLB:))
!>
!!    CALL ad_mp_assemble (dg, model, NptsWE, Fspv, F_west (JLB:))
!!    IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

!>    CALL mp_assemble (dg, model, NptsWE, Fspv, F_east (JLB:))
!>
!!    CALL ad_mp_assemble (dg, model, NptsWE, Fspv, F_east (JLB:))
!!    IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

!>    CALL mp_assemble (dg, model, NptsSN, Fspv, F_south(ILB:))
!>
!!    CALL ad_mp_assemble (dg, model, NptsSN, Fspv, F_south(ILB:))
!!    IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

!>    CALL mp_assemble (dg, model, NptsSN, Fspv, F_north(ILB:))
!>
!!    CALL ad_mp_assemble (dg, model, NptsSN, Fspv, F_north(ILB:))
!!    IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!  Northern boundary.
!
        IF (DOMAIN(dg)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
!>          F_north(i)=FE(i,Jend+1)*scale
!>
            FE(i,Jend+1)=FE(i,Jend+1)+scale*F_north(i)
            F_north(i)=0.0_r8
          END DO
        END IF
!
!  Southern boundary.
!
        IF (DOMAIN(dg)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
!>          F_south(i)=FE(i,Jstr)*scale
!>
            FE(i,Jstr)=FE(i,Jstr)+scale*F_south(i)
            F_south(i)=0.0_r8
          END DO
        END IF
!
!  Eastern boundary.
!
        IF (DOMAIN(dg)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
!>          F_east(j)=FX(Iend+1,j)*scale
!>
            FX(Iend+1,j)=FX(Iend+1,j)+scale*F_east(j)
            F_east(j)=0.0_r8
          END DO
        END IF
!
!  Western boundary.
!
        IF (DOMAIN(dg)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
!>          F_west(j)=FX(Istr,j)*scale
!>
            FX(Istr,j)=FX(Istr,j)+scale*F_west(j)
            F_west(j)=0.0_r8
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
!  Southern/Northern boundaries.
!
        Jb_south=J_bottom(dg)
        Jb_north=J_top(dg)
        DO i=Istr,Iend
          IF ((Jstr.le.Jb_south).and.(Jb_south.le.Jend)) THEN
!>          F_south(i)=FE(i,Jb_south)*scale
!>
            FE(i,Jb_south)=FE(i,Jb_south)+scale*F_south(i)
            F_south(i)=0.0_r8
          END IF
!
          IF ((Jstr.le.Jb_north).and.(Jb_north.le.Jend)) THEN
!>          F_north(i)=FE(i,Jb_north)*scale
!>
            FE(i,Jb_north)=FE(i,Jb_north)+scale*F_north(i)
            F_north(i)=0.0_r8
          END IF
        END DO
!
!  Western/Eastern boundaries.
!
        Ib_west=I_left(dg)
        Ib_east=I_right(dg)
        DO j=Jstr,Jend
          IF ((Istr.le.Ib_west).and.(Ib_west.le.Iend)) THEN
!>          F_west(j)=FX(Ib_west,j)*scale
!>
            FX(Ib_west,j)=FX(Ib_west,j)+scale*F_west(j)
            F_west(j)=0.0_r8
          END IF
!
          IF ((Istr.le.Ib_east).and.(Ib_east.le.Iend)) THEN
!>          F_east(j)=FX(Ib_east,j)*scale
!>
            FX(Ib_east,j)=FX(Ib_east,j)+scale*F_east(j)
            F_east(j)=0.0_r8
          END IF
        END DO
      END IF

# ifdef DISTRIBUTE
!
!  Initialize arrays to facilitate collective communications.
!
      NptsWE=JUB-JLB+1
      NptsSN=IUB-ILB+1
!
      F_west =0.0_r8
      F_east =0.0_r8
      F_south=0.0_r8
      F_north=0.0_r8
# endif

      RETURN
      END SUBROUTINE ad_bry_fluxes

# ifdef NESTING_DEBUG
!
      SUBROUTINE ad_check_massflux (ngf, model, tile)
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
      real(r8) :: ad_EastSum, ad_NorthSum, ad_SouthSum, ad_WestSum
# ifdef NESTING_DEBUG_NOT
      real(r8) :: MFratio
# endif
!
!  Clear adjoint constants.
!
      ad_EastSum=0.0_r8
      ad_NorthSum=0.0_r8
      ad_SouthSum=0.0_r8
      ad_WestSum=0.0_r8
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
# ifdef NESTING_DEBUG_NOT
          IF (DOMAIN(ngf)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (302,10) 'Western Boundary Mass Fluxes: ',          &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
              CALL my_flush (302)
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
              Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
              DO jsum=-half,half
                Jbf=Jedge+jsum
!>              tl_WestSum=tl_WestSum+                                  &
!>   &                     BRY_CONTACT(iwest,cr)%tl_Mflux(Jbf)
!>
                BRY_CONTACT(iwest,cr)%ad_Mflux(Jbf)=                    &
     &              BRY_CONTACT(iwest,cr)%ad_Mflux(Jbf)+ad_WestSum
              END DO
!>            tl_WestSum=0.0_r8
!>
              ad_WestSum=0.0_r8
              m=BRY_CONTACT(iwest,cr)%C2Bindex(Jbf)     ! pick last one
!
!  Load coarser grid western boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
!>            BRY_CONTACT(iwest,cjcr)%tl_Mflux(Jbc)=tl_WestSum
!>
              ad_WestSum=ad_WestSum+                                    &
     &                   BRY_CONTACT(iwest,cjcr)%ad_Mflux(Jbc)
              BRY_CONTACT(iwest,cjcr)%ad_Mflux(Jbc)=0.0_r8

# ifdef NESTING_DEBUG_NOT
              IF (WestSum.ne.0) THEN
                MFratio=REFINED(cr)%DU_avg2(1,m,tnew)/WestSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (302,30) Jbc, REFINED(cr)%DU_avg2(1,m,tnew),        &
     &                       WestSum, MFratio
              CALL my_flush (302)
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
# ifdef NESTING_DEBUG_NOT
          IF (DOMAIN(ngf)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (302,10) 'Eastern Boundary Mass Fluxes: ',          &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
              CALL my_flush (302)
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
              Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
              DO jsum=-half,half
                Jbf=Jedge+jsum
!>              tl_EastSum=tl_EastSum+                                  &
!>   &                     BRY_CONTACT(ieast,cr)%tl_Mflux(Jbf)
!>
                BRY_CONTACT(ieast,cr)%ad_Mflux(Jbf)=                    &
     &              BRY_CONTACT(ieast,cr)%ad_Mflux(Jbf)+ad_EastSum
              END DO
!>            tl_EastSum=0.0_r8
!>
              ad_EastSum=0.0_r8
              m=BRY_CONTACT(ieast,cr)%C2Bindex(Jbf)     ! pick last one
!
!  Load coarser grid eastern boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
!>            BRY_CONTACT(ieast,cjcr)%tl_Mflux(Jbc)=tl_EastSum
!>
              ad_EastSum=ad_EastSum+                                    &
     &                   BRY_CONTACT(ieast,cjcr)%ad_Mflux(Jbc)
              BRY_CONTACT(ieast,cjcr)%ad_Mflux(Jbc)=0.0_r8

# ifdef NESTING_DEBUG_NOT
              IF (EastSum.ne.0) THEN
                MFratio=REFINED(cr)%DU_avg2(1,m,tnew)/EastSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (302,30) Jbc, REFINED(cr)%DU_avg2(1,m,tnew),        &
     &                       EastSum, MFratio
              CALL my_flush (302)
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
# ifdef NESTING_DEBUG_NOT
          IF (DOMAIN(ngf)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (302,20) 'Southern Boundary Mass Fluxes: ',         &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
              CALL my_flush (302)
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
              Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
              DO isum=-half,half
                Ibf=Iedge+isum
!>              tl_SouthSum=tl_SouthSum+                                &
!>   &                      BRY_CONTACT(isouth,cr)%tl_Mflux(Ibf)
!>
                BRY_CONTACT(isouth,cr)%ad_Mflux(Ibf)=                   &
     &                BRY_CONTACT(isouth,cr)%ad_Mflux(Ibf)+ad_SouthSum
              END DO
!>            tl_SouthSum=0.0_r8
!>
              ad_SouthSum=0.0_r8
              m=BRY_CONTACT(isouth,cr)%C2Bindex(Ibf)    ! pick last one
!
!  Load coarser grid southern boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
!>            BRY_CONTACT(isouth,cjcr)%tl_Mflux(Ibc)=tl_SouthSum
!>
              ad_SouthSum=ad_SouthSum+                                  &
     &                    BRY_CONTACT(isouth,cjcr)%ad_Mflux(Ibc)
              BRY_CONTACT(isouth,cjcr)%ad_Mflux(Ibc)=0.0_r8

# ifdef NESTING_DEBUG_NOT
              IF (SouthSum.ne.0) THEN
                MFratio=REFINED(cr)%DV_avg2(1,m,tnew)/SouthSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (302,30) Ibc, REFINED(cr)%DV_avg2(1,m,tnew),        &
     &                       SouthSum, MFratio
              CALL my_flush (302)
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
# ifdef NESTING_DEBUG_NOT
          IF (DOMAIN(ngf)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              WRITE (302,20) 'Northern Boundary Mass Fluxes: ',         &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
              CALL my_flush (302)
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
              Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
              DO isum=-half,half
                Ibf=Iedge+isum
!>              tl_NorthSum=tl_NorthSum+                                &
!>   &                      BRY_CONTACT(inorth,cr)%tl_Mflux(Ibf)
!>
                BRY_CONTACT(inorth,cr)%ad_Mflux(Ibf)=                   &
     &              BRY_CONTACT(inorth,cr)%ad_Mflux(Ibf)+ad_NorthSum
                ad_NorthSum=0.0_r8
              END DO
!>            tl_NorthSum=0.0_r8
!>
              ad_NorthSum=0.0_r8
              m=BRY_CONTACT(inorth,cr)%C2Bindex(Ibf)    ! pick last one
!
!  Load coarser grid northern boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
!>            BRY_CONTACT(inorth,cjcr)%tl_Mflux(Ibc)=tl_NorthSum
!>
              ad_NorthSum=ad_NorthSum+                                  &
     &                    BRY_CONTACT(inorth,cjcr)%ad_Mflux(Ibc)
              BRY_CONTACT(inorth,cjcr)%ad_Mflux(Ibc)=0.0_r8

# ifdef NESTING_DEBUG_NOT
              IF (NorthSum.ne.0) THEN
                MFratio=REFINED(cr)%DV_avg2(1,m,tnew)/NorthSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (302,30) Ibc, REFINED(cr)%DV_avg2(1,m,tnew),        &
     &                       NorthSum, MFratio
# endif
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
          BRY_CONTACT(iwest ,cjcr)%ad_Mflux=0.0_r8
          BRY_CONTACT(ieast ,cjcr)%ad_Mflux=0.0_r8
          BRY_CONTACT(isouth,cjcr)%ad_Mflux=0.0_r8
          BRY_CONTACT(inorth,cjcr)%ad_Mflux=0.0_r8
# endif
# ifdef DISTRIBUTE
!
!  Collect data from all nodes.
!
!>        CALL mp_assemble (dg, model, NptsWE, spv,                     &
!>   &                      BRY_CONTACT(iwest ,cjcr)%tl_Mflux(JLB:))
!>
!>        CALL mp_assemble (dg, model, NptsWE, spv,                     &
!>   &                      BRY_CONTACT(ieast ,cjcr)%tl_Mflux(JLB:))
!>
!>        CALL mp_assemble (dg, model, NptsSN, spv,                     &
!>   &                      BRY_CONTACT(isouth,cjcr)%tl_Mflux(ILB:))
!>
!>        CALL mp_assemble (dg, model, NptsSN, spv,                     &
!>   &                      BRY_CONTACT(inorth,cjcr)%tl_Mflux(ILB:))
!>
# endif
        END IF
      END DO

# ifdef NESTING_DEBUG_NOT
!
      CALL my_flush (302)
!
  10  FORMAT (/,1x,a,/,/,4x,'cr = ',i2.2,4x,'dg = ',i2.2,4x,'rg = ',    &
     &        i2.2,4x,'iif(rg) = ',i3.3,4x,'iic(rg) = ',i10.10,4x,      &
     &        'time(rg) = ',i10.10,/,/,2x,'Coarse',6x,'Coarse Grid',8x, &
     &        'Fine Grid',11x,'Ratio',/,4x,'Jb',9x,'DU_avg2',9x,        &
     &        'SUM(DU_avg2)',/)
  20  FORMAT (/,1x,a,/,/,4x,'cr = ',i2.2,4x,'dg = ',i2.2,4x,'rg = ',    &
     &        i2.2,4x,'iif(rg) = ',i3.3,4x,'iic(rg) = ',i10.10,4x,      &
     &        'time(rg) = ',i10.10,/,/,2x,'Coarse',6x,'Coarse Grid',8x, &
     &        'Fine Grid',11x,'Ratio',/,4x,'Ib',9x,'DV_avg2',9x,        &
     &        'SUM(DV_avg2)',/)
  30  FORMAT (4x,i4.4,3(3x,1p,e15.8))
# endif

      RETURN
      END SUBROUTINE ad_check_massflux
# endif
#endif
      END MODULE ad_nesting_mod
