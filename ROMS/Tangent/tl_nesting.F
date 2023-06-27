#include "cppdefs.h"
      MODULE tl_nesting_mod

#if defined NESTING && defined TANGENT
!
!git $Id$
!svn $Id: tl_nesting.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several tangent linear routines to process the !
!  connectivity between nested grids. It process the contact region    !
!  points between data donor and data receiver grids.                  !
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
!  tl_nesting         Public interface to time-stepping kernel         !
!                                                                      !
!  tl_get_composite   Composite  grid, extracts contact points donor   !
!                       data                                           !
!  tl_get_refine      Refinement grid, extracts contact points donor   !
!                       data                                           !
!  tl_put_composite   Composite  grid, fills contact points (spatial   !
!                       interpolation)                                 !
!  tl_put_refine      Refinement grid, fills contact points (spatial   !
!                       and temporal interpolation)                    !
!                                                                      !
# ifdef NESTING_DEBUG
!  tl_check_massflux  If refinement, checks mass fluxes between coarse !
!                       and fine grids for volume conservation. It is  !
!                       use only for debugging and diagnostics.        !
# endif
!  tl_correct_tracer  Corrects coarse grid tracer at the refinement    !
!                       grid boundary with the refined accumulated     !
!                       fluxes                                         !
!  tl_fine2coarse     Replaces coarse grid state variables with the    !
!                       averaged fine grid values (two-way nesting)    !
!                                                                      !
!  tl_put_refine2d    Interpolates (space-time) 2D state variables     !
!  tl_put_refine3d    Interpolates (space-time) 3D state variables     !
!                                                                      !
!  tl_z_weights       Sets donor grid vertical indices (cell holding   !
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
      PUBLIC  :: tl_nesting
# ifdef NESTING_DEBUG
      PUBLIC  :: tl_check_massflux
# endif
# ifdef SOLVE3D
      PRIVATE :: tl_correct_tracer
      PRIVATE :: tl_correct_tracer_tile
# endif
      PRIVATE :: tl_fine2coarse
      PRIVATE :: tl_get_composite
      PRIVATE :: tl_get_refine
      PRIVATE :: tl_put_composite
      PRIVATE :: tl_put_refine
      PRIVATE :: tl_put_refine2d
# ifdef SOLVE3D
      PRIVATE :: tl_put_refine3d
      PRIVATE :: tl_z_weights
# endif
!
      CONTAINS
!
      SUBROUTINE tl_nesting (ng, model, isection)
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
!
      USE set_depth_mod,    ONLY : set_depth
      USE tl_set_depth_mod, ONLY : tl_set_depth
# endif
      USE nesting_mod,      ONLY : get_metrics
      USE nesting_mod,      ONLY : mask_hweights
      USE strings_mod,      ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, isection
!
!  Local variable declarations.
!
      logical :: LputFsur
!
      integer :: subs, tile, thread
      integer :: ngc
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 36, __LINE__, MyFile)
# endif
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
!^        CALL z_weights (ng, model, tile)
!^
          CALL tl_z_weights (ng, model, tile)
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
!^          CALL get_refine (ng, model, tile)
!^
            CALL tl_get_refine (ng, model, tile)
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
!^          CALL put_refine (ng, model, tile, LputFsur)
!^
            CALL tl_put_refine (ng, model, tile, LputFsur)
          END DO
!$OMP BARRIER
!
          LputFsur=.FALSE.
          DO tile=first_tile(ng),last_tile(ng),+1
!^          CALL put_refine (ng, model, tile, LputFsur)
!^
            CALL tl_put_refine (ng, model, tile, LputFsur)
          END DO
!$OMP BARRIER

# ifdef NESTING_DEBUG
!
!  If refinement, check mass flux conservation between coarser and
!  finer grids.
!
        ELSE IF (isection.eq.nmflx) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
!^          CALL check_massflux (ng, model, tile)
!^
            CALL tl_check_massflux (ng, model, tile)
          END DO
# endif

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
!^          CALL correct_tracer (ngc, ng, model, tile)
!^
            CALL tl_correct_tracer (ngc, ng, model, tile)
          END DO
!$OMP BARRIER
#  endif
!
!  Replace coarse grid 2D state variables with the averaged fine grid
!  values (two-way coupling).
!
          DO tile=last_tile(ngc),first_tile(ngc),-1
!^          CALL fine2coarse (ng, model, r2dvar, tile)
!^
            CALL tl_fine2coarse (ng, model, r2dvar, tile)
          END DO
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
!
!  Update coarse grid depth variables. We have a new coarse grid
!  adjusted free-surface, Zt_avg1.
!
          DO tile=first_tile(ngc),last_tile(ngc),+1
!^          CALL set_depth (ngc, tile, model)
!^
            CALL tl_set_depth (ngc, tile, model)
          END DO
!$OMP BARRIER
!
!  Replace coarse grid 3D state variables with the averaged fine grid
!  values (two-way coupling).
!
          DO tile=last_tile(ngc),first_tile(ngc),-1
!^          CALL fine2coarse (ng, model, r3dvar, tile)
!^
            CALL tl_fine2coarse (ng, model, r3dvar, tile)
          END DO
!$OMP BARRIER
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!^        CALL get_composite (ng, model, isection, tile)
!^
          CALL tl_get_composite (ng, model, isection, tile)
        END DO
!$OMP BARRIER
!
!  Fill composite grid contact points variables by interpolating from
!  extracted donor grid data.
!
        DO tile=last_tile(ng),first_tile(ng),-1
!^        CALL put_composite (ng, model, isection, tile)
!^
          CALL tl_put_composite (ng, model, isection, tile)
        END DO
!$OMP BARRIER

      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 36, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE tl_nesting
!
      SUBROUTINE tl_get_composite (ng, model, isection, tile)
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
!^          CALL get_contact2d (dg, model, tile,                        &
!^   &                          u2dvar, Vname(1,idUbms),                &
!^   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!^   &                          LBi, UBi, LBj, UBj,                     &
!^   &                          FORCES(dg) % bustr,                     &
!^   &                          COMPOSITE(cr) % bustr)
!^
            CALL get_contact2d (dg, model, tile,                        &
     &                          u2dvar, Vname(1,idUbms),                &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(dg) % tl_bustr,                  &
     &                          COMPOSITE(cr) % tl_bustr)
!^          CALL get_contact2d (dg, model, tile,                        &
!^   &                          v2dvar, Vname(1,idVbms),                &
!^   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!^   &                          LBi, UBi, LBj, UBj,                     &
!^   &                          FORCES(dg) % bvstr,                     &
!^   &                          COMPOSITE(cr) % bvstr)
!^
            CALL get_contact2d (dg, model, tile,                        &
     &                          v2dvar, Vname(1,idVbms),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(dg) % tl_bvstr,                  &
     &                          COMPOSITE(cr) % tl_bvstr)
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
!^            CALL get_contact2d (dg, model, tile,                      &
!^   &                            r2dvar, Vname(1,idFsur),              &
!^   &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
!^   &                            LBi, UBi, LBj, UBj,                   &
!^   &                            OCEAN(dg) % zeta(:,:,Tindex),         &
!^   &                            COMPOSITE(cr) % zeta(:,:,rec))
!^
              CALL get_contact2d (dg, model, tile,                      &
     &                            r2dvar, Vname(1,idFsur),              &
     &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OCEAN(dg) % tl_zeta(:,:,Tindex),      &
     &                            COMPOSITE(cr) % tl_zeta(:,:,rec))
            END DO
          END IF
!
!  Process free-surface equation rigth-hand-side (rzeta) term.
!
          IF (isection.eq.n2dPS) THEN
            Tindex=1
!^          CALL get_contact2d (dg, model, tile,                        &
!^   &                          r2dvar, Vname(1,idRzet),                &
!^   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!^   &                          LBi, UBi, LBj, UBj,                     &
!^   &                          OCEAN(dg) % rzeta(:,:,Tindex),          &
!^   &                          COMPOSITE(cr) % rzeta)
!^
            CALL get_contact2d (dg, model, tile,                        &
     &                          r2dvar, Vname(1,idRzet),                &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(dg) % tl_rzeta(:,:,Tindex),       &
     &                          COMPOSITE(cr) % tl_rzeta)
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
!^            CALL get_contact2d (dg, model, tile,                      &
!^   &                            u2dvar, Vname(1,idUbar),              &
!^   &                            cr, Ucontact(cr)%Npoints, Ucontact,   &
!^   &                            LBi, UBi, LBj, UBj,                   &
!^   &                            OCEAN(dg) % ubar(:,:,Tindex),         &
!^   &                            COMPOSITE(cr) % ubar(:,:,rec))
!^
              CALL get_contact2d (dg, model, tile,                      &
     &                            u2dvar, Vname(1,idUbar),              &
     &                            cr, Ucontact(cr)%Npoints, Ucontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OCEAN(dg) % tl_ubar(:,:,Tindex),      &
     &                            COMPOSITE(cr) % tl_ubar(:,:,rec))
!^            CALL get_contact2d (dg, model, tile,                      &
!^   &                            v2dvar, Vname(1,idVbar),              &
!^   &                            cr, Vcontact(cr)%Npoints, Vcontact,   &
!^   &                            LBi, UBi, LBj, UBj,                   &
!^   &                            OCEAN(dg) % vbar(:,:,Tindex),         &
!^   &                            COMPOSITE(cr) % vbar(:,:,rec))
!^
              CALL get_contact2d (dg, model, tile,                      &
     &                            v2dvar, Vname(1,idVbar),              &
     &                            cr, Vcontact(cr)%Npoints, Vcontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OCEAN(dg) % tl_vbar(:,:,Tindex),      &
     &                            COMPOSITE(cr) % tl_vbar(:,:,rec))
            END DO
          END IF

# ifdef SOLVE3D
!
!  Process time averaged free-surface (Zt_avg1) and 2D momentum fluxes
!  (DU_avg1, DV_avg1).
!
          IF (isection.eq.n2dfx) THEN
!^          CALL get_contact2d (dg, model, tile,                        &
!^   &                          r2dvar, 'Zt_avg1',                      &
!^   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!^   &                          LBi, UBi, LBj, UBj,                     &
!^   &                          COUPLING(dg) % Zt_avg1,                 &
!^   &                          COMPOSITE(cr) % Zt_avg1)
!^
            CALL get_contact2d (dg, model, tile,                        &
     &                          r2dvar, 'Zt_avg1',                      &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % tl_Zt_avg1,              &
     &                          COMPOSITE(cr) % tl_Zt_avg1)

            CALL get_contact2d (dg, model, tile,                        &
     &                          u2dvar, 'DU_avg1',                      &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % DU_avg1,                 &
     &                          COMPOSITE(cr) % DU_avg1)
            CALL get_contact2d (dg, model, tile,                        &
     &                          u2dvar, 'DU_avg1',                      &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % tl_DU_avg1,              &
     &                          COMPOSITE(cr) % tl_DU_avg1)

            CALL get_contact2d (dg, model, tile,                        &
     &                          v2dvar, 'DV_avg1',                      &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % DV_avg1,                 &
     &                          COMPOSITE(cr) % DV_avg1)
            CALL get_contact2d (dg, model, tile,                        &
     &                          v2dvar, 'DV_avg1',                      &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % tl_DV_avg1,              &
     &                          COMPOSITE(cr) % tl_DV_avg1)
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
!^            CALL get_contact3d (dg, model, tile,                      &
!^   &                            r3dvar, Vname(1,idTvar(itrc)),        &
!^   &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
!^   &                            LBi, UBi, LBj, UBj, 1, N(dg),         &
!^   &                            OCEAN(dg) % t(:,:,:,Tindex,itrc),     &
!^   &                            COMPOSITE(cr) % t(:,:,:,itrc))
!^
              CALL get_contact3d (dg, model, tile,                      &
     &                            r3dvar, Vname(1,idTvar(itrc)),        &
     &                            cr, Rcontact(cr)%Npoints, Rcontact,   &
     &                            LBi, UBi, LBj, UBj, 1, N(dg),         &
     &                            OCEAN(dg) % tl_t(:,:,:,Tindex,itrc),  &
     &                            COMPOSITE(cr) % tl_t(:,:,:,itrc))
            END DO
          END IF
#  endif
!
!  Process 3D momentum (u, v) at the appropriate time-index.
!
          IF ((isection.eq.n3dIC).or.                                   &
     &        (isection.eq.n3duv)) THEN
            Tindex=nnew(dg)
!^          CALL get_contact3d (dg, model, tile,                        &
!^   &                          u3dvar, Vname(1,idUvel),                &
!^   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!^   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!^   &                          OCEAN(dg) % u(:,:,:,Tindex),            &
!^   &                          COMPOSITE(cr) % u)
!^
            CALL get_contact3d (dg, model, tile,                        &
     &                          u3dvar, Vname(1,idUvel),                &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          OCEAN(dg) % tl_u(:,:,:,Tindex),         &
     &                          COMPOSITE(cr) % tl_u)
!^          CALL get_contact3d (dg, model, tile,                        &
!^   &                          v3dvar, Vname(1,idVvel),                &
!^   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!^   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!^   &                          OCEAN(dg) % v(:,:,:,Tindex),            &
!^   &                          COMPOSITE(cr) % v)
!^
            CALL get_contact3d (dg, model, tile,                        &
     &                          v3dvar, Vname(1,idVvel),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          OCEAN(dg) % tl_v(:,:,:,Tindex),         &
     &                          COMPOSITE(cr) % tl_v)
          END IF
!
!  Process 3D momentum fluxes (Huon, Hvom).
!
          IF (isection.eq.n3duv) THEN
!^          CALL get_contact3d (dg, model, tile,                        &
!^   &                          u3dvar, 'Huon',                         &
!^   &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
!^   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!^   &                          GRID(dg) % Huon,                        &
!^   &                          COMPOSITE(cr) % Huon)
!^
            CALL get_contact3d (dg, model, tile,                        &
     &                          u3dvar, 'Huon',                         &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          GRID(dg) % tl_Huon,                     &
     &                          COMPOSITE(cr) % tl_Huon)

!^          CALL get_contact3d (dg, model, tile,                        &
!^   &                          v3dvar, 'Hvom',                         &
!^   &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
!^   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!^   &                          GRID(dg) % Hvom,                        &
!^   &                          COMPOSITE(cr) % Hvom)
!^
            CALL get_contact3d (dg, model, tile,                        &
     &                          v3dvar, 'Hvom',                         &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          GRID(dg) % tl_Hvom,                     &
     &                          COMPOSITE(cr) % tl_Hvom)
          END IF
# endif

        END IF
      END DO
!
      RETURN
      END SUBROUTINE tl_get_composite
!
      SUBROUTINE tl_get_refine (ng, model, tile)
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
      USE nesting_mod, ONLY : get_contact2d
# ifdef SOLVE3D
      USE nesting_mod, ONLY : get_contact3d
# endif
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
!  If debugging, write information into Fortran unit 101 to check the
!  logic of processing donor grid data.
!
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            IF (Master) THEN
              IF (first) THEN
                first=.FALSE.
                WRITE (101,10)
              END IF
              WRITE (101,20) ng, cr, dg, rg, iic(dg), iic(rg),          &
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
              CALL my_flush (101)
            END IF
          END IF
# endif
!
!  Extract free-surface.
!
# ifdef SOLVE3D
!^        CALL get_contact2d (dg, model, tile,                          &
!^   &                        r2dvar, 'Zt_avg1',                        &
!^   &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
!^   &                        LBi, UBi, LBj, UBj,                       &
!^   &                        COUPLING(dg) % Zt_avg1,                   &
!^   &                        REFINED(cr) % zeta(:,:,tnew))
!^
          CALL get_contact2d (dg, model, tile,                          &
     &                        r2dvar, 'Zt_avg1',                        &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        COUPLING(dg) % tl_Zt_avg1,                &
     &                        REFINED(cr) % tl_zeta(:,:,tnew))
# else
!^        CALL get_contact2d (dg, model, tile,                          &
!^   &                        r2dvar, 'zeta',                           &
!^   &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
!^   &                        LBi, UBi, LBj, UBj,                       &
!^   &                        OCEAN(dg) % zeta(:,:,Tindex2d),           &
!^   &                        REFINED(cr) % zeta(:,:,tnew))
!^
          CALL get_contact2d (dg, model, tile,                          &
     &                        r2dvar, 'zeta',                           &
     &                        cr, Rcontact(cr)%Npoints, Rcontact,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(dg) % tl_zeta(:,:,Tindex2d),        &
     &                        REFINED(cr) % tl_zeta(:,:,tnew))
# endif
!
!  Extract 2D momentum components (ubar, vbar).
!
!^        CALL get_contact2d (dg, model, tile,                          &
!^   &                        u2dvar, Vname(1,idUbar),                  &
!^   &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
!^   &                        LBi, UBi, LBj, UBj,                       &
!^   &                        OCEAN(dg) % ubar(:,:,Tindex2d),           &
!^   &                        REFINED(cr) % ubar(:,:,tnew))
!^
          CALL get_contact2d (dg, model, tile,                          &
     &                        u2dvar, Vname(1,idUbar),                  &
     &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(dg) % tl_ubar(:,:,Tindex2d),        &
     &                        REFINED(cr) % tl_ubar(:,:,tnew))

!^        CALL get_contact2d (dg, model, tile,                          &
!^   &                        v2dvar, Vname(1,idVbar),                  &
!^   &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
!^   &                        LBi, UBi, LBj, UBj,                       &
!^   &                        OCEAN(dg) % vbar(:,:,Tindex2d),           &
!^   &                        REFINED(cr) % vbar(:,:,tnew))
!^
          CALL get_contact2d (dg, model, tile,                          &
     &                        v2dvar, Vname(1,idVbar),                  &
     &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OCEAN(dg) % tl_vbar(:,:,Tindex2d),        &
     &                        REFINED(cr) % tl_vbar(:,:,tnew))

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
          CALL get_persisted2d (dg, rg, model, tile,                    &
     &                          u2dvar, 'DU_avg2',                      &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % tl_DU_avg2,              &
     &                          REFINED(cr) % tl_DU_avg2(:,:,tnew))

          CALL get_persisted2d (dg, rg, model, tile,                    &
     &                          v2dvar, 'DV_avg2',                      &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % DV_avg2,                 &
     &                          REFINED(cr) % DV_avg2(:,:,tnew))
          CALL get_persisted2d (dg, rg, model, tile,                    &
     &                          v2dvar, 'DV_avg2',                      &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          COUPLING(dg) % tl_DV_avg2,              &
     &                          REFINED(cr) % tl_DV_avg2(:,:,tnew))
!
!  Tracer-type variables.
!
          DO itrc=1,NT(dg)
!^          CALL get_contact3d (dg, model, tile,                        &
!^   &                          r3dvar, Vname(1,idTvar(itrc)),          &
!^   &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
!^   &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
!^   &                          OCEAN(dg) % t(:,:,:,Tindex3d,itrc),     &
!^   &                          REFINED(cr) % t(:,:,:,tnew,itrc))
!^
            CALL get_contact3d (dg, model, tile,                        &
     &                          r3dvar, Vname(1,idTvar(itrc)),          &
     &                          cr, Rcontact(cr)%Npoints, Rcontact,     &
     &                          LBi, UBi, LBj, UBj, 1, N(dg),           &
     &                          OCEAN(dg) % tl_t(:,:,:,Tindex3d,itrc),  &
     &                          REFINED(cr) % tl_t(:,:,:,tnew,itrc))
          END DO
!
!  Extract 3D momentum components (u, v).
!
!^        CALL get_contact3d (dg, model, tile,                          &
!^   &                        u3dvar, Vname(1,idUvel),                  &
!^   &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
!^   &                        LBi, UBi, LBj, UBj, 1, N(dg),             &
!^   &                        OCEAN(dg) % u(:,:,:,Tindex3d),            &
!^   &                        REFINED(cr) % u(:,:,:,tnew))
!^
          CALL get_contact3d (dg, model, tile,                          &
     &                        u3dvar, Vname(1,idUvel),                  &
     &                        cr, Ucontact(cr)%Npoints, Ucontact,       &
     &                        LBi, UBi, LBj, UBj, 1, N(dg),             &
     &                        OCEAN(dg) % tl_u(:,:,:,Tindex3d),         &
     &                        REFINED(cr) % tl_u(:,:,:,tnew))

!^        CALL get_contact3d (dg, model, tile,                          &
!^   &                        v3dvar, Vname(1,idVvel),                  &
!^   &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
!^   &                        LBi, UBi, LBj, UBj, 1, N(dg),             &
!^   &                        OCEAN(dg) % v(:,:,:,Tindex3d),            &
!^   &                        REFINED(cr) % v(:,:,:,tnew))
!^
          CALL get_contact3d (dg, model, tile,                          &
     &                        v3dvar, Vname(1,idVvel),                  &
     &                        cr, Vcontact(cr)%Npoints, Vcontact,       &
     &                        LBi, UBi, LBj, UBj, 1, N(dg),             &
     &                        OCEAN(dg) % tl_v(:,:,:,Tindex3d),         &
     &                        REFINED(cr) % tl_v(:,:,:,tnew))
# endif
        END IF
      END DO
!
      RETURN
      END SUBROUTINE tl_get_refine
!
      SUBROUTINE tl_put_composite (ng, model, isection, tile)
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
!
      USE nesting_mod, ONLY : put_contact2d

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
     &                          COMPOSITE(cr) % tl_bustr,               &
     &                          FORCES(rg) % tl_bustr)
            CALL put_contact2d (rg, model, tile,                        &
     &                          v2dvar, Vname(1,idVbms),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
# ifdef MASKING
     &                          GRID(rg) % vmask,                       &
# endif
     &                          COMPOSITE(cr) % tl_bvstr,               &
     &                          FORCES(rg) % tl_bvstr)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (rg, tile, model, 2,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          FORCES(rg) % tl_bustr,                  &
     &                          FORCES(rg) % tl_bvstr)
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
     &                            COMPOSITE(cr) % tl_zeta(:,:,rec),     &
     &                            OCEAN(rg) % tl_zeta(:,:,Tindex))
# ifdef DISTRIBUTE
              CALL mp_exchange2d (rg, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(rg), NSperiodic(rg),       &
     &                            OCEAN(rg) % tl_zeta(:,:,Tindex))
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
     &                          COMPOSITE(cr) % tl_rzeta,               &
     &                          OCEAN(rg) % tl_rzeta(:,:,Tindex))
# ifdef DISTRIBUTE
            CALL mp_exchange2d (rg, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg) % tl_rzeta(:,:,Tindex))
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
     &                            COMPOSITE(cr) % tl_ubar(:,:,rec),     &
     &                            OCEAN(rg) % tl_ubar(:,:,Tindex))
              CALL put_contact2d (rg, model, tile,                      &
     &                            v2dvar, Vname(1,idVbar),              &
     &                            cr, Vcontact(cr)%Npoints, Vcontact,   &
     &                            LBi, UBi, LBj, UBj,                   &
# ifdef MASKING
     &                            GRID(rg) % vmask,                     &
# endif
     &                            COMPOSITE(cr) % tl_vbar(:,:,rec),     &
     &                            OCEAN(rg) % tl_vbar(:,:,Tindex))
# ifdef DISTRIBUTE
              CALL mp_exchange2d (rg, tile, model, 2,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(rg), NSperiodic(rg),       &
     &                            OCEAN(rg) % tl_ubar(:,:,Tindex),      &
     &                            OCEAN(rg) % tl_vbar(:,:,Tindex))
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
     &                          COMPOSITE(cr) % tl_Zt_avg1,             &
     &                          COUPLING(rg) % tl_Zt_avg1)
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
     &                          u2dvar, Vname(1,idUfx1),                &
     &                          cr, Ucontact(cr)%Npoints, Ucontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
#  ifdef MASKING
     &                          GRID(rg) % umask,                       &
#  endif
     &                          COMPOSITE(cr) % tl_DU_avg1,             &
     &                          COUPLING(rg) % tl_DU_avg1)
            CALL put_contact2d (rg, model, tile,                        &
     &                          v2dvar, Vname(1,idVfx1),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
#  ifdef MASKING
     &                          GRID(rg) % vmask,                       &
#  endif
     &                          COMPOSITE(cr) % DV_avg1,                &
     &                          COUPLING(rg) % DV_avg1)
            CALL put_contact2d (rg, model, tile,                        &
     &                          v2dvar, Vname(1,idVfx1),                &
     &                          cr, Vcontact(cr)%Npoints, Vcontact,     &
     &                          LBi, UBi, LBj, UBj,                     &
#  ifdef MASKING
     &                          GRID(rg) % vmask,                       &
#  endif
     &                          COMPOSITE(cr) % tl_DV_avg1,             &
     &                          COUPLING(rg) % tl_DV_avg1)
# ifdef DISTRIBUTE
            CALL mp_exchange2d (rg, tile, model, 2,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          COUPLING(rg) % DU_avg1,                 &
     &                          COUPLING(rg) % DV_avg1)
            CALL mp_exchange2d (rg, tile, model, 3,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          COUPLING(rg) % tl_Zt_avg1,              &
     &                          COUPLING(rg) % tl_DU_avg1,              &
     &                          COUPLING(rg) % tl_DV_avg1)
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
              CALL tl_put_contact3d (rg, model, tile,                   &
     &                               r3dvar, Vname(1,idTvar(itrc)),     &
     &                               cr, Rcontact(cr)%Npoints, Rcontact,&
     &                               LBi, UBi, LBj, UBj, 1, N(rg),      &
#   ifdef MASKING
     &                               GRID(rg) % rmask,                  &
#   endif
     &                               COMPOSITE(cr) % t(:,:,:,itrc),     &
     &                               COMPOSITE(cr) % tl_t(:,:,:,itrc),  &
     &                               OCEAN(rg) % tl_t(:,:,:,Tindex,     &
     &                                                itrc))
            END DO
#   ifdef DISTRIBUTE
            CALL mp_exchange4d (rg, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg), 1, NT(rg),&
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg) % tl_t(:,:,:,Tindex,:))
#   endif
          END IF
#  endif
!
!  Process 3D momentum (u, v) at the appropriate time-index.
!
          IF ((isection.eq.n3dIC).or.                                   &
     &        (isection.eq.n3duv)) THEN
            Tindex=nnew(rg)
            CALL tl_put_contact3d (rg, model, tile,                     &
     &                             u3dvar, Vname(1,idUvel),             &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
     &                             GRID(rg) % umask,                    &
#  endif
     &                             COMPOSITE(cr) % u,                   &
     &                             COMPOSITE(cr) % tl_u,                &
     &                             OCEAN(rg) % tl_u(:,:,:,Tindex))
            CALL tl_put_contact3d (rg, model, tile,                     &
     &                             v3dvar, Vname(1,idVvel),             &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
     &                             GRID(rg) % vmask,                    &
#  endif
     &                             COMPOSITE(cr) % v,                   &
     &                             COMPOSITE(cr) % tl_v,                &
     &                             OCEAN(rg) % tl_v(:,:,:,Tindex))
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (rg, tile, model, 2,                     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg) % tl_u(:,:,:,Tindex),         &
     &                          OCEAN(rg) % tl_v(:,:,:,Tindex))
#  endif
          END IF
!
!  Process 3D momentum fluxes (Huon, Hvom).
!
          IF (isection.eq.n3duv) THEN
            CALL tl_put_contact3d (rg, model, tile,                     &
     &                             u3dvar, 'Huon',                      &
     &                             cr, Ucontact(cr)%Npoints, Ucontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
     &                             GRID(rg) % umask,                    &
#  endif
     &                             COMPOSITE(cr) % Huon,                &
     &                             COMPOSITE(cr) % tl_Huon,             &
     &                             GRID(rg) % tl_Huon)
            CALL tl_put_contact3d (rg, model, tile,                     &
     &                             v3dvar, 'Hvom',                      &
     &                             cr, Vcontact(cr)%Npoints, Vcontact,  &
     &                             LBi, UBi, LBj, UBj, 1, N(rg),        &
#  ifdef MASKING
     &                             GRID(rg) % vmask,                    &
#  endif
     &                             COMPOSITE(cr) % Hvom,                &
     &                             COMPOSITE(cr) % tl_Hvom,             &
     &                             GRID(rg) % tl_Hvom)
#  ifdef DISTRIBUTE
            CALL mp_exchange3d (rg, tile, model, 2,                     &
     &                          LBi, UBi, LBj, UBj, 1, N(rg),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          GRID(rg) % tl_Huon,                     &
     &                          GRID(rg) % tl_Hvom)
#  endif
          END IF
# endif

        END IF
      END DO CR_LOOP
!
      RETURN
      END SUBROUTINE tl_put_composite
!
      SUBROUTINE tl_put_refine (ng, model, tile, LputFsur)
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
            CALL tl_put_refine2d (ng, dg, cr, model, tile, LputFsur,    &
     &                            LBi, UBi, LBj, UBj)
          ELSE
!
!  Fill other 2D state variables (like momentum) contact points.
!
            CALL tl_put_refine2d (ng, dg, cr, model, tile, LputFsur,    &
     &                            LBi, UBi, LBj, UBj)

# ifdef SOLVE3D
!
!  Fill 3D state variables contact points.
!
            CALL tl_put_refine3d (ng, dg, cr, model, tile,              &
     &                            LBi, UBi, LBj, UBj)
# endif
          END IF
        END IF
      END DO
!
      RETURN
      END SUBROUTINE tl_put_refine

# ifdef SOLVE3D
!
      SUBROUTINE tl_correct_tracer (ng, ngf, model, tile)
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
      CALL tl_correct_tracer_tile (ng, ngf, model, tile,                &
                                   LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS)
!
      RETURN
      END SUBROUTINE tl_correct_tracer
!
!***********************************************************************
      SUBROUTINE tl_correct_tracer_tile (ngc, ngf, model, tile,         &
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
      real(r8) :: tl_TFC, tl_TFF, tl_Tvalue, tl_cff

      real(r8) :: Dinv(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: tl_Dinv(IminS:ImaxS,JminS:JmaxS)
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
          tl_cff=GRID(ngc)%tl_Hz(i,j,1)
          DO k=2,N(rg)
            cff=cff+GRID(ngc)%Hz(i,j,k)
            tl_cff=tl_cff+GRID(ngc)%tl_Hz(i,j,k)
          END DO
          Dinv(i,j)=1.0_r8/cff
          tl_Dinv(i,j)=-tl_cff*Dinv(i,j)/cff
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
            tl_TFC=0.0_r8
            DO k=1,N(ngc)
              TFC=TFC+BRY_CONTACT(iwest,rgcr)%Tflux(Jbc,k,itrc)
              tl_TFC=tl_TFC+BRY_CONTACT(iwest,rgcr)%tl_Tflux(Jbc,k,itrc)
            END DO
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal J-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the J-edge.
!
            TFF=0.0_r8
            tl_TFF=0.0_r8
            Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
            DO jsum=-half,half
              Jbf=Jedge+jsum
              DO k=1,N(ngf)
                TFF=TFF+BRY_CONTACT(iwest,dgcr)%Tflux(Jbf,k,itrc)
                tl_TFF=tl_TFF+                                          &
     &                 BRY_CONTACT(iwest,dgcr)%tl_Tflux(Jbf,k,itrc)
              END DO
            END DO
!
!  Zeroth order correction to fine grid time integral (RIL, 2016).
!
            TFF=TFF*dt(ngc)/dt(ngf)
            tl_TFF=tl_TFF*dt(ngc)/dt(ngf)
!
!  Correct coarse grid tracer at the finer grid western boundary.
!
            cff=GRID(ngc)%pm(Ibc-1,Jbc)*                                &
     &          GRID(ngc)%pn(Ibc-1,Jbc)*                                &
     &          Dinv(Ibc-1,Jbc)
            tl_cff=GRID(ngc)%pm(Ibc-1,Jbc)*                             &
     &             GRID(ngc)%pn(Ibc-1,Jbc)*                             &
     &             tl_Dinv(Ibc-1,Jbc)
            DO k=1,N(ngc)
              Tvalue=MAX(0.0_r8,                                        &
     &                   OCEAN(ngc)%t(Ibc-1,Jbc,k,Tindex,itrc)-         &
     &                   cff*(TFF-TFC))
              tl_Tvalue=(0.5_r8-                                        &
     &             SIGN(0.5_r8,-(OCEAN(ngc)%t(Ibc-1,Jbc,k,Tindex,itrc)- &
     &                   cff*(TFF-TFC))))*                              &
     &                   (OCEAN(ngc)%tl_t(Ibc-1,Jbc,k,Tindex,itrc)-     &
     &                   tl_cff*(TFF-TFC)-cff*(tl_TFF-tl_TFC))
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
                Tvalue=Tvalue+                                          &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc-1,Jbc,k,itrc)*   &
     &                 (CLIMA(ngc)%tclm(Ibc-1,Jbc,k,ic)-Tvalue)
                tl_Tvalue=tl_Tvalue-                                    &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc-1,Jbc,k,itrc)*   &
     &                 tl_Tvalue
              END IF
#  ifdef MASKING
              Tvalue=Tvalue*GRID(ngc)%rmask(Ibc-1,Jbc)
              tl_Tvalue=tl_Tvalue*GRID(ngc)%rmask(Ibc-1,Jbc)
#  endif
!^            OCEAN(ngc)%t(Ibc-1,Jbc,k,Tindex,itrc)=Tvalue
!^
              OCEAN(ngc)%tl_t(Ibc-1,Jbc,k,Tindex,itrc)=tl_Tvalue
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
            tl_TFC=0.0_r8
            DO k=1,N(ngc)
              TFC=TFC+BRY_CONTACT(ieast,rgcr)%Tflux(Jbc,k,itrc)
              tl_TFC=tl_TFC+BRY_CONTACT(ieast,rgcr)%tl_Tflux(Jbc,k,itrc)
            END DO
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal J-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the J-edge.
!
            TFF=0.0_r8
            tl_TFF=0.0_r8
            Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
            DO jsum=-half,half
              Jbf=Jedge+jsum
              DO k=1,N(ngf)
                TFF=TFF+BRY_CONTACT(ieast,dgcr)%Tflux(Jbf,k,itrc)
                tl_TFF=tl_TFF+                                          &
     &                 BRY_CONTACT(ieast,dgcr)%tl_Tflux(Jbf,k,itrc)
              END DO
            END DO
!
!  Zeroth order correction to fine grid time integral (RIL, 2016).
!
            TFF=TFF*dt(ngc)/dt(ngf)
            tl_TFF=tl_TFF*dt(ngc)/dt(ngf)
!
!  Correct coarse grid tracer at the finer grid eastern boundary.
!
            cff=GRID(ngc)%pm(Ibc,Jbc)*                                  &
     &          GRID(ngc)%pn(Ibc,Jbc)*                                  &
     &          Dinv(Ibc,Jbc)
            tl_cff=GRID(ngc)%pm(Ibc,Jbc)*                               &
     &             GRID(ngc)%pn(Ibc,Jbc)*                               &
     &             tl_Dinv(Ibc,Jbc)
            DO k=1,N(ngc)
              Tvalue=MAX(0.0_r8,                                        &
     &                   OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-           &
     &                   cff*(TFF-TFC))
              tl_Tvalue=(0.5_r8-                                        &
     &             SIGN(0.5_r8,-(OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-   &
     &                   cff*(TFF-TFC))))*                              &
     &                   (OCEAN(ngc)%tl_t(Ibc,Jbc,k,Tindex,itrc)-       &
     &                   tl_cff*(TFF-TFC)-cff*(tl_TFF-tl_TFC))
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
                Tvalue=Tvalue+                                          &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc)*     &
     &                 (CLIMA(ngc)%tclm(Ibc,Jbc,k,ic)-Tvalue)
                tl_Tvalue=tl_Tvalue-                                    &
     &                    dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc)*  &
     &                    tl_Tvalue
              END IF
#  ifdef MASKING
              Tvalue=Tvalue*GRID(ngc)%rmask(Ibc,Jbc)
              tl_Tvalue=tl_Tvalue*GRID(ngc)%rmask(Ibc,Jbc)
#  endif
!^            OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)=Tvalue
!^
              OCEAN(ngc)%tl_t(Ibc,Jbc,k,Tindex,itrc)=tl_Tvalue
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
            tl_TFC=0.0_r8
            DO k=1,N(ngc)
             TFC=TFC+BRY_CONTACT(isouth,rgcr)%Tflux(Ibc,k,itrc)
             tl_TFC=tl_TFC+BRY_CONTACT(isouth,rgcr)%tl_Tflux(Ibc,k,itrc)
            END DO
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal I-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the I-edge.
!
            TFF=0.0_r8
            tl_TFF=0.0_r8
            Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
            DO isum=-half,half
              Ibf=Iedge+isum
              DO k=1,N(ngf)
                TFF=TFF+BRY_CONTACT(isouth,dgcr)%Tflux(Ibf,k,itrc)
                tl_TFF=tl_TFF+                                          &
     &                 BRY_CONTACT(isouth,dgcr)%tl_Tflux(Ibf,k,itrc)
              END DO
            END DO
!
!  Zeroth order correction to fine grid time integral (RIL, 2016).
!
            TFF=TFF*dt(ngc)/dt(ngf)
            tl_TFF=tl_TFF*dt(ngc)/dt(ngf)
!
!  Correct coarse grid tracer at the finer grid southern boundary.
!
            cff=GRID(ngc)%pm(Ibc,Jbc-1)*                                &
     &          GRID(ngc)%pn(Ibc,Jbc-1)*                                &
     &          Dinv(Ibc,Jbc-1)
            tl_cff=GRID(ngc)%pm(Ibc,Jbc-1)*                             &
     &             GRID(ngc)%pn(Ibc,Jbc-1)*                             &
     &             tl_Dinv(Ibc,Jbc-1)
            DO k=1,N(ngc)
              Tvalue=MAX(0.0_r8,                                        &
     &                   OCEAN(ngc)%t(Ibc,Jbc-1,k,Tindex,itrc)-         &
     &                   cff*(TFF-TFC))
              tl_Tvalue=(0.5_r8-                                        &
     &             SIGN(0.5_r8,-(OCEAN(ngc)%t(Ibc,Jbc-1,k,Tindex,itrc)- &
     &                   cff*(TFF-TFC))))*                              &
     &                   (OCEAN(ngc)%tl_t(Ibc,Jbc-1,k,Tindex,itrc)-     &
     &                   tl_cff*(TFF-TFC)-cff*(tl_TFF-tl_TFC))
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
                Tvalue=Tvalue+                                          &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc-1,k,itrc)*   &
     &                 (CLIMA(ngc)%tclm(Ibc,Jbc-1,k,ic)-Tvalue)
                tl_Tvalue=tl_Tvalue-                                    &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc-1,k,itrc)*   &
     &                 tl_Tvalue
              END IF
#  ifdef MASKING
              Tvalue=Tvalue*GRID(ngc)%rmask(Ibc,Jbc-1)
              tl_Tvalue=tl_Tvalue*GRID(ngc)%rmask(Ibc,Jbc-1)
#  endif
!^            OCEAN(ngc)%t(Ibc,Jbc-1,k,Tindex,itrc)=Tvalue
!^
              OCEAN(ngc)%tl_t(Ibc,Jbc-1,k,Tindex,itrc)=tl_Tvalue
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
            tl_TFC=0.0_r8
            DO k=1,N(ngc)
              TFC=TFC+BRY_CONTACT(inorth,rgcr)%Tflux(Ibc,k,itrc)
              tl_TFC=tl_TFC+                                            &
     &               BRY_CONTACT(inorth,rgcr)%tl_Tflux(Ibc,k,itrc)
            END DO
!
!  Sum vertically and horizontally finer grid advective tracer flux.
!  This is a vertical and horizontal I-integral because "RefineScale"
!  sub-divisions are done in the finer grid in each single coarse grid
!  at the I-edge.
!
            TFF=0.0_r8
            tl_TFF=0.0_r8
            Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
            DO isum=-half,half
              Ibf=Iedge+isum
              DO k=1,N(ngf)
                TFF=TFF+BRY_CONTACT(inorth,dgcr)%Tflux(Ibf,k,itrc)
                tl_TFF=tl_TFF+                                          &
     &                 BRY_CONTACT(inorth,dgcr)%tl_Tflux(Ibf,k,itrc)
              END DO
            END DO
!
!  Zeroth order correction to fine grid time integral.
!
            TFF=TFF*dt(ngc)/dt(ngf)
            tl_TFF=tl_TFF*dt(ngc)/dt(ngf)
!
!  Correct coarse grid tracer at the finer grid northern boundary.
!
            cff=GRID(ngc)%pm(Ibc,Jbc)*                                  &
     &          GRID(ngc)%pn(Ibc,Jbc)*                                  &
     &          Dinv(Ibc,Jbc)
            tl_cff=GRID(ngc)%pm(Ibc,Jbc)*                               &
     &             GRID(ngc)%pn(Ibc,Jbc)*                               &
     &             tl_Dinv(Ibc,Jbc)
            DO k=1,N(ngc)
              Tvalue=MAX(0.0_r8,                                        &
     &                   OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-           &
     &                   cff*(TFF-TFC))
              tl_Tvalue=(0.5_r8-                                        &
     &              SIGN(0.5_r8,-(OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)-  &
     &                   cff*(TFF-TFC))))*                              &
     &                   (OCEAN(ngc)%tl_t(Ibc,Jbc,k,Tindex,itrc)-       &
     &                   tl_cff*(TFF-TFC)-cff*(tl_TFF-tl_TFC))
              IF (LtracerCLM(itrc,ngc).and.LnudgeTCLM(itrc,ngc)) THEN
                Tvalue=Tvalue+                                          &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc)*     &
     &                 (CLIMA(ngc)%tclm(Ibc,Jbc,k,ic)-Tvalue)
                tl_Tvalue=tl_Tvalue-                                    &
     &                 dt(ngc)*CLIMA(ngc)%Tnudgcof(Ibc,Jbc,k,itrc)*     &
     &                 tl_Tvalue
              END IF
#  ifdef MASKING
              Tvalue=Tvalue*GRID(ngc)%rmask(Ibc,Jbc)
              tl_Tvalue=tl_Tvalue*GRID(ngc)%rmask(Ibc,Jbc)
#  endif
!^            OCEAN(ngc)%t(Ibc,Jbc,k,Tindex,itrc)=Tvalue
!^
              OCEAN(ngc)%tl_t(Ibc,Jbc,k,Tindex,itrc)=tl_Tvalue
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
!^    CALL mp_exchange4d (ngc, tile, model, 1,                          &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ngc),                &
!^   &                    1, NT(ngc),                                   &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ngc), NSperiodic(ngc),             &
!^   &                    OCEAN(ngc)%t(:,:,:,Tindex,:))
!^
      CALL mp_exchange4d (ngc, tile, model, 1,                          &
     &                    LBi, UBi, LBj, UBj, 1, N(ngc),                &
     &                    1, NT(ngc),                                   &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ngc), NSperiodic(ngc),             &
     &                    OCEAN(ngc)%tl_t(:,:,:,Tindex,:))
#  endif
!
      RETURN
      END SUBROUTINE tl_correct_tracer_tile
# endif
!
      SUBROUTINE tl_fine2coarse (ng, model, vtype, tile)
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

      USE nesting_mod, ONLY : fine2coarse2d
# ifdef SOLVE3D
      USE nesting_mod, ONLY : fine2coarse3d
# endif
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
      USE strings_mod,     ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, vtype, tile
!
!  Local variable declarations.
!
      logical :: AreaAvg
!
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
     &  __FILE__//", tl_fine2coarse"
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
  10          FORMAT (6x,'TL_FINE2COARSE - exchanging data between ',   &
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
     &                          GRID(rg)%rmask,                         &
#  endif
     &                          COUPLING(dg)%tl_Zt_avg1,                &
     &                          COUPLING(rg)%tl_Zt_avg1)
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
     &                          GRID(dg)%rmask,                         &
     &                          GRID(rg)%rmask,                         &
#  endif
     &                          OCEAN(dg)%tl_zeta(:,:,Dindex2d),        &
     &                          OCEAN(rg)%tl_zeta(:,:,Rindex2d))
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
     &                          OCEAN(dg)%tl_ubar(:,:,Dindex2d),        &
# ifdef SOLVE3D
     &                          OCEAN(rg)%tl_ubar(:,:,1),               &
     &                          OCEAN(rg)%tl_ubar(:,:,2))
# else
     &                          OCEAN(rg)%tl_ubar(:,:,Rindex2d))
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
     &                          OCEAN(dg)%tl_vbar(:,:,Dindex2d),        &
# ifdef SOLVE3D
     &                          OCEAN(rg)%tl_vbar(:,:,1),               &
     &                          OCEAN(rg)%tl_vbar(:,:,2))
# else
     &                          OCEAN(rg)%tl_vbar(:,:,Rindex2d))
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
     &                            GRID(dg)%rmask,                       &
     &                            GRID(rg)%rmask,                       &
#  endif
     &                            OCEAN(dg)%tl_t(:,:,:,Dindex3d,itrc),  &
     &                            OCEAN(rg)%tl_t(:,:,:,Rindex3d,itrc))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
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
     &                          OCEAN(dg)%tl_u(:,:,:,Dindex3d),         &
     &                          OCEAN(rg)%tl_u(:,:,:,Rindex3d))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
     &                          OCEAN(dg)%tl_v(:,:,:,Dindex3d),         &
     &                          OCEAN(rg)%tl_v(:,:,:,Rindex3d))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
          END IF
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
          IF (EWperiodic(rg).or.NSperiodic(rg)) THEN
            IF (vtype.eq.r2dvar) THEN
# ifdef SOLVE3D
!^            CALL exchange_r2d_tile (rg, tile,                         &
!^   &                                LBiR, UBiR, LBjR, UBjR,           &
!^   &                                COUPLING(rg)%Zt_avg1)
!^
              CALL exchange_r2d_tile (rg, tile,                         &
     &                                LBiR, UBiR, LBjR, UBjR,           &
     &                                COUPLING(rg)%tl_Zt_avg1)
              DO k=1,2
!^              CALL exchange_u2d_tile (rg, tile,                       &
!^   &                                  LBiR, UBiR, LBjR, UBjR,         &
!^   &                                  OCEAN(rg)%ubar(:,:,k))
!^
                CALL exchange_u2d_tile (rg, tile,                       &
     &                                  LBiR, UBiR, LBjR, UBjR,         &
     &                                  OCEAN(rg)%tl_ubar(:,:,k))
!^              CALL exchange_v2d_tile (rg, tile,                       &
!^   &                                  LBiR, UBiR, LBjR, UBjR,         &
!^   &                                  OCEAN(rg)%vbar(:,:,k))
!^
                CALL exchange_v2d_tile (rg, tile,                       &
     &                                  LBiR, UBiR, LBjR, UBjR,         &
     &                                  OCEAN(rg)%tl_vbar(:,:,k))
              END DO
# else
>             CALL exchange_r2d_tile (rg, tile,                         &
!^   &                                LBiR, UBiR, LBjR, UBjR,           &
!^   &                                OCEAN(rg)%zeta(:,:,Rindex2d))
!^
              CALL exchange_r2d_tile (rg, tile,                         &
     &                                LBiR, UBiR, LBjR, UBjR,           &
     &                                OCEAN(rg)%tl_zeta(:,:,Rindex2d))
!^            CALL exchange_u2d_tile (rg, tile,                         &
!^   &                                LBiR, UBiR, LBjR, UBjR,           &
!^   &                                OCEAN(rg)%ubar(:,:,Rindex2d))
!^
              CALL exchange_u2d_tile (rg, tile,                         &
     &                                LBiR, UBiR, LBjR, UBjR,           &
     &                                OCEAN(rg)%tl_ubar(:,:,Rindex2d))
!^            CALL exchange_v2d_tile (rg, tile,                         &
!^   &                                LBiR, UBiR, LBjR, UBjR,           &
!^   &                                OCEAN(rg)%vbar(:,:,Rindex2d))
!^
              CALL exchange_v2d_tile (rg, tile,                         &
     &                                LBiR, UBiR, LBjR, UBjR,           &
     &                                OCEAN(rg)%tl_vbar(:,:,Rindex2d))

# endif
# ifdef SOLVE3D
            ELSE IF (vtype.eq.r3dvar) THEN
!^            CALL exchange_u3d_tile (rg, tile,                         &
!^   &                                LBiR, UBiR, LBjR, UBjR, 1, N(rg), &
!^   &                                OCEAN(rg)%u(:,:,:,Rindex3d))
!^
              CALL exchange_u3d_tile (rg, tile,                         &
     &                                LBiR, UBiR, LBjR, UBjR, 1, N(rg), &
     &                                OCEAN(rg)%tl_u(:,:,:,Rindex3d))
!^            CALL exchange_v3d_tile (rg, tile,                         &
!^   &                                LBiR, UBiR, LBjR, UBjR, 1, N(rg), &
!^   &                                OCEAN(rg)%v(:,:,:,Rindex3d))
!^
              CALL exchange_v3d_tile (rg, tile,                         &
     &                                LBiR, UBiR, LBjR, UBjR, 1, N(rg), &
     &                                OCEAN(rg)%tl_v(:,:,:,Rindex3d))
              DO itrc=1,NT(rg)
!^              CALL exchange_r3d_tile (rg, tile,                       &
!^   &                                  LBiR, UBiR, LBjR, UBjR,         &
!^                                      1, N(rg),                       &
!^   &                                  OCEAN(rg)%t(:,:,:,Rindex3d,     &
!^   &                                              itrc))
!^
                CALL exchange_r3d_tile (rg, tile,                       &
     &                                  LBiR, UBiR, LBjR, UBjR,         &
     &                                  1, N(rg),                       &
     &                                  OCEAN(rg)%tl_t(:,:,:,Rindex3d,  &
     &                                                 itrc))
              END DO
# endif
            END IF
          END IF

# ifdef DISTRIBUTE
!
          IF (vtype.eq.r2dvar) THEN
#  ifdef SOLVE3D
!^          CALL mp_exchange2d (rg, tile, model, 1,                     &
!^   &                          LBiR, UBiR, LBjR, UBjR,                 &
!^   &                          NghostPoints,                           &
!^   &                          EWperiodic(rg), NSperiodic(rg),         &
!^   &                          COUPLING(rg)%Zt_avg1)
!^
            CALL mp_exchange2d (rg, tile, model, 1,                     &
     &                          LBiR, UBiR, LBjR, UBjR,                 &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          COUPLING(rg)%tl_Zt_avg1)
!^          CALL mp_exchange2d (rg, tile, model, 4,                     &
!^   &                          LBiR, UBiR, LBjR, UBjR,                 &
!^   &                          NghostPoints,                           &
!^   &                          EWperiodic(rg), NSperiodic(rg),         &
!^   &                          OCEAN(rg)%ubar(:,:,1),                  &
!^   &                          OCEAN(rg)%vbar(:,:,1),                  &
!^   &                          OCEAN(rg)%ubar(:,:,2),                  &
!^   &                          OCEAN(rg)%vbar(:,:,2))
!^
            CALL mp_exchange2d (rg, tile, model, 4,                     &
     &                          LBiR, UBiR, LBjR, UBjR,                 &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg)%tl_ubar(:,:,1),               &
     &                          OCEAN(rg)%tl_vbar(:,:,1),               &
     &                          OCEAN(rg)%tl_ubar(:,:,2),               &
     &                          OCEAN(rg)%tl_vbar(:,:,2))
#  else
!^          CALL mp_exchange2d (rg, tile, model, 3,                     &
!^   &                          LBiR, UBiR, LBjR, UBjR,                 &
!^   &                          NghostPoints,                           &
!^   &                          EWperiodic(rg), NSperiodic(rg),         &
!^   &                          OCEAN(rg)%zeta(:,:,Rindex2d),           &
!^   &                          OCEAN(rg)%ubar(:,:,Rindex2d),           &
!^   &                          OCEAN(rg)%vbar(:,:,Rindex2d))
!^
            CALL mp_exchange2d (rg, tile, model, 3,                     &
     &                          LBiR, UBiR, LBjR, UBjR,                 &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg)%tl_zeta(:,:,Rindex2d),        &
     &                          OCEAN(rg)%tl_ubar(:,:,Rindex2d),        &
     &                          OCEAN(rg)%tl_vbar(:,:,Rindex2d))
#  endif
#  ifdef SOLVE3D
          ELSE IF (vtype.eq.r3dvar) THEN
!^          CALL mp_exchange3d (rg, tile, model, 2,                     &
!^   &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
!^   &                          NghostPoints,                           &
!^   &                          EWperiodic(rg), NSperiodic(rg),         &
!^   &                          OCEAN(rg)%u(:,:,:,Rindex3d),            &
!^   &                          OCEAN(rg)%v(:,:,:,Rindex3d))
!^
            CALL mp_exchange3d (rg, tile, model, 2,                     &
     &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg)%tl_u(:,:,:,Rindex3d),         &
     &                          OCEAN(rg)%tl_v(:,:,:,Rindex3d))
!^          CALL mp_exchange4d (rg, tile, model, 1,                     &
!^   &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
!^   &                          1, NT(rg),                              &
!^   &                          NghostPoints,                           &
!^   &                          EWperiodic(rg), NSperiodic(rg),         &
!^   &                          OCEAN(rg)%t(:,:,:,Rindex3d,:))
!^
            CALL mp_exchange4d (rg, tile, model, 1,                     &
     &                          LBiR, UBiR, LBjR, UBjR, 1, N(rg),       &
     &                          1, NT(rg),                              &
     &                          NghostPoints,                           &
     &                          EWperiodic(rg), NSperiodic(rg),         &
     &                          OCEAN(rg)%tl_t(:,:,:,Rindex3d,:))
#  endif
          END IF
# endif
        END IF
      END DO
!
      RETURN
      END SUBROUTINE tl_fine2coarse
!
      SUBROUTINE tl_put_refine2d (ng, dg, cr, model, tile, LputFsur,    &
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
      USE mp_exchange_mod, ONLY : mp_exchange2d
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
!
# ifdef DISTRIBUTE
      integer :: ILB, IUB, JLB, JUB, NptsSN, NptsWE, my_tile
# endif
      integer :: NSUB, i, irec, j, m, tnew, told
      integer :: Idg, Jdg
!
# ifdef DISTRIBUTE
      real(r8), parameter :: spv = 0.0_r8
# endif
      real(dp) :: Wnew, Wold, SecScale, fac
      real(r8) :: cff, cff1, tl_cff
      real(r8) :: my_value, tl_my_value
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", tl_put_refined2d"

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
      IF (((Wold*Wnew).lt.0.0_dp).or.((Wold+Wnew).le.0.0_dp)) THEN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!^          my_value=Wold*                                              &
!^   &               (Rcontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%zeta(1,m,told)+              &
!^   &                Rcontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%zeta(2,m,told)+              &
!^   &                Rcontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%zeta(3,m,told)+              &
!^   &                Rcontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%zeta(4,m,told))+             &
!^   &               Wnew*                                              &
!^   &               (Rcontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%zeta(1,m,tnew)+              &
!^   &                Rcontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%zeta(2,m,tnew)+              &
!^   &                Rcontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%zeta(3,m,tnew)+              &
!^   &                Rcontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%zeta(4,m,tnew))
!^
            tl_my_value=Wold*                                           &
     &                  (Rcontact(cr)%Lweight(1,m)*                     &
     &                            REFINED(cr)%tl_zeta(1,m,told)+        &
     &                   Rcontact(cr)%Lweight(2,m)*                     &
     &                            REFINED(cr)%tl_zeta(2,m,told)+        &
     &                   Rcontact(cr)%Lweight(3,m)*                     &
     &                            REFINED(cr)%tl_zeta(3,m,told)+        &
     &                   Rcontact(cr)%Lweight(4,m)*                     &
     &                            REFINED(cr)%tl_zeta(4,m,told))+       &
     &                  Wnew*                                           &
     &                  (Rcontact(cr)%Lweight(1,m)*                     &
     &                            REFINED(cr)%tl_zeta(1,m,tnew)+        &
     &                   Rcontact(cr)%Lweight(2,m)*                     &
     &                            REFINED(cr)%tl_zeta(2,m,tnew)+        &
     &                   Rcontact(cr)%Lweight(3,m)*                     &
     &                            REFINED(cr)%tl_zeta(3,m,tnew)+        &
     &                   Rcontact(cr)%Lweight(4,m)*                     &
     &                            REFINED(cr)%tl_zeta(4,m,tnew))
# ifdef MASKING
!^          my_value=my_value*GRID(ng)%rmask(i,j)
!^
            tl_my_value=tl_my_value*GRID(ng)%rmask(i,j)
# endif
# ifdef WET_DRY
            IF (my_value.le.(Dcrit(ng)-GRID(ng)%h(i,j))) THEN
!^            my_value=Dcrit(ng)-GRID(ng)%h(i,j)
!^
              tl_my_value=-GRID(ng)%tl_h(i,j)
            END IF
# endif
!^          OCEAN(ng)%zeta(i,j,1)=my_value
!^
            OCEAN(ng)%tl_zeta(i,j,1)=tl_my_value
!^          OCEAN(ng)%zeta(i,j,2)=my_value
!^
            OCEAN(ng)%tl_zeta(i,j,2)=tl_my_value
!^          OCEAN(ng)%zeta(i,j,3)=my_value
!^
            OCEAN(ng)%tl_zeta(i,j,3)=tl_my_value
# ifdef SOLVE3D
!^          COUPLING(ng)%Zt_avg1(i,j)=my_value
!^
            COUPLING(ng)%tl_Zt_avg1(i,j)=tl_my_value
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
!^          my_value=Wold*                                              &
!^   &               (Ucontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%ubar(1,m,told)+              &
!^   &                Ucontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%ubar(2,m,told)+              &
!^   &                Ucontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%ubar(3,m,told)+              &
!^   &                Ucontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%ubar(4,m,told))+             &
!^   &               Wnew*                                              &
!^   &               (Ucontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%ubar(1,m,tnew)+              &
!^   &                Ucontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%ubar(2,m,tnew)+              &
!^   &                Ucontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%ubar(3,m,tnew)+              &
!^   &                Ucontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%ubar(4,m,tnew))
!^
            tl_my_value=Wold*                                           &
     &                  (Ucontact(cr)%Lweight(1,m)*                     &
     &                            REFINED(cr)%tl_ubar(1,m,told)+        &
     &                   Ucontact(cr)%Lweight(2,m)*                     &
     &                            REFINED(cr)%tl_ubar(2,m,told)+        &
     &                   Ucontact(cr)%Lweight(3,m)*                     &
     &                            REFINED(cr)%tl_ubar(3,m,told)+        &
     &                   Ucontact(cr)%Lweight(4,m)*                     &
     &                            REFINED(cr)%tl_ubar(4,m,told))+       &
     &                  Wnew*                                           &
     &                  (Ucontact(cr)%Lweight(1,m)*                     &
     &                            REFINED(cr)%tl_ubar(1,m,tnew)+        &
     &                   Ucontact(cr)%Lweight(2,m)*                     &
     &                            REFINED(cr)%tl_ubar(2,m,tnew)+        &
     &                   Ucontact(cr)%Lweight(3,m)*                     &
     &                            REFINED(cr)%tl_ubar(3,m,tnew)+        &
     &                   Ucontact(cr)%Lweight(4,m)*                     &
     &                            REFINED(cr)%tl_ubar(4,m,tnew))
# ifdef MASKING
!^          my_value=my_value*GRID(ng)%umask(i,j)
!^
            tl_my_value=tl_my_value*GRID(ng)%umask(i,j)
# endif
# ifdef WET_DRY
!^          my_value=my_value*GRID(ng)%umask_wet(i,j)
!^
            tl_my_value=tl_my_value*GRID(ng)%umask_wet(i,j)
# endif
            DO irec=1,3
# ifdef SOLVE3D
              Uboundary=(m.eq.BRY_CONTACT(iwest,cr)%C2Bindex(j)).or.    &
     &                  (m.eq.BRY_CONTACT(ieast,cr)%C2Bindex(j))
              IF(.not.(Uboundary.and.(irec.eq.indx1(ng)))) THEN
!^              OCEAN(ng)%ubar(i,j,irec)=my_value
!^
                OCEAN(ng)%tl_ubar(i,j,irec)=tl_my_value
!!            ELSE                                   ! for debugging
!!              OCEAN(ng)%ubar(i,j,irec)=0.0_r8      ! purposes
              END IF
# else
!^            OCEAN(ng)%ubar(i,j,irec)=my_value
!^
              OCEAN(ng)%tl_ubar(i,j,irec)=tl_my_value
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
!^          my_value=Wold*                                              &
!^   &               (Vcontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%vbar(1,m,told)+              &
!^   &                Vcontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%vbar(2,m,told)+              &
!^   &                Vcontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%vbar(3,m,told)+              &
!^   &                Vcontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%vbar(4,m,told))+             &
!^   &               Wnew*                                              &
!^   &               (Vcontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%vbar(1,m,tnew)+              &
!^   &                Vcontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%vbar(2,m,tnew)+              &
!^   &                Vcontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%vbar(3,m,tnew)+              &
!^   &                Vcontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%vbar(4,m,tnew))
!^
            tl_my_value=Wold*                                           &
     &                  (Vcontact(cr)%Lweight(1,m)*                     &
     &                            REFINED(cr)%tl_vbar(1,m,told)+        &
     &                   Vcontact(cr)%Lweight(2,m)*                     &
     &                            REFINED(cr)%tl_vbar(2,m,told)+        &
     &                   Vcontact(cr)%Lweight(3,m)*                     &
     &                            REFINED(cr)%tl_vbar(3,m,told)+        &
     &                   Vcontact(cr)%Lweight(4,m)*                     &
     &                            REFINED(cr)%tl_vbar(4,m,told))+       &
     &                  Wnew*                                           &
     &                  (Vcontact(cr)%Lweight(1,m)*                     &
     &                            REFINED(cr)%tl_vbar(1,m,tnew)+        &
     &                   Vcontact(cr)%Lweight(2,m)*                     &
     &                            REFINED(cr)%tl_vbar(2,m,tnew)+        &
     &                   Vcontact(cr)%Lweight(3,m)*                     &
     &                            REFINED(cr)%tl_vbar(3,m,tnew)+        &
     &                   Vcontact(cr)%Lweight(4,m)*                     &
     &                            REFINED(cr)%tl_vbar(4,m,tnew))
# ifdef MASKING
!^          my_value=my_value*GRID(ng)%vmask(i,j)
!^
            tl_my_value=tl_my_value*GRID(ng)%vmask(i,j)
# endif
# ifdef WET_DRY
!^          my_value=my_value*GRID(ng)%vmask_wet(i,j)
            tl_my_value=tl_my_value*GRID(ng)%vmask_wet(i,j)
# endif
            DO irec=1,3
# ifdef SOLVE3D
              Vboundary=(m.eq.BRY_CONTACT(isouth,cr)%C2Bindex(i)).or.   &
     &                  (m.eq.BRY_CONTACT(inorth,cr)%C2Bindex(i))
              IF(.not.(Vboundary.and.(irec.eq.indx1(ng)))) THEN
!^              OCEAN(ng)%vbar(i,j,irec)=my_value
!^
                OCEAN(ng)%tl_vbar(i,j,irec)=tl_my_value
!!            ELSE                                   ! for debugging
!!              OCEAN(ng)%vbar(i,j,irec)=0.0_r8      ! purposes
              END IF
# else
!^            OCEAN(ng)%vbar(i,j,irec)=my_value
!^
              OCEAN(ng)%tl_vbar(i,j,irec)=tl_my_value
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
            tl_cff=0.5_r8*GRID(ng)%on_u(Istr,j)*                        &
                (GRID(ng)%tl_h(Istr-1,j)+                               &
     &           OCEAN(ng)%tl_zeta(Istr-1,j,indx1(ng))+                 &
     &           GRID(ng)%tl_h(Istr  ,j)+                               &
     &           OCEAN(ng)%tl_zeta(Istr  ,j,indx1(ng)))
            cff1=GRID(ng)%on_u(Istr,j)/REFINED(cr)%on_u(m)
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DU_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DU_avg2(1,m,tnew))/cff
            tl_my_value=cff1*(Wold*REFINED(cr)%tl_DU_avg2(1,m,told)+    &
     &                     Wnew*REFINED(cr)%tl_DU_avg2(1,m,tnew))/cff-  &
     &                  tl_cff*my_value/cff
#  else
            my_value=cff1*REFINED(cr)%DU_avg2(1,m,tnew)/cff
            tl_my_value=cff1*REFINED(cr)%tl_DU_avg2(1,m,tnew)/cff-      &
     &                  tl_cff*my_value/cff
#  endif
#  ifdef MASKING
            my_value=my_value*GRID(ng)%umask(Istr,j)
            tl_my_value=tl_my_value*GRID(ng)%umask(Istr,j)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%umask_wet(Istr,j)
            tl_my_value=tl_my_value*GRID(ng)%umask_wet(Istr,j)
#  endif
#  ifdef NESTING_DEBUG
!^          BRY_CONTACT(iwest,cr)%Mflux(j)=cff*my_value
!^
            BRY_CONTACT(iwest,cr)%tl_Mflux(j)=cff*tl_my_value+          &
     &                                        tl_cff*my_value
#  endif
!^          OCEAN(ng)%ubar(Istr,j,indx1(ng))=my_value
!^
            OCEAN(ng)%tl_ubar(Istr,j,indx1(ng))=tl_my_value
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
            tl_cff=0.5_r8*GRID(ng)%on_u(Iend+1,j)*                      &
     &             (GRID(ng)%tl_h(Iend+1,j)+                            &
     &              OCEAN(ng)%tl_zeta(Iend+1,j,indx1(ng))+              &
     &              GRID(ng)%tl_h(Iend  ,j)+                            &
     &              OCEAN(ng)%tl_zeta(Iend  ,j,indx1(ng)))
            cff1=GRID(ng)%on_u(Iend+1,j)/REFINED(cr)%on_u(m)
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DU_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DU_avg2(1,m,tnew))/cff
            tl_my_value=cff1*(Wold*REFINED(cr)%tl_DU_avg2(1,m,told)+    &
     &                     Wnew*REFINED(cr)%tl_DU_avg2(1,m,tnew))/cff-  &
     &                     tl_cff*my_value/cff
#  else
            my_value=cff1*REFINED(cr)%DU_avg2(1,m,tnew)/cff
            tl_my_value=cff1*REFINED(cr)%tl_DU_avg2(1,m,tnew)/cff-      &
     &                  tl_cff*my_value/cff
#  endif
#  ifdef MASKING
            my_value=my_value*GRID(ng)%umask(Iend+1,j)
            tl_my_value=tl_my_value*GRID(ng)%umask(Iend+1,j)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%umask_wet(Iend+1,j)
            tl_my_value=tl_my_value*GRID(ng)%umask_wet(Iend+1,j)
#  endif
#  ifdef NESTING_DEBUG
!^          BRY_CONTACT(ieast,cr)%Mflux(j)=cff*my_value
!^
            BRY_CONTACT(ieast,cr)%tl_Mflux(j)=tl_cff*my_value+          &
     &                                        cff*tl_my_value
#  endif
!^          OCEAN(ng)%ubar(Iend+1,j,indx1(ng))=my_value
!^
            OCEAN(ng)%tl_ubar(Iend+1,j,indx1(ng))=tl_my_value
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
            tl_cff=0.5_r8*GRID(ng)%om_v(i,Jstr)*                        &
     &          (GRID(ng)%tl_h(i,Jstr-1)+                               &
     &           OCEAN(ng)%tl_zeta(i,Jstr-1,indx1(ng))+                 &
     &           GRID(ng)%tl_h(i,Jstr  )+                               &
     &           OCEAN(ng)%tl_zeta(i,Jstr  ,indx1(ng)))
            cff1=GRID(ng)%om_v(i,Jstr)/REFINED(cr)%om_v(m)
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DV_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DV_avg2(1,m,tnew))/cff
            tl_my_value=cff1*(Wold*REFINED(cr)%tl_DV_avg2(1,m,told)+    &
     &                     Wnew*REFINED(cr)%tl_DV_avg2(1,m,tnew))/cff-  &
     &                     tl_cff*my_value/cff
#  else
            my_value=cff1*REFINED(cr)%DV_avg2(1,m,tnew)/cff
            tl_my_value=cff1*REFINED(cr)%tl_DV_avg2(1,m,tnew)/cff-      &
     &                  tl_cff*my_value/cff
#  endif
#  ifdef MASKING
            my_value=my_value*GRID(ng)%vmask(i,Jstr)
            tl_my_value=tl_my_value*GRID(ng)%vmask(i,Jstr)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%vmask_wet(i,Jstr)
            tl_my_value=tl_my_value*GRID(ng)%vmask_wet(i,Jstr)
#  endif
#  ifdef NESTING_DEBUG
!^          BRY_CONTACT(isouth,cr)%Mflux(i)=cff*my_value
!^
            BRY_CONTACT(isouth,cr)%tl_Mflux(i)=tl_cff*my_value+         &
     &                                         cff*tl_my_value
#  endif
!^          OCEAN(ng)%vbar(i,Jstr,indx1(ng))=my_value
!^
            OCEAN(ng)%tl_vbar(i,Jstr,indx1(ng))=tl_my_value
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
            tl_cff=0.5_r8*GRID(ng)%om_v(i,Jend+1)*                      &
     &          (GRID(ng)%tl_h(i,Jend+1)+                               &
     &           OCEAN(ng)%tl_zeta(i,Jend+1,indx1(ng))+                 &
     &           GRID(ng)%tl_h(i,Jend  )+                               &
     &           OCEAN(ng)%tl_zeta(i,Jend  ,indx1(ng)))
            cff1=GRID(ng)%om_v(i,Jend+1)/REFINED(cr)%om_v(m)
#  ifdef TIME_INTERP_FLUX
            my_value=cff1*(Wold*REFINED(cr)%DV_avg2(1,m,told)+          &
     &                     Wnew*REFINED(cr)%DV_avg2(1,m,tnew))/cff
            tl_my_value=cff1*(Wold*REFINED(cr)%tl_DV_avg2(1,m,told)+    &
     &                     Wnew*REFINED(cr)%tl_DV_avg2(1,m,tnew))/cff-  &
     &                  tl_cff*my_value/cff
#  else
            my_value=cff1*REFINED(cr)%DV_avg2(1,m,tnew)/cff
            tl_my_value=cff1*REFINED(cr)%tl_DV_avg2(1,m,tnew)/cff-      &
     &                  tl_cff*my_value/cff
#  endif
#  ifdef MASKING
            my_value=my_value*GRID(ng)%vmask(i,Jend+1)
            tl_my_value=tl_my_value*GRID(ng)%vmask(i,Jend+1)
#  endif
#  ifdef WET_DRY
            my_value=my_value*GRID(ng)%vmask_wet(i,Jend+1)
            tl_my_value=tl_my_value*GRID(ng)%vmask_wet(i,Jend+1)
#  endif
#  ifdef NESTING_DEBUG
!^          BRY_CONTACT(inorth,cr)%Mflux(i)=cff*my_value
!^
            BRY_CONTACT(inorth,cr)%tl_Mflux(i)=tl_cff*my_value+         &
     &                                         cff*tl_my_value
#  endif
!^          OCEAN(ng)%vbar(i,Jend+1,indx1(ng))=my_value
!^
            OCEAN(ng)%tl_vbar(i,Jend+1,indx1(ng))=tl_my_value
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
!^      CALL mp_exchange2d (ng, tile, model,                            &
#  ifdef SOLVE3D
!^   &                      4,                                          &
#  else
!^   &                      3,                                          &
#  endif
!^   &                      LBi, UBi, LBj, UBj,                         &
!^   &                      NghostPoints,                               &
!^   &                      EWperiodic(ng), NSperiodic(ng),             &
#  ifdef SOLVE3D
!^   &                      COUPLING(ng)%Zt_avg1,                       &
#  endif
!^   &                      OCEAN(ng)%zeta(:,:,1),                      &
!^   &                      OCEAN(ng)%zeta(:,:,2),                      &
!^   &                      OCEAN(ng)%zeta(:,:,3))
!^
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
     &                      COUPLING(ng)%tl_Zt_avg1,                    &
#  endif
     &                      OCEAN(ng)%tl_zeta(:,:,1),                   &
     &                      OCEAN(ng)%tl_zeta(:,:,2),                   &
     &                      OCEAN(ng)%tl_zeta(:,:,3))
!
!  2D momentum.
!
      ELSE
!^      CALL mp_exchange2d (ng, tile, model, 3,                         &
!^   &                      LBi, UBi, LBj, UBj,                         &
!^   &                      NghostPoints,                               &
!^   &                      EWperiodic(ng), NSperiodic(ng),             &
!^   &                      OCEAN(ng)%ubar(:,:,1),                      &
!^   &                      OCEAN(ng)%ubar(:,:,2),                      &
!^   &                      OCEAN(ng)%ubar(:,:,3))
!^
        CALL mp_exchange2d (ng, tile, model, 3,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      OCEAN(ng)%tl_ubar(:,:,1),                   &
     &                      OCEAN(ng)%tl_ubar(:,:,2),                   &
     &                      OCEAN(ng)%tl_ubar(:,:,3))

!^      CALL mp_exchange2d (ng, tile, model, 3,                         &
!^   &                      LBi, UBi, LBj, UBj,                         &
!^   &                      NghostPoints,                               &
!^   &                      EWperiodic(ng), NSperiodic(ng),             &
!^   &                      OCEAN(ng)%vbar(:,:,1),                      &
!^   &                      OCEAN(ng)%vbar(:,:,2),                      &
!^   &                      OCEAN(ng)%vbar(:,:,3))
!^
        CALL mp_exchange2d (ng, tile, model, 3,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      OCEAN(ng)%tl_vbar(:,:,1),                   &
     &                      OCEAN(ng)%tl_vbar(:,:,2),                   &
     &                      OCEAN(ng)%tl_vbar(:,:,3))

#  ifdef NESTING_DEBUG
!
!^      CALL mp_assemble (ng, model, NptsWE, spv,                       &
!^   &                    BRY_CONTACT(iwest ,cr)%Mflux(JLB:))
!^
        CALL mp_assemble (ng, model, NptsWE, spv,                       &
     &                    BRY_CONTACT(iwest ,cr)%tl_Mflux(JLB:))
!^      CALL mp_assemble (ng, model, NptsWE, spv,                       &
!^   &                    BRY_CONTACT(ieast ,cr)%Mflux(JLB:))
!^
        CALL mp_assemble (ng, model, NptsWE, spv,                       &
     &                    BRY_CONTACT(ieast ,cr)%tl_Mflux(JLB:))
!^      CALL mp_assemble (ng, model, NptsSN, spv,                       &
!^   &                    BRY_CONTACT(isouth,cr)%Mflux(ILB:))
!^
        CALL mp_assemble (ng, model, NptsSN, spv,                       &
     &                    BRY_CONTACT(isouth,cr)%tl_Mflux(ILB:))
!^      CALL mp_assemble (ng, model, NptsSN, spv,                       &
!^   &                    BRY_CONTACT(inorth,cr)%Mflux(ILB:))
!^
        CALL mp_assemble (ng, model, NptsSN, spv,                       &
     &                    BRY_CONTACT(inorth,cr)%tl_Mflux(ILB:))
#  endif
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
!
      RETURN
      END SUBROUTINE tl_put_refine2d

# ifdef SOLVE3D
!
      SUBROUTINE tl_put_refine3d (ng, dg, cr, model, tile,              &
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
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
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
!
# endif
      integer :: i, itrc, j, k, m, tnew, told
!
      real(dp) :: Wnew, Wold, SecScale, fac
      real(r8) :: my_value, tl_my_value
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", tl_put_refine3d"

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
      IF (((Wold*Wnew).lt.0.0_dp).or.((Wold+Wnew).le.0.0_dp)) THEN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF

# ifdef NESTING_DEBUG
!
!  If debugging, write information into Fortran unit 201 to check the
!  logic of interpolating from donor grid data.
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        IF (Master) THEN
          IF (first) THEN
            first=.FALSE.
            WRITE (201,20)
          END IF
          WRITE (201,30) cr, dg, ng, iic(dg), iic(ng), told,  tnew,     &
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
          CALL my_flush (201)
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
!^            my_value=Wold*                                            &
!^   &                 (Rcontact(cr)%Lweight(1,m)*                      &
!^   &                           REFINED(cr)%t(1,k,m,told,itrc)+        &
!^   &                  Rcontact(cr)%Lweight(2,m)*                      &
!^   &                           REFINED(cr)%t(2,k,m,told,itrc)+        &
!^   &                  Rcontact(cr)%Lweight(3,m)*                      &
!^   &                           REFINED(cr)%t(3,k,m,told,itrc)+        &
!^   &                  Rcontact(cr)%Lweight(4,m)*                      &
!^   &                           REFINED(cr)%t(4,k,m,told,itrc))+       &
!^   &                 Wnew*                                            &
!^   &                 (Rcontact(cr)%Lweight(1,m)*                      &
!^   &                           REFINED(cr)%t(1,k,m,tnew,itrc)+        &
!^   &                  Rcontact(cr)%Lweight(2,m)*                      &
!^   &                           REFINED(cr)%t(2,k,m,tnew,itrc)+        &
!^   &                  Rcontact(cr)%Lweight(3,m)*                      &
!^   &                           REFINED(cr)%t(3,k,m,tnew,itrc)+        &
!^   &                  Rcontact(cr)%Lweight(4,m)*                      &
!^   &                           REFINED(cr)%t(4,k,m,tnew,itrc))
!^
              tl_my_value=Wold*                                         &
     &                    (Rcontact(cr)%Lweight(1,m)*                   &
     &                              REFINED(cr)%tl_t(1,k,m,told,itrc)+  &
     &                     Rcontact(cr)%Lweight(2,m)*                   &
     &                              REFINED(cr)%tl_t(2,k,m,told,itrc)+  &
     &                     Rcontact(cr)%Lweight(3,m)*                   &
     &                              REFINED(cr)%tl_t(3,k,m,told,itrc)+  &
     &                     Rcontact(cr)%Lweight(4,m)*                   &
     &                              REFINED(cr)%tl_t(4,k,m,told,itrc))+ &
     &                    Wnew*                                         &
     &                    (Rcontact(cr)%Lweight(1,m)*                   &
     &                              REFINED(cr)%tl_t(1,k,m,tnew,itrc)+  &
     &                     Rcontact(cr)%Lweight(2,m)*                   &
     &                              REFINED(cr)%tl_t(2,k,m,tnew,itrc)+  &
     &                     Rcontact(cr)%Lweight(3,m)*                   &
     &                              REFINED(cr)%tl_t(3,k,m,tnew,itrc)+  &
     &                     Rcontact(cr)%Lweight(4,m)*                   &
     &                              REFINED(cr)%tl_t(4,k,m,tnew,itrc))
#  ifdef MASKING
!^            my_value=my_value*GRID(ng)%rmask(i,j)
!^
              tl_my_value=tl_my_value*GRID(ng)%rmask(i,j)
#  endif
!^            OCEAN(ng)%t(i,j,k,1,itrc)=my_value
!^
              OCEAN(ng)%tl_t(i,j,k,1,itrc)=tl_my_value
!^            OCEAN(ng)%t(i,j,k,2,itrc)=my_value
!^
              OCEAN(ng)%tl_t(i,j,k,2,itrc)=tl_my_value
!^            OCEAN(ng)%t(i,j,k,3,itrc)=my_value
!^
              OCEAN(ng)%tl_t(i,j,k,3,itrc)=tl_my_value
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
!^          my_value=Wold*                                              &
!^   &               (Ucontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%u(1,k,m,told)+               &
!^   &                Ucontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%u(2,k,m,told)+               &
!^   &                Ucontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%u(3,k,m,told)+               &
!^   &                Ucontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%u(4,k,m,told))+              &
!^   &               Wnew*                                              &
!^   &               (Ucontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%u(1,k,m,tnew)+               &
!^   &                Ucontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%u(2,k,m,tnew)+               &
!^   &                Ucontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%u(3,k,m,tnew)+               &
!^   &                Ucontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%u(4,k,m,tnew))
!^
            tl_my_value=Wold*                                           &
     &               (Ucontact(cr)%Lweight(1,m)*                        &
     &                         REFINED(cr)%tl_u(1,k,m,told)+            &
     &                Ucontact(cr)%Lweight(2,m)*                        &
     &                         REFINED(cr)%tl_u(2,k,m,told)+            &
     &                Ucontact(cr)%Lweight(3,m)*                        &
     &                         REFINED(cr)%tl_u(3,k,m,told)+            &
     &                Ucontact(cr)%Lweight(4,m)*                        &
     &                         REFINED(cr)%tl_u(4,k,m,told))+           &
     &               Wnew*                                              &
     &               (Ucontact(cr)%Lweight(1,m)*                        &
     &                         REFINED(cr)%tl_u(1,k,m,tnew)+            &
     &                Ucontact(cr)%Lweight(2,m)*                        &
     &                         REFINED(cr)%tl_u(2,k,m,tnew)+            &
     &                Ucontact(cr)%Lweight(3,m)*                        &
     &                         REFINED(cr)%tl_u(3,k,m,tnew)+            &
     &                Ucontact(cr)%Lweight(4,m)*                        &
     &                         REFINED(cr)%tl_u(4,k,m,tnew))
#  ifdef MASKING
!^          my_value=my_value*GRID(ng)%umask(i,j)
!^
            tl_my_value=tl_my_value*GRID(ng)%umask(i,j)
#  endif
!^          OCEAN(ng)%u(i,j,k,1)=my_value
!^
            OCEAN(ng)%tl_u(i,j,k,1)=tl_my_value
!^          OCEAN(ng)%u(i,j,k,2)=my_value
!^
            OCEAN(ng)%tl_u(i,j,k,2)=tl_my_value
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
!^          my_value=Wold*                                              &
!^   &               (Vcontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%v(1,k,m,told)+               &
!^   &                Vcontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%v(2,k,m,told)+               &
!^   &                Vcontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%v(3,k,m,told)+               &
!^   &                Vcontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%v(4,k,m,told))+              &
!^   &               Wnew*                                              &
!^   &               (Vcontact(cr)%Lweight(1,m)*                        &
!^   &                         REFINED(cr)%v(1,k,m,tnew)+               &
!^   &                Vcontact(cr)%Lweight(2,m)*                        &
!^   &                         REFINED(cr)%v(2,k,m,tnew)+               &
!^   &                Vcontact(cr)%Lweight(3,m)*                        &
!^   &                         REFINED(cr)%v(3,k,m,tnew)+               &
!^   &                Vcontact(cr)%Lweight(4,m)*                        &
!^   &                         REFINED(cr)%v(4,k,m,tnew))
!^
            tl_my_value=Wold*                                           &
     &               (Vcontact(cr)%Lweight(1,m)*                        &
     &                         REFINED(cr)%tl_v(1,k,m,told)+            &
     &                Vcontact(cr)%Lweight(2,m)*                        &
     &                         REFINED(cr)%tl_v(2,k,m,told)+            &
     &                Vcontact(cr)%Lweight(3,m)*                        &
     &                         REFINED(cr)%tl_v(3,k,m,told)+            &
     &                Vcontact(cr)%Lweight(4,m)*                        &
     &                         REFINED(cr)%tl_v(4,k,m,told))+           &
     &               Wnew*                                              &
     &               (Vcontact(cr)%Lweight(1,m)*                        &
     &                         REFINED(cr)%tl_v(1,k,m,tnew)+            &
     &                Vcontact(cr)%Lweight(2,m)*                        &
     &                         REFINED(cr)%tl_v(2,k,m,tnew)+            &
     &                Vcontact(cr)%Lweight(3,m)*                        &
     &                         REFINED(cr)%tl_v(3,k,m,tnew)+            &
     &                Vcontact(cr)%Lweight(4,m)*                        &
     &                         REFINED(cr)%tl_v(4,k,m,tnew))
#  ifdef MASKING
!^          my_value=my_value*GRID(ng)%vmask(i,j)
!^
            tl_my_value=tl_my_value*GRID(ng)%vmask(i,j)
#  endif
!^          OCEAN(ng)%v(i,j,k,1)=my_value
!^
            OCEAN(ng)%tl_v(i,j,k,1)=tl_my_value
!^          OCEAN(ng)%v(i,j,k,2)=my_value
!^
            OCEAN(ng)%tl_v(i,j,k,2)=tl_my_value
          END DO
        END IF
      END DO

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange tile information.
!-----------------------------------------------------------------------
!
!^    CALL mp_exchange4d (ng, tile, model, 3,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    OCEAN(ng)%t(:,:,:,1,:),                       &
!^   &                    OCEAN(ng)%t(:,:,:,2,:),                       &
!^   &                    OCEAN(ng)%t(:,:,:,3,:))
!^
      CALL mp_exchange4d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    OCEAN(ng)%tl_t(:,:,:,1,:),                    &
     &                    OCEAN(ng)%tl_t(:,:,:,2,:),                    &
     &                    OCEAN(ng)%tl_t(:,:,:,3,:))
!^    CALL mp_exchange3d (ng, tile, model, 4,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    OCEAN(ng)%u(:,:,:,1),                         &
!^   &                    OCEAN(ng)%u(:,:,:,2),                         &
!^   &                    OCEAN(ng)%v(:,:,:,1),                         &
!^   &                    OCEAN(ng)%v(:,:,:,2))
!^
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    OCEAN(ng)%tl_u(:,:,:,1),                      &
     &                    OCEAN(ng)%tl_u(:,:,:,2),                      &
     &                    OCEAN(ng)%tl_v(:,:,:,1),                      &
     &                    OCEAN(ng)%tl_v(:,:,:,2))
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
!
      RETURN
      END SUBROUTINE tl_put_refine3d
# endif

# ifdef SOLVE3D
!
      SUBROUTINE tl_z_weights (ng, model, tile)
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
      integer :: cr, dg, rg, i, j, k, m
      integer :: Idg, Jdg, Kdg, IminD, ImaxD, JminD, JmaxD
      integer :: Irg, Jrg, Krg, IminR, ImaxR, JminR, JmaxR
      integer :: Idgm1, Idgp1, Jdgm1, Jdgp1
      integer :: Npoints
#  ifdef DISTRIBUTE
      integer :: Nkpts, Nwpts, Nzpts

      integer, parameter :: ispv = 0
#  endif
!
      real(r8), parameter :: spv = 0.0_r8

      real(r8) :: Zbot, Zr, Ztop, dz, r1, r2
      real(r8) :: tl_Zbot, tl_Zr, tl_Ztop, tl_dz, tl_r1, tl_r2

      real(r8), allocatable :: Zd(:,:,:)
      real(r8), allocatable :: tl_Zd(:,:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", tl_z_weights"
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
          Rcontact(cr)%tl_Vweight(1:2,1:N(rg),1:Npoints)=0.0_r8
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
                  Rcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
                  Rcontact(cr)%Vweight(2,Krg,m)=0.0_r8
                  Rcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
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
            IF (.not.allocated(tl_Zd)) THEN
              allocate ( tl_Zd(4,N(dg),Npoints) )
            END IF
            tl_Zd=0.0_r8
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
                  tl_Zd(1,Kdg,m)=GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg)
                  Zd(2,Kdg,m)=GRID(dg)%z_r(Idgp1,Jdg  ,Kdg)
                  tl_Zd(2,Kdg,m)=GRID(dg)%tl_z_r(Idgp1,Jdg  ,Kdg)
                  Zd(3,Kdg,m)=GRID(dg)%z_r(Idgp1,Jdgp1,Kdg)
                  tl_Zd(3,Kdg,m)=GRID(dg)%tl_z_r(Idgp1,Jdgp1,Kdg)
                  Zd(4,Kdg,m)=GRID(dg)%z_r(Idg  ,Jdgp1,Kdg)
                  tl_Zd(4,Kdg,m)=GRID(dg)%tl_z_r(Idg  ,Jdgp1,Kdg)
                END IF
              END DO
            END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
            CALL mp_assemble (dg, model, Nzpts, spv, Zd)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL mp_assemble (dg, model, Nzpts, spv, tl_Zd)
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
                  tl_Ztop=Rcontact(cr)%Lweight(1,m)*tl_Zd(1,N(dg),m)+   &
     &                    Rcontact(cr)%Lweight(2,m)*tl_Zd(2,N(dg),m)+   &
     &                    Rcontact(cr)%Lweight(3,m)*tl_Zd(3,N(dg),m)+   &
     &                    Rcontact(cr)%Lweight(4,m)*tl_Zd(4,N(dg),m)
                  Zbot=Rcontact(cr)%Lweight(1,m)*Zd(1,1    ,m)+         &
     &                 Rcontact(cr)%Lweight(2,m)*Zd(2,1    ,m)+         &
     &                 Rcontact(cr)%Lweight(3,m)*Zd(3,1    ,m)+         &
     &                 Rcontact(cr)%Lweight(4,m)*Zd(4,1    ,m)
                  tl_Zbot=Rcontact(cr)%Lweight(1,m)*tl_Zd(1,1    ,m)+   &
     &                    Rcontact(cr)%Lweight(2,m)*tl_Zd(2,1    ,m)+   &
     &                    Rcontact(cr)%Lweight(3,m)*tl_Zd(3,1    ,m)+   &
     &                    Rcontact(cr)%Lweight(4,m)*tl_Zd(4,1    ,m)
                  Zr=GRID(rg)%z_r(Irg,Jrg,Krg)
                  tl_Zr=GRID(rg)%tl_z_r(Irg,Jrg,Krg)
                  IF (Zr.ge.Ztop) THEN           ! If shallower, use top
                    Rcontact(cr)%Kdg(Krg,m)=N(dg)! donor grid cell value
                    Rcontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Rcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
                    Rcontact(cr)%Vweight(2,Krg,m)=1.0_r8
                    Rcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
                  ELSE IF (Zbot.ge.Zr) THEN      ! If deeper, use bottom
                    Rcontact(cr)%Kdg(Krg,m)=1    ! donor grid cell value
                    Rcontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Rcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
                    Rcontact(cr)%Vweight(2,Krg,m)=1.0_r8
                    Rcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
                  ELSE                           ! bounded, interpolate
                    DO Kdg=N(dg),2,-1
                      Ztop=Rcontact(cr)%Lweight(1,m)*Zd(1,Kdg  ,m)+     &
     &                     Rcontact(cr)%Lweight(2,m)*Zd(2,Kdg  ,m)+     &
     &                     Rcontact(cr)%Lweight(3,m)*Zd(3,Kdg  ,m)+     &
     &                     Rcontact(cr)%Lweight(4,m)*Zd(4,Kdg  ,m)
                      tl_Ztop=Rcontact(cr)%Lweight(1,m)*                &
     &                                            tl_Zd(1,Kdg  ,m)+     &
     &                        Rcontact(cr)%Lweight(2,m)*                &
                                                  tl_Zd(2,Kdg  ,m)+     &
     &                        Rcontact(cr)%Lweight(3,m)*                &
     &                                            tl_Zd(3,Kdg  ,m)+     &
     &                        Rcontact(cr)%Lweight(4,m)*tl_Zd(4,Kdg  ,m)
                      Zbot=Rcontact(cr)%Lweight(1,m)*Zd(1,Kdg-1,m)+     &
     &                     Rcontact(cr)%Lweight(2,m)*Zd(2,Kdg-1,m)+     &
     &                     Rcontact(cr)%Lweight(3,m)*Zd(3,Kdg-1,m)+     &
     &                     Rcontact(cr)%Lweight(4,m)*Zd(4,Kdg-1,m)
                      tl_Zbot=Rcontact(cr)%Lweight(1,m)*                &
     &                                            tl_Zd(1,Kdg-1,m)+     &
     &                     Rcontact(cr)%Lweight(2,m)*                   &
     &                                            tl_Zd(2,Kdg-1,m)+     &
     &                     Rcontact(cr)%Lweight(3,m)*                   &
     &                                            tl_Zd(3,Kdg-1,m)+     &
     &                     Rcontact(cr)%Lweight(4,m)*tl_Zd(4,Kdg-1,m)
                      IF ((Ztop.gt.Zr).and.(Zr.ge.Zbot)) THEN
                        dz=Ztop-Zbot
                        tl_dz=tl_Ztop-tl_Zbot
                        r2=(Zr-Zbot)/dz
                        tl_r2=(tl_Zr-tl_Zbot)/dz-tl_dz*r2/dz
                        r1=1.0_r8-r2
                        tl_r1=-tl_r2
                        Rcontact(cr)%Kdg(Krg,m)=Kdg
                        Rcontact(cr)%Vweight(1,Krg,m)=r1
                        Rcontact(cr)%tl_Vweight(1,Krg,m)=tl_r1
                        Rcontact(cr)%Vweight(2,Krg,m)=r2
                        Rcontact(cr)%tl_Vweight(2,Krg,m)=tl_r2
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL mp_assemble (rg, model, Nwpts, spv, Rcontact(cr)%Vweight)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL mp_assemble (rg, model, Nwpts, spv,                      &
     &                      Rcontact(cr)%tl_Vweight)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Deallocate local work arrays.
!
          IF (allocated(Zd)) THEN
            deallocate (Zd)
          END IF
          IF (allocated(tl_Zd)) THEN
            deallocate (tl_Zd)
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
          Ucontact(cr)%tl_Vweight(1:2,1:N(rg),1:Npoints)=0.0_r8
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
                  Ucontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
                  Ucontact(cr)%Vweight(2,Krg,m)=0.0_r8
                  Ucontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
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
            IF (.not.allocated(tl_Zd)) THEN
              allocate (tl_Zd(4,N(dg),Npoints))
            END IF
            tl_Zd=0.0_r8
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
                  tl_Zd(1,Kdg,m)=0.5_r8*                                &
     &                            (GRID(dg)%tl_z_r(Idgm1,Jdg  ,Kdg)+    &
     &                             GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg))
                  Zd(2,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idg  ,Jdg  ,Kdg)+    &
     &                                GRID(dg)%z_r(Idgp1,Jdg  ,Kdg))
                  tl_Zd(2,Kdg,m)=0.5_r8*                                &
     &                            (GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg)+    &
     &                             GRID(dg)%tl_z_r(Idgp1,Jdg  ,Kdg))
                  Zd(3,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idg  ,Jdgp1,Kdg)+    &
     &                                GRID(dg)%z_r(Idgp1,Jdgp1,Kdg))
                  tl_Zd(3,Kdg,m)=0.5_r8*                                &
     &                            (GRID(dg)%tl_z_r(Idg  ,Jdgp1,Kdg)+    &
     &                             GRID(dg)%tl_z_r(Idgp1,Jdgp1,Kdg))
                  Zd(4,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idgm1,Jdgp1,Kdg)+    &
     &                                GRID(dg)%z_r(Idg  ,Jdgp1,Kdg))
                  tl_Zd(4,Kdg,m)=0.5_r8*                                &
     &                            (GRID(dg)%tl_z_r(Idgm1,Jdgp1,Kdg)+    &
     &                             GRID(dg)%tl_z_r(Idg  ,Jdgp1,Kdg))
                END IF
              END DO
            END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
            CALL mp_assemble (dg, model, Nzpts, spv, Zd)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL mp_assemble (dg, model, Nzpts, spv, tl_Zd)
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
                  tl_Ztop=Ucontact(cr)%Lweight(1,m)*tl_Zd(1,N(dg),m)+   &
     &                    Ucontact(cr)%Lweight(2,m)*tl_Zd(2,N(dg),m)+   &
     &                    Ucontact(cr)%Lweight(3,m)*tl_Zd(3,N(dg),m)+   &
     &                    Ucontact(cr)%Lweight(4,m)*tl_Zd(4,N(dg),m)
                  Zbot=Ucontact(cr)%Lweight(1,m)*Zd(1,1    ,m)+         &
     &                 Ucontact(cr)%Lweight(2,m)*Zd(2,1    ,m)+         &
     &                 Ucontact(cr)%Lweight(3,m)*Zd(3,1    ,m)+         &
     &                 Ucontact(cr)%Lweight(4,m)*Zd(4,1    ,m)
                  tl_Zbot=Ucontact(cr)%Lweight(1,m)*tl_Zd(1,1    ,m)+   &
     &                    Ucontact(cr)%Lweight(2,m)*tl_Zd(2,1    ,m)+   &
     &                    Ucontact(cr)%Lweight(3,m)*tl_Zd(3,1    ,m)+   &
     &                    Ucontact(cr)%Lweight(4,m)*tl_Zd(4,1    ,m)
                  Zr=0.5_r8*(GRID(rg)%z_r(Irg  ,Jrg,Krg)+               &
     &                       GRID(rg)%z_r(Irg-1,Jrg,Krg))
                  tl_Zr=0.5_r8*(GRID(rg)%tl_z_r(Irg  ,Jrg,Krg)+         &
     &                          GRID(rg)%tl_z_r(Irg-1,Jrg,Krg))
                  IF (Zr.ge.Ztop) THEN           ! If shallower, use top
                    Ucontact(cr)%Kdg(Krg,m)=N(dg)! donor grid cell value
                    Ucontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Ucontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
                    Ucontact(cr)%Vweight(2,Krg,m)=1.0_r8
                    Ucontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
                  ELSE IF (Zbot.ge.Zr) THEN      ! If deeper, use bottom
                    Ucontact(cr)%Kdg(Krg,m)=1    ! donor grid cell value
                    Ucontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Ucontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
                    Ucontact(cr)%Vweight(2,Krg,m)=1.0_r8
                    Ucontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
                  ELSE                           ! bounded, interpolate
                    DO Kdg=N(dg),2,-1
                      Ztop=Ucontact(cr)%Lweight(1,m)*Zd(1,Kdg  ,m)+     &
     &                     Ucontact(cr)%Lweight(2,m)*Zd(2,Kdg  ,m)+     &
     &                     Ucontact(cr)%Lweight(3,m)*Zd(3,Kdg  ,m)+     &
     &                     Ucontact(cr)%Lweight(4,m)*Zd(4,Kdg  ,m)
                      tl_Ztop=Ucontact(cr)%Lweight(1,m)*                &
     &                                            tl_Zd(1,Kdg  ,m)+     &
     &                        Ucontact(cr)%Lweight(2,m)*                &
     &                                            tl_Zd(2,Kdg  ,m)+     &
     &                        Ucontact(cr)%Lweight(3,m)*                &
     &                                            tl_Zd(3,Kdg  ,m)+     &
     &                        Ucontact(cr)%Lweight(4,m)*tl_Zd(4,Kdg  ,m)
                      Zbot=Ucontact(cr)%Lweight(1,m)*Zd(1,Kdg-1,m)+     &
     &                     Ucontact(cr)%Lweight(2,m)*Zd(2,Kdg-1,m)+     &
     &                     Ucontact(cr)%Lweight(3,m)*Zd(3,Kdg-1,m)+     &
     &                     Ucontact(cr)%Lweight(4,m)*Zd(4,Kdg-1,m)
                      tl_Zbot=Ucontact(cr)%Lweight(1,m)*                &
     &                                            tl_Zd(1,Kdg-1,m)+     &
     &                        Ucontact(cr)%Lweight(2,m)*                &
     &                                            tl_Zd(2,Kdg-1,m)+     &
     &                        Ucontact(cr)%Lweight(3,m)*                &
     &                                            tl_Zd(3,Kdg-1,m)+     &
     &                        Ucontact(cr)%Lweight(4,m)*tl_Zd(4,Kdg-1,m)
                      IF ((Ztop.gt.Zr).and.(Zr.ge.Zbot)) THEN
                        dz=Ztop-Zbot
                        tl_dz=tl_Ztop-tl_Zbot
                        r2=(Zr-Zbot)/dz
                        tl_r2=(tl_Zr-tl_Zbot)/dz-tl_dz*r2/dz
                        r1=1.0_r8-r2
                        tl_r1=-tl_r2
                        Ucontact(cr)%Kdg(Krg,m)=Kdg
                        Ucontact(cr)%Vweight(1,Krg,m)=r1
                        Ucontact(cr)%tl_Vweight(1,Krg,m)=tl_r1
                        Ucontact(cr)%Vweight(2,Krg,m)=r2
                        Ucontact(cr)%tl_Vweight(2,Krg,m)=tl_r2
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL mp_assemble (rg, model, Nwpts, spv, Ucontact(cr)%Vweight)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL mp_assemble (rg, model, Nwpts, spv,                      &
     &                      Ucontact(cr)%tl_Vweight)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Deallocate local work arrays.
!
          IF (allocated(Zd)) THEN
            deallocate (Zd)
          END IF
          IF (allocated(tl_Zd)) THEN
            deallocate (tl_Zd)
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
          Vcontact(cr)%tl_Vweight(1:2,1:N(rg),1:Npoints)=0.0_r8
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
                  Vcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
                  Vcontact(cr)%Vweight(2,Krg,m)=0.0_r8
                  Vcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
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
            IF (.not.allocated(tl_Zd)) THEN
              allocate (tl_Zd(4,N(dg),Npoints))
            END IF
            tl_Zd=0.0_r8
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
                  tl_Zd(1,Kdg,m)=0.5_r8*                                &
     &                            (GRID(dg)%tl_z_r(Idg  ,Jdgm1,Kdg)+    &
     &                             GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg))
                  Zd(2,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idgp1,Jdgm1,Kdg)+    &
     &                                GRID(dg)%z_r(Idgp1,Jdg  ,Kdg))
                  tl_Zd(2,Kdg,m)=0.5_r8*                                &
     &                            (GRID(dg)%tl_z_r(Idgp1,Jdgm1,Kdg)+    &
     &                             GRID(dg)%tl_z_r(Idgp1,Jdg  ,Kdg))
                  Zd(3,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idgp1,Jdg  ,Kdg)+    &
     &                                GRID(dg)%z_r(Idgp1,Jdgp1,Kdg))
                  tl_Zd(3,Kdg,m)=0.5_r8*                                &
     &                            (GRID(dg)%tl_z_r(Idgp1,Jdg  ,Kdg)+    &
     &                             GRID(dg)%tl_z_r(Idgp1,Jdgp1,Kdg))
                  Zd(4,Kdg,m)=0.5_r8*(GRID(dg)%z_r(Idg  ,Jdg  ,Kdg)+    &
     &                                GRID(dg)%z_r(Idg  ,Jdgp1,Kdg))
                  tl_Zd(4,Kdg,m)=0.5_r8*                                &
     &                            (GRID(dg)%tl_z_r(Idg  ,Jdg  ,Kdg)+    &
     &                             GRID(dg)%tl_z_r(Idg  ,Jdgp1,Kdg))
                END IF
              END DO
            END DO

#  ifdef DISTRIBUTE
!
!  Exchange data between all parallel nodes.
!
            CALL mp_assemble (dg, model, Nzpts, spv, Zd)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL mp_assemble (dg, model, Nzpts, spv, tl_Zd)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
                  tl_Ztop=Vcontact(cr)%Lweight(1,m)*tl_Zd(1,N(dg),m)+   &
     &                    Vcontact(cr)%Lweight(2,m)*tl_Zd(2,N(dg),m)+   &
     &                    Vcontact(cr)%Lweight(3,m)*tl_Zd(3,N(dg),m)+   &
     &                    Vcontact(cr)%Lweight(4,m)*tl_Zd(4,N(dg),m)
                  Zbot=Vcontact(cr)%Lweight(1,m)*Zd(1,1    ,m)+         &
     &                 Vcontact(cr)%Lweight(2,m)*Zd(2,1    ,m)+         &
     &                 Vcontact(cr)%Lweight(3,m)*Zd(3,1    ,m)+         &
     &                 Vcontact(cr)%Lweight(4,m)*Zd(4,1    ,m)
                  tl_Zbot=Vcontact(cr)%Lweight(1,m)*tl_Zd(1,1    ,m)+   &
     &                    Vcontact(cr)%Lweight(2,m)*tl_Zd(2,1    ,m)+   &
     &                    Vcontact(cr)%Lweight(3,m)*tl_Zd(3,1    ,m)+   &
     &                    Vcontact(cr)%Lweight(4,m)*tl_Zd(4,1    ,m)
                  Zr=0.5_r8*(GRID(rg)%z_r(Irg,Jrg  ,Krg)+               &
     &                       GRID(rg)%z_r(Irg,Jrg-1,Krg))
                  tl_Zr=0.5_r8*(GRID(rg)%tl_z_r(Irg,Jrg  ,Krg)+         &
     &                       GRID(rg)%tl_z_r(Irg,Jrg-1,Krg))
                  IF (Zr.ge.Ztop) THEN           ! If shallower, use top
                    Vcontact(cr)%Kdg(Krg,m)=N(dg)! donor grid cell value
                    Vcontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Vcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
                    Vcontact(cr)%Vweight(2,Krg,m)=1.0_r8
                    Vcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
                  ELSE IF (Zbot.ge.Zr) THEN      ! If deeper, use bottom
                    Vcontact(cr)%Kdg(Krg,m)=1    ! donor grid cell value
                    Vcontact(cr)%Vweight(1,Krg,m)=0.0_r8
                    Vcontact(cr)%tl_Vweight(1,Krg,m)=0.0_r8
                    Vcontact(cr)%Vweight(2,Krg,m)=1.0_r8
                    Vcontact(cr)%tl_Vweight(2,Krg,m)=0.0_r8
                  ELSE                           ! bounded, interpolate
                    DO Kdg=N(dg),2,-1
                      Ztop=Vcontact(cr)%Lweight(1,m)*Zd(1,Kdg  ,m)+     &
     &                     Vcontact(cr)%Lweight(2,m)*Zd(2,Kdg  ,m)+     &
     &                     Vcontact(cr)%Lweight(3,m)*Zd(3,Kdg  ,m)+     &
     &                     Vcontact(cr)%Lweight(4,m)*Zd(4,Kdg  ,m)
                      tl_Ztop=Vcontact(cr)%Lweight(1,m)*                &
     &                                            tl_Zd(1,Kdg  ,m)+     &
     &                        Vcontact(cr)%Lweight(2,m)*                &
     &                                            tl_Zd(2,Kdg  ,m)+     &
     &                        Vcontact(cr)%Lweight(3,m)*                &
     &                                            tl_Zd(3,Kdg  ,m)+     &
     &                        Vcontact(cr)%Lweight(4,m)*tl_Zd(4,Kdg  ,m)
                      Zbot=Vcontact(cr)%Lweight(1,m)*Zd(1,Kdg-1,m)+     &
     &                     Vcontact(cr)%Lweight(2,m)*Zd(2,Kdg-1,m)+     &
     &                     Vcontact(cr)%Lweight(3,m)*Zd(3,Kdg-1,m)+     &
     &                     Vcontact(cr)%Lweight(4,m)*Zd(4,Kdg-1,m)
                      tl_Zbot=Vcontact(cr)%Lweight(1,m)*                &
     &                                            tl_Zd(1,Kdg-1,m)+     &
     &                        Vcontact(cr)%Lweight(2,m)*                &
     &                                            tl_Zd(2,Kdg-1,m)+     &
     &                        Vcontact(cr)%Lweight(3,m)*                &
     &                                            tl_Zd(3,Kdg-1,m)+     &
     &                        Vcontact(cr)%Lweight(4,m)*tl_Zd(4,Kdg-1,m)
                      IF ((Ztop.gt.Zr).and.(Zr.ge.Zbot)) THEN
                        dz=Ztop-Zbot
                        tl_dz=tl_Ztop-tl_Zbot
                        r2=(Zr-Zbot)/dz
                        tl_r2=(tl_Zr-tl_Zbot)/dz-tl_dz*r2/dz
                        r1=1.0_r8-r2
                        tl_r1=-tl_r2
                        Vcontact(cr)%Kdg(Krg,m)=Kdg
                        Vcontact(cr)%Vweight(1,Krg,m)=r1
                        Vcontact(cr)%tl_Vweight(1,Krg,m)=tl_r1
                        Vcontact(cr)%Vweight(2,Krg,m)=r2
                        Vcontact(cr)%tl_Vweight(2,Krg,m)=tl_r2
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL mp_assemble (rg, model, Nwpts, spv, Vcontact(cr)%Vweight)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL mp_assemble (rg, model, Nwpts, spv,                      &
     &                                          Vcontact(cr)%tl_Vweight)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Deallocate local work arrays.
!
          IF (allocated(Zd)) THEN
            deallocate (Zd)
          END IF
          IF (allocated(tl_Zd)) THEN
            deallocate (tl_Zd)
          END IF

        END IF
      END DO
!
      RETURN
      END SUBROUTINE tl_z_weights
# endif

# ifdef SOLVE3D
!
      SUBROUTINE tl_put_contact3d (rg, model, tile,                     &
     &                          gtype, svname,                          &
     &                          cr, Npoints, contact,                   &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
#  ifdef MASKING
     &                          Amask,                                  &
#  endif
     &                          Ac, tl_Ac, tl_Ar)
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
      real(r8), intent(in) :: tl_Ac(:,:,:)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: tl_Ar(LBi:,LBj:,LBk:)
#  else
      real(r8), intent(in) :: Ac(Npoints,LBk:UBk,4)
      real(r8), intent(in) :: tl_Ac(Npoints,LBk:UBk,4)
#   ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: tl_Ar(LBi:UBi,LBj:UBj,LBk:UBk)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, k, kdg, kdgm1, m
      integer :: Istr, Iend, Jstr, Jend, Kmin

      real(r8), dimension(8) :: cff
      real(r8), dimension(8) :: tl_cff
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
            tl_cff(1)=contact(cr)%Lweight(1,m)*                         &
     &                                   contact(cr)%tl_Vweight(1,k,m)
            cff(2)=contact(cr)%Lweight(2,m)*contact(cr)%Vweight(1,k,m)
            tl_cff(2)=contact(cr)%Lweight(2,m)*                         &
     &                                   contact(cr)%tl_Vweight(1,k,m)
            cff(3)=contact(cr)%Lweight(3,m)*contact(cr)%Vweight(1,k,m)
            tl_cff(3)=contact(cr)%Lweight(3,m)*                         &
     &                                   contact(cr)%tl_Vweight(1,k,m)
            cff(4)=contact(cr)%Lweight(4,m)*contact(cr)%Vweight(1,k,m)
            tl_cff(4)=contact(cr)%Lweight(4,m)*                         &
     &                                   contact(cr)%tl_Vweight(1,k,m)
            cff(5)=contact(cr)%Lweight(1,m)*contact(cr)%Vweight(2,k,m)
            tl_cff(5)=contact(cr)%Lweight(1,m)*                         &
     &                                   contact(cr)%tl_Vweight(2,k,m)
            cff(6)=contact(cr)%Lweight(2,m)*contact(cr)%Vweight(2,k,m)
            tl_cff(6)=contact(cr)%Lweight(2,m)*                         &
     &                                   contact(cr)%tl_Vweight(2,k,m)
            cff(7)=contact(cr)%Lweight(3,m)*contact(cr)%Vweight(2,k,m)
            tl_cff(7)=contact(cr)%Lweight(3,m)*                         &
     &                                   contact(cr)%tl_Vweight(2,k,m)
            cff(8)=contact(cr)%Lweight(4,m)*contact(cr)%Vweight(2,k,m)
            tl_cff(8)=contact(cr)%Lweight(4,m)*                         &
     &                                   contact(cr)%tl_Vweight(2,k,m)
!^          Ar(i,j,k)=cff(1)*Ac(1,kdgm1,m)+                             &
!^   &                cff(2)*Ac(2,kdgm1,m)+                             &
!^   &                cff(3)*Ac(3,kdgm1,m)+                             &
!^   &                cff(4)*Ac(4,kdgm1,m)+                             &
!^   &                cff(5)*Ac(1,kdg  ,m)+                             &
!^   &                cff(6)*Ac(2,kdg  ,m)+                             &
!^   &                cff(7)*Ac(3,kdg  ,m)+                             &
!^   &                cff(8)*Ac(4,kdg  ,m)
            tl_Ar(i,j,k)=tl_cff(1)*Ac(1,kdgm1,m)+                       &
     &                cff(1)*tl_Ac(1,kdgm1,m)+                          &
     &                tl_cff(2)*Ac(2,kdgm1,m)+                          &
     &                cff(2)*tl_Ac(2,kdgm1,m)+                          &
     &                tl_cff(3)*Ac(3,kdgm1,m)+                          &
     &                cff(3)*tl_Ac(3,kdgm1,m)+                          &
     &                tl_cff(4)*Ac(4,kdgm1,m)+                          &
     &                cff(4)*tl_Ac(4,kdgm1,m)+                          &
     &                tl_cff(5)*Ac(1,kdg  ,m)+                          &
     &                cff(5)*tl_Ac(1,kdg  ,m)+                          &
     &                tl_cff(6)*Ac(2,kdg  ,m)+                          &
     &                cff(6)*tl_Ac(2,kdg  ,m)+                          &
     &                tl_cff(7)*Ac(3,kdg  ,m)+                          &
     &                cff(7)*tl_Ac(3,kdg  ,m)+                          &
     &                tl_cff(8)*Ac(4,kdg  ,m)+                          &
     &                cff(8)*tl_Ac(4,kdg  ,m)
#  ifdef MASKING
!^          Ar(i,j,k)=Ar(i,j,k)*Amask(i,j)
            tl_Ar(i,j,k)=tl_Ar(i,j,k)*Amask(i,j)
#  endif
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE tl_put_contact3d
# endif
      SUBROUTINE tl_check_massflux (ngf, model, tile)
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
      real(r8) :: tl_EastSum, tl_NorthSum, tl_SouthSum, tl_WestSum
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
          BRY_CONTACT(iwest ,cjcr)%tl_Mflux=0.0_r8
          BRY_CONTACT(ieast ,cjcr)%tl_Mflux=0.0_r8
          BRY_CONTACT(isouth,cjcr)%tl_Mflux=0.0_r8
          BRY_CONTACT(inorth,cjcr)%tl_Mflux=0.0_r8
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
              WRITE (301,10) 'Western Boundary Mass Fluxes: ',          &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
              CALL my_flush (301)
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
              tl_WestSum=0.0_r8
              Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
              DO jsum=-half,half
                Jbf=Jedge+jsum
                WestSum=WestSum+BRY_CONTACT(iwest,cr)%Mflux(Jbf)
                tl_WestSum=tl_WestSum+                                  &
     &                       BRY_CONTACT(iwest,cr)%tl_Mflux(Jbf)
              END DO
              m=BRY_CONTACT(iwest,cr)%C2Bindex(Jbf)     ! pick last one
!
!  Load coarser grid western boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
              BRY_CONTACT(iwest,cjcr)%Mflux(Jbc)=WestSum
              BRY_CONTACT(iwest,cjcr)%tl_Mflux(Jbc)=tl_WestSum

# ifdef NESTING_DEBUG
              IF (WestSum.ne.0) THEN
                MFratio=REFINED(cr)%DU_avg2(1,m,tnew)/WestSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (301,30) Jbc, REFINED(cr)%DU_avg2(1,m,tnew),        &
     &                       WestSum, MFratio
              CALL my_flush (301)
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
              WRITE (301,10) 'Eastern Boundary Mass Fluxes: ',          &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
              CALL my_flush (301)
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
              tl_EastSum=0.0_r8
              Jedge=Jo+(Jbc-Jbc_min)*RefineScale(ngf)
              DO jsum=-half,half
                Jbf=Jedge+jsum
                EastSum=EastSum+BRY_CONTACT(ieast,cr)%Mflux(Jbf)
                tl_EastSum=tl_EastSum+                                  &
     &                       BRY_CONTACT(ieast,cr)%tl_Mflux(Jbf)
              END DO
              m=BRY_CONTACT(ieast,cr)%C2Bindex(Jbf)     ! pick last one
!
!  Load coarser grid eastern boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
              BRY_CONTACT(ieast,cjcr)%Mflux(Jbc)=EastSum
              BRY_CONTACT(ieast,cjcr)%tl_Mflux(Jbc)=tl_EastSum

# ifdef NESTING_DEBUG
              IF (EastSum.ne.0) THEN
                MFratio=REFINED(cr)%DU_avg2(1,m,tnew)/EastSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (301,30) Jbc, REFINED(cr)%DU_avg2(1,m,tnew),        &
     &                       EastSum, MFratio
              CALL my_flush (301)
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
              WRITE (301,20) 'Southern Boundary Mass Fluxes: ',         &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
              CALL my_flush (301)
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
              tl_SouthSum=0.0_r8
              Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
              DO isum=-half,half
                Ibf=Iedge+isum
                SouthSum=SouthSum+BRY_CONTACT(isouth,cr)%Mflux(Ibf)
                tl_SouthSum=tl_SouthSum+                                &
     &                         BRY_CONTACT(isouth,cr)%tl_Mflux(Ibf)
              END DO
              m=BRY_CONTACT(isouth,cr)%C2Bindex(Ibf)    ! pick last one
!
!  Load coarser grid southern boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
              BRY_CONTACT(isouth,cjcr)%Mflux(Ibc)=SouthSum
              BRY_CONTACT(isouth,cjcr)%tl_Mflux(Ibc)=tl_SouthSum

# ifdef NESTING_DEBUG
              IF (SouthSum.ne.0) THEN
                MFratio=REFINED(cr)%DV_avg2(1,m,tnew)/SouthSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (301,30) Ibc, REFINED(cr)%DV_avg2(1,m,tnew),        &
     &                       SouthSum, MFratio
              CALL my_flush (301)
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
              WRITE (301,20) 'Northern Boundary Mass Fluxes: ',         &
     &                       cr, dg, rg, iif(rg), iic(rg), INT(time(rg))
              CALL my_flush (301)
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
              tl_NorthSum=0.0_r8
              Iedge=Io+(Ibc-Ibc_min)*RefineScale(ngf)
              DO isum=-half,half
                Ibf=Iedge+isum
                NorthSum=NorthSum+BRY_CONTACT(inorth,cr)%Mflux(Ibf)
                tl_NorthSum=tl_NorthSum+                                &
     &                         BRY_CONTACT(inorth,cr)%tl_Mflux(Ibf)
              END DO
              m=BRY_CONTACT(inorth,cr)%C2Bindex(Ibf)    ! pick last one
!
!  Load coarser grid northern boundary mass flux that have been averaged
!  from finer grid. These values can be compared with the coarser grid
!  values REFINED(cr)%DU_avg2 to check if the mass flux between coarser
!  and finer grid is conserved.
!
              BRY_CONTACT(inorth,cjcr)%Mflux(Ibc)=NorthSum
              BRY_CONTACT(inorth,cjcr)%tl_Mflux(Ibc)=tl_NorthSum

# ifdef NESTING_DEBUG
              IF (NorthSum.ne.0) THEN
                MFratio=REFINED(cr)%DV_avg2(1,m,tnew)/NorthSum
              ELSE
                MFratio=1.0_r8
              END IF
              WRITE (301,30) Ibc, REFINED(cr)%DV_avg2(1,m,tnew),        &
     &                       NorthSum, MFratio
# endif
            END IF
          END DO

# ifdef DISTRIBUTE
!
!  Collect data from all nodes.
!
          CALL mp_assemble (dg, model, NptsWE, spv,                     &
     &                      BRY_CONTACT(iwest ,cjcr)%Mflux(JLB:))
          CALL mp_assemble (dg, model, NptsWE, spv,                     &
     &                      BRY_CONTACT(ieast ,cjcr)%Mflux(JLB:))
          CALL mp_assemble (dg, model, NptsSN, spv,                     &
     &                      BRY_CONTACT(isouth,cjcr)%Mflux(ILB:))
          CALL mp_assemble (dg, model, NptsSN, spv,                     &
     &                      BRY_CONTACT(inorth,cjcr)%Mflux(ILB:))
          CALL mp_assemble (dg, model, NptsWE, spv,                     &
     &                      BRY_CONTACT(iwest ,cjcr)%tl_Mflux(JLB:))
          CALL mp_assemble (dg, model, NptsWE, spv,                     &
     &                      BRY_CONTACT(ieast ,cjcr)%tl_Mflux(JLB:))
          CALL mp_assemble (dg, model, NptsSN, spv,                     &
     &                      BRY_CONTACT(isouth,cjcr)%tl_Mflux(ILB:))
          CALL mp_assemble (dg, model, NptsSN, spv,                     &
     &                      BRY_CONTACT(inorth,cjcr)%tl_Mflux(ILB:))
# endif
        END IF
      END DO

# ifdef NESTING_DEBUG
!
      CALL my_flush (301)
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
!
      RETURN
      END SUBROUTINE tl_check_massflux
#endif
      END MODULE tl_nesting_mod
