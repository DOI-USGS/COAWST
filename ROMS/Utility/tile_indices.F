#include "cppdefs.h"
      MODULE tile_indices_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module sets the application grid(s) tile decomposition bounds, !
!  indices, and switches.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     model        Numerical kernel descriptor (integer)               !
!                                                                      !
!     my_Im        Number of global grid points in the I-direction     !
!                    for each nested grid, [1:Ngrids].                 !
!                                                                      !
!     my_Jm        Number of global grid points in the J-direction     !
!                    for each nested grid, [1:Ngrids].                 !
!                                                                      !
!     my_Lm        Number of computational points in the I-direction   !
!                    for each nested grid, [1:Ngrids].                 !
!                                                                      !
!     my_Mm        Number of computational points in the J-direction.  !
!                    for each nested grid, [1:Ngrids].                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     my_BOUNDS    Lower and upper bounds indices structure per domain !
!                    partition for all grids.                          !
!                                                                      !
!     my_DOMAIN    Domain boundary edges switches and tile minimum     !
!                    and maximum fractional grid coordinates.          !
!                                                                      !
!     my_IOBOUNDS  I/O lower and upper bounds indices structure per    !
!                    domain partition for all grids.                   !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE get_bounds_mod,  ONLY : get_bounds,                           &
     &                            get_domain,                           &
     &                            get_domain_edges,                     &
     &                            get_iobounds,                         &
     &                            get_tile
!
      implicit none
!
      PUBLIC :: tile_indices
      PUBLIC :: tile_obs_bounds
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tile_indices (model,                                   &
     &                         my_Im, my_Jm,                            &
     &                         my_Lm, my_Mm,                            &
     &                         my_BOUNDS, my_DOMAIN, my_IOBOUNDS)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
      integer, intent(in) :: my_Im(1:Ngrids), my_Jm(1:Ngrids)
      integer, intent(in) :: my_Lm(1:Ngrids), my_Mm(1:Ngrids)
!
      TYPE (T_BOUNDS),   intent(out)   :: my_BOUNDS(1:Ngrids)
      TYPE (T_DOMAIN),   intent(inout) :: my_DOMAIN(1:Ngrids)
      TYPE (T_IOBOUNDS), intent(inout) :: my_IOBOUNDS(1:Ngrids)
!
!  Local variable declarations.
!
      integer :: Itile, Jtile, Nghost
      integer :: ng, tile
!
!-----------------------------------------------------------------------
!  Set lower and upper bounds indices per domain partition for all
!  nested grids.
!-----------------------------------------------------------------------
!
!  Set boundary edge I- or J-indices for each variable type.
!
      DO ng=1,Ngrids
        my_BOUNDS(ng) % edge(iwest ,p2dvar) = 1
        my_BOUNDS(ng) % edge(iwest ,r2dvar) = 0
        my_BOUNDS(ng) % edge(iwest ,u2dvar) = 1
        my_BOUNDS(ng) % edge(iwest ,v2dvar) = 0

        my_BOUNDS(ng) % edge(ieast ,p2dvar) = my_Lm(ng)+1
        my_BOUNDS(ng) % edge(ieast ,r2dvar) = my_Lm(ng)+1
        my_BOUNDS(ng) % edge(ieast ,u2dvar) = my_Lm(ng)+1
        my_BOUNDS(ng) % edge(ieast ,v2dvar) = my_Lm(ng)+1

        my_BOUNDS(ng) % edge(isouth,p2dvar) = 1
        my_BOUNDS(ng) % edge(isouth,u2dvar) = 0
        my_BOUNDS(ng) % edge(isouth,r2dvar) = 0
        my_BOUNDS(ng) % edge(isouth,v2dvar) = 1

        my_BOUNDS(ng) % edge(inorth,p2dvar) = my_Mm(ng)+1
        my_BOUNDS(ng) % edge(inorth,r2dvar) = my_Mm(ng)+1
        my_BOUNDS(ng) % edge(inorth,u2dvar) = my_Mm(ng)+1
        my_BOUNDS(ng) % edge(inorth,v2dvar) = my_Mm(ng)+1
      END DO
!
!  Set logical switches needed when processing variables in tiles
!  adjacent to the domain boundary edges or corners.  This needs to
!  be computed first since these switches are used in "get_tile".
!
      DO ng=1,Ngrids
        DO tile=-1,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain_edges (ng, tile,                              &
     &                           my_Lm(ng), my_Mm(ng),                  &
     &                           my_DOMAIN(ng)% Eastern_Edge    (tile), &
     &                           my_DOMAIN(ng)% Western_Edge    (tile), &
     &                           my_DOMAIN(ng)% Northern_Edge   (tile), &
     &                           my_DOMAIN(ng)% Southern_Edge   (tile), &
     &                           my_DOMAIN(ng)% NorthEast_Corner(tile), &
     &                           my_DOMAIN(ng)% NorthWest_Corner(tile), &
     &                           my_DOMAIN(ng)% SouthEast_Corner(tile), &
     &                           my_DOMAIN(ng)% SouthWest_Corner(tile), &
     &                           my_DOMAIN(ng)% NorthEast_Test  (tile), &
     &                           my_DOMAIN(ng)% NorthWest_Test  (tile), &
     &                           my_DOMAIN(ng)% SouthEast_Test  (tile), &
     &                           my_DOMAIN(ng)% SouthWest_Test  (tile))
        END DO
      END DO
!
!  Set tile computational indices and arrays allocation bounds
!
      Nghost=NghostPoints
      DO ng=1,Ngrids
        my_BOUNDS(ng) % LBij = 0
        my_BOUNDS(ng) % UBij = MAX(my_Lm(ng)+1,my_Mm(ng)+1)
        DO tile=-1,NtileI(ng)*NtileJ(ng)-1
          my_BOUNDS(ng) % tile(tile) = tile
          CALL get_tile (ng, tile,                                      &
     &                   my_Lm(ng), my_Mm(ng),                          &
     &                   Itile, Jtile,                                  &
     &                   my_BOUNDS(ng)% Istr   (tile),                  &
     &                   my_BOUNDS(ng)% Iend   (tile),                  &
     &                   my_BOUNDS(ng)% Jstr   (tile),                  &
     &                   my_BOUNDS(ng)% Jend   (tile),                  &
     &                   my_BOUNDS(ng)% IstrM  (tile),                  &
     &                   my_BOUNDS(ng)% IstrR  (tile),                  &
     &                   my_BOUNDS(ng)% IstrU  (tile),                  &
     &                   my_BOUNDS(ng)% IendR  (tile),                  &
     &                   my_BOUNDS(ng)% JstrM  (tile),                  &
     &                   my_BOUNDS(ng)% JstrR  (tile),                  &
     &                   my_BOUNDS(ng)% JstrV  (tile),                  &
     &                   my_BOUNDS(ng)% JendR  (tile),                  &
     &                   my_BOUNDS(ng)% IstrB  (tile),                  &
     &                   my_BOUNDS(ng)% IendB  (tile),                  &
     &                   my_BOUNDS(ng)% IstrP  (tile),                  &
     &                   my_BOUNDS(ng)% IendP  (tile),                  &
     &                   my_BOUNDS(ng)% IstrT  (tile),                  &
     &                   my_BOUNDS(ng)% IendT  (tile),                  &
     &                   my_BOUNDS(ng)% JstrB  (tile),                  &
     &                   my_BOUNDS(ng)% JendB  (tile),                  &
     &                   my_BOUNDS(ng)% JstrP  (tile),                  &
     &                   my_BOUNDS(ng)% JendP  (tile),                  &
     &                   my_BOUNDS(ng)% JstrT  (tile),                  &
     &                   my_BOUNDS(ng)% JendT  (tile),                  &
     &                   my_BOUNDS(ng)% Istrm3 (tile),                  &
     &                   my_BOUNDS(ng)% Istrm2 (tile),                  &
     &                   my_BOUNDS(ng)% Istrm1 (tile),                  &
     &                   my_BOUNDS(ng)% IstrUm2(tile),                  &
     &                   my_BOUNDS(ng)% IstrUm1(tile),                  &
     &                   my_BOUNDS(ng)% Iendp1 (tile),                  &
     &                   my_BOUNDS(ng)% Iendp2 (tile),                  &
     &                   my_BOUNDS(ng)% Iendp2i(tile),                  &
     &                   my_BOUNDS(ng)% Iendp3 (tile),                  &
     &                   my_BOUNDS(ng)% Jstrm3 (tile),                  &
     &                   my_BOUNDS(ng)% Jstrm2 (tile),                  &
     &                   my_BOUNDS(ng)% Jstrm1 (tile),                  &
     &                   my_BOUNDS(ng)% JstrVm2(tile),                  &
     &                   my_BOUNDS(ng)% JstrVm1(tile),                  &
     &                   my_BOUNDS(ng)% Jendp1 (tile),                  &
     &                   my_BOUNDS(ng)% Jendp2 (tile),                  &
     &                   my_BOUNDS(ng)% Jendp2i(tile),                  &
     &                   my_BOUNDS(ng)% Jendp3 (tile))
!
          CALL get_bounds (ng, tile, 0, Nghost,                         &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     Itile, Jtile,                                &
     &                     my_BOUNDS(ng)% LBi(tile),                    &
     &                     my_BOUNDS(ng)% UBi(tile),                    &
     &                     my_BOUNDS(ng)% LBj(tile),                    &
     &                     my_BOUNDS(ng)% UBj(tile))
        END DO
      END DO
!
!  Set I/O processing minimum (Imin, Jmax) and maximum (Imax, Jmax)
!  indices for non-overlapping (Nghost=0) and overlapping (Nghost>0)
!  tiles for each C-grid type variable.
!
      Nghost=NghostPoints
      DO ng=1,Ngrids
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_bounds (ng, tile, p2dvar, 0     ,                    &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     Itile, Jtile,                                &
     &                     my_BOUNDS(ng)% Imin(1,0,tile),               &
     &                     my_BOUNDS(ng)% Imax(1,0,tile),               &
     &                     my_BOUNDS(ng)% Jmin(1,0,tile),               &
     &                     my_BOUNDS(ng)% Jmax(1,0,tile))
          CALL get_bounds (ng, tile, p2dvar, Nghost,                    &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     Itile, Jtile,                                &
     &                     my_BOUNDS(ng)% Imin(1,1,tile),               &
     &                     my_BOUNDS(ng)% Imax(1,1,tile),               &
     &                     my_BOUNDS(ng)% Jmin(1,1,tile),               &
     &                     my_BOUNDS(ng)% Jmax(1,1,tile))
!
          CALL get_bounds (ng, tile, r2dvar, 0     ,                    &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     Itile, Jtile,                                &
     &                     my_BOUNDS(ng)% Imin(2,0,tile),               &
     &                     my_BOUNDS(ng)% Imax(2,0,tile),               &
     &                     my_BOUNDS(ng)% Jmin(2,0,tile),               &
     &                     my_BOUNDS(ng)% Jmax(2,0,tile))
          CALL get_bounds (ng, tile, r2dvar, Nghost,                    &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     Itile, Jtile,                                &
     &                     my_BOUNDS(ng)% Imin(2,1,tile),               &
     &                     my_BOUNDS(ng)% Imax(2,1,tile),               &
     &                     my_BOUNDS(ng)% Jmin(2,1,tile),               &
     &                     my_BOUNDS(ng)% Jmax(2,1,tile))
!
          CALL get_bounds (ng, tile, u2dvar, 0     ,                    &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     Itile, Jtile,                                &
     &                     my_BOUNDS(ng)% Imin(3,0,tile),               &
     &                     my_BOUNDS(ng)% Imax(3,0,tile),               &
     &                     my_BOUNDS(ng)% Jmin(3,0,tile),               &
     &                     my_BOUNDS(ng)% Jmax(3,0,tile))
          CALL get_bounds (ng, tile, u2dvar, Nghost,                    &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     Itile, Jtile,                                &
     &                     my_BOUNDS(ng)% Imin(3,1,tile),               &
     &                     my_BOUNDS(ng)% Imax(3,1,tile),               &
     &                     my_BOUNDS(ng)% Jmin(3,1,tile),               &
     &                     my_BOUNDS(ng)% Jmax(3,1,tile))
!
          CALL get_bounds (ng, tile, v2dvar, 0     ,                    &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     Itile, Jtile,                                &
     &                     my_BOUNDS(ng)% Imin(4,0,tile),               &
     &                     my_BOUNDS(ng)% Imax(4,0,tile),               &
     &                     my_BOUNDS(ng)% Jmin(4,0,tile),               &
     &                     my_BOUNDS(ng)% Jmax(4,0,tile))
          CALL get_bounds (ng, tile, v2dvar, Nghost,                    &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     Itile, Jtile,                                &
     &                     my_BOUNDS(ng)% Imin(4,1,tile),               &
     &                     my_BOUNDS(ng)% Imax(4,1,tile),               &
     &                     my_BOUNDS(ng)% Jmin(4,1,tile),               &
     &                     my_BOUNDS(ng)% Jmax(4,1,tile))
        END DO
      END DO
!
!  Set NetCDF IO bounds.
!
      DO ng=1,Ngrids
        CALL get_iobounds (ng, my_Lm(ng), my_Mm(ng),                    &
     &                     my_BOUNDS, my_IOBOUNDS)
      END DO
!
      RETURN
      END SUBROUTINE tile_indices
!
!***********************************************************************
      SUBROUTINE tile_obs_bounds (model,                               &
     &                            my_Im, my_Jm,                        &
     &                            my_Lm, my_Mm,                        &
     &                            my_DOMAIN)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
      integer, intent(in) :: my_Im(1:Ngrids), my_Jm(1:Ngrids)
      integer, intent(in) :: my_Lm(1:Ngrids), my_Mm(1:Ngrids)
!
      TYPE (T_DOMAIN), intent(inout) :: my_DOMAIN(1:Ngrids)
!
!  Local variable declarations.
!
      integer :: Itile, Jtile, Uoff, Voff
      integer :: ng, tile
!
      real(r8), parameter :: epsilon = 1.0E-8_r8
!
!-----------------------------------------------------------------------
!  Set minimum and maximum fractional coordinates for processing
!  observations. Either the full grid or only interior points will
!  be considered.  The strategy here is to add a small value (epsilon)
!  to the eastern and northern boundary values of Xmax and Ymax so
!  observations at such boundaries locations are processed. This
!  is needed because the .lt. operator in the following conditional:
!
!     IF (...
!    &    ((Xmin.le.Xobs(iobs)).and.(Xobs(iobs).lt.Xmax)).and.          &
!    &    ((Ymin.le.Yobs(iobs)).and.(Yobs(iobs).lt.Ymax))) THEN
!-----------------------------------------------------------------------
!
!  Set RHO-points domain lower and upper bounds (integer).
!
      DO ng=1,Ngrids
#ifdef DISTRIBUTE
        CALL get_bounds (ng, MyRank, r2dvar, 0,                         &
     &                   my_Im(ng), my_Jm(ng),                          &
     &                   my_Lm(ng), my_Mm(ng),                          &
     &                   Itile, Jtile,                                  &
     &                   rILB(ng), rIUB(ng), rJLB(ng), rJUB(ng))
# ifndef FULL_GRID
        IF (Itile.eq.0) THEN
          rILB(ng)=rILB(ng)+1
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          rIUB(ng)=rIUB(ng)-1
        END IF
        IF (Jtile.eq.0) THEN
          rJLB(ng)=rJLB(ng)+1
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          rJUB(ng)=rJUB(ng)-1
        END IF
# endif
#else
# ifdef FULL_GRID
        rILB(ng)=0
        rIUB(ng)=my_Lm(ng)+1
        rJLB(ng)=0
        rJUB(ng)=my_Mm(ng)+1
# else
        rILB(ng)=1
        rIUB(ng)=my_Lm(ng)
        rJLB(ng)=1
        rJUB(ng)=my_Mm(ng)
# endif
#endif
!
!  Minimum and maximum fractional coordinates for RHO-points.
!
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain (ng, tile, r2dvar, 0,                         &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     epsilon,                                     &
#ifdef FULL_GRID
     &                     .TRUE.,                                      &
#else
     &                     .FALSE.,                                     &
#endif
     &                     my_DOMAIN(ng)% Xmin_rho(tile),               &
     &                     my_DOMAIN(ng)% Xmax_rho(tile),               &
     &                     my_DOMAIN(ng)% Ymin_rho(tile),               &
     &                     my_DOMAIN(ng)% Ymax_rho(tile))
        END DO
#ifdef DISTRIBUTE
        rXmin(ng)=my_DOMAIN(ng)%Xmin_rho(MyRank)
        rXmax(ng)=my_DOMAIN(ng)%Xmax_rho(MyRank)
        rYmin(ng)=my_DOMAIN(ng)%Ymin_rho(MyRank)
        rYmax(ng)=my_DOMAIN(ng)%Ymax_rho(MyRank)
#else
        rXmin(ng)=my_DOMAIN(ng)%Xmin_rho(0)
        rXmax(ng)=my_DOMAIN(ng)%Xmax_rho(0)
        rYmin(ng)=my_DOMAIN(ng)%Ymin_rho(0)
        rYmax(ng)=my_DOMAIN(ng)%Ymax_rho(0)
#endif
      END DO
!
!  Set U-points domain lower and upper bounds (integer).
!
      DO ng=1,Ngrids
        IF (EWperiodic(ng)) THEN
          Uoff=0
        ELSE
          Uoff=1
        END IF
#ifdef DISTRIBUTE
        CALL get_bounds (ng, MyRank, u2dvar, 0,                         &
     &                   my_Im(ng), my_Jm(ng),                          &
     &                   my_Lm(ng), my_Mm(ng),                          &
     &                   Itile, Jtile,                                  &
     &                   uILB(ng), uIUB(ng), uJLB(ng), uJUB(ng))
# ifndef FULL_GRID
        IF (Itile.eq.0) THEN
          uILB(ng)=uILB(ng)+Uoff
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          uIUB(ng)=uIUB(ng)-1
        END IF
        IF (Jtile.eq.0) THEN
          uJLB(ng)=uJLB(ng)+1
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          uJUB(ng)=uJUB(ng)-1
        END IF
# endif
#else
# ifdef FULL_GRID
        uILB(ng)=1
        uIUB(ng)=my_Lm(ng)+1
        uJLB(ng)=0
        uJUB(ng)=my_Mm(ng)+1
# else
        uILB(ng)=1+Uoff
        uIUB(ng)=my_Lm(ng)
        uJLB(ng)=1
        uJUB(ng)=my_Mm(ng)
# endif
#endif
!
!  Minimum and maximum fractional coordinates for U-points.
!
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain (ng, tile, u2dvar, 0,                         &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     epsilon,                                     &
#ifdef FULL_GRID
     &                     .TRUE.,                                      &
#else
     &                     .FALSE.,                                     &
#endif
     &                     my_DOMAIN(ng)% Xmin_u(tile),                 &
     &                     my_DOMAIN(ng)% Xmax_u(tile),                 &
     &                     my_DOMAIN(ng)% Ymin_u(tile),                 &
     &                     my_DOMAIN(ng)% Ymax_u(tile))
        END DO
#ifdef DISTRIBUTE
        uXmin(ng)=my_DOMAIN(ng)%Xmin_u(MyRank)
        uXmax(ng)=my_DOMAIN(ng)%Xmax_u(MyRank)
        uYmin(ng)=my_DOMAIN(ng)%Ymin_u(MyRank)
        uYmax(ng)=my_DOMAIN(ng)%Ymax_u(MyRank)
#else
        uXmin(ng)=my_DOMAIN(ng)%Xmin_u(0)
        uXmax(ng)=my_DOMAIN(ng)%Xmax_u(0)
        uYmin(ng)=my_DOMAIN(ng)%Ymin_u(0)
        uYmax(ng)=my_DOMAIN(ng)%Ymax_u(0)
#endif
      END DO
!
!  Set V-points domain lower and upper bounds (integer).
!
      DO ng=1,Ngrids
        IF (NSperiodic(ng)) THEN
          Voff=0
        ELSE
          Voff=1
        END IF
#ifdef DISTRIBUTE
        CALL get_bounds (ng, MyRank, v2dvar, 0,                         &
     &                   my_Im(ng), my_Jm(ng),                          &
     &                   my_Lm(ng), my_Mm(ng),                          &
     &                   Itile, Jtile,                                  &
     &                   vILB(ng), vIUB(ng), vJLB(ng), vJUB(ng))
# ifndef FULL_GRID
        IF (Itile.eq.0) THEN
          vILB(ng)=vILB(ng)+1
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          vIUB(ng)=vIUB(ng)-1
        END IF
        IF (Jtile.eq.0) THEN
          vJLB(ng)=vJLB(ng)+Voff
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          vJUB(ng)=vJUB(ng)-1
        END IF
# endif
#else
#  ifdef FULL_GRID
        vILB(ng)=0
        vIUB(ng)=Lm(ng)+1
        vJLB(ng)=1
        vJUB(ng)=Mm(ng)+1
# else
        vILB(ng)=1
        vIUB(ng)=Lm(ng)
        vJLB(ng)=1+Voff
        vJUB(ng)=Mm(ng)
# endif
#endif
!
!  Minimum and maximum fractional coordinates for V-points.
!
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain (ng, tile, v2dvar, 0,                         &
     &                     my_Im(ng), my_Jm(ng),                        &
     &                     my_Lm(ng), my_Mm(ng),                        &
     &                     epsilon,                                     &
#ifdef FULL_GRID
     &                     .TRUE.,                                      &
#else
     &                     .FALSE.,                                     &
#endif
     &                     my_DOMAIN(ng)% Xmin_v(tile),                 &
     &                     my_DOMAIN(ng)% Xmax_v(tile),                 &
     &                     my_DOMAIN(ng)% Ymin_v(tile),                 &
     &                     my_DOMAIN(ng)% Ymax_v(tile))
        END DO
#ifdef DISTRIBUTE
        vXmin(ng)=my_DOMAIN(ng)%Xmin_v(MyRank)
        vXmax(ng)=my_DOMAIN(ng)%Xmax_v(MyRank)
        vYmin(ng)=my_DOMAIN(ng)%Ymin_v(MyRank)
        vYmax(ng)=my_DOMAIN(ng)%Ymax_v(MyRank)
#else
        vXmin(ng)=my_DOMAIN(ng)%Xmin_v(0)
        vXmax(ng)=my_DOMAIN(ng)%Xmax_v(0)
        vYmin(ng)=my_DOMAIN(ng)%Ymin_v(0)
        vYmax(ng)=my_DOMAIN(ng)%Ymax_v(0)
#endif
      END DO
!
      RETURN
      END SUBROUTINE tile_obs_bounds
!
      END MODULE tile_indices_mod
