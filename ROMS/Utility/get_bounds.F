#include "cppdefs.h"
      SUBROUTINE get_bounds (ng, tile, gtype, Nghost, Itile, Jtile,     &
     &                       LBi, UBi, LBj, UBj)
!
!svn $Id: get_bounds.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine compute grid bounds in the I- and J-directions.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     tile       Domain partition.                                     !
!     gtype      C-grid type. If zero, compute array allocation bounds.!
!                  Otherwise, compute bounds for IO processing.        !
!     Nghost     Number of ghost-points in the halo region:            !
!                  Nghost = 0,  compute non-overlapping bounds.        !
!                  Nghost > 0,  compute overlapping bounds.            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Itile      I-tile coordinate (a value from 0 to NtileI(ng)).     !
!     Jtile      J-tile coordinate (a value from 0 to NtileJ(ng)).     !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, gtype, Nghost
      integer, intent(out) :: Itile, Jtile, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
#ifdef DISTRIBUTE
      integer :: Iend, Istr, Jend, Jstr
      integer :: IstrM, IstrR, IstrU, IendR
      integer :: JstrM, JstrR, JstrV, JendR
      integer :: IstrB, IendB, IstrP, IendP, IstrT, IendT
      integer :: JstrB, JendB, JstrP, JendP, JstrT, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
      integer :: MyType
!
!-----------------------------------------------------------------------
!  Set array bounds in the I- and J-direction for distributed-memory
!  configurations.
!-----------------------------------------------------------------------
!
!  Set first and last grid-points according to staggered C-grid
!  classification.  If gtype = 0, it returns the values needed for
!  array allocation. Otherwise, it returns the values needed for IO
!  processing.
!
      MyType=ABS(gtype)
      IF (MyType.eq.0) THEN
# ifdef NESTING
        IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0) ) THEN
          Imin=-NghostPoints
          Imax=Lm(ng)+NghostPoints
          Jmin=-NghostPoints
          Jmax=Mm(ng)+NghostPoints
        ELSE
          IF (CompositeGrid(iwest,ng).or.EWperiodic(ng)) THEN
            Imin=-NghostPoints
          ELSE
            Imin=0
          END IF
          IF (CompositeGrid(ieast,ng).or.EWperiodic(ng)) THEN
            Imax=Lm(ng)+NghostPoints
          ELSE
            Imax=Im(ng)+1
          END IF
          IF (CompositeGrid(isouth,ng).or.NSperiodic(ng)) THEN
            Jmin=-NghostPoints
          ELSE
            Jmin=0
          END IF
          IF (CompositeGrid(inorth,ng).or.NSperiodic(ng)) THEN
            Jmax=Mm(ng)+NghostPoints
          ELSE
            Jmax=Jm(ng)+1
          END IF
        END IF
# else
        IF (EWperiodic(ng)) THEN
          IF (NSperiodic(ng)) THEN
            Imin=-NghostPoints
            Imax=Im(ng)+NghostPoints
            Jmin=-NghostPoints
            Jmax=Jm(ng)+NghostPoints
          ELSE
            Imin=-NghostPoints
            Imax=Im(ng)+NghostPoints
            Jmin=0
            Jmax=Jm(ng)+1
          END IF
        ELSE
          IF (NSperiodic(ng)) THEN
            Imin=0
            Imax=Im(ng)+1
            Jmin=-NghostPoints
            Jmax=Jm(ng)+NghostPoints
          ELSE
            Imin=0
            Imax=Im(ng)+1
            Jmin=0
            Jmax=Jm(ng)+1
          END IF
        END IF
# endif
      ELSE
        IF ((MyType.eq.p2dvar).or.(MyType.eq.u2dvar).or.                &
     &      (MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) THEN
          Imin=1
        ELSE
          Imin=0
        END IF
        Imax=Lm(ng)+1
        IF ((MyType.eq.p2dvar).or.(MyType.eq.v2dvar).or.                &
     &      (MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) THEN
          Jmin=1
        ELSE
          Jmin=0
        END IF
        Jmax=Mm(ng)+1
      END IF
!
!  Set physical, overlapping (Nghost>0) or non-overlapping (Nghost=0)
!  grid bounds according to tile rank.
!
      CALL get_tile (ng, tile,                                          &
     &               Itile, Jtile,                                      &
     &               Istr, Iend, Jstr, Jend,                            &
     &               IstrM, IstrR, IstrU, IendR,                        &
     &               JstrM, JstrR, JstrV, JendR,                        &
     &               IstrB, IendB, IstrP, IendP, IstrT, IendT,          &
     &               JstrB, JendB, JstrP, JendP, JstrT, JendT,          &
     &               Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1,          &
     &               Iendp1, Iendp2, Iendp2i, Iendp3,                   &
     &               Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1,          &
     &               Jendp1, Jendp2, Jendp2i, Jendp3)
!
      IF ((Itile.eq.-1).or.(Itile.eq.0)) THEN
        LBi=Imin
      ELSE
        LBi=Istr-Nghost
      END IF
      IF ((Itile.eq.-1).or.(Itile.eq.(NtileI(ng)-1))) THEN
        UBi=Imax
      ELSE
        UBi=Iend+Nghost
      END IF
      IF ((Jtile.eq.-1).or.(Jtile.eq.0)) THEN
        LBj=Jmin
      ELSE
        LBj=Jstr-Nghost
      END IF
      IF ((Jtile.eq.-1).or.(Jtile.eq.(NtileJ(ng)-1))) THEN
        UBj=Jmax
      ELSE
        UBj=Jend+Nghost
      END IF
#else
!
!-----------------------------------------------------------------------
!  Set array allocation bounds in the I- and J-direction for serial and
!  shared-memory configurations.
!-----------------------------------------------------------------------
!
      IF (tile.eq.-1) THEN
        Itile=-1
        Jtile=-1
      ELSE
        Jtile=tile/NtileI(ng)
        Itile=tile-Jtile*NtileI(ng)
      END IF
!
# ifdef NESTING
      IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0) ) THEN
        LBi=-NghostPoints
        UBi=Lm(ng)+NghostPoints
        LBj=-NghostPoints
        UBj=Mm(ng)+NghostPoints
      ELSE
        IF (CompositeGrid(iwest,ng).or.EWperiodic(ng)) THEN
          LBi=-NghostPoints
        ELSE
          LBi=0
        END IF
        IF (CompositeGrid(ieast,ng).or.EWperiodic(ng)) THEN
          UBi=Lm(ng)+NghostPoints
        ELSE
          UBi=Im(ng)+1
        END IF
        IF (CompositeGrid(isouth,ng).or.NSperiodic(ng)) THEN
          LBj=-NghostPoints
        ELSE
          LBj=0
        END IF
        IF (CompositeGrid(inorth,ng).or.NSperiodic(ng)) THEN
          UBj=Mm(ng)+NghostPoints
        ELSE
          UBj=Jm(ng)+1
        END IF
      END IF
# else

      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          LBi=-NghostPoints
          UBi=Im(ng)+NghostPoints
          LBj=-NghostPoints
          UBj=Jm(ng)+NghostPoints
        ELSE
          LBi=-NghostPoints
          UBi=Im(ng)+NghostPoints
          LBj=0
          UBj=Jm(ng)+1
        END IF
      ELSE
        IF (NSperiodic(ng)) THEN
          LBi=0
          UBi=Im(ng)+1
          LBj=-NghostPoints
          UBj=Jm(ng)+NghostPoints
        ELSE
          LBi=0
          UBi=Im(ng)+1
          LBj=0
          UBj=Jm(ng)+1
        END IF
      END IF
# endif
#endif

      RETURN
      END SUBROUTINE get_bounds

      SUBROUTINE get_domain (ng, tile, gtype, Nghost,                   &
     &                       epsilon, Lfullgrid,                        &
     &                       Xmin, Xmax, Ymin, Ymax)
!
!=======================================================================
!                                                                      !
!  This routine computes tile minimum and maximum fractional grid      !
!  coordinates.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     tile       Domain partition.                                     !
!     Nghost     Number of ghost-points in the halo region:            !
!                  Nghost = 0,  compute non-overlapping coordinates.   !
!                  Nghost > 0,  compute overlapping bounds.            !
!     gtype      C-grid type                                           !
!     epsilon    Small value to add to Xmax and Ymax when the tile     !
!                  is lying on the eastern and northern boundaries     !
!                  of the grid. This is usefull when processing        !
!                  observations.                                       !
!     Lfullgrid  Switch to include interior and boundaries points      !
!                  (TRUE) or just interior points (FALSE).             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Xmin       Minimum tile fractional X-coordinate.                 !
!     Xmax       Maximum tile fractional X-coordinate.                 !
!     Ymin       Minimum tile fractional Y-coordinate.                 !
!     Ymax       Maximum tile fractional Y-coordinate.                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lfullgrid

      integer, intent(in) :: ng, tile, gtype, Nghost

      real(r8), intent(in) :: epsilon
      real(r8), intent(out) :: Xmin, Xmax, Ymin, Ymax
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Itile, Jtile
!
!-----------------------------------------------------------------------
!  Computes tile minimum and maximum fractional-grid coordinates.
!-----------------------------------------------------------------------
!
      CALL get_bounds (ng, tile, gtype, Nghost, Itile, Jtile,           &
     &                 Imin, Imax, Jmin, Jmax)
!
!  Include interior and boundary points.
!
      IF (Lfullgrid) THEN
        IF ((Itile.eq.0).and.                                           &
     &      ((gtype.eq.r2dvar).or.(gtype.eq.r3dvar).or.                 &
     &       (gtype.eq.v2dvar).or.(gtype.eq.v3dvar))) THEN
          Xmin=REAL(Imin,r8)
        ELSE
          Xmin=REAL(Imin,r8)-0.5_r8
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          IF ((gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
            Xmax=REAL(Imax,r8)-0.5_r8
          ELSE
            Xmax=REAL(Imax,r8)
          END IF
        ELSE
          Xmax=REAL(Imax,r8)+0.5_r8
        END IF
        IF ((Jtile.eq.0).and.                                           &
     &      ((gtype.eq.r2dvar).or.(gtype.eq.r3dvar).or.                 &
     &       (gtype.eq.u2dvar).or.(gtype.eq.u3dvar))) THEN
          Ymin=REAL(Jmin,r8)
        ELSE
          Ymin=REAL(Jmin,r8)-0.5_r8
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          IF ((gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
            Ymax=REAL(Jmax,r8)-0.5_r8
          ELSE
            Ymax=REAL(Jmax,r8)
          END IF
        ELSE
          Ymax=REAL(Jmax,r8)+0.5_r8
        END IF
!
!   Include only interior points.
!
      ELSE
        IF (Itile.eq.0) THEN
          IF ((gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
             Xmin=REAL(Imin,r8)
          ELSE
             Xmin=REAL(Imin,r8)+0.5_r8
          END IF
        ELSE
          Xmin=REAL(Imin,r8)-0.5_r8
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          IF ((gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
            Xmax=REAL(Imax,r8)-1.0_r8
          ELSE
            Xmax=REAL(Imax,r8)-0.5_r8
          END IF
        ELSE
          Xmax=REAL(Imax,r8)+0.5_r8
        END IF
        IF (Jtile.eq.0) THEN
          IF ((gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
            Ymin=REAL(Jmin,r8)
          ELSE
            Ymin=REAL(Jmin,r8)+0.5
          END IF
        ELSE
          Ymin=REAL(Jmin,r8)-0.5_r8
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          IF ((gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
            Ymax=REAL(Jmax,r8)-1.0_r8
          ELSE
            Ymax=REAL(Jmax,r8)-0.5_r8
          END IF
        ELSE
          Ymax=REAL(Jmax,r8)+0.5_r8
        END IF
      END IF
!
!  If tile lie at the grid eastern or northen boundary, add provided
!  offset value to allow processing at those boundaries.
!
      IF (Itile.eq.(NtileI(ng)-1)) THEN
        Xmax=Xmax+epsilon
      END IF
      IF (Jtile.eq.(NtileJ(ng)-1)) THEN
        Ymax=Ymax+epsilon
      END IF

      RETURN

      END SUBROUTINE get_domain

      SUBROUTINE get_domain_edges (ng, tile,                            &
     &                             Eastern_Edge,                        &
     &                             Western_Edge,                        &
     &                             Northern_Edge,                       &
     &                             Southern_Edge,                       &
     &                             NorthEast_Corner,                    &
     &                             NorthWest_Corner,                    &
     &                             SouthEast_Corner,                    &
     &                             SouthWest_Corner,                    &
     &                             NorthEast_Test,                      &
     &                             NorthWest_Test,                      &
     &                             SouthEast_Test,                      &
     &                             SouthWest_Test)
!
!=======================================================================
!                                                                      !
!  This routine sets the logical switches (T/F) needed for processing  !
!  model variables in tiles adjacent to the domain boundary edges. It  !
!  facilitates complicated nesting configurations.                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     tile       Domain partition.                                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Eastern_Edge      tile next to the domain eastern  boundary      !
!     Western_Edge      tile next to the domain western  boundary      !
!     Northern_Edge     tile next to the domain northern boundary      !
!     Southern_Edge     tile next to the domain southern boundary      !
!                                                                      !
!     NorthEast_Corner  tile next to the domain northeastern corner    !
!     NorthWest_Corner  tile next to the domain northwestern corner    !
!     SouthEast_Corner  tile next to the domain southeastern corner    !
!     SouthWest_Corner  tile next to the domain southwestern corner    !
!                                                                      !
!     NorthEast_Test    test for tiles in the northeastern corner      !
!     NorthWest_Test    test for tiles in the northwestern corner      !
!     SouthEast_Test    test for tiles in the southeastern corner      !
!     SouthWest_Test    test for tiles in the southwestern corner      !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile

      logical, intent(out) :: Eastern_Edge
      logical, intent(out) :: Western_Edge
      logical, intent(out) :: Northern_Edge
      logical, intent(out) :: Southern_Edge
      logical, intent(out) :: NorthEast_Corner
      logical, intent(out) :: NorthWest_Corner
      logical, intent(out) :: SouthEast_Corner
      logical, intent(out) :: SouthWest_Corner
      logical, intent(out) :: NorthEast_Test
      logical, intent(out) :: NorthWest_Test
      logical, intent(out) :: SouthEast_Test
      logical, intent(out) :: SouthWest_Test
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: Itile, Jtile
!
!-----------------------------------------------------------------------
!  Compute Itile and Jtile.
!-----------------------------------------------------------------------
!
      IF (tile.eq.-1) THEN
        Itile=-1
        Jtile=-1
        Istr=1
        Iend=Lm(ng)
        Jstr=1
        Jend=Mm(ng)
      ELSE
        CALL tile_bounds_2d (ng, tile, Lm(ng), Mm(ng), Itile, Jtile,    &
     &                       Istr, Iend, Jstr, Jend)
      END IF
!
!-----------------------------------------------------------------------
!  Set the logical switches (T/F) needed for processing model variables
!  in tiles adjacent to the domain boundary edges.
!-----------------------------------------------------------------------
!
!  HGA:  Need to add the logic for composed grids.

      IF (tile.eq.-1) THEN
!
!  Set switches for the full grid (tile=-1) to TRUE, since it contains
!  all the boundary edges and corners.  This is a special case use for
!  other purposes and need only in routine "var_bounds".
!
        Western_Edge=.TRUE.
        Eastern_Edge=.TRUE.
        Southern_Edge=.TRUE.
        Northern_Edge=.TRUE.

        SouthWest_Test=.TRUE.
        SouthEast_Test=.TRUE.
        NorthWest_Test=.TRUE.
        NorthEast_Test=.TRUE.

        SouthWest_Corner=.TRUE.
        SouthEast_Corner=.TRUE.
        NorthWest_Corner=.TRUE.
        NorthEast_Corner=.TRUE.

      ELSE
!
!  Is the tile adjacent to the western or eastern domain
!  boundary edge?
!
        IF (Itile.eq.0) THEN
          Western_Edge=.TRUE.                   ! (Istr.eq.1)
        ELSE
          Western_Edge=.FALSE.                  ! elsewhere
        END IF

        IF (Itile.eq.(NtileI(ng)-1)) THEN
          Eastern_Edge=.TRUE.                   ! (Iend.eq.Lm(ng))
        ELSE
          Eastern_Edge=.FALSE.                  ! elsewhere
        END IF
!
!  Is the tile adjacent to the southern or northern domain
!  boundary edge?
!
        IF (Jtile.eq.0) THEN
          Southern_Edge=.TRUE.                  ! (Jstr.eq.1)
        ELSE
          Southern_Edge=.FALSE.                 ! elsewhere
        END IF

        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          Northern_Edge=.TRUE.                  ! (Jend.eq.Mm(ng))
        ELSE
          Northern_Edge=.FALSE.                 ! elsewhere
        END IF
!
!  Is the tile adjacent to the southwestern domain corner?
!
        IF ((Itile.eq.0).and.                                           &
     &      (Jtile.eq.0)) THEN
          SouthWest_Corner=.TRUE.               ! (Istr.eq.1).and.
          SouthWest_Test  =.TRUE.               ! (Jstr.eq.1)
        ELSE
          SouthWest_Corner=.FALSE.              ! elsewhere
#ifdef DISTRIBUTE
          SouthWest_Test  =.TRUE.
#else
          SouthWest_Test  =.FALSE.
#endif
        END IF
!
!  Is the tile adjacent to the southeastern domain corner?
!
        IF ((Itile.eq.(NtileI(ng)-1)).and.                              &
     &      (Jtile.eq.0)) THEN
          SouthEast_Corner=.TRUE.               ! (Iend.eq.Lm(ng)).and.
          SouthEast_Test  =.TRUE.               ! (Jstr.eq.1)
        ELSE
          SouthEast_Corner=.FALSE.              ! elsewhere
#ifdef DISTRIBUTE
          SouthEast_Test  =.TRUE.
#else
          SouthEast_Test  =.FALSE.
#endif
        END IF
!
!  Is the tile adjacent to the northwestern domain corner?
!
        IF ((Itile.eq.0).and.                                           &
     &      (Jtile.eq.(NtileJ(ng)-1))) THEN
          NorthWest_Corner=.TRUE.               ! (Istr.eq.1).and.
          NorthWest_Test  =.TRUE.               ! (Jend.eq.Mm(ng))
        ELSE
          NorthWest_Corner=.FALSE.              ! elsewhere
#ifdef DISTRIBUTE
          NorthWest_Test  =.TRUE.
#else
          NorthWest_Test  =.FALSE.
#endif
        END IF
!
!  Is the tile adjacent to the northeastern domain corner?
!
        IF ((Itile.eq.(NtileI(ng)-1)).and.                              &
     &      (Jtile.eq.(NtileJ(ng)-1))) THEN
          NorthEast_Corner=.TRUE.               ! (Iend.eq.Lm(ng)).and.
          NorthEast_Test  =.TRUE.               ! (Jend.eq.Mm(ng))
        ELSE
          NorthEast_Corner=.FALSE.              ! elsewhere
#ifdef DISTRIBUTE
          NorthEast_Test  =.TRUE.
#else
          NorthEast_Test  =.FALSE.
#endif
        END IF
      END IF

      RETURN
      END SUBROUTINE get_domain_edges

      SUBROUTINE get_iobounds (ng)
!
!=======================================================================
!                                                                      !
!  This routine computes the horizontal lower bound, upper bound, and  !
!  grid size for IO (NetCDF) variables.  Nested grids require special  !
!  attention due to their connetivity.                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!                                                                      !
!  On Output, the horizontal lower/upper bounds and grid size for      !
!  each variable type and nested grid number  are loaded into the      !
!  IOBOUNDS structure which is declared in module MOD_PARAM:           !
!                                                                      !
!   IOBOUNDS(ng) % ILB_psi     I-direction lower bound (PSI)           !
!   IOBOUNDS(ng) % IUB_psi     I-direction upper bound (PSI)           !
!   IOBOUNDS(ng) % JLB_psi     J-direction lower bound (PSI)           !
!   IOBOUNDS(ng) % JUB_psi     J-direction upper bound (PSI)           !
!                                                                      !
!   IOBOUNDS(ng) % ILB_rho     I-direction lower bound (RHO)           !
!   IOBOUNDS(ng) % IUB_rho     I-direction upper bound (RHO)           !
!   IOBOUNDS(ng) % JLB_rho     J-direction lower bound (RHO)           !
!   IOBOUNDS(ng) % JUB_rho     J-direction upper bound (RHO)           !
!                                                                      !
!   IOBOUNDS(ng) % ILB_u       I-direction lower bound (U)             !
!   IOBOUNDS(ng) % IUB_u       I-direction upper bound (U)             !
!   IOBOUNDS(ng) % JLB_u       J-direction lower bound (U)             !
!   IOBOUNDS(ng) % JUB_u       J-direction upper bound (U)             !
!                                                                      !
!   IOBOUNDS(ng) % ILB_v       I-direction lower bound (V)             !
!   IOBOUNDS(ng) % IUB_v       I-direction upper bound (V)             !
!   IOBOUNDS(ng) % JLB_v       J-direction lower bound (V)             !
!   IOBOUNDS(ng) % JUB_v       J-direction upper bound (V)             !
!                                                                      !
!   IOBOUNDS(ng) % xi_psi      Number of I-direction points (PSI)      !
!   IOBOUNDS(ng) % xi_rho      Number of I-direction points (RHO)      !
!   IOBOUNDS(ng) % xi_u        Number of I-direction points (U)        !
!   IOBOUNDS(ng) % xi_v        Number of I-direction points (V)        !
!                                                                      !
!   IOBOUNDS(ng) % eta_psi     Number of J-direction points (PSI)      !
!   IOBOUNDS(ng) % eta_rho     Number of J-direction points (RHO)      !
!   IOBOUNDS(ng) % eta_u       Number of I-direction points (U)        !
!   IOBOUNDS(ng) % eta_v       Number of I-direction points (V)        !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!-----------------------------------------------------------------------
!  Set IO lower/upper bounds and grid size for each C-grid type
!  variable.
!-----------------------------------------------------------------------
!
!  Recall that in non-nested applications the horizontal range,
!  including interior and boundary points, for all variable types
!  are:
!
!    PSI-type      [xi_psi, eta_psi] = [1:Lm(ng)+1, 1:Mm(ng)+1]
!    RHO-type      [xi_rho, eta_rho] = [0:Lm(ng)+1, 0:Mm(ng)+1]
!    U-type        [xi_u,   eta_u  ] = [1:Lm(ng)+1, 0:Mm(ng)+1]
!    V-type        [xi_v,   eta_v  ] = [0:Lm(ng)+1, 1:Mm(ng)+1]
!
      IOBOUNDS(ng) % ILB_psi = 1
      IOBOUNDS(ng) % IUB_psi = Lm(ng)+1
      IOBOUNDS(ng) % JLB_psi = 1
      IOBOUNDS(ng) % JUB_psi = Mm(ng)+1
!
      IOBOUNDS(ng) % ILB_rho = 0
      IOBOUNDS(ng) % IUB_rho = Lm(ng)+1
      IOBOUNDS(ng) % JLB_rho = 0
      IOBOUNDS(ng) % JUB_rho = Mm(ng)+1
!
      IOBOUNDS(ng) % ILB_u = 1
      IOBOUNDS(ng) % IUB_u = Lm(ng)+1
      IOBOUNDS(ng) % JLB_u = 0
      IOBOUNDS(ng) % JUB_u = Mm(ng)+1
!
      IOBOUNDS(ng) % ILB_v = 0
      IOBOUNDS(ng) % IUB_v = Lm(ng)+1
      IOBOUNDS(ng) % JLB_v = 1
      IOBOUNDS(ng) % JUB_v = Mm(ng)+1
!
!  Set IO NetCDF files horizontal dimension size. Recall that NetCDF
!  does not support arrays with zero index as an array element.
!
      IOBOUNDS(ng) % IorJ    = BOUNDS(ng) % UBij -                      &
     &                         BOUNDS(ng) % LBij + 1
!
      IOBOUNDS(ng) % xi_psi  = IOBOUNDS(ng) % IUB_psi -                 &
     &                         IOBOUNDS(ng) % ILB_psi + 1
      IOBOUNDS(ng) % xi_rho  = IOBOUNDS(ng) % IUB_rho -                 &
     &                         IOBOUNDS(ng) % ILB_rho + 1
      IOBOUNDS(ng) % xi_u    = IOBOUNDS(ng) % IUB_u   -                 &
     &                         IOBOUNDS(ng) % ILB_u   + 1
      IOBOUNDS(ng) % xi_v    = IOBOUNDS(ng) % IUB_v   -                 &
     &                         IOBOUNDS(ng) % ILB_v   + 1
!
      IOBOUNDS(ng) % eta_psi = IOBOUNDS(ng) % JUB_psi -                 &
     &                         IOBOUNDS(ng) % JLB_psi + 1
      IOBOUNDS(ng) % eta_rho = IOBOUNDS(ng) % JUB_rho -                 &
     &                         IOBOUNDS(ng) % JLB_rho + 1
      IOBOUNDS(ng) % eta_u   = IOBOUNDS(ng) % JUB_u   -                 &
     &                         IOBOUNDS(ng) % JLB_u   + 1
      IOBOUNDS(ng) % eta_v   = IOBOUNDS(ng) % JUB_v   -                 &
     &                         IOBOUNDS(ng) % JLB_v   + 1

      RETURN
      END SUBROUTINE get_iobounds

      SUBROUTINE get_tile (ng, tile,                                    &
     &                     Itile, Jtile,                                &
     &                     Istr, Iend, Jstr, Jend,                      &
     &                     IstrM, IstrR, IstrU, IendR,                  &
     &                     JstrM, JstrR, JstrV, JendR,                  &
     &                     IstrB, IendB, IstrP, IendP, IstrT, IendT,    &
     &                     JstrB, JendB, JstrP, JendP, JstrT, JendT,    &
     &                     Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1,    &
     &                     Iendp1, Iendp2, Iendp2i, Iendp3,             &
     &                     Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1,    &
     &                     Jendp1, Jendp2, Jendp2i, Jendp3)
!
!=======================================================================
!                                                                      !
!  This routine computes the starting and ending horizontal indices    !
!  for each sub-domain partition or tile.                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer).                         !
!     tile       Sub-domain partition.                                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Itile      I-tile coordinate (a value from 0 to NtileI(ng)).     !
!     Jtile      J-tile coordinate (a value from 0 to NtileJ(ng)).     !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!                                                                      !
!     IstrR      Starting tile index in the I-direction (RHO-points)   !
!     IstrU      Starting tile index in the I-direction (U-points)     !
!     IendR      Ending   tile index in the I-direction (RHO_points)   !
!                                                                      !
!     JstrR      Starting tile index in the J-direction (RHO-points)   !
!     JstrV      Starting tile index in the J-direction (V-points)     !
!     JendR      Ending   tile index in the J-direction (RHO_points)   !
!                                                                      !
!     IstrB      Starting nest tile in the I-direction (RHO-, V-obc)   !
!     IendB      Ending   nest tile in the I-direction (RHO-, V-obc)   !
!     IstrM      Starting nest tile in the I-direction (PSI-, U-obc)   !
!     IstrP      Starting nest tile in the I-direction (PSI, U-points) !
!     IendP      Ending   nest tile in the I-direction (PSI)           !
!     IstrT      Starting nest tile in the I-direction (RHO-points)    !
!     IendT      Ending   nest tile in the I-direction (RHO_points)    !
!                                                                      !
!     JstrB      Starting nest tile in the J-direction (RHO-, U-obc)   !
!     JendB      Ending   nest tile in the J-direction (RHO-, U-obc)   !
!     JstrM      Starting nest tile in the J-direction (PSI-, V-obc)   !
!     JstrP      Starting nest tile in the J-direction (PSI, V-points) !
!     JendP      Ending   nest tile in the J-direction (PSI)           !
!     JstrT      Starting nest tile in the J-direction (RHO-points)    !
!     JendT      Ending   nest tile in the J-direction (RHO-points)    !
!                                                                      !
!     Istrm3     Starting private I-halo computations, Istr-3          !
!     Istrm2     Starting private I-halo computations, Istr-2          !
!     Istrm1     Starting private I-halo computations, Istr-1          !
!     IstrUm2    Starting private I-halo computations, IstrU-2         !
!     IstrUm1    Starting private I-halo computations, IstrU-1         !
!     Iendp1     Ending   private I-halo computations, Iend+1          !
!     Iendp2     Ending   private I-halo computations, Iend+2          !
!     Iendp2i    Ending   private I-halo computations, Iend+2 interior !
!     Iendp3     Ending   private I-halo computations, Iend+3          !
!                                                                      !
!     Jstrm3     Starting private J-halo computations, Jstr-3          !
!     Jstrm2     Starting private J-halo computations, Jstr-2          !
!     Jstrm1     Starting private J-halo computations, Jstr-1          !
!     JstrVm2    Starting private J-halo computations, JstrV-2         !
!     JstrVm1    Starting private J-halo computations, JstrV-1         !
!     Jendp1     Ending   private J-halo computations, Jend+1          !
!     Jendp2     Ending   private J-halo computations, Jend+2          !
!     Jendp2i    Ending   private J-halo computations, Jend+2 interior !
!     Jendp3     Ending   private J-halo computations, Jend+3          !
!                                                                      !
!======================================================================!
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(out) :: Itile, Jtile
      integer, intent(out) :: Iend, Istr, Jend, Jstr
      integer, intent(out) :: IstrM, IstrR, IstrU, IendR
      integer, intent(out) :: JstrM, JstrR, JstrV, JendR
      integer, intent(out) :: IstrB, IendB, IstrP, IendP, IstrT, IendT
      integer, intent(out) :: JstrB, JendB, JstrP, JendP, JstrT, JendT
      integer, intent(out) :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer, intent(out) :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer, intent(out) :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer, intent(out) :: Jendp1, Jendp2, Jendp2i, Jendp3
!
!  Local variable declarations.
!
      integer :: my_Istr, my_Iend, my_Jstr, my_Jend
!
!-----------------------------------------------------------------------
!  Set physical non-overlapping grid bounds according to tile rank.
!-----------------------------------------------------------------------
!
!  Non-tiled grid bounds.  This is used in serial or shared-memory
!  modes to compute values in the full grid outside of parallel
!  regions.
!
      IF (tile.eq.-1) THEN
        Itile=-1
        Jtile=-1
        my_Istr=1
        my_Iend=Lm(ng)
        my_Jstr=1
        my_Jend=Mm(ng)
!
! Tiled grids bounds.
!
      ELSE
        CALL tile_bounds_2d (ng, tile, Lm(ng), Mm(ng), Itile, Jtile,    &
     &                       my_Istr, my_Iend, my_Jstr, my_Jend)
      END IF
!
!  Compute C-staggered variables bounds from tile bounds.
!
      CALL var_bounds (ng, tile, my_Istr, my_Iend, my_Jstr, my_Jend,    &
     &                 Istr, Iend, Jstr, Jend,                          &
     &                 IstrM, IstrR, IstrU, IendR,                      &
     &                 JstrM, JstrR, JstrV, JendR,                      &
     &                 IstrB, IendB, IstrP, IendP, IstrT, IendT,        &
     &                 JstrB, JendB, JstrP, JendP, JstrT, JendT,        &
     &                 Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1,        &
     &                 Iendp1, Iendp2, Iendp2i, Iendp3,                 &
     &                 Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1,        &
     &                 Jendp1, Jendp2, Jendp2i, Jendp3)

      RETURN
      END SUBROUTINE get_tile

      SUBROUTINE tile_bounds_1d (ng, tile, Imax, Istr, Iend)
!
!=======================================================================
!                                                                      !
!  This routine computes the starting and ending indices for the 1D    !
!  decomposition between all available threads or partitions.          !
!                                                                      !
!                    1 _____________________  Imax                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number                                    !
!     tile       Thread or partition                                   !
!     Imax       Global number of points                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Istr       Starting partition index                              !
!     Iend       Ending   partition index                              !
!                                                                      !
!======================================================================!
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Imax
      integer, intent(out) :: Iend, Istr
!
!  Local variable declarations.
!
      integer :: ChunkSize, Margin, Nnodes
!
!-----------------------------------------------------------------------
!  Compute 1D decomposition starting and ending indices.
!-----------------------------------------------------------------------
!
      Nnodes=NtileI(ng)*NtileJ(ng)
      ChunkSize=(Imax+Nnodes-1)/Nnodes
      Margin=(Nnodes*ChunkSize-Imax)/2

      IF (Imax.ge.Nnodes) THEN
        Istr=1+tile*ChunkSize-Margin
        Iend=Istr+ChunkSize-1
        Istr=MAX(Istr,1)
        Iend=MIN(Iend,Imax)
      ELSE
        Istr=1
        Iend=Imax
      END IF

      RETURN
      END SUBROUTINE tile_bounds_1d

      SUBROUTINE tile_bounds_2d (ng, tile, Imax, Jmax, Itile, Jtile,    &
     &                           Istr, Iend, Jstr, Jend)
!
!=======================================================================
!                                                                      !
!  This routine computes the starting and ending horizontal indices    !
!  for each sub-domain partition or tile for a grid bounded between    !
!  (1,1) and (Imax,Jmax):                                              !
!                                                                      !
!                      _________  (Imax,Jmax)                          !
!                     |         |                                      !
!                     |         |                                      !
!                     |_________|                                      !
!                (1,1)                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number                                    !
!     tile       Sub-domain partition                                  !
!     Imax       Global number of points in the I-direction            !
!     Jmax       Global number of points in the J-direction            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Itile      I-tile coordinate (a value from 0 to NtileI(ng))      !
!     Jtile      J-tile coordinate (a value from 0 to NtileJ(ng))      !
!     Istr       Starting tile index in the I-direction                !
!     Iend       Ending   tile index in the I-direction                !
!     Jstr       Starting tile index in the J-direction                !
!     Jend       Ending   tile index in the J-direction                !
!                                                                      !
!======================================================================!
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Imax, Jmax
      integer, intent(out) :: Itile, Jtile
      integer, intent(out) :: Iend, Istr, Jend, Jstr
!
!  Local variable declarations.
!
      integer :: ChunkSizeI, ChunkSizeJ, MarginI, MarginJ
!
!-----------------------------------------------------------------------
!  Compute tile decomposition for a horizontal grid bounded between
!  (1,1) and (Imax,Jmax).
!-----------------------------------------------------------------------
!
      ChunkSizeI=(Imax+NtileI(ng)-1)/NtileI(ng)
      ChunkSizeJ=(Jmax+NtileJ(ng)-1)/NtileJ(ng)
      MarginI=(NtileI(ng)*ChunkSizeI-Imax)/2
      MarginJ=(NtileJ(ng)*ChunkSizeJ-Jmax)/2
      Jtile=tile/NtileI(ng)
      Itile=tile-Jtile*NtileI(ng)
!
!  Tile bounds in the I-direction.
!
      Istr=1+Itile*ChunkSizeI-MarginI
      Iend=Istr+ChunkSizeI-1
      Istr=MAX(Istr,1)
      Iend=MIN(Iend,Imax)
!
!  Tile bounds in the J-direction.
!
      Jstr=1+Jtile*ChunkSizeJ-MarginJ
      Jend=Jstr+ChunkSizeJ-1
      Jstr=MAX(Jstr,1)
      Jend=MIN(Jend,Jmax)

      RETURN
      END SUBROUTINE tile_bounds_2d

      SUBROUTINE var_bounds (ng, tile,                                  &
     &                       my_Istr, my_Iend, my_Jstr, my_Jend,        &
     &                       Istr,  Iend,  Jstr,  Jend,                 &
     &                       IstrM, IstrR, IstrU, IendR,                &
     &                       JstrM, JstrR, JstrV, JendR,                &
     &                       IstrB, IendB, IstrP, IendP, IstrT, IendT,  &
     &                       JstrB, JendB, JstrP, JendP, JstrT, JendT,  &
     &                       Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1,  &
     &                       Iendp1, Iendp2, Iendp2i, Iendp3,           &
     &                       Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1,  &
     &                       Jendp1, Jendp2, Jendp2i, Jendp3)
!
!=======================================================================
!                                                                      !
!  This routine computes the computational grid starting and ending    !
!  horizontal indices for each C-staggered variable in terms of the    !
!  physical grid sub-domain partition or tile.                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number                                    !
!     tile       Domain partition                                      !
!     my_Istr    Physical grid starting tile index in the I-direction  !
!     my_Iend    Physical grid ending   tile index in the I-direction  !
!     my_Jstr    Physical grid starting tile index in the J-direction  !
!     my_Jend    Physical grid ending   tile index in the J-direction  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Istr       Starting tile index in the I-direction                !
!     Iend       Ending   tile index in the I-direction                !
!     Jstr       Starting tile index in the J-direction                !
!     Jend       Ending   tile index in the J-direction                !
!                                                                      !
!     IstrR      Starting tile index in the I-direction (RHO-points)   !
!     IstrU      Starting tile index in the I-direction (U-points)     !
!     IendR      Ending   tile index in the I-direction (RHO_points)   !
!                                                                      !
!     JstrR      Starting tile index in the J-direction (RHO-points)   !
!     JstrV      Starting tile index in the J-direction (V-points)     !
!     JendR      Ending   tile index in the J-direction (RHO_points)   !
!                                                                      !
!     IstrB      Starting nest tile in the I-direction (RHO-, V-obc)   !
!     IendB      Ending   nest tile in the I-direction (RHO-, V-obc)   !
!     IstrM      Starting nest tile in the I-direction (PSI-, U-obc)   !
!     IstrP      Starting nest tile in the I-direction (PSI, U-points) !
!     IendP      Ending   nest tile in the I-direction (PSI)           !
!     IstrT      Starting nest tile in the I-direction (RHO-points)    !
!     IendT      Ending   nest tile in the I-direction (RHO_points)    !
!                                                                      !
!     JstrB      Starting nest tile in the J-direction (RHO-, U-obc)   !
!     JendB      Ending   nest tile in the J-direction (RHO-, U-obc)   !
!     JstrM      Starting nest tile in the J-direction (PSI-, V-obc)   !
!     JstrP      Starting nest tile in the J-direction (PSI, V-points) !
!     JendP      Ending   nest tile in the J-direction (PSI)           !
!     JstrT      Starting nest tile in the J-direction (RHO-points)    !
!     JendT      Ending   nest tile in the J-direction (RHO-points)    !
!                                                                      !
!     Istrm3     Starting private I-halo computations, Istr-3          !
!     Istrm2     Starting private I-halo computations, Istr-2          !
!     Istrm1     Starting private I-halo computations, Istr-1          !
!     IstrUm2    Starting private I-halo computations, IstrU-2         !
!     IstrUm1    Starting private I-halo computations, IstrU-1         !
!     Iendp1     Ending   private I-halo computations, Iend+1          !
!     Iendp2     Ending   private I-halo computations, Iend+2          !
!     Iendp2i    Ending   private I-halo computations, Iend+2 interior !
!     Iendp3     Ending   private I-halo computations, Iend+3          !
!                                                                      !
!     Jstrm3     Starting private J-halo computations, Jstr-3          !
!     Jstrm2     Starting private J-halo computations, Jstr-2          !
!     Jstrm1     Starting private J-halo computations, Jstr-1          !
!     JstrVm2    Starting private J-halo computations, JstrV-2         !
!     JstrVm1    Starting private J-halo computations, JstrV-1         !
!     Jendp1     Ending   private J-halo computations, Jend+1          !
!     Jendp2     Ending   private J-halo computations, Jend+2          !
!     Jendp2i    Ending   private J-halo computations, Jend+2 interior !
!     Jendp3     Ending   private J-halo computations, Jend+3          !
!                                                                      !
!======================================================================!
!
      USE mod_param
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: my_Istr, my_Iend, my_Jstr, my_Jend
!
      integer, intent(out) :: Iend,  Istr,  Jend,  Jstr
      integer, intent(out) :: IstrM, IstrR, IstrU, IendR
      integer, intent(out) :: JstrM, JstrR, JstrV, JendR
      integer, intent(out) :: IstrB, IendB, IstrP, IendP, IstrT, IendT
      integer, intent(out) :: JstrB, JendB, JstrP, JendP, JstrT, JendT
      integer, intent(out) :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer, intent(out) :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer, intent(out) :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer, intent(out) :: Jendp1, Jendp2, Jendp2i, Jendp3
!
!=======================================================================
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables.
!=======================================================================
!
!  ROMS uses at staggered stencil:
!
!        -------v(i,j+1,k)-------               ------W(i,j,k)-------
!        |                      |               |                   |
!     u(i,j,k)   r(i,j,k)   u(i+1,j,k)          |     r(i,j,k)      |
!        |                      |               |                   |
!        --------v(i,j,k)--------               -----W(i,j,k-1)------
!
!            horizontal stencil                   vertical stencil
!                 C-grid
!
!
!  M   r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r
!      :                                                           :
!   M  v  p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p  v
!      :  +     |     |     |     |     |     |     |     |     +  :
!  Mm  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!      :  +     |     |     |     |     |     |     |     |     +  :
!   Mm v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v
!      :  +     |     |     |     |     |     |     |     |     +  :
!      r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!      :  +     |     |     |     |     |     |     |     |     +  :
!      v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v
!      :  +     |     |     |     |     |     |     |     |     +  :
!      r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!      :  +     |     |     |     |     |     |     |     |     +  :
!      v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v
!      :  +     |     |     |     |     |     |     |     |     +  :
!  2   r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!      :  +     |     |     |     |     |     |     |     |     +  :
!   2  v  p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p--v--p  v
!      :  +     |     |     |     |     |     |     |     |     +  :
!  1   r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!      :  +     |     |     |     |     |     |     |     |     +  :
!   1  v  p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p++v++p  v
!      :                                                           :
!  0   r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r..u..r
!         1     2                                         Lm    L
!      0     1     2                                         Lm    L
!
!                           interior       Boundary Conditions
!                         computations     W     E     S     N
!
!    RH0-type variables:  [1:Lm, 1:Mm]   [0,:] [L,:] [:,0] [:,M]
!    PSI-type variables:  [2:Lm, 2:Mm]   [1,:] [L,:] [:,1] [:,M]
!      U-type variables:  [2:Lm, 1:Mm]   [1,:] [L,:] [:,0] [:,M]
!      V-type variables:  [1:Lm, 2:Mm]   [0,:] [L,:] [:,1] [:,M]
!
!  Compute derived bounds for the loop indices over a subdomain tile.
!  The extended bounds (labelled by suffix R) are designed to cover
!  also the outer grid points (outlined above with :), if the subdomain
!  tile is adjacent to the physical boundary (outlined above with +).
!  Notice that IstrR, IendR, JstrR, JendR tile bounds computed here
!  DO NOT COVER ghost points (outlined below with *) associated with
!  periodic boundaries (if any) or the computational margins of MPI
!  subdomains.
!
!           Left/Top Tile                        Right/Top Tile
!
! JendR r..u..r..u..r..u..r..u  *  *      *  *  u..r..u..r..u..r..u..r
!       : Istr             Iend                Istr             Iend :
!       v  p++v++p++v++p++v++p  *  * Jend *  *  p++v++p++v++p++v++p  v
!       :  +     |     |     |                  |     |     |     +  :
!       r  u  r  u  r  u  r  u  *  *      *  *  u  r  u  r  u  r  u  r
!       :  +     |     |     |                  |     |     |     +  :
!       v  p--v--p--v--p--v--p  *  *      *  *  p--v--p--v--p--v--p  v
!       :  +     |     |     |                  |     |     |     +  :
!       r  u  r  u  r  u  r  u  *  *      *  *  u  r  u  r  u  r  u  r
!       :  +     |     |     |                  |     |     |     +  :
!       v  p--v--p--v--p--v--p  *  * Jstr *  *  p--v--p--v--p--v--p  v
!
!       *  *  *  *  *  *  *  *  *  *      *  *  *  *  *  *  *  *  *  *
!
!       *  *  *  *  *  *  *  *  *  *      *  *  *  *  *  *  *  *  *  *
!
!     IstrR    IstrU                                               IendR
!              IstrM
!
!                                                              IstrB=Istr
!                     *  *  *  *  *  *  *  *  *  *  *          IstrM=Istr
!                               Ghost Points                   IstrP=Istr
!                     *  *  *  *  *  *  *  *  *  *  *          IstrR=Istr
!                                                              IstrT=Istr
!                     *  *  p--v--p--v--p--v--p  *  *   Jend   IstrU=Istr
!                           |     |     |     |                IendB=Iend
!     Interior        *  *  u  r  u  r  u  r  u  *  *          IendP=Iend
!     Tile                  |     |     |     |                IendR=Iend
!                     *  *  p--v--p--v--p--v--p  *  *          IendT=Iend
!                           |     |     |     |                JstrB=Jstr
!                     *  *  u  r  u  r  u  r  u  *  *          JstrM=Jstr
!                           |     |     |     |                JstrP=Jstr
!                     *  *  p--v--p--v--p--v--p  *  *   Jstr   JstrR=Jstr
!                                                              JstrT=Jstr
!                     *  *  *  *  *  *  *  *  *  *  *          JstrV=Jstr
!                                                              JendB=Jend
!                     *  *  *  *  *  *  *  *  *  *  *          JendP=Jend
!                                                              JendR=Jend
!                          Istr              Iend              JendT=Jend
!
!
!
!       *  *  *  *  *  *  *  *  *  *      *  *  *  *  *  *  *  *  *  *
!
!       *  *  *  *  *  *  *  *  *  *      *  *  *  *  *  *  *  *  *  *
!         Istr             Iend
!       v  p--v--p--v--p--v--p  *  * Jend *  *  p--v--p--v--p--v--p  v
!       :  +     |     |     |                  |     |     |     +  :
!       r  u  r  u  r  u  r  u  *  *      *  *  u  r  u  r  u  r  u  r
!       :  +     |     |     |                  |     |     |     +  :
! JstrV v  p--v--p--v--p--v--p  *  *      *  *  p--v--p--v--p--v--p  v JstrM
!       :  +     |     |     |                  |     |     |     +  :
!       r  u  r  u  r  u  r  u  *  *      *  *  u  r  u  r  u  r  u  r
!       :  +     |     |     |                  |     |     |     +  :
!       v  p++v++p++v++p++v++p  *  * Jstr *  *  p++v++p++v++p++v++p  v
!       :                                                            :
!       r..u..r..u..r..u..r..u  *  *      *  *  u..r..u..r..u..r..u..r
!
!     IstrR    IstrU                                               IendR
!              IstrM
!
!           Left/Bottom Tile                    Right/Bottom Tile
!
!
!  It also computes loop-bounds for U- and V-type variables which
!  belong to the interior of the computational domain. These are
!  labelled by suffixes U,V and they step one grid point inward from
!  the side of the subdomain adjacent to the physical boundary.
!  Conversely, for an internal subdomain which does not include a
!  segment of the physical boundary, all bounds with suffixes R,U,V
!  are set to the same values of corresponding non-suffixed bounds.
!
!  In nested grids there are additional indices to process the
!  overlap regions and ranges when calling lateral boundary
!  conditions routines (see diagrams below):
!
!    IstrT    starting overlap I-direction (RHO-points)
!    IendT    ending   overlap I-direction (RHO-points)
!    JstrT    starting overlap J-direction (RHO-points)
!    JendT    ending   overlap J-direction (RHO-points)
!
!    IstrP    starting overlap I-direction (PSI-, U-points)
!    IendP    ending   overlap I-direction (PSI-points)
!    JstrP    starting overlap J-direction (PSI-, V-points)
!    JendP    ending   overlap J-direction (PSI-points)
!
!    IstrB    starting boundary I-direction (RHO-, V-points), IstrT+1
!    IendB    ending   boundary I-direction (RHO-, V-points), IendT-1
!    JstrB    starting boundary J-direction (RHO-, U-points), JstrT+1
!    JendB    ending   boundary J-direction (RHO-, U-points), JendT-1
!
!    IstrM    starting boundary I-direction (PSI-, U-points), IstrP+1
!    JstrM    starting boundary J-direction (PSI-, V-points), JstrP+1
!
!  If not nesting, these indices are set to:
!
!    IstrT = IstrR    IstrP = Istr    IstrB = Istr    IstrM = IstrU
!    IendT = IendR    IendP = Iend    IendB = Iend
!    JstrT = JstrR    JstrP = Jstr    JstrB = Jstr    JstrM = JstrV
!    JendT = JendR    JendP = Jend    JendB = Jend
!
!  The following diagram shows the lower left corner in nested grids:
!
!                             +
!        r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!           :     :     :     +     |     |     |     |     |
!        v..p..v..p..v..p..v  p--v--p--v--p--v--p--v--p--v--p--v
!           :     :     :     +     |     |     |     |     |
!   2    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!           :     :     :     +     |     |     |     |     |
!     2  v..p..v..p..v..p..v  p--v--p--v--p--v--p--v--p--v--p--v  JstrV
!           :     :     :     +     |     |     |     |     |
!   1    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  Jstr
!           :     :     :     +     |     |     |     |     |
!     1  v..p..v..p..v..p..v  p++v++p++v++p++v++p++v++p++v++p++v++
!           :     :     :
!   0    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  JstrR
!           :     :     :     :     :     :     :     :     :
!     0  v..p..v..p..v..p..v..p..v..p..v..p..v..p..v..p..v..p..v
!           :     :     :     :     :     :     :     :     :
!  -1    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!           :     :     :     :     :     :     :     :     :
!    -1  v..p..v..p..v..p..v..p..v..p..v..p..v..p..v..p..v..p..v  JstrM
!           :     :     :     :     :     :     :     :     :
!  -2    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  JstrB
!           :     :     :     :     :     :     :     :     :
!    -2  v..p..v..p..v..p..v..p..v..p..v..p..v..p..v..p..v..p..v  JstrP
!           :     :     :     :     :     :     :     :     :
!  -3    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  JstrT
!
!          -2    -1     0     1     2
!       -3    -2    -1     0     1     2
!
!         IstrP IstrM                IstrU
!            IstrB             Istr
!      IstrT             IstrR
!
!
!  The following diagram shows the upper right in nested grids:
!
!
!                                          IendR       IendT
!                                     Iend       IendB
!                                                   IendP
!
!                                   Lm    L    L+1   L+2
!                                      Lm    L    L+1   L+2
!
!  M+2   r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r   JendT
!           :     :     :     :     :     :     :     :
!    M+1 v  p..v..p..v..p..v..p..v..p..v..p..v..p..v..p..v   JendP
!           :     :     :     :     :     :     :     :
!  M+1   r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r   JendB
!           :     :     :     :     :     :     :     :
!    M+1 v  p..v..p..v..p..v..p..v..p..v..p..v..p..v..p..v
!           :     :     :     :     :     :     :     :
!  M     r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r   JendR
!                                               :     :
!    M ++v++p++v++p++v++p++v++p++v++p++v++p  v..p..v..p..v
!           |     |     |     |     |     +     :     :
!  Mm    r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r   Jend
!           |     |     |     |     |     +     :     :
!    Mm  v--p--v--p--v--p--v--p--v--p--v--p  v..p..v..p..v
!           |     |     |     |     |     +     :     :
!        r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!           |     |     |     |     |     +     :     :
!        v--p--v--p--v--p--v--p--v--p--v--p  v..p..v..p..v
!           |     |     |     |     |     +     :     :
!        r  u  r  u  r  u  r  u  r  u  r  u  r  u  r  u  r
!                                         +
!
!-----------------------------------------------------------------------
!  Starting I-tile indices.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          Istr =my_Istr
          IstrP=my_Istr
          IstrR=my_Istr
          IstrT=IstrR
          IstrU=my_Istr
          IstrB=my_Istr
          IstrM=IstrU
#ifdef NESTING
        ELSE IF (CompositeGrid(iwest,ng)) THEN
          Istr =my_Istr
          IstrP=-NghostPoints+1
          IstrR=my_Istr
          IstrT=-NghostPoints
          IstrU=my_Istr
          IstrB=IstrT                ! boundary conditions maybe avoided
          IstrM=IstrP                ! boundary conditions maybe avoided
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Istr =my_Istr
          IstrR=my_Istr-1
          IstrU=my_Istr+1
          IstrP=-NghostPoints+1
          IstrT=-NghostPoints
          IstrB=IstrT                ! boundary conditions are avoided
          IstrM=IstrP                ! boundary conditions are avoided
#endif
        ELSE
          Istr =my_Istr
          IstrP=my_Istr
          IstrR=my_Istr-1
          IstrT=IstrR
          IstrU=my_Istr+1
          IstrB=IstrT+1
          IstrM=IstrP+1
        END IF
      ELSE
        Istr =my_Istr
        IstrP=my_Istr
        IstrR=my_Istr
        IstrT=IstrR
        IstrU=my_Istr
        IstrB=my_Istr
        IstrM=IstrU
      END IF
!
!  Special case, Istrm3: used when MAX(0,Istr-3) is needed.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          Istrm3=my_Istr-3
#ifdef NESTING
        ELSE IF (CompositeGrid(iwest,ng)) THEN
          Istrm3=my_Istr-3
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Istrm3=my_Istr-3
#endif
        ELSE
          Istrm3=MAX(0,my_Istr-3)
        END IF
      ELSE
        Istrm3=my_Istr-3
      END IF
!
!  Special case, Istrm2: used when MAX(0,Istr-2) is needed.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          Istrm2=my_Istr-2
#ifdef NESTING
        ELSE IF (CompositeGrid(iwest,ng)) THEN
          Istrm2=my_Istr-2
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Istrm2=my_Istr-2
#endif
        ELSE
          Istrm2=MAX(0,my_Istr-2)
        END IF
      ELSE
        Istrm2=my_Istr-2
      END IF
!
!  Special case, IstrUm2: used when MAX(1,IstrU-2) is needed.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          IstrUm2=IstrU-2
#ifdef NESTING
        ELSE IF (CompositeGrid(iwest,ng)) THEN
          IstrUm2=IstrU-2
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          IstrUm2=IstrU-2
#endif
        ELSE
          IstrUm2=MAX(1,IstrU-2)
        END IF
      ELSE
        IstrUm2=IstrU-2
      END IF
!
!  Special case, Istrm1: used when MAX(1,Istr-1) is needed.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          Istrm1=my_Istr-1
#ifdef NESTING
        ELSE IF (CompositeGrid(iwest,ng)) THEN
          Istrm1=my_Istr-1
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Istrm1=my_Istr-1
#endif
        ELSE
          Istrm1=MAX(1,my_Istr-1)
        END IF
      ELSE
        Istrm1=my_Istr-1
      END IF
!
!  Special case, IstrUm1: used when MAX(2,IstrU-1) is needed.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          IstrUm1=IstrU-1
#ifdef NESTING
        ELSE IF (CompositeGrid(iwest,ng)) THEN
          IstrUm1=IstrU-1
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          IstrUm1=IstrU-1
#endif
        ELSE
          IstrUm1=MAX(2,IstrU-1)
        END IF
      ELSE
        IstrUm1=IstrU-1
      END IF
!
!-----------------------------------------------------------------------
!  Ending I-tile indices.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          Iend =my_Iend
          IendR=my_Iend
          IendP=IendR                ! check this one, same as IendR?
          IendT=IendR
          IendB=my_Iend
#ifdef NESTING
        ELSE IF (CompositeGrid(ieast,ng)) THEN
          Iend =my_Iend
          IendR=my_Iend
          IendP=my_Iend+NghostPoints ! check this one, same as IendT?
          IendT=my_Iend+NghostPoints
          IendB=IendT                ! boundary conditions maybe avoided
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Iend =my_Iend
          IendR=my_Iend+1
          IendP=my_Iend+NghostPoints ! check this one, same as IendT
          IendT=my_Iend+NghostPoints
          IendB=IendT                ! boundary conditions are avoided
#endif
        ELSE
          Iend =my_Iend
          IendR=my_Iend+1
          IendP=IendR
          IendT=IendR
          IendB=IendT-1
        END IF
      ELSE
        Iend =my_Iend
        IendR=my_Iend
        IendP=IendR
        IendT=IendR
        IendB=my_Iend
      END IF
!
!  Special case, Iendp1: used when MIN(Iend+1,Lm(ng)) is needed.
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          Iendp1=my_Iend+1
#ifdef NESTING
        ELSE IF (CompositeGrid(ieast,ng)) THEN
          Iendp1=my_Iend+1
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Iendp1=my_Iend+1
#endif
        ELSE
          Iendp1=MIN(my_Iend+1,Lm(ng))
        END IF
      ELSE
        Iendp1=my_Iend+1
      END IF
!
!  Special case, Iendp2i: used when MIN(Iend+2,Lm(ng)) is needed.
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          Iendp2i=my_Iend+2
#ifdef NESTING
        ELSE IF (CompositeGrid(ieast,ng)) THEN
          Iendp2i=my_Iend+2
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Iendp2i=my_Iend+2
#endif
        ELSE
          Iendp2i=MIN(my_Iend+2,Lm(ng))
        END IF
      ELSE
        Iendp2i=my_Iend+2
      END IF
!
!  Special case, Iendp2: used when MIN(Iend+2,Lm(ng)+1) is needed.
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          Iendp2=my_Iend+2
#ifdef NESTING
        ELSE IF (CompositeGrid(ieast,ng)) THEN
          Iendp2=my_Iend+2
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Iendp2=my_Iend+2
#endif
        ELSE
          Iendp2=MIN(my_Iend+2,Lm(ng)+1)
        END IF
      ELSE
        Iendp2=my_Iend+2
      END IF
!
!  Special case, Iendp3: used when MIN(Iend+3,Lm(ng)+1) is needed.
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        IF (EWperiodic(ng)) THEN
          Iendp3=my_Iend+3
#ifdef NESTING
        ELSE IF (CompositeGrid(ieast,ng)) THEN
          Iendp3=my_Iend+3
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Iendp3=my_Iend+3
#endif
        ELSE
          Iendp3=MIN(my_Iend+3,Lm(ng)+1)
        END IF
      ELSE
        Iendp3=my_Iend+3
      END IF
!
!-----------------------------------------------------------------------
!  Starting J-tile indices.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          Jstr =my_Jstr
          JstrP=my_Jstr
          JstrR=my_Jstr
          JstrT=JstrR
          JstrV=my_Jstr
          JstrB=my_Jstr
          JstrM=JstrV
#ifdef NESTING
        ELSE IF (CompositeGrid(isouth,ng)) THEN
          Jstr =my_Jstr
          JstrP=-NghostPoints+1
          JstrR=my_Jstr
          JstrT=-NghostPoints
          JstrV=my_Jstr
          JstrB=JstrT                ! boundary conditions maybe avoided
          JstrM=JstrP                ! boundary conditions maybe avoided
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Jstr =my_Jstr
          JstrR=my_Jstr-1
          JstrV=my_Jstr+1
          JstrP=-NghostPoints+1
          JstrT=-NghostPoints
          JstrB=JstrT                ! boundary conditions are avoided
          JstrM=JstrP                ! boundary conditions are avoided
#endif
        ELSE
          Jstr =my_Jstr
          JstrP=my_Jstr
          JstrR=my_Jstr-1
          JstrT=JstrR
          JstrV=my_Jstr+1
          JstrB=JstrT+1
          JstrM=JstrP+1
        END IF
      ELSE
        Jstr =my_Jstr
        JstrP=my_Jstr
        JstrR=my_Jstr
        JstrT=JstrR
        JstrV=my_Jstr
        JstrB=my_Jstr
        JstrM=JstrV
      END IF
!
!  Special case, Jstrm3: used when MAX(0,Jstr-3) is needed.
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          Jstrm3=my_Jstr-3
#ifdef NESTING
        ELSE IF (CompositeGrid(isouth,ng)) THEN
          Jstrm3=my_Jstr-3
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Jstrm3=my_Jstr-3
#endif
        ELSE
          Jstrm3=MAX(0,my_Jstr-3)
        END IF
      ELSE
        Jstrm3=my_Jstr-3
      END IF
!
!  Special case, Jstrm2: used when MAX(0,Jstr-2) is needed.
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          Jstrm2=my_Jstr-2
#ifdef NESTING
        ELSE IF (CompositeGrid(isouth,ng)) THEN
          Jstrm2=my_Jstr-2
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Jstrm2=my_Jstr-2
#endif
        ELSE
          Jstrm2=MAX(0,my_Jstr-2)
        END IF
      ELSE
        Jstrm2=my_Jstr-2
      END IF
!
!  Special case, JstrVm2: used when MAX(1,JstrV-2) is needed.
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          JstrVm2=JstrV-2
#ifdef NESTING
        ELSE IF (CompositeGrid(isouth,ng)) THEN
          JstrVm2=JstrV-2
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          JstrVm2=JstrV-2
#endif
        ELSE
          JstrVm2=MAX(1,JstrV-2)
        END IF
      ELSE
        JstrVm2=JstrV-2
      END IF
!
!  Special case, Jstrm1: used when MAX(1,Jstr-1) is needed.
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          Jstrm1=my_Jstr-1
#ifdef NESTING
        ELSE IF (CompositeGrid(isouth,ng)) THEN
          Jstrm1=my_Jstr-1
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Jstrm1=my_Jstr-1
#endif
        ELSE
          Jstrm1=MAX(1,my_Jstr-1)
        END IF
      ELSE
        Jstrm1=my_Jstr-1
      END IF
!
!  Special case, JstrVm1: used when MAX(2,JstrV-1) is needed.
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          JstrVm1=JstrV-1
#ifdef NESTING
        ELSE IF (CompositeGrid(isouth,ng)) THEN
          JstrVm1=JstrV-1
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          JstrVm1=JstrV-1
#endif
        ELSE
          JstrVm1=MAX(2,JstrV-1)
        END IF
      ELSE
        JstrVm1=JstrV-1
      END IF
!
!-----------------------------------------------------------------------
!  Ending J-tile indices.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          Jend =my_Jend
          JendR=my_Jend
          JendP=JendR                ! check this one, same as JendR?
          JendT=JendR
          JendB=my_Jend
#ifdef NESTING
        ELSE IF (CompositeGrid(inorth,ng)) THEN
          Jend =my_Jend
          JendR=my_Jend
          JendP=my_Jend+NghostPoints ! check this one, same as JendT?
          JendT=my_Jend+NghostPoints
          JendB=JendT                ! boundary conditions maybe avoided
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Jend =my_Jend
          JendR=my_Jend+1
          JendP=my_Jend+NghostPoints ! check this one, same as JendT
          JendT=my_Jend+NghostPoints
          JendB=JendT                ! boundary conditions are avoided
#endif
        ELSE
          Jend =my_Jend
          JendR=my_Jend+1
          JendP=JendR
          JendT=JendR
          JendB=JendT-1
        END IF
      ELSE
        Jend =my_Jend
        JendR=my_Jend
        JendP=JendR
        JendP=JendR
        JendT=JendR
        JendB=my_Jend
      END IF
!
!  Special case, Jendp1: used when MIN(Jend+1,Mm(ng)) is needed.
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          Jendp1=my_Jend+1
#ifdef NESTING
        ELSE IF (CompositeGrid(inorth,ng)) THEN
          Jendp1=my_Jend+1
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Jendp1=my_Jend+1
#endif
        ELSE
          Jendp1=MIN(my_Jend+1,Mm(ng))
        END IF
      ELSE
        Jendp1=my_Jend+1
      END IF
!
!  Special case, Jendp2i: used when MIN(Jend+2,Mm(ng)) is needed.
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          Jendp2i=my_Jend+2
#ifdef NESTING
        ELSE IF (CompositeGrid(inorth,ng)) THEN
          Jendp2i=my_Jend+2
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Jendp2i=my_Jend+2
#endif
        ELSE
          Jendp2i=MIN(my_Jend+2,Mm(ng))
        END IF
      ELSE
        Jendp2i=my_Jend+2
      END IF
!
!  Special case, Jendp2: used when MIN(Jend+2,Mm(ng)+1) is needed.
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          Jendp2=my_Jend+2
#ifdef NESTING
        ELSE IF (CompositeGrid(inorth,ng)) THEN
          Jendp2=my_Jend+2
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Jendp2=my_Jend+2
#endif
        ELSE
          Jendp2=MIN(my_Jend+2,Mm(ng)+1)
        END IF
      ELSE
        Jendp2=my_Jend+2
      END IF
!
!  Special case, Jendp3: used when MIN(Jend+3,Mm(ng)+1) is needed.
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        IF (NSperiodic(ng)) THEN
          Jendp3=my_Jend+3
#ifdef NESTING
        ELSE IF (CompositeGrid(inorth,ng)) THEN
          Jendp3=my_Jend+3
        ELSE IF (RefinedGrid(ng).and.(RefineScale(ng).gt.0)) THEN
          Jendp3=my_Jend+3
#endif
        ELSE
          Jendp3=MIN(my_Jend+3,Mm(ng)+1)
        END IF
      ELSE
        Jendp3=my_Jend+3
      END IF

      RETURN
      END SUBROUTINE var_bounds
