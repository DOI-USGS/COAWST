#include "cppdefs.h"
      SUBROUTINE get_bounds (ng, tile, gtype, Nghost, Itile, Jtile,     &
     &                       LBi, UBi, LBj, UBj)
!
!svn $Id: get_bounds.F 814 2008-10-29 01:42:17Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
!                  Nghost = 0,  compute non-overlaping bounds.         !
!                  Nghost > 0,  compute overlaping bounds.             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Itile      Tile coordinate in the I-direction, use only in       !
!                  distributed-memory applications.                    !
!     Jtile      Tile coordinate in the J-direction, use only in       !
!                  distributed-memory applications.                    !
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
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, gtype, Nghost
      integer, intent(out) :: Itile, Jtile, LBi, UBi, LBj, UBj
#ifdef DISTRIBUTE
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Istr, IstrR, IstrT, IstrU, Iend, IendR, IendT
      integer :: Jstr, JstrR, JstrT, JstrV, Jend, JendR, JendT
      integer :: IstrTU, JstrTV, IendTU, JendTV
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
        Imin=LOWER_BOUND_I
        Imax=UPPER_BOUND_I
        Jmin=LOWER_BOUND_J
        Jmax=UPPER_BOUND_J
      ELSE
# ifdef REFINED_GRID
        IF (ng.eq.1) THEN
          Imin=0
          Imax=Lm(ng)+1
          Jmin=0
          Jmax=Mm(ng)+1
          IF ((MyType.eq.p2dvar).or.(MyType.eq.u2dvar).or.              &
     &        (MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) Imin=1
          IF ((MyType.eq.p2dvar).or.(MyType.eq.v2dvar).or.              &
     &        (MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) Jmin=1
        ELSE
          Imin=LOWER_BOUND_I
          Imax=Lm(ng)+GHOST_POINTS
          Jmin=LOWER_BOUND_J
          Jmax=Mm(ng)+GHOST_POINTS
          IF ((MyType.eq.p2dvar).or.(MyType.eq.u2dvar).or.              &
     &        (MyType.eq.p3dvar).or.(MyType.eq.u3dvar))                 &
     &         Imin=LOWER_BOUND_I+1
          IF ((MyType.eq.p2dvar).or.(MyType.eq.v2dvar).or.              &
     &        (MyType.eq.p3dvar).or.(MyType.eq.v3dvar))                 &
     &         Jmin=LOWER_BOUND_J+1
        END IF
# else
        Imin=0
        Imax=Lm(ng)+1
        Jmin=0
        Jmax=Mm(ng)+1
        IF ((MyType.eq.p2dvar).or.(MyType.eq.u2dvar).or.                &
     &      (MyType.eq.p3dvar).or.(MyType.eq.u3dvar)) Imin=1
        IF ((MyType.eq.p2dvar).or.(MyType.eq.v2dvar).or.                &
     &      (MyType.eq.p3dvar).or.(MyType.eq.v3dvar)) Jmin=1
# endif
      END IF
!
!  Set physical, overlaping (Nghost>0) or non-overlaping (Nghost=0) 
!  grid bounds according to tile rank.
!
      CALL get_tile (ng, tile, Itile, Jtile,                            &
     &               Istr, Iend, Jstr, Jend,                            &
     &               IstrR, IstrT, IstrU, IendR, IendT,                 &
     &               JstrR, JstrT, JstrV, JendR, JendT,                 &
     &               IstrTU, JstrTV, IendTU, JendTV)
!
      IF (Itile.eq.0) THEN
        LBi=Imin
      ELSE
        LBi=Istr-Nghost
      END IF
      IF (Itile.eq.(NtileI(ng)-1)) THEN
        UBi=Imax
      ELSE
        UBi=Iend+Nghost
      END IF
      IF (Jtile.eq.0) THEN
        LBj=Jmin
      ELSE
        LBj=Jstr-Nghost
      END IF
      IF (Jtile.eq.(NtileJ(ng)-1)) THEN
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
      Itile=-1
      Jtile=-1
      LBi=LOWER_BOUND_I
      UBi=UPPER_BOUND_I
      LBj=LOWER_BOUND_J
      UBj=UPPER_BOUND_J
#endif
      RETURN
      END SUBROUTINE get_bounds

      SUBROUTINE get_tile (ng, tile, Itile, Jtile,                      &
     &                     Istr, Iend, Jstr, Jend,                      &
     &                     IstrR, IstrT, IstrU, IendR, IendT,           &
     &                     JstrR, JstrT, JstrV, JendR, JendT,           &
     &                     IstrTU, JstrTV, IendTU, JendTV)

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
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!                                                                      !
!     IstrR      Starting tile index in the I-direction (RHO-points).  !
!     IstrT      Starting nest tile  in the I-direction (RHO-points).  !
!     IstrU      Starting tile index in the I-direction (U-points).    !
!     IendR      Ending   tile index in the I-direction (RHO_points).  !
!     IendT      Ending   nest tile  in the I-direction (RHO_points).  !
!                                                                      !
!     JstrR      Starting tile index in the J-direction (RHO-points).  !
!     JstrT      Starting nest tile  in the J-direction (RHO-points).  !
!     JstrV      Starting tile index in the J-direction (V-points).    !
!     JendR      Ending   tile index in the J-direction (RHO_points).  !
!     JendT      Ending   nest tile  in the J-direction (RHO-points).  !
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
      integer, intent(out) :: IstrR, IstrT, IstrU, IendR, IendT
      integer, intent(out) :: JstrR, JstrT, JstrV, JendR, JendT
      integer, intent(out) :: IstrTU, JstrTV, IendTU, JendTV
!
!  Local variable declarations.
!
      integer :: ChunkSizeI, ChunkSizeJ, MarginI, MarginJ
!
!-----------------------------------------------------------------------
!  Set physical non-overlaping grid bounds according to tile rank.
!-----------------------------------------------------------------------
!
!  Non-tiled grid bounds.  This is used in serial or shared-memory
!  modes to compute values in the full grid outside of parallel
!  regions.
!
      IF (tile.eq.-1) THEN
        Itile=-1
        Jtile=-1
        Istr=1
        Iend=Lm(ng)
        Jstr=1
        Jend=Mm(ng)
!
! Tiled grids bounds.
!
      ELSE
        ChunkSizeI=(Lm(ng)+NtileI(ng)-1)/NtileI(ng)
        ChunkSizeJ=(Mm(ng)+NtileJ(ng)-1)/NtileJ(ng)
        MarginI=(NtileI(ng)*ChunkSizeI-Lm(ng))/2
        MarginJ=(NtileJ(ng)*ChunkSizeJ-Mm(ng))/2
        Jtile=tile/NtileI(ng)
        Itile=tile-Jtile*NtileI(ng)
!
!  Tile bounds in the I-direction.
!
        Istr=1+Itile*ChunkSizeI-MarginI
        Iend=Istr+ChunkSizeI-1
        Istr=MAX(Istr,1)
        Iend=MIN(Iend,Lm(ng))
!
!  Tile bounds in the J-direction.
!
        Jstr=1+Jtile*ChunkSizeJ-MarginJ
        Jend=Jstr+ChunkSizeJ-1
        Jstr=MAX(Jstr,1)
        Jend=MIN(Jend,Mm(ng))
      END IF
!
!  Compute C-staggered variables bounds from tile bounds.
!
      CALL var_bounds (ng, Istr, Iend, Jstr, Jend,                      &
     &                 IstrR, IstrT, IstrU, IendR, IendT,               &
     &                 JstrR, JstrT, JstrV, JendR, JendT,               &
     &                 IstrTU, JstrTV, IendTU, JendTV)

      RETURN
      END SUBROUTINE get_tile

      SUBROUTINE var_bounds (ng, Istr, Iend, Jstr, Jend,                &
     &                       IstrR, IstrT, IstrU, IendR, IendT,         &
     &                       JstrR, JstrT, JstrV, JendR, JendT,         &
     &                       IstrTU, JstrTV, IendTU, JendTV)
!
!=======================================================================
!                                                                      !
!  This routine computes the starting and ending horizontal indices    !
!  for each  C-staggered variable  from the sub-domain partition or    !
!  tile.                                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     IstrR      Starting tile index in the I-direction (RHO-points).  !
!     IstrT      Starting nest tile  in the I-direction (RHO-points).  !
!     IstrU      Starting tile index in the I-direction (U-points).    !
!     IendR      Ending   tile index in the I-direction (RHO_points).  !
!     IendT      Ending   nest tile  in the I-direction (RHO_points).  !
!                                                                      !
!     JstrR      Starting tile index in the J-direction (RHO-points).  !
!     JstrT      Starting nest tile  in the J-direction (RHO-points).  !
!     JstrV      Starting tile index in the J-direction (V-points).    !
!     JendR      Ending   tile index in the J-direction (RHO_points).  !
!     JendT      Ending   nest tile  in the J-direction (RHO-points).  !
!                                                                      !
!======================================================================!
!
      USE mod_param
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr
      integer, intent(out) :: IstrR, IstrT, IstrU, IendR, IendT
      integer, intent(out) :: JstrR, JstrT, JstrV, JendR, JendT
      integer, intent(out) :: IstrTU, JstrTV, IendTU, JendTV
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables.
!-----------------------------------------------------------------------
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
!         1     2                                        Lm     L
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
!
!
!
!                     *  *  *  *  *  *  *  *  *  *  *
!                               Ghost Points
!                     *  *  *  *  *  *  *  *  *  *  *
!
!                     *  *  p--v--p--v--p--v--p  *  *   Jend   IstrR=Istr
!                           |     |     |     |                IstrT=Istr
!     Interior        *  *  u  r  u  r  u  r  u  *  *          IstrU=Istr
!     Tile                  |     |     |     |                IendR=Iend
!                     *  *  p--v--p--v--p--v--p  *  *          IendT=Iend
!                           |     |     |     |                JstrR=Jstr
!                     *  *  u  r  u  r  u  r  u  *  *          JstrT=Jstr
!                           |     |     |     |                JstrV=Jstr
!                     *  *  p--v--p--v--p--v--p  *  *   Jstr   JendR=Jend
!                                                              JendT=Jend
!                     *  *  *  *  *  *  *  *  *  *  *
!
!                     *  *  *  *  *  *  *  *  *  *  *
!
!                          Istr              Iend
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
! JstrV v  p--v--p--v--p--v--p  *  *      *  *  p--v--p--v--p--v--p  v
!       :  +     |     |     |                  |     |     |     +  :
!       r  u  r  u  r  u  r  u  *  *      *  *  u  r  u  r  u  r  u  r
!       :  +     |     |     |                  |     |     |     +  :
!       v  p++v++p++v++p++v++p  *  * Jstr *  *  p++v++p++v++p++v++p  v
!       :                                                            :
!       r..u..r..u..r..u..r..u  *  *      *  *  u..r..u..r..u..r..u..r
!
!     IstrR    IstrU                                               IendR
!
!           Left/Bottom Tile                    Right/Bottom Tile
!
!
!  It also computes loop-bounds for U- and V-type variables which
!  belong to the interior of the computational domain. These are
!  labelled by suffixes U,V and they step one grid point inward from
!  the side of the subdomain adjacent to the physical boundary.
!  Conversely, for an internal subdomain which does not include a
!  segments of the physical boundary, all bounds with suffixes R,U,V
!  are set to the same values of corresponding non-suffixed bounds.
! 
!  Notice that the indices IstrT, IendT, JstrT, JendT are used during
!  nesting and include the overlap contact points between refined,
!  mosaic, or composed grids.  If not nesting, they are set to the
!  same values as IstrR, IendR, JstrR, JendR, respectively.
!
      IF (WESTERN_EDGE) THEN      ! Western Boundary Tile
#if defined EW_PERIODIC_REFINED
        IF (ng.eq.1) THEN
          IstrR=Istr
          IstrU=Istr
          IstrT=IstrR
          IstrTU=Istr-1
        ELSE
          IstrR=Istr-1
          IstrU=Istr+1
          IstrT=-NghostPoints
          IstrTU=IstrT
        END IF
#elif defined EW_PERIODIC
        IstrR=Istr
        IstrT=IstrR
        IstrU=Istr
#elif defined COMPOSED_GRID
        IstrR=Istr
        IstrT=-NghostPoints
        IstrU=Istr
        IstrTU=IstrT
#elif defined REFINED_GRID
        IstrR=Istr-1
        IstrU=Istr+1
        IF (ng.eq.1) THEN
          IstrT=IstrR
          IstrTU=IstrT
        ELSE
          IstrT=-NghostPoints
          IstrTU=IstrT
        END IF
#else
        IstrR=Istr-1
        IstrT=IstrR
        IstrU=Istr+1
        IstrTU=IstrT
#endif
      ELSE                        ! Interior Tile
        IstrR=Istr
        IstrT=IstrR
        IstrU=Istr
        IstrTU=IstrT-1
      END IF
      IF (EASTERN_EDGE) THEN      ! Eastern Boundary Tile
#if defined EW_PERIODIC_REFINED
        IF (ng.eq.1) THEN
          IendR=Iend
          IendT=IendR
        ELSE
          IendR=Iend+1
          IendT=Iend+NghostPoints
        END IF
#elif defined EW_PERIODIC
        IendR=Iend
        IendT=IendR
#elif defined COMPOSED_GRID
        IendR=Iend
        IendT=Iend+NghostPoints
#elif defined REFINED_GRID
        IendR=Iend+1
        IF (ng.eq.1) THEN
          IendT=IendR
          IendTU=IendT-1
        ELSE
          IendT=Iend+NghostPoints
          IendTU=IendT-1
        END IF
#else
        IendR=Iend+1
        IendT=IendR
#endif
      ELSE                        ! Interior Tile
        IendR=Iend
        IendT=IendR
        IendTU=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN     ! Southern Boundary Tile
#if defined NS_PERIODIC
        JstrR=Jstr
        JstrT=JstrR
        JstrV=Jstr
#elif defined COMPOSED_GRID
        JstrR=Jstr
        JstrT=-NghostPoints
        JstrV=Jstr
        JstrTV=JstrT
#elif defined REFINED_GRID
        JstrR=Jstr-1
        JstrV=Jstr+1
        IF (ng.eq.1) THEN
          JstrT=JstrR
          JstrTV=JstrT
        ELSE
          JstrT=-NghostPoints
          JstrTV=JstrT
        END IF
#else
        JstrR=Jstr-1
        JstrT=JstrR
        JstrV=Jstr+1
#endif
      ELSE                        ! Interior Tile
        JstrR=Jstr
        JstrT=JstrR
        JstrV=Jstr
        JstrTV=JstrT-1
      END IF
      IF (NORTHERN_EDGE) THEN     ! Northern Boundary Tile
#if defined NS_PERIODIC
        JendR=Jend
        JendT=JendR
#elif defined COMPOSED_GRID
        JendR=Jend
        JendT=Jend+NghostPoints
#elif defined REFINED_GRID
        JendR=Jend+1
        IF (ng.eq.1) THEN
          JendT=JendR
          JendTV=JendT-1
        ELSE
          JendT=Jend+NghostPoints
          JendTV=JendT-1
        END IF
#else
        JendR=Jend+1
        JendT=JendR
#endif
      ELSE                        ! Interior Tile
        JendR=Jend
        JendT=JendR
        JendTV=Jend
      END IF

      RETURN
      END SUBROUTINE var_bounds
