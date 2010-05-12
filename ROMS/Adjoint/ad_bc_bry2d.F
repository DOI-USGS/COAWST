#include "cppdefs.h"
      MODULE ad_bc_bry2d_mod
#ifdef ADJOINT
!
!svn $Id: ad_bc_bry2d.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package applies adjoint gradient conditions for generic 2D     !
!  boundary fields.                                                    !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    ad_bc_r2d_bry_tile    Boundary conditions for field at RHO-points !
!    ad_bc_u2d_bry_tile    Boundary conditions for field at U-points   !
!    ad_bc_v2d_bry_tile    Boundary conditions for field at V-points   !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_bc_r2d_bry_tile (ng, tile, boundary,                &
     &                               LBij, UBij,                        &
     &                               ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBij:)
# else
      real(r8), intent(inout) :: ad_A(LBij:UBij)
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.inorth) THEN
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>        tl_A(Iend+1)=tl_A(Iend)
!>
          ad_A(Iend)=ad_A(Iend)+ad_A(Iend+1)
          ad_A(Iend+1)=0.0_r8
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>        tl_A(Istr-1)=tl_A(Istr)
!>
          ad_A(Istr)=ad_A(Istr)+ad_A(Istr-1)
          ad_A(Istr-1)=0.0_r8
        END IF
      END IF

      IF (boundary.eq.isouth) THEN
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>        tl_A(Iend+1)=tl_A(Iend)
!>
          ad_A(Iend)=ad_A(Iend)+ad_A(Iend+1)
          ad_A(Iend+1)=0.0_r8
        END IF
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>        tl_A(Istr-1)=tl_A(Istr)
!>
          ad_A(Istr)=ad_A(Istr)+ad_A(Istr-1)
          ad_A(Istr-1)=0.0_r8
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.ieast) THEN
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Jend+1)=tl_A(Jend)
!>
          ad_A(Jend)=ad_A(Jend)+ad_A(Jend+1)
          ad_A(Jend+1)=0.0_r8
        END IF
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(Jstr-1)=tl_A(Jstr)
!>
          ad_A(Jstr)=ad_A(Jstr)+ad_A(Jstr-1)
          ad_A(Jstr-1)=0.0_r8
        END IF
      END IF

      IF (boundary.eq.iwest) THEN
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Jend+1)=tl_A(Jend)
!>
          ad_A(Jend)=ad_A(Jend)+ad_A(Jend+1)
          ad_A(Jend+1)=0.0_r8
        END IF
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(Jstr-1)=tl_A(Jstr)
!>
          ad_A(Jstr)=ad_A(Jstr)+ad_A(Jstr-1)
          ad_A(Jstr-1)=0.0_r8
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_bc_r2d_bry_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_u2d_bry_tile (ng, tile, boundary,                &
     &                               LBij, UBij,                        &
     &                               ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBij:)
# else
      real(r8), intent(inout) :: ad_A(LBij:UBij)
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.inorth) THEN
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>        tl_A(Iend+1)=tl_A(Iend)
!>
          ad_A(Iend)=ad_A(Iend)+ad_A(Iend+1)
          ad_A(Iend+1)=0.0_r8
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>        tl_A(IstrU-1)=tl_A(IstrU)
!>
          ad_A(IstrU)=ad_A(IstrU)+ad_A(IstrU-1)
          ad_A(IstrU-1)=0.0_r8
        END IF
      END IF

      IF (boundary.eq.isouth) THEN
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>        tl_A(Iend+1)=tl_A(Iend)
!>
          ad_A(Iend)=ad_A(Iend)+ad_A(Iend+1)
          ad_A(Iend+1)=0.0_r8
        END IF
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>        tl_A(IstrU-1)=tl_A(IstrU)
!>
          ad_A(IstrU)=ad_A(IstrU)+ad_A(IstrU-1)
          ad_A(IstrU-1)=0.0_r8
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.ieast) THEN
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Jend+1)=tl_A(Jend)
!>
          ad_A(Jend)=ad_A(Jend)+ad_A(Jend+1)
          ad_A(Jend+1)=0.0_r8
        END IF
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(Jstr-1)=tl_A(Jstr)
!>
          ad_A(Jstr)=ad_A(Jstr)+ad_A(Jstr-1)
          ad_A(Jstr-1)=0.0_r8
        END IF
      END IF

      IF (boundary.eq.iwest) THEN
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Jend+1)=tl_A(Jend)
!>
          ad_A(Jend)=ad_A(Jend)+ad_A(Jend+1)
          ad_A(Jend+1)=0.0_r8
        END IF
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(Jstr-1)=tl_A(Jstr)
!>
          ad_A(Jstr)=ad_A(Jstr)+ad_A(Jstr-1)
          ad_A(Jstr-1)=0.0_r8
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_bc_u2d_bry_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_v2d_bry_tile (ng, tile, boundary,                &
     &                               LBij, UBij,                        &
     &                               ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBij:)
# else
      real(r8), intent(inout) :: ad_A(LBij:UBij)
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.inorth) THEN
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>        tl_A(Iend+1)=tl_A(Iend)
!>
          ad_A(Iend)=ad_A(Iend)+ad_A(Iend+1)
          ad_A(Iend+1)=0.0_r8
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>        tl_A(Istr-1)=tl_A(Istr)
!>
          ad_A(Istr)=ad_A(Istr)+ad_A(Istr-1)
          ad_A(Istr-1)=0.0_r8
        END IF
      END IF

      IF (boundary.eq.isouth) THEN
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>        tl_A(Iend+1)=tl_A(Iend)
!>
          ad_A(Iend)=ad_A(Iend)+ad_A(Iend+1)
          ad_A(Iend+1)=0.0_r8
        END IF
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>        tl_A(Istr-1)=tl_A(Istr)
!>
          ad_A(Istr)=ad_A(Istr)+ad_A(Istr-1)
          ad_A(Istr-1)=0.0_r8
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.ieast) THEN
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Jend+1)=tl_A(Jend)
!>
          ad_A(Jend)=ad_A(Jend)+ad_A(Jend+1)
          ad_A(Jend+1)=0.0_r8
        END IF
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(JstrV-1)=tl_A(JstrV)
!>
          ad_A(JstrV)=ad_A(JstrV)+ad_A(JstrV-1)
          ad_A(JstrV-1)=0.0_r8
        END IF
      END IF

      IF (boundary.eq.iwest) THEN
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
!>        tl_A(Jend+1)=tl_A(Jend)
!>
          ad_A(Jend)=ad_A(Jend)+ad_A(Jend+1)
          ad_A(Jend+1)=0.0_r8
        END IF
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
!>        tl_A(JstrV-1)=tl_A(JstrV)
!>
          ad_A(JstrV)=ad_A(JstrV)+ad_A(JstrV-1)
          ad_A(JstrV-1)=0.0_r8
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_bc_v2d_bry_tile
#endif
      END MODULE ad_bc_bry2d_mod
