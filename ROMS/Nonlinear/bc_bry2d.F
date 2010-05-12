#include "cppdefs.h"
      MODULE bc_bry2d_mod
!
!svn $Id: bc_bry2d.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package applies gradient conditions for generic 2D boundary    !
!  fields.                                                             !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    bc_r2d_bry_tile    Boundary conditions for field at RHO-points    !
!    bc_u2d_bry_tile    Boundary conditions for field at U-points      !
!    bc_v2d_bry_tile    Boundary conditions for field at V-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bc_r2d_bry_tile (ng, tile, boundary,                   &
     &                            LBij, UBij,                           &
     &                            A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBij:)
#else
      real(r8), intent(inout) :: A(LBij:UBij)
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.iwest) THEN
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF

      IF (boundary.eq.ieast) THEN
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      IF (boundary.eq.inorth) THEN
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      RETURN
      END SUBROUTINE bc_r2d_bry_tile

!
!***********************************************************************
      SUBROUTINE bc_u2d_bry_tile (ng, tile, boundary,                   &
     &                            LBij, UBij,                           &
     &                            A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBij:)
#else
      real(r8), intent(inout) :: A(LBij:UBij)
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.iwest) THEN
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF

      IF (boundary.eq.ieast) THEN
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          A(IstrU-1)=A(IstrU)
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      IF (boundary.eq.inorth) THEN
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          A(IstrU-1)=A(IstrU)
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      RETURN
      END SUBROUTINE bc_u2d_bry_tile

!
!***********************************************************************
      SUBROUTINE bc_v2d_bry_tile (ng, tile, boundary,                   &
     &                            LBij, UBij,                           &
     &                            A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBij:)
#else
      real(r8), intent(inout) :: A(LBij:UBij)
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.iwest) THEN
        IF ((WESTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          A(JstrV-1)=A(JstrV)
        END IF
        IF ((WESTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF

      IF (boundary.eq.ieast) THEN
        IF ((EASTERN_EDGE).and.(SOUTHERN_EDGE)) THEN
          A(JstrV-1)=A(JstrV)
        END IF
        IF ((EASTERN_EDGE).and.(NORTHERN_EDGE)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      IF (boundary.eq.inorth) THEN
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      RETURN
      END SUBROUTINE bc_v2d_bry_tile

      END MODULE bc_bry2d_mod
