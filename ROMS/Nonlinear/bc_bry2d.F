#include "cppdefs.h"
      MODULE bc_bry2d_mod
!
!svn $Id: bc_bry2d.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF

      IF (boundary.eq.ieast) THEN
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      IF (boundary.eq.inorth) THEN
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
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
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF

      IF (boundary.eq.ieast) THEN
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(Jstr-1)=A(Jstr)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(IstrU-1)=A(IstrU)
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      IF (boundary.eq.inorth) THEN
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(IstrU-1)=A(IstrU)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
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
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(JstrV-1)=A(JstrV)
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF

      IF (boundary.eq.ieast) THEN
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(JstrV-1)=A(JstrV)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(Jend+1)=A(Jend)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.isouth) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      IF (boundary.eq.inorth) THEN
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Istr-1)=A(Istr)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(Iend+1)=A(Iend)
        END IF
      END IF

      RETURN
      END SUBROUTINE bc_v2d_bry_tile

      END MODULE bc_bry2d_mod
