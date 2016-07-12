#include "cppdefs.h"
      MODULE ad_bc_bry3d_mod
#ifdef ADJOINT
!
!svn $Id: ad_bc_bry3d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package applies adjoint gradient conditions for generic 3D     !
!  boundary fields.                                                    !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    ad_bc_r3d_bry_tile    Boundary conditions for field at RHO-points !
!    ad_bc_u3d_bry_tile    Boundary conditions for field at U-points   !
!    ad_bc_v3d_bry_tile    Boundary conditions for field at V-points   !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_bc_r3d_bry_tile (ng, tile, boundary,                &
     &                               LBij, UBij, LBk, UBk,              &
     &                               ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij, LBk, UBk

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBij:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.inorth) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Iend+1,k)=tl_A(Iend,k)
!>
            ad_A(Iend,k)=ad_A(Iend,k)+ad_A(Iend+1,k)
            ad_A(Iend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Istr-1,k)=tl_A(Istr,k)
!>
            ad_A(Istr,k)=ad_A(Istr,k)+ad_A(Istr-1,k)
            ad_A(Istr-1,k)=0.0_r8
          END DO
        END IF
      END IF

      IF (boundary.eq.isouth) THEN
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Iend+1,k)=tl_A(Iend,k)
!>
            ad_A(Iend,k)=ad_A(Iend,k)+ad_A(Iend+1,k)
            ad_A(Iend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Istr-1,k)=tl_A(Istr,k)
!>
            ad_A(Istr,k)=ad_A(Istr,k)+ad_A(Istr-1,k)
            ad_A(Istr-1,k)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.ieast) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jend+1,k)=tl_A(Jend,k)
!>
            ad_A(Jend,k)=ad_A(Jend,k)+ad_A(Jend+1,k)
            ad_A(Jend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jstr-1,k)=tl_A(Jstr,k)
!>
            ad_A(Jstr,k)=ad_A(Jstr,k)+ad_A(Jstr-1,k)
            ad_A(Jstr-1,k)=0.0_r8
          END DO
        END IF
      END IF

      IF (boundary.eq.iwest) THEN
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jend+1,k)=tl_A(Jend,k)
!>
            ad_A(Jend,k)=ad_A(Jend,k)+ad_A(Jend+1,k)
            ad_A(Jend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jstr-1,k)=tl_A(Jstr,k)
!>
            ad_A(Jstr,k)=ad_A(Jstr,k)+ad_A(Jstr-1,k)
            ad_A(Jstr-1,k)=0.0_r8
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_bc_r3d_bry_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_u3d_bry_tile (ng, tile, boundary,                &
     &                               LBij, UBij, LBk, UBk,              &
     &                               ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij, LBk, UBk

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBij:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.inorth) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Iend+1,k)=tl_A(Iend,k)
!>
            ad_A(Iend,k)=ad_A(Iend,k)+ad_A(Iend+1,k)
            ad_A(Iend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(IstrU-1,k)=tl_A(IstrU,k)
!>
            ad_A(IstrU,k)=ad_A(IstrU,k)+ad_A(IstrU-1,k)
            ad_A(IstrU-1,k)=0.0_r8
          END DO
        END IF
      END IF

      IF (boundary.eq.isouth) THEN
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Iend+1,k)=tl_A(Iend,k)
!>
            ad_A(Iend,k)=ad_A(Iend,k)+ad_A(Iend+1,k)
            ad_A(Iend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(IstrU-1,k)=tl_A(IstrU,k)
!>
            ad_A(IstrU,k)=ad_A(IstrU,k)+ad_A(IstrU-1,k)
            ad_A(IstrU-1,k)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.ieast) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jend+1,k)=tl_A(Jend,k)
!>
            ad_A(Jend,k)=ad_A(Jend,k)+ad_A(Jend+1,k)
            ad_A(Jend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jstr-1,k)=tl_A(Jstr,k)
!>
            ad_A(Jstr,k)=ad_A(Jstr,k)+ad_A(Jstr-1,k)
            ad_A(Jstr-1,k)=0.0_r8
          END DO
        END IF
      END IF

      IF (boundary.eq.iwest) THEN
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jend+1,k)=tl_A(Jend,k)
!>
            ad_A(Jend,k)=ad_A(Jend,k)+ad_A(Jend+1,k)
            ad_A(Jend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jstr-1,k)=tl_A(Jstr,k)
!>
            ad_A(Jstr,k)=ad_A(Jstr,k)+ad_A(Jstr-1,k)
            ad_A(Jstr-1,k)=0.0_r8
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_bc_u3d_bry_tile

!
!***********************************************************************
      SUBROUTINE ad_bc_v3d_bry_tile (ng, tile, boundary,                &
     &                               LBij, UBij, LBk, UBk,              &
     &                               ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, boundary
      integer, intent(in) :: LBij, UBij, LBk, UBk

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBij:,LBk:)
# else
      real(r8), intent(inout) :: ad_A(LBij:UBij,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint Southern and Northern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.inorth) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Iend+1,k)=tl_A(Iend,k)
!>
            ad_A(Iend,k)=ad_A(Iend,k)+ad_A(Iend+1,k)
            ad_A(Iend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Istr-1,k)=tl_A(Istr,k)
!>
            ad_A(Istr,k)=ad_A(Istr,k)+ad_A(Istr-1,k)
            ad_A(Istr-1,k)=0.0_r8
          END DO
        END IF
      END IF

      IF (boundary.eq.isouth) THEN
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Iend+1,k)=tl_A(Iend,k)
!>
            ad_A(Iend,k)=ad_A(Iend,k)+ad_A(Iend+1,k)
            ad_A(Iend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Istr-1,k)=tl_A(Istr,k)
!>
            ad_A(Istr,k)=ad_A(Istr,k)+ad_A(Istr-1,k)
            ad_A(Istr-1,k)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint Western and Eastern edges: gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (boundary.eq.ieast) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jend+1,k)=tl_A(Jend,k)
!>
            ad_A(Jend,k)=ad_A(Jend,k)+ad_A(Jend+1,k)
            ad_A(Jend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(JstrV-1,k)=tl_A(JstrV,k)
!>
            ad_A(JstrV,k)=ad_A(JstrV,k)+ad_A(JstrV-1,k)
            ad_A(JstrV-1,k)=0.0_r8
          END DO
        END IF
      END IF

      IF (boundary.eq.iwest) THEN
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(Jend+1,k)=tl_A(Jend,k)
!>
            ad_A(Jend,k)=ad_A(Jend,k)+ad_A(Jend+1,k)
            ad_A(Jend+1,k)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=LBk,UBk
!>          tl_A(JstrV-1,k)=tl_A(JstrV,k)
!>
            ad_A(JstrV,k)=ad_A(JstrV,k)+ad_A(JstrV-1,k)
            ad_A(JstrV-1,k)=0.0_r8
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_bc_v3d_bry_tile
#endif
      END MODULE ad_bc_bry3d_mod
