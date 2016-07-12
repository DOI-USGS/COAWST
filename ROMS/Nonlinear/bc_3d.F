#include "cppdefs.h"
      MODULE bc_3d_mod
#ifdef SOLVE3D
!
!svn $Id: bc_3d.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines apply close, gradient or periodic boundary           !
!  conditions to generic 3D fields.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                Nested grid number.                            !
!     tile              Domain partition.                              !
!     LBi               I-dimension Lower bound.                       !
!     UBi               I-dimension Upper bound.                       !
!     LBj               J-dimension Lower bound.                       !
!     UBj               J-dimension Upper bound.                       !
!     LBk               K-dimension Lower bound.                       !
!     UBk               K-dimension Upper bound.                       !
!     A                 3D field.                                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A                 Processed 3D field.                            !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     bc_r3d_tile       Boundary conditions for field at RHO-points    !
!     bc_u3d_tile       Boundary conditions for field at U-points      !
!     bc_v3d_tile       Boundary conditions for field at V-points      !
!     bc_w3d_tile       Boundary conditions for field at W-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bc_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                A(Iend+1,j,k)=A(Iend,j,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                A(Istr-1,j,k)=A(Istr,j,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                A(i,Jend+1,k)=A(i,Jend,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                A(i,Jstr-1,k)=A(i,Jstr,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO k=LBk,UBk
              A(Istr-1,Jstr-1,k)=0.5_r8*(A(Istr  ,Jstr-1,k)+            &
     &                                   A(Istr-1,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jstr-1,k)=0.5_r8*(A(Iend  ,Jstr-1,k)+            &
     &                                   A(Iend+1,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=LBk,UBk
              A(Istr-1,Jend+1,k)=0.5_r8*(A(Istr-1,Jend  ,k)+            &
     &                                   A(Istr  ,Jend+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+            &
     &                                   A(Iend  ,Jend+1,k))
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_r3d_tile

!
!***********************************************************************
      SUBROUTINE bc_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_u3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed or gradient.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (LBC(ieast,isBu3d,ng)%closed) THEN
            DO k=LBk,UBk
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%east(j)) THEN
                  A(Iend+1,j,k)=0.0_r8
                END IF
              END DO
            END DO
          ELSE
            DO k=LBk,UBk
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%east(j)) THEN
                  A(Iend+1,j,k)=A(Iend,j,k)
                END IF
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (LBC(iwest,isBu3d,ng)%closed) THEN
            DO k=LBk,UBk
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%west(j)) THEN
                  A(Istr,j,k)=0.0_r8
                END IF
              END DO
            END DO
          ELSE
            DO k=LBk,UBk
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%west(j)) THEN
                  A(Istr,j,k)=A(Istr+1,j,k)
                END IF
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (LBC(inorth,isBu3d,ng)%closed) THEN
            IF (EWperiodic(ng)) THEN
              Imin=IstrU
              Imax=Iend
            ELSE
              Imin=Istr
              Imax=IendR
            END IF
            DO k=LBk,UBk
              DO i=Imin,Imax
                IF (LBC_apply(ng)%north(i)) THEN
                  A(i,Jend+1,k)=gamma2(ng)*A(i,Jend,k)
# ifdef MASKING
                  A(i,Jend+1,k)=A(i,Jend+1,k)*GRID(ng)%umask(i,Jend+1)
# endif
                END IF
              END DO
            END DO
          ELSE
            DO k=LBk,UBk
              DO i=IstrU,Iend
                IF (LBC_apply(ng)%north(i)) THEN
                  A(i,Jend+1,k)=A(i,Jend,k)
                END IF
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (LBC(isouth,isBu3d,ng)%closed) THEN
            IF (EWperiodic(ng)) THEN
              Imin=IstrU
              Imax=Iend
            ELSE
              Imin=Istr
              Imax=IendR
            END IF
            DO k=LBk,UBk
              DO i=Imin,Imax
                IF (LBC_apply(ng)%south(i)) THEN
                  A(i,Jstr-1,k)=gamma2(ng)*A(i,Jstr,k)
# ifdef MASKING
                  A(i,Jstr-1,k)=A(i,Jstr-1,k)*GRID(ng)%umask(i,Jstr-1)
# endif
                END IF
              END DO
            END DO
          ELSE
            DO k=LBk,UBk
              DO i=IstrU,Iend
                IF (LBC_apply(ng)%south(i)) THEN
                  A(i,Jstr-1,k)=A(i,Jstr,k)
                END IF
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
             DO k=LBk,UBk
              A(Istr  ,Jstr-1,k)=0.5_r8*(A(Istr+1,Jstr-1,k)+            &
     &                                   A(Istr  ,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jstr-1,k)=0.5_r8*(A(Iend  ,Jstr-1,k)+            &
     &                                   A(Iend+1,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=LBk,UBk
              A(Istr  ,Jend+1,k)=0.5_r8*(A(Istr  ,Jend  ,k)+            &
     &                                   A(Istr+1,Jend+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+            &
     &                                   A(Iend  ,Jend+1,k))
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_u3d_tile

!
!***********************************************************************
      SUBROUTINE bc_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_v3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (LBC(ieast,isBv3d,ng)%closed) THEN
            IF (NSperiodic(ng)) THEN
              Jmin=JstrV
              Jmax=Jend
            ELSE
              Jmin=Jstr
              Jmax=JendR
            END IF
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                IF (LBC_apply(ng)%east(j)) THEN
                  A(Iend+1,j,k)=gamma2(ng)*A(Iend,j,k)
# ifdef MASKING
                  A(Iend+1,j,k)=A(Iend+1,j,k)*GRID(ng)%vmask(Iend+1,j)
# endif
                END IF
              END DO
            END DO
          ELSE
            DO k=LBk,UBk
              DO j=JstrV,Jend
                IF (LBC_apply(ng)%east(j)) THEN
                  A(Iend+1,j,k)=A(Iend,j,k)
                END IF
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (LBC(iwest,isBv3d,ng)%closed) THEN
            IF (NSperiodic(ng)) THEN
              Jmin=JstrV
              Jmax=Jend
            ELSE
              Jmin=Jstr
              Jmax=JendR
            END IF
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                IF (LBC_apply(ng)%west(j)) THEN
                  A(Istr-1,j,k)=gamma2(ng)*A(Istr,j,k)
# ifdef MASKING
                  A(Istr-1,j,k)=A(Istr-1,j,k)*GRID(ng)%vmask(Istr-1,j)
# endif
                END IF
              END DO
            END DO
          ELSE
            DO k=LBk,UBk
              DO j=JstrV,Jend
                IF (LBC_apply(ng)%west(j)) THEN
                  A(Istr-1,j,k)=A(Istr,j,k)
                END IF
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed or gradient.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (LBC(inorth,isBv3d,ng)%closed) THEN
            DO k=LBk,UBk
              DO i=Istr,Iend
                IF (LBC_apply(ng)%north(i)) THEN
                  A(i,Jend+1,k)=0.0_r8
                END IF
              END DO
            END DO
          ELSE
            DO k=LBk,UBk
              DO i=Istr,Iend
                IF (LBC_apply(ng)%north(i)) THEN
                  A(i,Jend+1,k)=A(i,Jend,k)
                END IF
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (LBC(isouth,isBv3d,ng)%closed) THEN
            DO k=LBk,UBk
              DO i=Istr,Iend
                IF (LBC_apply(ng)%south(i)) THEN
                  A(i,Jstr,k)=0.0_r8
                END IF
              END DO
            END DO
          ELSE
            DO k=LBk,UBk
              DO i=Istr,Iend
                IF (LBC_apply(ng)%south(i)) THEN
                  A(i,Jstr,k)=A(i,Jstr+1,k)
                END IF
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
            DO k=LBk,UBk
              A(Istr-1,Jstr  ,k)=0.5_r8*(A(Istr  ,Jstr  ,k)+            &
     &                                   A(Istr-1,Jstr+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            DO k=LBk,UBk
              A(Iend+1,Jstr  ,k)=0.5_r8*(A(Iend  ,Jstr  ,k)+            &
     &                                   A(Iend+1,Jstr+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=LBk,UBk
              A(Istr-1,Jend+1,k)=0.5_r8*(A(Istr-1,Jend  ,k)+            &
     &                                   A(Istr  ,Jend+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+            &
     &                                   A(Iend  ,Jend+1,k))
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_v3d_tile

!
!***********************************************************************
      SUBROUTINE bc_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                A(Iend+1,j,k)=A(Iend,j,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                A(Istr-1,j,k)=A(Istr,j,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                A(i,Jend+1,k)=A(i,Jend,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                A(i,Jstr-1,k)=A(i,Jstr,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO k=LBk,UBk
              A(Istr-1,Jstr-1,k)=0.5_r8*(A(Istr  ,Jstr-1,k)+            &
     &                                   A(Istr-1,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jstr-1,k)=0.5_r8*(A(Iend  ,Jstr-1,k)+            &
     &                                   A(Iend+1,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=LBk,UBk
              A(Istr-1,Jend+1,k)=0.5_r8*(A(Istr-1,Jend  ,k)+            &
     &                                   A(Istr  ,Jend+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+            &
     &                                   A(Iend  ,Jend+1,k))
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_w3d_tile

!
!***********************************************************************
      SUBROUTINE dabc_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                A(Iend+1,j,k)=A(Iend,j,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                A(Istr-1,j,k)=A(Istr,j,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                A(i,Jend+1,k)=A(i,Jend,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                A(i,Jstr-1,k)=A(i,Jstr,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO k=LBk,UBk
              A(Istr-1,Jstr-1,k)=0.5_r8*(A(Istr  ,Jstr-1,k)+            &
     &                                   A(Istr-1,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jstr-1,k)=0.5_r8*(A(Iend  ,Jstr-1,k)+            &
     &                                   A(Iend+1,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=LBk,UBk
              A(Istr-1,Jend+1,k)=0.5_r8*(A(Istr-1,Jend  ,k)+            &
     &                                   A(Istr  ,Jend+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+            &
     &                                   A(Iend  ,Jend+1,k))
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE dabc_r3d_tile

!
!***********************************************************************
      SUBROUTINE dabc_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_u3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                A(Iend+1,j,k)=A(Iend,j,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                A(Istr,j,k)=A(Istr+1,j,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                A(i,Jend+1,k)=A(i,Jend,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                A(i,Jstr-1,k)=A(i,Jstr,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
             DO k=LBk,UBk
              A(Istr  ,Jstr-1,k)=0.5_r8*(A(Istr+1,Jstr-1,k)+            &
     &                                   A(Istr  ,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jstr-1,k)=0.5_r8*(A(Iend  ,Jstr-1,k)+            &
     &                                   A(Iend+1,Jstr  ,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=LBk,UBk
              A(Istr  ,Jend+1,k)=0.5_r8*(A(Istr  ,Jend  ,k)+            &
     &                                   A(Istr+1,Jend+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+            &
     &                                   A(Iend  ,Jend+1,k))
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE dabc_u3d_tile

!
!***********************************************************************
      SUBROUTINE dabc_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_v3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j, k

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                A(Iend+1,j,k)=A(Iend,j,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO k=LBk,UBk
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                A(Istr-1,j,k)=A(Istr,j,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                A(i,Jend+1,k)=A(i,Jend,k)
              END IF
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO k=LBk,UBk
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                A(i,Jstr,k)=A(i,Jstr+1,k)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
            DO k=LBk,UBk
              A(Istr-1,Jstr  ,k)=0.5_r8*(A(Istr  ,Jstr  ,k)+            &
     &                                   A(Istr-1,Jstr+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            DO k=LBk,UBk
              A(Iend+1,Jstr  ,k)=0.5_r8*(A(Iend  ,Jstr  ,k)+            &
     &                                   A(Iend+1,Jstr+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=LBk,UBk
              A(Istr-1,Jend+1,k)=0.5_r8*(A(Istr-1,Jend  ,k)+            &
     &                                   A(Istr  ,Jend+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=LBk,UBk
              A(Iend+1,Jend+1,k)=0.5_r8*(A(Iend+1,Jend  ,k)+            &
     &                                   A(Iend  ,Jend+1,k))
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE dabc_v3d_tile
#endif
      END MODULE bc_3d_mod
