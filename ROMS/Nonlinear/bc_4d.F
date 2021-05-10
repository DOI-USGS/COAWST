#include "cppdefs.h"
      MODULE bc_4d_mod
#ifdef SOLVE3D
!
!svn $Id: bc_4d.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines apply close, gradient or periodic boundary           !
!  conditions to generic 4D fields.                                    !
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
!     LBl               L-dimension Lower bound.                       !
!     UBl               L-dimension Upper bound.                       !
!     A                 4D field.                                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A                 Processed 4D field.                            !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     bc_r4d_tile       Boundary conditions for field at RHO-points    !
!     bc_u4d_tile       Boundary conditions for field at U-points      !
!     bc_v4d_tile       Boundary conditions for field at V-points      !
!     bc_w4d_tile       Boundary conditions for field at W-points      !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bc_r4d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl,   &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_4d_mod, ONLY : exchange_r4d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,:,LBl:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBl:UBl)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, l

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO l=LBl,UBl
            DO k=LBk,UBk
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%east(j)) THEN
                  A(Iend+1,j,k,l)=A(Iend,j,k,l)
                END IF
              END DO
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO l=LBl,UBl
            DO k=LBk,UBk
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%west(j)) THEN
                  A(Istr-1,j,k,l)=A(Istr,j,k,l)
                END IF
              END DO
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
          DO l=LBl,UBl
            DO k=LBk,UBk
              DO i=Istr,Iend
                IF (LBC_apply(ng)%north(i)) THEN
                  A(i,Jend+1,k,l)=A(i,Jend,k,l)
                END IF
              END DO
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO l=LBl,UBl
            DO k=LBk,UBk
              DO i=Istr,Iend
                IF (LBC_apply(ng)%south(i)) THEN
                  A(i,Jstr-1,k,l)=A(i,Jstr,k,l)
                END IF
              END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Istr-1,Jstr-1,k,l)=0.5_r8*(A(Istr  ,Jstr-1,k,l)+      &
     &                                       A(Istr-1,Jstr  ,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Iend+1,Jstr-1,k,l)=0.5_r8*(A(Iend  ,Jstr-1,k,l)+      &
     &                                       A(Iend+1,Jstr  ,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Istr-1,Jend+1,k,l)=0.5_r8*(A(Istr-1,Jend  ,k,l)+      &
     &                                       A(Istr  ,Jend+1,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Iend+1,Jend+1,k,l)=0.5_r8*(A(Iend+1,Jend  ,k,l)+      &
     &                                       A(Iend  ,Jend+1,k,l))
              END DO
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
        CALL exchange_r4d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl, &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_r4d_tile

!
!***********************************************************************
      SUBROUTINE bc_u4d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl,   &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_4d_mod, ONLY : exchange_u4d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,:,LBl:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBl:UBl)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j, k, l

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed or gradient.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (LBC(ieast,isBu3d,ng)%closed) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jstr,Jend
                  IF (LBC_apply(ng)%east(j)) THEN
                    A(Iend+1,j,k,l)=0.0_r8
                  END IF
                END DO
              END DO
            END DO
          ELSE
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jstr,Jend
                  IF (LBC_apply(ng)%east(j)) THEN
                    A(Iend+1,j,k,l)=A(Iend,j,k,l)
                  END IF
                END DO
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (LBC(iwest,isBu3d,ng)%closed) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jstr,Jend
                  IF (LBC_apply(ng)%west(j)) THEN
                    A(Istr,j,k,l)=0.0_r8
                  END IF
                END DO
              END DO
            END DO
          ELSE
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jstr,Jend
                  IF (LBC_apply(ng)%west(j)) THEN
                    A(Istr,j,k,l)=A(Istr+1,j,k,l)
                  END IF
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  IF (LBC_apply(ng)%north(i)) THEN
                    A(i,Jend+1,k,l)=gamma2(ng)*A(i,Jend,k,l)
# ifdef MASKING
                    A(i,Jend+1,k,l)=A(i,Jend+1,k,l)*                    &
     &                              GRID(ng)%umask(i,Jend+1)
# endif
                  END IF
                END DO
              END DO
            END DO
          ELSE
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=IstrU,Iend
                  IF (LBC_apply(ng)%north(i)) THEN
                    A(i,Jend+1,k,l)=A(i,Jend,k,l)
                  END IF
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Imin,Imax
                  IF (LBC_apply(ng)%south(i)) THEN
                    A(i,Jstr-1,k,l)=gamma2(ng)*A(i,Jstr,k,l)
# ifdef MASKING
                    A(i,Jstr-1,k,l)=A(i,Jstr-1,k,l)*                    &
     &                              GRID(ng)%umask(i,Jstr-1)
# endif
                  END IF
                END DO
              END DO
            END DO
          ELSE
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=IstrU,Iend
                  IF (LBC_apply(ng)%south(i)) THEN
                    A(i,Jstr-1,k,l)=A(i,Jstr,k,l)
                  END IF
                END DO
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
             DO l=LBl,UBl
               DO k=LBk,UBk
                 A(Istr  ,Jstr-1,k,l)=0.5_r8*(A(Istr+1,Jstr-1,k,l)+     &
     &                                        A(Istr  ,Jstr  ,k,l))
               END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Iend+1,Jstr-1,k,l)=0.5_r8*(A(Iend  ,Jstr-1,k,l)+      &
     &                                       A(Iend+1,Jstr  ,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Istr  ,Jend+1,k,l)=0.5_r8*(A(Istr  ,Jend  ,k,l)+      &
     &                                       A(Istr+1,Jend+1,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Iend+1,Jend+1,k,l)=0.5_r8*(A(Iend+1,Jend  ,k,l)+      &
     &                                       A(Iend  ,Jend+1,k,l))
              END DO
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
        CALL exchange_u4d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl, &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_u4d_tile

!
!***********************************************************************
      SUBROUTINE bc_v4d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl,   &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_4d_mod, ONLY : exchange_v4d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,:,LBl:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBl:UBl)
# endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j, k, l

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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  IF (LBC_apply(ng)%east(j)) THEN
                    A(Iend+1,j,k,l)=gamma2(ng)*A(Iend,j,k,l)
# ifdef MASKING
                    A(Iend+1,j,k,l)=A(Iend+1,j,k,l)*                    &
     &                              GRID(ng)%vmask(Iend+1,j)
# endif
                  END IF
                END DO
              END DO
            END DO
          ELSE
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=JstrV,Jend
                  IF (LBC_apply(ng)%east(j)) THEN
                    A(Iend+1,j,k,l)=A(Iend,j,k,l)
                  END IF
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  IF (LBC_apply(ng)%west(j)) THEN
                    A(Istr-1,j,k,l)=gamma2(ng)*A(Istr,j,k,l)
# ifdef MASKING
                    A(Istr-1,j,k,l)=A(Istr-1,j,k,l)*                    &
     &                              GRID(ng)%vmask(Istr-1,j)
# endif
                  END IF
                END DO
              END DO
            END DO
          ELSE
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO j=JstrV,Jend
                  IF (LBC_apply(ng)%west(j)) THEN
                    A(Istr-1,j,k,l)=A(Istr,j,k,l)
                  END IF
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Istr,Iend
                  IF (LBC_apply(ng)%north(i)) THEN
                    A(i,Jend+1,k,l)=0.0_r8
                  END IF
                END DO
              END DO
            END DO
          ELSE
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Istr,Iend
                  IF (LBC_apply(ng)%north(i)) THEN
                    A(i,Jend+1,k,l)=A(i,Jend,k,l)
                  END IF
                END DO
              END DO
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (LBC(isouth,isBv3d,ng)%closed) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Istr,Iend
                  IF (LBC_apply(ng)%south(i)) THEN
                    A(i,Jstr,k,l)=0.0_r8
                  END IF
                END DO
              END DO
            END DO
          ELSE
            DO l=LBl,UBl
              DO k=LBk,UBk
                DO i=Istr,Iend
                  IF (LBC_apply(ng)%south(i)) THEN
                    A(i,Jstr,k,l)=A(i,Jstr+1,k,l)
                  END IF
                END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Istr-1,Jstr  ,k,l)=0.5_r8*(A(Istr  ,Jstr  ,k,l)+      &
     &                                       A(Istr-1,Jstr+1,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Iend+1,Jstr  ,k,l)=0.5_r8*(A(Iend  ,Jstr  ,k,l)+      &
     &                                       A(Iend+1,Jstr+1,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Istr-1,Jend+1,k,l)=0.5_r8*(A(Istr-1,Jend  ,k,l)+      &
     &                                       A(Istr  ,Jend+1,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Iend+1,Jend+1,k,l)=0.5_r8*(A(Iend+1,Jend  ,k,l)+      &
     &                                       A(Iend  ,Jend+1,k,l))
              END DO
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
        CALL exchange_v4d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl, &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_v4d_tile

!
!***********************************************************************
      SUBROUTINE bc_w4d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl,   &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_4d_mod, ONLY : exchange_w4d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:,LBl:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj,LBk:UBk,LBl:UBl)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, l

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO l=LBl,UBl
            DO k=LBk,UBk
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%east(j)) THEN
                  A(Iend+1,j,k,l)=A(Iend,j,k,l)
                END IF
              END DO
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO l=LBl,UBl
            DO k=LBk,UBk
              DO j=Jstr,Jend
                IF (LBC_apply(ng)%west(j)) THEN
                  A(Istr-1,j,k,l)=A(Istr,j,k,l)
                END IF
              END DO
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
          DO l=LBl,UBl
            DO k=LBk,UBk
              DO i=Istr,Iend
                IF (LBC_apply(ng)%north(i)) THEN
                  A(i,Jend+1,k,l)=A(i,Jend,k,l)
                END IF
              END DO
            END DO
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO l=LBl,UBl
            DO k=LBk,UBk
              DO i=Istr,Iend
                IF (LBC_apply(ng)%south(i)) THEN
                  A(i,Jstr-1,k,l)=A(i,Jstr,k,l)
                END IF
              END DO
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
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Istr-1,Jstr-1,k,l)=0.5_r8*(A(Istr  ,Jstr-1,k,l)+      &
     &                                       A(Istr-1,Jstr  ,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Iend+1,Jstr-1,k,l)=0.5_r8*(A(Iend  ,Jstr-1,k,l)+      &
     &                                       A(Iend+1,Jstr  ,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Istr-1,Jend+1,k,l)=0.5_r8*(A(Istr-1,Jend  ,k,l)+      &
     &                                       A(Istr  ,Jend+1,k,l))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO l=LBl,UBl
              DO k=LBk,UBk
                A(Iend+1,Jend+1,k,l)=0.5_r8*(A(Iend+1,Jend  ,k,l)+      &
     &                                       A(Iend  ,Jend+1,k,l))
              END DO
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
        CALL exchange_w4d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk, LBl, UBl, &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_w4d_tile
#endif
      END MODULE bc_4d_mod
