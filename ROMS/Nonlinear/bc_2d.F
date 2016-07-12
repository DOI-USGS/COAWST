#include "cppdefs.h"
      MODULE bc_2d_mod
!
!svn $Id: bc_2d.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines apply close, gradient or periodic boundary           !
!  conditions to generic 2D fields.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                Nested grid number.                            !
!     tile              Domain partition.                              !
!     LBi               I-dimension Lower bound.                       !
!     UBi               I-dimension Upper bound.                       !
!     LBj               J-dimension Lower bound.                       !
!     UBj               J-dimension Upper bound.                       !
!     A                 2D field.                                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A                 Processed 2D field.                            !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     bc_r2d_tile       Boundary conditions for field at RHO-points    !
!     bc_u2d_tile       Boundary conditions for field at U-points      !
!     bc_v2d_tile       Boundary conditions for field at V-points      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE bc_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              A(Iend+1,j)=A(Iend,j)
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              A(Istr-1,j)=A(Istr,j)
            END IF
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
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              A(i,Jend+1)=A(i,Jend)
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              A(i,Jstr-1)=A(i,Jstr)
            END IF
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
            A(Istr-1,Jstr-1)=0.5_r8*(A(Istr  ,Jstr-1)+                  &
     &                               A(Istr-1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                  &
     &                               A(Iend+1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                  &
     &                               A(Istr  ,Jend+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                  &
     &                               A(Iend  ,Jend+1))
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_r2d_tile

!
!***********************************************************************
      SUBROUTINE bc_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed or gradient
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (LBC(ieast,isBu2d,ng)%closed) THEN
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                A(Iend+1,j)=0.0_r8
              END IF
            END DO
          ELSE
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                A(Iend+1,j)=A(Iend,j)
              END IF
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (LBC(iwest,isBu2d,ng)%closed) THEN
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                A(Istr,j)=0.0_r8
              END IF
            END DO
          ELSE
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                A(Istr,j)=A(Istr+1,j)
              END IF
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
          IF (LBC(inorth,isBu2d,ng)%closed) THEN
            IF (EWperiodic(ng)) THEN
              Imin=IstrU
              Imax=Iend
            ELSE
              Imin=Istr
              Imax=IendR
            END IF
            DO i=Imin,Imax
              IF (LBC_apply(ng)%north(i)) THEN
                A(i,Jend+1)=gamma2(ng)*A(i,Jend)
#ifdef MASKING
                A(i,Jend+1)=A(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#endif
              END IF
            END DO
          ELSE
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                A(i,Jend+1)=A(i,Jend)
              END IF
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (LBC(isouth,isBu2d,ng)%closed) THEN
            IF (EWperiodic(ng)) THEN
              Imin=IstrU
              Imax=Iend
            ELSE
              Imin=Istr
              Imax=IendR
            END IF
            DO i=Imin,Imax
              IF (LBC_apply(ng)%south(i)) THEN
                A(i,Jstr-1)=gamma2(ng)*A(i,Jstr)
#ifdef MASKING
                A(i,Jstr-1)=A(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#endif
              END IF
            END DO
          ELSE
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                A(i,Jstr-1)=A(i,Jstr)
              END IF
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
            A(Istr  ,Jstr-1)=0.5_r8*(A(Istr+1,Jstr-1)+                  &
     &                               A(Istr  ,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                  &
     &                               A(Iend+1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            A(Istr  ,Jend+1)=0.5_r8*(A(Istr  ,Jend  )+                  &
     &                               A(Istr+1,Jend+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                  &
     &                               A(Iend  ,Jend+1))
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_u2d_tile

!
!***********************************************************************
      SUBROUTINE bc_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (LBC(ieast,isBv2d,ng)%closed) THEN
            IF (NSperiodic(ng)) THEN
              Jmin=JstrV
              Jmax=Jend
            ELSE
              Jmin=Jstr
              Jmax=JendR
            END IF
            DO j=Jmin,Jmax
              IF (LBC_apply(ng)%east(j)) THEN
                A(Iend+1,j)=gamma2(ng)*A(Iend,j)
#ifdef MASKING
                A(Iend+1,j)=A(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#endif
              END IF
            END DO
          ELSE
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                A(Iend+1,j)=A(Iend,j)
              END IF
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (LBC(iwest,isBv2d,ng)%closed) THEN
            IF (NSperiodic(ng)) THEN
              Jmin=JstrV
              Jmax=Jend
            ELSE
              Jmin=Jstr
              Jmax=JendR
            END IF
            DO j=Jmin,Jmax
              IF (LBC_apply(ng)%west(j)) THEN
                A(Istr-1,j)=gamma2(ng)*A(Istr,j)
#ifdef MASKING
                A(Istr-1,j)=A(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#endif
              END IF
            END DO
          ELSE
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                A(Istr-1,j)=A(Istr,j)
              END IF
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed or Gradient.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (LBC(inorth,isBv2d,ng)%closed) THEN
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                A(i,Jend+1)=0.0_r8
              END IF
            END DO
          ELSE
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                A(i,Jend+1)=A(i,Jend)
              END IF
            END DO
          END IF
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (LBC(isouth,isBv2d,ng)%closed) THEN
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                A(i,Jstr)=0.0_r8
              END IF
            END DO
          ELSE
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                A(i,Jstr)=A(i,Jstr+1)
              END IF
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
            A(Istr-1,Jstr  )=0.5_r8*(A(Istr  ,Jstr  )+                  &
     &                               A(Istr-1,Jstr+1))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            A(Iend+1,Jstr  )=0.5_r8*(A(Iend  ,Jstr  )+                  &
     &                               A(Iend+1,Jstr+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                  &
     &                               A(Istr  ,Jend+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                  &
     &                               A(Iend  ,Jend+1))
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE bc_v2d_tile

!
!***********************************************************************
      SUBROUTINE dabc_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              A(Iend+1,j)=A(Iend,j)
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              A(Istr-1,j)=A(Istr,j)
            END IF
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
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              A(i,Jend+1)=A(i,Jend)
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              A(i,Jstr-1)=A(i,Jstr)
            END IF
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
            A(Istr-1,Jstr-1)=0.5_r8*(A(Istr  ,Jstr-1)+                  &
     &                               A(Istr-1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                  &
     &                               A(Iend+1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                  &
     &                               A(Istr  ,Jend+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                  &
     &                               A(Iend  ,Jend+1))
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE dabc_r2d_tile

!
!***********************************************************************
      SUBROUTINE dabc_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
           IF (LBC_apply(ng)%east(j)) THEN
              A(Iend+1,j)=A(Iend,j)
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              A(Istr,j)=A(Istr+1,j)
            END IF
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
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              A(i,Jend+1)=A(i,Jend)
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              A(i,Jstr-1)=A(i,Jstr)
            END IF
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
            A(Istr  ,Jstr-1)=0.5_r8*(A(Istr+1,Jstr-1)+                  &
     &                               A(Istr  ,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                  &
     &                               A(Iend+1,Jstr  ))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            A(Istr  ,Jend+1)=0.5_r8*(A(Istr  ,Jend  )+                  &
     &                               A(Istr+1,Jend+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                  &
     &                               A(Iend  ,Jend+1))
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE dabc_u2d_tile

!
!***********************************************************************
      SUBROUTINE dabc_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              A(Iend+1,j)=A(Iend,j)
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              A(Istr-1,j)=A(Istr,j)
            END IF
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
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              A(i,Jend+1)=A(i,Jend)
            END IF
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              A(i,Jstr)=A(i,Jstr+1)
            END IF
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
            A(Istr-1,Jstr  )=0.5_r8*(A(Istr  ,Jstr  )+                  &
     &                               A(Istr-1,Jstr+1))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            A(Iend+1,Jstr  )=0.5_r8*(A(Iend  ,Jstr  )+                  &
     &                               A(Iend+1,Jstr+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                  &
     &                               A(Istr  ,Jend+1))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                  &
     &                               A(Iend  ,Jend+1))
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE dabc_v2d_tile

      END MODULE bc_2d_mod
