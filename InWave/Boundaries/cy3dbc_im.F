#include "cppdefs.h"
      MODULE cy3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: cy3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!                                                                      !
!  This subroutine sets lateral boundary conditions for 3D             !
!  etai component of the wave group velocity                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: cy3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE cy3dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_inwave_vars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL cy3dbc_tile (ng, tile,                                       &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 WAVEP(ng) % cy)
      RETURN
      END SUBROUTINE cy3dbc
!
!***********************************************************************
      SUBROUTINE cy3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       cy)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_grid
      USE mod_boundary
      USE mod_scalars      
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: cy(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: cy(LBi:UBi,LBj:UBj,ND)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, d

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, gradient boundary condition.
!
        IF (LBC(isouth,isCY3d,ng)%gradient) THEN
          DO d=1,ND
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                cy(i,Jstr,d)=cy(i,Jstr+1,d)
#   ifdef MASKING
                cy(i,Jstr,d)=cy(i,Jstr,d)*                              &
     &                       GRID(ng)%vmask(i,Jstr)
#   endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isCY3d,ng)%closed) THEN
          DO d=1,ND
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                cy(i,Jstr,d)=0.0_r8
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, gradient boundary condition.
!
        IF (LBC(inorth,isCY3d,ng)%gradient) THEN
          DO d=1,ND
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                cy(i,Jend+1,d)=cy(i,Jend,d)
#   ifdef MASKING
                cy(i,Jend+1,d)=cy(i,Jend+1,d)*                          &
     &                         GRID(ng)%vmask(i,Jend+1)
#   endif
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isCY3d,ng)%closed) THEN
          DO d=1,ND
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                cy(i,Jend+1,d)=0.0_r8
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, gradient boundary condition.
!
        IF (LBC(iwest,isCY3d,ng)%gradient) THEN
          DO d=1,ND
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                cy(Istr-1,j,d)=cy(Istr,j,d)
#   ifdef MASKING
                cy(Istr-1,j,d)=cy(Istr-1,j,d)*                          &
     &                         GRID(ng)%vmask(Istr-1,j)
#   endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isCY3d,ng)%closed) THEN
          DO d=1,ND
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                cy(Istr-1,j,d)=0.0_r8
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, gradient boundary condition.
!
        IF (LBC(ieast,isCY3d,ng)%gradient) THEN
          DO d=1,ND
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                cy(Iend+1,j,d)=cy(Iend,j,d)
#   ifdef MASKING
                cy(Iend+1,j,d)=cy(Iend+1,j,d)*                          &
     &                         GRID(ng)%vmask(Iend+1,j)
#   endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isCY3d,ng)%closed) THEN
          DO d=1,ND
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                cy(Iend+1,j,d)=0.0_r8
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
            DO d=1,ND
              cy(Istr-1,Jstr,d)=0.5_r8*(cy(Istr  ,Jstr  ,d)+            &
     &                                  cy(Istr-1,Jstr+1,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO d=1,ND
              cy(Iend+1,Jstr,d)=0.5_r8*(cy(Iend  ,Jstr  ,d)+            &
     &                                  cy(Iend+1,Jstr+1,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO d=1,ND
              cy(Istr-1,Jend+1,d)=0.5_r8*(cy(Istr-1,Jend  ,d)+          &
     &                                    cy(Istr  ,Jend+1,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO d=1,ND
              cy(Iend+1,Jend+1,d)=0.5_r8*(cy(Iend+1,Jend  ,d)+          &
     &                                    cy(Iend  ,Jend+1,d))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE cy3dbc_tile
#endif
      END MODULE cy3dbc_mod
