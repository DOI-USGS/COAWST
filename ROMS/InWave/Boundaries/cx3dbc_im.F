#include "cppdefs.h"
      MODULE cx3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: cx3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  This subroutine sets lateral boundary conditions for 3D             !
!  xi component of the wave group velocity                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: cx3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE cx3dbc (ng, tile)
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
      CALL cx3dbc_tile (ng, tile,                                       &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 WAVEP(ng) % cx)
      RETURN
      END SUBROUTINE cx3dbc
!
!***********************************************************************
      SUBROUTINE cx3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       cx)
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
      real(r8), intent(inout) :: cx(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: cx(LBi:UBi,LBj:UBj,ND)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, d

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, gradient boundary condition.
!
        IF (LBC(iwest,isCX3d,ng)%gradient) THEN
          DO d=1,ND
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                cx(Istr,j,d)=cx(Istr+1,j,d)
#   ifdef MASKING
                cx(Istr,j,d)=cx(Istr,j,d)*                              &
     &                       GRID(ng)%umask(Istr,j)
#   endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isCX3d,ng)%closed) THEN
          DO d=1,ND
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                cx(Istr,j,d)=0.0_r8
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
        IF (LBC(ieast,isCX3d,ng)%gradient) THEN
          DO d=1,ND
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                cx(Iend+1,j,d)=cx(Iend,j,d)
#   ifdef MASKING
                cx(Iend+1,j,d)=cx(Iend+1,j,d)*                          &
     &                         GRID(ng)%umask(Iend+1,j)
#   endif
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isCX3d,ng)%closed) THEN
          DO d=1,ND
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                cx(Iend+1,j,d)=0.0_r8
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, gradient boundary condition.
!
        IF (LBC(isouth,isCX3d,ng)%gradient) THEN
          DO d=1,ND
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                cx(i,Jstr-1,d)=cx(i,Jstr,d)
#   ifdef MASKING
                cx(i,Jstr-1,d)=cx(i,Jstr-1,d)*                          &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isCX3d,ng)%closed) THEN
          DO d=1,ND
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                cx(i,Jstr-1,d)=0.0_r8
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
        IF (LBC(inorth,isCX3d,ng)%gradient) THEN
          DO d=1,ND
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                cx(i,Jend+1,d)=cx(i,Jend,d)
#   ifdef MASKING
                cx(i,Jend+1,d)=cx(i,Jend+1,d)*                          &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isCX3d,ng)%closed) THEN
          DO d=1,ND
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                cx(i,Jend+1,d)=0.0_r8
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
              cx(Istr,Jstr-1,d)=0.5_r8*(cx(Istr+1,Jstr-1,d)+            &
     &                                  cx(Istr  ,Jstr  ,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO d=1,ND
              cx(Iend+1,Jstr-1,d)=0.5_r8*(cx(Iend  ,Jstr-1,d)+          &
     &                                    cx(Iend+1,Jstr  ,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO d=1,ND
              cx(Istr,Jend+1,d)=0.5_r8*(cx(Istr  ,Jend  ,d)+            &
     &                                  cx(Istr+1,Jend+1,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO d=1,ND
              cx(Iend+1,Jend+1,d)=0.5_r8*(cx(Iend+1,Jend  ,d)+          &
     &                                    cx(Iend  ,Jend+1,d))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE cx3dbc_tile
#endif
      END MODULE cx3dbc_mod
