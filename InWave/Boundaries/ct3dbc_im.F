#include "cppdefs.h"
      MODULE ct3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: ct3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  This subroutine sets lateral boundary conditions for the            !
!  celerity in the theta direction                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ct3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ct3dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_vars
      USE mod_inwave_params
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ct3dbc_tile (ng, tile,                                       &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 WAVEP(ng)% ct)
      RETURN
      END SUBROUTINE ct3dbc

!
!***********************************************************************
      SUBROUTINE ct3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ct)
!***********************************************************************
!
      USE mod_inwave_params
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ct(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: ct(LBi:UBi,LBj:UBj,ND)
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
        IF (LBC(iwest,isCT3d,ng)%gradient) THEN
          DO d=1,ND+1
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                ct(Istr-1,j,d)=ct(Istr,j,d)
#   ifdef MASKING
                ct(Istr-1,j,d)=ct(Istr-1,j,d)*                          &
     &                         GRID(ng)%rmask(Istr-1,j)
#   endif
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isCT3d,ng)%closed) THEN
          DO d=1,ND+1
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                ct(Istr-1,j,d)=0.0_r8
#   ifdef MASKING
                ct(Istr-1,j,d)=ct(Istr-1,j,d)*                          &
     &                         GRID(ng)%rmask(Istr-1,j)
#   endif
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
        IF (LBC(ieast,isCT3d,ng)%gradient) THEN
          DO d=1,ND+1
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                ct(Iend+1,j,d)=ct(Iend,j,d)
#   ifdef MASKING
                ct(Iend+1,j,d)=ct(Iend+1,j,d)*                          &
     &                         GRID(ng)%rmask(Iend+1,j)
#   endif
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isCT3d,ng)%closed) THEN
          DO d=1,ND+1
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                ct(Iend+1,j,d)=0.0_r8
#   ifdef MASKING
                ct(Iend+1,j,d)=ct(Iend+1,j,d)*                          &
     &                         GRID(ng)%rmask(Iend+1,j)
#   endif
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
        IF (LBC(isouth,isCT3d,ng)%gradient) THEN
          DO d=1,ND+1
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                ct(i,Jstr-1,d)=ct(i,Jstr,d)
#   ifdef MASKING
                ct(i,Jstr-1,d)=ct(i,Jstr-1,d)*                          &
     &                         GRID(ng)%rmask(i,Jstr-1)
#   endif
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isCT3d,ng)%closed) THEN
          DO d=1,ND+1
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                ct(i,Jstr-1,d)=0.0_r8
#   ifdef MASKING
                ct(i,Jstr-1,d)=ct(i,Jstr-1,d)*                          &
     &                         GRID(ng)%rmask(i,Jstr-1)
#   endif
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
        IF (LBC(inorth,isCT3d,ng)%gradient) THEN
          DO d=1,ND+1
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                ct(i,Jend+1,d)=ct(i,Jend,d)
#   ifdef MASKING
                ct(i,Jend+1,d)=ct(i,Jend+1,d)*                          &
     &                         GRID(ng)%rmask(i,Jend+1)
#   endif
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isCT3d,ng)%closed) THEN
          DO d=1,ND+1
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                ct(i,Jend+1,d)=0.0_r8
#   ifdef MASKING
                ct(i,Jend+1,d)=ct(i,Jend+1,d)*                          &
     &                         GRID(ng)%rmask(i,Jend+1)
#   endif
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
            DO d=1,ND+1
              ct(Istr-1,Jstr-1,d)=0.5_r8*                               &
     &                            (ct(Istr  ,Jstr-1,d)+                 &
     &                             ct(Istr-1,Jstr  ,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO d=1,ND+1
              ct(Iend+1,Jstr-1,d)=0.5_r8*                               &
     &                           (ct(Iend  ,Jstr-1,d)+                  &
     &                            ct(Iend+1,Jstr  ,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO d=1,ND+1
              ct(Istr-1,Jend+1,d)=0.5_r8*                               &
     &                            (ct(Istr-1,Jend  ,d)+                 &
     &                             ct(Istr  ,Jend+1,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO d=1,ND+1
              ct(Iend+1,Jend+1,d)=0.5_r8*                               &
     &                            (ct(Iend+1,Jend  ,d)+                 &
     &                             ct(Iend  ,Jend+1,d))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE ct3dbc_tile
#endif
      END MODULE ct3dbc_mod
