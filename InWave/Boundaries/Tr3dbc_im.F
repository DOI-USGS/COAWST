#include "cppdefs.h"
      MODULE Tr3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: Tr3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!                                                                      !
!  This subroutine sets lateral boundary conditions for the ITRC-th    !
!  tracer field.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: Tr3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE Tr3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj, ND,                    &
     &                       Tr)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
# ifdef REFINED_GRID
      USE mod_stepping
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, ND
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Tr(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: Tr(LBi:UBi,LBj:UBj,ND)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, d

# include "set_bounds.h"

!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO d=1,ND
            DO j=Jstr,Jend
              Tr(Istr-1,j,d)=Tr(Istr,j,d)
#   ifdef MASKING
              Tr(Istr-1,j,d)=Tr(Istr-1,j,d)*                            &
     &                       GRID(ng)%rmask(Istr-1,j)
#   endif
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO d=1,ND
            DO j=Jstr,Jend
              Tr(Iend+1,j,d)=Tr(Iend,j,d)
#   ifdef MASKING
              Tr(Iend+1,j,d)=Tr(Iend+1,j,d)*                            &
     &                       GRID(ng)%rmask(Iend+1,j)
#   endif
            END DO
          END DO
        END IF
      END IF
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO d=1,ND
            DO i=Istr,Iend
              Tr(i,Jstr-1,d)=Tr(i,Jstr,d)
#   ifdef MASKING
              Tr(i,Jstr-1,d)=Tr(i,Jstr-1,d)*                            &
     &                       GRID(ng)%rmask(i,Jstr-1)
#   endif
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO d=1,ND
            DO i=Istr,Iend
              Tr(i,Jend+1,d)=Tr(i,Jend,d)
#   ifdef MASKING
              Tr(i,Jend+1,d)=Tr(i,Jend+1,d)*                            &
     &                       GRID(ng)%rmask(i,Jend+1)
#   endif
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
              Tr(Istr-1,Jstr-1,d)=0.5_r8*                               &
     &                            (Tr(Istr,Jstr-1,d)+                   &
     &                             Tr(Istr-1,Jstr,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO d=1,ND
              Tr(Iend+1,Jstr-1,d)=0.5_r8*                               &
     &                            (Tr(Iend  ,Jstr-1,d)+                 &
     &                             Tr(Iend+1,Jstr  ,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO d=1,ND
              Tr(Istr-1,Jend+1,d)=0.5_r8*                               &
     &                            (Tr(Istr-1,Jend  ,d)+                 &
     &                             Tr(Istr  ,Jend+1,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO d=1,ND
              Tr(Iend+1,Jend+1,d)=0.5_r8*                               &
     &                            (Tr(Iend+1,Jend  ,d)+                 &
     &                             Tr(Iend  ,Jend+1,d))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE Tr3dbc_tile
#endif
      END MODULE Tr3dbc_mod
