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

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      IF (ng.eq.1) THEN
#endif

!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
       IF (WESTERN_EDGE) THEN
        DO d=1,ND
         DO j=Jstr,Jend
          Tr(Istr-1,j,d)=Tr(Istr,j,d)
#   ifdef MASKING
          Tr(Istr-1,j,d)=Tr(Istr-1,j,d)*                                &
     &                   GRID(ng)%rmask(Istr-1,j)
#   endif
         END DO
        END DO
       END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
       IF (EASTERN_EDGE) THEN
        DO d=1,ND
         DO j=Jstr,Jend
          Tr(Iend+1,j,d)=Tr(Iend,j,d)
#   ifdef MASKING
          Tr(Iend+1,j,d)=Tr(Iend+1,j,d)*                                &
     &                   GRID(ng)%rmask(Iend+1,j)
#   endif
         END DO
        END DO
       END IF

!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
       IF (SOUTHERN_EDGE) THEN
        DO d=1,ND
         DO i=Istr,Iend
          Tr(i,Jstr-1,d)=Tr(i,Jstr,d)
#   ifdef MASKING
          Tr(i,Jstr-1,d)=Tr(i,Jstr-1,d)*                                &
     &                   GRID(ng)%rmask(i,Jstr-1)
#   endif
         END DO
        END DO
       END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
       IF (NORTHERN_EDGE) THEN
        DO d=1,ND
         DO i=Istr,Iend
          Tr(i,Jend+1,d)=Tr(i,Jend,d)
#   ifdef MASKING
          Tr(i,Jend+1,d)=Tr(i,Jend+1,d)*                                &
     &                   GRID(ng)%rmask(i,Jend+1)
#   endif
         END DO
        END DO
       END IF


!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
       IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO d=1,ND
         Tr(Istr-1,Jstr-1,d)=0.5_r8*                                    &
     &                       (Tr(Istr,Jstr-1,d)+                        &
     &                        Tr(Istr-1,Jstr,d))
        END DO
       END IF
       IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND
         Tr(Iend+1,Jstr-1,d)=0.5_r8*                                    &
     &                       (Tr(Iend  ,Jstr-1,d)+                      &
     &                        Tr(Iend+1,Jstr  ,d))
        END DO
       END IF
       IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO d=1,ND
         Tr(Istr-1,Jend+1,d)=0.5_r8*                                    &
     &                       (Tr(Istr-1,Jend  ,d)+                      &
     &                        Tr(Istr  ,Jend+1,d))
        END DO
       END IF
       IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND
         Tr(Iend+1,Jend+1,d)=0.5_r8*                                    &
     &                       (Tr(Iend+1,Jend  ,d)+                      &
     &                        Tr(Iend  ,Jend+1,d))
        END DO
       END IF

# ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
# endif

      RETURN
      END SUBROUTINE Tr3dbc_tile
#endif
      END MODULE Tr3dbc_mod
