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
      SUBROUTINE cy3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj, ND,                    &
     &                       cy)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, ND
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

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      IF (ng.eq.1) THEN
#endif

# if !defined NS_C_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
       IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_CY_GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO d=1,ND
         DO i=Istr,Iend
          cy(i,Jstr,d)=cy(i,Jstr+1,d)
#   ifdef MASKING
          cy(i,Jstr,d)=cy(i,Jstr,d)*                                    &
     &                 GRID(ng)%vmask(i,Jstr)
#   endif
         END DO
        END DO

#  elif defined SOUTH_CY_WALL
!
!  Southern edge, closed boundary condition.
!
        DO d=1,ND
         DO i=Istr,Iend
          cy(i,Jstr,d)=0.0_r8
         END DO
        END DO
#  endif
       END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
       IF (NORTHERN_EDGE) THEN

#  if defined NORTH_CY_GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO d=1,ND
         DO i=Istr,Iend
          cy(i,Jend+1,d)=cy(i,Jend,d)
#   ifdef MASKING
          cy(i,Jend+1,d)=cy(i,Jend+1,d)*                                &
     &                   GRID(ng)%vmask(i,Jend+1)
#   endif
         END DO
        END DO

#  elif defined NORTH_CY_WALL
!
!  Northern edge, closed boundary condition.
!
        DO d=1,ND
         DO i=Istr,Iend
          cy(i,Jend+1,d)=0.0_r8
         END DO
        END DO
#  endif
       END IF
# endif

# if !defined EW_C_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
       IF (WESTERN_EDGE) THEN

#  if defined WEST_CY_GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO d=1,ND
         DO j=JstrV,Jend
          cy(Istr-1,j,d)=cy(Istr,j,d)
#   ifdef MASKING
          cy(Istr-1,j,d)=cy(Istr-1,j,d)*                                &
     &                   GRID(ng)%vmask(Istr-1,j)
#   endif
         END DO
        END DO

#  elif defined WEST_CY_WALL
!
!  Western edge, closed boundary condition.
!
        DO d=1,ND
         DO j=JstrV,Jend
          cy(Istr-1,j,d)=0.0_r8
         END DO
        END DO
#  endif
       END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
       IF (EASTERN_EDGE) THEN

#  if defined EAST_CY_GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO d=1,ND
         DO j=JstrV,Jend
          cy(Iend+1,j,d)=cy(Iend,j,d)
#   ifdef MASKING
          cy(Iend+1,j,d)=cy(Iend+1,j,d)*                                &
     &                   GRID(ng)%vmask(Iend+1,j)
#   endif
         END DO
        END DO

#  elif defined EAST_CY_WALL
!
!  Western edge, closed boundary condition.
!
        DO d=1,ND
         DO j=JstrV,Jend
          cy(Iend+1,j,d)=0.0_r8
         END DO
        END DO
#  endif
       END IF
# endif

# if !defined EW_C_PERIODIC && !defined NS_C_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
       IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO d=1,ND
         cy(Istr-1,Jstr,d)=0.5_r8*(cy(Istr  ,Jstr  ,d)+                 &
     &                             cy(Istr-1,Jstr+1,d))
        END DO
       END IF
       IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND
         cy(Iend+1,Jstr,d)=0.5_r8*(cy(Iend  ,Jstr  ,d)+                 &
     &                             cy(Iend+1,Jstr+1,d))
        END DO
       END IF
       IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO d=1,ND
         cy(Istr-1,Jend+1,d)=0.5_r8*(cy(Istr-1,Jend  ,d)+               &
     &                               cy(Istr  ,Jend+1,d))
        END DO
       END IF
       IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND
         cy(Iend+1,Jend+1,d)=0.5_r8*(cy(Iend+1,Jend  ,d)+               &
     &                               cy(Iend  ,Jend+1,d))
        END DO
       END IF
# endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif

      RETURN
      END SUBROUTINE cy3dbc_tile
#endif
      END MODULE cy3dbc_mod
