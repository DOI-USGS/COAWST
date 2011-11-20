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
     &                 LBi, UBi, LBj, UBj, ND(ng),                      &
     &                 WAVEP(ng) % cx)
      RETURN
      END SUBROUTINE cx3dbc
!
!***********************************************************************
      SUBROUTINE cx3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj, ND,                    &
     &                       cx)
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
      real(r8), intent(inout) :: cx(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: cx(LBi:UBi,LBj:UBj,ND)
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

# if !defined EW_C_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
       IF (WESTERN_EDGE) THEN

#  if defined WEST_CX_GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO d=1,ND
         DO j=Jstr,Jend
          cx(Istr,j,d)=cx(Istr+1,j,d)
#   ifdef MASKING
          cx(Istr,j,d)=cx(Istr,j,d)*                                    &
     &                 GRID(ng)%umask(Istr,j)
#   endif
         END DO
        END DO

#  elif defined WEST_CX_WALL
!
!  Western edge, closed boundary condition.
!
        DO d=1,ND
         DO j=Jstr,Jend
          cx(Istr,j,d)=0.0_r8
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

#  if defined EAST_CX_GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO d=1,ND
         DO j=Jstr,Jend
          cx(Iend+1,j,d)=cx(Iend,j,d)
#   ifdef MASKING
          cx(Iend+1,j,d)=cx(Iend+1,j,d)*                                &
     &                   GRID(ng)%umask(Iend+1,j)
#   endif
         END DO
        END DO

#  elif defined EAST_CX_WALL
!
!  Eastern edge, closed boundary condition.
!
        DO d=1,ND
         DO j=Jstr,Jend
          cx(Iend+1,j,d)=0.0_r8
         END DO
        END DO
#  endif
       END IF
# endif

# if !defined NS_C_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
       IF (SOUTHERN_EDGE) THEN


#  if defined SOUTH_CX_GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO d=1,ND
         DO i=IstrU,Iend
          cx(i,Jstr-1,d)=cx(i,Jstr,d)
#   ifdef MASKING
          cx(i,Jstr-1,d)=cx(i,Jstr-1,d)*                                &
     &                   GRID(ng)%umask(i,Jstr-1)
#   endif
         END DO
        END DO

#  elif defined SOUTH_CX_WALL
!
!  Southern edge, closed boundary condition.
!
        DO d=1,ND
         DO i=IstrU,Iend
          cx(i,Jstr-1,d)=0.0_r8
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

#  if defined NORTH_CX_GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO d=1,ND
         DO i=IstrU,Iend
          cx(i,Jend+1,d)=cx(i,Jend,d)
#   ifdef MASKING
          cx(i,Jend+1,d)=cx(i,Jend+1,d)*                                &
     &                   GRID(ng)%umask(i,Jend+1)
#   endif
         END DO
        END DO

#  elif defined NORTH_CX_WALL
!
!  Northern edge, closed boundary condition.
!
        DO d=1,ND
         DO i=IstrU,Iend
          cx(i,Jend+1,d)=0.0_r8
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
         cx(Istr,Jstr-1,d)=0.5_r8*(cx(Istr+1,Jstr-1,d)+                 &
     &                             cx(Istr  ,Jstr  ,d))
        END DO
       END IF
       IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND
         cx(Iend+1,Jstr-1,d)=0.5_r8*(cx(Iend  ,Jstr-1,d)+               &
     &                               cx(Iend+1,Jstr  ,d))
        END DO
       END IF
       IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO d=1,ND
         cx(Istr,Jend+1,d)=0.5_r8*(cx(Istr  ,Jend  ,d)+                 &
     &                             cx(Istr+1,Jend+1,d))
        END DO
       END IF
       IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND
         cx(Iend+1,Jend+1,d)=0.5_r8*(cx(Iend+1,Jend  ,d)+               &
     &                               cx(Iend  ,Jend+1,d))
        END DO
       END IF
# endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif

      RETURN
      END SUBROUTINE cx3dbc_tile
#endif
      END MODULE cx3dbc_mod
