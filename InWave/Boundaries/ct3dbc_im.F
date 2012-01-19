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
      USE mod_inwave_params
      USE mod_inwave_vars
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
     &                 LBi, UBi, LBj, UBj, ND(ng)+1,                    &
     &                 WAVEP(ng)% ct)
      RETURN
      END SUBROUTINE ct3dbc

!
!***********************************************************************
      SUBROUTINE ct3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj, ND,                    &
     &                       ct)
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
      real(r8), intent(inout) :: ct(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: ct(LBi:UBi,LBj:UBj,ND)
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

#  if defined WEST_CT_GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO d=1,ND
         DO j=Jstr,Jend
          ct(Istr-1,j,d)=ct(Istr,j,d)
#   ifdef MASKING
          ct(Istr-1,j,d)=ct(Istr-1,j,d)*                                &
     &                   GRID(ng)%rmask(Istr-1,j)
#   endif
         END DO
        END DO

#  elif defined WEST_CT_WALL
!
!  Western edge, closed boundary condition.
!
        DO d=1,ND
         DO j=Jstr,Jend
          ct(Istr-1,j,d)=0.0_r8
#   ifdef MASKING
          ct(Istr-1,j,d)=ct(Istr-1,j,d)*                                &
     &                   GRID(ng)%rmask(Istr-1,j)
#   endif
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

#  if defined EAST_CT_GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO d=1,ND
         DO j=Jstr,Jend
          ct(Iend+1,j,d)=ct(Iend,j,d)
#   ifdef MASKING
          ct(Iend+1,j,d)=ct(Iend+1,j,d)*                                &
     &                   GRID(ng)%rmask(Iend+1,j)
#   endif
         END DO
        END DO

#  elif defined EAST_CT_WALL
!
!  Eastern edge, closed boundary condition.
!
        DO d=1,ND
         DO j=Jstr,Jend
          ct(Iend+1,j,d)=0.0_r8
#   ifdef MASKING
          ct(Iend+1,j,d)=ct(Iend+1,j,d)*                                &
     &                   GRID(ng)%rmask(Iend+1,j)
#   endif
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

#  if defined SOUTH_CT_GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO d=1,ND
         DO i=Istr,Iend
          ct(i,Jstr-1,d)=ct(i,Jstr,d)
#   ifdef MASKING
          ct(i,Jstr-1,d)=ct(i,Jstr-1,d)*                                &
     &                   GRID(ng)%rmask(i,Jstr-1)
#   endif
         END DO
        END DO

#  elif defined SOUTH_CT_WALL
!
!  Southern edge, closed boundary condition.
!
        DO d=1,ND
         DO i=Istr,Iend
          ct(i,Jstr-1,d)=0.0_r8
#   ifdef MASKING
          ct(i,Jstr-1,d)=ct(i,Jstr-1,d)*                                &
     &                   GRID(ng)%rmask(i,Jstr-1)
#   endif
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

#  if defined NORTH_CT_GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO d=1,ND
         DO i=Istr,Iend
          ct(i,Jend+1,d)=ct(i,Jend,d)
#   ifdef MASKING
          ct(i,Jend+1,d)=ct(i,Jend+1,d)*                                &
     &                   GRID(ng)%rmask(i,Jend+1)
#   endif
         END DO
        END DO

#  elif defined NORTH_CT_WALL
!
!  Northern edge, closed boundary condition.
!
        DO d=1,ND
         DO i=Istr,Iend
          ct(i,Jend+1,d)=0.0_r8
#   ifdef MASKING
          ct(i,Jend+1,d)=ct(i,Jend+1,d)*                                &
     &                   GRID(ng)%rmask(i,Jend+1)
#   endif
         END DO
        END DO
#  endif
       END IF
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
       IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO d=1,ND
         ct(Istr-1,Jstr-1,d)=0.5_r8*                                    &
     &                       (ct(Istr  ,Jstr-1,d)+                      &
     &                        ct(Istr-1,Jstr  ,d))
        END DO
       END IF
       IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND
         ct(Iend+1,Jstr-1,d)=0.5_r8*                                    &
     &                       (ct(Iend  ,Jstr-1,d)+                      &
     &                        ct(Iend+1,Jstr  ,d))
        END DO
       END IF
       IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO d=1,ND
         ct(Istr-1,Jend+1,d)=0.5_r8*                                    &
     &                       (ct(Istr-1,Jend  ,d)+                      &
     &                        ct(Istr  ,Jend+1,d))
        END DO
       END IF
       IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND
         ct(Iend+1,Jend+1,d)=0.5_r8*                                    &
     &                       (ct(Iend+1,Jend  ,d)+                      &
     &                        ct(Iend  ,Jend+1,d))
        END DO
       END IF
# endif

# ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
# endif

      RETURN
      END SUBROUTINE ct3dbc_tile
#endif
      END MODULE ct3dbc_mod
