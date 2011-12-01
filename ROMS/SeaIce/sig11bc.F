#include "cppdefs.h"
      MODULE sig11bc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the sig11 ice stress.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC sig11bc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sig11bc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL sig11bc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   liold(ng), linew(ng),                          &
     &                   ICE(ng) % sig11)
      RETURN
      END SUBROUTINE sig11bc
!
!***********************************************************************
      SUBROUTINE sig11bc_tile (ng, tile,                                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew,                          &
     &                           sig11)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars

      implicit none

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: liold, linew

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: sig11(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: sig11(LBi:UBi,LBj:UBj,2)
# endif

!
!  Local variable declarations.
!
      integer :: i, j, know

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
        know=liold

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

# ifdef WEST_SIG11CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          sig11(0,j,linew)=BOUNDARY(ng)%sig11_west(j)
#  ifdef MASKING
          sig11(0,j,linew)=sig11(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sig11(0,j,linew)=sig11(0,j,linew)*                            &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# elif defined WEST_SIG11GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          sig11(0,j,linew)=sig11(1,j,linew)
#  ifdef MASKING
          sig11(0,j,linew)=sig11(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sig11(0,j,linew)=sig11(0,j,linew)*                            &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          sig11(0,j,linew)=sig11(1,j,linew)
#  ifdef MASKING
          sig11(0,j,linew)=sig11(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sig11(0,j,linew)=sig11(0,j,linew)*                            &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

# ifdef EAST_SIG11CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          sig11(Lm(ng)+1,j,linew)=BOUNDARY(ng)%sig11_east(j)
#  ifdef MASKING
          sig11(Lm(ng)+1,j,linew)=sig11(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sig11(Lm(ng)+1,j,linew)=sig11(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_SIG11GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          sig11(Lm(ng)+1,j,linew)=sig11(Lm(ng),j,linew)
#  ifdef MASKING
          sig11(Lm(ng)+1,j,linew)=sig11(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sig11(Lm(ng)+1,j,linew)=sig11(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          sig11(Lm(ng)+1,j,linew)=sig11(Lm(ng),j,linew)
#  ifdef MASKING
          sig11(Lm(ng)+1,j,linew)=sig11(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sig11(Lm(ng)+1,j,linew)=sig11(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# endif
      END IF
#endif
#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

# ifdef SOUTH_SIG11CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          sig11(i,0,linew)=BOUNDARY(ng)%sig11_south(i)
#  ifdef MASKING
          sig11(i,0,linew)=sig11(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sig11(i,0,linew)=sig11(i,0,linew)*                            &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# elif defined SOUTH_SIG11GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          sig11(i,0,linew)=sig11(i,1,linew)
#  ifdef MASKING
          sig11(i,0,linew)=sig11(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sig11(i,0,linew)=sig11(i,0,linew)*                            &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          sig11(i,0,linew)=sig11(i,1,linew)
#  ifdef MASKING
          sig11(i,0,linew)=sig11(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sig11(i,0,linew)=sig11(i,0,linew)*                            &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

# ifdef NORTH_SIG11CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          sig11(i,Mm(ng)+1,linew)=BOUNDARY(ng)%sig11_north(i)
#  ifdef MASKING
          sig11(i,Mm(ng)+1,linew)=sig11(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sig11(i,Mm(ng)+1,linew)=sig11(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_SIG11GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          sig11(i,Mm(ng)+1,linew)=sig11(i,Mm(ng),linew)
#  ifdef MASKING
          sig11(i,Mm(ng)+1,linew)=sig11(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sig11(i,Mm(ng)+1,linew)=sig11(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          sig11(i,Mm(ng)+1,linew)=sig11(i,Mm(ng),linew)
#  ifdef MASKING
          sig11(i,Mm(ng)+1,linew)=sig11(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sig11(i,Mm(ng)+1,linew)=sig11(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# endif
      END IF
#endif
#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE.and.WESTERN_EDGE) THEN
        sig11(0,0,linew)=0.5_r8*(sig11(1,0,linew)+                      &
     &                         sig11(0,1,linew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        sig11(Lm(ng)+1,0,linew)=0.5_r8*(sig11(Lm(ng)+1,1,linew)+        &
     &                                sig11(Lm(ng)  ,0,linew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        sig11(0,Mm(ng)+1,linew)=0.5_r8*(sig11(0,Mm(ng)  ,linew)+        &
     &                                sig11(1,Mm(ng)+1,linew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        sig11(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                          &
     &             (sig11(Lm(ng)+1,Mm(ng)  ,linew)+                     &
     &              sig11(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE sig11bc_tile
#endif

      END MODULE sig11bc_mod
