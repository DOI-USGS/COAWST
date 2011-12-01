#include "cppdefs.h"
      MODULE sig12bc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the sig12 ice stress.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC sig12bc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sig12bc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL sig12bc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   liold(ng), linew(ng),                          &
     &                   ICE(ng) % sig12)
      RETURN
      END SUBROUTINE sig12bc

!
!***********************************************************************
      SUBROUTINE sig12bc_tile (ng, tile,                                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew,                          &
     &                           sig12)
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
      real(r8), intent(inout) :: sig12(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: sig12(LBi:UBi,LBj:UBj,2)
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

# ifdef WEST_SIG12CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          sig12(0,j,linew)=BOUNDARY(ng)%sig12_west(j)
#  ifdef MASKING
          sig12(0,j,linew)=sig12(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sig12(0,j,linew)=sig12(0,j,linew)*                            &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# elif defined WEST_SIG12GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          sig12(0,j,linew)=sig12(1,j,linew)
#  ifdef MASKING
          sig12(0,j,linew)=sig12(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sig12(0,j,linew)=sig12(0,j,linew)*                            &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          sig12(0,j,linew)=sig12(1,j,linew)
#  ifdef MASKING
          sig12(0,j,linew)=sig12(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sig12(0,j,linew)=sig12(0,j,linew)*                            &
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

# ifdef EAST_SIG12CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          sig12(Lm(ng)+1,j,linew)=BOUNDARY(ng)%sig12_east(j)
#  ifdef MASKING
          sig12(Lm(ng)+1,j,linew)=sig12(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sig12(Lm(ng)+1,j,linew)=sig12(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_SIG12GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          sig12(Lm(ng)+1,j,linew)=sig12(Lm(ng),j,linew)
#  ifdef MASKING
          sig12(Lm(ng)+1,j,linew)=sig12(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sig12(Lm(ng)+1,j,linew)=sig12(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          sig12(Lm(ng)+1,j,linew)=sig12(Lm(ng),j,linew)
#  ifdef MASKING
          sig12(Lm(ng)+1,j,linew)=sig12(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sig12(Lm(ng)+1,j,linew)=sig12(Lm(ng)+1,j,linew)*              &
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

# ifdef SOUTH_SIG12CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          sig12(i,0,linew)=BOUNDARY(ng)%sig12_south(i)
#  ifdef MASKING
          sig12(i,0,linew)=sig12(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sig12(i,0,linew)=sig12(i,0,linew)*                            &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# elif defined SOUTH_SIG12GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          sig12(i,0,linew)=sig12(i,1,linew)
#  ifdef MASKING
          sig12(i,0,linew)=sig12(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sig12(i,0,linew)=sig12(i,0,linew)*                            &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          sig12(i,0,linew)=sig12(i,1,linew)
#  ifdef MASKING
          sig12(i,0,linew)=sig12(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sig12(i,0,linew)=sig12(i,0,linew)*                            &
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

# ifdef NORTH_SIG12CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          sig12(i,Mm(ng)+1,linew)=BOUNDARY(ng)%sig12_north(i)
#  ifdef MASKING
          sig12(i,Mm(ng)+1,linew)=sig12(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sig12(i,Mm(ng)+1,linew)=sig12(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_SIG12GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          sig12(i,Mm(ng)+1,linew)=sig12(i,Mm(ng),linew)
#  ifdef MASKING
          sig12(i,Mm(ng)+1,linew)=sig12(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sig12(i,Mm(ng)+1,linew)=sig12(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          sig12(i,Mm(ng)+1,linew)=sig12(i,Mm(ng),linew)
#  ifdef MASKING
          sig12(i,Mm(ng)+1,linew)=sig12(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sig12(i,Mm(ng)+1,linew)=sig12(i,Mm(ng)+1,linew)*              &
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
        sig12(0,0,linew)=0.5_r8*(sig12(1,0,linew)+                      &
     &                         sig12(0,1,linew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        sig12(Lm(ng)+1,0,linew)=0.5_r8*(sig12(Lm(ng)+1,1,linew)+        &
     &                                sig12(Lm(ng)  ,0,linew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        sig12(0,Mm(ng)+1,linew)=0.5_r8*(sig12(0,Mm(ng)  ,linew)+        &
     &                                sig12(1,Mm(ng)+1,linew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        sig12(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                          &
     &             (sig12(Lm(ng)+1,Mm(ng)  ,linew)+                     &
     &              sig12(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE sig12bc_tile
#endif

      END MODULE sig12bc_mod
