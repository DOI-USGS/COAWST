#include "cppdefs.h"
      MODULE sig22bc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the sig22 ice stress.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC sig22bc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sig22bc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL sig22bc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   liold(ng), linew(ng),                          &
     &                   ICE(ng) % sig22)
      RETURN
      END SUBROUTINE sig22bc

!
!***********************************************************************
      SUBROUTINE sig22bc_tile (ng, tile,                                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew,                          &
     &                           sig22)
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
      real(r8), intent(inout) :: sig22(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: sig22(LBi:UBi,LBj:UBj,2)
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

# ifdef WEST_SIG22CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          sig22(0,j,linew)=BOUNDARY(ng)%sig22_west(j)
#  ifdef MASKING
          sig22(0,j,linew)=sig22(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sig22(0,j,linew)=sig22(0,j,linew)*                            &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# elif defined WEST_SIG22GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          sig22(0,j,linew)=sig22(1,j,linew)
#  ifdef MASKING
          sig22(0,j,linew)=sig22(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sig22(0,j,linew)=sig22(0,j,linew)*                            &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          sig22(0,j,linew)=sig22(1,j,linew)
#  ifdef MASKING
          sig22(0,j,linew)=sig22(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sig22(0,j,linew)=sig22(0,j,linew)*                            &
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

# ifdef EAST_SIG22CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          sig22(Lm(ng)+1,j,linew)=BOUNDARY(ng)%sig22_east(j)
#  ifdef MASKING
          sig22(Lm(ng)+1,j,linew)=sig22(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sig22(Lm(ng)+1,j,linew)=sig22(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_SIG22GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          sig22(Lm(ng)+1,j,linew)=sig22(Lm(ng),j,linew)
#  ifdef MASKING
          sig22(Lm(ng)+1,j,linew)=sig22(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sig22(Lm(ng)+1,j,linew)=sig22(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          sig22(Lm(ng)+1,j,linew)=sig22(Lm(ng),j,linew)
#  ifdef MASKING
          sig22(Lm(ng)+1,j,linew)=sig22(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sig22(Lm(ng)+1,j,linew)=sig22(Lm(ng)+1,j,linew)*              &
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

# ifdef SOUTH_SIG22CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          sig22(i,0,linew)=BOUNDARY(ng)%sig22_south(i)
#  ifdef MASKING
          sig22(i,0,linew)=sig22(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sig22(i,0,linew)=sig22(i,0,linew)*                            &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# elif defined SOUTH_SIG22GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          sig22(i,0,linew)=sig22(i,1,linew)
#  ifdef MASKING
          sig22(i,0,linew)=sig22(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sig22(i,0,linew)=sig22(i,0,linew)*                            &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          sig22(i,0,linew)=sig22(i,1,linew)
#  ifdef MASKING
          sig22(i,0,linew)=sig22(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sig22(i,0,linew)=sig22(i,0,linew)*                            &
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

# ifdef NORTH_SIG22CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          sig22(i,Mm(ng)+1,linew)=BOUNDARY(ng)%sig22_north(i)
#  ifdef MASKING
          sig22(i,Mm(ng)+1,linew)=sig22(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sig22(i,Mm(ng)+1,linew)=sig22(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_SIG22GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          sig22(i,Mm(ng)+1,linew)=sig22(i,Mm(ng),linew)
#  ifdef MASKING
          sig22(i,Mm(ng)+1,linew)=sig22(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sig22(i,Mm(ng)+1,linew)=sig22(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          sig22(i,Mm(ng)+1,linew)=sig22(i,Mm(ng),linew)
#  ifdef MASKING
          sig22(i,Mm(ng)+1,linew)=sig22(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sig22(i,Mm(ng)+1,linew)=sig22(i,Mm(ng)+1,linew)*              &
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
        sig22(0,0,linew)=0.5_r8*(sig22(1,0,linew)+                      &
     &                         sig22(0,1,linew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        sig22(Lm(ng)+1,0,linew)=0.5_r8*(sig22(Lm(ng)+1,1,linew)+        &
     &                                sig22(Lm(ng)  ,0,linew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        sig22(0,Mm(ng)+1,linew)=0.5_r8*(sig22(0,Mm(ng)  ,linew)+        &
     &                                sig22(1,Mm(ng)+1,linew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        sig22(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                          &
     &             (sig22(Lm(ng)+1,Mm(ng)  ,linew)+                     &
     &              sig22(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE sig22bc_tile
#endif

      END MODULE sig22bc_mod
