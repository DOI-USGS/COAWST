#include "cppdefs.h"
      MODULE aibc_mod
#ifdef ICE_MODEL
! 
!***********************************************************************
!  Compute lateral boundary conditions for the ice concentration.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC aibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE aibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL aibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold(ng), linew(ng),                             &
     &                ICE(ng) % ui,                                     &
     &                ICE(ng) % vi,                                     &
     &                ICE(ng) % ai)
      RETURN
      END SUBROUTINE aibc

!
!***********************************************************************
      SUBROUTINE aibc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      liold, linew,                               &
     &                      ui, vi, ai)
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
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
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

# ifdef WEST_AICLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(1,j,linew).ge.0._r8) THEN
             ai(0,j,linew)=BOUNDARY(ng)%ai_west(j)
#  ifdef MASKING
          ai(0,j,linew)=ai(0,j,linew)*                               &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          ai(0,j,linew)=ai(0,j,linew)*                               &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
             ai(0,j,linew)=ai(1,j,liold)
#  ifdef MASKING
             ai(0,j,linew)=ai(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             ai(0,j,linew)=ai(0,j,linew)*                            &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ENDIF
        END DO
# elif defined WEST_AIGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ai(0,j,linew)=ai(1,j,linew)
#  ifdef MASKING
          ai(0,j,linew)=ai(0,j,linew)*                                  &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          ai(0,j,linew)=ai(0,j,linew)*                                  &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ai(0,j,linew)=ai(1,j,linew)
#  ifdef MASKING
          ai(0,j,linew)=ai(0,j,linew)*                                  &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          ai(0,j,linew)=ai(0,j,linew)*                                  &
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

# ifdef EAST_AICLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(Lm(ng)+1,j,linew).le.0._r8) THEN
             ai(Lm(ng)+1,j,linew)=BOUNDARY(ng)%ai_east(j)
#  ifdef MASKING
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,liold)
#  ifdef MASKING
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ENDIF
        END DO
# elif defined EAST_AIGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,linew)
#  ifdef MASKING
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,linew)
#  ifdef MASKING
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
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

# ifdef SOUTH_AICLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,1,linew).ge.0._r8) THEN
             ai(i,0,linew)=BOUNDARY(ng)%ai_south(i)
#  ifdef MASKING
             ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
          ai(i,0,linew)=ai(i,1,liold)
#  ifdef MASKING
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ENDIF
        END DO
# elif defined SOUTH_AIGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          ai(i,0,linew)=ai(i,1,linew)
#  ifdef MASKING
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          ai(i,0,linew)=ai(i,1,linew)
#  ifdef MASKING
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          ai(i,0,linew)=ai(i,0,linew)*                                  &
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

# ifdef NORTH_AICLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,Mm(ng)+1,linew).le.0._r8) THEN
             ai(i,Mm(ng)+1,linew)=BOUNDARY(ng)%ai_north(i)
#  ifdef MASKING
             ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ELSE
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),liold)
#  ifdef MASKING
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ENDIF
        END DO
# elif defined NORTH_AIGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),linew)
#  ifdef MASKING
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),linew)
#  ifdef MASKING
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
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
        ai(0,0,linew)=0.5_r8*(ai(1,0,linew)+                            &
     &                         ai(0,1,linew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        ai(Lm(ng)+1,0,linew)=0.5_r8*(ai(Lm(ng)+1,1,linew)+              &
     &                                ai(Lm(ng)  ,0,linew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        ai(0,Mm(ng)+1,linew)=0.5_r8*(ai(0,Mm(ng)  ,linew)+              &
     &                                ai(1,Mm(ng)+1,linew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        ai(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                             &
     &             (ai(Lm(ng)+1,Mm(ng)  ,linew)+                        &
     &              ai(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE aibc_tile
#endif

      END MODULE aibc_mod
