#include "cppdefs.h"
      MODULE hibc_mod
#ifdef ICE_MODEL
!
!***********************************************************************
!  Compute lateral boundary conditions for the ice thickness.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC hibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE hibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL hibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold(ng), linew(ng),                             &
     &                ICE(ng) % ui,                                     &
     &                ICE(ng) % vi,                                     &
     &                ICE(ng) % hi)
      RETURN
      END SUBROUTINE hibc

!
!***********************************************************************
      SUBROUTINE hibc_tile (ng, tile,                                   &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew,                          &
     &                           ui, vi, hi)
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
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj,2)
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

# ifdef WEST_HICLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(1,j,linew).ge.0._r8) THEN
             hi(0,j,linew)=BOUNDARY(ng)%hi_west(j)
#  ifdef MASKING
             hi(0,j,linew)=hi(0,j,linew)*                               &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             hi(0,j,linew)=hi(0,j,linew)*                               &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
             hi(0,j,linew)=hi(1,j,liold)
#  ifdef MASKING
             hi(0,j,linew)=hi(0,j,linew)*                               &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             hi(0,j,linew)=hi(0,j,linew)*                               &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ENDIF
        END DO
# elif defined WEST_HIGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          hi(0,j,linew)=hi(1,j,linew)
#  ifdef MASKING
          hi(0,j,linew)=hi(0,j,linew)*                                  &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hi(0,j,linew)=hi(0,j,linew)*                                  &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          hi(0,j,linew)=hi(1,j,linew)
#  ifdef MASKING
          hi(0,j,linew)=hi(0,j,linew)*                                  &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hi(0,j,linew)=hi(0,j,linew)*                                  &
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

# ifdef EAST_HICLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(Lm(ng)+1,j,linew).le.0._r8) THEN
             hi(Lm(ng)+1,j,linew)=BOUNDARY(ng)%hi_east(j)
#  ifdef MASKING
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,liold)
#  ifdef MASKING
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ENDIF
        END DO
# elif defined EAST_HIGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)
#  ifdef MASKING
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask(Lm(ng),j)
#  endif
#  ifdef WET_DRY
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask_wet(Lm(ng),j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)
#  ifdef MASKING
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask(Lm(ng),j)
#  endif
#  ifdef WET_DRY
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask_wet(Lm(ng),j)
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

# ifdef SOUTH_HICLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,1,linew).ge.0._r8) THEN
             hi(i,0,linew)=BOUNDARY(ng)%hi_south(i)
#  ifdef MASKING
             hi(i,0,linew)=hi(i,0,linew)*                               &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             hi(i,0,linew)=hi(i,0,linew)*                               &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
             hi(i,0,linew)=hi(i,1,linew)
#  ifdef MASKING
             hi(i,0,linew)=hi(i,0,linew)*                               &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             hi(i,0,linew)=hi(i,0,linew)*                               &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ENDIF
        END DO
# elif defined SOUTH_HIGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          hi(i,0,linew)=hi(i,1,linew)
#  ifdef MASKING
          hi(i,0,linew)=hi(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hi(i,0,linew)=hi(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          hi(i,0,linew)=hi(i,1,linew)
#  ifdef MASKING
          hi(i,0,linew)=hi(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hi(i,0,linew)=hi(i,0,linew)*                                  &
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

# ifdef NORTH_HICLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,Mm(ng)+1,linew).le.0._r8) THEN
             hi(i,Mm(ng)+1,linew)=BOUNDARY(ng)%hi_north(i)
#  ifdef MASKING
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
           ELSE
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),liold)
#  ifdef MASKING
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
           ENDIF
        END DO
# elif defined NORTH_HIGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)
#  ifdef WET_DRY
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)
#  ifdef MASKING
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                    &
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
        hi(0,0,linew)=0.5_r8*(hi(1,0,linew)+                            &
     &                         hi(0,1,linew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        hi(Lm(ng)+1,0,linew)=0.5_r8*(hi(Lm(ng)+1,1,linew)+              &
     &                                hi(Lm(ng)  ,0,linew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        hi(0,Mm(ng)+1,linew)=0.5_r8*(hi(0,Mm(ng)  ,linew)+              &
     &                                hi(1,Mm(ng)+1,linew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        hi(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                             &
     &     (hi(Lm(ng)+1,Mm(ng)  ,linew)+                                &
     &      hi(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE hibc_tile
#endif

      END MODULE hibc_mod
