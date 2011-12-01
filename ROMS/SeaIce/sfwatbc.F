#include "cppdefs.h"
      MODULE sfwatbc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the surface melt water.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC sfwatbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sfwatbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL sfwatbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   liold(ng), linew(ng),                          &
     &                   ICE(ng) % ui,                                  &
     &                   ICE(ng) % vi,                                  &
     &                   ICE(ng) % sfwat)
      RETURN
      END SUBROUTINE sfwatbc

!
!***********************************************************************
      SUBROUTINE sfwatbc_tile (ng, tile,                                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew,                          &
     &                           ui, vi, sfwat)
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
      real(r8), intent(inout) :: sfwat(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: sfwat(LBi:UBi,LBj:UBj,2)
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

# ifdef WEST_SFWATCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(1,j,linew).ge.0._r8) THEN
             sfwat(0,j,linew)=BOUNDARY(ng)%sfwat_west(j)
#  ifdef MASKING
             sfwat(0,j,linew)=sfwat(0,j,linew)*                         &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             sfwat(0,j,linew)=sfwat(0,j,linew)*                         &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
             sfwat(0,j,linew)=sfwat(1,j,liold)
#  ifdef MASKING
             sfwat(0,j,linew)=sfwat(0,j,linew)*                         &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             sfwat(0,j,linew)=sfwat(0,j,linew)*                         &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ENDIF
        END DO
# elif defined WEST_SFWATGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          sfwat(0,j,linew)=sfwat(1,j,linew)
#  ifdef MASKING
          sfwat(0,j,linew)=sfwat(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sfwat(0,j,linew)=sfwat(0,j,linew)*                            &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          sfwat(0,j,linew)=sfwat(1,j,linew)
#  ifdef MASKING
          sfwat(0,j,linew)=sfwat(0,j,linew)*                            &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          sfwat(0,j,linew)=sfwat(0,j,linew)*                            &
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

# ifdef EAST_SFWATCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(Lm(ng)+1,j,linew).ge.0._r8) THEN
             sfwat(Lm(ng)+1,j,linew)=BOUNDARY(ng)%sfwat_east(j)
#  ifdef MASKING
             sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng)+1,j,linew)*           &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng)+1,j,linew)*           &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
             sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng),j,liold)
#  ifdef MASKING
             sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng)+1,j,linew)*           &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng)+1,j,linew)*           &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ENDIF
        END DO
# elif defined EAST_SFWATGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng),j,linew)
#  ifdef MASKING
          sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng),j,linew)
#  ifdef MASKING
          sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng)+1,j,linew)*              &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          sfwat(Lm(ng)+1,j,linew)=sfwat(Lm(ng)+1,j,linew)*              &
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

# ifdef SOUTH_SFWATCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,1,linew).ge.0._r8) THEN
             sfwat(i,0,linew)=BOUNDARY(ng)%sfwat_south(i)
#  ifdef MASKING
             sfwat(i,0,linew)=sfwat(i,0,linew)*                         &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             sfwat(i,0,linew)=sfwat(i,0,linew)*                         &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
             sfwat(i,0,linew)=sfwat(i,1,liold)
#  ifdef MASKING
             sfwat(i,0,linew)=sfwat(i,0,linew)*                         &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             sfwat(i,0,linew)=sfwat(i,0,linew)*                         &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ENDIF
        END DO
# elif defined SOUTH_SFWATGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          sfwat(i,0,linew)=sfwat(i,1,linew)
#  ifdef MASKING
          sfwat(i,0,linew)=sfwat(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sfwat(i,0,linew)=sfwat(i,0,linew)*                            &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          sfwat(i,0,linew)=sfwat(i,1,linew)
#  ifdef MASKING
          sfwat(i,0,linew)=sfwat(i,0,linew)*                            &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          sfwat(i,0,linew)=sfwat(i,0,linew)*                            &
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

# ifdef NORTH_SFWATCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,Mm(ng)+1,linew).le.0._r8) THEN
             sfwat(i,Mm(ng)+1,linew)=BOUNDARY(ng)%sfwat_north(i)
#  ifdef MASKING
             sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng)+1,linew)*           &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng)+1,linew)*           &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ELSE
             sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng),liold)
#  ifdef MASKING
             sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng)+1,linew)*           &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng)+1,linew)*           &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ENDIF
        END DO
# elif defined NORTH_SFWATGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng),linew)
#  ifdef MASKING
          sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng),linew)
#  ifdef MASKING
          sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng)+1,linew)*              &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          sfwat(i,Mm(ng)+1,linew)=sfwat(i,Mm(ng)+1,linew)*              &
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
        sfwat(0,0,linew)=0.5_r8*(sfwat(1,0,linew)+                      &
     &                         sfwat(0,1,linew))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        sfwat(Lm(ng)+1,0,linew)=0.5_r8*(sfwat(Lm(ng)+1,1,linew)+        &
     &                                sfwat(Lm(ng)  ,0,linew))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        sfwat(0,Mm(ng)+1,linew)=0.5_r8*(sfwat(0,Mm(ng)  ,linew)+        &
     &                                sfwat(1,Mm(ng)+1,linew))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        sfwat(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                          &
     &             (sfwat(Lm(ng)+1,Mm(ng)  ,linew)+                     &
     &              sfwat(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE sfwatbc_tile
#endif

      END MODULE sfwatbc_mod
