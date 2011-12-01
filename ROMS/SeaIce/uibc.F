#include "cppdefs.h"
      MODULE uibc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the ice U-velocity.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC uibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE uibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL  uibc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 liuol(ng), liunw(ng),                            &
     &                 ICE(ng) % ui)
      RETURN
      END SUBROUTINE uibc
!
!***********************************************************************
      SUBROUTINE uibc_tile (ng, tile,                                   &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liuol, liunw,                          &
     &                           ui)
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
      integer, intent(in) :: liuol, liunw

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: ui(LBi:UBi,LBj:UBj,2)
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
        know=liuol
!
#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  ifdef WEST_MICLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ui(1,j,liunw)=BOUNDARY(ng)%ui_west(j)
#  ifdef MASKING
          ui(1,j,liunw)=ui(1,j,liunw)*                                  &
     &                   GRID(ng)%umask(1,j)
#  endif
#  ifdef WET_DRY
          ui(1,j,liunw)=ui(1,j,liunw)*                                  &
     &                   GRID(ng)%umask_wet(1,j)
#  endif
        END DO
# elif defined WEST_MIGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ui(1,j,liunw)=ui(2,j,liunw)
#  ifdef MASKING
          ui(1,j,liunw)=ui(1,j,liunw)*                                &
     &                   GRID(ng)%umask(1,j)
#  endif
#  ifdef WET_DRY
          ui(1,j,liunw)=ui(1,j,liunw)*                                &
     &                   GRID(ng)%umask_wet(1,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ui(1,j,liunw)=0.0_r8
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

# ifdef EAST_MICLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ui(Lm(ng)+1,j,liunw)=BOUNDARY(ng)%ui_east(j)
#  ifdef MASKING
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%umask_wet(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_MIGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng),j,liunw)
#  ifdef MASKING
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%umask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ui(Lm(ng)+1,j,liunw)=0.0_r8
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

# ifdef SOUTH_MICLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ui(i,0,liunw)=BOUNDARY(ng)%ui_south(i)
#  ifdef MASKING
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask(i,0)
#  endif
#  ifdef WET_DRY
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask_wet(i,0)
#  endif
        END DO
# elif defined SOUTH_MIGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ui(i,0,liunw)=ui(i,1,liunw)
#  ifdef MASKING
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask(i,0)
#  endif
#  ifdef WET_DRY
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ui(i,0,liunw)=gamma2(ng)*ui(i,1,liunw)
#  ifdef MASKING
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask(i,0)
#  endif
#  ifdef WET_DRY
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask_wet(i,0)
#  endif
        END DO
#  undef I_RANGE
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

# ifdef NORTH_MICLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ui(i,Mm(ng)+1,liunw)=BOUNDARY(ng)%ui_north(i)
#  ifdef MASKING
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask_wet(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_MIGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng),liunw)
#  ifdef MASKING
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ui(i,Mm(ng)+1,liunw)=gamma2(ng)*ui(i,Mm(ng),liunw)
#  ifdef MASKING
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask_wet(i,Mm(ng)+1)
#  endif
        END DO
#  undef I_RANGE
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
        ui(1,0,liunw)=0.5_r8*(ui(2,0,liunw)+                            &
     &                         ui(1,1,liunw))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        ui(Lm(ng)+1,0,liunw)=0.5_r8*(ui(Lm(ng)  ,0,liunw)+              &
     &                                ui(Lm(ng)+1,1,liunw))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        ui(1,Mm(ng)+1,liunw)=0.5_r8*(ui(2,Mm(ng)+1,liunw)+              &
     &                                ui(1,Mm(ng)  ,liunw))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        ui(Lm(ng)+1,Mm(ng)+1,liunw)=0.5_r8*                              &
     &                               (ui(Lm(ng)  ,Mm(ng)+1,liunw)+       &
     &                                ui(Lm(ng)+1,Mm(ng)  ,liunw))
      END IF
#endif
      RETURN
      END SUBROUTINE uibc_tile
#endif

      END MODULE uibc_mod
