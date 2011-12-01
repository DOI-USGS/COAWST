#include "cppdefs.h"
      MODULE vibc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the ice V-velocity.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC vibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL vibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                liuol(ng), liunw(ng),                             &
     &                ICE(ng) % vi)
      RETURN
      END SUBROUTINE vibc
!
!***********************************************************************
      SUBROUTINE vibc_tile (ng, tile,                                   &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liuol, liunw,                          &
     &                           vi)
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
      real(r8), intent(inout) :: vi(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: vi(LBi:UBi,LBj:UBj,2)
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

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

# ifdef WEST_MICLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=JstrV,Jend
          vi(0,j,liunw)=BOUNDARY(ng)%vi_west(j)
#  ifdef MASKING
          vi(0,j,liunw)=vi(0,j,liunw)*                                  &
     &                   GRID(ng)%vmask(0,j)
#  endif
#  ifdef WET_DRY
          vi(0,j,liunw)=vi(0,j,liunw)*                                  &
     &                   GRID(ng)%vmask_wet(0,j)
#  endif
        END DO
# elif defined WEST_MIGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=JstrV,Jend
          vi(0,j,liunw)=vi(1,j,liunw)
#  ifdef MASKING
          vi(0,j,liunw)=vi(0,j,liunw)*                                  &
     &                   GRID(ng)%vmask(0,j)
#  endif
#  ifdef WET_DRY
          vi(0,j,liunw)=vi(0,j,liunw)*                                  &
     &                   GRID(ng)%vmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
        DO j=J_RANGE
          vi(0,j,liunw)=gamma2(ng)*vi(0,j,liunw)
#  ifdef MASKING
          vi(0,j,liunw)=vi(0,j,liunw)*                                  &
     &                   GRID(ng)%vmask(0,j)
#  endif
#  ifdef WET_DRY
          vi(0,j,liunw)=vi(0,j,liunw)*                                  &
     &                   GRID(ng)%vmask_wet(0,j)
#  endif
        END DO
#  undef J_RANGE
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
        DO j=JstrV,Jend
          vi(Lm(ng)+1,j,liunw)=BOUNDARY(ng)%vi_east(j)
#  ifdef MASKING
          vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%vmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%vmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_MIGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=JstrV,Jend
          vi(Lm(ng)+1,j,liunw)=vi(Lm(ng),j,liunw)
#  ifdef MASKING
          vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                    &
     &                       GRID(ng)%vmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                    &
     &                       GRID(ng)%vmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend
#  else
#   define J_RANGE Jstr,JendR
#  endif
        DO j=J_RANGE
          vi(Lm(ng)+1,j,liunw)=gamma2(ng)*vi(Lm(ng),j,liunw)
#  ifdef MASKING
          vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%vmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          vi(Lm(ng)+1,j,liunw)=vi(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%vmask_wet(Lm(ng)+1,j)
#  endif
        END DO
#  undef J_RANGE
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
        DO i=Istr,Iend
          vi(i,1,liunw)=BOUNDARY(ng)%vi_south(i)
#  ifdef MASKING
          vi(i,1,liunw)=vi(i,1,liunw)*                                  &
     &                   GRID(ng)%vmask(i,1)
#  endif
#  ifdef WET_DRY
          vi(i,1,liunw)=vi(i,1,liunw)*                                  &
     &                   GRID(ng)%vmask_wet(i,1)
#  endif
        END DO
# elif defined SOUTH_MIGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          vi(i,1,liunw)=vi(i,2,liunw)
#  ifdef MASKING
          vi(i,1,liunw)=vi(i,1,liunw)*                                  &
     &                   GRID(ng)%vmask(i,1)
#  endif
#  ifdef WET_DRY
          vi(i,1,liunw)=vi(i,1,liunw)*                                  &
     &                   GRID(ng)%vmask_wet(i,1)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          vi(i,1,liunw)=0.0_r8
        END DO
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
        DO i=Istr,Iend
          vi(i,Mm(ng)+1,liunw)=BOUNDARY(ng)%vi_north(i)
#  ifdef MASKING
          vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%vmask(i,Mm(ng))
#  endif
#  ifdef WET_DRY
          vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%vmask_wet(i,Mm(ng))
#  endif
        END DO
# elif defined NORTH_MIGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng),liunw)
#  ifdef MASKING
          vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%vmask(i,Mm(ng))
#  endif
#  ifdef WET_DRY
          vi(i,Mm(ng)+1,liunw)=vi(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%vmask_wet(i,Mm(ng))
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          vi(i,Mm(ng)+1,liunw)=0.0_r8
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
        vi(0,1,liunw)=0.5_r8*(vi(0,2,liunw)+                            &
     &                         vi(1,1,liunw))
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        vi(Lm(ng)+1,1,liunw)=0.5_r8*(vi(Lm(ng)  ,1,liunw)+              &
     &                                vi(Lm(ng)+1,2,liunw))
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        vi(0,Mm(ng)+1,liunw)=0.5_r8*(vi(0,Mm(ng)  ,liunw)+              &
     &                                vi(1,Mm(ng)+1,liunw))
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        vi(Lm(ng)+1,Mm(ng)+1,liunw)=0.5_r8*                             &
     &                               (vi(Lm(ng)+1,Mm(ng)  ,liunw)+      &
     &                                vi(Lm(ng)  ,Mm(ng)+1,liunw))
      END IF
#endif
      RETURN
      END SUBROUTINE vibc_tile
#endif

      END MODULE vibc_mod
