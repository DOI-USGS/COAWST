#include "cppdefs.h"
      MODULE ageicebc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the ice age.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC ageicebc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ageicebc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
      USE mod_scalars
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL ageicebc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    liold(ng), linew(ng), min_h(ng),              &
     &                    ICE(ng) % ui,                                 &
     &                    ICE(ng) % vi,                                 &
     &                    ICE(ng) % hi,                                 &
     &                    ICE(ng) % ageice,                             &
     &                    ICE(ng) % hage)
      RETURN
      END SUBROUTINE ageicebc

!
!***********************************************************************
      SUBROUTINE ageicebc_tile (ng, tile,                               &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew, min_h,                   &
     &                           ui, vi, hi,                            &
     &                           ageice, hage)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid

      implicit none

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: liold, linew
      real(r8), intent(in) :: min_h

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(in)    :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ageice(LBi:,LBj:,:)
      real(r8), intent(inout) :: hage(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ageice(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hage(LBi:UBi,LBj:UBj,2)
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

# ifdef WEST_AGEICECLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(1,j,linew).ge.0._r8) THEN
             hage(0,j,linew)=BOUNDARY(ng)%hi_west(j)*               &
     &                           BOUNDARY(ng)%ageice_west(j)
#  ifdef MASKING
             hage(0,j,linew)=hage(0,j,linew)*                   &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             hage(0,j,linew)=hage(0,j,linew)*                   &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
             hage(0,j,linew)=hi(1,j,liold)*ageice(1,j,liold)
#  ifdef MASKING
             hage(0,j,linew)=hage(0,j,linew)*                   &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             hage(0,j,linew)=hage(0,j,linew)*                   &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ENDIF
             ageice(0,j,linew) = hage(0,j,linew)/                       &
     &                       MAX(hi(0,j,linew),1.0E-6_r8)
             IF(hi(0,j,linew).LE.min_h) THEN
                hage(0,j,linew) = 0.0_r8
                ageice(0,j,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(0,j,linew)=ageice(0,j,linew)*                       &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             ageice(0,j,linew)=ageice(0,j,linew)*                       &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# elif defined WEST_AGEICEGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          hage(0,j,linew)=hi(1,j,linew)*ageice(1,j,linew)
#  ifdef MASKING
          hage(0,j,linew)=hage(0,j,linew)*                      &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hage(0,j,linew)=hage(0,j,linew)*                      &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
             ageice(0,j,linew) = hage(0,j,linew)/                       &
     &                       MAX(hi(0,j,linew),1.0E-6_r8)
             IF(hi(0,j,linew).LE.min_h) THEN
                hage(0,j,linew) = 0.0_r8
                ageice(0,j,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(0,j,linew)=ageice(0,j,linew)*                       &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             ageice(0,j,linew)=ageice(0,j,linew)*                       &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          hage(0,j,linew)=hi(1,j,linew)*ageice(1,j,linew)
          ageice(0,j,linew) = hage(0,j,linew)/                          &
     &         MAX(hi(0,j,linew),1.0E-6_r8)
#  ifdef MASKING
          hage(0,j,linew)=hage(0,j,linew)*                      &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hage(0,j,linew)=hage(0,j,linew)*                      &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
             ageice(0,j,linew) = hage(0,j,linew)/                       &
     &                       MAX(hi(0,j,linew),1.0E-6_r8)
             IF(hi(0,j,linew).LE.min_h) THEN
                hage(0,j,linew) = 0.0_r8
                ageice(0,j,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(0,j,linew)=ageice(0,j,linew)*                       &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             ageice(0,j,linew)=ageice(0,j,linew)*                       &
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

# ifdef EAST_AGEICECLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(Lm(ng)+1,j,linew).le.0._r8) THEN
             hage(Lm(ng)+1,j,linew)=BOUNDARY(ng)%hi_east(j)*        &
     &                                  BOUNDARY(ng)%ageice_east(j)
#  ifdef MASKING
             hage(Lm(ng)+1,j,linew)=hage(Lm(ng)+1,j,linew)*     &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             hage(Lm(ng)+1,j,linew)=hage(Lm(ng)+1,j,linew)*     &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
             hage(Lm(ng)+1,j,linew)=hi(Lm(ng),j,liold)*             &
     &                                  ageice(Lm(ng),j,liold)
#  ifdef MASKING
             hage(Lm(ng)+1,j,linew)=hage(Lm(ng)+1,j,linew)*     &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             hage(Lm(ng)+1,j,linew)=hage(Lm(ng)+1,j,linew)*     &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ENDIF
             ageice(Lm(ng)+1,j,linew) = hage(Lm(ng)+1,j,linew)/         &
     &                       MAX(hi(Lm(ng)+1,j,linew),1.0E-6_r8)
             IF(hi(Lm(ng)+1,j,linew).LE.min_h) THEN
                hage(Lm(ng)+1,j,linew) = 0.0_r8
                ageice(Lm(ng)+1,j,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(Lm(ng)+1,j,linew)=ageice(Lm(ng)+1,j,linew)*         &
     &                   GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             ageice(Lm(ng)+1,j,linew)=ageice(Lm(ng)+1,j,linew)*         &
     &                   GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_AGEICEGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          hage(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)*                &
     &                               ageice(Lm(ng),j,linew)
#  ifdef MASKING
          hage(Lm(ng)+1,j,linew)=hage(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          hage(Lm(ng)+1,j,linew)=hage(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
             ageice(Lm(ng)+1,j,linew) = hage(Lm(ng)+1,j,linew)/         &
     &                       MAX(hi(Lm(ng)+1,j,linew),1.0E-6_r8)
             IF(hi(Lm(ng)+1,j,linew).LE.min_h) THEN
                hage(Lm(ng)+1,j,linew) = 0.0_r8
                ageice(Lm(ng)+1,j,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(Lm(ng)+1,j,linew)=ageice(Lm(ng)+1,j,linew)*          &
     &                   GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             ageice(Lm(ng)+1,j,linew)=ageice(Lm(ng)+1,j,linew)*          &
     &                   GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          hage(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)*                &
     &                               ageice(Lm(ng),j,linew)
#  ifdef MASKING
          hage(Lm(ng)+1,j,linew)=hage(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          hage(Lm(ng)+1,j,linew)=hage(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
             ageice(Lm(ng)+1,j,linew) = hage(Lm(ng)+1,j,linew)/         &
     &                       MAX(hi(Lm(ng)+1,j,linew),1.0E-6_r8)
             IF(hi(Lm(ng)+1,j,linew).LE.min_h) THEN
                hage(Lm(ng)+1,j,linew) = 0.0_r8
                ageice(Lm(ng)+1,j,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(Lm(ng)+1,j,linew)=ageice(Lm(ng)+1,j,linew)*         &
     &                   GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             ageice(Lm(ng)+1,j,linew)=ageice(Lm(ng)+1,j,linew)*         &
     &                   GRID(ng)%rmask_wet(Lm(ng)+1,j)
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

# ifdef SOUTH_AGEICECLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,1,linew).ge.0._r8) THEN
             hage(i,0,linew)=BOUNDARY(ng)%hi_south(i)*              &
     &                           BOUNDARY(ng)%ageice_south(i)
#  ifdef MASKING
             hage(i,0,linew)=hage(i,0,linew)*                   &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             hage(i,0,linew)=hage(i,0,linew)*                   &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
             hage(i,0,linew)=hi(i,1,liold)*ageice(i,1,liold)
#  ifdef MASKING
             hage(i,0,linew)=hage(i,0,linew)*                   &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             hage(i,0,linew)=hage(i,0,linew)*                   &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ENDIF
             ageice(i,0,linew) = hage(i,0,linew)/                       &
     &                       MAX(hi(i,0,linew),1.0E-6_r8)
             IF(hi(i,0,linew).LE.min_h) THEN
                hage(i,0,linew) = 0.0_r8
                ageice(i,0,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(i,0,linew)=ageice(i,0,linew)*                       &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             ageice(i,0,linew)=ageice(i,0,linew)*                       &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# elif defined SOUTH_AGEICEGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          hage(i,0,linew)=hi(i,1,linew)*ageice(i,1,linew)
#  ifdef MASKING
          hage(i,0,linew)=hage(i,0,linew)*                      &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hage(i,0,linew)=hage(i,0,linew)*                      &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
             ageice(i,0,linew) = hage(i,0,linew)/                       &
     &                       MAX(hi(i,0,linew),1.0E-6_r8)
             IF(hi(i,0,linew).LE.min_h) THEN
                hage(i,0,linew) = 0.0_r8
                ageice(i,0,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(i,0,linew)=ageice(i,0,linew)*                               &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             ageice(i,0,linew)=ageice(i,0,linew)*                               &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          hage(i,0,linew)=hi(i,1,linew)*ageice(i,1,linew)
#  ifdef MASKING
          hage(i,0,linew)=hage(i,0,linew)*                      &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hage(i,0,linew)=hage(i,0,linew)*                      &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
             ageice(i,0,linew) = hage(i,0,linew)/                     &
     &                       MAX(hi(i,0,linew),1.0E-6_r8)
             IF(hi(i,0,linew).LE.min_h) THEN
                hage(i,0,linew) = 0.0_r8
                ageice(i,0,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(i,0,linew)=ageice(i,0,linew)*                     &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             ageice(i,0,linew)=ageice(i,0,linew)*                     &
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

# ifdef NORTH_AGEICECLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF(vi(i,Mm(ng)+1,linew).le.0._r8) THEN
             hage(i,Mm(ng)+1,linew)=BOUNDARY(ng)%hi_north(i)*       &
     &                                  BOUNDARY(ng)%ageice_north(i)
#  ifdef MASKING
             hage(i,Mm(ng)+1,linew)=hage(i,Mm(ng)+1,linew)*     &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             hage(i,Mm(ng)+1,linew)=hage(i,Mm(ng)+1,linew)*     &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ELSE
             hage(i,Mm(ng)+1,linew)=hi(i,Mm(ng),liold)*             &
     &                                  ageice(i,Mm(ng),liold)
#  ifdef MASKING
             hage(i,Mm(ng)+1,linew)=hage(i,Mm(ng)+1,linew)*     &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             hage(i,Mm(ng)+1,linew)=hage(i,Mm(ng)+1,linew)*     &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ENDIF
             ageice(i,Mm(ng)+1,linew) = hage(i,Mm(ng)+1,linew)/         &
     &                       MAX(hi(i,Mm(ng)+1,linew),1.0E-6_r8)
             IF(hi(i,Mm(ng)+1,linew).LE.min_h) THEN
                hage(i,Mm(ng)+1,linew) = 0.0_r8
                ageice(i,Mm(ng)+1,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(i,Mm(ng)+1,linew)=ageice(i,Mm(ng)+1,linew)*                 &
     &                   GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             ageice(i,Mm(ng)+1,linew)=ageice(i,Mm(ng)+1,linew)*                 &
     &                   GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_AGEICEGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          hage(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)*                &
     &                               ageice(i,Mm(ng),linew)
#  ifdef MASKING
          hage(i,Mm(ng)+1,linew)=hage(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          hage(i,Mm(ng)+1,linew)=hage(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
             ageice(i,Mm(ng)+1,linew) = hage(i,Mm(ng)+1,linew)/         &
     &                       MAX(hi(i,Mm(ng)+1,linew),1.0E-6_r8)
             IF(hi(i,Mm(ng)+1,linew).LE.min_h) THEN
                hage(i,Mm(ng)+1,linew) = 0.0_r8
                ageice(i,Mm(ng)+1,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(i,Mm(ng)+1,linew)=ageice(i,Mm(ng)+1,linew)*         &
     &                   GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             ageice(i,Mm(ng)+1,linew)=ageice(i,Mm(ng)+1,linew)*         &
     &                   GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif

        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          hage(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)*                &
     &                               ageice(i,Mm(ng),linew)
#  ifdef MASKING
          hage(i,Mm(ng)+1,linew)=hage(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          hage(i,Mm(ng)+1,linew)=hage(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
             ageice(i,Mm(ng)+1,linew) = hage(i,Mm(ng)+1,linew)/         &
     &                       MAX(hi(i,Mm(ng)+1,linew),1.0E-6_r8)
             IF(hi(i,Mm(ng)+1,linew).LE.min_h) THEN
                hage(i,Mm(ng)+1,linew) = 0.0_r8
                ageice(i,Mm(ng)+1,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             ageice(i,Mm(ng)+1,linew)=ageice(i,Mm(ng)+1,linew)*         &
     &                   GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             ageice(i,Mm(ng)+1,linew)=ageice(i,Mm(ng)+1,linew)*         &
     &                   GRID(ng)%rmask_wet(i,Mm(ng)+1)
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
        hage(0,0,linew)=0.5_r8*(hage(1,0,linew)+                &
     &                         hage(0,1,linew))
             ageice(0,0,linew) = hage(0,0,linew)/                       &
     &                       MAX(hi(0,0,linew),1.0E-6_r8)
             IF(hi(0,0,linew).LE.min_h) THEN
                hage(0,0,linew) = 0.0_r8
                ageice(0,0,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             hage(0,0,linew)=hage(0,0,linew)*                   &
     &                   GRID(ng)%rmask(0,0)
             ageice(0,0,linew)=ageice(0,0,linew)*                               &
     &                   GRID(ng)%rmask(0,0)
#  endif
#  ifdef WET_DRY
             hage(0,0,linew)=hage(0,0,linew)*                   &
     &                   GRID(ng)%rmask_wet(0,0)
             ageice(0,0,linew)=ageice(0,0,linew)*                               &
     &                   GRID(ng)%rmask_wet(0,0)
#  endif
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        hage(Lm(ng)+1,0,linew)=0.5_r8*(hage(Lm(ng)+1,1,linew)+  &
     &                                hage(Lm(ng)  ,0,linew))
             ageice(Lm(ng)+1,0,linew) = hage(Lm(ng)+1,0,linew)/         &
     &                       MAX(hi(Lm(ng)+1,0,linew),1.0E-6_r8)
             IF(hi(Lm(ng)+1,0,linew).LE.min_h) THEN
                hage(Lm(ng)+1,0,linew) = 0.0_r8
                ageice(Lm(ng)+1,0,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             hage(Lm(ng)+1,0,linew)=hage(Lm(ng)+1,0,linew)*     &
     &                   GRID(ng)%rmask(Lm(ng)+1,0)
             ageice(Lm(ng)+1,0,linew)=ageice(Lm(ng)+1,0,linew)*         &
     &                   GRID(ng)%rmask(Lm(ng)+1,0)
#  endif
#  ifdef WET_DRY
             hage(Lm(ng)+1,0,linew)=hage(Lm(ng)+1,0,linew)*     &
     &                   GRID(ng)%rmask_wet(Lm(ng)+1,0)
             ageice(Lm(ng)+1,0,linew)=ageice(Lm(ng)+1,0,linew)*         &
     &                   GRID(ng)%rmask_wet(Lm(ng)+1,0)
#  endif
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        hage(0,Mm(ng)+1,linew)=0.5_r8*(hage(0,Mm(ng)  ,linew)+  &
     &                                hage(1,Mm(ng)+1,linew))

             ageice(0,Mm(ng)+1,linew) = hage(0,Mm(ng)+1,linew)/         &
     &                       MAX(hi(0,Mm(ng)+1,linew),1.0E-6_r8)
             IF(hi(0,Mm(ng)+1,linew).LE.min_h) THEN
                hage(0,Mm(ng)+1,linew) = 0.0_r8
                ageice(0,Mm(ng)+1,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             hage(0,Mm(ng)+1,linew)=hage(0,Mm(ng)+1,linew)*     &
     &                   GRID(ng)%rmask(0,Mm(ng)+1)
             ageice(0,Mm(ng)+1,linew)=ageice(0,Mm(ng)+1,linew)* &
     &                   GRID(ng)%rmask(0,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             hage(0,Mm(ng)+1,linew)=hage(0,Mm(ng)+1,linew)*     &
     &                   GRID(ng)%rmask_wet(0,Mm(ng)+1)
             ageice(0,Mm(ng)+1,linew)=ageice(0,Mm(ng)+1,linew)* &
     &                   GRID(ng)%rmask_wet(0,Mm(ng)+1)
#  endif
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        hage(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                       &
     &             (hage(Lm(ng)+1,Mm(ng)  ,linew)+                  &
     &              hage(Lm(ng)  ,Mm(ng)+1,linew))
             ageice(Lm(ng)+1,Mm(ng)+1,linew) =                      &
     &         hage(Lm(ng)+1,Mm(ng)+1,linew)/                       &
     &                       MAX(hi(Lm(ng)+1,Mm(ng)+1,linew),1.0E-6_r8)
             IF(hi(Lm(ng)+1,Mm(ng)+1,linew).LE.min_h) THEN
                hage(Lm(ng)+1,Mm(ng)+1,linew) = 0.0_r8
                ageice(Lm(ng)+1,Mm(ng)+1,linew) = 0.0_r8
             END IF
#  ifdef MASKING
             hage(Lm(ng)+1,Mm(ng)+1,linew)=                         &
     &             hage(Lm(ng)+1,Mm(ng)+1,linew)*                   &
     &                   GRID(ng)%rmask(Lm(ng)+1,Mm(ng)+1)
      ageice(Lm(ng)+1,Mm(ng)+1,linew)=ageice(Lm(ng)+1,Mm(ng)+1,linew)*  &
     &                   GRID(ng)%rmask(Lm(ng)+1,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             hage(Lm(ng)+1,Mm(ng)+1,linew)=                         &
     &             hage(Lm(ng)+1,Mm(ng)+1,linew)*                   &
     &                   GRID(ng)%rmask_wet(Lm(ng)+1,Mm(ng)+1)
      ageice(Lm(ng)+1,Mm(ng)+1,linew)=ageice(Lm(ng)+1,Mm(ng)+1,linew)*  &
     &                   GRID(ng)%rmask_wet(Lm(ng)+1,Mm(ng)+1)
#  endif
      END IF
#endif
      RETURN
      END SUBROUTINE ageicebc_tile
#endif

      END MODULE ageicebc_mod
