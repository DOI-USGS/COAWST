#include "cppdefs.h"
      MODULE tibc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the internal ice
!  temperature.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC tibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tibc (ng, tile)
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
      CALL tibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold(ng), linew(ng), min_h(ng),                  &
     &                ICE(ng) % ui,                                     &
     &                ICE(ng) % vi,                                     &
     &                ICE(ng) % hi,                                     &
     &                ICE(ng) % ti,                                     &
     &                ICE(ng) % enthalpi)
      RETURN
      END SUBROUTINE tibc

!
!***********************************************************************
      SUBROUTINE tibc_tile (ng, tile,                                   &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew, min_h,                   &
     &                           ui, vi, hi, ti, enthalpi)
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
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: enthalpi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: enthalpi(LBi:UBi,LBj:UBj,2)
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

# ifdef WEST_TICLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF(ui(1,j,linew).ge.0._r8) THEN
             enthalpi(0,j,linew)=BOUNDARY(ng)%hi_west(j)* &
     &                           BOUNDARY(ng)%ti_west(j)
#  ifdef MASKING
             enthalpi(0,j,linew)=enthalpi(0,j,linew)*                   &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             enthalpi(0,j,linew)=enthalpi(0,j,linew)*                   &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
             enthalpi(0,j,linew)=enthalpi(1,j,liold)
#  ifdef MASKING
             enthalpi(0,j,linew)=enthalpi(0,j,linew)*                   &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             enthalpi(0,j,linew)=enthalpi(0,j,linew)*                   &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif

#  ifdef WET_DRY
             enthalpi(0,j,linew)=enthalpi(0,j,linew)*                   &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
         ENDIF
          ti(0,j,linew) = enthalpi(0,j,linew)/                          &
     &                       MAX(hi(0,j,linew),1.0E-6_r8)
          IF (hi(0,j,linew).LE.min_h) THEN
                enthalpi(0,j,linew) = 0.0_r8
                ti(0,j,linew) = 0.0_r8
          END IF
        END DO
# elif defined WEST_TIGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          enthalpi(0,j,linew)=hi(1,j,linew)*ti(1,j,linew)
#  ifdef MASKING
          enthalpi(0,j,linew)=enthalpi(0,j,linew)*                      &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          enthalpi(0,j,linew)=enthalpi(0,j,linew)*                      &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ti(0,j,linew) = enthalpi(0,j,linew)/                          &
     &                       MAX(hi(0,j,linew),1.0E-6_r8)
          IF (hi(0,j,linew).LE.min_h) THEN
            enthalpi(0,j,linew) = 0.0_r8
            ti(0,j,linew) = 0.0_r8
          END IF
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          enthalpi(0,j,linew)=hi(1,j,linew)*ti(1,j,linew)
          ti(0,j,linew) = enthalpi(0,j,linew)/                          &
     &         MAX(hi(0,j,linew),1.0E-6_r8)
#  ifdef MASKING
          enthalpi(0,j,linew)=enthalpi(0,j,linew)*                      &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          enthalpi(0,j,linew)=enthalpi(0,j,linew)*                      &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
             ti(0,j,linew) = enthalpi(0,j,linew)/                       &
     &                       MAX(hi(0,j,linew),1.0E-6_r8)
             IF(hi(0,j,linew).LE.min_h) THEN
                enthalpi(0,j,linew) = 0.0_r8
                ti(0,j,linew) = 0.0_r8
             END IF
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

# ifdef EAST_TICLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          IF (ui(Lm(ng)+1,j,linew).le.0._r8) THEN
            enthalpi(Lm(ng)+1,j,linew)=BOUNDARY(ng)%hi_east(j)*         &
     &                                  BOUNDARY(ng)%ti_east(j)
#  ifdef MASKING
            enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*      &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
            enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*      &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
            enthalpi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,liold)*              &
     &                                  ti(Lm(ng),j,liold)
#  ifdef MASKING
            enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*      &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
            enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*      &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
            ti(Lm(ng)+1,j,linew) = enthalpi(Lm(ng)+1,j,linew)/          &
     &                       MAX(hi(Lm(ng)+1,j,linew),1.0E-6_r8)
            IF (hi(Lm(ng)+1,j,linew).LE.min_h) THEN
              enthalpi(Lm(ng)+1,j,linew) = 0.0_r8
              ti(Lm(ng)+1,j,linew) = 0.0_r8
            END IF
          ENDIF
        END DO
# elif defined EAST_TIGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          enthalpi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)*                &
     &                               ti(Lm(ng),j,linew)
#  ifdef MASKING
          enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ti(Lm(ng)+1,j,linew) = enthalpi(Lm(ng)+1,j,linew)/            &
     &                       MAX(hi(Lm(ng)+1,j,linew),1.0E-6_r8)
          IF (hi(Lm(ng)+1,j,linew).LE.min_h) THEN
            enthalpi(Lm(ng)+1,j,linew) = 0.0_r8
            ti(Lm(ng)+1,j,linew) = 0.0_r8
          END IF
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          enthalpi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)*                &
     &                               ti(Lm(ng),j,linew)
#  ifdef MASKING
          enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          enthalpi(Lm(ng)+1,j,linew)=enthalpi(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ti(Lm(ng)+1,j,linew) = enthalpi(Lm(ng)+1,j,linew)/            &
     &                       MAX(hi(Lm(ng)+1,j,linew),1.0E-6_r8)
          IF (hi(Lm(ng)+1,j,linew).LE.min_h) THEN
            enthalpi(Lm(ng)+1,j,linew) = 0.0_r8
            ti(Lm(ng)+1,j,linew) = 0.0_r8
          END IF
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

# ifdef SOUTH_TICLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF (vi(i,1,linew).ge.0._r8) THEN
            enthalpi(i,0,linew)=BOUNDARY(ng)%hi_south(i)*               &
     &                          BOUNDARY(ng)%ti_south(i)            
#  ifdef MASKING
            enthalpi(i,0,linew)=enthalpi(i,0,linew)*                    &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
            enthalpi(i,0,linew)=enthalpi(i,0,linew)*                    &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
            enthalpi(i,0,linew)=enthalpi(i,1,liold)
#  ifdef MASKING
            enthalpi(i,0,linew)=enthalpi(i,0,linew)*                    &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
            enthalpi(i,0,linew)=enthalpi(i,0,linew)*                    &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
            ti(i,0,linew) = enthalpi(i,0,linew)/                        &
     &                       MAX(hi(i,0,linew),1.0E-6_r8)
            IF (hi(i,0,linew).LE.min_h) THEN
              enthalpi(i,0,linew) = 0.0_r8
              ti(i,0,linew) = 0.0_r8
            END IF
          ENDIF
        END DO
# elif defined SOUTH_TIGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          enthalpi(i,0,linew)=hi(i,1,linew)*ti(i,1,linew)
#  ifdef MASKING
          enthalpi(i,0,linew)=enthalpi(i,0,linew)*                      &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          enthalpi(i,0,linew)=enthalpi(i,0,linew)*                      &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ti(i,0,linew) = enthalpi(i,0,linew)/                          &
     &                       MAX(hi(i,0,linew),1.0E-6_r8)
          IF (hi(i,0,linew).LE.min_h) THEN
            enthalpi(i,0,linew) = 0.0_r8
            ti(i,0,linew) = 0.0_r8
          END IF
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          enthalpi(i,0,linew)=enthalpi(i,1,linew)
#  ifdef MASKING
          enthalpi(i,0,linew)=enthalpi(i,0,linew)*                      &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          enthalpi(i,0,linew)=enthalpi(i,0,linew)*                      &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ti(i,0,linew) = enthalpi(i,0,linew)/                          &
     &                    MAX(hi(i,0,linew),1.0E-6_r8)
          IF (hi(i,0,linew).LE.min_h) THEN
            enthalpi(i,0,linew) = 0.0_r8
            ti(i,0,linew) = 0.0_r8
          END IF
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

# ifdef NORTH_TICLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF (vi(i,Mm(ng)+1,linew).le.0._r8) THEN
            enthalpi(i,Mm(ng)+1,linew)=BOUNDARY(ng)%hi_north(i)*        &
     &                                 BOUNDARY(ng)%ti_north(i)
#  ifdef MASKING
            enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*      &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
            enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*      &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ELSE
             enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng),liold)
#  ifdef MASKING
             enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*     &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*     &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ENDIF
          ti(i,Mm(ng)+1,linew) = enthalpi(i,Mm(ng)+1,linew)/            &
     &                       MAX(hi(i,Mm(ng)+1,linew),1.0E-6_r8)
          IF (hi(i,Mm(ng)+1,linew).LE.min_h) THEN
            enthalpi(i,Mm(ng)+1,linew) = 0.0_r8
            ti(i,Mm(ng)+1,linew) = 0.0_r8
          END IF
        END DO
# elif defined NORTH_TIGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          enthalpi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)*                &
     &                               ti(i,Mm(ng),linew)
#  ifdef MASKING
          enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ti(i,Mm(ng)+1,linew) = enthalpi(i,Mm(ng)+1,linew)/            &
     &                       MAX(hi(i,Mm(ng)+1,linew),1.0E-6_r8)
          IF (hi(i,Mm(ng)+1,linew).LE.min_h) THEN
            enthalpi(i,Mm(ng)+1,linew) = 0.0_r8
            ti(i,Mm(ng)+1,linew) = 0.0_r8
          END IF
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          enthalpi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)*                &
     &                               ti(i,Mm(ng),linew)
#  ifdef MASKING
          enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          enthalpi(i,Mm(ng)+1,linew)=enthalpi(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ti(i,Mm(ng)+1,linew) = enthalpi(i,Mm(ng)+1,linew)/            &
     &                       MAX(hi(i,Mm(ng)+1,linew),1.0E-6_r8)
          IF (hi(i,Mm(ng)+1,linew).LE.min_h) THEN
            enthalpi(i,Mm(ng)+1,linew) = 0.0_r8
            ti(i,Mm(ng)+1,linew) = 0.0_r8
          END IF
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
        enthalpi(0,0,linew)=0.5_r8*(enthalpi(1,0,linew)+                &
     &                         enthalpi(0,1,linew))
#  ifdef MASKING
             enthalpi(0,0,linew)=enthalpi(0,0,linew)* &
     &                   GRID(ng)%rmask(0,0)
#  endif
#  ifdef WET_DRY
             enthalpi(0,0,linew)=enthalpi(0,0,linew)* &
     &                   GRID(ng)%rmask_wet(0,0)
#  endif
        ti(0,0,linew) = enthalpi(0,0,linew)/ &
     &                       MAX(hi(0,0,linew),1.0E-6_r8)
        IF (hi(0,0,linew).LE.min_h) THEN
          enthalpi(0,0,linew) = 0.0_r8
          ti(0,0,linew) = 0.0_r8
        END IF
      END IF
      IF (SOUTHERN_EDGE.and.EASTERN_EDGE) THEN
        enthalpi(Lm(ng)+1,0,linew)=0.5_r8*(enthalpi(Lm(ng)+1,1,linew)+  &
     &                                enthalpi(Lm(ng)  ,0,linew))
#  ifdef MASKING
        enthalpi(Lm(ng)+1,0,linew)=enthalpi(Lm(ng)+1,0,linew)*          &
     &                   GRID(ng)%rmask(Lm(ng)+1,0)
#  endif
#  ifdef WET_DRY
        enthalpi(Lm(ng)+1,0,linew)=enthalpi(Lm(ng)+1,0,linew)*          &
     &                   GRID(ng)%rmask_wet(Lm(ng)+1,0)
#  endif
        ti(Lm(ng)+1,0,linew) = enthalpi(Lm(ng)+1,0,linew)/              &
     &                       MAX(hi(Lm(ng)+1,0,linew),1.0E-6_r8)
        IF (hi(Lm(ng)+1,0,linew).LE.min_h) THEN
          enthalpi(Lm(ng)+1,0,linew) = 0.0_r8
          ti(Lm(ng)+1,0,linew) = 0.0_r8
        END IF
      END IF
      IF (NORTHERN_EDGE.and.WESTERN_EDGE) THEN
        enthalpi(0,Mm(ng)+1,linew)=0.5_r8*(enthalpi(0,Mm(ng)  ,linew)+  &
     &                                enthalpi(1,Mm(ng)+1,linew))
#  ifdef MASKING
        enthalpi(0,Mm(ng)+1,linew)=enthalpi(0,Mm(ng)+1,linew)*          &
     &                   GRID(ng)%rmask(0,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
        enthalpi(0,Mm(ng)+1,linew)=enthalpi(0,Mm(ng)+1,linew)*          &
     &                   GRID(ng)%rmask_wet(0,Mm(ng)+1)
#  endif
        ti(0,Mm(ng)+1,linew) = enthalpi(0,Mm(ng)+1,linew)/              &
     &                       MAX(hi(0,Mm(ng)+1,linew),1.0E-6_r8)
        IF (hi(0,Mm(ng)+1,linew).LE.min_h) THEN
          enthalpi(0,Mm(ng)+1,linew) = 0.0_r8
          ti(0,Mm(ng)+1,linew) = 0.0_r8
        END IF
      END IF
      IF (NORTHERN_EDGE.and.EASTERN_EDGE) THEN
        enthalpi(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                       &
     &             (enthalpi(Lm(ng)+1,Mm(ng)  ,linew)+                  &
     &              enthalpi(Lm(ng)  ,Mm(ng)+1,linew))
#  ifdef MASKING
        enthalpi(Lm(ng)+1,Mm(ng)+1,linew)=                              &
     &             enthalpi(Lm(ng)+1,Mm(ng)+1,linew)*                   &
     &             GRID(ng)%rmask(Lm(ng)+1,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
        enthalpi(Lm(ng)+1,Mm(ng)+1,linew)=                              &
     &             enthalpi(Lm(ng)+1,Mm(ng)+1,linew)*                   &
     &             GRID(ng)%rmask_wet(Lm(ng)+1,Mm(ng)+1)
#  endif
        ti(Lm(ng)+1,Mm(ng)+1,linew) =                                   &
     &         enthalpi(Lm(ng)+1,Mm(ng)+1,linew)/                       &
     &         MAX(hi(Lm(ng)+1,Mm(ng)+1,linew),1.0E-6_r8)
        IF (hi(Lm(ng)+1,Mm(ng)+1,linew).LE.min_h) THEN
          enthalpi(Lm(ng)+1,Mm(ng)+1,linew) = 0.0_r8
          ti(Lm(ng)+1,Mm(ng)+1,linew) = 0.0_r8
        END IF
      END IF
#endif
      RETURN
      END SUBROUTINE tibc_tile
#endif

      END MODULE tibc_mod
