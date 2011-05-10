#include "cppdefs.h"
      MODULE vs2dbc_mod
#ifdef WEC
!
!svn $Id: vs2dbc_im.F 779 2008-10-04 23:15:47Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated V-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: vs2dbc, vs2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vs2dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL vs2dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  OCEAN(ng) % vbar_stokes)
      RETURN
      END SUBROUTINE vs2dbc
!
!***********************************************************************
      SUBROUTINE vs2dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       vbar_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_forces
      USE mod_grid
      USE mod_scalars
# ifdef WET_DRY
      USE mod_grid
# endif
# ifdef REFINED_GRID
      USE mod_refined
      USE mod_stepping
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: vbar_stokes(LBi:,LBj:)
# else
      real(r8), intent(inout) :: vbar_stokes(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff5
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8):: cff, cff1, cff2, dVde, dVdt, dVdx, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

# ifdef REFINED_GRID
      integer :: t1, t2, ngp
      real(r8) :: fac, fac1, fac2, my_val
# endif

#include "set_bounds.h"

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      IF (ng.eq.1) THEN
#endif

# if !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_M2SRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,Jstr  )=vbar_stokes(i  ,Jstr  )-                       &
     &                   vbar_stokes(i-1,Jstr  )
          grad(i,Jstr+1)=vbar_stokes(i  ,Jstr+1)-                       &
     &                   vbar_stokes(i-1,Jstr+1)
        END DO
        DO i=Istr,Iend
          dVdt=vbar_stokes(i,Jstr+1)-vbar_stokes(i,Jstr+1)
          dVde=vbar_stokes(i,Jstr+1)-vbar_stokes(i,Jstr+2)
          IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(i,Jstr+1)+grad(i+1,Jstr+1))).gt.0.0_r8) THEN
            dVdx=grad(i  ,Jstr+1)
          ELSE
            dVdx=grad(i+1,Jstr+1)
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
#   ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#   else
          Cx=0.0_r8
#   endif
          Ce=dVdt*dVde
          vbar_stokes(i,Jstr)=(cff*vbar_stokes(i,Jstr  )+               &
     &                        Ce *vbar_stokes(i,Jstr+1)-                &
     &                        MAX(Cx,0.0_r8)*grad(i  ,Jstr)-            &
     &                        MIN(Cx,0.0_r8)*grad(i+1,Jstr))/           &
     &                        (cff+Ce)
#   ifdef MASKING
          vbar_stokes(i,Jstr)=vbar_stokes(i,Jstr)*                      &
     &                        GRID(ng)%vmask(i,Jstr)
#   endif
#   ifdef WET_DRY
          vbar_stokes(i,Jstr)=vbar_stokes(i,Jstr)*                      &
     &                        GRID(ng)%vmask_wet(i,Jstr)
#   endif
        END DO

#  elif defined SOUTH_M2SCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          vbar_stokes(i,Jstr)=BOUNDARY(ng)%vbarstokes_south(i)
#   ifdef MASKING
          vbar_stokes(i,Jstr)=vbar_stokes(i,Jstr)*                      &
     &                        GRID(ng)%vmask(i,Jstr)
#   endif
#   ifdef WET_DRY
          vbar_stokes(i,Jstr)=vbar_stokes(i,Jstr)*                      &
     &                        GRID(ng)%vmask_wet(i,Jstr)
#   endif
        END DO

#  else
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          vbar_stokes(i,Jstr)=vbar_stokes(i,Jstr+1)
#   ifdef MASKING
          vbar_stokes(i,Jstr)=vbar_stokes(i,Jstr)*                      &
     &                        GRID(ng)%vmask(i,Jstr)
#   endif
#   ifdef WET_DRY
          vbar_stokes(i,Jstr)=vbar_stokes(i,Jstr)*                      &
     &                        GRID(ng)%vmask_wet(i,Jstr)
#   endif
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

#  if defined NORTH_M2SRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,Jend  )=vbar_stokes(i  ,Jend  )-                       &
     &                   vbar_stokes(i-1,Jend  )
          grad(i,Jend+1)=vbar_stokes(i  ,Jend+1)-                       &
     &                   vbar_stokes(i-1,Jend+1)
        END DO
        DO i=Istr,Iend
          dVdt=vbar_stokes(i,Jend)-vbar_stokes(i,Jend  )
          dVde=vbar_stokes(i,Jend)-vbar_stokes(i,Jend-1)
          IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
            dVdx=grad(i  ,Jend)
          ELSE
            dVdx=grad(i+1,Jend)
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
#   ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
#   else
          Cx=0.0_r8
#   endif
          Ce=dVdt*dVde
          vbar_stokes(i,Jend+1)=(cff*vbar_stokes(i,Jend+1)+             &
     &                          Ce *vbar_stokes(i,Jend  )-              &
     &                          MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/       &
     &                          (cff+Ce)
#   ifdef MASKING
          vbar_stokes(i,Jend+1)=vbar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%vmask(i,Jend+1)
#   endif
#   ifdef WET_DRY
          vbar_stokes(i,Jend+1)=vbar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%vmask_wet(i,Jend+1)
#   endif
        END DO

#  elif defined NORTH_M2SCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          vbar_stokes(i,Jend+1)=BOUNDARY(ng)%vbarstokes_north(i)
#   ifdef MASKING
          vbar_stokes(i,Jend+1)=vbar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%vmask(i,Jend+1)
#   endif
#   ifdef WET_DRY
          vbar_stokes(i,Jend+1)=vbar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%vmask_wet(i,Jend+1)
#   endif
        END DO

#  else
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          vbar_stokes(i,Jend+1)=vbar_stokes(i,Jend)
#   ifdef MASKING
          vbar_stokes(i,Jend+1)=vbar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%vmask(i,Jend+1)
#   endif
#   ifdef WET_DRY
          vbar_stokes(i,Jend+1)=vbar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%vmask_wet(i,Jend+1)
#   endif
        END DO
#  endif
      END IF
# endif

# if !defined EW_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_M2SRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO j=JstrV-1,Jend
          grad(Istr-1,j)=vbar_stokes(Istr-1,j+1)-                       &
     &                   vbar_stokes(Istr-1,j  )
          grad(Istr  ,j)=vbar_stokes(Istr  ,j+1)-                       &
     &                   vbar_stokes(Istr  ,j  )
        END DO
        DO j=JstrV,Jend
          dVdt=vbar_stokes(Istr,j)-vbar_stokes(Istr  ,j)
          dVdx=vbar_stokes(Istr,j)-vbar_stokes(Istr+1,j)
          IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(Istr,j-1)+grad(Istr,j))).gt.0.0_r8) THEN
            dVde=grad(Istr,j-1)
          ELSE
            dVde=grad(Istr,j  )
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
          Cx=dVdt*dVdx
#   ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#   else
          Ce=0.0_r8
#   endif
          vbar_stokes(Istr-1,j)=(cff*vbar_stokes(Istr-1,j)+             &
     &                          Cx *vbar_stokes(Istr  ,j)-              &
     &                          MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-        &
     &                          MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/       &
     &                          (cff+Cx)
#   ifdef MASKING
          vbar_stokes(Istr-1,j)=vbar_stokes(Istr-1,j)*                  &
     &                          GRID(ng)%vmask(Istr-1,j)
#   endif
#   ifdef WET_DRY
          vbar_stokes(Istr-1,j)=vbar_stokes(Istr-1,j)*                  &
     &                          GRID(ng)%vmask_wet(Istr-1,j)
#   endif
        END DO

#  elif defined WEST_M2SCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=JstrV,Jend
          vbar_stokes(Istr-1,j)=BOUNDARY(ng)%vbarstokes_west(j)
#   ifdef MASKING
          vbar_stokes(Istr-1,j)=vbar_stokes(Istr-1,j)*                  &
     &                          GRID(ng)%vmask(Istr-1,j)
#   endif
#   ifdef WET_DRY
          vbar_stokes(Istr-1,j)=vbar_stokes(Istr-1,j)*                  &
     &                          GRID(ng)%vmask_wet(Istr-1,j)
#   endif
        END DO

#  else
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO j=J_RANGE
          vbar_stokes(Istr-1,j)=gamma2(ng)*vbar_stokes(Istr,j)
#   ifdef MASKING
          vbar_stokes(Istr-1,j)=vbar_stokes(Istr-1,j)*                  &
     &                          GRID(ng)%vmask(Istr-1,j)
#   endif
#   ifdef WET_DRY
          vbar_stokes(Istr-1,j)=vbar_stokes(Istr-1,j)*                  &
     &                          GRID(ng)%vmask_wet(Istr-1,j)
#   endif
        END DO
#   undef J_RANGE
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

#  if defined EAST_M2SRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO j=JstrV-1,Jend
          grad(Iend  ,j)=vbar_stokes(Iend  ,j+1)-                       &
     &                   vbar_stokes(Iend  ,j  )
          grad(Iend+1,j)=vbar_stokes(Iend+1,j+1)-                       &
     &                   vbar_stokes(Iend+1,j  )
        END DO
        DO j=JstrV,Jend
          dVdt=vbar_stokes(Iend,j)-vbar_stokes(Iend  ,j)
          dVdx=vbar_stokes(Iend,j)-vbar_stokes(Iend-1,j)
          IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
          IF ((dVdt*(grad(Iend,j-1)+grad(Iend,j))).gt.0.0_r8) THEN
            dVde=grad(Iend,j-1)
          ELSE
            dVde=grad(Iend,j  )
          END IF
          cff=MAX(dVdx*dVdx+dVde*dVde,eps)
          Cx=dVdt*dVdx
#   ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dVdt*dVde,-cff))
#   else
          Ce=0.0_r8
#   endif
          vbar_stokes(Iend+1,j)=(cff*vbar_stokes(Iend+1,j)+             &
     &                          Cx *vbar_stokes(Iend  ,j)-              &
     &                          MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-        &
     &                          MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/       &
     &                          (cff+Cx)
#   ifdef MASKING
          vbar_stokes(Iend+1,j)=vbar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%vmask(Iend+1,j)
#   endif
#   ifdef WET_DRY
          vbar_stokes(Iend+1,j)=vbar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%vmask_wet(Iend+1,j)
#   endif
        END DO

#  elif defined EAST_M2SCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=JstrV,Jend
          vbar_stokes(Iend+1,j)=BOUNDARY(ng)%vbarstokes_east(j)
#   ifdef MASKING
          vbar_stokes(Iend+1,j)=vbar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%vmask(Iend+1,j)
#   endif
#   ifdef WET_DRY
          vbar_stokes(Iend+1,j)=vbar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%vmask_wet(Iend+1,j)
#   endif
        END DO

#  else
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO j=J_RANGE
          vbar_stokes(Iend+1,j)=gamma2(ng)*vbar_stokes(Iend,j)
#   ifdef MASKING
          vbar_stokes(Iend+1,j)=vbar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%vmask(Iend+1,j)
#   endif
#   ifdef WET_DRY
          vbar_stokes(Iend+1,j)=vbar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%vmask_wet(Iend+1,j)
#   endif
        END DO
#   undef J_RANGE
#  endif
      END IF
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        vbar_stokes(Istr-1,Jstr)=0.5_r8*(vbar_stokes(Istr  ,Jstr  )+    &
     &                                   vbar_stokes(Istr-1,Jstr+1))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        vbar_stokes(Iend+1,Jstr)=0.5_r8*(vbar_stokes(Iend  ,Jstr  )+    &
     &                                   vbar_stokes(Iend+1,Jstr+1))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        vbar_stokes(Istr-1,Jend+1)=0.5_r8*(vbar_stokes(Istr-1,Jend  )+  &
     &                                     vbar_stokes(Istr  ,Jend+1))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        vbar_stokes(Iend+1,Jend+1)=0.5_r8*(vbar_stokes(Iend+1,Jend  )+  &
     &                                     vbar_stokes(Iend  ,Jend+1))
      END IF
# endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif

      RETURN
      END SUBROUTINE vs2dbc_tile
#endif
      END MODULE vs2dbc_mod
