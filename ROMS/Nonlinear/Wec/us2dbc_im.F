#include "cppdefs.h"
      MODULE us2dbc_mod
#ifdef WEC
!
!svn $Id: u2sdbc_im.F 779 2008-10-04 23:15:47Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated Ustokes-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: us2dbc, us2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE us2dbc (ng, tile)
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
      CALL us2dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  OCEAN(ng) % ubar_stokes)
      RETURN
      END SUBROUTINE us2dbc
!
!***********************************************************************
      SUBROUTINE us2dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       ubar_stokes)
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
      real(r8), intent(inout) :: ubar_stokes(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ubar_stokes(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff5
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8) :: cff, cff1, cff2, dUde, dUdt, dUdx, tau

# ifdef REFINED_GRID
      integer :: t1, t2, ngp
      real(r8) :: fac, fac1, fac2, my_val
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

#include "set_bounds.h"

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      IF (ng.eq.1) THEN
#endif

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
        DO j=Jstr,Jend+1
          grad(Istr  ,j)=ubar_stokes(Istr  ,j  )-                       &
     &                   ubar_stokes(Istr  ,j-1)
          grad(Istr+1,j)=ubar_stokes(Istr+1,j  )-                       &
     &                   ubar_stokes(Istr+1,j-1)
        END DO
        DO j=Jstr,Jend
          dUdt=ubar_stokes(Istr+1,j)-ubar_stokes(Istr+1,j)
          dUdx=ubar(Istr+1,j)-ubar(Istr+2,j)
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.0.0_r8) THEN
            dUde=grad(Istr+1,j  )
          ELSE
            dUde=grad(Istr+1,j+1)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=dUdt*dUdx
#   ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#   else
          Ce=0.0_r8
#   endif
          ubar_stokes(Istr,j)=(cff*ubar_stokes(Istr  ,j)+               &
     &                        Cx *ubar_stokes(Istr+1,j)-                &
     &                        MAX(Ce,0.0_r8)*grad(Istr,j  )-            &
     &                        MIN(Ce,0.0_r8)*grad(Istr,j+1))/           &
     &                       (cff+Cx)
#   ifdef MASKING
          ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                      &
     &                        GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WET_DRY
          ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                      &
     &                        GRID(ng)%umask_wet(Istr,j)
#   endif
        END DO

#  elif defined WEST_M2SCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ubar_stokes(Istr,j)=BOUNDARY(ng)%ubarstokes_west(j)
#   ifdef MASKING
          ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                      &
     &                        GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WET_DRY
          ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                      &
     &                        GRID(ng)%umask_wet(Istr,j)
#   endif
        END DO

# elif defined WEST_M2SGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ubar_stokes(Istr,j)=ubar_stokes(Istr+1,j)
#   ifdef MASKING
          ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                      &
     &                        GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WET_DRY
          ubar_stokes(Istr,j)=ubar_stokes(Istr,j)*                      &
     &                        GRID(ng)%umask_wet(Istr,j)
#   endif
        END DO

#  else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ubar_stokes(Istr,j)=0.0_r8
        END DO
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
        DO j=Jstr,Jend+1
          grad(Iend  ,j)=ubar_stokes(Iend  ,j  )-                       &
     &                   ubar_stokes(Iend  ,j-1)
          grad(Iend+1,j)=ubar_stokes(Iend+1,j  )-                       &
     &                   ubar_stokes(Iend+1,j-1)
        END DO
        DO j=Jstr,Jend
          dUdt=ubar_stokes(Iend,j)-ubar_stokes(Iend  ,j)
          dUdx=ubar_stokes(Iend,j)-ubar_stokes(Iend-1,j)
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
            dUde=grad(Iend,j)
          ELSE
            dUde=grad(Iend,j+1)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=dUdt*dUdx
#   ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#   else
          Ce=0.0_r8
#   endif
          ubar_stokes(Iend+1,j)=(cff*ubar_stokes(Iend+1,j)+             &
     &                          Cx *ubar_stokes(Iend  ,j)-              &
     &                          MAX(Ce,0.0_r8)*grad(Iend+1,j  )-        &
     &                          MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/       &
     &                         (cff+Cx)
#   ifdef MASKING
          ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef WET_DRY
          ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%umask_wet(Iend+1,j)
#   endif
        END DO

#  elif defined EAST_M2SCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ubar_stokes(Iend+1,j)=BOUNDARY(ng)%ubarstokes_east(j)
#   ifdef MASKING
          ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef WET_DRY
          ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%umask_wet(Iend+1,j)
#   endif
        END DO

#  elif defined EAST_M2SGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ubar_stokes(Iend+1,j)=ubar_stokes(Iend,j)
#   ifdef MASKING
          ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef WET_DRY
          ubar_stokes(Iend+1,j)=ubar_stokes(Iend+1,j)*                  &
     &                          GRID(ng)%umask_wet(Iend+1,j)
#   endif
        END DO

#  else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ubar_stokes(Iend+1,j)=0.0_r8
        END DO
#  endif
      END IF
# endif
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
        DO i=IstrU-1,Iend
          grad(i,Jstr-1)=ubar_stokes(i+1,Jstr-1)-                       &
     &                   ubar_stokes(i  ,Jstr-1)
          grad(i,Jstr  )=ubar_stokes(i+1,Jstr  )-                       &
     &                   ubar_stokes(i  ,Jstr  )
        END DO
        DO i=IstrU,Iend
          dUdt=ubar_stokes(i,Jstr)-ubar_stokes(i,Jstr  )
          dUde=ubar_stokes(i,Jstr)-ubar_stokes(i,Jstr+1)
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,Jstr)+grad(i,Jstr))).gt.0.0_r8) THEN
            dUdx=grad(i-1,Jstr)
          ELSE
            dUdx=grad(i  ,Jstr)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#   ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#   else
          Cx=0.0_r8
#   endif
          Ce=dUdt*dUde
          ubar_stokes(i,Jstr-1)=(cff*ubar_stokes(i,Jstr-1)+             &
     &                          Ce *ubar_stokes(i,Jstr  )-              &
     &                          MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/       &
     &                          (cff+Ce)
#   ifdef MASKING
          ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*                  &
     &                          GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
          ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*                  &
     &                          GRID(ng)%umask_wet(i,Jstr-1)
#  endif
        END DO

#  elif defined SOUTH_M2SCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ubar_stokes(i,Jstr-1)=BOUNDARY(ng)%ubarstokes_south(i)
#   ifdef MASKING
          ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*                  &
     &                          GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
          ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*                  &
     &                          GRID(ng)%umask_wet(i,Jstr-1)
#   endif
        END DO

#  elif defined SOUTH_M2SGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr)
#   ifdef MASKING
          ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*                  &
     &                          GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
          ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*                  &
     &                          GRID(ng)%umask_wet(i,Jstr-1)
#   endif
        END DO

#  else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO i=I_RANGE
          ubar_stokes(i,Jstr-1)=gamma2(ng)*ubar_stokes(i,Jstr)
#   ifdef MASKING
          ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*                  &
     &                          GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
          ubar_stokes(i,Jstr-1)=ubar_stokes(i,Jstr-1)*                  &
     &                          GRID(ng)%umask_wet(i,Jstr-1)
#   endif
        END DO
#   undef I_RANGE
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

#  if defined NORTH_M2RADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=IstrU-1,Iend
          grad(i,Jend  )=ubar_stokes(i+1,Jend  )-                       &
     &                   ubar_stokes(i  ,Jend  )
          grad(i,Jend+1)=ubar_stokes(i+1,Jend+1)-                       &
     &                   ubar_stokes(i  ,Jend+1)
        END DO
        DO i=IstrU,Iend
          dUdt=ubar_stokes(i,Jend)-ubar_stokes(i,Jend  )
          dUde=ubar_stokes(i,Jend)-ubar_stokes(i,Jend-1)
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,Jend)+grad(i,Jend))).gt.0.0_r8) THEN
            dUdx=grad(i-1,Jend)
          ELSE
            dUdx=grad(i  ,Jend)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#   ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#   else
          Cx=0.0_r8
#   endif
          Ce=dUdt*dUde
          ubar_stokes(i,Jend+1)=(cff*ubar_stokes(i,Jend+1)+             &
     &                          Ce *ubar_stokes(i,Jend  )-              &
     &                          MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/       &
     &                          (cff+Ce)
#   ifdef MASKING
          ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
          ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%umask_wet(i,Jend+1)
#   endif
        END DO

#  elif defined NORTH_M2SCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ubar_stokes(i,Jend+1)=BOUNDARY(ng)%ubarstokes_north(i)
#   ifdef MASKING
          ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
          ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%umask_wet(i,Jend+1)
#   endif
        END DO

#  elif defined NORTH_M2SGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend)
#   ifdef MASKING
          ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
          ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%umask_wet(i,Jend+1)
#   endif
        END DO

#  else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   ifdef EW_PERIODIC
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO i=I_RANGE
          ubar_stokes(i,Jend+1)=gamma2(ng)*ubar_stokes(i,Jend)
#   ifdef MASKING
          ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
          ubar_stokes(i,Jend+1)=ubar_stokes(i,Jend+1)*                  &
     &                          GRID(ng)%umask_wet(i,Jend+1)
#   endif
        END DO
#   undef I_RANGE
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
        ubar_stokes(Istr,Jstr-1)=0.5_r8*(ubar_stokes(Istr+1,Jstr-1)+    &
     &                                   ubar_stokes(Istr  ,Jstr  ))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar_stokes(Iend+1,Jstr-1)=0.5_r8*(ubar_stokes(Iend  ,Jstr-1)+  &
     &                                     ubar_stokes(Iend+1,Jstr  ))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        ubar_stokes(Istr,Jend+1)=0.5_r8*(ubar_stokes(Istr  ,Jend  )+    &
     &                                   ubar_stokes(Istr+1,Jend+1))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar_stokes(Iend+1,Jend+1)=0.5_r8*(ubar_stokes(Iend+1,Jend  )+  &
     &                                     ubar_stokes(Iend  ,Jend+1))
      END IF
# endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif
!
      RETURN
      END SUBROUTINE us2dbc_tile
#endif
      END MODULE us2dbc_mod
