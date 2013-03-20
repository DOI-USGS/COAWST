#include "cppdefs.h"
      MODULE us3dbc_mod
#if defined WEC && defined SOLVE3D
!
!svn $Id: us3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  Ustokes-velocity.                                                   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: us3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE us3dbc (ng, tile)
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
      CALL us3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, N(ng),                      &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  OCEAN(ng) % u_stokes)
      RETURN
      END SUBROUTINE us3dbc
!
!***********************************************************************
      SUBROUTINE us3dbc_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, UBk,                  &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        u_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
# ifdef REFINED_GRID
      USE mod_refined
      USE mod_stepping
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: u_stokes(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: u_stokes(LBi:UBi,LBj:UBj,UBk)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dUde, dUdt, dUdx, tau

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad
#ifdef REFINED_GRID
      integer :: t1, t2, ngp
      real(r8) :: fac, fac1, fac2, my_val
#endif

# include "set_bounds.h"

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

#  if defined WEST_M3SRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=u_stokes(Istr  ,j  ,k)-                      &
     &                     u_stokes(Istr  ,j-1,k)
            grad(Istr+1,j)=u_stokes(Istr+1,j  ,k)-                      &
     &                     u_stokes(Istr+1,j-1,k)
          END DO
          DO j=Jstr,Jend
            dUdt=u_stokes(Istr+1,j,k)-u_stokes(Istr+1,j,k)
            dUdx=u_stokes(Istr+1,j,k)-u_stokes(Istr+2,j,k)
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
            u_stokes(Istr,j,k)=(cff*u_stokes(Istr  ,j,k)+               &
     &                         Cx *u_stokes(Istr+1,j,k)-                &
     &                         MAX(Ce,0.0_r8)*grad(Istr,j  )-           &
     &                         MIN(Ce,0.0_r8)*grad(Istr,j+1))/          &
     &                        (cff+Cx)
#   ifdef MASKING
            u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                      &
     &                         GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WET_DRY
            u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                      &
     &                         GRID(ng)%umask_wet(Istr,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3SCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u_stokes(Istr,j,k)=BOUNDARY(ng)%ustokes_west(j,k)
#   ifdef MASKING
            u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                      &
     &                         GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WET_DRY
            u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                      &
     &                         GRID(ng)%umask_wet(Istr,j)
#   endif
          END DO
        END DO

#  else
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u_stokes(Istr,j,k)=u_stokes(Istr+1,j,k)
#   ifdef MASKING
            u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                      &
     &                         GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WET_DRY
            u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                      &
     &                         GRID(ng)%umask_wet(Istr,j)
#   endif
          END DO
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

#  if defined EAST_M3SRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=u_stokes(Iend  ,j  ,k)-                      &
     &                     u_stokes(Iend  ,j-1,k)
            grad(Iend+1,j)=u_stokes(Iend+1,j  ,k)-                      &
     &                     u_stokes(Iend+1,j-1,k)
          END DO
          DO j=Jstr,Jend
            dUdt=u_stokes(Iend,j,k)-u_stokes(Iend  ,j,k)
            dUdx=u_stokes(Iend,j,k)-u_stokes(Iend-1,j,k)
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
              dUde=grad(Iend,j  )
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
            u_stokes(Iend+1,j,k)=(cff*u_stokes(Iend+1,j,k)+             &
     &                           Cx *u_stokes(Iend  ,j,k)-              &
     &                           MAX(Ce,0.0_r8)*grad(Iend+1,j  )-       &
     &                           MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/      &
     &                          (cff+Cx)
#   ifdef MASKING
            u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                  &
     &                           GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef WET_DRY
            u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                  &
     &                           GRID(ng)%umask_wet(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_M3SCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u_stokes(Iend+1,j,k)=BOUNDARY(ng)%ustokes_east(j,k)
#   ifdef MASKING
            u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                  &
     &                           GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef WET_DRY
            u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                   &
     &                           GRID(ng)%umask_wet(Iend+1,j)
#   endif
          END DO
        END DO

#  else
!
!  Eastern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u_stokes(Iend+1,j,k)=u_stokes(Iend,j,k)
#   ifdef MASKING
            u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                  &
     &                           GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef WET_DRY
            u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                  &
     &                           GRID(ng)%umask_wet(Iend+1,j)
#   endif
          END DO
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

#  if defined SOUTH_M3SRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=u_stokes(i+1,Jstr-1,k)-                      &
     &                     u_stokes(i  ,Jstr-1,k)
            grad(i,Jstr  )=u_stokes(i+1,Jstr  ,k)-                      &
     &                     u_stokes(i  ,Jstr  ,k)
          END DO
          DO i=IstrU,Iend
            dUdt=u_stokes(i,Jstr,k)-u_stokes(i,Jstr  ,k)
            dUde=u_stokes(i,Jstr,k)-u_stokes(i,Jstr+1,k)
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
            u_stokes(i,Jstr-1,k)=(cff*u_stokes(i,Jstr-1,k)+             &
     &                           Ce*u_stokes(i,Jstr  ,k)-               &
     &                           MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-       &
     &                           MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/      &
     &                          (cff+Ce)
#   ifdef MASKING
            u_stokes(i,Jstr-1,k)=u_stokes(i,Jstr-1,k)*                  &
     &                           GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
            u_stokes(i,Jstr-1,k)=u_stokes(i,Jstr-1,k)*                  &
     &                           GRID(ng)%umask_wet(i,Jstr-1)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3SCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u_stokes(i,Jstr-1,k)=BOUNDARY(ng)%ustokes_south(i,k)
#   ifdef MASKING
            u_stokes(i,Jstr-1,k)=u_stokes(i,Jstr-1,k)*                  &
     &                           GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
            u_stokes(i,Jstr-1,k)=u_stokes(i,Jstr-1,k)*                  &
     &                           GRID(ng)%umask_wet(i,Jstr-1)
#   endif
          END DO
        END DO

#  else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   if defined EW_PERIODIC || defined COMPOSED_GRID
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=1,N(ng)
          DO i=I_RANGE
            u_stokes(i,Jstr-1,k)=gamma2(ng)*u_stokes(i,Jstr,k)
#   ifdef MASKING
            u_stokes(i,Jstr-1,k)=u_stokes(i,Jstr-1,k)*                  &
     &                           GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
            u_stokes(i,Jstr-1,k)=u_stokes(i,Jstr-1,k)*                  &
     &                           GRID(ng)%umask_wet(i,Jstr-1)
#   endif
          END DO
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

#  if defined NORTH_M3SRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            grad(i,Jend  )=u_stokes(i+1,Jend  ,k)-                      &
     &                     u_stokes(i  ,Jend  ,k)
            grad(i,Jend+1)=u_stokes(i+1,Jend+1,k)-                      &
     &                     u_stokes(i  ,Jend+1,k,nstp)
          END DO
          DO i=IstrU,Iend
            dUdt=u_stokes(i,Jend,k)-u_stokes(i,Jend  ,k)
            dUde=u_stokes(i,Jend,k)-u_stokes(i,Jend-1,k)
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
            u_stokes(i,Jend+1,k)=(cff*u_stokes(i,Jend+1,k)+             &
     &                           Ce *u_stokes(i,Jend  ,k)-              &
     &                           MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-       &
     &                           MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/      &
     &                          (cff+Ce)
#   ifdef MASKING
            u_stokes(i,Jend+1,k)=u_stokes(i,Jend+1,k)*                  &
     &                           GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            u_stokes(i,Jend+1,k)=u_stokes(i,Jend+1,k)*                  &
     &                           GRID(ng)%umask_wet(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3SCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u_stokes(i,Jend+1,k)=BOUNDARY(ng)%ustokes_north(i,k)
#   ifdef MASKING
            u_stokes(i,Jend+1,k)=u_stokes(i,Jend+1,k)*                  &
     &                           GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            u_stokes(i,Jend+1,k)=u_stokes(i,Jend+1,k)*                  &
     &                           GRID(ng)%umask_wet(i,Jend+1)
#   endif
          END DO
        END DO

#  else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#   if defined EW_PERIODIC || defined COMPOSED_GRID
#    define I_RANGE IstrU,Iend
#   else
#    define I_RANGE Istr,IendR
#   endif
        DO k=1,N(ng)
          DO i=I_RANGE
            u_stokes(i,Jend+1,k)=gamma2(ng)*u_stokes(i,Jend,k)
#   ifdef MASKING
            u_stokes(i,Jend+1,k)=u_stokes(i,Jend+1,k)*                  &
     &                           GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            u_stokes(i,Jend+1,k)=u_stokes(i,Jend+1,k)*                  &
     &                           GRID(ng)%umask_wet(i,Jend+1)
#   endif
          END DO
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
        DO k=1,N(ng)
          u_stokes(Istr,Jstr-1,k)=0.5_r8*(u_stokes(Istr+1,Jstr-1,k)+    &
     &                                    u_stokes(Istr  ,Jstr  ,k))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          u_stokes(Iend+1,Jstr-1,k)=0.5_r8*(u_stokes(Iend  ,Jstr-1,k)+  &
     &                                      u_stokes(Iend+1,Jstr  ,k))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
          u_stokes(Istr,Jend+1,k)=0.5_r8*(u_stokes(Istr  ,Jend  ,k)+    &
     &                                    u_stokes(Istr+1,Jend+1,k))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          u_stokes(Iend+1,Jend+1,k)=0.5_r8*(u_stokes(Iend+1,Jend  ,k)+  &
     &                                      u_stokes(Iend  ,Jend+1,k))
        END DO
      END IF
# endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif

      RETURN
      END SUBROUTINE us3dbc_tile
#endif
      END MODULE us3dbc_mod
