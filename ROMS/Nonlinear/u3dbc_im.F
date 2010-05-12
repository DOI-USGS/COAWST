#include "cppdefs.h"
      MODULE u3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: u3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  U-velocity.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: u3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE u3dbc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL u3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng), nout,                                  &
     &                 OCEAN(ng) % u)
      RETURN
      END SUBROUTINE u3dbc
!
!***********************************************************************
      SUBROUTINE u3dbc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nout,                                &
     &                       u)
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
      integer, intent(in) :: nstp, nout
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,UBk,2)
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

#  if defined WEST_M3RADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=u(Istr  ,j  ,k,nstp)-                        &
     &                     u(Istr  ,j-1,k,nstp)
            grad(Istr+1,j)=u(Istr+1,j  ,k,nstp)-                        &
     &                     u(Istr+1,j-1,k,nstp)
          END DO
          DO j=Jstr,Jend
            dUdt=u(Istr+1,j,k,nstp)-u(Istr+1,j,k,nout)
            dUdx=u(Istr+1,j,k,nout)-u(Istr+2,j,k,nout)
#   ifdef WEST_M3NUDGING
            IF ((dUdt*dUdx).lt.0.0_r8) THEN
              tau=M3obc_in(ng,iwest)
            ELSE
              tau=M3obc_out(ng,iwest)
            END IF
            tau=tau*dt(ng)
#   endif
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
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
            BOUNDARY(ng)%u_west_Cx(j,k)=Cx
            BOUNDARY(ng)%u_west_Ce(j,k)=Ce
            BOUNDARY(ng)%u_west_C2(j,k)=cff
#   endif
            u(Istr,j,k,nout)=(cff*u(Istr  ,j,k,nstp)+                   &
     &                        Cx *u(Istr+1,j,k,nout)-                   &
     &                        MAX(Ce,0.0_r8)*grad(Istr,j  )-            &
     &                        MIN(Ce,0.0_r8)*grad(Istr,j+1))/           &
     &                       (cff+Cx)
#   ifdef WEST_M3NUDGING
            u(Istr,j,k,nout)=u(Istr,j,k,nout)+                          &
     &                       tau*(BOUNDARY(ng)%u_west(j,k)-             &
     &                            u(Istr,j,k,nstp))
#   endif
#   ifdef MASKING
            u(Istr,j,k,nout)=u(Istr,j,k,nout)*                          &
     &                       GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WET_DRY
            u(Istr,j,k,nout)=u(Istr,j,k,nout)*                          &
     &                       GRID(ng)%umask_wet(Istr,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Istr,j,k,nout)=BOUNDARY(ng)%u_west(j,k)
#   ifdef MASKING
            u(Istr,j,k,nout)=u(Istr,j,k,nout)*                          &
     &                       GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WET_DRY
            u(Istr,j,k,nout)=u(Istr,j,k,nout)*                          &
     &                       GRID(ng)%umask_wet(Istr,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Istr,j,k,nout)=u(Istr+1,j,k,nout)
#   ifdef MASKING
            u(Istr,j,k,nout)=u(Istr,j,k,nout)*                                &
     &                       GRID(ng)%umask(Istr,j)
#   endif
#   ifdef WET_DRY
            u(Istr,j,k,nout)=u(Istr,j,k,nout)*                                &
     &                       GRID(ng)%umask_wet(Istr,j)
#   endif
          END DO
        END DO

#  else
!
!  Western edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Istr,j,k,nout)=0.0_r8
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

#  if defined EAST_M3RADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=u(Iend  ,j  ,k,nstp)-                        &
     &                     u(Iend  ,j-1,k,nstp)
            grad(Iend+1,j)=u(Iend+1,j  ,k,nstp)-                        &
     &                     u(Iend+1,j-1,k,nstp)
          END DO
          DO j=Jstr,Jend
            dUdt=u(Iend,j,k,nstp)-u(Iend  ,j,k,nout)
            dUdx=u(Iend,j,k,nout)-u(Iend-1,j,k,nout)
#   ifdef EAST_M3NUDGING
            IF ((dUdt*dUdx).lt.0.0_r8) THEN
              tau=M3obc_in(ng,ieast)
            ELSE
              tau=M3obc_out(ng,ieast)
            END IF
            tau=tau*dt(ng)
#   endif
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
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
            BOUNDARY(ng)%u_east_Cx(j,k)=Cx
            BOUNDARY(ng)%u_east_Ce(j,k)=Ce
            BOUNDARY(ng)%u_east_C2(j,k)=cff
#   endif
            u(Iend+1,j,k,nout)=(cff*u(Iend+1,j,k,nstp)+                 &
     &                          Cx *u(Iend  ,j,k,nout)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Iend+1,j  )-        &
     &                          MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/       &
     &                         (cff+Cx)
#   ifdef EAST_M3NUDGING
            u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)+                      &
     &                         tau*(BOUNDARY(ng)%u_east(j,k)-           &
     &                              u(Iend+1,j,k,nstp))
#   endif
#   ifdef MASKING
            u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef WET_DRY
            u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%umask_wet(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_M3CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Iend+1,j,k,nout)=BOUNDARY(ng)%u_east(j,k)
#   ifdef MASKING
            u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef WET_DRY
            u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%umask_wet(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_M3GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          u(Iend+1,j,k,nout)=u(Iend,j,k,nout)
#   ifdef MASKING
          u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                        &
     &                       GRID(ng)%umask(Iend+1,j)
#   endif
#   ifdef WET_DRY
          u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                        &
     &                       GRID(ng)%umask_wet(Iend+1,j)
#   endif
        END DO
      END DO

#  else
!
!  Eastern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            u(Iend+1,j,k,nout)=0.0_r8
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

#  if defined SOUTH_M3RADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=u(i+1,Jstr-1,k,nstp)-                        &
     &                     u(i  ,Jstr-1,k,nstp)
            grad(i,Jstr  )=u(i+1,Jstr  ,k,nstp)-                        &
     &                     u(i  ,Jstr  ,k,nstp)
          END DO
          DO i=IstrU,Iend
            dUdt=u(i,Jstr,k,nstp)-u(i,Jstr  ,k,nout)
            dUde=u(i,Jstr,k,nout)-u(i,Jstr+1,k,nout)
#   ifdef SOUTH_M3NUDGING
            IF ((dUdt*dUde).lt.0.0_r8) THEN
              tau=M3obc_in(ng,isouth)
            ELSE
              tau=M3obc_out(ng,isouth)
            END IF
            tau=tau*dt(ng)
#   endif
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
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
            BOUNDARY(ng)%u_south_Cx(i,k)=Cx
            BOUNDARY(ng)%u_south_Ce(i,k)=Ce
            BOUNDARY(ng)%u_south_C2(i,k)=cff
#   endif
            u(i,Jstr-1,k,nout)=(cff*u(i,Jstr-1,k,nstp)+                 &
     &                          Ce *u(i,Jstr  ,k,nout)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/       &
     &                         (cff+Ce)
#   ifdef SOUTH_M3NUDGING
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)+                      &
     &                         tau*(BOUNDARY(ng)%u_south(i,k)-          &
     &                              u(i,Jstr-1,k,nstp))
#   endif
#   ifdef MASKING
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                      &
     &                         GRID(ng)%umask_wet(i,Jstr-1)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Jstr-1,k,nout)=BOUNDARY(ng)%u_south(i,k)
#   ifdef MASKING
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                      &
     &                         GRID(ng)%umask_wet(i,Jstr-1)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Jstr-1,k,nout)=u(i,Jstr,k,nout)
#   ifdef MASKING
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_MASK
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                      &
     &                         GRID(ng)%umask_wet(i,Jstr-1)
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
            u(i,Jstr-1,k,nout)=gamma2(ng)*u(i,Jstr,k,nout)
#   ifdef MASKING
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
#   ifdef WET_DRY
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                      &
     &                         GRID(ng)%umask_wet(i,Jstr-1)
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

#  if defined NORTH_M3RADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            grad(i,Jend  )=u(i+1,Jend  ,k,nstp)-                        &
     &                     u(i  ,Jend  ,k,nstp)
            grad(i,Jend+1)=u(i+1,Jend+1,k,nstp)-                        &
     &                     u(i  ,Jend+1,k,nstp)
          END DO
          DO i=IstrU,Iend
            dUdt=u(i,Jend,k,nstp)-u(i,Jend  ,k,nout)
            dUde=u(i,Jend,k,nout)-u(i,Jend-1,k,nout)
#   ifdef NORTH_M3NUDGING
            IF ((dUdt*dUde).lt.0.0_r8) THEN
              tau=M3obc_in(ng,inorth)
            ELSE
              tau=M3obc_out(ng,inorth)
            END IF
            tau=tau*dt(ng)
#   endif
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
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
            BOUNDARY(ng)%u_north_Cx(i,k)=Cx
            BOUNDARY(ng)%u_north_Ce(i,k)=Ce
            BOUNDARY(ng)%u_north_C2(i,k)=cff
#   endif
            u(i,Jend+1,k,nout)=(cff*u(i,Jend+1,k,nstp)+                 &
     &                          Ce *u(i,Jend  ,k,nout)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/       &
     &                         (cff+Ce)
#   ifdef NORTH_M3NUDGING
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)+                      &
     &                         tau*(BOUNDARY(ng)%u_north(i,k)-          &
     &                              u(i,Jend+1,k,nstp))
#    endif
#   ifdef MASKING
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%umask_wet(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Jend+1,k,nout)=BOUNDARY(ng)%u_north(i,k)
#   ifdef MASKING
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%umask_wet(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            u(i,Jend+1,k,nout)=u(i,Jend,k,nout)
#   ifdef MASKING
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%umask_wet(i,Jend+1)
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
            u(i,Jend+1,k,nout)=gamma2(ng)*u(i,Jend,k,nout)
#   ifdef MASKING
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%umask_wet(i,Jend+1)
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
          u(Istr,Jstr-1,k,nout)=0.5_r8*(u(Istr+1,Jstr-1,k,nout)+        &
     &                                  u(Istr  ,Jstr  ,k,nout))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          u(Iend+1,Jstr-1,k,nout)=0.5_r8*(u(Iend  ,Jstr-1,k,nout)+      &
     &                                    u(Iend+1,Jstr  ,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
          u(Istr,Jend+1,k,nout)=0.5_r8*(u(Istr  ,Jend  ,k,nout)+        &
     &                                  u(Istr+1,Jend+1,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          u(Iend+1,Jend+1,k,nout)=0.5_r8*(u(Iend+1,Jend  ,k,nout)+      &
     &                                    u(Iend  ,Jend+1,k,nout))
        END DO
      END IF
# endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif


# ifdef REFINED_GRID_BC
      IF (ng.ge.2) THEN
        ngp=myparent(ng)
        fac1=time(ngp)+dt(ngp)-time(ng)-dt(ng)
        fac2=time(ng)-time(ngp)+dt(ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
        IF (iic(1).eq.ntstart(1)) THEN
          t1=1
          t2=2
        ELSE
          t1=nstp
          t2=nnew(ng)
        END IF
        IF (WESTERN_EDGE) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=u(Istr  ,j  ,k,nstp)-                        &
     &                     u(Istr  ,j-1,k,nstp)
            grad(Istr+1,j)=u(Istr+1,j  ,k,nstp)-                        &
     &                     u(Istr+1,j-1,k,nstp)
          END DO
          DO j=Jstr,Jend
            dUdt=u(Istr+1,j,k,nstp)-u(Istr+1,j,k,nout)
            dUdx=u(Istr+1,j,k,nout)-u(Istr+2,j,k,nout)
!#   ifdef WEST_M3NUDGING
            IF ((dUdt*dUdx).lt.0.0_r8) THEN
              tau=M3obc_in(ng,iwest)
            ELSE
              tau=M3obc_out(ng,iwest)
            END IF
            tau=tau*dt(ng)
!#   endif
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.0.0_r8) THEN
              dUde=grad(Istr+1,j  )
            ELSE
              dUde=grad(Istr+1,j+1)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=dUdt*dUdx
!#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dUdt*dUde,-cff))
!#   else
!            Ce=0.0_r8
!#   endif
            u(Istr,j,k,nout)=(cff*u(Istr  ,j,k,nstp)+                   &
     &                        Cx *u(Istr+1,j,k,nout)-                   &
     &                        MAX(Ce,0.0_r8)*grad(Istr,j  )-            &
     &                        MIN(Ce,0.0_r8)*grad(Istr,j+1))/           &
     &                       (cff+Cx)
!#   ifdef WEST_M3NUDGING
            my_val=fac1*REFINED(ng)%uR_west(4,j,k,t1)+                  &
     &             fac2*REFINED(ng)%uR_west(4,j,k,t2)
            u(Istr,j,k,nout)=u(Istr,j,k,nout)+                          &
     &                       tau*(my_val-u(Istr,j,k,nstp))
!#   endif
#   ifdef MASKING
            u(Istr,j,k,nout)=u(Istr,j,k,nout)*                          &
     &                       GRID(ng)%umask(Istr,j)
#   endif
          END DO
        END DO
        END IF
        IF (EASTERN_EDGE) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=u(Iend  ,j  ,k,nstp)-                        &
     &                     u(Iend  ,j-1,k,nstp)
            grad(Iend+1,j)=u(Iend+1,j  ,k,nstp)-                        &
     &                     u(Iend+1,j-1,k,nstp)
          END DO
          DO j=Jstr,Jend
            dUdt=u(Iend,j,k,nstp)-u(Iend  ,j,k,nout)
            dUdx=u(Iend,j,k,nout)-u(Iend-1,j,k,nout)
!#   ifdef EAST_M3NUDGING
            IF ((dUdt*dUdx).lt.0.0_r8) THEN
              tau=M3obc_in(ng,ieast)
            ELSE
              tau=M3obc_out(ng,ieast)
            END IF
            tau=tau*dt(ng)
!#   endif
            IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
              dUde=grad(Iend,j  )
            ELSE
              dUde=grad(Iend,j+1)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
            Cx=dUdt*dUdx
!#   ifdef RADIATION_2D
            Ce=MIN(cff,MAX(dUdt*dUde,-cff))
!#   else
!            Ce=0.0_r8
!#   endif
            u(Iend+1,j,k,nout)=(cff*u(Iend+1,j,k,nstp)+                 &
     &                          Cx *u(Iend  ,j,k,nout)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Iend+1,j  )-        &
     &                          MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/       &
     &                         (cff+Cx)
!#   ifdef EAST_M3NUDGING
            my_val=fac1*REFINED(ng)%uR_east(1,j,k,t1)+                  &
     &             fac2*REFINED(ng)%uR_east(1,j,k,t2)
            u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)+                      &
     &                         tau*(my_val-u(Iend+1,j,k,nstp))
!#   endif
#   ifdef MASKING
            u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%umask(Iend+1,j)
#   endif
          END DO
        END DO
        END IF
        IF (SOUTHERN_EDGE) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=u(i+1,Jstr-1,k,nstp)-                        &
     &                     u(i  ,Jstr-1,k,nstp)
            grad(i,Jstr  )=u(i+1,Jstr  ,k,nstp)-                        &
     &                     u(i  ,Jstr  ,k,nstp)
          END DO
          DO i=IstrU,Iend
            dUdt=u(i,Jstr,k,nstp)-u(i,Jstr  ,k,nout)
            dUde=u(i,Jstr,k,nout)-u(i,Jstr+1,k,nout)
!#   ifdef SOUTH_M3NUDGING
            IF ((dUdt*dUde).lt.0.0_r8) THEN
              tau=M3obc_in(ng,isouth)
            ELSE
              tau=M3obc_out(ng,isouth)
            END IF
            tau=tau*dt(ng)
!#   endif
            IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(i-1,Jstr)+grad(i,Jstr))).gt.0.0_r8) THEN
              dUdx=grad(i-1,Jstr)
            ELSE
              dUdx=grad(i  ,Jstr)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
!#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
!#   else
!            Cx=0.0_r8
!#   endif
            Ce=dUdt*dUde
            u(i,Jstr-1,k,nout)=(cff*u(i,Jstr-1,k,nstp)+                 &
     &                          Ce *u(i,Jstr  ,k,nout)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/       &
     &                         (cff+Ce)
!#   ifdef SOUTH_M3NUDGING
            my_val=fac1*REFINED(ng)%uR_south(i,4,k,t1)+                 &
     &             fac2*REFINED(ng)%uR_south(i,4,k,t2)
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)+                      &
     &                         tau*(my_val-u(i,Jstr-1,k,nstp))
!#   endif
#   ifdef MASKING
            u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jstr-1)
#   endif
          END DO
        END DO
        END IF
        IF (NORTHERN_EDGE) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
            grad(i,Jend  )=u(i+1,Jend  ,k,nstp)-                        &
     &                     u(i  ,Jend  ,k,nstp)
            grad(i,Jend+1)=u(i+1,Jend+1,k,nstp)-                        &
     &                     u(i  ,Jend+1,k,nstp)
          END DO
          DO i=IstrU,Iend
            dUdt=u(i,Jend,k,nstp)-u(i,Jend  ,k,nout)
            dUde=u(i,Jend,k,nout)-u(i,Jend-1,k,nout)
!#   ifdef NORTH_M3NUDGING
            IF ((dUdt*dUde).lt.0.0_r8) THEN
              tau=M3obc_in(ng,inorth)
            ELSE
              tau=M3obc_out(ng,inorth)
            END IF
            tau=tau*dt(ng)
!#   endif
            IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
            IF ((dUdt*(grad(i-1,Jend)+grad(i,Jend))).gt.0.0_r8) THEN
              dUdx=grad(i-1,Jend)
            ELSE
              dUdx=grad(i  ,Jend)
            END IF
            cff=MAX(dUdx*dUdx+dUde*dUde,eps)
!#   ifdef RADIATION_2D
            Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
!#   else
!            Cx=0.0_r8
!#   endif
            Ce=dUdt*dUde
            u(i,Jend+1,k,nout)=(cff*u(i,Jend+1,k,nstp)+                 &
     &                          Ce *u(i,Jend  ,k,nout)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/       &
     &                         (cff+Ce)
!#   ifdef NORTH_M3NUDGING
            my_val=fac1*REFINED(ng)%uR_north(i,1,k,t1)+                 &
     &             fac2*REFINED(ng)%uR_north(i,1,k,t2)
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)+                      &
     &                         tau*(my_val-u(i,Jend+1,k,nstp))
!#    endif
#   ifdef MASKING
            u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%umask(i,Jend+1)
#   endif
          END DO
        END DO
        END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          DO k=1,N(ng)
            u(Istr,Jstr-1,k,nout)=0.5_r8*(u(Istr+1,Jstr-1,k,nout)+        &
     &                                    u(Istr  ,Jstr  ,k,nout))
          END DO
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          DO k=1,N(ng)
            u(Iend+1,Jstr-1,k,nout)=0.5_r8*(u(Iend  ,Jstr-1,k,nout)+      &
     &                                      u(Iend+1,Jstr  ,k,nout))
          END DO
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          DO k=1,N(ng)
            u(Istr,Jend+1,k,nout)=0.5_r8*(u(Istr  ,Jend  ,k,nout)+        &
     &                                    u(Istr+1,Jend+1,k,nout))
          END DO
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          DO k=1,N(ng)
            u(Iend+1,Jend+1,k,nout)=0.5_r8*(u(Iend+1,Jend  ,k,nout)+      &
     &                                      u(Iend  ,Jend+1,k,nout))
          END DO
        END IF

      END IF
# endif


      RETURN
      END SUBROUTINE u3dbc_tile
#endif
      END MODULE u3dbc_mod
