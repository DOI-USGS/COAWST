#include "cppdefs.h"
      MODULE v3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: v3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  V-velocity.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: v3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE v3dbc (ng, tile, nout)
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
      CALL v3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng), nout,                                  &
     &                 OCEAN(ng) % v)
      RETURN
      END SUBROUTINE v3dbc
!
!***********************************************************************
      SUBROUTINE v3dbc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nout,                                &
     &                       v)
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
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,UBk,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dVde, dVdt, dVdx, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
#ifdef REFINED_GRID
      integer :: t1, t2, ngp
      real(r8) :: fac, fac1, fac2, my_val
#endif

# include "set_bounds.h"

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

#  if defined SOUTH_M3RADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend+1
            grad(i,Jstr  )=v(i  ,Jstr  ,k,nstp)-                        &
     &                     v(i-1,Jstr  ,k,nstp)
            grad(i,Jstr+1)=v(i  ,Jstr+1,k,nstp)-                        &
     &                     v(i-1,Jstr+1,k,nstp)
          END DO
          DO i=Istr,Iend
            dVdt=v(i,Jstr+1,k,nstp)-v(i,Jstr+1,k,nout)
            dVde=v(i,Jstr+1,k,nout)-v(i,Jstr+2,k,nout)
#   ifdef SOUTH_M3NUDGING
            IF ((dVdt*dVde).lt.0.0_r8) THEN
              tau=M3obc_in(ng,isouth)
            ELSE
              tau=M3obc_out(ng,isouth)
            END IF
            tau=tau*dt(ng)
#   endif
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
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
            BOUNDARY(ng)%v_south_Cx(i,k)=Cx
            BOUNDARY(ng)%v_south_Ce(i,k)=Ce
            BOUNDARY(ng)%v_south_C2(i,k)=cff
#   endif
            v(i,Jstr,k,nout)=(cff*v(i,Jstr  ,k,nstp)+                   &
     &                        Ce *v(i,Jstr+1,k,nout)-                   &
     &                        MAX(Cx,0.0_r8)*grad(i  ,Jstr)-            &
     &                        MIN(Cx,0.0_r8)*grad(i+1,Jstr))/           &
     &                       (cff+Ce)
#   ifdef SOUTH_M3NUDGING
            v(i,Jstr,k,nout)=v(i,Jstr,k,nout)+                          &
     &                       tau*(BOUNDARY(ng)%v_south(i,k)-            &
     &                            v(i,Jstr,k,nstp))
#   endif
#   ifdef MASKING
            v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                          &
     &                       GRID(ng)%vmask(i,Jstr)
#   endif
#   ifdef WET_DRY
            v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                          &
     &                       GRID(ng)%vmask_wet(i,Jstr)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jstr,k,nout)=BOUNDARY(ng)%v_south(i,k)
#   ifdef MASKING
            v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                          &
     &                       GRID(ng)%vmask(i,Jstr)
#   endif
#   ifdef WET_DRY
            v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                          &
     &                       GRID(ng)%vmask_wet(i,Jstr)
#   endif
          END DO
        END DO

#  elif defined SOUTH_M3GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jstr,k,nout)=v(i,Jstr+1,k,nout)
#   ifdef MASKING
            v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                          &
     &                       GRID(ng)%vmask(i,Jstr)
#   endif
#   ifdef WET_DRY
            v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                          &
     &                       GRID(ng)%vmask_wet(i,Jstr)
#   endif
          END DO
        END DO
#  else
!
!  Southern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jstr,k,nout)=0.0_r8
          END DO
        END DO
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
          DO i=Istr,Iend+1
            grad(i,Jend  )=v(i  ,Jend  ,k,nstp)-                        &
     &                     v(i-1,Jend  ,k,nstp)
            grad(i,Jend+1)=v(i  ,Jend+1,k,nstp)-                        &
     &                     v(i-1,Jend+1,k,nstp)
          END DO
          DO i=Istr,Iend
            dVdt=v(i,Jend,k,nstp)-v(i,Jend  ,k,nout)
            dVde=v(i,Jend,k,nout)-v(i,Jend-1,k,nout)
#   ifdef NORTH_M3NUDGING
            IF ((dVdt*dVde).lt.0.0_r8) THEN
              tau=M3obc_in(ng,inorth)
            ELSE
              tau=M3obc_out(ng,inorth)
            END IF
            tau=tau*dt(ng)
#   endif
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
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
            BOUNDARY(ng)%v_north_Cx(i,k)=Cx
            BOUNDARY(ng)%v_north_Ce(i,k)=Ce
            BOUNDARY(ng)%v_north_C2(i,k)=cff
#   endif
            v(i,Jend+1,k,nout)=(cff*v(i,Jend+1,k,nstp)+                 &
     &                          Ce *v(i,Jend  ,k,nout)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-        &
     &                          MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/       &
     &                         (cff+Ce)
#   ifdef NORTH_M3NUDGING
            v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)+                      &
     &                         tau*(BOUNDARY(ng)%v_north(i,k)-          &
     &                              v(i,Jend+1,k,nstp))
#   endif
#   ifdef MASKING
            v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%vmask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jend+1,k,nout)=BOUNDARY(ng)%v_north(i,k)
#   ifdef MASKING
            v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%vmask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_M3GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jend+1,k,nout)=v(i,Jend,k,nout)
#   ifdef MASKING
            v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%vmask(i,Jend+1)
#   endif
#   ifdef WET_DRY
            v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(i,Jend+1)
#   endif
          END DO
        END DO

#  else
!
!  Northern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            v(i,Jend+1,k,nout)=0.0_r8
          END DO
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

#  if defined WEST_M3RADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend
            grad(Istr-1,j)=v(Istr-1,j+1,k,nstp)-                        &
     &                     v(Istr-1,j  ,k,nstp)
            grad(Istr  ,j)=v(Istr  ,j+1,k,nstp)-                        &
     &                     v(Istr  ,j  ,k,nstp)
          END DO
          DO j=JstrV,Jend
            dVdt=v(Istr,j,k,nstp)-v(Istr  ,j,k,nout)
            dVdx=v(Istr,j,k,nout)-v(Istr+1,j,k,nout)
#   ifdef WEST_M3NUDGING
            IF ((dVdt*dVdx).lt.0.0_r8) THEN
              tau=M3obc_in(ng,iwest)
            ELSE
              tau=M3obc_out(ng,iwest)
            END IF
            tau=tau*dt(ng)
#   endif
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
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
            BOUNDARY(ng)%v_west_Cx(j,k)=Cx
            BOUNDARY(ng)%v_west_Ce(j,k)=Ce
            BOUNDARY(ng)%v_west_C2(j,k)=cff
#   endif
            v(Istr-1,j,k,nout)=(cff*v(Istr-1,j,k,nstp)+                 &
     &                          Cx *v(Istr  ,j,k,nout)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-        &
     &                          MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/       &
     &                         (cff+Cx)
#   ifdef WEST_M3NUDGING
            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)+                      &
     &                         tau*(BOUNDARY(ng)%v_west(j,k)-           &
     &                              v(Istr-1,j,k,nstp))
#    endif
#   ifdef MASKING
            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
     &                         GRID(ng)%vmask(Istr-1,j)
#   endif
#   ifdef WET_DRY
            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Istr-1,j,k,nout)=BOUNDARY(ng)%v_west(j,k)
#   ifdef MASKING
            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
     &                         GRID(ng)%vmask(Istr-1,j)
#   endif
#   ifdef WET_DRY
            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_M3GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Istr-1,j,k,nout)=v(Istr,j,k,nout)
#   ifdef MASKING
            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
     &                         GRID(ng)%vmask(Istr-1,j)
#   endif
#   ifdef WET_DRY
            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(Istr-1,j)
#   endif
          END DO
        END DO

#  else
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   if defined NS_PERIODIC || defined COMPOSED_GRID
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=1,N(ng)
          DO j=J_RANGE
            v(Istr-1,j,k,nout)=gamma2(ng)*v(Istr,j,k,nout)
#   ifdef MASKING
            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
     &                         GRID(ng)%vmask(Istr-1,j)
#   endif
#   ifdef WET_DRY
            v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(Istr-1,j)
#   endif
          END DO
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

#  if defined EAST_M3RADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend
            grad(Iend  ,j)=v(Iend  ,j+1,k,nstp)-                        &
     &                     v(Iend  ,j  ,k,nstp)
            grad(Iend+1,j)=v(Iend+1,j+1,k,nstp)-                        &
     &                     v(Iend+1,j  ,k,nstp)
          END DO
          DO j=JstrV,Jend
            dVdt=v(Iend,j,k,nstp)-v(Iend  ,j,k,nout)
            dVdx=v(Iend,j,k,nout)-v(Iend-1,j,k,nout)
#   ifdef EAST_M3NUDGING
            IF ((dVdt*dVdx).lt.0.0_r8) THEN
              tau=M3obc_in(ng,ieast)
            ELSE
              tau=M3obc_out(ng,ieast)
            END IF
            tau=tau*dt(ng)
#   endif
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
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
            BOUNDARY(ng)%v_east_Cx(j,k)=Cx
            BOUNDARY(ng)%v_east_Ce(j,k)=Ce
            BOUNDARY(ng)%v_east_C2(j,k)=cff
#   endif
            v(Iend+1,j,k,nout)=(cff*v(Iend+1,j,k,nstp)+                 &
     &                          Cx *v(Iend  ,j,k,nout)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-        &
     &                          MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/       &
     &                         (cff+Cx)
#   ifdef EAST_M3NUDGING
            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)+                      &
     &                         tau*(BOUNDARY(ng)%v_east(j,k)-           &
     &                              v(Iend+1,j,k,nstp))
#   endif
#   ifdef MASKING
            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%vmask(Iend+1,j)
#   endif
#   ifdef WET_DRY
            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_M3CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Iend+1,j,k,nout)=BOUNDARY(ng)%v_east(j,k)
#   ifdef MASKING
            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%vmask(Iend+1,j)
#   endif
#   ifdef WET_DRY
            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_M3GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=JstrV,Jend
            v(Iend+1,j,k,nout)=v(Iend,j,k,nout)
#   ifdef MASKING
            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%vmask(Iend+1,j)
#   endif
#   ifdef WET_DRY
            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(Iend+1,j)
#   endif
          END DO
        END DO

#  else
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   if defined NS_PERIODIC || defined COMPOSED_GRID
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO k=1,N(ng)
          DO j=J_RANGE
            v(Iend+1,j,k,nout)=gamma2(ng)*v(Iend,j,k,nout)
#   ifdef MASKING
            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%vmask(Iend+1,j)
#   endif
#   ifdef WET_DRY
            v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                      &
     &                         GRID(ng)%vmask_wet(Iend+1,j)
#   endif
          END DO
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
        DO k=1,N(ng)
          v(Istr-1,Jstr,k,nout)=0.5_r8*(v(Istr  ,Jstr  ,k,nout)+        &
     &                                  v(Istr-1,Jstr+1,k,nout))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          v(Iend+1,Jstr,k,nout)=0.5_r8*(v(Iend  ,Jstr  ,k,nout)+        &
     &                                  v(Iend+1,Jstr+1,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
          v(Istr-1,Jend+1,k,nout)=0.5_r8*(v(Istr-1,Jend  ,k,nout)+      &
     &                                    v(Istr  ,Jend+1,k,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          v(Iend+1,Jend+1,k,nout)=0.5_r8*(v(Iend+1,Jend  ,k,nout)+      &
     &                                    v(Iend  ,Jend+1,k,nout))
        END DO
      END IF
# endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif

# ifdef REFINED_GRID_BC
!
! For child grid, impose bc's to reduce noise reflections of a 
! clamped condition.
!
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
        IF (SOUTHERN_EDGE) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jstr  )=v(i  ,Jstr  ,k,nstp)-                      &
     &                       v(i-1,Jstr  ,k,nstp)
#   ifdef MASKING
              grad(i,Jstr)=grad(i,Jstr)*GRID(ng)%vmask(i  ,Jstr)*       &
     &                                  GRID(ng)%vmask(i-1,Jstr)
#   endif
              grad(i,Jstr+1)=v(i  ,Jstr+1,k,nstp)-                      &
     &                       v(i-1,Jstr+1,k,nstp)
#   ifdef MASKING
              grad(i,Jstr+1)=grad(i,Jstr+1)*GRID(ng)%vmask(i  ,Jstr+1)* &
     &                                      GRID(ng)%vmask(i-1,Jstr+1)
#   endif
            END DO
            DO i=Istr,Iend
              dVdt=v(i,Jstr+1,k,nstp)-v(i,Jstr+1,k,nout)
              dVde=v(i,Jstr+1,k,nout)-v(i,Jstr+2,k,nout)
!
!  Set tau for nudging to parent data.
!
!             IF ((dVdt*dVde).lt.0.0_r8) THEN
!               tau=M3obc_in(ng,isouth)
!             ELSE
!               tau=M3obc_out(ng,isouth)
!             END IF
!             tau=tau*dt(ng)
              tau=1.0_r8
              IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(i,Jstr+1)+grad(i+1,Jstr+1))).gt.          &
     &             0.0_r8) THEN
                dVdx=grad(i  ,Jstr+1)
              ELSE
                dVdx=grad(i+1,Jstr+1)
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
              Ce=dVdt*dVde
              v(i,Jstr,k,nout)=(cff*v(i,Jstr  ,k,nstp)+                 &
     &                          Ce *v(i,Jstr+1,k,nout)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i  ,Jstr)-          &
     &                          MIN(Cx,0.0_r8)*grad(i+1,Jstr))/         &
     &                         (cff+Ce)
              my_val=fac1*REFINED(ng)%vR_south(i,4,k,t1)+               &
     &               fac2*REFINED(ng)%vR_south(i,4,k,t2)
              v(i,Jstr,k,nout)=v(i,Jstr,k,nout)+                        &
     &                         tau*(my_val-                             &
     &                              v(i,Jstr,k,nstp))
#   ifdef MASKING
              v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                        &
     &                         GRID(ng)%vmask(i,Jstr)
#   endif
#   ifdef WET_DRY
              v(i,Jstr,k,nout)=v(i,Jstr,k,nout)*                        &
     &                         GRID(ng)%vmask_wet(i,Jstr)
#   endif
            END DO
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jend  )=v(i  ,Jend  ,k,nstp)-                      &
     &                       v(i-1,Jend  ,k,nstp)
#   ifdef MASKING
              grad(i,Jend)=grad(i,Jend)*GRID(ng)%vmask(i  ,Jend)*       &
     &                                  GRID(ng)%vmask(i-1,Jend)
#   endif
              grad(i,Jend+1)=v(i  ,Jend+1,k,nstp)-                      &
     &                       v(i-1,Jend+1,k,nstp)
#   ifdef MASKING
              grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%vmask(i  ,Jend+1)* &
         &                                  GRID(ng)%vmask(i-1,Jend+1)
#   endif
            END DO
            DO i=Istr,Iend
              dVdt=v(i,Jend,k,nstp)-v(i,Jend  ,k,nout)
              dVde=v(i,Jend,k,nout)-v(i,Jend-1,k,nout)
!
!  Set tau for nudging to parent data.
!
!             IF ((dVdt*dVde).lt.0.0_r8) THEN
!               tau=M3obc_in(ng,inorth)
!             ELSE
!               tau=M3obc_out(ng,inorth)
!             END IF
!             tau=tau*dt(ng)
              tau=1.0_r8
              IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
                dVdx=grad(i  ,Jend)
              ELSE
                dVdx=grad(i+1,Jend)
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=MIN(cff,MAX(dVdt*dVdx,-cff))
              Ce=dVdt*dVde
              v(i,Jend+1,k,nout)=(cff*v(i,Jend+1,k,nstp)+               &
     &                            Ce *v(i,Jend  ,k,nout)-               &
     &                            MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-      &
     &                            MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/     &
     &                           (cff+Ce)
              my_val=fac1*REFINED(ng)%vR_north(i,1,k,t1)+               &
     &               fac2*REFINED(ng)%vR_north(i,1,k,t2)
              v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)+                    &
     &                           tau*(my_val-                           &
     &                                v(i,Jend+1,k,nstp))
#   ifdef MASKING
              v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                    &
     &                           GRID(ng)%vmask(i,Jend+1)
#   endif
#   ifdef WET_DRY
              v(i,Jend+1,k,nout)=v(i,Jend+1,k,nout)*                    &
     &                           GRID(ng)%vmask_wet(i,Jend+1)
#   endif
            END DO
          END DO
        END IF
        IF (WESTERN_EDGE) THEN
!
!  Western edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              grad(Istr-1,j)=v(Istr-1,j+1,k,nstp)-                      &
     &                       v(Istr-1,j  ,k,nstp)
              grad(Istr  ,j)=v(Istr  ,j+1,k,nstp)-                      &
     &                       v(Istr  ,j  ,k,nstp)
            END DO
            DO j=JstrV,Jend
              dVdt=v(Istr,j,k,nstp)-v(Istr  ,j,k,nout)
              dVdx=v(Istr,j,k,nout)-v(Istr+1,j,k,nout)
!
!  Set tau for nudging to parent data.
!
!             IF ((dVdt*dVdx).lt.0.0_r8) THEN
!               tau=M3obc_in(ng,iwest)
!             ELSE
!               tau=M3obc_out(ng,iwest)
!             END IF
!             tau=tau*dt(ng)
              tau=1.0_r8
              IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(Istr,j-1)+grad(Istr,j))).gt.0.0_r8) THEN
                dVde=grad(Istr,j-1)
              ELSE
                dVde=grad(Istr,j  )
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=dVdt*dVdx
              Ce=MIN(cff,MAX(dVdt*dVde,-cff))
              v(Istr-1,j,k,nout)=(cff*v(Istr-1,j,k,nstp)+               &
     &                            Cx *v(Istr  ,j,k,nout)-               &
     &                            MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-      &
     &                            MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/     &
     &                           (cff+Cx)
              my_val=fac1*REFINED(ng)%vR_west(4,j,k,t1)+                &
     &               fac2*REFINED(ng)%vR_west(4,j,k,t2)
              v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)+                    &
     &                           tau*(my_val-                           &
     &                                v(Istr-1,j,k,nstp))
#   ifdef MASKING
              v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                    &
     &                           GRID(ng)%vmask(Istr-1,j)
#   endif
#   ifdef WET_DRY
              v(Istr-1,j,k,nout)=v(Istr-1,j,k,nout)*                    &
     &                           GRID(ng)%vmask_wet(Istr-1,j)
#   endif
            END DO
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              grad(Iend  ,j)=v(Iend  ,j+1,k,nstp)-                      &
     &                       v(Iend  ,j  ,k,nstp)
              grad(Iend+1,j)=v(Iend+1,j+1,k,nstp)-                      &
     &                       v(Iend+1,j  ,k,nstp)
            END DO
            DO j=JstrV,Jend
              dVdt=v(Iend,j,k,nstp)-v(Iend  ,j,k,nout)
              dVdx=v(Iend,j,k,nout)-v(Iend-1,j,k,nout)
!
!  Set tau for nudging to parent data.
!
!             IF ((dVdt*dVdx).lt.0.0_r8) THEN
!               tau=M3obc_in(ng,ieast)
!             ELSE
!               tau=M3obc_out(ng,ieast)
!             END IF
!             tau=tau*dt(ng)
              tau=1.0_r8
              IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(Iend,j-1)+grad(Iend,j))).gt.0.0_r8) THEN
                dVde=grad(Iend,j-1)
              ELSE
                dVde=grad(Iend,j  )
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=dVdt*dVdx
              Ce=MIN(cff,MAX(dVdt*dVde,-cff))
              v(Iend+1,j,k,nout)=(cff*v(Iend+1,j,k,nstp)+               &
     &                            Cx *v(Iend  ,j,k,nout)-               &
     &                            MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-      &
     &                            MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/     &
     &                           (cff+Cx)
              my_val=fac1*REFINED(ng)%vR_east(1,j,k,t1)+                &
     &               fac2*REFINED(ng)%vR_east(1,j,k,t2)
              v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)+                    &
     &                           tau*(my_val-                           &
     &                                v(Iend+1,j,k,nstp))
#   ifdef MASKING
              v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                    &
     &                           GRID(ng)%vmask(Iend+1,j)
#   endif
#   ifdef WET_DRY
              v(Iend+1,j,k,nout)=v(Iend+1,j,k,nout)*                    &
     &                           GRID(ng)%vmask_wet(Iend+1,j)
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
            v(Istr-1,Jstr,k,nout)=0.5_r8*(v(Istr  ,Jstr  ,k,nout)+      &
     &                                    v(Istr-1,Jstr+1,k,nout))
          END DO
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          DO k=1,N(ng)
            v(Iend+1,Jstr,k,nout)=0.5_r8*(v(Iend  ,Jstr  ,k,nout)+      &
     &                                    v(Iend+1,Jstr+1,k,nout))
          END DO
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          DO k=1,N(ng)
            v(Istr-1,Jend+1,k,nout)=0.5_r8*(v(Istr-1,Jend  ,k,nout)+    &
     &                                      v(Istr  ,Jend+1,k,nout))
          END DO
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          DO k=1,N(ng)
            v(Iend+1,Jend+1,k,nout)=0.5_r8*(v(Iend+1,Jend  ,k,nout)+    &
     &                                      v(Iend  ,Jend+1,k,nout))
          END DO
        END IF
      END IF
# endif

      RETURN
      END SUBROUTINE v3dbc_tile
#endif
      END MODULE v3dbc_mod
