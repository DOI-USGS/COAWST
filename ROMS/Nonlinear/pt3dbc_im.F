#include "cppdefs.h"
      MODULE pt3dbc_mod
#ifdef SOLVE3D
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for the ITRC-th    !
!  tracer field.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: pt3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE pt3dbc (ng, tile, nout, itrc)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout, itrc
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL pt3dbc_tile (ng, tile, itrc,                                 &
     &                  LBi, UBi, LBj, UBj, N(ng), NT(ng),              &
     &                  IminS, ImaxS, JminS, JmaxS,                      &
     &                  nstp(ng), nout,                                 &
     &                  OCEAN(ng)% t)
      RETURN
      END SUBROUTINE pt3dbc

!
!***********************************************************************
      SUBROUTINE pt3dbc_tile (ng, tile, itrc,                           &
     &                        LBi, UBi, LBj, UBj, UBk, UBt,             &
     &                        IminS, ImaxS, JminS, JmaxS,                &
     &                        nstp, nout,                               &
     &                       t)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: eps =1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

# include "set_bounds.h"

      IF (.not.EWperiodic(ng)) THEN
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN

# if defined WEST_PTRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
              grad(Istr-1,j)=t(Istr-1,j  ,k,nstp,itrc)-                 &
     &                       t(Istr-1,j-1,k,nstp,itrc)
#   ifdef MASKING
              grad(Istr-1,j)=grad(Istr-1,j)*                            &
     &                       GRID(ng)%vmask(Istr-1,j)
#   endif
              grad(Istr  ,j)=t(Istr  ,j  ,k,nstp,itrc)-                 &
     &                       t(Istr  ,j-1,k,nstp,itrc)
#   ifdef MASKING
              grad(Istr  ,j)=grad(Istr  ,j)*                            &
     &                       GRID(ng)%vmask(Istr  ,j)
#   endif
            END DO
            DO j=Jstr,Jend
              dTdt=t(Istr,j,k,nstp,itrc)-t(Istr  ,j,k,nout,itrc)
              dTdx=t(Istr,j,k,nout,itrc)-t(Istr+1,j,k,nout,itrc)
#   ifdef WEST_PTNUDGING
              tau=Tobc_out(itrc,ng,iwest)
              IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(itrc,ng,iwest)
              tau=tau*dt(ng)
#   endif
              IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
              IF ((dTdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
                dTde=grad(Istr,j  )
              ELSE
                dTde=grad(Istr,j+1)
              END IF
              cff=MAX(dTdx*dTdx+dTde*dTde,eps)
              Cx=dTdt*dTdx
#   ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dTdt*dTde,-cff))
#   else
              Ce=0.0_r8
#   endif
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%t_west_Cx(j,k,itrc)=Cx
              BOUNDARY(ng)%t_west_Ce(j,k,itrc)=Ce
              BOUNDARY(ng)%t_west_C2(j,k,itrc)=cff
#   endif
              t(Istr-1,j,k,nout,itrc)=(cff*t(Istr-1,j,k,nstp,itrc)+     &
     &                                 Cx *t(Istr  ,j,k,nout,itrc)-     &
     &                                 MAX(Ce,0.0_r8)*                  &
     &                                    grad(Istr-1,j  )-             &
     &                                 MIN(Ce,0.0_r8)*                  &
     &                                    grad(Istr-1,j+1))/            &
     &                                (cff+Cx)
#   ifdef WEST_PTNUDGING
              t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)+          &
     &                              tau*(BOUNDARY(ng)%t_west(j,k,itrc)- &
     &                                   t(Istr-1,j,k,nstp,itrc))
#   endif
#   ifdef MASKING
              t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*          &
     &                                GRID(ng)%rmask(Istr-1,j)
#   endif
            END DO
          END DO

#  elif defined WEST_PTCLAMPED
!
!  Western edge, clamped boundary condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
              t(Istr-1,j,k,nout,itrc)=BOUNDARY(ng)%t_west(j,k,itrc)
#   ifdef MASKING
              t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*          &
     &                              GRID(ng)%rmask(Istr-1,j)
#   endif
            END DO
          END DO

#  elif defined WEST_PTGRADIENT
!
!  Western edge, gradient boundary condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
              t(Istr-1,j,k,nout,itrc)=t(Istr,j,k,nout,itrc)
#   ifdef MASKING
              t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*          &
     &                                GRID(ng)%rmask(Istr-1,j)
#   endif
            END DO
          END DO

#  else
!
!  Western edge, closed boundary condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
              t(Istr-1,j,k,nout,itrc)=t(Istr,j,k,nout,itrc)
#   ifdef MASKING
              t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*          &
     &                                GRID(ng)%rmask(Istr-1,j)
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
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN

#  if defined EAST_PTRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
             grad(Iend  ,j)=t(Iend  ,j  ,k,nstp,itrc)-                  &
     &                      t(Iend  ,j-1,k,nstp,itrc)
#   ifdef MASKING
             grad(Iend  ,j)=grad(Iend  ,j)*                             &
     &                      GRID(ng)%vmask(Iend  ,j)
#   endif
             grad(Iend+1,j)=t(Iend+1,j  ,k,nstp,itrc)-                  &
     &                      t(Iend+1,j-1,k,nstp,itrc)
#   ifdef MASKING
             grad(Iend+1,j)=grad(Iend+1,j)*                             &
     &                      GRID(ng)%vmask(Iend+1,j)
#   endif
            END DO
            DO j=Jstr,Jend
              dTdt=t(Iend,j,k,nstp,itrc)-t(Iend  ,j,k,nout,itrc)
              dTdx=t(Iend,j,k,nout,itrc)-t(Iend-1,j,k,nout,itrc)
#   ifdef EAST_PTNUDGING
              tau=Tobc_out(itrc,ng,ieast)
              IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(itrc,ng,ieast)
              tau=tau*dt(ng)
#   endif
              IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
              IF ((dTdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
                dTde=grad(Iend,j  )
              ELSE
                dTde=grad(Iend,j+1)
              END IF
              cff=MAX(dTdx*dTdx+dTde*dTde,eps)
              Cx=dTdt*dTdx
#   ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dTdt*dTde,-cff))
#   else
              Ce=0.0_r8
#   endif
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%t_east_Cx(j,k,itrc)=Cx
              BOUNDARY(ng)%t_east_Ce(j,k,itrc)=Ce
              BOUNDARY(ng)%t_east_C2(j,k,itrc)=cff
#   endif
              t(Iend+1,j,k,nout,itrc)=(cff*t(Iend+1,j,k,nstp,itrc)+     &
     &                                 Cx *t(Iend  ,j,k,nout,itrc)-     &
     &                                 MAX(Ce,0.0_r8)*                  &
     &                                    grad(Iend+1,j  )-             &
     &                                 MIN(Ce,0.0_r8)*                  &
     &                                    grad(Iend+1,j+1))/            &
     &                                (cff+Cx)
#   ifdef EAST_PTNUDGING
              t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)+          &
     &                              tau*(BOUNDARY(ng)%t_east(j,k,itrc)- &
     &                                   t(Iend+1,j,k,nstp,itrc))
#   endif
#   ifdef MASKING
              t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*          &
     &                                GRID(ng)%rmask(Iend+1,j)
#   endif
            END DO
          END DO

#  elif defined EAST_PTCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
              t(Iend+1,j,k,nout,itrc)=BOUNDARY(ng)%t_east(j,k,itrc)
#   ifdef MASKING
              t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*          &
     &                              GRID(ng)%rmask(Iend+1,j)
#   endif
            END DO
          END DO

#  elif defined EAST_PTGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
              t(Iend+1,j,k,nout,itrc)=t(Iend,j,k,nout,itrc)
#   ifdef MASKING
              t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*          &
     &                              GRID(ng)%rmask(Iend+1,j)
#   endif
            END DO
          END DO

#  else
!
!  Eastern edge, closed boundary condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend
              t(Iend+1,j,k,nout,itrc)=t(Iend,j,k,nout,itrc)
#   ifdef MASKING
              t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*          &
     &                                GRID(ng)%rmask(Iend+1,j)
#   endif
            END DO
          END DO
#  endif
        END IF
      END IF

!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN

#  if defined SOUTH_PTRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jstr  )=t(i  ,Jstr  ,k,nstp,itrc)-                   &
     &                       t(i-1,Jstr  ,k,nstp,itrc)
#   ifdef MASKING
              grad(i,Jstr  )=grad(i,Jstr  )*                              &
     &                       GRID(ng)%umask(i,Jstr  )
#   endif
              grad(i,Jstr-1)=t(i  ,Jstr-1,k,nstp,itrc)-                   &
     &                       t(i-1,Jstr-1,k,nstp,itrc)
#   ifdef MASKING
              grad(i,Jstr-1)=grad(i,Jstr-1)*                              &
     &                       GRID(ng)%umask(i,Jstr-1)
#   endif
            END DO
            DO i=Istr,Iend
              dTdt=t(i,Jstr,k,nstp,itrc)-t(i,Jstr  ,k,nout,itrc)
              dTde=t(i,Jstr,k,nout,itrc)-t(i,Jstr+1,k,nout,itrc)
#   ifdef SOUTH_PTNUDGING
              tau=Tobc_out(itrc,ng,isouth)
              IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(itrc,ng,isouth)
              tau=tau*dt(ng)
#   endif
              IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
              IF ((dTdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
                dTdx=grad(i  ,Jstr)
              ELSE
                dTdx=grad(i+1,Jstr)
              END IF
              cff=MAX(dTdx*dTdx+dTde*dTde,eps)
#   ifdef RADIATION_2D
              Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
#   else
              Cx=0.0_r8
#   endif
              Ce=dTdt*dTde
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%t_south_Cx(i,k,itrc)=Cx
              BOUNDARY(ng)%t_south_Ce(i,k,itrc)=Ce
              BOUNDARY(ng)%t_south_C2(i,k,itrc)=cff
#   endif
              t(i,Jstr-1,k,nout,itrc)=(cff*t(i,Jstr-1,k,nstp,itrc)+       &
     &                                 Ce *t(i,Jstr  ,k,nout,itrc )-      &
     &                                 MAX(Cx,0.0_r8)*                    &
     &                                    grad(i  ,Jstr-1)-               &
     &                                 MIN(Cx,0.0_r8)*                    &
     &                                    grad(i+1,Jstr-1))/              &
     &                                (cff+Ce)
#   ifdef SOUTH_PTNUDGING
              t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)+            &
     &                                tau*(BOUNDARY(ng)%t_south(i,k,itrc)-&
     &                                     t(i,Jstr-1,k,nstp,itrc))
#   endif
#   ifdef MASKING
              t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
     &                                GRID(ng)%rmask(i,Jstr-1)
#   endif
            END DO
          END DO

#  elif defined SOUTH_PTCLAMPED
!
!  Southern edge, clamped boundary condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,Jstr-1,k,nout,itrc)=BOUNDARY(ng)%t_south(i,k,itrc)
#   ifdef MASKING
              t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(i,Jstr-1)
#   endif
            END DO
          END DO

#  elif defined SOUTH_PTGRADIENT
!
!  Southern edge, gradient boundary condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,Jstr-1,k,nout,itrc)=t(i,Jstr,k,nout,itrc)
#   ifdef MASKING
              t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
     &                                GRID(ng)%rmask(i,Jstr-1)
#   endif
            END DO
          END DO

#  else
!
!  Southern edge, closed boundary condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,Jstr-1,k,nout,itrc)=t(i,Jstr,k,nout,itrc)
#   ifdef MASKING
              t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
     &                                GRID(ng)%rmask(i,Jstr-1)
#   endif
            END DO
          END DO
#  endif
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN

#  if defined NORTH_PTRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jend  )=t(i  ,Jend  ,k,nstp,itrc)-                   &
     &                       t(i-1,Jend  ,k,nstp,itrc)
#   ifdef MASKING
              grad(i,Jend  )=grad(i,Jend  )*                              &
     &                       GRID(ng)%umask(i,Jend  )
#   endif
              grad(i,Jend+1)=t(i  ,Jend+1,k,nstp,itrc)-                   &
     &                       t(i-1,Jend+1,k,nstp,itrc)
#   ifdef MASKING
              grad(i,Jend+1)=grad(i,Jend+1)*                              &
     &                       GRID(ng)%umask(i,Jend+1)
#   endif
            END DO
            DO i=Istr,Iend
              dTdt=t(i,Jend,k,nstp,itrc)-t(i,Jend  ,k,nout,itrc)
              dTde=t(i,Jend,k,nout,itrc)-t(i,Jend-1,k,nout,itrc)
#   ifdef NORTH_PTNUDGING
              tau=Tobc_out(itrc,ng,inorth)
              IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(itrc,ng,inorth)
              tau=tau*dt(ng)
#   endif
              IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
              IF ((dTdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
                dTdx=grad(i  ,Jend)
              ELSE
                dTdx=grad(i+1,Jend)
              END IF
              cff=MAX(dTdx*dTdx+dTde*dTde,eps)
#   ifdef RADIATION_2D
              Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
#   else
              Cx=0.0_r8
#   endif
              Ce=dTdt*dTde
#   if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%t_north_Cx(i,k,itrc)=Cx
              BOUNDARY(ng)%t_north_Ce(i,k,itrc)=Ce
              BOUNDARY(ng)%t_north_C2(i,k,itrc)=cff
#   endif
              t(i,Jend+1,k,nout,itrc)=(cff*t(i,Jend+1,k,nstp,itrc)+       &
     &                                 Ce *t(i,Jend  ,k,nout,itrc)-       &
     &                                 MAX(Cx,0.0_r8)*                    &
     &                                    grad(i  ,Jend+1)-               &
     &                                 MIN(Cx,0.0_r8)*                    &
     &                                    grad(i+1,Jend+1))/              &
     &                                (cff+Ce)
#   ifdef NORTH_PTNUDGING
              t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)+            &
     &                                tau*(BOUNDARY(ng)%t_north(i,k,itrc)-&
     &                                     t(i,Jend+1,k,nstp,itrc))
#   endif
#   ifdef MASKING
              t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
     &                                GRID(ng)%rmask(i,Jend+1)
#   endif
            END DO
          END DO

#  elif defined NORTH_PTCLAMPED
!
!  Northern edge, clamped boundary condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,Jend+1,k,nout,itrc)=BOUNDARY(ng)%t_north(i,k,itrc)
#   ifdef MASKING
              t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
     &                                GRID(ng)%rmask(i,Jend+1)
#   endif
            END DO
          END DO

#  elif defined NORTH_PTGRADIENT
!
!  Northern edge, gradient boundary condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,Jend+1,k,nout,itrc)=t(i,Jend,k,nout,itrc)
#   ifdef MASKING
              t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
     &                                GRID(ng)%rmask(i,Jend+1)
#   endif
            END DO
          END DO

#  else
!
!  Northern edge, closed boundary condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,Jend+1,k,nout,itrc)=t(i,Jend,k,nout,itrc)
#   ifdef MASKING
              t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
     &                                GRID(ng)%rmask(i,Jend+1)
#   endif
            END DO
          END DO
#  endif
        END IF
      END IF

!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng).and. .not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          DO k=1,N(ng)
            t(Istr-1,Jstr-1,k,nout,itrc)=0.5_r8*                        &
     &                                   (t(Istr  ,Jstr-1,k,nout,itrc)+ &
     &                                    t(Istr-1,Jstr  ,k,nout,itrc))
          END DO
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          DO k=1,N(ng)
            t(Iend+1,Jstr-1,k,nout,itrc)=0.5_r8*                        &
     &                                   (t(Iend  ,Jstr-1,k,nout,itrc)+ &
     &                                    t(Iend+1,Jstr  ,k,nout,itrc))
          END DO
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          DO k=1,N(ng)
            t(Istr-1,Jend+1,k,nout,itrc)=0.5_r8*                        &
     &                                   (t(Istr-1,Jend  ,k,nout,itrc)+ &
     &                                    t(Istr  ,Jend+1,k,nout,itrc))
          END DO
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          DO k=1,N(ng)
            t(Iend+1,Jend+1,k,nout,itrc)=0.5_r8*                        &
     &                                   (t(Iend+1,Jend  ,k,nout,itrc)+ &
     &                                    t(Iend  ,Jend+1,k,nout,itrc))
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE pt3dbc_tile
#endif
      END MODULE pt3dbc_mod
