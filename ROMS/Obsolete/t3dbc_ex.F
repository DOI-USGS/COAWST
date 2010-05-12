#include "cppdefs.h"
      MODULE t3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: t3dbc_ex.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
      PUBLIC  :: t3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE t3dbc (ng, tile, nout, itrc)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean

      integer, intent(in) :: ng, tile, nout, itrc

# include "tile.h"
!
      CALL t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,                &
     &                 LBi, UBi, LBj, UBj, N(ng), NT(ng),               &
     &                 nout, nstp(ng),                                  &
     &                 OCEAN(ng)% t)
      RETURN
      END SUBROUTINE t3dbc
!
!***********************************************************************
      SUBROUTINE t3dbc_tile (ng, Istr, Iend, Jstr, Jend, itrc,          &
     &                       LBi, UBi, LBj, UBj, UBk, UBt,              &
     &                       nout, nstp,                                &
     &                       t)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod mask
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Iend, Istr, Jend, Jstr, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: nout, nstp

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
# endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, k

      real(r8), parameter :: eps =1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, eps, tau

      real(r80, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_TRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
            grad(Istr-1,j)=t(Istr-1,j  ,k,nstp,itrc)-                   &
     &                     t(Istr-1,j-1,k,nstp,itrc)
#   ifdef MASKING
            grad(Istr-1,j)=grad(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#   endif
            grad(Istr  ,j)=t(Istr  ,j  ,k,nstp,itrc)-                   &
     &                     t(Istr  ,j-1,k,nstp,itrc)
#   ifdef MASKING
            grad(Istr  ,j)=grad(Istr  ,j)*GRID(ng)%vmask(Istr  ,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dTdt=t(Istr,j,k,nstp,itrc)-t(Istr  ,j,k,nout,itrc)
            dTdx=t(Istr,j,k,nstp,itrc)-t(Istr+1,j,k,nstp,itrc)
#   ifdef WEST_TNUDGING
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
            cff=dTdt/MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=MIN(1.0_r8,cff*dTdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dTde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            t(Istr-1,j,k,nout,itrc)=(1.0_r8-Cx)*t(Istr-1,j,k,nstp,itrc)+&
     &                              Cx*t(Istr,j,k,nstp,itrc )-          &
     &                              MAX(Ce,0.0_r8)*grad(Istr-1,j  )-    &
     &                              MIN(Ce,0.0_r8)*grad(Istr-1,j+1)
#   ifdef WEST_TNUDGING
     &      t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)+            &
     &                              tau*(BOUNDARY(ng)%t_west(j,k,itrc)- &
     &                                   t(Istr-1,j,k,nstp,itrc))
#   endif
#   ifdef MASKING
            t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_TCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(Istr-1,j,k,nout,itrc)=BOUNDARY(ng)%t_west(j,k,itrc)
#   ifdef MASKING
            t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_TGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(Istr-1,j,k,nout,itrc)=t(Istr,j,k,nout,itrc)
#   ifdef MASKING
            t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(Istr-1,j)
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
            t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(0,j)
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

#  if defined EAST_TRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend+1
           grad(Iend  ,j)=t(Iend  ,j  ,k,nstp,itrc)-                    &
     &                    t(Iend  ,j-1,k,nstp,itrc)
#   ifdef MASKING
           grad(Iend  ,j)=grad(Iend  ,j)*GRID(ng)%vmask(Iend  ,j)
#   endif
           grad(Iend+1,j)=t(Iend+1,j  ,k,nstp,itrc)-                    &
     &                    t(Iend+1,j-1,k,nstp,itrc)
#   ifdef MASKING
           grad(Iend+1,j)=grad(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
          DO j=Jstr,Jend
            dTdt=t(Iend,j,k,nstp,itrc)-t(Iend  ,j,k,nout,itrc)
            dTdx=t(Iend,j,k,nstp,itrc)-t(Iend-1,j,k,nstp,itrc)
#   ifdef EAST_TNUDGING
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
            cff=dTdt/MAX(dTdx*dTdx+dTde*dTde,eps)
            Cx=MIN(1.0_r8,cff*dTdx)
#   ifdef RADIATION_2D
            Ce=MIN(1.0_r8,MAX(cff*dTde,-1.0_r8))
#   else
            Ce=0.0_r8
#   endif
            t(Iend+1,j,k,nout,itrc)=(1.0_r8-Cx)*t(Iend+1,j,k,nstp,itrc)+&
     &                              Cx*t(Iend,j,k,nstp,itrc)-           &
     &                              MAX(Ce,0.0_r8)*grad(Iend+1,j  )-    &
     &                              MIN(Ce,0.0_r8)*grad(Iend+1,j+1)
#   ifdef EAST_TNUDGING
     &      t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)+            &
     &                              tau*(BOUNDARY(ng)%t_east(j,k,itrc)- &
     &                                   t(Iend+1,j,k,nstp,itrc))
#   endif
#   ifdef MASKING
            t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_TCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(Iend+1,j,k,nout,itrc)=BOUNDARY(ng)%t_east(j,k,itrc)
#   ifdef MASKING
            t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_TGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            t(Iend+1,j,k,nout,itrc)=t(Iend,j,k,nout,itrc)
#   ifdef MASKING
            t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*            &
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
            t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO
#  endif
      END IF
# endif

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_TRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend+1
            grad(i,Jstr  )=t(i  ,Jstr  ,k,nstp,itrc)-                   &
     &                     t(i-1,Jstr  ,k,nstp,itrc)
#   ifdef MASKING
            grad(i,Jstr  )=grad(i,Jstr  )*GRID(ng)%umask(i,Jstr  )
#   endif
            grad(i,Jstr-1)=t(i  ,Jstr-1,k,nstp,itrc)-                   &
     &                     t(i-1,Jstr-1,k,nstp,itrc)
#   ifdef MASKING
            grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#   endif
          END DO
          DO i=Istr,Iend
            dTdt=t(i,Jstr,k,nstp,itrc)-t(i,Jstr  ,k,nout,itrc)
            dTde=t(i,Jstr,k,nstp,itrc)-t(i,Jstr+1,k,nstp,itrc)
#   ifdef SOUTH_TNUDGING
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
            cff=dTdt/MAX(dTdx*dTdx+dTde*dTde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dTdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dTde)
            t(i,Jstr-1,k,nout,itrc)=(1.0_r8-Ce)*t(i,Jstr-1,k,nstp,itrc)+&
     &                              Ce*t(i,Jstr,k,nstp,itrc )-          &
     &                              MAX(Cx,0.0_r8)*grad(i  ,Jstr-1)-    &
     &                              MIN(Cx,0.0_r8)*grad(i+1,Jstr-1)
#   ifdef SOUTH_TNUDGING
            t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)+            &
     &                              tau*(BOUNDARY(ng)%t_south(i,k,itrc)-&
     &                                   t(i,Jstr-1,k,nstp,itrc))
#   endif
#   ifdef MASKING
            t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(i,0)
#   endif
          END DO
        END DO

#  elif defined SOUTH_TCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,Jstr-1,k,nout,itrc)=BOUNDARY(ng)%t_south(i,k,itrc)
#   ifdef MASKING
            t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(i,0)
#   endif
          END DO
        END DO

#  elif defined SOUTH_TGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,Jstr-1,k,nout,itrc)=t(i,Jstr,k,nout,itrc)
#   ifdef MASKING
            t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(i,0)
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
     &                              GRID(ng)%rmask(i,0)
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
      IF (NORTHERN_EDGE) THEN

#  if defined NORTH_TRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend+1
            grad(i,Jend  )=t(i  ,Jend  ,k,nstp,itrc)-                   &
     &                     t(i-1,Jend  ,k,nstp,itrc)
#   ifdef MASKING
            grad(i,Jend  )=grad(i,Jend  )*GRID(ng)%umask(i,Jend  )
#   endif
            grad(i,Jend+1)=t(i  ,Jend+1,k,nstp,itrc)-                   &
     &                     t(i-1,Jend+1,k,nstp,itrc)
#   ifdef MASKING
            grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#   endif
          END DO
          DO i=Istr,Iend
            dTdt=t(i,Jend,k,nstp,itrc)-t(i,Jend  ,k,nout,itrc)
            dTde=t(i,Jend,k,nstp,itrc)-t(i,Jend-1,k,nstp,itrc)
#   ifdef NORTH_TNUDGING
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
            cff=dTdt/MAX(dTdx*dTdx+dTde*dTde,eps)
#   ifdef RADIATION_2D
            Cx=MIN(1.0_r8,MAX(cff*dTdx,-1.0_r8))
#   else
            Cx=0.0_r8
#   endif
            Ce=MIN(1.0_r8,cff*dTde)
            t(i,Jend+1,k,nout,itrc)=(1.0_r8-Ce)*t(i,Jend+1,k,nstp,itrc)+&
     &                              Ce*t(i,Jend,k,nstp,itrc)-           &
     &                              MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-    &
     &                              MIN(Cx,0.0_r8)*grad(i+1,Jend+1)
#   ifdef NORTH_TNUDGING
            t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)+            &
     &                              tau*(BOUNDARY(ng)%t_north(i,k,itrc)-&
     &                                   t(i,Jend+1,k,nstp,itrc))
#   endif
#   ifdef MASKING
            t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_TCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,Jend+1,k,nout,itrc)=BOUNDARY(ng)%t_north(i,k,itrc)
#   ifdef MASKING
            t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_TGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            t(i,Jend+1,k,nout,itrc)=t(i,Jend,k,nout,itrc)
#   ifdef MASKING
            t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
     &                              GRID(ng)%rmask(i,Jend+1)
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
     &                              GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO
#  endif
      END IF
# endif

# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
          t(Istr-1,Jstr-1,k,nout,itrc)=0.5_r8*                          &
     &                                 (t(Istr  ,Jstr-1,k,nout,itrc)+   &
     &                                  t(Istr-1,Jstr  ,k,nout,itrc))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
          t(Iend+1,Jstr-1,k,nout,itrc)=0.5_r8*                          &
     &                                 (t(Iend  ,Jstr-1,k,nout,itrc)+   &
     &                                  t(Iend+1,Jstr  ,k,nout,itrc))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
          t(Istr-1,Jend+1,k,nout,itrc)=0.5_r8*                          &
     &                                 (t(Istr  ,Jend+1,k,nout,itrc)+   &
     &                                  t(Istr-1,Jend  ,k,nout,itrc))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
           t(Iend+1,Jend+1,k,nout,itrc)=0.5_r8*                         &
     &                                  (t(Iend  ,Jend+1,k,nout,itrc)+  &
     &                                   t(Iend+1,Jend  ,k,nout,itrc))
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE t3dbc_tile
#endif
      END MODULE t3dbc_mod
