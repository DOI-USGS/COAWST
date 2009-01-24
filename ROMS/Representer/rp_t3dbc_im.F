#include "cppdefs.h"
      MODULE rp_t3dbc_mod
#if defined TL_IOMS && defined SOLVE3D
!
!svn $Id: rp_t3dbc_im.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets representers tangent linear lateral boundary   !
!  conditions for the ITRC-th tracer field. It updates the specified   !
!  "nout" time index.                                                  !
!                                                                      !
!  BASIC STATE variables needed: t                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: rp_t3dbc, rp_t3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_t3dbc (ng, tile, nout, itrc)
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
      CALL rp_t3dbc_tile (ng, tile, itrc,                               &
     &                    LBi, UBi, LBj, UBj, N(ng), NT(ng),            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng), nout,                               &
     &                    OCEAN(ng)% tl_t)
      RETURN
      END SUBROUTINE rp_t3dbc
!
!***********************************************************************
      SUBROUTINE rp_t3dbc_tile (ng, tile, itrc,                         &
     &                          LBi, UBi, LBj, UBj, UBk, UBt,           &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp, nout,                             &
     &                          tl_t)
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
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
# else
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,UBk,3,UBt)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: Ce, Cx, cff, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_grad

# include "set_bounds.h"

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_TRADIATION

        IF (iic(ng).ne.0) THEN
!
!  Western edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
!>            grad(Istr-1,j)=t(Istr-1,j  ,k,nstp,itrc)-                 &
!>   &                       t(Istr-1,j-1,k,nstp,itrc)
!>
              tl_grad(Istr-1,j)=0.0_r8
            END DO
            DO j=Jstr,Jend
#   ifdef WEST_TNUDGING
              IF (BOUNDARY(ng)%t_west_Cx(j,k,itrc).eq.0.0_r8) THEN
                tau=Tobc_in(itrc,ng,iwest)
              ELSE
                tau=Tobc_out(itrc,ng,iwest)
              END IF
              tau=tau*dt(ng)
#   endif
              Cx=BOUNDARY(ng)%t_west_Cx(j,k,itrc)
#   ifdef RADIATION_2D
              Ce=BOUNDARY(ng)%t_west_Ce(j,k,itrc)
#   else
              Ce=0.0_r8
#   endif
              cff=BOUNDARY(ng)%t_west_C2(j,k,itrc)
!>            t(Istr-1,j,k,nout,itrc)=(cff*t(Istr-1,j,k,nstp,itrc)+     &
!>   &                                 Cx *t(Istr  ,j,k,nout,itrc)-     &
!>   &                                 MAX(Ce,0.0_r8)*                  &
!>   &                                    grad(Istr-1,j  )-             &
!>   &                                 MIN(Ce,0.0_r8)*                  &
!>   &                                    grad(Istr-1,j+1))/            &
!>   &                                (cff+Cx)
!>
              tl_t(Istr-1,j,k,nout,itrc)=(cff*                          &
     &                                    tl_t(Istr-1,j,k,nstp,itrc)+   &
     &                                    Cx *                          &
     &                                    tl_t(Istr  ,j,k,nout,itrc)-   &
     &                                    MAX(Ce,0.0_r8)*               &
     &                                       tl_grad(Istr-1,j  )-       &
     &                                    MIN(Ce,0.0_r8)*               &
     &                                       tl_grad(Istr-1,j+1))/      &
     &                                   (cff+Cx)
#   ifdef WEST_TNUDGING
!>            t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)+          &
!>   &                                tau*                              &
!>   &                                (BOUNDARY(ng)%t_west(j,k,itrc)-   &
!>   &                                 t(Istr-1,j,k,nstp,itrc))
!>
              tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,itrc)-    &
     &                                   tl_t(Istr-1,j,k,nstp,itrc)
#   endif
#   ifdef MASKING
!>            t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*          &
!>   &                                GRID(ng)%rmask(Istr-1,j)
!>
              tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,itrc)*    &
     &                                   GRID(ng)%rmask(Istr-1,j)
#   endif
            END DO
          END DO
        END IF

#  elif defined WEST_TCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          t(Istr-1,j,k,nout,itrc)=BOUNDARY(ng)%t_west(j,k,itrc)
#   ifdef MASKING
!>          t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(Istr-1,j)
#   endif
!>
            tl_t(Istr-1,j,k,nout,itrc)=BOUNDARY(ng)%t_west(j,k,itrc)
#   ifdef MASKING
            tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(Istr-1,j)
#   endif
          END DO
        END DO

#  elif defined WEST_TGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          t(Istr-1,j,k,nout,itrc)=t(Istr,j,k,nout,itrc)
!>
            tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr,j,k,nout,itrc)
#   ifdef MASKING
!>          t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(Istr-1,j)
!>
            tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(Istr-1,j)
#   endif
          END DO
        END DO

#  else
!
!  Western edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          t(Istr-1,j,k,nout,itrc)=t(Istr,j,k,nout,itrc)
!>
            tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr,j,k,nout,itrc)
#   ifdef MASKING
!>          t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(Istr-1,j)
!>
            tl_t(Istr-1,j,k,nout,itrc)=tl_t(Istr-1,j,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(Istr-1,j)
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

        IF (iic(ng).ne.0) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
!>            grad(Iend+1,j)=t(Iend+1,j  ,k,nstp,itrc)-                 &
!>   &                       t(Iend+1,j-1,k,nstp,itrc)
!>
              tl_grad(Iend+1,j)=0.0_r8
            END DO
            DO j=Jstr,Jend
#   ifdef EAST_TNUDGING
              IF (BOUNDARY(ng)%t_east_Cx(j,k,itrc).eq.0.0_r8) THEN
                tau=Tobc_in(itrc,ng,ieast)
              ELSE
                tau=Tobc_out(itrc,ng,ieast)
              END IF
              tau=tau*dt(ng)
#   endif
              Cx=BOUNDARY(ng)%t_east_Cx(j,k,itrc)
#   ifdef RADIATION_2D
              Ce=BOUNDARY(ng)%t_east_Ce(j,k,itrc)
#   else
              Ce=0.0_r8
#   endif
              cff=BOUNDARY(ng)%t_east_C2(j,k,itrc)
!>            t(Iend+1,j,k,nout,itrc)=(cff*t(Iend+1,j,k,nstp,itrc)+     &
!>   &                                 Cx *t(Iend  ,j,k,nout,itrc)-     &
!>   &                                 MAX(Ce,0.0_r8)*                  &
!>   &                                    grad(Iend+1,j  )-             &
!>   &                                 MIN(Ce,0.0_r8)*                  &
!>   &                                    grad(Iend+1,j+1))/            &
!>   &                                (cff+Cx)
!>
              tl_t(Iend+1,j,k,nout,itrc)=(cff*                          &
     &                                    tl_t(Iend+1,j,k,nstp,itrc)+   &
     &                                    Cx *                          &
     &                                    tl_t(Iend  ,j,k,nout,itrc)-   &
     &                                    MAX(Ce,0.0_r8)*               &
     &                                       tl_grad(Iend+1,j  )-       &
     &                                    MIN(Ce,0.0_r8)*               &
     &                                       tl_grad(Iend+1,j+1))/      &
     &                                   (cff+Cx)
#   ifdef EAST_TNUDGING
!>            t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)+          &
!>   &                                tau*                              &
!>   &                                (BOUNDARY(ng)%t_east(j,k,itrc)-   &
!>   &                                 t(Iend+1,j,k,nstp,itrc))
!>
              tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,itrc)-    &
     &                                   tau*tl_t(Iend+1,j,k,nstp,itrc)
#   endif
#   ifdef MASKING
!>            t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*          &
!>   &                                GRID(ng)%rmask(Iend+1,j)
!>
              tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,itrc)*    &
     &                                   GRID(ng)%rmask(Iend+1,j)
#   endif
            END DO
          END DO
        END IF

#  elif defined EAST_TCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          t(Iend+1,j,k,nout,itrc)=BOUNDARY(ng)%t_east(j,k,itrc)
#   ifdef MASKING
!>          t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(Iend+1,j)
#   endif
!>
            tl_t(Iend+1,j,k,nout,itrc)=BOUNDARY(ng)%t_east(j,k,itrc)
#   ifdef MASKING
            tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO

#  elif defined EAST_TGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          t(Iend+1,j,k,nout,itrc)=t(Iend,j,k,nout,itrc)
!>
            tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend,j,k,nout,itrc)
#   ifdef MASKING
!>          t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(Iend+1,j)
!>
            tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(Iend+1,j)
#   endif
          END DO
        END DO

#  else
!
!  Eastern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
!>          t(Iend+1,j,k,nout,itrc)=t(Iend,j,k,nout,itrc)
!>
            tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend,j,k,nout,itrc)
#   ifdef MASKING
!>          t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(Iend+1,j)
!>
            tl_t(Iend+1,j,k,nout,itrc)=tl_t(Iend+1,j,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(Iend+1,j)
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

        IF (iic(ng).ne.0) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend+1
!>            grad(i,Jstr-1)=t(i  ,Jstr-1,k,nstp,itrc)-                 &
!>   &                       t(i-1,Jstr-1,k,nstp,itrc)
!>
              tl_grad(i,Jstr-1)=0.0_r8
            END DO
            DO i=Istr,Iend
#   ifdef SOUTH_TNUDGING
              IF (BOUNDARY(ng)%t_south_Ce(i,k,itrc).eq.0.0_r8) THEN
                tau=Tobc_in(itrc,ng,isouth)
              ELSE
                tau=Tobc_out(itrc,ng,isouth)
              END IF
              tau=tau*dt(ng)
#   endif
#   ifdef RADIATION_2D
              Cx=BOUNDARY(ng)%t_south_Cx(i,k,itrc)
#   else
              Cx=0.0_r8
#   endif
              Ce=BOUNDARY(ng)%t_south_Ce(i,k,itrc)
              cff=BOUNDARY(ng)%t_south_C2(i,k,itrc)
!>            t(i,Jstr-1,k,nout,itrc)=(cff*t(i,Jstr-1,k,nstp,itrc)+     &
!>   &                                 Ce *t(i,Jstr  ,k,nout,itrc )-    &
!>   &                                 MAX(Cx,0.0_r8)*                  &
!>   &                                    grad(i  ,Jstr-1)-             &
!>   &                                 MIN(Cx,0.0_r8)*                  &
!>   &                                    grad(i+1,Jstr-1))/            &
!>   &                                (cff+Ce)
!>
              tl_t(i,Jstr-1,k,nout,itrc)=(cff*                          &
     &                                    tl_t(i,Jstr-1,k,nstp,itrc)+   &
     &                                    Ce *                          &
     &                                    tl_t(i,Jstr  ,k,nout,itrc)-   &
     &                                    MAX(Cx,0.0_r8)*               &
     &                                       tl_grad(i  ,Jstr-1)-       &
     &                                    MIN(Cx,0.0_r8)*               &
     &                                       tl_grad(i+1,Jstr-1))/      &
     &                                   (cff+Ce)
#   ifdef SOUTH_TNUDGING
!>            t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)+          &
!>   &                                tau*                              &
!>   &                                (BOUNDARY(ng)%t_south(i,k,itrc)-  &
!>   &                                 t(i,Jstr-1,k,nstp,itrc))
!>
              tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,itrc)-    &
     &                                   tau*tl_t(i,Jstr-1,k,nstp,itrc)
#   endif
#   ifdef MASKING
!>            t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*          &
!>   &                                GRID(ng)%rmask(i,Jstr-1)
!>
              tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,itrc)*    &
     &                                   GRID(ng)%rmask(i,Jstr-1)
#   endif
            END DO
          END DO
        END IF

#  elif defined SOUTH_TCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          t(i,Jstr-1,k,nout,itrc)=BOUNDARY(ng)%t_south(i,k,itrc)
#   ifdef MASKING
!>          t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(i,Jstr-1)
#   endif
!>
            tl_t(i,Jstr-1,k,nout,itrc)=BOUNDARY(ng)%t_south(i,k,itrc)
#   ifdef MASKING
            tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(i,Jstr-1)
#   endif
          END DO
        END DO

#  elif defined SOUTH_TGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          t(i,Jstr-1,k,nout,itrc)=t(i,Jstr,k,nout,itrc)
!>
            tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr,k,nout,itrc)
#   ifdef MASKING
!>          t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(i,Jstr-1)
!>
            tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(i,Jstr-1)
#   endif
          END DO
        END DO

#  else
!
!  Southern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          t(i,Jstr-1,k,nout,itrc)=t(i,Jstr,k,nout,itrc)
!>
            tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr,k,nout,itrc)
#   ifdef MASKING
!>          t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(i,Jstr-1)
!>
            tl_t(i,Jstr-1,k,nout,itrc)=tl_t(i,Jstr-1,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(i,Jstr-1)
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

        IF (iic(ng).ne.0) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
          DO k=1,N(ng)
            DO i=Istr,Iend+1
!>            grad(i,Jend+1)=t(i  ,Jend+1,k,nstp,itrc)-                 &
!>   &                       t(i-1,Jend+1,k,nstp,itrc)
!>
              tl_grad(i,Jend+1)=0.0_r8
            END DO
            DO i=Istr,Iend
#   ifdef NORTH_TNUDGING
              IF (BOUNDARY(ng)%t_north_Ce(i,k,itrc).eq.0.0_r8) THEN
                tau=Tobc_in(itrc,ng,inorth)
              ELSE
                tau=Tobc_out(itrc,ng,inorth)
              END IF
              tau=tau*dt(ng)
#   endif
#   ifdef RADIATION_2D
              Cx=BOUNDARY(ng)%t_north_Cx(i,k,itrc)
#   else
              Cx=0.0_r8
#   endif
              Ce=BOUNDARY(ng)%t_north_Ce(i,k,itrc)
              cff=BOUNDARY(ng)%t_north_C2(i,k,itrc)
!>            t(i,Jend+1,k,nout,itrc)=(cff*t(i,Jend+1,k,nstp,itrc)+     &
!>   &                                 Ce *t(i,Jend  ,k,nout,itrc)-     &
!>   &                                 MAX(Cx,0.0_r8)*                  &
!>   &                                    grad(i  ,Jend+1)-             &
!>   &                                 MIN(Cx,0.0_r8)*                  &
!>   &                                    grad(i+1,Jend+1))/            &
!>   &                                (cff+Ce)
!>
              tl_t(i,Jend+1,k,nout,itrc)=(cff*                          &
     &                                    tl_t(i,Jend+1,k,nstp,itrc)+   &
     &                                    Ce *                          &
     &                                    tl_t(i,Jend  ,k,nout,itrc)-   &
     &                                    MAX(Cx,0.0_r8)*               &
     &                                       tl_grad(i  ,Jend+1)-       &
     &                                    MIN(Cx,0.0_r8)*               &
     &                                       tl_grad(i+1,Jend+1))/      &
     &                                   (cff+Ce)
#   ifdef NORTH_TNUDGING
!>            t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)+          &
!>   &                                tau*                              & 
!>   &                                (BOUNDARY(ng)%t_north(i,k,itrc)-  &
!>   &                                 t(i,Jend+1,k,nstp,itrc))
!>
              tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,itrc)-    &
     &                                   tau*tl_t(i,Jend+1,k,nstp,itrc)
#   endif
#   ifdef MASKING
!>            t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*          &
!>   &                                GRID(ng)%rmask(i,Jend+1)
!>
              tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,itrc)*    &
     &                                   GRID(ng)%rmask(i,Jend+1)
#   endif
            END DO
          END DO
        END IF

#  elif defined NORTH_TCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          t(i,Jend+1,k,nout,itrc)=BOUNDARY(ng)%t_north(i,k,itrc)
#   ifdef MASKING
!>          t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(i,Jend+1)
#   endif
!>
            tl_t(i,Jend+1,k,nout,itrc)=BOUNDARY(ng)%t_north(i,k,itrc)
#   ifdef MASKING
            tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO

#  elif defined NORTH_TGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          t(i,Jend+1,k,nout,itrc)=t(i,Jend,k,nout,itrc)
!>
            tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend,k,nout,itrc)
#   ifdef MASKING
!>          t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(i,Jend+1)
!>
            tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(i,Jend+1)
#   endif
          END DO
        END DO

#  else
!
!  Northern edge, closed boundary condition.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          t(i,Jend+1,k,nout,itrc)=t(i,Jend,k,nout,itrc)
!>
            tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend,k,nout,itrc)
#   ifdef MASKING
!>          t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)*            &
!>   &                              GRID(ng)%rmask(i,Jend+1)
!>
            tl_t(i,Jend+1,k,nout,itrc)=tl_t(i,Jend+1,k,nout,itrc)*      &
     &                                 GRID(ng)%rmask(i,Jend+1)
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
!>        t(Istr-1,Jstr-1,k,nout,itrc)=0.5_r8*                          &
!>   &                                 (t(Istr  ,Jstr-1,k,nout,itrc)+   &
!>   &                                  t(Istr-1,Jstr  ,k,nout,itrc))
!>
          tl_t(Istr-1,Jstr-1,k,nout,itrc)=0.5_r8*                       &
     &                              (tl_t(Istr  ,Jstr-1,k,nout,itrc)+   &
     &                               tl_t(Istr-1,Jstr  ,k,nout,itrc))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        t(Iend+1,Jstr-1,k,nout,itrc)=0.5_r8*                          &
!>   &                                 (t(Iend  ,Jstr-1,k,nout,itrc)+   &
!>   &                                  t(Iend+1,Jstr  ,k,nout,itrc))
!>
          tl_t(Iend+1,Jstr-1,k,nout,itrc)=0.5_r8*                       &
     &                              (tl_t(Iend  ,Jstr-1,k,nout,itrc)+   &
     &                               tl_t(Iend+1,Jstr  ,k,nout,itrc))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        t(Istr-1,Jend+1,k,nout,itrc)=0.5_r8*                          &
!>   &                                 (t(Istr-1,Jend  ,k,nout,itrc)+   &
!>   &                                  t(Istr  ,Jend+1,k,nout,itrc))
!>
          tl_t(Istr-1,Jend+1,k,nout,itrc)=0.5_r8*                       &
     &                              (tl_t(Istr-1,Jend  ,k,nout,itrc)+   &
     &                               tl_t(Istr  ,Jend+1,k,nout,itrc))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO k=1,N(ng)
!>        t(Iend+1,Jend+1,k,nout,itrc)=0.5_r8*                          &
!>   &                                 (t(Iend+1,Jend  ,k,nout,itrc)+   &
!>   &                                  t(Iend  ,Jend+1,k,nout,itrc))
!>
          tl_t(Iend+1,Jend+1,k,nout,itrc)=0.5_r8*                       &
     &                              (tl_t(Iend+1,Jend  ,k,nout,itrc)+   &
     &                               tl_t(Iend  ,Jend+1,k,nout,itrc))
        END DO
      END IF
# endif
      RETURN
      END SUBROUTINE rp_t3dbc_tile
#endif
      END MODULE rp_t3dbc_mod
