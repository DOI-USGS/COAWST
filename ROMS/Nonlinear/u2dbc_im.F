#include "cppdefs.h"
      MODULE u2dbc_mod
!
!svn $Id: u2dbc_im.F 779 2008-10-04 23:15:47Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated U-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: u2dbc, u2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE u2dbc (ng, tile, kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kout
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL u2dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 krhs(ng), kstp(ng), kout,                        &
     &                 OCEAN(ng) % ubar,                                &
     &                 OCEAN(ng) % vbar,                                &
     &                 OCEAN(ng) % zeta)
      RETURN
      END SUBROUTINE u2dbc
!
!***********************************************************************
      SUBROUTINE u2dbc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       krhs, kstp, kout,                          &
     &                       ubar, vbar, zeta)
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
      USE mod_ocean
      USE mod_stepping
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: krhs, kstp, kout

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)

      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, know

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx, cff5
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8) :: cff, cff1, cff2, dt2d, dUde, dUdt, dUdx, tau

#ifdef REFINED_GRID
      integer :: t1, t2, ngp
      real(r8) :: fac, fac1, fac2, my_val
#endif

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
      IF (FIRST_2D_STEP) THEN
        know=krhs
        dt2d=dtfast(ng)
      ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
        know=krhs
        dt2d=2.0_r8*dtfast(ng)
      ELSE
        know=kstp
        dt2d=dtfast(ng)
      END IF

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      IF (ng.eq.1) THEN
#endif

#if !defined EW_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

# if defined WEST_M2RADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Istr  ,j)=ubar(Istr  ,j  ,know)-                         &
     &                   ubar(Istr  ,j-1,know)
          grad(Istr+1,j)=ubar(Istr+1,j  ,know)-                         &
     &                   ubar(Istr+1,j-1,know)
        END DO
        DO j=Jstr,Jend
          dUdt=ubar(Istr+1,j,know)-ubar(Istr+1,j,kout)
          dUdx=ubar(Istr+1,j,kout)-ubar(Istr+2,j,kout)
#  ifdef WEST_M2NUDGING
          IF ((dUdt*dUdx).lt.0.0_r8) THEN
            tau=M2obc_in(ng,iwest)
          ELSE
            tau=M2obc_out(ng,iwest)
          END IF
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.0.0_r8) THEN
            dUde=grad(Istr+1,j  )
          ELSE
            dUde=grad(Istr+1,j+1)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=dUdt*dUdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#  else
          Ce=0.0_r8
#  endif
#  if defined CELERITY_WRITE && defined FORWARD_WRITE
          BOUNDARY(ng)%ubar_west_Cx(j)=Cx
          BOUNDARY(ng)%ubar_west_Ce(j)=Ce
          BOUNDARY(ng)%ubar_west_C2(j)=cff
#  endif
          ubar(Istr,j,kout)=(cff*ubar(Istr  ,j,know)+                   &
     &                       Cx *ubar(Istr+1,j,kout)-                   &
     &                       MAX(Ce,0.0_r8)*grad(Istr,j  )-             &
     &                       MIN(Ce,0.0_r8)*grad(Istr,j+1))/            &
     &                      (cff+Cx)
#  ifdef WEST_M2NUDGING
         ubar(Istr,j,kout)=ubar(Istr,j,kout)+                           &
     &                     tau*(BOUNDARY(ng)%ubar_west(j)-              &
     &                          ubar(Istr,j,know))
#  endif
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# elif defined WEST_M2FLATHER
!
!  Western edge, Flather boundary condition.
!
        DO j=Jstr,Jend
#  if defined SSH_TIDES && !defined UV_TIDES
#   ifdef FSOBC_REDUCED
          bry_pgr=-g*(zeta(Istr,j,know)-                                &
     &                BOUNDARY(ng)%zeta_west(j))*                       &
     &            0.5_r8*GRID(ng)%pm(Istr,j)
#   else
          bry_pgr=-g*(zeta(Istr  ,j,know)-                              &
     &                zeta(Istr-1,j,know))*                             &
     &            0.5_r8*(GRID(ng)%pm(Istr-1,j)+                        &
     &                    GRID(ng)%pm(Istr  ,j))
#   endif
#   ifdef UV_COR
          bry_cor=0.125_r8*(vbar(Istr-1,j  ,know)+                      &
     &                      vbar(Istr-1,j+1,know)+                      &
     &                      vbar(Istr  ,j  ,know)+                      &
     &                      vbar(Istr  ,j+1,know))*                     &
     &                     (GRID(ng)%f(Istr-1,j)+                       &
     &                      GRID(ng)%f(Istr  ,j))
#   else
          bry_cor=0.0_r8
#   endif
          cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                    &
     &                         zeta(Istr-1,j,know)+                     &
     &                         GRID(ng)%h(Istr  ,j)+                    &
     &                         zeta(Istr  ,j,know)))
          bry_str=cff1*(FORCES(ng)%sustr(Istr,j)-                       &
     &                  FORCES(ng)%bustr(Istr,j))
          Cx=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(Istr-1,j)+                &
     &                             zeta(Istr-1,j,know)+                 &
     &                             GRID(ng)%h(Istr  ,j)+                &
     &                             zeta(Istr  ,j,know)))
          cff2=GRID(ng)%om_u(Istr,j)*Cx
!!        cff2=dt2d
          bry_val=ubar(Istr+1,j,know)+                                  &
     &            cff2*(bry_pgr+                                        &
     &                  bry_cor+                                        &
     &                  bry_str)
#  else
          bry_val=BOUNDARY(ng)%ubar_west(j)
#  endif
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                     &
     &                        zeta(Istr-1,j,know)+                      &
     &                        GRID(ng)%h(Istr  ,j)+                     &
     &                        zeta(Istr  ,j,know)))
          Cx=SQRT(g*cff)
          ubar(Istr,j,kout)=bry_val-                                    &
     &                      Cx*(0.5_r8*(zeta(Istr-1,j,know)+            &
     &                                  zeta(Istr  ,j,know))-           &
     &                          BOUNDARY(ng)%zeta_west(j))
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# elif defined WEST_M2CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ubar(Istr,j,kout)=BOUNDARY(ng)%ubar_west(j)
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# elif defined WEST_M2GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ubar(Istr,j,kout)=ubar(Istr+1,j,kout)
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# elif defined WEST_M2REDUCED
!
!  Western edge, reduced-physics boundary condition.
!
        DO j=Jstr,Jend
#  ifdef FSOBC_REDUCED
          bry_pgr=-g*(zeta(Istr,j,know)-                                &
     &                BOUNDARY(ng)%zeta_west(j))*                       &
     &            0.5_r8*GRID(ng)%pm(Istr,j)
#  else
          bry_pgr=-g*(zeta(Istr  ,j,know)-                              &
     &                zeta(Istr-1,j,know))*                             &
     &            0.5_r8*(GRID(ng)%pm(Istr-1,j)+                        &
     &                    GRID(ng)%pm(Istr  ,j))
#  endif
#  ifdef UV_COR
          bry_cor=0.125_r8*(vbar(Istr-1,j  ,know)+                      &
     &                      vbar(Istr-1,j+1,know)+                      &
     &                      vbar(Istr  ,j  ,know)+                      &
     &                      vbar(Istr  ,j+1,know))*                     &
     &                     (GRID(ng)%f(Istr-1,j)+                       &
     &                      GRID(ng)%f(Istr  ,j))
#  else
          bry_cor=0.0_r8
#  endif
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                     &
     &                        zeta(Istr-1,j,know)+                      &
     &                        GRID(ng)%h(Istr  ,j)+                     &
     &                        zeta(Istr  ,j,know)))
          bry_str=cff*(FORCES(ng)%sustr(Istr,j)-                        &
     &                 FORCES(ng)%bustr(Istr,j))
          ubar(Istr,j,kout)=ubar(Istr,j,know)+                          &
     &                      dt2d*(bry_pgr+                              &
     &                            bry_cor+                              &
     &                            bry_str)
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ubar(Istr,j,kout)=0.0_r8
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN

# if defined EAST_M2RADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Iend  ,j)=ubar(Iend  ,j  ,know)-                         &
     &                   ubar(Iend  ,j-1,know)
          grad(Iend+1,j)=ubar(Iend+1,j  ,know)-                         &
     &                   ubar(Iend+1,j-1,know)
        END DO
        DO j=Jstr,Jend
          dUdt=ubar(Iend,j,know)-ubar(Iend  ,j,kout)
          dUdx=ubar(Iend,j,kout)-ubar(Iend-1,j,kout)
#  ifdef EAST_M2NUDGING
          IF ((dUdt*dUdx).lt.0.0_r8) THEN
            tau=M2obc_in(ng,ieast)
          ELSE
            tau=M2obc_out(ng,ieast)
          END IF
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
            dUde=grad(Iend,j)
          ELSE
            dUde=grad(Iend,j+1)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=dUdt*dUdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#  else
          Ce=0.0_r8
#  endif
#  if defined CELERITY_WRITE && defined FORWARD_WRITE
          BOUNDARY(ng)%ubar_east_Cx(j)=Cx
          BOUNDARY(ng)%ubar_east_Ce(j)=Ce
          BOUNDARY(ng)%ubar_east_C2(j)=cff
#  endif
          ubar(Iend+1,j,kout)=(cff*ubar(Iend+1,j,know)+                 &
     &                         Cx *ubar(Iend  ,j,kout)-                 &
     &                         MAX(Ce,0.0_r8)*grad(Iend+1,j  )-         &
     &                         MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/        &
     &                        (cff+Cx)
#  ifdef EAST_M2NUDGING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)+                      &
     &                        tau*(BOUNDARY(ng)%ubar_east(j)-           &
     &                             ubar(Iend+1,j,know))
#  endif
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2FLATHER
!
!  Eastern edge, Flather boundary condition.
!
        DO j=Jstr,Jend
#  if defined SSH_TIDES && !defined UV_TIDES
#   ifdef FSOBC_REDUCED
          bry_pgr=-g*(BOUNDARY(ng)%zeta_east(j)-                        &
     &                zeta(Iend,j,know))*                               &
     &            0.5_r8*GRID(ng)%pm(Iend,j)
#   else
          bry_pgr=-g*(zeta(Iend+1,j,know)-                              &
     &                zeta(Iend  ,j,know))*                             &
     &            0.5_r8*(GRID(ng)%pm(Iend  ,j)+                        &
     &                    GRID(ng)%pm(Iend+1,j))
#   endif
#   ifdef UV_COR
          bry_cor=0.125_r8*(vbar(Iend  ,j  ,know)+                      &
     &                      vbar(Iend  ,j+1,know)+                      &
     &                      vbar(Iend+1,j  ,know)+                      &
     &                      vbar(Iend+1,j+1,know))*                     &
     &                     (GRID(ng)%f(Iend  ,j)+                       &
     &                      GRID(ng)%f(Iend+1,j))
#   else
          bry_cor=0.0_r8
#   endif
          cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                    &
     &                         zeta(Iend  ,j,know)+                     &
     &                         GRID(ng)%h(Iend+1,j)+                    &
     &                         zeta(Iend+1,j,know)))
          bry_str=cff1*(FORCES(ng)%sustr(Iend+1,j)-                     &
     &                  FORCES(ng)%bustr(Iend+1,j))
          Cx=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(Iend+1,j)+                &
     &                             zeta(Iend+1,j,know)+                 &
     &                             GRID(ng)%h(Iend  ,j)+                &
     &                             zeta(Iend  ,j,know)))
          cff2=GRID(ng)%om_u(Iend+1,j)*Cx
!!        cff2=dt2d
          bry_val=ubar(Iend,j,know)+                                    &
     &            cff2*(bry_pgr+                                        &
     &                  bry_cor+                                        &
     &                  bry_str)
#  else
          bry_val=BOUNDARY(ng)%ubar_east(j)
#  endif
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                     &
     &                        zeta(Iend  ,j,know)+                      &
     &                        GRID(ng)%h(Iend+1,j)+                     &
     &                        zeta(Iend+1,j,know)))
          Cx=SQRT(g*cff)
          ubar(Iend+1,j,kout)=bry_val+                                  &
     &                        Cx*(0.5_r8*(zeta(Iend  ,j,know)+          &
     &                                    zeta(Iend+1,j,know))-         &
     &                            BOUNDARY(ng)%zeta_east(j))
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ubar(Iend+1,j,kout)=BOUNDARY(ng)%ubar_east(j)
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ubar(Iend+1,j,kout)=ubar(Iend,j,kout)
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_M2REDUCED
!
!  Eastern edge, reduced-physics boundary condition.
!
        DO j=Jstr,Jend
#  ifdef FSOBC_REDUCED
          bry_pgr=-g*(BOUNDARY(ng)%zeta_east(j)-                        &
     &                zeta(Iend,j,know))*                               &
     &             0.5_r8*GRID(ng)%pm(Iend,j)
#  else
          bry_pgr=-g*(zeta(Iend+1,j,know)-                              &
     &                zeta(Iend  ,j,know))*                             &
     &            0.5_r8*(GRID(ng)%pm(Iend  ,j)+                        &
     &                    GRID(ng)%pm(Iend+1,j))
#  endif
#  ifdef UV_COR
          bry_cor=0.125_r8*(vbar(Iend  ,j  ,know)+                      &
     &                      vbar(Iend  ,j+1,know)+                      &
     &                      vbar(Iend+1,j  ,know)+                      &
     &                      vbar(Iend+1,j+1,know))*                     &
     &                     (GRID(ng)%f(Iend  ,j)+                       &
     &                      GRID(ng)%f(Iend+1,j))
#  else
          bry_cor=0.0_r8
#  endif
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                     &
     &                        zeta(Iend  ,j,know)+                      &
     &                        GRID(ng)%h(Iend+1,j)+                     &
     &                        zeta(Iend+1,j,know)))
          bry_str=cff*(FORCES(ng)%sustr(Iend+1,j)-                      &
     &                 FORCES(ng)%bustr(Iend+1,j))
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,know)+                      &
     &                        dt2d*(bry_pgr+                            &
     &                              bry_cor+                            &
     &                              bry_str)
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ubar(Iend+1,j,kout)=0.0_r8
        END DO
# endif
      END IF
#endif
#if !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

# if defined SOUTH_M2RADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO i=IstrU-1,Iend
          grad(i,Jstr-1)=ubar(i+1,Jstr-1,know)-                         &
     &                   ubar(i  ,Jstr-1,know)
          grad(i,Jstr  )=ubar(i+1,Jstr  ,know)-                         &
     &                   ubar(i  ,Jstr  ,know)
        END DO
        DO i=IstrU,Iend
          dUdt=ubar(i,Jstr,know)-ubar(i,Jstr  ,kout)
          dUde=ubar(i,Jstr,kout)-ubar(i,Jstr+1,kout)
#  ifdef SOUTH_M2NUDGING
          IF ((dUdt*dUde).lt.0.0_r8) THEN
            tau=M2obc_in(ng,isouth)
          ELSE
            tau=M2obc_out(ng,isouth)
          END IF
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,Jstr)+grad(i,Jstr))).gt.0.0_r8) THEN
            dUdx=grad(i-1,Jstr)
          ELSE
            dUdx=grad(i  ,Jstr)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dUdt*dUde
#  if defined CELERITY_WRITE && defined FORWARD_WRITE
          BOUNDARY(ng)%ubar_south_Cx(i)=Cx
          BOUNDARY(ng)%ubar_south_Ce(i)=Ce
          BOUNDARY(ng)%ubar_south_C2(i)=cff
#  endif
          ubar(i,Jstr-1,kout)=(cff*ubar(i,Jstr-1,know)+                 &
     &                         Ce *ubar(i,Jstr  ,kout)-                 &
     &                         MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/        &
     &                        (cff+Ce)
#  ifdef SOUTH_M2NUDGING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)+                      &
     &                        tau*(BOUNDARY(ng)%ubar_south(i)-          &
     &                             ubar(i,Jstr-1,know))
#  endif
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_M2FLATHER || defined SOUTH_M2REDUCED
!
!  Southern edge, Chapman boundary condition.
!
        DO i=IstrU,Iend
          cff=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jstr)+                       &
     &                     GRID(ng)%pn(i  ,Jstr))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jstr)+                     &
     &                        zeta(i-1,Jstr,know)+                      &
     &                        GRID(ng)%h(i  ,Jstr)+                     &
     &                        zeta(i  ,Jstr,know)))
          Ce=cff*cff1
          cff2=1.0_r8/(1.0_r8+Ce)
          ubar(i,Jstr-1,kout)=cff2*(ubar(i,Jstr-1,know)+                &
     &                              Ce*ubar(i,Jstr,kout))
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_M2CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Jstr-1,kout)=BOUNDARY(ng)%ubar_south(i)
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_M2GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Jstr-1,kout)=ubar(i,Jstr,kout)
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO

# else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ubar(i,Jstr-1,kout)=gamma2(ng)*ubar(i,Jstr,kout)
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO
#  undef I_RANGE
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

# if defined NORTH_M2RADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=IstrU-1,Iend
          grad(i,Jend  )=ubar(i+1,Jend  ,know)-                         &
     &                   ubar(i  ,Jend  ,know)
          grad(i,Jend+1)=ubar(i+1,Jend+1,know)-                         &
     &                   ubar(i  ,Jend+1,know)
        END DO
        DO i=IstrU,Iend
          dUdt=ubar(i,Jend,know)-ubar(i,Jend  ,kout)
          dUde=ubar(i,Jend,kout)-ubar(i,Jend-1,kout)
#  ifdef NORTH_M2NUDGING
          IF ((dUdt*dUde).lt.0.0_r8) THEN
            tau=M2obc_in(ng,inorth)
          ELSE
            tau=M2obc_out(ng,inorth)
          END IF
          tau=tau*dt2d
#  endif
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,Jend)+grad(i,Jend))).gt.0.0_r8) THEN
            dUdx=grad(i-1,Jend)
          ELSE
            dUdx=grad(i  ,Jend)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dUdt*dUde
#  if defined CELERITY_WRITE && defined FORWARD_WRITE
          BOUNDARY(ng)%ubar_north_Cx(i)=Cx
          BOUNDARY(ng)%ubar_north_Ce(i)=Ce
          BOUNDARY(ng)%ubar_north_C2(i)=cff
#  endif
          ubar(i,Jend+1,kout)=(cff*ubar(i,Jend+1,know)+                 &
     &                         Ce *ubar(i,Jend  ,kout)-                 &
     &                         MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/        &
     &                        (cff+Ce)
#  ifdef NORTH_M2NUDGING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)+                      &
     &                        tau*(BOUNDARY(ng)%ubar_north(i)-          &
     &                             ubar(i,Jend+1,know))
#  endif
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_M2FLATHER || defined NORTH_M2REDUCED
!
!  Northern edge, Chapman boundary condition.
!
        DO i=IstrU,Iend
          cff=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jend)+                       &
     &                     GRID(ng)%pn(i  ,Jend))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jend)+                     &
     &                        zeta(i-1,Jend,know)+                      &
     &                        GRID(ng)%h(i  ,Jend)+                     &
     &                        zeta(i  ,Jend,know)))
          Ce=cff*cff1
          cff2=1.0_r8/(1.0_r8+Ce)
          ubar(i,Jend+1,kout)=cff2*(ubar(i,Jend+1,know)+                &
     &                              Ce*ubar(i,Jend,kout))
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_M2CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Jend+1,kout)=BOUNDARY(ng)%ubar_north(i)
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_M2GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ubar(i,Jend+1,kout)=ubar(i,Jend,kout)
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO

# else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ubar(i,Jend+1,kout)=gamma2(ng)*ubar(i,Jend,kout)
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO
#  undef I_RANGE
# endif
      END IF
#endif

#if !defined EW_PERIODIC && !defined NS_PERIODIC && !defined COMPOSED_GRID
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        ubar(Istr,Jstr-1,kout)=0.5_r8*(ubar(Istr+1,Jstr-1,kout)+        &
     &                                 ubar(Istr  ,Jstr  ,kout))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar(Iend+1,Jstr-1,kout)=0.5_r8*(ubar(Iend  ,Jstr-1,kout)+      &
     &                                   ubar(Iend+1,Jstr  ,kout))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        ubar(Istr,Jend+1,kout)=0.5_r8*(ubar(Istr  ,Jend  ,kout)+        &
     &                                 ubar(Istr+1,Jend+1,kout))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar(Iend+1,Jend+1,kout)=0.5_r8*(ubar(Iend+1,Jend  ,kout)+      &
     &                                   ubar(Iend  ,Jend+1,kout))
      END IF
#endif

#ifdef REFINED_GRID
!set this for Parent Grid ONLY
      END IF
#endif

#ifdef REFINED_GRID_BC
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
# ifdef SOLVE3D
        t1=nstp(ng)
        t2=nnew(ng)
# else
        IF (MOD(iic(1)*3,2).eq.0) THEN
          t1=2
          t2=1
        ELSE
          t1=1
          t2=2
        ENDIF
# endif
      END IF
!!!!!!!!!!!!!!!!!!!!!!!
! bc's for refined grid
!!!!!!!!!!!!!!!!!!!!!!!
      IF (WESTERN_EDGE) THEN
!
!  Western edge, Flather boundary condition.
!
        DO j=Jstr,Jend
          bry_val=fac1*BOUNDARY(ng)%ubarR_west(4,j,t1)+               &
     &            fac2*BOUNDARY(ng)%ubarR_west(4,j,t2)
# ifdef MASKING
          bry_val=bry_val*GRID(ng)%umask(Istr,j)
# endif
!         bry_val=OCEAN(ng)%ubar(Istr,j,kout)
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                     &
     &                        zeta(Istr-1,j,know)+                      &
     &                        GRID(ng)%h(Istr  ,j)+                     &
     &                        zeta(Istr  ,j,know)))
          Cx=SQRT(g*cff)
          my_val=fac1*BOUNDARY(ng)%zetaR_west(4,j,t1)+                  &
     &           fac2*BOUNDARY(ng)%zetaR_west(4,j,t2)
# ifdef MASKING
          my_val=my_val*GRID(ng)%rmask(Istr-1,j)
# endif
!         my_val=OCEAN(ng)%zeta(Istr-1,j,know)
          ubar(Istr,j,kout)=bry_val-                                    &
     &                      Cx*(0.5_r8*(zeta(Istr-1,j,know)+            &
     &                                  zeta(Istr  ,j,know))-           &
     &                          my_val)
#  ifdef MASKING
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
!
!  Eastern edge, Flather boundary condition.
!
        DO j=Jstr,Jend
          bry_val=fac1*BOUNDARY(ng)%ubarR_east(1,j,t1)+               &
     &            fac2*BOUNDARY(ng)%ubarR_east(1,j,t2)
# ifdef MASKING
          bry_val=bry_val*GRID(ng)%umask(Iend+1,j)
# endif
!         bry_val=OCEAN(ng)%ubar(Iend+1,j,kout)
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                     &
     &                        zeta(Iend  ,j,know)+                      &
     &                        GRID(ng)%h(Iend+1,j)+                     &
     &                        zeta(Iend+1,j,know)))
          Cx=SQRT(g*cff)
          my_val=fac1*BOUNDARY(ng)%zetaR_east(1,j,t1)+                  &
     &           fac2*BOUNDARY(ng)%zetaR_east(1,j,t2)
# ifdef MASKING
          my_val=my_val*GRID(ng)%rmask(Iend+1,j)
# endif
!         my_val=OCEAN(ng)%zeta(Iend+1,j,know)
          ubar(Iend+1,j,kout)=bry_val+                                  &
     &                        Cx*(0.5_r8*(zeta(Iend  ,j,know)+          &
     &                                    zeta(Iend+1,j,know))-         &
     &                            my_val)
#  ifdef MASKING
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO 
      END IF
      IF (SOUTHERN_EDGE) THEN
!
!  Southern edge, Chapman boundary condition.
!
        DO i=IstrU,Iend
          cff=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jstr)+                       &
     &                     GRID(ng)%pn(i  ,Jstr))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jstr)+                     &
     &                        zeta(i-1,Jstr,know)+                      &
     &                        GRID(ng)%h(i  ,Jstr)+                     &
     &                        zeta(i  ,Jstr,know)))
          Ce=cff*cff1
          cff2=1.0_r8/(1.0_r8+Ce)
          ubar(i,Jstr-1,kout)=cff2*(ubar(i,Jstr-1,know)+                &
     &                              Ce*ubar(i,Jstr,kout))
#  ifdef MASKING
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
!
!  Northern edge, Chapman boundary condition.
!
        DO i=IstrU,Iend
          cff=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jend)+                       &
     &                     GRID(ng)%pn(i  ,Jend))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jend)+                     &
     &                        zeta(i-1,Jend,know)+                      &
     &                        GRID(ng)%h(i  ,Jend)+                     &
     &                        zeta(i  ,Jend,know)))
          Ce=cff*cff1
          cff2=1.0_r8/(1.0_r8+Ce)
          ubar(i,Jend+1,kout)=cff2*(ubar(i,Jend+1,know)+                &
     &                              Ce*ubar(i,Jend,kout))
#  ifdef MASKING
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO
      END IF
!
! Corners
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        ubar(Istr,Jstr-1,kout)=0.5_r8*(ubar(Istr+1,Jstr-1,kout)+        &
     &                                 ubar(Istr  ,Jstr  ,kout))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar(Iend+1,Jstr-1,kout)=0.5_r8*(ubar(Iend  ,Jstr-1,kout)+      &
     &                                   ubar(Iend+1,Jstr  ,kout))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        ubar(Istr,Jend+1,kout)=0.5_r8*(ubar(Istr  ,Jend  ,kout)+        &
     &                                 ubar(Istr+1,Jend+1,kout))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        ubar(Iend+1,Jend+1,kout)=0.5_r8*(ubar(Iend+1,Jend  ,kout)+      &
     &                                   ubar(Iend  ,Jend+1,kout))
      END IF

      END IF
#endif   !ref grid bcs
!
#if defined WET_DRY
!
!-----------------------------------------------------------------------
!  Impose wetting and drying conditions.
!-----------------------------------------------------------------------
!
# ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          cff5=GRID(ng)%umask_wet(Istr,j)
          ubar(Istr,j,kout)=ubar(Istr,j,kout)*                          &
     &                   (0.5_r8*cff5*ABS(ABS(cff5)-1.0_r8)+            &
     &                   (0.5_r8+cff5*DSIGN(0.5_r8,ubar(Istr,j,kout)))* &
     &                   (1.0_r8-ABS(ABS(cff5)-1.0_r8)))
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
          cff5=GRID(ng)%umask_wet(Iend+1,j)
          ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                      &
     &                 (0.5_r8*cff5*ABS(ABS(cff5)-1.0_r8)+              &
     &                 (0.5_r8+cff5*DSIGN(0.5_r8,ubar(Iend+1,j,kout)))* &
     &                 (1.0_r8-ABS(ABS(cff5)-1.0_r8)))
        END DO
      END IF
# endif
# ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrU,Iend
          cff5=GRID(ng)%umask_wet(i,Jstr-1)
          ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                      &
     &                 (0.5_r8*cff5*ABS(ABS(cff5)-1.0_r8)+              &
     &                 (0.5_r8+cff5*DSIGN(0.5_r8,ubar(i,Jstr-1,kout)))* &
     &                 (1.0_r8-ABS(ABS(cff5)-1.0_r8)))
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,Iend
          cff5=GRID(ng)%umask_wet(i,Jend+1)
          ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                      &
     &                 (0.5_r8*cff5*ABS(ABS(cff5)-1.0_r8)+              &
     &                 (0.5_r8+cff5*DSIGN(0.5_r8,ubar(i,Jend+1,kout)))* &
     &                 (1.0_r8-ABS(ABS(cff5)-1.0_r8)))
        END DO
      END IF
# endif
# if !defined EW_PERIODIC && !defined NS_PERIODIC
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          cff5=GRID(ng)%umask_wet(Istr,Jstr-1)
          ubar(Istr,Jstr-1,kout)=ubar(Istr,Jstr-1,kout)*                &
     &                   (0.5_r8*cff5*ABS(ABS(cff5)-1.0_r8)+            &
     &              (0.5_r8+cff5*DSIGN(0.5_r8,ubar(Istr,Jstr-1,kout)))* &
     &                   (1.0_r8-ABS(ABS(cff5)-1.0_r8)))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          cff5=GRID(ng)%umask_wet(Iend+1,Jstr-1)
          ubar(Iend+1,Jstr-1,kout)=ubar(Iend+1,Jstr-1,kout)*            &
     &                   (0.5_r8*cff5*ABS(ABS(cff5)-1.0_r8)+            &
     &            (0.5_r8+cff5*DSIGN(0.5_r8,ubar(Iend+1,Jstr-1,kout)))* &
     &                   (1.0_r8-ABS(ABS(cff5)-1.0_r8)))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          cff5=GRID(ng)%umask_wet(Istr,Jend+1)
          ubar(Istr,Jend+1,kout)=ubar(Istr,Jend+1,kout)*                &
     &                   (0.5_r8*cff5*ABS(ABS(cff5)-1.0_r8)+            &
     &            (0.5_r8+cff5*DSIGN(0.5_r8,ubar(Istr,Jend+1,kout)))*   &
     &                   (1.0_r8-ABS(ABS(cff5)-1.0_r8)))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          cff5=GRID(ng)%umask_wet(Iend+1,Jend+1)
          ubar(Iend+1,Jend+1,kout)=ubar(Iend+1,Jend+1,kout)*            &
     &                   (0.5_r8*cff5*ABS(ABS(cff5)-1.0_r8)+            &
     &            (0.5_r8+cff5*DSIGN(0.5_r8,ubar(Iend+1,Jend+1,kout)))* &
     &                   (1.0_r8-ABS(ABS(cff5)-1.0_r8)))
      END IF
# endif
#endif

      RETURN
      END SUBROUTINE u2dbc_tile
      END MODULE u2dbc_mod
