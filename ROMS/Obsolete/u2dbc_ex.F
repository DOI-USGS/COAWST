#include "cppdefs.h"
      MODULE u2dbc_mod
!
!svn $Id: u2dbc_ex.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
!
      PRIVATE
      PUBLIC  :: u2dbc, u2dbc_tile
!
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
     &                 IminS, ImaxS, JminS, JmaxS,                      &
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
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, kout,                          &
     &                       ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: krhs, kstp, kout
!
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
      integer :: Imin, Imax
      integer :: i, j, know

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8) :: cff, cff1, cff2, dt2d, dUde, dUdt, dUdx
      real(r8) :: obc_in, obc_out, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

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
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isUbar,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=ubar(Istr  ,j  ,know)-                       &
     &                     ubar(Istr  ,j-1,know)
            grad(Istr+1,j)=ubar(Istr+1,j  ,know)-                       &
     &                     ubar(Istr+1,j-1,know)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              dUdt=ubar(Istr+1,j,know)-ubar(Istr+1,j,kout)
              dUdx=ubar(Istr+1,j,know)-ubar(Istr+2,j,know)

              IF (LBC(iwest,isUbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(Istr-1,j)+               &
     &                     CLIMA(ng)%M2nudgcof(Istr  ,j))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,iwest)
                  obc_in =M2obc_in (ng,iwest)
                END IF
                IF ((dUdt*dUdx).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
                tau=tau*dt2d
              END IF

              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Istr+1,j  )+                              &
     &                   grad(Istr+1,j+1))).gt.0.0_r8) THEN
                dUde=grad(Istr+1,j  )
              ELSE
                dUde=grad(Istr+1,j+1)
              END IF
              cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=MIN(1.0_r8,cff*dUdx)
#ifdef RADIATION_2D
              Ce=MIN(1.0_r8,MAX(-1.0_r8,cff*dUde))
#else
              Ce=0.0_r8
#endif
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%ubar_west_Cx(j)=Cx
              BOUNDARY(ng)%ubar_west_Ce(j)=Ce
              BOUNDARY(ng)%ubar_west_C2(j)=cff
#endif
              ubar(Istr,j,kout)=(1.0_r8-Cx)*ubar(Istr,j,know)+          &
     &                          Cx*ubar(Istr+1,j,know)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Istr,j  )-          &
     &                          MIN(Ce,0.0_r8)*grad(Istr,j+1)

              IF (LBC(iwest,isUbar,ng)%nudging) THEN
                ubar(Istr,j,kout)=ubar(Istr,j,kout)+                    &
     &                            tau*(BOUNDARY(ng)%ubar_west(j)-       &
     &                                 ubar(Istr,j,know))
              END IF
#ifdef MASKING
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
#endif
            END IF
          END DO
!
!  Western edge, Flather boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%Flather) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
#if defined SSH_TIDES && !defined UV_TIDES
              IF (LBC(iwest,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(Istr,j,know)-                          &
     &                      BOUNDARY(ng)%zeta_west(j))*                 &
     &                  0.5_r8*GRID(ng)%pm(Istr,j)
              ELSE
                bry_pgr=-g*(zeta(Istr  ,j,know)-                        &
     &                      zeta(Istr-1,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Istr-1,j)+                  &
     &                          GRID(ng)%pm(Istr  ,j))
              END IF
# ifdef UV_COR
              bry_cor=0.125_r8*(vbar(Istr-1,j  ,know)+                  &
     &                          vbar(Istr-1,j+1,know)+                  &
     &                          vbar(Istr  ,j  ,know)+                  &
     &                          vbar(Istr  ,j+1,know))*                 &
     &                         (GRID(ng)%f(Istr-1,j)+                   &
     &                          GRID(ng)%f(Istr  ,j))
# else
              bry_cor=0.0_r8
# endif
              cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                &
     &                             zeta(Istr-1,j,know)+                 &
     &                             GRID(ng)%h(Istr  ,j)+                &
     &                             zeta(Istr  ,j,know)))
              bry_str=cff1*(FORCES(ng)%sustr(Istr,j)-                   &
     &                      FORCES(ng)%bustr(Istr,j))
              Cx=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(Istr-1,j)+            &
     &                                 zeta(Istr-1,j,know)+             &
     &                                 GRID(ng)%h(Istr  ,j)+            &
     &                                 zeta(Istr  ,j,know)))
              cff2=GRID(ng)%om_u(Istr,j)*Cx
!!            cff2=dt2d
              bry_val=ubar(Istr+1,j,know)+                              &
     &                cff2*(bry_pgr+                                    &
     &                      bry_cor+                                    &
     &                      bry_str)
#else
              bry_val=BOUNDARY(ng)%ubar_west(j)
#endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                 &
     &                            zeta(Istr-1,j,know)+                  &
     &                            GRID(ng)%h(Istr  ,j)+                 &
     &                            zeta(Istr  ,j,know)))
              Cx=SQRT(g*cff)
              ubar(Istr,j,kout)=bry_val-                                &
     &                          Cx*(0.5_r8*(zeta(Istr-1,j,know)+        &
     &                                      zeta(Istr  ,j,know))-       &
     &                              BOUNDARY(ng)%zeta_west(j))
#ifdef MASKING
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
#endif
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar(Istr,j,kout)=BOUNDARY(ng)%ubar_west(j)
#ifdef MASKING
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
#endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar(Istr,j,kout)=ubar(Istr+1,j,kout)
#ifdef MASKING
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
#endif
            END IF
          END DO
!
!  Western edge, reduced-physics boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%reduced) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              IF (LBC(iwest,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(Istr,j,know)-                          &
     &                      BOUNDARY(ng)%zeta_west(j))*                 &
     &                  0.5_r8*GRID(ng)%pm(Istr,j)
              ELSE
                bry_pgr=-g*(zeta(Istr  ,j,know)-                        &
     &                      zeta(Istr-1,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Istr-1,j)+                  &
     &                          GRID(ng)%pm(Istr  ,j))
              END IF
#ifdef UV_COR
              bry_cor=0.125_r8*(vbar(Istr-1,j  ,know)+                  &
     &                          vbar(Istr-1,j+1,know)+                  &
     &                          vbar(Istr  ,j  ,know)+                  &
     &                          vbar(Istr  ,j+1,know))*                 &
     &                         (GRID(ng)%f(Istr-1,j)+                   &
     &                          GRID(ng)%f(Istr  ,j))
#else
              bry_cor=0.0_r8
#endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                 &
     &                            zeta(Istr-1,j,know)+                  &
     &                            GRID(ng)%h(Istr  ,j)+                 &
     &                            zeta(Istr  ,j,know)))
              bry_str=cff*(FORCES(ng)%sustr(Istr,j)-                    &
     &                     FORCES(ng)%bustr(Istr,j))
              ubar(Istr,j,kout)=ubar(Istr,j,know)+                      &
     &                          dt2d*(bry_pgr+                          &
     &                                bry_cor+                          &
     &                                bry_str)
#ifdef MASKING
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
#endif
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar(Istr,j,kout)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (LBC(ieast,isUbar,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=ubar(Iend  ,j  ,know)-                       &
     &                     ubar(Iend  ,j-1,know)
            grad(Iend+1,j)=ubar(Iend+1,j  ,know)-                       &
     &                     ubar(Iend+1,j-1,know)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              dUdt=ubar(Iend,j,know)-ubar(Iend  ,j,kout)
              dUdx=ubar(Iend,j,know)-ubar(Iend-1,j,know)

              IF (LBC(ieast,isUbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(Iend  ,j)+               &
     &                     CLIMA(ng)%M2nudgcof(Iend+1,j))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,ieast)
                  obc_in =M2obc_in (ng,ieast)
                END IF
                IF ((dUdt*dUdx).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
                tau=tau*dt2d
              END IF

              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Iend,j  )+                                &
     &                   grad(Iend,j+1))).gt.0.0_r8) THEN
                dUde=grad(Iend,j)
              ELSE
                dUde=grad(Iend,j+1)
              END IF
              cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=MIN(1.0_r8,cff*dUdx)
#ifdef RADIATION_2D
              Ce=MIN(1.0_r8,MAX(-1.0_r8,cff*dUde))
#else
              Ce=0.0_r8
#endif
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%ubar_east_Cx(j)=Cx
              BOUNDARY(ng)%ubar_east_Ce(j)=Ce
              BOUNDARY(ng)%ubar_east_C2(j)=cff
#endif
              ubar(Iend+1,j,kout)=(1.0_r8-Cx)*ubar(Iend+1,j,know)+      &
     &                            Cx*ubar(Iend,j,know)-                 &
     &                            MAX(Ce,0.0_r8)*grad(Iend+1,j  )-      &
     &                            MIN(Ce,0.0_r8)*grad(Iend+1,j+1)

              IF (LBC(ieast,isUbar,ng)%nudging) THEN
                ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)+                &
     &                              tau*(BOUNDARY(ng)%ubar_east(j)-     &
     &                                   ubar(Iend+1,j,know))
              END IF
#ifdef MASKING
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, Flather boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%Flather) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
#if defined SSH_TIDES && !defined UV_TIDES
              IF (LBC(ieast,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_east(j)-                  &
     &                      zeta(Iend,j,know))*                         &
     &                  0.5_r8*GRID(ng)%pm(Iend,j)
              ELSE
                bry_pgr=-g*(zeta(Iend+1,j,know)-                        &
     &                      zeta(Iend  ,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Iend  ,j)+                  &
     &                          GRID(ng)%pm(Iend+1,j))
              END IF
# ifdef UV_COR
              bry_cor=0.125_r8*(vbar(Iend  ,j  ,know)+                  &
     &                          vbar(Iend  ,j+1,know)+                  &
     &                          vbar(Iend+1,j  ,know)+                  &
     &                          vbar(Iend+1,j+1,know))*                 &
     &                         (GRID(ng)%f(Iend  ,j)+                   &
     &                          GRID(ng)%f(Iend+1,j))
# else
              bry_cor=0.0_r8
# endif
              cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                &
     &                             zeta(Iend  ,j,know)+                 &
     &                             GRID(ng)%h(Iend+1,j)+                &
     &                             zeta(Iend+1,j,know)))
              bry_str=cff1*(FORCES(ng)%sustr(Iend+1,j)-                 &
     &                      FORCES(ng)%bustr(Iend+1,j))
              Cx=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(Iend+1,j)+            &
     &                                 zeta(Iend+1,j,know)+             &
     &                                 GRID(ng)%h(Iend  ,j)+            &
     &                                 zeta(Iend  ,j,know)))
              cff2=GRID(ng)%om_u(Iend+1,j)*Cx
!!            cff2=dt2d
              bry_val=ubar(Iend,j,know)+                                &
     &                cff2*(bry_pgr+                                    &
     &                      bry_cor+                                    &
     &                      bry_str)
#else
              bry_val=BOUNDARY(ng)%ubar_east(j)
#endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                 &
     &                            zeta(Iend  ,j,know)+                  &
     &                            GRID(ng)%h(Iend+1,j)+                 &
     &                            zeta(Iend+1,j,know)))
              Cx=SQRT(g*cff)
              ubar(Iend+1,j,kout)=bry_val+                              &
     &                            Cx*(0.5_r8*(zeta(Iend  ,j,know)+      &
     &                                        zeta(Iend+1,j,know))-     &
     &                                BOUNDARY(ng)%zeta_east(j))
#ifdef MASKING
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar(Iend+1,j,kout)=BOUNDARY(ng)%ubar_east(j)
#ifdef MASKING
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar(Iend+1,j,kout)=ubar(Iend,j,kout)
#ifdef MASKING
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, reduced-physics boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%reduced) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              IF (LBC(ieast,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_east(j)-                  &
     &                      zeta(Iend,j,know))*                         &
     &                  0.5_r8*GRID(ng)%pm(Iend,j)
              ELSE
                bry_pgr=-g*(zeta(Iend+1,j,know)-                        &
     &                      zeta(Iend  ,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Iend  ,j)+                  &
     &                          GRID(ng)%pm(Iend+1,j))
              END IF
#ifdef UV_COR
              bry_cor=0.125_r8*(vbar(Iend  ,j  ,know)+                  &
     &                          vbar(Iend  ,j+1,know)+                  &
     &                          vbar(Iend+1,j  ,know)+                  &
     &                          vbar(Iend+1,j+1,know))*                 &
     &                         (GRID(ng)%f(Iend  ,j)+                   &
     &                          GRID(ng)%f(Iend+1,j))
#else
              bry_cor=0.0_r8
#endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                 &
     &                            zeta(Iend  ,j,know)+                  &
     &                            GRID(ng)%h(Iend+1,j)+                 &
     &                            zeta(Iend+1,j,know)))
              bry_str=cff*(FORCES(ng)%sustr(Iend+1,j)-                  &
     &                     FORCES(ng)%bustr(Iend+1,j))
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,know)+                  &
     &                            dt2d*(bry_pgr+                        &
     &                                  bry_cor+                        &
     &                                  bry_str)
#ifdef MASKING
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar(Iend+1,j,kout)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isUbar,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=ubar(i+1,Jstr-1,know)-                       &
     &                     ubar(i  ,Jstr-1,know)
            grad(i,Jstr  )=ubar(i+1,Jstr  ,know)-                       &
     &                     ubar(i  ,Jstr  ,know)
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              dUdt=ubar(i,Jstr,know)-ubar(i,Jstr  ,kout)
              dUde=ubar(i,Jstr,know)-ubar(i,Jstr+1,know)

              IF (LBC(isouth,isUbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(i-1,Jstr-1)+             &
     &                     CLIMA(ng)%M2nudgcof(i  ,Jstr-1))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,isouth)
                  obc_in =M2obc_in (ng,isouth)
                END IF
                IF ((dUdt*dUde).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
                tau=tau*dt2d
              END IF

              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jstr)+                                &
     &                   grad(i  ,Jstr))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jstr)
              ELSE
                dUdx=grad(i  ,Jstr)
              END IF
              cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
#ifdef RADIATION_2D
              Cx=MIN(1.0_r8,MAX(-1.0_r8,cff*dUdx))
#else
              Cx=0.0_r8
#endif
              Ce=MIN(1.0_r8,cff*dUde)
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%ubar_south_Cx(i)=Cx
              BOUNDARY(ng)%ubar_south_Ce(i)=Ce
              BOUNDARY(ng)%ubar_south_C2(i)=cff
#endif
              ubar(i,Jstr-1,kout)=(1.0_r8-Ce)*ubar(i,Jstr-1,know)+      &
     &                            Ce*ubar(i,Jstr,know)-                 &
     &                            MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-      &
     &                            MIN(Cx,0.0_r8)*grad(i  ,Jstr-1)

              IF (LBC(isouth,isUbar,ng)%nudging) THEN
                ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)+                &
     &                              tau*(BOUNDARY(ng)%ubar_south(i)-    &
     &                                   ubar(i,Jstr-1,know))
              END IF
#ifdef MASKING
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
#endif
            END IF
          END DO
!
!  Southern edge, Chapman boundary condition.
!
        ELSE IF (LBC(isouth,isUbar,ng)%Flather.or.                      &
     &           LBC(isouth,isUbar,ng)%reduced) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jstr)+                   &
     &                         GRID(ng)%pn(i  ,Jstr))
              cff1=SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jstr)+                 &
     &                            zeta(i-1,Jstr,know)+                  &
     &                            GRID(ng)%h(i  ,Jstr)+                 &
     &                            zeta(i  ,Jstr,know)))
              Ce=cff*cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              ubar(i,Jstr-1,kout)=cff2*(ubar(i,Jstr-1,know)+            &
     &                                  Ce*ubar(i,Jstr,kout))
#ifdef MASKING
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
#endif
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isUbar,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ubar(i,Jstr-1,kout)=BOUNDARY(ng)%ubar_south(i)
#ifdef MASKING
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
#endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isUbar,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ubar(i,Jstr-1,kout)=ubar(i,Jstr,kout)
#ifdef MASKING
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
#endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(isouth,isUbar,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%south(i)) THEN
              ubar(i,Jstr-1,kout)=gamma2(ng)*ubar(i,Jstr,kout)
#ifdef MASKING
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
#endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (LBC(inorth,isUbar,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jend  )=ubar(i+1,Jend  ,know)-                       &
     &                     ubar(i  ,Jend  ,know)
            grad(i,Jend+1)=ubar(i+1,Jend+1,know)-                       &
     &                     ubar(i  ,Jend+1,know)
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              dUdt=ubar(i,Jend,know)-ubar(i,Jend  ,kout)
              dUde=ubar(i,Jend,know)-ubar(i,Jend-1,know)

              IF (LBC(inorth,isUbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(i-1,Jend+1)+             &
     &                     CLIMA(ng)%M2nudgcof(i  ,Jend+1))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,inorth)
                  obc_in =M2obc_in (ng,inorth)
                END IF
                IF ((dUdt*dUde).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
                tau=tau*dt2d
              END IF

              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jend)+                                &
     &                   grad(i  ,Jend))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jend)
              ELSE
                dUdx=grad(i  ,Jend)
              END IF
              cff=dUdt/MAX(dUdx*dUdx+dUde*dUde,eps)
#ifdef RADIATION_2D
              Cx=MIN(1.0_r8,MAX(-1.0_r8,cff*dUdx))
#else
              Cx=0.0_r8
#endif
              Ce=MIN(1.0_r8,cff*dUde)
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%ubar_north_Cx(i)=Cx
              BOUNDARY(ng)%ubar_north_Ce(i)=Ce
              BOUNDARY(ng)%ubar_north_C2(i)=cff
#endif
              ubar(i,Jend+1,kout)=(1.0_r8-Ce)*ubar(i,Jend+1,know)+      &
     &                            Ce*ubar(i,Jend,know)-                 &
     &                            MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-      &
     &                            MIN(Cx,0.0_r8)*grad(i  ,Jend+1)

              IF (LBC(inorth,isUbar,ng)%nudging) THEN
                ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)+                &
     &                              tau*(BOUNDARY(ng)%ubar_north(i)-    &
     &                                   ubar(i,Jend+1,know))
              END IF
#ifdef MASKING
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, Chapman boundary condition.
!
        ELSE IF (LBC(inorth,isUbar,ng)%Flather.or.                      &
     &           LBC(inorth,isUbar,ng)%reduced) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jend)+                   &
     &                         GRID(ng)%pn(i  ,Jend))
              cff1=SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jend)+                 &
     &                            zeta(i-1,Jend,know)+                  &
     &                            GRID(ng)%h(i  ,Jend)+                 &
     &                            zeta(i  ,Jend,know)))
              Ce=cff*cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              ubar(i,Jend+1,kout)=cff2*(ubar(i,Jend+1,know)+            &
     &                                  Ce*ubar(i,Jend,kout))
#ifdef MASKING
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isUbar,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ubar(i,Jend+1,kout)=BOUNDARY(ng)%ubar_north(i)
#ifdef MASKING
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isUbar,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ubar(i,Jend+1,kout)=ubar(i,Jend,kout)
#ifdef MASKING
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(inorth,isUbar,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%north(i)) THEN
              ubar(i,Jend+1,kout)=gamma2(ng)*ubar(i,Jend,kout)
#ifdef MASKING
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
#endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            ubar(Istr,Jstr-1,kout)=0.5_r8*(ubar(Istr+1,Jstr-1,kout)+    &
     &                                     ubar(Istr  ,Jstr  ,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            ubar(Iend+1,Jstr-1,kout)=0.5_r8*(ubar(Iend  ,Jstr-1,kout)+  &
     &                                       ubar(Iend+1,Jstr  ,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            ubar(Istr,Jend+1,kout)=0.5_r8*(ubar(Istr  ,Jend  ,kout)+    &
     &                                     ubar(Istr+1,Jend+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            ubar(Iend+1,Jend+1,kout)=0.5_r8*(ubar(Iend+1,Jend  ,kout)+  &
     &                                       ubar(Iend  ,Jend+1,kout))
          END IF
        END IF
      END IF

#if defined WET_DRY
!
!-----------------------------------------------------------------------
!  Impose wetting and drying conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff1=ABS(ABS(GRID(ng)%umask_wet(Istr,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,ubar(Istr,j,kout))*              &
     &                    GRID(ng)%umask_wet(Istr,j)
              cff=0.5_r8*GRID(ng)%umask_wet(Istr,j)*cff1+               &
     &            cff2*(1.0_r8-cff1)
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff1=ABS(ABS(GRID(ng)%umask_wet(Iend+1,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,ubar(Iend+1,j,kout))*            &
     &                    GRID(ng)%umask_wet(Iend+1,j)
              cff=0.5_r8*GRID(ng)%umask_wet(Iend+1,j)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*cff
            END IF
          END DO
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff1=ABS(ABS(GRID(ng)%umask_wet(i,Jstr-1))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,ubar(i,Jstr-1,kout))*            &
     &                    GRID(ng)%umask_wet(i,Jstr-1)
              cff=0.5_r8*GRID(ng)%umask_wet(i,Jstr-1)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff1=ABS(ABS(GRID(ng)%umask_wet(i,Jend+1))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,ubar(i,Jend+1,kout))*            &
     &                    GRID(ng)%umask_wet(i,Jend+1)
              cff=0.5_r8*GRID(ng)%umask_wet(i,Jend+1)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*cff
            END IF
          END DO
        END IF
      END IF

      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            cff1=ABS(ABS(GRID(ng)%umask_wet(Istr,Jstr-1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Istr,Jstr-1,kout))*           &
     &                  GRID(ng)%umask_wet(Istr,Jstr-1)
            cff=0.5_r8*GRID(ng)%umask_wet(Istr,Jstr-1)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            ubar(Istr,Jstr-1,kout)=ubar(Istr,Jstr-1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            cff1=ABS(ABS(GRID(ng)%umask_wet(Iend+1,Jstr-1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Iend+1,Jstr-1,kout))*         &
     &                  GRID(ng)%umask_wet(Iend+1,Jstr-1)
            cff=0.5_r8*GRID(ng)%umask_wet(Iend+1,Jstr-1)*cff1+          &
     &          cff2*(1.0_r8-cff1)
            ubar(Iend+1,Jstr-1,kout)=ubar(Iend+1,Jstr-1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            cff1=ABS(ABS(GRID(ng)%umask_wet(Istr,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Istr,Jend+1,kout))*           &
     &                  GRID(ng)%umask_wet(Istr,Jend+1)
            cff=0.5_r8*GRID(ng)%umask_wet(Istr,Jend+1)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            ubar(Istr,Jend+1,kout)=ubar(Istr,Jend+1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            cff1=ABS(ABS(GRID(ng)%umask_wet(Iend+1,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Iend+1,Jend+1,kout))*         &
     &                  GRID(ng)%umask_wet(Iend+1,Jend+1)
            cff=0.5_r8*GRID(ng)%umask_wet(Iend+1,Jend+1)+cff1+          &
     &          cff2*(1.0_r8-cff1)
            ubar(Iend+1,Jend+1,kout)=ubar(Iend+1,Jend+1,kout)*cff
          END IF
        END IF
      END IF
#endif

      RETURN
      END SUBROUTINE u2dbc_tile
      END MODULE u2dbc_mod
