#include "cppdefs.h"
      MODULE v2dbc_mod
!
!svn $Id: v2dbc_ex.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
!
      PRIVATE
      PUBLIC  :: v2dbc, v2dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE v2dbc (ng, tile, kout)
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
      CALL v2dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 krhs(ng), kstp(ng), kout,                        &
     &                 OCEAN(ng) % ubar,                                &
     &                 OCEAN(ng) % vbar,                                &
     &                 OCEAN(ng) % zeta)

      RETURN
      END SUBROUTINE v2dbc
!
!***********************************************************************
      SUBROUTINE v2dbc_tile (ng, tile,                                  &
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
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)

      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j, know

      real(r8), parameter :: eps = 1.0E-20_r8

      real(r8) :: Ce, Cx
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8):: cff, cff1, cff2, dt2d, dVde, dVdt, dVdx
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
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isVbar,ng)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jstr  )=vbar(i  ,Jstr  ,know)-                       &
     &                     vbar(i-1,Jstr  ,know)
            grad(i,Jstr+1)=vbar(i  ,Jstr+1,know)-                       &
     &                     vbar(i-1,Jstr+1,know)
          END DO
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              dVdt=vbar(i,Jstr+1,know)-vbar(i,Jstr+1,kout)
              dVde=vbar(i,Jstr+1,know)-vbar(i,Jstr+2,know)

              IF (LBC(isouth,isVbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(i,Jstr-1)+               &
     &                     CLIMA(ng)%M2nudgcof(i,Jstr  ))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,isouth)
                  obc_in =M2obc_in (ng,isouth)
                END IF
                IF ((dVdt*dVde).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
                tau=tau*dt2d
              END IF

              IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(i  ,Jstr+1)+                              &
     &                   grad(i+1,Jstr+1))).gt.0.0_r8) THEN
                dVdx=grad(i  ,Jstr+1)
              ELSE
                dVdx=grad(i+1,Jstr+1)
              END IF
              cff=dVdt/MAX(dVdx*dVdx+dVde*dVde,eps)
#ifdef RADIATION_2D
              Cx=MIN(1.0_r8,MAX(-1.0_r8,cff*dVdx))
#else
              Cx=0.0_r8
#endif
              Ce=MIN(1.0_r8,cff*dVde)
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%vbar_south_Cx(i)=Cx
              BOUNDARY(ng)%vbar_south_Ce(i)=Ce
              BOUNDARY(ng)%vbar_south_C2(i)=cff
#endif
              vbar(i,Jstr,kout)=(1.0_r8-Ce)*vbar(i,Jstr,know)+          &
     &                          Ce*vbar(i,Jstr+1,know)-                 &
     &                          MAX(Cx,0.0_r8)*grad(i  ,Jstr)-          &
     &                          MIN(Cx,0.0_r8)*grad(i+1,Jstr)

              IF (LBC(isouth,isVbar,ng)%nudging) THEN
                vbar(i,Jstr,kout)=vbar(i,Jstr,kout)+                    &
     &                            tau*(BOUNDARY(ng)%vbar_south(i)-      &
     &                                 vbar(i,Jstr,know))
              END IF
#ifdef MASKING
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
#endif
            END IF
          END DO
!
!  Southern edge, Flather boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
#if defined SSH_TIDES && !defined UV_TIDES
              IF (LBC(isouth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(i,Jstr,know)-                          &
     &                      BOUNDARY(ng)%zeta_south(i))*                &
     &                  0.5_r8*GRID(ng)%pn(i,Jstr)
              ELSE
                bry_pgr=-g*(zeta(i,Jstr  ,know)-                        &
     &                      zeta(i,Jstr-1,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                  &
     &                          GRID(ng)%pn(i,Jstr  ))
              END IF
# ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jstr-1,know)+                 &
     &                           ubar(i+1,Jstr-1,know)+                 &
     &                           ubar(i  ,Jstr  ,know)+                 &
     &                           ubar(i+1,Jstr  ,know))*                &
     &                          (GRID(ng)%f(i,Jstr-1)+                  &
     &                           GRID(ng)%f(i,Jstr  ))
# else
              bry_cor=0.0_r8
# endif
              cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                &
     &                             zeta(i,Jstr-1,know)+                 &
     &                             GRID(ng)%h(i,Jstr  )+                &
     &                             zeta(i,Jstr  ,know)))
              bry_str=cff1*(FORCES(ng)%svstr(i,Jstr)-                   &
     &                      FORCES(ng)%bvstr(i,Jstr))
              Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(i,Jstr-1)+            &
     &                                 zeta(i,Jstr-1,know)+             &
     &                                 GRID(ng)%h(i,Jstr  )+            &
     &                                 zeta(i,Jstr  ,know)))
              cff2=GRID(ng)%on_v(i,Jstr)*Ce
!!            cff2=dt2d
              bry_val=vbar(i,Jstr+1,know)+                              &
     &                cff2*(bry_pgr+                                    &
     &                      bry_cor+                                    &
     &                      bry_str)
#else
              bry_val=BOUNDARY(ng)%vbar_south(i)
#endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                 &
     &                            zeta(i,Jstr-1,know)+                  &
     &                            GRID(ng)%h(i,Jstr  )+                 &
     &                            zeta(i,Jstr  ,know)))
              Ce=SQRT(g*cff)
              vbar(i,Jstr,kout)=bry_val-                                &
     &                          Ce*(0.5_r8*(zeta(i,Jstr-1,know)+        &
     &                                      zeta(i,Jstr  ,know))-       &
     &                              BOUNDARY(ng)%zeta_south(i))
#ifdef MASKING
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
#endif
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              vbar(i,Jstr,kout)=BOUNDARY(ng)%vbar_south(i)
#ifdef MASKING
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
#endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              vbar(i,Jstr,kout)=vbar(i,Jstr+1,kout)
#ifdef MASKING
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
#endif
            END IF
          END DO
!
!  Southern edge, reduced-physics boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              IF (LBC(isouth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(i,Jstr,know)-                          &
     &                      BOUNDARY(ng)%zeta_south(i))*                &
     &                  0.5_r8*GRID(ng)%pn(i,Jstr)
              ELSE
                bry_pgr=-g*(zeta(i,Jstr  ,know)-                        &
     &                      zeta(i,Jstr-1,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                  &
     &                          GRID(ng)%pn(i,Jstr  ))
              END IF
#ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jstr-1,know)+                 &
     &                           ubar(i+1,Jstr-1,know)+                 &
     &                           ubar(i  ,Jstr  ,know)+                 &
     &                           ubar(i+1,Jstr  ,know))*                &
     &                          (GRID(ng)%f(i,Jstr-1)+                  &
     &                           GRID(ng)%f(i,Jstr  ))
#else
              bry_cor=0.0_r8
#endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                 &
     &                            zeta(i,Jstr-1,know)+                  &
     &                            GRID(ng)%h(i,Jstr  )+                 &
     &                            zeta(i,Jstr  ,know)))
              bry_str=cff*(FORCES(ng)%svstr(i,Jstr)-                    &
     &                     FORCES(ng)%bvstr(i,Jstr))
              vbar(i,Jstr,kout)=vbar(i,Jstr,know)+                      &
     &                          dt2d*(bry_pgr+                          &
     &                                bry_cor+                          &
     &                                bry_str)
#ifdef MASKING
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
#endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              vbar(i,Jstr,kout)=0.0_r8
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
        IF (LBC(inorth,isVbar,ng)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jend  )=vbar(i  ,Jend  ,know)-                       &
     &                     vbar(i-1,Jend  ,know)
            grad(i,Jend+1)=vbar(i  ,Jend+1,know)-                       &
     &                     vbar(i-1,Jend+1,know)
          END DO
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              dVdt=vbar(i,Jend,know)-vbar(i,Jend  ,kout)
              dVde=vbar(i,Jend,know)-vbar(i,Jend-1,know)

              IF (LBC(inorth,isVbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(i,Jend  )+               &
     &                     CLIMA(ng)%M2nudgcof(i,Jend+1))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,inorth)
                  obc_in =M2obc_in (ng,inorth)
                END IF
                IF ((dVdt*dVde).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
                tau=tau*dt2d
              END IF

              IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(i  ,Jend)+                                &
     &                   grad(i+1,Jend))).gt.0.0_r8) THEN
                dVdx=grad(i  ,Jend)
              ELSE
                dVdx=grad(i+1,Jend)
              END IF
              cff=dVdt/MAX(dVdx*dVdx+dVde*dVde,eps)
#ifdef RADIATION_2D
              Cx=MIN(1.0_r8,MAX(-1.0_r8,cff*dVdx))
#else
              Cx=0.0_r8
#endif
              Ce=MIN(1.0_r8,cff*dVde)
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%vbar_north_Cx(i)=Cx
              BOUNDARY(ng)%vbar_north_Ce(i)=Ce
              BOUNDARY(ng)%vbar_north_C2(i)=cff
#endif
              vbar(i,Jend+1,kout)=(1.0_r8-Ce)*vbar(i,Jend+1,know)+      &
     &                            Ce*vbar(i,Jend,know)-                 &
     &                            MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-      &
     &                            MIN(Cx,0.0_r8)*grad(i+1,Jend+1)

              IF (LBC(inorth,isVbar,ng)%nudging) THEN
                vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)+                &
     &                              tau*(BOUNDARY(ng)%vbar_north(i)-    &
     &                                   vbar(i,Jend+1,know))
              END IF
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, Flather boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
#if defined SSH_TIDES && !defined UV_TIDES
              IF (LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                 &
     &                      zeta(i,Jend,know))*                         &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
              ELSE
                bry_pgr=-g*(zeta(i,Jend+1,know)-                        &
     &                      zeta(i,Jend  ,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
              END IF
# ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jend  ,know)+                 &
     &                           ubar(i+1,Jend  ,know)+                 &
     &                           ubar(i  ,Jend+1,know)+                 &
     &                           ubar(i+1,Jend+1,know))*                &
     &                          (GRID(ng)%f(i,Jend  )+                  &
     &                           GRID(ng)%f(i,Jend+1))
# else
              bry_cor=0.0_r8
# endif
              cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                &
     &                             zeta(i,Jend  ,know)+                 &
     &                             GRID(ng)%h(i,Jend+1)+                &
     &                             zeta(i,Jend+1,know)))
              bry_str=cff1*(FORCES(ng)%svstr(i,Jend+1)-                 &
     &                      FORCES(ng)%bvstr(i,Jend+1))
              Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(i,Jend+1)+            &
     &                                 zeta(i,Jend+1,know)+             &
     &                                 GRID(ng)%h(i,Jend  )+            &
     &                                 zeta(i,Jend  ,know)))
              cff2=GRID(ng)%on_v(i,Jend+1)*Ce
!!            cff2=dt2d
              bry_val=vbar(i,Jend,know)+                                &
     &                cff2*(bry_pgr+                                    &
     &                      bry_cor+                                    &
     &                      bry_str)
#else
              bry_val=BOUNDARY(ng)%vbar_north(i)
#endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                 &
     &                            zeta(i,Jend  ,know)+                  &
     &                            GRID(ng)%h(i,Jend+1)+                 &
     &                            zeta(i,Jend+1,know)))
              Ce=SQRT(g*cff)
              vbar(i,Jend+1,kout)=bry_val+                              &
     &                            Ce*(0.5_r8*(zeta(i,Jend  ,know)+      &
     &                                        zeta(i,Jend+1,know))-     &
     &                                BOUNDARY(ng)%zeta_north(i))
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              vbar(i,Jend+1,kout)=BOUNDARY(ng)%vbar_north(i)
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              vbar(i,Jend+1,kout)=vbar(i,Jend,kout)
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, reduced-physics boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              IF (LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                 &
     &                      zeta(i,Jend,know))*                         &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
              ELSE
                bry_pgr=-g*(zeta(i,Jend+1,know)-                        &
     &                      zeta(i,Jend  ,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
              END IF
#ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jend  ,know)+                 &
     &                           ubar(i+1,Jend  ,know)+                 &
     &                           ubar(i  ,Jend+1,know)+                 &
     &                           ubar(i+1,Jend+1,know))*                &
     &                          (GRID(ng)%f(i,Jend  )+                  &
     &                           GRID(ng)%f(i,Jend+1))
#else
              bry_cor=0.0_r8
#endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                 &
     &                            zeta(i,Jend  ,know)+                  &
     &                            GRID(ng)%h(i,Jend+1)+                 &
     &                            zeta(i,Jend+1,know)))
              bry_str=cff*(FORCES(ng)%svstr(i,Jend+1)-                  &
     &                     FORCES(ng)%bvstr(i,Jend+1))
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,know)+                  &
     &                            dt2d*(bry_pgr+                        &
     &                                  bry_cor+                        &
     &                                  bry_str)
#ifdef MASKING
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
#endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              vbar(i,Jend+1,kout)=0.0_r8
            END IF
          END DO
        END IF
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
        IF (LBC(iwest,isVbar,ng)%radiation) THEN
          DO j=JstrV-1,Jend
            grad(Istr-1,j)=vbar(Istr-1,j+1,know)-                       &
     &                     vbar(Istr-1,j  ,know)
            grad(Istr  ,j)=vbar(Istr  ,j+1,know)-                       &
     &                     vbar(Istr  ,j  ,know)
          END DO
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              dVdt=vbar(Istr,j,know)-vbar(Istr  ,j,kout)
              dVdx=vbar(Istr,j,know)-vbar(Istr+1,j,know)

              IF (LBC(iwest,isVbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(Istr-1,j-1)+             &
     &                     CLIMA(ng)%M2nudgcof(Istr-1,j  ))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,iwest)
                  obc_in =M2obc_in (ng,iwest)
                END IF
                IF ((dVdt*dVdx).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
                tau=tau*dt2d
              END IF

              IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(Istr,j-1)+                                &
     &                   grad(Istr,j  ))).gt.0.0_r8) THEN
                dVde=grad(Istr,j-1)
              ELSE
                dVde=grad(Istr,j  )
              END IF
              cff=dVdt/MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=MIN(1.0_r8,cff*dVdx)
#ifdef RADIATION_2D
              Ce=MIN(1.0_r8,MAX(-1.0_r8,cff*dVde))
#else
              Ce=0.0_r8
#endif
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%vbar_west_Cx(j)=Cx
              BOUNDARY(ng)%vbar_west_Ce(j)=Ce
              BOUNDARY(ng)%vbar_west_C2(j)=cff
#endif
              vbar(Istr-1,j,kout)=(1.0_r8-Cx)*vbar(Istr-1,j,know)+      &
     &                            Cx*vbar(Istr,j,know)-                 &
     &                            MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-      &
     &                            MIN(Ce,0.0_r8)*grad(Istr-1,j  )

              IF (LBC(iwest,isVbar,ng)%nudging) THEN
                vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)+                &
     &                              tau*(BOUNDARY(ng)%vbar_west(j)-     &
     &                                   vbar(Istr-1,j,know))
              END IF
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, Chapman boundary condition.
!
        ELSE IF (LBC(iwest,isVbar,ng)%Flather.or.                       &
     &           LBC(iwest,isVbar,ng)%reduced) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pm(Istr,j-1)+                   &
     &                         GRID(ng)%pm(Istr,j  ))
              cff1=SQRT(g*0.5_r8*(GRID(ng)%h(Istr,j-1)+                 &
     &                            zeta(Istr,j-1,know)+                  &
     &                            GRID(ng)%h(Istr,j  )+                 &
     &                            zeta(Istr,j  ,know)))
              Cx=cff*cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              vbar(Istr-1,j,kout)=cff2*(vbar(Istr-1,j,know)+            &
     &                                  Cx*vbar(Istr,j,kout))
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              vbar(Istr-1,j,kout)=BOUNDARY(ng)%vbar_west(j)
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              vbar(Istr-1,j,kout)=vbar(Istr,j,kout)
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
            END IF
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(iwest,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%west(j)) THEN
              vbar(Istr-1,j,kout)=gamma2(ng)*vbar(Istr,j,kout)
#ifdef MASKING
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
#endif
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
        IF (LBC(ieast,isVbar,ng)%radiation) THEN
          DO j=JstrV-1,Jend
            grad(Iend  ,j)=vbar(Iend  ,j+1,know)-                       &
     &                     vbar(Iend  ,j  ,know)
            grad(Iend+1,j)=vbar(Iend+1,j+1,know)-                       &
     &                     vbar(Iend+1,j  ,know)
          END DO
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              dVdt=vbar(Iend,j,know)-vbar(Iend  ,j,kout)
              dVdx=vbar(Iend,j,know)-vbar(Iend-1,j,know)

              IF (LBC(ieast,isVbar,ng)%nudging) THEN
                IF (LnudgeM2CLM(ng)) THEN
                  obc_out=0.5_r8*                                       &
     &                    (CLIMA(ng)%M2nudgcof(Iend+1,j-1)+             &
     &                     CLIMA(ng)%M2nudgcof(Iend+1,j  ))
                  obc_in =obcfac(ng)*obc_out
                ELSE
                  obc_out=M2obc_out(ng,ieast)
                  obc_in =M2obc_in (ng,ieast)
                END IF
                IF ((dVdt*dVdx).lt.0.0_r8) THEN
                  tau=obc_in
                ELSE
                  tau=obc_out
                END IF
                tau=tau*dt2d
              END IF

              IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(Iend,j-1)+                                &
     &                   grad(Iend,j  ))).gt.0.0_r8) THEN
                dVde=grad(Iend,j-1)
              ELSE
                dVde=grad(Iend,j  )
              END IF
              cff=dVdt/MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=MIN(1.0_r8,cff*dVdx)
#ifdef RADIATION_2D
              Ce=MIN(1.0_r8,MAX(-1.0_r8,cff*dVde))
#else
              Ce=0.0_r8
#endif
#if defined CELERITY_WRITE && defined FORWARD_WRITE
              BOUNDARY(ng)%vbar_east_Cx(j)=Cx
              BOUNDARY(ng)%vbar_east_Ce(j)=Ce
              BOUNDARY(ng)%vbar_east_C2(j)=cff
#endif
              vbar(Iend+1,j,kout)=(1.0_r8-Cx)*vbar(Iend+1,j,know)+      &
     &                            Cx*vbar(Iend,j,know)-                 &
     &                            MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-      &
     &                            MIN(Ce,0.0_r8)*grad(Iend+1,j  )

              IF (LBC(ieast,isVbar,ng)%nudging) THEN
                vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)+                &
     &                              tau*(BOUNDARY(ng)%vbar_east(j)-     &
     &                                   vbar(Iend+1,j,know))
              END IF
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, Chapman boundary condition.
!
        ELSE IF (LBC(ieast,isVbar,ng)%Flather.or.                       &
     &           LBC(ieast,isVbar,ng)%reduced) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pm(Iend,j-1)+                   &
     &                         GRID(ng)%pm(Iend,j  ))
              cff1=SQRT(g*0.5_r8*(GRID(ng)%h(Iend,j-1)+                 &
     &                            zeta(Iend,j-1,know)+                  &
     &                            GRID(ng)%h(Iend,j  )+                 &
     &                            zeta(Iend,j  ,know)))
              Cx=cff*cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              vbar(Iend+1,j,kout)=cff2*(vbar(Iend+1,j,know)+            &
     &                                  Cx*vbar(Iend,j,kout))
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              vbar(Iend+1,j,kout)=BOUNDARY(ng)%vbar_east(j)
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              vbar(Iend+1,j,kout)=vbar(Iend,j,kout)
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
#endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(ieast,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%east(j)) THEN
              vbar(Iend+1,j,kout)=gamma2(ng)*vbar(Iend,j,kout)
#ifdef MASKING
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
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
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
            vbar(Istr-1,Jstr,kout)=0.5_r8*(vbar(Istr  ,Jstr  ,kout)+    &
     &                                     vbar(Istr-1,Jstr+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            vbar(Iend+1,Jstr,kout)=0.5_r8*(vbar(Iend  ,Jstr  ,kout)+    &
     &                                     vbar(Iend+1,Jstr+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            vbar(Istr-1,Jend+1,kout)=0.5_r8*(vbar(Istr-1,Jend  ,kout)+  &
     &                                       vbar(Istr  ,Jend+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            vbar(Iend+1,Jend+1,kout)=0.5_r8*(vbar(Iend+1,Jend  ,kout)+  &
     &                                      vbar(Iend  ,Jend+1,kout))
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
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,j,kout))*            &
     &                    GRID(ng)%vmask_wet(Istr-1,j)
              cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,j)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              vbar(Istr,j,kout)=vbar(Istr,j,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,j,kout))*            &
     &                    GRID(ng)%vmask_wet(Iend+1,j)
              cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,j)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*cff
            END IF
          END DO
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(i,Jstr))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(i,Jstr,kout))*              &
     &                    GRID(ng)%vmask_wet(i,Jstr)
              cff=0.5_r8*GRID(ng)%vmask_wet(i,Jstr)*cff1+               &
     &            cff2*(1.0_r8-cff1)
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(i,Jend+1))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(i,Jend+1,kout))*            &
     &                    GRID(ng)%vmask_wet(i,Jend+1)
              cff=0.5_r8*GRID(ng)%vmask_wet(i,Jend+1)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*cff
            END IF
          END DO
        END IF
      END IF

      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,Jstr))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,Jstr,kout))*           &
     &                  GRID(ng)%vmask_wet(Istr-1,Jstr)
            cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,Jstr)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            vbar(Istr-1,Jstr,kout)=vbar(Istr-1,Jstr,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,Jstr))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,Jstr,kout))*           &
     &                  GRID(ng)%vmask_wet(Iend+1,Jstr)
            cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,Jstr)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            vbar(Iend+1,Jstr,kout)=vbar(Iend+1,Jstr,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,Jend+1,kout))*         &
     &                  GRID(ng)%vmask_wet(Istr-1,Jend+1)
            cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,Jend+1)*cff1+          &
     &          cff2*(1.0_r8-cff1)
            vbar(Istr-1,Jend+1,kout)=vbar(Istr-1,Jend+1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,Jend+1,kout))*         &
     &                  GRID(ng)%vmask_wet(Iend+1,Jend+1)
            cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,Jend+1)*cff1+          &
     &          cff2*(1.0_r8-cff1)
            vbar(Iend+1,Jend+1,kout)=vbar(Iend+1,Jend+1,kout)*cff
          END IF
        END IF
      END IF
#endif

      RETURN

      END SUBROUTINE v2dbc_tile
      END MODULE v2dbc_mod
