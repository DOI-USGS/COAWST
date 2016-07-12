#include "cppdefs.h"
      MODULE tl_v2dbc_mod
#ifdef TANGENT
!
!svn $Id: tl_v2dbc_im.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets tangent linear lateral boundary conditions for !
!  vertically integrated  V-velocity.  It updates the specified "kout" !
!  index.                                                              !
!                                                                      !
!  BASIC STATE variables needed: vbar                                  !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_v2dbc, tl_v2dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_v2dbc (ng, tile, kout)
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
# include "tile.h"
!
      CALL tl_v2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs(ng), kstp(ng), kout,                     &
     &                    OCEAN(ng) % ubar,                             &
     &                    OCEAN(ng) % vbar,                             &
     &                    OCEAN(ng) % zeta,                             &
     &                    OCEAN(ng) % tl_ubar,                          &
     &                    OCEAN(ng) % tl_vbar,                          &
     &                    OCEAN(ng) % tl_zeta)

      RETURN
      END SUBROUTINE tl_v2dbc
!
!***********************************************************************
      SUBROUTINE tl_v2dbc_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          krhs, kstp, kout,                       &
     &                          ubar, vbar, zeta,                       &
     &                          tl_ubar, tl_vbar, tl_zeta)
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
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)

      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j, know

      real(r8) :: Ce, Cx, Ze
      real(r8) :: bry_pgr, bry_cor, bry_str
      real(r8) :: cff, cff1, cff2, cff3, dt2d
      real(r8) :: obc_in, obc_out, tau

      real(r8) :: tl_Ce, tl_Cx, tl_Ze
      real(r8) :: tl_bry_pgr, tl_bry_cor, tl_bry_str, tl_bry_val
      real(r8) :: tl_cff, tl_cff1, tl_cff2, tl_cff3

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_grad

# include "set_bounds.h"
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
        IF (tl_LBC(isouth,isVbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO i=Istr,Iend+1
!>            grad(i,Jstr)=vbar(i  ,Jstr,know)-                         &
!>   &                     vbar(i-1,Jstr,know)
!>
              tl_grad(i,Jstr)=0.0_r8
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(isouth,isVbar,ng)%nudging) THEN
                  IF (LnudgeM2CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M2nudgcof(i,Jstr-1)+             &
     &                       CLIMA(ng)%M2nudgcof(i,Jstr  ))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M2obc_out(ng,isouth)
                    obc_in =M2obc_in (ng,isouth)
                  END IF
                  IF (BOUNDARY(ng)%vbar_south_Ce(i).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt2d
                END IF
#  ifdef RADIATION_2D
                Cx=BOUNDARY(ng)%vbar_south_Cx(i)
#  else
                Cx=0.0_r8
#  endif
                Ce=BOUNDARY(ng)%vbar_south_Ce(i)
                cff=BOUNDARY(ng)%vbar_south_C2(i)
# endif
!>              vbar(i,Jstr,kout)=(cff*vbar(i,Jstr  ,know)+             &
!>   &                             Ce *vbar(i,Jstr+1,kout)-             &
!>   &                             MAX(Cx,0.0_r8)*grad(i  ,Jstr)-       &
!>   &                             MIN(Cx,0.0_r8)*grad(i+1,Jstr))/      &
!>   &                            (cff+Ce)
!>
                tl_vbar(i,Jstr,kout)=(cff*tl_vbar(i,Jstr  ,know)+       &
     &                                Ce *tl_vbar(i,Jstr+1,kout)-       &
     &                                MAX(Cx,0.0_r8)*                   &
     &                                   tl_grad(i  ,Jstr)-             &
     &                                MIN(Cx,0.0_r8)*                   &
     &                                   tl_grad(i+1,Jstr))/            &
     &                               (cff+Ce)

                IF (tl_LBC(isouth,isVbar,ng)%nudging) THEN
!>                vbar(i,Jstr,kout)=vbar(i,Jstr,kout)+                  &
!>   &                              tau*(BOUNDARY(ng)%vbar_south(i)-    &
!>   &                                   vbar(i,Jstr,know))
!>
                  tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)-            &
     &                                 tau*tl_vbar(i,Jstr,know)
                END IF
# ifdef MASKING
!>              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                    &
!>   &                            GRID(ng)%vmask(i,Jstr)
!>
                tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*              &
     &                               GRID(ng)%vmask(i,Jstr)
# endif
              END IF
            END DO
          END IF
!
!  Southern edge, Flather boundary condition.
!
        ELSE IF (tl_LBC(isouth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (tl_LBC(isouth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(i,Jstr,know)-                          &
     &                      BOUNDARY(ng)%zeta_south(i))*                &
     &                  0.5_r8*GRID(ng)%pn(i,Jstr)
                tl_bry_pgr=-g*tl_zeta(i,Jstr,know)*                     &
     &                     0.5_r8*GRID(ng)%pn(i,Jstr)
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(isouth,isVbar,ng)) THEN
                  tl_bry_pgr=tl_bry_pgr+                                &
     &                       g*BOUNDARY(ng)%tl_zeta_south(i)*           &
     &                       0.5_r8*GRID(ng)%pn(i,Jstr)
                END IF
#  endif
              ELSE
                bry_pgr=-g*(zeta(i,Jstr  ,know)-                        &
     &                      zeta(i,Jstr-1,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                  &
     &                          GRID(ng)%pn(i,Jstr  ))
                tl_bry_pgr=-g*(tl_zeta(i,Jstr  ,know)-                  &
     &                         tl_zeta(i,Jstr-1,know))*                 &
     &                     0.5_r8*(GRID(ng)%pn(i,Jstr-1)+               &
     &                             GRID(ng)%pn(i,Jstr  ))
              END IF
#  ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jstr-1,know)+                 &
     &                           ubar(i+1,Jstr-1,know)+                 &
     &                           ubar(i  ,Jstr  ,know)+                 &
     &                           ubar(i+1,Jstr  ,know))*                &
     &                          (GRID(ng)%f(i,Jstr-1)+                  &
     &                           GRID(ng)%f(i,Jstr  ))
              tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jstr-1,know)+           &
     &                              tl_ubar(i+1,Jstr-1,know)+           &
     &                              tl_ubar(i  ,Jstr  ,know)+           &
     &                              tl_ubar(i+1,Jstr  ,know))*          &
     &                             (GRID(ng)%f(i,Jstr-1)+               &
     &                              GRID(ng)%f(i,Jstr  ))
#  else
              bry_cor=0.0_r8
              tl_bry_cor=0.0_r8
#  endif
              cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                &
     &                             zeta(i,Jstr-1,know)+                 &
     &                             GRID(ng)%h(i,Jstr  )+                &
     &                             zeta(i,Jstr  ,know)))
              tl_cff1=-cff1*cff1*(0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+      &
     &                                    tl_zeta(i,Jstr-1,know)+       &
     &                                    GRID(ng)%tl_h(i,Jstr  )+      &
     &                                    tl_zeta(i,Jstr  ,know)))
              bry_str=cff1*(FORCES(ng)%svstr(i,Jstr)-                   &
     &                      FORCES(ng)%bvstr(i,Jstr))
              tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jstr)-             &
     &                            FORCES(ng)%bvstr(i,Jstr))+            &
     &                   cff1*(FORCES(ng)%tl_svstr(i,Jstr)-             &
     &                         FORCES(ng)%tl_bvstr(i,Jstr))
              Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(i,Jstr-1)+            &
     &                                 zeta(i,Jstr-1,know)+             &
     &                                 GRID(ng)%h(i,Jstr  )+            &
     &                                 zeta(i,Jstr  ,know)))
              tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jstr-1)+       &
     &                                   tl_zeta(i,Jstr-1,know)+        &
     &                                   GRID(ng)%tl_h(i,Jstr  )+       &
     &                                   tl_zeta(i,Jstr  ,know))
              cff2=GRID(ng)%on_v(i,Jstr)*Ce
              tl_cff2=GRID(ng)%on_v(i,Jstr)*tl_Ce
!>            bry_val=vbar(i,Jstr+1,know)+                              &
!>   &                cff2*(bry_pgr+                                    &
!>   &                      bry_cor+                                    &
!>   &                      bry_str)
!>
              tl_bry_val=tl_vbar(i,Jstr+1,know)+                        &
     &                   tl_cff2*(bry_pgr+                              &
     &                            bry_cor+                              &
     &                            bry_str)+                             &
     &                   cff2*(tl_bry_pgr+                              &
     &                         tl_bry_cor+                              &
     &                         tl_bry_str)
# else
!>            bry_val=BOUNDARY(ng)%vbar_south(i)
!>
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
                tl_bry_val=BOUNDARY(ng)%tl_vbar_south(i)
              ELSE
                tl_bry_val=0.0_r8
              END IF
#  else
              tl_bry_val=0.0_r8
#  endif
# endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                 &
     &                            zeta(i,Jstr-1,know)+                  &
     &                            GRID(ng)%h(i,Jstr  )+                 &
     &                            zeta(i,Jstr  ,know)))
              tl_cff=-cff*cff*(0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+         &
     &                                 tl_zeta(i,Jstr-1,know)+          &
     &                                 GRID(ng)%tl_h(i,Jstr  )+         &
     &                                 tl_zeta(i,Jstr  ,know)))
              Ce=SQRT(g*cff)
              tl_Ce=0.5_r8*g*tl_cff/Ce
!>            vbar(i,Jstr,kout)=bry_val-                                &
!>   &                          Ce*(0.5_r8*(zeta(i,Jstr-1,know)+        &
!>   &                                      zeta(i,Jstr  ,know))-       &
!>   &                              BOUNDARY(ng)%zeta_south(i))
!>
              tl_vbar(i,Jstr,kout)=tl_bry_val-                          &
     &                             tl_Ce*                               &
     &                             (0.5_r8*(zeta(i,Jstr-1,know)+        &
     &                                      zeta(i,Jstr  ,know))-       &
     &                              BOUNDARY(ng)%zeta_south(i))-        &
     &                             Ce*                                  &
     &                             (0.5_r8*(tl_zeta(i,Jstr-1,know)+     &
     &                                      tl_zeta(i,Jstr  ,know)))
# ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
                tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)+              &
     &                               Ce*BOUNDARY(ng)%tl_zeta_south(i)
              END IF
# endif
# ifdef MASKING
!>            vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
!>   &                          GRID(ng)%vmask(i,Jstr)
!>
              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
            END IF
          END DO
!
!  Southern edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (tl_LBC(isouth,isVbar,ng)%Shchepetkin) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (tl_LBC(isouth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(i,Jstr,know)-                          &
     &                      BOUNDARY(ng)%zeta_south(i))*                &
     &                  0.5_r8*GRID(ng)%pn(i,Jstr)
                tl_bry_pgr=-g*tl_zeta(i,Jstr,know)*                     &
     &                     0.5_r8*GRID(ng)%pn(i,Jstr)
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(isouth,isVbar,ng)) THEN
                  tl_bry_pgr=tl_bry_pgr+                                &
     &                       g*BOUNDARY(ng)%tl_zeta_south(i)*           &
     &                       0.5_r8*GRID(ng)%pn(i,Jstr)
                END IF
#  endif
              ELSE
                bry_pgr=-g*(zeta(i,Jstr  ,know)-                        &
     &                      zeta(i,Jstr-1,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                  &
     &                          GRID(ng)%pn(i,Jstr  ))
                tl_bry_pgr=-g*(tl_zeta(i,Jstr  ,know)-                  &
     &                         tl_zeta(i,Jstr-1,know))*                 &
     &                     0.5_r8*(GRID(ng)%pn(i,Jstr-1)+               &
     &                             GRID(ng)%pn(i,Jstr  ))
              END IF
#  ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jstr-1,know)+                 &
     &                           ubar(i+1,Jstr-1,know)+                 &
     &                           ubar(i  ,Jstr  ,know)+                 &
     &                           ubar(i+1,Jstr  ,know))*                &
     &                          (GRID(ng)%f(i,Jstr-1)+                  &
     &                           GRID(ng)%f(i,Jstr  ))
              tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jstr-1,know)+           &
     &                              tl_ubar(i+1,Jstr-1,know)+           &
     &                              tl_ubar(i  ,Jstr  ,know)+           &
     &                              tl_ubar(i+1,Jstr  ,know))*          &
     &                             (GRID(ng)%f(i,Jstr-1)+               &
     &                              GRID(ng)%f(i,Jstr  ))
#  else
              bry_cor=0.0_r8
              tl_bry_cor=0.0_r8
#  endif
              cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                &
     &                             zeta(i,Jstr-1,know)+                 &
     &                             GRID(ng)%h(i,Jstr  )+                &
     &                             zeta(i,Jstr  ,know)))
              tl_cff1=-cff1*cff1*(0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+      &
     &                                    tl_zeta(i,Jstr-1,know)+       &
     &                                    GRID(ng)%tl_h(i,Jstr  )+      &
     &                                    tl_zeta(i,Jstr  ,know)))
              bry_str=cff1*(FORCES(ng)%svstr(i,Jstr)-                   &
     &                      FORCES(ng)%bvstr(i,Jstr))
              tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jstr)-             &
     &                            FORCES(ng)%bvstr(i,Jstr))+            &
     &                   cff1*(FORCES(ng)%tl_svstr(i,Jstr)-             &
     &                         FORCES(ng)%tl_bvstr(i,Jstr))
              Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(i,Jstr-1)+            &
     &                                 zeta(i,Jstr-1,know)+             &
     &                                 GRID(ng)%h(i,Jstr  )+            &
     &                                 zeta(i,Jstr  ,know)))
              tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jstr-1)+       &
     &                                   tl_zeta(i,Jstr-1,know)+        &
     &                                   GRID(ng)%tl_h(i,Jstr  )+       &
     &                                   tl_zeta(i,Jstr  ,know))
              cff2=GRID(ng)%on_v(i,Jstr)*Ce
              tl_cff2=GRID(ng)%on_v(i,Jstr)*tl_Ce
!>            bry_val=vbar(i,Jstr+1,know)+                              &
!>   &                cff2*(bry_pgr+                                    &
!>   &                      bry_cor+                                    &
!>   &                      bry_str)
!>
              tl_bry_val=tl_vbar(i,Jstr+1,know)+                        &
     &                   tl_cff2*(bry_pgr+                              &
     &                            bry_cor+                              &
     &                            bry_str)+                             &
     &                   cff2*(tl_bry_pgr+                              &
     &                         tl_bry_cor+                              &
     &                         tl_bry_str)
# else
!>            bry_val=BOUNDARY(ng)%vbar_south(i)
!>
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
                tl_bry_val=BOUNDARY(ng)%tl_vbar_south(i)
              ELSE
                tl_bry_val=0.0_r8
              END IF
#  else
              tl_bry_val=0.0_r8
#  endif
# endif
              cff=0.5_r8*(GRID(ng)%h(i,Jstr-1)+                         &
     &                    GRID(ng)%h(i,Jstr  ))
              tl_cff=0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+                   &
     &                       GRID(ng)%tl_h(i,Jstr  ))
              cff1=SQRT(g/cff)
              tl_cff1=-0.5_r8*cff1*tl_cff/cff
              Ce=dt2d*cff1*cff*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+           &
     &                                 GRID(ng)%pn(i,Jstr  ))
              tl_Ce=dt2d*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                 &
     &                           GRID(ng)%pn(i,Jstr  ))*                &
     &              (cff1*tl_cff+                                       &
     &               tl_cff1*cff)
              Ze=(0.5_r8+Ce)*zeta(i,Jstr  ,know)+                       &
     &           (0.5_r8-Ce)*zeta(i,Jstr-1,know)
              tl_Ze=(0.5_r8+Ce)*tl_zeta(i,Jstr  ,know)+                 &
     &              (0.5_r8-Ce)*tl_zeta(i,Jstr-1,know)+                 &
     &              tl_Ce*(zeta(i,Jstr  ,know)-                         &
     &                     zeta(i,Jstr-1,know))
              IF (Ce.gt.Co) THEN
                cff2=(1.0_r8-Co/Ce)**2
                tl_cff2=2.0_r8*cff2*Co*tl_Ce/(Ce*Ce)
                cff3=zeta(i,Jstr,kout)+                                 &
     &               Ce*zeta(i,Jstr-1,know)-                            &
     &               (1.0_r8+Ce)*zeta(i,Jstr,know)
                tl_cff3=tl_zeta(i,Jstr,kout)+                           &
     &                  Ce*tl_zeta(i,Jstr-1,know)+                      &
     &                  tl_Ce*(zeta(i,Jstr-1,know)+                     &
     &                         zeta(i,Jstr  ,know))-                    &
     &                  (1.0_r8+Ce)*tl_zeta(i,Jstr,know)
                Ze=Ze+cff2*cff3
                tl_Ze=tl_Ze+cff2*tl_cff3+                               &
     &                      tl_cff2*cff3
              END IF
!>            vbar(i,Jstr,kout)=0.5_r8*                                 &
!>   &                          ((1.0_r8-Ce)*vbar(i,Jstr,know)+         &
!>   &                           Ce*vbar(i,Jstr+1,know)+                &
!>   &                           bry_val-                               &
!>   &                           cff1*(Ze-BOUNDARY(ng)%zeta_south(i)))
!>
              tl_vbar(i,Jstr,kout)=0.5_r8*                              &
     &                             ((1.0_r8-Ce)*                        &
     &                              tl_vbar(i,Jstr,know)-               &
     &                              tl_Ce*(vbar(i,Jstr  ,know)-         &
     &                                     vbar(i,Jstr+1,know))+        &
     &                              Ce*tl_vbar(i,Jstr+1,know)+          &
     &                              tl_bry_val-                         &
     &                              tl_cff1*                            &
     &                              (Ze-BOUNDARY(ng)%zeta_south(i))-    &
     &                              cff1*tl_Ze)
# ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
                tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)+              &
     &                               0.5_r8*cff1*                       &
     &                               BOUNDARY(ng)%tl_zeta_south(i)
              END IF
# endif
# ifdef MASKING
!>            vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
!>   &                          GRID(ng)%vmask(i,Jstr)
!>
              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(isouth,isVbar,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            vbar(i,Jstr,kout)=BOUNDARY(ng)%vbar_south(i)
!>
# ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
                tl_vbar(i,Jstr,kout)=BOUNDARY(ng)%tl_vbar_south(i)
              ELSE
                tl_vbar(i,Jstr,kout)=0.0_r8
              END IF
# else
              tl_vbar(i,Jstr,kout)=0.0_r8
# endif
# ifdef MASKING
!>            vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
!>   &                          GRID(ng)%vmask(i,Jstr)
!>
              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(isouth,isVbar,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            vbar(i,Jstr,kout)=vbar(i,Jstr+1,kout)
!>
              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr+1,kout)
# ifdef MASKING
!>            vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
!>   &                          GRID(ng)%vmask(i,Jstr)
!>
              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
            END IF
          END DO
!
!  Southern edge, reduced-physics boundary condition.
!
        ELSE IF (tl_LBC(isouth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              IF (tl_LBC(isouth,isFsur,ng)%acquire) THEN
!>              bry_pgr=-g*(zeta(i,Jstr,know)-                          &
!>   &                      BOUNDARY(ng)%zeta_south(i))*                &
!>   &                  0.5_r8*GRID(ng)%pn(i,Jstr)
!>
                tl_bry_pgr=-g*tl_zeta(i,Jstr,know)*                     &
     &                     0.5_r8*GRID(ng)%pn(i,Jstr)
# ifdef ADJUST_BOUNDARY
                IF (Lobc(inorth,isVbar,ng)) THEN
                  tl_bry_pgr=tl_bry_pgr+                                &
     &                       g*BOUNDARY(ng)%tl_zeta_south(i)*           &
     &                       0.5_r8*GRID(ng)%pn(i,Jstr)
                END IF
# endif
              ELSE
!>              bry_pgr=-g*(zeta(i,Jstr,know)-                          &
!>   &                      zeta(i,Jstr-1,know))*                       &
!>   &                  0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                  &
!>   &                          GRID(ng)%pn(i,Jstr  ))
!>
                tl_bry_pgr=-g*(tl_zeta(i,Jstr  ,know)-                  &
     &                         tl_zeta(i,Jstr-1,know))*                 &
     &                     0.5_r8*(GRID(ng)%pn(i,Jstr-1)+               &
     &                             GRID(ng)%pn(i,Jstr  ))
              END IF
# ifdef UV_COR
!>            bry_cor=-0.125_r8*(ubar(i  ,Jstr-1,know)+                 &
!>   &                           ubar(i+1,Jstr-1,know)+                 &
!>   &                           ubar(i  ,Jstr  ,know)+                 &
!>   &                           ubar(i+1,Jstr  ,know))*                &
!>   &                          (GRID(ng)%f(i,Jstr-1)+                  &
!>   &                           GRID(ng)%f(i,Jstr  ))
!>
              tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jstr-1,know)+           &
     &                              tl_ubar(i+1,Jstr-1,know)+           &
     &                              tl_ubar(i  ,Jstr  ,know)+           &
     &                              tl_ubar(i+1,Jstr  ,know))*          &
     &                             (GRID(ng)%f(i,Jstr-1)+               &
     &                              GRID(ng)%f(i,Jstr  ))
# else
!>            bry_cor=0.0_r8
!>
              tl_bry_cor=0.0_r8
# endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                 &
     &                            zeta(i,Jstr-1,know)+                  &
     &                            GRID(ng)%h(i,Jstr  )+                 &
     &                            zeta(i,Jstr  ,know)))
              tl_cff=-cff*cff*0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+          &
     &                                tl_zeta(i,Jstr-1,know)+           &
     &                                GRID(ng)%tl_h(i,Jstr  )+          &
     &                                tl_zeta(i,Jstr  ,know))
!>            bry_str=cff*(FORCES(ng)%svstr(i,Jstr)-                    &
!>   &                     FORCES(ng)%bvstr(i,Jstr))
!>
              tl_bry_str=tl_cff*(FORCES(ng)%svstr(i,Jstr)-              &
     &                           FORCES(ng)%bvstr(i,Jstr))+             &
     &                   cff*(FORCES(ng)%tl_svstr(i,Jstr)-              &
     &                        FORCES(ng)%tl_bvstr(i,Jstr))

!>            vbar(i,Jstr,kout)=vbar(i,Jstr,know)+                      &
!>   &                          dt2d*(bry_pgr+                          &
!>   &                                bry_cor+                          &
!>   &                                bry_str)
!>
              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,know)+                &
     &                             dt2d*(tl_bry_pgr+                    &
     &                                   tl_bry_cor+                    &
     &                                   tl_bry_str)
# ifdef MASKING
!>            vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
!>   &                          GRID(ng)%vmask(i,Jstr)
!>
              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(isouth,isVbar,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            vbar(i,Jstr,kout)=0.0_r8
!>
              tl_vbar(i,Jstr,kout)=0.0_r8
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
        IF (tl_LBC(inorth,isVbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO i=Istr,Iend+1
!>            grad(i,Jend+1)=vbar(i  ,Jend+1,know)-                     &
!>   &                       vbar(i-1,Jend+1,know)
!>
              tl_grad(i,Jend+1)=0.0_r8
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(inorth,isVbar,ng)%nudging) THEN
                  IF (LnudgeM2CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M2nudgcof(i,Jend  )+             &
     &                       CLIMA(ng)%M2nudgcof(i,Jend+1))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M2obc_out(ng,inorth)
                    obc_in =M2obc_in (ng,inorth)
                  END IF
                  IF (BOUNDARY(ng)%vbar_north_Ce(i).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt2d
                END IF
#  ifdef RADIATION_2D
                Cx=BOUNDARY(ng)%vbar_north_Cx(i)
#  else
                Cx=0.0_r8
#  endif
                Ce=BOUNDARY(ng)%vbar_north_Ce(i)
                cff=BOUNDARY(ng)%vbar_north_C2(i)
# endif
!>              vbar(i,Jend+1,kout)=(cff*vbar(i,Jend+1,know)+           &
!>   &                               Ce *vbar(i,Jend  ,kout)-           &
!>   &                               MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-   &
!>   &                               MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/  &
!>   &                              (cff+Ce)
!>
                tl_vbar(i,Jend+1,kout)=(cff*tl_vbar(i,Jend+1,know)+     &
     &                                  Ce *tl_vbar(i,Jend  ,kout)-     &
     &                                  MAX(Cx,0.0_r8)*                 &
     &                                     tl_grad(i  ,Jend+1)-         &
     &                                  MIN(Cx,0.0_r8)*                 &
     &                                     tl_grad(i+1,Jend+1))/        &
     &                                 (cff+Ce)

                IF (tl_LBC(inorth,isVbar,ng)%nudging) THEN
!>                vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)+              &
!>   &                                tau*(BOUNDARY(ng)%vbar_north(i)-  &
!>   &                                     vbar(i,Jend+1,know))
!>
                  tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)-        &
     &                                   tau*tl_vbar(i,Jend+1,know)
                END IF
# ifdef MASKING
!>              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                &
!>   &                              GRID(ng)%vmask(i,Jend+1)
!>
                tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*          &
     &                                 GRID(ng)%vmask(i,Jend+1)
# endif
              END IF
            END DO
          END IF
!
!  Northern edge, Flather boundary condition.
!
        ELSE IF (tl_LBC(inorth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (tl_LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                 &
     &                      zeta(i,Jend,know))*                         &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
                tl_bry_pgr=g*tl_zeta(i,Jend,know)*                      &
     &                     0.5_r8*GRID(ng)%pn(i,Jend)
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(inorth,isVbar,ng)) THEN
                  tl_bry_pgr=tl_bry_pgr-                                &
     &                       g*BOUNDARY(ng)%tl_zeta_north(i)*           &
     &                       0.5_r8*GRID(ng)%pn(i,Jend)
                END IF
#  endif
              ELSE
                bry_pgr=-g*(zeta(i,Jend+1,know)-                        &
     &                      zeta(i,Jend  ,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
                tl_bry_pgr=-g*(tl_zeta(i,Jend+1,know)-                  &
     &                         tl_zeta(i,Jend  ,know))*                 &
     &                     0.5_r8*(GRID(ng)%pn(i,Jend  )+               &
     &                             GRID(ng)%pn(i,Jend+1))
              END IF
#  ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jend  ,know)+                 &
     &                           ubar(i+1,Jend  ,know)+                 &
     &                           ubar(i  ,Jend+1,know)+                 &
     &                           ubar(i+1,Jend+1,know))*                &
     &                          (GRID(ng)%f(i,Jend  )+                  &
     &                           GRID(ng)%f(i,Jend+1))
              tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jend  ,know)+           &
     &                              tl_ubar(i+1,Jend  ,know)+           &
     &                              tl_ubar(i  ,Jend+1,know)+           &
     &                              tl_ubar(i+1,Jend+1,know))*          &
     &                             (GRID(ng)%f(i,Jend  )+               &
     &                              GRID(ng)%f(i,Jend+1))
#  else
              bry_cor=0.0_r8
              tl_bry_cor=0.0_r8
#  endif
              cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                &
     &                             zeta(i,Jend  ,know)+                 &
     &                             GRID(ng)%h(i,Jend+1)+                &
     &                             zeta(i,Jend+1,know)))
              tl_cff1=-cff1*cff1*0.5_r8*(GRID(ng)%tl_h(i,Jend  )+       &
     &                                   tl_zeta(i,Jend  ,know)+        &
     &                                   GRID(ng)%tl_h(i,Jend+1)+       &
     &                                   tl_zeta(i,Jend+1,know))
              bry_str=cff1*(FORCES(ng)%svstr(i,Jend+1)-                 &
     &                      FORCES(ng)%bvstr(i,Jend+1))
              tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jend+1)-           &
     &                            FORCES(ng)%bvstr(i,Jend+1))+          &
     &                   cff1*(FORCES(ng)%tl_svstr(i,Jend+1)-           &
     &                         FORCES(ng)%tl_bvstr(i,Jend+1))
              Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(i,Jend+1)+            &
     &                                 zeta(i,Jend+1,know)+             &
     &                                 GRID(ng)%h(i,Jend  )+            &
     &                                 zeta(i,Jend  ,know)))
              tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jend+1)+       &
     &                                   tl_zeta(i,Jend+1,know)+        &
     &                                   GRID(ng)%tl_h(i,Jend  )+       &
     &                                   tl_zeta(i,Jend  ,know))
              cff2=GRID(ng)%on_v(i,Jend+1)*Ce
              tl_cff2=GRID(ng)%on_v(i,Jend+1)*tl_Ce
!>            bry_val=vbar(i,Jend,know)+                                &
!>   &                cff2*(bry_pgr+                                    &
!>   &                      bry_cor+                                    &
!>   &                      bry_str)
!>
              tl_bry_val=tl_vbar(i,Jend,know)+                          &
     &                   tl_cff2*(bry_pgr+                              &
     &                            bry_cor+                              &
     &                            bry_str)+                             &
     &                   cff2*(tl_bry_pgr+                              &
     &                         tl_bry_cor+                              &
     &                         tl_bry_str)
# else
!>            bry_val=BOUNDARY(ng)%vbar_north(i)
!>
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
                tl_bry_val=BOUNDARY(ng)%tl_vbar_north(i)
              ELSE
                tl_bry_val=0.0_r8
              END IF
#  else
              tl_bry_val=0.0_r8
#  endif
# endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                 &
     &                            zeta(i,Jend  ,know)+                  &
     &                            GRID(ng)%h(i,Jend+1)+                 &
     &                            zeta(i,Jend+1,know)))
              tl_cff=-cff*cff*(0.5_r8*(GRID(ng)%tl_h(i,Jend  )+         &
     &                                 tl_zeta(i,Jend  ,know)+          &
     &                                 GRID(ng)%tl_h(i,Jend+1)+         &
     &                                 tl_zeta(i,Jend+1,know)))
              Ce=SQRT(g*cff)
              tl_Ce=0.5_r8*g*tl_cff/Ce
!>            vbar(i,Jend+1,kout)=bry_val+                              &
!>   &                            Ce*(0.5_r8*(zeta(i,Jend  ,know)+      &
!>   &                                        zeta(i,Jend+1,know))-     &
!>   &                                BOUNDARY(ng)%zeta_north(i))
!>
              tl_vbar(i,Jend+1,kout)=tl_bry_val+                        &
     &                               tl_Ce*                             &
     &                               (0.5_r8*(zeta(i,Jend  ,know)+      &
     &                                        zeta(i,Jend+1,know))-     &
     &                                BOUNDARY(ng)%zeta_north(i))+      &
     &                               Ce*                                &
     &                               (0.5_r8*(tl_zeta(i,Jend  ,know)+   &
     &                                        tl_zeta(i,Jend+1,know)))
# ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
                tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)-          &
     &                                 Ce*BOUNDARY(ng)%tl_zeta_north(i)
              END IF
# endif
# ifdef MASKING
!>            vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
!>   &                            GRID(ng)%vmask(i,Jend+1)
!>
              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (LBC(inorth,isVbar,ng)%Shchepetkin) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (tl_LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                 &
     &                      zeta(i,Jend,know))*                         &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
                tl_bry_pgr=g*tl_zeta(i,Jend,know)*                      &
     &                     0.5_r8*GRID(ng)%pn(i,Jend)
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(inorth,isVbar,ng)) THEN
                  tl_bry_pgr=tl_bry_pgr-                                &
     &                       g*BOUNDARY(ng)%tl_zeta_north(i)*           &
     &                       0.5_r8*GRID(ng)%pn(i,Jend)
                END IF
#  endif
              ELSE
                bry_pgr=-g*(zeta(i,Jend+1,know)-                        &
     &                      zeta(i,Jend  ,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
                tl_bry_pgr=-g*(tl_zeta(i,Jend+1,know)-                  &
     &                         tl_zeta(i,Jend  ,know))*                 &
     &                     0.5_r8*(GRID(ng)%pn(i,Jend  )+               &
     &                             GRID(ng)%pn(i,Jend+1))
              END IF
#  ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jend  ,know)+                 &
     &                           ubar(i+1,Jend  ,know)+                 &
     &                           ubar(i  ,Jend+1,know)+                 &
     &                           ubar(i+1,Jend+1,know))*                &
     &                          (GRID(ng)%f(i,Jend  )+                  &
     &                           GRID(ng)%f(i,Jend+1))
              tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jend  ,know)+           &
     &                              tl_ubar(i+1,Jend  ,know)+           &
     &                              tl_ubar(i  ,Jend+1,know)+           &
     &                              tl_ubar(i+1,Jend+1,know))*          &
     &                             (GRID(ng)%f(i,Jend  )+               &
     &                              GRID(ng)%f(i,Jend+1))
#  else
              bry_cor=0.0_r8
              tl_bry_cor=0.0_r8
#  endif
              cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                &
     &                             zeta(i,Jend  ,know)+                 &
     &                             GRID(ng)%h(i,Jend+1)+                &
     &                             zeta(i,Jend+1,know)))
              tl_cff1=-cff1*cff1*0.5_r8*(GRID(ng)%tl_h(i,Jend  )+       &
     &                                   tl_zeta(i,Jend  ,know)+        &
     &                                   GRID(ng)%tl_h(i,Jend+1)+       &
     &                                   tl_zeta(i,Jend+1,know))
              bry_str=cff1*(FORCES(ng)%svstr(i,Jend+1)-                 &
     &                      FORCES(ng)%bvstr(i,Jend+1))
              tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jend+1)-           &
     &                            FORCES(ng)%bvstr(i,Jend+1))+          &
     &                   cff1*(FORCES(ng)%tl_svstr(i,Jend+1)-           &
     &                         FORCES(ng)%tl_bvstr(i,Jend+1))
              Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(i,Jend+1)+            &
     &                                 zeta(i,Jend+1,know)+             &
     &                                 GRID(ng)%h(i,Jend  )+            &
     &                                 zeta(i,Jend  ,know)))
              tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jend+1)+       &
     &                                   tl_zeta(i,Jend+1,know)+        &
     &                                   GRID(ng)%tl_h(i,Jend  )+       &
     &                                   tl_zeta(i,Jend  ,know))
              cff2=GRID(ng)%on_v(i,Jend+1)*Ce
              tl_cff2=GRID(ng)%on_v(i,Jend+1)*tl_Ce
!>            bry_val=vbar(i,Jend,know)+                                &
!>   &                cff2*(bry_pgr+                                    &
!>   &                      bry_cor+                                    &
!>   &                      bry_str)
!>
              tl_bry_val=tl_vbar(i,Jend,know)+                          &
     &                   tl_cff2*(bry_pgr+                              &
     &                            bry_cor+                              &
     &                            bry_str)+                             &
     &                   cff2*(tl_bry_pgr+                              &
     &                         tl_bry_cor+                              &
     &                         tl_bry_str)
# else
!>            bry_val=BOUNDARY(ng)%vbar_north(i)
!>
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
                tl_bry_val=BOUNDARY(ng)%tl_vbar_north(i)
              ELSE
                tl_bry_val=0.0_r8
              END IF
#  else
              tl_bry_val=0.0_r8
#  endif
# endif
              cff=0.5_r8*(GRID(ng)%h(i,Jend  )+                         &
     &                    GRID(ng)%h(i,Jend+1))
              tl_cff=0.5_r8*(GRID(ng)%tl_h(i,Jend  )+                   &
     &                       GRID(ng)%tl_h(i,Jend+1))
              cff1=SQRT(g/cff)
              tl_cff1=-0.5_r8*cff1*tl_cff/cff
              Ce=dt2d*cff1*cff*0.5_r8*(GRID(ng)%pn(i,Jend  )+           &
     &                                 GRID(ng)%pn(i,Jend+1))
              tl_Ce=dt2d*0.5_r8*(GRID(ng)%pn(i,Jend  )+                 &
     &                           GRID(ng)%pn(i,Jend+1))*                &
     &              (cff1*tl_cff+                                       &
     &               tl_cff1*cff)
              Ze=(0.5_r8+Ce)*zeta(i,Jend  ,know)+                       &
     &           (0.5_r8-Ce)*zeta(i,Jend+1,know)
              tl_Ze=(0.5_r8+Ce)*tl_zeta(i,Jend  ,know)+                 &
     &              (0.5_r8-Ce)*tl_zeta(i,Jend+1,know)+                 &
     &              tl_Ce*(zeta(i,Jend  ,know)-                         &
     &                     zeta(i,Jend+1,know))
              IF (Ce.gt.Co) THEN
                cff2=(1.0_r8-Co/Ce)**2
                tl_cff2=2.0_r8*cff2*Co*tl_Ce/(Ce*Ce)
                cff3=zeta(i,Jend,kout)+                                 &
     &               Ce*zeta(i,Jend+1,know)-                            &
     &               (1.0_r8+Ce)*zeta(i,Jend,know)
                tl_cff3=tl_zeta(i,Jend,kout)+                           &
     &                  Ce*tl_zeta(i,Jend+1,know)+                      &
     &                  tl_Ce*(zeta(i,Jend  ,know)+                     &
     &                         zeta(i,Jend+1,know))-                    &
     &                  (1.0_r8+Ce)*tl_zeta(i,Jend,know)
                Ze=Ze+cff2*cff3
                tl_Ze=tl_Ze+cff2*tl_cff3+                               &
     &                      tl_cff2*cff3
              END IF
!>            vbar(i,Jend+1,kout)=0.5_r8*                               &
!>   &                            ((1.0_r8-Ce)*vbar(i,Jend+1,know)+     &
!>   &                             Ce*vbar(i,Jend,know)+                &
!>   &                             bry_val+                             &
!>   &                             cff1*(Ze-BOUNDARY(ng)%zeta_north(i)))
!>
              tl_vbar(i,Jend+1,kout)=0.5_r8*                            &
     &                               ((1.0_r8-Ce)*                      &
     &                                tl_vbar(i,Jend+1,know)+           &
     &                                tl_Ce*(vbar(i,Jend  ,know)-       &
     &                                       vbar(i,Jend+1,know))+      &
     &                                Ce*tl_vbar(i,Jend,know)+          &
     &                                tl_bry_val+                       &
     &                                tl_cff1*                          &
     &                                (Ze-BOUNDARY(ng)%zeta_north(i))-  &
     &                                 cff1*tl_ze)
# ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
                tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)-          &
     &                                 0.5_r8*cff1*                     &
     &                                 Ce*BOUNDARY(ng)%tl_zeta_north(i)
              END IF
# endif
# ifdef MASKING
!>            vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
!>   &                            GRID(ng)%vmask(i,Jend+1)
!>
              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(inorth,isVbar,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            vbar(i,Jend+1,kout)=BOUNDARY(ng)%vbar_north(i)
!>
# ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
                tl_vbar(i,Jend+1,kout)=BOUNDARY(ng)%tl_vbar_north(i)
              ELSE
                tl_vbar(i,Jend+1,kout)=0.0_r8
              END IF
# else
              tl_vbar(i,Jend+1,kout)=0.0_r8
# endif
# ifdef MASKING
!>            vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
!>   &                            GRID(ng)%vmask(i,Jend+1)
!>
              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(inorth,isVbar,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            vbar(i,Jend+1,kout)=vbar(i,Jend,kout)
!>
              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend,kout)
# ifdef MASKING
!>            vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
!>   &                            GRID(ng)%vmask(i,Jend+1)
!>
              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, reduced-physics boundary condition.
!
        ELSE IF (tl_LBC(inorth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              IF (tl_LBC(inorth,isFsur,ng)%acquire) THEN
!>              bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                 &
!>   &                      zeta(i,Jend,know))*                         &
!>   &                  0.5_r8*GRID(ng)%pn(i,Jend)
!>
                tl_bry_pgr=g*tl_zeta(i,Jend,know)*                      &
     &                     0.5_r8*GRID(ng)%pn(i,Jend)
# ifdef ADJUST_BOUNDARY
                IF (Lobc(inorth,isVbar,ng)) THEN
                  tl_bry_pgr=tl_bry_pgr-                                &
     &                       g*BOUNDARY(ng)%tl_zeta_north(i)*           &
     &                       0.5_r8*GRID(ng)%pn(i,Jend)
                END IF
# endif
              ELSE
!>              bry_pgr=-g*(zeta(i,Jend+1,know)-                        &
!>   &                      zeta(i,Jend  ,know))*                       &
!>   &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
!>   &                          GRID(ng)%pn(i,Jend+1))
!>
                tl_bry_pgr=-g*(tl_zeta(i,Jend+1,know)-                  &
     &                         tl_zeta(i,Jend  ,know))*                 &
     &                     0.5_r8*(GRID(ng)%pn(i,Jend  )+               &
     &                             GRID(ng)%pn(i,Jend+1))
              END IF
# ifdef UV_COR
!>            bry_cor=-0.125_r8*(ubar(i  ,Jend  ,know)+                 &
!>   &                           ubar(i+1,Jend  ,know)+                 &
!>   &                           ubar(i  ,Jend+1,know)+                 &
!>   &                           ubar(i+1,Jend+1,know))*                &
!>   &                          (GRID(ng)%f(i,Jend  )+                  &
!>   &                           GRID(ng)%f(i,Jend+1))
!>
              tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jend  ,know)+           &
     &                              tl_ubar(i+1,Jend  ,know)+           &
     &                              tl_ubar(i  ,Jend+1,know)+           &
     &                              tl_ubar(i+1,Jend+1,know))*          &
     &                             (GRID(ng)%f(i,Jend  )+               &
     &                              GRID(ng)%f(i,Jend+1))
# else
!>            bry_cor=0.0_r8
!>
              tl_bry_cor=0.0_r8
# endif
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                 &
     &                            zeta(i,Jend  ,know)+                  &
     &                            GRID(ng)%h(i,Jend+1)+                 &
     &                            zeta(i,Jend+1,know)))
              tl_cff=-cff*cff*0.5_r8*(GRID(ng)%tl_h(i,Jend  )+          &
     &                                tl_zeta(i,Jend  ,know)+           &
     &                                GRID(ng)%tl_h(i,Jend+1)+          &
     &                                tl_zeta(i,Jend+1,know))
!>            bry_str=cff*(FORCES(ng)%svstr(i,Jend+1)-                  &
!>   &                     FORCES(ng)%bvstr(i,Jend+1))
!>
              tl_bry_str=tl_cff*(FORCES(ng)%svstr(i,Jend+1)-            &
     &                           FORCES(ng)%bvstr(i,Jend+1))+           &
     &                   cff*(FORCES(ng)%tl_svstr(i,Jend+1)-            &
     &                        FORCES(ng)%tl_bvstr(i,Jend+1))
!>            vbar(i,Jend+1,kout)=vbar(i,Jend+1,know)+                  &
!>   &                            dt2d*(bry_pgr+                        &
!>   &                                  bry_cor+                        &
!>   &                                  bry_str)
!>
              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,know)+            &
     &                               dt2d*(tl_bry_pgr+                  &
     &                                     tl_bry_cor+                  &
     &                                     tl_bry_str)
# ifdef MASKING
!>            vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
!>   &                            GRID(ng)%vmask(i,Jend+1)
!>
              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (tl_LBC(inorth,isVbar,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            vbar(i,Jend+1,kout)=0.0_r8
!>
              tl_vbar(i,Jend+1,kout)=0.0_r8
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
        IF (tl_LBC(iwest,isVbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO j=JstrV-1,Jend
!>            grad(Istr-1,j)=vbar(Istr-1,j+1,know)-                     &
!>   &                       vbar(Istr-1,j  ,know)
!>
              tl_grad(Istr-1,j)=0.0_r8
            END DO
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(iwest,isVbar,ng)%nudging) THEN
                  IF (LnudgeM2CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M2nudgcof(Istr-1,j-1)+           &
     &                       CLIMA(ng)%M2nudgcof(Istr-1,j  ))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M2obc_out(ng,iwest)
                    obc_in =M2obc_in (ng,iwest)
                  END IF
                  IF (BOUNDARY(ng)%vbar_west_Cx(j).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt2d
                END IF
                Cx=BOUNDARY(ng)%vbar_west_Cx(j)
#  ifdef RADIATION_2D
                Ce=BOUNDARY(ng)%vbar_west_Ce(j)
#  else
                Ce=0.0_r8
#  endif
                cff=BOUNDARY(ng)%vbar_west_C2(j)
# endif
!>              vbar(Istr-1,j,kout)=(cff*vbar(Istr-1,j,know)+           &
!>   &                               Cx *vbar(Istr  ,j,kout)-           &
!>   &                               MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-   &
!>   &                               MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/  &
!>   &                              (cff+Cx)
!>
                tl_vbar(Istr-1,j,kout)=(cff*tl_vbar(Istr-1,j,know)+     &
     &                                  Cx *tl_vbar(1,j,kout)-          &
     &                                  MAX(Ce,0.0_r8)*                 &
     &                                     tl_grad(Istr-1,j-1)-         &
     &                                  MIN(Ce,0.0_r8)*                 &
     &                                     tl_grad(Istr-1,j  ))/        &
     &                                 (cff+Cx)

                IF (tl_LBC(iwest,isVbar,ng)%nudging) THEN
!>                vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)+              &
!>   &                                tau*(BOUNDARY(ng)%vbar_west(j)-   &
!>   &                                     vbar(Istr-1,j,know))
!>
                  tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)-        &
     &                                   tau*tl_vbar(1,j,know)
                END IF
# ifdef MASKING
!>              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                &
!>   &                              GRID(ng)%vmask(Istr-1,j)
!>
                tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*          &
     &                                 GRID(ng)%vmask(Istr-1,j)
# endif
              END IF
            END DO
          END IF
!
!  Western edge, Chapman boundary condition.
!
        ELSE IF (tl_LBC(iwest,isVbar,ng)%Flather.or.                    &
     &           tl_LBC(iwest,isVbar,ng)%reduced.or.                    &
     &           tl_LBC(iwest,isVbar,ng)%Shchepetkin) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pm(Istr,j-1)+                   &
     &                         GRID(ng)%pm(Istr,j  ))
              cff1=SQRT(g*0.5_r8*(GRID(ng)%h(Istr,j-1)+                 &
     &                            zeta(Istr,j-1,know)+                  &
     &                            GRID(ng)%h(Istr,j  )+                 &
     &                            zeta(Istr,j  ,know)))
              tl_cff1=0.25_r8*g*(GRID(ng)%tl_h(Istr,j-1)+               &
     &                           tl_zeta(Istr,j-1,know)+                &
     &                           GRID(ng)%tl_h(Istr,j  )+               &
     &                           tl_zeta(Istr,j  ,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              tl_cff2=-cff2*cff2*tl_Cx
!>            vbar(Istr-1,j,kout)=cff2*(vbar(Istr-1,j,know)+            &
!>   &                                  Cx*vbar(Istr,j,kout))
!>
              tl_vbar(Istr-1,j,kout)=tl_cff2*(vbar(Istr-1,j,know)+      &
     &                                        Cx*vbar(Istr,j,kout))+    &
     &                               cff2*(tl_vbar(Istr-1,j,know)+      &
     &                                     tl_Cx*vbar(Istr,j,kout)+     &
     &                                     Cx*tl_vbar(Istr,j,kout))
# ifdef MASKING
!>            vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
!>   &                            GRID(ng)%vmask(Istr-1,j)
!>
              tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*            &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(iwest,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!>            vbar(Istr-1,j,kout)=BOUNDARY(ng)%vbar_west(j)
!>
# ifdef ADJUST_BOUNDARY
              IF (Lobc(iwest,isVbar,ng)) THEN
                tl_vbar(Istr-1,j,kout)=BOUNDARY(ng)%tl_vbar_west(j)
              ELSE
                tl_vbar(Istr-1,j,kout)=0.0_r8
              END IF
# else
              tl_vbar(Istr-1,j,kout)=0.0_r8
# endif
# ifdef MASKING
!>            vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
!>   &                            GRID(ng)%vmask(Istr-1,j)
!>
              tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*            &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(iwest,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!>            vbar(Istr-1,j,kout)=vbar(Istr,j,kout)
!>
              tl_vbar(Istr-1,j,kout)=tl_vbar(Istr,j,kout)
# ifdef MASKING
!>            vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
!>   &                            GRID(ng)%vmask(Istr-1,j)
!>
              tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*            &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
            END IF
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (tl_LBC(iwest,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%west(j)) THEN
!>            vbar(Istr-1,j,kout)=gamma2(ng)*vbar(Istr,j,kout)
!>
              tl_vbar(Istr-1,j,kout)=gamma2(ng)*tl_vbar(Istr,j,kout)
# ifdef MASKING
!>            vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
!>   &                            GRID(ng)%vmask(Istr-1,j)
!>
              tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*            &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
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
        IF (tl_LBC(ieast,isVbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO j=JstrV-1,Jend
!>            grad(Iend+1,j)=vbar(Iend+1,j+1,know)-                     &
!>   &                       vbar(Iend+1,j  ,know)
!>
              tl_grad(Iend+1,j)=0.0_r8
            END DO
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (tl_LBC(ieast,isVbar,ng)%nudging) THEN
                  IF (LnudgeM2CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M2nudgcof(Iend+1,j-1)+           &
     &                       CLIMA(ng)%M2nudgcof(Iend+1,j  ))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M2obc_out(ng,ieast)
                    obc_in =M2obc_in (ng,ieast)
                  END IF
                  IF (BOUNDARY(ng)%vbar_east_Cx(j).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt2d
                END IF
                Cx=BOUNDARY(ng)%vbar_east_Cx(j)
#  ifdef RADIATION_2D
                Ce=BOUNDARY(ng)%vbar_east_Ce(j)
#  else
                Ce=0.0_r8
#  endif
                cff=BOUNDARY(ng)%vbar_east_C2(j)
# endif
!>              vbar(Iend+1,j,kout)=(cff*vbar(Iend+1,j,know)+           &
!>   &                               Cx *vbar(Iend  ,j,kout)-           &
!>   &                               MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-   &
!>   &                               MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/  &
!>   &                              (cff+Cx)
!>
                tl_vbar(Iend+1,j,kout)=(cff*tl_vbar(Iend+1,j,know)+     &
     &                                  Cx *tl_vbar(Iend  ,j,kout)-     &
     &                                  MAX(Ce,0.0_r8)*                 &
     &                                     tl_grad(Iend+1,j-1)-         &
     &                                  MIN(Ce,0.0_r8)*                 &
     &                                     tl_grad(Iend+1,j  ))/        &
     &                                 (cff+Cx)

                IF (tl_LBC(ieast,isVbar,ng)%nudging) THEN
!>                vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)+              &
!>   &                                tau*(BOUNDARY(ng)%vbar_east(j)-   &
!>   &                                     vbar(Iend+1,j,know))
!>
                  tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)-        &
     &                                   tau*tl_vbar(Iend+1,j,know)
                END IF
# ifdef MASKING
!>              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                &
!>   &                              GRID(ng)%vmask(Iend+1,j)
!>
                tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*          &
     &                                 GRID(ng)%vmask(Iend+1,j)
# endif
              END IF
            END DO
          END IF
!
!  Eastern edge, Chapman boundary condition.
!
        ELSE IF (tl_LBC(ieast,isVbar,ng)%Flather.or.                    &
     &           tl_LBC(ieast,isVbar,ng)%reduced.or.                    &
     &           tl_LBC(ieast,isVbar,ng)%Shchepetkin) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pm(Iend,j-1)+                   &
     &                         GRID(ng)%pm(Iend,j  ))
              cff1=SQRT(g*0.5_r8*(GRID(ng)%h(Iend,j-1)+                 &
     &                            zeta(Iend,j-1,know)+                  &
     &                            GRID(ng)%h(Iend,j  )+                 &
     &                            zeta(Iend,j  ,know)))
              tl_cff1=0.25_r8*g*(GRID(ng)%tl_h(Iend,j-1)+               &
     &                           tl_zeta(Iend,j-1,know)+                &
     &                           GRID(ng)%tl_h(Iend,j  )+               &
     &                           tl_zeta(Iend,j  ,know))/cff1
              Cx=cff*cff1
              tl_Cx=cff*tl_cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              tl_cff2=-cff2*cff2*tl_Cx
!>            vbar(Iend+1,j,kout)=cff2*(vbar(Iend+1,j,know)+            &
!>   &                                  Cx*vbar(Iend,j,kout))
!>
              tl_vbar(Iend+1,j,kout)=tl_cff2*(vbar(Iend+1,j,know)+      &
     &                                        Cx*vbar(Iend,j,kout))+    &
     &                               cff2*(tl_vbar(Iend+1,j,know)+      &
     &                                     tl_Cx*vbar(Iend,j,kout)+     &
     &                                     Cx*tl_vbar(Iend,j,kout))
# ifdef MASKING
!>            vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
!>   &                            GRID(ng)%vmask(Iend+1,j)
!>
              tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*            &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (tl_LBC(ieast,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!>            vbar(Iend+1,j,kout)=BOUNDARY(ng)%vbar_east(j)
!>
# ifdef ADJUST_BOUNDARY
              IF (Lobc(ieast,isVbar,ng)) THEN
                tl_vbar(Iend+1,j,kout)=BOUNDARY(ng)%tl_vbar_east(j)
              ELSE
                tl_vbar(Iend+1,j,kout)=0.0_r8
              END IF
# else
              tl_vbar(Iend+1,j,kout)=0.0_r8
# endif
# ifdef MASKING
!>            vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
!>   &                            GRID(ng)%vmask(Iend+1,j)
!>
              tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*            &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (tl_LBC(ieast,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!>            vbar(Iend+1,j,kout)=vbar(Iend,j,kout)
!>
              tl_vbar(Iend+1,j,kout)=tl_vbar(Iend,j,kout)
# ifdef MASKING
!>            vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
!>   &                            GRID(ng)%vmask(Iend+1,j)
!>
              tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*            &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (tl_LBC(ieast,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%east(j)) THEN
!>            vbar(Iend+1,j,kout)=gamma2(ng)*vbar(Iend,j,kout)
!>
              tl_vbar(Iend+1,j,kout)=gamma2(ng)*tl_vbar(Iend,j,kout)
# ifdef MASKING
!>            vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
!>   &                            GRID(ng)%vmask(Iend+1,j)
!>
              tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*            &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
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
!>          vbar(Istr-1,Jstr,kout)=0.5_r8*(vbar(Istr  ,Jstr  ,kout)+    &
!>   &                                     vbar(Istr-1,Jstr+1,kout))
!>
            tl_vbar(Istr-1,Jstr,kout)=0.5_r8*                           &
     &                                (tl_vbar(Istr  ,Jstr  ,kout)+     &
     &                                 tl_vbar(Istr-1,Jstr+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
!>          vbar(Iend+1,Jstr,kout)=0.5_r8*(vbar(Iend  ,Jstr  ,kout)+    &
!>   &                                     vbar(Iend+1,Jstr+1,kout))
!>
            tl_vbar(Iend+1,Jstr,kout)=0.5_r8*                           &
     &                                (tl_vbar(Iend  ,Jstr  ,kout)+     &
     &                                 tl_vbar(Iend+1,Jstr+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          vbar(Istr-1,Jend+1,kout)=0.5_r8*(vbar(Istr-1,Jend  ,kout)+  &
!>   &                                       vbar(Istr  ,Jend+1,kout))
!>
            tl_vbar(Istr-1,Jend+1,kout)=0.5_r8*                         &
     &                                  (tl_vbar(Istr-1,Jend  ,kout)+   &
     &                                   tl_vbar(Istr  ,Jend+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          vbar(Iend+1,Jend+1,kout)=0.5_r8*(vbar(Iend+1,Jend  ,kout)+  &
!>   &                                       vbar(Iend  ,Jend+1,kout))
!>
            tl_vbar(Iend+1,Jend+1,kout)=0.5_r8*                         &
     &                                  (tl_vbar(Iend+1,Jend  ,kout)+   &
     &                                   tl_vbar(Iend  ,Jend+1,kout))
          END IF
        END IF
      END IF

# if defined WET_DRY_NOT_YET
!
!-----------------------------------------------------------------------
!  Impose wetting and drying conditions.
!
!  HGA:  need TLM code here.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!>            cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,j))-1.0_r8)
!>            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,j,kout))*            &
!>   &                    GRID(ng)%vmask_wet(Istr-1,j)
!>            cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,j)*cff1+             &
!>   &            cff2*(1.0_r8-cff1)
!>            vbar(Istr,j,kout)=vbar(Istr,j,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!>            cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,j))-1.0_r8)
!>            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,j,kout))*            &
!>   &                    GRID(ng)%vmask_wet(Iend+1,j)
!>            cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,j)*cff1+             &
!>   &            cff2*(1.0_r8-cff1)
!>            vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*cff
            END IF
          END DO
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            cff1=ABS(ABS(GRID(ng)%vmask_wet(i,Jstr))-1.0_r8)
!>            cff2=0.5_r8+DSIGN(0.5_r8,vbar(i,Jstr,kout))*              &
!>   &                    GRID(ng)%vmask_wet(i,Jstr)
!>            cff=0.5_r8*GRID(ng)%vmask_wet(i,Jstr)*cff1+               &
!>   &            cff2*(1.0_r8-cff1)
!>            vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            cff1=ABS(ABS(GRID(ng)%vmask_wet(i,Jend+1))-1.0_r8)
!>            cff2=0.5_r8+DSIGN(0.5_r8,vbar(i,Jend+1,kout))*            &
!>   &                    GRID(ng)%vmask_wet(i,Jend+1)
!>            cff=0.5_r8*GRID(ng)%vmask_wet(i,Jend+1)*cff1+             &
!>   &            cff2*(1.0_r8-cff1)
!>            vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*cff
            END IF
          END DO
        END IF
      END IF

      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
!>          cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,Jstr))-1.0_r8)
!>          cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,Jstr,kout))*           &
!>   &                  GRID(ng)%vmask_wet(Istr-1,Jstr)
!>          cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,Jstr)*cff1+            &
!>   &          cff2*(1.0_r8-cff1)
!>          vbar(Istr-1,Jstr,kout)=vbar(Istr-1,Jstr,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
!>          cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,Jstr))-1.0_r8)
!>          cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,Jstr,kout))*           &
!>   &                  GRID(ng)%vmask_wet(Iend+1,Jstr)
!>          cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,Jstr)*cff1+            &
!>   &          cff2*(1.0_r8-cff1)
!>          vbar(Iend+1,Jstr,kout)=vbar(Iend+1,Jstr,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,Jend+1))-1.0_r8)
!>          cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,Jend+1,kout))*         &
!>   &                  GRID(ng)%vmask_wet(Istr-1,Jend+1)
!>          cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,Jend+1)*cff1+          &
!>   &          cff2*(1.0_r8-cff1)
!>          vbar(Istr-1,Jend+1,kout)=vbar(Istr-1,Jend+1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,Jend+1))-1.0_r8)
!>          cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,Jend+1,kout))*         &
!>   &                  GRID(ng)%vmask_wet(Iend+1,Jend+1)
!>          cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,Jend+1)*cff1+          &
!>   &          cff2*(1.0_r8-cff1)
!>          vbar(Iend+1,Jend+1,kout)=vbar(Iend+1,Jend+1,kout)*cff
          END IF
        END IF
      END IF
# endif

      RETURN
      END SUBROUTINE tl_v2dbc_tile
#endif
      END MODULE tl_v2dbc_mod
