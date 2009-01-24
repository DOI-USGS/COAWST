#include "cppdefs.h"
      MODULE tl_v2dbc_mod
#ifdef TANGENT
!
!svn $Id: tl_v2dbc_im.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
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

      PRIVATE
      PUBLIC  :: tl_v2dbc, tl_v2dbc_tile

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
      USE mod_forces
      USE mod_grid
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
      integer :: i, j, know

      real(r8) :: Ce, Cx
      real(r8) :: bry_pgr, bry_cor, bry_str
      real(r8) :: cff, cff1, cff2, dt2d, tau

      real(r8) :: tl_Ce, tl_Cx
      real(r8) :: tl_bry_pgr, tl_bry_cor, tl_bry_str, tl_bry_val
      real(r8) :: tl_cff, tl_cff1, tl_cff2

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

# ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN

#  if defined SOUTH_M2RADIATION

        IF (iic(ng).ne.0) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
          DO i=Istr,Iend+1
!>          grad(i,Jstr)=vbar(i  ,Jstr,know)-                           &
!>   &                   vbar(i-1,Jstr,know)
!>
            tl_grad(i,Jstr)=0.0_r8
          END DO
          DO i=Istr,Iend
#   ifdef SOUTH_M2NUDGING
            IF (BOUNDARY(ng)%vbar_south_Ce(i).eq.0.0_r8) THEN
              tau=M2obc_in(ng,isouth)
            ELSE
              tau=M2obc_out(ng,isouth)
            END IF
            tau=tau*dt2d
#   endif
#   ifdef RADIATION_2D
            Cx=BOUNDARY(ng)%vbar_south_Cx(i)
#   else
            Cx=0.0_r8
#   endif
            Ce=BOUNDARY(ng)%vbar_south_Ce(i)
            cff=BOUNDARY(ng)%vbar_south_C2(i)
!>          vbar(i,Jstr,kout)=(cff*vbar(i,Jstr  ,know)+                 &
!>   &                         Ce *vbar(i,Jstr+1,kout)-                 &
!>   &                         MAX(Cx,0.0_r8)*grad(i  ,Jstr)-           &
!>   &                         MIN(Cx,0.0_r8)*grad(i+1,Jstr))/          &
!>   &                        (cff+Ce)
!>
            tl_vbar(i,Jstr,kout)=(cff*tl_vbar(i,Jstr  ,know)+           &
     &                            Ce *tl_vbar(i,Jstr+1,kout)-           &
     &                            MAX(Cx,0.0_r8)*tl_grad(i  ,Jstr)-     &
     &                            MIN(Cx,0.0_r8)*tl_grad(i+1,Jstr))/    &
     &                           (cff+Ce)
#   ifdef SOUTH_M2NUDGING
!>          vbar(i,Jstr,kout)=vbar(i,Jstr,kout)+                        &
!>   &                        tau*(BOUNDARY(ng)%vbar_south(i)-          &
!>   &                             vbar(i,Jstr,know))
!>
            tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)-                  &
     &                           tau*tl_vbar(i,Jstr,know)
#   endif
#   ifdef MASKING
!>          vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                        &
!>   &                        GRID(ng)%vmask(i,Jstr)
!>
            tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                  &
     &                           GRID(ng)%vmask(i,Jstr)
#   endif
          END DO
        END IF

#  elif defined SOUTH_M2FLATHER
!
!  Southern edge, Flather boundary condition.
!
        DO i=Istr,Iend
#   if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
#    ifdef FSOBC_REDUCED
          bry_pgr=-g*(zeta(i,Jstr,know)-                                &
     &                BOUNDARY(ng)%zeta_south(i))*                      &
     &            0.5_r8*GRID(ng)%pn(i,Jstr)
          tl_bry_pgr=-g*tl_zeta(i,Jstr,know)*                           &
     &               0.5_r8*GRID(ng)%pn(i,Jstr)
#    else
          bry_pgr=-g*(zeta(i,Jstr  ,know)-                              &
     &                zeta(i,Jstr-1,know))*                             &
     &            0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                        &
     &                    GRID(ng)%pn(i,Jstr  ))
          tl_bry_pgr=-g*(tl_zeta(i,Jstr  ,know)-                        &
     &                   tl_zeta(i,Jstr-1,know))*                       &
     &               0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                     &
     &                       GRID(ng)%pn(i,Jstr  ))
#    endif
#    ifdef UV_COR
          bry_cor=-0.125_r8*(ubar(i  ,Jstr-1,know)+                     &
     &                       ubar(i+1,Jstr-1,know)+                     &
     &                       ubar(i  ,Jstr  ,know)+                     &
     &                       ubar(i+1,Jstr  ,know))*                    &
     &                      (GRID(ng)%f(i,Jstr-1)+                      &
     &                       GRID(ng)%f(i,Jstr  ))
          tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jstr-1,know)+               &
     &                          tl_ubar(i+1,Jstr-1,know)+               &
     &                          tl_ubar(i  ,Jstr  ,know)+               &
     &                          tl_ubar(i+1,Jstr  ,know))*              &
     &                         (GRID(ng)%f(i,Jstr-1)+                   &
     &                          GRID(ng)%f(i,Jstr  ))
#    else
          bry_cor=0.0_r8
          tl_bry_cor=0.0_r8
#    endif
          cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                    &
     &                         zeta(i,Jstr-1,know)+                     &
     &                         GRID(ng)%h(i,Jstr  )+                    &
     &                         zeta(i,Jstr  ,know)))
          tl_cff1=-cff1*cff1*(0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+          &
     &                                tl_zeta(i,Jstr-1,know)+           &
     &                                GRID(ng)%tl_h(i,Jstr  )+          &
     &                                tl_zeta(i,Jstr  ,know)))
          bry_str=cff1*(FORCES(ng)%svstr(i,Jstr)-                       &
     &                  FORCES(ng)%bvstr(i,Jstr))
          tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jstr)-                 &
     &                        FORCES(ng)%bvstr(i,Jstr))+                &
     &               cff1*(FORCES(ng)%tl_svstr(i,Jstr)-                 &
     &                     FORCES(ng)%tl_bvstr(i,Jstr))
          Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(i,Jstr-1)+                &
     &                             zeta(i,Jstr-1,know)+                 &
     &                             GRID(ng)%h(i,Jstr  )+                &
     &                             zeta(i,Jstr  ,know)))
          tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jstr-1)+           &
     &                               tl_zeta(i,Jstr-1,know)+            &
     &                               GRID(ng)%tl_h(i,Jstr  )+           &
     &                               tl_zeta(i,Jstr  ,know))
          cff2=GRID(ng)%on_v(i,Jstr)*Ce
          tl_cff2=GRID(ng)%on_v(i,Jstr)*tl_Ce
!>        bry_val=vbar(i,Jstr+1,know)+                                  &
!>   &            cff2*(bry_pgr+                                        &
!>   &                  bry_cor+                                        &
!>   &                  bry_str)
!>
          tl_bry_val=tl_vbar(i,Jstr+1,know)+                            &
     &               tl_cff2*(bry_pgr+                                  &
     &                        bry_cor+                                  &
     &                        bry_str)+                                 &
     &               cff2*(tl_bry_pgr+                                  &
     &                     tl_bry_cor+                                  &
     &                     tl_bry_str)
#   else
!>        bry_val=BOUNDARY(ng)%vbar_south(i)
!>
          tl_bry_val=0.0_r8
#   endif
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                     &
     &                        zeta(i,Jstr-1,know)+                      &
     &                        GRID(ng)%h(i,Jstr  )+                     &
     &                        zeta(i,Jstr  ,know)))
          tl_cff=-cff*cff*(0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+             &
     &                             tl_zeta(i,Jstr-1,know)+              &
     &                             GRID(ng)%tl_h(i,Jstr  )+             &
     &                             tl_zeta(i,Jstr  ,know)))
          Ce=SQRT(g*cff)
          tl_Ce=0.5_r8*g*tl_cff/Ce
!>        vbar(i,Jstr,kout)=bry_val-                                    &
!>   &                      Ce*(0.5_r8*(zeta(i,Jstr-1,know)+            &
!>   &                                  zeta(i,Jstr  ,know))-           &
!>   &                          BOUNDARY(ng)%zeta_south(i))
!>
          tl_vbar(i,Jstr,kout)=tl_bry_val-                              &
     &                         tl_Ce*(0.5_r8*(zeta(i,Jstr-1,know)+      &
     &                                        zeta(i,Jstr  ,know))-     &
     &                                BOUNDARY(ng)%zeta_south(i))-      &
     &                         Ce*(0.5_r8*(tl_zeta(i,Jstr-1,know)+      &
     &                                     tl_zeta(i,Jstr  ,know)))
#   ifdef MASKING
!>        vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                          &
!>   &                      GRID(ng)%vmask(i,Jstr)
!>
          tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                    &
     &                         GRID(ng)%vmask(i,Jstr)
#   endif
        END DO

#  elif defined SOUTH_M2CLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
!>        vbar(i,Jstr,kout)=BOUNDARY(ng)%vbar_south(i)
#   ifdef MASKING
!>        vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                          &
!>    &                     GRID(ng)%vmask(i,Jstr)
#   endif
!>
          tl_vbar(i,Jstr,kout)=0.0_r8
        END DO

#  elif defined SOUTH_M2GRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
!>        vbar(i,Jstr,kout)=vbar(i,Jstr+1,kout)
!>
          tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr+1,kout)
#   ifdef MASKING
!>        vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                          &
!>   &                      GRID(ng)%vmask(i,Jstr)
!>
          tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                    &
     &                         GRID(ng)%vmask(i,Jstr)
#   endif
        END DO

# elif defined SOUTH_M2REDUCED
!
!  Southern edge, reduced-physics boundary condition.
!
        DO i=Istr,Iend
#  ifdef FSOBC_REDUCED
!>        bry_pgr=-g*(zeta(i,Jstr,know)-                                &
!>   &                BOUNDARY(ng)%zeta_south(i))*                      &
!>   &            0.5_r8*GRID(ng)%pn(i,Jstr)
!>
          tl_bry_pgr=-g*tl_zeta(i,Jstr,know)*                           &
     &               0.5_r8*GRID(ng)%pn(i,Jstr)
#  else
!>        bry_pgr=-g*(zeta(i,Jstr,know)-                                &
!>   &                zeta(i,Jstr-1,know))*                             &
!>   &            0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                        &
!>   &                    GRID(ng)%pn(i,Jstr  ))
!>
          tl_bry_pgr=-g*(tl_zeta(i,Jstr  ,know)-                        &
     &                   tl_zeta(i,Jstr-1,know))*                       &
     &               0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                     &
     &                       GRID(ng)%pn(i,Jstr  ))
#  endif
#  ifdef UV_COR
!>        bry_cor=-0.125_r8*(ubar(i  ,Jstr-1,know)+                     &
!>   &                       ubar(i+1,Jstr-1,know)+                     &
!>   &                       ubar(i  ,Jstr  ,know)+                     &
!>   &                       ubar(i+1,Jstr  ,know))*                    &
!>   &                      (GRID(ng)%f(i,Jstr-1)+                      &
!>   &                       GRID(ng)%f(i,Jstr  ))
!>
          tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jstr-1,know)+               &
     &                          tl_ubar(i+1,Jstr-1,know)+               &
     &                          tl_ubar(i  ,Jstr  ,know)+               &
     &                          tl_ubar(i+1,Jstr  ,know))*              &
     &                         (GRID(ng)%f(i,Jstr-1)+                   &
     &                          GRID(ng)%f(i,Jstr  ))
#  else
!>        bry_cor=0.0_r8
!>
          tl_bry_cor=0.0_r8
#  endif
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                     &
     &                        zeta(i,Jstr-1,know)+                      &
     &                        GRID(ng)%h(i,Jstr  )+                     &
     &                        zeta(i,Jstr  ,know)))
          tl_cff=-cff*cff*0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+              &
     &                            tl_zeta(i,Jstr-1,know)+               &
     &                            GRID(ng)%tl_h(i,Jstr  )+              &
     &                            tl_zeta(i,Jstr  ,know))
!>        bry_str=cff*(FORCES(ng)%svstr(i,Jstr)-                        &
!>   &                 FORCES(ng)%bvstr(i,Jstr))
!>
          tl_bry_str=tl_cff*(FORCES(ng)%svstr(i,Jstr)-                  &
     &                       FORCES(ng)%bvstr(i,Jstr))+                 &
     &               cff*(FORCES(ng)%tl_svstr(i,Jstr)-                  &
     &                    FORCES(ng)%tl_bvstr(i,Jstr))

!>        vbar(i,Jstr,kout)=vbar(i,Jstr,know)+                          &
!>   &                      dt2d*(bry_pgr+                              &
!>   &                            bry_cor+                              &
!>   &                            bry_str)
!>
          tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,know)+                    &
     &                         dt2d*(tl_bry_pgr+                        &
     &                               tl_bry_cor+                        &
     &                               tl_bry_str)
#  ifdef MASKING
!>        vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                          &
!>   &                      GRID(ng)%vmask(i,Jstr)
!>
          tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                    &
     &                         GRID(ng)%vmask(i,Jstr)
#  endif
        END DO

#  else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
!>        vbar(i,Jstr,kout)=0.0_r8
!>
          tl_vbar(i,Jstr,kout)=0.0_r8
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN

#  if defined NORTH_M2RADIATION

        IF (iic(ng).ne.0) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
          DO i=Istr,Iend+1
!>          grad(i,Jend+1)=vbar(i  ,Jend+1,know)-                       &
!>   &                     vbar(i-1,Jend+1,know)
!>
            tl_grad(i,Jend+1)=0.0_r8
          END DO
          DO i=Istr,Iend
#   ifdef NORTH_M2NUDGING
            IF (BOUNDARY(ng)%vbar_north_Ce(i).eq.0.0_r8) THEN
              tau=M2obc_in(ng,inorth)
            ELSE
               tau=M2obc_out(ng,inorth)
            END IF
            tau=tau*dt2d
#   endif
#   ifdef RADIATION_2D
            Cx=BOUNDARY(ng)%vbar_north_Cx(i)
#   else
            Cx=0.0_r8
#   endif
            Ce=BOUNDARY(ng)%vbar_north_Ce(i)
            cff=BOUNDARY(ng)%vbar_north_C2(i)
!>          vbar(i,Jend+1,kout)=(cff*vbar(i,Jend+1,know)+               &
!>   &                           Ce *vbar(i,Jend  ,kout)-               &
!>   &                           MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-       &
!>   &                           MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/      &
!>   &                          (cff+Ce)
!>
            tl_vbar(i,Jend+1,kout)=(cff*tl_vbar(i,Jend+1,know)+         &
     &                              Ce *tl_vbar(i,Jend  ,kout)-         &
     &                              MAX(Cx,0.0_r8)*tl_grad(i  ,Jend+1)- &
     &                              MIN(Cx,0.0_r8)*tl_grad(i+1,Jend+1))/&
     &                             (cff+Ce)
#   ifdef NORTH_M2NUDGING
!>          vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)+                    &
!>   &                           tau*(BOUNDARY(ng)%vbar_north(i)-       &
!>   &                                vbar(i,Jend+1,know))
!>
            tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)-              &
     &                             tau*tl_vbar(i,Jend+1,know)
#   endif
#   ifdef MASKING
!>          vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                    &
!>   &                          GRID(ng)%vmask(i,Jend+1)
!>
            tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*              &
     &                             GRID(ng)%vmask(i,Jend+1)
#   endif
          END DO
        END IF

#  elif defined NORTH_M2FLATHER
!
!  Northern edge, Flather boundary condition.
!
        DO i=Istr,Iend
#   if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
#    ifdef FSOBC_REDUCED
          bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                       &
     &                zeta(i,Jend,know))*                               &
     &            0.5_r8*GRID(ng)%pn(i,Jend)
          tl_bry_pgr=g*tl_zeta(i,Jend,know)*                            &
     &               0.5_r8*GRID(ng)%pn(i,Jend)
#    else
          bry_pgr=-g*(zeta(i,Jend+1,know)-                              &
     &                zeta(i,Jend  ,know))*                             &
     &            0.5_r8*(GRID(ng)%pn(i,Jend  )+                        &
     &                    GRID(ng)%pn(i,Jend+1))
          tl_bry_pgr=-g*(tl_zeta(i,Jend+1,know)-                        &
     &                   tl_zeta(i,Jend  ,know))*                       &
     &               0.5_r8*(GRID(ng)%pn(i,Jend  )+                     &
     &                       GRID(ng)%pn(i,Jend+1))
#    endif
#    ifdef UV_COR
          bry_cor=-0.125_r8*(ubar(i  ,Jend  ,know)+                     &
     &                       ubar(i+1,Jend  ,know)+                     &
     &                       ubar(i  ,Jend+1,know)+                     &
     &                       ubar(i+1,Jend+1,know))*                    &
     &                      (GRID(ng)%f(i,Jend  )+                      &
     &                       GRID(ng)%f(i,Jend+1))
          tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jend  ,know)+               &
     &                          tl_ubar(i+1,Jend  ,know)+               &
     &                          tl_ubar(i  ,Jend+1,know)+               &
     &                          tl_ubar(i+1,Jend+1,know))*              &
     &                         (GRID(ng)%f(i,Jend  )+                   &
     &                          GRID(ng)%f(i,Jend+1))
#    else
          bry_cor=0.0_r8
          tl_bry_cor=0.0_r8
#    endif
          cff1=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                    &
     &                         zeta(i,Jend  ,know)+                     &
     &                         GRID(ng)%h(i,Jend+1)+                    &
     &                         zeta(i,Jend+1,know)))
          tl_cff1=-cff1*cff1*0.5_r8*(GRID(ng)%tl_h(i,Jend  )+           &
     &                               tl_zeta(i,Jend  ,know)+            &
     &                               GRID(ng)%tl_h(i,Jend+1)+           &
     &                               tl_zeta(i,Jend+1,know))
          bry_str=cff1*(FORCES(ng)%svstr(i,Jend+1)-                     &
     &                  FORCES(ng)%bvstr(i,Jend+1))
          tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jend+1)-               &
     &                        FORCES(ng)%bvstr(i,Jend+1))+              &
     &               cff1*(FORCES(ng)%tl_svstr(i,Jend+1)-               &
     &                     FORCES(ng)%tl_bvstr(i,Jend+1))
          Ce=1.0_r8/SQRT(g*0.5_r8*(GRID(ng)%h(i,Jend+1)+                &
     &                             zeta(i,Jend+1,know)+                 &
     &                             GRID(ng)%h(i,Jend  )+                &
     &                             zeta(i,Jend  ,know)))
          tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jend+1)+           &
     &                               tl_zeta(i,Jend+1,know)+            &
     &                               GRID(ng)%tl_h(i,Jend  )+           &
     &                               tl_zeta(i,Jend  ,know))
          cff2=GRID(ng)%on_v(i,Jend+1)*Ce
          tl_cff2=GRID(ng)%on_v(i,Jend+1)*tl_Ce
!>        bry_val=vbar(i,Jend,know)+                                    &
!>   &            cff2*(bry_pgr+                                        &
!>   &                  bry_cor+                                        &
!>   &                  bry_str)
!>
          tl_bry_val=tl_vbar(i,Jend,know)+                              &
     &               tl_cff2*(bry_pgr+                                  &
     &                        bry_cor+                                  &
     &                        bry_str)+                                 &
     &               cff2*(tl_bry_pgr+                                  &
     &                     tl_bry_cor+                                  &
     &                     tl_bry_str)
#   else
!>        bry_val=BOUNDARY(ng)%vbar_north(i)
!>
          tl_bry_val=0.0_r8
#   endif
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                     &
     &                        zeta(i,Jend  ,know)+                      &
     &                        GRID(ng)%h(i,Jend+1)+                     &
     &                        zeta(i,Jend+1,know)))
          tl_cff=-cff*cff*(0.5_r8*(GRID(ng)%tl_h(i,Jend  )+             &
     &                             tl_zeta(i,Jend  ,know)+              &
     &                             GRID(ng)%tl_h(i,Jend+1)+             &
     &                             tl_zeta(i,Jend+1,know)))
          Ce=SQRT(g*cff)
          tl_Ce=0.5_r8*g*tl_cff/Ce
!>        vbar(i,Jend+1,kout)=bry_val+                                  &
!>   &                        Ce*(0.5_r8*(zeta(i,Jend  ,know)+          &
!>   &                                    zeta(i,Jend+1,know))-         &
!>   &                            BOUNDARY(ng)%zeta_north(i))
!>
          tl_vbar(i,Jend+1,kout)=tl_bry_val+                            &
     &                           tl_Ce*(0.5_r8*(zeta(i,Jend  ,know)+    &
     &                                          zeta(i,Jend+1,know))-   &
     &                                  BOUNDARY(ng)%zeta_north(i))+    &
     &                           Ce*(0.5_r8*(tl_zeta(i,Jend  ,know)+    &
     &                                       tl_zeta(i,Jend+1,know)))
#   ifdef MASKING
!>        vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                      &
!>   &                        GRID(ng)%vmask(i,Jend+1)
!>
          tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*                &
     &                           GRID(ng)%vmask(i,Jend+1)
#   endif
        END DO

#  elif defined NORTH_M2CLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
!>        vbar(i,Jend+1,kout)=BOUNDARY(ng)%vbar_north(i)
#   ifdef MASKING
!>        vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                      &
!>   &                        GRID(ng)%vmask(i,Jend+1)
#   endif
!>
          tl_vbar(i,Jend+1,kout)=0.0_r8
        END DO

#  elif defined NORTH_M2GRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
!>        vbar(i,Jend+1,kout)=vbar(i,Jend,kout)
!>
          tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend,kout)
#   ifdef MASKING
!>        vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                      &
!>   &                        GRID(ng)%vmask(i,Jend+1)
!>
          tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*                &
     &                           GRID(ng)%vmask(i,Jend+1)
#   endif
        END DO

# elif defined NORTH_M2REDUCED
!
!  Northern edge, reduced-physics boundary condition.
!
        DO i=Istr,Iend
#  ifdef FSOBC_REDUCED
!>        bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                       &
!>   &                zeta(i,Jend,know))*                               &
!>   &            0.5_r8*GRID(ng)%pn(i,Jend)
!>
          tl_bry_pgr=g*zeta(i,Jend,know))*                              &
     &               0.5_r8*GRID(ng)%pn(i,Jend)
#  else
!>        bry_pgr=-g*(zeta(i,Jend+1,know)-                              &
!>   &                zeta(i,Jend  ,know))*                             &
!>   &            0.5_r8*(GRID(ng)%pn(i,Jend  )+                        &
!>   &                    GRID(ng)%pn(i,Jend+1))
!>
          tl_bry_pgr=-g*(tl_zeta(i,Jend+1,know)-                        &
     &                   tl_zeta(i,Jend  ,know))*                       &
     &               0.5_r8*(GRID(ng)%pn(i,Jend  )+                     &
     &                       GRID(ng)%pn(i,Jend+1))
#  endif
#  ifdef UV_COR
!>        bry_cor=-0.125_r8*(ubar(i  ,Jend  ,know)+                     &
!>   &                       ubar(i+1,Jend  ,know)+                     &
!>   &                       ubar(i  ,Jend+1,know)+                     &
!>   &                       ubar(i+1,Jend+1,know))*                    &
!>   &                      (GRID(ng)%f(i,Jend  )+                      &
!>   &                       GRID(ng)%f(i,Jend+1))
!>
          tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jend  ,know)+               &
     &                          tl_ubar(i+1,Jend  ,know)+               &
     &                          tl_ubar(i  ,Jend+1,know)+               &
     &                          tl_ubar(i+1,Jend+1,know))*              &
     &                         (GRID(ng)%f(i,Jend  )+                   &
     &                          GRID(ng)%f(i,Jend+1))
#  else
!>        bry_cor=0.0_r8
!>
          tl_bry_cor=0.0_r8
#  endif
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                     &
     &                        zeta(i,Jend  ,know)+                      &
     &                        GRID(ng)%h(i,Jend+1)+                     &
     &                        zeta(i,Jend+1,know)))
          tl_cff=-cff*cff*0.5_r8*(GRID(ng)%tl_h(i,Jend  )+              &
     &                            tl_zeta(i,Jend  ,know)+               &
     &                            GRID(ng)%tl_h(i,Jend+1)+              &
     &                            tl_zeta(i,Jend+1,know))
!>        bry_str=cff*(FORCES(ng)%svstr(i,Jend+1)-                      &
!>   &                 FORCES(ng)%bvstr(i,Jend+1))
!>
          tl_bry_str=tl_cff*(FORCES(ng)%svstr(i,Jend+1)-                &
     &                       FORCES(ng)%bvstr(i,Jend+1))+               &
     &               cff*(FORCES(ng)%tl_svstr(i,Jend+1)-                &
     &                    FORCES(ng)%tl_bvstr(i,Jend+1))
!>        vbar(i,Jend+1,kout)=vbar(i,Jend+1,know)+                      &
!>   &                      dt2d*(bry_pgr+                              &
!>   &                            bry_cor+                              &
!>   &                            bry_str)
!>
          tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,know)+                &
     &                           dt2d*(tl_bry_pgr+                      &
     &                                 tl_bry_cor+                      &
     &                                 tl_bry_str)
#  ifdef MASKING
!>        vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                      &
!>   &                        GRID(ng)%vmask(i,Jend+1)
!>
          tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*                &
     &                           GRID(ng)%vmask(i,Jend+1)
#  endif
        END DO

#  else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
!>        vbar(i,Jend+1,kout)=0.0_r8
!>
          tl_vbar(i,Jend+1,kout)=0.0_r8
        END DO
#  endif
      END IF
# endif

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

#  if defined WEST_M2RADIATION

        IF (iic(ng).ne.0) THEN
!
!  Western edge, implicit upstream radiation condition.
!
          DO j=JstrV-1,Jend
!>          grad(Istr-1,j)=vbar(Istr-1,j+1,know)-                       &
!>   &                     vbar(Istr-1,j  ,know)
!>
            tl_grad(Istr-1,j)=0.0_r8
          END DO
          DO j=JstrV,Jend
#   ifdef WEST_M2NUDGING
            IF (BOUNDARY(ng)%vbar_west_Cx(j).eq.0.0_r8) THEN
     &        tau=M2obc_in(ng,iwest)
            ELSE
              tau=M2obc_out(ng,iwest)
            END IF
            tau=tau*dt2d
#   endif
            Cx=BOUNDARY(ng)%vbar_west_Cx(j)
#   ifdef RADIATION_2D
            Ce=BOUNDARY(ng)%vbar_west_Ce(j)
#   else
            Ce=0.0_r8
#   endif
            cff=BOUNDARY(ng)%vbar_west_C2(j)
!>          vbar(Istr-1,j,kout)=(cff*vbar(Istr-1,j,know)+               &
!>   &                           Cx *vbar(Istr  ,j,kout)-               &
!>   &                           MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-       &
!>   &                           MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/      &
!>   &                          (cff+Cx)
!>
            tl_vbar(Istr-1,j,kout)=(cff*tl_vbar(Istr-1,j,know)+         &
     &                              Cx *tl_vbar(1,j,kout)-              &
     &                              MAX(Ce,0.0_r8)*tl_grad(Istr-1,j-1)- &
     &                              MIN(Ce,0.0_r8)*tl_grad(Istr-1,j  ))/&
     &                             (cff+Cx)
#   ifdef WEST_M2NUDGING
!>          vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)+                    &
!>   &                          tau*(BOUNDARY(ng)%vbar_west(j)-         &
!>   &                               vbar(Istr-1,j,know))
!>
            tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)-              &
     &                             tau*tl_vbar(1,j,know)
#   endif
#   ifdef MASKING
!>           vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                   &
!>   &                           GRID(ng)%vmask(Istr-1,j)
!>
             tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*             &
     &                              GRID(ng)%vmask(Istr-1,j)
#   endif
          END DO
        END IF

#  elif defined WEST_M2FLATHER || defined WEST_M2REDUCED
!
!  Western edge, Chapman boundary condition.
!
        DO j=JstrV,Jend
          cff=dt2d*0.5_r8*(GRID(ng)%pm(Istr,j-1)+                       &
     &                     GRID(ng)%pm(Istr,j  ))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(Istr,j-1)+                     &
     &                        zeta(Istr,j-1,know)+                      &
     &                        GRID(ng)%h(Istr,j  )+                     &
     &                        zeta(Istr,j  ,know)))
          tl_cff1=0.25_r8*g*(GRID(ng)%tl_h(Istr,j-1)+                   &
     &                       tl_zeta(Istr,j-1,know)+                    &
     &                       GRID(ng)%tl_h(Istr,j  )+                   &
     &                       tl_zeta(Istr,j  ,know))/cff1
          Cx=cff*cff1
          tl_Cx=cff*tl_cff1
          cff2=1.0_r8/(1.0_r8+Cx)
          tl_cff2=-cff2*cff2*tl_Cx
!>        vbar(Istr-1,j,kout)=cff2*(vbar(Istr-1,j,know)+                &
!>   &                              Cx*vbar(Istr,j,kout))
!>
          tl_vbar(Istr-1,j,kout)=tl_cff2*(vbar(Istr-1,j,know)+          &
     &                                    Cx*vbar(Istr,j,kout))+        &
     &                           cff2*(tl_vbar(Istr-1,j,know)+          &
     &                                 tl_Cx*vbar(Istr,j,kout)+         &
     &                                 Cx*tl_vbar(Istr,j,kout))
#   ifdef MASKING
!>        vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                      &
!>   &                        GRID(ng)%vmask(Istr-1,j)
!>
          tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*                &
     &                           GRID(ng)%vmask(Istr-1,j)
#   endif
        END DO

#  elif defined WEST_M2CLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=JstrV,Jend
!>        vbar(Istr-1,j,kout)=BOUNDARY(ng)%vbar_west(j)
#   ifdef MASKING
!>        vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                      &
!>   &                        GRID(ng)%vmask(Istr-1,j)
#   endif
!>
          tl_vbar(Istr-1,j,kout)=0.0_r8
        END DO

#  elif defined WEST_M2GRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=JstrV,Jend
!>        vbar(Istr-1,j,kout)=vbar(Istr,j,kout)
!>
          tl_vbar(Istr-1,j,kout)=tl_vbar(Istr,j,kout)
#   ifdef MASKING
!>        vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                      &
!>   &                        GRID(ng)%vmask(Istr-1,j)
!>
          tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*                &
     &                           GRID(ng)%vmask(Istr-1,j)
#   endif
        END DO

#  else
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO j=J_RANGE
!>        vbar(Istr-1,j,kout)=gamma2(ng)*vbar(Istr,j,kout)
!>
          tl_vbar(Istr-1,j,kout)=gamma2(ng)*tl_vbar(Istr,j,kout)
#   ifdef MASKING
!>        vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                      &
!>   &                        GRID(ng)%vmask(Istr-1,j)
!>
          tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*                &
     &                           GRID(ng)%vmask(Istr-1,j)
#   endif
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

#  if defined EAST_M2RADIATION

        IF (iic(ng).ne.0) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
          DO j=JstrV-1,Jend
!>          grad(Iend+1,j)=vbar(Iend+1,j+1,know)-                       &
!>   &                     vbar(Iend+1,j  ,know)
!>
            tl_grad(Iend+1,j)=0.0_r8
          END DO
          DO j=JstrV,Jend
#   ifdef EAST_M2NUDGING
            IF (BOUNDARY(ng)%vbar_east_Cx(j).eq.0.0_r8) THEN
              tau=M2obc_in(ng,ieast)
            ELSE
              tau=M2obc_out(ng,ieast)
            END IF
            tau=tau*dt2d
#   endif
            Cx=BOUNDARY(ng)%vbar_east_Cx(j)
#   ifdef RADIATION_2D
            Ce=BOUNDARY(ng)%vbar_east_Ce(j)
#   else
            Ce=0.0_r8
#   endif
            cff=BOUNDARY(ng)%vbar_east_C2(j)
!>          vbar(Iend+1,j,kout)=(cff*vbar(Iend+1,j,know)+               &
!>   &                           Cx *vbar(Iend  ,j,kout)-               &
!>   &                           MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-       &
!>   &                           MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/      &
!>   &                          (cff+Cx)
!>
            tl_vbar(Iend+1,j,kout)=(cff*tl_vbar(Iend+1,j,know)+         &
     &                              Cx *tl_vbar(Iend  ,j,kout)-         &
     &                              MAX(Ce,0.0_r8)*tl_grad(Iend+1,j-1)- &
     &                              MIN(Ce,0.0_r8)*tl_grad(Iend+1,j  ))/&
     &                             (cff+Cx)
#   ifdef EAST_M2NUDGING
!>          vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)+                    &
!>   &                          tau*(BOUNDARY(ng)%vbar_east(j)-         &
!>   &                               vbar(Iend+1,j,know))
!>
            tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)-              &
     &                             tau*tl_vbar(Iend+1,j,know)
#   endif
#   ifdef MASKING
!>          vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                    &
!>   &                          GRID(ng)%vmask(Iend+1,j)
!>
            tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*              &
     &                             GRID(ng)%vmask(Iend+1,j)
#   endif
          END DO
        END IF

#  elif defined EAST_M2FLATHER || defined EAST_M2REDUCED
!
!  Eastern edge, Chapman boundary condition.
!
        DO j=JstrV,Jend
          cff=dt2d*0.5_r8*(GRID(ng)%pm(Iend,j-1)+                       &
     &                     GRID(ng)%pm(Iend,j  ))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(Iend,j-1)+                     &
     &                        zeta(Iend,j-1,know)+                      &
     &                        GRID(ng)%h(Iend,j  )+                     &
     &                        zeta(Iend,j  ,know)))
          tl_cff1=0.25_r8*g*(GRID(ng)%tl_h(Iend,j-1)+                   &
     &                       tl_zeta(Iend,j-1,know)+                    &
     &                       GRID(ng)%tl_h(Iend,j  )+                   &
     &                       tl_zeta(Iend,j  ,know))/cff1
          Cx=cff*cff1
          tl_Cx=cff*tl_cff1
          cff2=1.0_r8/(1.0_r8+Cx)
          tl_cff2=-cff2*cff2*tl_Cx
!>        vbar(Iend+1,j,kout)=cff2*(vbar(Iend+1,j,know)+                &
!>   &                              Cx*vbar(Iend,j,kout))
!>
          tl_vbar(Iend+1,j,kout)=tl_cff2*(vbar(Iend+1,j,know)+          &
     &                                    Cx*vbar(Iend,j,kout))+        &
     &                           cff2*(tl_vbar(Iend+1,j,know)+          &
     &                                 tl_Cx*vbar(Iend,j,kout)+         &
     &                                 Cx*tl_vbar(Iend,j,kout))
#   ifdef MASKING
!>        vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                      &
!>   &                        GRID(ng)%vmask(Iend+1,j)
!>
          tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*                &
     &                           GRID(ng)%vmask(Iend+1,j)
#   endif
        END DO

#  elif defined EAST_M2CLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=JstrV,Jend
!>        vbar(Iend+1,j,kout)=BOUNDARY(ng)%vbar_east(j)
#   ifdef MASKING
!>        vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                      &
!>   &                        GRID(ng)%vmask(Iend+1,j)
#   endif
!>
          tl_vbar(Iend+1,j,kout)=0.0_r8
        END DO

#  elif defined EAST_M2GRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=JstrV,Jend
!>        vbar(Iend+1,j,kout)=vbar(Iend,j,kout)
!>
          tl_vbar(Iend+1,j,kout)=tl_vbar(Iend,j,kout)
#   ifdef MASKING
!>        vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                      &
!>   &                        GRID(ng)%vmask(Iend+1,j)
!>
          tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*                &
     &                           GRID(ng)%vmask(Iend+1,j)
#   endif
        END DO

#  else
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
#   ifdef NS_PERIODIC
#    define J_RANGE JstrV,Jend
#   else
#    define J_RANGE Jstr,JendR
#   endif
        DO j=J_RANGE
!>        vbar(Iend+1,j,kout)=gamma2(ng)*vbar(Iend,j,kout)
!>
          tl_vbar(Iend+1,j,kout)=gamma2(ng)*tl_vbar(Iend,j,kout)
#   ifdef MASKING
!>        vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                      &
!>   &                        GRID(ng)%vmask(Iend+1,j)
!>
          tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*                &
     &                           GRID(ng)%vmask(Iend+1,j)
#   endif
        END DO
#   undef J_RANGE
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
!>      vbar(Istr-1,Jstr,kout)=0.5_r8*(vbar(Istr  ,Jstr  ,kout)+        &
!>   &                                 vbar(Istr-1,Jstr+1,kout))
!>
        tl_vbar(Istr-1,Jstr,kout)=0.5_r8*                               &
     &                            (tl_vbar(Istr  ,Jstr  ,kout)+         &
     &                             tl_vbar(Istr-1,Jstr+1,kout))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>      vbar(Iend+1,Jstr,kout)=0.5_r8*(vbar(Iend  ,Jstr  ,kout)+        &
!>   &                                 vbar(Iend+1,Jstr+1,kout))
!>
        tl_vbar(Iend+1,Jstr,kout)=0.5_r8*                               &
     &                            (tl_vbar(Iend  ,Jstr  ,kout)+         &
     &                             tl_vbar(Iend+1,Jstr+1,kout))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>      vbar(Istr-1,Jend+1,kout)=0.5_r8*(vbar(Istr-1,Jend  ,kout)+      &
!>   &                                   vbar(Istr  ,Jend+1,kout))
!>
        tl_vbar(Istr-1,Jend+1,kout)=0.5_r8*                             &
     &                              (tl_vbar(Istr-1,Jend  ,kout)+       &
     &                               tl_vbar(Istr  ,Jend+1,kout))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>      vbar(Iend+1,Jend+1,kout)=0.5_r8*(vbar(Iend+1,Jend  ,kout)+      &
!>   &                                   vbar(Iend  ,Jend+1,kout))
!>
        tl_vbar(Iend+1,Jend+1,kout)=0.5_r8*                             &
     &                              (tl_vbar(Iend+1,Jend  ,kout)+       &
     &                               tl_vbar(Iend  ,Jend+1,kout))
      END IF
# endif
      RETURN
      END SUBROUTINE tl_v2dbc_tile
#endif
      END MODULE tl_v2dbc_mod
