#include "cppdefs.h"
      MODULE ad_u2dbc_mod
#ifdef ADJOINT
!
!svn $Id: ad_u2dbc_im.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets adjoint lateral boundary conditions for        !
!  vertically integrated  U-velocity.  It updates the specified        !
!  "kout" index.                                                       !
!                                                                      !
!  BASIC STATE variables needed: zeta                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_u2dbc, ad_u2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_u2dbc (ng, tile, kout)
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
      CALL ad_u2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs(ng), kstp(ng), kout,                     &
     &                    OCEAN(ng) % ubar,                             &
     &                    OCEAN(ng) % vbar,                             &
     &                    OCEAN(ng) % zeta,                             &
     &                    OCEAN(ng) % ad_ubar,                          &
     &                    OCEAN(ng) % ad_ubar,                          &
     &                    OCEAN(ng) % ad_zeta)

      RETURN
      END SUBROUTINE ad_u2dbc

!
!***********************************************************************
      SUBROUTINE ad_u2dbc_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          krhs, kstp, kout,                       &
     &                          ubar, vbar, zeta,                       &
     &                          ad_ubar, ad_vbar, ad_zeta)
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

      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j, know

      real(r8) :: Ce, Cx, Zx
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8) :: cff, cff1, cff2, cff3, dt2d
      real(r8) :: obc_in, obc_out, tau

      real(r8) :: ad_Ce, ad_Cx
      real(r8) :: ad_bry_pgr, ad_bry_cor, ad_bry_str, ad_bry_val, ad_Zx
      real(r8) :: ad_cff, ad_cff1, ad_cff2, ad_cff3
      real(r8) :: adfac

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_grad

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_Ce=0.0_r8
      ad_Cx=0.0_r8
      ad_Zx=0.0_r8
      ad_cff=0.0_r8
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_cff3=0.0_r8
      ad_bry_pgr=0.0_r8
      ad_bry_cor=0.0_r8
      ad_bry_str=0.0_r8
      ad_bry_val=0.0_r8

      ad_grad(LBi:UBi,LBj:UBj)=0.0_r8
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

# if defined WET_DRY_NOT_YET
!
!-----------------------------------------------------------------------
!  Impose wetting and drying conditions.
!
!  HGA:  need ADM code here for the NLM code below.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          cff1=ABS(ABS(GRID(ng)%umask_wet(Iend+1,Jend+1))-1.0_r8)
!>          cff2=0.5_r8+DSIGN(0.5_r8,ubar(Iend+1,Jend+1,kout))*         &
!>   &                  GRID(ng)%umask_wet(Iend+1,Jend+1)
!>          cff=0.5_r8*GRID(ng)%umask_wet(Iend+1,Jend+1)+cff1+          &
!>   &          cff2*(1.0_r8-cff1)
!>          ubar(Iend+1,Jend+1,kout)=ubar(Iend+1,Jend+1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          cff1=ABS(ABS(GRID(ng)%umask_wet(Istr,Jend+1))-1.0_r8)
!>          cff2=0.5_r8+DSIGN(0.5_r8,ubar(Istr,Jend+1,kout))*           &
!>   &                  GRID(ng)%umask_wet(Istr,Jend+1)
!>          cff=0.5_r8*GRID(ng)%umask_wet(Istr,Jend+1)*cff1+            &
!>   &          cff2*(1.0_r8-cff1)
!>          ubar(Istr,Jend+1,kout)=ubar(Istr,Jend+1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
!>          cff1=ABS(ABS(GRID(ng)%umask_wet(Iend+1,Jstr-1))-1.0_r8)
!>          cff2=0.5_r8+DSIGN(0.5_r8,ubar(Iend+1,Jstr-1,kout))*         &
!>   &                  GRID(ng)%umask_wet(Iend+1,Jstr-1)
!>          cff=0.5_r8*GRID(ng)%umask_wet(Iend+1,Jstr-1)*cff1+          &
!>   &          cff2*(1.0_r8-cff1)
!>          ubar(Iend+1,Jstr-1,kout)=ubar(Iend+1,Jstr-1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
!>          cff1=ABS(ABS(GRID(ng)%umask_wet(Istr,Jstr-1))-1.0_r8)
!>          cff2=0.5_r8+DSIGN(0.5_r8,ubar(Istr,Jstr-1,kout))*           &
!>   &                  GRID(ng)%umask_wet(Istr,Jstr-1)
!>          cff=0.5_r8*GRID(ng)%umask_wet(Istr,Jstr-1)*cff1+            &
!>   &          cff2*(1.0_r8-cff1)
!>          ubar(Istr,Jstr-1,kout)=ubar(Istr,Jstr-1,kout)*cff
          END IF
        END IF
      END IF

      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            cff1=ABS(ABS(GRID(ng)%umask_wet(i,Jend+1))-1.0_r8)
!>            cff2=0.5_r8+DSIGN(0.5_r8,ubar(i,Jend+1,kout))*            &
!>   &                    GRID(ng)%umask_wet(i,Jend+1)
!>            cff=0.5_r8*GRID(ng)%umask_wet(i,Jend+1)*cff1+             &
!>   &            cff2*(1.0_r8-cff1)
!>            ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            cff1=ABS(ABS(GRID(ng)%umask_wet(i,Jstr-1))-1.0_r8)
!>            cff2=0.5_r8+DSIGN(0.5_r8,ubar(i,Jstr-1,kout))*            &
!>   &                    GRID(ng)%umask_wet(i,Jstr-1)
!>            cff=0.5_r8*GRID(ng)%umask_wet(i,Jstr-1)*cff1+             &
!>   &            cff2*(1.0_r8-cff1)
!>            ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*cff
            END IF
          END DO
        END IF
      END IF

      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!>            cff1=ABS(ABS(GRID(ng)%umask_wet(Iend+1,j))-1.0_r8)
!>            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Iend+1,j,kout))*            &
!>   &                    GRID(ng)%umask_wet(Iend+1,j)
!>            cff=0.5_r8*GRID(ng)%umask_wet(Iend+1,j)*cff1+             &
!>   &            cff2*(1.0_r8-cff1)
!>            ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!>            cff1=ABS(ABS(GRID(ng)%umask_wet(Istr,j))-1.0_r8)
!>            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Istr,j,kout))*              &
!>   &                    GRID(ng)%umask_wet(Istr,j)
!>            cff=0.5_r8*GRID(ng)%umask_wet(Istr,j)*cff1+               &
!>   &            cff2*(1.0_r8-cff1)
!>            ubar(Istr,j,kout)=ubar(Istr,j,kout)*cff
            END IF
          END DO
        END IF
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
!>          tl_ubar(Iend+1,Jend+1,kout)=0.5_r8*                         &
!>   &                                  (tl_ubar(Iend+1,Jend  ,kout)+   &
!>   &                                   tl_ubar(Iend  ,Jend+1,kout))
!>
            adfac=0.5_r8*ad_ubar(Iend+1,Jend+1,kout)
            ad_ubar(Iend+1,Jend  ,kout)=ad_ubar(Iend+1,Jend  ,kout)+    &
     &                                  adfac
            ad_ubar(Iend  ,Jend+1,kout)=ad_ubar(Iend  ,Jend+1,kout)+    &
     &                                  adfac
            ad_ubar(Iend+1,Jend+1,kout)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          tl_ubar(Istr,Jend+1,kout)=0.5_r8*                           &
!>   &                                (tl_ubar(Istr  ,Jend  ,kout)+     &
!>   &                                 tl_ubar(Istr+1,Jend+1,kout))
!>
            adfac=0.5_r8*ad_ubar(Istr,Jend+1,kout)
            ad_ubar(Istr  ,Jend  ,kout)=ad_ubar(Istr  ,Jend  ,kout)+    &
     &                                  adfac
            ad_ubar(Istr+1,Jend+1,kout)=ad_ubar(Istr+1,Jend+1,kout)+    &
     &                                  adfac
            ad_ubar(Istr  ,Jend+1,kout)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
!>          tl_ubar(Iend+1,Jstr-1,kout)=0.5_r8*                         &
!>   &                                  (tl_ubar(Iend  ,Jstr-1,kout)+   &
!>   &                                   tl_ubar(Iend+1,Jstr  ,kout))
!>
            adfac=0.5_r8*ad_ubar(Iend+1,Jstr-1,kout)
            ad_ubar(Iend  ,Jstr-1,kout)=ad_ubar(Iend  ,Jstr-1,kout)+    &
     &                                  adfac
            ad_ubar(Iend+1,Jstr  ,kout)=ad_ubar(Iend+1,Jstr  ,kout)+    &
     &                                  adfac
            ad_ubar(Iend+1,Jstr-1,kout)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
!>          ubar(Istr,Jstr-1,kout)=0.5_r8*                              &
!>   &                             (ubar(Istr+1,Jstr-1,kout)+           &
!>   &                              ubar(Istr  ,Jstr  ,kout))
!>
            adfac=0.5_r8*ad_ubar(Istr,Jstr-1,kout)
            ad_ubar(Istr+1,Jstr-1,kout)=ad_ubar(Istr+1,Jstr-1,kout)+    &
     &                                  adfac
            ad_ubar(Istr  ,Jstr  ,kout)=ad_ubar(Istr  ,Jstr  ,kout)+    &
     &                                  adfac
            ad_ubar(Istr  ,Jstr-1,kout)=0.0_r8
          END IF
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
        IF (ad_LBC(inorth,isUbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (ad_LBC(inorth,isUbar,ng)%nudging) THEN
                  IF (LnudgeM2CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M2nudgcof(i-1,Jend+1)+           &
     &                       CLIMA(ng)%M2nudgcof(i  ,Jend+1))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M2obc_out(ng,inorth)
                    obc_in =M2obc_in (ng,inorth)
                  END IF
                  IF (BOUNDARY(ng)%ubar_north_Ce(i).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt2d
                END IF
#  ifdef RADIATION_2D
                Cx=BOUNDARY(ng)%ubar_north_Cx(i)
#  else
                Cx=0.0_r8
#  endif
                Ce=BOUNDARY(ng)%ubar_north_Ce(i)
                cff=BOUNDARY(ng)%ubar_north_C2(i)
# endif
# ifdef MASKING
!>              tl_ubar(i,Jend+1,kout)=tl_ubar(i,Jend+1,kout)*          &
!>   &                                 GRID(ng)%umask(i,Jend+1)
!>
                ad_ubar(i,Jend+1,kout)=ad_ubar(i,Jend+1,kout)*          &
     &                                 GRID(ng)%umask(i,Jend+1)
# endif
                IF (ad_LBC(inorth,isUbar,ng)%nudging) THEN
!>                tl_ubar(i,Jend+1,kout)=tl_ubar(i,Jend+1,kout)-        &
!>   &                                   tau*tl_ubar(i,Jend+1,know)
!>
                  ad_ubar(i,Jend+1 ,know)=ad_ubar(i,Jend+1 ,know)-      &
     &                                    tau*ad_ubar(i,Jend+1,kout)
                END IF
!>              tl_ubar(i,Jend+1,kout)=(cff*tl_ubar(i,Jend+1,know)+     &
!>   &                                  Ce *tl_ubar(i,Jend  ,kout)-     &
!>   &                                  MAX(Cx,0.0_r8)*                 &
!>   &                                     tl_grad(i-1,Jend+1)-         &
!>   &                                  MIN(Cx,0.0_r8)*                 &
!>   &                                     tl_grad(i  ,Jend+1))/        &
!>   &                                 (cff+Ce)
!>
                adfac=ad_ubar(i,Jend+1,kout)/(cff+Ce)
                ad_grad(i-1,Jend+1)=ad_grad(i-1,Jend+1)-                &
     &                              MAX(Cx,0.0_r8)*adfac
                ad_grad(i  ,Jend+1)=ad_grad(i  ,Jend+1)-                &
     &                              MIN(Cx,0.0_r8)*adfac
                ad_ubar(i,Jend  ,kout)=ad_ubar(i,Jend  ,kout)+Ce *adfac
                ad_ubar(i,Jend+1,know)=ad_ubar(i,Jend+1,know)+cff*adfac
                ad_ubar(i,Jend+1,kout)=0.0_r8
              END IF
            END DO
          END IF
!
!  Northern edge, Chapman boundary condition.
!
        ELSE IF (ad_LBC(inorth,isUbar,ng)%Flather.or.                   &
     &           ad_LBC(inorth,isUbar,ng)%Shchepetkin.or.               &
     &           ad_LBC(inorth,isUbar,ng)%reduced) THEN
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
# ifdef MASKING
!>            tl_ubar(i,Jend+1,kout)=tl_ubar(i,Jend+1,kout)*            &
!>   &                               GRID(ng)%umask(i,Jend+1)
!>
              ad_ubar(i,Jend+1,kout)=ad_ubar(i,Jend+1,kout)*            &
     &                               GRID(ng)%umask(i,Jend+1)
# endif
!>            tl_ubar(i,Jend+1,kout)=tl_cff2*(ubar(i,Jend+1,know)+      &
!>   &                                        Ce*ubar(i,Jend,kout))+    &
!>   &                               cff2*(tl_ubar(i,Jend+1,know)+      &
!>   &                                     tl_Ce*ubar(i,Jend,kout)+     &
!>   &                                     Ce*tl_ubar(i,Jend,kout))
!>
              adfac=cff2*ad_ubar(i,Jend+1,kout)
              ad_ubar(i,Jend  ,kout)=ad_ubar(i,Jend  ,kout)+Ce*adfac
              ad_ubar(i,Jend+1,know)=ad_ubar(i,Jend+1,know)+adfac
              ad_Ce=ad_Ce+ubar(i,Jend,kout)*adfac
              ad_cff2=ad_cff2+                                          &
     &                (ubar(i,Jend+1,know)+                             &
     &                 Ce*ubar(i,Jend,kout))*ad_ubar(i,Jend+1,kout)
              ad_ubar(i,Jend+1,kout)=0.0_r8
!>            tl_cff2=-cff2*cff2*tl_Ce
!>
              ad_Ce=ad_Ce-cff2*cff2*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Ce=cff*tl_cff1
!>
              ad_cff1=ad_cff1+cff*ad_Ce
              ad_Ce=0.0_r8
!>            tl_cff1=0.25_r8*g*(GRID(ng)%tl_h(i-1,Jend)+               &
!>   &                           tl_zeta(i-1,Jend,know)+                &
!>   &                           GRID(ng)%tl_h(i  ,Jend)+               &
!>   &                           tl_zeta(i  ,Jend,know))/cff1
!>
              adfac=0.25_r8*g*ad_cff1/cff1
              GRID(ng)%ad_h(i-1,Jend)=GRID(ng)%ad_h(i-1,Jend)+adfac
              GRID(ng)%ad_h(i  ,Jend)=GRID(ng)%ad_h(i  ,Jend)+adfac
              ad_zeta(i-1,Jend,know)=ad_zeta(i-1,Jend,know)+adfac
              ad_zeta(i  ,Jend,know)=ad_zeta(i  ,Jend,know)+adfac
              ad_cff1=0.0_r8
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(inorth,isUbar,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
# ifdef MASKING
!>            tl_ubar(i,Jend+1,kout)=tl_ubar(i,Jend+1,kout)*            &
!>   &                               GRID(ng)%umask(i,Jend+1)
!>
              ad_ubar(i,Jend+1,kout)=ad_ubar(i,Jend+1,kout)*            &
     &                               GRID(ng)%umask(i,Jend+1)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isUbar,ng)) THEN
!>              tl_ubar(i,Jend+1,kout)=BOUNDARY(ng)%tl_ubar_north(i)
!>
                BOUNDARY(ng)%ad_ubar_north(i)=BOUNDARY(ng)%             &
     &                                              ad_ubar_north(i)+   &
     &                                        ad_ubar(i,Jend+1,kout)
                ad_ubar(i,Jend+1,kout)=0.0_r8
              ELSE
!>              tl_ubar(i,Jend+1,kout)=0.0_r8
!>
                ad_ubar(i,Jend+1,kout)=0.0_r8
              END IF
# else
!>            tl_ubar(i,Jend+1,kout)=0.0_r8
!>
              ad_ubar(i,Jend+1,kout)=0.0_r8
# endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(inorth,isUbar,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
# ifdef MASKING
!>            tl_ubar(i,Jend+1,kout)=tl_ubar(i,Jend+1,kout)*            &
!>   &                               GRID(ng)%umask(i,Jend+1)
!>
              ad_ubar(i,Jend+1,kout)=ad_ubar(i,Jend+1,kout)*            &
     &                               GRID(ng)%umask(i,Jend+1)
# endif
!>            tl_ubar(i,Jend+1,kout)=tl_ubar(i,Jend,kout)
!>
              ad_ubar(i,Jend,kout)=ad_ubar(i,Jend,kout)+                &
     &                             ad_ubar(i,Jend+1,kout)
              ad_ubar(i,Jend+1,kout)=0.0_r8
            END IF
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (ad_LBC(inorth,isUbar,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%north(i)) THEN
# ifdef MASKING
!>            tl_ubar(i,Jend+1,kout)=tl_ubar(i,Jend+1,kout)*            &
!>   &                               GRID(ng)%umask(i,Jend+1)
!>
              ad_ubar(i,Jend+1,kout)=ad_ubar(i,Jend+1,kout)*            &
     &                               GRID(ng)%umask(i,Jend+1)
# endif
!>            tl_ubar(i,Jend+1,kout)=gamma2(ng)*tl_ubar(i,Jend,kout)
!>
              ad_ubar(i,Jend  ,kout)=ad_ubar(i,Jend,kout)+              &
     &                               gamma2(ng)*ad_ubar(i,Jend+1,kout)
              ad_ubar(i,Jend+1,kout)=0.0_r8
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
        IF (ad_LBC(isouth,isUbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (ad_LBC(isouth,isUbar,ng)%nudging) THEN
                  IF (LnudgeM2CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M2nudgcof(i-1,Jstr-1)+           &
     &                       CLIMA(ng)%M2nudgcof(i  ,Jstr-1))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M2obc_out(ng,isouth)
                    obc_in =M2obc_in (ng,isouth)
                  END IF
                  IF (BOUNDARY(ng)%ubar_south_Ce(i).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt2d
                END IF
#  ifdef RADIATION_2D
                Cx=BOUNDARY(ng)%ubar_south_Cx(i)
#  else
                Cx=0.0_r8
#  endif
                Ce=BOUNDARY(ng)%ubar_south_Ce(i)
                cff=BOUNDARY(ng)%ubar_south_C2(i)
# endif
# ifdef MASKING
!>              tl_ubar(i,Jstr-1,kout)=tl_ubar(i,Jstr-1,kout)*          &
!>   &                                 GRID(ng)%umask(i,Jstr-1)
!>
                ad_ubar(i,Jstr-1,kout)=ad_ubar(i,Jstr-1,kout)*          &
     &                                 GRID(ng)%umask(i,Jstr-1)
# endif
                IF (ad_LBC(isouth,isUbar,ng)%nudging) THEN
!>                tl_ubar(i,Jstr-1,kout)=tl_ubar(i,Jstr-1,kout)-        &
!>   &                                   tau*tl_ubar(i,Jstr-1,know)
!>
                  ad_ubar(i,Jstr-1,know)=ad_ubar(i,Jstr-1,know)-        &
     &                                 tau*ad_ubar(i,Jstr-1,kout)
                END IF
!>              tl_ubar(i,Jstr-1,kout)=(cff*tl_ubar(i,Jstr-1,know)+     &
!>   &                                  Ce *tl_ubar(i,Jstr  ,kout)-     &
!>   &                                  MAX(Cx,0.0_r8)*                 &
!>   &                                     tl_grad(i-1,Jstr-1)-         &
!>   &                                  MIN(Cx,0.0_r8)*                 &
!>   &                                     tl_grad(i  ,Jstr-1))/        &
!>   &                                 (cff+Ce)
!>
                adfac=ad_ubar(i,Jstr-1,kout)/(cff+Ce)
                ad_grad(i-1,Jstr-1)=ad_grad(i-1,Jstr-1)-                &
     &                              MAX(Cx,0.0_r8)*adfac
                ad_grad(i  ,Jstr-1)=ad_grad(i  ,Jstr-1)-                &
     &                              MIN(Cx,0.0_r8)*adfac
                ad_ubar(i,Jstr-1,know)=ad_ubar(i,Jstr-1,know)+cff*adfac
                ad_ubar(i,Jstr  ,kout)=ad_ubar(i,Jstr  ,kout)+Ce *adfac
                ad_ubar(i,Jstr-1,kout)=0.0_r8
              END IF
            END DO
          END IF
!
!  Southern edge, Chapman boundary condition.
!
        ELSE IF (ad_LBC(isouth,isUbar,ng)%Flather.or.                   &
     &           ad_LBC(isouth,isUbar,ng)%Shchepetkin.or.               &
     &           ad_LBC(isouth,isUbar,ng)%reduced) THEN
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
# ifdef MASKING
!>            tl_ubar(i,Jstr-1,kout)=tl_ubar(i,Jstr-1,kout)*            &
!>   &                               GRID(ng)%umask(i,Jstr-1)
!>
              ad_ubar(i,Jstr-1,kout)=ad_ubar(i,Jstr-1,kout)*            &
     &                               GRID(ng)%umask(i,Jstr-1)
# endif
!>            tl_ubar(i,Jstr-1,kout)=tl_cff2*(ubar(i,Jstr-1,know)+      &
!>   &                                        Ce*ubar(i,Jstr,kout))+    &
!>   &                               cff2*(tl_ubar(i,Jstr-1,know)+      &
!>   &                                     tl_Ce*ubar(i,Jstr,kout)+     &
!>   &                                     Ce*tl_ubar(i,Jstr,kout))
!>
              adfac=cff2*ad_ubar(i,Jstr-1,kout)
              ad_ubar(i,Jstr-1,know)=ad_ubar(i,Jstr-1,know)+adfac
              ad_ubar(i,Jstr  ,kout)=ad_ubar(i,Jstr  ,kout)+Ce*adfac
              ad_Ce=ad_Ce+ubar(i,Jstr,kout)*adfac
              ad_cff2=ad_cff2+                                          &
     &                (ubar(i,Jstr-1,know)+                             &
     &                 Ce*ubar(i,Jstr,kout))*ad_ubar(i,Jstr-1,kout)
              ad_ubar(i,Jstr-1,kout)=0.0_r8
!>            tl_cff2=-cff2*cff2*tl_Ce
!>
              ad_Ce=ad_Ce-cff2*cff2*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Ce=cff*tl_cff1
!>
              ad_cff1=ad_cff1+cff*ad_Ce
              ad_Ce=0.0_r8
!>            tl_cff1=0.25_r8*g*(GRID(ng)%tl_h(i-1,Jstr)+               &
!>   &                           tl_zeta(i-1,Jstr,know)+                &
!>   &                           GRID(ng)%tl_h(i  ,Jstr)+               &
!>   &                           tl_zeta(i  ,Jstr,know))/cff1
!>
              adfac=0.25_r8*g*ad_cff1/cff1
              GRID(ng)%ad_h(i-1,Jstr)=GRID(ng)%ad_h(i-1,Jstr)+adfac
              GRID(ng)%ad_h(i  ,Jstr)=GRID(ng)%ad_h(i  ,Jstr)+adfac
              ad_zeta(i-1,Jstr,know)=ad_zeta(i-1,Jstr,know)+adfac
              ad_zeta(i  ,Jstr,know)=ad_zeta(i  ,Jstr,know)+adfac
              ad_cff1=0.0_r8
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(isouth,isUbar,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
# ifdef MASKING
!>            tl_ubar(i,Jstr-1,kout)=tl_ubar(i,Jstr-1,kout)*            &
!>   &                               GRID(ng)%umask(i,Jstr-1)
!>
              ad_ubar(i,Jstr-1,kout)=ad_ubar(i,Jstr-1,kout)*            &
     &                               GRID(ng)%umask(i,Jstr-1)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isUbar,ng)) THEN
!>              tl_ubar(i,Jstr-1,kout)=BOUNDARY(ng)%tl_ubar_south(i)
!>
                BOUNDARY(ng)%ad_ubar_south(i)=BOUNDARY(ng)%             &
     &                                              ad_ubar_south(i)+   &
     &                                        ad_ubar(i,Jstr-1,kout)
                ad_ubar(i,Jstr-1,kout)=0.0_r8
              ELSE
!>              tl_ubar(i,Jstr-1,kout)=0.0_r8
!>
                ad_ubar(i,Jstr-1,kout)=0.0_r8
              END IF
# else
!>            tl_ubar(i,Jstr-1,kout)=0.0_r8
!>
              ad_ubar(i,Jstr-1,kout)=0.0_r8
# endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(isouth,isUbar,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
# ifdef MASKING
!>            tl_ubar(i,Jstr-1,kout)=tl_ubar(i,Jstr-1,kout)*            &
!>   &                               GRID(ng)%umask(i,Jstr-1)
!>
              ad_ubar(i,Jstr-1,kout)=ad_ubar(i,Jstr-1,kout)*            &
     &                               GRID(ng)%umask(i,Jstr-1)
# endif
!>            tl_ubar(i,Jstr-1,kout)=tl_ubar(i,Jstr,kout)
!>
              ad_ubar(i,Jstr  ,kout)=ad_ubar(i,Jstr  ,kout)+            &
     &                               ad_ubar(i,Jstr-1,kout)
              ad_ubar(i,Jstr-1,kout)=0.0_r8
            END IF
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (ad_LBC(isouth,isUbar,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%south(i)) THEN
# ifdef MASKING
!>            tl_ubar(i,Jstr-1,kout)=tl_ubar(i,Jstr-1,kout)*            &
!>   &                               GRID(ng)%umask(i,Jstr-1)
!>
              ad_ubar(i,Jstr-1,kout)=ad_ubar(i,Jstr-1,kout)*            &
     &                               GRID(ng)%umask(i,Jstr-1)
# endif
!>            tl_ubar(i,Jstr-1,kout)=gamma2(ng)*tl_ubar(i,Jstr,kout)
!>
              ad_ubar(i,Jstr  ,kout)=ad_ubar(i,Jstr,kout)+              &
     &                               gamma2(ng)*ad_ubar(i,Jstr-1,kout)
              ad_ubar(i,Jstr-1,kout)=0.0_r8
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
        IF (ad_LBC(ieast,isUbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (ad_LBC(ieast,isUbar,ng)%nudging) THEN
                  IF (LnudgeM2CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M2nudgcof(Iend  ,j)+             &
     &                       CLIMA(ng)%M2nudgcof(Iend+1,j))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M2obc_out(ng,ieast)
                    obc_in =M2obc_in (ng,ieast)
                  END IF
                  IF (BOUNDARY(ng)%ubar_east_Cx(j).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt2d
                END IF
                Cx=BOUNDARY(ng)%ubar_east_Cx(j)
#  ifdef RADIATION_2D
                Ce=BOUNDARY(ng)%ubar_east_Ce(j)
#  else
                Ce=0.0_r8
#  endif
                cff=BOUNDARY(ng)%ubar_east_C2(j)
# endif
# ifdef MASKING
!>              tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,kout)*          &
!>   &                                 GRID(ng)%umask(Iend+1,j)
!>
                ad_ubar(Iend+1,j,kout)=ad_ubar(Iend+1,j,kout)*          &
     &                                 GRID(ng)%umask(Iend+1,j)
# endif
                IF (ad_LBC(ieast,isUbar,ng)%nudging) THEN
!>                tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,kout)-        &
!>   &                                   tau*tl_ubar(Iend+1,j,know)
!>
                  ad_ubar(Iend+1,j,know)=ad_ubar(Iend+1,j,know)-        &
     &                                   tau*ad_ubar(Iend+1,j,kout)
                END IF
!>              tl_ubar(Iend+1,j,kout)=(cff*tl_ubar(Iend+1,j,know)+     &
!>   &                                  Cx *tl_ubar(Iend  ,j,kout)-     &
!>   &                                  MAX(Ce,0.0_r8)*                 &
!>   &                                     tl_grad(Iend+1,j  )-         &
!>   &                                  MIN(Ce,0.0_r8)*                 &
!>   &                                     tl_grad(Iend+1,j+1))/        &
!>   &                                 (cff+Cx)
!>
                adfac=ad_ubar(Iend+1,j,kout)/(cff+Cx)
                ad_grad(Iend+1,j  )=ad_grad(Iend+1,j  )-                &
     &                              MAX(Ce,0.0_r8)*adfac
                ad_grad(Iend+1,j+1)=ad_grad(Iend+1,j+1)-                &
     &                              MIN(Ce,0.0_r8)*adfac
                ad_ubar(Iend  ,j,kout)=ad_ubar(Iend  ,j,kout)+Cx *adfac
                ad_ubar(Iend+1,j,know)=ad_ubar(Iend+1,j,know)+cff*adfac
                ad_ubar(Iend+1,j,kout)=0.0_r8
              END IF
            END DO
          END IF
!
!  Eastern edge, Flather boundary condition.
!
        ELSE IF (ad_LBC(ieast,isUbar,ng)%Flather) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                 &
     &                            zeta(Iend  ,j,know)+                  &
     &                            GRID(ng)%h(Iend+1,j)+                 &
     &                            zeta(Iend+1,j,know)))
              Cx=SQRT(g*cff)
# ifdef MASKING
!>            tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,kout)*            &
!>   &                               GRID(ng)%umask(Iend+1,j)
!>
              ad_ubar(Iend+1,j,kout)=ad_ubar(Iend+1,j,kout)*            &
     &                               GRID(ng)%umask(Iend+1,j)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(ieast,isUbar,ng)) THEN
!>              tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,kout)-          &
!>   &                                 Cx*BOUNDARY(ng)%tl_zeta_east(j)
!>
                BOUNDARY(ng)%ad_zeta_east(j)=                           &
     &                                   BOUNDARY(ng)%ad_zeta_east(j)-  &
     &                                       Cx*ad_ubar(Iend+1,j,kout)
              END IF
# endif
!>            tl_ubar(Iend+1,j,kout)=tl_bry_val+                        &
!>   &                               tl_Cx*                             &
!>   &                               (0.5_r8*(zeta(Iend  ,j,know)+      &
!>   &                                        zeta(Iend+1,j,know))-     &
!>   &                                BOUNDARY(ng)%zeta_east(j))+       &
!>   &                               Cx*                                &
!>                                   (0.5_r8*(tl_zeta(Iend  ,j,know)+   &
!>   &                                        tl_zeta(Iend+1,j,know)))
!>
              adfac=Cx*0.5_r8*ad_ubar(Iend+1,j,kout)
              ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)+adfac
              ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
              ad_Cx=ad_Cx+                                              &
     &              (0.5_r8*(zeta(Iend  ,j,know)+                       &
     &                       zeta(Iend+1,j,know))-                      &
     &               BOUNDARY(ng)%zeta_east(j))*ad_ubar(Iend+1,j,kout)
              ad_bry_val=ad_bry_val+ad_ubar(Iend+1,j,kout)
              ad_ubar(Iend+1,j,kout)=0.0_r8
!>            tl_Cx=0.5_r8*g*tl_cff/Cx
!>
              ad_cff=ad_cff+0.5_r8*g*ad_Cx/Cx
              ad_Cx=0.0_r8
!>            tl_cff=-cff*cff*(0.5_r8*(GRID(ng)%tl_h(Iend  ,j)+         &
!>   &                                 tl_zeta(Iend  ,j,know)+          &
!>   &                                 GRID(ng)%tl_h(Iend+1,j)+         &
!>   &                                 tl_zeta(Iend+1,j,know)))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              GRID(ng)%ad_h(Iend  ,j)=GRID(ng)%ad_h(Iend  ,j)+adfac
              GRID(ng)%ad_h(Iend+1,j)=GRID(ng)%ad_h(Iend+1,j)+adfac
              ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)+adfac
              ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
              ad_cff=0.0_r8

# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (ad_LBC(ieast,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_east(j)-                  &
     &                      zeta(Iend,j,know))*                         &
     &                  0.5_r8*GRID(ng)%pm(Iend,j)
              ELSE
                bry_pgr=-g*(zeta(Iend+1,j,know)-                        &
     &                      zeta(Iend  ,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Iend  ,j)+                  &
     &                          GRID(ng)%pm(Iend+1,j))
              END IF
#  ifdef UV_COR
              bry_cor=0.125_r8*(vbar(Iend  ,j  ,know)+                  &
     &                          vbar(Iend  ,j+1,know)+                  &
     &                          vbar(Iend+1,j  ,know)+                  &
     &                          vbar(Iend+1,j+1,know))*                 &
     &                         (GRID(ng)%f(Iend  ,j)+                   &
     &                          GRID(ng)%f(Iend+1,j))
#  else
              bry_cor=0.0_r8
#  endif
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
!>            tl_bry_val=tl_ubar(Iend,j,know)+                          &
!>   &                   tl_cff2*(bry_pgr+                              &
!>   &                            bry_cor+                              &
!>   &                            bry_str)+                             &
!>   &                   cff2*(tl_bry_pgr+                              &
!>   &                         tl_bry_cor+                              &
!>   &                         tl_bry_str)
!>
              adfac=cff2*ad_bry_val
              ad_bry_pgr=ad_bry_pgr+adfac
              ad_bry_cor=ad_bry_cor+adfac
              ad_bry_str=ad_bry_str+adfac
              ad_cff2=ad_cff2+(bry_pgr+                                 &
     &                         bry_cor+                                 &
     &                         bry_str)*ad_bry_val
              ad_ubar(Iend,j,know)=ad_ubar(Iend,j,know)+ad_bry_val
              ad_bry_val=0.0_r8
!>            tl_cff2=GRID(ng)%om_u(Iend+1,j)*tl_Cx
!>
              ad_Cx=ad_Cx+GRID(ng)%om_u(Iend+1,j)*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Cx=-Cx*Cx*Cx*0.25_r8*g*(GRID(ng)%tl_h(Iend+1,j)+       &
!>   &                                   tl_zeta(Iend+1,j,know)+        &
!>   &                                   GRID(ng)%tl_h(Iend  ,j)+       &
!>   &                                   tl_zeta(Iend  ,j,know))
!>
              adfac=-Cx*Cx*Cx*0.25_r8*g*ad_Cx
              ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)+adfac
              ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
              GRID(ng)%ad_h(Iend  ,j)=GRID(ng)%ad_h(Iend  ,j)+adfac
              GRID(ng)%ad_h(Iend+1,j)=GRID(ng)%ad_h(Iend+1,j)+adfac
              ad_Cx=0.0_r8
!>            tl_bry_str=tl_cff1*(FORCES(ng)%sustr(Iend+1,j)-           &
!>   &                            FORCES(ng)%bustr(Iend+1,j))+          &
!>   &                   cff1*(FORCES(ng)%tl_sustr(Iend+1,j)-           &
!>   &                         FORCES(ng)%tl_bustr(Iend+1,j))
!>
              adfac=cff1*ad_bry_str
              FORCES(ng)%ad_sustr(Iend+1,j)=                            &
     &                                   FORCES(ng)%ad_sustr(Iend+1,j)+ &
     &                                      adfac
              FORCES(ng)%tl_bustr(Iend+1,j)=                            &
     &                                   FORCES(ng)%tl_bustr(Iend+1,j)- &
     &                                      adfac
              ad_cff1=ad_cff1+(FORCES(ng)%sustr(Iend+1,j)-              &
     &                         FORCES(ng)%bustr(Iend+1,j))*ad_bry_str
              ad_bry_str=0.0_r8
!>            tl_cff1=-cff1*cff1*0.5_r8*(GRID(ng)%tl_h(Iend  ,j)+       &
!>   &                                   tl_zeta(Iend  ,j,know)+        &
!>   &                                   GRID(ng)%tl_h(Iend+1,j)+       &
!>   &                                   tl_zeta(Iend+1,j,know))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)+adfac
              ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
              GRID(ng)%ad_h(Iend  ,j)=GRID(ng)%ad_h(Iend  ,j)+adfac
              GRID(ng)%ad_h(Iend+1,j)=GRID(ng)%ad_h(Iend+1,j)+adfac
              ad_cff1=0.0_r8
#  ifdef UV_COR
!>            tl_bry_cor=0.125_r8*(tl_vbar(Iend  ,j  ,know)+            &
!>   &                             tl_vbar(Iend  ,j+1,know)+            &
!>   &                             tl_vbar(Iend+1,j  ,know)+            &
!>   &                             tl_vbar(Iend+1,j+1,know))*           &
!>   &                            (GRID(ng)%f(Iend  ,j)+                &
!>   &                             GRID(ng)%f(Iend+1,j))
!>
              adfac=0.125_r8*(GRID(ng)%f(Iend  ,j)+                     &
     &                        GRID(ng)%f(Iend+1,j))*ad_bry_cor
              ad_vbar(Iend  ,j  ,know)=ad_vbar(Iend  ,j  ,know)+adfac
              ad_vbar(Iend+1,j  ,know)=ad_vbar(Iend+1,j  ,know)+adfac
              ad_vbar(Iend  ,j+1,know)=ad_vbar(Iend  ,j+1,know)+adfac
              ad_vbar(Iend+1,j+1,know)=ad_vbar(Iend+1,j+1,know)+adfac
              ad_bry_cor=0.0_r8
#  else
!>            tl_bry_cor=0.0_r8
!>
              ad_bry_cor=0.0_r8
#  endif
              IF (ad_LBC(ieast,isFsur,ng)%acquire) THEN
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(ieast,isUbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr-                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_east(j)*            &
!>   &                       0.5_r8*GRID(ng)%pm(Iend,j)
!>
                  BOUNDARY(ng)%ad_zeta_east(j)=BOUNDARY(ng)%            &
     &                                               ad_zeta_east(j)-   &
     &                                         g*0.5_r8*                &
     &                                         GRID(ng)%pm(Iend,j)*     &
     &                                         ad_bry_pgr
                END IF
#  endif
!>              tl_bry_pgr=g*tl_zeta(Iend,j,know)*                      &
!>   &                     0.5_r8*GRID(ng)%pm(Iend,j)
!>
                ad_zeta(Iend,j,know)=ad_zeta(Iend,j,know)+              &
     &                               g*0.5_r8*GRID(ng)%pm(Iend,j)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(Iend+1,j,know)-                  &
!>   &                         tl_zeta(Iend  ,j,know))*                 &
!>   &                      0.5_r8*(GRID(ng)%pm(Iend  ,j)+              &
!>   &                              GRID(ng)%pm(Iend+1,j))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pm(Iend  ,j)+                 &
     &                           GRID(ng)%pm(Iend+1,j))*ad_bry_pgr
                ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)-adfac
                ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
# else
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(ieast,isUbar,ng)) THEN
!>              tl_bry_val=BOUNDARY(ng)%tl_ubar_east(j)
!>
                BOUNDARY(ng)%ad_ubar_east(j)=                           &
     &                                   BOUNDARY(ng)%ad_ubar_east(j)+  &
     &                                       ad_bry_val
                ad_bry_val=0.0_r8
              ELSE
!>              tl_bry_val=0.0_r8
!>
                ad_bry_val=0.0_r8
              END IF
#  else
!>            tl_bry_val=0.0_r8
!>
              ad_bry_val=0.0_r8
#  endif
# endif
            END IF
          END DO
!
!  Eastern edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (ad_LBC(ieast,isUbar,ng)%Shchepetkin) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=0.5_r8*(GRID(ng)%h(Iend  ,j)+                         &
     &                    GRID(ng)%h(Iend+1,j))
              cff1=SQRT(g/cff)
              Cx=dt2d*cff1*cff*0.5_r8*(GRID(ng)%pm(Iend  ,j)+           &
     &                                 GRID(ng)%pm(Iend+1,j))
              Zx=(0.5_r8+Cx)*zeta(Iend  ,j,know)+                       &
     &           (0.5_r8-Cx)*zeta(Iend+1,j,know)
              IF (Cx.gt.Co) THEN
                cff2=(1.0_r8-Co/Cx)**2
                cff3=zeta(Iend,j,kout)+                                 &
     &               Cx*zeta(Iend+1,j,know)-                            &
     &               (1.0_r8+Cx)*zeta(Iend,j,know)
                Zx=Zx+cff2*cff3
              END IF
# ifdef MASKING
!>            tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,kout)*            &
!>   &                               GRID(ng)%umask(Iend+1,j)
!>
              ad_ubar(Iend+1,j,kout)=ad_ubar(Iend+1,j,kout)*            &
     &                               GRID(ng)%umask(Iend+1,j)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(ieast,isUbar,ng)) THEN
!>              tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,kout)-          &
!>   &                                 0.5_r8*cff1*                     &
!>   &                                 BOUNDARY(ng)%tl_zeta_east(j)
!>
                BOUNDARY(ng)%ad_zeta_east(j)=                           &
     &                                   BOUNDARY(ng)%ad_zeta_east(j)-  &
     &                                       0.5_r8*cff1*               &
     &                                       ad_ubar(Iend+1,j,kout)
              END IF
# endif
!>            tl_ubar(Iend+1,j,kout)=0.5_r8*                            &
!>   &                               ((1.0_r8-Cx)*                      &
!>   &                                tl_ubar(Iend+1,j,know)+           &
!>   &                                tl_Cx*(ubar(Iend  ,j,know)-       &
!>   &                                       ubar(Iend+1,j,know))+      &
!>   &                                Cx*tl_ubar(Iend,j,know)+          &
!>   &                                tl_bry_val+                       &
!>   &                                tl_cff1*                          &
!>   &                                (Zx-BOUNDARY(ng)%zeta_east(j))-   &
!>   &                                cff1*tl_Zx)
!>
              adfac=0.5_r8*ad_ubar(Iend+1,j,kout)
              ad_ubar(Iend+1,j,know)=ad_ubar(Iend+1,j,know)+            &
     &                               (1.0_r8-Cx)*adfac
              ad_ubar(Iend  ,j,know)=ad_ubar(Iend  ,j,know)+            &
     &                               Cx*adfac
              ad_Cx=ad_Cx+                                              &
     &              (ubar(Iend  ,j,know)-                               &
     &               ubar(Iend+1,j,know))*adfac
              ad_bry_val=ad_bry_val+adfac
              ad_cff1=ad_cff1+                                          &
     &                (Zx-BOUNDARY(ng)%zeta_east(j))*adfac
              ad_Zx=ad_Zx-cff1*adfac
              ad_ubar(Iend+1,j,kout)=0.0_r8
              IF (Cx.gt.Co) THEN
!>              tl_Zx=tl_Zx+cff2*tl_cff3+                               &
!>   &                      tl_cff2*cff3
!>
                ad_cff2=ad_cff2+cff3*ad_Zx
                ad_cff3=ad_cff3+cff2*ad_Zx
!>              tl_cff3=tl_zeta(Iend,j,kout)+                           &
!>   &                  Cx*tl_zeta(Iend+1,j,know)+                      &
!>   &                  tl_Cx*(zeta(Iend  ,j,know)+                     &
!>   &                         zeta(Iend+1,j,know))-                    &
!>   &                  (1.0_r8+Cx)*tl_zeta(Iend,j,know)
!>
                ad_zeta(Iend  ,j,kout)=ad_zeta(Iend  ,j,kout)+          &
     &                                 ad_cff3
                ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)-          &
     &                                 (1.0_r8+Cx)*ad_cff3
                ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+          &
     &                                 Cx*ad_cff3
                ad_Cx=ad_Cx+                                            &
     &                (zeta(Iend  ,j,know)+                             &
     &                 zeta(Iend+1,j,know))*ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=2.0_r8*cff2*Co*tl_Cx/(Cx*Cx)
!>
                ad_Cx=ad_Cx+                                            &
     &                2.0_r8*cff2*Co*ad_cff2/(Cx*Cx)
                ad_cff2=0.0_r8
              END IF
!>            tl_Zx=(0.5_r8+Cx)*tl_zeta(Iend  ,j,know)+                 &
!>   &              (0.5_r8-Cx)*tl_zeta(Iend+1,j,know)+                 &
!>   &              tl_Cx*(zeta(Iend  ,j,know)-                         &
!>   &                     zeta(Iend+1,j,know))
!>
              ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)+            &
     &                               (0.5_r8+Cx)*ad_Zx
              ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+            &
     &                               (0.5_r8-Cx)*ad_Zx
              ad_Cx=ad_Cx+                                              &
     &              (zeta(Iend  ,j,know)-                               &
     &               zeta(Iend+1,j,know))*ad_Zx
              ad_Zx=0.0_r8
!>            tl_Cx=dt2d*0.5_r8*(GRID(ng)%pm(Iend  ,j)+                 &
!>   &                           GRID(ng)%pm(Iend+1,j))*                &
!>   &              (cff1*tl_cff+                                       &
!>   &               tl_cff1*cff)
!>
              adfac=dt2d*0.5_r8*(GRID(ng)%pm(Iend  ,j)+                 &
     &                           GRID(ng)%pm(Iend+1,j))*ad_Cx
              ad_cff=ad_cff+cff1*adfac
              ad_cff1=ad_cff1+cff*adfac
              ad_Cx=0.0_r8
!>            tl_cff1=-0.5_r8*cff1*tl_cff/cff
!>
              ad_cff=ad_cff-                                            &
     &               0.5_r8*cff1*ad_cff/cff
              ad_cff1=0.0_r8
!>            tl_cff=0.5_r8*(GRID(ng)%tl_h(Iend  ,j)+                   &
!>   &                       GRID(ng)%tl_h(Iend+1,j))
!>
              adfac=0.5_r8*ad_cff
              GRID(ng)%ad_h(Iend  ,j)=GRID(ng)%ad_h(Iend  ,j)+adfac
              GRID(ng)%ad_h(Iend+1,j)=GRID(ng)%ad_h(Iend+1,j)+adfac
              ad_cff=0.0_r8

# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (ad_LBC(ieast,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_east(j)-                  &
     &                      zeta(Iend,j,know))*                         &
     &                  0.5_r8*GRID(ng)%pm(Iend,j)
              ELSE
                bry_pgr=-g*(zeta(Iend+1,j,know)-                        &
     &                      zeta(Iend  ,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Iend  ,j)+                  &
     &                          GRID(ng)%pm(Iend+1,j))
              END IF
#  ifdef UV_COR
              bry_cor=0.125_r8*(vbar(Iend  ,j  ,know)+                  &
     &                          vbar(Iend  ,j+1,know)+                  &
     &                          vbar(Iend+1,j  ,know)+                  &
     &                          vbar(Iend+1,j+1,know))*                 &
     &                         (GRID(ng)%f(Iend  ,j)+                   &
     &                          GRID(ng)%f(Iend+1,j))
#  else
              bry_cor=0.0_r8
#  endif
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
!>            tl_bry_val=tl_ubar(Iend,j,know)+                          &
!>   &                   tl_cff2*(bry_pgr+                              &
!>   &                            bry_cor+                              &
!>   &                            bry_str)+                             &
!>   &                   cff2*(tl_bry_pgr+                              &
!>   &                         tl_bry_cor+                              &
!>   &                         tl_bry_str)
!>
              adfac=cff2*ad_bry_val
              ad_bry_pgr=ad_bry_pgr+adfac
              ad_bry_cor=ad_bry_cor+adfac
              ad_bry_str=ad_bry_str+adfac
              ad_cff2=ad_cff2+(bry_pgr+                                 &
     &                         bry_cor+                                 &
     &                         bry_str)*ad_bry_val
              ad_ubar(Iend,j,know)=ad_ubar(Iend,j,know)+ad_bry_val
              ad_bry_val=0.0_r8
!>            tl_cff2=GRID(ng)%om_u(Iend+1,j)*tl_Cx
!>
              ad_Cx=ad_Cx+GRID(ng)%om_u(Iend+1,j)*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Cx=-Cx*Cx*Cx*0.25_r8*g*(GRID(ng)%tl_h(Iend+1,j)+       &
!>   &                                   tl_zeta(Iend+1,j,know)+        &
!>   &                                   GRID(ng)%tl_h(Iend  ,j)+       &
!>   &                                   tl_zeta(Iend  ,j,know))
!>
              adfac=-Cx*Cx*Cx*0.25_r8*g*ad_Cx
              ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)+adfac
              ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
              GRID(ng)%ad_h(Iend  ,j)=GRID(ng)%ad_h(Iend  ,j)+adfac
              GRID(ng)%ad_h(Iend+1,j)=GRID(ng)%ad_h(Iend+1,j)+adfac
              ad_Cx=0.0_r8
!>            tl_bry_str=tl_cff1*(FORCES(ng)%sustr(Iend+1,j)-           &
!>   &                            FORCES(ng)%bustr(Iend+1,j))+          &
!>   &                   cff1*(FORCES(ng)%tl_sustr(Iend+1,j)-           &
!>   &                         FORCES(ng)%tl_bustr(Iend+1,j))
!>
              adfac=cff1*ad_bry_str
              FORCES(ng)%ad_sustr(Iend+1,j)=                            &
     &                                   FORCES(ng)%ad_sustr(Iend+1,j)+ &
     &                                      adfac
              FORCES(ng)%tl_bustr(Iend+1,j)=                            &
     &                                   FORCES(ng)%tl_bustr(Iend+1,j)- &
     &                                      adfac
              ad_cff1=ad_cff1+(FORCES(ng)%sustr(Iend+1,j)-              &
     &                         FORCES(ng)%bustr(Iend+1,j))*ad_bry_str
              ad_bry_str=0.0_r8
!>            tl_cff1=-cff1*cff1*0.5_r8*(GRID(ng)%tl_h(Iend  ,j)+       &
!>   &                                   tl_zeta(Iend  ,j,know)+        &
!>   &                                   GRID(ng)%tl_h(Iend+1,j)+       &
!>   &                                   tl_zeta(Iend+1,j,know))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)+adfac
              ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
              GRID(ng)%ad_h(Iend  ,j)=GRID(ng)%ad_h(Iend  ,j)+adfac
              GRID(ng)%ad_h(Iend+1,j)=GRID(ng)%ad_h(Iend+1,j)+adfac
              ad_cff1=0.0_r8
#  ifdef UV_COR
!>            tl_bry_cor=0.125_r8*(tl_vbar(Iend  ,j  ,know)+            &
!>   &                             tl_vbar(Iend  ,j+1,know)+            &
!>   &                             tl_vbar(Iend+1,j  ,know)+            &
!>   &                             tl_vbar(Iend+1,j+1,know))*           &
!>   &                            (GRID(ng)%f(Iend  ,j)+                &
!>   &                             GRID(ng)%f(Iend+1,j))
!>
              adfac=0.125_r8*(GRID(ng)%f(Iend  ,j)+                     &
     &                        GRID(ng)%f(Iend+1,j))*ad_bry_cor
              ad_vbar(Iend  ,j  ,know)=ad_vbar(Iend  ,j  ,know)+adfac
              ad_vbar(Iend+1,j  ,know)=ad_vbar(Iend+1,j  ,know)+adfac
              ad_vbar(Iend  ,j+1,know)=ad_vbar(Iend  ,j+1,know)+adfac
              ad_vbar(Iend+1,j+1,know)=ad_vbar(Iend+1,j+1,know)+adfac
              ad_bry_cor=0.0_r8
#  else
!>            tl_bry_cor=0.0_r8
!>
              ad_bry_cor=0.0_r8
#  endif
              IF (ad_LBC(ieast,isFsur,ng)%acquire) THEN
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(ieast,isUbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr-                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_east(j)*            &
!>   &                       0.5_r8*GRID(ng)%pm(Iend,j)
!>
                  BOUNDARY(ng)%ad_zeta_east(j)=BOUNDARY(ng)%            &
     &                                               ad_zeta_east(j)-   &
     &                                         g*0.5_r8*                &
     &                                         GRID(ng)%pm(Iend,j)*     &
     &                                         ad_bry_pgr
                END IF
#  endif
!>              tl_bry_pgr=g*tl_zeta(Iend,j,know)*                      &
!>   &                     0.5_r8*GRID(ng)%pm(Iend,j)
!>
                ad_zeta(Iend,j,know)=ad_zeta(Iend,j,know)+              &
     &                               g*0.5_r8*GRID(ng)%pm(Iend,j)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(Iend+1,j,know)-                  &
!>   &                         tl_zeta(Iend  ,j,know))*                 &
!>   &                      0.5_r8*(GRID(ng)%pm(Iend  ,j)+              &
!>   &                              GRID(ng)%pm(Iend+1,j))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pm(Iend  ,j)+                 &
     &                           GRID(ng)%pm(Iend+1,j))*ad_bry_pgr
                ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)-adfac
                ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
# else
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(ieast,isUbar,ng)) THEN
!>              tl_bry_val=BOUNDARY(ng)%tl_ubar_east(j)
!>
                BOUNDARY(ng)%ad_ubar_east(j)=                           &
     &                                   BOUNDARY(ng)%ad_ubar_east(j)+  &
     &                                       ad_bry_val
                ad_bry_val=0.0_r8
              ELSE
!>              tl_bry_val=0.0_r8
!>
                ad_bry_val=0.0_r8
              END IF
#  else
!>            tl_bry_val=0.0_r8
!>
              ad_bry_val=0.0_r8
#  endif
# endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(ieast,isUbar,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
# ifdef MASKING
!>            tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,kout)*            &
!>   &                               GRID(ng)%umask(Iend+1,j)
!>
              ad_ubar(Iend+1,j,kout)=ad_ubar(Iend+1,j,kout)*            &
     &                               GRID(ng)%umask(Iend+1,j)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(ieast,isUbar,ng)) THEN
!>              tl_ubar(Iend+1,j,kout)=BOUNDARY(ng)%tl_ubar_east(j)
!>
                BOUNDARY(ng)%ad_ubar_east(j)=                           &
     &                                   BOUNDARY(ng)%ad_ubar_east(j)+  &
     &                                       ad_ubar(Iend+1,j,kout)
                ad_ubar(Iend+1,j,kout)=0.0_r8
              ELSE
!>              tl_ubar(Iend+1,j,kout)=0.0_r8
!>
                ad_ubar(Iend+1,j,kout)=0.0_r8
              END IF
# else
!>            tl_ubar(Iend+1,j,kout)=0.0_r8
!>
              ad_ubar(Iend+1,j,kout)=0.0_r8
# endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(ieast,isUbar,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
# ifdef MASKING
!>            tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,kout)*            &
!>   &                               GRID(ng)%umask(Iend+1,j)
!>
              ad_ubar(Iend+1,j,kout)=ad_ubar(Iend+1,j,kout)*            &
     &                               GRID(ng)%umask(Iend+1,j)
# endif
!>            tl_ubar(Iend+1,j,kout)=tl_ubar(Iend,j,kout)
!>
              ad_ubar(Iend,j,kout)=ad_ubar(Iend,j,kout)+                &
     &                             ad_ubar(Iend+1,j,kout)
              ad_ubar(Iend+1,j,kout)=0.0_r8
            END IF
          END DO
!
!  Eastern edge, reduced-physics boundary condition.
!
        ELSE IF (ad_LBC(ieast,isUbar,ng)%reduced) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                 &
     &                            zeta(Iend  ,j,know)+                  &
     &                            GRID(ng)%h(Iend+1,j)+                 &
     &                            zeta(Iend+1,j,know)))
# ifdef MASKING
!>            tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,kout)*            &
!>   &                               GRID(ng)%umask(Iend+1,j)
!>
              ad_ubar(Iend+1,j,kout)=ad_ubar(Iend+1,j,kout)*            &
     &                               GRID(ng)%umask(Iend+1,j)
# endif
!>            tl_ubar(Iend+1,j,kout)=tl_ubar(Iend+1,j,know)+            &
!>   &                               dt2d*(tl_bry_pgr+                  &
!>   &                                     tl_bry_cor+                  &
!>   &                                     tl_bry_str)
!>
              adfac=dt2d*ad_ubar(Iend+1,j,kout)
              ad_bry_pgr=ad_bry_pgr+adfac
              ad_bry_cor=ad_bry_cor+adfac
              ad_bry_str=ad_bry_str+adfac
              ad_ubar(Iend+1,j,know)=ad_ubar(Iend+1,j,know)+            &
     &                               ad_ubar(Iend+1,j,kout)
              ad_ubar(Iend+1,j,kout)=0.0_r8
!>            tl_bry_str=tl_cff*(FORCES(ng)%sustr(Iend+1,j)-            &
!>   &                           FORCES(ng)%bustr(Iend+1,j))+           &
!>   &                   cff*(FORCES(ng)%tl_sustr(Iend+1,j)-            &
!>   &                        FORCES(ng)%tl_bustr(Iend+1,j))
!>
              adfac=cff*ad_bry_str
              FORCES(ng)%ad_sustr(Iend+1,j)=                            &
     &                                   FORCES(ng)%ad_sustr(Iend+1,j)+ &
     &                                      adfac
              FORCES(ng)%ad_bustr(Iend+1,j)=                            &
     &                                   FORCES(ng)%ad_bustr(Iend+1,j)- &
     &                                      adfac
              ad_cff=ad_cff+(FORCES(ng)%sustr(Iend+1,j)-                &
     &                       FORCES(ng)%bustr(Iend+1,j))*ad_bry_str
              ad_bry_str=0.0_r8
# ifdef UV_COR
!>            tl_bry_cor=0.125_r8*(tl_vbar(Iend,  j  ,know)+            &
!>   &                             tl_vbar(Iend  ,j+1,know)+            &
!>   &                             tl_vbar(Iend+1,j  ,know)+            &
!>   &                             tl_vbar(Iend+1,j+1,know))*           &
!>   &                            (GRID(ng)%f(Iend  ,j)+                &
!>   &                             GRID(ng)%f(Iend+1,j))
!>
              adfac=0.125_r8*(GRID(ng)%f(Iend,j)+                       &
     &                        GRID(ng)%f(Iend+1,j))*ad_bry_cor
              ad_vbar(Iend,  j  ,know)=ad_vbar(Iend  ,j  ,know)+adfac
              ad_vbar(Iend+1,j  ,know)=ad_vbar(Iend+1,j  ,know)+adfac
              ad_vbar(Iend  ,j+1,know)=ad_vbar(Iend  ,j+1,know)+adfac
              ad_vbar(Iend+1,j+1,know)=ad_vbar(Iend+1,j+1,know)+adfac
              ad_bry_cor=0.0_r8
# else
!>            tl_bry_cor=0.0_r8
!>
              ad_bry_cor=0.0_r8
# endif
              IF (ad_LBC(ieast,isFsur,ng)%acquire) THEN
# ifdef ADJUST_BOUNDARY
                IF (Lobc(ieast,isUbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr-                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_east(j)*            &
!>   &                       0.5_r8*GRID(ng)%pm(Iend,j)
!>
                  BOUNDARY(ng)%ad_zeta_east(j)=BOUNDARY(ng)%            &
     &                                               ad_zeta_east(j)-   &
     &                                         g*0.5_r8*                &
     &                                         GRID(ng)%pm(Iend,j)*     &
     &                                         ad_bry_pgr
                END IF
# endif
!>              tl_bry_pgr=g*tl_zeta(Iend,j,know)*                      &
!>   &                     0.5_r8*GRID(ng)%pm(Iend,j)
!>
                ad_zeta(Iend,j,know)=ad_zeta(Iend,j,know)+              &
     &                               g*0.5_r8*GRID(ng)%pm(Iend,j)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(Iend+1,j,know)-                  &
!>   &                         tl_zeta(Iend  ,j,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pm(Iend  ,j)+               &
!>   &                             GRID(ng)%pm(Iend+1,j))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pm(Iend  ,j)+                 &
     &                           GRID(ng)%pm(Iend+1,j))*ad_bry_pgr
                ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)-adfac
                ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (ad_LBC(ieast,isUbar,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
!>            tl_ubar(Iend+1,j,kout)=0.0_r8
!>
              ad_ubar(Iend+1,j,kout)=0.0_r8
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
        IF (ad_LBC(iwest,isUbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (ad_LBC(iwest,isUbar,ng)%nudging) THEN
                  IF (LnudgeM2CLM(ng)) THEN
                    obc_out=0.5_r8*                                     &
     &                      (CLIMA(ng)%M2nudgcof(Istr-1,j)+             &
     &                       CLIMA(ng)%M2nudgcof(Istr  ,j))
                    obc_in =obcfac(ng)*obc_out
                  ELSE
                    obc_out=M2obc_out(ng,iwest)
                    obc_in =M2obc_in (ng,iwest)
                  END IF
                  IF (BOUNDARY(ng)%ubar_west_Cx(j).lt.0.0_r8) THEN
                    tau=obc_in
                  ELSE
                    tau=obc_out
                  END IF
                  tau=tau*dt2d
                END IF
                Cx=BOUNDARY(ng)%ubar_west_Cx(j)
#  ifdef RADIATION_2D
                Ce=BOUNDARY(ng)%ubar_west_Ce(j)
#  else
                Ce=0.0_r8
#  endif
                cff=BOUNDARY(ng)%ubar_west_C2(j)
# endif
# ifdef MASKING
!>              tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,kout)*              &
!>   &                               GRID(ng)%umask(Istr,j)
!>
                ad_ubar(Istr,j,kout)=ad_ubar(Istr,j,kout)*              &
     &                               GRID(ng)%umask(Istr,j)
# endif
                IF (ad_LBC(iwest,isUbar,ng)%nudging) THEN
!>                tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,kout)-            &
!>   &                                 tau*tl_ubar(Istr,j,know)
!>
                  ad_ubar(Istr,j,know)=ad_ubar(Istr,j,know)-            &
     &                                 tau*ad_ubar(Istr,j,kout)
                END IF
!>              tl_ubar(Istr,j,kout)=(cff*tl_ubar(Istr  ,j,know)+       &
!>   &                                Cx *tl_ubar(Istr+1,j,kout)-       &
!>   &                                MAX(Ce,0.0_r8)*                   &
!>   &                                   tl_grad(Istr,j  )-             &
!>   &                                MIN(Ce,0.0_r8)*                   &
!>   &                                   tl_grad(Istr,j+1))/            &
!>   &                               (cff+Cx)
!>
                adfac=ad_ubar(Istr,j,kout)/(cff+Cx)
                ad_grad(Istr,j  )=ad_grad(Istr,j  )-MAX(Ce,0.0_r8)*adfac
                ad_grad(Istr,j+1)=ad_grad(Istr,j+1)-MIN(Ce,0.0_r8)*adfac
                ad_ubar(Istr  ,j,know)=ad_ubar(Istr  ,j,know)+cff*adfac
                ad_ubar(Istr+1,j,kout)=ad_ubar(Istr+1,j,kout)+Cx *adfac
                ad_ubar(Istr  ,j,kout)=0.0_r8
              END IF
            END DO
          END IF
!
!  Western edge, Flather boundary condition.
!
        ELSE IF (ad_LBC(iwest,isUbar,ng)%Flather) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                 &
     &                            zeta(Istr-1,j,know)+                  &
     &                            GRID(ng)%h(Istr  ,j)+                 &
     &                            zeta(Istr  ,j,know)))
              Cx=SQRT(g*cff)
# ifdef MASKING
!>            tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,kout)*                &
!>   &                             GRID(ng)%umask(Istr,j)
!>
              ad_ubar(Istr,j,kout)=ad_ubar(Istr,j,kout)*                &
     &                             GRID(ng)%umask(Istr,j)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(iwest,isUbar,ng)) THEN
!>              tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,kout)+              &
!>   &                               Cx*BOUNDARY(ng)%tl_zeta_west(j)
!>
                BOUNDARY(ng)%ad_zeta_west(j)=                           &
     &                                   BOUNDARY(ng)%ad_zeta_west(j)+  &
     &                                       Cx*ad_ubar(Istr,j,kout)
              END IF
# endif
!>            tl_ubar(Istr,j,kout)=tl_bry_val-                          &
!>   &                             tl_Cx*                               &
!>   &                             (0.5_r8*(zeta(Istr-1,j,know)+        &
!>   &                                      zeta(Istr  ,j,know))-       &
!>   &                              BOUNDARY(ng)%zeta_west(j))-         &
!>   &                             Cx*                                  &
!>   &                             (0.5_r8*(tl_zeta(Istr-1,j,know)+     &
!>   &                                      tl_zeta(Istr  ,j,know)))
!>
              adfac=Cx*0.5_r8*ad_ubar(Istr,j,kout)
              ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)-adfac
              ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)-adfac
              ad_Cx=ad_Cx-                                              &
     &              (0.5_r8*(zeta(Istr-1,j,know)+                       &
     &                       zeta(Istr  ,j,know))-                      &
     &               BOUNDARY(ng)%zeta_west(j))*ad_ubar(Istr,j,kout)
              ad_bry_val=ad_bry_val+ad_ubar(Istr,j,kout)
              ad_ubar(Istr,j,kout)=0.0_r8
!>            tl_Cx=0.5_r8*g*tl_cff/Cx
!>
              ad_cff=ad_cff+0.5_r8*g*ad_Cx/Cx
              ad_Cx=0.0_r8
!>            tl_cff=-cff*cff*(0.5_r8*(GRID(ng)%tl_h(Istr-1,j)+         &
!>   &                                 tl_zeta(Istr-1,j,know)+          &
!>   &                                 GRID(ng)%tl_h(Istr  ,j)+         &
!>   &                                 tl_zeta(Istr  ,j,know)))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              GRID(ng)%ad_h(Istr-1,j)=GRID(ng)%ad_h(Istr-1,j)+adfac
              GRID(ng)%ad_h(Istr  ,j)=GRID(ng)%ad_h(Istr  ,j)+adfac
              ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+adfac
              ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+adfac
              ad_cff=0.0_r8

# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (ad_LBC(iwest,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(Istr,j,know)-                          &
     &                      BOUNDARY(ng)%zeta_west(j))*                 &
     &                  0.5_r8*GRID(ng)%pm(Istr,j)
              ELSE
                bry_pgr=-g*(zeta(Istr  ,j,know)-                        &
     &                      zeta(Istr-1,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Istr-1,j)+                  &
     &                          GRID(ng)%pm(Istr,j))
              END IF
#  ifdef UV_COR
              bry_cor=0.125_r8*(vbar(Istr-1,j  ,know)+                  &
     &                          vbar(Istr-1,j+1,know)+                  &
     &                          vbar(Istr  ,j  ,know)+                  &
     &                          vbar(Istr  ,j+1,know))*                 &
     &                         (GRID(ng)%f(Istr-1,j)+                   &
     &                          GRID(ng)%f(Istr,j))
#  else
              bry_cor=0.0_r8
#  endif
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
!>            tl_bry_val=tl_ubar(Istr+1,j,know)+                        &
!>   &                   tl_cff2*(bry_pgr+                              &
!>   &                            bry_cor+                              &
!>   &                            bry_str)+                             &
!>   &                   cff2*(tl_bry_pgr+                              &
!>   &                         tl_bry_cor+                              &
!>   &                         tl_bry_str)
!>
              adfac=cff2*ad_bry_val
              ad_bry_pgr=ad_bry_pgr+adfac
              ad_bry_cor=ad_bry_cor+adfac
              ad_bry_str=ad_bry_str+adfac
              ad_cff2=ad_cff2+(bry_pgr+                                 &
     &                         bry_cor+                                 &
     &                         bry_str)*ad_bry_val
              ad_ubar(Istr+1,j,know)=ad_ubar(Istr+1,j,know)+ad_bry_val
              ad_bry_val=0.0_r8
!>            tl_cff2=GRID(ng)%om_u(Istr,j)*tl_Cx
!>
              ad_Cx=ad_Cx+GRID(ng)%om_u(Istr,j)*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Cx=-Cx*Cx*Cx*0.25_r8*g*(GRID(ng)%tl_h(Istr-1,j)+       &
!>   &                                   tl_zeta(Istr-1,j,know)+        &
!>   &                                   GRID(ng)%tl_h(Istr  ,j)+       &
!>   &                                   tl_zeta(Istr  ,j,know))
!>
              adfac=-Cx*Cx*Cx*0.25_r8*g*ad_Cx
              ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+adfac
              ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+adfac
              GRID(ng)%ad_h(Istr-1,j)=GRID(ng)%ad_h(Istr-1,j)+adfac
              GRID(ng)%ad_h(Istr  ,j)=GRID(ng)%ad_h(Istr  ,j)+adfac
              ad_Cx=0.0_r8
!>            tl_bry_str=tl_cff1*(FORCES(ng)%sustr(Istr,j)-             &
!>   &                            FORCES(ng)%bustr(Istr,j))+            &
!>   &                   cff1*(FORCES(ng)%tl_sustr(Istr,j)-             &
!>   &                         FORCES(ng)%tl_bustr(Istr,j))
!>
              adfac=cff1*ad_bry_str
              FORCES(ng)%ad_sustr(Istr,j)=FORCES(ng)%ad_sustr(Istr,j)+  &
     &                                    adfac
              FORCES(ng)%ad_bustr(Istr,j)=FORCES(ng)%ad_bustr(Istr,j)-  &
     &                                    adfac
              ad_cff1=ad_cff1+(FORCES(ng)%sustr(Istr,j)-                &
     &                         FORCES(ng)%bustr(Istr,j))*ad_bry_str
              ad_bry_str=0.0_r8
!>            tl_cff1=-cff1*cff1*0.5_r8*(GRID(ng)%tl_h(Istr-1,j)+       &
!>   &                                   tl_zeta(Istr-1,j,know)+        &
!>   &                                   GRID(ng)%tl_h(Istr  ,j)+       &
!>   &                                   tl_zeta(Istr  ,j,know))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+adfac
              ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+adfac
              GRID(ng)%ad_h(Istr-1,j)=GRID(ng)%ad_h(Istr-1,j)+adfac
              GRID(ng)%ad_h(Istr  ,j)=GRID(ng)%ad_h(Istr  ,j)+adfac
              ad_cff1=0.0_r8
#  ifdef UV_COR
!>            tl_bry_cor=0.125_r8*(tl_vbar(Istr-1,j  ,know)+            &
!>   &                             tl_vbar(Istr-1,j+1,know)+            &
!>   &                             tl_vbar(Istr  ,j  ,know)+            &
!>   &                             tl_vbar(Istr  ,j+1,know))*           &
!>   &                            (GRID(ng)%f(Istr-1,j)+                &
!>   &                             GRID(ng)%f(Istr  ,j))
!>
              adfac=0.125_r8*(GRID(ng)%f(Istr-1,j)+                     &
     &                        GRID(ng)%f(Istr  ,j))*ad_bry_cor
              ad_vbar(Istr-1,j  ,know)=ad_vbar(Istr-1,j  ,know)+adfac
              ad_vbar(Istr  ,j  ,know)=ad_vbar(Istr  ,j  ,know)+adfac
              ad_vbar(Istr-1,j+1,know)=ad_vbar(Istr-1,j+1,know)+adfac
              ad_vbar(Istr  ,j+1,know)=ad_vbar(Istr  ,j+1,know)+adfac
              ad_bry_cor=0.0_r8
#  else
!>            tl_bry_cor=0.0_r8
!>
              ad_bry_cor=0.0_r8
#  endif
              IF (ad_LBC(iwest,isFsur,ng)%acquire) THEN
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(iwest,isUbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr+                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_west(j)*            &
!>   &                       0.5_r8*GRID(ng)%pm(Istr,j)
!>
                  BOUNDARY(ng)%ad_zeta_west(j)=BOUNDARY(ng)%            &
     &                                               ad_zeta_west(j)+   &
     &                                         g*0.5_r8*                &
     &                                         GRID(ng)%pm(Istr,j)*     &
     &                                         ad_bry_pgr
                END IF
#  endif
!>              tl_bry_pgr=-g*tl_zeta(Istr,j,know)*                     &
!>   &                     0.5_r8*GRID(ng)%pm(Istr,j)
!>
                ad_zeta(Istr,j,know)=ad_zeta(Istr,j,know)-              &
     &                               g*0.5_r8*GRID(ng)%pm(Istr,j)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(Istr  ,j,know)-                  &
!>   &                         tl_zeta(Istr-1,j,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pm(Istr-1,j)+               &
!>   &                             GRID(ng)%pm(Istr  ,j))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pm(Istr-1,j)+                 &
     &                           GRID(ng)%pm(Istr  ,j))*ad_bry_pgr
                ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)-adfac
                ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
# else
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(iwest,isUbar,ng)) THEN
!>              tl_bry_val=BOUNDARY(ng)%tl_ubar_west(j)
!>
                BOUNDARY(ng)%ad_ubar_west(j)=                           &
     &                                   BOUNDARY(ng)%ad_ubar_west(j)+  &
     &                                       ad_bry_val
                ad_bry_val=0.0_r8
              ELSE
!>              tl_bry_val=0.0_r8
!>
                ad_bry_val=0.0_r8
              END IF
#  else
!>            tl_bry_val=0.0_r8
!>
              ad_bry_val=0.0_r8
#  endif
# endif
            END IF
          END DO
!
!  Western edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (ad_LBC(iwest,isUbar,ng)%Shchepetkin) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=0.5_r8*(GRID(ng)%h(Iend  ,j)+                         &
     &                    GRID(ng)%h(Iend+1,j))
              cff1=SQRT(g/cff)
              Cx=dt2d*cff1*cff*0.5_r8*(GRID(ng)%pm(Iend  ,j)+           &
     &                                 GRID(ng)%pm(Iend+1,j))
              Zx=(0.5_r8+Cx)*zeta(Iend  ,j,know)+                       &
     &           (0.5_r8-Cx)*zeta(Iend+1,j,know)
              IF (Cx.gt.Co) THEN
                cff2=(1.0_r8-Co/Cx)**2
                cff3=zeta(Iend,j,kout)+                                 &
     &               Cx*zeta(Iend+1,j,know)-                            &
     &               (1.0_r8+Cx)*zeta(Iend,j,know)
                Zx=Zx+cff2*cff3
              END IF
# ifdef MASKING
!>            tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,kout)*                &
!>   &                             GRID(ng)%umask(Istr,j)
!>
              ad_ubar(Istr,j,kout)=ad_ubar(Istr,j,kout)*                &
     &                             GRID(ng)%umask(Istr,j)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(iwest,isUbar,ng)) THEN
!>              tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,kout)+              &
!>   &                               Cx*BOUNDARY(ng)%tl_zeta_west(j)
!>
                BOUNDARY(ng)%ad_zeta_west(j)=                           &
     &                                   BOUNDARY(ng)%ad_zeta_west(j)+  &
     &                                       Cx*ad_ubar(Istr,j,kout)
              END IF
# endif
!>            tl_ubar(Iend+1,j,kout)=0.5_r8*                            &
!>   &                               ((1.0_r8-Cx)*                      &
!>   &                                tl_ubar(Iend+1,j,know)+           &
!>   &                                tl_Cx*(ubar(Iend  ,j,know)-       &
!>   &                                       ubar(Iend+1,j,know))+      &
!>   &                                Cx*tl_ubar(Iend,j,know)+          &
!>   &                                tl_bry_val+                       &
!>   &                                tl_cff1*                          &
!>   &                                (Zx-BOUNDARY(ng)%zeta_east(j))-   &
!>   &                                cff1*tl_Zx)
!>
              adfac=0.5_r8*ad_ubar(Iend+1,j,kout)
              ad_ubar(Iend+1,j,know)=ad_ubar(Iend+1,j,know)+            &
     &                               (1.0_r8-Cx)*adfac
              ad_ubar(Iend  ,j,know)=ad_ubar(Iend  ,j,know)+            &
     &                               Cx*adfac
              ad_Cx=ad_Cx+                                              &
     &              (ubar(Iend  ,j,know)-                               &
     &               ubar(Iend+1,j,know))*adfac
              ad_bry_val=ad_bry_val+adfac
              ad_cff1=ad_cff1+                                          &
     &                (Zx-BOUNDARY(ng)%zeta_east(j))*adfac
              ad_Zx=ad_Zx-cff1*adfac
              ad_ubar(Iend+1,j,kout)=0.0_r8
              IF (Cx.gt.Co) THEN
!>              tl_Zx=tl_Zx+cff2*tl_cff3+                               &
!>   &                      tl_cff2*cff3
!>
                ad_cff2=ad_cff2+cff3*ad_Zx
                ad_cff3=ad_cff3+cff2*ad_Zx
!>              tl_cff3=tl_zeta(Iend,j,kout)+                           &
!>   &                  Cx*tl_zeta(Iend+1,j,know)+                      &
!>   &                  tl_Cx*(zeta(Iend  ,j,know)+                     &
!>   &                         zeta(Iend+1,j,know))-                    &
!>   &                  (1.0_r8+Cx)*tl_zeta(Iend,j,know)
!>
                ad_zeta(Iend  ,j,kout)=ad_zeta(Iend,j,kout)+            &
     &                                 ad_cff3
                ad_zeta(Iend  ,j,know)=ad_zeta(Iend,j,know)-            &
     &                                 (1.0_r8+Cx)*ad_cff3
                ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+          &
     &                                 Cx*ad_cff3
                ad_Cx=ad_Cx+                                            &
     &                (zeta(Iend  ,j,know)+                             &
     &                 zeta(Iend+1,j,know))*ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=2.0_r8*cff2*Co*tl_Cx/(Cx*Cx)
!>
                ad_Cx=ad_Cx+                                            &
     &                2.0_r8*cff2*Co*ad_cff2/(Cx*Cx)
                ad_cff2=0.0_r8
              END IF
!>            tl_Zx=(0.5_r8+Cx)*tl_zeta(Iend  ,j,know)+                 &
!>   &              (0.5_r8-Cx)*tl_zeta(Iend+1,j,know)+                 &
!>   &              tl_Cx*(zeta(Iend  ,j,know)-                         &
!>   &                     zeta(Iend+1,j,know))
!>
              ad_zeta(Iend  ,j,know)=ad_zeta(Iend  ,j,know)+            &
     &                               (0.5_r8+Cx)*ad_Zx
              ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+            &
     &                               (0.5_r8-Cx)*ad_Zx
              ad_Cx=ad_Cx+                                              &
     &              (zeta(Iend  ,j,know)-                               &
     &               zeta(Iend+1,j,know))*ad_Zx
              ad_Zx=0.0_r8
!>            tl_Cx=dt2d*0.5_r8*(GRID(ng)%pm(Iend  ,j)+                 &
!>   &                           GRID(ng)%pm(Iend+1,j))*                &
!>   &              (cff1*tl_cff+                                       &
!>   &               tl_cff1*cff)
!>
              adfac=dt2d*0.5_r8*(GRID(ng)%pm(Iend  ,j)+                 &
     &                           GRID(ng)%pm(Iend+1,j))*ad_Cx
              ad_cff=ad_cff+cff1*adfac
              ad_cff1=ad_cff1+cff*adfac
              ad_Cx=0.0_r8
!>            tl_cff1=-0.5_r8*cff1*tl_cff/cff
!>
              ad_cff=ad_cff-                                            &
     &               0.5_r8*cff1*ad_cff/cff
              ad_cff1=0.0_r8
!>            tl_cff=0.5_r8*(GRID(ng)%tl_h(Iend  ,j)+                   &
!>   &                       GRID(ng)%tl_h(Iend+1,j))
!>
              adfac=0.5_r8*ad_cff
              GRID(ng)%ad_h(Iend  ,j)=GRID(ng)%ad_h(Iend  ,j)+adfac
              GRID(ng)%ad_h(Iend+1,j)=GRID(ng)%ad_h(Iend+1,j)+adfac
              ad_cff=0.0_r8

# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (ad_LBC(iwest,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(Istr,j,know)-                          &
     &                      BOUNDARY(ng)%zeta_west(j))*                 &
     &                  0.5_r8*GRID(ng)%pm(Istr,j)
              ELSE
                bry_pgr=-g*(zeta(Istr  ,j,know)-                        &
     &                      zeta(Istr-1,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Istr-1,j)+                  &
     &                          GRID(ng)%pm(Istr,j))
              END IF
#  ifdef UV_COR
              bry_cor=0.125_r8*(vbar(Istr-1,j  ,know)+                  &
     &                          vbar(Istr-1,j+1,know)+                  &
     &                          vbar(Istr  ,j  ,know)+                  &
     &                          vbar(Istr  ,j+1,know))*                 &
     &                         (GRID(ng)%f(Istr-1,j)+                   &
     &                          GRID(ng)%f(Istr,j))
#  else
              bry_cor=0.0_r8
#  endif
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
!>            tl_bry_val=tl_ubar(Istr+1,j,know)+                        &
!>   &                   tl_cff2*(bry_pgr+                              &
!>   &                            bry_cor+                              &
!>   &                            bry_str)+                             &
!>   &                   cff2*(tl_bry_pgr+                              &
!>   &                         tl_bry_cor+                              &
!>   &                         tl_bry_str)
!>
              adfac=cff2*ad_bry_val
              ad_bry_pgr=ad_bry_pgr+adfac
              ad_bry_cor=ad_bry_cor+adfac
              ad_bry_str=ad_bry_str+adfac
              ad_cff2=ad_cff2+(bry_pgr+                                 &
     &                         bry_cor+                                 &
     &                         bry_str)*ad_bry_val
              ad_ubar(Istr+1,j,know)=ad_ubar(Istr+1,j,know)+ad_bry_val
              ad_bry_val=0.0_r8
!>            tl_cff2=GRID(ng)%om_u(Istr,j)*tl_Cx
!>
              ad_Cx=ad_Cx+GRID(ng)%om_u(Istr,j)*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Cx=-Cx*Cx*Cx*0.25_r8*g*(GRID(ng)%tl_h(Istr-1,j)+       &
!>   &                                   tl_zeta(Istr-1,j,know)+        &
!>   &                                   GRID(ng)%tl_h(Istr  ,j)+       &
!>   &                                   tl_zeta(Istr  ,j,know))
!>
              adfac=-Cx*Cx*Cx*0.25_r8*g*ad_Cx
              ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+adfac
              ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+adfac
              GRID(ng)%ad_h(Istr-1,j)=GRID(ng)%ad_h(Istr-1,j)+adfac
              GRID(ng)%ad_h(Istr  ,j)=GRID(ng)%ad_h(Istr  ,j)+adfac
              ad_Cx=0.0_r8
!>            tl_bry_str=tl_cff1*(FORCES(ng)%sustr(Istr,j)-             &
!>   &                            FORCES(ng)%bustr(Istr,j))+            &
!>   &                   cff1*(FORCES(ng)%tl_sustr(Istr,j)-             &
!>   &                         FORCES(ng)%tl_bustr(Istr,j))
!>
              adfac=cff1*ad_bry_str
              FORCES(ng)%ad_sustr(Istr,j)=FORCES(ng)%ad_sustr(Istr,j)+  &
     &                                    adfac
              FORCES(ng)%ad_bustr(Istr,j)=FORCES(ng)%ad_bustr(Istr,j)-  &
     &                                    adfac
              ad_cff1=ad_cff1+(FORCES(ng)%sustr(Istr,j)-                &
     &                         FORCES(ng)%bustr(Istr,j))*ad_bry_str
              ad_bry_str=0.0_r8
!>            tl_cff1=-cff1*cff1*0.5_r8*(GRID(ng)%tl_h(Istr-1,j)+       &
!>   &                                   tl_zeta(Istr-1,j,know)+        &
!>   &                                   GRID(ng)%tl_h(Istr  ,j)+       &
!>   &                                   tl_zeta(Istr  ,j,know))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+adfac
              ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+adfac
              GRID(ng)%ad_h(Istr-1,j)=GRID(ng)%ad_h(Istr-1,j)+adfac
              GRID(ng)%ad_h(Istr  ,j)=GRID(ng)%ad_h(Istr  ,j)+adfac
              ad_cff1=0.0_r8
#  ifdef UV_COR
!>            tl_bry_cor=0.125_r8*(tl_vbar(Istr-1,j  ,know)+            &
!>   &                             tl_vbar(Istr-1,j+1,know)+            &
!>   &                             tl_vbar(Istr  ,j  ,know)+            &
!>   &                             tl_vbar(Istr  ,j+1,know))*           &
!>   &                            (GRID(ng)%f(Istr-1,j)+                &
!>   &                             GRID(ng)%f(Istr  ,j))
!>
              adfac=0.125_r8*(GRID(ng)%f(Istr-1,j)+                     &
     &                        GRID(ng)%f(Istr  ,j))*ad_bry_cor
              ad_vbar(Istr-1,j  ,know)=ad_vbar(Istr-1,j  ,know)+adfac
              ad_vbar(Istr  ,j  ,know)=ad_vbar(Istr  ,j  ,know)+adfac
              ad_vbar(Istr-1,j+1,know)=ad_vbar(Istr-1,j+1,know)+adfac
              ad_vbar(Istr  ,j+1,know)=ad_vbar(Istr  ,j+1,know)+adfac
              ad_bry_cor=0.0_r8
#  else
!>            tl_bry_cor=0.0_r8
!>
              ad_bry_cor=0.0_r8
#  endif
              IF (ad_LBC(iwest,isFsur,ng)%acquire) THEN
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(iwest,isUbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr+                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_west(j)*            &
!>   &                       0.5_r8*GRID(ng)%pm(Istr,j)
!>
                  BOUNDARY(ng)%ad_zeta_west(j)=BOUNDARY(ng)%            &
     &                                               ad_zeta_west(j)+   &
     &                                         g*0.5_r8*                &
     &                                         GRID(ng)%pm(Istr,j)*     &
     &                                         ad_bry_pgr
                END IF
#  endif
!>              tl_bry_pgr=-g*tl_zeta(Istr,j,know)*                     &
!>   &                     0.5_r8*GRID(ng)%pm(Istr,j)
!>
                ad_zeta(Istr,j,know)=ad_zeta(Istr,j,know)-              &
     &                               g*0.5_r8*GRID(ng)%pm(Istr,j)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(Istr  ,j,know)-                  &
!>   &                         tl_zeta(Istr-1,j,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pm(Istr-1,j)+               &
!>   &                             GRID(ng)%pm(Istr  ,j))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pm(Istr-1,j)+                 &
     &                           GRID(ng)%pm(Istr  ,j))*ad_bry_pgr
                ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)-adfac
                ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
# else
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(iwest,isUbar,ng)) THEN
!>              tl_bry_val=BOUNDARY(ng)%tl_ubar_west(j)
!>
                BOUNDARY(ng)%ad_ubar_west(j)=                           &
     &                                   BOUNDARY(ng)%ad_ubar_west(j)+  &
     &                                       ad_bry_val
                ad_bry_val=0.0_r8
              ELSE
!>              tl_bry_val=0.0_r8
!>
                ad_bry_val=0.0_r8
              END IF
#  else
!>            tl_bry_val=0.0_r8
!>
              ad_bry_val=0.0_r8
#  endif
# endif
            END IF
          END DO









!
!  Western edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(iwest,isUbar,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
# ifdef MASKING
!>            tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,kout)*                &
!>   &                             GRID(ng)%umask(Istr,j)
!>
              ad_ubar(Istr,j,kout)=ad_ubar(Istr,j,kout)*                &
     &                             GRID(ng)%umask(Istr,j)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(iwest,isUbar,ng)) THEN
!>              tl_ubar(Istr,j,kout)=BOUNDARY(ng)%tl_ubar_west(j)
!>
                BOUNDARY(ng)%ad_ubar_west(j)=                           &
     &                                   BOUNDARY(ng)%ad_ubar_west(j)+  &
     &                                       ad_ubar(Istr,j,kout)
                ad_ubar(Istr,j,kout)=0.0_r8
              ELSE
!>              tl_ubar(Istr,j,kout)=0.0_r8
!>
                ad_ubar(Istr,j,kout)=0.0_r8
              END IF
# else
!>            tl_ubar(Istr,j,kout)=0.0_r8
!>
              ad_ubar(Istr,j,kout)=0.0_r8
# endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(iwest,isUbar,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
# ifdef MASKING
!>            tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,kout)*                &
!>   &                             GRID(ng)%umask(Istr,j)
!>
              ad_ubar(Istr,j,kout)=ad_ubar(Istr,j,kout)*                &
     &                             GRID(ng)%umask(Istr,j)
# endif
!>            tl_ubar(Istr  ,j,kout)=tl_ubar(Istr+1,j,kout)
!>
              ad_ubar(Istr+1,j,kout)=ad_ubar(Istr+1,j,kout)+            &
     &                               ad_ubar(Istr,j,kout)
              ad_ubar(Istr  ,j,kout)=0.0_r8
            END IF
          END DO
!
!  Western edge, reduced-physics boundary condition.
!
        ELSE IF (ad_LBC(iwest,isUbar,ng)%reduced) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                 &
     &                            zeta(Iend  ,j,know)+                  &
     &                            GRID(ng)%h(Iend+1,j)+                 &
     &                            zeta(Iend+1,j,know)))
# ifdef MASKING
!>            tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,kout)*                &
!>   &                             GRID(ng)%umask(Istr,j)
!>
              ad_ubar(Istr,j,kout)=ad_ubar(Istr,j,kout)*                &
     &                             GRID(ng)%umask(Istr,j)
# endif
!>            tl_ubar(Istr,j,kout)=tl_ubar(Istr,j,know)+                &
!>   &                             dt2d*(tl_bry_pgr+                    &
!>   &                                   tl_bry_cor+                    &
!>   &                                   tl_bry_str)
!>
              adfac=dt2d*ad_ubar(Istr,j,kout)
              ad_bry_pgr=ad_bry_pgr+adfac
              ad_bry_cor=ad_bry_cor+adfac
              ad_bry_str=ad_bry_str+adfac
              ad_ubar(Istr,j,know)=ad_ubar(Istr,j,know)+                &
     &                             ad_ubar(Istr,j,kout)
              ad_ubar(Istr,j,kout)=0.0_r8
!>            tl_bry_str=tl_cff*(FORCES(ng)%sustr(Istr,j)-              &
!>   &                           FORCES(ng)%bustr(Istr,j))+             &
!>   &                   cff*(FORCES(ng)%tl_sustr(Istr,j)-              &
!>   &                        FORCES(ng)%tl_bustr(Istr,j))
!>
              adfac=cff*ad_bry_str
              FORCES(ng)%ad_sustr(Istr,j)=FORCES(ng)%ad_sustr(Istr,j)+  &
     &                                    adfac
              FORCES(ng)%ad_bustr(Istr,j)=FORCES(ng)%ad_bustr(Istr,j)-  &
     &                                    adfac
              ad_cff=ad_cff+(FORCES(ng)%sustr(Istr,j)-                  &
     &                       FORCES(ng)%bustr(Istr,j))*ad_bry_str
              ad_bry_str=0.0_r8
!>            tl_cff=-cff*cff*0.5_r8*(GRID(ng)%tl_h(Istr-1,j)+          &
!>   &                                tl_zeta(Istr-1,j,know)+           &
!>   &                                GRID(ng)%tl_h(Istr  ,j)+          &
!>   &                                tl_zeta(Istr  ,j,know))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+adfac
              ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+adfac
              GRID(ng)%ad_h(Istr-1,j)=GRID(ng)%ad_h(Istr-1,j)+adfac
              GRID(ng)%ad_h(Istr  ,j)=GRID(ng)%ad_h(Istr  ,j)+adfac
              ad_cff=0.0_r8
# ifdef UV_COR
!>            tl_bry_cor=0.125_r8*(tl_vbar(Istr-1,j  ,know)+            &
!>   &                             tl_vbar(Istr-1,j+1,know)+            &
!>   &                             tl_vbar(Istr  ,j  ,know)+            &
!>   &                             tl_vbar(Istr  ,j+1,know))*           &
!>   &                            (GRID(ng)%f(Istr-1,j)+                &
!>   &                             GRID(ng)%f(Istr  ,j))
!>
              adfac=0.125_r8*(GRID(ng)%f(Istr-1,j)+                     &
     &                        GRID(ng)%f(Istr  ,j))*ad_bry_cor
              ad_vbar(Istr-1,j  ,know)=ad_vbar(Istr-1,j  ,know)+adfac
              ad_vbar(Istr  ,j  ,know)=ad_vbar(Istr  ,j  ,know)+adfac
              ad_vbar(Istr-1,j+1,know)=ad_vbar(Istr-1,j+1,know)+adfac
              ad_vbar(Istr  ,j+1,know)=ad_vbar(Istr  ,j+1,know)+adfac
              ad_bry_cor=0.0_r8
# else
!>            tl_bry_cor=0.0_r8
!>
              ad_bry_cor=0.0_r8
# endif
              IF (ad_LBC(iwest,isFsur,ng)%acquire) THEN
# ifdef ADJUST_BOUNDARY
                IF (Lobc(iwest,isUbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr+                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_west(j)*            &
!>   &                       0.5_r8*GRID(ng)%pm(Istr,j)
!>
                  BOUNDARY(ng)%ad_zeta_west(j)=BOUNDARY(ng)%            &
     &                                               ad_zeta_west(j)+   &
     &                                         g*0.5_r8*                &
     &                                         GRID(ng)%pm(Istr,j)*     &
     &                                         ad_bry_pgr
                END IF
# endif
!>              tl_bry_pgr=-g*(zeta(Istr,j,know)*                       &
!>   &                     0.5_r8*GRID(ng)%pm(Istr,j)
!>
                ad_zeta(Istr,j,know)=ad_zeta(Istr,j,know)-              &
     &                               g*0.5_r8*GRID(ng)%pm(Istr,j)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(Istr  ,j,know)-                  &
!>   &                         tl_zeta(Istr-1,j,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pm(Istr-1,j)+               &
!>   &                             GRID(ng)%pm(Istr  ,j))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pm(Istr-1,j)+                 &
     &                           GRID(ng)%pm(Istr  ,j))*ad_bry_pgr
                ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)-adfac
                ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (ad_LBC(iwest,isUbar,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
!>            tl_ubar(Istr,j,kout)=0.0_r8
!>
              ad_ubar(Istr,j,kout)=0.0_r8
            END IF
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_u2dbc_tile
#endif
      END MODULE ad_u2dbc_mod
