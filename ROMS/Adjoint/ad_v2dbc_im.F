#include "cppdefs.h"
      MODULE ad_v2dbc_mod
#ifdef ADJOINT
!
!svn $Id: ad_v2dbc_im.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets adjoint lateral boundary conditions for        !
!  vertically integrated  V-velocity.  It updates the specified        !
!  "kout" index.                                                       !
!                                                                      !
!  BASIC STATE variables needed: zeta                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_v2dbc, ad_v2dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_v2dbc (ng, tile, kout)
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
      CALL ad_v2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs(ng), kstp(ng), kout,                     &
     &                    OCEAN(ng) % ubar,                             &
     &                    OCEAN(ng) % vbar,                             &
     &                    OCEAN(ng) % zeta,                             &
     &                    OCEAN(ng) % ad_ubar,                          &
     &                    OCEAN(ng) % ad_vbar,                          &
     &                    OCEAN(ng) % ad_zeta)

      RETURN
      END SUBROUTINE ad_v2dbc

!
!***********************************************************************
      SUBROUTINE ad_v2dbc_tile (ng, tile,                               &
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
      integer :: Jmin, Jmax
      integer :: i, j, know

      real(r8) :: Ce, Cx, Ze
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8) :: cff, cff1, cff2, cff3, dt2d
      real(r8) :: obc_in, obc_out, tau

      real(r8) :: ad_Ce, ad_Cx
      real(r8) :: ad_bry_pgr, ad_bry_cor, ad_bry_str, ad_bry_val, ad_Ze
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
      ad_Ze=0.0_r8
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
!>          cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,Jend+1))-1.0_r8)
!>          cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,Jend+1,kout))*         &
!>   &                  GRID(ng)%vmask_wet(Iend+1,Jend+1)
!>          cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,Jend+1)*cff1+          &
!>   &          cff2*(1.0_r8-cff1)
!>          vbar(Iend+1,Jend+1,kout)=vbar(Iend+1,Jend+1,kout)*cff
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
      END IF

      IF (.not.NSperiodic(ng)) THEN
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
      END IF

      IF (.not.EWperiodic(ng)) THEN
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
      ENDIF
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
!>          tl_vbar(Iend+1,Jend+1,kout)=0.5_r8*                         &
!>   &                                  (tl_vbar(Iend+1,Jend  ,kout)+   &
!>   &                                   tl_vbar(Iend  ,Jend+1,kout))
!>
            adfac=0.5_r8*ad_vbar(Iend+1,Jend+1,kout)
            ad_vbar(Iend+1,Jend  ,kout)=ad_vbar(Iend+1,Jend  ,kout)+    &
     &                                  adfac
            ad_vbar(Iend  ,Jend+1,kout)=ad_vbar(Iend  ,Jend+1,kout)+    &
     &                                  adfac
            ad_vbar(Iend+1,Jend+1,kout)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
!>          tl_vbar(Istr-1,Jend+1,kout)=0.5_r8*                         &
!>   &                                  (tl_vbar(Istr-1,Jend  ,kout)+   &
!>   &                                   tl_vbar(Istr  ,Jend+1,kout))
!>
            adfac=0.5_r8*ad_vbar(Istr-1,Jend+1,kout)
            ad_vbar(Istr-1,Jend  ,kout)=ad_vbar(Istr-1,Jend  ,kout)+    &
     &                                  adfac
            ad_vbar(Istr  ,Jend+1,kout)=ad_vbar(Istr  ,Jend+1,kout)+    &
     &                                  adfac
            ad_vbar(Istr-1,Jend+1,kout)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
!>          tl_vbar(Iend+1,Jstr,kout)=0.5_r8*                           &
!>   &                                (tl_vbar(Iend  ,Jstr  ,kout)+     &
!>   &                                 tl_vbar(Iend+1,Jstr+1,kout))
!>
            adfac=0.5_r8*ad_vbar(Iend+1,Jstr,kout)
            ad_vbar(Iend  ,Jstr  ,kout)=ad_vbar(Iend  ,Jstr  ,kout)+    &
     &                                  adfac
            ad_vbar(Iend+1,Jstr+1,kout)=ad_vbar(Iend+1,Jstr+1,kout)+    &
     &                                  adfac
            ad_vbar(Iend+1,Jstr  ,kout)=0.0_r8
          END IF
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
!>          tl_vbar(Istr-1,Jstr,kout)=0.5_r8*                           &
!>   &                                (tl_vbar(Istr  ,Jstr  ,kout)+     &
!>   &                                 tl_vbar(Istr-1,Jstr+1,kout))
!>
            adfac=0.5_r8*ad_vbar(Istr-1,Jstr,kout)
            ad_vbar(Istr  ,Jstr  ,kout)=ad_vbar(Istr  ,Jstr  ,kout)+    &
     &                                  adfac
            ad_vbar(Istr-1,Jstr+1,kout)=ad_vbar(Istr-1,Jstr+1,kout)+    &
     &                                  adfac
            ad_vbar(Istr-1,Jstr  ,kout)=0.0_r8
          END IF
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
        IF (ad_LBC(ieast,isVbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (ad_LBC(ieast,isVbar,ng)%nudging) THEN
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
# ifdef MASKING
!>              tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*          &
!>   &                                 GRID(ng)%vmask(Iend+1,j)
!>
                ad_vbar(Iend+1,j,kout)=ad_vbar(Iend+1,j,kout)*          &
     &                                 GRID(ng)%vmask(Iend+1,j)
# endif
                IF (ad_LBC(ieast,isVbar,ng)%nudging) THEN
!>                tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)-        &
!>   &                                   tau*tl_vbar(Iend+1,j,know)
!>
                  ad_vbar(Iend+1 ,j,know)=ad_vbar(Iend+1 ,j,know)-      &
     &                                    tau*ad_vbar(Iend+1,j,kout)
                END IF
!>              tl_vbar(Iend+1,j,kout)=(cff*tl_vbar(Iend+1,j,know)+     &
!>   &                                  Cx *tl_vbar(Iend  ,j,kout)-     &
!>   &                                  MAX(Ce,0.0_r8)*                 &
!>   &                                     tl_grad(Iend+1,j-1)-         &
!>   &                                  MIN(Ce,0.0_r8)*                 &
!>   &                                     tl_grad(Iend+1,j  ))/        &
!>   &                                 (cff+Cx)
!>
                adfac=ad_vbar(Iend+1,j,kout)/(cff+Cx)
                ad_grad(Iend+1,j-1)=ad_grad(Iend+1,j-1)-                &
     &                              MAX(Ce,0.0_r8)*adfac
                ad_grad(Iend+1,j  )=ad_grad(Iend+1,j  )-                &
     &                              MIN(Ce,0.0_r8)*adfac
                ad_vbar(Iend  ,j,kout)=ad_vbar(Iend  ,j,kout)+Cx* adfac
                ad_vbar(Iend+1,j,know)=ad_vbar(Iend+1,j,know)+cff*adfac
                ad_vbar(Iend+1,j,kout)=0.0_r8
              END IF
            END DO
          END IF
!
!  Eastern edge, Chapman boundary condition.
!
        ELSE IF (ad_LBC(ieast,isVbar,ng)%Flather.or.                    &
     &           ad_LBC(ieast,isVbar,ng)%reduced.or.                    &
     &           ad_LBC(ieast,isVbar,ng)%Shchepetkin) THEN
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
# ifdef MASKING
!>            tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*            &
!>   &                               GRID(ng)%vmask(Iend+1,j)
!>
              ad_vbar(Iend+1,j,kout)=ad_vbar(Iend+1,j,kout)*            &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
!>            tl_vbar(Iend+1,j,kout)=tl_cff2*(vbar(Iend+1,j,know)+      &
!>   &                                        Cx*vbar(Iend,j,kout))+    &
!>   &                               cff2*(tl_vbar(Iend+1,j,know)+      &
!>   &                                     tl_Cx*vbar(Iend,j,kout)+     &
!>   &                                     Cx*tl_vbar(Iend,j,kout))
!>
              adfac=cff2*ad_vbar(Iend+1,j,kout)
              ad_vbar(Iend  ,j,kout)=ad_vbar(Iend  ,j,kout)+Cx*adfac
              ad_vbar(Iend+1,j,know)=ad_vbar(Iend+1,j,know)+adfac
              ad_Cx=ad_Cx+vbar(Iend,j,kout)*adfac
              ad_cff2=ad_cff2+                                          &
     &                (vbar(Iend+1,j,know)+                             &
     &                 Cx*vbar(Iend,j,kout))*ad_vbar(Iend+1,j,kout)
              ad_vbar(Iend+1,j,kout)=0.0_r8
!>            tl_cff2=-cff2*cff2*tl_Cx
!>
              ad_Cx=ad_Cx-cff2*cff2*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Cx=cff*tl_cff1
!>
              ad_cff1=ad_cff1+cff*ad_Cx
              ad_Cx=0.0_r8
!>            tl_cff1=0.25_r8*g*(GRID(ng)%tl_h(Iend,j-1)+               &
!>   &                           tl_zeta(Iend,j-1,know)+                &
!>   &                           GRID(ng)%tl_h(Iend,j  )+               &
!>   &                           tl_zeta(Iend,j  ,know))/cff1
!>
              adfac=0.25_r8*g*ad_cff1/cff1
              GRID(ng)%ad_h(Iend,j-1)=GRID(ng)%ad_h(Iend,j-1)+adfac
              GRID(ng)%ad_h(Iend,j  )=GRID(ng)%ad_h(Iend,j  )+adfac
              ad_zeta(Iend,j-1,know)=ad_zeta(Iend,j-1,know)+adfac
              ad_zeta(Iend,j  ,know)=ad_zeta(Iend,j  ,know)+adfac
              ad_cff1=0.0_r8
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(ieast,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
# ifdef MASKING
!>            tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*            &
!>   &                               GRID(ng)%vmask(Iend+1,j)
!>
              ad_vbar(Iend+1,j,kout)=ad_vbar(Iend+1,j,kout)*            &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(ieast,isVbar,ng)) THEN
!>              tl_vbar(Iend+1,j,kout)=BOUNDARY(ng)%tl_vbar_east(j)
!>
                BOUNDARY(ng)%ad_vbar_east(j)=                           &
     &                                   BOUNDARY(ng)%ad_vbar_east(j)+  &
     &                                       ad_vbar(Iend+1,j,kout)
                ad_vbar(Iend+1,j,kout)=0.0_r8
              ELSE
!>              tl_vbar(Iend+1,j,kout)=0.0_r8
!>
                ad_vbar(Iend+1,j,kout)=0.0_r8
              END IF
# else
!>            tl_vbar(Iend+1,j,kout)=0.0_r8
!>
              ad_vbar(Iend+1,j,kout)=0.0_r8
# endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(ieast,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
# ifdef MASKING
!>            tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*            &
!>   &                               GRID(ng)%vmask(Iend+1,j)
!>
              ad_vbar(Iend+1,j,kout)=ad_vbar(Iend+1,j,kout)*            &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
!>            tl_vbar(Iend+1,j,kout)=tl_vbar(Iend,j,kout)
!>
              ad_vbar(Iend  ,j,kout)=ad_vbar(Iend,j,kout)+              &
     &                               ad_vbar(Iend+1,j,kout)
              ad_vbar(Iend+1,j,kout)=0.0_r8
            END IF
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (ad_LBC(ieast,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%east(j)) THEN
# ifdef MASKING
!>            tl_vbar(Iend+1,j,kout)=tl_vbar(Iend+1,j,kout)*            &
!>   &                               GRID(ng)%vmask(Iend+1,j)
!>
              ad_vbar(Iend+1,j,kout)=ad_vbar(Iend+1,j,kout)*            &
     &                               GRID(ng)%vmask(Iend+1,j)
# endif
!>            tl_vbar(Iend+1,j,kout)=gamma2(ng)*tl_vbar(Iend,j,kout)
!>
              ad_vbar(Iend  ,j,kout)=ad_vbar(Iend,j,kout)+              &
     &                               gamma2(ng)*ad_vbar(Iend+1,j,kout)
              ad_vbar(Iend+1,j,kout)=0.0_r8
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
        IF (ad_LBC(iwest,isVbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (ad_LBC(iwest,isVbar,ng)%nudging) THEN
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
# ifdef MASKING
!>              tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*          &
!>   &                                 GRID(ng)%vmask(Istr-1,j)
!>
                ad_vbar(Istr-1,j,kout)=ad_vbar(Istr-1,j,kout)*          &
     &                                  GRID(ng)%vmask(Istr-1,j)
# endif
                IF (ad_LBC(iwest,isVbar,ng)%nudging) THEN
!>                tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)-        &
!>   &                                   tau*tl_vbar(1,j,know)
!>
                  ad_vbar(Istr,j,know)=ad_vbar(Istr,j,know)-            &
     &                               tau*ad_vbar(Istr-1,j,kout)
                END IF
!>              tl_vbar(Istr-1,j,kout)=(cff*tl_vbar(Istr-1,j,know)+     &
!>   &                                  Cx *tl_vbar(1,j,kout)-          &
!>   &                                  MAX(Ce,0.0_r8)*                 &
!>   &                                     tl_grad(Istr-1,j-1)-         &
!>   &                                  MIN(Ce,0.0_r8)*                 &
!>   &                                     tl_grad(Istr-1,j  ))/        &
!>   &                                 (cff+Cx)
!>
                adfac=ad_vbar(Istr-1,j,kout)/(cff+Cx)
                ad_grad(Istr-1,j-1)=ad_grad(Istr-1,j-1)-                &
     &                              MAX(Ce,0.0_r8)*adfac
                ad_grad(Istr-1,j  )=ad_grad(Istr-1,j  )-                &
     &                              MIN(Ce,0.0_r8)*adfac
                ad_vbar(Istr-1,j,know)=ad_vbar(Istr-1,j,know)+cff*adfac
                ad_vbar(Istr  ,j,kout)=ad_vbar(Istr  ,j,kout)+Cx *adfac
                ad_vbar(Istr-1,j,kout)=0.0_r8
              END IF
            END DO
          END IF
!
!  Western edge, Chapman boundary condition.
!
        ELSE IF (ad_LBC(iwest,isVbar,ng)%Flather.or.                    &
     &           ad_LBC(iwest,isVbar,ng)%reduced.or.                    &
     &           ad_LBC(iwest,isVbar,ng)%Shchepetkin) THEN
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
# ifdef MASKING
!>            tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*            &
!>   &                               GRID(ng)%vmask(Istr-1,j)
!>
              ad_vbar(Istr-1,j,kout)=ad_vbar(Istr-1,j,kout)*            &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
!>            tl_vbar(Istr-1,j,kout)=tl_cff2*(vbar(Istr-1,j,know)+      &
!>   &                                        Cx*vbar(Istr,j,kout))+    &
!>   &                               cff2*(tl_vbar(Istr-1,j,know)+      &
!>   &                                     tl_Cx*vbar(Istr,j,kout)+     &
!>   &                                     Cx*tl_vbar(Istr,j,kout))
!>
              adfac=cff2*ad_vbar(Istr-1,j,kout)
              ad_vbar(Istr-1,j,know)=ad_vbar(Istr-1,j,know)+adfac
              ad_vbar(Istr  ,j,kout)=ad_vbar(Istr  ,j,kout)+Cx*adfac
              ad_Cx=ad_Cx+vbar(Istr,j,kout)*adfac
              ad_cff2=ad_cff2+                                          &
     &                (vbar(Istr-1,j,know)+                             &
     &                 Cx*vbar(Istr,j,kout))*ad_vbar(Istr-1,j,kout)
              ad_vbar(Istr-1,j,kout)=0.0_r8
!>            tl_cff2=-cff2*cff2*tl_Cx
!>
              ad_Cx=ad_Cx-cff2*cff2*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Cx=cff*tl_cff1
!>
              ad_cff1=ad_cff1+cff*ad_Cx
              ad_Cx=0.0_r8
!>            tl_cff1=0.25_r8*g*(GRID(ng)%tl_h(Istr,j-1)+               &
!>   &                           tl_zeta(Istr,j-1,know)+                &
!>   &                           GRID(ng)%tl_h(Istr,j  )+               &
!>   &                           tl_zeta(Istr,j  ,know))/cff1
!>
              adfac=0.25_r8*g*ad_cff1/cff1
              GRID(ng)%ad_h(Istr,j-1)=GRID(ng)%ad_h(Istr,j-1)+adfac
              GRID(ng)%ad_h(Istr,j  )=GRID(ng)%ad_h(Istr,j  )+adfac
              ad_zeta(Istr,j-1,know)=ad_zeta(Istr,j-1,know)+adfac
              ad_zeta(Istr,j  ,know)=ad_zeta(Istr,j  ,know)+adfac
              ad_cff1=0.0_r8
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(iwest,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
# ifdef MASKING
!>            tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*            &
!>   &                               GRID(ng)%vmask(Istr-1,j)
!>
              ad_vbar(Istr-1,j,kout)=ad_vbar(Istr-1,j,kout)*            &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(iwest,isVbar,ng)) THEN
!>              tl_vbar(Istr-1,j,kout)=BOUNDARY(ng)%tl_vbar_west(j)
!>
                BOUNDARY(ng)%ad_vbar_west(j)=                           &
     &                                   BOUNDARY(ng)%ad_vbar_west(j)+  &
     &                                       ad_vbar(Istr-1,j,kout)
                ad_vbar(Istr-1,j,kout)=0.0_r8
              ELSE
!>              tl_vbar(Istr-1,j,kout)=0.0_r8
!>
                ad_vbar(Istr-1,j,kout)=0.0_r8
              END IF
# else
!>            tl_vbar(Istr-1,j,kout)=0.0_r8
!>
              ad_vbar(Istr-1,j,kout)=0.0_r8
# endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(iwest,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
# ifdef MASKING
!>            tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*            &
!>   &                               GRID(ng)%vmask(Istr-1,j)
!>
              ad_vbar(Istr-1,j,kout)=ad_vbar(Istr-1,j,kout)*            &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
!>            tl_vbar(Istr-1,j,kout)=tl_vbar(Istr,j,kout)
!>
              ad_vbar(Istr  ,j,kout)=ad_vbar(Istr,j,kout)+              &
     &                               ad_vbar(Istr-1,j,kout)
              ad_vbar(Istr-1,j,kout)=0.0_r8
            END IF
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (ad_LBC(iwest,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%west(j)) THEN
# ifdef MASKING
!>            tl_vbar(Istr-1,j,kout)=tl_vbar(Istr-1,j,kout)*            &
!>   &                               GRID(ng)%vmask(Istr-1,j)
!>
              ad_vbar(Istr-1,j,kout)=ad_vbar(Istr-1,j,kout)*            &
     &                               GRID(ng)%vmask(Istr-1,j)
# endif
!>            tl_vbar(Istr-1,j,kout)=gamma2(ng)*tl_vbar(Istr,j,kout)
!>
              ad_vbar(Istr  ,j,kout)=ad_vbar(Istr,j,kout)+              &
     &                               gamma2(ng)*ad_vbar(Istr-1,j,kout)
              ad_vbar(Istr-1,j,kout)=0.0_r8
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
        IF (ad_LBC(inorth,isVbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (ad_LBC(inorth,isVbar,ng)%nudging) THEN
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
# ifdef MASKING
!>              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*          &
!>   &                                 GRID(ng)%vmask(i,Jend+1)
!>
                ad_vbar(i,Jend+1,kout)=ad_vbar(i,Jend+1,kout)*          &
     &                                 GRID(ng)%vmask(i,Jend+1)
# endif
                IF (ad_LBC(inorth,isVbar,ng)%nudging) THEN
!>                tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)-        &
!>   &                                   tau*tl_vbar(i,Jend+1,know)
!>
                  ad_vbar(i,Jend+1 ,know)=ad_vbar(i,Jend+1 ,know)-      &
     &                                    tau*ad_vbar(i,Jend+1,kout)
                END IF
!>              tl_vbar(i,Jend+1,kout)=(cff*tl_vbar(i,Jend+1,know)+     &
!>   &                                  Ce *tl_vbar(i,Jend  ,kout)-     &
!>   &                                  MAX(Cx,0.0_r8)*                 &
!>   &                                     tl_grad(i  ,Jend+1)-         &
!>   &                                  MIN(Cx,0.0_r8)*                 &
!>   &                                     tl_grad(i+1,Jend+1))/        &
!>   &                                 (cff+Ce)
!>
                adfac=ad_vbar(i,Jend+1,kout)/(cff+Ce)
                ad_grad(i  ,Jend+1)=ad_grad(i  ,Jend+1)-                &
     &                              MAX(Cx,0.0_r8)*adfac
                ad_grad(i+1,Jend+1)=ad_grad(i+1,Jend+1)-                &
     &                              MIN(Cx,0.0_r8)*adfac
                ad_vbar(i,Jend  ,kout)=ad_vbar(i,Jend  ,kout)+Ce*adfac
                ad_vbar(i,Jend+1,know)=ad_vbar(i,Jend+1,know)+cff*adfac
                ad_vbar(i,Jend+1,kout)=0.0_r8
              END IF
            END DO
          END IF
!
!  Northern edge, Flather boundary condition.
!
        ELSE IF (ad_LBC(inorth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                 &
     &                            zeta(i,Jend  ,know)+                  &
     &                            GRID(ng)%h(i,Jend+1)+                 &
     &                            zeta(i,Jend+1,know)))
              Ce=SQRT(g*cff)
# ifdef MASKING
!>            tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
!>   &                               GRID(ng)%vmask(i,Jend+1)
!>
              ad_vbar(i,Jend+1,kout)=ad_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
!>              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)-          &
!>   &                                 Ce*BOUNDARY(ng)%tl_zeta_north(i)
!>
                BOUNDARY(ng)%ad_zeta_north(i)=BOUNDARY(ng)%             &
     &                                              ad_zeta_north(i)-   &
     &                                        Ce*ad_vbar(i,Jend+1,kout)
              END IF
# endif
!>            tl_vbar(i,Jend+1,kout)=tl_bry_val+                        &
!>   &                               tl_Ce*                             &
!>   &                               (0.5_r8*(zeta(i,Jend  ,know)+      &
!>   &                                        zeta(i,Jend+1,know))-     &
!>   &                                BOUNDARY(ng)%zeta_north(i))+      &
!>   &                               Ce*                                &
!>   &                               (0.5_r8*(tl_zeta(i,Jend  ,know)+   &
!>   &                                        tl_zeta(i,Jend+1,know)))
!>
              adfac=Ce*0.5_r8*ad_vbar(i,Jend+1,kout)
              ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)+adfac
              ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
              ad_Ce=ad_Ce+                                              &
     &              (0.5_r8*(zeta(i,Jend  ,know)+                       &
     &                       zeta(i,Jend+1,know))-                      &
     &               BOUNDARY(ng)%zeta_north(i))*ad_vbar(i,Jend+1,kout)
              ad_bry_val=ad_bry_val+ad_vbar(i,Jend+1,kout)
              ad_vbar(i,Jend+1,kout)=0.0_r8
!>            tl_Ce=0.5_r8*g*tl_cff/Ce
!>
              ad_cff=ad_cff+0.5_r8*g*ad_Ce/Ce
              ad_Ce=0.0_r8
!>            tl_cff=-cff*cff*(0.5_r8*(GRID(ng)%tl_h(i,Jend  )+         &
!>   &                                 tl_zeta(i,Jend  ,know)+          &
!>   &                                 GRID(ng)%tl_h(i,Jend+1)+         &
!>   &                                 tl_zeta(i,Jend+1,know)))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              GRID(ng)%ad_h(i,Jend  )=GRID(ng)%ad_h(i,Jend  )+adfac
              GRID(ng)%ad_h(i,Jend+1)=GRID(ng)%ad_h(i,Jend+1)+adfac
              ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)+adfac
              ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
              ad_cff=0.0_r8

# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (ad_LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                 &
     &                      zeta(i,Jend,know))*                         &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
              ELSE
                bry_pgr=-g*(zeta(i,Jend+1,know)-                        &
     &                      zeta(i,Jend  ,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
              END IF
#  ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jend  ,know)+                 &
     &                           ubar(i+1,Jend  ,know)+                 &
     &                           ubar(i  ,Jend+1,know)+                 &
     &                           ubar(i+1,Jend+1,know))*                &
     &                          (GRID(ng)%f(i,Jend  )+                  &
     &                           GRID(ng)%f(i,Jend+1))
#  else
              bry_cor=0.0_r8
#  endif
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
!>            tl_bry_val=tl_vbar(i,Jend,know)+                          &
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
              ad_vbar(i,Jend,know)=ad_vbar(i,Jend,know)+ad_bry_val
              ad_bry_val=0.0_r8
!>            tl_cff2=GRID(ng)%on_v(i,Jend+1)*tl_Ce
!>
              ad_Ce=ad_Ce+GRID(ng)%on_v(i,Jend+1)*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jend+1)+       &
!>   &                                   tl_zeta(i,Jend+1,know)+        &
!>   &                                   GRID(ng)%tl_h(i,Jend  )+       &
!>   &                                   tl_zeta(i,Jend  ,know))
!>
              adfac=-Ce*Ce*Ce*0.25_r8*g*ad_Ce
              ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)+adfac
              ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
              GRID(ng)%ad_h(i,Jend  )=GRID(ng)%ad_h(i,Jend  )+adfac
              GRID(ng)%ad_h(i,Jend+1)=GRID(ng)%ad_h(i,Jend+1)+adfac
              ad_Ce=0.0_r8
!>            tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jend+1)-           &
!>   &                            FORCES(ng)%bvstr(i,Jend+1))+          &
!>   &                   cff1*(FORCES(ng)%tl_svstr(i,Jend+1)-           &
!>   &                         FORCES(ng)%tl_bvstr(i,Jend+1))
!>
              adfac=cff1*ad_bry_str
              FORCES(ng)%ad_svstr(i,Jend+1)=                            &
     &                                   FORCES(ng)%ad_svstr(i,Jend+1)+ &
     &                                      adfac
              FORCES(ng)%ad_bvstr(i,Jend+1)=                            &
     &                                   FORCES(ng)%ad_bvstr(i,Jend+1)- &
     &                                    adfac
              ad_cff1=ad_cff1+(FORCES(ng)%svstr(i,Jend+1)-              &
     &                         FORCES(ng)%bvstr(i,Jend+1))*ad_bry_str

              ad_bry_str=0.0_r8
!>            tl_cff1=-cff1*cff1*0.5_r8*(GRID(ng)%tl_h(i,Jend  )+       &
!>   &                                   tl_zeta(i,Jend  ,know)+        &
!>   &                                   GRID(ng)%tl_h(i,Jend+1)+       &
!>   &                                   tl_zeta(i,Jend+1,know))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)+adfac
              ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
              GRID(ng)%ad_h(i,Jend  )=GRID(ng)%ad_h(i,Jend  )+adfac
              GRID(ng)%ad_h(i,Jend+1)=GRID(ng)%ad_h(i,Jend+1)+adfac
              ad_cff1=0.0_r8
#  ifdef UV_COR
!>            tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jend  ,know)+           &
!>   &                              tl_ubar(i+1,Jend  ,know)+           &
!>   &                              tl_ubar(i  ,Jend+1,know)+           &
!>   &                              tl_ubar(i+1,Jend+1,know))*          &
!>   &                             (GRID(ng)%f(i,Jend  )+               &
!>   &                              GRID(ng)%f(i,Jend+1))
!>
              adfac=0.125_r8*(GRID(ng)%f(i,Jend  )+                     &
     &                        GRID(ng)%f(i,Jend+1))*ad_bry_cor
              ad_ubar(i  ,Jend  ,know)=ad_ubar(i  ,Jend  ,know)-adfac
              ad_ubar(i+1,Jend  ,know)=ad_ubar(i+1,Jend  ,know)-adfac
              ad_ubar(i  ,Jend+1,know)=ad_ubar(i  ,Jend+1,know)-adfac
              ad_ubar(i+1,Jend+1,know)=ad_ubar(i+1,Jend+1,know)-adfac
              ad_bry_cor=0.0_r8
#  else
!>            tl_bry_cor=0.0_r8
!>
#  endif
              IF (ad_LBC(inorth,isFsur,ng)%acquire) THEN
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(inorth,isVbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr-                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_north(i)*           &
!>   &                       0.5_r8*GRID(ng)%pn(i,Jend)
!>
                  BOUNDARY(ng)%ad_zeta_north(i)=BOUNDARY(ng)%           &
     &                                                ad_zeta_north(i)- &
     &                                          g*0.5_r8*               &
     &                                          GRID(ng)%pn(i,Jend)*    &
     &                                          ad_bry_pgr
                END IF
#  endif
!>              tl_bry_pgr=g*tl_zeta(i,Jend,know)*                      &
!>   &                     0.5_r8*GRID(ng)%pn(i,Jend)
!>
                ad_zeta(i,Jend,know)=ad_zeta(i,Jend,know)+              &
     &                               g*0.5_r8*GRID(ng)%pn(i,Jend)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(i,Jend+1,know)-                  &
!>   &                         tl_zeta(i,Jend  ,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pn(i,Jend  )+               &
!>   &                             GRID(ng)%pn(i,Jend+1))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pn(i,Jend  )+                 &
     &                           GRID(ng)%pn(i,Jend+1))*ad_bry_pgr
                ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)-adfac
                ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
# else
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
!>              tl_bry_val=BOUNDARY(ng)%tl_vbar_north(i)
!>
                BOUNDARY(ng)%ad_vbar_north(i)=BOUNDARY(ng)%             &
     &                                              ad_vbar_north(i)+   &
     &                                        ad_bry_val
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
!  Northern edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (ad_LBC(inorth,isVbar,ng)%Shchepetkin) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=0.5_r8*(GRID(ng)%h(i,Jend  )+                         &
     &                    GRID(ng)%h(i,Jend+1))
              cff1=SQRT(g/cff)
              Ce=dt2d*cff1*cff*0.5_r8*(GRID(ng)%pn(i,Jend  )+           &
     &                                 GRID(ng)%pn(i,Jend+1))
              Ze=(0.5_r8+Ce)*zeta(i,Jend  ,know)+                       &
     &           (0.5_r8-Ce)*zeta(i,Jend+1,know)
              IF (Ce.gt.Co) THEN
                cff2=(1.0_r8-Co/Ce)**2
                cff3=zeta(i,Jend,kout)+                                 &
     &               Ce*zeta(i,Jend+1,know)-                            &
     &               (1.0_r8+Ce)*zeta(i,Jend,know)
                Ze=Ze+cff2*cff3
              END IF
# ifdef MASKING
!>            tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
!>   &                               GRID(ng)%vmask(i,Jend+1)
!>
              ad_vbar(i,Jend+1,kout)=ad_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
!>              tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)-          &
!>   &                                 Ce*BOUNDARY(ng)%tl_zeta_north(i)
!>
                BOUNDARY(ng)%ad_zeta_north(i)=BOUNDARY(ng)%             &
     &                                              ad_zeta_north(i)-   &
     &                                        Ce*ad_vbar(i,Jend+1,kout)
              END IF
# endif
!>            tl_vbar(i,Jend+1,kout)=0.5_r8*                            &
!>   &                               ((1.0_r8-Ce)*                      &
!>   &                                tl_vbar(i,Jend+1,know)+           &
!>   &                                tl_Ce*(vbar(i,Jend  ,know)-       &
!>   &                                       vbar(i,Jend+1,know))+      &
!>   &                                Ce*tl_vbar(i,Jend,know)+          &
!>   &                                tl_bry_val+                       &
!>   &                                tl_cff1*                          &
!>   &                                (Ze-BOUNDARY(ng)%zeta_north(i))-  &
!>   &                                 cff1*tl_ze)
!>
              adfac=0.5_r8*ad_vbar(i,Jend+1,kout)
              ad_vbar(i,Jend+1,know)=ad_vbar(i,Jend+1,know)+            &
     &                               (1.0_r8-Ce)*adfac
              ad_vbar(i,Jend  ,know)=ad_vbar(i,Jend  ,know)+            &
     &                               Ce*adfac
              ad_Ce=ad_Ce+                                              &
     &              (vbar(i,Jend  ,know)-                               &
     &               vbar(i,Jend+1,know))*adfac
              ad_cff1=ad_cff1+                                          &
     &                (Ze-BOUNDARY(ng)%zeta_north(i))*adfac
              ad_ze=ad_ze-cff1*adfac
              ad_vbar(i,Jend+1,kout)=0.0_r8

              IF (Ce.gt.Co) THEN
!>              tl_Ze=tl_Ze+cff2*tl_cff3+                               &
!>   &                      tl_cff2*cff3
!>
                ad_cff2=ad_cff2+cff3*ad_Ze
                ad_cff3=ad_cff3+cff2*ad_Ze
!>              tl_cff3=tl_zeta(i,Jend,kout)+                           &
!>   &                  Ce*tl_zeta(i,Jend+1,know)+                      &
!>   &                  tl_Ce*(zeta(i,Jend  ,know)+                     &
!>   &                         zeta(i,Jend+1,know))-                    &
!>   &                  (1.0_r8+Ce)*tl_zeta(i,Jend,know)
!>
                ad_zeta(i,Jend  ,kout)=ad_zeta(i,Jend  ,kout)+          &
     &                                 ad_cff3
                ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)-          &
     &                                 (1.0_r8+Ce)*ad_cff3
                ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+          &
     &                                 Ce*ad_cff3
                ad_Ce=ad_Ce+                                            &
     &                (zeta(i,Jend  ,know)+                             &
     &                 zeta(i,Jend+1,know))*ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=2.0_r8*cff2*Co*tl_Ce/(Ce*Ce)
!>
                ad_Ce=ad_Ce+                                            &
     &                2.0_r8*cff2*Co*ad_cff2/(Ce*Ce)
                ad_cff2=0.0_r8
              END IF
!>            tl_Ze=(0.5_r8+Ce)*tl_zeta(i,Jend  ,know)+                 &
!>   &              (0.5_r8-Ce)*tl_zeta(i,Jend+1,know)+                 &
!>   &              tl_Ce*(zeta(i,Jend  ,know)-                         &
!>   &                     zeta(i,Jend+1,know))
!>
              ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)+            &
     &                               (0.5_r8+Ce)*ad_Ze
              ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+            &
     &                               (0.5_r8-Ce)*ad_Ze
              ad_Ce=ad_Ce+                                              &
     &              (zeta(i,Jend  ,know)-                               &
     &               zeta(i,Jend+1,know))*ad_Ze
              ad_Ze=0.0_r8
!>            tl_Ce=dt2d*0.5_r8*(GRID(ng)%pn(i,Jend  )+                 &
!>   &                           GRID(ng)%pn(i,Jend+1))*                &
!>   &              (cff1*tl_cff+                                       &
!>   &               tl_cff1*cff)
!>
              adfac=dt2d*0.5_r8*(GRID(ng)%pn(i,Jend  )+                 &
     &                           GRID(ng)%pn(i,Jend+1))*ad_Ce
              ad_cff=ad_cff+cff1*adfac
              ad_cff1=ad_cff1+cff*adfac
              ad_Ce=0.0_r8
!>            tl_cff1=-0.5_r8*cff1*tl_cff/cff
!>
              ad_cff=ad_cff-                                            &
     &               0.5_r8*cff1*ad_cff/cff
              ad_cff1=0.0_r8
!>            tl_cff=0.5_r8*(GRID(ng)%tl_h(i,Jend  )+                   &
!>   &                       GRID(ng)%tl_h(i,Jend+1))
!>
              adfac=0.5_r8*ad_cff
              GRID(ng)%ad_h(i,Jend  )=GRID(ng)%ad_h(i,Jend  )+adfac
              GRID(ng)%ad_h(i,Jend+1)=GRID(ng)%ad_h(i,Jend+1)+adfac
              ad_cff=0.0_r8

# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (ad_LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                 &
     &                      zeta(i,Jend,know))*                         &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
              ELSE
                bry_pgr=-g*(zeta(i,Jend+1,know)-                        &
     &                      zeta(i,Jend  ,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
              END IF
#  ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jend  ,know)+                 &
     &                           ubar(i+1,Jend  ,know)+                 &
     &                           ubar(i  ,Jend+1,know)+                 &
     &                           ubar(i+1,Jend+1,know))*                &
     &                          (GRID(ng)%f(i,Jend  )+                  &
     &                           GRID(ng)%f(i,Jend+1))
#  else
              bry_cor=0.0_r8
#  endif
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
!>            tl_bry_val=tl_vbar(i,Jend,know)+                          &
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
              ad_vbar(i,Jend,know)=ad_vbar(i,Jend,know)+ad_bry_val
              ad_bry_val=0.0_r8
!>            tl_cff2=GRID(ng)%on_v(i,Jend+1)*tl_Ce
!>
              ad_Ce=ad_Ce+GRID(ng)%on_v(i,Jend+1)*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jend+1)+       &
!>   &                                   tl_zeta(i,Jend+1,know)+        &
!>   &                                   GRID(ng)%tl_h(i,Jend  )+       &
!>   &                                   tl_zeta(i,Jend  ,know))
!>
              adfac=-Ce*Ce*Ce*0.25_r8*g*ad_Ce
              ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)+adfac
              ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
              GRID(ng)%ad_h(i,Jend  )=GRID(ng)%ad_h(i,Jend  )+adfac
              GRID(ng)%ad_h(i,Jend+1)=GRID(ng)%ad_h(i,Jend+1)+adfac
              ad_Ce=0.0_r8
!>            tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jend+1)-           &
!>   &                            FORCES(ng)%bvstr(i,Jend+1))+          &
!>   &                   cff1*(FORCES(ng)%tl_svstr(i,Jend+1)-           &
!>   &                         FORCES(ng)%tl_bvstr(i,Jend+1))
!>
              adfac=cff1*ad_bry_str
              FORCES(ng)%ad_svstr(i,Jend+1)=                            &
     &                                   FORCES(ng)%ad_svstr(i,Jend+1)+ &
     &                                      adfac
              FORCES(ng)%ad_bvstr(i,Jend+1)=                            &
     &                                   FORCES(ng)%ad_bvstr(i,Jend+1)- &
     &                                    adfac
              ad_cff1=ad_cff1+(FORCES(ng)%svstr(i,Jend+1)-              &
     &                         FORCES(ng)%bvstr(i,Jend+1))*ad_bry_str

              ad_bry_str=0.0_r8
!>            tl_cff1=-cff1*cff1*0.5_r8*(GRID(ng)%tl_h(i,Jend  )+       &
!>   &                                   tl_zeta(i,Jend  ,know)+        &
!>   &                                   GRID(ng)%tl_h(i,Jend+1)+       &
!>   &                                   tl_zeta(i,Jend+1,know))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)+adfac
              ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
              GRID(ng)%ad_h(i,Jend  )=GRID(ng)%ad_h(i,Jend  )+adfac
              GRID(ng)%ad_h(i,Jend+1)=GRID(ng)%ad_h(i,Jend+1)+adfac
              ad_cff1=0.0_r8
#  ifdef UV_COR
!>            tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jend  ,know)+           &
!>   &                              tl_ubar(i+1,Jend  ,know)+           &
!>   &                              tl_ubar(i  ,Jend+1,know)+           &
!>   &                              tl_ubar(i+1,Jend+1,know))*          &
!>   &                             (GRID(ng)%f(i,Jend  )+               &
!>   &                              GRID(ng)%f(i,Jend+1))
!>
              adfac=0.125_r8*(GRID(ng)%f(i,Jend  )+                     &
     &                        GRID(ng)%f(i,Jend+1))*ad_bry_cor
              ad_ubar(i  ,Jend  ,know)=ad_ubar(i  ,Jend  ,know)-adfac
              ad_ubar(i+1,Jend  ,know)=ad_ubar(i+1,Jend  ,know)-adfac
              ad_ubar(i  ,Jend+1,know)=ad_ubar(i  ,Jend+1,know)-adfac
              ad_ubar(i+1,Jend+1,know)=ad_ubar(i+1,Jend+1,know)-adfac
              ad_bry_cor=0.0_r8
#  else
!>            tl_bry_cor=0.0_r8
!>
#  endif
              IF (ad_LBC(inorth,isFsur,ng)%acquire) THEN
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(inorth,isVbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr-                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_north(i)*           &
!>   &                       0.5_r8*GRID(ng)%pn(i,Jend)
!>
                  BOUNDARY(ng)%ad_zeta_north(i)=BOUNDARY(ng)%           &
     &                                                ad_zeta_north(i)- &
     &                                          g*0.5_r8*               &
     &                                          GRID(ng)%pn(i,Jend)*    &
     &                                          ad_bry_pgr
                END IF
#  endif
!>              tl_bry_pgr=g*tl_zeta(i,Jend,know)*                      &
!>   &                     0.5_r8*GRID(ng)%pn(i,Jend)
!>
                ad_zeta(i,Jend,know)=ad_zeta(i,Jend,know)+              &
     &                               g*0.5_r8*GRID(ng)%pn(i,Jend)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(i,Jend+1,know)-                  &
!>   &                         tl_zeta(i,Jend  ,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pn(i,Jend  )+               &
!>   &                             GRID(ng)%pn(i,Jend+1))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pn(i,Jend  )+                 &
     &                           GRID(ng)%pn(i,Jend+1))*ad_bry_pgr
                ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)-adfac
                ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
# else
#  ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
!>              tl_bry_val=BOUNDARY(ng)%tl_vbar_north(i)
!>
                BOUNDARY(ng)%ad_vbar_north(i)=BOUNDARY(ng)%             &
     &                                              ad_vbar_north(i)+   &
     &                                        ad_bry_val
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
!  Northern edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(inorth,isVbar,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
# ifdef MASKING
!>            tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
!>   &                               GRID(ng)%vmask(i,Jend+1)
!>
              ad_vbar(i,Jend+1,kout)=ad_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(inorth,isVbar,ng)) THEN
!>              tl_vbar(i,Jend+1,kout)=BOUNDARY(ng)%tl_vbar_north(i)
!>
                BOUNDARY(ng)%ad_vbar_north(i)=BOUNDARY(ng)%             &
     &                                              ad_vbar_north(i)+   &
     &                                        ad_vbar(i,Jend+1,kout)
                ad_vbar(i,Jend+1,kout)=0.0_r8
              ELSE
!>              tl_vbar(i,Jend+1,kout)=0.0_r8
!>
                ad_vbar(i,Jend+1,kout)=0.0_r8
              END IF
# else
!>            tl_vbar(i,Jend+1,kout)=0.0_r8
!>
              ad_vbar(i,Jend+1,kout)=0.0_r8
# endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(inorth,isVbar,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
# ifdef MASKING
!>            tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
!>   &                               GRID(ng)%vmask(i,Jend+1)
!>
              ad_vbar(i,Jend+1,kout)=ad_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
!>            tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend,kout)
!>
              ad_vbar(i,Jend  ,kout)=ad_vbar(i,Jend,kout)+              &
     &                               ad_vbar(i,Jend+1,kout)
              ad_vbar(i,Jend+1,kout)=0.0_r8
            END IF
          END DO
!
!  Northern edge, reduced-physics boundary condition.
!
        ELSE IF (ad_LBC(inorth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                 &
     &                            zeta(i,Jend  ,know)+                  &
     &                            GRID(ng)%h(i,Jend+1)+                 &
     &                            zeta(i,Jend+1,know)))
# ifdef MASKING
!>            tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,kout)*            &
!>   &                               GRID(ng)%vmask(i,Jend+1)
!>
              ad_vbar(i,Jend+1,kout)=ad_vbar(i,Jend+1,kout)*            &
     &                               GRID(ng)%vmask(i,Jend+1)
# endif
!>            tl_vbar(i,Jend+1,kout)=tl_vbar(i,Jend+1,know)+            &
!>   &                               dt2d*(tl_bry_pgr+                  &
!>   &                                     tl_bry_cor+                  &
!>   &                                     tl_bry_str)
!>
              adfac=dt2d*ad_vbar(i,Jend+1,kout)
              ad_bry_pgr=ad_bry_pgr+adfac
              ad_bry_cor=ad_bry_cor+adfac
              ad_bry_str=ad_bry_str+adfac
              ad_vbar(i,Jend+1,know)=ad_vbar(i,Jend+1,know)+            &
     &                               ad_vbar(i,Jend+1,kout)
              ad_vbar(i,Jend+1,kout)=0.0_r8
!>            tl_bry_str=tl_cff*(FORCES(ng)%svstr(i,Jend+1)-            &
!>   &                           FORCES(ng)%bvstr(i,Jend+1))+           &
!>   &                   cff*(FORCES(ng)%tl_svstr(i,Jend+1)-            &
!>   &                        FORCES(ng)%tl_bvstr(i,Jend+1))
!>
              adfac=cff*ad_bry_str
              FORCES(ng)%ad_svstr(i,Jend+1)=                            &
     &                                   FORCES(ng)%ad_svstr(i,Jend+1)+ &
     &                                      adfac
              FORCES(ng)%ad_bvstr(i,Jend+1)=                            &
     &                                   FORCES(ng)%ad_bvstr(i,Jend+1)- &
     &                                      adfac
              ad_cff=ad_cff+(FORCES(ng)%svstr(i,Jend+1)-                &
     &                       FORCES(ng)%bvstr(i,Jend+1))*ad_bry_str
              ad_bry_str=0.0_r8
!>            tl_cff=-cff*cff*0.5_r8*(GRID(ng)%tl_h(i,Jend  )+          &
!>   &                                tl_zeta(i,Jend  ,know)+           &
!>   &                                GRID(ng)%tl_h(i,Jend+1)+          &
!>   &                                tl_zeta(i,Jend+1,know))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)+adfac
              ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
              GRID(ng)%ad_h(i,Jend  )=GRID(ng)%ad_h(i,Jend  )+adfac
              GRID(ng)%ad_h(i,Jend+1)=GRID(ng)%ad_h(i,Jend+1)+adfac
              ad_cff=0.0_r8
# ifdef UV_COR
!>            tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jend  ,know)+           &
!>   &                              tl_ubar(i+1,Jend  ,know)+           &
!>   &                              tl_ubar(i  ,Jend+1,know)+           &
!>   &                              tl_ubar(i+1,Jend+1,know))*          &
!>   &                             (GRID(ng)%f(i,Jend  )+               &
!>   &                              GRID(ng)%f(i,Jend+1))
!>
              adfac=-0.125_r8*(GRID(ng)%f(i,Jend  )+                    &
     &                         GRID(ng)%f(i,Jend+1))*ad_bry_cor
              ad_ubar(i  ,Jend  ,know)=ad_ubar(i  ,Jend  ,know)+adfac
              ad_ubar(i+1,Jend  ,know)=ad_ubar(i+1,Jend  ,know)+adfac
              ad_ubar(i  ,Jend+1,know)=ad_ubar(i  ,Jend+1,know)+adfac
              ad_ubar(i+1,Jend+1,know)=ad_ubar(i+1,Jend+1,know)+adfac
              ad_bry_cor=0.0_r8
# else
!>            tl_bry_cor=0.0_r8
!>
              ad_bry_cor=0.0_r8
# endif
              IF (ad_LBC(inorth,isFsur,ng)%acquire) THEN
# ifdef ADJUST_BOUNDARY
                IF (Lobc(inorth,isVbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr-                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_north(i)*           &
!>   &                      0.5_r8*GRID(ng)%pn(i,Jend)
!>
                  BOUNDARY(ng)%ad_zeta_north(i)=BOUNDARY(ng)%           &
     &                                                ad_zeta_north(i)- &
     &                                          g*0.5_r8*               &
     &                                          GRID(ng)%pn(i,Jend)*    &
     &                                          ad_bry_pgr
                END IF
# endif
!>              tl_bry_pgr=g*tl_zeta(i,Jend,know)*                      &
!>   &                     0.5_r8*GRID(ng)%pn(i,Jend)
!>
                ad_zeta(i,Jend,know)=ad_zeta(i,Jend,know)+              &
     &                               g*0.5_r8*GRID(ng)%pn(i,Jend)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(i,Jend+1,know)-                  &
!>   &                         tl_zeta(i,Jend  ,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pn(i,Jend  )+               &
!>   &                             GRID(ng)%pn(i,Jend+1))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pn(i,Jend  )+                 &
     &                           GRID(ng)%pn(i,Jend+1))*ad_bry_pgr
                ad_zeta(i,Jend  ,know)=ad_zeta(i,Jend  ,know)-adfac
                ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (ad_LBC(inorth,isVbar,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
!>            tl_vbar(i,Jend+1,kout)=0.0_r8
!>
              ad_vbar(i,Jend+1,kout)=0.0_r8
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
        IF (ad_LBC(isouth,isVbar,ng)%radiation) THEN
          IF (iic(ng).ne.0) THEN
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
# if defined CELERITY_READ && defined FORWARD_READ
                IF (ad_LBC(isouth,isVbar,ng)%nudging) THEN
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
# ifdef MASKING
!>              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*              &
!>   &                               GRID(ng)%vmask(i,Jstr)
!>
                ad_vbar(i,Jstr,kout)=ad_vbar(i,Jstr,kout)*              &
     &                               GRID(ng)%vmask(i,Jstr)
# endif
                IF (ad_LBC(isouth,isVbar,ng)%nudging) THEN
!>                tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)-            &
!>   &                                 tau*tl_vbar(i,Jstr,know)
!>
                  ad_vbar(i,Jstr,know)=ad_vbar(i,Jstr,know)-            &
     &                                 tau*ad_vbar(i,Jstr,kout)
                END IF
!>              tl_vbar(i,Jstr,kout)=(cff*tl_vbar(i,Jstr  ,know)+       &
!>   &                                Ce *tl_vbar(i,Jstr+1,kout)-       &
!>   &                                MAX(Cx,0.0_r8)*                   &
!>   &                                   tl_grad(i  ,Jstr)-             &
!>   &                                MIN(Cx,0.0_r8)*                   &
!>   &                                   tl_grad(i+1,Jstr))/            &
!>   &                               (cff+Ce)
!>
                adfac=ad_vbar(i,Jstr,kout)/(cff+Ce)
                ad_grad(i  ,Jstr)=ad_grad(i  ,Jstr)-MAX(Cx,0.0_r8)*adfac
                ad_grad(i+1,Jstr)=ad_grad(i+1,Jstr)-MIN(Cx,0.0_r8)*adfac
                ad_vbar(i,Jstr  ,know)=ad_vbar(i,Jstr  ,know)+cff*adfac
                ad_vbar(i,Jstr-1,kout)=ad_vbar(i,Jstr-1,kout)+Ce *adfac
                ad_vbar(i,Jstr  ,kout)=0.0_r8
              END IF
            END DO
          END IF
!
!  Southern edge, Flather boundary condition.
!
        ELSE IF (ad_LBC(isouth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                 &
     &                            zeta(i,Jstr-1,know)+                  &
     &                            GRID(ng)%h(i,Jstr  )+                 &
     &                            zeta(i,Jstr  ,know)))
              Ce=SQRT(g*cff)
# ifdef MASKING
!>            tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
!>   &                             GRID(ng)%vmask(i,Jstr)
!>
              ad_vbar(i,Jstr,kout)=ad_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
!>              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)+              &
!>   &                               Ce*BOUNDARY(ng)%tl_zeta_south(i)
!>
                BOUNDARY(ng)%ad_zeta_south(i)=BOUNDARY(ng)%             &
     &                                              ad_zeta_south(i)+   &
     &                                        Ce*ad_vbar(i,Jstr,kout)
              END IF
# endif
!>            tl_vbar(i,Jstr,kout)=tl_bry_val-                          &
!>   &                             tl_Ce*                               &
!>   &                             (0.5_r8*(zeta(i,Jstr-1,know)+        &
!>   &                                      zeta(i,Jstr  ,know))-       &
!>   &                              BOUNDARY(ng)%zeta_south(i))-        &
!>   &                             Ce*                                  &
!>   &                             (0.5_r8*(tl_zeta(i,Jstr-1,know)+     &
!>   &                                     tl_zeta(i,Jstr  ,know)))
!>
              adfac=Ce*0.5_r8*ad_vbar(i,Jstr,kout)
              ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)-adfac
              ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)-adfac
              ad_Ce=ad_Ce-                                              &
     &              (0.5_r8*(zeta(i,Jstr-1,know)+                       &
     &                       zeta(i,Jstr  ,know))-                      &
     &               BOUNDARY(ng)%zeta_south(i))*ad_vbar(i,Jstr,kout)
              ad_bry_val=ad_bry_val+ad_vbar(i,Jstr,kout)
              ad_vbar(i,Jstr,kout)=0.0_r8
!>            tl_Ce=0.5_r8*g*tl_cff/Ce
!>
              ad_cff=ad_cff+0.5_r8*g*ad_Ce/Ce
              ad_Ce=0.0_r8
!>            tl_cff=-cff*cff*(0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+         &
!>   &                                 tl_zeta(i,Jstr-1,know)+          &
!>   &                                 GRID(ng)%tl_h(i,Jstr  )+         &
!>   &                                 tl_zeta(i,Jstr  ,know)))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              GRID(ng)%ad_h(i,Jstr-1)=GRID(ng)%ad_h(i,Jstr-1)+adfac
              GRID(ng)%ad_h(i,Jstr  )=GRID(ng)%ad_h(i,Jstr  )+adfac
              ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+adfac
              ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+adfac
              ad_cff=0.0_r8

# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (ad_LBC(isouth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(i,Jstr,know)-                          &
     &                      BOUNDARY(ng)%zeta_south(i))*                &
     &                  0.5_r8*GRID(ng)%pn(i,Jstr)
              ELSE
                bry_pgr=-g*(zeta(i,Jstr  ,know)-                        &
     &                      zeta(i,Jstr-1,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                  &
     &                          GRID(ng)%pn(i,Jstr  ))
              END IF
#  ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jstr-1,know)+                 &
     &                           ubar(i+1,Jstr-1,know)+                 &
     &                           ubar(i  ,Jstr  ,know)+                 &
     &                           ubar(i+1,Jstr  ,know))*                &
     &                          (GRID(ng)%f(i,Jstr-1)+                  &
     &                           GRID(ng)%f(i,Jstr  ))
#  else
              bry_cor=0.0_r8
#  endif
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
!>            tl_bry_val=tl_vbar(i,Jstr+1,know)+                        &
!>   &                   tl_cff2*(bry_pgr+                              &
!>   &                            bry_cor+                              &
!>   &                            bry_str)+                             &
!>   &                   cff2*(tl_bry_pgr+                              &
!>   &                         tl_bry_cor+                              &
!>   &                         tl_bry_str)
!>
              adfac=cff2*ad_bry_val
              tl_bry_pgr=tl_bry_pgr+adfac
              tl_bry_cor=tl_bry_cor+adfac
              tl_bry_str=tl_bry_str+adfac
              ad_cff2=ad_cff2+(bry_pgr+                                 &
     &                         bry_cor+                                 &
     &                         bry_str)*ad_bry_val
              ad_vbar(i,Jstr+1,know)=ad_vbar(i,Jstr+1,know)+ad_bry_val
              ad_bry_val=0.0_r8
!>            tl_cff2=GRID(ng)%on_v(i,Jstr)*tl_Ce
!>
              ad_Ce=ad_Ce+GRID(ng)%on_v(i,Jstr)*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jstr-1)+       &
!>   &                                   tl_zeta(i,Jstr-1,know)+        &
!>   &                                   GRID(ng)%tl_h(i,Jstr  )+       &
!>   &                                   tl_zeta(i,Jstr  ,know))
!>
              adfac=-Ce*Ce*Ce*0.25_r8*g*ad_Ce
              ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+adfac
              ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+adfac
              GRID(ng)%ad_h(i,Jstr-1)=GRID(ng)%ad_h(i,Jstr-1)+adfac
              GRID(ng)%ad_h(i,Jstr  )=GRID(ng)%ad_h(i,Jstr  )+adfac
              ad_Ce=0.0_r8
!>            tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jstr)-             &
!>   &                            FORCES(ng)%bvstr(i,Jstr))+            &
!>   &                   cff1*(FORCES(ng)%tl_svstr(i,Jstr)-             &
!>   &                         FORCES(ng)%tl_bvstr(i,Jstr))
!>
              adfac=cff1*ad_bry_str
              FORCES(ng)%ad_svstr(i,Jstr)=FORCES(ng)%ad_svstr(i,Jstr)+  &
     &                                    adfac
              FORCES(ng)%ad_bvstr(i,Jstr)=FORCES(ng)%ad_bvstr(i,Jstr)-  &
     &                                    adfac
              ad_cff1=ad_cff1+(FORCES(ng)%svstr(i,Jstr)-                &
     &                         FORCES(ng)%bvstr(i,Jstr))*ad_bry_str
              ad_bry_str=0.0_r8
!>            tl_cff1=-cff1*cff1*(0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+      &
!>   &                                    tl_zeta(i,Jstr-1,know)+       &
!>   &                                    GRID(ng)%tl_h(i,Jstr  )+      &
!>   &                                    tl_zeta(i,Jstr  ,know)))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+adfac
              ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+adfac
              GRID(ng)%ad_h(i,Jstr-1)=GRID(ng)%ad_h(i,Jstr-1)+adfac
              GRID(ng)%ad_h(i,Jstr  )=GRID(ng)%ad_h(i,Jstr  )+adfac
              ad_cff1=0.0_r8
#  ifdef UV_COR
!>            tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jstr-1,know)+           &
!>   &                              tl_ubar(i+1,Jstr-1,know)+           &
!>   &                              tl_ubar(i  ,Jstr  ,know)+           &
!>   &                              tl_ubar(i+1,Jstr  ,know))*          &
!>   &                             (GRID(ng)%f(i,Jstr-1)+               &
!>   &                              GRID(ng)%f(i,Jstr  ))
!>
              adfac=-0.125_r8*(GRID(ng)%f(i,Jstr-1)+                    &
     &                         GRID(ng)%f(i,Jstr  ))*ad_bry_cor
              ad_ubar(i  ,Jstr-1,know)=ad_ubar(i  ,Jstr-1,know)+adfac
              ad_ubar(i+1,Jstr-1,know)=ad_ubar(i+1,Jstr-1,know)+adfac
              ad_ubar(i  ,Jstr  ,know)=ad_ubar(i  ,Jstr  ,know)+adfac
              ad_ubar(i+1,Jstr  ,know)=ad_ubar(i+1,Jstr  ,know)+adfac
              ad_bry_cor=0.0_r8
#  else
!>            tl_bry_cor=0.0_r8
!>
              ad_bry_cor=0.0_r8
#  endif
              IF (ad_LBC(isouth,isFsur,ng)%acquire) THEN
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(isouth,isVbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr+                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_south(i)*           &
!>   &                       0.5_r8*GRID(ng)%pn(i,Jstr)
!>
                  BOUNDARY(ng)%ad_zeta_south(i)=BOUNDARY(ng)%           &
     &                                                ad_zeta_south(i)+ &
     &                                          g*0.5_r8*               &
     &                                          GRID(ng)%pn(i,Jstr)*    &
     &                                          ad_bry_pgr
                END IF
#  endif
!>              tl_bry_pgr=-g*tl_zeta(i,Jstr,know)*                     &
!>   &                     0.5_r8*GRID(ng)%pn(i,Jstr)
!>
                tl_zeta(i,Jstr,know)=tl_zeta(i,Jstr,know)-              &
     &                               g*0.5_r8*GRID(ng)%pn(i,Jstr)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(i,Jstr  ,know)-                  &
!>   &                         tl_zeta(i,Jstr-1,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pn(i,Jstr-1)+               &
!>   &                             GRID(ng)%pn(i,Jstr  ))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                 &
     &                           GRID(ng)%pn(i,Jstr  ))*ad_bry_pgr
                ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)-adfac
                ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
# else
#   ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
!>              tl_bry_val=BOUNDARY(ng)%tl_vbar_south(i)
!>
                BOUNDARY(ng)%ad_vbar_south(i)=BOUNDARY(ng)%             &
     &                                              ad_vbar_south(i)+   &
     &                                        ad_bry_val
                ad_bry_val=0.0_r8
              ELSE
!>              tl_bry_val=0.0_r8
!>
                ad_bry_val=0.0_r8
              END IF
#   else
!>            tl_bry_val=0.0_r8
!>
              ad_bry_val=0.0_r8
#   endif
# endif
            END IF
          END DO
!
!  Southern edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (ad_LBC(isouth,isVbar,ng)%Shchepetkin) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=0.5_r8*(GRID(ng)%h(i,Jstr-1)+                         &
     &                    GRID(ng)%h(i,Jstr  ))
              cff1=SQRT(g/cff)
              Ce=dt2d*cff1*cff*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+           &
     &                                 GRID(ng)%pn(i,Jstr  ))
              Ze=(0.5_r8+Ce)*zeta(i,Jstr  ,know)+                       &
     &           (0.5_r8-Ce)*zeta(i,Jstr-1,know)
              IF (Ce.gt.Co) THEN
                cff2=(1.0_r8-Co/Ce)**2
                cff3=zeta(i,Jstr,kout)+                                 &
     &               Ce*zeta(i,Jstr-1,know)-                            &
     &               (1.0_r8+Ce)*zeta(i,Jstr,know)
                Ze=Ze+cff2*cff3
              END IF
# ifdef MASKING
!>            tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
!>   &                             GRID(ng)%vmask(i,Jstr)
!>
              ad_vbar(i,Jstr,kout)=ad_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
!>              tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)+              &
!>   &                               Ce*BOUNDARY(ng)%tl_zeta_south(i)
!>
                BOUNDARY(ng)%ad_zeta_south(i)=BOUNDARY(ng)%             &
     &                                              ad_zeta_south(i)+   &
     &                                        Ce*ad_vbar(i,Jstr,kout)
              END IF
# endif
!>            tl_vbar(i,Jstr,kout)=0.5_r8*                              &
!>   &                             ((1.0_r8-Ce)*                        &
!>   &                              tl_vbar(i,Jstr,know)-               &
!>   &                              tl_Ce*(vbar(i,Jstr  ,know)-         &
!>   &                                     vbar(i,Jstr+1,know))+        &
!>   &                              Ce*tl_vbar(i,Jstr+1,know)+          &
!>   &                              tl_bry_val-                         &
!>   &                              tl_cff1*                            &
!>   &                              (Ze-BOUNDARY(ng)%zeta_south(i))-    &
!>   &                              cff1*tl_Ze)
!>
              adfac=0.5_r8*ad_vbar(i,Jstr,kout)
              ad_vbar(i,Jstr  ,know)=ad_vbar(i,Jstr,know)+              &
     &                               (1.0_r8-Ce)*adfac
              ad_vbar(i,Jstr+1,know)=ad_vbar(i,Jstr+1,know)+            &
     &                               Ce*adfac
              ad_Ce=ad_Ce-                                              &
     &              (vbar(i,Jstr  ,know)-                               &
     &               vbar(i,Jstr+1,know))*adfac
              ad_bry_val=ad_bry_val+adfac
              ad_cff1=ad_cff1-                                          &
     &                (Ze-BOUNDARY(ng)%zeta_south(i))*adfac
              ad_Ze=ad_Ze+cff1*adfac
              ad_vbar(i,Jstr,kout)=0.0_r8
              IF (Ce.gt.Co) THEN
!>              tl_Ze=tl_Ze+cff2*tl_cff3+                               &
!>   &                      tl_cff2*cff3
!>
                ad_cff2=ad_cff2+cff3*ad_Ze
                ad_cff3=ad_cff3+cff2*ad_Ze
!>              tl_cff3=tl_zeta(i,Jstr,kout)+                           &
!>   &                  Ce*tl_zeta(i,Jstr-1,know)+                      &
!>   &                  tl_Ce*(zeta(i,Jstr-1,know)+                     &
!>   &                         zeta(i,Jstr  ,know))-                    &
!>   &                  (1.0_r8+Ce)*tl_zeta(i,Jstr,know)
!>
                ad_zeta(i,Jstr  ,kout)=ad_zeta(i,Jstr  ,kout)+          &
     &                                 ad_cff3
                ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)-          &
     &                                 (1.0_r8+Ce)*ad_cff3
                ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+          &
     &                                 Ce*ad_cff3
                ad_Ce=ad_Ce+                                            &
     &                (zeta(i,Jstr-1,know)+                             &
     &                 zeta(i,Jstr  ,know))*ad_cff3
                ad_cff3=0.0_r8
!>              tl_cff2=2.0_r8*cff2*Co*tl_Ce/(Ce*Ce)
!>
                ad_Ce=ad_Ce+                                            &
     &                2.0_r8*cff2*Co*ad_cff2/(Ce*Ce)
                ad_cff2=0.0_r8
              END IF
!>            tl_Ze=(0.5_r8+Ce)*tl_zeta(i,Jstr  ,know)+                 &
!>   &              (0.5_r8-Ce)*tl_zeta(i,Jstr-1,know)+                 &
!>   &              tl_Ce*(zeta(i,Jstr  ,know)-                         &
!>   &                     zeta(i,Jstr-1,know))
!>
              ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+            &
     &                               (0.5_r8+Ce)*ad_Ze
              ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+            &
     &                               (0.5_r8-Ce)*ad_Ze
              ad_Ce=ad_Ce+                                              &
     &              (zeta(i,Jstr  ,know)-                               &
     &               zeta(i,Jstr-1,know))*ad_Ze
              ad_Ze=0.0_r8
!>            tl_Ce=dt2d*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                 &
!>   &                           GRID(ng)%pn(i,Jstr  ))*                &
!>   &              (cff1*tl_cff+                                       &
!>   &               tl_cff1*cff)
!>
              adfac=dt2d*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                 &
     &                           GRID(ng)%pn(i,Jstr  ))*ad_Ce
              ad_cff=ad_cff+cff1*adfac
              ad_cff1=ad_cff1+cff*adfac
              ad_Ce=0.0_r8
!>            tl_cff1=-0.5_r8*cff1*tl_cff/cff
!>
              ad_cff=ad_cff-                                            &
     &               0.5_r8*cff1*ad_cff/cff
              ad_cff1=0.0_r8
!>            tl_cff=0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+                   &
!>   &                       GRID(ng)%tl_h(i,Jstr  ))
!>
              adfac=0.5_r8*ad_cff
              GRID(ng)%ad_h(i,Jstr-1)=GRID(ng)%ad_h(i,Jstr-1)+adfac
              GRID(ng)%ad_h(i,Jstr  )=GRID(ng)%ad_h(i,Jstr  )+adfac
              ad_cff=0.0_r8

# if defined SSH_TIDES_NOT_YET && !defined UV_TIDES_NOT_YET
              IF (ad_LBC(isouth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(i,Jstr,know)-                          &
     &                      BOUNDARY(ng)%zeta_south(i))*                &
     &                  0.5_r8*GRID(ng)%pn(i,Jstr)
              ELSE
                bry_pgr=-g*(zeta(i,Jstr  ,know)-                        &
     &                      zeta(i,Jstr-1,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                  &
     &                          GRID(ng)%pn(i,Jstr  ))
              END IF
#  ifdef UV_COR
              bry_cor=-0.125_r8*(ubar(i  ,Jstr-1,know)+                 &
     &                           ubar(i+1,Jstr-1,know)+                 &
     &                           ubar(i  ,Jstr  ,know)+                 &
     &                           ubar(i+1,Jstr  ,know))*                &
     &                          (GRID(ng)%f(i,Jstr-1)+                  &
     &                           GRID(ng)%f(i,Jstr  ))
#  else
              bry_cor=0.0_r8
#  endif
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
!>            tl_bry_val=tl_vbar(i,Jstr+1,know)+                        &
!>   &                   tl_cff2*(bry_pgr+                              &
!>   &                            bry_cor+                              &
!>   &                            bry_str)+                             &
!>   &                   cff2*(tl_bry_pgr+                              &
!>   &                         tl_bry_cor+                              &
!>   &                         tl_bry_str)
!>
              adfac=cff2*ad_bry_val
              tl_bry_pgr=tl_bry_pgr+adfac
              tl_bry_cor=tl_bry_cor+adfac
              tl_bry_str=tl_bry_str+adfac
              ad_cff2=ad_cff2+(bry_pgr+                                 &
     &                         bry_cor+                                 &
     &                         bry_str)*ad_bry_val
              ad_vbar(i,Jstr+1,know)=ad_vbar(i,Jstr+1,know)+ad_bry_val
              ad_bry_val=0.0_r8
!>            tl_cff2=GRID(ng)%on_v(i,Jstr)*tl_Ce
!>
              ad_Ce=ad_Ce+GRID(ng)%on_v(i,Jstr)*ad_cff2
              ad_cff2=0.0_r8
!>            tl_Ce=-Ce*Ce*Ce*0.25_r8*g*(GRID(ng)%tl_h(i,Jstr-1)+       &
!>   &                                   tl_zeta(i,Jstr-1,know)+        &
!>   &                                   GRID(ng)%tl_h(i,Jstr  )+       &
!>   &                                   tl_zeta(i,Jstr  ,know))
!>
              adfac=-Ce*Ce*Ce*0.25_r8*g*ad_Ce
              ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+adfac
              ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+adfac
              GRID(ng)%ad_h(i,Jstr-1)=GRID(ng)%ad_h(i,Jstr-1)+adfac
              GRID(ng)%ad_h(i,Jstr  )=GRID(ng)%ad_h(i,Jstr  )+adfac
              ad_Ce=0.0_r8
!>            tl_bry_str=tl_cff1*(FORCES(ng)%svstr(i,Jstr)-             &
!>   &                            FORCES(ng)%bvstr(i,Jstr))+            &
!>   &                   cff1*(FORCES(ng)%tl_svstr(i,Jstr)-             &
!>   &                         FORCES(ng)%tl_bvstr(i,Jstr))
!>
              adfac=cff1*ad_bry_str
              FORCES(ng)%ad_svstr(i,Jstr)=FORCES(ng)%ad_svstr(i,Jstr)+  &
     &                                    adfac
              FORCES(ng)%ad_bvstr(i,Jstr)=FORCES(ng)%ad_bvstr(i,Jstr)-  &
     &                                    adfac
              ad_cff1=ad_cff1+(FORCES(ng)%svstr(i,Jstr)-                &
     &                         FORCES(ng)%bvstr(i,Jstr))*ad_bry_str
              ad_bry_str=0.0_r8
!>            tl_cff1=-cff1*cff1*(0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+      &
!>   &                                    tl_zeta(i,Jstr-1,know)+       &
!>   &                                    GRID(ng)%tl_h(i,Jstr  )+      &
!>   &                                    tl_zeta(i,Jstr  ,know)))
!>
              adfac=-cff1*cff1*0.5_r8*ad_cff1
              ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+adfac
              ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+adfac
              GRID(ng)%ad_h(i,Jstr-1)=GRID(ng)%ad_h(i,Jstr-1)+adfac
              GRID(ng)%ad_h(i,Jstr  )=GRID(ng)%ad_h(i,Jstr  )+adfac
              ad_cff1=0.0_r8
#  ifdef UV_COR
!>            tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jstr-1,know)+           &
!>   &                              tl_ubar(i+1,Jstr-1,know)+           &
!>   &                              tl_ubar(i  ,Jstr  ,know)+           &
!>   &                              tl_ubar(i+1,Jstr  ,know))*          &
!>   &                             (GRID(ng)%f(i,Jstr-1)+               &
!>   &                              GRID(ng)%f(i,Jstr  ))
!>
              adfac=-0.125_r8*(GRID(ng)%f(i,Jstr-1)+                    &
     &                         GRID(ng)%f(i,Jstr  ))*ad_bry_cor
              ad_ubar(i  ,Jstr-1,know)=ad_ubar(i  ,Jstr-1,know)+adfac
              ad_ubar(i+1,Jstr-1,know)=ad_ubar(i+1,Jstr-1,know)+adfac
              ad_ubar(i  ,Jstr  ,know)=ad_ubar(i  ,Jstr  ,know)+adfac
              ad_ubar(i+1,Jstr  ,know)=ad_ubar(i+1,Jstr  ,know)+adfac
              ad_bry_cor=0.0_r8
#  else
!>            tl_bry_cor=0.0_r8
!>
              ad_bry_cor=0.0_r8
#  endif
              IF (ad_LBC(isouth,isFsur,ng)%acquire) THEN
#  ifdef ADJUST_BOUNDARY
                IF (Lobc(isouth,isVbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr+                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_south(i)*           &
!>   &                       0.5_r8*GRID(ng)%pn(i,Jstr)
!>
                  BOUNDARY(ng)%ad_zeta_south(i)=BOUNDARY(ng)%           &
     &                                                ad_zeta_south(i)+ &
     &                                          g*0.5_r8*               &
     &                                          GRID(ng)%pn(i,Jstr)*    &
     &                                          ad_bry_pgr
                END IF
#  endif
!>              tl_bry_pgr=-g*tl_zeta(i,Jstr,know)*                     &
!>   &                     0.5_r8*GRID(ng)%pn(i,Jstr)
!>
                tl_zeta(i,Jstr,know)=tl_zeta(i,Jstr,know)-              &
     &                               g*0.5_r8*GRID(ng)%pn(i,Jstr)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(i,Jstr  ,know)-                  &
!>   &                         tl_zeta(i,Jstr-1,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pn(i,Jstr-1)+               &
!>   &                             GRID(ng)%pn(i,Jstr  ))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                 &
     &                           GRID(ng)%pn(i,Jstr  ))*ad_bry_pgr
                ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)-adfac
                ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
# else
#   ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
!>              tl_bry_val=BOUNDARY(ng)%tl_vbar_south(i)
!>
                BOUNDARY(ng)%ad_vbar_south(i)=BOUNDARY(ng)%             &
     &                                              ad_vbar_south(i)+   &
     &                                        ad_bry_val
                ad_bry_val=0.0_r8
              ELSE
!>              tl_bry_val=0.0_r8
!>
                ad_bry_val=0.0_r8
              END IF
#   else
!>            tl_bry_val=0.0_r8
!>
              ad_bry_val=0.0_r8
#   endif
# endif
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (ad_LBC(isouth,isVbar,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
# ifdef MASKING
!>            tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
!>   &                             GRID(ng)%vmask(i,Jstr)
!>
              ad_vbar(i,Jstr,kout)=ad_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
# ifdef ADJUST_BOUNDARY
              IF (Lobc(isouth,isVbar,ng)) THEN
!>              tl_vbar(i,Jstr,kout)=BOUNDARY(ng)%tl_vbar_south(i)
!>
                BOUNDARY(ng)%ad_vbar_south(i)=BOUNDARY(ng)%             &
     &                                              ad_vbar_south(i)+   &
     &                                        ad_vbar(i,Jstr,kout)
                ad_vbar(i,Jstr,kout)=0.0_r8
              ELSE
!>              tl_vbar(i,Jstr,kout)=0.0_r8
!>
                ad_vbar(i,Jstr,kout)=0.0_r8
              END IF
# else
!>            tl_vbar(i,Jstr,kout)=0.0_r8
!>
              ad_vbar(i,Jstr,kout)=0.0_r8
# endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (ad_LBC(isouth,isVbar,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
# ifdef MASKING
!>            tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
!>   &                             GRID(ng)%vmask(i,Jstr)
!>
              ad_vbar(i,Jstr,kout)=ad_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
!>            tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr+1,kout)
!>
              ad_vbar(i,Jstr+1,kout)=ad_vbar(i,Jstr+1,kout)+            &
     &                               ad_vbar(i,Jstr,kout)
              ad_vbar(i,Jstr  ,kout)=0.0_r8
            END IF
          END DO
!
!  Southern edge, reduced-physics boundary condition.
!
        ELSE IF (ad_LBC(isouth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                 &
     &                            zeta(i,Jstr-1,know)+                  &
     &                            GRID(ng)%h(i,Jstr  )+                 &
     &                            zeta(i,Jstr  ,know)))
# ifdef MASKING
!>            tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,kout)*                &
!>   &                             GRID(ng)%vmask(i,Jstr)
!>
              ad_vbar(i,Jstr,kout)=ad_vbar(i,Jstr,kout)*                &
     &                             GRID(ng)%vmask(i,Jstr)
# endif
!>            tl_vbar(i,Jstr,kout)=tl_vbar(i,Jstr,know)+                &
!>   &                             dt2d*(tl_bry_pgr+                    &
!>   &                                   tl_bry_cor+                    &
!>   &                                   tl_bry_str)
!>
              adfac=dt2d*ad_vbar(i,Jstr,kout)
              ad_bry_pgr=ad_bry_pgr+adfac
              ad_bry_cor=ad_bry_cor+adfac
              ad_bry_str=ad_bry_str+adfac
              ad_vbar(i,Jstr,know)=ad_vbar(i,Jstr,know)+                &
     &                             ad_vbar(i,Jstr,kout)
              ad_vbar(i,Jstr,kout)=0.0_r8
!>            tl_bry_str=tl_cff*(FORCES(ng)%svstr(i,Jstr)-              &
!>   &                           FORCES(ng)%bvstr(i,Jstr))+             &
!>   &                   cff*(FORCES(ng)%tl_svstr(i,Jstr)-              &
!>   &                        FORCES(ng)%tl_bvstr(i,Jstr))
!>
              adfac=cff*ad_bry_str
              FORCES(ng)%ad_svstr(i,Jstr)=FORCES(ng)%ad_svstr(i,Jstr)+  &
     &                                    adfac
              FORCES(ng)%ad_bvstr(i,Jstr)=FORCES(ng)%ad_bvstr(i,Jstr)-  &
     &                                    adfac
              ad_cff=ad_cff+(FORCES(ng)%svstr(i,Jstr)-                  &
     &                       FORCES(ng)%bvstr(i,Jstr))*ad_bry_str
              ad_bry_str=0.0_r8
!>            tl_cff=-cff*cff*0.5_r8*(GRID(ng)%tl_h(i,Jstr-1)+          &
!>   &                                tl_zeta(i,Jstr-1,know)+           &
!>   &                                GRID(ng)%tl_h(i,Jstr  )+          &
!>   &                                tl_zeta(i,Jstr  ,know))
!>
              adfac=-cff*cff*0.5_r8*ad_cff
              ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+adfac
              ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+adfac
              GRID(ng)%ad_h(i,Jstr-1)=GRID(ng)%ad_h(i,Jstr-1)+adfac
              GRID(ng)%ad_h(i,Jstr  )=GRID(ng)%ad_h(i,Jstr  )+adfac
              ad_cff=0.0_r8
# ifdef UV_COR
!>            tl_bry_cor=-0.125_r8*(tl_ubar(i  ,Jstr-1,know)+           &
!>   &                              tl_ubar(i+1,Jstr-1,know)+           &
!>   &                              tl_ubar(i  ,Jstr  ,know)+           &
!>   &                              tl_ubar(i+1,Jstr  ,know))*          &
!>   &                             (GRID(ng)%f(i,Jstr-1)+               &
!>   &                              GRID(ng)%f(i,Jstr  ))
!>
              adfac=-0.125_r8*(GRID(ng)%f(i,Jstr-1)+                    &
     &                         GRID(ng)%f(i,Jstr  ))*ad_bry_cor
              ad_ubar(i  ,Jstr-1,know)=ad_ubar(i  ,Jstr-1,know)+adfac
              ad_ubar(i+1,Jstr-1,know)=ad_ubar(i+1,Jstr-1,know)+adfac
              ad_ubar(i  ,Jstr  ,know)=ad_ubar(i  ,Jstr  ,know)+adfac
              ad_ubar(i+1,Jstr  ,know)=ad_ubar(i+1,Jstr  ,know)+adfac
              ad_bry_cor=0.0_r8
# else
!>            bry_cor=0.0_r8
!>
              tl_bry_cor=0.0_r8
# endif
              IF (ad_LBC(isouth,isFsur,ng)%acquire) THEN
# ifdef ADJUST_BOUNDARY
                IF (Lobc(isouth,isVbar,ng)) THEN
!>                tl_bry_pgr=tl_bry_pgr+                                &
!>   &                       g*BOUNDARY(ng)%tl_zeta_south(i)*           &
!>   &                       0.5_r8*GRID(ng)%pn(i,Jstr)
!>
                  BOUNDARY(ng)%ad_zeta_south(i)=BOUNDARY(ng)%           &
     &                                                ad_zeta_south(i)+ &
     &                                          g*0.5_r8*               &
     &                                          GRID(ng)%pn(i,Jstr)*    &
     &                                          ad_bry_pgr
                END IF
# endif
!>              tl_bry_pgr=-g*tl_zeta(i,Jstr,know)*                     &
!>   &                     0.5_r8*GRID(ng)%pn(i,Jstr)
!>
                ad_zeta(i,Jstr,know)=ad_zeta(i,Jstr,know)-              &
     &                               g*0.5_r8*GRID(ng)%pn(i,Jstr)*      &
     &                               ad_bry_pgr
                ad_bry_pgr=0.0_r8
              ELSE
!>              tl_bry_pgr=-g*(tl_zeta(i,Jstr  ,know)-                  &
!>   &                         tl_zeta(i,Jstr-1,know))*                 &
!>   &                     0.5_r8*(GRID(ng)%pn(i,Jstr-1)+               &
!>   &                             GRID(ng)%pn(i,Jstr  ))
!>
                adfac=-g*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                 &
     &                           GRID(ng)%pn(i,Jstr  ))*ad_bry_pgr
                ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)-adfac
                ad_zeta(i,Jstr  ,know)=ad_zeta(i,Jstr  ,know)+adfac
                ad_bry_pgr=0.0_r8
              END IF
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (ad_LBC(isouth,isVbar,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
!>            tl_vbar(i,Jstr,kout)=0.0_r8
!>
              ad_vbar(i,Jstr,kout)=0.0_r8
            END IF
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_v2dbc_tile
#endif
      END MODULE ad_v2dbc_mod
