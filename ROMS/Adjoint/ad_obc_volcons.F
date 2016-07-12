#include "cppdefs.h"
      MODULE ad_obc_volcons_mod
#ifdef ADJOINT
!
!svn $Id: ad_obc_volcons.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routines computes adjoint integral mass flux  "obc_flux"       !
!  across the open boundaries, which is needed to enforce global       !
!  mass conservation constraint.                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_obc_flux_tile, ad_set_DUV_bc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_obc_flux (ng, tile, kinp)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kinp
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_obc_flux_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       kinp,                                      &
# ifdef MASKING
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % ad_h,                           &
     &                       GRID(ng) % om_v,                           &
     &                       GRID(ng) % on_u,                           &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
     &                       OCEAN(ng) % ad_zeta)
      RETURN
      END SUBROUTINE ad_obc_flux
!
!***********************************************************************
      SUBROUTINE ad_obc_flux_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             kinp,                                &
# ifdef MASKING
     &                             umask, vmask,                        &
# endif
     &                             h, ad_h,                             &
     &                             om_v, on_u,                          &
     &                             ubar, vbar, zeta,                    &
     &                             ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kinp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_h(LBi:,LBj:)
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(inout) :: ad_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff, my_area, my_flux
      real(r8) :: adfac, ad_cff, ad_my_area, ad_my_flux

# ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
# endif
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8
!
!-----------------------------------------------------------------------
!  Perform adjoint global summation and compute correction velocity.
!-----------------------------------------------------------------------
!
      IF (ANY(ad_VolCons(:,ng))) THEN
!>      tl_ubar_xs=tl_bc_flux/bc_area-                                  &
!>   &             tl_bc_area*ubar_xs/bc_area
!>
        ad_bc_area=-ad_ubar_xs*ubar_xs/bc_area
        ad_bc_flux=ad_ubar_xs/bc_area
        ad_ubar_xs=0.0_r8
!>      tl_bc_flux=tl_bc_flux+tl_my_flux
!>
        ad_my_flux=ad_bc_flux
!>      tl_bc_area=tl_bc_area+tl_my_area
!>
        ad_my_area=ad_bc_area
      END IF
!
!-----------------------------------------------------------------------
!  Compute open segments cross-section area and mass flux.
!-----------------------------------------------------------------------
!
      IF (ad_VolCons(inorth,ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            cff=0.5_r8*om_v(i,Jend+1)*                                  &
     &          (zeta(i,Jend  ,kinp)+h(i,Jend  )+                       &
     &           zeta(i,Jend+1,kinp)+h(i,Jend+1))
# ifdef MASKING
            cff=cff*vmask(i,Jend+1)
# endif
!>          tl_my_flux=tl_my_flux-                                      &
!>   &                 tl_cff*vbar(i,Jend+1,kinp)-                      &
!>   &                 cff*tl_vbar(i,Jend+1,kinp)
!>
            ad_vbar(i,Jend+1,kinp)=ad_vbar(i,Jend+1,kinp)-              &
     &                             cff*ad_my_flux
            ad_cff=ad_cff-ad_my_flux*vbar(i,Jend+1,kinp)
!>          tl_my_area=tl_my_area+tl_cff
!>
            ad_cff=ad_cff+ad_my_area
# ifdef MASKING
!>          tl_cff=tl_cff*vmask(i,Jend+1)
!>
            ad_cff=ad_cff*vmask(i,Jend+1)
# endif
!>          tl_cff=0.5_r8*om_v(i,Jend+1)*                               &
!>   &             (tl_zeta(i,Jend  ,kinp)+tl_h(i,Jend  )+              &
!>   &              tl_zeta(i,Jend+1,kinp)+tl_h(i,Jend+1))
!>
            adfac=0.5_r8*om_v(i,Jend+1)*ad_cff
            ad_zeta(i,Jend  ,kinp)=ad_zeta(i,Jend  ,kinp)+adfac
            ad_zeta(i,Jend+1,kinp)=ad_zeta(i,Jend+1,kinp)+adfac
            ad_h(i,Jend  )=ad_h(i,Jend  )+adfac
            ad_h(i,Jend+1)=ad_h(i,Jend+1)+adfac
            ad_cff=0.0_r8
          END DO
        END IF
      END IF

      IF (ad_VolCons(isouth,ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            cff=0.5_r8*om_v(i,Jstr)*                                    &
     &          (zeta(i,Jstr-1,kinp)+h(i,Jstr-1)+                       &
     &           zeta(i,Jstr  ,kinp)+h(i,Jstr  ))
# ifdef MASKING
            cff=cff*vmask(i,Jstr)
# endif
!>          tl_my_flux=tl_my_flux+                                      &
!>   &                 tl_cff*vbar(i,JstrV-1,kinp)+                     &
!>   &                 cff*tl_vbar(i,JstrV-1,kinp)
!>
            ad_vbar(i,JstrV-1,kinp)=ad_vbar(i,JstrV-1,kinp)+            &
     &                              cff*ad_my_flux
            ad_cff=ad_cff+ad_my_flux*vbar(i,JstrV-1,kinp)
!>          tl_my_area=tl_my_area+tl_cff
!>
            ad_cff=ad_cff+ad_my_area
# ifdef MASKING
!>          tl_cff=tl_cff*vmask(i,Jstr)
!>
            ad_cff=ad_cff*vmask(i,Jstr)
# endif
!>          tl_cff=0.5_r8*om_v(i,Jstr)*                                 &
!>   &             (tl_zeta(i,Jstr-1,kinp)+tl_h(i,Jstr-1)+              &
!>   &              tl_zeta(i,Jstr  ,kinp)+tl_h(i,Jstr  ))
!>
            adfac=0.5_r8*om_v(i,Jstr)*ad_cff
            ad_zeta(i,Jstr-1,kinp)=ad_zeta(i,Jstr-1,kinp)+adfac
            ad_zeta(i,Jstr  ,kinp)=ad_zeta(i,Jstr  ,kinp)+adfac
            ad_h(i,Jstr-1)=ad_h(i,Jstr-1)+adfac
            ad_h(i,Jstr  )=ad_h(i,Jstr  )+adfac
            ad_cff=0.0_r8
          END DO
        END IF
      END IF

      IF (ad_VolCons(ieast,ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            cff=0.5_r8*on_u(Iend+1,j)*                                  &
     &          (zeta(Iend  ,j,kinp)+h(Iend  ,j)+                       &
     &           zeta(Iend+1,j,kinp)+h(Iend+1,j))
# ifdef MASKING
            cff=cff*umask(Iend+1,j)
# endif
!>          tl_my_flux=tl_my_flux-                                      &
!>   &                 tl_cff*ubar(Iend+1,j,kinp)-                      &
!>   &                 cff*tl_ubar(Iend+1,j,kinp)
!>
            ad_ubar(Iend+1,j,kinp)=ad_ubar(Iend+1,j,kinp)-              &
     &                             cff*ad_my_flux
            ad_cff=ad_cff-ad_my_flux*ubar(Iend+1,j,kinp)
!>          tl_my_area=tl_my_area+tl_cff
!>
            ad_cff=ad_cff+ad_my_area
# ifdef MASKING
!>          tl_cff=tl_cff*umask(Iend+1,j)
!>
            ad_cff=ad_cff*umask(Iend+1,j)
# endif
!>          tl_cff=0.5_r8*on_u(Iend+1,j)*                               &
!>   &             (tl_zeta(Iend  ,j,kinp)+tl_h(Iend  ,j)+              &
!>   &              tl_zeta(Iend+1,j,kinp)+tl_h(Iend+1,j))
!>
            adfac=0.5_r8*on_u(Iend+1,j)*ad_cff
            ad_zeta(Iend  ,j,kinp)=ad_zeta(Iend  ,j,kinp)+adfac
            ad_zeta(Iend+1,j,kinp)=ad_zeta(Iend+1,j,kinp)+adfac
            ad_h(Iend  ,j)=ad_h(Iend  ,j)+adfac
            ad_h(Iend+1,j)=ad_h(Iend+1,j)+adfac
            ad_cff=0.0_r8
          END DO
        END IF
      END IF

      IF (ad_VolCons(iwest,ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            cff=0.5_r8*on_u(Istr,j)*                                    &
     &          (zeta(Istr-1,j,kinp)+h(Istr-1,j)+                       &
     &           zeta(Istr  ,j,kinp)+h(Istr  ,j))
# ifdef MASKING
            cff=cff*umask(Istr,j)
# endif
!>          tl_my_flux=tl_my_flux+                                      &
!>   &                 tl_cff*ubar(Istr,j,kinp)+                        &
!>   &                 cff*tl_ubar(Istr,j,kinp)
!>
            ad_ubar(Istr,j,kinp)=ad_ubar(Istr,j,kinp)+                  &
     &                           cff*ad_my_flux
            ad_cff=ad_cff+ad_my_flux*ubar(Istr,j,kinp)
!>          tl_my_area=tl_my_area+tl_cff
!>
            ad_cff=ad_cff+ad_my_area
# ifdef MASKING
!>          tl_cff=tl_cff*umask(Istr,j)
!>
            ad_cff=ad_cff*umask(Istr,j)
# endif
!>          tl_cff=0.5_r8*on_u(Istr,j)*                                 &
!>   &             (tl_zeta(Istr-1,j,kinp)+tl_h(Istr-1,j)+              &
!>   &              tl_zeta(Istr  ,j,kinp)+tl_h(Istr  ,j))
!>
            adfac=0.5_r8*on_u(Istr,j)*ad_cff
            ad_zeta(Istr-1,j,kinp)=ad_zeta(Istr-1,j,kinp)+adfac
            ad_zeta(Istr  ,j,kinp)=ad_zeta(Istr  ,j,kinp)+adfac
            ad_h(Istr-1,j)=ad_h(Istr-1,j)+adfac
            ad_h(Istr-1,j)=ad_h(Istr-1,j)+adfac
            ad_cff=0.0_r8
          END DO
        END IF
      END IF

!>    tl_my_area=0.0_r8
!>    tl_my_flux=0.0_r8
!>
      ad_my_area=0.0_r8
      ad_my_flux=0.0_r8

      RETURN
      END SUBROUTINE ad_obc_flux_tile
!
!***********************************************************************
      SUBROUTINE ad_set_DUV_bc_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               kinp,                              &
# ifdef MASKING
     &                               umask, vmask,                      &
# endif
     &                               om_v, on_u,                        &
     &                               ubar, vbar,                        &
     &                               ad_ubar, ad_vbar,                  &
     &                               Drhs, Duon, Dvom,                  &
     &                               ad_Drhs, ad_Duon, ad_Dvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_parallel

# ifdef DISTRIBUTE
!
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kinp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: Drhs(IminS:,JminS:)
      real(r8), intent(in) :: Duon(IminS:,JminS:)
      real(r8), intent(in) :: Dvom(IminS:,JminS:)

      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Drhs(IminS:,JminS:)
      real(r8), intent(inout) :: ad_Duon(IminS:,JminS:)
      real(r8), intent(inout) :: ad_Dvom(IminS:,JminS:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: Drhs(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: Duon(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: Dvom(IminS:ImaxS,JminS:JmaxS)

      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_Drhs(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(inout) :: ad_Duon(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(inout) :: ad_Dvom(IminS:ImaxS,JminS:JmaxS)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j

      real(r8) :: adfac, adfac1, adfac2, adfac3
      real(r8) :: ad_my_ubar_xs
# ifdef DISTRIBUTE
      real(r8) :: buffer
      character (len=3) :: op_handle
# endif

# include "set_bounds.h"

# ifdef DISTRIBUTE
#  define I_RANGE IstrU,MIN(Iend+1,Lm(ng))
#  define J_RANGE JstrV,MIN(Jend+1,Mm(ng))
# else
#  define I_RANGE MAX(2,IstrU-1),MIN(Iend+1,Lm(ng))
#  define J_RANGE MAX(2,JstrV-1),MIN(Jend+1,Mm(ng))
# endif
!
!-----------------------------------------------------------------------
!  Set vertically integrated mass fluxes "Duon" and "Dvom" along
!  the open boundaries in such a way that the integral volume is
!  conserved.  This is done by applying "ubar_xs" correction to
!  the velocities.
!-----------------------------------------------------------------------
!
      ad_my_ubar_xs=0.0_r8
!
# ifdef DISTRIBUTE

! Do a special exchange to avoid having three ghost points for high
! order numerical stencil.
!
      IF (ad_VolCons(isouth,ng).or.ad_VolCons(inorth,ng)) THEN
!>      CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      tl_Dvom)
!>
        CALL ad_mp_exchange2d (ng, tile, iADM, 1,                       &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_Dvom)
      END IF

      IF (ad_VolCons(iwest,ng).or.ad_VolCons(ieast,ng)) THEN
!>      CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      tl_Duon)
!>
        CALL ad_mp_exchange2d (ng, tile, iADM, 1,                       &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_Duon)
      END IF
# endif

      IF (ad_VolCons(inorth,ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=-2+I_RANGE+1
# ifdef MASKING
!>          tl_Dvom(i,Jend+1)=tl_Dvom(i,Jend+1)*vmask(i,Jend+1)
!>
            ad_Dvom(i,Jend+1)=ad_Dvom(i,Jend+1)*vmask(i,Jend+1)
# endif
!>          tl_Dvom(i,Jend+1)=0.5_r8*                                   &
!>   &                        ((tl_Drhs(i,Jend+1)+tl_Drhs(i,Jend))*     &
!>   &                         (vbar(i,Jend+1,kinp)+ubar_xs)+           &
!>   &                         (Drhs(i,Jend+1)+Drhs(i,Jend))*           &
!>   &                         (tl_vbar(i,Jend+1,kinp)+tl_ubar_xs))*    &
!>   &                        om_v(i,Jend+1)
!>
            adfac=0.5_r8*om_v(i,Jend+1)*ad_Dvom(i,Jend+1)
            adfac1=adfac*(vbar(i,Jend+1,kinp)+ubar_xs)
            adfac2=adfac*(Drhs(i,Jend+1)+Drhs(i,Jend))
            ad_Drhs(i,Jend+1)=ad_Drhs(i,Jend+1)+adfac1
            ad_Drhs(i,Jend  )=ad_Drhs(i,Jend  )+adfac1
            ad_vbar(i,Jend+1,kinp)=ad_vbar(i,Jend+1,kinp)+adfac2
            ad_my_ubar_xs=ad_my_ubar_xs+adfac2
            ad_Dvom(i,Jend+1)=0.0_r8
          END DO
        END IF
      END IF

      IF (ad_VolCons(isouth,ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=-2+I_RANGE+1
# ifdef MASKING
!>          tl_Dvom(i,Jstr)=tl_Dvom(i,Jstr)*vmask(i,Jstr)
!>
            ad_Dvom(i,Jstr)=ad_Dvom(i,Jstr)*vmask(i,Jstr)
# endif
!>          tl_Dvom(i,Jstr)=0.5_r8*                                     &
!>   &                      ((tl_Drhs(i,Jstr)+tl_Drhs(i,Jstr-1))*       &
!>   &                       (vbar(i,Jstr,kinp)-ubar_xs)+               &
!>   &                       (Drhs(i,Jstr)+Drhs(i,Jstr-1))*             &
!>   &                       (tl_vbar(i,Jstr,kinp)-tl_ubar_xs))*        &
!>   &                      om_v(i,Jstr)
!>
            adfac=0.5_r8*om_v(i,Jstr)*ad_Dvom(i,Jstr)
            adfac1=adfac*(vbar(i,Jstr,kinp)-ubar_xs)
            adfac2=adfac*(Drhs(i,Jstr)+Drhs(i,Jstr-1))
            ad_Drhs(i,Jstr-1)=ad_Drhs(i,Jstr-1)+adfac1
            ad_Drhs(i,Jstr  )=ad_Drhs(i,Jstr  )+adfac1
            ad_vbar(i,Jstr,kinp)=ad_vbar(i,Jstr,kinp)+adfac2
            ad_my_ubar_xs=ad_my_ubar_xs-adfac2
            ad_Dvom(i,Jstr)=0.0_r8
          END DO
        END IF
      END IF

      IF (ad_VolCons(ieast,ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=-2+J_RANGE+1
# ifdef MASKING
!>          tl_Duon(Iend+1,j)=tl_Duon(Iend+1,j)*umask(Iend+1,j)
!>
            ad_Duon(Iend+1,j)=ad_Duon(Iend+1,j)*umask(Iend+1,j)
# endif
!>          tl_Duon(Iend+1,j)=0.5_r8*                                   &
!>   &                        ((tl_Drhs(Iend+1,j)+tl_Drhs(Iend,j))*     &
!>   &                         (ubar(Iend+1,j,kinp)+ubar_xs)+           &
!>   &                         (Drhs(Iend+1,j)+Drhs(Iend,j))*           &
!>   &                         (tl_ubar(Iend+1,j,kinp)+tl_ubar_xs))*    &
!>   &                        on_u(Iend+1,j)
!>
            adfac=0.5_r8*on_u(Iend+1,j)*ad_Duon(Iend+1,j)
            adfac1=adfac*(ubar(Iend+1,j,kinp)+ubar_xs)
            adfac2=adfac*(Drhs(Iend+1,j)+Drhs(Iend,j))
            ad_Drhs(Iend  ,j)=ad_Drhs(Iend  ,j)+adfac1
            ad_Drhs(Iend+1,j)=ad_Drhs(Iend+1,j)+adfac1
            ad_ubar(Iend+1,j,kinp)=ad_ubar(Iend+1,j,kinp)+adfac2
            ad_my_ubar_xs=ad_my_ubar_xs+adfac2
            ad_Duon(Iend+1,j)=0.0_r8
          END DO
        END IF
      END IF

      IF (ad_VolCons(iwest,ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=-2+J_RANGE+1
# ifdef MASKING
!>          tl_Duon(Istr,j)=tl_Duon(Istr,j)*umask(Istr,j)
!>
            ad_Duon(Istr,j)=ad_Duon(Istr,j)*umask(Istr,j)
# endif
!>          tl_Duon(Istr,j)=0.5_r8*                                     &
!>   &                      ((tl_Drhs(Istr,j)+tl_Drhs(Istr-1,j))*       &
!>   &                       (ubar(Istr,j,kinp)-ubar_xs)+               &
!>   &                       (Drhs(Istr,j)+Drhs(Istr-1,j))*             &
!>   &                       (tl_ubar(Istr,j,kinp)-tl_ubar_xs))*        &
!>   &                      on_u(Istr,j)
!>
            adfac=0.5_r8*on_u(Istr,j)*ad_Duon(Istr,j)
            adfac1=adfac*(ubar(Istr,j,kinp)-ubar_xs)
            adfac2=adfac*(Drhs(Istr,j)+Drhs(Istr-1,j))
            ad_Drhs(Istr-1,j)=ad_Drhs(Istr-1,j)+adfac1
            ad_Drhs(Istr  ,j)=ad_Drhs(Istr  ,j)+adfac1
            ad_ubar(Istr,j,kinp)=ad_ubar(Istr,j,kinp)+adfac2
            ad_my_ubar_xs=ad_my_ubar_xs-adfac2
            ad_Duon(Istr,j)=0.0_r8
          END DO
        END IF
      END IF

# undef I_RANGE
# undef J_RANGE
!
!-----------------------------------------------------------------------
!  Perform global summation and compute correction velocity.
!-----------------------------------------------------------------------
!
      IF (ANY(ad_VolCons(:,ng))) THEN
# ifdef DISTRIBUTE
        NSUB=1                           ! distributed-memory
# else
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
# endif
!$OMP CRITICAL (AD_OBC_VOLUME)
        IF (tile_count.eq.0) THEN
          adfac3=0.0_r8
        END IF
        adfac3=adfac3+ad_my_ubar_xs
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
# ifdef DISTRIBUTE
          buffer=adfac3
          op_handle='SUM'
          CALL mp_reduce (ng, iADM, 1, buffer, op_handle)
          adfac3=buffer
# endif
          IF (iif(ng).eq.nfast(ng)+1) THEN
            ad_ubar_xs=ad_ubar_xs+adfac3
          ELSE
            ad_ubar_xs=adfac3
          ENDIF
        END IF
!$OMP END CRITICAL (AD_OBC_VOLUME)
      END IF

      RETURN
      END SUBROUTINE ad_set_DUV_bc_tile
!
!***********************************************************************
      SUBROUTINE ad_conserve_mass_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  kinp,                           &
# ifdef MASKING
     &                                  umask, vmask,                   &
# endif
     &                                  ad_ubar, ad_vbar)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kinp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Corrects velocities across the open boundaries to enforce global
!  mass conservation constraint.
!-----------------------------------------------------------------------
!
      IF (ad_VolCons(inorth,ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
# ifdef MASKING
!>          vbar(i,Jend+1,kinp)=vbar(i,Jend+1,kinp)*                    &
!>   &                          vmask(i,Jend+1)
!>
            ad_vbar(i,Jend+1,kinp)=ad_vbar(i,Jend+1,kinp)*              &
     &                             vmask(i,Jend+1)
# endif
!>          tl_vbar(i,Jend+1,kinp)=(tl_vbar(i,Jend+1,kinp)+tl_ubar_xs)
!>
            ad_ubar_xs=ad_ubar_xs+ad_vbar(i,Jend+1,kinp)
          END DO
        END IF
      END IF

      IF (ad_VolCons(isouth,ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
# ifdef MASKING
!>          tl_vbar(i,Jstr,kinp)=tl_vbar(i,Jstr,kinp)*                  &
!>   &                           vmask(i,Jstr)
!>
            ad_vbar(i,Jstr,kinp)=ad_vbar(i,Jstr,kinp)*                  &
     &                           vmask(i,Jstr)
# endif
!>          tl_vbar(i,Jstr,kinp)=(tl_vbar(i,Jstr,kinp)-tl_ubar_xs)
!>
            ad_ubar_xs=ad_ubar_xs-ad_vbar(i,Jstr,kinp)
          END DO
        END IF
      END IF

      IF (ad_VolCons(ieast,ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
# ifdef MASKING
!>          tl_ubar(Iend+1,j,kinp)=tl_ubar(Iend+1,j,kinp)*              &
!>   &                             umask(Iend+1,j)
!>
            ad_ubar(Iend+1,j,kinp)=ad_ubar(Iend+1,j,kinp)*              &
     &                             umask(Iend+1,j)
# endif
!>          tl_ubar(Iend+1,j,kinp)=tl_ubar(Iend+1,j,kinp)+tl_ubar_xs
!>
            ad_ubar_xs=ad_ubar_xs+ad_ubar(Iend+1,j,kinp)
          END DO
        END IF
      END IF

      IF (ad_VolCons(iwest,ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
# ifdef MASKING
!>          tl_ubar(Istr,j,kinp)=tl_ubar(Istr,j,kinp)*                  &
!>   &                           umask(Istr,j)
!>
            ad_ubar(Istr,j,kinp)=ad_ubar(Istr,j,kinp)*                  &
     &                           umask(Istr,j)
# endif
!>          tl_ubar(Istr,j,kinp)=tl_ubar(Istr,j,kinp)-tl_ubar_xs
!>
            ad_ubar_xs=ad_ubar_xs-ad_ubar(Istr,j,kinp)
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_conserve_mass_tile
#endif
      END MODULE ad_obc_volcons_mod
