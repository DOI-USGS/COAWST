#undef DEBUG
#ifdef ADJOINT
      SUBROUTINE ad_step2d (ng, tile)
!
!svn $Id: ad_step2d_LF_AM3.h 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine performs a fast (predictor or corrector) time-step     !
!  for the free-surface and 2D  momentum  adjoint  equations.  The     !
!  predictor step  is  Leap-Frog  whereas  the  corrector step  is     !
!  trapezoidal, Adams-Moulton.  If applicable,  it also calculates     !
!  time filtering variables over all fast-time steps  to damp high     !
!  frequency signals in 3D applications.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef CLIMATOLOGY
      USE mod_clima
# endif
# ifdef SOLVE3D
      USE mod_coupling
# endif
# ifdef DIAGNOSTICS_UV
!!    USE mod_diags
# endif
      USE mod_forces
      USE mod_grid
# if defined UV_VIS2 || defined UV_VIS4 || defined NEARSHORE_MELLOR
      USE mod_mixing
# endif
      USE mod_ocean
# if defined UV_PSOURCE || defined Q_PSOURCE
      USE mod_sources
# endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"

# ifdef PROFILE
      CALL wclock_on (ng, iADM, 9)
# endif
      CALL ad_step2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, N(ng),                   &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs(ng), kstp(ng), knew(ng),                &
# ifdef SOLVE3D
     &                     nstp(ng), nnew(ng),                          &
# endif
# if defined UV_PSOURCE || defined Q_PSOURCE
     &                     Nsrc(ng),                                    &
     &                     SOURCES(ng) % Isrc,     SOURCES(ng) % Jsrc,  &
     &                     SOURCES(ng) % Dsrc,     SOURCES(ng) % Qbar,  &
# endif
# ifdef MASKING
     &                     GRID(ng) % pmask,       GRID(ng) % rmask,    &
     &                     GRID(ng) % umask,       GRID(ng) % vmask,    &
# endif
# ifdef WET_DRY_NOT_YET
     &                     GRID(ng) % rmask_wet,                        &
     &                     GRID(ng) % umask_wet,   GRID(ng) % vmask_wet,&
# endif
# ifdef SOLVE3D
#  ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                     GRID(ng) % ad_bed_thick,                     &
#  endif
     &                     GRID(ng) % ad_Hz,                            &
     &                     GRID(ng) % ad_z_r,      GRID(ng) % ad_z_w,   &
# endif
     &                     GRID(ng) % fomn,                             &
     &                     GRID(ng) % h,           GRID(ng) % ad_h,     &
     &                     GRID(ng) % om_u,        GRID(ng) % om_v,     &
     &                     GRID(ng) % on_u,        GRID(ng) % on_v,     &
     &                     GRID(ng) % omn,                              &
     &                     GRID(ng) % pm,          GRID(ng) % pn,       &
# if defined CURVGRID && defined UV_ADV
     &                     GRID(ng) % dndx,        GRID(ng) % dmde,     &
# endif
# if defined UV_VIS2 || defined UV_VIS4
     &                     GRID(ng) % pmon_r,      GRID(ng) % pnom_r,   &
     &                     GRID(ng) % pmon_p,      GRID(ng) % pnom_p,   &
     &                     GRID(ng) % om_r,        GRID(ng) % on_r,     &
     &                     GRID(ng) % om_p,        GRID(ng) % on_p,     &
#  ifdef UV_VIS2
     &                     MIXING(ng) % visc2_p,   MIXING(ng) % visc2_r,&
#  endif
#  ifdef UV_VIS4
     &                     MIXING(ng) % visc4_p,   MIXING(ng) % visc4_r,&
#  endif
# endif
# ifdef NEARSHORE_MELLOR
     &                     MIXING(ng) % ad_rustr2d,                     &
     &                     MIXING(ng) % ad_rvstr2d,                     &
     &                     OCEAN(ng) % ad_rulag2d,                      &
     &                     OCEAN(ng) % ad_rvlag2d,                      &
     &                     OCEAN(ng) % ubar_stokes,                     &
     &                     OCEAN(ng) % ad_ubar_stokes,                  &
     &                     OCEAN(ng) % vbar_stokes,                     &
     &                     OCEAN(ng) % ad_vbar_stokes,                  &
# endif
# ifdef M2CLIMATOLOGY
     &                     CLIMA(ng) % ubarclm,    CLIMA(ng) % vbarclm, &
#  ifdef M2CLM_NUDGING
     &                     CLIMA(ng) % M2nudgcof,                       &
#  endif
# endif
# ifndef SOLVE3D
     &                     FORCES(ng) % ad_sustr,                       &
     &                     FORCES(ng) % ad_svstr,                       &
     &                     FORCES(ng) % ad_bustr,                       &
     &                     FORCES(ng) % ad_bvstr,                       &
# else
#  ifdef VAR_RHO_2D
     &                     COUPLING(ng) % rhoA,                         &
     &                     COUPLING(ng) % ad_rhoA,                      &
     &                     COUPLING(ng) % rhoS,                         &
     &                     COUPLING(ng) % ad_rhoS,                      &
#  endif
     &                     COUPLING(ng) % ad_DU_avg1,                   &
     &                     COUPLING(ng) % ad_DU_avg2,                   &
     &                     COUPLING(ng) % ad_DV_avg1,                   &
     &                     COUPLING(ng) % ad_DV_avg2,                   &
     &                     COUPLING(ng) % Zt_avg1,                      &
     &                     COUPLING(ng) % ad_Zt_avg1,                   &
     &                     COUPLING(ng) % ad_rufrc,                     &
     &                     COUPLING(ng) % ad_rvfrc,                     &
     &                     OCEAN(ng) % ad_ru,                           &
     &                     OCEAN(ng) % ad_rv,                           &
# endif
# ifdef DIAGNOSTICS_UV
!!   &                     DIAGS(ng) % DiaU2wrk,   DIAGS(ng) % DiaV2wrk,&
!!   &                     DIAGS(ng) % DiaRUbar,   DIAGS(ng) % DiaRVbar,&
#  ifdef SOLVE3D
!!   &                     DIAGS(ng) % DiaU2int,   DIAGS(ng) % DiaV2int,&
!!   &                     DIAGS(ng) % DiaRUfrc,   DIAGS(ng) % DiaRVfrc,&
#  endif
# endif
# ifndef SOLVE3D
     &                     OCEAN(ng) % ad_ubar_sol,                     &
     &                     OCEAN(ng) % ad_vbar_sol,                     &
     &                     OCEAN(ng) % ad_zeta_sol,                     &
# endif
     &                     OCEAN(ng) % rubar,      OCEAN(ng) % ad_rubar,&
     &                     OCEAN(ng) % rvbar,      OCEAN(ng) % ad_rvbar,&
     &                     OCEAN(ng) % rzeta,      OCEAN(ng) % ad_rzeta,&
     &                     OCEAN(ng) % ubar,       OCEAN(ng) % ad_ubar, &
     &                     OCEAN(ng) % vbar,       OCEAN(ng) % ad_vbar, &
     &                     OCEAN(ng) % zeta,       OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 9)
# endif
      RETURN
      END SUBROUTINE ad_step2d
!
!***********************************************************************
      SUBROUTINE ad_step2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, UBk,               &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           krhs, kstp, knew,                      &
# ifdef SOLVE3D
     &                           nstp, nnew,                            &
# endif
# if defined UV_PSOURCE || defined Q_PSOURCE
     &                           Nsrc, Isrc, Jsrc, Dsrc, Qbar,          &
# endif
# ifdef MASKING
     &                           pmask, rmask, umask, vmask,            &
# endif
# ifdef WET_DRY_NOT_YET
     &                           rmask_wet, umask_wet, vmask_wet,       &
# endif
# ifdef SOLVE3D
#  ifdef ICESHELF
     &                           zice,                                  &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                           ad_bed_thick,                          &
#  endif
     &                           ad_Hz, ad_z_r, ad_z_w,                 &
# endif
     &                           fomn,                                  &
     &                           h, ad_h,                               &
     &                           om_u, om_v, on_u, on_v, omn, pm, pn,   &
# if defined CURVGRID && defined UV_ADV
     &                           dndx, dmde,                            &
# endif
# if defined UV_VIS2 || defined UV_VIS4
     &                           pmon_r, pnom_r, pmon_p, pnom_p,        &
     &                           om_r, on_r, om_p, on_p,                &
#  ifdef UV_VIS2
     &                           visc2_p, visc2_r,                      &
#  endif
#  ifdef UV_VIS4
     &                           visc4_p, visc4_r,                      &
#  endif
# endif
# ifdef NEARSHORE_MELLOR
     &                           ad_rustr2d, ad_rvstr2d,                &
     &                           ad_rulag2d, ad_rvlag2d,                &
     &                           ubar_stokes, ad_ubar_stokes,           &
     &                           vbar_stokes, ad_vbar_stokes,           &
# endif
# ifdef M2CLIMATOLOGY
     &                           ubarclm, vbarclm,                      &
#  ifdef M2CLM_NUDGING
     &                           M2nudgcof,                             &
#  endif
# endif
# ifndef SOLVE3D
     &                           ad_sustr, ad_svstr,                    &
     &                           ad_bustr, ad_bvstr,                    &
# else
#  ifdef VAR_RHO_2D
     &                           rhoA, ad_rhoA, rhoS, ad_rhoS,          &
#  endif
     &                           ad_DU_avg1, ad_DU_avg2,                &
     &                           ad_DV_avg1, ad_DV_avg2,                &
     &                           Zt_avg1, ad_Zt_avg1,                   &
     &                           ad_rufrc, ad_rvfrc,                    &
     &                           ad_ru, ad_rv,                          &
# endif
# ifdef DIAGNOSTICS_UV
!!   &                           DiaU2wrk, DiaV2wrk,                    &
!!   &                           DiaRUbar, DiaRVbar,                    &
#  ifdef SOLVE3D
!!   &                           DiaU2int, DiaV2int,                    &
!!   &                           DiaRUfrc, DiaRVfrc,                    &
#  endif
# endif
# ifndef SOLVE3D 
     &                           ad_ubar_sol, ad_vbar_sol,              &
     &                           ad_zeta_sol,                           &
# endif
     &                           rubar, ad_rubar,                       &
     &                           rvbar, ad_rvbar,                       &
     &                           rzeta, ad_rzeta,                       &
     &                           ubar, ad_ubar,                         &
     &                           vbar, ad_vbar,                         &
     &                           zeta, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
!
# ifdef WET_DRY_NOT_YET
      USE bc_2d_mod
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE ad_exchange_2d_mod
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
# ifdef OBC_VOLCONS
      USE obc_volcons_mod
      USE ad_obc_volcons_mod
# endif
# ifdef SOLVE3D
      USE ad_set_depth_mod, ONLY : ad_set_depth_tile
# endif
      USE ad_u2dbc_mod, ONLY : ad_u2dbc_tile
      USE ad_v2dbc_mod, ONLY : ad_v2dbc_tile
      USE ad_zetabc_mod, ONLY : ad_zetabc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: krhs, kstp, knew
# ifdef SOLVE3D
      integer, intent(in) :: nstp, nnew
# endif
!
# ifdef ASSUMED_SHAPE
#  if defined UV_PSOURCE || defined Q_PSOURCE
      integer, intent(in) :: Nsrc
      integer, intent(in) :: Isrc(:)
      integer, intent(in) :: Jsrc(:)
      real(r8), intent(in) :: Dsrc(:)
      real(r8), intent(in) :: Qbar(:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout):: ad_bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
#  endif
      real(r8), intent(in) :: fomn(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  if defined CURVGRID && defined UV_ADV
      real(r8), intent(in) :: dndx(LBi:,LBj:)
      real(r8), intent(in) :: dmde(LBi:,LBj:)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
#   ifdef UV_VIS2
      real(r8), intent(in) :: visc2_p(LBi:,LBj:)
      real(r8), intent(in) :: visc2_r(LBi:,LBj:)
#   endif
#   ifdef UV_VIS4
      real(r8), intent(in) :: visc4_p(LBi:,LBj:)
      real(r8), intent(in) :: visc4_r(LBi:,LBj:)
#   endif
#  endif
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: ubar_stokes(LBi:,LBj:)
      real(r8), intent(in) :: vbar_stokes(LBi:,LBj:)
#  endif
#  ifdef M2CLIMATOLOGY
      real(r8), intent(in) :: ubarclm(LBi:,LBj:)
      real(r8), intent(in) :: vbarclm(LBi:,LBj:)
#   ifdef M2CLM_NUDGING
      real(r8), intent(in) :: M2nudgcof(LBi:,LBj:)
#   endif
#  endif
      real(r8), intent(in) :: rubar(LBi:,LBj:,:)
      real(r8), intent(in) :: rvbar(LBi:,LBj:,:)
      real(r8), intent(in) :: rzeta(LBi:,LBj:,:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
#   if defined VAR_RHO_2D
      real(r8), intent(in) :: rhoA(LBi:,LBj:)
      real(r8), intent(in) :: rhoS(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)

      real(r8), intent(inout) :: ad_DU_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DU_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DV_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
#   if defined VAR_RHO_2D
      real(r8), intent(inout) :: ad_rhoA(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rhoS(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: ad_rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: ad_ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: ad_rv(LBi:,LBj:,0:,:)
#  else 
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_bustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_bvstr(LBi:,LBj:)
#  endif
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: ad_rustr2d(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rvstr2d(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rulag2d(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rvlag2d(LBi:,LBj:)
      real(r8), intent(inout) :: ad_ubar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: ad_vbar_stokes(LBi:,LBj:)
#  endif
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(inout) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: umask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_wet(LBi:,LBj:)
#  endif
#  ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaU2wrk(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaV2wrk(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaRUbar(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRVbar(LBi:,LBj:,:,:)
#   ifdef SOLVE3D
!!    real(r8), intent(inout) :: DiaU2int(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaV2int(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
#   endif
#  endif
      real(r8), intent(inout) :: ad_h(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_rvbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_rzeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifndef SOLVE3D
      real(r8), intent(out) :: ad_ubar_sol(LBi:,LBj:)
      real(r8), intent(out) :: ad_vbar_sol(LBi:,LBj:)
      real(r8), intent(out) :: ad_zeta_sol(LBi:,LBj:)
#  endif

# else

#  if defined UV_PSOURCE || defined Q_PSOURCE
      integer, intent(in) :: Nsrc
      integer, intent(in) :: Isrc(Nsrc)
      integer, intent(in) :: Jsrc(Nsrc)

      real(r8), intent(in) :: Dsrc(Nsrc)
      real(r8), intent(in) :: Qbar(Nsrc)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout):: ad_bed_thick(LBi:UBi,LBj:UBi,2)
#   endif
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:UBk)
#  endif
      real(r8), intent(in) :: fomn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  if defined CURVGRID && defined UV_ADV
      real(r8), intent(in) :: dndx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: dmde(LBi:UBi,LBj:UBj)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
#   ifdef UV_VIS2
      real(r8), intent(in) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc2_r(LBi:UBi,LBj:UBj)
#   endif
#   ifdef UV_VIS4
      real(r8), intent(in) :: visc4_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc4_r(LBi:UBi,LBj:UBj)
#   endif
#  endif
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbar_stokes(LBi:UBi,LBj:UBj)
#  endif
#  ifdef M2CLIMATOLOGY
      real(r8), intent(in) :: ubarclm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbarclm(LBi:UBi,LBj:UBj)
#   ifdef M2CLM_NUDGING
      real(r8), intent(in) :: M2nudgcof(LBi:UBi,LBj:UBj)
#   endif
#  endif
      real(r8), intent(in) :: rubar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: rvbar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: rzeta(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rhoS(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: ad_DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
#   if defined VAR_RHO_2D
      real(r8), intent(inout) :: ad_rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rhoS(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: ad_rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_ru(LBi:UBi,LBj:UBj,0:UBk,2)
      real(r8), intent(inout) :: ad_rv(LBi:UBi,LBj:UBj,0:UBk,2)
#  else
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bvstr(LBi:UBi,LBj:UBj)
#  endif
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: ad_rustr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rvstr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rulag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rvlag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_vbar_stokes(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(inout) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaU2wrk(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaV2wrk(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaRUbar(LBi:UBi,LBj:UBj,2,NDM2d-1)
!!    real(r8), intent(inout) :: DiaRVbar(LBi:UBi,LBj:UBj,2,NDM2d-1)
#   ifdef SOLVE3D
!!    real(r8), intent(inout) :: DiaU2int(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaV2int(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
#   endif
#  endif
      real(r8), intent(inout) :: ad_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rubar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ad_rvbar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ad_rzeta(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifndef SOLVE3D
      real(r8), intent(out) :: ad_ubar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: ad_vbar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: ad_zeta_sol(LBi:UBi,LBj:UBj)
#  endif
# endif
!
!  Local variable declarations.
!
      logical :: CORRECTOR_2D_STEP
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j, ptsk
# if defined UV_PSOURCE || defined Q_PSOURCE
      integer :: is
# endif
# ifdef DIAGNOSTICS_UV
!!    integer :: idiag
# endif

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5
      real(r8) :: fac, fac1, fac2
      real(r8) :: adfac, adfac1, adfac2, adfac3, adfac4
      real(r8) :: ad_cff, ad_cff1, ad_fac, ad_fac1

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dgrad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dnew
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Drhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Drhs_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: DUon
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: DVom
# ifdef NEARSHORE_MELLOR
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: DUSon
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: DVSom
# endif
# ifdef UV_VIS4
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LapU
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LapV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gzeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gzeta2
# if defined VAR_RHO_2D && defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gzetaSA
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rhs_ubar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rhs_vbar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rhs_zeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zeta_new
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zwrk
# ifdef WET_DRY_NOT_YET
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wetdry
# endif
# ifdef DIAGNOSTICS_UV
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,
!!   &                                     NDM2d-1) :: DiaU2rhs
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,
!!   &                                     NDM2d-1) :: DiaV2rhs
# endif

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Dgrad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Dnew
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Drhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Drhs_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_DUon
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_DVom
# ifdef NEARSHORE_MELLOR
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_DUSon
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_DVSom
# endif
# ifdef UV_VIS4
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_LapU
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_LapV
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_grad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_gzeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_gzeta2
# if defined VAR_RHO_2D && defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_gzetaSA
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_rhs_ubar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_rhs_vbar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_rhs_zeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_zeta_new
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_zwrk

# include "set_bounds.h"
!
      ptsk=3-kstp
      CORRECTOR_2D_STEP=.not.PREDICTOR_2D_STEP(ng)
# ifdef DEBUG
      WRITE (21,20) iic(ng), CORRECTOR_2D_STEP,                         &
     &              kstp, krhs, knew, ptsk
 20   FORMAT (' iic = ',i5.5,' corrector = ',l1,' kstp = ',i1,          &
     &        ' krhs = ',i1,' knew = ',i1,' ptsk = ',i1)
# endif
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8
      ad_cff1=0.0_r8
      ad_fac=0.0_r8
      ad_fac1=0.0_r8
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          ad_Dgrad(i,j)=0.0_r8
          ad_Dnew(i,j)=0.0_r8
          ad_Drhs(i,j)=0.0_r8
          ad_Drhs_p(i,j)=0.0_r8
          ad_Dstp(i,j)=0.0_r8
          ad_DUon(i,j)=0.0_r8
          ad_DVom(i,j)=0.0_r8
# ifdef NEARSHORE_MELLOR
          ad_DUSon(i,j)=0.0_r8
          ad_DVSom(i,j)=0.0_r8
# endif
# ifdef UV_VIS4
          ad_LapU(i,j)=0.0_r8
          ad_LapV(i,j)=0.0_r8
# endif
          ad_UFe(i,j)=0.0_r8
          ad_UFx(i,j)=0.0_r8
          ad_VFe(i,j)=0.0_r8
          ad_VFx(i,j)=0.0_r8
          ad_grad(i,j)=0.0_r8
          ad_gzeta(i,j)=0.0_r8 
          ad_gzeta2(i,j)=0.0_r8
# if defined VAR_RHO_2D && defined SOLVE3D
          ad_gzetaSA(i,j)=0.0_r8
# endif
          ad_rhs_ubar(i,j)=0.0_r8
          ad_rhs_vbar(i,j)=0.0_r8
          ad_rhs_zeta(i,j)=0.0_r8
          ad_rhs_zeta(i,j)=0.0_r8
          ad_zeta_new(i,j)=0.0_r8
          ad_zwrk(i,j)=0.0_r8
          ad_DUon(i,j)=0.0_r8
          ad_DVom(i,j)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE total depth (m) arrays and vertically
!  integerated mass fluxes.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend+1
#  else
#   define I_RANGE IstrU,MIN(Iend+1,Lm(ng))
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend+1
#  else
#   define J_RANGE JstrV,MIN(Jend+1,Mm(ng))
#  endif
# else
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU-1,Iend+1
#  else
#   define I_RANGE MAX(2,IstrU-1),MIN(Iend+1,Lm(ng))
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV-1,Jend+1
#  else
#   define J_RANGE MAX(2,JstrV-1),MIN(Jend+1,Mm(ng))
#  endif
# endif
      DO j=-2+J_RANGE+1
        DO i=-2+I_RANGE+1
          Dnew(i,j)=zeta(i,j,knew)+h(i,j)
          Drhs(i,j)=zeta(i,j,krhs)+h(i,j)
          Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
        END DO
      END DO
      DO j=-2+J_RANGE+1
        DO i=-1+I_RANGE+1
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
          DUon(i,j)=ubar(i,j,krhs)*cff1
# ifdef NEARSHORE_MELLOR
          DUSon(i,j)=ubar_stokes(i,j)*cff1
          DUon(i,j)=DUon(i,j)+DUSon(i,j)
# endif
        END DO
      END DO
      DO j=-1+J_RANGE+1
        DO i=-2+I_RANGE+1
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
          DVom(i,j)=vbar(i,j,krhs)*cff1
# ifdef NEARSHORE_MELLOR
          DVSom(i,j)=vbar_stokes(i,j)*cff1
          DVom(i,j)=DVom(i,j)+DVSom(i,j)
# endif
        END DO
      END DO
# ifdef DISTRIBUTE
!
!  Do a special exchange to avoid having three ghost points for
!  high order numerical stencil. Notice that a private array is
!  passed to the exchange routine.  It will also apply periodic
!  boundary conditions if no partitions in I- or J-directions.
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        DUon)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        DVom)
#  endif
      CALL mp_exchange2d (ng, tile, iADM, 2,                            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    DUon, DVom)
# endif
# undef I_RANGE
# undef J_RANGE
# ifdef OBC_VOLCONS
!
!  Compute integral mass flux across open boundaries and adjust
!  for volume conservation. Compute BASIC STATE value.
!  This must be computed here instead of below.
!
        CALL obc_flux_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      knew,                                       &
#  ifdef MASKING
     &                      umask, vmask,                               &
#  endif
     &                      h, om_v, on_u,                              &
     &                      ubar, vbar, zeta)
!
!  Set vertically integrated mass fluxes DUon and DVom along the open
!  boundaries in such a way that the integral volume is conserved.
!
      CALL set_DUV_bc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs,                                       &
#  ifdef MASKING
     &                      umask, vmask,                               &
#  endif
     &                      om_v, on_u,                                 &
     &                      ubar, vbar,                                 &
     &                      Drhs, DUon, DVom)
# endif
# ifdef UV_VIS4
!
!  Compute BASIC state depths at PSI-points for viscosity.
!
#  ifdef EW_PERIODIC
#   define IV_RANGE Istr-1,Iend+1
#   define IU_RANGE Istr-1,Iend+1
#  else
#   define IV_RANGE MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
#   define IU_RANGE MAX(2,IstrU-1),MIN(Iend+1,Lm(ng))
#  endif
#  ifdef NS_PERIODIC
#   define JU_RANGE Jstr-1,Jend+1
#   define JV_RANGE Jstr-1,Jend+1
#  else
#   define JU_RANGE MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
#   define JV_RANGE MAX(2,JstrV-1),MIN(Jend+1,Mm(ng))
#  endif
!
      DO j=JU_RANGE+1
        DO i=IV_RANGE+1
# else
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
# endif
          Drhs_p(i,j)=0.25_r8*(Drhs(i,j  )+Drhs(i-1,j  )+               &
     &                         Drhs(i,j-1)+Drhs(i-1,j-1))
        END DO
      END DO
# undef IU_RANGE
# undef IV_RANGE
# undef JU_RANGE
# undef JV_RANGE
!!
!! Since the BASIC STATE is not recomputed, set right-hand-side
!! terms.
!!
!!    DO j=Jstr,Jend
!!        DO i=IstrU,Iend
!!          rhs_ubar(i,j)=rubar(i,j,1)
!!        END DO
!!    END DO
!!    DO j=JstrV,Jend
!!      DO i=Istr,Iend
!!        rhs_vbar(i,j)=rvbar(i,j,1)
!!        END DO
!!      END DO
!
!  Initialize BASIC STATE right-hand-side terms.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          rhs_ubar(i,j)=0.0_r8
          rhs_vbar(i,j)=0.0_r8
        END DO
      END DO
!
!  Do not perform the actual time stepping during the auxiliary
!  (nfast(ng)+1) time step.
!
      STEP_LOOP : IF (iif(ng).le.nfast(ng)) THEN
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_ubar(:,:,knew),                            &
!>   &                    tl_vbar(:,:,knew))
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_ubar(:,:,knew),                         &
     &                       ad_vbar(:,:,knew))
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_vbar(:,:,knew))
!>
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,knew))
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_ubar(:,:,knew))
!>
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,knew))
#  endif
# endif
# ifdef UV_PSOURCE
!
!-----------------------------------------------------------------------
!  Apply adjoint of mass point sources.
!-----------------------------------------------------------------------
!
        DO is=1,Nsrc
          i=Isrc(is)
          j=Jsrc(is)
          IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
            IF (INT(Dsrc(is)).eq.0) THEN
              cff=1.0_r8/(on_u(i,j)*0.5_r8*(Dnew(i-1,j)+Dnew(i,j)))
#  ifdef SOLVE3D
!>            tl_DU_avg1(i,j)=0.0_r8
!>
              ad_DU_avg1(i,j)=0.0_r8
#  endif
!>            tl_ubar(i,j,knew)=Qbar(is)*tl_cff
!>
              ad_cff=ad_cff+Qbar(is)*ad_ubar(i,j,knew)
              ad_ubar(i,j,knew)=0.0_r8
!>            tl_cff=-cff*cff*on_u(i,j)*                                &
!>   &               0.5_r8*(tl_Dnew(i-1,j)+tl_Dnew(i,)
!>
              adfac=-cff*cff*on_u(i,j)*0.5_r8*ad_cff
              ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac
              ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac
              ad_cff=0.0_r8
            ELSE
              cff=1.0_r8/(om_v(i,j)*0.5_r8*(Dnew(i,j-1)+Dnew(i,j)))
#  ifdef SOLVE3D
!>            tl_DV_avg1(i,j)=0.08
!>
              ad_DV_avg1(i,j)=0.08
#  endif
!>            tl_vbar(i,j,knew)=Qbar(is)*tl_cff
!>
              ad_cff=ad_cff+Qbar(is)*ad_vbar(i,j,knew)
              ad_vbar(i,j,knew)=0.0_r8
!>            tl_cff=-cff*cff*om_v(i,j)*                                &
!>   &               0.5_r8*(tl_Dnew(i,j-1)+tl_Dnew(i,j))
!>
              adfac=-cff*cff*om_v(i,j)*0.5_r8*ad_cff
              ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac
              ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac
              ad_cff=0.0_r8
            END IF
          END IF
        END DO
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
!>          tl_Dnew(i,j)=tl_zeta(i,j,knew)+tl_h(i,j)
!>
            ad_zeta(i,j,knew)=ad_zeta(i,j,knew)+ad_Dnew(i,j)
            ad_h(i,j)=ad_h(i,j)+ad_Dnew(i,j)
            ad_Dnew(i,j)=0.0_r8
          END DO
        END DO
# endif
# ifdef OBC_VOLCONS
!>      CALL tl_obc_flux_tile (ng, tile,                                &
!>   &                         LBi, UBi, LBj, UBj,                      &
!>   &                         IminS, ImaxS, JminS, JmaxS,              &
!>   &                         knew,                                    &
#  ifdef MASKING
!>   &                         umask, vmask,                            &
#  endif
!>   &                         h, tl_h, om_v, on_u,                     &
!>   &                         ubar, vbar, zeta,                        &
!>   &                         tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_obc_flux_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         knew,                                    &
#  ifdef MASKING
     &                         umask, vmask,                            &
#  endif
     &                         h, ad_h, om_v, on_u,                     &
     &                         ubar, vbar, zeta,                        &
     &                         ad_ubar, ad_vbar, ad_zeta)
# endif
!
!-----------------------------------------------------------------------
!  Apply adjoint lateral boundary conditions.
!-----------------------------------------------------------------------
!
!>      CALL tl_v2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, knew,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
!>      CALL tl_u2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, knew,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
!
!  If predictor step, load right-side-term into shared arrays for
!  future use during the subsequent corrector step.
!
        IF (PREDICTOR_2D_STEP(ng)) THEN
# ifdef DIAGNOSTICS_UV
!!        DO idiag=1,NDM2d-1
!!          DO j=Jstr,Jend
!!            DO i=IstrU,Iend
!!              DiaRUbar(i,j,krhs,idiag)=DiaU2rhs(i,j,idiag)
!!            END DO
!!          END DO
!!          DO j=JstrV,Jend
!!            DO i=Istr,Iend
!!              DiaRVbar(i,j,krhs,idiag)=DiaV2rhs(i,j,idiag)
!!            END DO
!!          END DO
!!        END DO
# endif
          DO j=JstrV,Jend
            DO i=Istr,Iend
!>            tl_rvbar(i,j,krhs)=tl_rhs_vbar(i,j)
!>
              ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)+ad_rvbar(i,j,krhs)
              ad_rvbar(i,j,krhs)=0.0_r8
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
!>            tl_rubar(i,j,krhs)=tl_rhs_ubar(i,j)
!>
              ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)+ad_rubar(i,j,krhs)
              ad_rubar(i,j,krhs)=0.0_r8
            END DO
          END DO
        END IF
# ifdef DIAGNOSTICS_UV
!!
!!-----------------------------------------------------------------------
!!  Time step 2D momentum diagnostic terms.
!!-----------------------------------------------------------------------
!!
#  ifdef SOLVE3D
!!
!!  The arrays "DiaU2rhs" and "DiaV2rhs" contain the contributions of
!!  each of the 2D right-hand-side terms for the momentum equations.
!!
!!  These values are integrated, time-stepped and converted to mass flux
!!  units (m3 s-1) for coupling with the 3D diagnostic terms.
!!
!!      fac=weight(1,iif(ng),ng)
!!      IF (FIRST_2D_STEP.and.CORRECTOR_2D_STEP) THEN
!!        cff1=0.5_r8*dtfast(ng)
!!        DO idiag=1,NDM2d-1
!!          DO j=JstrV,Jend
!!            DO i=Istr,Iend
!!              DiaV2int(i,j,idiag)=cff1*DiaV2rhs(i,j,idiag)
!!              DiaV2wrk(i,j,idiag)=DiaV2int(i,j,idiag)*                &
!!   &                              (pn(i,j)+pn(i,j-1))*fac
!!            END DO
!!          END DO
!!          DO j=Jstr,Jend
!!            DO i=IstrU,Iend
!!              DiaU2int(i,j,idiag)=cff1*DiaU2rhs(i,j,idiag)
!!              DiaU2wrk(i,j,idiag)=DiaU2int(i,j,idiag)*                &
!!   &                              (pm(i-1,j)+pm(i,j))*fac
!!            END DO
!!          END DO
!!        END DO
!!      ELSE IF (CORRECTOR_2D_STEP) THEN
!!        cff1=0.5_r8*dtfast(ng)*5.0_r8/12.0_r8
!!        cff2=0.5_r8*dtfast(ng)*8.0_r8/12.0_r8
!!        cff3=0.5_r8*dtfast(ng)*1.0_r8/12.0_r8
!!        DO idiag=1,NDM2d-1
!!          DO j=JstrV,Jend
!!            DO i=Istr,Iend
!!              DiaV2int(i,j,idiag)=DiaV2int(i,j,idiag)+                &
!!   &                              (cff1*DiaV2rhs(i,j,idiag)+          &
!!   &                               cff2*DiaRVbar(i,j,kstp,idiag)-     &
!!   &                               cff3*DiaRVbar(i,j,ptsk,idiag))
!!              DiaV2wrk(i,j,idiag)=DiaV2wrk(i,j,idiag)+                &
!!   &                              DiaV2int(i,j,idiag)*                &
!!   &                              (pn(i,j)+pn(i,j-1))*fac
!!            END DO
!!          END DO
!!          DO j=Jstr,Jend
!!            DO i=IstrU,Iend
!!              DiaU2int(i,j,idiag)=DiaU2int(i,j,idiag)+                &
!!   &                              (cff1*DiaU2rhs(i,j,idiag)+          &
!!   &                               cff2*DiaRUbar(i,j,kstp,idiag)-     &
!!   &                               cff3*DiaRUbar(i,j,ptsk,idiag))
!!              DiaU2wrk(i,j,idiag)=DiaU2wrk(i,j,idiag)+                &
!!   &                              DiaU2int(i,j,idiag)*                &
!!   &                              (pm(i-1,j)+pm(i,j))*fac
!!            END DO
!!          END DO
!!        END DO
!!      END IF
#  else
!!
!!  Time-step the diagnostic terms.
!!
!!      IF (FIRST_2D_STEP.and.CORRECTOR_2D_STEP) THEN
!!        cff1=0.5_r8*dtfast(ng)
!!        DO j=JstrV,Jend
!!          DO i=Istr,Iend
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!            DiaV2wrk(i,j,M2rate)=vbar(i,j,knew)-vbar(i,j,kstp)*       &
!!   &                             (Dstp(i,j)+Dstp(i,j-1))*fac
!!          END DO
!!        END DO
!!        DO j=Jstr,Jend
!!          DO i=IstrU,Iend
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!            DiaU2wrk(i,j,M2rate)=ubar(i,j,knew)-ubar(i,j,kstp)*       &
!!   &                             (Dstp(i,j)+Dstp(i-1,j))*fac
!!          END DO
!!        END DO
!!        DO idiag=1,NDM2d-1
!!          DO j=JstrV,Jend
!!            DO i=Istr,Iend
!!              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!!              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!              DiaV2wrk(i,j,idiag)=cff*cff1*DiaV2rhs(i,j,idiag)*fac
!!            END DO
!!          END DO
!!          DO j=Jstr,Jend
!!            DO i=IstrU,Iend
!!              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
!!              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!              DiaU2wrk(i,j,idiag)=cff*cff1*DiaU2rhs(i,j,idiag)*fac
!!            END DO
!!          END DO
!!        END DO
!!      ELSE IF (CORRECTOR_2D_STEP) THEN
!!        cff1=0.5_r8*dtfast(ng)*5.0_r8/12.0_r8
!!        cff2=0.5_r8*dtfast(ng)*8.0_r8/12.0_r8
!!        cff3=0.5_r8*dtfast(ng)*1.0_r8/12.0_r8
!!        DO j=JstrV,Jend
!!          DO i=Istr,Iend
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!            DiaV2wrk(i,j,M2rate)=vbar(i,j,knew)-                      &
!!   &                             vbar(i,j,kstp)*                      &
!!   &                             (Dstp(i,j)+Dstp(i,j-1))*fac
!!          END DO
!!        END DO
!!        DO j=Jstr,Jend
!!          DO i=IstrU,Iend
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!            DiaU2wrk(i,j,M2rate)=ubar(i,j,knew)-                      &
!!   &                             ubar(i,j,kstp)*                      &
!!   &                             (Dstp(i,j)+Dstp(i-1,j))*fac
!!          END DO
!!        END DO
!!        DO idiag=1,NDM2d-1
!!          DO j=JstrV,Jend
!!            DO i=Istr,Iend
!!              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!!              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!              DiaV2wrk(i,j,idiag)=cff*(cff1*DiaV2rhs(i,j,idiag)+      &
!!   &                                   cff2*DiaRVbar(i,j,kstp,idiag)- &
!!   &                                   cff3*DiaRVbar(i,j,ptsk,idiag))*&
!!   &                                  fac
!!            END DO
!!          END DO
!!          DO j=Jstr,Jend
!!            DO i=IstrU,Iend
!!              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
!!              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!              DiaU2wrk(i,j,idiag)=cff*(cff1*DiaU2rhs(i,j,idiag)+      &
!!   &                                   cff2*DiaRUbar(i,j,kstp,idiag)- &
!!   &                                   cff3*DiaRUbar(i,j,ptsk,idiag))*&
!!   &                                  fac
!!            END DO
!!          END DO
!!        END DO
!!      END IF
#  endif
# endif
# if defined WET_DRY_NOT_YET
!
!-----------------------------------------------------------------------
! If wet/drying, compute new masks for cells with depth < Dcrit.
!-----------------------------------------------------------------------
!
! HGA:  This option is not adjointed. The code below is for the NLM.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          wetdry(i,j)=1.0_r8
          IF (zwrk(i,j).le.(Dcrit(ng)-h(i,j))) THEN
            wetdry(i,j)=0.0_r8
          END IF
#  ifdef MASKING
          wetdry(i,j)=wetdry(i,j)*rmask(i,j)
#  endif
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          rmask_wet(i,j)=wetdry(i,j)
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          umask_wet(i,j)=wetdry(i-1,j)+wetdry(i,j)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          vmask_wet(i,j)=wetdry(i,j-1)+wetdry(i,j)
        END DO
      END DO
!
!  Apply boundary conditions
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rmask_wet)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  umask_wet)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  vmask_wet)

#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rmask_wet)
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        umask_wet)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vmask_wet)
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rmask_wet, umask_wet, vmask_wet)
#  endif
# endif
!
!=======================================================================
!  Time step adjoint 2D momentum equations.
!=======================================================================

# ifndef SOLVE3D
!
!  Save 2D momentum adjoint solution for IO purposes.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_ubar_sol(i,j)=ad_ubar(i,j,knew)
          END DO
          IF (j.ge.Jstr) THEN
            DO i=IstrR,IendR
              ad_vbar_sol(i,j)=ad_vbar(i,j,knew)
            END DO
          END IF
        END DO
# endif
!
!  During the first time-step, the predictor step is Forward-Euler
!  and the corrector step is Backward-Euler. Otherwise, the predictor
!  step is Leap-frog and the corrector step is Adams-Moulton.
# ifdef WET_DRY_NOT_YET
!  HGA:  This option is not fully adjointed yet.  We need to resolve
!        the issued time-dependent wet/dry mask arrays.
# endif
!
        IF (FIRST_2D_STEP) THEN
          cff1=0.5_r8*dtfast(ng)
# ifdef WET_DRY_NOT_YET
          cff2=1.0_r8/cff1
# endif
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
# ifdef WET_DRY_NOT_YET
              fac1=cff2/cff
!>            tl_rhs_vbar(i,j)=(tl_vbar(i,j,knew)*                      &
!>   &                          (Dnew(i,j)+Dnew(i,j-1))+                &
!>   &                          vbar(i,j,knew)*                         &
!>   &                          (tl_Dnew(i,j)+tl_Dnew(i,j-1))-          &
!>   &                          tl_vbar(i,j,kstp)*                      &
!>   &                          (Dstp(i,j)+Dstp(i,j-1))-                &
!>   &                          vbar(i,j,kstp)*                         &
!>   &                          (tl_Dstp(i,j)+tl_Dstp(i,j-1)))*fac1
!>
              adfac=fac1*ad_rhs_vbar(i,j)
              adfac1=adfac*vbar(i,j,knew)
              adfac2=adfac*vbar(i,j,kstp)
              ad_vbar(i,j,knew)=ad_vbar(i,j,knew)+                      &
     &                          (Dnew(i,j)+Dnew(i,j-1))*adfac
              ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)-                      &
     &                          (Dstp(i,j)+Dstp(i,j-1))*adfac
              ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac1
              ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac1
              ad_Dstp(i,j-1)=ad_Dstp(i,j-1)-adfac2
              ad_Dstp(i,j  )=ad_Dstp(i,j  )-adfac2
              ad_rhs_vbar(i,j)=0.0_r8
              IF (vmask_wet(i,j).eq.1.0_r8) THEN
                IF (rmask_wet(i,j-1).eq.1.0_r8) THEN
!>                tl_vbar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8, vbar(i,j,knew)))*     &
!>   &                              tl_vbar(i,j,knew)
!>
                  ad_vbar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8, vbar(i,j,knew)))*     &
     &                              ad_vbar(i,j,knew)
                ELSE
!>                tl_vbar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8,-vbar(i,j,knew)))*     &
!>   &                              tl_vbar(i,j,knew)
!>
                  ad_vbar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8,-vbar(i,j,knew)))*     &
     &                              ad_vbar(i,j,knew)
                END IF
              ELSE
!>              tl_vbar(i,j,knew)=0.5_r8*tl_vbar(i,j,knew)*             &
!>   &                            vmask_wet(i,j)
                ad_vbar(i,j,knew)=0.5_r8*ad_vbar(i,j,knew)*             &
     &                            vmask_wet(i,j)
              END IF
# endif
# ifdef MASKING
!>            tl_vbar(i,j,knew)=tl_vbar(i,j,knew)*vmask(i,j)
!>
              ad_vbar(i,j,knew)=ad_vbar(i,j,knew)*vmask(i,j)
# endif
!>            tl_vbar(i,j,knew)=(tl_vbar(i,j,kstp)*                     &
!>   &                           (Dstp(i,j)+Dstp(i,j-1))+               &
!>   &                           vbar(i,j,kstp)*                        &
!>   &                           (tl_Dstp(i,j)+tl_Dstp(i,j-1))+         &
!>   &                           cff*cff1*tl_rhs_vbar(i,j))*fac+        &
!>   &                          (vbar(i,j,kstp)*                        &
!>   &                           (Dstp(i,j)+Dstp(i,j-1))+               &
!>   &                           cff*cff1*rhs_vbar(i,j))*tl_fac
!>
              adfac=fac*ad_vbar(i,j,knew)
              adfac1=adfac*(Dstp(i,j)+Dstp(i,j-1))
              adfac2=adfac*cff*cff1
              adfac3=adfac*vbar(i,j,kstp)
              ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)+adfac1
              ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)+adfac2
              ad_Dstp(i,j-1)=ad_Dstp(i,j-1)+adfac3
              ad_Dstp(i,j  )=ad_Dstp(i,j  )+adfac3
              ad_fac=ad_fac+                                            &
     &               (vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+           &
     &                cff*cff1*rhs_vbar(i,j))*ad_vbar(i,j,knew)
              ad_vbar(i,j,knew)=0.0_r8
!>            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i,j-1))
!>
              adfac=-fac*fac*ad_fac
              ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac
              ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac
              ad_fac=0.0_r8
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
# ifdef WET_DRY_NOT_YET
              fac1=cff2/cff
!>            tl_rhs_ubar(i,j)=(tl_ubar(i,j,knew)*                      &
!>   &                          (Dnew(i,j)+Dnew(i-1,j))+                &
!>   &                          ubar(i,j,knew)*                         &
!>   &                          (tl_Dnew(i,j)+tl_Dnew(i-1,j))-          &
!>   &                          tl_ubar(i,j,kstp)*                      &
!>   &                          (Dstp(i,j)+Dstp(i-1,j))-                &
!>   &                          ubar(i,j,kstp)*                         &
!>   &                          (tl_Dstp(i,j)+tl_Dstp(i-1,j)))*fac1
!>
              adfac=fac1*ad_rhs_ubar(i,j)
              adfac1=adfac*ubar(i,j,knew)
              adfac2=adfac*ubar(i,j,kstp)
              ad_ubar(i,j,knew)=ad_ubar(i,j,knew)+                      &
     &                          (Dnew(i,j)+Dnew(i-1,j))*adfac
              ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)-                      &
     &                          (Dstp(i,j)+Dstp(i-1,j))*adfac
              ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac1
              ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac1
              ad_Dstp(i-1,j)=ad_Dstp(i-1,j)-adfac2
              ad_Dstp(i  ,j)=ad_Dstp(i  ,j)-adfac2
              ad_rhs_ubar(i,j)=0.0_r8
              IF (umask_wet(i,j).eq.1.0_r8) THEN
                IF (rmask_wet(i-1,j).eq.1.0_r8) THEN
!>                tl_ubar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8, ubar(i,j,knew)))*     &
!>   &                              tl_ubar(i,j,knew)
!>
                  ad_ubar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8, ubar(i,j,knew)))*     &
     &                              ad_ubar(i,j,knew)
                ELSE
!>                tl_ubar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8,-ubar(i,j,knew)))*     &
!>   &                              tl_ubar(i,j,knew)
!>
                  ad_ubar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8,-ubar(i,j,knew)))*     &
     &                              ad_ubar(i,j,knew)

                END IF
              ELSE
!>              tl_ubar(i,j,knew)=0.5_r8*tl_ubar(i,j,knew)*             &
!>   &                            umask_wet(i,j)
!>
                ad_ubar(i,j,knew)=0.5_r8*ad_ubar(i,j,knew)*             &
     &                            umask_wet(i,j)
              END IF
# endif
# ifdef MASKING
!>            tl_ubar(i,j,knew)=tl_ubar(i,j,knew)*umask(i,j)
!>
              ad_ubar(i,j,knew)=ad_ubar(i,j,knew)*umask(i,j)
# endif
!>            tl_ubar(i,j,knew)=(tl_ubar(i,j,kstp)*                     &
!>   &                           (Dstp(i,j)+Dstp(i-1,j))+               &
!>   &                           ubar(i,j,kstp)*                        &
!>   &                           (tl_Dstp(i,j)+tl_Dstp(i-1,j))+         &
!>   &                           cff*cff1*tl_rhs_ubar(i,j))*fac+        &
!>   &                          (ubar(i,j,kstp)*                        &
!>   &                           (Dstp(i,j)+Dstp(i-1,j))+               &
!>   &                           cff*cff1*rhs_ubar(i,j))*tl_fac
!>
              adfac=fac*ad_ubar(i,j,knew)
              adfac1=adfac*(Dstp(i,j)+Dstp(i-1,j))
              adfac2=adfac*cff*cff1
              adfac3=adfac*ubar(i,j,kstp)
              ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)+adfac1
              ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)+adfac2
              ad_Dstp(i-1,j)=ad_Dstp(i-1,j)+adfac3
              ad_Dstp(i  ,j)=ad_Dstp(i  ,j)+adfac3
              ad_fac=ad_fac+                                            &
     &               (ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+           &
     &                cff*cff1*rhs_ubar(i,j))*ad_ubar(i,j,knew)
              ad_ubar(i,j,knew)=0.0_r8
!>            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i-1,j))
!>
              adfac=-fac*fac*ad_fac
              ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac
              ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac
              ad_fac=0.0_r8
            END DO
          END DO
        ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
          cff1=dtfast(ng)
# ifdef WET_DRY_NOT_YET
          cff2=1.0_r8/cff1
# endif
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
# ifdef WET_DRY_NOT_YET
              fac1=cff2/cff
!>            tl_rhs_vbar(i,j)=(tl_vbar(i,j,knew)*                      &
!>   &                          (Dnew(i,j)+Dnew(i,j-1))+                &
!>   &                          vbar(i,j,knew)*                         &
!>   &                          (tl_Dnew(i,j)+tl_Dnew(i,j-1))-          &
!>   &                          tl_vbar(i,j,kstp)*                      &
!>   &                          (Dstp(i,j)+Dstp(i,j-1))-                &
!>   &                          vbar(i,j,kstp)*                         &
!>   &                          (tl_Dstp(i,j)+tl_Dstp(i,j-1)))*fac1
!>
              adfac=fac1*ad_rhs_vbar(i,j)
              adfac1=adfac*vbar(i,j,knew)
              adfac2=adfac*vbar(i,j,kstp)
              ad_vbar(i,j,knew)=ad_vbar(i,j,knew)+                      &
     &                          (Dnew(i,j)+Dnew(i,j-1))*adfac
              ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)-                      &
     &                          (Dstp(i,j)+Dstp(i,j-1))*adfac
              ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac1
              ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac1
              ad_Dstp(i,j-1)=ad_Dstp(i,j-1)-adfac2
              ad_Dstp(i,j  )=ad_Dstp(i,j  )-adfac2
              ad_rhs_vbar(i,j)=0.0_r8
              IF (vmask_wet(i,j).eq.1.0_r8) THEN
                IF (rmask_wet(i,j-1).eq.1.0_r8) THEN
!>                tl_vbar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8, vbar(i,j,knew)))*     &
!>   &                              tl_vbar(i,j,knew)
!>
                  ad_vbar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8, vbar(i,j,knew)))*     &
     &                              ad_vbar(i,j,knew)
                ELSE
!>                tl_vbar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8,-vbar(i,j,knew)))*     &
!>   &                              tl_vbar(i,j,knew)
!>
                  ad_vbar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8,-vbar(i,j,knew)))*     &
     &                              ad_vbar(i,j,knew)
                END IF
              ELSE
!>              tl_vbar(i,j,knew)=0.5_r8*tl_vbar(i,j,knew)*             &
!>   &                            vmask_wet(i,j)
                ad_vbar(i,j,knew)=0.5_r8*ad_vbar(i,j,knew)*             &
     &                            vmask_wet(i,j)
              END IF
# endif
# ifdef MASKING
!>            tl_vbar(i,j,knew)=tl_vbar(i,j,knew)*vmask(i,j)
!>
              ad_vbar(i,j,knew)=ad_vbar(i,j,knew)*vmask(i,j)
# endif
!>            tl_vbar(i,j,knew)=(tl_vbar(i,j,kstp)*                     &
!>   &                           (Dstp(i,j)+Dstp(i,j-1))+               &
!>   &                           vbar(i,j,kstp)*                        &
!>   &                           (tl_Dstp(i,j)+tl_Dstp(i,j-1))+         &
!>   &                           cff*cff1*tl_rhs_vbar(i,j))*fac+        &
!>   &                          (vbar(i,j,kstp)*                        &
!>   &                           (Dstp(i,j)+Dstp(i,j-1))+               &
!>   &                           cff*cff1*rhs_vbar(i,j))*tl_fac
!>
              adfac=fac*ad_vbar(i,j,knew)
              adfac1=adfac*(Dstp(i,j)+Dstp(i,j-1))
              adfac2=adfac*cff*cff1
              adfac3=adfac*vbar(i,j,kstp)
              ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)+adfac1
              ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)+adfac2
              ad_Dstp(i,j-1)=ad_Dstp(i,j-1)+adfac3
              ad_Dstp(i,j  )=ad_Dstp(i,j  )+adfac3
              ad_fac=ad_fac+                                            &
     &               (vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+           &
     &                cff*cff1*rhs_vbar(i,j))*ad_vbar(i,j,knew)
              ad_vbar(i,j,knew)=0.0_r8
!>            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i,j-1))
!>
              adfac=-fac*fac*ad_fac              
              ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac
              ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac
              ad_fac=0.0_r8
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
# ifdef WET_DRY_NOT_YET
              fac1=cff2/cff
!>            tl_rhs_ubar(i,j)=(tl_ubar(i,j,knew)*                      &
!>   &                          (Dnew(i,j)+Dnew(i-1,j))+                &
!>   &                          ubar(i,j,knew)*                         &
!>   &                          (tl_Dnew(i,j)+tl_Dnew(i-1,j))-          &
!>   &                          tl_ubar(i,j,kstp)*                      &
!>   &                          (Dstp(i,j)+Dstp(i-1,j))-                &
!>   &                          ubar(i,j,kstp)*                         &
!>   &                          (tl_Dstp(i,j)+tl_Dstp(i-1,j)))*fac1
!>
              adfac=fac1*ad_rhs_ubar(i,j)
              adfac1=adfac*ubar(i,j,knew)
              adfac2=adfac*ubar(i,j,kstp)
              ad_ubar(i,j,knew)=ad_ubar(i,j,knew)+                      &
     &                          (Dnew(i,j)+Dnew(i-1,j))*adfac
              ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)-                      &
                                (Dstp(i,j)+Dstp(i-1,j))*adfac
              ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac1
              ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac1
              ad_Dstp(i-1,j)=ad_Dstp(i-1,j)-adfac2
              ad_Dstp(i  ,j)=ad_Dstp(i  ,j)-adfac2
              ad_rhs_ubar(i,j)=0.0_r8
              IF (umask_wet(i,j).eq.1.0_r8) THEN
                IF (rmask_wet(i-1,j).eq.1.0_r8) THEN
!>                tl_ubar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8, ubar(i,j,knew)))*     &
!>   &                              tl_ubar(i,j,knew)
!>
                  ad_ubar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8, ubar(i,j,knew)))*     &
     &                              ad_ubar(i,j,knew)
                ELSE
!>                tl_ubar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8,-ubar(i,j,knew)))*     &
!>   &                              tl_ubar(i,j,knew)
!>
                  ad_ubar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8,-ubar(i,j,knew)))*     &
     &                              ad_ubar(i,j,knew)
                END IF
              ELSE
!>              tl_ubar(i,j,knew)=0.5_r8*tl_ubar(i,j,knew)*             &
!>   &                            umask_wet(i,j)
!>
                ad_ubar(i,j,knew)=0.5_r8*ad_ubar(i,j,knew)*             &
     &                            umask_wet(i,j)
              END IF
# endif
# ifdef MASKING
!>            tl_ubar(i,j,knew)=tl_ubar(i,j,knew)*umask(i,j)
!>
              ad_ubar(i,j,knew)=ad_ubar(i,j,knew)*umask(i,j)
# endif
!>            tl_ubar(i,j,knew)=(tl_ubar(i,j,kstp)*                     &
!>   &                           (Dstp(i,j)+Dstp(i-1,j))+               &
!>   &                           ubar(i,j,kstp)*                        &
!>   &                           (tl_Dstp(i,j)+tl_Dstp(i-1,j))+         &
!>   &                           cff*cff1*tl_rhs_ubar(i,j))*fac+        &
!>   &                          (ubar(i,j,kstp)*                        &
!>   &                           (Dstp(i,j)+Dstp(i-1,j))+               &
!>   &                           cff*cff1*rhs_ubar(i,j))*tl_fac
!>
              adfac=fac*ad_ubar(i,j,knew)
              adfac1=adfac*(Dstp(i,j)+Dstp(i-1,j))
              adfac2=adfac*cff*cff1
              adfac3=adfac*ubar(i,j,kstp)
              ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)+adfac1
              ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)+adfac2
              ad_Dstp(i-1,j)=ad_Dstp(i-1,j)+adfac3
              ad_Dstp(i  ,j)=ad_Dstp(i  ,j)+adfac3
              ad_fac=ad_fac+                                            &
     &               (ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+           &
     &                cff*cff1*rhs_ubar(i,j))*ad_ubar(i,j,knew)
              ad_ubar(i,j,knew)=0.0_r8
!>            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i-1,j))
!>
              adfac=-fac*fac*ad_fac
              ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac
              ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac
              ad_fac=0.0_r8
            END DO
          END DO
        ELSE IF (CORRECTOR_2D_STEP) THEN
          cff1=0.5_r8*dtfast(ng)*5.0_r8/12.0_r8
          cff2=0.5_r8*dtfast(ng)*8.0_r8/12.0_r8
          cff3=0.5_r8*dtfast(ng)*1.0_r8/12.0_r8
# ifdef WET_DRY_NOT_YET
          cff4=1.0_r8/cff1
# endif
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
# ifdef WET_DRY_NOT_YET
              fac1=1.0_r8/cff
!>            tl_rhs_vbar(i,j)=((tl_vbar(i,j,knew)*                     &
!>   &                           (Dnew(i,j)+Dnew(i,j-1))+               &
!>   &                           vbar(i,j,knew)*                        &
!>   &                           (tl_Dnew(i,j)+tl_Dnew(i,j-1))-         &
!>   &                           tl_vbar(i,j,kstp)*                     &
!>   &                           (Dstp(i,j)+Dstp(i,j-1))-               &
!>   &                           vbar(i,j,kstp)*                        &
!>   &                           (tl_Dstp(i,j)+tl_Dstp(i,j-1)))*fac1-   &
!>   &                          cff2*tl_rvbar(i,j,kstp)+                &
!>   &                          cff3*tl_rvbar(i,j,ptsk))*cff4
!>
              adfac=cff4*ad_rhs_vbar(i,j)
              adfac1=adfac*fac1*vbar(i,j,knew)
              adfac2=adfac*fac1*vbar(i,j,kstp)
              ad_vbar(i,j,knew)=ad_vbar(i,j,knew)+                      &
     &                          (Dnew(i,j)+Dnew(i,j-1))*adfac
              ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)-                      &
     &                          (Dstp(i,j)+Dstp(i,j-1))*adfac
              ad_rvbar(i,j,kstp)=ad_rvbar(i,j,kstp)-cff2*adfac
              ad_rvbar(i,j,ptsk)=ad_rvbar(i,j,ptsk)+cff3*adfac
              ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac1
              ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac1
              ad_Dstp(i,j-1)=ad_Dstp(i,j-1)-adfac2
              ad_Dstp(i,j  )=ad_Dstp(i,j  )-adfac2
              ad_rhs_vbar(i,j)=0.0_r8              
              IF (vmask_wet(i,j).eq.1.0_r8) THEN
                IF (rmask_wet(i,j-1).eq.1.0_r8) THEN
!>                tl_vbar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8, vbar(i,j,knew)))*     &
!>   &                              tl_vbar(i,j,knew)
!>
                  ad_vbar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8, vbar(i,j,knew)))*     &
     &                              ad_vbar(i,j,knew)
                ELSE
!>                tl_vbar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8,-vbar(i,j,knew)))*     &
!>   &                              tl_vbar(i,j,knew)
!>
                  ad_vbar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8,-vbar(i,j,knew)))*     &
     &                              ad_vbar(i,j,knew)
                END IF
              ELSE
!>              tl_vbar(i,j,knew)=0.5_r8*tl_vbar(i,j,knew)*             &
!>   &                            vmask_wet(i,j)
!>
                ad_vbar(i,j,knew)=0.5_r8*ad_vbar(i,j,knew)*             &
     &                            vmask_wet(i,j)
              END IF
# endif
# ifdef MASKING
!>            tl_vbar(i,j,knew)=tl_vbar(i,j,knew)*vmask(i,j)
!>
              ad_vbar(i,j,knew)=ad_vbar(i,j,knew)*vmask(i,j)
# endif
!>            tl_vbar(i,j,knew)=(tl_vbar(i,j,kstp)*                     &
!>   &                           (Dstp(i,j)+Dstp(i,j-1))+               &
!>   &                           vbar(i,j,kstp)*                        &
!>   &                           (tl_Dstp(i,j)+tl_Dstp(i,j-1))+         &
!>   &                           cff*(cff1*tl_rhs_vbar(i,j)+            &
!>   &                                cff2*tl_rvbar(i,j,kstp)-          &
!>   &                                cff3*tl_rvbar(i,j,ptsk)))*fac+    &
!>   &                          (vbar(i,j,kstp)*                        &
!>   &                           (Dstp(i,j)+Dstp(i,j-1))+               &
!>   &                           cff*(cff1*rhs_vbar(i,j)+               &
!>   &                                cff2*rvbar(i,j,kstp)-             &
!>   &                                cff3*rvbar(i,j,ptsk)))*tl_fac
!>
              adfac=fac*ad_vbar(i,j,knew)
              adfac1=adfac*(Dstp(i,j)+Dstp(i,j-1))
              adfac2=adfac*cff
              adfac3=adfac*vbar(i,j,kstp)
              ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)+adfac1
              ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)+cff1*adfac2
              ad_rvbar(i,j,kstp)=ad_rvbar(i,j,kstp)+cff2*adfac2
              ad_rvbar(i,j,ptsk)=-cff3*adfac2
              ad_Dstp(i,j-1)=ad_Dstp(i,j-1)+adfac3
              ad_Dstp(i,j  )=ad_Dstp(i,j  )+adfac3
              ad_fac=ad_fac+                                            &
     &               (vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+           &
     &                cff*(cff1*rhs_vbar(i,j)+                          &
     &                     cff2*rvbar(i,j,kstp)-                        &
     &                     cff3*rvbar(i,j,ptsk)))*ad_vbar(i,j,knew)
              ad_vbar(i,j,knew)=0.0_r8
!>            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i,j-1))
!>
              adfac=-fac*fac*ad_fac
              ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac
              ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac
              ad_fac=0.0_r8
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
# ifdef WET_DRY_NOT_YET
              fac1=1.0_r8/cff
!>            tl_rhs_ubar(i,j)=((tl_ubar(i,j,knew)*                     &
!>   &                           (Dnew(i,j)+Dnew(i-1,j))+               &
!>   &                           ubar(i,j,knew)*                        &
!>   &                           (tl_Dnew(i,j)+tl_Dnew(i-1,j))-         &
!>   &                           tl_ubar(i,j,kstp)*                     &
!>   &                           (Dstp(i,j)+Dstp(i-1,j))-               &
!>   &                           ubar(i,j,kstp)*                        &
!>   &                           (tl_Dstp(i,j)+tl_Dstp(i-1,j)))*fac1-   &
!>   &                          cff2*tl_rubar(i,j,kstp)+                &
!>   &                          cff3*tl_rubar(i,j,ptsk))*cff4
!>
              adfac=cff4*ad_rhs_ubar(i,j)
              adfac1=adfac*fac1*ubar(i,j,knew)
              adfac2=adfac*fac1*ubar(i,j,kstp)
              ad_ubar(i,j,knew)=ad_ubar(i,j,knew)+                      &
     &                          (Dnew(i,j)+Dnew(i-1,j))*adfac
              ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)-                      &
     &                          (Dstp(i,j)+Dstp(i-1,j))*adfac
              ad_rubar(i,j,kstp)=ad_rubar(i,j,kstp)-cff2*adfac
              ad_rubar(i,j,ptsk)=ad_rubar(i,j,ptsk)+cff3*adfac
              ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac1
              ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac1
              ad_Dstp(i-1,j)=ad_Dstp(i-1,j)-adfac2
              ad_Dstp(i  ,j)=ad_Dstp(i  ,j)-adfac2
              ad_rhs_ubar(i,j)=0.0_r8
              IF (umask_wet(i,j).eq.1.0_r8) THEN
                IF (rmask_wet(i-1,j).eq.1.0_r8) THEN
!>                tl_ubar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8, ubar(i,j,knew)))*     &
!>   &                              tl_ubar(i,j,knew)
!>
                  ad_ubar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8, ubar(i,j,knew)))*     &
     &                              ad_ubar(i,j,knew)
                ELSE
!>                tl_ubar(i,j,knew)=(0.5_r8+                            &
!>   &                               SIGN(0.5_r8,-ubar(i,j,knew)))*     &
!>   &                              tl_ubar(i,j,knew)
!>
                  ad_ubar(i,j,knew)=(0.5_r8+                            &
     &                               SIGN(0.5_r8,-ubar(i,j,knew)))*     &
     &                              ad_ubar(i,j,knew)
                END IF
              ELSE
!>              tl_ubar(i,j,knew)=0.5_r8*tl_ubar(i,j,knew)*             &
!>   &                            umask_wet(i,j)
!>
                ad_ubar(i,j,knew)=0.5_r8*ad_ubar(i,j,knew)*             &
     &                            umask_wet(i,j)
              END IF
# endif
# ifdef MASKING
!>            tl_ubar(i,j,knew)=tl_ubar(i,j,knew)*umask(i,j)
!>
              ad_ubar(i,j,knew)=ad_ubar(i,j,knew)*umask(i,j)
# endif
!>            tl_ubar(i,j,knew)=(tl_ubar(i,j,kstp)*                     &
!>   &                           (Dstp(i,j)+Dstp(i-1,j))+               &
!>   &                           ubar(i,j,kstp)*                        &
!>   &                           (tl_Dstp(i,j)+tl_Dstp(i-1,j))+         &
!>   &                           cff*(cff1*tl_rhs_ubar(i,j)+            &
!>   &                                cff2*tl_rubar(i,j,kstp)-          &
!>   &                                cff3*tl_rubar(i,j,ptsk)))*fac+    &
!>   &                          (ubar(i,j,kstp)*                        &
!>   &                           (Dstp(i,j)+Dstp(i-1,j))+               &
!>   &                           cff*(cff1*rhs_ubar(i,j)+               &
!>   &                                cff2*rubar(i,j,kstp)-             &
!>   &                                cff3*rubar(i,j,ptsk)))*tl_fac
!>
              adfac=fac*ad_ubar(i,j,knew)
              adfac1=adfac*(Dstp(i,j)+Dstp(i-1,j))
              adfac2=adfac*cff
              adfac3=adfac*ubar(i,j,kstp)
              ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)+adfac1
              ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)+cff1*adfac2
              ad_rubar(i,j,kstp)=ad_rubar(i,j,kstp)+cff2*adfac2
              ad_rubar(i,j,ptsk)=-cff3*adfac2
              ad_Dstp(i-1,j)=ad_Dstp(i-1,j)+adfac3
              ad_Dstp(i  ,j)=ad_Dstp(i  ,j)+adfac3
              ad_fac=ad_fac+                                            &
     &               (ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+           &
     &                cff*(cff1*rhs_ubar(i,j)+                          &
     &                     cff2*rubar(i,j,kstp)-                        &
     &                     cff3*rubar(i,j,ptsk)))*ad_ubar(i,j,knew)
              ad_ubar(i,j,knew)=0.0_r8
!>            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i-1,j))
!>
              adfac=-fac*fac*ad_fac
              ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac
              ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac
              ad_fac=0.0_r8
            END DO
          END DO
        END IF
!
!  Compute adjoint total water column depth.
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
!>          tl_Dstp(i,j)=tl_zeta(i,j,kstp)+tl_h(i,j)
!>
            ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_Dstp(i,j)
            ad_h(i,j)=ad_h(i,j)+ad_Dstp(i,j)
            ad_Dstp(i,j)=0.0_r8
          END DO
        END DO
!
!=======================================================================
!  Compute right-hand-side for the 2D momentum equations.
!=======================================================================

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint Coupling between 2D and 3D equations.
!-----------------------------------------------------------------------
!
!  Before the predictor step of the first barotropic time-step,
!  arrays "rufrc" and "rvfrc" contain the vertical integrals of
!  the 3D right-hand-side terms for momentum equations (including
!  surface and bottom stresses, if so prescribed).
!
!  Convert them into forcing terms by subtracting the fast time
!  "rhs_ubar" and "rhs_vbar" from them; Also, immediately apply
!  these forcing terms "rhs_ubar" and "rhs_vbar".
!
!  From now on, these newly computed forcing terms will remain
!  constant during the fast time stepping and will added to
!  "rhs_ubar" and "rhs_vbar" during all subsequent time steps.
!
        IF (FIRST_2D_STEP.and.PREDICTOR_2D_STEP(ng)) THEN
          IF (iic(ng).eq.ntfirst(ng)) THEN
            DO j=JstrV,Jend
              DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef NEARSHORE_MELLOR
!!            DiaRVfrc(i,j,nstp,M2hrad)=DiaRVfrc(i,j,3,M2hrad)
!!            DiaV2rhs(i,j,M2hrad)=DiaV2rhs(i,j,M2hrad)+                &
!!   &                             DiaRVfrc(i,j,3,M2hrad)
!!            DiaRVfrc(i,j,3,M2hrad)=DiaRVfrc(i,j,3,M2hrad)-            &
!!   &                               DiaV2rhs(i,j,M2hrad)
#   endif
#   ifdef UV_ADV
!!              DiaRVfrc(i,j,nstp,M2hadv)=DiaRVfrc(i,j,3,M2hadv)
!!              DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)+              &
!!   &                               DiaRVfrc(i,j,3,M2hadv)
!!              DiaRVfrc(i,j,3,M2hadv)=DiaRVfrc(i,j,3,M2hadv)-          &
!!   &                                 DiaV2rhs(i,j,M2hadv)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!              DiaRVfrc(i,j,nstp,M2hvis)=DiaRVfrc(i,j,3,M2hvis)
!!              DiaV2rhs(i,j,M2hvis)=DiaV2rhs(i,j,M2hvis)+              &
!!   &                               DiaRVfrc(i,j,3,M2hvis)
!!              DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)-          &
!!   &                                 DiaV2rhs(i,j,M2hvis)
#   endif
#   ifdef UV_COR
!!              DiaRVfrc(i,j,nstp,M2fcor)=DiaRVfrc(i,j,3,M2fcor)
!!              DiaV2rhs(i,j,M2fcor)=DiaV2rhs(i,j,M2fcor)+              &
!!   &                               DiaRVfrc(i,j,3,M2fcor)
!!              DiaRVfrc(i,j,3,M2fcor)=DiaRVfrc(i,j,3,M2fcor)-          &
!!   &                                 DiaV2rhs(i,j,M2fcor)
#   endif
!!              DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!              DiaV2rhs(i,j,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!              DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!              DiaV2rhs(i,j,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!              DiaRVfrc(i,j,nstp,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)
!!              DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+              &
!!   &                               DiaRVfrc(i,j,3,M2pgrd)
!!              DiaRVfrc(i,j,3,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)-          &
!!   &                                 DiaV2rhs(i,j,M2pgrd)
#  endif
!>              tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
!>
                ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rv(i,j,0,nstp)
                ad_rv(i,j,0,nstp)=0.0_r8
!>              tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
!>
                ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rhs_vbar(i,j)
!>              tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
!>
                ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)-ad_rvfrc(i,j)
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef NEARSHORE_MELLOR
!!            DiaRUfrc(i,j,nstp,M2hrad)=DiaRUfrc(i,j,3,M2hrad)
!!            DiaU2rhs(i,j,M2hrad)=DiaU2rhs(i,j,M2hrad)+                &
!!   &                             DiaRUfrc(i,j,3,M2hrad)
!!            DiaRUfrc(i,j,3,M2hrad)=DiaRUfrc(i,j,3,M2hrad)-            &
!!   &                               DiaU2rhs(i,j,M2hrad)
#   endif
#   ifdef UV_ADV
!!              DiaRUfrc(i,j,nstp,M2hadv)=DiaRUfrc(i,j,3,M2hadv)
!!              DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+              &
!!   &                               DiaRUfrc(i,j,3,M2hadv)
!!              DiaRUfrc(i,j,3,M2hadv)=DiaRUfrc(i,j,3,M2hadv)-          &
!!   &                                 DiaU2rhs(i,j,M2hadv)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!              DiaRUfrc(i,j,nstp,M2hvis)=DiaRUfrc(i,j,3,M2hvis)
!!              DiaU2rhs(i,j,M2hvis)=DiaU2rhs(i,j,M2hvis)+              &
!!   &                               DiaRUfrc(i,j,3,M2hvis)
!!              DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)-          &
!!   &                                 DiaU2rhs(i,j,M2hvis)
#   endif
#   ifdef UV_COR
!!              DiaRUfrc(i,j,nstp,M2fcor)=DiaRUfrc(i,j,3,M2fcor)
!!              DiaU2rhs(i,j,M2fcor)=DiaU2rhs(i,j,M2fcor)+              &
!!   &                               DiaRUfrc(i,j,3,M2fcor)
!!              DiaRUfrc(i,j,3,M2fcor)=DiaRUfrc(i,j,3,M2fcor)-          &
!!   &                                 DiaU2rhs(i,j,M2fcor)
#   endif
!!              DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!              DiaU2rhs(i,j,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!              DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!              DiaU2rhs(i,j,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!              DiaRUfrc(i,j,nstp,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)
!!              DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+              &
!!   &                               DiaRUfrc(i,j,3,M2pgrd)
!!              DiaRUfrc(i,j,3,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)-          &
!!   &                                 DiaU2rhs(i,j,M2pgrd)
#  endif
!>              tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
!>
                ad_rufrc(i,j)=ad_rufrc(i,j)+ad_ru(i,j,0,nstp)
                ad_ru(i,j,0,nstp)=0.0_r8
!>              tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_rufrc(i,j)
!>
                ad_rufrc(i,j)=ad_rufrc(i,j)+ad_rhs_ubar(i,j)
!>              tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
!>
                ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)-ad_rufrc(i,j)
              END DO
            END DO
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
            DO j=JstrV,Jend
              DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef NEARSHORE_MELLOR
!!            DiaRVfrc(i,j,nstp,M2hrad)=DiaRVfrc(i,j,3,M2hrad)
!!            DiaV2rhs(i,j,M2hrad)=DiaV2rhs(i,j,M2hrad)+                &
!!   &                             1.5_r8*DiaRVfrc(i,j,3,M2hrad)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2hrad)
!!            DiaRVfrc(i,j,3,M2hrad)=DiaRVfrc(i,j,3,M2hrad)-            &
!!   &                               DiaV2rhs(i,j,M2hrad)
#   endif
#   ifdef UV_ADV
!!              DiaRVfrc(i,j,nstp,M2hadv)=DiaRVfrc(i,j,3,M2hadv)
!!              DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)+              &
!!   &                               1.5_r8*DiaRVfrc(i,j,3,M2hadv)-     &
!!   &                               0.5_r8*DiaRVfrc(i,j,nnew,M2hadv)
!!              DiaRVfrc(i,j,3,M2hadv)=DiaRVfrc(i,j,3,M2hadv)-          &
!!   &                                 DiaV2rhs(i,j,M2hadv)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!              DiaRVfrc(i,j,nstp,M2hvis)=DiaRVfrc(i,j,3,M2hvis)
!!              DiaV2rhs(i,j,M2hvis)=DiaV2rhs(i,j,M2hvis)+              &
!!   &                               1.5_r8*DiaRVfrc(i,j,3,M2hvis)-     &
!!   &                               0.5_r8*DiaRVfrc(i,j,nnew,M2hvis)
!!              DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)-          &
!!   &                                 DiaV2rhs(i,j,M2hvis)
#   endif
#   ifdef UV_COR
!!              DiaRVfrc(i,j,nstp,M2fcor)=DiaRVfrc(i,j,3,M2fcor)
!!              DiaV2rhs(i,j,M2fcor)=DiaV2rhs(i,j,M2fcor)+              &
!!   &                               1.5_r8*DiaRVfrc(i,j,3,M2fcor)-     &
!!   &                               0.5_r8*DiaRVfrc(i,j,nnew,M2fcor)
!!              DiaRVfrc(i,j,3,M2fcor)=DiaRVfrc(i,j,3,M2fcor)-          &
!!   &                                 DiaV2rhs(i,j,M2fcor)
#   endif
!!              DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!              DiaV2rhs(i,j,M2sstr)=1.5_r8*DiaRVfrc(i,j,3,M2sstr)-     &
!!   &                               0.5_r8*DiaRVfrc(i,j,nnew,M2sstr)
!!              DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!              DiaV2rhs(i,j,M2bstr)=1.5_r8*DiaRVfrc(i,j,3,M2bstr)-     &
!!   &                               0.5_r8*DiaRVfrc(i,j,nnew,M2bstr)
!!              DiaRVfrc(i,j,nstp,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)
!!              DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+              &
!!   &                               1.5_r8*DiaRVfrc(i,j,3,M2pgrd)-     &
!!   &                               0.5_r8*DiaRVfrc(i,j,nnew,M2pgrd)
!!              DiaRVfrc(i,j,3,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)-          &
!!   &                                 DiaV2rhs(i,j,M2pgrd)
#  endif
!>              tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
!>
                ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rv(i,j,0,nstp)
                ad_rv(i,j,0,nstp)=0.0_r8
!>              tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                      &
!>   &                           1.5_r8*tl_rvfrc(i,j)-                  &
!>   &                           0.5_r8*tl_rv(i,j,0,nnew)
!>
                ad_rvfrc(i,j)=ad_rvfrc(i,j)+1.5_r8*ad_rhs_vbar(i,j)
                ad_rv(i,j,0,nnew)=ad_rv(i,j,0,nnew)-                    &
     &                            0.5_r8*ad_rhs_vbar(i,j)
!>              tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
!>
                ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)-ad_rvfrc(i,j)
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef NEARSHORE_MELLOR
!!            DiaRUfrc(i,j,nstp,M2hrad)=DiaRUfrc(i,j,3,M2hrad)
!!            DiaU2rhs(i,j,M2hrad)=DiaU2rhs(i,j,M2hrad)+                &
!!   &                             1.5_r8*DiaRUfrc(i,j,3,M2hrad)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2hrad)
!!            DiaRUfrc(i,j,3,M2hrad)=DiaRUfrc(i,j,3,M2hrad)-            &
!!   &                               DiaU2rhs(i,j,M2hrad)
#   endif
#   ifdef UV_ADV
!!              DiaRUfrc(i,j,nstp,M2hadv)=DiaRUfrc(i,j,3,M2hadv)
!!              DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+              &
!!   &                               1.5_r8*DiaRUfrc(i,j,3,M2hadv)-     &
!!   &                               0.5_r8*DiaRUfrc(i,j,nnew,M2hadv)
!!              DiaRUfrc(i,j,3,M2hadv)=DiaRUfrc(i,j,3,M2hadv)-          &
!!   &                                 DiaU2rhs(i,j,M2hadv)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!              DiaRUfrc(i,j,nstp,M2hvis)=DiaRUfrc(i,j,3,M2hvis)
!!              DiaU2rhs(i,j,M2hvis)=DiaU2rhs(i,j,M2hvis)+              &
!!   &                               1.5_r8*DiaRUfrc(i,j,3,M2hvis)-     &
!!   &                               0.5_r8*DiaRUfrc(i,j,nnew,M2hvis)
!!              DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)-          &
!!   &                                 DiaU2rhs(i,j,M2hvis)
#   endif
#   ifdef UV_COR
!!              DiaRUfrc(i,j,nstp,M2fcor)=DiaRUfrc(i,j,3,M2fcor)
!!              DiaU2rhs(i,j,M2fcor)=DiaU2rhs(i,j,M2fcor)+              &
!!   &                               1.5_r8*DiaRUfrc(i,j,3,M2fcor)-     &
!!   &                               0.5_r8*DiaRUfrc(i,j,nnew,M2fcor)
!!              DiaRUfrc(i,j,3,M2fcor)=DiaRUfrc(i,j,3,M2fcor)-          &
!!   &                                 DiaU2rhs(i,j,M2fcor)
#   endif
!!              DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!              DiaU2rhs(i,j,M2sstr)=1.5_r8*DiaRUfrc(i,j,3,M2sstr)-     &
!!   &                               0.5_r8*DiaRUfrc(i,j,nnew,M2sstr)
!!              DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!              DiaU2rhs(i,j,M2bstr)=1.5_r8*DiaRUfrc(i,j,3,M2bstr)-     &
!!   &                               0.5_r8*DiaRUfrc(i,j,nnew,M2bstr)
!!              DiaRUfrc(i,j,nstp,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)
!!              DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+              &
!!   &                               1.5_r8*DiaRUfrc(i,j,3,M2pgrd)-     &
!!   &                               0.5_r8*DiaRUfrc(i,j,nnew,M2pgrd)
!!              DiaRUfrc(i,j,3,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)-          &
!!   &                                 DiaU2rhs(i,j,M2pgrd)
#  endif
!>              tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
!>
                ad_rufrc(i,j)=ad_rufrc(i,j)+ad_ru(i,j,0,nstp)
                ad_ru(i,j,0,nstp)=0.0_r8
!>              tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                      &
!>   &                           1.5_r8*tl_rufrc(i,j)-                  &
!>   &                           0.5_r8*tl_ru(i,j,0,nnew)
!>
                ad_rufrc(i,j)=ad_rufrc(i,j)+1.5_r8*ad_rhs_ubar(i,j)
                ad_ru(i,j,0,nnew)=ad_ru(i,j,0,nnew)-                    &
     &                            0.5_r8*ad_rhs_ubar(i,j)
!>              tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
!>
                ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)-ad_rufrc(i,j)
              END DO
            END DO
          ELSE
            cff1=23.0_r8/12.0_r8
            cff2=16.0_r8/12.0_r8
            cff3= 5.0_r8/12.0_r8
            DO j=JstrV,Jend
              DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef NEARSHORE_MELLOR
!!            DiaRVfrc(i,j,nstp,M2hrad)=DiaRVfrc(i,j,3,M2hrad)
!!            DiaV2rhs(i,j,M2hrad)=DiaV2rhs(i,j,M2hrad)+                &
!!   &                             cff1*DiaRVfrc(i,j,3,M2hrad)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2hrad)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2hrad)
!!            DiaRVfrc(i,j,3,M2hrad)=DiaRVfrc(i,j,3,M2hrad)-            &
!!   &                               DiaV2rhs(i,j,M2hrad)
#   endif
#   ifdef UV_ADV
!!              DiaRVfrc(i,j,nstp,M2hadv)=DiaRVfrc(i,j,3,M2hadv)
!!              DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)+              &
!!   &                               cff1*DiaRVfrc(i,j,3,M2hadv)-       &
!!   &                               cff2*DiaRVfrc(i,j,nnew,M2hadv)+    &
!!   &                               cff3*DiaRVfrc(i,j,nstp,M2hadv)
!!              DiaRVfrc(i,j,3,M2hadv)=DiaRVfrc(i,j,3,M2hadv)-          &
!!   &                                 DiaV2rhs(i,j,M2hadv)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!              DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)-          &
!!   &                                 DiaV2rhs(i,j,M2hvis)
!!              DiaV2rhs(i,j,M2hvis)=DiaV2rhs(i,j,M2hvis)+              &
!!   &                               cff1*DiaRVfrc(i,j,3,M2hvis)-       &
!!   &                               cff2*DiaRVfrc(i,j,nnew,M2hvis)+    &
!!   &                               cff3*DiaRVfrc(i,j,nstp,M2hvis)
!!              DiaRVfrc(i,j,nstp,M2hvis)=DiaRVfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_COR
!!              DiaRVfrc(i,j,nstp,M2fcor)=DiaRVfrc(i,j,3,M2fcor)
!!              DiaV2rhs(i,j,M2fcor)=DiaV2rhs(i,j,M2fcor)+              &
!!   &                               cff1*DiaRVfrc(i,j,3,M2fcor)-       &
!!   &                               cff2*DiaRVfrc(i,j,nnew,M2fcor)+    &
!!   &                               cff3*DiaRVfrc(i,j,nstp,M2fcor)
!!              DiaRVfrc(i,j,3,M2fcor)=DiaRVfrc(i,j,3,M2fcor)-          &
!!   &                                 DiaV2rhs(i,j,M2fcor)
#   endif
!!              DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!              DiaV2rhs(i,j,M2sstr)=cff1*DiaRVfrc(i,j,3,M2sstr)-       &
!!   &                               cff2*DiaRVfrc(i,j,nnew,M2sstr)+    &
!!   &                               cff3*DiaRVfrc(i,j,nstp,M2sstr)
!!              DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!              DiaV2rhs(i,j,M2bstr)=cff1*DiaRVfrc(i,j,3,M2bstr)-       &
!!   &                               cff2*DiaRVfrc(i,j,nnew,M2bstr)+    &
!!   &                               cff3*DiaRVfrc(i,j,nstp,M2bstr)
!!              DiaRVfrc(i,j,nstp,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)
!!              DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+              &
!!   &                               cff1*DiaRVfrc(i,j,3,M2pgrd)-       &
!!   &                               cff2*DiaRVfrc(i,j,nnew,M2pgrd)+    &
!!   &                               cff3*DiaRVfrc(i,j,nstp,M2pgrd)
!!              DiaRVfrc(i,j,3,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)-          &
!!   &                                 DiaV2rhs(i,j,M2pgrd)
#  endif
!>              tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
!>
                ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rv(i,j,0,nstp)
                ad_rv(i,j,0,nstp)=0.0_r8
!>              tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                      &
!>   &                           cff1*tl_rvfrc(i,j)-                    &
!>   &                           cff2*tl_rv(i,j,0,nnew)+                &
!>   &                           cff3*tl_rv(i,j,0,nstp)
!>
                ad_rvfrc(i,j)=ad_rvfrc(i,j)+cff1*ad_rhs_vbar(i,j)
                ad_rv(i,j,0,nnew)=ad_rv(i,j,0,nnew)-                    &
     &                            cff2*ad_rhs_vbar(i,j)
                ad_rv(i,j,0,nstp)=ad_rv(i,j,0,nstp)+                    &
     &                            cff3*ad_rhs_vbar(i,j)
!>              tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
!>
                ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)-ad_rvfrc(i,j)
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef NEARSHORE_MELLOR
!!            DiaRUfrc(i,j,nstp,M2hrad)=DiaRUfrc(i,j,3,M2hrad)
!!            DiaU2rhs(i,j,M2hrad)=DiaU2rhs(i,j,M2hrad)+                &
!!   &                             cff1*DiaRUfrc(i,j,3,M2hrad)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2hrad)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2hrad)
!!            DiaRUfrc(i,j,3,M2hrad)=DiaRUfrc(i,j,3,M2hrad)-            &
!!   &                               DiaU2rhs(i,j,M2hrad)
#   endif
#   ifdef UV_ADV
!!              DiaRUfrc(i,j,nstp,M2hadv)=DiaRUfrc(i,j,3,M2hadv)
!!              DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+              &
!!   &                               cff1*DiaRUfrc(i,j,3,M2hadv)-       &
!!   &                               cff2*DiaRUfrc(i,j,nnew,M2hadv)+    &
!!   &                               cff3*DiaRUfrc(i,j,nstp,M2hadv)
!!              DiaRUfrc(i,j,3,M2hadv)=DiaRUfrc(i,j,3,M2hadv)-          &
!!   &                                 DiaU2rhs(i,j,M2hadv)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!              DiaRUfrc(i,j,nstp,M2hvis)=DiaRUfrc(i,j,3,M2hvis)
!!              DiaU2rhs(i,j,M2hvis)=DiaU2rhs(i,j,M2hvis)+              &
!!   &                               cff1*DiaRUfrc(i,j,3,M2hvis)-       &
!!   &                               cff2*DiaRUfrc(i,j,nnew,M2hvis)+    &
!!   &                               cff3*DiaRUfrc(i,j,nstp,M2hvis)
!!              DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)-          &
!!   &                                 DiaU2rhs(i,j,M2hvis)
#   endif
#   ifdef UV_COR
!!              DiaRUfrc(i,j,nstp,M2fcor)=DiaRUfrc(i,j,3,M2fcor)
!!              DiaU2rhs(i,j,M2fcor)=DiaU2rhs(i,j,M2fcor)+              &
!!   &                               cff1*DiaRUfrc(i,j,3,M2fcor)-       &
!!   &                               cff2*DiaRUfrc(i,j,nnew,M2fcor)+    &
!!   &                               cff3*DiaRUfrc(i,j,nstp,M2fcor)
!!              DiaRUfrc(i,j,3,M2fcor)=DiaRUfrc(i,j,3,M2fcor)-          &
!!   &                                 DiaU2rhs(i,j,M2fcor)
#   endif
!!              DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!              DiaU2rhs(i,j,M2sstr)=cff1*DiaRUfrc(i,j,3,M2sstr)-       &
!!   &                               cff2*DiaRUfrc(i,j,nnew,M2sstr)+    &
!!   &                               cff3*DiaRUfrc(i,j,nstp,M2sstr)
!!              DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!              DiaU2rhs(i,j,M2bstr)=cff1*DiaRUfrc(i,j,3,M2bstr)-       &
!!   &                               cff2*DiaRUfrc(i,j,nnew,M2bstr)+    &
!!   &                               cff3*DiaRUfrc(i,j,nstp,M2bstr)
!!              DiaRUfrc(i,j,nstp,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)
!!              DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+              &
!!   &                               cff1*DiaRUfrc(i,j,3,M2pgrd)-       &
!!   &                               cff2*DiaRUfrc(i,j,nnew,M2pgrd)+    &
!!   &                               cff3*DiaRUfrc(i,j,nstp,M2pgrd)
!!              DiaRUfrc(i,j,3,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)-          &
!!   &                                 DiaU2rhs(i,j,M2pgrd)
#  endif
!>              tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
!>
                ad_rufrc(i,j)=ad_rufrc(i,j)+ad_ru(i,j,0,nstp)
                ad_ru(i,j,0,nstp)=0.0_r8
!>              tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                      &
!>   &                           cff1*tl_rufrc(i,j)-                    &
!>   &                           cff2*tl_ru(i,j,0,nnew)+                &
!>   &                           cff3*tl_ru(i,j,0,nstp)
!>
                ad_rufrc(i,j)=ad_rufrc(i,j)+cff1*ad_rhs_ubar(i,j)
                ad_ru(i,j,0,nnew)=ad_ru(i,j,0,nnew)-                    &
     &                            cff2*ad_rhs_ubar(i,j)
                ad_ru(i,j,0,nstp)=ad_ru(i,j,0,nstp)+                    &
     &                            cff3*ad_rhs_ubar(i,j)
!>              tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
!>
                ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)-ad_rufrc(i,j)
              END DO
            END DO
          END IF
        ELSE
          DO j=JstrV,Jend
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef NEARSHORE_MELLOR
!!          DiaV2rhs(i,j,M2hrad)=DiaV2rhs(i,j,M2hrad)+                  &
!!   &                           DiaRVfrc(i,j,3,M2hrad)
#   endif
#   ifdef UV_ADV
!!            DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)+                &
!!   &                             DiaRVfrc(i,j,3,M2hadv)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!            DiaV2rhs(i,j,M2hvis)=DiaV2rhs(i,j,M2hvis)+                &
!!   &                             DiaRVfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_COR
!!            DiaV2rhs(i,j,M2fcor)=DiaV2rhs(i,j,M2fcor)+                &
!!   &                             DiaRVfrc(i,j,3,M2fcor)
#   endif
!!            DiaV2rhs(i,j,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DiaV2rhs(i,j,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+                &
!!   &                             DiaRVfrc(i,j,3,M2pgrd)
#  endif
!>            tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
!>
              ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rhs_vbar(i,j)
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef NEARSHORE_MELLOR
!!          DiaU2rhs(i,j,M2hrad)=DiaU2rhs(i,j,M2hrad)+                  &
!!   &                           DiaRUfrc(i,j,3,M2hrad)
#   endif
#   ifdef UV_ADV
!!            DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+                &
!!   &                             DiaRUfrc(i,j,3,M2hadv)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!            DiaU2rhs(i,j,M2hvis)=DiaU2rhs(i,j,M2hvis)+                &
!!   &                             DiaRUfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_COR
!!            DiaU2rhs(i,j,M2fcor)=DiaU2rhs(i,j,M2fcor)+                &
!!   &                             DiaRUfrc(i,j,3,M2fcor)
#   endif
!!            DiaU2rhs(i,j,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DiaU2rhs(i,j,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+                &
!!   &                             DiaRUfrc(i,j,3,M2pgrd)
#  endif
!>            tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_rufrc(i,j)
!>
              ad_rufrc(i,j)=ad_rufrc(i,j)+ad_rhs_ubar(i,j)
            END DO
          END DO
        END IF
# else
!>
!>----------------------------------------------------------------------
!>  Add in surface momentum stress.
!>----------------------------------------------------------------------
!>
      DO j=Jstr,Jend
        DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2sstr)=fac
#  endif
!>        tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
!>
          ad_fac=ad_fac+ad_rhs_ubar(i,j)
!>        tl_fac=tl_sustr(i,j)*om_u(i,j)*on_u(i,j)
!>
          ad_sustr(i,j)=ad_sustr(i,j)+om_u(i,j)*on_u(i,j)*ad_fac
          ad_fac=0.0_r8
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2sstr)=fac
#  endif
!>        tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
!>
          ad_fac=ad_fac+ad_rhs_vbar(i,j)
!>        tl_fac=tl_svstr(i,j)*om_v(i,j)*on_v(i,j)
!>
          ad_svstr(i,j)=ad_svstr(i,j)+om_v(i,j)*on_v(i,j)*ad_fac
          ad_fac=0.0_r8
        END DO
      END DO
# endif
# ifdef M2CLM_NUDGING
!
!-----------------------------------------------------------------------
!  Add in adjoint nudging of 2D momentum climatology.
!-----------------------------------------------------------------------
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.25_r8*(M2nudgcof(i,j-1)+M2nudgcof(i,j))*              &
       &        om_v(i,j)*on_v(i,j)
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                          &
!>   &                       cff*((Drhs(i,j-1)+Drhs(i,j))*              &
!>   &                            (-tl_vbar(i,j,krhs))+                 &
!>   &                            (tl_Drhs(i,j-1)+tl_Drhs(i,j))*        &
!>   &                            (vbarclm(i,j)-vbar(i,j,krhs)))
!>
            adfac=cff*ad_rhs_vbar(i,j)
            adfac1=adfac*(Drhs(i,j-1)+Drhs(i,j))
            adfac2=adfac*(vbarclm(i,j)-vbar(i,j,krhs))
            ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)-adfac1
            ad_Drhs(i,j-1)=ad_Drhs(i,j-1)+adfac2
            ad_Drhs(i,j  )=ad_Drhs(i,j  )+adfac2
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.25_r8*(M2nudgcof(i-1,j)+M2nudgcof(i,j))*              &
     &          om_u(i,j)*on_u(i,j)
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                          &
!>   &                       cff*((Drhs(i-1,j)+Drhs(i,j))*              &
!>   &                            (-tl_ubar(i,j,krhs))+                 &
!>   &                            (tl_Drhs(i-1,j)+tl_Drhs(i,j))*        &
!>   &                            (ubarclm(i,j)-ubar(i,j,krhs)))
!>
            adfac=cff*ad_rhs_ubar(i,j)
            adfac1=adfac*(Drhs(i-1,j)+Drhs(i,j))
            adfac2=adfac*(ubarclm(i,j)-ubar(i,j,krhs))
            ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)-adfac1
            ad_Drhs(i-1,j)=ad_Drhs(i-1,j)+adfac2
            ad_Drhs(i  ,j)=ad_Drhs(i  ,j)+adfac2
          END DO
        END DO
# endif
# ifndef SOLVE3D
!
!-----------------------------------------------------------------------
!  Add in bottom stress.
!-----------------------------------------------------------------------
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2bstr)=-fac
#  endif
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
!>
            ad_fac=ad_fac-ad_rhs_vbar(i,j)
!>          tl_fac=tl_bvstr(i,j)*om_v(i,j)*on_v(i,j)
!>
            ad_bvstr(i,j)=ad_bvstr(i,j)+om_v(i,j)*on_v(i,j)*ad_fac
            ad_fac=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2bstr)=-fac
#  endif
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
!>
            ad_fac=ad_fac-ad_rhs_ubar(i,j)
!>          tl_fac=tl_bustr(i,j)*om_u(i,j)*on_u(i,j)
!>
            ad_bustr(i,j)=ad_bustr(i,j)+om_u(i,j)*on_u(i,j)*ad_fac
            ad_fac=0.0_r8
          END DO
        END DO
# else
#  ifdef DIAGNOSTICS_UV
!!
!!  Initialize the stress term if no bottom friction is defined.
!!
!!      DO j=Jstr,Jend
!!        DO i=IstrU,Iend
!!          DiaU2rhs(i,j,M2bstr)=0.0_r8
!!        END DO
!!      END DO
!!      DO j=JstrV,Jend
!!        DO i=Istr,Iend
!!          DiaV2rhs(i,j,M2bstr)=0.0_r8
!!        END DO
!!      END DO
#  endif
# endif
# if defined NEARSHORE_MELLOR && \
    (!defined SOLVE3D         || defined DIAGNOSTICS_UV)
!
!-----------------------------------------------------------------------
!  Add in radiation stress terms.
!-----------------------------------------------------------------------
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2hrad)=-cff1
#  endif
#  ifndef SOLVE3D
!>        tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_cff1-tl_cff2
!>
          ad_cff2=ad_cff2-ad_rhs_vbar(i,j)
          ad_cff1=ad_cff1-ad_rhs_vbar(i,j)
#  endif
!>        tl_cff2=tl_rvlag2d(i,j)
!>
          ad_rvlag2d(i,j)=ad_rvlag2d(i,j)+ad_cff2
          ad_cff2=0.0_r8
!>        tl_cff1=tl_rvstr2d(i,j)*om_v(i,j)*on_v(i,j)
!>
          ad_rvstr2d(i,j)=ad_rvstr2d(i,j)+                              &
     &                    om_v(i,j)*on_v(i,j)*ad_cff1
          ad_cff1=0.0_r8
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2hrad)=-cff1
#  endif
#  ifndef SOLVE3D
!>        tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_cff1-tl_cff2
!>
          ad_cff2=ad_cff2-ad_rhs_ubar(i,j)
          ad_cff1=ad_cff1-ad_rhs_ubar(i,j)
#  endif
!>        tl_cff2=tl_rulag2d(i,j)
!>
          ad_rulag2d(i,j)=ad_rulag2d(i,j)+ad_cff2
          ad_cff2=0.0_r8
!>        tl_cff1=tl_rustr2d(i,j)*om_u(i,j)*on_u(i,j)
!>
          ad_rustr2d(i,j)=ad_rustr2d(i,j)+                              &
     &                    om_u(i,j)*on_u(i,j)*ad_cff1
          ad_cff1=0.0_r8
        END DO
      END DO
# endif
# if defined UV_VIS2 || defined UV_VIS4
!
!-----------------------------------------------------------------------
!  Compute basic state total depth at PSI-points.
!-----------------------------------------------------------------------
!
#  ifdef UV_VIS4
#   ifdef EW_PERIODIC
#    define IV_RANGE Istr-1,Iend+1
#    define IU_RANGE Istr-1,Iend+1
#   else
#    define IV_RANGE MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
#    define IU_RANGE MAX(2,IstrU-1),MIN(Iend+1,Lm(ng))
#   endif
#   ifdef NS_PERIODIC
#    define JU_RANGE Jstr-1,Jend+1
#    define JV_RANGE Jstr-1,Jend+1
#   else
#    define JU_RANGE MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
#    define JV_RANGE MAX(2,JstrV-1),MIN(Jend+1,Mm(ng))
#   endif
!
        DO j=JU_RANGE+1
          DO i=IV_RANGE+1
#  else
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
#  endif
            Drhs_p(i,j)=0.25_r8*(Drhs(i,j  )+Drhs(i-1,j  )+             &
     &                           Drhs(i,j-1)+Drhs(i-1,j-1))
          END DO
        END DO 
# endif
# ifdef UV_VIS4
!
!-----------------------------------------------------------------------
!  Add in adjoint horizontal biharmonic viscosity. The biharmonic
!  operator is computed by applying the harmonic operator twice.
!-----------------------------------------------------------------------
!
!  Compute flux-components of the horizontal divergence of the
!  BASIC STATE stress tensor (m4 s^-3/2) in XI- and ETA-directions.
!
        DO j=-1+JV_RANGE
          DO i=-1+IU_RANGE
            cff=visc4_r(i,j)*0.5_r8*                                    &
     &          (pmon_r(i,j)*                                           &
     &           ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,krhs)-               &
     &            (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,krhs))-              &
     &           pnom_r(i,j)*                                           &
     &           ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,krhs)-               &
     &            (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,krhs)))
            UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
            VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
          END DO
        END DO
        DO j=JU_RANGE+1
          DO i=IV_RANGE+1
            cff=visc4_p(i,j)*0.5_r8*                                    &
     &          (pmon_p(i,j)*                                           &
     &           ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,krhs)-             &
     &            (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,krhs))+            &
     &           pnom_p(i,j)*                                           &
     &           ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,krhs)-             &
     &            (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,krhs)))
#  ifdef MASKING
            cff=cff*pmask(i,j)
#  endif
            UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
            VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
          END DO
        END DO
!
!  Compute BASIC STATE first harmonic operator (m s^-3/2).
!
        DO j=JU_RANGE
          DO i=IU_RANGE
            LapU(i,j)=0.125_r8*                                         &
     &                (pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*          &
     &                ((pn(i-1,j)+pn(i,j))*                             &
     &                 (UFx(i,j  )-UFx(i-1,j))+                         &
     &                 (pm(i-1,j)+pm(i,j))*                             &
     &                 (UFe(i,j+1)-UFe(i  ,j)))
          END DO
        END DO
        DO j=JV_RANGE
          DO i=IV_RANGE
            LapV(i,j)=0.125_r8*                                         &
     &                (pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))*          &
     &                ((pn(i,j-1)+pn(i,j))*                             &
     &                 (VFx(i+1,j)-VFx(i,j  ))-                         &
     &                 (pm(i,j-1)+pm(i,j))*                             &
     &                 (VFe(i  ,j)-VFe(i,j-1)))
          END DO
        END DO
!
!  Apply boundary conditions (other than periodic) to the first
!  BASIC STATE harmonic operator. These are gradient or closed
!  (free slip or no slip) boundary conditions.
!
#  ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=JU_RANGE
#   ifdef WESTERN_WALL
            LapU(IstrU-1,j)=0.0_r8
#   else
            LapU(IstrU-1,j)=LapU(IstrU,j)
#   endif
          END DO
          DO j=JV_RANGE
#   ifdef WESTERN_WALL
            LapV(Istr-1,j)=gamma2(ng)*LapV(Istr,j)
#   else
            LapV(Istr-1,j)=0.0_r8
#   endif
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=JU_RANGE
#   ifdef EASTERN_WALL
            LapU(Iend+1,j)=0.0_r8
#   else
            LapU(Iend+1,j)=LapU(Iend,j)
#   endif
          END DO
          DO j=JV_RANGE
#   ifdef EASTERN_WALL
            LapV(Iend+1,j)=gamma2(ng)*LapV(Iend,j)
#   else
            LapV(Iend+1,j)=0.0_r8
#   endif
          END DO
        END IF
#  endif
#  ifndef NS_PERIODIC
        IF (SOUTHERN_EDGE) THEN
          DO i=IU_RANGE
#   ifdef SOUTHERN_WALL
            LapU(i,Jstr-1)=gamma2(ng)*LapU(i,Jstr)
#   else
            LapU(i,Jstr-1)=0.0_r8
#   endif
          END DO
          DO i=IV_RANGE
#   ifdef SOUTHERN_WALL
            LapV(i,JstrV-1)=0.0_r8
#   else
            LapV(i,JstrV-1)=LapV(i,JstrV)
#   endif
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=IU_RANGE
#   ifdef NORTHERN_WALL
            LapU(i,Jend+1)=gamma2(ng)*LapU(i,Jend)
#   else
            LapU(i,Jend+1)=0.0_r8
#   endif
          END DO
          DO i=IV_RANGE
#   ifdef NORTHERN_WALL
            LapV(i,Jend+1)=0.0_r8
#   else
            LapV(i,Jend+1)=LapV(i,Jend)
#   endif
          END DO
        END IF
#  endif
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          LapU(Istr  ,Jstr-1)=0.5_r8*(LapU(Istr+1,Jstr-1)+              &
     &                                LapU(Istr  ,Jstr  ))
          LapV(Istr-1,Jstr  )=0.5_r8*(LapV(Istr-1,Jstr+1)+              &
     &                                LapV(Istr  ,Jstr  ))
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          LapU(Iend+1,Jstr-1)=0.5_r8*(LapU(Iend  ,Jstr-1)+              &
     &                                LapU(Iend+1,Jstr  ))
          LapV(Iend+1,Jstr  )=0.5_r8*(LapV(Iend  ,Jstr  )+              &
     &                                LapV(Iend+1,Jstr+1))
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
          LapU(Istr  ,Jend+1)=0.5_r8*(LapU(Istr+1,Jend+1)+              &
     &                                LapU(Istr  ,Jend  ))
          LapV(Istr-1,Jend+1)=0.5_r8*(LapV(Istr  ,Jend+1)+              &
     &                                LapV(Istr-1,Jend  ))
        END IF
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
          LapU(Iend+1,Jend+1)=0.5_r8*(LapU(Iend  ,Jend+1)+              &
     &                                LapU(Iend+1,Jend  ))
          LapV(Iend+1,Jend+1)=0.5_r8*(LapV(Iend  ,Jend+1)+              &
     &                                LapV(Iend+1,Jend  ))
        END IF
#  endif
!
!  Add in adjoint biharmocnic viscosity
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  if defined DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2hvis)=-fac
#  endif
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
!>
            ad_fac=ad_fac-ad_rhs_vbar(i,j)
!>          tl_fac=0.5_r8*((pn(i,j-1)+pn(i,j))*                         &
!>   &                     (tl_VFx(i+1,j)-tl_VFx(i,j  ))-               &
!>   &                     (pm(i,j-1)+pm(i,j))*                         &
!>   &                     (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
!>
            adfac=0.5_r8*ad_fac
            adfac1=adfac*(pn(i,j-1)+pn(i,j))
            adfac2=adfac*(pm(i,j-1)+pm(i,j))
            ad_VFx(i,j  )=ad_VFx(i,j  )-adfac1
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+adfac1
            ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac2
            ad_VFe(i  ,j)=ad_VFe(i  ,j)-adfac2
            ad_fac=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  if defined DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2hvis)=-fac
#  endif
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
!>
            ad_fac=ad_fac-ad_rhs_ubar(i,j)
!>          tl_fac=0.5_r8*((pn(i-1,j)+pn(i,j))*                         &
!>   &                     (tl_UFx(i,j  )-tl_UFx(i-1,j))+               &
!>   &                     (pm(i-1,j)+pm(i,j))*                         &
!>   &                     (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
!>
            adfac=0.5_r8*ad_fac
            adfac1=adfac*(pn(i-1,j)+pn(i,j))
            adfac2=adfac*(pm(i-1,j)+pm(i,j))
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-adfac1
            ad_UFx(i,j  )=ad_UFx(i  ,j)+adfac1
            ad_UFe(i  ,j)=ad_UFe(i  ,j)-adfac2
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+adfac2
            ad_fac=0.0_r8
          END DO
        END DO
!
!  Compute flux-components of the horizontal divergence of the
!  adjoint biharmonic stress tensor (m4/s2) in XI- and ETA-directions.
!
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
!>          tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
!>          tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
!>
            ad_cff=ad_cff+                                              &
     &             on_p(i,j)*on_p(i,j)*ad_VFx(i,j)+                     &
     &             om_p(i,j)*om_p(i,j)*ad_UFe(i,j)
            ad_VFx(i,j)=0.0_r8
            ad_UFe(i,j)=0.0_r8
#  ifdef MASKING
!>          tl_cff=tl_cff*pmask(i,j)
!>
            ad_cff=ad_cff*pmask(i,j)
#  endif
!>          tl_cff=visc4_p(i,j)*0.5_r8*                                 &
!>   &             (tl_Drhs_p(i,j)*                                     &
!>   &              (pmon_p(i,j)*                                       &
!>   &               ((pn(i  ,j-1)+pn(i  ,j))*LapV(i  ,j)-              &
!>   &                (pn(i-1,j-1)+pn(i-1,j))*LapV(i-1,j))+             &
!>   &               pnom_p(i,j)*                                       &
!>   &               ((pm(i-1,j  )+pm(i,j  ))*LapU(i,j  )-              &
!>   &                (pm(i-1,j-1)+pm(i,j-1))*LapU(i,j-1)))+            &
!>   &              Drhs_p(i,j)*                                        &
!>   &              (pmon_p(i,j)*                                       &
!>   &               ((pn(i  ,j-1)+pn(i  ,j))*tl_LapV(i  ,j)-           &
!>   &                (pn(i-1,j-1)+pn(i-1,j))*tl_LapV(i-1,j))+          &
!>   &               pnom_p(i,j)*                                       &
!>   &               ((pm(i-1,j  )+pm(i,j  ))*tl_LapU(i,j  )-           &
!>   &                (pm(i-1,j-1)+pm(i,j-1))*tl_LapU(i,j-1))))
!>
            adfac=visc4_p(i,j)*0.5_r8*ad_cff
            adfac1=adfac*Drhs_p(i,j)*pmon_p(i,j)
            adfac2=adfac*Drhs_p(i,j)*pnom_p(i,j)
            ad_Drhs_p(i,j)=ad_Drhs_p(i,j)+                              &
     &                     (pmon_p(i,j)*                                &
     &                      ((pn(i  ,j-1)+pn(i  ,j))*LapV(i  ,j)-       &
     &                       (pn(i-1,j-1)+pn(i-1,j))*LapV(i-1,j))+      &
     &                      pnom_p(i,j)*                                &
     &                      ((pm(i-1,j  )+pm(i,j  ))*LapU(i,j  )-       &
     &                       (pm(i-1,j-1)+pm(i,j-1))*LapU(i,j-1)))*     &
     &                     adfac
            ad_LapV(i  ,j)=ad_LapV(i  ,j)+                              &
     &                     (pn(i  ,j-1)+pn(i  ,j))*adfac1
            ad_LapV(i-1,j)=ad_LapV(i-1,j)-                              &
     &                     (pn(i-1,j-1)+pn(i-1,j))*adfac1
            ad_LapU(i,j  )=ad_LapU(i,j  )+                              &
     &                     (pm(i-1,j  )+pm(i,j  ))*adfac2
            ad_LapU(i,j-1)=ad_LapU(i,j-1)-                              &
     &                     (pm(i-1,j-1)+pm(i,j-1))*adfac2
            ad_cff=0.0_r8
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
!>          tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
!>          tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
!>
            ad_cff=ad_cff+                                              &
     &             om_r(i,j)*om_r(i,j)*ad_VFe(i,j)+                     &
     &             on_r(i,j)*on_r(i,j)*ad_UFx(i,j)
            ad_VFe(i,j)=0.0_r8
            ad_UFx(i,j)=0.0_r8
!>          tl_cff=visc4_r(i,j)*0.5_r8*                                 &
!>   &             (tl_Drhs(i,j)*                                       &
!>   &              (pmon_r(i,j)*                                       &
!>   &               ((pn(i  ,j)+pn(i+1,j))*LapU(i+1,j)-                &
!>   &                (pn(i-1,j)+pn(i  ,j))*LapU(i  ,j))-               &
!>   &               pnom_r(i,j)*                                       &
!>   &               ((pm(i,j  )+pm(i,j+1))*LapV(i,j+1)-                &
!>   &                (pm(i,j-1)+pm(i,j  ))*LapV(i,j  )))+              &
!>   &              Drhs(i,j)*                                          &
!>   &              (pmon_r(i,j)*                                       &
!>   &               ((pn(i  ,j)+pn(i+1,j))*tl_LapU(i+1,j)-             &
!>   &                (pn(i-1,j)+pn(i  ,j))*tl_LapU(i  ,j))-            &
!>   &               pnom_r(i,j)*                                       &
!>   &               ((pm(i,j  )+pm(i,j+1))*tl_LapV(i,j+1)-             &
!>   &                (pm(i,j-1)+pm(i,j  ))*tl_LapV(i,j  ))))
!>
            adfac=visc4_r(i,j)*0.5_r8*ad_cff
            adfac1=adfac*Drhs(i,j)*pmon_r(i,j)
            adfac2=adfac*Drhs(i,j)*pnom_r(i,j)
            ad_Drhs(i,j)=ad_Drhs(i,j)+                                  &
     &                   (pmon_r(i,j)*                                  &
     &                    ((pn(i  ,j)+pn(i+1,j))*LapU(i+1,j)-           &
     &                     (pn(i-1,j)+pn(i  ,j))*LapU(i  ,j))-          &
     &                    pnom_r(i,j)*                                  &
     &                    ((pm(i,j  )+pm(i,j+1))*LapV(i,j+1)-           &
     &                     (pm(i,j-1)+pm(i,j  ))*LapV(i,j  )))*adfac
            ad_LapU(i+1,j)=ad_LapU(i+1,j)+                              &
     &                     (pn(i  ,j)+pn(i+1,j))*adfac1
            ad_LapU(i  ,j)=ad_LapU(i  ,j)-                              &
     &                     (pn(i-1,j)+pn(i  ,j))*adfac1
            ad_LapV(i,j+1)=ad_LapV(i,j+1)-                              &
     &                     (pm(i,j  )+pm(i,j+1))*adfac2
            ad_LapV(i,j  )=ad_LapV(i,j  )+                              &
     &                     (pm(i,j-1)+pm(i,j  ))*adfac2
            ad_cff=0.0_r8
          END DO
        END DO
!
!  Apply boundary conditions (other than periodic) to the first
!  adjoint harmonic operator. These are gradient or closed (free
!  slip or no slip) boundary conditions.
!
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
        IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>        tl_LapV(Iend+1,Jend+1)=0.5_r8*(tl_LapV(Iend  ,Jend+1)+        &
!>   &                                   tl_LapV(Iend+1,Jend  ))
!>
          adfac=0.5_r8*ad_LapV(Iend+1,Jend+1)
          ad_LapV(Iend+1,Jend  )=ad_LapV(Iend+1,Jend  )+adfac
          ad_LapV(Iend  ,Jend+1)=ad_LapV(Iend  ,Jend+1)+adfac
          ad_LapV(Iend+1,Jend+1)=0.0_r8
!>        tl_LapU(Iend+1,Jend+1)=0.5_r8*(tl_LapU(Iend  ,Jend+1)+        &
!>   &                                   tl_LapU(Iend+1,Jend  ))
!>
          adfac=0.5_r8*ad_LapU(Iend+1,Jend+1)
          ad_LapU(Iend+1,Jend  )=ad_LapU(Iend+1,Jend  )+adfac
          ad_LapU(Iend  ,Jend+1)=ad_LapU(Iend  ,Jend+1)+adfac
          ad_LapU(Iend+1,Jend+1)=0.0_r8
        END IF
        IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>        tl_LapV(Istr-1,Jend+1)=0.5_r8*(tl_LapV(Istr  ,Jend+1)+        &
!>   &                                   tl_LapV(Istr-1,Jend ))
!>
          adfac=0.5_r8*ad_LapV(Istr-1,Jend+1)
          ad_LapV(Istr-1,Jend  )=ad_LapV(Istr-1,Jend  )+adfac
          ad_LapV(Istr  ,Jend+1)=ad_LapV(Istr  ,Jend+1)+adfac
          ad_LapV(Istr-1,Jend+1)=0.0_r8
!>        tl_LapU(Istr  ,Jend+1)=0.5_r8*(tl_LapU(Istr+1,Jend+1)+        &
!>   &                                   tl_LapU(Istr  ,Jend  ))
!>
          adfac=0.5_r8*ad_LapU(Istr  ,Jend+1)
          ad_LapU(Istr  ,Jend  )=ad_LapU(Istr  ,Jend  )+adfac
          ad_LapU(Istr+1,Jend+1)=ad_LapU(Istr+1,Jend+1)+adfac
          ad_LapU(Istr  ,Jend+1)=0.0_r8
        END IF
        IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
!>        tl_LapV(Iend+1,Jstr  )=0.5_r8*(tl_LapV(Iend  ,Jstr  )+        &
!>   &                                   tl_LapV(Iend+1,Jstr+1))
!>
          adfac=0.5_r8*ad_LapV(Iend+1,Jstr  )
          ad_LapV(Iend  ,Jstr  )=ad_LapV(Iend  ,Jstr  )+adfac
          ad_LapV(Iend+1,Jstr+1)=ad_LapV(Iend+1,Jstr+1)+adfac
          ad_LapV(Iend+1,Jstr  )=0.0_r8
!>        tl_LapU(Iend+1,Jstr-1)=0.5_r8*(tl_LapU(Iend  ,Jstr-1)+        &
!>   &                                   tl_LapU(Iend+1,Jstr  ))
!>
          adfac=0.5_r8*ad_LapU(Iend+1,Jstr-1)
          ad_LapU(Iend  ,Jstr-1)=ad_LapU(Iend  ,Jstr-1)+adfac
          ad_LapU(Iend+1,Jstr  )=ad_LapU(Iend+1,Jstr  )+adfac
          ad_LapU(Iend+1,Jstr-1)=0.0_r8
        END IF
        IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
!>        tl_LapV(Istr-1,Jstr  )=0.5_r8*(tl_LapV(Istr-1,Jstr+1)+        &
!>   &                                   tl_LapV(Istr  ,Jstr  ))
!>
          adfac=0.5_r8*ad_LapV(Istr-1,Jstr  )
          ad_LapV(Istr  ,Jstr  )=ad_LapV(Istr  ,Jstr  )+adfac
          ad_LapV(Istr-1,Jstr+1)=ad_LapV(Istr-1,Jstr+1)+adfac
          ad_LapV(Istr-1,Jstr  )=0.0_r8
!>        tl_LapU(Istr  ,Jstr-1)=0.5_r8*(tl_LapU(Istr+1,Jstr-1)+        &
!>   &                                   tl_LapU(Istr  ,Jstr  ))
!>
          adfac=0.5_r8*ad_LapU(Istr  ,Jstr-1)
          ad_LapU(Istr+1,Jstr-1)=ad_LapU(Istr+1,Jstr-1)+adfac
          ad_LapU(Istr  ,Jstr  )=ad_LapU(Istr  ,Jstr  )+adfac
          ad_LapU(Istr  ,Jstr-1)=0.0_r8
        END IF
#  endif
#  ifndef NS_PERIODIC
        IF (NORTHERN_EDGE) THEN
          DO i=IV_RANGE
#   ifdef NORTHERN_WALL
!>          tl_LapV(i,Jend+1)=0.0_r8
!>
            ad_LapV(i,Jend+1)=0.0_r8
#   else
!>          tl_LapV(i,Jend+1)=tl_LapV(i,Jend)
!>
            ad_LapV(i,Jend)=ad_LapV(i,Jend)+ad_LapV(i,Jend+1)
            ad_LapV(i,Jend+1)=0.0_r8
#   endif
          END DO
          DO i=IU_RANGE
#   ifdef NORTHERN_WALL
!>          tl_LapU(i,Jend+1)=gamma2(ng)*tl_LapU(i,Jend)
!>
            ad_LapU(i,Jend)=ad_LapU(i,Jend)+gamma2(ng)*ad_LapU(i,Jend+1)
            ad_LapU(i,Jend+1)=0.0_r8
#   else
!>          tl_LapU(i,Jend+1)=0.0_r8
!>
            ad_LapU(i,Jend+1)=0.0_r8
#   endif
          END DO
        END IF
        IF (SOUTHERN_EDGE) THEN
          DO i=IV_RANGE
#   ifdef SOUTHERN_WALL
!>          tl_LapV(i,JstrV-1)=0.0_r8
!>
            ad_LapV(i,JstrV-1)=0.0_r8
#   else
!>          tl_LapV(i,JstrV-1)=tl_LapV(i,JstrV)
!>
            ad_LapV(i,JstrV)=ad_LapV(i,JstrV)+ad_LapV(i,JstrV-1)
            ad_LapV(i,JstrV-1)=0.0_r8
#   endif
          END DO
          DO i=IU_RANGE
#   ifdef SOUTHERN_WALL
!>          tl_LapU(i,Jstr-1)=gamma2(ng)*tl_LapU(i,Jstr)
!>
            ad_LapU(i,Jstr)=ad_LapU(i,Jstr)+gamma2(ng)*ad_LapU(i,Jstr-1)
            ad_LapU(i,Jstr-1)=0.0_r8
#   else
!>          tl_LapU(i,Jstr-1)=0.0_r8
!>
            ad_LapU(i,Jstr-1)=0.0_r8
#   endif
          END DO
        END IF
#  endif
#  ifndef EW_PERIODIC
        IF (EASTERN_EDGE) THEN
          DO j=JV_RANGE
#   ifdef EASTERN_WALL
!>          tl_LapV(Iend+1,j)=gamma2(ng)*tl_LapV(Iend,j)
!>
            ad_LapV(Iend,j)=ad_LapV(Iend,j)+gamma2(ng)*ad_LapV(Iend+1,j)
            ad_LapV(Iend+1,j)=0.0_r8
#   else
!>          tl_LapV(Iend+1,j)=0.0_r8
!>
            ad_LapV(Iend+1,j)=0.0_r8
#   endif
          END DO
          DO j=JU_RANGE
#   ifdef EASTERN_WALL
!>          tl_LapU(Iend+1,j)=0.0_r8
!>
            ad_LapU(Iend+1,j)=0.0_r8
#   else
!>          tl_LapU(Iend+1,j)=tl_LapU(Iend,j)
!>
            ad_LapU(Iend,j)=ad_LapU(Iend,j)+ad_LapU(Iend+1,j)
            ad_LapU(Iend+1,j)=0.0_r8
#   endif
          END DO
        END IF
        IF (WESTERN_EDGE) THEN
          DO j=JV_RANGE
#   ifdef WESTERN_WALL
!>          tl_LapV(Istr-1,j)=gamma2(ng)*tl_LapV(Istr,j)
!>
            ad_LapV(Istr,j)=ad_LapV(Istr,j)+gamma2(ng)*ad_LapV(Istr-1,j)
            ad_LapV(Istr-1,j)=0.0_r8
#   else
!>          tl_LapV(Istr-1,j)=0.0_r8
!>
            ad_LapV(Istr-1,j)=0.0_r8
#   endif
          END DO
          DO j=JU_RANGE
#   ifdef WESTERN_WALL
!>          tl_LapU(IstrU-1,j)=0.0_r8
!>
            ad_LapU(IstrU-1,j)=0.0_r8
#   else
!>          tl_LapU(IstrU-1,j)=tl_LapU(IstrU,j)
!>
            ad_LapU(IstrU,j)=ad_LapU(IstrU,j)+ad_LapU(IstrU-1,j)
            ad_LapU(IstrU-1,j)=0.0_r8
#   endif
          END DO
        END IF
#  endif
!
!  Compute adjoint first harmonic operator (m s^-3/2).
!
        DO j=JV_RANGE
          DO i=IV_RANGE
!>          tl_LapV(i,j)=0.125_r8*                                      &
!>   &                   (pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))*       &
!>   &                   ((pn(i,j-1)+pn(i,j))*                          &
!>   &                    (tl_VFx(i+1,j)-tl_VFx(i,j  ))-                &
!>   &                    (pm(i,j-1)+pm(i,j))*                          &
!>   &                    (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
!>
            adfac=0.125_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))*     &
     &            ad_LapV(i,j)
            adfac1=adfac*(pn(i,j-1)+pn(i,j))
            adfac2=adfac*(pm(i,j-1)+pm(i,j))
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-adfac1
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+adfac1
            ad_VFe(i,j  )=ad_VFe(i,j  )-adfac2
            ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac2
            ad_LapV(i,j)=0.0_r8
          END DO
        END DO
!
        DO j=JU_RANGE
          DO i=IU_RANGE
!>          tl_LapU(i,j)=0.125_r8*                                      &
!>   &                   (pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*       &
!>   &                   ((pn(i-1,j)+pn(i,j))*                          &
!>   &                    (tl_UFx(i,j  )-tl_UFx(i-1,j))+                &
!>   &                    (pm(i-1,j)+pm(i,j))*                          &
!>   &                    (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
!>
            adfac=0.125_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*     &
     &            ad_LapU(i,j)
            adfac1=adfac*(pn(i-1,j)+pn(i,j))
            adfac2=adfac*(pm(i-1,j)+pm(i,j))
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-adfac1
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac1
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+adfac2
            ad_UFe(i,j  )=ad_UFe(i,j  )-adfac2
            ad_LapU(i,j)=0.0_r8
          END DO
        END DO
!
!  Compute flux-components of the adjoint horizontal divergence of the
!  stress tensor (m4 s^-3/2) in XI- and ETA-directions. It is assumed
!  here that "visc4_r" and "visc4_p" are the squared root of the
!  biharmonic viscosity coefficient.  For momentum balance purposes,
!  the total thickness "D" appears only when computing the second
!  harmonic operator.
!
        DO j=JU_RANGE+1
          DO i=IV_RANGE+1
!>          tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
!>          tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
!>
            ad_cff=ad_cff+                                              &
     &             on_p(i,j)*on_p(i,j)*ad_VFx(i,j)+                     &
     &             om_p(i,j)*om_p(i,j)*ad_UFe(i,j)
            ad_VFx(i,j)=0.0_r8
            ad_UFe(i,j)=0.0_r8
#  ifdef MASKING
!>          tl_cff=tl_cff*pmask(i,j)
!>
            ad_cff=ad_cff*pmask(i,j)
#  endif
!>          tl_cff=visc4_p(i,j)*0.5_r8*                                 &
!>   &             (pmon_p(i,j)*                                        &
!>   &              ((pn(i  ,j-1)+pn(i  ,j))*tl_vbar(i  ,j,krhs)-       &
!>   &               (pn(i-1,j-1)+pn(i-1,j))*tl_vbar(i-1,j,krhs))+      &
!>   &              pnom_p(i,j)*                                        &
!>   &              ((pm(i-1,j  )+pm(i,j  ))*tl_ubar(i,j  ,krhs)-       &
!>   &               (pm(i-1,j-1)+pm(i,j-1))*tl_ubar(i,j-1,krhs)))
!>
            adfac=visc4_p(i,j)*0.5_r8*ad_cff
            adfac1=adfac*pmon_p(i,j)
            adfac2=adfac*pnom_p(i,j)
            ad_vbar(i-1,j,krhs)=ad_vbar(i-1,j,krhs)-                    &
     &                          (pn(i-1,j-1)+pn(i-1,j))*adfac1
            ad_vbar(i  ,j,krhs)=ad_vbar(i  ,j,krhs)+                    &
     &                          (pn(i  ,j-1)+pn(i  ,j))*adfac1
            ad_ubar(i,j-1,krhs)=ad_ubar(i,j-1,krhs)-                    &
     &                          (pm(i-1,j-1)+pm(i,j-1))*adfac2
            ad_ubar(i,j  ,krhs)=ad_ubar(i,j  ,krhs)+                    &
     &                          (pm(i-1,j  )+pm(i,j  ))*adfac2
            ad_cff=0.0_r8
          END DO
        END DO
        DO j=-1+JV_RANGE
          DO i=-1+IU_RANGE
!>          tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
!>          tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
!>
            ad_cff=ad_cff+                                              &
     &             om_r(i,j)*om_r(i,j)*ad_VFe(i,j)+                     &
     &             on_r(i,j)*on_r(i,j)*ad_UFx(i,j)
            ad_VFe(i,j)=0.0_r8
            ad_UFx(i,j)=0.0_r8
!>          tl_cff=visc4_r(i,j)*0.5_r8*                                 &
!>   &             (pmon_r(i,j)*                                        &
!>   &              ((pn(i  ,j)+pn(i+1,j))*tl_ubar(i+1,j,krhs)-         &
!>   &               (pn(i-1,j)+pn(i  ,j))*tl_ubar(i  ,j,krhs))-        &
!>   &              pnom_r(i,j)*                                        &
!>   &              ((pm(i,j  )+pm(i,j+1))*tl_vbar(i,j+1,krhs)-         &
!>   &               (pm(i,j-1)+pm(i,j  ))*tl_vbar(i,j  ,krhs)))
!>
            adfac=visc4_r(i,j)*0.5_r8*ad_cff
            adfac1=adfac*pmon_r(i,j)
            adfac2=adfac*pnom_r(i,j)
            ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+                    &
     &                          (pn(i  ,j)+pn(i+1,j))*adfac1
            ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)-                    &
     &                          (pn(i-1,j)+pn(i  ,j))*adfac1
            ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)-                    &
     &                          (pm(i,j  )+pm(i,j+1))*adfac2
            ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+                    &
     &                          (pm(i,j-1)+pm(i,j  ))*adfac2
            ad_cff=0.0_r8
          END DO
        END DO
# endif
# ifdef UV_VIS2
!
!-----------------------------------------------------------------------
!  Add in adjoint horizontal harmonic viscosity.
!-----------------------------------------------------------------------
!
!  Add in harmonic viscosity.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  if defined DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2hvis)=fac
#  endif
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
!>
            ad_fac=ad_fac+ad_rhs_vbar(i,j)
!>          tl_fac=0.5_r8*((pn(i,j-1)+pn(i,j))*                         &
!>   &                     (tl_VFx(i+1,j)-tl_VFx(i,j  ))-               &
!>   &                     (pm(i,j-1)+pm(i,j))*                         &
!>   &                     (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
!>
            adfac=0.5_r8*ad_fac
            adfac1=adfac*(pn(i,j-1)+pn(i,j))
            adfac2=adfac*(pm(i,j-1)+pm(i,j))
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-adfac1
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+adfac1
            ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac2
            ad_VFe(i,j  )=ad_VFe(i,j  )-adfac2
            ad_fac=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  if defined DIAGNOSTICS_UV
!>          DiaU2rhs(i,j,M2hvis)=fac
#  endif
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
!>
            ad_fac=ad_fac+ad_rhs_ubar(i,j)
!>          tl_fac=0.5_r8*((pn(i-1,j)+pn(i,j))*                         &
!>   &                     (tl_UFx(i,j  )-tl_UFx(i-1,j))+               &
!>   &                     (pm(i-1,j)+pm(i,j))*                         &
!>   &                     (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
!>
            adfac=0.5_r8*ad_fac
            adfac1=adfac*(pn(i-1,j)+pn(i,j))
            adfac2=adfac*(pm(i-1,j)+pm(i,j))
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-adfac1
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac1
            ad_UFe(i,j  )=ad_UFe(i,j  )-adfac2
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+adfac2
            ad_fac=0.0_r8
          END DO
        END DO
!
!  Compute flux-components of the adjoint horizontal divergence of the
!  stress tensor (m5/s2) in XI- and ETA-directions.
!
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
!>          tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
!>          tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
!>
            ad_cff=ad_cff+                                              &
     &             on_p(i,j)*on_p(i,j)*ad_VFx(i,j)+                     &
     &             om_p(i,j)*om_p(i,j)*ad_UFe(i,j)
            ad_VFx(i,j)=0.0_r8
            ad_UFe(i,j)=0.0_r8
#  ifdef MASKING
!>          tl_cff=tl_cff*pmask(i,j)
!>
            ad_cff=ad_cff*pmask(i,j)
#  endif
!>          tl_cff=visc2_p(i,j)*0.5_r8*                                 &
!>   &             (tl_Drhs_p(i,j)*                                     &
!>   &              (pmon_p(i,j)*                                       &
!>   &               ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,krhs)-         &
!>   &                (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,krhs))+        &
!>   &               pnom_p(i,j)*                                       &
!>   &               ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,krhs)-         &
!>   &                (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,krhs)))+       &
!>   &              Drhs_p(i,j)*                                        &
!>   &              (pmon_p(i,j)*                                       &
!>   &               ((pn(i  ,j-1)+pn(i  ,j))*tl_vbar(i  ,j,krhs)-      &
!>   &                (pn(i-1,j-1)+pn(i-1,j))*tl_vbar(i-1,j,krhs))+     &
!>   &               pnom_p(i,j)*                                       &
!>   &               ((pm(i-1,j  )+pm(i,j  ))*tl_ubar(i,j  ,krhs)-      &
!>   &                (pm(i-1,j-1)+pm(i,j-1))*tl_ubar(i,j-1,krhs))))
!>
            adfac=visc2_p(i,j)*0.5_r8*ad_cff
            adfac1=adfac*Drhs_p(i,j)
            adfac2=adfac1*pmon_p(i,j)
            adfac3=adfac1*pnom_p(i,j)
            ad_Drhs_p(i,j)=ad_Drhs_p(i,j)+                              &
     &                     (pmon_p(i,j)*                                &
     &                      ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,krhs)-  &
     &                       (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,krhs))+ &
     &                      pnom_p(i,j)*                                &
     &                      ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,krhs)-  &
     &                       (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,krhs)))*&
     &                     adfac
            ad_vbar(i-1,j,krhs)=ad_vbar(i-1,j,krhs)-                    &
     &                          (pn(i-1,j-1)+pn(i-1,j))*adfac2
            ad_vbar(i  ,j,krhs)=ad_vbar(i  ,j,krhs)+                    &
     &                          (pn(i  ,j-1)+pn(i  ,j))*adfac2
            ad_ubar(i,j-1,krhs)=ad_ubar(i,j-1,krhs)-                    &
     &                          (pm(i-1,j-1)+pm(i,j-1))*adfac3
            ad_ubar(i,j  ,krhs)=ad_ubar(i,j  ,krhs)+                    &
     &                          (pm(i-1,j  )+pm(i,j  ))*adfac3
            ad_cff=0.0_r8
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
!>          tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
!>          tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
!>
            ad_cff=ad_cff+                                              &
     &             om_r(i,j)*om_r(i,j)*ad_VFe(i,j)+                     &
     &             on_r(i,j)*on_r(i,j)*ad_UFx(i,j)
            ad_VFe(i,j)=0.0_r8
            ad_UFx(i,j)=0.0_r8
!>          tl_cff=visc2_r(i,j)*0.5_r8*                                 &
!>   &             (tl_Drhs(i,j)*                                       &
!>   &              (pmon_r(i,j)*                                       &
!>   &               ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,krhs)-           &
!>   &                (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,krhs))-          &
!>   &               pnom_r(i,j)*                                       &
!>   &               ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,krhs)-           &
!>   &                (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,krhs)))+         &
!>   &              Drhs(i,j)*                                          &
!>   &              (pmon_r(i,j)*                                       &
!>   &               ((pn(i  ,j)+pn(i+1,j))*tl_ubar(i+1,j,krhs)-        &
!>   &                (pn(i-1,j)+pn(i  ,j))*tl_ubar(i  ,j,krhs))-       &
!>   &               pnom_r(i,j)*                                       &
!>   &               ((pm(i,j  )+pm(i,j+1))*tl_vbar(i,j+1,krhs)-        &
!>   &                (pm(i,j-1)+pm(i,j  ))*tl_vbar(i,j  ,krhs))))
!>
            adfac=visc2_r(i,j)*0.5_r8*ad_cff
            adfac1=adfac*Drhs(i,j)
            adfac2=adfac1*pmon_r(i,j)
            adfac3=adfac1*pnom_r(i,j)
            ad_Drhs(i,j)=ad_Drhs(i,j)+                                  &
     &                   (pmon_r(i,j)*                                  &
     &                    ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,krhs)-      &
     &                     (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,krhs))-     &
     &                    pnom_r(i,j)*                                  &
     &                    ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,krhs)-      &
     &                     (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,krhs)))*    &
     &                   adfac
            ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)-                    &
     &                          (pn(i-1,j)+pn(i  ,j))*adfac2
            ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+                    &
     &                          (pn(i  ,j)+pn(i+1,j))*adfac2
            ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+                    &
     &                          (pm(i,j-1)+pm(i,j  ))*adfac3
            ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)-                    &
     &                          (pm(i,j  )+pm(i,j+1))*adfac3
            ad_cff=0.0_r8
          END DO
        END DO
# endif
# if defined UV_VIS2 || defined UV_VIS4
!
!-----------------------------------------------------------------------
!  If horizontal mixing, compute adjoint total depth at PSI-points.
!-----------------------------------------------------------------------
!
#  ifdef UV_VIS4
        DO j=JU_RANGE+1
          DO i=IV_RANGE+1
#  else
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
#  endif
            Drhs_p(i,j)=0.25_r8*(Drhs(i,j  )+Drhs(i-1,j  )+             &
     &                           Drhs(i,j-1)+Drhs(i-1,j-1))
!>          tl_Drhs_p(i,j)=0.25_r8*(tl_Drhs(i,j  )+tl_Drhs(i-1,j  )+    &
!>   &                              tl_Drhs(i,j-1)+tl_Drhs(i-1,j-1))
!>
            adfac=0.25_r8*ad_Drhs_p(i,j)
            ad_Drhs(i-1,j  )=ad_Drhs(i-1,j  )+adfac
            ad_Drhs(i  ,j  )=ad_Drhs(i  ,j  )+adfac
            ad_Drhs(i-1,j-1)=ad_Drhs(i-1,j-1)+adfac
            ad_Drhs(i  ,j-1)=ad_Drhs(i  ,j-1)+adfac
            ad_Drhs_p(i,j)=0.0_r8
          END DO
        END DO
#  undef IU_RANGE
#  undef IV_RANGE
#  undef JU_RANGE
#  undef JV_RANGE
# endif
# if defined UV_COR || (defined CURVGRID && defined UV_ADV)
!
!-----------------------------------------------------------------------
!  Add in adjoint Coriolis and adjoint curvilinear transformation
!  terms, if any.
!-----------------------------------------------------------------------
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  if defined DIAGNOSTICS_UV
!!          fac2=0.5_r8*(Vwrk(i,j)+Vwrk(i,j-1))
#   if (defined CURVGRID && defined UV_ADV)
!!          DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)-fac1+fac2
#   endif
#   ifdef UV_COR
!!          DiaV2rhs(i,j,M2fcor)=-fac2
#   endif
#  endif
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac1
!>
            ad_fac1=ad_fac1-ad_rhs_vbar(i,j)
!>          tl_fac1=0.5_r8*(tl_VFe(i,j)+tl_VFe(i,j-1))
!>
            adfac=0.5_r8*ad_fac1
            ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac
            ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
            ad_fac1=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  if defined DIAGNOSTICS_UV
!!          fac2=0.5_r8*(Uwrk(i,j)+Uwrk(i-1,j))
#   if (defined CURVGRID && defined UV_ADV)
!!          DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+fac1-fac2
#   endif
#   ifdef UV_COR
!!          DiaU2rhs(i,j,M2fcor)=fac2
#   endif
#  endif
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac1
!>
            ad_fac1=ad_fac1+ad_rhs_ubar(i,j)
!>          tl_fac1=0.5_r8*(tl_UFx(i,j)+tl_UFx(i-1,j))
!>
            adfac=0.5_r8*ad_fac1
            ad_UFx(i-1,j)=ad_UFx(i-1,j)+adfac
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
            ad_fac1=0.0_r8
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=0.5_r8*Drhs(i,j)*(                                      &
#  ifdef UV_COR
     &          fomn(i,j)                                               &
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     &         +0.5_r8*((vbar(i,j  ,krhs)+                              &
#   ifdef NEARSHORE_MELLOR
     &                   vbar_stokes(i,j  )+                            &
     &                   vbar_stokes(i,j+1)+                            &
#   endif
     &                   vbar(i,j+1,krhs))*dndx(i,j)-                   &
#   ifdef NEARSHORE_MELLOR
     &                   ubar_stokes(i  ,j)+                            &
     &                   ubar_stokes(i+1,j)+                            &
#   endif
     &                  (ubar(i  ,j,krhs)+                              &
     &                   ubar(i+1,j,krhs))*dmde(i,j))                   &
#  endif
     &          )
#  if defined DIAGNOSTICS_UV
#   ifdef UV_COR
!!          Vwrk(i,j)=0.5_r8*Drhs(i,j)*fomn(i,j)*                       &
!!   &                (ubar(i,j,krhs)+ubar(i+1,j  ,krhs))
!!          Uwrk(i,j)=0.5_r8*Drhs(i,j)*fomn(i,j)*                       &
!!   &                (vbar(i,j,krhs)+vbar(i  ,j+1,krhs))
#   else
!!          Vwrk(i,j)=0.0_r8
!!          Uwrk(i,j)=0.0_r8
#   endif
#  endif
!>          tl_VFe(i,j)=tl_cff*(ubar(i  ,j,krhs)+                       &
#  ifdef NEARSHORE_MELLOR
!>   &                          ubar_stokes(i  ,j)+                     &
!>   &                          ubar_stokes(i+1,j)+                     &
#  endif
!>   &                          ubar(i+1,j,krhs))+                      &
!>   &                  cff*(tl_ubar(i  ,j,krhs)+                       &
#  ifdef NEARSHORE_MELLOR
!>   &                       tl_ubar_stokes(i  ,j)+                     &
!>   &                       tl_ubar_stokes(i+1,j)+                     &
#  endif
!>   &                       tl_ubar(i+1,j,krhs))
!>
            adfac=cff*ad_VFe(i,j)
            ad_cff=ad_cff+                                              &
     &             ad_VFe(i,j)*(ubar(i  ,j,krhs)+                       &
#  ifdef NEARSHORE_MELLOR
     &                          ubar_stokes(i  ,j)+                     &
     &                          ubar_stokes(i+1,j)+                     &
#  endif
     &                          ubar(i+1,j,krhs))
            ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)+adfac
            ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+adfac
#  ifdef NEARSHORE_MELLOR
            ad_ubar_stokes(i  ,j)=ad_ubar_stokes(i  ,j)+adfac
            ad_ubar_stokes(i+1,j)=ad_ubar_stokes(i+1,j)+adfac
#  endif
            ad_VFe(i,j)=0.0_r8
!>          tl_UFx(i,j)=tl_cff*(vbar(i,j  ,krhs)+                       &
#  ifdef NEARSHORE_MELLOR
!>   &                          vbar_stokes(i,j  )+                     &
!>   &                          vbar_stokes(i,j+1)+                     &
#  endif
!>   &                          vbar(i,j+1,krhs))+                      &
!>   &                  cff*(tl_vbar(i,j  ,krhs)+                       &
#  ifdef NEARSHORE_MELLOR
!>   &                       tl_vbar_stokes(i,j  )+                     &
!>   &                       tl_vbar_stokes(i,j+1)+                     &
#  endif
!>   &                       tl_vbar(i,j+1,krhs))
!>
            adfac=cff*ad_UFx(i,j)
            ad_cff=ad_cff+                                              &
     &             ad_UFx(i,j)*(vbar(i,j  ,krhs)+                       &
#  ifdef NEARSHORE_MELLOR
     &                          vbar_stokes(i,j  )+                     &
     &                          vbar_stokes(i,j+1)+                     &
#  endif
     &                          vbar(i,j+1,krhs))
            ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+adfac
            ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+adfac
#  ifdef NEARSHORE_MELLOR
            ad_vbar_stokes(i,j  )=ad_vbar_stokes(i,j  )+adfac
            ad_vbar_stokes(i,j+1)=ad_vbar_stokes(i,j+1)+adfac
#  endif
            ad_UFx(i,j)=0.0_r8
#  if (defined CURVGRID && defined UV_ADV)
!>          tl_cff=tl_cff+                                              &
!>   &             0.25_r8*Drhs(i,j)*                                   &
!>   &             ((tl_vbar(i,j  ,krhs)+                               &
#   ifdef NEARSHORE_MELLOR
!>   &               tl_vbar_stokes(i,j  )+                             &
!>   &               tl_vbar_stokes(i,j+1)+                             &
#   endif
!>   &               tl_vbar(i,j+1,krhs))*dndx(i,j)-                    &
!>   &              (tl_ubar(i,j,krhs)+
#   ifdef NEARSHORE_MELLOR
!>   &               tl_ubar_stokes(i  ,j)+                             &
!>   &               tl_ubar_stokes(i+1,j)+                             &
#   endif
!>   &               tl_ubar(i+1,j,krhs))*dmde(i,j))
!>
            adfac=0.25_r8*Drhs(i,j)*ad_cff
            adfac1=adfac*dndx(i,j)
            adfac2=adfac*dmde(i,j)
            ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+adfac1
            ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+adfac1
#   ifdef NEARSHORE_MELLOR
            ad_vbar_stokes(i,j  )=ad_vbar_stokes(i,j  )+adfac1
            ad_vbar_stokes(i,j+1)=ad_vbar_stokes(i,j+1)+adfac1
#   endif
            ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)-adfac2
            ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)-adfac2
#   ifdef NEARSHORE_MELLOR
            ad_ubar_stokes(i  ,j)=ad_ubar_stokes(i  ,j)-adfac2
            ad_ubar_stokes(i+1,j)=ad_ubar_stokes(i+1,j)-adfac2
#   endif
#  endif
!>          tl_cff=0.5_r8*tl_Drhs(i,j)*(                                &
#  ifdef UV_COR
!>   &             fomn(i,j)                                            &
#  endif
#  if (defined CURVGRID && defined UV_ADV)
!>   &            +0.5_r8*((vbar(i,j  ,krhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      vbar_stokes(i,j  )+                         &
!>   &                      vbar_stokes(i,j+1)+                         &
#   endif
!>   &                      vbar(i,j+1,krhs))*                          &
!>   &                     dndx(i,j)-                                   &
!>   &                     (ubar(i  ,j,krhs)+                           &
#  ifdef NEARSHORE_MELLOR
!>   &                      ubar_stokes(i  ,j)+                         &
!>   &                      ubar_stokes(i+1,j)+                         &
#  endif
!>   &                      ubar(i+1,j,krhs))*                          &
!>   &                     dmde(i,j))                                   &
#  endif
!>   &            )
!>
            adfac=0.5_r8*ad_cff
#  ifdef UV_COR
            ad_Drhs(i,j)=ad_Drhs(i,j)+adfac*fomn(i,j)
#  endif
#  if (defined CURVGRID && defined UV_ADV)
            ad_Drhs(i,j)=ad_Drhs(i,j)+                                  &
     &                   adfac*0.5_r8*((vbar(i  ,j  ,krhs)+             &
#   ifdef NEARSHORE_MELLOR
     &                                  vbar_stokes(i,j  )+             &
     &                                  vbar_stokes(i,j+1)+             &
#   endif
     &                                  vbar(i  ,j+1,krhs))*            &
     &                                 dndx(i,j)-                       &
     &                                 (ubar(i  ,j  ,krhs)+             &
#  ifdef NEARSHORE_MELLOR
     &                                  ubar_stokes(i  ,j)+             &
     &                                  ubar_stokes(i+1,j)+             &
#  endif
     &                                  ubar(i+1,j  ,krhs))*            &
     &                                 dmde(i,j))
#  endif
            ad_cff=0.0_r8
          END DO
        END DO
# endif
# ifdef UV_ADV
!
!-----------------------------------------------------------------------
!  Add in adjoint horizontal advection of momentum.
!-----------------------------------------------------------------------
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  if defined DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2hadv)=-fac
#  endif
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
!>
            ad_fac=ad_fac-ad_rhs_vbar(i,j)
!>          tl_fac=(tl_VFx(i+1,j)-tl_VFx(i,j))+                         &
!>   &             (tl_VFe(i,j)-tl_VFe(i,j-1))
!>
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-ad_fac
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+ad_fac
            ad_VFe(i,j-1)=ad_VFe(i,j-1)-ad_fac
            ad_VFe(i,j  )=ad_VFe(i,j  )+ad_fac
            ad_fac=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  if defined DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2hadv)=-fac
#  endif
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
!>
            ad_fac=ad_fac-ad_rhs_ubar(i,j)
!>          tl_fac=(tl_UFx(i,j)-tl_UFx(i-1,j))+                        &
!>   &             (tl_UFe(i,j+1)-tl_UFe(i,j))
!>
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-ad_fac
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+ad_fac
            ad_UFe(i,j  )=ad_UFe(i,j  )-ad_fac
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+ad_fac
            ad_fac=0.0_r8
          END DO
        END DO
#  ifdef UV_C2ADVECTION
!
!  Second-order, centered differences advection.
!
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
!>          tl_VFe(i,j)=0.25_r8*                                        &
!>   &                  ((tl_DVom(i,j)+tl_DVom(i,j+1))*                 &
!>   &                   (vbar(i,j  ,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
!>   &                    vbar_stokes(i,j  )+                           &
!>   &                    vbar_stokes(i,j+1)+                           &
#   endif
!>   &                    vbar(i,j+1,krhs))+                            &
!>   &                   (DVom(i,j)+DVom(i,j+1))*                       &
!>   &                   (tl_vbar(i,j  ,krhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_vbar_stokes(i,j  )+                        &
!>   &                    tl_vbar_stokes(i,j+1)+                        &
#   endif
!>   &                    tl_vbar(i,j+1,krhs)))
!>
            adfac=0.25_r8*ad_VFe(i,j)
            adfac1=adfac*(vbar(i,j  ,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                    vbar_stokes(i,j  )+                           &
     &                    vbar_stokes(i,j+1)+                           &
#   endif
     &                    vbar(i,j+1,krhs))
            adfac2=adfac*(DVom(i,j)+DVom(i,j+1))
            ad_DVom(i,j  )=ad_DVom(i,j  )+adfac1
            ad_DVom(i,j+1)=ad_DVom(i,j+1)+adfac1
            ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+adfac2
            ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+adfac2
#   ifdef NEARSHORE_MELLOR
            ad_vbar_stokes(i,j  )=ad_vbar_stokes(i,j  )+adfac2
            ad_vbar_stokes(i,j+1)=ad_vbar_stokes(i,j+1)+adfac2
#   endif
            ad_VFe(i,j)=0.0_r8
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
!>          tl_VFx(i,j)=0.25_r8*                                        &
!>   &                  ((tl_DUon(i,j)+tl_DUon(i,j-1))*                 &
!>   &                   (vbar(i  ,j,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
!>   &                    vbar_stokes(i  ,j)+                           &
!>   &                    vbar_stokes(i-1,j)+                           &
#   endif
!>   &                    vbar(i-1,j,krhs))+                            &
!>   &                   (DUon(i,j)+DUon(i,j-1))*                       &
!>   &                   (tl_vbar(i  ,j,krhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_vbar_stokes(i  ,j)+                        &
!>   &                    tl_vbar_stokes(i-1,j)+                        &
#   endif
!>   &                    tl_vbar(i-1,j,krhs)))
!>
            adfac=0.25_r8*ad_VFx(i,j)
            adfac1=adfac*(vbar(i  ,j,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                    vbar_stokes(i  ,j)+                           &
     &                    vbar_stokes(i-1,j)+                           &
#   endif
     &                    vbar(i-1,j,krhs))
            adfac2=adfac*(DUon(i,j)+DUon(i,j-1))
            ad_DUon(i,j  )=ad_DUon(i,j  )+adfac1
            ad_DUon(i,j-1)=ad_DUon(i,j-1)+adfac1
            ad_vbar(i  ,j,krhs)=ad_vbar(i  ,j,krhs)+adfac2
            ad_vbar(i-1,j,krhs)=ad_vbar(i-1,j,krhs)+adfac2
#   ifdef NEARSHORE_MELLOR
            ad_vbar_stokes(i-1,j)=ad_vbar_stokes(i-1,j)+adfac2
            ad_vbar_stokes(i  ,j)=ad_vbar_stokes(i  ,j)+adfac2
#   endif
            ad_VFx(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
!>          tl_UFe(i,j)=0.25_r8*                                        &
!>   &                  ((tl_DVom(i,j)+tl_DVom(i-1,j))*                 &
!>   &                   (ubar(i,j  ,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
!>   &                    ubar_stokes(i,j  )+                           &
!>   &                    ubar_stokes(i,j-1)+                           &
#   endif
!>   &                    ubar(i,j-1,krhs))+                            &
!>   &                   (DVom(i,j)+DVom(i-1,j))*                       &
!>   &                   (tl_ubar(i,j  ,krhs)+
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_ubar_stokes(i,j  )+                        &
!>   &                    tl_ubar_stokes(i,j-1)+                        &
#   endif
!>   &                    tl_ubar(i,j-1,krhs)))
!>
            adfac=0.25_r8*ad_UFe(i,j)
            adfac1=adfac*(ubar(i,j  ,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                    ubar_stokes(i,j  )+                           &
     &                    ubar_stokes(i,j-1)+                           &
#   endif
     &                    ubar(i,j-1,krhs))
            adfac2=adfac*(DVom(i,j)+DVom(i-1,j))
            ad_DVom(i  ,j)=ad_DVom(i  ,j)+adfac1
            ad_DVom(i-1,j)=ad_DVom(i-1,j)+adfac1
            ad_ubar(i,j  ,krhs)=ad_ubar(i,j  ,krhs)+adfac2
            ad_ubar(i,j-1,krhs)=ad_ubar(i,j-1,krhs)+adfac2
#   ifdef NEARSHORE_MELLOR
            ad_ubar_stokes(i,j-1)=ad_ubar_stokes(i,j-1)+adfac2
            ad_ubar_stokes(i,j  )=ad_ubar_stokes(i,j  )+adfac2
#   endif
            ad_UFe(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
!>          tl_UFx(i,j)=0.25_r8*                                        &
!>   &                  ((tl_DUon(i,j)+tl_DUon(i+1,j))*                 &
!>   &                   (ubar(i  ,j,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
!>   &                    ubar_stokes(i  ,j)+                           &
!>   &                    ubar_stokes(i+1,j)+                           &
#   endif
!>   &                    ubar(i+1,j,krhs))+                            &
!>   &                   (DUon(i,j)+DUon(i+1,j))*                       &
!>   &                   (tl_ubar(i  ,j,krhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_ubar_stokes(i  ,j)+                        &
!>   &                    tl_ubar_stokes(i+1,j)+                        &
#   endif
!>   &                    tl_ubar(i+1,j,krhs)))
!>
            adfac=0.25_r8*ad_UFx(i,j)
            adfac1=adfac*(ubar(i  ,j,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                    ubar_stokes(i  ,j)+                           &
     &                    ubar_stokes(i+1,j)+                           &
#   endif
     &                    ubar(i+1,j,krhs))
            adfac2=adfac*(DUon(i,j)+DUon(i+1,j))
            ad_DUon(i  ,j)=ad_DUon(i  ,j)+adfac1
            ad_DUon(i+1,j)=ad_DUon(i+1,j)+adfac1
            ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)+adfac2
            ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+adfac2
#   ifdef NEARSHORE_MELLOR
            ad_ubar_stokes(i  ,j)=ad_ubar_stokes(i  ,j)+adfac2
            ad_ubar_stokes(i+1,j)=ad_ubar_stokes(i+1,j)+adfac2
#   endif
            ad_UFx(i,j)=0.0_r8
          END DO
        END DO
#  else
!
!  Fourth-order, centered differences advection.
!
#   ifdef NS_PERIODIC
#    define JV_RANGE JstrV-1,Jend+1
#   else
#    define JV_RANGE MAX(JstrV-1,2),MIN(Jend+1,Mm(ng))
#   endif
        DO j=JV_RANGE
          DO i=Istr,Iend
            grad (i,j)=vbar(i,j-1,krhs)-2.0_r8*vbar(i,j,krhs)+          &
#   ifdef NEARSHORE_MELLOR
     &                 vbar_stokes(i,j-1)-2.0_r8*vbar_stokes(i,j)+      &
     &                 vbar_stokes(i,j+1)+                              &
#   endif
     &                 vbar(i,j+1,krhs)
            Dgrad(i,j)=DVom(i,j-1)-2.0_r8*DVom(i,j)+DVom(i,j+1)
          END DO
        END DO
#   undef JV_RANGE
#   ifndef NS_PERIODIC
        IF (NORTHERN_EDGE) THEN
          DO i=Istr,Iend
            grad (i,Jend+1)=grad (i,Jend)
            Dgrad(i,Jend+1)=Dgrad(i,Jend)
          END DO
        END IF
        IF (SOUTHERN_EDGE) THEN
          DO i=Istr,Iend
            grad (i,Jstr)=grad (i,Jstr+1)
            Dgrad(i,Jstr)=Dgrad(i,Jstr+1)
          END DO
        END IF
#   endif
        cff=1.0_r8/6.0_r8
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
!>          tl_VFe(i,j)=0.25_r8*                                        &
!>   &                  ((tl_vbar(i,j  ,krhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_vbar_stokes(i,j  )+                        &
!>   &                    tl_vbar_stokes(i,j+1)+                        &
#   endif
!>   &                    tl_vbar(i,j+1,krhs)-                          &
!>   &                    cff*(tl_grad (i,j)+tl_grad (i,j+1)))*         &
!>   &                   (DVom(i,j)+DVom(i,j+1)-                        &
!>   &                    cff*(Dgrad(i,j)+Dgrad(i,j+1)))+               &
!>   &                   (vbar(i,j  ,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
!>   &                    vbar_stokes(i,j  )+                           &
!>   &                    vbar_stokes(i,j+1)+                           &
#   endif
!>   &                    vbar(i,j+1,krhs)-                             &
!>   &                    cff*(grad (i,j)+grad (i,j+1)))*               &
!>   &                   (tl_DVom(i,j)+tl_DVom(i,j+1)-                  &
!>   &                    cff*(tl_Dgrad(i,j)+tl_Dgrad(i,j+1))))
!>
            adfac=0.25_r8*ad_VFe(i,j)
            adfac1=adfac*(DVom(i,j)+DVom(i,j+1)-                        &
     &                    cff*(Dgrad(i,j)+Dgrad(i,j+1)))
            adfac2=adfac1*cff
            adfac3=adfac*(vbar(i,j  ,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                    vbar_stokes(i,j  )+                           &
     &                    vbar_stokes(i,j+1)+                           &
#   endif
     &                    vbar(i,j+1,krhs)-                             &
     &                    cff*(grad (i,j)+grad (i,j+1)))
            adfac4=adfac3*cff
            ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+adfac1
            ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+adfac1
#   ifdef NEARSHORE_MELLOR
            ad_vbar_stokes(i,j  )=ad_vbar_stokes(i,j  )+adfac1
            ad_vbar_stokes(i,j+1)=ad_vbar_stokes(i,j+1)+adfac1
#   endif
            ad_grad (i,j  )=ad_grad (i,j  )-adfac2
            ad_grad (i,j+1)=ad_grad (i,j+1)-adfac2
            ad_DVom(i,j  )=ad_DVom(i,j  )+adfac3
            ad_DVom(i,j+1)=ad_DVom(i,j+1)+adfac3
            ad_Dgrad(i,j  )=ad_Dgrad(i,j  )-adfac4
            ad_Dgrad(i,j+1)=ad_Dgrad(i,j+1)-adfac4
            ad_VFe(i,j)=0.0_r8
          END DO
        END DO
#   ifndef NS_PERIODIC
        IF (NORTHERN_EDGE) THEN
          DO i=Istr,Iend
!>          tl_Dgrad(i,Jend+1)=tl_Dgrad(i,Jend)
!>
            ad_Dgrad(i,Jend)=ad_Dgrad(i,Jend)+ad_Dgrad(i,Jend+1)
            ad_Dgrad(i,Jend+1)=0.0_r8
!>          tl_grad (i,Jend+1)=tl_grad (i,Jend)
!>
            ad_grad (i,Jend)=ad_grad (i,Jend)+ad_grad (i,Jend+1)
            ad_grad (i,Jend+1)=0.0_r8
          END DO
        END IF
        IF (SOUTHERN_EDGE) THEN
          DO i=Istr,Iend
!>          tl_Dgrad(i,Jstr)=tl_Dgrad(i,Jstr+1)
!>
            ad_Dgrad(i,Jstr+1)=ad_Dgrad(i,Jstr+1)+ad_Dgrad(i,Jstr)
            ad_Dgrad(i,Jstr)=0.0_r8
!>          tl_grad (i,Jstr)=tl_grad (i,Jstr+1)
!>
            ad_grad (i,Jstr+1)=ad_grad (i,Jstr+1)+ad_grad (i,Jstr)
            ad_grad (i,Jstr)=0.0_r8
          END DO
        END IF
#   endif
#   ifdef NS_PERIODIC
#    define JV_RANGE JstrV-1,Jend+1
#   else
#    define JV_RANGE MAX(JstrV-1,2),MIN(Jend+1,Mm(ng))
#   endif
        DO j=JV_RANGE
          DO i=Istr,Iend
!>          tl_Dgrad(i,j)=tl_DVom(i,j-1)-2.0_r8*tl_DVom(i,j)+           &
!>   &                    tl_DVom(i,j+1)
!>
            ad_DVom(i,j-1)=ad_DVom(i,j-1)+ad_Dgrad(i,j)
            ad_DVom(i,j  )=ad_DVom(i,j  )-2.0_r8*ad_Dgrad(i,j)
            ad_DVom(i,j+1)=ad_DVom(i,j+1)+ad_Dgrad(i,j)
            ad_Dgrad(i,j)=0.0_r8
!>          tl_grad (i,j)=tl_vbar(i,j-1,krhs)-2.0_r8*tl_vbar(i,j,krhs)+ &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_vbar_stokes(i,j-1)-                        &
!>   &                    2.0_r8*tl_vbar_stokes(i,j)+                   &
!>   &                    tl_vbar_stokes(i,j+1)+                        &
#   endif
!>   &                    tl_vbar(i,j+1,krhs)
!>
            ad_vbar(i,j-1,krhs)=ad_vbar(i,j-1,krhs)+ad_grad(i,j)
            ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)-                    &
     &                          2.0_r8*ad_grad(i,j)
            ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+ad_grad(i,j)
#   ifdef NEARSHORE_MELLOR
            ad_vbar_stokes(i,j-1)=ad_vbar_stokes(i,j-1)+ad_grad(i,j)
            ad_vbar_stokes(i,j  )=ad_vbar_stokes(i,j  )-                &
     &                            2.0_r8*ad_grad(i,j)
            ad_vbar_stokes(i,j+1)=ad_vbar_stokes(i,j+1)+ad_grad(i,j)
#   endif
            ad_grad(i,j)=0.0_r8
          END DO
        END DO
#   undef JV_RANGE
#   ifdef EW_PERIODIC
#    define IV_RANGE Istr-1,Iend+1
#   else
#    define IV_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#   endif
        DO j=JstrV,Jend
          DO i=IV_RANGE
            grad(i,j)=vbar(i-1,j,krhs)-2.0_r8*vbar(i,j,krhs)+           &
#   ifdef NEARSHORE_MELLOR
     &                vbar_stokes(i-1,j)-2.0_r8*vbar_stokes(i,j)+       &
     &                vbar_stokes(i+1,j)+                               &
#   endif
     &                vbar(i+1,j,krhs)
          END DO
        END DO
#   undef IV_RANGE
#   ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=JstrV,Jend
            grad(Istr-1,j)=grad(Istr,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=JstrV,Jend
            grad(Iend+1,j)=grad(Iend,j)
          END DO
        END IF
#   endif
        DO j=JstrV-1,Jend
          DO i=Istr,Iend+1
            Dgrad(i,j)=DUon(i,j-1)-2.0_r8*DUon(i,j)+DUon(i,j+1)
          END DO
        END DO
        cff=1.0_r8/6.0_r8
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
!>          tl_VFx(i,j)=0.25_r8*                                        &
!>   &                  ((tl_vbar(i  ,j,krhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_vbar_stokes(i  ,j)+                        &
!>   &                    tl_vbar_stokes(i-1,j)+                        &
#   endif
!>   &                    tl_vbar(i-1,j,krhs)-                          &
!>   &                    cff*(tl_grad (i,j)+tl_grad (i-1,j)))*         &
!>   &                   (DUon(i,j)+DUon(i,j-1)-                        &
!>   &                    cff*(Dgrad(i,j)+Dgrad(i,j-1)))+               &
!>   &                   (vbar(i  ,j,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
!>   &                    vbar_stokes(i  ,j)+                           &
!>   &                    vbar_stokes(i-1,j)+                           &
#   endif
!>   &                    vbar(i-1,j,krhs)-                             &
!>   &                    cff*(grad (i,j)+grad (i-1,j)))*               &
!>   &                   (tl_DUon(i,j)+tl_DUon(i,j-1)-                  &
!>   &                    cff*(tl_Dgrad(i,j)+tl_Dgrad(i,j-1))))
!>
            adfac=0.25_r8*ad_VFx(i,j)
            adfac1=adfac*(DUon(i,j)+DUon(i,j-1)-                        &
     &                    cff*(Dgrad(i,j)+Dgrad(i,j-1)))
            adfac2=adfac1*cff
            adfac3=adfac*(vbar(i  ,j,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                    vbar_stokes(i  ,j)+                           &
     &                    vbar_stokes(i-1,j)+                           &
#   endif
     &                    vbar(i-1,j,krhs)-                             &
     &                    cff*(grad (i,j)+grad (i-1,j)))
            adfac4=adfac3*cff
            ad_vbar(i-1,j,krhs)=ad_vbar(i-1,j,krhs)+adfac1
            ad_vbar(i  ,j,krhs)=ad_vbar(i  ,j,krhs)+adfac1
#   ifdef NEARSHORE_MELLOR
            ad_vbar_stokes(i-1,j)=ad_vbar_stokes(i-1,j)+adfac1
            ad_vbar_stokes(i  ,j)=ad_vbar_stokes(i  ,j)+adfac1
#   endif
            ad_grad (i-1,j)=ad_grad (i-1,j)-adfac2
            ad_grad (i  ,j)=ad_grad (i  ,j)-adfac2
            ad_DUon(i,j-1)=ad_DUon(i,j-1)+adfac3
            ad_DUon(i,j  )=ad_DUon(i,j  )+adfac3
            ad_Dgrad(i,j-1)=ad_Dgrad(i,j-1)-adfac4
            ad_Dgrad(i,j  )=ad_Dgrad(i,j  )-adfac4
            ad_VFx(i,j)=0.0_r8
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=Istr,Iend+1
!>          tl_Dgrad(i,j)=tl_DUon(i,j-1)-2.0_r8*tl_DUon(i,j)+           &
!>   &                    tl_DUon(i,j+1)
!>
            ad_DUon(i,j-1)=ad_DUon(i,j-1)+ad_Dgrad(i,j)
            ad_DUon(i,j  )=ad_DUon(i,j  )-2.0_r8*ad_Dgrad(i,j)
            ad_DUon(i,j+1)=ad_DUon(i,j+1)+ad_Dgrad(i,j)
            ad_Dgrad(i,j)=0.0_r8
          END DO
        END DO
#   ifndef EW_PERIODIC
        IF (EASTERN_EDGE) THEN
          DO j=JstrV,Jend
!>          tl_grad(Iend+1,j)=tl_grad(Iend,j)
!>
            ad_grad(Iend,j)=ad_grad(Iend,j)+ad_grad(Iend+1,j)
            ad_grad(Iend+1,j)=0.0_r8
          END DO
        END IF
        IF (WESTERN_EDGE) THEN
          DO j=JstrV,Jend
!>          tl_grad(Istr-1,j)=tl_grad(Istr,j)
!>
            ad_grad(Istr,j)=ad_grad(Istr,j)+ad_grad(Istr-1,j)
            ad_grad(Istr-1,j)=0.0_r8
          END DO
        END IF
#   endif
#   ifdef EW_PERIODIC
#    define IV_RANGE Istr-1,Iend+1
#   else
#    define IV_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#   endif
        DO j=JstrV,Jend
          DO i=IV_RANGE
!>          tl_grad(i,j)=tl_vbar(i-1,j,krhs)-2.0_r8*tl_vbar(i,j,krhs)+  &
#   ifdef NEARSHORE_MELLOR
!>   &                   tl_vbar_stokes(i-1,j)-                         &
!>   &                   2.0_r8*tl_vbar_stokes(i,j)+                    &
!>   &                   tl_vbar_stokes(i+1,j)+                         &
#   endif
!>   &                   tl_vbar(i+1,j,krhs)
!>
            ad_vbar(i-1,j,krhs)=ad_vbar(i-1,j,krhs)+ad_grad(i,j)
            ad_vbar(i  ,j,krhs)=ad_vbar(i  ,j,krhs)-                    &
     &                          2.0_r8*ad_grad(i,j)
            ad_vbar(i+1,j,krhs)=ad_vbar(i+1,j,krhs)+ad_grad(i,j)
#   ifdef NEARSHORE_MELLOR
            ad_vbar_stokes(i-1,j)=ad_vbar_stokes(i-1,j)+ad_grad(i,j)
            ad_vbar_stokes(i  ,j)=ad_vbar_stokes(i  ,j)-                &
     &                            2.0_r8*ad_grad(i,j)
            ad_vbar_stokes(i+1,j)=ad_vbar_stokes(i+1,j)+ad_grad(i,j)
#   endif
            ad_grad(i,j)=0.0_r8
          END DO
        END DO
#   undef IV_RANGE
#   ifdef NS_PERIODIC
#    define JU_RANGE Jstr-1,Jend+1
#   else
#    define JU_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#   endif
        DO j=JU_RANGE
          DO i=IstrU,Iend
            grad(i,j)=ubar(i,j-1,krhs)-2.0_r8*ubar(i,j,krhs)+           &
#   ifdef NEARSHORE_MELLOR
     &                ubar_stokes(i,j-1)-2.0_r8*ubar_stokes(i,j)+       &
     &                ubar_stokes(i,j+1)+                               &
#   endif
     &                ubar(i,j+1,krhs)
          END DO
        END DO
#   undef JU_RANGE
#   ifndef NS_PERIODIC
        IF (SOUTHERN_EDGE) THEN
          DO i=IstrU,Iend
            grad(i,Jstr-1)=grad(i,Jstr)
          END DO
        END IF
        IF (NORTHERN_EDGE) THEN
          DO i=IstrU,Iend
            grad(i,Jend+1)=grad(i,Jend)
          END DO
        END IF
#   endif
        DO j=Jstr,Jend+1
          DO i=IstrU-1,Iend
            Dgrad(i,j)=DVom(i-1,j)-2.0_r8*DVom(i,j)+DVom(i+1,j)
          END DO
        END DO
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
!>          tl_UFe(i,j)=0.25_r8*                                        &
!>   &                  ((tl_ubar(i,j  ,krhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_ubar_stokes(i,j  )+                        &
!>   &                    tl_ubar_stokes(i,j-1)+                        &
#   endif
!>   &                    tl_ubar(i,j-1,krhs)-                          &
!>   &                    cff*(tl_grad (i,j)+tl_grad (i,j-1)))*         &
!>   &                   (DVom(i,j)+DVom(i-1,j)-                        &
!>   &                    cff*(Dgrad(i,j)+Dgrad(i-1,j)))+               &
!>   &                   (ubar(i,j  ,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
!>   &                    ubar_stokes(i,j  )+                           &
!>   &                    ubar_stokes(i,j-1)+                           &
#   endif
!>   &                    ubar(i,j-1,krhs)-                             &
!>   &                    cff*(grad (i,j)+grad (i,j-1)))*               &
!>   &                   (tl_DVom(i,j)+tl_DVom(i-1,j)-                  &
!>   &                    cff*(tl_Dgrad(i,j)+tl_Dgrad(i-1,j))))
!>
            adfac=0.25_r8*ad_UFe(i,j)
            adfac1=adfac*(DVom(i,j)+DVom(i-1,j)-                        &
     &                    cff*(Dgrad(i,j)+Dgrad(i-1,j)))
            adfac2=adfac1*cff
            adfac3=adfac*(ubar(i,j  ,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                    ubar_stokes(i,j  )+                           &
     &                    ubar_stokes(i,j-1)+                           &
#   endif
     &                    ubar(i,j-1,krhs)-                             &
     &                    cff*(grad (i,j)+grad (i,j-1)))
            adfac4=adfac3*cff
            ad_ubar(i,j-1,krhs)=ad_ubar(i,j-1,krhs)+adfac1
            ad_ubar(i,j  ,krhs)=ad_ubar(i,j  ,krhs)+adfac1
#   ifdef NEARSHORE_MELLOR
            ad_ubar_stokes(i,j-1)=ad_ubar_stokes(i,j-1)+adfac1
            ad_ubar_stokes(i,j  )=ad_ubar_stokes(i,j  )+adfac1
#   endif
            ad_grad (i,j-1)=ad_grad (i,j-1)-adfac2
            ad_grad (i,j  )=ad_grad (i,j  )-adfac2
            ad_DVom(i-1,j)=ad_DVom(i-1,j)+adfac3
            ad_DVom(i  ,j)=ad_DVom(i  ,j)+adfac3
            ad_Dgrad(i-1,j)=ad_Dgrad(i-1,j)-adfac4
            ad_Dgrad(i  ,j)=ad_Dgrad(i  ,j)-adfac4
            ad_UFe(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=IstrU-1,Iend
!>          tl_Dgrad(i,j)=tl_DVom(i-1,j)-2.0_r8*tl_DVom(i,j)+           &
!>   &                    tl_DVom(i+1,j)
!>
            ad_DVom(i-1,j)=ad_DVom(i-1,j)+ad_Dgrad(i,j)
            ad_DVom(i  ,j)=ad_DVom(i  ,j)-2.0_r8*ad_Dgrad(i,j)
            ad_DVom(i+1,j)=ad_DVom(i+1,j)+ad_Dgrad(i,j)
            ad_Dgrad(i,j)=0.0_r8
          END DO
        END DO
#   ifndef NS_PERIODIC
        IF (NORTHERN_EDGE) THEN
          DO i=IstrU,Iend
!>          tl_grad(i,Jend+1)=tl_grad(i,Jend)
!>
            ad_grad(i,Jend)=ad_grad(i,Jend)+ad_grad(i,Jend+1)
            ad_grad(i,Jend+1)=0.0_r8
          END DO
        END IF
        IF (SOUTHERN_EDGE) THEN
          DO i=IstrU,Iend
!>          tl_grad(i,Jstr-1)=tl_grad(i,Jstr)
!>
            ad_grad(i,Jstr)=ad_grad(i,Jstr)+ad_grad(i,Jstr-1)
            ad_grad(i,Jstr-1)=0.0_r8
          END DO
        END IF
#   endif
#   ifdef NS_PERIODIC
#    define JU_RANGE Jstr-1,Jend+1
#   else
#    define JU_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#   endif
        DO j=JU_RANGE
          DO i=IstrU,Iend
!>          tl_grad(i,j)=tl_ubar(i,j-1,krhs)-2.0_r8*tl_ubar(i,j,krhs)+  &
#   ifdef NEARSHORE_MELLOR
!>   &                   tl_ubar_stokes(i,j-1)-                         &
!>   &                   2.0_r8*tl_ubar_stokes(i,j)+                    &
!>   &                   tl_ubar_stokes(i,j+1)+                         &
#   endif
!>   &                   tl_ubar(i,j+1,krhs)
!>
            ad_ubar(i,j-1,krhs)=ad_ubar(i,j-1,krhs)+ad_grad(i,j)
            ad_ubar(i,j  ,krhs)=ad_ubar(i,j  ,krhs)-                    &
     &                          2.0_r8*ad_grad(i,j)
            ad_ubar(i,j+1,krhs)=ad_ubar(i,j+1,krhs)+ad_grad(i,j)
#   ifdef NEARSHORE_MELLOR
            ad_ubar_stokes(i,j-1)=ad_ubar_stokes(i,j-1)+ad_grad(i,j)
            ad_ubar_stokes(i,j  )=ad_ubar_stokes(i,j)-                  &
     &                            2.0_r8*ad_grad(i,j)
            ad_ubar_stokes(i,j+1)=ad_ubar_stokes(i,j+1)+ad_grad(i,j)
#   endif
            ad_grad(i,j)=0.0_r8
          END DO
        END DO
#   undef JU_RANGE
#   ifdef EW_PERIODIC
#    define IU_RANGE IstrU-1,Iend+1
#   else
#    define IU_RANGE MAX(IstrU-1,2),MIN(Iend+1,Lm(ng))
#   endif
        DO j=Jstr,Jend
          DO i=IU_RANGE
            grad (i,j)=ubar(i-1,j,krhs)-2.0_r8*ubar(i,j,krhs)+          &
#    ifdef NEARSHORE_MELLOR
     &                 ubar_stokes(i-1,j)-2.0_r8*ubar_stokes(i,j)+      &
     &                 ubar_stokes(i+1,j)+                              &
#    endif
     &                 ubar(i+1,j,krhs)
            Dgrad(i,j)=DUon(i-1,j)-2.0_r8*DUon(i,j)+DUon(i+1,j)
          END DO
        END DO
#   undef IU_RANGE
#   ifndef EW_PERIODIC
        IF (WESTERN_EDGE) THEN
          DO j=Jstr,Jend
            grad (Istr,j)=grad (Istr+1,j)
            Dgrad(Istr,j)=Dgrad(Istr+1,j)
          END DO
        END IF
        IF (EASTERN_EDGE) THEN
          DO j=Jstr,Jend
            grad (Iend+1,j)=grad (Iend,j)
            Dgrad(Iend+1,j)=Dgrad(Iend,j)
          END DO
        END IF
#   endif
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
!>          tl_UFx(i,j)=0.25_r8*                                        &
!>   &                  ((ubar(i  ,j,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
!>   &                    ubar_stokes(i  ,j)+                           &
!>   &                    ubar_stokes(i+1,j)+                           &
#   endif
!>   &                    ubar(i+1,j,krhs)-                             &
!>   &                    cff*(grad (i,j)+grad (i+1,j)))*               &
!>   &                   (tl_DUon(i,j)+tl_DUon(i+1,j)-                  &
!>   &                    cff*(tl_Dgrad(i,j)+tl_Dgrad(i+1,j)))+         &
!>   &                   (tl_ubar(i  ,j,krhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_ubar_stokes(i  ,j)+                        &
!>   &                    tl_ubar_stokes(i+1,j)+                        &
#   endif
!>   &                    tl_ubar(i+1,j,krhs)-                          &
!>   &                    cff*(tl_grad (i,j)+tl_grad (i+1,j)))*         &
!>   &                   (DUon(i,j)+DUon(i+1,j)-                        &
!>   &                    cff*(Dgrad(i,j)+Dgrad(i+1,j))))
!>
            adfac=0.25_r8*ad_UFx(i,j)
            adfac1=adfac*(DUon(i,j)+DUon(i+1,j)-                        &
     &                    cff*(Dgrad(i,j)+Dgrad(i+1,j)))
            adfac2=adfac1*cff
            adfac3=adfac*(ubar(i  ,j,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                    ubar_stokes(i  ,j)+                           &
     &                    ubar_stokes(i+1,j)+                           &
#   endif
     &                    ubar(i+1,j,krhs)-                             &
     &                    cff*(grad (i,j)+grad (i+1,j)))
            adfac4=adfac3*cff
            ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)+adfac1
            ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+adfac1
#   ifdef NEARSHORE_MELLOR
            ad_ubar_stokes(i  ,j)=ad_ubar_stokes(i  ,j)+adfac1
            ad_ubar_stokes(i+1,j)=ad_ubar_stokes(i+1,j)+adfac1
#   endif
            ad_grad (i  ,j)=ad_grad (i  ,j)-adfac2
            ad_grad (i+1,j)=ad_grad (i+1,j)-adfac2
            ad_DUon(i  ,j)=ad_DUon(i  ,j)+adfac3
            ad_DUon(i+1,j)=ad_DUon(i+1,j)+adfac3
            ad_Dgrad(i  ,j)=ad_Dgrad(i  ,j)-adfac4
            ad_Dgrad(i+1,j)=ad_Dgrad(i+1,j)-adfac4
            ad_UFx(i,j)=0.0_r8
          END DO
        END DO
#   ifndef EW_PERIODIC
        IF (EASTERN_EDGE) THEN
          DO j=Jstr,Jend
!>          tl_Dgrad(Iend+1,j)=tl_Dgrad(Iend,j)
!>
            ad_Dgrad(Iend,j)=ad_Dgrad(Iend,j)+ad_Dgrad(Iend+1,j)
            ad_Dgrad(Iend+1,j)=0.0_r8
!>          tl_grad (Iend+1,j)=tl_grad (Iend,j)
!>
            ad_grad (Iend,j)=ad_grad (Iend,j)+ad_grad (Iend+1,j)
            ad_grad (Iend+1,j)=0.0_r8
          END DO
        END IF
        IF (WESTERN_EDGE) THEN
          DO j=Jstr,Jend
!>          tl_Dgrad(Istr,j)=tl_Dgrad(Istr+1,j)
!>
            ad_Dgrad(Istr+1,j)=ad_Dgrad(Istr+1,j)+ad_Dgrad(Istr,j)
            ad_Dgrad(Istr,j)=0.0_r8
!>          tl_grad (Istr,j)=tl_grad (Istr+1,j)
!>
            ad_grad (Istr+1,j)=ad_grad (Istr+1,j)+ad_grad (Istr,j)
            ad_grad (Istr,j)=0.0_r8
          END DO
        END IF
#   endif
#   ifdef EW_PERIODIC
#    define IU_RANGE IstrU-1,Iend+1
#   else
#    define IU_RANGE MAX(IstrU-1,2),MIN(Iend+1,Lm(ng))
#   endif
        DO j=Jstr,Jend
          DO i=IU_RANGE
!>          tl_Dgrad(i,j)=tl_DUon(i-1,j)-2.0_r8*tl_DUon(i,j)+           &
!>   &                    tl_DUon(i+1,j)
!>
            ad_DUon(i-1,j)=ad_DUon(i-1,j)+ad_Dgrad(i,j)
            ad_DUon(i  ,j)=ad_DUon(i  ,j)-2.0_r8*ad_Dgrad(i,j)
            ad_DUon(i+1,j)=ad_DUon(i+1,j)+ad_Dgrad(i,j)
            ad_Dgrad(i,j)=0.0_r8
!>          tl_grad (i,j)=tl_ubar(i-1,j,krhs)-2.0_r8*tl_ubar(i,j,krhs)+ &
#   ifdef NEARHSORE_MELLOR
!>   &                    tl_ubar_stokes(i-1,j)-                        &
!>   &                    2.0_r8*tl_ubar_stokes(i,j)+                   &
!>   &                    tl_ubar_stokes(i+1,j)+                        &
#   endif
!>   &                    tl_ubar(i+1,j,krhs)
!>
            ad_ubar(i-1,j,krhs)=ad_ubar(i-1,j,krhs)+ad_grad (i,j)
            ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)-                    &
     &                          2.0_r8*ad_grad (i,j)
            ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+ad_grad (i,j)
#   ifdef NEARHSORE_MELLOR
            ad_ubar_stokes(i-1,j)=ad_ubar_stokes(i-1,j)+ad_grad (i,j)
            ad_ubar_stokes(i  ,j)=ad_ubar_stokes(i  ,j)-                &
     &                            2.0_r8*ad_grad (i,j)
            ad_ubar_stokes(i+1,j)=ad_ubar_stokes(i+1,j)+ad_grad (i,j)
#   endif
            ad_grad(i,j)=0.0_r8
          END DO
        END DO
#   undef IU_RANGE
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute adjoint pressure gradient terms.
!-----------------------------------------------------------------------
!      
!  Compute BASIC STATE fields associated with pressure gradient and
!  time-stepping of adjoint free-surface.
!
        fac=1000.0_r8/rho0
        IF (FIRST_2D_STEP) THEN
          cff1=dtfast(ng)
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              rhs_zeta(i,j)=(DUon(i,j)-DUon(i+1,j))+                    &
     &                      (DVom(i,j)-DVom(i,j+1))
              zeta_new(i,j)=zeta(i,j,kstp)+                             &
     &                      pm(i,j)*pn(i,j)*cff1*rhs_zeta(i,j)
# ifdef MASKING
              zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
# endif
              zwrk(i,j)=0.5_r8*(zeta(i,j,kstp)+zeta_new(i,j))
# if defined VAR_RHO_2D && defined SOLVE3D
              gzeta(i,j)=(fac+rhoS(i,j))*zwrk(i,j)
              gzeta2(i,j)=gzeta(i,j)*zwrk(i,j)
              gzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
              gzeta(i,j)=zwrk(i,j)
              gzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
# endif
            END DO
          END DO
        ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
          cff1=2.0_r8*dtfast(ng)
          cff4=4.0_r8/25.0_r8
          cff5=1.0_r8-2.0_r8*cff4
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              rhs_zeta(i,j)=(DUon(i,j)-DUon(i+1,j))+                    &
     &                      (DVom(i,j)-DVom(i,j+1))
              zeta_new(i,j)=zeta(i,j,kstp)+                             &
     &                      pm(i,j)*pn(i,j)*cff1*rhs_zeta(i,j)
# ifdef MASKING
              zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
# endif
              zwrk(i,j)=cff5*zeta(i,j,krhs)+                            &
     &                  cff4*(zeta(i,j,kstp)+zeta_new(i,j))
# if defined VAR_RHO_2D && defined SOLVE3D
              gzeta(i,j)=(fac+rhoS(i,j))*zwrk(i,j)
              gzeta2(i,j)=gzeta(i,j)*zwrk(i,j)
              gzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
              gzeta(i,j)=zwrk(i,j)
              gzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
# endif
            END DO
          END DO
        ELSE IF (CORRECTOR_2D_STEP) THEN
          cff1=dtfast(ng)*5.0_r8/12.0_r8
          cff2=dtfast(ng)*8.0_r8/12.0_r8
          cff3=dtfast(ng)*1.0_r8/12.0_r8
          cff4=2.0_r8/5.0_r8
          cff5=1.0_r8-cff4
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff=cff1*((DUon(i,j)-DUon(i+1,j))+                        &
     &                  (DVom(i,j)-DVom(i,j+1)))
              zeta_new(i,j)=zeta(i,j,kstp)+                             &
     &                      pm(i,j)*pn(i,j)*(cff+                       &
     &                                       cff2*rzeta(i,j,kstp)-      &
     &                                       cff3*rzeta(i,j,ptsk))
# ifdef MASKING
              zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
# endif
              zwrk(i,j)=cff5*zeta_new(i,j)+cff4*zeta(i,j,krhs)
# if defined VAR_RHO_2D && defined SOLVE3D
              gzeta(i,j)=(fac+rhoS(i,j))*zwrk(i,j)
              gzeta2(i,j)=gzeta(i,j)*zwrk(i,j)
              gzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
              gzeta(i,j)=zwrk(i,j)
              gzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
# endif
            END DO
          END DO
        END IF
!
!  Compute adjoint pressure gradient.
!
        cff1=0.5_r8*g
        cff2=1.0_r8/3.0_r8
        DO j=Jstr,Jend
          IF (j.ge.JstrV) THEN
            DO i=Istr,Iend
# ifdef DIAGNOSTICS_UV
!!            DiaV2rhs(i,j,M2pgrd)=rhs_vbar(i,j)
# endif
!>            tl_rhs_vbar(i,j)=cff1*om_v(i,j)*                          &
!>   &                         ((tl_h(i,j-1)+                           &
!>   &                           tl_h(i,j  ))*                          &
!>   &                          (gzeta(i,j-1)-                          &
!>   &                           gzeta(i,j  ))+                         &
!>   &                          (h(i,j-1)+                              &
!>   &                           h(i,j  ))*                             &
!>   &                          (tl_gzeta(i,j-1)-                       &
!>   &                           tl_gzeta(i,j  ))+                      &
# if defined VAR_RHO_2D && defined SOLVE3D
!>   &                          (tl_h(i,j-1)-                           &
!>   &                           tl_h(i,j  ))*                          &
!>   &                          (gzetaSA(i,j-1)+                        &
!>   &                           gzetaSA(i,j  )+                        &
!>   &                           cff2*(rhoA(i,j-1)-                     &
!>   &                                 rhoA(i,j  ))*                    &
!>   &                                (zwrk(i,j-1)-                     &
!>   &                                 zwrk(i,j  )))+                   &
!>   &                          (h(i,j-1)-                              &
!>   &                           h(i,j  ))*                             &
!>   &                          (tl_gzetaSA(i,j-1)+                     &
!>   &                           tl_gzetaSA(i,j  )+                     &
!>   &                           cff2*((tl_rhoA(i,j-1)-                 &
!>   &                                  tl_rhoA(i,j  ))*                &
!>   &                                 (zwrk(i,j-1)-                    &
!>   &                                  zwrk(i,j  ))+                   &
!>   &                                 (rhoA(i,j-1)-                    &
!>   &                                  rhoA(i,j  ))*                   &
!>   &                                 (tl_zwrk(i,j-1)-                 &
!>   &                                  tl_zwrk(i,j  ))))+              &
# endif
!>   &                          (tl_gzeta2(i,j-1)-                      &
!>   &                           tl_gzeta2(i,j  )
!>
              adfac=cff1*om_v(i,j)*ad_rhs_vbar(i,j)
              adfac1=adfac*(gzeta(i,j-1)-gzeta(i,j  ))
              adfac2=adfac*(h(i,j-1)+h(i,j  ))
              ad_h(i,j-1)=ad_h(i,j-1)+adfac1
              ad_h(i,j  )=ad_h(i,j  )+adfac1
              ad_gzeta(i,j-1)=ad_gzeta(i,j-1)+adfac2
              ad_gzeta(i,j  )=ad_gzeta(i,j  )-adfac2
              ad_gzeta2(i,j-1)=ad_gzeta2(i,j-1)+adfac
              ad_gzeta2(i,j  )=ad_gzeta2(i,j  )-adfac
# if defined VAR_RHO_2D && defined SOLVE3D
              adfac1=adfac*(gzetaSA(i,j-1)+                             &
     &                      gzetaSA(i,j  )+                             &
     &                      cff2*(rhoA(i,j-1)-                          &
     &                            rhoA(i,j  ))*                         &
     &                           (zwrk(i,j-1)-                          &
     &                            zwrk(i,j  )))
              adfac2=adfac*(h(i,j-1)-h(i,j))
              adfac3=adfac2*cff2*(zwrk(i,j-1)-zwrk(i,j))
              adfac4=adfac2*cff2*(rhoA(i,j-1)-rhoA(i,j))
              ad_h(i,j-1)=ad_h(i,j-1)+adfac1
              ad_h(i,j  )=ad_h(i,j  )-adfac1
              ad_gzetaSA(i,j-1)=ad_gzetaSA(i,j-1)+adfac2
              ad_gzetaSA(i,j  )=ad_gzetaSA(i,j  )+adfac2
              ad_rhoA(i,j-1)=ad_rhoA(i,j-1)+adfac3
              ad_rhoA(i,j  )=ad_rhoA(i,j  )-adfac3
              ad_zwrk(i,j-1)=ad_zwrk(i,j-1)+adfac4
              ad_zwrk(i,j  )=ad_zwrk(i,j  )-adfac4
# endif
              ad_rhs_vbar(i,j)=0.0_r8
            END DO
          END IF
          DO i=IstrU,Iend
# ifdef DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2pgrd)=rhs_ubar(i,j)
# endif
!>          tl_rhs_ubar(i,j)=cff1*on_u(i,j)*                            &
!>   &                       ((tl_h(i-1,j)+                             &
!>   &                         tl_h(i ,j))*                             &
!>   &                        (gzeta(i-1,j)-                            &
!>   &                         gzeta(i  ,j))+                           &
!>   &                        (h(i-1,j)+                                &
!>   &                         h(i  ,j))*                               &
!>   &                        (tl_gzeta(i-1,j)-                         &
!>   &                         tl_gzeta(i  ,j))+                        &
# if defined VAR_RHO_2D && defined SOLVE3D
!>   &                        (tl_h(i-1,j)-                             &
!>   &                         tl_h(i  ,j))*                            &
!>   &                        (gzetaSA(i-1,j)+                          &
!>   &                         gzetaSA(i  ,j)+                          &
!>   &                         cff2*(rhoA(i-1,j)-                       &
!>   &                               rhoA(i  ,j))*                      &
!>   &                              (zwrk(i-1,j)-                       &
!>   &                               zwrk(i  ,j)))+                     &
!>   &                        (h(i-1,j)-                                &
!>   &                         h(i  ,j))*                               &
!>   &                        (tl_gzetaSA(i-1,j)+                       &
!>   &                         tl_gzetaSA(i  ,j)+                       &
!>   &                         cff2*((tl_rhoA(i-1,j)-                   &
!>   &                                tl_rhoA(i  ,j))*                  &
!>   &                               (zwrk(i-1,j)-                      &
!>   &                                zwrk(i  ,j))+                     &
!>   &                               (rhoA(i-1,j)-                      &
!>   &                                rhoA(i  ,j))*                     &
!>   &                               (tl_zwrk(i-1,j)-                   &
!>   &                                tl_zwrk(i  ,j))))+                &
# endif
!>   &                        (tl_gzeta2(i-1,j)-                        &
!>   &                         tl_gzeta2(i  ,j)))
!>
            adfac=cff1*on_u(i,j)*ad_rhs_ubar(i,j)
            adfac1=adfac*(gzeta(i-1,j)-gzeta(i  ,j))
            adfac2=adfac*(h(i-1,j)+h(i  ,j))
            ad_h(i-1,j)=ad_h(i-1,j)+adfac1
            ad_h(i  ,j)=ad_h(i  ,j)+adfac1
            ad_gzeta(i-1,j)=ad_gzeta(i-1,j)+adfac2
            ad_gzeta(i  ,j)=ad_gzeta(i  ,j)-adfac2
            ad_gzeta2(i-1,j)=ad_gzeta2(i-1,j)+adfac
            ad_gzeta2(i  ,j)=ad_gzeta2(i  ,j)-adfac
# if defined VAR_RHO_2D && defined SOLVE3D
            adfac1=adfac*(gzetaSA(i-1,j)+                               &
     &                    gzetaSA(i  ,j)+                               &
     &                    cff2*(rhoA(i-1,j)-                            &
     &                          rhoA(i  ,j))*                           &
     &                         (zwrk(i-1,j)-                            &
     &                          zwrk(i  ,j)))
            adfac2=adfac*(h(i-1,j)-h(i  ,j))
            adfac3=adfac2*cff2*(zwrk(i-1,j)-zwrk(i,j))
            adfac4=adfac2*cff2*(rhoA(i-1,j)-rhoA(i,j))
            ad_h(i-1,j)=ad_h(i-1,j)+adfac1
            ad_h(i  ,j)=ad_h(i  ,j)-adfac1
            ad_gzetaSA(i-1,j)=ad_gzetaSA(i-1,j)+adfac2
            ad_gzetaSA(i  ,j)=ad_gzetaSA(i  ,j)+adfac2
            ad_rhoA(i-1,j)=ad_rhoA(i-1,j)+adfac3
            ad_rhoA(i  ,j)=ad_rhoA(i  ,j)-adfac3
            ad_zwrk(i-1,j)=ad_zwrk(i-1,j)+adfac4
            ad_zwrk(i  ,j)=ad_zwrk(i  ,j)-adfac4
# endif
            ad_rhs_ubar(i,j)=0.0_r8
          END DO
        END DO
!
!  Set adjoint free-surface lateral boundary conditions.
!
# ifdef DISTRIBUTE
!>      CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints, EWperiodic, NSperiodic,       &
!>   &                      tl_zeta(:,:,knew))
!>
        CALL ad_mp_exchange2d (ng, tile, iADM, 1,                       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints, EWperiodic, NSperiodic,    &
     &                         ad_zeta(:,:,knew))
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_zeta(:,:,knew))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,knew))
# endif
!>      CALL tl_zetabc_tile (ng, tile,                                  &
!>   &                       LBi, UBi, LBj, UBj,                        &
!>   &                       IminS, ImaxS, JminS, JmaxS,                &
!>   &                       krhs, kstp, knew,                          &
!>   &                       zeta, tl_zeta)
!>
        CALL ad_zetabc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, knew,                          &
     &                       zeta, ad_zeta)
# ifdef Q_PSOURCE
!
!  Apply mass point sources - Volume influx.
!
      DO is=1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        IF (((IstrR.le.i).and.(i.le.IendR)).and.                        &
     &      ((JstrR.le.j).and.(j.le.JendR))) THEN
!>        tl_zeta(i,j,knew)=tl_zeta(i,j,knew)+0.0_r8
!>
        END IF
      END DO
# endif
!
!  If adjoint predictor step, load right-side-term into shared array.
!
        IF (PREDICTOR_2D_STEP(ng)) THEN
# ifdef DISTRIBUTE
!>        CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        NghostPoints, EWperiodic, NSperiodic,     &
!>   &                        tl_rzeta(:,:,krhs))
!>
          CALL ad_mp_exchange2d (ng, tile, iADM, 1,                     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints, EWperiodic, NSperiodic,  &
     &                           ad_rzeta(:,:,krhs))
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
!>        CALL exchange_r2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_rzeta(:,:,krhs))
!>
          CALL ad_exchange_r2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_rzeta(:,:,krhs))
# endif
          DO j=Jstr,Jend
            DO i=Istr,Iend
!>            tl_rzeta(i,j,krhs)=tl_rhs_zeta(i,j)
!>
              ad_rhs_zeta(i,j)=ad_rhs_zeta(i,j)+ad_rzeta(i,j,krhs)
              ad_rzeta(i,j,krhs)=0.0
            END DO
          END DO
        END IF

# ifndef SOLVE3D
!
!  Save free-surface adjoint solution for IO purposes.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta_sol(i,j)=ad_zeta(i,j,knew)
          END DO
        END DO
# endif
!
!  Load new adjoint free-surface values into shared array at both
!  predictor and corrector steps.
# ifdef WET_DRY_NOT_YET
!  Modify new free-surface to Ensure that depth is > Dcrit for masked
!  cells.
# endif
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
# if defined WET_DRY_NOT_YET && defined MASKING
!>          tl_zeta(i,j,knew)=tl_zeta(i,j,knew)-                          &
!>    &                       tl_h(i,j)*(1.0_r8-rmask(i,j))
!>
            ad_h(i,j)=ad_h(i,j)+(1.0_r8-rmask(i,j))*ad_zeta(i,j,knew)
# endif
!>          tl_zeta(i,j,knew)=tl_zeta_new(i,j)
!>
            ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_zeta(i,j,knew)
            ad_zeta(i,j,knew)=0.0_r8
          END DO
        END DO
!
!=======================================================================
!  Time step adjoint free-surface equation.
!=======================================================================
!
!  During the first time-step, the predictor step is Forward-Euler
!  and the corrector step is Backward-Euler. Otherwise, the predictor
!  step is Leap-frog and the corrector step is Adams-Moulton.
!
# if defined VAR_RHO_2D && defined SOLVE3D
        fac=1000.0_r8/rho0
# endif
        IF (FIRST_2D_STEP) THEN
          cff1=dtfast(ng)
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
# if defined VAR_RHO_2D && defined SOLVE3D
!>            tl_gzetaSA(i,j)=tl_zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))+       &
!>   &                        zwrk(i,j)*(tl_rhoS(i,j)-tl_rhoA(i,j))
!>
              adfac=zwrk(i,j)*ad_gzetaSA(i,j)
              ad_zwrk(i,j)=ad_zwrk(i,j)+                                &
     &                     (rhoS(i,j)-rhoA(i,j))*ad_gzetaSA(i,j)
              ad_rhoS(i,j)=ad_rhoS(i,j)+adfac
              ad_rhoA(i,j)=ad_rhoA(i,j)-adfac
              ad_gzetaSA(i,j)=0.0_r8
!>            tl_gzeta2(i,j)=tl_gzeta(i,j)*zwrk(i,j)+                   &
!>   &                       gzeta(i,j)*tl_zwrk(i,j)
!>
              ad_gzeta(i,j)=ad_gzeta(i,j)+zwrk(i,j)*ad_gzeta2(i,j)
              ad_zwrk(i,j)=ad_zwrk(i,j)+gzeta(i,j)*ad_gzeta2(i,j)
              ad_gzeta2(i,j)=0.0_r8
!>            tl_gzeta(i,j)=(fac+rhoS(i,j))*tl_zwrk(i,j)+               &
!>   &                      tl_rhoS(i,j)*zwrk(i,j)
!>
              ad_rhoS(i,j)=ad_rhoS(i,j)+zwrk(i,j)*ad_gzeta(i,j)
              ad_zwrk(i,j)=ad_zwrk(i,j)+(fac+rhoS(i,j))*ad_gzeta(i,j)
              ad_gzeta(i,j)=0.0_r8
# else
!>            tl_gzeta2(i,j)=2.0_r8*tl_zwrk(i,j)*zwrk(i,j)
!>            tl_gzeta(i,j)=tl_zwrk(i,j)
!>
              ad_zwrk(i,j)=ad_zwrk(i,j)+                                &
     &                     2.0_r8*zwrk(i,j)*ad_gzeta2(i,j)+             &
     &                     ad_gzeta(i,j)
              ad_gzeta2(i,j)=0.0_r8
              ad_gzeta(i,j)=0.0_r8
# endif
!>            tl_zwrk(i,j)=0.5_r8*(tl_zeta(i,j,kstp)+tl_zeta_new(i,j))
!>
              adfac=0.5_r8*ad_zwrk(i,j)
              ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+adfac
              ad_zeta_new(i,j)=ad_zeta_new(i,j)+adfac
              ad_zwrk(i,j)=0.0_r8
!>            tl_Dnew(i,j)=tl_zeta_new(i,j)+tl_h(i,j)
!>
              ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_Dnew(i,j)
              ad_h(i,j)=ad_h(i,j)+ad_Dnew(i,j)
              ad_Dnew(i,j)=0.0_r8
# ifdef MASKING
!>            tl_zeta_new(i,j)=tl_zeta_new(i,j)*rmask(i,j)
!>
              ad_zeta_new(i,j)=ad_zeta_new(i,j)*rmask(i,j)
# endif
!>            tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+                       &
!>   &                         pm(i,j)*pn(i,j)*cff1*tl_rhs_zeta(i,j)
!>
              ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_zeta_new(i,j)
              ad_rhs_zeta(i,j)=ad_rhs_zeta(i,j)+                        &
     &                         pm(i,j)*pn(i,j)*cff1*ad_zeta_new(i,j)
              ad_zeta_new(i,j)=0.0_r8
!>            tl_rhs_zeta(i,j)=(tl_DUon(i,j)-tl_DUon(i+1,j))+           &
!>   &                         (tl_DVom(i,j)-tl_DVom(i,j+1))
!>
              ad_DUon(i  ,j  )=ad_DUon(i  ,j  )+ad_rhs_zeta(i,j)
              ad_DUon(i+1,j  )=ad_DUon(i+1,j  )-ad_rhs_zeta(i,j)
              ad_DVom(i  ,j  )=ad_DVom(i  ,j  )+ad_rhs_zeta(i,j)
              ad_DVom(i  ,j+1)=ad_DVom(i  ,j+1)-ad_rhs_zeta(i,j)
              ad_rhs_zeta(i,j)=0.0_r8
            END DO
          END DO
        ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
          cff1=2.0_r8*dtfast(ng)
          cff4=4.0_r8/25.0_r8
          cff5=1.0_r8-2.0_r8*cff4
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
# if defined VAR_RHO_2D && defined SOLVE3D
!>            tl_gzetaSA(i,j)=tl_zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))+       &
!>   &                        zwrk(i,j)*(tl_rhoS(i,j)-tl_rhoA(i,j))
!>
              adfac=zwrk(i,j)*ad_gzetaSA(i,j)
              ad_zwrk(i,j)=ad_zwrk(i,j)+                                &
     &                     (rhoS(i,j)-rhoA(i,j))*ad_gzetaSA(i,j)
              ad_rhoS(i,j)=ad_rhoS(i,j)+adfac
              ad_rhoA(i,j)=ad_rhoA(i,j)-adfac
              ad_gzetaSA(i,j)=0.0_r8
!>            tl_gzeta2(i,j)=tl_gzeta(i,j)*zwrk(i,j)+                   &
!>   &                       gzeta(i,j)*tl_zwrk(i,j)
!>
              ad_gzeta(i,j)=ad_gzeta(i,j)+zwrk(i,j)*ad_gzeta2(i,j)
              ad_zwrk(i,j)=ad_zwrk(i,j)+gzeta(i,j)*ad_gzeta2(i,j)
              ad_gzeta2(i,j)=0.0_r8
!>            tl_gzeta(i,j)=(fac+rhoS(i,j))*tl_zwrk(i,j)+               &
!>   &                      tl_rhoS(i,j)*zwrk(i,j)
!>
              ad_zwrk(i,j)=ad_zwrk(i,j)+(fac+rhoS(i,j))*ad_gzeta(i,j)
              ad_rhoS(i,j)=ad_rhoS(i,j)+zwrk(i,j)*ad_gzeta(i,j)
              ad_gzeta(i,j)=0.0_r8
# else
!>            tl_gzeta2(i,j)=2.0_r8*tl_zwrk(i,j)*zwrk(i,j)
!>            tl_gzeta(i,j)=tl_zwrk(i,j)
!>
              ad_zwrk(i,j)=ad_zwrk(i,j)+                                &
     &                     2.0_r8*zwrk(i,j)*ad_gzeta2(i,j)+             &
     &                     ad_gzeta(i,j)
              ad_gzeta2(i,j)=0.0_r8
              ad_gzeta(i,j)=0.0_r8
# endif
!>            tl_zwrk(i,j)=cff5*tl_zeta(i,j,krhs)+                      &
!>   &                     cff4*(tl_zeta(i,j,kstp)+tl_zeta_new(i,j))
!>
              adfac=cff4*ad_zwrk(i,j)
              ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+cff5*ad_zwrk(i,j)
              ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+adfac
              ad_zeta_new(i,j)=ad_zeta_new(i,j)+adfac
              ad_zwrk(i,j)=0.0_r8
!>            tl_Dnew(i,j)=tl_zeta_new(i,j)+tl_h(i,j)
!>
              ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_Dnew(i,j)
              ad_h(i,j)=ad_h(i,j)+ad_Dnew(i,j)
              ad_Dnew(i,j)=0.0_r8
# ifdef MASKING
!>            tl_zeta_new(i,j)=tl_zeta_new(i,j)*rmask(i,j)
!>
              ad_zeta_new(i,j)=ad_zeta_new(i,j)*rmask(i,j)
# endif
!>            tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+                       &
!>   &                         pm(i,j)*pn(i,j)*cff1*tl_rhs_zeta(i,j)
!>
              ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_zeta_new(i,j)
              ad_rhs_zeta(i,j)=ad_rhs_zeta(i,j)+                        &
     &                         pm(i,j)*pn(i,j)*cff1*ad_zeta_new(i,j)
              ad_zeta_new(i,j)=0.0_r8
!>            tl_rhs_zeta(i,j)=(tl_DUon(i,j)-tl_DUon(i+1,j))+           &
!>   &                         (tl_DVom(i,j)-tl_DVom(i,j+1))
!>
              ad_DUon(i  ,j  )=ad_DUon(i  ,j  )+ad_rhs_zeta(i,j)
              ad_DUon(i+1,j  )=ad_DUon(i+1,j  )-ad_rhs_zeta(i,j)
              ad_DVom(i  ,j  )=ad_DVom(i  ,j  )+ad_rhs_zeta(i,j)
              ad_DVom(i  ,j+1)=ad_DVom(i  ,j+1)-ad_rhs_zeta(i,j)
              ad_rhs_zeta(i,j)=0.0_r8
            END DO
          END DO
        ELSE IF (CORRECTOR_2D_STEP) THEN
          cff1=dtfast(ng)*5.0_r8/12.0_r8
          cff2=dtfast(ng)*8.0_r8/12.0_r8
          cff3=dtfast(ng)*1.0_r8/12.0_r8
          cff4=2.0_r8/5.0_r8
          cff5=1.0_r8-cff4
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
# if defined VAR_RHO_2D && defined SOLVE3D
!>            tl_gzetaSA(i,j)=tl_zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))+       &
!>   &                        zwrk(i,j)*(tl_rhoS(i,j)-tl_rhoA(i,j))
!>
              adfac=zwrk(i,j)*ad_gzetaSA(i,j)
              ad_zwrk(i,j)=ad_zwrk(i,j)+                                &
     &                     (rhoS(i,j)-rhoA(i,j))*ad_gzetaSA(i,j)
              ad_rhoS(i,j)=ad_rhoS(i,j)+adfac
              ad_rhoA(i,j)=ad_rhoA(i,j)-adfac
              ad_gzetaSA(i,j)=0.0_r8
!>            tl_gzeta2(i,j)=tl_gzeta(i,j)*zwrk(i,j)+                   &
!>   &                       gzeta(i,j)*tl_zwrk(i,j)
!>
              ad_zwrk(i,j)=ad_zwrk(i,j)+gzeta(i,j)*ad_gzeta2(i,j)
              ad_gzeta(i,j)=ad_gzeta(i,j)+zwrk(i,j)*ad_gzeta2(i,j)
              ad_gzeta2(i,j)=0.0_r8
!>            tl_gzeta(i,j)=(fac+rhoS(i,j))*tl_zwrk(i,j)+               &
!>   &                      tl_rhoS(i,j)*zwrk(i,j)
!>
              ad_zwrk(i,j)=ad_zwrk(i,j)+(fac+rhoS(i,j))*ad_gzeta(i,j)
              ad_rhoS(i,j)=ad_rhoS(i,j)+zwrk(i,j)*ad_gzeta(i,j)
              ad_gzeta(i,j)=0.0_r8
# else
!>            tl_gzeta(i,j)=tl_zwrk(i,j)
!>            tl_gzeta2(i,j)=2.0_r8*tl_zwrk(i,j)*zwrk(i,j)
!>
              ad_zwrk(i,j)=ad_zwrk(i,j)+                                &
     &                     2.0_r8*zwrk(i,j)*ad_gzeta2(i,j)+             &
     &                     ad_gzeta(i,j)
              ad_gzeta2(i,j)=0.0_r8
              ad_gzeta(i,j)=0.0_r8
# endif
!>            tl_zwrk(i,j)=cff5*tl_zeta_new(i,j)+cff4*tl_zeta(i,j,krhs)
!>
              ad_zeta_new(i,j)=ad_zeta_new(i,j)+cff5*ad_zwrk(i,j)
              ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+cff4*ad_zwrk(i,j)
              ad_zwrk(i,j)=0.0_r8
!>            tl_Dnew(i,j)=tl_zeta_new(i,j)+tl_h(i,j)
!>
              ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_Dnew(i,j)
              ad_h(i,j)=ad_h(i,j)+ad_Dnew(i,j)
              ad_Dnew(i,j)=0.0_r8
# ifdef MASKING
!>            tl_zeta_new(i,j)=tl_zeta_new(i,j)*rmask(i,j)
!>
              ad_zeta_new(i,j)=ad_zeta_new(i,j)*rmask(i,j)
# endif
!>            tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+                       &
!>   &                         pm(i,j)*pn(i,j)*(tl_cff+                 &
!>   &                                          cff2*tl_rzeta(i,j,kstp)-&
!>   &                                          cff3*tl_rzeta(i,j,ptsk))
!>
              adfac=pm(i,j)*pn(i,j)*ad_zeta_new(i,j)
              ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_zeta_new(i,j)
              ad_cff=ad_cff+adfac
              ad_rzeta(i,j,kstp)=ad_rzeta(i,j,kstp)+adfac*cff2
              ad_rzeta(i,j,ptsk)=-adfac*cff3
              ad_zeta_new(i,j)=0.0_r8
!>            tl_cff=cff1*((tl_DUon(i,j)-tl_DUon(i+1,j))+               &
!>   &                     (tl_DVom(i,j)-tl_DVom(i,j+1)))
!>
              adfac=cff1*ad_cff
              ad_DUon(i  ,j  )=ad_DUon(i  ,j  )+adfac
              ad_DUon(i+1,j  )=ad_DUon(i+1,j  )-adfac
              ad_DVom(i  ,j  )=ad_DVom(i  ,j  )+adfac
              ad_DVom(i  ,j+1)=ad_DVom(i  ,j+1)-adfac
              ad_cff=0.0_r8
            END DO
          END DO
        END IF
      END IF STEP_LOOP
# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute adjoint time averaged fields over all short timesteps.
!-----------------------------------------------------------------------
!
!  After all fast time steps are completed, recompute S-coordinate
!  surfaces according to the new free surface field. Apply boundary
!  conditions to time averaged fields.
!
      IF ((iif(ng).eq.(nfast(ng)+1)).and.PREDICTOR_2D_STEP(ng)) THEN
!>      CALL tl_set_depth_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          IminS, ImaxS, JminS, JmaxS,             &
!>   &                          nstp, nnew,                             &
!>   &                          h, tl_h,                                &
#  ifdef ICESHELF
!>   &                          zice,                                   &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>   &                          tl_bed_thick,                           &
#  endif
!>   &                          Zt_avg1, tl_Zt_avg1,                    &
!>   &                          tl_Hz, tl_z_r, tl_z_w)
!>
        CALL ad_set_depth_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp, nnew,                             &
     &                          h, ad_h,                                &
#  ifdef ICESHELF
     &                          zice,                                   &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                          ad_bed_thick,                           &
#  endif
     &                          Zt_avg1, ad_Zt_avg1,                    &
     &                          ad_Hz, ad_z_r, ad_z_w)
#  ifdef DISTRIBUTE
!>      CALL mp_exchange2d (ng, tile, iTLM, 3,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints, EWperiodic, NSperiodic,       &
!>   &                      tl_Zt_avg1, tl_DU_avg1, tl_DV_avg1)
!>
        CALL ad_mp_exchange2d (ng, tile, iADM, 3,                       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints, EWperiodic, NSperiodic,    &
     &                         ad_Zt_avg1, ad_DU_avg1, ad_DV_avg1)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_DV_avg1)
!>
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_DV_avg1)
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_DU_avg1)
!>
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_DU_avg1)
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_Zt_avg1)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_Zt_avg1)
#  endif
      END IF
!
!  Compute time-averaged fields.
!
      IF (PREDICTOR_2D_STEP(ng)) THEN
        IF (FIRST_2D_STEP) THEN
!
!  Reset arrays for 2D fields averaged within the short time-steps.
!
          cff2=(-1.0_r8/12.0_r8)*weight(2,iif(ng)+1,ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
!>            tl_Zt_avg1(i,j)=0.0_r8
!>
              ad_Zt_avg1(i,j)=0.0_r8
            END DO
            DO i=Istr,IendR
!>            tl_DU_avg2(i,j)=cff2*tl_DUon(i,j)
!>
              ad_DUon(i,j)=ad_DUon(i,j)+cff2*ad_DU_avg2(i,j)
              ad_DU_avg2(i,j)=0.0_r8
!>            tl_DU_avg1(i,j)=0.0_r8
!>
              ad_DU_avg1(i,j)=0.0_r8
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
!>            tl_DV_avg2(i,j)=cff2*tl_DVom(i,j)
!>
              ad_DVom(i,j)=ad_DVom(i,j)+cff2*ad_DV_avg2(i,j)
              ad_DV_avg2(i,j)=0.0_r8
!>            tl_DV_avg1(i,j)=0.0_r8
!>
              ad_DV_avg1(i,j)=0.0_r8
            END DO
          END DO
        ELSE
!
!  Accumulate field averages of previous time-step after they are
!  computed in the previous corrector step, updated their boundaries,
!  and synchronized.
!
          cff1=weight(1,iif(ng)-1,ng)
          cff2=(8.0_r8/12.0_r8)*weight(2,iif(ng)  ,ng)-                 &
     &         (1.0_r8/12.0_r8)*weight(2,iif(ng)+1,ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
!>            tl_Zt_avg1(i,j)=tl_Zt_avg1(i,j)+cff1*tl_zeta(i,j,krhs)
!>
              ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+cff1*ad_Zt_avg1(i,j)
            END DO
            DO i=Istr,IendR
!>            tl_DU_avg2(i,j)=tl_DU_avg2(i,j)+cff2*tl_DUon(i,j)
!>
              ad_DUon(i,j)=ad_DUon(i,j)+                                &
     &                     cff2*ad_DU_avg2(i,j)
#  ifdef NEARSHORE_MELLOR
!>            tl_DU_avg1(i,j)=tl_DU_avg1(i,j)-cff1*tl_DUSon(i,j)
!>
              ad_DUSon(i,j)=ad_DUSon(i,j)-                              &
     &                      cff1*ad_DU_avg1(i,j)
#  endif
!>            tl_DU_avg1(i,j)=tl_DU_avg1(i,j)+cff1*tl_DUon(i,j)
!>
              ad_DUon(i,j)=ad_DUon(i,j)+                                &
     &                     cff1*ad_DU_avg1(i,j)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
!>            tl_DV_avg2(i,j)=tl_DV_avg2(i,j)+cff2*tl_DVom(i,j)
!>
              ad_DVom(i,j)=ad_DVom(i,j)+                                &
     &                     cff2*ad_DV_avg2(i,j)
#  ifdef NEARSHORE_MELLOR
!>            tl_DV_avg1(i,j)=tl_DV_avg1(i,j)-cff1*tl_DVSom(i,j)
!>
              ad_DVSom(i,j)=ad_DVSom(i,j)-                              &
     &                      cff1*ad_DV_avg1(i,j)
#  endif
!>            tl_DV_avg1(i,j)=tl_DV_avg1(i,j)+cff1*tl_DVom(i,j)
!>
              ad_DVom(i,j)=ad_DVom(i,j)+                                &
     &                     cff1*ad_DV_avg1(i,j)
            END DO
          END DO
        END IF
      ELSE
        IF (FIRST_2D_STEP) THEN
          cff2=weight(2,iif(ng),ng)
        ELSE
          cff2=(5.0_r8/12.0_r8)*weight(2,iif(ng),ng)
        END IF
        DO j=JstrR,JendR
          DO i=Istr,IendR
!>          tl_DV_avg2(i,j)=tl_DV_avg2(i,j)+cff2*tl_DVom(i,j)
!>
            ad_DVom(i,j)=ad_DVom(i,j)+cff2*ad_DV_avg2(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
!>          tl_DU_avg2(i,j)=tl_DU_avg2(i,j)+cff2*tl_DUon(i,j)
!>
            ad_DUon(i,j)=ad_DUon(i,j)+cff2*ad_DU_avg2(i,j)
          END DO
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute total depth (m) and vertically integrated mass fluxes.
!-----------------------------------------------------------------------

# ifdef OBC_VOLCONS
!
!  Set vertically integrated mass fluxes DUon and DVom along the open
!  boundaries in such a way that the integral volume is conserved.
!
!>    CALL tl_set_DUV_bc_tile (ng, tile,                                &
!>   &                         LBi, UBi, LBj, UBj,                      &
!>   &                         IminS, ImaxS, JminS, JmaxS,              &
!>   &                         krhs,                                    &
#  ifdef MASKING
!>   &                         umask, vmask,                            &
#  endif
!>   &                         om_v, on_u, ubar, vbar,                  &
!>   &                         tl_ubar, tl_vbar,                        &
!>   &                         Drhs, DUon, DVom,                        &
!>   &                         tl_Drhs, tl_DUon, tl_DVom)
!>
      CALL ad_set_DUV_bc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         krhs,                                    &
#  ifdef MASKING
     &                         umask, vmask,                            &
#  endif
     &                         om_v, on_u, ubar, vbar,                  &
     &                         ad_ubar, ad_vbar,                        &
     &                         Drhs, DUon, DVom,                        &
     &                         ad_Drhs, ad_DUon, ad_DVom)
# endif
# ifdef DISTRIBUTE
!
!  Do a special exchange to avoid having three ghost points for
!  high order numerical stencil. Notice that a private array is
!  passed to the "exchange" routine.  It will also apply periodic
!  boundary conditions if no partitions in I- or J-directions.
!
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_DUon, tl_DVom)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_DUon, ad_DVom)
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_u2d_tile (ng, tile,                                 &
!>   &                        IminS, ImaxS, JminS, JmaxS,               &
!>   &                        tl_DUon)
!>
      CALL ad_exchange_u2d_tile (ng, tile,                              &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ad_DUon)
!>    CALL exchange_v2d_tile (ng, tile,                                 &
!>   &                        IminS, ImaxS, JminS, JmaxS,               &
!>   &                        tl_DVom)
!>
      CALL ad_exchange_v2d_tile (ng, tile,                              &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ad_DVom)
#  endif
# endif
!
!  Compute adjoint adjoint vertically integrated mass fluxes.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend+1
#  else
#   define I_RANGE IstrU,MIN(Iend+1,Lm(ng))
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV,Jend+1
#  else
#   define J_RANGE JstrV,MIN(Jend+1,Mm(ng))
#  endif
# else
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU-1,Iend+1
#  else
#   define I_RANGE MAX(2,IstrU-1),MIN(Iend+1,Lm(ng))
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE JstrV-1,Jend+1
#  else
#   define J_RANGE MAX(2,JstrV-1),MIN(Jend+1,Mm(ng))
#  endif
# endif
      DO j=-1+J_RANGE+1
        DO i=-2+I_RANGE+1
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
# ifdef NEARSHORE_MELLOR
!>        tl_DVom(i,j)=tl_DVom(i,j)+tl_DVSom(i,j)
!>
          ad_DVSom(i,j)=ad_DVSom(i,j)+ad_DVom(i,j)
!>        tl_DVSom(i,j)=tl_vbar_stokes(i,j)*cff1+                       &
!>   &                  vbar_stokes(i,j)*tl_cff1
!>
          ad_cff1=ad_cff1+vbar_stokes(i,j)*ad_DVSom(i,j)
          ad_vbar_stokes(i,j)=ad_vbar_stokes(i,j)+cff1*ad_DVSom(i,j)
          ad_DVSom(i,j)=0.0_r8
# endif
!>        tl_DVom(i,j)=tl_vbar(i,j,krhs)*cff1+                          &
!>   &                 vbar(i,j,krhs)*tl_cff1
!>
          ad_cff1=ad_cff1+vbar(i,j,krhs)*ad_DVom(i,j)
          ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)+cff1*ad_DVom(i,j)
          ad_DVom(i,j)=0.0_r8
!>        tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i,j-1))
!>
          adfac=cff*ad_cff1
          ad_Drhs(i,j-1)=ad_Drhs(i,j-1)+adfac
          ad_Drhs(i,j  )=ad_Drhs(i,j  )+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
      DO j=-2+J_RANGE+1
        DO i=-1+I_RANGE+1
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
# ifdef NEARSHORE_MELLOR
!>        tl_DUon(i,j)=tl_DUon(i,j)+tl_DUSon(i,j)
!>
          ad_DUSon(i,j)=ad_DUSon(i,j)+ad_DUon(i,j)
!>        tl_DUSon(i,j)=tl_ubar_stokes(i,j)*cff1+                       &
!>   &                  ubar_stokes(i,j)*tl_cff1
!>
          ad_cff1=ad_cff1+ubar_stokes(i,j)*ad_DUSon(i,j)
          ad_ubar_stokes(i,j)=ad_ubar_stokes(i,j)+cff1*ad_DUSon(i,j)
          ad_DUSon(i,j)=0.0_r8
# endif
!>        tl_DUon(i,j)=tl_ubar(i,j,krhs)*cff1+                          &
!>   &                 ubar(i,j,krhs)*tl_cff1
!>
          ad_cff1=ad_cff1+ubar(i,j,krhs)*ad_DUon(i,j)
          ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)+cff1*ad_DUon(i,j)
          ad_DUon(i,j)=0.0_r8
!>        tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i-1,j))
!>
          adfac=cff*ad_cff1
          ad_Drhs(i-1,j)=ad_Drhs(i-1,j)+adfac
          ad_Drhs(i  ,j)=ad_Drhs(i  ,j)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
!
!  Compute adjoint total depth.
!
      DO j=-2+J_RANGE+1
        DO i=-2+I_RANGE+1
!>        tl_Drhs(i,j)=tl_zeta(i,j,krhs)+tl_h(i,j)
!>
          ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+ad_Drhs(i,j)
          ad_h(i,j)=ad_h(i,j)+ad_Drhs(i,j)
          ad_Drhs(i,j)=0.0_r8
        END DO
      END DO
# undef I_RANGE
# undef J_RANGE
      RETURN
      END SUBROUTINE ad_step2d_tile
#else
      SUBROUTINE ad_step2d
      END SUBROUTINE ad_step2d
#endif
