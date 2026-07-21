#define DEBUG

      MODULE ad_step2d_mod
!
!git $Id$
!=======================================================================
!                                                                      !
!  Solves adjoint model shallow-water primitive equations (barotropic  !
!  mode) using the Generalized Forward-Backward 3rd-order Adams-       !
!  Bashforth / 4th-order Adams-Moulton (FB AB3-AM4) time stepping      !
!  algorithm (Shchepetkin and McWilliams, 2005). In 3D applications,   !
!  it perform fast-time averaging to interact with 3D momentum         !
!  equations (baroclinic mode).                                        !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic   !
!     modeling system (ROMS): a split-explicit, free-surface,          !
!     topography-following-coordinate oceanic model, Ocean Modelling,  !
!     9, 347-404, doi:10.1016/j.ocemod.2004.08.002.                    !
!                                                                      !
!  Shchepetkin, A.F., and J.C. McWilliams, 2009: Computational kernel  !
!     algorithms for fine-scale, multiprocess, longtime oceanic        !
!     simulations, pp 121-183. In 'Handbook of Numerical Analysis:     !
!     Computational Methods for the Atmosphere and Oceans', R.M. Teman !
!     and J.J. Tribbia, eds, Elsevier Science.                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef SOLVE3D
      USE mod_coupling
#endif
#ifdef DIAGNOSTICS_UV
!!    USE mod_diags
#endif
      USE mod_forces
      USE mod_grid
#if defined UV_VIS2 || defined UV_VIS4
      USE mod_mixing
#endif
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET && \
    defined SOLVE3D
      USE mod_sedbed
#endif
      USE mod_sources
      USE mod_stepping
!
      USE ad_exchange_2d_mod
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod,    ONLY : ad_mp_exchange2d
      USE mp_exchange_mod,    ONLY : mp_exchange2d
#endif
      USE obc_volcons_mod,    ONLY : obc_flux_tile,                     &
     &                               set_DUV_bc_tile
      USE ad_obc_volcons_mod, ONLY : ad_obc_flux_tile,                  &
     &                               ad_set_DUV_bc_tile
#ifdef SOLVE3D
      USE ad_set_depth_mod,   ONLY : ad_set_depth
#endif
      USE ad_u2dbc_mod,       ONLY : ad_u2dbc_tile
      USE ad_v2dbc_mod,       ONLY : ad_v2dbc_tile
      USE ad_zetabc_mod,      ONLY : ad_zetabc_local
#ifdef WET_DRY_NOT_YET
      USE wetdry_mod,         ONLY : wetdry_tile
#endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ad_step2d
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_step2d (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iADM, 9, __LINE__, MyFile)
#endif
      CALL ad_step2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj, N(ng),                   &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs(ng), kstp(ng), knew(ng),                &
#ifdef SOLVE3D
     &                     nstp(ng), nnew(ng),                          &
#endif
#ifdef MASKING
     &                     GRID(ng) % pmask,     GRID(ng) % rmask,      &
     &                     GRID(ng) % umask,     GRID(ng) % vmask,      &
#endif
#ifdef WET_DRY_NOT_YET
     &                     GRID(ng) % pmask_wet, GRID(ng) % pmask_full, &
     &                     GRID(ng) % rmask_wet, GRID(ng) % rmask_full, &
     &                     GRID(ng) % umask_wet, GRID(ng) % umask_full, &
     &                     GRID(ng) % vmask_wet, GRID(ng) % vmask_full, &
# ifdef SOLVE3D
     &                     GRID(ng) % rmask_wet_avg,                    &
# endif
#endif
#if (defined UV_COR && !defined SOLVE3D) || defined STEP2D_CORIOLIS
     &                     GRID(ng) % fomn,                             &
#endif
     &                     GRID(ng) % h,         GRID(ng) % ad_h,       &
     &                     GRID(ng) % om_u,      GRID(ng) % om_v,       &
     &                     GRID(ng) % on_u,      GRID(ng) % on_v,       &
     &                     GRID(ng) % pm,        GRID(ng) % pn,         &
#if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
     &                     GRID(ng) % dndx,      GRID(ng) % dmde,       &
#endif
     &                     GRID(ng) % rdrag,                            &
#if defined UV_QDRAG && !defined SOLVE3D
     &                     GRID(ng) % rdrag2,                           &
#endif
#if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
     &                     GRID(ng) % pmon_r,    GRID(ng) % pnom_r,     &
     &                     GRID(ng) % pmon_p,    GRID(ng) % pnom_p,     &
     &                     GRID(ng) % om_r,      GRID(ng) % on_r,       &
     &                     GRID(ng) % om_p,      GRID(ng) % on_p,       &
# ifdef UV_VIS2
     &                     MIXING(ng) % visc2_p,                        &
     &                     MIXING(ng) % visc2_r,                        &
# endif
# ifdef UV_VIS4
     &                     MIXING(ng) % visc4_p,                        &
     &                     MIXING(ng) % visc4_r,                        &
# endif
#endif
#if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                     SEDBED(ng) % ad_bed_thick,                   &
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
     &                     OCEAN(ng) % eq_tide,                         &
     &                     OCEAN(ng) % ad_eq_tide,                      &
#endif
#ifndef SOLVE3D
     &                     FORCES(ng) % sustr,   FORCES(ng) % ad_sustr, &
     &                     FORCES(ng) % svstr,   FORCES(ng) % ad_svstr, &
# ifdef ATM_PRESS
     &                     FORCES(ng) % Pair,                           &
# endif
#else
# ifdef VAR_RHO_2D
     &                     COUPLING(ng) % rhoA,                         &
     &                     COUPLING(ng) % ad_rhoA,                      &
     &                     COUPLING(ng) % rhoS,                         &
     &                     COUPLING(ng) % ad_rhoS,                      &
# endif
     &                     COUPLING(ng) % ad_DU_avg1,                   &
     &                     COUPLING(ng) % ad_DU_avg2,                   &
     &                     COUPLING(ng) % ad_DV_avg1,                   &
     &                     COUPLING(ng) % ad_DV_avg2,                   &
     &                     COUPLING(ng) % ad_Zt_avg1,                   &
     &                     COUPLING(ng) % rufrc,                        &
     &                     COUPLING(ng) % ad_rufrc,                     &
     &                     COUPLING(ng) % rvfrc,                        &
     &                     COUPLING(ng) % ad_rvfrc,                     &
     &                     COUPLING(ng) % ad_rufrc_bak,                 &
     &                     COUPLING(ng) % ad_rvfrc_bak,                 &
#endif
#if defined NESTING && !defined SOLVE3D
     &                     OCEAN(ng) % ad_DU_flux,                      &
     &                     OCEAN(ng) % ad_DV_flux,                      &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                     DIAGS(ng) % DiaU2wrk, DIAGS(ng) % DiaV2wrk,  &
!!   &                     DIAGS(ng) % DiaRUbar, DIAGS(ng) % DiaRVbar,  &
# ifdef SOLVE3D
!!   &                     DIAGS(ng) % DiaU2int, DIAGS(ng) % DiaV2int,  &
!!   &                     DIAGS(ng) % DiaRUfrc, DIAGS(ng) % DiaRVfrc,  &
# endif
#endif
     &                     OCEAN(ng) % ad_ubar_sol,                     &
     &                     OCEAN(ng) % ad_vbar_sol,                     &
     &                     OCEAN(ng) % ad_zeta_sol,                     &
     &                     OCEAN(ng) % ubar,     OCEAN(ng) % ad_ubar,   &
     &                     OCEAN(ng) % vbar,     OCEAN(ng) % ad_vbar,   &
     &                     OCEAN(ng) % zeta,     OCEAN(ng) % ad_zeta)
#ifdef PROFILE
      CALL wclock_off (ng, iADM, 9, __LINE__, MyFile)
#endif
!
      RETURN
      END SUBROUTINE ad_step2d
!
!***********************************************************************
      SUBROUTINE ad_step2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, UBk,               &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           krhs, kstp, knew,                      &
#ifdef SOLVE3D
     &                           nstp, nnew,                            &
#endif
#ifdef MASKING
     &                           pmask, rmask, umask, vmask,            &
#endif
#ifdef WET_DRY_NOT_YET
     &                           pmask_wet, pmask_full,                 &
     &                           rmask_wet, rmask_full,                 &
     &                           umask_wet, umask_full,                 &
     &                           vmask_wet, vmask_full,                 &
# ifdef SOLVE3D
     &                           rmask_wet_avg,                         &
# endif
#endif
#if (defined UV_COR && !defined SOLVE3D) || defined STEP2D_CORIOLIS
     &                           fomn,                                  &
#endif
     &                           h, ad_h,                               &
     &                           om_u, om_v, on_u, on_v, pm, pn,        &
#if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
     &                           dndx, dmde,                            &
#endif
     &                           rdrag,                                 &
#if defined UV_QDRAG && !defined SOLVE3D
     &                           rdrag2,                                &
#endif
#if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
     &                           pmon_r, pnom_r, pmon_p, pnom_p,        &
     &                           om_r, on_r, om_p, on_p,                &
# ifdef UV_VIS2
     &                           visc2_p, visc2_r,                      &
# endif
# ifdef UV_VIS4
     &                           visc4_p, visc4_r,                      &
# endif
#endif
#if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                           ad_bed_thick,                          &
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
     &                           eq_tide, ad_eq_tide,                   &
#endif
#ifndef SOLVE3D
     &                           sustr, ad_sustr,                       &
     &                           svstr, ad_svstr,                       &
# ifdef ATM_PRESS
     &                           Pair,                                  &
# endif
#else
# ifdef VAR_RHO_2D
     &                           rhoA, ad_rhoA,                         &
     &                           rhoS, ad_rhoS,                         &
# endif
     &                           ad_DU_avg1, ad_DU_avg2,                &
     &                           ad_DV_avg1, ad_DV_avg2,                &
     &                           ad_Zt_avg1,                            &
     &                           rufrc, ad_rufrc,                       &
     &                           rvfrc, ad_rvfrc,                       &
     &                           ad_rufrc_bak, ad_rvfrc_bak,            &
#endif
#if defined NESTING && !defined SOLVE3D
     &                           ad_DU_flux, ad_DV_flux,                &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                           DiaU2wrk, DiaV2wrk,                    &
!!   &                           DiaRUbar, DiaRVbar,                    &
# ifdef SOLVE3D
!!   &                           DiaU2int, DiaV2int,                    &
!!   &                           DiaRUfrc, DiaRVfrc,                    &
# endif
#endif
     &                           ad_ubar_sol,                           &
     &                           ad_vbar_sol,                           &
     &                           ad_zeta_sol,                           &
     &                           ubar, ad_ubar,                         &
     &                           vbar, ad_vbar,                         &
     &                           zeta, ad_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in    ) :: ng, tile
      integer, intent(in    ) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in    ) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in    ) :: krhs, kstp, knew
#ifdef SOLVE3D
      integer, intent(in    ) :: nstp, nnew
#endif
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in   ) :: pmask(LBi:,LBj:)
      real(r8), intent(in   ) :: rmask(LBi:,LBj:)
      real(r8), intent(in   ) :: umask(LBi:,LBj:)
      real(r8), intent(in   ) :: vmask(LBi:,LBj:)
# endif
# if (defined UV_COR && !defined SOLVE3D) || defined STEP2D_CORIOLIS
      real(r8), intent(in   ) :: fomn(LBi:,LBj:)
# endif
      real(r8), intent(in   ) :: h(LBi:,LBj:)
      real(r8), intent(in   ) :: om_u(LBi:,LBj:)
      real(r8), intent(in   ) :: om_v(LBi:,LBj:)
      real(r8), intent(in   ) :: on_u(LBi:,LBj:)
      real(r8), intent(in   ) :: on_v(LBi:,LBj:)
      real(r8), intent(in   ) :: pm(LBi:,LBj:)
      real(r8), intent(in   ) :: pn(LBi:,LBj:)
# if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
      real(r8), intent(in   ) :: dndx(LBi:,LBj:)
      real(r8), intent(in   ) :: dmde(LBi:,LBj:)
# endif
      real(r8), intent(in   ) :: rdrag(LBi:,LBj:)
# if defined UV_QDRAG && !defined SOLVE3D
      real(r8), intent(in   ) :: rdrag2(LBi:,LBj:)
# endif
      real(r8), intent(in   ) :: rufrc(LBi:,LBj:)
      real(r8), intent(in   ) :: rvfrc(LBi:,LBj:)
# if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
      real(r8), intent(in   ) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in   ) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in   ) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in   ) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in   ) :: om_r(LBi:,LBj:)
      real(r8), intent(in   ) :: on_r(LBi:,LBj:)
      real(r8), intent(in   ) :: om_p(LBi:,LBj:)
      real(r8), intent(in   ) :: on_p(LBi:,LBj:)
#  ifdef UV_VIS2
      real(r8), intent(in   ) :: visc2_p(LBi:,LBj:)
      real(r8), intent(in   ) :: visc2_r(LBi:,LBj:)
#  endif
#  ifdef UV_VIS4
      real(r8), intent(in   ) :: visc4_p(LBi:,LBj:)
      real(r8), intent(in   ) :: visc4_r(LBi:,LBj:)
#  endif
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout) :: ad_bed_thick(LBi:,LBj:,:)
# endif
# if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
      real(r8), intent(in   ) :: eq_tide(LBi:,LBj:)
      real(r8), intent(inout) :: ad_eq_tide(LBi:,LBj:)
# endif
      real(r8), intent(in   ) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in   ) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in   ) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_h(LBi:,LBj:)
# ifndef SOLVE3D
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
#  ifdef ATM_PRESS
      real(r8), intent(in   ) :: Pair(LBi:,LBj:)
#  endif
# else
#  ifdef VAR_RHO_2D
      real(r8), intent(in   ) :: rhoA(LBi:,LBj:)
      real(r8), intent(in   ) :: rhoS(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rhoA(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rhoS(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_DU_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DU_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DV_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rufrc_bak(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_rvfrc_bak(LBi:,LBj:,:)
# endif
# ifdef WET_DRY_NOT_YET
      real(r8), intent(inout) :: pmask_full(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_full(LBi:,LBj:)
      real(r8), intent(inout) :: umask_full(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_full(LBi:,LBj:)

      real(r8), intent(inout) :: pmask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: umask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_wet(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: rmask_wet_avg(LBi:,LBj:)
#  endif
# endif
# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaU2wrk(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaV2wrk(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaRUbar(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRVbar(LBi:,LBj:,:,:)
#  ifdef SOLVE3D
!!    real(r8), intent(inout) :: DiaU2int(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaV2int(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
#  endif
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# if defined NESTING && !defined SOLVE3D
      real(r8), intent(inout) :: ad_DU_flux(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DV_flux(LBi:,LBj:)
# endif
      real(r8), intent(out  ) :: ad_ubar_sol(LBi:,LBj:)
      real(r8), intent(out  ) :: ad_vbar_sol(LBi:,LBj:)
      real(r8), intent(out  ) :: ad_zeta_sol(LBi:,LBj:)

#else

# ifdef MASKING
      real(r8), intent(in   ) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: vmask(LBi:UBi,LBj:UBj)
# endif
# if (defined UV_COR && !defined SOLVE3D) || defined STEP2D_CORIOLIS
      real(r8), intent(in   ) :: fomn(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in   ) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pn(LBi:UBi,LBj:UBj)
# if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
      real(r8), intent(in   ) :: dndx(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: dmde(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in   ) :: rdrag(LBi:UBi,LBj:UBj)
# if defined UV_QDRAG && !defined SOLVE3D
      real(r8), intent(in   ) :: rdrag2(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in   ) :: rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: rvfrc(LBi:UBi,LBj:UBj)
# if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
      real(r8), intent(in   ) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: on_p(LBi:UBi,LBj:UBj)
#  ifdef UV_VIS2
      real(r8), intent(in   ) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: visc2_r(LBi:UBi,LBj:UBj)
#  endif
#  ifdef UV_VIS4
      real(r8), intent(in   ) :: visc4_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: visc4_r(LBi:UBi,LBj:UBj)
#  endif
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in   ) :: ad_bed_thick(LBi:UBi,LBj:UBj,3)
# endif
# if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
      real(r8), intent(in   ) :: eq_tide(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_eq_tide(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in   ) :: ubar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(in   ) :: vbar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(in   ) :: zeta(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: ad_h(LBi:UBi,LBj:UBj)
# ifndef SOLVE3D
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
#  ifdef ATM_PRESS
      real(r8), intent(in   ) :: Pair(LBi:UBi,LBj:UBj)
#  endif
# else
#  ifdef VAR_RHO_2D
      real(r8), intent(in   ) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: rhoS(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rhoS(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ad_DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rufrc_bak(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ad_rvfrc_bak(LBi:UBi,LBj:UBj,2)
# endif
# ifdef WET_DRY_NOT_YET
      real(r8), intent(inout) :: pmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_full(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: pmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_wet(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: rmask_wet_avg(LBi:UBi,LBj:UBj)
#  endif
# endif
# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaU2wrk(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaV2wrk(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaRUbar(LBi:UBi,LBj:UBj,2,NDM2d-1)
!!    real(r8), intent(inout) :: DiaRVbar(LBi:UBi,LBj:UBj,2,NDM2d-1)
#  ifdef SOLVE3D
!!    real(r8), intent(inout) :: DiaU2int(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaV2int(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
#  endif
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,:)
# if defined NESTING && !defined SOLVE3D
      real(r8), intent(inout) :: ad_DU_flux(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DV_flux(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out  ) :: ad_ubar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: ad_vbar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: ad_zeta_sol(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, is, j
      integer :: kbak, kold
#ifdef DIAGNOSTICS_UV
      integer :: idiag
#endif
!
      real(r8) :: bkw0, bkw1, bkw2, bkw_new
      real(r8) :: fwd0, fwd1, fwd2
#ifdef SOLVE3D
      real(r8) :: cfwd0, cfwd1, cfwd2
#endif
      real(r8) :: cff,  cff1, cff2, cff3, cff4
#ifdef WET_DRY_NOT_YET
      real(r8) :: cff5, cff6, cff7
#endif
      real(r8) :: fac, fac1, fac2
      real(r8) :: ad_cff,  ad_cff1, ad_cff2, ad_cff3, ad_cff4
#ifdef WET_DRY_NOT_YET
      real(r8) :: ad_cff5, ad_cff6, ad_cff7
#endif
      real(r8) :: ad_fac, ad_fac1, ad_fac2
      real(r8) :: adfac, adfac1, adfac2, adfac3, adfac4, adfac5
!
      real(r8), parameter :: IniVal = 0.0_r8
!
#if defined UV_C4ADVECTION && !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dgrad
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dnew
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dnew_rd
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Drhs
#if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Drhs_p
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: DUon
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: DVom
#if defined UV_C4ADVECTION && !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rzeta2
#if defined VAR_RHO_2D && defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rzetaSA
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rubar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rvbar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rzeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: urhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vrhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zeta_new
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zwrk
#ifdef WET_DRY_NOT_YET
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wetdry
#endif
#ifdef DIAGNOSTICS_UV
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,NDM2d-1) :: DiaU2rhs
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,NDM2d-1) :: DiaV2rhs
#endif
!
#if defined UV_C4ADVECTION && !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Dgrad
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Dnew
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Dnew_rd
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Drhs
#if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Drhs_p
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_DUon
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_DVom
#if defined STEP2D_CORIOLIS || !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFe
#endif
#if !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFx
#endif
#if defined UV_C4ADVECTION && !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_grad
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_rzeta2
#if defined VAR_RHO_2D && defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_rzetaSA
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_rzeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_rubar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_rvbar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_urhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_vrhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_zwrk
!
      real(r8), allocatable :: ad_zeta_new(:,:)

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=IniVal
      ad_cff1=IniVal
      ad_cff2=IniVal
      ad_cff3=IniVal
      ad_cff4=IniVal
      ad_fac=IniVal
      ad_fac1=IniVal
      ad_fac2=IniVal
!
#if defined UV_C4ADVECTION && !defined SOLVE3D
      ad_Dgrad=IniVal
#endif
      ad_Dnew=IniVal
      ad_Dnew_rd=IniVal
      ad_Drhs=IniVal
#if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
      ad_Drhs_p=IniVal
#endif
      ad_Dstp=IniVal
      ad_DUon=IniVal
      ad_DVom=IniVal
#if defined STEP2D_CORIOLIS || !defined SOLVE3D
      ad_UFx=IniVal
      ad_VFe=IniVal
#endif
#if !defined SOLVE3D
      ad_UFe=IniVal
      ad_VFx=IniVal
#endif
#if defined UV_C4ADVECTION && !defined SOLVE3D
      ad_grad=IniVal
#endif
      ad_rzeta2=IniVal
#if defined VAR_RHO_2D && defined SOLVE3D
      ad_rzetaSA=IniVal
#endif
      ad_rzeta=IniVal
      ad_rubar=IniVal
      ad_rvbar=IniVal
      ad_urhs=IniVal
      ad_vrhs=IniVal
      ad_zwrk=IniVal
!
!-----------------------------------------------------------------------
!  Set coefficients for AB3-AM4 forward-backward algorithm.
!-----------------------------------------------------------------------
!
!  Because the Forward Euler step is used to update "zeta" during the
!  first barotropic step, the pressure-gradient term in the momentum
!  equation must be computed via the Backward step to keep it
!  numerically stable. However, this interferes with the computation
!  of forcing terms "rufrc" and "rvfrc" because the free surface in
!  pressure gradient computation in 3D is exactly at the time
!  corresponding to baroclinic step "nstp" (rather than ahead by one
!  barotropic step after it updated by a normal forward-backward step).
!  To resolve this conflict, the pressure gradient term is computed in
!  two stages during the first barotropic step. It uses zeta(:,:,kstp)
!  at first to ensure exact consistency with 3D model. Then, after
!  vertical integrals of 3D right-hand-side "rufrc" and "rvfrc" are
!  converted into forcing terms, add correction based on the difference
!  zeta_new(:,:)-zeta(:,:,kstp) to "rubar" and "rvbar" to make them
!  consistent with the Backward step for pressure gradient.
!  For pressure gradient terms, search for the label PGF_FB_CORRECTION
!  below.
!
      IF (FIRST_2D_STEP) THEN            ! Meaning of time indices
        kbak=kstp                        !------------------------
        kold=kstp                        ! m-2   m-1   m     m+1
        fwd0=1.0_r8                      ! kold  kbak  kstp  knew
        fwd1=0.0_r8                      ! fwd2  fwd1  fwd0
        fwd2=0.0_r8                      ! bkw2  bkw1  bkw0  bkw_new
#ifdef SOLVE3D
        bkw_new=0.0_r8
        bkw0=1.0_r8
#else
        bkw_new=1.0_r8
        bkw0=0.0_r8
#endif
        bkw1=0.0_r8
        bkw2=0.0_r8
      ELSE IF (FIRST_2D_STEP+1) THEN
        kbak=kstp-1
        IF (kbak.lt.1) kbak=4
        kold=kbak
        fwd0=1.0_r8                      ! Logically AB2-AM3 forward-
        fwd1=0.0_r8                      ! backward scheme with maximum
        fwd2=0.0_r8                      ! stability coefficients while
        bkw_new=1.0833333333333_r8       ! maintaining third-order
        bkw0=-0.1666666666666_r8         ! accuracy, alpha_max=1.73
        bkw1= 0.0833333333333_r8
        bkw2= 0.0_r8
      ELSE
        kbak=kstp-1
        IF (kbak.lt.1) kbak=4
        kold=kbak-1
        IF (kold.lt.1) kold=4
        fwd0=1.781105_r8
        fwd1=-1.06221_r8
        fwd2=0.281105_r8
        bkw_new=0.614_r8
        bkw0=0.285_r8
        bkw1=0.0880_r8
        bkw2=0.013_r8
      END IF

#ifdef DEBUG
!
      IF (Master) THEN
        WRITE (20,10) iic(ng), iif(ng), kold, kbak, kstp, knew
 10     FORMAT (' iic = ',i5.5,' iif = ',i3.3,                          &
     &          ' kold = ',i1,' kbak = ',i1,' kstp = ',i1,' knew = ',i1)
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE total depth (m) arrays and vertically
!  integerated mass fluxes.
!-----------------------------------------------------------------------
!
#if defined DISTRIBUTE && !defined NESTING
# define IR_RANGE IstrUm2-1,Iendp2
# define JR_RANGE JstrVm2-1,Jendp2
# define IU_RANGE IstrUm1-1,Iendp2
# define JU_RANGE Jstrm1-1,Jendp2
# define IV_RANGE Istrm1-1,Iendp2
# define JV_RANGE JstrVm1-1,Jendp2
#else
# define IR_RANGE IstrUm2-1,Iendp2
# define JR_RANGE JstrVm2-1,Jendp2
# define IU_RANGE IstrUm2,Iendp2
# define JU_RANGE JstrVm2-1,Jendp2
# define IV_RANGE IstrUm2-1,Iendp2
# define JV_RANGE JstrVm2,Jendp2
#endif

      DO j=JR_RANGE
        DO i=IR_RANGE
!^        Drhs(i,j)=h(i,j)+fwd0*zeta(i,j,kstp)+                         &
!^   &                     fwd1*zeta(i,j,kbak)+                         &
!^   &                     fwd2*zeta(i,j,kold)
!^                                             using background instead
          Drhs(i,j)=h(i,j)+zeta(i,j,kstp)
        END DO
      END DO
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
!^        urhs(i,j)=fwd0*ubar(i,j,kstp)+                                &
!^   &              fwd1*ubar(i,j,kbak)+                                &
!^   &              fwd2*ubar(i,j,kold)
!^                                             using background instead
          urhs(i,j)=ubar(i,j,kstp)
          DUon(i,j)=urhs(i,j)*cff1
        END DO
      END DO
!
      DO j=JV_RANGE
        DO i=IV_RANGE
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
!^        vrhs(i,j)=fwd0*vbar(i,j,kstp)+                                &
!^   &              fwd1*vbar(i,j,kbak)+                                &
!^   &              fwd2*vbar(i,j,kold)
!^                                             using background instead
          vrhs(i,j)=vbar(i,j,kstp)
          DVom(i,j)=vrhs(i,j)*cff1
        END DO
      END DO

#undef IU_RANGE
#undef IV_RANGE
#undef JU_RANGE
#undef JV_RANGE

#if defined DISTRIBUTE && \
    defined UV_ADV     && defined UV_C4ADVECTION && !defined SOLVE3D
!
!  In distributed-memory, the I- and J-ranges are different and a
!  special exchange is done here to avoid having three ghost points
!  for high-order numerical stencils. Notice that a private array is
!  passed below to the exchange routine. It also applies periodic
!  boundary conditions, if appropriate and no partitions in I- or
!  J-directions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          DUon)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          DVom)
      END IF
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    DUon, DVom)
#endif
!
!  Compute integral mass flux across open boundaries and adjust
!  for volume conservation. Compute BASIC STATE value.
!  HGA: Need to resolve 'krhs' index here.
!
      IF (ANY(VolCons(:,ng))) THEN
        CALL obc_flux_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      knew,                                       &
#ifdef MASKING
     &                      umask, vmask,                               &
#endif
     &                      h, om_v, on_u,                              &
     &                      ubar, vbar, zeta)
!
!  Set vertically integrated mass fluxes DUon and DVom along the open
!  boundaries in such a way that the integral volume is conserved.
!
        CALL set_DUV_bc_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        krhs,                                     &
#ifdef MASKING
     &                        umask, vmask,                             &
#endif
     &                        om_v, on_u,                               &
     &                        ubar, vbar,                               &
     &                        Drhs, DUon, DVom)
      END IF
!
!-----------------------------------------------------------------------
!  Compute BASIC STATE fields associated with pressure gradient and
!  time-stepping of adjoint free-surface, "zeta_new".
!-----------------------------------------------------------------------
!
!  Notice that the new local free-surface is allocated so it can be
!  passed as an argumment to "zetabc_local". An automatic array cannot
!  be used here because of weird memory problems.
!
      allocate ( ad_zeta_new(IminS:ImaxS,JminS:JmaxS) )
      ad_zeta_new = 0.0_r8
!
!  Compute "zeta_new" at new time step and interpolate it half-step
!  backward, "zwrk" for the subsequent computation of the adjoint
!  barotropic pressure gradient. Here, we use the BASIC STATE values.
!  Thus, the nonlinear correction to the pressure-gradient term from
!  "kstp" to "knew" is not needed for Forward-Euler to Forward-Backward
!  steps (PGF_FB_CORRECTION method).
!
!  Get background zeta_new from BASIC state. Notice the I- and J-range
!  used to avoid calling nonlinear 'zetabc_local' routine.
!
      DO j=JR_RANGE
        DO i=IR_RANGE
          zeta_new(i,j)=zeta(i,j,knew)
#ifdef MASKING
          zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
# ifdef WET_DRY_NOT_YET
!^        zeta_new(i,j)=zeta_new(i,j)+                                  &
!^   &                  (Dcrit(ng)-h(i,j))*(1.0_r8-rmask(i,j))
# endif
#endif
          Dnew(i,j)=h(i,j)+zeta_new(i,j)
          Dnew_rd(i,j)=Dnew(i,j)
          Dstp(i,j)=h(i,j)+zeta(i,j,kstp)
        END DO
      END DO

#undef IR_RANGE
#undef JR_RANGE
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!^        zeta_new(i,j)=zeta(i,j,kstp)+                                 &
!^   &                  dtfast(ng)*pm(i,j)*pn(i,j)*                     &
!^   &                  (DUon(i,j)-DUon(i+1,j)+                         &
!^   &                   DVom(i,j)-DVom(i,j+1))
#ifdef MASKING
!^        zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
# ifdef WET_DRY_NOT_YET
!!        zeta_new(i,j)=zeta_new(i,j)+                                  &
!!   &                  (Dcrit(ng)-h(i,j))*(1.0_r8-rmask(i,j))
# endif
#endif
!^
!                                              using background instead
          zwrk(i,j)=bkw_new*zeta_new(i,j)+                              &
     &              bkw0*zeta(i,j,kstp)+                                &
     &              bkw1*zeta(i,j,kbak)+                                &
     &              bkw2*zeta(i,j,kold)
#if defined VAR_RHO_2D && defined SOLVE3D
          rzeta(i,j)=(1.0_r8+rhoS(i,j))*zwrk(i,j)
          rzeta2(i,j)=rzeta(i,j)*zwrk(i,j)
          rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
#else
          rzeta(i,j)=zwrk(i,j)
          rzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
#endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Save adjoint 2D solution at knew index for IO purposes.
!-----------------------------------------------------------------------
!
#ifdef SOLVE3D
      IF (iif(ng).eq.nfast(ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta_sol(i,j)=ad_zeta(i,j,knew)
          END DO
          DO i=Istr,IendR
            ad_ubar_sol(i,j)=ad_ubar(i,j,knew)
          END DO
          IF (j.ge.Jstr) THEN
            DO i=IstrR,IendR
              ad_vbar_sol(i,j)=ad_vbar(i,j,knew)
            END DO
          END IF
        END DO
      END IF
#else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta_sol(i,j)=ad_zeta(i,j,knew)
        END DO
        DO i=Istr,IendR
          ad_ubar_sol(i,j)=ad_ubar(i,j,knew)
        END DO
        IF (j.ge.Jstr) THEN
          DO i=IstrR,IendR
            ad_vbar_sol(i,j)=ad_vbar(i,j,knew)
          END DO
        END IF
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Adjoint of exchange halo tile information.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange2d (ng, tile, iTLM, 3,                            &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_zeta(:,:,knew),                            &
!^   &                    tl_ubar(:,:,knew),                            &
!^   &                    tl_vbar(:,:,knew))
!^
      CALL ad_mp_exchange2d (ng, tile, iADM, 3,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,knew),                         &
     &                       ad_ubar(:,:,knew),                         &
     &                       ad_vbar(:,:,knew))
#endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!^      CALL exchange_v2d_tile (ng, tile,                               &
!^   &                          LBi, UBi, LBj, UBj,                     &
!^   &                          tl_vbar(:,:,knew))
!^
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,knew))
!^      CALL exchange_u2d_tile (ng, tile,                               &
!^   &                          LBi, UBi, LBj, UBj,                     &
!^   &                          tl_ubar(:,:,knew))
!^
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,knew))
!^      CALL exchange_r2d_tile (ng, tile,                               &
!^   &                          LBi, UBi, LBj, UBj,                     &
!^   &                          tl_zeta(:,:,knew))
!^
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,knew))
      END IF

#ifdef NESTING
# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  If nesting and after all fast time steps are completed, adjoint of
!  exchange halo information to time averaged fields.
!-----------------------------------------------------------------------
!
      IF (iif(ng).eq.nfast(ng)) THEN
!
!  In nesting applications with refinement grids, we need to exchange
!  the DU_avg2 and DV_avg2 fluxes boundary information for the case
!  that a contact point is at a tile partition. Notice that in such
!  cases, we need i+1 and j+1 values for spatial/temporal interpolation.
!
#  ifdef DISTRIBUTE
!^      CALL mp_exchange2d (ng, tile, iTLM, 2,                          &
!^   &                      LBi, UBi, LBj, UBj,                         &
!^   &                      NghostPoints,                               &
!^   &                      EWperiodic(ng), NSperiodic(ng),             &
!^   &                      tl_DU_avg2, tl_DV_avg2)
!^
        CALL ad_mp_exchange2d (ng, tile, iADM, 2,                       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_DU_avg2, ad_DV_avg2)
!^      CALL mp_exchange2d (ng, tile, iTLM, 3,                          &
!^   &                      LBi, UBi, LBj, UBj,                         &
!^   &                      NghostPoints,                               &
!^   &                      EWperiodic(ng), NSperiodic(ng),             &
!^   &                      tl_Zt_avg1, tl_DU_avg1, tl_DV_avg1)
!^
        CALL ad_mp_exchange2d (ng, tile, iADM, 3,                       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_Zt_avg1, ad_DU_avg1, ad_DV_avg1)
!
#  endif
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!^        CALL exchange_v2d_tile (ng, tile,                             &
!^   &                            LBi, UBi, LBj, UBj,                   &
!^   &                            tl_DV_avg2)
!^
          CALL ad_exchange_v2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_DV_avg2)
!^        CALL exchange_u2d_tile (ng, tile,                             &
!^   &                            LBi, UBi, LBj, UBj,                   &
!^   &                            tl_DU_avg2)
!^
          CALL ad_exchange_u2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_DU_avg2)
!^        CALL exchange_v2d_tile (ng, tile,                             &
!^   &                            LBi, UBi, LBj, UBj,                   &
!^   &                            tl_DV_avg1)
!^
          CALL ad_exchange_v2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_DV_avg1)
!^        CALL exchange_u2d_tile (ng, tile,                             &
!^   &                            LBi, UBi, LBj, UBj,                   &
!^   &                            tl_DU_avg1)
!^
          CALL ad_exchange_u2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_DU_avg1)
!^        CALL exchange_r2d_tile (ng, tile,                             &
!^   &                            LBi, UBi, LBj, UBj,                   &
!^   &                            tl_Zt_avg1)
!^
          CALL ad_exchange_r2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_Zt_avg1)
        END IF

      END IF
# else
!
!  In nesting applications with refinement grids, we need to exchange
!  the DU_flux and DV_flux fluxes boundary information for the case
!  that a contact point is at a tile partition. Notice that in such
!  cases, we need i+1 and j+1 values for spatial/temporal interpolation.
!
#  ifdef DISTRIBUTE
!^    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_DU_flux, tl_DV_flux)
!^
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_DU_flux, ad_DV_flux)
!
#  endif
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!^      CALL exchange_v2d_tile (ng, tile,                               &
!^   &                          LBi, UBi, LBj, UBj,                     &
!^   &                          tl_DV_flux)
!^
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_DV_flux)
!^      CALL exchange_u2d_tile (ng, tile,                               &
!^   &                          LBi, UBi, LBj, UBj,                     &
!^   &                          tl_DU_flux)
!^
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_DU_flux)
      END IF
# endif
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint replace the new free-surface zeta(:,:,knew) with it fast
!  time-averaged value, Zt_avg1 at the of the last 2D time step. Recall
!  this is state variable is the one that communicates with the 3D
!  kernel.
!-----------------------------------------------------------------------
!
      IF (iif(ng).eq.nfast(ng)) THEN
!^      CALL tl_set_depth (ng, tile, iTLM)
!^
        CALL ad_set_depth (ng, tile, iADM)
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
!^          tl_zeta(i,j,knew)=tl_Zt_avg1(i,j)
!^
            ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+ad_zeta(i,j,knew)
            ad_zeta(i,j,knew)=0.0_r8
          END DO
        END DO
      END IF
#endif

#ifdef WET_DRY_NOT_YET
!
!-----------------------------------------------------------------------
!  Adjoint of compute new wet/dry masks.
!-----------------------------------------------------------------------
!
!^    CALL wetdry_tile (ng, tile,                                       &
!^   &                  LBi, UBi, LBj, UBj,                             &
!^   &                  IminS, ImaxS, JminS, JmaxS,                     &
# ifdef MASKING
!^   &                  pmask, rmask, umask, vmask,                     &
# endif
!^   &                  h, zeta(:,:,knew),                              &
# ifdef SOLVE3D
!^   &                  DU_avg1, DV_avg1,                               &
!^   &                  rmask_wet_avg,                                  &
# endif
!^   &                  pmask_wet, pmask_full,                          &
!^   &                  rmask_wet, rmask_full,                          &
!^   &                  umask_wet, umask_full,                          &
!^   &                  vmask_wet, vmask_full)
!^
!^  HGA: Need the ADM code here.
!^
#endif
!
!-----------------------------------------------------------------------
!  Apply adjoint momentum transport point sources (like river runoff),
!  if any.
!
!    Dsrc(is) = 0,  flow across grid cell u-face (positive or negative)
!    Dsrc(is) = 1,  flow across grid cell v-face (positive or negative)
!-----------------------------------------------------------------------
!
      IF (LuvSrc(ng)) THEN
        DO is=1,Nsrc(ng)
          i=SOURCES(ng)%Isrc(is)
          j=SOURCES(ng)%Jsrc(is)
          IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
#if defined NESTING && !defined SOLVE3D
!^            tl_DU_flux(i,j)=SOURCES(ng)%tl_Qbar(is)
!^
              SOURCES(ng)%ad_Qbar(is)=SOURCES(ng)%ad_Qbar(is)+          &
     &                                ad_DU_flux(i,j)
              ad_DU_flux(i,j)=0.0_r8
#endif
#ifdef SOLVE3D
!^            tl_DU_avg1(i,j)=SOURCES(ng)%tl_Qbar(is)
!^
              SOURCES(ng)%ad_Qbar(is)=SOURCES(ng)%ad_Qbar(is)+          &
     &                                ad_DU_avg1(i,j)
              ad_DU_avg1(i,j)=0.0_r8
#endif
              cff=1.0_r8/(on_u(i,j)*                                    &
     &                    0.5_r8*(Dnew(i-1,j)+Dnew(i,j)))
!^            tl_ubar(i,j,knew)=SOURCES(ng)%tl_Qbar(is)*cff+            &
!^   &                          SOURCES(ng)%Qbar(is)*tl_cff
!^
              SOURCES(ng)%ad_Qbar(is)=SOURCES(ng)%ad_Qbar(is)+          &
     &                                cff*ad_ubar(i,j,knew)
              ad_cff=ad_cff+                                            &
     &               SOURCES(ng)%Qbar(is)*ad_ubar(i,j,knew)

              ad_ubar(i,j,knew)=0.0_r8
!^            tl_cff=-cff*cff*on_u(i,j)*                                &
!^   &               0.5_r8*(tl_Dnew(i-1,j)+tl_Dnew(i,j))
!^
              adfac=-cff*cff*on_u(i,j)*0.5_r8*ad_cff
              ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac
              ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac
              ad_cff=0.0_r8
            ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
#if defined NESTING && !defined SOLVE3D
!^            tl_DV_flux(i,j)=SOURCES(ng)%tl_Qbar(is)
!^
              SOURCES(ng)%ad_Qbar(is)=SOURCES(ng)%ad_Qbar(is)+          &
     &                                ad_DV_flux(i,j)
              ad_DV_flux(i,j)=0.0_r8
#endif
#ifdef SOLVE3D
!^            tl_DV_avg1(i,j)=SOURCES(ng)%tl_Qbar(is)
!^
              SOURCES(ng)%ad_Qbar(is)=SOURCES(ng)%ad_Qbar(is)+          &
     &                                ad_DV_avg1(i,j)
              ad_DV_avg1(i,j)=0.0_r8
#endif
              cff=1.0_r8/(om_v(i,j)*                                    &
     &                    0.5_r8*(Dnew(i,j-1)+Dnew(i,j)))
!^            tl_vbar(i,j,knew)=SOURCES(ng)%tl_Qbar(is)*cff+            &
!^   &                          SOURCES(ng)%Qbar(is)*tl_cff
!^
              SOURCES(ng)%ad_Qbar(is)=SOURCES(ng)%ad_Qbar(is)+          &
     &                                cff*ad_vbar(i,j,knew)
              ad_cff=ad_cff+                                            &
     &               SOURCES(ng)%Qbar(is)*ad_vbar(i,j,knew)
              ad_vbar(i,j,knew)=0.0_r8
!^            tl_cff=-cff*cff*om_v(i,j)*                                &
!^   &               0.5_r8*(tl_Dnew(i,j-1)+tl_Dnew(i,j))
!^
              adfac=-cff*cff*om_v(i,j)*0.5_r8*ad_cff
              ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac
              ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac
              ad_cff=0.0_r8
            END IF
          END IF
        END DO
      END IF

#if defined SOLVE3D || (defined NESTING && !defined SOLVE3D)
!
!  Set adjoint barotropic fluxes along physical boundaries.
!
# ifdef SOLVE3D
      cff1=0.5*weight(1,iif(ng),ng)
!
# endif
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrR,IendR
# if defined NESTING && !defined SOLVE3D
!^          tl_DV_flux(i,Jend+1)=0.5_r8*om_v(i,Jend+1)*                 &
!^   &                           ((Dnew(i,Jend+1)+                      &
!^   &                             Dnew(i,Jend  ))*                     &
!^   &                            tl_vbar(i,Jend+1,knew)+               &
!^   &                            (tl_Dnew(i,Jend+1)+                   &
!^   &                             tl_Dnew(i,Jend  ))*                  &
!^   &                            vbar(i,Jend+1,knew))
!^
            adfac=0.5_r8*om_v(i,Jend+1)*ad_DV_flux(i,Jend+1)
            adfac1=adfac1*vbar(i,Jend+1,knew)
            ad_vbar(i,Jend+1,knew)=ad_vbar(i,Jend+1,knew)+              &
     &                             (Dnew(i,Jend+1)+                     &
     &                              Dnew(i,Jend  ))*adfac
            ad_Dnew(i,Jend  )=ad_Dnew(i,Jend  )+adfac1
            ad_Dnew(i,Jend+1)=ad_Dnew(i,Jend+1)+adfac1
            ad_DV_flux(i,Jend+1)=0.0_r8
# else
!^          tl_DV_avg1(i,Jend+1)=tl_DV_avg1(i,Jend+1)+                  &
!^   &                           cff1*om_v(i,Jend+1)*                   &
!^   &                           ((Dnew(i,Jend+1)+                      &
!^   &                             Dnew(i,Jend  ))*                     &
!^   &                            tl_vbar(i,Jend+1,knew)+               &
!^   &                            (tl_Dnew(i,Jend+1)+                   &
!^   &                             tl_Dnew(i,Jend  ))*                  &
!^   &                            vbar(i,Jend+1,knew))
!^
            adfac=cff1*om_v(i,Jend+1)*ad_DV_avg1(i,Jend+1)
            adfac1=adfac*vbar(i,Jend+1,knew)
            ad_vbar(i,Jend+1,knew)=ad_vbar(i,Jend+1,knew)+              &
     &                             (Dnew(i,Jend+1)+                     &
     &                              Dnew(i,Jend  ))*adfac
            ad_Dnew(i,Jend  )=ad_Dnew(i,Jend  )+adfac1
            ad_Dnew(i,Jend+1)=ad_Dnew(i,Jend+1)+adfac1
# endif
          END DO
          DO i=IstrU,Iend
# if defined NESTING && !defined SOLVE3D
!^          tl_DU_flux(i,Jend+1)=0.5_r8*on_u(i,Jend+1)*                 &
!^   &                           ((Dnew(i  ,Jend+1)+                    &
!^   &                             Dnew(i-1,Jend+1))*                   &
!^   &                            tl_ubar(i,Jend+1,knew)+               &
!^   &                            (tl_Dnew(i  ,Jend+1)+                 &
!^   &                             tl_Dnew(i-1,Jend+1))*                &
!^   &                            ubar(i,Jend+1,knew))
!^
            adfac=0.5_r8*on_u(i,Jend+1)*ad_DU_flux(i,Jend+1)
            adfac1=adfac*ubar(i,Jend+1,knew)
            ad_ubar(i,Jend+1,knew)=ad_ubar(i,Jend+1,knew)+              &
     &                             (Dnew(i  ,Jend+1)+                   &
     &                              Dnew(i-1,Jend+1))*adfac
            ad_Dnew(i-1,Jend+1)=ad_Dnew(i-1,Jend+1)+adfac1
            ad_Dnew(i  ,Jend+1)=ad_Dnew(i  ,Jend+1)+adfac1
            ad_DU_flux(i,Jend+1)=0.0_r8
# else
!^          tl_DU_avg1(i,Jend+1)=tl_DU_avg1(i,Jend+1)+                  &
!^   &                           cff1*on_u(i,Jend+1)*                   &
!^   &                           ((Dnew(i  ,Jend+1)+                    &
!^   &                             Dnew(i-1,Jend+1))*                   &
!^   &                            tl_ubar(i,Jend+1,knew)+               &
!^   &                            (tl_Dnew(i  ,Jend+1)+                 &
!^   &                             tl_Dnew(i-1,Jend+1))*                &
!^   &                            ubar(i,Jend+1,knew))
!^
            adfac=cff1*on_u(i,Jend+1)*ad_DU_avg1(i,Jend+1)
            adfac1=adfac*ubar(i,Jend+1,knew)
            ad_ubar(i,Jend+1,knew)=ad_ubar(i,Jend+1,knew)+              &
     &                             (Dnew(i  ,Jend+1)+                   &
     &                              Dnew(i-1,Jend+1))*adfac
            ad_Dnew(i-1,Jend+1)=ad_Dnew(i-1,Jend+1)+adfac1
            ad_Dnew(i  ,Jend+1)=ad_Dnew(i  ,Jend+1)+adfac1
# endif
          END DO
        END IF
      END IF

      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrR,IendR
# if defined NESTING && !defined SOLVE3D
!^          tl_DV_flux(i,JstrV-1)=0.5_r8*om_v(i,JstrV-1)*               &
!^   &                            ((Dnew(i,JstrV-1)+                    &
!^   &                              Dnew(i,JstrV-2))*                   &
!^   &                             tl_vbar(i,JstrV-1,knew)+             &
!^   &                             (tl_Dnew(i,JstrV-1)+                 &
!^   &                              tl_Dnew(i,JstrV-2))*                &
!^   &                             vbar(i,JstrV-1,knew))
!^
            adfac=0.5_r8*om_v(i,JstrV-1)*ad_DV_flux(i,JstrV-1)
            adfac1=adfac*vbar(i,JstrV-1,knew)
            ad_vbar(i,JstrV-1,knew)=ad_vbar(i,JstrV-1,knew)+            &
     &                              (Dnew(i,JstrV-1)+                   &
     &                               Dnew(i,JstrV-2))*adfac
            ad_Dnew(i,JstrV-2)=ad_Dnew(i,JstrV-2)+adfac1
            ad_Dnew(i,JstrV-1)=ad_Dnew(i,JstrV-1)+adfac1
            ad_DV_flux(i,JstrV-1)=0.0_r8
# else
!^          tl_DV_avg1(i,JstrV-1)=tl_DV_avg1(i,JstrV-1)+                &
!^   &                            cff1*om_v(i,JstrV-1)*                 &
!^   &                            ((Dnew(i,JstrV-1)+                    &
!^   &                              Dnew(i,JstrV-2))*                   &
!^   &                             tl_vbar(i,JstrV-1,knew)+             &
!^   &                             (tl_Dnew(i,JstrV-1)+                 &
!^   &                              tl_Dnew(i,JstrV-2))*                &
!^   &                             vbar(i,JstrV-1,knew))
!^
            adfac=cff1*om_v(i,JstrV-1)*ad_DV_avg1(i,JstrV-1)
            adfac1=adfac*vbar(i,JstrV-1,knew)
            ad_vbar(i,JstrV-1,knew)=ad_vbar(i,JstrV-1,knew)+            &
     &                              (Dnew(i,JstrV-1)+                   &
     &                               Dnew(i,JstrV-2))*adfac
            ad_Dnew(i,JstrV-2)=ad_Dnew(i,JstrV-2)+adfac1
            ad_Dnew(i,JstrV-1)=ad_Dnew(i,JstrV-1)+adfac1
# endif
          END DO
          DO i=IstrU,Iend
# if defined NESTING && !defined SOLVE3D
!^          tl_DU_flux(i,Jstr-1)=0.5_r8*on_u(i,Jstr-1)*                 &
!^   &                           ((Dnew(i  ,Jstr-1)+                    &
!^   &                             Dnew(i-1,Jstr-1))*                   &
!^   &                            tl_ubar(i,Jstr-1,knew)+               &
!^   &                            (tl_Dnew(i  ,Jstr-1)+                 &
!^   &                             tl_Dnew(i-1,Jstr-1))*                &
!^   &                            ubar(i,Jstr-1,knew))
!^
            adfac=0.5_r8*on_u(i,Jstr-1)*ad_DU_flux(i,Jstr-1)
            adfac1=adfac*ubar(i,Jstr-1,knew)
            ad_ubar(i,Jstr-1,knew)=ad_ubar(i,Jstr-1,knew)+              &
     &                             (Dnew(i  ,Jstr-1)+                   &
     &                              Dnew(i-1,Jstr-1))*adfac
            ad_Dnew(i-1,Jstr-1)=ad_Dnew(i-1,Jstr-1)+adfac1
            ad_Dnew(i  ,Jstr-1)=ad_Dnew(i  ,Jstr-1)+adfac1
            ad_DU_flux(i,Jstr-1)=0.0_r8
# else
!^          tl_DU_avg1(i,Jstr-1)=tl_DU_avg1(i,Jstr-1)+                  &
!^   &                           cff1*on_u(i,Jstr-1)*                   &
!^   &                           ((Dnew(i  ,Jstr-1)+                    &
!^   &                             Dnew(i-1,Jstr-1))*                   &
!^   &                            tl_ubar(i,Jstr-1,knew)+               &
!^   &                            (tl_Dnew(i  ,Jstr-1)+                 &
!^   &                             tl_Dnew(i-1,Jstr-1))*                &
!^   &                            ubar(i,Jstr-1,knew))
!^
# endif
            adfac=cff1*on_u(i,Jstr-1)*ad_DU_avg1(i,Jstr-1)
            adfac1=adfac*ubar(i,Jstr-1,knew)
            ad_ubar(i,Jstr-1,knew)=ad_ubar(i,Jstr-1,knew)+              &
     &                             (Dnew(i  ,Jstr-1)+                   &
     &                              Dnew(i-1,Jstr-1))*adfac
            ad_Dnew(i-1,Jstr-1)=ad_Dnew(i-1,Jstr-1)+adfac1
            ad_Dnew(i  ,Jstr-1)=ad_Dnew(i  ,Jstr-1)+adfac1
          END DO
        END IF
      END IF

      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
# if defined NESTING && !defined SOLVE3D
!^          tl_DV_flux(Iend+1,j)=0.5_r8*om_v(Iend+1,j)*                 &
!^   &                           ((Dnew(Iend+1,j  )+                    &
!^   &                             Dnew(Iend+1,j-1))*                   &
!^   &                            tl_vbar(Iend+1,j,knew)+               &
!^   &                            (tl_Dnew(Iend+1,j  )+                 &
!^   &                             tl_Dnew(Iend+1,j-1))*                &
!^   &                            vbar(Iend+1,j,knew))
!^
            adfac=0.5_r8*om_v(Iend+1,j)*ad_DV_flux(Iend+1,j)
            adfac1=adfac*vbar(Iend+1,j,knew)
            ad_vbar(Iend+1,j,knew)=ad_vbar(Iend+1,j,knew)+              &
     &                             (Dnew(Iend+1,j  )+                   &
     &                              Dnew(Iend+1,j-1))*adfac
            ad_Dnew(Iend+1,j-1)=ad_Dnew(Iend+1,j-1)+adfac1
            ad_Dnew(Iend+1,j  )=ad_Dnew(Iend+1,j  )+adfac1
            ad_DV_flux(Iend+1,j)=0.0_r8
# else
!^          tl_DV_avg1(Iend+1,j)=tl_DV_avg1(Iend+1,j)+                  &
!^   &                           cff1*om_v(Iend+1,j)*                   &
!^   &                           ((Dnew(Iend+1,j  )+                    &
!^   &                             Dnew(Iend+1,j-1))*                   &
!^   &                            tl_vbar(Iend+1,j,knew)+               &
!^   &                            (tl_Dnew(Iend+1,j  )+                 &
!^   &                             tl_Dnew(Iend+1,j-1))*                &
!^   &                            vbar(Iend+1,j,knew))
!^
            adfac=cff1*om_v(Iend+1,j)*ad_DV_avg1(Iend+1,j)
            adfac1=adfac*vbar(Iend+1,j,knew)
            ad_vbar(Iend+1,j,knew)=ad_vbar(Iend+1,j,knew)+              &
     &                             (Dnew(Iend+1,j  )+                   &
     &                              Dnew(Iend+1,j-1))*adfac
            ad_Dnew(Iend+1,j-1)=ad_Dnew(Iend+1,j-1)+adfac1
            ad_Dnew(Iend+1,j  )=ad_Dnew(Iend+1,j  )+adfac1
# endif
          END DO
          DO j=JstrR,JendR
# if defined NESTING && !defined SOLVE3D
!^          tl_DU_flux(Iend+1,j)=0.5_r8*on_u(Iend+1,j)*                 &
!^   &                           ((Dnew(Iend+1,j)+                      &
!^   &                             Dnew(Iend  ,j))*                     &
!^   &                            tl_ubar(Iend+1,j,knew)+               &
!^   &                            (tl_Dnew(Iend+1,j)+                   &
!^   &                             tl_Dnew(Iend  ,j))*                  &
!^   &                            ubar(Iend+1,j,knew))
!^
            adfac=0.5_r8*on_u(Iend+1,j)*ad_DU_flux(Iend+1,j)
            adfac1=adfac*ubar(Iend+1,j,knew)
            ad_ubar(Iend+1,j,knew)=ad_ubar(Iend+1,j,knew)+              &
     &                             (Dnew(Iend+1,j)+                     &
     &                              Dnew(Iend  ,j))*adfac
            ad_Dnew(Iend  ,j)=ad_Dnew(Iend  ,j)+adfac1
            ad_Dnew(Iend+1,j)=ad_Dnew(Iend+1,j)+adfac1
            ad_DU_flux(Iend+1,j)=0.0_r8
# else
!^          tl_DU_avg1(Iend+1,j)=tl_DU_avg1(Iend+1,j)+                  &
!^   &                           cff1*on_u(Iend+1,j)*                   &
!^   &                           ((Dnew(Iend+1,j)+                      &
!^   &                             Dnew(Iend  ,j))*                     &
!^   &                            tl_ubar(Iend+1,j,knew)+               &
!^   &                            (tl_Dnew(Iend+1,j)+                   &
!^   &                             tl_Dnew(Iend  ,j))*                  &
!^   &                            ubar(Iend+1,j,knew))
!^
            adfac=cff1*on_u(Iend+1,j)*ad_DU_avg1(Iend+1,j)
            adfac1=adfac*ubar(Iend+1,j,knew)
            ad_ubar(Iend+1,j,knew)=ad_ubar(Iend+1,j,knew)+              &
     &                             (Dnew(Iend+1,j)+                     &
     &                              Dnew(Iend  ,j))*adfac
            ad_Dnew(Iend  ,j)=ad_Dnew(Iend  ,j)+adfac1
            ad_Dnew(Iend+1,j)=ad_Dnew(Iend+1,j)+adfac1
# endif
          END DO
        END IF
      END IF

      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
# if defined NESTING && !defined SOLVE3D
!^          tl_DV_flux(Istr-1,j)=0.5_r8*om_v(Istr-1,j)*                 &
!^   &                           ((Dnew(Istr-1,j  )+                    &
!^   &                             Dnew(Istr-1,j-1))*                   &
!^   &                            tl_vbar(Istr-1,j,knew)+               &
!^   &                            (tl_Dnew(Istr-1,j  )+                 &
!^   &                             tl_Dnew(Istr-1,j-1))*                &
!^   &                            vbar(Istr-1,j,knew))
!^
            adfac=0.5_r8*om_v(Istr-1,j)*ad_DV_flux(Istr-1,j)
            adfac1=adfac*vbar(Istr-1,j,knew)
            ad_vbar(Istr-1,j,knew)=ad_vbar(Istr-1,j,knew)+              &
     &                             (Dnew(Istr-1,j  )+                   &
     &                              Dnew(Istr-1,j-1))*adfac
            ad_Dnew(Istr-1,j-1)=ad_Dnew(Istr-1,j-1)+adfac1
            ad_Dnew(Istr-1,j  )=ad_Dnew(Istr-1,j  )+adfac1
            ad_DV_flux(Istr-1,j)=0.0_r8
# else
!^          tl_DV_avg1(Istr-1,j)=tl_DV_avg1(Istr-1,j)+                  &
!^   &                           cff1*om_v(Istr-1,j)*                   &
!^   &                           ((Dnew(Istr-1,j  )+                    &
!^   &                             Dnew(Istr-1,j-1))*                   &
!^   &                            tl_vbar(Istr-1,j,knew)+               &
!^   &                            (tl_Dnew(Istr-1,j  )+                 &
!^   &                             tl_Dnew(Istr-1,j-1))*                &
!^   &                            vbar(Istr-1,j,knew))
!^
            adfac=cff1*om_v(Istr-1,j)*ad_DV_avg1(Istr-1,j)
            adfac1=adfac*vbar(Istr-1,j,knew)
            ad_vbar(Istr-1,j,knew)=ad_vbar(Istr-1,j,knew)+              &
     &                             (Dnew(Istr-1,j  )+                   &
     &                              Dnew(Istr-1,j-1))*adfac
            ad_Dnew(Istr-1,j-1)=ad_Dnew(Istr-1,j-1)+adfac1
            ad_Dnew(Istr-1,j  )=ad_Dnew(Istr-1,j  )+adfac1
# endif
          END DO
          DO j=JstrR,JendR
# if defined NESTING && !defined SOLVE3D
!^          tl_DU_flux(IstrU-1,j)=0.5_r8*on_u(IstrU-1,j)*               &
!^   &                            ((Dnew(IstrU-1,j)+                    &
!^   &                              Dnew(IstrU-2,j))*                   &
!^   &                             tl_ubar(IstrU-1,j,knew)+             &
!^   &                             (tl_Dnew(IstrU-1,j)+                 &
!^   &                              tl_Dnew(IstrU-2,j))*                &
!^   &                             ubar(IstrU-1,j,knew))
!^
            adfac=0.5_r8*on_u(IstrU-1,j)*ad_DU_flux(IstrU-1,j)
            adfac1=adfac*ubar(IstrU-1,j,knew)
            ad_ubar(IstrU-1,j,knew)=ad_ubar(IstrU-1,j,knew)+            &
     &                              (Dnew(IstrU-1,j)+                   &
     &                               Dnew(IstrU-2,j))*adfac
            ad_Dnew(IstrU-2,j)=ad_Dnew(IstrU-2,j)+adfac1
            ad_Dnew(IstrU-1,j)=ad_Dnew(IstrU-1,j)+adfac1
            ad_DU_flux(IstrU-1,j)=0.0_r8
# else
!^          tl_DU_avg1(IstrU-1,j)=tl_DU_avg1(IstrU-1,j)+                &
!^   &                            cff1*on_u(IstrU-1,j)*                 &
!^   &                            ((Dnew(IstrU-1,j)+                    &
!^   &                              Dnew(IstrU-2,j))*                   &
!^   &                             tl_ubar(IstrU-1,j,knew)+             &
!^   &                             (tl_Dnew(IstrU-1,j)+                 &
!^   &                              tl_Dnew(IstrU-2,j))*                &
!^   &                             ubar(IstrU-1,j,knew))
!^
            adfac=cff1*on_u(IstrU-1,j)*ad_DU_avg1(IstrU-1,j)
            adfac1=adfac*ubar(IstrU-1,j,knew)
            ad_ubar(IstrU-1,j,knew)=ad_ubar(IstrU-1,j,knew)+            &
     &                              (Dnew(IstrU-1,j)+                   &
     &                               Dnew(IstrU-2,j))*adfac
            ad_Dnew(IstrU-2,j)=ad_Dnew(IstrU-2,j)+adfac1
            ad_Dnew(IstrU-1,j)=ad_Dnew(IstrU-1,j)+adfac1
# endif
          END DO
        END IF
      END IF
!
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr-1,IendR
!^          tl_Dnew(i,Jend+1)=tl_h(i,Jend+1)+tl_zeta_new(i,Jend+1)
!^
            ad_h(i,Jend+1)=ad_h(i,Jend+1)+                              &
     &                     ad_Dnew(i,Jend+1)
            ad_zeta_new(i,Jend+1)=ad_zeta_new(i,Jend+1)+                &
     &                            ad_Dnew(i,Jend+1)
            ad_Dnew(i,Jend+1)=0.0_r8
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr-1,IendR
!^          tl_Dnew(i,Jstr-1)=tl_h(i,Jstr-1)+tl_zeta_new(i,Jstr-1)
!^
            ad_h(i,Jstr-1)=ad_h(i,Jstr-1)+                              &
     &                     ad_Dnew(i,Jstr-1)
            ad_zeta_new(i,Jstr-1)=ad_zeta_new(i,Jstr-1)+                &
     &                            ad_Dnew(i,Jstr-1)
            ad_Dnew(i,Jstr-1)=0.0_r8
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr-1,JendR
!^          tl_Dnew(Iend+1,j)=tl_h(Iend+1,j)+tl_zeta_new(Iend+1,j)
!^
            ad_h(Iend+1,j)=ad_h(Iend+1,j)+                              &
     &                     ad_Dnew(Iend+1,j)
            ad_zeta_new(Iend+1,j)=ad_zeta_new(Iend+1,j)+                &
     &                            ad_Dnew(Iend+1,j)
            ad_Dnew(Iend+1,j)=0.0_r8
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr-1,JendR
!^          tl_Dnew(Istr-1,j)=tl_h(Istr-1,j)+tl_zeta_new(Istr-1,j)
!^
            ad_h(Istr-1,j)=ad_h(Istr-1,j)+                              &
     &                     ad_Dnew(Istr-1,j)
            ad_zeta_new(Istr-1,j)=ad_zeta_new(Istr-1,j)+                &
     &                            ad_Dnew(Istr-1,j)
            ad_Dnew(Istr-1,j)=0.0_r8
          END DO
        END IF
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Adjoint of time step 2D momentum equations.
!-----------------------------------------------------------------------
!
!  Compute adjoint integral mass flux across open boundaries and adjust
!  for volume conservation.
!
      IF (ANY(VolCons(:,ng))) THEN
!^      CALL tl_obc_flux_tile (ng, tile,                                &
!^   &                         LBi, UBi, LBj, UBj,                      &
!^   &                         IminS, ImaxS, JminS, JmaxS,              &
!^   &                         knew,                                    &
#ifdef MASKING
!^   &                         umask, vmask,                            &
#endif
!^   &                         h, tl_h, om_v, on_u,                     &
!^   &                         ubar, vbar, zeta,                        &
!^   &                         tl_ubar, tl_vbar, tl_zeta)
!^
        CALL ad_obc_flux_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         knew,                                    &
#ifdef MASKING
     &                         umask, vmask,                            &
#endif
     &                         h, ad_h, om_v, on_u,                     &
     &                         ubar, vbar, zeta,                        &
     &                         ad_ubar, ad_vbar, ad_zeta)

      END IF
!
!  Apply adjoint ateral boundary conditions.
!
!^    CALL tl_v2dbc_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    IminS, ImaxS, JminS, JmaxS,                   &
!^   &                    krhs, kstp, knew,                             &
!^   &                    ubar, vbar, zeta,                             &
!^   &                    tl_ubar, tl_vbar, tl_zeta)
!^
      CALL ad_v2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, knew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
!^    CALL tl_u2dbc_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    IminS, ImaxS, JminS, JmaxS,                   &
!^   &                    krhs, kstp, knew,                             &
!^   &                    ubar, vbar, zeta,                             &
!^   &                    tl_ubar, tl_vbar, tl_zeta)
!^
      CALL ad_u2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, knew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
!
!  Advance 2D momentum components while simultaneously adding them to
!  accumulate fast-time averages to compute barotropic fluxes. Doing so
!  straight away yields a more computationally dense code. However, the
!  fast-time averaged fluxes (DU_avg1 and DV_avg1) are needed both at
!  the interior and physical boundary points. Thus, we need separate
!  loops along the domain boundaries after setting "ubar" and "vbar"
!  lateral boundary conditions. Also, note that bottom drag is treated
!  implicitly:
!
!    Dnew*ubar(:,:,m+1) = Dold*ubar(:,:,m) +
!                         dtfast(ng)*rhs2D(:,:) -
!                         dtfast(ng)*rdrag(:,:)*ubar(:,:,m+1)
!  hence
!
!    ubar(:,:,m+1)=[Dold * ubar(..,m) + dtfast(ng) * rhs2D(:,:)] /
!                  [Dnew + dtfast(ng) * rdrag(:,:)]
!
!    DU_avg1 = DU_avg1 +
!              weight(m+1) * Dnew * ubar(:,:,m+1) * on_u(:,:)
!
!  where it should be noted that Dnew .ne. Dnew + dtfast * rdrag
!
      cff=0.5_r8*dtfast(ng)
#ifdef SOLVE3D
      cff1=0.5_r8*weight(1,iif(ng),ng)
#else
      cff2=2.0_r8*dtfast(ng)
#endif
      DO j=JstrV,Jend
        DO i=Istr,Iend
#if defined NESTING && !defined SOLVE3D
!^        tl_DV_flux(i,j)=0.5_r8*om_v(i,j)*                             &
!^   &                    ((Dnew(i,j)+Dnew(i,j-1))*                     &
!^   &                     tl_vbar(i,j,knew)+                           &
!^   &                     (tl_Dnew(i,j)+tl_Dnew(i,j-1))*               &
!^   &                     vbar(i,j,knew))
!^
          adfac=0.5_r8*om_v(i,j)*ad_DV_flux(i,j)
          adfac1=adfac*vbar(i,j,knew)
          ad_vbar(i,j,knew)=ad_vbar(i,j,knew)+                          &
     &                      (Dnew(i,j)+Dnew(i,j-1))*adfac
          ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac1
          ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac1
          ad_DV_flux(i,j)=0.0_r8
#endif
#ifdef SOLVE3D
!^        tl_DV_avg1(i,j)=tl_DV_avg1(i,j)+                              &
!^   &                    cff1*om_v(i,j)*                               &
!^   &                    ((Dnew(i,j)+Dnew(i,j-1))*                     &
!^   &                     tl_vbar(i,j,knew)+                           &
!^   &                     (tl_Dnew(i,j)+tl_Dnew(i,j-1))*               &
!^   &                     vbar(i,j,knew))
!^
          adfac=cff1*om_v(i,j)*ad_DV_avg1(i,j)
          adfac1=adfac*vbar(i,j,knew)
          ad_vbar(i,j,knew)=ad_vbar(i,j,knew)+                          &
     &                      (Dnew(i,j)+Dnew(i,j-1))*adfac
          ad_Dnew(i,j-1)=ad_Dnew(i,j-1)+adfac1
          ad_Dnew(i,j  )=ad_Dnew(i,j  )+adfac1
#endif
#ifdef WET_DRY_NOT_YET
!^        cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
!^        cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
!^        cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!^        vbar(i,j,knew)=vbar(i,j,knew)*cff7
!^
!^  HGA: ADM code needed here.
!^
#endif
#ifdef MASKING
!^        tl_vbar(i,j,knew)=tl_vbar(i,j,knew)*vmask(i,j)
!^
          ad_vbar(i,j,knew)=ad_vbar(i,j,knew)*vmask(i,j)
#endif
          cff3=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          fac2=1.0_r8/(Dnew_rd(i,j)+Dnew_rd(i,j-1))
!^        tl_vbar(i,j,knew)=tl_fac2*                                    &
!^   &                      ((Dstp(i,j)+Dstp(i,j-1))*vbar(i,j,kstp)+    &
#ifdef SOLVE3D
!^   &                       cff3*(rvbar(i,j)+rvfrc(i,j)))+             &
#else
!^   &                       cff3*rvbar(i,j)+cff2*svstr(i,j))+          &
#endif
!^   &                      fac2*                                       &
!^   &                      ((Dstp(i,j)+Dstp(i,j-1))*                   &
!^   &                       tl_vbar(i,j,kstp)+                         &
!^   &                       (tl_Dstp(i,j)+tl_Dstp(i,j-1))*             &
!^   &                       vbar(i,j,kstp)+                            &
#ifdef SOLVE3D
!^   &                       cff3*(tl_rvbar(i,j)+tl_rvfrc(i,j)))
#else
!^   &                       cff3*tl_rvbar(i,j)+cff2*tl_svstr(i,j))
#endif
!^
          adfac=fac2*ad_vbar(i,j,knew)
          adfac1=adfac*(Dstp(i,j)+Dstp(i,j-1))
          adfac2=adfac*cff3
          adfac3=adfac*vbar(i,j,kstp)
          ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)+adfac1
#ifdef SOLVE3D
          ad_rvbar(i,j)=ad_rvbar(i,j)+adfac2
          ad_rvfrc(i,j)=ad_rvfrc(i,j)+adfac2
#else
          ad_rvbar(i,j)=ad_rvbar(i,j)+adfac2
          ad_svstr(i,j)=ad_svstr(i,j)+cff2*adfac
#endif
          ad_Dstp(i,j-1)=ad_Dstp(i,j-1)+adfac3
          ad_Dstp(i,j  )=ad_Dstp(i,j  )+adfac3
          ad_fac2=ad_fac2+                                              &
     &            ad_vbar(i,j,knew)*                                    &
     &            ((Dstp(i,j)+Dstp(i,j-1))*vbar(i,j,kstp)+              &
#ifdef SOLVE3D
     &             cff3*(rvbar(i,j)+rvfrc(i,j)))
#else
    &              cff3*rvbar(i,j)+cff2*svstr(i,j))
#endif
          ad_vbar(i,j,knew)=0.0_r8
!^        tl_fac2=-fac2*fac2*(tl_Dnew_rd(i,j)+tl_Dnew_rd(i,j-1))
!^
          adfac=-fac2*fac2*ad_fac2
          ad_Dnew_rd(i,j-1)=ad_Dnew_rd(i,j-1)+adfac
          ad_Dnew_rd(i,j  )=ad_Dnew_rd(i,j  )+adfac
          ad_fac2=0.0_r8
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
#if defined NESTING && !defined SOLVE3D
!^        tl_DU_flux(i,j)=0.5_r8*on_u(i,j)*                             &
!^   &                    ((Dnew(i,j)+Dnew(i-1,j))*                     &
!^   &                     tl_ubar(i,j,knew)+                           &
!^   &                     (tl_Dnew(i,j)+tl_Dnew(i-1,j))*               &
!^   &                     ubar(i,j,knew))
!^
          adfac=0.5_r8*on_u(i,j)*ad_DU_flux(i,j)
          adfac1=adfac*ubar(i,j,knew)
          ad_ubar(i,j,knew)=ad_ubar(i,j,knew)+                          &
     &                      (Dnew(i,j)+Dnew(i-1,j))*adfac
          ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac1
          ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac1
          ad_DU_flux(i,j)=0.0_r8
#endif
#ifdef SOLVE3D
!^        tl_DU_avg1(i,j)=tl_DU_avg1(i,j)+                              &
!^   &                    cff1*on_u(i,j)*                               &
!^   &                    ((Dnew(i,j)+Dnew(i-1,j))*                     &
!^   &                     tl_ubar(i,j,knew)+                           &
!^   &                     (tl_Dnew(i,j)+tl_Dnew(i-1,j))*               &
!^   &                     ubar(i,j,knew))
!^
          adfac=cff1*on_u(i,j)*ad_DU_avg1(i,j)
          adfac1=adfac*ubar(i,j,knew)
          ad_ubar(i,j,knew)=ad_ubar(i,j,knew)+                          &
     &                      (Dnew(i,j)+Dnew(i-1,j))*adfac
          ad_Dnew(i-1,j)=ad_Dnew(i-1,j)+adfac1
          ad_Dnew(i  ,j)=ad_Dnew(i  ,j)+adfac1
#endif
#ifdef WET_DRY_NOT_YET
!^        cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
!^        cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
!^        cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!^        ubar(i,j,knew)=ubar(i,j,knew)*cff7
!^
!^  HGA: TLM code needed here.
!^
#endif
#ifdef MASKING
!^        tl_ubar(i,j,knew)=tl_ubar(i,j,knew)*umask(i,j)
!^
          ad_ubar(i,j,knew)=ad_ubar(i,j,knew)*umask(i,j)
#endif
          cff3=cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
          fac1=1.0_r8/(Dnew_rd(i,j)+Dnew_rd(i-1,j))
!^        tl_ubar(i,j,knew)=tl_fac1*                                    &
!^   &                      ((Dstp(i,j)+Dstp(i-1,j))*ubar(i,j,kstp)+    &
#ifdef SOLVE3D
!^   &                       cff3*(rubar(i,j)+rufrc(i,j)))+             &
#else
!^   &                       cff3*rubar(i,j)+cff2*sustr(i,j))+          &
#endif
!^   &                      fac1*                                       &
!^   &                      ((Dstp(i,j)+Dstp(i-1,j))*                   &
!^   &                       tl_ubar(i,j,kstp)+                         &
!^   &                       (tl_Dstp(i,j)+tl_Dstp(i-1,j))*             &
!^   &                       ubar(i,j,kstp)+                            &
#ifdef SOLVE3D
!^   &                       cff3*(tl_rubar(i,j)+tl_rufrc(i,j)))
#else
!^   &                       cff3*tl_rubar(i,j)+cff2*tl_sustr(i,j))
#endif
!>
          adfac=fac1*ad_ubar(i,j,knew)
          adfac1=adfac*(Dstp(i,j)+Dstp(i-1,j))
          adfac2=adfac*cff3
          adfac3=adfac*ubar(i,j,kstp)
          ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)+adfac1
#ifdef SOLVE3D
          ad_rubar(i,j)=ad_rubar(i,j)+adfac2
          ad_rufrc(i,j)=ad_rufrc(i,j)+adfac2
#else
          ad_rubar(i,j)=ad_rubar(i,j)+adfac2
          ad_sustr(i,j)=ad_sustr(i,j)+cff2*adfac
#endif
          ad_Dstp(i-1,j)=ad_Dstp(i-1,j)+adfac3
          ad_Dstp(i  ,j)=ad_Dstp(i  ,j)+adfac3
          ad_fac1=ad_fac1+                                              &
     &            ad_ubar(i,j,knew)*                                    &
     &            ((Dstp(i,j)+Dstp(i-1,j))*ubar(i,j,kstp)+              &
#ifdef SOLVE3D
     &             cff3*(rubar(i,j)+rufrc(i,j)))
#else
     &             cff3*rubar(i,j)+cff2*sustr(i,j))
#endif
          ad_ubar(i,j,knew)=0.0_r8
!^        tl_fac1=-fac1*fac1*(tl_Dnew_rd(i,j)+tl_Dnew_rd(i-1,j))
!^
          adfac=-fac1*fac1*ad_fac1
          ad_Dnew_rd(i-1,j)=ad_Dnew_rd(i-1,j)+adfac
          ad_Dnew_rd(i  ,j)=ad_Dnew_rd(i  ,j)+adfac
          ad_fac1=0.0_r8
        END DO
      END DO

#if defined UV_QDRAG && !defined SOLVE3D
!
!  Adjoint of add quadratic drag term associated in shallow-water
!  applications.
!
!  Here, the SQRT(3) is due to a linear interpolation with second order
!  accuaracy that ensures positive and negative values of the velocity
!  components:
!
!   u^2(i+1/2) = (1/3)*[u(i)*u(i) + u(i)*u(i+1) + u(i+1)*u(i+1)]
!
!  If u(i)=1 and u(i+1)=-1, then u^2(i+1/2)=1/3 as it should be.
!
      cff=dtfast(ng)/SQRT(3.0_r8)
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=ubar(i  ,j,kstp)**2+                                     &
     &         ubar(i+1,j,kstp)**2+                                     &
     &         ubar(i  ,j,kstp)*ubar(i+1,j,kstp)+                       &
     &         vbar(i,j  ,kstp)**2+                                     &
     &         vbar(i,j+1,kstp)**2+                                     &
     &         vbar(i,j  ,kstp)*vbar(i,j+1,kstp)
          cff2=SQRT(cff1)
!^        tl_Dnew_rd(i,j)=tl_Dnew_rd(i,j)+                              &
!^   &                    cff*rdrag2(i,j)*tl_cff2
!^
          ad_cff2=ad_cff2+                                              &
     &            cff*rdrag2(i,j)*ad_Dnew_rd(i,j)
!^        tl_cff2=0.5_r8*tl_cff1/cff2
!^
          ad_cff1=ad_cff1+0.5_r8*ad_cff2/cff2
          ad_cff2=0.0_r8
!^        tl_cff1=2.0_r8*ubar(i  ,j,kstp)*tl_ubar(i  ,j,kstp)+          &
!^   &            2.0_r8*ubar(i+1,j,kstp)*tl_ubar(i+1,j,kstp)+          &
!^   &            tl_ubar(i  ,j,kstp)*ubar(i+1,j,kstp)+                 &
!^   &            tl_ubar(i+1,j,kstp)*ubar(i  ,j,kstp)+                 &
!^   &            2.0_r8*vbar(i,j  ,kstp)*tl_vbar(i,j  ,kstp)+          &
!^   &            2.0_r8*vbar(i,j+1,kstp)*tl_vbar(i,j+1,kstp)+          &
!^   &            tl_vbar(i,j  ,kstp)*vbar(i,j+1,kstp)+                 &
!^   &            tl_vbar(i,j+1,kstp)*vbar(i,j  ,kstp)
!^
          adfac=2.0_r8*ad_cff1
          ad_ubar(i  ,j,kstp)=ad_ubar(i  ,j,kstp)+                      &
     &                        ubar(i  ,j,kstp)*adfac+                   &
     &                        ubar(i+1,j,kstp)*ad_cff1
          ad_ubar(i+1,j,kstp)=ad_ubar(i+1,j,kstp)+                      &
     &                        ubar(i+1,j,kstp)*adfac+                   &
     &                        ubar(i  ,j,kstp)*ad_cff1
          ad_vbar(i,j  ,kstp)=ad_vbar(i,j  ,kstp)+                      &
     &                        vbar(i,j  ,kstp)*adfac+                   &
     &                        vbar(i,j+1,kstp)*ad_cff1
          ad_vbar(i,j+1,kstp)=ad_vbar(i,j+1,kstp)+                      &
     &                        vbar(i,j+1,kstp)*adfac+                   &
     &                        vbar(i,j  ,kstp)*ad_cff1
          ad_cff1=0.0_r8
        END DO
      END DO
#endif
!
!  Adjoint of compute depths.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!^        tl_Dstp(i,j)=tl_h(i,j)+tl_zeta(i,j,kstp)
!^
          ad_h(i,j)=ad_h(i,j)+ad_Dstp(i,j)
          ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_Dstp(i,j)
          ad_Dstp(i,j)=0.0_r8
!^        tl_Dnew_rd(i,j)=tl_Dnew(i,j)
!^
          ad_Dnew(i,j)=ad_Dnew(i,j)+ad_Dnew_rd(i,j)
          ad_Dnew_rd(i,j)=0.0_r8
!^        tl_Dnew(i,j)=tl_h(i,j)+tl_zeta_new(i,j)
!^
          ad_h(i,j)=ad_h(i,j)+ad_Dnew(i,j)
          ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_Dnew(i,j)
          ad_Dnew(i,j)=0.0_r8
        END DO
      END DO

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Coupling between 2D and 3D equations.
!-----------------------------------------------------------------------
!
!  Before the first barotropic time step, arrays "rufrc" and "rvfrc"
!  contain vertical integrals of the 3D right-hand-side terms for the
!  momentum equations (including surface and bottom stresses). During
!  the first barotropic time step, convert them into forcing terms by
!  subtracting the fast-time "rubar" and "rvbar" from them.
!
!  In the predictor-coupled mode, the resultant forcing terms "rufrc"
!  and "rvfrc" are extrapolated forward in time, so they become
!  centered effectively at time n+1/2. This is done using optimized
!  Adams-Bashforth weights. In the code below, rufrc_bak(:,:,nstp) is
!  at (n-1)time step, while rufrc_bak(:,:,3-nstp) is at (n-2). After
!  its use as input, the latter is overwritten by the value at time
!  step "nstp" (mathematically "n") during the next step.
!
!  From now on, the computed forcing terms "rufrc" and "rvfrc" will
!  remain constant during  the fast-time stepping and will be added
!  to "rubar" and "rvbar" during all subsequent barotropic steps.
!
      COUPLED_STEP : IF (FIRST_2D_STEP) THEN
!
!  Predictor coupled barotropic mode: Set coefficients for AB3-like
!  forward-in-time extrapolation of 3D to 2D forcing terms "rufrc" and
!  "rvfrc".
!
        IF (iic(ng).eq.ntstart(ng)) THEN
          cfwd0=1.0_r8
          cfwd1=0.0_r8
          cfwd2=0.0_r8
        ELSE IF (iic(ng).eq.ntstart(ng)+1) THEN
          cfwd0=1.5_r8
          cfwd1=-0.5_r8
          cfwd2=0.0_r8
        ELSE
          cfwd2=0.281105_r8
          cfwd1=-0.5_r8-2.0_r8*cfwd2
          cfwd0=1.5_r8+cfwd2
        END IF
!
        cff1=0.5*g
# if defined VAR_RHO_2D && defined SOLVE3D
        cff2=0.333333333333_r8
# endif
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (j.ge.JstrV) THEN
# ifdef DIAGNOSTICS_UV
!!            DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+                &
!!   &                             rvbar(i,j)
# endif
!^            tl_rvbar(i,j)=tl_rvbar(i,j)+                              &
!^   &                      cff1*om_v(i,j)*                             &
!^   &                      ((tl_h(i,j-1)+                              &
!^   &                        tl_h(i,j  ))*                             &
!^   &                       (rzeta(i,j-1)-                             &
!^   &                        rzeta(i,j  ))+                            &
!^   &                       (h(i,j-1)+                                 &
!^   &                        h(i,j  ))*                                &
!^   &                       (tl_rzeta(i,j-1)-                          &
!^   &                        tl_rzeta(i,j  ))+                         &
# if defined VAR_RHO_2D && defined SOLVE3D
!^   &                       (tl_h(i,j-1)-                              &
!^   &                        tl_h(i,j  ))*                             &
!^   &                       (rzetaSA(i,j-1)+                           &
!^   &                        rzetaSA(i,j  )+                           &
!^   &                        cff2*(rhoA(i,j-1)-                        &
!^   &                              rhoA(i,j  ))*                       &
!^   &                             (zwrk(i,j-1)-                        &
!^   &                              zwrk(i,j  )))+                      &
!^   &                       (h(i,j-1)-                                 &
!^   &                        h(i,j  ))*                                &
!^   &                       (tl_rzetaSA(i,j-1)+                        &
!^   &                        tl_rzetaSA(i,j  )+                        &
!^   &                        cff2*((tl_rhoA(i,j-1)-                    &
!^   &                               tl_rhoA(i,j  ))*                   &
!^   &                              (zwrk(i,j-1)-                       &
!^   &                               zwrk(i,j  ))+                      &
!^   &                              (rhoA(i,j-1)-                       &
!^   &                               rhoA(i,j  ))*                      &
!^   &                              (tl_zwrk(i,j-1)-                    &
!^   &                               tl_zwrk(i,j  ))))+                 &
# endif
!^   &                       (tl_rzeta2(i,j-1)-                         &
!^   &                        tl_rzeta2(i,j  )))
!^
              adfac=cff1*om_v(i,j)*ad_rvbar(i,j)
              adfac1=adfac*(rzeta(i,j-1)-rzeta(i,j  ))
              adfac2=adfac*(h(i,j-1)-h(i,j  ))
              ad_h(i,j-1)=ad_h(i,j-1)+adfac1
              ad_h(i,j  )=ad_h(i,j  )+adfac1
              ad_rzeta(i,j-1)=ad_rzeta(i,j-1)+adfac2
              ad_rzeta(i,j  )=ad_rzeta(i,j  )-adfac2
              ad_rzeta2(i,j-1)=ad_rzeta2(i,j-1)+adfac
              ad_rzeta2(i,j  )=ad_rzeta2(i,j  )-adfac
# if defined VAR_RHO_2D && defined SOLVE3D
              adfac3=adfac*(rzetaSA(i,j-1)+                             &
     &                      rzetaSA(i,j  )+                             &
     &                      cff2*(rhoA(i,j-1)-                          &
     &                            rhoA(i,j  ))*                         &
     &                           (zwrk(i,j-1)-                          &
     &                            zwrk(i,j  )))
              adfac4=adfac2*cff2*(zwrk(i,j-1)-zwrk(i,j))
              adfac5=adfac2*cff2*(rhoA(i,j-1)-rhoA(i,j))
              ad_h(i,j-1)=ad_h(i,j-1)+adfac3
              ad_h(i,j  )=ad_h(i,j  )-adfac3
              ad_rzetaSA(i,j-1)=ad_rzetaSA(i,j-1)+adfac2
              ad_rzetaSA(i,j  )=ad_rzetaSA(i,j  )+adfac2
              ad_rhoA(i,j-1)=ad_rhoA(i,j-1)+adfac4
              ad_rhoA(i,j  )=ad_rhoA(i,j  )-adfac4
              ad_zwrk(i,j-1)=ad_zwrk(i,j-1)+adfac5
              ad_zwrk(i,j  )=ad_zwrk(i,j  )-adfac5
# endif
            END IF
!
            IF (i.ge.IstrU) THEN
# ifdef DIAGNOSTICS_UV
!!            DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+                &
!!   &                             rubar(i,j)
# endif
!^            tl_rubar(i,j)=tl_rubar(i,j)+                              &
!^   &                      cff1*on_u(i,j)*                             &
!^   &                      ((tl_h(i-1,j)+                              &
!^   &                        tl_h(i ,j))*                              &
!^   &                       (rzeta(i-1,j)-                             &
!^   &                        rzeta(i  ,j))+                            &
!^   &                       (h(i-1,j)+                                 &
!^   &                        h(i  ,j))*                                &
!^   &                       (tl_rzeta(i-1,j)-                          &
!^   &                        tl_rzeta(i  ,j))+                         &
# if defined VAR_RHO_2D && defined SOLVE3D
!^   &                       (tl_h(i-1,j)-                              &
!^   &                        tl_h(i  ,j))*                             &
!^   &                       (rzetaSA(i-1,j)+                           &
!^   &                        rzetaSA(i  ,j)+                           &
!^   &                        cff2*(rhoA(i-1,j)-                        &
!^   &                              rhoA(i  ,j))*                       &
!^   &                             (zwrk(i-1,j)-                        &
!^   &                              zwrk(i  ,j)))+                      &
!^   &                       (h(i-1,j)-                                 &
!^   &                        h(i  ,j))*                                &
!^   &                       (tl_rzetaSA(i-1,j)+                        &
!^   &                        tl_rzetaSA(i  ,j)+                        &
!^   &                        cff2*((tl_rhoA(i-1,j)-                    &
!^   &                               tl_rhoA(i  ,j))*                   &
!^   &                              (zwrk(i-1,j)-                       &
!^   &                               zwrk(i  ,j))+                      &
!^   &                              (rhoA(i-1,j)-                       &
!^   &                               rhoA(i  ,j))*                      &
!^   &                              (tl_zwrk(i-1,j)-                    &
!^   &                               tl_zwrk(i  ,j))))+                 &
# endif
!^   &                       (tl_rzeta2(i-1,j)-                         &
!^   &                        tl_rzeta2(i  ,j)))
!^
              adfac=cff1*on_u(i,j)*ad_rubar(i,j)
              adfac1=adfac*(rzeta(i-1,j)-rzeta(i  ,j))
              adfac2=adfac*(h(i-1,j)+h(i  ,j))
              ad_h(i-1,j)=ad_h(i-1,j)+adfac1
              ad_h(i  ,j)=ad_h(i  ,j)+adfac1
              ad_rzeta(i-1,j)=ad_rzeta(i-1,j)+adfac2
              ad_rzeta(i  ,j)=ad_rzeta(i  ,j)-adfac2
              ad_rzeta2(i-1,j)=ad_rzeta2(i-1,j)+adfac
              ad_rzeta2(i  ,j)=ad_rzeta2(i  ,j)-adfac
# if defined VAR_RHO_2D && defined SOLVE3D
              adfac3=adfac*(rzetaSA(i-1,j)+                             &
     &                      rzetaSA(i  ,j)+                             &
     &                        cff2*(rhoA(i-1,j)-                        &
     &                              rhoA(i  ,j))*                       &
     &                             (zwrk(i-1,j)-                        &
     &                              zwrk(i  ,j)))
              adfac4=adfac2*cff2*(zwrk(i-1,j)-zwrk(i,j))
              adfac5=adfac2*cff2*(rhoA(i-1,j)-rhoA(i,j))
              ad_h(i-1,j)=ad_h(i-1,j)+adfac3
              ad_h(i  ,j)=ad_h(i  ,j)-adfac3
              ad_rzetaSA(i-1,j)=ad_rzetaSA(i-1,j)+adfac2
              ad_rzetaSA(i  ,j)=ad_rzetaSA(i  ,j)+adfac2
              ad_rhoA(i-1,j)=ad_rhoA(i-1,j)+adfac4
              ad_rhoA(i  ,j)=ad_rhoA(i  ,j)-adfac4
              ad_zwrk(i-1,j)=ad_zwrk(i-1,j)+adfac5
              ad_zwrk(i  ,j)=ad_zwrk(i  ,j)-adfac5
# endif
            END IF
          END DO
        END DO
!
!  Add correction term to shift pressure-gradient terms from "kstp" to
!  "knew". That is, it converts the first 2D step from Forward-Euler
!  to Forward-Backward (this is PGF_FB_CORRECTION mentioned above).
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
# if defined VAR_RHO_2D && defined SOLVE3D
!^          tl_rzetaSA(i,j)=tl_zwrk(i,j)*                               &
!^   &                      (rhoS(i,j)-rhoA(i,j))+                      &
!^   &                      zwrk(i,j)*                                  &
!^   &                      (tl_rhoS(i,j)-tl_rhoA(i,j))
!^
            adfac=zwrk(i,j)*ad_rzetaSA(i,j)
            ad_zwrk(i,j)=ad_zwrk(i,j)+                                  &
     &                   (rhoS(i,j)-rhoA(i,j))*ad_rzetaSA(i,j)
            ad_rhoS(i,j)=ad_rhoS(i,j)+adfac
            ad_rhoA(i,j)=ad_rhoA(i,j)-adfac
            ad_rzetaSA(i,j)=0.0_r8
!^          tl_rzeta2(i,j)=tl_rzeta(i,j)*                               &
!^   &                     (zeta_new(i,j)+zeta(i,j,kstp))+              &
!^   &                     rzeta(i,j)*                                  &
!^   &                     (tl_zeta_new(i,j)+tl_zeta(i,j,kstp))
!^
            adfac=rzeta(i,j)*ad_rzeta2(i,j)
            ad_rzeta(i,j)=ad_rzeta(i,j)+                                &
     &                    (zeta_new(i,j)+zeta(i,j,kstp))*ad_rzeta2(i,j)
            ad_zeta_new(i,j)=ad_zeta_new(i,j)+adfac
            ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+adfac
            ad_rzeta2(i,j)=0.0_r8
!^          tl_rzeta(i,j)=(1.0_r8+rhoS(i,j))*tl_zwrk(i,j)+              &
!^   &                    tl_rhoS(i,j)*zwrk(i,j)
!^
            ad_zwrk(i,j)=ad_zwrk(i,j)+(1.0_r8+rhoS(i,j))*ad_rzeta(i,j)
            ad_rhoS(i,j)=ad_rhoS(i,j)+zwrk(i,j)*ad_rzeta(i,j)
            ad_rzeta(i,j)=0.0_r8
# else
!^          tl_rzeta2(i,j)=tl_zwrk(i,j)*                                &
!^   &                     (zeta_new(i,j)+zeta(i,j,kstp))+              &
!^   &                     zwrk(i,j)*                                   &
!^   &                     (tl_zeta_new(i,j)+tl_zeta(i,j,kstp))
!^
            adfac=zwrk(i,j)*ad_rzeta2(i,j)
            ad_zwrk(i,j)=ad_zwrk(i,j)+                                  &
     &                   (zeta_new(i,j)+zeta(i,j,kstp))*ad_rzeta2(i,j)
            ad_zeta_new(i,j)=ad_zeta_new(i,j)+adfac
            ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+adfac
            ad_rzeta2(i,j)=0.0_r8
!^          tl_rzeta(i,j)=tl_zwrk(i,j)
!^
            ad_zwrk(i,j)=ad_zwrk(i,j)+ad_rzeta(i,j)
            ad_rzeta(i,j)=0.0_r8
# endif
!^          tl_zwrk(i,j)=tl_zeta_new(i,j)-tl_zeta(i,j,kstp)
!^
            ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_zwrk(i,j)
            ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)-ad_zwrk(i,j)
            ad_zwrk(i,j)=0.0_r8
          END DO
        END DO
!
!  Barotropic mode running predictor stage: forward extrapolation.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (j.ge.JstrV) THEN
!^            tl_rvfrc_bak(i,j,3-nstp)=tl_cff2
!^
              ad_cff2=ad_cff2+ad_rvfrc_bak(i,j,3-nstp)
              ad_rvfrc_bak(i,j,3-nstp)=0.0_r8
!^            tl_rvfrc(i,j)=cfwd0*tl_cff2+                              &
!^   &                      cfwd1*tl_rvfrc_bak(i,j,  nstp)+             &
!^   &                      cfwd2*tl_rvfrc_bak(i,j,3-nstp)
!^
              ad_cff2=ad_cff2+cfwd0*ad_rvfrc(i,j)
              ad_rvfrc_bak(i,j,  nstp)=ad_rvfrc_bak(i,j,  nstp)+        &
     &                                 cfwd1*ad_rvfrc(i,j)
              ad_rvfrc_bak(i,j,3-nstp)=ad_rvfrc_bak(i,j,3-nstp)+        &
     &                                 cfwd2*ad_rvfrc(i,j)
              ad_rvfrc(i,j)=0.0_r8
!^            tl_cff2=tl_rvfrc(i,j)-tl_rvbar(i,j)
!^
              ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_cff2
              ad_rvbar(i,j)=ad_rvbar(i,j)-ad_cff2
              ad_cff2=0.0_r8
            END IF
!
            IF (i.ge.IstrU) THEN
!^            tl_rufrc_bak(i,j,3-nstp)=tl_cff1
!^
              ad_cff1=ad_cff1+ad_rufrc_bak(i,j,3-nstp)
              ad_rufrc_bak(i,j,3-nstp)=0.0_r8
!^            tl_rufrc(i,j)=cfwd0*tl_cff1+                              &
!^   &                      cfwd1*tl_rufrc_bak(i,j,  nstp)+             &
!^   &                      cfwd2*tl_rufrc_bak(i,j,3-nstp)
!^
              ad_cff1=ad_cff1+cfwd0*ad_rufrc(i,j)
              ad_rufrc_bak(i,j,  nstp)=ad_rufrc_bak(i,j,  nstp)+        &
     &                                 cfwd1*ad_rufrc(i,j)
              ad_rufrc_bak(i,j,3-nstp)=ad_rufrc_bak(i,j,3-nstp)+        &
     &                                 cfwd2*ad_rufrc(i,j)
              ad_rufrc(i,j)=0.0_r8
!^            tl_cff1=tl_rufrc(i,j)-tl_rubar(i,j)
!^
              ad_rufrc(i,j)=ad_rufrc(i,j)+ad_cff1
              ad_rubar(i,j)=ad_rubar(i,j)-ad_cff1
              ad_cff1=0.0_r8
            END IF
!
!  Compensate for (cancel) bottom drag terms: at input into step2d
!  "rufrc" and "rvfrc" contain bottom drag terms computed by 3D mode.
!  However, there are no 2D counterparts in "rubar" and "rvbar" because
!  2D bottom drag will be computed implicitly during the final stage of
!  updating ubar(:,:,knew) and vbar(:,:,knew) below.  Note that unlike
!  the other terms, the bottom drag should not be extrapolated forward,
!  if "rufrc" and "rvfrc" are, so this cancelation needs to be done
!  right now rather than at the bottom of this loop.
!
            IF (j.ge.JstrV) THEN
!^            tl_rvfrc(i,j)=tl_rvfrc(i,j)+                              &
!^   &                      0.5_r8*(rdrag(i,j)+rdrag(i,j-1))*           &
!^   &                      om_v(i,j)*on_v(i,j)*tl_vbar(i,j,kstp)
!^
              ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)+                      &
     &                          0.5_r8*(rdrag(i,j)+rdrag(i,j-1))*       &
     &                          om_v(i,j)*on_v(i,j)*ad_rvfrc(i,j)
            END IF

            IF (i.ge.IstrU) THEN
!^            tl_rufrc(i,j)=tl_rufrc(i,j)+                              &
!^   &                      0.5_r8*(rdrag(i,j)+rdrag(i-1,j))*           &
!^   &                      om_u(i,j)*on_u(i,j)*tl_ubar(i,j,kstp)
!^
              ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)+                      &
     &                          0.5_r8*(rdrag(i,j)+rdrag(i-1,j))*       &
     &                          om_u(i,j)*on_u(i,j)*ad_rufrc(i,j)
            END IF
          END DO
        END DO
!
      END IF COUPLED_STEP
#endif
!
!=======================================================================
!  Adjoint of compute right-hand-side for the 2D momentum equations.
!=======================================================================
#ifdef SOLVE3D
!
!  Notice that we are suppressing the computation of momentum advection,
!  Coriolis, and lateral viscosity terms in 3D Applications because
!  these terms are already included in the baroclinic-to-barotropic
!  forcing arrays "rufrc" and "rvfrc". It does not mean we are entirely
!  omitting them, but it is a choice between recomputing them at every
!  barotropic step or keeping them "frozen" during the fast-time
!  stepping.
# ifdef STEP2D_CORIOLIS
!  However, in some coarse grid applications with larger baroclinic
!  timestep (say, DT around 20 minutes or larger), adding the Coriolis
!  term in the barotropic equations is useful since f*DT is no longer
!  small.
# endif
#endif

#if defined UV_VIS2 && !defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint of Add in horizontal harmonic viscosity.
!-----------------------------------------------------------------------
!
!  Compute BASIC STATE total depth at PSI-points.
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
          Drhs_p(i,j)=0.25_r8*(Drhs(i,j  )+Drhs(i-1,j  )+               &
     &                         Drhs(i,j-1)+Drhs(i-1,j-1))
        END DO
      END DO
!
!  Add in harmonic viscosity.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (j.ge.JstrV) THEN
# if defined DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2hvis)=fac1
!!          DiaV2rhs(i,j,M2xvis)= cff1
!!          DiaV2rhs(i,j,M2yvis)=-cff2
# endif
!^          tl_rvbar(i,j)=tl_rvbar(i,j)+tl_fac1
!^
            ad_fac1=ad_fac1+ad_rvbar(i,j)
!^          tl_fac1=tl_cff1-tl_cff2
!^
            ad_cff1=ad_cff1+ad_fac1
            ad_cff2=ad_cff2-ad_fac1
            ad_fac1=0.0_r8
!^          tl_cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*                         &
!^   &              (tl_VFe(i  ,j)-tl_VFe(i,j-1))
!^
            adfac=0.5_r8*(pm(i,j-1)+pm(i,j))*ad_cff2
            ad_VFe(i,j-1)=ad_VFe(i,j-1)-adfac
            ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
            ad_cff2=0.0_r8
!^          tl_cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*                         &
!^   &              (tl_VFx(i+1,j)-tl_VFx(i,j  ))
!^
            adfac=0.5_r8*(pn(i,j-1)+pn(i,j))*ad_cff1
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-adfac
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+adfac
            ad_cff1=0.0_r8
          END IF
!
          IF (i.ge.IstrU) THEN
# if defined DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2hvis)=fac1
!!          DiaU2rhs(i,j,M2xvis)=cff1
!!          DiaU2rhs(i,j,M2yvis)=cff2
# endif
!^          tl_rubar(i,j)=tl_rubar(i,j)+tl_fac1
!^
            ad_fac1=ad_fac1+ad_rubar(i,j)
!^          tl_fac1=tl_cff1+tl_cff2
!^
            ad_cff1=ad_cff1+ad_fac1
            ad_cff2=ad_cff2+ad_fac1
            ad_fac1=0.0_r8
!^          tl_cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*                         &
!^   &              (tl_UFe(i,j+1)-tl_UFe(i  ,j))
!^
            adfac=0.5_r8*(pm(i-1,j)+pm(i,j))*ad_cff2
            ad_UFe(i,j  )=ad_UFe(i,j  )-adfac
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+adfac
            ad_cff2=0.0_r8
!^          tl_cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*                         &
!^   &              (tl_UFx(i,j  )-tl_UFx(i-1,j))
!^
            adfac=0.5_r8*(pn(i-1,j)+pn(i,j))*ad_cff1
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-adfac
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
            ad_cff1=0.0_r8
          END IF
        END DO
      END DO
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
!^        tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
!^        tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
!^
          ad_cff=ad_cff+                                                &
     &           on_p(i,j)*on_p(i,j)*ad_VFx(i,j)+                       &
     &           om_p(i,j)*om_p(i,j)*ad_UFe(i,j)
          ad_VFx(i,j)=0.0_r8
          ad_UFe(i,j)=0.0_r8
#  ifdef WET_DRY_NOT_YET
!^        tl_cff=tl_cff*pmask_wet(i,j)
!^
          ad_cff=ad_cff*pmask_wet(i,j)
#  endif
#  ifdef MASKING
!^        tl_cff=tl_cff*pmask(i,j)
!^
          ad_cff=ad_cff*pmask(i,j)
#  endif
!^        tl_cff=visc2_p(i,j)*0.5_r8*                                   &
!^   &           (tl_Drhs_p(i,j)*                                       &
!^   &            (pmon_p(i,j)*                                         &
!^   &             ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,kstp)-           &
!^   &              (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,kstp))+          &
!^   &             pnom_p(i,j)*                                         &
!^   &             ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,kstp)-           &
!^   &              (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,kstp)))+         &
!^   &            Drhs_p(i,j)*                                          &
!^   &            (pmon_p(i,j)*                                         &
!^   &             ((pn(i  ,j-1)+pn(i  ,j))*tl_vbar(i  ,j,kstp)-        &
!^   &              (pn(i-1,j-1)+pn(i-1,j))*tl_vbar(i-1,j,kstp))+       &
!^   &             pnom_p(i,j)*                                         &
!^   &             ((pm(i-1,j  )+pm(i,j  ))*tl_ubar(i,j  ,kstp)-        &
!^   &              (pm(i-1,j-1)+pm(i,j-1))*tl_ubar(i,j-1,kstp))))
!^
          adfac=visc2_p(i,j)*0.5_r8*ad_cff
          adfac1=adfac*Drhs_p(i,j)
          adfac2=adfac1*pmon_p(i,j)
          adfac3=adfac1*pnom_p(i,j)
          ad_Drhs_p(i,j)=ad_Drhs_p(i,j)+                                &
     &                   (pmon_p(i,j)*                                  &
     &                    ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,kstp)-    &
     &                     (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,kstp))+   &
     &                    pnom_p(i,j)*                                  &
     &                    ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,kstp)-    &
     &                     (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,kstp)))*  &
     &                   adfac
          ad_vbar(i-1,j,kstp)=ad_vbar(i-1,j,kstp)-                      &
     &                        (pn(i-1,j-1)+pn(i-1,j))*adfac2
          ad_vbar(i  ,j,kstp)=ad_vbar(i  ,j,kstp)+                      &
     &                        (pn(i  ,j-1)+pn(i  ,j))*adfac2
          ad_ubar(i,j-1,kstp)=ad_ubar(i,j-1,kstp)-                      &
     &                        (pm(i-1,j-1)+pm(i,j-1))*adfac3
          ad_ubar(i,j  ,kstp)=ad_ubar(i,j  ,kstp)+                      &
     &                        (pm(i-1,j  )+pm(i,j  ))*adfac3
          ad_cff=0.0_r8
        END DO
      END DO
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!^        tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
!^        tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
!^
          ad_cff=ad_cff+                                                &
     &           om_r(i,j)*om_r(i,j)*ad_VFe(i,j)+                       &
     &           on_r(i,j)*on_r(i,j)*ad_UFx(i,j)
          ad_VFe(i,j)=0.0_r8
          ad_UFx(i,j)=0.0_r8
!^        tl_cff=visc2_r(i,j)*0.5_r8*                                   &
!^   &           (tl_Drhs(i,j)*                                         &
!^   &            (pmon_r(i,j)*                                         &
!^   &             ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,kstp)-             &
!^   &              (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,kstp))-            &
!^   &             pnom_r(i,j)*                                         &
!^   &             ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,kstp)-             &
!^   &              (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,kstp)))+           &
!^   &            Drhs(i,j)*                                            &
!^   &            (pmon_r(i,j)*                                         &
!^   &             ((pn(i  ,j)+pn(i+1,j))*tl_ubar(i+1,j,kstp)-          &
!^   &              (pn(i-1,j)+pn(i  ,j))*tl_ubar(i  ,j,kstp))-         &
!^   &             pnom_r(i,j)*                                         &
!^   &             ((pm(i,j  )+pm(i,j+1))*tl_vbar(i,j+1,kstp)-          &
!^   &              (pm(i,j-1)+pm(i,j  ))*tl_vbar(i,j  ,kstp))))
!^
          adfac=visc2_r(i,j)*0.5_r8*ad_cff
          adfac1=adfac*Drhs(i,j)
          adfac2=adfac1*pmon_r(i,j)
          adfac3=adfac1*pnom_r(i,j)
          ad_Drhs(i,j)=ad_Drhs(i,j)+                                    &
     &                 (pmon_r(i,j)*                                    &
     &                  ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,kstp)-        &
     &                   (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,kstp))-       &
     &                  pnom_r(i,j)*                                    &
     &                  ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,kstp)-        &
     &                   (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,kstp)))*      &
     &                 adfac
          ad_ubar(i  ,j,kstp)=ad_ubar(i  ,j,kstp)-                      &
     &                        (pn(i-1,j)+pn(i  ,j))*adfac2
          ad_ubar(i+1,j,kstp)=ad_ubar(i+1,j,kstp)+                      &
     &                        (pn(i  ,j)+pn(i+1,j))*adfac2
          ad_vbar(i,j  ,kstp)=ad_vbar(i,j  ,kstp)+                      &
     &                        (pm(i,j-1)+pm(i,j  ))*adfac3
          ad_vbar(i,j+1,kstp)=ad_vbar(i,j+1,kstp)-                      &
     &                        (pm(i,j  )+pm(i,j+1))*adfac3
          ad_cff=0.0_r8
        END DO
      END DO
!
!  Compute total depth at PSI-points.
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
!^        tl_Drhs_p(i,j)=0.25_r8*(tl_Drhs(i,j  )+tl_Drhs(i-1,j  )+      &
!^   &                            tl_Drhs(i,j-1)+tl_Drhs(i-1,j-1))
!^
          adfac=0.25_r8*ad_Drhs_p(i,j)
          ad_Drhs(i-1,j-1)=ad_Drhs(i-1,j-1)+adfac
          ad_Drhs(i-1,j  )=ad_Drhs(i-1,j  )+adfac
          ad_Drhs(i,  j-1)=ad_Drhs(i  ,j-1)+adfac
          ad_Drhs(i  ,j  )=ad_Drhs(i  ,j  )+adfac
          ad_Drhs_p(i,j)=0.0_r8
        END DO
      END DO
#endif

#if (defined CURVGRID && defined UV_ADV) && !defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint of add in curvilinear transformation terms.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (j.ge.JstrV) THEN
# if defined DIAGNOSTICS_UV
!!          fac2=0.5_r8*(Vwrk(i,j)+Vwrk(i,j-1))
!!          DiaV2rhs(i,j,M2xadv)=DiaV2rhs(i,j,M2xadv)-fac1+fac2
!!          DiaV2rhs(i,j,M2yadv)=DiaV2rhs(i,j,M2yadv)-fac2
!!          DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)-fac1
# endif
!^          tl_rvbar(i,j)=tl_rvbar(i,j)-tl_fac1
!^
            ad_fac1=ad_fac1-ad_rvbar(i,j)
!^          tl_fac1=0.5_r8*(tl_VFe(i,j)+tl_VFe(i,j-1))
!^
            adfac=0.5_r8*ad_fac1
            ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac
            ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
            ad_fac1=0.0_r8
          END IF
!
          IF (i.ge.IstrU) THEN
# if defined DIAGNOSTICS_UV
!!          fac2=0.5_r8*(Uwrk(i,j)+Uwrk(i-1,j))
!!          DiaU2rhs(i,j,M2xadv)=DiaU2rhs(i,j,M2xadv)+fac1-fac2
!!          DiaU2rhs(i,j,M2yadv)=DiaU2rhs(i,j,M2yadv)+fac2
!!          DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+fac1
# endif
!^          tl_rubar(i,j)=tl_rubar(i,j)+tl_fac1
!^
            ad_fac1=ad_fac1+ad_rubar(i,j)
!^          tl_fac1=0.5_r8*(tl_UFx(i,j)+tl_UFx(i-1,j))
!^
            adfac=0.5_r8*ad_fac1
            ad_UFx(i-1,j)=ad_UFx(i-1,j)+adfac
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
            ad_fac1=0.0_r8
          END IF
        END DO
      END DO
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=0.5_r8*(vrhs(i,j  )+                                     &
     &                 vrhs(i,j+1))
          cff2=0.5_r8*(urhs(i  ,j)+
     &                 urhs(i+1,j))
          cff3=cff1*dndx(i,j)
          cff4=cff2*dmde(i,j)
          cff=Drhs(i,j)*(cff3-cff4)
# if defined DIAGNOSTICS_UV
!!        cff=Drhs(i,j)*cff4
!!        Uwrk(i,j)=-cff*cff1                  ! ubar equation, ETA-term
!!        Vwrk(i,j)=-cff*cff2                  ! vbar equation, ETA-term
# endif
!^        tl_VFe(i,j)=tl_cff*cff2+cff*tl_cff2
!^        tl_UFx(i,j)=tl_cff*cff1+cff*tl_cff1
!^
          ad_cff=ad_cff+                                                &
     &           cff1*ad_UFx(i,j)+                                      &
     &           cff2*ad_VFe(i,j)
          ad_cff1=ad_cff1+cff*ad_UFx(i,j)
          ad_cff2=ad_cff2+cff*ad_VFe(i,j)
          ad_UFx(i,j)=0.0_r8
          ad_VFe(i,j)=0.0_r8
!^        tl_cff=tl_Drhs(i,j)*(cff3-cff4)+                              &
!^   &           Drhs(i,j)*(tl_cff3-tl_cff4)
!^
          adfac=Drhs(i,j)*ad_cff
          ad_cff4=ad_cff4-adfac
          ad_cff3=ad_cff3+adfac
          ad_Drhs(i,j)=ad_Drhs(i,j)+(cff3-cff4)*ad_cff
          ad_cff=0.0_r8
!^        tl_cff4=tl_cff2*dmde(i,j)
!^
          ad_cff2=ad_cff2+dmde(i,j)*ad_cff4
          ad_cff4=0.0_r8
!^        tl_cff3=tl_cff1*dndx(i,j)
!^
          ad_cff1=ad_cff1+dndx(i,j)*ad_cff3
          ad_cff3=0.0_r8
!^        tl_cff2=0.5_r8*(tl_urhs(i  ,j)+                               &
!^   &                    tl_urhs(i+1,j))
!^
          adfac=0.5_r8*ad_cff2
          ad_urhs(i  ,j)=ad_urhs(i  ,j)+adfac
          ad_urhs(i+1,j)=ad_urhs(i+1,j)+adfac
          ad_cff2=0.0_r8
!^        tl_cff1=0.5_r8*(tl_vrhs(i,j  )+                               &
!^   &                    tl_vrhs(i,j+1))
!^
          adfac=0.5_r8*ad_cff1
          ad_vrhs(i,j  )=ad_vrhs(i,j  )+adfac
          ad_vrhs(i,j+1)=ad_vrhs(i,j+1)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
#endif

#if (defined UV_COR && !defined SOLVE3D) || defined STEP2D_CORIOLIS
!
!-----------------------------------------------------------------------
!  Adjoint of add in Coriolis term.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (j.ge.JstrV) THEN
# if defined DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2fcor)=-fac2
# endif
!^          tl_rvbar(i,j)=tl_rvbar(i,j)-tl_fac2
!^
            ad_fac2=ad_fac2-ad_rvbar(i,j)
!^          tl_fac2=0.5_r8*(tl_VFe(i,j)+tl_VFe(i,j-1))
!^
            adfac=0.5_r8*ad_fac2
            ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac
            ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
            ad_fac2=0.0_r8
          END IF
!
          IF (i.ge.IstrU) THEN
# if defined DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2fcor)=fac1
# endif
!^          tl_rubar(i,j)=tl_rubar(i,j)+tl_fac1
!^
            ad_fac1=ad_fac1+ad_rubar(i,j)
!^          tl_fac1=0.5_r8*(tl_UFx(i,j)+tl_UFx(i-1,j))
!^
            adfac=0.5_r8*ad_fac1
            ad_UFx(i-1,j)=ad_UFx(i-1,j)+adfac
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
            ad_fac1=0.0_r8
          END IF
        END DO
      END DO
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff=0.5_r8*Drhs(i,j)*fomn(i,j)
!^        tl_VFe(i,j)=tl_cff*(urhs(i  ,j)+                              &
!^   &                        urhs(i+1,j))+                             &
!^   &                cff*(tl_urhs(i  ,j)+                              &
!^   &                     tl_urhs(i+1,j))
!^
          adfac=cff*ad_VFe(i,j)
          ad_urhs(i  ,j)=ad_urhs(i  ,j)+adfac
          ad_urhs(i+1,j)=ad_urhs(i+1,j)+adfac
          ad_cff=ad_cff+                                                &
     &           (urhs(i  ,j)+                                          &
     &            urhs(i+1,j))*ad_VFe(i,j)
          ad_VFe(i,j)=0.0_r8
!
!^        tl_UFx(i,j)=tl_cff*(vrhs(i,j  )+                              &
!^   &                        vrhs(i,j+1))+                             &
!^   &                cff*(tl_vrhs(i,j  )+                              &
!^   &                     tl_vrhs(i,j+1))
!^
          adfac=cff*ad_UFx(i,j)
          ad_vrhs(i,j  )=ad_vrhs(i,j  )+adfac
          ad_vrhs(i,j+1)=ad_vrhs(i,j+1)+adfac
          ad_cff=ad_cff+                                                &
     &           (vrhs(i,j  )+                                          &
     &            vrhs(i,j+1))*ad_UFx(i,j)
          ad_UFx(i,j)=0.0_r8
!^        tl_cff=0.5_r8*tl_Drhs(i,j)*fomn(i,j)
!^
          ad_Drhs(i,j)=ad_Drhs(i,j)+0.5_r8*fomn(i,j)*ad_cff
          ad_cff=0.0_r8
        END DO
      END DO
#endif

#if defined UV_ADV && !defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint of add in horizontal advection of momentum.
!-----------------------------------------------------------------------
!
!  Add advection to RHS terms.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (j.ge.JstrV) THEN
# if defined DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2xadv)=-cff3
!!          DiaV2rhs(i,j,M2yadv)=-cff4
!!          DiaV2rhs(i,j,M2hadv)=-fac2
# endif
!^          tl_rvbar(i,j)=tl_rvbar(i,j)-tl_fac2
!^
            ad_fac2=ad_fac2-ad_rvbar(i,j)
!^          tl_fac2=tl_cff3+tl_cff4
!^
            ad_cff3=ad_cff3+ad_fac2
            ad_cff4=ad_cff4+ad_fac2
            ad_fac2=0.0_r8
!^          tl_cff4=tl_VFe(i,j)-tl_VFe(i,j-1)
!^
            ad_VFe(i,j-1)=ad_VFe(i,j-1)-ad_cff4
            ad_VFe(i,j  )=ad_VFe(i,j  )+ad_cff4
            ad_cff4=0.0_r8
!^          tl_cff3=tl_VFx(i+1,j)-tl_VFx(i,j)
!^
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-ad_cff3
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+ad_cff3
            ad_cff3=0.0_r8
          END IF
!
          IF (i.ge.IstrU) THEN
# if defined DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2xadv)=-cff1
!!          DiaU2rhs(i,j,M2yadv)=-cff2
!!          DiaU2rhs(i,j,M2hadv)=-fac1
# endif
!^          tl_rubar(i,j)=tl_rubar(i,j)-tl_fac1
!^
            ad_fac1=ad_fac1-ad_rubar(i,j)
!^          tl_fac1=tl_cff1+tl_cff2
!^
            ad_cff1=ad_cff1+ad_fac1
            ad_cff2=ad_cff2+ad_fac1
            ad_fac1=0.0_r8
!^          tl_cff2=tl_UFe(i,j+1)-tl_UFe(i,j)
!^
            ad_UFe(i,j  )=ad_UFe(i,j  )-ad_cff2
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+ad_cff2
            ad_cff2=0.0_r8
!^          tl_cff1=tl_UFx(i,j)-tl_UFx(i-1,j)
!^
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-ad_cff1
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+ad_cff1
            ad_cff1=0.0_r8
          END IF
        END DO
      END DO

# ifdef UV_C2ADVECTION
!
!  Second-order, centered differences advection fluxes.
!
      DO j=Jstr-1,Jend
        DO i=Istr,Iend
!^        tl_UFe(i,j+1)=0.25_r8*                                        &
#  ifdef MASKING
!^   &                  pmask(i,j+1)*                                   &
#  endif
!^   &                  ((tl_DVom(i,j+1)+tl_DVom(i-1,j+1))*             &
!^   &                   (urhs(i,j+1)+                                  &
!^   &                    urhs(i,j  ))+                                 &
!^   &                   (DVom(i,j+1)+DVom(i-1,j+1))*                   &
!^   &                   (tl_urhs(i,j+1)+                               &
!^   &                    tl_urhs(i,j  )))
!^
          adfac=0.25_r8*ad_UFe(i,j+1)
          adfac1=adfac*(urhs(i,j+1)+                                    &
     &                  urhs(i,j  ))
          adfac2=adfac*(DVom(i,j+1)+DVom(i-1,j+1))
          ad_DVom(i-1,j+1)=ad_DVom(i-1,j+1)+adfac1
          ad_DVom(i  ,j+1)=ad_DVom(i,j+1)+adfac1
          ad_urhs(i,j  )=ad_urhs(i,j  )+adfac2
          ad_urhs(i,j+1)=ad_urhs(i,j+1)+adfac2
          ad_UFe(i,j+1)=0.0_r8
!
          IF (j.ge.JstrV-1) THEN
!^          tl_VFe(i,j)=0.25_r8*                                        &
!^   &                  ((tl_DVom(i,j)+tl_DVom(i,j+1))*                 &
!^   &                   (vrhs(i,j  )+                                  &
!^   &                    vrhs(i,j+1))+                                 &
!^   &                   (DVom(i,j)+DVom(i,j+1))*                       &
!^   &                   (tl_vrhs(i,j  )+                               &
!^   &                    tl_vrhs(i,j+1)))
!^
            adfac=0.25_r8*ad_VFe(i,j)
            adfac1=adfac*(vrhs(i,j  )+                                  &
     &                    vrhs(i,j+1))
            adfac2=adfac*(DVom(i,j)+DVom(i,j+1))
            ad_DVom(i,j  )=ad_DVom(i,j  )+adfac1
            ad_DVom(i,j+1)=ad_DVom(i,j+1)+adfac1
            ad_vrhs(i,j  )=ad_vrhs(i,j  )+adfac2
            ad_vrhs(i,j+1)=ad_vrhs(i,j+1)+adfac2
            ad_VFe(i,j)=0.0_r8
          END IF
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend
!^        tl_VFx(i+1,j)=0.25_r8*                                        &
#  ifdef MASKING
!^   &                  pmask(i+1,j)*                                   &
#  endif
!^   &                  ((tl_DUon(i+1,j)+tl_DUon(i+1,j-1))*             &
!^   &                   (vrhs(i+1,j)+                                  &
!^   &                    vrhs(i  ,j))+                                 &
!^   &                    (DUon(i+1,j)+DUon(i+1,j-1))*                  &
!^   &                    (tl_vrhs(i+1,j)+                              &
!^   &                     tl_vrhs(i  ,j)))
!^
          adfac=0.25_r8*                                                &
#  ifdef MASKING
     &          pmask(i+1,j)*                                           &
#  endif
     &          ad_VFx(i+1,j)
          adfac1=adfac*(vrhs(i+1,j)+                                    &
     &                  vrhs(i  ,j))
          adfac2=adfac*(DUon(i+1,j)+DUon(i+1,j-1))
          ad_DUon(i+1,j-1)=ad_DUon(i+1,j-1)+adfac1
          ad_DUon(i+1,j  )=ad_DUon(i+1,j  )+adfac1
          ad_vrhs(i  ,j)=ad_vrhs(i  ,j)+adfac2
          ad_vrhs(i+1,j)=ad_vrhs(i+1,j)+adfac2
          ad_VFx(i+1,j)=0.0_r8
!
          IF (i.ge.IstrU-1) THEN
!^          tl_UFx(i,j)=0.25_r8*                                        &
!^   &                  ((tl_DUon(i,j)+tl_DUon(i+1,j))*                 &
!^   &                   (urhs(i  ,j)+                                  &
!^   &                    urhs(i+1,j))+                                 &
!^   &                   (DUon(i,j)+DUon(i+1,j))*                       &
!^   &                   (tl_urhs(i  ,j)+                               &
!^   &                    tl_urhs(i+1,j)))
!^
            adfac=0.25_r8*ad_UFx(i,j
            adfac1=adfac*(urhs(i  ,j)+                                  &
     &                    urhs(i+1,j))
            adfac2=adfac*(DUon(i,j)+DUon(i+1,j))
            ad_DUon(i  ,j)=ad_DUon(i  ,j)+adfac1
            ad_DUon(i+1,j)=ad_DUon(i+1,j)+adfac1
            ad_urhs(i  ,j)=ad_urhs(i  ,j)+adfac2
            ad_urhs(i+1,j)=ad_urhs(i+1,j)+adfac2
            ad_UFx(i,j)=0.0_r8
          END IF
        END DO
      END DO
!

# elif defined UV_C4ADVECTION
!
!  Fourth-order, centered differences v-momentum advection fluxes.
!
      DO j=JstrVm1,Jendp1                               ! BASIC STATE
        DO i=Istr,Iend
          grad (i,j)=vrhs(i,j-1)-2.0_r8*vrhs(i,j)+                      &
     &               vrhs(i,j+1)
          Dgrad(i,j)=DVom(i,j-1)-2.0_r8*DVom(i,j)+DVom(i,j+1)
        END DO
      END DO
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            grad (i,Jend+1)=grad (i,Jend)
            Dgrad(i,Jend+1)=Dgrad(i,Jend)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            grad (i,Jstr)=grad (i,Jstr+1)
            Dgrad(i,Jstr)=Dgrad(i,Jstr+1)
          END DO
        END IF
      END IF
!                                                         d/dy(Dvv/m)
      cff=1.0_r8/6.0_r8
      DO j=JstrV-1,Jend
        DO i=Istr,Iend
!^        tl_VFe(i,j)=0.25_r8*                                          &
!^   &                ((tl_vrhs(i,j  )+                                 &
!^   &                  tl_vrhs(i,j+1)-                                 &
!^   &                  cff*(tl_grad (i,j)+tl_grad (i,j+1)))*           &
!^   &                 (DVom(i,j)+DVom(i,j+1)-                          &
!^   &                  cff*(Dgrad(i,j)+Dgrad(i,j+1)))+                 &
!^   &                 (vrhs(i,j  )+                                    &
!^   &                  vrhs(i,j+1)-                                    &
!^   &                  cff*(grad (i,j)+grad (i,j+1)))*                 &
!^   &                 (tl_DVom(i,j)+tl_DVom(i,j+1)-                    &
!^   &                  cff*(tl_Dgrad(i,j)+tl_Dgrad(i,j+1))))
!^
          adfac=0.25_r8*ad_VFe(i,j)
          adfac1=adfac*(DVom(i,j)+DVom(i,j+1)-                          &
     &                  cff*(Dgrad(i,j)+Dgrad(i,j+1)))
          adfac2=adfac1*cff
          adfac3=adfac*(vrhs(i,j  )+                                    &
     &                  vrhs(i,j+1,krhs)-                               &
     &                  cff*(grad (i,j)+grad (i,j+1)))
          adfac4=adfac3*cff
          ad_vrhs(i,j  )=ad_vrhs(i,j  )+adfac1
          ad_vrhs(i,j+1)=ad_vrhs(i,j+1)+adfac1
          ad_grad (i,j  )=ad_grad (i,j  )-adfac2
          ad_grad (i,j+1)=ad_grad (i,j+1)-adfac2
          ad_DVom(i,j  )=ad_DVom(i,j  )+adfac3
          ad_DVom(i,j+1)=ad_DVom(i,j+1)+adfac3
          ad_Dgrad(i,j  )=ad_Dgrad(i,j  )-adfac4
          ad_Dgrad(i,j+1)=ad_Dgrad(i,j+1)-adfac4
          ad_VFe(i,j)=0.0_r8
        END DO
      END DO
!
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
!^          tl_Dgrad(i,Jend+1)=tl_Dgrad(i,Jend)
!^
            ad_Dgrad(i,Jend)=ad_Dgrad(i,Jend)+ad_Dgrad(i,Jend+1)
            ad_Dgrad(i,Jend+1)=0.0_r8
!^          tl_grad (i,Jend+1)=tl_grad (i,Jend)
!^
            ad_grad (i,Jend)=ad_grad (i,Jend)+ad_grad (i,Jend+1)
            ad_grad (i,Jend+1)=0.0_r8
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
!^          tl_Dgrad(i,Jstr)=tl_Dgrad(i,Jstr+1)
!^
            ad_Dgrad(i,Jstr+1)=ad_Dgrad(i,Jstr+1)+ad_Dgrad(i,Jstr)
            ad_Dgrad(i,Jstr)=0.0_r8
!^          tl_grad (i,Jstr)=tl_grad (i,Jstr+1)
!^
            ad_grad (i,Jstr+1)=ad_grad (i,Jstr+1)+ad_grad (i,Jstr)
            ad_grad (i,Jstr)=0.0_r8
          END DO
       END IF
      END IF
      DO j=JstrVm1,Jendp1
        DO i=Istr,Iend
!^        tl_Dgrad(i,j)=tl_DVom(i,j-1)-2.0_r8*tl_DVom(i,j)+             &
!^   &                  tl_DVom(i,j+1)
!^
          ad_DVom(i,j-1)=ad_DVom(i,j-1)+ad_Dgrad(i,j)
          ad_DVom(i,j  )=ad_DVom(i,j  )-2.0_r8*ad_Dgrad(i,j)
          ad_DVom(i,j+1)=ad_DVom(i,j+1)+ad_Dgrad(i,j)
          ad_Dgrad(i,j)=0.0_r8
!^        tl_grad(i,j)=tl_vrhs(i,j-1)-2.0_r8*tl_vrhs(i,j)+              &
!^   &                 tl_vrhs(i,j+1)
!^
          ad_vrhs(i,j-1)=ad_vrhs(i,j-1)+ad_grad(i,j)
          ad_vrhs(i,j  )=ad_vrhs(i,j  )-2.0_r8*ad_grad(i,j)
          ad_vrhs(i,j+1)=ad_vrhs(i,j+1)+ad_grad(i,j)
          ad_grad(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=JstrV,Jend                                   ! BASIC STATE
        DO i=Istrm1,Iendp1
          grad(i,j)=vrhs(i-1,j)-2.0_r8*vrhs(i,j,krhs)+                  &
     &              vrhs(i+1,j)
        END DO
      END DO
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
            grad(Istr-1,j)=grad(Istr,j)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            grad(Iend+1,j)=grad(Iend,j)
          END DO
        END IF
      END IF
      DO j=JstrV-1,Jend
        DO i=Istr,Iend+1
          Dgrad(i,j)=DUon(i,j-1)-2.0_r8*DUon(i,j)+DUon(i,j+1)
        END DO
      END DO
!                                                         d/dx(Duv/n)
      cff=1.0_r8/6.0_r8
      DO j=JstrV,Jend
        DO i=Istr,Iend+1
!^        tl_VFx(i,j)=0.25_r8*                                          &
!^   &                ((tl_vrhs(i  ,j)+                                 &
!^   &                  tl_vrhs(i-1,j)-                                 &
!^   &                  cff*(tl_grad (i,j)+tl_grad (i-1,j)))*           &
!^   &                 (DUon(i,j)+DUon(i,j-1)-                          &
!^   &                  cff*(Dgrad(i,j)+Dgrad(i,j-1)))+                 &
!^   &                 (vrhs(i  ,j)+                                    &
!^   &                  vrhs(i-1,j)-                                    &
!^   &                  cff*(grad (i,j)+grad (i-1,j)))*                 &
!^   &                 (tl_DUon(i,j)+tl_DUon(i,j-1)-                    &
!^   &                  cff*(tl_Dgrad(i,j)+tl_Dgrad(i,j-1))))
!^
          adfac=0.25_r8*ad_VFx(i,j)
          adfac1=adfac*(DUon(i,j)+DUon(i,j-1)-                          &
     &                  cff*(Dgrad(i,j)+Dgrad(i,j-1)))
          adfac2=adfac1*cff
          adfac3=adfac*(vrhs(i  ,j)+                                    &
     &                  vrhs(i-1,j)-                                    &
     &                  cff*(grad (i,j)+grad (i-1,j)))
          adfac4=adfac3*cff
          ad_vrhs(i-1,j)=ad_vrhs(i-1,j)+adfac1
          ad_vrhs(i  ,j)=ad_vrhs(i  ,j)+adfac1
          ad_grad (i-1,j)=ad_grad (i-1,j)-adfac2
          ad_grad (i  ,j)=ad_grad (i  ,j)-adfac2
          ad_DUon(i,j-1)=ad_DUon(i,j-1)+adfac3
          ad_DUon(i,j  )=ad_DUon(i,j  )+adfac3
          ad_Dgrad(i,j-1)=ad_Dgrad(i,j-1)-adfac4
          ad_Dgrad(i,j  )=ad_Dgrad(i,j  )-adfac4
          ad_VFx(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=JstrV-1,Jend
        DO i=Istr,Iend+1
!^        tl_Dgrad(i,j)=tl_DUon(i,j-1)-2.0_r8*tl_DUon(i,j)+             &
!^   &                  tl_DUon(i,j+1)
!^
          ad_DUon(i,j-1)=ad_DUon(i,j-1)+ad_Dgrad(i,j)
          ad_DUon(i,j  )=ad_DUon(i,j  )-2.0_r8*ad_Dgrad(i,j)
          ad_DUon(i,j+1)=ad_DUon(i,j+1)+ad_Dgrad(i,j)
          ad_Dgrad(i,j)=0.0_r8
        END DO
      END DO
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
!^          tl_grad(Iend+1,j)=tl_grad(Iend,j)
!^
            ad_grad(Iend,j)=ad_grad(Iend,j)+ad_grad(Iend+1,j)
            ad_grad(Iend+1,j)=0.0_r8
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
!^          tl_grad(Istr-1,j)=tl_grad(Istr,j)
!^
            ad_grad(Istr,j)=ad_grad(Istr,j)+ad_grad(Istr-1,j)
            ad_grad(Istr-1,j)=0.0_r8
          END DO
        END IF
      END IF
      DO j=JstrV,Jend
        DO i=Istrm1,Iendp1
!^        tl_grad(i,j)=tl_vrhs(i-1,j)-2.0_r8*tl_vrhs(i,j)+              &
!^   &                 tl_vrhs(i+1,j)
!^
          ad_vrhs(i-1,j)=ad_vrhs(i-1,j)+ad_grad(i,j)
          ad_vrhs(i  ,j)=ad_vrhs(i  ,j)-2.0_r8*ad_grad(i,j)
          ad_vrhs(i+1,j)=ad_vrhs(i+1,j)+ad_grad(i,j)
          ad_grad(i,j)=0.0_r8
        END DO
      END DO
!
!  Fourth-order, centered differences u-momentum advection fluxes.
!
      DO j=Jstrm1,Jendp1                                ! BASIC STATE
        DO i=IstrU,Iend
          grad(i,j)=urhs(i,j-1)-2.0_r8*urhs(i,j)+                       &
     &              urhs(i,j+1)
        END DO
      END DO
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            grad(i,Jstr-1)=grad(i,Jstr)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
            grad(i,Jend+1)=grad(i,Jend)
          END DO
        END IF
      END IF
      DO j=Jstr,Jend+1
        DO i=IstrU-1,Iend
          Dgrad(i,j)=DVom(i-1,j)-2.0_r8*DVom(i,j)+DVom(i+1,j)
        END DO
      END DO
!                                                         d/dy(Duv/m)
      cff=1.0_r8/6.0_r8
      DO j=Jstr,Jend+1
        DO i=IstrU,Iend
!^        tl_UFe(i,j)=0.25_r8*                                          &
!^   &                ((tl_urhs(i,j  )+                                 &
!^   &                  tl_urhs(i,j-1)-                                 &
!^   &                  cff*(tl_grad (i,j)+tl_grad (i,j-1)))*           &
!^   &                 (DVom(i,j)+DVom(i-1,j)-                          &
!^   &                  cff*(Dgrad(i,j)+Dgrad(i-1,j)))+                 &
!^   &                 (urhs(i,j  )+                                    &
!^   &                  urhs(i,j-1)-                                    &
!^   &                  cff*(grad (i,j)+grad (i,j-1)))*                 &
!^   &                 (tl_DVom(i,j)+tl_DVom(i-1,j)-                    &
!^   &                  cff*(tl_Dgrad(i,j)+tl_Dgrad(i-1,j))))
!^
          adfac=0.25_r8*ad_UFe(i,j)
          adfac1=adfac*(DVom(i,j)+DVom(i-1,j)-                          &
     &                  cff*(Dgrad(i,j)+Dgrad(i-1,j)))
          adfac2=adfac1*cff
          adfac3=adfac*(urhs(i,j  )+                                    &
     &                  urhs(i,j-1,krhs)-                               &
     &                  cff*(grad (i,j)+grad (i,j-1)))
          adfac4=adfac3*cff
          ad_urhs(i,j-1)=ad_urhs(i,j-1)+adfac1
          ad_urhs(i,j  )=ad_urhs(i,j  )+adfac1
          ad_grad (i,j-1)=ad_grad (i,j-1)-adfac2
          ad_grad (i,j  )=ad_grad (i,j  )-adfac2
          ad_DVom(i-1,j)=ad_DVom(i-1,j)+adfac3
          ad_DVom(i  ,j)=ad_DVom(i  ,j)+adfac3
          ad_Dgrad(i-1,j)=ad_Dgrad(i-1,j)-adfac4
          ad_Dgrad(i  ,j)=ad_Dgrad(i  ,j)-adfac4
          ad_UFe(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=Jstr,Jend+1
        DO i=IstrU-1,Iend
!^        tl_Dgrad(i,j)=tl_DVom(i-1,j)-2.0_r8*tl_DVom(i,j)+             &
!^   &                  tl_DVom(i+1,j)
!^
          ad_DVom(i-1,j)=ad_DVom(i-1,j)+ad_Dgrad(i,j)
          ad_DVom(i  ,j)=ad_DVom(i  ,j)-2.0_r8*ad_Dgrad(i,j)
          ad_DVom(i+1,j)=ad_DVom(i+1,j)+ad_Dgrad(i,j)
          ad_Dgrad(i,j)=0.0_r8
        END DO
      END DO
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
!^          tl_grad(i,Jend+1)=tl_grad(i,Jend)
!^
            ad_grad(i,Jend)=ad_grad(i,Jend)+ad_grad(i,Jend+1)
            ad_grad(i,Jend+1)=0.0_r8
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
!^          tl_grad(i,Jstr-1)=tl_grad(i,Jstr)
!^
            ad_grad(i,Jstr)=ad_grad(i,Jstr)+ad_grad(i,Jstr-1)
            ad_grad(i,Jstr-1)=0.0_r8
          END DO
        END IF
      END IF
      DO j=Jstrm1,Jendp1
        DO i=IstrU,Iend
!^        tl_grad(i,j)=tl_urhs(i,j-1)-2.0_r8*tl_urhs(i,j)+              &
!^   &                 tl_urhs(i,j+1)
!^
          ad_urhs(i,j-1)=ad_urhs(i,j-1)+ad_grad(i,j)
          ad_urhs(i,j  )=ad_urhs(i,j  )-2.0_r8*ad_grad(i,j)
          ad_urhs(i,j+1)=ad_urhs(i,j+1)+ad_grad(i,j)
          ad_grad(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=Jstr,Jend                                    ! BASIC STATE
        DO i=IstrUm1,Iendp1
          grad (i,j)=urhs(i-1,j)-2.0_r8*urhs(i,j)+                      &
     &               urhs(i+1,j,krhs)
          Dgrad(i,j)=DUon(i-1,j)-2.0_r8*DUon(i,j)+DUon(i+1,j)
        END DO
      END DO
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            grad (Istr,j)=grad (Istr+1,j)
            Dgrad(Istr,j)=Dgrad(Istr+1,j)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            grad (Iend+1,j)=grad (Iend,j)
            Dgrad(Iend+1,j)=Dgrad(Iend,j)
          END DO
        END IF
      END IF
!                                                         d/dx(Duu/n)
      cff=1.0_r8/6.0_r8
      DO j=Jstr,Jend
        DO i=IstrU-1,Iend
!^        tl_UFx(i,j)=0.25_r8*                                          &
!^   &                ((urhs(i  ,j)+                                    &
!^   &                  urhs(i+1,j)-                                    &
!^   &                  cff*(grad (i,j)+grad (i+1,j)))*                 &
!^   &                 (tl_DUon(i,j)+tl_DUon(i+1,j)-                    &
!^   &                  cff*(tl_Dgrad(i,j)+tl_Dgrad(i+1,j)))+           &
!^   &                 (tl_urhs(i  ,j)+                                 &
!^   &                  tl_urhs(i+1,j)-                                 &
!^   &                  cff*(tl_grad (i,j)+tl_grad (i+1,j)))*           &
!^   &                 (DUon(i,j)+DUon(i+1,j)-                          &
!^   &                  cff*(Dgrad(i,j)+Dgrad(i+1,j))))
!^
          adfac=0.25_r8*ad_UFx(i,j)
          adfac1=adfac*(DUon(i,j)+DUon(i+1,j)-                          &
     &                  cff*(Dgrad(i,j)+Dgrad(i+1,j)))
          adfac2=adfac1*cff
          adfac3=adfac*(urhs(i  ,j)+                                    &
     &                  urhs(i+1,j)-                                    &
     &                  cff*(grad (i,j)+grad (i+1,j)))
          adfac4=adfac3*cff
          ad_urhs(i  ,j)=ad_urhs(i  ,j)+adfac1
          ad_urhs(i+1,j)=ad_urhs(i+1,j)+adfac1
          ad_grad (i  ,j)=ad_grad (i  ,j)-adfac2
          ad_grad (i+1,j)=ad_grad (i+1,j)-adfac2
          ad_DUon(i  ,j)=ad_DUon(i  ,j)+adfac3
          ad_DUon(i+1,j)=ad_DUon(i+1,j)+adfac3
          ad_Dgrad(i  ,j)=ad_Dgrad(i  ,j)-adfac4
          ad_Dgrad(i+1,j)=ad_Dgrad(i+1,j)-adfac4
          ad_UFx(i,j)=0.0_r8
        END DO
      END DO
!
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
!^          tl_Dgrad(Iend+1,j)=tl_Dgrad(Iend,j)
!^
            ad_Dgrad(Iend,j)=ad_Dgrad(Iend,j)+ad_Dgrad(Iend+1,j)
            ad_Dgrad(Iend+1,j)=0.0_r8
!^          tl_grad (Iend+1,j)=tl_grad (Iend,j)
!^
            ad_grad (Iend,j)=ad_grad (Iend,j)+ad_grad (Iend+1,j)
            ad_grad (Iend+1,j)=0.0_r8
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
!^          tl_Dgrad(Istr,j)=tl_Dgrad(Istr+1,j)
!^
            ad_Dgrad(Istr+1,j)=ad_Dgrad(Istr+1,j)+ad_Dgrad(Istr,j)
            ad_Dgrad(Istr,j)=0.0_r8
!^          tl_grad (Istr,j)=tl_grad (Istr+1,j)
!^
            ad_grad (Istr+1,j)=ad_grad (Istr+1,j)+ad_grad (Istr,j)
            ad_grad (Istr,j)=0.0_r8
          END DO
        END IF
      END IF
      DO j=Jstr,Jend
        DO i=IstrUm1,Iendp1
!^        tl_Dgrad(i,j)=tl_DUon(i-1,j)-2.0_r8*tl_DUon(i,j)+             &
!^   &                  tl_DUon(i+1,j)
!^
          ad_DUon(i-1,j)=ad_DUon(i-1,j)+ad_Dgrad(i,j)
          ad_DUon(i  ,j)=ad_DUon(i  ,j)-2.0_r8*ad_Dgrad(i,j)
          ad_DUon(i+1,j)=ad_DUon(i+1,j)+ad_Dgrad(i,j)
          ad_Dgrad(i,j)=0.0_r8
!^        tl_grad(i,j)=tl_urhs(i-1,j)-2.0_r8*tl_urhs(i,j)+              &
!^   &                 tl_urhs(i+1,j)
!^
          ad_urhs(i-1,j)=ad_urhs(i-1,j)+ad_grad (i,j)
          ad_urhs(i  ,j)=ad_urhs(i  ,j)-2.0_r8*ad_grad (i,j)
          ad_urhs(i+1,j)=ad_urhs(i+1,j)+ad_grad (i,j)
          ad_grad(i,j)=0.0_r8
        END DO
      END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  Adjoint of compute pressure-gradient terms.
!-----------------------------------------------------------------------
!
      cff1=0.5_r8*g
#if defined VAR_RHO_2D && defined SOLVE3D
      cff2=0.333333333333_r8
#endif
#if defined ATM_PRESS && !defined SOLVE3D
      cff3=0.5_r8*100.0_r8/rho0
#endif
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (j.ge.JstrV) THEN
#ifdef DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2pgrd)=rvbar(i,j)
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
!^          tl_rvbar(i,j)=tl_rvbar(i,j)-                                &
!^   &                    cff1*om_v(i,j)*                               &
!^   &                    ((tl_h(i,j-1)+tl_h(i,j)+                      &
!^   &                      tl_rzeta(i,j-1)+tl_rzeta(i,j))*             &
!^   &                     (eq_tide(i,j)-eq_tide(i,j-1))+               &
!^   &                     (h(i,j-1)+h(i,j)+                            &
!^   &                      rzeta(i,j-1)+rzeta(i,j))*                   &
!^   &                     (tl_eq_tide(i,j)-tl_eq_tide(i,j-1)))
!^
            adfac=cff1*om_v(i,j)*ad_rvbar(i,j)
            adfac1=adfac*(eq_tide(i,j)-eq_tide(i,j-1))
            adfac2=adfac*(h(i,j-1)+h(i,j)+                            &
     &                    rzeta(i,j-1)+rzeta(i,j))
            ad_h(i,j-1)=ad_h(i,j-1)-adfac1
            ad_h(i,j  )=ad_h(i,j  )-adfac1
            ad_rzeta(i,j-1)=ad_rzeta(i,j-1)-adfac1
            ad_rzeta(i,j  )=ad_rzeta(i,j  )-adfac1
            ad_eq_tide(i,j-1)=ad_eq_tide(i,j-1)+adfac2
            ad_eq_tide(i,j  )=ad_eq_tide(i,j  )-adfac2
#endif
#if defined ATM_PRESS && !defined SOLVE3D
!^          tl_rvbar(i,j)=tl_rvbar(i,j)-                                &
!^   &                    cff3*om_v(i,j)*                               &
!^   &                    (tl_h(i,j-1)+tl_h(i,j)+                       &
!^   &                     tl_rzeta(i,j-1)+tl_rzeta(i,j))*              &
!^   &                    (Pair(i,j)-Pair(i,j-1))
!^
            adfac=-cff3*om_v(i,j)*(Pair(i,j)-Pair(i,j-1)*ad_rvbar(i,j)
            ad_h(i,j-1)=ad_h(i,j-1)+adfac
            ad_h(i,j  )=ad_h(i,j  )+adfac
            ad_rzeta(i,j-1)=ad_rzeta(i,j-1)+adfac
            ad_rzeta(i,j  )=ad_rzeta(i,j  )+adfac
#endif
!^          tl_rvbar(i,j)=cff1*om_v(i,j)*                               &
!^   &                    ((tl_h(i,j-1)+                                &
!^   &                      tl_h(i,j  ))*                               &
!^   &                     (rzeta(i,j-1)-                               &
!^   &                      rzeta(i,j  ))+                              &
!^   &                     (h(i,j-1)+                                   &
!^   &                      h(i,j  ))*                                  &
!^   &                     (tl_rzeta(i,j-1)-                            &
!^   &                      tl_rzeta(i,j  ))+                           &
#if defined VAR_RHO_2D && defined SOLVE3D
!^   &                     (tl_h(i,j-1)-                                &
!^   &                      tl_h(i,j  ))*                               &
!^   &                     (rzetaSA(i,j-1)+                             &
!^   &                      rzetaSA(i,j  )+                             &
!^   &                      cff2*(rhoA(i,j-1)-                          &
!^   &                            rhoA(i,j  ))*                         &
!^   &                           (zwrk(i,j-1)-                          &
!^   &                            zwrk(i,j  )))+                        &
!^   &                     (h(i,j-1)-                                   &
!^   &                      h(i,j  ))*                                  &
!^   &                     (tl_rzetaSA(i,j-1)+                          &
!^   &                      tl_rzetaSA(i,j  )+                          &
!^   &                      cff2*((tl_rhoA(i,j-1)-                      &
!^   &                             tl_rhoA(i,j  ))*                     &
!^   &                            (zwrk(i,j-1)-                         &
!^   &                             zwrk(i,j  ))+                        &
!^   &                            (rhoA(i,j-1)-                         &
!^   &                             rhoA(i,j  ))*                        &
!^   &                            (tl_zwrk(i,j-1)-                      &
!^   &                             tl_zwrk(i,j  ))))+                   &
#endif
!^   &                     (tl_rzeta2(i,j-1)-                           &
!^   &                      tl_rzeta2(i,j  )))
!^
            adfac=cff1*om_v(i,j)*ad_rvbar(i,j)
            adfac1=adfac*(rzeta(i,j-1)-rzeta(i,j  ))
            adfac2=adfac*(h(i,j-1)+h(i,j  ))
            ad_h(i,j-1)=ad_h(i,j-1)+adfac1
            ad_h(i,j  )=ad_h(i,j  )+adfac1
            ad_rzeta(i,j-1)=ad_rzeta(i,j-1)+adfac2
            ad_rzeta(i,j  )=ad_rzeta(i,j  )-adfac2
            ad_rzeta2(i,j-1)=ad_rzeta2(i,j-1)+adfac
            ad_rzeta2(i,j  )=ad_rzeta2(i,j  )-adfac
#if defined VAR_RHO_2D && defined SOLVE3D
            adfac3=adfac*(rzetaSA(i,j-1)+                               &
     &                    rzetaSA(i,j  )+                               &
     &                    cff2*(rhoA(i,j-1)-                            &
     &                          rhoA(i,j  ))*                           &
     &                         (zwrk(i,j-1)-                            &
     &                          zwrk(i,j  )))
            adfac4=adfac2*cff2*(zwrk(i,j-1)-zwrk(i,j))
            adfac5=adfac2*cff2*(rhoA(i,j-1)-rhoA(i,j))
            ad_h(i,j-1)=ad_h(i,j-1)+adfac3
            ad_h(i,j  )=ad_h(i,j  )-adfac3
            ad_rzetaSA(i,j-1)=ad_rzetaSA(i,j-1)+adfac2
            ad_rzetaSA(i,j  )=ad_rzetaSA(i,j  )+adfac2
            ad_rhoA(i,j-1)=ad_rhoA(i,j-1)+adfac4
            ad_rhoA(i,j  )=ad_rhoA(i,j  )-adfac4
            ad_zwrk(i,j-1)=ad_zwrk(i,j-1)+adfac5
            ad_zwrk(i,j  )=ad_zwrk(i,j  )-adfac5
#endif
            ad_rvbar(i,j)=0.0_r8
          END IF
!
          IF (i.ge.IstrU) THEN
#ifdef DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2pgrd)=rubar(i,j)
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
!^          tl_rubar(i,j)=tl_rubar(i,j)-                                &
!^   &                    cff1*on_u(i,j)*                               &
!^   &                    ((tl_h(i-1,j)+tl_h(i,j)+                      &
!^   &                      tl_rzeta(i-1,j)+tl_rzeta(i,j))*             &
!^   &                     (eq_tide(i,j)-eq_tide(i-1,j))+               &
!^   &                     (h(i-1,j)+h(i,j)+                            &
!^   &                      rzeta(i-1,j)+rzeta(i,j))*                   &
!^   &                     (tl_eq_tide(i,j)-tl_eq_tide(i-1,j)))
!^
            adfac=cff1*on_u(i,j)*ad_rubar(i,j)
            adfac1=adfac*(eq_tide(i,j)-eq_tide(i-1,j))
            adfac2=adfac*(h(i-1,j)+h(i,j)+                              &
     &                    rzeta(i-1,j)+rzeta(i,j))
            ad_h(i-1,j)=ad_h(i-1,j)-adfac1
            ad_h(i  ,j)=ad_h(i  ,j)-adfac1
            ad_rzeta(i-1,j)=ad_rzeta(i-1,j)-adfac1
            ad_rzeta(i  ,j)=ad_rzeta(i  ,j)-adfac1
            ad_eq_tide(i-1,j)=ad_eq_tide(i-1,j)+adfac2
            ad_eq_tide(i  ,j)=ad_eq_tide(i  ,j)-adfac2
#endif
#if defined ATM_PRESS && !defined SOLVE3D
!^          tl_rubar(i,j)=tl_rubar(i,j)-                                &
!^   &                    cff3*on_u(i,j)*                               &
!^   &                    (tl_h(i-1,j)+tl_h(i,j)+                       &
!^   &                     tl_rzeta(i-1,j)+tl_rzeta(i,j))*              &
!^   &                    (Pair(i,j)-Pair(i-1,j))
!^
            adfac=-cff3*on_u(i,j)*(Pair(i,j)-Pair(i-1,j))*ad_rubar(i,j)
            ad_h(i-1,j)=ad_h(i-1,j)+adfac
            ad_h(i  ,j)=ad_h(i  ,j)+adfac
            ad_rzeta(i-1,j)=ad_rzeta(i-1,j)+adfac
            ad_rzeta(i  ,j)=ad_rzeta(i  ,j)+adfac
#endif
!^          tl_rubar(i,j)=cff1*on_u(i,j)*                               &
!^   &                    ((tl_h(i-1,j)+                                &
!^   &                      tl_h(i  ,j))*                               &
!^   &                     (rzeta(i-1,j)-                               &
!^   &                      rzeta(i  ,j))+                              &
!^   &                     (h(i-1,j)+                                   &
!^   &                      h(i  ,j))*                                  &
!^   &                     (tl_rzeta(i-1,j)-                            &
!^   &                      tl_rzeta(i  ,j))+                           &
#if defined VAR_RHO_2D && defined SOLVE3D
!^   &                     (tl_h(i-1,j)-                                &
!^   &                      tl_h(i  ,j))*                               &
!^   &                     (rzetaSA(i-1,j)+                             &
!^   &                      rzetaSA(i  ,j)+                             &
!^   &                      cff2*(rhoA(i-1,j)-                          &
!^   &                            rhoA(i  ,j))*                         &
!^   &                           (zwrk(i-1,j)-                          &
!^   &                            zwrk(i  ,j)))+                        &
!^   &                     (h(i-1,j)-                                   &
!^   &                      h(i  ,j))*                                  &
!^   &                     (tl_rzetaSA(i-1,j)+                          &
!^   &                      tl_rzetaSA(i  ,j)+                          &
!^   &                      cff2*((tl_rhoA(i-1,j)-                      &
!^   &                             tl_rhoA(i  ,j))*                     &
!^   &                            (zwrk(i-1,j)-                         &
!^   &                             zwrk(i  ,j))+                        &
!^   &                            (rhoA(i-1,j)-                         &
!^   &                             rhoA(i  ,j))*                        &
!^   &                            (tl_zwrk(i-1,j)-                      &
!^   &                             tl_zwrk(i  ,j))))+                   &
#endif
!^   &                     (tl_rzeta2(i-1,j)-                           &
!^   &                      tl_rzeta2(i  ,j)))
!^
            adfac=cff1*on_u(i,j)*ad_rubar(i,j)
            adfac1=adfac*(rzeta(i-1,j)-rzeta(i  ,j))
            adfac2=adfac*(h(i-1,j)+h(i  ,j))
            ad_h(i-1,j)=ad_h(i-1,j)+adfac1
            ad_h(i  ,j)=ad_h(i  ,j)+adfac1
            ad_rzeta(i-1,j)=ad_rzeta(i-1,j)+adfac2
            ad_rzeta(i  ,j)=ad_rzeta(i  ,j)-adfac2
            ad_rzeta2(i-1,j)=ad_rzeta2(i-1,j)+adfac
            ad_rzeta2(i  ,j)=ad_rzeta2(i  ,j)-adfac
#if defined VAR_RHO_2D && defined SOLVE3D
            adfac3=adfac*(rzetaSA(i-1,j)+                               &
     &                    rzetaSA(i  ,j)+                               &
     &                    cff2*(rhoA(i-1,j)-                            &
     &                          rhoA(i  ,j))*                           &
     &                         (zwrk(i-1,j)-                            &
     &                          zwrk(i  ,j)))
            adfac4=adfac2*cff2*(zwrk(i-1,j)-zwrk(i,j))
            adfac5=adfac2*cff2*(rhoA(i-1,j)-rhoA(i,j))
            ad_h(i-1,j)=ad_h(i-1,j)+adfac3
            ad_h(i  ,j)=ad_h(i  ,j)-adfac3
            ad_rzetaSA(i-1,j)=ad_rzetaSA(i-1,j)+adfac2
            ad_rzetaSA(i  ,j)=ad_rzetaSA(i  ,j)+adfac2
            ad_rhoA(i-1,j)=ad_rhoA(i-1,j)+adfac4
            ad_rhoA(i  ,j)=ad_rhoA(i  ,j)-adfac4
            ad_zwrk(i-1,j)=ad_zwrk(i-1,j)+adfac5
            ad_zwrk(i  ,j)=ad_zwrk(i  ,j)-adfac5
#endif
            ad_rubar(i,j)=0.0_r8
          END IF
        END DO
      END DO

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute fast-time-averaged fields over all barotropic time steps.
!-----------------------------------------------------------------------
!
!  Reset/initialize arrays for averaged fields during the first
!  barotropic time step. Then, accumulate it time average. Include
!  physical boundary points, but not periodic ghost points or
!  computation distributed-memory computational margins.
!
      cff1=weight(1,iif(ng),ng)
      cff2=weight(2,iif(ng),ng)
!
      IF (FIRST_2D_STEP) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            IF (j.ge.Jstr) THEN
!^            tl_DV_avg2(i,j)=cff2*tl_DVom(i,j)
!^
              ad_DVom(i,j)=ad_DVom(i,j)+cff2*ad_DV_avg2(i,j)
              ad_DV_avg2(i,j)=0.0_r8
!^            tl_DV_avg1(i,j)=0.0_r8
!^
              ad_DV_avg1(i,j)=0.0_r8
            END IF
            IF (i.ge.Istr) THEN
!^            tl_DU_avg2(i,j)=cff2*tl_DUon(i,j)
!^
              ad_DUon(i,j)=ad_DUon(i,j)+cff2*ad_DU_avg2(i,j)
              ad_DU_avg2(i,j)=0.0_r8
!^            tl_DU_avg1(i,j)=0.0_r8
!^
              ad_DU_avg1(i,j)=0.0_r8
            END IF
!^          tl_Zt_avg1(i,j)=cff1*tl_zeta(i,j,knew)
!^
            ad_zeta(i,j,knew)=ad_zeta(i,j,knew)+cff1*ad_Zt_avg1(i,j)
            ad_Zt_avg1(i,j)=0.0_r8
          END DO
        END DO
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            IF (j.ge.Jstr) THEN
!^            tl_DV_avg2(i,j)=tl_DV_avg2(i,j)+cff2*tl_DVom(i,j)
!^
              ad_DVom(i,j)=ad_DVom(i,j)+cff2*ad_DV_avg2(i,j)
            END IF
            IF (i.ge.Istr) THEN
!^            tl_DU_avg2(i,j)=tl_DU_avg2(i,j)+cff2*tl_DUon(i,j)
!^
              ad_DUon(i,j)=ad_DUon(i,j)+cff2*ad_DU_avg2(i,j)
            END IF
!^          tl_Zt_avg1(i,j)=tl_Zt_avg1(i,j)+cff1*tl_zeta(i,j,knew)
!^
            ad_zeta(i,j,knew)=ad_zeta(i,j,knew)+cff1*ad_Zt_avg1(i,j)
          END DO
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Adjoint of advance free-surface.
!-----------------------------------------------------------------------

#ifndef SOLVE3D
!
!  Save free-surface adjoint solution for IO purposes.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ad_zeta_sol(i,j)=ad_zeta(i,j,knew)
          END DO
        END DO
#endif
!
!  Load new computed free-surface into global state array.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!^        tl_zeta(i,j,knew)=tl_zeta_new(i,j)
!^
          ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_zeta(i,j,knew)
          ad_zeta(i,j,knew)=0.0_r8
        END DO
      END DO
!
!  Apply boundary conditions to newly computed free-surface "zeta_new".
!  Notice that we are using the local routine, which passes the private
!  "zeta_new" array as argument.
!
!  Here, we use the local "zetabc" since the private array "zeta_new"
!  is passed as an argument to allow computing the lateral boundary
!  conditions on the range IstrU-1:Iend and JstrV-1:Jend, so parallel
!  tile exchanges are avoided.
!
!^    CALL tl_zetabc_local (ng, tile,                                   &
!^   &                      LBi, UBi, LBj, UBj,                         &
!^   &                      IminS, ImaxS, JminS, JmaxS,                 &
!^   &                      kstp,                                       &
!^   &                      zeta, tl_zeta,                              &
!^   &                      zeta_new, tl_zeta_new)
!^
      CALL ad_zetabc_local (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      kstp,                                       &
     &                      zeta, ad_zeta,                              &
     &                      zeta_new, ad_zeta_new)
!
!  Apply mass point sources (volume vertical influx), if any.
!
!    Dsrc(is) = 2,  flow across grid cell w-face (positive or negative)
!
      IF (LwSrc(ng)) THEN
        DO is=1,Nsrc(ng)
          IF (INT(SOURCES(ng)%Dsrc(is)).eq.2) THEN
            i=SOURCES(ng)%Isrc(is)
            j=SOURCES(ng)%Jsrc(is)
            IF (((IstrR.le.i).and.(i.le.IendR)).and.                    &
     &          ((JstrR.le.j).and.(j.le.JendR))) THEN
!^            tl_zeta_new(i,j)=tl_zeta_new(i,j)+0.0_r8
!^
!             ad_zeta_new(i,j)=ad_zeta_new(i,j)+0.0_r8
            END IF
          END IF
        END DO
      END IF
!
!  Compute "zeta_new" at new time step and interpolate it half-step
!  backward, "zwrk" for the subsequent computation of barotropic
!  pressure gradient.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          fac=dtfast(ng)*pm(i,j)*pn(i,j)
#if defined VAR_RHO_2D && defined SOLVE3D
!^        tl_rzetaSA(i,j)=tl_zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))+           &
!^   &                    zwrk(i,j)*(tl_rhoS(i,j)-tl_rhoA(i,j))
!^
          adfac=zwrk(i,j)*ad_rzetaSA(i,j)
          ad_zwrk(i,j)=ad_zwrk(i,j)+                                    &
     &                 (rhoS(i,j)-rhoA(i,j))*ad_rzetaSA(i,j)
          ad_rhoS(i,j)=ad_rhoS(i,j)+adfac
          ad_rhoA(i,j)=ad_rhoA(i,j)-adfac
          ad_rzetaSA(i,j)=0.0_r8
!^        tl_rzeta2(i,j)=tl_rzeta(i,j)*zwrk(i,j)+                       &
!^   &                   rzeta(i,j)*tl_zwrk(i,j)
!^
          ad_rzeta(i,j)=ad_rzeta(i,j)+zwrk(i,j)*ad_rzeta2(i,j)
          ad_zwrk(i,j)=ad_zwrk(i,j)+rzeta(i,j)*ad_rzeta2(i,j)
          ad_rzeta2(i,j)=0.0_r8
!^        tl_rzeta(i,j)=(1.0_r8+rhoS(i,j))*tl_zwrk(i,j)+                &
!^   &                  tl_rhoS(i,j)*zwrk(i,j)
!^
          ad_rhoS(i,j)=ad_rhoS(i,j)+zwrk(i,j)*ad_rzeta(i,j)
          ad_zwrk(i,j)=ad_zwrk(i,j)+(1.0_r8+rhoS(i,j))*ad_rzeta(i,j)
          ad_rzeta(i,j)=0.0_r8
#else
!^        tl_rzeta2(i,j)=2.0_r8*tl_zwrk(i,j)*zwrk(i,j)
!^        tl_rzeta(i,j)=tl_zwrk(i,j)
!^
          ad_zwrk(i,j)=ad_zwrk(i,j)+                                    &
     &                 2.0_r8*zwrk(i,j)*ad_rzeta2(i,j)+                 &
     &                 ad_rzeta(i,j)
          ad_rzeta2(i,j)=0.0_r8
          ad_rzeta(i,j)=0.0_r8
#endif
!^        tl_zwrk(i,j)=bkw_new*tl_zeta_new(i,j)+                        &
!^   &                 bkw0*tl_zeta(i,j,kstp)+                          &
!^   &                 bkw1*tl_zeta(i,j,kbak)+                          &
!^   &                 bkw2*tl_zeta(i,j,kold)
!^
          ad_zeta_new(i,j)=ad_zeta_new(i,j)+bkw_new*ad_zwrk(i,j)
          ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+bkw0*ad_zwrk(i,j)
          ad_zeta(i,j,kbak)=ad_zeta(i,j,kbak)+bkw1*ad_zwrk(i,j)
          ad_zeta(i,j,kold)=ad_zeta(i,j,kold)+bkw2*ad_zwrk(i,j)
          ad_zwrk(i,j)=0.0_r8
#ifdef MASKING
# ifdef WET_DRY_NOT_YET
!^        tl_zeta_new(i,j)=tl_zeta_new(i,j)-                            &
!^   &                     tl_h(i,j)*(1.0_r8-rmask(i,j))
!^
# endif
!^        tl_zeta_new(i,j)=tl_zeta_new(i,j)*rmask(i,j)
!^
          ad_zeta_new(i,j)=ad_zeta_new(i,j)*rmask(i,j)
#endif
!^        tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+                           &
!^   &                     fac*(tl_DUon(i,j)-tl_DUon(i+1,j)+            &
!^   &                          tl_DVom(i,j)-tl_DVom(i,j+1))
!^
          adfac=fac*ad_zeta_new(i,j)
          ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_zeta_new(i,j)
          ad_DUon(i  ,j)=ad_DUon(i  ,j)+adfac
          ad_DUon(i+1,j)=ad_DUon(i+1,j)-adfac
          ad_DVom(i,j  )=ad_DVom(i,j  )+adfac
          ad_DVom(i,j+1)=ad_DVom(i,j+1)-adfac
          ad_zeta_new(i,j)=0.0_r8            ! HGA comment for debugging
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Adjoint of preliminary steps.
!-----------------------------------------------------------------------
!
!  Set vertically integrated mass fluxes DUon and DVom along the open
!  boundaries in such a way that the integral volume is conserved.
!
      IF (ANY(VolCons(:,ng))) THEN
!^      CALL tl_set_DUV_bc_tile (ng, tile,                              &
!^   &                           LBi, UBi, LBj, UBj,                    &
!^   &                           IminS, ImaxS, JminS, JmaxS,            &
!^   &                           krhs,                                  &
#ifdef MASKING
!^   &                           umask, vmask,                          &
#endif
!^   &                           om_v, on_u,                            &
!^   &                           ubar, vbar,                            &
!^   &                           tl_ubar, tl_vbar,                      &
!^   &                           Drhs, DUon, DVom,                      &
!^   &                           tl_Drhs, tl_DUon, tl_DVom)
!^
        CALL ad_set_DUV_bc_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           krhs,                                  &
#ifdef MASKING
     &                           umask, vmask,                          &
#endif
     &                           om_v, on_u,                            &
     &                           ubar, vbar,                            &
     &                           ad_ubar, ad_vbar,                      &
     &                           Drhs, DUon, DVom,                      &
     &                           ad_Drhs, ad_DUon, ad_DVom)
!
!  Compute integral mass flux across open boundaries and adjust
!  for volume conservation. Compute BASIC STATE value.
!  HGA: Need to resolve 'krhs' index here.
!
!^        CALL tl_obc_flux_tile (ng, tile,                              &
!^   &                           LBi, UBi, LBj, UBj,                    &
!^   &                           IminS, ImaxS, JminS, JmaxS,            &
!^   &                           knew,                                  &
# ifdef MASKING
!^   &                           umask, vmask,                          &
# endif
!^   &                           h, tl_h, om_v, on_u,                   &
!^   &                           ubar, vbar, zeta,                      &
!^   &                           tl_ubar, tl_vbar, tl_zeta)
!^
          CALL ad_obc_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           knew,                                  &
# ifdef MASKING
     &                           umask, vmask,                          &
# endif
     &                           h, ad_h, om_v, on_u,                   &
     &                           ubar, vbar, zeta,                      &
     &                           ad_ubar, ad_vbar, ad_zeta)
      END IF

#if defined DISTRIBUTE && \
    defined UV_ADV     && defined UV_C4ADVECTION && !defined SOLVE3D
!
!  In distributed-memory, the I- and J-ranges are different and a
!  special exchange is done here to avoid having three ghost points
!  for high-order numerical stencils. Notice that a private array is
!  passed below to the exchange routine. It also applies periodic
!  boundary conditions, if appropriate and no partitions in I- or
!  J-directions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!^      CALL exchange_v2d_tile (ng, tile,                               &
!^   &                          IminS, ImaxS, JminS, JmaxS,             &
!^   &                          tl_DVom)
!^
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             ad_DVom)
!^      CALL exchange_u2d_tile (ng, tile,                               &
!^   &                          IminS, ImaxS, JminS, JmaxS,             &
!^   &                          tl_DUon)
!^
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             ad_DUon)
      END IF
!^    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!^   &                    IminS, ImaxS, JminS, JmaxS,                   &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_DUon, tl_DVom)
!^
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_DUon, ad_DVom)
#endif
!
!  Compute total depth of water column and vertically integrated fluxes
!  needed for computing horizontal divergence to advance free surface
!  and for nonlinear advection terms for the barotropic momentum
!  equations.
!
#if defined DISTRIBUTE && !defined NESTING
# define IR_RANGE IstrUm2-1,Iendp2
# define JR_RANGE JstrVm2-1,Jendp2
# define IU_RANGE IstrUm1-1,Iendp2
# define JU_RANGE Jstrm1-1,Jendp2
# define IV_RANGE Istrm1-1,Iendp2
# define JV_RANGE JstrVm1-1,Jendp2
#else
# define IR_RANGE IstrUm2-1,Iendp2
# define JR_RANGE JstrVm2-1,Jendp2
# define IU_RANGE IstrUm2,Iendp2
# define JU_RANGE JstrVm2-1,Jendp2
# define IV_RANGE IstrUm2-1,Iendp2
# define JV_RANGE JstrVm2,Jendp2
#endif

      DO j=JV_RANGE
        DO i=IV_RANGE
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
!^        tl_DVom(i,j)=tl_vrhs(i,j)*cff1+                               &
!^   &                 vrhs(i,j)*tl_cff1
!^
          ad_vrhs(i,j)=ad_vrhs(i,j)+cff1*ad_DVom(i,j)
          ad_cff1=ad_cff1+vrhs(i,j)*ad_DVom(i,j)
          ad_DVom(i,j)=0.0_r8
!^        tl_vrhs(i,j)=fwd0*tl_vbar(i,j,kstp)+                          &
!^   &                 fwd1*tl_vbar(i,j,kbak)+                          &
!^   &                 fwd2*tl_vbar(i,j,kold)
!^
          ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)+fwd0*ad_vrhs(i,j)
          ad_vbar(i,j,kbak)=ad_vbar(i,j,kbak)+fwd1*ad_vrhs(i,j)
          ad_vbar(i,j,kold)=ad_vbar(i,j,kold)+fwd2*ad_vrhs(i,j)
          ad_vrhs(i,j)=0.0_r8
!^        tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i,j-1))
!^
          adfac=cff*ad_cff1
          ad_Drhs(i,j-1)=ad_Drhs(i,j-1)+adfac
          ad_Drhs(i,j  )=ad_Drhs(i,j  )+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
!^        tl_DUon(i,j)=tl_urhs(i,j)*cff1+                               &
!^   &                 urhs(i,j)*tl_cff1
!^
          ad_urhs(i,j)=ad_urhs(i,j)+cff1*ad_DUon(i,j)
          ad_cff1=ad_cff1+urhs(i,j)*ad_DUon(i,j)
          ad_DUon(i,j)=0.0_r8
!^        tl_urhs(i,j)=fwd0*tl_ubar(i,j,kstp)+                          &
!^   &                 fwd1*tl_ubar(i,j,kbak)+                          &
!^   &                 fwd2*tl_ubar(i,j,kold)
!^
          ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)+fwd0*ad_urhs(i,j)
          ad_ubar(i,j,kbak)=ad_ubar(i,j,kbak)+fwd1*ad_urhs(i,j)
          ad_ubar(i,j,kold)=ad_ubar(i,j,kold)+fwd2*ad_urhs(i,j)
          ad_urhs(i,j)=0.0_r8
!^        tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i-1,j))
!^
          adfac=cff*ad_cff1
          ad_Drhs(i-1,j)=ad_Drhs(i-1,j)+adfac
          ad_Drhs(i  ,j)=ad_Drhs(i  ,j)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
!
      DO j=JR_RANGE
        DO i=IR_RANGE
!^        tl_Drhs(i,j)=tl_h(i,j)+fwd0*tl_zeta(i,j,kstp)+                &
!^   &                           fwd1*tl_zeta(i,j,kbak)+                &
!^   &                           fwd2*tl_zeta(i,j,kold)
!^
          ad_h(i,j)=ad_h(i,j)+ad_Drhs(i,j)
          ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+fwd0*ad_Drhs(i,j)
          ad_zeta(i,j,kbak)=ad_zeta(i,j,kbak)+fwd1*ad_Drhs(i,j)
          ad_zeta(i,j,kold)=ad_zeta(i,j,kold)+fwd2*ad_Drhs(i,j)
          ad_Drhs(i,j)=0.0_r8
        END DO
      END DO

#undef IR_RANGE
#undef IU_RANGE
#undef IV_RANGE
#undef JR_RANGE
#undef JU_RANGE
#undef JV_RANGE
!
!  Deallocate local new free-surface.
!
      deallocate ( ad_zeta_new )
!
      RETURN
      END SUBROUTINE ad_step2d_tile

      END MODULE ad_step2d_mod
