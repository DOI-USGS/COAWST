#define DEBUG

      MODULE step2d_mod
!
!git $Id$
!=======================================================================
!                                                                      !
!  Solves nonlinear shallow-water primitive equations (barotropic mode)!
!  using the Generalized Forward-Backward 3rd-order Adams-Bashforth /  !
!  4th-order Adams-Moulton (FB AB3-AM4) time stepping algorithm        !
!  (Shchepetkin and McWilliams, 2005); see section 2.3 starting with   !
!  equation 2.49. In 3D applications, it perform fast-time averaging   !
!  to interact with 3D momentum equations (baroclinic mode).           !
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
!  Adapted from A.F. Shchepetkin routine "step2d_FB.F" (07-26-2022)    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef SOLVE3D
      USE mod_coupling
#endif
#ifdef DIAGNOSTICS_UV
      USE mod_diags
#endif
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
#endif
      USE mod_sources
      USE mod_stepping
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod,    ONLY : mp_exchange2d
#endif
      USE obc_volcons_mod,    ONLY : obc_flux_tile, set_DUV_bc_tile
#ifdef SOLVE3D
      USE set_depth_mod,      ONLY : set_depth
#endif
      USE u2dbc_mod,          ONLY : u2dbc_tile
      USE v2dbc_mod,          ONLY : v2dbc_tile
#ifdef WET_DRY
      USE wetdry_mod,         ONLY : wetdry_tile
#endif
      USE zetabc_mod,         ONLY : zetabc_local
!
      implicit none
!
      PRIVATE
      PUBLIC  :: step2d
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step2d (ng, tile)
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
      CALL wclock_on (ng, iNLM, 9, __LINE__, MyFile)
#endif
      CALL step2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, N(ng),                      &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  krhs(ng), kstp(ng), knew(ng),                   &
#ifdef SOLVE3D
     &                  nstp(ng), nnew(ng),                             &
#endif
#ifdef MASKING
     &                  GRID(ng) % pmask,       GRID(ng) % rmask,       &
     &                  GRID(ng) % umask,       GRID(ng) % vmask,       &
#endif
#ifdef WET_DRY
     &                  GRID(ng) % pmask_wet,   GRID(ng) % pmask_full,  &
     &                  GRID(ng) % rmask_wet,   GRID(ng) % rmask_full,  &
     &                  GRID(ng) % umask_wet,   GRID(ng) % umask_full,  &
     &                  GRID(ng) % vmask_wet,   GRID(ng) % vmask_full,  &
# ifdef SOLVE3D
     &                  GRID(ng) % rmask_wet_avg,                       &
# endif
#endif
#if (defined UV_COR && !defined SOLVE3D) || defined STEP2D_CORIOLIS
     &                  GRID(ng) % fomn,                                &
#endif
     &                  GRID(ng) % h,                                   &
     &                  GRID(ng) % om_u,        GRID(ng) % om_v,        &
     &                  GRID(ng) % on_u,        GRID(ng) % on_v,        &
     &                  GRID(ng) % pm,          GRID(ng) % pn,          &
#if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
     &                  GRID(ng) % dndx,        GRID(ng) % dmde,        &
#endif
     &                  GRID(ng) % rdrag,                               &
#if defined UV_QDRAG && !defined SOLVE3D
     &                  GRID(ng) % rdrag2,                              &
#endif
#if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
     &                  GRID(ng) % pmon_r,      GRID(ng) % pnom_r,      &
     &                  GRID(ng) % pmon_p,      GRID(ng) % pnom_p,      &
     &                  GRID(ng) % om_r,        GRID(ng) % on_r,        &
     &                  GRID(ng) % om_p,        GRID(ng) % on_p,        &
# ifdef UV_VIS2
     &                  MIXING(ng) % visc2_p,   MIXING(ng) % visc2_r,   &
# endif
# ifdef UV_VIS4
     &                  MIXING(ng) % visc4_p,   MIXING(ng) % visc4_r,   &
# endif
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
     &                  OCEAN(ng) % eq_tide,                            &
#endif
#ifndef SOLVE3D
     &                  FORCES(ng) % sustr,     FORCES(ng) % svstr,     &
# ifdef ATM_PRESS
     &                  FORCES(ng) % Pair,                              &
# endif
#else
# ifdef VAR_RHO_2D
     &                  COUPLING(ng) % rhoA,    COUPLING(ng) % rhoS,    &
# endif
     &                  COUPLING(ng) % DU_avg1, COUPLING(ng) % DU_avg2, &
     &                  COUPLING(ng) % DV_avg1, COUPLING(ng) % DV_avg2, &
     &                  COUPLING(ng) % Zt_avg1,                         &
     &                  COUPLING(ng) % rufrc,                           &
     &                  COUPLING(ng) % rvfrc,                           &
     &                  COUPLING(ng) % rufrc_bak,                       &
     &                  COUPLING(ng) % rvfrc_bak,                       &
#endif
#if defined NESTING && !defined SOLVE3D
     &                  OCEAN(ng) % DU_flux,    OCEAN(ng) % DV_flux,    &
#endif
     &                  OCEAN(ng) % ubar,       OCEAN(ng) % vbar,       &
     &                  OCEAN(ng) % zeta)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 9, __LINE__, MyFile)
#endif
!
      RETURN
      END SUBROUTINE step2d
!
!***********************************************************************
      SUBROUTINE step2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, UBk,                  &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        krhs, kstp, knew,                         &
#ifdef SOLVE3D
     &                        nstp, nnew,                               &
#endif
#ifdef MASKING
     &                        pmask, rmask, umask, vmask,               &
#endif
#ifdef WET_DRY
     &                        pmask_wet, pmask_full,                    &
     &                        rmask_wet, rmask_full,                    &
     &                        umask_wet, umask_full,                    &
     &                        vmask_wet, vmask_full,                    &
# ifdef SOLVE3D
     &                        rmask_wet_avg,                            &
# endif
#endif
#if (defined UV_COR && !defined SOLVE3D) || defined STEP2D_CORIOLIS
     &                        fomn,                                     &
#endif
     &                        h,                                        &
     &                        om_u, om_v, on_u, on_v, pm, pn,           &
#if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
     &                        dndx, dmde,                               &
#endif
     &                        rdrag,                                    &
#if defined UV_QDRAG && !defined SOLVE3D
     &                        rdrag2,                                   &
#endif
#if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
     &                        pmon_r, pnom_r, pmon_p, pnom_p,           &
     &                        om_r, on_r, om_p, on_p,                   &
# ifdef UV_VIS2
     &                        visc2_p, visc2_r,                         &
# endif
# ifdef UV_VIS4
     &                        visc4_p, visc4_r,                         &
# endif
#endif
#if defined SEDIMENT && defined SED_MORPH
     &                        bed_thick,                                &
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
     &                        eq_tide,                                  &
#endif
#ifndef SOLVE3D
     &                        sustr, svstr,
# ifdef ATM_PRESS
     &                        Pair,                                     &
# endif
#else
# ifdef VAR_RHO_2D
     &                        rhoA, rhoS,                               &
# endif
     &                        DU_avg1, DU_avg2,                         &
     &                        DV_avg1, DV_avg2,                         &
     &                        Zt_avg1,                                  &
     &                        rufrc, rvfrc,                             &
     &                        rufrc_bak, rvfrc_bak,                     &
#endif
#ifdef DIAGNOSTICS_UV
     &                        DiaU2wrk, DiaV2wrk,                       &
     &                        DiaRUbar, DiaRVbar,                       &
# ifdef SOLVE3D
     &                        DiaU2int, DiaV2int,                       &
     &                        DiaRUfrc, DiaRVfrc,                       &
# endif
#endif
#if defined NESTING && !defined SOLVE3D
     &                        DU_flux, DV_flux,                         &
#endif
     &                        ubar, vbar, zeta)
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
# if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: h(LBi:,LBj:)
# else
      real(r8), intent(in   ) :: h(LBi:,LBj:)
# endif
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
# if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in   ) :: bed_thick(LBi:,LBj:,:)
# endif
# if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
      real(r8), intent(in   ) :: eq_tide(LBi:,LBj:)
# endif
# ifndef SOLVE3D
      real(r8), intent(in   ) :: sustr(LBi:,LBj:)
      real(r8), intent(in   ) :: svstr(LBi:,LBj:)
#  ifdef ATM_PRESS
      real(r8), intent(in   ) :: Pair(LBi:,LBj:)
#  endif
# else
#  ifdef VAR_RHO_2D
      real(r8), intent(in   ) :: rhoA(LBi:,LBj:)
      real(r8), intent(in   ) :: rhoS(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: DU_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: DV_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: rufrc_bak(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvfrc_bak(LBi:,LBj:,:)
# endif
# ifdef WET_DRY
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
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
# if defined NESTING && !defined SOLVE3D
      real(r8), intent(out  ) :: DU_flux(LBi:,LBj:)
      real(r8), intent(out  ) :: DV_flux(LBi:,LBj:)
# endif

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
# if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in   ) :: h(LBi:UBi,LBj:UBj)
# endif
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
# if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in   ) :: bed_thick(LBi:UBi,LBj:UBj,1:3)
# endif
# if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
      real(r8), intent(in   ) :: eq_tide(LBi:UBi,LBj:UBj)
# endif
# ifndef SOLVE3D
      real(r8), intent(in   ) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: svstr(LBi:UBi,LBj:UBj)
#  ifdef ATM_PRESS
      real(r8), intent(in   ) :: Pair(LBi:UBi,LBj:UBj)
#  endif
# else
#  ifdef VAR_RHO_2D
      real(r8), intent(in   ) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: rhoS(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rufrc_bak(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: rvfrc_bak(LBi:UBi,LBj:UBj,2)
# endif
# ifdef WET_DRY
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
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,:)
# if defined NESTING && !defined SOLVE3D
      real(r8), intent(out  ) :: DU_flux(LBi:UBi,LBj:UBj)
      real(r8), intent(out  ) :: DV_flux(LBi:UBi,LBj:UBj)
# endif
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
#ifdef WET_DRY
      real(r8) :: cff5, cff6, cff7
#endif
      real(r8) :: fac, fac1, fac2
#ifdef DEBUG
      real(r8), parameter :: IniVal = 0.0_r8
#endif
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
#if defined STEP2D_CORIOLIS || !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
#endif
#if !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
#endif
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
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zwrk
#ifdef WET_DRY
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wetdry
#endif
#ifdef DIAGNOSTICS_UV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,NDM2d-1) :: DiaU2rhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,NDM2d-1) :: DiaV2rhs
#endif
!
      real(r8), allocatable :: zeta_new(:,:)

#include "set_bounds.h"

#ifdef DEBUG
!
!-----------------------------------------------------------------------
!  Initialize private arrays for debugging.
!-----------------------------------------------------------------------
!
# if defined UV_C4ADVECTION && !defined SOLVE3D
      Dgrad=IniVal
# endif
      Dnew=IniVal
      Dnew_rd=IniVal
      Drhs=IniVal
# if (defined UV_VIS2 || defined UV_VIS4) && !defined SOLVE3D
      Drhs_p=IniVal
# endif
      Dstp=IniVal
      DUon=IniVal
      DVom=IniVal
# if defined STEP2D_CORIOLIS || !defined SOLVE3D
      UFx=IniVal
      VFe=IniVal
# endif
# if !defined SOLVE3D
      UFe=IniVal
      VFx=IniVal
# endif
# if defined UV_C4ADVECTION && !defined SOLVE3D
      grad=IniVal
# endif
      rzeta2=IniVal
# if defined VAR_RHO_2D && defined SOLVE3D
      rzetaSA=IniVal
# endif
      rzeta=IniVal
      rubar=IniVal
      rvbar=IniVal
      urhs=IniVal
      vrhs=IniVal
      zwrk=IniVal
# ifdef WET_DRY
      wetdry=IniVal
# endif
# ifdef DIAGNOSTICS_UV
      Uwrk=IniVal
      Vwrk=IniVal
      DiaU2rhs=IniVal
      DiaV2rhs=IniVal
# endif
#endif
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
!  Preliminary steps.
!-----------------------------------------------------------------------
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

      DO j=JR_RANGE
        DO i=IR_RANGE
          Drhs(i,j)=h(i,j)+fwd0*zeta(i,j,kstp)+                         &
     &                     fwd1*zeta(i,j,kbak)+                         &
     &                     fwd2*zeta(i,j,kold)
        END DO
      END DO
!
      DO j=JU_RANGE
        DO i=IU_RANGE
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
          urhs(i,j)=fwd0*ubar(i,j,kstp)+                                &
     &              fwd1*ubar(i,j,kbak)+                                &
     &              fwd2*ubar(i,j,kold)
          DUon(i,j)=urhs(i,j)*cff1
        END DO
      END DO
!
      DO j=JV_RANGE
        DO i=IV_RANGE
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
          vrhs(i,j)=fwd0*vbar(i,j,kstp)+                                &
     &              fwd1*vbar(i,j,kbak)+                                &
     &              fwd2*vbar(i,j,kold)
          DVom(i,j)=vrhs(i,j)*cff1
        END DO
      END DO

#undef IR_RANGE
#undef IU_RANGE
#undef IV_RANGE
#undef JR_RANGE
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
!  Set vertically integrated mass fluxes DUon and DVom along the open
!  boundaries in such a way that the integral volume is conserved.
!  HGA: Need to resolve 'krhs' index here.
!
      IF (ANY(VolCons(:,ng))) THEN
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
!  Advance free-surface.
!-----------------------------------------------------------------------
!
!  Notice that the new local free-surface is allocated so it can be
!  passed as an argumment to "zetabc_local". An automatic array cannot
!  be used here because of weird memory problems.
!
      allocate ( zeta_new(IminS:ImaxS,JminS:JmaxS) )
      zeta_new = 0.0_r8
!
!  Compute "zeta_new" at new time step and interpolate it half-step
!  backward, "zwrk" for the subsequent computation of the tangent
!  linear barotropic pressure gradient. Here, we use the BASIC STATE
!  values. Thus, the nonlinear correction to the pressure-gradient
!  term from "kstp" to "knew" is not needed for Forward-Euler to
!  Forward-Backward steps (PGF_FB_CORRECTION method).
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          fac=dtfast(ng)*pm(i,j)*pn(i,j)
          zeta_new(i,j)=zeta(i,j,kstp)+                                 &
     &                  fac*(DUon(i,j)-DUon(i+1,j)+                     &
     &                       DVom(i,j)-DVom(i,j+1))
#ifdef MASKING
          zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
# ifdef WET_DRY
          zeta_new(i,j)=zeta_new(i,j)+                                  &
     &                  (Dcrit(ng)-h(i,j))*(1.0_r8-rmask(i,j))
# endif
#endif
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
              zeta_new(i,j)=zeta_new(i,j)+                              &
     &                      SOURCES(ng)%Qbar(is)*                       &
     &                      pm(i,j)*pn(i,j)*dtfast(ng)
            END IF
          END IF
        END DO
      END IF
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
      CALL zetabc_local (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   kstp,                                          &
     &                   zeta,                                          &
     &                   zeta_new)
!
!  Load new computed free-surface into global state array.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          zeta(i,j,knew)=zeta_new(i,j)
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
            Zt_avg1(i,j)=cff1*zeta(i,j,knew)
            IF (i.ge.Istr) THEN
              DU_avg1(i,j)=0.0_r8
              DU_avg2(i,j)=cff2*DUon(i,j)
            END IF
            IF (j.ge.Jstr) THEN
              DV_avg1(i,j)=0.0_r8
              DV_avg2(i,j)=cff2*DVom(i,j)
            END IF
          END DO
        END DO
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
            IF (i.ge.Istr) THEN
              DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
            END IF
            IF (j.ge.Jstr) THEN
              DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
            END IF
          END DO
        END DO
      END IF
#endif
!
!=======================================================================
!  Compute right-hand-side for the 2D momentum equations.
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
!
!-----------------------------------------------------------------------
!  Compute pressure-gradient terms.
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
          IF (i.ge.IstrU) THEN
            rubar(i,j)=cff1*on_u(i,j)*                                  &
     &                 ((h(i-1,j)+                                      &
     &                   h(i  ,j))*                                     &
     &                  (rzeta(i-1,j)-                                  &
     &                   rzeta(i  ,j))+                                 &
#if defined VAR_RHO_2D && defined SOLVE3D
     &                  (h(i-1,j)-                                      &
     &                   h(i  ,j))*                                     &
     &                  (rzetaSA(i-1,j)+                                &
     &                   rzetaSA(i  ,j)+                                &
     &                   cff2*(rhoA(i-1,j)-                             &
     &                         rhoA(i  ,j))*                            &
     &                        (zwrk(i-1,j)-                             &
     &                         zwrk(i,j)))+                             &
#endif
     &                  (rzeta2(i-1,j)-                                 &
     &                   rzeta2(i  ,j)))
#if defined ATM_PRESS && !defined SOLVE3D
            rubar(i,j)=rubar(i,j)-                                      &
     &                 cff3*on_u(i,j)*                                  &
     &                 (h(i-1,j)+h(i,j)+                                &
     &                  rzeta(i-1,j)+rzeta(i,j))*                       &
     &                 (Pair(i,j)-Pair(i-1,j))
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
            rubar(i,j)=rubar(i,j)-                                      &
     &                 cff1*on_u(i,j)*                                  &
     &                 (h(i-1,j)+h(i,j)+                                &
     &                  rzeta(i-1,j)+rzeta(i,j))*                       &
     &                 (eq_tide(i,j)-eq_tide(i-1,j))
#endif
#ifdef DIAGNOSTICS_UV
            DiaU2rhs(i,j,M2pgrd)=rubar(i,j)
#endif
          END IF
!
          IF (j.ge.JstrV) THEN
            rvbar(i,j)=cff1*om_v(i,j)*                                  &
     &                 ((h(i,j-1)+                                      &
     &                   h(i,j  ))*                                     &
     &                  (rzeta(i,j-1)-                                  &
     &                   rzeta(i,j  ))+                                 &
#if defined VAR_RHO_2D && defined SOLVE3D
     &                  (h(i,j-1)-                                      &
     &                   h(i,j  ))*                                     &
     &                  (rzetaSA(i,j-1)+                                &
     &                   rzetaSA(i,j  )+                                &
     &                   cff2*(rhoA(i,j-1)-                             &
     &                         rhoA(i,j  ))*                            &
     &                        (zwrk(i,j-1)-                             &
     &                         zwrk(i,j  )))+                           &
#endif
     &                  (rzeta2(i,j-1)-                                 &
     &                   rzeta2(i,j  )))
#if defined ATM_PRESS && !defined SOLVE3D
            rvbar(i,j)=rvbar(i,j)-                                      &
     &                 cff3*om_v(i,j)*                                  &
     &                 (h(i,j-1)+h(i,j)+                                &
     &                  rzeta(i,j-1)+rzeta(i,j))*                       &
     &                 (Pair(i,j)-Pair(i,j-1))
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
            rvbar(i,j)=rvbar(i,j)-                                      &
     &                 cff1*om_v(i,j)*                                  &
     &                 (h(i,j-1)+h(i,j)+                                &
     &                  rzeta(i,j-1)+rzeta(i,j))*                       &
     &                 (eq_tide(i,j)-eq_tide(i,j-1))
#endif
#ifdef DIAGNOSTICS_UV
            DiaV2rhs(i,j,M2pgrd)=rvbar(i,j)
#endif
          END IF
        END DO
      END DO

#if defined UV_ADV && !defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Add in horizontal advection of momentum.
!-----------------------------------------------------------------------

# ifdef UV_C2ADVECTION
!
!  Second-order, centered differences advection fluxes.
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend
          IF (i.ge.IstrU-1) THEN
            UFx(i,j)=0.25_r8*                                           &
     &               (DUon(i,j)+DUon(i+1,j))*                           &
     &               (urhs(i  ,j)+                                      &
     &                urhs(i+1,j))
          END IF
!
          VFx(i+1,j)=0.25_r8*                                           &
#  ifdef MASKING
     &               pmask(i+1,j)*                                      &
#  endif
     &               (DUon(i+1,j)+DUon(i+1,j-1))*                       &
     &               (vrhs(i+1,j)+                                      &
     &                vrhs(i  ,j))
        END DO
      END DO
!
      DO j=Jstr-1,Jend
        DO i=Istr,Iend
          IF (j.ge.JstrV-1) THEN
            VFe(i,j)=0.25_r8*                                           &
     &               (DVom(i,j)+DVom(i,j+1))*                           &
     &               (vrhs(i,j  )+                                      &
     &                vrhs(i,j+1))
          END IF
!
          UFe(i,j+1)=0.25_r8*                                           &
#  ifdef MASKING
     &               pmask(i,j+1)*                                      &
#  endif
     &               (DVom(i,j+1)+DVom(i-1,j+1))*                       &
     &               (urhs(i,j+1)+                                      &
     &                urhs(i,j  ))
        END DO
      END DO

# elif defined UV_C4ADVECTION
!
!  Fourth-order, centered differences u-momentum advection fluxes.
!
      DO j=Jstr,Jend
        DO i=IstrUm1,Iendp1
          grad (i,j)=urhs(i-1,j)-2.0_r8*urhs(i,j)+                      &
     &               urhs(i+1,j)
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
          UFx(i,j)=0.25_r8*(urhs(i  ,j)+                                &
     &                      urhs(i+1,j)-                                &
     &                      cff*(grad (i,j)+grad (i+1,j)))*             &
     &                     (DUon(i,j)+DUon(i+1,j)-                      &
     &                      cff*(Dgrad(i,j)+Dgrad(i+1,j)))
        END DO
      END DO
!
      DO j=Jstrm1,Jendp1
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
          UFe(i,j)=0.25_r8*(urhs(i,j  )+                                &
     &                      urhs(i,j-1)-                                &
     &                      cff*(grad (i,j)+grad (i,j-1)))*             &
     &                     (DVom(i,j)+DVom(i-1,j)-                      &
     &                      cff*(Dgrad(i,j)+Dgrad(i-1,j)))
        END DO
      END DO
!
!  Fourth-order, centered differences v-momentum advection fluxes.
!
      DO j=JstrV,Jend
        DO i=Istrm1,Iendp1
          grad(i,j)=vrhs(i-1,j)-2.0_r8*vrhs(i,j)+                       &
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
          VFx(i,j)=0.25_r8*(vrhs(i  ,j)+                                &
     &                      vrhs(i-1,j)-                                &
     &                      cff*(grad (i,j)+grad (i-1,j)))*             &
     &                     (DUon(i,j)+DUon(i,j-1)-                      &
     &                      cff*(Dgrad(i,j)+Dgrad(i,j-1)))
        END DO
      END DO
!
      DO j=JstrVm1,Jendp1
        DO i=Istr,Iend
          grad(i,j)=vrhs(i,j-1)-2.0_r8*vrhs(i,j)+                       &
     &              vrhs(i,j+1)
          Dgrad(i,j)=DVom(i,j-1)-2.0_r8*DVom(i,j)+DVom(i,j+1)
        END DO
      END DO
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            grad (i,Jstr)=grad (i,Jstr+1)
            Dgrad(i,Jstr)=Dgrad(i,Jstr+1)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            grad (i,Jend+1)=grad (i,Jend)
            Dgrad(i,Jend+1)=Dgrad(i,Jend)
          END DO
        END IF
      END IF
!                                                         d/dy(Dvv/m)
      cff=1.0_r8/6.0_r8
      DO j=JstrV-1,Jend
        DO i=Istr,Iend
          VFe(i,j)=0.25_r8*(vrhs(i,j  )+                                &
     &                      vrhs(i,j+1)-                                &
     &                      cff*(grad (i,j)+grad (i,j+1)))*             &
     &                     (DVom(i,j)+DVom(i,j+1)-                      &
     &                      cff*(Dgrad(i,j)+Dgrad(i,j+1)))
        END DO
      END DO
# endif
!
!  Add advection to RHS terms.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (i.ge.IstrU) THEN
            cff1=UFx(i,j)-UFx(i-1,j)
            cff2=UFe(i,j+1)-UFe(i,j)
            fac1=cff1+cff2
            rubar(i,j)=rubar(i,j)-fac1
# if defined DIAGNOSTICS_UV
            DiaU2rhs(i,j,M2xadv)=-cff1
            DiaU2rhs(i,j,M2yadv)=-cff2
            DiaU2rhs(i,j,M2hadv)=-fac1
# endif
          END IF
!
          IF (j.ge.JstrV) THEN
            cff3=VFx(i+1,j)-VFx(i,j)
            cff4=VFe(i,j)-VFe(i,j-1)
            fac2=cff3+cff4
            rvbar(i,j)=rvbar(i,j)-fac2
# if defined DIAGNOSTICS_UV
            DiaV2rhs(i,j,M2xadv)=-cff3
            DiaV2rhs(i,j,M2yadv)=-cff4
            DiaV2rhs(i,j,M2hadv)=-fac2
# endif
          END IF
        END DO
      END DO
#endif

#if (defined UV_COR && !defined SOLVE3D) || defined STEP2D_CORIOLIS
!
!-----------------------------------------------------------------------
!  Add in Coriolis term.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff=0.5_r8*Drhs(i,j)*fomn(i,j)
          UFx(i,j)=cff*(vrhs(i,j  )+                                    &
     &                  vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i  ,j)+                                    &
     &                  urhs(i+1,j))
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (i.ge.IstrU) THEN
            fac1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            rubar(i,j)=rubar(i,j)+fac1
# if defined DIAGNOSTICS_UV
            DiaU2rhs(i,j,M2fcor)=fac1
# endif

          END IF
!
          IF (j.ge.JstrV) THEN
            fac2=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rvbar(i,j)=rvbar(i,j)-fac2
# if defined DIAGNOSTICS_UV
            DiaV2rhs(i,j,M2fcor)=-fac2
# endif
          END IF
        END DO
      END DO
#endif

#if (defined CURVGRID && defined UV_ADV) && !defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Add in curvilinear transformation terms.
!-----------------------------------------------------------------------
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
          UFx(i,j)=cff*cff1
          VFe(i,j)=cff*cff2
# if defined DIAGNOSTICS_UV
          cff=Drhs(i,j)*cff4
          Uwrk(i,j)=-cff*cff1                  ! ubar equation, ETA-term
          Vwrk(i,j)=-cff*cff2                  ! vbar equation, ETA-term
# endif
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (i.ge.IstrU) THEN
            fac1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            rubar(i,j)=rubar(i,j)+fac1
# if defined DIAGNOSTICS_UV
            fac2=0.5_r8*(Uwrk(i,j)+Uwrk(i-1,j))
            DiaU2rhs(i,j,M2xadv)=DiaU2rhs(i,j,M2xadv)+fac1-fac2
            DiaU2rhs(i,j,M2yadv)=DiaU2rhs(i,j,M2yadv)+fac2
            DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+fac1
# endif
          END IF
!
          IF (j.ge.JstrV) THEN
            fac1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rvbar(i,j)=rvbar(i,j)-fac1
# if defined DIAGNOSTICS_UV
            fac2=0.5_r8*(Vwrk(i,j)+Vwrk(i,j-1))
            DiaV2rhs(i,j,M2xadv)=DiaV2rhs(i,j,M2xadv)-fac1+fac2
            DiaV2rhs(i,j,M2yadv)=DiaV2rhs(i,j,M2yadv)-fac2
            DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)-fac1
# endif
          END IF
        END DO
      END DO
#endif

#if defined UV_VIS2 && !defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Add in horizontal harmonic viscosity.
!-----------------------------------------------------------------------
!
!  Compute total depth at PSI-points.
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
          Drhs_p(i,j)=0.25_r8*(Drhs(i,j  )+Drhs(i-1,j  )+               &
     &                         Drhs(i,j-1)+Drhs(i-1,j-1))
        END DO
      END DO
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff=visc2_r(i,j)*Drhs(i,j)*0.5_r8*                            &
     &        (pmon_r(i,j)*                                             &
     &         ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,kstp)-                 &
     &          (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,kstp))-                &
     &         pnom_r(i,j)*                                             &
     &         ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,kstp)-                 &
     &          (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,kstp)))
          UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
          VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
        END DO
      END DO
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
          cff=visc2_p(i,j)*Drhs_p(i,j)*0.5_r8*                          &
     &        (pmon_p(i,j)*                                             &
     &         ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,kstp)-               &
     &          (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,kstp))+              &
     &         pnom_p(i,j)*                                             &
     &         ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,kstp)-               &
     &          (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,kstp)))
#  ifdef MASKING
          cff=cff*pmask(i,j)
#  endif
#  ifdef WET_DRY
          cff=cff*pmask_wet(i,j)
#  endif
          UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
          VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
        END DO
      END DO
!
!  Add in harmonic viscosity.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (i.ge.IstrU) THEN
            cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*(UFx(i,j  )-UFx(i-1,j))
            cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i  ,j))
            fac1=cff1+cff2
            rubar(i,j)=rubar(i,j)+fac1
# if defined DIAGNOSTICS_UV
            DiaU2rhs(i,j,M2hvis)=fac1
            DiaU2rhs(i,j,M2xvis)=cff1
            DiaU2rhs(i,j,M2yvis)=cff2
# endif
          END IF
!
          IF (j.ge.JstrV) THEN
            cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j  ))
            cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*(VFe(i  ,j)-VFe(i,j-1))
            fac1=cff1-cff2
            rvbar(i,j)=rvbar(i,j)+fac1
# if defined DIAGNOSTICS_UV
            DiaV2rhs(i,j,M2hvis)=fac1
            DiaV2rhs(i,j,M2xvis)= cff1
            DiaV2rhs(i,j,M2yvis)=-cff2
# endif
          END IF
        END DO
      END DO
#endif

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
        DO j=Jstr,Jend
          DO i=Istr,Iend
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
            IF (i.ge.IstrU) THEN
              rufrc(i,j)=rufrc(i,j)+                                    &
     &                   0.5_r8*(rdrag(i,j)+rdrag(i-1,j))*              &
     &                   om_u(i,j)*on_u(i,j)*ubar(i,j,kstp)
            END IF
!
            IF (j.ge.JstrV) THEN
              rvfrc(i,j)=rvfrc(i,j)+                                    &
     &                   0.5_r8*(rdrag(i,j)+rdrag(i,j-1))*              &
     &                   om_v(i,j)*on_v(i,j)*vbar(i,j,kstp)
            END IF
!
!  Barotropic mode running predictor stage: forward extrapolation.
!
            IF (i.ge.IstrU) THEN
              cff1=rufrc(i,j)-rubar(i,j)
              rufrc(i,j)=cfwd0*cff1+                                    &
     &                   cfwd1*rufrc_bak(i,j,  nstp)+                   &
     &                   cfwd2*rufrc_bak(i,j,3-nstp)
              rufrc_bak(i,j,3-nstp)=cff1
            END IF
!
            IF (j.ge.JstrV) THEN
              cff2=rvfrc(i,j)-rvbar(i,j)
              rvfrc(i,j)=cfwd0*cff2+                                    &
     &                   cfwd1*rvfrc_bak(i,j,  nstp)+                   &
     &                   cfwd2*rvfrc_bak(i,j,3-nstp)
              rvfrc_bak(i,j,3-nstp)=cff2
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
            zwrk(i,j)=zeta_new(i,j)-zeta(i,j,kstp)
# if defined VAR_RHO_2D && defined SOLVE3D
            rzeta(i,j)=(1.0_r8+rhoS(i,j))*zwrk(i,j)
            rzeta2(i,j)=rzeta(i,j)*(zeta_new(i,j)+zeta(i,j,kstp))
            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
            rzeta(i,j)=zwrk(i,j)
            rzeta2(i,j)=zwrk(i,j)*(zeta_new(i,j)+zeta(i,j,kstp))
# endif
          END DO
        END DO
!
        cff1=0.5*g
# if defined VAR_RHO_2D && defined SOLVE3D
        cff2=0.333333333333_r8
# endif
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (i.ge.IstrU) THEN
              rubar(i,j)=rubar(i,j)+                                    &
     &                   cff1*on_u(i,j)*                                &
     &                   ((h(i-1,j)+                                    &
     &                     h(i  ,j))*                                   &
     &                    (rzeta(i-1,j)-                                &
     &                     rzeta(i  ,j))+                               &
# if defined VAR_RHO_2D && defined SOLVE3D
     &                    (h(i-1,j)-                                    &
     &                     h(i  ,j))*                                   &
     &                    (rzetaSA(i-1,j)+                              &
     &                     rzetaSA(i  ,j)+                              &
     &                     cff2*(rhoA(i-1,j)-                           &
     &                           rhoA(i  ,j))*                          &
     &                          (zwrk(i-1,j)-                           &
     &                           zwrk(i  ,j)))+                         &
# endif
     &                    (rzeta2(i-1,j)-                               &
     &                     rzeta2(i  ,j)))
# ifdef DIAGNOSTICS_UV
              DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+                &
     &                             rubar(i,j)
# endif
            END IF
!
            IF (j.ge.JstrV) THEN
              rvbar(i,j)=rvbar(i,j)+                                    &
     &                   cff1*om_v(i,j)*                                &
     &                   ((h(i,j-1)+                                    &
     &                     h(i,j  ))*                                   &
     &                    (rzeta(i,j-1)-                                &
     &                     rzeta(i,j  ))+                               &
# if defined VAR_RHO_2D && defined SOLVE3D
     &                    (h(i,j-1)-                                    &
     &                     h(i,j  ))*                                   &
     &                    (rzetaSA(i,j-1)+                              &
     &                     rzetaSA(i,j  )+                              &
     &                     cff2*(rhoA(i,j-1)-                           &
     &                           rhoA(i,j  ))*                          &
     &                          (zwrk(i,j-1)-                           &
     &                           zwrk(i,j  )))+                         &
# endif
     &                    (rzeta2(i,j-1)-                               &
     &                     rzeta2(i,j  )))
# ifdef DIAGNOSTICS_UV
              DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+                &
     &                             rvbar(i,j)
# endif
            END IF
          END DO
        END DO
      END IF COUPLED_STEP
#endif
!
!-----------------------------------------------------------------------
!  Time step 2D momentum equations.
!-----------------------------------------------------------------------
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
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          Dnew(i,j)=h(i,j)+zeta_new(i,j)
          Dnew_rd(i,j)=Dnew(i,j)+dtfast(ng)*rdrag(i,j)
          Dstp(i,j)=h(i,j)+zeta(i,j,kstp)
        END DO
      END DO

#if defined UV_QDRAG && !defined SOLVE3D
!
!  Add quadratic drag term associated in shallow-water applications.
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
          Dnew_rd(i,j)=Dnew_rd(i,j)+                                    &
     &                 cff*rdrag2(i,j)*cff2
        END DO
      END DO
#endif
!
!  Step 2D momentum equations.
!
      cff=0.5_r8*dtfast(ng)
#ifdef SOLVE3D
      cff1=0.5_r8*weight(1,iif(ng),ng)
#else
      cff2=2.0_r8*dtfast(ng)
#endif
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff3=cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
          fac1=1.0_r8/(Dnew_rd(i,j)+Dnew_rd(i-1,j))
          ubar(i,j,knew)=fac1*((Dstp(i,j)+Dstp(i-1,j))*ubar(i,j,kstp)+  &
#ifdef SOLVE3D
     &                         cff3*(rubar(i,j)+rufrc(i,j)))
#else
     &                         cff3*rubar(i,j)+cff2*sustr(i,j))
#endif
#ifdef MASKING
          ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
#endif
#ifdef WET_DRY
          cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
          cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
          cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
          ubar(i,j,knew)=ubar(i,j,knew)*cff7
#endif
#ifdef SOLVE3D
          DU_avg1(i,j)=DU_avg1(i,j)+                                    &
     &                 cff1*on_u(i,j)*                                  &
     &                 (Dnew(i,j)+Dnew(i-1,j))*ubar(i,j,knew)
#endif
#if defined NESTING && !defined SOLVE3D
          DU_flux(i,j)=0.5_r8*on_u(i,j)*                                &
     &                 (Dnew(i,j)+Dnew(i-1,j))*ubar(i,j,knew)
#endif
        END DO
      END DO
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff3=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          fac2=1.0_r8/(Dnew_rd(i,j)+Dnew_rd(i,j-1))
          vbar(i,j,knew)=fac2*((Dstp(i,j)+Dstp(i,j-1))*vbar(i,j,kstp)+  &
#ifdef SOLVE3D
     &                         cff3*(rvbar(i,j)+rvfrc(i,j)))
#else
     &                         cff3*rvbar(i,j)+cff2*svstr(i,j))
#endif
#ifdef MASKING
          vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
#endif
#ifdef WET_DRY
          cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
          cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
          cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
          vbar(i,j,knew)=vbar(i,j,knew)*cff7
#endif
#ifdef SOLVE3D
          DV_avg1(i,j)=DV_avg1(i,j)+                                    &
     &                 cff1*om_v(i,j)*                                  &
     &                 (Dnew(i,j)+Dnew(i,j-1))*vbar(i,j,knew)
#endif
#if defined NESTING && !defined SOLVE3D
          DV_flux(i,j)=0.5_r8*om_v(i,j)*                                &
     &                 (Dnew(i,j)+Dnew(i,j-1))*vbar(i,j,knew)
#endif
        END DO
      END DO
!
!  Apply lateral boundary conditions.
!
      CALL u2dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 krhs, kstp, knew,                                &
     &                 ubar, vbar, zeta)
      CALL v2dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 krhs, kstp, knew,                                &
     &                 ubar, vbar, zeta)
!
!  Compute integral mass flux across open boundaries and adjust
!  for volume conservation.
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
      END IF

#if defined SOLVE3D || (defined NESTING && !defined SOLVE3D)
!
!  Set barotropic fluxes along physical boundaries.
!
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr-1,JendR
            Dnew(Istr-1,j)=h(Istr-1,j)+zeta_new(Istr-1,j)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr-1,JendR
            Dnew(Iend+1,j)=h(Iend+1,j)+zeta_new(Iend+1,j)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr-1,IendR
            Dnew(i,Jstr-1)=h(i,Jstr-1)+zeta_new(i,Jstr-1)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr-1,IendR
            Dnew(i,Jend+1)=h(i,Jend+1)+zeta_new(i,Jend+1)
          END DO
        END IF
      END IF

# ifdef SOLVE3D
!
      cff1=0.5*weight(1,iif(ng),ng)
# endif
!
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrR,JendR
# if defined NESTING && !defined SOLVE3D
            DU_flux(IstrU-1,j)=0.5_r8*on_u(IstrU-1,j)*                  &
     &                         (Dnew(IstrU-1,j)+Dnew(IstrU-2,j))*       &
     &                         ubar(IstrU-1,j,knew)
# else
            DU_avg1(IstrU-1,j)=DU_avg1(IstrU-1,j)+                      &
     &                         cff1*on_u(IstrU-1,j)*                    &
     &                         (Dnew(IstrU-1,j)+Dnew(IstrU-2,j))*       &
     &                         ubar(IstrU-1,j,knew)
# endif
          END DO
          DO j=JstrV,Jend
# if defined NESTING && !defined SOLVE3D
            DV_flux(Istr-1,j)=0.5_r8*om_v(Istr-1,j)*                    &
     &                        (Dnew(Istr-1,j)+Dnew(Istr-1,j-1))*        &
     &                        vbar(Istr-1,j,knew)
# else
            DV_avg1(Istr-1,j)=DV_avg1(Istr-1,j)+                        &
     &                        cff1*om_v(Istr-1,j)*                      &
     &                        (Dnew(Istr-1,j)+Dnew(Istr-1,j-1))*        &
     &                        vbar(Istr-1,j,knew)
# endif
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrR,JendR
# if defined NESTING && !defined SOLVE3D
            DU_flux(Iend+1,j)=0.5_r8*on_u(Iend+1,j)*                    &
     &                        (Dnew(Iend+1,j)+Dnew(Iend,j))*            &
     &                        ubar(Iend+1,j,knew)
# else
            DU_avg1(Iend+1,j)=DU_avg1(Iend+1,j)+                        &
     &                        cff1*on_u(Iend+1,j)*                      &
     &                        (Dnew(Iend+1,j)+Dnew(Iend,j))*            &
     &                        ubar(Iend+1,j,knew)
# endif
          END DO
          DO j=JstrV,Jend
# if defined NESTING && !defined SOLVE3D
            DV_flux(Iend+1,j)=0.5_r8*om_v(Iend+1,j)*                    &
     &                        (Dnew(Iend+1,j)+Dnew(Iend+1,j-1))*        &
     &                        vbar(Iend+1,j,knew)
# else
            DV_avg1(Iend+1,j)=DV_avg1(Iend+1,j)+                        &
     &                        cff1*om_v(Iend+1,j)*                      &
     &                        (Dnew(Iend+1,j)+Dnew(Iend+1,j-1))*        &
     &                        vbar(Iend+1,j,knew)
# endif
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
# if defined NESTING && !defined SOLVE3D
            DU_flux(i,Jstr-1)=0.5_r8*on_u(i,Jstr-1)*                    &
     &                        (Dnew(i,Jstr-1)+Dnew(i-1,Jstr-1))*        &
     &                        ubar(i,Jstr-1,knew)
# else
            DU_avg1(i,Jstr-1)=DU_avg1(i,Jstr-1)+                        &
     &                        cff1*on_u(i,Jstr-1)*                      &
     &                        (Dnew(i,Jstr-1)+Dnew(i-1,Jstr-1))*        &
     &                        ubar(i,Jstr-1,knew)
# endif
          END DO
          DO i=IstrR,IendR
# if defined NESTING && !defined SOLVE3D
            DV_flux(i,JstrV-1)=0.5_r8*om_v(i,JstrV-1)*                  &
     &                         (Dnew(i,JstrV-1)+Dnew(i,JstrV-2))*       &
     &                         vbar(i,JstrV-1,knew)
# else
            DV_avg1(i,JstrV-1)=DV_avg1(i,JstrV-1)+                      &
     &                         cff1*om_v(i,JstrV-1)*                    &
     &                         (Dnew(i,JstrV-1)+Dnew(i,JstrV-2))*       &
     &                         vbar(i,JstrV-1,knew)
# endif
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
# if defined NESTING && !defined SOLVE3D
            DU_flux(i,Jend+1)=0.5_r8*on_u(i,Jend+1)*                    &
     &                        (Dnew(i,Jend+1)+Dnew(i-1,Jend+1))*        &
     &                        ubar(i,Jend+1,knew)
# else
            DU_avg1(i,Jend+1)=DU_avg1(i,Jend+1)+                        &
     &                        cff1*on_u(i,Jend+1)*                      &
     &                        (Dnew(i,Jend+1)+Dnew(i-1,Jend+1))*        &
     &                        ubar(i,Jend+1,knew)
# endif
          END DO
          DO i=IstrR,IendR
# if defined NESTING && !defined SOLVE3D
            DV_flux(i,Jend+1)=0.5_r8*om_v(i,Jend+1)*                    &
     &                        (Dnew(i,Jend+1)+Dnew(i,Jend))*            &
     &                        vbar(i,Jend+1,knew)
# else
            DV_avg1(i,Jend+1)=DV_avg1(i,Jend+1)+                        &
     &                        cff1*om_v(i,Jend+1)*                      &
     &                        (Dnew(i,Jend+1)+Dnew(i,Jend))*            &
     &                        vbar(i,Jend+1,knew)
# endif
          END DO
        END IF
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Apply momentum transport point sources (like river runoff), if any.
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
              cff=1.0_r8/(on_u(i,j)*                                    &
     &                    0.5_r8*(Dnew(i-1,j)+Dnew(i,j)))
              ubar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
#ifdef SOLVE3D
              DU_avg1(i,j)=SOURCES(ng)%Qbar(is)
#endif
#if defined NESTING && !defined SOLVE3D
              DU_flux(i,j)=SOURCES(ng)%Qbar(is)
#endif
            ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
              cff=1.0_r8/(om_v(i,j)*                                    &
     &                    0.5_r8*(Dnew(i,j-1)+Dnew(i,j)))
              vbar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
#ifdef SOLVE3D
              DV_avg1(i,j)=SOURCES(ng)%Qbar(is)
#endif
#if defined NESTING && !defined SOLVE3D
              DV_flux(i,j)=SOURCES(ng)%Qbar(is)
#endif
            END IF
          END IF
        END DO
      END IF
!
!  Deallocate local new free-surface.
!
      deallocate ( zeta_new )

#ifdef WET_DRY
!
!-----------------------------------------------------------------------
!  Compute new wet/dry masks.
!-----------------------------------------------------------------------
!
      CALL wetdry_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
# ifdef MASKING
     &                  pmask, rmask, umask, vmask,                     &
# endif
     &                  h, zeta(:,:,knew),                              &
# ifdef SOLVE3D
     &                  DU_avg1, DV_avg1,                               &
     &                  rmask_wet_avg,                                  &
# endif
     &                  pmask_wet, pmask_full,                          &
     &                  rmask_wet, rmask_full,                          &
     &                  umask_wet, umask_full,                          &
     &                  vmask_wet, vmask_full)
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  At the end of the last 2D time step replace the new free-surface
!  zeta(:,:,knew) with it fast time-averaged value, Zt_avg1. Recall
!  this is state variable is the one that communicates with the 3D
!  kernel. Then, compute time-dependent depths.
!-----------------------------------------------------------------------
!
      IF (iif(ng).eq.nfast(ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            zeta(i,j,knew)=Zt_avg1(i,j)
          END DO
        END DO
        CALL set_depth (ng, tile, iNLM)
      END IF
#endif

#ifdef NESTING
# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  If nesting and after all fast time steps are completed, exchange
!  halo information to time averaged fields.
!-----------------------------------------------------------------------
!
      IF (iif(ng).eq.nfast(ng)) THEN
!
!  In nesting applications with refinement grids, we need to exchange
!  the DU_avg2 and DV_avg2 fluxes boundary information for the case
!  that a contact point is at a tile partition. Notice that in such
!  cases, we need i+1 and j+1 values for spatial/temporal interpolation.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Zt_avg1)
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            DU_avg1)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            DV_avg1)
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            DU_avg2)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            DV_avg2)
        END IF

#  ifdef DISTRIBUTE
!
        CALL mp_exchange2d (ng, tile, iNLM, 3,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Zt_avg1, DU_avg1, DV_avg1)
        CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      DU_avg2, DV_avg2)
#  endif
      END IF
# else
!
!  In nesting applications with refinement grids, we need to exchange
!  the DU_flux and DV_flux fluxes boundary information for the case
!  that a contact point is at a tile partition. Notice that in such
!  cases, we need i+1 and j+1 values for spatial/temporal interpolation.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          DU_flux)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          DV_flux)
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    DU_flux, DV_flux)
#  endif
# endif
#endif
!
!-----------------------------------------------------------------------
!  Exchange halo tile information.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zeta(:,:,knew))
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar(:,:,knew))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vbar(:,:,knew))
      END IF

#ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zeta(:,:,knew),                               &
     &                    ubar(:,:,knew),                               &
     &                    vbar(:,:,knew))
#endif
!
      RETURN
      END SUBROUTINE step2d_tile

      END MODULE step2d_mod
