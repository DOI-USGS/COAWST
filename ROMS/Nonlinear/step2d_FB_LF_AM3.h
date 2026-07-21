#define PREDICTOR_2D_STEP (knew.eq.3)

      MODULE step2d_mod
!
!git $Id$
!=======================================================================
!                                                                      !
!  Solves Nonlinear (NLM) shallow-water primitive equations predictor  !
!  (Leap-frog) and corrector (Adams-Moulton) time-stepping engine with !
!  a Forward Backward feedback.                                        !
!                                                                      !
!  The kernel formulation is based on Shchepetkin and McWilliams       !
!  (2005), equations (2.38)-(2.39) and (2.40)-(2.41).                  !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Shchepetkin, A.F. and J.C. McWilliams, 2005: The regional oceanic   !
!     modeling system (ROMS): a split-explicit, free-surface,          !
!     topography-following-coordinate oceanic model, Ocean Modelling,  !
!     9, 347-404, doi:10.1016/j.ocemod.2004.08.002.                    !
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
#if defined UV_VIS2
      USE mod_mixing
#endif
      USE mod_ncparam
      USE mod_scalars
      USE mod_ocean
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
      USE obc_volcons_mod,    ONLY : obc_flux_tile,                     &
     &                               set_DUV_bc_tile
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
!************************************************************************
      SUBROUTINE step2d (ng, tile)
!************************************************************************
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
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  kstp(ng), knew(ng),                             &
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
     &                  GRID(ng) % omn,                                 &
     &                  GRID(ng) % pm,          GRID(ng) % pn,          &
#if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
     &                  GRID(ng) % dndx,        GRID(ng) % dmde,        &
#endif
#if defined UV_VIS2 && !defined SOLVE3D
     &                  GRID(ng) % pmon_r,      GRID(ng) % pnom_r,      &
     &                  GRID(ng) % pmon_p,      GRID(ng) % pnom_p,      &
     &                  GRID(ng) % om_r,        GRID(ng) % on_r,        &
     &                  GRID(ng) % om_p,        GRID(ng) % on_p,        &
     &                  MIXING(ng) % visc2_p,                           &
     &                  MIXING(ng) % visc2_r,                           &
#endif
#if defined SEDIMENT && defined SED_MORPH
     &                  SEDBED(ng) % bed_thick,                         &
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
     &                  OCEAN(ng) % eq_tide,                            &
#endif
#ifndef SOLVE3D
     &                  FORCES(ng) % sustr,     FORCES(ng) % svstr,     &
     &                  FORCES(ng) % bustr,     FORCES(ng) % bvstr,     &
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
#ifdef DIAGNOSTICS_UV
     &                  DIAGS(ng) % DiaU2wrk,   DIAGS(ng) % DiaV2wrk,   &
     &                  DIAGS(ng) % DiaRUbar,   DIAGS(ng) % DiaRVbar,   &
# ifdef SOLVE3D
     &                  DIAGS(ng) % DiaU2int,   DIAGS(ng) % DiaV2int,   &
     &                  DIAGS(ng) % DiaRUfrc,   DIAGS(ng) % DiaRVfrc,   &
# endif
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
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        kstp, knew,                               &
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
     &                        om_u, om_v, on_u, on_v, omn, pm, pn,      &
#if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
     &                        dndx, dmde,                               &
#endif
#if defined UV_VIS2 && !defined SOLVE3D
     &                        pmon_r, pnom_r, pmon_p, pnom_p,           &
     &                        om_r, on_r, om_p, on_p,                   &
     &                        visc2_p, visc2_r,                         &
#endif
#if defined SEDIMENT && defined SED_MORPH
     &                        bed_thick,                                &
#endif
#if defined TIDE_GENERATING_FORCES && !defined SOLVE3D
     &                        eq_tide,                                  &
#endif
#ifndef SOLVE3D
     &                        sustr, svstr, bustr, bvstr,               &
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
     &                        ubar,  vbar, zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in    ) :: ng, tile
      integer, intent(in    ) :: LBi, UBi, LBj, UBj
      integer, intent(in    ) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in    ) :: kstp, knew
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
      real(r8), intent(in   ) :: omn(LBi:,LBj:)
      real(r8), intent(in   ) :: pm(LBi:,LBj:)
      real(r8), intent(in   ) :: pn(LBi:,LBj:)
# if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
      real(r8), intent(in   ) :: dndx(LBi:,LBj:)
      real(r8), intent(in   ) :: dmde(LBi:,LBj:)
# endif
# if defined UV_VIS2 && !defined SOLVE3D
      real(r8), intent(in   ) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in   ) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in   ) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in   ) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in   ) :: om_r(LBi:,LBj:)
      real(r8), intent(in   ) :: on_r(LBi:,LBj:)
      real(r8), intent(in   ) :: om_p(LBi:,LBj:)
      real(r8), intent(in   ) :: on_p(LBi:,LBj:)
      real(r8), intent(in   ) :: visc2_p(LBi:,LBj:)
      real(r8), intent(in   ) :: visc2_r(LBi:,LBj:)
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
      real(r8), intent(in   ) :: bustr(LBi:,LBj:)
      real(r8), intent(in   ) :: bvstr(LBi:,LBj:)
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
# ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaU2wrk(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaV2wrk(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaRUbar(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRVbar(LBi:,LBj:,:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: DiaU2int(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaV2int(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
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
      real(r8), intent(in   ) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pn(LBi:UBi,LBj:UBj)
# if defined CURVGRID && defined UV_ADV && !defined SOLVE3D
      real(r8), intent(in   ) :: dndx(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: dmde(LBi:UBi,LBj:UBj)
# endif
# if defined UV_VIS2 && !defined SOLVE3D
      real(r8), intent(in   ) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: visc2_r(LBi:UBi,LBj:UBj)
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
      real(r8), intent(in   ) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: bvstr(LBi:UBi,LBj:UBj)
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
# ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaU2wrk(LBi:UBi,LBj:UBj,NDM2d)
      real(r8), intent(inout) :: DiaV2wrk(LBi:UBi,LBj:UBj,NDM2d)
      real(r8), intent(inout) :: DiaRUbar(LBi:UBi,LBj:UBj,2,NDM2d-1)
      real(r8), intent(inout) :: DiaRVbar(LBi:UBi,LBj:UBj,2,NDM2d-1)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: DiaU2int(LBi:UBi,LBj:UBj,NDM2d)
      real(r8), intent(inout) :: DiaV2int(LBi:UBi,LBj:UBj,NDM2d)
      real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
      real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
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
      integer :: krhs, kbak
#ifdef DIAGNOSTICS_UV
      integer :: idiag
#endif
!
      real(r8) :: cff, cff1, cff2, cff3, cff4
#ifdef WET_DRY
      real(r8) :: cff5, cff6, cff7
#endif
      real(r8) :: fac, fac1, fac2
!
#if defined UV_C4ADVECTION && !defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dgrad
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dnew
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Drhs
#if defined UV_VIS2 && !defined SOLVE3D
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
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rubar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rvbar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rzeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rzeta2
#if defined VAR_RHO_2D && defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rzetaSA
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: urhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vrhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zwrk
#ifdef DIAGNOSTICS_UV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,NDM2d-1) :: DiaU2rhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,NDM2d-1) :: DiaV2rhs
#endif
!
      real(r8), allocatable :: zeta_new(:,:)
!
!  The following stability limits are obtained empirically using 3/4
!  degree Atlantic model configuration. In all these cases barotropic
!  mode time step is about 180...250 seconds, which is much less than
!  the inertial period. The maximum stability coefficients turned out
!  to be slightly different than predicted by linear theory, although
!  all theoretical tendencies agree with the practice. Note the nearly
!  70% gain in stability compared with LF-TR for appropriate
!  coefficients (linear theory predicts beta=0.166, epsil=0.84).
!
!!    real(r8), parameter :: gamma=0.0_r8,                              &
!!   &                       beta =0.0_r8,  epsil=0.0_r8  !--> Cu=0.818
!!    real(r8), parameter :: gamma=1.0_r8/12.0_r8,                      &
!!   &                       beta =0.0_r8,  epsil=0.0_r8  !--> Cu=0.878
!!    real(r8), parameter :: gamma=1./12.,                              &
!!                           beta =0.1_r8,  epsil=0.6_r8  !--> Cu=1.050
      real(r8), parameter :: gamma=0.0_r8,                              &
     &                       beta =0.14_r8, epsil=0.74_r8 !==> Cu=1.341

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Timestep vertically integrated (barotropic) equations.
!-----------------------------------------------------------------------
!
!  In the code below it is assumed that variables with time index "krhs"
!  are time-centered at step "n" in barotropic time during predictor
!  sub-step and "n+1/2" during corrector.
!
      IF (PREDICTOR_2D_STEP) THEN
        krhs=kstp
      ELSE
        krhs=3
      END IF
      IF (FIRST_2D_STEP) THEN
        kbak=kstp                      ! "kbak" is used as "from"
      ELSE                             ! time index for LF timestep
        kbak=3-kstp
      END IF
!
!-----------------------------------------------------------------------
!  Preliminary steps.
!-----------------------------------------------------------------------
!
!  Compute total depth of the water column and vertically integrated
!  mass fluxes, which are used in computation of free-surface elevation
!  time tendency and advection terms for the barotropic momentum
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
          Drhs(i,j)=zeta(i,j,krhs)+h(i,j)
        END DO
      END DO
      DO j=JU_RANGE
        DO i=IU_RANGE
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
          DUon(i,j)=ubar(i,j,krhs)*cff1
        END DO
      END DO
      DO j=JV_RANGE
        DO i=IV_RANGE
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
          DVom(i,j)=vbar(i,j,krhs)*cff1
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

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Fields averaged over all barotropic time steps.
!-----------------------------------------------------------------------
!
!  Notice that the index ranges here are designed to include physical
!  boundaries only. Periodic ghost points and internal mpi computational
!  margins are NOT included.
!
!  Reset all barotropic mode time-averaged arrays during the first
!  predictor step. At all subsequent time steps, accumulate averages
!  of the first kind using the DELAYED way. For example, "Zt_avg1" is
!  not summed immediately after the corrector step when computed but
!  during the subsequent predictor substep. It allows saving operations
!  because "DUon" and "DVom" are calculated anyway. The last time step
!  has a special code to add all three barotropic variables after the
!  last corrector substep.
!
      IF (PREDICTOR_2D_STEP) THEN                ! PREDICTOR STEP
        IF (FIRST_2D_STEP) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              Zt_avg1(i,j)=0.0_r8
              DU_avg1(i,j)=0.0_r8
              DV_avg1(i,j)=0.0_r8
              DU_avg2(i,j)=0.0_r8
              DV_avg2(i,j)=0.0_r8
            END DO
          END DO
        ELSE
          cff=weight(1,iif(ng)-1,ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              Zt_avg1(i,j)=Zt_avg1(i,j)+cff*zeta(i,j,krhs)
              IF (i.ge.Istr) THEN
                DU_avg1(i,j)=DU_avg1(i,j)+cff*DUon(i,j)
              END IF
              IF (j.ge.Jstr) THEN
                DV_avg1(i,j)=DV_avg1(i,j)+cff*DVom(i,j)
              END IF
            END DO
          END DO
        END IF
      ELSE                                       ! CORRECTOR STEP
        cff=weight(2,iif(ng),ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            IF (i.ge.Istr) THEN
              DU_avg2(i,j)=DU_avg2(i,j)+cff*DUon(i,j)
            END IF
            IF (j.ge.Jstr) THEN
              DV_avg2(i,j)=DV_avg2(i,j)+cff*DVom(i,j)
            END IF
          END DO
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Advance free-surface.
!-----------------------------------------------------------------------
!
!  Notice that the new local free-surface is allocated so it can be
!  passed as an argumment to "zetabc" to avoid memory issues.
!
      allocate ( zeta_new(IminS:ImaxS,JminS:JmaxS) )
      zeta_new = 0.0_r8
!
!  Compute "zeta_new" at the new time step and interpolate backward for
!  the subsequent computation of barotropic pressure-gradient terms.
!  Notice that during the predictor of the first 2D step in 3D mode,
!  the pressure gradient terms are computed using just zeta(:,:,kstp),
!  i.e., like in the Forward Euler step, rather than the more accurate
!  predictor of generalized RK2. This is to keep it consistent with the
!  computation of pressure gradient in 3D mode, which uses precisely
!  the initial value of "zeta" rather than the value changed by the
!  first barotropic predictor step. Later in this code, just after
!  "rufrc, rvfrc" are finalized, a correction term based on the
!  difference zeta_new(:,:)-zeta(:,:,kstp) to "rubar, rvbar" to make
!  them consistent with generalized RK2 stepping for pressure gradient
!  terms.
!
      IF (PREDICTOR_2D_STEP) THEN
        IF (FIRST_2D_STEP) THEN     ! Modified RK2 time step (with
          cff=dtfast(ng)            ! Forward-Backward feedback with
#ifdef SOLVE3D
          cff1=0.0_r8               !==> Forward Euler
          cff2=1.0_r8
#else
          cff1=0.333333333333_r8    ! optimally chosen beta=1/3 and
          cff2=0.666666666667_r8    ! epsilon=2/3, see below) is used
#endif
          cff3=0.0_r8               ! here for the start up.
        ELSE
          cff=2.0_r8*dtfast(ng)     ! In the code below "zwrk" is
          cff1=beta                 ! time-centered at time step "n"
          cff2=1.0_r8-2.0_r8*beta   ! in the case of LF (for all but
          cff3=beta                 ! the first time step)
        END IF
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            fac=cff*pm(i,j)*pn(i,j)
            zeta_new(i,j)=zeta(i,j,kbak)+                               &
     &                    fac*(DUon(i,j)-DUon(i+1,j)+                   &
     &                         DVom(i,j)-DVom(i,j+1))
#ifdef MASKING
            zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
# ifdef WET_DRY
            zeta_new(i,j)=zeta_new(i,j)+                                &
     &                    (Dcrit(ng)-h(i,j))*(1.0_r8-rmask(i,j))
# endif
#endif
            Dnew(i,j)=zeta_new(i,j)+h(i,j)

            zwrk(i,j)=cff1*zeta_new(i,j)+                               &
     &                cff2*zeta(i,j,kstp)+                              &
     &                cff3*zeta(i,j,kbak)
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
      ELSE                                     !--> CORRECTOR STEP
        IF (FIRST_2D_STEP) THEN
          cff =0.333333333333_r8               ! Modified RK2 weighting:
          cff1=0.333333333333_r8               ! here "zwrk" is time-
          cff2=0.333333333333_r8               ! centered at "n+1/2".
          cff3=0.0_r8
        ELSE
          cff =1.0_r8-epsil                    ! zwrk is always time-
          cff1=(0.5_r8-gamma)*epsil            ! centered at n+1/2
          cff2=(0.5_r8+2.0_r8*gamma)*epsil     ! during corrector sub-
          cff3=-gamma *epsil                   ! step.
        END IF
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            fac=dtfast(ng)*pm(i,j)*pn(i,j)
            zeta_new(i,j)=zeta(i,j,kstp)+                               &
     &                    fac*(DUon(i,j)-DUon(i+1,j)+                   &
     &                         DVom(i,j)-DVom(i,j+1))
#ifdef MASKING
            zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
#endif
            Dnew(i,j)=zeta_new(i,j)+h(i,j)

            zwrk(i,j)=cff *zeta(i,j,krhs)+                              &
     &                cff1*zeta_new(i,j)+                               &
     &                cff2*zeta(i,j,kstp)+                              &
     &                cff3*zeta(i,j,kbak)
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
      END IF
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
!  Apply boundary conditions to newly computed free-surface "zeta_new"
!  and load into global state array. Notice that "zeta_new" is always
!  centered at time step "m+1", while zeta(:,:,knew) should be centered
!  either at "m+1/2" after predictor step and at "m+1" after corrector.
!  Chosing it to be this way makes it possible avoid storing RHS for
!  zeta, ubar, and vbar between predictor and corrector sub-steps.
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
      IF (PREDICTOR_2D_STEP) THEN
        IF (FIRST_2D_STEP) THEN
          cff1=0.5_r8
          cff2=0.5_r8
          cff3=0.0_r8
        ELSE
          cff1=0.5_r8-gamma
          cff2=0.5_r8+2.0_r8*gamma
          cff3=-gamma
        END IF
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            zeta(i,j,knew)=cff1*zeta_new(i,j)+                          &
     &                     cff2*zeta(i,j,kstp)+                         &
     &                     cff3*zeta(i,j,kbak)
          END DO
        END DO
      ELSE
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            zeta(i,j,knew)=zeta_new(i,j)
          END DO
        END DO
      END IF
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
!  Notice that "rubar" and "rvbar" are computed within the same to allow
!  shared references to array elements (i,j), which increases the
!  computational density by almost a factor of 1.5 resulting in overall
!  more efficient code.
!
      cff1=0.5*g
      cff2=0.333333333333_r8
#if !defined SOLVE3D && defined ATM_PRESS
      fac=0.5_r8*100.0_r8/rho0
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
     &                         zwrk(i  ,j)))+                           &
#endif
     &                  (rzeta2(i-1,j)-                                 &
     &                   rzeta2(i  ,j)))
#if defined ATM_PRESS && !defined SOLVE3D
            rubar(i,j)=rubar(i,j)-                                      &
     &                 fac*on_u(i,j)*                                   &
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
     &                 fac*om_v(i,j)*                                   &
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
        DO i=IstrU-1,Iend
          UFx(i,j)=0.25_r8*                                             &
     &             (DUon(i,j)+DUon(i+1,j))*                             &
     &             (ubar(i  ,j,krhs)+                                   &
     &              ubar(i+1,j,krhs))
        END DO
      END DO
!
      DO j=Jstr,Jend+1
        DO i=IstrU,Iend
          UFe(i,j)=0.25_r8*                                             &
     &             (DVom(i,j)+DVom(i-1,j))*                             &
     &             (ubar(i,j  ,krhs)+                                   &
     &              ubar(i,j-1,krhs))
        END DO
      END DO
!
      DO j=JstrV,Jend
        DO i=Istr,Iend+1
          VFx(i,j)=0.25_r8*                                             &
     &             (DUon(i,j)+DUon(i,j-1))*                             &
     &             (vbar(i  ,j,krhs)+                                   &
     &              vbar(i-1,j,krhs))
        END DO
      END DO
!
      DO j=JstrV-1,Jend
        DO i=Istr,Iend
          VFe(i,j)=0.25_r8*                                             &
     &             (DVom(i,j)+DVom(i,j+1))*                             &
     &             (vbar(i,j  ,krhs)+                                   &
     &              vbar(i,j+1,krhs))
        END DO
      END DO

# elif defined UV_C4ADVECTION
!
!  Fourth-order, centered differences u-momentum advection fluxes.
!
      DO j=Jstr,Jend
        DO i=IstrUm1,Iendp1
          grad (i,j)=ubar(i-1,j,krhs)-2.0_r8*ubar(i,j,krhs)+            &
     &               ubar(i+1,j,krhs)
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
          UFx(i,j)=0.25_r8*(ubar(i  ,j,krhs)+                           &
     &                      ubar(i+1,j,krhs)-                           &
     &                      cff*(grad (i,j)+grad (i+1,j)))*             &
     &                     (DUon(i,j)+DUon(i+1,j)-                      &
     &                      cff*(Dgrad(i,j)+Dgrad(i+1,j)))
        END DO
      END DO
!
      DO j=Jstrm1,Jendp1
        DO i=IstrU,Iend
          grad(i,j)=ubar(i,j-1,krhs)-2.0_r8*ubar(i,j,krhs)+             &
     &              ubar(i,j+1,krhs)
        END DO
      END DO
!
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
          UFe(i,j)=0.25_r8*(ubar(i,j  ,krhs)+                           &
     &                      ubar(i,j-1,krhs)-                           &
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
          grad(i,j)=vbar(i-1,j,krhs)-2.0_r8*vbar(i,j,krhs)+             &
     &              vbar(i+1,j,krhs)
        END DO
      END DO
!
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
          VFx(i,j)=0.25_r8*(vbar(i  ,j,krhs)+                           &
     &                      vbar(i-1,j,krhs)-                           &
     &                      cff*(grad (i,j)+grad (i-1,j)))*             &
     &                     (DUon(i,j)+DUon(i,j-1)-                      &
     &                      cff*(Dgrad(i,j)+Dgrad(i,j-1)))
        END DO
      END DO
!
      DO j=JstrVm1,Jendp1
        DO i=Istr,Iend
          grad(i,j)=vbar(i,j-1,krhs)-2.0_r8*vbar(i,j,krhs)+             &
     &              vbar(i,j+1,krhs)
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
          VFe(i,j)=0.25_r8*(vbar(i,j  ,krhs)+                           &
     &                      vbar(i,j+1,krhs)-                           &
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
            fac=cff1+cff2
            rubar(i,j)=rubar(i,j)-fac
# if defined DIAGNOSTICS_UV
            DiaU2rhs(i,j,M2xadv)=-cff1
            DiaU2rhs(i,j,M2yadv)=-cff2
            DiaU2rhs(i,j,M2hadv)=-fac
# endif
          END IF
!
          IF (j.ge.JstrV) THEN
            cff1=VFx(i+1,j)-VFx(i,j)
            cff2=VFe(i,j)-VFe(i,j-1)
            fac=cff1+cff2
            rvbar(i,j)=rvbar(i,j)-fac
# if defined DIAGNOSTICS_UV
            DiaV2rhs(i,j,M2xadv)=-cff1
            DiaV2rhs(i,j,M2yadv)=-cff2
            DiaV2rhs(i,j,M2hadv)=-fac
# endif
          END IF
        END DO
      END DO
#endif

#if (defined UV_COR & !defined SOLVE3D) || defined STEP2D_CORIOLIS
!
!-----------------------------------------------------------------------
!  Add in Coriolis term.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff=0.5_r8*Drhs(i,j)*fomn(i,j)
          UFx(i,j)=cff*(vbar(i,j  ,krhs)+                               &
     &                  vbar(i,j+1,krhs))
          VFe(i,j)=cff*(ubar(i  ,j,krhs)+                               &
     &                  ubar(i+1,j,krhs))
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
          cff1=0.5_r8*(vbar(i,j  ,krhs)+                                &
     &                 vbar(i,j+1,krhs))
          cff2=0.5_r8*(ubar(i  ,j,krhs)+                                &
     &                 ubar(i+1,j,krhs))
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
     &         ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,krhs)-                 &
     &          (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,krhs))-                &
     &         pnom_r(i,j)*                                             &
     &         ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,krhs)-                 &
     &          (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,krhs)))
          UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
          VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
        END DO
      END DO
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
          cff=visc2_p(i,j)*Drhs_p(i,j)*0.5_r8*                          &
     &        (pmon_p(i,j)*                                             &
     &         ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,krhs)-               &
     &          (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,krhs))+              &
     &         pnom_p(i,j)*                                             &
     &         ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,krhs)-               &
     &          (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,krhs)))
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
            fac=cff1+cff2
            rubar(i,j)=rubar(i,j)+fac
# if defined DIAGNOSTICS_UV
            DiaU2rhs(i,j,M2hvis)=fac
            DiaU2rhs(i,j,M2xvis)=cff1
            DiaU2rhs(i,j,M2yvis)=cff2
# endif
          END IF
!
          IF (j.ge.JstrV) THEN
            cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j  ))
            cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*(VFe(i  ,j)-VFe(i,j-1))
            fac=cff1-cff2
            rvbar(i,j)=rvbar(i,j)+fac
# if defined DIAGNOSTICS_UV
            DiaV2rhs(i,j,M2hvis)=fac
            DiaV2rhs(i,j,M2xvis)= cff1
            DiaV2rhs(i,j,M2yvis)=-cff2
# endif
          END IF
        END DO
      END DO
#endif

#ifndef SOLVE3D
!
!-----------------------------------------------------------------------
!  Add in bottom stress.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          fac=bustr(i,j)*om_u(i,j)*on_u(i,j)
          rubar(i,j)=rubar(i,j)-fac
# ifdef DIAGNOSTICS_UV
          DiaU2rhs(i,j,M2bstr)=-fac
# endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          fac=bvstr(i,j)*om_v(i,j)*on_v(i,j)
          rvbar(i,j)=rvbar(i,j)-fac
# ifdef DIAGNOSTICS_UV
          DiaV2rhs(i,j,M2bstr)=-fac
# endif
        END DO
      END DO
#else
# ifdef DIAGNOSTICS_UV
!
!  Initialize the stress term if no bottom friction is defined.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          DiaU2rhs(i,j,M2bstr)=0.0_r8
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          DiaV2rhs(i,j,M2bstr)=0.0_r8
        END DO
      END DO
# endif
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Coupling between 2D and 3D equations.
!-----------------------------------------------------------------------
!
!  Before the predictor step of the first barotropic time step, arrays
!  "rufrc" and "rvfrc" contain vertical integrals of the 3D RHS terms
!  for the momentum equations (including surface and bottom stresses,
!  if so prescribed). During the first barotropic time step, convert
!  them into forcing terms by subtracting the fast-time "rubar" and
!  "rvbar" from them.
!
!  These forcing terms are then extrapolated forward in time using
!  optimized Adams-Bashforth weights, so that the resultant "rufrc"
!  and "rvfrc" are centered effectively at time n+1/2 in baroclinic
!  time.
!
!  From now on, these newly computed forcing terms remain unchanged
!  during the fast time stepping and will be added to "rubar" and
!  "rvbar" during all subsequent barotropic time steps.
!
!  Thus, the algorithm below is designed for coupling during the 3D
!  predictor sub-step. The forcing terms "rufrc" and "rvfrc" are
!  computed as instantaneous values at 3D time index "nstp" first and
!  then extrapolated half-step forward using AM3-like weights optimized
!  for maximum stability (with particular care for startup).
!
      IF (FIRST_2D_STEP.and.PREDICTOR_2D_STEP) THEN
        IF (FIRST_TIME_STEP) THEN
          cff3=0.0_r8
          cff2=0.0_r8
          cff1=1.0_r8
        ELSE IF (FIRST_TIME_STEP+1) THEN
          cff3=0.0_r8
          cff2=-0.5_r8
          cff1=1.5_r8
        ELSE
          cff3=0.281105_r8
          cff2=-0.5_r8-2.0_r8*cff3
          cff1=1.5_r8+cff3
        END IF
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=rufrc(i,j)-rubar(i,j)
            rufrc(i,j)=cff1*cff+                                        &
     &                 cff2*rufrc_bak(i,j,3-nstp)+                      &
     &                 cff3*rufrc_bak(i,j,nstp  )
            rufrc_bak(i,j,nstp)=cff
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=rvfrc(i,j)-rvbar(i,j)
            rvfrc(i,j)=cff1*cff+                                        &
     &                 cff2*rvfrc_bak(i,j,3-nstp)+                      &
     &                 cff3*rvfrc_bak(i,j,nstp  )
            rvfrc_bak(i,j,nstp)=cff
          END DO
        END DO
!
!  Since coupling requires that the pressure gradient term is computed
!  using zeta(:,:,kstp) instead of 1/3 toward zeta_new(:,:) as needed
!  by generalized RK2 scheme, apply compensation to shift pressure
!  gradient terms from "kstp" to 1/3 toward "knew".
!
        cff1=0.5_r8*g
        cff2=0.333333333333_r8
        cff3=1.666666666666_r8

        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            zwrk(i,j)=cff2*(zeta_new(i,j)-zeta(i,j,kstp))
# if defined VAR_RHO_2D && defined SOLVE3D
            rzeta(i,j)=(1.0_r8+rhoS(i,j))*zwrk(i,j)
            rzeta2(i,j)=rzeta(i,j)*                                     &
     &                  (cff2*zeta_new(i,j)+                            &
     &                   cff3*zeta(i,j,kstp))
            rzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
# else
            rzeta(i,j)=zwrk(i,j)
            rzeta2(i,j)=zwrk(i,j)*                                      &
     &                 (cff2*zeta_new(i,j)+                             &
     &                  cff3*zeta(i,j,kstp))
# endif
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (i.ge.IstrU) THEN
              rubar(i,j)=rubar(i,j)+                                    &
     &                   cff1*on_u(i,j)*                                &
     &                   ((h(i-1,j)+                                    &
     &                     h(i  ,j))*                                   &
     &                   (rzeta(i-1,j)-                                 &
     &                    rzeta(i  ,j))+                                &
# if defined VAR_RHO_2D && defined SOLVE3D
     &                   (h(i-1,j)-                                     &
     &                    h(i  ,j))*                                    &
     &                   (rzetaSA(i-1,j)+                               &
     &                    rzetaSA(i  ,j)+                               &
     &                    cff2*(rhoA(i-1,j)-                            &
     &                          rhoA(i  ,j))*                           &
     &                         (zwrk(i-1,j)-                            &
     &                          zwrk(i  ,j)))+                          &
# endif
     &                   (rzeta2(i-1,j)-                                &
     &                    rzeta2(i  ,j)))
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
      END IF
#endif
!
!=======================================================================
!  Time step 2D momentum equations.
!=======================================================================
!
!  Compute total water column depth.
!
      IF (FIRST_2D_STEP.or.(.not.PREDICTOR_2D_STEP)) THEN
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            Dstp(i,j)=h(i,j)+zeta(i,j,kstp)
          END DO
        END DO
      ELSE
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            Dstp(i,j)=h(i,j)+zeta(i,j,kbak)
          END DO
        END DO
      END IF
!
!  During the predictor sub-step, once newly computed "ubar" and "vbar"
!  become available, interpolate them half-step backward in barotropic
!  time (i.e., they end up time-centered at n+1/2) in order to use it
!  during subsequent corrector sub-step.
!
      IF (PREDICTOR_2D_STEP) THEN
        IF (FIRST_2D_STEP) THEN
          cff1=0.5_r8*dtfast(ng)
          cff2=0.5_r8
          cff3=0.5_r8
          cff4=0.0_r8
        ELSE
          cff1=dtfast(ng)
          cff2=0.5_r8-gamma
          cff3=0.5_r8+2.0_r8*gamma
          cff4=-gamma
        ENDIF
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=cff1*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
            fac1=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
            ubar(i,j,knew)=fac1*                                        &
     &                     (ubar(i,j,kbak)*                             &
     &                      (Dstp(i,j)+Dstp(i-1,j))+                    &
#ifdef SOLVE3D
     &                      cff*(rubar(i,j)+rufrc(i,j)))
#else
     &                      cff*rubar(i,j)+4.0_r8*cff1*sustr(i,j))
#endif
#ifdef MASKING
            ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
#endif
            ubar(i,j,knew)=cff2*ubar(i,j,knew)+                         &
     &                     cff3*ubar(i,j,kstp)+                         &
     &                     cff4*ubar(i,j,kbak)
#ifdef WET_DRY
            cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
            cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
            cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
            ubar(i,j,knew)=ubar(i,j,knew)*cff7
#endif
#if defined NESTING && !defined SOLVE3D
            DU_flux(i,j)=0.5_r8*on_u(i,j)*                              &
     &                   (Dnew(i,j)+Dnew(i-1,j))*ubar(i,j,knew)
#endif
          END DO
        END DO
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=cff1*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            fac2=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
            vbar(i,j,knew)=fac2*                                        &
     &                     (vbar(i,j,kbak)*                             &
     &                      (Dstp(i,j)+Dstp(i,j-1))+                    &
#ifdef SOLVE3D
     &                      cff*(rvbar(i,j)+rvfrc(i,j)))
#else
     &                      cff*rvbar(i,j)+4.0_r8*cff1*svstr(i,j))
#endif
#ifdef MASKING
            vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
#endif
            vbar(i,j,knew)=cff2*vbar(i,j,knew)+                         &
     &                     cff3*vbar(i,j,kstp)+                         &
     &                     cff4*vbar(i,j,kbak)
#ifdef WET_DRY
            cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
            cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
            cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
            vbar(i,j,knew)=vbar(i,j,knew)*cff7
#endif
#if defined NESTING && !defined SOLVE3D
            DV_flux(i,j)=0.5_r8*om_v(i,j)*                              &
     &                   (Dnew(i,j)+Dnew(i,j-1))*vbar(i,j,knew)
#endif
          END DO
        END DO

      ELSE                                    !--> CORRECTOR_2D_STEP

        cff1=0.5_r8*dtfast(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=cff1*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
            fac1=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
            ubar(i,j,knew)=fac1*                                        &
     &                     (ubar(i,j,kstp)*                             &
     &                      (Dstp(i,j)+Dstp(i-1,j))+                    &
#ifdef SOLVE3D
     &                      cff*(rubar(i,j)+rufrc(i,j)))
#else
     &                      cff*rubar(i,j)+4.0_r8*cff1*sustr(i,j))
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
#if defined NESTING && !defined SOLVE3D
            DU_flux(i,j)=0.5_r8*on_u(i,j)*                              &
     &                   (Dnew(i,j)+Dnew(i-1,j))*ubar(i,j,knew)
#endif
          END DO
        END DO
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=cff1*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            fac2=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
            vbar(i,j,knew)=fac2*                                        &
     &                     (vbar(i,j,kstp)*                             &
     &                      (Dstp(i,j)+Dstp(i,j-1))+                    &
#ifdef SOLVE3D
     &                      cff*(rvbar(i,j)+rvfrc(i,j)))
#else
     &                      cff*rvbar(i,j)+4.0_r8*cff1*svstr(i,j))
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
#if defined NESTING && !defined SOLVE3D
            DV_flux(i,j)=0.5_r8*om_v(i,j)*                              &
     &                   (Dnew(i,j)+Dnew(i,j-1))*vbar(i,j,knew)
#endif
          END DO
        END DO
      END IF
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

#if defined NESTING && !defined SOLVE3D
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
!
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrR,JendR
            DU_flux(IstrU-1,j)=0.5_r8*on_u(IstrU-1,j)*                  &
     &                         (Dnew(IstrU-1,j)+Dnew(IstrU-2,j))*       &
     &                         ubar(IstrU-1,j,knew)
          END DO
          DO j=JstrV,Jend
            DV_flux(Istr-1,j)=0.5_r8*om_v(Istr-1,j)*                    &
     &                        (Dnew(Istr-1,j)+Dnew(Istr-1,j-1))*        &
     &                        vbar(Istr-1,j,knew)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrR,JendR
            DU_flux(Iend+1,j)=0.5_r8*on_u(Iend+1,j)*                    &
     &                        (Dnew(Iend+1,j)+Dnew(Iend,j))*            &
     &                        ubar(Iend+1,j,knew)
          END DO
          DO j=JstrV,Jend
            DV_flux(Iend+1,j)=0.5_r8*om_v(Iend+1,j)*                    &
     &                        (Dnew(Iend+1,j)+Dnew(Iend+1,j-1))*        &
     &                        vbar(Iend+1,j,knew)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            DU_flux(i,Jstr-1)=0.5_r8*on_u(i,Jstr-1)*                    &
     &                        (Dnew(i,Jstr-1)+Dnew(i-1,Jstr-1))*        &
     &                        ubar(i,Jstr-1,knew)
          END DO
          DO i=IstrR,IendR
            DV_flux(i,JstrV-1)=0.5_r8*om_v(i,JstrV-1)*                  &
     &                         (Dnew(i,JstrV-1)+Dnew(i,JstrV-2))*       &
     &                         vbar(i,JstrV-1,knew)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
            DU_flux(i,Jend+1)=0.5_r8*on_u(i,Jend+1)*                    &
     &                        (Dnew(i,Jend+1)+Dnew(i-1,Jend+1))*        &
     &                        ubar(i,Jend+1,knew)
          END DO
          DO i=IstrR,IendR
            DV_flux(i,Jend+1)=0.5_r8*om_v(i,Jend+1)*                    &
     &                        (Dnew(i,Jend+1)+Dnew(i,Jend))*            &
     &                        vbar(i,Jend+1,knew)
          END DO
        END IF
      END IF
#endif
!
!  Apply momentum transport point sources (like river runoff), if any.
!
!    Dsrc(is) = 0,  flow across grid cell u-face (positive or negative)
!    Dsrc(is) = 1,  flow across grid cell v-face (positive or negative)
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

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Finalize computation of barotropic mode averages.
!-----------------------------------------------------------------------
!
!  This procedure starts with filling in boundary rows of total depths
!  at the new time step, which is needed to be done only during the
!  last barotropic time step, Normally, the computation of averages
!  occurs at the beginning of the next predictor step because "DUon"
!  and "DVom" are being computed anyway. Strictly speaking, the filling
!  the boundaries are necessary only in the case of open boundaries,
!  otherwise, the associated fluxes are all zeros.
!
      IF ((iif(ng).eq.nfast(ng)).and.(knew.lt.3)) THEN
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
!
!  At the end of the last 2D time step replace the new free-surface
!  zeta(:,:,knew) with it fast time-averaged value, Zt_avg1. Recall
!  this is state variable is the one that communicates with the 3D
!  kernel. Then, compute time-dependent depths.
!
        cff=weight(1,iif(ng),ng)
        cff1=0.5*cff
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Zt_avg1(i,j)=Zt_avg1(i,j)+                                  &
     &                   cff*zeta(i,j,knew)
            IF (i.ge.Istr) THEN
              DU_avg1(i,j)=DU_avg1(i,j)+                                &
     &                     cff1*on_u(i,j)*                              &
     &                     (Dnew(i,j)+Dnew(i-1,j))*ubar(i,j,knew)
            END IF
            IF (j.ge.Jstr) THEN
              DV_avg1(i,j)=DV_avg1(i,j)+                                &
     &                     cff1*om_v(i,j)*                              &
     &                     (Dnew(i,j)+Dnew(i,j-1))*vbar(i,j,knew)
            END IF
            zeta(i,j,knew)=Zt_avg1(i,j)
          END DO
        END DO
        CALL set_depth (ng, tile, iNLM)

# ifdef NESTING
!
!  After all fast time steps are completed, apply boundary conditions
!  to time averaged fields.
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
# endif
      END IF
#endif
#if defined NESTING && !defined SOLVE3D
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

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    DU_flux, DV_flux)
# endif
#endif
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
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
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
!
      END MODULE step2d_mod
