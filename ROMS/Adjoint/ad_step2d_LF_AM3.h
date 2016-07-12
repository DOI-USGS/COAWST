#undef DEBUG
#ifdef ADJOINT
      SUBROUTINE ad_step2d (ng, tile)
!
!svn $Id: ad_step2d_LF_AM3.h 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!                                                                      !
!  Adjoint shallow-water primitive equations predictor (Leap-frog)     !
!  and corrector (Adams-Moulton) time-stepping engine.                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
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
# if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
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
# ifdef MASKING
     &                     GRID(ng) % pmask,       GRID(ng) % rmask,    &
     &                     GRID(ng) % umask,       GRID(ng) % vmask,    &
# endif
# ifdef WET_DRY_NOT_YET
     &                     GRID(ng) % pmask_wet, GRID(ng) % pmask_full, &
     &                     GRID(ng) % rmask_wet, GRID(ng) % rmask_full, &
     &                     GRID(ng) % umask_wet, GRID(ng) % umask_full, &
     &                     GRID(ng) % vmask_wet, GRID(ng) % vmask_full, &
#  ifdef SOLVE3D
     &                     GRID(ng) % rmask_wet_avg,                    &
#  endif
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
# ifndef SOLVE3D
     &                     FORCES(ng) % ad_sustr,                       &
     &                     FORCES(ng) % ad_svstr,                       &
     &                     FORCES(ng) % ad_bustr,                       &
     &                     FORCES(ng) % ad_bvstr,                       &
#  ifdef ATM_PRESS
     &                     FORCES(ng) % Pair,                           &
#  endif
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
# ifdef MASKING
     &                           pmask, rmask, umask, vmask,            &
# endif
# ifdef WET_DRY_NOT_YET
     &                           pmask_wet, pmask_full,                 &
     &                           rmask_wet, rmask_full,                 &
     &                           umask_wet, umask_full,                 &
     &                           vmask_wet, vmask_full,                 &
#  ifdef SOLVE3D
     &                           rmask_wet_avg,                         &
#  endif
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
# ifndef SOLVE3D
     &                           ad_sustr, ad_svstr,                    &
     &                           ad_bustr, ad_bvstr,                    &
#  ifdef ATM_PRESS
     &                           Pair,                                  &
#  endif
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
      USE mod_clima
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
      USE mod_sources
!
      USE ad_exchange_2d_mod
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE obc_volcons_mod
      USE ad_obc_volcons_mod
      USE ad_u2dbc_mod, ONLY : ad_u2dbc_tile
      USE ad_v2dbc_mod, ONLY : ad_v2dbc_tile
      USE ad_zetabc_mod, ONLY : ad_zetabc_tile
# ifdef WET_DRY_NOT_YET
!>    USE wetdry_mod, ONLY : wetdry_tile
# endif
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
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
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
      real(r8), intent(in) :: rubar(LBi:,LBj:,:)
      real(r8), intent(in) :: rvbar(LBi:,LBj:,:)
      real(r8), intent(in) :: rzeta(LBi:,LBj:,:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  if !defined SOLVE3D && defined ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:,LBj:)
#  endif
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
      real(r8), intent(inout) :: pmask_full(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_full(LBi:,LBj:)
      real(r8), intent(inout) :: umask_full(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_full(LBi:,LBj:)

      real(r8), intent(inout) :: pmask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: umask_wet(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_wet(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: rmask_wet_avg(LBi:,LBj:)
#   endif
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

#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
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
      real(r8), intent(in) :: rubar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: rvbar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: rzeta(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  if !defined SOLVE3D && defined ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
#  endif
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
      real(r8), intent(inout) :: pmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_full(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: pmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_wet(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(inout) :: rmask_wet_avg(LBi:UBi,LBj:UBj)
#   endif
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

      integer :: i, is, j, ptsk
# ifdef DIAGNOSTICS_UV
!!    integer :: idiag
# endif

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6, cff7
      real(r8) :: fac, fac1, fac2, fac3
      real(r8) :: ad_cff, ad_cff1, ad_cff2, ad_cff3, ad_cff4
      real(r8) :: ad_fac, ad_fac1
      real(r8) :: adfac, adfac1, adfac2, adfac3, adfac4

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
!>    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wetdry
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
# ifdef WET_DRY_NOT_YET
!>    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_wetdry
# endif

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
      ad_cff2=0.0_r8
      ad_cff3=0.0_r8
      ad_cff4=0.0_r8
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

!  In distributed-memory, the I- and J-ranges are different and a
!  special exchange is done to avoid having three ghost points for
!  high order numerical stencils. Notice that a private array is
!  passed below to the exchange routine. It also applies periodic
!  boundary conditions, if appropriate and no partitions in I- or
!  J-directions.
!
      DO j=JstrV-2,Jendp2
        DO i=IstrU-2,Iendp2
          Dnew(i,j)=zeta(i,j,knew)+h(i,j)
          Drhs(i,j)=zeta(i,j,krhs)+h(i,j)
          Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
        END DO
      END DO
      DO j=JstrV-2,Jendp2
        DO i=IstrU-1,Iendp2
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
          DUon(i,j)=ubar(i,j,krhs)*cff1
#  ifdef NEARSHORE_MELLOR
          DUSon(i,j)=ubar_stokes(i,j)*cff1
          DUon(i,j)=DUon(i,j)+DUSon(i,j)
#  endif
        END DO
      END DO
      DO j=JstrV-1,Jendp2
        DO i=IstrU-2,Iendp2
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
          DVom(i,j)=vbar(i,j,krhs)*cff1
#  ifdef NEARSHORE_MELLOR
          DVSom(i,j)=vbar_stokes(i,j)*cff1
          DVom(i,j)=DVom(i,j)+DVSom(i,j)
#  endif
        END DO
      END DO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          DUon)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          DVom)
      END IF
      CALL mp_exchange2d (ng, tile, iADM, 2,                            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    DUon, DVom)

# else

      DO j=JstrVm2-1,Jendp2
        DO i=IstrUm2-1,Iendp2
          Dnew(i,j)=zeta(i,j,knew)+h(i,j)
          Drhs(i,j)=zeta(i,j,krhs)+h(i,j)
          Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
        END DO
      END DO
      DO j=JstrVm2-1,Jendp2
        DO i=IstrUm2,Iendp2
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
          DUon(i,j)=ubar(i,j,krhs)*cff1
#  ifdef NEARSHORE_MELLOR
          DUSon(i,j)=ubar_stokes(i,j)*cff1
          DUon(i,j)=DUon(i,j)+DUSon(i,j)
#  endif
        END DO
      END DO
      DO j=JstrVm2,Jendp2
        DO i=IstrUm2-1,Iendp2
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
          DVom(i,j)=vbar(i,j,krhs)*cff1
#  ifdef NEARSHORE_MELLOR
          DVSom(i,j)=vbar_stokes(i,j)*cff1
          DVom(i,j)=DVom(i,j)+DVSom(i,j)
#  endif
        END DO
      END DO
# endif
!
!  Compute integral mass flux across open boundaries and adjust
!  for volume conservation. Compute BASIC STATE value.
!  This must be computed here instead of below.
!
      IF (ANY(ad_VolCons(:,ng))) THEN
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
        CALL set_DUV_bc_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        krhs,                                     &
#  ifdef MASKING
     &                        umask, vmask,                             &
#  endif
     &                        om_v, on_u,                               &
     &                        ubar, vbar,                               &
     &                        Drhs, DUon, DVom)
      END IF
# if defined UV_VIS2 || defined UV_VIS4
!
!  Compute BASIC state depths at PSI-points for viscosity.
!
#  ifdef UV_VIS4
      DO j=Jstrm1,Jendp2
        DO i=Istrm1,Iendp2
#  else
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
#  endif
          Drhs_p(i,j)=0.25_r8*(Drhs(i,j  )+Drhs(i-1,j  )+               &
     &                         Drhs(i,j-1)+Drhs(i-1,j-1))
        END DO
      END DO
# endif
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
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>      CALL mp_exchange2d (ng, tile, iTLM, 2,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      tl_ubar(:,:,knew),                          &
!>   &                      tl_vbar(:,:,knew))
!>
        CALL ad_mp_exchange2d (ng, tile, iADM, 2,                       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_ubar(:,:,knew),                       &
     &                         ad_vbar(:,:,knew))
!
# endif

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>        CALL exchange_v2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_vbar(:,:,knew))
!>
          CALL ad_exchange_v2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_vbar(:,:,knew))
!>        CALL exchange_u2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_ubar(:,:,knew))
!>
          CALL ad_exchange_u2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_ubar(:,:,knew))
        END IF
!
!-----------------------------------------------------------------------
!  Apply adjoint momentum transport point sources (like river runoff),
!  if any.
!-----------------------------------------------------------------------
!
        IF (LuvSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            i=SOURCES(ng)%Isrc(is)
            j=SOURCES(ng)%Jsrc(is)
            IF (((IstrR.le.i).and.(i.le.IendR)).and.                    &
     &          ((JstrR.le.j).and.(j.le.JendR))) THEN
              IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
                cff=1.0_r8/(on_u(i,j)*                                  &
     &                      0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+          &
     &                              zeta(i  ,j,knew)+h(i  ,j)))
!>              tl_ubar(i,j,knew)=SOURCES(ng)%Qbar(is)*tl_cff
!>
                ad_cff=ad_cff+SOURCES(ng)%Qbar(is)*ad_ubar(i,j,knew)
                ad_ubar(i,j,knew)=0.0_r8
!>              tl_cff=-cff*cff*on_u(i,j)*                              &
!>   &                 0.5_r8*(tl_zeta(i-1,j,knew)+tl_h(i-1,j)+         &
!>   &                         tl_zeta(i  ,j,knew)+tl_h(i  ,j))
!>
                adfac=-cff*cff*on_u(i,j)*0.5_r8*ad_cff
                ad_h(i-1,j)=ad_h(i-1,j)+adfac
                ad_h(i  ,j)=ad_h(i  ,j)+adfac
                ad_zeta(i-1,j,knew)=ad_zeta(i-1,j,knew)+adfac
                ad_zeta(i  ,j,knew)=ad_zeta(i  ,j,knew)+adfac
                ad_cff=0.0_r8
              ELSE
                cff=1.0_r8/(om_v(i,j)*                                  &
     &                      0.5_r8*(zeta(i,j-1,knew)+h(i,j-1)+          &
     &                              zeta(i,j  ,knew)+h(i,j  )))
!>              tl_vbar(i,j,knew)=SOURCES(ng)%Qbar(is)*tl_cff
!>
                ad_cff=ad_cff+SOURCES(ng)%Qbar(is)*ad_vbar(i,j,knew)
                ad_vbar(i,j,knew)=0.0_r8
!>              tl_cff=-cff*cff*om_v(i,j)*                              &
!>   &                 0.5_r8*(tl_zeta(i,j-1,knew)+tl_h(i,j-1)+         &
!>   &                         tl_zeta(i,j  ,knew)+tl_h(i,j  ))
!>
                adfac=-cff*cff*om_v(i,j)*0.5_r8*ad_cff
                ad_h(i,j-1)=ad_h(i,j-1)+adfac
                ad_h(i,j  )=ad_h(i,j  )+adfac
                ad_zeta(i,j-1,knew)=ad_zeta(i,j-1,knew)+adfac
                ad_zeta(i,j  ,knew)=ad_zeta(i,j  ,knew)+adfac
                ad_cff=0.0_r8
              END IF
            END IF
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Apply adjoint lateral boundary conditions.
!-----------------------------------------------------------------------
!
!  Compute integral mass flux across open boundaries and adjust
!  for adjoint volume conservation.
!
        IF (ANY(ad_VolCons(:,ng))) THEN
!>        CALL tl_obc_flux_tile (ng, tile,                              &
!>   &                           LBi, UBi, LBj, UBj,                    &
!>   &                           IminS, ImaxS, JminS, JmaxS,            &
!>   &                           knew,                                  &
#  ifdef MASKING
!>   &                           umask, vmask,                          &
#  endif
!>   &                           h, tl_h, om_v, on_u,                   &
!>   &                           ubar, vbar, zeta,                      &
!>   &                           tl_ubar, tl_vbar, tl_zeta)
!>
          CALL ad_obc_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           knew,                                  &
#  ifdef MASKING
     &                           umask, vmask,                          &
#  endif
     &                           h, ad_h, om_v, on_u,                   &
     &                           ubar, vbar, zeta,                      &
     &                           ad_ubar, ad_vbar, ad_zeta)
        END IF
!
!  Adjoint lateral boundary conditons.
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
!!              DiaV2wrk(i,j,idiag)=DiaV2int(i,j,idiag)*                &
!!   &                              (pn(i,j)+pn(i,j-1))*fac
!!              DiaV2int(i,j,idiag)=cff1*DiaV2rhs(i,j,idiag)
!!            END DO
!!          END DO
!!          DO j=Jstr,Jend
!!            DO i=IstrU,Iend
!!              DiaU2wrk(i,j,idiag)=DiaU2int(i,j,idiag)*                &
!!   &                              (pm(i-1,j)+pm(i,j))*fac
!!              DiaU2int(i,j,idiag)=cff1*DiaU2rhs(i,j,idiag)
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
!!              DiaV2wrk(i,j,idiag)=DiaV2wrk(i,j,idiag)+                &
!!   &                              DiaV2int(i,j,idiag)*                &
!!   &                              (pn(i,j)+pn(i,j-1))*fac
!!              DiaV2int(i,j,idiag)=DiaV2int(i,j,idiag)+                &
!!   &                              (cff1*DiaV2rhs(i,j,idiag)+          &
!!   &                               cff2*DiaRVbar(i,j,kstp,idiag)-     &
!!   &                               cff3*DiaRVbar(i,j,ptsk,idiag))
!!            END DO
!!          END DO
!!          DO j=Jstr,Jend
!!            DO i=IstrU,Iend
!!              DiaU2wrk(i,j,idiag)=DiaU2wrk(i,j,idiag)+                &
!!   &                              DiaU2int(i,j,idiag)*                &
!!   &                              (pm(i-1,j)+pm(i,j))*fac
!!              DiaU2int(i,j,idiag)=DiaU2int(i,j,idiag)+                &
!!   &                              (cff1*DiaU2rhs(i,j,idiag)+          &
!!   &                               cff2*DiaRUbar(i,j,kstp,idiag)-     &
!!   &                               cff3*DiaRUbar(i,j,ptsk,idiag))
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
!!            DiaV2wrk(i,j,M2rate)=vbar(i,j,knew)-vbar(i,j,kstp)*       &
!!   &                             (Dstp(i,j)+Dstp(i,j-1))*fac
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!          END DO
!!        END DO
!!        DO j=Jstr,Jend
!!          DO i=IstrU,Iend
!!            DiaU2wrk(i,j,M2rate)=ubar(i,j,knew)-ubar(i,j,kstp)*       &
!!   &                             (Dstp(i,j)+Dstp(i-1,j))*fac
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!          END DO
!!        END DO
!!        DO idiag=1,NDM2d-1
!!          DO j=JstrV,Jend
!!            DO i=Istr,Iend
!!              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!!              DiaV2wrk(i,j,idiag)=cff*cff1*DiaV2rhs(i,j,idiag)*fac
!!              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!            END DO
!!          END DO
!!          DO j=Jstr,Jend
!!            DO i=IstrU,Iend
!!              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
!!              DiaU2wrk(i,j,idiag)=cff*cff1*DiaU2rhs(i,j,idiag)*fac
!!              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!            END DO
!!          END DO
!!        END DO
!!      ELSE IF (CORRECTOR_2D_STEP) THEN
!!        cff1=0.5_r8*dtfast(ng)*5.0_r8/12.0_r8
!!        cff2=0.5_r8*dtfast(ng)*8.0_r8/12.0_r8
!!        cff3=0.5_r8*dtfast(ng)*1.0_r8/12.0_r8
!!        DO j=JstrV,Jend
!!          DO i=Istr,Iend
!!            DiaV2wrk(i,j,M2rate)=vbar(i,j,knew)-                      &
!!   &                             vbar(i,j,kstp)*                      &
!!   &                             (Dstp(i,j)+Dstp(i,j-1))*fac
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!          END DO
!!        END DO
!!        DO j=Jstr,Jend
!!          DO i=IstrU,Iend
!!            DiaU2wrk(i,j,M2rate)=ubar(i,j,knew)-                      &
!!   &                             ubar(i,j,kstp)*                      &
!!   &                             (Dstp(i,j)+Dstp(i-1,j))*fac
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!          END DO
!!        END DO
!!        DO idiag=1,NDM2d-1
!!          DO j=JstrV,Jend
!!            DO i=Istr,Iend
!!              cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!!              DiaV2wrk(i,j,idiag)=cff*(cff1*DiaV2rhs(i,j,idiag)+      &
!!   &                                   cff2*DiaRVbar(i,j,kstp,idiag)- &
!!   &                                   cff3*DiaRVbar(i,j,ptsk,idiag))*&
!!   &                                  fac
!!              fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!            END DO
!!          END DO
!!          DO j=Jstr,Jend
!!            DO i=IstrU,Iend
!!              cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
!!              DiaU2wrk(i,j,idiag)=cff*(cff1*DiaU2rhs(i,j,idiag)+      &
!!   &                                   cff2*DiaRUbar(i,j,kstp,idiag)- &
!!   &                                   cff3*DiaRUbar(i,j,ptsk,idiag))*&
!!   &                                  fac
!!              fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!            END DO
!!          END DO
!!        END DO
!!      END IF
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
!>
!>            cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
!>            cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
!>            cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>            vbar(i,j,knew)=vbar(i,j,knew)*cff7
!>
!>  HGA: ADM code needed here for the above NLM code.
!>
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
!>
!>            cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
!>            cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
!>            cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>            ubar(i,j,knew)=ubar(i,j,knew)*cff7
!>
!>  HGA: ADM code needed here for the above NLM code.
!>
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
!>
!>            cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
!>            cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
!>            cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>            vbar(i,j,knew)=vbar(i,j,knew)*cff7
!>
!>  HGA: ADM code needed here for the above NLM code.
!>
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
!>
!>            cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
!>            cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
!>            cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>            ubar(i,j,knew)=ubar(i,j,knew)*cff7
!>
!>  HGA: ADM code needed here for the above NLM code.
!>
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
!>
!>            cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
!>            cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
!>            cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>            vbar(i,j,knew)=vbar(i,j,knew)*cff7
!>
!>  HGA: ADM code needed here for the above NLM code.
!>
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
!>
!>            cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
!>            cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
!>            cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>            ubar(i,j,knew)=ubar(i,j,knew)*cff7
!>
!>  HGA: ADM code needed here for the above NLM code.
!>
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
!!            DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaV2rhs(i,j,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DiaV2rhs(i,j,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DO idiag=1,M2pgrd
!!              DiaRVfrc(i,j,nstp,idiag)=DiaRVfrc(i,j,3,idiag)
!!              DiaV2rhs(i,j,idiag)=DiaV2rhs(i,j,idiag)+                &
!!   &                              DiaRVfrc(i,j,3,idiag)
!!              DiaRVfrc(i,j,3,idiag)=DiaRVfrc(i,j,3,idiag)-            &
!!   &                                DiaV2rhs(i,j,idiag)
!!            END DO
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
!!            DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaU2rhs(i,j,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DiaU2rhs(i,j,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DO idiag=1,M2pgrd
!!              DiaRUfrc(i,j,nstp,idiag)=DiaRUfrc(i,j,3,idiag)
!!              DiaU2rhs(i,j,idiag)=DiaU2rhs(i,j,idiag)+                &
!!   &                              DiaRUfrc(i,j,3,idiag)
!!              DiaRUfrc(i,j,3,idiag)=DiaRUfrc(i,j,3,idiag)-            &
!!   &                                DiaU2rhs(i,j,idiag)
!!            END DO
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
!!            DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaV2rhs(i,j,M2bstr)=1.5_r8*DiaRVfrc(i,j,3,M2bstr)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2bstr)
!!            DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DiaV2rhs(i,j,M2sstr)=1.5_r8*DiaRVfrc(i,j,3,M2sstr)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2sstr)
!!            DO idiag=1,M2pgrd
!!              DiaRVfrc(i,j,nstp,idiag)=DiaRVfrc(i,j,3,idiag)
!!              DiaV2rhs(i,j,idiag)=DiaV2rhs(i,j,idiag)+                &
!!   &                              1.5_r8*DiaRVfrc(i,j,3,idiag)-       &
!!   &                              0.5_r8*DiaRVfrc(i,j,nnew,idiag)
!!              DiaRVfrc(i,j,3,idiag)=DiaRVfrc(i,j,3,idiag)-            &
!!   &                                DiaV2rhs(i,j,idiag)
!!            END DO
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
!!            DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaU2rhs(i,j,M2bstr)=1.5_r8*DiaRUfrc(i,j,3,M2bstr)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2bstr)
!!            DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DiaU2rhs(i,j,M2sstr)=1.5_r8*DiaRUfrc(i,j,3,M2sstr)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2sstr)
!!            DO idiag=1,M2pgrd
!!              DiaRUfrc(i,j,nstp,idiag)=DiaRUfrc(i,j,3,idiag)
!!              DiaU2rhs(i,j,idiag)=DiaU2rhs(i,j,idiag)+                &
!!   &                              1.5_r8*DiaRUfrc(i,j,3,idiag)-       &
!!   &                              0.5_r8*DiaRUfrc(i,j,nnew,idiag)
!!              DiaRUfrc(i,j,3,idiag)=DiaRUfrc(i,j,3,idiag)-            &
!!   &                                DiaU2rhs(i,j,idiag)
!!            END DO
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
!!            DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaV2rhs(i,j,M2bstr)=cff1*DiaRVfrc(i,j,3,M2bstr)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2bstr)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2bstr)
!!            DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DiaV2rhs(i,j,M2sstr)=cff1*DiaRVfrc(i,j,3,M2sstr)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2sstr)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2sstr)
!!            DO idiag=1,M2pgrd
!!              DiaRVfrc(i,j,nstp,idiag)=DiaRVfrc(i,j,3,idiag)
!!              DiaV2rhs(i,j,idiag)=DiaV2rhs(i,j,idiag)+                &
!!   &                              cff1*DiaRVfrc(i,j,3,idiag)-         &
!!   &                              cff2*DiaRVfrc(i,j,nnew,idiag)+      &
!!   &                              cff3*DiaRVfrc(i,j,nstp,idiag)
!!              DiaRVfrc(i,j,3,idiag)=DiaRVfrc(i,j,3,idiag)-            &
!!   &                                DiaV2rhs(i,j,idiag)
!!            END DO
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
!!            DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaU2rhs(i,j,M2bstr)=cff1*DiaRUfrc(i,j,3,M2bstr)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2bstr)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2bstr)
!!            DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DiaU2rhs(i,j,M2sstr)=cff1*DiaRUfrc(i,j,3,M2sstr)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2sstr)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2sstr)
!!            DO idiag=1,M2pgrd
!!              DiaRUfrc(i,j,nstp,idiag)=DiaRUfrc(i,j,3,idiag)
!!              DiaU2rhs(i,j,idiag)=DiaU2rhs(i,j,idiag)+                &
!!   &                              cff1*DiaRUfrc(i,j,3,idiag)-         &
!!   &                              cff2*DiaRUfrc(i,j,nnew,idiag)+      &
!!   &                              cff3*DiaRUfrc(i,j,nstp,idiag)
!!              DiaRUfrc(i,j,3,idiag)=DiaRUfrc(i,j,3,idiag)-            &
!!   &                                DiaU2rhs(i,j,idiag)
!!            END DO
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
!!          DiaV2rhs(i,j,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!          DiaV2rhs(i,j,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!          DO idiag=1,M2pgrd
!!            DiaV2rhs(i,j,idiag)=DiaV2rhs(i,j,idiag)+                  &
!!   &                            DiaRVfrc(i,j,3,idiag)
!!          END DO
#  endif
!>            tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
!>
              ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rhs_vbar(i,j)
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!          DiaU2rhs(i,j,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!          DO idiag=1,M2pgrd
!!            DiaU2rhs(i,j,idiag)=DiaU2rhs(i,j,idiag)+                  &
!!   &                            DiaRUfrc(i,j,3,idiag)
!!          END DO
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
!!          DiaU2rhs(i,j,M2sstr)=fac
#  endif
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
!>
            ad_fac=ad_fac+ad_rhs_ubar(i,j)
!>          tl_fac=tl_sustr(i,j)*om_u(i,j)*on_u(i,j)
!>
            ad_sustr(i,j)=ad_sustr(i,j)+om_u(i,j)*on_u(i,j)*ad_fac
            ad_fac=0.0_r8
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2sstr)=fac
#  endif
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
!>
            ad_fac=ad_fac+ad_rhs_vbar(i,j)
!>          tl_fac=tl_svstr(i,j)*om_v(i,j)*on_v(i,j)
!>
            ad_svstr(i,j)=ad_svstr(i,j)+om_v(i,j)*on_v(i,j)*ad_fac
            ad_fac=0.0_r8
          END DO
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Add in adjoint nudging of 2D momentum climatology.
!-----------------------------------------------------------------------
!
        IF (LnudgeM2CLM(ng)) THEN
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=0.25_r8*(CLIMA(ng)%M2nudgcof(i,j-1)+                  &
     &                     CLIMA(ng)%M2nudgcof(i,j  ))*                 &
     &            om_v(i,j)*on_v(i,j)
!>            tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                        &
!>   &                         cff*((Drhs(i,j-1)+Drhs(i,j))*            &
!>   &                              (-tl_vbar(i,j,krhs))+               &
!>   &                              (tl_Drhs(i,j-1)+tl_Drhs(i,j))*      &
!>   &                              (CLIMA(ng)%vbarclm(i,j)-
!>   &                               vbar(i,j,krhs)))
!>
              adfac=cff*ad_rhs_vbar(i,j)
              adfac1=adfac*(Drhs(i,j-1)+Drhs(i,j))
              adfac2=adfac*(CLIMA(ng)%vbarclm(i,j)-vbar(i,j,krhs))
              ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)-adfac1
              ad_Drhs(i,j-1)=ad_Drhs(i,j-1)+adfac2
              ad_Drhs(i,j  )=ad_Drhs(i,j  )+adfac2
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=0.25_r8*(CLIMA(ng)%M2nudgcof(i-1,j)+                  &
     &                     CLIMA(ng)%M2nudgcof(i  ,j))*                 &
     &            om_u(i,j)*on_u(i,j)
!>            tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                        &
!>   &                         cff*((Drhs(i-1,j)+Drhs(i,j))*            &
!>   &                              (-tl_ubar(i,j,krhs))+               &
!>   &                              (tl_Drhs(i-1,j)+tl_Drhs(i,j))*      &
!>   &                              (CLIMA(ng)%ubarclm(i,j)-
!>   &                               ubar(i,j,krhs)))
!>
              adfac=cff*ad_rhs_ubar(i,j)
              adfac1=adfac*(Drhs(i-1,j)+Drhs(i,j))
              adfac2=adfac*(CLIMA(ng)%ubarclm(i,j)-ubar(i,j,krhs))
              ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)-adfac1
              ad_Drhs(i-1,j)=ad_Drhs(i-1,j)+adfac2
              ad_Drhs(i  ,j)=ad_Drhs(i  ,j)+adfac2
            END DO
          END DO
        END IF

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
!!          DiaV2rhs(i,j,M2hrad)=-cff1
#  endif
#  ifndef SOLVE3D
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_cff1-tl_cff2
!>
            ad_cff2=ad_cff2-ad_rhs_vbar(i,j)
            ad_cff1=ad_cff1-ad_rhs_vbar(i,j)
#  endif
!>          tl_cff2=tl_rvlag2d(i,j)
!>
            ad_rvlag2d(i,j)=ad_rvlag2d(i,j)+ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=tl_rvstr2d(i,j)*om_v(i,j)*on_v(i,j)
!>
            ad_rvstr2d(i,j)=ad_rvstr2d(i,j)+                            &
     &                      om_v(i,j)*on_v(i,j)*ad_cff1
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2hrad)=-cff1
#  endif
#  ifndef SOLVE3D
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_cff1-tl_cff2
!>
            ad_cff2=ad_cff2-ad_rhs_ubar(i,j)
            ad_cff1=ad_cff1-ad_rhs_ubar(i,j)
#  endif
!>          tl_cff2=tl_rulag2d(i,j)
!>
            ad_rulag2d(i,j)=ad_rulag2d(i,j)+ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=tl_rustr2d(i,j)*om_u(i,j)*on_u(i,j)
!>
            ad_rustr2d(i,j)=ad_rustr2d(i,j)+                            &
     &                      om_u(i,j)*on_u(i,j)*ad_cff1
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
        DO j=Jstrm1,Jendp2
          DO i=Istrm1,Iendp2
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
        DO j=JstrVm2,Jendp1
          DO i=IstrUm2,Iendp1
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
        DO j=Jstrm1,Jendp2
          DO i=Istrm1,Iendp2
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
        DO j=Jstrm1,Jendp1
          DO i=IstrUm1,Iendp1
            LapU(i,j)=0.125_r8*                                         &
     &                (pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*          &
     &                ((pn(i-1,j)+pn(i,j))*                             &
     &                 (UFx(i,j  )-UFx(i-1,j))+                         &
     &                 (pm(i-1,j)+pm(i,j))*                             &
     &                 (UFe(i,j+1)-UFe(i  ,j)))
          END DO
        END DO
        DO j=JstrVm1,Jendp1
          DO i=Istrm1,Iendp1
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
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (ad_LBC(iwest,isUbar,ng)%closed) THEN
              DO j=Jstrm1,Jendp1
                LapU(IstrU-1,j)=0.0_r8
              END DO
            ELSE
              DO j=Jstrm1,Jendp1
                LapU(IstrU-1,j)=LapU(IstrU,j)
              END DO
            END IF
            IF (ad_LBC(iwest,isVbar,ng)%closed) THEN
              DO j=JstrVm1,Jendp1
                LapV(Istr-1,j)=gamma2(ng)*LapV(Istr,j)
              END DO
            ELSE
              DO j=JstrVm1,Jendp1
                LapV(Istr-1,j)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (ad_LBC(ieast,isUbar,ng)%closed) THEN
              DO j=Jstrm1,Jendp1
                LapU(Iend+1,j)=0.0_r8
              END DO
            ELSE
              DO j=Jstrm1,Jendp1
                LapU(Iend+1,j)=LapU(Iend,j)
              END DO
            END IF
            IF (ad_LBC(ieast,isVbar,ng)%closed) THEN
              DO j=JstrVm1,Jendp1
                LapV(Iend+1,j)=gamma2(ng)*LapV(Iend,j)
              END DO
            ELSE
              DO j=JstrVm1,Jendp1
                LapV(Iend+1,j)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (ad_LBC(isouth,isUbar,ng)%closed) THEN
              DO i=IstrUm1,Iendp1
                LapU(i,Jstr-1)=gamma2(ng)*LapU(i,Jstr)
              END DO
            ELSE
              DO i=IstrUm1,Iendp1
                LapU(i,Jstr-1)=0.0_r8
              END DO
            END IF
            IF (ad_LBC(isouth,isVbar,ng)%closed) THEN
              DO i=Istrm1,Iendp1
                LapV(i,JstrV-1)=0.0_r8
              END DO
            ELSE
              DO i=Istrm1,Iendp1
                LapV(i,JstrV-1)=LapV(i,JstrV)
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (ad_LBC(inorth,isUbar,ng)%closed) THEN
              DO i=IstrUm1,Iendp1
                LapU(i,Jend+1)=gamma2(ng)*LapU(i,Jend)
              END DO
            ELSE
              DO i=IstrUm1,Iendp1
                LapU(i,Jend+1)=0.0_r8
              END DO
            END IF
            IF (ad_LBC(inorth,isVbar,ng)%closed) THEN
              DO i=Istrm1,Iendp1
                LapV(i,Jend+1)=0.0_r8
              END DO
            ELSE
              DO i=Istrm1,Iendp1
                LapV(i,Jend+1)=LapV(i,Jend)
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            LapU(Istr  ,Jstr-1)=0.5_r8*(LapU(Istr+1,Jstr-1)+            &
     &                                  LapU(Istr  ,Jstr  ))
            LapV(Istr-1,Jstr  )=0.5_r8*(LapV(Istr-1,Jstr+1)+            &
     &                                  LapV(Istr  ,Jstr  ))
          END IF
        END IF

        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            LapU(Iend+1,Jstr-1)=0.5_r8*(LapU(Iend  ,Jstr-1)+            &
     &                                  LapU(Iend+1,Jstr  ))
            LapV(Iend+1,Jstr  )=0.5_r8*(LapV(Iend  ,Jstr  )+            &
     &                                  LapV(Iend+1,Jstr+1))
          END IF
        END IF

        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            LapU(Istr  ,Jend+1)=0.5_r8*(LapU(Istr+1,Jend+1)+            &
     &                                  LapU(Istr  ,Jend  ))
            LapV(Istr-1,Jend+1)=0.5_r8*(LapV(Istr  ,Jend+1)+            &
     &                                  LapV(Istr-1,Jend  ))
          END IF
        END IF

        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            LapU(Iend+1,Jend+1)=0.5_r8*(LapU(Iend  ,Jend+1)+            &
     &                                  LapU(Iend+1,Jend  ))
            LapV(Iend+1,Jend+1)=0.5_r8*(LapV(Iend  ,Jend+1)+            &
     &                                  LapV(Iend+1,Jend  ))
          END IF
        END IF
!
!  Add in adjoint biharmocnic viscosity
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  if defined DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2yvis)=-cff2
!!          DiaV2rhs(i,j,M2xvis)= cff1
!!          DiaV2rhs(i,j,M2hvis)=fac
#  endif
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
!>
            ad_fac=ad_fac+ad_rhs_vbar(i,j)
!>          tl_fac=tl_cff1-tl_cff2
!>
            ad_cff1=ad_cff1+ad_fac
            ad_cff2=ad_cff2-ad_fac
            ad_fac=0.0_r8
!>          tl_cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*                         &
!>   &              (tl_VFe(i  ,j)-tl_VFe(i,j-1))
!>
            adfac=0.5_r8*(pm(i,j-1)+pm(i,j))*ad_cff2
            ad_VFe(i,j-1)=ad_VFe(i,j-1)-adfac
            ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
            ad_cff2=0.0_r8
!>          tl_cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*                         &
!>   &              (tl_VFx(i+1,j)-tl_VFx(i,j  ))
!>
            adfac=0.5_r8*(pn(i,j-1)+pn(i,j))*ad_cff1
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-adfac
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+adfac
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  if defined DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2yvis)=cff2
!!          DiaU2rhs(i,j,M2xvis)=cff1
!!          DiaU2rhs(i,j,M2hvis)=fac
#  endif
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
!>
            ad_fac=ad_fac+ad_rhs_ubar(i,j)
!>          tl_fac=tl_cff1+tl_cff2
!>
            ad_cff1=ad_cff1+ad_fac
            ad_cff2=ad_cff2+ad_fac
            ad_fac=0.0_r8
!>          tl_cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*                         &
!>   &              (UFe(i,j+1)-UFe(i  ,j))
!>
            adfac=0.5_r8*(pm(i-1,j)+pm(i,j))*ad_cff2
            UFe(i,j  )=UFe(i,j  )-adfac
            UFe(i,j+1)=UFe(i,j+1)+adfac
            ad_cff2=0.0_r8
!>          tl_cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*                         &
!>   &              (tl_UFx(i,j  )-tl_UFx(i-1,j))
!>
            adfac=0.5_r8*(pn(i-1,j)+pn(i,j))*ad_cff1
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-adfac
            ad_UFx(i,j  )=ad_UFx(i,j  )+adfac
            ad_cff1=0.0_r8
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
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
!>          tl_LapV(Iend+1,Jend+1)=0.5_r8*(tl_LapV(Iend  ,Jend+1)+      &
!>   &                                     tl_LapV(Iend+1,Jend  ))
!>
            adfac=0.5_r8*ad_LapV(Iend+1,Jend+1)
            ad_LapV(Iend+1,Jend  )=ad_LapV(Iend+1,Jend  )+adfac
            ad_LapV(Iend  ,Jend+1)=ad_LapV(Iend  ,Jend+1)+adfac
            ad_LapV(Iend+1,Jend+1)=0.0_r8
!>          tl_LapU(Iend+1,Jend+1)=0.5_r8*(tl_LapU(Iend  ,Jend+1)+      &
!>   &                                     tl_LapU(Iend+1,Jend  ))
!>
            adfac=0.5_r8*ad_LapU(Iend+1,Jend+1)
            ad_LapU(Iend+1,Jend  )=ad_LapU(Iend+1,Jend  )+adfac
            ad_LapU(Iend  ,Jend+1)=ad_LapU(Iend  ,Jend+1)+adfac
            ad_LapU(Iend+1,Jend+1)=0.0_r8
          END IF
        END IF

        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
!>          tl_LapV(Istr-1,Jend+1)=0.5_r8*(tl_LapV(Istr  ,Jend+1)+      &
!>   &                                     tl_LapV(Istr-1,Jend ))
!>
            adfac=0.5_r8*ad_LapV(Istr-1,Jend+1)
            ad_LapV(Istr-1,Jend  )=ad_LapV(Istr-1,Jend  )+adfac
            ad_LapV(Istr  ,Jend+1)=ad_LapV(Istr  ,Jend+1)+adfac
            ad_LapV(Istr-1,Jend+1)=0.0_r8
!>          tl_LapU(Istr  ,Jend+1)=0.5_r8*(tl_LapU(Istr+1,Jend+1)+      &
!>   &                                     tl_LapU(Istr  ,Jend  ))
!>
            adfac=0.5_r8*ad_LapU(Istr  ,Jend+1)
            ad_LapU(Istr  ,Jend  )=ad_LapU(Istr  ,Jend  )+adfac
            ad_LapU(Istr+1,Jend+1)=ad_LapU(Istr+1,Jend+1)+adfac
            ad_LapU(Istr  ,Jend+1)=0.0_r8
          END IF
        END IF

        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
!>          tl_LapV(Iend+1,Jstr  )=0.5_r8*(tl_LapV(Iend  ,Jstr  )+      &
!>   &                                     tl_LapV(Iend+1,Jstr+1))
!>
            adfac=0.5_r8*ad_LapV(Iend+1,Jstr  )
            ad_LapV(Iend  ,Jstr  )=ad_LapV(Iend  ,Jstr  )+adfac
            ad_LapV(Iend+1,Jstr+1)=ad_LapV(Iend+1,Jstr+1)+adfac
            ad_LapV(Iend+1,Jstr  )=0.0_r8
!>          tl_LapU(Iend+1,Jstr-1)=0.5_r8*(tl_LapU(Iend  ,Jstr-1)+      &
!>   &                                     tl_LapU(Iend+1,Jstr  ))
!>
            adfac=0.5_r8*ad_LapU(Iend+1,Jstr-1)
            ad_LapU(Iend  ,Jstr-1)=ad_LapU(Iend  ,Jstr-1)+adfac
            ad_LapU(Iend+1,Jstr  )=ad_LapU(Iend+1,Jstr  )+adfac
            ad_LapU(Iend+1,Jstr-1)=0.0_r8
          END IF
        END IF

        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
!>          tl_LapV(Istr-1,Jstr  )=0.5_r8*(tl_LapV(Istr-1,Jstr+1)+      &
!>   &                                     tl_LapV(Istr  ,Jstr  ))
!>
            adfac=0.5_r8*ad_LapV(Istr-1,Jstr  )
            ad_LapV(Istr  ,Jstr  )=ad_LapV(Istr  ,Jstr  )+adfac
            ad_LapV(Istr-1,Jstr+1)=ad_LapV(Istr-1,Jstr+1)+adfac
            ad_LapV(Istr-1,Jstr  )=0.0_r8
!>          tl_LapU(Istr  ,Jstr-1)=0.5_r8*(tl_LapU(Istr+1,Jstr-1)+      &
!>   &                                     tl_LapU(Istr  ,Jstr  ))
!>
            adfac=0.5_r8*ad_LapU(Istr  ,Jstr-1)
            ad_LapU(Istr+1,Jstr-1)=ad_LapU(Istr+1,Jstr-1)+adfac
            ad_LapU(Istr  ,Jstr  )=ad_LapU(Istr  ,Jstr  )+adfac
            ad_LapU(Istr  ,Jstr-1)=0.0_r8
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (ad_LBC(inorth,isVbar,ng)%closed) THEN
              DO i=Istrm1,Iendp1
!>              tl_LapV(i,Jend+1)=0.0_r8
!>
                ad_LapV(i,Jend+1)=0.0_r8
              END DO
            ELSE
              DO i=Istrm1,Iendp1
!>              tl_LapV(i,Jend+1)=tl_LapV(i,Jend)
!>
                ad_LapV(i,Jend)=ad_LapV(i,Jend)+ad_LapV(i,Jend+1)
                ad_LapV(i,Jend+1)=0.0_r8
              END DO
            END IF
            IF (ad_LBC(inorth,isUbar,ng)%closed) THEN
              DO i=IstrUm1,Iendp1
!>              tl_LapU(i,Jend+1)=gamma2(ng)*tl_LapU(i,Jend)
!>
                ad_LapU(i,Jend)=ad_LapU(i,Jend)+                        &
     &                          gamma2(ng)*ad_LapU(i,Jend+1)
                ad_LapU(i,Jend+1)=0.0_r8
              END DO
            ELSE
              DO i=IstrUm1,Iendp1
!>              tl_LapU(i,Jend+1)=0.0_r8
!>
                ad_LapU(i,Jend+1)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (ad_LBC(isouth,isVbar,ng)%closed) THEN
              DO i=Istrm1,Iendp1
!>              tl_LapV(i,JstrV-1)=0.0_r8
!>
                ad_LapV(i,JstrV-1)=0.0_r8
              END DO
            ELSE
              DO i=Istrm1,Iendp1
!>              tl_LapV(i,JstrV-1)=tl_LapV(i,JstrV)
!>
                ad_LapV(i,JstrV)=ad_LapV(i,JstrV)+ad_LapV(i,JstrV-1)
                ad_LapV(i,JstrV-1)=0.0_r8
              END DO
            END IF
            IF (ad_LBC(isouth,isUbar,ng)%closed) THEN
              DO i=IstrUm1,Iendp1
!>              tl_LapU(i,Jstr-1)=gamma2(ng)*tl_LapU(i,Jstr)
!>
                ad_LapU(i,Jstr)=ad_LapU(i,Jstr)+                        &
     &                          gamma2(ng)*ad_LapU(i,Jstr-1)
                ad_LapU(i,Jstr-1)=0.0_r8
              END DO
            ELSE
              DO i=IstrUm1,Iendp1
!>              tl_LapU(i,Jstr-1)=0.0_r8
!>
                ad_LapU(i,Jstr-1)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (ad_LBC(ieast,isVbar,ng)%closed) THEN
              DO j=JstrVm1,Jendp1
!>              tl_LapV(Iend+1,j)=gamma2(ng)*tl_LapV(Iend,j)
!>
                ad_LapV(Iend,j)=ad_LapV(Iend,j)+                        &
     &                          gamma2(ng)*ad_LapV(Iend+1,j)
                ad_LapV(Iend+1,j)=0.0_r8
              END DO
            ELSE
              DO j=JstrVm1,Jendp1
!>              tl_LapV(Iend+1,j)=0.0_r8
!>
                ad_LapV(Iend+1,j)=0.0_r8
              END DO
            END IF
            IF (ad_LBC(ieast,isUbar,ng)%closed) THEN
              DO j=Jstrm1,Jendp1
!>              tl_LapU(Iend+1,j)=0.0_r8
!>
                ad_LapU(Iend+1,j)=0.0_r8
              END DO
            ELSE
              DO j=Jstrm1,Jendp1
!>              tl_LapU(Iend+1,j)=tl_LapU(Iend,j)
!>
                ad_LapU(Iend,j)=ad_LapU(Iend,j)+ad_LapU(Iend+1,j)
                ad_LapU(Iend+1,j)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (ad_LBC(iwest,isVbar,ng)%closed) THEN
              DO j=JstrVm1,Jendp1
!>              tl_LapV(Istr-1,j)=gamma2(ng)*tl_LapV(Istr,j)
!>
                ad_LapV(Istr,j)=ad_LapV(Istr,j)+                        &
     &                          gamma2(ng)*ad_LapV(Istr-1,j)
                ad_LapV(Istr-1,j)=0.0_r8
              END DO
            ELSE
              DO j=JstrVm1,Jendp1
!>              tl_LapV(Istr-1,j)=0.0_r8
!>
                ad_LapV(Istr-1,j)=0.0_r8
              END DO
            END IF
            IF (ad_LBC(iwest,isUbar,ng)%closed) THEN
              DO j=Jstrm1,Jendp1
!>              tl_LapU(IstrU-1,j)=0.0_r8
!>
                ad_LapU(IstrU-1,j)=0.0_r8
              END DO
            ELSE
              DO j=Jstrm1,Jendp1
!>              tl_LapU(IstrU-1,j)=tl_LapU(IstrU,j)
!>
                ad_LapU(IstrU,j)=ad_LapU(IstrU,j)+ad_LapU(IstrU-1,j)
                ad_LapU(IstrU-1,j)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
!  Compute adjoint first harmonic operator (m s^-3/2).
!
        DO j=JstrVm1,Jendp1
          DO i=Istrm1,Iendp1
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
        DO j=Jstrm1,Jendp1
          DO i=IstrUm1,Iendp1
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
        DO j=Jstrm1,Jendp2
          DO i=Istrm1,Iendp2
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
        DO j=JstrVm2,Jendp1
          DO i=IstrUm2,Iendp1
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
!!          DiaV2rhs(i,j,M2yvis)=-cff2
!!          DiaV2rhs(i,j,M2xvis)= cff1
!!          DiaV2rhs(i,j,M2hvis)=fac
#  endif
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
!>
            ad_fac=ad_fac+ad_rhs_vbar(i,j)
!>          tl_fac=tl_cff1-tl_cff2
!>
            ad_cff1=ad_cff1+ad_fac
            ad_cff2=ad_cff2-ad_fac
            ad_fac=0.0_r8
!>          tl_cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*                         &
!>   &              (tl_VFe(i  ,j)-tl_VFe(i,j-1))
!>
            adfac=0.5_r8*(pm(i,j-1)+pm(i,j))*ad_cff2
            ad_VFe(i,j-1)=ad_VFe(i,j-1)-adfac
            ad_VFe(i  ,j)=ad_VFe(i  ,j)+adfac
            ad_cff2=0.0_r8
!>          tl_cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*                         &
!>   &              (tl_VFx(i+1,j)-tl_VFx(i,j  ))
!>
            adfac=0.5_r8*(pn(i,j-1)+pn(i,j))*ad_cff1
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-adfac
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+adfac
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  if defined DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2yvis)=cff2
!!          DiaU2rhs(i,j,M2xvis)=cff1
!!          DiaU2rhs(i,j,M2hvis)=fac
#  endif
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
!>
            ad_fac=ad_fac+ad_rhs_ubar(i,j)
!>          tl_fac=tl_cff1+tl_cff2
!>
            ad_cff1=ad_cff1+ad_fac
            ad_cff2=ad_cff2+ad_fac
            ad_fac=0.0_r8
!>          tl_cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*                         &
!>   &              (tl_UFe(i,j+1)-tl_UFe(i  ,j))
!>
            adfac=0.5_r8*(pm(i-1,j)+pm(i,j))*ad_cff2
            ad_UFe(i,j  )=ad_UFe(i,j  )-adfac
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+adfac
            ad_cff2=0.0_r8
!>          tl_cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*                         &
!>   &              (tl_UFx(i,j  )-tl_UFx(i-1,j))
!>
            adfac=0.5_r8*(pn(i-1,j)+pn(i,j))*ad_cff1
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-adfac
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
            ad_cff1=0.0_r8
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
        DO j=Jstrm1,Jendp2
          DO i=Istrm1,Iendp2
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
# endif
# if defined CURVGRID && defined UV_ADV
!
!-----------------------------------------------------------------------
!  Add in curvilinear transformation terms.
!-----------------------------------------------------------------------
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
#  if defined DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)-fac1
!!        DiaV2rhs(i,j,M2yadv)=DiaV2rhs(i,j,M2yadv)-fac2
!!        DiaV2rhs(i,j,M2xadv)=DiaV2rhs(i,j,M2xadv)-fac1+fac2
!!        fac2=0.5_r8*(Vwrk(i,j)+Vwrk(i,j-1))
#  endif
!>        tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac1
!>
          ad_fac1=ad_fac1-ad_rhs_vbar(i,j)
!>        tl_fac1=0.5_r8*(tl_VFe(i,j)+tl_VFe(i,j-1))
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
!!        DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+fac1
!!        DiaU2rhs(i,j,M2yadv)=DiaU2rhs(i,j,M2yadv)+fac2
!!        DiaU2rhs(i,j,M2xadv)=DiaU2rhs(i,j,M2xadv)+fac1-fac2
!!        fac2=0.5_r8*(Uwrk(i,j)+Uwrk(i-1,j))
#  endif
!>        tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac1
!>
          ad_fac1=ad_fac1+ad_rhs_ubar(i,j)
!>        tl_fac1=0.5_r8*(tl_UFx(i,j)+tl_UFx(i-1,j))
!>
          adfac=0.5_r8*ad_fac1
          ad_UFx(i-1,j)=ad_UFx(i-1,j)+adfac
          ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
          ad_fac1=0.0_r8
        END DO
      END DO
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=0.5_r8*(vbar(i,j  ,krhs)+                                &
#  ifdef NEARSHORE_MELLOR
     &                 vbar_stokes(i,j  )+                              &
     &                 vbar_stokes(i,j+1)+                              &
#  endif
     &                 vbar(i,j+1,krhs))
          cff2=0.5_r8*(ubar(i  ,j,krhs)+                                &
#  ifdef NEARSHORE_MELLOR
     &                 ubar_stokes(i  ,j)+                              &
     &                 ubar_stokes(i+1,j)+                              &
#  endif
     &                 ubar(i+1,j,krhs))
          cff3=cff1*dndx(i,j)
          cff4=cff2*dmde(i,j)
          cff=Drhs(i,j)*(cff3-cff4)
#  if defined DIAGNOSTICS_UV
!!        Vwrk(i,j)=-cff*cff2                  ! vbar equation, ETA-term
!!        Uwrk(i,j)=-cff*cff1                  ! ubar equation, ETA-term
!!        cff=Drhs(i,j)*cff4
#  endif
!>        tl_VFe(i,j)=tl_cff*cff2+cff*tl_cff2
!>        tl_UFx(i,j)=tl_cff*cff1+cff*tl_cff1
!>
          ad_cff=ad_cff+                                                &
     &           cff1*ad_UFx(i,j)+                                      &
     &           cff2*ad_VFe(i,j)
          ad_cff1=ad_cff1+cff*ad_UFx(i,j)
          ad_cff2=ad_cff2+cff*ad_VFe(i,j)
          ad_UFx(i,j)=0.0_r8
          ad_VFe(i,j)=0.0_r8
!>        tl_cff=tl_Drhs(i,j)*(cff3-cff4)+                              &
!>   &           Drhs(i,j)*(tl_cff3-tl_cff4)
!>
          adfac=Drhs(i,j)*ad_cff
          ad_cff4=ad_cff4-adfac
          ad_cff3=ad_cff3+adfac
          ad_Drhs(i,j)=ad_Drhs(i,j)+(cff3-cff4)*ad_cff
          ad_cff=0.0_r8
!>        tl_cff4=tl_cff2*dmde(i,j)
!>
          ad_cff2=ad_cff2+dmde(i,j)*ad_cff4
          ad_cff4=0.0_r8
!>        tl_cff3=tl_cff1*dndx(i,j)
!>
          ad_cff1=ad_cff1+dndx(i,j)*ad_cff3
          ad_cff3=0.0_r8
!>        tl_cff2=0.5_r8*(tl_ubar(i  ,j,krhs)+                          &
#  ifdef NEARSHORE_MELLOR
!>   &                    tl_ubar_stokes(i  ,j)+                        &
!>   &                    tl_ubar_stokes(i+1,j)+                        &
#  endif
!>   &                    tl_ubar(i+1,j,krhs))
!>
          adfac=0.5_r8*ad_cff2
          ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)+adfac
          ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+adfac
#  ifdef NEARSHORE_MELLOR
          ad_ubar_stokes(i  ,j)=ad_ubar_stokes(i  ,j)+adfac
          ad_ubar_stokes(i+1,j)=ad_ubar_stokes(i+1,j)+adfac
#  endif
          ad_cff2=0.0_r8
!>        tl_cff1=0.5_r8*(tl_vbar(i,j  ,krhs)+                          &
#  ifdef NEARSHORE_MELLOR
!>   &                    tl_vbar_stokes(i,j  )+                        &
!>   &                    tl_vbar_stokes(i,j+1)+                        &
#  endif
!>   &                    tl_vbar(i,j+1,krhs))
!>
          adfac=0.5_r8*ad_cff1
          ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+adfac
          ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+adfac
#  ifdef NEARSHORE_MELLOR
          ad_vbar_stokes(i,j  )=ad_vbar_stokes(i,j  )+adfac
          ad_vbar_stokes(i,j+1)=ad_vbar_stokes(i,j+1)+adfac
#  endif
          ad_cff1=0.0_r8
        END DO
      END DO
# endif
# ifdef UV_COR
!
!-----------------------------------------------------------------------
!  Add in Coriolis term.
!-----------------------------------------------------------------------
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
#  if defined DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2fcor)=-fac1
#  endif
!>        tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac1
!>
          ad_fac1=ad_fac1-ad_rhs_vbar(i,j)
!>        tl_fac1=0.5_r8*(tl_VFe(i,j)+tl_VFe(i,j-1))
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
!!        DiaU2rhs(i,j,M2fcor)=fac1
#  endif
!>        tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac1
!>
          ad_fac1=ad_fac1+ad_rhs_ubar(i,j)
!>        tl_fac1=0.5_r8*(tl_UFx(i,j)+tl_UFx(i-1,j))
!>
          adfac=0.5_r8*ad_fac1
          ad_UFx(i-1,j)=ad_UFx(i-1,j)+adfac
          ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
          ad_fac1=0.0_r8
        END DO
      END DO
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff=0.5_r8*Drhs(i,j)*fomn(i,j)
!>        tl_VFe(i,j)=tl_cff*(ubar(i  ,j,krhs)+                         &
#  ifdef NEARSHORE_MELLOR
!>   &                        ubar_stokes(i  ,j)+                       &
!>   &                        ubar_stokes(i+1,j)+                       &
#  endif
!>   &                        ubar(i+1,j,krhs))+                        &
!>   &                cff*(tl_ubar(i  ,j,krhs)+                         &
#  ifdef NEARSHORE_MELLOR
!>   &                     tl_ubar_stokes(i  ,j)+                       &
!>   &                     tl_ubar_stokes(i+1,j)+                       &
#  endif
!>   &                     tl_ubar(i+1,j,krhs))
!>
          adfac=cff*ad_VFe(i,j)
          ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)+adfac
          ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+adfac
#  ifdef NEARSHORE_MELLOR
          ad_ubar_stokes(i  ,j)=ad_ubar_stokes(i  ,j)+adfac
          ad_ubar_stokes(i+1,j)=ad_ubar_stokes(i+1,j)+adfac
#  endif
          ad_cff=ad_cff+(ubar(i  ,j,krhs)+                              &
#  ifdef NEARSHORE_MELLOR
     &                   ubar_stokes(i  ,j)+                            &
     &                   ubar_stokes(i+1,j)+                            &
#  endif
     &                   ubar(i+1,j,krhs))*ad_VFe(i,j)
          ad_VFe(i,j)=0.0_r8
!>        tl_UFx(i,j)=tl_cff*(vbar(i,j  ,krhs)+                         &
#  ifdef NEARSHORE_MELLOR
!>   &                        vbar_stokes(i,j  )+                       &
!>   &                        vbar_stokes(i,j+1)+                       &
#  endif
!>   &                        vbar(i,j+1,krhs))+                        &
!>   &                cff*(tl_vbar(i,j  ,krhs)+                         &
#  ifdef NEARSHORE_MELLOR
!>   &                     tl_vbar_stokes(i,j  )+                       &
!>   &                     tl_vbar_stokes(i,j+1)+                       &
#  endif
!>   &                     tl_vbar(i,j+1,krhs))
!>
          adfac=cff*ad_UFx(i,j)
          ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+adfac
          ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+adfac
#  ifdef NEARSHORE_MELLOR
          ad_vbar_stokes(i,j  )=ad_vbar_stokes(i,j  )+adfac
          ad_vbar_stokes(i,j+1)=ad_vbar_stokes(i,j+1)+adfac
#  endif
          ad_cff=ad_cff+(vbar(i,j  ,krhs)+                              &
#  ifdef NEARSHORE_MELLOR
     &                   vbar_stokes(i,j  )+                            &
     &                   vbar_stokes(i,j+1)+                            &
#  endif
     &                   vbar(i,j+1,krhs))*ad_UFx(i,j)
          ad_UFx(i,j)=0.0_r8
!>        tl_cff=0.5_r8*tl_Drhs(i,j)*fomn(i,j)
!>
          ad_Drhs(i,j)=ad_Drhs(i,j)+0.5_r8*fomn(i,j)*ad_cff
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
!!          DiaV2rhs(i,j,M2yadv)=-cff2
!!          DiaV2rhs(i,j,M2xadv)=-cff1
#  endif
!>          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
!>
            ad_fac=ad_fac-ad_rhs_vbar(i,j)
!>          tl_fac=tl_cff1+tl_cff2
!>
            ad_cff1=ad_cff1+ad_fac
            ad_cff2=ad_cff2+ad_fac
            ad_fac=0.0_r8
!>          tl_cff2=tl_VFe(i,j)-tl_VFe(i,j-1)
!>
            ad_VFe(i,j-1)=ad_VFe(i,j-1)-ad_cff2
            ad_VFe(i,j  )=ad_VFe(i,j  )+ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=tl_VFx(i+1,j)-tl_VFx(i,j)
!>
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-ad_cff1
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+ad_cff1
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  if defined DIAGNOSTICS_UV
!!          DiaU2rhs(i,j,M2xadv)=-cff1
!!          DiaU2rhs(i,j,M2yadv)=-cff2
!!          DiaU2rhs(i,j,M2hadv)=-fac
#  endif
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
!>
            ad_fac=ad_fac-ad_rhs_ubar(i,j)
!>          tl_fac=tl_cff1+tl_cff2
!>
            ad_cff1=ad_cff1+ad_fac
            ad_cff2=ad_cff2+ad_fac
            ad_fac=0.0_r8
!>          tl_cff2=tl_UFe(i,j+1)-tl_UFe(i,j)
!>
            ad_UFe(i,j  )=ad_UFe(i,j  )-ad_cff2
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=tl_UFx(i,j)-tl_UFx(i-1,j)
!>
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-ad_cff1
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+ad_cff1
            ad_cff1=0.0_r8
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
!
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
!
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
!
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
        DO j=JstrVm1,Jendp1
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
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
!>            tl_Dgrad(i,Jend+1)=tl_Dgrad(i,Jend)
!>
              ad_Dgrad(i,Jend)=ad_Dgrad(i,Jend)+ad_Dgrad(i,Jend+1)
              ad_Dgrad(i,Jend+1)=0.0_r8
!>            tl_grad (i,Jend+1)=tl_grad (i,Jend)
!>
              ad_grad (i,Jend)=ad_grad (i,Jend)+ad_grad (i,Jend+1)
              ad_grad (i,Jend+1)=0.0_r8
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
!>            tl_Dgrad(i,Jstr)=tl_Dgrad(i,Jstr+1)
!>
              ad_Dgrad(i,Jstr+1)=ad_Dgrad(i,Jstr+1)+ad_Dgrad(i,Jstr)
              ad_Dgrad(i,Jstr)=0.0_r8
!>            tl_grad (i,Jstr)=tl_grad (i,Jstr+1)
!>
              ad_grad (i,Jstr+1)=ad_grad (i,Jstr+1)+ad_grad (i,Jstr)
              ad_grad (i,Jstr)=0.0_r8
            END DO
          END IF
        END IF

        DO j=JstrVm1,Jendp1
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
        DO j=JstrV,Jend
          DO i=Istrm1,Iendp1
            grad(i,j)=vbar(i-1,j,krhs)-2.0_r8*vbar(i,j,krhs)+           &
#   ifdef NEARSHORE_MELLOR
     &                vbar_stokes(i-1,j)-2.0_r8*vbar_stokes(i,j)+       &
     &                vbar_stokes(i+1,j)+                               &
#   endif
     &                vbar(i+1,j,krhs)
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
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrV,Jend
!>            tl_grad(Iend+1,j)=tl_grad(Iend,j)
!>
               ad_grad(Iend,j)=ad_grad(Iend,j)+ad_grad(Iend+1,j)
              ad_grad(Iend+1,j)=0.0_r8
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrV,Jend
!>            tl_grad(Istr-1,j)=tl_grad(Istr,j)
!>
              ad_grad(Istr,j)=ad_grad(Istr,j)+ad_grad(Istr-1,j)
              ad_grad(Istr-1,j)=0.0_r8
            END DO
           END IF
        END IF
        DO j=JstrV,Jend
          DO i=Istrm1,Iendp1
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
        DO j=Jstrm1,Jendp1
          DO i=IstrU,Iend
            grad(i,j)=ubar(i,j-1,krhs)-2.0_r8*ubar(i,j,krhs)+           &
#   ifdef NEARSHORE_MELLOR
     &                ubar_stokes(i,j-1)-2.0_r8*ubar_stokes(i,j)+       &
     &                ubar_stokes(i,j+1)+                               &
#   endif
     &                ubar(i,j+1,krhs)
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
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrU,Iend
!>            tl_grad(i,Jend+1)=tl_grad(i,Jend)
!>
              ad_grad(i,Jend)=ad_grad(i,Jend)+ad_grad(i,Jend+1)
              ad_grad(i,Jend+1)=0.0_r8
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrU,Iend
!>            tl_grad(i,Jstr-1)=tl_grad(i,Jstr)
!>
              ad_grad(i,Jstr)=ad_grad(i,Jstr)+ad_grad(i,Jstr-1)
               ad_grad(i,Jstr-1)=0.0_r8
            END DO
          END IF
        END IF
        DO j=Jstrm1,Jendp1
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
        DO j=Jstr,Jend
          DO i=IstrUm1,Iendp1
            grad (i,j)=ubar(i-1,j,krhs)-2.0_r8*ubar(i,j,krhs)+          &
#    ifdef NEARSHORE_MELLOR
     &                 ubar_stokes(i-1,j)-2.0_r8*ubar_stokes(i,j)+      &
     &                 ubar_stokes(i+1,j)+                              &
#    endif
     &                 ubar(i+1,j,krhs)
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
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
!>            tl_Dgrad(Iend+1,j)=tl_Dgrad(Iend,j)
!>
              ad_Dgrad(Iend,j)=ad_Dgrad(Iend,j)+ad_Dgrad(Iend+1,j)
              ad_Dgrad(Iend+1,j)=0.0_r8
!>            tl_grad (Iend+1,j)=tl_grad (Iend,j)
!>
              ad_grad (Iend,j)=ad_grad (Iend,j)+ad_grad (Iend+1,j)
              ad_grad (Iend+1,j)=0.0_r8
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
!>            tl_Dgrad(Istr,j)=tl_Dgrad(Istr+1,j)
!>
              ad_Dgrad(Istr+1,j)=ad_Dgrad(Istr+1,j)+ad_Dgrad(Istr,j)
              ad_Dgrad(Istr,j)=0.0_r8
!>            tl_grad (Istr,j)=tl_grad (Istr+1,j)
!>
              ad_grad (Istr+1,j)=ad_grad (Istr+1,j)+ad_grad (Istr,j)
              ad_grad (Istr,j)=0.0_r8
            END DO
          END IF
        END IF
        DO j=Jstr,Jend
          DO i=IstrUm1,Iendp1
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
# if !defined SOLVE3D && defined ATM_PRESS
        fac3=0.5_r8*100.0_r8/rho0
# endif
        DO j=Jstr,Jend
          IF (j.ge.JstrV) THEN
            DO i=Istr,Iend
# ifdef DIAGNOSTICS_UV
!!            DiaV2rhs(i,j,M2pgrd)=rhs_vbar(i,j)
# endif
# if defined ATM_PRESS && !defined SOLVE3D
!>            tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                        &
!>   &                         fac3*om_v(i,j)*                          &
!>   &                         (tl_h(i,j-1)+tl_h(i,j)+                  &
!>   &                          tl_gzeta(i,j-1)+tl_gzeta(i,j))*         &
!>   &                         (Pair(i,j-1)-Pair(i,j))
!>
              adfac=fac3*om_v(i,j)*(Pair(i,j-1)-Pair(i,j)*              &
     &              ad_rhs_vbar(i,j)
              ad_h(i,j-1)=ad_h(i,j-1)+adfac
              ad_h(i,j  )=ad_h(i,j  )+adfac
              ad_gzeta(i,j-1)=ad_gzeta(i,j-1)+adfac
              ad_gzeta(i,j  )=ad_gzeta(i,j  )+adfac
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
# if defined ATM_PRESS && !defined SOLVE3D
!>          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                          &
!>   &                       fac3*on_u(i,j)*                            &
!>   &                       (tl_h(i-1,j)+tl_h(i,j)+                    &
!>   &                        tl_gzeta(i-1,j)+tl_gzeta(i,j))*           &
!>   &                       (Pair(i-1,j)-Pair(i,j))
!>
            adfac=fac3*on_u(i,j)*(Pair(i-1,j)-Pair(i,j))*               &
     &            ad_rhs_ubar(i,j)
            ad_h(i-1,j)=ad_h(i-1,j)+adfac
            ad_h(i  ,j)=ad_h(i  ,j)+adfac
            ad_gzeta(i-1,j)=ad_gzeta(i-1,j)+adfac
            ad_gzeta(i  ,j)=ad_gzeta(i  ,j)+adfac
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
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      tl_zeta(:,:,knew))
!>
        CALL ad_mp_exchange2d (ng, tile, iADM, 1,                       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_zeta(:,:,knew))
# endif
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>        CALL exchange_r2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_zeta(:,:,knew))
!>
          CALL ad_exchange_r2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_zeta(:,:,knew))
        END IF
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
!
!  Apply adjoint mass point sources (volume vertical influx), if any.
!
        IF (LwSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            i=SOURCES(ng)%Isrc(is)
            j=SOURCES(ng)%Jsrc(is)
            IF (((IstrR.le.i).and.(i.le.IendR)).and.                    &
     &          ((JstrR.le.j).and.(j.le.JendR))) THEN
!>            tl_zeta(i,j,knew)=tl_zeta(i,j,knew)+0.0_r8
!>
            END IF
          END DO
        END IF
!
!  If adjoint predictor step, load right-side-term into shared array.
!
        IF (PREDICTOR_2D_STEP(ng)) THEN
# ifdef DISTRIBUTE
!>        CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        NghostPoints,                             &
!>   &                        EWperiodic(ng), NSperiodic(ng),           &
!>   &                        tl_rzeta(:,:,krhs))
!>
          CALL ad_mp_exchange2d (ng, tile, iADM, 1,                     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_rzeta(:,:,krhs))
# endif
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>          CALL exchange_r2d_tile (ng, tile,                           &
!>   &                              LBi, UBi, LBj, UBj,                 &
!>   &                              tl_rzeta(:,:,krhs))
!>
            CALL ad_exchange_r2d_tile (ng, tile,                        &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 ad_rzeta(:,:,krhs))
          END IF
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

# ifdef WET_DRY_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute new wet/dry masks.
!-----------------------------------------------------------------------
!
!>    CALL wetdry_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  IminS, ImaxS, JminS, JmaxS,                     &
#  ifdef MASKING
!>   &                  pmask, rmask, umask, vmask,                     &
#  endif
!>   &                  h, zeta(:,:,kstp),                              &
#  ifdef SOLVE3D
!>   &                  DU_avg1, DV_avg1,                               &
!>   &                  rmask_wet_avg,                                  &
#  endif
!>   &                  pmask_wet, pmask_full,                          &
!>   &                  rmask_wet, rmask_full,                          &
!>   &                  umask_wet, umask_full,                          &
!>   &                  vmask_wet, vmask_full)
!>
!>  HGA: Need the ADM code here for the above NLM code.
!>
# endif
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
#  ifdef NESTING
!  In nesting applications with refinement grids, we need to exchange
!  the DU_avg2 and DV_avg2 fluxes boundary information for the case
!  that a contact point is at a tile partition. Notice that in such
!  cases, we need i+1 and j+1 values for spatial/temporal interpolation.
#  endif
!
      IF ((iif(ng).eq.(nfast(ng)+1)).and.PREDICTOR_2D_STEP(ng)) THEN

#  ifdef DISTRIBUTE
#   ifdef NESTING
!>      CALL mp_exchange2d (ng, tile, iTLM, 2,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      tl_DU_avg2, tl_DV_avg2)
!>
        CALL ad_mp_exchange2d (ng, tile, iADM, 2,                       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_DU_avg2, ad_DV_avg2)
#   endif
!>      CALL mp_exchange2d (ng, tile, iTLM, 3,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      tl_Zt_avg1, tl_DU_avg1, tl_DV_avg1)
!>
        CALL ad_mp_exchange2d (ng, tile, iADM, 3,                       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_Zt_avg1, ad_DU_avg1, ad_DV_avg1)
#  endif
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
#  ifdef NESTING
!>        CALL exchange_v2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_DV_avg2)
!>
          CALL ad_exchange_v2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_DV_avg2)
!>        CALL exchange_u2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_DU_avg2)
          CALL ad_exchange_u2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_DU_avg2)
#  endif
!>        CALL exchange_v2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_DV_avg1)
!>
          CALL ad_exchange_v2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_DV_avg1)
!>        CALL exchange_u2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_DU_avg1)
!>
          CALL ad_exchange_u2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_DU_avg1)
!>        CALL exchange_r2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_Zt_avg1)
!>
          CALL ad_exchange_r2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_Zt_avg1)
        END IF
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
!
!  Set vertically integrated mass fluxes DUon and DVom along the open
!  boundaries in such a way that the integral volume is conserved.
!
      IF (ANY(ad_VolCons(:,ng))) THEN
!>      CALL tl_set_DUV_bc_tile (ng, tile,                              &
!>   &                           LBi, UBi, LBj, UBj,                    &
!>   &                           IminS, ImaxS, JminS, JmaxS,            &
!>   &                           krhs,                                  &
# ifdef MASKING
!>   &                           umask, vmask,                          &
# endif
!>   &                           om_v, on_u, ubar, vbar,                &
!>   &                           tl_ubar, tl_vbar,                      &
!>   &                           Drhs, DUon, DVom,                      &
!>   &                           tl_Drhs, tl_DUon, tl_DVom)
!>
        CALL ad_set_DUV_bc_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           krhs,                                  &
# ifdef MASKING
     &                           umask, vmask,                          &
# endif
     &                           om_v, on_u, ubar, vbar,                &
     &                           ad_ubar, ad_vbar,                      &
     &                           Drhs, DUon, DVom,                      &
     &                           ad_Drhs, ad_DUon, ad_DVom)
      END IF

# ifdef DISTRIBUTE
!
!  In distributed-memory, the I- and J-ranges are different and a
!  special exchange is done to avoid having three ghost points for
!  high order numerical stencils. Notice that a private array is
!  passed below to the exchange routine. It also applies periodic
!  boundary conditions, if appropriate and no partitions in I- or
!  J-directions.
!
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_DUon, tl_DVom)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_DUon, ad_DVom)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          IminS, ImaxS, JminS, JmaxS,             &
!>   &                          tl_DUon)
!>
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             ad_DUon)
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          IminS, ImaxS, JminS, JmaxS,             &
!>   &                          tl_DVom)
!>
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             ad_DVom)
      END IF
# endif
# if defined DISTRIBUTE && !defined NESTING
!
!  Compute adjoint adjoint vertically integrated mass fluxes.
!
      DO j=JstrV-1,Jendp2
        DO i=IstrU-2,Iendp2
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
#  ifdef NEARSHORE_MELLOR
!>        tl_DVom(i,j)=tl_DVom(i,j)+tl_DVSom(i,j)
!>
          ad_DVSom(i,j)=ad_DVSom(i,j)+ad_DVom(i,j)
!>        tl_DVSom(i,j)=tl_vbar_stokes(i,j)*cff1+                       &
!>   &                  vbar_stokes(i,j)*tl_cff1
!>
          ad_cff1=ad_cff1+vbar_stokes(i,j)*ad_DVSom(i,j)
          ad_vbar_stokes(i,j)=ad_vbar_stokes(i,j)+cff1*ad_DVSom(i,j)
          ad_DVSom(i,j)=0.0_r8
#  endif
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
      DO j=JstrV-2,Jendp2
        DO i=IstrU-1,Iendp2
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
#  ifdef NEARSHORE_MELLOR
!>        tl_DUon(i,j)=tl_DUon(i,j)+tl_DUSon(i,j)
!>
          ad_DUSon(i,j)=ad_DUSon(i,j)+ad_DUon(i,j)
!>        tl_DUSon(i,j)=tl_ubar_stokes(i,j)*cff1+                       &
!>   &                  ubar_stokes(i,j)*tl_cff1
!>
          ad_cff1=ad_cff1+ubar_stokes(i,j)*ad_DUSon(i,j)
          ad_ubar_stokes(i,j)=ad_ubar_stokes(i,j)+cff1*ad_DUSon(i,j)
          ad_DUSon(i,j)=0.0_r8
#  endif
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
      DO j=JstrV-2,Jendp2
        DO i=IstrU-2,Iendp2
!>        tl_Drhs(i,j)=tl_zeta(i,j,krhs)+tl_h(i,j)
!>
          ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+ad_Drhs(i,j)
          ad_h(i,j)=ad_h(i,j)+ad_Drhs(i,j)
          ad_Drhs(i,j)=0.0_r8
        END DO
      END DO

# else

      DO j=JstrVm2,Jendp2
        DO i=IstrUm2-1,Iendp2
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
#  ifdef NEARSHORE_MELLOR
!>        tl_DVom(i,j)=tl_DVom(i,j)+tl_DVSom(i,j)
!>
          ad_DVSom(i,j)=ad_DVSom(i,j)+ad_DVom(i,j)
!>        tl_DVSom(i,j)=tl_vbar_stokes(i,j)*cff1+                       &
!>   &                  vbar_stokes(i,j)*tl_cff1
!>
          ad_cff1=ad_cff1+vbar_stokes(i,j)*ad_DVSom(i,j)
          ad_vbar_stokes(i,j)=ad_vbar_stokes(i,j)+cff1*ad_DVSom(i,j)
          ad_DVSom(i,j)=0.0_r8
#  endif
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
      DO j=JstrVm2-1,Jendp2
        DO i=IstrUm2,Iendp2
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
#  ifdef NEARSHORE_MELLOR
!>        tl_DUon(i,j)=tl_DUon(i,j)+tl_DUSon(i,j)
!>
          ad_DUSon(i,j)=ad_DUSon(i,j)+ad_DUon(i,j)
!>        tl_DUSon(i,j)=tl_ubar_stokes(i,j)*cff1+                       &
!>   &                  ubar_stokes(i,j)*tl_cff1
!>
          ad_cff1=ad_cff1+ubar_stokes(i,j)*ad_DUSon(i,j)
          ad_ubar_stokes(i,j)=ad_ubar_stokes(i,j)+cff1*ad_DUSon(i,j)
          ad_DUSon(i,j)=0.0_r8
#  endif
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
      DO j=JstrVm2-1,Jendp2
        DO i=IstrUm2-1,Iendp2
!>        tl_Drhs(i,j)=tl_zeta(i,j,krhs)+tl_h(i,j)
!>
          ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+ad_Drhs(i,j)
          ad_h(i,j)=ad_h(i,j)+ad_Drhs(i,j)
          ad_Drhs(i,j)=0.0_r8
        END DO
      END DO
# endif
      RETURN
      END SUBROUTINE ad_step2d_tile
#else
      SUBROUTINE ad_step2d
      END SUBROUTINE ad_step2d
#endif
