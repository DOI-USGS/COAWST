#undef DEBUG
#ifdef TANGENT
      SUBROUTINE tl_step2d (ng, tile)
!
!svn $Id: tl_step2d_LF_AM3.h 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!                                                                      !
!  Tangent linear model shallow-water primitive equations predictor    !
!  (Leap-frog) and corrector (Adams-Moulton) time-stepping engine.     !
!                                                                      !
!=======================================================================
!
      USE mod_param
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
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 9)
# endif
      CALL tl_step2d_tile (ng, tile,                                    &
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
     &                     GRID(ng) % h,           GRID(ng) % tl_h,     &
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
     &                     MIXING(ng) % tl_rustr2d,                     &
     &                     MIXING(ng) % tl_rvstr2d,                     &
     &                     OCEAN(ng) % tl_rulag2d,                      &
     &                     OCEAN(ng) % tl_rvlag2d,                      &
     &                     OCEAN(ng) % ubar_stokes,                     &
     &                     OCEAN(ng) % tl_ubar_stokes,                  &
     &                     OCEAN(ng) % vbar_stokes,                     &
     &                     OCEAN(ng) % tl_vbar_stokes,                  &
# endif
# ifndef SOLVE3D
     &                     FORCES(ng) % tl_bustr,                       &
     &                     FORCES(ng) % tl_bvstr,                       &
#  ifdef ATM_PRESS
     &                     FORCES(ng) % Pair,                           &
#  endif
# else
#  ifdef VAR_RHO_2D
     &                     COUPLING(ng) % rhoA,                         &
     &                     COUPLING(ng) % tl_rhoA,                      &
     &                     COUPLING(ng) % rhoS,                         &
     &                     COUPLING(ng) % tl_rhoS,                      &
#  endif
     &                     COUPLING(ng) % tl_DU_avg1,                   &
     &                     COUPLING(ng) % tl_DU_avg2,                   &
     &                     COUPLING(ng) % tl_DV_avg1,                   &
     &                     COUPLING(ng) % tl_DV_avg2,                   &
     &                     COUPLING(ng) % Zt_avg1,                      &
     &                     COUPLING(ng) % tl_Zt_avg1,                   &
     &                     COUPLING(ng) % tl_rufrc,                     &
     &                     COUPLING(ng) % tl_rvfrc,                     &
     &                     OCEAN(ng) % tl_ru,                           &
     &                     OCEAN(ng) % tl_rv,                           &
# endif
# ifdef DIAGNOSTICS_UV
!!   &                     DIAGS(ng) % DiaU2wrk,   DIAGS(ng) % DiaV2wrk,&
!!   &                     DIAGS(ng) % DiaRUbar,   DIAGS(ng) % DiaRVbar,&
#  ifdef SOLVE3D
!!   &                     DIAGS(ng) % DiaU2int,   DIAGS(ng) % DiaV2int,&
!!   &                     DIAGS(ng) % DiaRUfrc,   DIAGS(ng) % DiaRVfrc,&
#  endif
# endif
     &                     OCEAN(ng) % rubar,      OCEAN(ng) % tl_rubar,&
     &                     OCEAN(ng) % rvbar,      OCEAN(ng) % tl_rvbar,&
     &                     OCEAN(ng) % rzeta,      OCEAN(ng) % tl_rzeta,&
     &                     OCEAN(ng) % ubar,       OCEAN(ng) % tl_ubar, &
     &                     OCEAN(ng) % vbar,       OCEAN(ng) % tl_vbar, &
     &                     OCEAN(ng) % zeta,       OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 9)
# endif

      RETURN
      END SUBROUTINE tl_step2d
!
!***********************************************************************
      SUBROUTINE tl_step2d_tile (ng, tile,                              &
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
     &                           h, tl_h,                               &
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
     &                           tl_rustr2d, tl_rvstr2d,                &
     &                           tl_rulag2d, tl_rvlag2d,                &
     &                           ubar_stokes, tl_ubar_stokes,           &
     &                           vbar_stokes, tl_vbar_stokes,           &
# endif
# ifndef SOLVE3D
     &                           tl_bustr, tl_bvstr,                    &
#  ifdef ATM_PRESS
     &                           Pair,                                  &
#  endif
# else
#  ifdef VAR_RHO_2D
     &                           rhoA, tl_rhoA, rhoS, tl_rhoS,          &
#  endif
     &                           tl_DU_avg1, tl_DU_avg2,                &
     &                           tl_DV_avg1, tl_DV_avg2,                &
     &                           Zt_avg1, tl_Zt_avg1,                   &
     &                           tl_rufrc, tl_rvfrc,                    &
     &                           tl_ru, tl_rv,                          &
# endif
# ifdef DIAGNOSTICS_UV
!!   &                           DiaU2wrk, DiaV2wrk,                    &
!!   &                           DiaRUbar, DiaRVbar,                    &
#  ifdef SOLVE3D
!!   &                           DiaU2int, DiaV2int,                    &
!!   &                           DiaRUfrc, DiaRVfrc,                    &
#  endif
# endif
     &                           rubar, tl_rubar,                       &
     &                           rvbar, tl_rvbar,                       &
     &                           rzeta, tl_rzeta,                       &
     &                           ubar, tl_ubar,                         &
     &                           vbar, tl_vbar,                         &
     &                           zeta, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
      USE mod_sources
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE obc_volcons_mod
      USE tl_obc_volcons_mod
      USE tl_u2dbc_mod, ONLY : tl_u2dbc_tile
      USE tl_v2dbc_mod, ONLY : tl_v2dbc_tile
      USE tl_zetabc_mod, ONLY : tl_zetabc_tile
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

      real(r8), intent(in) :: tl_h(LBi:,LBj:)
#  ifndef SOLVE3D
      real(r8), intent(in) :: tl_bustr(LBi:,LBj:)
      real(r8), intent(in) :: tl_bvstr(LBi:,LBj:)
#   ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:,LBj:)
#   endif
#  else
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: rhoA(LBi:,LBj:)
      real(r8), intent(in) :: rhoS(LBi:,LBj:)
      real(r8), intent(in) :: tl_rhoA(LBi:,LBj:)
      real(r8), intent(in) :: tl_rhoS(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)

      real(r8), intent(inout) :: tl_DU_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: tl_DU_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: tl_DV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: tl_DV_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: tl_rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: tl_rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: tl_ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tl_rv(LBi:,LBj:,0:,:)
#  endif
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: tl_rustr2d(LBi:,LBj:)
      real(r8), intent(inout) :: tl_rvstr2d(LBi:,LBj:)
      real(r8), intent(inout) :: tl_rulag2d(LBi:,LBj:)
      real(r8), intent(inout) :: tl_rvlag2d(LBi:,LBj:)
      real(r8), intent(inout) :: tl_ubar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: tl_vbar_stokes(LBi:,LBj:)
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
      real(r8), intent(inout) :: tl_rubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_rvbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_rzeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

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

      real(r8), intent(in) :: tl_h(LBi:UBi,LBj:UBj)
#  ifndef SOLVE3D
      real(r8), intent(in) :: tl_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_bvstr(LBi:UBi,LBj:UBj)
#   ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
#   endif
#  else
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rhoS(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_rhoS(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: tl_DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_ru(LBi:UBi,LBj:UBj,0:UBk,2)
      real(r8), intent(inout) :: tl_rv(LBi:UBi,LBj:UBj,0:UBk,2)
#  endif
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: tl_rustr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_rvstr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_rulag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_rvlag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_vbar_stokes(LBi:UBi,LBj:UBj)
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
      real(r8), intent(inout) :: tl_rubar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: tl_rvbar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: tl_rzeta(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
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
      real(r8) :: tl_cff, tl_cff1, tl_cff2, tl_cff3, tl_cff4
      real(r8) :: tl_fac, tl_fac1

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
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,NDM2d-1) :: DiaU2rhs
!!    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,NDM2d-1) :: DiaV2rhs
# endif

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Dgrad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Dnew
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Drhs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Drhs_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_DUon
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_DVom
# ifdef NEARSHORE_MELLOR
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_DUSon
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_DVSom
# endif
# ifdef UV_VIS4
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_LapU
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_LapV
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_grad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_gzeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_gzeta2
# if defined VAR_RHO_2D && defined SOLVE3D
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_gzetaSA
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_rhs_ubar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_rhs_vbar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_rhs_zeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_zeta_new
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_zwrk
# ifdef WET_DRY_NOT_YET
!>    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_wetdry
# endif

# include "set_bounds.h"
!
      ptsk=3-kstp
      CORRECTOR_2D_STEP=.not.PREDICTOR_2D_STEP(ng)
# ifdef DEBUG
      IF (Master) THEN
        WRITE (20,20) iic(ng), PREDICTOR_2D_STEP(ng),                   &
     &                kstp, krhs, knew, ptsk
 20     FORMAT (' iic = ',i5.5,' predictor = ',l1,' kstp = ',i1,        &
     &          ' krhs = ',i1,' knew = ',i1,' ptsk = ',i1)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute total depth (m) and vertically integrated mass fluxes.
!-----------------------------------------------------------------------
!
# if defined DISTRIBUTE && !defined NESTING

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
          tl_Drhs(i,j)=tl_zeta(i,j,krhs)+tl_h(i,j)
        END DO
      END DO
      DO j=JstrV-2,Jendp2
        DO i=IstrU-1,Iendp2
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
          tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i-1,j))
          DUon(i,j)=ubar(i,j,krhs)*cff1
          tl_DUon(i,j)=tl_ubar(i,j,krhs)*cff1+                          &
     &                 ubar(i,j,krhs)*tl_cff1
#  ifdef NEARSHORE_MELLOR
          DUSon(i,j)=ubar_stokes(i,j)*cff1
          tl_DUSon(i,j)=tl_ubar_stokes(i,j)*cff1+                       &
     &                  ubar_stokes(i,j)*tl_cff1
          DUon(i,j)=DUon(i,j)+DUSon(i,j)
          tl_DUon(i,j)=tl_DUon(i,j)+tl_DUSon(i,j)
#  endif
        END DO
      END DO
      DO j=JstrV-1,Jendp2
        DO i=IstrU-2,Iendp2
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
          tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i,j-1))
          DVom(i,j)=vbar(i,j,krhs)*cff1
          tl_DVom(i,j)=tl_vbar(i,j,krhs)*cff1+                          &
     &                 vbar(i,j,krhs)*tl_cff1
#  ifdef NEARSHORE_MELLOR
          DVSom(i,j)=vbar_stokes(i,j)*cff1
          tl_DVSom(i,j)=tl_vbar_stokes(i,j)*cff1+                       &
     &                  vbar_stokes(i,j)*tl_cff1
          DVom(i,j)=DVom(i,j)+DVSom(i,j)
          tl_DVom(i,j)=tl_DVom(i,j)+tl_DVSom(i,j)
#  endif
        END DO
      END DO

# else

      DO j=JstrVm2-1,Jendp2
        DO i=IstrUm2-1,Iendp2
          Dnew(i,j)=zeta(i,j,knew)+h(i,j)
          Drhs(i,j)=zeta(i,j,krhs)+h(i,j)
          tl_Drhs(i,j)=tl_zeta(i,j,krhs)+tl_h(i,j)
        END DO
      END DO
      DO j=JstrVm2-1,Jendp2
        DO i=IstrUm2,Iendp2
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
          tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i-1,j))
          DUon(i,j)=ubar(i,j,krhs)*cff1
          tl_DUon(i,j)=tl_ubar(i,j,krhs)*cff1+                          &
     &                 ubar(i,j,krhs)*tl_cff1
#  ifdef NEARSHORE_MELLOR
          DUSon(i,j)=ubar_stokes(i,j)*cff1
          tl_DUSon(i,j)=tl_ubar_stokes(i,j)*cff1+                       &
     &                  ubar_stokes(i,j)*tl_cff1
          DUon(i,j)=DUon(i,j)+DUSon(i,j)
          tl_DUon(i,j)=tl_DUon(i,j)+tl_DUSon(i,j)
#  endif
        END DO
      END DO
      DO j=JstrVm2,Jendp2
        DO i=IstrUm2-1,Iendp2
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
          tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i,j-1))
          DVom(i,j)=vbar(i,j,krhs)*cff1
          tl_DVom(i,j)=tl_vbar(i,j,krhs)*cff1+                          &
     &                 vbar(i,j,krhs)*tl_cff1
#  ifdef NEARSHORE_MELLOR
          DVSom(i,j)=vbar_stokes(i,j)*cff1
          tl_DVSom(i,j)=tl_vbar_stokes(i,j)*cff1+                       &
     &                  vbar_stokes(i,j)*tl_cff1
          DVom(i,j)=DVom(i,j)+DVSom(i,j)
          tl_DVom(i,j)=tl_DVom(i,j)+tl_DVSom(i,j)
#  endif
        END DO
      END DO
# endif
# ifdef DISTRIBUTE
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          DUon)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          tl_DUon)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          DVom)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          tl_DVom)
      END IF

      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    DUon, DVom)
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_DUon, tl_DVom)
# endif
# if !defined FORWARD_RHS
!
!  Initialize right-hand-side terms.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          rhs_ubar(i,j)=0.0_r8
          rhs_vbar(i,j)=0.0_r8
        END DO
      END DO
# endif
!
!  Compute integral mass flux across open boundaries and adjust
!  for volume conservation. Compute BASIC STATE value.
!  This needs to be computed here instead of below.
!
      IF (ANY(tl_VolCons(:,ng))) THEN
        CALL obc_flux_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      knew,                                       &
# ifdef MASKING
     &                      umask, vmask,                               &
# endif
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
# ifdef MASKING
     &                        umask, vmask,                             &
# endif
     &                        om_v, on_u,                               &
     &                        ubar, vbar,                               &
     &                        Drhs, DUon, DVom)
        CALL tl_set_DUV_bc_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           krhs,                                  &
# ifdef MASKING
     &                           umask, vmask,                          &
# endif
     &                           om_v, on_u, ubar, vbar,                &
     &                           tl_ubar, tl_vbar,                      &
     &                           Drhs, DUon, DVom,                      &
     &                           tl_Drhs, tl_DUon, tl_DVom)
      END IF
# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time averaged fields over all short time-steps.
!-----------------------------------------------------------------------
!
      IF (PREDICTOR_2D_STEP(ng)) THEN
        IF (FIRST_2D_STEP) THEN
!
!  Reset arrays for 2D fields averaged within the short time-steps.
!
          cff2=(-1.0_r8/12.0_r8)*weight(2,iif(ng)+1,ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
!>            Zt_avg1(i,j)=0.0_r8
!>
              tl_Zt_avg1(i,j)=0.0_r8
            END DO
            DO i=Istr,IendR
!>            DU_avg1(i,j)=0.0_r8
!>
              tl_DU_avg1(i,j)=0.0_r8
!>            DU_avg2(i,j)=cff2*DUon(i,j)
!>
              tl_DU_avg2(i,j)=cff2*tl_DUon(i,j)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
!>            DV_avg1(i,j)=0.0_r8
!>
              tl_DV_avg1(i,j)=0.0_r8
!>            DV_avg2(i,j)=cff2*DVom(i,j)
!>
              tl_DV_avg2(i,j)=cff2*tl_DVom(i,j)
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
!>            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,krhs)
!>
              tl_Zt_avg1(i,j)=tl_Zt_avg1(i,j)+cff1*tl_zeta(i,j,krhs)
            END DO
            DO i=Istr,IendR
!>            DU_avg1(i,j)=DU_avg1(i,j)+cff1*DUon(i,j)
!>
              tl_DU_avg1(i,j)=tl_DU_avg1(i,j)+cff1*tl_DUon(i,j)
#  ifdef NEARSHORE_MELLOR
!>            DU_avg1(i,j)=DU_avg1(i,j)-cff1*DUSon(i,j)
!>
              tl_DU_avg1(i,j)=tl_DU_avg1(i,j)-cff1*tl_DUSon(i,j)
#  endif
!>            DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
!>
              tl_DU_avg2(i,j)=tl_DU_avg2(i,j)+cff2*tl_DUon(i,j)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
!>            DV_avg1(i,j)=DV_avg1(i,j)+cff1*DVom(i,j)
!>
              tl_DV_avg1(i,j)=tl_DV_avg1(i,j)+cff1*tl_DVom(i,j)
#  ifdef NEARSHORE_MELLOR
!>            DV_avg1(i,j)=DV_avg1(i,j)-cff1*DVSom(i,j)
!>
              tl_DV_avg1(i,j)=tl_DV_avg1(i,j)-cff1*tl_DVSom(i,j)
#  endif
!>            DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
!>
              tl_DV_avg2(i,j)=tl_DV_avg2(i,j)+cff2*tl_DVom(i,j)
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
!>          DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
!>
            tl_DU_avg2(i,j)=tl_DU_avg2(i,j)+cff2*tl_DUon(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
!>          DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
!>
            tl_DV_avg2(i,j)=tl_DV_avg2(i,j)+cff2*tl_DVom(i,j)
          END DO
        END DO
      END IF
!
!  After all fast time steps are completed, apply boundary conditions
!  to time averaged fields.
#  ifdef NESTING
!  In nesting applications with refinement grids, we need to exchange
!  the DU_avg2 and DV_avg2 fluxes boundary information for the case
!  that a contact point is at a tile partition. Notice that in such
!  cases, we need i+1 and j+1 values for spatial/temporal interpolation.
#  endif
!
      IF ((iif(ng).eq.(nfast(ng)+1)).and.PREDICTOR_2D_STEP(ng)) THEN
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>        CALL exchange_r2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            Zt_avg1)
!>
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_Zt_avg1)
!>        CALL exchange_u2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            DU_avg1)
!>
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_DU_avg1)
!>        CALL exchange_v2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            DV_avg1)
!>
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_DV_avg1)
#  ifdef NESTING
!>        CALL exchange_u2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            DU_avg2)
!>
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_DU_avg2)
!>        CALL exchange_v2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            DV_avg2)
!>
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_DV_avg2)
#  endif
        END IF

#  ifdef DISTRIBUTE
!>      CALL mp_exchange2d (ng, tile, iNLM, 3,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      Zt_avg1, DU_avg1, DV_avg1)
!>
        CALL mp_exchange2d (ng, tile, iTLM, 3,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_Zt_avg1, tl_DU_avg1, tl_DV_avg1)
#   ifdef NESTING
!>      CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      DU_avg2, DV_avg2)
!>
        CALL mp_exchange2d (ng, tile, iTLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_DU_avg2, tl_DV_avg2)
#   endif
#  endif
      END IF
# endif
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
!>  HGA: Need the TLM code here.
!>
# endif
!
!  Do not perform the actual time stepping during the auxiliary
!  (nfast(ng)+1) time step.
!
      IF (iif(ng).gt.nfast(ng)) RETURN
!
!=======================================================================
!  Time step free-surface equation.
!=======================================================================
!
!  During the first time-step, the predictor step is Forward-Euler
!  and the corrector step is Backward-Euler. Otherwise, the predictor
!  step is Leap-frog and the corrector step is Adams-Moulton.
!
# if defined VAR_RHO_2D && defined SOLVE3D
      fac=1000.0_r8/rho0
# endif
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE && \
    !defined SOLVE3D
      IF (FIRST_2D_STEP.and.SOinitial(ng)) THEN
# else
      IF (FIRST_2D_STEP) THEN
# endif
        cff1=dtfast(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            rhs_zeta(i,j)=(DUon(i,j)-DUon(i+1,j))+                      &
     &                    (DVom(i,j)-DVom(i,j+1))
            tl_rhs_zeta(i,j)=(tl_DUon(i,j)-tl_DUon(i+1,j))+             &
     &                       (tl_DVom(i,j)-tl_DVom(i,j+1))
            zeta_new(i,j)=zeta(i,j,kstp)+                               &
     &                    pm(i,j)*pn(i,j)*cff1*rhs_zeta(i,j)
            tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+                         &
     &                       pm(i,j)*pn(i,j)*cff1*tl_rhs_zeta(i,j)
# ifdef MASKING
            zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
            tl_zeta_new(i,j)=tl_zeta_new(i,j)*rmask(i,j)
# endif
!>          Dnew(i,j)=zeta_new(i,j)+h(i,j)
!>
            tl_Dnew(i,j)=tl_zeta_new(i,j)+tl_h(i,j)
!
            zwrk(i,j)=0.5_r8*(zeta(i,j,kstp)+zeta_new(i,j))
            tl_zwrk(i,j)=0.5_r8*(tl_zeta(i,j,kstp)+tl_zeta_new(i,j))
# if defined VAR_RHO_2D && defined SOLVE3D
            gzeta(i,j)=(fac+rhoS(i,j))*zwrk(i,j)
            tl_gzeta(i,j)=(fac+rhoS(i,j))*tl_zwrk(i,j)+                 &
     &                    tl_rhoS(i,j)*zwrk(i,j)
            gzeta2(i,j)=gzeta(i,j)*zwrk(i,j)
            tl_gzeta2(i,j)=tl_gzeta(i,j)*zwrk(i,j)+                     &
     &                     gzeta(i,j)*tl_zwrk(i,j)
            gzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
            tl_gzetaSA(i,j)=tl_zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))+         &
     &                      zwrk(i,j)*(tl_rhoS(i,j)-tl_rhoA(i,j))
# else
            gzeta(i,j)=zwrk(i,j)
            tl_gzeta(i,j)=tl_zwrk(i,j)
            gzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
            tl_gzeta2(i,j)=2.0_r8*tl_zwrk(i,j)*zwrk(i,j)
# endif
          END DO
        END DO
      ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
        cff1=2.0_r8*dtfast(ng)
        cff4=4.0_r8/25.0_r8
        cff5=1.0_r8-2.0_r8*cff4
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            rhs_zeta(i,j)=(DUon(i,j)-DUon(i+1,j))+                      &
     &                    (DVom(i,j)-DVom(i,j+1))
            tl_rhs_zeta(i,j)=(tl_DUon(i,j)-tl_DUon(i+1,j))+             &
     &                       (tl_DVom(i,j)-tl_DVom(i,j+1))
            zeta_new(i,j)=zeta(i,j,kstp)+                               &
     &                    pm(i,j)*pn(i,j)*cff1*rhs_zeta(i,j)
            tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+                         &
     &                       pm(i,j)*pn(i,j)*cff1*tl_rhs_zeta(i,j)
# ifdef MASKING
            zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
            tl_zeta_new(i,j)=tl_zeta_new(i,j)*rmask(i,j)
# endif
!>          Dnew(i,j)=zeta_new(i,j)+h(i,j)
!>
            tl_Dnew(i,j)=tl_zeta_new(i,j)+tl_h(i,j)
!
            zwrk(i,j)=cff5*zeta(i,j,krhs)+                              &
     &                cff4*(zeta(i,j,kstp)+zeta_new(i,j))
            tl_zwrk(i,j)=cff5*tl_zeta(i,j,krhs)+                        &
     &                   cff4*(tl_zeta(i,j,kstp)+tl_zeta_new(i,j))
# if defined VAR_RHO_2D && defined SOLVE3D
            gzeta(i,j)=(fac+rhoS(i,j))*zwrk(i,j)
            tl_gzeta(i,j)=(fac+rhoS(i,j))*tl_zwrk(i,j)+                 &
     &                    tl_rhoS(i,j)*zwrk(i,j)
            gzeta2(i,j)=gzeta(i,j)*zwrk(i,j)
            tl_gzeta2(i,j)=tl_gzeta(i,j)*zwrk(i,j)+                     &
     &                     gzeta(i,j)*tl_zwrk(i,j)
            gzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
            tl_gzetaSA(i,j)=tl_zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))+         &
     &                      zwrk(i,j)*(tl_rhoS(i,j)-tl_rhoA(i,j))
# else
            gzeta(i,j)=zwrk(i,j)
            tl_gzeta(i,j)=tl_zwrk(i,j)
            gzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
            tl_gzeta2(i,j)=2.0_r8*tl_zwrk(i,j)*zwrk(i,j)
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
            cff=cff1*((DUon(i,j)-DUon(i+1,j))+                          &
     &                (DVom(i,j)-DVom(i,j+1)))
            tl_cff=cff1*((tl_DUon(i,j)-tl_DUon(i+1,j))+                 &
     &                   (tl_DVom(i,j)-tl_DVom(i,j+1)))
            zeta_new(i,j)=zeta(i,j,kstp)+                               &
     &                    pm(i,j)*pn(i,j)*(cff+                         &
     &                                     cff2*rzeta(i,j,kstp)-        &
     &                                     cff3*rzeta(i,j,ptsk))
            tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+                         &
     &                       pm(i,j)*pn(i,j)*(tl_cff+                   &
     &                                        cff2*tl_rzeta(i,j,kstp)-  &
     &                                        cff3*tl_rzeta(i,j,ptsk))
# ifdef MASKING
            zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
            tl_zeta_new(i,j)=tl_zeta_new(i,j)*rmask(i,j)
# endif
!>          Dnew(i,j)=zeta_new(i,j)+h(i,j)
!>
            tl_Dnew(i,j)=tl_zeta_new(i,j)+tl_h(i,j)
!
            zwrk(i,j)=cff5*zeta_new(i,j)+cff4*zeta(i,j,krhs)
            tl_zwrk(i,j)=cff5*tl_zeta_new(i,j)+cff4*tl_zeta(i,j,krhs)
# if defined VAR_RHO_2D && defined SOLVE3D
            gzeta(i,j)=(fac+rhoS(i,j))*zwrk(i,j)
            tl_gzeta(i,j)=(fac+rhoS(i,j))*tl_zwrk(i,j)+                 &
     &                    tl_rhoS(i,j)*zwrk(i,j)
            gzeta2(i,j)=gzeta(i,j)*zwrk(i,j)
            tl_gzeta2(i,j)=tl_gzeta(i,j)*zwrk(i,j)+                     &
     &                     gzeta(i,j)*tl_zwrk(i,j)
            gzetaSA(i,j)=zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))
            tl_gzetaSA(i,j)=tl_zwrk(i,j)*(rhoS(i,j)-rhoA(i,j))+         &
     &                      zwrk(i,j)*(tl_rhoS(i,j)-tl_rhoA(i,j))
# else
            gzeta(i,j)=zwrk(i,j)
            tl_gzeta(i,j)=tl_zwrk(i,j)
            gzeta2(i,j)=zwrk(i,j)*zwrk(i,j)
            tl_gzeta2(i,j)=2.0_r8*tl_zwrk(i,j)*zwrk(i,j)
# endif
          END DO
        END DO
      END IF
!
!  Load new free-surface values into shared array at both predictor
!  and corrector steps.
# ifdef WET_DRY_NOT_YET
!  Modify new free-surface to Ensure that depth is > Dcrit for masked
!  cells.
# endif
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!>        zeta(i,j,knew)=zeta_new(i,j)
!>
          tl_zeta(i,j,knew)=tl_zeta_new(i,j)
# if defined WET_DRY_NOT_YET && defined MASKING
!>        zeta(i,j,knew)=zeta(i,j,knew)+                                &
!>   &                   (Dcrit(ng)-h(i,j))*(1.0_r8-rmask(i,j))
!>
          tl_zeta(i,j,knew)=tl_zeta(i,j,knew)-                          &
     &                      tl_h(i,j)*(1.0_r8-rmask(i,j))
# endif
        END DO
      END DO
!
!  If predictor step, load right-side-term into shared array.
!
      IF (PREDICTOR_2D_STEP(ng)) THEN
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          rzeta(i,j,krhs)=rhs_zeta(i,j)
!>
            tl_rzeta(i,j,krhs)=tl_rhs_zeta(i,j)
          END DO
        END DO

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>        CALL exchange_r2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            rzeta(:,:,krhs))
!>
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_rzeta(:,:,krhs))
        END IF

# ifdef DISTRIBUTE
!>      CALL mp_exchange2d (ng, tile, iNLM, 1,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      rzeta(:,:,krhs))
!>
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_rzeta(:,:,krhs))
# endif
      END IF
!
!  Apply mass point sources (volume vertical influx), if any.
!
      IF (LwSrc(ng)) THEN
        DO is=1,Nsrc(ng)
          i=SOURCES(ng)%Isrc(is)
          j=SOURCES(ng)%Jsrc(is)
          IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
!>          zeta(i,j,knew)=zeta(i,j,knew)+                              &
!>   &                     SOURCES(ng)%Qbar(is)*                        &
!>   &                     pm(i,j)*pn(i,j)*dtfast(ng)
!>
!!          tl_zeta(i,j,knew)=tl_zeta(i,j,knew)+0.0_r8
          END IF
        END DO
      END IF
!
!  Set free-surface lateral boundary conditions.
!
!>    CALL zetabc_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  IminS, ImaxS, JminS, JmaxS,                     &
!>   &                  krhs, kstp, knew,                               &
!>   &                  zeta)
!>
      CALL tl_zetabc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, knew,                            &
     &                     zeta, tl_zeta)

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          zeta(:,:,knew))
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_zeta(:,:,knew))
      END IF

# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    zeta(:,:,knew))
!>
      CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_zeta(:,:,knew))
# endif
!
!=======================================================================
!  Compute right-hand-side for the 2D momentum equations.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Compute pressure gradient terms.
!-----------------------------------------------------------------------
!
      cff1=0.5_r8*g
      cff2=1.0_r8/3.0_r8
# if !defined SOLVE3D && defined ATM_PRESS
      fac3=0.5_r8*100.0_r8/rho0
# endif
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        rhs_ubar(i,j)=cff1*on_u(i,j)*                                 &
!>   &                  ((h(i-1,j)+                                     &
!>   &                    h(i  ,j))*                                    &
!>   &                   (gzeta(i-1,j)-                                 &
!>   &                    gzeta(i  ,j))+                                &
# if defined VAR_RHO_2D && defined SOLVE3D
!>   &                   (h(i-1,j)-                                     &
!>   &                    h(i  ,j))*                                    &
!>   &                   (gzetaSA(i-1,j)+                               &
!>   &                    gzetaSA(i  ,j)+                               &
!>   &                    cff2*(rhoA(i-1,j)-                            &
!>   &                          rhoA(i  ,j))*                           &
!>   &                         (zwrk(i-1,j)-                            &
!>   &                          zwrk(i  ,j)))+                          &
# endif
!>   &                   (gzeta2(i-1,j)-                                &
!>   &                    gzeta2(i  ,j)))
!>
          tl_rhs_ubar(i,j)=cff1*on_u(i,j)*                              &
     &                     ((tl_h(i-1,j)+                               &
     &                       tl_h(i ,j))*                               &
     &                      (gzeta(i-1,j)-                              &
     &                       gzeta(i  ,j))+                             &
     &                      (h(i-1,j)+                                  &
     &                       h(i ,j))*                                  &
     &                      (tl_gzeta(i-1,j)-                           &
     &                       tl_gzeta(i  ,j))+                          &
# if defined VAR_RHO_2D && defined SOLVE3D
     &                      (tl_h(i-1,j)-                               &
     &                       tl_h(i  ,j))*                              &
     &                      (gzetaSA(i-1,j)+                            &
     &                       gzetaSA(i  ,j)+                            &
     &                       cff2*(rhoA(i-1,j)-                         &
     &                             rhoA(i  ,j))*                        &
     &                            (zwrk(i-1,j)-                         &
     &                             zwrk(i  ,j)))+                       &
     &                      (h(i-1,j)-                                  &
     &                       h(i  ,j))*                                 &
     &                      (tl_gzetaSA(i-1,j)+                         &
     &                       tl_gzetaSA(i  ,j)+                         &
     &                       cff2*((tl_rhoA(i-1,j)-                     &
     &                              tl_rhoA(i  ,j))*                    &
     &                             (zwrk(i-1,j)-                        &
     &                              zwrk(i  ,j))+                       &
     &                             (rhoA(i-1,j)-                        &
     &                              rhoA(i  ,j))*                       &
     &                             (tl_zwrk(i-1,j)-                     &
     &                              tl_zwrk(i  ,j))))+                  &
# endif
     &                      (tl_gzeta2(i-1,j)-                          &
     &                       tl_gzeta2(i  ,j)))
# if defined ATM_PRESS && !defined SOLVE3D
!>        rhs_ubar(i,j)=rhs_ubar(i,j)+                                  &
!>   &                  fac3*on_u(i,j)*                                 &
!>   &                  (h(i-1,j)+h(i,j)+                               &
!>   &                   gzeta(i-1,j)+gzeta(i,j))*                      &
!>   &                  (Pair(i-1,j)-Pair(i,j))
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                            &
     &                     fac3*on_u(i,j)*                              &
     &                     (tl_h(i-1,j)+tl_h(i,j)+                      &
     &                      tl_gzeta(i-1,j)+tl_gzeta(i,j))*             &
     &                     (Pair(i-1,j)-Pair(i,j))
# endif
# ifdef DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2pgrd)=rhs_ubar(i,j)
# endif
        END DO
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
!>          rhs_vbar(i,j)=cff1*om_v(i,j)*                               &
!>   &                    ((h(i,j-1)+                                   &
!>   &                      h(i,j  ))*                                  &
!>   &                     (gzeta(i,j-1)-                               &
!>   &                      gzeta(i,j  ))+                              &
# if defined VAR_RHO_2D && defined SOLVE3D
!>   &                     (h(i,j-1)-                                   &
!>   &                      h(i,j  ))*                                  &
!>   &                     (gzetaSA(i,j-1)+                             &
!>   &                      gzetaSA(i,j  )+                             &
!>   &                      cff2*(rhoA(i,j-1)-                          &
!>   &                            rhoA(i,j  ))*                         &
!>   &                           (zwrk(i,j-1)-                          &
!>   &                            zwrk(i,j  )))+                        &
# endif
!>   &                     (gzeta2(i,j-1)-                              &
!>   &                      gzeta2(i,j  )))
!>
            tl_rhs_vbar(i,j)=cff1*om_v(i,j)*                            &
     &                       ((tl_h(i,j-1)+                             &
     &                         tl_h(i,j  ))*                            &
     &                        (gzeta(i,j-1)-                            &
     &                         gzeta(i,j  ))+                           &
     &                        (h(i,j-1)+                                &
     &                         h(i,j  ))*                               &
     &                        (tl_gzeta(i,j-1)-                         &
     &                         tl_gzeta(i,j  ))+                        &
# if defined VAR_RHO_2D && defined SOLVE3D
     &                        (tl_h(i,j-1)-                             &
     &                         tl_h(i,j  ))*                            &
     &                        (gzetaSA(i,j-1)+                          &
     &                         gzetaSA(i,j  )+                          &
     &                         cff2*(rhoA(i,j-1)-                       &
     &                               rhoA(i,j  ))*                      &
     &                              (zwrk(i,j-1)-                       &
     &                               zwrk(i,j  )))+                     &
     &                        (h(i,j-1)-                                &
     &                         h(i,j  ))*                               &
     &                        (tl_gzetaSA(i,j-1)+                       &
     &                         tl_gzetaSA(i,j  )+                       &
     &                         cff2*((tl_rhoA(i,j-1)-                   &
     &                                tl_rhoA(i,j  ))*                  &
     &                               (zwrk(i,j-1)-                      &
     &                                zwrk(i,j  ))+                     &
     &                               (rhoA(i,j-1)-                      &
     &                                rhoA(i,j  ))*                     &
     &                               (tl_zwrk(i,j-1)-                   &
     &                                tl_zwrk(i,j  ))))+                &
# endif
     &                        (tl_gzeta2(i,j-1)-                        &
     &                         tl_gzeta2(i,j  )))
# if defined ATM_PRESS && !defined SOLVE3D
!>          rhs_vbar(i,j)=rhs_vbar(i,j)+                                &
!>   &                    fac3*om_v(i,j)*                               &
!>   &                    (h(i,j-1)+h(i,j)+                             &
!>   &                     gzeta(i,j-1)+gzeta(i,j))*                    &
!>   &                    (Pair(i,j-1)-Pair(i,j))
!>
            tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                          &
     &                       fac3*om_v(i,j)*                            &
     &                       (tl_h(i,j-1)+tl_h(i,j)+                    &
     &                        tl_gzeta(i,j-1)+tl_gzeta(i,j))*           &
     &                       (Pair(i,j-1)-Pair(i,j))
# endif
# ifdef DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2pgrd)=rhs_vbar(i,j)
# endif
          END DO
        END IF
      END DO
# ifdef UV_ADV
!
!-----------------------------------------------------------------------
!  Add in horizontal advection of momentum.
!-----------------------------------------------------------------------

#  ifdef UV_C2ADVECTION
!
!  Second-order, centered differences advection.
!
      DO j=Jstr,Jend
        DO i=IstrU-1,Iend
!>        UFx(i,j)=0.25_r8*(DUon(i,j)+DUon(i+1,j))*                     &
!>   &                     (ubar(i  ,j,krhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      ubar_stokes(i  ,j)+                         &
!>   &                      ubar_stokes(i+1,j)+                         &
#   endif
!>   &                      ubar(i+1,j,krhs))
!>
          tl_UFx(i,j)=0.25_r8*                                          &
     &                ((tl_DUon(i,j)+tl_DUon(i+1,j))*                   &
     &                 (ubar(i  ,j,krhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                  ubar_stokes(i  ,j)+                             &
     &                  ubar_stokes(i+1,j)+                             &
#   endif
     &                  ubar(i+1,j,krhs))+                              &
     &                 (DUon(i,j)+DUon(i+1,j))*                         &
     &                 (tl_ubar(i  ,j,krhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  tl_ubar_stokes(i  ,j)+                          &
     &                  tl_ubar_stokes(i+1,j)+                          &
#   endif
     &                  tl_ubar(i+1,j,krhs)))
        END DO
      END DO
!
      DO j=Jstr,Jend+1
        DO i=IstrU,Iend
!>        UFe(i,j)=0.25_r8*(DVom(i,j)+DVom(i-1,j))*                     &
!>   &                     (ubar(i,j  ,krhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      ubar_stokes(i,j  )+                         &
!>   &                      ubar_stokes(i,j-1)+                         &
#   endif
!>   &                      ubar(i,j-1,krhs))
!>
          tl_UFe(i,j)=0.25_r8*                                          &
     &                ((tl_DVom(i,j)+tl_DVom(i-1,j))*                   &
     &                 (ubar(i,j  ,krhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                  ubar_stokes(i,j  )+                             &
     &                  ubar_stokes(i,j-1)+                             &
#   endif
     &                  ubar(i,j-1,krhs))+                              &
     &                 (DVom(i,j)+DVom(i-1,j))*                         &
     &                 (tl_ubar(i,j  ,krhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  tl_ubar_stokes(i,j  )+                          &
     &                  tl_ubar_stokes(i,j-1)+                          &
#   endif
     &                  tl_ubar(i,j-1,krhs)))
        END DO
      END DO
!
      DO j=JstrV,Jend
        DO i=Istr,Iend+1
!>        VFx(i,j)=0.25_r8*(DUon(i,j)+DUon(i,j-1))*                     &
!>   &                     (vbar(i  ,j,krhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      vbar_stokes(i  ,j)+                         &
!>   &                      vbar_stokes(i-1,j)+                         &
#   endif
!>   &                      vbar(i-1,j,krhs))
!>
          tl_VFx(i,j)=0.25_r8*                                          &
     &                ((tl_DUon(i,j)+tl_DUon(i,j-1))*                   &
     &                 (vbar(i  ,j,krhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                  vbar_stokes(i  ,j)+                             &
     &                  vbar_stokes(i-1,j)+                             &
#   endif
     &                  vbar(i-1,j,krhs))+                              &
     &                 (DUon(i,j)+DUon(i,j-1))*                         &
     &                 (tl_vbar(i  ,j,krhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  tl_vbar_stokes(i  ,j)+                          &
     &                  tl_vbar_stokes(i-1,j)+                          &
#   endif
     &                  tl_vbar(i-1,j,krhs)))
        END DO
      END DO
!
      DO j=JstrV-1,Jend
        DO i=Istr,Iend
!>        VFe(i,j)=0.25_r8*(DVom(i,j)+DVom(i,j+1))*                     &
!>   &                     (vbar(i,j  ,krhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      vbar_stokes(i,j  )+                         &
!>   &                      vbar_stokes(i,j+1)+                         &
#   endif
!>   &                      vbar(i,j+1,krhs))
!>
          tl_VFe(i,j)=0.25_r8*                                          &
     &                ((tl_DVom(i,j)+tl_DVom(i,j+1))*                   &
     &                 (vbar(i,j  ,krhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                  vbar_stokes(i,j  )+                             &
     &                  vbar_stokes(i,j+1)+                             &
#   endif
     &                  vbar(i,j+1,krhs))+                              &
     &                 (DVom(i,j)+DVom(i,j+1))*                         &
     &                 (tl_vbar(i,j  ,krhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  tl_vbar_stokes(i,j  )+                          &
     &                  tl_vbar_stokes(i,j+1)+                          &
#   endif
     &                  tl_vbar(i,j+1,krhs)))
        END DO
      END DO
#  else
!
!  Fourth-order, centered differences advection.
!
      DO j=Jstr,Jend
        DO i=IstrUm1,Iendp1
          grad (i,j)=ubar(i-1,j,krhs)-2.0_r8*ubar(i,j,krhs)+            &
#   ifdef NEARSHORE_MELLOR
     &               ubar_stokes(i-1,j)-2.0_r8*ubar_stokes(i,j)+        &
     &               ubar_stokes(i+1,j)+                                &
#   endif
     &               ubar(i+1,j,krhs)
          tl_grad(i,j)=tl_ubar(i-1,j,krhs)-2.0_r8*tl_ubar(i,j,krhs)+    &
#   ifdef NEARSHORE_MELLOR
     &                 tl_ubar_stokes(i-1,j)-2.0_r8*tl_ubar_stokes(i,j)+&
     &                 tl_ubar_stokes(i+1,j)+                           &
#   endif
     &                 tl_ubar(i+1,j,krhs)
          Dgrad(i,j)=DUon(i-1,j)-2.0_r8*DUon(i,j)+DUon(i+1,j)
          tl_Dgrad(i,j)=tl_DUon(i-1,j)-2.0_r8*tl_DUon(i,j)+             &
     &                  tl_DUon(i+1,j)
        END DO
      END DO
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            grad (Istr,j)=grad (Istr+1,j)
            tl_grad (Istr,j)=tl_grad (Istr+1,j)
            Dgrad(Istr,j)=Dgrad(Istr+1,j)
            tl_Dgrad(Istr,j)=tl_Dgrad(Istr+1,j)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            grad (Iend+1,j)=grad (Iend,j)
            tl_grad (Iend+1,j)=tl_grad (Iend,j)
            Dgrad(Iend+1,j)=Dgrad(Iend,j)
            tl_Dgrad(Iend+1,j)=tl_Dgrad(Iend,j)
          END DO
        END IF
      END IF

      cff=1.0_r8/6.0_r8
      DO j=Jstr,Jend
        DO i=IstrU-1,Iend
!>        UFx(i,j)=0.25_r8*(ubar(i  ,j,krhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      ubar_stokes(i  ,j)+                         &
!>   &                      ubar_stokes(i+1,j)+                         &
#   endif
!>   &                      ubar(i+1,j,krhs)-                           &
!>   &                      cff*(grad (i,j)+grad (i+1,j)))*             &
!>   &                     (DUon(i,j)+DUon(i+1,j)-                      &
!>   &                      cff*(Dgrad(i,j)+Dgrad(i+1,j)))
!>
          tl_UFx(i,j)=0.25_r8*                                          &
     &                ((ubar(i  ,j,krhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                  ubar_stokes(i  ,j)+                             &
     &                  ubar_stokes(i+1,j)+                             &
#   endif
     &                  ubar(i+1,j,krhs)-                               &
     &                  cff*(grad (i,j)+grad (i+1,j)))*                 &
     &                 (tl_DUon(i,j)+tl_DUon(i+1,j)-                    &
     &                  cff*(tl_Dgrad(i,j)+tl_Dgrad(i+1,j)))+           &
     &                 (tl_ubar(i  ,j,krhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  tl_ubar_stokes(i  ,j)+                          &
     &                  tl_ubar_stokes(i+1,j)+                          &
#   endif
     &                  tl_ubar(i+1,j,krhs)-                            &
     &                  cff*(tl_grad (i,j)+tl_grad (i+1,j)))*           &
     &                 (DUon(i,j)+DUon(i+1,j)-                          &
     &                  cff*(Dgrad(i,j)+Dgrad(i+1,j))))
        END DO
      END DO
!
      DO j=Jstrm1,Jendp1
        DO i=IstrU,Iend
          grad(i,j)=ubar(i,j-1,krhs)-2.0_r8*ubar(i,j,krhs)+             &
#   ifdef NEARSHORE_MELLOR
     &              ubar_stokes(i,j-1)-2.0_r8*ubar_stokes(i,j)+         &
     &              ubar_stokes(i,j+1)+                                 &
#   endif
     &              ubar(i,j+1,krhs)
          tl_grad(i,j)=tl_ubar(i,j-1,krhs)-2.0_r8*tl_ubar(i,j,krhs)+    &
#   ifdef NEARSHORE_MELLOR
     &                 tl_ubar_stokes(i,j-1)-2.0_r8*tl_ubar_stokes(i,j)+&
     &                 tl_ubar_stokes(i,j+1)+                           &
#   endif
     &                 tl_ubar(i,j+1,krhs)
        END DO
      END DO
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            grad(i,Jstr-1)=grad(i,Jstr)
            tl_grad(i,Jstr-1)=tl_grad(i,Jstr)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
            grad(i,Jend+1)=grad(i,Jend)
            tl_grad(i,Jend+1)=tl_grad(i,Jend)
          END DO
        END IF
      END IF
      DO j=Jstr,Jend+1
        DO i=IstrU-1,Iend
          Dgrad(i,j)=DVom(i-1,j)-2.0_r8*DVom(i,j)+DVom(i+1,j)
          tl_Dgrad(i,j)=tl_DVom(i-1,j)-2.0_r8*tl_DVom(i,j)+             &
     &                  tl_DVom(i+1,j)
        END DO
      END DO

      cff=1.0_r8/6.0_r8
      DO j=Jstr,Jend+1
        DO i=IstrU,Iend
!>        UFe(i,j)=0.25_r8*(ubar(i,j  ,krhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      ubar_stokes(i,j  )+                         &
!>   &                      ubar_stokes(i,j-1)+                         &
#   endif
!>   &                      ubar(i,j-1,krhs)-                           &
!>   &                      cff*(grad (i,j)+grad (i,j-1)))*             &
!>   &                     (DVom(i,j)+DVom(i-1,j)-                      &
!>   &                      cff*(Dgrad(i,j)+Dgrad(i-1,j)))
!>
          tl_UFe(i,j)=0.25_r8*                                          &
     &                ((tl_ubar(i,j  ,krhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  tl_ubar_stokes(i,j  )+                          &
     &                  tl_ubar_stokes(i,j-1)+                          &
#   endif
     &                  tl_ubar(i,j-1,krhs)-                            &
     &                  cff*(tl_grad (i,j)+tl_grad (i,j-1)))*           &
     &                 (DVom(i,j)+DVom(i-1,j)-                          &
     &                  cff*(Dgrad(i,j)+Dgrad(i-1,j)))+                 &
     &                 (ubar(i,j  ,krhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                  ubar_stokes(i,j  )+                             &
     &                  ubar_stokes(i,j-1)+                             &
#   endif
     &                  ubar(i,j-1,krhs)-                               &
     &                  cff*(grad (i,j)+grad (i,j-1)))*                 &
     &                 (tl_DVom(i,j)+tl_DVom(i-1,j)-                    &
     &                  cff*(tl_Dgrad(i,j)+tl_Dgrad(i-1,j))))
        END DO
      END DO
!
      DO j=JstrV,Jend
        DO i=Istrm1,Iendp1
          grad(i,j)=vbar(i-1,j,krhs)-2.0_r8*vbar(i,j,krhs)+             &
#   ifdef NEARSHORE_MELLOR
     &              vbar_stokes(i-1,j)-2.0_r8*vbar_stokes(i,j)+         &
     &              vbar_stokes(i+1,j)+                                 &
#   endif
     &              vbar(i+1,j,krhs)
          tl_grad(i,j)=tl_vbar(i-1,j,krhs)-2.0_r8*tl_vbar(i,j,krhs)+    &
#   ifdef NEARSHORE_MELLOR
     &                 tl_vbar_stokes(i-1,j)-2.0_r8*tl_vbar_stokes(i,j)+&
     &                 tl_vbar_stokes(i+1,j)+                           &
#   endif
     &                 tl_vbar(i+1,j,krhs)
        END DO
      END DO
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
            grad(Istr-1,j)=grad(Istr,j)
            tl_grad(Istr-1,j)=tl_grad(Istr,j)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            grad(Iend+1,j)=grad(Iend,j)
            tl_grad(Iend+1,j)=tl_grad(Iend,j)
          END DO
        END IF
      END IF
      DO j=JstrV-1,Jend
        DO i=Istr,Iend+1
          Dgrad(i,j)=DUon(i,j-1)-2.0_r8*DUon(i,j)+DUon(i,j+1)
          tl_Dgrad(i,j)=tl_DUon(i,j-1)-2.0_r8*tl_DUon(i,j)+             &
     &                  tl_DUon(i,j+1)
        END DO
      END DO

      cff=1.0_r8/6.0_r8
      DO j=JstrV,Jend
        DO i=Istr,Iend+1
!>        VFx(i,j)=0.25_r8*(vbar(i  ,j,krhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      vbar_stokes(i  ,j)+                         &
!>   &                      vbar_stokes(i-1,j)+                         &
#   endif
!>   &                      vbar(i-1,j,krhs)-                           &
!>   &                      cff*(grad (i,j)+grad (i-1,j)))*             &
!>   &                     (DUon(i,j)+DUon(i,j-1)-                      &
!>   &                      cff*(Dgrad(i,j)+Dgrad(i,j-1)))
!>
          tl_VFx(i,j)=0.25_r8*                                          &
     &                ((tl_vbar(i  ,j,krhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  tl_vbar_stokes(i  ,j)+                          &
     &                  tl_vbar_stokes(i-1,j)+                          &
#   endif
     &                  tl_vbar(i-1,j,krhs)-                            &
     &                  cff*(tl_grad (i,j)+tl_grad (i-1,j)))*           &
     &                 (DUon(i,j)+DUon(i,j-1)-                          &
     &                  cff*(Dgrad(i,j)+Dgrad(i,j-1)))+                 &
     &                 (vbar(i  ,j,krhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                  vbar_stokes(i  ,j)+                             &
     &                  vbar_stokes(i-1,j)+                             &
#   endif
     &                  vbar(i-1,j,krhs)-                               &
     &                  cff*(grad (i,j)+grad (i-1,j)))*                 &
     &                 (tl_DUon(i,j)+tl_DUon(i,j-1)-                    &
     &                  cff*(tl_Dgrad(i,j)+tl_Dgrad(i,j-1))))
        END DO
      END DO
!
      DO j=JstrVm1,Jendp1
        DO i=Istr,Iend
          grad(i,j)=vbar(i,j-1,krhs)-2.0_r8*vbar(i,j,krhs)+             &
#   ifdef NEARSHORE_MELLOR
     &              vbar_stokes(i,j-1)-2.0_r8*vbar_stokes(i,j)+         &
     &              vbar_stokes(i,j+1)+                                 &
#   endif
     &              vbar(i,j+1,krhs)
          tl_grad(i,j)=tl_vbar(i,j-1,krhs)-2.0_r8*tl_vbar(i,j,krhs)+    &
#   ifdef NEARSHORE_MELLOR
     &                 tl_vbar_stokes(i,j-1)-2.0_r8*tl_vbar_stokes(i,j)+&
     &                 tl_vbar_stokes(i,j+1)+                           &
#   endif
     &                 tl_vbar(i,j+1,krhs)
          Dgrad(i,j)=DVom(i,j-1)-2.0_r8*DVom(i,j)+DVom(i,j+1)
          tl_Dgrad(i,j)=tl_DVom(i,j-1)-2.0_r8*tl_DVom(i,j)+             &
     &                  tl_DVom(i,j+1)
        END DO
      END DO
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            grad (i,Jstr)=grad (i,Jstr+1)
            tl_grad (i,Jstr)=tl_grad (i,Jstr+1)
            Dgrad(i,Jstr)=Dgrad(i,Jstr+1)
            tl_Dgrad(i,Jstr)=tl_Dgrad(i,Jstr+1)
          END DO
        END IF
      END IF
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            grad (i,Jend+1)=grad (i,Jend)
            tl_grad (i,Jend+1)=tl_grad (i,Jend)
            Dgrad(i,Jend+1)=Dgrad(i,Jend)
            tl_Dgrad(i,Jend+1)=tl_Dgrad(i,Jend)
          END DO
        END IF
      END IF

      cff=1.0_r8/6.0_r8
      DO j=JstrV-1,Jend
        DO i=Istr,Iend
!>        VFe(i,j)=0.25_r8*(vbar(i,j  ,krhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      vbar_stokes(i,j  )+                         &
!>   &                      vbar_stokes(i,j+1)+                         &
#   endif
!>   &                      vbar(i,j+1,krhs)-                           &
!>   &                      cff*(grad (i,j)+grad (i,j+1)))*             &
!>   &                     (DVom(i,j)+DVom(i,j+1)-                      &
!>   &                      cff*(Dgrad(i,j)+Dgrad(i,j+1)))
!>
          tl_VFe(i,j)=0.25_r8*                                          &
     &                ((tl_vbar(i,j  ,krhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  tl_vbar_stokes(i,j  )+                          &
     &                  tl_vbar_stokes(i,j+1)+                          &
#   endif
     &                  tl_vbar(i,j+1,krhs)-                            &
     &                  cff*(tl_grad (i,j)+tl_grad (i,j+1)))*           &
     &                 (DVom(i,j)+DVom(i,j+1)-                          &
     &                  cff*(Dgrad(i,j)+Dgrad(i,j+1)))+                 &
     &                 (vbar(i,j  ,krhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                  vbar_stokes(i,j  )+                             &
     &                  vbar_stokes(i,j+1)+                             &
#   endif
     &                  vbar(i,j+1,krhs)-                               &
     &                  cff*(grad (i,j)+grad (i,j+1)))*                 &
     &                 (tl_DVom(i,j)+tl_DVom(i,j+1)-                    &
     &                  cff*(tl_Dgrad(i,j)+tl_Dgrad(i,j+1))))
        END DO
      END DO
#  endif
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        cff1=UFx(i,j)-UFx(i-1,j)
!>
          tl_cff1=tl_UFx(i,j)-tl_UFx(i-1,j)
!>        cff2=UFe(i,j+1)-UFe(i,j)
!>
          tl_cff2=tl_UFe(i,j+1)-tl_UFe(i,j)
!>        fac=cff1+cff2
!>
          tl_fac=tl_cff1+tl_cff2
!>        rhs_ubar(i,j)=rhs_ubar(i,j)-fac
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2xadv)=-cff1
!!        DiaU2rhs(i,j,M2yadv)=-cff2
!!        DiaU2rhs(i,j,M2hadv)=-fac
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        cff1=VFx(i+1,j)-VFx(i,j)
!>
          tl_cff1=tl_VFx(i+1,j)-tl_VFx(i,j)
!>        cff2=VFe(i,j)-VFe(i,j-1)
!>
          tl_cff2=tl_VFe(i,j)-tl_VFe(i,j-1)
!>        fac=cff1+cff2
!>
          tl_fac=tl_cff1+tl_cff2
!>        rhs_vbar(i,j)=rhs_vbar(i,j)-fac
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2xadv)=-cff1
!!        DiaV2rhs(i,j,M2yadv)=-cff2
!!        DiaV2rhs(i,j,M2hadv)=-fac
#  endif
        END DO
      END DO
# endif
# ifdef UV_COR
!
!-----------------------------------------------------------------------
!  Add in Coriolis term.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff=0.5_r8*Drhs(i,j)*fomn(i,j)
          tl_cff=0.5_r8*tl_Drhs(i,j)*fomn(i,j)
!>        UFx(i,j)=cff*(vbar(i,j  ,krhs)+                               &
#  ifdef NEARSHORE_MELLOR
!>   &                  vbar_stokes(i,j  )+                             &
!>   &                  vbar_stokes(i,j+1)+                             &
#  endif
!>   &                  vbar(i,j+1,krhs))
!>
          tl_UFx(i,j)=tl_cff*(vbar(i,j  ,krhs)+                         &
#  ifdef NEARSHORE_MELLOR
     &                        vbar_stokes(i,j  )+                       &
     &                        vbar_stokes(i,j+1)+                       &
#  endif
     &                        vbar(i,j+1,krhs))+                        &
     &                cff*(tl_vbar(i,j  ,krhs)+                         &
#  ifdef NEARSHORE_MELLOR
     &                     tl_vbar_stokes(i,j  )+                       &
     &                     tl_vbar_stokes(i,j+1)+                       &
#  endif
     &                     tl_vbar(i,j+1,krhs))
!>        VFe(i,j)=cff*(ubar(i  ,j,krhs)+                               &
#  ifdef NEARSHORE_MELLOR
!>   &                  ubar_stokes(i  ,j)+                             &
!>   &                  ubar_stokes(i+1,j)+                             &
#  endif
!>   &                  ubar(i+1,j,krhs))
!>
          tl_VFe(i,j)=tl_cff*(ubar(i  ,j,krhs)+                         &
#  ifdef NEARSHORE_MELLOR
     &                        ubar_stokes(i  ,j)+                       &
     &                        ubar_stokes(i+1,j)+                       &
#  endif
     &                        ubar(i+1,j,krhs))+                        &
     &                cff*(tl_ubar(i  ,j,krhs)+                         &
#  ifdef NEARSHORE_MELLOR
     &                     tl_ubar_stokes(i  ,j)+                       &
     &                     tl_ubar_stokes(i+1,j)+                       &
#  endif
     &                     tl_ubar(i+1,j,krhs))
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        fac1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
!>
          tl_fac1=0.5_r8*(tl_UFx(i,j)+tl_UFx(i-1,j))
!>        rhs_ubar(i,j)=rhs_ubar(i,j)+fac1
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac1
#  if defined DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2fcor)=fac1
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        fac1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
!>
          tl_fac1=0.5_r8*(tl_VFe(i,j)+tl_VFe(i,j-1))
!>        rhs_vbar(i,j)=rhs_vbar(i,j)-fac1
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac1
#  if defined DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2fcor)=-fac1
#  endif
        END DO
      END DO
# endif
# if defined CURVGRID && defined UV_ADV
!
!-----------------------------------------------------------------------
!  Add in curvilinear transformation terms.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=0.5_r8*(vbar(i,j  ,krhs)+                                &
#  ifdef NEARSHORE_MELLOR
     &                 vbar_stokes(i,j  )+                              &
     &                 vbar_stokes(i,j+1)+                              &
#  endif
     &                 vbar(i,j+1,krhs))
          tl_cff1=0.5_r8*(tl_vbar(i,j  ,krhs)+                          &
#  ifdef NEARSHORE_MELLOR
     &                    tl_vbar_stokes(i,j  )+                        &
     &                    tl_vbar_stokes(i,j+1)+                        &
#  endif
     &                    tl_vbar(i,j+1,krhs))
          cff2=0.5_r8*(ubar(i  ,j,krhs)+                                &
#  ifdef NEARSHORE_MELLOR
     &                 ubar_stokes(i  ,j)+                              &
     &                 ubar_stokes(i+1,j)+                              &
#  endif
     &                 ubar(i+1,j,krhs))
          tl_cff2=0.5_r8*(tl_ubar(i  ,j,krhs)+                          &
#  ifdef NEARSHORE_MELLOR
     &                    tl_ubar_stokes(i  ,j)+                        &
     &                    tl_ubar_stokes(i+1,j)+                        &
#  endif
     &                    tl_ubar(i+1,j,krhs))
          cff3=cff1*dndx(i,j)
          tl_cff3=tl_cff1*dndx(i,j)
          cff4=cff2*dmde(i,j)
          tl_cff4=tl_cff2*dmde(i,j)
          cff=Drhs(i,j)*(cff3-cff4)
          tl_cff=tl_Drhs(i,j)*(cff3-cff4)+                              &
     &           Drhs(i,j)*(tl_cff3-tl_cff4)
!>        UFx(i,j)=cff*cff1
!>
          tl_UFx(i,j)=tl_cff*cff1+cff*tl_cff1
!>        VFe(i,j)=cff*cff2
!>
          tl_VFe(i,j)=tl_cff*cff2+cff*tl_cff2
#  if defined DIAGNOSTICS_UV
!!        cff=Drhs(i,j)*cff4
!!        Uwrk(i,j)=-cff*cff1                  ! ubar equation, ETA-term
!!        Vwrk(i,j)=-cff*cff2                  ! vbar equation, ETA-term
#  endif
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        fac1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
!>
          tl_fac1=0.5_r8*(tl_UFx(i,j)+tl_UFx(i-1,j))
!>        rhs_ubar(i,j)=rhs_ubar(i,j)+fac1
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac1
#  if defined DIAGNOSTICS_UV
!!        fac2=0.5_r8*(Uwrk(i,j)+Uwrk(i-1,j))
!!        DiaU2rhs(i,j,M2xadv)=DiaU2rhs(i,j,M2xadv)+fac1-fac2
!!        DiaU2rhs(i,j,M2yadv)=DiaU2rhs(i,j,M2yadv)+fac2
!!        DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+fac1
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        fac1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
!>
          tl_fac1=0.5_r8*(tl_VFe(i,j)+tl_VFe(i,j-1))
!>        rhs_vbar(i,j)=rhs_vbar(i,j)-fac1
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac1
#  if defined DIAGNOSTICS_UV
!!        fac2=0.5_r8*(Vwrk(i,j)+Vwrk(i,j-1))
!!        DiaV2rhs(i,j,M2xadv)=DiaV2rhs(i,j,M2xadv)-fac1+fac2
!!        DiaV2rhs(i,j,M2yadv)=DiaV2rhs(i,j,M2yadv)-fac2
!!        DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)-fac1
#  endif
        END DO
      END DO
# endif
# if defined UV_VIS2 || defined UV_VIS4
!
!-----------------------------------------------------------------------
!  If horizontal mixing, compute total depth at PSI-points.
!-----------------------------------------------------------------------
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
          tl_Drhs_p(i,j)=0.25_r8*(tl_Drhs(i,j  )+tl_Drhs(i-1,j  )+      &
     &                            tl_Drhs(i,j-1)+tl_Drhs(i-1,j-1))
        END DO
      END DO
# endif
# ifdef UV_VIS2
!
!-----------------------------------------------------------------------
!  Add in horizontal harmonic viscosity.
!-----------------------------------------------------------------------
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!>        cff=visc2_r(i,j)*Drhs(i,j)*0.5_r8*                            &
!>   &        (pmon_r(i,j)*                                             &
!>   &         ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,krhs)-                 &
!>   &          (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,krhs))-                &
!>   &         pnom_r(i,j)*                                             &
!>   &         ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,krhs)-                 &
!>   &          (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,krhs)))
!>
          tl_cff=visc2_r(i,j)*0.5_r8*                                   &
     &           (tl_Drhs(i,j)*                                         &
     &            (pmon_r(i,j)*                                         &
     &             ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,krhs)-             &
     &              (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,krhs))-            &
     &             pnom_r(i,j)*                                         &
     &             ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,krhs)-             &
     &              (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,krhs)))+           &
     &            Drhs(i,j)*                                            &
     &            (pmon_r(i,j)*                                         &
     &             ((pn(i  ,j)+pn(i+1,j))*tl_ubar(i+1,j,krhs)-          &
     &              (pn(i-1,j)+pn(i  ,j))*tl_ubar(i  ,j,krhs))-         &
     &             pnom_r(i,j)*                                         &
     &             ((pm(i,j  )+pm(i,j+1))*tl_vbar(i,j+1,krhs)-          &
     &              (pm(i,j-1)+pm(i,j  ))*tl_vbar(i,j  ,krhs))))
!>        UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
!>
          tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
!>        VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
!>
          tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
        END DO
      END DO
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
!>        cff=visc2_p(i,j)*Drhs_p(i,j)*0.5_r8*                          &
!>   &        (pmon_p(i,j)*                                             &
!>   &         ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,krhs)-               &
!>   &          (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,krhs))+              &
!>   &         pnom_p(i,j)*                                             &
!>   &         ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,krhs)-               &
!>   &          (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,krhs)))
!>
          tl_cff=visc2_p(i,j)*0.5_r8*                                   &
     &           (tl_Drhs_p(i,j)*                                       &
     &            (pmon_p(i,j)*                                         &
     &             ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,krhs)-           &
     &              (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,krhs))+          &
     &             pnom_p(i,j)*                                         &
     &             ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,krhs)-           &
     &              (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,krhs)))+         &
     &            Drhs_p(i,j)*                                          &
     &            (pmon_p(i,j)*                                         &
     &             ((pn(i  ,j-1)+pn(i  ,j))*tl_vbar(i  ,j,krhs)-        &
     &              (pn(i-1,j-1)+pn(i-1,j))*tl_vbar(i-1,j,krhs))+       &
     &             pnom_p(i,j)*                                         &
     &             ((pm(i-1,j  )+pm(i,j  ))*tl_ubar(i,j  ,krhs)-        &
     &              (pm(i-1,j-1)+pm(i,j-1))*tl_ubar(i,j-1,krhs))))
#  ifdef MASKING
!>        cff=cff*pmask(i,j)
!>
          tl_cff=tl_cff*pmask(i,j)
#  endif
!>        UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
!>
          tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
!>        VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
!>
          tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
        END DO
      END DO
!
!  Add in harmonic viscosity.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*(UFx(i,j  )-UFx(i-1,j))
!>
          tl_cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*                           &
     &            (tl_UFx(i,j  )-tl_UFx(i-1,j))
!>        cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i  ,j))
!>
          tl_cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*                           &
     &            (tl_UFe(i,j+1)-tl_UFe(i  ,j))
!>        fac=cff1+cff2
!>
          tl_fac=tl_cff1+tl_cff2
!>        rhs_ubar(i,j)=rhs_ubar(i,j)+fac
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2hvis)=fac
!!        DiaU2rhs(i,j,M2xvis)=cff1
!!        DiaU2rhs(i,j,M2yvis)=cff2
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j  ))
!>
          tl_cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*                           &
     &            (tl_VFx(i+1,j)-tl_VFx(i,j  ))
!>        cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*(VFe(i  ,j)-VFe(i,j-1))
!>
          tl_cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*                           &
     &            (tl_VFe(i  ,j)-tl_VFe(i,j-1))
!>        fac=cff1-cff2
!>
          tl_fac=tl_cff1-tl_cff2
!>        rhs_vbar(i,j)=rhs_vbar(i,j)+fac
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2hvis)=fac
!!        DiaV2rhs(i,j,M2xvis)= cff1
!!        DiaV2rhs(i,j,M2yvis)=-cff2
#  endif
        END DO
      END DO
# endif
# ifdef UV_VIS4
!
!-----------------------------------------------------------------------
!  Add in horizontal biharmonic viscosity. The biharmonic operator
!  is computed by applying the harmonic operator twice.
!-----------------------------------------------------------------------
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m4 s^-3/2) in XI- and ETA-directions.  It is assumed here
!  that "visc4_r" and "visc4_p" are the squared root of the biharmonic
!  viscosity coefficient.  For momentum balance purposes, the total
!  thickness "D" appears only when computing the second harmonic
!  operator.
!
      DO j=JstrVm2,Jendp1
        DO i=IstrUm2,Iendp1
          cff=visc4_r(i,j)*0.5_r8*                                      &
     &        (pmon_r(i,j)*                                             &
     &         ((pn(i  ,j)+pn(i+1,j))*ubar(i+1,j,krhs)-                 &
     &          (pn(i-1,j)+pn(i  ,j))*ubar(i  ,j,krhs))-                &
     &         pnom_r(i,j)*                                             &
     &         ((pm(i,j  )+pm(i,j+1))*vbar(i,j+1,krhs)-                 &
     &          (pm(i,j-1)+pm(i,j  ))*vbar(i,j  ,krhs)))
          tl_cff=visc4_r(i,j)*0.5_r8*                                   &
     &           (pmon_r(i,j)*                                          &
     &            ((pn(i  ,j)+pn(i+1,j))*tl_ubar(i+1,j,krhs)-           &
     &             (pn(i-1,j)+pn(i  ,j))*tl_ubar(i  ,j,krhs))-          &
     &            pnom_r(i,j)*                                          &
     &            ((pm(i,j  )+pm(i,j+1))*tl_vbar(i,j+1,krhs)-           &
     &             (pm(i,j-1)+pm(i,j  ))*tl_vbar(i,j  ,krhs)))
          UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
          tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
          VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
          tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
        END DO
      END DO
      DO j=Jstrm1,Jendp2
        DO i=Istrm1,Iendp2
          cff=visc4_p(i,j)*0.5_r8*                                      &
     &        (pmon_p(i,j)*                                             &
     &         ((pn(i  ,j-1)+pn(i  ,j))*vbar(i  ,j,krhs)-               &
     &          (pn(i-1,j-1)+pn(i-1,j))*vbar(i-1,j,krhs))+              &
     &         pnom_p(i,j)*                                             &
     &         ((pm(i-1,j  )+pm(i,j  ))*ubar(i,j  ,krhs)-               &
     &          (pm(i-1,j-1)+pm(i,j-1))*ubar(i,j-1,krhs)))
          tl_cff=visc4_p(i,j)*0.5_r8*                                   &
     &           (pmon_p(i,j)*                                          &
     &            ((pn(i  ,j-1)+pn(i  ,j))*tl_vbar(i  ,j,krhs)-         &
     &             (pn(i-1,j-1)+pn(i-1,j))*tl_vbar(i-1,j,krhs))+        &
     &            pnom_p(i,j)*                                          &
     &            ((pm(i-1,j  )+pm(i,j  ))*tl_ubar(i,j  ,krhs)-         &
     &             (pm(i-1,j-1)+pm(i,j-1))*tl_ubar(i,j-1,krhs)))
#  ifdef MASKING
          cff=cff*pmask(i,j)
          tl_cff=tl_cff*pmask(i,j)
#  endif
          UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
          tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
          VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
          tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
        END DO
      END DO
!
!  Compute first harmonic operator (m s^-3/2).
!
      DO j=Jstrm1,Jendp1
        DO i=IstrUm1,Iendp1
          LapU(i,j)=0.125_r8*                                           &
     &              (pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*            &
     &              ((pn(i-1,j)+pn(i,j))*                               &
     &               (UFx(i,j  )-UFx(i-1,j))+                           &
     &               (pm(i-1,j)+pm(i,j))*                               &
     &               (UFe(i,j+1)-UFe(i  ,j)))
          tl_LapU(i,j)=0.125_r8*                                        &
     &                 (pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*         &
     &                 ((pn(i-1,j)+pn(i,j))*                            &
     &                  (tl_UFx(i,j  )-tl_UFx(i-1,j))+                  &
     &                  (pm(i-1,j)+pm(i,j))*                            &
     &                  (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
        END DO
      END DO
      DO j=JstrVm1,Jendp1
        DO i=Istrm1,Iendp1
          LapV(i,j)=0.125_r8*                                           &
     &              (pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))*            &
     &              ((pn(i,j-1)+pn(i,j))*                               &
     &               (VFx(i+1,j)-VFx(i,j  ))-                           &
     &               (pm(i,j-1)+pm(i,j))*                               &
     &               (VFe(i  ,j)-VFe(i,j-1)))
          tl_LapV(i,j)=0.125_r8*                                        &
     &                 (pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))*         &
     &                 ((pn(i,j-1)+pn(i,j))*                            &
     &                  (tl_VFx(i+1,j)-tl_VFx(i,j  ))-                  &
     &                  (pm(i,j-1)+pm(i,j))*                            &
     &                  (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
        END DO
      END DO
!
!  Apply boundary conditions (other than periodic) to the first
!  harmonic operator. These are gradient or closed (free slip or
!  no slip) boundary conditions.
!
      IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (tl_LBC(iwest,isUbar,ng)%closed) THEN
            DO j=Jstrm1,Jendp1
              LapU(IstrU-1,j)=0.0_r8
              tl_LapU(IstrU-1,j)=0.0_r8
            END DO
          ELSE
            DO j=Jstrm1,Jendp1
              LapU(IstrU-1,j)=LapU(IstrU,j)
              tl_LapU(IstrU-1,j)=tl_LapU(IstrU,j)
            END DO
          END IF
          IF (tl_LBC(iwest,isVbar,ng)%closed) THEN
            DO j=JstrVm1,Jendp1
              LapV(Istr-1,j)=gamma2(ng)*LapV(Istr,j)
              tl_LapV(Istr-1,j)=gamma2(ng)*tl_LapV(Istr,j)
            END DO
          ELSE
            DO j=JstrVm1,Jendp1
              LapV(Istr-1,j)=0.0_r8
              tl_LapV(Istr-1,j)=0.0_r8
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (tl_LBC(ieast,isUbar,ng)%closed) THEN
            DO j=Jstrm1,Jendp1
              LapU(Iend+1,j)=0.0_r8
              tl_LapU(Iend+1,j)=0.0_r8
            END DO
          ELSE
            DO j=Jstrm1,Jendp1
              LapU(Iend+1,j)=LapU(Iend,j)
              tl_LapU(Iend+1,j)=tl_LapU(Iend,j)
            END DO
          END IF
          IF (tl_LBC(ieast,isVbar,ng)%closed) THEN
            DO j=JstrVm1,Jendp1
              LapV(Iend+1,j)=gamma2(ng)*LapV(Iend,j)
              tl_LapV(Iend+1,j)=gamma2(ng)*tl_LapV(Iend,j)
            END DO
          ELSE
            DO j=JstrVm1,Jendp1
              LapV(Iend+1,j)=0.0_r8
              tl_LapV(Iend+1,j)=0.0_r8
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (tl_LBC(isouth,isUbar,ng)%closed) THEN
            DO i=IstrUm1,Iendp1
              LapU(i,Jstr-1)=gamma2(ng)*LapU(i,Jstr)
              tl_LapU(i,Jstr-1)=gamma2(ng)*tl_LapU(i,Jstr)
            END DO
          ELSE
            DO i=IstrUm1,Iendp1
              LapU(i,Jstr-1)=0.0_r8
              tl_LapU(i,Jstr-1)=0.0_r8
            END DO
          END IF
          IF (tl_LBC(isouth,isVbar,ng)%closed) THEN
            DO i=Istrm1,Iendp1
              LapV(i,JstrV-1)=0.0_r8
              tl_LapV(i,JstrV-1)=0.0_r8
            END DO
          ELSE
            DO i=Istrm1,Iendp1
              LapV(i,JstrV-1)=LapV(i,JstrV)
              tl_LapV(i,JstrV-1)=tl_LapV(i,JstrV)
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (tl_LBC(inorth,isUbar,ng)%closed) THEN
            DO i=IstrUm1,Iendp1
              LapU(i,Jend+1)=gamma2(ng)*LapU(i,Jend)
              tl_LapU(i,Jend+1)=gamma2(ng)*tl_LapU(i,Jend)
            END DO
          ELSE
            DO i=IstrUm1,Iendp1
              LapU(i,Jend+1)=0.0_r8
              tl_LapU(i,Jend+1)=0.0_r8
            END DO
          END IF
          IF (tl_LBC(inorth,isVbar,ng)%closed) THEN
            DO i=Istrm1,Iendp1
              LapV(i,Jend+1)=0.0_r8
              tl_LapV(i,Jend+1)=0.0_r8
            END DO
          ELSE
            DO i=Istrm1,Iendp1
              LapV(i,Jend+1)=LapV(i,Jend)
              tl_LapV(i,Jend+1)=tl_LapV(i,Jend)
            END DO
          END IF
        END IF
      END IF
!
      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          LapU(Istr  ,Jstr-1)=0.5_r8*(LapU(Istr+1,Jstr-1)+              &
     &                                LapU(Istr  ,Jstr  ))
          tl_LapU(Istr  ,Jstr-1)=0.5_r8*(tl_LapU(Istr+1,Jstr-1)+        &
     &                                   tl_LapU(Istr  ,Jstr  ))
          LapV(Istr-1,Jstr  )=0.5_r8*(LapV(Istr-1,Jstr+1)+              &
     &                                LapV(Istr  ,Jstr  ))
          tl_LapV(Istr-1,Jstr  )=0.5_r8*(tl_LapV(Istr-1,Jstr+1)+        &
     &                                   tl_LapV(Istr  ,Jstr  ))
        END IF
      END IF

      IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          LapU(Iend+1,Jstr-1)=0.5_r8*(LapU(Iend  ,Jstr-1)+              &
     &                                LapU(Iend+1,Jstr  ))
          tl_LapU(Iend+1,Jstr-1)=0.5_r8*(tl_LapU(Iend  ,Jstr-1)+        &
     &                                   tl_LapU(Iend+1,Jstr  ))
          LapV(Iend+1,Jstr  )=0.5_r8*(LapV(Iend  ,Jstr  )+              &
     &                                LapV(Iend+1,Jstr+1))
          tl_LapV(Iend+1,Jstr  )=0.5_r8*(tl_LapV(Iend  ,Jstr  )+        &
     &                                   tl_LapV(Iend+1,Jstr+1))
        END IF
      END IF

      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          LapU(Istr  ,Jend+1)=0.5_r8*(LapU(Istr+1,Jend+1)+              &
     &                                LapU(Istr  ,Jend  ))
          tl_LapU(Istr  ,Jend+1)=0.5_r8*(tl_LapU(Istr+1,Jend+1)+        &
     &                                   tl_LapU(Istr  ,Jend  ))
          LapV(Istr-1,Jend+1)=0.5_r8*(LapV(Istr  ,Jend+1)+              &
     &                                LapV(Istr-1,Jend  ))
          tl_LapV(Istr-1,Jend+1)=0.5_r8*(tl_LapV(Istr  ,Jend+1)+        &
     &                                   tl_LapV(Istr-1,Jend  ))
        END IF
      END IF

      IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.          &
     &          CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          LapU(Iend+1,Jend+1)=0.5_r8*(LapU(Iend  ,Jend+1)+              &
     &                                LapU(Iend+1,Jend  ))
          tl_LapU(Iend+1,Jend+1)=0.5_r8*(tl_LapU(Iend  ,Jend+1)+        &
     &                                   tl_LapU(Iend+1,Jend  ))
          LapV(Iend+1,Jend+1)=0.5_r8*(LapV(Iend  ,Jend+1)+              &
     &                                LapV(Iend+1,Jend  ))
          tl_LapV(Iend+1,Jend+1)=0.5_r8*(tl_LapV(Iend  ,Jend+1)+        &
     &                                   tl_LapV(Iend+1,Jend  ))
        END IF
      END IF
!
!  Compute flux-components of the horizontal divergence of the
!  biharmonic stress tensor (m4/s2) in XI- and ETA-directions.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!>        cff=visc4_r(i,j)*Drhs(i,j)*0.5_r8*                            &
!>   &        (pmon_r(i,j)*                                             &
!>   &         ((pn(i  ,j)+pn(i+1,j))*LapU(i+1,j)-                      &
!>   &          (pn(i-1,j)+pn(i  ,j))*LapU(i  ,j))-                     &
!>   &         pnom_r(i,j)*                                             &
!>   &         ((pm(i,j  )+pm(i,j+1))*LapV(i,j+1)-                      &
!>   &          (pm(i,j-1)+pm(i,j  ))*LapV(i,j  )))
!>
          tl_cff=visc4_r(i,j)*0.5_r8*                                   &
     &           (tl_Drhs(i,j)*                                         &
     &            (pmon_r(i,j)*                                         &
     &             ((pn(i  ,j)+pn(i+1,j))*LapU(i+1,j)-                  &
     &              (pn(i-1,j)+pn(i  ,j))*LapU(i  ,j))-                 &
     &             pnom_r(i,j)*                                         &
     &             ((pm(i,j  )+pm(i,j+1))*LapV(i,j+1)-                  &
     &              (pm(i,j-1)+pm(i,j  ))*LapV(i,j  )))+                &
     &            Drhs(i,j)*                                            &
     &            (pmon_r(i,j)*                                         &
     &             ((pn(i  ,j)+pn(i+1,j))*tl_LapU(i+1,j)-               &
     &              (pn(i-1,j)+pn(i  ,j))*tl_LapU(i  ,j))-              &
     &             pnom_r(i,j)*                                         &
     &             ((pm(i,j  )+pm(i,j+1))*tl_LapV(i,j+1)-               &
     &              (pm(i,j-1)+pm(i,j  ))*tl_LapV(i,j  ))))
!>        UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
!>
          tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
!>        VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
!>
          tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
        END DO
      END DO
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
!>        cff=visc4_p(i,j)*Drhs_p(i,j)*0.5_r8*                          &
!>   &        (pmon_p(i,j)*                                             &
!>   &         ((pn(i  ,j-1)+pn(i  ,j))*LapV(i  ,j)-                    &
!>   &          (pn(i-1,j-1)+pn(i-1,j))*LapV(i-1,j))+                   &
!>   &         pnom_p(i,j)*                                             &
!>   &         ((pm(i-1,j  )+pm(i,j  ))*LapU(i,j  )-                    &
!>   &          (pm(i-1,j-1)+pm(i,j-1))*LapU(i,j-1)))
!>
          tl_cff=visc4_p(i,j)*0.5_r8*                                   &
     &           (tl_Drhs_p(i,j)*                                       &
     &            (pmon_p(i,j)*                                         &
     &             ((pn(i  ,j-1)+pn(i  ,j))*LapV(i  ,j)-                &
     &              (pn(i-1,j-1)+pn(i-1,j))*LapV(i-1,j))+               &
     &             pnom_p(i,j)*                                         &
     &             ((pm(i-1,j  )+pm(i,j  ))*LapU(i,j  )-                &
     &              (pm(i-1,j-1)+pm(i,j-1))*LapU(i,j-1)))+              &
     &            Drhs_p(i,j)*                                          &
     &            (pmon_p(i,j)*                                         &
     &             ((pn(i  ,j-1)+pn(i  ,j))*tl_LapV(i  ,j)-             &
     &              (pn(i-1,j-1)+pn(i-1,j))*tl_LapV(i-1,j))+            &
     &             pnom_p(i,j)*                                         &
     &             ((pm(i-1,j  )+pm(i,j  ))*tl_LapU(i,j  )-             &
     &              (pm(i-1,j-1)+pm(i,j-1))*tl_LapU(i,j-1))))
#  ifdef MASKING
!>        cff=cff*pmask(i,j)
!>
          tl_cff=tl_cff*pmask(i,j)
#  endif
!>        UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
!>
          tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
!>        VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
!>
          tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
        END DO
      END DO
!
!  Add in biharmonic viscosity.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*(UFx(i,j  )-UFx(i-1,j))
!>
          tl_cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*                           &
     &            (tl_UFx(i,j  )-tl_UFx(i-1,j))
!>        cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i  ,j))

          tl_cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*                           &
     &            (UFe(i,j+1)-UFe(i  ,j))
!>        fac=cff1+cff2
!>
          tl_fac=tl_cff1+tl_cff2
!>        rhs_ubar(i,j)=rhs_ubar(i,j)+fac
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2hvis)=fac
!!        DiaU2rhs(i,j,M2xvis)=cff1
!!        DiaU2rhs(i,j,M2yvis)=cff2
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j  ))
!>
          tl_cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*                           &
     &            (tl_VFx(i+1,j)-tl_VFx(i,j  ))
!>        cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*(VFe(i  ,j)-VFe(i,j-1))
!>
          tl_cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*                           &
     &            (tl_VFe(i  ,j)-tl_VFe(i,j-1))
!>        fac=cff1-cff2
!>
          tl_fac=tl_cff1-tl_cff2
!>        rhs_vbar(i,j)=rhs_vbar(i,j)+fac
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2hvis)=fac
!!        DiaV2rhs(i,j,M2xvis)= cff1
!!        DiaV2rhs(i,j,M2yvis)=-cff2
#  endif
        END DO
      END DO
# endif
# if defined NEARSHORE_MELLOR && \
    (!defined SOLVE3D         || defined DIAGNOSTICS_UV)
!
!-----------------------------------------------------------------------
!  Add in radiation stress terms.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        cff1=rustr2d(i,j)*om_u(i,j)*on_u(i,j)
!>
          tl_cff1=tl_rustr2d(i,j)*om_u(i,j)*on_u(i,j)
!>        cff2=rulag2d(i,j)
!>
          tl_cff2=tl_rulag2d(i,j)
#  ifndef SOLVE3D
!>        rhs_ubar(i,j)=rhs_ubar(i,j)-cff1-cff2
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_cff1-tl_cff2
#  endif
#  ifdef DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2hrad)=-cff1
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        cff1=rvstr2d(i,j)*om_v(i,j)*on_v(i,j)
!>
          tl_cff1=tl_rvstr2d(i,j)*om_v(i,j)*on_v(i,j)
!>        cff2=rvlag2d(i,j)
!>
          tl_cff2=tl_rvlag2d(i,j)
#  ifndef SOLVE3D
!>        rhs_vbar(i,j)=rhs_vbar(i,j)-cff1-cff2
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_cff1-tl_cff2
#  endif
#  ifdef DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2hrad)=-cff1
#  endif
        END DO
      END DO
# endif
# ifndef SOLVE3D
!
!-----------------------------------------------------------------------
!  Add in bottom stress.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        fac=bustr(i,j)*om_u(i,j)*on_u(i,j)
!>
          tl_fac=tl_bustr(i,j)*om_u(i,j)*on_u(i,j)
!>        rhs_ubar(i,j)=rhs_ubar(i,j)-fac
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
#  ifdef DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2bstr)=-fac
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        fac=bvstr(i,j)*om_v(i,j)*on_v(i,j)
!>
          tl_fac=tl_bvstr(i,j)*om_v(i,j)*on_v(i,j)
!>        rhs_vbar(i,j)=rhs_vbar(i,j)-fac
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
#  ifdef DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2bstr)=-fac
#  endif
        END DO
      END DO
# else
#  ifdef DIAGNOSTICS_UV
!!
!!  Initialize the stress term if no bottom friction is defined.
!!
!!    DO j=Jstr,Jend
!!      DO i=IstrU,Iend
!!        DiaU2rhs(i,j,M2bstr)=0.0_r8
!!      END DO
!!    END DO
!!    DO j=JstrV,Jend
!!      DO i=Istr,Iend
!!        DiaV2rhs(i,j,M2bstr)=0.0_r8
!!      END DO
!!    END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Add in nudging of 2D momentum climatology.
!-----------------------------------------------------------------------
!
      IF (LnudgeM2CLM(ng)) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.25_r8*(CLIMA(ng)%M2nudgcof(i-1,j)+                    &
     &                   CLIMA(ng)%M2nudgcof(i  ,j))*                   &
     &          om_u(i,j)*on_u(i,j)
!>          rhs_ubar(i,j)=rhs_ubar(i,j)+                                &
!>   &                    cff*(Drhs(i-1,j)+Drhs(i,j))*                  &
!>   &                        (CLIMA(ng)%ubarclm(i,j)-                  &
!>   &                         ubar(i,j,krhs))
!>
            tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                          &
     &                       cff*((Drhs(i-1,j)+Drhs(i,j))*              &
     &                            (-tl_ubar(i,j,krhs))+                 &
     &                            (tl_Drhs(i-1,j)+tl_Drhs(i,j))*        &
     &                            (CLIMA(ng)%ubarclm(i,j)-              &
     &                             ubar(i,j,krhs)))
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.25_r8*(CLIMA(ng)%M2nudgcof(i,j-1)+                    &
     &                   CLIMA(ng)%M2nudgcof(i,j  ))*                   &
     &          om_v(i,j)*on_v(i,j)
!>          rhs_vbar(i,j)=rhs_vbar(i,j)+                                &
!>   &                    cff*(Drhs(i,j-1)+Drhs(i,j))*                  &
!>   &                        (CLIMA(ng)%vbarclm(i,j)-                  &
!>   &                         vbar(i,j,krhs))
!>
            tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                          &
     &                       cff*((Drhs(i,j-1)+Drhs(i,j))*              &
     &                            (-tl_vbar(i,j,krhs))+                 &
     &                            (tl_Drhs(i,j-1)+tl_Drhs(i,j))*        &
     &                            (CLIMA(ng)%vbarclm(i,j)-              &
     &                             vbar(i,j,krhs)))
          END DO
        END DO
      END IF

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Coupling between 2D and 3D equations.
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
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
        IF (iic(ng).eq.ntfirst(ng).and.SOinitial(ng)) THEN
# else
        IF (iic(ng).eq.ntfirst(ng)) THEN
# endif
          DO j=Jstr,Jend
            DO i=IstrU,Iend
!>            rufrc(i,j)=rufrc(i,j)-rhs_ubar(i,j)
!>
              tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
!>            rhs_ubar(i,j)=rhs_ubar(i,j)+rufrc(i,j)
!>
              tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_rufrc(i,j)
!>            ru(i,j,0,nstp)=rufrc(i,j)
!>
              tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
#  ifdef DIAGNOSTICS_UV
!!            DO idiag=1,M2pgrd
!!              DiaRUfrc(i,j,3,idiag)=DiaRUfrc(i,j,3,idiag)-            &
!!   &                                DiaU2rhs(i,j,idiag)
!!              DiaU2rhs(i,j,idiag)=DiaU2rhs(i,j,idiag)+                &
!!   &                              DiaRUfrc(i,j,3,idiag)
!!              DiaRUfrc(i,j,nstp,idiag)=DiaRUfrc(i,j,3,idiag)
!!            END DO
!!            DiaU2rhs(i,j,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DiaU2rhs(i,j,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
#  endif
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
!>            rvfrc(i,j)=rvfrc(i,j)-rhs_vbar(i,j)
!>
              tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
!>            rhs_vbar(i,j)=rhs_vbar(i,j)+rvfrc(i,j)
!>
              tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
!>            rv(i,j,0,nstp)=rvfrc(i,j)
!>
              tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
#  ifdef DIAGNOSTICS_UV
!!            DO idiag=1,M2pgrd
!!              DiaRVfrc(i,j,3,idiag)=DiaRVfrc(i,j,3,idiag)-            &
!!   &                                DiaV2rhs(i,j,idiag)
!!              DiaV2rhs(i,j,idiag)=DiaV2rhs(i,j,idiag)+                &
!!   &                              DiaRVfrc(i,j,3,idiag)
!!              DiaRVfrc(i,j,nstp,idiag)=DiaRVfrc(i,j,3,idiag)
!!            END DO
!!            DiaV2rhs(i,j,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DiaV2rhs(i,j,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
#  endif
            END DO
          END DO
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1).and.SOinitial(ng)) THEN
# else
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
# endif
          DO j=Jstr,Jend
            DO i=IstrU,Iend
!>            rufrc(i,j)=rufrc(i,j)-rhs_ubar(i,j)
!>
              tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
!>            rhs_ubar(i,j)=rhs_ubar(i,j)+                              &
!>   &                      1.5_r8*rufrc(i,j)-0.5_r8*ru(i,j,0,nnew)
!>
              tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                        &
     &                         1.5_r8*tl_rufrc(i,j)-                    &
     &                         0.5_r8*tl_ru(i,j,0,nnew)
!>            ru(i,j,0,nstp)=rufrc(i,j)
!>
              tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
#  ifdef DIAGNOSTICS_UV
!!            DO idiag=1,M2pgrd
!!              DiaRUfrc(i,j,3,idiag)=DiaRUfrc(i,j,3,idiag)-            &
!!   &                                DiaU2rhs(i,j,idiag)
!!              DiaU2rhs(i,j,idiag)=DiaU2rhs(i,j,idiag)+                &
!!   &                              1.5_r8*DiaRUfrc(i,j,3,idiag)-       &
!!   &                              0.5_r8*DiaRUfrc(i,j,nnew,idiag)
!!              DiaRUfrc(i,j,nstp,idiag)=DiaRUfrc(i,j,3,idiag)
!!            END DO
!!            DiaU2rhs(i,j,M2sstr)=1.5_r8*DiaRUfrc(i,j,3,M2sstr)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2sstr)
!!            DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DiaU2rhs(i,j,M2bstr)=1.5_r8*DiaRUfrc(i,j,3,M2bstr)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2bstr)
!!            DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
#  endif
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
!>            rvfrc(i,j)=rvfrc(i,j)-rhs_vbar(i,j)
!>
              tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
!>            rhs_vbar(i,j)=rhs_vbar(i,j)+                              &
!>   &                      1.5_r8*rvfrc(i,j)-0.5_r8*rv(i,j,0,nnew)
!>
              tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                        &
     &                         1.5_r8*tl_rvfrc(i,j)-                    &
     &                         0.5_r8*tl_rv(i,j,0,nnew)
!>            rv(i,j,0,nstp)=rvfrc(i,j)
!>
              tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
#  ifdef DIAGNOSTICS_UV
!!            DO idiag=1,M2pgrd
!!              DiaRVfrc(i,j,3,idiag)=DiaRVfrc(i,j,3,idiag)-            &
!!   &                                DiaV2rhs(i,j,idiag)
!!              DiaV2rhs(i,j,idiag)=DiaV2rhs(i,j,idiag)+                &
!!   &                              1.5_r8*DiaRVfrc(i,j,3,idiag)-       &
!!   &                              0.5_r8*DiaRVfrc(i,j,nnew,idiag)
!!              DiaRVfrc(i,j,nstp,idiag)=DiaRVfrc(i,j,3,idiag)
!!            END DO
!!            DiaV2rhs(i,j,M2sstr)=1.5_r8*DiaRVfrc(i,j,3,M2sstr)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2sstr)
!!            DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DiaV2rhs(i,j,M2bstr)=1.5_r8*DiaRVfrc(i,j,3,M2bstr)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2bstr)
!!            DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
#  endif
            END DO
          END DO
        ELSE
          cff1=23.0_r8/12.0_r8
          cff2=16.0_r8/12.0_r8
          cff3= 5.0_r8/12.0_r8
          DO j=Jstr,Jend
            DO i=IstrU,Iend
!>            rufrc(i,j)=rufrc(i,j)-rhs_ubar(i,j)
!>
              tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
!>            rhs_ubar(i,j)=rhs_ubar(i,j)+                              &
!>   &                      cff1*rufrc(i,j)-                            &
!>   &                      cff2*ru(i,j,0,nnew)+                        &
!>   &                      cff3*ru(i,j,0,nstp)
!>
              tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                        &
     &                         cff1*tl_rufrc(i,j)-                      &
     &                         cff2*tl_ru(i,j,0,nnew)+                  &
     &                         cff3*tl_ru(i,j,0,nstp)
!>            ru(i,j,0,nstp)=rufrc(i,j)
!>
              tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
#  ifdef DIAGNOSTICS_UV
!!            DO idiag=1,M2pgrd
!!              DiaRUfrc(i,j,3,idiag)=DiaRUfrc(i,j,3,idiag)-            &
!!   &                                DiaU2rhs(i,j,idiag)
!!              DiaU2rhs(i,j,idiag)=DiaU2rhs(i,j,idiag)+                &
!!   &                              cff1*DiaRUfrc(i,j,3,idiag)-         &
!!   &                              cff2*DiaRUfrc(i,j,nnew,idiag)+      &
!!   &                              cff3*DiaRUfrc(i,j,nstp,idiag)
!!              DiaRUfrc(i,j,nstp,idiag)=DiaRUfrc(i,j,3,idiag)
!!            END DO
!!            DiaU2rhs(i,j,M2sstr)=cff1*DiaRUfrc(i,j,3,M2sstr)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2sstr)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2sstr)
!!            DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DiaU2rhs(i,j,M2bstr)=cff1*DiaRUfrc(i,j,3,M2bstr)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2bstr)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2bstr)
!!            DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
#  endif
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
!>            rvfrc(i,j)=rvfrc(i,j)-rhs_vbar(i,j)
!>
              tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
!>            rhs_vbar(i,j)=rhs_vbar(i,j)+                              &
!>   &                      cff1*rvfrc(i,j)-                            &
!>   &                      cff2*rv(i,j,0,nnew)+                        &
!>   &                      cff3*rv(i,j,0,nstp)
!>
              tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                        &
     &                         cff1*tl_rvfrc(i,j)-                      &
     &                         cff2*tl_rv(i,j,0,nnew)+                  &
     &                         cff3*tl_rv(i,j,0,nstp)
!>            rv(i,j,0,nstp)=rvfrc(i,j)
!>
              tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
#  ifdef DIAGNOSTICS_UV
!!            DO idiag=1,M2pgrd
!!              DiaRVfrc(i,j,3,idiag)=DiaRVfrc(i,j,3,idiag)-            &
!!   &                                DiaV2rhs(i,j,idiag)
!!              DiaV2rhs(i,j,idiag)=DiaV2rhs(i,j,idiag)+                &
!!   &                              cff1*DiaRVfrc(i,j,3,idiag)-         &
!!   &                              cff2*DiaRVfrc(i,j,nnew,idiag)+      &
!!   &                              cff3*DiaRVfrc(i,j,nstp,idiag)
!!              DiaRVfrc(i,j,nstp,idiag)=DiaRVfrc(i,j,3,idiag)
!!            END DO
!!            DiaV2rhs(i,j,M2sstr)=cff1*DiaRVfrc(i,j,3,M2sstr)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2sstr)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2sstr)
!!            DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DiaV2rhs(i,j,M2bstr)=cff1*DiaRVfrc(i,j,3,M2bstr)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2bstr)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2bstr)
!!            DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
#  endif
            END DO
          END DO
        END IF
      ELSE
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!>          rhs_ubar(i,j)=rhs_ubar(i,j)+rufrc(i,j)
!>
            tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_rufrc(i,j)
#  ifdef DIAGNOSTICS_UV
!!          DO idiag=1,M2pgrd
!!            DiaU2rhs(i,j,idiag)=DiaU2rhs(i,j,idiag)+                  &
!!   &                            DiaRUfrc(i,j,3,idiag)
!!          END DO
!!          DiaU2rhs(i,j,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!          DiaU2rhs(i,j,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
!>          rhs_vbar(i,j)=rhs_vbar(i,j)+rvfrc(i,j)
!>
            tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
#  ifdef DIAGNOSTICS_UV
!!          DO idiag=1,M2pgrd
!!            DiaV2rhs(i,j,idiag)=DiaV2rhs(i,j,idiag)+                  &
!!   &                            DiaRVfrc(i,j,3,idiag)
!!          END DO
!!          DiaV2rhs(i,j,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!          DiaV2rhs(i,j,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
#  endif
          END DO
        END DO
      END IF
# else
!>
!>----------------------------------------------------------------------
!>  Add in surface momentum stress.
!>----------------------------------------------------------------------
!>
!>    DO j=Jstr,Jend
!>      DO i=IstrU,Iend
!>        fac=sustr(i,j)*om_u(i,j)*on_u(i,j)
!>        rhs_ubar(i,j)=rhs_ubar(i,j)+fac
#  ifdef DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2sstr)=fac
#  endif
!>      END DO
!>    END DO
!>    DO j=JstrV,Jend
!>      DO i=Istr,Iend
!>        fac=svstr(i,j)*om_v(i,j)*on_v(i,j)
!>        rhs_vbar(i,j)=rhs_vbar(i,j)+fac
#  ifdef DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2sstr)=fac
#  endif
!>      END DO
!>    END DO
# endif
!
!=======================================================================
!  Time step 2D momentum equations.
!=======================================================================
!
!  Compute total water column depth.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
          tl_Dstp(i,j)=tl_zeta(i,j,kstp)+tl_h(i,j)
        END DO
      END DO
!
!  During the first time-step, the predictor step is Forward-Euler
!  and the corrector step is Backward-Euler. Otherwise, the predictor
!  step is Leap-frog and the corrector step is Adams-Moulton.
# ifdef WET_DRY_NOT_YET
!  HGA:  We need to think more about TLM of the wet/dry mask arrays
!        since they are time-dependent.
# endif
!
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE && \
    !defined SOLVE3D
      IF (FIRST_2D_STEP.and.SOinitial(ng)) THEN
# else
      IF (FIRST_2D_STEP) THEN
# endif
        cff1=0.5_r8*dtfast(ng)
# ifdef WET_DRY_NOT_YET
        cff2=1.0_r8/cff1
# endif
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i-1,j))
!>          ubar(i,j,knew)=(ubar(i,j,kstp)*                             &
!>   &                      (Dstp(i,j)+Dstp(i-1,j))+                    &
!>   &                      cff*cff1*rhs_ubar(i,j))*fac
!>
            tl_ubar(i,j,knew)=(tl_ubar(i,j,kstp)*                       &
     &                         (Dstp(i,j)+Dstp(i-1,j))+                 &
     &                         ubar(i,j,kstp)*                          &
     &                         (tl_Dstp(i,j)+tl_Dstp(i-1,j))+           &
     &                         cff*cff1*tl_rhs_ubar(i,j))*fac+          &
     &                        (ubar(i,j,kstp)*                          &
     &                         (Dstp(i,j)+Dstp(i-1,j))+                 &
     &                         cff*cff1*rhs_ubar(i,j))*tl_fac
# ifdef MASKING
!>          ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
!>
            tl_ubar(i,j,knew)=tl_ubar(i,j,knew)*umask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>          cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
!>          cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
!>          cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>          ubar(i,j,knew)=ubar(i,j,knew)*cff7
!>
!>  HGA: TLM code needed here.
!>
            fac1=cff2/cff
!>          rhs_ubar(i,j)=(ubar(i,j,knew)*(Dnew(i,j)+Dnew(i-1,j))-      &
!>   &                     ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j)))*     &
!>   &                    fac1
!>
            tl_rhs_ubar(i,j)=(tl_ubar(i,j,knew)*                        &
     &                        (Dnew(i,j)+Dnew(i-1,j))+                  &
     &                        ubar(i,j,knew)*                           &
     &                        (tl_Dnew(i,j)+tl_Dnew(i-1,j))-            &
     &                        tl_ubar(i,j,kstp)*                        &
     &                        (Dstp(i,j)+Dstp(i-1,j))-                  &
     &                        ubar(i,j,kstp)*                           &
     &                        (tl_Dstp(i,j)+tl_Dstp(i-1,j)))*fac1
# endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i,j-1))
!>          vbar(i,j,knew)=(vbar(i,j,kstp)*                             &
!>   &                      (Dstp(i,j)+Dstp(i,j-1))+                    &
!>   &                      cff*cff1*rhs_vbar(i,j))*fac
!>
            tl_vbar(i,j,knew)=(tl_vbar(i,j,kstp)*                       &
     &                         (Dstp(i,j)+Dstp(i,j-1))+                 &
     &                         vbar(i,j,kstp)*                          &
     &                         (tl_Dstp(i,j)+tl_Dstp(i,j-1))+           &
     &                         cff*cff1*tl_rhs_vbar(i,j))*fac+          &
     &                        (vbar(i,j,kstp)*                          &
     &                         (Dstp(i,j)+Dstp(i,j-1))+                 &
     &                         cff*cff1*rhs_vbar(i,j))*tl_fac
# ifdef MASKING
!>          vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
!>
            tl_vbar(i,j,knew)=tl_vbar(i,j,knew)*vmask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>          cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
!>          cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
!>          cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>          vbar(i,j,knew)=vbar(i,j,knew)*cff7
!>
!>  HGA: TLM code needed here.
!>
            fac1=cff2/cff
!>          rhs_vbar(i,j)=(vbar(i,j,knew)*(Dnew(i,j)+Dnew(i,j-1))-      &
!>   &                     vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1)))*     &
!>   &                    fac1
!>
            tl_rhs_vbar(i,j)=(tl_vbar(i,j,knew)*                        &
     &                        (Dnew(i,j)+Dnew(i,j-1))+                  &
     &                        vbar(i,j,knew)*                           &
     &                        (tl_Dnew(i,j)+tl_Dnew(i,j-1))-            &
     &                        tl_vbar(i,j,kstp)*                        &
     &                        (Dstp(i,j)+Dstp(i,j-1))-                  &
     &                        vbar(i,j,kstp)*                           &
     &                        (tl_Dstp(i,j)+tl_Dstp(i,j-1)))*fac1
# endif
          END DO
        END DO
      ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
        cff1=dtfast(ng)
# ifdef WET_DRY_NOT_YET
        cff2=1.0_r8/cff1
# endif
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i-1,j))
!>          ubar(i,j,knew)=(ubar(i,j,kstp)*                             &
!>   &                      (Dstp(i,j)+Dstp(i-1,j))+                    &
!>   &                      cff*cff1*rhs_ubar(i,j))*fac
!>
            tl_ubar(i,j,knew)=(tl_ubar(i,j,kstp)*                       &
     &                         (Dstp(i,j)+Dstp(i-1,j))+                 &
     &                         ubar(i,j,kstp)*                          &
     &                         (tl_Dstp(i,j)+tl_Dstp(i-1,j))+           &
     &                         cff*cff1*tl_rhs_ubar(i,j))*fac+          &
     &                        (ubar(i,j,kstp)*                          &
     &                         (Dstp(i,j)+Dstp(i-1,j))+                 &
     &                         cff*cff1*rhs_ubar(i,j))*tl_fac
# ifdef MASKING
!>          ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
!>
            tl_ubar(i,j,knew)=tl_ubar(i,j,knew)*umask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>          cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
!>          cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
!>          cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>          ubar(i,j,knew)=ubar(i,j,knew)*cff7
!>
!>  HGA: TLM code needed here.
!>
            fac1=cff2/cff
!>          rhs_ubar(i,j)=(ubar(i,j,knew)*(Dnew(i,j)+Dnew(i-1,j))-      &
!>   &                     ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j)))*     &
!>   &                    fac1
!>
            tl_rhs_ubar(i,j)=(tl_ubar(i,j,knew)*                        &
     &                        (Dnew(i,j)+Dnew(i-1,j))+                  &
     &                        ubar(i,j,knew)*                           &
     &                        (tl_Dnew(i,j)+tl_Dnew(i-1,j))-            &
     &                        tl_ubar(i,j,kstp)*                        &
     &                        (Dstp(i,j)+Dstp(i-1,j))-                  &
     &                        ubar(i,j,kstp)*                           &
     &                        (tl_Dstp(i,j)+tl_Dstp(i-1,j)))*fac1
# endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i,j-1))
!>          vbar(i,j,knew)=(vbar(i,j,kstp)*                             &
!>   &                      (Dstp(i,j)+Dstp(i,j-1))+                    &
!>   &                      cff*cff1*rhs_vbar(i,j))*fac
!>
            tl_vbar(i,j,knew)=(tl_vbar(i,j,kstp)*                       &
     &                         (Dstp(i,j)+Dstp(i,j-1))+                 &
     &                         vbar(i,j,kstp)*                          &
     &                         (tl_Dstp(i,j)+tl_Dstp(i,j-1))+           &
     &                         cff*cff1*tl_rhs_vbar(i,j))*fac+          &
     &                        (vbar(i,j,kstp)*                          &
     &                         (Dstp(i,j)+Dstp(i,j-1))+                 &
     &                         cff*cff1*rhs_vbar(i,j))*tl_fac
# ifdef MASKING
!>          vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
!>
            tl_vbar(i,j,knew)=tl_vbar(i,j,knew)*vmask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>          cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
!>          cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
!>          cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>          vbar(i,j,knew)=vbar(i,j,knew)*cff7
!>
!>  HGA: TLM code needed here.
!>
            fac1=cff2/cff
!>          rhs_vbar(i,j)=(vbar(i,j,knew)*(Dnew(i,j)+Dnew(i,j-1))-      &
!>   &                     vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1)))*     &
!>   &                    fac1
!>
            tl_rhs_vbar(i,j)=(tl_vbar(i,j,knew)*                        &
     &                        (Dnew(i,j)+Dnew(i,j-1))+                  &
     &                        vbar(i,j,knew)*                           &
     &                        (tl_Dnew(i,j)+tl_Dnew(i,j-1))-            &
     &                        tl_vbar(i,j,kstp)*                        &
     &                        (Dstp(i,j)+Dstp(i,j-1))-                  &
     &                        vbar(i,j,kstp)*                           &
     &                        (tl_Dstp(i,j)+tl_Dstp(i,j-1)))*fac1
# endif
          END DO
        END DO
      ELSE IF (CORRECTOR_2D_STEP) THEN
        cff1=0.5_r8*dtfast(ng)*5.0_r8/12.0_r8
        cff2=0.5_r8*dtfast(ng)*8.0_r8/12.0_r8
        cff3=0.5_r8*dtfast(ng)*1.0_r8/12.0_r8
# ifdef WET_DRY_NOT_YET
        cff4=1.0_r8/cff1
# endif
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i-1,j))
!>          ubar(i,j,knew)=(ubar(i,j,kstp)*                             &
!>   &                      (Dstp(i,j)+Dstp(i-1,j))+                    &
!>   &                      cff*(cff1*rhs_ubar(i,j)+                    &
!>   &                           cff2*rubar(i,j,kstp)-                  &
!>   &                           cff3*rubar(i,j,ptsk)))*fac
!>
            tl_ubar(i,j,knew)=(tl_ubar(i,j,kstp)*                       &
     &                         (Dstp(i,j)+Dstp(i-1,j))+                 &
     &                         ubar(i,j,kstp)*                          &
     &                         (tl_Dstp(i,j)+tl_Dstp(i-1,j))+           &
     &                         cff*(cff1*tl_rhs_ubar(i,j)+              &
     &                              cff2*tl_rubar(i,j,kstp)-            &
     &                              cff3*tl_rubar(i,j,ptsk)))*fac+      &
     &                        (ubar(i,j,kstp)*                          &
     &                         (Dstp(i,j)+Dstp(i-1,j))+                 &
     &                         cff*(cff1*rhs_ubar(i,j)+                 &
     &                              cff2*rubar(i,j,kstp)-               &
     &                              cff3*rubar(i,j,ptsk)))*tl_fac
# ifdef MASKING
!>          ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
!>
            tl_ubar(i,j,knew)=tl_ubar(i,j,knew)*umask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>          cff5=ABS(ABS(umask_wet(i,j))-1.0_r8)
!>          cff6=0.5_r8+DSIGN(0.5_r8,ubar(i,j,knew))*umask_wet(i,j)
!>          cff7=0.5_r8*umask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>          ubar(i,j,knew)=ubar(i,j,knew)*cff7
!>
!>  HGA: TLM code needed here.
!>
            fac1=1.0_r8/cff
!>          rhs_ubar(i,j)=((ubar(i,j,knew)*(Dnew(i,j)+Dnew(i-1,j))-     &
!>   &                      ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j)))*    &
!>   &                     fac1-                                        &
!>   &                     cff2*rubar(i,j,kstp)+                        &
!>   &                     cff3*rubar(i,j,ptsk))*cff4
!>
            tl_rhs_ubar(i,j)=((tl_ubar(i,j,knew)*                       &
     &                         (Dnew(i,j)+Dnew(i-1,j))+                 &
     &                         ubar(i,j,knew)*                          &
     &                         (tl_Dnew(i,j)+tl_Dnew(i-1,j))-           &
     &                         tl_ubar(i,j,kstp)*                       &
     &                         (Dstp(i,j)+Dstp(i-1,j))-                 &
     &                         ubar(i,j,kstp)*                          &
     &                         (tl_Dstp(i,j)+tl_Dstp(i-1,j)))*fac1-     &
     &                        cff2*tl_rubar(i,j,kstp)+                  &
     &                        cff3*tl_rubar(i,j,ptsk))*cff4
# endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
            tl_fac=-fac*fac*(tl_Dnew(i,j)+tl_Dnew(i,j-1))
!>          vbar(i,j,knew)=(vbar(i,j,kstp)*                             &
!>   &                      (Dstp(i,j)+Dstp(i,j-1))+                    &
!>   &                      cff*(cff1*rhs_vbar(i,j)+                    &
!>   &                           cff2*rvbar(i,j,kstp)-                  &
!>   &                           cff3*rvbar(i,j,ptsk)))*fac
!>
            tl_vbar(i,j,knew)=(tl_vbar(i,j,kstp)*                       &
     &                         (Dstp(i,j)+Dstp(i,j-1))+                 &
     &                         vbar(i,j,kstp)*                          &
     &                         (tl_Dstp(i,j)+tl_Dstp(i,j-1))+           &
     &                         cff*(cff1*tl_rhs_vbar(i,j)+              &
     &                              cff2*tl_rvbar(i,j,kstp)-            &
     &                              cff3*tl_rvbar(i,j,ptsk)))*fac+      &
     &                        (vbar(i,j,kstp)*                          &
     &                         (Dstp(i,j)+Dstp(i,j-1))+                 &
     &                         cff*(cff1*rhs_vbar(i,j)+                 &
     &                              cff2*rvbar(i,j,kstp)-               &
     &                              cff3*rvbar(i,j,ptsk)))*tl_fac
# ifdef MASKING
!>          vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
!>
            tl_vbar(i,j,knew)=tl_vbar(i,j,knew)*vmask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>          cff5=ABS(ABS(vmask_wet(i,j))-1.0_r8)
!>          cff6=0.5_r8+DSIGN(0.5_r8,vbar(i,j,knew))*vmask_wet(i,j)
!>          cff7=0.5_r8*vmask_wet(i,j)*cff5+cff6*(1.0_r8-cff5)
!>          vbar(i,j,knew)=vbar(i,j,knew)*cff7
!>
!>  HGA: TLM code needed here.
!>
            fac1=1.0_r8/cff
!>          rhs_vbar(i,j)=((vbar(i,j,knew)*(Dnew(i,j)+Dnew(i,j-1))-     &
!>   &                      vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1)))*    &
!>   &                     fac1-                                        &
!>   &                     cff2*rvbar(i,j,kstp)+                        &
!>   &                     cff3*rvbar(i,j,ptsk))*cff4
!>
            tl_rhs_vbar(i,j)=((tl_vbar(i,j,knew)*                       &
     &                         (Dnew(i,j)+Dnew(i,j-1))+                 &
     &                         vbar(i,j,knew)*                          &
     &                         (tl_Dnew(i,j)+tl_Dnew(i,j-1))-           &
     &                         tl_vbar(i,j,kstp)*                       &
     &                         (Dstp(i,j)+Dstp(i,j-1))-                 &
     &                         vbar(i,j,kstp)*                          &
     &                         (tl_Dstp(i,j)+tl_Dstp(i,j-1)))*fac1-     &
     &                        cff2*tl_rvbar(i,j,kstp)+                  &
     &                        cff3*tl_rvbar(i,j,ptsk))*cff4
# endif
          END DO
        END DO
      END IF
# ifdef DIAGNOSTICS_UV
!!
!!-----------------------------------------------------------------------
!!  Time step 2D momentum diagnostic terms.
!!-----------------------------------------------------------------------

#  ifdef MASKING
!!
!!  Apply land/sea mask.
!!
!!    DO idiag=1,NDM2d-1
!!      DO j=Jstr,Jend
!!        DO i=IstrU,Iend
!!          DiaU2rhs(i,j,idiag)=DiaU2rhs(i,j,idiag)*umask(i,j)
!!        END DO
!!      END DO
!!      DO j=JstrV,Jend
!!        DO i=Istr,Iend
!!          DiaV2rhs(i,j,idiag)=DiaV2rhs(i,j,idiag)*vmask(i,j)
!!        END DO
!!      END DO
!!    END DO
#  endif
#  ifdef SOLVE3D
!!
!!  The arrays "DiaU2rhs" and "DiaV2rhs" contain the contributions of
!!  each of the 2D right-hand-side terms for the momentum equations.
!!
!!  These values are integrated, time-stepped and converted to mass flux
!!  units (m3 s-1) for coupling with the 3D diagnostic terms.
!!
!!    fac=weight(1,iif(ng),ng)
!!    IF (FIRST_2D_STEP.and.CORRECTOR_2D_STEP) THEN
!!      cff1=0.5_r8*dtfast(ng)
!!      DO idiag=1,NDM2d-1
!!        DO j=Jstr,Jend
!!          DO i=IstrU,Iend
!!            DiaU2int(i,j,idiag)=cff1*DiaU2rhs(i,j,idiag)
!!            DiaU2wrk(i,j,idiag)=DiaU2int(i,j,idiag)*                  &
!!   &                            (pm(i-1,j)+pm(i,j))*fac
!!          END DO
!!        END DO
!!        DO j=JstrV,Jend
!!          DO i=Istr,Iend
!!            DiaV2int(i,j,idiag)=cff1*DiaV2rhs(i,j,idiag)
!!            DiaV2wrk(i,j,idiag)=DiaV2int(i,j,idiag)*                  &
!!   &                            (pn(i,j)+pn(i,j-1))*fac
!!          END DO
!!        END DO
!!      END DO
!!    ELSE IF (CORRECTOR_2D_STEP) THEN
!!      cff1=0.5_r8*dtfast(ng)*5.0_r8/12.0_r8
!!      cff2=0.5_r8*dtfast(ng)*8.0_r8/12.0_r8
!!      cff3=0.5_r8*dtfast(ng)*1.0_r8/12.0_r8
!!      DO idiag=1,NDM2d-1
!!        DO j=Jstr,Jend
!!          DO i=IstrU,Iend
!!            DiaU2int(i,j,idiag)=DiaU2int(i,j,idiag)+                  &
!!   &                            (cff1*DiaU2rhs(i,j,idiag)+            &
!!   &                             cff2*DiaRUbar(i,j,kstp,idiag)-       &
!!   &                             cff3*DiaRUbar(i,j,ptsk,idiag))
!!            DiaU2wrk(i,j,idiag)=DiaU2wrk(i,j,idiag)+                  &
!!   &                            DiaU2int(i,j,idiag)*                  &
!!   &                            (pm(i-1,j)+pm(i,j))*fac
!!          END DO
!!        END DO
!!        DO j=JstrV,Jend
!!          DO i=Istr,Iend
!!            DiaV2int(i,j,idiag)=DiaV2int(i,j,idiag)+                  &
!!   &                            (cff1*DiaV2rhs(i,j,idiag)+            &
!!   &                             cff2*DiaRVbar(i,j,kstp,idiag)-       &
!!   &                             cff3*DiaRVbar(i,j,ptsk,idiag))
!!            DiaV2wrk(i,j,idiag)=DiaV2wrk(i,j,idiag)+                  &
!!   &                            DiaV2int(i,j,idiag)*                  &
!!   &                            (pn(i,j)+pn(i,j-1))*fac
!!          END DO
!!        END DO
!!      END DO
!!    END IF
#  else
!!
!!  Time-step the diagnostic terms.
!!
!!    IF (FIRST_2D_STEP.and.CORRECTOR_2D_STEP) THEN
!!      cff1=0.5_r8*dtfast(ng)
!!      DO idiag=1,NDM2d-1
!!        DO j=Jstr,Jend
!!          DO i=IstrU,Iend
!!            cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!            DiaU2wrk(i,j,idiag)=cff*cff1*DiaU2rhs(i,j,idiag)*fac
!!          END DO
!!        END DO
!!        DO j=JstrV,Jend
!!          DO i=Istr,Iend
!!            cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!            DiaV2wrk(i,j,idiag)=cff*cff1*DiaV2rhs(i,j,idiag)*fac
!!          END DO
!!        END DO
!!      END DO
!!      DO j=Jstr,Jend
!!        DO i=IstrU,Iend
!!          fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!          DiaU2wrk(i,j,M2rate)=ubar(i,j,knew)-ubar(i,j,kstp)*         &
!!   &                           (Dstp(i,j)+Dstp(i-1,j))*fac
!!        END DO
!!      END DO
!!      DO j=JstrV,Jend
!!        DO i=Istr,Iend
!!          fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!          DiaV2wrk(i,j,M2rate)=vbar(i,j,knew)-vbar(i,j,kstp)*         &
!!   &                           (Dstp(i,j)+Dstp(i,j-1))*fac
!!        END DO
!!      END DO
!!    ELSE IF (CORRECTOR_2D_STEP) THEN
!!      cff1=0.5_r8*dtfast(ng)*5.0_r8/12.0_r8
!!      cff2=0.5_r8*dtfast(ng)*8.0_r8/12.0_r8
!!      cff3=0.5_r8*dtfast(ng)*1.0_r8/12.0_r8
!!      DO idiag=1,NDM2d-1
!!        DO j=Jstr,Jend
!!          DO i=IstrU,Iend
!!            cff=(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!            DiaU2wrk(i,j,idiag)=cff*(cff1*DiaU2rhs(i,j,idiag)+        &
!!   &                                 cff2*DiaRUbar(i,j,kstp,idiag)-   &
!!   &                                 cff3*DiaRUbar(i,j,ptsk,idiag))*  &
!!   &                                fac
!!          END DO
!!        END DO
!!        DO j=JstrV,Jend
!!          DO i=Istr,Iend
!!            cff=(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!!            fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!            DiaV2wrk(i,j,idiag)=cff*(cff1*DiaV2rhs(i,j,idiag)+        &
!!   &                                 cff2*DiaRVbar(i,j,kstp,idiag)-   &
!!   &                                 cff3*DiaRVbar(i,j,ptsk,idiag))*
!!   &                                fac
!!          END DO
!!        END DO
!!      END DO
!!      DO j=Jstr,Jend
!!        DO i=IstrU,Iend
!!          fac=1.0_r8/(Dnew(i,j)+Dnew(i-1,j))
!!          DiaU2wrk(i,j,M2rate)=ubar(i,j,knew)-                        &
!!   &                           ubar(i,j,kstp)*                        &
!!   &                           (Dstp(i,j)+Dstp(i-1,j))*fac
!!        END DO
!!      END DO
!!      DO j=JstrV,Jend
!!        DO i=Istr,Iend
!!          fac=1.0_r8/(Dnew(i,j)+Dnew(i,j-1))
!!          DiaV2wrk(i,j,M2rate)=vbar(i,j,knew)-                        &
!!   &                           vbar(i,j,kstp)*                        &
!!   &                           (Dstp(i,j)+Dstp(i,j-1))*fac
!!        END DO
!!      END DO
!!    END IF
#  endif
# endif
!
!  If predictor step, load right-side-term into shared arrays for
!  future use during the subsequent corrector step.
!
      IF (PREDICTOR_2D_STEP(ng)) THEN
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!>          rubar(i,j,krhs)=rhs_ubar(i,j)
!>
            tl_rubar(i,j,krhs)=tl_rhs_ubar(i,j)
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
!>          rvbar(i,j,krhs)=rhs_vbar(i,j)
!>
            tl_rvbar(i,j,krhs)=tl_rhs_vbar(i,j)
          END DO
        END DO
# ifdef DIAGNOSTICS_UV
!!      DO idiag=1,NDM2d-1
!!        DO j=Jstr,Jend
!!          DO i=IstrU,Iend
!!            DiaRUbar(i,j,krhs,idiag)=DiaU2rhs(i,j,idiag)
!!          END DO
!!        END DO
!!        DO j=JstrV,Jend
!!          DO i=Istr,Iend
!!            DiaRVbar(i,j,krhs,idiag)=DiaV2rhs(i,j,idiag)
!!          END DO
!!        END DO
!!      END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions.
!-----------------------------------------------------------------------
!
!>    CALL u2dbc_tile (ng, tile,                                        &
!>   &                 LBi, UBi, LBj, UBj,                              &
!>   &                 IminS, ImaxS, JminS, JmaxS,                      &
!>   &                 krhs, kstp, knew,                                &
!>   &                 ubar, vbar, zeta)
!>
      CALL tl_u2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, knew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
!>    CALL v2dbc_tile (ng, tile,                                        &
!>   &                 LBi, UBi, LBj, UBj,                              &
!>   &                 IminS, ImaxS, JminS, JmaxS,                      &
!>   &                 krhs, kstp, knew,                                &
!>   &                 ubar, vbar, zeta)
!>
      CALL tl_v2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, knew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    tl_ubar, tl_vbar, tl_zeta)
!
!  Compute integral mass flux across open boundaries and adjust
!  for volume conservation.
!
      IF (ANY(tl_VolCons(:,ng))) THEN
        CALL tl_obc_flux_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         knew,                                    &
# ifdef MASKING
     &                         umask, vmask,                            &
# endif
     &                         h, tl_h, om_v, on_u,                     &
     &                         ubar, vbar, zeta,                        &
     &                         tl_ubar, tl_vbar, tl_zeta)
      END IF
!
!-----------------------------------------------------------------------
!  Apply momentum transport point sources (like river runoff), if any.
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
     &                    0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+            &
     &                            zeta(i  ,j,knew)+h(i  ,j)))
              tl_cff=-cff*cff*on_u(i,j)*                                &
     &               0.5_r8*(tl_zeta(i-1,j,knew)+tl_h(i-1,j)+           &
     &                       tl_zeta(i  ,j,knew)+tl_h(i  ,j))
!>            ubar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
!>
              tl_ubar(i,j,knew)=SOURCES(ng)%Qbar(is)*tl_cff
            ELSE
              cff=1.0_r8/(om_v(i,j)*                                    &
     &                    0.5_r8*(zeta(i,j-1,knew)+h(i,j-1)+            &
     &                            zeta(i,j  ,knew)+h(i,j  )))
              tl_cff=-cff*cff*om_v(i,j)*                                &
     &               0.5_r8*(tl_zeta(i,j-1,knew)+tl_h(i,j-1)+           &
     &                       tl_zeta(i,j  ,knew)+tl_h(i,j  ))
!>            vbar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
!>
              tl_vbar(i,j,knew)=SOURCES(ng)%Qbar(is)*tl_cff
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          ubar(:,:,knew))
!>
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_ubar(:,:,knew))
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          vbar(:,:,knew))
!>
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_vbar(:,:,knew))
      END IF

# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    ubar(:,:,knew),                               &
!>   &                    vbar(:,:,knew))
!>
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_ubar(:,:,knew),                            &
     &                    tl_vbar(:,:,knew))
# endif

      RETURN
      END SUBROUTINE tl_step2d_tile
#else
      SUBROUTINE tl_step2d
      END SUBROUTINE tl_step2d
#endif
