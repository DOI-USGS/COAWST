#undef DEBUG
#ifdef TANGENT
      SUBROUTINE tl_step2d (ng, tile)
!
!svn $Id: tl_step2d_LF_AM3.h 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine performs a fast (predictor or corrector) time-step     !
!  for the free-surface and 2D momentum  tangent linear equations.     !
!  The predictor step is  Leap-Frog  whereas the corrector step is     !
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
     &                     GRID(ng) % tl_bed_thick,                     &
#  endif
     &                     GRID(ng) % tl_Hz,                            &
     &                     GRID(ng) % tl_z_r,      GRID(ng) % tl_z_w,   &
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
# ifdef M2CLIMATOLOGY
     &                     CLIMA(ng) % ubarclm,    CLIMA(ng) % vbarclm, &
#  ifdef M2CLM_NUDGING
     &                     CLIMA(ng) % M2nudgcof,                       &
#  endif
# endif
# ifndef SOLVE3D
     &                     FORCES(ng) % tl_bustr,                       &
     &                     FORCES(ng) % tl_bvstr,                       &
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
     &                           tl_bed_thick,                          &
#  endif
     &                           tl_Hz, tl_z_r, tl_z_w,                 &
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
# ifdef M2CLIMATOLOGY
     &                           ubarclm, vbarclm,                      &
#  ifdef M2CLM_NUDGING
     &                           M2nudgcof,                             &
#  endif
# endif
# ifndef SOLVE3D
     &                           tl_bustr, tl_bvstr,                    &
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
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
!
# ifdef WET_DRY_NOT_YET
      USE bc_2d_mod
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
# ifdef OBC_VOLCONS
      USE obc_volcons_mod
      USE tl_obc_volcons_mod
# endif
# ifdef SOLVE3D
      USE tl_set_depth_mod, ONLY : tl_set_depth_tile
# endif
      USE tl_u2dbc_mod, ONLY : tl_u2dbc_tile
      USE tl_v2dbc_mod, ONLY : tl_v2dbc_tile
      USE tl_zetabc_mod, ONLY : tl_zetabc_tile
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
      real(r8), intent(in) :: tl_bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(out) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_z_w(LBi:,LBj:,0:)
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
#  ifndef SOLVE3D
      real(r8), intent(in) :: tl_bustr(LBi:,LBj:)
      real(r8), intent(in) :: tl_bvstr(LBi:,LBj:)
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
      real(r8), intent(inout) :: tl_h(LBi:,LBj:)
      real(r8), intent(inout) :: tl_rubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_rvbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_rzeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

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
      real(r8), intent(in):: tl_bed_thick(LBi:UBi,LBj:UBi,2)
#   endif
      real(r8), intent(out) :: tl_Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(out) :: tl_z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(out) :: tl_z_w(LBi:UBi,LBj:UBj,0:UBk)
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
#  ifndef SOLVE3D
      real(r8), intent(in) :: tl_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_bvstr(LBi:UBi,LBj:UBj)
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
      real(r8), intent(inout) :: tl_h(LBi:UBi,LBj:UBj)
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
      real(r8) :: tl_cff, tl_cff1, tl_fac, tl_fac1

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
          tl_Drhs(i,j)=tl_zeta(i,j,krhs)+tl_h(i,j)
        END DO
      END DO
      DO j=-2+J_RANGE+1
        DO i=-1+I_RANGE+1
          cff=0.5_r8*on_u(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i-1,j))
          tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i-1,j))
          DUon(i,j)=ubar(i,j,krhs)*cff1
          tl_DUon(i,j)=tl_ubar(i,j,krhs)*cff1+                          &
     &                 ubar(i,j,krhs)*tl_cff1
# ifdef NEARSHORE_MELLOR
          DUSon(i,j)=ubar_stokes(i,j)*cff1
          tl_DUSon(i,j)=tl_ubar_stokes(i,j)*cff1+                       &
     &                  ubar_stokes(i,j)*tl_cff1
          DUon(i,j)=DUon(i,j)+DUSon(i,j)
          tl_DUon(i,j)=tl_DUon(i,j)+tl_DUSon(i,j)
# endif
        END DO
      END DO
      DO j=-1+J_RANGE+1
        DO i=-2+I_RANGE+1
          cff=0.5_r8*om_v(i,j)
          cff1=cff*(Drhs(i,j)+Drhs(i,j-1))
          tl_cff1=cff*(tl_Drhs(i,j)+tl_Drhs(i,j-1))
          DVom(i,j)=vbar(i,j,krhs)*cff1
          tl_DVom(i,j)=tl_vbar(i,j,krhs)*cff1+                          &
     &                 vbar(i,j,krhs)*tl_cff1
# ifdef NEARSHORE_MELLOR
          DVSom(i,j)=vbar_stokes(i,j)*cff1
          tl_DVSom(i,j)=tl_vbar_stokes(i,j)*cff1+                       &
     &                  vbar_stokes(i,j)*tl_cff1
          DVom(i,j)=DVom(i,j)+DVSom(i,j)
          tl_DVom(i,j)=tl_DVom(i,j)+tl_DVSom(i,j)
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
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        tl_DUon)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        DVom)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        tl_DVom)
#  endif
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    DUon, DVom)
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_DUon, tl_DVom)
# endif
# undef I_RANGE
# undef J_RANGE
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
# ifdef OBC_VOLCONS
!
!  Compute integral mass flux across open boundaries and adjust
!  for volume conservation. Compute BASIC STATE value.
!  This needs to be computed here instead of below.
!
      CALL obc_flux_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    knew,                                         &
#  ifdef MASKING
     &                    umask, vmask,                                 &
#  endif
     &                    h, om_v, on_u,                                &
     &                    ubar, vbar, zeta)
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
      CALL tl_set_DUV_bc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         krhs,                                    &
#  ifdef MASKING
     &                         umask, vmask,                            &
#  endif
     &                         om_v, on_u, ubar, vbar,                  &
     &                         tl_ubar, tl_vbar,                        &
     &                         Drhs, DUon, DVom,                        &
     &                         tl_Drhs, tl_DUon, tl_DVom)
# endif
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
!  After all fast time steps are completed, recompute S-coordinate
!  surfaces according to the new free surface field.  Apply boundary
!  conditions to time averaged fields.
!
      IF ((iif(ng).eq.(nfast(ng)+1)).and.PREDICTOR_2D_STEP(ng)) THEN
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          Zt_avg1)
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_Zt_avg1)
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          DU_avg1)
!>
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_DU_avg1)
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          DV_avg1)
!>
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_DV_avg1)
#  endif
#  ifdef DISTRIBUTE
!>      CALL mp_exchange2d (ng, tile, iNLM, 3,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints, EWperiodic, NSperiodic,       &
!>   &                      Zt_avg1, DU_avg1, DV_avg1)
!>
        CALL mp_exchange2d (ng, tile, iTLM, 3,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      tl_Zt_avg1, tl_DU_avg1, tl_DV_avg1)
#  endif
!>      CALL set_depth_tile (ng, tile,                                  &
!>   &                       LBi, UBi, LBj, UBj,                        &
!>   &                       IminS, ImaxS, JminS, JmaxS,                &
!>   &                       nstp, nnew,                                &
!>   &                       h,                                         &
#  ifdef ICESHELF
!>   &                       zice,                                      &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>   &                       bed_thick,                                 &
#  endif
!>   &                       Zt_avg1, Hz, z_r, z_w)
!>
        CALL tl_set_depth_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp, nnew,                             &
     &                          h, tl_h,                                &
#  ifdef ICESHELF
     &                          zice,                                   &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                          tl_bed_thick,                           &
#  endif
     &                          Zt_avg1, tl_Zt_avg1,                    &
     &                          tl_Hz, tl_z_r, tl_z_w)
      END IF
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
      IF (FIRST_2D_STEP) THEN
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
# if defined EW_PERIODIC || defined NS_PERIODIC
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          rzeta(:,:,krhs))
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_rzeta(:,:,krhs))
# endif
# ifdef DISTRIBUTE
!>      CALL mp_exchange2d (ng, tile, iNLM, 1,                          &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      NghostPoints, EWperiodic, NSperiodic,       &
!>   &                      rzeta(:,:,krhs))
!>
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      tl_rzeta(:,:,krhs))
# endif
      END IF

# ifdef Q_PSOURCE
!
!  Apply mass point sources - Volume influx.
!
      DO is=1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        IF (((IstrR.le.i).and.(i.le.IendR)).and.                        &
     &      ((JstrR.le.j).and.(j.le.JendR))) THEN
!>        zeta(i,j,knew)=zeta(i,j,knew)+Qbar(is)*pm(i,j)*pn(i,j)*       &
!>   &                   dtfast(ng)
!>
!!        tl_zeta(i,j,knew)=tl_zeta(i,j,knew)+0.0_r8
        END IF
      END DO
# endif
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
# if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        zeta(:,:,knew))
!>
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_zeta(:,:,knew))
# endif
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    zeta(:,:,knew))
!>
      CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
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
!
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
#   ifdef EW_PERIODIC
#    define IU_RANGE IstrU-1,Iend+1
#   else
#    define IU_RANGE MAX(IstrU-1,2),MIN(Iend+1,Lm(ng))
#   endif
      DO j=Jstr,Jend
        DO i=IU_RANGE
          grad (i,j)=ubar(i-1,j,krhs)-2.0_r8*ubar(i,j,krhs)+            &
#    ifdef NEARSHORE_MELLOR
     &               ubar_stokes(i-1,j)-2.0_r8*ubar_stokes(i,j)+        &
     &               ubar_stokes(i+1,j)+                                &
#    endif
     &               ubar(i+1,j,krhs)
          tl_grad(i,j)=tl_ubar(i-1,j,krhs)-2.0_r8*tl_ubar(i,j,krhs)+    &
#    ifdef NEARSHORE_MELLOR
     &                 tl_ubar_stokes(i-1,j)-2.0_r8*tl_ubar_stokes(i,j)+&
     &                 tl_ubar_stokes(i+1,j)+                           &
#    endif
     &                 tl_ubar(i+1,j,krhs)
          Dgrad(i,j)=DUon(i-1,j)-2.0_r8*DUon(i,j)+DUon(i+1,j)
          tl_Dgrad(i,j)=tl_DUon(i-1,j)-2.0_r8*tl_DUon(i,j)+             &
     &                  tl_DUon(i+1,j)
        END DO
      END DO
#   undef IU_RANGE
#   ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN
        DO j=Jstr,Jend
          grad (Istr,j)=grad (Istr+1,j)
          tl_grad (Istr,j)=tl_grad (Istr+1,j)
          Dgrad(Istr,j)=Dgrad(Istr+1,j)
          tl_Dgrad(Istr,j)=tl_Dgrad(Istr+1,j)
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=Jstr,Jend
          grad (Iend+1,j)=grad (Iend,j)
          tl_grad (Iend+1,j)=tl_grad (Iend,j)
          Dgrad(Iend+1,j)=Dgrad(Iend,j)
          tl_Dgrad(Iend+1,j)=tl_Dgrad(Iend,j)
        END DO
      END IF
#   endif
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
#   ifdef NS_PERIODIC
#    define JU_RANGE Jstr-1,Jend+1
#   else
#    define JU_RANGE MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
#   endif
      DO j=JU_RANGE
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
#   undef JU_RANGE
#   ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
        DO i=IstrU,Iend
          grad(i,Jstr-1)=grad(i,Jstr)
          tl_grad(i,Jstr-1)=tl_grad(i,Jstr)
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=IstrU,Iend
          grad(i,Jend+1)=grad(i,Jend)
          tl_grad(i,Jend+1)=tl_grad(i,Jend)
        END DO
      END IF
#   endif
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
#   ifdef EW_PERIODIC
#    define IV_RANGE Istr-1,Iend+1
#   else
#    define IV_RANGE MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
#   endif
      DO j=JstrV,Jend
        DO i=IV_RANGE
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
#   undef IV_RANGE
#   ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN
        DO j=JstrV,Jend
          grad(Istr-1,j)=grad(Istr,j)
          tl_grad(Istr-1,j)=tl_grad(Istr,j)
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JstrV,Jend
          grad(Iend+1,j)=grad(Iend,j)
          tl_grad(Iend+1,j)=tl_grad(Iend,j)
        END DO
      END IF
#   endif
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
#   ifdef NS_PERIODIC
#    define JV_RANGE JstrV-1,Jend+1
#   else
#    define JV_RANGE MAX(JstrV-1,2),MIN(Jend+1,Mm(ng))
#   endif
      DO j=JV_RANGE
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
#   undef JV_RANGE
#   ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
        DO i=Istr,Iend
          grad (i,Jstr)=grad (i,Jstr+1)
          tl_grad (i,Jstr)=tl_grad (i,Jstr+1)
          Dgrad(i,Jstr)=Dgrad(i,Jstr+1)
          tl_Dgrad(i,Jstr)=tl_Dgrad(i,Jstr+1)
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=Istr,Iend
          grad (i,Jend+1)=grad (i,Jend)
          tl_grad (i,Jend+1)=tl_grad (i,Jend)
          Dgrad(i,Jend+1)=Dgrad(i,Jend)
          tl_Dgrad(i,Jend+1)=tl_Dgrad(i,Jend)
        END DO
      END IF
#   endif
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
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        fac=(UFx(i,j)-UFx(i-1,j))+                                    &
!>   &        (UFe(i,j+1)-UFe(i,j))
!>
          tl_fac=(tl_UFx(i,j)-tl_UFx(i-1,j))+                           &
     &           (tl_UFe(i,j+1)-tl_UFe(i,j))
!>        rhs_ubar(i,j)=rhs_ubar(i,j)-fac
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2hadv)=-fac
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        fac=(VFx(i+1,j)-VFx(i,j))+                                    &
!>   &        (VFe(i,j)-VFe(i,j-1))
!>
          tl_fac=(tl_VFx(i+1,j)-tl_VFx(i,j))+                           &
     &           (tl_VFe(i,j)-tl_VFe(i,j-1))
!>        rhs_vbar(i,j)=rhs_vbar(i,j)-fac
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2hadv)=-fac
#  endif
        END DO
      END DO
# endif
# if defined UV_COR || (defined CURVGRID && defined UV_ADV)
!
!-----------------------------------------------------------------------
!  Add in Coriolis and curvilinear transformation terms, if any.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff=0.5_r8*Drhs(i,j)*(                                        &
#  ifdef UV_COR
     &        fomn(i,j)                                                 &
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     &       +0.5_r8*((vbar(i,j  ,krhs)+                                &
#   ifdef NEARSHORE_MELLOR
     &                 vbar_stokes(i,j  )+                              &
     &                 vbar_stokes(i,j+1)+                              &
#   endif
     &                 vbar(i,j+1,krhs))*dndx(i,j)-                     &
     &                (ubar(i  ,j,krhs)+                                &
#   ifdef NEARSHORE_MELLOR
     &                 ubar_stokes(i  ,j)+                              &
     &                 ubar_stokes(i+1,j)+                              &
#   endif
     &                 ubar(i+1,j,krhs))*dmde(i,j))                     &
#  endif
     &        )
          tl_cff=0.5_r8*tl_Drhs(i,j)*(                                  &
#  ifdef UV_COR
     &           fomn(i,j)                                              &
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     &          +0.5_r8*((vbar(i,j  ,krhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                    vbar_stokes(i,j  )+                           &
     &                    vbar_stokes(i,j+1)+                           &
#   endif
     &                    vbar(i,j+1,krhs))*dndx(i,j)-                  &
     &                   (ubar(i  ,j,krhs)+                             &
#  ifdef NEARSHORE_MELLOR
     &                    ubar_stokes(i  ,j)+                           &
     &                    ubar_stokes(i+1,j)+                           &
#  endif
     &                    ubar(i+1,j,krhs))*dmde(i,j))                  &
#  endif
     &          )
#  if (defined CURVGRID && defined UV_ADV)
          tl_cff=tl_cff+                                                &
     &           0.25_r8*Drhs(i,j)*                                     &
     &           ((tl_vbar(i,j  ,krhs)+                                 &
#   ifdef NEARSHORE_MELLOR
     &             tl_vbar_stokes(i,j  )+                               &
     &             tl_vbar_stokes(i,j+1)+                               &
#   endif
     &             tl_vbar(i,j+1,krhs))*dndx(i,j)-                      &
     &            (tl_ubar(i  ,j,krhs)+                                 &
#   ifdef NEARSHORE_MELLOR
     &             tl_ubar_stokes(i  ,j)+                               &
     &             tl_ubar_stokes(i+1,j)+                               &
#   endif
     &             tl_ubar(i+1,j,krhs))*dmde(i,j))
#  endif
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
#  if defined DIAGNOSTICS_UV
#   ifdef UV_COR
!!        Uwrk(i,j)=0.5_r8*Drhs(i,j)*fomn(i,j)*                         &
!!   &              (vbar(i,j  ,krhs)+                                  &
#    ifdef NEARHSORE_MELLOR
!!   &               vbar_stokes(i,j  )+                                &
!!   &               vbar_stokes(i,j+1)+                                &
#    endif
!!                   vbar(i,j+1,krhs))
!!        Vwrk(i,j)=0.5_r8*Drhs(i,j)*fomn(i,j)*                         &
!!   &              (ubar(i  ,j,krhs)+                                  &
#    ifdef NEARSHORE_MELLOR
!!   &               ubar_stokes(i  ,j)+                                &
!!   &               ubar_stokes(i+1,j)+                                &
#    endif
!!   &               ubar(i+1,j,krhs))
#   else
!!        Uwrk(i,j)=0.0_r8
!!        Vwrk(i,j)=0.0_r8
#   endif
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
#   ifdef UV_COR
!!        DiaU2rhs(i,j,M2fcor)=fac2
#   endif
#   if (defined CURVGRID && defined UV_ADV)
!!        DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+fac1-fac2
#   endif
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
#   ifdef UV_COR
!!        DiaV2rhs(i,j,M2fcor)=-fac2
#   endif
#   if (defined CURVGRID && defined UV_ADV)
!!        DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)-fac1+fac2
#   endif
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
!>        fac=0.5_r8*((pn(i-1,j)+pn(i,j))*                              &
!>   &                (UFx(i,j  )-UFx(i-1,j))+                          &
!>   &                (pm(i-1,j)+pm(i,j))*                              &
!>   &                (UFe(i,j+1)-UFe(i  ,j)))
!>
          tl_fac=0.5_r8*((pn(i-1,j)+pn(i,j))*                           &
     &                   (tl_UFx(i,j  )-tl_UFx(i-1,j))+                 &
     &                   (pm(i-1,j)+pm(i,j))*                           &
     &                   (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
!>        rhs_ubar(i,j)=rhs_ubar(i,j)+fac
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2hvis)=fac
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        fac=0.5_r8*((pn(i,j-1)+pn(i,j))*                              &
!>   &                (VFx(i+1,j)-VFx(i,j  ))-                          &
!>   &                (pm(i,j-1)+pm(i,j))*                              &
!>   &                (VFe(i  ,j)-VFe(i,j-1)))
!>
          tl_fac=0.5_r8*((pn(i,j-1)+pn(i,j))*                           &
     &                   (tl_VFx(i+1,j)-tl_VFx(i,j  ))-                 &
     &                   (pm(i,j-1)+pm(i,j))*                           &
     &                   (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
!>        rhs_vbar(i,j)=rhs_vbar(i,j)+fac
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2hvis)=fac
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
      DO j=-1+JV_RANGE
        DO i=-1+IU_RANGE
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
      DO j=JU_RANGE+1
        DO i=IV_RANGE+1
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
#   ifdef MASKING
          cff=cff*pmask(i,j)
          tl_cff=tl_cff*pmask(i,j)
#   endif
          UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
          tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
          VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
          tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
        END DO
      END DO
!
!  Compute first harmonic operator (m s^-3/2).
!
      DO j=JU_RANGE
        DO i=IU_RANGE
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
      DO j=JV_RANGE
        DO i=IV_RANGE
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
#  ifndef EW_PERIODIC
      IF (WESTERN_EDGE) THEN
        DO j=JU_RANGE
#   ifdef WESTERN_WALL
          LapU(IstrU-1,j)=0.0_r8
          tl_LapU(IstrU-1,j)=0.0_r8
#   else
          LapU(IstrU-1,j)=LapU(IstrU,j)
          tl_LapU(IstrU-1,j)=tl_LapU(IstrU,j)
#   endif
        END DO
        DO j=JV_RANGE
#   ifdef WESTERN_WALL
          LapV(Istr-1,j)=gamma2(ng)*LapV(Istr,j)
          tl_LapV(Istr-1,j)=gamma2(ng)*tl_LapV(Istr,j)
#   else
          LapV(Istr-1,j)=0.0_r8
          tl_LapV(Istr-1,j)=0.0_r8
#   endif
        END DO
      END IF
      IF (EASTERN_EDGE) THEN
        DO j=JU_RANGE
#   ifdef EASTERN_WALL
          LapU(Iend+1,j)=0.0_r8
          tl_LapU(Iend+1,j)=0.0_r8
#   else
          LapU(Iend+1,j)=LapU(Iend,j)
          tl_LapU(Iend+1,j)=tl_LapU(Iend,j)
#   endif
        END DO
        DO j=JV_RANGE
#   ifdef EASTERN_WALL
          LapV(Iend+1,j)=gamma2(ng)*LapV(Iend,j)
          tl_LapV(Iend+1,j)=gamma2(ng)*tl_LapV(Iend,j)
#   else
          LapV(Iend+1,j)=0.0_r8
          tl_LapV(Iend+1,j)=0.0_r8
#   endif
        END DO
      END IF
#  endif
#  ifndef NS_PERIODIC
      IF (SOUTHERN_EDGE) THEN
        DO i=IU_RANGE
#   ifdef SOUTHERN_WALL
          LapU(i,Jstr-1)=gamma2(ng)*LapU(i,Jstr)
          tl_LapU(i,Jstr-1)=gamma2(ng)*tl_LapU(i,Jstr)
#   else
          LapU(i,Jstr-1)=0.0_r8
          tl_LapU(i,Jstr-1)=0.0_r8
#   endif
        END DO
        DO i=IV_RANGE
#   ifdef SOUTHERN_WALL
          LapV(i,JstrV-1)=0.0_r8
          tl_LapV(i,JstrV-1)=0.0_r8
#   else
          LapV(i,JstrV-1)=LapV(i,JstrV)
          tl_LapV(i,JstrV-1)=tl_LapV(i,JstrV)
#   endif
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO i=IU_RANGE
#   ifdef NORTHERN_WALL
          LapU(i,Jend+1)=gamma2(ng)*LapU(i,Jend)
          tl_LapU(i,Jend+1)=gamma2(ng)*tl_LapU(i,Jend)
#   else
          LapU(i,Jend+1)=0.0_r8
          tl_LapU(i,Jend+1)=0.0_r8
#   endif
        END DO
        DO i=IV_RANGE
#   ifdef NORTHERN_WALL
          LapV(i,Jend+1)=0.0_r8
          tl_LapV(i,Jend+1)=0.0_r8
#   else
          LapV(i,Jend+1)=LapV(i,Jend)
          tl_LapV(i,Jend+1)=tl_LapV(i,Jend)
#   endif
        END DO
      END IF
#  endif
#  if !defined EW_PERIODIC && !defined NS_PERIODIC
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        LapU(Istr  ,Jstr-1)=0.5_r8*(LapU(Istr+1,Jstr-1)+                &
     &                              LapU(Istr  ,Jstr  ))
        tl_LapU(Istr  ,Jstr-1)=0.5_r8*(tl_LapU(Istr+1,Jstr-1)+          &
     &                                 tl_LapU(Istr  ,Jstr  ))
        LapV(Istr-1,Jstr  )=0.5_r8*(LapV(Istr-1,Jstr+1)+                &
     &                              LapV(Istr  ,Jstr  ))
        tl_LapV(Istr-1,Jstr  )=0.5_r8*(tl_LapV(Istr-1,Jstr+1)+          &
     &                                 tl_LapV(Istr  ,Jstr  ))
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        LapU(Iend+1,Jstr-1)=0.5_r8*(LapU(Iend  ,Jstr-1)+                &
     &                              LapU(Iend+1,Jstr  ))
        tl_LapU(Iend+1,Jstr-1)=0.5_r8*(tl_LapU(Iend  ,Jstr-1)+          &
     &                                 tl_LapU(Iend+1,Jstr  ))
        LapV(Iend+1,Jstr  )=0.5_r8*(LapV(Iend  ,Jstr  )+                &
     &                              LapV(Iend+1,Jstr+1))
        tl_LapV(Iend+1,Jstr  )=0.5_r8*(tl_LapV(Iend  ,Jstr  )+          &
     &                                 tl_LapV(Iend+1,Jstr+1))
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        LapU(Istr  ,Jend+1)=0.5_r8*(LapU(Istr+1,Jend+1)+                &
     &                              LapU(Istr  ,Jend  ))
        tl_LapU(Istr  ,Jend+1)=0.5_r8*(tl_LapU(Istr+1,Jend+1)+          &
     &                                 tl_LapU(Istr  ,Jend  ))
        LapV(Istr-1,Jend+1)=0.5_r8*(LapV(Istr  ,Jend+1)+                &
     &                              LapV(Istr-1,Jend  ))
        tl_LapV(Istr-1,Jend+1)=0.5_r8*(tl_LapV(Istr  ,Jend+1)+          &
     &                                 tl_LapV(Istr-1,Jend  ))
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        LapU(Iend+1,Jend+1)=0.5_r8*(LapU(Iend  ,Jend+1)+                &
     &                              LapU(Iend+1,Jend  ))
        tl_LapU(Iend+1,Jend+1)=0.5_r8*(tl_LapU(Iend  ,Jend+1)+          &
     &                                 tl_LapU(Iend+1,Jend  ))
        LapV(Iend+1,Jend+1)=0.5_r8*(LapV(Iend  ,Jend+1)+                &
     &                              LapV(Iend+1,Jend  ))
        tl_LapV(Iend+1,Jend+1)=0.5_r8*(tl_LapV(Iend  ,Jend+1)+          &
     &                                 tl_LapV(Iend+1,Jend  ))
      END IF
#  endif
#  undef IU_RANGE
#  undef IV_RANGE
#  undef JU_RANGE
#  undef JV_RANGE
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
!>        fac=0.5_r8*((pn(i-1,j)+pn(i,j))*                              &
!>   &                (UFx(i,j  )-UFx(i-1,j))+                          &
!>   &                (pm(i-1,j)+pm(i,j))*                              &
!>   &                (UFe(i,j+1)-UFe(i  ,j)))
!>
          tl_fac=0.5_r8*((pn(i-1,j)+pn(i,j))*                           &
     &                   (tl_UFx(i,j  )-tl_UFx(i-1,j))+                 &
     &                   (pm(i-1,j)+pm(i,j))*                           &
     &                   (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
!>        rhs_ubar(i,j)=rhs_ubar(i,j)-fac
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaU2rhs(i,j,M2hvis)=-fac
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        fac=0.5_r8*((pn(i,j-1)+pn(i,j))*                              &
!>   &                (VFx(i+1,j)-VFx(i,j  ))-                          &
!>   &                (pm(i,j-1)+pm(i,j))*                              &
!>   &                (VFe(i  ,j)-VFe(i,j-1)))
!>
          tl_fac=0.5_r8*((pn(i,j-1)+pn(i,j))*                           &
     &                   (tl_VFx(i+1,j)-tl_VFx(i,j  ))-                 &
     &                   (pm(i,j-1)+pm(i,j))*                           &
     &                   (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
!>        rhs_vbar(i,j)=rhs_vbar(i,j)-fac
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
#  if defined DIAGNOSTICS_UV
!!        DiaV2rhs(i,j,M2hvis)=-fac
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
# ifdef M2CLM_NUDGING
!
!-----------------------------------------------------------------------
!  Add in nudging of 2D momentum climatology.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=0.25_r8*(M2nudgcof(i-1,j)+M2nudgcof(i,j))*                &
     &        om_u(i,j)*on_u(i,j)
!>        rhs_ubar(i,j)=rhs_ubar(i,j)+                                  &
!>   &                  cff*(Drhs(i-1,j)+Drhs(i,j))*                    &
!>   &                      (ubarclm(i,j)-ubar(i,j,krhs))
!>
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+                            &
     &                     cff*((Drhs(i-1,j)+Drhs(i,j))*                &
     &                          (-tl_ubar(i,j,krhs))+                   &
     &                          (tl_Drhs(i-1,j)+tl_Drhs(i,j))*          &
     &                          (ubarclm(i,j)-ubar(i,j,krhs)))
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff=0.25_r8*(M2nudgcof(i,j-1)+M2nudgcof(i,j))*                &
     &        om_v(i,j)*on_v(i,j)
!>        rhs_vbar(i,j)=rhs_vbar(i,j)+                                  &
!>   &                  cff*(Drhs(i,j-1)+Drhs(i,j))*                    &
!>   &                      (vbarclm(i,j)-vbar(i,j,krhs))
!>
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+                            &
     &                     cff*((Drhs(i,j-1)+Drhs(i,j))*                &
     &                          (-tl_vbar(i,j,krhs))+                   &
     &                          (tl_Drhs(i,j-1)+tl_Drhs(i,j))*          &
     &                          (vbarclm(i,j)-vbar(i,j,krhs)))
        END DO
      END DO
# endif
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
        IF (iic(ng).eq.ntfirst(ng)) THEN
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
!!            DiaRUfrc(i,j,3,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)-            &
!!   &                               DiaU2rhs(i,j,M2pgrd)
!!            DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+                &
!!   &                             DiaRUfrc(i,j,3,M2pgrd)
!!            DiaRUfrc(i,j,nstp,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)
!!            DiaU2rhs(i,j,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaU2rhs(i,j,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
!!            DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
#   ifdef UV_COR
!!            DiaRUfrc(i,j,3,M2fcor)=DiaRUfrc(i,j,3,M2fcor)-            &
!!   &                               DiaU2rhs(i,j,M2fcor)
!!            DiaU2rhs(i,j,M2fcor)=DiaU2rhs(i,j,M2fcor)+                &
!!   &                             DiaRUfrc(i,j,3,M2fcor)
!!            DiaRUfrc(i,j,nstp,M2fcor)=DiaRUfrc(i,j,3,M2fcor)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!            DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)-            &
!!   &                               DiaU2rhs(i,j,M2hvis)
!!            DiaU2rhs(i,j,M2hvis)=DiaU2rhs(i,j,M2hvis)+                &
!!   &                             DiaRUfrc(i,j,3,M2hvis)
!!            DiaRUfrc(i,j,nstp,M2hvis)=DiaRUfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_ADV
!!            DiaRUfrc(i,j,3,M2hadv)=DiaRUfrc(i,j,3,M2hadv)-            &
!!   &                               DiaU2rhs(i,j,M2hadv)
!!            DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+                &
!!   &                             DiaRUfrc(i,j,3,M2hadv)
!!            DiaRUfrc(i,j,nstp,M2hadv)=DiaRUfrc(i,j,3,M2hadv)
#   endif
#   ifdef NEARSHORE_MELLOR
!!            DiaRUfrc(i,j,3,M2hrad)=DiaRUfrc(i,j,3,M2hrad)-            &
!!   &                               DiaU2rhs(i,j,M2hrad)
!!            DiaU2rhs(i,j,M2hrad)=DiaU2rhs(i,j,M2hrad)+                &
!!   &                             DiaRUfrc(i,j,3,M2hrad)
!!            DiaRUfrc(i,j,nstp,M2hrad)=DiaRUfrc(i,j,3,M2hrad)
#   endif
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
!!            DiaRVfrc(i,j,3,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)-            &
!!   &                               DiaV2rhs(i,j,M2pgrd)
!!            DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+                &
!!   &                             DiaRVfrc(i,j,3,M2pgrd)
!!            DiaRVfrc(i,j,nstp,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)
!!            DiaV2rhs(i,j,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaV2rhs(i,j,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
!!            DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
#   ifdef UV_COR
!!            DiaRVfrc(i,j,3,M2fcor)=DiaRVfrc(i,j,3,M2fcor)-            &
!!   &                               DiaV2rhs(i,j,M2fcor)
!!            DiaV2rhs(i,j,M2fcor)=DiaV2rhs(i,j,M2fcor)+                &
!!   &                             DiaRVfrc(i,j,3,M2fcor)
!!            DiaRVfrc(i,j,nstp,M2fcor)=DiaRVfrc(i,j,3,M2fcor)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!            DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)-            &
!!   &                               DiaV2rhs(i,j,M2hvis)
!!            DiaV2rhs(i,j,M2hvis)=DiaV2rhs(i,j,M2hvis)+                &
!!   &                             DiaRVfrc(i,j,3,M2hvis)
!!            DiaRVfrc(i,j,nstp,M2hvis)=DiaRVfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_ADV
!!            DiaRVfrc(i,j,3,M2hadv)=DiaRVfrc(i,j,3,M2hadv)-            &
!!   &                               DiaV2rhs(i,j,M2hadv)
!!            DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)+                &
!!   &                             DiaRVfrc(i,j,3,M2hadv)
!!            DiaRVfrc(i,j,nstp,M2hadv)=DiaRVfrc(i,j,3,M2hadv)
#   endif
#   ifdef NEARSHORE_MELLOR
!!            DiaRVfrc(i,j,3,M2hrad)=DiaRVfrc(i,j,3,M2hrad)-            &
!!   &                               DiaV2rhs(i,j,M2hrad)
!!            DiaV2rhs(i,j,M2hrad)=DiaV2rhs(i,j,M2hrad)+                &
!!   &                             DiaRVfrc(i,j,3,M2hrad)
!!            DiaRVfrc(i,j,nstp,M2hrad)=DiaRVfrc(i,j,3,M2hrad)
#   endif
#  endif
            END DO
          END DO
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
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
!!            DiaRUfrc(i,j,3,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)-            &
!!   &                               DiaU2rhs(i,j,M2pgrd)
!!            DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+                &
!!   &                             1.5_r8*DiaRUfrc(i,j,3,M2pgrd)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2pgrd)
!!            DiaRUfrc(i,j,nstp,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)
!!            DiaU2rhs(i,j,M2bstr)=1.5_r8*DiaRUfrc(i,j,3,M2bstr)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2bstr)
!!            DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaU2rhs(i,j,M2sstr)=1.5_r8*DiaRUfrc(i,j,3,M2sstr)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2sstr)
!!            DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
#   ifdef UV_COR
!!            DiaRUfrc(i,j,3,M2fcor)=DiaRUfrc(i,j,3,M2fcor)-            &
!!   &                               DiaU2rhs(i,j,M2fcor)
!!            DiaU2rhs(i,j,M2fcor)=DiaU2rhs(i,j,M2fcor)+                &
!!   &                             1.5_r8*DiaRUfrc(i,j,3,M2fcor)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2fcor)
!!            DiaRUfrc(i,j,nstp,M2fcor)=DiaRUfrc(i,j,3,M2fcor)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!            DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)-            &
!!   &                               DiaU2rhs(i,j,M2hvis)
!!            DiaU2rhs(i,j,M2hvis)=DiaU2rhs(i,j,M2hvis)+                &
!!   &                             1.5_r8*DiaRUfrc(i,j,3,M2hvis)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2hvis)
!!            DiaRUfrc(i,j,nstp,M2hvis)=DiaRUfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_ADV
!!            DiaRUfrc(i,j,3,M2hadv)=DiaRUfrc(i,j,3,M2hadv)-            &
!!   &                               DiaU2rhs(i,j,M2hadv)
!!            DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+                &
!!   &                             1.5_r8*DiaRUfrc(i,j,3,M2hadv)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2hadv)
!!            DiaRUfrc(i,j,nstp,M2hadv)=DiaRUfrc(i,j,3,M2hadv)
#   endif
#   ifdef NEARSHORE_MELLOR
!!            DiaRUfrc(i,j,3,M2hrad)=DiaRUfrc(i,j,3,M2hrad)-            &
!!   &                               DiaU2rhs(i,j,M2hrad)
!!            DiaU2rhs(i,j,M2hrad)=DiaU2rhs(i,j,M2hrad)+                &
!!   &                             1.5_r8*DiaRUfrc(i,j,3,M2hrad)-       &
!!   &                             0.5_r8*DiaRUfrc(i,j,nnew,M2hrad)
!!            DiaRUfrc(i,j,nstp,M2hrad)=DiaRUfrc(i,j,3,M2hrad)
#   endif
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
!!            DiaRVfrc(i,j,3,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)-            &
!!   &                               DiaV2rhs(i,j,M2pgrd)
!!            DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+                &
!!   &                             1.5_r8*DiaRVfrc(i,j,3,M2pgrd)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2pgrd)
!!            DiaRVfrc(i,j,nstp,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)
!!            DiaV2rhs(i,j,M2bstr)=1.5_r8*DiaRVfrc(i,j,3,M2bstr)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2bstr)
!!            DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaV2rhs(i,j,M2sstr)=1.5_r8*DiaRVfrc(i,j,3,M2sstr)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2sstr)
!!            DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
#   ifdef UV_COR
!!            DiaRVfrc(i,j,3,M2fcor)=DiaRVfrc(i,j,3,M2fcor)-            &
!!   &                               DiaV2rhs(i,j,M2fcor)
!!            DiaV2rhs(i,j,M2fcor)=DiaV2rhs(i,j,M2fcor)+                &
!!   &                             1.5_r8*DiaRVfrc(i,j,3,M2fcor)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2fcor)
!!            DiaRVfrc(i,j,nstp,M2fcor)=DiaRVfrc(i,j,3,M2fcor)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!            DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)-            &
!!   &                               DiaV2rhs(i,j,M2hvis)
!!            DiaV2rhs(i,j,M2hvis)=DiaV2rhs(i,j,M2hvis)+                &
!!   &                             1.5_r8*DiaRVfrc(i,j,3,M2hvis)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2hvis)
!!            DiaRVfrc(i,j,nstp,M2hvis)=DiaRVfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_ADV
!!            DiaRVfrc(i,j,3,M2hadv)=DiaRVfrc(i,j,3,M2hadv)-            &
!!   &                               DiaV2rhs(i,j,M2hadv)
!!            DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)+                &
!!   &                             1.5_r8*DiaRVfrc(i,j,3,M2hadv)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2hadv)
!!            DiaRVfrc(i,j,nstp,M2hadv)=DiaRVfrc(i,j,3,M2hadv)
#   endif
#   ifdef NEARSHORE_MELLOR
!!            DiaRVfrc(i,j,3,M2hrad)=DiaRVfrc(i,j,3,M2hrad)-            &
!!   &                               DiaV2rhs(i,j,M2hrad)
!!            DiaV2rhs(i,j,M2hrad)=DiaV2rhs(i,j,M2hrad)+                &
!!   &                             1.5_r8*DiaRVfrc(i,j,3,M2hrad)-       &
!!   &                             0.5_r8*DiaRVfrc(i,j,nnew,M2hrad)
!!            DiaRVfrc(i,j,nstp,M2hrad)=DiaRVfrc(i,j,3,M2hrad)
#   endif
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
!!            DiaRUfrc(i,j,3,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)-            &
!!   &                               DiaU2rhs(i,j,M2pgrd)
!!            DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+                &
!!   &                             cff1*DiaRUfrc(i,j,3,M2pgrd)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2pgrd)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2pgrd)
!!            DiaRUfrc(i,j,nstp,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)
!!            DiaU2rhs(i,j,M2bstr)=cff1*DiaRUfrc(i,j,3,M2bstr)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2bstr)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2bstr)
!!            DiaRUfrc(i,j,nstp,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!            DiaU2rhs(i,j,M2sstr)=cff1*DiaRUfrc(i,j,3,M2sstr)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2sstr)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2sstr)
!!            DiaRUfrc(i,j,nstp,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
#   ifdef UV_COR
!!            DiaRUfrc(i,j,3,M2fcor)=DiaRUfrc(i,j,3,M2fcor)-            &
!!   &                               DiaU2rhs(i,j,M2fcor)
!!            DiaU2rhs(i,j,M2fcor)=DiaU2rhs(i,j,M2fcor)+                &
!!   &                             cff1*DiaRUfrc(i,j,3,M2fcor)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2fcor)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2fcor)
!!            DiaRUfrc(i,j,nstp,M2fcor)=DiaRUfrc(i,j,3,M2fcor)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!            DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)-            &
!!   &                               DiaU2rhs(i,j,M2hvis)
!!            DiaU2rhs(i,j,M2hvis)=DiaU2rhs(i,j,M2hvis)+                &
!!   &                             cff1*DiaRUfrc(i,j,3,M2hvis)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2hvis)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2hvis)
!!            DiaRUfrc(i,j,nstp,M2hvis)=DiaRUfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_ADV
!!            DiaRUfrc(i,j,3,M2hadv)=DiaRUfrc(i,j,3,M2hadv)-            &
!!   &                               DiaU2rhs(i,j,M2hadv)
!!            DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+                &
!!   &                             cff1*DiaRUfrc(i,j,3,M2hadv)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2hadv)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2hadv)
!!            DiaRUfrc(i,j,nstp,M2hadv)=DiaRUfrc(i,j,3,M2hadv)
#   endif
#   ifdef NEARSHORE_MELLOR
!!            DiaRUfrc(i,j,3,M2hrad)=DiaRUfrc(i,j,3,M2hrad)-            &
!!   &                               DiaU2rhs(i,j,M2hrad)
!!            DiaU2rhs(i,j,M2hrad)=DiaU2rhs(i,j,M2hrad)+                &
!!   &                             cff1*DiaRUfrc(i,j,3,M2hrad)-         &
!!   &                             cff2*DiaRUfrc(i,j,nnew,M2hrad)+      &
!!   &                             cff3*DiaRUfrc(i,j,nstp,M2hrad)
!!            DiaRUfrc(i,j,nstp,M2hrad)=DiaRUfrc(i,j,3,M2hrad)
#   endif
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
!!            DiaRVfrc(i,j,3,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)-            &
!!   &                               DiaV2rhs(i,j,M2pgrd)
!!            DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+                &
!!   &                             cff1*DiaRVfrc(i,j,3,M2pgrd)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2pgrd)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2pgrd)
!!            DiaRVfrc(i,j,nstp,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)
!!            DiaV2rhs(i,j,M2bstr)=cff1*DiaRVfrc(i,j,3,M2bstr)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2bstr)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2bstr)
!!            DiaRVfrc(i,j,nstp,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!            DiaV2rhs(i,j,M2sstr)=cff1*DiaRVfrc(i,j,3,M2sstr)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2sstr)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2sstr)
!!            DiaRVfrc(i,j,nstp,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
#   ifdef UV_COR
!!            DiaRVfrc(i,j,3,M2fcor)=DiaRVfrc(i,j,3,M2fcor)-            &
!!   &                               DiaV2rhs(i,j,M2fcor)
!!            DiaV2rhs(i,j,M2fcor)=DiaV2rhs(i,j,M2fcor)+                &
!!   &                             cff1*DiaRVfrc(i,j,3,M2fcor)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2fcor)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2fcor)
!!            DiaRVfrc(i,j,nstp,M2fcor)=DiaRVfrc(i,j,3,M2fcor)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!            DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)-            &
!!   &                               DiaV2rhs(i,j,M2hvis)
!!            DiaV2rhs(i,j,M2hvis)=DiaV2rhs(i,j,M2hvis)+                &
!!   &                             cff1*DiaRVfrc(i,j,3,M2hvis)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2hvis)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2hvis)
!!            DiaRVfrc(i,j,nstp,M2hvis)=DiaRVfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_ADV
!!            DiaRVfrc(i,j,3,M2hadv)=DiaRVfrc(i,j,3,M2hadv)-            &
!!   &                               DiaV2rhs(i,j,M2hadv)
!!            DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)+                &
!!   &                             cff1*DiaRVfrc(i,j,3,M2hadv)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2hadv)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2hadv)
!!            DiaRVfrc(i,j,nstp,M2hadv)=DiaRVfrc(i,j,3,M2hadv)
#   endif
#   ifdef NEARSHORE_MELLOR
!!            DiaRVfrc(i,j,3,M2hrad)=DiaRVfrc(i,j,3,M2hrad)-            &
!!   &                               DiaV2rhs(i,j,M2hrad)
!!            DiaV2rhs(i,j,M2hrad)=DiaV2rhs(i,j,M2hrad)+                &
!!   &                             cff1*DiaRVfrc(i,j,3,M2hrad)-         &
!!   &                             cff2*DiaRVfrc(i,j,nnew,M2hrad)+      &
!!   &                             cff3*DiaRVfrc(i,j,nstp,M2hrad)
!!            DiaRVfrc(i,j,nstp,M2hrad)=DiaRVfrc(i,j,3,M2hrad)
#   endif
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
!!          DiaU2rhs(i,j,M2pgrd)=DiaU2rhs(i,j,M2pgrd)+                  &
!!   &                           DiaRUfrc(i,j,3,M2pgrd)
!!          DiaU2rhs(i,j,M2bstr)=DiaRUfrc(i,j,3,M2bstr)
!!          DiaU2rhs(i,j,M2sstr)=DiaRUfrc(i,j,3,M2sstr)
#   ifdef UV_COR
!!          DiaU2rhs(i,j,M2fcor)=DiaU2rhs(i,j,M2fcor)+                  &
!!   &                           DiaRUfrc(i,j,3,M2fcor)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!          DiaU2rhs(i,j,M2hvis)=DiaU2rhs(i,j,M2hvis)+                  &
!!   &                           DiaRUfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_ADV
!!          DiaU2rhs(i,j,M2hadv)=DiaU2rhs(i,j,M2hadv)+                  &
!!   &                           DiaRUfrc(i,j,3,M2hadv)
#   endif
#   ifdef NEARSHORE_MELLOR
!!          DiaU2rhs(i,j,M2hrad)=DiaU2rhs(i,j,M2hrad)+                  &
!!   &                           DiaRUfrc(i,j,3,M2hrad)
#   endif
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
!>          rhs_vbar(i,j)=rhs_vbar(i,j)+rvfrc(i,j)
!>
            tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
#  ifdef DIAGNOSTICS_UV
!!          DiaV2rhs(i,j,M2pgrd)=DiaV2rhs(i,j,M2pgrd)+                  &
!!   &                           DiaRVfrc(i,j,3,M2pgrd)
!!          DiaV2rhs(i,j,M2bstr)=DiaRVfrc(i,j,3,M2bstr)
!!          DiaV2rhs(i,j,M2sstr)=DiaRVfrc(i,j,3,M2sstr)
#   ifdef UV_COR
!!          DiaV2rhs(i,j,M2fcor)=DiaV2rhs(i,j,M2fcor)+                  &
!!   &                           DiaRVfrc(i,j,3,M2fcor)
#   endif
#   if defined UV_VIS2 || defined UV_VIS4
!!          DiaV2rhs(i,j,M2hvis)=DiaV2rhs(i,j,M2hvis)+                  &
!!   &                           DiaRVfrc(i,j,3,M2hvis)
#   endif
#   ifdef UV_ADV
!!          DiaV2rhs(i,j,M2hadv)=DiaV2rhs(i,j,M2hadv)+                  &
!!   &                           DiaRVfrc(i,j,3,M2hadv)
#   endif
#   ifdef NEARSHORE_MELLOR
!!          DiaV2rhs(i,j,M2hrad)=DiaV2rhs(i,j,M2hrad)+                  &
!!   &                           DiaRVfrc(i,j,3,M2hrad)
#   endif
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
      IF (FIRST_2D_STEP) THEN
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
            IF (umask_wet(i,j).eq.1.0_r8) THEN
              IF (rmask_wet(i-1,j).eq.1.0_r8) THEN
!>              ubar(i,j,knew)=MAX(ubar(i,j,knew),0.0_r8)
!>
                tl_ubar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8, ubar(i,j,knew)))*       &
     &                            tl_ubar(i,j,knew)
              ELSE
!>              ubar(i,j,knew)=MIN(ubar(i,j,knew),0.0_r8)
!>
                tl_ubar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8,-ubar(i,j,knew)))*       &
     &                            tl_ubar(i,j,knew)
              END IF
            ELSE
!>            ubar(i,j,knew)=0.5_r8*ubar(i,j,knew)*umask_wet(i,j)
!>
              tl_ubar(i,j,knew)=0.5_r8*tl_ubar(i,j,knew)*umask_wet(i,j)
            END IF
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
            IF (vmask_wet(i,j).eq.1.0_r8) THEN
              IF (rmask_wet(i,j-1).eq.1.0_r8) THEN
!>              vbar(i,j,knew)=MAX(vbar(i,j,knew),0.0_r8)
!>
                tl_vbar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8, vbar(i,j,knew)))*       &
     &                            tl_vbar(i,j,knew)
              ELSE
!>              vbar(i,j,knew)=MIN(vbar(i,j,knew),0.0_r8)
!>
                tl_vbar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8,-vbar(i,j,knew)))*       &
     &                            tl_vbar(i,j,knew)
              END IF
            ELSE
!>            vbar(i,j,knew)=0.5_r8*vbar(i,j,knew)*vmask_wet(i,j)
!>
              tl_vbar(i,j,knew)=0.5_r8*tl_vbar(i,j,knew)*vmask_wet(i,j)
            END IF
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
            IF (umask_wet(i,j).eq.1.0_r8) THEN
              IF (rmask_wet(i-1,j).eq.1.0_r8) THEN
!>              ubar(i,j,knew)=MAX(ubar(i,j,knew),0.0_r8)
!>
                tl_ubar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8, ubar(i,j,knew)))*       &
     &                            tl_ubar(i,j,knew)
              ELSE
!>              ubar(i,j,knew)=MIN(ubar(i,j,knew),0.0_r8)
!>
                tl_ubar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8,-ubar(i,j,knew)))*       &
     &                            tl_ubar(i,j,knew)
              END IF
            ELSE
!>            ubar(i,j,knew)=0.5_r8*ubar(i,j,knew)*umask_wet(i,j)
!>
              tl_ubar(i,j,knew)=0.5_r8*tl_ubar(i,j,knew)*umask_wet(i,j)
            END IF
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
            IF (vmask_wet(i,j).eq.1.0_r8) THEN
              IF (rmask_wet(i,j-1).eq.1.0_r8) THEN
!>              vbar(i,j,knew)=MAX(vbar(i,j,knew),0.0_r8)
!>
                tl_vbar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8, vbar(i,j,knew)))*       &
     &                            tl_vbar(i,j,knew)
              ELSE
!>              vbar(i,j,knew)=MIN(vbar(i,j,knew),0.0_r8)
!>
                tl_vbar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8,-vbar(i,j,knew)))*       &
     &                            tl_vbar(i,j,knew)
              END IF
            ELSE
!>            vbar(i,j,knew)=0.5_r8*vbar(i,j,knew)*vmask_wet(i,j)
!>
              tl_vbar(i,j,knew)=0.5_r8*tl_vbar(i,j,knew)*vmask_wet(i,j)
            END IF
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
            IF (umask_wet(i,j).eq.1.0_r8) THEN
              IF (rmask_wet(i-1,j).eq.1.0_r8) THEN
!>              ubar(i,j,knew)=MAX(ubar(i,j,knew),0.0_r8)
!>
                tl_ubar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8, ubar(i,j,knew)))*       &
     &                            tl_ubar(i,j,knew)
              ELSE
!>              ubar(i,j,knew)=MIN(ubar(i,j,knew),0.0_r8)
!>
                tl_ubar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8,-ubar(i,j,knew)))*       &
     &                            tl_ubar(i,j,knew)
              END IF
            ELSE
!>            ubar(i,j,knew)=0.5_r8*ubar(i,j,knew)*umask_wet(i,j)
!>
              tl_ubar(i,j,knew)=0.5_r8*tl_ubar(i,j,knew)*umask_wet(i,j)
            END IF
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
            IF (vmask_wet(i,j).eq.1.0_r8) THEN
              IF (rmask_wet(i,j-1).eq.1.0_r8) THEN
!>              vbar(i,j,knew)=MAX(vbar(i,j,knew),0.0_r8)
!>
                tl_vbar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8, vbar(i,j,knew)))*       &
     &                            tl_vbar(i,j,knew)
              ELSE
!>              vbar(i,j,knew)=MIN(vbar(i,j,knew),0.0_r8)
!>
                tl_vbar(i,j,knew)=(0.5_r8+                              &
     &                             SIGN(0.5_r8,-vbar(i,j,knew)))*       &
     &                            tl_vbar(i,j,knew)
              END IF
            ELSE
!>            vbar(i,j,knew)=0.5_r8*vbar(i,j,knew)*vmask_wet(i,j)
!>
              tl_vbar(i,j,knew)=0.5_r8*tl_vbar(i,j,knew)*vmask_wet(i,j)
            END IF
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

# if defined WET_DRY_NOT_YET
!
!-----------------------------------------------------------------------
! If wet/drying, compute new masks for cells with depth < Dcrit.
!-----------------------------------------------------------------------
!
! HGA:  We need to think more about TLM of the wet/dry mask arrays
!       since they are time-dependent.

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
# ifdef OBC_VOLCONS
!
!  Compute integral mass flux across open boundaries and adjust
!  for volume conservation.
!
      CALL tl_obc_flux_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       knew,                                      &
#  ifdef MASKING
     &                       umask, vmask,                              &
#  endif
     &                       h, tl_h, om_v, on_u,                       &
     &                       ubar, vbar, zeta,                          &
     &                       tl_ubar, tl_vbar, tl_zeta)
# endif
# ifdef UV_PSOURCE
!
!-----------------------------------------------------------------------
!  Apply mass point sources.
!-----------------------------------------------------------------------
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Dnew(i,j)=zeta(i,j,knew)+h(i,j)
          tl_Dnew(i,j)=tl_zeta(i,j,knew)+tl_h(i,j)
        END DO
      END DO
      DO is=1,Nsrc
        i=Isrc(is)
        j=Jsrc(is)
        IF (((IstrR.le.i).and.(i.le.IendR)).and.                        &
     &      ((JstrR.le.j).and.(j.le.JendR))) THEN
          IF (INT(Dsrc(is)).eq.0) THEN
            cff=1.0_r8/(on_u(i,j)*0.5_r8*(Dnew(i-1,j)+Dnew(i,j)))
            tl_cff=-cff*cff*on_u(i,j)*                                  &
     &             0.5_r8*(tl_Dnew(i-1,j)+tl_Dnew(i,j))
!>          ubar(i,j,knew)=Qbar(is)*cff
!>
            tl_ubar(i,j,knew)=Qbar(is)*tl_cff
#  ifdef SOLVE3D
!>          DU_avg1(i,j)=Qbar(is)
!>
            tl_DU_avg1(i,j)=0.0_r8
#  endif
          ELSE
            cff=1.0_r8/(om_v(i,j)*0.5_r8*(Dnew(i,j-1)+Dnew(i,j)))
            tl_cff=-cff*cff*om_v(i,j)*                                  &
     &             0.5_r8*(tl_Dnew(i,j-1)+tl_Dnew(i,j))
!>          vbar(i,j,knew)=Qbar(is)*cff
!>
            tl_vbar(i,j,knew)=Qbar(is)*tl_cff
#  ifdef SOLVE3D
!>          DV_avg1(i,j)=Qbar(is)
!>
            tl_DV_avg1(i,j)=0.0_r8
#  endif
          END IF
        END IF
      END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_u2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        ubar(:,:,knew))
!>
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_ubar(:,:,knew))
!>    CALL exchange_v2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        vbar(:,:,knew))
!>
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_vbar(:,:,knew))
#  endif
#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    ubar(:,:,knew),                               &
!>   &                    vbar(:,:,knew))
!>
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_ubar(:,:,knew),                            &
     &                    tl_vbar(:,:,knew))
#  endif
# endif
      RETURN
      END SUBROUTINE tl_step2d_tile
#else
      SUBROUTINE tl_step2d
      END SUBROUTINE tl_step2d
#endif
