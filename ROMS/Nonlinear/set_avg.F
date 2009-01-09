#include "cppdefs.h"
      MODULE set_avg_mod
#if defined AVERAGES && (!defined ADJOINT && defined NONLINEAR)
!
!svn $Id: set_avg.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_avg

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_average
# if defined FORWARD_WRITE && defined SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
# ifdef SOLVE3D
      USE mod_grid
      USE mod_mixing
# endif
      USE mod_ocean
      USE mod_stepping
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_tides
# endif
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
      CALL wclock_on (ng, iNLM, 5)
# endif
      CALL set_avg_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
     &                   NTC(ng),                                       &
# endif
     &                   KOUT,                                          &
# ifdef SOLVE3D
     &                   NOUT,                                          &
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
#  ifdef AVERAGES_QUADRATIC
     &                   GRID(ng) % Huon,                               &
     &                   GRID(ng) % Hvom,                               &
#  endif
#  ifdef FORWARD_WRITE
     &                   COUPLING(ng) % DU_avg1,                        &
     &                   COUPLING(ng) % DU_avg2,                        &
     &                   COUPLING(ng) % DV_avg1,                        &
     &                   COUPLING(ng) % DV_avg2,                        &
#  endif
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   OCEAN(ng) % W,                                 &
     &                   OCEAN(ng) % wvel,                              &
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % rho,                               &
#  ifdef AVERAGES_NEARSHORE
     &                   OCEAN(ng) % u_stokes,                          &
     &                   OCEAN(ng) % v_stokes,                          &
     &                   MIXING(ng) % rustr3d,                          &
     &                   MIXING(ng) % rvstr3d,                          &
     &                   MIXING(ng) % Sxx,                              &
     &                   MIXING(ng) % Sxy,                              &
     &                   MIXING(ng) % Syy,                              &
     &                   MIXING(ng) % Szx,                              &
     &                   MIXING(ng) % Szy,                              &
#  endif
#  ifdef LMD_SKPP
     &                   MIXING(ng) % hsbl,                             &
#  endif
#  ifdef LMD_BKPP
     &                   MIXING(ng) % hbbl,                             &
#  endif
#  ifdef AVERAGES_AKV
     &                   MIXING(ng) % Akv,                              &
#  endif
#  if defined AVERAGES_AKT || defined AVERAGES_AKS
     &                   MIXING(ng) % Akt,                              &
#  endif
#  ifdef AVERAGES_FLUXES
     &                   FORCES(ng) % stflx,                            &
#   ifdef BULK_FLUXES
     &                   FORCES(ng) % lhflx,                            &
     &                   FORCES(ng) % shflx,                            &
     &                   FORCES(ng) % lrflx,                            &
#    ifdef EMINUSP
     &                   FORCES(ng) % evap,                             &
     &                   FORCES(ng) % rain,                             &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                   FORCES(ng) % srflx,                            &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
# endif
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
     &                   OCEAN(ng) % zeta,                              &
# ifdef AVERAGES_NEARSHORE
     &                   OCEAN(ng) % ubar_stokes,                       &
     &                   OCEAN(ng) % vbar_stokes,                       &
     &                   MIXING(ng) % rustr2d,                          &
     &                   MIXING(ng) % rvstr2d,                          &
     &                   MIXING(ng) % Sxx_bar,                          &
     &                   MIXING(ng) % Sxy_bar,                          &
     &                   MIXING(ng) % Syy_bar,                          &
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
     &                   TIDES(ng) % CosOmega, TIDES(ng) % SinOmega,    &
     &                   TIDES(ng) % CosW_avg, TIDES(ng) % CosW_sum,    &
     &                   TIDES(ng) % SinW_avg, TIDES(ng) % SinW_sum,    &
     &                   TIDES(ng) % CosWCosW, TIDES(ng) % SinWSinW,    &
     &                   TIDES(ng) % SinWCosW,                          &
     &                   TIDES(ng) % ubar_tide,                         &
     &                   TIDES(ng) % ubar_detided,                      &
     &                   TIDES(ng) % vbar_tide,                         &
     &                   TIDES(ng) % vbar_detided,                      &
     &                   TIDES(ng) % zeta_tide,                         &
     &                   TIDES(ng) % zeta_detided,                      &
# endif
# ifdef SOLVE3D
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
     &                   TIDES(ng) % u_tide,                            &
     &                   TIDES(ng) % u_detided,                         &
     &                   TIDES(ng) % v_tide,                            &
     &                   TIDES(ng) % v_detided,                         &
#  endif
#  ifdef FORWARD_WRITE
     &                   AVERAGE(ng) % avgDU_avg1,                      &
     &                   AVERAGE(ng) % avgDU_avg2,                      &
     &                   AVERAGE(ng) % avgDV_avg1,                      &
     &                   AVERAGE(ng) % avgDV_avg2,                      &
#  endif
     &                   AVERAGE(ng) % avgu3d,                          &
     &                   AVERAGE(ng) % avgv3d,                          &
     &                   AVERAGE(ng) % avgw3d,                          &
     &                   AVERAGE(ng) % avgwvel,                         &
     &                   AVERAGE(ng) % avgt,                            &
     &                   AVERAGE(ng) % avgrho,                          &
#  ifdef AVERAGES_NEARSHORE
     &                   AVERAGE(ng) % avgu3Sd,                         &
     &                   AVERAGE(ng) % avgv3Sd,                         &
     &                   AVERAGE(ng) % avgu3RS,                         &
     &                   AVERAGE(ng) % avgv3RS,                         &
     &                   AVERAGE(ng) % avgSxx3d,                        &
     &                   AVERAGE(ng) % avgSxy3d,                        &
     &                   AVERAGE(ng) % avgSyy3d,                        &
     &                   AVERAGE(ng) % avgSzx3d,                        &
     &                   AVERAGE(ng) % avgSzy3d,                        &
#  endif
#  ifdef LMD_SKPP
     &                   AVERAGE(ng) % avghsbl,                         &
#  endif
#  ifdef LMD_BKPP
     &                   AVERAGE(ng) % avghbbl,                         &
#  endif
#  ifdef AVERAGES_AKV
     &                   AVERAGE(ng) % avgAKv,                          &
#  endif
#  ifdef AVERAGES_AKT
     &                   AVERAGE(ng) % avgAKt,                          &
#  endif
#  ifdef AVERAGES_AKS
     &                   AVERAGE(ng) % avgAKs,                          &
#  endif
#  ifdef AVERAGES_FLUXES
     &                   AVERAGE(ng) % avgstf,                          &
     &                   AVERAGE(ng) % avgswf,                          &
#   ifdef BULK_FLUXES
     &                   AVERAGE(ng) % avglhf,                          &
     &                   AVERAGE(ng) % avgshf,                          &
     &                   AVERAGE(ng) % avglrf,                          &
#    ifdef EMINUSP
     &                   AVERAGE(ng) % avgevap,                         &
     &                   AVERAGE(ng) % avgrain,                         &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                   AVERAGE(ng) % avgsrf,                          &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                   AVERAGE(ng) % avgsus,                          &
     &                   AVERAGE(ng) % avgsvs,                          &
     &                   AVERAGE(ng) % avgbus,                          &
     &                   AVERAGE(ng) % avgbvs,                          &
# endif
# ifdef AVERAGES_QUADRATIC
#  ifdef SOLVE3D
     &                   AVERAGE(ng) % avgHuon,                         &
     &                   AVERAGE(ng) % avgHvom,                         &
     &                   AVERAGE(ng) % avgUU,                           &
     &                   AVERAGE(ng) % avgUV,                           &
     &                   AVERAGE(ng) % avgVV,                           &
     &                   AVERAGE(ng) % avgHuonT,                        &
     &                   AVERAGE(ng) % avgHvomT,                        &
     &                   AVERAGE(ng) % avgUT,                           &
     &                   AVERAGE(ng) % avgVT,                           &
     &                   AVERAGE(ng) % avgTT,                           &
#  endif     
     &                   AVERAGE(ng) % avgU2,                           &
     &                   AVERAGE(ng) % avgV2,                           &
     &                   AVERAGE(ng) % avgZZ,                           &
# endif
# ifdef AVERAGES_NEARSHORE
     &                   AVERAGE(ng) % avgu2Sd,                         &
     &                   AVERAGE(ng) % avgv2Sd,                         &
     &                   AVERAGE(ng) % avgu2RS,                         &
     &                   AVERAGE(ng) % avgv2RS,                         &
     &                   AVERAGE(ng) % avgSxx2d,                        &
     &                   AVERAGE(ng) % avgSxy2d,                        &
     &                   AVERAGE(ng) % avgSyy2d,                        &
# endif
     &                   AVERAGE(ng) % avgu2d,                          &
     &                   AVERAGE(ng) % avgv2d,                          &
     &                   AVERAGE(ng) % avgzeta)

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_avg
!
!***********************************************************************
      SUBROUTINE set_avg_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
     &                         NTC,                                     &
# endif
     &                         Kout,                                    &
# ifdef SOLVE3D
     &                         Nout,                                    &
     &                         pm, pn,                                  &
#  ifdef AVERAGES_QUADRATIC
     &                         Huon, Hvom,                              &
#  endif
#  ifdef FORWARD_WRITE
     &                         DU_avg1, DU_avg2,                        &
     &                         DV_avg1, DV_avg2,                        &
#  endif
     &                         u, v,                                    &
     &                         W, wvel, t, rho,                         &
#  ifdef AVERAGES_NEARSHORE
     &                         u_stokes, v_stokes,                      &
     &                         rustr3d, rvstr3d,                        &
     &                         Sxx, Sxy, Syy, Szx, Szy,                 &
#  endif
#  ifdef LMD_SKPP
     &                         hsbl,                                    &
#  endif
#  ifdef LMD_BKPP
     &                         hbbl,                                    &
#  endif
#  ifdef AVERAGES_AKV
     &                         Akv,                                     &
#  endif
#  if defined AVERAGES_AKT || defined AVERAGES_AKS
     &                         Akt,                                     &
#  endif
#  ifdef AVERAGES_FLUXES
     &                         stflx,                                   &
#   ifdef BULK_FLUXES
     &                         lhflx, shflx, lrflx,                     &
#    ifdef EMINUSP
     &                         evap, rain,                              &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                         srflx,                                   &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                         sustr, svstr, bustr, bvstr,              &
# endif
     &                         ubar, vbar,                              &
     &                         zeta,                                    &
# ifdef AVERAGES_NEARSHORE
     &                         ubar_stokes, vbar_stokes,                &
     &                         rustr2d, rvstr2d,                        &
     &                         Sxx_bar, Sxy_bar, Syy_bar,               &
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
     &                         CosOmega, SinOmega,                      &
     &                         CosW_avg, CosW_sum,                      &
     &                         SinW_avg, SinW_sum,                      &
     &                         CosWCosW, SinWSinW, SinWCosW,            &
     &                         ubar_tide, ubar_detided,                 &
     &                         vbar_tide, vbar_detided,                 &
     &                         zeta_tide, zeta_detided,                 &
# endif
# ifdef SOLVE3D
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
     &                         u_tide, u_detided,                       &
     &                         v_tide, v_detided,                       &
#  endif
#  ifdef FORWARD_WRITE
     &                         avgDU_avg1, avgDU_avg2,                  &
     &                         avgDV_avg1, avgDV_avg2,                  &
#  endif
     &                         avgu3d, avgv3d,                          &
     &                         avgw3d, avgwvel,                         &
     &                         avgt, avgrho,                            &
#  ifdef AVERAGES_NEARSHORE
     &                         avgu3Sd, avgv3Sd,                        &
     &                         avgu3RS, avgv3RS,                        &
     &                         avgSxx3d, avgSxy3d, avgSyy3d,            &
     &                         avgSzx3d, avgSzy3d,                      &
#  endif
#  ifdef LMD_SKPP
     &                         avghsbl,                                 &
#  endif
#  ifdef LMD_BKPP
     &                         avghbbl,                                 &
#  endif
#  ifdef AVERAGES_AKV
     &                         avgAKv,                                  &
#  endif
#  ifdef AVERAGES_AKT
     &                         avgAKt,                                  &
#  endif
#  ifdef AVERAGES_AKS
     &                         avgAKs,                                  &
#  endif
#  ifdef AVERAGES_FLUXES
     &                         avgstf, avgswf,                          &
#   ifdef BULK_FLUXES
     &                         avglhf, avgshf, avglrf,                  &
#    ifdef EMINUSP
     &                         avgevap, avgrain,                        &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                         avgsrf,                                  &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                         avgsus, avgsvs, avgbus, avgbvs,          &
# endif
# ifdef AVERAGES_QUADRATIC
#  ifdef SOLVE3D
     &                         avgHuon, avgHvom,                        &
     &                         avgUU, avgUV, avgVV,                     &
     &                         avgHuonT, avgHvomT,                      &
     &                         avgUT, avgVT, avgTT,                     &
#  endif     
     &                         avgU2, avgV2, avgZZ,                     &
# endif
# ifdef AVERAGES_NEARSHORE
     &                         avgu2Sd, avgv2Sd,                        &
     &                         avgu2RS, avgv2RS,                        &
     &                         avgSxx2d, avgSxy2d, avgSyy2d,            &
# endif
     &                         avgu2d, avgv2d,                          &
     &                         avgzeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Kout
# ifdef SOLVE3D
      integer, intent(in) :: Nout
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      integer, intent(in) :: NTC
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#   ifdef AVERAGES_QUADRATIC
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
#   endif
#   ifdef FORWARD_WRITE
      real(r8), intent(in) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DU_avg2(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg2(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(in) :: wvel(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#   ifdef AVERAGES_NEARSHORE
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: rustr3d(LBi:,LBj:,:)
      real(r8), intent(in) :: rvstr3d(LBi:,LBj:,:)
      real(r8), intent(in) :: Sxx(LBi:,LBj:,:)
      real(r8), intent(in) :: Sxy(LBi:,LBj:,:)
      real(r8), intent(in) :: Syy(LBi:,LBj:,:)
      real(r8), intent(in) :: Szx(LBi:,LBj:,:)
      real(r8), intent(in) :: Szy(LBi:,LBj:,:)
#   endif
#   ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:,LBj:)
#   endif
#   ifdef LMD_BKPP
      real(r8), intent(in) :: hbbl(LBi:,LBj:)
#   endif
#   ifdef AVERAGES_AKV
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
#   endif
#   if defined AVERAGES_AKT || defined AVERAGES_AKS
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
#    ifdef BULK_FLUXES
      real(r8), intent(in) :: lhflx(LBi:,LBj:)
      real(r8), intent(in) :: shflx(LBi:,LBj:)
      real(r8), intent(in) :: lrflx(LBi:,LBj:)
#     ifdef EMINUSP
      real(r8), intent(in) :: evap(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(in) :: srflx(LBi:,LBj:)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef AVERAGES_NEARSHORE
      real(r8), intent(in) :: ubar_stokes(LBi:,LBj:)
      real(r8), intent(in) :: vbar_stokes(LBi:,LBj:)
      real(r8), intent(in) :: rustr2d(LBi:,LBj:)
      real(r8), intent(in) :: rvstr2d(LBi:,LBj:)
      real(r8), intent(in) :: Sxx_bar(LBi:,LBj:)
      real(r8), intent(in) :: Sxy_bar(LBi:,LBj:)
      real(r8), intent(in) :: Syy_bar(LBi:,LBj:)
#  endif
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      real(r8), intent(in) :: CosOmega(:)
      real(r8), intent(in) :: SinOmega(:)
      real(r8), intent(inout) :: CosW_avg(:)
      real(r8), intent(inout) :: CosW_sum(:)
      real(r8), intent(inout) :: SinW_avg(:)
      real(r8), intent(inout) :: SinW_sum(:)
      real(r8), intent(inout) :: CosWCosW(:,:)
      real(r8), intent(inout) :: SinWSinW(:,:)
      real(r8), intent(inout) :: SinWCosW(:,:)
      real(r8), intent(inout) :: ubar_tide(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ubar_detided(LBi:,LBj:)
      real(r8), intent(inout) :: vbar_tide(LBi:,LBj:,0:)
      real(r8), intent(inout) :: vbar_detided(LBi:,LBj:)
      real(r8), intent(inout) :: zeta_tide(LBi:,LBj:,0:)
      real(r8), intent(inout) :: zeta_detided(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      real(r8), intent(inout) :: u_tide(LBi:,LBj:,:,0:)
      real(r8), intent(inout) :: u_detided(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_tide(LBi:,LBj:,:,0:)
      real(r8), intent(inout) :: v_detided(LBi:,LBj:,:)
#   endif
#   ifdef FORWARD_WRITE
      real(r8), intent(inout) :: avgDU_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: avgDU_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: avgDV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: avgDV_avg2(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: avgu3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgv3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgw3d(LBi:,LBj:,0:)
      real(r8), intent(inout) :: avgwvel(LBi:,LBj:,0:)
      real(r8), intent(inout) :: avgt(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgrho(LBi:,LBj:,:)
#   ifdef AVERAGES_NEARSHORE
      real(r8), intent(inout) :: avgu3Sd(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgv3Sd(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgu3RS(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgv3RS(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgSxx3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgSxy3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgSyy3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgSzx3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgSzy3d(LBi:,LBj:,:)
#   endif
#   ifdef LMD_SKPP
      real(r8), intent(inout) :: avghsbl(LBi:,LBj:)
#   endif
#   ifdef LMD_BKPP
      real(r8), intent(inout) :: avghbbl(LBi:,LBj:)
#   endif
#   ifdef AVERAGES_AKV
      real(r8), intent(inout) :: avgAKv(LBi:,LBj:,0:)
#   endif
#   ifdef AVERAGES_AKT
      real(r8), intent(inout) :: avgAKt(LBi:,LBj:,0:)
#   endif
#   ifdef AVERAGES_AKS
      real(r8), intent(inout) :: avgAKs(LBi:,LBj:,0:)
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(inout) :: avgstf(LBi:,LBj:)
      real(r8), intent(inout) :: avgswf(LBi:,LBj:)
#    ifdef BULK_FLUXES
      real(r8), intent(inout) :: avglhf(LBi:,LBj:)
      real(r8), intent(inout) :: avgshf(LBi:,LBj:)
      real(r8), intent(inout) :: avglrf(LBi:,LBj:)
#     ifdef EMINUSP
      real(r8), intent(inout) :: avgevap(LBi:,LBj:)
      real(r8), intent(inout) :: avgrain(LBi:,LBj:)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: avgsrf(LBi:,LBj:)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(inout) :: avgsus(LBi:,LBj:)
      real(r8), intent(inout) :: avgsvs(LBi:,LBj:)
      real(r8), intent(inout) :: avgbus(LBi:,LBj:)
      real(r8), intent(inout) :: avgbvs(LBi:,LBj:)
#  endif
#  ifdef AVERAGES_QUADRATIC
#   ifdef SOLVE3D
      real(r8), intent(inout) :: avgHuon(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgHvom(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgUU(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgUV(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgVV(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgHuonT(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgHvomT(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgUT(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgVT(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgTT(LBi:,LBj:,:,:)
#   endif     
      real(r8), intent(inout) :: avgU2(LBi:,LBj:)
      real(r8), intent(inout) :: avgV2(LBi:,LBj:)
      real(r8), intent(inout) :: avgZZ(LBi:,LBj:)
#  endif
#  ifdef AVERAGES_NEARSHORE
      real(r8), intent(inout) :: avgu2Sd(LBi:,LBj:)
      real(r8), intent(inout) :: avgv2Sd(LBi:,LBj:)
      real(r8), intent(inout) :: avgu2RS(LBi:,LBj:)
      real(r8), intent(inout) :: avgv2RS(LBi:,LBj:)
      real(r8), intent(inout) :: avgSxx2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgSxy2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgSyy2d(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: avgu2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgv2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgzeta(LBi:,LBj:)

# else

#  ifdef SOLVE3D
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#   ifdef AVERAGES_QUADRATIC
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
#   endif
#   ifdef FORWARD_WRITE
      real(r8), intent(in) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg2(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: wvel(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
#   ifdef AVERAGES_NEARSHORE
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rustr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rvstr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Sxx(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Sxy(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Syy(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Szx(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Szy(LBi:UBi,LBj:UBj,N(ng))
#   endif
#   ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef LMD_BKPP
      real(r8), intent(in) :: hbbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef AVERAGES_AKV
      real(r8), intent(in) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined AVERAGES_AKT || defined AVERAGES_AKS
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#    ifdef BULK_FLUXES
      real(r8), intent(in) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: lrflx(LBi:UBi,LBj:UBj)
#     ifdef EMINUSP
      real(r8), intent(in) :: evap(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef AVERAGES_NEARSHORE
      real(r8), intent(in) :: ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rustr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rvstr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Sxx_bar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Sxy_bar(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Syy_bar(LBi:UBi,LBj:UBj)
#  endif
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      real(r8), intent(in) :: CosOmega(NTC)
      real(r8), intent(in) :: SinOmega(NTC)
      real(r8), intent(inout) :: CosW_avg(NTC)
      real(r8), intent(inout) :: CosW_sum(NTC)
      real(r8), intent(inout) :: SinW_avg(NTC)
      real(r8), intent(inout) :: SinW_sum(NTC)
      real(r8), intent(inout) :: CosWCosW(NTC,NTC)
      real(r8), intent(inout) :: SinWSinW(NTC,NTC)
      real(r8), intent(inout) :: SinWCosW(NTC,NTC)
      real(r8), intent(inout) :: ubar_tide(LBi:UBi,LBj:UBj,0:2*NTC)
      real(r8), intent(inout) :: ubar_detided(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vbar_tide(LBi:UBi,LBj:UBj,0:2*NTC)
      real(r8), intent(inout) :: vbar_detided(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: zeta_tide(LBi:UBi,LBj:UBj,0:2*NTC)
      real(r8), intent(inout) :: zeta_detided(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      real(r8), intent(inout) :: u_tide(LBi:UBi,LBj:UBj,N(ng),0:2*NTC)
      real(r8), intent(inout) :: u_detided(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v_tide(LBi:UBi,LBj:UBj,N(ng),0:2*NTC)
      real(r8), intent(inout) :: v_detided(LBi:UBi,LBj:UBj,N(ng))
#   endif
#   ifdef FORWARD_WRITE
      real(r8), intent(inout) :: avgDU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgDU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgDV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgDV_avg2(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: avgu3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgv3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgw3d(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: avgwvel(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: avgt(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: avgrho(LBi:UBi,LBj:UBj,N(ng))
#   ifdef AVERAGES_NEARSHORE
      real(r8), intent(inout) :: avgu3Sd(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgv3Sd(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgu3RS(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgv3RS(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgSxx3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgSxy3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgSyy3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgSzx3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgSzy3d(LBi:UBi,LBj:UBj,N(ng))
#   endif
#   ifdef LMD_SKPP
      real(r8), intent(inout) :: avghsbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef LMD_BKPP
      real(r8), intent(inout) :: avghbbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef AVERAGES_AKV
      real(r8), intent(inout) :: avgAKv(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   ifdef AVERAGES_AKT
      real(r8), intent(inout) :: avgAKt(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   ifdef AVERAGES_AKS
      real(r8), intent(inout) :: avgAKs(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(inout) :: avgstf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgswf(LBi:UBi,LBj:UBj)
#    ifdef BULK_FLUXES
      real(r8), intent(inout) :: avglhf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgshf(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avglrf(LBi:UBi,LBj:UBj)
#     ifdef EMINUSP
      real(r8), intent(inout) :: avgevap(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgrain(LBi:UBi,LBj:UBj)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: avgsrf(LBi:UBi,LBj:UBj)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(inout) :: avgsus(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgsvs(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgbus(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgbvs(LBi:UBi,LBj:UBj)
#  endif
#  ifdef AVERAGES_QUADRATIC
#   ifdef SOLVE3D
      real(r8), intent(inout) :: avgHuon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgHvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgUU(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgUV(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgVV(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgHuonT(LBi:UBi,LBj:UBj,N(ng),NAT)
      real(r8), intent(inout) :: avgHvomT(LBi:UBi,LBj:UBj,N(ng),NAT)
      real(r8), intent(inout) :: avgUT(LBi:UBi,LBj:UBj,N(ng),NAT)
      real(r8), intent(inout) :: avgVT(LBi:UBi,LBj:UBj,N(ng),NAT)
      real(r8), intent(inout) :: avgTT(LBi:UBi,LBj:UBj,N(ng),NAT)
#   endif     
      real(r8), intent(inout) :: avgU2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgV2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgZZ(LBi:UBi,LBj:UBj)
#  endif
#  ifdef AVERAGES_NEARSHORE
      real(r8), intent(inout) :: avgu2Sd(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgv2Sd(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgu2RS(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgv2RS(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgSxx2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgSxy2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgSyy2d(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: avgu2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgv2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgzeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      integer :: mk, nk
      integer, dimension(2*NTC+1) :: indx
# endif

      real(r8) :: fac, fac1
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      real(r8) :: d, ubar_sum, vbar_sum, zeta_sum
#  ifdef SOLVE3D
      real(r8) :: u_sum, v_sum
#  endif
      real(r8), dimension(0:2*NTC) :: Ak_zeta
      real(r8), dimension(0:2*NTC) :: Ak_ubar
      real(r8), dimension(0:2*NTC) :: Ak_vbar
#  ifdef SOLVE3D
      real(r8), dimension(0:2*NTC) :: Ak_u
      real(r8), dimension(0:2*NTC) :: Ak_v
#  endif
      real(r8), dimension(0:2*NTC) :: tide_harmonics
      real(r8), dimension(0:2*NTC,0:2*NTC) :: C, Y
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Compute least-squares coefficients to detide time-averaged fields.
!  Notice that the coefficients are always accumulated and not
!  re-initialized.  This allows better tidal fit as the simulation
!  progresses.
!
        IF (SOUTH_WEST_TEST) THEN
          Hcount(ng)=Hcount(ng)+1
          DO nk=1,NTC
            SinW_avg(nk)=SinOmega(nk)
            CosW_avg(nk)=CosOmega(nk)
            SinW_sum(nk)=SinW_sum(nk)+SinOmega(nk)
            CosW_sum(nk)=CosW_sum(nk)+CosOmega(nk)
            DO mk=1,NTC
              SinWSinW(mk,nk)=SinWSinW(mk,nk)+SinOmega(mk)*SinOmega(nk)
              CosWCosW(mk,nk)=CosWCosW(mk,nk)+CosOmega(mk)*CosOmega(nk)
              SinWCosW(mk,nk)=SinWCosW(mk,nk)+SinOmega(mk)*CosOmega(nk)
            END DO
          END DO
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinOmega(nk)
            tide_harmonics(nk+NTC)=CosOmega(nk)
          END DO
        END IF
# endif
!
!  Initialize 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
            DO nk=0,2*NTC
              zeta_tide(i,j,nk)=zeta_tide(i,j,nk)+                      &
     &                          zeta(i,j,Kout)*tide_harmonics(nk)
              ubar_tide(i,j,nk)=ubar_tide(i,j,nk)+                      &
     &                          ubar(i,j,Kout)*tide_harmonics(nk)
              vbar_tide(i,j,nk)=vbar_tide(i,j,nk)+                      &
     &                          vbar(i,j,Kout)*tide_harmonics(nk)
            END DO
# endif
            avgzeta(i,j)=zeta(i,j,Kout)
            avgu2d (i,j)=ubar(i,j,Kout)
            avgv2d (i,j)=vbar(i,j,Kout)
# if defined FORWARD_WRITE && defined SOLVE3D
            avgDU_avg1(i,j)=DU_avg1(i,j)
            avgDU_avg2(i,j)=DU_avg2(i,j)
            avgDV_avg1(i,j)=DV_avg1(i,j)
            avgDV_avg2(i,j)=DV_avg2(i,j)
# endif
# ifdef AVERAGES_NEARSHORE
            avgu2Sd(i,j)=ubar_stokes(i,j)
            avgv2Sd(i,j)=vbar_stokes(i,j)
            avgu2RS(i,j)=rustr2d(i,j)
            avgv2RS(i,j)=rvstr2d(i,j)
            avgSxx2d(i,j)=Sxx_bar(i,j)
            avgSxy2d(i,j)=Sxy_bar(i,j)
            avgSyy2d(i,j)=Syy_bar(i,j)
# endif
# ifdef AVERAGES_QUADRATIC
            avgZZ(i,j)=zeta(i,j,Kout)*zeta(i,j,Kout)
            avgU2(i,j)=ubar(i,j,Kout)*ubar(i,j,Kout)
            avgV2(i,j)=vbar(i,j,Kout)*vbar(i,j,Kout)
# endif
# ifdef SOLVE3D
#  ifdef LMD_SKPP
            avghsbl(i,j)=hsbl(i,j)
#  endif
#  ifdef LMD_BKPP
            avghbbl(i,j)=hbbl(i,j)
#  endif
#  ifdef AVERAGES_FLUXES
            avgstf(i,j)=stflx(i,j,itemp)
            avgswf(i,j)=stflx(i,j,isalt)
#   ifdef BULK_FLUXES
            avglhf(i,j)=lhflx(i,j)
            avgshf(i,j)=shflx(i,j)
            avglrf(i,j)=lrflx(i,j)
#    ifdef EMINUSP
            avgevap(i,j)=evap(i,j)
            avgrain(i,j)=rain(i,j)
#    endif
#   endif
#   ifdef SHORTWAVE
            avgsrf(i,j)=srflx(i,j)
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
            avgsus(i,j)=sustr(i,j)
            avgsvs(i,j)=svstr(i,j)
            avgbus(i,j)=bustr(i,j)
            avgbvs(i,j)=bvstr(i,j)
# endif
          END DO
        END DO

# ifdef SOLVE3D
!
!  Initialize fields associated with 3D horizontal momentum.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
              DO nk=0,2*NTC
                u_tide(i,j,k,nk)=u_tide(i,j,k,nk)+                      &
     &                           u(i,j,k,Nout)*tide_harmonics(nk)
                v_tide(i,j,k,nk)=v_tide(i,j,k,nk)+                      &
     &                           v(i,j,k,Nout)*tide_harmonics(nk)
              END DO
#  endif
              avgu3d(i,j,k)=u(i,j,k,Nout)
              avgv3d(i,j,k)=v(i,j,k,Nout)
              avgrho(i,j,k)=rho(i,j,k)
#  ifdef AVERAGES_NEARSHORE
              avgu3Sd(i,j,k)=u_stokes(i,j,k)
              avgv3Sd(i,j,k)=v_stokes(i,j,k)
              avgu3RS(i,j,k)=rustr3d(i,j,k)
              avgv3RS(i,j,k)=rvstr3d(i,j,k)
              avgSxx3d(i,j,k)=Sxx(i,j,k)
              avgSxy3d(i,j,k)=Sxy(i,j,k)
              avgSyy3d(i,j,k)=Syy(i,j,k)
              avgSzx3d(i,j,k)=Szx(i,j,k)
              avgSzy3d(i,j,k)=Szy(i,j,k)
#  endif
#  ifdef AVERAGES_QUADRATIC
              avgHuon(i,j,k)=Huon(i,j,k)
              avgHvom(i,j,k)=Hvom(i,j,k)
              avgUU(i,j,k)=u(i,j,k,Nout)*u(i,j,k,Nout)
              avgVV(i,j,k)=v(i,j,k,Nout)*v(i,j,k,Nout)
#  endif
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          DO j=Jstr,Jend
            DO i=Istr,Iend
              avgUV(i,j,k)=0.25_r8*(u(i  ,j  ,k,Nout)+                  &
     &                              u(i+1,j  ,k,Nout))*                 &
     &                             (v(i  ,j  ,k,Nout)+                  &
     &                              v(i  ,j+1,k,Nout))
            END DO
          END DO
#  endif
        END DO
!
!  Initialized fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                avgt(i,j,k,itrc)=t(i,j,k,Nout,itrc)
              END DO
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          IF (itrc.le.NAT) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  avgTT(i,j,k,itrc)=t(i,j,k,Nout,itrc)*                 &
     &                              t(i,j,k,Nout,itrc)
                END DO
                DO i=Istr,Iend
                  avgUT(i,j,k,itrc)=0.5_r8*u(i,j,k,Nout)*               &
     &                              (t(i-1,j,k,Nout,itrc)+              &
     &                               t(i  ,j,k,Nout,itrc))
                  avgHuonT(i,j,k,itrc)=0.5_r8*Huon(i,j,k)*              &
     &                                 (t(i-1,j,k,Nout,itrc)+           &
     &                                  t(i  ,j,k,Nout,itrc))
                END DO
              END DO
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  avgVT(i,j,k,itrc)=0.5_r8*v(i,j,k,Nout)*               &
     &                              (t(i,j-1,k,Nout,itrc)+              &
     &                               t(i,j  ,k,Nout,itrc))
                  avgHvomT(i,j,k,itrc)=0.5_r8*Hvom(i,j,k)*              &
     &                                 (t(i,j-1,k,Nout,itrc)+           &
     &                                  t(i,j  ,k,Nout,itrc))
                END DO
              END DO
            END DO
          END IF
#  endif
        END DO
!
!  Initialize fields at W-points.
!
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgw3d(i,j,k)=W(i,j,k)*pm(i,j)*pn(i,j)
              avgwvel(i,j,k)=wvel(i,j,k)
#  ifdef AVERAGES_AKV
              avgAKv(i,j,k)=Akv(i,j,k)
#  endif
#  ifdef AVERAGES_AKT
              avgAKt(i,j,k)=Akt(i,j,k,itemp)
#  endif
#  ifdef AVERAGES_AKS
              avgAKs(i,j,k)=Akt(i,j,k,isalt)
#  endif
            END DO
          END DO
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Accumukate Detide least-squares coefficients. They only vary in time
!  since omega (as computed in set_tides) used model time coordinate.
!
        IF (SOUTH_WEST_TEST) THEN
          Hcount(ng)=Hcount(ng)+1
          DO nk=1,NTC
            SinW_avg(nk)=SinW_avg(nk)+SinOmega(nk)
            CosW_avg(nk)=CosW_avg(nk)+CosOmega(nk)
            SinW_sum(nk)=SinW_sum(nk)+SinOmega(nk)
            CosW_sum(nk)=CosW_sum(nk)+CosOmega(nk)
            DO mk=1,NTC
              SinWSinW(mk,nk)=SinWSinW(mk,nk)+SinOmega(mk)*SinOmega(nk)
              CosWCosW(mk,nk)=CosWCosW(mk,nk)+CosOmega(mk)*CosOmega(nk)
              SinWCosW(mk,nk)=SinWCosW(mk,nk)+SinOmega(mk)*CosOmega(nk)
            END DO
          END DO
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinOmega(nk)
            tide_harmonics(nk+NTC)=CosOmega(nk)
          END DO
        END IF
# endif
!
!  Accumulate 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
            DO nk=0,2*NTC
              zeta_tide(i,j,nk)=zeta_tide(i,j,nk)+                      &
     &                          zeta(i,j,Kout)*tide_harmonics(nk)
              ubar_tide(i,j,nk)=ubar_tide(i,j,nk)+                      &
     &                          ubar(i,j,Kout)*tide_harmonics(nk)
              vbar_tide(i,j,nk)=vbar_tide(i,j,nk)+                      &
     &                          vbar(i,j,Kout)*tide_harmonics(nk)
            END DO
# endif
            avgzeta(i,j)=avgzeta(i,j)+zeta(i,j,Kout)
            avgu2d (i,j)=avgu2d (i,j)+ubar(i,j,Kout)
            avgv2d (i,j)=avgv2d (i,j)+vbar(i,j,Kout)
# if defined FORWARD_WRITE && defined SOLVE3D
            avgDU_avg1(i,j)=avgDU_avg1(i,j)+DU_avg1(i,j)
            avgDU_avg2(i,j)=avgDU_avg2(i,j)+DU_avg2(i,j)
            avgDV_avg1(i,j)=avgDV_avg1(i,j)+DV_avg1(i,j)
            avgDV_avg2(i,j)=avgDV_avg2(i,j)+DV_avg2(i,j)
# endif
# ifdef AVERAGES_NEARSHORE
            avgu2Sd(i,j)=avgu2Sd(i,j)+ubar_stokes(i,j)
            avgv2Sd(i,j)=avgv2Sd(i,j)+vbar_stokes(i,j)
            avgu2RS(i,j)=avgu2RS(i,j)+rustr2d(i,j)
            avgv2RS(i,j)=avgv2RS(i,j)+rvstr2d(i,j)
            avgSxx2d(i,j)=avgSxx2d(i,j)+Sxx_bar(i,j)
            avgSxy2d(i,j)=avgSxy2d(i,j)+Sxy_bar(i,j)
            avgSyy2d(i,j)=avgSyy2d(i,j)+Syy_bar(i,j)
# endif
# ifdef AVERAGES_QUADRATIC
            avgZZ(i,j)=avgZZ(i,j)+zeta(i,j,Kout)*zeta(i,j,Kout)
            avgU2(i,j)=avgU2(i,j)+ubar(i,j,Kout)*ubar(i,j,Kout)
            avgV2(i,j)=avgU2(i,j)+vbar(i,j,Kout)*vbar(i,j,Kout)
# endif
# ifdef SOLVE3D
#  ifdef LMD_SKPP
            avghsbl(i,j)=avghsbl(i,j)+hsbl(i,j)
#  endif
#  ifdef LMD_BKPP
            avghbbl(i,j)=avghbbl(i,j)+hbbl(i,j)
#  endif
#  ifdef AVERAGES_FLUXES
            avgstf(i,j)=avgstf(i,j)+stflx(i,j,itemp)
            avgswf(i,j)=avgswf(i,j)+stflx(i,j,isalt)
#   ifdef BULK_FLUXES
            avglhf(i,j)=avglhf(i,j)+lhflx(i,j)
            avgshf(i,j)=avgshf(i,j)+shflx(i,j)
            avglrf(i,j)=avglrf(i,j)+lrflx(i,j)
#    ifdef EMINUSP
            avgevap(i,j)=avgevap(i,j)+evap(i,j)
            avgrain(i,j)=avgrain(i,j)+rain(i,j)
#    endif
#   endif
#   ifdef SHORTWAVE
            avgsrf(i,j)=avgsrf(i,j)+srflx(i,j)
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
            avgsus(i,j)=avgsus(i,j)+sustr(i,j)
            avgsvs(i,j)=avgsvs(i,j)+svstr(i,j)
            avgbus(i,j)=avgbus(i,j)+bustr(i,j)
            avgbvs(i,j)=avgbvs(i,j)+bvstr(i,j)
# endif
          END DO
        END DO

# ifdef SOLVE3D
!
!  Accumulate fields associated with 3D horizontal momentum.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
              DO nk=0,2*NTC
                u_tide(i,j,k,nk)=u_tide(i,j,k,nk)+                      &
     &                           u(i,j,k,Nout)*tide_harmonics(nk)
                v_tide(i,j,k,nk)=v_tide(i,j,k,nk)+                      &
     &                           v(i,j,k,Nout)*tide_harmonics(nk)
              END DO
#  endif
              avgu3d(i,j,k)=avgu3d(i,j,k)+u(i,j,k,Nout)
              avgv3d(i,j,k)=avgv3d(i,j,k)+v(i,j,k,Nout)
              avgrho(i,j,k)=avgrho(i,j,k)+rho(i,j,k)
#  ifdef AVERAGES_NEARSHORE
              avgu3Sd(i,j,k)=avgu3Sd(i,j,k)+u_stokes(i,j,k)
              avgv3Sd(i,j,k)=avgv3Sd(i,j,k)+v_stokes(i,j,k)
              avgu3RS(i,j,k)=avgu3RS(i,j,k)+rustr3d(i,j,k)
              avgv3RS(i,j,k)=avgv3RS(i,j,k)+rvstr3d(i,j,k)
              avgSxx3d(i,j,k)=avgSxx3d(i,j,k)+Sxx(i,j,k)
              avgSxy3d(i,j,k)=avgSxy3d(i,j,k)+Sxy(i,j,k)
              avgSyy3d(i,j,k)=avgSyy3d(i,j,k)+Syy(i,j,k)
              avgSzx3d(i,j,k)=avgSzx3d(i,j,k)+Szx(i,j,k)
              avgSzy3d(i,j,k)=avgSzy3d(i,j,k)+Szy(i,j,k)
#  endif
#  ifdef AVERAGES_QUADRATIC
              avgHuon(i,j,k)=avgHuon(i,j,k)+Huon(i,j,k)
              avgHvom(i,j,k)=avgHvom(i,j,k)+Hvom(i,j,k)
              avgUU(i,j,k)=avgUU(i,j,k)+u(i,j,k,Nout)*u(i,j,k,Nout)
              avgVV(i,j,k)=avgVV(i,j,k)+v(i,j,k,Nout)*v(i,j,k,Nout)
#  endif
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          DO j=Jstr,Jend
            DO i=Istr,Iend
              avgUV(i,j,k)=avgUV(i,j,k)+                                &
     &                     0.25_r8*(u(i  ,j  ,k,Nout)+                  &
     &                              u(i+1,j  ,k,Nout))*                 &
     &                             (v(i  ,j  ,k,Nout)+                  &
     &                              v(i  ,j+1,k,Nout))
            END DO
          END DO
#  endif
        END DO
!
!  Accumulate fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                avgt(i,j,k,itrc)=avgt(i,j,k,itrc)+t(i,j,k,Nout,itrc)
              END DO
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          IF (itrc.le.NAT) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  avgTT(i,j,k,itrc)=avgTT(i,j,k,itrc)+                  &
     &                              t(i,j,k,Nout,itrc)*                 &
     &                              t(i,j,k,Nout,itrc)
                END DO
                DO i=Istr,Iend
                  avgUT(i,j,k,itrc)=avgUT(i,j,k,itrc)+                  &
     &                              0.5_r8*u(i,j,k,Nout)*               &
     &                              (t(i-1,j,k,Nout,itrc)+              &
     &                               t(i  ,j,k,Nout,itrc))
                  avgHuonT(i,j,k,itrc)=avgHuonT(i,j,k,itrc)+            &
     &                                 0.5_r8*Huon(i,j,k)*              &
     &                                 (t(i-1,j,k,Nout,itrc)+           &
     &                                  t(i  ,j,k,Nout,itrc))
                END DO
              END DO
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  avgVT(i,j,k,itrc)=avgVT(i,j,k,itrc)+                  &
     &                              0.5_r8*v(i,j,k,Nout)*               &
     &                              (t(i,j-1,k,Nout,itrc)+              &
     &                               t(i,j  ,k,Nout,itrc))
                  avgHvomT(i,j,k,itrc)=avgVT(i,j,k,itrc)+               &
     &                                 0.5_r8*Hvom(i,j,k)*              &
     &                                 (t(i,j-1,k,Nout,itrc)+           &
     &                                  t(i,j  ,k,Nout,itrc))
                END DO
              END DO
            END DO
          END IF
#  endif
        END DO
!
!  Accumulate fields at W-points.
!
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgw3d(i,j,k)=avgw3d(i,j,k)+W(i,j,k)*pm(i,j)*pn(i,j)
              avgwvel(i,j,k)=avgwvel(i,j,k)+wvel(i,j,k)
#  ifdef AVERAGES_AKV
              avgAKv(i,j,k)=avgAKv(i,j,k)+Akv(i,j,k)
#  endif
#  ifdef AVERAGES_AKT
              avgAKt(i,j,k)=avgAKt(i,j,k)+Akt(i,j,k,itemp)
#  endif
#  ifdef AVERAGES_AKS
              avgAKs(i,j,k)=avgAKs(i,j,k)+Akt(i,j,k,isalt)
#  endif
            END DO
          END DO
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        fac=1.0_r8/REAL(nAVG(ng),r8)
        IF (SOUTH_WEST_TEST) THEN
          AVGtime(ng)=AVGtime(ng)+REAL(nAVG(ng),r8)*dt(ng)
        END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Compute detide least-squares coefficients.  Build coefficient squared
!  matrix C(0:2*NTC,0:2*NTC) to invert. It is 2*NTC because we are
!  solving for real components Ak and Bk. The zero rows and column is
!  for the coefficients associated with the time mean.
!
!        F(t) = Fmean + SUM [ Ak sin(omega(k)*t) ]
!                     + SUM [ Bk cos(omega(k)*t) ]   for k=1:NTC
!
!  In the code below, all the arrays are collapsed into a single
!  dimension index such that:
!
!            k=0               mean term
!            k=1:NTC           sine terms
!            k=NTC+1:2*NTC     cosine terms
!
        IF (SOUTH_WEST_TEST) THEN
          C(0,0)=1.0_r8                     ! time-averaged coefficient
          fac1=1.0_r8/REAL(Hcount(ng),r8)   ! global summation factor
          DO nk=1,NTC
            C(0,nk    )=fac1*SinW_sum(nk)
            C(0,nk+NTC)=fac1*CosW_sum(nk)
            C(nk,0    )=C(0,nk)             ! symmetric
            C(nk+NTC,0)=C(0,nk+NTC)         ! symmetric
            DO mk=1,NTC
              C(mk,nk)=fac1*SinWSinW(mk,nk)
              C(mk,nk+NTC)=fac1*SinWCosW(mk,nk)
              C(mk+NTC,nk)=fac1*SinWCosW(nk,mk)
              C(mk+NTC,nk+NTC)=fac1*CosWCosW(mk,nk)
            END DO
          END DO
          DO nk=0,2*NTC
            DO mk=0,2*NTC
              C(nk,mk)=C(mk,nk)
            END DO
          END DO
!
!  Invert least-squares coefficient matrix by LU decomposition.
!
          DO mk=0,2*NTC
            DO nk=0,2*NTC
              Y(mk,nk)=0.0_r8
            END DO
            Y(mk,mk)=1.0_r8                ! identity matrix
          END DO
          CALL ludcmp (C(0,0), 2*NTC+1, 2*NTC+1, indx, d)
!
!  Find inverse by columns. The matrix Y will now contain the inverse
!  of the least-squares coefficient matrix C, which will have been
!  destroyed.
!
          DO nk=0,2*NTC
            CALL lubksb (C(0,0), 2*NTC+1, 2*NTC+1, indx, Y(0,nk))
          END DO
!
!  Compute time-averaged harmonics for current field average window.
!
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=fac*SinW_avg(nk)
            tide_harmonics(nk+NTC)=fac*cosW_avg(nk)
          END DO
!
!  Scale inverse by the global summation factor.
!
          DO nk=0,2*NTC
            DO mk=0,2*NTC
              Y(nk,mk)=fac1*Y(nk,mk)
            END DO
          END DO
        END IF
# endif
!
!  Process 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
            zeta_sum=0.0_r8
            ubar_sum=0.0_r8
            vbar_sum=0.0_r8
            DO nk=0,2*NTC
              Ak_zeta(nk)=0.0_r8
              Ak_ubar(nk)=0.0_r8
              Ak_vbar(nk)=0.0_r8
              DO mk=0,2*NTC
!!              Ak_zeta(nk)=Ak_zeta(nk)+Y(nk,mk)*zeta_tide(i,j,mk)
!!              Ak_ubar(nk)=Ak_ubar(nk)+Y(nk,mk)*ubar_tide(i,j,mk)
!!              Ak_vbar(nk)=Ak_vbar(nk)+Y(nk,mk)*vbar_tide(i,j,mk)
                Ak_zeta(nk)=Ak_zeta(nk)+Y(mk,nk)*zeta_tide(i,j,mk)
                Ak_ubar(nk)=Ak_ubar(nk)+Y(mk,nk)*ubar_tide(i,j,mk)
                Ak_vbar(nk)=Ak_vbar(nk)+Y(mk,nk)*vbar_tide(i,j,mk)
              END DO
              zeta_sum=zeta_sum+Ak_zeta(nk)*tide_harmonics(nk)
              ubar_sum=ubar_sum+Ak_ubar(nk)*tide_harmonics(nk)
              vbar_sum=vbar_sum+Ak_vbar(nk)*tide_harmonics(nk)
            END DO
            avgzeta(i,j)=fac*avgzeta(i,j)
            avgu2d (i,j)=fac*avgu2d (i,j)
            avgv2d (i,j)=fac*avgv2d (i,j)
            IF (.TRUE.) THEN 
              zeta_detided(i,j)=avgzeta(i,j)-zeta_sum
              ubar_detided(i,j)=avgu2d (i,j)-ubar_sum
              vbar_detided(i,j)=avgv2d (i,j)-vbar_sum
            ELSE
              zeta_detided(i,j)=avgzeta(i,j)
              ubar_detided(i,j)=avgu2d (i,j)
              vbar_detided(i,j)=avgv2d (i,j)
            END IF
# else
            avgzeta(i,j)=fac*avgzeta(i,j)
            avgu2d (i,j)=fac*avgu2d (i,j)
            avgv2d (i,j)=fac*avgv2d (i,j)
# endif
# if defined FORWARD_WRITE && defined SOLVE3D
            avgDU_avg1(i,j)=fac*avgDU_avg1(i,j)
            avgDU_avg2(i,j)=fac*avgDU_avg2(i,j)
            avgDV_avg1(i,j)=fac*avgDV_avg1(i,j)
            avgDV_avg2(i,j)=fac*avgDV_avg2(i,j)
# endif
# ifdef AVERAGES_NEARSHORE
            avgu2Sd(i,j)=fac*avgu2Sd(i,j)
            avgv2Sd(i,j)=fac*avgv2Sd(i,j)
            avgu2RS(i,j)=fac*avgu2RS(i,j)
            avgv2RS(i,j)=fac*avgv2RS(i,j)
            avgSxx2d(i,j)=fac*avgSxx2d(i,j)
            avgSxy2d(i,j)=fac*avgSxy2d(i,j)
            avgSyy2d(i,j)=fac*avgSyy2d(i,j)
# endif
# ifdef AVERAGES_QUADRATIC
            avgZZ(i,j)=fac*avgZZ(i,j)
            avgU2(i,j)=fac*avgU2(i,j)
            avgV2(i,j)=fac*avgU2(i,j)
# endif
# ifdef SOLVE3D
#  ifdef LMD_SKPP
            avghsbl(i,j)=fac*avghsbl(i,j)
#  endif
#  ifdef LMD_BKPP
            avghbbl(i,j)=fac*avghbbl(i,j)
#  endif
#  ifdef AVERAGES_FLUXES
            avgstf(i,j)=fac*avgstf(i,j)
            avgswf(i,j)=fac*avgswf(i,j)
#   ifdef BULK_FLUXES
            avglhf(i,j)=fac*avglhf(i,j)
            avgshf(i,j)=fac*avgshf(i,j)
            avglrf(i,j)=fac*avglrf(i,j)
#    ifdef EMINUSP
            avgevap(i,j)=fac*avgevap(i,j)
            avgrain(i,j)=fac*avgrain(i,j)
#    endif
#   endif
#   ifdef SHORTWAVE
            avgsrf(i,j)=fac*avgsrf(i,j)
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
            avgsus(i,j)=fac*avgsus(i,j)
            avgsvs(i,j)=fac*avgsvs(i,j)
            avgbus(i,j)=fac*avgbus(i,j)
            avgbvs(i,j)=fac*avgbvs(i,j)
# endif
          END DO
        END DO

# ifdef SOLVE3D
!
!  Process fields associated with 3D horizontal momentum.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
              u_sum=0.0_r8
              v_sum=0.0_r8
              DO nk=0,2*NTC
                Ak_u(nk)=0.0_r8
                Ak_v(nk)=0.0_r8
                DO mk=0,2*NTC
!!                Ak_u(nk)=Ak_u(nk)+Y(nk,mk)*u_tide(i,j,k,mk)
!!                Ak_v(nk)=Ak_v(nk)+Y(nk,mk)*v_tide(i,j,k,mk)
                  Ak_u(nk)=Ak_u(nk)+Y(mk,nk)*u_tide(i,j,k,mk)
                  Ak_v(nk)=Ak_v(nk)+Y(mk,nk)*v_tide(i,j,k,mk)
                END DO
                u_sum=u_sum+Ak_u(nk)*tide_harmonics(nk)
                v_sum=v_sum+Ak_v(nk)*tide_harmonics(nk)
              END DO
              avgu3d(i,j,k)=fac*avgu3d(i,j,k)
              avgv3d(i,j,k)=fac*avgv3d(i,j,k)
              IF (.TRUE.) THEN
                u_detided(i,j,k)=avgu3d(i,j,k)-u_sum
                v_detided(i,j,k)=avgv3d(i,j,k)-v_sum
              ELSE
                u_detided(i,j,k)=avgu3d(i,j,k)
                v_detided(i,j,k)=avgv3d(i,j,k)
              END IF
#  else
              avgu3d(i,j,k)=fac*avgu3d(i,j,k)
              avgv3d(i,j,k)=fac*avgv3d(i,j,k)
#  endif
              avgrho(i,j,k)=fac*avgrho(i,j,k)
#  ifdef AVERAGES_NEARSHORE
              avgu3Sd(i,j,k)=fac*avgu3Sd(i,j,k)
              avgv3Sd(i,j,k)=fac*avgv3Sd(i,j,k)
              avgu3RS(i,j,k)=fac*avgu3RS(i,j,k)
              avgv3RS(i,j,k)=fac*avgv3RS(i,j,k)
              avgSxx3d(i,j,k)=fac*avgSxx3d(i,j,k)
              avgSxy3d(i,j,k)=fac*avgSxy3d(i,j,k)
              avgSyy3d(i,j,k)=fac*avgSyy3d(i,j,k)
              avgSzx3d(i,j,k)=fac*avgSzx3d(i,j,k)
              avgSzy3d(i,j,k)=fac*avgSzy3d(i,j,k)
#  endif
#  ifdef AVERAGES_QUADRATIC
              avgHuon(i,j,k)=fac*avgHuon(i,j,k)
              avgHvom(i,j,k)=fac*avgHvom(i,j,k)
              avgUU(i,j,k)=fac*avgUU(i,j,k)
              avgVV(i,j,k)=fac*avgVV(i,j,k)
#  endif
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          DO j=Jstr,Jend
            DO i=Istr,Iend
              avgUV(i,j,k)=fac*avgUV(i,j,k)
            END DO
          END DO
#  endif
        END DO
!
!  Process fields associated with tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                avgt(i,j,k,itrc)=fac*avgt(i,j,k,itrc)
              END DO
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          IF (itrc.le.NAT) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  avgTT(i,j,k,itrc)=fac*avgTT(i,j,k,itrc)
                END DO
                DO i=Istr,Iend
                  avgUT(i,j,k,itrc)=fac*avgUT(i,j,k,itrc)
                  avgHuonT(i,j,k,itrc)=fac*avgHuonT(i,j,k,itrc)
                END DO
              END DO
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  avgVT(i,j,k,itrc)=fac*avgVT(i,j,k,itrc)
                  avgHvomT(i,j,k,itrc)=fac*avgHvomT(i,j,k,itrc)
                END DO
              END DO
            END DO
          END IF
#  endif
        END DO
!
!  Process fields at W-points.
!
        DO k=0,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              avgw3d(i,j,k)=fac*avgw3d(i,j,k)
              avgwvel(i,j,k)=fac*avgwvel(i,j,k)
#  ifdef AVERAGES_AKV
              avgAKv(i,j,k)=fac*avgAKv(i,j,k)
#  endif
#  ifdef AVERAGES_AKT
              avgAKt(i,j,k)=fac*avgAKt(i,j,k)
#  endif
#  ifdef AVERAGES_AKS
              avgAKs(i,j,k)=fac*avgAKs(i,j,k)
#  endif
            END DO
          END DO
        END DO
# endif
      END IF

      RETURN
      END SUBROUTINE set_avg_tile
#endif
      END MODULE set_avg_mod
