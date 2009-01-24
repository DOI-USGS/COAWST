#include "cppdefs.h"
      MODULE ad_set_avg_mod
#if defined AVERAGES && defined ADJOINT
!
!svn $Id: ad_set_avg.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  adjoint fields.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: ad_set_avg

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_average
      USE mod_coupling
      USE mod_forces
# ifdef SOLVE3D
      USE mod_grid
      USE mod_mixing
# endif
      USE mod_ocean
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
      CALL wclock_on (ng, iADM, 5)
# endif
      CALL ad_set_avg_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
# ifdef SOLVE3D
     &                      nstp(ng),                                   &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      OCEAN(ng) % ad_u,                           &
     &                      OCEAN(ng) % ad_v,                           &
     &                      OCEAN(ng) % ad_W,                           &
     &                      OCEAN(ng) % ad_t,                           &
     &                      OCEAN(ng) % ad_rho,                         &
#  ifdef LMD_SKPP_NOT_YET
     &                      MIXING(ng) % ad_hsbl,                       &
#  endif
#  ifdef LMD_BKPP_NOT_YET
     &                      MIXING(ng) % ad_hbbl,                       &
#  endif
#  ifdef AVERAGES_AKV
     &                      MIXING(ng) % ad_Akv,                        &
#  endif
#  if defined AVERAGES_AKT || defined AVERAGES_AKS
     &                      MIXING(ng) % ad_Akt,                        &
#  endif
#  ifdef AVERAGES_FLUXES
     &                      FORCES(ng) % ad_stflx,                      &
#   ifdef BULK_FLUXES
     &                      FORCES(ng) % ad_lhflx,                      &
     &                      FORCES(ng) % ad_shflx,                      &
     &                      FORCES(ng) % ad_lrflx,                      &
#    ifdef EMINUSP
     &                      FORCES(ng) % ad_evap,                       &
!!   &                      FORCES(ng) % ad_rain,                       &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                      FORCES(ng) % ad_srflx,                      &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                      FORCES(ng) % ad_sustr,                      &
     &                      FORCES(ng) % ad_svstr,                      &
     &                      FORCES(ng) % ad_bustr,                      &
     &                      FORCES(ng) % ad_bvstr,                      &
# endif
     &                      OCEAN(ng) % ad_ubar_sol,                    &
     &                      OCEAN(ng) % ad_vbar_sol,                    &
     &                      OCEAN(ng) % ad_zeta_sol,                    &
# ifdef SOLVE3D
     &                      AVERAGE(ng) % avgu3d,                       &
     &                      AVERAGE(ng) % avgv3d,                       &
     &                      AVERAGE(ng) % avgw3d,                       &
     &                      AVERAGE(ng) % avgt,                         &
     &                      AVERAGE(ng) % avgrho,                       &
#  ifdef LMD_SKPP_NOT_YET
     &                      AVERAGE(ng) % avghsbl,                      &
#  endif
#  ifdef LMD_BKPP_NOT_YET
     &                      AVERAGE(ng) % avghbbl,                      &
#  endif
#  ifdef AVERAGES_AKV
     &                      AVERAGE(ng) % avgAKv,                       &
#  endif
#  ifdef AVERAGES_AKT
     &                      AVERAGE(ng) % avgAKt,                       &
#  endif
#  ifdef AVERAGES_AKS
     &                      AVERAGE(ng) % avgAKs,                       &
#  endif
#  ifdef AVERAGES_FLUXES
     &                      AVERAGE(ng) % avgstf,                       &
     &                      AVERAGE(ng) % avgswf,                       &
#   ifdef BULK_FLUXES
     &                      AVERAGE(ng) % avglhf,                       &
     &                      AVERAGE(ng) % avgshf,                       &
     &                      AVERAGE(ng) % avglrf,                       &
#    ifdef EMINUSP
     &                      AVERAGE(ng) % avgevap,                      &
!!   &                      AVERAGE(ng) % avgrain,                      &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                      AVERAGE(ng) % avgsrf,                       &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                      AVERAGE(ng) % avgsus,                       &
     &                      AVERAGE(ng) % avgsvs,                       &
     &                      AVERAGE(ng) % avgbus,                       &
     &                      AVERAGE(ng) % avgbvs,                       &
# endif
# ifdef AVERAGES_QUADRATIC
#  ifdef SOLVE3D
     &                      AVERAGE(ng) % avgUU,                        &
     &                      AVERAGE(ng) % avgUV,                        &
     &                      AVERAGE(ng) % avgVV,                        &
     &                      AVERAGE(ng) % avgUT,                        &
     &                      AVERAGE(ng) % avgVT,                        &
     &                      AVERAGE(ng) % avgTT,                        &
#  endif     
     &                      AVERAGE(ng) % avgU2,                        &
     &                      AVERAGE(ng) % avgV2,                        &
     &                      AVERAGE(ng) % avgZZ,                        &
# endif
     &                      AVERAGE(ng) % avgu2d,                       &
     &                      AVERAGE(ng) % avgv2d,                       &
     &                      AVERAGE(ng) % avgzeta)

# ifdef PROFILE
      CALL wclock_off (ng, iADM, 5)
# endif
      RETURN
      END SUBROUTINE ad_set_avg
!
!***********************************************************************
      SUBROUTINE ad_set_avg_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
# ifdef SOLVE3D
     &                            nstp,                                 &
     &                            pm, pn,                               &
     &                            ad_u, ad_v, ad_W,                     &
     &                            ad_t, ad_rho,                         &
#  ifdef LMD_SKPP_NOT_YET
     &                            ad_hsbl,                              &
#  endif
#  ifdef LMD_BKPP_NOT_YET
     &                            ad_hbbl,                              &
#  endif
#  ifdef AVERAGES_AKV
     &                            ad_Akv,                               &
#  endif
#  if defined AVERAGES_AKT || defined AVERAGES_AKS
     &                            ad_Akt,                               &
#  endif
#  ifdef AVERAGES_FLUXES
     &                            ad_stflx,                             &
#   ifdef BULK_FLUXES
     &                            ad_lhflx, ad_shflx, ad_lrflx,         &
#    ifdef EMINUSP
     &                            ad_evap,                              &
!!   &                            ad_rain,                              &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                            ad_srflx,                             &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                            ad_sustr, ad_svstr,                   &
     &                            ad_bustr, ad_bvstr,                   &
# endif
     &                            ad_ubar_sol, ad_vbar_sol,             &
     &                            ad_zeta_sol,                          &
# ifdef SOLVE3D
     &                            avgu3d, avgv3d, avgw3d,               &
     &                            avgt, avgrho,                         &
#  ifdef LMD_SKPP_NOT_YET
     &                            avghsbl,                              &
#  endif
#  ifdef LMD_BKPP_NOT_YET
     &                            avghbbl,                              &
#  endif
#  ifdef AVERAGES_AKV
     &                            avgAKv,                               &
#  endif
#  ifdef AVERAGES_AKT
     &                            avgAKt,                               &
#  endif
#  ifdef AVERAGES_AKS
     &                            avgAKs,                               &
#  endif
#  ifdef AVERAGES_FLUXES
     &                            avgstf, avgswf,                       &
#   ifdef BULK_FLUXES
     &                            avglhf, avgshf, avglrf,               &
#    ifdef EMINUSP
     &                            avgevap,                              &
!!   &                            avgrain,                              &
#    endif
#   endif
#   ifdef SHORTWAVE
     &                            avgsrf,                               &
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
     &                            avgsus, avgsvs, avgbus, avgbvs,       &
# endif
# ifdef AVERAGES_QUADRATIC
#  ifdef SOLVE3D
     &                            avgUU, avgUV, avgVV,                  &
     &                            avgUT, avgVT, avgTT,                  &
#  endif     
     &                            avgU2, avgV2, avgZZ,                  &
# endif
     &                            avgu2d, avgv2d, avgzeta)
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
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef SOLVE3D
      integer, intent(in) :: nstp
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_W(LBi:,LBj:,0:)
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_rho(LBi:,LBj:,:)
#   ifdef LMD_SKPP_NOT_YET
      real(r8), intent(in) :: ad_hsbl(LBi:,LBj:)
#   endif
#   ifdef LMD_BKPP_NOT_YET
      real(r8), intent(in) :: hbbl(LBi:,LBj:)
#   endif
#   ifdef AVERAGES_AKV
      real(r8), intent(in) :: ad_Akv(LBi:,LBj:,0:)
#   endif
#   if defined AVERAGES_AKT || defined AVERAGES_AKS
      real(r8), intent(in) :: ad_Akt(LBi:,LBj:,0:,:)
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: ad_stflx(LBi:,LBj:,:)
#    ifdef BULK_FLUXES
      real(r8), intent(in) :: ad_lhflx(LBi:,LBj:)
      real(r8), intent(in) :: ad_shflx(LBi:,LBj:)
      real(r8), intent(in) :: ad_lrflx(LBi:,LBj:)
#     ifdef EMINUSP
      real(r8), intent(in) :: ad_evap(LBi:,LBj:)
!!    real(r8), intent(in) :: ad_rain(LBi:,LBj:)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(in) :: ad_srflx(LBi:,LBj:)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(in) :: ad_svstr(LBi:,LBj:)
      real(r8), intent(in) :: ad_bustr(LBi:,LBj:)
      real(r8), intent(in) :: ad_bvstr(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ad_ubar_sol(LBi:,LBj:)
      real(r8), intent(in) :: ad_vbar_sol(LBi:,LBj:)
      real(r8), intent(in) :: ad_zeta_sol(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: avgu3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgv3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgw3d(LBi:,LBj:,0:)
      real(r8), intent(inout) :: avgt(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgrho(LBi:,LBj:,:)
#   ifdef LMD_SKPP_NOT_YET
      real(r8), intent(inout) :: avghsbl(LBi:,LBj:)
#   endif
#   ifdef LMD_BKPP_NOT_YET
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
!!    real(r8), intent(inout) :: avgrain(LBi:,LBj:)
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
      real(r8), intent(inout) :: avgUU(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgUV(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgVV(LBi:,LBj:,:)
      real(r8), intent(inout) :: avgUT(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgVT(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: avgTT(LBi:,LBj:,:,:)
#   endif     
      real(r8), intent(inout) :: avgU2(LBi:,LBj:)
      real(r8), intent(inout) :: avgV2(LBi:,LBj:)
      real(r8), intent(inout) :: avgZZ(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: avgu2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgv2d(LBi:,LBj:)
      real(r8), intent(inout) :: avgzeta(LBi:,LBj:)

# else

#  ifdef SOLVE3D
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_rho(LBi:UBi,LBj:UBj,N(ng))
#   ifdef LMD_SKPP_NOT_YET
      real(r8), intent(in) :: ad_hsbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef LMD_BKPP_NOT_YET
      real(r8), intent(in) :: ad_hbbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef AVERAGES_AKV
      real(r8), intent(in) :: ad_Akv(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined AVERAGES_AKT || defined AVERAGES_AKS
      real(r8), intent(in) :: ad_Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#   endif
#   ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
#    ifdef BULK_FLUXES
      real(r8), intent(in) :: ad_lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_lrflx(LBi:UBi,LBj:UBj)
#     ifdef EMINUSP
      real(r8), intent(in) :: ad_evap(LBi:UBi,LBj:UBj)
!!    real(r8), intent(in) :: ad_rain(LBi:UBi,LBj:UBj)
#     endif
#    endif
#    ifdef SHORTWAVE
      real(r8), intent(in) :: ad_srflx(LBi:UBi,LBj:UBj)
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
      real(r8), intent(in) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_bvstr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ad_ubar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_vbar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_zeta_sol(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: avgu3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgv3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgw3d(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: avgt(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: avgrho(LBi:UBi,LBj:UBj,N(ng))
#   ifdef LMD_SKPP_NOT_YET
      real(r8), intent(inout) :: avghsbl(LBi:UBi,LBj:UBj)
#   endif
#   ifdef LMD_BKPP_NOT_YET
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
!!    real(r8), intent(inout) :: avgrain(LBi:UBi,LBj:UBj)
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
      real(r8), intent(inout) :: avgUU(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgUV(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgVV(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: avgUT(LBi:UBi,LBj:UBj,N(ng),NAT)
      real(r8), intent(inout) :: avgVT(LBi:UBi,LBj:UBj,N(ng),NAT)
      real(r8), intent(inout) :: avgTT(LBi:UBi,LBj:UBj,N(ng),NAT)
#   endif     
      real(r8), intent(inout) :: avgU2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgV2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgZZ(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: avgu2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgv2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: avgzeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8) :: fac

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
      IF (((iic(ng).lt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng),nAVG(ng)).eq.0)).or.                            &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Initialize 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgzeta(i,j)=ad_zeta_sol(i,j)
            avgu2d (i,j)=ad_ubar_sol(i,j)
            avgv2d (i,j)=ad_vbar_sol(i,j)
# ifdef AVERAGES_QUADRATIC
            avgZZ(i,j)=ad_zeta_sol(i,j)*ad_zeta_sol(i,j)
            avgU2(i,j)=ad_ubar_sol(i,j)*ad_ubar_sol(i,j)
            avgV2(i,j)=ad_vbar_sol(i,j)*ad_vbar_sol(i,j)
# endif
# ifdef SOLVE3D
#  ifdef LMD_SKPP_NOT_YET
            avghsbl(i,j)=ad_hsbl(i,j)
#  endif
#  ifdef LMD_BKPP_NOT_YET
            avghbbl(i,j)=ad_hbbl(i,j)
#  endif
#  ifdef AVERAGES_FLUXES
            avgstf(i,j)=ad_stflx(i,j,itemp)
            avgswf(i,j)=ad_stflx(i,j,isalt)
#   ifdef BULK_FLUXES
            avglhf(i,j)=ad_lhflx(i,j)
            avgshf(i,j)=ad_shflx(i,j)
            avglrf(i,j)=ad_lrflx(i,j)
#    ifdef EMINUSP
            avgevap(i,j)=ad_evap(i,j)
!!          avgrain(i,j)=ad_rain(i,j)
#    endif
#   endif
#   ifdef SHORTWAVE
            avgsrf(i,j)=ad_srflx(i,j)
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
            avgsus(i,j)=ad_sustr(i,j)
            avgsvs(i,j)=ad_svstr(i,j)
            avgbus(i,j)=ad_bustr(i,j)
            avgbvs(i,j)=ad_bvstr(i,j)
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
              avgu3d(i,j,k)=ad_u(i,j,k,nstp)
              avgv3d(i,j,k)=ad_v(i,j,k,nstp)
              avgrho(i,j,k)=ad_rho(i,j,k)
#  ifdef AVERAGES_QUADRATIC
              avgUU(i,j,k)=ad_u(i,j,k,nstp)*ad_u(i,j,k,nstp)
              avgVV(i,j,k)=ad_v(i,j,k,nstp)*ad_v(i,j,k,nstp)
#  endif
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          DO j=Jstr,Jend
            DO i=Istr,Iend
              avgUV(i,j,k)=0.25_r8*(ad_u(i  ,j  ,k,nstp)+               &
     &                              ad_u(i+1,j  ,k,nstp))*              &
     &                             (ad_v(i  ,j  ,k,nstp)+               &
     &                              ad_v(i  ,j+1,k,nstp))
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
                avgt(i,j,k,itrc)=ad_t(i,j,k,nstp,itrc)
              END DO
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          IF (itrc.le.NAT) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  avgTT(i,j,k,itrc)=ad_t(i,j,k,nstp,itrc)*              &
     &                              ad_t(i,j,k,nstp,itrc)
                END DO
                DO i=Istr,Iend
                  avgUT(i,j,k,itrc)=0.5_r8*ad_u(i,j,k,nstp)*            &
     &                              (ad_t(i-1,j,k,nstp,itrc)+           &
     &                               ad_t(i  ,j,k,nstp,itrc))
                END DO
              END DO
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  avgVT(i,j,k,itrc)=0.5_r8*ad_v(i,j,k,nstp)*            &
     &                              (ad_t(i,j-1,k,nstp,itrc)+           &
     &                               ad_t(i,j  ,k,nstp,itrc))
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
              avgw3d(i,j,k)=ad_W(i,j,k)*pm(i,j)*pn(i,j)
#  ifdef AVERAGES_AKV
              avgAKv(i,j,k)=ad_Akv(i,j,k)
#  endif
#  ifdef AVERAGES_AKT
              avgAKt(i,j,k)=ad_Akt(i,j,k,itemp)
#  endif
#  ifdef AVERAGES_AKS
              avgAKs(i,j,k)=ad_Akt(i,j,k,isalt)
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
      ELSE IF ((iic(ng)-1).le.ntsAVG(ng)) THEN
!
!  Accumulate 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgzeta(i,j)=avgzeta(i,j)+ad_zeta_sol(i,j)
            avgu2d (i,j)=avgu2d (i,j)+ad_ubar_sol(i,j)
            avgv2d (i,j)=avgv2d (i,j)+ad_vbar_sol(i,j)
# ifdef AVERAGES_QUADRATIC
            avgZZ(i,j)=avgZZ(i,j)+ad_zeta_sol(i,j)*ad_zeta_sol(i,j)
            avgU2(i,j)=avgU2(i,j)+ad_ubar_sol(i,j)*ad_ubar_sol(i,j)
            avgV2(i,j)=avgV2(i,j)+ad_vbar_sol(i,j)*ad_vbar_sol(i,j)
# endif
# ifdef SOLVE3D
#  ifdef LMD_SKPP_NOT_YET
            avghsbl(i,j)=avghsbl(i,j)+ad_hsbl(i,j)
#  endif
#  ifdef LMD_BKPP_NOT_YET
            avghbbl(i,j)=avghbbl(i,j)+ad_hbbl(i,j)
#  endif
#  ifdef AVERAGES_FLUXES
            avgstf(i,j)=avgstf(i,j)+ad_stflx(i,j,itemp)
            avgswf(i,j)=avgswf(i,j)+ad_stflx(i,j,isalt)
#   ifdef BULK_FLUXES
            avglhf(i,j)=avglhf(i,j)+ad_lhflx(i,j)
            avgshf(i,j)=avgshf(i,j)+ad_shflx(i,j)
            avglrf(i,j)=avglrf(i,j)+ad_lrflx(i,j)
#    ifdef EMINUSP
            avgevap(i,j)=avgevap(i,j)+ad_evap(i,j)
!!          avgrain(i,j)=avgrain(i,j)+ad_rain(i,j)
#    endif
#   endif
#   ifdef SHORTWAVE
            avgsrf(i,j)=avgsrf(i,j)+ad_srflx(i,j)
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
            avgsus(i,j)=avgsus(i,j)+ad_sustr(i,j)
            avgsvs(i,j)=avgsvs(i,j)+ad_svstr(i,j)
            avgbus(i,j)=avgbus(i,j)+ad_bustr(i,j)
            avgbvs(i,j)=avgbvs(i,j)+ad_bvstr(i,j)
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
              avgu3d(i,j,k)=avgu3d(i,j,k)+ad_u(i,j,k,nstp)
              avgv3d(i,j,k)=avgv3d(i,j,k)+ad_v(i,j,k,nstp)
              avgrho(i,j,k)=avgrho(i,j,k)+ad_rho(i,j,k)
#  ifdef AVERAGES_QUADRATIC
              avgUU(i,j,k)=avgUU(i,j,k)+                                &
     &                     ad_u(i,j,k,nstp)*ad_u(i,j,k,nstp)
              avgVV(i,j,k)=avgVV(i,j,k)+                                &
     &                     ad_v(i,j,k,nstp)*ad_v(i,j,k,nstp)
#  endif
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          DO j=Jstr,Jend
            DO i=Istr,Iend
              avgUV(i,j,k)=avgUV(i,j,k)+                                &
     &                     0.25_r8*(ad_u(i  ,j  ,k,nstp)+               &
     &                              ad_u(i+1,j  ,k,nstp))*              &
     &                             (ad_v(i  ,j  ,k,nstp)+               &
     &                              ad_v(i  ,j+1,k,nstp))
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
                avgt(i,j,k,itrc)=avgt(i,j,k,itrc)+ad_t(i,j,k,nstp,itrc)
              END DO
            END DO
          END DO
#  ifdef AVERAGES_QUADRATIC
          IF (itrc.le.NAT) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  avgTT(i,j,k,itrc)=avgTT(i,j,k,itrc)+                  &
     &                              ad_t(i,j,k,nstp,itrc)*              &
     &                              ad_t(i,j,k,nstp,itrc)
                END DO
                DO i=Istr,Iend
                  avgUT(i,j,k,itrc)=avgUT(i,j,k,itrc)+                  &
     &                              0.5_r8*ad_u(i,j,k,nstp)*            &
     &                              (ad_t(i-1,j,k,nstp,itrc)+           &
     &                               ad_t(i  ,j,k,nstp,itrc))
                END DO
              END DO
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  avgVT(i,j,k,itrc)=avgVT(i,j,k,itrc)+                  &
     &                              0.5_r8*ad_v(i,j,k,nstp)*            &
     &                              (ad_t(i,j-1,k,nstp,itrc)+           &
     &                               ad_t(i,j  ,k,nstp,itrc))
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
              avgw3d(i,j,k)=avgw3d(i,j,k)+ad_W(i,j,k)*pm(i,j)*pn(i,j)
#  ifdef AVERAGES_AKV
              avgAKv(i,j,k)=avgAKv(i,j,k)+ad_Akv(i,j,k)
#  endif
#  ifdef AVERAGES_AKT
              avgAKt(i,j,k)=avgAKt(i,j,k)+ad_Akt(i,j,k,itemp)
#  endif
#  ifdef AVERAGES_AKS
              avgAKs(i,j,k)=avgAKs(i,j,k)+ad_Akt(i,j,k,isalt)
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
      IF ((iic(ng).lt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
        fac=1.0_r8
# else
        fac=1.0_r8/REAL(nAVG(ng),r8)
# endif
        IF (SOUTH_WEST_TEST) THEN
          AVGtime(ng)=AVGtime(ng)-REAL(nAVG(ng),r8)*dt(ng)
        END IF
!
!  Process 2D fields.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            avgzeta(i,j)=fac*avgzeta(i,j)
            avgu2d (i,j)=fac*avgu2d (i,j)
            avgv2d (i,j)=fac*avgv2d (i,j)
# ifdef AVERAGES_QUADRATIC
            avgZZ(i,j)=fac*avgZZ(i,j)
            avgU2(i,j)=fac*avgU2(i,j)
            avgV2(i,j)=fac*avgU2(i,j)
# endif
# ifdef SOLVE3D
#  ifdef LMD_SKPP_NOT_YET
            avghsbl(i,j)=fac*avghsbl(i,j)
#  endif
#  ifdef LMD_BKPP_NOT_YET
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
!!          avgrain(i,j)=fac*avgrain(i,j)
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
              avgu3d(i,j,k)=fac*avgu3d(i,j,k)
              avgv3d(i,j,k)=fac*avgv3d(i,j,k)
              avgrho(i,j,k)=fac*avgrho(i,j,k)
#  ifdef AVERAGES_QUADRATIC
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
                END DO
              END DO
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  avgVT(i,j,k,itrc)=fac*avgVT(i,j,k,itrc)
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
      END SUBROUTINE ad_set_avg_tile
#endif
      END MODULE ad_set_avg_mod
