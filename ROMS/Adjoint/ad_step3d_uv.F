#include "cppdefs.h"
      MODULE ad_step3d_uv_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_step3d_uv.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine time-steps  the  adjoint  horizontal  momentum      !
!  equations. The vertical viscosity terms are time-stepped using      !
!  an implicit algorithm.                                              !
!                                                                      !
!  BASIC STATE variables needed:  u, v, ru, rv, Akv,                   !
!                                 Hz, Huon, Hvom, z_r, z_w,            !
!                                 DU_avg1, DU_avg2, DV_avg1, DV_avg2,  !
!                                 Qsrc                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_step3d_uv

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_step3d_uv (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
# ifdef DIAGNOSTICS_UV
!!    USE mod_diags
# endif
      USE mod_forces
      USE mod_grid
      USE mod_mixing
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
      CALL wclock_on (ng, iADM, 34)
# endif
      CALL ad_step3d_uv_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        knew(ng), nrhs(ng), nstp(ng), nnew(ng),   &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef WET_DRY_NOT_YET
     &                        GRID(ng) % umask_wet,                     &
     &                        GRID(ng) % vmask_wet,                     &
# endif
     &                        GRID(ng) % om_v,                          &
     &                        GRID(ng) % on_u,                          &
     &                        GRID(ng) % pm,                            &
     &                        GRID(ng) % pn,                            &
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % ad_Hz,                         &
     &                        GRID(ng) % z_r,                           &
     &                        GRID(ng) % ad_z_r,                        &
     &                        GRID(ng) % z_w,                           &
     &                        GRID(ng) % ad_z_w,                        &
     &                        MIXING(ng) % Akv,                         &
     &                        MIXING(ng) % ad_Akv,                      &
     &                        COUPLING(ng) % DU_avg1,                   &
     &                        COUPLING(ng) % ad_DU_avg1,                &
     &                        COUPLING(ng) % DV_avg1,                   &
     &                        COUPLING(ng) % ad_DV_avg1,                &
     &                        COUPLING(ng) % DU_avg2,                   &
     &                        COUPLING(ng) % ad_DU_avg2,                &
     &                        COUPLING(ng) % DV_avg2,                   &
     &                        COUPLING(ng) % ad_DV_avg2,                &
     &                        OCEAN(ng) % ad_ru,                        &
     &                        OCEAN(ng) % ad_rv,                        &
# ifdef DIAGNOSTICS_UV
!!   &                        DIAGS(ng) % DiaU2wrk,                     &
!!   &                        DIAGS(ng) % DiaV2wrk,                     &
!!   &                        DIAGS(ng) % DiaU2int,                     &
!!   &                        DIAGS(ng) % DiaV2int,                     &
!!   &                        DIAGS(ng) % DiaU3wrk,                     &
!!   &                        DIAGS(ng) % DiaV3wrk,                     &
!!   &                        DIAGS(ng) % DiaRU,                        &
!!   &                        DIAGS(ng) % DiaRV,                        &
# endif
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % ad_u,                         &
     &                        OCEAN(ng) % v,                            &
     &                        OCEAN(ng) % ad_v,                         &
     &                        OCEAN(ng) % ad_ubar,                      &
     &                        OCEAN(ng) % ad_vbar,                      &
     &                        OCEAN(ng) % ad_ubar_sol,                  &
     &                        OCEAN(ng) % ad_vbar_sol,                  &
# ifdef NEARSHORE_MELLOR
     &                       OCEAN(ng) % ubar_stokes,                   &
     &                       OCEAN(ng) % tl_ubar_stokes,                &
     &                       OCEAN(ng) % vbar_stokes,                   &
     &                       OCEAN(ng) % tl_vbar_stokes,                &
     &                       OCEAN(ng) % u_stokes,                      &
     &                       OCEAN(ng) % tl_u_stokes,                   &
     &                       OCEAN(ng) % v_stokes,                      &
     &                       OCEAN(ng) % tl_v_stokes,                   &
# endif
     &                        GRID(ng) % Huon,                          &
     &                        GRID(ng) % ad_Huon,                       &
     &                        GRID(ng) % Hvom,                          &
     &                        GRID(ng) % ad_Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 34)
# endif
      RETURN
      END SUBROUTINE ad_step3d_uv
!
!***********************************************************************
      SUBROUTINE ad_step3d_uv_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              knew, nrhs, nstp, nnew,             &
# ifdef MASKING
     &                              umask, vmask,                       &
# endif
# ifdef WET_DRY_NOT_YET
     &                              umask_wet, vmask_wet,               &
# endif
     &                              om_v, on_u, pm, pn,                 &
     &                              Hz, ad_Hz,                          &
     &                              z_r, ad_z_r,                        &
     &                              z_w, ad_z_w,                        &
     &                              Akv, ad_Akv,                        &
     &                              DU_avg1, ad_DU_avg1,                &
     &                              DV_avg1, ad_DV_avg1,                &
     &                              DU_avg2, ad_DU_avg2,                &
     &                              DV_avg2, ad_DV_avg2,                &
     &                              ad_ru, ad_rv,                       &
# ifdef DIAGNOSTICS_UV
!!   &                              DiaU2wrk, DiaV2wrk,                 &
!!   &                              DiaU2int, DiaV2int,                 &
!!   &                              DiaU3wrk, DiaV3wrk,                 &
!!   &                              DiaRU,    DiaRV,                    &
# endif
     &                              u, ad_u,                            &
     &                              v, ad_v,                            &
     &                              ad_ubar, ad_vbar,                   &
     &                              ad_ubar_sol, ad_vbar_sol,           &
# ifdef NEARSHORE_MELLOR
     &                              ubar_stokes, ad_ubar_stokes,        &
     &                              vbar_stokes, ad_vbar_stokes,        &
     &                              u_stokes, ad_u_stokes,              &
     &                              v_stokes, ad_v_stokes,              &
# endif
     &                              Huon, ad_Huon,                      &
     &                              Hvom, ad_Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sources
!
      USE ad_exchange_2d_mod
      USE ad_exchange_3d_mod
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d, ad_mp_exchange3d
# endif
      USE ad_u3dbc_mod, ONLY : ad_u3dbc_tile
      USE ad_v3dbc_mod, ONLY : ad_v3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: knew, nrhs, nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(in) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DU_avg2(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg2(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: ubar_stokes(LBi:,LBj:)
      real(r8), intent(in) :: vbar_stokes(LBi:,LBj:)
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
#  endif
#  ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaU2wrk(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaV2wrk(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaU2int(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaV2int(LBi:,LBj:,:)
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
!!    real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_DU_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DU_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: ad_DV_avg2(LBi:,LBj:)
      real(r8), intent(inout) :: ad_ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: ad_rv(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: ad_ubar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: ad_vbar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: ad_u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Hvom(LBi:,LBj:,:)

      real(r8), intent(out) :: ad_ubar_sol(LBi:,LBj:)
      real(r8), intent(out) :: ad_vbar_sol(LBi:,LBj:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vbar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
#  ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaU2wrk(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaV2wrk(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaU2int(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaV2int(LBi:UBi,LBj:UBj,NDM2d)
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
!!    real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
!!    real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif
      real(r8), intent(inout) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: ad_rv(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: ad_ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_vbar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Hvom(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(out) :: ad_ubar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: ad_vbar_sol(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, idiag, is, j, k

      real(r8) :: cff, cff1, cff2
      real(r8) :: ad_cff, ad_cff1, ad_cff2
      real(r8) :: adfac, adfac1, adfac2

      real(r8), dimension(IminS:ImaxS) :: CF1
      real(r8), dimension(IminS:ImaxS) :: FC1
# ifdef NEARSHORE_MELLOR
      real(r8), dimension(IminS:ImaxS) :: CFs1
      real(r8), dimension(IminS:ImaxS) :: DCs1
# endif
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: AK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC1
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
# ifdef NEARSHORE_MELLOR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CFs
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DCs
# endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hzk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: oHz
# ifdef DIAGNOSTICS_UV
!!    real(r8), dimension(IminS:ImaxS,0:N(ng)) :: wrk
!!    real(r8), dimension(IminS:ImaxS,1:NDM2d-1) :: Dwrk
# endif
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_AK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_FC
# ifdef NEARSHORE_MELLOR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_CFs
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_DCs
# endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Hzk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_oHz
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8
      ad_cff1=0.0_r8
      DO k=0,N(ng)
        DO i=IminS,ImaxS
          ad_AK(i,k)=0.0_r8
          ad_BC(i,k)=0.0_r8
          ad_CF(i,k)=0.0_r8
          ad_DC(i,k)=0.0_r8
# ifdef NEARSHORE_MELLOR
          ad_CFs(i,k)=0.0_r8
          ad_DCs(i,k)=0.0_r8
# endif
          ad_FC(i,k)=0.0_r8
        END DO
      END DO
      DO k=1,N(ng)
        DO i=IminS,ImaxS
          ad_Hzk(i,k)=0.0_r8
          ad_oHz(i,k)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iTLM, 4,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_ubar(:,:,1), tl_vbar(:,:,1),               &
!>   &                    tl_ubar(:,:,2), tl_vbar(:,:,2))
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 4,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ubar(:,:,1), ad_vbar(:,:,1),            &
     &                       ad_ubar(:,:,2), ad_vbar(:,:,2))

!>    CALL mp_exchange3d (ng, tile, iTLM, 4,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_u(:,:,:,nnew), tl_v(:,:,:,nnew),           &
!>   &                    tl_Huon, tl_Hvom)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 4,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_u(:,:,:,nnew), ad_v(:,:,:,nnew),        &
     &                       ad_Huon, ad_Hvom)
!
# endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO k=1,2
!>        CALL exchange_v2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_vbar(:,:,k))
!>
          CALL ad_exchange_v2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_vbar(:,:,k))
!>        CALL exchange_u2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            tl_ubar(:,:,k))
!>
          CALL ad_exchange_u2d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               ad_ubar(:,:,k))
        END DO
!
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Hvom)
!>      CALL exchange_v3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_Hvom)
!>
        CALL ad_exchange_v3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_Hvom)

        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Huon)
!>      CALL exchange_u3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_Huon)
!>
        CALL ad_exchange_u3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_Huon)

!>      CALL exchange_v3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_v(:,:,:,nnew))
!>
        CALL ad_exchange_v3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_v(:,:,:,nnew))
!>      CALL exchange_u3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_u(:,:,:,nnew))
!>
        CALL ad_exchange_u3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_u(:,:,:,nnew))
      END IF
!
!-----------------------------------------------------------------------
!  Couple 2D and 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Couple adjoint velocity component in the ETA-direction.
!
      J_LOOP1 : DO j=JstrT,JendT
        IF (j.ge.Jstr) THEN
!
!  First, compute BASIC STATE quantities. This section was computed
!  three times in the original code. This can avoided by saving the
!  intermediate values in scratch arrays DC1, CF1, and FC1 (HGA).
!
          DO i=IstrT,IendT
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
# ifdef NEARSHORE_MELLOR
            CFs(i,0)=0.0_r8
# endif
            FC(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IstrT,IendT
              cff=0.5_r8*om_v(i,j)
              DC(i,k)=cff*(Hz(i,j,k)+Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+                                          &
     &                DC(i,k)*v(i,j,k,nnew)
# ifdef NEARSHORE_MELLOR
              CFs(i,0)=CFs(i,0)+                                        &
     &                 DC(i,k)*v_stokes(i,j,k)
# endif
            END DO
          END DO
          DO i=IstrT,IendT
            DC1(i,0)=DC(i,0)                            ! intermediate
# ifdef NEARSHORE_MELLOR
            cff2=DC(i,0)*vbar_stokes(i,j)
# endif
            DC(i,0)=1.0_r8/DC(i,0)                      ! recursive
            CF1(i)=CF(i,0)                              ! intermediate
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg1(i,j))      ! recursive
# ifdef NEARSHORE_MELLOR
            CFs1(i)=CFs(i,0)
            CFs(i,0)=DC(i,0)*(CFs(i,0)-cff2)            ! recursive
# endif
          END DO
          DO k=N(ng),1,-1
            DO i=IstrT,IendT
!>            Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
# ifdef NEARSHORE_MELLOR
!>            Hvom(i,j,k)=Hvom(i,j,k)+0.5_r8*v_stokes(i,j,k)*DC(i,k)
# endif
!>            FC(i,0)=FC(i,0)+Hvom(i,j,k)
!>
              FC(i,0)=FC(i,0)+                                          &
     &                0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
# ifdef NEARSHORE_MELLOR
              FC(i,0)=FC(i,0)+                                          &
     &                0.5_r8*v_stokes(i,j,k)*DC(i,k)
# endif
            END DO
          END DO
          DO i=IstrT,IendT
            FC1(i)=FC(i,0)                              ! intermediate
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))      ! recursive
          END DO
!
!  Compute correct mass flux, Hz*v/m.
!
          DO k=1,N(ng)
            DO i=IstrT,IendT
!>            tl_Hvom(i,j,k)=tl_Hvom(i,j,k)-                            &
!>   &                       tl_DC(i,k)*FC(i,0)-                        &
!>   &                       DC(i,k)*tl_FC(i,0)
!>
              ad_DC(i,k)=ad_DC(i,k)-ad_Hvom(i,j,k)*FC(i,0)
              ad_FC(i,0)=ad_FC(i,0)-ad_Hvom(i,j,k)*DC(i,k)
            END DO
          END DO
          DO i=IstrT,IendT
!>          tl_FC(i,0)=tl_DC(i,0)*(FC1(i)-DV_avg2(i,j))+                &
!>   &                 DC(i,0)*(tl_FC(i,0)-tl_DV_avg2(i,j))
!>
            adfac=DC(i,0)*ad_FC(i,0)
            ad_DC(i,0)=ad_DC(i,0)+ad_FC(i,0)*(FC1(i)-DV_avg2(i,j))
            ad_DV_avg2(i,j)=ad_DV_avg2(i,j)-adfac
            ad_FC(i,0)=adfac
          END DO
!>        DO k=N(ng),1,-1
!>
          DO k=1,N(ng)
            DO i=IstrT,IendT
# ifdef DIAGNOSTICS_UV
!!            DiaV3wrk(i,j,k,M3rate)=v(i,j,k,nnew)-                     &
!!   &                               DiaV3wrk(i,j,k,M3rate)
# endif
!>            tl_FC(i,0)=tl_FC(i,0)+tl_Hvom(i,j,k)
!>
              ad_Hvom(i,j,k)=ad_Hvom(i,j,k)+ad_FC(i,0)
# ifdef NEARSHORE_MELLOR
!>            tl_Hvom(i,j,k)=tl_Hvom(i,j,k)+                            &
!>   &                       0.5_r8*(tl_v_stokes(i,j,k)*DC(i,k)+        &
!>   &                               v_stokes(i,j,k)*tl_DC(i,k))
!>
              adfac=0.5_r8*ad_Hvom(i,j,k)
              ad_DC(i,k)=ad_DC(i,k)+v_stokes(i,j,k)*adfac
              ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)+DC(i,k)*adfac
# endif
!>            tl_Hvom(i,j,k)=0.5_r8*(tl_Hvom(i,j,k)+                    &
!>   &                               tl_v(i,j,k,nnew)*DC(i,k)+          &
!>   &                               v(i,j,k,nnew)*tl_DC(i,k))
!>
              adfac=0.5_r8*ad_Hvom(i,j,k)
              ad_DC(i,k)=ad_DC(i,k)+v(i,j,k,nnew)*adfac
              ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)+DC(i,k)*adfac
              ad_Hvom(i,j,k)=adfac
            END DO
          END DO
!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, vbar=DV_avg1/(D*om_v).  Recall that,
!  D=CF(:,0). The replacement is avoided if the boundary edge is
!  periodic. The J-loop is pipelined, so we need to do a special
!  test for the southern and northern domain edges.
!
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (j.eq.Mm(ng)+1) THEN                  ! northern boundary
              DO k=1,N(ng)                           ! J-loop pipelined
                DO i=Istr,Iend
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>                tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*                &
!>   &                               vmask_wet(i,j)
!>
                  ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)*                &
     &                               vmask_wet(i,j)
#  endif
#  ifdef MASKING
!>                tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*                &
!>   &                               vmask(i,j)
!>
                  ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)*                &
     &                               vmask(i,j)
#  endif
!>                tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)-                &
!>   &                               tl_CFs(i,0)
!>
                  ad_CFs(i,0)=ad_CFs(i,0)-ad_v_stokes(i,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>                tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*                    &
!>   &                               vmask_wet(i,j)
!>
                  ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*                    &
     &                             vmask_wet(i,j)
# endif
# ifdef MASKING
!>                tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*                    &
!>   &                             vmask(i,j)
!>
                  ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*                    &
     &                             vmask(i,j)
# endif
!>                tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-                    &
!>   &                             tl_CF(i,0)
!>
                  ad_CF(i,0)=ad_CF(i,0)-ad_v(i,j,k,nnew)
                END DO
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (j.eq.1) THEN                         ! southern boundary
              DO k=1,N(ng)                           ! J-loop pipelined
                DO i=Istr,Iend
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>                tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*                &
!>   &                               vmask_wet(i,j)
!>
                  ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)*                &
     &                               vmask_wet(i,j)
#  endif
#  ifdef MASKING
!>                tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*                &
!>   &                               vmask(i,j)
!>
                  ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)*                &
     &                               vmask(i,j)
#  endif
!>                tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)-                &
!>   &                               tl_CFs(i,0)
!>
                  ad_CFs(i,0)=ad_CFs(i,0)-ad_v_stokes(i,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>                tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*                     &
!>   &                             vmask_wet(i,j)
!>
                  ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*                     &
     &                             vmask_wet(i,j)
# endif
# ifdef MASKING
!>                tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*                     &
!>   &                             vmask(i,j)
!>
                  ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*                     &
     &                             vmask(i,j)
# endif
!>                tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-                     &
!>   &                             tl_CF(i,0)
!>
                  ad_CF(i,0)=ad_CF(i,0)-ad_v(i,j,k,nnew)
                END DO
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO k=1,N(ng)
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>              tl_v_stokes(Iend+1,j,k)=tl_v_stokes(Iend+1,j,k)*        &
!>   &                                  vmask_wet(Iend+1,j)
!>
                ad_v_stokes(Iend+1,j,k)=ad_v_stokes(Iend+1,j,k)*        &
     &                                  vmask_wet(Iend+1,j)
#  endif
#  ifdef MASKING
!>              tl_v_stokes(Iend+1,j,k)=tl_v_stokes(Iend+1,j,k)*        &
!>   &                                  vmask(Iend+1,j)
!>
                ad_v_stokes(Iend+1,j,k)=ad_v_stokes(Iend+1,j,k)*        &
     &                                  vmask(Iend+1,j)
#  endif
!>              tl_v_stokes(Iend+1,j,k)=tl_v_stokes(Iend+1,j,k)-        &
!>   &                                  tl_CFs(Iend+1,0)
!>
                ad_CFs(Iend+1,0)=ad_CFs(Iend+1,0)-                      &
     &                           ad_v_stokes(Iend+1,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>              tl_v(Iend+1,j,k,nnew)=tl_v(Iend+1,j,k,nnew)*            &
!>   &                                vmask_wet(Iend+1,j)
!>
                ad_v(Iend+1,j,k,nnew)=ad_v(Iend+1,j,k,nnew)*            &
     &                                vmask_wet(Iend+1,j)
# endif
# ifdef MASKING
!>              tl_v(Iend+1,j,k,nnew)=tl_v(Iend+1,j,k,nnew)*            &
!>   &                                vmask(Iend+1,j)
!>
                ad_v(Iend+1,j,k,nnew)=ad_v(Iend+1,j,k,nnew)*            &
     &                                vmask(Iend+1,j)
# endif
!>              tl_v(Iend+1,j,k,nnew)=tl_v(Iend+1,j,k,nnew)-            &
!>   &                                    tl_CF(Iend+1,0)
!>
                ad_CF(Iend+1,0)=ad_CF(Iend+1,0)-                        &
     &                          ad_v(Iend+1,j,k,nnew)
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO k=1,N(ng)
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>              tl_v_stokes(Istr-1,j,k)=tl_v_stokes(Istr-1,j,k)*        &
!>   &                                  vmask_wet(Istr-1,j)
!>
                ad_v_stokes(Istr-1,j,k)=ad_v_stokes(Istr-1,j,k)*        &
     &                                  vmask_wet(Istr-1,j)
#  endif
#  ifdef MASKING
!>              tl_v_stokes(Istr-1,j,k)=tl_v_stokes(Istr-1,j,k)*        &
!>   &                                  vmask(Istr-1,j)
!>
                ad_v_stokes(Istr-1,j,k)=ad_v_stokes(Istr-1,j,k)*        &
     &                                  vmask(Istr-1,j)
#  endif
!>              tl_v_stokes(Istr-1,j,k)=tl_v_stokes(Istr-1,j,k)-        &
!>   &                                  tl_CFs(Istr-1,0)
!>
                ad_CFs(Istr-1,0)=ad_CFs(Istr-1,0)-                      &
     &                           ad_v_stokes(Istr-1,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>              tl_v(Istr-1,j,k,nnew)=tl_v(Istr-1,j,k,nnew)*            &
!>   &                                vmask_wet(Istr-1,j)
!>
                ad_v(Istr-1,j,k,nnew)=ad_v(Istr-1,j,k,nnew)*            &
     &                                vmask_wet(Istr-1,j)
# endif
# ifdef MASKING
!>              tl_v(Istr-1,j,k,nnew)=tl_v(Istr-1,j,k,nnew)*            &
!>   &                                vmask(Istr-1,j)
!>
                ad_v(Istr-1,j,k,nnew)=ad_v(Istr-1,j,k,nnew)*            &
     &                                vmask(Istr-1,j)
# endif
!>              tl_v(Istr-1,j,k,nnew)=tl_v(Istr-1,j,k,nnew)-            &
!>   &                                tl_CF(Istr-1,0)
!>
                ad_CF(Istr-1,0)=ad_CF(Istr-1,0)-                        &
     &                          ad_v(Istr-1,j,k,nnew)
              END DO
            END IF
          END IF

# ifdef DIAGNOSTICS_UV
!!
!! Convert the units of the fast-time integrated change in mass flux
!! diagnostic values to change in velocity (m s-1).
!!
!!        DO idiag=1,NDM2d-1
!!          DO i=IstrT,IendT
#  ifdef MASKING
!!            DiaV2wrk(i,j,idiag)=DiaV2wrk(i,j,idiag)*vmask(i,j)
#  endif
!!            DiaV2wrk(i,j,idiag)=DC(i,0)*DiaV2wrk(i,j,idiag)
!!          END DO
!!        END DO
# endif
!
!  Save adjoint solution for time-step iic(ng)-1 as the sum of time
!  records 1 and 2.
!
          DO i=IstrT,IendT
            ad_vbar_sol(i,j)=ad_vbar(i,j,1)+ad_vbar(i,j,2)
          END DO
!
!  Compute adjoint thicknesses of V-boxes DC(i,1:N), total depth of the
!  water column DC(i,0), and incorrect vertical mean CF(i,0).  Notice
!  that barotropic component is replaced with its fast-time averaged
!  values.
!
          DO i=IstrT,IendT
# ifdef DIAGNOSTICS_UV
!!          DiaV2int(i,j,M2rate)=vbar(i,j,1)               ! HGA check
!!          DiaV2wrk(i,j,M2rate)=vbar(i,j,1)-DiaV2int(i,j,M2rate)
!!          DiaV2int(i,j,M2rate)=vbar(i,j,1)*DC1(i,0)
!!          DiaV2wrk(i,j,M2rate)=vbar(i,j,1)-                           &
!!   &                           DiaV2int(i,j,M2rate)*DC(i,0)
# endif
!>          tl_vbar(i,j,2)=tl_vbar(i,j,1)
!>
            ad_vbar(i,j,1)=ad_vbar(i,j,1)+ad_vbar(i,j,2)
            ad_vbar(i,j,2)=0.0_r8
!>          tl_vbar(i,j,1)=tl_DC(i,0)*DV_avg1(i,j)+                     &
!>   &                     DC(i,0)*tl_DV_avg1(i,j)
!>
            ad_DC(i,0)=ad_DC(i,0)+ad_vbar(i,j,1)*DV_avg1(i,j)
            ad_DV_avg1(i,j)=ad_DV_avg1(i,j)+ad_vbar(i,j,1)*DC(i,0)
            ad_vbar(i,j,1)=0.0_r8
# ifdef NEARSHORE_MELLOR
!>          tl_CFs(i,0)=tl_DC(i,0)*(CFs1(i)-cff2)+                      &
!>   &                  DC(i,0)*(tl_CFs(i,0)-tl_cff2)
!>
            adfac=DC(i,0)*ad_CFs(i,0)
            ad_DC(i,0)=ad_DC(i,0)+(CFs1(i)-cff2)*ad_CFs(i,0)
            ad_cff2=ad_cff2-adfac
            ad_CFs(i,0)=ad_CFs(i,0)+adfac
            ad_CFs(i,0)=0.0_r8
# endif
!>          tl_CF(i,0)=tl_DC(i,0)*(CF1(i)-DV_avg1(i,j))+                &
!>   &                 DC(i,0)*(tl_CF(i,0)-tl_DV_avg1(i,j))
!>
            adfac=DC(i,0)*ad_CF(i,0)
            ad_DC(i,0)=ad_DC(i,0)+ad_CF(i,0)*(CF1(i)-DV_avg1(i,j))
            ad_DV_avg1(i,j)=ad_DV_avg1(i,j)-adfac
            ad_CF(i,0)=adfac
!>          tl_DC(i,0)=-tl_DC(i,0/(DC1(i,0)*DC1(i,0))
!>
            ad_DC(i,0)=-ad_DC(i,0)/(DC1(i,0)*DC1(i,0))
# ifdef NEARSHORE_MELLOR
!>          tl_cff2=tl_DC(i,0)*vbar_stokes(i,j)+                        &
!>   &              DC(i,0)*tl_vbar_stokes(i,j))
!>
            ad_vbar_stokes(i,j)=ad_vbar_stokes(i,j)+DC(i,0)*ad_cff2
            ad_DC(i,0)=ad_DC(i,0)+vbar_stokes(i,j)*ad_cff2
            ad_cff2=0.0_r8
# endif
          END DO
          DO i=IstrT,IendT
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
# ifdef NEARSHORE_MELLOR
            CFs(i,0)=0.0_r8
# endif
            FC(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IstrT,IendT
              cff=0.5_r8*om_v(i,j)
              DC(i,k)=cff*(Hz(i,j,k)+Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+                                          &
     &                DC(i,k)*v(i,j,k,nnew)
# ifdef NEARSHORE_MELLOR
              CFs(i,0)=CFs(i,0)+DC(i,k)*v_stokes(i,j,k)
!>            tl_CFs(i,0)=tl_CFs(i,0)+                                  &
!>   &                    tl_DC(i,k)*v_stokes(i,j,k)+                   &
!>   &                    DC(i,k)*tl_v_stokes(i,j,k)
!>
              ad_DC(i,k)=ad_DC(i,k)+v_stokes(i,j,k)*ad_CFs(i,0)
              ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)+DC(i,k)*ad_CFs(i,0)
# endif
!>            tl_CF(i,0)=tl_CF(i,0)+                                    &
!>   &                   tl_DC(i,k)*v(i,j,k,nnew)+                      &
!>   &                   DC(i,k)*tl_v(i,j,k,nnew)
!>
              ad_DC(i,k)=ad_DC(i,k)+ad_CF(i,0)*v(i,j,k,nnew)
              ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)+ad_CF(i,0)*DC(i,k)
!>            tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
!>
              ad_DC(i,k)=ad_DC(i,k)+ad_DC(i,0)
!>            tl_DC(i,k)=cff*(tl_Hz(i,j,k)+tl_Hz(i,j-1,k))
!>
              adfac=cff*ad_DC(i,k)
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
          DO i=IstrT,IendT
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
# ifdef NEARSHORE_MELLOR
!>          tl_CFs(i,0)=0.0_r8
!>
            ad_CFs(i,0)=0.0_r8
# endif
!>          tl_CF(i,0)=0.0_r8
!>
            ad_CF(i,0)=0.0_r8
!>          tl_DC(i,0)=0.0_r8
!>
            ad_DC(i,0)=0.0_r8
          END DO
        END IF
!
!  Couple adjoint velocity component in the XI-direction.  First,
!  compute BASIC STATE quantities. This section was computed three
!  times in the original code. This can avoided by saving the
!  intermediate values in scratch arrays DC1, CF1, and FC1 (HGA).
!
        DO i=IstrP,IendT
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
# ifdef NEARSHORE_MELLOR
          CFs(i,0)=0.0_r8
# endif
          FC(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrP,IendT
            cff=0.5_r8*on_u(i,j)
            DC(i,k)=cff*(Hz(i,j,k)+Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+                                            &
     &              DC(i,k)*u(i,j,k,nnew)
# ifdef NEARSHORE_MELLOR
            CFs(i,0)=CFs(i,0)+                                          &
     &               DC(i,k)*u_stokes(i,j,k)
# endif
          END DO
        END DO
        DO i=IstrP,IendT
          DC1(i,0)=DC(i,0)                              ! intermediate
# ifdef NEARSHORE_MELLOR
          cff2=DC(i,0)*ubar_stokes(i,j)
# endif
          DC(i,0)=1.0_r8/DC(i,0)                        ! recursive
          CF1(i)=CF(i,0)                                ! intermediate
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j))        ! recursive
# ifdef NEARSHORE_MELLOR
          CFs1(i)=CFs(i,0)                              ! intermediate
          CFs(i,0)=DC(i,0)*(CFs(i,0)-cff2)              ! recursive
# endif
        END DO
        DO k=N(ng),1,-1
          DO i=IstrP,IendT
!>          Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))
# ifdef NEARSHORE_MELLOR
!>          Huon(i,j,k)=Huon(i,j,k)+0.5_r8*u_stokes(i,j,k)*DC(i,k)
# endif
!>          FC(i,0)=FC(i,0)+Huon(i,j,k)
!>
            FC(i,0)=FC(i,0)+                                            &
     &              0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))
# ifdef NEARSHORE_MELLOR
            FC(i,0)=FC(i,0)+                                            &
     &              0.5_r8*u_stokes(i,j,k)*DC(i,k)
# endif
          END DO
        END DO
        DO i=IstrP,IendT
          FC1(i)=FC(i,0)                                ! intermediate
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))        ! recursive
        END DO
!
!  Compute correct adjoint mass flux, Hz*u/n.
!
        DO k=1,N(ng)
          DO i=IstrP,IendT
!>          tl_Huon(i,j,k)=tl_Huon(i,j,k)-                              &
!>   &                     tl_DC(i,k)*FC(i,0)-                          &
!>   &                     DC(i,k)*tl_FC(i,0)
!>
            ad_DC(i,k)=ad_DC(i,k)-ad_Huon(i,j,k)*FC(i,0)
            ad_FC(i,0)=ad_FC(i,0)-ad_Huon(i,j,k)*DC(i,k)
          END DO
        END DO
        DO i=IstrP,IendT
!>        tl_FC(i,0)=tl_DC(i,0)*(FC1(i)-DU_avg2(i,j))+                  &
!>   &               DC(i,0)*(tl_FC(i,0)-tl_DU_avg2(i,j))
!>
          adfac=DC(i,0)*ad_FC(i,0)
          ad_DC(i,0)=ad_DC(i,0)+ad_FC(i,0)*(FC1(i)-DU_avg2(i,j))
          ad_DU_avg2(i,j)=ad_DU_avg2(i,j)-adfac
          ad_FC(i,0)=adfac
        END DO
!>      DO k=N(ng),1,-1
!>
        DO k=1,N(ng)
          DO i=IstrP,IendT
# ifdef DIAGNOSTICS_UV
!!          DiaU3wrk(i,j,k,M3rate)=u(i,j,k,nnew)-DiaU3wrk(i,j,k,M3rate)
# endif
!>          tl_FC(i,0)=tl_FC(i,0)+tl_Huon(i,j,k)
!>
            ad_Huon(i,j,k)=ad_Huon(i,j,k)+ad_FC(i,0)
# ifdef NEARSHORE_MELLOR
!>          tl_Huon(i,j,k)=tl_Huon(i,j,k)+                              &
!>   &                     0.5_r8*(tl_u_stokes(i,j,k)*DC(i,k)+          &
!>   &                             u_stokes(i,j,k)*tl_DC(i,k))
!>
            adfac=0.5_r8*ad_Huon(i,j,k)
            ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)+DC(i,k)*adfac
            ad_DC(i,k)=ad_DC(i,k)+u_stokes(i,j,k)*adfac
# endif
!>          tl_Huon(i,j,k)=0.5_r8*(tl_Huon(i,j,k)+                      &
!>   &                             tl_u(i,j,k,nnew)*DC(i,k)+            &
!>   &                             u(i,j,k,nnew)*tl_DC(i,k))
!>
            adfac=0.5_r8*ad_Huon(i,j,k)
            ad_DC(i,k)=ad_DC(i,k)+u(i,j,k,nnew)*adfac
            ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)+DC(i,k)*adfac
            ad_Huon(i,j,k)=adfac
          END DO
        END DO
!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, ubar=DU_avg1/(D*on_u). Recall that,
!  D=CF(:,0). The replacement is avoided if the boundary edge is
!  periodic. The J-loop is pipelined, so we need to do a special
!  test for the southern and northern domain edges.
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (j.eq.Mm(ng)+1) THEN                    ! northern boundary
            DO k=1,N(ng)                             ! J-loop piplined
              DO i=IstrU,Iend
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>              tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*                  &
!>   &                             umask_wet(i,j)
!>
                ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)*                  &
     &                             umask_wet(i,j)
#  endif
#  ifdef MASKING
!>              tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*                  &
!>   &                             umask(i,j)
!>
                ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)*                  &
     &                             umask(i,j)
#  endif
!>              tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)-                  &
!>   &                             tl_CFs(i,0)
!>
                ad_CFs(i,0)=ad_CFs(i,0)-ad_u_stokes(i,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>              tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*                      &
!>   &                           umask_wet(i,j)
!>
                ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*                      &
     &                           umask_wet(i,j)
# endif
# ifdef MASKING
!>              tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*                      &
!>   &                           umask(i,j)
!>
                ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*                      &
     &                           umask(i,j)
# endif
!>              tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_CF(i,0)
!>
                ad_CF(i,0)=ad_CF(i,0)-ad_u(i,j,k,nnew)
              END DO
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (j.eq.0) THEN                           ! southern boundary
            DO k=1,N(ng)                             ! J-loop pipelined
              DO i=IstrU,Iend
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>              tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*                  &
!>   &                             umask_wet(i,j)
!>
                ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)*                  &
     &                             umask_wet(i,j)
#  endif
#  ifdef MASKING
!>              tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*                  &
!>   &                             umask(i,j)
!>
                ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)*                  &
     &                             umask(i,j)
#  endif
!>              tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)-                  &
!>   &                             tl_CFs(i,0)
!>
                ad_CFs(i,0)=ad_CFs(i,0)-ad_u_stokes(i,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>              tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*                      &
!>   &                           umask_wet(i,j)
!>
                ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*                      &
     &                           umask_wet(i,j)
# endif
# ifdef MASKING
!>              tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*                      &
!>   &                           umask(i,j)
!>
                ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*                      &
     &                           umask(i,j)
# endif
!>              tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_CF(i,0)
!>
                ad_CF(i,0)=ad_CF(i,0)-ad_u(i,j,k,nnew)
              END DO
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=1,N(ng)
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>            tl_u_stokes(Iend+1,j,k)=tl_u_stokes(Iend+1,j,k)*          &
!>   &                                umask_wet(Iend+1,j)
!>
              ad_u_stokes(Iend+1,j,k)=ad_u_stokes(Iend+1,j,k)*          &
     &                                umask_wet(Iend+1,j)
#  endif
#  ifdef MASKING
!>            tl_u_stokes(Iend+1,j,k)=tl_u_stokes(Iend+1,j,k)*          &
!>   &                                umask(Iend+1,j)
!>
              ad_u_stokes(Iend+1,j,k)=ad_u_stokes(Iend+1,j,k)*          &
     &                                umask(Iend+1,j)
#  endif
!>            tl_u_stokes(Iend+1,j,k)=tl_u_stokes(Iend+1,j,k)-          &
!>   &                                tl_CFs(Iend+1,0)
!>
              ad_CFs(Iend+1,0)=ad_CFs(Iend+1,0)-                        &
     &                         ad_u_stokes(Iend+1,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>            tl_u(Iend+1,j,k,nnew)=tl_u(Iend+1,j,k,nnew)*              &
!>   &                              umask_wet(Iend+1,j)
!>
              ad_u(Iend+1,j,k,nnew)=ad_u(Iend+1,j,k,nnew)*              &
     &                              umask_wet(Iend+1,j)
# endif
# ifdef MASKING
!>            tl_u(Iend+1,j,k,nnew)=tl_u(Iend+1,j,k,nnew)*              &
!>   &                              umask(Iend+1,j)
!>
              ad_u(Iend+1,j,k,nnew)=ad_u(Iend+1,j,k,nnew)*              &
     &                              umask(Iend+1,j)
# endif
!>            tl_u(Iend+1,j,k,nnew)=tl_u(Iend+1,j,k,nnew)-              &
!>   &                              tl_CF(Iend+1,0)
!>
              ad_CF(Iend+1,0)=ad_CF(Iend+1,0)-                          &
     &                        ad_u(Iend+1,j,k,nnew)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=1,N(ng)
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>            tl_u_stokes(Istr,j,k)=tl_u_stokes(Istr,j,k)*              &
!>   &                              umask_wet(Istr,j)
!>
              ad_u_stokes(Istr,j,k)=ad_u_stokes(Istr,j,k)*              &
     &                              umask_wet(Istr,j)
#  endif
#  ifdef MASKING
!>            tl_u_stokes(Istr,j,k)=tl_u_stokes(Istr,j,k)*              &
!>   &                              umask(Istr,j)
!>
              ad_u_stokes(Istr,j,k)=ad_u_stokes(Istr,j,k)*              &
     &                              umask(Istr,j)
#  endif
!>            tl_u_stokes(Istr,j,k)=tl_u_stokes(Istr,j,k)-              &
!>   &                              tl_CFs(Istr,0)
!>
              ad_CFs(Istr,0)=ad_CFs(Istr,0)-                            &
     &                       ad_u_stokes(Istr,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>            tl_u(Istr,j,k,nnew)=tl_u(Istr,j,k,nnew)*                  &
!>   &                            umask_wet(Istr,j)
!>
              ad_u(Istr,j,k,nnew)=ad_u(Istr,j,k,nnew)*                  &
     &                            umask_wet(Istr,j)
# endif
# ifdef MASKING
!>            tl_u(Istr,j,k,nnew)=tl_u(Istr,j,k,nnew)*                  &
!>   &                            umask(Istr,j)
!>
              ad_u(Istr,j,k,nnew)=ad_u(Istr,j,k,nnew)*                  &
     &                            umask(Istr,j)
# endif
!>            tl_u(Istr,j,k,nnew)=tl_u(Istr,j,k,nnew)-                  &
!>   &                            tl_CF(Istr,0)
!>
              ad_CF(Istr,0)=ad_CF(Istr,0)-                              &
     &                      ad_u(Istr,j,k,nnew)
            END DO
          END IF
        END IF

# ifdef DIAGNOSTICS_UV
!!
!! Convert the units of the fast-time integrated change in mass flux
!! diagnostic values to change in velocity (m s-1).
!!
!!      DO idiag=1,NDM2d-1
!!        DO i=IstrP,IendT
#  ifdef MASKING
!!          DiaU2wrk(i,j,idiag)=DiaU2wrk(i,j,idiag)*umask(i,j)
#  endif
!!          DiaU2wrk(i,j,idiag)=DC(i,0)*DiaU2wrk(i,j,idiag)
!!        END DO
!!      END DO
# endif
!
!  Save adjoint solution for time-step iic(ng)-1 as the sum of time
!  records 1 and 2.
!
        DO i=IstrP,IendT
          ad_ubar_sol(i,j)=ad_ubar(i,j,1)+ad_ubar(i,j,2)
        END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
        DO i=IstrP,IendT
# ifdef DIAGNOSTICS_UV
!!        DiaU2int(i,j,M2rate)=ubar(i,j,1)*DC1(i,0)
!!        DiaU2wrk(i,j,M2rate)=ubar(i,j,1)-DiaU2int(i,j,M2rate)*DC(i,0)
# endif
!>        tl_ubar(i,j,2)=tl_ubar(i,j,1)
!>
          ad_ubar(i,j,1)=ad_ubar(i,j,1)+ad_ubar(i,j,2)
          ad_ubar(i,j,2)=0.0_r8
!>        tl_ubar(i,j,1)=tl_DC(i,0)*DU_avg1(i,j)+                       &
!>   &                   DC(i,0)*tl_DU_avg1(i,j)
!>
          ad_DC(i,0)=ad_DC(i,0)+ad_ubar(i,j,1)*DU_avg1(i,j)
          ad_DU_avg1(i,j)=ad_DU_avg1(i,j)+ad_ubar(i,j,1)*DC(i,0)
          ad_ubar(i,j,1)=0.0_r8
# ifdef NEARSHORE_MELLOR
!>        tl_CFs(i,0)=tl_DC(i,0)*(CFs1(i)-cff2)+                        &
!>   &                DC(i,0)*(tl_CFs(i,0)-tl_cff2)
!>
          adfac=DC(i,0)*ad_CFs(i,0)
          ad_DC(i,0)=ad_DC(i,0)+(CFs1(i)-cff2)*ad_CFs(i,0)
          ad_CFs(i,0)=ad_CFs(i,0)+adfac
          ad_cff2=tl_cff2-adfac
          ad_CFs(i,0)=0.0_r8
# endif
!>        tl_CF(i,0)=tl_DC(i,0)*(CF1(i)-DU_avg1(i,j))+                  &
!>   &               DC(i,0)*(tl_CF(i,0)-tl_DU_avg1(i,j))
!>
          adfac=DC(i,0)*ad_CF(i,0)
          ad_DC(i,0)=ad_DC(i,0)+ad_CF(i,0)*(CF1(i)-DU_avg1(i,j))
          ad_DU_avg1(i,j)=ad_DU_avg1(i,j)-adfac
          ad_CF(i,0)=adfac
!>        tl_DC(i,0)=-tl_DC(i,0)/(DC1(i,0)*DC1(i,0))
!>
          ad_DC(i,0)=-ad_DC(i,0)/(DC1(i,0)*DC1(i,0))
# ifdef NEARSHORE_MELLOR
!>        tl_cff2=tl_DC(i,0)*ubar_stokes(i,j)+                          &
!>   &            DC(i,0)*tl_ubar_stokes(i,j)
!>
          ad_DC(i,0)=ad_DC(i,0)+ubar_stokes(i,j)*ad_cff2
          ad_ubar_stokes(i,j)=ad_ubar_stokes(i,j)+DC(i,0)*ad_cff2
          ad_cff2=0.0_r8
# endif
        END DO
        DO i=IstrP,IendT
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
# ifdef NEARSHORE_MELLOR
          CFs(i,0)=0.0_r8
# endif
          FC(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrP,IendT
            cff=0.5_r8*on_u(i,j)
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+                                            &
     &              DC(i,k)*u(i,j,k,nnew)
# ifdef NEARSHORE_MELLOR
            CFs(i,0)=CFs(i,0)+                                          &
     &               DC(i,k)*u_stokes(i,j,k)
!>          tl_CFs(i,0)=tl_CFs(i,0)+                                    &
!>   &                  tl_DC(i,k)*u_stokes(i,j,k)+                     &
!>   &                  DC(i,k)*tl_u_stokes(i,j,k)
!>
            ad_DC(i,k)=ad_DC(i,k)+u_stokes(i,j,k)*ad_CFs(i,0)
            ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)+DC(i,k)*ad_CFs(i,0)
# endif
!>          tl_CF(i,0)=tl_CF(i,0)+                                      &
!>   &                 tl_DC(i,k)*u(i,j,k,nnew)+                        &
!>   &                 DC(i,k)*tl_u(i,j,k,nnew)
!>
            ad_DC(i,k)=ad_DC(i,k)+ad_CF(i,0)*u(i,j,k,nnew)
            ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)+ad_CF(i,0)*DC(i,k)
!>          tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
!>
            ad_DC(i,k)=ad_DC(i,k)+ad_DC(i,0)
!>          tl_DC(i,k)=cff*(tl_Hz(i,j,k)+tl_Hz(i-1,j,k))
!>
            adfac=cff*ad_DC(i,k)
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac
            ad_DC(i,k)=0.0_r8
          END DO
        END DO
        DO i=IstrP,IendT
!>        tl_FC(i,0)=0.0_r8
!>
          ad_FC(i,0)=0.0_r8
# ifdef NEARSHORE_MELLOR
!>        tl_CFs(i,0)=0.0_r8
!>
          ad_CFs(i,0)=0.0_r8
# endif
!>        tl_CF(i,0)=0.0_r8
!>
          ad_CF(i,0)=0.0_r8
!>        tl_DC(i,0)=0.0_r8
!>
          ad_DC(i,0)=0.0_r8
        END DO
      END DO J_LOOP1
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
          IF (((IstrR.le.i).and.(i.le.IendR)).and.                      &
     &        ((JstrR.le.j).and.(j.le.JendR))) THEN
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
              DO k=1,N(ng)
                cff1=1.0_r8/(on_u(i,j)*                                 &
     &                       0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+       &
     &                               z_w(i  ,j,k)-z_w(i  ,j,k-1)))
!>              tl_u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*tl_cff1
!>
                ad_cff1=ad_cff1+SOURCES(ng)%Qsrc(is,k)*ad_u(i,j,k,nnew)
                ad_u(i,j,k,nnew)=0.0_r8
!>              tl_cff1=-cff1*cff1*on_u(i,j)*                           &
!>   &                  0.5_r8*(tl_z_w(i-1,j,k)-tl_z_w(i-1,j,k-1)+      &
!>   &                          tl_z_w(i  ,j,k)-tl_z_w(i  ,j,k-1))
!>
                adfac=-0.5_r8*cff1*cff1*on_u(i,j)*ad_cff1
                ad_z_w(i-1,j,k-1)=ad_z_w(i-1,j,k-1)-adfac
                ad_z_w(i  ,j,k-1)=ad_z_w(i  ,j,k-1)-adfac
                ad_z_w(i-1,j,k  )=ad_z_w(i-1,j,k  )+adfac
                ad_z_w(i  ,j,k  )=ad_z_w(i  ,j,k  )+adfac
                ad_cff1=0.0_r8
              END DO
            ELSE
              DO k=1,N(ng)
                cff1=1.0_r8/(om_v(i,j)*                                 &
     &                       0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+       &
     &                               z_w(i,j  ,k)-z_w(i,j  ,k-1)))
!>              tl_v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*tl_cff1
!>
                ad_cff1=ad_cff1+SOURCES(ng)%Qsrc(is,k)*ad_v(i,j,k,nnew)
                ad_v(i,j,k,nnew)=0.0_r8
!>              tl_cff1=-cff1*cff1*om_v(i,j)*                           &
!>   &                  0.5_r8*(tl_z_w(i,j-1,k)-tl_z_w(i,j-1,k-1)+      &
!>   &                          tl_z_w(i,j  ,k)-tl_z_w(i,j  ,k-1))
!>
                adfac=-0.5_r8*cff1*cff1*om_v(i,j)*ad_cff1
                ad_z_w(i,j-1,k-1)=ad_z_w(i,j-1,k-1)-adfac
                ad_z_w(i,j  ,k-1)=ad_z_w(i,j  ,k-1)-adfac
                ad_z_w(i,j-1,k  )=ad_z_w(i,j-1,k  )+adfac
                ad_z_w(i,j  ,k  )=ad_z_w(i,j  ,k  )+adfac
                ad_cff1=0.0_r8
              END DO
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Set lateral adjoint boundary conditions.
!-----------------------------------------------------------------------
!
!>    CALL tl_v3dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    nstp, nnew,                                   &
!>   &                    v)
!>
      CALL ad_v3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nnew,                                   &
     &                    ad_v)
!>    CALL tl_u3dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    nstp, nnew,                                   &
!>   &                    tl_u)
!>
      CALL ad_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nnew,                                   &
     &                    ad_u)
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the ETA-direction.
!-----------------------------------------------------------------------
!
      J_LOOP2 : DO j=Jstr,Jend
        IF (j.ge.JstrV) THEN
!
!  Compute BASIC STATE quantities.
!
          DO i=Istr,Iend
            AK(i,0)=0.5_r8*(Akv(i,j-1,0)+                               &
     &                      Akv(i,j  ,0))
            DO k=1,N(ng)
              AK(i,k)=0.5_r8*(Akv(i,j-1,k)+                             &
     &                        Akv(i,j  ,k))
              Hzk(i,k)=0.5_r8*(Hz(i,j-1,k)+                             &
     &                         Hz(i,j  ,k))
# if defined SPLINES_VVISC || defined DIAGNOSTICS_UV
              oHz(i,k)=1.0_r8/Hzk(i,k)
# endif
            END DO
          END DO
!
!  Couple and update new adjoint solution.
!
# if defined DIAGNOSTICS_UV && defined MASKING
!!        DO k=1,N(ng)
!!          DO i=Istr,Iend
!!            DO idiag=1,NDM3d
!!              DiaV3wrk(i,j,k,idiag)=DiaV3wrk(i,j,k,idiag)*vmask(i,j)
!!            END DO
!!          END DO
!!        END DO
# endif

          DO k=1,N(ng)
            DO i=Istr,Iend
# ifdef DIAGNOSTICS_UV
#  ifdef NEARSHORE_MELLOR
!!            DiaV3wrk(i,j,k,M3hrad)=DiaV3wrk(i,j,k,M3hrad)-            &
!!   &                               Dwrk(i,M2hrad)
#  endif
#  ifdef UV_ADV
!!            DiaV3wrk(i,j,k,M3hadv)=DiaV3wrk(i,j,k,M3hadv)-            &
!!   &                               Dwrk(i,M2hadv)
!!            DiaV3wrk(i,j,k,M3yadv)=DiaV3wrk(i,j,k,M3yadv)-            &
!!   &                               Dwrk(i,M2yadv)
!!            DiaV3wrk(i,j,k,M3xadv)=DiaV3wrk(i,j,k,M3xadv)-            &
!!   &                               Dwrk(i,M2xadv)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!            DiaV3wrk(i,j,k,M3hvis)=DiaV3wrk(i,j,k,M3hvis)-            &
!!   &                               Dwrk(i,M2hvis)
!!            DiaV3wrk(i,j,k,M3yvis)=DiaV3wrk(i,j,k,M3yvis)-            &
!!   &                               Dwrk(i,M2yvis)
!!            DiaV3wrk(i,j,k,M3xvis)=DiaV3wrk(i,j,k,M3xvis)-            &
!!   &                               Dwrk(i,M2xvis)
#  endif
#  ifdef UV_COR
!!            DiaV3wrk(i,j,k,M3fcor)=DiaV3wrk(i,j,k,M3fcor)-            &
!!   &                               Dwrk(i,M2fcor)
#  endif
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)-            &
!!   &                               Dwrk(i,M2bstr)
!!            DiaV3wrk(i,j,k,M3pgrd)=DiaV3wrk(i,j,k,M3pgrd)-            &
!!   &                               Dwrk(i,M2pgrd)
# endif
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>            tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*vmask_wet(i,j)
!>
              ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)*vmask_wet(i,j)
#  endif
#  ifdef MASKING
!>            tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*vmask(i,j)
!>
              ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)*vmask(i,j)
#  endif
!>            tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)-tl_DCs(i,0)
!>
              ad_DCs(i,0)=ad_DCs(i,0)-ad_v_stokes(i,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*vmask_wet(i,j)
!>
              ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*vmask_wet(i,j)
# endif
# ifdef MASKING
!>            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*vmask(i,j)
!>
              ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*vmask(i,j)
# endif
!>            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_DC(i,0)
!>
              ad_DC(i,0)=ad_DC(i,0)-ad_v(i,j,k,nnew)
            END DO
          END DO
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, vbar=DV_avg1/(D*om_v). Recall that, D=CF(:,0).
!
          DO i=Istr,Iend
            CF(i,0)=Hzk(i,1)
            DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              CF(i,0)=CF(i,0)+Hzk(i,k)
              DC(i,0)=DC(i,0)+v(i,j,k,nnew)*Hzk(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            DC1(i,0)=DC(i,0)                              ! intermediate
            cff1=1.0_r8/(CF(i,0)*om_v(i,j))
            DC(i,0)=(DC(i,0)*om_v(i,j)-DV_avg1(i,j))*cff1    ! recursive
# ifdef NEARSHORE_MELLOR
            DCs1(i)=DCs(i,0)                              ! intermediate
            cff2=1.0_r8/CF(i,0)
            DCs(i,0)=DCs(i,0)*cff2-vbar_stokes(i,j)          ! recursive
# endif
# ifdef DIAGNOSTICS_UV
!!          Dwrk(i,M2bstr)=(Dwrk(i,M2bstr)*om_v(i,j)-                   &
!!   &                      DiaV2wrk(i,j,M2bstr)-DiaV2wrk(i,j,M2sstr))* &
!!   &                     cff1
!!          DO idiag=1,M2pgrd
!!            Dwrk(i,idiag)=(Dwrk(i,idiag)*om_v(i,j)-                   &
!!   &                       DiaV2wrk(i,j,idiag))*cff1
!!          END DO
# endif
# ifdef NEARSHORE_MELLOR
!>          tl_DCs(i,0)=tl_DCs(i,0)*cff2+                               &
!>   &                  DCs1(i,0)*tl_cff2-tl_vbar_stokes(i,j)
!>
            ad_vbar_stokes(i,j)=ad_vbar_stokes(i,j)-ad_DCs(i,0)
            ad_cff2=ad_cff2+DCs1(i,0)*ad_DCs(i,0)
            ad_DCs(i,0)=cff2*ad_DCs(i,0)
!>          tl_cff2=-cff2*cff2*tl_CF(i,0)
!>
            ad_CF(i,0)=ad_CF(i,0)-cff2*cff2*tl_cff2
            ad_cff2=0.0_r8
# endif
!>          tl_DC(i,0)=(tl_DC(i,0)*om_v(i,j)-tl_DV_avg1(i,j))*cff1+     &
!>   &                 (DC1(i,0)*om_v(i,j)-DV_avg1(i,j))*tl_cff1
!>
            adfac=cff1*ad_DC(i,0)
            ad_DV_avg1(i,j)=ad_DV_avg1(i,j)-adfac
            ad_cff1=ad_cff1+                                            &
     &              (DC1(i,0)*om_v(i,j)-DV_avg1(i,j))*ad_DC(i,0)
            ad_DC(i,0)=om_v(i,j)*adfac
!>          tl_cff1=-cff1*cff1*tl_CF(i,0)*om_v(i,j)
!>
            ad_CF(i,0)=ad_CF(i,0)-cff1*cff1*om_v(i,j)*ad_cff1
            ad_cff1=0.0_r8
          END DO
          DO i=Istr,Iend
            CF(i,0)=Hzk(i,1)
            DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              CF(i,0)=CF(i,0)+Hzk(i,k)
              DC(i,0)=DC(i,0)+v(i,j,k,nnew)*Hzk(i,k)
# ifdef DIAGNOSTICS_UV
#  ifdef NEARSHORE_MELLOR
!!            Dwrk(i,M2hrad)=Dwrk(i,M2hrad)+                            &
!!   &                       DiaV3wrk(i,j,k,M3hrad)*Hzk(i,k)
#  endif
#  ifdef UV_ADV
!!            Dwrk(i,M2hadv)=Dwrk(i,M2hadv)+                            &
!!   &                       DiaV3wrk(i,j,k,M3hadv)*Hzk(i,k)
!!            Dwrk(i,M2yadv)=Dwrk(i,M2yadv)+                            &
!!   &                       DiaV3wrk(i,j,k,M3yadv)*Hzk(i,k)
!!            Dwrk(i,M2xadv)=Dwrk(i,M2xadv)+                            &
!!   &                       DiaV3wrk(i,j,k,M3xadv)*Hzk(i,k)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!            Dwrk(i,M2hvis)=Dwrk(i,M2hvis)+                            &
!!   &                       DiaV3wrk(i,j,k,M3hvis)*Hzk(i,k)
!!            Dwrk(i,M2yvis)=Dwrk(i,M2yvis)+                            &
!!   &                       DiaV3wrk(i,j,k,M3yvis)*Hzk(i,k)
!!            Dwrk(i,M2xvis)=Dwrk(i,M2xvis)+                            &
!!   &                       DiaV3wrk(i,j,k,M3xvis)*Hzk(i,k)
#  endif
#  ifdef UV_COR
!!            Dwrk(i,M2fcor)=Dwrk(i,M2fcor)+                            &
!!   &                       DiaV3wrk(i,j,k,M3fcor)*Hzk(i,k)
#  endif
!!            Dwrk(i,M2bstr)=Dwrk(i,M2bstr)+                            &
!!   &                       DiaV3wrk(i,j,k,M3vvis)*Hzk(i,k)
!!            Dwrk(i,M2pgrd)=Dwrk(i,M2pgrd)+                            &
!!   &                       DiaV3wrk(i,j,k,M3pgrd)*Hzk(i,k)
# endif
# ifdef NEARSHORE_MELLOR
!>            tl_DCs(i,0)=tl_DCs(i,0)+                                  &
!>   &                    tl_v_stokes(i,j,k)*Hzk(i,k)+                  &
!>   &                    v_stokes(i,j,k)*tl_Hzk(i,k)
!>
              ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)+Hzk(i,k)*ad_DCs(i,0)
              ad_Hzk(i,k)=ad_Hzk(i,k)+v_stokes(i,j,k)*ad_DCs(i,0)
# endif
!>            tl_DC(i,0)=tl_DC(i,0)+                                    &
!>   &                   tl_v(i,j,k,nnew)*Hzk(i,k)+                     &
!>   &                   v(i,j,k,nnew)*tl_Hzk(i,k)
!>
              ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)+ad_DC(i,0)*Hzk(i,k)
              ad_Hzk(i,k)=ad_Hzk(i,k)+v(i,j,k,nnew)*ad_DC(i,0)
!>            tl_CF(i,0)=tl_CF(i,0)+tl_Hzk(i,k)
!>
              ad_Hzk(i,k)=ad_Hzk(i,k)+ad_CF(i,0)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=Hzk(i,1)
            DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)
# ifdef DIAGNOSTICS_UV
#  ifdef NEARSHORE_MELLOR
!!          Dwrk(i,M2hrad)=DiaV3wrk(i,j,1,M3hrad)*Hzk(i,1)
#  endif
#  ifdef UV_ADV
!!          Dwrk(i,M2hadv)=DiaV3wrk(i,j,1,M3hadv)*Hzk(i,1)
!!          Dwrk(i,M2yadv)=DiaV3wrk(i,j,1,M3yadv)*Hzk(i,1)
!!          Dwrk(i,M2xadv)=DiaV3wrk(i,j,1,M3xadv)*Hzk(i,1)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!          Dwrk(i,M2hvis)=DiaV3wrk(i,j,1,M3hvis)*Hzk(i,1)
!!          Dwrk(i,M2yvis)=DiaV3wrk(i,j,1,M3yvis)*Hzk(i,1)
!!          Dwrk(i,M2xvis)=DiaV3wrk(i,j,1,M3xvis)*Hzk(i,1)
#  endif
#  ifdef UV_COR
!!          Dwrk(i,M2fcor)=DiaV3wrk(i,j,1,M3fcor)*Hzk(i,1)
#  endif
!!          Dwrk(i,M2bstr)=DiaV3wrk(i,j,1,M3vvis)*Hzk(i,1)
!!          Dwrk(i,M2pgrd)=DiaV3wrk(i,j,1,M3pgrd)*Hzk(i,1)
# endif
# ifdef NEARSHORE_MELLOR
!>          tl_DCs(i,0)=tl_v_stokes(i,j,1)*Hzk(i,1)+                    &
!>   &                  v_stokes(i,j,1)*tl_Hzk(i,1)
!>
            ad_v_stokes(i,j,1)=ad_v_stokes(i,j,1)+Hzk(i,1)*ad_DCs(i,0)
            ad_Hzk(i,1)=ad_Hzk(i,1)+v_stokes(i,j,1)*ad_DCs(i,0)
            ad_DCs(i,0)=0.0_r8
# endif
!>          tl_DC(i,0)=tl_v(i,j,1,nnew)*Hzk(i,1)+                       &
!>   &                 v(i,j,1,nnew)*tl_Hzk(i,1)
!>
            ad_v(i,j,1,nnew)=ad_v(i,j,1,nnew)+ad_DC(i,0)*Hzk(i,1)
            ad_Hzk(i,1)=ad_Hzk(i,1)+v(i,j,1,nnew)*ad_DC(i,0)
            ad_DC(i,0)=0.0_r8
!>          tl_CF(i,0)=tl_Hzk(i,1)
!>
            ad_Hzk(i,1)=ad_Hzk(i,1)+ad_CF(i,0)
            ad_CF(i,0)=0.0_r8
          END DO

# if defined SPLINES_VVISC
!
!  Apply spline algorithm to BASIC STATE "v" which should be
!  in units of velocity (m/s). DC will be used in the tangent
!  linear spline algorithm below. Solve BASIC STATE tridiagonal
!  system.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hzk(i,k  )-dt(ng)*AK(i,k-1)*oHz(i,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-dt(ng)*AK(i,k+1)*oHz(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(v(i,j,k+1,nnew)-v(i,j,k,nnew)-               &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution. Save DC for the adjoint below.
!  DC is scaled later by AK.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
!  Multiply DC (BASIC STATE spline gradients) by AK and save it in DC1.
!
          DO k=0,N(ng)
            DO i=Istr,Iend
              DC1(i,k)=DC(i,k)*AK(i,k)
            END DO
          END DO
!>        DO k=1,N(ng)
!>
          DO k=N(ng),1,-1
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+cff
#  endif
!>            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+tl_cff
!>
              ad_cff=ad_cff+ad_v(i,j,k,nnew)
!>            tl_cff=dt(ng)*(tl_oHz(i,k)*(DC(i,k)-DC(i,k-1))+           &
!>   &                       oHz(i,k)*(tl_DC(i,k)-tl_DC(i,k-1)))
!>                                                               use DC1
              adfac=dt(ng)*ad_cff
              adfac1=adfac*oHz(i,k)
              ad_DC(i,k-1)=ad_DC(i,k-1)-adfac1
              ad_DC(i,k  )=ad_DC(i,k  )+adfac1
              ad_oHz(i,k)=ad_oHz(i,k)+(DC1(i,k)-DC1(i,k-1))*adfac
              ad_cff=0.0_r8
!>            tl_DC(i,k)=tl_DC(i,k)*AK(i,k)+DC(i,k)*tl_AK(i,k)
!>                                                               use DC
              ad_AK(i,k)=ad_AK(i,k)+DC(i,k)*ad_DC(i,k)
              ad_DC(i,k)=AK(i,k)*ad_DC(i,k)
            END DO
          END DO
!
!  Adjoint back substitution.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
!>            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>
              ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
            END DO
          END DO
          DO i=Istr,Iend
!>          tl_DC(i,N(ng))=0.0_r8
!>
            ad_DC(i,N(ng))=0.0_r8
          END DO
!
!  Adjoint LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!>            tl_DC(i,k)=cff*(tl_v(i,j,k+1,nnew)-                       &
!>   &                        tl_v(i,j,k  ,nnew)-                       &
!>   &                        (tl_FC(i,k)*DC(i,k-1)+                    &
!>   &                         tl_BC(i,k)*DC(i,k  )+                    &
!>   &                         tl_CF(i,k)*DC(i,k+1))-                   &
!>   &                        FC(i,k)*tl_DC(i,k-1))
!>
              adfac=cff*ad_DC(i,k)
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
              ad_CF(i,k)=ad_CF(i,k)-DC(i,k+1)*adfac
              ad_BC(i,k)=ad_BC(i,k)-DC(i,k  )*adfac
              ad_FC(i,k)=ad_FC(i,k)-DC(i,k-1)*adfac
              ad_v(i,j,k,  nnew)=ad_v(i,j,k  ,nnew)-adfac
              ad_v(i,j,k+1,nnew)=ad_v(i,j,k+1,nnew)+adfac
              ad_DC(i,k)=0.0_r8
!>            tl_BC(i,k)=cff1*(tl_Hzk(i,k)+tl_Hzk(i,k+1))+              &
!>   &                   dt(ng)*(tl_AK(i,k)*(oHz(i,k)+oHz(i,k+1))+      &
!>   &                           AK(i,k)*(tl_oHz(i,k)+tl_oHz(i,k+1)))
!>
              adfac=dt(ng)*ad_BC(i,k)
              adfac1=adfac*AK(i,k)
              adfac2=cff1*ad_BC(i,k)
              ad_oHz(i,k  )=ad_oHz(i,k  )+adfac1
              ad_oHz(i,k+1)=ad_oHz(i,k+1)+adfac1
              ad_AK(i,k)=ad_AK(i,k)+(oHz(i,k)+oHz(i,k+1))*adfac
              ad_Hzk(i,k  )=ad_Hzk(i,k  )+adfac2
              ad_Hzk(i,k+1)=ad_Hzk(i,k+1)+adfac2
              ad_BC(i,k)=0.0_r8
            END DO
          END DO
!
!  Use conservative, parabolic spline reconstruction of adjoint
!  vertical viscosity derivatives.  Then, time step vertical
!  viscosity term implicitly by solving a tridiagonal system.
!
          DO i=Istr,Iend
!>          tl_DC(i,0)=0.0_r8
!>
            ad_DC(i,0)=0.0_r8
!>          tl_CF(i,0)=0.0_r8
!>
            ad_CF(i,0)=0.0_r8
          END DO
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
!>            tl_CF(i,k)=cff1*tl_Hzk(i,k+1)-                            &
!>   &                   dt(ng)*(tl_AK(i,k+1)*oHz(i,k+1)+               &
!>   &                           AK(i,k+1)*tl_oHz(i,k+1))
!>
              adfac=dt(ng)*ad_CF(i,k)
              ad_oHz(i,k+1)=ad_oHz(i,k+1)-AK(i,k+1)*adfac
              ad_AK(i,k+1)=ad_AK(i,k+1)-oHz(i,k+1)*adfac
              ad_Hzk(i,k+1)=ad_Hzk(i,k+1)+cff1*ad_CF(i,k)
              ad_CF(i,k)=0.0_r8
!>            tl_FC(i,k)=cff1*tl_Hzk(i,k  )-
!>   &                   dt(ng)*(tl_AK(i,k-1)*oHz(i,k  )+               &
!>   &                           AK(i,k-1)*tl_oHz(i,k  ))
!>
              adfac=dt(ng)*ad_FC(i,k)
              ad_oHz(i,k)=ad_oHz(i,k)-AK(i,k-1)*adfac
              ad_AK(i,k-1)=ad_AK(i,k-1)-oHz(i,k)*adfac
              ad_Hzk(i,k)=ad_Hzk(i,k)+cff1*ad_FC(i,k)
              ad_FC(i,k)=0.0_r8
            END DO
          END DO
# else
!
!  Compute BASIC STATE off-diagonal coefficients [lambda*dt*Akv/Hz]
!  for the implicit vertical viscosity term at horizontal V-points
!  and vertical W-points.
!
!  NOTE: The original code solves the tridiagonal system A*v=r where
!        A is a matrix and v and r are vectors. We need to solve the
!        tangent linear of this system which is A*tl_v+tl_A*v=tl_r.
!        Here, tl_A*v and tl_r are known, so we must solve for tl_v
!        by inverting A*tl_v=tl_r-tl_A*v.
!
          cff=-lambda*dt(ng)/0.5_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-                 &
     &                     z_r(i,j,k  )-z_r(i,j-1,k  ))
              FC(i,k)=cff*cff1*AK(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
            END DO
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
            END DO
          END DO
!
!  Compute new adjoint solution by back substitution.
!  (DC is a tangent linear variable here).
!
!>        DO k=N(ng)-1,1,-1
!>
          DO k=1,N(ng)-1
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+            &
!!   &                               v(i,j,k,nnew)-wrk(i,k)
#  endif
!>            tl_v(i,j,k,nnew)=DC(i,k)
!>
              ad_DC(i,k)=ad_DC(i,k)+ad_v(i,j,k,nnew)
              ad_v(i,j,k,nnew)=0.0_r8
!>            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>
              ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
#  ifdef DIAGNOSTICS_UV
!!            wrk(i,k)=v(i,j,k,nnew)*oHz(i,k)
#  endif
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaV3wrk(i,j,N(ng),M3vvis)=DiaV3wrk(i,j,N(ng),M3vvis)+      &
!!   &                                 v(i,j,N(ng),nnew)-wrk(i,N(ng))
#  endif
!>          tl_v(i,j,N(ng),nnew)=DC(i,N(ng))
!>
            ad_DC(i,N(ng))=ad_DC(i,N(ng))+ad_v(i,j,N(ng),nnew)
            ad_v(i,j,N(ng),nnew)=0.0_r8
!>          DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/      &
!>   &                  (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
!>
            adfac=ad_DC(i,N(ng))/                                       &
     &            (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
            ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-FC(i,N(ng)-1)*adfac
            ad_DC(i,N(ng)  )=adfac
          END DO
!
!  Solve the adjoint tridiagonal system.
!  (DC is a tangent linear variable here).
!
!>        DO k=2,N(ng)-1
!>
          DO k=N(ng)-1,2,-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
!>            DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
!>
              adfac=cff*ad_DC(i,k)
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k-1)*adfac
              ad_DC(i,k  )=adfac
            END DO
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
!>          DC(i,1)=cff*DC(i,1)
!>
            ad_DC(i,1)=cff*ad_DC(i,1)
          END DO
          DO i=Istr,Iend
!>          DC(i,N(ng))=tl_v(i,j,N(ng),nnew)-                           &
!>   &                  (tl_FC(i,N(ng)-1)*v(i,j,N(ng)-1,nnew)+          &
!>   &                   tl_BC(i,N(ng)  )*v(i,j,N(ng)  ,nnew))
!>
            ad_BC(i,N(ng)  )=-v(i,j,N(ng)  ,nnew)*ad_DC(i,N(ng))
            ad_FC(i,N(ng)-1)=-v(i,j,N(ng)-1,nnew)*ad_DC(i,N(ng))
            ad_v(i,j,N(ng),nnew)=ad_DC(i,N(ng))
            ad_DC(i,N(ng))=0.0_r8
!>          DC(i,1)=tl_v(i,j,1,nnew)-                                   &
!>   &              (tl_BC(i,1)*v(i,j,1,nnew)+                          &
!>   &               tl_FC(i,1)*v(i,j,2,nnew))
!>
            ad_FC(i,1)=-v(i,j,2,nnew)*ad_DC(i,1)
            ad_BC(i,1)=-v(i,j,1,nnew)*ad_DC(i,1)
            ad_v(i,j,1,nnew)=ad_DC(i,1)
            ad_DC(i,1)=0.0_r8
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
!>            DC(i,k)=tl_v(i,j,k,nnew)-                                 &
!>   &                (tl_FC(i,k-1)*v(i,j,k-1,nnew)+                    &
!>   &                 tl_BC(i,k  )*v(i,j,k  ,nnew)+                    &
!>   &                 tl_FC(i,k  )*v(i,j,k+1,nnew))
!>
              ad_FC(i,k-1)=ad_FC(i,k-1)-v(i,j,k-1,nnew)*ad_DC(i,k)
              ad_FC(i,k  )=ad_FC(i,k  )-v(i,j,k+1,nnew)*ad_DC(i,k)
              ad_BC(i,k  )=ad_BC(i,k  )-v(i,j,k  ,nnew)*ad_DC(i,k)
              ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)+ad_DC(i,k)
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_BC(i,k)=tl_Hzk(i,k)-tl_FC(i,k)-tl_FC(i,k-1)
!>
              ad_Hzk(i,k)=ad_Hzk(i,k)+ad_BC(i,k)
              ad_FC(i,k-1)=ad_FC(i,k-1)-ad_BC(i,k)
              ad_FC(i,k  )=ad_FC(i,k  )-ad_BC(i,k)
              ad_BC(i,k)=0.0_r8
            END DO
          END DO
!
!  Compute adjoint off-diagonal coefficients [lambda*dt*Akv/Hz] for
!  the implicit vertical viscosity term at horizontal V-points and
!  vertical W-points.
!
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
          END DO
          cff=-lambda*dt(ng)/0.5_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-                 &
     &                     z_r(i,j,k  )-z_r(i,j-1,k  ))
!>            tl_FC(i,k)=cff*(tl_cff1*AK(i,k)+cff1*tl_AK(i,k))
!>
              adfac=cff*ad_FC(i,k)
              ad_AK(i,k)=ad_AK(i,k)+cff1*adfac
              ad_cff1=ad_cff1+AK(i,k)*adfac
              ad_FC(i,k)=0.0_r8
!>            tl_cff1=-cff1*cff1*(tl_z_r(i,j,k+1)+tl_z_r(i,j-1,k+1)-    &
!>   &                            tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  ))
!>
              adfac=-cff1*cff1*ad_cff1
              ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-adfac
              ad_z_r(i,j  ,k  )=ad_z_r(i,j  ,k  )-adfac
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)+adfac
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+adfac
              ad_cff1=0.0_r8
            END DO
          END DO
# endif
!
!  Time step adjoint right-hand-side terms.
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.25_r8*dt(ng)
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
            cff=0.25_r8*dt(ng)*3.0_r8/2.0_r8
          ELSE
            cff=0.25_r8*dt(ng)*23.0_r8/12.0_r8
          END IF
          DO i=Istr,Iend
            DC(i,0)=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          END DO
!
!  The BASIC STATE "v" used below must be in transport units, but "v"
!  retrived is in m/s so we multiply by "Hzk".
!
          DO k=1,N(ng)
            DO i=Istr,Iend
# ifdef DIAGNOSTICS_UV
!!            DiaV3wrk(i,j,k,M3rate)=DiaV3wrk(i,j,k,M3rate)*oHz(i,k)
#  ifdef BODYFORCE
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+            &
!!   &                               DC(i,0)*DiaRV(i,j,k,nrhs,M3vvis)*  &
!!   &                               oHz(i,k)
#  endif
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)*oHz(i,k)
#  if defined UV_VIS2 || defined UV_VIS4
!!            DiaV3wrk(i,j,k,M3hvis)=DiaV3wrk(i,j,k,M3hvis)*oHz(i,k)
!!            DiaV3wrk(i,j,k,M3yvis)=DiaV3wrk(i,j,k,M3yvis)*oHz(i,k)
!!            DiaV3wrk(i,j,k,M3xvis)=DiaV3wrk(i,j,k,M3xvis)*oHz(i,k)
#  endif
!!            DO idiag=1,M3pgrd
!!              DiaV3wrk(i,j,k,idiag)=(DiaV3wrk(i,j,k,idiag)+           &
!!   &                                 DC(i,0)*                         &
!!   &                                 DiaRV(i,j,k,nrhs,idiag))*        &
!!   &                                oHz(i,k)
!!            END DO
# endif
# ifdef SPLINES_VVISC
!>            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*oHz(i,k)+               &
!>   &                         (v(i,j,k,nnew)*Hzk(i,k))*tl_oHz(i,k)
!>
              ad_oHz(i,k)=ad_oHz(i,k)+                                  &
     &                    (v(i,j,k,nnew)*Hzk(i,k))*ad_v(i,j,k,nnew)
              ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*oHz(i,k)
# endif
!>            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+                        &
!>   &                         DC(i,0)*tl_rv(i,j,k,nrhs)
!>
              ad_rv(i,j,k,nrhs)=ad_rv(i,j,k,nrhs)+                      &
     &                          DC(i,0)*ad_v(i,j,k,nnew)
            END DO
          END DO
          DO i=Istr,Iend
            DO k=1,N(ng)
# ifdef SPLINES_VVISC
              oHz(i,k)=1.0_r8/Hzk(i,k)
!>            tl_oHz(i,k)=-oHz(i,k)*oHz(i,k)*tl_Hzk(i,k)
!>
              ad_Hzk(i,k)=ad_Hzk(i,k)-oHz(i,k)*oHz(i,k)*ad_oHz(i,k)
              ad_oHz(i,k)=0.0_r8
# endif
!>            tl_Hzk(i,k)=0.5_r8*(tl_Hz(i,j-1,k)+                       &
!>   &                            tl_Hz(i,j  ,k))
!>
              adfac=0.5_r8*ad_Hzk(i,k)
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac
              ad_Hzk(i,k)=0.0_r8
!>            tl_AK(i,k)=0.5_r8*(tl_Akv(i,j-1,k)+                       &
!>   &                           tl_Akv(i,j  ,k))
!>
              adfac=0.5_r8*ad_AK(i,k)
              ad_Akv(i,j-1,k)=ad_Akv(i,j-1,k)+adfac
              ad_Akv(i,j  ,k)=ad_Akv(i,j  ,k)+adfac
              ad_AK(i,k)=0.0_r8
            END DO
!>          tl_AK(i,0)=0.5_r8*(tl_Akv(i,j-1,0)+                         &
!>   &                         tl_Akv(i,j  ,0))
!>
            adfac=0.5_r8*ad_AK(i,0)
            ad_Akv(i,j-1,0)=ad_Akv(i,j-1,0)+adfac
            ad_Akv(i,j  ,0)=ad_Akv(i,j  ,0)+adfac
            ad_AK(i,0)=0.0_r8
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the XI-direction.
!-----------------------------------------------------------------------
!
!  Compute BASIC STATE quatities.
!
        DO i=IstrU,Iend
          AK(i,0)=0.5_r8*(Akv(i-1,j,0)+                                 &
     &                    Akv(i  ,j,0))
          DO k=1,N(ng)
            AK(i,k)=0.5_r8*(Akv(i-1,j,k)+                               &
     &                      Akv(i  ,j,k))
            Hzk(i,k)=0.5_r8*(Hz(i-1,j,k)+                               &
     &                       Hz(i  ,j,k))
# ifdef SPLINES_VVISC
            oHz(i,k)=1.0_r8/Hzk(i,k)
# endif
          END DO
        END DO
!
!  Couple and update new adjoint solution.
!
# if defined DIAGNOSTICS_UV && defined MASKING
!!      DO k=1,N(ng)
!!        DO i=IstrU,Iend
!!          DO idiag=1,NDM3d
!!            DiaU3wrk(i,j,k,idiag)=DiaU3wrk(i,j,k,idiag)*umask(i,j)
!!          END DO
!!        END DO
!!      END DO
# endif

        DO k=1,N(ng)
          DO i=IstrU,Iend
# ifdef DIAGNOSTICS_UV
#  ifdef NEARSHORE_MELLOR
!!          DiaU3wrk(i,j,k,M3hrad)=DiaU3wrk(i,j,k,M3hrad)-              &
!!   &                             Dwrk(i,M2hrad)
#  endif
#  ifdef UV_ADV
!!          DiaU3wrk(i,j,k,M3hadv)=DiaU3wrk(i,j,k,M3hadv)-              &
!!   &                             Dwrk(i,M2hadv)
!!          DiaU3wrk(i,j,k,M3yadv)=DiaU3wrk(i,j,k,M3yadv)-              &
!!   &                             Dwrk(i,M2yadv)
!!          DiaU3wrk(i,j,k,M3xadv)=DiaU3wrk(i,j,k,M3xadv)-              &
!!   &                             Dwrk(i,M2xadv)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!          DiaU3wrk(i,j,k,M3hvis)=DiaU3wrk(i,j,k,M3hvis)-              &
!!   &                             Dwrk(i,M2hvis)
!!          DiaU3wrk(i,j,k,M3yvis)=DiaU3wrk(i,j,k,M3yvis)-              &
!!   &                             Dwrk(i,M2yvis)
!!          DiaU3wrk(i,j,k,M3xvis)=DiaU3wrk(i,j,k,M3xvis)-              &
!!   &                             Dwrk(i,M2xvis)
#  endif
#  ifdef UV_COR
!!          DiaU3wrk(i,j,k,M3fcor)=DiaU3wrk(i,j,k,M3fcor)-              &
!!   &                             Dwrk(i,M2fcor)
#  endif
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)-              &
!!   &                             Dwrk(i,M2bstr)
!!          DiaU3wrk(i,j,k,M3pgrd)=DiaU3wrk(i,j,k,M3pgrd)-              &
!!   &                             Dwrk(i,M2pgrd)
# endif
# ifdef NEARSHORE_MELLOR
#  ifdef WET_DRY_NOT_YET
!>          tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*umask_wet(i,j)
!>
            ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)*umask_wet(i,j)
#  endif
#  ifdef MASKING
!>          tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*umask(i,j)
!>
            ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)*umask(i,j)
#  endif
!>          tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)-tl_DCs(i,0)
!>
            ad_DCs(i,0)=ad_DCs(i,0)-ad_u_stokes(i,j,k)
# endif
# ifdef WET_DRY_NOT_YET
!>          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*umask_wet(i,j)
!>
            ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*umask_wet(i,j)
# endif
# ifdef MASKING
!>          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*umask(i,j)
!>
            ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*umask(i,j)
# endif
!>          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_DC(i,0)
!>
            ad_DC(i,0)=ad_DC(i,0)-ad_u(i,j,k,nnew)
          END DO
        END DO
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, ubar=DU_avg1/(D*on_u). Recall that, D=CF(:,0).
!
        DO i=IstrU,Iend
          CF(i,0)=Hzk(i,1)
          DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            CF(i,0)=CF(i,0)+Hzk(i,k)
            DC(i,0)=DC(i,0)+u(i,j,k,nnew)*Hzk(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
          DC1(i,0)=DC(i,0)                                ! intermediate
          cff1=1.0/(CF(i,0)*on_u(i,j))
          DC(i,0)=(DC(i,0)*on_u(i,j)-DU_avg1(i,j))*cff1      ! recursive
# ifdef NEARSHORE_MELLOR
          DCs1(i)=DCs(i,0)                                ! intermediate
          cff2=1.0_r8/CF(i,0)
          DCs(i,0)=DCs(i,0)*cff2-ubar_stokes(i,j)            ! recursive
# endif
# ifdef DIAGNOSTICS_UV
!!        Dwrk(i,M2bstr)=(Dwrk(i,M2bstr)*on_u(i,j)-                     &
!!   &                    DiaU2wrk(i,j,M2bstr)-DiaU2wrk(i,j,M2sstr))*   &
!!   &                   cff1
!!        DO idiag=1,M2pgrd
!!          Dwrk(i,idiag)=(Dwrk(i,idiag)*on_u(i,j)-                     &
!!   &                     DiaU2wrk(i,j,idiag))*cff1
!!        END DO
# endif
# ifdef NEARSHORE_MELLOR
!>        tl_DCs(i,0)=tl_DCs(i,0)*cff2+                                 &
!>   &                DCs1(i)*tl_cff2-tl_ubar_stokes(i,j)
!>
          ad_cff2=ad_cff2+DCs1(i)*ad_DCs(i,0)
          ad_ubar_stokes(i,j)=ad_ubar_stokes(i,j)-ad_DCs(i,0)
          ad_DCs(i,0)=cff2*ad_DCs(i,0)
!>        tl_cff2=-cff2*cff2*tl_CF(i,0)
!>
          ad_CF(i,0)=ad_CF(i,0)-cff2*cff2*ad_cff2
          ad_cff2=0.0_r8
# endif
!>        tl_DC(i,0)=(tl_DC(i,0)*on_u(i,j)-tl_DU_avg1(i,j))*cff1+       &
!>   &               (DC1(i,0)*on_u(i,j)-DU_avg1(i,j))*tl_cff1
!>
          adfac=cff1*ad_DC(i,0)
          ad_DU_avg1(i,j)=ad_DU_avg1(i,j)-adfac
          ad_cff1=ad_cff1+                                              &
     &            (DC1(i,0)*on_u(i,j)-DU_avg1(i,j))*ad_DC(i,0)
          ad_DC(i,0)=on_u(i,j)*adfac
!>        tl_cff1=-cff1*cff1*tl_CF(i,0)*on_u(i,j)
!>
          ad_CF(i,0)=ad_CF(i,0)-cff1*cff1*on_u(i,j)*ad_cff1
          ad_cff1=0.0_r8
        END DO
        DO i=IstrU,Iend
          CF(i,0)=Hzk(i,1)
          DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            CF(i,0)=CF(i,0)+Hzk(i,k)
            DC(i,0)=DC(i,0)+u(i,j,k,nnew)*Hzk(i,k)
# ifdef NEARSHORE_MELLOR
            DCs(i,0)=DCs(i,0)+u_stokes(i,j,k)*Hzk(i,k)
# endif
# ifdef DIAGNOSTICS_UV
#  ifdef NEARSHORE_MELLOR
!!          Dwrk(i,M2hrad)=Dwrk(i,M2hrad)+                              &
!!   &                     DiaU3wrk(i,j,k,M3hrad)*Hzk(i,k)
#  endif
#  ifdef UV_ADV
!!          Dwrk(i,M2hadv)=Dwrk(i,M2hadv)+                              &
!!   &                     DiaU3wrk(i,j,k,M3hadv)*Hzk(i,k)
!!          Dwrk(i,M2yadv)=Dwrk(i,M2yadv)+                              &
!!   &                     DiaU3wrk(i,j,k,M3yadv)*Hzk(i,k)
!!          Dwrk(i,M2xadv)=Dwrk(i,M2xadv)+                              &
!!   &                     DiaU3wrk(i,j,k,M3xadv)*Hzk(i,k)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!          Dwrk(i,M2hvis)=Dwrk(i,M2hvis)+                              &
!!   &                     DiaU3wrk(i,j,k,M3hvis)*Hzk(i,k)
!!          Dwrk(i,M2yvis)=Dwrk(i,M2yvis)+                              &
!!   &                     DiaU3wrk(i,j,k,M3yvis)*Hzk(i,k)
!!          Dwrk(i,M2xvis)=Dwrk(i,M2xvis)+                              &
!!   &                     DiaU3wrk(i,j,k,M3xvis)*Hzk(i,k)
#  endif
#  ifdef UV_COR
!!          Dwrk(i,M2fcor)=Dwrk(i,M2fcor)+                              &
!!   &                     DiaU3wrk(i,j,k,M3fcor)*Hzk(i,k)
#  endif
!!          Dwrk(i,M2bstr)=Dwrk(i,M2bstr)+                              &
!!   &                     DiaU3wrk(i,j,k,M3vvis)*Hzk(i,k)
!!          Dwrk(i,M2pgrd)=Dwrk(i,M2pgrd)+                              &
!!   &                     DiaU3wrk(i,j,k,M3pgrd)*Hzk(i,k)
# endif
# ifdef NEARSHORE_MELLOR
!>          tl_DCs(i,0)=tl_DCs(i,0)+                                    &
!>   &                  tl_u_stokes(i,j,k)*Hzk(i,k)+                    &
!>   &                  u_stokes(i,j,k)*tl_Hzk(i,k)
!>
            ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)+Hzk(i,k)*ad_DCs(i,0)
            ad_Hzk(i,k)=ad_Hzk(i,k)+u_stokes(i,j,k)*ad_DCs(i,0)
# endif
!>          tl_DC(i,0)=tl_DC(i,0)+                                      &
!>   &                 tl_u(i,j,k,nnew)*Hzk(i,k)+                       &
!>   &                 u(i,j,k,nnew)*tl_Hzk(i,k)
!>
            ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)+ad_DC(i,0)*Hzk(i,k)
            ad_Hzk(i,k)=ad_Hzk(i,k)+u(i,j,k,nnew)*ad_DC(i,0)
!>          tl_CF(i,0)=tl_CF(i,0)+tl_Hzk(i,k)
!>
            ad_Hzk(i,k)=ad_Hzk(i,k)+ad_CF(i,0)
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,0)=Hzk(i,1)
          DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)
# ifdef NEARSHORE_MELLOR
          DCs(i,0)=u_stokes(i,j,1)*Hzk(i,1)
# endif
# ifdef DIAGNOSTICS_UV
#  ifdef NEARSHORE_MELLOR
!!        Dwrk(i,M2hrad)=DiaU3wrk(i,j,1,M3hrad)*Hzk(i,1)
#  endif
#  ifdef UV_ADV
!!        Dwrk(i,M2hadv)=DiaU3wrk(i,j,1,M3hadv)*Hzk(i,1)
!!        Dwrk(i,M2yadv)=DiaU3wrk(i,j,1,M3yadv)*Hzk(i,1)
!!        Dwrk(i,M2xadv)=DiaU3wrk(i,j,1,M3xadv)*Hzk(i,1)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!        Dwrk(i,M2hvis)=DiaU3wrk(i,j,1,M3hvis)*Hzk(i,1)
!!        Dwrk(i,M2yvis)=DiaU3wrk(i,j,1,M3yvis)*Hzk(i,1)
!!        Dwrk(i,M2xvis)=DiaU3wrk(i,j,1,M3xvis)*Hzk(i,1)
#  endif
#  ifdef UV_COR
!!        Dwrk(i,M2fcor)=DiaU3wrk(i,j,1,M3fcor)*Hzk(i,1)
#  endif
!!        Dwrk(i,M2bstr)=DiaU3wrk(i,j,1,M3vvis)*Hzk(i,1)
!!        Dwrk(i,M2pgrd)=DiaU3wrk(i,j,1,M3pgrd)*Hzk(i,1)
# endif
# ifdef NEARSHORE_MELLOR
!>        tl_DCs(i,0)=tl_u_stokes(i,j,1)*Hzk(i,1)+                      &
!>   &                u_stokes(i,j,1)*tl_Hzk(i,1)
!>
          ad_u_stokes(i,j,1)=ad_u_stokes(i,j,1)+Hzk(i,1)*ad_DCs(i,0)
          ad_Hzk(i,1)=ad_Hzk(i,1)+u_stokes(i,j,1)*ad_DCs(i,0)
          ad_DCs(i,0)=0.0_r8
# endif
!>        tl_DC(i,0)=tl_u(i,j,1,nnew)*Hzk(i,1)+                         &
!>   &               u(i,j,1,nnew)*tl_Hzk(i,1)
!>
          ad_u(i,j,1,nnew)=ad_u(i,j,1,nnew)+ad_DC(i,0)*Hzk(i,1)
          ad_Hzk(i,1)=ad_Hzk(i,1)+u(i,j,1,nnew)*ad_DC(i,0)
          ad_DC(i,0)=0.0_r8
!>        tl_CF(i,0)=tl_Hzk(i,1)
!>
          ad_Hzk(i,1)=ad_Hzk(i,1)+ad_CF(i,0)
          ad_CF(i,0)=0.0_r8
        END DO

# if defined SPLINES_VVISC
!
!  Apply spline algorithm to BASIC STATE "u" which should be
!  in units of velocity (m/s). DC will be used in the tangent
!  linear spline algorithm below.  Solve BASIC STATE tridiagonal
!  system.
!
        cff1=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=cff1*Hzk(i,k  )-dt(ng)*AK(i,k-1)*oHz(i,k  )
            CF(i,k)=cff1*Hzk(i,k+1)-dt(ng)*AK(i,k+1)*oHz(i,k+1)
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,0)=0.0_r8
          DC(i,0)=0.0_r8
        END DO
!
!  LU decomposition and forward substitution.
!
        cff1=1.0_r8/3.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                         &
     &              dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
            cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
            CF(i,k)=cff*CF(i,k)
            DC(i,k)=cff*(u(i,j,k+1,nnew)-u(i,j,k,nnew)-                 &
     &                   FC(i,k)*DC(i,k-1))
          END DO
        END DO
!
!  Backward substitution. Save DC for the adjoint below.
!  DC is scaled later by AK.
!
        DO i=IstrU,Iend
          DC(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
          END DO
        END DO
!
!  Multiply DC (BASIC STATE spline gradients) by AK and save it in DC1.
!
        DO k=0,N(ng)
          DO i=IstrU,Iend
            DC1(i,k)=DC(i,k)*AK(i,k)
          END DO
        END DO
!>      DO k=1,N(ng)
!>
        DO k=N(ng),1,-1
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+cff
#  endif
!>          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+tl_cff
!>
            ad_cff=ad_cff+ad_u(i,j,k,nnew)
!>          tl_cff=dt(ng)*(tl_oHz(i,k)*(DC(i,k)-DC(i,k-1))+             &
!>   &                     oHz(i,k)*(tl_DC(i,k)-tl_DC(i,k-1)))
!>                                                               use DC1
            adfac=dt(ng)*ad_cff
            adfac1=adfac*oHz(i,k)
            ad_DC(i,k-1)=ad_DC(i,k-1)-adfac1
            ad_DC(i,k  )=ad_DC(i,k  )+adfac1
            ad_oHz(i,k)=ad_oHz(i,k)+(DC1(i,k)-DC1(i,k-1))*adfac
            ad_cff=0.0_r8
!>          tl_DC(i,k)=tl_DC(i,k)*AK(i,k)+DC(i,k)*tl_AK(i,k)
!>                                                               use DC
            ad_AK(i,k)=ad_AK(i,k)+DC(i,k)*ad_DC(i,k)
            ad_DC(i,k)=AK(i,k)*ad_DC(i,k)
          END DO
        END DO
!
!  Adjoint back substitution.
!
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
!>          tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!>
            ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
!>        tl_DC(i,N(ng))=0.0_r8
!>
          ad_DC(i,N(ng))=0.0_r8
        END DO
!
!  Adjoint LU decomposition and forward substitution.
!
        cff1=1.0_r8/3.0_r8
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                         &
     &              dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
            cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!>          tl_DC(i,k)=cff*(tl_u(i,j,k+1,nnew)-                         &
!>   &                      tl_u(i,j,k  ,nnew)-                         &
!>   &                      (tl_FC(i,k)*DC(i,k-1)+                      &
!>   &                       tl_BC(i,k)*DC(i,k  )+                      &
!>   &                       tl_CF(i,k)*DC(i,k+1))-                     &
!>   &                      FC(i,k)*tl_DC(i,k-1))
!>
            adfac=cff*ad_DC(i,k)
            ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
            ad_CF(i,k)=ad_CF(i,k)-DC(i,k+1)*adfac
            ad_BC(i,k)=ad_BC(i,k)-DC(i,k  )*adfac
            ad_FC(i,k)=ad_FC(i,k)-DC(i,k-1)*adfac
            ad_u(i,j,k  ,nnew)=ad_u(i,j,k  ,nnew)-adfac
            ad_u(i,j,k+1,nnew)=ad_u(i,j,k+1,nnew)+adfac
            ad_DC(i,k)=0.0_r8
!>          tl_BC(i,k)=cff1*(tl_Hzk(i,k)+tl_Hzk(i,k+1))+                &
!>   &                 dt(ng)*(tl_AK(i,k)*(oHz(i,k)+oHz(i,k+1))+        &
!>   &                         AK(i,k)*(tl_oHz(i,k)+tl_oHz(i,k+1)))
!>
            adfac=dt(ng)*ad_BC(i,k)
            adfac1=adfac*AK(i,k)
            adfac2=cff1*ad_BC(i,k)
            ad_oHz(i,k  )=ad_oHz(i,k  )+adfac1
            ad_oHz(i,k+1)=ad_oHz(i,k+1)+adfac1
            ad_AK(i,k)=ad_AK(i,k)+(oHz(i,k)+oHz(i,k+1))*adfac
            ad_Hzk(i,k  )=ad_Hzk(i,k  )+adfac2
            ad_Hzk(i,k+1)=ad_Hzk(i,k+1)+adfac2
            ad_BC(i,k)=0.0_r8
          END DO
        END DO
!
!  Use conservative, parabolic spline reconstruction of adjoint
!  vertical viscosity derivatives.  Then, time step vertical
!  viscosity term implicitly by solving a tridiagonal system.
!
        DO i=IstrU,Iend
!>        tl_DC(i,0)=0.0_r8
!>
          ad_DC(i,0)=0.0_r8
!>        tl_CF(i,0)=0.0_r8
!>
          ad_CF(i,0)=0.0_r8
        END DO
        cff1=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
!>          tl_CF(i,k)=cff1*tl_Hzk(i,k+1)-                              &
!>   &                 dt(ng)*(tl_AK(i,k+1)*oHz(i,k+1)+                 &
!>   &                         AK(i,k+1)*tl_oHz(i,k+1))
!>
            adfac=dt(ng)*ad_CF(i,k)
            ad_oHz(i,k+1)=ad_oHz(i,k+1)-AK(i,k+1)*adfac
            ad_AK(i,k+1)=ad_AK(i,k+1)-oHz(i,k+1)*adfac
            ad_Hzk(i,k+1)=ad_Hzk(i,k+1)+cff1*ad_CF(i,k)
            ad_CF(i,k)=0.0_r8
!>          tl_FC(i,k)=cff1*tl_Hzk(i,k  )-                              &
!>   &                 dt(ng)*(tl_AK(i,k-1)*oHz(i,k  )+                 &
!>   &                         AK(i,k-1)*tl_oHz(i,k  ))
!>
            adfac=dt(ng)*ad_FC(i,k)
            ad_oHz(i,k)=ad_oHz(i,k)-AK(i,k-1)*adfac
            ad_AK(i,k-1)=ad_AK(i,k-1)-oHz(i,k)*adfac
            ad_Hzk(i,k)=ad_Hzk(i,k)+cff1*ad_FC(i,k)
            ad_FC(i,k)=0.0_r8
          END DO
        END DO
# else
!
!  Compute BASIC STATE off-diagonal coefficients [lambda*dt*Akv/Hz]
!  for the implicit vertical viscosity term at horizontal U-points
!  and vertical W-points.
!
!  NOTE: The original code solves the tridiagonal system A*u=r where
!        A is a matrix and u and r are vectors. We need to solve the
!        tangent linear of this system which is A*tl_u+tl_A*u=tl_r.
!        Here, tl_A*u and tl_r are known, so we must solve for tl_u
!        by inverting A*tl_u=tl_r-tl_A*u.
!
        cff=-lambda*dt(ng)/0.5_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff1=1.0_r8/(z_r(i,j,k+1)+z_r(i-1,j,k+1)-                   &
     &                   z_r(i,j,k  )-z_r(i-1,j,k  ))
            FC(i,k)=cff*cff1*AK(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          FC(i,N(ng))=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
          END DO
        END DO
        DO i=IstrU,Iend
          cff=1.0_r8/BC(i,1)
          CF(i,1)=cff*FC(i,1)
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
            CF(i,k)=cff*FC(i,k)
          END DO
        END DO
!
!  Compute new adjoint solution by back substitution.
!  (DC is a tangent linear variable here).
!
!>      DO k=N(ng)-1,1,-1
!>
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+              &
!!   &                             u(i,j,k,nnew)-wrk(i,k)
#  endif
!>          tl_u(i,j,k,nnew)=DC(i,k)
!>
            ad_DC(i,k)=ad_DC(i,k)+ad_u(i,j,k,nnew)
            ad_u(i,j,k,nnew)=0.0_r8
!>          DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>
            ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
#  ifdef DIAGNOSTICS_UV
!!          wrk(i,k)=u(i,j,k,nnew)*oHz(i,k)
#  endif
          END DO
        END DO
        DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!        DiaU3wrk(i,j,N(ng),M3vvis)=DiaU3wrk(i,j,N(ng),M3vvis)+        &
!!   &                               u(i,j,N(ng),nnew)-wrk(i,N(ng))
#  endif
!>        tl_u(i,j,N(ng),nnew)=DC(i,N(ng))
!>
          ad_DC(i,N(ng))=ad_DC(i,N(ng))+ad_u(i,j,N(ng),nnew)
          ad_u(i,j,N(ng),nnew)=0.0_r8
!>        DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/        &
!>   &                (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
!>
          adfac=ad_DC(i,N(ng))/                                         &
     &          (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
          ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-FC(i,N(ng)-1)*adfac
          ad_DC(i,N(ng)  )=adfac
        END DO
!
!  Solve the adjoint tridiagonal system.
!  (DC is a tangent linear variable here).
!
!>      DO k=2,N(ng)-1
!>
        DO k=N(ng)-1,2,-1
          DO i=IstrU,Iend
            cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
!>          DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
!>
            adfac=cff*ad_DC(i,k)
            ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k-1)*adfac
            ad_DC(i,k  )=adfac
          END DO
        END DO
        DO i=IstrU,Iend
          cff=1.0_r8/BC(i,1)
!>        DC(i,1)=cff*DC(i,1)
!>
          ad_DC(i,1)=cff*ad_DC(i,1)
        END DO
        DO i=IstrU,Iend
!>        DC(i,N(ng))=tl_u(i,j,N(ng),nnew)-                             &
!>   &                (tl_FC(i,N(ng)-1)*u(i,j,N(ng)-1,nnew)+            &
!>   &                 tl_BC(i,N(ng)  )*u(i,j,N(ng)  ,nnew))
!>
          ad_BC(i,N(ng)  )=-u(i,j,N(ng)  ,nnew)*ad_DC(i,N(ng))
          ad_FC(i,N(ng)-1)=-u(i,j,N(ng)-1,nnew)*ad_DC(i,N(ng))
          ad_u(i,j,N(ng),nnew)=ad_DC(i,N(ng))
          ad_DC(i,N(ng))=0.0_r8
!>        DC(i,1)=tl_u(i,j,1,nnew)-                                     &
!>   &            (tl_BC(i,1)*u(i,j,1,nnew)+                            &
!>   &             tl_FC(i,1)*u(i,j,2,nnew))
!>
          ad_FC(i,1)=-u(i,j,2,nnew)*ad_DC(i,1)
          ad_BC(i,1)=-u(i,j,1,nnew)*ad_DC(i,1)
          ad_u(i,j,1,nnew)=ad_DC(i,1)
          ad_DC(i,1)=0.0_r8
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU,Iend
!>          DC(i,k)=tl_u(i,j,k,nnew)-                                   &
!>   &              (tl_FC(i,k-1)*u(i,j,k-1,nnew)+                      &
!>   &               tl_BC(i,k  )*u(i,j,k  ,nnew)+                      &
!>   &               tl_FC(i,k  )*u(i,j,k+1,nnew))
!>
            ad_FC(i,k-1)=ad_FC(i,k-1)-u(i,j,k-1,nnew)*ad_DC(i,k)
            ad_FC(i,k  )=ad_FC(i,k  )-u(i,j,k+1,nnew)*ad_DC(i,k)
            ad_BC(i,k  )=ad_BC(i,k  )-u(i,j,k  ,nnew)*ad_DC(i,k)
            ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)+ad_DC(i,k)
            ad_DC(i,k)=0.0_r8
          END DO
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
!>          tl_BC(i,k)=tl_Hzk(i,k)-tl_FC(i,k)-tl_FC(i,k-1)
!>
            ad_Hzk(i,k)=ad_Hzk(i,k)+ad_BC(i,k)
            ad_FC(i,k-1)=ad_FC(i,k-1)-ad_BC(i,k)
            ad_FC(i,k  )=ad_FC(i,k  )-ad_BC(i,k)
            ad_BC(i,k)=0.0_r8
          END DO
        END DO
!
!  Compute adjoint off-diagonal coefficients [lambda*dt*Akv/Hz] for
!  the implicit vertical viscosity term at horizontal U-points and
!  vertical W-points.
!
        DO i=IstrU,Iend
!>        tl_FC(i,N(ng))=0.0_r8
!>
          ad_FC(i,N)=0.0_r8
!>        tl_FC(i,0)=0.0_r8
!>
          ad_FC(i,0)=0.0_r8
        END DO
        cff=-lambda*dt(ng)/0.5_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff1=1.0_r8/(z_r(i,j,k+1)+z_r(i-1,j,k+1)-                   &
     &                   z_r(i,j,k  )-z_r(i-1,j,k  ))
!>          tl_FC(i,k)=cff*(tl_cff1*AK(i,k)+cff1*tl_AK(i,k))
!>
            adfac=cff*ad_FC(i,k)
            ad_AK(i,k)=ad_AK(i,k)+cff1*adfac
            ad_cff1=ad_cff1+AK(i,k)*adfac
            ad_FC(i,k)=0.0_r8
!>          tl_cff1=-cff1*cff1*(tl_z_r(i,j,k+1)+tl_z_r(i-1,j,k+1)-      &
!>   &                          tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  ))
!>
            adfac=-cff1*cff1*ad_cff1
            ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-adfac
            ad_z_r(i  ,j,k  )=ad_z_r(i  ,j,k  )-adfac
            ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)+adfac
            ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+adfac
            ad_cff1=0.0_r8
          END DO
        END DO
# endif
!
!  Time step adjoint right-hand-side terms.
!
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff=0.25_r8*dt(ng)
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
          cff=0.25_r8*dt(ng)*3.0_r8/2.0_r8
        ELSE
          cff=0.25_r8*dt(ng)*23.0_r8/12.0_r8
        END IF
        DO i=IstrU,Iend
          DC(i,0)=cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        END DO
!
!  The BASIC STATE "u" used below must be in transport units, but "u"
!  retrived is in m/s so we multiply by "Hzk".
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
# ifdef DIAGNOSTICS_UV
!!          DiaU3wrk(i,j,k,M3rate)=DiaU3wrk(i,j,k,M3rate)*oHz(i,k)
#  ifdef BODYFORCE
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+              &
!!   &                             DC(i,0)*DiaRU(i,j,k,nrhs,M3vvis)*    &
!!   &                             oHz(i,k)
#  endif
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)*oHz(i,k)
#  if defined UV_VIS2 || defined UV_VIS4
!!          DiaU3wrk(i,j,k,M3hvis)=DiaU3wrk(i,j,k,M3hvis)*oHz(i,k)
!!          DiaU3wrk(i,j,k,M3yvis)=DiaU3wrk(i,j,k,M3yvis)*oHz(i,k)
!!          DiaU3wrk(i,j,k,M3xvis)=DiaU3wrk(i,j,k,M3xvis)*oHz(i,k)
#  endif
!!          DO idiag=1,M3pgrd
!!            DiaU3wrk(i,j,k,idiag)=(DiaU3wrk(i,j,k,idiag)+             &
!!   &                               DC(i,0)*DiaRU(i,j,k,nrhs,idiag))*  &
!!   &                              oHz(i,k)
!!          END DO
# endif
# ifdef SPLINES_VVISC
!>          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*oHz(i,k)+                 &
!>   &                       (u(i,j,k,nnew)*Hzk(i,k))*tl_oHz(i,k)
!>
            ad_oHz(i,k)=ad_oHz(i,k)+                                    &
     &                  (u(i,j,k,nnew)*Hzk(i,k))*ad_u(i,j,k,nnew)
            ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*oHz(i,k)
# endif
!>          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+                          &
!>   &                       DC(i,0)*tl_ru(i,j,k,nrhs)
!>
            ad_ru(i,j,k,nrhs)=ad_ru(i,j,k,nrhs)+                        &
     &                        DC(i,0)*ad_u(i,j,k,nnew)
          END DO
        END DO
        DO i=IstrU,Iend
          DO k=1,N(ng)
# if defined SPLINES_VVISC || defined DIAGNOSTICS_UV
            oHz(i,k)=1.0_r8/Hzk(i,k)
!>          tl_oHz(i,k)=-oHz(i,k)*oHz(i,k)*tl_Hzk(i,k)
!>
            ad_Hzk(i,k)=ad_Hzk(i,k)-oHz(i,k)*oHz(i,k)*ad_oHz(i,k)
            ad_oHz(i,k)=0.0_r8
# endif
!>          tl_Hzk(i,k)=0.5_r8*(tl_Hz(i-1,j,k)+                         &
!>   &                          tl_Hz(i  ,j,k))
!>
            adfac=0.5_r8*ad_Hzk(i,k)
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac
            ad_Hzk(i,k)=0.0_r8
!>          tl_AK(i,k)=0.5_r8*(tl_Akv(i-1,j,k)+                         &
!>   &                         tl_Akv(i  ,j,k))
!>
            adfac=0.5_r8*ad_AK(i,k)
            ad_Akv(i-1,j,k)=ad_Akv(i-1,j,k)+adfac
            ad_Akv(i  ,j,k)=ad_Akv(i  ,j,k)+adfac
            ad_AK(i,k)=0.0_r8
          END DO
!>        tl_AK(i,0)=0.5_r8*(tl_Akv(i-1,j,0)+                           &
!>   &                       tl_Akv(i  ,j,0))
!>
          adfac=0.5_r8*ad_AK(i,0)
          ad_Akv(i-1,j,0)=ad_Akv(i-1,j,0)+adfac
          ad_Akv(i  ,j,0)=ad_Akv(i  ,j,0)+adfac
          ad_AK(i,0)=0.0_r8
        END DO
      END DO J_LOOP2
      RETURN
      END SUBROUTINE ad_step3d_uv_tile
#endif
      END MODULE ad_step3d_uv_mod
