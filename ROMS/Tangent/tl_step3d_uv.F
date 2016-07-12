#include "cppdefs.h"
      MODULE tl_step3d_uv_mod
#if defined TANGENT && defined SOLVE3D
!
!svn $Id: tl_step3d_uv.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-steps the tangent linear horizontal momentum      !
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
!
      PRIVATE
      PUBLIC  :: tl_step3d_uv
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_step3d_uv (ng, tile)
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
      CALL wclock_on (ng, iTLM, 34)
# endif
      CALL tl_step3d_uv_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs(ng), nstp(ng), nnew(ng),             &
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
     &                        GRID(ng) % tl_Hz,                         &
     &                        GRID(ng) % z_r,                           &
     &                        GRID(ng) % tl_z_r,                        &
     &                        GRID(ng) % z_w,                           &
     &                        GRID(ng) % tl_z_w,                        &
     &                        MIXING(ng) % Akv,                         &
     &                        MIXING(ng) % tl_Akv,                      &
     &                        COUPLING(ng) % DU_avg1,                   &
     &                        COUPLING(ng) % tl_DU_avg1,                &
     &                        COUPLING(ng) % DV_avg1,                   &
     &                        COUPLING(ng) % tl_DV_avg1,                &
     &                        COUPLING(ng) % DU_avg2,                   &
     &                        COUPLING(ng) % tl_DU_avg2,                &
     &                        COUPLING(ng) % DV_avg2,                   &
     &                        COUPLING(ng) % tl_DV_avg2,                &
     &                        OCEAN(ng) % tl_ru,                        &
     &                        OCEAN(ng) % tl_rv,                        &
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
     &                        OCEAN(ng) % tl_u,                         &
     &                        OCEAN(ng) % v,                            &
     &                        OCEAN(ng) % tl_v,                         &
     &                        OCEAN(ng) % tl_ubar,                      &
     &                        OCEAN(ng) % tl_vbar,                      &
# ifdef NEARSHORE_MELLOR
     &                        OCEAN(ng) % ubar_stokes,                  &
     &                        OCEAN(ng) % tl_ubar_stokes,               &
     &                        OCEAN(ng) % vbar_stokes,                  &
     &                        OCEAN(ng) % tl_vbar_stokes,               &
     &                        OCEAN(ng) % u_stokes,                     &
     &                        OCEAN(ng) % tl_u_stokes,                  &
     &                        OCEAN(ng) % v_stokes,                     &
     &                        OCEAN(ng) % tl_v_stokes,                  &
# endif
     &                        GRID(ng) % Huon,                          &
     &                        GRID(ng) % tl_Huon,                       &
     &                        GRID(ng) % Hvom,                          &
     &                        GRID(ng) % tl_Hvom)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 34)
# endif

      RETURN
      END SUBROUTINE tl_step3d_uv
!
!***********************************************************************
      SUBROUTINE tl_step3d_uv_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nrhs, nstp, nnew,                   &
# ifdef MASKING
     &                              umask, vmask,                       &
# endif
# ifdef WET_DRY_NOT_YET
     &                              umask_wet, vmask_wet,               &
# endif
     &                              om_v, on_u, pm, pn,                 &
     &                              Hz, tl_Hz,                          &
     &                              z_r, tl_z_r,                        &
     &                              z_w, tl_z_w,                        &
     &                              Akv, tl_Akv,                        &
     &                              DU_avg1, tl_DU_avg1,                &
     &                              DV_avg1, tl_DV_avg1,                &
     &                              DU_avg2, tl_DU_avg2,                &
     &                              DV_avg2, tl_DV_avg2,                &
     &                              tl_ru, tl_rv,                       &
# ifdef DIAGNOSTICS_UV
!!   &                              DiaU2wrk, DiaV2wrk,                 &
!!   &                              DiaU2int, DiaV2int,                 &
!!   &                              DiaU3wrk, DiaV3wrk,                 &
!!   &                              DiaRU,    DiaRV,                    &
# endif
     &                              u, tl_u,                            &
     &                              v, tl_v,                            &
     &                              tl_ubar, tl_vbar,                   &
# ifdef NEARSHORE_MELLOR
     &                              ubar_stokes, tl_ubar_stokes,        &
     &                              vbar_stokes, tl_vbar_stokes,        &
     &                              u_stokes, tl_u_stokes,              &
     &                              v_stokes, tl_v_stokes,              &
# endif
     &                              Huon, tl_Huon,                      &
     &                              Hvom, tl_Hvom)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_sources
!
      USE exchange_2d_mod
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
      USE tl_u3dbc_mod, ONLY : tl_u3dbc_tile
      USE tl_v3dbc_mod, ONLY : tl_v3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
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
      real(r8), intent(in) :: tl_ubar_stokes(LBi:,LBj:)
      real(r8), intent(in) :: tl_vbar_stokes(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_Akv(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: tl_DV_avg1(LBi:,LBj:)
      real(r8), intent(in) :: tl_DU_avg2(LBi:,LBj:)
      real(r8), intent(in) :: tl_DV_avg2(LBi:,LBj:)
      real(r8), intent(in) :: tl_ru(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: tl_rv(LBi:,LBj:,0:,:)

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
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_v_stokes(LBi:,LBj:,:)
#  endif
      real(r8), intent(out) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_Huon(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_Hvom(LBi:,LBj:,:)

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
      real(r8), intent(in) :: tl_ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_vbar_stokes(LBi:UBi,LBj:UBj)
#  endif

      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_DU_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_DV_avg2(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(in) :: tl_rv(LBi:UBi,LBj:UBj,0:N(ng),2)
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
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(out) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: tl_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_Hvom(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, idiag, is, j, k

      real(r8) :: cff, cff1, cff2
      real(r8) :: tl_cff, tl_cff1, tl_cff2

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
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_AK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_FC
# ifdef NEARSHORE_MELLOR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_CFs
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_DCs
# endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Hzk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_oHz
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the XI-direction.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          AK(i,0)=0.5_r8*(Akv(i-1,j,0)+                                 &
     &                    Akv(i  ,j,0))
          tl_AK(i,0)=0.5_r8*(tl_Akv(i-1,j,0)+                           &
     &                       tl_Akv(i  ,j,0))
          DO k=1,N(ng)
            AK(i,k)=0.5_r8*(Akv(i-1,j,k)+                               &
     &                      Akv(i  ,j,k))
            tl_AK(i,k)=0.5_r8*(tl_Akv(i-1,j,k)+                         &
     &                         tl_Akv(i  ,j,k))
            Hzk(i,k)=0.5_r8*(Hz(i-1,j,k)+                               &
     &                       Hz(i  ,j,k))
            tl_Hzk(i,k)=0.5_r8*(tl_Hz(i-1,j,k)+                         &
     &                          tl_Hz(i  ,j,k))
# if defined SPLINES_VVISC || defined DIAGNOSTICS_UV
            oHz(i,k)=1.0_r8/Hzk(i,k)
            tl_oHz(i,k)=-oHz(i,k)*oHz(i,k)*tl_Hzk(i,k)
# endif
          END DO
        END DO
!
!  Time step right-hand-side terms.
!
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
        IF (iic(ng).eq.ntfirst(ng).and.SOinitial(ng)) THEN
# else
        IF (iic(ng).eq.ntfirst(ng)) THEN
# endif
          cff=0.25_r8*dt(ng)
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1).and.SOinitial(ng)) THEN
# else
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
# endif
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
!>          u(i,j,k,nnew)=u(i,j,k,nnew)+                                &
!>   &                    DC(i,0)*ru(i,j,k,nrhs)
!>
            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+                          &
     &                       DC(i,0)*tl_ru(i,j,k,nrhs)
# ifdef SPLINES_VVISC
!>          u(i,j,k,nnew)=u(i,j,k,nnew)*oHz(i,k)
!>
            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*oHz(i,k)+                 &
     &                       (u(i,j,k,nnew)*Hzk(i,k))*tl_oHz(i,k)
# endif
# ifdef DIAGNOSTICS_UV
!!          DO idiag=1,M3pgrd
!!            DiaU3wrk(i,j,k,idiag)=(DiaU3wrk(i,j,k,idiag)+             &
!!   &                               DC(i,0)*DiaRU(i,j,k,nrhs,idiag))*  &
!!   &                              oHz(i,k)
!!          END DO
#  if defined UV_VIS2 || defined UV_VIS4
!!          DiaU3wrk(i,j,k,M3xvis)=DiaU3wrk(i,j,k,M3xvis)*oHz(i,k)
!!          DiaU3wrk(i,j,k,M3yvis)=DiaU3wrk(i,j,k,M3yvis)*oHz(i,k)
!!          DiaU3wrk(i,j,k,M3hvis)=DiaU3wrk(i,j,k,M3hvis)*oHz(i,k)
#  endif
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)*oHz(i,k)
#  ifdef BODYFORCE
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+              &
!!   &                             DC(i,0)*DiaRU(i,j,k,nrhs,M3vvis)*    &
!!   &                             oHz(i,k)
#  endif
!!          DiaU3wrk(i,j,k,M3rate)=DiaU3wrk(i,j,k,M3rate)*oHz(i,k)
# endif
          END DO
        END DO

# ifdef SPLINES_VVISC
!
!  Apply spline algorithim to BASIC STATE "u" which should be
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
!  Backward substitution. Save DC for the tangent linearization.
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
!  Use conservative, parabolic spline reconstruction of tangent
!  linear vertical viscosity derivatives.  Then, time step vertical
!  viscosity term implicitly by solving a tridiagonal system.
!
        cff1=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=cff1*Hzk(i,k  )-                                    &
     &              dt(ng)*AK(i,k-1)*oHz(i,k  )
            tl_FC(i,k)=cff1*tl_Hzk(i,k  )-                              &
     &                 dt(ng)*(tl_AK(i,k-1)*oHz(i,k  )+                 &
     &                         AK(i,k-1)*tl_oHz(i,k  ))
            CF(i,k)=cff1*Hzk(i,k+1)-                                    &
     &               dt(ng)*AK(i,k+1)*oHz(i,k+1)
            tl_CF(i,k)=cff1*tl_Hzk(i,k+1)-                              &
     &                 dt(ng)*(tl_AK(i,k+1)*oHz(i,k+1)+                 &
     &                         AK(i,k+1)*tl_oHz(i,k+1))
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,0)=0.0_r8
          tl_CF(i,0)=0.0_r8
          tl_DC(i,0)=0.0_r8
        END DO
!
!  Tangent linear LU decomposition and forward substitution.
!
        cff1=1.0_r8/3.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                         &
     &              dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
            tl_BC(i,k)=cff1*(tl_Hzk(i,k)+tl_Hzk(i,k+1))+                &
     &                 dt(ng)*(tl_AK(i,k)*(oHz(i,k)+oHz(i,k+1))+        &
     &                         AK(i,k)*(tl_oHz(i,k)+tl_oHz(i,k+1)))
            cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
            CF(i,k)=cff*CF(i,k)
            tl_DC(i,k)=cff*(tl_u(i,j,k+1,nnew)-                         &
     &                      tl_u(i,j,k  ,nnew)-                         &
     &                      (tl_FC(i,k)*DC(i,k-1)+                      &
     &                       tl_BC(i,k)*DC(i,k  )+                      &
     &                       tl_CF(i,k)*DC(i,k+1))-                     &
     &                      FC(i,k)*tl_DC(i,k-1))
          END DO
        END DO
!
!  Tangent linear backward substitution.
!
        DO i=IstrU,Iend
          tl_DC(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
          END DO
        END DO
!
!  Compute tl_DC before multiplying BASIC STATE spline gradients
!  DC by AK.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            tl_DC(i,k)=tl_DC(i,k)*AK(i,k)+DC(i,k)*tl_AK(i,k)
            DC(i,k)=DC(i,k)*AK(i,k)
!>          cff=dt(ng)*oHz(i,k)*(DC(i,k)-DC(i,k-1))
!>
            tl_cff=dt(ng)*(tl_oHz(i,k)*(DC(i,k)-DC(i,k-1))+             &
     &                     oHz(i,k)*(tl_DC(i,k)-tl_DC(i,k-1)))
!>          u(i,j,k,nnew)=u(i,j,k,nnew)+cff
!>
            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+tl_cff
#  ifdef DIAGNOSTICS_UV
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+cff
#  endif
          END DO
        END DO
# else
!
!  Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal U-points and
!  vertical W-points.
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
            tl_cff1=-cff1*cff1*(tl_z_r(i,j,k+1)+tl_z_r(i-1,j,k+1)-      &
     &                          tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  ))
            FC(i,k)=cff*cff1*AK(i,k)
            tl_FC(i,k)=cff*(tl_cff1*AK(i,k)+cff1*tl_AK(i,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          tl_FC(i,0)=0.0_r8
          FC(i,N(ng))=0.0_r8
          tl_FC(i,N(ng))=0.0_r8
        END DO
!
!  Solve the tangent linear tridiagonal system.
!  (DC is a tangent linear variable here).
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
            BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
            tl_BC(i,k)=tl_Hzk(i,k)-tl_FC(i,k)-tl_FC(i,k-1)
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU,Iend
            DC(i,k)=tl_u(i,j,k,nnew)-                                   &
     &              (tl_FC(i,k-1)*u(i,j,k-1,nnew)+                      &
     &               tl_BC(i,k  )*u(i,j,k  ,nnew)+                      &
     &               tl_FC(i,k  )*u(i,j,k+1,nnew))
          END DO
        END DO
        DO i=IstrU,Iend
          DC(i,1)=tl_u(i,j,1,nnew)-                                     &
     &            (tl_BC(i,1)*u(i,j,1,nnew)+                            &
     &             tl_FC(i,1)*u(i,j,2,nnew))
          DC(i,N(ng))=tl_u(i,j,N(ng),nnew)-                             &
     &                (tl_FC(i,N(ng)-1)*u(i,j,N(ng)-1,nnew)+            &
     &                 tl_BC(i,N(ng)  )*u(i,j,N(ng)  ,nnew))
        END DO
        DO i=IstrU,Iend
          cff=1.0_r8/BC(i,1)
          CF(i,1)=cff*FC(i,1)
          DC(i,1)=cff*DC(i,1)
        END DO
        DO k=2,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
            CF(i,k)=cff*FC(i,k)
            DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
          END DO
        END DO
!
!  Compute new solution by back substitution.
!  (DC is a tangent linear variable here).
!
        DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!        wrk(i,N(ng))=u(i,j,N(ng),nnew)*oHz(i,N(ng))
#  endif
          DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/        &
     &                (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
!>        u(i,j,N(ng),nnew)=DC(i,N(ng))
!>
          tl_u(i,j,N(ng),nnew)=DC(i,N(ng))
#  ifdef DIAGNOSTICS_UV
!!        DiaU3wrk(i,j,N(ng),M3vvis)=DiaU3wrk(i,j,N(ng),M3vvis)+        &
!!   &                               u(i,j,N(ng),nnew)-wrk(i,N(ng))
#  endif
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          wrk(i,k)=u(i,j,k,nnew)*oHz(i,k)
#  endif
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>          u(i,j,k,nnew)=DC(i,k)
!>
            tl_u(i,j,k,nnew)=DC(i,k)
#  ifdef DIAGNOSTICS_UV
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+              &
!!   &                             u(i,j,k,nnew)-wrk(i,k)
#  endif
          END DO
        END DO
# endif
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, ubar=DU_avg1/(D*on_u). Recall that, D=CF(:,0).
!
        DO i=IstrU,Iend
          CF(i,0)=Hzk(i,1)
          tl_CF(i,0)=tl_Hzk(i,1)
          DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)
          tl_DC(i,0)=tl_u(i,j,1,nnew)*Hzk(i,1)+                         &
     &               u(i,j,1,nnew)*tl_Hzk(i,1)
# ifdef NEARSHORE_MELLOR
          DCs(i,0)=u_stokes(i,j,1)*Hzk(i,1)
          tl_DCs(i,0)=tl_u_stokes(i,j,1)*Hzk(i,1)+                      &
     &                u_stokes(i,j,1)*tl_Hzk(i,1)
# endif
# ifdef DIAGNOSTICS_UV
!!        Dwrk(i,M2pgrd)=DiaU3wrk(i,j,1,M3pgrd)*Hzk(i,1)
!!        Dwrk(i,M2bstr)=DiaU3wrk(i,j,1,M3vvis)*Hzk(i,1)
#  ifdef UV_COR
!!        Dwrk(i,M2fcor)=DiaU3wrk(i,j,1,M3fcor)*Hzk(i,1)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!        Dwrk(i,M2xvis)=DiaU3wrk(i,j,1,M3xvis)*Hzk(i,1)
!!        Dwrk(i,M2yvis)=DiaU3wrk(i,j,1,M3yvis)*Hzk(i,1)
!!        Dwrk(i,M2hvis)=DiaU3wrk(i,j,1,M3hvis)*Hzk(i,1)
#  endif
#  ifdef UV_ADV
!!        Dwrk(i,M2xadv)=DiaU3wrk(i,j,1,M3xadv)*Hzk(i,1)
!!        Dwrk(i,M2yadv)=DiaU3wrk(i,j,1,M3yadv)*Hzk(i,1)
!!        Dwrk(i,M2hadv)=DiaU3wrk(i,j,1,M3hadv)*Hzk(i,1)
#  endif
#  ifdef NEARSHORE_MELLOR
!!        Dwrk(i,M2hrad)=DiaU3wrk(i,j,1,M3hrad)*Hzk(i,1)
#  endif
# endif
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            CF(i,0)=CF(i,0)+Hzk(i,k)
            tl_CF(i,0)=tl_CF(i,0)+tl_Hzk(i,k)
            DC(i,0)=DC(i,0)+u(i,j,k,nnew)*Hzk(i,k)
            tl_DC(i,0)=tl_DC(i,0)+                                      &
     &                 tl_u(i,j,k,nnew)*Hzk(i,k)+                       &
     &                 u(i,j,k,nnew)*tl_Hzk(i,k)
# ifdef NEARSHORE_MELLOR
            DCs(i,0)=DCs(i,0)+u_stokes(i,j,k)*Hzk(i,k)
            tl_DCs(i,0)=tl_DCs(i,0)+                                    &
     &                  tl_u_stokes(i,j,k)*Hzk(i,k)+                    &
     &                  u_stokes(i,j,k)*tl_Hzk(i,k)
# endif
# ifdef DIAGNOSTICS_UV
!!          Dwrk(i,M2pgrd)=Dwrk(i,M2pgrd)+                              &
!!   &                     DiaU3wrk(i,j,k,M3pgrd)*Hzk(i,k)
!!          Dwrk(i,M2bstr)=Dwrk(i,M2bstr)+                              &
!!   &                     DiaU3wrk(i,j,k,M3vvis)*Hzk(i,k)
#  ifdef UV_COR
!!          Dwrk(i,M2fcor)=Dwrk(i,M2fcor)+                              &
!!   &                     DiaU3wrk(i,j,k,M3fcor)*Hzk(i,k)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!          Dwrk(i,M2xvis)=Dwrk(i,M2xvis)+                              &
!!   &                     DiaU3wrk(i,j,k,M3xvis)*Hzk(i,k)
!!          Dwrk(i,M2yvis)=Dwrk(i,M2yvis)+                              &
!!   &                     DiaU3wrk(i,j,k,M3yvis)*Hzk(i,k)
!!          Dwrk(i,M2hvis)=Dwrk(i,M2hvis)+                              &
!!   &                     DiaU3wrk(i,j,k,M3hvis)*Hzk(i,k)
#  endif
#  ifdef UV_ADV
!!          Dwrk(i,M2xadv)=Dwrk(i,M2xadv)+                              &
!!   &                     DiaU3wrk(i,j,k,M3xadv)*Hzk(i,k)
!!          Dwrk(i,M2yadv)=Dwrk(i,M2yadv)+                              &
!!   &                     DiaU3wrk(i,j,k,M3yadv)*Hzk(i,k)
!!          Dwrk(i,M2hadv)=Dwrk(i,M2hadv)+                              &
!!   &                     DiaU3wrk(i,j,k,M3hadv)*Hzk(i,k)
#  endif
#  ifdef NEARSHORE_MELLOR
!!          Dwrk(i,M2hrad)=Dwrk(i,M2hrad)+                              &
!!   &                     DiaU3wrk(i,j,k,M3hrad)*Hzk(i,k)
#  endif
# endif
          END DO
        END DO
        DO i=IstrU,Iend
          DC1(i,0)=DC(i,0)                                ! intermediate
          cff1=1.0_r8/(CF(i,0)*on_u(i,j))
          tl_cff1=-cff1*cff1*tl_CF(i,0)*on_u(i,j)
          DC(i,0)=(DC(i,0)*on_u(i,j)-DU_avg1(i,j))*cff1      ! recursive
          tl_DC(i,0)=(tl_DC(i,0)*on_u(i,j)-tl_DU_avg1(i,j))*cff1+       &
     &               (DC1(i,0)*on_u(i,j)-DU_avg1(i,j))*tl_cff1
# ifdef NEARSHORE_MELLOR
          DCs1(i)=DCs(i,0)                                ! intermediate
          cff2=1.0_r8/CF(i,0)
          tl_cff2=-cff2*cff2*tl_CF(i,0)
          DCs(i,0)=DCs(i,0)*cff2-ubar_stokes(i,j)            ! recursive
          tl_DCs(i,0)=tl_DCs(i,0)*cff2+                                 &
     &                DCs1(i)*tl_cff2-tl_ubar_stokes(i,j)
# endif
# ifdef DIAGNOSTICS_UV
!!        DO idiag=1,M2pgrd
!!          Dwrk(i,idiag)=(Dwrk(i,idiag)*on_u(i,j)-                     &
!!   &                     DiaU2wrk(i,j,idiag))*cff1
!!        END DO
!!        Dwrk(i,M2bstr)=(Dwrk(i,M2bstr)*on_u(i,j)-                     &
!!   &                    DiaU2wrk(i,j,M2bstr)-DiaU2wrk(i,j,M2sstr))*   &
!!   &                   cff1
# endif
        END DO
!
!  Couple and update new solution.
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
!>          u(i,j,k,nnew)=u(i,j,k,nnew)-DC(i,0)
!>
            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_DC(i,0)
# ifdef MASKING
!>          u(i,j,k,nnew)=u(i,j,k,nnew)*umask(i,j)
!>
            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*umask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>          u(i,j,k,nnew)=u(i,j,k,nnew)*umask_wet(i,j)
!>
            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*umask_wet(i,j)
# endif
# ifdef NEARSHORE_MELLOR
!>          u_stokes(i,j,k)=u_stokes(i,j,k)-DCs(i,0)
!>
            tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)-tl_DCs(i,0)
#  ifdef MASKING
!>          u_stokes(i,j,k)=u_stokes(i,j,k)*umask(i,j)
!>
            tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*umask(i,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>          u_stokes(i,j,k)=u_stokes(i,j,k)*umask_wet(i,j)
!>
            tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*umask_wet(i,j)
#  endif
# endif
# ifdef DIAGNOSTICS_UV
!!          DiaU3wrk(i,j,k,M3pgrd)=DiaU3wrk(i,j,k,M3pgrd)-              &
!!   &                             Dwrk(i,M2pgrd)
!!          DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)-              &
!!   &                             Dwrk(i,M2bstr)
#  ifdef UV_COR
!!          DiaU3wrk(i,j,k,M3fcor)=DiaU3wrk(i,j,k,M3fcor)-              &
!!   &                             Dwrk(i,M2fcor)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!          DiaU3wrk(i,j,k,M3xvis)=DiaU3wrk(i,j,k,M3xvis)-              &
!!   &                             Dwrk(i,M2xvis)
!!          DiaU3wrk(i,j,k,M3yvis)=DiaU3wrk(i,j,k,M3yvis)-              &
!!   &                             Dwrk(i,M2yvis)
!!          DiaU3wrk(i,j,k,M3hvis)=DiaU3wrk(i,j,k,M3hvis)-              &
!!   &                             Dwrk(i,M2hvis)
#  endif
#  ifdef UV_ADV
!!          DiaU3wrk(i,j,k,M3xadv)=DiaU3wrk(i,j,k,M3xadv)-              &
!!   &                             Dwrk(i,M2xadv)
!!          DiaU3wrk(i,j,k,M3yadv)=DiaU3wrk(i,j,k,M3yadv)-              &
!!   &                             Dwrk(i,M2yadv)
!!          DiaU3wrk(i,j,k,M3hadv)=DiaU3wrk(i,j,k,M3hadv)-              &
!!   &                             Dwrk(i,M2hadv)
#  endif
#  ifdef NEARSHORE_MELLOR
!!          DiaU3wrk(i,j,k,M3hrad)=DiaU3wrk(i,j,k,M3hrad)-              &
!!   &                             Dwrk(i,M2hrad)
#  endif
# endif
          END DO
        END DO

# if defined DIAGNOSTICS_UV && defined MASKING
!!      DO k=1,N(ng)
!!        DO i=IstrU,Iend
!!          DO idiag=1,NDM3d
!!            DiaU3wrk(i,j,k,idiag)=DiaU3wrk(i,j,k,idiag)*umask(i,j)
!!          END DO
!!        END DO
!!      END DO
# endif
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the ETA-direction.
!-----------------------------------------------------------------------
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            AK(i,0)=0.5_r8*(Akv(i,j-1,0)+                               &
     &                      Akv(i,j  ,0))
            tl_AK(i,0)=0.5_r8*(tl_Akv(i,j-1,0)+                         &
     &                         tl_Akv(i,j  ,0))
            DO k=1,N(ng)
              AK(i,k)=0.5_r8*(Akv(i,j-1,k)+                             &
     &                        Akv(i,j  ,k))
              tl_AK(i,k)=0.5_r8*(tl_Akv(i,j-1,k)+                       &
     &                           tl_Akv(i,j  ,k))
              Hzk(i,k)=0.5_r8*(Hz(i,j-1,k)+                             &
     &                         Hz(i,j  ,k))
              tl_Hzk(i,k)=0.5_r8*(tl_Hz(i,j-1,k)+                       &
     &                            tl_Hz(i,j  ,k))
# if defined SPLINES_VVISC || defined DIAGNOSTICS_UV
              oHz(i,k)=1.0_r8/Hzk(i,k)
              tl_oHz(i,k)=-oHz(i,k)*oHz(i,k)*tl_Hzk(i,k)
# endif
            END DO
          END DO
!
!  Time step right-hand-side terms.
!
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
          IF (iic(ng).eq.ntfirst(ng).and.SOinitial(ng)) THEN
# else
          IF (iic(ng).eq.ntfirst(ng)) THEN
# endif
            cff=0.25_r8*dt(ng)
# if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1).and.SOinitial(ng)) THEN
# else
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
# endif
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
!>            v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k,nrhs)
!>
              tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+                        &
     &                         DC(i,0)*tl_rv(i,j,k,nrhs)
# ifdef SPLINES_VVISC
!>            v(i,j,k,nnew)=v(i,j,k,nnew)*oHz(i,k)
!>
              tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*oHz(i,k)+               &
     &                         (v(i,j,k,nnew)*Hzk(i,k))*tl_oHz(i,k)
# endif
# ifdef DIAGNOSTICS_UV
!!            DO idiag=1,M3pgrd
!!              DiaV3wrk(i,j,k,idiag)=(DiaV3wrk(i,j,k,idiag)+           &
!!   &                                 DC(i,0)*                         &
!!   &                                 DiaRV(i,j,k,nrhs,idiag))*        &
!!   &                                oHz(i,k)
!!            END DO
#  if defined UV_VIS2 || defined UV_VIS4
!!            DiaV3wrk(i,j,k,M3xvis)=DiaV3wrk(i,j,k,M3xvis)*oHz(i,k)
!!            DiaV3wrk(i,j,k,M3yvis)=DiaV3wrk(i,j,k,M3yvis)*oHz(i,k)
!!            DiaV3wrk(i,j,k,M3hvis)=DiaV3wrk(i,j,k,M3hvis)*oHz(i,k)
#  endif
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)*oHz(i,k)
#  ifdef BODYFORCE
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+            &
!!   &                               DC(i,0)*DiaRV(i,j,k,nrhs,M3vvis)*  &
!!   &                               oHz(i,k)
#  endif
!!            DiaV3wrk(i,j,k,M3rate)=DiaV3wrk(i,j,k,M3rate)*oHz(i,k)
# endif
            END DO
          END DO

# ifdef SPLINES_VVISC
!
!  Apply spline algorithim to BASIC STATE "v" which should be
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
!  Backward substitution. Save DC for the tangent linearization.
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
!  Use conservative, parabolic spline reconstruction of tangent
!  linear vertical viscosity derivatives.  Then, time step vertical
!  viscosity term implicitly by solving a tridiagonal system.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hzk(i,k  )-                                  &
     &                dt(ng)*AK(i,k-1)*oHz(i,k  )
              tl_FC(i,k)=cff1*tl_Hzk(i,k  )-                            &
     &                   dt(ng)*(tl_AK(i,k-1)*oHz(i,k  )+               &
     &                           AK(i,k-1)*tl_oHz(i,k  ))
              CF(i,k)=cff1*Hzk(i,k+1)-                                  &
     &                dt(ng)*AK(i,k+1)*oHz(i,k+1)
              tl_CF(i,k)=cff1*tl_Hzk(i,k+1)-                            &
     &                   dt(ng)*(tl_AK(i,k+1)*oHz(i,k+1)+               &
     &                           AK(i,k+1)*tl_oHz(i,k+1))
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            tl_CF(i,0)=0.0_r8
            tl_DC(i,0)=0.0_r8
          END DO
!
!  Tangent linear LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                dt(ng)*AK(i,k)*(oHz(i,k)+oHz(i,k+1))
              tl_BC(i,k)=cff1*(tl_Hzk(i,k)+tl_Hzk(i,k+1))+              &
     &                   dt(ng)*(tl_AK(i,k)*(oHz(i,k)+oHz(i,k+1))+      &
     &                           AK(i,k)*(tl_oHz(i,k)+tl_oHz(i,k+1)))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              tl_DC(i,k)=cff*(tl_v(i,j,k+1,nnew)-                       &
     &                        tl_v(i,j,k  ,nnew)-                       &
     &                        (tl_FC(i,k)*DC(i,k-1)+                    &
     &                         tl_BC(i,k)*DC(i,k  )+                    &
     &                         tl_CF(i,k)*DC(i,k+1))-                   &
     &                        FC(i,k)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Tangent linear backward substitution.
!
          DO i=Istr,Iend
            tl_DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
            END DO
          END DO
!
!  Compute tl_DC before multiplying BASIC STATE spline gradients
!  DC by AK.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              tl_DC(i,k)=tl_DC(i,k)*AK(i,k)+DC(i,k)*tl_AK(i,k)
              DC(i,k)=DC(i,k)*AK(i,k)
!>            cff=dt(ng)*oHz(i,k)*(DC(i,k)-DC(i,k-1))
!>
              tl_cff=dt(ng)*(tl_oHz(i,k)*(DC(i,k)-DC(i,k-1))+           &
     &                       oHz(i,k)*(tl_DC(i,k)-tl_DC(i,k-1)))
!>            v(i,j,k,nnew)=v(i,j,k,nnew)+cff
!>
              tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+tl_cff
#  ifdef DIAGNOSTICS_UV
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+cff
#  endif
            END DO
          END DO
# else
!
!  Compute off-diagonal coefficients [lambda*dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal V-points and
!  vertical W-points.
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
              tl_cff1=-cff1*cff1*(tl_z_r(i,j,k+1)+tl_z_r(i,j-1,k+1)-    &
     &                            tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  ))
              FC(i,k)=cff*cff1*AK(i,k)
              tl_FC(i,k)=cff*(tl_cff1*AK(i,k)+cff1*tl_AK(i,k))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            tl_FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
            tl_FC(i,N(ng))=0.0_r8
          END DO
!
!  Solve the tangent linear tridiagonal system.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
              tl_BC(i,k)=tl_Hzk(i,k)-tl_FC(i,k)-tl_FC(i,k-1)
            END DO
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              DC(i,k)=tl_v(i,j,k,nnew)-                                 &
     &                (tl_FC(i,k-1)*v(i,j,k-1,nnew)+                    &
     &                 tl_BC(i,k  )*v(i,j,k  ,nnew)+                    &
     &                 tl_FC(i,k  )*v(i,j,k+1,nnew))
            END DO
          END DO
          DO i=Istr,Iend
            DC(i,1)=tl_v(i,j,1,nnew)-                                   &
     &              (tl_BC(i,1)*v(i,j,1,nnew)+                          &
     &               tl_FC(i,1)*v(i,j,2,nnew))
            DC(i,N(ng))=tl_v(i,j,N(ng),nnew)-                           &
     &                  (tl_FC(i,N(ng)-1)*v(i,j,N(ng)-1,nnew)+          &
     &                   tl_BC(i,N(ng)  )*v(i,j,N(ng)  ,nnew))
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!  (DC is a tangent linear variable here).
!
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!          wrk(i,N(ng))=v(i,j,N(ng),nnew)*oHz(i,N(ng))
#  endif
            DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/      &
     &                  (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
!>          v(i,j,N(ng),nnew)=DC(i,N(ng))
!>
            tl_v(i,j,N(ng),nnew)=DC(i,N(ng))
#  ifdef DIAGNOSTICS_UV
!!          DiaV3wrk(i,j,N(ng),M3vvis)=DiaV3wrk(i,j,N(ng),M3vvis)+      &
!!   &                                 v(i,j,N(ng),nnew)-wrk(i,N(ng))
#  endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!            wrk(i,k)=v(i,j,k,nnew)*oHz(i,k)
#  endif
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!>            v(i,j,k,nnew)=DC(i,k)
!>
              tl_v(i,j,k,nnew)=DC(i,k)
#  ifdef DIAGNOSTICS_UV
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+            &
!!   &                               v(i,j,k,nnew)-wrk(i,k)
#  endif
            END DO
          END DO
# endif
!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, vbar=DV_avg1/(D*om_v). Recall that, D=CF(:,0).
!
          DO i=Istr,Iend
            CF(i,0)=Hzk(i,1)
            tl_CF(i,0)=tl_Hzk(i,1)
            DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)
            tl_DC(i,0)=tl_v(i,j,1,nnew)*Hzk(i,1)+                       &
     &                 v(i,j,1,nnew)*tl_Hzk(i,1)
# ifdef NEARSHORE_MELLOR
            DCs(i,0)=v_stokes(i,j,1)*Hzk(i,1)
            tl_DCs(i,0)=tl_v_stokes(i,j,1)*Hzk(i,1)+                    &
     &                  v_stokes(i,j,1)*tl_Hzk(i,1)
# endif
# ifdef DIAGNOSTICS_UV
!!          Dwrk(i,M2pgrd)=DiaV3wrk(i,j,1,M3pgrd)*Hzk(i,1)
!!          Dwrk(i,M2bstr)=DiaV3wrk(i,j,1,M3vvis)*Hzk(i,1)
#  ifdef UV_COR
!!          Dwrk(i,M2fcor)=DiaV3wrk(i,j,1,M3fcor)*Hzk(i,1)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!          Dwrk(i,M2xvis)=DiaV3wrk(i,j,1,M3xvis)*Hzk(i,1)
!!          Dwrk(i,M2yvis)=DiaV3wrk(i,j,1,M3yvis)*Hzk(i,1)
!!          Dwrk(i,M2hvis)=DiaV3wrk(i,j,1,M3hvis)*Hzk(i,1)
#  endif
#  ifdef UV_ADV
!!          Dwrk(i,M2xadv)=DiaV3wrk(i,j,1,M3xadv)*Hzk(i,1)
!!          Dwrk(i,M2yadv)=DiaV3wrk(i,j,1,M3yadv)*Hzk(i,1)
!!          Dwrk(i,M2hadv)=DiaV3wrk(i,j,1,M3hadv)*Hzk(i,1)
#  endif
#  ifdef NEARSHORE_MELLOR
!!          Dwrk(i,M2hrad)=DiaV3wrk(i,j,1,M3hrad)*Hzk(i,1)
#  endif
# endif
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              CF(i,0)=CF(i,0)+Hzk(i,k)
              tl_CF(i,0)=tl_CF(i,0)+tl_Hzk(i,k)
              DC(i,0)=DC(i,0)+v(i,j,k,nnew)*Hzk(i,k)
              tl_DC(i,0)=tl_DC(i,0)+                                    &
     &                   tl_v(i,j,k,nnew)*Hzk(i,k)+                     &
     &                   v(i,j,k,nnew)*tl_Hzk(i,k)
# ifdef NEARSHORE_MELLOR
              DCs(i,0)=DCs(i,0)+v_stokes(i,j,k)*Hzk(i,k)
              tl_DCs(i,0)=tl_DCs(i,0)+                                  &
     &                    tl_v_stokes(i,j,k)*Hzk(i,k)+                  &
     &                    v_stokes(i,j,k)*tl_Hzk(i,k)
# endif
# ifdef DIAGNOSTICS_UV
!!            Dwrk(i,M2pgrd)=Dwrk(i,M2pgrd)+                            &
!!   &                       DiaV3wrk(i,j,k,M3pgrd)*Hzk(i,k)
!!            Dwrk(i,M2bstr)=Dwrk(i,M2bstr)+                            &
!!   &                       DiaV3wrk(i,j,k,M3vvis)*Hzk(i,k)
#  ifdef UV_COR
!!            Dwrk(i,M2fcor)=Dwrk(i,M2fcor)+                            &
!!   &                       DiaV3wrk(i,j,k,M3fcor)*Hzk(i,k)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!            Dwrk(i,M2xvis)=Dwrk(i,M2xvis)+                            &
!!   &                       DiaV3wrk(i,j,k,M3xvis)*Hzk(i,k)
!!            Dwrk(i,M2yvis)=Dwrk(i,M2yvis)+                            &
!!   &                       DiaV3wrk(i,j,k,M3yvis)*Hzk(i,k)
!!            Dwrk(i,M2hvis)=Dwrk(i,M2hvis)+                            &
!!   &                       DiaV3wrk(i,j,k,M3hvis)*Hzk(i,k)
#  endif
#  ifdef UV_ADV
!!            Dwrk(i,M2xadv)=Dwrk(i,M2xadv)+                            &
!!   &                       DiaV3wrk(i,j,k,M3xadv)*Hzk(i,k)
!!            Dwrk(i,M2yadv)=Dwrk(i,M2yadv)+                            &
!!   &                       DiaV3wrk(i,j,k,M3yadv)*Hzk(i,k)
!!            Dwrk(i,M2hadv)=Dwrk(i,M2hadv)+                            &
!!   &                       DiaV3wrk(i,j,k,M3hadv)*Hzk(i,k)
#  endif
#  ifdef NEARSHORE_MELLOR
!!            Dwrk(i,M2hrad)=Dwrk(i,M2hrad)+                            &
!!   &                       DiaV3wrk(i,j,k,M3hrad)*Hzk(i,k)
#  endif
# endif
            END DO
          END DO
          DO i=Istr,Iend
            DC1(i,0)=DC(i,0)                              ! intermediate
            cff1=1.0_r8/(CF(i,0)*om_v(i,j))
            tl_cff1=-cff1*cff1*tl_CF(i,0)*om_v(i,j)
            DC(i,0)=(DC(i,0)*om_v(i,j)-DV_avg1(i,j))*cff1    ! recursive
            tl_DC(i,0)=(tl_DC(i,0)*om_v(i,j)-tl_DV_avg1(i,j))*cff1+     &
     &                 (DC1(i,0)*om_v(i,j)-DV_avg1(i,j))*tl_cff1
# ifdef NEARSHORE_MELLOR
            DCs1(i)=DCs(i,0)                              ! intermediate
            cff2=1.0_r8/CF(i,0)
            tl_cff2=-cff2*cff2*tl_CF(i,0)
            DCs(i,0)=DCs(i,0)*cff2-vbar_stokes(i,j)          ! recursive
            tl_DCs(i,0)=tl_DCs(i,0)*cff2+                               &
     &                  DCs1(i,0)*tl_cff2-tl_vbar_stokes(i,j)
# endif
# ifdef DIAGNOSTICS_UV
!!          DO idiag=1,M2pgrd
!!            Dwrk(i,idiag)=(Dwrk(i,idiag)*om_v(i,j)-                   &
!!   &                       DiaV2wrk(i,j,idiag))*cff1
!!          END DO
!!          Dwrk(i,M2bstr)=(Dwrk(i,M2bstr)*om_v(i,j)-                   &
!!   &                      DiaV2wrk(i,j,M2bstr)-DiaV2wrk(i,j,M2sstr))* &
!!   &                     cff1
# endif
          END DO
!
!  Couple and update new solution.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            v(i,j,k,nnew)=v(i,j,k,nnew)-DC(i,0)
!>
              tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_DC(i,0)
# ifdef MASKING
!>            v(i,j,k,nnew)=v(i,j,k,nnew)*vmask(i,j)
!>
              tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*vmask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>            v(i,j,k,nnew)=v(i,j,k,nnew)*vmask_wet(i,j)
!>
              tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*vmask_wet(i,j)
# endif
# ifdef NEARSHORE_MELLOR
!>            v_stokes(i,j,k)=v_stokes(i,j,k)-DCs(i,0)
!>
              tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)-tl_DCs(i,0)
#  ifdef MASKING
!>            v_stokes(i,j,k)=v_stokes(i,j,k)*vmask(i,j)
!>
              tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*vmask(i,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>            v_stokes(i,j,k)=v_stokes(i,j,k)*vmask_wet(i,j)
!>
              tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*vmask_wet(i,j)
#  endif
# endif
# ifdef DIAGNOSTICS_UV
!!            DiaV3wrk(i,j,k,M3pgrd)=DiaV3wrk(i,j,k,M3pgrd)-            &
!!   &                               Dwrk(i,M2pgrd)
!!            DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)-            &
!!   &                               Dwrk(i,M2bstr)
#  ifdef UV_COR
!!            DiaV3wrk(i,j,k,M3fcor)=DiaV3wrk(i,j,k,M3fcor)-            &
!!   &                               Dwrk(i,M2fcor)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!            DiaV3wrk(i,j,k,M3xvis)=DiaV3wrk(i,j,k,M3xvis)-            &
!!   &                               Dwrk(i,M2xvis)
!!            DiaV3wrk(i,j,k,M3yvis)=DiaV3wrk(i,j,k,M3yvis)-            &
!!   &                               Dwrk(i,M2yvis)
!!            DiaV3wrk(i,j,k,M3hvis)=DiaV3wrk(i,j,k,M3hvis)-            &
!!   &                               Dwrk(i,M2hvis)
#  endif
#  ifdef UV_ADV
!!            DiaV3wrk(i,j,k,M3xadv)=DiaV3wrk(i,j,k,M3xadv)-            &
!!   &                               Dwrk(i,M2xadv)
!!            DiaV3wrk(i,j,k,M3yadv)=DiaV3wrk(i,j,k,M3yadv)-            &
!!   &                               Dwrk(i,M2yadv)
!!            DiaV3wrk(i,j,k,M3hadv)=DiaV3wrk(i,j,k,M3hadv)-            &
!!   &                               Dwrk(i,M2hadv)
#  endif
#  ifdef NEARSHORE_MELLOR
!!            DiaV3wrk(i,j,k,M3hrad)=DiaV3wrk(i,j,k,M3hrad)-            &
!!   &                               Dwrk(i,M2hrad)
#  endif
# endif
            END DO
          END DO

# if defined DIAGNOSTICS_UV && defined MASKING
!!        DO k=1,N(ng)
!!          DO i=Istr,Iend
!!            DO idiag=1,NDM3d
!!              DiaV3wrk(i,j,k,idiag)=DiaV3wrk(i,j,k,idiag)*vmask(i,j)
!!            END DO
!!          END DO
!!        END DO
# endif
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
!>    CALL u3dbc_tile (ng, tile,                                        &
!>   &                 LBi, UBi, LBj, UBj, N(ng),                       &
!>   &                 IminS, ImaxS, JminS, JmaxS,                      &
!>   &                 nstp, nnew,                                      &
!>   &                 u)
!>
      CALL tl_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nnew,                                   &
     &                    tl_u)
!>    CALL v3dbc_tile (ng, tile,                                        &
!>   &                 LBi, UBi, LBj, UBj, N(ng),                       &
!>   &                 IminS, ImaxS, JminS, JmaxS,                      &
!>   &                 nstp, nnew,                                      &
!>   &                 v)
!>
      CALL tl_v3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nnew,                                   &
     &                    tl_v)
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
              DO k=1,N(ng)
                cff1=1.0_r8/(on_u(i,j)*                                 &
     &                       0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+       &
     &                               z_w(i  ,j,k)-z_w(i  ,j,k-1)))
                tl_cff1=-cff1*cff1*on_u(i,j)*                           &
     &                  0.5_r8*(tl_z_w(i-1,j,k)-tl_z_w(i-1,j,k-1)+      &
     &                          tl_z_w(i  ,j,k)-tl_z_w(i  ,j,k-1))
!>              u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
!>
                tl_u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*tl_cff1
              END DO
            ELSE
              DO k=1,N(ng)
                cff1=1.0_r8/(om_v(i,j)*                                 &
     &                       0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+       &
     &                               z_w(i,j  ,k)-z_w(i,j  ,k-1)))
                tl_cff1=-cff1*cff1*om_v(i,j)*                           &
     &                  0.5_r8*(tl_z_w(i,j-1,k)-tl_z_w(i,j-1,k-1)+      &
     &                          tl_z_w(i,j  ,k)-tl_z_w(i,j  ,k-1))
!>              v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
!>
                tl_v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*tl_cff1
              END DO
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Couple 2D and 3D momentum equations.
!-----------------------------------------------------------------------
!
!  Couple velocity component in the XI-direction.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          DC(i,0)=0.0_r8
          tl_DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
          tl_CF(i,0)=0.0_r8
# ifdef NEARSHORE_MELLOR
          CFs(i,0)=0.0_r8
          tl_CFs(i,0)=0.0_r8
# endif
          FC(i,0)=0.0_r8
          tl_FC(i,0)=0.0_r8
        END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
        DO k=1,N(ng)
          DO i=IstrP,IendT
            cff=0.5_r8*on_u(i,j)
            DC(i,k)=cff*(Hz(i,j,k)+Hz(i-1,j,k))
            tl_DC(i,k)=cff*(tl_Hz(i,j,k)+tl_Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
            tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
            CF(i,0)=CF(i,0)+                                            &
     &              DC(i,k)*u(i,j,k,nnew)
            tl_CF(i,0)=tl_CF(i,0)+                                      &
     &                 tl_DC(i,k)*u(i,j,k,nnew)+                        &
     &                 DC(i,k)*tl_u(i,j,k,nnew)
# ifdef NEARSHORE_MELLOR
            CFs(i,0)=CFs(i,0)+                                          &
     &               DC(i,k)*u_stokes(i,j,k)
            tl_CFs(i,0)=tl_CFs(i,0)+                                    &
     &                  tl_DC(i,k)*u_stokes(i,j,k)+                     &
     &                  DC(i,k)*tl_u_stokes(i,j,k)
# endif
          END DO
        END DO
        DO i=IstrP,IendT
          DC1(i,0)=DC(i,0)                              ! intermediate
# ifdef NEARSHORE_MELLOR
          cff2=DC(i,0)*ubar_stokes(i,j)
          tl_cff2=tl_DC(i,0)*ubar_stokes(i,j)+                          &
     &            DC(i,0)*tl_ubar_stokes(i,j)
# endif
          DC(i,0)=1.0_r8/DC(i,0)                        ! recursive
          tl_DC(i,0)=-tl_DC(i,0)/(DC1(i,0)*DC1(i,0))
          CF1(i)=CF(i,0)                                ! intermediate
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j))        ! recursive
          tl_CF(i,0)=tl_DC(i,0)*(CF1(i)-DU_avg1(i,j))+                  &
     &               DC(i,0)*(tl_CF(i,0)-tl_DU_avg1(i,j))
# ifdef NEARSHORE_MELLOR
          CFs1(i)=CFs(i,0)                              ! intermediate
          CFs(i,0)=DC(i,0)*(CFs(i,0)-cff2)              ! recursive
          tl_CFs(i,0)=tl_DC(i,0)*(CFs1(i)-cff2)+                        &
     &                DC(i,0)*(tl_CFs(i,0)-tl_cff2)
# endif
!>        ubar(i,j,1)=DC(i,0)*DU_avg1(i,j)
!>
          tl_ubar(i,j,1)=tl_DC(i,0)*DU_avg1(i,j)+                       &
     &                   DC(i,0)*tl_DU_avg1(i,j)
!>        ubar(i,j,2)=ubar(i,j,1)
!>
          tl_ubar(i,j,2)=tl_ubar(i,j,1)
# ifdef DIAGNOSTICS_UV
!!        DiaU2wrk(i,j,M2rate)=ubar(i,j,1)-DiaU2int(i,j,M2rate)*DC(i,0)
!!        DiaU2int(i,j,M2rate)=ubar(i,j,1)*DC1(i,0)
# endif
        END DO
# ifdef DIAGNOSTICS_UV
!!
!! Convert the units of the fast-time integrated change in mass flux
!! diagnostic values to change in velocity (m s-1).
!!
!!      DO idiag=1,NDM2d-1
!!        DO i=IstrP,IendT
!!          DiaU2wrk(i,j,idiag)=DC(i,0)*DiaU2wrk(i,j,idiag)
#  ifdef MASKING
!!          DiaU2wrk(i,j,idiag)=DiaU2wrk(i,j,idiag)*umask(i,j)
#  endif
!!        END DO
!!      END DO
# endif
!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, ubar=DU_avg1/(D*on_u). Recall that,
!  D=CF(:,0).
!
!  NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant
!         update in the interior again for computational purposes which
!         will not affect the nonlinear code.  However, the adjoint
!         code is wrong because the interior solution is corrected
!         twice. The replacement is avoided if the boundary edge is
!         periodic. The J-loop is pipelined, so we need to do a special
!         test for the southern and northern domain edges.
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=1,N(ng)
!>            u(Istr,j,k,nnew)=u(Istr,j,k,nnew)-CF(Istr,0)
!>
              tl_u(Istr,j,k,nnew)=tl_u(Istr,j,k,nnew)-                  &
     &                            tl_CF(Istr,0)
# ifdef MASKING
!>            u(Istr,j,k,nnew)=u(Istr,j,k,nnew)*                        &
!>   &                         umask(Istr,j)
!>
              tl_u(Istr,j,k,nnew)=tl_u(Istr,j,k,nnew)*                  &
     &                            umask(Istr,j)
# endif
# ifdef WET_DRY_NOT_YET
!>            u(Istr,j,k,nnew)=u(Istr,j,k,nnew)*                        &
!>   &                         umask_wet(Istr,j)
!>
              tl_u(Istr,j,k,nnew)=tl_u(Istr,j,k,nnew)*                  &
     &                            umask_wet(Istr,j)
# endif
# ifdef NEARSHORE_MELLOR
!>            u_stokes(Istr,j,k)=u_stokes(Istr,j,k)-CFs(Istr,0)
!>
              tl_u_stokes(Istr,j,k)=tl_u_stokes(Istr,j,k)-              &
     &                              tl_CFs(Istr,0)
#  ifdef MASKING
!>            u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                    &
!>   &                           umask(Istr,j)
!>
              tl_u_stokes(Istr,j,k)=tl_u_stokes(Istr,j,k)*              &
     &                              umask(Istr,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>            u_stokes(Istr,j,k)=u_stokes(Istr,j,k)*                    &
!>   &                           umask_wet(Istr,j)
!>
              tl_u_stokes(Istr,j,k)=tl_u_stokes(Istr,j,k)*              &
     &                              umask_wet(Istr,j)
#  endif
# endif
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=1,N(ng)
!>            u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)-CF(Iend+1,0)
!>
              tl_u(Iend+1,j,k,nnew)=tl_u(Iend+1,j,k,nnew)-              &
     &                              tl_CF(Iend+1,0)
# ifdef MASKING
!>            u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)*                    &
!>   &                           umask(Iend+1,j)
!>
              tl_u(Iend+1,j,k,nnew)=tl_u(Iend+1,j,k,nnew)*              &
     &                              umask(Iend+1,j)
# endif
# ifdef WET_DRY_NOT_YET
!>            u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)*                    &
!>   &                           umask_wet(Iend+1,j)
!>
              tl_u(Iend+1,j,k,nnew)=tl_u(Iend+1,j,k,nnew)*              &
     &                              umask_wet(Iend+1,j)
# endif
# ifdef NEARSHORE_MELLOR
!>            u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)-CFs(Iend+1,0)
!>
              tl_u_stokes(Iend+1,j,k)=tl_u_stokes(Iend+1,j,k)-          &
     &                                tl_CFs(Iend+1,0)
#  ifdef MASKING
!>            u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                &
!>   &                             umask(Iend+1,j)
!>
              tl_u_stokes(Iend+1,j,k)=tl_u_stokes(Iend+1,j,k)*          &
     &                                umask(Iend+1,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>            u_stokes(Iend+1,j,k)=u_stokes(Iend+1,j,k)*                &
!>   &                             umask_wet(Iend+1,j)
!>
              tl_u_stokes(Iend+1,j,k)=tl_u_stokes(Iend+1,j,k)*          &
     &                                umask_wet(Iend+1,j)
#  endif
# endif
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (j.eq.0) THEN                           ! southern boundary
            DO k=1,N(ng)                             ! J-loop pipelined
              DO i=IstrU,Iend
!>              u(i,j,k,nnew)=u(i,j,k,nnew)-CF(i,0)
!>
                tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-                      &
     &                           tl_CF(i,0)
# ifdef MASKING
!>              u(i,j,k,nnew)=u(i,j,k,nnew)*umask(i,j)
!>
                tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*                      &
     &                           umask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>              u(i,j,k,nnew)=u(i,j,k,nnew)*                            &
!>   &                        umask_wet(i,j)
!>
                tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*                      &
     &                           umask_wet(i,j)
# endif
# ifdef NEARSHORE_MELLOR
!>              u_stokes(i,j,k)=u_stokes(i,j,k)-CFs(i,0)
!>
                tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)-                  &
     &                             tl_CFs(i,0)
#  ifdef MASKING
!>              u_stokes(i,j,k)=u_stokes(i,j,k)*                        &
!>   &                          umask(i,j)
!>
                tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*                  &
     &                             umask(i,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>              u_stokes(i,j,k)=u_stokes(i,j,k)*                        &
!>   &                          umask_wet(i,j)
!>
                tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*                  &
     &                             umask_wet(i,j)
#  endif
# endif
              END DO
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (j.eq.Mm(ng)+1) THEN                    ! northern boundary
            DO k=1,N(ng)                             ! J-loop pipelined
              DO i=IstrU,Iend
!>              u(i,j,k,nnew)=u(i,j,k,nnew)-CF(i,0)
!>
                tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-                      &
     &                           tl_CF(i,0)
# ifdef MASKING
!>              u(i,j,k,nnew)=u(i,j,k,nnew)*umask(i,j)
!>
                tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*                      &
     &                           umask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>              u(i,j,k,nnew)=u(i,j,k,nnew)*umask_wet(i,j)
!>
                tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*                      &
     &                           umask_wet(i,j)
# endif
# ifdef NEARSHORE_MELLOR
!>              u_stokes(i,j,k)=u_stokes(i,j,k)-CFs(i,0)
!>
                tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)-                  &
     &                             tl_CFs(i,0)
#  ifdef MASKING
!>              u_stokes(i,j,k)=u_stokes(i,j,k)*                        &
!>   &                          umask(i,j)
!>
                tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*                  &
     &                             umask(i,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>              u_stokes(i,j,k)=u_stokes(i,j,k)*                        &
!>   &                          umask_wet(i,j)
!>
                tl_u_stokes(i,j,k)=tl_u_stokes(i,j,k)*                  &
     &                             umask_wet(i,j)
#  endif
# endif
              END DO
            END DO
          END IF
        END IF
!
!  Compute correct mass flux, Hz*u/n.
!
        DO k=N(ng),1,-1
          DO i=IstrP,IendT
            Huon(i,j,k)=0.5_r8*(Huon(i,j,k)+u(i,j,k,nnew)*DC(i,k))
            tl_Huon(i,j,k)=0.5_r8*(tl_Huon(i,j,k)+                      &
     &                             tl_u(i,j,k,nnew)*DC(i,k)+            &
     &                             u(i,j,k,nnew)*tl_DC(i,k))
# ifdef NEARSHORE_MELLOR
            Huon(i,j,k)=Huon(i,j,k)+0.5_r8*u_stokes(i,j,k)*DC(i,k)
            tl_Huon(i,j,k)=tl_Huon(i,j,k)+                              &
     &                     0.5_r8*(tl_u_stokes(i,j,k)*DC(i,k)+          &
     &                             u_stokes(i,j,k)*tl_DC(i,k))
# endif
            FC(i,0)=FC(i,0)+Huon(i,j,k)
            tl_FC(i,0)=tl_FC(i,0)+tl_Huon(i,j,k)
# ifdef DIAGNOSTICS_UV
!!          DiaU3wrk(i,j,k,M3rate)=u(i,j,k,nnew)-DiaU3wrk(i,j,k,M3rate)
# endif
          END DO
        END DO
        DO i=IstrP,IendT
          FC1(i)=FC(i,0)                                ! intermediate
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))        ! recursive
          tl_FC(i,0)=tl_DC(i,0)*(FC1(i)-DU_avg2(i,j))+                  &
     &               DC(i,0)*(tl_FC(i,0)-tl_DU_avg2(i,j))
        END DO
        DO k=1,N(ng)
          DO i=IstrP,IendT
            Huon(i,j,k)=Huon(i,j,k)-DC(i,k)*FC(i,0)
            tl_Huon(i,j,k)=tl_Huon(i,j,k)-                              &
     &                     tl_DC(i,k)*FC(i,0)-                          &
     &                     DC(i,k)*tl_FC(i,0)
          END DO
        END DO
!
!  Couple velocity component in the ETA-direction.
!
        IF (j.ge.Jstr) THEN
          DO i=IstrT,IendT
            DC(i,0)=0.0_r8
            tl_DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
            tl_CF(i,0)=0.0_r8
# ifdef NEARSHORE_MELLOR
            CFs(i,0)=0.0_r8
            tl_CFs(i,0)=0.0_r8
# endif
            FC(i,0)=0.0_r8
            tl_FC(i,0)=0.0_r8
          END DO
!
!  Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!
          DO k=1,N(ng)
            DO i=IstrT,IendT
              cff=0.5_r8*om_v(i,j)
              DC(i,k)=cff*(Hz(i,j,k)+Hz(i,j-1,k))
              tl_DC(i,k)=cff*(tl_Hz(i,j,k)+tl_Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
              CF(i,0)=CF(i,0)+                                          &
     &                DC(i,k)*v(i,j,k,nnew)
              tl_CF(i,0)=tl_CF(i,0)+                                    &
     &                   tl_DC(i,k)*v(i,j,k,nnew)+                      &
     &                   DC(i,k)*tl_v(i,j,k,nnew)
# ifdef NEARSHORE_MELLOR
              CFs(i,0)=CFs(i,0)+                                        &
     &                 DC(i,k)*v_stokes(i,j,k)
              tl_CFs(i,0)=tl_CFs(i,0)+                                  &
     &                    tl_DC(i,k)*v_stokes(i,j,k)+                   &
     &                    DC(i,k)*tl_v_stokes(i,j,k)
# endif
            END DO
          END DO
          DO i=IstrT,IendT
            DC1(i,0)=DC(i,0)                            ! intermediate
# ifdef NEARSHORE_MELLOR
            cff2=DC(i,0)*vbar_stokes(i,j)
            tl_cff2=tl_DC(i,0)*vbar_stokes(i,j)+                        &
     &              DC(i,0)*tl_vbar_stokes(i,j)
# endif
            DC(i,0)=1.0_r8/DC(i,0)                      ! recursive
            tl_DC(i,0)=-tl_DC(i,0)/(DC1(i,0)*DC1(i,0))
            CF1(i)=CF(i,0)                              ! intermediate
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg1(i,j))      ! recursive
            tl_CF(i,0)=tl_DC(i,0)*(CF1(i)-DV_avg1(i,j))+                &
     &                 DC(i,0)*(tl_CF(i,0)-tl_DV_avg1(i,j))
# ifdef NEARSHORE_MELLOR
            CFs1(i)=CFs(i,0)
            CFs(i,0)=DC(i,0)*(CFs(i,0)-cff2)            ! recursive
            tl_CFs(i,0)=tl_DC(i,0)*(CFs1(i)-cff2)+                      &
     &                  DC(i,0)*(tl_CFs(i,0)-tl_cff2)
# endif
!>          vbar(i,j,1)=DC(i,0)*DV_avg1(i,j)
!>
            tl_vbar(i,j,1)=tl_DC(i,0)*DV_avg1(i,j)+                     &
     &                     DC(i,0)*tl_DV_avg1(i,j)
!>          vbar(i,j,2)=vbar(i,j,1)
!>
            tl_vbar(i,j,2)=tl_vbar(i,j,1)
# ifdef DIAGNOSTICS_UV
!!          DiaV2wrk(i,j,M2rate)=vbar(i,j,1)-                           &
!!   &                           DiaV2int(i,j,M2rate)*DC(i,0)
!!          DiaV2int(i,j,M2rate)=vbar(i,j,1)*DC1(i,0)
!!          DiaV2wrk(i,j,M2rate)=vbar(i,j,1)-DiaV2int(i,j,M2rate)
!!          DiaV2int(i,j,M2rate)=vbar(i,j,1)
# endif
          END DO
# ifdef DIAGNOSTICS_UV
!!
!! Convert the units of the fast-time integrated change in mass flux
!! diagnostic values to change in velocity (m s-1).
!!
!!        DO idiag=1,NDM2d-1
!!          DO i=IstrT,IendT
!!            DiaV2wrk(i,j,idiag)=DC(i,0)*DiaV2wrk(i,j,idiag)
#  ifdef MASKING
!!            DiaV2wrk(i,j,idiag)=DiaV2wrk(i,j,idiag)*vmask(i,j)
#  endif
!!          END DO
!!        END DO
# endif
!
!  Replace only BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, vbar=DV_avg1/(D*om_v).  Recall that,
!  D=CF(:,0).
!
!  NOTE:  Only the BOUNDARY POINTS need to be replaced. Avoid redundant
!         update in the interior again for computational purposes which
!         will not affect the nonlinear code.  However, the adjoint
!         code is wrong because the interior solution is corrected
!         twice. The replacement is avoided if the boundary edge is
!         periodic. The J-loop is pipelined, so we need to do a special
!         test for the southern and northern domain edges.
!
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO k=1,N(ng)
!>              v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)-CF(Istr-1,0)
!>
                tl_v(Istr-1,j,k,nnew)=tl_v(Istr-1,j,k,nnew)-            &
     &                                tl_CF(Istr-1,0)
# ifdef MASKING
!>              v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)*                  &
!>   &                             vmask(Istr-1,j)
!>
                tl_v(Istr-1,j,k,nnew)=tl_v(Istr-1,j,k,nnew)*            &
     &                                vmask(Istr-1,j)
# endif
# ifdef WET_DRY_NOT_YET
!>              v(Istr-1,j,k,nnew)=v(Istr-1,j,k,nnew)*                  &
!>   &                             vmask_wet(Istr-1,j)
!>
                tl_v(Istr-1,j,k,nnew)=tl_v(Istr-1,j,k,nnew)*            &
     &                                vmask_wet(Istr-1,j)
# endif
# ifdef NEARSHORE_MELLOR
!>              v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)-              &
!>   &                               CFs(Istr-1,0)
!>
                tl_v_stokes(Istr-1,j,k)=tl_v_stokes(Istr-1,j,k)-        &
     &                                  tl_CFs(Istr-1,0)
#  ifdef MASKING
!>              v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
!>   &                               vmask(Istr-1,j)
!>
                tl_v_stokes(Istr-1,j,k)=tl_v_stokes(Istr-1,j,k)*        &
     &                                  vmask(Istr-1,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>              v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
!>   &                               vmask_wet(Istr-1,j)
!>
                tl_v_stokes(Istr-1,j,k)=tl_v_stokes(Istr-1,j,k)*        &
     &                                  vmask_wet(Istr-1,j)
#  endif
# endif
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO k=1,N(ng)
!>              v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)-CF(Iend+1,0)
!>
                tl_v(Iend+1,j,k,nnew)=tl_v(Iend+1,j,k,nnew)-            &
     &                                tl_CF(Iend+1,0)
# ifdef MASKING
!>              v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)*                  &
!>   &                             vmask(Iend+1,j)
!>
                tl_v(Iend+1,j,k,nnew)=tl_v(Iend+1,j,k,nnew)*            &
     &                                vmask(Iend+1,j)
# endif
# ifdef WET_DRY_NOT_YET
!>              v(Iend+1,j,k,nnew)=v(Iend+1,j,k,nnew)*                  &
!>   &                             vmask_wet(Iend+1,j)
!>
                tl_v(Iend+1,j,k,nnew)=tl_v(Iend+1,j,k,nnew)*            &
     &                                vmask_wet(Iend+1,j)
# endif
# ifdef NEARSHORE_MELLOR
!>              v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)-              &
!>   &                               CFs(Iend+1,0)
!>
                tl_v_stokes(Iend+1,j,k)=tl_v_stokes(Iend+1,j,k)-        &
     &                                  tl_CFs(Iend+1,0)
#  ifdef MASKING
!>              v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
!>   &                               vmask(Iend+1,j)
!>
                tl_v_stokes(Iend+1,j,k)=tl_v_stokes(Iend+1,j,k)*        &
     &                                  vmask(Iend+1,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>              v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
!>   &                               vmask_wet(Iend+1,j)
!>
                tl_v_stokes(Iend+1,j,k)=tl_v_stokes(Iend+1,j,k)*        &
     &                                  vmask_wet(Iend+1,j)
#  endif
# endif
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (j.eq.1) THEN                         ! southern boundary
              DO k=1,N(ng)                           ! J-loop pipelined
                DO i=Istr,Iend
!>                v(i,j,k,nnew)=v(i,j,k,nnew)-CF(i,0)
!>
                  tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-                    &
     &                             tl_CF(i,0)
# ifdef MASKING
!>                v(i,j,k,nnew)=v(i,j,k,nnew)*
!>   &                            vmask(i,j)
!>
                  tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*                    &
     &                             vmask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>                v(i,j,k,nnew)=v(i,j,k,nnew)*                          &
!>   &                          vmask_wet(i,j)
!>
                  tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*                    &
     &                             vmask_wet(i,j)
# endif
# ifdef NEARSHORE_MELLOR
!>                v_stokes(i,j,k)=v_stokes(i,j,k)-CFs(i,0)
!>
                  tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)-                &
     &                               tl_CFs(i,0)
#  ifdef MASKING
!>                v_stokes(i,j,k)=v_stokes(i,j,k)*vmask(i,j)
!>
                  tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*                &
     &                               vmask(i,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>                v_stokes(i,j,k)=v_stokes(i,j,k)*                      &
!>   &                            vmask_wet(i,j)
!>
                  tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*                &
     &                               vmask_wet(i,j)
#  endif
# endif
                END DO
              END DO
            END IF
          END IF
!
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (j.eq.Mm(ng)+1) THEN                  ! northern boundary
              DO k=1,N(ng)                           ! J-loop pipelined
                DO i=Istr,Iend
!>                v(i,j,k,nnew)=v(i,j,k,nnew)-CF(i,0)
!>
                  tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-                    &
     &                             tl_CF(i,0)
# ifdef MASKING
!>                v(i,j,k,nnew)=v(i,j,k,nnew)*                          &
!>   &                          vmask(i,j)
!>
                  tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*                    &
     &                             vmask(i,j)
# endif
# ifdef WET_DRY_NOT_YET
!>                v(i,j,k,nnew)=v(i,j,k,nnew)*                          &
!>   &                          vmask_wet(i,j)
!>
                  tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*                    &
     &                             vmask_wet(i,j)
# endif
# ifdef NEARSHORE_MELLOR
!>                v_stokes(i,j,k)=v_stokes(i,j,k)-CFs(i,0)
!>
                  tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)-                &
     &                               tl_CFs(i,0)
#  ifdef MASKING
!>                v_stokes(i,j,k)=v_stokes(i,j,k)*                      &
!>   &                            vmask(i,j)
!>
                  tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*                &
     &                               vmask(i,j)
#  endif
#  ifdef WET_DRY_NOT_YET
!>                v_stokes(i,j,k)=v_stokes(i,j,k)*                      &
!>   &                            vmask_wet(i,j)
!>
                  tl_v_stokes(i,j,k)=tl_v_stokes(i,j,k)*                &
     &                               vmask_wet(i,j)
#  endif
# endif
                END DO
              END DO
            END IF
          END IF
!
!  Compute correct mass flux, Hz*v/m.
!
          DO k=N(ng),1,-1
            DO i=IstrT,IendT
              Hvom(i,j,k)=0.5_r8*(Hvom(i,j,k)+v(i,j,k,nnew)*DC(i,k))
              tl_Hvom(i,j,k)=0.5_r8*(tl_Hvom(i,j,k)+                    &
     &                               tl_v(i,j,k,nnew)*DC(i,k)+          &
     &                               v(i,j,k,nnew)*tl_DC(i,k))
# ifdef NEARSHORE_MELLOR
              Hvom(i,j,k)=Hvom(i,j,k)+0.5_r8*v_stokes(i,j,k)*DC(i,k)
              tl_Hvom(i,j,k)=tl_Hvom(i,j,k)+                            &
     &                       0.5_r8*(tl_v_stokes(i,j,k)*DC(i,k)+        &
     &                               v_stokes(i,j,k)*tl_DC(i,k))
# endif
              FC(i,0)=FC(i,0)+Hvom(i,j,k)
              tl_FC(i,0)=tl_FC(i,0)+tl_Hvom(i,j,k)
# ifdef DIAGNOSTICS_UV
!!            DiaV3wrk(i,j,k,M3rate)=v(i,j,k,nnew)-                     &
!!   &                               DiaV3wrk(i,j,k,M3rate)
# endif
            END DO
          END DO
          DO i=IstrT,IendT
            FC1(i)=FC(i,0)                              ! intermediate
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))      ! recursive
            tl_FC(i,0)=tl_DC(i,0)*(FC1(i)-DV_avg2(i,j))+                &
     &                 DC(i,0)*(tl_FC(i,0)-tl_DV_avg2(i,j))
          END DO
          DO k=1,N(ng)
            DO i=IstrT,IendT
              Hvom(i,j,k)=Hvom(i,j,k)-DC(i,k)*FC(i,0)
              tl_Hvom(i,j,k)=tl_Hvom(i,j,k)-                            &
     &                       tl_DC(i,k)*FC(i,0)-                        &
     &                       DC(i,k)*tl_FC(i,0)
            END DO
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_u3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          u(:,:,:,nnew))
!>
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_u(:,:,:,nnew))
!>      CALL exchange_v3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          v(:,:,:,nnew))
!>
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_v(:,:,:,nnew))

        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Huon)
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_Huon)

        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Hvom)
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_Hvom)
!
        DO k=1,2
!>        CALL exchange_u2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            ubar(:,:,k))
!>
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_ubar(:,:,k))
!>        CALL exchange_v2d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            vbar(:,:,k))
!>
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_vbar(:,:,k))
        END DO
      END IF

# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iNLM, 4,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    u(:,:,:,nnew), v(:,:,:,nnew),                 &
!>   &                    Huon, Hvom)
!>
      CALL mp_exchange3d (ng, tile, iTLM, 4,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_u(:,:,:,nnew), tl_v(:,:,:,nnew),           &
     &                    tl_Huon, tl_Hvom)
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Huon, Hvom)
!
!>    CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    ubar(:,:,1), vbar(:,:,1),                     &
!>   &                    ubar(:,:,2), vbar(:,:,2))
!>
      CALL mp_exchange2d (ng, tile, iTLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_ubar(:,:,1), tl_vbar(:,:,1),               &
     &                    tl_ubar(:,:,2), tl_vbar(:,:,2))
# endif

      RETURN
      END SUBROUTINE tl_step3d_uv_tile
#endif
      END MODULE tl_step3d_uv_mod
