#include "cppdefs.h"
      MODULE ad_pre_step3d_mod

#if defined ADJOINT && defined SOLVE3D
# define NEUMANN
!
!svn $Id: ad_pre_step3d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine initialize computations for new time step of the    !
!  adjoint 3D primitive variables.                                     !
!                                                                      !
!  Both n-1 and n-2 time-step contributions of the  Adams/Bashforth    !
!  scheme are added here to u and v at time index "nnew", since the    !
!  right-hand-side  arrays ru and rv at  n-2  will be overwriten in    !
!  subsequent calls to routines within the 3D engine.                  !
!                                                                      !
!  It also computes the time  "n"  vertical viscosity and diffusion    !
!  contributions of the Crank-Nicholson implicit scheme because the    !
!  thicknesses "Hz" will be overwriten at the end of the  2D engine    !
!  (barotropic mode) computations.                                     !
!                                                                      !
!  The actual time step will be carried out in routines "step3d_uv"    !
!  and "step3d_t".                                                     !
!                                                                      !
!  BASIC STATE variables needed: AKt, AKv, Huon, Hvom, Hz, Tsrc, W,    !
!                                bustr, bvstr, ghats, srflx, sustr,    !
!                                svstr,  t, z_r, z_w                   !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ad_pre_step3d
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_pre_step3d (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef DIAGNOSTICS
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
      CALL wclock_on (ng, iADM, 22)
# endif
      CALL ad_pre_step3d_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
# ifdef FOUR_DVAR
     &                         kstp(ng), knew(ng),                      &
# endif
     &                         nrhs(ng), nstp(ng), nnew(ng),            &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
# endif
# ifdef FOUR_DVAR
     &                         GRID(ng) % om_v,                         &
     &                         GRID(ng) % on_u,                         &
# endif
     &                         GRID(ng) % pm,                           &
     &                         GRID(ng) % pn,                           &
     &                         GRID(ng) % Hz,                           &
     &                         GRID(ng) % ad_Hz,                        &
     &                         GRID(ng) % Huon,                         &
     &                         GRID(ng) % ad_Huon,                      &
     &                         GRID(ng) % Hvom,                         &
     &                         GRID(ng) % ad_Hvom,                      &
     &                         GRID(ng) % z_r,                          &
     &                         GRID(ng) % ad_z_r,                       &
     &                         GRID(ng) % z_w,                          &
     &                         GRID(ng) % ad_z_w,                       &
     &                         FORCES(ng) % ad_btflx,                   &
     &                         FORCES(ng) % ad_bustr,                   &
     &                         FORCES(ng) % ad_bvstr,                   &
     &                         FORCES(ng) % ad_stflx,                   &
     &                         FORCES(ng) % ad_sustr,                   &
     &                         FORCES(ng) % ad_svstr,                   &
# ifdef SOLAR_SOURCE
     &                         FORCES(ng) % srflx,                      &
# endif
     &                         MIXING(ng) % Akt,                        &
     &                         MIXING(ng) % ad_Akt,                     &
     &                         MIXING(ng) % Akv,                        &
     &                         MIXING(ng) % ad_Akv,                     &
# ifdef LMD_NONLOCAL_NOT_YET
     &                         MIXING(ng) % ad_ghats,                   &
# endif
# ifdef FOUR_DVAR
     &                         OCEAN(ng) % ad_ubar,                     &
     &                         OCEAN(ng) % ad_vbar,                     &
# endif
     &                         OCEAN(ng) % W,                           &
     &                         OCEAN(ng) % ad_W,                        &
     &                         OCEAN(ng) % ad_ru,                       &
     &                         OCEAN(ng) % ad_rv,                       &
# ifdef DIAGNOSTICS_TS
!!   &                         DIAGS(ng) % DiaTwrk,                     &
# endif
# ifdef DIAGNOSTICS_UV
!!   &                         DIAGS(ng) % DiaU3wrk,                    &
!!   &                         DIAGS(ng) % DiaV3wrk,                    &
!!   &                         DIAGS(ng) % DiaRU,                       &
!!   &                         DIAGS(ng) % DiaRV,                       &
# endif
     &                         OCEAN(ng) % t,                           &
     &                         OCEAN(ng) % ad_t,                        &
     &                         OCEAN(ng) % u,                           &
     &                         OCEAN(ng) % ad_u,                        &
     &                         OCEAN(ng) % v,                           &
     &                         OCEAN(ng) % ad_v)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 22)
# endif
      RETURN
      END SUBROUTINE ad_pre_step3d
!
!***********************************************************************
      SUBROUTINE ad_pre_step3d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
# ifdef FOUR_DVAR
     &                               kstp, knew,                        &
# endif
     &                               nrhs, nstp, nnew,                  &
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
# ifdef FOUR_DVAR
     &                               om_v, on_u,                        &
# endif
     &                               pm, pn,                            &
     &                               Hz, ad_Hz,                         &
     &                               Huon, ad_Huon,                     &
     &                               Hvom, ad_Hvom,                     &
     &                               z_r, ad_z_r,                       &
     &                               z_w, ad_z_w,                       &
     &                               ad_btflx, ad_bustr, ad_bvstr,      &
     &                               ad_stflx, ad_sustr, ad_svstr,      &
# ifdef SOLAR_SOURCE
     &                               srflx,                             &
# endif
     &                               Akt, ad_Akt,                       &
     &                               Akv, ad_Akv,                       &
# ifdef LMD_NONLOCAL_NOT_YET
     &                               ad_ghats,                          &
# endif
# ifdef FOUR_DVAR
     &                               ad_ubar, ad_vbar,                  &
# endif
     &                               W, ad_W,                           &
     &                               ad_ru, ad_rv,                      &
# ifdef DIAGNOSTICS_TS
!!   &                               DiaTwrk,                           &
# endif
# ifdef DIAGNOSTICS_UV
!!   &                               DiaU3wrk, DiaV3wrk,                &
!!   &                               DiaRU, DiaRV,                      &
# endif
     &                               t, ad_t,                           &
     &                               u, ad_u,                           &
     &                               v, ad_v)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sources
!
      USE ad_exchange_3d_mod, ONLY : ad_exchange_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange4d
# endif
      USE ad_t3dbc_mod, ONLY : ad_t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef FOUR_DVAR
      integer, intent(in) :: kstp, knew
# endif
      integer, intent(in) :: nrhs, nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef FOUR_DVAR
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  ifdef SOLAR_SOURCE
      real(r8), intent(in) :: srflx(LBi:,LBj:)
#  endif
#  ifdef SUN
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  else
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
#  endif
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
#  ifdef FOUR_DVAR
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
#  ifdef SUN
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)

#  ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
#  endif
#  ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
!!    real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Hvom(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_btflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_bustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_bvstr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: ad_rv(LBi:,LBj:,0:,:)
#  ifdef LMD_NONLOCAL_NOT_YET
      real(r8), intent(inout) :: ad_ghats(LBi:,LBj:,0:,:)
#  endif
#  ifdef SUN
      real(r8), intent(inout) :: ad_Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  else
      real(r8), intent(inout) :: ad_Akt(LBi:,LBj:,0:,:)
#  endif
      real(r8), intent(inout) :: ad_Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_W(LBi:,LBj:,0:)
#  ifdef SUN
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef FOUR_DVAR
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef SOLAR_SOURCE
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(in) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef FOUR_DVAR
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)

#  ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
!!   &                                   NDT)
#  endif
#  ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
!!    real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
!!    real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_btflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: ad_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: ad_rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#  ifdef LMD_NONLOCAL_NOT_YET
      real(r8), intent(inout) :: ad_ghats(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  endif
      real(r8), intent(inout) :: ad_Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: ad_Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: Isrc, Jsrc
      integer :: i, ic, indx, is, itrc, j, k, ltrc
# if defined AGE_MEAN && defined T_PASSIVE
      integer :: iage
# endif
# if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_UV
      integer :: idiag
# endif
# ifdef TS_MPDATA_NOT_YET
      real(r8), parameter :: Gamma = 0.5_r8
# else
      real(r8), parameter :: Gamma = 1.0_r8/6.0_r8
# endif
      real(r8), parameter :: eps = 1.0E-16_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: ad_cff, ad_cff1, ad_cff2, ad_cff3, ad_cff4
      real(r8) :: adfac, adfac1, adfac2, adfac3

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
     real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_FC

# ifdef SOLAR_SOURCE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: ad_swdk
# endif

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_grad

# include "set_bounds.h"
!
!  In the tangent model, we have:
!
!>                          ! nnew=3-nstp
!>    indx=3-nrhs           ! nrhs=nstp
!>
      indx=nnew
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
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          ad_FE(i,j)=0.0_r8
          ad_FX(i,j)=0.0_r8
          ad_curv(i,j)=0.0_r8
          ad_grad(i,j)=0.0_r8
        END DO
# ifdef SOLAR_SOURCE
        DO k=0,N(ng)
          DO i=IminS,ImaxS
            ad_swdk(i,j,k)=0.0_r8
          END DO
        END DO
# endif
      END DO
      DO k=0,N(ng)
        DO i=IminS,ImaxS
          ad_CF(i,k)=0.0_r8
          ad_DC(i,k)=0.0_r8
          ad_FC(i,k)=0.0_r8
        END DO
      END DO

# ifndef TS_FIXED
!
!=======================================================================
!  Apply tracers lateral boundary conditions.
!=======================================================================
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange4d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_t(:,:,:,3,:))
!>
      CALL ad_mp_exchange4d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_t(:,:,:,3,:))
!
#  endif
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>        CALL exchange_r3d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
!>   &                            tl_t(:,:,:,3,itrc))
!>
          CALL ad_exchange_r3d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               ad_t(:,:,:,3,itrc))
        END IF
!>      CALL tl_t3dbc_tile (ng, tile, itrc, ic,                         &
!>   &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      nstp, 3,                                    &
!>   &                      tl_t)
!>
        CALL ad_t3dbc_tile (ng, tile, itrc, ic,                         &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, 3,                                    &
     &                      ad_t)
      END DO
# endif
!
!=======================================================================
!  3D adjoint momentum equation in the ETA-direction.
!=======================================================================
!
      J_LOOP2 : DO j=Jstr,Jend
        IF (j.ge.JstrV) THEN
!
!  Compute new V-momentum (m m/s).
!
          cff=dt(ng)*0.25_r8
          DO i=Istr,Iend
            DC(i,0)=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          END DO
          IF (iic(ng).eq.ntfirst(ng)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
# ifdef DIAGNOSTICS_UV
!!              DiaV3wrk(i,j,k,M3rate)=cff1
!!              DiaV3wrk(i,j,k,M3vvis)=cff2
!!              DO idiag=1,M3pgrd
!!                DiaV3wrk(i,j,k,idiag)=0.0_r8
!!              END DO
# endif
!>              tl_v(i,j,k,nnew)=tl_cff1+tl_cff2
!>
                ad_cff1=ad_cff1+ad_v(i,j,k,nnew)
                ad_cff2=ad_cff2+ad_v(i,j,k,nnew)
                ad_v(i,j,k,nnew)=0.0_r8
!>              tl_cff2=tl_FC(i,k)-tl_FC(i,k-1)
!>
                ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff2
                ad_FC(i,k  )=ad_FC(i,k  )+ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=0.5_r8*(tl_v(i,j,k,nstp)*                       &
!>   &                          (Hz(i,j,k)+Hz(i,j-1,k))+                &
!>   &                          v(i,j,k,nstp)*                          &
!>   &                          (tl_Hz(i,j,k)+tl_Hz(i,j-1,k)))
!>
                adfac=0.5_r8*ad_cff1
                adfac1=adfac*v(i,j,k,nstp)
                ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+                      &
     &                           (Hz(i,j,k)+Hz(i,j-1,k))*adfac
                ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
                ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
                ad_cff1=0.0_r8
              END DO
            END DO
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff3=0.5_r8*DC(i,0)
# ifdef DIAGNOSTICS_UV
!!              DiaV3wrk(i,j,k,M3rate)=cff1
#  ifdef BODYFORCE
!!              DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)-          &
!!   &                                 cff3*DiaRV(i,j,k,indx,M3vvis)
#  endif
!!              DiaV3wrk(i,j,k,M3vvis)=cff2
!!              DO idiag=1,M3pgrd
!!                DiaV3wrk(i,j,k,idiag)=-cff3*DiaRV(i,j,k,indx,idiag)
!!              END DO
# endif
!>              tl_v(i,j,k,nnew)=tl_cff1-                               &
!>   &                           cff3*tl_rv(i,j,k,indx)+                &
!>   &                           tl_cff2
!
                ad_rv(i,j,k,indx)=ad_rv(i,j,k,indx)-                    &
     &                            cff3*ad_v(i,j,k,nnew)
                ad_cff1=ad_cff1+ad_v(i,j,k,nnew)
                ad_cff2=ad_cff2+ad_v(i,j,k,nnew)
                ad_v(i,j,k,nnew)=0.0_r8
!>              tl_cff2=tl_FC(i,k)-tl_FC(i,k-1)
!>
                ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff2
                ad_FC(i,k  )=ad_FC(i,k  )+ad_cff2
                ad_cff2=0.0_r8
!>              tl_cff1=0.5_r8*(tl_v(i,j,k,nstp)*                       &
!>   &                          (Hz(i,j,k)+Hz(i,j-1,k))+                &
!>   &                          v(i,j,k,nstp)*                          &
!>   &                          (tl_Hz(i,j,k)+tl_Hz(i,j-1,k)))
!>
                adfac=0.5_r8*ad_cff1
                adfac1=adfac*v(i,j,k,nstp)
                ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
                ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
                ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+                      &
     &                           (Hz(i,j,k)+Hz(i,j-1,k))*adfac
                ad_cff1=0.0_r8
              END DO
            END DO
          ELSE
            cff1= 5.0_r8/12.0_r8
            cff2=16.0_r8/12.0_r8
            DO k=1,N(ng)
              DO i=Istr,Iend
# ifdef DIAGNOSTICS_UV
!!              DiaV3wrk(i,j,k,M3rate)=cff3
#  ifdef BODYFORCE
!!              DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+          &
!!   &                                 DC(i,0)*                         &
!!   &                                 (cff1*DiaRV(i,j,k,nrhs,M3vvis)-  &
!!   &                                  cff2*DiaRV(i,j,k,indx,M3vvis))
#  endif
!!              DiaV3wrk(i,j,k,M3vvis)=cff4
!!              DO idiag=1,M3pgrd
!!                DiaV3wrk(i,j,k,idiag)=DC(i,0)*                        &
!!   &                                  (cff1*DiaRV(i,j,k,nrhs,idiag)-  &
!!   &                                   cff2*DiaRV(i,j,k,indx,idiag))
!!              END DO
# endif
!>              tl_v(i,j,k,nnew)=tl_cff3+                               &
!>   &                           DC(i,0)*(cff1*tl_rv(i,j,k,nrhs)-       &
!>   &                                    cff2*tl_rv(i,j,k,indx))+      &
!>   &                           tl_cff4
!>
                adfac=DC(i,0)*ad_v(i,j,k,nnew)
                ad_rv(i,j,k,nrhs)=ad_rv(i,j,k,nrhs)+cff1*adfac
                ad_rv(i,j,k,indx)=ad_rv(i,j,k,indx)-cff2*adfac
                ad_cff3=ad_cff3+ad_v(i,j,k,nnew)
                ad_cff4=ad_cff4+ad_v(i,j,k,nnew)
                ad_v(i,j,k,nnew)=0.0_r8
!>              tl_cff4=tl_FC(i,k)-tl_FC(i,k-1)
!>
                ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff4
                ad_FC(i,k  )=ad_FC(i,k  )+ad_cff4
                ad_cff4=0.0_r8
!>              tl_cff3=0.5_r8*(tl_v(i,j,k,nstp)*                       &
!>   &                          (Hz(i,j,k)+Hz(i,j-1,k))+                &
!>   &                          v(i,j,k,nstp)*                          &
!>   &                          (tl_Hz(i,j,k)+tl_Hz(i,j-1,k)))
!>
                adfac=0.5_r8*ad_cff3
                adfac1=adfac*v(i,j,k,nstp)
                ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
                ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
                ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+                      &
     &                           (Hz(i,j,k)+Hz(i,j-1,k))*adfac
                ad_cff3=0.0_r8
              END DO
            END DO
          END IF
!
!  Apply bottom and surface stresses, if so is prescribed.
!
          DO i=Istr,Iend
# ifdef BODYFORCE
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
# else
!>          tl_FC(i,N(ng))=dt(ng)*tl_svstr(i,j)
!>
            ad_svstr(i,j)=ad_svstr(i,j)+dt(ng)*ad_FC(i,N(ng))
            ad_FC(i,N(ng))=0.0_r8
!>          tl_FC(i,0)=dt(ng)*tl_bvstr(i,j)
!>
            ad_bvstr(i,j)=ad_bvstr(i,j)+dt(ng)*ad_FC(i,0)
            ad_FC(i,0)=0.0_r8
# endif
          END DO
!
!  Compute adjoint V-component viscous vertical momentum fluxes "FC"
!  at current time-step n, and at horizontal V-points and vertical
!  W-points.
!
          cff3=dt(ng)*(1.0_r8-lambda)
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-                  &
     &                    z_r(i,j,k  )-z_r(i,j-1,k  ))
!>            tl_FC(i,k)=cff3*                                          &
!>   &                   (cff*((tl_v(i,j,k+1,nstp)-tl_v(i,j,k,nstp))*   &
!>   &                         (Akv(i,j,k)+Akv(i,j-1,k))+               &
!>   &                         (v(i,j,k+1,nstp)-v(i,j,k,nstp))*         &
!>   &                         (tl_Akv(i,j,k)+tl_Akv(i,j-1,k)))+        &
!>   &                    tl_cff*(v(i,j,k+1,nstp)-v(i,j,k,nstp))*       &
!>   &                           (Akv(i,j,k)+Akv(i,j-1,k)))
!>
              adfac=cff3*ad_FC(i,k)
              adfac1=adfac*cff
              adfac2=adfac1*(Akv(i,j,k)+Akv(i,j-1,k))
              adfac3=adfac1*(v(i,j,k+1,nstp)-v(i,j,k,nstp))
              ad_v(i,j,k  ,nstp)=ad_v(i,j,k  ,nstp)-adfac2
              ad_v(i,j,k+1,nstp)=ad_v(i,j,k+1,nstp)+adfac2
              ad_Akv(i,j-1,k)=ad_Akv(i,j-1,k)+adfac3
              ad_Akv(i,j  ,k)=ad_Akv(i,j  ,k)+adfac3
              ad_cff=ad_cff+                                            &
     &               (v(i,j,k+1,nstp)-v(i,j,k,nstp))*                   &
     &               (Akv(i,j,k)+Akv(i,j-1,k))*adfac
              ad_FC(i,k)=0.0_r8
!>            tl_cff=-cff*cff*(tl_z_r(i,j,k+1)+tl_z_r(i,j-1,k+1)-       &
!>   &                         tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  ))
!>
              adfac=-cff*cff*ad_cff
              ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-adfac
              ad_z_r(i,j  ,k  )=ad_z_r(i,j  ,k  )-adfac
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)+adfac
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+adfac
              ad_cff=0.0_r8
            END DO
          END DO
        END IF
!
!=====================================================================
!  3D adjoint momentum equation in the XI-direction.
!=====================================================================
!
!  Compute new U-momentum (m m/s).
!
        cff=dt(ng)*0.25_r8
        DO i=IstrU,Iend
          DC(i,0)=cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        END DO
        IF (iic(ng).eq.ntfirst(ng)) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
# ifdef DIAGNOSTICS_UV
!!            DiaU3wrk(i,j,k,M3rate)=cff1
!!            DiaU3wrk(i,j,k,M3vvis)=cff2
!!            DO idiag=1,M3pgrd
!!              DiaU3wrk(i,j,k,idiag)=0.0_r8
!!            END DO
# endif
!>            tl_u(i,j,k,nnew)=tl_cff1+tl_cff2
!>
              ad_cff1=ad_cff1+ad_u(i,j,k,nnew)
              ad_cff2=ad_cff2+ad_u(i,j,k,nnew)
              ad_u(i,j,k,nnew)=0.0_r8
!>            tl_cff2=tl_FC(i,k)-tl_FC(i,k-1)
!>
              ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff2
              ad_FC(i,k  )=ad_FC(i,k  )+ad_cff2
              ad_cff2=0.0_r8
!>            tl_cff1=0.5_r8*(tl_u(i,j,k,nstp)*                         &
!>   &                        (Hz(i,j,k)+Hz(i-1,j,k))+                  &
!>   &                        u(i,j,k,nstp)*                            &
!>   &                        (tl_Hz(i,j,k)+tl_Hz(i-1,j,k)))
!>
              adfac=0.5_r8*ad_cff1
              adfac1=adfac*u(i,j,k,nstp)
              ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
              ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
              ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+                        &
     &                         (Hz(i,j,k)+Hz(i-1,j,k))*adfac
              ad_cff1=0.0_r8
            END DO
          END DO
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff3=0.5_r8*DC(i,0)
# ifdef DIAGNOSTICS_UV
!!            DiaU3wrk(i,j,k,M3rate)=cff1
#  ifdef BODYFORCE
!!            DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)-            &
!!   &                               cff3*DiaRU(i,j,k,indx,M3vvis)
#  endif
!!            DO idiag=1,M3pgrd
!!              DiaU3wrk(i,j,k,idiag)=-cff3*DiaRU(i,j,k,indx,idiag)
!!            END DO
!!            DiaU3wrk(i,j,k,M3vvis)=cff2
# endif
!>            tl_u(i,j,k,nnew)=tl_cff1-                                 &
!>   &                         cff3*tl_ru(i,j,k,indx)+                  &
!>   &                         tl_cff2
!>
              ad_ru(i,j,k,indx)=ad_ru(i,j,k,indx)-                      &
     &                          cff3*ad_u(i,j,k,nnew)
              ad_cff1=ad_cff1+ad_u(i,j,k,nnew)
              ad_cff2=ad_cff2+ad_u(i,j,k,nnew)
              ad_u(i,j,k,nnew)=0.0_r8
!>            tl_cff2=tl_FC(i,k)-tl_FC(i,k-1)
!>
              ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff2
              ad_FC(i,k  )=ad_FC(i,k  )+ad_cff2
              ad_cff2=0.0_r8
!>            tl_cff1=0.5_r8*(tl_u(i,j,k,nstp)*                         &
!>   &                        (Hz(i,j,k)+Hz(i-1,j,k))+                  &
!>   &                        u(i,j,k,nstp)*                            &
!>   &                        (tl_Hz(i,j,k)+tl_Hz(i-1,j,k)))
!>
              adfac=0.5_r8*ad_cff1
              adfac1=adfac*u(i,j,k,nstp)
              ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
              ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
              ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+                        &
     &                         (Hz(i,j,k)+Hz(i-1,j,k))*adfac
              ad_cff1=0.0_r8
            END DO
          END DO
        ELSE
          cff1= 5.0_r8/12.0_r8
          cff2=16.0_r8/12.0_r8
          DO k=1,N(ng)
            DO i=IstrU,Iend
# ifdef DIAGNOSTICS_UV
!!            DiaU3wrk(i,j,k,M3rate)=cff3
#  ifdef BODYFORCE
!!            DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+            &
!!   &                               DC(i,0)*                           &
!!   &                               (cff1*DiaRU(i,j,k,nrhs,M3vvis)-    &
!!   &                                cff2*DiaRU(i,j,k,indx,M3vvis))
#  endif
!!            DiaU3wrk(i,j,k,M3vvis)=cff4
!!            DO idiag=1,M3pgrd
!!              DiaU3wrk(i,j,k,idiag)=DC(i,0)*                          &
!!   &                                (cff1*DiaRU(i,j,k,nrhs,idiag)-    &
!!   &                                 cff2*DiaRU(i,j,k,indx,idiag))
!!            END DO
# endif
!>            tl_u(i,j,k,nnew)=tl_cff3+                                 &
!>   &                         DC(i,0)*(cff1*tl_ru(i,j,k,nrhs)-         &
!>   &                                  cff2*tl_ru(i,j,k,indx))+        &
!>   &                         tl_cff4
!>
              adfac=DC(i,0)*ad_u(i,j,k,nnew)
              ad_ru(i,j,k,nrhs)=ad_ru(i,j,k,nrhs)+cff1*adfac
              ad_ru(i,j,k,indx)=ad_ru(i,j,k,indx)-cff2*adfac
              ad_cff3=ad_cff3+ad_u(i,j,k,nnew)
              ad_cff4=ad_cff4+ad_u(i,j,k,nnew)
              ad_u(i,j,k,nnew)=0.0_r8
!>            tl_cff4=tl_FC(i,k)-tl_FC(i,k-1)
!>
              ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff4
              ad_FC(i,k  )=ad_FC(i,k  )+ad_cff4
              ad_cff4=0.0_r8
!>            tl_cff3=0.5_r8*(tl_u(i,j,k,nstp)*                         &
!>   &                        (Hz(i,j,k)+Hz(i-1,j,k))+                  &
!>   &                        u(i,j,k,nstp)*                            &
!>   &                        (tl_Hz(i,j,k)+tl_Hz(i-1,j,k)))
!>
              adfac=0.5_r8*ad_cff3
              adfac1=adfac*u(i,j,k,nstp)
              ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
              ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
              ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+                        &
     &                         (Hz(i,j,k)+Hz(i-1,j,k))*adfac
              ad_cff3=0.0_r8
            END DO
          END DO
        END IF
!
!  Apply bottom and surface stresses, if so is prescribed.
!
        DO i=IstrU,Iend
# ifdef BODYFORCE
!>        tl_FC(i,N(ng))=0.0_r8
!>
          ad_FC(i,N(ng))=0.0_r8
!>        tl_FC(i,0)=0.0_r8
!>
          ad_FC(i,0)=0.0_r8
# else
!>        tl_FC(i,N(ng))=dt(ng)*tl_sustr(i,j)
!>
          ad_sustr(i,j)=ad_sustr(i,j)+dt(ng)*ad_FC(i,N(ng))
          ad_FC(i,N(ng))=0.0_r8
!>        tl_FC(i,0)=dt(ng)*tl_bustr(i,j)
!>
          ad_bustr(i,j)=ad_bustr(i,j)+dt(ng)*ad_FC(i,0)
          ad_FC(i,0)=0.0_r8
# endif
        END DO
!
!  Compute adjoint U-component viscous vertical momentum fluxes "FC"
!  at current time-step n, and at horizontal U-points and vertical
!  W-points.
!
        cff3=dt(ng)*(1.0_r8-lambda)
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0/(z_r(i,j,k+1)+z_r(i-1,j,k+1)-                       &
     &               z_r(i,j,k  )-z_r(i-1,j,k  ))
!>          tl_FC(i,k)=cff3*                                            &
!>   &                 (cff*((tl_u(i,j,k+1,nstp)-tl_u(i,j,k,nstp))*     &
!>   &                       (Akv(i,j,k)+Akv(i-1,j,k))+                 &
!>   &                       (u(i,j,k+1,nstp)-u(i,j,k,nstp))*           &
!>   &                       (tl_Akv(i,j,k)+tl_Akv(i-1,j,k)))+          &
!>   &                  tl_cff*(u(i,j,k+1,nstp)-u(i,j,k,nstp))*         &
!>   &                         (Akv(i,j,k)+Akv(i-1,j,k)))
!>
            adfac=cff3*ad_FC(i,k)
            adfac1=adfac*cff
            adfac2=adfac1*(Akv(i,j,k)+Akv(i-1,j,k))
            adfac3=adfac1*(u(i,j,k+1,nstp)-u(i,j,k,nstp))
            ad_u(i,j,k  ,nstp)=ad_u(i,j,k  ,nstp)-adfac2
            ad_u(i,j,k+1,nstp)=ad_u(i,j,k+1,nstp)+adfac2
            ad_Akv(i-1,j,k)=ad_Akv(i-1,j,k)+adfac3
            ad_Akv(i  ,j,k)=ad_Akv(i  ,j,k)+adfac3
            ad_cff=ad_cff+                                              &
     &             (u(i,j,k+1,nstp)-u(i,j,k,nstp))*                     &
     &             (Akv(i,j,k)+Akv(i-1,j,k))*adfac
            ad_FC(i,k)=0.0_r8
!>          tl_cff=-cff*cff*(tl_z_r(i,j,k+1)+tl_z_r(i-1,j,k+1)-         &
!>   &                       tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  ))
!>
            adfac=-cff*cff*ad_cff
            ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-adfac
            ad_z_r(i  ,j,k  )=ad_z_r(i  ,j,k  )-adfac
            ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)+adfac
            ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+adfac
            ad_cff=0.0_r8
          END DO
        END DO
      END DO J_LOOP2

# ifndef TS_FIXED
!
!=======================================================================
!  Adjoint tracer equation(s).
!=======================================================================
!
!-----------------------------------------------------------------------
!  Start computation of tracers at n+1 time-step, t(i,j,k,nnew,itrc).
!-----------------------------------------------------------------------
!
!  Compute vertical diffusive fluxes "FC" of the tracer fields at
!  current time step n, and at horizontal RHO-points and vertical
!  W-points.
!
      DO j=Jstr,Jend
        DO itrc=1,NT(ng)
!
!  Compute new tracer field (m Tunits).
!
          DO k=1,N(ng)
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_TS
!!            DiaTwrk(i,j,k,itrc,iTvdif)=cff2
!!            DiaTwrk(i,j,k,itrc,iTrate)=cff1
#  endif
!>            tl_t(i,j,k,nnew,itrc)=tl_cff1+tl_cff2
!>
              ad_cff1=ad_cff1+ad_t(i,j,k,nnew,itrc)
              ad_cff2=ad_cff2+ad_t(i,j,k,nnew,itrc)
              ad_t(i,j,k,nnew,itrc)=0.0_r8
!>            tl_cff2=tl_FC(i,k)-tl_FC(i,k-1)
!>
              ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff2
              ad_FC(i,k  )=ad_FC(i,k  )+ad_cff2
              ad_cff2=0.0_r8
!>            tl_cff1=tl_Hz(i,j,k)*t(i,j,k,nstp,itrc)+                  &
!>   &                Hz(i,j,k)*tl_t(i,j,k,nstp,itrc)
!>
              ad_t(i,j,k,nstp,itrc)=ad_t(i,j,k,nstp,itrc)+              &
     &                              Hz(i,j,k)*ad_cff1
              ad_Hz(i,j,k)=ad_Hz(i,j,k)+t(i,j,k,nstp,itrc)*ad_cff1
              ad_cff1=0.0_r8
            END DO
          END DO
!
!  Apply bottom and surface tracer flux conditions.
!
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=dt(ng)*tl_stflx(i,j,itrc)
!>
            ad_stflx(i,j,itrc)=ad_stflx(i,j,itrc)+dt(ng)*ad_FC(i,N(ng))
            ad_FC(i,N(ng))=0.0_r8
!>          tl_FC(i,0)=dt(ng)*tl_btflx(i,j,itrc)
!>
            ad_btflx(i,j,itrc)=ad_btflx(i,j,itrc)+dt(ng)*ad_FC(i,0)
            ad_FC(i,0)=0.0_r8
          END DO

#  ifdef SOLAR_SOURCE
!
!  Add in incoming solar radiation at interior W-points using decay
!  decay penetration function based on Jerlov water type.
!
          IF (itrc.eq.itemp) THEN
            DO k=1,N(ng)-1
              DO i=Istr,Iend
!>              tl_FC(i,k)=tl_FC(i,k)+dt(ng)*srflx(i,j)*tl_swdk(i,j,k)
!>
                ad_swdk(i,j,k)=ad_swdk(i,j,k)+                          &
     &                         dt(ng)*srflx(i,j)*ad_FC(i,k)
              END DO
            END DO
          END IF
#  endif
#  ifdef LMD_NONLOCAL_NOT_YET
!
!  Add in the nonlocal transport flux for unstable (convective)
!  forcing conditions into matrix FC when using the Large et al.
!  KPP scheme. The nonlocal transport is only applied to active
!  tracers.
!
          IF (itrc.le.NAT) THEN
            DO k=1,N(ng)-1
              DO i=Istr,Iend
!>              tl_FC(i,k)=tl_FC(i,k)-                                  &
!>   &                     dt(ng)*(tl_Akt(i,j,k,itrc)*                  &
!>   &                             ghats(i,j,k,itrc)+                   &
!>   &                             Akt(i,j,k,itrc)*                     &
!>   &                             tl_ghats(i,j,k,itrc))
!>
                adfac=dt(ng)*ad_FC(i,k)
                ad_ghats(i,j,k,itrc)=ad_ghats(i,j,k,itrc)-              &
     &                               Akt(i,j,k,itrc)*adfac
                ad_Akt(i,j,k,itrc)=ad_Akt(i,j,k,itrc)-                  &
     &                             ghats(i,j,k,itrc)*adfac
              END DO
            END DO
          END IF
#  endif
!
!  Compute adjoint vertical diffusive fluxes "FC" of the tracer fields
!  at current time step n, and at horizontal RHO-points and vertical
!  W-points. Notice that the vertical diffusion coefficients for
!  passive tracers is the same as that for salinity (ltrc=NAT).
!
          cff3=dt(ng)*(1.0_r8-lambda)
          ltrc=MIN(NAT,itrc)
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
!>            tl_FC(i,k)=cff3*                                          &
!>   &                   (cff*(tl_Akt(i,j,k,ltrc)*                      &
!>   &                         (t(i,j,k+1,nstp,itrc)-                   &
!>   &                          t(i,j,k  ,nstp,itrc))+                  &
!>   &                         Akt(i,j,k,ltrc)*                         &
!>   &                         (tl_t(i,j,k+1,nstp,itrc)-                &
!>   &                          tl_t(i,j,k  ,nstp,itrc)))+              &
!>   &                    tl_cff*(Akt(i,j,k,ltrc)*                      &
!>   &                            (t(i,j,k+1,nstp,itrc)-                &
!>   &                             t(i,j,k,nstp,itrc))))
!>
              adfac=cff3*ad_FC(i,k)
              adfac1=adfac*cff
              adfac2=adfac1*Akt(i,j,k,ltrc)
              ad_Akt(i,j,k,ltrc)=ad_Akt(i,j,k,ltrc)+                    &
     &                           (t(i,j,k+1,nstp,itrc)-                 &
     &                            t(i,j,k  ,nstp,itrc))*adfac1
              ad_t(i,j,k  ,nstp,itrc)=ad_t(i,j,k  ,nstp,itrc)-adfac2
              ad_t(i,j,k+1,nstp,itrc)=ad_t(i,j,k+1,nstp,itrc)+adfac2
              ad_cff=ad_cff+                                            &
     &               (Akt(i,j,k,ltrc)*                                  &
     &                (t(i,j,k+1,nstp,itrc)-                            &
     &                 t(i,j,k  ,nstp,itrc)))*adfac
              ad_FC(i,k)=0.0_r8
!>            tl_cff=-cff*cff*(tl_z_r(i,j,k+1)-tl_z_r(i,j,k))
!>
              adfac=-cff*cff*ad_cff
              ad_z_r(i,j,k  )=ad_z_r(i,j,k  )-adfac
              ad_z_r(i,j,k+1)=ad_z_r(i,j,k+1)+adfac
              ad_cff=0.0_r8
            END DO
          END DO
        END DO
      END DO
!
!------------------------------------------------------------------------
!  Time-step adjoint vertical advection of tracers.  Impose artificial
!  continuity equation.
!------------------------------------------------------------------------
!
      J_LOOP1 : DO j=Jstr,Jend
!
!  First, compute BASIC STATE artificial continuity equation and BASIC
!  STATE tracer vertical advection flux needed for time-stepping.
!
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff=0.5_r8*dt(ng)
        ELSE
          cff=(1.0_r8-Gamma)*dt(ng)
        END IF
        DO k=1,N(ng)
          DO i=Istr,Iend
            DC(i,k)=1.0_r8/(Hz(i,j,k)-cff*pm(i,j)*pn(i,j)*              &
     &                      (Huon(i+1,j,k)-Huon(i,j,k)+                 &
     &                       Hvom(i,j+1,k)-Hvom(i,j,k)+                 &
     &                      (W(i,j,k)-W(i,j,k-1))))
          END DO
        END DO

        T_LOOP2 : DO itrc=1,NT(ng)

#  ifdef TS_SVADVECTION_TL
!
!  Build BASIC conservative parabolic splines for the vertical
!  derivatives "FC" of the tracer.  Then, the interfacial "FC" values
!  are converted to vertical advective flux.
!
          DO i=Istr,Iend
#   ifdef NEUMANN
            FC(i,0)=1.5_r8*t(i,j,1,nstp,itrc)
            CF(i,1)=0.5_r8
#   else
            FC(i,0)=2.0_r8*t(i,j,1,nstp,itrc)
            CF(i,1)=1.0_r8
#   endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              FC(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,nstp,itrc)+    &
     &                             Hz(i,j,k+1)*t(i,j,k  ,nstp,itrc))-   &
     &                     Hz(i,j,k+1)*FC(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
#   ifdef NEUMANN
            FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),nstp,itrc)-                 &
     &                   FC(i,N(ng)-1))/(2.0_r8-CF(i,N(ng)))
#   else
            FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),nstp,itrc)-                 &
     &                   FC(i,N(ng)-1))/(1.0_r8-CF(i,N(ng)))
#   endif
          END DO
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
              FC(i,k+1)=W(i,j,k+1)*FC(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO

#  elif defined TS_A4VADVECTION_TL
!
!  Fourth-order, BASIC STATE Akima vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=t(i,j,k+1,nstp,itrc)-                             &
     &                t(i,j,k  ,nstp,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=FC(i,1)
            FC(i,N(ng))=FC(i,N(ng)-1)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=2.0_r8*FC(i,k)*FC(i,k-1)
              IF (cff.gt.eps) THEN
                CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
              ELSE
                CF(i,k)=0.0_r8
              END IF
            END DO
          END DO
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,nstp,itrc)+                     &
     &                        t(i,j,k+1,nstp,itrc)-                     &
     &                        cff1*(CF(i,k+1)-CF(i,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO

#  elif defined TS_C2VADVECTION_TL
!
!  Second-order, BASIC STATE central differences vertical advective
!  flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,nstp,itrc)+                     &
     &                        t(i,j,k+1,nstp,itrc))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO

#  elif defined TS_MPDATA_NOT_YET
!
!  First_order, BASIC STATE upstream differences vertical advective
!  flux.
!
          DO i=Istr,Iend
            DO k=1,N(ng)-1
              cff1=MAX(W(i,j,k),0.0_r8)
              cff2=MIN(W(i,j,k),0.0_r8)
              FC(i,k)=cff1*t(i,j,k  ,nstp,itrc)+                        &
     &                cff2*t(i,j,k+1,nstp,itrc)
            END DO
#   ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,nstp,itrc)
#   else
            FC(i,0)=0.0_r8
#   endif
            FC(i,N(ng))=0.0_r8
          END DO

#  else
!
!  Fourth-order, BASIC STATE central differences vertical advective
!  flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                (cff2*(t(i,j,k  ,nstp,itrc)+                      &
     &                       t(i,j,k+1,nstp,itrc))-                     &
     &                 cff3*(t(i,j,k-1,nstp,itrc)+                      &
     &                       t(i,j,k+2,nstp,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,1)=W(i,j,1)*                                           &
     &              (cff1*t(i,j,1,nstp,itrc)+                           &
     &               cff2*t(i,j,2,nstp,itrc)-                           &
     &               cff3*t(i,j,3,nstp,itrc))
            FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                               &
     &                    (cff1*t(i,j,N(ng)  ,nstp,itrc)+               &
     &                     cff2*t(i,j,N(ng)-1,nstp,itrc)-               &
     &                     cff3*t(i,j,N(ng)-2,nstp,itrc))
            FC(i,N(ng))=0.0_r8
          END DO
#  endif
!
! Time-step vertical advection of tracers (Tunits).
!
! WARNING: t(:,:,:,3,itrc) at this point should be in units of
! =======  m Tunits. But, t(:,:,:,3,itrc) is read from a BASIC
!          STATE file and is in Tunits, so we need to multiply
!          by level thickness (Hz).
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.5_r8*dt(ng)
          ELSE
            cff=(1.0_r8-Gamma)*dt(ng)
          END IF
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=cff*pm(i,j)*pn(i,j)
!>            tl_t(i,j,k,3,itrc)=tl_DC(i,k)*                            &
!>   &                           (t(i,j,k,3,itrc)*Hz(i,j,k)-            &
!>   &                            cff1*(FC(i,k)-FC(i,k-1)))+            &
!>   &                           DC(i,k)*                               &
!>   &                           (tl_t(i,j,k,3,itrc)-                   &
!>   &                            cff1*(tl_FC(i,k)-tl_FC(i,k-1)))
!>
              adfac=DC(i,k)*ad_t(i,j,k,3,itrc)
              adfac1=adfac*cff1
              ad_DC(i,k)=ad_DC(i,k)+                                    &
     &                   (t(i,j,k,3,itrc)*Hz(i,j,k)-                    &
     &                    cff1*(FC(i,k)-FC(i,k-1)))*                    &
     &                    ad_t(i,j,k,3,itrc)
              ad_FC(i,k-1)=ad_FC(i,k-1)+adfac1
              ad_FC(i,k  )=ad_FC(i,k  )-adfac1
              ad_t(i,j,k,3,itrc)=adfac
            END DO
          END DO
!
!  Apply mass point sources (volume vertical influx), if any.
!  (Notice the j-loop is pipelined above).
!
          IF (LwSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (LtracerSrc(itrc,ng).and.                              &
     &            ((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.            &
     &            (j.eq.Jsrc)) THEN
                DO k=1,N(ng)-1
!>                tl_FC(Isrc,k)=tl_FC(Isrc,k)+0.0_r8
!>
                END DO
              END IF
            END DO
          END IF
!
!-----------------------------------------------------------------------
!  Compute adjoint time rate of change of intermediate tracer due to
!  vertical advection.
!-----------------------------------------------------------------------

#  ifdef TS_SVADVECTION_TL
!
!  Construct BASIC STATE conservative parabolic splines for the vertical
!  derivatives "CF" of the tracer.
!
          DO i=Istr,Iend
#   ifdef NEUMANN
            FC(i,0)=1.5_r8*t(i,j,1,nstp,itrc)
            CF(i,1)=0.5_r8
#   else
            FC(i,0)=2.0_r8*t(i,j,1,nstp,itrc)
            CF(i,1)=1.0_r8
#   endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              FC(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,nstp,itrc)+    &
     &                             Hz(i,j,k+1)*t(i,j,k  ,nstp,itrc))-   &
     &                     Hz(i,j,k+1)*FC(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
#   ifdef NEUMANN
            FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),nstp,itrc)-                 &
     &                   FC(i,N(ng)-1))/(2.0_r8-CF(i,N(ng)))
#   else
            FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),nstp,itrc)-                 &
     &                   FC(i,N(ng)-1))/(1.0_r8-CF(i,N(ng)))
#   endif
          END DO
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
            END DO
          END DO
!
!  Now the adjoint spline code.
!
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
          END DO
!
!   Adjoint back substitution.
!
          DO k=0,N(ng)-1
            DO i=Istr,Iend
!>            tl_FC(i,k+1)=tl_W(i,j,k+1)*FC(i,k+1)+                     &
!>   &                     W(i,j,k+1)*tl_FC(i,k+1)
!>
              ad_W(i,j,k+1)=ad_W(i,j,k+1)+FC(i,k+1)*ad_FC(i,k+1)
              ad_FC(i,k+1)=W(i,j,k+1)*ad_FC(i,k+1)
!>            tl_FC(i,k)=tl_FC(i,k)-CF(i,k+1)*tl_FC(i,k+1)
!>
              ad_FC(i,k+1)=ad_FC(i,k+1)-CF(i,k+1)*ad_FC(i,k)
            END DO
          END DO
!
          DO i=Istr,Iend
#   ifdef NEUMANN
!>          tl_FC(i,N(ng))=(3.0_r8*tl_t(i,j,N(ng),nstp,itrc)-           &
!>     &                    tl_FC(i,N(ng)-1))/                          &
!>     &                   (2.0_r8-CF(i,N(ng)))
!>
            adfac=ad_FC(i,N(ng))/(2.0_r8-CF(i,N(ng)))
            ad_t(i,j,N(ng),nstp,itrc)=ad_t(i,j,N(ng),nstp,itrc)+        &
     &                                3.0_r8*adfac
            ad_FC(i,N(ng)-1)=ad_FC(i,N(ng)-1)-adfac
            ad_FC(i,N(ng))=0.0_r8
#   else
!>          tl_FC(i,N(ng))=(2.0_r8*tl_t(i,j,N(ng),nstp,itrc)-           &
!>   &                      tl_FC(i,N(ng)-1))/                          &
!>   &                     (1.0_r8-CF(i,N(ng)))
!>
            adfac=ad_FC(i,N(ng))/(1.0_r8-CF(i,N(ng)))
            ad_t(i,j,N(ng),nstp,itrc)=ad_t(i,j,N(ng),nstp,itrc)+        &
     &                                2.0_r8*adfac
            ad_FC(i,N(ng)-1))=ad_FC(i,N(ng)-1))-adfac
            ad_FC(i,N(ng))=0.0_r8
#   endif
          END DO
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
!>            tl_FC(i,k)=cff*                                           &
!>   &                   (3.0_r8*(Hz(i,j,k  )*tl_t(i,j,k+1,nstp,itrc)+  &
!>   &                            Hz(i,j,k+1)*tl_t(i,j,k  ,nstp,itrc)+  &
!>   &                            tl_Hz(i,j,k  )*t(i,j,k+1,nstp,itrc)+  &
!>   &                            tl_Hz(i,j,k+1)*t(i,j,k  ,nstp,itrc))- &
!>   &                    (tl_Hz(i,j,k+1)*FC(i,k-1)+                    &
!>   &                     2.0_r8*(tl_Hz(i,j,k  )+                      &
!>   &                             tl_Hz(i,j,k+1))*FC(i,k)+             &
!>   &                     tl_Hz(i,j,k  )*FC(i,k+1))-                   &
!>   &                    Hz(i,j,k+1)*tl_FC(i,k-1))
!>
              adfac=cff*ad_FC(i,k)
              adfac1=3.0_r8*adfac
              adfac2=2.0_r8*adfac
              ad_t(i,j,k  ,nstp,itrc)=ad_t(i,j,k  ,nstp,itrc)+          &
     &                                Hz(i,j,k+1)*adfac1
              ad_t(i,j,k+1,nstp,itrc)=ad_t(i,j,k+1,nstp,itrc)+          &
     &                                Hz(i,j,k  )*adfac1
              ad_Hz(i,j,k  )=ad_Hz(i,j,k  )+                            &
     &                       t(i,j,k+1,nstp,itrc)*adfac1-               &
     &                       FC(i,k  )*adfac2-                          &
     &                       FC(i,k+1)*adfac
              ad_Hz(i,j,k+1)=ad_Hz(i,j,k+1)+                            &
     &                       t(i,j,k  ,nstp,itrc)*adfac1-               &
     &                       FC(i,k-1)*adfac-                           &
     &                       FC(i,k  )*adfac2
              ad_FC(i,k-1)=ad_FC(i,k-1)-Hz(i,j,k+1)*adfac
              ad_FC(i,k)=0.0_r8
            END DO
          END DO
!
          DO i=Istr,Iend
#   ifdef NEUMANN
!>          tl_FC(i,0)=1.5_r8*tl_t(i,j,1,nstp,itrc)
!>
            ad_t(i,j,1,nstp,itrc)=ad_t(i,j,1,nstp,itrc)+                &
     &                            1.5_r8*ad_FC(i,0)
            ad_FC(i,0)=0.0_r8
#   else
!>          tl_FC(i,0)=2.0_r8*tl_t(i,j,1,nstp,itrc)
!>
            ad_t(i,j,1,nstp,itrc)=ad_t(i,j,1,nstp,itrc)+                &
     &                            2.0_r8*ad_FC(i,0)
            ad_FC(i,0)=0.0_r8
#   endif
          END DO

#  elif defined TS_A4VADVECTION_TL
!
!  Fourth-order, adjoint Akima vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=t(i,j,k+1,nstp,itrc)-                             &
     &                t(i,j,k  ,nstp,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=FC(i,1)
            FC(i,N(ng))=FC(i,N(ng)-1)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=2.0_r8*FC(i,k)*FC(i,k-1)
              IF (cff.gt.eps) THEN
                CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
              ELSE
                CF(i,k)=0.0_r8
              END IF
            END DO
          END DO
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
#   ifdef SED_MORPH
!>          tl_FC(i,0)=tl_W(i,j,0)*t(i,j,1,nstp,itrc)+                  &
!>   &                 W(i,j,0)*tl_t(i,j,1,nstp,itrc)
!>
            ad_W(i,j,0)=ad_W(i,j,0)+                                    &
     &                  t(i,j,1,nstp,itrc)*ad_FC(i,0)
            ad_t(i,j,1,nstp,itrc)=ad_t(i,j,1,nstp,itrc)+                &
     &                            W(i,j,0)*ad_FC(i,0)
            ad_FC(i,0)=0.0_r8
#   else
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
#   endif
          END DO
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
!>            tl_FC(i,k)=0.5_r8*                                        &
!>   &                   (tl_W(i,j,k)*                                  &
!>   &                    (t(i,j,k  ,nstp,itrc)+                        &
!>   &                     t(i,j,k+1,nstp,itrc)-                        &
!>   &                     cff1*(CF(i,k+1)-CF(i,k)))+                   &
!>   &                    W(i,j,k)*                                     &
!>   &                    (tl_t(i,j,k  ,nstp,itrc)+                     &
!>   &                     tl_t(i,j,k+1,nstp,itrc)-                     &
!>   &                     cff1*(tl_CF(i,k+1)-tl_CF(i,k))))
!>
              adfac=0.5_r8*ad_FC(i,k)
              adfac1=adfac*W(i,j,k)
              adfac2=adfac1*cff1
              ad_CF(i,k  )=ad_CF(i,k  )+adfac2
              ad_CF(i,k+1)=ad_CF(i,k+1)-adfac2
              ad_t(i,j,k  ,nstp,itrc)=ad_t(i,j,k  ,nstp,itrc)+adfac1
              ad_t(i,j,k+1,nstp,itrc)=ad_t(i,j,k+1,nstp,itrc)+adfac1
              ad_W(i,j,k)=ad_W(i,j,k)+                                  &
     &                    (t(i,j,k  ,nstp,itrc)+                        &
     &                     t(i,j,k+1,nstp,itrc)-                        &
     &                     cff1*(CF(i,k+1)-CF(i,k)))*adfac
              ad_FC(i,k)=0.0_r8
            END DO
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=2.0_r8*FC(i,k)*FC(i,k-1)
              IF (cff.gt.eps) THEN
!>              tl_CF(i,k)=((FC(i,k)+FC(i,k-1))*tl_cff-                 &
!>   &                      cff*(tl_FC(i,k)+tl_FC(i,k-1)))/             &
!>   &                     ((FC(i,k)+FC(i,k-1))*(FC(i,k)+FC(i,k-1)))
!>
                adfac=ad_CF(i,k)/                                       &
     &                ((FC(i,k)+FC(i,k-1))*(FC(i,k)+FC(i,k-1)))
                adfac1=adfac*cff
                ad_FC(i,k-1)=ad_FC(i,k-1)-adfac1
                ad_FC(i,k  )=ad_FC(i,k  )-adfac1
                ad_cff=ad_cff+(FC(i,k)+FC(i,k-1))*adfac
                ad_CF(i,k)=0.0_r8
              ELSE
!>              tl_CF(i,k)=0.0_r8
!>
                ad_CF(i,k)=0.0_r8
              END IF
!>            tl_cff=2.0_r8*(tl_FC(i,k)*FC(i,k-1)+                      &
!>   &                       FC(i,k)*tl_FC(i,k-1))
!>
              adfac=2.0_r8*ad_cff
              ad_FC(i,k-1)=ad_FC(i,k-1)+FC(i,k  )*adfac
              ad_FC(i,k  )=ad_FC(i,k  )+FC(i,k-1)*adfac
              ad_cff=0.0_r8
            END DO
          END DO
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=tl_FC(i,N(ng)-1)
!>
            ad_FC(i,N(ng)-1)=ad_FC(i,N(ng)-1)+ad_FC(i,N(ng))
            ad_FC(i,N(ng))=0.0_r8
!>          tl_FC(i,0)=tl_FC(i,1)
!>
            ad_FC(i,1)=ad_FC(i,1)+ad_FC(i,0)
            ad_FC(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
!>            tl_FC(i,k)=tl_t(i,j,k+1,nstp,itrc)-                       &
!>   &                   tl_t(i,j,k  ,nstp,itrc)
!>
              ad_t(i,j,k  ,nstp,itrc)=ad_t(i,j,k  ,nstp,itrc)-ad_FC(i,k)
              ad_t(i,j,k+1,nstp,itrc)=ad_t(i,j,k+1,nstp,itrc)+ad_FC(i,k)
              ad_FC(i,k)=0.0_r8
            END DO
          END DO

#  elif defined TS_C2VADVECTION_TL
!
!  Second-order, central differences vertical advective flux.
!
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
#   ifdef SED_MORPH
!>          tl_FC(i,0)=tl_W(i,j,0)*t(i,j,1,nstp,itrc)                   &
!>   &                 W(i,j,0)*tl_t(i,j,1,nstp,itrc)
!>
            ad_W(i,j,0)=ad_W(i,j,0)+                                    &
     &                  t(i,j,1,nstp,itrc)*ad_FC(i,0)
            ad_t(i,j,1,nstp,itrc)=ad_t(i,j,1,nstp,itrc)+                &
     &                            W(i,j,0)*ad_FC(i,0)
            ad_FC(i,0)=0.0_r8
#   else
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
#   endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
!>            tl_FC(i,k)=0.5_r8*                                        &
!>   &                   (tl_W(i,j,k)*                                  &
!>   &                    (t(i,j,k  ,nstp,itrc)+                        &
!>   &                     t(i,j,k+1,nstp,itrc))+                       &
!>   &                    W(i,j,k)*                                     &
!>   &                    (tl_t(i,j,k  ,nstp,itrc)+                     &
!>   &                     tl_t(i,j,k+1,nstp,itrc)))
!>
              adfac=0.5_r8*ad_FC(i,k)
              adfac1=adfac*W(i,j,k)
              ad_W(i,j,k)=ad_W(i,j,k)+                                  &
     &                    (t(i,j,k  ,nstp,itrc)+                        &
     &                     t(i,j,k+1,nstp,itrc))*adfac
              ad_t(i,j,k  ,nstp,itrc)=ad_t(i,j,k  ,nstp,itrc)+adfac1
              ad_t(i,j,k+1,nstp,itrc)=ad_t(i,j,k+1,nstp,itrc)+adfac1
              ad_FC(i,k)=0.0_r8
            END DO
          END DO

#  elif defined TS_MPDATA_NOT_YET
!
!  First_order, upstream differences vertical advective flux.
!
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
#   ifdef SED_MORPH
!>          tl_FC(i,0)=tl_W(i,j,0)*t(i,j,1,nstp,itrc)+                  &
!>   &                 W(i,j,0)*tl_t(i,j,1,nstp,itrc)
!>
            ad_t(i,j,1,nstp,itrc)=ad_t(i,j,1,nstp,itrc)+                &
     &                            W(i,j,0)*ad_FC(i,0)
            ad_W(i,j,0)=ad_W(i,j,0)+                                    &
     &                  t(i,j,1,nstp,itrc)+ad_FC(i,0)
            ad_FC(i,0)=0.0_r8
#   else
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
#   endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=MAX(W(i,j,k),0.0_r8)
              cff2=MIN(W(i,j,k),0.0_r8)
!>            tl_FC(i,k)=tl_cff1*t(i,j,k  ,nstp,itrc)+                  &
!>   &                   cff1*tl_t(i,j,k  ,nstp,itrc)+                  &
!>   &                   tl_cff2*t(i,j,k+1,nstp,itrc)+                  &
!>   &                   cff2*tl_t(i,j,k+1,nstp,itrc)
!>
              ad_t(i,j,k  ,nstp,itrc)=ad_t(i,j,k  ,nstp,itrc)+          &
     &                                cff1*ad_FC(i,k)
              ad_t(i,j,k+1,nstp,itrc)=ad_t(i,j,k+1,nstp,itrc)+          &
     &                                cff2*ad_FC(i,k)
              ad_cff1=ad_cff1+t(i,j,k  ,nstp,itrc)*ad_FC(i,k)
              ad_cff2=ad_cff2+t(i,j,k+1,nstp,itrc)*ad_FC(i,k)
              ad_FC(i,k)=0.0_r8
!>            tl_cff1=(0.5_r8+SIGN(0.5_r8, W(i,j,k)))*tl_W(i,j,k)
!>            tl_cff2=(0.5_r8+SIGN(0.5_r8,-W(i,j,k)))*tl_W(i,j,k)
!>
              ad_W(i,j,k)=ad_W(i,j,k)+                                  &
     &                    (0.5_r8+SIGN(0.5_r8,-W(i,j,k)))*ad_cff2+      &
     &                    (0.5_r8+SIGN(0.5_r8, W(i,j,k)))*ad_cff1
              ad_cff2=0.0_r8
              ad_cff1=0.0_r8
            END DO
          END DO

#  else
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
!>          tl_FC(i,N(ng)-1)=tl_W(i,j,N(ng)-1)*                         &
!>   &                       (cff1*t(i,j,N(ng)  ,nstp,itrc)+            &
!>   &                        cff2*t(i,j,N(ng)-1,nstp,itrc)-            &
!>   &                        cff3*t(i,j,N(ng)-2,nstp,itrc))+           &
!>   &                       W(i,j,N(ng)-1)*                            &
!>   &                       (cff1*tl_t(i,j,N(ng)  ,nstp,itrc)+         &
!>   &                        cff2*tl_t(i,j,N(ng)-1,nstp,itrc)-         &
!>   &                        cff3*tl_t(i,j,N(ng)-2,nstp,itrc))
!>
            adfac=W(i,j,N(ng)-1)*ad_FC(i,N(ng)-1)
            ad_W(i,j,N(ng)-1)=ad_W(i,j,N(ng)-1)+                        &
     &                        (cff1*t(i,j,N(ng)  ,nstp,itrc)+           &
     &                         cff2*t(i,j,N(ng)-1,nstp,itrc)-           &
     &                         cff3*t(i,j,N(ng)-2,nstp,itrc))*          &
     &                        ad_FC(i,N(ng)-1)
            ad_t(i,j,N(ng)-2,nstp,itrc)=ad_t(i,j,N(ng)-2,nstp,itrc)-    &
     &                                  cff3*adfac
            ad_t(i,j,N(ng)-1,nstp,itrc)=ad_t(i,j,N(ng)-1,nstp,itrc)+    &
     &                                  cff2*adfac
            ad_t(i,j,N(ng)  ,nstp,itrc)=ad_t(i,j,N(ng)  ,nstp,itrc)+    &
     &                                  cff1*adfac
            ad_FC(i,N(ng)-1)=0.0_r8
!>          tl_FC(i,1)=tl_W(i,j,1)*                                     &
!>   &                 (cff1*t(i,j,1,nstp,itrc)+                        &
!>   &                  cff2*t(i,j,2,nstp,itrc)-                        &
!>   &                  cff3*t(i,j,3,nstp,itrc))+                       &
!>   &                 W(i,j,1)*                                        &
!>   &                 (cff1*tl_t(i,j,1,nstp,itrc)+                     &
!>   &                  cff2*tl_t(i,j,2,nstp,itrc)-                     &
!>   &                  cff3*tl_t(i,j,3,nstp,itrc))
!>
            adfac=W(i,j,1)*ad_FC(i,1)
            ad_W(i,j,1)=ad_W(i,j,1)+                                    &
     &                  (cff1*t(i,j,1,nstp,itrc)+                       &
     &                   cff2*t(i,j,2,nstp,itrc)-                       &
     &                   cff3*t(i,j,3,nstp,itrc))*ad_FC(i,1)
            ad_t(i,j,1,nstp,itrc)=ad_t(i,j,1,nstp,itrc)+cff1*adfac
            ad_t(i,j,2,nstp,itrc)=ad_t(i,j,2,nstp,itrc)+cff2*adfac
            ad_t(i,j,3,nstp,itrc)=ad_t(i,j,3,nstp,itrc)-cff3*adfac
            ad_FC(i,1)=0.0_r8
#   ifdef SED_MORPH
!>          tl_FC(i,0)=2.0_r8*(tl_W(i,j,0)*                             &
!>   &                         (cff2*t(i,j,1,nstp,itrc)-                &
!>   &                          cff3*t(i,j,2,nstp,itrc))+               &
!>   &                         W(i,j,0)*                                &
!>   &                         (cff2*tl_t(i,j,1,nstp,itrc)-             &
!>   &                          cff3*tl_t(i,j,2,nstp,itrc)))
!>
            adfac=2.0_r8*ad_FC(i,0)
            adfac1=adfac*W(i,j,0)
            ad_W(i,j,0)=ad_W(i,j,0)+                                    &
     &                  (cff2*t(i,j,1,nstp,itrc)-                       &
     &                   cff3*t(i,j,2,nstp,itrc))*adfac
            ad_t(i,j,1,nstp,itrc)=ad_t(i,j,1,nstp,itrc)+cff2*adfac1
            ad_t(i,j,2,nstp,itrc)=ad_t(i,j,2,nstp,itrc)-cff3*adfac1
            ad_FC(i,0)=0.0_r8
#   else
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
#   endif
          END DO
          DO k=2,N(ng)-2
            DO i=Istr,Iend
!>            tl_FC(i,k)=tl_W(i,j,k)*                                   &
!>   &                   (cff2*(t(i,j,k  ,nstp,itrc)+                   &
!>   &                          t(i,j,k+1,nstp,itrc))-                  &
!>   &                    cff3*(t(i,j,k-1,nstp,itrc)+                   &
!>   &                          t(i,j,k+2,nstp,itrc)))+                 &
!>   &                   W(i,j,k)*                                      &
!>   &                   (cff2*(tl_t(i,j,k  ,nstp,itrc)+                &
!>   &                          tl_t(i,j,k+1,nstp,itrc))-               &
!>   &                    cff3*(tl_t(i,j,k-1,nstp,itrc)+                &
!>   &                          tl_t(i,j,k+2,nstp,itrc)))
!>
              adfac=W(i,j,k)*ad_FC(i,k)
              adfac1=adfac*cff2
              adfac2=adfac*cff3
              ad_W(i,j,k)=ad_W(i,j,k)+                                  &
     &                    (cff2*(t(i,j,k  ,nstp,itrc)+                  &
     &                           t(i,j,k+1,nstp,itrc))-                 &
     &                     cff3*(t(i,j,k-1,nstp,itrc)+                  &
     &                           t(i,j,k+2,nstp,itrc)))*ad_FC(i,k)
              ad_t(i,j,k-1,nstp,itrc)=ad_t(i,j,k-1,nstp,itrc)-adfac2
              ad_t(i,j,k  ,nstp,itrc)=ad_t(i,j,k  ,nstp,itrc)+adfac1
              ad_t(i,j,k+1,nstp,itrc)=ad_t(i,j,k+1,nstp,itrc)+adfac1
              ad_t(i,j,k+2,nstp,itrc)=ad_t(i,j,k+2,nstp,itrc)-adfac2
              ad_FC(i,k)=0.0_r8
            END DO
          END DO
#  endif
        END DO T_LOOP2
!
!-----------------------------------------------------------------------
!  Compute intermediate adjoint tracer at n+1/2 time-step,
!  t(i,j,k,3,itrc).
!-----------------------------------------------------------------------
!
!  Compute adjoint artificial continuity equation (same for all tracers)
!  and load it into private array DC (1/m). Notice pipelined J-loop.
!
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff=0.5_r8*dt(ng)
        ELSE
          cff=(1.0_r8-Gamma)*dt(ng)
        END IF
        DO k=1,N(ng)
          DO i=Istr,Iend
            DC(i,k)=1.0_r8/(Hz(i,j,k)-cff*pm(i,j)*pn(i,j)*              &
     &                      (Huon(i+1,j,k)-Huon(i,j,k)+                 &
     &                       Hvom(i,j+1,k)-Hvom(i,j,k)+                 &
     &                      (W(i,j,k)-W(i,j,k-1))))
!>          tl_DC(i,k)=-DC(i,k)*DC(i,k)*                                &
!>   &                  (tl_Hz(i,j,k)-                                  &
!>   &                   cff*pm(i,j)*pn(i,j)*                           &
!>   &                   (tl_Huon(i+1,j,k)-tl_Huon(i,j,k)+              &
!>   &                    tl_Hvom(i,j+1,k)-tl_Hvom(i,j,k)+              &
!>   &                   (tl_W(i,j,k)-tl_W(i,j,k-1))))
!>
            adfac=-DC(i,k)*DC(i,k)*ad_DC(i,k)
            adfac1=adfac*cff*pm(i,j)*pn(i,j)
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+adfac
            ad_Huon(i  ,j,k)=ad_Huon(i  ,j,k)+adfac1
            ad_Huon(i+1,j,k)=ad_Huon(i+1,j,k)-adfac1
            ad_Hvom(i,j  ,k)=ad_Hvom(i,j  ,k)+adfac1
            ad_Hvom(i,j+1,k)=ad_Hvom(i,j+1,k)-adfac1
            ad_W(i,j,k-1)=ad_W(i,j,k-1)+adfac1
            ad_W(i,j,k  )=ad_W(i,j,k  )-adfac1
            ad_DC(i,k)=0.0_r8
          END DO
        END DO
      END DO J_LOOP1

#  if (defined AGE_MEAN && defined T_PASSIVE) && !defined TS_MPDATA
!
!  If inert passive tracer and Mean Age, compute age concentration (even
!  inert index) forced by the right-hand-side term that is concentration
!  of an associated conservative passive tracer (odd inert index).
!  Implemented and tested by W.G. Zhang and J. Wilkin. (m Tunits)
!
      IF (iic(ng).eq.ntfirst(ng)) THEN
        cff=0.5_r8*dt(ng)
      ELSE
        cff=(1.0_r8-Gamma)*dt(ng)
      END IF
!
      DO itrc=1,NPT,2
        iage=inert(itrc+1)                     ! even inert tracer index
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
!>            tl_t(i,j,k,3,iage)=tl_t(i,j,k,3,iage)+                    &
!>   &                           cff*                                   &
!>   &                           (Hz(i,j,k)*                            &
!>   &                            tl_t(i,j,k,nnew,inert(itrc))+         &
!>   &                            tl_Hz(i,j,k)*                         &
!>   &                            t(i,j,k,nnew,inert(itrc)))
!>
              adfac=cff*ad_t(i,j,k,3,iage)
              ad_t(i,j,k,nnew,inert(itrc))=ad_t(i,j,k,nnew,inert(itrc))+&
     &                                     Hz(i,j,k)*adfac
              ad_Hz(i,j,k)=ad_Hz(i,j,k)+                                &
     &                     t(i,j,k,nnew,inert(itrc))*adfac
            END DO
          END DO
        END DO
      END DO
#  endif
!
!  Compute adjoint time rate of change of intermediate tracer due to
!  horizontal advection.
!
      T_LOOP1 :DO itrc=1,NT(ng)
        K_LOOP: DO k=1,N(ng)
!
!  Time-step horizontal advection (m Tunits).
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.5_r8*dt(ng)
            cff1=1.0_r8
            cff2=0.0_r8
          ELSE
            cff=(1.0_r8-Gamma)*dt(ng)
            cff1=0.5_r8+Gamma
            cff2=0.5_r8-Gamma
          END IF
          DO j=Jstr,Jend
            DO i=Istr,Iend
!>            tl_t(i,j,k,3,itrc)=tl_Hz(i,j,k)*                          &
!>   &                           (cff1*t(i,j,k,nstp,itrc)+              &
!>   &                            cff2*t(i,j,k,nnew,itrc))+             &
!>   &                           Hz(i,j,k)*                             &
!>   &                           (cff1*tl_t(i,j,k,nstp,itrc)+           &
!>   &                            cff2*tl_t(i,j,k,nnew,itrc))-          &
!>   &                           cff*pm(i,j)*pn(i,j)*                   &
!>   &                           (tl_FX(i+1,j)-tl_FX(i,j)+              &
!>   &                            tl_FE(i,j+1)-tl_FE(i,j))
!>
              adfac1=Hz(i,j,k)*ad_t(i,j,k,3,itrc)
              adfac2=cff*pm(i,j)*pn(i,j)*ad_t(i,j,k,3,itrc)
              ad_Hz(i,j,k)=ad_Hz(i,j,k)+                                &
     &                     (cff1*t(i,j,k,nstp,itrc)+                    &
     &                      cff2*t(i,j,k,nnew,itrc))*ad_t(i,j,k,3,itrc)
              ad_t(i,j,k,nstp,itrc)=ad_t(i,j,k,nstp,itrc)+cff1*adfac1
              ad_t(i,j,k,nnew,itrc)=ad_t(i,j,k,nnew,itrc)+cff2*adfac1
              ad_FE(i,j  )=ad_FE(i,j  )+adfac2
              ad_FE(i,j+1)=ad_FE(i,j+1)-adfac2
              ad_FX(i  ,j)=ad_FX(i  ,j)+adfac2
              ad_FX(i+1,j)=ad_FX(i+1,j)-adfac2
              ad_t(i,j,k,3,itrc)=0.0_r8
            END DO
          END DO
!
!  Apply tracers point sources to the horizontal advection terms,
!  if any.
!
          IF (LuvSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (LtracerSrc(itrc,ng).and.                              &
     &            ((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.            &
     &            ((Jstr.le.Jsrc).and.(Jsrc.le.Jend+1))) THEN
                IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
!>                tl_FX(Isrc,Jsrc)=tl_Huon(Isrc,Jsrc,k)*                &
!>   &                             SOURCES(ng)%Tsrc(is,k,itrc)
!>
                  ad_Huon(Isrc,Jsrc,k)=ad_Huon(Isrc,Jsrc,k)+            &
     &                                 SOURCES(ng)%Tsrc(is,k,itrc)*     &
     &                                 ad_FX(Isrc,Jsrc)
                  ad_FX(Isrc,Jsrc)=0.0_r8
                ELSE
!>                tl_FE(Isrc,Jsrc)=tl_Hvom(Isrc,Jsrc,k)*                &
!>   &                             SOURCES(ng)%Tsrc(is,k,itrc)
!>
                  ad_Hvom(Isrc,Jsrc,k)=ad_Hvom(Isrc,Jsrc,k)+            &
     &                                 SOURCES(ng)%Tsrc(is,k,itrc)*     &
     &                                 ad_FE(Isrc,Jsrc)
                  ad_FE(Isrc,Jsrc)=0.0_r8
                END IF
              END IF
            END DO
          END IF

#  if defined TS_C2HADVECTION_TL
!
!  Second-order, centered differences horizontal advective fluxes.
!
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
!>            tl_FE(i,j)=0.5_r8*                                        &
!>   &                   (tl_Hvom(i,j,k)*                               &
!>   &                    (t(i,j-1,k,nstp,itrc)+                        &
!>   &                     t(i,j  ,k,nstp,itrc))+                       &
!>   &                    Hvom(i,j,k)*                                  &
!>   &                    (tl_t(i,j-1,k,nstp,itrc)+                     &
!>   &                     tl_t(i,j  ,k,nstp,itrc)))
!>
              adfac=0.5_r8*ad_FE(i,j)
              adfac1=adfac*Hvom(i,j,k)
              ad_Hvom(i,j,k)=ad_Hvom(i,j,k)+                            &
     &                       adfac*(t(i,j  ,k,nstp,itrc)+               &
     &                              t(i,j-1,k,nstp,itrc))
              ad_t(i,j-1,k,nstp,itrc)=ad_t(i,j-1,k,nstp,itrc)+adfac1
              ad_t(i,j  ,k,nstp,itrc)=ad_t(i,j  ,k,nstp,itrc)+adfac1
              ad_FE(i,j)=0.0_r8
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
!>            tl_FX(i,j)=0.5_r8*                                        &
!>   &                   (tl_Huon(i,j,k)*                               &
!>   &                    (t(i-1,j,k,nstp,itrc)+                        &
!>   &                     t(i  ,j,k,nstp,itrc))+                       &
!>   &                    Huon(i,j,k)*                                  &
!>   &                    (tl_t(i-1,j,k,nstp,itrc)+                     &
!>   &                     tl_t(i  ,j,k,nstp,itrc)))
!>
              adfac=0.5_r8*ad_FX(i,j)
              adfac1=adfac*Huon(i,j,k)
              ad_Huon(i,j,k)=ad_Huon(i,j,k)+                            &
     &                       adfac*(t(i-1,j,k,nstp,itrc)+               &
     &                              t(i  ,j,k,nstp,itrc))
              ad_t(i-1,j,k,nstp,itrc)=ad_t(i-1,j,k,nstp,itrc)+adfac1
              ad_t(i  ,j,k,nstp,itrc)=ad_t(i  ,j,k,nstp,itrc)+adfac1
              ad_FX(i,j)=0.0_r8
            END DO
          END DO

#  elif defined TS_MPDATA_NOT_YET
!
!  First-order, upstream differences horizontal advective fluxes.
!
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff1=MAX(Hvom(i,j,k),0.0_r8)
              cff2=MIN(Hvom(i,j,k),0.0_r8)
!>            tl_FE(i,j)=tl_cff1*t(i,j-1,k,nstp,itrc)+                  &
!>   &                   cff1*tl_t(i,j-1,k,nstp,itrc)+                  &
!>   &                   tl_cff2*t(i,j  ,k,nstp,itrc)+                  &
!>   &                   cff2*tl_t(i,j  ,k,nstp,itrc)
!>
              ad_t(i,j-1,k,nstp,itrc)=ad_t(i,j-1,k,nstp,itrc)+          &
     &                                cff1*ad_FE(i,j)
              ad_t(i,j  ,k,nstp,itrc)=ad_t(i,j  ,k,nstp,itrc)+          &
     &                                cff2*ad_FE(i,j)
              ad_cff1=ad_cff1+t(i,j-1,k,nstp,itrc)*ad_FE(i,j)
              ad_cff2=ad_cff2+t(i,j  ,k,nstp,itrc)*ad_FE(i,j)
              ad_FE(i,j)=0.0_r8
!>            tl_cff1=(0.5_r8+SIGN(0.5_r8, Hvom(i,j,k)))*tl_Hvom(i,j,k)
!>            tl_cff2=(0.5_r8+SIGN(0.5_r8,-Hvom(i,j,k)))*tl_Hvom(i,j,k)
!>
              ad_Hvom(i,j,k)=ad_Hvom(i,j,k)+                            &
     &                       (0.5_r8+SIGN(0.5_r8,-Hvom(i,j,k)))*        &
     &                       ad_cff2+                                   &
     &                       (0.5_r8+SIGN(0.5_r8, Hvom(i,j,k)))*        &
     &                       ad_cff1
              ad_cff2=0.0_r8
              ad_cff1=0.0_r8
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff1=MAX(Huon(i,j,k),0.0_r8)
              cff2=MIN(Huon(i,j,k),0.0_r8)
!>            tl_FX(i,j)=tl_cff1*t(i-1,j,k,nstp,itrc)+                  &
!>   &                   cff1*tl_t(i-1,j,k,nstp,itrc)+                  &
!>   &                   tl_cff2*t(i  ,j,k,nstp,itrc)+                  &
!>   &                   cff2*tl_t(i  ,j,k,nstp,itrc)
!>
              ad_t(i  ,j,k,nstp,itrc)=ad_t(i  ,j,k,nstp,itrc)+          &
     &                                cff2*ad_FX(i,j)
              ad_t(i-1,j,k,nstp,itrc)=ad_t(i-1,j,k,nstp,itrc)+          &
     &                                cff1*ad_FX(i,j)
              ad_cff1=ad_cff1+t(i-1,j,k,nstp,itrc)*ad_FX(i,j)
              ad_cff2=ad_cff2+t(i  ,j,k,nstp,itrc)*ad_FX(i,j)
              ad_FX(i,j)=0.0_r8
!>            tl_cff1=(0.5_r8+SIGN(0.5_r8, Huon(i,j,k)))*tl_Huon(i,j,k)
!>            tl_cff2=(0.5_r8+SIGN(0.5_r8,-Huon(i,j,k)))*tl_Huon(i,j,k)
!>
              ad_Huon(i,j,k)=ad_Huon(i,j,k)+                            &
     &                       (0.5_r8+SIGN(0.5_r8,-Huon(i,j,k)))*        &
     &                       ad_cff2+                                   &
     &                       (0.5_r8+SIGN(0.5_r8, Huon(i,j,k)))*        &
     &                       ad_cff1
              ad_cff2=0.0_r8
              ad_cff1=0.0_r8
            END DO
          END DO

#  else
!
#   if defined TS_U3HADVECTION_TL
!  Third-order, uptream-biased horizontal advective fluxes.
#  elif defined TS_A4HADVECTION_TL
!  Fourth-order, Akima horizontal advective fluxes.
#   else
!  Fourth-order, centered differences horizontal advective fluxes.
#   endif
!  Compute BASIC STATE "curv" and "grad" scratch arrays.
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=t(i,j  ,k,nstp,itrc)-                             &
     &                t(i,j-1,k,nstp,itrc)
#   ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#   endif
            END DO
          END DO
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jstr-1)=FE(i,Jstr)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jend+2)=FE(i,Jend+1)
              END DO
            END IF
          END IF
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
#   if defined TS_U3HADVECTION_TL
              curv(i,j)=FE(i,j+1)-FE(i,j)
#   elif defined TS_A4HADVECTION_TL
              cff=2.0_r8*FE(i,j+1)*FE(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#   else
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
#   endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#   ifdef TS_U3HADVECTION_TL
!>            tl_FE(i,j)=0.5_r8*                                        &
!>   &                   (tl_Hvom(i,j,k)*                               &
!>   &                    (t(i,j-1,k,nstp,itrc)+                        &
!>   &                     t(i,j  ,k,nstp,itrc))+                       &
!>   &                    Hvom(i,j,k)*                                  &
!>   &                    (tl_t(i,j-1,k,nstp,itrc)+                     &
!>   &                     tl_t(i,j  ,k,nstp,itrc)))-                   &
!>   &                   cff1*                                          &
!>   &                   (tl_curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+       &
!>   &                    curv(i,j-1)*                                  &
!>   &                    (0.5_r8+SIGN(0.5_r8, Hvom(i,j,k)))*           &
!>   &                    tl_Hvom(i,j,k)+                               &
!>   &                    tl_curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8)+       &
!>   &                    curv(i,j  )*                                  &
!>   &                    (0.5_r8+SIGN(0.5_r8,-Hvom(i,j,k)))*           &
!>   &                    tl_Hvom(i,j,k))
!>
              adfac=0.5_r8*ad_FE(i,j)
              adfac1=adfac*Hvom(i,j,k)
              adfac2=cff1*ad_FE(i,j)
              ad_Hvom(i,j,k)=ad_Hvom(i,j,k)+                            &
     &                       (t(i,j-1,k,nstp,itrc)+                     &
     &                        t(i,j  ,k,nstp,itrc))*adfac-              &
     &                       (curv(i,j-1)*                              &
     &                        (0.5_r8+SIGN(0.5_r8, Hvom(i,j,k)))+       &
     &                        curv(i,j  )*                              &
     &                        (0.5_r8+SIGN(0.5_r8,-Hvom(i,j,k))))*      &
     &                       adfac2
              ad_t(i,j-1,k,nstp,itrc)=ad_t(i,j-1,k,nstp,itrc)+adfac1
              ad_t(i,j  ,k,nstp,itrc)=ad_t(i,j  ,k,nstp,itrc)+adfac1
              ad_curv(i,j-1)=ad_curv(i,j-1)-                            &
     &                       MAX(Hvom(i,j,k),0.0_r8)*adfac2
              ad_curv(i,j  )=ad_curv(i,j  )-                            &
     &                       MIN(Hvom(i,j,k),0.0_r8)*adfac2
              ad_FE(i,j)=0.0_r8
#   else
!>            tl_FE(i,j)=0.5_r8*                                        &
!>   &                   (tl_Hvom(i,j,k)*                               &
!>   &                    (t(i,j-1,k,nstp,itrc)+                        &
!>   &                     t(i,j  ,k,nstp,itrc)-                        &
!>   &                     cff2*(grad(i,j  )-                           &
!>   &                           grad(i,j-1)))+                         &
!>   &                    Hvom(i,j,k)*                                  &
!>   &                    (tl_t(i,j-1,k,nstp,itrc)+                     &
!>   &                     tl_t(i,j  ,k,nstp,itrc)-                     &
!>   &                     cff2*(tl_grad(i,j  )-                        &
!>   &                           tl_grad(i,j-1))))
!>
              adfac=0.5_r8*ad_FE(i,j)
              adfac1=adfac*Hvom(i,j,k)
              adfac2=adfac1*cff2
              ad_Hvom(i,j,k)=ad_Hvom(i,j,k)+                            &
     &                       adfac*(t(i,j-1,k,nstp,itrc)+               &
     &                              t(i,j  ,k,nstp,itrc)-               &
     &                              cff2*(grad(i,j  )-                  &
     &                                    grad(i,j-1)))
              ad_t(i,j-1,k,nstp,itrc)=ad_t(i,j-1,k,nstp,itrc)+adfac1
              ad_t(i,j  ,k,nstp,itrc)=ad_t(i,j  ,k,nstp,itrc)+adfac1
              ad_grad(i,j-1)=ad_grad(i,j-1)+adfac2
              ad_grad(i,j  )=ad_grad(i,j  )-adfac2
              ad_FE(i,j)=0.0_r8
#   endif
            END DO
          END DO
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
#   if defined TS_U3HADVECTION_TL
!>            tl_curv(i,j)=tl_FE(i,j+1)-tl_FE(i,j)
!>
              ad_FE(i,j  )=ad_FE(i,j  )-ad_curv(i,j)
              ad_FE(i,j+1)=ad_FE(i,j+1)+ad_curv(i,j)
              ad_curv(i,j)=0.0_r8
#   elif defined TS_A4HADVECTION_TL
              cff=2.0_r8*FE(i,j+1)*FE(i,j)
              IF (cff.gt.eps) THEN
!>              tl_grad(i,j)=((FE(i,j+1)+FE(i,j))*tl_cff-               &
!>   &                        cff*(tl_FE(i,j+1)+tl_FE(i,j)))/           &
!>   &                       ((FE(i,j+1)+FE(i,j))*(FE(i,j+1)+FE(i,j)))
!>
                adfac=ad_grad(i,j)/                                     &
     &                ((FE(i,j+1)+FE(i,j))*(FE(i,j+1)+FE(i,j)))
                adfac1=adfac*cff
                ad_FE(i,j  )=ad_FE(i,j  )-adfac1
                ad_FE(i,j+1)=ad_FE(i,j+1)-adfac1
                ad_cff=ad_cff+(FE(i,j+1)+FE(i,j))*adfac
                ad_grad(i,j)=0.0_r8
              ELSE
!>              tl_grad(i,j)=0.0_r8
!>
                ad_grad(i,j)=0.0_r8
              END IF
!>            tl_cff=2.0_r8*(tl_FE(i,j+1)*FE(i,j)+                      &
!>   &                       FE(i,j+1)*tl_FE(i,j))
!>
              adfac=2.0_r8*ad_cff
              ad_FE(i,j  )=ad_FE(i,j  )+FE(i,j+1)*adfac
              ad_FE(i,j+1)=ad_FE(i,j+1)+FE(i,j  )*adfac
              ad_cff=0.0_r8
#   else
!>            tl_grad(i,j)=0.5_r8*(tl_FE(i,j+1)+tl_FE(i,j))
!>
              adfac=0.5_r8*ad_grad(i,j)
              ad_FE(i,j  )=ad_FE(i,j  )+adfac
              ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
              ad_grad(i,j)=0.0_r8
#   endif
            END DO
          END DO
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO i=Istr,Iend
!>              tl_FE(i,Jend+2)=tl_FE(i,Jend+1)
!>
                ad_FE(i,Jend+1)=ad_FE(i,Jend+1)+ad_FE(i,Jend+2)
                ad_FE(i,Jend+2)=0.0_r8
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO i=Istr,Iend
!>              tl_FE(i,Jstr-1)=tl_FE(i,Jstr)
!>
                ad_FE(i,Jstr)=ad_FE(i,Jstr)+ad_FE(i,Jstr-1)
                ad_FE(i,Jstr-1)=0.0_r8
              END DO
            END IF
          END IF
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
#   ifdef MASKING
!>            tl_FE(i,j)=tl_FE(i,j)*vmask(i,j)
!>
              ad_FE(i,j)=ad_FE(i,j)*vmask(i,j)
#   endif
!>            tl_FE(i,j)=tl_t(i,j  ,k,nstp,itrc)-                       &
!>   &                   tl_t(i,j-1,k,nstp,itrc)
!>
              ad_t(i,j-1,k,nstp,itrc)=ad_t(i,j-1,k,nstp,itrc)-ad_FE(i,j)
              ad_t(i,j  ,k,nstp,itrc)=ad_t(i,j  ,k,nstp,itrc)+ad_FE(i,j)
              ad_FE(i,j)=0.0_r8
            END DO
          END DO
!
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=t(i  ,j,k,nstp,itrc)-                             &
     &                t(i-1,j,k,nstp,itrc)
#   ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#   endif
            END DO
          END DO
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Istr-1,j)=FX(Istr,j)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Iend+2,j)=FX(Iend+1,j)
              END DO
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
#   if defined TS_U3HADVECTION_TL
              curv(i,j)=FX(i+1,j)-FX(i,j)
#   elif defined TS_A4HADVECTION_TL
              cff=2.0_r8*FX(i+1,j)*FX(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#   else
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
#   endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#   ifdef TS_U3HADVECTION_TL
!>            tl_FX(i,j)=0.5_r8*                                        &
!>   &                   (tl_Huon(i,j,k)*                               &
!>   &                    (t(i-1,j,k,nstp,itrc)+                        &
!>   &                     t(i  ,j,k,nstp,itrc))+                       &
!>   &                    Huon(i,j,k)*                                  &
!>   &                    (tl_t(i-1,j,k,nstp,itrc)+                     &
!>   &                     tl_t(i  ,j,k,nstp,itrc)))-                   &
!>   &                   cff1*                                          &
!>   &                   (tl_curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+       &
!>   &                    curv(i-1,j)*                                  &
!>   &                    (0.5_r8+SIGN(0.5_r8, Huon(i,j,k)))*           &
!>   &                    tl_Huon(i,j,k)+                               &
!>   &                    tl_curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8)+       &
!>   &                    curv(i  ,j)*                                  &
!>   &                    (0.5_r8+SIGN(0.5_r8,-Huon(i,j,k)))*           &
!>   &                    tl_Huon(i,j,k))
!>
              adfac=0.5_r8*ad_FX(i,j)
              adfac1=adfac*Huon(i,j,k)
              adfac2=cff1*ad_FX(i,j)
              ad_Huon(i,j,k)=ad_Huon(i,j,k)+                            &
     &                       (t(i-1,j,k,nstp,itrc)+                     &
     &                        t(i  ,j,k,nstp,itrc))*adfac-              &
     &                       (curv(i-1,j)*                              &
     &                        (0.5_r8+SIGN(0.5_r8, Huon(i,j,k)))+       &
     &                        curv(i  ,j)*                              &
     &                        (0.5_r8+SIGN(0.5_r8,-Huon(i,j,k))))*      &
     &                       adfac2
              ad_t(i-1,j,k,nstp,itrc)=ad_t(i-1,j,k,nstp,itrc)+adfac1
              ad_t(i  ,j,k,nstp,itrc)=ad_t(i  ,j,k,nstp,itrc)+adfac1
              ad_curv(i-1,j)=ad_curv(i-1,j)-                            &
     &                       MAX(Huon(i,j,k),0.0_r8)*adfac2
              ad_curv(i  ,j)=ad_curv(i  ,j)-                            &
     &                       MIN(Huon(i,j,k),0.0_r8)*adfac2
              ad_FX(i,j)=0.0_r8
#   else
!>            tl_FX(i,j)=0.5_r8*                                        &
!>   &                   (tl_Huon(i,j,k)*                               &
!>   &                    (t(i-1,j,k,nstp,itrc)+                        &
!>   &                     t(i  ,j,k,nstp,itrc)-                        &
!>   &                     cff2*(grad(i  ,j)-                           &
!>   &                           grad(i-1,j)))+                         &
!>   &                    Huon(i,j,k)*                                  &
!>   &                    (tl_t(i-1,j,k,nstp,itrc)+                     &
!>   &                     tl_t(i  ,j,k,nstp,itrc)-                     &
!>   &                     cff2*(tl_grad(i  ,j)-                        &
!>   &                           tl_grad(i-1,j))))
!>
              adfac=0.5_r8*ad_FX(i,j)
              adfac1=adfac*Huon(i,j,k)
              adfac2=adfac1*cff2
              ad_Huon(i,j,k)=ad_Huon(i,j,k)+                            &
     &                       adfac*(t(i-1,j,k,nstp,itrc)+               &
     &                              t(i  ,j,k,nstp,itrc)-               &
     &                              cff2*(grad(i  ,j)-                  &
     &                                    grad(i-1,j)))
              ad_t(i-1,j,k,nstp,itrc)=ad_t(i-1,j,k,nstp,itrc)+adfac1
              ad_t(i  ,j,k,nstp,itrc)=ad_t(i  ,j,k,nstp,itrc)+adfac1
              ad_grad(i-1,j)=ad_grad(i-1,j)+adfac2
              ad_grad(i  ,j)=ad_grad(i  ,j)-adfac2
              ad_FX(i,j)=0.0_r8
#   endif
            END DO
          END DO

          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
#   if defined TS_U3HADVECTION_TL
!>            tl_curv(i,j)=tl_FX(i+1,j)-tl_FX(i,j)
!>
              ad_FX(i  ,j)=ad_FX(i  ,j)-ad_curv(i,j)
              ad_FX(i+1,j)=ad_FX(i+1,j)+ad_curv(i,j)
              ad_curv(i,j)=0.0_r8
#   elif defined TS_A4HADVECTION_TL
              cff=2.0_r8*FX(i+1,j)*FX(i,j)
              IF (cff.gt.eps) THEN
!>              tl_grad(i,j)=((FX(i+1,j)+FX(i,j))*tl_cff-               &
!>   &                        cff*(tl_FX(i+1,j)+tl_FX(i,j)))/           &
!>   &                       ((FX(i+1,j)+FX(i,j))*(FX(i+1,j)+FX(i,j)))
!>
                adfac=ad_grad(i,j)/                                     &
     &                ((FX(i+1,j)+FX(i,j))*(FX(i+1,j)+FX(i,j)))
                adfac1=adfac*cff
                ad_FX(i  ,j)=ad_FX(i  ,j)-adfac1
                ad_FX(i+1,j)=ad_FX(i+1,j)-adfac1
                ad_cff=ad_cff+(FX(i+1,j)+FX(i,j))*adfac
                ad_grad(i,j)=0.0_r8
              ELSE
!>              tl_grad(i,j)=0.0_r8
!>
                ad_grad(i,j)=0.0_r8
              END IF
#   else
!>            tl_grad(i,j)=0.5_r8*(tl_FX(i+1,j)+tl_FX(i,j))
!>
              adfac=0.5_r8*ad_grad(i,j)
              ad_FX(i  ,j)=ad_FX(i  ,j)+adfac
              ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
              ad_grad(i,j)=0.0_r8
#   endif
!>            tl_cff=2.0_r8*(tl_FX(i+1,j)*FX(i,j)+                      &
!>   &                       FX(i+1,j)*tl_FX(i,j))
!>
              adfac=2.0_r8*ad_cff
              ad_FX(i  ,j)=ad_FX(i  ,j)+FX(i+1,j)*adfac
              ad_FX(i+1,j)=ad_FX(i+1,j)+FX(i  ,j)*adfac
              ad_cff=0.0_r8
            END DO
          END DO
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO j=Jstr,Jend
!>              tl_FX(Iend+2,j)=tl_FX(Iend+1,j)
!>
                ad_FX(Iend+1,j)=ad_FX(Iend+1,j)+ad_FX(Iend+2,j)
                ad_FX(Iend+2,j)=0.0_r8
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO j=Jstr,Jend
!>              tl_FX(Istr-1,j)=tl_FX(Istr,j)
!>
                ad_FX(Istr,j)=ad_FX(Istr,j)+ad_FX(Istr-1,j)
                ad_FX(Istr-1,j)=0.0_r8
              END DO
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
#   ifdef MASKING
!>            tl_FX(i,j)=tl_FX(i,j)*umask(i,j)
!>
              ad_FX(i,j)=ad_FX(i,j)*umask(i,j)
#   endif
!>            tl_FX(i,j)=tl_t(i  ,j,k,nstp,itrc)-                       &
!>   &                   tl_t(i-1,j,k,nstp,itrc)
!>
              ad_t(i-1,j,k,nstp,itrc)=ad_t(i-1,j,k,nstp,itrc)-ad_FX(i,j)
              ad_t(i  ,j,k,nstp,itrc)=ad_t(i  ,j,k,nstp,itrc)+ad_FX(i,j)
              ad_FX(i,j)=0.0_r8
            END DO
          END DO
#  endif
        END DO K_LOOP
      END DO T_LOOP1

#  ifdef SOLAR_SOURCE
!
!  Compute fraction of the solar shortwave radiation, "swdk"
!  (at vertical W-points) penetrating water column. First, compute
!  BASIC STATE arrays.
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            FX(i,j)=z_w(i,j,N(ng))-z_w(i,j,k)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_swdk(i,j,k)=tl_FE(i,j)
!>
            ad_FE(i,j)=ad_FE(i,j)+ad_swdk(i,j,k)
            ad_swdk(i,j,k)=0.0_r8
          END DO
        END DO
!>      CALL tl_lmd_swfrac_tile (ng, tile,                              &
!>   &                           LBi, UBi, LBj, UBj,                    &
!>   &                           IminS, ImaxS, JminS, JmaxS,            &
!>   &                           -1.0_r8, FX, tl_FX, tl_FE)
!>
        CALL ad_lmd_swfrac_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           -1.0_r8, FX, ad_FX, ad_FE)
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_FX(i,j)=tl_z_w(i,j,N(ng))-tl_z_w(i,j,k)
!>
            ad_z_w(i,j,k    )=ad_z_w(i,j,k    )-ad_FX(i,j)
            ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))+ad_FX(i,j)
            ad_FX(i,j)=0.0_r8
          END DO
        END DO
      END DO
#  endif
# endif /* !TS_FIXED */

      RETURN
      END SUBROUTINE ad_pre_step3d_tile
#endif
      END MODULE ad_pre_step3d_mod
