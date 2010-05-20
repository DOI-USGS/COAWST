#include "cppdefs.h"
      MODULE wec_mod
#ifdef NEARSHORE_WEC
!
!svn $Id: wec.F 1428 2008-03-12 13:07:21Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!================================================== John C. Warner ====!
!                                                                      !
!  This routine computes the terms corresponding to vortex forces in   !
!  momentum equations.                                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Uchiyama, Y., McWilliams, J.C., and Shchepetkin, A.F. (2010).       !
!  Wave current interacation in an oceanic circulation model with a    !
!  vortex-force formalism: Applications to surf zone, Ocean Modeling,  !
!  xx, xx-yy.
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: wec
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec (ng, tile)
!***********************************************************************
!
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
      USE mod_coupling
# if defined DIAGNOSTICS_UV
      USE mod_diags
# endif
!
      integer, intent(in) :: ng, tile
# include "tile.h"
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 21)
# endif
      CALL wec_tile (ng, tile,    LBi, UBi, LBj, UBj, N(ng),            &
     &                            IminS, ImaxS, JminS, JmaxS,           &
# ifdef SOLVE3D
     &                            nrhs(ng),                             &
# endif
# ifdef MASKING
     &                            GRID(ng) % pmask,                     &
     &                            GRID(ng) % rmask,                     &
     &                            GRID(ng) % umask,                     &
     &                            GRID(ng) % vmask,                     &
# endif
# ifdef WET_DRY
     &                            GRID(ng) % umask_wet,                 &
     &                            GRID(ng) % vmask_wet,                 &
# endif
     &                            GRID(ng) % om_u,                      &
     &                            GRID(ng) % om_v,                      &
     &                            GRID(ng) % on_u,                      &
     &                            GRID(ng) % on_v,                      &
     &                            GRID(ng) % pm,                        &
     &                            GRID(ng) % pn,                        &
     &                            GRID(ng) % angler,                    &
# if defined CURVGRID
     &                            GRID(ng) % dndx,                      &
     &                            GRID(ng) % dmde,                      &
# endif
     &                            GRID(ng) % h,                         &
     &                            OCEAN(ng) % zeta,                     &
     &                            OCEAN(ng) % zetat,                    &
# ifdef SOLVE3D
     &                            GRID(ng) % Hz,                        &
     &                            GRID(ng) % z_r,                       &
     &                            GRID(ng) % z_w,                       &
# endif
     &                            FORCES(ng) % Hwave,                   &
     &                            FORCES(ng) % Dwave,                   &
     &                            FORCES(ng) % Lwave,                   &
     &                            FORCES(ng) % Pair,                    &
     &                            FORCES(ng) % Pwave_top,               &
# ifdef SOLVE3D
#  ifdef DIAGNOSTICS_UV
     &                            DIAGS(ng) % DiaRU,                    &
     &                            DIAGS(ng) % DiaRV,                    &
#  endif
     &                            MIXING(ng) % rustr3d,                 &
     &                            MIXING(ng) % rvstr3d,                 &
     &                            OCEAN(ng) % rulag3d,                  &
     &                            OCEAN(ng) % rvlag3d,                  &
     &                            OCEAN(ng) % u,                        &
     &                            OCEAN(ng) % v,                        &
     &                            OCEAN(ng) % u_stokes,                 &
     &                            OCEAN(ng) % v_stokes,                 &
     &                            OCEAN(ng) % w_stokes,                 &
# endif
     &                            MIXING(ng) % rustr2d,                 &
     &                            MIXING(ng) % rvstr2d,                 &
     &                            OCEAN(ng) % rulag2d,                  &
     &                            OCEAN(ng) % rvlag2d,                  &
     &                            OCEAN(ng) % ubar,                     &
     &                            OCEAN(ng) % vbar,                     &
     &                            OCEAN(ng) % ubar_stokes,              &
     &                            OCEAN(ng) % vbar_stokes)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 21)
# endif
      RETURN
      END SUBROUTINE wec
!
!***********************************************************************
      SUBROUTINE wec_tile (ng, tile,    LBi, UBi, LBj, UBj, UBk,        &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
# ifdef SOLVE3D
     &                                  nrhs,                           &
# endif
# ifdef MASKING
     &                                  pmask, rmask, umask, vmask,     &
# endif
# ifdef WET_DRY
     &                                  umask_wet, vmask_wet,           &
# endif
     &                                  om_u, om_v, on_u, on_v,         &
     &                                  pm, pn,                         &
     &                                  angler,                         &
# if defined CURVGRID
     &                                  dndx, dmde,                     &
# endif
     &                                  h, zeta, zetat,                 &
# ifdef SOLVE3D
     &                                  Hz, z_r, z_w,                   &
# endif
     &                                  Hwave, Dwave, Lwave, Pair,      &
     &                                  Pwave_top,                      &
# ifdef SOLVE3D
#  ifdef DIAGNOSTICS_UV
     &                                  DiaRU, DiaRV,                   &
#  endif
     &                                  rustr3d,  rvstr3d,              &
     &                                  rulag3d, rvlag3d,               &
     &                                  u, v,                           &
     &                                  u_stokes, v_stokes,w_stokes,    &
# endif
     &                                  rustr2d, rvstr2d,               &
     &                                  rulag2d, rvlag2d, ubar, vbar,   &
     &                                  ubar_stokes, vbar_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
      USE bc_2d_mod
      USE bc_3d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef SOLVE3D
      integer, intent(in) :: nrhs
# endif
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
#  if defined CURVGRID
      real(r8), intent(in) :: dndx(LBi:,LBj:)
      real(r8), intent(in) :: dmde(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: zetat(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  endif
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_top(LBi:,LBj:)
#  ifdef SOLVE3D
#   ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: rustr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvstr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rulag3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvlag3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: w_stokes(LBi:,LBj:,0:)
#  endif
      real(r8), intent(inout) :: rustr2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvstr2d(LBi:,LBj:)
      real(r8), intent(inout) :: rulag2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvlag2d(LBi:,LBj:)
      real(r8), intent(in)    :: ubar(LBi:,LBj:,:)
      real(r8), intent(in)    :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ubar_stokes(LBi:,LBj:)
      real(r8), intent(inout) :: vbar_stokes(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
#  if defined CURVGRID
      real(r8), intent(in) :: dndx(LBi:UBI,LBj:UBj)
      real(r8), intent(in) :: dmde(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zetat(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
#  endif
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pwave_top(LBi:UBi,LBj:UBj)

#  ifdef SOLVE3D
#   ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#   endif
      real(r8), intent(inout) :: rustr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rvstr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rulag3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rvlag3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: w_stokes(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(inout) :: rustr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvstr2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rulag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvlag2d(LBi:UBi,LBj:UBj)
      real(r8), intent(in)    :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in)    :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ubar_stokes(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vbar_stokes(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
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
      integer :: i, j, k, numits, it
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6
      real(r8) :: fac1, fac2, RB1, RB2, EBAR
      real(r8), parameter :: alphaR=0.0_r8
      real(r8), parameter :: sinb=0.1_r8
      real(r8), parameter :: gammaw=0.31_r8
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavec
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owaven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveAA
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveRE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: WD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigmat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zetaw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rollA
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rollAA
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gamr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EER
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBX_bar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBY_bar
                 
# ifdef SOLVE3D
      real(r8) :: vi1, vi2
      real(r8), dimension(IminS:ImaxS) :: wrk
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCS
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FSS
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FSC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF1
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: oroller
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBXL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBYL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huons
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvoms
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: z_psi
	  
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dU
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dku
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dkv
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ddku
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ddkv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vf2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: bh
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: qsp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: ku
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: kv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shearu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shearv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shear1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shear2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: vfvert
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: ustux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: ustuy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: vstvx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: vstvy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: BBX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: BBY      
#  ifdef CURVGRID
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
#  endif
# endif
# include "set_bounds.h"
      fac1=1.0_r8/dt(ng)
      fac2=1.0_r8/(g*rho0)
# ifdef SOLVE3D
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!
!  Compute total depth
!
          Dstp(i,j)=zeta(i,j,1)+h(i,j)
!
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!
          waven(i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          owaven(i,j)=1.0_r8/waven(i,j)
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff)
          waveny(i,j)=waven(i,j)*SIN(cff)
          cff=0.25_r8*SQRT(2.0_r8)*Hwave(i,j)
          waveAA(i,j)=cff*cff
          sigma(i,j)=SQRT(g*waven(i,j)*TANH(waven(i,j)*Dstp(i,j)))
          osigma(i,j)=1.0_r8/sigma(i,j)
          waveE(i,j)=0.0625_r8*g*Hwave(i,j)*Hwave(i,j)
!
!  Compute wave celerity and nonlinear water depth
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j)+eps,kDmax)
          wavec(i,j)=SQRT(g*owaven(i,j)*TANH(kD(i,j)))          
!
!  Compute quasi static sea surface elevation, Eq(7).
!
          zetaw(i,j)=-fac2*Pair(i,j)-(0.5_r8*waveAA(i,j)*waven(i,j)/    &
     &               SINH(2.0_r8*kD(i,j)))
!
!  Compute metrics for vertical roller distribution.
!
          oroller(i,j)=0.0_r8
          gamr(i,j)=MIN(0.707_r8*Dstp(i,j)/                             &
     &              (1.25_r8*Hwave(i,j)+eps),1.0_r8)
          DO k=1,N(ng)
            cff2=(1.0_r8+SCALARS(ng)%Cs_r(k))*gamr(i,j)
            oroller(i,j)=oroller(i,j)+Hz(i,j,k)*                        &
     &                COSH(2.0_r8*pi*cff2)
          END DO
          oroller(i,j)=1.0_r8/(oroller(i,j)+eps)
!
!  Calculate wave dissipation using empirical parameters of
!  Church Thornton 1993.
! 
          RB1= 1.0_r8+tanh(8.0_r8*((0.707_r8*Hwave(i,j)/                &
     &          (gammaw*Dstp(i,j)))-1.0_r8))
          RB2= 1.0_r8-(1.0_r8+(0.707_r8*Hwave(i,j)/                     &
     &          (gammaw*Dstp(i,j)))**2.0_r8)**(-2.5_r8) 
          sigmat(i,j)=MIN(1.0_r8/Pwave_top(i,j),1.0_r8)
          WD(i,j)=(3.0_r8*sqrt(pi)/16.0_r8)*g*sigmat(i,j)*              &
     &                (0.64_r8**3.0_r8/Dstp(i,j))*                      &
     &                 (0.707_r8*Hwave(i,j))**3.0_r8*RB1*RB2
!
!  Initialize depth independent arrays for summation
!
        END DO
        IF ((j.ge.Jstr) .and. (j.le.Jend)) THEN
          DO i=Istr,Iend
            BBX_bar(i,j)=0.0_r8
            BBY_bar(i,j)=0.0_r8
            w_stokes(i,j,0)=0.0_r8
          END DO
        END IF
      END DO
!
!  Solve roller evolution equation
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          rollAA(i,j)=rollA(i,j)
        END DO
      END DO
      numits=30
      DO it=1,numits
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
             cff3=(ubar(i,j,nrhs)+wavenx(i,j)*owaven(i,j)*wavec(i,j))*  &
     &          on_u(i,j)

            cff1=MAX(cff3,0.0_r8)
            cff2=MIN(cff3,0.0_r8)
            FX(i,j)=cff1*rollA(i-1,j)+cff2*rollA(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff3=(vbar(i,j,nrhs)+waveny(i,j)*owaven(i,j)*wavec(i,j))*   &
     &           om_v(i,j)*1.0_r8
            cff1=MAX(cff3,0.0_r8)
            cff2=MIN(cff3,0.0_r8)
            FE(i,j)=cff1*rollA(i,j-1)+cff2*rollA(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)/REAL(numits,r8)
            cff1=cff*(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))
           rollAA(i,j)=rollAA(i,j)-cff1
          END DO
        END DO
!
! Roller source / sink term.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)/REAL(numits,r8)
            EBAR=-alphaR*WD(i,j)/sigma(i,j)+g*sinb*rollA(i,j)/wavec(i,j)
            rollAA(i,j)=rollAA(i,j)-cff*EBAR
            rollA(i,j)=rollAA(i,j)
            waveRE(i,j)=rollA(i,j)*sigma(i,j)
!
! Compute roller dissipation rate here (Eqn. 41 U10)
!
            EER(i,j)=g*sinb*waveRE(i,j)/wavec(i,j) 
          END DO
        END DO
        DO j=Jstr-1,Jend+1
          rollA(Istr-1,j)=rollA(Istr,j)
          waveRE(Istr-1,j)=waveRE(Istr,j)
          EER(Istr-1,j)=EER(Istr,j)
          rollA(Iend+1,j)=rollA(Iend,j)
          waveRE(Iend+1,j)=waveRE(Iend,j)
          EER(Iend+1,j)=EER(Iend,j)
        END DO
        DO i=Istr-1,Iend+1       
          rollA(i,Jstr-1)=rollA(i,Jstr)
          waveRE(i,Jstr-1)=waveRE(i,Jstr)
          EER(i,Jstr-1)=EER(i,Jstr)
          rollA(i,Jend+1)=rollA(i,Jend)
          waveRE(i,Jend+1)=waveRE(i,Jend)
          EER(i,Jend+1)=EER(i,Jend)
        END DO
      END DO
!
!  Compute Bernoulli's Head (BH)
!  BH has been represented as the symbol cursive K in the U10 (see Eqn. 5)
!
!  Use splines method to calculate shear of velocity
!
#  ifdef SPLINES_Notyet
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          CF(i,0)=0.0_r8
          dU(i,0)=0.0_r8
          dV(i,0)=0.0_r8
        END DO

        DO k=1,N(ng)-1
          DO i=Istr-1,Iend+1
             ku(i,j,k)= 0.5_r8*(u(i,j,k,nrhs)+u(i-1,j,k,nrhs))*         &
     &                  wavenx(i,j)
             kv(i,j,k)= 0.5_r8*(v(i,j,k,nrhs)+v(i,j-1,k,nrhs))*         &
     &                  waveny(i,j)
             cff= 1.0_r8/(2.0_r8*Hz(i,j,k+1)+                           &
     &                    Hz(i,j,k)*(2.0_r8-CF(i,k-1)))
             CF(i,k)= cff*Hz(i,j,k+1)
             dU(i,k)= cff*(3.0_r8*(u(i  ,j,k+1,nrhs)-u(i,  j,k,nrhs)+   &
     &                             u(i+1,j,k+1,nrhs)-u(i+1,j,k,nrhs))-  &
     &                Hz(i,j,k)*dU(i,k-1))
             dV(i,k)= cff*(3.0_r8*(v(i,j  ,k+1,nrhs)-v(i,j  ,k,nrhs)+   &
     &                             v(i,j+1,k+1,nrhs)-v(i,j+1,k,nrhs))-  &
     &                Hz(i,j,k)*dV(i,k-1))              
          END DO
        END DO            

        DO i=Istr-1,Iend+1
          dU(i,N(ng))=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO

        DO k=N(ng)-1,1,-1
          DO i=Istr-1,Iend+1
            dU(i,k)=dU(i,k)-CF(i,k)*dU(i,k+1)
            dV(i,k)=dV(i,k)-CF(i,k)*dV(i,k+1)
          END DO
        END DO
		
        DO k=1,N(ng)-1
          DO i=Istr-1,Iend+1
!
! Calculate du/dz (i.e. (1/Hz)*(du/ds))
!
            shearu(i,j,k)=dU(i,k)
!
! Calculate dv/dz (i.e. (1/Hz)*(dv/ds))
            shearv(i,j,k)=dV(i,k)

          END DO
        END DO
      END DO

!	  
! Apply boundary conditions here:	  
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
           ku(i,j,N(ng))=ku(i,j,N(ng)-1)
           kv(i,j,N(ng))=kv(i,j,N(ng)-1)
        END DO
      END DO
	  
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
           ku(Iend+1,j,k)=ku(Iend,j,k)
        END DO

        DO i=Istr-1,Iend+1
           kv(i,Jend+1,k)=kv(i,Jend,k)
        END DO
      END DO

! 
! Calculate Shear for KU and KV
!

      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          CF(i,0)=0.0_r8
          dku(i,j,0)=0.0_r8
          dkv(i,j,0)=0.0_r8
        END DO

        DO k=1,N(ng)-1
          DO i=Istr-1,Iend+1
             cff= 1.0_r8/(2.0_r8*Hz(i,j,k+1)+                           &
     &                    Hz(i,j,k)*(2.0_r8-CF(i,k-1)))
             CF(i,k)= cff*Hz(i,j,k+1)

             dku(i,j,k)= cff*(3.0_r8*(ku(i  ,j,k+1)-ku(i,  j,k)+        &
     &                              ku(i+1,j,k+1)-ku(i+1,j,k))-         &
     &                 Hz(i,j,k)*dku(i,j,k-1))
             dkv(i,j,k)= cff*(3.0_r8*(kv(i,j  ,k+1)-kv(i,j  ,k)+        &
     &                              kv(i,j+1,k+1)-kv(i,j+1,k))-         &
     &                 Hz(i,j,k)*dkv(i,j,k-1))              
          END DO
        END DO            

        DO i=Istr-1,Iend+1
          dku(i,j,N(ng))=0.0_r8
          dkv(i,j,N(ng))=0.0_r8
        END DO

        DO k=N(ng)-1,1,-1
          DO i=Istr-1,Iend+1
            dku(i,j,k)=dku(i,j,k)-CF(i,k)*dku(i,j,k+1)
            dkv(i,j,k)=dkv(i,j,k)-CF(i,k)*dkv(i,j,k+1)
          END DO
        END DO
		
        DO k=1,N(ng)-1
          DO i=Istr-1,Iend+1
!
! Calculate d(k_x.u)/dz + d(k_y.v)/dz (i.e. (1/Hz)*(d(k_x.u)/ds+d(k_y.v/ds))
!
            shear1(i,j,k)=dku(i,j,k)+dkv(i,j,k)

          END DO
        END DO
      END DO


!	  
! Apply boundary conditions here:	  
!
      DO i=Istr-1,Iend+1
        DO j=Jstr-1,Jend+1    
         dku(i,j,N(ng))=dku(i,j,N(ng)-1)
         dkv(i,j,N(ng))=dkv(i,j,N(ng)-1)
         shear1(i,j,N(ng))=shear1(i,j,N(ng)-1)
        END DO
      END DO
      
      DO k=1,N(ng)	  
        DO j=Jstr-1,Jend+1
         dku(Iend+1,j,k)=dku(Iend,j,k)
         dkv(Iend+1,j,k)=dkv(Iend,j,k)
        END DO

        DO i=Istr-1,Iend+1
         dku(i,Jend+1,k)=dku(i,Jend,k)
         dkv(i,Jend+1,k)=dkv(i,Jend,k)
        END DO
      END DO
      

!
!Calculate second order shear terms
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          CF(i,0)=0.0_r8
          ddku(i,0)=0.0_r8
          ddkv(i,0)=0.0_r8
        END DO

        DO k=1,N(ng)-1
          DO i=Istr-1,Iend+1
             cff= 1.0_r8/(2.0_r8*Hz(i,j,k+1)+                           &
     &                    Hz(i,j,k)*(2.0_r8-CF(i,k-1)))
             CF(i,k)= cff*Hz(i,j,k+1)

             ddku(i,k)= cff*(3.0_r8*(dku(i  ,j,k+1)-dku(i,  j,k)+       &
     &                               dku(i+1,j,k+1)-dku(i+1,j,k))-      &
     &                  Hz(i,j,k)*ddku(i,k-1))
             ddkv(i,k)= cff*(3.0_r8*(dkv(i,j  ,k+1)-dkv(i,j  ,k)+       &
     &                               dkv(i,j+1,k+1)-dkv(i,j+1,k))-      &
     &                  Hz(i,j,k)*ddkv(i,k-1))              
          END DO
        END DO            

        DO i=Istr-1,Iend+1
          ddku(i,N(ng))=0.0_r8
          ddkv(i,N(ng))=0.0_r8
        END DO

        DO k=N(ng)-1,1,-1
          DO i=Istr-1,Iend+1
            ddku(i,k)=ddku(i,k)-CF(i,k)*ddku(i,k+1)
            ddkv(i,k)=ddkv(i,k)-CF(i,k)*ddkv(i,k+1)
          END DO
        END DO
		
        DO k=1,N(ng)-1
          DO i=Istr-1,Iend+1
!
! Calculate d^2(k_x.u)/dz^2 + d^2(k_y.v)/dz^2 
!
            shear2(i,j,k)=ddku(i,k)+ddkv(i,k)

          END DO
        END DO
      END DO

#  else      
!
!  Compute Bernoulli head, Eq 5.
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
!
! Compute k*u and k*v. This will be at rho pts.
!
            ku(i,j,k)=0.5_r8*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))            &
                      *wavenx(i,j)
            kv(i,j,k)= 0.5_r8*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))           &
                      *waveny(i,j)
!
!  Calculates du/dz (i.e. (1/Hz)*(du/ds)). This will be at rho and w pts.
! 
            shearu(i,j,k)=cff*(u(i  ,j,k+1,nrhs)-u(i  ,j,k,nrhs)+       &
     &                         u(i+1 ,j,k+1,nrhs)-u(i+1 ,j,k,nrhs))
!
!  Calculates dv/dz (i.e. (1/Hz)*(dv/ds)).
!
            shearv(i,j,k)=cff*(v(i,j  ,k+1,nrhs)-v(i,j  ,k,nrhs)+       &
     &                         v(i,j+1,k+1,nrhs)-v(i,j+1,k,nrhs))
          END DO
        END DO
      END DO
!
!  Applying boundary conditions.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ku(i,j,N(ng))=ku(i,j,N(ng)-1)
          kv(i,j,N(ng))=kv(i,j,N(ng)-1)
          shearu(i,j,N(ng))=shearu(i,j,N(ng)-1)
          shearu(i,j,0)=shearu(i,j,1)
          shearv(i,j,N(ng))=shearv(i,j,N(ng)-1)
          shearv(i,j,0)=shearv(i,j,1)
        END DO
      END DO
!
!  Compute gradients for shear.
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
!
!  Calculates d(k_x*u)/dz+d(k_y*v)/dz (i.e. (1/Hz)*(d(k_x*u)/ds+d(k_y*v)/ds)
!
!  Shear1 is at w points and rho horiz points.

            shear1(i,j,k)=cff*((ku(i,j,k+1)-ku(i,j,k  ))+               &
     &                         (kv(i,j,k+1)-kv(i,j,k  )))
          END DO
        END DO
      END DO
!
! Applying boundary conditions here.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          shear1(i,j,N(ng))=shear1(i,j,N(ng)-1)
          shear1(i,j,0)=shear1(i,j,1)
        END DO
      END DO
!
! Compute second derivative.
!
      DO k=0,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=0.5_r8/(z_w(i,j,k+1)-z_w(i,j,k))
!
! Calculates d(k_x*u)^2/dz^2+d(k_y*v)^2/dz^2 
! in the form  (1/Hz)*(d/ds)*(1/Hz*d(k_x*u)/ds+1/Hz*d(k_y*v)/ds)
!
!           at k=0 shear2=0 because shear1(1)=shear1(0).
            shear2(i,j,k)=(cff*(shear1(i,j,k+1)-shear1(i,j,k)))
          END DO
        END DO
      END DO
!
! Applying boundary conditions here.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          shear2(i,j,N(ng))=shear2(i,j,N(ng)-1)
        END DO
      END DO
#  endif
!
!  Perform integration on the shear2 term calcuated above.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
! there is no Hz(,,0)
!           cff=Hz(i,j,0)
!           fac2=(1.0_r8+SCALARS(ng)%Cs_r(N(ng)))-                       &
!                (1.0_r8+SCALARS(ng)%Cs_r(0))
!           vi1=cff*shear2(i,j,0)*SINH(2*kD(i,j)*fac2)                   
           vi1=0.0_r8
!
          DO k=1,N(ng)
             fac2=(1.0_r8+SCALARS(ng)%Cs_r(N(ng)))-                     &
     &            (1.0_r8+SCALARS(ng)%Cs_r(k))
             cff=Hz(i,j,k)
             vi1=vi1+cff*shear2(i,j,k)*SINH(2*kD(i,j)*fac2)             
          END DO
!
!  Calculate Bernoulli's Head at the surface
!
          bh(i,j)=0.25_r8*vi1*sigma(i,j)*waveAA(i,j)/                   & 
     &                    (waven(i,j)*(SINH(kD(i,j)))**2)

        END DO
      END DO

!
! Compute quasi-static components of pressure terms (i.e. Cursive P)
! as in Eqn. 9 in U10
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!          fac2=1.0_r8+SCALARS(ng)%Cs_r(0)
!          cff=Hz(i,j,0)
!          vi2=cff*shear2(i,j,0)*COSH(2*kD(i,j)*fac2)
          vi2=0.0_r8
          DO k=1,N(ng)
             fac2=1.0_r8+SCALARS(ng)%Cs_r(k)
             cff=Hz(i,j,k)
             vi2=vi2+cff*shear2(i,j,k)*COSH(2*kD(i,j))*fac2	  
          END DO
          cff1=0.5_r8*waveAA(i,j)*osigma(i,j)
          cff2=TANH(kD(i,j))/SINH(2*kD(i,j))
          qsp(i,j)=cff1*(cff2*(shear1(i,j,N(ng))+                       &
     &                         COSH(2*kD(i,j))*shear1(i,j,0)+           &
     &                         vi2)-                                    &
     &                   2*waven(i,j)*TANH(kD(i,j))*                    &
     &                   (ku(i,j,N(ng))+kv(i,j,N(ng))))
        END DO
      END DO
!
! Compute Vertical Vortex Force terms (denoted as K in Eqn. 5, U10):
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            vfvert(i,j,k)=0.25_r8*(u_stokes(i,j,k)+u_stokes(i-1,j,k))*   &
     &	                  (shearu(i,j,k)+shearu(i,j,k+1))+               &
     &                    0.25_r8*(v_stokes(i,j,k)+v_stokes(i,j-1,k))*   &
     &                    (shearv(i,j,k)+shearv(i,j,k+1))
          END DO
        END DO
      END DO
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          vfvert(Istr-1,j,k)=vfvert(Istr,j,k)
        END DO
        
        DO i=Istr,Iend
          vfvert(i,Jstr-1,k)=vfvert(i,Jstr,k)
        END DO
      END DO		
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
          rustr3d(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*                & 
     &                   (vfvert(i  ,j,k)*Hz(i  ,j,k)-                  &
     &                    vfvert(i-1,j,k)*Hz(i-1,j,k))*on_u(i,j)
          rvstr3d(i,j,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*                & 
     &                   (vfvert(i,j  ,k)*Hz(i,j  ,k)-                  &
     &                    vfvert(i,j-1,k)*Hz(i,j-1,k))*on_v(i,j)       
          END DO
        END DO
      END DO 
!
!---------------------------------------------------------------------------
! Evaluate geopotential function (Eqn. 33, Uchiyama et al., 2010)
!---------------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zetat(i,j)=Dstp(i,j)*                                         &
     &                      (-g*zetaw(i,j)-(g*qsp(i,j)-bh(i,j)))
        END DO
      END DO
!
!------------------------------------------------------------------------
! Compute non conservative wave acceleration terms here
!------------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
        fac2=1.0_r8+SCALARS(ng)%Cs_r(k)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            cff2=fac2*gamr(i,j)
            cff3=COSH(2.0_r8*pi*cff2)
            cff4=(((1-alphaR)*WD(i,j)+EER(i,j))/sigma(i,j))
            BBXL(i,j)=cff3*cff4*wavenx(i,j)*oroller(i,j)
            BBYL(i,j)=cff3*cff4*waveny(i,j)*oroller(i,j) 
          END DO
!
! Compute depth integrated terms and copy stress terms from local arrays
! to state variables.
!
          IF ((j.ge.Jstr) .and. (j.le.Jend)) THEN
            DO i=Istr,Iend
               BBX(i,j,k)=BBXL(i,j)
               BBY(i,j,k)=BBYL(i,j)
               BBX_bar(i,j)=BBX_bar(i,j)+Hz(i,j,k)*BBX(i,j,k)
               BBY_bar(i,j)=BBY_bar(i,j)+Hz(i,j,k)*BBY(i,j,k)
            END DO
          END IF
        END DO 
!
! Compute contribution to U-momentum
! 
        DO j=Jstr,Jend
          DO i=IstrU,Iend
             cff=(BBXL(i  ,j)*Hz(i  ,j,k)-                              &
     &            BBXL(i-1,j)*Hz(i-1,j,k))*on_u(i,j)
             rustr3d(i,j,k)=rustr3d(i,j,k)+cff
# ifdef DIAGNOSTICS_UV
             DiaRU(i,j,k,nrhs,M3hrad)=-cff
# endif
          END DO
        END DO
!
! Compute contribution to V-momentum
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
             cff=(BBYL(i,j  )*Hz(i,j  ,k)-                              &
     &            BBYL(i,j-1)*Hz(i,j-1,k))*om_v(i,j)
             rvstr3d(i,j,k)=rvstr3d(i,j,k)+cff
# ifdef DIAGNOSTICS_UV
             DiaRV(i,j,k,nrhs,M3hrad)=-cff
# endif
          END DO
        END DO
!	  
! Compute u-stokes velocities. See Eqn. 2 in U10
!
        fac2=1.0_r8+SCALARS(ng)%Cs_r(k)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
             cff2=(waveE(i-1,j)+waveE(i,j))
!  NEED to ADD BBR ?? to wave E for here and vstokes
             cff3=(kD(i-1,j)+kD(i,j))
             u_stokes(i,j,k)=cff2*                                      &
     &                      (wavenx(i-1,j)+wavenx(i,j))/                &
     &                      (wavec (i-1,j)+wavec (i,j))*                & 
     &                      COSH(cff3*fac2)/SINH(cff3)

#  ifdef MASKING
             u_stokes(i,j,k)=u_stokes(i,j,k)*umask(i,j)
#  endif
#  ifdef WET_DRY
             u_stokes(i,j,k)=u_stokes(i,j,k)*umask_wet(i,j)
#  endif
          END DO
        END DO
!
!  Compute v-stokes velocity.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
             cff2=(waveE(i,j-1)+waveE(i,j))
             cff3=(kD(i,j-1)+kD(i,j))
             v_stokes(i,j,k)=cff2*                                      &
     &                      (waveny(i-1,j)+waveny(i,j))/                &
     &                      (wavec (i,j-1)+wavec (i,j))*                &
     &                      COSH(cff3*fac2)/SINH(cff3)
#  ifdef MASKING
             v_stokes(i,j,k)=v_stokes(i,j,k)*vmask(i,j)
#  endif
#  ifdef WET_DRY
             v_stokes(i,j,k)=v_stokes(i,j,k)*vmask_wet(i,j)
#  endif
          END DO
        END DO
!	  
!  Calculate vertical stokes velocity terms.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
             Huons(i,j)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*                 &
     &                  u_stokes(i,j,k)*on_u(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
             Hvoms(i,j)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*                 &
     &                  v_stokes(i,j,k)*om_v(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
             w_stokes(i,j,k)=w_stokes(i,j,k-1)-                         &
     &                      (Huons(i+1,j)-Huons(i,j)+                   &
     &                       Hvoms(i,j+1)-Hvoms(i,j))*                  &
     &                       pm(i,j)*pn(i,j)
          END DO
        END DO
      END DO K_LOOP
!
! Applying boundary condition for w_stokes
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wrk(i)=w_stokes(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
        END DO
!
!  In order to insure zero vertical velocity at the free-surface,
!  subtract the vertical velocities of the moving S-coordinates
!  isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
!  The proportionally coefficients are a linear function of the
!  S-coordinate with zero value at the bottom (k=0) and unity at
!  the free-surface (k=N).
!
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            w_stokes(i,j,k)=w_stokes(i,j,k)-                            &
     &                      wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          END DO
        END DO
        DO i=Istr,Iend
          w_stokes(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
! Contribution of a term corresponding to product of Stokes and Eulerian Velocity
! Eqn. 26 and 27, Uchiyama et al. (2010)
!
      DO k=1,N(ng)      
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            ustux(i,j,k)=0.5_r8*                                        &
     &                   (u_stokes(i,j,k)+u_stokes(i-1,j,k))*           &
     &		         (u(i+1,j,k,nrhs)-u(i,j,k,nrhs))*on_u(i,j)
            ustuy(i,j,k)=0.5_r8*                                        &
                           (u_stokes(i,j,k)+u_stokes(i,j-1,k))*         &
     &			   (u(i,j+1,k,nrhs)-u(i,j,k,nrhs))*om_v(i,j)
          END DO
        END DO
        DO i=Istr,Iend
          DO j=JstrV,Jend
            cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
            vstvx(i,j,k)=0.5_r8*                                        &
     &                   (v_stokes(i,j,k)+v_stokes(i,j-1,k))*           &
     &	                 (v(i,j+1,k,nrhs)-v(i,j,k,nrhs))*on_u(i,j)
            vstvy(i,j,k)=0.5_r8*                                        &
     &                     (v_stokes(i,j,k)+v_stokes(i,j-1,k))*         &
     &	                   (v(i,j+1,k,nrhs)-v(i,j,k,nrhs))*om_v(i,j)
          END DO
        END DO		
      END DO
      DO k=1,N(ng)		
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=(vstvx(i,j,k)*Hz(i,j,k)-                                &
     &           vstvx(i-1,j,k)*Hz(i-1,j,k))*on_u(i,j)  
            rustr3d(i,j,k)=rustr3d(i,j,k)+cff
          END DO
        END DO
        DO i=Istr,Iend
          DO j=JstrV,Jend 
            cff=(vstvy(i,j,k)*Hz(i,j,k)-                                &
     &           vstvx(i,j-1,k)*Hz(i,j-1,k))*om_v(i,j)
            rvstr3d(i,j,k)=rvstr3d(i,j,k)+cff
          END DO
        END DO
      END DO     
!
!  For a 3D application, compute associated 2D fields by taking the
!  vertical integral of 3D fields.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff5=0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          ubar_stokes(i,j)=cff5*u_stokes(i,j,1)
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            ubar_stokes(i,j)=ubar_stokes(i,j)+cff5*u_stokes(i,j,k)
          END DO
          cff3=0.25_r8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
          cff4=2.0_r8/(Dstp(i-1,j)+Dstp(i,j))
          ubar_stokes(i,j)=ubar_stokes(i,j)*cff4
        END DO
      END DO
      DO i=Istr,Iend
        DO j=JstrV,Jend
          cff5=0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          vbar_stokes(i,j)=cff5*v_stokes(i,j,1)
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
            vbar_stokes(i,j)=vbar_stokes(i,j)+cff5*v_stokes(i,j,k)
          END DO
          cff3=0.25_r8*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
          cff4=2.0_r8/(Dstp(i,j-1)+Dstp(i,j))
          vbar_stokes(i,j)=vbar_stokes(i,j)*cff4
        END DO
      END DO
#else
!PUT IN BAROTROPIC MOMENTUM BALANCE HERE
#endif
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  ubar_stokes)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  vbar_stokes)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rustr2d)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rvstr2d)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  zetat)
# ifdef SOLVE3D
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  u_stokes)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  v_stokes)
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  w_stokes)
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rustr3d)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rvstr3d)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rustr2d, rvstr2d, rulag2d, rvlag2d)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    zetat, ubar_stokes, vbar_stokes)
#  ifdef SOLVE3D
      CALL mp_exchange3d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rustr3d,  rvstr3d, rulag3d,  rvlag3d)
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    u_stokes, v_stokes)
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    w_stokes)
#  endif
# endif
      RETURN
      END SUBROUTINE wec_tile
#endif
      END MODULE wec_mod
