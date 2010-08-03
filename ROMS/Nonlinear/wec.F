#include "cppdefs.h"
#define WD_THORCHUR
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
!  34, 16-35.
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
# ifdef ATM_PRESS
     &                            FORCES(ng) % Pair,                    &
# endif
     &                            FORCES(ng) % Pwave_top,               &
     &                            FORCES(ng) % rollA,                   &
# ifdef SOLVE3D
#  ifdef DIAGNOSTICS_UV
     &                            DIAGS(ng) % DiaRU,                    &
     &                            DIAGS(ng) % DiaRV,                    &
#  endif
     &                            MIXING(ng) % rustr3d,                 &
     &                            MIXING(ng) % rvstr3d,                 &
#  ifdef WAVE_INDUCED_MIXING
     &                            MIXING(ng) % Akv,                     &
#  endif
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
     &                                  Hwave, Dwave, Lwave,            &
# ifdef ATM_PRESS
     &                                  Pair,                           &
# endif
     &                                  Pwave_top, rollA,               &
# ifdef SOLVE3D
#  ifdef DIAGNOSTICS_UV
     &                                  DiaRU, DiaRV,                   &
#  endif
     &                                  rustr3d,  rvstr3d,              &
#  ifdef WAVE_INDUCED_MIXING
     &                                  Akv,                            &
#  endif
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
#  ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Pwave_top(LBi:,LBj:)
      real(r8), intent(inout) :: rollA(LBi:,LBj:)
#  ifdef SOLVE3D
#   ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: rustr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvstr3d(LBi:,LBj:,:)
#   ifdef WAVE_INDUCED_MIXING
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
#   endif
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
      real(r8), intent(inout) :: zetat(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
#  endif
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
#  ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Pwave_top(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rollA(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
#   ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#   endif
      real(r8), intent(inout) :: rustr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rvstr3d(LBi:UBi,LBj:UBj,N(ng))
#   ifdef WAVE_INDUCED_MIXING
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
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
      real(r8) :: fac1, fac2, RB1, RB2
      real(r8), parameter :: alphaR=0.0_r8
      real(r8), parameter :: sinb=0.1_r8
      real(r8), parameter :: gammaw=0.31_r8
!     real(r8), parameter :: gammaw=0.4_r8
# ifdef BOTTOM_STREAMING 
      real(r8), parameter :: ks=0.03
      real(r8), parameter :: awd=3.0_r8
      real(r8), parameter :: KWDmax=200.0_r8
# endif
# ifdef WAVE_INDUCED_MIXING
      real(r8), parameter :: Cb=0.07_r8
# endif
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
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: WD
#  ifdef BOTTOM_STREAMING
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: KWD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owd
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EWD
#  endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigmat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zetaw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gamr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EER
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EBAR
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBX_bar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBY_bar
                 
# ifdef SOLVE3D
      real(r8) :: vi1, vi2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: oroller
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBXL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BBYL
#  ifdef BOTTOM_STREAMING
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BWDXL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BWDYL
#  endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huons
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvoms
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: bh
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: qsp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ku
      real(r8), dimension(IminS:ImaxS,N(ng)) :: kv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: VF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: shear1
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: shear2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ustux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ustuy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vstvx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vstvy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huston
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvston
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hustom
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvstom
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: BBX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: BBY      
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shearuv
# endif
# include "set_bounds.h"
      fac1=1.0_r8/dt(ng)
# ifdef ATM_PRESS
      fac2=1.0_r8/(g*rho0)
# endif
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
          zetaw(i,j)=-(0.5_r8*waveAA(i,j)*waven(i,j)/                   &
     &               SINH(2.0_r8*kD(i,j)))
# ifdef ATM_PRESS
          zetaw(i,j)=zetaw(i,j)-fac2*Pair(i,j)
# endif
!
!  Compute metrics for vertical roller distribution.
!
          oroller(i,j)=0.0_r8
          gamr(i,j)=MIN(0.707_r8*Dstp(i,j)/                             &
     &              (1.25_r8*Hwave(i,j)+eps),1.0_r8)
          DO k=1,N(ng)
            cff2=(1.0_r8+SCALARS(ng)%Cs_r(k))*gamr(i,j)
            oroller(i,j)=oroller(i,j)+Hz(i,j,k)*                        &
     &                   COSH(2.0_r8*pi*cff2)
          END DO
          oroller(i,j)=1.0_r8/(oroller(i,j)+eps)

# ifdef WD_THORGUZA
!
!  Calcualate wave dissipation using empirical parameters of 
!  Thornton and Guza, 1986
!
          cff1=0.707_r8*Hwave(i,j)
          cff2=1.0_r8/((gammaw*Dstp(i,j))**4)
          WD(i,j)=(3.0_r8*sqrt(pi)/16.0_r8)*g*sigmat(i,j)*              &
     &            (0.64_r8**3.0_r8/Dstp(i,j))*                          &
     &            (cff1**7.0_r8)*cff2    
# endif
# ifdef WD_THORCHUR
!
!  Calculate wave dissipation using empirical parameters of
!  Church and Thornton, 1993.
! 
          cff1=0.707_r8*Hwave(i,j)
          cff2=1.0_r8/(gammaw*Dstp(i,j))
          RB1=1.0_r8+tanh(8.0_r8*((cff1*cff2)-1.0_r8))
          RB2=1.0_r8-(1.0_r8+(cff1*cff2)**2.0_r8)**(-2.5_r8) 
          sigmat(i,j)=MIN(1.0_r8/Pwave_top(i,j),1.0_r8)
          WD(i,j)=(3.0_r8*sqrt(pi)/16.0_r8)*g*sigmat(i,j)*              &
     &            (0.64_r8**3.0_r8/Dstp(i,j))*                          &
     &            (cff1**3.0_r8)*RB1*RB2
# endif

# ifdef BOTTOM_STREAMING
!
!  Compute metrics for vertical bottom streaming distribution
!          
          owd(i,j)=0.0_r8
!  Bottom Orbital Velocity
          cff=0.25_r8*SQRT(2.0_r8)*sigma(i,j)*Hwave(i,j)/               &
     &                (SINH(kD(i,j))+eps)
          cff1=awd*0.09_r8*ks*(cff/                                     &
     &             (ks*sigma(i,j)))**0.82_r8
          KWD(i,j)=MIN(Dstp(i,j)/(cff1+eps),KWDmax)
          DO k=1,N(ng)
            cff2=(1.0_r8+SCALARS(ng)%Cs_r(N(ng)))-                      &
     &           (1.0_r8+SCALARS(ng)%Cs_r(k))   
            owd(i,j)=owd(i,j)+Hz(i,j,k)*                                &
     &               COSH(KWD(i,j)*cff2)
          END DO
          owd(i,j)=1.0_r8/(owd(i,j)+eps)
!  Wave friction factor (Soulsby, 1997)
          cff3=MIN(1.39_r8*(sigma(i,j)*(ks/30.0_r8)/                    &
     &               cff)**0.52_r8,0.2_r8)
!
!  Wave dissipation rate due to wave bottom drag Reniers et al. (2004b)
!
          EWD(i,j)=(0.5_r8/sqrt(pi))*cff3*(cff**3.0_r8) 
# endif       
!
!  Initialize depth independent arrays for summation
!
        END DO
        IF ((j.ge.Jstr) .and. (j.le.Jend)) THEN
          DO i=Istr,Iend
            BBX_bar(i,j)=0.0_r8
            BBY_bar(i,j)=0.0_r8
          END DO
        END IF
      END DO
!
!  Solve roller evolution equation for rollA.
!
      numits=30
      DO it=1,numits
!
!  Computer roller breaking source term and roller disspation sink term:
!  Roller breaking source (EBAR) Eqn 40 source.
!  Roller disspation rate (EER) Eqn 41.
!
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            EBAR(i,j)=-alphaR*WD(i,j)
            EER(i,j)=g*sinb*rollA(i,j)*sigma(i,j)/wavec(i,j) 
          END DO
        END DO
!
!  Time step roller action density.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff3=(ubar(i,j,nrhs)+wavenx(i,j)*owaven(i,j)*               &
     &           wavec(i,j))*on_u(i,j)
            cff1=MAX(cff3,0.0_r8)
            cff2=MIN(cff3,0.0_r8)
            FX(i,j)=cff1*rollA(i-1,j)+cff2*rollA(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff3=(vbar(i,j,nrhs)+waveny(i,j)*owaven(i,j)*               &
     &           wavec(i,j))*om_v(i,j)
            cff1=MAX(cff3,0.0_r8)
            cff2=MIN(cff3,0.0_r8)
            FE(i,j)=cff1*rollA(i,j-1)+cff2*rollA(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)/REAL(numits,r8)
            cff1=cff*(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))
            rollA(i,j)=rollA(i,j)-cff1
          END DO
        END DO
!
!  Add roller source / sink term.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)/REAL(numits,r8)
            rollA(i,j)=rollA(i,j)-cff*osigma(i,j)*(EBAR(i,j)+EER(i,j))
          END DO
        END DO
!
!  Call bc's.
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    rollA)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      rollA)
# endif
      END DO
!
!---------------------------------------------------------------------------
!  Compute Bernoulli's Head (BH), represented as the 
!  symbol cursive K in Eqn. 5.
!---------------------------------------------------------------------------
!
      J_LOOP1 : DO j=Jstr,Jend
!
! Compute k*u and k*v at rho pts.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            ku(i,k)=0.5_r8*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*             &
                    wavenx(i,j)
            kv(i,k)=0.5_r8*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*             &
                    waveny(i,j)
          END DO
        END DO
!
!  Compute gradients for shear d(k_x*u)/dz at cell faces.
!
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            shear1(i,k)=cff*((ku(i,k+1)-ku(i,k))+                       &
     &                       (kv(i,k+1)-kv(i,k)))
          END DO
        END DO
!
! Applying boundary conditions.
!
        DO i=Istr,Iend
          shear1(i,N(ng))=shear1(i,N(ng)-1)
          shear1(i,0)=shear1(i,1)
        END DO
!
! Compute second derivative d(k_x*u)^2/dz^2
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            cff=0.5_r8/(z_w(i,j,k)-z_w(i,j,k-1))
            shear2(i,k)=cff*(shear1(i,k)-shear1(i,k-1))
          END DO
        END DO
!
!  Perform integration on the shear2 term calcuated above.
!
        DO i=Istr,Iend
          cff=Hz(i,j,1)
          fac2=(1.0_r8+SCALARS(ng)%Cs_r(N(ng)))-                        &
     &         (1.0_r8+SCALARS(ng)%Cs_r(1))
          vi1=cff*shear2(i,1)*SINH(2.0_r8*kD(i,j)*fac2)
          DO k=2,N(ng)
            cff=Hz(i,j,k)
            fac2=(1.0_r8+SCALARS(ng)%Cs_r(N(ng)))-                      &
     &           (1.0_r8+SCALARS(ng)%Cs_r(k))
            vi1=vi1+cff*shear2(i,k)*SINH(2.0_r8*kD(i,j)*fac2)             
          END DO
!
!  Calculate Bernoulli's Head at the surface
!
          bh(i,j)=0.25_r8*vi1*sigma(i,j)*waveAA(i,j)/                   & 
     &                    (waven(i,j)*(SINH(kD(i,j)))**2)
        END DO
!
!---------------------------------------------------------------------------
!  Compute quasi-static components of pressure terms (i.e. Cursive P)
!  in Eqn. 9.
!---------------------------------------------------------------------------
!
        DO i=Istr,Iend
          fac2=1.0_r8+SCALARS(ng)%Cs_r(1)
          cff=Hz(i,j,1)
          vi2=cff*shear2(i,1)*COSH(2.0_r8*kD(i,j)*fac2)
          DO k=2,N(ng)
             fac2=1.0_r8+SCALARS(ng)%Cs_r(k)
             cff=Hz(i,j,k)
             vi2=vi2+cff*shear2(i,k)*COSH(2.0_r8*kD(i,j))*fac2	  
          END DO
          cff1=0.5_r8*waveAA(i,j)*osigma(i,j)
          cff2=TANH(kD(i,j))/SINH(2.0_r8*kD(i,j))
          qsp(i,j)=cff1*(cff2*(-shear1(i,N(ng))+                        &
     &                         COSH(2.0_r8*kD(i,j))*shear1(i,0)+        &
     &                         vi2)-                                    &
     &                   2.0_r8*waven(i,j)*TANH(kD(i,j))*               &
     &                   (ku(i,N(ng))+kv(i,N(ng))))
        END DO
      END DO J_LOOP1
!
!---------------------------------------------------------------------------
! Evaluate geopotential function (Eqn. 33, Uchiyama et al., 2010)
!---------------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zetat(i,j)=-Dstp(i,j)*                                        &
     &                      (-g*zetaw(i,j)-(g*qsp(i,j)-bh(i,j)))
        END DO
      END DO
!
!---------------------------------------------------------------------------
! Stokes velocities.
!---------------------------------------------------------------------------
      DO k=1,N(ng)
!
!  Compute u-stokes velocities, Eqn. 2.
!
        fac2=1.0_r8+SCALARS(ng)%Cs_r(k)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
             cff2=(waveE(i-1,j)+waveE(i,j))
             cff2=cff2+(rollA(i-1,j)*sigma(i-1,j)+                      &
     &                  rollA(i,j)*sigma(i,j))
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
!  Compute v-stokes velocity, Eqn. 2.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
             cff2=(waveE(i,j-1)+waveE(i,j))
             cff2=cff2+(rollA(i,j-1)*sigma(i,j-1)+                      &
     &                  rollA(i,j)*sigma(i,j))
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
      END DO
!
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  u_stokes)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  v_stokes)
!	  
!  Compute vertical stokes velocity, Eqn. 31.
!
      DO k=1,N(ng)
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
     &                       Hvoms(i,j+1)-Hvoms(i,j))
          END DO
        END DO
      END DO
!
!---------------------------------------------------------------------------
! Compute Vertical Vortex Force terms, denoted as K in Eqn. 5.
!---------------------------------------------------------------------------
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          DO k=1,N(ng)-1
!            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            cff=1.0_r8
            cff1=0.5_r8*(u_stokes(i,j,k+1)+u_stokes(i,j,k))
            cff2=0.5_r8*(u_stokes(i+1,j,k+1)+u_stokes(i+1,j,k))
            cff3=0.5_r8*(v_stokes(i,j,k+1)+v_stokes(i,j,k))
            cff4=0.5_r8*(v_stokes(i,j+1,k+1)+v_stokes(i,j+1,k))
            shearuv(i,j,k)=cff*                                         &
     &                      (cff1*(u(i  ,j,k+1,nrhs)-u(i  ,j,k,nrhs))+  &
     &                       cff2*(u(i+1,j,k+1,nrhs)-u(i+1,j,k,nrhs)))+ &
     &                     cff*                                         &
     &                      (cff3*(v(i,j  ,k+1,nrhs)-v(i,j  ,k,nrhs))+  &
     &                       cff4*(v(i,j+1,k+1,nrhs)-v(i,j+1,k,nrhs)))
          END DO
          shearuv(i,j,0)=shearuv(i,j,1)
          shearuv(i,j,N(ng))=shearuv(i,j,N(ng)-1)
        END DO
      END DO
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
            VF(i,k)=0.5_r8*((shearuv(i  ,j,k)+shearuv(i  ,j,k-1))-      &
     &                      (shearuv(i-1,j,k)+shearuv(i-1,j,k-1)))
            rustr3d(i,j,k)=cff*VF(i,k)*on_u(i,j)
          END DO
        END DO
        IF (j.ge.JstrV) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
              VF(i,k)=0.5_r8*((shearuv(i,j  ,k)+shearuv(i,j  ,k-1))-    &
     &                        (shearuv(i,j-1,k)+shearuv(i,j-1,k-1)))
              rvstr3d(i,j,k)=cff*VF(i,k)*om_v(i,j)
            END DO
          END DO
        END IF
      END DO J_LOOP
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
            cff4=((1.0_r8-alphaR)*WD(i,j)+EER(i,j))*osigma(i,j)
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
            cff=0.5_r8*(BBXL(i  ,j)*Hz(i  ,j,k)+                        &
     &                  BBXL(i-1,j)*Hz(i-1,j,k))*                       &
     &                  om_u(i,j)*on_u(i,j)
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
# ifdef DIAGNOSTICS_UV
!             DiaRU(i,j,k,nrhs,M3hrad)=cff
# endif
          END DO
        END DO
!
! Compute contribution to V-momentum
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(BBYL(i,j  )*Hz(i,j  ,k)+                        &
     &                  BBYL(i,j-1)*Hz(i,j-1,k))*                       &
     &                  om_v(i,j)*on_v(i,j)
            rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
# ifdef DIAGNOSTICS_UV
!           DiaRV(i,j,k,nrhs,M3hrad)=cff
# endif
          END DO
        END DO

# ifdef BOTTOM_STREAMING
!
! Compute bottom streaming based acceleration terms
!
        fac2=SCALARS(ng)%Cs_r(N(ng))-                                   &
     &       SCALARS(ng)%Cs_r(k)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            cff2=COSH(fac2*KWD(i,j))
            cff3=EWD(i,j)*osigma(i,j)
            BWDXL(i,j)=cff2*cff3*wavenx(i,j)*owd(i,j)
            BWDYL(i,j)=cff2*cff3*waveny(i,j)*owd(i,j) 
          END DO
        END DO
!
! Compute contribution to U-momentum
! 
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.5_r8*(BWDXL(i  ,j)*Hz(i  ,j,k)+                       &
     &                  BWDXL(i-1,j)*Hz(i-1,j,k))*                      &
     &                  om_u(i,j)*on_u(i,j)
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
# ifdef DIAGNOSTICS_UV
!             DiaRU(i,j,k,nrhs,M3hrad)=cff
# endif
          END DO
        END DO
!
! Compute contribution to V-momentum
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(BWDYL(i,j  )*Hz(i,j  ,k)+                       &
     &                  BWDYL(i,j-1)*Hz(i,j-1,k))*                      &
     &                  om_v(i,j)*on_v(i,j)
            rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
# ifdef DIAGNOSTICS_UV
!           DiaRV(i,j,k,nrhs,M3hrad)=cff
# endif
          END DO
        END DO
# endif       
# ifdef WAVE_INDUCED_MIXING
! 
! Compute contribution of wave breaking induced mixing to Akv
!
        fac2=1.0_r8+SCALARS(ng)%Cs_r(k)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff2=fac2*gamr(i,j)
            cff3=COSH(2.0_r8*pi*cff2)
            cff4=(EER(i,j)*sigma(i,j))**(1.0_r8/3.0_r8)
            Akv(i,j,k)=Akv(i,j,k)+                                      &
     &                 Cb*cff4*SQRT(2.0_r8)*Hwave(i,j)*Dstp(i,j)*       &
     &                 cff3*oroller(i,j)
          END DO
        END DO   
# endif
!---------------------------------------------------------------------------
! Contribution of a term corresponding to product of 
! Stokes and Eulerian Velocity Eqn. 26 and 27.
!---------------------------------------------------------------------------
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            Huston(i,j)=cff*on_u(i,j)*u_stokes(i,j,k)
            Hvston(i,j)=0.25_r8*cff*on_u(i,j)*                          &
     &                  (v_stokes(i  ,j  ,k)+                           &
     &                   v_stokes(i  ,j+1,k)+                           &
     &                   v_stokes(i-1,j  ,k)+                           &
     &                   v_stokes(i-1,j+1,k))
          END DO
          DO i=IstrU-1,Iend
            UFx(i,j)=0.5_r8*(u(i  ,j  ,k,nrhs)+                         &
                             u(i+1,j  ,k,nrhs))
            VFx(i,j)=0.5_r8*(v(i  ,j  ,k,nrhs)+                         &
     &                       v(i  ,j+1,k,nrhs))           
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
            Hustom(i,j)=cff*0.25_r8*om_v(i,j)*                          &
     &                  (u_stokes(i  ,j  ,k)+                           &
     &                   u_stokes(i+1,j  ,k)+                           &
     &                   u_stokes(i  ,j-1,k)+                           &
     &                   u_stokes(i+1,j-1,k))
            Hvstom(i,j)=cff*om_v(i,j)*v_stokes(i,j,k)
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
            UFe(i,j)=0.5_r8*(u(i+1,j  ,k,nrhs)+                         &
     &                       u(i  ,j  ,k,nrhs))
            VFe(i,j)=0.5_r8*(v(i  ,j  ,k,nrhs)+                         &
     &                       v(i  ,j+1,k,nrhs))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=UFx(i,j)-UFx(i-1,j)
            cff2=VFx(i,j)-VFx(i-1,j)
            cff=Huston(i,j)*cff1+Hvston(i,j)*cff2
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
!
!  Conversion to m2/s2.
!
            cff1=0.25_r8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
            rustr3d(i,j,k)=rustr3d(i,j,k)*cff1
          END DO
        END DO
        DO i=Istr,Iend
          DO j=JstrV,Jend 
            cff1=UFe(i,j)-UFe(i,j-1)
            cff2=VFe(i,j)-VFe(i,j-1)
            cff=Hustom(i,j)*cff1+Hvstom(i,j)*cff2 
            rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
!
!  Conversion to m2/s2.
!
            cff1=0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            rvstr3d(i,j,k)=rvstr3d(i,j,k)*cff1
          END DO
        END DO
      END DO K_LOOP
!
!  For a 3D application, compute associated 2D fields by taking the
!  vertical integral of 3D fields.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff5=0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          ubar_stokes(i,j)=cff5*u_stokes(i,j,1)
          rustr2d(i,j)=cff5*rustr3d(i,j,1)
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))
            ubar_stokes(i,j)=ubar_stokes(i,j)+cff5*u_stokes(i,j,k)
            rustr2d(i,j)=rustr2d(i,j)+cff5*rustr3d(i,j,k)
          END DO
          cff3=0.25_r8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
          cff4=2.0_r8/(Dstp(i-1,j)+Dstp(i,j))
          rustr2d(i,j)=rustr2d(i,j)*cff3*cff4
          ubar_stokes(i,j)=ubar_stokes(i,j)*cff4
        END DO
      END DO
      DO i=Istr,Iend
        DO j=JstrV,Jend
          cff5=0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          vbar_stokes(i,j)=cff5*v_stokes(i,j,1)
          rvstr2d(i,j)=cff5*rvstr3d(i,j,1)
          DO k=2,N(ng)
            cff5=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))
            vbar_stokes(i,j)=vbar_stokes(i,j)+cff5*v_stokes(i,j,k)
            rvstr2d(i,j)=rvstr2d(i,j)+cff5*rvstr3d(i,j,k)
          END DO
          cff3=0.25_r8*(pm(i,j-1)+pm(i,j))*(pn(i,j-1)+pn(i,j))
          cff4=2.0_r8/(Dstp(i,j-1)+Dstp(i,j))
          vbar_stokes(i,j)=vbar_stokes(i,j)*cff4
          rvstr2d(i,j)=rvstr2d(i,j)*cff3*cff4
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
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  w_stokes)
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rustr3d)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rvstr3d)
#  ifdef WAVE_INDUCED_MIXING
#   if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        Akv)
#   endif
#  endif
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
#   ifdef WAVE_INDUCED_MIXING
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Akv)
#   endif
#  endif
# endif
      RETURN
      END SUBROUTINE wec_tile
#endif
      END MODULE wec_mod
