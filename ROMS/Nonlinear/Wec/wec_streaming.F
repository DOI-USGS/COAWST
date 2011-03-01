#include "cppdefs.h"
!#define BOTTOM_STREAMING_YU
#undef  BOTTOM_STREAMING_XU_BOWEN
#undef  SURFACE_STREAMING
      MODULE wec_streaming_mod
#if defined SOLVE3D && defined BOTTOM_STREAMING
!
!svn $Id: wec_streaming.F 1428 2008-03-12 13:07:21Z jcwarner $
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
      PUBLIC  :: wec_streaming
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_streaming (ng, tile)
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
      CALL wec_streaming_tile (ng, tile,    LBi, UBi, LBj, UBj, N(ng),  &
     &                         IminS, ImaxS, JminS, JmaxS,              &
# ifdef SOLVE3D
     &                         nrhs(ng),                                &
# endif
     &                         GRID(ng) % angler,                       &
     &                         GRID(ng) % h,                            &
# ifdef SOLVE3D
     &                         GRID(ng) % Hz,                           &
# endif
     &                         GRID(ng) % om_u,                         &
     &                         GRID(ng) % om_v,                         &
     &                         GRID(ng) % on_u,                         &
     &                         GRID(ng) % on_v,                         &
# ifdef BOTTOM_STREAMING_XU_BOWEN
     &                         GRID(ng) % fomn,
# endif
     &                         OCEAN(ng) % zeta,                        &
     &                         FORCES(ng) % Hwave,                      &
     &                         FORCES(ng) % Dwave,                      &
     &                         FORCES(ng) % Lwave,                      &
     &                         FORCES(ng) % Pwave_top,                  &
# ifdef WAVES_OCEAN
     &                         FORCES(ng) % Dissip_fric,                &
# endif
# ifdef DIAGNOSTICS_UV
     &                         DIAGS(ng) % DiaRU,                       &
     &                         DIAGS(ng) % DiaRV,                       &
# endif
# ifdef BOTTOM_STREAMING_XU_BOWEN
     &                         MIXING(ng) % visc2_r,                    &
# endif
# ifdef SURFACE_STREAMING
     &                         MIXING(ng) & Akv,                        &
# endif
     &                         MIXING(ng) % rustr3d,                    &
     &                         MIXING(ng) % rvstr3d)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 21)
# endif
      RETURN
      END SUBROUTINE wec_streaming
!
!***********************************************************************
      SUBROUTINE wec_streaming_tile (ng, tile, LBi, UBi, LBj, UBj, UBk, &
     &                               IminS, ImaxS, JminS, JmaxS,        &
# ifdef SOLVE3D
     &                               nrhs,                              &
# endif
     &                               angler, h,                         &
# ifdef SOLVE3D
     &                               Hz,                                &
# endif
     &                               om_u, om_v, on_u, on_v,            &
# ifdef BOTTOM_STREAMING_XU_BOWEN
     &                               fomn,                              &
# endif
     &                               zeta,                              &
     &                               Hwave, Dwave, Lwave,               &
     &                               Pwave_top,                         &
# ifdef WAVES_OCEAN
     &                               Dissip_fric,                       &
# endif
# ifdef DIAGNOSTICS_UV
     &                               DiaRU, DiaRV,                      &
# endif
# ifdef BOTTOM_STREAMING_XU_BOWEN
     &                               visc2_r,                           &
# endif
# ifdef SURFACE_STREAMING
     &                               Akv,                               &
# endif
     &                               rustr3d,  rvstr3d)
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
      integer, intent(in) :: nrhs
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
#  ifdef BOTTOM_STREAMING_XU_BOWEN
      real(r8), intent(in) :: fomn(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_top(LBi:,LBj:)
#  ifdef WAVES_OCEAN
      real(r8), intent(in) :: Dissip_fric(LBi:,LBj:)
#  endif
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif
#  ifdef BOTTOM_STREAMING_XU_BOWEN
      real(r8), intent(in) :: visc2_r(LBi:,LBj:)
#  endif
#  ifdef SURFACE_STREAMING
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
#  endif
      real(r8), intent(inout) :: rustr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvstr3d(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
#  endif
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
#  ifdef BOTTOM_STREAMING_XU_BOWEN
      real(r8), intent(in) :: fomn(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pwave_top(LBi:UBi,LBj:UBj)
#  ifdef WAVES_OCEAN
      real(r8), intent(in) :: Dissip_fric(LBi:UBi,LBj:UBj)
#  endif
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif
#  ifdef BOTTOM_STREAMING_XU_BOWEN
      real(r8), intent(in) :: visc2_r(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SURFACE_STREAMING
      real(r8), intent(in) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(inout) :: rustr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rvstr3d(LBi:UBi,LBj:UBj,N(ng))
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
      integer :: i, j, k
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: fac1, fac2
      real(r8), parameter :: ks=0.03
      real(r8), parameter :: awd=3.0_r8
      real(r8), parameter :: KWDmax=200.0_r8
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
# ifdef BOTTOM_STREAMING
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: KWD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owd
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EWD
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osigma
# ifdef BOTTOM_STREAMING_XU_BOWEN
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: WaveCo
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Kobeta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Fosigma
# endif
# ifdef SURFACE_STREAMING
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Surst
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Sstopx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Sstopy
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BWDXL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BWDYL
# include "set_bounds.h"
      fac1=1.0_r8/dt(ng)
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
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff)
          waveny(i,j)=waven(i,j)*SIN(cff)
          sigma(i,j)=MIN(SQRT(g*waven(i,j)*TANH(waven(i,j)*Dstp(i,j))), &
     &                   2.0_r8)
          osigma(i,j)=1.0_r8/sigma(i,j)
!
!  Compute wave celerity and nonlinear water depth
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j)+eps,kDmax)

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
!  Hold this constant for now.  Need to add logic for 
!  zo if no bbl on, or zo_apparent if use a bbl.
         cff3=MIN(1.39_r8*(sigma(i,j)*(ks/30.0_r8)/                    &
     &               cff)**0.52_r8,0.2_r8)
!
!  Wave dissipation rate due to wave bottom drag Reniers et al. (2004b)
!
#  ifdef WAVES_OCEAN
          EWD(i,j)=Dissip_fric(i,j)
#  else
          EWD(i,j)=(0.5_r8/sqrt(pi))*cff3*(cff**3.0_r8)
#  endif
# endif
# ifdef BOTTOM_STREAMING_XU_BOWEN
          Kobeta(i,j)=waven(i,j)/((SQRT(0.5_r8*sigma(i,j)/              &
     &                visc2_r(i,j)))+eps)
          WaveCo(i,j)=(0.5_r8*(0.5_r8*Hwave(i,j)*                       &
     &                sigma(i,j))**2.0_r8)/                             &
     &                (4.0_r8*(SINH(kD(i,j)))**2.0_r8+eps)
          Fosigma(i,j)=fomn(i,j)*osigma(i,j)
# endif
# ifdef SURFACE_STREAMING
          Surst(i,j)=0.25_r8*(Hwave(i,j)**2.0_r8)*g*                    &
     &               (waven(i,j)**2.0_r8)*osigma(i,j)          
# endif
        END DO
      END DO
!
! Compute bottom streaming based acceleration terms
!
      K_LOOP : DO k=1,N(ng)
# ifdef BOTTOM_STREAMING
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
# endif
# ifdef BOTTOM_STREAMING_XU_BOWEN
        fac2=1+SCALARS(ng)%Cs_r(k)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            cff =Dstp(i,j)-(Dstp(i,j)*fac2)
            cff3=Beta(i,j)*cff
            cff1=(2.0_r8*cff3)*                                         &
     &           EXP(-cff3)*SIN(cff3)+                                  &
     &           2.0_r8*EXP(-cff3)*COS(cff3)-1.0_r8-                    &
     &           EXP(-2.0_r8*cff3)
            cff2=SINH(2.0_r8*waven(i,j)*cff)
            BWDXL(i,j)=cff1*Kobeta(i,j)*WaveCo(i,j)
            BWDYL(i,j)=cff2*Fosigma(i,j)*WaveCo(i,j)
          END DO
        END DO
# endif
!
! Compute contribution to U-momentum
! 
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.5_r8*(BWDXL(i  ,j)*Hz(i  ,j,k)+                       &
   &                    BWDXL(i-1,j)*Hz(i-1,j,k))
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
# ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3bstm)=cff
# endif
          END DO
        END DO
!
! Compute contribution to V-momentum
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(BWDYL(i,j  )*Hz(i,j  ,k)+                       &
   &                    BWDYL(i,j-1)*Hz(i,j-1,k))
            rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
# ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3bstm)=cff
# endif
          END DO
        END DO
      END DO K_LOOP
!
# ifdef SURFACE_STREAMING
!
! Compute surface streaming based acceleration terms
!
         DO j=Jstr-1,Jend+1
           DO i=Istr-1,Iend+1
             IF (k.eq.N(ng)) THEN
               Sstopx(i,j)= Akv(i,j,k)*Surst(i,j)*wavenx(i,j)
               Sstopy(i,j)= Akv(i,j,k)*Surst(i,j)*waveny(i,j)
             ELSE
               Sstopx(i,j)= 0.0_r8
               Sstopy(i,j)= 0.0_r8
             ENDIF
           END DO
         END DO
!
! Compute contribution to U-momentum
! 
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.5_r8*(Sstopx(i  ,j)*Hz(i  ,j,k)+                      &
     &                  Sstopx(i-1,j)*Hz(i-1,j,k))
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3sstm)=cff
#  endif
          END DO
        END DO
!
! Compute contribution to V-momentum
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(Sstopy(i,j  )*Hz(i,j  ,k)+                      &
     &                  Sstopy(i,j-1)*Hz(i,j-1,k))
            rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3sstm)=cff
#  endif
          END DO
        END DO
# endif
!
!  Apply boundary conditions.
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rustr3d)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rvstr3d)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rustr3d,  rvstr3d)
# endif
      RETURN
      END SUBROUTINE wec_streaming_tile
#endif
      END MODULE wec_streaming_mod
