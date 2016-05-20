#include "cppdefs.h"
#define BOTTOM_STREAMING_YU
#undef  BOTTOM_STREAMING_XU_BOWEN
      MODULE wec_streaming_mod
#if defined SOLVE3D && (defined BOTTOM_STREAMING || \
    defined SURFACE_STREAMING)
!
!svn $Id: wec_streaming.F 1428 2008-03-12 13:07:21Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!                                                   Nirnimesh Kumar    !
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
      USE mod_stepping
# if defined DIAGNOSTICS_UV
      USE mod_diags
# endif
# if defined VEGETATION && defined VEG_STREAMING 
      USE mod_vegarr
      USE vegetation_stream_mod, ONLY : vegetation_stream_cal
# endif 
!
      integer, intent(in) :: ng, tile
# include "tile.h"

# if defined VEGETATION && defined VEG_STREAMING 
      CALL vegetation_stream_cal (ng, tile)
# endif 

# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 21)
# endif
      CALL wec_streaming_tile (ng, tile, LBi, UBi, LBj, UBj, N(ng),     &
     &                         IminS, ImaxS, JminS, JmaxS,              &
# ifdef SOLVE3D
     &                         nrhs(ng),                                &
# endif
     &                         GRID(ng) % angler,                       &
     &                         GRID(ng) % om_u,                         &
     &                         GRID(ng) % om_v,                         &
     &                         GRID(ng) % on_u,                         &
     &                         GRID(ng) % on_v,                         &
# ifdef SOLVE3D
     &                         GRID(ng) % Hz,                           &
     &                         GRID(ng) % z_r,                          &
     &                         GRID(ng) % z_w,                          &
# endif
# ifdef BOTTOM_STREAMING_XU_BOWEN
     &                         GRID(ng) % f,                            &
# endif
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
# if defined BOTTOM_STREAMING_XU_BOWEN || defined SURFACE_STREAMING
     &                         MIXING(ng) % Akv,                        &
# endif
# ifdef BOTTOM_STREAMING
     &                         MIXING(ng) % rubst2d,                    &
     &                         MIXING(ng) % rvbst2d,                    &
# endif
# ifdef SURFACE_STREAMING
     &                         MIXING(ng) % russt2d,                    &
     &                         MIXING(ng) % rvsst2d,                    &
# endif
# if defined VEGETATION && defined VEG_STREAMING
     &                         VEG(ng) % BWDXL_veg,                     &
     &                         VEG(ng) % BWDYL_veg,                     &
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
     &                               angler,                            &
     &                               om_u, om_v, on_u, on_v,            &
# ifdef SOLVE3D
     &                               Hz, z_r, z_w,                      &
# endif
# ifdef BOTTOM_STREAMING_XU_BOWEN
     &                               f,                                 &
# endif
     &                               Hwave, Dwave, Lwave,               &
     &                               Pwave_top,                         &
# ifdef WAVES_OCEAN
     &                               Dissip_fric,                       &
# endif
# ifdef DIAGNOSTICS_UV
     &                               DiaRU, DiaRV,                      &
# endif
# if defined BOTTOM_STREAMING_XU_BOWEN || defined SURFACE_STREAMING
     &                               Akv,                               &
# endif
# ifdef BOTTOM_STREAMING
     &                               rubst2d, rvbst2d,                  &
# endif
# ifdef SURFACE_STREAMING
     &                               russt2d, rvsst2d,                  &
# endif
# if defined VEGETATION && defined VEG_STREAMING
     &                               BWDXL_veg, BWDYL_veg,              &
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
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  endif
#  ifdef BOTTOM_STREAMING_XU_BOWEN
      real(r8), intent(in) :: f(LBi:,LBj:)
#  endif
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
#  if defined BOTTOM_STREAMING_XU_BOWEN || defined SURFACE_STREAMING
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
#  endif
#  ifdef BOTTOM_STREAMING
      real(r8), intent(inout) :: rubst2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvbst2d(LBi:,LBj:)
#  endif
#  ifdef SURFACE_STREAMING
      real(r8), intent(inout) :: russt2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvsst2d(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: rustr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvstr3d(LBi:,LBj:,:)
#  if defined VEGETATION && defined VEG_STREAMING 
      real(r8), intent(inout) :: BWDXL_veg(LBi:,LBj:,:)
      real(r8), intent(inout) :: BWDYL_veg(LBi:,LBj:,:)
#  endif
# else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
#  endif
#  ifdef BOTTOM_STREAMING_XU_BOWEN
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
#  endif
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
#  if defined BOTTOM_STREAMING_XU_BOWEN || defined SURFACE_STREAMING
      real(r8), intent(in) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
#  ifdef BOTTOM_STREAMING
      real(r8), intent(inout) :: rubst2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvbst2d(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SURFACE_STREAMING
      real(r8), intent(inout) :: russt2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rvsst2d(LBi:UBi,LBj:UBj)
#  endif
#  if defined VEGETATION && defined VEG_STREAMING 
      real(r8), intent(inout) :: BWDXL_veg(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: BWDYL_veg(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: rustr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rvstr3d(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: fac2, sqrt2
# ifdef BOTTOM_STREAMING_XU_BOWEN
      real(r8) ::  Beta, fosigma
# endif
      real(r8), parameter :: ks=0.03_r8
      real(r8), parameter :: awd=1.0_r8
      real(r8), parameter :: KWDmax=200.0_r8
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: oDstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
# ifdef BOTTOM_STREAMING_YU
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: KWD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owd
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EWD
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osigma
# ifdef SURFACE_STREAMING
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Surst
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Sstopx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Sstopy
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BWDXL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BWDYL
# include "set_bounds.h"
!
      sqrt2=SQRT(2.0_r8)
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!
!  Compute total depth
!
          Dstp(i,j)=z_w(i,j,N(ng))-z_w(i,j,0)
          oDstp(i,j)=1.0_r8/Dstp(i,j)
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

# ifdef BOTTOM_STREAMING_YU
!
!  Compute metrics for vertical bottom streaming distribution
!          
          owd(i,j)=0.0_r8
!  Bottom Orbital Velocity
          cff=0.25_r8*sqrt2*sigma(i,j)*Hwave(i,j)/                      &
     &        (SINH(kD(i,j))+eps)
          cff1=awd*0.09_r8*ks*(cff/                                     &
     &         (ks*sigma(i,j)))**0.82_r8
          KWD(i,j)=MIN(Dstp(i,j)/(cff1+eps),KWDmax)
          DO k=1,N(ng)
            cff2=(z_r(i,j,N(ng))-z_r(i,j,k))*oDstp(i,j)
            owd(i,j)=owd(i,j)+                                          &
     &               Hz(i,j,k)*COSH(KWD(i,j)*cff2)
          END DO
          owd(i,j)=1.0_r8/owd(i,j)
!
!  Wave friction factor (based on Soulsby, 1997).
!  Hold this constant for now.  Need to add logic for 
!  zo if no bbl on, or zo_apparent if use a bbl.
         cff3=MIN(1.39_r8*(sigma(i,j)*(ks/30.0_r8)/                     &
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
# ifdef SURFACE_STREAMING
          Surst(i,j)=0.25_r8*(Hwave(i,j)**2.0_r8)*g*                    &
     &               (waven(i,j)**2.0_r8)*osigma(i,j)          
# endif
        END DO
      END DO
!
!  Initialize depth-avg forcing terms for subsequent computation.
!  For now, these are used for the diagnostics.  This logic will change
!  when 2D method is implemented.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
# ifdef BOTTOM_STREAMING
          rubst2d(i,j)=0.0_r8
          rvbst2d(i,j)=0.0_r8
# endif
# ifdef SURFACE_STREAMING
          russt2d(i,j)=0.0_r8
          rvsst2d(i,j)=0.0_r8
# endif
        END DO
      END DO
!
! Compute bottom streaming based acceleration terms
!
      K_LOOP : DO k=1,N(ng)
# ifdef BOTTOM_STREAMING_YU
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            fac2=(z_r(i,j,N(ng))-z_r(i,j,k))*oDstp(i,j)
            cff2=COSH(fac2*KWD(i,j))
            cff3=EWD(i,j)*osigma(i,j)
            BWDXL(i,j)=cff2*cff3*wavenx(i,j)*owd(i,j)
            BWDYL(i,j)=cff2*cff3*waveny(i,j)*owd(i,j)
          END DO
        END DO
# endif
# ifdef BOTTOM_STREAMING_XU_BOWEN
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            fac2=(z_r(i,j,k)-z_w(i,j,0))*oDstp(i,j)
            cff=0.0313_r8*(Hwave(i,j)*sigma(i,j)/                     &
     &                     SINH(kD(i,j)))**2.0_r8
            Beta=MAX(SQRT(0.5_r8*sigma(i,j)/                          &
     &               MAX(Akv(i,j,k),eps)),40.0_r8)
            fosigma=f(i,j)*osigma(i,j)
            cff1=(Dstp(i,j)-Dstp(i,j)*fac2)
            cff2=(2.0_r8*Beta*cff1)*                                  &
     &           EXP(-cff1*Beta)*SIN(cff1*Beta)+                      &
     &           2.0_r8*EXP(-cff1*Beta)*COS(cff1*Beta)-1.0_r8-        &
     &           EXP(-2.0_r8*cff1*Beta)
            cff2=MIN(cff2,0.0_r8)
            cff3=SINH(2.0_r8*waven(i,j)*cff1)
            BWDXL(i,j)=cff*cff2*waven(i,j)/Beta
            BWDYL(i,j)=fosigma*cff*cff3
          END DO
        END DO
!
! Bottom Streaming in Lentz et al. (2008)
!          DO j=Jstr-1,Jend+1
!            DO i=Istr-1,Iend+1
!              fac2=(z_r(i,j,k)-z_w(i,j,0))*oDstp(i,j)
!              cff=0.0625_r8*waven(i,j)*(Hwave(i,j)*sigma(i,j)/          &
!     &                       SINH(kD(i,j)))**2.0_r8
!              Beta=MAX(SQRT(0.5_r8*sigma(i,j)/                          &
!     &                 MAX(Akv(i,j,k),eps)),40.0_r8)
!              cff1=(Dstp(i,j)-Dstp(i,j)*fac2)
!              cff2=(-Beta*cff1*SIN(cff1*Beta)+                          &
!     &               Beta*cff1*COS(cff1*Beta)-                          &
!     &               COS(Beta*cff1))*EXP(-cff1*Beta)-                   &
!     &               EXP(-2.0_r8*cff1*Beta)
!              cff3=SINH(2.0_r8*waven(i,j)*cff1)
!              fosigma=f(i,j)*osigma(i,j)
!              BWDXL(i,j)=cff*cff2*waven(i,j)/Beta
!              BWDYL(i,j)=fosigma*cff*cff3*0.0_r8
!          END DO
!        END DO

# endif
# if defined VEGETATION && defined VEG_STREAMING
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            BWDXL(i,j)=BWDXL_veg(i,j,k)
            BWDYL(i,j)=BWDYL_veg(i,j,k)
          END DO
        END DO
# endif
!
! Compute contribution to U-momentum
! 
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=0.5_r8*(Hz(i-1,j,k)+                                   &
     &                   Hz(i  ,j,k))
            cff=0.5_r8*(BWDXL(i  ,j)*Hz(i  ,j,k)+                       &
     &                  BWDXL(i-1,j)*Hz(i-1,j,k))
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
            rubst2d(i,j)=rubst2d(i,j)+cff*cff1*om_u(i,j)*on_u(i,j)
# ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3bstm)=cff*om_u(i,j)*on_u(i,j)
# endif
          END DO
        END DO
!
! Compute contribution to V-momentum
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=0.5_r8*(Hz(i,j  ,k)+                                   &
     &                   Hz(i,j-1,k))
            cff=0.5_r8*(BWDYL(i,j  )*Hz(i,j  ,k)+                       &
     &                  BWDYL(i,j-1)*Hz(i,j-1,k))
            rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
            rvbst2d(i,j)=rvbst2d(i,j)+cff*cff1*om_v(i,j)*on_v(i,j)
# ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3bstm)=cff*om_v(i,j)*on_v(i,j)
# endif
          END DO
        END DO
!
# ifdef SURFACE_STREAMING
!
! Compute surface streaming based acceleration terms
!
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            IF (k.eq.N(ng)) THEN
              Sstopx(i,j)= 0.5_r8*(Akv(i,j,k)+Akv(i,j,k-1))*            &
     &                             Surst(i,j)*wavenx(i,j)
              Sstopy(i,j)= 0.5_r8*(Akv(i,j,k)+Akv(i,j,k-1))*            &
     &                             Surst(i,j)*waveny(i,j)
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
            cff1=0.5_r8*(Hz(i-1,j,k)+                                   &
     &                   Hz(i  ,j,k))
            cff=0.5_r8*(Sstopx(i  ,j)*Hz(i  ,j,k)+                      &
     &                  Sstopx(i-1,j)*Hz(i-1,j,k))
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
            russt2d(i,j)=russt2d(i,j)+cff*cff1*om_u(i,j)*on_u(i,j)
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3sstm)=cff*om_u(i,j)*on_u(i,j)
#  endif
          END DO
        END DO
!
! Compute contribution to V-momentum
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=0.5_r8*(Hz(i,j  ,k)+                                   &
     &                   Hz(i,j-1,k))
            cff=0.5_r8*(Sstopy(i,j  )*Hz(i,j  ,k)+                      &
     &                  Sstopy(i,j-1)*Hz(i,j-1,k))
            rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
            rvsst2d(i,j)=rvsst2d(i,j)+cff*cff1*om_v(i,j)*on_v(i,j)
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3sstm)=cff*om_v(i,j)*on_v(i,j)
#  endif
          END DO
        END DO
# endif
      END DO K_LOOP

      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=0.5_r8*(Dstp(i,j)+Dstp(i-1,j))
# ifdef BOTTOM_STREAMING
          rubst2d(i,j)=rubst2d(i,j)/cff
# endif
# ifdef SURFACE_STREAMING
          russt2d(i,j)=russt2d(i,j)/cff
# endif
!
          IF (j.ge.JstrV) THEN
            cff=0.5_r8*(Dstp(i,j)+Dstp(i,j-1))
# ifdef BOTTOM_STREAMING
            rvbst2d(i,j)=rvbst2d(i,j)/cff
# endif
# ifdef SURFACE_STREAMING
            rvsst2d(i,j)=rvsst2d(i,j)/cff
# endif
          END IF
        END DO
      END DO
!
!  Apply boundary conditions.
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rustr3d)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rvstr3d)
# ifdef BOTTOM_STREAMING
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rubst2d)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rvbst2d)
# endif
# ifdef SURFACE_STREAMING
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  russt2d)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rvsst2d)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rustr3d,  rvstr3d)

#  ifdef BOTTOM_STREAMING
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rubst2d, rvbst2d)
#  endif
#  ifdef SURFACE_STREAMING
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    russt2d, rvsst2d)
#  endif
# endif
      RETURN
      END SUBROUTINE wec_streaming_tile
#endif
      END MODULE wec_streaming_mod
