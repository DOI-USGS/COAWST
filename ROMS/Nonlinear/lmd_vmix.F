#include "cppdefs.h"
      MODULE lmd_vmix_mod
#if defined NONLINEAR && defined LMD_MIXING && defined SOLVE3D
!
!svn $Id: lmd_vmix.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine computes the vertical mixing coefficients for       !
!  momentum and tracers  at the ocean surface boundary layer and       !
!  interior using the Large, McWilliams and Doney  (1994) mixing       !
!  scheme.                                                             !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Large, W.G., J.C. McWilliams, and S.C. Doney, 1994: A Review      !
!      and model with a nonlocal boundary layer parameterization,      !
!      Reviews of Geophysics, 32,363-403.                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: lmd_vmix

      CONTAINS
!
!***********************************************************************
      SUBROUTINE lmd_vmix (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
# ifdef LMD_SKPP
      USE lmd_skpp_mod, ONLY : lmd_skpp
# endif
# ifdef LMD_BKPP
      USE lmd_bkpp_mod, ONLY : lmd_bkpp
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
      CALL wclock_on (ng, iNLM, 18)
# endif
      IF (.not.PerfectRST(ng).or.iic(ng).ne.ntstart(ng)) THEN
        CALL lmd_vmix_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng),                                   &
     &                      GRID(ng) % Hz,                              &
# ifndef RI_SPLINES
     &                      GRID(ng) % z_r,                             &
# endif
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
# ifdef LMD_DDMIX
     &                      OCEAN(ng) % t,                              &
     &                      MIXING(ng) % alfaobeta,                     &
# endif
     &                      MIXING(ng) % bvf,                           &
     &                      MIXING(ng) % Akt,                           &
     &                      MIXING(ng) % Akv)
# ifdef LMD_SKPP
        CALL lmd_skpp (ng, tile)
# endif
# ifdef LMD_BKPP
        CALL lmd_bkpp (ng, tile)
# endif
        CALL lmd_finish (ng, tile)
      END IF
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 18)
# endif
      RETURN
      END SUBROUTINE lmd_vmix
!
!***********************************************************************
      SUBROUTINE lmd_vmix_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp,                                   &
     &                          Hz,                                     &
# ifndef RI_SPLINES
     &                          z_r,                                    &
# endif
     &                          rho, u, v,                              &
# ifdef LMD_DDMIX
     &                          t, alfaobeta,                           &
# endif
     &                          bvf, Akt, Akv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      Integer, intent(in) :: nstp
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#  ifndef RI_SPLINES
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef LMD_DDMIX
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: alfaobeta(LBi:,LBj:,0:)
#  endif
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)

      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#  ifndef RI_SPLINES
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),3)
#  ifdef LMD_DDMIX
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8), parameter :: eps = 1.0E-14_r8

      real(r8) :: cff, lmd_iwm, lmd_iws, nu_sx, shear2
# ifdef LMD_DDMIX
      real(r8) :: Rrho, ddDS, ddDT, nu_dds, nu_ddt
# endif

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: Rig

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dU
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dV

# include "set_bounds.h"

# ifdef LMD_RIMIX
!
!-----------------------------------------------------------------------
! Compute gradient Richardson number.
!-----------------------------------------------------------------------
!
!  Compute gradient Richardson number at horizontal RHO-points and
!  vertical W-points.  If zero or very small velocity shear, bound
!  computation by a large negative value.
!
#  ifdef RI_SPLINES
      DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
        DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
          FC(i,0)=0.0_r8
          dR(i,0)=0.0_r8
          dU(i,0)=0.0_r8
          dV(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                             &
     &                  Hz(i,j,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*Hz(i,j,k+1)
            dR(i,k)=cff*(6.0_r8*(rho(i,j,k+1)-rho(i,j,k))-              &
     &                   Hz(i,j,k)*dR(i,k-1))
            dU(i,k)=cff*(3.0_r8*(u(i  ,j,k+1,nstp)-u(i  ,j,k,nstp)+     &
     &                           u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))-    &
     &                   Hz(i,j,k)*dU(i,k-1))
            dV(i,k)=cff*(3.0_r8*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+     &
     &                           v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp))-    &
     &                   Hz(i,j,k)*dV(i,k-1))
          END DO
        END DO
        DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
          dR(i,N(ng))=0.0_r8
          dU(i,N(ng))=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            dR(i,k)=dR(i,k)-FC(i,k)*dR(i,k+1)
            dU(i,k)=dU(i,k)-FC(i,k)*dU(i,k+1)
            dV(i,k)=dV(i,k)-FC(i,k)*dV(i,k+1)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            shear2=dU(i,k)*dU(i,k)+dV(i,k)*dV(i,k)
            Rig(i,j,k)=bvf(i,j,k)/max(shear2,eps)
!!          Rig(i,j,k)=-gorho0*dR(i,k)/(shear2+eps)
          END DO
        END DO
      END DO
#  else
      DO k=1,N(ng)-1
        DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            shear2=(cff*(u(i  ,j,k+1,nstp)-u(i  ,j,k,nstp)+             &
     &                   u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp)))**2+        &
     &             (cff*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+             &
     &                   v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))**2
            Rig(i,j,k)=bvf(i,j,k)/max(shear2,eps)
          END DO
        END DO
      END DO
#  endif
#  ifdef RI_HORAVG
      DO k=1,N(ng)-1
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
            Rig(Istr-1,j,k)=Rig(Istr,j,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
            Rig(Iend+1,j,k)=Rig(Iend,j,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            Rig(i,Jstr-1,k)=Rig(i,Jstr,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            Rig(i,Jend+1,k)=Rig(i,Jend,k)
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          Rig(Istr-1,Jstr-1,k)=Rig(Istr,Jstr,k)
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          Rig(Istr-1,Jend+1,k)=Rig(Istr,Jend,k)
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          Rig(Iend+1,Jstr-1,k)=Rig(Iend,Jstr,k)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          Rig(Iend+1,Jend+1,k)=Rig(Iend,Jend,k)
        END IF
!
!  Smooth gradient Richardson number horizontally.  Use Rig(:,:,0)
!  as scratch utility array.
!
        DO j=Jstr-1,Jend
          DO i=Istr-1,Iend
            Rig(i,j,0)=0.25_r8*(Rig(i,j  ,k)+Rig(i+1,j  ,k)+            &
     &                          Rig(i,j+1,k)+Rig(i+1,j+1,k))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Rig(i,j,k)=0.25_r8*(Rig(i,j  ,0)+Rig(i-1,j  ,0)+            &
     &                          Rig(i,j-1,0)+Rig(i-1,j-1,0))
          END DO
        END DO
      END DO
#  endif
#  ifdef RI_VERAVG
!
!  Smooth gradient Richardson number vertically at the interior points.
!
      DO k=N(ng)-2,2,-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Rig(i,j,k)=0.25_r8*Rig(i,j,k-1)+                            &
     &                 0.50_r8*Rig(i,j,k  )+                            &
     &                 0.25_r8*Rig(i,j,k+1)
          END DO
        END DO
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute "interior" viscosities and diffusivities everywhere as
!  the superposition of three processes: local Richardson number
!  instability due to resolved vertical shear, internal wave
!  breaking, and double diffusion.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
!
!  Compute interior diffusivity due to shear instability mixing.
!
# ifdef LMD_RIMIX
            cff=MIN(1.0_r8,MAX(0.0_r8,Rig(i,j,k))/lmd_Ri0)
            nu_sx=1.0_r8-cff*cff
            nu_sx=nu_sx*nu_sx*nu_sx
!
!  The shear mixing should be also a function of the actual magnitude
!  of the shear, see Polzin (1996, JPO, 1409-1425).
!
            shear2=bvf(i,j,k)/(Rig(i,j,k)+eps)
            cff=shear2*shear2/(shear2*shear2+16.0E-10_r8)
            nu_sx=cff*nu_sx
# else
            nu_sx=0.0_r8
# endif
!
!  Compute interior diffusivity due to wave breaking (Gargett and
!  Holloway.
!
            cff=1.0_r8/SQRT(MAX(bvf(i,j,k),1.0E-7_r8))
# ifdef ISOMIP
            lmd_iwm=1.0E-3_r8
            lmd_iws=5.0E-5_r8
            if (bvf(i,j,k) .lt. 0.0_r8) then
              lmd_iwm=lmd_iwm*100.0_r8
              lmd_iws=lmd_iws*100.0_r8
            endif
# else
            lmd_iwm=1.0E-6_r8*cff
            lmd_iws=1.0E-7_r8*cff
!           lmd_iwm=lmd_nuwm
!           lmd_iws=lmd_nuws
# endif
!
! Sum contributions due to internal wave breaking, shear instability
! and convective diffusivity due to shear instability.
!
# ifdef ISOMIP
            Akv(i,j,k)=lmd_iwm
            Akt(i,j,k,itemp)=lmd_iws
# else
            Akv(i,j,k)=lmd_iwm+lmd_nu0m*nu_sx
            Akt(i,j,k,itemp)=lmd_iws+lmd_nu0s*nu_sx
# endif
# ifdef SALINITY
            Akt(i,j,k,isalt)=Akt(i,j,k,itemp)
# endif
          END DO
        END DO
# ifdef LMD_DDMIX
!
!-----------------------------------------------------------------------
!  Compute double-diffusive mixing.  It can occur when vertical
!  gradient of density is stable but the vertical gradient of
!  salinity (salt figering) or temperature (diffusive convection)
!  is unstable.
!-----------------------------------------------------------------------
!
!  Compute double-diffusive density ratio, Rrho.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            ddDT=t(i,j,k+1,nstp,itemp)-t(i,j,k,nstp,itemp)
            ddDS=t(i,j,k+1,nstp,isalt)-t(i,j,k,nstp,isalt)
            ddDS=SIGN(1.0_r8,ddDS)*MAX(ABS(ddDS),1.0E-14_r8)
            Rrho=alfaobeta(i,j,k)*ddDT/ddDS
!
!  Salt fingering case.
!
            IF ((Rrho.gt.1.0_r8).and.(ddDS.gt.0.0_r8)) THEN
!
!  Compute interior diffusivity for double diffusive mixing of
!  salinity.  Upper bound "Rrho" by "Rrho0"; (lmd_Rrho0=1.9,
!  lmd_nuf=0.001).
!
              Rrho=MIN(Rrho,lmd_Rrho0)
              nu_dds=1.0_r8-((Rrho-1.0_r8)/(lmd_Rrho0-1.0_r8))**2
              nu_dds=lmd_nuf*nu_dds*nu_dds*nu_dds
!
!  Compute interior diffusivity for double diffusive mixing
!  of temperature (lmd_fdd=0.7).
!
              nu_ddt=lmd_fdd*nu_dds
!
!  Diffusive convection case.
!
            ELSE IF ((0.0_r8.lt.Rrho).and.(Rrho.lt.1.0_r8).and.         &
     &              (ddDS.lt.0.0_r8)) THEN
!
!  Compute interior diffusivity for double diffusive mixing of
!  temperature (Marmorino and Caldwell, 1976); (lmd_nu=1.5e-6,
!  lmd_tdd1=0.909, lmd_tdd2=4.6, lmd_tdd3=0.54).
!
              nu_ddt=lmd_nu*lmd_tdd1*                                   &
     &               EXP(lmd_tdd2*                                      &
     &                   EXP(-lmd_tdd3*((1.0_r8/Rrho)-1.0_r8)))
!
!  Compute interior diffusivity for double diffusive mixing
!  of salinity (lmd_sdd1=0.15, lmd_sdd2=1.85, lmd_sdd3=0.85).
!
              IF (Rrho.lt.0.5_r8) THEN
                nu_dds=nu_ddt*lmd_sdd1*Rrho
              ELSE
                nu_dds=nu_ddt*(lmd_sdd2*Rrho-lmd_sdd3)
              END IF
            ELSE
              nu_ddt=0.0_r8
              nu_dds=0.0_r8
            END IF
!
!  Add double diffusion contribution to temperature and salinity
!  mixing coefficients.
!
            Akt(i,j,k,itemp)=Akt(i,j,k,itemp)+nu_ddt
#  ifdef SALINITY
            Akt(i,j,k,isalt)=Akt(i,j,k,isalt)+nu_dds
#  endif
          END DO
        END DO
# endif
      END DO

      RETURN
      END SUBROUTINE lmd_vmix_tile
!
!***********************************************************************
      SUBROUTINE lmd_finish (ng, tile)
!***********************************************************************
!
      USE mod_param
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
      CALL lmd_finish_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      nstp(ng),                                   &
     &                      GRID(ng) % Hz,                              &
# ifdef M2TIDE_DIFF
     &                      GRID(ng) % z_w,                             &
     &                      MIXING(ng) % m2_amp,                        &
# endif
     &                      MIXING(ng) % bvf,                           &
     &                      MIXING(ng) % Akt,                           &
     &                      MIXING(ng) % Akv)
      RETURN
      END SUBROUTINE lmd_finish
!
!***********************************************************************
      SUBROUTINE lmd_finish_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            nstp,                                 &
     &                            Hz,                                   &
# ifdef M2TIDE_DIFF
     &                            z_w,                                  &
     &                            m2_amp,                               &
# endif
     &                            bvf, Akt, Akv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_3d_mod, ONLY : bc_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      Integer, intent(in) :: nstp
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
# ifdef M2TIDE_DIFF
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: m2_amp(LBi:,LBj:)
# endif
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)

      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
# ifdef M2TIDE_DIFF
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: m2_amp(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8), parameter :: eps = 1.0E-14_r8

      real(r8) :: cff, nu_sx, nu_sxc, shear2
# ifdef M2TIDE_DIFF
      real(r8), parameter :: lmd_k0 = 5.0e-3_r8  ! m2 s-1
      real(r8), parameter :: lmd_sigma = 3.0_r8
      real(r8), parameter :: lmd_p = 0.25_r8
      real(r8), parameter :: lmd_drag = 2.4e-3_r8
      real(r8), parameter :: eps = 1.e-20_r8
      real(r8) :: lmd_kz, lmd_ri
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute "interior" viscosities and diffusivities everywhere as
!  the superposition of three processes: local Richardson number
!  instability due to resolved vertical shear, internal wave
!  breaking, and double diffusion.
# ifdef M2TIDE_DIFF
!  Add in contribution of simulated tidal mixing.
!   Lee, H.-C. et al., 2006 Ocean Modelling, 11, 464-477.
# endif
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
!
!  Compute interior convective diffusivity due to static instability
!  mixing.
!
# ifdef LMD_CONVEC
            cff=MAX(bvf(i,j,k),lmd_bvfcon)
            cff=MIN(1.0_r8,(lmd_bvfcon-cff)/lmd_bvfcon)
            nu_sxc=1.0_r8-cff*cff
            nu_sxc=nu_sxc*nu_sxc*nu_sxc
# else
            nu_sxc=0.0_r8
# endif
!
! Sum contributions due to internal wave breaking, shear instability
! and convective diffusivity due to shear instability.
!
# ifndef ISOMIP
            Akv(i,j,k)=Akv(i,j,k)+lmd_nu0c*nu_sxc
            Akt(i,j,k,itemp)=Akt(i,j,k,itemp)+lmd_nu0c*nu_sxc
#  ifdef SALINITY
            Akt(i,j,k,isalt)=Akt(i,j,k,isalt)+lmd_nu0c*nu_sxc
#  endif
# endif
          END DO
        END DO
      END DO
!
# ifdef M2TIDE_DIFF
!  Add in contribution of simulated tidal mixing.
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            shear2 = 0.5 * lmd_drag *                                   &
     &               (m2_amp(i,j)/(vonKar*(z_w(i,j,0)-z_w(i,j,k)) ))**2
            lmd_ri = max(0.0,bvf(i,j,k)) / (shear2 + eps)
            lmd_kz = lmd_k0 * (1 + lmd_sigma*lmd_ri)**(-lmd_p)
            Akt(i,j,k,itemp)=Akt(i,j,k,itemp)+lmd_kz
#  ifdef SALINITY
            Akt(i,j,k,isalt)=Akt(i,j,k,isalt)+lmd_kz
#  endif
          END DO
        END DO
      END DO
# endif
!
!  Apply boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  Akv)
      DO itrc=1,NAT
        CALL bc_w3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    Akt(:,:,:,itrc))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akv)
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng), 1, NAT,         &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akt)
# endif
      RETURN
      END SUBROUTINE lmd_finish_tile
#endif
      END MODULE lmd_vmix_mod
