#include "cppdefs.h"
      MODULE gls_corstep_mod
#if defined NONLINEAR && defined GLS_MIXING && defined SOLVE3D
!
!svn $Id: gls_corstep.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine perfoms the corrector step for turbulent kinetic       !
!  energy and generic length scale prognostic variables, tke and       !
!  gls.                                                                !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Umlauf, L. and H. Burchard, 2003:  A generic length-scale           !
!    Equation for geophysical turbulence models, J. Marine Res.,       !
!    61, 235-265.                                                      !
!                                                                      !
!  Warner, J. C., C. R. Sherwood, H. G. Arango and R. P. Signell,      !
!    2005:  Performance of four turbulence closure models              !
!    implemented using a generic length scale method, Ocean            !
!    Modelling, 8, 81-113.                                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: gls_corstep

      CONTAINS
!
!***********************************************************************
      SUBROUTINE gls_corstep (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
# if defined VEGETATION && defined VEG_TURB
      USE mod_vegarr
      USE vegetation_turb_mod, ONLY : vegetation_turb_cal
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
# if defined VEGETATION && defined VEG_TURB
! 
!-----------------------------------------------------------------------
!  Add the effect of vegetation on the turbulence model.
!-----------------------------------------------------------------------
!
      CALL vegetation_turb_cal (ng, tile)
# endif

# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 19)
# endif
      CALL gls_corstep_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp(ng), nnew(ng),                        &
# ifdef MASKING
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
     &                       GRID(ng) % Huon,                           &
     &                       GRID(ng) % Hvom,                           &
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % z_w,                            &
     &                       GRID(ng) % ZoBot,                          &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % W,                             &
# ifdef WEC_VF
     &                       OCEAN(ng) % W_stokes,                      &
# endif
# if defined VEGETATION && defined VEG_TURB 
     &                       VEG(ng) % tke_veg,                         &
     &                       VEG(ng) % gls_veg,                         &
# endif
     &                       FORCES(ng) % bustr,                        &
     &                       FORCES(ng) % bvstr,                        &
     &                       FORCES(ng) % sustr,                        &
     &                       FORCES(ng) % svstr,                        &
# ifdef ZOS_HSIG
     &                       FORCES(ng) % Hwave,                        &
# endif
# ifdef TKE_WAVEDISS
     &                       FORCES(ng) % Dissip_break,                 &
     &                       FORCES(ng) % Dissip_wcap,                  &
# endif
     &                       MIXING(ng) % Akt,                          &
     &                       MIXING(ng) % Akv,                          &
     &                       MIXING(ng) % bvf,                          &
     &                       MIXING(ng) % Akk,                          &
     &                       MIXING(ng) % Akp,                          &
     &                       MIXING(ng) % Lscale,                       &
     &                       MIXING(ng) % gls,                          &
     &                       MIXING(ng) % tke)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 19)
# endif
      RETURN
      END SUBROUTINE gls_corstep
!
!***********************************************************************
      SUBROUTINE gls_corstep_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew,                          &
# ifdef MASKING
     &                             umask, vmask,                        &
# endif
     &                             Huon, Hvom, Hz, pm, pn, z_r, z_w,    &
     &                             ZoBot,                               &
     &                             u, v, W,                             &
# ifdef WEC_VF
     &                             W_stokes,                            &
# endif
# if defined VEGETATION && defined VEG_TURB
     &                             tke_veg, gls_veg,                    &
# endif
     &                             bustr, bvstr, sustr, svstr,          &
# ifdef ZOS_HSIG
     &                             Hwave,                               &
# endif
# ifdef TKE_WAVEDISS
     &                             Dissip_break, Dissip_wcap,           &
# endif
     &                             Akt, Akv, bvf,                       &
     &                             Akk, Akp, Lscale, gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
      USE tkebc_mod, ONLY : tkebc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
# ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:,LBj:,0:)
# endif
# if defined VEGETATION && defined VEG_TURB
      real(r8), intent(in) :: tke_veg(LBi:,LBj:,0:)
      real(r8), intent(in) :: gls_veg(LBi:,LBj:,0:)
# endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  ifdef ZOS_HSIG
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
#  endif
#  ifdef TKE_WAVEDISS
      real(r8), intent(in) :: Dissip_break(LBi:,LBj:)
      real(r8), intent(in) :: Dissip_wcap(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)

      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akk(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akp(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Lscale(LBi:,LBj:,0:)
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: ZoBot(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
# ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:UBi,LBj:UBj,0:N(ng))
# endif
# if defined VEGETATION && defined VEG_TURB
      real(r8), intent(in) :: tke_veg(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: gls_veg(LBi:UBi,LBj:UBj,0:N(ng))
# endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  if defined ZOS_HSIG
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
#  endif
#  ifdef TKE_WAVEDISS
      real(r8), intent(in) :: Dissip_break(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dissip_wcap(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: Akk(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: Akp(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: Lscale(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: gls(LBi:UBi,LBj:UBj,0:N(ng),3)
      real(r8), intent(inout) :: tke(LBi:UBi,LBj:UBj,0:N(ng),3)
# endif
!
!  Local variable declarations.
!
      logical :: Lmy25
      integer :: i, itrc, j, k

#ifdef AKLIMIT
      real(r8), parameter :: Ak_max = 10.0_r8
#endif
      real(r8), parameter :: Gadv = 1.0_r8/3.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8) :: Zos_min

      real(r8) :: Gh, Gm, Kprod, Ls_unlmt, Ls_lmt, Pprod, Sh, Sm
      real(r8) :: cff, cff1, cff2, cff3
      real(r8) :: cmu_fac1, cmu_fac2, cmu_fac3, cmu_fac4
      real(r8) :: gls_c3, gls_exp1, gls_fac1, gls_fac2, gls_fac3
      real(r8) :: gls_fac4, gls_fac5, gls_fac6, ql, sqrt2, strat2
      real(r8) :: tke_exp1, tke_exp2, tke_exp3, tke_exp4, wall_fac
      real(r8) :: gls_d, gls_sigp_cb, ogls_sigp, sig_eff
      real(r8) :: L_sft
# if defined CRAIG_BANNER || defined TKE_WAVEDISS
      real(r8) :: cb_wallE
# endif

      real(r8), dimension(IminS:ImaxS) :: tke_fluxt
      real(r8), dimension(IminS:ImaxS) :: tke_fluxb
      real(r8), dimension(IminS:ImaxS) :: gls_fluxt
      real(r8), dimension(IminS:ImaxS) :: gls_fluxb
      real(r8), dimension(IminS:ImaxS) :: Zos_eff

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BCK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BCP
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCP
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dU
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dV

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shear2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: buoy2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zob_min
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curvK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curvP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradP

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute several constants.
!-----------------------------------------------------------------------
!
      Zos_min=MAX(Zos(ng),0.0001_r8)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Zob_min(i,j)=MAX(ZoBot(i,j),0.0001_r8)
        END DO
      END DO
!
      Lmy25=.FALSE.
      IF ((gls_p(ng).eq.0.0_r8).and.                                    &
     &    (gls_n(ng).eq.1.0_r8).and.                                    &
     &    (gls_m(ng).eq.1.0_r8)) THEN
        Lmy25=.TRUE.
      END IF

# if defined CRAIG_BANNER || defined TKE_WAVEDISS
      IF (Lmy25) THEN
!!      cb_wallE=1.0_r8+gls_E2
        cb_wallE=1.25_r8
      ELSE
        cb_wallE=1.0_r8
      END IF
      L_sft=vonKar
      cff1=sqrt(1.5_r8*gls_sigk(ng))*gls_cmu0(ng)/L_sft
      gls_sigp_cb=L_sft**2/(gls_cmu0(ng)**2*gls_c2(ng)*cb_wallE)*       &
     &            (gls_n(ng)**2-                                        &
     &             cff1*gls_n(ng)/3.0_r8*(4.0_r8*gls_m(ng)+1.0_r8)+     &
     &             cff1**2*gls_m(ng)/9.0_r8*(2.0_r8+4.0_r8*gls_m(ng)))
# else
      L_sft=vonKar
      gls_sigp_cb=gls_sigp(ng)
# endif
      ogls_sigp=1.0_r8/gls_sigp_cb
!
      sqrt2=SQRT(2.0_r8)
      cmu_fac1=gls_cmu0(ng)**(-gls_p(ng)/gls_n(ng))
      cmu_fac2=gls_cmu0(ng)**(3.0_r8+gls_p(ng)/gls_n(ng))
      cmu_fac3=1.0_r8/gls_cmu0(ng)**2.0_r8
      cmu_fac4=(1.5_r8*gls_sigk(ng))**(1.0_r8/3.0_r8)/                  &
     &         (gls_cmu0(ng)**(4.0_r8/3.0_r8))
      gls_fac1=gls_n(ng)*gls_cmu0(ng)**(gls_p(ng)+1.0_r8)
      gls_fac2=gls_cmu0(ng)**(gls_p(ng))*gls_n(ng)*                     &
     &         vonKar**(gls_n(ng))
      gls_fac3=gls_cmu0(ng)**(gls_p(ng))*gls_n(ng)
      gls_fac4=gls_cmu0(ng)**(gls_p(ng))
      gls_fac5=0.56_r8**(0.5_r8*gls_n(ng))*gls_cmu0(ng)**gls_p(ng)
      gls_fac6=8.0_r8/gls_cmu0(ng)**6.0_r8
!
      gls_exp1=1.0_r8/gls_n(ng)
      tke_exp1=gls_m(ng)/gls_n(ng)
      tke_exp2=0.5_r8+gls_m(ng)/gls_n(ng)
      tke_exp3=0.5_r8+gls_m(ng)
      tke_exp4=gls_m(ng)+0.5_r8*gls_n(ng)
!
!-----------------------------------------------------------------------
!  Compute vertical velocity shear at W-points.
!-----------------------------------------------------------------------
!
# ifdef RI_SPLINES
      DO j=Jstrm1,Jendp1
        DO i=Istrm1,Iendp1
          CF(i,0)=0.0_r8
          dU(i,0)=0.0_r8
          dV(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=Istrm1,Iendp1
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                             &
     &                  Hz(i,j,k)*(2.0_r8-CF(i,k-1)))
            CF(i,k)=cff*Hz(i,j,k+1)
            dU(i,k)=cff*(3.0_r8*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp)+     &
     &                           u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))-    &
     &                   Hz(i,j,k)*dU(i,k-1))
            dV(i,k)=cff*(3.0_r8*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+     &
     &                           v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp))-    &
     &                   Hz(i,j,k)*dV(i,k-1))
          END DO
        END DO
        DO i=Istrm1,Iendp1
          dU(i,N(ng))=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=Istrm1,Iendp1
            dU(i,k)=dU(i,k)-CF(i,k)*dU(i,k+1)
            dV(i,k)=dV(i,k)-CF(i,k)*dV(i,k+1)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istrm1,Iendp1
            shear2(i,j,k)=dU(i,k)*dU(i,k)+dV(i,k)*dV(i,k)
          END DO
        END DO
      END DO
# else
      DO k=1,N(ng)-1
        DO j=Jstrm1,Jendp1
          DO i=Istrm1,Iendp1
            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            shear2(i,j,k)=(cff*(u(i  ,j,k+1,nstp)-u(i  ,j,k,nstp)+      &
     &                          u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp)))**2+ &
     &                    (cff*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+      &
     &                          v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))**2
          END DO
        END DO
      END DO
# endif
!
! Load Brunt-Vaisala frequency.
!
      DO k=1,N(ng)-1
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            buoy2(i,j,k)=bvf(i,j,k)
          END DO
        END DO
      END DO
# ifdef N2S2_HORAVG
!
!-----------------------------------------------------------------------
!  Smooth horizontally buoyancy and shear.  Use buoy2(:,:,0) and
!  shear2(:,:,0) as scratch utility array.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
            shear2(Istr-1,j,k)=shear2(Istr,j,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
            shear2(Iend+1,j,k)=shear2(Iend,j,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            shear2(i,Jstr-1,k)=shear2(i,Jstr,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            shear2(i,Jend+1,k)=shear2(i,Jend,k)
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          shear2(Istr-1,Jstr-1,k)=shear2(Istr,Jstr,k)
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          shear2(Istr-1,Jend+1,k)=shear2(Istr,Jend,k)
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          shear2(Iend+1,Jstr-1,k)=shear2(Iend,Jstr,k)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          shear2(Iend+1,Jend+1,k)=shear2(Iend,Jend,k)
        END IF
!
!  Average horizontally.
!
        DO j=Jstr-1,Jend
          DO i=Istr-1,Iend
            buoy2(i,j,0)=0.25_r8*(buoy2(i,j  ,k)+buoy2(i+1,j  ,k)+      &
     &                            buoy2(i,j+1,k)+buoy2(i+1,j+1,k))
            shear2(i,j,0)=0.25_r8*(shear2(i,j  ,k)+shear2(i+1,j  ,k)+   &
     &                             shear2(i,j+1,k)+shear2(i+1,j+1,k))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            buoy2(i,j,k)=0.25_r8*(buoy2(i,j  ,0)+buoy2(i-1,j  ,0)+      &
     &                            buoy2(i,j-1,0)+buoy2(i-1,j-1,0))
            shear2(i,j,k)=0.25_r8*(shear2(i,j  ,0)+shear2(i-1,j  ,0)+   &
     &                             shear2(i,j-1,0)+shear2(i-1,j-1,0))
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Time-step advective terms.
!-----------------------------------------------------------------------
!
!  At entry, it is assumed that the turbulent kinetic energy fields
!  "tke" and "gls", at time level "nnew", are set to its values at
!  time level "nstp" multiplied by the grid box thicknesses Hz
!  (from old time step and at W-points).
!
      DO k=1,N(ng)-1
# ifdef K_C2ADVECTION
!
!  Second-order, centered differences advection.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.25_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            FXK(i,j)=cff*(tke(i,j,k,3)+tke(i-1,j,k,3))
            FXP(i,j)=cff*(gls(i,j,k,3)+gls(i-1,j,k,3))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.25_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FEK(i,j)=cff*(tke(i,j,k,3)+tke(i,j-1,k,3))
            FEP(i,j)=cff*(gls(i,j,k,3)+gls(i,j-1,k,3))
          END DO
        END DO
# else
        DO j=Jstr,Jend
          DO i=Istrm1,Iendp2
            gradK(i,j)=(tke(i,j,k,3)-tke(i-1,j,k,3))
#  ifdef MASKING
            gradK(i,j)=gradK(i,j)*umask(i,j)
#  endif
            gradP(i,j)=(gls(i,j,k,3)-gls(i-1,j,k,3))
#  ifdef MASKING
            gradP(i,j)=gradP(i,j)*umask(i,j)
#  endif
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              gradK(Istr-1,j)=gradK(Istr,j)
              gradP(Istr-1,j)=gradP(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              gradK(Iend+2,j)=gradK(Iend+1,j)
              gradP(Iend+2,j)=gradP(Iend+1,j)
            END DO
          END IF
        END IF
#  ifdef K_C4ADVECTION
!
!  Fourth-order, centered differences advection.
!
        cff1=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            FXK(i,j)=cff*0.5_r8*(tke(i-1,j,k,3)+tke(i,j,k,3)-           &
     &                           cff1*(gradK(i+1,j)-gradK(i-1,j)))
            FXP(i,j)=cff*0.5_r8*(gls(i-1,j,k,3)+gls(i,j,k,3)-           &
     &                           cff1*(gradP(i+1,j)-gradP(i-1,j)))
          END DO
        END DO
#  else
!
!  Third-order, upstream bias advection with velocity dependent
!  hyperdiffusion.
!
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
            curvK(i,j)=gradK(i+1,j)-gradK(i,j)
            curvP(i,j)=gradP(i+1,j)-gradP(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            IF (cff.gt.0.0_r8) THEN
              cff1=curvK(i-1,j)
              cff2=curvP(i-1,j)
            ELSE
              cff1=curvK(i,j)
              cff2=curvP(i,j)
            END IF
            FXK(i,j)=cff*0.5_r8*(tke(i-1,j,k,3)+tke(i,j,k,3)-           &
     &                           Gadv*cff1)
            FXP(i,j)=cff*0.5_r8*(gls(i-1,j,k,3)+gls(i,j,k,3)-           &
     &                           Gadv*cff2)
          END DO
        END DO
#  endif
        DO j=Jstrm1,Jendp2
          DO i=Istr,Iend
            gradK(i,j)=(tke(i,j,k,3)-tke(i,j-1,k,3))
#  ifdef MASKING
            gradK(i,j)=gradK(i,j)*vmask(i,j)
#  endif
            gradP(i,j)=(gls(i,j,k,3)-gls(i,j-1,k,3))
#  ifdef MASKING
            gradP(i,j)=gradP(i,j)*vmask(i,j)
#  endif
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              gradK(i,Jstr-1)=gradK(i,Jstr)
              gradP(i,Jstr-1)=gradP(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              gradK(i,Jend+2)=gradK(i,Jend+1)
              gradP(i,Jend+2)=gradP(i,Jend+1)
            END DO
          END IF
        END IF
#  ifdef K_C4ADVECTION
        cff1=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.5_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FEK(i,j)=cff*0.5_r8*(tke(i,j-1,k,3)+tke(i,j,k,3)-           &
     &                           cff1*(gradK(i,j+1)-gradK(i,j-1)))
            FEP(i,j)=cff*0.5_r8*(gls(i,j-1,k,3)+gls(i,j,k,3)-           &
     &                           cff1*(gradP(i,j+1)-gradP(i,j-1)))
          END DO
        END DO
#  else
        DO j=Jstr-1,Jend+1
          DO i=Istr,Iend
            curvK(i,j)=gradK(i,j+1)-gradK(i,j)
            curvP(i,j)=gradP(i,j+1)-gradP(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.5_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            IF (cff.gt.0.0_r8) THEN
              cff1=curvK(i,j-1)
              cff2=curvP(i,j-1)
            ELSE
              cff1=curvK(i,j)
              cff2=curvP(i,j)
            END IF
            FEK(i,j)=cff*0.5_r8*(tke(i,j-1,k,3)+tke(i,j,k,3)-           &
     &                           Gadv*cff1)
            FEP(i,j)=cff*0.5_r8*(gls(i,j-1,k,3)+gls(i,j,k,3)-           &
     &                           Gadv*cff2)
          END DO
        END DO
#  endif
# endif
!
!  Time-step horizontal advection.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-                            &
     &                      cff*(FXK(i+1,j)-FXK(i,j)+                   &
     &                           FEK(i,j+1)-FEK(i,j))
            tke(i,j,k,nnew)=MAX(tke(i,j,k,nnew),gls_Kmin(ng))
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-                            &
     &                      cff*(FXP(i+1,j)-FXP(i,j)+                   &
     &                           FEP(i,j+1)-FEP(i,j))
            gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),gls_Pmin(ng))
          END DO
        END DO
      END DO
!
! Compute vertical advection.
!
      DO j=Jstr,Jend
# ifdef K_C2ADVECTION
        DO k=1,N(ng)
          DO i=Istr,Iend
            cff=0.25_r8*(W(i,j,k)+W(i,j,k-1))
#  ifdef WEC_VF
            cff=cff+0.25_r8*(W_stokes(i,j,k)+W_stokes(i,j,k-1))
#  endif
            FCK(i,k)=cff*(tke(i,j,k,3)+tke(i,j,k-1,3))
            FCP(i,k)=cff*(gls(i,j,k,3)+gls(i,j,k-1,3))
          END DO
        END DO
# else
        cff1=7.0_r8/12.0_r8
        cff2=1.0_r8/12.0_r8
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            cff=0.5*(W(i,j,k)+W(i,j,k-1))
#  ifdef WEC_VF
            cff=cff+0.5_r8*(W_stokes(i,j,k)+W_stokes(i,j,k-1))
#  endif
            FCK(i,k)=cff*(cff1*(tke(i,j,k-1,3)+                         &
     &                          tke(i,j,k  ,3))-                        &
     &                    cff2*(tke(i,j,k-2,3)+                         &
     &                          tke(i,j,k+1,3)))
            FCP(i,k)=cff*(cff1*(gls(i,j,k-1,3)+                         &
     &                          gls(i,j,k  ,3))-                        &
     &                    cff2*(gls(i,j,k-2,3)+                         &
     &                          gls(i,j,k+1,3)))
          END DO
        END DO
        cff1=1.0_r8/3.0_r8
        cff2=5.0_r8/6.0_r8
        cff3=1.0_r8/6.0_r8
         DO i=Istr,Iend
          cff=0.5_r8*(W(i,j,0)+W(i,j,1))
#  ifdef WEC_VF
          cff=cff+0.5_r8*(W_stokes(i,j,0)+W_stokes(i,j,1))
#  endif
          FCK(i,1)=cff*(cff1*tke(i,j,0,3)+                              &
     &                  cff2*tke(i,j,1,3)-                              &
     &                  cff3*tke(i,j,2,3))
          FCP(i,1)=cff*(cff1*gls(i,j,0,3)+                              &
     &                  cff2*gls(i,j,1,3)-                              &
     &                  cff3*gls(i,j,2,3))
          cff=0.5_r8*(W(i,j,N(ng))+W(i,j,N(ng)-1))
#  ifdef WEC_VF
          cff=cff+0.5_r8*(W_stokes(i,j,N(ng))+W_stokes(i,j,N(ng)-1))
#  endif
          FCK(i,N(ng))=cff*(cff1*tke(i,j,N(ng)  ,3)+                    &
     &                      cff2*tke(i,j,N(ng)-1,3)-                    &
     &                      cff3*tke(i,j,N(ng)-2,3))
          FCP(i,N(ng))=cff*(cff1*gls(i,j,N(ng)  ,3)+                    &
     &                      cff2*gls(i,j,N(ng)-1,3)-                    &
     &                      cff3*gls(i,j,N(ng)-2,3))
        END DO
# endif
!
!  Time-step vertical advection term.
!
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-                            &
     &                      cff*(FCK(i,k+1)-FCK(i,k))
            tke(i,j,k,nnew)=MAX(tke(i,j,k,nnew),gls_Kmin(ng))
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-                            &
     &                      cff*(FCP(i,k+1)-FCP(i,k))
            gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),gls_Pmin(ng))
          END DO
        END DO
!
!----------------------------------------------------------------------
!  Compute vertical mixing, turbulent production and turbulent
!  dissipation terms.
!----------------------------------------------------------------------
!
!  Set term for vertical mixing of turbulent fields.
!
        cff=-0.5_r8*dt(ng)
        DO i=Istr,Iend
          DO k=2,N(ng)-1
            FCK(i,k)=cff*(Akk(i,j,k)+Akk(i,j,k-1))/Hz(i,j,k)
            FCP(i,k)=cff*(Akp(i,j,k)+Akp(i,j,k-1))/Hz(i,j,k)
            CF(i,k)=0.0_r8
          END DO
          FCP(i,1)=0.0_r8
          FCP(i,N(ng))=0.0_r8
          FCK(i,1)=0.0_r8
          FCK(i,N(ng))=0.0_r8
        END DO
!
!  Compute production and dissipation terms.
!
        DO i=Istr,Iend
          DO k=1,N(ng)-1
!
!  Compute shear and bouyant production of turbulent energy (m3/s3)
!  at W-points (ignore small negative values of buoyancy).
!
            strat2=buoy2(i,j,k)
            IF (strat2.gt.0.0_r8) THEN
              gls_c3=gls_c3m(ng)
            ELSE
              gls_c3=gls_c3p(ng)
            END IF
            Kprod=shear2(i,j,k)*(Akv(i,j,k)-Akv_bak(ng))-               &
     &            strat2*(Akt(i,j,k,itemp)-Akt_bak(itemp,ng))
            Pprod=gls_c1(ng)*shear2(i,j,k)*(Akv(i,j,k)-Akv_bak(ng))-    &
     &            gls_c3*strat2*(Akt(i,j,k,itemp)-Akt_bak(itemp,ng))
!
!  If negative production terms, then add buoyancy to dissipation terms
!  (BCK and BCP) below, using "cff1" and "cff2" as the on/off switch.
!
            cff1=1.0_r8
            IF (Kprod.lt.0.0_r8) THEN
              Kprod=Kprod+strat2*(Akt(i,j,k,itemp)-Akt_bak(itemp,ng))
              cff1=0.0_r8
            END IF
            cff2=1.0_r8
            IF (Pprod.lt.0.0_r8) THEN
              Pprod=Pprod+gls_c3*strat2*(Akt(i,j,k,itemp)-              &
     &                                   Akt_bak(itemp,ng))
              cff2=0.0_r8
            END IF
!
!  Time-step shear and buoyancy production terms.
!
            cff=0.5_r8*(Hz(i,j,k)+Hz(i,j,k+1))
            tke(i,j,k,nnew)=tke(i,j,k,nnew)+                            &
     &                      dt(ng)*cff*Kprod
            gls(i,j,k,nnew)=gls(i,j,k,nnew)+                            &
     &                      dt(ng)*cff*Pprod*gls(i,j,k,nstp)/           &
     &                      MAX(tke(i,j,k,nstp),gls_Kmin(ng))

# if defined VEGETATION && defined VEG_TURB
!
!  Add the effect of vegetation on tke and gls 
!          
            tke(i,j,k,nnew)=tke(i,j,k,nnew)+dt(ng)*tke_veg(i,j,k)
            gls(i,j,k,nnew)=gls(i,j,k,nnew)+dt(ng)*gls_veg(i,j,k)
# endif 
!
!  Compute dissipation of turbulent energy (m3/s3).
!
            wall_fac=1.0_r8
            IF (Lmy25) THEN
!
!  Parabolic wall function,  L = ds db / (ds + db).
!
!!            wall_fac=1.0_r8+gls_E2/(vonKar*vonKar)*                   &
!!   &                 (gls(i,j,k,nstp)**( gls_exp1)*cmu_fac1*          &
!!   &                  tke(i,j,k,nstp)**(-tke_exp1)*                   &
!!   &                  (1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,k))+            &
!!   &                   1.0_r8/(z_w(i,j,k)-z_w(i,j,0))))**2
!!
!! Triangular wall function, L = min (ds, db).
!!
!!            wall_fac=1.0_r8+gls_E2/(vonKar*vonKar)*                   &
!!   &                 (gls(i,j,k,nstp)**( gls_exp1)*cmu_fac1*          &
!!   &                  tke(i,j,k,nstp)**(-tke_exp1)*                   &
!!   &                  (1.0_r8/MIN((z_w(i,j,N(ng))-z_w(i,j,k)),        &
!!   &                              (z_w(i,j,k)-z_w(i,j,0)))))**2
!!
!! Linear wall function for , L = ds (=dist to surface).
!!
!!            wall_fac=1.0_r8+gls_E2/(vonKar*vonKar)*                   &
!!   &                 (gls(i,j,k,nstp)**( gls_exp1)*cmu_fac1*          &
!!   &                  tke(i,j,k,nstp)**(-tke_exp1)*                   &
!!   &                 (1.0_r8/ (z_w(i,j,N(ng))-z_w(i,j,k))))**2
!!
!! Parabolic wall function + free surface correction
!
            wall_fac=1.0_r8+gls_E2/(vonKar*vonKar)*                     &
     &                (gls(i,j,k,nstp)**( gls_exp1)*cmu_fac1*           &
     &                 tke(i,j,k,nstp)**(-tke_exp1)*                    &
     &                 (1.0_r8/ (z_w(i,j,k)-z_w(i,j,0))))**2+           &
     &                0.25_r8/(vonKar*vonKar)*                          &
     &                (gls(i,j,k,nstp)**( gls_exp1)*cmu_fac1*           &
     &                 tke(i,j,k,nstp)**(-tke_exp1)*                    &
     &                 (1.0_r8/ (z_w(i,j,N(ng))-z_w(i,j,k))))**2
            END IF
!
            BCK(i,k)=cff*(1.0_r8+dt(ng)*                                &
     &                    gls(i,j,k,nstp)**(-gls_exp1)*cmu_fac2*        &
     &                    tke(i,j,k,nstp)**( tke_exp2)+                 &
     &                    dt(ng)*(1.0_r8-cff1)*strat2*                  &
     &                    (Akt(i,j,k,itemp)-Akt_bak(itemp,ng))/         &
     &                    tke(i,j,k,nstp))-                             &
     &                    FCK(i,k)-FCK(i,k+1)
            BCP(i,k)=cff*(1.0_r8+dt(ng)*gls_c2(ng)*wall_fac*            &
     &                    gls(i,j,k,nstp)**(-gls_exp1)*cmu_fac2*        &
     &                    tke(i,j,k,nstp)**( tke_exp2)+                 &
     &                    dt(ng)*(1.0_r8-cff2)*gls_c3*strat2*           &
     &                    (Akt(i,j,k,itemp)-Akt_bak(itemp,ng))/         &
     &                    tke(i,j,k,nstp))-                             &
     &                    FCP(i,k)-FCP(i,k+1)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step dissipation and vertical diffusion terms implicitly.
!-----------------------------------------------------------------------
!
!  Set Dirichlet surface and bottom boundary conditions. Compute
!  surface roughness from wind stress (Charnok) and set Craig and
!  Banner wave breaking surface flux, if appropriate.
!
        DO i=Istr,Iend
# if defined CRAIG_BANNER
          tke(i,j,N(ng),nnew)=MAX(cmu_fac4*0.5_r8*                      &
     &                            SQRT((sustr(i,j)+sustr(i+1,j))**2+    &
     &                                 (svstr(i,j)+svstr(i,j+1))**2)*   &
     &                            crgban_cw(ng)**(2.0_r8/3.0_r8),       &
     &                            gls_Kmin(ng))
# elif defined TKE_WAVEDISS
          tke(i,j,N(ng),nnew)=MAX(cmu_fac4*(sz_alpha(ng)*               &
     &                        (Dissip_break(i,j)+Dissip_wcap(i,j)))**   &
     &                        (2.0_r8/3.0_r8), gls_Kmin(ng))
# else
          tke(i,j,N(ng),nnew)=MAX(cmu_fac3*0.5_r8*                      &
     &                            SQRT((sustr(i,j)+sustr(i+1,j))**2+    &
     &                                 (svstr(i,j)+svstr(i,j+1))**2),   &
     &                            gls_Kmin(ng))
# endif
          tke(i,j,0,nnew)=MAX(cmu_fac3*0.5_r8*                          &
     &                        SQRT((bustr(i,j)+bustr(i+1,j))**2+        &
     &                             (bvstr(i,j)+bvstr(i,j+1))**2),       &
     &                        gls_Kmin(ng))
# if defined CHARNOK
          Zos_eff(i)=MAX(charnok_alpha(ng)/g*0.5_r8*                    &
     &                   SQRT((sustr(i,j)+sustr(i+1,j))**2+             &
     &                        (svstr(i,j)+svstr(i,j+1))**2),            &
     &                   Zos_min)
# elif defined ZOS_HSIG
          Zos_eff(i)=MAX(zos_hsig_alpha(ng)*Hwave(i,j), Zos_min)
# else
          Zos_eff(i)=Zos_min
# endif
          gls(i,j,N(ng),nnew)=MAX(gls_cmu0(ng)**gls_p(ng)*              &
     &                            tke(i,j,N(ng),nnew)**gls_m(ng)*       &
     &                            (L_sft*Zos_eff(i))**gls_n(ng),        &
     &                            gls_Pmin(ng))
          cff=gls_fac4*(vonKar*Zob_min(i,j))**(gls_n(ng))
          gls(i,j,0,nnew)=MAX(cff*tke(i,j,0,nnew)**(gls_m(ng)),         &
     &                        gls_Pmin(ng))
        END DO
!
!  Solve tri-diagonal system for turbulent kinetic energy.
!
        DO i=Istr,Iend
# if defined CRAIG_BANNER
          tke_fluxt(i)=dt(ng)*crgban_cw(ng)*                            &
     &                 (0.50_r8*                                        &
     &                  SQRT((sustr(i,j)+sustr(i+1,j))**2+              &
     &                       (svstr(i,j)+svstr(i,j+1))**2))**1.5_r8
# elif defined TKE_WAVEDISS
          tke_fluxt(i)=dt(ng)*sz_alpha(ng)*(Dissip_break(i,j)+          &
     &                                      Dissip_wcap(i,j))
# else
          tke_fluxt(i)=0.0_r8
# endif
          tke_fluxb(i)=0.0_r8
!
          cff=1.0_r8/BCK(i,N(ng)-1)
          CF(i,N(ng)-1)=cff*FCK(i,N(ng)-1)
          tke(i,j,N(ng)-1,nnew)=cff*(tke(i,j,N(ng)-1,nnew)+tke_fluxt(i))
        END DO
        DO i=Istr,Iend
          DO k=N(ng)-2,1,-1
            cff=1.0_r8/(BCK(i,k)-CF(i,k+1)*FCK(i,k+1))
            CF(i,k)=cff*FCK(i,k)
            tke(i,j,k,nnew)=cff*(tke(i,j,k,nnew)-                       &
     &                           FCK(i,k+1)*tke(i,j,k+1,nnew))
          END DO
          tke(i,j,1,nnew)=tke(i,j,1,nnew)-cff*tke_fluxb(i)
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-CF(i,k)*tke(i,j,k-1,nnew)
          END DO
        END DO
!
!  Solve tri-diagonal system for generic statistical field.
!
        DO i=Istr,Iend
          cff=0.5_r8*(tke(i,j,N(ng),nnew)+tke(i,j,N(ng)-1,nnew))
          gls_fluxt(i)=dt(ng)*gls_fac3*cff**gls_m(ng)*                  &
     &                 L_sft**(gls_n(ng))*                              &
     &                 (Zos_eff(i)+0.5_r8*Hz(i,j,N(ng)))**              &
     &                 (gls_n(ng)-1.0_r8)*                              &
     &                 0.5_r8*(Akp(i,j,N(ng))+Akp(i,j,N(ng)-1))
# ifdef CRAIG_BANNER
          gls_fluxt(i)=gls_fluxt(i)-dt(ng)*                             &
     &                 gls_m(ng)*(gls_cmu0(ng)**gls_p(ng))*             &
     &                 cff**(gls_m(ng)-1.0_r8)*                         &
     &                 ((Zos_eff(i)+0.5_r8*Hz(i,j,N(ng)))*L_sft)**      &
     &                 gls_n(ng)*                                       &
     &                 gls_sigk(ng)*ogls_sigp*crgban_cw(ng)*            &
     &                 (0.5_r8*                                         &
     &                  SQRT((sustr(i,j)+sustr(i+1,j))**2+              &
     &                       (svstr(i,j)+svstr(i,j+1))**2))**1.5_r8
# elif defined TKE_WAVEDISS
          gls_fluxt(i)=gls_fluxt(i)-dt(ng)*                             &
     &                 gls_m(ng)*(gls_cmu0(ng)**gls_p(ng))*             &
     &                 cff**(gls_m(ng)-1.0_r8)*                         &
     &                 ((Zos_eff(i)+0.5_r8*Hz(i,j,N(ng)))*L_sft)**      &
     &                 gls_n(ng)*                                       &
     &                 gls_sigk(ng)*ogls_sigp*sz_alpha(ng)*             &
     &                 (Dissip_break(i,j)+Dissip_wcap(i,j))
# endif
          cff=0.5_r8*(tke(i,j,0,nnew)+tke(i,j,1,nnew))
          gls_fluxb(i)=dt(ng)*gls_fac2*(cff**gls_m(ng))*                &
     &                 (0.5_r8*Hz(i,j,1)+Zob_min(i,j))**                &
     &                 (gls_n(ng)-1.0_r8)*                              &
     &                 0.5_r8*(Akp(i,j,0)+Akp(i,j,1))
!
          cff=1.0_r8/BCP(i,N(ng)-1)
          CF(i,N(ng)-1)=cff*FCP(i,N(ng)-1)
          gls(i,j,N(ng)-1,nnew)=cff*(gls(i,j,N(ng)-1,nnew)-gls_fluxt(i))
        END DO
        DO i=Istr,Iend
          DO k=N(ng)-2,1,-1
            cff=1.0_r8/(BCP(i,k)-CF(i,k+1)*FCP(i,k+1))
            CF(i,k)=cff*FCP(i,k)
            gls(i,j,k,nnew)=cff*(gls(i,j,k,nnew)-                       &
     &                           FCP(i,k+1)*gls(i,j,k+1,nnew))
          END DO
          gls(i,j,1,nnew)=gls(i,j,1,nnew)-cff*gls_fluxb(i)
!!        gls(i,j,1,nnew)=MAX(gls(i,j,1,nnew), gls_Pmin(ng))
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-CF(i,k)*gls(i,j,k-1,nnew)
!!          gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew), gls_Pmin(ng))
          END DO
        END DO
!
!---------------------------------------------------------------------
!  Compute vertical mixing coefficients (m2/s).
!---------------------------------------------------------------------
!
        DO i=Istr,Iend
          DO k=1,N(ng)-1
!
!  Compute turbulent length scale (m).
!
            tke(i,j,k,nnew)=MAX(tke(i,j,k,nnew),gls_Kmin(ng))
            gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),gls_Pmin(ng))
            IF (gls_n(ng).ge.0.0_r8) THEN
              gls(i,j,k,nnew)=MIN(gls(i,j,k,nnew),gls_fac5*             &
     &                            tke(i,j,k,nnew)**(tke_exp4)*          &
     &                            (SQRT(MAX(0.0_r8,                     &
     &                                  buoy2(i,j,k)))+eps)**           &
     &                            (-gls_n(ng)))
            ELSE
              gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),gls_fac5*             &
     &                            tke(i,j,k,nnew)**(tke_exp4)*          &
     &                            (SQRT(MAX(0.0_r8,                     &
     &                                  buoy2(i,j,k)))+eps)**           &
     &                            (-gls_n(ng)))
            END IF
            Ls_unlmt=MAX(eps,                                           &
     &                   gls(i,j,k,nnew)**( gls_exp1)*cmu_fac1*         &
     &                   tke(i,j,k,nnew)**(-tke_exp1))
            IF (buoy2(i,j,k).gt.0.0_r8) THEN
              Ls_lmt=MIN(Ls_unlmt,                                      &
     &                 SQRT(0.56_r8*tke(i,j,k,nnew)/                    &
     &                      (MAX(0.0_r8,buoy2(i,j,k))+eps)))
            ELSE
              Ls_lmt=Ls_unlmt
            END IF
!
! Recompute gls based on limited length scale
!
            gls(i,j,k,nnew)=MAX(gls_cmu0(ng)**gls_p(ng)*                &
     &                          tke(i,j,k,nnew)**gls_m(ng)*             &
     &                          Ls_lmt**gls_n(ng), gls_Pmin(ng))
!
!  Compute nondimensional stability functions for tracers (Sh) and
!  momentum (Sm).
!
            Gh=MIN(gls_Gh0,-buoy2(i,j,k)*Ls_lmt*Ls_lmt/                 &
     &                    (2.0_r8*tke(i,j,k,nnew)))
            Gh=MIN(Gh,Gh-(Gh-gls_Ghcri)**2/                             &
     &                    (Gh+gls_Gh0-2.0_r8*gls_Ghcri))
            Gh=MAX(Gh,gls_Ghmin)
# if defined CANUTO_A || defined CANUTO_B
!
!  Compute shear number.
!
            Gm=(gls_b0/gls_fac6-gls_b1*Gh+gls_b3*gls_fac6*(Gh**2))/     &
     &         (gls_b2-gls_b4*gls_fac6*Gh)
            Gm=MIN(Gm,shear2(i,j,k)*Ls_lmt*Ls_lmt/                      &
     &                    (2.0_r8*tke(i,j,k,nnew)))
!!          Gm=MIN(Gm,(gls_s1*gls_fac6*Gh-gls_s0)/(gls_s2*gls_fac6))
!
!  Compute stability functions
!
            cff=gls_b0-gls_b1*gls_fac6*Gh+gls_b2*gls_fac6*Gm+           &
     &          gls_b3*gls_fac6**2*Gh**2-gls_b4*gls_fac6**2*Gh*Gm+      &
     &          gls_b5*gls_fac6**2*Gm*Gm
            Sm=(gls_s0-gls_s1*gls_fac6*Gh+gls_s2*gls_fac6*Gm)/cff
            Sh=(gls_s4-gls_s5*gls_fac6*Gh+gls_s6*gls_fac6*Gm)/cff
            Sm=MAX(Sm,0.0_r8)
            Sh=MAX(Sh,0.0_r8)
!
!  Relate Canuto stability to ROMS notation
!
            Sm=Sm*sqrt2/gls_cmu0(ng)**3
            Sh=Sh*sqrt2/gls_cmu0(ng)**3
# elif defined KANTHA_CLAYSON
            cff=1.0_r8-my_Sh2*Gh
            Sh=my_Sh1/cff
            Sm=(my_B1pm1o3+my_Sm4*Sh*Gh)/(1.0_r8-my_Sm2*Gh)
# else  /* Galperin */
            cff=1.0_r8-my_Sh2*Gh
            Sh=my_Sh1/cff
            Sm=(my_Sm3+Sh*Gh*my_Sm4)/(1.0_r8-my_Sm2*Gh)
# endif
!
!  Compute vertical mixing (m2/s) coefficients of momentum and
!  tracers.  Average ql over the two timesteps rather than using
!  the new Lscale and just averaging tke.
!
            ql=sqrt2*0.5_r8*(Ls_lmt*SQRT(tke(i,j,k,nnew))+              &
     &                       Lscale(i,j,k)*SQRT(tke(i,j,k,nstp)))
            Akv(i,j,k)=Akv_bak(ng)+Sm*ql
#ifdef AKLIMIT
            Akv(i,j,k)=MIN(Ak_max,Akv(i,j,k))
#endif
            DO itrc=1,NAT
              Akt(i,j,k,itrc)=Akt_bak(itrc,ng)+Sh*ql
#ifdef AKLIMIT
              Akt(i,j,k,itrc)=MIN(Ak_max,Akt(i,j,k,itrc))
#endif
            END DO
!
!  Compute vertical mixing (m2/s) coefficents of turbulent kinetic
!  energy and generic statistical field.
!
            Akk(i,j,k)=Akk_bak(ng)+                                     &
     &                 Sm*ql/gls_sigk(ng)

# if defined CRAIG_BANNER || defined TKE_WAVEDISS
!
!  If wave breaking, modify surface boundary condition for
!  gls diffusivity Schmidt number.
!
            Pprod=gls_c1(ng)*shear2(i,j,k)*Akv(i,j,k)
            cff=cmu_fac2*tke(i,j,k,nnew)**(1.5_r8+tke_exp1)*            &
     &          gls(i,j,k,nnew)**(-1.0_r8/gls_n(ng))
            cff2=MIN(Pprod/cff, 1.0_r8)
            sig_eff=cff2*gls_sigp(ng)+(1.0_r8-cff2)*gls_sigp_cb
            Akp(i,j,k)=Akp_bak(ng)+Sm*ql/sig_eff
# else
            Akp(i,j,k)=Akp_bak(ng)+Sm*ql*ogls_sigp
# endif
!
!  Save limited length scale.
!
            Lscale(i,j,k)=Ls_lmt
          END DO
!
!  Compute vertical mixing coefficients at the surface and bottom.
!
!
          Akv(i,j,N(ng))=Akv_bak(ng)+L_sft*Zos_eff(i)*gls_cmu0(ng)*     &
     &                   SQRT(tke(i,j,N(ng),nnew))
#ifdef AKLIMIT
          Akv(i,j,N(ng))=MIN(Ak_max,Akv(i,j,N(ng)))
#endif
          Akv(i,j,0)=Akv_bak(ng)+vonKar*Zob_min(i,j)*gls_cmu0(ng)*      &
     &               SQRT(tke(i,j,0,nnew))
#ifdef AKLIMIT
          Akv(i,j,0)=MIN(Ak_max,Akv(i,j,0))
#endif
!
          Akk(i,j,N(ng))=Akk_bak(ng)+Akv(i,j,N(ng))/gls_sigk(ng)
          Akk(i,j,0)=Akk_bak(ng)+Akv(i,j,0)/gls_sigk(ng)
          Akp(i,j,N(ng))=Akp_bak(ng)+Akv(i,j,N(ng))*ogls_sigp
          Akp(i,j,0)=Akp_bak(ng)+Akv(i,j,0)/gls_sigp(ng)
!
          DO itrc=1,NAT
            Akt(i,j,N(ng),itrc)=Akt_bak(itrc,ng)
            Akt(i,j,0,itrc)=Akt_bak(itrc,ng)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL tkebc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nnew, nstp,                                      &
     &                 gls, tke)

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tke(:,:,:,nnew))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          gls(:,:,:,nnew))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akv)
        DO itrc=1,NAT
          CALL exchange_w3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            Akt(:,:,:,itrc))
        END DO
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tke(:,:,:,nnew),                              &
     &                    gls(:,:,:,nnew),                              &
     &                    Akv)
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng), 1, NAT,         &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akt)
# endif

      RETURN
      END SUBROUTINE gls_corstep_tile
#endif
      END MODULE gls_corstep_mod
