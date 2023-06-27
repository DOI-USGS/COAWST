#include "cppdefs.h"
#ifdef BBL_MODEL
# undef BSTRESS_UPWIND
#else
# define BSTRESS_UPWIND
#endif
#define SED_WENO
#undef  SED_UPWIND

      MODULE sed_bedload_mod

#if defined NONLINEAR && defined SEDIMENT && \
   (defined BEDLOAD_SOULSBY || defined BEDLOAD_MPM) 
!
!git $Id$
!svn $Id: sed_bedload.F 1151 2023-02-09 03:08:53Z arango $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes sediment bedload transport using the Meyer-   !
!  Peter and Muller (1948) formulation  for unidirectional flow and    !
!  Souksby and Damgaard (2005) algorithm that accounts for combined    !
!  effect of currents and waves.                                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Meyer-Peter, E. and R. Muller, 1948: Formulas for bedload transport !
!    In: Report on the 2nd Meeting International Association Hydraulic !
!    Research, Stockholm, Sweden, pp 39-64.                            !
!                                                                      !
!  Soulsby, R.L. and J.S. Damgaard, 2005: Bedload sediment transport   !
!    in coastal waters, Coastal Engineering, 52 (8), 673-689.          !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!  Udated sed bed evolution scheme to the WENO method of:              !
!  Wen Long, James T. Kirby, Zhiyu Shao,                               !
!  A numerical scheme for morphological bed level calculations,        !
!  Coastal Engineering,55, Issue 2, 2008, 167-180.                     !
!  https://doi.org/10.1016/j.coastaleng.2007.09.009.                   !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: sed_bedload
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_bedload (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
      USE mod_stepping
# ifdef BBL_MODEL
      USE mod_bbl
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 16, __LINE__, MyFile)
# endif
      CALL sed_bedload_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp(ng), nnew(ng),                        &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef WET_DRY
     &                       GRID(ng) % rmask_wet,                      &
     &                       GRID(ng) % umask_wet,                      &
     &                       GRID(ng) % vmask_wet,                      &
# endif
     &                       GRID(ng) % z_w,                            &
# ifdef BBL_MODEL
     &                       BBL(ng) % bustrc,                          &
     &                       BBL(ng) % bvstrc,                          &
     &                       BBL(ng) % bustrw,                          &
     &                       BBL(ng) % bvstrw,                          &
     &                       BBL(ng) % bustrcwmax,                      &
     &                       BBL(ng) % bvstrcwmax,                      &
     &                       FORCES(ng) % Dwave,                        &
     &                       FORCES(ng) % Pwave_bot,                    &
# endif
     &                       FORCES(ng) % bustr,                        &
     &                       FORCES(ng) % bvstr,                        &
# if defined BEDLOAD_SOULSBY
     &                       FORCES(ng) % Hwave,                        &
     &                       FORCES(ng) % Lwave,                        &
     &                       GRID(ng) % angler,                         &
# endif
# if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % om_r,                           &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % on_r,                           &
     &                       GRID(ng) % on_v,                           &
     &                       SEDBED(ng) % bedldu,                       &
     &                       SEDBED(ng) % bedldv,                       &
# endif
     &                       SEDBED(ng) % bed,                          &
     &                       SEDBED(ng) % bed_frac,                     &
     &                       SEDBED(ng) % bed_mass,                     &
     &                       SEDBED(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE sed_bedload
!
!***********************************************************************
      SUBROUTINE sed_bedload_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew,                          &
     &                             pm, pn,                              &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef WET_DRY
     &                             rmask_wet, umask_wet, vmask_wet,     &
# endif
     &                             z_w,                                 &
# ifdef BBL_MODEL
     &                             bustrc, bvstrc,                      &
     &                             bustrw, bvstrw,                      &
     &                             bustrcwmax, bvstrcwmax,              &
     &                             Dwave, Pwave_bot,                    &
# endif
     &                             bustr, bvstr,                        &
# if defined BEDLOAD_SOULSBY
     &                             Hwave, Lwave,                        &
     &                             angler,                              &
# endif
# if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY
     &                             h, om_r, om_u, on_r, on_v,           &
     &                             bedldu, bedldv,                      &
# endif
     &                             bed, bed_frac, bed_mass,             &
     &                             bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_sediment
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
# ifdef BEDLOAD
      USE exchange_2d_mod, ONLY : exchange_u2d_tile, exchange_v2d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrc(LBi:,LBj:)
      real(r8), intent(in) :: bvstrc(LBi:,LBj:)
      real(r8), intent(in) :: bustrw(LBi:,LBj:)
      real(r8), intent(in) :: bvstrw(LBi:,LBj:)
      real(r8), intent(in) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: bvstrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
#  if defined BEDLOAD_SOULSBY
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
#  endif
#  if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(inout) :: bedldu(LBi:,LBj:,:)
      real(r8), intent(inout) :: bedldv(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustrw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pwave_bot(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
#  if defined BEDLOAD_SOULSBY
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
#  endif
#  if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bedldu(LBi:UBi,LBj:UBj,NST)
      real(r8), intent(inout) :: bedldv(LBi:UBi,LBj:UBj,NST)
#  endif
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
# endif
!
!  Local variable declarations.
!
      integer :: i, ised, j, k

      real(r8), parameter :: eps = 1.0E-14_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, fac1, fac2
      real(r8) :: Dstp, bed_change, dz, roll

# if defined BEDLOAD_MPM
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_w
!     nondimensional critical erosion stress for MPM
      real(r8), parameter :: tau_mpmc = 0.047_r8
# endif
# ifdef BSTRESS_UPWIND
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_wX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_wE
# endif
# ifdef SED_SLUMP
      real(r8) :: slopefac_wet_local, slopefac_dry_local
# endif
# ifdef SED_WENO
      real(r8) :: S1m, S2m, S3m, S1p, S2p, S3p
      real(r8) :: alpha1m, alpha2m, alpha3m
      real(r8) :: alpha1p, alpha2p, alpha3p, alpham, alphap
      real(r8) :: w1m, w2m, w3m, w1p, w2p, w3p
      real(r8) :: q1m, q2m, q3m, q1p, q2p, q3p 
      real(r8) :: signa, FXm, FXp, FEm, FEp
      real(r8), parameter :: thirtotwelv = 13.0_r8/12.0_r8
      real(r8), parameter :: elevenosix  = 11.0_r8/6.0_r8
      real(r8), parameter :: sevenosix   = 7.0_r8/6.0_r8
      real(r8), parameter :: fiveosix    = 5.0_r8/6.0_r8
      real(r8), parameter :: oneosix     = 1.0_r8/6.0_r8
      real(r8), parameter :: oneothree   = 1.0_r8/3.0_r8
# endif
# ifdef BEDLOAD
      real(r8) :: bedld, bedld_mass, dzdx, dzdy, dzdxdy
      real(r8) :: a_slopex, a_slopey, sed_angle
      real(r8) :: smgd, smgdr, osmgd, Umag
      real(r8) :: rhs_bed, Ua, Ra, phi, Clim

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE_r
# endif
# if defined BEDLOAD_MPM && !defined BSTRESS_UPWIND

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: angleu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: anglev
# endif
# if defined BEDLOAD_SOULSBY
      real(r8) :: theta_mean, theta_wav, w_asym
      real(r8) :: theta_max, theta_max1, theta_max2
      real(r8) :: phi_x1, phi_x2, phi_x, phi_y
      real(r8) :: bedld_x, bedld_y, tau_cur, waven, wavec

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phic
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phicw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_wav
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_mean
      real(r8), parameter :: kdmax = 100.0_r8
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! Compute maximum bottom stress for MPM bedload or suspended load.
!-----------------------------------------------------------------------
!
# if defined BEDLOAD_MPM
#  ifdef BBL_MODEL
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          tau_w(i,j)=SQRT(bustrcwmax(i,j)*bustrcwmax(i,j)+              &
     &                    bvstrcwmax(i,j)*bvstrcwmax(i,j))
        END DO
      END DO
#  else
      DO j=Jstrm1,Jendp1
        DO i=Istrm1,Iendp1
#   ifdef BSTRESS_UPWIND
          cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,bustr(i+1,j)))
          cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,bustr(i+1,j)))
          cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,bustr(i  ,j)))
          cff4=0.5_r8*(1.0_r8-SIGN(1.0_r8,bustr(i  ,j)))
          tau_wX(i,j)=cff3*(cff1*bustr(i,j)+                            &
     &                cff2*0.5_r8*(bustr(i,j)+bustr(i+1,j)))+           &
     &                cff4*(cff2*bustr(i+1,j)+                          &
     &                cff1*0.5_r8*(bustr(i,j)+bustr(i+1,j)))
          cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,bvstr(i,j+1)))
          cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,bvstr(i,j+1)))
          cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,bvstr(i,j)))
          cff4=0.5_r8*(1.0_r8-SIGN(1.0_r8,bvstr(i,j)))
          tau_wE(i,j)=cff3*(cff1*bvstr(i,j)+                            &
     &                cff2*0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))+           &
     &                cff4*(cff2*bvstr(i,j+1)+                          &
     &                cff1*0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))
#   endif
          tau_w(i,j)=0.5_r8*SQRT((bustr(i,j)+bustr(i+1,j))*             &
     &                           (bustr(i,j)+bustr(i+1,j))+             &
     &                           (bvstr(i,j)+bvstr(i,j+1))*             &
     &                           (bvstr(i,j)+bvstr(i,j+1)))
        END DO
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute bedload sediment transport.
!-----------------------------------------------------------------------
!
! Compute some constant bed slope parameters.
!
      sed_angle=DTAN(33.0_r8*pi/180.0_r8)
!
!  Compute angle between currents and waves (radians).
!
# if defined SED_WENO
      DO j=Jstrm3,Jendp2i
        DO i=Istrm3,Iendp2i
# else
      DO j=Jstrm1,Jendp1
        DO i=Istrm1,Iendp1
# endif
# if defined BEDLOAD_SOULSBY
!
! Compute angle between currents and waves, measure CCW from current
! direction toward wave vector.
!
          IF (bustrc(i,j).eq.0.0_r8) THEN
            phic(i,j)=0.5_r8*pi*SIGN(1.0_r8,bvstrc(i,j))
          ELSE
            phic(i,j)=ATAN2(bvstrc(i,j),bustrc(i,j))
          ENDIF
          phicw(i,j)=1.5_r8*pi-Dwave(i,j)-phic(i,j)-angler(i,j)
!
! Compute stress components at rho points.
!
          tau_cur=SQRT(bustrc(i,j)*bustrc(i,j)+                         &
     &                 bvstrc(i,j)*bvstrc(i,j))
          tau_wav(i,j)=MIN(SQRT(bustrw(i,j)*bustrw(i,j)+                &
     &                          bvstrw(i,j)*bvstrw(i,j)),20.0_r8)
          tau_mean(i,j)=tau_cur*(1.0_r8+1.2_r8*((tau_wav(i,j)/          &
     &                  (tau_cur+tau_wav(i,j)+eps))**3.2_r8))
!
# elif defined BEDLOAD_MPM && !defined BSTRESS_UPWIND
          cff1=0.5_r8*(bustr(i,j)+bustr(i+1,j))
          cff2=0.5_r8*(bvstr(i,j)+bvstr(i,j+1))
          Umag=SQRT(cff1*cff1+cff2*cff2)+eps
          angleu(i,j)=cff1/Umag
          anglev(i,j)=cff2/Umag
# endif
        END DO
      END DO
!
      DO ised=NCS+1,NST
        smgd=(Srho(ised,ng)/rho0-1.0_r8)*g*Sd50(ised,ng)
        osmgd=1.0_r8/smgd
        smgdr=SQRT(smgd)*Sd50(ised,ng)*Srho(ised,ng)
!
# if defined SED_WENO
        DO j=Jstrm3,Jendp2i
          DO i=Istrm3,Iendp2i
# else
        DO j=Jstrm1,Jendp1
          DO i=Istrm1,Iendp1
# endif
# ifdef BEDLOAD_SOULSBY
!
! Compute wave asymmetry factor, based on Fredosoe and Deigaard.
!
            Dstp=z_w(i,j,N(ng))+h(i,j)
            waven=2.0_r8*pi/(Lwave(i,j)+eps)
!           wavec=SQRT(g/waven*tanh(waven*Dstp))
            cff4=MIN(waven*Dstp,kdmax)
!           cff1=-0.1875_r8*wavec*(waven*Dstp)**2/(SINH(cff4))**4
!           cff2=0.125_r8*g*Hwave(i,j)**2/(wavec*Dstp+eps)
!           cff3=pi*Hwave(i,j)/(Pwave_bot(i,j)*SINH(cff4)+eps)
!
! Compute wave asymmetry factor, based on the note of Soulsby
!
            cff1=MIN(0.375_r8*(Hwave(i,j)/Dstp)*                        &
     &                    ((waven*Dstp)/(SINH(cff4))**3),0.15_r8)
            w_asym=2.0_r8*cff1/(1.0_r8+cff1**2.0_r8)
!
! Compute nondimensional stresses.
!
            theta_wav=tau_wav(i,j)*osmgd+eps
            theta_mean=tau_mean(i,j)*osmgd
!
            cff1=theta_wav*(1.0_r8+w_asym)
            cff2=theta_wav*(1.0_r8-w_asym)
            theta_max1=SQRT((theta_mean+                                &
     &                       cff1*COS(phicw(i,j)))**2+                  &
     &                      (cff1*SIN(phicw(i,j)))**2)
            theta_max2=SQRT((theta_mean+                                &
     &                       cff2*COS(phicw(i,j)+pi))**2+               &
     &                      (cff2*SIN(phicw(i,j)+pi))**2)
            theta_max=MAX(theta_max1,theta_max2)
!
! Motion initiation factor.
!
#  if defined COHESIVE_BED || defined MIXED_BED
! Check that both sediment class and bed critical stresses are exceeded
            cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,                            &
     &                   theta_max/                                     &
     &          (max(tau_ce(ised,ng),bed(i,j,1,ibtcr))*osmgd)-1.0_r8))
	
#  else
! Only sediment class critical stresses needs to be exceeded
            cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,                            &
     &                   theta_max/(tau_ce(ised,ng)*osmgd)-1.0_r8))
#  endif
!
! Calculate bed loads in direction of current and perpendicular
! direction.
!
            phi_x1=12.0_r8*SQRT(theta_mean)*                            &
     &             MAX((theta_mean-tau_ce(ised,ng)*osmgd),0.0_r8)
            phi_x2=12.0_r8*(0.9534_r8+0.1907*COS(2.0_r8*phicw(i,j)))*   &
     &             SQRT(theta_wav)*theta_mean+                          &
     &             12.0_r8*(0.229_r8*w_asym*theta_wav**1.5_r8*          &
     &                      COS(phicw(i,j)))
!           phi_x=MAX(phi_x1,phi_x2) !  <- original
            IF (ABS(phi_x2).gt.phi_x1) THEN
              phi_x=phi_x2
            ELSE
              phi_x=phi_x1
            END IF
            bedld_x=phi_x*smgdr*cff3
!
            cff5=theta_wav**1.5_r8+1.5_r8*(theta_mean**1.5_r8)
            phi_y=12.0_r8*0.1907_r8*theta_wav*theta_wav*                &
     &            (theta_mean*SIN(2.0_r8*phicw(i,j))+1.2_r8*w_asym*     &
     &            theta_wav*SIN(phicw(i,j)))/cff5*cff3
            bedld_y=phi_y*smgdr
!
! Partition bedld into xi and eta directions, still at rho points.
! (FX_r and FE_r have dimensions of kg).
!
            FX_r(i,j)=(bedld_x*COS(phic(i,j))-bedld_y*SIN(phic(i,j)))*  &
     &                on_r(i,j)*dt(ng)
            FE_r(i,j)=(bedld_x*SIN(phic(i,j))+bedld_y*COS(phic(i,j)))*  &
     &                om_r(i,j)*dt(ng)
# elif defined BEDLOAD_MPM
#  ifdef BSTRESS_UPWIND
!
! Magnitude of bed load at rho points. Meyer-Peter Muller formulation.
! bedld has dimensions of kg m-1 s-1. Use partitions of stress
! from upwind direction, still at rho points.
! (FX_r and FE_r have dimensions of kg).
!
#   if defined COHESIVE_BED || defined MIXED_BED
            IF (tau_wX(i,j).gt.bed(i,j,1,ibtcr)) THEN
               cff=1.0_r8
            ELSE
               cff=0.0_r8
            END
            bedld=8.0_r8*(MAX((ABS(tau_wX(i,j))*osmgd-tau_mpmc),        &
     &                        0.0_r8)**1.5_r8)*smgdr*                   &
     &                        SIGN(1.0_r8,tau_wX(i,j))*cff
            FX_r(i,j)=bedld*on_r(i,j)*dt(ng)
            IF (tau_wE(i,j).gt.bed(i,j,1,ibtcr)) THEN
               cff=1.0_r8
            ELSE
               cff=0.0_r8
            END
            bedld=8.0_r8*(MAX((ABS(tau_wE(i,j))*osmgd-tau_mpmc),        &
     &                        0.0_r8)**1.5_r8)*smgdr*                   &
     &                        SIGN(1.0_r8,tau_wE(i,j))*cff
            FE_r(i,j)=bedld*om_r(i,j)*dt(ng)
#   else
            bedld=8.0_r8*(MAX((ABS(tau_wX(i,j))*osmgd-tau_mpmc),        &
     &                        0.0_r8)**1.5_r8)*smgdr*                   &
     &                        SIGN(1.0_r8,tau_wX(i,j))
            FX_r(i,j)=bedld*on_r(i,j)*dt(ng)
            bedld=8.0_r8*(MAX((ABS(tau_wE(i,j))*osmgd-tau_mpmc),        &
     &                        0.0_r8)**1.5_r8)*smgdr*                   &
     &                        SIGN(1.0_r8,tau_wE(i,j))
            FE_r(i,j)=bedld*om_r(i,j)*dt(ng)
#   endif
#  else
!
! Magnitude of bed load at rho points. Meyer-Peter Muller formulation.
! (BEDLD has dimensions of kg m-1 s-1).
!
#   if defined COHESIVE_BED || defined MIXED_BED
            IF (tau_w(i,j).gt.bed(i,j,1,ibtcr)) THEN
              bedld=8.0_r8*(MAX((tau_w(i,j)*osmgd-tau_mpmc),            &
     &                           0.0_r8)**1.5_r8)*smgdr
            ELSE
              bedld=0.0_r8
            END IF
#   else
            bedld=8.0_r8*(MAX((tau_w(i,j)*osmgd-tau_mpmc),              &
     &                         0.0_r8)**1.5_r8)*smgdr
#   endif
!
! Partition bedld into xi and eta directions, still at rho points.
! (FX_r and FE_r have dimensions of kg).
!
            FX_r(i,j)=angleu(i,j)*bedld*on_r(i,j)*dt(ng)
            FE_r(i,j)=anglev(i,j)*bedld*om_r(i,j)*dt(ng)
#  endif
# endif
!
! Correct for along-direction slope. Limit slope to 0.9*sed angle.
!
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i,j)))
            cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,FE_r(i,j)))
            cff4=0.5_r8*(1.0_r8-SIGN(1.0_r8,FE_r(i,j)))
# if defined SLOPE_NEMETH
            dzdx=(h(i+1,j)-h(i,j))/om_u(i+1,j)*cff1+                    &
     &           (h(i-1,j)-h(i,j))/om_u(i  ,j)*cff2
            dzdy=(h(i,j+1)-h(i,j))/on_v(i,j+1)*cff3+                    &
     &           (h(i,j-1)-h(i,j))/on_v(i  ,j)*cff4
            a_slopex=1.7_r8*dzdx
            a_slopey=1.7_r8*dzdy
!
! Add contriubiton of bed slope to bed load transport fluxes.
!
            FX_r(i,j)=FX_r(i,j)*(1.0_r8+a_slopex)
            FE_r(i,j)=FE_r(i,j)*(1.0_r8+a_slopey)
!
# elif defined SLOPE_LESSER
            dzdx=MIN(((h(i+1,j)-h(i  ,j))/om_u(i+1,j)*cff1+             &
     &                (h(i  ,j)-h(i-1,j))/om_u(i  ,j)*cff2),0.52_r8)*   &
     &                SIGN(1.0_r8,FX_r(i,j))
            dzdy=MIN(((h(i,j+1)-h(i,j  ))/on_v(i,j+1)*cff3+             &
     &                (h(i,j  )-h(i,j-1))/on_v(i  ,j)*cff4),0.52_r8)*   &
     &                SIGN(1.0_r8,FE_r(i,j))
            cff=DATAN(dzdx)
            a_slopex=sed_angle/(COS(cff)*(sed_angle-dzdx))
            cff=DATAN(dzdy)
            a_slopey=sed_angle/(COS(cff)*(sed_angle-dzdy))
!
! Add contriubiton of bed slope to bed load transport fluxes.
!
            FX_r(i,j)=FX_r(i,j)*a_slopex
            FE_r(i,j)=FE_r(i,j)*a_slopey
# endif
# ifdef SED_MORPH
!
! Apply morphology factor.
!
            FX_r(i,j)=FX_r(i,j)*morph_fac(ised,ng)
            FE_r(i,j)=FE_r(i,j)*morph_fac(ised,ng)
# endif
!
! Apply bedload transport rate coefficient. Also limit
! bedload to the fraction of each sediment class.
!
            FX_r(i,j)=FX_r(i,j)*bedload_coeff(ng)*bed_frac(i,j,1,ised)
            FE_r(i,j)=FE_r(i,j)*bedload_coeff(ng)*bed_frac(i,j,1,ised)
!
! Limit bed load to available bed mass.
!
            bedld_mass=ABS(FX_r(i,j))+ABS(FE_r(i,j))+eps
            FX_r(i,j)=MIN(ABS(FX_r(i,j)),                               &
     &                    bed_mass(i,j,1,nstp,ised)*                    &
     &                    om_r(i,j)*on_r(i,j)*ABS(FX_r(i,j))/           &
     &                    bedld_mass)*                                  &
     &                    SIGN(1.0_r8,FX_r(i,j))
            FE_r(i,j)=MIN(ABS(FE_r(i,j)),                               &
     &                    bed_mass(i,j,1,nstp,ised)*                    &
     &                    om_r(i,j)*on_r(i,j)*ABS(FE_r(i,j))/           &
     &                    bedld_mass)*                                  &
     &                    SIGN(1.0_r8,FE_r(i,j))
# ifdef MASKING
#  ifdef WET_DRY
            FX_r(i,j)=FX_r(i,j)*rmask_wet(i,j)
            FE_r(i,j)=FE_r(i,j)*rmask_wet(i,j)
#  else
            FX_r(i,j)=FX_r(i,j)*rmask(i,j)
            FE_r(i,j)=FE_r(i,j)*rmask(i,j)
#  endif
# endif
          END DO
        END DO
!
!  Apply boundary conditions (gradient).
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstrm1,Jendp1
              FX_r(Istr-1,j)=FX_r(Istr,j)
              FE_r(Istr-1,j)=FE_r(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstrm1,Jendp1
              FX_r(Iend+1,j)=FX_r(Iend,j)
              FE_r(Iend+1,j)=FE_r(Iend,j)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istrm1,Iendp1
              FX_r(i,Jstr-1)=FX_r(i,Jstr)
              FE_r(i,Jstr-1)=FE_r(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istrm1,Iendp1
              FX_r(i,Jend+1)=FX_r(i,Jend)
              FE_r(i,Jend+1)=FE_r(i,Jend)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            FX_r(Istr-1,Jstr-1)=0.5_r8*(FX_r(Istr  ,Jstr-1)+            &
     &                                  FX_r(Istr-1,Jstr  ))
            FE_r(Istr-1,Jstr-1)=0.5_r8*(FE_r(Istr  ,Jstr-1)+            &
     &                                  FE_r(Istr-1,Jstr  ))
          END IF
        END IF

        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            FX_r(Iend+1,Jstr-1)=0.5_r8*(FX_r(Iend  ,Jstr-1)+            &
     &                                  FX_r(Iend+1,Jstr  ))
            FE_r(Iend+1,Jstr-1)=0.5_r8*(FE_r(Iend  ,Jstr-1)+            &
     &                                  FE_r(Iend+1,Jstr  ))
          END IF
        END IF

        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            FX_r(Istr-1,Jend+1)=0.5_r8*(FX_r(Istr-1,Jend  )+            &
     &                                  FX_r(Istr  ,Jend+1))
            FE_r(Istr-1,Jend+1)=0.5_r8*(FE_r(Istr-1,Jend  )+            &
     &                                  FE_r(Istr  ,Jend+1))
          END IF
        END IF

        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            FX_r(Iend+1,Jend+1)=0.5_r8*(FX_r(Iend+1,Jend  )+            &
     &                                  FX_r(Iend  ,Jend+1))
            FE_r(Iend+1,Jend+1)=0.5_r8*(FE_r(Iend+1,Jend  )+            &
     &                                  FE_r(Iend  ,Jend+1))
          END IF
        END IF
!
! Compute face fluxes at u and v points before taking divergence.
!
!       DO j=JstrR,JendR
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
# if defined SED_UPWIND
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i,j)))
            FX(i,j)=0.5_r8*(1.0_r8+SIGN(1.0_r8,FX_r(i-1,j)))*           &
     &              (cff1*FX_r(i-1,j)+                                  &
     &               cff2*0.5_r8*(FX_r(i-1,j)+FX_r(i,j)))+              &
     &              0.5_r8*(1.0_r8-SIGN(1.0_r8,FX_r(i-1,j)))*           &
     &              (cff2*FX_r(i  ,j)+                                  &
     &               cff1*0.5_r8*(FX_r(i-1,j)+FX_r(i,j)))
# elif defined SED_WENO
!
!  Long et al. (2008). Coastal Engr, 55, 167-180.
!
            S1m=thirtotwelv*                                            &
     &          (FX_r(i-3,j)-2.0_r8*FX_r(i-2,j)+FX_r(i-1,j))**2+        &
     &          0.25_r8*                                                &
     &          (FX_r(i-3,j)-4.0_r8*FX_r(i-2,j)+3.0_r8*FX_r(i-1,j))**2
            S2m=thirtotwelv*                                            &
     &          (FX_r(i-2,j)-2.0_r8*FX_r(i-1,j)+FX_r(i,j))**2+          &
     &          0.25_r8*                                                &
     &          (FX_r(i-2,j)-FX_r(i,j))**2
            S3m=thirtotwelv*                                            &
     &          (FX_r(i-1,j)-2.0_r8*FX_r(i,j)+FX_r(i+1,j))**2+          &
     &          0.25_r8*                                                &
     &          (3.0_r8*FX_r(i-1,j)-4.0_r8*FX_r(i,j)+FX_r(i+1,j))**2
!
            S1p=thirtotwelv*                                            &
     &          (FX_r(i-2,j)-2.0_r8*FX_r(i-1,j)+FX_r(i,j))**2+          &
     &          0.25_r8*                                                &
     &          (FX_r(i-2,j)-4.0_r8*FX_r(i-1,j)+3.0_r8*FX_r(i,j))**2
            S2p=thirtotwelv*                                            &
     &          (FX_r(i-1,j)-2.0_r8*FX_r(i,j)+FX_r(i+1,j))**2+          &
     &          0.25_r8*                                                &
     &          (FX_r(i-1,j)-FX_r(i+1,j))**2
!     &          (-FX_r(i-1,j)+FX_r(i+1,j))**2
            S3p=thirtotwelv*                                            &
     &          (FX_r(i,j)-2.0_r8*FX_r(i+1,j)+FX_r(i+2,j))**2+          &
     &          0.25_r8*                                                &
     &          (3.0_r8*FX_r(i,j)-4.0_r8*FX_r(i+1,j)+FX_r(i+2,j))**2


!
            alpha1m=0.1_r8/(S1m+eps)**2
            alpha2m=0.6_r8/(S2m+eps)**2
            alpha3m=0.3_r8/(S3m+eps)**2
!
            alpha1p=0.3_r8/(S1p+eps)**2
            alpha2p=0.6_r8/(S2p+eps)**2
            alpha3p=0.1_r8/(S3p+eps)**2
!
            alpham=alpha1m+alpha2m+alpha3m
            alphap=alpha1p+alpha2p+alpha3p
!
            w1m=alpha1m/alpham
            w2m=alpha2m/alpham
            w3m=alpha3m/alpham
            w1p=alpha1p/alphap
            w2p=alpha2p/alphap
            w3p=alpha3p/alphap
!
            q1m=oneothree*FX_r(i-3,j)-sevenosix*FX_r(i-2,j)+            &
     &          elevenosix*FX_r(i-1,j)
            q2m=-oneosix*FX_r(i-2,j)+fiveosix*FX_r(i-1,j)+              &
     &           oneothree*FX_r(i,j)
            q3m=oneothree*FX_r(i-1,j)+fiveosix*FX_r(i,j)-               &
     &          oneosix*FX_r(i+1,j)
!
            q1p=-oneosix*FX_r(i-2,j)+fiveosix*FX_r(i-1,j)+              &
     &           oneothree*FX_r(i,j)
            q2p=oneothree*FX_r(i-1,j)+fiveosix*FX_r(i,j)-               &
     &          oneosix*FX_r(i+1,j)
            q3p=elevenosix*FX_r(i,j)-sevenosix*FX_r(i+1,j)+             &
     &          oneothree*FX_r(i+2,j)
!
!           signa=(FX_r(i,j)-FX_r(i-1,j))*(h(i,j)-h(i-1,j)) 
            signa=FX_r(i,j)
            cff=SIGN(1.0_r8,signa)
            FXm=0.5_r8*(1.0_r8+cff)*(w1m*q1m+w2m*q2m+w3m*q3m)
            FXp=0.5_r8*(1.0_r8-cff)*(w1p*q1p+w2p*q2p+w3p*q3p)
!
            FX(i,j)=FXm+FXp
# else
            FX(i,j)=0.5_r8*(FX_r(i-1,j)+FX_r(i,j))
# endif
# ifdef SLOPE_KIRWAN
!           cff1=30.0_r8
            cff1=10.0_r8
            dzdx=(h(i,j)-h(i-1  ,j))/om_u(i,j)
            a_slopex=(MAX(0.0_r8,abs(dzdx)-0.05_r8)                     &
     &               *SIGN(1.0_r8,dzdx)*cff1)                           &
     &               *om_r(i,j)*dt(ng)
#  ifdef SED_MORPH
            a_slopex=a_slopex*morph_fac(ised,ng)
#  endif
            FX(i,j)=FX(i,j)+a_slopex
# endif
# ifdef MASKING
            FX(i,j)=FX(i,j)*umask(i,j)
#  ifdef WET_DRY
            FX(i,j)=FX(i,j)*umask_wet(i,j)
#  endif
# endif
          END DO
        END DO
        DO j=Jstr,Jend+1
!         DO i=IstrR,IendR
          DO i=Istr,Iend
# ifdef SED_UPWIND
            cff1=0.5_r8*(1.0_r8+SIGN(1.0_r8,FE_r(i,j)))
            cff2=0.5_r8*(1.0_r8-SIGN(1.0_r8,FE_r(i,j)))
            FE(i,j)=0.5_r8*(1.0_r8+SIGN(1.0_r8,FE_r(i,j-1)))*           &
     &              (cff1*FE_r(i,j-1)+                                  &
     &               cff2*0.5_r8*(FE_r(i,j-1)+FE_r(i,j)))+              &
     &              0.5_r8*(1.0_r8-SIGN(1.0_r8,FE_r(i,j-1)))*           &
     &              (cff2*FE_r(i  ,j)+                                  &
     &               cff1*0.5_r8*(FE_r(i,j-1)+FE_r(i,j)))
# elif defined SED_WENO
            S1m=thirtotwelv*                                            &
     &          (FE_r(i,j-3)-2.0_r8*FE_r(i,j-2)+FE_r(i,j-1))**2+        &
     &          0.25_r8*                                                &
     &          (FE_r(i,j-3)-4.0_r8*FE_r(i,j-2)+3.0_r8*FE_r(i,j-1))**2
            S2m=thirtotwelv*                                            &
     &          (FE_r(i,j-2)-2.0_r8*FE_r(i,j-1)+FE_r(i,j))**2+          &
     &          0.25_r8*                                                &
     &          (FE_r(i,j-2)-FE_r(i,j))**2
            S3m=thirtotwelv*                                            &
     &          (FE_r(i,j-1)-2.0_r8*FE_r(i,j)+FE_r(i,j+1))**2+          &
     &          0.25_r8*                                                &
     &          (3.0_r8*FE_r(i,j-1)-4.0_r8*FE_r(i,j)+FE_r(i,j+1))**2
!
            S1p=thirtotwelv*                                            &
     &          (FE_r(i,j-2)-2.0_r8*FE_r(i,j-1)+FE_r(i,j))**2+          &
     &          0.25_r8*                                                &
     &          (FE_r(i,j-2)-4.0_r8*FE_r(i,j-1)+3.0_r8*FE_r(i,j))**2
            S2p=thirtotwelv*                                            &
     &          (FE_r(i,j-1)-2.0_r8*FE_r(i,j)+FE_r(i,j+1))**2+          &
     &          0.25_r8*                                                &
     &          (FE_r(i,j-1)-FE_r(i,j+1))**2
!     &          (-FE_r(i,j-1)+FE_r(i,j+1))**2
            S3p=thirtotwelv*                                            &
     &          (FE_r(i,j)-2.0_r8*FE_r(i,j+1)+FE_r(i,j+2))**2+          &
     &          0.25_r8*                                                &
     &          (3.0_r8*FE_r(i,j)-4.0_r8*FE_r(i,j+1)+FE_r(i,j+2))**2


            alpha1m=0.1_r8/(S1m+eps)**2
            alpha2m=0.6_r8/(S2m+eps)**2
            alpha3m=0.3_r8/(S3m+eps)**2
!
            alpha1p=0.3_r8/(S1p+eps)**2
            alpha2p=0.6_r8/(S2p+eps)**2
            alpha3p=0.1_r8/(S3p+eps)**2
!
            alpham=alpha1m+alpha2m+alpha3m
            alphap=alpha1p+alpha2p+alpha3p
!
            w1m=alpha1m/alpham
            w2m=alpha2m/alpham
            w3m=alpha3m/alpham
            w1p=alpha1p/alphap
            w2p=alpha2p/alphap
            w3p=alpha3p/alphap
!
            q1m=oneothree*FE_r(i,j-3)-sevenosix*FE_r(i,j-2)+            &
     &          elevenosix*FE_r(i,j-1)
            q2m=-oneosix*FE_r(i,j-2)+fiveosix*FE_r(i,j-1)+              &
     &          oneothree*FE_r(i,j)
            q3m=oneothree*FE_r(i,j-1)+fiveosix*FE_r(i,j)-               &
     &          oneosix*FE_r(i,j+1)
!
            q1p=-oneosix*FE_r(i,j-2)+fiveosix*FE_r(i,j-1)+              &
     &          oneothree*FE_r(i,j)
            q2p=oneothree*FE_r(i,j-1)+fiveosix*FE_r(i,j)-               &
     &          oneosix*FE_r(i,j+1)
            q3p=elevenosix*FE_r(i,j)-sevenosix*FE_r(i,j+1)+             &
     &          oneothree*FE_r(i,j+2)
!
!           signa=(FE_r(i,j)-FE_r(i,j-1))*(h(i,j)-h(i,j-1))
            signa=FE_r(i,j)
            cff=SIGN(1.0_r8,signa)
            FEm=0.5_r8*(1.0_r8+cff)*(w1m*q1m+w2m*q2m+w3m*q3m)
            FEp=0.5_r8*(1.0_r8-cff)*(w1p*q1p+w2p*q2p+w3p*q3p)
!
            FE(i,j)=FEm+FEp

# else
            FE(i,j)=0.5_r8*(FE_r(i,j-1)+FE_r(i,j))
# endif
# ifdef SLOPE_KIRWAN
!           cff1=30.0_r8
            cff1=10.0_r8
            dzdy=(h(i,j)-h(i  ,j-1))/on_v(i,j)
            a_slopey=(MAX(0.0_r8,abs(dzdy)-0.05_r8)                     &
     &               *SIGN(1.0_r8,dzdy)*cff1)                           &
     &               *on_r(i,j)*dt(ng)
#  ifdef SED_MORPH
            a_slopey=a_slopey*morph_fac(ised,ng)
#  endif
            FE(i,j)=FE(i,j)+a_slopey
# endif
# ifdef MASKING
            FE(i,j)=FE(i,j)*vmask(i,j)
#  ifdef WET_DRY
            FE(i,j)=FE(i,j)*vmask_wet(i,j)
#  endif
# endif
          END DO
        END DO
# ifdef SED_SLUMP
!
!  Sed slump computation to allow slumping for wet areas everywhere
!  and  at the wet/dry interface.
!
!  sedslopes are the critical slopes to allow slumping.
!  slopefac are the scale factors for sediment movement.
!
!  U-direction slumping
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff2=Srho(ised,ng)*(1.0_r8-bed(i,j,1,iporo))
            slopefac_dry_local=cff2*dt(ng)*slopefac_dry(ng)
            slopefac_wet_local=cff2*dt(ng)*slopefac_wet(ng)
#  ifdef SED_MORPH
            slopefac_wet_local=slopefac_wet_local*morph_fac(ised,ng)
            slopefac_dry_local=slopefac_dry_local*morph_fac(ised,ng)
#  endif
            dzdx=(h(i,j)-h(i-1,j))/om_u(i,j)
            dzdy=(h(i,j)-h(i,j-1))/on_v(i,j)
            dzdxdy=sqrt(dzdx**2.0_r8+dzdy**2.0_r8)
!  For the wet part
            cff=dzdxdy-sedslope_crit_wet(ng)
            cff1=(0.5_r8+SIGN(0.5_r8,cff))
!jcw
!           cff=0.5_r8*cff*cff1/(pm(i,j)*pn(i,j))
!           cff2=cff*slopefac_wet_local*SIGN(1.0_r8,dzdx)
!mai
            cff=ABS(dzdxdy-sedslope_crit_wet(ng))*dzdx/(dzdxdy+eps)
            cff2=cff*cff1/(pm(i,j)*pn(i,j))*slopefac_wet_local
#  ifdef MASKING
             cff2=cff2*umask(i,j)
#   ifdef WET_DRY
            cff2=cff2*umask_wet(i,j)
#   endif
#  endif
            FX(i,j)=FX(i,j)+cff2
!  For the dry part
!jcw
!           cff=ABS(dzdx)-sedslope_crit_dry(ng)
!           cff1=(0.5_r8+SIGN(0.5_r8,cff))
!           cff=0.5_r8*cff*cff1/(pm(i,j)*pn(i,j))
!           cff2=cff*slopefac_dry_local*SIGN(1.0_r8,dzdx)
!mai
            cff=dzdxdy-sedslope_crit_dry(ng)
            cff1=(0.5_r8+SIGN(0.5_r8,cff))
            cff=ABS(dzdxdy-sedslope_crit_dry(ng))*dzdx/(dzdxdy+eps)
            cff2=cff*cff1/(pm(i,j)*pn(i,j))*slopefac_dry_local
#  ifdef MASKING
            cff2=cff2*umask(i,j)
#   ifdef WET_DRY
            ii=MAX(i-1,1)
            ip=MIN(i+1,Lm(ng)+1)
            cff2=cff2*(1.0_r8-umask_wet(i,j))*                          &
     &                ((1.0_r8-umask_wet(ii,j))*umask_wet(ip,j)+        &
     &                 (1.0_r8-umask_wet(ip,j))*umask_wet(ii,j))
#   endif
#  endif
            FX(i,j)=FX(i,j)+cff2
          END DO
        END DO
!
!  V-direction slumping
!
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff2=Srho(ised,ng)*(1.0_r8-bed(i,j,1,iporo))
            slopefac_dry_local=cff2*dt(ng)*slopefac_dry(ng)
            slopefac_wet_local=cff2*dt(ng)*slopefac_wet(ng)
#  ifdef SED_MORPH
            slopefac_wet_local=slopefac_wet_local*morph_fac(ised,ng)
            slopefac_dry_local=slopefac_dry_local*morph_fac(ised,ng)
#  endif
            dzdy=(h(i,j)-h(i,j-1))/on_v(i,j)   ! positive is downhill
            dzdx=(h(i,j)-h(i-1,j))/om_u(i,j)
            dzdxdy=sqrt(dzdx**2.0_r8+dzdy**2.0_r8)
! For the wet part
            cff=ABS(dzdxdy)-sedslope_crit_wet(ng)
            cff1=(0.5_r8+SIGN(0.5_r8,cff))
!jcw
!           cff=0.5_r8*cff*cff1/(pm(i,j)*pn(i,j))
!           cff2=cff*slopefac_wet_local*SIGN(1.0_r8,dzdy)
!mai
            cff=ABS(dzdxdy-sedslope_crit_wet(ng))*dzdy/(dzdxdy+eps)
            cff2=cff*cff1/(pm(i,j)*pn(i,j))*slopefac_wet_local
#  ifdef MASKING
            cff2=cff2*vmask(i,j)
#   ifdef WET_DRY
            cff2=cff2*vmask_wet(i,j)
#   endif
#  endif
            FE(i,j)=FE(i,j)+cff2
! For the dry part
!jcw
!           cff=ABS(dzdy)-sedslope_crit_dry(ng)
!           cff1=(0.5_r8+SIGN(0.5_r8,cff))
!           cff=0.5_r8*cff*cff1/(pm(i,j)*pn(i,j))
!           cff2=cff*slopefac_dry_local*SIGN(1.0_r8,dzdy)
!mai
            cff=dzdxdy-sedslope_crit_dry(ng)
            cff1=(0.5_r8+SIGN(0.5_r8,cff))
            cff=ABS(dzdxdy-sedslope_crit_dry(ng))*dzdy/(dzdxdy+eps)
            cff2=cff*cff1/(pm(i,j)*pn(i,j))*slopefac_dry_local
#  ifdef MASKING
            cff2=cff2*vmask(i,j)
#   ifdef WET_DRY
            jj=MAX(j-1,1)
            jp=MIN(j+1,Mm(ng)+1)
            cff2=cff2*(1.0_r8-vmask_wet(i,j))*                          &
     &                ((1.0_r8-vmask_wet(i,jj))*vmask_wet(i,jp)+        &
     &                 (1.0_r8-vmask_wet(i,jp))*vmask_wet(i,jj))
#   endif
#  endif
            FE(i,j)=FE(i,j)+cff2
          END DO
        END DO
# endif
!
! Limit fluxes to prevent bottom from breaking thru water surface.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
!
! Compute Total thickness available and change.
!
            IF (FX(i,j).ge.0.0_r8) THEN
              Dstp=z_w(i,j,1)-z_w(i,j,0)
!             Dstp=z_w(i,j,N(ng))-z_w(i,j,0)
              rhs_bed=FX(i,j)*pm(i,j)*pn(i,j)
              bed_change=rhs_bed/(Srho(ised,ng)*                        &
     &                   (1.0_r8-bed(i,j,1,iporo)))
            ELSE
              Dstp=z_w(i-1,j,1)-z_w(i-1,j,0)
!             Dstp=z_w(i-1,j,N(ng))-z_w(i-1,j,0)
              rhs_bed=ABS(FX(i,j))*pm(i-1,j)*pn(i-1,j)
              bed_change=rhs_bed/(Srho(ised,ng)*                        &
     &                   (1.0_r8-bed(i-1,j,1,iporo)))
            END IF
!
! Limit that change to be less than available.
!
!           cff=MAX(bed_change-0.75_r8*Dstp,0.0_r8)
            cff=MAX(bed_change-1.00_r8*Dstp,0.0_r8)
            cff1=cff/ABS(bed_change+eps)
            FX(i,j)=FX(i,j)*(1.0_r8-cff1)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
!
! Compute Total thickness available and change.
!
            IF (FE(i,j).ge.0.0_r8) THEN
              Dstp=z_w(i,j,1)-z_w(i,j,0)
!             Dstp=z_w(i,j,N(ng))-z_w(i,j,0)
              rhs_bed=FE(i,j)*pm(i,j)*pn(i,j)
              bed_change=rhs_bed/(Srho(ised,ng)*                        &
     &                   (1.0_r8-bed(i,j,1,iporo)))
            ELSE
              Dstp=z_w(i,j-1,1)-z_w(i,j-1,0)
!             Dstp=z_w(i,j-1,N(ng))-z_w(i,j-1,0)
              rhs_bed=ABS(FE(i,j))*pm(i,j-1)*pn(i,j-1)
              bed_change=rhs_bed/(Srho(ised,ng)*                        &
     &                   (1.0_r8-bed(i,j-1,1,iporo)))
            END IF
!
! Limit that change to be less than available.
!
!           cff=MAX(bed_change-0.75_r8*Dstp,0.0_r8)
            cff=MAX(bed_change-1.00_r8*Dstp,0.0_r8)
            cff1=cff/ABS(bed_change+eps)
            FE(i,j)=FE(i,j)*(1.0_r8-cff1)
          END DO
        END DO
!
!  Apply boundary conditions (gradient).
!
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,IendR
              FX(i,Jstr-1)=FX(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,IendR
              FX(i,Jend+1)=FX(i,Jend)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (LBC(iwest,isTvar(idsed(ised)),ng)%closed) THEN
              DO j=Jstr-1,Jend+1
                FX(Istr,j)=0.0_r8
              END DO
            END IF
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (LBC(ieast,isTvar(idsed(ised)),ng)%closed) THEN
              DO j=Jstr-1,Jend+1
                FX(Iend+1,j)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,JendR
              FE(Istr-1,j)=FE(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,JendR
              FE(Iend+1,j)=FE(Iend,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (LBC(isouth,isTvar(idsed(ised)),ng)%closed) THEN
              DO i=Istr-1,Iend+1
                FE(i,Jstr)=0.0_r8
              END DO
            END IF
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (LBC(inorth,isTvar(idsed(ised)),ng)%closed) THEN
              DO i=Istr-1,Iend+1
                FE(i,Jend+1)=0.0_r8
              END DO
            END IF
          END IF
        END IF
!
!  Determine flux divergence and evaluate change in bed properties.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=(FX(i+1,j)-FX(i,j)+                                     &
     &           FE(i,j+1)-FE(i,j))*pm(i,j)*pn(i,j)
            bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,nstp,ised)-    &
     &                                    cff,0.0_r8)
# if !defined SUSPLOAD
            DO k=2,Nbed
              bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            END DO
# endif
            bed(i,j,1,ithck)=MAX(bed(i,j,1,ithck)-                      &
     &                           cff/(Srho(ised,ng)*                    &
     &                                (1.0_r8-bed(i,j,1,iporo))),       &
     &                           0.0_r8)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Output bedload fluxes.
!-----------------------------------------------------------------------
!
        cff=0.5_r8/dt(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            bedldu(i,j,ised)=FX(i,j)*(pn(i-1,j)+pn(i,j))*cff
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            bedldv(i,j,ised)=FE(i,j)*(pm(i,j-1)+pm(i,j))*cff
          END DO
        END DO
      END DO
!
!  Need to update bed mass for the non-cohesive sediment types, becasue 
!  they did not partake in the bedload transport.
!
      DO ised=1,NCS
        DO j=Jstr,Jend
          DO i=Istr,Iend
            bed_mass(i,j,1,nnew,ised)=bed_mass(i,j,1,nstp,ised)
# if !defined SUSPLOAD
            DO k=2,Nbed
              bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            END DO
# endif
          END DO
        END DO
      END DO
!
!  Update mean surface properties.
!  Sd50 must be positive definite, due to BBL routines.
!  Srho must be >1000, due to (s-1) in BBL routines.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          cff3=0.0_r8
          DO ised=1,NST
            cff3=cff3+bed_mass(i,j,1,nnew,ised)
          END DO
          DO ised=1,NST
            bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/MAX(cff3,eps)
          END DO
!
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,NST
            cff1=cff1*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
            cff2=cff2*Sd50(ised,ng)**bed_frac(i,j,1,ised)
            cff3=cff3*(wsed(ised,ng)+eps)**bed_frac(i,j,1,ised)
            cff4=cff4*Srho(ised,ng)**bed_frac(i,j,1,ised)
          END DO
          bottom(i,j,itauc)=cff1
          bottom(i,j,isd50)=MIN(cff2,Zob(ng))
          bottom(i,j,iwsed)=cff3
          bottom(i,j,idens)=MAX(cff4,1050.0_r8)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
      DO ised=1,NST
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_frac(:,:,:,ised))
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_mass(:,:,:,nnew,ised))
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            bedldu(:,:,ised))
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            bedldv(:,:,ised))
        END IF
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed_frac,                                     &
     &                    bed_mass(:,:,:,nnew,:))
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL mp_exchange3d (ng, tile, iNLM, 2,                          &
     &                      LBi, UBi, LBj, UBj, 1, NST,                 &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      bedldu, bedldv)
      END IF
# endif
      DO i=1,MBEDP
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed(:,:,:,i))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, MBEDP,        &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed)
# endif
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, MBOTP,                   &
     &                  bottom)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, MBOTP,                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom)
# endif
!
      RETURN
      END SUBROUTINE sed_bedload_tile
#endif
      END MODULE sed_bedload_mod
