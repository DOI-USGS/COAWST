#include "cppdefs.h"
#define SED_WENO
#undef  SED_UPWIND
      MODULE sed_bedload_vandera_mod

#if defined NONLINEAR && defined SEDIMENT && defined BEDLOAD_VANDERA  
!
!svn $Id: sed_bedload_vandera.F 429 2009-12-20 17:30:26Z arango $
!======================================================================!
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!----------------------------------------------Tarandeep S. Kalra-------
!------------------------------------------------Chris Sherwood --------
!----------------------------------------------- John C. Warner---------
!-----------------------------------------------------------------------
!  This routine computes sediment bedload transport using              !
!  Van der A et al.(2013) formulation  for unidirectional flow and     !
!  accounts for wave asymmetry leading to differential sediment        !
!  transport for crest and trough cycles.                              ! 
!                                                                      !
!  References:                                                         !
!
!  van der A, D.A., Ribberink, J.S., van der Werf, J.J.,O'Donoghue, T.,!
!  Buijsrogge, R.H., Kranenburg, W.M., (2013). Practical sand transport!
!  formula for non-breaking waves and currents. Coastal Engineering,   !
!  76, pp.26-42
!
!  Kalra, T.S., Suttles, S., Sherwood, C.R., Warner, J.C.              !
!  Aretxabaleta,A.L., and Leavitt, G.R. (2022). Shoaling Wave Shape    !
!  Estimates from Field Observations and derived Bedload Sediment Rates!
!  J. Mar.Sci.Eng. 2022, 10, 223. https://doi.org/10.3390/jmse10020223 !
!
!  Udated sed bed evolution scheme to the WENO method of:              !
!  Wen Long, James T. Kirby, Zhiyu Shao,                               !
!  A numerical scheme for morphological bed level calculations,        !
!  Coastal Engineering,55, Issue 2, 2008, 167-180.                     !
!  https://doi.org/10.1016/j.coastaleng.2007.09.009.                   !
!                                                                      !
!----------------------------------------------------------------------!
!======================================================================!
!
      implicit none

      PRIVATE
      PUBLIC  :: sed_bedload_vandera

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_bedload_vandera (ng, tile)
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL sed_bedload_vandera_tile (ng, tile,                          &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp(ng), nnew(ng),                        &
     &                       BBL(ng)    % bustrc,                       &
     &                       BBL(ng)    % bvstrc,                       &
     &                       FORCES(ng) % Hwave,                        &
     &                       FORCES(ng) % Lwave,                        &
     &                       FORCES(ng) % Dwave,                        &
     &                       FORCES(ng) % Pwave_bot,                    &
     &                       FORCES(ng) % Uwave_rms,                    &
     &                       GRID(ng) % angler,                         &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % om_r,                           &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % on_r,                           &
     &                       GRID(ng) % on_v,                           &
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
     &                       OCEAN(ng) % zeta,                          &
     &                       SEDBED(ng) % ursell_no,                    &
     &                       SEDBED(ng) % RR_asymwave,                  &
     &                       SEDBED(ng) % beta_asymwave,                &
     &                       SEDBED(ng) % ucrest_r,                     &
     &                       SEDBED(ng) % utrough_r,                    &
     &                       SEDBED(ng) % T_crest,                      &
     &                       SEDBED(ng) % T_trough,                     &
     &                       SEDBED(ng) % bedldu,                       &
     &                       SEDBED(ng) % bedldv,                       &
     &                       SEDBED(ng) % bed,                          &
     &                       SEDBED(ng) % bed_frac,                     &
     &                       SEDBED(ng) % bed_mass,                     &
     &                       SEDBED(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sed_bedload_vandera
!
!***********************************************************************
      SUBROUTINE sed_bedload_vandera_tile (ng, tile,                    &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew,                          &
     &                             bustrc, bvstrc,                      &
     &                             Hwave, Lwave,                        &
     &                             Dwave, Pwave_bot,                    &
     &                             Uwave_rms,                           &
     &                             angler,                              &
     &                             h, om_r, om_u, on_r, on_v,           &
     &                             pm, pn,                              &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef WET_DRY
     &                             rmask_wet, umask_wet, vmask_wet,     &
# endif
     &                             z_w,                                 &
     &                             zeta,                                &
     &                             ursell_no,                           &
     &                             RR_asymwave, beta_asymwave,          &
     &                             ucrest_r, utrough_r,                 &
     &                             T_crest, T_trough,                   &
     &                             bedldu, bedldv,                      &
     &                             bed, bed_frac, bed_mass,             &
     &                             bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_sediment
      USE mod_vandera_funcs
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE bc_3d_mod, ONLY : bc_r3d_tile
      USE exchange_2d_mod, ONLY : exchange_u2d_tile, exchange_v2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d,         &
     &                            mp_exchange4d
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
      real(r8), intent(in) :: bustrc(LBi:,LBj:)
      real(r8), intent(in) :: bvstrc(LBi:,LBj:)
!
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
      real(r8), intent(in) :: Uwave_rms(LBi:,LBj:)
!
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
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
!
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
!
      real(r8), intent(inout) :: ursell_no(LBi:,LBj:)
      real(r8), intent(inout) :: RR_asymwave(LBi:,LBj:)
      real(r8), intent(inout) :: beta_asymwave(LBi:,LBj:)
      real(r8), intent(inout) :: ucrest_r(LBi:,LBj:)
      real(r8), intent(inout) :: utrough_r(LBi:,LBj:)
      real(r8), intent(inout) :: T_crest(LBi:,LBj:)
      real(r8), intent(inout) :: T_trough(LBi:,LBj:)
!
      real(r8), intent(inout) :: bedldu(LBi:,LBj:,:)
      real(r8), intent(inout) :: bedldv(LBi:,LBj:,:)
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: bustrc(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrc(LBi:UBi,LBj:UBj)
!
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pwave_bot(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uwave_rms(LBi:UBi,LBj:UBj)
!
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
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
!
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
!
      real(r8), intent(inout) :: ursell_no(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: RR_asymwave(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: beta_asymwave(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ucrest_r((LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: utrough_r(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: T_crest(LBi:UBi,LBj:UBj) 
      real(r8), intent(inout) :: T_trough(LBi:UBi,LBj:UBj)
!
      real(r8), intent(inout) :: bedldu(LBi:UBi,LBj:UBj,NST)
      real(r8), intent(inout) :: bedldv(LBi:UBi,LBj:UBj,NST)
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
# endif
!
!  Local variable declarations.
!
      integer :: i, ii, ip, ised, j, jj, jp, k
!
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, fac1, fac2
      real(r8) :: Dstp, bed_change, dz, roll
      real(r8) :: a_slopex, a_slopey, sed_angle
      real(r8) :: bedld, bedld_mass, dzdx, dzdy, dzdxdy
      real(r8) :: rhs_bed, Ua, Ra, Clim, phi_cw
!
      real(r8) :: Hs, Td, depth
      real(r8) :: d50, d50_mix, d90, rhos
      real(r8) :: urms, umag_curr, phi_curwave, udelta
      real(r8) :: y, uhat, ahat
      real(r8) :: k_wn, c_w
      real(r8) :: smgd, osmgd
!
      real(r8) :: r, phi, Su, Au
      real(r8) :: Sk, Ak
      real(r8) :: T_cu, T_tu
      real(r8) :: umax, umin
!
      real(r8) :: uhat_c, uhat_t
      real(r8) :: mag_uc, mag_ut
!
      real(r8) :: theta
      real(r8) :: fd, ksw, eta, alpha, tau_wRe
      real(r8) :: dsf_c, dsf_t
      real(r8) :: theta_c, theta_t
      real(r8) :: theta_cx, theta_cy, theta_tx, theta_ty
      real(r8) :: mag_theta_c, mag_theta_t
      real(r8) :: mag_bstrc
      real(r8) :: om_cc, om_tt, om_ct, om_tc
      real(r8) :: smgd_3
!
      real(r8) :: bedld_cx, bedld_cy
      real(r8) :: bedld_tx, bedld_ty
      real(r8) :: bedld_x, bedld_y
!
      real(r8) :: wavecycle, alphac, alphaw
      real(r8) :: twopi, otwopi, sqrt2
!
# ifdef SED_SLUMP
      real(r8) :: slopefac_wet_local, slopefac_dry_local
# endif
      real(r8), parameter :: eps = 1.0E-14_r8

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
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phic
!
!  Need local arrays for the global vars because we fill the local
!  arrays larger than the standard stencil.
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ursell_nol
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: RR_asymwavel
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: beta_asymwavel
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ksd_wbll
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: fd_wbll
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phi_wcl
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ucrest_rl
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: utrough_rl
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: T_crestl
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: T_troughl

# include "set_bounds.h"
!
      twopi=2.0_r8*pi
      otwopi=1.0_r8/twopi
      sqrt2=SQRT(2.0_r8)
      FX_r=0.0_r8
      FE_r=0.0_r8
      fd=0.0_r8
!
!-----------------------------------------------------------------------
!  Compute bedload sediment transport.
!-----------------------------------------------------------------------
!
      sed_angle=DTAN(33.0_r8*pi/180.0_r8)
      alphac=bedload_vandera_alphac(ng)
      alphaw=bedload_vandera_alphaw(ng)
!
# if defined SED_WENO
      DO j=Jstrm3,Jendp2i
        DO i=Istrm3,Iendp2i
# else
      DO j=Jstrm1,Jendp1
        DO i=Istrm1,Iendp1
# endif
!
! Compute angle between currents and waves, measure CCW from current
! direction toward wave vector.
!
          phi_cw=1.5_r8*pi-Dwave(i,j)-bottom(i,j,idpcx)-angler(i,j)
!
! Compute angle between waves and current, measure CCW from wave 
! towards current vector
!
          phi_wcl(i,j)=2.0_r8*pi-phi_cw
!
        END DO
      END DO
!
      DO ised=NCS+1,NST
        rhos=Srho(ised,ng)                      ! (kg/m3)
        d50=sd50(ised,ng)                       ! (m)
        d90=1.3_r8*d50                          ! (m)
        IF(NST>1) THEN 
          d50_mix=0.0003                          ! 0.3 mm 
        ELSE
          d50_mix=d50
        ENDIF 
!
        cff=rhos/rho0
        smgd=(cff-1.0_r8)*g*d50
        osmgd=1.0_r8/smgd
!
        smgd_3=sqrt((cff-1.0_r8)*g*d50**3.0_r8)
!
# if defined SED_WENO
      DO j=Jstrm3,Jendp2i
        DO i=Istrm3,Iendp2i
# else
        DO j=Jstrm1,Jendp1
          DO i=Istrm1,Iendp1
# endif
!
            Hs=Hwave(i,j)                             ! (m)
            depth=MAX(h(i,j)+zeta(i,j,1),Dcrit(ng))   ! (m)
            Td=MAX(Pwave_bot(i,j),1.0_r8)             ! (s)
            urms=Uwave_rms(i,j)                       ! (m/s)
            phi_curwave=phi_wcl(i,j)
            udelta=bottom(i,j,idubl)
# if defined BEDLOAD_VANDERA_MADSEN_UDELTA
            ksd_wbll(i,j)=bottom(i,j,idksd)
# endif
! 
! Compute magnitude of stress for computing current velocity 
! at the wave boundary layer
!
            mag_bstrc=SQRT(bustrc(i,j)*bustrc(i,j)+                     &
     &                     bvstrc(i,j)*bvstrc(i,j))
!
            uhat=urms*sqrt2
            ahat=uhat*Td*otwopi
            k_wn=kh(Td,depth)/depth                ! Wave number 
            c_w=2.0_r8*pi/(k_wn*Td)                ! Wave speed
!
! VA-2013 equation 1 is solved in 3 sub-steps
!
!----------------------------------------------------------------------
! Ruessink et al. provides equations for calculating skewness parameters
! Uses Malarkey and Davies equations to get "r" and "phi"
! common to both crest and trough.
!-----------------------------------------------------------------------
!
            CALL skewness_params(Hs, Td, depth, r, phi, ursell_nol(i,j))
!        
!-----------------------------------------------------------------------
! Abreu et al. use skewness params to get representative critical orbital
! velocity for crest and trough cycles , use r and phi.
!-----------------------------------------------------------------------
! 
            CALL abreu_points(r, phi, uhat, Td,                         &
     &                        T_crestl(i,j), T_troughl(i,j),            &
     &                        T_cu, T_tu, umax, umin,                   &
     &                        RR_asymwavel(i,j), beta_asymwavel(i,j))
!
!-----------------------------------------------------------------------
!                         Crest half cycle
!-----------------------------------------------------------------------
! Step 1. Representative crest half cycle water particle velocity
! as well as full cycle orbital velocity and excursion.
!-----------------------------------------------------------------------
!
            uhat_c=umax
            uhat_t=umin
!
!-----------------------------------------------------------------------
! VA2013 Equation 10, 11.
!-----------------------------------------------------------------------
!
            ucrest_rl(i,j)=0.5_r8*sqrt2*uhat_c
            utrough_rl(i,j)=0.5_r8*sqrt2*uhat_t
!
            smgd=(rhos/rho0-1.0_r8)*g*d50
            osmgd=1.0_r8/smgd
!
! Full wave cycle
!
            CALL full_wave_cycle_stress_factors(ng, d50, d90, osmgd,    &
     &                                                 Td, depth,       &
     &                                    umag_curr, phi_curwave,       &
     &                                    RR_asymwavel(i,j), uhat, ahat,&
     &                                                umax, umin,       &
     &                                                 mag_bstrc,       &
     &                                           alphac, alphaw,        &
     &                                   T_crestl(i,j), T_troughl(i,j), &
     &                                                T_cu, T_tu,       &
     &                                         ksd_wbll(i,j),  udelta,  &
     &                                                    fd_wbll(i,j), &
                                           alpha, eta, ksw, tau_wRe )
!
!-----------------------------------------------------------------------
! 2. Bed shear stress (Shields parameter) for Crest half cycle 
!    alpha VA2013 Eqn. 19  
!-----------------------------------------------------------------------
!
!    alpha VA2013 Eqn. 19  
!-----------------------------------------------------------------------
!
            CALL half_wave_cycle_stress_factors( T_cu, T_crestl(i,j),   &
     &                                           ahat, ksw,             &
     &                                          fd_wbll(i,j), alpha,    &
     &                                           alphac, alphaw,        &
     &                                           d50, osmgd,            &
     &            ucrest_rl(i,j), uhat_c, udelta, phi_curwave,          &
     &                                           tau_wRe,               &
     &                          dsf_c, theta_cx, theta_cy, mag_theta_c )
!
!-----------------------------------------------------------------------
! 3. Compute sediment load entrained during each crest half cycle
!-----------------------------------------------------------------------
!
            wavecycle=1.0_r8
            CALL sandload_vandera( wavecycle,                           &
     &                              Hs, Td,  depth, RR_asymwavel(i,j),  &
     &                              d50, d50_mix, rhos, c_w,            &
     &                              eta, dsf_c,                         &
     &                              T_crestl(i,j), T_cu, uhat_c,        &
     &                              mag_theta_c, om_cc, om_ct )
!
!-----------------------------------------------------------------------
! 2. Bed shear stress (Shields parameter) for Trough half cycle 
!    alpha VA2013 Eqn. 19  
!-----------------------------------------------------------------------
!
            CALL half_wave_cycle_stress_factors( T_tu, T_troughl(i,j),  &
     &                                           ahat, ksw,             &
     &                                          fd_wbll(i,j), alpha,    &
     &                                           alphac, alphaw,        &
     &                                           d50, osmgd,            &
     &           utrough_rl(i,j), uhat_t, udelta, phi_curwave,          &
     &                                           tau_wRe,               &
     &                          dsf_t, theta_tx, theta_ty, mag_theta_t )
!
!-----------------------------------------------------------------------
! 3. Compute sediment load entrained during each trough half cycle
!-----------------------------------------------------------------------
!
            wavecycle=-1.0_r8
            CALL sandload_vandera( wavecycle,                           &
     &                              Hs, Td,  depth, RR_asymwavel(i,j),  &
     &                              d50, d50_mix, rhos, c_w,            &
     &                              eta, dsf_t,                         &
     &                              T_troughl(i,j), T_tu, uhat_t,       &
     &                              mag_theta_t, om_tt, om_tc )
!
!-----------------------------------------------------------------------
! 3. Compute sediment load entrained during each trough half cycle
!-----------------------------------------------------------------------
!
            cff1=MAX(0.5_r8*T_crestl(i,j)/(T_cu+eps), 0.0_r8)
!
            cff2=sqrt(mag_theta_c)*(om_cc+cff1*om_tc)  
            cff3=(theta_cx/mag_theta_c)
            bedld_cx=cff2*cff3
!
            cff3=(theta_cy/mag_theta_c)
            bedld_cy=cff2*cff3
!
            cff1=MAX(0.5_r8*T_troughl(i,j)/(T_tu+eps), 0.0_r8)
!
            cff2=sqrt(mag_theta_t)*(om_tt+cff1*om_ct)
            cff3=(theta_tx/mag_theta_t)
            bedld_tx=cff2*cff3
!
            cff3=(theta_ty/mag_theta_t)
            bedld_ty=cff2*cff3
!
!-----------------------------------------------------------------------
! VA2013  Use the velocity-load equation 1. 
! Units of smgd_3 are m2-s-1
!-----------------------------------------------------------------------
!
            smgd_3=sqrt((rhos/rho0-1.0_r8)*g*d50**3.0_r8)
!
            bedld_x=smgd_3*( bedld_cx*T_crestl(i,j)+                     &
     &                       bedld_tx*T_troughl(i,j) )/(Td+eps)
            bedld_y=smgd_3*( bedld_cy*T_crestl(i,j)+                     &
     &                       bedld_ty*T_troughl(i,j) )/(Td+eps)
!
! The units of these are kg m-1 sec-1
! COMMENTED FOR NOW 
!
            bedld_x=rhos*bedld_x*bed_frac(i,j,1,ised)
            bedld_y=rhos*bedld_y*bed_frac(i,j,1,ised)
!           
! Convert bedload from the wave aligned axis to xi and eta directions
! 
            theta=1.5_r8*pi-Dwave(i,j)-angler(i,j)
!
! Partition bedld into xi and eta directions, still at rho points.
! (FX_r and FE_r have dimensions of kg).
!
            FX_r(i,j)=(bedld_x*COS(theta)-bedld_y*SIN(theta))*          &
     &                on_r(i,j)*dt(ng)
            FE_r(i,j)=(bedld_x*SIN(theta)+bedld_y*COS(theta))*          &
     &                om_r(i,j)*dt(ng)
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
! Add contribution of bed slope to bed load transport fluxes.
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
! Add contribution of bed slope to bed load transport fluxes.
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
     &               cff2*(FX_r(i-1,j)+FX_r(i,j)))+                     &
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
     &               cff2*(FE_r(i,j-1)+FE_r(i,j)))+                     &
     &              0.5_r8*(1.0_r8-SIGN(1.0_r8,FE_r(i,j-1)))*           &
     &              (cff2*FE_r(i  ,j)+                                  &
     &               cff1*0.5_r8*(FE_r(i,j-1)+FE_r(i,j)))
# elif defined SED_WENO
!
!  Long et al. (2008). Coastal Engr, 55, 167-180.
!
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
#  ifdef WET_DRY
            cff2=cff2*umask_wet(i,j)
#   endif
#  endif
            FX(i,j)=FX(i,j)+cff2*bed_frac(i,j,1,ised)
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
            FX(i,j)=FX(i,j)+cff2*bed_frac(i,j,1,ised)
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
#   ifdef SED_MORPH
            slopefac_wet_local=slopefac_wet_local*morph_fac(ised,ng)
            slopefac_dry_local=slopefac_dry_local*morph_fac(ised,ng)
#   endif
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
            FE(i,j)=FE(i,j)+cff2*bed_frac(i,j,1,ised)
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
            FE(i,j)=FE(i,j)+cff2*bed_frac(i,j,1,ised)
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
              DO j=JstrR,JendR
                FX(Istr,j)=0.0_r8
              END DO
            END IF
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (LBC(ieast,isTvar(idsed(ised)),ng)%closed) THEN
              DO j=JstrR,JendR
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
              DO i=IstrR,IendR
                FE(i,Jstr)=0.0_r8
              END DO
            END IF
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (LBC(inorth,isTvar(idsed(ised)),ng)%closed) THEN
              DO i=IstrR,IendR
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
     &                                                   0.0_r8)
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
!  Fill global arrays with local work array data.
!  We cant use a global array outside of the tile.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ursell_no(i,j)=ursell_nol(i,j)
          RR_asymwave(i,j)=RR_asymwavel(i,j)
          beta_asymwave(i,j)=beta_asymwavel(i,j)
          ucrest_r(i,j)=ucrest_rl(i,j)
          utrough_r(i,j)=utrough_rl(i,j)
          T_crest(i,j)=T_crestl(i,j)
          T_trough(i,j)=T_troughl(i,j)
# if defined BEDLOAD_VANDERA_DIRECT_UDELTA
          bottom(i,j,idksd)=ksd_wbll(i,j)
# endif
          bottom(i,j,idpwc)=phi_wcl(i,j)
          bottom(i,j,idfdw)=fd_wbll(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   ursell_no)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   RR_asymwave)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   beta_asymwave)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   ucrest_r)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   utrough_r)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   T_crest)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   T_trough)
# if defined BEDLOAD_VANDERA_DIRECT_UDELTA
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,idksd))
# endif
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,idpwc))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,idfdw))
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
!     DO i=1,1 !dont do all MBEDP, we only changed ithck in this routine
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed(:,:,:,ithck))
!     END DO
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed(:,:,:,ithck))
# endif
!      CALL bc_r3d_tile (ng, tile,                                       &
!     &                  LBi, UBi, LBj, UBj, 1, MBOTP,                   &
!     &                  bottom)
!# ifdef DISTRIBUTE
!# if defined BEDLOAD_VANDERA_DIRECT_UDELTA
!      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
!     &                    LBi, UBi, LBj, UBj,                           &
!     &                    NghostPoints,                                 &
!     &                    EWperiodic(ng), NSperiodic(ng),               &
!     &                    bottom(:,:,idksd))
!# endif
!      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!     &                    LBi, UBi, LBj, UBj,                           &
!     &                    NghostPoints,                                 &
!     &                    EWperiodic(ng), NSperiodic(ng),               &
!     &                    bottom(:,:,idpwc),                            &
!     &                    bottom(:,:,idfdw))
!# endif

      RETURN
      END SUBROUTINE sed_bedload_vandera_tile
! 
! Subroutines and functions required for Van der A formulation.
! 
      SUBROUTINE sandload_vandera( wavecycle,                           &
     &                              Hs, Td,  depth, RR,                 &
     &                              d50, d50_mix, rhos, c_w,            &
     &                              eta, dsf,                           &
     &                              T_i, T_iu, uhat_i, mag_theta_i,     &
     &                              om_ii, om_iy )
!
      USE mod_kinds
      USE mod_scalars
      USE mod_vandera_funcs
!
      implicit none
!
      real(r8), intent(in) :: wavecycle
      real(r8), intent(in) :: Hs, Td, depth, RR
      real(r8), intent(in) :: d50, d50_mix, rhos, c_w
      real(r8), intent(in) :: eta, dsf
      real(r8), intent(in) :: T_i, T_iu
      real(r8), intent(in) :: uhat_i, mag_theta_i
      real(r8), intent(out):: om_ii, om_iy
!
! local variables
! 
      real(r8), parameter :: m_fac=11.0_r8, n_fac=1.2_r8
      real(r8), parameter :: alpha_fac=8.2_r8
      real(r8), parameter :: xi=1.7_r8 ! Based on Santoss_core.m
      real(r8), parameter :: eps=1.0E-14_r8
      real(r8) :: eps_eff
      real(r8) :: om_i
      real(r8) :: theta_diff, theta_ieff, theta_cr
      real(r8) :: w_s
      real(r8) :: ws_eta, ws_dsf
      real(r8) :: w_sc_eta, w_sc_dsf
      real(r8) :: cff, cff1_eta, cff1_dsf
      real(r8) :: P
! 
! Find settling velocity based on Soulsby (1997). 
! VA2013 Use 0.8*d50 for settling velocity (text under equation 28).
!
      w_s=w_s_calc(0.8_r8*d50, rhos)
!
! VA2013 Equation 29, for crest cycle
!
      ws_eta=w_sc_calc(Hs, Td, depth, RR, w_s, eta)
      ws_dsf=w_sc_calc(Hs, Td, depth, RR, w_s, dsf)
      IF(wavecycle.eq.1.0_r8) THEN
        w_sc_eta=MAX(w_s+ws_eta,0.0_r8)
        w_sc_dsf=MAX(w_s+ws_dsf,0.0_r8)
      ENDIF
!
! VA2013 Equation 30, for trough cycle
!
      IF(wavecycle.eq.-1.0_r8) THEN
!        w_sc_eta=(w_s-ws_eta)
!        w_sc_dsf=(w_s-ws_dsf)
        w_sc_eta=MAX(w_s-ws_eta,0.36*w_s)
        w_sc_dsf=MAX(w_s-ws_dsf,0.36*w_s)
!        w_sc_eta=MIN(w_s-ws_eta,0.0_r8)
!        w_sc_dsf=MIN(w_s-ws_dsf,0.0_r8)
      ENDIF
!
! VA2013 Equation 33, Phase lag parameter
!
      cff=1.0_r8-(wavecycle*xi*uhat_i/c_w)
!
      IF( (T_i-T_iu).eq.0.0_r8 ) THEN 
        cff1_eta=0.0_r8
        cff1_dsf=0.0_r8
      ELSE
        cff1_eta=(1.0_r8/(2.0_r8*(T_i-T_iu)*w_sc_eta))
        cff1_dsf=(1.0_r8/(2.0_r8*(T_i-T_iu)*w_sc_dsf))
      ENDIF 
!
      IF(eta.gt.0.0_r8) THEN
!
! For ripple regime 
!
        P=alpha_fac*eta*cff*cff1_eta
      ELSEIF(eta.eq.0.0_r8)THEN
!
! For sheet flow regime 
!
        P=alpha_fac*dsf*cff*cff1_dsf
      ENDIF
!
      eps_eff=(d50/d50_mix)**0.25_r8
!
! CRS for multiple sed types
!
!      eps_eff=1.0_r8 
      theta_ieff=eps_eff*mag_theta_i
! 
! Find critical Shields parameters based on Soulsby (1997).
!
      theta_cr=theta_cr_calc(d50, rhos)
!
! Sand load entrained in the flow during each half-cycle
!
      theta_diff=MAX((theta_ieff-theta_cr),0.0_r8)
      om_i=m_fac*(theta_diff)**n_fac
!
! VA2013 Equation 23-26, Sandload entrained during half cycle 
 
      IF(P.lt.eps) THEN
 
! This is Taran's addition if there are no waves then phase lag=0.0
! 
        om_ii=1.0_r8 
        om_iy=0.0_r8
      ELSEIF(P.gt.eps.and.P.lt.1.0_r8) THEN
        om_ii=om_i
        om_iy=0.0_r8
      ELSE
        om_ii=om_i/P
        cff=1.0_r8/P
        om_iy=om_i*(1.0_r8-cff)
      ENDIF
!
      RETURN
      END SUBROUTINE sandload_vandera
!
      SUBROUTINE full_wave_cycle_stress_factors( ng, d50, d90, osmgd,   &
     &                                                 Td, depth,       &
     &                                    umag_curr, phi_curwave,       &
     &                                            RR, uhat, ahat,       &
     &                                                umax, umin,       &
     &                                                 mag_bstrc,       &
     &                                           alphac, alphaw,        &
     &                                      T_c, T_t, T_cu, T_tu,       &
     &                                                ksd,              &
     &                                                udelta, fd,       &
     &                                   alpha, eta, ksw, tau_wRe )
!
!**********************************************************************
!  This subroutine returns the following:
!  eta                 : ripple height
!  udelta              : current velocity at the wave boundary layer
!  fd                  : current friction factor
!  tau_wRe             : Wave averaged Reynolds stress
!  T_c, T_t, T_cu, T_tu: Updated time periods in half cycles
!                        based on current velocity
!**********************************************************************
!
      USE mod_grid
      USE mod_kinds
      USE mod_scalars
      USE mod_sediment
      USE mod_sedbed
      USE mod_vandera_funcs
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
! Input the crest or trough half cycle velocity
! d50 -- grain size in meters
! Different for crest and trough half cycles
!
      real(r8), intent(in) :: d50, d90, osmgd
      real(r8), intent(in) :: Td, depth
      real(r8), intent(in) :: umag_curr, phi_curwave
      real(r8), intent(in) :: RR, uhat, ahat
      real(r8), intent(in) :: umax, umin
      real(r8), intent(in) :: mag_bstrc
      real(r8), intent(in) :: alphac, alphaw
      real(r8), intent(inout) :: T_c, T_t, T_cu, T_tu, fd, ksd
      real(r8), intent(in) :: udelta
      real(r8), intent(inout) :: alpha, eta, ksw, tau_wRe
!
!  Local variables
! 
      integer  :: iter
      integer,  parameter :: total_iters=15
      real(r8), parameter :: tol=0.001_r8, von_k=0.41_r8
      real(r8), parameter :: eps=1.0E-14_r8
      real(r8), parameter :: crs_fac=1.0_r8
      real(r8) :: theta_timeavg_old, theta_timeavg, theta_hat_i
      real(r8) :: k_wn
      real(r8) :: psi  ! maximum mobility number
      real(r8) :: rlen ! ripple length
      real(r8) :: omega
!     real(r8) :: ksd
      real(r8) :: fw
      real(r8) :: alpha_w, fwd, c_w
      real(r8) :: ustarw
!
! Iterative solution to obtain current and wave related bed roughness
! VA2013 Apendix A, Shields parameter (Stress) depends on bed roughness
! Bed roughness computed from converged Shields parameter
!
! Maximum mobility number at crest and trough
! For irregular waves, use Rayleigh distributed maximum value
! VA, text under equation Appendix B.4
!
      psi=(1.27_r8*uhat)**2*osmgd
!
! Use Appendix B eqn B.1 and B.2 to get ripple height and length
!
      CALL ripple_dim(psi, d50, eta, rlen)
!
      eta=eta*ahat
      rlen=rlen*ahat
!
      omega=2.0_r8*pi/Td
!
! Initiliaze with theta_timeavg=0 and theta_hat_i=theta_timeavg
!
      theta_timeavg=0.0_r8
      theta_timeavg_old=0.0_r8
      fd=0.0_r8
!
# if defined BEDLOAD_VANDERA_MADSEN_UDELTA
      fd=fd_calc_madsen(udelta, mag_bstrc)
# endif
!
      DO iter=1,total_iters
!
! Calculate wave related bed roughness from VA2013 A.5
!
        ksw=ksw_calc(d50, mu_calc(d50), theta_timeavg, eta, rlen)
!
! Calculate full-cycle wave friction factor VA2013 Appendix Eqn. A.4
!
        fw=fw_calc(ahat, ksw)
!
# if defined BEDLOAD_VANDERA_DIRECT_UDELTA
!
! Calculate current-related bed roughness from VA2013 Appendix A.1
! dont need ksd if fd is directly computed from mag_bustrc
        ksd=ksd_calc(d50, d90, mu_calc(d50), theta_timeavg, eta, rlen)
!
! Calculate friction factor from ksd, udelta, delta
!
        fd=fd_calc_santoss(udelta, sg_zwbl(ng), ksd)
!
# endif
!
! Calculate Time-averaged absolute Shields stress VA2013 Appendix Eq. A.3
!
!        theta_timeavg=osmgd*(0.5_r8*fd*udelta**2.0_r8+                  &
!     &                       0.25_r8*fw*uhat**2.0_r8)
        theta_timeavg=osmgd*(0.5_r8*fd*alphac*udelta**2.0_r8+            &
     &                       0.25_r8*fw*alphaw*uhat**2.0_r8)
!
        IF(ABS(theta_timeavg-theta_timeavg_old).lt.tol) THEN
          EXIT
        ENDIF
        theta_timeavg_old=theta_timeavg
      END DO
!
! Calculate full-cycle current friction factor from VA2013 Eqn. 20
! use the stress from COAWST and corresponding current velocity 
!
      alpha=udelta/(udelta+uhat+eps)
!     fwd=alpha*fd+(1.0-alpha)*fw
      fwd=alpha*fd*alphac+(1.0-alpha)*fw*alphaw
!
      k_wn=kh(Td,depth)/depth          ! Wave number 
      c_w=2.0_r8*pi/(k_wn*Td)          ! Wave speed
      alpha_w=0.424_r8
!
# ifdef BEDLOAD_VANDERA_WAVE_AVGD_STRESS
      tau_wRe=MAX((rho0*fwd*alpha_w*uhat**3.0_r8/(2.0_r8*c_w)),eps)
# else 
      tau_wRe=0.0_r8
# endif
!
! Compute the change in time period based on converged udelta
! (current velocity at wave boundary layer)
!
      CALL current_timeperiod(udelta, phi_curwave, umax, umin, RR,      &
     &                        T_c, T_t, Td)
!
# ifdef BEDLOAD_VANDERA_SURFACE_WAVE
!
! Calculate the effect of surface waves
!
      CALL surface_wave_mod(Td, depth, uhat, T_c, T_cu, T_t, T_tu)
# endif
!
      END SUBROUTINE full_wave_cycle_stress_factors
!
      SUBROUTINE half_wave_cycle_stress_factors( T_iu, T_i, ahat, ksw,  &
     &                                           fd, alpha,             &
     &                                           alphac, alphaw,        &
     &                                           d50, osmgd,            &
     &                                ui_r, uhat_i, udelta, phi_curwave,&
     &                                           tau_wRe,               &
     &                            dsf, theta_ix, theta_iy, mag_theta_i )
!
!**********************************************************************
!  This subroutine returns the following:
!  dsf                 : sheetflow thickness
!  theta_ix, theta_iy  : Shields parameter in x and y dir.
!  mag_theta_i         : Magnitude of Shields parameter for half cycle
!**********************************************************************
!
      USE mod_kinds
      USE mod_scalars
      USE mod_vandera_funcs
!
      implicit none
!
! Input the crest or trough half cycle velocity
! d50 -- grain size in meters
! Different for crest and trough half cycles 
!       
      real(r8), intent(in) :: T_iu, T_i, ahat, ksw
      real(r8), intent(in) :: fd, alpha
      real(r8), intent(in) :: alphac, alphaw
      real(r8), intent(in) :: d50, osmgd
      real(r8), intent(in) :: ui_r, uhat_i, udelta, phi_curwave
      real(r8), intent(in) :: tau_wRe
      real(r8), intent(inout) :: dsf, theta_ix, theta_iy, mag_theta_i
!
!  Local variables
!
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8) :: fw_i, fwd_i
      real(r8) :: alpha_w, fwd, k, c_w
      real(r8) :: theta_hat_i
      real(r8) :: ui_rx, ui_ry, mag_ui
!        
! Wave friction factor for wave and crest half cycle VA2013 Eqn. 21
! 
      fw_i=fwi_calc(T_iu, T_i, ahat, ksw)
!
! Wave current friction factor (Madsen and Grant) VA2013 Eqn. 18
! Different for crest and trough 
!
!     fwd_i=alpha*fd+(1.0_r8-alpha)*fw_i
      fwd_i=alpha*fd*alphac+(1.0_r8-alpha)*fw_i*alphaw
!
! VA2013-Magnitude of Shields parameter Eqn. 17
! 
      theta_hat_i=0.5_r8*fwd_i*uhat_i**2*osmgd
!
! Sheet flow thickness VA2013 Appendix C C.1 
! Update from initial value 
!
      dsf=dsf_calc(d50, theta_hat_i) !this dsf is in m 
!
! Calculated the velocity magnitude based on representative velocities
! equation 12 from Van der A, 2013
!
!-----------------------------------------------------------------------
! Get the representative trough half cycle water particle velocity
!    as well as full cycle orbital velocity and excursion
!-----------------------------------------------------------------------
!
      ui_rx=udelta*COS(phi_curwave)*alphac+ui_r*alphaw
      ui_ry=udelta*SIN(phi_curwave)*alphac
!
! mag_ui is set to a min value to avoid non-zero division
!
      mag_ui=MAX( SQRT(ui_rx*ui_rx+ui_ry*ui_ry), eps )
!
! VA2013-Magnitude of Shields parameter Eqn. 17
! 
      mag_theta_i=MAX(0.5_r8*fwd_i*osmgd*(mag_ui**2), 0.0_r8)
!
!-----------------------------------------------------------------------
! Shields parameter in crest cycle
! rho0 is required for non-dimensionalizing 
!-----------------------------------------------------------------------
!
      theta_ix=ABS(mag_theta_i)*ui_rx/(mag_ui)+tau_wRe*osmgd/rho0
      theta_iy=ABS(mag_theta_i)*ui_ry/(mag_ui)
!
! mag_theta_i is set to a min value to avoid non-zero division
!
      mag_theta_i=MAX( sqrt(theta_ix*theta_ix+theta_iy*theta_iy),eps )
!
!
      END SUBROUTINE half_wave_cycle_stress_factors
!
      SUBROUTINE current_timeperiod(unet, phi_curwave, umax, umin,      &
     &                              RR, T_c, T_t, Td)
!
!**********************************************************************
!  This subroutine returns the following:
!  T_c, T_t  : Time period in crest and trough cycle
!**********************************************************************
!
! Modify the crest and trough time periods based on current velocity
! This function was developed by Chris Sherwood and Tarandeep Kalra
!
! The basis for these formulations are formed from Appendix A.3 
! in SANTOSS report. 
! Report number: SANTOSS_UT_IR3 Date: January 2010
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none
!  
      real(r8), intent(in) :: unet, phi_curwave
      real(r8), intent(in) :: umax, umin
      real(r8), intent(in) :: RR, Td
      real(r8), intent(inout) :: T_c, T_t
!
!  Local variables
!
      real(r8) :: unet_xdir, udelta, delt
!
      unet_xdir=unet*cos(phi_curwave)

      IF(RR.eq.0.5_r8) THEN
        T_c=0.5_r8*Td
        T_t=0.5_r8*Td
        IF(unet_xdir.ge.umax) THEN
          T_c=Td
          T_t=0.0_r8
        ELSEIF(unet_xdir.le.umin) THEN
          T_c=0.0_r8
          T_t=Td
        ELSEIF(unet_xdir.lt.0.0_r8.and.unet_xdir.gt.umin) THEN
          delt=ASIN(-unet/umin)/pi
          T_t=T_t*(1.0_r8-2.0_r8*delt)
          T_c=Td-T_t
        ELSEIF(unet_xdir.gt.0.0_r8.and.unet_xdir.lt.umax) THEN
          delt=ASIN(unet_xdir/(-umax))/pi
          T_c=T_c*(1.0_r8-2.0_r8*delt)
          T_t=Td-T_c
        ELSEIF(unet_xdir.eq.0.0_r8) THEN
          T_c=T_c
          T_t=T_t
        ENDIF
      ELSEIF(RR.gt.0.5_r8) THEN
        T_c=T_c
        T_t=T_t
        IF(unet_xdir.ge.umax) THEN
          T_c=Td
          T_t=0.0_r8
        ELSEIF(unet_xdir<=umin) THEN
          T_c=0.0_r8
          T_t=Td
        ELSEIF(unet_xdir.lt.0.0_r8.and.unet_xdir.gt.umin) THEN
          delt=ASIN(-unet_xdir/umin)/pi
          T_t=T_t*(1.0_r8-2.0_r8*delt)
          T_c=Td-T_t
        ELSEIF(unet_xdir.gt.0.0_r8.and.unet_xdir.lt.umax) THEN
          delt=ASIN(unet_xdir/(-umax))/pi
          T_c=T_c*(1.0_r8-2.0_r8*delt)
          T_t=Td-T_c
        ELSEIF(unet_xdir.eq.0.0_r8) THEN
          T_c=T_c
          T_t=T_t
        ENDIF
      ENDIF
      T_c=MAX(T_c,0.0_r8)
      T_t=MAX(T_t,0.0_r8)
!
      END SUBROUTINE current_timeperiod
!
      SUBROUTINE surface_wave_mod(Td, depth, uhat,                      &
     &                            T_c, T_cu, T_t, T_tu)
! 
!**********************************************************************
!  This subroutine returns the following:
!  T_c, T_cu, T_t, T_tu  : Change in time period in crest and 
!                          trough cycle due to particle displacement
!                          under surface waves. 
!**********************************************************************
!
! Crest period extension for horizontal particle displacement.
! Tuning factor eps_eff = 0.55 from measurements GWK Schretlen 2010         
! Equations in Appendix B of SANTOSS Report 
! Report number: SANTOSS_UT_IR3 Date: January 2010 
!
      USE mod_kinds
      USE mod_scalars
      USE mod_vandera_funcs
!
      implicit none
! 
      real(r8), intent(in) :: Td, depth, uhat
      real(r8), intent(inout) :: T_c, T_cu, T_t, T_tu
!
!  Local variables
!
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8) :: k_wn, eps_eff, c
      real(r8) :: delta_Tc, delta_Tt
      real(r8) :: T_c_new, T_cu_new
      real(r8) :: T_t_new, T_tu_new
!
      k_wn=kh(Td,depth)/depth
      c=2.0_r8*pi/(k_wn*Td)
!
      eps_eff=0.55_r8 
!
      delta_Tc=eps_eff*uhat/(c*pi-eps_eff*2.0*uhat)
      T_c_new=T_c+delta_Tc
!
! Avoid non zero values for time periods 
!
      T_c_new=MAX( T_c_new, 0.0_r8)
      T_cu_new=MAX( T_cu*T_c_new/(T_c+eps), 0.0_r8 )
!
      delta_Tt=eps_eff*uhat/(c*pi+eps_eff*2.0*uhat)
      T_t_new=T_t-delta_Tt
      T_t_new=MAX( T_t_new, 0.0_r8)
      T_tu_new=MAX( T_tu*T_t_new/(T_t+eps), 0.0_r8 )
!
      T_c=T_c_new
      T_cu=T_cu_new
      T_t=T_t_new
      T_tu=T_tu_new
!
      END SUBROUTINE surface_wave_mod
!
      SUBROUTINE ripple_dim(psi, d50, eta, rlen)
! 
!**********************************************************************
!  This subroutine returns the following:
!  eta, rlen : Ripple dimensions: (height and length) 
!**********************************************************************
!
! Calculate ripple dimensions of O'Donoghue et al. 2006
! based on VA2013 Appendix B
!        
      USE mod_kinds
      USE mod_scalars 

      implicit none 
!
      real(r8), intent(in)  :: psi, d50
      real(r8), intent(out) :: eta, rlen
!
      real(r8) :: d50_mm 
      real(r8) :: m_eta, m_lambda, n_eta, n_lambda 
      real(r8), parameter :: eps=1.0E-14_r8
!     
      d50_mm=0.001_r8*d50
      IF(d50_mm.lt.0.22_r8) THEN
        m_eta=0.55_r8
        m_lambda=0.73_r8
      ELSEIF(d50_mm.ge.0.22_r8.and.d50_mm.lt.0.30_r8) THEN
        m_eta=0.55_r8+(0.45_r8*(d50_mm-0.22_r8)/(0.30_r8-0.22_r8))
        m_lambda=0.73_r8+(0.27_r8*(d50_mm-0.22_r8)/(0.30_r8-0.22_r8))
      ELSE
        m_eta=1.0_r8
        m_lambda=1.0_r8
      ENDIF
! 
! Smooth transition between ripple regime and bed sheet flow regime 
!
      IF(psi.le.190.0_r8) THEN
        n_eta=1.0_r8
      ELSEIF(psi.gt.190.0_r8.and.psi.lt.240.0_r8) THEN
        n_eta=0.5_r8*(1.0_r8+cos(pi*(psi-190.0_r8)/(50.0_r8)))
      ELSEIF(psi.ge.240.0_r8) THEN
        n_eta=0.0_r8
      ENDIF
      n_lambda=n_eta
!
      eta=MAX(0.0_r8,m_eta*n_eta*(0.275_r8-0.022*psi**0.42_r8))
!      rlen=MAX(0.0_r8,m_lambda*n_lambda*                                &
!     &                             (1.97_r8-0.44_r8*psi**0.21_r8))
      rlen=MAX(eps,m_lambda*n_lambda*                                   &
     &                             (1.97_r8-0.44_r8*psi**0.21_r8))
!
      RETURN
      END SUBROUTINE ripple_dim
!
      SUBROUTINE skewness_params( H_s, T, depth, r, phi, Ur )
!        
! Ruessink et al. provides equations for calculating skewness parameters
! Uses Malarkey and Davies equations to get "bb" and "r"
! Given input of H_s, T and depth 
! r     - skewness/asymmetry parameter r=2b/(1+b^2)            [value]
! phi   - skewness/asymmetry parameter                         [value]
! Su     - umax/(umax-umin)                                    [value]
! Au   - amax/(amax-amin)                                      [value]
! alpha - tmax/pi                                              [value]
!
      USE mod_kinds
      USE mod_scalars
      USE mod_vandera_funcs
!
      implicit none
!
      real(r8), intent(in)    :: H_s, T, depth
      real(r8), intent(inout) :: Ur
      real(r8), intent(out)   :: r, phi
!
! Local variables 
! 
      real(r8), parameter :: p1=0.0_r8
      real(r8), parameter :: p2=0.857_r8
      real(r8), parameter :: p3=-0.471_r8
      real(r8), parameter :: p4=0.297_r8
      real(r8), parameter :: p5=0.815_r8
      real(r8), parameter :: p6=0.672_r8
      real(r8), parameter :: eps=1.0E-14_r8
      real(r8) :: a_w
      real(r8) :: B, psi, bb
      real(r8) :: k_wn, cff
!      real(r8) :: kh_calc
      real(r8) :: Su, Au
!
! Ruessink et al., 2012, Coastal Engineering 65:56-63.
!
! k is the local wave number computed with linear wave theory.
!
      k_wn=kh(T,depth)/depth    
!
      a_w=0.5_r8*H_s 
      Ur=0.75_r8*a_w*k_wn/((k_wn*depth)**3.0_r8)+eps
!
! Ruessink et al., 2012 Equation 9.
!
      cff=EXP((p3-log10(Ur))/p4)
      B=p1+((p2-p1)/(1.0_r8+cff))
!
      psi=-90.0_r8*deg2rad*(1.0_r8-TANH(p5/Ur**p6))
!
# ifdef BEDLOAD_VANDERA_ASYM_LIMITS
      B=MIN(B,0.8554_r8) ! according to fig.2, max values w.r.t Ur=24
      psi=MAX(psi,-1.4233_r8) 
# endif     
! 
! Markaley and Davies, equation provides bb which is "b" in paper
! Check from where CRS found these equations
! 
      bb=sqrt(2.0_r8)*B/(sqrt(2.0_r8*B**2.0_r8+9.0_r8))
      r=2.0_r8*bb/(bb**2.0_r8+1.0_r8)
!
! Ruessink et al., 2012 under Equation 12.
!
      phi=-psi-0.5_r8*pi
!
! Where are these asymmetry Su, Au utilized 
! recreate the asymetry 
!          
      Su=B*cos(psi)
      Au=B*sin(psi)
!
      RETURN
      END SUBROUTINE skewness_params

      SUBROUTINE abreu_points( r, phi, Uw, T, T_c, T_t,                 &
     &                            T_cu, T_tu, umax, umin, RR, beta )
! 
!  Calculate umax, umin, and phases of asymmetrical wave orbital velocity 
!
!  Use the asymmetry parameters from Ruessink et al, 2012
!  to get the umax, umin and phases of asymettrical wave 
!  orbital velocity to be used by Van Der A. 
!  T_c is duration of crest
!  T_cu Duration of accerating flow within crest half cycle
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none
!
      real(r8), intent(in)  :: r, phi, Uw, T
      real(r8), intent(out) :: T_c, T_t, T_cu, T_tu
      real(r8), intent(out) :: umax, umin, RR, beta
!
! Local variables 
! 
      real(r8) :: b, c, ratio, tmt, tmc, tzd, tzu
      real(r8) :: omega, w, phi_new
      real(r8) :: P, F0, betar_0
!     real(r8) :: T_tu, T_cu, T_c, T_t
      real(r8) :: cff1, cff2, cff
      real(r8) :: Sk, Ak
!
      omega=2.0_r8*pi/T
!
      phi_new=-phi

! Malarkey and Davies (Under equation 16b) 
      P=SQRT(1.0_r8-r*r)
!
! Malarkey and Davies (Under equation 16b) 
!
      b=r/(1.0_r8+P)
!
! Appendix E of Malarkey and Davies 
!
      c=b*SIN(phi_new)
!
      cff1=4.0_r8*c*(b*b-c*c)+(1.0_r8-b*b)*(1.0_r8+b*b-2.0_r8*c*c)
      cff2=(1.0_r8+b*b)**2.0_r8-4.0_r8*c*c
      ratio=cff1/cff2
!
! These if conditionals prevent ASIN to be between [-1,1] and prevent NaNs
! Not a problem in the MATLAB code
!
      IF(ratio.gt.1.0_r8)THEN
        ratio=1.0_r8
      ENDIF
      IF(ratio.lt.-1.0_r8)THEN
        ratio=-1.0_r8
      ENDIF
      tmc=ASIN(ratio)
!
!
      cff1=4.0_r8*c*(b*b-c*c)-(1.0_r8-b*b)*(1.0_r8+b*b-2.0_r8*c*c)
      cff2=(1.0_r8+b*b)**2.0_r8-4.0_r8*c*c
      ratio=cff1/cff2
      IF(ratio.gt.1.0_r8)THEN
        ratio=1.0_r8
      ENDIF
      IF(ratio.lt.-1.0_r8)THEN
        ratio=-1.0_r8
      ENDIF
      tmt=ASIN(ratio)
!       
      IF(tmc.lt.0.0_r8) THEN
        tmc=tmc+2.0_r8*pi
      ENDIF
      IF(tmt.lt.0.0_r8) THEN
        tmt=tmt+2.0_r8*pi
      ENDIF
! 
! Non dimensional umax and umin, under E5 in Malarkey and Davies 
! 
      umax=1.0_r8+c
      umin=umax-2.0_r8
!
!     Dimensionalize
!
      umax=umax*Uw
      umin=umin*Uw
!
! phase of zero upcrossing and downcrossing (radians)
!
      tzu=ASIN(b*SIN(phi_new))
      tzd=2.0_r8*ACOS(c)+tzu
! 
! MD, equation 17
!
      RR=0.5_r8*(1.0_r8+b*SIN(phi_new))
! 
! MD, under equation 18
! 
      IF(r.le.0.5_r8) THEN
        F0=1.0_r8-0.27_r8*(2.0_r8*r)**(2.1_r8)
      ELSE
        F0=0.59_r8+0.14_r8*(2.0_r8*r)**(-6.2_r8)
      ENDIF
!
! MD, Equation 15a,b 
!
      IF(r.ge.0.0_r8.and.r.lt.0.5)THEN
        betar_0=0.5_r8*(1.0_r8+r)
      ELSEIF(r.gt.0.5_r8.and.r.lt.1.0_r8)THEN
        cff1=4.0_r8*r*(1.0_r8+r)
        cff2=cff1+1.0_r8
        betar_0=cff1/cff2
      ENDIF
!
! MD, Equation 18
!
      cff=SIN((0.5_r8*pi-ABS(phi_new))*F0)/SIN(0.5_r8*pi*F0)
      beta=0.5_r8+(betar_0-0.5_r8)*cff
!
! MD, Table 1, get asymmetry parameterization
! using GSSO (10a,b)
!
      cff=SQRT(2.0_r8*(1.0_r8+b*b)**3.0_r8)
      Sk=3.0_r8*SIN(phi_new)/cff
      Ak=-3.0_r8*COS(phi_new)/cff
!
! These are the dimensional fractions of wave periods needed by Van der A eqn.
!
      w=1.0_r8/omega
      T_c=(tzd-tzu)*w
      T_t=T-T_c
      T_cu=(tmc-tzu)*w
      T_tu=(tmt-tzd)*w
!
      RETURN
      END SUBROUTINE abreu_points
!
#endif
      END MODULE sed_bedload_vandera_mod



