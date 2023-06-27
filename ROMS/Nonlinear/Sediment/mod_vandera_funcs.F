#include "cppdefs.h"
      Module MOD_VANDERA_FUNCS
#if defined SEDIMENT && defined BEDLOAD_VANDERA
!
!svn $Id: mod_bedload_vandera_funcs.F 429 2018-04-10 17:30:26Z arango $
!======================================================================!
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!----------------------------------------------Tarandeep S. Kalra------!
!------------------------------------------------Chris Sherwood -------!
!----------------------------------------------- John C. Warner--------!
!----------------------------------------------------------------------!
!  This module contains several functions that are required for the    !
!  sediment bedload calculations using Van der A's formulations.       !     
!                                                                      !
!----------------------------------------------------------------------!
!======================================================================!
!
      implicit none
      
      CONTAINS 
!
      REAL(r8) FUNCTION kh(Td,depth)
!
!  Calculate wave number from Wave period and depth 
!
! RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
! HR Wallingford Report TR 155, February 2006
!
      USE mod_scalars
!
      implicit none
! 
      real(r8) :: Td, depth
      real(r8) :: cff
      real(r8) :: x, y, t, omega
!
      omega=2.0_r8*pi/Td
        
!
! depth (i.e.x) cannot go negative to avoid negative wave number
!
      IF(depth.lt.0.0_r8) THEN 
        x=0.0_r8
      ELSE 
        x=omega**2.0_r8*depth/g
      ENDIF 
!
      IF(x.lt.1.0_r8) THEN
        y=SQRT(x)
      ELSE
        y=x
      ENDIF
!
! Iteratively solving 3 times for eqn.7 of Soulsby 1997 by using 
! eqns. (12a-14)
!      
      t=TANH(y)
      cff=(y*t-x)/(t+y*(1.0_r8-t*t))
      y=y-cff
!
      t=TANH(y)
      cff=(y*t-x)/(t+y*(1.0_r8-t*t))
      y=y-cff

      t=TANH(y)
      cff=(y*t-x)/(t+y*(1.0_r8-t*t))
      y=y-cff
      kh=y
!
      RETURN
      END FUNCTION
! 
      REAL(r8) FUNCTION w_s_calc(d50, rhos)
!
! Critical Shields parameter from Soulsby (1997).
! Dynamics of Marine Sands 
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none
! 
      real(r8), parameter :: nu=1.36E-6_r8
      real(r8) :: d50, rhos
      real(r8) :: s, dstar
      real(r8) :: cff, cff1
!
      s=rhos/rho0
      dstar=(g*(s-1)/(nu*nu))**(1.0_r8/3.0_r8)*d50
      cff=nu/d50
      cff1=10.36_r8
      w_s_calc=cff*(sqrt(cff1*cff1+1.049_r8*dstar**3.0_r8)-cff1)
!
      RETURN
      END FUNCTION
!
      REAL(r8) FUNCTION w_sc_calc(Hs, Td, depth, RR, w_s, zws) 
! 
! Second order Stokes theory to get vertical velocity of water particle
! at a given elevation based on santoss_core.m
! 
      USE mod_kinds
      USE mod_scalars
!
      implicit none 
!     
      real(r8), parameter :: eps_inv=1.0E14_r8
      real(r8) :: Hs, Td, depth, RR, zws, w_s
      real(r8) :: cff, worb1, worb2, worb 
      real(r8), parameter :: eps = 1.0E-14_r8
!
      worb1=pi*Hs*zws/(Td*depth+eps)
      worb2=worb1*2.0_r8*(RR+RR-1.0_r8)
!
!  Using the SANTOSS model formulation 
! 
      cff=1.0_r8/8.0_r8
      worb=cff*worb1*SQRT(64.0_r8-(-worb1+                              &
     &          SQRT(worb1**2+32.0_r8*                                  &
     &               worb2**2+eps))**2/(worb2**2+eps))+                 &
     &               worb2*SIN(2.0_r8*ACOS(cff*(-worb1+                 &
     &          SQRT(worb1**2+32.0_r8*worb2**2))/                       &
     &              (worb2+eps)))

!
! Prevent worb from going to Infinity when worb2=0.0 
! 
      worb=MIN(worb, eps_inv)
      w_sc_calc=worb
!
      RETURN  
      END FUNCTION w_sc_calc
!
      REAL(r8) FUNCTION mu_calc(d50)
!
! Calculate bed roughness factor based on grain size
! VA2013 Appendix A., required for current related bed roughness
! and wave related bed roughness. 
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none
!
      real(r8) :: d50, d50_mm
!
      d50_mm=d50*1000.0_r8
!
      IF(d50_mm.le.0.15_r8) THEN
        mu_calc=6.0_r8
      ELSEIF(d50_mm.gt.0.15_r8.and.d50_mm.lt.0.20_r8) THEN
        mu_calc=6.0_r8-5.0_r8*((d50_mm-0.15_r8)/(0.2_r8-0.15_r8))
      ELSEIF(d50_mm.ge.0.20_r8) THEN
        mu_calc=1.0_r8
      ENDIF
!
      RETURN
      END FUNCTION mu_calc
!
      REAL(r8) FUNCTION ksd_calc(d50, d90, mu, theta_timeavg,           &
     &                         eta, rlen)
!
! Calculate current-related bed roughness from VA2013 Appendix A.1.
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none
!
      real(r8) :: d50, d90, mu, theta_timeavg, eta, rlen
      real(r8) :: ripple_fac
!
      rlen=MAX(rlen,d50)
      ripple_fac=0.4_r8*eta**2.0_r8/rlen
      ksd_calc=MAX( 3.0_r8*d90,                                         &
     &             d50*(mu+6.0_r8*(theta_timeavg-1.0_r8)) )+            &
     &             ripple_fac
!
      RETURN
      END FUNCTION ksd_calc
!
      REAL(r8) FUNCTION ksw_calc(d50, mu, theta_timeavg, eta, rlen)
!
! Calculate wave related bed roughness from VA2013 Eqn. A.5.
! 
      USE mod_kinds
      USE mod_scalars
!
      implicit none
      real(r8) :: d50, mu, theta_timeavg, eta, rlen
      real(r8) :: ripple_fac, ksw
!
      rlen=MAX(rlen,d50)
!
      ripple_fac=0.4_r8*eta**2.0_r8/rlen
      ksw_calc=MAX( d50,                                                &
     &                d50*(mu+6.0_r8*(theta_timeavg-1.0_r8)) )          &
     &               +ripple_fac
!
      RETURN
      END FUNCTION ksw_calc
!
      REAL(r8) FUNCTION fw_calc(ahat, ksw)
!
! Calculate full-cycle wave friction factor from VA2013 Eqn. A.4. 
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none
      real(r8) :: ahat, ksw, ratio, fw
!
      ratio=ahat/ksw
      IF(ratio.gt.1.587_r8) THEN
        fw_calc=0.00251_r8*EXP(5.21_r8*(ratio)**(-0.19_r8))
      ELSE
        fw_calc=0.3_r8
      ENDIF
!
      RETURN
      END FUNCTION fw_calc
!
! This function is not getting used at the time because
! we compute current frction factor from bottom current stress 
! 
      REAL(r8) FUNCTION fd_calc_santoss(udelta, delta, ksd)

      USE mod_kinds
      USE mod_scalars

      implicit none

! Calculate current related friction factor VA2013 Eqn. 20
! Assuming logarithmic velocity profile.
 
      real(r8) :: udelta, delta, ksd
      real(r8), parameter :: min_udelta=1.0E-4_r8
      real(r8), parameter :: von_k=0.41_r8
 
      IF(udelta.lt.min_udelta) THEN
        fd_calc_santoss=0.0_r8
      ELSE
        fd_calc_santoss=2.0_r8*(von_k/LOG(30.0_r8*delta/ksd))**2.0_r8
      ENDIF
     
      RETURN
      END FUNCTION fd_calc_santoss
!
      REAL(r8) FUNCTION fd_calc_madsen(udelta, mag_bstrc)
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none
!
! Calculate current related friction factor 
! directly from the current stresses. 
! 
      real(r8), parameter :: eps=1.0E-14_r8
      real(r8), parameter :: min_udelta=1.0E-4_r8
      real(r8) :: udelta, mag_bstrc
! 
      IF(udelta.lt.min_udelta) THEN
        fd_calc_madsen=0.0_r8
      ELSE
        !fd_calc_new=MAX((mag_bstrc/(0.5_r8*udelta*udelta*rho0)),eps)
        fd_calc_madsen=MAX((mag_bstrc/(0.5_r8*udelta*udelta)),eps)
!        fd_calc_new=MIN(fd_calc_new,2.0_r8) 
        fd_calc_madsen=fd_calc_madsen
      ENDIF
!     
      RETURN
      END FUNCTION fd_calc_madsen
!
      REAL(r8) FUNCTION fwi_calc(T_iu, T_i, ahat, ksw)
!
! Wave friction factor for wave and crest half cycle VA2013 Eqn. 21.
!       
      USE mod_kinds
      USE mod_scalars
!
      implicit none
!
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8) :: T_iu, T_i, ahat, ksw
      real(r8) :: c1, ratio, fwi
      real(r8) :: cff1, cff2, cff3
!
      fwi_calc=0.3_r8
!      
      c1=2.6_r8
      ratio=ahat/ksw
      IF(ratio.gt.1.587_r8) THEN
        cff1=MAX( (T_iu/(T_i+eps)),0.0_r8 )
        cff2=(2.0_r8*cff1)**c1
        cff3=cff2*ratio
!
! These if condition prevents arithematic overflow error
! when 0.0**-0.19_r8, think about that 
!
        IF(cff3.le.0.0_r8) THEN 
          fwi_calc=0.0_r8
        ELSE 
          fwi_calc=0.00251_r8*EXP(5.21_r8*(cff3)**(-0.19_r8))
        END IF
      END IF 
!
      RETURN
      END FUNCTION fwi_calc
!
      REAL(r8) FUNCTION dsf_calc(d50, theta_i)
!
! Sheet flow thickness VA2013 Appendix C.1.
!       
      USE mod_kinds
      USE mod_scalars
!
      implicit none
!
      real(r8) :: d50, theta_i
      real(r8) :: d50_mm
      real(r8) :: cff
!
      d50_mm=d50*1000.0_r8
      IF(d50_mm.le.0.15_r8)THEN
        cff=25.0_r8*theta_i
      ELSEIF(d50_mm.gt.0.15_r8.and.d50_mm.lt.0.20_r8)THEN
        cff=25.0_r8-(12.0_r8*(d50_mm-0.15_r8)/0.05_r8)
      ELSEIF(d50_mm.ge.0.20_r8)THEN
        cff=13.0_r8*theta_i
      ENDIF
      dsf_calc=MAX(d50*cff,d50)
!
      RETURN
      END FUNCTION dsf_calc
!
      REAL(r8) FUNCTION theta_cr_calc(d50, rhos)
!
! Critical Shields parameter from Soulsby (1997).
!
      USE mod_kinds
      USE mod_scalars

      implicit none
! 
      real(r8), parameter :: nu=1.36E-6_r8
      real(r8) :: d50, rhos
      real(r8) :: s, dstar
      real(r8) :: cff1, cff2
!
      s=rhos/rho0
      dstar=(g*(s-1)/(nu*nu))**(1.0_r8/3.0_r8)*d50
      cff1=0.30_r8/(1.0_r8+1.2_r8*dstar)
      cff2=0.055_r8*(1.0_r8-EXP(-0.020_r8*dstar))
      theta_cr_calc=cff1+cff2
!
      RETURN
      END FUNCTION theta_cr_calc
!
#endif
      END MODULE mod_vandera_funcs
