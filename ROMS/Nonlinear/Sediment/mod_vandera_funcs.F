#include "cppdefs.h"
      MODULE MOD_VANDERA_FUNCS
#if defined SEDIMENT && defined BEDLOAD_VANDERA
!
!git $Id$
!================================================ Tarandeep S. Kalra ===
!  Copyright (c) 2002-2026 The ROMS Group             Chris Sherwood   !
!    Licensed under a MIT/X style license             John C. Warner   !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module contains several functions that are required for the    !
!  sediment bedload calculations using Van der A's formulations.       !     
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Soulsby, R.L., 1997: Dynamics of Marine Sands: A Manual for         !
!    Practical Applications, Thomas Telford, London.                   !
!                                                                      !
!  Soulsby, R.L., 2006: Simplified calculation of wave orbital         !
!    velocities, Technical Report, KR Hallingford.                     !
!                                                                      !
!  VA2013 Appendix A.                                                  !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_scalars
!
      implicit none
!      
      PUBLIC  :: kh              ! wave number
      PUBLIC  :: w_s_calc        ! critical shields parameter
      PUBLIC  :: w_sc_calc       ! water particle vertical velocity
      PUBLIC  :: mu_calc         ! bed roughness factor from grain size
      PUBLIC  :: ksd_calc        ! current-related bed roughness
      PUBLIC  :: ksw_calc        ! wave-related bed roughness
      PUBLIC  :: fw_calc         ! full-cycle wave friction factor
      PUBLIC  :: fd_calc_santoss ! current-related frction factor
      PUBLIC  :: fd_calc_madsen  ! current-related frction factor
      PUBLIC  :: fwi_calc        ! wave friction factor
      PUBLIC  :: dsf_calc        ! sheet flow thickness
      PUBLIC  :: theta_cr_calc   ! critical shilds parameter, Soulsby
      PRIVATE
!
      CONTAINS 
!
      FUNCTION kh (Td, depth) RESULT (wn)
!
!=======================================================================
!                                                                      !
!  Calculate Wave Number from Wave period and depth (Soulsby, 2006).   !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: Td, depth
!
!  Local variable declarations.
!
      integer  :: iter
!
      real(r8) :: cff
      real(r8) :: x, y, t, omega
      real(r8) :: wn
!
!-----------------------------------------------------------------------
!  Compute wave number.
!-----------------------------------------------------------------------
!
      omega=2.0_r8*pi/Td
!
!  Depth (i.e. x) cannot go negative to avoid negative wave number.
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
!  Iteratively solving 3 times for equation 7 of Soulsby (1997) by
!  using equations 12a-14.
!      
      DO iter=1,3
        t=TANH(y)
        cff=(y*t-x)/(t+y*(1.0_r8-t*t))
        y=y-cff
      END DO
      wn=y
!
      RETURN
      END FUNCTION kh
! 
      FUNCTION w_s_calc (d50, rhos) RESULT (w_s)
!
!=======================================================================
!                                                                      !
!  Computes Critical Shields parameter (Soulsby, 1997).                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: d50, rhos
!
!  Local variable declarations.
!
      real(r8), parameter :: nu=1.36E-6_r8
!
      real(r8) :: cff, cff1, cff2
      real(r8) :: s, dstar
      real(r8) :: w_s
!
!-----------------------------------------------------------------------
!  Compute critical shields parameter.
!-----------------------------------------------------------------------
!
      s=rhos/rho0
      dstar=(g*(s-1)/(nu*nu))**(1.0_r8/3.0_r8)*d50
      cff=nu/d50
      cff1=10.36_r8
      cff2=dstar*dstar*dstar
      w_s=cff*(SQRT(cff1*cff1+1.049_r8*cff2)-cff1)
!
      RETURN
      END FUNCTION w_s_calc
!
      FUNCTION w_sc_calc (Hs, Td, depth, RR, w_s, zws) RESULT (worb)
! 
!=======================================================================
!                                                                      !
!  Uses second order Stokes theory to get vertical velocity of water   !
!  particle at a given elevation based on "santoss_core.m" script.     !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
!
      real(r8), intent(in) :: Hs, Td, depth, RR, zws, w_s
!
!  Local variable declarations.
!
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: eps_inv = 1.0E14_r8
!
      real(r8) :: cff, worb1, worb2, worb 
!
!-----------------------------------------------------------------------
!  Compute water particle vertical velocity at given elevation.
!-----------------------------------------------------------------------
!
      worb1=pi*Hs*zws/(Td*depth+eps)
      worb2=worb1*2.0_r8*(RR+RR-1.0_r8)
!
!  Using the SANTOSS model formulation.
! 
      cff=1.0_r8/8.0_r8
      worb=cff*worb1*SQRT(64.0_r8-(-worb1+                              &
     &                             SQRT(worb1**2+32.0_r8*               &
     &                                  worb2**2+eps))**2/              &
     &                            (worb2**2+eps))+                      &
     &         worb2*SIN(2.0_r8*ACOS(cff*(-worb1+                       &
     &                                    SQRT(worb1**2+                &
     &                                         32.0_r8*worb2**2))/      &
     &                                   (worb2+eps)))
!
! Prevent worb from going to Infinity when worb2=0.0 
! 
      worb=MIN(worb, eps_inv)
!
      RETURN  
      END FUNCTION w_sc_calc
!
      FUNCTION mu_calc (d50) RESULT (mu)
! 
!=======================================================================
!                                                                      !
!  It calculates bed roughness factor based on grain size using        !
!  VA2013 Appendix A formulation.                                      !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
!
      real(r8), intent(in) :: d50
!
!  Local variable declarations.
!
      real(r8) :: d50_mm
      real(r8) :: mu
!
!-----------------------------------------------------------------------
!  Compute bed roughness factor.
!-----------------------------------------------------------------------
!
      d50_mm=d50*1000.0_r8
!
      IF (d50_mm.le.0.15_r8) THEN
        mu=6.0_r8
      ELSEIF(d50_mm.gt.0.15_r8.and.d50_mm.lt.0.20_r8) THEN
        mu=6.0_r8-5.0_r8*((d50_mm-0.15_r8)/(0.2_r8-0.15_r8))
      ELSEIF(d50_mm.ge.0.20_r8) THEN
        mu=1.0_r8
      END IF
!
      RETURN
      END FUNCTION mu_calc
!
      FUNCTION ksd_calc (d50, d90, mu, theta_timeavg, eta, rlen)        &
     &           RESULT (ksd)
! 
!=======================================================================
!                                                                      !
!  It calculates current-related bed roughness from VA2013             !
!  Appendix A.1.                                                       !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
!
      real(r8), intent(in)    :: d50, d90, mu, theta_timeavg, eta
      real(r8), intent(inout) :: rlen       ! HGA weird, why overwrite?
!                                             seems like a bug
!  Local variable declarations.
!
      real(r8) :: ripple_fac
      real(r8) :: ksd
!
!-----------------------------------------------------------------------
!  Compute current-related bed roughness.
!-----------------------------------------------------------------------
!
      rlen=MAX(rlen,d50)                   ! notice overwritting
      ripple_fac=0.4_r8*eta**2.0_r8/rlen
      ksd=MAX(3.0_r8*d90,                                              &
     &        d50*(mu+6.0_r8*(theta_timeavg-1.0_r8)))+                 &
     &    ripple_fac
!
      RETURN
      END FUNCTION ksd_calc
!
      FUNCTION ksw_calc (d50, mu, theta_timeavg, eta, rlen)             &
     &           RESULT (ksw)
! 
!=======================================================================
!                                                                      !
!  It calculates wave related bed roughness from VA2013 Appendix A,    !
!  Equation A.5                                                        !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
!
      real(r8), intent(in)    :: d50, mu, theta_timeavg, eta
      real(r8), intent(inout) :: rlen       ! HGA weird, why overwrite?
!                                             seems like a bug
!  Local variable declarations.
!
      real(r8) :: ripple_fac
      real(r8) :: ksw
!
!-----------------------------------------------------------------------
!  Compute wave-related bed roughness.
!-----------------------------------------------------------------------
!
      rlen=MAX(rlen,d50)                   ! notice overwritting
      ripple_fac=0.4_r8*eta**2.0_r8/rlen
      ksw=MAX(d50,                                                      &
     &        d50*(mu+6.0_r8*(theta_timeavg-1.0_r8)))+                  &
     &    ripple_fac
!
      RETURN
      END FUNCTION ksw_calc
!
      FUNCTION fw_calc (ahat, ksw) RESULT (fw)
! 
!=======================================================================
!                                                                      !
!  It calculates full-cycle wave friction factor from from VA2013      !
!  Appendix A, Equation A.4.                                           !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
!
      real(r8), intent(in) :: ahat, ksw
!
!  Local variable declarations.
!
      real(r8) :: ratio
      real(r8) :: fw
!
!-----------------------------------------------------------------------
!  Compute full-cycle wave friction factor.
!-----------------------------------------------------------------------
!
      ratio=ahat/ksw
      IF (ratio.gt.1.587_r8) THEN
        fw=0.00251_r8*EXP(5.21_r8*(ratio)**(-0.19_r8))
      ELSE
        fw=0.3_r8
      ENDIF
!
      RETURN
      END FUNCTION fw_calc
!
      FUNCTION fd_calc_santoss (udelta, delta, ksd) RESULT (fd)
! 
!=======================================================================
!                                                                      !
!  It calculates current-related friction factor by Assuming           !
!  logarithmic velocity profile (VA2013, Equation 20).                 !
!                                                                      !
!  It is not used since it is currently computed from  bottom          !
!  current stress.                                                     !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
! 
      real(r8) :: udelta, delta, ksd
!
!  Local variable declarations.
!
      real(r8), parameter :: min_udelta=1.0E-4_r8
      real(r8), parameter :: von_k=0.41_r8
!
      real(r8) :: fd
!
!-----------------------------------------------------------------------
!  Compute current-related friction factor.
!-----------------------------------------------------------------------
!
      IF (udelta.lt.min_udelta) THEN
        fd=0.0_r8
      ELSE
        fd=2.0_r8*(von_k/LOG(30.0_r8*delta/ksd))**2.0_r8
      ENDIF
!     
      RETURN
      END FUNCTION fd_calc_santoss
!
      FUNCTION fd_calc_madsen (udelta, mag_bstrc) RESULT (fd)
! 
!=======================================================================
!                                                                      !
!  It calculates current-related friction factor directly from the     !
!  current stresses.                                                   !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
! 
      real(r8), intent(in) :: udelta, mag_bstrc
!
!  Local variable declarations.
!
      real(r8), parameter :: eps=1.0E-14_r8
      real(r8), parameter :: min_udelta=1.0E-4_r8
!
      real(r8) :: fd
!
!-----------------------------------------------------------------------
!  Compute current related friction factor 
!-----------------------------------------------------------------------
! 
      IF (udelta.lt.min_udelta) THEN
        fd=0.0_r8
      ELSE
!!      fd=MAX((mag_bstrc/(0.5_r8*udelta*udelta*rho0)),eps)
!!      fd=MIN(fd, 2.0_r8) 
        fd=MAX((mag_bstrc/(0.5_r8*udelta*udelta)),eps)
      END IF
!     
      RETURN
      END FUNCTION fd_calc_madsen
!
      FUNCTION fwi_calc (T_iu, T_i, ahat, ksw) RESULT (fwi)
! 
!=======================================================================
!                                                                      !
!  It calculates wave friction factor for wave and crest half cycle,   !
!  VA2013 Equation 21.                                                 !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
!
      real(r8), intent(in) :: T_iu, T_i, ahat, ksw
!
!  Local variable declarations.
!
      real(r8), parameter :: eps = 1.0E-14_r8
!
      real(r8) :: c1, ratio
      real(r8) :: cff1, cff2, cff3
      real(r8) :: fwi
!
!-----------------------------------------------------------------------
!  Compute wave friction factor.
!-----------------------------------------------------------------------
!
      fwi=0.3_r8
!      
      c1=2.6_r8
      ratio=ahat/ksw
      IF (ratio.gt.1.587_r8) THEN
        cff1=MAX( (T_iu/(T_i+eps)),0.0_r8 )
        cff2=(2.0_r8*cff1)**c1
        cff3=cff2*ratio
!
! Prevents arithematic overflow error when 0.0**(-0.19_r8) 
!
        IF (cff3.le.0.0_r8) THEN 
          fwi=0.0_r8
        ELSE 
          fwi=0.00251_r8*EXP(5.21_r8*(cff3)**(-0.19_r8))
        END IF
      END IF 
!
      RETURN
      END FUNCTION fwi_calc
!
      FUNCTION dsf_calc (d50, theta_i) RESULT (dsf)
! 
!=======================================================================
!                                                                      !
!  It calculates sheet flow thickness, VA2013 Appendix C.1.            !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
!
      real(r8), intent(in) :: d50, theta_i
!
!  Local variable declarations.
!
      real(r8) :: d50_mm
      real(r8) :: cff
      real(r8) :: dsf
!
!-----------------------------------------------------------------------
!  Compute sheet flow thickness.
!-----------------------------------------------------------------------
!
      d50_mm=d50*1000.0_r8
      IF (d50_mm.le.0.15_r8) THEN
        cff=25.0_r8*theta_i
      ELSE IF((d50_mm.gt.0.15_r8).and.(d50_mm.lt.0.20_r8)) THEN
        cff=25.0_r8-(12.0_r8*(d50_mm-0.15_r8)/0.05_r8)
      ELSE IF (d50_mm.ge.0.20_r8) THEN
        cff=13.0_r8*theta_i
      ENDIF
      dsf=MAX(d50*cff,d50)
!
      RETURN
      END FUNCTION dsf_calc
!
      FUNCTION theta_cr_calc (d50, rhos) RESULT (theta)
!
!=======================================================================
!                                                                      !
!  It calculates critical chields parameter from Soulsby (1997).       !
!                                                                      !
!=======================================================================
!   
!  Imported variable declarations.
!
      real(r8), intent(in) :: d50, rhos
!
!  Local variable declarations.
!
      real(r8), parameter :: nu=1.36E-6_r8
!
      real(r8) :: s, dstar
      real(r8) :: cff1, cff2
      real(r8) :: theta
!
!-----------------------------------------------------------------------
!  Compute critical shields parameter.
!-----------------------------------------------------------------------
!
      s=rhos/rho0
      dstar=(g*(s-1)/(nu*nu))**(1.0_r8/3.0_r8)*d50
      cff1=0.30_r8/(1.0_r8+1.2_r8*dstar)
      cff2=0.055_r8*(1.0_r8-EXP(-0.020_r8*dstar))
      theta=cff1+cff2
!
      RETURN
      END FUNCTION theta_cr_calc
!
#endif
      END MODULE mod_vandera_funcs
