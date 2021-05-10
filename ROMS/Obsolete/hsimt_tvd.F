#include "cppdefs.h"

      MODULE hsimt_tvd_mod
#if defined NONLINEAR && defined SOLVE3D
!
!=======================================================================
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                     Hui Wu   !
!    See License_ROMS.txt                            Tarandeep Kalra   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine computes horizontal tracer advection fluxes FX and FE  !
!  using tthe third High-order Spatial Interpolation at the Middle     !
!  Temporal level (HSIMT; Wu and Zhu, 2010) with a Total Variation     !
!  Diminishing limiter scheme.                                         !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Hui Wu and Jianrong Zhu, 2010: Advection scheme with 3rd          !
!      high-order spatial interpolation at the middle temporal         !
!      level and its application to saltwater intrusion in the         !
!      Changjiang Estuary, Ocean Modelling 33, 33-51,                  !
!      doi:10.1016/j.ocemod.2009.12.001                                !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC :: hsimt_tvd_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE hsimt_tvd_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
# ifdef MASKING
   &                             rmask, umask, vmask,                   &
# endif
# ifdef WET_DRY
     &                           rmask_wet, umask_wet, vmask_wet,       &
# endif
     &                           pm, pn,                                &
     &                           Huon, Hvom, oHz, t,                    &
     &                           FX, FE)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
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
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Huon(LBi:,LBj:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:)
      real(r8), intent(in) :: oHz(IminS:,JminS:)
      real(r8), intent(in) :: t(LBi:,LBj:)
      real(r8), intent(out) :: FX(IminS:,JminS:)
      real(r8), intent(out) :: FE(IminS:,JminS:)
# else
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
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: oHz(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: FX(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(out) :: FE(IminS:ImaxS,JminS:JmaxS)
# endif
!
!  Local variable declarations.
!
      integer  :: i, is, j, k, ii, jj

      real(r8) :: cc1 = 0.25_r8
      real(r8) :: cc2 = 0.5_r8
      real(r8) :: cc3 = 1.0_r8/12.0_r8
      real(r8) :: eps1 = 1.0E-12_r8

      real(r8) :: cff, cff1
      real(r8) :: betaL, betaR, betaD, betaU
      real(r8) :: rL, rR, rD, rU, rkaL, rkaR, rkaD, rkaU
      real(r8) :: a1, b1, sw_eta, sw_xi

      real(r8), dimension(IminS:ImaxS) :: gradX, KaX, oKaX
      real(r8), dimension(JminS:JmaxS) :: gradE, KaE, oKaE

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute tracer horizontal aadvective fluxes using the HSIMT scheme
!  (Wu and Zhu, 2010).
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU-1,Iendp2
          cff=0.125_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*dt(ng)
          cff1=cff*(oHz(i-1,j)+oHz(i,j))
          gradX(i)=t(i,j)-t(i-1,j)
          KaX(i)=1.0_r8-ABS(Huon(i,j)*cff1)
# ifdef MASKING
          gradX(i)=gradX(i)*umask(i,j)
          KaX(i)=KaX(i)*umask(i,j)
# endif
        END DO
        IF (.not.EWperiodic(ng)) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (Huon(Istr,j).ge.0.0_r8) THEN
              gradX(Istr-1)=0.0_r8
              KaX(Istr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (Huon(Iend+1,j).lt.0.0_r8) THEN
              gradX(Iend+2)=0.0_r8
              KaX(Iend+2)=0.0_r8
            END IF
          END IF
        END IF
        DO i=Istr,Iend+1
          IF (KaX(i).le.eps1) THEN
            oKaX(i)=0.0_r8
          ELSE
            oKaX(i)=1.0_r8/MAX(KaX(i),eps1)
          END IF
          IF (Huon(i,j).ge.0.0_r8) THEN
            IF (ABS(gradX(i)).le.eps1) THEN
              rL=0.0_r8
              rkaL=0.0_r8
            ELSE
              rL=gradX(i-1)/gradX(i)
              rkaL=KaX(i-1)*oKaX(i)
            END IF
            a1= cc1*KaX(i)+cc2-cc3*oKaX(i)
            b1=-cc1*KaX(i)+cc2+cc3*oKaX(i)
            betaL=a1+b1*rL
            cff=0.5_r8*MAX(0.0_r8,                                      &
     &                     MIN(2.0_r8, 2.0_r8*rL*rkaL, betaL))*         &
     &          gradX(i)*KaX(i)
# ifdef MASKING
            ii=MAX(i-2,0)
            cff=cff*rmask(ii,j)
# endif
            sw_xi=t(i-1,j)+cff
          ELSE
            IF (ABS(gradX(i)).le.eps1) THEN
              rR=0.0_r8
              rkaR=0.0_r8
            ELSE
              rR=gradX(i+1)/gradX(i)
              rkaR=KaX(i+1)*oKaX(i)
            END IF
            a1= cc1*KaX(i)+cc2-cc3*oKaX(i)
            b1=-cc1*KaX(i)+cc2+cc3*oKaX(i)
            betaR=a1+b1*rR
            cff=0.5_r8*MAX(0.0_r8,                                      &
     &                     MIN(2.0_r8, 2.0_r8*rR*rkaR, betaR))*         &
     &          gradX(i)*KaX(i)
# ifdef MASKING
            ii=MIN(i+1,Lm(ng)+1)
            cff=cff*rmask(ii,j)
# endif
            sw_xi=t(i,j)-cff
          END IF
          FX(i,j)=sw_xi*Huon(i,j)
        END DO
      END DO
!
      DO i=Istr,Iend
        DO j=JstrV-1,Jendp2
          cff=0.125_r8*(pn(i,j)+pn(i,j-1))*(pm(i,j)+pm(i,j-1))*dt(ng)
          cff1=cff*(oHz(i,j)+oHz(i,j-1))
          gradE(j)=t(i,j)-t(i,j-1)
          KaE(j)=1.0_r8-ABS(Hvom(i,j)*cff1)
# ifdef MASKING
          gradE(j)=gradE(j)*vmask(i,j)
          KaE(j)=KaE(j)*vmask(i,j)
# endif
        END DO
        IF (.not.NSperiodic(ng)) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (Hvom(i,Jstr).ge.0.0_r8) THEN
              gradE(Jstr-1)=0.0_r8
              KaE(Jstr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (Hvom(i,Jend+1).lt.0.0_r8) THEN
              gradE(Jend+2)=0.0_r8
              KaE(Jend+2)=0.0_r8
            END IF
          END IF
        END IF
        DO j=Jstr,Jend+1
          IF (KaE(j).le.eps1) THEN
            oKaE(j)=0.0_r8
          ELSE
            oKaE(j)=1.0_r8/MAX(KaE(j),eps1)
          END IF
          IF (Hvom(i,j).ge.0.0_r8) THEN
            IF (ABS(gradE(j)).le.eps1) THEN
              rD=0.0_r8
              rkaD=0.0_r8
            ELSE
              rD=gradE(j-1)/gradE(j)
              rkaD=KaE(j-1)*oKaE(j)
            END IF
            a1= cc1*KaE(j)+cc2-cc3*oKaE(j)
            b1=-cc1*KaE(j)+cc2+cc3*oKaE(j)
            betaD=a1+b1*rD
            cff=0.5_r8*MAX(0.0_r8,                                      &
     &                     MIN(2.0_r8, 2.0_r8*rD*rkaD, betaD))*         &
     &          gradE(j)*KaE(j)
# ifdef MASKING
            jj=MAX(j-2,0)
            cff=cff*rmask(i,jj)
# endif
            sw_eta=t(i,j-1)+cff
          ELSE
            IF (ABS(gradE(j)).le.eps1) THEN
              rU=0.0_r8
              rkaU=0.0_r8
            ELSE
              rU=gradE(j+1)/gradE(j)
              rkaU=KaE(j+1)*oKaE(j)
            END IF
            a1= cc1*KaE(j)+cc2-cc3*oKaE(j)
            b1=-cc1*KaE(j)+cc2+cc3*oKaE(j)
            betaU=a1+b1*rU
            cff=0.5*MAX(0.0_r8,                                         &
     &                  MIN(2.0_r8, 2.0_r8*rU*rkaU, betaU))*            &
     &          gradE(j)*KaE(j)
# ifdef MASKING
            jj=MIN(j+1,Mm(ng)+1)
            cff=cff*rmask(i,jj)
# endif
            sw_eta=t(i,j)-cff
          END IF
          FE(i,j)=sw_eta*Hvom(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE hsimt_tvd_tile
#endif
      END MODULE hsimt_tvd_mod
