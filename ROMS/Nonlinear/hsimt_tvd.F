#include "cppdefs.h"

      MODULE hsimt_tvd_mod
#if defined NONLINEAR && defined TS_HSIMT && defined SOLVE3D
!
!======================= Hui Wu, Tarandeep Kalra, and John C. Warner====
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!This routine computes anti-diffusive tracer flux based on HSIMT-TVD   !
!by Wu and Zhu (2010). This routine is for personal test only currently!
!                                                                      !
!  On Output: FX, FE                                                   !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Hui Wu and Jianrong Zhu (2010), Advection scheme with 3rd         !
!    high-order spatial interpolation at the middle temporal level     !
!    and its application to saltwater intrusion in the Changjiang      !
!    Estuary, Ocean Modelling 33, 33-51.                               !
!   Please contact Hui Wu (hwusklec@gmail.com) if have any questions   !
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
     &                           Huon_k, Hvom_k, oHz_k, t_k,            &
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
      real(r8), intent(in) :: Huon_k(LBi:,LBj:)
      real(r8), intent(in) :: Hvom_k(LBi:,LBj:)
      real(r8), intent(in) :: oHz_k(IminS:,JminS:)
      real(r8), intent(in) :: t_k(LBi:,LBj:)
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
      real(r8), intent(in) :: Huon_k(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hvom_k(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: oHz_k(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: t_k(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: FX(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(out) :: FE(IminS:ImaxS,JminS:JmaxS)
# endif
!
!  Local variable declarations.
!
      integer  :: i, is, j, k
      real(r8) :: cc1, cc2, cc3
      real(r8) :: sw_xi, rl, rkal, a1, b1, betal, rt, rkar, betar
      real(r8) :: sw_eta, rd, rkad, betad, ru, rkau, betau
      real(r8) :: cff, cff1, cff2, epson
      real(r8), dimension(IminS:ImaxS) :: kax, kax_inverse
      real(r8), dimension(IminS:ImaxS) :: grad_x
      real(r8), dimension(JminS:JmaxS) :: grad_y
      real(r8), dimension(JminS:JmaxS) :: kay, kay_inverse

# include "set_bounds.h"

!************Declare some constants locally***************************
      cc1=0.25_r8
      cc2=0.5_r8
      cc3=1.0_r8/12.0_r8
      epson=1.0E-12_r8
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+2
          grad_x(i)=(t_k(i,j)-t_k(i-1,j))
          cff=0.125_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))*         &
     &        (oHz_k(i-1,j)+oHz_k(i,j))
          kax(i)=(1.0_r8-abs(Huon_k(i,j)*dt(ng)*cff))
# ifdef MASKING
          grad_x(i)=grad_x(i)*umask(i,j)
          kax(i)=kax(i)*umask(i,j)
# endif
        END DO
        IF (.not.EWperiodic(ng)) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (Huon_k(Istr,j).ge.0.0_r8) THEN
              grad_x(Istr-1)=0.0_r8
              kax(Istr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (Huon_k(Iend+1,j).lt.0.0_r8) THEN
              grad_x(Iend+2)=0.0_r8
              kax(Iend+2)=0.0_r8
            END IF
          END IF
        END IF
        DO i=Istr,Iend+1
          IF (kax(i).le.epson) THEN
            kax_inverse(i)=0.0_r8
          ELSE
            kax_inverse(i)=1.0_r8/MAX(kax(i),epson)
          END IF
          IF (Huon_k(i,j).ge.0.0_r8) THEN
            IF (abs(grad_x(i)).le.epson) THEN
              rl=0.0_r8
              rkal=0.0_r8
            ELSE
              rl=grad_x(i-1)/(grad_x(i))
              rkal=kax(i-1)*kax_inverse(i)
!             cff1=Huon_k(i-1,j)*0.25_r8*(oHz_k(i-2,j)+oHz_k(i-1,j))*(pn(i-2,j)+pn(i-1,j))
!             cff2=Huon_k(i,j)*0.25_r8*(oHz_k(i-1,j)+oHz_k(i,j))*(pn(i-1,j)+pn(i,j))
!             rkal=kax(i-1)*kax_inverse(i)*cff1/cff2
            END IF
            a1= cc1*kax(i)+cc2-cc3*kax_inverse(i)
            b1=-cc1*kax(i)+cc2+cc3*kax_inverse(i)
            betal=a1+b1*rl
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rl*rkal,betal))*    &
     &                                  grad_x(i)*kax(i)
# ifdef MASKING
            cff=cff*rmask(i-2,j)
# endif
            sw_xi=t_k(i-1,j)+cff
          ELSE
            IF (abs(grad_x(i)).le.epson) THEN
              rt=0.0_r8
              rkar=0.0_r8
            ELSE
              rt=grad_x(i+1)/(grad_x(i))
              rkar=kax(i+1)*kax_inverse(i)
!             cff1=Huon_k(i+1,j)*0.25_r8*(oHz_k(i,j)+oHz_k(i+1,j))*(pn(i,j)+pn(i+1,j))
!             cff2=Huon_k(i,j)*0.25_r8*(oHz_k(i-1,j)+oHz_k(i,j))*(pn(i-1,j)+pn(i,j))
!             rkar=kax(i+1)*kax_inverse(i)*cff1/cff2
            END IF
            a1= cc1*kax(i)+cc2-cc3*kax_inverse(i)
            b1=-cc1*kax(i)+cc2+cc3*kax_inverse(i)
            betar=a1+b1*rt
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rt*rkar,betar))*    &
     &                                grad_x(i)*kax(i)
# ifdef MASKING
            cff=cff*rmask(i+1,j)
# endif
            sw_xi=t_k(i,j)-cff
          END IF
          FX(i,j)=sw_xi*huon_k(i,j)
        END DO
      END DO
!
      DO i=Istr,Iend
        DO j=Jstr-1,Jend+2
          grad_y(j)=(t_k(i,j)-t_k(i,j-1))
          cff=0.125_r8*(pn(i,j)+pn(i,j-1))*(pm(i,j)+pm(i,j-1))*         &
     &        (oHz_k(i,j)+oHz_k(i,j-1))
          kay(j)=(1.0_r8-abs(Hvom_k(i,j)*dt(ng)*cff))
# ifdef MASKING
          grad_y(j)=grad_y(j)*vmask(i,j)
          kay(j)=kay(j)*vmask(i,j)
# endif
        END DO
        IF (.not.NSperiodic(ng)) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (Hvom_k(i,Jstr).ge.0.0_r8) THEN
              grad_y(Jstr-1)=0.0_r8
              kay(Jstr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (Hvom_k(i,Jend+1).lt.0.0_r8) THEN
              grad_y(Jend+2)=0.0_r8
              kay(Jend+2)=0.0_r8
            END IF
          END IF
        END IF
        DO j=Jstr,Jend+1
          IF (kay(j).le.epson) THEN
            kay_inverse(j)=0.0_r8
          ELSE
            kay_inverse(j)=1.0_r8/MAX(kay(j),epson)
          END IF
          IF (Hvom_k(i,j).ge.0.0_r8) THEN
            IF (abs(grad_y(j)).le.epson) THEN
              rd=0.0_r8
              rkad=0.0_r8
            ELSE
              rd=grad_y(j-1)/grad_y(j)
              rkad=kay(j-1)*kay_inverse(j)
!             cff1=Hvom_k(i,j-1)*0.25_r8*(oHz_k(i,j-2)+oHz_k(i,j-1))*(pm(i,j-2)+pm(i,j-1))
!             cff2=Hvom_k(i,j)*0.25_r8*(oHz_k(i,j-1)+oHz_k(i,j))*(pm(i,j-1)+pm(i,j))
!             rkad=kay(j-1)*kay_inverse(j)*cff1/cff2
            END IF
            a1= cc1*kay(j)+cc2-cc3*kay_inverse(j)
            b1=-cc1*kay(j)+cc2+cc3*kay_inverse(j)
            betad=a1+b1*rd
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rd*rkad,betad))*    &
     &                              grad_y(j)*kay(j)
# ifdef MASKING
            cff=cff*rmask(i,j-2)
# endif
            sw_eta=t_k(i,j-1)+cff
          ELSE
            IF (abs(grad_y(j)).le.epson) THEN
              ru=0.0_r8
              rkau=0.0_r8
            ELSE
              ru=grad_y(j+1)/(grad_y(j))
              rkau=kay(j+1)*kay_inverse(j)
!             cff1=Hvom_k(i,j+1)*0.25_r8*(oHz_k(i,j+1)+oHz_k(i,j))*(pm(i,j+1)+pm(i,j))
!             cff2=Hvom_k(i,j)*0.25_r8*(oHz_k(i,j-1)+oHz_k(i,j))*(pm(i,j-1)+pm(i,j))
!             rkau=kay(j+1)*kay_inverse(j)*cff1/cff2
            END IF
            a1= cc1*kay(j)+cc2-cc3*kay_inverse(j)
            b1=-cc1*kay(j)+cc2+cc3*kay_inverse(j)
            betau=a1+b1*ru
            cff=0.5*max(0.0_r8,min(2.0_r8,2.0_r8*ru*rkau,betau))*       &
     &                            grad_y(j)*kay(j)
# ifdef MASKING
            cff=cff*rmask(i,j+1)
# endif
            sw_eta=t_k(i,j)-cff
          END IF
          FE(i,j)=sw_eta*hvom_k(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE hsimt_tvd_tile
#endif
      END MODULE hsimt_tvd_mod
