#include "cppdefs.h"
#define HSIMT_HOT

      MODULE hsimt_tvd_mod
#if defined NONLINEAR && defined TS_HSIMT && defined SOLVE3D
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2014 The ROMS/TOMS Group        John C. Warner   !
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
     &                           pm, pn , omn, z_r,                     &
     &                           Huon_k, Hvom_k, u_k, v_k, t_k,         &
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
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon_k(LBi:,LBj:)
      real(r8), intent(in) :: Hvom_k(LBi:,LBj:)
      real(r8), intent(in) :: u_k(IminS:,JminS:)
      real(r8), intent(in) :: v_k(IminS:,JminS:)
      real(r8), intent(in) :: t_k(LBi:,LBj:)
      real(r8), intent(out) :: FE(IminS:,JminS:)
      real(r8), intent(out) :: FX(IminS:,JminS:)
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
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon_k(LBi:,LBj:)
      real(r8), intent(in) :: Hvom_k(LBi:,LBj:)
      real(r8), intent(in) :: u_k(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: v_k(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: t_k(LBi:,LBj:)
      real(r8), intent(out) :: FE(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(out) :: FX(IminS:ImaxS,JminS:JmaxS)
# endif
!
!  Local variable declarations.
!
      integer  :: i, is, j, k
      real(r8) :: cc1, cc2, cc3
      real(r8) :: sw_xi, rl, rkal, a1, b1, betal, rt, rkar, betar
      real(r8) :: sw_eta, rd, rkad, betad, ru, rkau, betau
      real(r8) :: cff, epson
      real(r8), dimension(IminS:ImaxS) :: kax, kax_inverse
      real(r8), dimension(IminS:ImaxS) :: grad_x
      real(r8), dimension(JminS:JmaxS) :: grad_y
      real(r8), dimension(JminS:JmaxS) :: kay, kay_inverse

# include "set_bounds.h"

!************Declare some constants locally***************************
      cc1=0.25_r8
      cc2=0.5_r8
      cc3=1.0_r8/12.0_r8
      epson=1.e-10_r8
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+2
          grad_x(i)=t_k(i,j)-t_k(i-1,j)
          cff=0.5_r8*(pm(i-1,j)+pm(i,j))
          kax(i)=1.0_r8-abs(u_k(i,j)*dt(ng)*cff)+epson
        END DO 
        IF (.not.EWperiodic(ng)) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (u_k(Istr,j).ge.0.0_r8) THEN
              grad_x(Istr-1)=0.0_r8
              kax(Istr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (u_k(Iend+1,j).lt.0.0_r8) THEN
              grad_x(Iend+2)=0.0_r8 
              kax(Iend+2)=0.0_r8 
            END IF 
          END IF 
        END IF
        DO i=Istr,Iend+1
!         kax(i)=1.0_r8-abs(u_k(i,j)*dt(ng)*om_u_inverse)+epson
          kax_inverse(i)=1.0_r8/MAX(kax(i),epson)
          IF (u_k(i,j).ge.0.0_r8) THEN
            IF (abs(grad_x(i)).le.epson) THEN
              rl=0.0_r8
              rkal=0.0_r8
            ELSE                                           
              rl=grad_x(i-1)/(grad_x(i)+epson)
              rkal=kax(i-1)*kax_inverse(i)
            END IF
            a1=cc1*kax(i)+cc2-cc3*kax_inverse(i)
            b1=-cc1*kax(i)+cc2+cc3*kax_inverse(i)
            betal=a1+b1*rl
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rl*rkal,betal))*    &
     &                                  grad_x(i)*kax(i)
            
#  ifdef MASKING
            sw_xi=t_k(i-1,j)+rmask(i-2,j)*rmask(i,j)*cff
#  else
            sw_xi=t_k(i-1,j)+cff
#  endif
#  ifdef MASKING
            FX(i,j)=sw_xi*huon_k(i,j)*umask(i,j)
#  else
            FX(i,j)=sw_xi*huon_k(i,j)
#  endif
          ELSE
            IF (abs(grad_x(i)).le.epson) THEN
              rt=0.0_r8
              rkar=0.0_r8
            ELSE
              rt=grad_x(i+1)/(grad_x(i)+epson)
              rkar=kax(i+1)*kax_inverse(i)
            END IF
            a1=cc1*kax(i)+cc2-cc3*kax_inverse(i)
            b1=-cc1*kax(i)+cc2+cc3*kax_inverse(i)
            betar=a1+b1*rt
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rt*rkar,betar))*    &
     &                                grad_x(i)*kax(i)
#  ifdef MASKING
            sw_xi=t_k(i,j)-rmask(i+1,j)*rmask(i-1,j)*cff
#  else
            sw_xi=t_k(i,j)-cff
#  endif  
#  ifdef MASKING
            FX(i,j)=sw_xi*huon_k(i,j)*umask(i,j)
#  else
            FX(i,j)=sw_xi*huon_k(i,j)
#  endif
          END IF
        END DO
      END DO
!
      DO i=Istr,Iend
        DO j=Jstr-1,Jend+2
          grad_y(j)=t_k(i,j)-t_k(i,j-1)
          cff=0.5_r8*(pn(i,j-1)+pn(i,j))
          kay(j)=1.0_r8-abs(v_k(i,j)*dt(ng)*cff)+epson
        END DO
        IF (.not.NSperiodic(ng)) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (v_k(i,Jstr).ge.0.0_r8) THEN
              grad_y(Jstr-1)=0.0_r8
              kay(Jstr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (v_k(i,Jend+1).lt.0.0_r8) THEN
              grad_y(Jend+2)=0.0_r8
              kay(Jend+2)=0.0_r8
            END IF
          END IF
        END IF 
        DO j=Jstr,Jend+1
!         kay(j)=1.0_r8-abs(v_k(i,j)*dt(ng)*on_v_inverse)+epson
          kay_inverse(j)=1.0_r8/MAX(kay(j),epson)
          IF (v_k(i,j).ge.0.0_r8) THEN
            IF (abs(grad_y(j)).le.epson) THEN
              rd=0.0_r8
              rkad=0.0_r8
            ELSE
              rd=grad_y(j-1)/(grad_y(j)+epson)
              rkad=kay(j-1)*kay_inverse(j)
            END IF
            a1=cc1*kay(j)+cc2-cc3*kay_inverse(j)
            b1=-cc1*kay(j)+cc2+cc3*kay_inverse(j)
            betad=a1+b1*rd
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rd*rkad,betad))*    &
     &                              grad_y(j)*kay(j)
#  ifdef MASKING
            sw_eta=t_k(i,j-1)+rmask(i,j-2)*rmask(i,j)*cff
#  else
            sw_eta=t_k(i,j-1)+cff
#  endif
#  ifdef MASKING
            FE(i,j)=sw_eta*hvom_k(i,j)*vmask(i,j)
#  else
            FE(i,j)=sw_eta*hvom_k(i,j)
#  endif
          ELSE
            IF (abs(grad_y(j)).le.epson) THEN
              ru=0.0_r8
              rkau=0.0_r8
            ELSE
              ru=grad_y(j+1)/(grad_y(j)+epson)
              rkau=kay(j+1)*kay_inverse(j)
            END IF
            a1=cc1*kay(j)+cc2-cc3*kay_inverse(j)
            b1=-cc1*kay(j)+cc2+cc3*kay_inverse(j)
            betau=a1+b1*ru
            cff=0.5*max(0.0_r8,min(2.0_r8,2.0_r8*ru*rkau,betau))*       &
     &                            grad_y(j)*kay(j)
#  ifdef MASKING
            sw_eta=t_k(i,j)-rmask(i,j+1)*rmask(i,j-1)*cff
#  else
            sw_eta=t_k(i,j)-cff
#  endif
#  ifdef MASKING
            FE(i,j)=sw_eta*hvom_k(i,j)*vmask(i,j)
#  else
            FE(i,j)=sw_eta*hvom_k(i,j)
#  endif
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE hsimt_tvd_tile
#endif
      END MODULE hsimt_tvd_mod
