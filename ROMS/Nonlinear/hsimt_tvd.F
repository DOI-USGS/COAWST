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
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Hui Wu and Jianrong Zhu (2010), Advection scheme with 3rd         !
!		high-order spatial interpolation at the middle temporal level  !
!		and its application to saltwater intrusion in the Changjiang   !
!		Estuary, Ocean Modelling 33, 33-51.                            !
!   Please contact Hui Wu (hwusklec@gmail.com) if have any questions   !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC :: hsimt_tvd_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE hsimt_tvd_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
# ifdef MASKING
     &                              rmask, umask, vmask,                &
# endif
# ifdef WET_DRY
     &                              rmask_wet, umask_wet, vmask_wet,    &
# endif
!     &                              pm, pn, omn, om_u, on_v,u_k,v_k,   &
     &                              pm, pn , omn, u_k, v_k,             &
     &                              z_r,                                &
     &                              Huon_k, Hvom_k, t_k,                &
     &                              FX,FE)
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
!      real(r8), intent(in) :: om_u(LBi:,LBj:)
!      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon_k(LBi:,LBj:)
      real(r8), intent(in) :: Hvom_k(LBi:,LBj:)
      real(r8), intent(in) :: t_k(LBi:,LBj:)
      real(r8), intent(in) :: u_k(LBi:,LBj:)
      real(r8), intent(in) :: v_k(LBi:,LBj:)
      real(r8), intent(out) :: FE(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(out) :: FX(IminS:ImaxS,JminS:JmaxS)
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
!      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
!      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon_k(LBi:,LBj:)
      real(r8), intent(in) :: Hvom_k(LBi:,LBj:)
      real(r8), intent(in) :: t_k(LBi:,LBj:)
      real(r8), intent(in) :: u_k(LBi:,LBj:)
      real(r8), intent(in) :: v_k(LBi:,LBj:)
      real(r8), intent(out) :: FE(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(out) :: FX(IminS:ImaxS,JminS:JmaxS)
# endif
!
!  Local variable declarations.
!
      integer  :: i, is, j, k
      real(r8) :: cc1,cc2,cc3
      real(r8) :: om_u_inverse,on_v_inverse
      real(r8) :: sw_xi,rl,rkal,a1,b1,betal,rtr,rkar,betar
      real(r8) :: sw_eta,rd,rkad,betad,ru,rkau,betau,epson
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kax,kax_inverse
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kay,kay_inverse

# include "set_bounds.h"

!************Declare some constants locally***************************
      cc1=0.25_r8
      cc2=0.5_r8
      cc3=1.0_r8/12.0_r8
      epson=1.e-10_r8

      DO j=JstrV-1,Jendp1
        DO i=IstrU-1,Iendp2
          om_u_inverse = 0.5_r8*(pm(i-1,j) + pm(i,j))
          kax(i,j)=1.0_r8-abs(u_k(i,j)*dt(ng)*om_u_inverse)+epson
          kax_inverse(i,j)=1.0_r8/kax(i,j)
        END DO
      END DO

      DO j=JstrV-1,Jendp2
        DO i=IstrU-1,Iendp1
          on_v_inverse = 0.5_r8*(pn(i,j-1) + pn(i,j))
          kay(i,j)=1.0_r8-abs(v_k(i,j)*dt(ng)*on_v_inverse)+epson
          kay_inverse(i,j)=1.0_r8/kay(i,j)
        END DO
      END DO

      DO j=JstrV-1,Jendp1
        DO i=IstrU-1,Iendp2
            IF (u_k(i,j).ge.0.0_r8) THEN
              IF (i.le.3) THEN
                sw_xi=t_k(i-1,j)
              ELSE
#  ifdef MASKING
                IF (rmask(i-2,j).eq.0.0_r8) THEN
                  sw_xi=t_k(i-1,j)
                ELSE
                  IF (rmask(i,j).eq.0.0_r8) THEN
                    sw_xi=t_k(i-1,j)
                  ELSE
#   endif
                    IF (abs(t_k(i,j)-t_k(i-1,j)).le.epson) THEN
                      rl=0.0_r8
                      rkal=0.0_r8
                    ELSE
                      rl=(t_k(i-1,j)-t_k(i-2,j))/                           &
     &                   (t_k(i,j)-t_k(i-1,j))
                      rkal=kax(i-1,j)*kax_inverse(i,j)
                    ENDIF
                    a1=cc1*kax(i,j)+cc2-cc3*kax_inverse(i,j)
                    b1=-cc1*kax(i,j)+cc2+cc3*kax_inverse(i,j)
                    betal=a1+b1*rl
                    sw_xi=t_k(i-1,j)+                                       &
     &                 0.5_r8*max(0.0_r8,                                   &
     &                                 min(2.0_r8,2.0_r8*rl*rkal,betal))*   &
     &                                 (t_k(i,j)-t_k(i-1,j))*kax(i,j)
#  ifdef MASKING
                  END IF
                END IF
#  endif
              END IF
              FX(i,j)=sw_xi*huon_k(i,j)
#  ifdef MASKING
              FX(i,j)=sw_xi*huon_k(i,j)*umask(i,j)
#  endif

            ELSE
              IF (i.ge.Lm(ng)-2) THEN
                sw_xi=t_k(i,j)
              ELSE
#  ifdef MASKING
                IF (rmask(i+1,j).eq.0.0_r8) THEN
                  sw_xi=t_k(i,j)
                ELSE
                  IF (rmask(i-1,j).eq.0.0_r8) THEN
                    sw_xi=t_k(i,j)
                  ELSE
#   endif
                    IF (abs(t_k(i,j)-t_k(i-1,j)).le.epson) THEN
                      rtr=0.0_r8
                      rkar=0.0_r8
                    ELSE
                      rtr=(t_k(i+1,j)-t_k(i,j))/                            &
     &                    (t_k(i,j)-t_k(i-1,j))
                      rkar=kax(i+1,j)*kax_inverse(i,j)
                    ENDIF
                    a1=cc1*kax(i,j)+cc2-cc3*kax_inverse(i,j)
                    b1=-cc1*kax(i,j)+cc2+cc3*kax_inverse(i,j)
                    betar=a1+b1*rtr
                    sw_xi=t_k(i,j)-                                         &
     &                 0.5_r8*max(0.0_r8,                                   &
     &                            min(2.0_r8,2.0_r8*rtr*rkar,betar))*       &
     &                            (t_k(i,j)-t_k(i-1,j))*kax(i,j)
#  ifdef MASKING
                  END IF  
                END IF
#   endif
              END IF
              FX(i,j)=sw_xi*huon_k(i,j)
#  ifdef MASKING
              FX(i,j)=sw_xi*huon_k(i,j)*umask(i,j)
#  endif
            END IF
        END DO
      END DO

      DO j=JstrV-1,Jendp2
        DO i=IstrU-1,Iendp1
          IF (v_k(i,j).ge.0.0_r8) THEN
            IF (j.le.3) THEN
              sw_eta=t_k(i,j-1)
            ELSE
#  ifdef MASKING
              IF (rmask(i,j-2).eq.0.0_r8) THEN
                sw_eta=t_k(i,j-1)
              ELSE
                IF (rmask(i,j).eq.0.0_r8) THEN
                  sw_eta=t_k(i,j-1)
                ELSE
#   endif   
                  IF (abs(t_k(i,j)-t_k(i,j-1)).le.epson) THEN
                    rd=0.0_r8
                    rkad=0.0_r8
                  ELSE
                    rd=(t_k(i,j-1)-t_k(i,j-2))/                       &
     &                 (t_k(i,j)-t_k(i,j-1))
                    rkad=kay(i,j-1)*kay_inverse(i,j)
                  END IF
                  a1=cc1*kay(i,j)+cc2-cc3*kay_inverse(i,j)
                  b1=-cc1*kay(i,j)+cc2+cc3*kay_inverse(i,j)
                  betad=a1+b1*rd
                  sw_eta=t_k(i,j-1)+                                  &
     &                   0.5_r8*max(0.0_r8,                           &
     &                          min(2.0_r8,2.0_r8*rd*rkad,betad))*    &
     &                           (t_k(i,j)-t_k(i,j-1))*kay(i,j)
#  ifdef MASKING
                END IF
              END IF
#  endif
            END IF
              FE(i,j)=sw_eta*hvom_k(i,j)
#  ifdef MASKING
              FE(i,j)=sw_eta*hvom_k(i,j)*vmask(i,j)
#  endif          
          ELSE
            IF (j.ge.Mm(ng)-2) THEN
              sw_eta=t_k(i,j)
            ELSE
#  ifdef MASKING
              IF (rmask(i,j+1).eq.0.0_r8) THEN
                  sw_eta=t_k(i,j)
              ELSE
                IF (rmask(i,j-1).eq.0.0_r8) THEN
                  sw_eta=t_k(i,j)
                ELSE
#  endif
                  IF (abs(t_k(i,j)-t_k(i,j-1)).le.epson) THEN
                    ru=0.0_r8
                    rkau=0.0_r8
                  ELSE
                    ru=(t_k(i,j+1)-t_k(i,j))/                        &
     &                 (t_k(i,j)-t_k(i,j-1))
                    rkau=kay(i,j+1)*kay_inverse(i,j)
                  END IF
                  a1=cc1*kay(i,j)+cc2-cc3*kay_inverse(i,j)
                  b1=-cc1*kay(i,j)+cc2+cc3*kay_inverse(i,j)
                  betau=a1+b1*ru
                  sw_eta=t_k(i,j)-                                   &
     &                   0.5*max(0.0_r8,                             &
     &                   min(2.0_r8,2.0_r8*ru*rkau,betau))* &
     &                   (t_k(i,j)-t_k(i,j-1))*kay(i,j)
#  ifdef MASKING  
                END IF
              END IF
#  endif
            END IF
              FE(i,j)=sw_eta*hvom_k(i,j)
#  ifdef MASKING
              FE(i,j)=sw_eta*hvom_k(i,j)*vmask(i,j)
#  endif
          END IF 
        END DO
      END DO

      RETURN
      END SUBROUTINE hsimt_tvd_tile
#endif
      END MODULE hsimt_tvd_mod
