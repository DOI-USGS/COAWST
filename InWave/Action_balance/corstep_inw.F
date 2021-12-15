#include "cppdefs.h"
#undef  AC_HSIMT
#define AC_U3HADVECTION
#undef  AC_CENTERED2
#undef  AC_CENTERED4
#undef  AC_AKIMA4
      MODULE corstep_inw_mod
#if defined INWAVE_MODEL
!
!=======================================================================
!                                                                      !
!  This routine time-steps action density equations.                   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: corstep_inw

      CONTAINS
!
!***********************************************************************
      SUBROUTINE corstep_inw (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_inwave_vars
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
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 35)
# endif

      CALL corstep_inw_tile (ng, tile,                                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), nstp(ng), nnew(ng),               &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
# endif
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % pmon_u,                          &
     &                      GRID(ng) % pnom_v,                          &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
     &                      WAVEP(ng) % AC,                             &
     &                      WAVEP(ng) % cx,                             &
     &                      WAVEP(ng) % cy,                             &
     &                      WAVEP(ng) % ct,                             &
     &                      WAVEP(ng) % Tr,                             &
     &                      WAVEP(ng) % kwc,                            &
     &                      WAVEG(ng) % pd)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 35)
# endif
      RETURN
      END SUBROUTINE corstep_inw
!
!***********************************************************************
      SUBROUTINE corstep_inw_tile (ng, tile,                            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, nstp, nnew,                     &
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# ifdef WET_DRY
     &                            rmask_wet, umask_wet, vmask_wet,      &
# endif
     &                            pm, pn, on_u, om_u, om_v,             &
     &                            pmon_u, pnom_v,                       &
     &                            u, v,                                 &
     &                            AC, cx, cy, ct, Tr, kwc, pd)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_coupling
      USE mod_scalars
      USE mod_forces
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_inwave_swan
      USE exchange_3d_mod, ONLY : exchange_AC3d_tile

# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
      USE AC3dbc_mod, ONLY : AC3dbc_tile
# if defined WDISS_ROELVINK || defined WDISS_GAMMA
      USE dissip_inw_mod, ONLY : dissip_inw_tile
# endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
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
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
      real(r8), intent(in) :: cx(LBi:,LBj:,:)
      real(r8), intent(in) :: cy(LBi:,LBj:,:)
      real(r8), intent(in) :: ct(LBi:,LBj:,:)
      real(r8), intent(in) :: Tr(LBi:,LBj:,:)
      real(r8), intent(in) :: kwc(LBi:,LBj:,:)
      real(r8), intent(in) :: pd
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND,3)
      real(r8), intent(in) :: cx(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: cy(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: ct(LBi:UBi,LBj:UBj,0:ND)
      real(r8), intent(in) :: Tr(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: kwc(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: pd
# endif
!
!  Local variable declarations.
!
      integer :: i, is, itrc, j, k, d, ii, jj
      real(r8) :: eps = 1.0E-16_r8


      real(r8) :: cff, cff1, cff2, cff3, opd
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
      real(r8), dimension(IminS:ImaxS,0:ND+1) :: curvd
      real(r8), dimension(IminS:ImaxS,0:ND+2) :: FD
# ifdef AC_HSIMT
      real(r8) :: epson, sw
      real(r8) :: sw_xi, rl, rkal, a1, b1, betal, rt, rkar, betar
      real(r8) :: sw_eta, rd, rkad, betad, ru, rkau, betau
      real(r8), dimension(IminS:ImaxS) :: grad_x
      real(r8), dimension(JminS:JmaxS) :: grad_y
      real(r8), dimension(IminS:ImaxS) :: kax, kax_inverse
      real(r8), dimension(JminS:JmaxS) :: kay, kay_inverse
      real(r8), dimension(0:N(ng))     :: kaz, kaz_inverse, grad_k
      real(r8), dimension(0:ND+1)      :: grad_d, kad, kad_inverse
# endif

# include "set_bounds.h"

!  Compute horizontal AC advection fluxes.
!
      D_LOOP: DO d=1,ND
# if defined AC_CENTERED2
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=cx(i,j,d)*on_u(i,j)
            FX(i,j)=cff*0.5_r8*                                         &
     &              (AC(i-1,j,d,3)+                                     &
     &               AC(i  ,j,d,3))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=cy(i,j,d)*om_v(i,j)
            FE(i,j)=cff*0.5_r8*                                         &
     &              (AC(i,j-1,d,3)+                                     &
     &               AC(i,j  ,d,3))
          END DO
        END DO
# elif defined AC_U3HADVECTION || defined AC_CENTERED4 || defined AC_AKIMA4
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=AC(i  ,j,d,3)-                                    &
     &                AC(i-1,j,d,3)
#  ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Istr-1,j)=FX(Istr,j)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Iend+2,j)=FX(Iend+1,j)
              END DO
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
#  if defined AC_U3HADVECTION
              curv(i,j)=FX(i+1,j)-FX(i,j)
#  elif defined AC_AKIMA4
              cff=2.0_r8*FX(i+1,j)*FX(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  elif defined AC_CENTERED4
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=cx(i,j,d)*on_u(i,j)
#  if defined AC_U3HADVECTION
              FX(i,j)=cff*0.5_r8*                                       &
     &                (AC(i-1,j,d,3)+                                   &
     &                 AC(i  ,j,d,3))-                                  &
     &                cff1*(curv(i-1,j)*MAX(cff,0.0_r8)+                &
     &                      curv(i  ,j)*MIN(cff,0.0_r8))
#  elif defined AC_CENTERED4 || defined AC_AKIMA4
              FX(i,j)=cff*0.5_r8*                                       &
     &                (AC(i-1,j,d,3)+                                   &
     &                 AC(i  ,j,d,3)-                                   &
     &                 cff2*(grad(i  ,j)-                               &
     &                       grad(i-1,j)))
#  endif
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=AC(i,j  ,d,3)-                                    &
     &                AC(i,j-1,d,3)
#  ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jstr-1)=FE(i,Jstr)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jend+2)=FE(i,Jend+1)
              END DO
            END IF
          END IF
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
#  if defined AC_U3HADVECTION
              curv(i,j)=FE(i,j+1)-FE(i,j)
#  elif defined AC_AKIMA4
              cff=2.0_r8*FE(i,j+1)*FE(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  elif defined AC_CENTERED4
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=cy(i,j,d)*om_v(i,j)
#  if defined AC_U3HADVECTION
              FE(i,j)=cff*0.5_r8*                                       &
     &                (AC(i,j-1,d,3)+                                   &
     &                 AC(i,j  ,d,3))-                                  &
     &                cff1*(curv(i,j-1)*MAX(cff,0.0_r8)+                &
     &                      curv(i,j  )*MIN(cff,0.0_r8))
#  elif defined AC_CENTERED4 || defined AC_AKIMA4
              FE(i,j)=cff*0.5_r8*                                       &
     &                (AC(i,j-1,d,3)+                                   &
     &                 AC(i,j  ,d,3)-                                   &
     &                 cff2*(grad(i,j  )-                               &
     &                       grad(i,j-1)))
#  endif
            END DO
          END DO
# elif defined AC_HSIMT
      epson=1.0E-12_r8
!
      DO j=Jstr,Jend
        DO i=Istrm1,Iendp2
          grad_x(i)=(AC(i,j,d,3)-AC(i-1,j,d,3))
          cff=0.5_r8*(pm(i-1,j)+pm(i,j))
          kax(i)=1.0_r8-abs(cx(i,j,d))*dt(ng)*cff*WAVES(ng)%dtifac
#  ifdef MASKING
          grad_x(i)=grad_x(i)*umask(i,j)
          kax(i)=kax(i)*umask(i,j)
#  endif
        END DO
        IF (.not.EWperiodic(ng)) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            IF (cx(Istr,j,d).ge.0.0_r8) THEN
              grad_x(Istr-1)=0.0_r8
              kax(Istr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            IF (cx(Iend+1,j,d).lt.0.0_r8) THEN
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
          IF (cx(i,j,d).ge.0.0_r8) THEN
            IF (abs(grad_x(i)).le.epson) THEN
              rl=0.0_r8
              rkal=0.0_r8
            ELSE
              rl=grad_x(i-1)/(grad_x(i))
              rkal=kax(i-1)*kax_inverse(i)
            END IF
            a1= cc1*kax(i)+cc2-cc3*kax_inverse(i)
            b1=-cc1*kax(i)+cc2+cc3*kax_inverse(i)
            betal=a1+b1*rl
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rl*rkal,betal))*    &
     &                                  grad_x(i)*kax(i)
#  ifdef MASKING
            ii=MAX(i-2,0)
            cff=cff*rmask(ii,j)
#  endif
            sw_xi=AC(i-1,j,d,3)+cff
          ELSE
            IF (abs(grad_x(i)).le.epson) THEN
              rt=0.0_r8
              rkar=0.0_r8
            ELSE
              rt=grad_x(i+1)/(grad_x(i))
              rkar=kax(i+1)*kax_inverse(i)
            END IF
            a1= cc1*kax(i)+cc2-cc3*kax_inverse(i)
            b1=-cc1*kax(i)+cc2+cc3*kax_inverse(i)
            betar=a1+b1*rt
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rt*rkar,betar))*    &
     &                                grad_x(i)*kax(i)
#  ifdef MASKING
            ii=MIN(i+1,Lm(ng)+1)
            cff=cff*rmask(ii,j)
#  endif
            sw_xi=AC(i,j,d,3)-cff
          END IF
          FX(i,j)=sw_xi*cx(i,j,d)*on_u(i,j)
        END DO
      END DO
!
      DO i=Istr,Iend
        DO j=Jstrm1,Jendp2
          grad_y(j)=(AC(i,j,d,3)-AC(i,j-1,d,3))
          cff=0.5_r8*(pn(i,j)+pn(i,j-1))
          kay(j)=1.0_r8-abs(cy(i,j,d))*dt(ng)*cff*WAVES(ng)%dtifac
#  ifdef MASKING
          grad_y(j)=grad_y(j)*vmask(i,j)
          kay(j)=kay(j)*vmask(i,j)
#  endif
        END DO
        IF (.not.NSperiodic(ng)) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            IF (cy(i,Jstr,d).ge.0.0_r8) THEN
              grad_y(Jstr-1)=0.0_r8
              kay(Jstr-1)=0.0_r8
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            IF (cy(i,Jend+1,d).lt.0.0_r8) THEN
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
          IF (cy(i,j,d).ge.0.0_r8) THEN
            IF (abs(grad_y(j)).le.epson) THEN
              rd=0.0_r8
              rkad=0.0_r8
            ELSE
              rd=grad_y(j-1)/grad_y(j)
              rkad=kay(j-1)*kay_inverse(j)
            END IF
            a1= cc1*kay(j)+cc2-cc3*kay_inverse(j)
            b1=-cc1*kay(j)+cc2+cc3*kay_inverse(j)
            betad=a1+b1*rd
            cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rd*rkad,betad))*    &
     &                              grad_y(j)*kay(j)
#  ifdef MASKING
            jj=MAX(j-2,0)
            cff=cff*rmask(i,jj)
#  endif
            sw_eta=AC(i,j-1,d,3)+cff
          ELSE
            IF (abs(grad_y(j)).le.epson) THEN
              ru=0.0_r8
              rkau=0.0_r8
            ELSE
              ru=grad_y(j+1)/(grad_y(j))
              rkau=kay(j+1)*kay_inverse(j)
            END IF
            a1= cc1*kay(j)+cc2-cc3*kay_inverse(j)
            b1=-cc1*kay(j)+cc2+cc3*kay_inverse(j)
            betau=a1+b1*ru
            cff=0.5*max(0.0_r8,min(2.0_r8,2.0_r8*ru*rkau,betau))*       &
     &                            grad_y(j)*kay(j)
#  ifdef MASKING
            jj=MIN(j+1,Mm(ng)+1)
            cff=cff*rmask(i,jj)
#  endif
            sw_eta=AC(i,j,d,3)-cff
          END IF
          FE(i,j)=sw_eta*cy(i,j,d)*om_v(i,j)
        END DO
      END DO
# endif
!
!  Time-step horizontal advection term.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)*WAVES(ng)%dtifac
            cff1=cff*(FX(i+1,j)-FX(i,j)+                                &
     &                FE(i,j+1)-FE(i,j))
            AC(i,j,d,nnew)=AC(i,j,d,nstp)-cff1
          END DO
        END DO
# if defined AC_CENTERED2 || defined AC_CENTERED4 || defined AC_AKIMA4
!
! here lets try horiz mix
!
         cff2=2.0_r8 ! for a4test chan   
!        cff2=0.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=0.5_r8*cff2*pmon_u(i,j)
              FX(i,j)=cff*(AC(i,j,d,nnew)-AC(i-1,j,d,nnew))
#  ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#  endif
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=0.5_r8*cff2*pnom_v(i,j)
              FE(i,j)=cff*(AC(i,j,d,nnew)-AC(i,j-1,d,nnew))
#  ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j  )-FX(i,j))
              cff2=cff*(FE(i  ,j+1)-FE(i,j))
              cff3=cff1+cff2
              AC(i,j,d,nnew)=AC(i,j,d,nnew)+cff3
            END DO
          END DO
# endif
      END DO D_LOOP
!
!-----------------------------------------------------------------------
!  Time-step directional advection term.
!-----------------------------------------------------------------------
!
      opd=1.0_r8/pd
      J_LOOP: DO j=Jstr,Jend
# if defined AC_U3HADVECTION || defined AC_CENTERED2 || defined AC_CENTERED4 || \
     defined AC_AKIMA4
        DO i=Istr,Iend
#  if defined THETA_AC_PERIODIC
            FD(i,0)=AC(i,j,ND  ,3)-                                     &
     &              AC(i,j,ND-1,3)
            FD(i,1)=AC(i,j,1   ,3)-                                     &
     &              AC(i,j,ND  ,3)
#  else
          FD(i,0)=0.0_r8
          FD(i,1)=0.0_r8
#  endif
          DO d=2,ND
            FD(i,d)=AC(i,j,d  ,3)-                                      &
     &              AC(i,j,d-1,3)
          END DO
#  if defined THETA_AC_PERIODIC
          FD(i,ND+1)=FD(i,1)
          FD(i,ND+2)=FD(i,2)
#  else
          FD(i,ND+1)=0.0_r8
          FD(i,ND+2)=0.0_r8
#  endif
        END DO
!
        DO i=Istr,Iend
          DO d=0,ND+1
            curvd(i,d)=FD(i,d+1)-FD(i,d)
          END DO
        END DO
!
        cff1=1.0_r8/6.0_r8
        cff2=1.0_r8/3.0_r8
        DO i=Istr,Iend
          DO d=1,1
#  if defined THETA_AC_PERIODIC
            cff=ct(i,j,d)*opd
#  else
#   if defined THETA_AC_WALL
            cff=0.0_r8
#   else
            cff=ct(i,j,d)*opd
#   endif
#  endif
            FD(i,d)=cff*0.5_r8*                                         &
#  if defined THETA_AC_PERIODIC
     &              (AC(i,j,ND,3)+                                      &
     &               AC(i,j,d  ,3))-                                    &
#  else
     &              (AC(i,j,d     ,3)+                                  &
     &               AC(i,j,d  ,3))-                                    &
#  endif
     &               cff1*(curvd(i,d-1)*MAX(cff,0.0_r8)+                &
     &                     curvd(i,d  )*MIN(cff,0.0_r8))
          END DO
          DO d=2,ND
            cff=ct(i,j,d)*opd
            FD(i,d)=cff*0.5_r8*                                         &
     &              (AC(i,j,d-1,3)+                                     &
     &               AC(i,j,d  ,3))-                                    &
     &               cff1*(curvd(i,d-1)*MAX(cff,0.0_r8)+                &
     &                     curvd(i,d  )*MIN(cff,0.0_r8))
          END DO
          DO d=ND+1,ND+1
#  if defined THETA_AC_PERIODIC
            cff=ct(i,j,d)*opd
#  else
#   if defined THETA_AC_WALL
            cff=0.0_r8
#   else
            cff=ct(i,j,d)*opd
#   endif
#  endif
            FD(i,d)=cff*0.5_r8*                                         &
#  if defined THETA_AC_PERIODIC
     &              (AC(i,j,ND,3)+                                      &
     &               AC(i,j,1  ,3))-                                    &
#  else
     &              (AC(i,j,ND,3)+                                      &
     &               AC(i,j,ND,3))-                                     &
#  endif
     &               cff1*(curvd(i,d-1)*MAX(cff,0.0_r8)+                &
     &                     curvd(i,d  )*MIN(cff,0.0_r8))
          END DO
        END DO
# elif defined AC_HSIMT
        epson=1.0E-12_r8
!
        DO i=Istr,Iend
          DO d=2,ND
            grad_d(d)=(AC(i,j,d,3)-AC(i,j,d-1,3))
            kad(d)=1.0_r8-abs(ct(i,j,d))*dt(ng)*WAVES(ng)%dtifac
          END DO

          DO d=0,0
#  if defined THETA_AC_PERIODIC
            grad_d(d)=(AC(i,j,ND,3)-AC(i,j,ND-1,3))
#  else
            grad_d(d)=0.0_r8
#  endif
            kad(d)=0.0_r8
          END DO
          DO d=1,1
#  if defined THETA_AC_PERIODIC
            grad_d(d)=(AC(i,j,d,3)-AC(i,j,ND,3))
#  else
            grad_d(d)=0.0_r8
#  endif
            kad(d)=1.0_r8-abs(ct(i,j,d))*dt(ng)*WAVES(ng)%dtifac
          END DO
          DO d=ND+1,ND+1
#  if defined THETA_AC_PERIODIC
            grad_d(d)=(AC(i,j,1,3)-AC(i,j,d-1,3))
#  else
            grad_d(d)=0.0_r8
#  endif
            kad(d)=0.0_r8
          END DO
          DO d=1,1
            IF (kad(d).le.epson) THEN
              kad_inverse(d)=0.0_r8
            ELSE
              kad_inverse(d)=1.0_r8/MAX(kad(d),epson)
            END IF
#  if defined THETA_AC_WALL
            FD(i,d)=0.0_r8
#  else
            IF (ct(i,j,d).ge.0.0_r8) THEN
              IF (abs(grad_d(d)).le.epson) THEN
                rl=0.0_r8
                rkal=0.0_r8
              ELSE
                rl=grad_d(d-1)/(grad_d(d))
                rkal=kad(d-1)*kad_inverse(d)
              END IF
              a1= cc1*kad(d)+cc2-cc3*kad_inverse(d)
              b1=-cc1*kad(d)+cc2+cc3*kad_inverse(d)
              betal=a1+b1*rl
              cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rl*rkal,betal))*    &
     &                                    grad_d(d)*kad(d)
#   if defined THETA_AC_PERIODIC
              sw_xi=AC(i,j,ND,3)+cff
#   else
              sw_xi=0.0_r8
#   endif
            ELSE
              IF (abs(grad_d(d)).le.epson) THEN
                rt=0.0_r8
                rkar=0.0_r8
              ELSE
                rt=grad_d(d+1)/(grad_d(d))
                rkar=kad(d+1)*kad_inverse(d)
              END IF
              a1= cc1*kad(d)+cc2-cc3*kad_inverse(d)
              b1=-cc1*kad(d)+cc2+cc3*kad_inverse(d)
              betar=a1+b1*rt
              cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rt*rkar,betar))*    &
     &                                  grad_d(d)*kad(d)
              sw_xi=AC(i,j,d,3)-cff
            END IF
            FD(i,d)=sw_xi*ct(i,j,d)*opd
#  endif
          END DO

          DO d=2,ND
            IF (kad(d).le.epson) THEN
              kad_inverse(d)=0.0_r8
            ELSE
              kad_inverse(d)=1.0_r8/MAX(kad(d),epson)
            END IF
            IF (ct(i,j,d).ge.0.0_r8) THEN
              IF (abs(grad_d(d)).le.epson) THEN
                rl=0.0_r8
                rkal=0.0_r8
              ELSE
                rl=grad_d(d-1)/(grad_d(d))
                rkal=kad(d-1)*kad_inverse(d)
              END IF
              a1= cc1*kad(d)+cc2-cc3*kad_inverse(d)
              b1=-cc1*kad(d)+cc2+cc3*kad_inverse(d)
              betal=a1+b1*rl
              cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rl*rkal,betal))*    &
     &                                    grad_d(d)*kad(d)
              sw_xi=AC(i,j,d-1,3)+cff
            ELSE
              IF (abs(grad_d(d)).le.epson) THEN
                rt=0.0_r8
                rkar=0.0_r8
              ELSE
                rt=grad_d(d+1)/(grad_d(d))
                rkar=kad(d+1)*kad_inverse(d)
              END IF
              a1= cc1*kad(d)+cc2-cc3*kad_inverse(d)
              b1=-cc1*kad(d)+cc2+cc3*kad_inverse(d)
              betar=a1+b1*rt
              cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rt*rkar,betar))*    &
     &                                  grad_d(d)*kad(d)
              sw_xi=AC(i,j,d,3)-cff
            END IF
            FD(i,d)=sw_xi*ct(i,j,d)*opd
          END DO

          DO d=ND+1,ND+1
            IF (kad(d).le.epson) THEN
              kad_inverse(d)=0.0_r8
            ELSE
              kad_inverse(d)=1.0_r8/MAX(kad(d),epson)
            END IF
#  if defined THETA_AC_WALL
            FD(i,d)=0.0_r8
#  else
            IF (ct(i,j,d).ge.0.0_r8) THEN
              IF (abs(grad_d(d)).le.epson) THEN
                rl=0.0_r8
                rkal=0.0_r8
              ELSE
                rl=grad_d(d-1)/(grad_d(d))
                rkal=kad(d-1)*kad_inverse(d)
              END IF
              a1= cc1*kad(d)+cc2-cc3*kad_inverse(d)
              b1=-cc1*kad(d)+cc2+cc3*kad_inverse(d)
              betal=a1+b1*rl
              cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rl*rkal,betal))*    &
     &                                    grad_d(d)*kad(d)
              sw_xi=AC(i,j,d-1,3)+cff
            ELSE
              IF (abs(grad_d(d)).le.epson) THEN
                rt=0.0_r8
                rkar=0.0_r8
              ELSE
                rt=grad_d(d+1)/(grad_d(d))
                rkar=kad(d+1)*kad_inverse(d)
              END IF
              a1= cc1*kad(d)+cc2-cc3*kad_inverse(d)
              b1=-cc1*kad(d)+cc2+cc3*kad_inverse(d)
              betar=a1+b1*rt
              cff=0.5_r8*max(0.0_r8,min(2.0_r8,2.0_r8*rt*rkar,betar))*    &
     &                                  grad_d(d)*kad(d)
#   if defined THETA_AC_PERIODIC
              sw_xi=AC(i,j,1,3)-cff
#   else
              sw_xi=0.0_r8
#   endif
            END IF
            FD(i,d)=sw_xi*ct(i,j,d)*opd
#  endif
          END DO
        END DO
# endif
!
!  Time-step directional advection (m Tunits).
!
        DO d=1,ND
          DO i=Istr,Iend
            cff=dt(ng)*pd*WAVES(ng)%dtifac
            AC(i,j,d,nnew)=AC(i,j,d,nnew)-                              &
     &                     cff*(FD(i,d+1)-FD(i,d))
            AC(i,j,d,nnew)=MAX(0.0_r8,AC(i,j,d,nnew))
          END DO
        END DO
      END DO J_LOOP
# if defined WDISS_ROELVINK || defined WDISS_GAMMA
      CALL dissip_inw_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, nnew,                                 &
#  ifdef MASKING
     &                      rmask,                                      &
#  endif
#  ifdef WET_DRY
     &                      rmask_wet, umask_wet, vmask_wet,            &
#  endif
     &                      FORCES(ng) % Dissip_break,                  &
     &                      FORCES(ng) % Dissip_wcap,                   &
     &                      FORCES(ng) % Lwave,                         &
#  ifdef WAVES_BOT_PERIOD
     &                      FORCES(ng) % Pwave_bot,                     &
#  endif
     &                      WAVEP(ng) % h_tot,                          &
     &                      AC, Tr, kwc)
# endif
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions.
!-----------------------------------------------------------------------
      CALL AC3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  nstp, nnew,                                     &
     &                  AC)
!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_AC3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, ND,             &
     &                           AC(:,:,:,nnew))
      END IF
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND,                    &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    AC(:,:,:,nnew))
# endif
      RETURN
      END SUBROUTINE corstep_inw_tile
#endif
      END MODULE corstep_inw_mod

