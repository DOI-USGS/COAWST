      SUBROUTINE bblm (ng, tile)
!
!svn $Id: sg_bbl.h 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group        Richard Styles   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine compute bottom stresses for the case when the wave     !
!  solution in the wave boundary layer is based on a  2-layer eddy     !
!  viscosity that is linear increasing above Zo and constant above     !
!  Z1.                                                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Styles, R. and S.M. glenn,  2000: Modeling stratified wave and      !
!    current bottom boundary layers in the continental shelf, JGR,     !
!    105, 24119-24139.                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_bbl
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
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
      CALL wclock_on (ng, iNLM, 37)
# endif
      CALL bblm_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                nrhs(ng),                                         &
     &                GRID(ng) % h,                                     &
     &                GRID(ng) % z_r,                                   &
     &                GRID(ng) % z_w,                                   &
     &                GRID(ng) % angler,                                &
     &                GRID(ng) % ZoBot,                                 &
# if defined SG_CALC_UB
     &                FORCES(ng) % Hwave,                               &
# else
     &                FORCES(ng) % Uwave_rms,                           &
# endif
     &                FORCES(ng) % Dwave,                               &
     &                FORCES(ng) % Pwave_bot,                           &
     &                OCEAN(ng) % rho,                                  &
     &                OCEAN(ng) % u,                                    &
     &                OCEAN(ng) % v,                                    &
     &                SEDBED(ng) % bottom,                              &
     &                BBL(ng) % Iconv,                                  &
     &                BBL(ng) % Ubot,                                   &
     &                BBL(ng) % Vbot,                                   &
     &                BBL(ng) % Ur,                                     &
     &                BBL(ng) % Vr,                                     &
     &                BBL(ng) % bustrc,                                 &
     &                BBL(ng) % bvstrc,                                 &
     &                BBL(ng) % bustrw,                                 &
     &                BBL(ng) % bvstrw,                                 &
     &                BBL(ng) % bustrcwmax,                             &
     &                BBL(ng) % bvstrcwmax,                             &
     &                FORCES(ng) % bustr,                               &
     &                FORCES(ng) % bvstr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 37)
# endif

      RETURN
      END SUBROUTINE bblm
!
!***********************************************************************
      SUBROUTINE bblm_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs,                                       &
     &                      h, z_r, z_w, angler, ZoBot,                 &
# if defined SG_CALC_UB
     &                      Hwave,                                      &
# else
     &                      Uwave_rms,                                  &
# endif
     &                      Dwave, Pwave_bot,                           &
     &                      rho, u, v,                                  &
     &                      bottom,                                     &
     &                      Iconv,                                      &
     &                      Ubot, Vbot, Ur, Vr,                         &
     &                      bustrc, bvstrc,                             &
     &                      bustrw, bvstrw,                             &
     &                      bustrcwmax, bvstrcwmax,                     &
     &                      bustr, bvstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      integer, intent(inout) :: Iconv(LBi:,LBj:)

      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
#  if defined SG_CALC_UB
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
#  else
      real(r8), intent(in) :: Uwave_rms(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)

      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)

      real(r8), intent(out) :: Ubot(LBi:,LBj:)
      real(r8), intent(out) :: Vbot(LBi:,LBj:)
      real(r8), intent(out) :: Ur(LBi:,LBj:)
      real(r8), intent(out) :: Vr(LBi:,LBj:)
      real(r8), intent(out) :: bustrc(LBi:,LBj:)
      real(r8), intent(out) :: bvstrc(LBi:,LBj:)
      real(r8), intent(out) :: bustrw(LBi:,LBj:)
      real(r8), intent(out) :: bvstrw(LBi:,LBj:)
      real(r8), intent(out) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(out) :: bvstrcwmax(LBi:,LBj:)
      real(r8), intent(out) :: bustr(LBi:,LBj:)
      real(r8), intent(out) :: bvstr(LBi:,LBj:)
# else
      integer, intent(inout) :: Iconv(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ZoBot(LBi:UBi,LBj:UBj)
#  if defined SG_CALC_UB
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: Uwave_rms(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pwave_bot(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)

      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)

      real(r8), intent(out) :: Ubot(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Vbot(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Ur(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Vr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: bustrc(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: bvstrc(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: bustrw(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: bvstrw(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: bustrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: bvstrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: bvstr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      logical :: ITERATE

      integer :: Iter, i, j, k

      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8) :: Fwave_bot, Kb, Kbh, KboKb0, Kb0, Kdelta, Ustr
      real(r8) :: anglec, anglew
      real(r8) :: cff, cff1, cff2, cff3, og, fac, fac1, fac2
      real(r8) :: sg_ab, sg_abokb, sg_a1, sg_b1, sg_chi, sg_c1, sg_dd
      real(r8) :: sg_epsilon, sg_eta, sg_fofa, sg_fofb, sg_fofc, sg_fwm
      real(r8) :: sg_kbs, sg_lambda, sg_mu, sg_phicw, sg_ro, sg_row
      real(r8) :: sg_shdnrm, sg_shld, sg_shldcr, sg_scf, sg_ss, sg_star
      real(r8) :: sg_ub, sg_ubokur, sg_ubouc, sg_ubouwm, sg_ur
      real(r8) :: sg_ustarc, sg_ustarcw, sg_ustarwm, sg_znot, sg_znotp
      real(r8) :: sg_zr, sg_zrozn, sg_z1, sg_z1ozn, sg_z2, twopi, z1, z2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ab
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauc
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taucwmax
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ur_sg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vr_sg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ub
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ucur
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Umag
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vcur
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phic
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: phicw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rheight
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: rlength
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: u100
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: znot
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: znotc

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initalize to default values.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          Tauc(i,j)=0.0_r8
          Tauw(i,j)=0.0_r8
          Taucwmax(i,j)=0.0_r8
          u100(i,j)=0.0_r8
          rheight(i,j)=0.0_r8
          rlength(i,j)=0.0_r8
          znot(i,j)=ZoBot(i,j)
          znotc(i,j)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set currents above bed.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend+1
        DO i=IstrU-1,Iend+1
          Zr(i,j)=z_r(i,j,1)-z_w(i,j,0)
          Ur_sg(i,j)=u(i,j,1,nrhs)
          Vr_sg(i,j)=v(i,j,1,nrhs)
# ifdef SG_LOGINT
!
!  If current height is less than z1ur, interpolate logarithmically
!  to z1ur.
!
          IF (Zr(i,j).lt.sg_z1min) THEN
            DO k=2,N(ng)
              z1=z_r(i,j,k-1)-z_w(i,j,0)
              z2=z_r(i,j,k  )-z_w(i,j,0)
              IF ((z1.lt.sg_z1min).and.(sg_z1min.lt.z2)) THEN
                fac=1.0_r8/LOG(z2/z1)
                fac1=fac*LOG(z2/sg_z1min)
                fac2=fac*LOG(sg_z1min/z1)
                Ur_sg(i,j)=fac1*u(i,j,k-1,nrhs)+fac2*u(i,j,k,nrhs)
                Vr_sg(i,j)=fac1*v(i,j,k-1,nrhs)+fac2*v(i,j,k,nrhs)
                Zr(i,j)=sg_z1min
              END IF
            END DO
          END IF
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute bed wave orbital velocity (m/s) and excursion amplitude
!  (m) from wind-induced waves.  Use linear wave theory dispersion
!  relation for wave number.
!-----------------------------------------------------------------------
!
      twopi=2.0_r8*pi
      og=1.0_r8/g
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!
!  Compute first guess for wavenumber, Kb0.  Use deep water (Kb0*h>1)
!  and shallow water (Kb0*H<1) approximations.
!
          Fwave_bot=twopi/MAX(Pwave_bot(i,j),0.05_r8)
!
!  Compute bed wave orbital velocity and excursion amplitude.
!
# ifdef SG_CALC_UB
          Kb0=Fwave_bot*Fwave_bot*og
          IF (Kb0*h(i,j).ge.1.0_r8) THEN
            Kb=Kb0
          ELSE
            Kb=Fwave_bot/SQRT(g*h(i,j))
          END IF
!
!  Compute bottom wave number via Newton-Raphson method.
!
          ITERATE=.TRUE.
          DO Iter=1,sg_n
            IF (ITERATE) THEN
              Kbh=Kb*h(i,j)
              KboKb0=Kb/Kb0
              Kdelta=(1.0_r8-KboKb0*TANH(Kbh))/                         &
     &               (1.0_r8+Kbh*(KboKb0-1.0_r8/KboKb0))
              ITERATE=ABS(Kb*Kdelta) .ge. sg_tol
              Kb=Kb*(1.0_r8+Kdelta)
            END IF
          END DO
          Ab(i,j)=0.5_r8*Hwave(i,j)/SINH(Kb*h(i,j))+eps
          Ub(i,j)=Fwave_bot*Ab(i,j)+eps
# else
          Ub(i,j)=ABS(Uwave_rms(i,j))+eps
          Ab(i,j)=Ub(i,j)/Fwave_bot+eps
# endif
!
!  Compute bottom current magnitude at RHO-points.
!
          Ucur(i,j)=0.5_r8*(Ur_sg(i,j)+Ur_sg(i+1,j))
          Vcur(i,j)=0.5_r8*(Vr_sg(i,j)+Vr_sg(i,j+1))
          Umag(i,j)=SQRT(Ucur(i,j)*Ucur(i,j)+Vcur(i,j)*Vcur(i,j))+eps
!
!  Compute angle between currents and waves (radians)
!
          IF (Ucur(i,j).eq.0.0_r8) THEN
            phic(i,j)=0.5_r8*pi*SIGN(1.0_r8,Vcur(i,j))
          ELSE
            phic(i,j)=ATAN2(Vcur(i,j),Ucur(i,j))
          ENDIF
          phicw(i,j)=1.5_r8*pi-Dwave(i,j)-phic(i,j)-angler(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set default logarithmic profile.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          IF (Umag(i,j).gt.0.0_r8) THEN
!!          Ustr=MIN(sg_ustarcdef,Umag(i,j)*vonKar/                     &
!!   &                            LOG(Zr(i,j)/ZoBot(i,j)))
!!          Tauc(i,j)=Ustr*Ustr
            cff1=vonKar/LOG(Zr(i,j)/ZoBot(i,j))
            cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
            Tauc(i,j)=cff2*Umag(i,j)*Umag(i,j)
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Wave-current interaction case.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          sg_dd=bottom(i,j,isd50)
          sg_ss=bottom(i,j,idens)/(rho(i,j,1)+1000.0_r8)
          sg_ab=Ab(i,j)
          sg_ub=Ub(i,j)
          sg_phicw=phicw(i,j)
          sg_ur=Umag(i,j)
          sg_zr=Zr(i,j)
!
!  Compute hydraulic roughness "Znot" (m), ripple height "eta" (m),
!  and ripple length "lambda" (m).
!
# ifdef SG_CALC_ZNOT
          sg_star=sg_dd/(4.0_r8*sg_nu)*SQRT((sg_ss-1.0_r8)*sg_g*sg_dd)
!
!  Compute critical shield parameter based on grain diameter.
!  (sg_scf is a correction factor).
!
          sg_scf=1.0_r8
          IF (sg_star.le.1.5_r8) THEN
            sg_shldcr=sg_scf*0.0932_r8*sg_star**(-0.707_r8)
          ELSE IF ((1.5_r8.lt.sg_star).and.(sg_star.lt.4.0_r8)) THEN
            sg_shldcr=sg_scf*0.0848_r8*sg_star**(-0.473_r8)
          ELSE IF ((4.0_r8.le.sg_star).and.(sg_star.lt.10.0_r8)) THEN
            sg_shldcr=sg_scf*0.0680_r8*sg_star**(-0.314_r8)
          ELSE IF ((10.0_r8.le.sg_star).and.(sg_star.lt.34.0_r8)) THEN
            sg_shldcr=sg_scf*0.033_r8
          ELSE IF ((34.0_r8.le.sg_star).and.(sg_star.lt.270.0_r8)) THEN
            sg_shldcr=sg_scf*0.0134_r8*sg_star**(0.255_r8)
          ELSE
            sg_shldcr=sg_scf*0.056_r8
          END IF
!
!  Calculate skin friction shear stress based on Ole Madsen (1994)
!  empirical formula. Check initiation of sediment motion criteria,
!  to see if we compute sg_znot based on the wave-formed ripples.
!  If the skin friction calculation indicates that sediment is NOT
!  in motion, the ripple model is invalid and take the default value,
!  ZoBot.
!
          sg_abokb=sg_ab/sg_dd
          IF (sg_abokb.le.100.0_r8) THEN
            sg_fwm=EXP(7.02_r8*sg_abokb**(-0.078_r8)-8.82_r8)
          ELSE
            sg_fwm=EXP(5.61_r8*sg_abokb**(-0.109_r8)-7.30_r8)
          END IF
          sg_ustarwm=SQRT(0.5_r8*sg_fwm)*sg_ub
          sg_shdnrm=(sg_ss-1.0_r8)*sg_dd*sg_g
          sg_shld=sg_ustarwm*sg_ustarwm/sg_shdnrm
          IF ((sg_shld/sg_shldcr).le.1.0_r8) THEN
            sg_znot=ZoBot(i,j)
            sg_eta=0.0_r8
            sg_lambda=0.0_r8
          ELSE
!
!  Calculate ripple height and length and bottom roughness
!
            sg_chi=4.0_r8*sg_nu*sg_ub*sg_ub/                            &
     &             (sg_dd*((sg_ss-1.0_r8)*sg_g*sg_dd)**1.5_r8)
            IF (sg_chi.le.2.0_r8) THEN
              sg_eta=sg_ab*0.30_r8*sg_chi**(-0.39_r8)
              sg_lambda=sg_ab*1.96_r8*sg_chi**(-0.28_r8)
            ELSE
              sg_eta=sg_ab*0.45_r8*sg_chi**(-0.99_r8)
              sg_lambda=sg_ab*2.71_r8*sg_chi**(-0.75_r8)
            END IF
            sg_kbs=sg_ab*0.0655_r8*                                     &
     &             (sg_ub*sg_ub/((sg_ss-1.0_r8)*sg_g*sg_ab))**1.4_r8
            sg_znot=(sg_dd+2.3_r8*sg_eta+sg_kbs)/30.0_r8
          END IF
# else
          sg_znot=ZoBot(i,j)
          sg_chi=4.0_r8*sg_nu*sg_ub*sg_ub/                              &
     &           (sg_dd*((sg_ss-1.0_r8)*sg_g*sg_dd)**1.5_r8)
          IF (sg_chi.le.2.0_r8) THEN
            sg_eta=sg_ab*0.32_r8*sg_chi**(-0.34_r8)
            sg_lambda=sg_ab*2.04_r8*sg_chi**(-0.23_r8)
          ELSE
            sg_eta=sg_ab*0.52_r8*sg_chi**(-1.01_r8)
            sg_lambda=sg_ab*2.7_r8*sg_chi**(-0.78_r8)
          END IF
# endif
          znot(i,j)=sg_znot
          rheight(i,j)=sg_eta
          rlength(i,j)=sg_lambda
!
!  Compute only when nonzero currents and waves.
!
          sg_zrozn=sg_zr/sg_znot
          IF ((sg_ur.gt.0.0_r8).and.(sg_ub.gt.0.0_r8).and.              &
     &        (sg_zrozn.gt.1.0_r8)) THEN
!
! Compute bottom stress based on ripple roughness.
!
            sg_ubokur=sg_ub/(sg_kappa*sg_ur)
            sg_row=sg_ab/sg_znot
            sg_a1=1.0E-6_r8
            CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur,     &
     &                       sg_a1, sg_mu, sg_epsilon, sg_ro, sg_fofa)
            sg_abokb=sg_ab/(30.0_r8*sg_znot)
            IF (sg_abokb.le.100.0_r8) THEN
              sg_fwm=EXP(-8.82_r8+7.02_r8*sg_abokb**(-0.078_r8))
            ELSE
              sg_fwm=EXP(-7.30_r8+5.61_r8*sg_abokb**(-0.109_r8))
            END IF
            sg_ubouwm=SQRT(2.0_r8/sg_fwm)
!
!  Determine the maximum ratio of wave over combined shear stresses,
!  sg_ubouwm (ub/ustarwm).
!
            CALL sg_purewave (sg_row, sg_ubouwm, sg_znotp, sg_ro)
!
!  Set initial guess of the ratio of wave over shear stress, sg_c1
!  (ub/ustarc).
!
            sg_b1=sg_ubouwm
            sg_fofb=-sg_fofa
            sg_c1=0.5_r8*(sg_a1+sg_b1)
            CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur,     &
     &                       sg_c1, sg_mu, sg_epsilon, sg_ro, sg_fofc)
!
!  Solve PDE via bi-section method.
!
            ITERATE=.TRUE.
            DO Iter=1,sg_n
              IF (ITERATE) THEN
                IF ((sg_fofb*sg_fofc).lt.0.0_r8) THEN
                  sg_a1=sg_c1
                ELSE
                  sg_b1=sg_c1
                END IF
                sg_c1=0.5_r8*(sg_a1+sg_b1)
                CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur, &
     &                           sg_c1, sg_mu, sg_epsilon, sg_ro,       &
     &                           sg_fofc)
                ITERATE=(sg_b1-sg_c1) .ge. sg_tol
                IF (ITERATE) Iconv(i,j)=Iter
              END IF
            END DO
            sg_ubouc=sg_c1
!
!  Compute bottom shear stress magnitude (m/s).
!
            sg_ustarcw=sg_ub/sg_ubouc
            sg_ustarwm=sg_mu*sg_ustarcw
!!          sg_ustarc=MIN(sg_ustarcdef,sg_epsilon*sg_ustarcw)
!!          sg_ustarc=MIN(SQRT(Tauc(i,j)),sg_epsilon*sg_ustarcw)
            sg_ustarc=MAX(SQRT(Tauc(i,j)),sg_epsilon*sg_ustarcw)
            Tauc(i,j)=sg_ustarc*sg_ustarc
            Tauw(i,j)=sg_ustarwm*sg_ustarwm
            Taucwmax(i,j)=SQRT((Tauc(i,j)+                              &
     &                          Tauw(i,j)*COS(phicw(i,j)))**2+          &
     &                         (Tauw(i,j)*SIN(phicw(i,j)))**2)
!
!  Compute apparent hydraulic roughness (m).
!
            IF (sg_epsilon.gt.0.0_r8) THEN
              sg_z1=sg_alpha*sg_kappa*sg_ab/sg_ubouc
              sg_z2=sg_z1/sg_epsilon
              sg_z1ozn=sg_z1/sg_znot
              znotc(i,j)=sg_z2*                                         &
     &                   EXP(-(1.0_r8-sg_epsilon+                       &
     &                         sg_epsilon*LOG(sg_z1ozn)))
!
!  Compute mean (m/s) current at 100 cm above the bottom.
!
              IF (sg_z100.gt.sg_z2) THEN
                u100(i,j)=sg_ustarc*                                    &
     &                    (LOG(sg_z100/sg_z2)+1.0_r8-sg_epsilon+        &
     &                    sg_epsilon*LOG(sg_z1ozn))/sg_kappa
              ELSE IF ((sg_z100.le.sg_z2).and.(sg_zr.gt.sg_z1)) THEN
                u100(i,j)=sg_ustarc*sg_epsilon*                         &
     &                    (sg_z100/sg_z1-1.0_r8+LOG(sg_z1ozn))/sg_kappa
              ELSE
                u100(i,j)=sg_ustarc*sg_epsilon*                         &
     &                    LOG(sg_z100/sg_znot)/sg_kappa
              END IF
            END IF
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute kinematic bottom stress components due current and wind-
!  induced waves.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          anglec=Ur_sg(i,j)/(0.5*(Umag(i-1,j)+Umag(i,j)))
          bustr(i,j)=0.5_r8*(Tauc(i-1,j)+Tauc(i,j))*anglec
#  ifdef WET_DRY
          cff2=0.75_r8*0.5_r8*(z_w(i-1,j,1)+z_w(i,j,1)-                 &
     &                         z_w(i-1,j,0)-z_w(i,j,0))
          bustr(i,j)=SIGN(1.0_r8,bustr(i,j))*MIN(ABS(bustr(i,j)),       &
     &               ABS(u(i,j,1,nrhs))*cff2/dt(ng))
#  endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          anglec=Vr_sg(i,j)/(0.5_r8*(Umag(i,j-1)+Umag(i,j)))
          bvstr(i,j)=0.5_r8*(Tauc(i,j-1)+Tauc(i,j))*anglec
#  ifdef WET_DRY
          cff2=0.75_r8*0.5_r8*(z_w(i,j-1,1)+z_w(i,j,1)-                 &
     &                         z_w(i,j-1,0)-z_w(i,j,0))
          bvstr(i,j)=SIGN(1.0_r8,bvstr(i,j))*MIN(ABS(bvstr(i,j)),       &
     &               ABS(v(i,j,1,nrhs))*cff2/dt(ng))
#  endif
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          anglec=Ucur(i,j)/Umag(i,j)
          anglew=COS(1.5_r8*pi-Dwave(i,j)-angler(i,j))
          bustrc(i,j)=Tauc(i,j)*anglec
          bustrw(i,j)=Tauw(i,j)*anglew
          bustrcwmax(i,j)=Taucwmax(i,j)*anglew
          Ubot(i,j)=Ub(i,j)*anglew
          Ur(i,j)=Ucur(i,j)
!
          anglec=Vcur(i,j)/Umag(i,j)
          anglew=SIN(1.5_r8*pi-Dwave(i,j)-angler(i,j))
          bvstrc(i,j)=Tauc(i,j)*anglec
          bvstrw(i,j)=Tauw(i,j)*anglew
          bvstrcwmax(i,j)=Taucwmax(i,j)*anglew
          Vbot(i,j)=Ub(i,j)*anglew
          Vr(i,j)=Vcur(i,j)
!
          bottom(i,j,ibwav)=Ab(i,j)
          bottom(i,j,irhgt)=rheight(i,j)
          bottom(i,j,irlen)=rlength(i,j)
          bottom(i,j,izdef)=znot(i,j)
          bottom(i,j,izapp)=znotc(i,j)
        END DO
      END DO
!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustrc)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstrc)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustrw)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstrw)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustrcwmax)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstrcwmax)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Ubot)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Vbot)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Ur)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Vr)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,ibwav))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,irhgt))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,irlen))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izdef))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izapp))
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bustr, bvstr, bustrc, bvstrc)
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bustrw, bvstrw, bustrcwmax, bvstrcwmax)
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Ubot, Vbot, Ur, Vr)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom(:,:,ibwav),                            &
     &                    bottom(:,:,irhgt),                            &
     &                    bottom(:,:,irlen))
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom(:,:,izdef),                            &
     &                    bottom(:,:,izapp))
#endif

      RETURN
      END SUBROUTINE bblm_tile

      SUBROUTINE sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur,     &
     &                       sg_ubouc, sg_mu, sg_epsilon, sg_ro,        &
     &                       sg_fofx)
!
!=======================================================================
!                                                                      !
!  This routine computes bottom stresses via bottom boundary layer     !
!  formulation of Styles and Glenn (1999).                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     sg_row      Ratio of wave excursion amplitude over roughness.    !
!     sg_zrozn    Ratio of height of current over roughness.           !
!     sg_phiwc    Angle between wave and currents (radians).           !
!     sg_ubokur   Ratio of wave over current velocity:                 !
!                   ub/(vonKar*ur)                                     !
!     sg_ubouc    Ratio of bed wave orbital over bottom shear stress   !
!                   (ub/ustarc), first guess.                          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     sg_ubouc    Ratio of bed wave orbital over bottom shear stress   !
!                   (ub/ustarc), iterated value.                       !
!     sg_mu       Ratio between wave and current bottom shear          !
!                   stresses (ustarwm/ustarc).                         !
!     sg_epsilon  Ratio between combined (wave and current) and        !
!                   current bottom shear stresses (ustarc/ustarcw).    !
!     sg_ro       Internal friction Rossby number:                     !
!                   ustarc/(omega*znot)                                !
!     sg_fofx     Root of PDE used for convergence.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: sg_row, sg_zrozn, sg_phicw, sg_ubokur

      real(r8), intent(inout) :: sg_ubouc

      real(r8), intent(out) :: sg_mu, sg_epsilon, sg_ro, sg_fofx
!
!  Local variable declarations.
!
      logical :: ITERATE

      integer :: Iter

      real(r8) :: cff, sg_bei, sg_beip, sg_ber, sg_berp, sg_cosphi
      real(r8) :: sg_eps2, sg_kei, sg_keip, sg_ker, sg_kerp, sg_mu2
      real(r8) :: sg_phi, sg_ror, sg_x, sg_z2p, sg_znotp, sg_zroz1
      real(r8) :: sg_zroz2, sg_z1ozn, sg_z2ozn

      complex(c8) :: sg_argi, sg_bnot, sg_bnotp, sg_b1, sg_b1p
      complex(c8) :: sg_gammai, sg_knot, sg_knotp, sg_k1, sg_k1p
      complex(c8) :: sg_ll, sg_nn
!
!-----------------------------------------------------------------------
!  Compute bottom stresses.
!-----------------------------------------------------------------------
!
!  Compute nondimensional bottom wave shear, phi.  Iterate to make
!  sure that there is an upper limit in "ubouc".  It usually requires
!  only one pass.
!
      ITERATE=.TRUE.
      DO Iter=1,sg_n
        IF (ITERATE) THEN
          sg_ro=sg_row/sg_ubouc
          sg_znotp=1.0_r8/(sg_kappa*sg_ro)
          IF ((sg_z1p/sg_znotp).gt.1.0_r8) THEN
            sg_x=2.0_r8*SQRT(sg_znotp)
            IF (sg_x.le.8.0_r8) THEN
              CALL sg_kelvin8m (sg_x, sg_ber, sg_bei, sg_ker, sg_kei,   &
     &                          sg_berp, sg_beip, sg_kerp, sg_keip)
            ELSE
              CALL sg_kelvin8p (sg_x, sg_ker, sg_kei, sg_ber, sg_bei,   &
     &                          sg_kerp, sg_keip, sg_berp, sg_beip)
            END IF
            cff=1.0_r8/SQRT(sg_znotp)
            sg_bnot =CMPLX(sg_ber,sg_bei)
            sg_knot =CMPLX(sg_ker,sg_kei)
            sg_bnotp=CMPLX(sg_berp,sg_beip)*cff
            sg_knotp=CMPLX(sg_kerp,sg_keip)*cff
!
            sg_x=2.0_r8*SQRT(sg_z1p)
            IF (sg_x.le.8.0_r8) THEN
              CALL sg_kelvin8m (sg_x, sg_ber, sg_bei, sg_ker, sg_kei,   &
     &                          sg_berp, sg_beip, sg_kerp, sg_keip)
            ELSE
              CALL sg_kelvin8p (sg_x, sg_ker, sg_kei, sg_ber, sg_bei,   &
     &                          sg_kerp, sg_keip, sg_berp, sg_beip)
            END IF
            cff=1.0_r8/SQRT(sg_z1p)
            sg_b1 =CMPLX(sg_ber,sg_bei)
            sg_k1 =CMPLX(sg_ker,sg_kei)
            sg_b1p=CMPLX(sg_berp,sg_beip)*cff
            sg_k1p=CMPLX(sg_kerp,sg_keip)*cff
!
            sg_ll=sg_mp*sg_b1+sg_b1p
            sg_nn=sg_mp*sg_k1+sg_k1p
            sg_argi=sg_bnotp*sg_nn/(sg_bnot*sg_nn-sg_knot*sg_ll)+       &
     &              sg_knotp*sg_ll/(sg_knot*sg_ll-sg_bnot*sg_nn)
            sg_gammai=-sg_kappa*sg_znotp*sg_argi
            sg_phi=CABS(sg_gammai)
          ELSE
            sg_gammai=-sg_kappa*sg_z1p*sg_mp
            sg_phi=CABS(sg_gammai)
          END IF
!
          IF (sg_ubouc.gt.(1.0_r8/sg_phi)) THEN
            sg_ubouc=1.0_r8/sg_phi
          ELSE
            ITERATE=.FALSE.
          END IF
        END IF
      END DO
!
!  Compute ratio of wave over current bottom shear stresses.
!
      sg_mu=SQRT(sg_ubouc*sg_phi)
!
!  Compute ratio of current over combined bottom shear stresses.
!
      IF (sg_mu.eq.1.0_r8) THEN
        sg_epsilon=0.0_r8
      ELSE
        sg_mu2=sg_mu*sg_mu
        sg_cosphi=ABS(COS(sg_phicw))
        sg_eps2=-sg_mu2*sg_cosphi+                                      &
     &          SQRT(1.0_r8+sg_mu2*sg_mu2*(sg_cosphi*sg_cosphi-1.0_r8))
        sg_epsilon=SQRT(sg_eps2)
      END IF
!
!  Determine root of PDE used for convergence.
!
      IF (sg_epsilon.ne.0.0_r8) THEN
        sg_z2p=sg_z1p/sg_epsilon
        sg_ror=sg_ro/sg_zrozn
        sg_zroz1=1.0_r8/(sg_alpha*sg_kappa*sg_ror)
        sg_zroz2=sg_epsilon*sg_zroz1
        sg_z1ozn=sg_alpha*sg_kappa*sg_ro
        sg_z2ozn=sg_z1ozn/sg_epsilon
!
        IF ((sg_zroz2.gt.1.0_r8).and.(sg_z1ozn.gt.1.0_r8)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*                       &
     &                      (LOG(sg_zroz2)+1.0_r8-sg_epsilon+           &
     &                       sg_epsilon*LOG(sg_z1ozn))
        ELSE IF ((sg_zroz2.le.1.0_r8).and.(sg_zroz1.gt.1.0_r8).and.     &
     &          (sg_z1ozn.gt.1.0_r8)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*sg_epsilon*            &
     &                      (sg_zroz1-1.0_r8+LOG(sg_z1ozn))
        ELSE IF ((sg_zroz1.le.1.0_r8).and.(sg_z1ozn.gt.1.0_r8)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*sg_epsilon*            &
     &                      LOG(sg_zrozn)
        ELSE IF ((sg_zroz2.gt.1.0_r8).and.(sg_z1ozn.le.1.0_r8).and.     &
     &          (sg_z2ozn.gt.1.0_r8)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*                       &
     &                      (LOG(sg_zroz2)+1.0_r8-1.0_r8/sg_z2ozn)
        ELSE IF ((sg_zroz2.le.1.0_r8).and.(sg_zroz1.gt.1.0_r8).and.     &
     &          (sg_z1ozn.le.1.0_r8).and.(sg_z2ozn.gt.1.0_r8)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*sg_epsilon*            &
     &                      (sg_zroz1-1.0_r8/sg_z1ozn)
        ELSE IF ((sg_zroz2.gt.1.0_r8).and.(sg_z2ozn.le.1.0_r8)) THEN
          sg_fofx=-sg_ubouc+sg_ubokur*sg_epsilon*LOG(sg_zrozn)
        END IF
      END IF
      RETURN
      END SUBROUTINE sg_bstress

      SUBROUTINE  sg_purewave (sg_row, sg_ubouwm, sg_znotp, sg_ro)
!
!=======================================================================
!                                                                      !
!  This routine determines the maximum ratio of waves over combined    !
!  bottom shear stress.                                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     sg_row      Ratio of wave excursion amplitude over roughness.    !
!     sg_ubouwm   Maximum ratio of waves over combined bottom shear    !
!                   stress.                                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     sg_ubouwm   Maximum ratio of waves over combined bottom shear    !
!                   stress.                                            !
!     sg_znotp    Ratio of hydraulic roughness over scaled height      !
!                   of bottom boundary layer.                          !
!     sg_ro       Internal friction Rossby number:                     !
!                   ustarc/(omega*znot)                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: sg_row

      real(r8), intent(inout) :: sg_ubouwm

      real(r8), intent(out) :: sg_znotp, sg_ro
!
!  Local variable declarations.
!
      integer :: Iter

      real(r8) :: cff, sg_bei, sg_beip, sg_ber, sg_berp, sg_kei
      real(r8) :: sg_keip, sg_ker, sg_kerp, sg_phi, sg_ubouwmn, sg_x

      complex(c8) :: sg_argi, sg_bnot, sg_bnotp, sg_b1, sg_b1p
      complex(c8) :: sg_gammai, sg_knot, sg_knotp, sg_k1, sg_k1p
      complex(c8) :: sg_ll, sg_nn
!
!-----------------------------------------------------------------------
!  Compute wind-induced wave stress.
!-----------------------------------------------------------------------
!
      DO Iter=1,sg_n
        sg_ro=sg_row/sg_ubouwm
        sg_znotp=1.0_r8/(sg_kappa*sg_ro)
        IF (sg_z1p/sg_znotp.gt.1.0_r8) THEN
          sg_x=2.0_r8*SQRT(sg_znotp)
          IF (sg_x.le.8.0_r8) THEN
            CALL sg_kelvin8m (sg_x, sg_ber, sg_bei, sg_ker, sg_kei,     &
     &                        sg_berp, sg_beip, sg_kerp, sg_keip)
          ELSE
            CALL sg_kelvin8p (sg_x, sg_ker, sg_kei, sg_ber, sg_bei,     &
     &                        sg_kerp, sg_keip, sg_berp, sg_beip)
          END IF
          cff=1.0_r8/SQRT(sg_znotp)
          sg_bnot =CMPLX(sg_ber,sg_bei)
          sg_knot =CMPLX(sg_ker,sg_kei)
          sg_bnotp=CMPLX(sg_berp,sg_beip)*cff
          sg_knotp=CMPLX(sg_kerp,sg_keip)*cff
!
          sg_x=2.0*SQRT(sg_z1p)
          IF (sg_x.le.8.0_r8) THEN
            CALL sg_kelvin8m (sg_x, sg_ber, sg_bei, sg_ker, sg_kei,     &
     &                        sg_berp, sg_beip, sg_kerp, sg_keip)
          ELSE
            CALL sg_kelvin8p (sg_x, sg_ker, sg_kei, sg_ber, sg_bei,     &
     &                        sg_kerp, sg_keip, sg_berp, sg_beip)
          END IF
          cff=1.0_r8/SQRT(sg_z1p)
          sg_b1 =CMPLX(sg_ber,sg_bei)
          sg_k1 =CMPLX(sg_ker,sg_kei)
          sg_b1p=CMPLX(sg_berp,sg_beip)*cff
          sg_k1p=CMPLX(sg_kerp,sg_keip)*cff
!
          sg_ll=sg_mp*sg_b1+sg_b1p
          sg_nn=sg_mp*sg_k1+sg_k1p
          sg_argi=sg_bnotp*sg_nn/(sg_bnot*sg_nn-sg_knot*sg_ll)+         &
     &            sg_knotp*sg_ll/(sg_knot*sg_ll-sg_bnot*sg_nn)
          sg_gammai=-sg_kappa*sg_znotp*sg_argi
          sg_phi=CABS(sg_gammai)
        ELSE
          sg_gammai=-sg_kappa*sg_z1p*sg_mp
          sg_phi=CABS(sg_gammai)
        END IF
!
        sg_ubouwmn=1.0_r8/sg_phi
        IF (abs((sg_ubouwmn-sg_ubouwm)/sg_ubouwmn).le.sg_tol) THEN
          sg_ubouwm=sg_ubouwmn
          RETURN
        ELSE
          sg_ubouwm=sg_ubouwmn
        END IF
      END DO
      RETURN
      END SUBROUTINE  sg_purewave

      SUBROUTINE sg_kelvin8m (x, ber, bei, ker, kei, berp, beip,        &
     &                        kerp, keip)
!
!=======================================================================
!                                                                      !
! This rotuine computes the Kelvin functions for arguments less        !
! than eight.                                                          !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: x
      real(r8), intent(out) :: ber, bei, ker, kei
      real(r8), intent(out) :: berp, beip, kerp, keip
!
!  Local variable declarations.
!
      integer :: i

      real(r8) :: cff, xhalf

      real(r8), dimension(28) :: xp
!
!-----------------------------------------------------------------------
!  Compute Kelvin functions.
!-----------------------------------------------------------------------
!
      cff=0.125_r8*x
      xp(1)=cff
      DO i=2,28
        xp(i)=xp(i-1)*cff
      END DO
      xhalf=0.5_r8*x
!
      ber=1.0_r8-                                                       &
     &    64.0_r8*xp(4)+113.77777774_r8*xp(8)-                          &
     &    32.36345652_r8*xp(12)+2.64191397_r8*xp(16)-                   &
     &    0.08349609_r8*xp(20)+0.00122552_r8*xp(24)-                    &
     &    0.00000901_r8*xp(28)
      bei=16.0_r8*xp(2)-113.77777774_r8*xp(6)+                          &
     &    72.81777742*xp(10)-10.56765779_r8*xp(14)+                     &
     &    0.52185615_r8*xp(18)-0.01103667_r8*xp(22)+                    &
     &    0.00011346*xp(26)
!
      ker=-ber*LOG(xhalf)+0.25_r8*pi*bei-                               &
     &    0.57721566_r8-59.05819744*xp(4)+                              &
     &    171.36272133_r8*xp(8)-60.60977451_r8*xp(12)+                  &
     &    5.65539121_r8*xp(16)-0.19636347_r8*xp(20)+                    &
     &    0.00309699_r8*xp(24)-0.00002458_r8*xp(28)
      kei=-bei*LOG(xhalf)-0.25_r8*pi*ber+                               &
     &    6.76454936_r8*xp(2)-142.91827687_r8*xp(6)+                    &
     &    124.23569650_r8*xp(10)-21.30060904_r8*xp(14)+                 &
     &    1.17509064_r8*xp(18)-0.02695875_r8*xp(22)+                    &
     &    0.00029532_r8*xp(26)
!
      berp=x*(-4.0_r8*xp(2)+14.22222222_r8*xp(6)-                       &
     &        6.06814810_r8*xp(10)+0.66047849_r8*xp(14)-                &
     &        0.02609253_r8*xp(18)+0.00045957_r8*xp(22)-                &
     &        0.00000394_r8*xp(26))
      beip=x*(0.5_r8-10.66666666_r8*xp(4)+11.37777772_r8*xp(8)-         &
     &        2.31167514_r8*xp(12)+0.14677204_r8*xp(16)-                &
     &        0.00379386_r8*xp(20)+0.00004609_r8*xp(24))
!
      kerp=-berp*LOG(xhalf)-ber/x+0.25*pi*beip+                         &
     &     x*(-3.69113734_r8*xp(2)+21.42034017_r8*xp(6)-                &
     &        11.36433272_r8*xp(10)+1.41384780_r8*xp(14)-               &
     &        0.06136358_r8*xp(18)+0.00116137_r8*xp(22)-                &
     &        0.00001075*xp(26))
      keip=-beip*LOG(xhalf)-bei/x-0.25_r8*pi*berp+                      &
     &     x*(0.21139217_r8-13.39858846_r8*xp(4)+                       &
     &        19.41182758_r8*xp(8)-4.65950823_r8*xp(12)+                &
     &        0.33049424_r8*xp(16)-0.00926707_r8*xp(20)+                &
     &        0.00011997_r8*xp(24))
      RETURN
      END SUBROUTINE sg_kelvin8m

      SUBROUTINE sg_kelvin8p (x, ker, kei, ber, bei, kerp, keip,        &
     &                        berp, beip)
!
!=======================================================================
!                                                                      !
! This rotuine computes the Kelvin functions for arguments greater     !
! than eight.                                                          !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: x
      real(r8), intent(out) :: ker, kei, ber, bei
      real(r8), intent(out) :: kerp, keip, berp, beip
!
!  Local variable declarations.
!
      integer :: i

      real(r8) :: cff, xhalf

      real(r8), dimension(6) :: xm, xp

      complex(c8) :: argm, argp, fofx, gofx, phim, phip, thetam, thetap
!
!-----------------------------------------------------------------------
!  Compute Kelvin functions.
!-----------------------------------------------------------------------
!
      cff=8.0_r8/x
      xp(1)=cff
      xm(1)=-cff
      DO i=2,6
        xp(i)=xp(i-1)*cff
        xm(i)=-xm(i-1)*cff
      END DO
!
      thetap=CMPLX(0.0_r8,-0.3926991_r8)+                               &
     &       CMPLX(0.0110486_r8,-0.0110485_r8)*xp(1)+                   &
     &       CMPLX(0.0_r8,-0.0009765_r8)*xp(2)+                         &
     &       CMPLX(-0.0000906_r8,-0.0000901_r8)*xp(3)+                  &
     &       CMPLX(-0.0000252_r8,0.0_r8)*xp(4)+                         &
     &       CMPLX(-0.0000034_r8,0.0000051_r8)*xp(5)+                   &
     &       CMPLX(0.0000006,0.0000019)*xp(6)
      thetam=CMPLX(0.0_r8,-0.3926991_r8)+                               &
     &       CMPLX(0.0110486_r8,-0.0110485_r8)*xm(1)+                   &
     &       CMPLX(0.0_r8,-0.0009765_r8)*xm(2)+                         &
     &       CMPLX(-0.0000906_r8,-0.0000901_r8)*xm(3)+                  &
     &       CMPLX(-0.0000252_r8,0.0_r8)*xm(4)+                         &
     &       CMPLX(-0.0000034_r8,0.0000051_r8)*xm(5)+                   &
     &       CMPLX(0.0000006_r8,0.0000019_r8)*xm(6)
!
      phip=CMPLX(0.7071068_r8,0.7071068_r8)+                            &
     &     CMPLX(-0.0625001_r8,-0.0000001_r8)*xp(1)+                    &
     &     CMPLX(-0.0013813_r8,0.0013811_r8)*xp(2)+                     &
     &     CMPLX(0.0000005_r8,0.0002452_r8)*xp(3)+                      &
     &     CMPLX(0.0000346_r8,0.0000338_r8)*xp(4)+                      &
     &     CMPLX(0.0000117_r8,-0.0000024_r8)*xp(5)+                     &
     &     CMPLX(0.0000016_r8,-0.0000032_r8)*xp(6)
      phim=CMPLX(0.7071068_r8,0.7071068_r8)+                            &
     &     CMPLX(-0.0625001_r8,-0.0000001_r8)*xm(1)+                    &
     &     CMPLX(-0.0013813_r8,0.0013811_r8)*xm(2)+                     &
     &     CMPLX(0.0000005_r8,0.0002452_r8)*xm(3)+                      &
     &     CMPLX(0.0000346_r8,0.0000338_r8)*xm(4)+                      &
     &     CMPLX(0.0000117_r8,-0.0000024_r8)*xm(5)+                     &
     &     CMPLX(0.0000016_r8,-0.0000032_r8)*xm(6)
!
      cff=x/SQRT(2.0_r8)
      argm=-cff*CMPLX(1.0_r8,1.0_r8)+thetam
      fofx=SQRT(pi/(2.0_r8*x))*CEXP(argm)
      ker=REAL(fofx)
      kei=AIMAG(fofx)
!
      argp=cff*CMPLX(1.0_r8,1.0_r8)+thetap
      gofx=1.0_r8/SQRT(2.0_r8*pi*x)*CEXP(argp)
      ber=REAL(gofx)-kei/pi
      bei=AIMAG(gofx)+ker/pi
!
      kerp=REAL(-fofx*phim)
      keip=AIMAG(-fofx*phim)
!
      berp=REAL(gofx*phip)-keip/pi
      beip=AIMAG(gofx*phip)+kerp/pi
      RETURN
      END SUBROUTINE sg_kelvin8p
