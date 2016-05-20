#define M94WC
#undef SGWC
#define N92_RIPRUF

      SUBROUTINE bblm (ng, tile)
!
!svn $Id: ssw_bbl.h 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group        Chris Sherwood   !
!    Licensed under a MIT/X style license               Rich Signell   !
!    See License_ROMS.txt                             John C. Warner   !
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
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 37)
#endif
      CALL bblm_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                nrhs(ng),                                         &
     &                GRID(ng) % h,                                     &
     &                GRID(ng) % z_r,                                   &
     &                GRID(ng) % z_w,                                   &
     &                GRID(ng) % angler,                                &
     &                GRID(ng) % ZoBot,                                 &
#if defined SSW_CALC_UB
     &                FORCES(ng) % Hwave,                               &
#else
     &                FORCES(ng) % Uwave_rms,                           &
#endif
     &                FORCES(ng) % Dwave,                               &
     &                FORCES(ng) % Pwave_bot,                           &
#ifdef BEDLOAD
     &                SEDBED(ng) % bedldu,                              &
     &                SEDBED(ng) % bedldv,                              &
#endif
     &                SEDBED(ng) % bottom,                              &
     &                OCEAN(ng) % rho,                                  &
     &                OCEAN(ng) % u,                                    &
     &                OCEAN(ng) % v,                                    &
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
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 37)
#endif
      RETURN
      END SUBROUTINE bblm
!
!***********************************************************************
      SUBROUTINE bblm_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs,                                       &
     &                      h, z_r, z_w, angler, ZoBot,                 &
#if defined SSW_CALC_UB
     &                      Hwave,                                      &
#else
     &                      Uwave_rms,                                  &
#endif
     &                      Dwave, Pwave_bot,                           &
#ifdef BEDLOAD
     &                      bedldu, bedldv,                             &
#endif
     &                      bottom, rho, u, v,                          &
     &                      Iconv,                                      &
     &                      Ubot, Vbot, Ur, Vr,                         &
     &                      bustrc, bvstrc,                             &
     &                      bustrw, bvstrw,                             &
     &                      bustrcwmax, bvstrcwmax,                     &
     &                      bustr, bvstr)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_sediment
!
      USE bc_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs

#ifdef ASSUMED_SHAPE
      integer, intent(inout) :: Iconv(LBi:,LBj:)

      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
# if defined SSW_CALC_UB
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
# else
      real(r8), intent(in) :: Uwave_rms(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
# ifdef BEDLOAD
      real(r8), intent(in) :: bedldu(LBi:,LBj:,:)
      real(r8), intent(in) :: bedldv(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
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
#else
      integer, intent(inout) :: Iconv(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ZoBot(LBi:UBi,LBj:UBj)
# if defined SSW_CALC_UB
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: Uwave_rms(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pwave_bot(LBi:UBi,LBj:UBj)
# ifdef BEDLOAD
      real(r8), intent(in) :: bedldu(LBi:UBi,LBj:UBj,1:NST)
      real(r8), intent(in) :: bedldv(LBi:UBi,LBj:UBj,1:NST)
# endif
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
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
#endif
!
!  Local variable declarations.
!
      logical :: ITERATE

      integer :: Iter, i, j, k

      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8) :: Kbh, Kbh2, Kdh
      real(r8) :: taucr, wsedr, tstar, coef_st
      real(r8) :: coef_b1, coef_b2, coef_b3, d0
      real(r8) :: dolam, dolam1, doeta1, doeta2, fdo_etaano
      real(r8) :: lamorb, lamanorb
      real(r8) :: m_ubr, m_wr, m_ucr, m_zr, m_phicw, m_kb
      real(r8) :: m_ustrc, m_ustrwm, m_ustrr, m_fwc, m_zoa
      real(r8) :: zo
      real(r8) :: Kb, Kdelta, Ustr
      real(r8) :: anglec, anglew
      real(r8) :: cff, cff1, cff2, cff3, og, fac, fac1, fac2
      real(r8) :: sg_ab, sg_abokb, sg_a1, sg_b1, sg_chi, sg_c1, d50
      real(r8) :: sg_epsilon, ssw_eta, sg_fofa, sg_fofb, sg_fofc, sg_fwm
      real(r8) :: sg_kbs, ssw_lambda, sg_mu, sg_phicw, sg_ro, sg_row
      real(r8) :: sg_shdnrm, sg_shld, sg_shldcr, sg_scf, rhos, sg_star
      real(r8) :: sg_ub, sg_ubokur, sg_ubouc, sg_ubouwm, sg_ur
      real(r8) :: sg_ustarc, sg_ustarcw, sg_ustarwm, sg_znot, sg_znotp
      real(r8) :: sg_zr, sg_zrozn, sg_z1, sg_z1ozn, sg_z2, twopi, z1, z2
      real(r8) :: zoMIN, zoMAX
      real(r8) :: coef_fd

      real(r8), parameter :: absolute_zoMIN = 5.0d-5  ! in Harris-Wiberg
!!    real(r8), parameter :: absolute_zoMIN = 5.0d-8  ! in Harris-Wiberg
      real(r8), parameter ::  Cd_fd = 0.5_r8

      real(r8), parameter :: K1 = 0.6666666666_r8     ! Coefficients for
      real(r8), parameter :: K2 = 0.3555555555_r8     ! explicit
      real(r8), parameter :: K3 = 0.1608465608_r8     ! wavenumber
      real(r8), parameter :: K4 = 0.0632098765_r8     ! calculation
      real(r8), parameter :: K5 = 0.0217540484_r8     ! (Dean and
      real(r8), parameter :: K6 = 0.0065407983_r8     !  Dalrymple, 1991)

      real(r8), parameter :: coef_a1=0.095_r8         ! Coefficients for
      real(r8), parameter :: coef_a2=0.442_r8         ! ripple predictor
      real(r8), parameter :: coef_a3=2.280_r8         ! (Wiberg-Harris)

#if defined GM82_RIPRUF
      real(r8), parameter :: ar = 27.7_r8/30.0_r8  ! Grant-Madsen (1982)
#elif defined N92_RIPRUF
      real(r8), parameter :: ar = 0.267_r8         ! Nielsen (1992)
#elif defined R88_RIPRUF
      real(r8), parameter :: ar = 0.533_r8         ! Raudkivi (1988)
#else
      No ripple roughness coeff. chosen
#endif

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ab
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Fwave_bot
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
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoN
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoST
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoBF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoDEF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: zoBIO

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set currents above the bed.
!-----------------------------------------------------------------------
!
      twopi=2.0_r8*pi

      DO j=JstrV-1,Jend+1
        DO i=IstrU-1,Iend+1
          Zr(i,j)=z_r(i,j,1)-z_w(i,j,0)
          Ur_sg(i,j)=u(i,j,1,nrhs)
          Vr_sg(i,j)=v(i,j,1,nrhs)
#ifdef SSW_LOGINT
!
!  If current height is less than z1ur, interpolate logarithmically
!  to z1ur. (This has not been updated wrt ssw...uses outdated zo defs)
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
#endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute bottom stresses.
!-----------------------------------------------------------------------
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!
!  Set bed wave orbital velocity and excursion amplitude.  Use data
!  from wave models (SWAN) or use Dean and Dalrymple (1991) 6th-degree
!  polynomial to approximate wave number on shoaling water.

          Fwave_bot(i,j)=twopi/MAX(Pwave_bot(i,j),0.05_r8)
#ifdef SSW_CALC_UB
          Kdh=h(i,j)*Fwave_bot(i,j)**2/g
          Kbh2=Kdh*Kdh+                                                 &
     &         Kdh/(1.0_r8+Kdh*(K1+Kdh*(K2+Kdh*(K3+Kdh*(K4+             &
     &              Kdh*(K5+K6*Kdh))))))
          Kbh=SQRT(Kbh2)
          Ab(i,j)=0.5_r8*Hwave(i,j)/SINH(Kbh)+eps
          Ub(i,j)=Fwave_bot(i,j)*Ab(i,j)+eps
#else
          Ub(i,j)=MAX(Uwave_rms(i,j),0.0_r8)+eps
          Ab(i,j)=Ub(i,j)/Fwave_bot(i,j)+eps
#endif
!
!  Compute bottom current magnitude at RHO-points.
!
          Ucur(i,j)=0.5_r8*(Ur_sg(i,j)+Ur_sg(i+1,j))
          Vcur(i,j)=0.5_r8*(Vr_sg(i,j)+Vr_sg(i,j+1))
          Umag(i,j)=SQRT(Ucur(i,j)*Ucur(i,j)+Vcur(i,j)*Vcur(i,j)+eps)
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
!  Loop over RHO points.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!
!  Load sediment properties, stresses, and roughness from previous time
!  step (stresses are in m2 s-2).
!
          d50=bottom(i,j,isd50)
          rhos=bottom(i,j,idens)/(rho(i,j,1)+1000.0_r8)
          wsedr=bottom(i,j,iwsed)
          Taucr=bottom(i,j,itauc)
          Tauc(i,j)=SQRT(bustrc(i,j)**2+bvstrc(i,j)**2)
          Tauw(i,j)=SQRT(bustrw(i,j)**2+bvstrw(i,j)**2)
          Taucwmax(i,j)=SQRT( bustrcwmax(i,j)**2+bvstrcwmax(i,j)**2)
!
          rheight(i,j)=bottom(i,j,irhgt)
          rlength(i,j)=bottom(i,j,irlen)
          zoMAX=0.9_r8*Zr(i,j)
          zoMIN=MAX(absolute_zoMIN,2.5_r8*d50/30.0_r8)
!
!  Initialize arrays.
!
          zoN(i,j)=MIN(MAX(2.5_r8*d50/30.0_r8, zoMIN ),zoMAX)
          zoST(i,j)=0.0_r8
          zoBF(i,j)=0.0_r8
          zoBIO(i,j)=0.0_r8
!!        zoDEF(i,j)=bottom(i,j,izdef)
          zoDEF(i,j)=ZoBot(i,j)

#ifdef SSW_CALC_ZNOT
!
!  Calculate components of roughness and sum: zo = zoN + zoST + zoBF
!  Determine whether sediment is in motion. Use Shields criterion to
!  determine if sediment is mobile.
!
          tstar=Taucwmax(i,j)/(Taucr+eps)
          IF (tstar.lt.1.0_r8) THEN                         ! no motion
            zoST(i,j)=0.0_r8
            zoBF(i,j)=ar*rheight(i,j)**2/rlength(i,j)
          ELSE
!
!  Threshold of motion exceeded - calculate new zoST and zoBF
!  Calculate saltation roughness according to Wiberg & Rubin (1989)
!  (Eqn. 11 in Harris & Wiberg, 2001)
!  (d50 is in m, but this formula needs cm)
!
             coef_st=0.0204_r8*LOG(100.0_r8*d50+eps)**2+                &
     &               0.0220_r8*LOG(100.0_r8*d50+eps)+0.0709_r8
             zoST(i,j)=0.056_r8*d50*0.68_r8*tstar/                      &
     &                 (1.0_r8+coef_st*tstar)
             IF (zoST(i,j).lt.0.0_r8) THEN
               IF (Master) THEN
                 PRINT *, ' Warning: zoST<0  tstar, d50, coef_st:'
                 PRINT *, tstar,d50,coef_st
               END IF
             END IF
!
!  Calculate ripple height and wavelength.
!  Use Malarkey & Davies (2003) explict version of Wiberg & Harris.
!
             coef_b1=1.0_r8/coef_a1
             coef_b2=0.5_r8*(1.0_r8 + coef_a2)*coef_b1
             coef_b3=coef_b2**2-coef_a3*coef_b1
             d0=2.0_r8*Ab(i,j)
             IF ((d0/d50).gt.13000.0_r8) THEN              ! sheet flow
               rheight(i,j)=0.0_r8
               rlength(i,j)=535.0_r8*d50        ! does not matter since
             ELSE                               ! rheight=0
               dolam1=d0/(535.0_r8*d50)
               doeta1=EXP(coef_b2-SQRT(coef_b3-coef_b1*LOG(dolam1)))
               lamorb=0.62_r8*d0
               lamanorb=535.0_r8*d50
               IF (doeta1.lt.20.0_r8) THEN
                 dolam=1.0_r8/0.62_r8
               ELSE IF (doeta1.gt.100.0_r8) THEN
                 dolam=dolam1
               ELSE
                 fdo_etaano=-LOG(lamorb/lamanorb)*                      &
     &                       LOG(0.01_r8*doeta1)/LOG(5.0_r8)
                 dolam=dolam1*EXP(-fdo_etaano)
               END IF
               doeta2=EXP(coef_b2-SQRT(coef_b3-coef_b1*LOG(dolam)))
               rheight(i,j)=d0/doeta2
               rlength(i,j)=d0/dolam
             END IF
!
!  Value of ar can range from 0.3 to 3 (Soulsby, 1997, p. 124)
!
             zoBF(i,j)=ar*rheight(i,j)**2/rlength(i,j)
          END IF
          zo=zoN(i,j)
# ifdef SSW_ZOBL
          zo=zo+zoST(i,j)
# endif
# ifdef SSW_ZORIP
          zo=zo+zoBF(i,j)
# endif
# ifdef SSW_ZOBIO
          zo=zo+zoBIO(i,j)
# endif
#else
          IF (zoDEF(i,j).lt.absolute_zoMIN) THEN
            zoDEF(i,j)=absolute_zoMIN
            IF (Master) THEN
              PRINT *, ' Warning: default zo < 0.05 mm, replaced with:',&
     &                 zoDEF
            END IF
          END IF
          zo=zoDEF(i,j)
#endif
!
!  Compute stresses.
!
!  Default stress calcs for pure currents
!
          zo=MIN(MAX(zo,zoMIN),zoMAX)
!
          cff1=vonKar/LOG(Zr(i,j)/zo)
          cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
          Tauc(i,j)=cff2*Umag(i,j)*Umag(i,j)
          Tauw(i,j)=0.0_r8
          Taucwmax(i,j)=Tauc(i,j)
          znot(i,j)=zo
          znotc(i,j)=zo
!
          IF ((Umag(i,j).le.eps).and.(Ub(i,j).gt.eps)) THEN
!
!  Pure waves - use wave friction factor approach from Madsen
!  (1994, eqns 32-33).
!
            sg_abokb=Ab(i,j)/(30.0_r8*zo)
            sg_fwm=0.3_r8
            IF ((sg_abokb.gt.0.2_r8).and.(sg_abokb.le.100.0_r8)) THEN
              sg_fwm=EXP(-8.82_r8+7.02_r8*sg_abokb**(-0.078_r8))
            ELSE IF (sg_abokb.gt.100.0_r8)THEN
              sg_fwm=EXP(-7.30_r8+5.61_r8*sg_abokb**(-0.109_r8))
            END IF
            Tauc(i,j)= 0.0_r8
            Tauw(i,j)= 0.5_r8*sg_fwm*Ub(i,j)*Ub(i,j)
            Taucwmax(i,j)=Tauw(i,j)
            znot(i,j)=zo
            znotc(i,j)=zo
          ELSE IF ((Umag(i,j).gt.0.0_r8).and.(Ub(i,j).gt.eps).and.      &
     &             ((Zr(i,j)/zo).le.1.0_r8)) THEN
!
!  Waves and currents, but zr <= zo.
!
            IF (Master) THEN
              PRINT *,' Warning: w-c calcs ignored because zr <= zo'
            END IF
          ELSE IF ((Umag(i,j).gt.0.0_r8).and.(Ub(i,j).gt.eps).and.      &
     &             ((Zr(i,j)/zo).gt.1.0_r8)) THEN
!
!  Waves and currents, zr > zo.
!
#if defined SGWC
            sg_zrozn=Zr(i,j)/zo
            sg_ubokur=Ub(i,j)/(sg_kappa*Umag(i,j))
            sg_row=Ab(i,j)/zo
            sg_a1=1.0d-6
            sg_phicw=phicw(i,j)
            CALL sg_bstress (sg_row, sg_zrozn, sg_phicw, sg_ubokur,     &
     &                       sg_a1, sg_mu, sg_epsilon, sg_ro, sg_fofa)
            sg_abokb=Ab(i,j)/(30.0_r8*zo)
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
            ITERATE=.true.
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
            sg_ustarcw=Ub(i,j)/sg_ubouc
            sg_ustarwm=sg_mu*sg_ustarcw
!!          sg_ustarc=MIN(sg_ustarcdef,sg_epsilon*sg_ustarcw)  !original
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
              sg_z1=sg_alpha*sg_kappa*Ab(i,j)/sg_ubouc
              sg_z2=sg_z1/sg_epsilon
              sg_z1ozn=sg_z1/zo
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
              ELSE IF ((sg_z100.le.sg_z2).and.(Zr(i,j).gt.sg_z1)) THEN
                u100(i,j)=sg_ustarc*sg_epsilon*                         &
     &                    (sg_z100/sg_z1-1.0_r8+LOG(sg_z1ozn))/sg_kappa
              ELSE
                u100(i,j)=sg_ustarc*sg_epsilon*                         &
     &                    LOG(sg_z100/zo)/sg_kappa
              END IF
            END IF
#elif defined M94WC
            m_ubr=Ub(i,j)
            m_wr=Fwave_bot(i,j)
            m_ucr=Umag(i,j)
            m_zr=Zr(i,j)
            m_phicw=phicw(i,j)
            m_kb=30.0_r8*zo
            CALL madsen94 (m_ubr, m_wr, m_ucr,                          &
     &                     m_zr, m_phicw, m_kb,                         &
     &                     m_ustrc, m_ustrwm, m_ustrr, m_fwc, m_zoa)
            Tauc(i,j)=m_ustrc*m_ustrc
            Tauw(i,j)=m_ustrwm*m_ustrwm
            Taucwmax(i,j)=m_ustrr*m_ustrr
            znotc(i,j)=min( m_zoa, zoMAX )
            u100(i,j)=(m_ustrc/vonKar)*LOG(1.0_r8/m_zoa)
#endif
#if defined SSW_FORM_DRAG_COR
            IF (rheight(i,j).gt.(zoN(i,j)+zoST(i,j))) THEN
              coef_fd=0.5_r8*Cd_fd*(rheight(i,j)/rlength(i,j))*         &
     &                (1.0_r8/(vonKar*vonKar))*                         &
     &                (LOG(rheight(i,j)/                                &
     &                 (zoN(i,j)+zoST(i,j))-1.0_r8))**2
               Taucwmax(i,j)=Taucwmax(i,j)/(1.0_r8+coef_fd)
               Taucwmax(i,j)=Taucwmax(i,j)*(1.0_r8+8.0_r8*              &
     &                       rheight(i,j)/rlength(i,j))
            END IF
#endif
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
          bottom(i,j,izdef)=zoDEF(i,j)
          bottom(i,j,izapp)=znotc(i,j)
          bottom(i,j,izNik)=zoN(i,j)
          bottom(i,j,izbio)=zoBIO(i,j)
          bottom(i,j,izbfm)=zoBF(i,j)
          bottom(i,j,izbld)=zoST(i,j)
          bottom(i,j,izwbl)=znot(i,j)
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
      CALL bc_r2d_tile (ng, tile,                                       &
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
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izNik))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izbio))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izbfm))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izbld))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bottom(:,:,izwbl))
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
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom(:,:,izdef),                            &
     &                    bottom(:,:,izapp))
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom(:,:,izNik),                            &
     &                    bottom(:,:,izbio),                            &
     &                    bottom(:,:,izbfm),                            &
     &                    bottom(:,:,izwbl))
#endif

      RETURN
      END SUBROUTINE bblm_tile

#ifdef SGWC
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
      ITERATE=.true.
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
            sg_bnot =CMPLX(sg_ber,sg_bei,8)
            sg_knot =CMPLX(sg_ker,sg_kei,8)
            sg_bnotp=CMPLX(sg_berp,sg_beip,8)*cff
            sg_knotp=CMPLX(sg_kerp,sg_keip,8)*cff
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
            sg_b1 =CMPLX(sg_ber,sg_bei,8)
            sg_k1 =CMPLX(sg_ker,sg_kei,8)
            sg_b1p=CMPLX(sg_berp,sg_beip,8)*cff
            sg_k1p=CMPLX(sg_kerp,sg_keip,8)*cff
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
            ITERATE=.false.
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
          sg_bnot =CMPLX(sg_ber,sg_bei,8)
          sg_knot =CMPLX(sg_ker,sg_kei,8)
          sg_bnotp=CMPLX(sg_berp,sg_beip,8)*cff
          sg_knotp=CMPLX(sg_kerp,sg_keip,8)*cff
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
          sg_b1 =CMPLX(sg_ber,sg_bei,8)
          sg_k1 =CMPLX(sg_ker,sg_kei,8)
          sg_b1p=CMPLX(sg_berp,sg_beip,8)*cff
          sg_k1p=CMPLX(sg_kerp,sg_keip,8)*cff
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
! than eight (p 384 Abram and Stegun).                                 !
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
      thetap=CMPLX(0.0_r8,-0.3926991_r8,8)+                             &
     &       CMPLX(0.0110486_r8,-0.0110485_r8,8)*xp(1)+                 &
     &       CMPLX(0.0_r8,-0.0009765_r8,8)*xp(2)+                       &
     &       CMPLX(-0.0000906_r8,-0.0000901_r8,8)*xp(3)+                &
     &       CMPLX(-0.0000252_r8,0.0_r8,8)*xp(4)+                       &
     &       CMPLX(-0.0000034_r8,0.0000051_r8,8)*xp(5)+                 &
     &       CMPLX(0.0000006,0.0000019,8)*xp(6)
      thetam=CMPLX(0.0_r8,-0.3926991_r8,8)+                             &
     &       CMPLX(0.0110486_r8,-0.0110485_r8,8)*xm(1)+                 &
     &       CMPLX(0.0_r8,-0.0009765_r8,8)*xm(2)+                       &
     &       CMPLX(-0.0000906_r8,-0.0000901_r8,8)*xm(3)+                &
     &       CMPLX(-0.0000252_r8,0.0_r8,8)*xm(4)+                       &
     &       CMPLX(-0.0000034_r8,0.0000051_r8,8)*xm(5)+                 &
     &       CMPLX(0.0000006_r8,0.0000019_r8,8)*xm(6)
!
      phip=CMPLX(0.7071068_r8,0.7071068_r8,8)+                          &
     &     CMPLX(-0.0625001_r8,-0.0000001_r8,8)*xp(1)+                  &
     &     CMPLX(-0.0013813_r8,0.0013811_r8,8)*xp(2)+                   &
     &     CMPLX(0.0000005_r8,0.0002452_r8,8)*xp(3)+                    &
     &     CMPLX(0.0000346_r8,0.0000338_r8,8)*xp(4)+                    &
     &     CMPLX(0.0000117_r8,-0.0000024_r8,8)*xp(5)+                   &
     &     CMPLX(0.0000016_r8,-0.0000032_r8,8)*xp(6)
      phim=CMPLX(0.7071068_r8,0.7071068_r8,8)+                          &
     &     CMPLX(-0.0625001_r8,-0.0000001_r8,8)*xm(1)+                  &
     &     CMPLX(-0.0013813_r8,0.0013811_r8,8)*xm(2)+                   &
     &     CMPLX(0.0000005_r8,0.0002452_r8,8)*xm(3)+                    &
     &     CMPLX(0.0000346_r8,0.0000338_r8,8)*xm(4)+                    &
     &     CMPLX(0.0000117_r8,-0.0000024_r8,8)*xm(5)+                   &
     &     CMPLX(0.0000016_r8,-0.0000032_r8,8)*xm(6)
!
      cff=x/SQRT(2.0_r8)
      argm=-cff*CMPLX(1.0_r8,1.0_r8,8)+thetam
      fofx=SQRT(pi/(2.0_r8*x))*CEXP(argm)
      ker=REAL(fofx)
      kei=AIMAG(fofx)
!
      argp=cff*CMPLX(1.0_r8,1.0_r8,8)+thetap
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
#endif

#ifdef M94WC
      SUBROUTINE madsen94 (ubr, wr, ucr, zr, phiwc, kN,                 &
     &                     ustrc, ustrwm, ustrr, fwc, zoa)
!
!=======================================================================
!                                                                      !
!  Grant-Madsen model from Madsen (1994).                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ubr     Rep. wave-orbital velocity amplitude outside WBL (m/s).  !
!     wr      Rep. angular wave frequency,  2* pi/T (rad/s).           !
!     ucr     Current velocity at height zr (m/s).                     !
!     zr      Reference height for current velocity (m).               !
!     phiwc   Angle between currents and waves at zr (radians).        !
!     kN      Bottom roughness height, like Nikuradse k, (m).          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ustrc   Current friction velocity, u*c (m/s).                    !
!     ustrwm  Wave maximum friction velocity, u*wm (m/s).              !
!     ustrr   Wave-current combined friction velocity, u*r (m/s).      !
!     fwc     Wave friction factor (nondimensional).                   !
!     zoa     Apparent bottom roughness (m).                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      real(r8), intent(in) ::  ubr, wr, ucr, zr, phiwc, kN
      real(r8), intent(out) ::  ustrc, ustrwm, ustrr, fwc, zoa
!
!  Local variable declarations.
!
      integer, parameter :: MAXIT = 20
      integer :: i, nit
      integer :: iverbose = 1

      real(r8) :: bigsqr, cosphiwc, cukw, diff, lndw, lnln, lnzr
      real(r8) :: phicwc, zo
      real(r8) :: dval = 99.99_r8

      real(r8), dimension(MAXIT) :: Cmu
      real(r8), dimension(MAXIT) :: dwc
      real(r8), dimension(MAXIT) :: fwci
      real(r8), dimension(MAXIT) :: rmu
      real(r8), dimension(MAXIT) :: ustrci
      real(r8), dimension(MAXIT) :: ustrr2
      real(r8), dimension(MAXIT) :: ustrwm2
!
!-----------------------------------------------------------------------
!  Compute bottom friction velocities and roughness.
!-----------------------------------------------------------------------
!
!  Set special default values.
!
      ustrc=dval
      ustrwm=dval
      ustrr=dval
      fwc=0.4_r8
      zoa=kN/30.0_r8
      phicwc=phiwc

      zo = kN/30.0_r8

      IF (ubr.le.0.01_r8) THEN
        IF (ucr.le. 0.01_r8) THEN          ! no waves or currents
          ustrc=0.0_r8
          ustrwm=0.0_r8
          ustrr=0.0_r8
          RETURN
        END IF
        ustrc=ucr*vonKar/LOG(zr/zo)        ! no waves
        ustrwm=0.0_r8
        ustrr=ustrc
        RETURN
      END IF
!
!  Iterate to compute friction velocities, roughness, and wave friction
!  factor.  Notice that the computation of the wave friction factor
!  has been inlined for efficiency.
!
      cosphiwc=ABS(COS(phiwc))
      rmu(1)=0.0_r8
      Cmu(1)=1.0_r8

      cukw=Cmu(1)*ubr/(kN*wr)
      IF ((cukw.gt.0.2_r8).and.(cukw.le.100.0_r8)) THEN       ! Eq 32/33
        fwci(1)=Cmu(1)*EXP(7.02_r8*cukw**(-0.078_r8)-8.82_r8)
      ELSE IF ((cukw.gt.100.).and.(cukw.le.10000.0_r8)) THEN
        fwci(1)=Cmu(1)*EXP(5.61_r8*cukw**(-0.109_r8)-7.30_r8)
      ELSE IF (cukw.gt.10000.0_r8 ) THEN
        fwci(1)=Cmu(1)*EXP(5.61_r8*10000.0_r8**(-0.109_r8)-7.30_r8)
      ELSE
        fwci(1)=Cmu(1)*0.43_r8
      END IF
      ustrwm2(1)=0.5_r8*fwci(1)*ubr*ubr                       ! Eq 29
      ustrr2(1)=Cmu(1)*ustrwm2(1)                             ! Eq 26
      ustrr=SQRT(ustrr2(1))
      IF (cukw.ge.8.0_r8) THEN
        dwc(1)=2.0_r8*vonKar*ustrr/wr                         ! Eq 36
      ELSE
        dwc(1)=kN
      END IF
      lnzr=LOG(zr/dwc(1))
      lndw=LOG(dwc(1)/zo)
      lnln=lnzr/lndw
      bigsqr=-1.0_r8+SQRT(1.0_r8+((4.0_r8*vonKar*lndw)/                 &
     &                            (lnzr*lnzr))*ucr/ustrr)
      ustrci(1)=0.5_r8*ustrr*lnln*bigsqr
!
      i=1
      diff=1.0_r8
      DO WHILE ((i.lt.MAXIT).and.(diff.gt.0.000005_r8))
        i=i+1
        rmu(i)=ustrci(i-1)*ustrci(i-1)/ustrwm2(i-1)
        Cmu(i)=SQRT(1.0_r8+                                             &
     &              2.0_r8*rmu(i)*cosphiwc+rmu(i)*rmu(i))     ! Eq 27
        cukw=Cmu(i)*ubr/(kN*wr)
        IF ((cukw.gt.0.2_r8).and.(cukw.le.100.0_r8)) THEN     ! Eq 32/33
          fwci(i)=Cmu(i)*EXP(7.02_r8*cukw**(-0.078_r8)-8.82_r8)
        ELSE IF ((cukw.gt.100.).and.(cukw.le.10000.0_r8)) THEN
          fwci(i)=Cmu(i)*EXP(5.61_r8*cukw**(-0.109_r8)-7.30_r8)
        ELSE IF (cukw.gt.10000.0_r8 ) THEN
          fwci(i)=Cmu(i)*EXP(5.61_r8*10000.0_r8**(-0.109_r8)-7.30_r8)
        ELSE
          fwci(i)=Cmu(i)*0.43_r8
        END IF
        ustrwm2(i)=0.5_r8*fwci(i)*ubr*ubr                     ! Eq 29
        ustrr2(i)=Cmu(i)*ustrwm2(i)                           ! Eq 26
        ustrr=SQRT(ustrr2(i))
!!      IF ((Cmu(1)*ubr/(kN*wr)).ge.8.0_r8) THEN  ! HGA Why 1?
        IF (cukw.ge.8.0_r8) THEN
          dwc(i)=2.0_r8*vonKar*ustrr/wr                       ! Eq 36
        ELSE
          dwc(i)=kN
        END IF
        lnzr=LOG(zr/dwc(i))
        lndw=LOG(dwc(i)/zo)
        lnln=lnzr/lndw
        bigsqr=-1.0_r8+SQRT(1.0_r8+((4.0_r8*vonKar*lndw)/               &
     &                              (lnzr*lnzr))*ucr/ustrr)
        ustrci(i)=0.5_r8*ustrr*lnln*bigsqr                    ! Eq 38
        diff=ABS((fwci(i)-fwci(i-1))/fwci(i))
      END DO
      ustrwm=SQRT(ustrwm2(i))
      ustrc=ustrci(i)
      ustrr=SQRT(ustrr2(i))
      phicwc=phiwc
      zoa=EXP(LOG(dwc(i))-(ustrc/ustrr)*LOG(dwc(i)/zo))       ! Eq 11
      fwc=fwci(i)

      RETURN
      END SUBROUTINE madsen94
#endif
