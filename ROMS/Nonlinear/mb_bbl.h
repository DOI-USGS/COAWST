      SUBROUTINE bblm (ng,tile)
!
!svn $Id: mb_bbl.h 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group          Meinte Blaas   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine computes bottom stresses for combined waves and     !
!  currents using the parametric approximation by Soulsby 1997:        !
!                                                                      !
!        tauCW=tauC*[1+1.2*(tauW/(tauC+tauW))^3.2]                     !
!                                                                      !
!  and                                                                 !
!                                                                      !
!        tauCWmax=SQRT([tauCW+tauW*COS(phiCW)]^2+[tauW*SIN(phiCW)]^2)  !
!                                                                      !
!  where                                                               !
!                                                                      !
!  tauCW      Combined wave-averaged stress (in current direction).    !
!  tauC       Stress due to currents, if waves would be absent.        !
!  tauW       Amplitude of stress due to waves without currents.       !
!  tauCWmax   Maximum combined wave-averaged stress.                   !
!  phiCW      Angle between current and waves.                         !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Dyer 1986, Coastal and Estuarine Sediment Dynamics, Wiley,        !
!      342 pp.                                                         !
!                                                                      !
!    Harris and Wiberg 2001, Comp. and Geosci. 27, 675-690.            !
!                                                                      !
!    Li and Amos 2001, Comp. and Geosci. 27, 619-645.                  !
!                                                                      !
!    Soulsby 1997, Dynamics of Marine Sands, Telford  Publ., 249 pp.   !
!                                                                      !
!    Soulsby 1995, Bed shear-stresses due to combined waves and        !
!      currents, in: Stive et al: Advances in Coastal Morphodynamics,  !
!      Wiley, 4.20-4.23                                                !
!                                                                      !
!    Wiberg and Harris 1994, J. Geophys. Res. 99(C4), 775-789.         !
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
# ifdef MB_CALC_UB
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
# ifdef MB_CALC_UB
     &                      Hwave,                                      &
# else
     &                      Uwave_rms,                                  &
# endif
     &                      Dwave, Pwave_bot,                           &
     &                      rho, u, v, bottom,                          &
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
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
#  ifdef MB_CALC_UB
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
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ZoBot(LBi:UBi,LBj:UBj)
#  ifdef MB_CALC_UB
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
      integer :: i, ised, j

      real(r8), parameter :: K1 = 0.6666666666_r8
      real(r8), parameter :: K2 = 0.3555555555_r8
      real(r8), parameter :: K3 = 0.1608465608_r8
      real(r8), parameter :: K4 = 0.0632098765_r8
      real(r8), parameter :: K5 = 0.0217540484_r8
      real(r8), parameter :: K6 = 0.0065407983_r8

      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8), parameter :: scf1 = 0.5_r8 * 1.39_r8
      real(r8), parameter :: scf2 = 0.52_r8
      real(r8), parameter :: scf3 = 2.0_r8 - scf2
      real(r8), parameter :: scf4 = 1.2_r8
      real(r8), parameter :: scf5 = 3.2_r8

      real(r8) :: twopi = 2.0_r8 * pi

      real(r8) :: RHbiomax = 0.006_r8    ! maximum biogenic ripple height
      real(r8) :: RHmin = 0.001_r8       ! arbitrary minimum ripple height
      real(r8) :: RLmin = 0.01_r8        ! arbitrary minimum ripple length

      real(r8) :: Ab, Fwave_bot, Kbh, Kbh2, Kdh
      real(r8) :: angleC, angleW, phiC, phiCW
      real(r8) :: cff, cff1, cff2, d50, viscosity, wset
      real(r8) :: RHbio, RHbiofac, RHmax, RLbio
      real(r8) :: rhoWater, rhoSed
      real(r8) :: rhgt, rlen, thetw
      real(r8) :: Znot, ZnotC, Znot_bl, Znot_rip
      real(r8) :: tau_bf, tau_ex,  tau_en, tau_up, tau_w, tau_wb
      real(r8) :: tau_c, tau_cb, tau_cs, tau_cw, tau_cwb, tau_cws

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ub
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ucur
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vcur
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Ur_mb
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vr_mb
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Umag
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tauC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tauW
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tauCW
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tauCWmax

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initalize stresses due to currents and waves.
!-----------------------------------------------------------------------
!
      RHbiofac=1.0_r8/EXP(4.11_r8)
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          tauC(i,j)=0.0_r8
          tauCW(i,j)=0.0_r8
          tauCWmax(i,j)=0.0_r8
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
          Ur_mb(i,j)=u(i,j,1,nrhs)
          Vr_mb(i,j)=v(i,j,1,nrhs)
        END DO
      END DO
!
!=======================================================================
!  Compute bottom stresses.
!=======================================================================
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          RHbio=0.0_r8
          RLbio=0.0_r8
          rlen=bottom(i,j,irlen)
          rhgt=bottom(i,j,irhgt)
          rhoWater=rho(i,j,1)+1000.0_r8
          viscosity=0.0013_r8/rhoWater         ! kinematic viscosity
!
!-----------------------------------------------------------------------
!  Compute bed wave orbital velocity (m/s) and excursion amplitude (m)
!  from wind-induced waves.  Use Dean and Dalrymple (1991) 6th-degree
!  polynomial to approximate wave number on shoaling water.
!-----------------------------------------------------------------------
!
          Fwave_bot=twopi/MAX(Pwave_bot(i,j),0.05_r8)
# ifdef MB_CALC_UB
          Kdh=h(i,j)*Fwave_bot*Fwave_bot/g
          Kbh2=Kdh*Kdh+                                                 &
     &         Kdh/(1.0_r8+Kdh*(K1+Kdh*(K2+Kdh*(K3+Kdh*(K4+             &
     &              Kdh*(K5+K6*Kdh))))))
          Kbh=SQRT(Kbh2)
!
!  Compute bed wave orbital velocity and excursion amplitude.
!
          Ab=0.5_r8*Hwave(i,j)/SINH(Kbh)+eps
          Ub(i,j)=Fwave_bot*Ab
# else
          Ub(i,j)=Uwave_rms(i,j)
          Ab=Ub(i,j)/Fwave_bot+eps
# endif
!
!  Compute bottom current magnitude at RHO-points.
!
          Ucur(i,j)=0.5_r8*(Ur_mb(i,j)+Ur_mb(i+1,j))
          Vcur(i,j)=0.5_r8*(Vr_mb(i,j)+Vr_mb(i,j+1))
          Umag(i,j)=SQRT(Ucur(i,j)*Ucur(i,j)+Vcur(i,j)*Vcur(i,j))+eps
!
!  Compute angle between currents and waves (radians)
!
          IF (Ucur(i,j).ne.0.0_r8) THEN
            phiC=ATAN2(Vcur(i,j),Ucur(i,j))
          ELSE
            phiC=0.5_r8*pi*SIGN(1.0_r8,Vcur(i,j))
          END IF
          phiCW=1.5_r8*pi-Dwave(i,j)-phiC-angler(i,j)
!
!-----------------------------------------------------------------------
!  Determine skin roughness from sediment size. Set default logarithmic
!  profile consistent with current-only case.
!  Establish local median grain size for all calculations here and
!  determine local values of critical stresses.
!  Since most parameterizations have been derived ignoring multiple
!  grain sizes, we apply this single d50 also in the case of mixed beds.
!-----------------------------------------------------------------------
!
          d50=bottom(i,j,isd50)              ! (m)
          tau_cb=bottom(i,j,itauc)           ! (m2/s2)
          wset=bottom(i,j,iwsed)             ! (m/s)
          rhoSed=bottom(i,j,idens)/rhoWater  ! (nondimensional)
!
! Critical stress (m2/s2) for transition to sheet flow
! (Li and Amos, 2001, eq. 11).
!
          tau_up=0.172_r8*(rhoSed-1.0_r8)*g*(d50**0.624_r8)
!
! Threshold stress (m2/s2) for break off (Grant and Madsen,1982).
!
          tau_bf=0.79_r8*(viscosity**(-0.6_r8))*                        &
     &           (((rhoSed-1.0_r8)*g)**0.3_r8)*(d50**0.9_r8)*tau_cb
!
!  Set Znot for currents as maximum of user value or grain roughness.
!
           ZnotC=d50/12.0_r8
           Znot=MAX(ZoBot(i,j),ZnotC)
!
!-----------------------------------------------------------------------
!  Set tauC (m2/s2) acc. to current-only case (skin friction) [m/s]
!-----------------------------------------------------------------------
!
          cff1=vonKar/LOG(Zr(i,j)/Znot)
          cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
          tauC(i,j)=cff2*Umag(i,j)*Umag(i,j)

          cff1=vonKar/LOG(Zr(i,j)/ZnotC)
          tau_cs=cff1*cff1*Umag(i,j)*Umag(i,j)
!
!-----------------------------------------------------------------------
!  If significant waves (Ub > 0.01 m/s), wave-current interaction case
!  according to Soulsby (1995). Otherwise, tauCWmax = tauC for sediment
!  purposes.
!-----------------------------------------------------------------------
!
          IF (Ub(i,j).gt.0.01_r8) THEN
!
!  Determine skin stresses (m2/s2) for pure waves and combined flow
!  using Soulsby (1995) approximation of the wave friction factor,
!  fw=2*scf1*(Znot/Ab)**scf2 so tauW=0.5fw*Ub**2.
!
            tau_w=scf1*((ZnotC*Fwave_bot)**scf2)*(Ub(i,j)**scf3)
!
!  Wave-averaged, combined wave-current stress.(Eqn. 69, Soulsby, 1997).
!
            tau_cw=tau_cs*                                              &
     &             (1.0_r8+scf4*((tau_w/(tau_w+tau_cs))**scf5))
!
!  Maximum of combined wave-current skin stress (m2/s2) component for
!  sediment.(Eqn. 70, Soulsby, 1997).
!
            tau_cws=SQRT((tau_cw+tau_w*COS(phiCW))**2+                  &
     &                   (tau_w*SIN(phiCW))**2)
!!          tauw(i,j)=tau_cws
            tauCWmax(i,j)=tau_cws
            tauW(i,j)=tau_w
!
!  Set combined stress for Znot.
!
            tau_w=scf1*((Znot*Fwave_bot)**scf2)*(Ub(i,j)**scf3)
            tau_cw=tauC(i,j)*                                           &
     &             (1.0_r8+scf4*((tau_w/(tau_w+tauC(i,j)))**scf5))

# ifdef MB_Z0BL
!
!-----------------------------------------------------------------------
!  Compute bedload roughness for ripple predictor and sediment purposes.
!  At high transport stages, friction depends on thickness of bedload
!  layer. Shear stress due to combined grain and bedload roughness is
!  used to predict ripples and onset of suspension (Li and Amos, 2001).
!-----------------------------------------------------------------------
!
#  ifdef MB_CALC_ZNOT
            tau_ex=MAX((tau_cws-tau_cb),0.0_r8)
            cff=(1.0_r8/((rhoSed-1.0_r8)*g*d50))
            Znot_bl=17.4_r8*d50*(cff*tau_ex)**0.75_r8
            ZnotC=ZnotC+Znot_bl
#  endif
!
!-----------------------------------------------------------------------
!  Compute stresses (m2/s2)for sediment purposes, using grain and
!  bedload roughness.
!-----------------------------------------------------------------------
!
            cff1=vonKar/LOG(Zr(i,j)/ZnotC)
            tau_c=cff1*cff1*Umag(i,j)*Umag(i,j)
            tau_wb=scf1*((ZnotC*Fwave_bot)**scf2)*(Ub(i,j)**scf3)
            tau_cw=tau_c*(1.0_r8+scf4*((tau_wb/(tau_wb+tau_c))**scf5))
!
!  Maximum of combined wave-current stress (m2/s2) component for
!  sediment purposes.
!
            tau_cwb=SQRT((tau_cw+tau_wb*COS(phiCW))**2+                 &
     &                   (tau_wb*SIN(phiCW))**2)
            tauCWmax(i,j)=tau_cwb
            tauW(i,j)=tau_wb
# endif

# ifdef MB_Z0RIP
!
!-----------------------------------------------------------------------
!  Determine bedform roughness ripple height (m) and ripple length (m)
!  for sandy beds. Use structure according to Li and Amos (2001).
!-----------------------------------------------------------------------
!
!  Check median grain diameter
!
            IF (d50.ge.0.000063_r8) THEN
!
!  Enhanced skin stress if pre-exisiting ripples (Nielsen, 1986).
!
              RHmax=0.8_r8*rlen/pi
              rhgt=MAX(MIN(RHmax,rhgt),RHmin)
              tau_en=MAX(tau_cws,tau_cws*(rlen/(rlen-pi*rhgt))**2)
!
              IF ((tau_cws.lt.tau_cb).and.(tau_en.ge.tau_cb)) THEN
                 rhgt=(19.6_r8*(SQRT(tau_cws/tau_cb))+20.9_r8)*d50
                 rlen=rhgt/0.12_r8        ! local transport
              ELSE
                IF ((tau_cws.ge.tau_cb).and.(tau_cwb.lt.tau_bf)) THEN
                  rhgt=(22.15_r8*(SQRT(tau_cwb/tau_cb))+6.38_r8)*d50
                  rlen=rhgt/0.12_r8       ! bed load in eq. range
                ELSE
                  IF ((tau_cwb.ge.tau_bf).and.(tau_cwb.lt.tau_up)) THEN
                    rlen=535.0_r8*d50     ! break off regime
                    rhgt=0.15_r8*rlen*(SQRT(tau_up)-SQRT(tau_cwb))/     &
     &                                (SQRT(tau_up)-SQRT(tau_bf ))
                  ELSE IF (tau_cwb.ge.tau_up) THEN
                    rlen=0.0_r8           ! sheet flow, plane bed
                    rhgt=0.0_r8
                  ELSE
                    rlen=bottom(i,j,irlen)
                    rhgt=bottom(i,j,irhgt)
                  END IF
                END IF
              END IF
            END IF
# endif

# ifdef MB_Z0BIO
!
!-----------------------------------------------------------------------
!  Determine (biogenic) bedform roughness ripple height (m) and ripple
!  length (m) for silty beds, using Harris and Wiberg (2001).
!-----------------------------------------------------------------------
!
!  Use 10 cm default biogenic ripple length, RLbio (Wheatcroft 1994).
!
!  NOTE:  For mixed beds take average of silty and sandy bedform
!         roughnesses weighted according to silt and sand fractions.
!
            IF (d50.lt.0.000063_r8) THEN
              RLbio=0.1_r8
              thetw=tau_cws*(1.0_r8/((rhoSed-1.0_r8)*g*d50))
              RHbio=(thetw**(-1.67_r8))*RLbio*RHbiofac
              rhgt=MIN(RHbio,RHbiomax)
              rlen=RLbio
            END IF
# endif

# if defined MB_Z0RIP || defined MB_Z0BIO
!
!  Ripple roughness using Grant and Madsen (1982) roughness length.
!
#  ifdef MB_CALC_ZNOT
            Znot_rip=0.92_r8*rhgt*rhgt/(MAX(rlen,RLmin))
            ZnotC=ZnotC+Znot_rip
#  endif
!
!-----------------------------------------------------------------------
! Compute bottom stress (m2/s2) components based on total roughnes.
!-----------------------------------------------------------------------
!
            cff1=vonKar/LOG(Zr(i,j)/ZnotC)
            tau_c=cff1*cff1*Umag(i,j)*Umag(i,j)
            tau_w=scf1*((ZnotC*Fwave_bot)**scf2)*(Ub(i,j)**scf3)
            tau_cw=tau_c*                                               &
     &             (1.0_r8+scf4*((tau_w/(tau_w+tau_c))**scf5))
!!          tau_cwb=SQRT((tau_cw+tau_w*COS(phiCW))**2+                  &
!!   &                   (tau_w*SIN(phiCW))**2)
!!          tauCWmax(i,j)=tau_cwb
# endif
!
!-----------------------------------------------------------------------
!  Compute effective bottom shear velocity (m/s) relevant for flow and
!  eddy-diffusivities/viscosity.
!-----------------------------------------------------------------------
!
!!          tauC(i,j)=tau_cw
            tauCW(i,j)=tau_cw
            tauW(i,j)=tau_w
          ELSE IF (Ub(i,j).le.0.01_r8) THEN
!
!  If current-only, tauCWmax=tauC(skin) for use in sediment model; tauC
!  is determined using roughness due to current ripples.
!  In this limit, tauCWmax=tauCW=tauC
!
            tauCWmax(i,j)=tauC(i,j)
            tauW(i,j)=0.0_r8
!!          tauW(i,j)=tau_cs
!!          tauCW(i,j)=tau_cs

# ifdef MB_Z0RIP
!!          IF (tauC(i,j).gt.tau_up) THEN
            IF (tau_cs(i,j).gt.tau_up) THEN
              rhgt=0.0_r8
              rlen=0.0_r8
!!          ELSE IF (tauC(i,j).lt.tau_cb) THEN
            ELSE IF (tau_cs(i,j).lt.tau_cb) THEN
              rlen=bottom(i,j,irlen)
              rhgt=bottom(i,j,irhgt)
            ELSE
              rlen=1000.0_r8*d50                       ! Yalin (1964)
              rhgt=0.0308_r8*(rlen**1.19_r8)
!!            rhgt=100.0_r8*0.074_r8*                                   &
!!   &             (0.01_r8*rlen)**1.19_r8             ! Allen (1970)
            END IF
#  ifdef MB_CALC_ZNOT
            ZnotC=ZnotC+0.92_r8*rhgt*rhgt/(MAX(rlen,RLmin))
#  endif
# endif
            cff1=vonKar/LOG(Zr(i,j)/ZnotC)
            cff2=MIN(Cdb_max,MAX(Cdb_min,cff1*cff1))
!!          tauc(i,j)=cff2*Umag(i,j)*Umag(i,j)
            tauCW(i,j)=cff2*Umag(i,j)*Umag(i,j)
          END IF
!
!-----------------------------------------------------------------------
!  Load variables for output purposes.
!-----------------------------------------------------------------------
!
          bottom(i,j,ibwav)=Ab
          bottom(i,j,irlen)=rlen
          bottom(i,j,irhgt)=rhgt
          bottom(i,j,izdef)=Znot                       ! Zob(ng)
          bottom(i,j,izapp)=ZnotC
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute kinematic bottom stress (m2/s2) components for flow due to
!  combined current and wind-induced waves.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          angleC=Ur_mb(i,j)/(0.5*(Umag(i-1,j)+Umag(i,j)))
          bustr(i,j)=0.5_r8*(tauCW(i-1,j)+tauCW(i,j))*angleC
# ifdef WET_DRY
          cff2=0.75_r8*0.5_r8*(z_w(i-1,j,1)+z_w(i,j,1)-                 &
     &                         z_w(i-1,j,0)-z_w(i,j,0))
          bustr(i,j)=SIGN(1.0_r8,bustr(i,j))*MIN(ABS(bustr(i,j)),       &
     &               ABS(u(i,j,1,nrhs))*cff2/dt(ng))
# endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          angleC=Vr_mb(i,j)/(0.5_r8*(Umag(i,j-1)+Umag(i,j)))
          bvstr(i,j)=0.5_r8*(tauCW(i,j-1)+tauCW(i,j))*angleC
# ifdef WET_DRY
          cff2=0.75_r8*0.5_r8*(z_w(i,j-1,1)+z_w(i,j,1)-                 &
     &                         z_w(i,j-1,0)-z_w(i,j,0))
          bvstr(i,j)=SIGN(1.0_r8,bvstr(i,j))*MIN(ABS(bvstr(i,j)),       &
     &               ABS(v(i,j,1,nrhs))*cff2/dt(ng))
# endif
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=Istr,Iend
          angleC=Ucur(i,j)/Umag(i,j)
          angleW=COS(1.5_r8*pi-Dwave(i,j)-angler(i,j))
          bustrc(i,j)=tauCW(i,j)*angleC
          bustrw(i,j)=tauW(i,j)*angleW
          bustrcwmax(i,j)=tauCWmax(i,j)*angleW
          Ubot(i,j)=Ub(i,j)*angleW
          Ur(i,j)=Ucur(i,j)
!
          angleC=Vcur(i,j)/Umag(i,j)
          angleW=SIN(1.5_r8*pi-Dwave(i,j)-angler(i,j))
          bvstrc(i,j)=tauCW(i,j)*angleC
          bvstrw(i,j)=tauW(i,j)*angleW
          bvstrcwmax(i,j)=tauCWmax(i,j)*angleW
          Vbot(i,j)=Ub(i,j)*angleW
          Vr(i,j)=Vcur(i,j)
        END DO
      END DO
!
! Apply periodic or gradient boundary conditions for output purposes.
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
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bottom(:,:,izdef),                            &
     &                    bottom(:,:,izapp))
#endif

      RETURN
      END SUBROUTINE bblm_tile
