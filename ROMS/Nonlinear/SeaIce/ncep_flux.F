#include "cppdefs.h"
      MODULE ncep_flux_mod
#ifdef NCEP_FLUXES

!
!=====================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group
!================================================ Hernan G. Arango ===
!                                                                    !
!  This routine computes the bulk parameterization of surface wind   !
!  stress and surface net heat fluxes.                               !
!                                                                    !
!  References:                                                       !
!                                                                    !
!    Fairall, C.W., E.F. Bradley, D.P. Rogers, J.B. Edson and G.S.   !
!      Young, 1996:  Bulk parameterization of air-sea fluxes for     !
!      tropical ocean-global atmosphere Coupled-Ocean Atmosphere     !
!      Response Experiment, JGR, 101, 3747-3764.                     !
!                                                                    !
!    Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B.      !
!      Edson, and G.S. Young, 1996:  Cool-skin and warm-layer        !
!      effects on sea surface temperature, JGR, 101, 1295-1308.      !
!                                                                    !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk        !
!        parameterization of the air-sea exchange of heat and        !
!        water vapor including the molecular constraints at          !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                !
!                                                                    !
!  Adapted from Bergen Climate Model code written by Mats Bentsen    !
!                                                                    !
!                                                                    !
!=====================================================================
!
      implicit none

      PRIVATE
      PUBLIC ncep_flux

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ncep_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 17)
# endif
      CALL ncep_flux_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
# ifdef ICE_MODEL
     &                     liold(ng),                                   &
# endif
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
# ifdef WET_DRY
     &                     GRID(ng) % rmask_wet,                        &
# endif
# ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
# endif
     &                     FORCES(ng) % rhoa_n,                         &
     &                     FORCES(ng) % wg2_d, FORCES(ng) % wg2_m,      &
     &                     FORCES(ng) % cd_d,  FORCES(ng) % cd_m,       &
     &                     FORCES(ng) % ch_d,  FORCES(ng) % ch_m,       &
     &                     FORCES(ng) % ce_d,  FORCES(ng) % ce_m,       &
     &                     FORCES(ng) % nustr, FORCES(ng) % nvstr,      &
     &                     FORCES(ng) % srflx, FORCES(ng) % lrflx,      &
#ifdef RUNOFF
     &                     FORCES(ng) % runoff,                         &
#endif
     &                     FORCES(ng) % rain,                           &
     &                     FORCES(ng) % cloud,  FORCES(ng) % cawdir,    &
     &                     FORCES(ng) % icec,   FORCES(ng) % lhflx,     &
     &                     FORCES(ng) % shflx,  FORCES(ng) % Pair,      &
     &                     FORCES(ng) % skt,                            &
     &                     FORCES(ng) % tau_awx_n,                      &
     &                     FORCES(ng) % tau_awy_n,                      &
# ifdef ICE_MODEL
     &                     ICE(ng) % ai,    ICE(ng) % hi,               &
# ifdef ICE_THERMO
     &                     FORCES(ng) % qao_n,                          &
     &                     FORCES(ng) % qai_n,                          &
     &                     FORCES(ng) % p_e_n,                          &
     &                     ICE(ng) % hsn,   ICE(ng) % tis,              &
#  ifdef MELT_PONDS
     &                     ICE(ng) % apond,                             &
#  endif
     &                     ICE(ng) % coef_ice_heat,                     &
     &                     ICE(ng) % rhs_ice_heat,                      &
     &                     FORCES(ng) % snow_n,                         &
# endif
     &                     FORCES(ng) % tau_aix_n,                      &
     &                     FORCES(ng) % tau_aiy_n,                      &
     &                     FORCES(ng) % sustr_aw,                       &
     &                     FORCES(ng) % svstr_aw,                       &
# endif
     &                     FORCES(ng) % stflx,                          &
# if !defined ICE_MODEL
     &                     FORCES(ng) % sustr,                          &
     &                     FORCES(ng) % svstr,                          &
# endif
     &                     OCEAN(ng) % t,                               &
     &                     OCEAN(ng) % rho                              &
     &                     )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17)
# endif
      RETURN
      END SUBROUTINE ncep_flux
!
!********************************************************************
      SUBROUTINE ncep_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
# ifdef ICE_MODEL
     &                           liold,                                 &
# endif
# ifdef MASKING
     &                           rmask,                                 &
# endif
# ifdef WET_DRY
     &                           rmask_wet,                             &
# endif
# ifdef ICESHELF
     &                           zice,                                  &
# endif
     &                           rhoa_n, wg2_d, wg2_m, cd_d, cd_m,      &
     &                           ch_d, ch_m, ce_d, ce_m,                &
     &                           nustr, nvstr,                          &
     &                           srflx, lrflx,                          &
#ifdef RUNOFF
     &                           runoff,                                &
#endif
     &                           rain,                                  &
     &                           cloud, cawdir, icec, lhflx, shflx,     &
     &                           Pair, skt,                             &
     &                           tau_awx_n, tau_awy_n,                  &
# ifdef ICE_MODEL
     &                           ai, hi,                                &
# ifdef ICE_THERMO
     &                           qao_n,                                 &
     &                           qai_n,                                 &
     &                           p_e_n,                                 &
     &                           hsn, tis,                              &
#  ifdef MELT_PONDS
     &                           apond,                                 &
#  endif
     &                           coef_ice_heat, rhs_ice_heat, snow_n,   &
# endif
     &                           tau_aix_n, tau_aiy_n,                  &
     &                           sustr_aw, svstr_aw,                    &
# endif
     &                           stflx,                                 &
# if !defined ICE_MODEL
     &                           sustr,                                 &
     &                           svstr,                                 &
# endif
     &                           t,                                     &
     &                           rho)
!********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
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
# ifdef ICE_MODEL
      integer, intent(in) :: liold
# endif

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: rhoa_n(LBi:,LBj:)
      real(r8), intent(inout) :: wg2_d(LBi:,LBj:)
      real(r8), intent(inout) :: cd_d(LBi:,LBj:)
      real(r8), intent(inout) :: ch_d(LBi:,LBj:)
      real(r8), intent(inout) :: ce_d(LBi:,LBj:)
      real(r8), intent(inout) :: wg2_m(LBi:,LBj:)
      real(r8), intent(inout) :: cd_m(LBi:,LBj:)
      real(r8), intent(inout) :: ch_m(LBi:,LBj:)
      real(r8), intent(inout) :: ce_m(LBi:,LBj:)
      real(r8), intent(in) :: nustr(LBi:,LBj:)
      real(r8), intent(in) :: nvstr(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: lrflx(LBi:,LBj:)
#ifdef RUNOFF
      real(r8), intent(in) :: runoff(LBi:,LBj:)
#endif
      real(r8), intent(in) :: rain(LBi:,LBj:)
      real(r8), intent(in) :: cloud(LBi:,LBj:)
      real(r8), intent(in) :: cawdir(LBi:,LBj:)
      real(r8), intent(in) :: icec(LBi:,LBj:)
      real(r8), intent(in) :: lhflx(LBi:,LBj:)
      real(r8), intent(in) :: shflx(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: skt(LBi:,LBj:)
      real(r8), intent(out) :: tau_awx_n(LBi:,LBj:)
      real(r8), intent(out) :: tau_awy_n(LBi:,LBj:)
#  ifdef ICE_MODEL
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
#   ifdef ICE_THERMO
      real(r8), intent(out) :: qao_n(LBi:,LBj:)
      real(r8), intent(out) :: qai_n(LBi:,LBj:)
      real(r8), intent(out) :: p_e_n(LBi:,LBj:)
      real(r8), intent(in) :: hsn(LBi:,LBj:,:)
      real(r8), intent(in) :: tis(LBi:,LBj:)
#    ifdef MELT_PONDS
      real(r8), intent(in) :: apond(LBi:,LBj:,:)
#    endif
      real(r8), intent(out) :: coef_ice_heat(LBi:,LBj:)
      real(r8), intent(out) :: rhs_ice_heat(LBi:,LBj:)
      real(r8), intent(out) :: snow_n(LBi:,LBj:)
#    endif
      real(r8), intent(out) :: tau_aix_n(LBi:,LBj:)
      real(r8), intent(out) :: tau_aiy_n(LBi:,LBj:)
      real(r8), intent(out) :: sustr_aw(LBi:,LBj:)
      real(r8), intent(out) :: svstr_aw(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: stflx(LBi:,LBj:,:)
#  if !defined ICE_MODEL
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: rhoa_n(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wg2_d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: cd_d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ch_d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ce_d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wg2_m(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: cd_m(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ch_m(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ce_m(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: nustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: nvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: lrflx(LBi:UBi,LBj:UBj)
#ifdef RUNOFF
      real(r8), intent(in) :: runoff(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: cawdir(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: icec(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: skt(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tau_awx_n(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tau_awy_n(LBi:UBi,LBj:UBj)
#  ifdef ICE_MODEL
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
#   ifdef ICE_THERMO
      real(r8), intent(out) :: qao_n(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qai_n(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: p_e_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: tis(LBi:UBi,LBj:UBj)
#    ifdef MELT_PONDS
      real(r8), intent(in) :: apond(LBi:UBi,LBj:UBj,2)
#    endif
      real(r8), intent(out) :: coef_ice_heat(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhs_ice_heat(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: snow_n(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(out) :: tau_aix_n(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tau_aiy_n(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: sustr_aw(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr_aw(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#  if !defined ICE_MODEL
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: rhoswi, zu, zt, zq, slp, TseaK, tsrf, sfc_temp
      real(r8) :: tsrf_d, tsrf_m, cpair, fice, omice, rice
      real(r8) :: lhtfl, shtfl, taux_d, tauy_d, tau_d
      real(r8) :: hice_n, ua, sa, ta, qa, qsrf_d, qsh_i
      real(r8) :: le_i, dq_i, qlh_i, qlww, qlwi, qsww, qswi
      real(r8) :: albo, albi, albsn, albs, sigmaw, sigmai
      real(r8) :: emisw, emisi, rhocp, rhocpr, cawdif, nlwrs
      real(r8) :: albw_d, dswrf, cc, evap, qlh, dq, fqlat1
      real(r8) :: qsh, taux_m, tauy_m, f1, ustara, ustari
      real(r8) :: taufac, qsrf_m, tml, le, cd_i, cd_w
      real(r8) :: facw, rain_n, qsatw, qsati, rhoair, rhoSea
      real(r8) :: fac_thin, tai_lim, taw_lim
      real(r8) :: hsn_n
      real(r8) :: thki_n
      real(r8) :: thksn_n
      real(r8) :: apond_n
      real(r8) :: rhoO
      real(r8), parameter :: alb_i_dry=0.65_r8
      real(r8), parameter :: alb_i_wet=0.60_r8
      real(r8), parameter :: alb_s_dry=0.85_r8
      real(r8), parameter :: alb_s_wet=0.72_r8
# if !defined ICE_MODEL
      real(r8), parameter :: cdai=3.0_r8
# endif
      real(r8) :: cff1, cff2
!
# include "set_bounds.h"
!
!=====================================================================
!  Atmosphere-Ocean bulk fluxes parameterization.
!=====================================================================
!

      rhocp = 4.09374E+6_r8
      rhocpr = 0.2442754E-6_r8

      zu=10._r8 ! measurement height of wind [m]
      zt=10._r8 ! measurement height of temperature [m]
      zq=10._r8 ! measurement height of specific humidity [m]
      cpair=1004.7_r8 ! specific heat of dry air (J/K/kg)

      DO j=Jstr-1,JendT
        DO i=Istr-1,IendT
# ifdef ICE_THERMO
             IF(ai(i,j,liold).gt.min_a(ng)) THEN
                sfc_temp = ai(i,j,liold)*tis(i,j) +                     &
     &                    (1._r8-ai(i,j,liold))*t(i,j,N(ng),nrhs,itemp)
             ELSE
                sfc_temp = t(i,j,N(ng),nrhs,itemp)
             ENDIF
#undef DIAG_WPB
#ifdef DIAG_WPB
      if (i.eq.155.and.j.eq.482) then
        print *, 'ice4', sfc_temp, ai(i,j,liold), tis(i,j), &
              t(i,j,N(ng),nrhs,itemp)
      endif
#endif
# else
          sfc_temp = skt(i,j) - 273.16_r8
# endif
          TseaK=t(i,j,N(ng),nrhs,itemp)+273.16
          rhoSea = rho(i,j,N(ng))+1000.0_r8
          rhoswi=1._r8/rhoSea
! Convert slp from mb to Pascals
          slp = Pair(i,j)*100.0_r8
          tsrf_d = skt(i,j)
          tsrf = sfc_temp + 273.16
          tml = TseaK

#ifdef ICE_MODEL
          fice = max(ai(i,j,liold), 0.0_r8)
#else
          fice = icec(i,j)
#endif
          omice = 1.0_r8-fice
          rice = icec(i,j)
#if defined ICE_MODEL && defined ICE_THERMO
          hice_n = hi(i,j,liold)
          hsn_n = hsn(i,j,liold)
# ifdef MELT_PONDS
          apond_n = apond(i,j,liold)
# else
          apond_n = 0._r8
# endif
#else
          hice_n = 3._r8*fice
          hsn_n = 0.2_r8*fice
          apond_n = 0._r8
#endif
          thki_n = hice_n/(MAX(fice,0.01_r8))
          thksn_n = hsn_n/(MAX(fice,0.01_r8))
!
          lhtfl = lhflx(i,j)
          shtfl = shflx(i,j)
          taux_d = nustr(i,j)
          tauy_d = nvstr(i,j)
          tau_d = sqrt(taux_d*taux_d + tauy_d*tauy_d)
!

!
! --- ------------------------------------------------------------------
! --- compute the atmospheric state by using the prescribed momentum
! --- heat fluxes and the prescribed sea state
! --- ------------------------------------------------------------------
!
! --- first guess on the atmospheric state by using the transfer
! --- coefficients and density from the previous time step
!
      ua=sqrt(.5*(-wg2_d(i,j)                                           &
     &            +sqrt(wg2_d(i,j)*wg2_d(i,j)                           &
     &                 +4.*(tau_d/(rhoa_n(i,j)*cd_d(i,j)))**2)))
      sa=sqrt(ua*ua+wg2_d(i,j))
      ta=tsrf_d-.0098*zt-shtfl/(rhoa_n(i,j)*cpair*ch_d(i,j)*sa)
      ta = max(ta,203._r8)
      qsrf_d=rice*qsati(tsrf_d,slp)+(1.-rice)*qsatw(tsrf_d,slp)
      le=(2.501-0.00237*(tsrf_d-273.16))*1.e6
      qa=qsrf_d-lhtfl/(rhoa_n(i,j)*le*ce_d(i,j)*sa)
      rhoa_n(i,j)=rhoair(ta,qa,slp)
!
! --- update the transfer coefficients and gustiness
!
      CALL bulktf(ua,zu,ta,zt,qa,zq,tsrf_d,qsrf_d,rice,                 &
     &            cd_d(i,j),ch_d(i,j),ce_d(i,j),wg2_d(i,j),             &
     &            cd_i,cd_w)
!
! --- update the atmospheric state
      ua=sqrt(.5*(-wg2_d(i,j)                                           &
     &            +sqrt(wg2_d(i,j)*wg2_d(i,j)                           &
     &                 +4.*(tau_d/(rhoa_n(i,j)*cd_d(i,j)))**2)))
      sa=sqrt(ua*ua+wg2_d(i,j))
      ta=tsrf_d-.0098*zt-shtfl/(rhoa_n(i,j)*cpair*ch_d(i,j)*sa)
      qa=qsrf_d-lhtfl/(rhoa_n(i,j)*le*ce_d(i,j)*sa)
      rhoa_n(i,j)=rhoair(ta,qa,slp)
!
! --- ------------------------------------------------------------------
! --- update transfer coefficients and gustiness with the computed
! --- atmospheric state the models ocean state
! --- ------------------------------------------------------------------
!
      tsrf_m=fice*tsrf+omice*tml
      qsrf_m=fice*qsati(tsrf,slp)+omice*qsatw(tml,slp)
!
      CALL bulktf(ua,zu,ta,zt,qa,zq,tsrf_m,qsrf_m,fice,                 &
     &            cd_m(i,j), ch_m(i,j), ce_m(i,j), wg2_m(i,j),          &
     &            cd_i, cd_w)
!
! Endelig beregning av vindstresskorreksjon og frisksjonshastighet til
! modell (hice_n er istykkelse i m):
!
! --- ------------------------------------------------------------------
! --- compute correction of the wind stress on the surface and friction
! --- velocity
! --- ------------------------------------------------------------------
!
      sa=sqrt(ua*ua+wg2_m(i,j))
      taufac = rhoa_n(i,j)*cd_m(i,j)*sa*ua/tau_d

!
! --- friction velocity (m/s) is modified in the presence
! --- of ice
      ustara=sqrt(cd_m(i,j)*sa*ua*rhoa_n(i,j)*rhoswi)
!      ustari=sqrt(tauice(i,j)*rhoswi)
# if !defined ICE_MODEL
      ustari=sqrt(cdai*sa*ua*rhoa_n(i,j)*rhoswi)
# else
      ustari=sqrt(cdai(ng)*sa*ua*rhoa_n(i,j)*rhoswi)
# endif
      f1=fice*min(1._r8,hice_n)

! Stresskorreksjonen benyttes senere slik at:
      taux_m=taufac*taux_d
      tauy_m=taufac*tauy_d
# ifdef ICE_MODEL
      facw = fice*(cd_i/cd_w) + omice
      fac_thin = 0.5_r8*(1.0_r8-COS(2.0_r8*pi*MIN(                      &
     &  hi(i,j,liold)/(ai(i,j,liold)+0.02_r8),0.5_r8)))
#   ifdef ICE_MOM_BULK
      cd_i = cdai(ng)
      cd_i = 0.5_r8*(1.0_r8-COS(2.0_r8*pi*MIN(                          &
     &  hi(i,j,liold),0.5_r8)))*cdai(ng)
#   else
      cd_i = cd_w + fac_thin*(cd_i-cd_w)
#   endif
      tau_awx_n(i,j) = taux_m/facw
      tau_awy_n(i,j) = tauy_m/facw
      tau_aix_n(i,j) = (cd_i/cd_w)*tau_awx_n(i,j)
      tau_aiy_n(i,j) = (cd_i/cd_w)*tau_awy_n(i,j)
# else
      tau_awx_n(i,j) = taux_m
      tau_awy_n(i,j) = tauy_m
# endif
!
! -- sensible heat flux for open water (W/m**2)
!
      taw_lim = ta
      qsh =rhoa_n(i,j)*cpair*ch_m(i,j)*sa*(taw_lim+0.0098*zt-tml)
!
! --- latent heat for open water (W/m**2)
      fqlat1=rhoa_n(i,j)*ce_m(i,j)*sa
      le=(2.501-0.00237*(tml-273.16))*1.e6
      dq=qa-qsatw(tml,slp)
      qlh =fqlat1*le*dq
!
! --- evaporation (m/s)
      evap = -fqlat1*dq*1.e-3

! Short-wave radiation to ocean
      cc = cloud(i,j)
      cc = max(cc, 0.0_r8)
      dswrf = srflx(i,j)
      albw_d=0.065_r8
      cawdif = 1._r8-albw_d
      qsww =dswrf*(cawdir(i,j)*(1._r8-cc)+cawdif*cc)

! Net long-wave radiation to ocean
      nlwrs = lrflx(i,j)
      sigmaw =  5.67E-8_r8
      emisw = 0.97_r8
      qlww = nlwrs + 4._r8*sigmaw*emisw*ta**3*(tml-tsrf_d)

#ifdef ICE_THERMO
! Net upward oceanic heat flux without ice, in W/m**2
      qao_n(i,j) = (-qsh - qlh - qsww + qlww)
! Net precipitation - evaporation in m/s
      IF(ta.ge.273.16) THEN
         rain_n = rain(i,j)/1000._r8
      ELSE
         rain_n = omice*rain(i,j)/1000._r8
         snow_n(i,j) = rain(i,j)/rhosnow_dry(ng)
      ENDIF
      p_e_n(i,j) = (evap*omice-rain_n)
      stflx(i,j,isalt) = p_e_n(i,j)
#ifdef RUNOFF
      stflx(i,j,isalt) = stflx(i,j,isalt) - runoff(i,j)/1000._r8
#endif
#undef DIAG_WPB
#ifdef DIAG_WPB
      IF(i.eq.354.and.j.eq.101) THEN
        write(9,100) tdays,qao_n(i,j),-qsh,-qlh,-qsww,qlww
        write(8,100) tdays,ta,skt(i,j),tsrf_d,tml
        write(10,100) tdays,evap,rain_n,runoff(i,j)*0.001_r8,stflx(i,j,isalt)
 100    format(f16.5,5e16.4)
!        write(*,*) '  From ncep_flux: '
!        write(*,*) 'qao_n = ',qao_n(i,j)
!        write(*,*) '-qsh = ',-qsh
!        write(*,*) '-qlh = ',-qlh
!        write(*,*) '-qsww = ',-qsww
!        write(*,*) 'qlww = ',qlww
!        write(*,*) ' '
      ENDIF
#endif


#else
! Net downward oceanic heat flux in kinematic units
      stflx(i,j,itemp) = (qsh + qlh + qsww - qlww) *rhocpr
! Evaporation-precipitation salt flux in m/s (kinematic units)
      stflx(i,j,isalt) = evap-rain(i,j)*0.001_r8
# ifdef RUNOFF
      stflx(i,j,isalt) = stflx(i,j,isalt) - runoff(i,j)*0.001_r8
# endif
#endif

!  Correct ocean net short-wave radiation for ice cover and convert to
!  kinematic units
#ifdef ICE_THERMO
      srflx(i,j) = omice*qsww*rhocpr
#else
      srflx(i,j) = qsww*rhocpr
#endif

#ifdef ICE_THERMO
      coef_ice_heat(i,j) = 0._r8
      rhs_ice_heat(i,j) = 0._r8

! Sensible heat flux over ice
!
      tai_lim = ta
      qsh_i = rhoa_n(i,j)*cpair*ch_m(i,j)*sa*(tai_lim+0.0098*zt-tsrf)

      coef_ice_heat(i,j) = coef_ice_heat(i,j) +                         &
     &      rhoa_n(i,j)*cpair*ch_m(i,j)*sa
      rhs_ice_heat(i,j) = rhs_ice_heat(i,j) +                           &
     &      rhoa_n(i,j)*cpair*ch_m(i,j)*sa*(tai_lim+0.0098*zt)
!
! Latent heat flux over ice
      le_i=2.834E+6_r8
      dq_i=qa-qsati(tsrf,slp)
      qlh_i =fqlat1*le_i*dq_i
!
! Correct NCEP latent heat flux over ice to Jordan et al (1999):
!      qlh_i = 0.52_r8*qlh_i - 1.01_r8 ! Reduce cooling in summer
      rhs_ice_heat(i,j) = rhs_ice_heat(i,j) +                           &
     &     qlh_i
!
! Calculate the ice/snow albedo
# ifdef ICE_ALB_EC92
! Ice and snow albedo is calculated from Ebert and Curry,1992b
          albi=0._r8
          albsn=0._r8
          albi=0.082409_r8*LOG(thki_n)+0.485472_r8
          IF (thki_n.GE.1._r8) albi=0.07616_r8*thki_n+0.414492_r8
          IF (thki_n.GE.2._r8) albi=0.561632_r8
!approximated values for albsn(depends on COSZ, but small variation)
          albsn=0.83_r8
          IF (apond_n.GT.0._r8) THEN
            IF (thksn_n.GE.0.1_r8) THEN
              albsn=0.701009_r8
            ELSE
              albsn=albi+(0.701009_r8-albi)*thksn_n/0.1
            ENDIF
          ENDIF
          albs=albi
          IF (hsn_n.GT.0._r8) albs=albsn
!          IF (apond_n.GT.0._r8) albs=0.100737_r8                        &
!     &         +0.518_r8*EXP(-8.1_r8 *apond_n-0.47_r8)                  &
!     &         +0.341_r8*EXP(-31.8_r8*apond_n-0.94_r8)                  &
!     &         +0.131_r8*EXP(-2.6_r8 *apond_n-3.82_r8)
!
# else
      cff1 = alb_s_wet - alb_s_dry
      cff2 = alb_i_wet - alb_i_dry
          IF(hsn(i,j,liold).gt.0._r8) THEN
            IF (tsrf-273.16_r8.gt. -1._r8) THEN
                albs = cff1*(tsrf-272.16_r8)+alb_s_dry
            ELSE
                albs = alb_s_dry
            ENDIF
          ELSE
            IF (tsrf-273.16_r8.gt. -1._r8) THEN
                albs = cff2*(tsrf-272.16_r8)+alb_i_dry
            ELSE
                albs = alb_i_dry
            ENDIF
          ENDIF
#undef DIAG_WPB
#ifdef DIAG_WPB
      if(i.eq.155.and.j.eq.482) then
        print *, 'ice3', albs, tsrf, fice, tsrf_d, tsrf_m, sfc_temp, &
                  tis(i,j), skt(i,j)
      endif
#endif

# endif
!
! Short-wave radiation to ice
      qswi = (1._r8-albs)*dswrf

      rhs_ice_heat(i,j) = rhs_ice_heat(i,j) +                           &
     &     qswi

! Net long-wave radiation from ice
      sigmai = 5.67E-8_r8 ! 5.512_e8-8
      emisi = 0.97_r8
      qlwi = nlwrs + 4._r8*sigmai*emisi*ta**3*(tsrf-tsrf_d)

      coef_ice_heat(i,j) = coef_ice_heat(i,j) +                         &
     &     4._r8*sigmai*emisi*ta**3
      rhs_ice_heat(i,j) = rhs_ice_heat(i,j) -                           &
     &     nlwrs + 4._r8*sigmai*emisi*ta**3*tsrf_d

! Correct for Kelvin to Celsius in RHS
      rhs_ice_heat(i,j) = rhs_ice_heat(i,j) -                           &
     &     coef_ice_heat(i,j)*273.16_r8

! Net heat flux from ice to atmosphere
       qai_n(i,j) =-qsh_i - qlh_i - qswi + qlwi

#undef DIAG_WPB
#ifdef DIAG_WPB
      if(i.eq.155.and.j.eq.482) then
        write(*,*) 'ice',tdays,qai_n(i,j),-qsh_i,-qlh_i,-qswi,qlwi,     &
     &              tsrf,tsrf_d
        print *, 'ice2', dswrf, albs, tsrf, alb_s_wet, alb_i_wet
      END IF
#endif

#endif
      END DO
      END DO

# ifdef ICE_MODEL
      DO j=JstrT,JendT
        DO i=Istr,IendT
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i-1,j,N(ng)))
          sustr_aw(i,j) = 0.5_r8*(tau_awx_n(i-1,j)+tau_awx_n(i,j))/rhoO
        END DO
      END DO

      DO j=Jstr,JendT
        DO i=IstrT,IendT
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i,j-1,N(ng)))
          svstr_aw(i,j) = 0.5_r8*(tau_awy_n(i,j-1)+tau_awy_n(i,j))/rhoO
        END DO
      END DO
# else
      DO j=JstrT,JendT
        DO i=Istr,IendT
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i-1,j,N(ng)))
          sustr(i,j) = 0.5_r8*(tau_awx_n(i-1,j)+tau_awx_n(i,j))/rhoO
        END DO
      END DO

      DO j=Jstr,JendT
        DO i=IstrT,IendT
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i,j-1,N(ng)))
          svstr(i,j) = 0.5_r8*(tau_awy_n(i,j-1)+tau_awy_n(i,j))/rhoO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoa_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wg2_d)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          cd_d)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ch_d)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ce_d)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wg2_m)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          cd_m)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ch_m)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ce_m)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tau_awx_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tau_awy_n)
# ifdef ICE_MODEL
#  ifdef ICE_THERMO
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          qao_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          qai_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          p_e_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          coef_ice_heat)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhs_ice_heat)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          snow_n)
#  endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tau_aix_n)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tau_aiy_n)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr_aw)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr_aw)
# endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,isalt))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,itemp))
# ifndef ICE_MODEL
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
# endif
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    rhoa_n, wg2_d, cd_d, ch_d)
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ce_d, wg2_m, cd_m, ch_m)
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ce_m, srflx, tau_awx_n, tau_awy_n)
#  ifdef ICE_MODEL
#   ifdef ICE_THERMO
       CALL mp_exchange2d (ng, tile, iNLM, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    qao_n, qai_n, p_e_n, coef_ice_heat)
       CALL mp_exchange2d (ng, tile, iNLM, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    rhs_ice_heat, snow_n)
#   endif
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    tau_aix_n, tau_aiy_n, sustr_aw, svstr_aw)
#  endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    stflx(:,:,itemp),                             &
     &                    stflx(:,:,isalt))
#  ifndef ICE_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    sustr, svstr)
#  endif
# endif
      RETURN
      END SUBROUTINE ncep_flux_tile
#endif
      END MODULE ncep_flux_mod

! -----------------------------------------------------------------
      SUBROUTINE bulktf(du,zu,ta,zt,qa,zq,ts,qs,cice,cd,ch,ce,wg2,      &
     &                  cd_i,cd_w)
!
! --- this routine computes momentum, sensible heat, and latent heat
! --- transfer coefficients by one interation of a bulk flux algorithm
! --- (Fairall et al.,1996)
!
      USE mod_kinds

      implicit none

! --- input variables:
! ---   du     - velocity difference (wind speed - current) [m/s]
! ---   zu     - measurement height of wind speed [m]
! ---   ta     - air temperature [K]
! ---   zt     - measurement height of air temperature [m]
! ---   qa     - specific humidity of air [kg/kg]
! ---   zq     - measurement height of specific humidity [m]
! ---   ts     - sea surface temperature [K]
! ---   qs     - specific humidity of air at the surface [kg/kg]
! ---   cice   - ice concentration (0 open water, 1 sea ice cover)
! ---   cd     - momentum transfer coeffisient
! ---   ch     - sensible heat transfer coeffisient
! ---   ce     - latent heat transfer coeffisient
! ---   wg2    - gustiness squared [m^2/s^2]
!
      real(r8), intent(in) :: du
      real(r8), intent(in) :: zu
      real(r8), intent(in) :: ta
      real(r8), intent(in) :: zt
      real(r8), intent(in) :: qa
      real(r8), intent(in) :: zq
      real(r8), intent(in) :: ts
      real(r8), intent(in) :: qs
      real(r8), intent(in) :: cice
      real(r8), intent(inout) :: cd
      real(r8), intent(inout) :: ch
      real(r8), intent(inout) :: ce
      real(r8), intent(inout) :: wg2
      real(r8), intent(out) :: cd_i
      real(r8), intent(out) :: cd_w

!
! --- parameters:
! ---   eps    - molecular weight ratio of dry air and water vapour
! ---   t0     - freezing point of water [K]
! ---   zi     - inversion height [m]
! ---   g      - acceleration due to gravity [m/s^2]
! ---   beta   - empirical constant relating gustiness to convective
! ---            scaling velocity
! ---   alpha  - Charnock constant
! ---   k      - von Karman constant
!
      real(r8), parameter :: eps=.62197_r8, cv=1._r8/eps-1._r8
      real(r8), parameter :: t0=273.16_r8
      real(r8), parameter :: zi=600._r8, g=9.8_r8, beta=1.2_r8
      real(r8), parameter :: athird=1._r8/3._r8
      real(r8), parameter :: alpha=.011_r8, gi=1._r8/g
      real(r8), parameter :: k=.4_r8, ki=2.5_r8
!
      real(r8) :: tv, tac, visca, dt, dq, du1, du2, s, ustar2, ustar
      real(r8) :: fac, tstar, qstar
      real(r8) :: tvstar, li, w3, wg, zetau, zetat, zetaq, z0, cd2
      real(r8) :: psiu, reu, ret, req
      real(r8) :: z0t, z0q, ct2, psitq, cq2, zice, zwat
      real(r8) :: cd2i, cd2w
!
! --- virtual temperature
      tv=ta*(1.+cv*qa)
!
! --- kinematic viscosity of dry air - Andreas (1989) CRREL Rep. 89-11
! --- [m^2/s]
      tac=ta-t0
      visca=1.326E-5_r8*(1._r8+tac*(6.542E-3_r8+tac*                    &
     &                            (8.301E-6_r8-tac*4.84E-9_r8)))
!
! --- temperature difference
      dt=ta-ts+.0098_r8*zt
!
! --- humidity difference
      dq=qa-qs
!
! --- initial values
      du1=max(du,1.E-2_r8)
      du2=du1*du1
      s=SQRT(du2+wg2)
      ustar2=cd*s*du1
      ustar=SQRT(ustar2)
      fac=ustar/(cd*du1)
      tstar=fac*ch*dt
      qstar=fac*ce*dq
!
! --- inverse Monin-Obukov lenght
      tvstar=tstar*(1+cv*qa)+cv*ta*qstar
      li=min(3._r8/zu,g*k*tvstar/(ustar2*tv))
!
! --- gustiness
      w3=-zi*g*ustar*tvstar/ta
      wg=max(.1_r8,beta*max(0._r8,w3)**athird)
!
! --- wind speed included gustiness
      s=sqrt(du2+wg*wg)
!
! --- nondimensional heights
      zetau=zu*li
      zetat=zt*li
      zetaq=zq*li
!
! --- roughness length - choose roughness lenght for sea ice
! --- corresponding to smooth multiyear ice, Guest and Davidson (1991)
! --- JGR 96, 4709-4721
      z0=cice*2.E-3_r8+(1._r8-cice)*                                    &
     &   (0.11_r8*visca/ustar+alpha*ustar2*gi)
!
! --- square root of the momentum transfer coefficient
      cd2=k/max(7._r8,LOG(zu/z0)-psiu(zetau))
!
! --- update friction velocity
      ustar=cd2*SQRT(s*du1)
!
! --- roughness Reynolds numbers
      reu=ustar*z0/visca
      CALL lkb(reu,ret,req)
!
! --- temperature and humidity roughness scales
      fac=visca/ustar
      z0t=fac*ret
      z0q=fac*req
!
! --- transfer coeffisient components for temperature and humidity
      ct2=k/max(7._r8,LOG(zt/z0t)-psitq(zetat))
      cq2=k/max(7._r8,LOG(zq/z0q)-psitq(zetaq))
!
! --- momentum transfer coeffisient
      cd=cd2*cd2
!
! --- heat transfer coeffisients
      ch=cd2*ct2
      ce=cd2*cq2
!
! --- gustiness squared
      wg2=wg*wg

! Get drag coefficients for open water, sea ice cover
      zwat=(0.11_r8*visca/ustar+alpha*ustar2*gi)
      zice=2.E-3_r8
      cd2i=k/MAX(7._r8,LOG(zu/zice)-psiu(zetau))
      cd2w=k/MAX(7._r8,LOG(zu/zwat)-psiu(zetau))
      cd_i = cd2i*cd2i
      cd_w = cd2w*cd2w

      RETURN
      END SUBROUTINE bulktf
!
! --- ------------------------------------------------------------------
!
      SUBROUTINE lkb(reu, ret, req)
!
! --- computes the roughness Reynolds for temperature and humidity
! --- following Liu, Katsaros and Businger (1979), J. Atmos.  Sci., 36,
! --- 1722-1735.
!
      USE mod_kinds

      real(r8), intent(in) :: reu
      real(r8), intent(out) :: ret
      real(r8), intent(out) :: req
!
      integer :: i
!
      real(r8), dimension(8,2) ::  a
      real(r8), dimension(8) :: re
      real(r8), dimension(8,2) ::  b

      data a /                                                          &
     &       0.177_r8, 1.376_r8, 1.026_r8, 1.625_r8,                    &
     &       4.661_r8, 34.904_r8, 1667.19_r8, 5.88E5_r8,                &
     &       0.292_r8, 1.808_r8, 1.393_r8, 1.956_r8,                    &
     &       4.994_r8, 30.709_r8, 1448.68_r8, 2.98E5_r8 /

      data b                                                            &
     &      /0._r8, 0.929_r8, -0.599_r8, -1.018_r8,                     &
     &      -1.475_r8, -2.067_r8, -2.907_r8, -3.935_r8,                 &
     &       0._r8, 0.826_r8, -0.528_r8, -0.870_r8,                     &
     &      -1.297_r8, -1.845_r8, -2.682_r8, -3.616_r8/

      data re                                                           &
     &      /0.11_r8, 0.825_r8, 3.0_r8, 10.0_r8,                        &
     &       30.0_r8, 100._r8, 300._r8, 1000._r8/
!
      i=1
 10   IF (reu.gt.re(i).and.i.lt.8) THEN
        i=i+1
        GOTO 10
      ENDIF
!
      ret=a(i,1)*reu**b(i,1)
      req=a(i,2)*reu**b(i,2)
!
      RETURN
      END SUBROUTINE lkb

!
! --- ------------------------------------------------------------------
!
      FUNCTION psiu(zeta)
!
! --- Monin-Obukhov similarity velocity profile function
!
      USE mod_kinds

      real(r8) :: psiu, zeta
!
      real(r8), parameter :: athird=1._r8/3._r8
      real(r8), parameter :: pi=3.141592653589793_r8
      real(r8), parameter :: sqrt3=1.732050807568877_r8
      real(r8), parameter :: sqrt3i=.5773502691896258_r8
!
      real(r8) :: x, psik, y, psic, f
!
      IF (zeta.eq.0._r8) THEN
        psiu=0.
      ELSEIF (zeta.gt.0._r8) THEN
        psiu=-4.7_r8*zeta
      ELSE
        x=(1._r8-16._r8*zeta)**.25_r8
        psik=2._r8*LOG((1._r8+x)*.5_r8)+LOG((1._r8+x*x)*.5_r8)          &
     &       -2._r8*ATAN(x)+pi*.5_r8
        y=(1.-12.87*zeta)**athird
        psic=1.5_r8*LOG((y*y+y+1._r8)*athird)                           &
     &      -sqrt3*ATAN((2._r8*y+1._r8)*sqrt3i)                         &
     &      +pi*sqrt3i
        f=1._r8/(1._r8+zeta*zeta)
        psiu=f*psik+(1._r8-f)*psic
      ENDIF
!
      RETURN
      END FUNCTION psiu
!
! --- ------------------------------------------------------------------
!
      FUNCTION psitq(zeta)
!
! --- Monin-Obukhov similarity temperature and humidity profile function
!
      USE mod_kinds

      real(r8) :: psitq, zeta
!
      real(r8), parameter :: athird=1._r8/3._r8
      real(r8), parameter :: pi=3.141592653589793_r8
      real(r8), parameter :: sqrt3=1.732050807568877_r8
      real(r8), parameter :: sqrt3i=.5773502691896258_r8
!
      real(r8) :: x, psik, y, psic, f
!
      IF (zeta.eq.0._r8) THEN
        psitq=0._r8
      ELSEIF (zeta.gt.0._r8) THEN
        psitq=-4.7_r8*zeta
      ELSE
        x=(1._r8-16._r8*zeta)**.25_r8
        psik=2._r8*LOG((1._r8+x*x)*.5_r8)
        y=(1._r8-12.87_r8*zeta)**athird
        psic=1.5_r8*LOG((y*y+y+1._r8)*athird)                           &
     &             -sqrt3*ATAN((2._r8*y+1._r8)*sqrt3i)                  &
     &      +pi*sqrt3i
        f=1./(1._r8+zeta*zeta)
        psitq=f*psik+(1._r8-f)*psic
      ENDIF
!
      RETURN
      END FUNCTION psitq

      FUNCTION qsatw(t,p)
!
! --- computes saturation specific humidity [kg/kg] over open water,
! --- Buck (1981) JAM 20, 1527-1532
!
! --- input variables:
! ---   t      - air temperature [K]
! ---   p      - atmospheric pressure [Pa]
!
      USE mod_kinds

      real(r8) :: qsatw, t, p
!
! --- parameters:
! ---   eps    - molecular weight ratio of dry air and water vapour
!
      real(r8), parameter :: eps=0.62197_r8
!
      real(r8) :: e
!
      e=611.21_r8*(1.0007_r8+3.46E-8_r8*p)*EXP(17.502_r8*(t-273.16_r8)/ &
     &            (t-32.19_r8))
      qsatw=eps*e/(p-(1._r8-eps)*e)
!
      RETURN
      END FUNCTION qsatw
!
! --- ------------------------------------------------------------------
!
      FUNCTION qsati(t,p)
!
! --- computes saturation specific humidity [kg/kg] over sea ice,
! --- Parkinson and Washington (1979) JGR 84, 311-337
!
! --- input variables:
! ---   t      - air temperature [K]
! ---   p      - atmospheric pressure [Pa]
!
      USE mod_kinds

      real(r8) :: qsati, t, p
!
! --- parameters:
! ---   eps    - molecular weight ratio of dry air and water vapour
!
      real(r8), parameter :: eps=0.62197_r8
!
      real(r8) :: e
!
      e=611._r8*10._r8**(9.5_r8*(t-273.16_r8)/(t-7.66_r8))
      qsati=eps*e/(p-(1._r8-eps)*e)
!
      RETURN
      END FUNCTION qsati
!
! --- ------------------------------------------------------------------
!
      FUNCTION rhoair(t,q,p)
!
! --- computes air density [kg/m^3]
!
! --- input variables:
! ---   t      - air temperature [K]
! ---   q      - specific humidity of air [kg/kg]
! ---   p      - atmospheric pressure [Pa]
!

      USE mod_kinds

      real(r8) :: rhoair, t, q, p
!
! --- parameters:
! ---   eps    - molecular weight ratio of dry air and water vapour
! ---   rgas   - gas constant for dry air [J/(kg*K)]
!
      real(r8), parameter :: eps=0.62197_r8
      real(r8), parameter :: cv=1.0_r8/eps-1.0_r8
      real(r8), parameter :: rgas=287.04_r8
!
      real(r8) :: tv
!
! --- virtual temperature
      tv=t*(1.0_r8+cv*q)
!
! --- moist air density [kg/m^3]
      rhoair=p/(rgas*tv)
!
      RETURN
      END FUNCTION rhoair


