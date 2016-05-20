     SUBROUTINE ice_thermo (ng, tile)
!
!*************************************************** W. Paul Budgell ***
!  Copyright (c) 2002-2015 ROMS/TOMS Group                             !
!************************************************** Hernan G. Arango ***
!                                                                      !
!  This subroutine evaluates the ice thermodynamic growth and decay    !
!  term based on  Mellor and Kantha (1989) and Parkinson and           !
!  Washington (1979)                                                   !
!                                                                      !
!***********************************************************************
!

      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_ice
      USE mod_forces
      USE mod_stepping
#ifdef ICE_SHOREFAST
      USE mod_coupling
#endif
#ifdef AICLM_NUDGING
      USE mod_clima
#endif

      implicit none

      integer, intent(in) :: ng, tile

#include "tile.h"

#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 51)
#endif

      CALL ice_thermo_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), liold(ng), linew(ng),             &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
#ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
#endif
#ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#endif
#ifdef ICE_SHOREFAST
     &                      GRID(ng) % h,                               &
     &                      COUPLING(ng) % Zt_avg1,                     &
#endif
#ifdef AICLM_NUDGING
     &                      CLIMA(ng) % aiclm,                          &
     &                      CLIMA(ng) % hiclm,                          &
     &                      CLIMA(ng) % AInudgcof,                      &
#endif
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      OCEAN(ng) % t,                              &
     &                      ICE(ng) % wfr,                              &
     &                      ICE(ng) % wai,                              &
     &                      ICE(ng) % wao,                              &
     &                      ICE(ng) % wio,                              &
     &                      ICE(ng) % wro,                              &
     &                      ICE(ng) % ai,                               &
     &                      ICE(ng) % hi,                               &
     &                      ICE(ng) % hsn,                              &
     &                      ICE(ng) % ageice,                           &
#ifdef MELT_PONDS
     &                      ICE(ng) % apond,                            &
     &                      ICE(ng) % hpond,                            &
#endif
     &                      ICE(ng) % tis,                              &
     &                      ICE(ng) % ti,                               &
     &                      ICE(ng) % enthalpi,                         &
     &                      ICE(ng) % hage,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % coef_ice_heat,                    &
     &                      ICE(ng) % rhs_ice_heat,                     &
     &                      ICE(ng) % s0mk,                             &
     &                      ICE(ng) % t0mk,                             &
     &                      ICE(ng) % io_mflux,                         &
#if defined ICE_BIO && defined BERING_10K
     &                      ICE(ng) % IcePhL,                           &
     &                      ICE(ng) % IceNO3,                           &
     &                      ICE(ng) % IceNH4,                           &
#endif
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
     &                      FORCES(ng) % qai_n,                         &
     &                      FORCES(ng) % qi_o_n,                        &
     &                      FORCES(ng) % qao_n,                         &
     &                      FORCES(ng) % snow_n,                        &
     &                      FORCES(ng) % rain,                          &
     &                      FORCES(ng) % stflx)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 51)
#endif
      RETURN
      END SUBROUTINE ice_thermo
!
!***********************************************************************
      SUBROUTINE ice_thermo_tile (ng, tile,                             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs, liold, linew,                       &
#ifdef MASKING
     &                        rmask,                                    &
#endif
#ifdef WET_DRY
     &                        rmask_wet,                                &
#endif
#ifdef ICESHELF
     &                        zice,                                     &
#endif
#ifdef ICE_SHOREFAST
     &                        h, Zt_avg1,                               &
#endif
#ifdef AICLM_NUDGING
     &                        aiclm, hiclm, AInudgcof,                  &
#endif
     &                        z_r, z_w, t,                              &
     &                        wfr, wai, wao, wio, wro,                  &
     &                        ai, hi, hsn, ageice,                      &
#ifdef MELT_PONDS
     &                        apond, hpond,                             &
#endif
     &                        tis, ti, enthalpi, hage,                  &
     &                        ui, vi, coef_ice_heat, rhs_ice_heat,      &
     &                        s0mk, t0mk, io_mflux,                     &
#if defined ICE_BIO && defined BERING_10K
     &                        IcePhL, IceNO3, IceNH4,                   &
#endif
     &                        sustr, svstr,                             &
     &                        qai_n, qi_o_n, qao_n,                     &
     &                        snow_n,                                   &
     &                        rain,                                     &
     &                        stflx)
!***********************************************************************
!
!  Original comment:
!      beregner varmefluxer og produskjonsrater
!      og oppdaterer tis (t3 i mellor et.al.)
!
!  means compute heat fluxes and ice production rates:
!
!      wai(i,j)=-(qai(i,j) -qi2(i,j)) /(hfus1(i,j)*rhosw)
!
!  and up date the internal ice temperature (t3 in Mellor et all).
!
!        tis(i,j)=tis(i,j)+del
!
!        the following global arrays are calculated:
!        (description is given below)
!
!        apond
!        hpond
!        ageice
!        qai
!        qio
!        qi2
!        wsm
!        wai
!        wro
!        tis
!        t2
!        hfus1
!        coa
!
!     D1 = BULK SENSIBLE HEAT TRANSFER COEFFICIENT          [J/(K*m**3)]
!     D2 = LATENT HEAT TRANSFER COEFFICIENT,                [J/(K*m**3)]
!          D2I FOR OVER ICE, D2W FOR OVER WATER
!     D3 = STEFAN-BOLTZMAN CONST. * SURFACE EMISSIVITY      [W/(K**4*m**2)]
!
!        parameters:
!
!         inp from atmosphere model:
!               wind_speed(im,jm)     -  abs(wind_10_meter)
!               Tair(im,jm)       -  atmos. temperature
!               rh(im,jm)       -  atmosphere specific humidity
!               snow(im,jm)      -  snow fall rate
!
!
!         inp from ocean model:
!
!             t0mk(im,jm)       -  sea surface temperature
!             rmask(im,jm)      -  pointer (land/ocean)
!
!               dtice             -  time step
!                     dtice=float(isplitc)*dti
!
!       global variables transfered by module:
!
!            ----  needs to be initiated elswhere ----
!
!   (rads)   sw_flux(i,j)      -  incoming short wave radiation
!   (rads)   lwrad(i,j)        -  incoming long wave radiation
!
!            qi2(i,j)        -  heat flux in ice
!            hi(i,j,linew)   -  ice mass (divided by area)
!            ai(i,j,linew)   -  ice concentration
!            hsn(i,j,linew)  -  mass snow (pr. area) ai*snow_thick
!            ti(i,j,linew)   -  temperature in middle of ice
!                               (t1 in mellor ...)
!            enthalpi(i,j,linew) -  scaled perturbation ice heat content
!            tis(i,j)        -  temperature at snow/atmos. interface
!                               (t3 in Mellor..)
!            brnfr(i,j)      -  brine fraction
!            wsm(i,j)        -  snow melting rate
!            wai(i,j)        -  melt rate at atmos./ice
!            apond(i,j,linew)-  melt water fraction
!            hpond(i,j,linew)-  melt water depth
!            ageice(i,j,linew)- ice age
!
!
!            ----  initiated in this routine ----
!
!            qai(i,j)        -  heat flux atmosphere/ice
!                               (positive from ice to atm.)
!            qio(i,j)        -  heat flux ice/oceam (possitive from ocean)
!            hfus1(i,j)      -  heat of fusion (L_o or L_3)
!            wro(i,j)        -  production rate of surface runoff
!            t2(i,j)         -  temperature at ice/snow interface
!
!***********************************************************************

      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE mod_boundary
!
      USE i2d_bc_mod
      USE tibc_mod, ONLY : tibc_tile
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, liold, linew
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
# endif
# ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
# endif
# ifdef AICLM_NUDGING
      real(r8), intent(in) :: aiclm(LBi:,LBj:)
      real(r8), intent(in) :: hiclm(LBi:,LBj:)
      real(r8), intent(in) :: AInudgcof(LBi:,LBj:)
# endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: wfr(LBi:,LBj:)
      real(r8), intent(inout) :: wai(LBi:,LBj:)
      real(r8), intent(inout) :: wao(LBi:,LBj:)
      real(r8), intent(inout) :: wio(LBi:,LBj:)
      real(r8), intent(inout) :: wro(LBi:,LBj:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hsn(LBi:,LBj:,:)
      real(r8), intent(inout) :: ageice(LBi:,LBj:,:)
#ifdef MELT_PONDS
      real(r8), intent(inout) :: apond(LBi:,LBj:,:)
      real(r8), intent(inout) :: hpond(LBi:,LBj:,:)
#endif
      real(r8), intent(inout) :: tis(LBi:,LBj:)
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: enthalpi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hage(LBi:,LBj:,:)
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: coef_ice_heat(LBi:,LBj:)
      real(r8), intent(inout) :: rhs_ice_heat(LBi:,LBj:)
      real(r8), intent(inout) :: s0mk(LBi:,LBj:)
      real(r8), intent(inout) :: t0mk(LBi:,LBj:)
      real(r8), intent(out) :: io_mflux(LBi:,LBj:)
#if defined ICE_BIO && defined BERING_10K
      real(r8), intent(inout) :: IcePhL(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNO3(LBi:,LBj:,:)
      real(r8), intent(inout) :: IceNH4(LBi:,LBj:,:)
#endif
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: qai_n(LBi:,LBj:)
      real(r8), intent(in) :: qi_o_n(LBi:,LBj:)
      real(r8), intent(in) :: qao_n(LBi:,LBj:)
      real(r8), intent(in) :: snow_n(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)
      real(r8), intent(out) :: stflx(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
# endif
# ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
# endif
# ifdef AICLM_NUDGING
      real(r8), intent(in) :: aiclm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: hiclm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: AInudgcof(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: wfr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wai(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wao(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wio(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wro(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ageice(LBi:UBi,LBj:UBj,2)
#ifdef MELT_PONDS
      real(r8), intent(inout) :: apond(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hpond(LBi:UBi,LBj:UBj,2)
#endif
      real(r8), intent(inout) :: tis(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: enthalpi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hage(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: coef_ice_heat(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rhs_ice_heat(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: s0mk(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t0mk(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: io_mflux(LBi:UBi,LBj:UBj)
#if defined ICE_BIO && defined BERING_10K
      real(r8), intent(inout) :: IcePhL(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNO3(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: IceNH4(LBi:UBi,LBj:UBj,2)
# endif
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: qai_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: qi_o_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: qao_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: snow_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#endif

! Local variable definitions
!

      integer :: i, j
      integer :: iday, month, year
      real(r8) :: hour, yday

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: b2d

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: alph
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ws

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: temp_top
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: salt_top
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sice
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: brnfr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: hfus1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: qi2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: qai
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: qio
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wsm
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: utau
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dztop
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: coa
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: t2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: cht
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: chs
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ai_old
!      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: enthal

#ifdef AICLM_NUDGING
      real(r8) :: cff
#endif
      real(r8) :: tfrz
      real(r8) :: cot
      real(r8) :: ai_tmp
#ifdef MELT_PONDS
      real(r8) :: vpond
      real(r8) :: vpond_new
      real(r8) :: pond_r, apond_old
#endif
      real(r8) :: pmelt

      real(r8), parameter :: eps = 1.E-4_r8
      real(r8), parameter :: prt = 13._r8
      real(r8), parameter :: prs = 2432._r8
      real(r8), parameter :: tpr = 0.85_r8
      real(r8), parameter :: nu = 1.8E-6_r8
      real(r8), parameter :: z0ii = 0.02_r8
      real(r8), parameter :: kappa = 0.4_r8
      real(r8), parameter :: rhosw = 1026._r8           ! [kg m-3]
      real(r8), parameter :: frln = -0.0543_r8          ! [psu C-1]
      real(r8), parameter :: sice_ref = 3.2_r8          ! [psu]
      real(r8), parameter :: alphic = 2.034_r8          ! [W m-1 K-1]
      real(r8), parameter :: alphsn = 0.31_r8           ! [W m-1 K-1]
      real(r8), parameter :: hfus = 3.347E+5_r8         ! [J kg-1]
      real(r8), parameter :: cpi = 2093.0_r8            ! [J kg-1 K-1]
      real(r8), parameter :: cpw = 3990.0_r8            ! [J kg-1 K-1]
      real(r8), parameter :: rhocpr = 0.2442754E-6_r8   ! [m s2 K kg-1]
      real(r8), parameter :: ykf = 3.14_r8
#ifdef MELT_PONDS
      real(r8), parameter :: pond_Tp = -2.0_r8          ! [C]
      real(r8), parameter :: pond_delta = 0.8_r8
      real(r8), parameter :: pond_rmin = 0.15_r8
      real(r8), parameter :: pond_rmax = 0.7_r8
#endif

      real(r8) :: corfac
      real(r8) :: hicehinv  ! 1./(0.5*ice_thick)
      real(r8) :: z0
      real(r8) :: zdz0
      real(r8) :: rno
      real(r8) :: termt
      real(r8) :: terms
      real(r8) :: tfz
      real(r8) :: xtot
      real(r8) :: phi
      real(r8) :: d1
      real(r8) :: d2i
      real(r8) :: d3

      real(r8) :: fac_shflx

#ifdef ICE_SHOREFAST
      real(r8) :: hh
      real(r8) :: clear
      real(r8) :: fac_sf
#endif
#ifdef ICE_CONVSNOW
      real(r8) :: hstar
#endif

#include "set_bounds.h"

      CALL caldate(r_date, tdays(ng), year, yday, month, iday, hour)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          temp_top(i,j)=t(i,j,N(ng),nrhs,itemp)
          salt_top(i,j)=t(i,j,N(ng),nrhs,isalt)
          salt_top(i,j) = MIN(MAX(0.0_r8,salt_top(i,j)),40.0_r8)
          dztop(i,j)=z_w(i,j,N(ng))-z_r(i,j,N(ng))
          stflx(i,j,isalt) = stflx(i,j,isalt)*                          &
     &          MIN(MAX(t(i,j,N(ng),nrhs,isalt),0.0_r8),60.0_r8)
#  if defined WET_DRY && defined CASPIAN
          stflx(i,j,isalt) = stflx(i,j,isalt)*rmask_wet(i,j)
#  endif
        END DO
      END DO

      d1 = rho_air(ng) * spec_heat_air * trans_coeff
      d2i = rho_air(ng) * sublim_latent_heat * trans_coeff
      d3 = StefBo * ice_emiss

      DO j=Jstr,Jend
        DO i=Istr,Iend
          utau(i,j) = sqrt(sqrt(                                        &
     &             (0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2                &
     &           + (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2                &
     &                   )    )
          utau(i,j) = max(utau(i,j),1.E-4_r8)
        END DO
      END DO

!------------------------------------------------------!
!   Get incoming long and shortwave radiation
!------------------------------------------------------!
!
! *** all rho's 0n 1026 kg m-3. cp's on 4.186e+6 j m-3 c-1.
! *** sigma=5.67e-8 w m-2 k-4 : sigma=sigma*epsilon. m&u used 5.78e-8
! *** compute sw ,lw & back radiation (see p&w 1979)
!
!-----------------------------------------------------------------
!  calculate those parts of the energy balance which do not depend
!  on the surface temperature.
!-----------------------------------------------------------------
! *** ignore snow effects except change albi to albsn value

!      beregner sno- og is-tykkelse
!   (compute snow and ice thicknesses)
      DO j = Jstr,Jend
        DO i = Istr,Iend
          sice(i,j) = MIN(sice_ref,salt_top(i,j))
          ice_thick(i,j) = MAX(0.05_r8,                                 &
     &                hi(i,j,linew)/MAX(ai(i,j,linew),eps))
          snow_thick(i,j) = hsn(i,j,linew)/MAX(ai(i,j,linew),eps)
          ai_old(i,j) = ai(i,j,linew)
          brnfr(i,j) = frln*sice(i,j)/MIN(ti(i,j,linew),-eps)
          brnfr(i,j) = MIN(brnfr(i,j),0.2_r8)
          brnfr(i,j) = MAX(brnfr(i,j),0.0_r8)
!      alph - thermal conductivity of ice
          alph(i,j) = alphic*(1._r8-1.2_r8*brnfr(i,j))
#ifndef ICE_BOX
          corfac = 1._r8/(0.5_r8*(1._r8+EXP(-(hi(i,j,linew)/1._r8)**2)))
          alph(i,j) = alph(i,j)*corfac
#endif
          coa(i,j) = 2.0_r8*alph(i,j)*snow_thick(i,j)/                  &
     &                   (alphsn*ice_thick(i,j))
        END DO
      END DO

! *** compute ice thermodynamic variables
!*    specify snow fall rate and snow thickness
!*    compute net ice atmos. surface heat transfer
!*    zero if temp. is below freezing.
!     frysepunktspemp. (t=-0.27 c)



!-----------------------------------------------------------------------
!     SOLVE FOR TEMPERATURE AT THE TOP OF THE ICE LAYER
!-----------------------------------------------------------------------
      DO j = Jstr,Jend
        DO i = Istr,Iend
! gradient coefficient for heat conductivity part
          b2d(i,j) = 2.0_r8*alph(i,j)/(ice_thick(i,j)*(1._r8+coa(i,j)))
          coef_ice_heat(i,j) = coef_ice_heat(i,j) + b2d(i,j)

          IF (ai(i,j,linew) .gt. min_a(ng)) THEN

! downward conductivity term, assuming the ocean at the freezing point
            rhs_ice_heat(i,j) = rhs_ice_heat(i,j) +                     &
     &              b2d(i,j)*ti(i,j,linew)
            tis(i,j) = rhs_ice_heat(i,j)/coef_ice_heat(i,j)
            tis(i,j) = MAX(tis(i,j),-45._r8)
            qai(i,j) = qai_n(i,j)
          ELSE
            tis(i,j) = temp_top(i,j)
            qai(i,j) = qai_n(i,j)
          END IF
        END DO
      END DO

      DO j = Jstr,Jend
        DO i = Istr,Iend
!**** calculate interior ice temp and heat fluxes
!       new temperature in ice
          IF (ai(i,j,linew) .gt. min_a(ng)) THEN
            cot = cpi - frln*sice(i,j)*hfus/(ti(i,j,linew)-eps)**2
!            enthal(i,j,1) = brnfr(i,j) * (hfus + cpw*ti(i,j,linew)) +  &
!     &                (1 - brnfr(i,j)) * cpi * ti(i,j,linew)
#ifdef ICE_I_O
            ti(i,j,linew) = ti(i,j,linew) +                             &
     &               dtice(ng)/(rhoice(ng)*ice_thick(i,j)*cot)*         &
     &        (2._r8*alph(i,j)/ice_thick(i,j)*                          &
     &         (t0mk(i,j) + (tis(i,j) - (2._r8+coa(i,j))*ti(i,j,linew)) &
     &                      /(1._r8+coa(i,j))) + qi_o_n(i,j))
#else
            ti(i,j,linew) = ti(i,j,linew) + dtice(ng)*(                 &
     &        2._r8*alph(i,j)/(rhoice(ng)*ice_thick(i,j)**2*cot)        &
     &        *(t0mk(i,j) + (tis(i,j) - (2._r8+coa(i,j))*ti(i,j,linew)) &
     &                                        /(1._r8+coa(i,j))))
#endif
            ti(i,j,linew) = max(ti(i,j,linew),-35._r8)
            ti(i,j,linew) = min(ti(i,j,linew),-eps)
!            brnfr(i,j) = frln*sice(i,j)/MIN(ti(i,j,linew),-eps)
!            enthal(i,j,2) = brnfr(i,j) * (hfus + cpw*ti(i,j,linew)) +  &
!     &                (1 - brnfr(i,j)) * cpi * ti(i,j,linew)
          ELSE
            ti(i,j,linew) = temp_top(i,j)
          END IF
        END DO
      END DO

      DO j = Jstr,Jend
        DO i = Istr,Iend
          IF (ai(i,j,linew) .gt. min_a(ng)) THEN
            t2(i,j) = (tis(i,j)+coa(i,j)*ti(i,j,linew))/(1._r8+coa(i,j))
            hicehinv = 2._r8/ice_thick(i,j)
            qi2(i,j) = alph(i,j)*(ti(i,j,linew)-t2(i,j))*hicehinv
            qio(i,j) = alph(i,j)*(t0mk(i,j)-ti(i,j,linew))*hicehinv
          END IF

!  Compute net heat flux from ice to atmosphere - Mellor and Kantha (7)
        END DO
      END DO

!**** for open water ice fluxes set to zero
      DO j = Jstr,Jend
        DO i = Istr,Iend
          IF (ai(i,j,linew) .le. min_a(ng)) THEN
#ifdef MASKING
# ifdef WET_DRY
            tis(i,j) = t0mk(i,j)*rmask(i,j)*rmask_wet(i,j)
            t2(i,j) = t0mk(i,j)*rmask(i,j)*rmask_wet(i,j)
            ti(i,j,linew) = -2.0_r8*rmask(i,j)*rmask_wet(i,j)
# else
            tis(i,j) = t0mk(i,j)*rmask(i,j)
            t2(i,j) = t0mk(i,j)*rmask(i,j)
            ti(i,j,linew) = -2.0_r8*rmask(i,j)
# endif
#elif defined WET_DRY
            tis(i,j) = t0mk(i,j)*rmask_wet(i,j)
            t2(i,j) = t0mk(i,j)*rmask_wet(i,j)
            ti(i,j,linew) = -2.0_r8*rmask_wet(i,j)
#else
            tis(i,j) = t0mk(i,j)
            t2(i,j) = t0mk(i,j)
            ti(i,j,linew) = -2.0_r8
#endif
#ifdef ICESHELF
            IF (zice(i,j).ne.0.0_r8) THEN
              tis(i,j) = 0.0_r8
              t2(i,j) = 0.0_r8
              ti(i,j,linew) = 0.0_r8
            END IF
#endif
            qi2(i,j) = 0._r8
            qai(i,j) = 0._r8
            qio(i,j) = 0._r8
            hsn(i,j,linew) = 0._r8
#ifdef MELT_PONDS
            apond(i,j,linew) = 0._r8
            hpond(i,j,linew) = 0._r8
#endif
          END IF
        END DO
      END DO

! Set snow fall rate to value derived from precipitation rate

      DO j = Jstr,Jend
        DO i = Istr,Iend
          snow(i,j) = max(snow_n(i,j),0._r8)
          ws(i,j) = snow(i,j)
        END DO
      END DO

      DO j = Jstr,Jend
        DO i = Istr,Iend
          tfrz = frln*sice(i,j)
          wsm(i,j) = 0._r8
          wai(i,j) = 0._r8
          wro(i,j) = 0._r8
          IF (ai(i,j,linew) .gt. min_a(ng)) THEN
! Melt ice or freeze surface water in the fall if there is no snow
            IF (hsn(i,j,linew) .le. 0.0_r8) THEN
#ifdef MELT_PONDS
              IF (tis(i,j) .gt. tfrz .or. hpond(i,j,linew) .gt. 0._r8)  &
     &                                                            THEN
#else
              IF (tis(i,j) .gt. tfrz) THEN
#endif
!   ice warmer than freezing point
                tis(i,j) = tfrz
		t2(i,j) = tfrz
!   ice warmer than freezing point
                hfus1(i,j) = hfus*(1._r8-brnfr(i,j))+tis(i,j)*cpw       &
     &         -((1._r8-brnfr(i,j))*cpi+brnfr(i,j)*cpw)*ti(i,j,linew)
                qai(i,j) = qai_n(i,j)
                qi2(i,j) = b2d(i,j)*(ti(i,j,linew)-tis(i,j))
!    compute ice production rate (negative here) from atmosphere-ice exchange
!    Means wai is positive for melt
                wai(i,j) = -(qai(i,j)-qi2(i,j)) /(hfus1(i,j)*rhosw)
!    compute production rate for melt water (melting rate)
                wsm(i,j) = ws(i,j)
              END IF
            ELSE
!      there is snow cover
!***** to melt snow or
!***********  freeze surface water under snow
              IF (tis(i,j) .gt. 0.0_r8) THEN
!          ice temperature warmer than the freezing point
                tis(i,j) = 0._r8
                qai(i,j) = qai_n(i,j)
                qi2(i,j) = b2d(i,j)*(ti(i,j,linew)-tis(i,j))
                t2(i,j) = (tis(i,j)+coa(i,j)*ti(i,j,linew))/            &
     &                (1._r8+coa(i,j))
!          snow melting
! When does snow get denser???
                wsm(i,j) = max(0.0_r8,-(qai(i,j)-qi2(i,j))/             &
     &                    (rhosnow_dry(ng)*hfus)) + ws(i,j)
!     &                    (rhosnow_wet(ng)*hfus)) + ws(i,j)
              END IF

#ifdef MELT_PONDS
              IF (tis(i,j) < 0.0_r8 .and.                               &
     &            hpond(i,j,linew) > 0.0_r8) THEN
!          colder than the freezing point
                tis(i,j) = 0._r8
                qai(i,j) = qai_n(i,j)
                qi2(i,j) = b2d(i,j)*(ti(i,j,linew)-tis(i,j))
                wai(i,j) = -(qai(i,j)-qi2(i,j))/(hfus*rhosw)
              END IF
#endif
            END IF
!
!***** compute snow thickness
!       hsn - snow thickness
#ifdef NO_SNOW
            hsn(i,j,linew) = 0.0_r8
#else
            hsn(i,j,linew) = hsn(i,j,linew)+(ai(i,j,linew)              &
     &                             *(-wsm(i,j)+ws(i,j)))*dtice(ng)
            hsn(i,j,linew) = max(0.0_r8,hsn(i,j,linew))
#endif
          END IF
#ifdef MELT_PONDS
!
!  Update melt ponds
!  This all comes from CICE's cesm pond scheme.
!
          IF (ai(i,j,linew) > min_a(ng)) THEN
            vpond = apond(i,j,linew)*hpond(i,j,linew)*ai(i,j,linew)
!  pond growth (should have rain...)
            pmelt = MAX(0._r8,wai(i,j)+wsm(i,j))
            pond_r = pond_rmin+(pond_rmax-pond_rmin)*ai(i,j,linew)
            vpond = vpond + pmelt*pond_r*dtice(ng)
            wro(i,j) = (1.0_r8-pond_r)*pmelt
!  pond contraction
            vpond = vpond*exp(0.01_r8*MAX((pond_Tp-tis(i,j)),0._r8)/    &
     &                         pond_Tp)
!  New pond shape
            apond(i,j,linew) = MIN(1.0_r8,                              &
     &                   sqrt(vpond/(pond_delta*ai(i,j,linew))))
            hpond(i,j,linew) = pond_delta*apond(i,j,linew)
            IF (hi(i,j,linew) < 0.01_r8) THEN
              hpond(i,j,linew) = 0.0_r8
              apond(i,j,linew) = 0.0_r8
            ELSE IF (hpond(i,j,linew) .gt. 0.9_r8*hi(i,j,linew)) THEN
              hpond(i,j,linew) = 0.9_r8*hi(i,j,linew)
              apond(i,j,linew) = hi(i,j,linew)/pond_delta
            END IF
            vpond_new = apond(i,j,linew)*hpond(i,j,linew)*ai(i,j,linew)
            wro(i,j) = wro(i,j) + (vpond-vpond_new)/dtice(ng)
          ELSE
            vpond = apond(i,j,linew)*hpond(i,j,linew)*ai(i,j,linew)
            wro(i,j) = vpond/dtice(ng)
            apond(i,j,linew) = 0.0_r8
            hpond(i,j,linew) = 0.0_r8
          END IF
#else
          pmelt = MAX(0._r8,wai(i,j)+wsm(i,j))
          wro(i,j) = pmelt
#endif
        END DO
      END DO

      DO j = Jstr,Jend
        DO i = Istr,Iend
          z0 = max(z0ii*ice_thick(i,j),0.01_r8)
          z0 = min(z0,0.1_r8)
!
!     *** Yaglom and Kader formulation for z0t and z0s
!
          zdz0 = dztop(i,j)/z0   !WPB
          zdz0 = MAX(zdz0,3._r8)

          rno = utau(i,j)*0.09_r8/nu
          termt = ykf*sqrt(rno)*prt**0.666667_r8
          terms = ykf*sqrt(rno)*prs**0.666667_r8
          cht(i,j) = utau(i,j)/(tpr*log(zdz0)/kappa+termt)
          chs(i,j) = utau(i,j)/(tpr*log(zdz0)/kappa+terms)
        END DO
      END DO

      DO j = Jstr,Jend
        DO i = Istr,Iend
          tfz = frln*salt_top(i,j)
          wao(i,j) = 0._r8
          wio(i,j) = 0._r8

          hfus1(i,j) = hfus*(1.0_r8-brnfr(i,j))+t0mk(i,j)*cpw           &
     &         -((1.0_r8-brnfr(i,j))*cpi+brnfr(i,j)*cpw)*ti(i,j,linew)
          IF (temp_top(i,j) .le. tfz)                                   &
     &             wao(i,j) = qao_n(i,j)/(hfus1(i,j)*rhosw)
          IF (ai(i,j,linew) .le. min_a(ng) .or.                         &
     &        hi(i,j,linew) .le. min_h(ng)) THEN
            s0mk(i,j) = salt_top(i,j)
            t0mk(i,j) = temp_top(i,j)
            wai(i,j) = 0._r8
            xtot = (1._r8-ai(i,j,linew))*wao(i,j)
          ELSE

! MK89 version
#ifdef ICE_BOX
!  F_t set to 2 W/m^2
            wio(i,j) = (qio(i,j) - 2.0_r8)/(rhosw*hfus1(i,j))
            xtot = ai(i,j,linew)*wio(i,j)                               &
     &              +(1._r8-ai(i,j,linew))*wao(i,j)
#else
            wio(i,j) = (qio(i,j)/rhosw +                                &
     &               cpw*cht(i,j)*(t0mk(i,j)-temp_top(i,j)))/hfus1(i,j)

            xtot = ai(i,j,linew)*wio(i,j)                               &
     &              +(1._r8-ai(i,j,linew))*wao(i,j)

            s0mk(i,j) =                                                 &
     &            (chs(i,j)*salt_top(i,j)+(wro(i,j)-wio(i,j))*sice(i,j))&
     &                /(chs(i,j)+wro(i,j)-wio(i,j))
            s0mk(i,j) = max(s0mk(i,j),0._r8)
            s0mk(i,j) = min(s0mk(i,j),40._r8)

            t0mk(i,j) = frln*s0mk(i,j)
#endif

          END IF

#ifdef CASPIAN_XXX
          hh = h(i,j)+Zt_avg1(i,j)
          IF (hh.LT.1.0_r8) THEN
            fac_shflx = hh
!            fac_shflx = 0.0_r8
          ELSE
            fac_shflx = 1.0_r8
          END IF
#else
          fac_shflx = 1.0_r8
#endif
#ifdef ICESHELF
          IF (zice(i,j).eq.0.0_r8) THEN
#endif
            IF(ai(i,j,linew).LE.min_a(ng)) THEN
               stflx(i,j,itemp) = qao_n(i,j)*fac_shflx
            ELSE
#ifdef ICE_SHOREFAST
              hh = h(i,j)+Zt_avg1(i,j)
              clear = hh-0.9_r8*hi(i,j,liold)
              clear = MAX(clear,0.0_r8)
              IF (clear.lt.1.5_r8) THEN
                fac_sf = MAX(clear-0.5_r8,0.0_r8)/1.0_r8
              ELSE
                fac_sf = 1.0_r8
              END IF
              stflx(i,j,itemp) = (1.0_r8-ai(i,j,linew))*qao_n(i,j)      &
     &                          *fac_shflx                              &
     &                   +(ai(i,j,linew)*qio(i,j)                       &
     &                   -xtot*hfus1(i,j))*fac_sf
#else
              stflx(i,j,itemp) = (1.0_r8-ai(i,j,linew))*qao_n(i,j)      &
     &                   +ai(i,j,linew)*qio(i,j)                        &
     &                   -xtot*hfus1(i,j)
#endif
#if defined WET_DRY && defined CASPIAN
          stflx(i,j,itemp) = stflx(i,j,itemp)*rmask_wet(i,j)
#endif
            END IF

! Change stflx(i,j,itemp) back to ROMS convention
            stflx(i,j,itemp) = -stflx(i,j,itemp) * rhocpr

#ifdef MASKING
            stflx(i,j,itemp) = stflx(i,j,itemp)*rmask(i,j)
#endif
#ifdef WET_DRY
!            stflx(i,j,itemp) = stflx(i,j,itemp)*rmask_wet(i,j)
#endif
#ifdef ICE_SHOREFAST
            stflx(i,j,isalt) = stflx(i,j,isalt) +                       &
     &        (- (xtot-ai(i,j,linew)*wro(i,j))*                         &
     &          (sice(i,j)-MIN(MAX(s0mk(i,j),0.0_r8),60.0_r8)) )        &
     &                     *fac_sf
#else
            stflx(i,j,isalt) = stflx(i,j,isalt)                         &
     &          - (xtot-ai(i,j,linew)*wro(i,j))*(sice(i,j)-s0mk(i,j))
#endif

! Test for case of rainfall on snow/ice and assume 100% drainage
#ifndef NCEP_FLUXES
            IF (rain(i,j).gt.0._r8.AND.snow_n(i,j).EQ.0._r8) THEN
              stflx(i,j,isalt) = stflx(i,j,isalt) -                       &
     &                           ai(i,j,linew)*rain(i,j)*0.001_r8
            END IF
#endif
!  io_mflux is ice production rate (+ve for growth)
            io_mflux(i,j) = xtot - ai(i,j,linew)*wro(i,j) + wfr(i,j)
#ifdef MASKING
            stflx(i,j,isalt) = stflx(i,j,isalt)*rmask(i,j)
            io_mflux(i,j) = io_mflux(i,j)*rmask(i,j)
#endif
#ifdef WET_DRY
            stflx(i,j,isalt) = stflx(i,j,isalt)*rmask_wet(i,j)
            io_mflux(i,j) = io_mflux(i,j)*rmask_wet(i,j)
#endif
#ifdef ICESHELF
          ELSE
            io_mflux(i,j) = 0.0_r8
          END IF
#endif
        END DO
      END DO

!********************************

      DO j = Jstr,Jend
        DO i = Istr,Iend
          phi = 4._r8
          if (wao(i,j) .lt. 0.0_r8 ) phi = 0.5_r8
          hi(i,j,linew) = hi(i,j,linew)+dtice(ng)                       &
     &             *(ai(i,j,linew)                                      &
     &             *(wio(i,j)-wai(i,j))                                 &
     &        +(1.0_r8-ai(i,j,linew))*wao(i,j) + wfr(i,j))

          ai_tmp = ai(i,j,linew)
          ai(i,j,linew) = ai(i,j,linew) +                               &
     &             dtice(ng)*(1.0_r8-ai(i,j,linew))                     &
     &                      *(phi*wao(i,j)+wfr(i,j))
          ai(i,j,linew) = min(ai(i,j,linew),max_a(ng))

#ifndef NO_SNOW
! adjust snow volume when ice melting out from under it
          IF (ai(i,j,linew) .lt. ai_tmp)                                &
     &        hsn(i,j,linew) =                                          &
     &           hsn(i,j,linew)*ai(i,j,linew)/max(ai_tmp,eps)

# ifdef ICE_CONVSNOW
!
! If snow base is below sea level, then raise the snow base to sea level
!  by converting some snow to ice (N.B. hstar is also weighted by ai
!  like hsn and hi)
!
          hstar = hsn(i,j,linew) - (rhosw - rhoice(ng)) *               &
     &             hi(i,j,linew) / rhosnow_dry(ng)
          IF (hstar .gt. 0.0_r8) THEN
            hsn(i,j,linew) = hsn(i,j,linew) - rhoice(ng)*hstar/rhosw
            hi(i,j,linew) = hi(i,j,linew) + rhosnow_dry(ng)*hstar/rhosw
          ENDIF
# endif
#endif
#ifdef AICLM_NUDGING
          cff = AInudgcof(i,j)
          ai(i,j,linew)=ai(i,j,linew)+                                  &
     &                  dtice(ng)*cff*(aiclm(i,j)-ai(i,j,linew))
          hi(i,j,linew)=hi(i,j,linew)+                                  &
     &                  dtice(ng)*cff*(hiclm(i,j)-hi(i,j,linew))
#endif
#ifdef MASKING
          ai(i,j,linew) = ai(i,j,linew)*rmask(i,j)
          hi(i,j,linew) = hi(i,j,linew)*rmask(i,j)
#endif
#ifdef WET_DRY
!          ai(i,j,linew) = ai(i,j,linew)*rmask_wet(i,j)
!          hi(i,j,linew) = hi(i,j,linew)*rmask_wet(i,j)
#endif
#ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) THEN
            ai(i,j,linew) = 0.0_r8
            hi(i,j,linew) = 0.0_r8
          END IF
#endif
#ifdef MELT_PONDS
!
! Adjust ponds for changes in ai
!
          IF (ai(i,j,linew) > min_a(ng) .and.                           &
     &                       apond(i,j,linew) > 0._r8) THEN
            apond_old = apond(i,j,linew)
            apond(i,j,linew) = apond(i,j,linew)*ai_old(i,j)/            &
     &                      ai(i,j,linew)
            hpond(i,j,linew) = hpond(i,j,linew)*apond_old*ai_old(i,j)   &
     &                 /(ai(i,j,linew)*apond(i,j,linew))
          END IF
#endif

! determine age of the sea ice
! Case 1 - new ice
          IF (ageice(i,j,linew).le.0.0_r8                               &
     &                .and.hi(i,j,linew).gt.min_h(ng)) THEN
            ageice(i,j,linew)=dtice(ng)/86400._r8
! Case 2 - existing ice gets older
          ELSEIF(ageice(i,j,linew).gt.0.0_r8                            &
     &                .and.hi(i,j,linew).gt.min_h(ng)) THEN
            ageice(i,j,linew) = ageice(i,j,linew)+dtice(ng)/86400._r8
! Case 3 - all ice in cell has melted or is open water and stays open water
          ELSE
            ageice(i,j,linew) = 0.0_r8
          ENDIF

#undef DIAG_WPB
#ifdef DIAG_WPB
      IF (i.eq.1.and.j.eq.1) THEN
         write(*,*) tdays,wio(i,j),wai(i,j),wao(i,j),wfr(i,j),          &
     &              ai(i,j,linew),hi(i,j,linew),tis(i,j),               &
#ifdef MELT_PONDS
     &              apond(i,j,linew), hpond(i,j,linew),                 &
#endif
     &              temp_top(i,j),t0mk(i,j),stflx(i,j,itemp),           &
     &              salt_top(i,j),s0mk(i,j),stflx(i,j,isalt),           &
     &              qio(i,j), ti(i,j,linew), brnfr(i,j),                &
     &              t2(i,j), qao_n(i,j), qi2(i,j), qai_n(i,j)
        print *
      END IF
#endif
#ifdef ICE_BOX
!      IF (i.eq.1.and.j.eq.1) THEN
!         write(*,*) tdays,enthal(i,j,1),enthal(i,j,2),                  &
!     &              hi(i,j,linew),hsn(i,j,linew),tis(i,j),              &
!     &              ti(i,j,linew), t2(i,j),                             &
!     &              qio(i,j), qi2(i,j), qi_o_n(i,j),                    &
!     &              (qio(i,j) - qi2(i,j) + qi_o_n(i,j))*dtice(ng)/      &
!     &              (ice_thick(i,j)*rhoice(ng))
!        print *
!      END IF
      IF (i.eq.1.and.j.eq.1.and.iday==15.and.int(hour)==0) THEN
         write(*,*) tdays,wio(i,j),wai(i,j),wro(i,j),                   &
     &              hi(i,j,linew),hsn(i,j,linew),tis(i,j),              &
#ifdef MELT_PONDS
     &              apond(i,j,linew), hpond(i,j,linew),                 &
#endif
     &              ti(i,j,linew), t2(i,j),                             &
     &              qio(i,j), qi2(i,j), qai_n(i,j),                     &
     &              alph(i,j), coa(i,j), qi_o_n(i,j), cot, t0mk(i,j)
        print *
      END IF
#endif

        ENDDO
      ENDDO

!********************************
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ai(i,j,linew) = MIN(ai(i,j,linew),max_a(ng))
          ai(i,j,linew) = MAX(ai(i,j,linew),0.0_r8)
          hi(i,j,linew) = MAX(hi(i,j,linew),0.0_r8)
          hsn(i,j,linew) = MAX(hsn(i,j,linew),0.0_r8)
          ti(i,j,linew) = MAX(ti(i,j,linew),-70.0_r8)
          if (hi(i,j,linew) .le. 0.0_r8) ai(i,j,linew) = 0.0_r8
          if (ai(i,j,linew) .le. 0.0_r8) hi(i,j,linew) = 0.0_r8
        ENDDO
      ENDDO

      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tis)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  coef_ice_heat)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rhs_ice_heat)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,isalt))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,itemp))

      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold, linew,                                   &
     &                  BOUNDARY(ng)%ai_west(LBj:UBj),                  &
     &                  BOUNDARY(ng)%ai_east(LBj:UBj),                  &
     &                  BOUNDARY(ng)%ai_north(LBi:UBi),                 &
     &                  BOUNDARY(ng)%ai_south(LBi:UBi),                 &
     &                  ui, vi, ai, LBC(:,isAice,ng))
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold, linew,                                   &
     &                  BOUNDARY(ng)%hi_west(LBj:UBj),                  &
     &                  BOUNDARY(ng)%hi_east(LBj:UBj),                  &
     &                  BOUNDARY(ng)%hi_north(LBi:UBi),                 &
     &                  BOUNDARY(ng)%hi_south(LBi:UBi),                 &
     &                  ui, vi, hi, LBC(:,isHice,ng))
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold, linew,                                   &
     &                  BOUNDARY(ng)%hsn_west(LBj:UBj),                 &
     &                  BOUNDARY(ng)%hsn_east(LBj:UBj),                 &
     &                  BOUNDARY(ng)%hsn_north(LBi:UBi),                &
     &                  BOUNDARY(ng)%hsn_south(LBi:UBi),                &
     &                  ui, vi, hsn, LBC(:,isHsno,ng))
      CALL tibc_tile (ng, tile, iNLM,                                   &
     &                          LBi, UBi, LBj, UBj, liold, linew,       &
     &                          ui, vi, hi, ti, enthalpi)
#ifdef MELT_PONDS
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold, linew,                                   &
     &                  BOUNDARY(ng)%apond_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%apond_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%apond_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%apond_south(LBi:UBi),              &
     &                  ui, vi, apond, LBC(:,isApond,ng))
      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold, linew,                                   &
     &                  BOUNDARY(ng)%hpond_west(LBj:UBj),               &
     &                  BOUNDARY(ng)%hpond_east(LBj:UBj),               &
     &                  BOUNDARY(ng)%hpond_north(LBi:UBi),              &
     &                  BOUNDARY(ng)%hpond_south(LBi:UBi),              &
     &                  ui, vi, hpond, LBC(:,isHpond,ng))
#endif
!      CALL i2d_bc_tile (ng, tile, iNLM,                                 &
!     &                  LBi, UBi, LBj, UBj,                             &
!     &                  IminS, ImaxS, JminS, JmaxS,                     &
!     &                  liold, linew,                                   &
!     &                  BOUNDARY(ng)%ageice_west,                       &
!     &                  BOUNDARY(ng)%ageice_east,                       &
!     &                  BOUNDARY(ng)%ageice_north,                      &
!     &                  BOUNDARY(ng)%ageice_south,                      &
!     &                  ui, vi, ageice, LBC(:,isAgeice,ng))
!     CALL ageicebc_tile (ng, tile,                                     &
!    &                          LBi, UBi, LBj, UBj, liold, linew,       &
!    &                          min_h(ng), ui, vi, hi, ageice, hage)
#if defined ICE_BIO && defined BERING_10K
FOOO
! Convert these too.
      CALL IcePhLbc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, IcePhL)
      CALL IceNO3bc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, IceNO3)
      CALL IceNH4bc_tile (ng, tile,                                     &
     &                LBi, UBi, LBj, UBj,                               &
     &                liold, linew,                                     &
     &                ui, vi, IceNH4)
#endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ai(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hi(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hsn(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ti(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          enthalpi(:,:,linew))
#ifdef MELT_PONDS
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          apond(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hpond(:,:,linew))
#endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ageice(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hage(:,:,linew))
# if defined ICE_BIO && defined BERING_10K
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IcePhL(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IceNO3(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IceNH4(:,:,linew))
# endif
      END IF
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ai(:,:,linew), hi(:,:,linew),                 &
     &                    hsn(:,:,linew), ti(:,:,linew))
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    enthalpi(:,:,linew),                          &
     &                    ageice(:,:,linew), hage(:,:,linew))
# ifdef MELT_PONDS
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    apond(:,:,linew), hpond(:,:,linew))
# endif
# if defined ICE_BIO && defined BERING_10K
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    IcePhL(:,:,linew), IceNO3(:,:,linew),         &
     &                    IceNH4(:,:,linew))

# endif
#endif

      RETURN
      END SUBROUTINE ice_thermo_tile
