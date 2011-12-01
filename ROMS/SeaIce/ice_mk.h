     SUBROUTINE ice_thermo (ng, tile)
!
!*************************************************** W. Paul Budgell ***
!  Copyright (c) 2009 ROMS/TOMS Group                                  !
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

      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"

# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 44)
# endif

      CALL ice_thermo_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), liold(ng), linew(ng),             &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      OCEAN(ng) % t,                              &
     &                      ICE(ng) % wfr,                              &
     &                      ICE(ng) % wai,                              &
     &                      ICE(ng) % wao,                              &
     &                      ICE(ng) % wio,                              &
     &                      ICE(ng) % wro,                              &
     &                      ICE(ng) % ai,                               &
     &                      ICE(ng) % hi,                               &
     &                      ICE(ng) % hsn,                              &
     &                      ICE(ng) % sfwat,                            &
     &                      ICE(ng) % ageice,                           &
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
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
     &                      FORCES(ng) % qai_n,                         &
     &                      FORCES(ng) % qao_n,                         &
     &                      FORCES(ng) % p_e_n,                         &
     &                      FORCES(ng) % snow_n,                        &
     &                      FORCES(ng) % rain,                          &
     &                      FORCES(ng) % stflx)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 44)
# endif
      RETURN
      END SUBROUTINE ice_thermo
!
!***********************************************************************
      SUBROUTINE ice_thermo_tile (ng, tile,                             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs, liold, linew,                       &
# ifdef MASKING
     &                        rmask,                                    &
# endif
# ifdef WET_DRY
     &                        rmask_wet,                                &
# endif
# ifdef ICESHELF
     &                        zice,                                     &
# endif
     &                        z_r, z_w, pm, pn, t,                      &
     &                        wfr, wai, wao, wio, wro,                  &
     &                        ai, hi, hsn, sfwat, ageice, tis, ti,      &
     &                        enthalpi, hage,                           &
     &                        ui, vi, coef_ice_heat, rhs_ice_heat,      &
     &                        s0mk, t0mk, io_mflux,                     &
     &                        sustr, svstr,                             &
     &                        qai_n, qao_n,                             &
     &                        p_e_n,                                    &
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
!        sfcalb
!        sfwat
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
!            wsm(i,j)        -  melting rate
!            wai(i,j)        -  production rate at atmos./ice
!            sfwat(i,j,linew)-  melt water
!            ageice(i,j,linew)- ice age
!
!
!            ----  initiated in this routine ----
!
!            qai(i,j)        -  heat flux atmosphere/ice
!                               (positive from ice to atm.)
!            qio(i,j)        -  heat flux ice/oceam (possitive from ocean)
!            sfcalb(i,j)     -  global albedo
!            hfus1(i,j)      -  heat of fusion (L_o or L_3)
!            wro(i,j)        -  production rate of surface runoff
!            t2(i,j)         -  temperature at ice/snow interface
!            hsn(i,j,linew)  -  snow avg. thickness next time step
!
!
!***********************************************************************

      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
!
      USE aibc_mod, ONLY : aibc_tile
      USE hibc_mod, ONLY : hibc_tile
      USE hsnbc_mod, ONLY : hsnbc_tile
      USE tibc_mod, ONLY : tibc_tile
      USE sfwatbc_mod, ONLY : sfwatbc_tile
      USE ageicebc_mod, ONLY : ageicebc_tile
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      implicit none

!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, liold, linew
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
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: wfr(LBi:,LBj:)
      real(r8), intent(inout) :: wai(LBi:,LBj:)
      real(r8), intent(inout) :: wao(LBi:,LBj:)
      real(r8), intent(inout) :: wio(LBi:,LBj:)
      real(r8), intent(inout) :: wro(LBi:,LBj:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hsn(LBi:,LBj:,:)
      real(r8), intent(inout) :: sfwat(LBi:,LBj:,:)
      real(r8), intent(inout) :: ageice(LBi:,LBj:,:)
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
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: qai_n(LBi:,LBj:)
      real(r8), intent(in) :: qao_n(LBi:,LBj:)
      real(r8), intent(in) :: p_e_n(LBi:,LBj:)
      real(r8), intent(in) :: snow_n(LBi:,LBj:)
      real(r8), intent(in) :: rain(LBi:,LBj:)
      real(r8), intent(out) :: stflx(LBi:,LBj:,:)
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
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: wfr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wai(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wao(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wio(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: wro(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: sfwat(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ageice(LBi:UBi,LBj:UBj,2)
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
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: qai_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: qao_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: p_e_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: snow_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
# endif

! Local variable definitions
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif

      integer :: i, j, it

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: b2d
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: c2d
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: f2d
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: g2d
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: h2d

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: lathi
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: alph
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ws
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: qa

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: temp_top
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: salt_top
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
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: w0
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: cht
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: chs


      real(r8) :: tfrz
      real(r8) :: cot
      real(r8) :: xmelt
      real(r8) :: ai_tmp

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
      real(r8), parameter :: sfwatmx = 0.1_r8           ! [m]
      real(r8), parameter :: rhocpr = 0.2442754E-6_r8   ! [m s2 K kg-1]

      real(r8) :: sice
      real(r8) :: thic
      real(r8) :: corfac
      real(r8) :: hicehinv  ! 1./(0.5*thic)
      real(r8) :: z0
      real(r8) :: zdz0
      real(r8) :: rno
      real(r8) :: termt
      real(r8) :: terms
      real(r8) :: tfz
      real(r8) :: xwai
      real(r8) :: xtot
      real(r8) :: phi
      real(r8) :: ykf
      real(r8) :: d1
      real(r8) :: d2i
      real(r8) :: d3

# include "set_bounds.h"

      DO j=Jstr,Jend
        DO i=Istr,Iend
          temp_top(i,j)=t(i,j,N(ng),nrhs,itemp)
          salt_top(i,j)=t(i,j,N(ng),nrhs,isalt)
          salt_top(i,j) = MIN(MAX(0.0_r8,salt_top(i,j)),40.0_r8)
          dztop(i,j)=z_w(i,j,N(ng))-z_r(i,j,N(ng))
          stflx(i,j,isalt) = stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)
        END DO
      END DO

      ykf = 3.14_r8
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
! *** all rho's 0n 1026 kg m-3. cp's on 4.186e+6 j m-3 c-1.
! *** sigma=5.67e-8 w m-2 k-4 : sigma=sigma*epsilon. m&u used 5.78e-8
! *** ignore snow effects except change albi to albsn value

!      beregner sno- og is-tykkelse
!   (compute snow and ice thicknesses)
      DO j = Jstr,Jend
        DO i = Istr,Iend
          sice = MIN(sice_ref,salt_top(i,j))
          ice_thick(i,j) = 0.05_r8+hi(i,j,linew)/                       &
     &                    (ai(i,j,linew)+eps)
          snow_thick(i,j) = hsn(i,j,linew)/(ai(i,j,linew)+eps)
          brnfr(i,j) = frln*sice/(ti(i,j,linew)-eps)
          brnfr(i,j) = min(brnfr(i,j),0.2_r8)
          brnfr(i,j) = max(brnfr(i,j),0.0_r8)
!      alph - thermal conductivity of ice
          alph(i,j) = alphic*MAX(1._r8-1.2_r8*brnfr(i,j),0.25_r8)
          corfac = 1._r8/(0.5_r8*(1._r8+EXP(-(hi(i,j,linew)/1._r8)**2)))
          alph(i,j) = alph(i,j)*corfac
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
          coef_ice_heat(i,j) = coef_ice_heat(i,j) +                     &
     &          b2d(i,j)

          IF (ai(i,j,linew) .gt. min_a(ng)) THEN

! downward conductivity term, assuming the ocean at the freezing point
              rhs_ice_heat(i,j) = rhs_ice_heat(i,j) +                   &
     &              b2d(i,j)*ti(i,j,linew)
              tis(i,j) = rhs_ice_heat(i,j)/coef_ice_heat(i,j)
            if (tis(i,j) .lt. -45._r8) tis(i,j) = -45._r8
          ELSE
            tis(i,j) = temp_top(i,j)
          END IF
        END DO
      END DO

      DO j = Jstr,Jend
        DO i = Istr,Iend
          sice = MIN(sice_ref,salt_top(i,j))
!**** calculate interior ice temp and heat fluxes
!       new temperature in ice
          IF (ai(i,j,linew) .gt. min_a(ng)) THEN
            cot = -frln*sice*hfus/(ti(i,j,linew)-eps)**2 + cpi
            ti(i,j,linew) = ti(i,j,linew) + dtice(ng)*(                 &
     &      2._r8*alph(i,j)/(rhoice(ng)*ice_thick(i,j)**2*cot)          &
     &         *(t0mk(i,j) + (tis(i,j) - (2._r8+coa(i,j))*ti(i,j,linew))&
     &                                        /(1._r8+coa(i,j))))
            ti(i,j,linew) = max(ti(i,j,linew),-35._r8)
            ti(i,j,linew) = min(ti(i,j,linew),-eps)
          ELSE
            ti(i,j,linew) = temp_top(i,j)
          END IF
        END DO
      END DO

      DO j = Jstr,Jend
        DO i = Istr,Iend
          IF (ai(i,j,linew) .gt. min_a(ng)) THEN
            t2(i,j) = (tis(i,j)+coa(i,j)*ti(i,j,linew))/(1._r8+coa(i,j))
            hicehinv = 1._r8/(0.5_r8*ice_thick(i,j))
            qi2(i,j) = alph(i,j)*(ti(i,j,linew)-t2(i,j))*hicehinv
            qio(i,j) = alph(i,j)*(t0mk(i,j)-ti(i,j,linew))*hicehinv
          END IF

!  Compute net heat flux from ice to atmosphere - Mellor and Kantha (7)
          qai(i,j) = qai_n(i,j)
        END DO
      END DO

!**** for open water ice fluxes set to zero
      DO j = Jstr,Jend
        DO i = Istr,Iend
          IF (ai(i,j,linew) .le. min_a(ng)) THEN
# ifdef MASKING
#  ifdef WET_DRY
            tis(i,j) = t0mk(i,j)*rmask(i,j)*rmask_wet(i,j)
            t2(i,j) = t0mk(i,j)*rmask(i,j)*rmask_wet(i,j)
            ti(i,j,linew) = -2.0_r8*rmask(i,j)*rmask_wet(i,j)
#  else
            tis(i,j) = t0mk(i,j)*rmask(i,j)
            t2(i,j) = t0mk(i,j)*rmask(i,j)
            ti(i,j,linew) = -2.0_r8*rmask(i,j)
#  endif
# elif defined WET_DRY
            tis(i,j) = t0mk(i,j)*rmask_wet(i,j)
            t2(i,j) = t0mk(i,j)*rmask_wet(i,j)
            ti(i,j,linew) = -2.0_r8*rmask_wet(i,j)
# else
            tis(i,j) = t0mk(i,j)
            t2(i,j) = t0mk(i,j)
            ti(i,j,linew) = -2.0_r8
# endif
# ifdef ICESHELF
            IF (zice(i,j).ne.0.0_r8) THEN
              tis(i,j) = 0.0_r8
              t2(i,j) = 0.0_r8
              ti(i,j,linew) = 0.0_r8
            END IF
# endif
            qi2(i,j) = 0._r8
            qai(i,j) = 0._r8
            qio(i,j) = 0._r8
            hsn(i,j,linew) = 0._r8
            sfwat(i,j,linew) = 0._r8
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
          sice = MIN(sice_ref,salt_top(i,j))
          tfrz = frln*sice
          wsm(i,j) = 0._r8
          wai(i,j) = 0._r8
          wro(i,j) = 0._r8
          IF (ai(i,j,linew) .gt. min_a(ng)) THEN
! Melt ice or freeze surface water in the fall if there is no snow
            IF (hsn(i,j,linew) .le. 0.0_r8) THEN
              IF (tis(i,j) .gt. tfrz .or. sfwat(i,j,linew) .gt. 0._r8)   &
     &                                                            THEN
!   ice warmer than freezing point
                tis(i,j) = tfrz
                hfus1(i,j) = hfus*(1._r8-brnfr(i,j))+tis(i,j)*cpw        &
     &         -((1._r8-brnfr(i,j))*cpi+brnfr(i,j)*cpw)*ti(i,j,linew)
                qai(i,j) = qai_n(i,j)
                qi2(i,j) = b2d(i,j)*(ti(i,j,linew)-tis(i,j))
!    compute ice production rate (negative here) from atmosphere-ice exchange
                wai(i,j) = -(qai(i,j) -qi2(i,j)) /(hfus1(i,j)*rhosw)
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
!          snow melting
                wsm(i,j) = max(0.0_r8,-(qai(i,j)-qi2(i,j))/             &
     &                    (rhosnow_wet(ng)*hfus)) + ws(i,j)
              END IF

              IF (tis(i,j) .lt. 0.0_r8 .and.                            &
     &            sfwat(i,j,linew) .gt. 0.0_r8) THEN
!          colder than the freezing point
                tis(i,j) = 0._r8
                qai(i,j) = qai_n(i,j)
                qi2(i,j) = b2d(i,j)*(ti(i,j,linew)-tis(i,j))
                wai(i,j) = -(qai(i,j)-qi2(i,j))/(hfus*rhosw)
              END IF
            END IF
!
!***** compute snow thickness and surface melt water
!       hsn - snow thickness
!       sfwat - surface meltwater
            hsn(i,j,linew) = hsn(i,j,linew)+(ai(i,j,linew)              &
     &                             *(-wsm(i,j)+ws(i,j)))*dtice(ng)
            hsn(i,j,linew) = max(0.0_r8,hsn(i,j,linew))
            sfwat(i,j,linew) = sfwat(i,j,linew)+                        &
     &                        (wai(i,j)+wsm(i,j))*dtice(ng)
            sfwat(i,j,linew) = max(0.0_r8,sfwat(i,j,linew))
!
            IF (sfwat(i,j,linew) .gt. sfwatmx) THEN
              wro(i,j) = (sfwat(i,j,linew)-sfwatmx)/dtice(ng)
              sfwat(i,j,linew) = sfwatmx
            END IF
          END IF
        END DO
      END DO

      DO j = Jstr,Jend
        DO i = Istr,Iend
          thic = ice_thick(i,j)

          z0 = max(z0ii*thic,0.01_r8)
          z0 = min(z0,0.1_r8)
!
!     *** Yaglom and Kader formulation for z0t and z0s
!
          zdz0 = dztop(i,j)/z0   !WPB
          if (zdz0 .lt. 3._r8) zdz0 = 3._r8

          rno = utau(i,j)*0.09_r8/nu
          termt = ykf*sqrt(rno)*prt**0.666667_r8
          terms = ykf*sqrt(rno)*prs**0.666667_r8
          cht(i,j) = utau(i,j)/(tpr*(log(zdz0)/kappa+termt))
          chs(i,j) = utau(i,j)/(tpr*(log(zdz0)/kappa+terms))
        END DO
      END DO

      DO j = Jstr,Jend
        DO i = Istr,Iend
          sice = MIN(sice_ref,salt_top(i,j))
          tfz = frln*salt_top(i,j)
          wao(i,j) = 0._r8
          wio(i,j) = 0._r8

          xwai = max(0._r8,wai(i,j))
          hfus1(i,j) = hfus*(1.0_r8-brnfr(i,j))+t0mk(i,j)*cpw           &
     &         -((1.0_r8-brnfr(i,j))*cpi+brnfr(i,j)*cpw)*ti(i,j,linew)
          IF (temp_top(i,j) .le. tfz)                                   &
     &             wao(i,j) = qao_n(i,j)/(hfus1(i,j)*rhosw)
          IF (ai(i,j,linew) .le. min_a(ng)) THEN
            s0mk(i,j) = salt_top(i,j)
            t0mk(i,j) = temp_top(i,j)
            wai(i,j) = 0._r8
            xtot = (1._r8-ai(i,j,linew))*wao(i,j)
          ELSE

! MK89 version
            wio(i,j) = (qio(i,j)/rhosw +                                &
     &               cpw*cht(i,j)*(t0mk(i,j)-temp_top(i,j)))/hfus1(i,j)

            xtot = ai(i,j,linew)*wio(i,j)                               &
     &              +(1._r8-ai(i,j,linew))*wao(i,j)

            s0mk(i,j) = (chs(i,j)*salt_top(i,j)+(xwai-wio(i,j))*sice)   &
     &                /(chs(i,j)+xwai+wro(i,j)-wio(i,j))
            s0mk(i,j) = max(s0mk(i,j),0._r8)
            s0mk(i,j) = min(s0mk(i,j),40._r8)

            t0mk(i,j) = frln*s0mk(i,j)

          END IF

          w0(i,j) = xtot-ai(i,j,linew)*wai(i,j)
          IF(ai(i,j,linew).LE.min_a(ng)) THEN
             stflx(i,j,itemp) = qao_n(i,j)
          ELSE
             stflx(i,j,itemp) = (1.0_r8-ai(i,j,linew))*qao_n(i,j)       &
     &                   +ai(i,j,linew)*qio(i,j)                        &
     &                   -xtot*hfus1(i,j)
          END IF

! Change stflx(i,j,itemp) back to ROMS convention
          stflx(i,j,itemp) = -stflx(i,j,itemp) * rhocpr

# ifdef MASKING
          stflx(i,j,itemp) = stflx(i,j,itemp)*rmask(i,j)
# endif
# ifdef WET_DRY
          stflx(i,j,itemp) = stflx(i,j,itemp)*rmask_wet(i,j)
# endif
# ifdef ICESHELF
          IF(zice(i,j).ne.0.0_r8) THEN
              stflx(i,j,itemp) = 0.0_r8
          END IF
# endif
          stflx(i,j,isalt) = stflx(i,j,isalt)                           &
     &        - (xtot-ai(i,j,linew)*xwai)*(sice-s0mk(i,j))              &
     &        - ai(i,j,linew)*wro(i,j)*s0mk(i,j)

! Test for case of rainfall on snow/ice and assume 100% drainage
#ifndef NCEP_FLUXES
          IF (rain(i,j).gt.0._r8.AND.snow_n(i,j).EQ.0._r8) THEN
            stflx(i,j,isalt) = stflx(i,j,isalt) -                       &
     &                         ai(i,j,linew)*rain(i,j)*0.001_r8
          END IF
#endif
!  io_mflux is ice production rate (+ve for growth)
          io_mflux(i,j) = xtot -ai(i,j,linew)*xwai -                    &
     &                          ai(i,j,linew)*wro(i,j) + wfr(i,j)
# ifdef MASKING
          stflx(i,j,isalt) = stflx(i,j,isalt)*rmask(i,j)
          io_mflux(i,j) = io_mflux(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          stflx(i,j,isalt) = stflx(i,j,isalt)*rmask_wet(i,j)
          io_mflux(i,j) = io_mflux(i,j)*rmask_wet(i,j)
# endif
# ifdef ICESHELF
         IF (zice(i,j).ne.0.0_r8) THEN
           stflx(i,j,isalt) = 0.0_r8
           io_mflux(i,j) = 0.0_r8
         END IF
# endif
        END DO
      END DO

!********************************

      DO j = Jstr,Jend
        DO i = Istr,Iend
          phi = 4._r8
          if (wao(i,j) .lt. 0.0_r8 ) phi = 0.5_r8
          xmelt = min((wio(i,j)-wai(i,j)),0.0_r8)
          hi(i,j,linew) = hi(i,j,linew)+dtice(ng)                       &
     &             *(ai(i,j,linew)                                      &
     &             *(wio(i,j)-wai(i,j))                                 &
     &        +(1.0_r8-ai(i,j,linew))*wao(i,j) + wfr(i,j))

          ai_tmp = ai(i,j,linew)
          ai(i,j,linew) = ai(i,j,linew) +                               &
     &             dtice(ng)*(1.0_r8-ai(i,j,linew))                     &
     &                      *(phi*wao(i,j)+wfr(i,j))
          ai(i,j,linew) = min(ai(i,j,linew),max_a(ng))

! adjust snow volume when ice melting out from under it
          IF (ai(i,j,linew) .lt. ai_tmp)                                &
     &        hsn(i,j,linew) =                                          &
     &           hsn(i,j,linew)*ai(i,j,linew)/max(ai_tmp,eps)
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
#ifdef MASKING
          ai(i,j,linew) = ai(i,j,linew)*rmask(i,j)
          hi(i,j,linew) = hi(i,j,linew)*rmask(i,j)
#endif
#ifdef WET_DRY
          ai(i,j,linew) = ai(i,j,linew)*rmask_wet(i,j)
          hi(i,j,linew) = hi(i,j,linew)*rmask_wet(i,j)
#endif
#ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) THEN
            ai(i,j,linew) = 0.0_r8
            hi(i,j,linew) = 0.0_r8
          END IF
#endif

#undef DIAG_WPB
#ifdef DIAG_WPB
      IF (i.eq.156.and.j.eq.481) THEN
         write(*,*) tdays,wio(i,j),wai(i,j),wao(i,j),wfr(i,j),          &
     &              xmelt,ai(i,j,linew),tis(i,j),                       &
     &                                     sfwat(i,j,linew),            &
     &              temp_top(i,j),t0mk(i,j),stflx(i,j,itemp),           &
     &              salt_top(i,j),s0mk(i,j),stflx(i,j,isalt),           &
     &              qio(i,j), ti(i,j,linew), brnfr(i,j),                &
     &              t2(i,j)
      END IF
#endif

        ENDDO
      ENDDO

!********************************
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ai(i,j,linew) = MIN(ai(i,j,linew),max_a(ng))
          ai(i,j,linew) = MAX(ai(i,j,linew),min_a(ng))
          hi(i,j,linew) = MAX(hi(i,j,linew),min_h(ng))
          hsn(i,j,linew) = MAX(hsn(i,j,linew),0.0_r8)
          sfwat(i,j,linew) = MAX(sfwat(i,j,linew),0.0_r8)
          ti(i,j,linew) = MAX(ti(i,j,linew),-70.0_r8)
          if (hi(i,j,linew) .le. min_h(ng)) ai(i,j,linew) = min_a(ng)
          if (ai(i,j,linew) .le. min_a(ng)) hi(i,j,linew) = min_h(ng)
        ENDDO
      ENDDO

        CALL bc_r2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tis)
        CALL bc_r2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          coef_ice_heat)
        CALL bc_r2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhs_ice_heat)
        CALL bc_r2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,isalt))
        CALL bc_r2d_tile (ng, tile,                                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,itemp))

        CALL aibc_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj, liold, linew,       &
     &                          ui, vi, ai)
        CALL hibc_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj, liold, linew,       &
     &                          ui, vi, hi)
        CALL hsnbc_tile (ng, tile,                                      &
     &                          LBi, UBi, LBj, UBj, liold, linew,       &
     &                          ui, vi, hsn)
        CALL tibc_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj, liold, linew,       &
     &                          min_h(ng), ui, vi, hi, ti, enthalpi)
        CALL sfwatbc_tile (ng, tile,                                    &
     &                        LBi, UBi, LBj, UBj, liold, linew,         &
     &                        ui, vi, sfwat)
        CALL ageicebc_tile (ng, tile,                                   &
     &                          LBi, UBi, LBj, UBj, liold, linew,       &
     &                          min_h(ng), ui, vi, hi, ageice, hage)
#if defined EW_PERIODIC || defined NS_PERIODIC
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
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sfwat(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ageice(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hage(:,:,linew))

#endif
#ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      ai(:,:,linew), hi(:,:,linew),               &
     &                      hsn(:,:,linew), ti(:,:,linew))
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      enthalpi(:,:,linew), sfwat(:,:,linew),      &
     &                      ageice(:,:,linew), hage(:,:,linew))
#endif

      RETURN
      END SUBROUTINE ice_thermo_tile
