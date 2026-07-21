module module_mp_tempo_main
  !! main tempo microphysics code
  use module_mp_tempo_cfgs, only : ty_tempo_cfgs
  use module_mp_tempo_params, only : wp, sp, dp, &
    min_qv, roverrv, rdry, r1, r2, nt_c_max, t0, nrhg, rho_g, &
    meters3_to_liters, eps, aero_max, nwfa_default, nifa_default
  use module_mp_tempo_utils, only : get_nuc, get_constant_cloud_number, snow_moments, calc_rslf, calc_rsif
  use module_mp_tempo_diags, only : reflectivity_10cm, effective_radius, max_hail_diam, &
    freezing_rain
  use module_mp_tempo_aerosols, only : init_ice_friendly_aerosols, init_water_friendly_aerosols, &
    aerosol_collection_efficiency
  use module_mp_tempo_ml, only : tempo_ml_predict_cloud_number
  implicit none
  private

  public :: tempo_main, ty_tempo_main_diags

#ifdef FV3
  public :: cloud_check_and_update, ice_check_and_update, snow_check_and_update
#endif
  
#ifdef unit_testing
  public :: get_cloud_table_index, get_snow_table_index, &
    get_temperature_table_index, get_rain_table_index, &
    get_graupel_table_index, get_ice_table_index
#endif
 
  type :: ty_tempo_main_diags
    real(wp) :: rain_precip
    real(wp) :: cloud_precip
    real(wp) :: ice_liquid_equiv_precip
    real(wp) :: snow_liquid_equiv_precip
    real(wp) :: graupel_liquid_equiv_precip
    real(wp) :: frozen_fraction
    real(wp) :: frz_rain_precip
    real(wp), dimension(:), allocatable :: rain_med_vol_diam
    real(wp), dimension(:), allocatable :: graupel_med_vol_diam
    real(wp), dimension(:), allocatable :: refl10cm
    real(wp), dimension(:), allocatable :: re_cloud
    real(wp), dimension(:), allocatable :: re_ice
    real(wp), dimension(:), allocatable :: re_snow
    real(wp), dimension(:), allocatable :: max_hail_diameter
    real(wp), dimension(:), allocatable :: cloud_number_mixing_ratio
  end type

  type :: ty_tend
    real(dp), pointer, contiguous, dimension(:) :: &
      prr_wau, pnr_wau, pnc_wau, prr_rcw, pnc_rcw, pnr_rcr, & ! warm rain
      prs_scw, pnc_scw, png_scw, pbg_scw, prg_gcw, pnc_gcw, pbg_gcw, & ! riming
      pri_ihm, pni_ihm, prs_ihm, prg_ihm, prg_scw, & ! riming
      prr_rcs, pnr_rcs, prg_rcs, png_rcs, prs_rcs, pbg_rcs, & ! rain-snow
      prr_rcg, pnr_rcg, prg_rcg, png_rcg, pbg_rcg, & ! rain-graupel
      pri_inu, pni_inu, pri_iha, pni_iha, & ! ice nucleation
      pri_wfz, pni_wfz, & ! water freezing
      prg_rfz, png_rfz, pnr_rfz, pri_rfz, pni_rfz, pbg_rfz, & ! rain freezing
      prs_sde, pri_ide, pni_ide, prs_ide, prg_gde, png_gde, & ! depositional growth
      pni_iau, prs_iau, & ! ice-snow conversion
      prr_sml, prr_gml, pbg_sml, pbg_gml, pnr_sml, pnr_gml, & ! melting
      prr_rci, pnr_rci, pri_rci, pni_rci, prg_rci, png_rci, pbg_rci, & ! rain-ice
      pni_sci, prs_sci, & ! snow-ice
      prw_vcd, pnc_wcd, prv_rev, pnr_rev, & ! condensation/evaporation
      pna_rca, pna_sca, pna_gca, pnd_rcd, pnd_scd, pnd_gcd ! aerosol
  end type

  interface get_cloud_number
    module procedure tempo_ml_predict_cloud_number
    module procedure get_constant_cloud_number
  end interface

  contains

  subroutine tempo_main(tempo_cfgs, &
    qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, qb1d, ni1d, nr1d, nc1d, ng1d, &
    nwfa1d, nifa1d, t1d, p1d, w1d, dz1d, &
    qcfrac1d, qifrac1d, qc_bl1d, qcfrac_bl1d, &
    thten_bl1d, qvten_bl1d, qcten_bl1d, qiten_bl1d, &
    thten_lwrad1d, thten_swrad1d, &
    kts, kte, dt, ii, jj, tempo_main_diags)
    !! tempo main

    type(ty_tend) :: tend
    type(ty_tempo_main_diags), intent(out) :: tempo_main_diags

    type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
    integer, intent(in) :: kts, kte, ii, jj
    real(wp), intent(in) :: dt
    real(wp), dimension(kts:kte), intent(inout) :: t1d !! 1D temperature \([K]\)
    real(wp), dimension(kts:kte), intent(in) :: p1d !! 1D pressure \([Pa]\)
    real(wp), dimension(kts:kte), intent(inout) :: qv1d !! 1D water vapor mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte), intent(inout) :: qc1d !! 1D cloud water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte), intent(inout) :: qr1d !! 1D rain water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte), intent(inout) :: qi1d !! 1D cloud ice mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte), intent(inout) :: qs1d !! 1D snow mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte), intent(inout) :: qg1d !! 1D graupel mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte), intent(inout) :: ni1d !! 1D cloud ice number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(kts:kte), intent(inout) :: nr1d !! 1D rain water number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(:), intent(inout), optional :: nc1d !! 1D cloud water number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(:), intent(inout), optional :: nwfa1d !! 1D water-friendly aerosol number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(:), intent(inout), optional :: nifa1d !! 1D ice-friendly aerosol number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(:), intent(inout), optional :: qb1d !! 1D graupel volume mixing ratio \([m^{-3}\; kg^{-1}]\)
    real(wp), dimension(:), intent(inout), optional :: ng1d !! 1D graupel number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(kts:kte), intent(in) :: w1d !! 1D vertical velocity \(m\; s^{-1}]\)
    real(wp), dimension(kts:kte), intent(in) :: dz1d !! 1D vertical grid spacing \([m]\)

    ! additional optional arrays
    real(wp), dimension(:), intent(inout), optional :: qcfrac1d !! cloud fraction
    real(wp), dimension(:), intent(inout), optional :: qifrac1d !! cloud ice fraction
    real(wp), dimension(:), intent(in), optional :: qc_bl1d !! cloud water mixing ratio from boundary layer scheme
    real(wp), dimension(:), intent(in), optional :: qcfrac_bl1d !! cloud fraction from boundary layer scheme
    real(wp), dimension(:), intent(in), optional :: thten_bl1d !! potential temperature tendency from boundary layer scheme
    real(wp), dimension(:), intent(in), optional :: qvten_bl1d !! water vapor mixing ratio tendency from boundary layer scheme
    real(wp), dimension(:), intent(in), optional :: qcten_bl1d !! cloud water mixing ratio tendency from boundary layer scheme
    real(wp), dimension(:), intent(in), optional :: qiten_bl1d !! cloud ice mixing ratio from boundary layer scheme
    real(wp), dimension(:), intent(in), optional :: thten_lwrad1d !! potential temperature tendency from longwave radiation scheme
    real(wp), dimension(:), intent(in), optional :: thten_swrad1d !! potential temperature tendency from shortwave radiation scheme
  
    real(wp), dimension(kts:kte) :: tten, qvten, qcten, qiten, qrten, qsten, &
      qgten, qbten, niten, nrten, ncten, ngten, nwfaten, nifaten !! tendencies

    logical, dimension(kts:kte) :: l_qc, l_qi, l_qr, l_qs, l_qg !! hydrometeor existence logicals
    integer, dimension(kts:kte) :: idx_bg !! graupel density index

    ! thermodynamic variables
    real(wp), dimension(kts:kte) :: temp, pres, qv !! thermodynamic variables
    real(wp), dimension(kts:kte) :: rho, rhof, rhof2 !! thermodynamic variables
    real(wp), dimension(kts:kte) :: qvs, qvsi, delqvs !! thermodynamic variables
    real(wp), dimension(kts:kte) :: satw, sati, ssatw, ssati !! thermodynamic variables
    real(wp), dimension(kts:kte) :: diffu, visco, vsc2, tcond, lvap, ocp, lvt2 !! thermodynamic variables
 
    real(wp), dimension(kts:kte) :: rc, ri, rr, rs, rg, rb !! local microphysical variables
    real(wp), dimension(kts:kte) :: ni, nr, nc, ng, nwfa, nifa !! local microphysics variables

    real(dp), dimension(kts:kte) :: ilamc, ilami, ilamr, ilamg !! inverse lambda
    real(wp), dimension(kts:kte) :: mvd_r, mvd_c, mvd_g !! median volume diameter
    real(dp), dimension(kts:kte) :: smob, smo2, smo1, smo0, smoc, smoe, smof, smog, ns, smoz !! snow moments
    
    real(wp), dimension(kts:kte) :: xrx, xnx !! temporary arrays
    real(wp), dimension(:), allocatable :: xncx, xngx, xqbx, ncsave !! temporary arrays

    real(wp), dimension(kts:kte+1) :: vtrr, vtnr, vtrs, vtri, vtni, vtrg, vtng, vtrc, vtnc !! fallspeeds
    real(wp), dimension(kts:kte) :: vtboost !! snow fallspeed boost factor
    integer :: substeps_sedi, ktop_sedi, n !! sedimentation substepping variables
    real(wp) :: semi_sedi_factor !! semi-lagrangian sedimentation factor

    real(dp), target, dimension(kts:kte, 74) :: tend_work !! array to store tendencies
    
    ! local variables
    real(wp) :: tempc, tc0, odt
    logical :: do_micro, supersaturated
    logical, save :: first_call_main = .true.
    integer :: k, nz

    ! --------------------------------------------------------------------------------------------
    do_micro = .false.
    supersaturated = .false.
    odt = 1._wp / dt

    nz = size(qv1d)

    ! map pointers to the contiguous stack workspace
    ! warm rain
    tend%prr_wau => tend_work(:, 1)
    tend%pnr_wau => tend_work(:, 2)
    tend%pnc_wau => tend_work(:, 3)
    tend%prr_rcw => tend_work(:, 4)
    tend%pnc_rcw => tend_work(:, 5)
    tend%pnr_rcr => tend_work(:, 6)

    ! riming
    tend%prs_scw => tend_work(:, 7)
    tend%pnc_scw => tend_work(:, 8)
    tend%png_scw => tend_work(:, 9)
    tend%pbg_scw => tend_work(:, 10)
    tend%prg_gcw => tend_work(:, 11)
    tend%pnc_gcw => tend_work(:, 12)
    tend%pbg_gcw => tend_work(:, 13)
    tend%pri_ihm => tend_work(:, 14)
    tend%pni_ihm => tend_work(:, 15)
    tend%prs_ihm => tend_work(:, 16)
    tend%prg_ihm => tend_work(:, 17)
    tend%prg_scw => tend_work(:, 18)

    ! rain-snow
    tend%prr_rcs => tend_work(:, 19)
    tend%pnr_rcs => tend_work(:, 20)
    tend%prg_rcs => tend_work(:, 21)
    tend%png_rcs => tend_work(:, 22)
    tend%prs_rcs => tend_work(:, 23)
    tend%pbg_rcs => tend_work(:, 24)

    ! rain-graupel
    tend%prr_rcg => tend_work(:, 25)
    tend%pnr_rcg => tend_work(:, 26)
    tend%prg_rcg => tend_work(:, 27)
    tend%png_rcg => tend_work(:, 28)
    tend%pbg_rcg => tend_work(:, 29)

    ! ice nucleation
    tend%pri_inu => tend_work(:, 30)
    tend%pni_inu => tend_work(:, 31)
    tend%pri_iha => tend_work(:, 32)
    tend%pni_iha => tend_work(:, 33)

    ! water freezing
    tend%pri_wfz => tend_work(:, 34)
    tend%pni_wfz => tend_work(:, 35)

    ! rain freezing
    tend%prg_rfz => tend_work(:, 36)
    tend%png_rfz => tend_work(:, 37)
    tend%pnr_rfz => tend_work(:, 38)
    tend%pri_rfz => tend_work(:, 39)
    tend%pni_rfz => tend_work(:, 40)
    tend%pbg_rfz => tend_work(:, 41)

    ! depositional growth
    tend%prs_sde => tend_work(:, 42)
    tend%pri_ide => tend_work(:, 43)
    tend%pni_ide => tend_work(:, 44)
    tend%prs_ide => tend_work(:, 45)
    tend%prg_gde => tend_work(:, 46)
    tend%png_gde => tend_work(:, 47)

    ! ice-snow conversion
    tend%pni_iau => tend_work(:, 48)
    tend%prs_iau => tend_work(:, 49)

    ! melting
    tend%prr_sml => tend_work(:, 50)
    tend%prr_gml => tend_work(:, 51)
    tend%pbg_sml => tend_work(:, 52)
    tend%pbg_gml => tend_work(:, 53)
    tend%pnr_sml => tend_work(:, 54)
    tend%pnr_gml => tend_work(:, 55)

    ! rain-ice
    tend%prr_rci => tend_work(:, 56)
    tend%pnr_rci => tend_work(:, 57)
    tend%pri_rci => tend_work(:, 58)
    tend%pni_rci => tend_work(:, 59)
    tend%prg_rci => tend_work(:, 60)
    tend%png_rci => tend_work(:, 61)
    tend%pbg_rci => tend_work(:, 62)

    ! snow-ice
    tend%pni_sci => tend_work(:, 63)
    tend%prs_sci => tend_work(:, 64)

    ! condensation/evaporation
    tend%prw_vcd => tend_work(:, 65)
    tend%pnc_wcd => tend_work(:, 66)
    tend%prv_rev => tend_work(:, 67)
    tend%pnr_rev => tend_work(:, 68)

    ! aerosol
    tend%pna_rca => tend_work(:, 69)
    tend%pna_sca => tend_work(:, 70)
    tend%pna_gca => tend_work(:, 71)
    tend%pnd_rcd => tend_work(:, 72)
    tend%pnd_scd => tend_work(:, 73)
    tend%pnd_gcd => tend_work(:, 74)

    ! zero out all mp tendencies
    tend_work = 0._dp
    
    ! zero tendencies
    do k = 1, nz
      tten(k) = 0._wp
      qvten(k) = 0._wp
      qcten(k) = 0._wp
      qiten(k) = 0._wp
      qrten(k) = 0._wp
      qsten(k) = 0._wp
      qgten(k) = 0._wp
      ngten(k) = 0._wp
      qbten(k) = 0._wp
      niten(k) = 0._wp
      nrten(k) = 0._wp
      ncten(k) = 0._wp
      nwfaten(k) = 0._wp
      nifaten(k) = 0._wp
      smo0(k) = 0._dp
      smo1(k) = 0._dp
      smo2(k) = 0._dp
      smob(k) = 0._dp
      smoc(k) = 0._dp
      smoe(k) = 0._dp
      smof(k) = 0._dp
      smog(k) = 0._dp
      smoz(k) = 0._dp
      ns(k) = 0._dp
      vtboost(k) = 1._wp
    enddo

    ! fallspeeds and sedimentation
    do k = 1, nz+1
      vtrr(k) = 0._wp
      vtnr(k) = 0._wp
      vtrs(k) = 0._wp
      vtri(k) = 0._wp
      vtni(k) = 0._wp
      vtrc(k) = 0._wp
      vtnc(k) = 0._wp
      vtrg(k) = 0._wp
      vtng(k) = 0._wp
    enddo

    tempo_main_diags%rain_precip = 0._wp
    tempo_main_diags%cloud_precip = 0._wp
    tempo_main_diags%ice_liquid_equiv_precip = 0._wp
    tempo_main_diags%snow_liquid_equiv_precip = 0._wp
    tempo_main_diags%graupel_liquid_equiv_precip = 0._wp
    tempo_main_diags%frz_rain_precip = 0._wp
  
    ! initialization -----------------------------------------------------------------------------
    do k = 1, nz
      temp(k) = t1d(k)
      qv(k) = max(min_qv, qv1d(k))
      pres(k) = p1d(k)
      rho(k) = roverrv*pres(k)/(rdry*temp(k)*(qv(k)+roverrv))
    enddo

    if (present(nwfa1d)) then
      if (first_call_main) then
        if (sum(nwfa1d) < eps) call init_water_friendly_aerosols(dz1d, nwfa)
      endif
    else
      call init_water_friendly_aerosols(dz1d, nwfa)
    endif
    
    if (present(nifa1d)) then
      if (first_call_main) then
        if (sum(nifa1d) < eps) call init_ice_friendly_aerosols(dz1d, nifa)
      endif
    else
      call init_ice_friendly_aerosols(dz1d, nifa)    
    endif

    call aerosol_check_and_update(rho=rho, nwfa1d=nwfa1d, nifa1d=nifa1d, &
      nwfa=nwfa, nifa=nifa, nwfaten=nwfaten, nifaten=nifaten, dt=dt)

    call rain_check_and_update(rho, l_qr, qr1d, nr1d, rr, nr, qrten, nrten, ilamr, mvd_r, dt, odt)
  
    call ice_check_and_update(rho, l_qi, qi1d, ni1d, ri, ni, qiten, niten, ilami, dt, odt)
  
    call snow_check_and_update(rho, l_qs, qs1d, rs, qsten, dt, odt)
    ! snow moments
    do k = 1, nz
      if (l_qs(k)) then
        tc0 = min(-0.1, temp(k)-t0)
        call snow_moments(rs=rs(k), tc=tc0, &
          smob=smob(k), smoc=smoc(k), ns=ns(k), &
          smo0=smo0(k), smo1=smo1(k), smo2=smo2(k), &
          smoe=smoe(k), smof=smof(k), smog=smog(k))
      endif 
    enddo

    ! set one-moment cloud number concentration
    if (.not. present(nc1d)) then
      if (tempo_cfgs%ml_for_nc_flag) then
        xrx = qc1d
        where(xrx <= 1.e-12_wp) xrx = 0._wp
        ! ml prediction
        call get_cloud_number(xrx, qr1d, qi1d, qs1d, pres, temp, w1d, xnx)
        nc = xnx * rho
      else
        ! single modment constant value
        call get_cloud_number(nc=nc)
      endif
      allocate(ncsave(nz), source=nc)
    endif

    call cloud_check_and_update(rho=rho, l_qc=l_qc, qc1d=qc1d, nc1d=nc1d, &
      ncsave=ncsave, rc=rc, nc=nc, qcten=qcten, ncten=ncten, ilamc=ilamc, mvd_c=mvd_c, &
      dt=dt, odt=odt)

    ! init ng and qb
    if (first_call_main) then
      if (present(ng1d) .and. present(qb1d)) then
        if (sum(qg1d) > r1 .and. sum(ng1d) < eps .and. sum(qb1d) < eps) then
          call graupel_init(rho, qg1d, ng1d, qb1d)
        endif 
      endif
    endif 
    call graupel_check_and_update(rho=rho, l_qg=l_qg, qg1d=qg1d, ng1d=ng1d, &
      qb1d=qb1d, rg=rg, ng=ng, rb=rb, idx=idx_bg, qgten=qgten, ngten=ngten, &
      qbten=qbten, ilamg=ilamg, mvd_g=mvd_g, dt=dt, odt=odt)

    ! re-zero tendencies after initial check zero tendencies
    do k = 1, nz
      qcten(k) = 0._wp
      qiten(k) = 0._wp
      qrten(k) = 0._wp
      qsten(k) = 0._wp
      qgten(k) = 0._wp
      ngten(k) = 0._wp
      qbten(k) = 0._wp
      niten(k) = 0._wp
      nrten(k) = 0._wp
      ncten(k) = 0._wp
    enddo

    call thermo_vars(qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
      satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, &
      supersaturated)

    if (first_call_main) first_call_main = .false.

    ! check for hydrometeors or supersaturation --------------------------------------------------
    do_micro = any(l_qc) .or. any(l_qr) .or. any(l_qi) .or. any(l_qs) .or. any(l_qg)
    if (.not. do_micro .and. .not. supersaturated) return

    ! main microphysical processes ---------------------------------------------------------------
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call warm_rain(rhof, l_qc, rc, nc, ilamc, mvd_c, l_qr, rr, nr, mvd_r, tend, odt)
    endif 
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call rain_snow_rain_graupel(temp, l_qr, rr, nr, ilamr, l_qs, rs, &
        l_qg, rg, ng, ilamg, idx_bg, tend, odt)
    endif 
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call ice_nucleation(temp, rho, w1d, qv, qvsi, ssati, ssatw, &
        nifa, nwfa, ni, smo0, rc, nc, rr, nr, ilamr, tend, dt, odt)
    endif 
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call ice_processes(rhof, rhof2, rho, w1d, temp, qv, qvsi, tcond, diffu, &
      vsc2, ssati, l_qi, ri, ni, ilami, l_qs, rs, smoe, smof, smo1, rr, nr, &
      ilamr, mvd_r, l_qg, rg, ng, ilamg, idx_bg, tend, odt)
    endif
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call riming(temp, rhof, visco, l_qc, rc, nc, ilamc, mvd_c, l_qs, rs, &
        smo0, smob, smoc, smoe, vtboost, l_qg, rg, ng, ilamg, idx_bg, tend, odt)
    endif 
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call melting(rhof2, rho, temp, qvsi, tcond, diffu, vsc2, ssati, delqvs, &
      l_qs, rs, smof, smo0, smo1, l_qg, rg, ng, ilamg, idx_bg, tend, dt, odt)
    endif
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call aerosol_scavenging(temp, rho, rhof, visco, nwfa, nifa, l_qr, nr, ilamr, &
      mvd_r, l_qs, rs, smob, smoc, smoe, l_qg, rg, ng, ilamg, idx_bg, tend, odt)
    endif 

    ! check and sum tendencies -------------------------------------------------------------------
    call check_over_depletion(rho, temp, qvsi, qv, l_qc, rc, l_qi, ri, &
      l_qr, rr, l_qs, rs, l_qg, rg, tend, odt)

    call sum_tendencies(rho, temp, idx_bg, lvap, ocp, tend, tten, qvten, qcten, &
      ncten, qiten, niten, qsten, qrten, nrten, qgten, ngten, qbten)

    ! update after tendencies applied ------------------------------------------------------------
    do k = 1, nz
      temp(k) = t1d(k) + tten(k)*dt
      tempc = temp(k) - t0
      qv(k) = max(min_qv, qv1d(k) + qvten(k)*dt)
      rho(k) = roverrv*pres(k)/(rdry*temp(k)*(qv(k)+roverrv))
      nwfaten(k) = nwfaten(k) - (tend%pna_rca(k) + tend%pna_sca(k) + tend%pna_gca(k) + &
        tend%pni_iha(k)) / rho(k)
      nifaten(k) = nifaten(k) - (tend%pnd_rcd(k) + tend%pnd_scd(k) + tend%pnd_gcd(k)) / rho(k)
    enddo 
    
    ! only updates nwfa, nifa
    call aerosol_check_and_update(rho=rho, nwfa1d=nwfa1d, nifa1d=nifa1d, &
      nwfa=nwfa, nifa=nifa, nwfaten=nwfaten, nifaten=nifaten, dt=dt)

    ! send temporary arrays to avoid updates to 1d variables at this point
    xrx = qc1d
    if (present(nc1d)) then
      if (.not. allocated(xncx)) allocate(xncx(nz), source=0._wp)
      xncx = nc1d
    endif 
    call cloud_check_and_update(rho=rho, l_qc=l_qc, qc1d=xrx, nc1d=xncx, &
      ncsave=ncsave, rc=rc, nc=nc, qcten=qcten, ncten=ncten, ilamc=ilamc, mvd_c=mvd_c, &
      dt=dt, odt=odt)

    xrx = qr1d
    xnx = nr1d
    call rain_check_and_update(rho, l_qr, xrx, xnx, rr, nr, qrten, nrten, ilamr, mvd_r, dt, odt)

    xrx = qi1d
    xnx = ni1d
    call ice_check_and_update(rho, l_qi, xrx, xnx, ri, ni, qiten, niten, ilami, dt, odt)

    xrx = qs1d
    call snow_check_and_update(rho, l_qs, xrx, rs, qsten, dt, odt)
    ! snow moments
    do k = 1, nz
      smo0(k) = 0._dp
      smo1(k) = 0._dp
      smo2(k) = 0._dp
      smob(k) = 0._dp
      smoc(k) = 0._dp
      smoe(k) = 0._dp
      smof(k) = 0._dp
      smog(k) = 0._dp
      ns(k) = 0._dp
      if (l_qs(k)) then
        tc0 = min(-0.1, temp(k)-t0)
        call snow_moments(rs=rs(k), tc=tc0, &
          smob=smob(k), smoc=smoc(k), &
          smo2=smo2(k))
      endif 
    enddo

    xrx = qg1d
    if (present(ng1d) .and. present(qb1d)) then
      if (.not. allocated(xngx)) allocate(xngx(nz), source=0._wp)
      if (.not. allocated(xqbx)) allocate(xqbx(nz), source=0._wp)
      xngx = ng1d
      xqbx = qb1d
    endif 
    call graupel_check_and_update(rho=rho, l_qg=l_qg, qg1d=xrx, ng1d=xngx, &
      qb1d=xqbx, rg=rg, ng=ng, rb=rb, idx=idx_bg, qgten=qgten, ngten=ngten, &
      qbten=qbten, ilamg=ilamg, mvd_g=mvd_g, dt=dt, odt=odt)

    call thermo_vars(qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
      satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, &
      supersaturated)

    ! after update do cloud condensation / rain evaporation --------------------------------------
    ! cloud condensation
    if (.not. tempo_cfgs%turn_off_micro_flag .and. tempo_cfgs%cloud_condensation_flag) then
      call cloud_condensation(rho, temp, w1d, ssatw, lvap, tcond, diffu, lvt2, &
        nwfa, qv, qvs, l_qc, rc, nc, tend, dt, odt)

      do k = 1, nz
        qvten(k) = qvten(k) - tend%prw_vcd(k)
        qcten(k) = qcten(k) + tend%prw_vcd(k)
        ncten(k) = ncten(k) + tend%pnc_wcd(k)
        nwfaten(k) = nwfaten(k) - tend%pnc_wcd(k)
        tten(k) = tten(k) + lvap(k)*ocp(k)*tend%prw_vcd(k)
      enddo 

      xrx = qc1d
      if (present(nc1d)) then
        if (.not. allocated(xncx)) allocate(xncx(nz), source=0._wp)
        xncx = nc1d
      endif 
      call cloud_check_and_update(rho=rho, l_qc=l_qc, qc1d=xrx, nc1d=xncx, &
        ncsave=ncsave, rc=rc, nc=nc, qcten=qcten, ncten=ncten, ilamc=ilamc, mvd_c=mvd_c, &
        dt=dt, odt=odt)

      do k = 1, nz
        qv(k) = max(min_qv, qv1d(k) + qvten(k)*dt)
        temp(k) = t1d(k) + tten(k)*dt
        rho(k) = roverrv*pres(k)/(rdry*temp(k)*(qv(k)+roverrv))
      enddo 
      
      call thermo_vars(qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
      satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, &
      supersaturated)
    endif 

    ! rain evaporation
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call rain_evaporation(rho, temp, ssatw, lvap, tcond, diffu, vsc2, rhof2, &
        qv, qvs, l_qr, rr, nr, ilamr, tend, odt)

      do k = 1, nz
        qrten(k) = qrten(k) - tend%prv_rev(k)
        qvten(k) = qvten(k) + tend%prv_rev(k)
        nrten(k) = nrten(k) - tend%pnr_rev(k)
        nwfaten(k) = nwfaten(k) + tend%pnr_rev(k)
        tten(k) = tten(k) - lvap(k)*ocp(k)*tend%prv_rev(k)
      enddo

      xrx = qr1d
      xnx = nr1d
      call rain_check_and_update(rho, l_qr, xrx, xnx, rr, nr, qrten, nrten, ilamr, mvd_r, dt, odt)

      do k = 1, nz
        qv(k) = max(min_qv, qv1d(k) + qvten(k)*dt)
        temp(k) = t1d(k) + tten(k)*dt
        rho(k) = roverrv*pres(k)/(rdry*temp(k)*(qv(k)+roverrv))
      enddo 
        
      call thermo_vars(qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
        satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, &
        supersaturated)
    endif 

    ! sedimentation ------------------------------------------------------------------------------
    
    ! rain
    ktop_sedi = 1
    substeps_sedi = 1
    semi_sedi_factor = 10._wp
    if (any(l_qr)) then
      call rain_fallspeed(rhof=rhof, l_qr=l_qr, rr=rr, ilamr=ilamr, dz1d=dz1d, &
        vt=vtrr, vtn=vtnr, substeps_sedi=substeps_sedi, ktop_sedi=ktop_sedi, dt=dt)

      if (tempo_cfgs%semi_sedi_flag) then
        substeps_sedi = max(int(substeps_sedi/semi_sedi_factor) + 1, 1)
        do n = 1, substeps_sedi
          call semilagrangian_sedimentation(dz1d=dz1d, rho=rho, xr=rr, xten=qrten, &
            vt=vtrr, steps=substeps_sedi, limit=r1, precip=tempo_main_diags%rain_precip, dt=dt, odt=odt)
          call semilagrangian_sedimentation(dz1d=dz1d, rho=rho, xr=nr, xten=nrten, &
            vt=vtnr, steps=substeps_sedi, limit=r2, dt=dt, odt=odt)
          vtrr = 0._wp
          vtnr = 0._wp
          xrx = qr1d
          xnx = nr1d
          call rain_check_and_update(rho, l_qr, xrx, xnx, rr, nr, qrten, nrten, ilamr, mvd_r, dt, odt) 
          call rain_fallspeed(rhof=rhof, l_qr=l_qr, rr=rr, ilamr=ilamr, dz1d=dz1d, &
            vt=vtrr, vtn=vtnr, dt=dt)
        enddo 
      else
        do n = 1, substeps_sedi
          call sedimentation(xr=rr, vt=vtrr, dz1d=dz1d, rho=rho, xten=qrten, limit=r1, &
            steps=substeps_sedi, ktop_sedi=ktop_sedi, precip=tempo_main_diags%rain_precip, dt=dt)
          call sedimentation(xr=nr, vt=vtnr, dz1d=dz1d, rho=rho, xten=nrten, limit=r2, &
            steps=substeps_sedi, ktop_sedi=ktop_sedi, dt=dt)
          ! vtrr = 0._wp
          ! vtnr = 0._wp
          ! xrx = qr1d
          ! xnx = nr1d
          ! call rain_check_and_update(rho, l_qr, xrx, xnx, rr, nr, qrten, nrten, ilamr, mvd_r, dt, odt)
          ! call rain_fallspeed(rhof=rhof, l_qr=l_qr, rr=rr, ilamr=ilamr, dz1d=dz1d, &
          !   vt=vtrr, vtn=vtnr, dt=dt)
        enddo
      endif 
    endif

    ! graupel
    ktop_sedi = 1
    substeps_sedi = 1
    semi_sedi_factor = 10._wp
    if (any(l_qg)) then
      call graupel_fallspeed(rhof=rhof, rho=rho, visco=visco, &
        l_qg=l_qg, rg=rg, rb=rb, qb1d=qb1d, idx=idx_bg, ilamg=ilamg, dz1d=dz1d, &
        vt=vtrg, vtn=vtng, substeps_sedi=substeps_sedi, ktop_sedi=ktop_sedi, dt=dt)

      if (tempo_cfgs%semi_sedi_flag) then 
        substeps_sedi = max(int(substeps_sedi/semi_sedi_factor) + 1, 1)
        do n = 1, substeps_sedi
          call semilagrangian_sedimentation(dz1d=dz1d, rho=rho, xr=rg, xten=qgten, &
            vt=vtrg, steps=substeps_sedi, limit=r1, precip=tempo_main_diags%graupel_liquid_equiv_precip, &
            dt=dt, odt=odt)
           call semilagrangian_sedimentation(dz1d=dz1d, rho=rho, xr=ng, xten=ngten, &
            vt=vtng, steps=substeps_sedi, limit=r2, dt=dt, odt=odt)
          call semilagrangian_sedimentation(dz1d=dz1d, rho=rho, xr=rb, xten=qbten, &
            vt=vtrg, steps=substeps_sedi, limit=meters3_to_liters*r1/rho_g(nrhg), dt=dt, odt=odt)
          vtrg = 0._wp
          vtng = 0._wp
          xrx = qg1d
          if (present(ng1d) .and. present(qb1d)) then
            if (.not. allocated(xngx)) allocate(xngx(nz), source=0._wp)
            if (.not. allocated(xqbx)) allocate(xqbx(nz), source=0._wp)
            xngx = ng1d
            xqbx = qb1d
          endif 
          call graupel_check_and_update(rho=rho, l_qg=l_qg, qg1d=xrx, ng1d=xngx, &
            qb1d=xqbx, rg=rg, ng=ng, rb=rb, idx=idx_bg, qgten=qgten, ngten=ngten, &
            qbten=qbten, ilamg=ilamg, mvd_g=mvd_g, dt=dt, odt=odt)
          call graupel_fallspeed(rhof=rhof, rho=rho, visco=visco, &
            l_qg=l_qg, rg=rg, rb=rb, qb1d=qb1d, idx=idx_bg, ilamg=ilamg, dz1d=dz1d, &
            vt=vtrg, vtn=vtng, dt=dt)
        enddo 
      else
        do n = 1, substeps_sedi
          call sedimentation(xr=rg, vt=vtrg, dz1d=dz1d, rho=rho, xten=qgten, limit=r1, &
            steps=substeps_sedi, ktop_sedi=ktop_sedi, precip=tempo_main_diags%graupel_liquid_equiv_precip, dt=dt)
          call sedimentation(xr=ng, vt=vtng, dz1d=dz1d, rho=rho, xten=ngten, limit=r2, &
            steps=substeps_sedi, ktop_sedi=ktop_sedi, dt=dt)
          call sedimentation(xr=rb, vt=vtrg, dz1d=dz1d, rho=rho, xten=qbten, &
            limit=meters3_to_liters*r1/rho_g(nrhg), steps=substeps_sedi, ktop_sedi=ktop_sedi, dt=dt)
          ! vtrg = 0._wp
          ! vtng = 0._wp
          ! xrx = qg1d
          ! if (present(ng1d) .and. present(qb1d)) then
          !   if (.not. allocated(xngx)) allocate(xngx(nz), source=0._wp)
          !   if (.not. allocated(xqbx)) allocate(xqbx(nz), source=0._wp)
          !   xngx = ng1d
          !   xqbx = qb1d
          ! endif 
          ! call graupel_check_and_update(rho=rho, l_qg=l_qg, qg1d=xrx, ng1d=xngx, &
          !   qb1d=xqbx, rg=rg, ng=ng, rb=rb, idx=idx_bg, qgten=qgten, ngten=ngten, &
          !   qbten=qbten, ilamg=ilamg, mvd_g=mvd_g, dt=dt, odt=odt)
          ! call graupel_fallspeed(rhof=rhof, rho=rho, visco=visco, &
          !   l_qg=l_qg, rg=rg, rb=rb, qb1d=qb1d, idx=idx_bg, ilamg=ilamg, dz1d=dz1d, &
          !   vt=vtrg, vtn=vtng, dt=dt)
        enddo
      endif 
    endif

    ! snow
    ktop_sedi = 1
    substeps_sedi = 1
    if (any(l_qs)) then
      call snow_fallspeed(rhof=rhof, l_qs=l_qs, rs=rs, prr_sml=tend%prr_sml, smob=smob, smoc=smoc, &
        rr=rr, vtrr=vtrr, dz1d=dz1d, vt=vtrs, vtboost=vtboost, substeps_sedi=substeps_sedi, ktop_sedi=ktop_sedi, dt=dt)
      do n = 1, substeps_sedi
        call sedimentation(xr=rs, vt=vtrs, dz1d=dz1d, rho=rho, xten=qsten, limit=r1, &
        steps=substeps_sedi, ktop_sedi=ktop_sedi, precip=tempo_main_diags%snow_liquid_equiv_precip, dt=dt)
      enddo
    endif 

    ! ice
    ktop_sedi = 1
    substeps_sedi = 1
    if (any(l_qi)) then
      call ice_fallspeed(rhof, l_qi, ri, ilami, dz1d, vtri, vtni, &
        substeps_sedi, ktop_sedi, dt=dt)
      call sedimentation(xr=ri, vt=vtri, dz1d=dz1d, rho=rho, xten=qiten, limit=r1, &
        steps=substeps_sedi, ktop_sedi=ktop_sedi, precip=tempo_main_diags%ice_liquid_equiv_precip, dt=dt)
      call sedimentation(xr=ni, vt=vtni, dz1d=dz1d, rho=rho, xten=niten, limit=r2, &
        steps=substeps_sedi, ktop_sedi=ktop_sedi, dt=dt)
    endif 

    ! cloud
    ktop_sedi = 1
    substeps_sedi = 1
    if (any(l_qc)) then
      call cloud_fallspeed(rhof, w1d, l_qc, rc, nc, ilamc, dz1d, vtrc, vtnc, ktop_sedi)
      call sedimentation(xr=rc, vt=vtrc, dz1d=dz1d, rho=rho, xten=qcten, limit=r1, &
        steps=substeps_sedi, ktop_sedi=ktop_sedi, precip=tempo_main_diags%cloud_precip, dt=dt)
      call sedimentation(xr=nc, vt=vtnc, dz1d=dz1d, rho=rho, xten=ncten, limit=r2, &
        steps=substeps_sedi, ktop_sedi=ktop_sedi, dt=dt)
    endif 

    ! after sedimentation freeze all cloud water below hgfrz temperature
    ! and melt all cloud ice above freezing
    if (.not. tempo_cfgs%turn_off_micro_flag) then
      call freeze_cloud_melt_ice(temp=temp, rho=rho, ocp=ocp, lvap=lvap, &
        qi1d=qi1d, ni1d=ni1d, qiten=qiten, niten=niten, qc1d=qc1d, nc1d=nc1d, &
        qcten=qcten, ncten=ncten, tten=tten, ncsave=ncsave, dt=dt, odt=odt)
    endif 

    ! final update -------------------------------------------------------------------------------
    do k = 1, nz
      t1d(k)  = t1d(k) + tten(k)*dt
      qv1d(k) = max(min_qv, (qv1d(k) + qvten(k)*dt))
      rho(k) = roverrv*pres(k)/(rdry*temp(k)*(qv(k)+roverrv))
      if (present(nwfa1d)) then
        nwfa1d(k) = max(nwfa_default, min(aero_max, (nwfa1d(k)+nwfaten(k)*dt)))
      endif 
      if (present(nifa1d)) then
        nifa1d(k) = max(nifa_default, min(aero_max, (nifa1d(k)+nifaten(k)*dt)))
      endif 
    enddo
    
    call cloud_check_and_update(rho=rho, l_qc=l_qc, qc1d=qc1d, nc1d=nc1d, &
      ncsave=ncsave, rc=rc, nc=nc, qcten=qcten, ncten=ncten, ilamc=ilamc, mvd_c=mvd_c, &
      dt=dt, odt=odt)

    call rain_check_and_update(rho, l_qr, qr1d, nr1d, rr, nr, qrten, nrten, ilamr, mvd_r, dt, odt)
    
    call ice_check_and_update(rho, l_qi, qi1d, ni1d, ri, ni, qiten, niten, ilami, dt, odt)
    
    call snow_check_and_update(rho, l_qs, qs1d, rs, qsten, dt, odt)
    ! snow moments
    do k = 1, nz
      smo0(k) = 0._dp
      smo1(k) = 0._dp
      smo2(k) = 0._dp
      smob(k) = 0._dp
      smoc(k) = 0._dp
      smoe(k) = 0._dp
      smof(k) = 0._dp
      smog(k) = 0._dp
      smoz(k) = 0._dp
      ns(k) = 0._dp
      if (l_qs(k)) then
        tc0 = min(-0.1, temp(k)-t0)
        call snow_moments(rs=rs(k), tc=tc0, &
          smob=smob(k), smoc=smoc(k), &
          smo2=smo2(k), smoz=smoz(k))
      endif 
    enddo
  
    call graupel_check_and_update(rho=rho, l_qg=l_qg, qg1d=qg1d, ng1d=ng1d, &
      qb1d=qb1d, rg=rg, ng=ng, rb=rb, idx=idx_bg, qgten=qgten, ngten=ngten, &
      qbten=qbten, ilamg=ilamg, mvd_g=mvd_g, dt=dt, odt=odt)

    ! diagnostic output --------------------------------------------------------------------------
    ! frozen fraction
    tempo_main_diags%frozen_fraction = &
      (tempo_main_diags%ice_liquid_equiv_precip + tempo_main_diags%snow_liquid_equiv_precip + &
      tempo_main_diags%graupel_liquid_equiv_precip) / &
      (tempo_main_diags%ice_liquid_equiv_precip + tempo_main_diags%snow_liquid_equiv_precip + &
      tempo_main_diags%graupel_liquid_equiv_precip + tempo_main_diags%rain_precip + r1)

    ! freezing rain
    call freezing_rain(temp=temp(1), rain_precip=tempo_main_diags%rain_precip, &
      cloud_precip=tempo_main_diags%cloud_precip, &
      frz_rain=tempo_main_diags%frz_rain_precip)

    if (tempo_cfgs%cloud_number_mixing_ratio_flag) then
      allocate(tempo_main_diags%cloud_number_mixing_ratio(nz), source=nc*rho)
    endif 

    ! median volume diameter of rain and graupel
    if (tempo_cfgs%rain_med_vol_diam_flag) then
      allocate(tempo_main_diags%rain_med_vol_diam(nz), source=0._wp)
      tempo_main_diags%rain_med_vol_diam = mvd_r
    endif 
    if (tempo_cfgs%graupel_med_vol_diam_flag) then
      allocate(tempo_main_diags%graupel_med_vol_diam(nz), source=0._wp)
      tempo_main_diags%graupel_med_vol_diam = mvd_g
    endif 

    ! max hail diameter
    if (tempo_cfgs%max_hail_diameter_flag) then
      allocate(tempo_main_diags%max_hail_diameter(nz), source=0._wp)
      call max_hail_diam(rho, rg, ng, ilamg, idx_bg, &
        tempo_main_diags%max_hail_diameter)
    endif

    ! 10-cm reflectivity
    if (tempo_cfgs%refl10cm_flag) then
      allocate(tempo_main_diags%refl10cm(nz), source=-35._wp)
      call reflectivity_10cm(tempo_cfgs%refl10cm_from_melting_flag, &
        temp, l_qr, rr, nr, ilamr, &
        l_qs, rs, smoc, smob, smoz, l_qg, rg, ng, idx_bg, ilamg, &
        tempo_main_diags%refl10cm)
    endif 

    ! effective radii
    if ((tempo_cfgs%re_cloud_flag) .and. (tempo_cfgs%re_ice_flag) .and. (tempo_cfgs%re_snow_flag)) then
      allocate(tempo_main_diags%re_cloud(nz), source=0._wp)
      allocate(tempo_main_diags%re_ice(nz), source=0._wp)
      allocate(tempo_main_diags%re_snow(nz), source=0._wp)

      ! this next code block optionally adds pbl clouds to resolved clouds 
      ! for the effective radius calculation
      ! thus qc -> qc + qc_bl, a new value of nc is predicted with ML, 
      ! and then ilamc is updated along with rc and nc before re is calculated
      ! please output any cloud diagnostics before this calculation
      ! because rc, nc, ilamc, and mvd_c will include contributions from 
      ! resolved and explicit clouds and qcten and ncten are zeroed
      if (present(qc_bl1d) .and. present(qcfrac_bl1d)) then
        xrx = qc1d
        if (.not. allocated(xncx)) allocate(xncx(nz), source=0._wp)
        xncx = nc1d
        do k = 1, nz
          if ((xrx(k) <= r1) .and. & 
            (qc_bl1d(k) > 1.e-9_wp) .and. (qcfrac_bl1d(k) > 0._wp)) then
            xrx(k) = xrx(k) + qc_bl1d(k) / qcfrac_bl1d(k) ! use in-cloud PBL mass
          endif 
        enddo 
        where(xrx <= 1.e-12_wp) xrx = 0._wp
        ! ml prediction
        call get_cloud_number(xrx, qr1d, qi1d, qs1d, pres, temp, w1d, xncx)
      
        ! xrx and xncx have been updated to include pbl contribution -> update ilamc and nc 
        ! for effective radius calculation
        qcten = 0._wp
        ncten = 0._wp
        call cloud_check_and_update(rho=rho, l_qc=l_qc, qc1d=xrx, nc1d=xncx, &
          rc=rc, nc=nc, qcten=qcten, ncten=ncten, ilamc=ilamc, mvd_c=mvd_c, &
          dt=dt, odt=odt)
      endif
      call effective_radius(temp, l_qc, nc, ilamc, l_qi, ilami, l_qs, rs, &
        tempo_main_diags%re_cloud, tempo_main_diags%re_ice, tempo_main_diags%re_snow)
    endif
  end subroutine tempo_main


  subroutine aerosol_check_and_update(rho, nwfa1d, nifa1d, nwfa, nifa, nwfaten, nifaten, dt)
    !! sets aerosol number concentrations and checks bounds
    use module_mp_tempo_params, only : nwfa_default, aero_max, nifa_default

    real(wp), intent(in) :: dt
    real(wp), dimension(:), intent(in) :: rho
    real(wp), dimension(:), intent(inout) :: nwfa, nifa, nwfaten, nifaten
    real(wp), dimension(:), intent(inout), optional :: nwfa1d, nifa1d
    integer :: k, nz

    nz = size(rho)
    do k = 1, nz
      if (present(nwfa1d)) then
        nwfa(k) = (nwfa1d(k)+nwfaten(k)*dt)*rho(k)
      endif 
      if (present(nifa1d)) then
        nifa(k) = (nifa1d(k)+nifaten(k)*dt)*rho(k)
      endif
      nwfa(k) = max(nwfa_default*rho(k), min(aero_max*rho(k), nwfa(k)))
      nifa(k) = max(nifa_default*rho(k), min(aero_max*rho(k), nifa(k)))
    enddo 
  end subroutine aerosol_check_and_update


  subroutine cloud_check_and_update(rho, l_qc, qc1d, nc1d, ncsave, rc, nc, &
    qcten, ncten, ilamc, mvd_c, dt, odt)
    !! computes cloud water contents, ilamc, and mvd_c and checks bounds
    use module_mp_tempo_params, only : r1, nt_c_max, nt_c_min, nu_c_scale, &
      am_r, bm_r, cce, ccg, d0c, d0r, ocg1, ocg2, obmr, nt_c_l, d0r, nt_c_l

    real(wp), intent(in) :: dt, odt
    real(wp), dimension(:), intent(in) :: rho
    real(wp), dimension(:), intent(inout) :: qc1d, qcten, ncten, rc, nc
    real(wp), dimension(:), intent(out) :: mvd_c
    real(dp), dimension(:), intent(out) :: ilamc
    real(wp), dimension(:), intent(inout), optional :: nc1d, ncsave
    logical, dimension(:), intent(inout) :: l_qc
    integer :: k, nz, nu_c
    real(dp) :: lamc, xdc
    logical :: hit_limit

    nz = size(qc1d)
    do k = 1, nz
      hit_limit = .false.
      if (qc1d(k)+qcten(k)*dt > r1) then
        l_qc(k) = .true.
        ! update mass
        rc(k) = (qc1d(k)+qcten(k)*dt)*rho(k)
        qc1d(k) = qc1d(k)+qcten(k)*dt

        ! update number
        if (present(nc1d)) then
          nc(k) = max(nt_c_min, (nc1d(k)+ncten(k)*dt)*rho(k))
          
          ! number check
          if (nc(k) <= nt_c_min) then
            hit_limit = .true.
            nc(k) = nt_c_min
          endif 
          if (nc(k) > nt_c_max) then
            hit_limit = .true.
            nc(k) = nt_c_max
          endif 
        
          ! size check
          nu_c = get_nuc(nc(k))
          lamc = (nc(k)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
          xdc = (bm_r + nu_c + 1._dp) / lamc
          if (xdc < d0c) then
            lamc = cce(2,nu_c)/d0c
            hit_limit = .true.
          elseif (xdc > d0r*2._dp) then
            lamc = cce(2,nu_c)/(d0r*2._dp)
            hit_limit = .true.
          endif
          ! update number to be consistent with lamc
          nc(k) = ccg(1,nu_c)*ocg2(nu_c)*rc(k) / am_r*lamc**bm_r

          if (hit_limit) ncten(k) = (nc(k)/rho(k) - nc1d(k)) * odt
          nc1d(k) = max(nt_c_min/rho(k), &
            min(ccg(1,nu_c)*ocg2(nu_c)*qc1d(k)/am_r*lamc**bm_r, nt_c_max/rho(k)))
        else
          if (present(ncsave)) then
            nc(k) = ncsave(k)
          else
            nc(k) = nt_c_l
          endif
        endif
        nu_c = get_nuc(nc(k))
        lamc = (nc(k)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
        ilamc(k) = 1._dp / lamc
        mvd_c(k) = max(min((3.0_wp + nu_c + 0.672_wp) * ilamc(k), d0r), d0c)
      else
        l_qc(k) = .false.
        rc(k) = r1
        nc(k) = nt_c_min
        mvd_c(k) = d0c
        ilamc(k) = 0._dp
        qcten(k) = -qc1d(k) * odt
        qc1d(k) = 0.0_wp
        if (present(nc1d)) then
          ncten(k) = -nc1d(k) * odt
          nc1d(k) = 0.0_wp
        endif 
      endif
    enddo
  end subroutine cloud_check_and_update


  subroutine rain_check_and_update(rho, l_qr, qr1d, nr1d, rr, nr, &
      qrten, nrten, ilamr, mvd_r, dt, odt)
    !! computes rain water contents, ilamr, and mvd_r and checks bounds
    use module_mp_tempo_params, only : r1, r2, mu_r, crg, &
      am_r, bm_r, obmr, d0r, d0r_max, org2, org3, rho_g

    real(wp), intent(in) :: dt, odt
    real(wp), dimension(:), intent(in) :: rho
    real(wp), dimension(:), intent(inout) :: qr1d, nr1d, qrten, nrten, rr, nr
    real(wp), dimension(:), intent(out) :: mvd_r
    real(dp), dimension(:), intent(out) :: ilamr
    logical, dimension(:), intent(inout) :: l_qr
    integer :: k, nz
    real(dp) :: lamr
    logical :: hit_limit

    nz = size(qr1d)
    do k = 1, nz
      hit_limit = .false.
      if (qr1d(k)+qrten(k)*dt > r1) then
        l_qr(k) = .true.
        ! update mass
        rr(k) = (qr1d(k)+qrten(k)*dt)*rho(k)
        qr1d(k) = qr1d(k)+qrten(k)*dt

        ! update number
        nr(k) = max(r2, (nr1d(k)+nrten(k)*dt)*rho(k))

        ! number check
        if (nr(k) <= r2) then
          hit_limit = .true.
          mvd_r(k) = 1.0e-3_wp
          lamr = (3.0_dp + mu_r + 0.672_dp) / mvd_r(k)
          nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
        endif
        
        ! size check
        lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
        mvd_r(k) = (3.0_wp + mu_r + 0.672_wp) / lamr
        if (mvd_r(k) > d0r_max) then
          hit_limit = .true.
          mvd_r(k) = d0r_max
          lamr = (3.0_dp + mu_r + 0.672_dp) / mvd_r(k)
          nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
        elseif (mvd_r(k) < d0r*0.75_wp) then
          hit_limit = .true.
          mvd_r(k) = d0r*0.75_wp
          lamr = (3.0_dp + mu_r + 0.672_dp) / mvd_r(k)
          nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
        endif
        ! update number to be consistent with lamc
        if (hit_limit) nrten(k) = (nr(k)/rho(k) - nr1d(k))*odt
        nr1d(k) = crg(2)*org3*qr1d(k)*lamr**bm_r / am_r
        ilamr(k) = 1._dp / lamr
      else
        l_qr(k) = .false.
        rr(k) = r1
        nr(k) = r2
        mvd_r(k) = d0r
        ilamr(k) = 0._dp
        qrten(k) = -qr1d(k) * odt
        nrten(k) = -nr1d(k) * odt
        qr1d(k) = 0.0_wp
        nr1d(k) = 0.0_wp
      endif
    enddo
  end subroutine rain_check_and_update


 subroutine ice_check_and_update(rho, l_qi, qi1d, ni1d, ri, ni, &
      qiten, niten, ilami, dt, odt)
    !! computes ice contents, ilami and checks bounds
    use module_mp_tempo_params, only : max_ni, r1, r2, cie, cig, &
      mu_i, am_i, bm_i, oig1, oig2, obmi , d0s

    real(wp), intent(in) :: dt, odt
    real(wp), dimension(:), intent(in) :: rho
    real(wp), dimension(:), intent(inout) :: qi1d, ni1d, qiten, niten, ri, ni
    real(dp), dimension(:), intent(out) :: ilami
    logical, dimension(:), intent(inout) :: l_qi
    integer :: k, nz
    real(dp) :: lami, xdi
    logical :: hit_limit

    nz = size(qi1d)
    do k = 1, nz 
      hit_limit = .false.
      if (qi1d(k)+qiten(k)*dt > r1) then
        l_qi(k) = .true.
        !update mass
        ri(k) = (qi1d(k)+qiten(k)*dt)*rho(k)
        qi1d(k) = qi1d(k)+qiten(k)*dt

        !update number
        ni(k) = max(r2, (ni1d(k)+niten(k)*dt)*rho(k))
    
        ! check number
        if (ni(k) <= r2) then
          hit_limit = .true.
          lami = cie(2)/5.e-6_dp
          ni(k) = min(max_ni, cig(1)*oig2*ri(k)/am_i*lami**bm_i)
        endif

        ! check size
        lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
        xdi = (bm_i + mu_i + 1._dp) / lami
        if (xdi < 5.e-6_dp) then
          hit_limit = .true.
          lami = cie(2)/5.e-6_dp
          ni(k) = min(max_ni, cig(1)*oig2*ri(k)/am_i*lami**bm_i)
        elseif (xdi > d0s) then
          hit_limit = .true.
          lami = cie(2)/d0s
          ni(k) = cig(1)*oig2*ri(k)/am_i*lami**bm_i
        endif
        
        if (hit_limit) niten(k) = (ni(k)/rho(k) - ni1d(k))*odt
        ni1d(k) = max(r2/rho(k), &
          min(cig(1)*oig2*qi1d(k)/am_i*lami**bm_i, max_ni/rho(k)))
        ilami(k) = 1._dp / lami
      else
        l_qi(k) = .false.
        ri(k) = r1
        ni(k) = r2
        ilami(k) = 0._dp
        qiten(k) = -qi1d(k) * odt
        niten(k) = -ni1d(k) * odt
        qi1d(k) = 0.0_wp
        ni1d(k) = 0.0_wp
      endif
    enddo
  end subroutine ice_check_and_update


  subroutine snow_check_and_update(rho, l_qs, qs1d, rs, qsten, dt, odt)
    !! computes snow mass
    use module_mp_tempo_params, only : max_ni, r1, r2

    real(wp), intent(in) :: dt, odt
    real(wp), dimension(:), intent(in) :: rho
    real(wp), dimension(:), intent(inout) :: qs1d, rs, qsten
    logical, dimension(:), intent(inout) :: l_qs
    integer :: k, nz

    nz = size(qs1d)
    do k = 1, nz
      if (qs1d(k)+qsten(k)*dt > r1) then
        l_qs(k) = .true.
        ! update mass
        rs(k) = (qs1d(k)+qsten(k)*dt)*rho(k)
        qs1d(k) = qs1d(k)+qsten(k)*dt
      else
        l_qs(k) = .false.
        rs(k) = r1
        qsten(k) = -qs1d(k) * odt
        qs1d(k) = 0.0_wp
      endif
    enddo 
  end subroutine snow_check_and_update


 subroutine graupel_check_and_update(rho, l_qg, qg1d, ng1d, qb1d, rg, ng, rb, &
      idx, qgten, ngten, qbten, ilamg, mvd_g, dt, odt)
    !! computes graupel contents, ilamg, and mvd_g and checks bounds
    use module_mp_tempo_params, only : r1, r2, nrhg, rho_g, mu_g, &
      am_g, bm_g, ogg3, cgg, ogg2, obmg, d0r, idx_bg1, gonv_max, &
      gonv_min, oge1, ogg1, d0g, meters3_to_liters

    real(wp), intent(in) :: dt, odt    
    real(wp), dimension(:), intent(in) :: rho
    real(wp), dimension(:), intent(inout) :: qg1d, qgten, rg, ng, rb, ngten, qbten
    real(wp), dimension(:), intent(inout), optional :: ng1d, qb1d
    real(dp), dimension(:), intent(out) :: ilamg
    real(wp), dimension(:), intent(out) :: mvd_g
    logical, dimension(:), intent(inout) :: l_qg
    integer, dimension(:), intent(inout) :: idx
    integer :: k, nz
    real(dp) :: lamg, ygra1, zans1, n0_exp, lam_exp
    logical :: hit_limit

    nz = size(qg1d)
    do k = 1, nz
      hit_limit = .false.
      if (qg1d(k)+qgten(k)*dt > r1) then
        l_qg(k) = .true.
        !update mass
        rg(k) = (qg1d(k)+qgten(k)*dt)*rho(k)
        qg1d(k) = qg1d(k)+qgten(k)*dt

        !update number and density
        if (present(ng1d) .and. present(qb1d)) then
          ng(k) = max(r2, (ng1d(k)+ngten(k)*dt)*rho(k))
          ! qb1d is L/kg and rb is L/m^3
          rb(k) = min(max(rg(k)*meters3_to_liters/rho_g(nrhg), &
            (qb1d(k)+qbten(k)*dt)*rho(k)), rg(k)*meters3_to_liters/rho_g(1))
          idx(k) = max(1, min(nint(10._wp*rg(k)/rb(k))+1, nrhg))

          ! check number
          if (ng(k) <= r2) then
            hit_limit = .true.
            mvd_g(k) = 1.5e-3_wp
            lamg = (3.0_dp + mu_g + 0.672_dp) / mvd_g(k)
            ng(k) = cgg(2,1)*ogg3*rg(k)*lamg**bm_g / am_g(idx(k))
          endif

          ! check size
          lamg = (am_g(idx(k))*cgg(3,1)*ogg2*ng(k)/rg(k))**obmg
          mvd_g(k) = (3.0_wp + mu_g + 0.672_wp) / lamg
          if (mvd_g(k) > 25.4e-3_wp) then
            hit_limit = .true.
            mvd_g(k) = 25.4e-3_wp
            lamg = (3.0_dp + mu_g + 0.672_dp) / mvd_g(k)
            ng(k) = cgg(2,1)*ogg3*rg(k)*lamg**bm_g / am_g(idx(k))
          elseif (mvd_g(k) < d0r) then
            hit_limit = .true.
            mvd_g(k) = d0r
            lamg = (3.0_dp + mu_g + 0.672_dp) / mvd_g(k)
            ng(k) = cgg(2,1)*ogg3*rg(k)*lamg**bm_g / am_g(idx(k))
          endif

          if (hit_limit) ngten(k) = (ng(k)/rho(k) - ng1d(k)) * odt
          ng1d(k) = cgg(2,1)*ogg3*qg1d(k)*lamg**bm_g / am_g(idx(k))
          qb1d(k) = min(max(qg1d(k)*meters3_to_liters/rho_g(nrhg), &
            qb1d(k)+qbten(k)*dt), meters3_to_liters*qg1d(k)/rho_g(1))
          idx(k) = max(1, min(nint(10._wp*qg1d(k)/qb1d(k))+1, nrhg))
        else
          idx(k) = idx_bg1
          ygra1 = log10(max(1.e-9_dp, real(rg(k), kind=dp)))
          zans1 = 3.4_dp + 2._dp/7._dp*(ygra1+8._dp)
          n0_exp = max(gonv_min, min(10._dp**(zans1), gonv_max))
          lam_exp = (n0_exp*am_g(idx(k))*cgg(1,1)/rg(k))**oge1
          lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
          ng(k) = cgg(2,1)*ogg3*rg(k)*lamg**bm_g / am_g(idx(k))
          rb(k) = meters3_to_liters*rg(k)/rho_g(idx(k))
        endif
        ilamg(k) = 1._dp / lamg
        mvd_g(k) = (3.0_wp + mu_g + 0.672_wp) * ilamg(k)
      else
        l_qg(k) = .false.
        rg(k) = r1
        ng(k) = r2
        mvd_g(k) = d0g
        ilamg(k) = 0._dp
        idx(k) = idx_bg1
        rb(k) = meters3_to_liters*r1/rho_g(idx(k))
        qgten(k) = -qg1d(k) * odt
        qg1d(k) = 0.0_wp
        if (present(ng1d) .and. present(qb1d)) then
          ngten(k) = -ng1d(k) * odt
          qbten(k) = -qb1d(k) * odt
          ng1d(k) = 0.0_wp
          qb1d(k) = 0.0_wp
        endif 
      endif
    enddo 
  end subroutine graupel_check_and_update


  subroutine graupel_init(rho, qg1d, ng1d, qb1d)
    !! initializes graupel number and volume if both are zero
    !! and hail-aware = true
    use module_mp_tempo_params, only : r1, meters3_to_liters, &
      idx_bg1, am_g, bm_g, mu_g, ogg3, cgg

    real(wp), dimension(:), intent(in) :: rho
    real(wp), dimension(:), intent(in) :: qg1d
    real(wp), dimension(:), intent(inout) :: ng1d, qb1d
    integer :: k, nz, idx
    real(dp) :: lamg
    real(wp) :: mvd_g, rg, ng, rb

    nz = size(rho)
    do k = 1, nz
      if (qg1d(k) > r1) then
        rg = qg1d(k)*rho(k)
        rb = rg*meters3_to_liters/rho_g(idx_bg1)
        idx = max(1, min(nint(10._wp*rg/rb)+1, nrhg))
        mvd_g = 5.e-3_wp
        lamg = (3.0_dp + mu_g + 0.672_dp) / mvd_g
        ng = cgg(2,1)*ogg3*rg*lamg**bm_g / am_g(idx)
        ng1d(k) = ng/rho(k)
        qb1d(k) = rb/rho(k)
      endif
    enddo 
  end subroutine graupel_init


  subroutine thermo_vars(qv, temp, pres, rho, rhof, rhof2, qvs, delqvs, qvsi, &
    satw, sati, ssatw, ssati, diffu, visco, vsc2, ocp, lvap, tcond, lvt2, &
    supersaturated)
    !! computes thermodynamic variables
    use module_mp_tempo_params, only : t0, rho_not, eps, cp, lvap0, orv

    real(wp), dimension(:), intent(in) :: qv, temp, pres, rho
    real(wp), dimension(:), intent(out) :: rhof, rhof2, qvs, &
      delqvs, qvsi, satw, sati, ssatw, ssati, diffu, visco, vsc2, &
      ocp, lvap, tcond, lvt2
    logical, intent(inout) :: supersaturated
    integer :: k, nz
    real(wp) :: tempc, otemp

    nz = size(temp)
    supersaturated = .false.
    do k = 1, nz
      otemp = 1._wp / temp(k)
      tempc = temp(k) - t0
      rhof(k) = sqrt(rho_not/rho(k))
      rhof2(k) = sqrt(rhof(k))
      qvs(k) = calc_rslf(pres(k), temp(k))
      delqvs(k) = max(0._wp, calc_rslf(pres(k), t0)-qv(k))
      if (tempc <= 0._wp) then
        qvsi(k) = calc_rsif(pres(k), temp(k))
      else
        qvsi(k) = qvs(k)
      endif
      satw(k) = qv(k)/qvs(k)
      sati(k) = qv(k)/qvsi(k)
      ssatw(k) = satw(k) - 1._wp
      ssati(k) = sati(k) - 1._wp
      if (abs(ssatw(k)) < eps) ssatw(k) = 0._wp
      if (abs(ssati(k)) < eps) ssati(k) = 0._wp
      if (ssati(k) > 0._wp) supersaturated = .true.
      diffu(k) = 2.11e-5_wp*(temp(k)/t0)**1.94_wp * (101325._wp/pres(k))
      if (tempc >= 0._wp) then
        visco(k) = (1.718_wp+0.0049_wp*tempc)*1.0e-5_wp
      else
        visco(k) = (1.718_wp+0.0049_wp*tempc-1.2e-5_wp*tempc*tempc)*1.0e-5_wp
      endif
      ocp(k) = 1._wp/(cp*(1._wp+0.887_wp*qv(k)))
      vsc2(k) = sqrt(rho(k)/visco(k))
      lvap(k) = lvap0 + (2106.0_wp - 4218.0_wp)*tempc
      tcond(k) = (5.69_wp + 0.0168_wp*tempc)*1.0e-5_wp * 418.936_wp
      lvt2(k) = lvap(k)*lvap(k)*ocp(k)*orv*otemp*otemp
    enddo
  end subroutine thermo_vars


  subroutine check_over_depletion(rho, temp, qvsi, qv, l_qc, rc, l_qi, ri, &
    l_qr, rr, l_qs, rs, l_qg, rg, tend, odt)
    !! check to ensure that loss terms don't over-deplete a category and
    !! adjusts tendencies if needed
    use module_mp_tempo_params, only : eps, rho_i, t0, meters3_to_liters

    real(wp), intent(in) :: odt
    type(ty_tend), intent(inout) :: tend
    logical, dimension(:), intent(in) :: l_qc, l_qi, l_qr, l_qs, l_qg
    real(wp), dimension(:), intent(in) :: rho, temp, qvsi, qv, rc, ri, rr, rs, rg
    real(wp) :: sump, rate_max, ratio
    integer :: k, nz

    nz = size(qv)
    do k = 1, nz
      ! losses to vapor include, deposition nucleation; ice, snow, and graupel
      ! depositional growth; Koop nucleation
      sump = tend%pri_inu(k) + tend%pri_ide(k) + tend%prs_ide(k) + &
        tend%prs_sde(k) + tend%prg_gde(k) + tend%pri_iha(k)
      rate_max = (qv(k)-qvsi(k))*rho(k)*odt*0.999_wp
      if ((sump > eps .and. sump > rate_max) .or. &
        (sump < -eps .and. sump < rate_max)) then
        ratio = rate_max/sump
        tend%pri_inu(k) = tend%pri_inu(k) * ratio
        tend%pri_ide(k) = tend%pri_ide(k) * ratio
        tend%pni_ide(k) = tend%pni_ide(k) * ratio ! scales with pri_ide
        tend%prs_ide(k) = tend%prs_ide(k) * ratio
        tend%prs_sde(k) = tend%prs_sde(k) * ratio
        tend%prg_gde(k) = tend%prg_gde(k) * ratio
        tend%pri_iha(k) = tend%pri_iha(k) * ratio
      endif

      ! losses to cloud water include conversion to rain;
      ! freezing; collection by rain, snow, and graupel.
      sump = -tend%prr_wau(k) - tend%pri_wfz(k) - tend%prr_rcw(k) - &
        tend%prs_scw(k) - tend%prg_scw(k) - tend%prg_gcw(k)
      rate_max = -rc(k)*odt
      if (l_qc(k)) then
        if (sump < rate_max) then
          ratio = rate_max/sump
          tend%prr_wau(k) = tend%prr_wau(k) * ratio
          tend%pri_wfz(k) = tend%pri_wfz(k) * ratio
          tend%prr_rcw(k) = tend%prr_rcw(k) * ratio
          tend%prs_scw(k) = tend%prs_scw(k) * ratio
          tend%prg_scw(k) = tend%prg_scw(k) * ratio
          tend%prg_gcw(k) = tend%prg_gcw(k) * ratio
        endif
      endif 

      ! losses to cloud ice include sublimation; conversion to snow;
      ! collection by snow and rain
      sump = tend%pri_ide(k) - tend%prs_iau(k) - tend%prs_sci(k) - tend%pri_rci(k)
      rate_max = -ri(k)*odt
      if (l_qi(k)) then
        if (sump < rate_max) then
          ratio = rate_max/sump
          tend%pri_ide(k) = tend%pri_ide(k) * ratio
          tend%prs_iau(k) = tend%prs_iau(k) * ratio
          tend%prs_sci(k) = tend%prs_sci(k) * ratio
          tend%pri_rci(k) = tend%pri_rci(k) * ratio
        endif
      endif 

      ! losses to rain include freezing; collection by ice, snow, and graupel
      ! resulting in freezing
      sump = -tend%prg_rfz(k) - tend%pri_rfz(k) - tend%prr_rci(k) + &
        tend%prr_rcs(k) + tend%prr_rcg(k)
      rate_max = -rr(k)*odt
      if (l_qr(k)) then
        if (sump < rate_max) then
          ratio = rate_max/sump
          tend%prg_rfz(k) = tend%prg_rfz(k) * ratio
          tend%pbg_rfz(k) = tend%pbg_rfz(k) * ratio ! scales with prg_rfz
          tend%pri_rfz(k) = tend%pri_rfz(k) * ratio
          tend%prr_rci(k) = tend%prr_rci(k) * ratio
          tend%prr_rcs(k) = tend%prr_rcs(k) * ratio
          tend%prr_rcg(k) = tend%prr_rcg(k) * ratio
        endif
      endif

      ! losses to snow include sublimation; melting; collection by rain; 
      ! rime splintering
      sump = tend%prs_sde(k) - tend%prs_ihm(k) - tend%prr_sml(k) + &
        tend%prs_rcs(k)
      rate_max = -rs(k)*odt
      if (l_qs(k)) then 
        if (sump < rate_max) then
          ratio = rate_max/sump
          tend%prs_sde(k) = tend%prs_sde(k) * ratio
          tend%prs_ihm(k) = tend%prs_ihm(k) * ratio
          tend%prr_sml(k) = tend%prr_sml(k) * ratio
          tend%prs_rcs(k) = tend%prs_rcs(k) * ratio
        endif
      endif

      ! losses to graupel include sublimation; melting; rime splintering; 
      ! collection by rain
      sump = tend%prg_gde(k) - tend%prg_ihm(k) - tend%prr_gml(k) + tend%prg_rcg(k)
      rate_max = -rg(k)*odt
      if (l_qg(k)) then
        if (sump < rate_max) then
          ratio = rate_max/sump
          tend%prg_gde(k) = tend%prg_gde(k) * ratio
          tend%prg_ihm(k) = tend%prg_ihm(k) * ratio
          tend%prr_gml(k) = tend%prr_gml(k) * ratio
          tend%prg_rcg(k) = tend%prg_rcg(k) * ratio
          tend%pbg_rcg(k) = tend%pbg_rcg(k) * ratio ! scales with prg_rcg
        endif
      endif 

      ! updates after adjustments
      ! reset sum of rime splintering
      tend%pri_ihm(k) = tend%prs_ihm(k) + tend%prg_ihm(k)
      ! reset total rain-graupel collection amount if reduced
      ratio = min(abs(tend%prr_rcg(k)), abs(tend%prg_rcg(k)))
      tend%prr_rcg(k) = ratio * sign(1.0_dp, tend%prr_rcg(k))
      tend%prg_rcg(k) = -tend%prr_rcg(k)
      tend%pbg_rcg(k) = meters3_to_liters*tend%prg_rcg(k)/rho_i ! scale density change
      ! reset total rain-snow collection amount if reduced
      if (temp(k) > t0) then
        ratio = min(abs(tend%prr_rcs(k)), abs(tend%prs_rcs(k)))
        tend%prr_rcs(k) = ratio * sign(1.0_dp, tend%prr_rcs(k))
        tend%prs_rcs(k) = -tend%prr_rcs(k)
      endif
    enddo
  end subroutine check_over_depletion


  subroutine sum_tendencies(rho, temp, idx, lvap, ocp, tend, tten, qvten, qcten, &
    ncten, qiten, niten, qsten, qrten, nrten, qgten, ngten, qbten)
    !! sums tendencies for each hydrometeor category and temperature and moisture
    use module_mp_tempo_params, only : lsub, rho_g, t0, lfus, meters3_to_liters

    type(ty_tend), intent(in) :: tend
    real(wp), dimension(:), intent(in) :: rho, temp, lvap, ocp
    integer, dimension(:), intent(in) :: idx
    real(wp), dimension(:), intent(inout) :: qvten, qcten, ncten, qiten, niten, &
      qsten, qrten, nrten, qgten, ngten, qbten, tten
    real(wp) :: orho, lfus2
    integer :: k, nz

    nz = size(temp)
    do k = 1, nz
      orho = 1./rho(k)
      lfus2 = lsub - lvap(k)

      qvten(k) = qvten(k) + (-tend%pri_inu(k) - tend%pri_iha(k) - tend%pri_ide(k) - &
        tend%prs_ide(k) - tend%prs_sde(k) - tend%prg_gde(k)) * orho

      qcten(k) = qcten(k) + (-tend%prr_wau(k) - tend%pri_wfz(k) - tend%prr_rcw(k) - &
        tend%prs_scw(k) - tend%prg_scw(k) - tend%prg_gcw(k)) * orho

      ncten(k) = ncten(k) + (-tend%pnc_wau(k) - tend%pnc_rcw(k) - tend%pni_wfz(k) - &
        tend%pnc_scw(k) - tend%pnc_gcw(k)) * orho

      qiten(k) = qiten(k) + (tend%pri_inu(k) + tend%pri_iha(k) + tend%pri_ihm(k) + &
        tend%pri_wfz(k) + tend%pri_rfz(k) + tend%pri_ide(k) - tend%prs_iau(k) - &
        tend%prs_sci(k) - tend%pri_rci(k)) * orho

      niten(k) = niten(k) + (tend%pni_inu(k) + tend%pni_iha(k) + tend%pni_ihm(k) + &
        tend%pni_wfz(k) + tend%pni_rfz(k) + tend%pni_ide(k) - tend%pni_iau(k) - &
        tend%pni_sci(k) - tend%pni_rci(k)) * orho

      qrten(k) = qrten(k) + (tend%prr_wau(k) + tend%prr_rcw(k) + tend%prr_sml(k) + &
        tend%prr_gml(k) + tend%prr_rcs(k) + tend%prr_rcg(k) - tend%prg_rfz(k) - &
        tend%pri_rfz(k) - tend%prr_rci(k)) * orho

      nrten(k) = nrten(k) + (tend%pnr_wau(k) + tend%pnr_sml(k) + tend%pnr_gml(k) - &
        (tend%pnr_rfz(k) + tend%pnr_rcr(k) + tend%pnr_rcg(k) + tend%pnr_rcs(k) + &
        tend%pnr_rci(k) + tend%pni_rfz(k))) * orho

      qsten(k) = qsten(k) + (tend%prs_iau(k) + tend%prs_sde(k) + tend%prs_sci(k) + &
        tend%prs_scw(k) + tend%prs_rcs(k) + tend%prs_ide(k) - tend%prs_ihm(k) - &
        tend%prr_sml(k)) * orho

      qgten(k) = qgten(k) + (tend%prg_scw(k) + tend%prg_rfz(k) + tend%prg_gde(k) + &
        tend%prg_rcg(k) + tend%prg_gcw(k) + tend%prg_rci(k) + tend%prg_rcs(k) - &
        tend%prg_ihm(k) - tend%prr_gml(k)) * orho

      ngten(k) = ngten(k) + (tend%png_scw(k) + tend%png_rfz(k) - tend%png_rcg(k) + &
        tend%png_rci(k) + tend%png_rcs(k) + tend%png_gde(k) - tend%pnr_gml(k)) * orho

      qbten(k) = qbten(k) + (tend%pbg_scw(k) + tend%pbg_rfz(k) + tend%pbg_gcw(k) + &
        tend%pbg_rci(k) + tend%pbg_rcs(k) + tend%pbg_rcg(k) + tend%pbg_sml(k) - &
        tend%pbg_gml(k) + meters3_to_liters * (tend%prg_gde(k) - tend%prg_ihm(k)) / rho_g(idx(k))) * orho

      if (temp(k) < t0) then
        tten(k) = tten(k) + &
          (lsub*ocp(k)*(tend%pri_inu(k) + tend%pri_ide(k) + &
          tend%prs_ide(k) + tend%prs_sde(k) + tend%prg_gde(k) + tend%pri_iha(k)) + &
          lfus2*ocp(k)*(tend%pri_wfz(k) + tend%pri_rfz(k) + tend%prg_rfz(k) + &
          tend%prs_scw(k) + tend%prg_scw(k) + tend%prg_gcw(k) + tend%prg_rcs(k) + &
          tend%prs_rcs(k) + tend%prr_rci(k) + tend%prg_rcg(k)))*orho
      else
        tten(k) = tten(k) + &
          (lfus*ocp(k)*(-tend%prr_sml(k) - tend%prr_gml(k) - &
          tend%prr_rcg(k) - tend%prr_rcs(k)) + &
          lsub*ocp(k)*(tend%prs_sde(k) + tend%prg_gde(k)))*orho
      endif
    enddo
  end subroutine sum_tendencies


  subroutine sedimentation(xr, vt, dz1d, rho, xten, limit, steps, ktop_sedi, precip, dt)
    !! computes sedimentation fluxes, adds fluxes to tendencies, and updates hydrometeor
    !! mass (and number and volume)
    use module_mp_tempo_params, only : low_limit_mass_for_precip

    real(wp), intent(in) :: dt
    integer, intent(in) :: steps
    integer, intent(in), optional :: ktop_sedi
    real(wp), dimension(:), intent(inout) :: xr, xten
    real(wp), dimension(:), intent(in) :: dz1d, rho
    real(wp), dimension(:), intent(in) :: vt
    real(wp), intent(in) :: limit
    real(wp), intent(inout), optional :: precip
    real(wp) :: odz, orho
    real(wp), allocatable, dimension(:) :: sed_r
    integer :: k, nz, ktop

    nz = size(xr)
    allocate(sed_r(nz), source=0._wp)
    ktop = nz-1
    if (present(ktop_sedi)) ktop = ktop_sedi

    do k = nz, 1, -1
      sed_r(k) = vt(k)*xr(k)
    enddo
    k = nz
    odz = 1._wp/dz1d(k)
    orho = 1._wp/rho(k)
    xten(k) = xten(k) - sed_r(k)*odz*(1._wp/real(steps, kind=wp))*orho
    xr(k) = max(limit, xr(k) - sed_r(k)*odz*dt*(1._wp/real(steps, kind=wp)))

    do k = ktop, 1, -1
      odz = 1._wp/dz1d(k)
      orho = 1._wp/rho(k)
      xten(k) = xten(k) + (sed_r(k+1)-sed_r(k))*odz*(1._wp/real(steps, kind=wp))*orho
      xr(k) = max(limit, xr(k) + (sed_r(k+1)-sed_r(k))*odz*dt*(1._wp/real(steps, kind=wp)))
    enddo

    if (present(precip)) then 
      if (xr(1) > low_limit_mass_for_precip) then
        precip = precip + sed_r(1)*dt*(1._wp/real(steps, kind=wp))
      endif 
    endif 
  end subroutine sedimentation


  subroutine semilagrangian_sedimentation(dz1d, rho, xr, xten, vt, steps, limit, precip, dt, odt)
    !! semi-lagrangian sedimentation scheme from 
    !! [Juang and Hong (2010)](https://doi.org/10.1175/2009MWR3109.1)
    !!
    !! original author: hann-ming henry juang <henry.juang@noaa.gov>
    !! original implemented by: song-you hong

    real(wp), intent(in) :: dt, odt
    integer, intent(in) :: steps
    real(wp), intent(in) :: limit
    real(wp), dimension(:), intent(in) :: dz1d, rho
    real(wp), dimension(:), intent(inout) :: xr, xten
    real(wp), dimension(:), intent(in) :: vt
    real(wp) :: fa1, fa2, con1, decfl, dip, dim, tl, th, tl2, th2, qqd, qqh, qql, zsum, qsum, &
      orho, dql, dqh
    real(wp), dimension(:), allocatable :: zi, wi, za, dza, qa, qmi, qpi, net_flx, precip_flx, rr_save
    real(wp), intent(out), optional :: precip
    integer :: k, nz, kk, kb, kt, m

    nz = size(dz1d)
    allocate(zi(nz+1), source=0._wp)
    allocate(wi(nz+1), source=0._wp)
    allocate(dza(nz+1), source=0._wp)
    allocate(za(nz+2), source=0._wp)
    allocate(qa(nz+1), source=0._wp)
    allocate(qmi(nz+1), source=0._wp)
    allocate(qpi(nz+1), source=0._wp)
    allocate(net_flx(nz), source=0._wp)
    allocate(precip_flx(nz), source=0._wp)
    allocate(rr_save(nz), source=xr)

    zi(1) = 0._wp ! zi(1) needs to be zero so zero out explicitly
    do k = 1, nz
      zi(k+1) = zi(k) + dz1d(k)
    enddo
    
    ! 3rd order interpolation to get wi
    fa1 = 9._wp/16._wp
    fa2 = 1._wp/16._wp
    wi(1) = vt(1)
    wi(2) = 0.5_wp*(vt(2)+vt(1))
    do k = 3, nz-1
      wi(k) = fa1*(vt(k)+vt(k-1))-fa2*(vt(k+1)+vt(k-2))
    enddo
    wi(nz) = 0.5_wp*(vt(nz)+vt(nz-1))
    wi(nz+1) = vt(nz+1)

    do k = 2, nz
      if(vt(k) == 0._wp) wi(k) = vt(k-1)
    enddo

    ! diffusivity of wi
    ! con1 should be > 0 and < 1
    con1 = 0.05_wp
    do k = nz, 1, -1
      decfl = (wi(k+1)-wi(k))*dt*(1._wp/real(steps, kind=wp))/dz1d(k)
      if(decfl > con1) then
        wi(k) = wi(k+1) - con1*dz1d(k)*odt*real(steps, kind=wp)
      endif
    enddo

    ! compute arrival point
    do k = 1, nz+1
      za(k) = zi(k) - wi(k)*dt*(1._wp/real(steps, kind=wp))
    enddo
    za(nz+2) = zi(nz+1)
    do k = 1, nz+1
      dza(k) = za(k+1)-za(k)
    enddo

  ! computer deformation at arrival point
    do k = 1, nz
      qa(k) = xr(k)*dz1d(k)/dza(k)
    enddo
    qa(nz+1) = 0._wp

    ! estimate values at arrival cell interface with monotone
    do k = 2, nz
      dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
      dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
      if(dip*dim <= 0._wp) then
        qmi(k)=qa(k)
        qpi(k)=qa(k)
      else
        qpi(k)=qa(k)+0.5_wp*(dip+dim)*dza(k)
        qmi(k)=2._wp*qa(k)-qpi(k)
        if(qpi(k) > 0._wp .or. qmi(k) < 0._wp) then
          qpi(k) = qa(k)
          qmi(k) = qa(k)
        endif
      endif
    enddo
    qpi(1) = qa(1)
    qmi(1) = qa(1)
    qmi(nz+1) = qa(nz+1)
    qpi(nz+1) = qa(nz+1)

    ! interpolation to regular point
    kb = 1
    kt = 1
    intp : do k = 1, nz
      kb = max(kb-1,1)
      kt = max(kt-1,1)
      ! find kb and kt
      if(zi(k) >= za(nz+1)) then
        exit intp
      else
        find_kb : do kk = kb, nz
          if(zi(k) <= za(kk+1)) then
            kb = kk
            exit find_kb
          else
            cycle find_kb
          endif
        enddo find_kb
        find_kt : do kk = kt, nz+2
          if(zi(k+1) <= za(kk)) then
            kt = kk
            exit find_kt
          else  
            cycle find_kt
          endif
        enddo find_kt
        kt = kt - 1
        
        ! compute q with piecewise constant method
        if(kt == kb) then
          tl = (zi(k)-za(kb))/dza(kb)
          th = (zi(k+1)-za(kb))/dza(kb)
          tl2 = tl*tl
          th2 = th*th
          qqd = 0.5_wp*(qpi(kb)-qmi(kb))
          qqh = qqd*th2+qmi(kb)*th
          qql = qqd*tl2+qmi(kb)*tl
          xr(k) = (qqh-qql)/(th-tl)
        elseif(kt > kb) then
          tl = (zi(k)-za(kb))/dza(kb)
          tl2 = tl*tl
          qqd = 0.5_wp*(qpi(kb)-qmi(kb))
          qql = qqd*tl2+qmi(kb)*tl
          dql = qa(kb)-qql
          zsum = (1._wp-tl)*dza(kb)
          qsum = dql*dza(kb)
          if(kt-kb > 1) then
            do m = kb+1, kt-1
              zsum = zsum + dza(m)
              qsum = qsum + qa(m) * dza(m)
            enddo
          endif
          th = (zi(k+1)-za(kt))/dza(kt)
          th2 = th*th
          qqd = 0.5_wp*(qpi(kt)-qmi(kt))
          dqh = qqd*th2+qmi(kt)*th
          zsum = zsum + th*dza(kt)
          qsum = qsum + dqh*dza(kt)
          xr(k) = qsum/zsum
        endif
      endif
      orho = 1._wp / rho(k)
      xr(k) = max(xr(k), limit)
      xten(k) = xten(k) + (xr(k) - rr_save(k)) * &
        orho*odt
    enddo intp

    precip_loop: do k = 1, nz
      if(za(k) < 0._wp .and. za(k+1) <= 0.0_wp) then
        if (present(precip)) precip = precip + qa(k)*dza(k)
        net_flx(k) =  qa(k)*dza(k)
      elseif (za(k) < 0._wp .and. za(k+1) > 0._wp) then
        th = (0._wp-za(k))/dza(k)
        th2 = th*th
        qqd = 0.5_wp*(qpi(k)-qmi(k))
        qqh = qqd*th2+qmi(k)*th
        if (present(precip)) precip = precip + qqh*dza(k)
        net_flx(k) = qqh*dza(k)
        exit precip_loop
      endif
    enddo precip_loop

    ! calculating precipitation fluxes
    do k = nz, 1, -1
      if(k == nz) then
        precip_flx(k) = net_flx(k)
      else
        precip_flx(k) = precip_flx(k+1) + net_flx(k)
      end if
    enddo
  end subroutine semilagrangian_sedimentation


  subroutine rain_fallspeed(rhof, l_qr, rr, ilamr, dz1d, vt, vtn, substeps_sedi, ktop_sedi, dt)
    !! calculates mass and number weighted fall speeds for rain
    !! and optionally the substepping required and the top k-level of sedimentation
    use module_mp_tempo_params, only : crg, av_r, org3, fv_r, cre

    real(wp), intent(in) :: dt
    real(wp), dimension(:), intent(in) :: rhof, dz1d, rr
    real(dp), dimension(:), intent(in) :: ilamr
    logical, dimension(:), intent(in) :: l_qr
    real(wp), dimension(:), intent(inout) :: vt, vtn
    integer, intent(out), optional :: substeps_sedi, ktop_sedi
    real(wp) :: dz_by_vt
    real(dp) :: lamr
    integer :: k, nz

    nz = size(l_qr)
    if (present(ktop_sedi)) ktop_sedi = 1
    if (present(substeps_sedi)) substeps_sedi = 1
    do k = nz, 1, -1
      if (rr(k) > r1) then
        lamr = 1._dp / ilamr(k)
        vt(k) = rhof(k)*av_r*crg(6)*org3 * lamr**cre(3) *((lamr+fv_r)**(-cre(6)))
        vtn(k) = rhof(k)*av_r*crg(7)/crg(12) * lamr**cre(12)*((lamr+fv_r)**(-cre(7)))
      else
        vt(k) = vt(k+1)
        vtn(k) = vtn(k+1)
      endif
      if (max(vt(k), vtn(k)) > 1.e-3_wp) then
        if (present(ktop_sedi)) ktop_sedi = max(ktop_sedi, k)
        dz_by_vt = dz1d(k) / (max(vt(k), vtn(k)))
        if (present(substeps_sedi)) then
          substeps_sedi = max(substeps_sedi, int(dt/dz_by_vt + 1._wp))
        endif 
      endif
    enddo
    if (present(ktop_sedi)) then
      if (ktop_sedi == nz) ktop_sedi = nz-1 
    endif 
  end subroutine rain_fallspeed


  subroutine graupel_fallspeed(rhof, rho, visco, l_qg, rg, rb, qb1d, idx, ilamg, &
      dz1d, vt, vtn, substeps_sedi, ktop_sedi, dt)
    !! calculates mass and number weighted fall speeds for graupel
    !! and optionally the substepping required and the top k-level of sedimentation
    use module_mp_tempo_params, only : nrhg, rho_g, av_g_old, bv_g_old, &
      cgg, t0, mu_g, ogg2, ogg3, a_coeff, b_coeff, meters3_to_liters, earth_gravity

    real(wp), intent(in) :: dt    
    real(wp), dimension(:), intent(in) :: rhof, rho, visco, dz1d, rg, rb
    real(wp), dimension(:), intent(in), optional :: qb1d
    real(dp), dimension(:), intent(in) :: ilamg
    logical, dimension(:), intent(in) :: l_qg
    integer, dimension(:), intent(in) :: idx
    real(wp), dimension(:), intent(inout) :: vt, vtn
    integer, intent(out), optional :: substeps_sedi, ktop_sedi
    real(wp) :: dz_by_vt, dens_g, afall, bfall
    integer :: k, nz

    nz = size(l_qg)
    if (present(ktop_sedi)) ktop_sedi = 1
    if (present(substeps_sedi)) substeps_sedi = 1
    do k = nz, 1, -1
      if (rg(k) > r1) then
        if (present(qb1d)) then
          dens_g = max(rho_g(1), min(meters3_to_liters*rg(k)/rb(k), rho_g(nrhg)))
          afall = a_coeff*((4._wp*dens_g*earth_gravity)/(3._wp*rho(k)))**b_coeff
          afall = afall * visco(k)**(1._wp-2._wp*b_coeff)
          bfall = 3._wp*b_coeff - 1._wp
        else
          afall = av_g_old
          bfall = bv_g_old
        endif
        vt(k) = rhof(k)*afall*cgg(6,idx(k))*ogg3 * ilamg(k)**bfall
        ! idea: if (temp(k) > t0) vt(k) = max(vt(k), vtrr(k))

        if (mu_g == 0) then
          vtn(k) = rhof(k)*afall*cgg(7,idx(k))/cgg(12,idx(k)) * ilamg(k)**bfall
        else
          vtn(k) = rhof(k)*afall*cgg(8,idx(k))*ogg2 * ilamg(k)**bfall
        endif
      else
        vt(k) = vt(k+1)
        vtn(k) = vtn(k+1)
      endif
      if (vt(k) > 1.e-3_wp) then
        if (present(ktop_sedi)) ktop_sedi = max(ktop_sedi, k)
        dz_by_vt = dz1d(k) / vt(k)
        if (present(substeps_sedi)) then
          substeps_sedi = max(substeps_sedi, int(dt/dz_by_vt + 1._wp))
        endif 
      endif
    enddo
    if (present(ktop_sedi)) then
      if (ktop_sedi == nz) ktop_sedi = nz-1 
    endif 
  end subroutine graupel_fallspeed


  subroutine snow_fallspeed(rhof, l_qs, rs, prr_sml, smob, smoc, &
      rr, vtrr, dz1d, vt, vtboost, substeps_sedi, ktop_sedi, dt)
    !! calcules mass weighted fall speeds for snow
    !! and optionally the substepping required and the top k-level of sedimentation
    use module_mp_tempo_params, only : lam0, lam1, fv_s, kap0, kap1, mu_s, &
      cse, csg, av_s

    real(wp), intent(in) :: dt
    real(wp), dimension(:), intent(in) :: rhof, dz1d, rs, rr, vtboost
    real(wp), dimension(:), intent(in) :: vtrr
    logical, dimension(:), intent(in) :: l_qs
    real(wp), dimension(:), intent(inout) :: vt
    real(dp), dimension(:), intent(in) :: smob, smoc, prr_sml
    integer, intent(out), optional :: substeps_sedi, ktop_sedi
    real(wp) :: dz_by_vt, vts, sr
    real(dp) :: xds, mrat, ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts
    integer :: k, nz

    nz = size(l_qs)
    if (present(ktop_sedi)) ktop_sedi = 1
    if (present(substeps_sedi)) substeps_sedi = 1
    do k = nz, 1, -1
      if (rs(k) > r1) then
        xds = smoc(k) / smob(k)
        mrat = 1._dp/xds
        ils1 = 1._dp/(mrat*lam0 + fv_s)
        ils2 = 1._dp/(mrat*lam1 + fv_s)
        t1_vts = kap0*csg(4)*ils1**cse(4)
        t2_vts = kap1*mrat**mu_s*csg(10)*ils2**cse(10)
        ils1 = 1._dp/(mrat*lam0)
        ils2 = 1._dp/(mrat*lam1)
        t3_vts = kap0*csg(1)*ils1**cse(1)
        t4_vts = kap1*mrat**mu_s*csg(7)*ils2**cse(7)
        vts = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)

        if (prr_sml(k) > 0._dp) then
          sr = rs(k)/(rs(k)+rr(k))
          vt(k) = vts*sr + (1._wp-sr)*vtrr(k)
        else
          vt(k) = vts*vtboost(k)
        endif
      else
        vt(k) = vt(k+1)
      endif
      if (vt(k) > 1.e-3_wp) then
        if (present(ktop_sedi)) ktop_sedi = max(ktop_sedi, k)
        dz_by_vt = dz1d(k) / vt(k)
        if (present(substeps_sedi)) then
          substeps_sedi = max(substeps_sedi, int(dt/dz_by_vt + 1._wp))
        endif 
      endif
    enddo
    if (present(ktop_sedi)) then
      if (ktop_sedi == nz) ktop_sedi = nz-1
    endif 
  end subroutine snow_fallspeed


  subroutine ice_fallspeed(rhof, l_qi, ri, ilami, dz1d, vt, vtn, &
    substeps_sedi, ktop_sedi, dt)
    !! calculates mass and number weighted fall speeds for ice
    !! and the substepping required and the top k-level of sedimentation
    use module_mp_tempo_params, only : av_i, cig, oig2, bv_i

    real(wp), intent(in) :: dt    
    real(wp), dimension(:), intent(in) :: rhof, dz1d, ri
    real(dp), dimension(:), intent(in) :: ilami
    logical, dimension(:), intent(in) :: l_qi
    real(wp), dimension(:), intent(inout) :: vt, vtn
    integer, intent(out) :: substeps_sedi, ktop_sedi
    real(wp) :: dz_by_vt
    integer :: k, nz

    nz = size(l_qi)
    ktop_sedi = 1
    substeps_sedi = 1
    do k = nz, 1, -1
      if (ri(k) > r1) then
        vt(k) = rhof(k)*av_i*cig(3)*oig2 * ilami(k)**bv_i
        vtn(k) = rhof(k)*av_i*cig(6)/cig(7) * ilami(k)**bv_i
      else
        vt(k) = vt(k+1)
        vtn(k) = vtn(k+1)
      endif
      if (vt(k) > 1.e-3_wp) then
        ktop_sedi = max(ktop_sedi, k)
        dz_by_vt = dz1d(k) / vt(k)
        substeps_sedi = max(substeps_sedi, int(dt/dz_by_vt + 1._wp))
      endif
    enddo
    if (ktop_sedi == nz) ktop_sedi = nz-1 
  end subroutine ice_fallspeed


  subroutine cloud_fallspeed(rhof, w1d, l_qc, rc, nc, ilamc, dz1d, vt, vtn, &
    ktop_sedi)
    !! calculates mass and number weighted fall speeds for cloud
    !! and the top k-level of sedimentation
    use module_mp_tempo_params, only : av_c, ccg, ocg1, ocg2, bv_c, r2

    real(wp), dimension(:), intent(in) :: rhof, w1d, dz1d, rc, nc
    real(dp), dimension(:), intent(in) :: ilamc
    logical, dimension(:), intent(in) :: l_qc
    real(wp), dimension(:), intent(inout) :: vt, vtn
    integer, intent(out) :: ktop_sedi
    real(wp) :: hgt
    integer :: k, nz, nu_c

    nz = size(l_qc)
    ktop_sedi = 1

    ! clouds/fog settle below 500 m agl
    hgt = 0._wp
    hgt_loop : do k = 1, nz-1
      if (rc(k) > r2) ktop_sedi = k
      hgt = hgt + dz1d(k)
      if (hgt > 500._wp) exit hgt_loop
    enddo hgt_loop

    do k = ktop_sedi, 1, -1
      if (rc(k) > r1 .and. w1d(k) < 0.1_wp) then
        nu_c = get_nuc(nc(k))
        vt(k) = rhof(k)*av_c*ccg(5,nu_c)*ocg2(nu_c) * ilamc(k)**bv_c
        vtn(k) = rhof(k)*av_c*ccg(4,nu_c)*ocg1(nu_c) * ilamc(k)**bv_c
      endif 
    enddo
  end subroutine cloud_fallspeed


  subroutine cloud_condensation(rho, temp, w1d, ssatw, lvap, tcond, diffu, lvt2, &
    nwfa, qv, qvs, l_qc, rc, nc, tend, dt, odt)
    !! cloud condensation and evaporation
    use module_mp_tempo_params, only : eps, r1, t0, orv, pi, rho_w, nbc, &
      tnc_wev, nt_c_min

    real(wp), intent(in) :: dt, odt
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(:), intent(in) :: rho, temp, w1d, ssatw, lvap, tcond, diffu, lvt2, &
      nwfa, qv, qvs, rc, nc
    logical, dimension(:), intent(in) :: l_qc
    real(wp) :: clap, fcd, dfcd, xrc, xnc, orho, tempc, otemp, &
      rvs, rvs_p, rvs_pp, gamsc, alphsc, xsat, t1_evap
    real(dp) :: dc_star
    integer :: k, nz, n, idx_d, idx_n, idx_c

    nz = size(qv)
    do k = 1, nz
      if (abs(ssatw(k)) < eps) cycle ! RH = 100%

      orho = 1._wp/rho(k)
      clap = (qv(k)-qvs(k))/(1._wp + lvt2(k)*qvs(k))
      do n = 1, 3
        fcd = qvs(k)*exp(lvt2(k)*clap) - qv(k) + clap
        dfcd = qvs(k)*lvt2(k)*exp(lvt2(k)*clap) + 1._wp
        clap = clap - fcd/dfcd
      enddo
      xrc = rc(k) + clap*rho(k)
      xnc = 0._wp

      if (xrc > r1) then
        ! mass tendency
        tend%prw_vcd(k) = clap*odt

        if (clap > eps) then ! condensation
          xnc = max(nt_c_min, activate_cloud_number(temp(k), w1d(k), nwfa(k)))
          tend%pnc_wcd(k) = 0.5_wp*(xnc-nc(k) + abs(xnc-nc(k)))*odt*orho
        elseif (l_qc(k) .and. ssatw(k) < -1.e-6_wp .and. clap < -eps) then ! evaporation
          tempc = temp(k) - t0
          otemp = 1._wp/temp(k)
          rvs = rho(k)*qvs(k)
          rvs_p = rvs*otemp*(lvap(k)*otemp*orv - 1._wp)
          rvs_pp = rvs * (otemp*(lvap(k)*otemp*orv - 1._wp) * &
            otemp*(lvap(k)*otemp*orv - 1._wp) + &
            (-2._wp*lvap(k)*otemp*otemp*otemp*orv) + otemp*otemp)
          gamsc = lvap(k)*diffu(k)/tcond(k) * rvs_p
          alphsc = 0.5_wp*(gamsc/(1._wp+gamsc))*(gamsc/(1._wp+gamsc)) * &
            rvs_pp/rvs_p * rvs/rvs_p
          alphsc = max(1.e-9_wp, alphsc)
          xsat = ssatw(k)
          if (abs(xsat) < 1.e-9_wp) xsat = 0._wp
          t1_evap = 2._wp*pi*(1.0_wp - alphsc*xsat + 2._wp*alphsc*alphsc*xsat*xsat - &
            5._wp*alphsc*alphsc*alphsc*xsat*xsat*xsat) / (1._wp+gamsc)

          dc_star = sqrt(-2._dp*dt * t1_evap/(2._dp*pi) * &
            4._dp*diffu(k)*ssatw(k)*rvs/rho_w)
          idx_d = max(1, min(int(1.e6_wp*dc_star), nbc))
          call get_cloud_table_index(rc(k), nc(k), idx_c, idx_n)

          tend%prw_vcd(k) = max(real(-rc(k)*0.99_wp*orho*odt, kind=dp), &
            tend%prw_vcd(k))
          tend%pnc_wcd(k) = max(real(-nc(k)*0.99_wp*orho*odt, kind=dp), &
            -tnc_wev(idx_d, idx_c, idx_n)*orho*odt)
        endif 
      else
        tend%prw_vcd(k) = -rc(k)*orho*odt
        tend%pnc_wcd(k) = -nc(k)*orho*odt
      endif 
    enddo 
  end subroutine cloud_condensation


  subroutine rain_evaporation(rho, temp, ssatw, lvap, tcond, diffu, &
    vsc2, rhof2, qv, qvs, l_qr, rr, nr, ilamr, tend, odt)
    !! rain evaporation that includes reduction in the evaporation rate
    !! in the presence of melting graupel
    use module_mp_tempo_params, only : eps, r1, t0, orv, pi, rho_w, &
      org2, cre, t1_qr_ev, t2_qr_ev, fv_r

    real(wp), intent(in) :: odt
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(:), intent(in) :: rho, temp, ssatw, lvap, tcond, diffu, vsc2, rhof2, &
      qv, qvs, rr, nr
    logical, dimension(:), intent(in) :: l_qr
    real(dp), dimension(:), intent(in) :: ilamr
    real(wp) :: orho, tempc, otemp, &
      rvs, rvs_p, rvs_pp, gamsc, alphsc, xsat, t1_evap, rate_max, eva_factor
    real(dp) :: lamr, n0_r
    integer :: k, nz

    nz = size(qv)
    do k = 1, nz
      if(l_qr(k)) then
        if ((ssatw(k) < -eps) .and. tend%prw_vcd(k) <= 0._dp) then
          orho = 1._wp/rho(k)
          tempc = temp(k) - t0
          otemp = 1._wp/temp(k)
          rvs = rho(k)*qvs(k)
          rvs_p = rvs*otemp*(lvap(k)*otemp*orv - 1._wp)
          rvs_pp = rvs * (otemp*(lvap(k)*otemp*orv - 1._wp) * &
            otemp*(lvap(k)*otemp*orv - 1._wp) + &
            (-2._wp*lvap(k)*otemp*otemp*otemp*orv) + otemp*otemp)
          gamsc = lvap(k)*diffu(k)/tcond(k) * rvs_p
          alphsc = 0.5_wp*(gamsc/(1._wp+gamsc))*(gamsc/(1._wp+gamsc)) * &
            rvs_pp/rvs_p * rvs/rvs_p
          alphsc = max(1.e-9_wp, alphsc)
          xsat = min(-1.e-9_wp, ssatw(k))
          t1_evap = 2._wp*pi*(1.0_wp - alphsc*xsat + 2._wp*alphsc*alphsc*xsat*xsat - &
            5._wp*alphsc*alphsc*alphsc*xsat*xsat*xsat) / (1._wp+gamsc)

          !> @note
          !> rain evaporation rapidly eliminates near zero values when low humidity (<95%)
          !> @endnotes
          if (qv(k)/qvs(k) < 0.95_wp .and. rr(k)*orho <= 1.e-8_wp) then
            tend%prv_rev(k) = rr(k)*orho*odt
          else
            lamr = 1._dp/ilamr(k)
            n0_r = nr(k)*org2*lamr**cre(2)
            tend%prv_rev(k) = t1_evap*diffu(k)*(-ssatw(k))*n0_r*rvs * &
              (t1_qr_ev*ilamr(k)**cre(10) + t2_qr_ev*vsc2(k)*rhof2(k)* &
              ((lamr+0.5*fv_r)**(-cre(11))))
            rate_max = min((rr(k)*orho*odt), &
              (qvs(k)-qv(k))*odt)
            tend%prv_rev(k) = min(real(rate_max, kind=dp), tend%prv_rev(k)*orho)
            if (tend%prr_gml(k) > 0._dp) then
              eva_factor = min(1._wp, 0.01_wp+(0.99_wp-0.01_wp)*(tempc/20._wp))
              tend%prv_rev(k) = tend%prv_rev(k)*eva_factor
            endif
          endif
          tend%pnr_rev(k) = min(real(nr(k)*0.99*orho*odt, kind=dp),  &
              tend%prv_rev(k) * nr(k)/rr(k))
         endif  
        endif 
    enddo 
  end subroutine rain_evaporation


  subroutine freeze_cloud_melt_ice(temp, rho, ocp, lvap, qi1d, ni1d, qiten, niten, &
    qc1d, nc1d, qcten, ncten, tten, ncsave, dt, odt)
    ! freezes all cloud water and melts all cloud ice instantly given the temperature
    use module_mp_tempo_params, only : t0, lfus, lsub, hgfrz, nt_c_l

    real(wp), intent(in) :: dt, odt
    real(wp), dimension(:), intent(in) :: temp, rho, ocp, lvap, qi1d, ni1d, qc1d
    real(wp), dimension(:), intent(inout) :: qiten, niten, qcten, ncten, tten
    real(wp), dimension(:), intent(in), optional :: nc1d, ncsave
    real(wp) :: xri, xrc, lfus2, xnc
    integer :: k, nz

    nz = size(temp)
    do k = 1, nz
      ! instantly melt all cloud ice
      xri = max(0._wp, qi1d(k)+qiten(k)*dt)
      if ((temp(k) > t0) .and. (xri > 0._wp)) then
        qcten(k) = qcten(k) + xri*odt
        ncten(k) = ncten(k) + ni1d(k)*odt
        qiten(k) = qiten(k) - xri*odt
        niten(k) = -ni1d(k)*odt
        tten(k) = tten(k) - lfus*ocp(k)*xri*odt
      endif
      ! instantly freeze all cloud water
      xrc = max(0._wp, qc1d(k)+qcten(k)*dt)
      if ((temp(k) < hgfrz) .and. (xrc > 0._wp)) then
        lfus2 = lsub - lvap(k)
        if (present(nc1d)) then
          xnc = nc1d(k) + ncten(k)*dt
        elseif (present(ncsave)) then
          xnc = ncsave(k)/rho(k) + ncten(k)*dt
        else
          xnc = nt_c_l/rho(k) + ncten(k)*dt
        endif
        qiten(k) = qiten(k) + xrc*odt
        niten(k) = niten(k) + xnc*odt
        qcten(k) = qcten(k) - xrc*odt
        ncten(k) = ncten(k) - xnc*odt
        tten(k) = tten(k) + lfus2*ocp(k)*xrc*odt
      endif
    enddo
  end subroutine freeze_cloud_melt_ice


  function koop_nucleation(temp, satw, naero, dt) result(nuc)
    !! aqueous solution freezing of water from 
    !! [Koop et al. (2000)](https://doi.org/10.1038/35020537)
    !! newer research suggests that the freezing rate should be lower 
    !! than original paper, so J_rate is reduced by two orders of magnitude
    use module_mp_tempo_params, only : r_uni, ar_volume

    real(wp), intent(in) :: temp, satw, naero, dt
    real(wp) :: xni, mu_diff, a_w_i, delta_aw, log_j_rate, j_rate, prob_h
    real(wp) :: nuc

    xni = 0.0_wp

    mu_diff = 210368._wp + (131.438_wp*temp) - &
      (3.32373e6_wp/temp) - (41729.1_wp*log(temp))
    a_w_i = exp(mu_diff/(r_uni*temp))
    delta_aw = satw - a_w_i

    log_j_rate = -906.7_wp + (8502._wp*delta_aw) - &
      (26924._wp*delta_aw*delta_aw) + (29180._wp*delta_aw*delta_aw*delta_aw)
    log_j_rate = min(20._wp, log_j_rate)
    j_rate = 10._wp**log_j_rate ! cm-3 s-1
    prob_h = min(1._wp-exp(-j_rate*ar_volume*dt), 1._wp)
    if (prob_h > 0._wp) then
      xni = min(prob_h*naero, 1000.e3_wp)
    endif
    nuc = max(0._wp, xni)
  end function koop_nucleation


  function activate_cloud_number(temp, w1d, nwfa, land) result(activ)
    !! calculations numer of cloud droplets activated
    use module_mp_tempo_params, only : ta_na, ntb_arc, ta_ww, ntb_arw, &
      ta_tk, ntb_art, tnccn_act

    real(wp), intent(in) :: temp, w1d, nwfa
    integer, intent(in), optional :: land
    real(wp) :: n_local, w_local
    real(wp) :: a, b, c, d, t, u, x1, x2, y1, y2, nx, wy, fraction
    real(wp) :: lower_lim_nuc_frac
    integer :: i, j, k, l, m, n
    real(wp) :: activ

    ! index for number of aerosols
    n_local = nwfa * 1.e-6_wp
    if (n_local >= ta_na(ntb_arc)) then
      n_local = ta_na(ntb_arc) - 1.0_wp
    elseif (n_local <= ta_na(1)) then
      n_local = ta_na(1) + 1.0_wp
    endif
    nindex: do n = 2, ntb_arc
      if (n_local >= ta_na(n-1) .and. n_local < ta_na(n)) exit nindex
    enddo nindex
    i = n
    x1 = log(ta_na(i-1))
    x2 = log(ta_na(i))

    ! index for vertical velocity
    w_local = w1d
    if (w_local >= ta_ww(ntb_arw)) then
        w_local = ta_ww(ntb_arw) - 1.0_wp
    elseif (w_local <= ta_ww(1)) then
        w_local = ta_ww(1) + 0.001_wp
    endif
    windex: do n = 2, ntb_arw
      if (w_local >= ta_ww(n-1) .and. w_local < ta_ww(n)) exit windex
    enddo windex
    j = n
    y1 = log(ta_ww(j-1))
    y2 = log(ta_ww(j))

    k = max(1, min(nint((temp - ta_tk(1))*0.1_wp) + 1, ntb_art))

    ! the next two values are indexes of mean aerosol radius and
    ! hygroscopicity and are currently constant 
    !> @todo
    !> separation tiny size sulfates from larger sea salts
    !> @endtodo
    l = 3
    m = 2

    !> @note
    !> there is a lower limit set for activation over water to improve cloud coverage
    !> @endnote
    lower_lim_nuc_frac = 0.
    if (present(land)) then
      if (land == 1) then ! land
        lower_lim_nuc_frac = 0.
      elseif (land == 0) then ! not land (water/ice)
        lower_lim_nuc_frac = 0.15
      else
        lower_lim_nuc_frac = 0.15 ! catch-all for anything else
      endif
    endif
    
    a = tnccn_act(i-1,j-1,k,l,m)
    b = tnccn_act(i,j-1,k,l,m)
    c = tnccn_act(i,j,k,l,m)
    d = tnccn_act(i-1,j,k,l,m)
    nx = log(n_local)
    wy = log(w_local)
    t = (nx-x1)/(x2-x1)
    u = (wy-y1)/(y2-y1)

    fraction = (1.0_wp-t)*(1.0_wp-u)*a + t*(1.0_wp-u)*b + t*u*c + (1.0_wp-t)*u*d
    fraction = max(fraction, lower_lim_nuc_frac)

    activ = nwfa*fraction
  end function activate_cloud_number


  subroutine warm_rain(rhof, l_qc, rc, nc, ilamc, mvd_c, l_qr, rr, nr, mvd_r, tend, odt)
    !! computes warm-rain process rates -- condensation/evaporation happen later
    use module_mp_tempo_params, only : d0r, d0c, r1, nbr, t_efrw, &
      t1_qr_qc, mu_r, am_r, ccg, obmr, ocg2, dr, org2, cre, fv_r, &
      autocon_nr_factor

    real(wp), intent(in) :: odt
    real(wp), dimension(:), intent(in) :: rhof, mvd_r, mvd_c, rr, nr, rc, nc
    real(dp), dimension(:), intent(in) :: ilamc
    logical, dimension(:), intent(in) :: l_qc, l_qr
    type(ty_tend), intent(inout) :: tend

    real(dp) :: lamr, lamc, n0_r
    real(wp) :: ef_rr, dc_g, dc_b, xdc, zeta1, zeta, taud, tau, ef_rw
    integer :: k, nz, nu_c, idx

    nz = size(l_qc)
    !> @note
    !> rain self-collection is from
    !> [Seifert and Beheng (2001)](https://doi.org/10.1016/S0169-8095(01)00126-0)
    !> and drop break-up follows
    !> [Verlinde and Cotton (1993)](https://doi.org/10.1175/1520-0493(1993)121<2776:FMOONC>2.0.CO;2)
    do k = 1, nz
      if (l_qr(k)) then
        if (mvd_r(k) > d0r) then
          ef_rr = max(-0.1_wp, 1.0_wp - exp(2300.0_wp*(mvd_r(k)-1950.0e-6_wp)))
          tend%pnr_rcr(k) = ef_rr * 2.0_wp*nr(k)*rr(k)
        endif
      endif

      !>
      !> autoconversion follows 
      !> [Berry and Reinhardt (1974)](https://doi.org/10.1175/1520-0469(1974)031<1814:AAOCDG>2.0.CO;2)
      !> with characteristic diameters correctly computed from gamma distribution of cloud droplets
      if (l_qc(k)) then
        if (rc(k) > 0.01e-3_wp) then
          nu_c = get_nuc(nc(k))
          lamc = 1._dp / ilamc(k)       
          xdc = max(d0c*1.e6_wp, ((rc(k)/(am_r*nc(k)))**obmr) * 1.e6_wp)
          dc_g = ((ccg(3,nu_c)*ocg2(nu_c))**obmr / lamc) * 1.e6_wp
          dc_b = (xdc*xdc*xdc*dc_g*dc_g*dc_g - xdc*xdc*xdc*xdc*xdc*xdc) &
              **(1._wp/6._wp)
          zeta1 = 0.5_wp*((6.25e-6_wp*xdc*dc_b*dc_b*dc_b - 0.4_wp) &
              + abs(6.25e-6_wp*xdc*dc_b*dc_b*dc_b - 0.4_wp))
          zeta = 0.027_wp*rc(k)*zeta1
          taud = 0.5_wp*((0.5_wp*dc_b - 7.5_wp) + abs(0.5_wp*dc_b - 7.5_wp)) + r1
          tau = 3.72_wp/(rc(k)*taud)
          tend%prr_wau(k) = zeta/tau
          tend%prr_wau(k) = min(real(rc(k)*odt, kind=dp), &
            tend%prr_wau(k))
          tend%pnr_wau(k) = tend%prr_wau(k) / (am_r*nu_c*autocon_nr_factor*d0r*d0r*d0r) 
          tend%pnc_wau(k) = min(real(nc(k)*odt, kind=dp), &
            tend%prr_wau(k) / (am_r*mvd_c(k)*mvd_c(k)*mvd_c(k)))
        endif
      endif

      !>
      !> rain collecting cloud water - assumes dc << dr and vtc \(\approx 0\)
      !> @endnote
      if (l_qr(k) .and. l_qc(k)) then  
        if (mvd_r(k) > d0r .and. mvd_c(k) > d0c) then
          lamr = (3.0_dp + mu_r + 0.672_dp) / mvd_r(k)
          idx = 1 + int(nbr*log(real(mvd_r(k)/dr(1), kind=dp)) / &
            log(real(dr(nbr)/dr(1), kind=dp)))
          idx = min(idx, nbr)
          ef_rw = t_efrw(idx, int(mvd_c(k)*1.e6_wp))
          n0_r = nr(k)*org2*lamr**cre(2)
          tend%prr_rcw(k) = rhof(k)*t1_qr_qc*ef_rw*rc(k)*n0_r * &
            ((lamr+fv_r)**(-cre(9)))
          tend%prr_rcw(k) = min(real(rc(k)*odt, kind=dp), tend%prr_rcw(k))
          tend%pnc_rcw(k) = rhof(k)*t1_qr_qc*ef_rw*nc(k)*n0_r * &
            ((lamr+fv_r)**(-cre(9)))
          tend%pnc_rcw(k) = min(real(nc(k)*odt, kind=dp), tend%pnc_rcw(k))
        endif
      endif 
    enddo
  end subroutine warm_rain


  subroutine riming(temp, rhof, visco, l_qc, rc, nc, ilamc, mvd_c, &
    l_qs, rs, smo0, smob, smoc, smoe, vtboost, l_qg, rg, ng, ilamg, idx, tend, odt)
    !! snow and graupel riming
    use module_mp_tempo_params, only : d0c, d0s, nbs, ds, t_efsw, t1_qs_qc, &
      r_g, bm_g, mu_g, av_g, cgg, ogg3, bv_g, rho_w, t0, d0g, pi, cge, ogg2, &
      rime_threshold, rime_conversion, av_s, bv_s, rho_s, xm0i, eps, fv_s, &
      meters3_to_liters

    real(wp), intent(in) :: odt    
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(:), intent(in) :: rhof, visco, temp, rc, nc, rs, rg, ng
    real(wp), dimension(:), intent(in) :: mvd_c
    real(dp), dimension(:), intent(in) :: smo0, smob, smoc, smoe, ilamg, ilamc
    logical, dimension(:), intent(in) :: l_qc, l_qs, l_qg
    integer, dimension(:), intent(in) :: idx
    real(wp), dimension(:), intent(out) :: vtboost

    real(dp) :: xds, xdg, n0_g, lamc
    real(wp) :: ef_sw, vtg, stoke_g, const_ri, tempc, rime_dens, ef_gw
    real(wp) :: t1_qg_qc, r_frac, g_frac, vts, tf, snow_dens_frac
    integer :: k, nz, idxs, nu_c

    nz = size(l_qc)
    do k = 1, nz
      tempc = temp(k) - t0
      vtboost(k) = 1._wp
      if (l_qc(k) .and. l_qs(k)) then
        nu_c = get_nuc(nc(k))
        lamc = 1._dp / ilamc(k)
        xds = smoc(k) / smob(k)
        if ((mvd_c(k) > d0c) .and. (xds > d0s)) then
          !> @note
          !> snow collecting cloud water - assume dc << ds and vtc \(\approx 0\)
          idxs = 1 + int(nbs*log(real(xds/ds(1), kind=dp)) / log(real(ds(nbs)/ds(1), kind=dp)))
          idxs = min(idxs, nbs)
          ef_sw = t_efsw(idxs, int(mvd_c(k)*1.e6_wp))
          tend%prs_scw(k) = rhof(k)*t1_qs_qc*ef_sw*rc(k)*smoe(k)
          tend%prs_scw(k) = min(real(rc(k)*odt, kind=dp), tend%prs_scw(k))
          tend%pnc_scw(k) = rhof(k)*t1_qs_qc*ef_sw*nc(k)*smoe(k)
          tend%pnc_scw(k) = min(real(nc(k)*odt, kind=dp), tend%pnc_scw(k))

          !>
          !> at temperatures below melting, if the riming rate is greater than the depositional
          !> growth rate for snow by a factor rime_threshold, convert a portion of rimed snow
          !> to graupel
          if (temp(k) < t0) then
            if (tend%prs_scw(k) > rime_threshold*tend%prs_sde(k) .and. &
              tend%prs_sde(k) > eps) then
              r_frac = min(30.0_dp, tend%prs_scw(k)/tend%prs_sde(k))
              g_frac = min(rime_conversion, 0.15_wp + (r_frac-2._wp)*.028_wp)
              vtboost(k) = min(1.5_wp, 1.1_wp + (r_frac-2.)*.014_wp)
              tend%prg_scw(k) = g_frac*tend%prs_scw(k)
              tend%png_scw(k) = tend%prg_scw(k)*smo0(k)/rs(k)
              vts = av_s*xds**bv_s * exp(-fv_s*xds)
              const_ri = -1._wp*(mvd_c(k)*0.5e6_wp)*vts/min(-0.1_wp,tempc)
              const_ri = max(0.1_wp, min(const_ri, 10._wp))
              rime_dens = (0.051_wp + 0.114_wp*const_ri - 0.0055_wp*const_ri*const_ri)*1000._wp
              if(rime_dens < 150._wp) then
                g_frac = 0._wp
                tend%prg_scw(k) = 0._dp
                tend%png_scw(k) = 0._dp
              endif
              snow_dens_frac = min(1._wp, max(0._wp, rs(k)*odt / &
                (rs(k)*odt + tend%prg_scw(k))))
              tend%pbg_scw(k) = meters3_to_liters*tend%prg_scw(k) / &
                (rho_s * snow_dens_frac + rime_dens * (1._wp-snow_dens_frac))
              ! tend%pbg_scw(k) = meters3_to_liters*tend%prg_scw(k) / &
              !  (0.5_wp*(rho_s+rime_dens))
              tend%prs_scw(k) = (1._wp - g_frac)*tend%prs_scw(k)
            endif
          endif 
        endif
      endif 

      if (l_qc(k) .and. l_qg(k)) then
        !>
        !> graupel collecting cloud water - assume dc << dg and vtc \(\approx 0\)
        if (rg(k) >= r_g(1) .and. mvd_c(k) > d0c) then
          xdg = (bm_g + mu_g + 1._wp) * ilamg(k)
          vtg = rhof(k)*av_g(idx(k))*cgg(6,idx(k))*ogg3 * ilamg(k)**bv_g(idx(k))
          stoke_g = mvd_c(k)*mvd_c(k)*vtg*rho_w/(9._wp*visco(k)*xdg)
          !>
          !> rime density formula is from 
          !> [Cober and List (1993)](https://doi.org/10.1175/1520-0469(1993)050<1591:MOTHAM>2.0.CO;2)
          const_ri = -1._wp*(mvd_c(k)*0.5e6_wp)*vtg/min(-0.1_wp, tempc)
          const_ri = max(0.1_wp, min(const_ri, 10._wp))
          rime_dens = (0.051_wp + 0.114_wp*const_ri - 0.0055_wp*const_ri*const_ri)*1000._wp
          if (xdg > d0g) then
            if (stoke_g >= 0.4_wp .and. stoke_g <= 10._wp) then
              ef_gw = 0.55_wp*log10(2.51_wp*stoke_g)
            elseif (stoke_g < 0.4_wp) then
              ef_gw = 0.0_wp
            elseif (stoke_g > 10._wp) then
              ef_gw = 0.77_wp
            endif
            !>
            !> hail size increases below the melting level so the collection efficiency
            !> is reduced (proxy for shedding of collected cloud water)
            if (temp(k) > t0) ef_gw = ef_gw*0.1_wp
            t1_qg_qc = pi*.25_wp*av_g(idx(k)) * cgg(9,idx(k))
            n0_g = ng(k)*ogg2*(1._wp/ilamg(k))**cge(2,1)
            tend%prg_gcw(k) = rhof(k)*t1_qg_qc*ef_gw*rc(k)* &
              n0_g*ilamg(k)**cge(9,idx(k))
            tend%pnc_gcw(k) = rhof(k)*t1_qg_qc*ef_gw*nc(k)* &
              n0_g*ilamg(k)**cge(9,idx(k))
            tend%pnc_gcw(k) = min(real(nc(k)*odt, kind=dp), tend%pnc_gcw(k))
            if (temp(k) < t0) tend%pbg_gcw(k) = meters3_to_liters*tend%prg_gcw(k)/rime_dens

            if (temp(k) < t0) then
              !>
              !> rime splintering is from
              !> [Hallet and Mossop (1974)](https://doi.org/10.1038/249026a0)
              !> @endnote
              if (tend%prg_gcw(k) > eps .and. tempc > -8._wp) then
                tf = 0._wp
                if (tempc >= -5._wp .and. tempc < -3._wp) then
                  tf = 0.5_wp*(-3.0_wp - tempc)
                elseif (tempc > -8._wp .and. tempc < -5._wp) then
                  tf = 0.33333333_wp*(8._wp + tempc)
                endif
                tend%pni_ihm(k) = 3.5e8_wp*tf*tend%prg_gcw(k)
                tend%pri_ihm(k) = xm0i*tend%pni_ihm(k)
                tend%prs_ihm(k) = tend%prs_scw(k)/(tend%prs_scw(k)+tend%prg_gcw(k)) * &
                  tend%pri_ihm(k)
                tend%prg_ihm(k) = tend%prg_gcw(k)/(tend%prs_scw(k)+tend%prg_gcw(k)) * &
                  tend%pri_ihm(k)
              endif
            endif
          endif 
        endif
      endif 
    end do 
  end subroutine riming


  subroutine get_snow_table_index(rs, idx_s)
    !! get snow table index from snow mass
    use module_mp_tempo_params, only : ntb_s, nis2

    real(wp), intent(in) :: rs
    integer :: nis, nn, n
    integer, intent(out) :: idx_s

    nis = nint(log10(rs))
    do_loop_rs: do nn = nis-1, nis+1
      n = nn
      if ((rs/10._wp**nn) >= 1._wp .and. (rs/10._wp**nn) < 10._wp) exit do_loop_rs
    enddo do_loop_rs
    idx_s = int(rs/10._wp**n) + 10*(n-nis2) - (n-nis2)
    idx_s = max(1, min(idx_s, ntb_s))
  end subroutine get_snow_table_index


  subroutine get_temperature_table_index(tempk, idx_t)
    !! get temperature table index
    use module_mp_tempo_params, only : t0, ntb_t
    
    real(wp), intent(in) :: tempk
    real(wp) :: tempc
    integer, intent(out) :: idx_t

    tempc = tempk - t0
    idx_t = int((tempc-2.5_wp)/5._wp) - 1
    idx_t = max(1, -idx_t)
    idx_t = min(idx_t, ntb_t)
  end subroutine get_temperature_table_index


  subroutine get_rain_table_index(rr, ilamr, idx_r, idx_r1)
    !! get rain table indices from rain mass and lambda
    use module_mp_tempo_params, only : nir2, nir3, ntb_r, ntb_r1, &
      org2, org1, bm_r, am_r, crg, cre

    real(wp), intent(in) :: rr
    real(dp), intent(in) :: ilamr
    real(dp) :: lamr, lam_exp, n0_exp
    integer :: nir, nn, n
    integer, intent(out) :: idx_r, idx_r1

    nir = nint(log10(rr))
    do_loop_rr: do nn = nir-1, nir+1
      n = nn
      if ((rr/10._wp**nn) >= 1._wp .and. (rr/10._wp**nn) < 10._wp) exit do_loop_rr
    enddo do_loop_rr
    idx_r = int(rr/10._wp**n) + 10*(n-nir2) - (n-nir2)
    idx_r = max(1, min(idx_r, ntb_r))

    lamr = 1./ilamr
    lam_exp = lamr * (crg(3)*org2*org1)**bm_r
    n0_exp = org1*rr/am_r * lam_exp**cre(1)
    nir = nint(log10(real(n0_exp, kind=dp)))
    do_loop_nr: do nn = nir-1, nir+1
      n = nn
      if ((n0_exp/10._wp**nn) >= 1._wp .and. (n0_exp/10._wp**nn) < 10._wp) exit do_loop_nr
    enddo do_loop_nr
    idx_r1 = int(n0_exp/10._wp**n) + 10*(n-nir3) - (n-nir3)
    idx_r1 = max(1, min(idx_r1, ntb_r1))
  end subroutine get_rain_table_index


  subroutine get_graupel_table_index(rg, ilamg, idx, idx_g, idx_g1)
    !! get graupel table indices from graupel mass, lambda, and density index
    use module_mp_tempo_params, only : nig2, ntb_g, ntb_g1, ogg2, ogg1, &
      bm_g, cgg, ogg1, am_g, cge, nig3

    real(wp), intent(in) :: rg
    real(dp), intent(in) :: ilamg
    integer, intent(in) :: idx
    real(dp) :: lamg, lam_exp, n0_exp
    integer :: nig, nn, n
    integer, intent(out) :: idx_g, idx_g1

    nig = nint(log10(rg))
    do_loop_rg: do nn = nig-1, nig+1
      n = nn
      if ( (rg/10._wp**nn) >= 1._wp .and. (rg/10._wp**nn).lt.10._wp) exit do_loop_rg
    enddo do_loop_rg
    idx_g = int(rg/10._wp**n) + 10*(n-nig2) - (n-nig2)
    idx_g = max(1, min(idx_g, ntb_g))

    lamg = 1./ilamg
    lam_exp = lamg * (cgg(3,1)*ogg2*ogg1)**bm_g
    n0_exp = ogg1*rg/am_g(idx) * lam_exp**cge(1,1)
    nig = nint(log10(real(n0_exp, kind=dp)))
    do_loop_ng: do nn = nig-1, nig+1
      n = nn
      if ( (n0_exp/10._wp**nn) >= 1._wp .and. (n0_exp/10._wp**nn) < 10._wp) exit do_loop_ng
    enddo do_loop_ng
    idx_g1 = int(n0_exp/10._wp**n) + 10*(n-nig3) - (n-nig3)
    idx_g1 = max(1, min(idx_g1, ntb_g1))
  end subroutine get_graupel_table_index


  subroutine get_cloud_table_index(rc, nc, idx_c, idx_n)
    !! get cloud table index from mass and number
    use module_mp_tempo_params, only : nbc, ntb_c, r_c, nic2, t_nc, nic1

    real(wp), intent(in) :: rc, nc
    integer, intent(out) :: idx_c, idx_n
    integer :: nic, nn, n

    nic = nint(log10(rc))
    do_loop_rc: do nn = nic-1, nic+1
      n = nn
      if ( (rc/10._wp**nn) >= 1._wp .and. (rc/10._wp**nn) < 10._wp) exit do_loop_rc
    enddo do_loop_rc
    idx_c = int(rc/10._wp**n) + 10*(n-nic2) - (n-nic2)
    idx_c = max(1, min(idx_c, ntb_c))
          
    idx_n = nint(1._wp + real(nbc, kind=wp) * log(real(nc/t_nc(1), kind=dp)) / nic1)
    idx_n = max(1, min(idx_n, nbc))
  end subroutine get_cloud_table_index


  subroutine get_ice_table_index(ri, ni, idx_i, idx_i1)
    !! get ice table index from mass and number
    use module_mp_tempo_params, only : ntb_i, ntb_i1, nii2, nii3

    real(wp), intent(in) :: ri, ni
    integer, intent(out) :: idx_i, idx_i1
    integer :: nii, nn, n

    nii = nint(log10(ri))
    do_loop_ri: do nn = nii-1, nii+1
      n = nn
      if ( (ri/10._wp**nn) >= 1._wp .and. (ri/10._wp**nn) < 10._wp) exit do_loop_ri
    enddo do_loop_ri
    idx_i = int(ri/10._wp**n) + 10*(n-nii2) - (n-nii2)
    idx_i = max(1, min(idx_i, ntb_i))
  
    nii = nint(log10(ni))
    do_loop_ni: do nn = nii-1, nii+1
      n = nn
      if ( (ni/10._wp**nn) >= 1._wp .and. (ni/10._wp**nn) < 10._wp) exit do_loop_ni
    enddo do_loop_ni
    idx_i1 = int(ni/10._wp**n) + 10*(n-nii3) - (n-nii3)
    idx_i1 = max(1, min(idx_i1, ntb_i1))
  end subroutine get_ice_table_index


  subroutine rain_snow_rain_graupel(temp, l_qr, rr, nr, ilamr, l_qs, rs, &
      l_qg, rg, ng, ilamg, idx, tend, odt)
    !! calculates rain-snow and rain-graupel collection
    use module_mp_tempo_params, only : t0, r_r, r_s, r_g, rho_i, rho_g, meters3_to_liters, &
      tmr_racs2, tcr_sacr2, tmr_racs1, tcr_sacr1, tms_sacr1, tcs_racs1, &
      tnr_sacr1, tnr_sacr2, tnr_racs1, tnr_racs2, &
      tcr_gacr, tmr_racg, tcg_racg, tnr_gacr, tnr_racg

    real(wp), intent(in) :: odt    
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(:), intent(in) :: temp, rr, rs, rg, nr, ng
    real(dp), dimension(:), intent(in) :: ilamr, ilamg
    logical, dimension(:), intent(in) :: l_qr, l_qs, l_qg
    integer, dimension(:), intent(in) :: idx
    integer :: k, nz, idx_r, idx_r1, idx_s, idx_t, &
      idx_g, idx_g1

    nz = size(l_qg)
    do k = 1, nz
      if (l_qr(k) .and. l_qs(k)) then
        if (rr(k) >= r_r(1) .and. rs(k) >= r_s(1)) then
          call get_temperature_table_index(temp(k), idx_t)
          call get_rain_table_index(rr(k), ilamr(k), idx_r, idx_r1)
          call get_snow_table_index(rs(k), idx_s)
          if (temp(k) < t0) then
            tend%prr_rcs(k) = -(tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
              + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
              + tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
              + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r))
            tend%prs_rcs(k) = tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
              + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
              - tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
              - tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
            tend%prg_rcs(k) = tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
              + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r) &
              + tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
              + tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
            tend%prr_rcs(k) = max(real(-rr(k)*odt, kind=dp), tend%prr_rcs(k))
            tend%prs_rcs(k) = max(real(-rs(k)*odt, kind=dp), tend%prs_rcs(k))
            tend%prg_rcs(k) = min(real((rr(k)+rs(k))*odt, kind=dp), &
              tend%prg_rcs(k))
            tend%pnr_rcs(k) = tnr_racs1(idx_s,idx_t,idx_r1,idx_r) &
              + tnr_racs2(idx_s,idx_t,idx_r1,idx_r) &
              + tnr_sacr1(idx_s,idx_t,idx_r1,idx_r) &
              + tnr_sacr2(idx_s,idx_t,idx_r1,idx_r)
            tend%pnr_rcs(k) = min(real(nr(k)*odt, kind=dp), tend%pnr_rcs(k))
            tend%png_rcs(k) = tend%pnr_rcs(k)
            tend%pbg_rcs(k) = meters3_to_liters*tend%prg_rcs(k)/rho_i
          else
            tend%prs_rcs(k) = -tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
              - tms_sacr1(idx_s,idx_t,idx_r1,idx_r) &
              + tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
              + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r)
            tend%prs_rcs(k) = max(real(-rs(k)*odt, kind=dp), tend%prs_rcs(k))
            tend%prr_rcs(k) = -tend%prs_rcs(k)
          endif
        endif
      endif 

      if (l_qr(k) .and. l_qg(k)) then
        if (rr(k) >= r_r(1) .and. rg(k) >= r_g(1)) then
          call get_temperature_table_index(temp(k), idx_t)        
          call get_rain_table_index(rr(k), ilamr(k), idx_r, idx_r1)
          call get_graupel_table_index(rg(k), ilamg(k), idx(k), idx_g, idx_g1)
          if (temp(k) < t0) then
            tend%prg_rcg(k) = tmr_racg(idx_g1,idx_g,idx(k),idx_r1,idx_r) &
              + tcr_gacr(idx_g1,idx_g,idx(k),idx_r1,idx_r)
            tend%prg_rcg(k) = min(real(rr(k)*odt, kind=dp), tend%prg_rcg(k))
            tend%prr_rcg(k) = -tend%prg_rcg(k)
            tend%pnr_rcg(k) = tnr_racg(idx_g1,idx_g,idx(k),idx_r1,idx_r) &
              + tnr_gacr(idx_g1,idx_g,idx(k),idx_r1,idx_r)
            tend%pnr_rcg(k) = min(real(nr(k)*odt, kind=dp), tend%pnr_rcg(k))
            tend%pbg_rcg(k) = meters3_to_liters*tend%prg_rcg(k)/rho_i
          else
            tend%prr_rcg(k) = tcg_racg(idx_g1,idx_g,idx(k),idx_r1,idx_r)
            tend%prr_rcg(k) = min(real(rg(k)*odt, kind=dp), tend%prr_rcg(k))
            tend%prg_rcg(k) = -tend%prr_rcg(k)
            tend%png_rcg(k) = tnr_racg(idx_g1,idx_g,idx(k),idx_r1,idx_r)
            tend%png_rcg(k) = min(real(ng(k)*odt, kind=dp), tend%png_rcg(k))
            tend%pbg_rcg(k) = meters3_to_liters*tend%prg_rcg(k)/rho_g(idx(k))
            !> @note
            !> adds explicit rain drop break-up due to collisions with graupel
            !> at temperatures above melting
            !> @endnote
            tend%pnr_rcg(k) = -1.5_wp*tnr_gacr(idx_g1,idx_g,idx(k),idx_r1,idx_r)
          endif
        endif
      endif 
    enddo 
  end subroutine rain_snow_rain_graupel


  subroutine ice_nucleation(temp, rho, w1d, qv, qvsi, ssati, ssatw, &
      nifa, nwfa, ni, smo0, rc, nc, rr, nr, ilamr, tend, dt, odt)
    !! ice nulceation
    use module_mp_tempo_params, only : r_r, r_c, hgfrz, rho_i, xm0i, &
      tpg_qrfz, tpi_qrfz, tni_qrfz, tnr_qrfz, tpi_qcfz, tni_qcfz, &
      demott_nuc_ssati, eps, icenuc_max, tno, ato, max_ni, meters3_to_liters

    real(wp), intent(in) :: dt, odt
    type(ty_tend), intent(inout) :: tend
    real(wp), dimension(:), intent(in) :: qv, temp, rho, qvsi, rr, nr, rc, nc, w1d, &
      ssati, ssatw, ni
    real(dp), dimension(:), intent(in) :: ilamr, smo0
    real(wp), dimension(:), intent(in), optional :: nifa, nwfa
    real(wp) :: rate_max, tempc, xni, xnc
    integer :: k, nz, idx_in, idx_r, idx_r1, idx_tc, idx_c, idx_n

    nz = size(qv)
    do k = 1, nz
      if (temp(k) < t0) then
        tempc = temp(k) - t0
        idx_tc = max(1, min(nint(-tempc), 45))
        rate_max = (qv(k)-qvsi(k))*rho(k)*odt*0.999_wp
        if (present(nifa)) then
          xni = demott_nucleation(tempc, rho(k), nifa(k))
        else  
          xni = 1._wp * 1000._wp ! 1 / Liter
        endif 
        call get_in_table_index(xni, idx_in)

        !> @note
        !> freezing of water drops into either cloud ice or graupel is from 
        !> [Bigg (1953)](https://doi.org/10.1002/qj.49707934207)
        if (rr(k) > r_r(1)) then
          call get_rain_table_index(rr(k), ilamr(k), idx_r, idx_r1)
          tend%prg_rfz(k) = tpg_qrfz(idx_r,idx_r1,idx_tc,idx_in)*odt
          tend%pri_rfz(k) = tpi_qrfz(idx_r,idx_r1,idx_tc,idx_in)*odt
          tend%pni_rfz(k) = tni_qrfz(idx_r,idx_r1,idx_tc,idx_in)*odt
          tend%pnr_rfz(k) = tnr_qrfz(idx_r,idx_r1,idx_tc,idx_in)*odt
          tend%prg_rfz(k) = min(real(rr(k)*odt, kind=dp), tend%prg_rfz(k))
          tend%pnr_rfz(k) = min(real(nr(k)*odt, kind=dp), tend%pnr_rfz(k))
          ! reduce number of graupel particles created at higher vertical velocities
          tend%png_rfz(k) = tend%pnr_rfz(k) * &
            max(min((10._wp**(-0.1_wp*w1d(k)) + 0.1_wp), 1._wp), 0.1_wp)
          ! tend%png_rfz(k) = tend%pnr_rfz(k)
        elseif (rr(k) > r1 .and. temp(k) < hgfrz) then
          tend%pri_rfz(k) = rr(k)*odt
          tend%pni_rfz(k) = nr(k)*odt
        endif
        tend%pbg_rfz(k) = meters3_to_liters*tend%prg_rfz(k)/rho_i

        if (rc(k) > r_c(1)) then
          call get_cloud_table_index(rc(k), nc(k), idx_c, idx_n)
          tend%pri_wfz(k) = tpi_qcfz(idx_c,idx_n,idx_tc,idx_in)*odt
          tend%pri_wfz(k) = min(real(rc(k)*odt, kind=dp), tend%pri_wfz(k))
          tend%pni_wfz(k) = tni_qcfz(idx_c,idx_n,idx_tc,idx_in)*odt
          tend%pni_wfz(k) = min(real(nc(k)*odt, kind=dp), &
            tend%pri_wfz(k)/(2.0_dp*xm0i), tend%pni_wfz(k))
        elseif (rc(k) > r1 .and. temp(k) < hgfrz) then
          tend%pri_wfz(k) = rc(k)*odt
          tend%pni_wfz(k) = nc(k)*odt
        endif
        !>
        !> deposition nucleation from dust is from
        !> [DeMott et al. (2010)](https://doi.org/10.1073/pnas.0910818107)
        if ((ssati(k) >= demott_nuc_ssati) .or. (ssatw(k) > eps &
            .and. tempc < -20._wp)) then
          if (present(nifa)) then
            xnc = demott_nucleation(tempc, rho(k), nifa(k))
          else
            xnc = min(icenuc_max, tno*exp(ato*(t0-temp(k))))
          endif
          xni = ni(k) + (tend%pni_rfz(k)+tend%pni_wfz(k))*dt
          tend%pni_inu(k) = 0.5_wp*(xnc-xni + abs(xnc-xni))*odt
          tend%pri_inu(k) = min(real(rate_max, kind=dp), xm0i*tend%pni_inu(k))
          tend%pni_inu(k) = tend%pri_inu(k)/xm0i
        endif
        !>
        !> freezing of aqueous aerosols is based on [Koop et al. (2000)](https://doi.org/10.1038/35020537)
        xni = smo0(k)+ni(k) + (tend%pni_rfz(k)+tend%pni_wfz(k)+tend%pni_inu(k))*dt
        if (present(nwfa)) then
          if ((xni <= max_ni) .and.(temp(k) < 238._wp) .and. (ssati(k) >= 0.4_wp)) then
            xnc = koop_nucleation(temp(k), ssatw(k), nwfa(k), dt)
            tend%pni_iha(k) = xnc*odt
            tend%pri_iha(k) = min(real(rate_max, kind=dp), xm0i*0.1_wp*tend%pni_iha(k))
            tend%pni_iha(k) = tend%pri_iha(k)/(xm0i*0.1_wp)
          endif
        endif 
      endif 
    enddo 
  end subroutine ice_nucleation


  function demott_nucleation(tempc, rho, nifa) result(nuc)
    !! DeMott nucleation
    use module_mp_tempo_params, only : rho_not0

    real(wp), intent(in) :: tempc, rho, nifa
    real(wp) :: xni, nifa_cc
    real(wp) :: nuc

    xni = 0._wp
    nifa_cc = max(0.5_wp, nifa*rho_not0*1.e-6_wp/rho)
    xni = (5.94e-5_wp*(-tempc)**3.33_wp) * (nifa_cc**((-0.0264_wp*(tempc))+0.0033_wp))
    xni = xni*rho/rho_not0 * 1000._wp
    nuc = max(0._wp, xni)
  end function demott_nucleation


  subroutine get_in_table_index(xni, idx_in)
    !! get ice nuclei table index
    use module_mp_tempo_params, only : ntb_in, nt_in, niin2

    real(wp), intent(in) :: xni
    integer, intent(out) :: idx_in
    integer :: niin, nn, n

    if (xni >  nt_in(1)) then
        niin = nint(log10(xni))
        do_loop_xni: do nn = niin-1, niin+1
          n = nn
          if ( (xni/10._wp**nn) >= 1._wp .and. (xni/10._wp**nn) < 10._wp) exit do_loop_xni
        enddo do_loop_xni
        idx_in = int(xni/10._wp**n) + 10*(n-niin2) - (n-niin2)
        idx_in = max(1, min(idx_in, ntb_in))
    else
        idx_in = 1
    endif
  end subroutine get_in_table_index


  subroutine get_t1_subl(rho, temp, qvsi, tcond, diffu, ssati, t1_subl)
    !! calculations thermodynamic term used in depositional growth and melting
    use module_mp_tempo_params, only : t0, bm_i, mu_i, lsub, orv, &
      pi, c_sqrd, c_cube, d0s, ntb_i, r_s, ef_si, r_r, fv_r, ef_ri, &
      rho_i, eps, rho_w, rho_g

    real(wp), dimension(:), intent(in) :: rho, temp, qvsi, tcond, diffu, ssati
    real(wp) :: otemp, rvs, rvs_p, rvs_pp, gamsc, alphsc, xsat
    real(wp), dimension(:), intent(out) :: t1_subl
    integer :: k, nz

    nz = size(rho)
    do k = 1, nz
      otemp = 1._wp/temp(k)
      rvs = rho(k)*qvsi(k)
      rvs_p = rvs*otemp*(lsub*otemp*orv - 1._wp)
      rvs_pp = rvs * (otemp*(lsub*otemp*orv - 1._wp) * otemp*(lsub*otemp*orv - 1._wp) + &
        (-2.*lsub*otemp*otemp*otemp*orv) + otemp*otemp)
      gamsc = lsub*diffu(k)/tcond(k) * rvs_p
      alphsc = 0.5_wp*(gamsc/(1._wp+gamsc))*(gamsc/(1._wp+gamsc)) * &
        rvs_pp/rvs_p * rvs/rvs_p
      alphsc = max(1.e-9_wp, alphsc)
      xsat = ssati(k)
      if (abs(xsat) < 1.e-9_wp) xsat = 0._wp
      t1_subl(k) = 4._wp*pi*(1._wp - alphsc*xsat + 2._wp*alphsc*alphsc*xsat*xsat - &
        5._wp*alphsc*alphsc*alphsc*xsat*xsat*xsat) / (1._wp+gamsc)
    end do  
  end subroutine get_t1_subl


  subroutine ice_processes(rhof, rhof2, rho, w1d, temp, qv, qvsi, tcond, diffu, &
    vsc2, ssati, l_qi, ri, ni, ilami, l_qs, rs, smoe, smof, smo1, rr, nr, ilamr, &
    mvd_r, l_qg, rg, ng, ilamg, idx, tend, odt)
    !! ice processes including cloud ice depositional growth, conversion of cloud ice
    !! to snow, snow collecting cloud ice, rain collecting cloud ice, snow depositional growth, 
    !! and graupel sublimation
    use module_mp_tempo_params, only : t0, d0i, bm_i, mu_i, am_i, &
      c_sqrd, c_cube, oig1, cig, d0s, ntb_i, tpi_ide, tps_iaus, tni_iaus, &
      obmi, r_s, ef_si, t1_qs_qi, r_r, org2, cre, t1_qr_qi, t2_qr_qi, &
      fv_r, ef_ri, rho_i, t1_qs_sd, t2_qs_sd, eps, t1_qg_sd, &
      sc3, ogg2, cge, cgg, av_g, rho_w, rho_g

    real(wp), intent(in) :: odt    
    type(ty_tend), intent(inout) :: tend
    logical, dimension(:), intent(in) :: l_qi, l_qs, l_qg
    real(wp), dimension(:), intent(in) :: rhof, rhof2, rho, w1d, ri, ni, rs, rr, nr, &
      temp, qv, qvsi, tcond, diffu, ssati, vsc2, mvd_r, rg, ng
    real(dp), dimension(:), intent(in) :: ilami, smoe, smof, smo1, ilamr, ilamg
    integer, dimension(:), intent(in) :: idx
    real(wp) :: xdi, xmi, oxmi, c_snow, rate_max, otemp, rvs, t2_qg_sd
    real(dp) :: lami, lamr, n0_r, n0_g
    integer :: k, nz, idx_i, idx_i1
    real(wp), dimension(:), allocatable :: t1_subl

    nz = size(l_qi)
    allocate(t1_subl(nz), source=0._wp)
    call get_t1_subl(rho, temp, qvsi, tcond, diffu, ssati, t1_subl)

    do k = 1, nz
      otemp = 1._wp/temp(k)
      rvs = rho(k)*qvsi(k)
      rate_max = (qv(k)-qvsi(k))*rho(k)*odt*0.999_wp

      if (temp(k) < t0) then
        if (l_qi(k)) then
          call get_ice_table_index(ri(k), ni(k), idx_i, idx_i1)
          lami = 1._dp/ilami(k)
          xdi = max(real(d0i, kind=dp), (bm_i + mu_i + 1.) * ilami(k))
          xmi = am_i*xdi**bm_i
          oxmi = 1._wp/xmi
          tend%pri_ide(k) = c_cube*t1_subl(k)*diffu(k)*ssati(k)*rvs &
            *oig1*cig(5)*ni(k)*ilami(k)
          if (tend%pri_ide(k) < 0._dp) then
            tend%pri_ide(k) = max(real(-ri(k)*odt, kind=dp), &
              tend%pri_ide(k), real(rate_max, kind=dp))
            tend%pni_ide(k) = tend%pri_ide(k)*oxmi
            tend%pni_ide(k) = max(real(-ni(k)*odt, kind=dp), tend%pni_ide(k))
          else
            tend%pri_ide(k) = min(tend%pri_ide(k), real(rate_max, kind=dp))
            tend%prs_ide(k) = (1.0_dp-tpi_ide(idx_i,idx_i1))*tend%pri_ide(k)
            tend%pri_ide(k) = tpi_ide(idx_i,idx_i1)*tend%pri_ide(k)
          endif

          ! conversion of cloud ice to snow
          if ((idx_i == ntb_i) .or. (xdi >  5.0_wp*d0s)) then
            tend%prs_iau(k) = ri(k)*.99_wp*odt
            tend%pni_iau(k) = ni(k)*.95_wp*odt
          elseif (xdi < 0.1_wp*d0s) then
            tend%prs_iau(k) = 0._dp
            tend%pni_iau(k) = 0._dp
          else
            tend%prs_iau(k) = tps_iaus(idx_i,idx_i1)*odt
            tend%prs_iau(k) = min(real(ri(k)*.99_wp*odt, kind=dp), tend%prs_iau(k))
            tend%pni_iau(k) = tni_iaus(idx_i,idx_i1)*odt
            tend%pni_iau(k) = min(real(ni(k)*.95_wp*odt, kind=dp), tend%pni_iau(k))
          endif

          ! snow collecting cloud ice assumes di << ds and vti ~ 0
          lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
          xdi = max(real(d0i, kind=dp), (bm_i + mu_i + 1.) * ilami(k))
          xmi = am_i*xdi**bm_i
          oxmi = 1./xmi
          if (rs(k) >= r_s(1)) then
            tend%prs_sci(k) = t1_qs_qi*rhof(k)*ef_si*ri(k)*smoe(k)
            tend%pni_sci(k) = tend%prs_sci(k) * oxmi
          endif

          ! rain collecting cloud ice assumes di << dr and vti= ~ 0
          if (rr(k) >= r_r(1) .and. mvd_r(k) > 4._wp*xdi) then
            lamr = 1._wp/ilamr(k)
            n0_r = nr(k)*org2*lamr**cre(2)
            tend%pri_rci(k) = rhof(k)*t1_qr_qi*ef_ri*ri(k)*n0_r * &
              ((lamr+fv_r)**(-cre(9)))
            tend%pnr_rci(k) = rhof(k)*t1_qr_qi*ef_ri*ni(k)*n0_r * &
              ((lamr+fv_r)**(-cre(9)))
            tend%pnr_rci(k) = min(real(nr(k)*odt, kind=dp), tend%pnr_rci(k))
            tend%png_rci(k) = tend%pnr_rci(k) * &
              max(min((10._wp**(-0.1*w1d(k)) + 0.1_wp), 1._wp), 0.1_wp)
            tend%pni_rci(k) = tend%pri_rci(k) * oxmi
            tend%prr_rci(k) = rhof(k)*t2_qr_qi*ef_ri*ni(k)*n0_r * &
              ((lamr+fv_r)**(-cre(8)))
            tend%prr_rci(k) = min(real(rr(k)*odt, kind=dp), tend%prr_rci(k))
            tend%prg_rci(k) = tend%pri_rci(k) + tend%prr_rci(k)
            tend%pbg_rci(k) = tend%prg_rci(k)/rho_i
          endif
        endif

        if (l_qs(k)) then
          c_snow = c_sqrd + (temp(k)-t0+1.5_wp)*(c_cube-c_sqrd)/(-30._wp+1.5_wp)
          c_snow = max(c_sqrd, min(c_snow, c_cube))
          tend%prs_sde(k) = c_snow*t1_subl(k)*diffu(k)*ssati(k)*rvs * (t1_qs_sd*smo1(k) + &
            t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
          if (tend%prs_sde(k) < 0._dp) then
            tend%prs_sde(k) = max(real(-rs(k)*odt, kind=dp), &
              tend%prs_sde(k), real(rate_max, kind=dp))
          else
            tend%prs_sde(k) = min(tend%prs_sde(k), real(rate_max, kind=dp))
          endif
        endif
        if (l_qg(k)) then
          if (ssati(k) < -eps) then
            n0_g = ng(k)*ogg2*(1._wp/ilamg(k))**cge(2,1)
            t2_qg_sd = 0.28_wp*sc3*sqrt(av_g(idx(k))) * cgg(11,idx(k))
            tend%prg_gde(k) = c_cube*t1_subl(k)*diffu(k)*ssati(k)*rvs &
                * n0_g * (t1_qg_sd*ilamg(k)**cge(10,1) &
                + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11,idx(k)))
            if (tend%prg_gde(k) < 0._wp) then
                tend%prg_gde(k) = max(real(-rg(k)*odt, kind=dp), &
                  tend%prg_gde(k), real(rate_max, kind=dp))
                tend%png_gde(k) = tend%prg_gde(k) * ng(k)/rg(k)
            else
                tend%prg_gde(k) = min(tend%prg_gde(k), real(rate_max, kind=dp))
            endif
          endif 
        endif
      endif 
    enddo
  end subroutine ice_processes


  subroutine melting(rhof2, rho, temp, qvsi, tcond, diffu, vsc2, ssati, &
    delqvs, l_qs, rs, smof, smo0, smo1, l_qg, rg, ng, ilamg, idx, tend, dt, odt)
    !! melting of snow and graupel  
    use module_mp_tempo_params, only : t0, bm_i, mu_i, pi, c_sqrd, c_cube, d0s, &
      ntb_i, r_s, ef_si, r_r, fv_r, ef_ri, rho_i, t1_qs_sd, t2_qs_sd, eps, &
      t1_qg_sd, sc3, ogg2, cge, cgg, av_g, t1_qs_me, t2_qs_me, lvap0, olfus, &
      t1_qg_me, rho_w, rho_g, meters3_to_liters, timestep_conversion_rime_to_rain

    real(wp), intent(in) :: dt, odt
    type(ty_tend), intent(inout) :: tend
    logical, dimension(:), intent(in) :: l_qs, l_qg
    real(wp), dimension(:), intent(in) :: rhof2, rho, rs, &
      temp, qvsi, tcond, diffu, ssati, delqvs, vsc2, rg, ng
    real(dp), dimension(:), intent(in) :: smof, smo0, smo1, ilamg
    integer, dimension(:), intent(in) :: idx
    real(wp) :: tempc, otemp, rvs, melt_f, t2_qg_me, t2_qg_sd
    real(dp) :: n0_g, n0_melt, lamg
    integer :: k, nz
    real(wp), dimension(:), allocatable :: t1_subl

    nz = size(l_qs)
    allocate(t1_subl(nz), source=0._wp)
    call get_t1_subl(rho, temp, qvsi, tcond, diffu, ssati, t1_subl)

    do k = 1, nz
      otemp = 1._wp/temp(k)
      tempc = temp(k) - t0
      rvs = rho(k)*qvsi(k)

      if (temp(k) > t0) then
        if(l_qs(k)) then  
          tend%prr_sml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delqvs(k)) * &
            (t1_qs_me*smo1(k) + t2_qs_me*rhof2(k)*vsc2(k)*smof(k))

          if (tend%prr_sml(k) > 0._dp) then
            tend%prr_sml(k) = tend%prr_sml(k) + 4218._wp*olfus*tempc * &
              (tend%prr_rcs(k)+tend%prs_scw(k))
            tend%prr_sml(k) = min(real(rs(k)*odt, kind=dp), &
              max(0._dp, tend%prr_sml(k)))
            tend%pnr_sml(k) = smo0(k)/rs(k)*tend%prr_sml(k) * 10.0_wp**(-0.25_wp*tempc) 
            tend%pnr_sml(k) = min(real(smo0(k)*odt, kind=dp), tend%pnr_sml(k))
          else
            tend%prr_sml(k) = 0._dp
            tend%pnr_sml(k) = 0._dp
            if (ssati(k) < 0._wp) then
              tend%prs_sde(k) = c_cube*t1_subl(k)*diffu(k)*ssati(k)*rvs * &
                (t1_qs_sd*smo1(k) + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
              tend%prs_sde(k) = max(real(-rs(k)*odt, kind=dp), tend%prs_sde(k))
            endif
          endif
        endif 

        if (l_qg(k)) then
          n0_g = ng(k)*ogg2*(1._wp/ilamg(k))**cge(2,1)
          n0_melt = ng(k)*ogg2*(1._dp/ilamg(k))**cge(2,1)
          if ((rg(k)*ng(k)) < 1.e-4_wp) then
            lamg = 1./ilamg(k)
            n0_melt = (1.e-4_wp/rg(k))*ogg2*lamg**cge(2,1)
          endif
          t2_qg_me = pi*4._wp * c_cube*olfus * &
            0.28_wp*sc3*sqrt(av_g(idx(k))) * cgg(11,idx(k))
          tend%prr_gml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delqvs(k)) * &
            n0_melt*(t1_qg_me*ilamg(k)**cge(10,1) + &
            t2_qg_me*rhof2(k)*vsc2(k)*ilamg(k)**cge(11,idx(k)))
          tend%prr_gml(k) = min(real(rg(k)*odt, kind=dp), max(0._dp, tend%prr_gml(k)))
          if (tend%prr_gml(k) > 0._dp) then
            melt_f = max(0.05_wp, min(tend%prr_gml(k)*dt/rg(k),1._wp))
            ! 1000 is density water, 50 is lower limit (max ice density is 800)
            tend%pbg_gml(k) = meters3_to_liters*tend%prr_gml(k) / &
              max(min(melt_f*rho_g(idx(k)), rho_w), 50._wp)
            tend%pnr_gml(k) = tend%prr_gml(k)*ng(k)/rg(k) * 10.0_wp**(-0.33_wp*(temp(k)-t0))
          else
            tend%prr_gml(k) = 0._dp
            tend%pnr_gml(k) = 0._dp
            tend%pbg_gml(k) = 0._dp
            if (ssati(k) < 0._wp) then
              t2_qg_sd = 0.28_wp*sc3*sqrt(av_g(idx(k))) * cgg(11,idx(k))
              tend%prg_gde(k) = c_cube*t1_subl(k)*diffu(k)*ssati(k)*rvs * n0_g * &
                (t1_qg_sd*ilamg(k)**cge(10,1) + &
                t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11,idx(k)))
              tend%prg_gde(k) = max(real(-rg(k)*odt, kind=dp), tend%prg_gde(k))
              tend%png_gde(k) = tend%prg_gde(k) * ng(k)/rg(k)
            endif 
          endif 
        endif 
        !> @note
        !> if snow and graupel collected cloud water at temperatures above melting
        !> and if the timestep is too long, the melting process should have converted
        !> everything to rain
        !> 
        !> credit to Bjorn-Egil Nygaard for this find
        if (dt > timestep_conversion_rime_to_rain) then
          tend%prr_rcw(k) = tend%prr_rcw(k)+tend%prs_scw(k)+tend%prg_gcw(k)
          tend%prs_scw(k) = 0._dp
          tend%prg_gcw(k) = 0._dp
        endif
      endif 
    enddo
  end subroutine melting


  subroutine aerosol_scavenging(temp, rho, rhof, visco, nwfa, nifa, &
    l_qr, nr, ilamr, mvd_r, l_qs, rs, smob, smoc, smoe, &
    l_qg, rg, ng, ilamg, idx, tend, odt)
    !! scavenging of aerosols by rain, snow, and graupel
    use module_mp_tempo_params, only : d0r, t1_qr_qc, fv_r, cre, &
      org2, r_s, t1_qs_qc, r_g, bm_g, mu_g, av_g, cge, cgg, pi, ogg2

    real(wp), intent(in) :: odt    
    real(wp), dimension(:), intent(in) :: temp, rho,rhof, visco, nr, mvd_r, &
      nwfa, nifa, rs, rg, ng
    real(dp), dimension(:), intent(in) :: smob, smoc, smoe, ilamg, ilamr
    integer, dimension(:), intent(in) :: idx
    logical, dimension(:), intent(in) :: l_qr, l_qs, l_qg
    type(ty_tend), intent(inout) :: tend
    real(wp) :: ef_ra, ef_sa, ef_ga, t1_qg_qc 
    real(dp) :: n0_r, xds, xdg, n0_g, lamr
    real(wp), parameter :: wf_aerosol_size = 0.04e-6_wp
    real(wp), parameter :: if_aerosol_size = 0.8e-6_wp
    integer :: k, nz

    nz = size(l_qr)
    do k = 1, nz
      if (l_qr(k) .and. mvd_r(k).gt. d0r) then
        ef_ra = aerosol_collection_efficiency(real(mvd_r(k), kind=dp), &
          wf_aerosol_size, visco(k), rho(k), temp(k), 'r')
        lamr = 1._dp/ilamr(k)
        n0_r = nr(k)*org2*lamr**cre(2)
        tend%pna_rca(k) = rhof(k)*t1_qr_qc*ef_ra*nwfa(k)*n0_r * &
          ((lamr+fv_r)**(-cre(9)))
        tend%pna_rca(k) = min(real(nwfa(k)*odt, kind=dp), &
          tend%pna_rca(k))
        ef_ra = aerosol_collection_efficiency(real(mvd_r(k), kind=dp), &
          if_aerosol_size, visco(k), rho(k), temp(k), 'r')
        tend%pnd_rcd(k) = rhof(k)*t1_qr_qc*ef_ra*nifa(k)*n0_r * &
          ((lamr+fv_r)**(-cre(9)))
        tend%pnd_rcd(k) = min(real(nifa(k)*odt, kind=dp), &
          tend%pnd_rcd(k))
      endif

      if (l_qs(k) .and. rs(k) > r_s(1)) then
        xds = smoc(k) / smob(k)
        ef_sa = aerosol_collection_efficiency(xds,wf_aerosol_size, &
          visco(k), rho(k), temp(k), 's')
        tend%pna_sca(k) = rhof(k)*t1_qs_qc*ef_sa*nwfa(k)*smoe(k)
        tend%pna_sca(k) = min(real(nwfa(k)*odt, kind=dp), &
          tend%pna_sca(k))
        ef_sa = aerosol_collection_efficiency(xds, if_aerosol_size, &
          visco(k), rho(k), temp(k), 's')
        tend%pnd_scd(k) = rhof(k)*t1_qs_qc*ef_sa*nifa(k)*smoe(k)
        tend%pnd_scd(k) = min(real(nifa(k)*odt, kind=dp), &
          tend%pnd_scd(k))
      endif

      if (l_qg(k) .and. rg(k) > r_g(1)) then
        xdg = (bm_g + mu_g + 1._dp) * ilamg(k)
        ef_ga = aerosol_collection_efficiency(xdg, wf_aerosol_size, &
          visco(k), rho(k), temp(k), 'g')
        t1_qg_qc = pi*.25_wp*av_g(idx(k)) * cgg(9,idx(k))
        n0_g = ng(k)*ogg2*(1._wp/ilamg(k))**cge(2,1)
        tend%pna_gca(k) = rhof(k)*t1_qg_qc*ef_ga*nwfa(k)*n0_g * &
          ilamg(k)**cge(9,idx(k))
        tend%pna_gca(k) = min(real(nwfa(k)*odt, kind=dp), &
          tend%pna_gca(k))
        ef_ga = aerosol_collection_efficiency(xdg, if_aerosol_size, &
          visco(k), rho(k), temp(k), 'g')
        tend%pnd_gcd(k) = rhof(k)*t1_qg_qc*ef_ga*nifa(k)*n0_g * &
          ilamg(k)**cge(9,idx(k))
        tend%pnd_gcd(k) = min(real(nifa(k)*odt, kind=dp), &
          tend%pnd_gcd(k))
      endif
    enddo
  end subroutine aerosol_scavenging

end module module_mp_tempo_main
