module module_mp_tempo_diags
  !! diagnostic output
  use module_mp_tempo_params, only : wp, sp, dp, create_bins, r1, pi
  use module_mp_tempo_utils, only : get_nuc, snow_moments
  
  implicit none
  private
  
  public :: reflectivity_10cm, effective_radius, max_hail_diam, freezing_rain

  contains 

  subroutine effective_radius(temp, l_qc, nc, ilamc, l_qi, ilami, l_qs, rs, &
    re_qc, re_qi, re_qs)
    !! effective radius values for cloud water, cloud ice and snow
    !!
    !! \(r_{e} = 0.5\frac{\int_0^\infty D^{3}n(D)dD}{\int_0^\infty D^2n(D)dD}\)
    use module_mp_tempo_params, only : mu_i, t0

    real(wp), dimension(:), intent(in) :: temp, nc, rs
    real(dp), dimension(:), intent(in) :: ilamc, ilami
    logical, dimension(:), intent(in) :: l_qc, l_qi, l_qs
    real(wp), dimension(:), intent(out) :: re_qc, re_qi, re_qs
    real(wp), dimension(15), parameter :: g_ratio = &
      [24._wp,60._wp,120._wp,210._wp,336._wp,504._wp,720._wp,990._wp, &
      1320._wp,1716._wp,2184._wp,2730._wp,3360._wp,4080._wp,4896._wp]
    real(dp) :: smob, smoc
    real(wp) :: tc0
    integer :: k, nz, nu_c

    nz = size(l_qc)
    do k = 1, nz
       re_qc(k) = 0._wp
       re_qi(k) = 0._wp
       re_qs(k) = 0._wp

      !> @note
      !> limiting values of \(2.51-50 \mu m\) for cloud water, \(2.51-125 \mu m\) for
      !> cloud ice, and \(5.01-999 \mu m\) for snow are consistent with RRTMG radiation
      !> @endnote
      if (l_qc(k)) then
        nu_c = get_nuc(nc(k))
        re_qc(k) = max(2.51e-6_wp, &
          min(real(0.5_dp*(3._dp+real(nu_c, kind=dp))*ilamc(k), kind=wp), 50.e-6_wp))
      endif
      if (l_qi(k)) then
        re_qi(k) = max(2.51e-6_wp, &
          min(real(0.5_dp*(3._dp+real(mu_i, kind=dp))*ilami(k), kind=wp), 125.e-6_wp))
      endif
      if (l_qs(k)) then
        tc0 = min(-0.1, temp(k)-t0)
        call snow_moments(rs=rs(k), tc=tc0, smob=smob, smoc=smoc)
        re_qs(k) = max(5.01e-6_wp, min(0.5_wp*(smoc/smob), 999.e-6_wp))
      endif 
    enddo 
  end subroutine effective_radius


  subroutine reflectivity_10cm(refl10cm_from_melting_flag, &
    temp, l_qr, rr, nr, ilamr, l_qs, rs, smoc, smob, smoz, &
    l_qg, rg, ng, idx, ilamg, dbz)
    !! 10-cm radar reflectivity
    !! 
    !! contributions from melting snow and graupel are optionally included
    !!
    !! \(Z_{e} = \int_0^\infty D^{6}n(D)dD\) and \(dbz = 10*log10(Z_{e}*1\times 10^{18})\)
    use module_mp_tempo_params, only : pi, org2, cre, crg, am_s, am_g, cge, cgg, ogg2

    logical, intent(in) :: refl10cm_from_melting_flag
    logical, dimension(:), intent(in) :: l_qr, l_qs, l_qg
    real(wp), dimension(:), intent(in) :: temp, rg, ng, rr, nr, rs
    real(dp), dimension(:), intent(in) :: ilamr, smoc, smob, smoz, ilamg
    integer, dimension(:), intent(in) :: idx
    real(wp), dimension(:), intent(out) :: dbz
    real(wp) :: ze_rain(size(temp)), ze_snow(size(temp)), ze_graupel(size(temp))
    real(dp) :: n0_r, lamr, n0_g
    integer :: k, nz, k_melt

    nz = size(temp)
    ze_rain = 1.e-22_wp
    ze_snow = 1.e-22_wp
    ze_graupel = 1.e-22_wp

    if (refl10cm_from_melting_flag) then
      k_melt = find_melting_level(temp, l_qr, l_qs, l_qg)
    endif

    do k = nz, 1, -1
      dbz(k) = -35._wp
       
      if (l_qr(k)) then
        lamr = 1._dp/ilamr(k)
        n0_r = nr(k)*org2*lamr**cre(2)
        ze_rain(k) = n0_r*crg(4)*ilamr(k)**cre(4)
      endif
      if (l_qs(k)) then
        ze_snow(k) = (0.176_wp/0.93_wp) * (6._wp/pi)*(6._wp/pi) * &
          (am_s/900._wp)*(am_s/900._wp)*smoz(k)
        ! include melting
        if (refl10cm_from_melting_flag) then
          if (k_melt > 2 .and. k < k_melt-1) then
            ze_snow(k) = reflectivity_from_melting_snow(rs(k), &
              smob(k), smoc(k), rr(k))
          endif
        endif 
      endif
      if (l_qg(k)) then
        n0_g = ng(k)*ogg2*(1._dp/ilamg(k))**cge(2,1)
        ze_graupel(k) = (0.176_wp/0.93_wp) * (6._wp/pi)*(6._wp/pi) * &
          (am_g(idx(k))/900._wp)*(am_g(idx(k))/900._wp) * n0_g*cgg(4,1)*ilamg(k)**cge(4,1)
        ! include melting
        if (refl10cm_from_melting_flag) then
          if (k_melt > 2 .and. k < k_melt-1) then
            ze_graupel(k) = reflectivity_from_melting_graupel(rg(k), ng(k), &
              ilamg(k), idx(k), rr(k))
          endif
        endif 
      endif
      dbz(k) = max(-35._wp, 10._wp*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.e18_dp))
    enddo
  end subroutine reflectivity_10cm


  function find_melting_level(temp, l_qr, l_qs, l_qg) result(k_melt)
    !! finds the melting level
    use module_mp_tempo_params, only : t0

    real(wp), dimension(:), intent(in) :: temp
    logical, dimension(:), intent(in) :: l_qr, l_qs, l_qg
    integer :: k, nz
    integer :: k_melt

    nz = size(l_qr)
    k_melt = 1
    kloop: do k = nz-1, 1, -1
      if ((temp(k) > t0) .and. l_qr(k) .and. (l_qs(k+1) .or. l_qg(k+1))) then
        k_melt = max(k+1, k_melt)
        exit kloop
      endif
    enddo kloop
  end function find_melting_level

  
  function complex_water_ray(lambda, t) result(refractive_index)
    use module_mp_tempo_params, only : pi
    !! complex refractive index of water
    !! from [Ray (1972)](https://doi.org/10.1364/AO.11.001836)
    !! calculated as function of temperature t [Celsius] (valid from -10 to 30)
    !! and radar wavelength lambda [m] (valid from 0.001 to 1)
    !!
    !! original credit: Ulrich Blahak and G. Thompson
    real(dp), intent(in) :: lambda, t
    real(dp) :: epsinf,epss,epsr,epsi,alpha,lambdas,nenner
    complex(dp), parameter :: i = (0._dp, 1._dp)
    complex(dp) :: refractive_index

    epsinf = 5.27137_dp + 0.02164740d0*t - 0.00131198_dp*t*t
    epss = 78.54_dp * (1.0_dp - 4.579e-3_dp * (t - 25.0_dp) + &
      1.190e-5_dp * (t - 25.0_dp)*(t - 25.0_dp) - 2.800e-8_dp * &
      (t - 25.0_dp)*(t - 25.0_dp)*(t - 25.0_dp))
    alpha = -16.8129_dp/(t+273.16_dp) + 0.0609265_dp
    lambdas = 0.00033836_dp * exp(2513.98_dp/(t+273.16_dp)) * 1e-2_dp

    nenner = 1._dp+2._dp*(lambdas/lambda)**(1_dp-alpha)*sin(alpha*pi*0.5_dp) + &
      (lambdas/lambda)**(2._dp-2._dp*alpha)
    epsr = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1_dp-alpha) * &
      sin(alpha*pi*0.5_dp)+1._dp)) / nenner
    epsi = ((epss-epsinf) * ((lambdas/lambda)**(1_dp-alpha) * &
      cos(alpha*pi*0.5_dp)+0._dp)) / nenner + lambda*1.25664_dp/1.88496_dp
    
    refractive_index = sqrt(cmplx(epsr,-epsi, kind=dp))
  end function complex_water_ray


  function complex_ice_maetzler(lambda, t) result(refractive_index)
    !! complex refractive index of ice from
    !! [Maetzler (1998)](https://doi.org/10.1007/978-94-011-5252-5_10)
    !! calculated as function of temperature t [Celsius] (valid from -250 to 0)
    !! and radar wavelength lambda [m] (valid from 0.0001 to 30)
    !!
    !! original credit: Ulrich Blahak and G. Thompson
    real(dp), intent(in) :: lambda, t
    real(dp) :: f,c,tk,b1,b2,b,deltabeta,betam,beta,theta,alfa
    complex(dp) :: refractive_index

    c = 2.99e8_dp
    tk = t + 273.16_dp
    f = c / lambda * 1e-9_dp
    b1 = 0.0207_dp
    b2 = 1.16e-11_dp
    b = 335._dp
    deltabeta = exp(-10.02_dp + 0.0364_dp*(tk-273.16_dp))
    betam = (b1/tk) * (exp(b/tk) / ((exp(b/tk)-1._dp)**2_dp) ) + b2*f*f
    beta = betam + deltabeta
    theta = 300._dp / tk - 1._dp
    alfa = (0.00504_dp + 0.0062_dp*theta) * exp(-22.1_dp*theta)

    refractive_index = 3.1884_dp + 9.1e-4_dp*(tk-273.16_dp)
    refractive_index = refractive_index+ cmplx(0.0_dp, (alfa/f + beta*f), kind=dp) 
    refractive_index = sqrt(conjg(refractive_index))     
  end function complex_ice_maetzler


  function reflectivity_from_melting_graupel(rg, ng, ilamg, idx, rr) result(ze_graupel)
    !! calculates radar reflectivity from melting graupel using binned approach
    !!
    !! original credit: Ulrich Blahak and G. Thompson
    use module_mp_tempo_params, only : am_g, bm_g, obmg, ocmg, mu_g, ogg2, cge, &
      gbins_radar, dgbins_radar, radar_bins

    real(wp), intent(in) :: rg, ng, rr 
    real(dp), intent(in) :: ilamg
    integer, intent(in) :: idx
    real(dp), parameter :: melt_outside = 0.9_dp
    real(dp), parameter :: lambda_radar = 0.10_dp  ! in meters
    complex(dp) :: m_w_0, m_i_0
    real(dp), dimension(radar_bins+1) :: simpson
    real(dp), dimension(3), parameter :: basis = [1._dp/3._dp, 4._dp/3._dp, 1._dp/3._dp]
    real(dp) :: sr, fmelt, eta, lamg, backscatter, f_d, n0_g, mass, k_w
    integer :: n
    real(wp) :: ze_graupel

    do n = 1, radar_bins+1
      simpson(n) = 0._dp
    enddo
    do n = 1, radar_bins-1, 2
      simpson(n) = simpson(n) + basis(1)
      simpson(n+1) = simpson(n+1) + basis(2)
      simpson(n+2) = simpson(n+2) + basis(3)
    enddo

    m_w_0 = complex_water_ray(lambda_radar, 0._dp)
    m_i_0 = complex_ice_maetzler(lambda_radar, 0._dp)
    k_w = (abs((m_w_0*m_w_0 - 1._dp) /(m_w_0*m_w_0 + 2._dp)))**2

    sr = max(0.01_dp, min(1.0_dp - rg/max(rg + rr, r1), 0.99_dp))
    fmelt = real(sr*sr, kind=dp)
    eta = 0._dp
    lamg = 1._dp/ilamg
    n0_g = ng*ogg2*lamg**cge(2,1)

    do n = 1, radar_bins
      mass = am_g(idx) * gbins_radar(n)**bm_g
      call rayleigh_soak_wetgraupel(x_g=mass, a_geo=real(ocmg(idx), kind=dp), b_geo=real(obmg, kind=dp), &
        fmelt=fmelt, lambda_radar=lambda_radar, meltratio_outside=melt_outside, &
        m_w=m_w_0, m_i=m_i_0, backscatter=backscatter)
      f_d = n0_g*gbins_radar(n)**mu_g * exp(-lamg*gbins_radar(n))
      eta = eta + f_d * backscatter * simpson(n) * dgbins_radar(n)
    enddo
    ze_graupel = lambda_radar*lambda_radar*lambda_radar*lambda_radar / &
      (pi*pi*pi*pi*pi * k_w) * eta
  end function reflectivity_from_melting_graupel


  function reflectivity_from_melting_snow(rs, smob, smoc, rr) result(ze_snow)
    !! calculates radar reflectivity from melting snow using binned approach
    !!
    !! original credit: Ulrich Blahak and G. Thompson
    use module_mp_tempo_params, only : am_s, bm_s, lam0, lam1, obms, ocms, kap0, kap1, mu_s, &
      sbins_radar, dsbins_radar, radar_bins

    real(wp), intent(in) :: rs, rr 
    real(dp), intent(in) :: smob, smoc
    real(dp), parameter :: melt_outside = 0.9_dp
    real(dp), parameter :: lambda_radar = 0.10_dp  ! in meters
    complex(dp) :: m_w_0, m_i_0
    real(dp), dimension(radar_bins+1) :: simpson
    real(dp), dimension(3), parameter :: basis = [1._dp/3._dp, 4._dp/3._dp, 1._dp/3._dp]
    real(dp) :: sr, fmelt, eta, backscatter, f_d, mass, k_w, om3, m0, mrat, slam1, slam2
    integer :: n
    real(wp) :: ze_snow

    do n = 1, radar_bins+1
      simpson(n) = 0._dp
    enddo
    do n = 1, radar_bins-1, 2
      simpson(n) = simpson(n) + basis(1)
      simpson(n+1) = simpson(n+1) + basis(2)
      simpson(n+2) = simpson(n+2) + basis(3)
    enddo

    m_w_0 = complex_water_ray(lambda_radar, 0._dp)
    m_i_0 = complex_ice_maetzler(lambda_radar, 0._dp)
    k_w = (abs((m_w_0*m_w_0 - 1._dp) /(m_w_0*m_w_0 + 2._dp)))**2

    sr = max(0.01_dp, min(1.0_dp - rs/(rs + rr), 0.99_dp))
    fmelt = real(sr*sr, kind=dp)
    eta = 0._dp

    om3 = 1._dp/smoc
    m0 = (smob*om3)
    mrat = smob*m0*m0*m0
    slam1 = m0 * lam0
    slam2 = m0 * lam1
    do n = 1, radar_bins
      mass = am_s * sbins_radar(n)**bm_s

      call rayleigh_soak_wetgraupel(x_g=mass, a_geo=real(ocms, kind=dp), b_geo=real(obms, kind=dp), &
        fmelt=fmelt, lambda_radar=lambda_radar, meltratio_outside=melt_outside, &
        m_w=m_w_0, m_i=m_i_0, backscatter=backscatter)
      f_d = mrat*(kap0*exp(-slam1*sbins_radar(n)) + kap1*(m0*sbins_radar(n))**mu_s * exp(-slam2*sbins_radar(n)))
      eta = eta + f_d * backscatter * simpson(n) * dsbins_radar(n)
    enddo
    ze_snow = lambda_radar*lambda_radar*lambda_radar*lambda_radar / &
      (pi*pi*pi*pi*pi * k_w) * eta
  end function reflectivity_from_melting_snow


  subroutine rayleigh_soak_wetgraupel(x_g, a_geo, b_geo, fmelt, lambda_radar, &
    meltratio_outside, m_w, m_i, backscatter)
    !! calculates backscatter cross section of wet snow or graupel
    !! using Maxwell-Garnett mixing formula and Rayleigh approximation
    !!
    !! original credit: Ulrich Blahak and G. Thompson

    real(dp), intent(in) :: x_g, a_geo, b_geo, fmelt, lambda_radar, meltratio_outside
    complex(dp), intent(in) :: m_w, m_i
    real(dp), intent(out) :: backscatter
    real(dp) :: fm, mra, x_w, d_g, vg, rhog, meltratio_outside_grenz, &
      volg, fmgrenz, d_large, volwater, volice, volair, volmix1, volmix2, &
      vol1, vol2, vol3
    complex(dp) :: m_core, m_air, m_tmp, m1t, m2t, m3t, beta2, beta3

    m_air = (1._dp, 0._dp)
    fm = max(min(fmelt, 1._dp), 0._dp) ! limit melt fraction between 0 and 1
    mra = max(min(meltratio_outside, 1._dp), 0._dp) ! ratio of (melt on outside) / (melt on indside)
    mra = mra + (1._dp-mra)*fm ! force mra to 1 when fm = 1
    
    x_w = x_g * fm
    d_g = a_geo * x_g**b_geo
    if (d_g >= r1) then
      vg = pi/6._dp * d_g**3
      rhog = max(min(x_g / vg, 900._dp), 10._dp)
      vg = x_g / rhog
      meltratio_outside_grenz = 1._dp - rhog / 1000._dp

      if (mra <= meltratio_outside_grenz) then
        volg = vg * (1._dp - mra * fm)
      else
        ! at some value of fm all air gets filled with meltwater
        fmgrenz=(900._dp-rhog)/(mra*900._dp-rhog+900._dp*rhog/1000._dp)
        if (fm <= fmgrenz) then
          ! not all air is filled with meltwater
          volg = (1.0_dp - mra * fm) * vg
        else
          ! all air is filled with meltwater
          volg = (x_g - x_w) / 900._dp + x_w / 1000._dp
        endif
      endif

      ! ice-air-water mixture volumes
      d_large = (6._dp / pi * volg) ** (1._dp/3._dp)
      volice = (x_g - x_w) / (volg * 900._dp)
      volwater = x_w / (1000._dp * volg)
      volair = 1.0_dp - volice - volwater    

      volmix1 = volice / max(volice+volwater,1e-10_dp)
      volmix2 = 1._dp - volmix1

      ! Maxwell-Garnett mixing for ice-water mixture
      m1t = m_w**2
      m2t = m_air**2
      m3t = m_i**2
      vol1 = volmix2
      vol2 = 0._dp
      vol3 = volmix1

      beta2 = 2._dp*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*log(m2t/m1t)-1._dp)
      beta3 = 2._dp*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*log(m3t/m1t)-1._dp)

      m_tmp = sqrt(((1._dp-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
       (1._dp-vol2-vol3+vol2*beta2+vol3*beta3))

      ! Maxwell-Garnett mixing (including air)
      m1t = m_tmp**2
      m2t = m_air**2
      m3t = (2._dp*m_air)**2
      vol1 = (1._dp-volair)
      vol2 = volair
      vol3 = 0._dp

      beta2 = 2._dp*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*log(m2t/m1t)-1._dp)
      beta3 = 2._dp*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*log(m3t/m1t)-1._dp)

      ! complex index of refraction for ice-air-water mixture
      m_core = sqrt(((1._dp-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
       (1._dp-vol2-vol3+vol2*beta2+vol3*beta3))

      ! Rayleigh backscattering coefficient of melting particle
      backscatter = (abs((m_core**2-1._dp)/(m_core**2+2._dp)))**2 * pi*pi*pi*pi*pi * d_large**6 / &
        (lambda_radar*lambda_radar*lambda_radar*lambda_radar)
    else 
      backscatter = 0._dp
    endif
  end subroutine rayleigh_soak_wetgraupel


  subroutine max_hail_diam(rho, rg, ng, ilamg, idx, max_hail_diameter)
    !! estimates maximmum hail diameter [mm] using a binned approach
    !! 
    !! see [Jensen et al. (2023)](https://doi.org/10.1175/MWR-D-21-0319.1)
    use module_mp_tempo_params, only : hbins, dhbins, rho_g, ogg2, cge, nhbins, mu_g

    real(wp), dimension(:), intent(in) :: rho, rg, ng
    real(dp), dimension(:), intent(in) :: ilamg
    integer, dimension(:), intent(in) :: idx
    real(wp), dimension(:), intent(out) :: max_hail_diameter
    real(dp) :: lamg, n0_g, sum_nh, sum_t, f_d, hail_max
    integer :: k, nz, n
    real(dp), parameter :: threshold_conc = 0.0005

    nz = size(rho)
    do k = 1, nz
      max_hail_diameter(k) = 0._wp
      if(rg(k)/rho(k) >= 1.e-6_wp) then
        if (rho_g(idx(k)) < 350._wp) cycle ! density too low to be hail/ice pellets
        lamg = 1._dp / ilamg(k)
        n0_g = ng(k)*ogg2*lamg**cge(2,1)

        sum_nh = 0._dp
        sum_t = 0._dp
        do n = nhbins, 1, -1
          f_d = n0_g*hbins(n)**mu_g * exp(-lamg*hbins(n)) * dhbins(n)
          sum_nh = sum_nh + f_d
          if (sum_nh > threshold_conc) exit
          sum_t = sum_nh
        enddo
        if (n >= nhbins) then
          hail_max = hbins(nhbins)
        elseif (hbins(n+1) > 1.e-3_wp) then
          hail_max = hbins(n) - (sum_nh-threshold_conc)/(sum_nh-sum_t) * (hbins(n)-hbins(n+1))
        else
          hail_max = 1.e-4_wp
        endif
        max_hail_diameter(k) = 1000._wp * hail_max ! convert to mm
      endif
    enddo
  end subroutine max_hail_diam


  subroutine freezing_rain(temp, rain_precip, cloud_precip, frz_rain)
    !! estimates freezing rain/drizzle accumulation
    use module_mp_tempo_params, only : t0

    real(wp), intent(in) :: temp, rain_precip
    real(wp), intent(in), optional :: cloud_precip
    real(wp), intent(out) :: frz_rain

    frz_rain = 0._wp
    if ((temp-t0) < -0.5_wp) then
      frz_rain = rain_precip
      if (present(cloud_precip)) then
        frz_rain = frz_rain + cloud_precip
      endif 
    endif   
  end subroutine freezing_rain

end module module_mp_tempo_diags
