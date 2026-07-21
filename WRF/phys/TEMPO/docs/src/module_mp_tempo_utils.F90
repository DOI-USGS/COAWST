module module_mp_tempo_utils
  !! utilities for tempo microphysics
  use module_mp_tempo_params, only : wp, sp, dp
  implicit none
  private
  
  public :: snow_moments, calc_gamma_p, get_nuc, get_constant_cloud_number, &
    calc_rslf, calc_rsif, compute_efrw, compute_efsw, compute_drop_evap, qi_aut_qs

  contains

  subroutine get_constant_cloud_number(land, nc)
    !! returns land-specific value of cloud droplet number concentration
    !! when aerosol-aware = false if land = 1, else returns ocean-specific value
    use module_mp_tempo_params, only : nt_c_l, nt_c_o
    
    integer, intent(in), optional :: land
    real(wp), dimension(:), intent(out) :: nc

    nc = nt_c_l
    if (present(land)) then
      if (land /= 1) nc = nt_c_o
    endif 
  end subroutine get_constant_cloud_number


  function calc_gamma_p(a, x) result(gamma_p)
    !! normalized lower gamma function calculated either with a 
    !! series expansion or continued fraction method
    !!
    !! input: a = gamma function argument, x = upper limit of integration
    !!
    !! output: gamma_p = \(\gamma(a, x) / \Gamma(a)\)
    real(wp), intent(in) :: a, x
    real(wp) :: gamma_p

    if ((x < 0.0_wp) .or. (a <= 0.0_wp)) stop "Invalid arguments for function gamma_p"
    if (x < (a+1.0_wp)) then
      gamma_p = calc_gamma_series(a, x)
    else
      ! gammma_cf computes the upper series
      gamma_p = 1.0_wp - calc_gamma_cf(a, x)
    endif
  end function calc_gamma_p


  function calc_gamma_series(a, x) result(gamma_series)
    !! solves the normalized lower gamma function
    !!
    !! \(\gamma(a,x) / \Gamma(a) = x^{a} * \gamma*(a,x)\)
    !!
    !! see [Equation 8.7.1](https://dlmf.nist.gov/8.7)
    !! \(\gamma(a,x) = exp(-x) * \sum_{k=0}^{\infty} \frac{x^{k}}{\Gamma(a+k+1)}\)
    !!
    !! input: a = gamma function argument, x = upper limit of integration
    !!
    !! output: normalized lower gamma function \(\gamma(a, x) / \Gamma(a)\)
    real(wp), intent(in) :: a, x
    integer :: k
    integer, parameter :: it_max = 100
    real(wp), parameter :: smallvalue = 1.e-7_wp
    real(wp) :: ap1, sum_term, sum
    real(wp) :: gamma_series

    if (x <= 0.0_wp) stop "Invalid arguments for function gamma_series"
    ! k = 0 summation term is 1 / Gamma(a+1)
    ap1 = a
    sum_term = 1.0_wp / gamma(ap1+1.0_wp)
    sum = sum_term
    do k = 1, it_max
      ap1 = ap1 + 1.0_wp
      sum_term = sum_term * x / ap1
      sum = sum + sum_term
      if (abs(sum_term) < (abs(sum) * smallvalue)) exit
    enddo
    if (k == it_max) stop "gamma_series solution did not converge"
    gamma_series = sum * x**a * exp(-x)
  end function calc_gamma_series
    

  function calc_gamma_cf(a, x) result(gamma_cf)
  !! solves the normalized upper gamma function \(\gamma(a,x) / \Gamma(a)\)
  !! using a continued fractions method
  !! [(modified Lentz Algorithm)](http://functions.wolfram.com/06.06.10.0003.01)
  !!
  !!input: a = gamma function argument, x = lower limit of integration
  !!
  !!output: normalized upper gamma function: \(\gamma(a, x) / \Gamma(a)\)
    real(wp), intent(in) :: a, x
    integer :: k
    integer, parameter :: it_max = 100
    real(wp), parameter :: smallvalue = 1.e-7_wp
    real(wp), parameter :: offset = 1.e-30_wp
    real(wp) :: b, d, h0, c, delta, h, aj
    real(wp) :: gamma_cf
  
    b = 1.0_wp - a + x
    d = 1.0_wp / b
    h0 = offset
    c = b + (1.0_wp/offset)
    delta = c * d
    h = h0 * delta

    do k = 1, it_max
      aj = k * (a-k)
      b = b + 2.0_wp
      d = b + aj*d
      if(abs(d) < offset) d = offset
      c = b + aj/c
      if(abs(c) < offset) c = offset
      d = 1.0_wp / d
      delta = c * d
      h = h * delta
      if (abs(delta-1.0_wp) < smallvalue) exit
    enddo
    if (k == it_max) stop "gamma_cf solution did not converge"
    gamma_cf = exp(-x+a*log(x)) * h / gamma(a)
  end function calc_gamma_cf   


  subroutine snow_moments(rs, tc, smob, smoc, ns, smo0, smo1, smo2, smoe, smof, smog, smoz)
    !! computes snow moments from
    !! [Field et al. (2005)](https://doi.org/10.1256/qj.04.134)
    ! smo0 = 0th moment
    ! smo1 = 1st moment
    ! smo2 = 2nd moment
    ! ns   = total number concentration (not smo0)
    ! smob = rs*oams (2nd moment when bm_s=2)
    ! smoc = (bm_s+1)th moment
    ! smoe = (bv_s+2)th moment
    ! smof = (1+(bv_s+1)/2)th moment
    ! smog = (bm_s+bv_s+2)th moment
    ! somz = (bm**2)th moment for reflectivity
    use module_mp_tempo_params, only : bm_s, sa, sb, &
      oams, cse, csg, lam0, lam1, kap0, kap1, mu_s
    
    real(wp), intent(in) :: rs, tc
    real(dp) :: loga_, a_, b_, smo2_, m0, mrat, slam1, slam2
    real(dp), intent(out) :: smob, smoc
    real(dp), intent(out), optional :: ns, smo0, smo1, smo2, smoe, smof, smog, smoz

    ! Second moment and smob
    smob = real(rs*oams, kind=dp)
    if (bm_s > 2.0_wp-1.e-3_wp .and. bm_s < 2.0_wp+1.e-3_wp) then
      loga_ = sa(1) + sa(2)*tc + sa(3)*bm_s &
        + sa(4)*tc*bm_s + sa(5)*tc*tc &
        + sa(6)*bm_s*bm_s + sa(7)*tc*tc*bm_s &
        + sa(8)*tc*bm_s*bm_s + sa(9)*tc*tc*tc &
        + sa(10)*bm_s*bm_s*bm_s
      a_ = 10.0_wp**loga_
      b_ = sb(1) + sb(2)*tc + sb(3)*bm_s &
        + sb(4)*tc*bm_s + sb(5)*tc*tc &
        + sb(6)*bm_s*bm_s + sb(7)*tc*tc*bm_s &
        + sb(8)*tc*bm_s*bm_s + sb(9)*tc*tc*tc &
        + sb(10)*bm_s*bm_s*bm_s
      smo2_ = (smob/a_)**(1._wp/b_)
    else
      smo2_ = smob
    endif
    if (present(smo2)) smo2 = smo2_

    ! bm+1 moment. Useful for diameter calcs.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(1) &
      + sa(4)*tc*cse(1) + sa(5)*tc*tc &
      + sa(6)*cse(1)*cse(1) + sa(7)*tc*tc*cse(1) &
      + sa(8)*tc*cse(1)*cse(1) + sa(9)*tc*tc*tc &
      + sa(10)*cse(1)*cse(1)*cse(1)
    a_ = 10.0_dp**loga_
    b_ = sb(1)+sb(2)*tc+sb(3)*cse(1) + sb(4)*tc*cse(1) &
      + sb(5)*tc*tc + sb(6)*cse(1)*cse(1) &
      + sb(7)*tc*tc*cse(1) + sb(8)*tc*cse(1)*cse(1) &
      + sb(9)*tc*tc*tc+sb(10)*cse(1)*cse(1)*cse(1)
    smoc = a_ * smo2_**b_

    ! 0th moment. Represents snow number concentration.
    loga_ = sa(1) + sa(2)*tc + sa(5)*tc*tc + sa(9)*tc*tc*tc
    a_ = 10.0**loga_
    b_ = sb(1) + sb(2)*tc + sb(5)*tc*tc + sb(9)*tc*tc*tc
    if (present(smo0)) smo0 = a_ * smo2_**b_

    ! 1st moment. Useful for depositional growth and melting.
    loga_ = sa(1) + sa(2)*tc + sa(3) &
      + sa(4)*tc + sa(5)*tc*tc &
      + sa(6) + sa(7)*tc*tc &
      + sa(8)*tc + sa(9)*tc*tc*tc &
      + sa(10)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3) + sb(4)*tc &
      + sb(5)*tc*tc + sb(6) &
      + sb(7)*tc*tc + sb(8)*tc &
      + sb(9)*tc*tc*tc + sb(10)
    if (present(smo1)) smo1 = a_ * smo2_**b_

    ! snow number concentration (explicit integral, not smo0)
    m0 = smob/smoc
    mrat = smob*m0*m0*m0
    slam1 = m0 * lam0
    slam2 = m0 * lam1
    if (present(ns)) then
      ns = mrat*kap0/slam1 + mrat*kap1*m0**mu_s*csg(15)/slam2**cse(15)
    endif 

    ! bv_s+2 (th) moment. Useful for riming.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(13) &
      + sa(4)*tc*cse(13) + sa(5)*tc*tc &
      + sa(6)*cse(13)*cse(13) + sa(7)*tc*tc*cse(13) &
      + sa(8)*tc*cse(13)*cse(13) + sa(9)*tc*tc*tc &
      + sa(10)*cse(13)*cse(13)*cse(13)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3)*cse(13) + sb(4)*tc*cse(13) &
      + sb(5)*tc*tc + sb(6)*cse(13)*cse(13) &
      + sb(7)*tc*tc*cse(13) + sb(8)*tc*cse(13)*cse(13) &
      + sb(9)*tc*tc*tc + sb(10)*cse(13)*cse(13)*cse(13)
    if (present(smoe)) smoe = a_ * smo2_**b_

    ! 1+(bv_s+1)/2 (th) moment.  Useful for depositional growth.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(16) &
      + sa(4)*tc*cse(16) + sa(5)*tc*tc &
      + sa(6)*cse(16)*cse(16) + sa(7)*tc*tc*cse(16) &
      + sa(8)*tc*cse(16)*cse(16) + sa(9)*tc*tc*tc &
      + sa(10)*cse(16)*cse(16)*cse(16)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3)*cse(16) + sb(4)*tc*cse(16) &
      + sb(5)*tc*tc + sb(6)*cse(16)*cse(16) &
      + sb(7)*tc*tc*cse(16) + sb(8)*tc*cse(16)*cse(16) &
      + sb(9)*tc*tc*tc + sb(10)*cse(16)*cse(16)*cse(16)
    if (present(smof)) smof = a_ * smo2_**b_

    ! bm_s + bv_s+2 (th) moment.  Useful for riming into graupel.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(17) &
        + sa(4)*tc*cse(17) + sa(5)*tc*tc &
        + sa(6)*cse(17)*cse(17) + sa(7)*tc*tc*cse(17) &
        + sa(8)*tc*cse(17)*cse(17) + sa(9)*tc*tc*tc &
        + sa(10)*cse(17)*cse(17)*cse(17)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3)*cse(17) + sb(4)*tc*cse(17) &
        + sb(5)*tc*tc + sb(6)*cse(17)*cse(17) &
        + sb(7)*tc*tc*cse(17) + sb(8)*tc*cse(17)*cse(17) &
        + sb(9)*tc*tc*tc + sb(10)*cse(17)*cse(17)*cse(17)
    if (present(smog)) smog = a_ * smo2_**b_

    !..Calculate bm_s*2 (th) moment.  Useful for reflectivity.
    loga_ = sa(1) + sa(2)*tc + sa(3)*cse(3) &
      + sa(4)*tc*cse(3) + sa(5)*tc*tc &
      + sa(6)*cse(3)*cse(3) + sa(7)*tc*tc*cse(3) &
      + sa(8)*tc*cse(3)*cse(3) + sa(9)*tc*tc*tc &
      + sa(10)*cse(3)*cse(3)*cse(3)
    a_ = 10.0**loga_
    b_ = sb(1)+ sb(2)*tc + sb(3)*cse(3) + sb(4)*tc*cse(3) &
      + sb(5)*tc*tc + sb(6)*cse(3)*cse(3) &
      + sb(7)*tc*tc*cse(3) + sb(8)*tc*cse(3)*cse(3) &
      + sb(9)*tc*tc*tc + sb(10)*cse(3)*cse(3)*cse(3)
    if (present(smoz)) smoz = a_ * smo2**b_

  end subroutine snow_moments


  function calc_rslf(p, t) result(rslf)
    !! calculates liquid saturation vapor mixing ratio
    real(wp), intent(in) :: p, t
    real(wp) :: esl, x
    real(wp), parameter :: c0 = .611583699E03_wp
    real(wp), parameter :: c1 = .444606896E02_wp
    real(wp), parameter :: c2 = .143177157e01_wp
    real(wp), parameter :: c3 = .264224321e-1_wp
    real(wp), parameter :: c4 = .299291081e-3_wp
    real(wp), parameter :: c5 = .203154182e-5_wp
    real(wp), parameter :: c6 = .702620698e-8_wp
    real(wp), parameter :: c7 = .379534310e-11_wp
    real(wp), parameter :: c8 = -.321582393e-13_wp
    real(wp) :: rslf

    x = max(-80._wp, t-273.16_wp)
    esl = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    esl = min(esl, p*0.15) 
    !! @note
    !! even with p = 1050 mb and t = 55 C, the saturation vapor 
    !! pressure only contributes to 15% of the total pressure,
    !! thus the limit on the saturation vapor pressure is set to 15%
    !! @endnote
    rslf = .622*esl/(p-esl)
  end function calc_rslf


  function calc_rsif(p, t) result(rsif)
    !! calculates liquid saturation vapor mixing ratio
    real(wp), intent(in) :: p, t
    real(wp) :: esi, x
    real(wp), parameter :: c0 = .609868993e03_wp
    real(wp), parameter :: c1 = .499320233e02_wp
    real(wp), parameter :: c2 = .184672631e01_wp
    real(wp), parameter :: c3 = .402737184e-1_wp
    real(wp), parameter :: c4 = .565392987e-3_wp
    real(wp), parameter :: c5 = .521693933e-5_wp
    real(wp), parameter :: c6 = .307839583e-7_wp
    real(wp), parameter :: c7 = .105785160e-9_wp
    real(wp), parameter :: c8 = .161444444e-12_wp
    real(wp) :: rsif

    x = max(-80._wp, t-273.16_wp)
    esi = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    esi = min(esi, p*0.15)
    rsif = .622*esi/max(1.e-4_wp,(p-esi))
  end function calc_rsif


  function get_nuc(nc) result(nu_c)
    !! returns nu_c for cloud water (values from 2-15)
    use module_mp_tempo_params, only : nu_c_scale

    real(wp), intent(in) :: nc
    integer :: nu_c

    if ((nu_c_scale/nc) >= 12.5_wp) then
      nu_c = 15
    else
      nu_c = min(15, nint(nu_c_scale/nc) + 2)
    endif 
  end function get_nuc


  subroutine compute_efrw()
    !! collision efficiency for rain collecting cloud water from Beard and Grover (1974)
    !! if a/A < 0.25
    !! https://doi.org/10.1175/1520-0469(1974)031<0543:NCEFSR>2.0.CO;2
    !! otherwise uses polynomials to get close match of Pruppacher and Klett Fig. 14-9
    use module_mp_tempo_params, only : nbc, nbr, dc, dr, t_efrw, rho_w, pi

    real(dp) :: vtr, stokes, reynolds, ef_rw
    real(dp) :: p, yc0, f, g, h, z, k0, x
    integer :: i, j

    do j = 1, nbc
      do i = 1, nbr
        ef_rw = 0.0_dp
        p = dc(j) / dr(i)
        if (dr(i) < 50.e-6_dp .or. dc(j) < 3.e-6_dp) then
          t_efrw(i,j) = 0.0_dp
        elseif (p > 0.25_dp) then
          x = dc(j) * 1.e6_dp
          if (dr(i) < 75.e-6_dp) then
            ef_rw = 0.026794_dp*x - 0.20604_dp
          elseif (dr(i) < 125.e-6_dp) then
            ef_rw = -0.00066842_dp*x*x + 0.061542_dp*x - 0.37089_dp
          elseif (dr(i) < 175.e-6_dp) then
            ef_rw = 4.091e-06_dp*x*x*x*x - 0.00030908_dp*x*x*x + &
              0.0066237_dp*x*x - 0.0013687_dp*x - 0.073022_dp
          elseif (dr(i) < 250.e-6_dp) then
            ef_rw = 9.6719e-5_dp*x*x*x - 0.0068901_dp*x*x + 0.17305_dp*x - 0.65988_dp
          elseif (dr(i) < 350.e-6_dp) then
            ef_rw = 9.0488e-5_dp*x*x*x - 0.006585_dp*x*x + 0.16606_dp*x - 0.56125_dp
          else
            ef_rw = 0.00010721_dp*x*x*x - 0.0072962_dp*x*x + 0.1704_dp*x - 0.46929_dp
          endif
        else
          vtr = -0.1021_dp + 4.932e3_dp*dr(i) - 0.9551e6_dp*dr(i)*dr(i) + 0.07934e9_dp*dr(i)*dr(i)*dr(i) &
            - 0.002362e12_dp*dr(i)*dr(i)*dr(i)*dr(i)
          stokes = dc(j) * dc(j) * vtr * rho_w / (9._dp*1.718e-5_dp*dr(i))
          reynolds = 9._dp * stokes / (p*p*rho_w)

          f = log(reynolds)
          g = -0.1007_dp - 0.358_dp*f + 0.0261_dp*f*f
          k0 = exp(g)
          z = log(stokes / (k0+1.e-15_dp))
          h = 0.1465_dp + 1.302_dp*z - 0.607_dp*z*z + 0.293_dp*z*z*z
          yc0 = 2.0_dp / pi * atan(h)
          ef_rw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))
        endif

        t_efrw(i,j) = max(0.0_dp, min(ef_rw, 0.95_dp))
      enddo
    enddo
  end subroutine compute_efrw


  subroutine compute_efsw()
    !! collision efficiency for snow collecting cloud water from Wang and Ji (2000)
    !! https://doi.org/10.1175/1520-0469(2000)057<1001:CEOICA>2.0.CO;2
    !! equating melted snow diameter to effective collision cross-section
    use module_mp_tempo_params, only : wp, sp, dp, &
      nbc, dc, av_s, bv_s, ds, nbs, am_s, bm_s, am_r, obmr, fv_s, &
      t_efsw, d0s, rho_w, pi

    real(dp) :: ds_m, vts, vtc, stokes, reynolds, ef_sw
    real(dp) :: p, yc0, f, g, h, z, k0
    integer :: i, j

    do j = 1, nbc
      vtc = 1.19e4_dp * (1.0e4_dp*dc(j)*dc(j)*0.25_dp)
      do i = 1, nbs
        vts = av_s*ds(i)**bv_s * exp(real(-fv_s*ds(i), kind=dp)) - vtc
        ds_m = (am_s*ds(i)**bm_s / am_r)**obmr
        p = dc(j) / ds_m

        if (p > 0.25_dp .or. ds(i) < d0s .or. dc(j) < 6.e-6_dp .or. vts < 1.e-3_dp) then
          t_efsw(i,j) = 0.0_dp
        else
          stokes = dc(j) * dc(j) * vts * rho_w / (9.*1.718e-5_dp*ds_m)
          reynolds = 9._dp * stokes / (p*p*rho_w)

          f = log(reynolds)
          g = -0.1007_dp - 0.358_dp*f + 0.0261_dp*f*f
          k0 = exp(g)
          z = log(stokes / (k0+1.e-15_dp))
          h = 0.1465_dp + 1.302_dp*z - 0.607_dp*z*z + 0.293_dp*z*z*z
          yc0 = 2.0_dp / pi * atan(h)
          ef_sw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))

          t_efsw(i,j) = max(0.0_dp, min(ef_sw, 0.95_dp))
        endif
      enddo
    enddo
  end subroutine compute_efsw


  subroutine qi_aut_qs()
    !! calculates cloud ice conversion to snow
    !! and depositional growth by binning cloud ice distributions and 
    !! determining both the size bins > d0s (that are converted to snow) and the 
    !! depositional growth up to d0s (for cloud ice) and > d0s for snow 
    !! following Harrington et al. (1995)
    !! https://doi.org/10.1175/1520-0469(1995)052<4344:POICCP>2.0.CO;2
    use module_mp_tempo_params, only : wp, sp, dp, &
      nbi, ntb_i, ntb_i1, &
      am_i, cie, cig, oig1, nt_i, r_i, obmi, bm_i, mu_i, d0s, d0i, &
      tpi_ide, tps_iaus, tni_iaus, di, dti

    integer :: i, j, n2
    real(dp), dimension(nbi) :: n_i
    real(dp) :: n0_i, lami, di_mean, t1, t2
    real(wp) :: xlimit_intg

    do j = 1, ntb_i1
      do i = 1, ntb_i
        lami = (am_i*cig(2)*oig1*nt_i(j)/r_i(i))**obmi
        di_mean = (bm_i + mu_i + 1.) / lami
        n0_i = nt_i(j)*oig1 * lami**cie(1)
        t1 = 0.
        t2 = 0.
        if (real(di_mean, kind=wp) > 5.*d0s) then
          t1 = r_i(i)
          t2 = nt_i(j)
          tpi_ide(i,j) = 0.
        elseif (real(di_mean, kind=wp) < d0i) then
          t1 = 0.
          t2 = 0.
          tpi_ide(i,j) = 1.
        else
          xlimit_intg = lami*d0s
          tpi_ide(i,j) = real(calc_gamma_p(mu_i+2.0, xlimit_intg), kind=dp)
          do n2 = 1, nbi
            n_i(n2) = n0_i*di(n2)**mu_i * exp(-lami*di(n2))*dti(n2)
            if (di(n2) >= d0s) then
              t1 = t1 + n_i(n2) * am_i*di(n2)**bm_i
              t2 = t2 + n_i(n2)
            endif
          enddo
        endif
        tps_iaus(i,j) = t1
        tni_iaus(i,j) = t2
      enddo
    enddo
  end subroutine qi_aut_qs


  subroutine compute_drop_evap()
    !! calculates droplet evaporation data
    use module_mp_tempo_params, only: wp, sp, dp, &
      nbc, am_r, dc, bm_r, nu_c_scale, nu_c_max, nu_c_min, &
      t_nc, ntb_c, cce, ccg, ocg1, r_c, obmr, dtc, &
      tpc_wev, tnc_wev

    integer :: i, j, k, n
    real(dp), dimension(nbc) :: n_c, massc
    real(dp) :: summ, summ2, lamc, n0_c
    integer :: nu_c

    do n = 1, nbc
      massc(n) = am_r*dc(n)**bm_r
    enddo

    do k = 1, nbc
      nu_c = get_nuc(real(t_nc(k), kind=wp))
      do j = 1, ntb_c
        lamc = (t_nc(k)*am_r* ccg(2,nu_c)*ocg1(nu_c) / r_c(j))**obmr
        n0_c = t_nc(k)*ocg1(nu_c) * lamc**cce(1,nu_c)
        do i = 1, nbc
          n_c(i) = n0_c* dc(i)**nu_c*exp(-lamc*dc(i))*dtc(i)
          summ = 0._dp
          summ2 = 0._dp
          do n = 1, i
            summ = summ + massc(n)*n_c(n)
            summ2 = summ2 + n_c(n)
          enddo
          tpc_wev(i,j,k) = summ
          tnc_wev(i,j,k) = summ2
        enddo
      enddo
    enddo
  end subroutine compute_drop_evap

end module module_mp_tempo_utils
