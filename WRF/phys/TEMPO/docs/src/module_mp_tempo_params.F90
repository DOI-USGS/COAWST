module module_mp_tempo_params
  !! parameters and variables used in tempo microphysics

! define machine precision
#if defined(tempo_mpas)
  use mpas_kind_types, only : wp => RKIND, sp => R4KIND, dp => R8KIND
#elif defined(tempo_ccpp)
  use machine, only : wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#elif defined(wrfmodel)
  use ccpp_kind_types, only : wp => kind_phys, sp => kind_phys
#else
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
#endif 
  use iso_fortran_env, only : real32, real64 ! for machine-independent lookup table precisions

  implicit none

  public

#if defined(wrfmodel)
  integer, parameter :: dp = selected_real_kind(15,307)
#endif

  character(len=11) :: tempo_version !! tempo version string (max is xxx.xxx.xxx)

  ! parameters that can be changed ------------------------------------------------------------------------
  integer, parameter :: idx_bg1 = 6 !! index from rho_g when hail_aware = false: density = 500 \(kg\, m^{-3}\)
 
    !> @note
    !> av_r is the rain fallspeed power-law coefficient
    !>
    !> fallspeed power law relations are
    !>
    !> \[ v =(a_{v}D^{b_{v}})\exp\left(-f_{v}D\right), \textrm{where}\,fv = 0\, \textrm{for graupel/ice} \]
    !>
    !> and coefficients are from from [Ferrier (1994)](https://doi.org/10.1175/1520-0469(1994)051<0249:ADMMPF>2.0.CO;2) for rain and 
    !> [Thompson et al. (2008)](https://doi.org/10.1175/2008MWR2387.1) for ice, snow, and graupel
    !> @endnote
  real(wp), parameter :: av_r = 4854._wp
  real(wp), parameter :: bv_r = 1.0_wp !! rain fallspeed power-law coefficient
  real(wp), parameter :: av_s = 40._wp !! snow fallspeed power-law coefficient
  real(wp), parameter :: bv_s = 0.55_wp !! snow fallspeed power-law coefficient
  real(wp), parameter :: fv_s = 100._wp !! snow fallspeed power-law coefficient
  real(wp), parameter :: bv_c = 2.0_wp !! cloud fallspeed power-law coefficient
  real(wp), parameter :: bv_i = 1.0_wp !! ice fallspeed power-law coefficient
  real(wp), parameter :: av_g_old = 442._wp !! graupel fallspeed power-law coefficient (hail_aware = false)
  real(wp), parameter :: bv_g_old = 0.89_wp !! graupel fallspeed power-law coefficient (hail_aware = false)
  real(wp), parameter :: av_g_new = 161.794724_wp !! graupel fallspeed power-law coefficient HARDCODED for idx_bg1 = 6 (hail_aware = true)
  real(wp), parameter :: bv_g_new = 0.640961647_wp !! graupel fallspeed power-law coefficient HARDCODED for idx_bg1 = 6 (hail_aware = true)
  real(wp), parameter :: fv_r = 195.0_wp !! rain fallspeed power-law coefficient
  real(wp), parameter :: av_c = 0.316946e8_wp !! cloud fallspeed power-law coefficient
  real(wp), parameter :: a_coeff = 0.47244157_wp !! graupel fallspeed power-law coefficient
  real(wp), parameter :: b_coeff = 0.54698726_wp !! grapuel fallspeed power-law coefficient
  real(wp), parameter :: av_i = 1493.9 !! ice fallspeed power-law coefficient
    !> @note
    !> am_s is the snow mass power-law coefficient
    !>
    !> mass power law relations are
    !> \[ m = a_{m}D^{b_{m}} \]
    !> and coefficients for snow are from
    !> [Field et al. (2005)](https://doi.org/10.1256/qj.04.134)
    !> and others assume a spherical form
    !> @endnote
  real(wp), parameter :: am_s = 0.069_wp 
  real(wp), parameter :: bm_s = 2.0_wp !! snow mass power-law coefficient
  real(wp), parameter :: bm_g = 3.0_wp !! graupel mass power-law coefficient
  real(wp), parameter :: bm_i = 3.0_wp !! ice mass power-law coefficient
  real(wp), parameter :: bm_r = 3.0_wp !! rain mass power-law coefficient

  real(wp), parameter :: rho_i = 890._wp !! density of cloud ice \([kg\, m^{-3}]\)
  real(wp), parameter :: xm0i = 1.e-12_wp !! ice initiates with this mass \([kg]\)
  real(wp), parameter :: d0c = 1.e-6_wp !! minimum diameter of cloud droplets \([m]\)
  real(wp), parameter :: d0r = 50.e-6_wp !! minimum diameter of raindrops \([m]\)
  real(wp), parameter :: d0s = 300.e-6_wp !! minimum diameter of snow \([m]\)
  real(wp), parameter :: d0g = 350.e-6_wp !! minimum diameter of graupel \([m]\)
  real(wp), parameter :: d0r_max = 2.5e-3_wp !! maximum diameter of raindrops \([m]\)

  real(wp), parameter :: c_cube = 0.5_wp !! capacitance of a sphere \(\left(D^{3}\right)\)
  real(wp), parameter :: c_sqrd = 0.15_wp !! capacitance of plates/aggregates \(\left(D^{2}\right)\)

    !> @note
    !> mu_r is the shape parameter for rain
    !> generalized gamma distributions for rain, graupel and cloud ice have the form
    !>
    !> \[ n\left(D\right) = n_{0} D^{\mu}\exp(-\lambda D) \]
    !> \[ \textrm{where}\\ \mu = 0 \\ \textrm{is exponential} \]
    !> @endnote
  real(wp), parameter :: mu_r = 0.0_wp 
  real(wp), parameter :: mu_s = 0.6357_wp !! shape parameter for snow
  real(wp), parameter :: mu_g = 0.0_wp !! shape parameter for graupel
  real(wp), parameter :: mu_i = 0.0_wp !! shape parameter for cloud ice

  real(wp), parameter :: nu_c_scale = 1000.e6_wp !! scaling parameter for nu_c
  integer, parameter :: nu_c_max = 15 !! maximum value for nu_c
  integer, parameter :: nu_c_min = 2 !! minimum value for nu_c

  real(wp), parameter :: naccn0 = 300.0e6_wp !! used for water-friendly aerosol initialization
  real(wp), parameter :: naccn1 = 50.0e6_wp !! used for water-friendly aerosol initialization
  real(wp), parameter :: nain0 = 1.5e6_wp !! used for ice-friendly aerosol initialization
  real(wp), parameter :: nain1 = 0.5e6_wp !! used for ice-friendly aerosol initialization
  real(wp), parameter :: nwfa_default = 11.1e6_wp !! default value for water-friendly aerosols
  real(wp), parameter :: nifa_default = nain1*0.01_wp !! default value for ice-friendly aerosols
  real(wp), parameter :: aero_max = 9999.e6_wp !! maximum aerosol value
  real(wp), parameter :: hgfrz = 235.16_wp !! temperature to freeze all liquid \([K]\)
  real(wp), parameter :: nt_c_o = 50.e6_wp !! cloud number concentration over ocean (non-aerosol aware) \([m^{-3}]\)
  real(wp), parameter :: nt_c_l = 100.e6_wp !! cloud number concentration over land (non-aerosol aware) \([m^{-3}]\)
  real(wp), parameter :: nt_c_max = 1999.e6_wp !! maximum cloud number concentration \([m^{-3}]\)
  real(wp), parameter :: nt_c_min = 2._wp !! minimum cloud number concentration \([m^{-3}]\)

  real(wp), parameter :: tno = 5.0_wp !! constant in the [Cooper](https://doi.org/10.1007/978-1-935704-17-1_4) curve for ice nucleation
  real(wp), parameter :: ato = 0.304_wp !! constant in the [Cooper](https://doi.org/10.1007/978-1-935704-17-1_4) curve for ice nucleation
  real(wp) :: rho_s = 100.0_wp !! density of snow \([kg\, m^{-3}]\)

  real(wp), parameter :: demott_nuc_ssati = 0.25_wp !! ice supersaturation threshold for [DeMott](https://doi.org/10.1073/pnas.0910818107) nucleation
  real(dp), parameter :: max_ni = 4999.e3_wp !! maximum ice number concentration \([m^{-3}]\)
  real(wp), parameter :: icenuc_max = 1000.e3_wp !! maximum ice nucleation number \([m^{-3}]\)
  real(wp), parameter :: rime_threshold = 2.0_wp !! snow to graupel rime threshold parameter
  real(wp), parameter :: rime_conversion = 0.95_wp !! snow to graupel rime conversion parameter
  real(wp), parameter :: ef_si = 0.05_wp !! snow-ice collection efficiency
  real(wp), parameter :: ef_rs = 0.95_wp !! rain-snow collection efficiency
  real(wp), parameter :: ef_rg = 0.75_wp !! rain-graupel collection efficiency
  real(wp), parameter :: ef_ri = 0.95_wp !! rain-ice collection efficiency
  real(wp), parameter :: autocon_nr_factor = 10._wp !! factor controlling rain number tendency from autconversion (larger produces few drops)
  real(wp), parameter :: timestep_conversion_rime_to_rain = 120._wp !! timestep above which rime above freezing becomes rain (this timestep should be at least 120s)

  ! parameters that should NOT be changed -----------------------------------------------------------------
  integer, parameter :: table_sp = real32 !! precision for lookup tables (machine independent)
  integer, parameter :: table_dp = real64 !! precision for lookup tables (machine independent)

  integer, parameter :: nrhg = 9 !! graupel density array size when hail_aware = true
  integer, parameter :: nrhg1 = 1 !! graupel density array size when hail_aware = false

  real(wp), parameter :: min_qv = 1.e-10_wp !! minimum value of water vapor mixing ratio \([kg\, kg^{-1}]\)
  real(wp), parameter :: r1 = 1.e-12_wp !! minimum hydrometeor mass \([kg\, m^{-3}]\) 
  real(wp), parameter :: r2 = 1.e-6_wp !! minimum hydrometeor number \([kg\, m^{-3}]\)
  real(wp), parameter :: low_limit_mass_for_precip = 1.e-9_wp !! minimum hydrometor mass needed in the lowest-model level for precipitation
  real(wp), parameter :: eps = 1.e-15_wp !! small non-zero number
  real(wp), parameter :: meters3_to_liters = 1000._wp !! number of liters in 1 \(m^{3}\)
  real(dp), parameter :: gonv_min = 1.e2_dp !! minimum graupel y-intercept \([m^{-4}]\)
  real(dp), parameter :: gonv_max = 1.e6_dp !! maximum graupel y-intercept \([m^{-4}]\)

  real(wp), parameter :: t0 = 273.15_wp !! melting point of ice \([K]\)
  real(wp), parameter :: rho_w = 1000._wp !! density of liquid water \([kg\, m^{-3}]\)

  real(wp), dimension(nrhg), parameter :: rho_g = [50._wp, 100._wp, 200._wp, 300._wp, 400._wp, &
    500._wp, 600._wp, 700._wp, 800._wp] !! !! densities of graupel when hail_aware = true \([kg\, m^{-3}]\)

  real(wp), parameter :: sc = 0.632_wp !! [schmidt number](https://glossary.ametsoc.org/wiki/Schmidt_number)
  real(wp), parameter :: earth_gravity = 9.8_wp !! gravity of Earth \([m\, s^{-2}]\)

  ! these can be overwritten by a host model and don't have a parameter attribute
  real(wp) :: pi = 3.1415926536_wp !! pi is approximately 355/113
  real(wp) :: lsub = 2.834e6_wp !! enthalpy of sublimation \([J\, kg^{-1}]\)
  real(wp) :: lvap0 = 2.5e6_wp !! enthalpy of vaporization \([J\, kg^{-1}]\)
  real(wp) :: rv = 461.5_wp !! gas constant for water vapor \([J\, K^{-1}\, kg^{-1}]\)
  real(wp) :: rdry = 287.04_wp !! gas constant for dry air \([J\, K^{-1}\, kg^{-1}]\)
  real(wp) :: roverrv = 0.622_wp !! dry gas constant divided by water vapor gas constant
  real(wp) :: r = 287.04_wp !! gas constant for dry air \([J\, K^{-1}\, kg^{-1}]\)
  real(wp) :: rho_not !! density constant \([kg\, m^{-3}]\)
  real(wp) :: rho_not0 !! density constant \([kg\, m^{-3}]\)
  real(wp) :: cp = 1004.0_wp !! heat capacity of air at constant pressure \([J\, K^{-1}\, kg^{-1}]\)
  real(wp) :: r_uni = 8.314  !! gas constant \([J\, K^{-1}\, mol^{-1}]\)

  real(wp), parameter :: kap0 = 490.6_wp !! snow parameter from [Field et al. (2005)](https://doi.org/10.1256/qj.04.134)
  real(wp), parameter :: kap1 = 17.46_wp !! snow parameter from [Field et al. (2005)](https://doi.org/10.1256/qj.04.134)
  real(wp), parameter :: lam0 = 20.78_wp !! snow parameter from [Field et al. (2005)](https://doi.org/10.1256/qj.04.134)
  real(wp), parameter :: lam1 = 3.29_wp !! snow parameter from [Field et al. (2005)](https://doi.org/10.1256/qj.04.134)
  
  ! lookup table dimensions
  integer, parameter :: nbins = 100 !! lookup table dimension (number of bins)
  integer, parameter :: nbc = nbins !! lookup table dimension for cloud water
  integer, parameter :: nbr = nbins !! lookup table dimension for rain
  integer, parameter :: nbs = nbins !! lookup table dimension for snow
  integer, parameter :: nbi = nbins !! lookup table dimension
  integer, parameter :: nbg = nbins !! lookup table dimension
  integer, parameter :: ntb_i = 64 !! lookup table dimension for cloud ice
  integer, parameter :: ntb_i1 = 55 !! lookup table dimension for cloud ice
  integer, parameter :: ntb_c = 37 !! lookup table dimension for cloud water
  integer, parameter :: ntb_t = 9 !! lookup table dimension for temperature
  integer, parameter :: ntb_g1 = 37 !! lookup table dimension for graupel
  integer, parameter :: ntb_s = 37 !! lookup table dimension for snow
  integer, parameter :: ntb_g = 37 !! lookup table dimension for graupel
  integer, parameter :: ntb_r = 37 !! lookup table dimension for rain
  integer, parameter :: ntb_r1 = 37 !! lookup table dimension for rain
  integer, parameter :: ntb_t1 = 45 !! lookup table dimension for temperature
  integer, parameter :: ntb_in = 55 !! lookup table dimension for IN
  integer, parameter :: ntb_arc = 7 !! lookup table dimension for CCN activation
  integer, parameter :: ntb_arw = 9 !! lookup table dimension for CCN activation
  integer, parameter :: ntb_art = 7 !! lookup table dimension for CCN activation
  integer, parameter :: ntb_arr = 5 !! lookup table dimension for CCN activation
  integer, parameter :: ntb_ark = 4 !! lookup table dimension for CCN activation

  ! lookup table data
  real(wp), dimension(ntb_c), parameter :: &
    r_c = [1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp,2.e-3_wp,3.e-3_wp,4.e-3_wp,5.e-3_wp,6.e-3_wp,7.e-3_wp,8.e-3_wp,9.e-3_wp, &
    1.e-2_wp] !! mass bins for cloud water \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_i), parameter :: &
    r_i = [1.e-10_wp,2.e-10_wp,3.e-10_wp,4.e-10_wp, &
    5.e-10_wp,6.e-10_wp,7.e-10_wp,8.e-10_wp,9.e-10_wp, &
    1.e-9_wp,2.e-9_wp,3.e-9_wp,4.e-9_wp,5.e-9_wp,6.e-9_wp,7.e-9_wp,8.e-9_wp,9.e-9_wp, &
    1.e-8_wp,2.e-8_wp,3.e-8_wp,4.e-8_wp,5.e-8_wp,6.e-8_wp,7.e-8_wp,8.e-8_wp,9.e-8_wp, &
    1.e-7_wp,2.e-7_wp,3.e-7_wp,4.e-7_wp,5.e-7_wp,6.e-7_wp,7.e-7_wp,8.e-7_wp,9.e-7_wp, &
    1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp] !! mass bins for ice water \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_r), parameter :: &
    r_r = [1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp,2.e-3_wp,3.e-3_wp,4.e-3_wp,5.e-3_wp,6.e-3_wp,7.e-3_wp,8.e-3_wp,9.e-3_wp, &
    1.e-2_wp] !! mass bins for rain \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_s), parameter :: &
    r_s = [1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp,2.e-3_wp,3.e-3_wp,4.e-3_wp,5.e-3_wp,6.e-3_wp,7.e-3_wp,8.e-3_wp,9.e-3_wp, &
    1.e-2_wp] !! mass bins for snow \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_g), parameter :: &
    r_g = [1.e-6_wp,2.e-6_wp,3.e-6_wp,4.e-6_wp,5.e-6_wp,6.e-6_wp,7.e-6_wp,8.e-6_wp,9.e-6_wp, &
    1.e-5_wp,2.e-5_wp,3.e-5_wp,4.e-5_wp,5.e-5_wp,6.e-5_wp,7.e-5_wp,8.e-5_wp,9.e-5_wp, &
    1.e-4_wp,2.e-4_wp,3.e-4_wp,4.e-4_wp,5.e-4_wp,6.e-4_wp,7.e-4_wp,8.e-4_wp,9.e-4_wp, &
    1.e-3_wp,2.e-3_wp,3.e-3_wp,4.e-3_wp,5.e-3_wp,6.e-3_wp,7.e-3_wp,8.e-3_wp,9.e-3_wp, &
    1.e-2_wp] !! mass bins for graupel \([kg\, m^{-3}]\)

  real(wp), dimension(ntb_r1), parameter :: &
    n0r_exp = [1.e6_wp,2.e6_wp,3.e6_wp,4.e6_wp,5.e6_wp,6.e6_wp,7.e6_wp,8.e6_wp,9.e6_wp, &
    1.e7_wp,2.e7_wp,3.e7_wp,4.e7_wp,5.e7_wp,6.e7_wp,7.e7_wp,8.e7_wp,9.e7_wp, &
    1.e8_wp,2.e8_wp,3.e8_wp,4.e8_wp,5.e8_wp,6.e8_wp,7.e8_wp,8.e8_wp,9.e8_wp, &
    1.e9_wp,2.e9_wp,3.e9_wp,4.e9_wp,5.e9_wp,6.e9_wp,7.e9_wp,8.e9_wp,9.e9_wp, &
    1.e10_wp] !! y-intercept bins for rain \([m^{-4}]\)

  real(wp), dimension(ntb_g1), parameter :: &
    n0g_exp = [1.e2_wp,2.e2_wp,3.e2_wp,4.e2_wp,5.e2_wp,6.e2_wp,7.e2_wp,8.e2_wp,9.e2_wp, &
    1.e3_wp,2.e3_wp,3.e3_wp,4.e3_wp,5.e3_wp,6.e3_wp,7.e3_wp,8.e3_wp,9.e3_wp, &
    1.e4_wp,2.e4_wp,3.e4_wp,4.e4_wp,5.e4_wp,6.e4_wp,7.e4_wp,8.e4_wp,9.e4_wp, &
    1.e5_wp,2.e5_wp,3.e5_wp,4.e5_wp,5.e5_wp,6.e5_wp,7.e5_wp,8.e5_wp,9.e5_wp, &
    1.e6_wp] !! y-intercept bins for graupel \([m^{-4}]\)

  real(wp), dimension(ntb_i1), parameter :: &
    nt_i = [1.0_wp,2.0_wp,3.0_wp,4.0_wp,5.0_wp,6.0_wp,7.0_wp,8.0_wp,9.0_wp, &
    1.e1_wp,2.e1_wp,3.e1_wp,4.e1_wp,5.e1_wp,6.e1_wp,7.e1_wp,8.e1_wp,9.e1_wp, &
    1.e2_wp,2.e2_wp,3.e2_wp,4.e2_wp,5.e2_wp,6.e2_wp,7.e2_wp,8.e2_wp,9.e2_wp, &
    1.e3_wp,2.e3_wp,3.e3_wp,4.e3_wp,5.e3_wp,6.e3_wp,7.e3_wp,8.e3_wp,9.e3_wp, &
    1.e4_wp,2.e4_wp,3.e4_wp,4.e4_wp,5.e4_wp,6.e4_wp,7.e4_wp,8.e4_wp,9.e4_wp, &
    1.e5_wp,2.e5_wp,3.e5_wp,4.e5_wp,5.e5_wp,6.e5_wp,7.e5_wp,8.e5_wp,9.e5_wp, &
    1.e6_wp] !! number bins for ice \([m^{-3}]\)

  real(wp), dimension(ntb_in), parameter :: &
    nt_in = [1.0_wp,2.0_wp,3.0_wp,4.0_wp,5.0_wp,6.0_wp,7.0_wp,8.0_wp,9.0_wp, &
    1.e1_wp,2.e1_wp,3.e1_wp,4.e1_wp,5.e1_wp,6.e1_wp,7.e1_wp,8.e1_wp,9.e1_wp, &
    1.e2_wp,2.e2_wp,3.e2_wp,4.e2_wp,5.e2_wp,6.e2_wp,7.e2_wp,8.e2_wp,9.e2_wp, &
    1.e3_wp,2.e3_wp,3.e3_wp,4.e3_wp,5.e3_wp,6.e3_wp,7.e3_wp,8.e3_wp,9.e3_wp, &
    1.e4_wp,2.e4_wp,3.e4_wp,4.e4_wp,5.e4_wp,6.e4_wp,7.e4_wp,8.e4_wp,9.e4_wp, &
    1.e5_wp,2.e5_wp,3.e5_wp,4.e5_wp,5.e5_wp,6.e5_wp,7.e5_wp,8.e5_wp,9.e5_wp, &
    1.e6_wp] !! number bins for IN concentration from \(0.001-1000\, L^{-1}\) \([m^{-3}]\)

  real(wp), dimension(ntb_arc), parameter :: &
    ta_na = [10._wp, 31.6_wp, 100._wp, 316._wp, &
      1000._wp, 3160._wp, 10000._wp] !! aerosol lookup table bins for number concentration
  real(wp), dimension(ntb_arw), parameter :: &
    ta_ww = [0.01_wp, 0.0316_wp, 0.1_wp, 0.316_wp, &
      1._wp, 3.16_wp, 10._wp, 31.6_wp, 100._wp] !! aerosol lookup table bins for vertical velocity
  real(wp), dimension(ntb_art), parameter :: &
    ta_tk = [243.15_wp, 253.15_wp, 263.15_wp, &
      273.15_wp, 283.15_wp, 293.15_wp, 303.15_wp] !! aerosol lookup table bins for temperature
  real(wp), dimension(ntb_arr), parameter :: &
    ta_ra = [0.01_wp, 0.02_wp, 0.04_wp, 0.08_wp, 0.16_wp] !! aerosol lookup table bins for radius
  real(wp), dimension(ntb_ark), parameter :: &
    ta_ka = [0.2_wp, 0.4_wp, 0.6_wp, 0.8_wp] !! aerosol lookup table bins for hygroscopicity

  real(wp), dimension(10), parameter :: &
    sa = [5.065339_wp, -0.062659_wp, -3.032362_wp, 0.029469_wp, -0.000285_wp, &
    0.31255_wp, 0.000204_wp, 0.003199_wp, 0._wp, -0.015952_wp] !! snow moment data from [Field et al. (2005)](https://doi.org/10.1256/qj.04.134)
  real(wp), dimension(10), parameter :: &
    sb = [0.476221_wp, -0.015896_wp, 0.165977_wp, 0.007468_wp, -0.000141_wp, &
    0.060366_wp, 0.000079_wp, 0.000594_wp, 0._wp, -0.003577_wp] !! snow moment data from [Field et al. (2005)](https://doi.org/10.1256/qj.04.134)

  real(wp), dimension(ntb_t), parameter :: &
    tc = [-0.01_wp, -5._wp, -10._wp, -15._wp, &
      -20._wp, -25._wp, -30._wp, -35._wp, -40._wp] !! temperature lookup table data

  ! variables ---------------------------------------------------------------------------------------------
  integer, protected :: dim_nrhg !! number of dimensions for graupel density

    !> @note
    !> av_g contains graupel fallspeed power-law coefficients (hail_aware = true)
    !>
    !> av_g and bv_g values from A. Heymsfield: Best - Reynolds relationship
    !> @endnote
  real(wp), protected, dimension(nrhg) :: av_g = [45.9173813_wp, 67.0867386_wp, 98.0158463_wp, &
    122.353378_wp, 143.204224_wp, 161.794724_wp, &
    178.762115_wp, 194.488785_wp, 209.225876_wp]
  real(wp), protected, dimension(nrhg) :: bv_g = [0.640961647_wp, 0.640961647_wp, 0.640961647_wp, &
    0.640961647_wp, 0.640961647_wp, 0.640961647_wp, &
    0.640961647_wp, 0.640961647_wp, 0.640961647_wp] !! graupel fallspeed power-law coefficients (hail_aware = true)

  real(wp), protected :: am_i !! ice mass-diameter power-law coefficient
  real(wp), protected :: am_r !! rain mass-diameter power-law coefficient
  real(wp), protected, dimension (nrhg) :: am_g !! graupel mass-diameter power-law coefficient
  real(wp), protected :: lfus !! enthalpy of fusion \([J\, kg^{-1}]\)
  real(wp), protected :: olfus !! 1 / lfus \([kg\, J^{-1}]\)
  real(wp), protected :: orv !! 1 / rv \([K\, kg\, J^{-1}]\)
  real(wp), protected :: ar_volume !! volume for Koop nucleation
  real(wp), protected :: sc3 !! schmidt number to the 1/3 power
  real(wp), protected :: d0i !! minimum diameter of cloud ice \([m]\)
  real(wp), protected :: xm0s !! minimum mass of snow \([kg]\)
  real(wp), protected :: xm0g !! minimum mass of graupel \([kg]\)
  real(wp), protected :: obmi !! 1 / bm_i
  real(wp), protected :: obmr !! 1 / bm_r
  real(wp), protected :: oams !! 1 / am_s
  real(wp), protected :: obms !! 1 / bm_s
  real(wp), protected :: ocms !! oams ^ obms
  real(wp), protected, dimension(nrhg) :: oamg !! 1 / am_g
  real(wp), protected, dimension(nrhg) :: ocmg !! oamg ^ obmg
  real(wp), protected :: obmg !! 1 / bm_g

  ! various gamma calculations used throughout the microphysics
  real(wp), protected, dimension(5,15) :: cce, ccg !! for \(ccg = \Gamma(x)\), cce is x for cloud water
  real(wp), protected, dimension(15) :: ocg1, ocg2 !! inverse of specific ccg values
  real(wp), protected, dimension(7) :: cie, cig !! for \(cig = \Gamma(x)\), cie is x for cloud ice
  real(wp), protected :: oig1, oig2 !! inverse of specific cig values
  real(wp), protected, dimension(13) :: cre, crg !! for \(crg = \Gamma(x)\), cre is x for rain
  real(wp), protected :: ore1, org1, org2, org3 !! inverse of specific cre and crg values
  real(wp), protected, dimension(17) :: cse, csg !! for \(csg = \Gamma(x)\), cse is x for snow
  real(wp), protected, dimension(12,nrhg) :: cge, cgg !! for \(cgg = \Gamma(x)\), cge is x for graupel
  real(wp), protected :: oge1, ogg1, ogg2, ogg3 !! inverse of specific cge and cgg values

  ! precomputed constants in various rate equations
  real(wp), protected :: t1_qr_qc, t1_qr_qi, t2_qr_qi !! terms for rain collecting cloud water and cloud ice equations
  real(wp), protected :: t1_qs_qc, t1_qs_qi !! terms for snow collecting cloud water and cloud ice equations
  real(wp), protected :: t1_qr_ev, t2_qr_ev !! terms for rain evaporation equation
  real(wp), protected :: t1_qs_sd, t2_qs_sd !! terms for deposition/sublimation of snow equation
  real(wp), protected :: t1_qs_me, t2_qs_me !! terms for melting snow equation
  real(wp), protected :: t1_qg_sd !! term for deposition/sublimation of graupel equation
  real(wp), protected :: t1_qg_me !! term for melting graupel equation

  integer :: nic2, nii2, nii3, nir2, nir3, nis2, nig2, nig3, niin2 !! lookup table indexes
  !> @history
  !> nic1 is used for cloud droplet number concentration lookup table
  !>
  !> nic1 was changed from integer in previous code versions
  !> @endhistory
  real(dp) :: nic1

  real(dp), protected, dimension(nbc) :: dc, dtc !! diameter and bin space for cloud water bins \([m]\)
  real(dp), protected, dimension(nbi) :: di, dti !! diameter and bin space for ice bins \([m]\)
  real(dp), protected, dimension(nbr) :: dr, dtr !! diameter and bin space for rain bins \([m]\)
  real(dp), protected, dimension(nbs) :: ds, dts !! diameter and bin space for snow bins \([m]\)
  real(dp), protected, dimension(nbg) :: dg, dtg !! diameter and bin space for graupel bins \([m]\)
  real(dp), protected, dimension(nbc) :: t_nc !! cloud droplet number concentration bins \([cm^{-3}]\)

  integer, parameter :: nhbins = 50 !! used for hail size calculation
  real(dp), protected, dimension(:), allocatable :: hbins, dhbins !! diameter and bin space for hail size bins \([m]\)

  integer, parameter :: radar_bins = 50 !! used for radar caculation
  real(dp), protected, dimension(:), allocatable :: sbins_radar, dsbins_radar !! diameter and bin space for snow used in radar calculation \([m]\)
  real(dp), protected, dimension(:), allocatable :: gbins_radar, dgbins_radar !! diameter and bin space for graupel used in radar calculation \([m]\)

  ! lookup table data set in module_mp_tempo_init and cannot be protected
  real(dp), allocatable, dimension(:,:) :: t_efrw, t_efsw !! collection efficiency data arrays
  real(dp), allocatable, dimension(:,:,:) :: tpc_wev, tnc_wev !! evaporation data arrays
  real(table_sp), allocatable, dimension(:,:,:,:,:) :: tnccn_act !! cloud condensation nuclei data arrays
  real(table_dp), allocatable, dimension(:,:,:,:,:) :: tcg_racg, tmr_racg, tcr_gacr, &
    tnr_racg, tnr_gacr !! rain-graupel collection data arrays
  real(table_dp), allocatable, dimension(:,:,:,:) :: tcs_racs1, tmr_racs1, tcs_racs2, &
    tmr_racs2, tcr_sacr1, tms_sacr1, tcr_sacr2, tms_sacr2, &
    tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2 !! rain-snow collection data arrays
  real(table_dp), allocatable, dimension(:,:,:,:) :: tpi_qcfz, tni_qcfz !! cloud droplet freezing data arrays
  real(table_dp), allocatable, dimension(:,:,:,:) :: tpi_qrfz, tpg_qrfz, tni_qrfz, tnr_qrfz !! rain freezing data arrays
  real(dp), allocatable, dimension(:,:) :: tps_iaus, tni_iaus, tpi_ide !! cloud ice depositional growth and conversion to snow data array
    
  ! -------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------
  contains

  subroutine get_version(version)
    !! returns the tempo version string from the README.md file
    !! or returns empty string if not found
  
    character(len=*), intent(inout) :: version
    character(len=100) :: first_line, filename
    integer :: io_unit
    logical :: fileexists

    filename = 'README.md'
    inquire(file=trim(filename), exist=fileexists)
    if (.not. fileexists) then
      version = ''
      ! write(*,'(A)') 'Unable to determine TEMPO Microphysics Version'
      return
    endif
  
    open(newunit=io_unit, file=filename, status='old', action='read')
    read(io_unit, '(A)') first_line
    close(io_unit)

    ! format is tempo-vX.X.X
    version = trim(first_line(8:))
    write(*,'(A)') 'TEMPO Microphysics Version: '//trim(version)
  end subroutine get_version


  subroutine initialize_graupel_vars(hail_flag)
    !! initialize graupel variables based on hail-aware configuration flag

    logical, intent(in) :: hail_flag

    if (hail_flag) then
      ! in case these were previously set to av_g_old and bv_g_old reset
      av_g(idx_bg1) = av_g_new
      bv_g(idx_bg1) = bv_g_new
      dim_nrhg = nrhg
    else
      av_g(idx_bg1) = av_g_old
      bv_g(idx_bg1) = bv_g_old
      dim_nrhg = nrhg1
    endif
  end subroutine initialize_graupel_vars


  subroutine initialize_parameters()
    !! initialize tempo parameters and variables
    
    integer :: m, n

    ! pi could be set by a host model, thus these parameters need to be calculated here
    am_i = pi * rho_i / 6.0_wp
    am_r = pi * rho_w / 6.0_wp
    am_g = [pi * rho_g(1) / 6.0_wp, &
      pi * rho_g(2) / 6.0_wp, &
      pi * rho_g(3) / 6.0_wp, &
      pi * rho_g(4) / 6.0_wp, &
      pi * rho_g(5) / 6.0_wp, &
      pi * rho_g(6) / 6.0_wp, &
      pi * rho_g(7) / 6.0_wp, &
      pi * rho_g(8) / 6.0_wp, &
      pi * rho_g(9) / 6.0_wp]
      ! av_i = av_s * d0s ** (bv_s - bv_i)
    ar_volume = 4.0_wp / 3.0_wp * pi * (2.5e-6_wp)**3

    lfus = lsub - lvap0
    olfus = 1.0_wp / lfus
    orv = 1.0_wp / rv
    rho_not = 101325.0_wp / (rdry*298.0_wp)
    rho_not0 = 101325.0_wp / (rdry*t0)

    ! Schmidt number to one-third used numerous times
    sc3 = sc**(1.0_wp/3.0_wp)

    ! compute minimum ice diameter from mass and minimum snow/graupel mass from diameter
    d0i = (xm0i/am_i)**(1.0_wp/bm_i)

    ! pre-compute various constants used in the microphysics equations
    xm0s = am_s * d0s**bm_s
    xm0g = am_g(nrhg) * d0g**bm_g
    obmi = 1.0_wp / bm_i
    obmr = 1.0_wp / bm_r
    oams = 1.0_wp / am_s
    obms = 1.0_wp / bm_s
    ocms = oams**obms
    obmg = 1.0_wp / bm_g
    do m = 1, nrhg
      oamg(m) = 1.0_wp / am_g(m)
      ocmg(m) = oamg(m)**obmg
    enddo

    ! gamma functions for cloud water
    do n = 1, 15
      cce(1,n) = n + 1._wp
      cce(2,n) = bm_r + n + 1._wp
      cce(3,n) = bm_r + n + 4._wp
      cce(4,n) = n + bv_c + 1._wp
      cce(5,n) = bm_r + n + bv_c + 1._wp
      ccg(1,n) = gamma(cce(1,n))
      ccg(2,n) = gamma(cce(2,n))
      ccg(3,n) = gamma(cce(3,n))
      ccg(4,n) = gamma(cce(4,n))
      ccg(5,n) = gamma(cce(5,n))
      ocg1(n) = 1.0_wp / ccg(1,n)
      ocg2(n) = 1.0_wp / ccg(2,n)
    enddo

    ! gamma functions for cloud ice
    cie(1) = mu_i + 1._wp
    cie(2) = bm_i + mu_i + 1._wp
    cie(3) = bm_i + mu_i + bv_i + 1._wp
    cie(4) = mu_i + bv_i + 1._wp
    cie(5) = mu_i + 2._wp
    cie(6) = bm_i*0.5_wp + mu_i + bv_i + 1._wp
    cie(7) = bm_i*0.5_wp + mu_i + 1._wp
    cig(1) = gamma(cie(1))
    cig(2) = gamma(cie(2))
    cig(3) = gamma(cie(3))
    cig(4) = gamma(cie(4))
    cig(5) = gamma(cie(5))
    cig(6) = gamma(cie(6))
    cig(7) = gamma(cie(7))
    oig1 = 1.0_wp / cig(1)
    oig2 = 1.0_wp / cig(2)

    ! gamma functions for rain
    cre(1) = bm_r + 1._wp
    cre(2) = mu_r + 1._wp
    cre(3) = bm_r + mu_r + 1._wp
    cre(4) = bm_r*2._wp + mu_r + 1._wp
    cre(5) = mu_r + bv_r + 1._wp
    cre(6) = bm_r + mu_r + bv_r + 1._wp
    cre(7) = bm_r*0.5_wp + mu_r + bv_r + 1._wp
    cre(8) = bm_r + mu_r + bv_r + 3._wp
    cre(9) = mu_r + bv_r + 3._wp
    cre(10) = mu_r + 2._wp
    cre(11) = 0.5_wp*(bv_r + 5._wp + 2._wp*mu_r)
    cre(12) = bm_r*0.5_wp + mu_r + 1._wp
    cre(13) = bm_r*2._wp + mu_r + bv_r + 1._wp

    do n = 1, 13
      crg(n) = gamma(cre(n))
    enddo

    ore1 = 1.0_wp / cre(1)
    org1 = 1.0_wp / crg(1)
    org2 = 1.0_wp / crg(2)
    org3 = 1.0_wp / crg(3)

    ! gamma functions for snow
    cse(1) = bm_s + 1._wp
    cse(2) = bm_s + 2._wp
    cse(3) = bm_s*2._wp
    cse(4) = bm_s + bv_s + 1._wp
    cse(5) = bm_s*2._wp + bv_s + 1._wp
    cse(6) = bm_s*2._wp + 1._wp
    cse(7) = bm_s + mu_s + 1._wp
    cse(8) = bm_s + mu_s + 2._wp
    cse(9) = bm_s + mu_s + 3._wp
    cse(10) = bm_s + mu_s + bv_s + 1._wp
    cse(11) = bm_s*2._wp + mu_s + bv_s + 1._wp
    cse(12) = bm_s*2._wp + mu_s + 1._wp
    cse(13) = bv_s + 2._wp
    cse(14) = bm_s + bv_s
    cse(15) = mu_s + 1._wp
    cse(16) = 1.0_wp + (1.0_wp + bv_s)/2._wp
    cse(17) = bm_s + bv_s + 2._wp

    do n = 1, 17
      csg(n) = gamma(cse(n))
    enddo

    ! gamma functions for graupel
    cge(1,:) = bm_g + 1._wp
    cge(2,:) = mu_g + 1._wp
    cge(3,:) = bm_g + mu_g + 1._wp
    cge(4,:) = bm_g*2. + mu_g + 1._wp
    cge(10,:) = mu_g + 2._wp
    cge(12,:) = bm_g*0.5_wp + mu_g + 1._wp

    do m = 1, nrhg
      cge(5,m) = bm_g*2._wp + mu_g + bv_g(m) + 1._wp
      cge(6,m) = bm_g + mu_g + bv_g(m) + 1._wp
      cge(7,m) = bm_g*0.5_wp + mu_g + bv_g(m) + 1._wp
      cge(8,m) = mu_g + bv_g(m) + 1._wp
      cge(9,m) = mu_g + bv_g(m) + 3._wp
      cge(11,m) = 0.5_wp*(bv_g(m) + 5._wp + 2._wp*mu_g)
    enddo

    do m = 1, nrhg
      do n = 1, 12
        cgg(n,m) = gamma(cge(n,m))
      enddo
    enddo
    oge1 = 1.0_wp / cge(1,1)
    ogg1 = 1.0_wp / cgg(1,1)
    ogg2 = 1.0_wp / cgg(2,1)
    ogg3 = 1.0_wp / cgg(3,1)

    ! rain collecting cloud water and cloud ice
    t1_qr_qc = pi * 0.25_wp * av_r * crg(9)
    t1_qr_qi = pi * 0.25_wp * av_r * crg(9)
    t2_qr_qi = pi * 0.25_wp * am_r*av_r * crg(8)

    ! snow collecting cloud water and cloud ice
    t1_qs_qc = pi * 0.25_wp * av_s
    t1_qs_qi = pi * 0.25_wp * av_s

    ! evaporation of rain; ignore depositional growth of rain.
    t1_qr_ev = 0.78_wp * crg(10)
    t2_qr_ev = 0.308_wp * sc3 * sqrt(av_r) * crg(11)

    ! sublimation/depositional growth of snow
    t1_qs_sd = 0.86_wp
    t2_qs_sd = 0.28_wp * sc3 * sqrt(av_s)

    ! melting of snow
    t1_qs_me = pi * 4._wp *c_sqrd * olfus * 0.86_wp
    t2_qs_me = pi * 4._wp *c_sqrd * olfus * 0.28_wp * sc3 * sqrt(av_s)

    ! sublimation/depositional growth of graupel
    t1_qg_sd = 0.86_wp * cgg(10,1)

    ! melting of graupel
    t1_qg_me = pi * 4._wp * c_cube * olfus * 0.86_wp * cgg(10,1)
  end subroutine initialize_parameters
    

  subroutine initialize_bins_for_tables()
    !! initialize log-spaced bins of hydrometer quantities used for lookup tables

    integer :: n

    ! constants for helping find lookup table indexes
    nic2 = nint(log10(r_c(1)))
    nii2 = nint(log10(r_i(1)))
    nir2 = nint(log10(r_r(1)))
    nis2 = nint(log10(r_s(1)))
    nig2 = nint(log10(r_g(1)))
    nii3 = nint(log10(nt_i(1)))
    nir3 = nint(log10(n0r_exp(1)))
    nig3 = nint(log10(n0g_exp(1)))
    niin2 = nint(log10(nt_in(1)))

    ! bins of cloud water (from minimum diameter to 100 microns)
    dc(1) = real(d0c, kind=dp)
    dtc(1) = real(d0c, kind=dp)
    do n = 2, nbc
      dc(n) = dc(n-1) + 1.0e-6_dp
      dtc(n) = (dc(n) - dc(n-1))
    enddo

    ! bins of cloud ice (from min diameter up to 2x min snow size)
    call create_bins(numbins=nbi, lowbin=real(d0i, kind=dp), &
      highbin=2.0_dp*d0s, bins=di, deltabins=dti)

    ! bins of rain (from min diameter up to 5 mm)
    call create_bins(numbins=nbr, lowbin=real(d0r, kind=dp), &
      highbin=0.005_dp, bins=dr, deltabins=dtr)

    ! bins of snow (from min diameter up to 2 cm)
    call create_bins(numbins=nbs, lowbin=real(d0s, kind=dp), &
      highbin=0.02_dp, bins=ds, deltabins=dts)

    ! bins of graupel (from min diameter up to 5 cm)
    call create_bins(numbins=nbg, lowbin=real(d0g, kind=dp), &
      highbin=0.05_dp, bins=dg, deltabins=dtg)

    ! bins of cloud droplet number concentration (1 to 3000 per cc)
    call create_bins(numbins=nbc, lowbin=1.0_dp, &
      highbin=3000.0_dp, bins=t_nc)
    t_nc = t_nc * 1.0e6_dp
    nic1 = real(log(t_nc(nbc)/t_nc(1)), kind=dp)
  end subroutine initialize_bins_for_tables


  subroutine initialize_bins_for_hail_size()
    !! initialize log-spaced bins for hail size calculation

    real(dp), parameter :: lowbin = 500.e-6_dp
    real(dp), parameter :: highbin = 0.075_dp
  
    if (.not. allocated(hbins)) allocate(hbins(nhbins), source=0._dp)
    if (.not. allocated(dhbins)) allocate(dhbins(nhbins), source=0._dp)
    call create_bins(numbins=nhbins, lowbin=lowbin, highbin=highbin, &
      bins=hbins, deltabins=dhbins)
  end subroutine initialize_bins_for_hail_size


  subroutine initialize_bins_for_radar()
    !! initialize log-spaced bins for radar calculation

    real(dp), parameter :: lowbin = 100.e-6_dp
    real(dp), parameter :: s_highbin = 0.02_dp
    real(dp), parameter :: g_highbin = 0.05_dp
  
    if (.not. allocated(sbins_radar)) allocate(sbins_radar(radar_bins), source=0._dp)
    if (.not. allocated(dsbins_radar)) allocate(dsbins_radar(radar_bins), source=0._dp)
    ! bins of snow (from 100 microns up to 2 cm)
    call create_bins(numbins=radar_bins, lowbin=lowbin, &
      highbin=s_highbin, bins=sbins_radar, deltabins=dsbins_radar)

    if (.not. allocated(gbins_radar)) allocate(gbins_radar(radar_bins), source=0._dp)
    if (.not. allocated(dgbins_radar)) allocate(dgbins_radar(radar_bins), source=0._dp)
      ! bins of graupel (from 100 microns up to 5 cm)
    call create_bins(numbins=radar_bins, lowbin=lowbin, &
      highbin=g_highbin, bins=gbins_radar, deltabins=dgbins_radar)
  end subroutine initialize_bins_for_radar


  subroutine create_bins(numbins, lowbin, highbin, bins, deltabins)
    !! calculates log-spaced bins of hydrometer sizes to simplify calculations later
  
    integer, intent(in) :: numbins
    real(dp), intent(in) :: lowbin, highbin

    real(dp), dimension(:), intent(out) :: bins
    real(dp), dimension(:), intent(out), optional :: deltabins
    
    integer :: n
    real(dp), dimension(numbins+1) :: xdx
  
    xdx(1) = lowbin
    xdx(numbins+1) = highbin
    do  n = 2, numbins
      xdx(n) = exp(real(n-1, kind=dp)/real(numbins, kind=dp) * log(xdx(numbins+1)/xdx(1)) + log(xdx(1)))
    enddo

    do n = 1, numbins
      bins(n) = sqrt(xdx(n)*xdx(n+1))
    enddo

    if (present(deltabins)) then
      do n = 1, numbins
        deltabins(n) = xdx(n+1) - xdx(n)
      enddo
    endif
  end subroutine create_bins


  subroutine initialize_arrays_freezewater(table_size)
    !! initialize data arrays for [Bigg (1953)](https://doi.org/10.1002/qj.49707934207) freezing of cloud water and rain
    !! @warning
    !! six lookup table arrays are allocated and set here
    !! @endwarning

    integer, intent(out), optional :: table_size

    ! cloud water freezing
    if (.not. allocated(tpi_qcfz)) allocate(tpi_qcfz(ntb_c,nbc,ntb_t1,ntb_in), source=0._table_dp)
    if (.not. allocated(tni_qcfz)) allocate(tni_qcfz(ntb_c,nbc,ntb_t1,ntb_in), source=0._table_dp)

    ! rain freezing
    if (.not. allocated(tpi_qrfz)) allocate(tpi_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_in), source=0._table_dp)
    if (.not. allocated(tpg_qrfz)) allocate(tpg_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_in), source=0._table_dp)
    if (.not. allocated(tni_qrfz)) allocate(tni_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_in), source=0._table_dp)
    if (.not. allocated(tnr_qrfz)) allocate(tnr_qrfz(ntb_r,ntb_r1,ntb_t1,ntb_in), source=0._table_dp)

    ! table size is precision * entries * dimensions
    if (present(table_size)) then
      table_size = (table_dp * 2 * (ntb_c*nbc*ntb_t1*ntb_in)) + &
        (table_dp * 4 * (ntb_r*ntb_r1*ntb_t1*ntb_in))
    endif 
  end subroutine initialize_arrays_freezewater


  subroutine initialize_arrays_qr_acr_qs(table_size)
    !! initialize data arrays for rain-snow collection
    !! @warning
    !! twelve lookup table arrays are allocated and set here
    !! @endwarning

    integer, intent(out), optional :: table_size

    if (.not. allocated(tcs_racs1)) allocate(tcs_racs1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tmr_racs1)) allocate(tmr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tcs_racs2)) allocate(tcs_racs2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tmr_racs2)) allocate(tmr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tcr_sacr1)) allocate(tcr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tms_sacr1)) allocate(tms_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tcr_sacr2)) allocate(tcr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tms_sacr2)) allocate(tms_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_racs1)) allocate(tnr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_racs2)) allocate(tnr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_sacr1)) allocate(tnr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_sacr2)) allocate(tnr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r), source=0._table_dp)

    ! table size is precision * entries * dimensions
    if (present(table_size)) then
      table_size = table_dp * 12 * (ntb_s*ntb_t*ntb_r1*ntb_r)
    endif 
  end subroutine initialize_arrays_qr_acr_qs


  subroutine initialize_arrays_qr_acr_qg(table_size)
    !! initialize data arrays for rain-graupel collection
    !! @warning
    !! five lookup table arrays are allocated and set here
    !! @endwarning

    integer, intent(out), optional :: table_size

     ! rain-graupel
    if (.not. allocated(tcg_racg)) allocate(tcg_racg(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tmr_racg)) allocate(tmr_racg(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tcr_gacr)) allocate(tcr_gacr(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_racg)) allocate(tnr_racg(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)
    if (.not. allocated(tnr_gacr)) allocate(tnr_gacr(ntb_g1,ntb_g,nrhg,ntb_r1,ntb_r), source=0._table_dp)

    ! table size is precision * entries * dimensions 
    if (present(table_size)) then
      table_size = table_dp * 5 * (ntb_g1*ntb_g*nrhg*ntb_r1*ntb_r)
    endif 
  end subroutine initialize_arrays_qr_acr_qg


  subroutine initialize_arrays_ccn(table_size)
    !! initialize data arrays for ccn lookup table
    !! @warning
    !! one lookup table array is allocated and set here
    !! @endwarning

    integer, intent(out), optional :: table_size

    if (.not. allocated(tnccn_act)) &
      allocate(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark), source=0._table_sp)

    ! table size is precision * entries * dimensions
    if (present(table_size)) then
      table_size = table_sp * 1 * (ntb_arc*ntb_arw*ntb_art*ntb_arr*ntb_ark) + (table_sp + table_sp) * 1
    endif
  end subroutine initialize_arrays_ccn


  subroutine initialize_arrays_drop_evap()
    !! initialize data arrays for drop evaporation data
    !! @warning
    !! two lookup table arrays are allocated and set here
    !! @endwarning

      if (.not. allocated(tpc_wev)) allocate(tpc_wev(nbc,ntb_c,nbc), source=0._dp)
      if (.not. allocated(tnc_wev)) allocate(tnc_wev(nbc,ntb_c,nbc), source=0._dp)
  end subroutine initialize_arrays_drop_evap

  
  subroutine initialize_array_efsw()
    !! initializes the collision efficiency data array for snow collecting cloud water
    !! @warning
    !! one lookup table array is allocated and set here
    !! @endwarning

    if (.not. allocated(t_efsw)) allocate(t_efsw(nbs,nbc), source=0._dp)
  end subroutine initialize_array_efsw


  subroutine initialize_array_efrw()
    !! initializes the collision efficiency data array for rain collecting cloud water
    !! @warning
    !! one lookup table array is allocated and set here
    !! @endwarning

    if (.not. allocated(t_efrw)) allocate(t_efrw(nbr,nbc), source=0._dp)
  end subroutine initialize_array_efrw


  subroutine initialize_arrays_qi_aut_qs()
    !! initializes data arrays for cloud ice to snow conversion and growth
    !! @warning
    !! three lookup table arrays are allocated and set here
    !! @endwarning

    if (.not. allocated(tps_iaus)) allocate(tps_iaus(ntb_i,ntb_i1), source=0._dp)
    if (.not. allocated(tni_iaus)) allocate(tni_iaus(ntb_i,ntb_i1), source=0._dp)
    if (.not. allocated(tpi_ide)) allocate(tpi_ide(ntb_i,ntb_i1), source=0._dp)
  end subroutine initialize_arrays_qi_aut_qs

end module module_mp_tempo_params

