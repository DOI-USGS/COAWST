module module_mp_tempo_aerosols
  !! contains produces used when aerosol-aware = true
  use module_mp_tempo_params, only : wp, sp, dp, &
    naccn0, naccn1, nain0, nain1, nwfa_default, aero_max

  implicit none
  private

  public :: init_water_friendly_aerosols, init_ice_friendly_aerosols, &
    aerosol_collection_efficiency

  contains

  subroutine init_water_friendly_aerosols(dz1d, nwfa)
    !! sets water-friendly aerosols to an exponential profile
    !! if aerosol-aware = true and no initial condition
    !! is provided by the host model
    real(wp), dimension(:), intent(in) :: dz1d
    real(wp), dimension(:), intent(inout) :: nwfa
    real(wp) :: hgt(size(dz1d))
    real(wp) :: h_01, niccn3
    integer :: k, nz
    
    nz = size(dz1d)
    hgt = 0._wp
    do k = 2, nz
      hgt(k) = hgt(k-1) + dz1d(k)
    enddo

    if(hgt(1) <= 1000.0_wp) then
      h_01 = 0.8_wp
    elseif(hgt(1) >= 2500.0_wp) then
      h_01 = 0.01_wp
    else
      h_01 = 0.8_wp*cos(hgt(1)*0.001_wp - 1.0_wp)
    endif
    niccn3 = -1.0_wp*log(naccn1/naccn0)/h_01
    nwfa(1) = naccn1+naccn0*exp(-((hgt(2)-hgt(1))/1000._wp)*niccn3)
    do k = 2, nz
      nwfa(k) = naccn1+naccn0*exp(-((hgt(k)-hgt(1))/1000._wp)*niccn3)
    enddo
  end subroutine init_water_friendly_aerosols


  subroutine init_ice_friendly_aerosols(dz1d, nifa)
    !! sets ice-friendly aerosols to an exponential profile
    !! if aerosol-aware = true and no initial condition
    !! is provided by the host model
    real(wp), dimension(:), intent(in) :: dz1d
    real(wp), dimension(:), intent(inout) :: nifa
    real(wp), dimension(:), allocatable :: hgt
    real(wp) :: h_01, niin3
    integer :: k, nz
    
    nz = size(dz1d)
    allocate(hgt(nz), source=0._wp)
    do k = 2, nz
      hgt(k) = hgt(k-1) + dz1d(k)
    enddo

    if(hgt(1) <= 1000.0_wp) then
      h_01 = 0.8_wp
    elseif(hgt(1) >= 2500.0_wp) then
      h_01 = 0.01_wp
    else
      h_01 = 0.8_wp*cos(hgt(1)*0.001_wp - 1.0_wp)
    endif
    niin3 = -1.0_wp*log(nain1/nain0)/h_01
    nifa(1) = nain1+nain0*exp(-((hgt(2)-hgt(1))/1000._wp)*niin3)
    do k = 2, nz
      nifa(k) = nain1+nain0*exp(-((hgt(k)-hgt(1))/1000._wp)*niin3)
    enddo
  end subroutine init_ice_friendly_aerosols


  function aerosol_collection_efficiency(d, da, visc, rhoa, temp, species) result(eff_a)
    !! computes aerosol collection efficiency for precipitation scavenging
    !! from [Wang et al. (2010)](https://doi.org/10.5194/acp-10-5685-2010)
    use module_mp_tempo_params, only : rho_w, rho_s, av_s, bv_s, &
      idx_bg1, pi, av_g, bv_g, rho_g

    real(dp), intent(in) :: d
    real(wp), intent(in) :: da, visc, rhoa, temp
    character(len=1), intent(in) :: species
    real(wp) :: aval, cc, diff, re, sc, st, st2, vt, eff, rho_p
    real(wp), parameter :: boltzman = 1.3806503e-23_wp
    real(wp), parameter :: meanpath = 0.0256e-6_wp
    real(wp) :: eff_a

    vt = 1._wp
    rho_p = rho_w
    ! rain
    if (species == 'r') then
      vt = -0.1021_wp + 4.932e3_wp*d - 0.9551e6_wp*d*d + &
        0.07934e9_wp*d*d*d - 0.002362e12_wp*d*d*d*d
      rho_p = rho_w
    ! snow
    elseif (species == 's') then
      vt = av_s*d**bv_s
      rho_p = rho_s
    ! graupel
    elseif (species .eq. 'g') then
      vt = av_g(idx_bg1)*d**bv_g(idx_bg1)
      rho_p = rho_g(idx_bg1)
    endif
    cc = 1._wp + 2._wp*meanpath/da *(1.257_wp+0.4_wp*exp(-0.55_wp*da/meanpath))
    diff = boltzman*temp*cc/(3._wp*pi*visc*da)
    re = 0.5_wp*rhoa*d*vt/visc
    sc = visc/(rhoa*diff)
    st = (rho_p-rhoa)*da*da*vt*cc/(9._wp*visc*d)
    aval = log(1._wp + re)
    st2 = (1.2_wp + 1._wp/12._wp*aval)/(1._wp+aval)
    eff = 4._wp/(re*sc) * (1._wp + 0.4_wp*sqrt(re)*sc**0.3333_wp + &
      0.16_wp*sqrt(re)*sqrt(sc)) + 4._wp*da/d * &
      (0.02_wp + da/d*(1._wp+2._wp*sqrt(re)))
    if (st > st2) eff = eff + ((st-st2)/(st-st2+0.666667_wp))**1.5_wp
    eff_a = max(1.e-5_wp, min(eff, 1._wp))
  end function aerosol_collection_efficiency

end module module_mp_tempo_aerosols
