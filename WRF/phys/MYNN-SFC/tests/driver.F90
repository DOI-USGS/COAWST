! ###########################################################################################
! CI test driver for MYNN SFC scheme
! ###########################################################################################
program driver
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
  use module_sf_mynnsfc_common,only: kind_phys,cp,lsm_ruc => ruclsmscheme
  use module_sf_mynnsfc_ccpp_tests
  use module_sf_mynnsfc_wrf_tests

  implicit none

  call ccpp_test(saveoutput=.true.)
  
  call wrf_test(case='wat', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=1, saveoutput=.true.)
  call wrf_test(case='wat', sf_mynn_sfcflux_water=0, sf_mynn_sfcflux_land=1, saveoutput=.false.)
  call wrf_test(case='wat', sf_mynn_sfcflux_water=2, sf_mynn_sfcflux_land=1, saveoutput=.false.)
  call wrf_test(case='wat', sf_mynn_sfcflux_water=3, sf_mynn_sfcflux_land=1, saveoutput=.false.)
  call wrf_test(case='wat', sf_mynn_sfcflux_water=4, sf_mynn_sfcflux_land=1, saveoutput=.false.)

  call wrf_test(case='lnd', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=1, saveoutput=.true.)
  call wrf_test(case='lnd', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=0, saveoutput=.false.)
  call wrf_test(case='lnd', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=2, saveoutput=.false.)

  call wrf_test(case='icy', sf_mynn_sfcflux_water=1, sf_mynn_sfcflux_land=0, saveoutput=.true.)

end program driver
