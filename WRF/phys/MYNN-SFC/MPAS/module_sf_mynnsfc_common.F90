!====================================================================

 module module_sf_mynnsfc_common

!------------------------------------------
!Define Model-specific constants/parameters.
!This module will be used at the initialization stage
!where all model-specific constants are read and saved into
!memory. This module is then used again in the MYNNSFC_* scheme. All
!MYNN-specific constants are declared globally in the main
!module (module_sf_mynnsfc) further below:
!------------------------------------------
!
! For MPAS:
 use mpas_kind_types,only: kind_phys => RKIND

!WRF CONSTANTS
  use mpas_atmphys_constants, only:          &
    & karman, grav => gravity, p1000mb => P0,&
    & cp, r_d, r_v, rcp, xlv, xlf, xls,      &
    & svp1, svp2, svp3, ep_1, ep_2,          &
    & rvovrd => rvord, cpv, cliq, cice,      &
    & svpt0

!WRF CONFIGURATIONS
!  use module_state_description, only:       &
!    & ruclsmscheme
  
 implicit none
 save

! For reference:
! real:: cp           != 7.*r_d/2. (J/kg/K)
! real:: cpv          != 4.*r_v    (J/kg/K) Spec heat H2O gas
! real:: cice         != 2106.     (J/kg/K) Spec heat H2O ice
! real:: cliq         != 4190.     (J/kg/K) Spec heat H2O liq
! real:: p608         != R_v/R_d-1.
! real:: ep_2         != R_d/R_v
! real:: grav         != accel due to gravity
! real:: karman       != von Karman constant
! real:: t0c          != temperature of water at freezing, 273.15 K
! real:: rcp          != r_d/cp
! real:: r_d          != 287.  (J/kg/K) gas const dry air
! real:: r_v          != 461.6 (J/kg/K) gas const water
! real:: xlf          != 0.35E6 (J/kg) fusion at 0 C
! real:: xlv          != 2.50E6 (J/kg) vaporization at 0 C
! real:: xls          != 2.85E6 (J/kg) sublimation
! real:: rvovrd       != r_v/r_d != 1.608
! real:: p1000mb      != 100000.0
! real:: svp1         != 0.6112 !(kPa)
! real:: svp2         != 17.67  !(dimensionless)
! real:: svp3         != 29.65  !(K)
 
! Specified locally
 real(kind_phys),parameter:: tref   = 300.0   !=reference temperature (K)
 real(kind_phys),parameter:: p608   = ep_1
 real(kind_phys),parameter:: t0c    = svpt0   != 273.15
 
! To be derived in the init routine
 real(kind_phys),parameter:: ep_3   = 1.-ep_2 != 0.378
 real(kind_phys),parameter:: gtr    = grav/tref
 real(kind_phys),parameter:: rk     = cp/r_d
 real(kind_phys),parameter:: tv0    =  p608*tref
 real(kind_phys),parameter:: tv1    = (1.+p608)*tref
 real(kind_phys),parameter:: xlscp  = (xlv+xlf)/cp
 real(kind_phys),parameter:: xlvcp  = xlv/cp
 real(kind_phys),parameter:: g_inv  = 1./grav

! !input arguments:
! type(mpas_pool_type),intent(in):: configs
! !local pointers:
! character(len=StrKIND),pointer:: config_lsm_scheme
!
! call mpas_pool_get_config(configs,'config_lsm_scheme',config_lsm_scheme)

!---------------------------------------------------

end module module_sf_mynnsfc_common
