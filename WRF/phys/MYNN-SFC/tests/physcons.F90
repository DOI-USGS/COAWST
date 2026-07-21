!> \file physcons.F90
!! This file contains module physcons

!  ==========================================================  !!!!!
!                 module  'physcons' description               !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!   this module contains some the most frequently used math and        !
!   physics constatns for gcm models.                                  !
!                                                                      !
!   references:                                                        !
!     as set in NMC handbook from Smithsonian tables.                  !
!                                                                      !
!   modification history:                                              !
!                                                                      !
!     1990-04-30  g and rd are made consistent with NWS usage          !
!     2001-10-22  g made consistent with SI usage                      !
!     2005-04-13  added molecular weights for gases          - y-t hou !
!     2013-07-12  added temperature for homogen. nuc. for ice. - R.sun !
!                                                                      !
!   external modules referenced:                                       !
!                                                                      !
!       'module machine'                    in 'machine.f'             !
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!

!> \defgroup physcons GFS Physics Constants Module
!> This module contains some of the most frequently used math and physics
!! constants for GCM models.

!> This module contains some of the most frequently used math and physics
!! constants for GCM models.
          module physcons
!
  use machine, only: kind_phys, kind_dyn
!
  implicit none
!
  public

!> \name Math constants
! real(kind=kind_phys),parameter:: con_pi     =3.1415926535897931        !< pi
  real(kind=kind_phys),parameter:: con_pi     =4.0d0*atan(1.0d0)         !< pi
  real(kind=kind_phys),parameter:: con_sqrt2  =1.414214e+0_kind_phys               !< square root of 2
  real(kind=kind_phys),parameter:: con_sqrt3  =1.732051e+0_kind_phys               !< quare root of 3

!> \name Geophysics/Astronomy constants
  real(kind=kind_phys),parameter:: con_rerth  =6.3712e+6_kind_phys                 !< radius of earth (\f$m\f$)
  real(kind=kind_phys),parameter:: con_g      =9.80665e+0_kind_phys                !< gravity (\f$m/s^{2}\f$)
  real(kind=kind_phys),parameter:: con_omega  =7.2921e-5_kind_phys                 !< ang vel of earth (\f$s^{-1}\f$)
  real(kind=kind_phys),parameter:: con_p0     =1.01325e5_kind_phys                 !< standard atmospheric pressure (\f$Pa\f$)
! real(kind=kind_phys),parameter:: con_solr   =1.36822e+3_kind_phys                ! solar constant    (W/m2)-aer(2001)
  real(kind=kind_phys),parameter:: con_solr_2002 =1.3660e+3_kind_phys              !< solar constant (\f$W/m^{2}\f$)-Liu(2002)
  real(kind=kind_phys),parameter:: con_solr_2008 =1.3608e+3_kind_phys              !< solar constant (\f$W/m^{2}\f$)-nasa-sorce Tim(2008)
! real(kind=kind_phys),parameter:: con_solr   =1.36742732e+3_kind_phys             ! solar constant    (W/m2)-gfdl(1989) - OPR as of Jan 2006
  ! Selected geophysics/astronomy constants with kind=kind_dyn
  real(kind=kind_dyn), parameter:: con_g_dyn  =9.80665e+0_kind_dyn                 !< gravity (\f$m/s^{2}\f$)

!> \name Thermodynamics constants
  real(kind=kind_phys),parameter:: con_rgas   =8.314472_kind_phys                  !< molar gas constant (\f$J/mol/K\f$)
  real(kind=kind_phys),parameter:: con_rd     =2.8705e+2_kind_phys                 !< gas constant air (\f$J/kg/K\f$)
  real(kind=kind_phys),parameter:: con_rv     =4.6150e+2_kind_phys                 !< gas constant H2O (\f$J/kg/K\f$)
  real(kind=kind_phys),parameter:: con_cp     =1.0046e+3_kind_phys                 !< spec heat air at p (\f$J/kg/K\f$)
  real(kind=kind_phys),parameter:: con_cv     =7.1760e+2_kind_phys                 !< spec heat air at v (\f$J/kg/K\f$)
  real(kind=kind_phys),parameter:: con_cvap   =1.8460e+3_kind_phys                 !< spec heat H2O gas (\f$J/kg/K\f$)
  real(kind=kind_phys),parameter:: con_cliq   =4.1855e+3_kind_phys                 !< spec heat H2O liq (\f$J/kg/K\f$)
  real(kind=kind_phys),parameter:: con_csol   =2.1060e+3_kind_phys                 !< spec heat H2O ice (\f$J/kg/K\f$)
  real(kind=kind_phys),parameter:: con_hvap   =2.5000e+6_kind_phys                 !< lat heat H2O cond (\f$J/kg\f$)
! real(kind=kind_phys),parameter:: con_hvap   =2.5010e+6_kind_phys                 ! from AMS
  real(kind=kind_phys),parameter:: con_hfus   =3.3358e+5_kind_phys                 !< lat heat H2O fusion (\f$J/kg\f$)
! real(kind=kind_phys),parameter:: con_hfus   =3.3370e+5_kind_phys                 ! from AMS
  real(kind=kind_phys),parameter:: con_psat   =6.1078e+2_kind_phys                 !< pres at H2O 3pt (\f$Pa\f$)
  real(kind=kind_phys),parameter:: con_t0c    =2.7315e+2_kind_phys                 !< temp at 0C (K)
  real(kind=kind_phys),parameter:: con_ttp    =2.7316e+2_kind_phys                 !< temp at H2O 3pt (K)
  real(kind=kind_phys),parameter:: con_tice   =2.7120e+2_kind_phys                 !< temp freezing sea (K)
  real(kind=kind_phys),parameter:: con_jcal   =4.1855E+0_kind_phys                 !< joules per calorie
  real(kind=kind_phys),parameter:: con_rhw0   =1022.0_kind_phys                    !< sea water reference density (\f$kg/m^{3}\f$)
  real(kind=kind_phys),parameter:: con_epsq   =1.0E-12_kind_phys                   !< min q for computing precip type
  real(kind=kind_phys),parameter:: con_epsqs  =1.0E-10_kind_phys
  ! Selected thermodynamics constants with kind=kind_dyn
  real(kind=kind_dyn), parameter:: con_rd_dyn   =2.8705e+2_kind_dyn                !< gas constant air (\f$J/kg/K\f$)
  real(kind=kind_dyn), parameter:: con_rv_dyn   =4.6150e+2_kind_dyn                !< gas constant H2O (\f$J/kg/K\f$)
  real(kind=kind_dyn), parameter:: con_cp_dyn   =1.0046e+3_kind_dyn                !< spec heat air at p (\f$J/kg/K\f$)
  real(kind=kind_dyn), parameter:: con_hvap_dyn =2.5000e+6_kind_dyn                !< lat heat H2O cond (\f$J/kg\f$)
  real(kind=kind_dyn), parameter:: con_hfus_dyn =3.3358e+5_kind_dyn                !< lat heat H2O fusion (\f$J/kg\f$)

!> \name Secondary constants
  real(kind=kind_phys),parameter:: con_rocp   =con_rd/con_cp
  real(kind=kind_phys),parameter:: con_cpor   =con_cp/con_rd
  real(kind=kind_phys),parameter:: con_rog    =con_rd/con_g
  real(kind=kind_phys),parameter:: con_fvirt  =con_rv/con_rd-1.
  real(kind=kind_phys),parameter:: con_eps    =con_rd/con_rv
  real(kind=kind_phys),parameter:: con_epsm1  =con_rd/con_rv-1.
  real(kind=kind_phys),parameter:: con_dldt   =con_cvap-con_cliq
  real(kind=kind_phys),parameter:: con_xpona  =-con_dldt/con_rv
  real(kind=kind_phys),parameter:: con_xponb  =-con_dldt/con_rv+con_hvap/(con_rv*con_ttp)
  real(kind=kind_phys),parameter:: con_1ovg   = 1._kind_phys/con_g

!> \name Other Physics/Chemistry constants (source: 2002 CODATA)
  real(kind=kind_phys),parameter:: con_c      =2.99792458e+8_kind_phys             !< speed of light (\f$m/s\f$)
  real(kind=kind_phys),parameter:: con_plnk   =6.6260693e-34_kind_phys             !< planck constant (\f$J/s\f$)
  real(kind=kind_phys),parameter:: con_boltz  =1.3806505e-23_kind_phys             !< boltzmann constant (\f$J/K\f$)
  real(kind=kind_phys),parameter:: con_sbc    =5.670400e-8_kind_phys               !< stefan-boltzmann (\f$W/m^{2}/K^{4}\f$)
  real(kind=kind_phys),parameter:: con_avgd   =6.0221415e23_kind_phys              !< avogadro constant (\f$mol^{-1}\f$)
  real(kind=kind_phys),parameter:: con_gasv   =22413.996e-6_kind_phys              !< vol of ideal gas at 273.15K, 101.325kPa (\f$m^{3}/mol\f$)
! real(kind=kind_phys),parameter:: con_amd    =28.970_kind_phys                    !< molecular wght of dry air (g/mol)
  real(kind=kind_phys),parameter:: con_amd    =28.9644_kind_phys                   !< molecular wght of dry air (\f$g/mol\f$)
  real(kind=kind_phys),parameter:: con_amw    =18.0154_kind_phys                   !< molecular wght of water vapor (\f$g/mol\f$)
  real(kind=kind_phys),parameter:: con_amo3   =47.9982_kind_phys                   !< molecular wght of o3 (\f$g/mol\f$)
! real(kind=kind_phys),parameter:: con_amo3   =48.0_kind_phys                      !< molecular wght of o3  (g/mol)
  real(kind=kind_phys),parameter:: con_amco2  =44.011_kind_phys                    !< molecular wght of co2 (\f$g/mol\f$)
  real(kind=kind_phys),parameter:: con_amo2   =31.9999_kind_phys                   !< molecular wght of o2 (\f$g/mol\f$)
  real(kind=kind_phys),parameter:: con_amch4  =16.043_kind_phys                    !< molecular wght of ch4 (\f$g/mol\f$)
  real(kind=kind_phys),parameter:: con_amn2o  =44.013_kind_phys                    !< molecular wght of n2o (\f$g/mol\f$)
  real(kind=kind_phys),parameter:: con_thgni  =-38.15_kind_phys                    !< temperature the H.G.Nuc. ice starts
  real(kind=kind_phys),parameter:: karman     =0.4_kind_phys                       !< Von Karman constant
  real(kind=kind_phys),parameter:: con_runiver=con_avgd*con_boltz

!> minimum ice concentration
  real(kind=kind_phys),parameter:: cimin      =0.15                                !< minimum ice concentration

!> minimum aerosol concentration
  real(kind=kind_phys),parameter:: qamin = 1.e-16_kind_phys
!> minimum rain amount
  real(kind=kind_phys),parameter:: rainmin = 1.e-13_kind_phys
!> \name Miscellaneous physics related constants (For WSM6; Moorthi - Jul 2014)
! integer, parameter :: max_lon=16000, max_lat=8000, min_lon=192, min_lat=94
! integer, parameter :: max_lon=5000,  max_lat=2500, min_lon=192, min_lat=94 ! current opr
! integer, parameter :: max_lon=5000,  max_lat=2000, min_lon=192, min_lat=94 ! current opr
! integer, parameter :: max_lon=8000,  max_lat=4000, min_lon=192, min_lat=94 ! current opr
! real(kind=kind_phys), parameter:: rlapse  = 0.65e-2, rhc_max = 0.9999      ! current opr
! real(kind=kind_phys), parameter:: rlapse  = 0.65e-2, rhc_max = 0.9999999   ! new
! real(kind=kind_phys), parameter:: rlapse  = 0.65e-2, rhc_max = 0.9900

  real(kind=kind_phys), parameter:: rlapse  = 0.65e-2_kind_phys
  real(kind=kind_phys), parameter:: cb2mb   = 10.0_kind_phys, pa2mb   = 0.01_kind_phys
! for wsm6
  real(kind=kind_phys),parameter:: rhowater   = 1000._kind_phys                    !< density of water (kg/m^3)
  real(kind=kind_phys),parameter:: rhosnow    = 100._kind_phys                     !< density of snow (kg/m^3)
  real(kind=kind_phys),parameter:: rhoair     = 1.28_kind_phys                     !< density of air near surface (kg/m^3)
  real(kind=kind_phys),parameter:: rholakeice = 0.917e3_kind_phys                  !< density of ice on lake (kg/m^3)
  real(kind=kind_phys),parameter:: rhoair_IFS = 1._kind_phys                       !< reference air density (kg/m^3), ref: IFS

! Decorrelation length constant (km) for iovr = 4 or 5 and idcor = 0
  real(kind=kind_phys),parameter:: decorr_con = 2.50_kind_phys

! for gfdlmp v3
  real(kind=kind_phys), parameter :: visd  = 1.717e-5 ! dynamics viscosity of air at 0 deg C and 1000 hPa (Mason, 1971) (kg/m/s)
  real(kind=kind_phys), parameter :: visk  = 1.35e-5  ! kinematic viscosity of air at 0 deg C  and 1000 hPa (Mason, 1971) (m^2/s)
  real(kind=kind_phys), parameter :: vdifu = 2.25e-5  ! diffusivity of water vapor in air at 0 deg C  and 1000 hPa (Mason, 1971) (m^2/s)
  real(kind=kind_phys), parameter :: tcond = 2.40e-2  ! thermal conductivity of air at 0 deg C  and 1000 hPa (Mason, 1971) (J/m/s/K)
  real(kind=kind_phys), parameter :: cdg   = 3.15121  ! drag coefficient of graupel (Locatelli and Hobbs, 1974)
  real(kind=kind_phys), parameter :: cdh   = 0.5      ! drag coefficient of hail (Heymsfield and Wright, 2014)
  real(kind=kind_phys), parameter :: rhocw = 1.0e3    ! density of cloud water (kg/m^3)
  real(kind=kind_phys), parameter :: rhoci = 9.17e2   ! density of cloud ice (kg/m^3)
  real(kind=kind_phys), parameter :: rhocr = 1.0e3    ! density of rain (Lin et al. 1983) (kg/m^3)
  real(kind=kind_phys), parameter :: rhocg = 4.0e2    ! density of graupel (Rutledge and Hobbs 1984) (kg/m^3)
  real(kind=kind_phys), parameter :: rhoch = 9.17e2   ! density of hail (Lin et al. 1983) (kg/m^3)
  real(kind=kind_phys), parameter :: qcmin = 1.0e-15  ! min value for cloud condensates (kg/kg)
  real(kind=kind_phys), parameter :: qfmin = 1.0e-8   ! min value for sedimentation (kg/kg)
  real(kind=kind_phys), parameter :: con_one      = 1_kind_phys
  real(kind=kind_phys), parameter :: con_p001     = 0.001_kind_phys
  real(kind=kind_phys), parameter :: con_secinday = 86400._kind_phys

!........................................!
      end module physcons                !
!========================================!
