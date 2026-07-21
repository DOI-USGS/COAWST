!WRF:MODEL_LAYER:CONSTANTS
!

MODULE module_model_constants

  implicit none

  !-----------------------------------------
  ! Numeric kind
  !-----------------------------------------
  integer, parameter :: wp = kind(1.0d0)

  !-----------------------------------------
  ! 2. Real number bounds
  !-----------------------------------------
  real(wp), parameter :: epsilon = 1.0e-15_wp

  !-----------------------------------------
  ! 4. Physical constants
  !-----------------------------------------

  ! JM NOTE -- can we name this grav instead?
  real(wp), parameter :: g   = 9.81_wp        ! gravity (m s^-2)

  real(wp), parameter :: r_d = 287.0_wp       ! gas constant dry air (J K^-1 kg^-1)
  real(wp), parameter :: cp  = 7.0_wp*r_d/2.0_wp

  real(wp), parameter :: r_v  = 461.6_wp
  real(wp), parameter :: cv   = cp - r_d
  real(wp), parameter :: cpv  = 4.0_wp*r_v
  real(wp), parameter :: cvv  = cpv - r_v
  real(wp), parameter :: cvpm = -cv/cp

  real(wp), parameter :: cliq = 4190.0_wp
  real(wp), parameter :: cice = 2106.0_wp
  real(wp), parameter :: psat = 610.78_wp

  real(wp), parameter :: rcv  = r_d/cv
  real(wp), parameter :: rcp  = r_d/cp
  real(wp), parameter :: rovg = r_d/g
  real(wp), parameter :: c2   = cp*rcv

  real(wp), parameter :: mwdry = 28.966_wp    ! g mol^-1

  real(wp), parameter :: p1000mb = 100000.0_wp
  real(wp), parameter :: t0      = 300.0_wp
  real(wp), parameter :: p0      = p1000mb

  real(wp), parameter :: cpovcv = cp/(cp-r_d)
  real(wp), parameter :: cvovcp = 1.0_wp/cpovcv
  real(wp), parameter :: rvovrd = r_v/r_d

  real(wp), parameter :: reradius = 1.0_wp/6370.0e03_wp

  real(wp), parameter :: asselin = 0.025_wp
  real(wp), parameter :: cb      = 25.0_wp

  ! Latent heats
  real(wp), parameter :: XLV0 = 3.15e6_wp
  real(wp), parameter :: XLV1 = 2370.0_wp
  real(wp), parameter :: XLS0 = 2.905e6_wp
  real(wp), parameter :: XLS1 = 259.532_wp

  real(wp), parameter :: XLS = 2.85e6_wp
  real(wp), parameter :: XLV = 2.5e6_wp
  real(wp), parameter :: XLF = 3.50e5_wp

  ! Densities
  real(wp), parameter :: rhowater = 1000.0_wp
  real(wp), parameter :: rhosnow  = 100.0_wp
  real(wp), parameter :: rhoair0  = 1.28_wp

  ! Effective radii
  real(wp), parameter :: RE_QC_BG  = 2.49e-6_wp
  real(wp), parameter :: RE_QI_BG  = 4.99e-6_wp
  real(wp), parameter :: RE_QS_BG  = 9.99e-6_wp
  real(wp), parameter :: RE_QC_MAX = 50.0e-6_wp
  real(wp), parameter :: RE_QI_MAX = 125.0e-6_wp
  real(wp), parameter :: RE_QS_MAX = 999.0e-6_wp

  ! Math / geometry
  real(wp), parameter :: piconst = 3.14159265358979323846_wp
  real(wp), parameter :: DEGRAD  = piconst/180.0_wp
  real(wp), parameter :: DPD     = 360.0_wp/365.0_wp

  ! Thermodynamics
  real(wp), parameter :: SVP1  = 0.6112_wp
  real(wp), parameter :: SVP2  = 17.67_wp
  real(wp), parameter :: SVP3  = 29.65_wp
  real(wp), parameter :: SVPT0 = 273.15_wp

  real(wp), parameter :: EP_1 = r_v/r_d - 1.0_wp
  real(wp), parameter :: EP_2 = r_d/r_v

  real(wp), parameter :: KARMAN = 0.4_wp
  real(wp), parameter :: EOMEG  = 7.2921e-5_wp
  real(wp), parameter :: STBOLT = 5.67051e-8_wp

  ! Turbulence / damping
  real(wp), parameter :: prandtl = 1.0_wp/3.0_wp
  real(wp), parameter :: w_alpha = 0.3_wp
  real(wp), parameter :: w_beta  = 1.0_wp

  ! Misc scheme constants
  real(wp), parameter :: pq0    = 379.90516_wp
  real(wp), parameter :: epsq2  = 0.2_wp
  real(wp), parameter :: a2     = 17.2693882_wp
  real(wp), parameter :: a3     = 273.16_wp
  real(wp), parameter :: a4     = 35.86_wp
  real(wp), parameter :: epsq   = 1.0e-12_wp
  real(wp), parameter :: p608   = rvovrd - 1.0_wp
  real(wp), parameter :: climit = 1.0e-20_wp

  real(wp), parameter :: cm1 = 2937.4_wp
  real(wp), parameter :: cm2 = 4.9283_wp
  real(wp), parameter :: cm3 = 23.5518_wp

  real(wp), parameter :: defc = 0.0_wp
  real(wp), parameter :: defm = 99999.0_wp

  real(wp), parameter :: epsfc = 1.0_wp/1.05_wp
  real(wp), parameter :: epswet = 0.0_wp

  real(wp), parameter :: fcdif = 1.0_wp/3.0_wp
  real(wp), parameter :: fcm   = 3.0e-5_wp

  real(wp), parameter :: gma   = -r_d*(1.0_wp-rcp)*0.5_wp

  real(wp), parameter :: p400  = 40000.0_wp
  real(wp), parameter :: phitp = 15000.0_wp

  real(wp), parameter :: pi1 = 3.1415926_wp
  real(wp), parameter :: pi2 = 2.0_wp*pi1

  real(wp), parameter :: plbtm = 105000.0_wp
  real(wp), parameter :: plomd = 64200.0_wp
  real(wp), parameter :: pmdhi = 35000.0_wp

  real(wp), parameter :: q2ini = 0.50_wp
  real(wp), parameter :: rfcp  = 0.25_wp/cp

  real(wp), parameter :: rhcrit_land = 0.75_wp
  real(wp), parameter :: rhcrit_sea  = 0.80_wp

  real(wp), parameter :: rlag  = 14.8125_wp
  real(wp), parameter :: rlx   = 0.90_wp
  real(wp), parameter :: scq2  = 50.0_wp
  real(wp), parameter :: slopht = 0.001_wp
  real(wp), parameter :: tlc   = 2.0_wp*0.703972477_wp

  real(wp), parameter :: wa   = 0.15_wp
  real(wp), parameter :: wght = 0.35_wp
  real(wp), parameter :: wpc  = 0.075_wp

  real(wp), parameter :: z0land = 0.10_wp
  real(wp), parameter :: z0max  = 0.008_wp
  real(wp), parameter :: z0sea  = 0.001_wp

  !-----------------------------------------
  ! Orbital / Earth parameters
  !-----------------------------------------
  integer, parameter :: PLANET_YEAR = 365
  real(wp), parameter :: OBLIQUITY  = 23.5_wp
  real(wp), parameter :: ECCENTRICITY = 0.014_wp
  real(wp), parameter :: SEMIMAJORAXIS = 1.0_wp
  real(wp), parameter :: zero_date = 0.0_wp
  real(wp), parameter :: EQUINOX_FRACTION = 0.0_wp

#if (EM_CORE == 1)
  integer, parameter :: ZONE_SOLVE_EM = 1
  integer, parameter :: ZONE_SFS      = 2
#endif

CONTAINS

  subroutine init_module_model_constants
    implicit none
  end subroutine init_module_model_constants

END MODULE module_model_constants
