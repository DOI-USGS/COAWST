module module_mp_tempo_cfgs
  !! tempo configs

  implicit none
  private

  public :: ty_tempo_cfgs, ty_tempo_table_cfgs

!> \section arg_table_ty_tempo_cfgs Argument Table
!! \htmlinclude ty_tempo_cfgs.html
!!
  
  ! tempo configuration flags for init
  type :: ty_tempo_cfgs
    logical :: aerosolaware_flag = .true. !! flag to run aerosol-aware microphysics
    logical :: hailaware_flag = .true. !! flag to run hail-aware microphysics
    logical :: ml_for_bl_nc_flag = .false. !! flag to run machine-learning prediction for subgrid cloud number concentration
    logical :: ml_for_nc_flag = .false. !! flag to run machine-learning prediction for tempo cloud number concentration
    logical :: semi_sedi_flag = .false. !! flag for semi-lagrangian sedimentation
    logical :: refl10cm_from_melting_flag = .false. !! flag to calculate reflectivity for melting snow and graupel
    logical :: turn_off_micro_flag = .false. !! flag to turn off all microphysical processes
    logical :: cloud_condensation_flag = .true. !! flag to control cloud condensation
    logical :: verbose = .false. !! flag to turn on verbose print statements
    ! flags to turn on/off diagnostic output
    logical :: refl10cm_flag = .true. !! flag to output 10cm reflectivity
    logical :: re_cloud_flag = .true. !! flag to output cloud effective radius
    logical :: re_ice_flag = .true. !! flag to output ice effective radius
    logical :: re_snow_flag = .true. !! flag to output snow effective radius
    logical :: max_hail_diameter_flag = .true. !! flag to output maximum hail diameter
    logical :: rain_med_vol_diam_flag = .false. !! flag to output median volume diameter for rain
    logical :: graupel_med_vol_diam_flag = .false. !! flag to output median volume diameter for graupel
    logical :: cloud_number_mixing_ratio_flag = .false. !! flag to output cloud number mixing ratio
  end type

  ! tempo lookup table filenames
  type :: ty_tempo_table_cfgs
    character(len=100) :: ccn_table_name = 'ccn_activate.bin' !! ccn table name
    character(len=100) :: qrqg_table_name = 'qr_acr_qg_data_tempo_v3' !! rain-graupel collection table name
    character(len=100) :: qrqs_table_name = 'qr_acr_qs_data_tempo_v3' !! rain-snow collection table name
    character(len=100) :: freezewater_table_name = 'freeze_water_data_tempo_v3' !! freeze water collection table name
  end type

end module module_mp_tempo_cfgs
