module module_mp_tempo_driver
  !! tempo driver that converts 3d model input to 1d arrays used by the main code
  !! also allocates and fills diagnostic arrays
  use module_mp_tempo_cfgs, only : ty_tempo_cfgs, ty_tempo_table_cfgs
  use module_mp_tempo_params, only : wp, sp, dp
  use module_mp_tempo_main, only : tempo_main, ty_tempo_main_diags
  use module_mp_tempo_utils, only : compute_efrw, compute_efsw, compute_drop_evap, qi_aut_qs
  use module_mp_tempo_ml, only : ty_tempo_ml_data, nc_ml_nodes, nc_ml_input, nc_ml_output, &
    nc_ml_trans_mean, nc_ml_trans_var, nc_ml_w00, nc_ml_w01, nc_ml_b00, nc_ml_b01, save_or_read_ml_data

  implicit none
  private

  public :: tempo_init, tempo_run, ty_tempo_driver_diags, tempo_aerosol_surface_emissions
  
  type(ty_tempo_table_cfgs) :: tempo_table_cfgs

  type :: ty_tempo_driver_diags
    real(wp), dimension(:,:), allocatable :: rain_precip
    real(wp), dimension(:,:), allocatable :: ice_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: snow_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: graupel_liquid_equiv_precip
    real(wp), dimension(:,:), allocatable :: frozen_fraction
    real(wp), dimension(:,:), allocatable :: frz_rain_precip
    real(wp), dimension(:,:), allocatable :: max_hail_diameter_sfc
    real(wp), dimension(:,:), allocatable :: max_hail_diameter_column
    real(wp), dimension(:,:,:), allocatable :: refl10cm
    real(wp), dimension(:,:,:), allocatable :: re_cloud
    real(wp), dimension(:,:,:), allocatable :: re_ice
    real(wp), dimension(:,:,:), allocatable :: re_snow
    real(wp), dimension(:,:,:), allocatable :: rain_med_vol_diam
    real(wp), dimension(:,:,:), allocatable :: graupel_med_vol_diam
    real(wp), dimension(:,:,:), allocatable :: cloud_number_mixing_ratio
  end type

  contains

!> initialize tempo microphysics
!! \section arg_table_tempo_init Argument Table
!! \htmlinclude tempo_init.html
!!
  subroutine tempo_init(aerosolaware_flag, hailaware_flag, semi_sedi_flag, cloud_condensation_flag, &
    refl10cm_from_melting_flag, ml_for_bl_nc_flag, ml_for_nc_flag, force_init_flag, tempo_cfgs)
    !! initialize tempo microphysics
    use module_mp_tempo_params, only : get_version, tempo_version, t_efrw, &
      initialize_graupel_vars, initialize_parameters, initialize_bins_for_tables, &
      initialize_array_efrw, initialize_array_efsw, initialize_arrays_drop_evap, &
      initialize_arrays_ccn, initialize_arrays_qi_aut_qs, &
      initialize_arrays_qr_acr_qs, initialize_arrays_qr_acr_qg, initialize_arrays_freezewater, &
      initialize_bins_for_hail_size, initialize_bins_for_radar

    logical, intent(in), optional :: aerosolaware_flag, hailaware_flag, refl10cm_from_melting_flag, &
      ml_for_bl_nc_flag, ml_for_nc_flag, force_init_flag, semi_sedi_flag, cloud_condensation_flag
    type(ty_tempo_cfgs), intent(inout) :: tempo_cfgs

    character(len=100) :: table_filename
    integer :: table_size
    logical :: initialize_mp_vars, force_init

    ! get tempo version from readme file
    call get_version(tempo_version) 

    ! check an allocatable array (t_efrw) to see if initialization can be skipped
    ! but allow for force initialization useful for testing
    force_init = .false.
    if (present(force_init_flag)) force_init = force_init_flag

    initialize_mp_vars = .true.
    if (allocated(t_efrw)) initialize_mp_vars = .false.
    if (force_init) initialize_mp_vars = .true.

    if (initialize_mp_vars) then
      if (present(aerosolaware_flag)) tempo_cfgs%aerosolaware_flag = aerosolaware_flag
      if (present(hailaware_flag)) tempo_cfgs%hailaware_flag = hailaware_flag
      if (present(ml_for_bl_nc_flag)) tempo_cfgs%ml_for_bl_nc_flag = ml_for_bl_nc_flag
      if (present(ml_for_nc_flag)) tempo_cfgs%ml_for_nc_flag = ml_for_nc_flag
      if (present(semi_sedi_flag)) tempo_cfgs%semi_sedi_flag = semi_sedi_flag
      if (present(cloud_condensation_flag)) tempo_cfgs%cloud_condensation_flag = cloud_condensation_flag
      if (present(refl10cm_from_melting_flag)) tempo_cfgs%refl10cm_from_melting_flag = refl10cm_from_melting_flag

      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- TEMPO microphysics configuration options: '
        write(*,'(A,L)') 'tempo_init() --- aerosol aware = ', tempo_cfgs%aerosolaware_flag
        write(*,'(A,L)') 'tempo_init() --- hail aware = ', tempo_cfgs%hailaware_flag
        write(*,'(A,L)') 'tempo_init() --- ML for subgrid cloud number = ', tempo_cfgs%ml_for_bl_nc_flag
        write(*,'(A,L)') 'tempo_init() --- ML for cloud number = ', tempo_cfgs%ml_for_nc_flag
        write(*,'(A,L)') 'tempo_init() --- reflectivity from melting snow/graupel = ', tempo_cfgs%refl10cm_from_melting_flag
        write(*,'(A,L)') 'tempo_init() --- semi-lagrangian sedimentation = ', tempo_cfgs%semi_sedi_flag
      endif 

      ! set graupel variables from hail_aware_flag
      call initialize_graupel_vars(tempo_cfgs%hailaware_flag) 
      if (tempo_cfgs%verbose) then
        write(*,'(A,L)') 'tempo_init() --- initialized graupel variables using hail aware = ', tempo_cfgs%hailaware_flag
      endif 

      ! set parameters that can depend on the host model
      call initialize_parameters() 
      if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized parameters'
      
      ! creates log-spaced bins of hydrometers for tables
      call initialize_bins_for_tables() 
      if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized bins for lookup tables'

      ! collision efficiencies between rain/snow and cloud water.
      call initialize_array_efrw()
      call compute_efrw()
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized collision efficiency data for rain collecting cloud water'
      endif 
      call initialize_array_efsw()
      call compute_efsw()
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized collision efficiency data for snow collecting cloud water'
      endif 

      ! drop evaporation
      call initialize_arrays_drop_evap()
      call compute_drop_evap()
      if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized drop evaporation data'

      ! cloud ice to snow and depositional growth
      call initialize_arrays_qi_aut_qs()
      call qi_aut_qs()

      ! CCN activation table
      table_filename = tempo_table_cfgs%ccn_table_name
      call initialize_arrays_ccn(table_size)
      call read_table_ccn(trim(table_filename), table_size)
      if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized data for ccn lookup table'

      ! freeze water collection lookup table
      table_filename = tempo_table_cfgs%freezewater_table_name
      call initialize_arrays_freezewater(table_size)
      call read_table_freezewater(trim(table_filename), table_size)
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized data for frozen cloud water and rain lookup table'
      endif 

      ! rain-snow collection lookup table
      table_filename = tempo_table_cfgs%qrqs_table_name
      call initialize_arrays_qr_acr_qs(table_size)
      call read_table_qr_acr_qs(trim(table_filename), table_size)
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized data for rain-snow collection lookup table'
      endif 

      ! rain-graupel collection lookup table
      table_filename = tempo_table_cfgs%qrqg_table_name
      call initialize_arrays_qr_acr_qg(table_size)
      call read_table_qr_acr_qg(trim(table_filename), table_size)
      if (tempo_cfgs%verbose) then
        write(*,'(A)') 'tempo_init() --- initialized data for rain-graupel collection lookup table'
      endif 

      ! bins used for optional refl10cm calculation with melting
      if (tempo_cfgs%refl10cm_from_melting_flag) then
        call initialize_bins_for_radar()
        if (tempo_cfgs%verbose) then
          write(*,'(A,L)') 'tempo_init() ---  flag to calcuate reflectivity with contributions from melting snow and graupel = ', &
            tempo_cfgs%refl10cm_from_melting_flag
          write(*,'(A)') 'tempo_init() --- initialized bins for reflectivity calcuation with meting snow and graupel'
        endif 
      endif

      ! bins used for optional hail size calculation
      if (tempo_cfgs%max_hail_diameter_flag) then
        call initialize_bins_for_hail_size()
        if (tempo_cfgs%verbose) then
          write(*,'(A,L)') 'tempo_init() ---  flag to calculate max hail diameter = ', &
            tempo_cfgs%max_hail_diameter_flag
          write(*,'(A)') 'tempo_init() --- initialized bins for hail size calculation'
        endif
      endif

      ! data for machine learning
      if(tempo_cfgs%ml_for_bl_nc_flag .or. tempo_cfgs%ml_for_nc_flag) then
        call init_ml_data()
        if (tempo_cfgs%verbose) write(*,'(A)') 'tempo_init() --- initialized data for cloud number machine learning'
      endif 
    endif
  end subroutine tempo_init

!> \section arg_table_tempo_run Argument Table
!! \htmlinclude tempo_run.html
!!
  subroutine tempo_run(tempo_cfgs, dt, itimestep, &
    t, th, pii, p, w, dz, &
    qv, qc, qr, qi, qs, qg, ni, nr, &
    nc, nwfa, nifa, ng, qb, &
    qc_bl, qcfrac_bl, &
    qcfrac, qifrac, &
    thten_bl, qvten_bl, qcten_bl, qiten_bl, &
    thten_lwrad, thten_swrad, &
    ids, ide, jds, jde, kds, kde, &
    ims, ime, jms, jme, kms, kme, &
    its, ite, jts, jte, kts, kte, tempo_diags)

    type(ty_tempo_cfgs), intent(in) :: tempo_cfgs
    real(wp), intent(in) :: dt !! timestep \([s]]\)
    integer, intent(in) :: itimestep !! integer timestep = integration time / dt
    integer, intent(in) :: ids, ide, jds, jde, kds, kde !! domain locations
    integer, intent(in) :: ims, ime, jms, jme, kms, kme !! memory locations
    integer, intent(in) :: its, ite, jts, jte, kts, kte !! tile locations

    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: t !! temperature \([K]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: th !! theta \([K]\)

    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in) :: p !! pressure \([Pa]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in) :: w !! vertical velocity \([m\; s^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in) :: dz !! vertical grid spacing \([m]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: pii !! exner function

    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qv !! 3D water vapor mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qc !! 3D cloud water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qr !! 3D rain water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qi !! 3D cloud ice mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qs !! 3D snow mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qg !! 3D graupel mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: ni !! 3D cloud ice number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: nr !! 3D rain water number mixing ratio \([kg^{-1}]\)

    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: nc !! 3D cloud water number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: nwfa !! 3D water-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: nifa !! 3D ice-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: qb !! 3D graupel volume mixing ratio \([m^{-3}\; kg^{-1}]\) (hail-aware)
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: ng !! 3D graupel number mixing ratio \([kg^{-1}]\) (hail-aware)

    ! additional optional arguments
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: qcfrac
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: qifrac
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: qc_bl
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: qcfrac_bl
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: thten_bl
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: qvten_bl
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: qcten_bl
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: qiten_bl
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: thten_lwrad
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(in), optional :: thten_swrad

    real(wp), dimension(kts:kte) :: t1d !! 1D temperature \([K]\)
    real(wp), dimension(kts:kte) :: p1d !! 1D pressure \([Pa]\)
    real(wp), dimension(kts:kte) :: qv1d !! 1D water vapor mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qc1d !! 1D cloud water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qr1d !! 1D rain water mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qi1d !! 1D cloud ice mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qs1d !! 1D snow mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: qg1d !! 1D graupel mass mixing ratio \([kg\; kg^{-1}]\)
    real(wp), dimension(kts:kte) :: ni1d !! 1D cloud ice number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(kts:kte) :: nr1d !! 1D rain water number mixing ratio \([kg^{-1}]\)
    real(wp), dimension(kts:kte) :: w1d !! 1D vertical velocity \(m\; s^{-1}]\)
    real(wp), dimension(kts:kte) :: dz1d !! 1D vertical grid spacing \([m]\)

    real(wp), dimension(:), allocatable :: nc1d !! 1D cloud water number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(:), allocatable :: nwfa1d !! 1D water-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(:), allocatable :: nifa1d !! 1D ice-friendly aerosol number mixing ratio \([kg^{-1}]\) (aerosol-aware)
    real(wp), dimension(:), allocatable :: qb1d !! 1D graupel volume mixing ratio \([m^{-3}\; kg^{-1}]\) (hail-aware)
    real(wp), dimension(:), allocatable :: ng1d !! 1D graupel number mixing ratio \([kg^{-1}]\) (hail-aware)

    ! additional optional 1d arrays
    real(wp), dimension(:), allocatable :: qcfrac1d
    real(wp), dimension(:), allocatable :: qifrac1d
    real(wp), dimension(:), allocatable :: qc_bl1d
    real(wp), dimension(:), allocatable :: qcfrac_bl1d
    real(wp), dimension(:), allocatable :: thten_bl1d
    real(wp), dimension(:), allocatable :: qvten_bl1d
    real(wp), dimension(:), allocatable :: qcten_bl1d
    real(wp), dimension(:), allocatable :: qiten_bl1d
    real(wp), dimension(:), allocatable :: thten_lwrad1d
    real(wp), dimension(:), allocatable :: thten_swrad1d

    integer :: i, j, k, nz
    logical :: use_temperature 

    type(ty_tempo_main_diags) :: tempo_main_diags
    type(ty_tempo_driver_diags), intent(inout) :: tempo_diags

    nz = kte - kts + 1
    ! allocate 1d arrays if 3d arrays are present
    if (present(nwfa)) allocate(nwfa1d(nz), source=0._wp)
    if (present(nifa)) allocate(nifa1d(nz), source=0._wp)
    if (present(nc)) allocate(nc1d(nz), source=0._wp)
    if (present(ng)) allocate(ng1d(nz), source=0._wp)
    if (present(qb)) allocate(qb1d(nz), source=0._wp)

    ! additional optional 1d arrays
    if (present(qcfrac)) allocate(qcfrac1d(nz), source=0._wp)
    if (present(qifrac)) allocate(qifrac1d(nz), source=0._wp)
    if (present(qc_bl)) allocate(qc_bl1d(nz), source=0._wp)
    if (present(qcfrac_bl)) allocate(qcfrac_bl1d(nz), source=0._wp)
    if (present(thten_bl)) allocate(thten_bl1d(nz), source=0._wp)
    if (present(qvten_bl)) allocate(qvten_bl1d(nz), source=0._wp)
    if (present(qcten_bl)) allocate(qcten_bl1d(nz), source=0._wp)  
    if (present(qiten_bl)) allocate(qiten_bl1d(nz), source=0._wp)
    if (present(thten_lwrad)) allocate(thten_lwrad1d(nz), source=0._wp)
    if (present(thten_swrad)) allocate(thten_swrad1d(nz), source=0._wp) 

    ! allocate diagnostics
    ! 3d diagnostics have configuration flags
    if (tempo_cfgs%cloud_number_mixing_ratio_flag) then
      if (.not. allocated(tempo_diags%cloud_number_mixing_ratio)) then
        allocate(tempo_diags%cloud_number_mixing_ratio(its:ite, kts:kte, jts:jte), source=0._wp)
      else
        tempo_diags%cloud_number_mixing_ratio = 0._wp
      endif
    endif

    if (tempo_cfgs%rain_med_vol_diam_flag) then
      if (.not. allocated(tempo_diags%rain_med_vol_diam)) then
        allocate(tempo_diags%rain_med_vol_diam(its:ite, kts:kte, jts:jte), source=0._wp)
      else
        tempo_diags%rain_med_vol_diam = 0._wp
      endif
    endif

    if (tempo_cfgs%graupel_med_vol_diam_flag) then
      if (.not. allocated(tempo_diags%graupel_med_vol_diam)) then
        allocate(tempo_diags%graupel_med_vol_diam(its:ite, kts:kte, jts:jte), source=0._wp)
      else
        tempo_diags%graupel_med_vol_diam = 0._wp
      endif
    endif

    if (tempo_cfgs%refl10cm_flag) then
      if (.not. allocated(tempo_diags%refl10cm)) then
        allocate(tempo_diags%refl10cm(its:ite, kts:kte, jts:jte), source=-35._wp)
      else
        tempo_diags%refl10cm = -35._wp
      endif
    endif

    if (tempo_cfgs%re_cloud_flag) then
      if (.not. allocated(tempo_diags%re_cloud)) then
        allocate(tempo_diags%re_cloud(its:ite, kts:kte, jts:jte), source=0._wp)
      else
        tempo_diags%re_cloud = 0._wp
      endif
    endif

    if (tempo_cfgs%re_ice_flag) then
      if (.not. allocated(tempo_diags%re_ice)) then
        allocate(tempo_diags%re_ice(its:ite, kts:kte, jts:jte), source=0._wp)
      else
        tempo_diags%re_ice = 0._wp
      endif
    endif

    if (tempo_cfgs%re_snow_flag) then
      if (.not. allocated(tempo_diags%re_snow)) then
        allocate(tempo_diags%re_snow(its:ite, kts:kte, jts:jte), source=0._wp)
      else
        tempo_diags%re_snow = 0._wp
      endif
    endif

    ! 2d diagnostics
    if (tempo_cfgs%max_hail_diameter_flag) then
      if (.not. allocated(tempo_diags%max_hail_diameter_sfc)) then
        allocate(tempo_diags%max_hail_diameter_sfc(its:ite, jts:jte), source=0._wp)
        allocate(tempo_diags%max_hail_diameter_column(its:ite, jts:jte), source=0._wp)
      else
        tempo_diags%max_hail_diameter_sfc = 0._wp
        tempo_diags%max_hail_diameter_column = 0._wp
      endif
    endif

    ! precipitation
    if (.not. allocated(tempo_diags%rain_precip)) then
      allocate(tempo_diags%rain_precip(its:ite, jts:jte), source=0._wp)
      allocate(tempo_diags%ice_liquid_equiv_precip(its:ite, jts:jte), source=0._wp)
      allocate(tempo_diags%snow_liquid_equiv_precip(its:ite, jts:jte), source=0._wp)
      allocate(tempo_diags%graupel_liquid_equiv_precip(its:ite, jts:jte), source=0._wp)
      allocate(tempo_diags%frozen_fraction(its:ite, jts:jte), source=0._wp)
      allocate(tempo_diags%frz_rain_precip(its:ite, jts:jte), source=0._wp)
    else
      tempo_diags%rain_precip = 0._wp
      tempo_diags%ice_liquid_equiv_precip = 0._wp
      tempo_diags%snow_liquid_equiv_precip = 0._wp
      tempo_diags%graupel_liquid_equiv_precip = 0._wp
      tempo_diags%frozen_fraction = 0._wp
      tempo_diags%frz_rain_precip = 0._wp
    endif

    ! temperature or theta and exner
    if (present(t)) then
      use_temperature = .true.
    elseif (present(th) .and. present(pii)) then
      use_temperature = .false.
    else  
      error stop "tempo_run() --- requires either temperature or theta and Exner function"
    endif 

    ! tempo driver code
    do j = jts, jte
      do i = its, ite
        do k = kts, kte
          if (use_temperature) then
            t1d(k) = t(i,k,j)
          else  
            t1d(k) = th(i,k,j) * pii(i,k,j)
          endif 
          p1d(k) = p(i,k,j)
          w1d(k) = w(i,k,j)
          dz1d(k) = dz(i,k,j)
          qv1d(k) = qv(i,k,j)
          qc1d(k) = qc(i,k,j)
          qi1d(k) = qi(i,k,j)
          qr1d(k) = qr(i,k,j)
          qs1d(k) = qs(i,k,j)
          qg1d(k) = qg(i,k,j)
          ni1d(k) = ni(i,k,j)
          nr1d(k) = nr(i,k,j)

          ! nwfa, nifa, and nc are optional aerosol-aware variables
          if (present(nwfa)) nwfa1d(k) = nwfa(i,k,j)
          if (present(nifa)) nifa1d(k) = nifa(i,k,j)
          if (present(nc)) nc1d(k) = nc(i,k,j)

          ! ng and qb are optional hail-aware variables
          if ((present(ng)) .and. (present(qb))) then
            ng1d(k) = ng(i,k,j)
            qb1d(k) = qb(i,k,j)
          endif 

          ! machine learning for pbl clouds
          if (present(qc_bl) .and. present(qcfrac_bl)) then
            qc_bl1d(k) = qc_bl(i,k,j)
            qcfrac_bl1d(k) = qcfrac_bl(i,k,j)
          endif 
        enddo

        ! main call to the 1d tempo microphysics
        call tempo_main(tempo_cfgs=tempo_cfgs, &
          qv1d=qv1d, qc1d=qc1d, qi1d=qi1d, qr1d=qr1d, qs1d=qs1d, qg1d=qg1d, qb1d=qb1d, &
          ni1d=ni1d, nr1d=nr1d, nc1d=nc1d, ng1d=ng1d, nwfa1d=nwfa1d, nifa1d=nifa1d, t1d=t1d, p1d=p1d, &
          w1d=w1d, dz1d=dz1d, &
          qcfrac1d=qcfrac1d, qifrac1d=qifrac1d, qc_bl1d=qc_bl1d, qcfrac_bl1d=qcfrac_bl1d, &
          thten_bl1d=thten_bl1d, qvten_bl1d=qvten_bl1d, qcten_bl1d=qcten_bl1d, qiten_bl1d=qiten_bl1d, &
          thten_lwrad1d=thten_lwrad1d, thten_swrad1d=thten_swrad1d, &
          kts=kts, kte=kte, dt=dt, ii=i, jj=j, tempo_main_diags=tempo_main_diags)
          
        ! precipitation
        tempo_diags%rain_precip(i,j) = tempo_main_diags%rain_precip
        tempo_diags%ice_liquid_equiv_precip(i,j) = tempo_main_diags%ice_liquid_equiv_precip
        tempo_diags%snow_liquid_equiv_precip(i,j) = tempo_main_diags%snow_liquid_equiv_precip
        tempo_diags%graupel_liquid_equiv_precip(i,j) = tempo_main_diags%graupel_liquid_equiv_precip
        tempo_diags%frozen_fraction(i,j) = tempo_main_diags%frozen_fraction
        tempo_diags%frz_rain_precip(i,j) = tempo_main_diags%frz_rain_precip

        ! 3d diagnostics
        if (allocated(tempo_diags%cloud_number_mixing_ratio) .and. allocated(tempo_main_diags%cloud_number_mixing_ratio)) then
          tempo_diags%cloud_number_mixing_ratio(i,:,j) = tempo_main_diags%cloud_number_mixing_ratio
        endif 
        if (allocated(tempo_diags%rain_med_vol_diam) .and. allocated(tempo_main_diags%rain_med_vol_diam)) then
          tempo_diags%rain_med_vol_diam(i,:,j) = tempo_main_diags%rain_med_vol_diam
        endif 
        if (allocated(tempo_diags%graupel_med_vol_diam) .and. allocated(tempo_main_diags%graupel_med_vol_diam)) then
          tempo_diags%graupel_med_vol_diam(i,:,j) = tempo_main_diags%graupel_med_vol_diam
        endif 
        if (allocated(tempo_diags%re_cloud) .and. allocated(tempo_main_diags%re_cloud)) then
          tempo_diags%re_cloud(i,:,j) = tempo_main_diags%re_cloud
        endif 
        if (allocated(tempo_diags%re_ice) .and. allocated(tempo_main_diags%re_ice)) then
          tempo_diags%re_ice(i,:,j) = tempo_main_diags%re_ice
        endif 
        if (allocated(tempo_diags%re_snow) .and. allocated(tempo_main_diags%re_snow)) then
          tempo_diags%re_snow(i,:,j) = tempo_main_diags%re_snow
        endif 
        if (allocated(tempo_diags%refl10cm) .and. allocated(tempo_main_diags%refl10cm)) then
          tempo_diags%refl10cm(i,:,j) = tempo_main_diags%refl10cm
        endif
        if (allocated(tempo_diags%max_hail_diameter_sfc) .and. allocated(tempo_diags%max_hail_diameter_column) .and. allocated(tempo_main_diags%max_hail_diameter)) then
          tempo_diags%max_hail_diameter_sfc(i,j) = tempo_main_diags%max_hail_diameter(1)
          tempo_diags%max_hail_diameter_column(i,j) = maxval(tempo_main_diags%max_hail_diameter)
        endif

        ! return variables to model
        do k = kts, kte
          if (present(nc)) nc(i,k,j) = nc1d(k)
          if (present(nwfa)) nwfa(i,k,j) = nwfa1d(k)
          if (present(nifa)) nifa(i,k,j) = nifa1d(k)
          if ((present(ng)) .and. (present(qb))) then
            ng(i,k,j) = ng1d(k)
            qb(i,k,j) = qb1d(k)
          endif 
          qv(i,k,j) = qv1d(k)
          qc(i,k,j) = qc1d(k)
          qi(i,k,j) = qi1d(k)
          qr(i,k,j) = qr1d(k)
          qs(i,k,j) = qs1d(k)
          qg(i,k,j) = qg1d(k)
          ni(i,k,j) = ni1d(k)
          nr(i,k,j) = nr1d(k)
          
          ! tempo main returns temperature (t1d), so convert to theta if needed
          if (use_temperature) then
            t(i,k,j) = t1d(k)
          else  
            th(i,k,j) = t1d(k) / pii(i,k,j)
          endif 
        enddo
      enddo
    enddo
  end subroutine tempo_run


  subroutine tempo_aerosol_surface_emissions(dt, nwfa, nwfa2d, ims, ime, jms, jme, kms, kme, kts)
    !! adds aerosol surface emissions to the 3D field
    real(wp), intent(in) :: dt
    real(wp), dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: nwfa 
    real(wp), dimension(ims:ime, jms:jme), intent(in) :: nwfa2d
    integer, intent(in) :: ims, ime, jms, jme, kms, kme, kts
    integer :: i, j

    do j = jms, jme
      do i = ims, ime
        nwfa(i,kts,j) = nwfa(i,kts,j) + nwfa2d(i,j) * dt
      enddo
    enddo
  end subroutine tempo_aerosol_surface_emissions


  subroutine read_table_freezewater(filename, table_size)
    !! read lookup table for frozen cloud and rain water
    use module_mp_tempo_params, only : tpi_qrfz, tni_qrfz, &
      tpg_qrfz, tnr_qrfz, tpi_qcfz, tni_qcfz

    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size
    
    integer :: mp_unit, istat

    mp_unit = 11
    call check_before_table_read(filename, table_size)
    open(unit=mp_unit, file=filename, form='unformatted', status='old', access='stream', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tpi_qrfz
    read(mp_unit) tni_qrfz
    read(mp_unit) tpg_qrfz
    read(mp_unit) tnr_qrfz
    read(mp_unit) tpi_qcfz
    read(mp_unit) tni_qcfz
    close(unit=mp_unit)
  end subroutine read_table_freezewater


  subroutine read_table_qr_acr_qs(filename, table_size)
    !! read lookup table for rain-snow collection
    use module_mp_tempo_params, only : tcs_racs1, tmr_racs1, &
      tcs_racs2, tmr_racs2, tcr_sacr1, tms_sacr1, tcr_sacr2, &
      tms_sacr2, tnr_racs1, tnr_racs2, tnr_sacr1, tnr_sacr2

    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size
    
    integer :: mp_unit, istat

    mp_unit = 11
    call check_before_table_read(filename, table_size)
    open(unit=mp_unit, file=filename, form='unformatted', status='old', access='stream', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tcs_racs1
    read(mp_unit) tmr_racs1
    read(mp_unit) tcs_racs2
    read(mp_unit) tmr_racs2
    read(mp_unit) tcr_sacr1
    read(mp_unit) tms_sacr1
    read(mp_unit) tcr_sacr2
    read(mp_unit) tms_sacr2
    read(mp_unit) tnr_racs1
    read(mp_unit) tnr_racs2
    read(mp_unit) tnr_sacr1
    read(mp_unit) tnr_sacr2
    close(unit=mp_unit)
  end subroutine read_table_qr_acr_qs


  subroutine read_table_qr_acr_qg(filename, table_size)
    !! read lookup table for rain-graupel collection
    use module_mp_tempo_params, only : tcg_racg, tmr_racg, &
      tcr_gacr, tnr_racg, tnr_gacr

    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size
    
    integer :: mp_unit, istat

    mp_unit = 11
    call check_before_table_read(filename, table_size)
    open(unit=mp_unit, file=filename, form='unformatted', status='old', access='stream', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tcg_racg
    read(mp_unit) tmr_racg
    read(mp_unit) tcr_gacr
    read(mp_unit) tnr_racg
    read(mp_unit) tnr_gacr
    close(unit=mp_unit)
  end subroutine read_table_qr_acr_qg


  subroutine read_table_ccn(filename, table_size)
    !! read static file containing CCN activation of aerosols;
    !! the data were created from a parcel model by Feingold and Heymsfield (1992)
    !! https://doi.org/10.1175/1520-0469(1992)049<2325:POCGOD>2.0.CO;2
    !! with further changes by Eidhammer and Kreidenweis
    use module_mp_tempo_params, only : tnccn_act
  
    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size
    
    integer :: mp_unit, istat

    call check_before_table_read(filename=filename, table_size=table_size)

    mp_unit = 11
    open(unit=mp_unit, file=filename, form='unformatted', status='old', &
      action='read', iostat=istat, convert='big_endian')
    read(mp_unit) tnccn_act
    close(unit=mp_unit)
  end subroutine read_table_ccn


  subroutine check_before_table_read(filename, table_size)
    !! checks that lookup tables exist and are the correct size
    !! before attempting to read them

    character(len=*), intent(in) :: filename
    integer, intent(in) :: table_size

    logical :: fileexists
    integer :: filesize
    character(len=100) :: int_to_str1, int_to_str2

    inquire(file=filename, size=filesize, exist=fileexists)
    if (.not. fileexists) then
      write(*,'(3A)') 'tempo_init() --- *** FATAL *** file "', filename, &
        '" was not found in this directory.'
      write(*,'(A)') ''
      write(*,'(A)') 'How to fix issues with tables (datasets stored in files):'
      write(*,'(3A)') '   (1) The table ', trim(tempo_table_cfgs%ccn_table_name), &
        ' is located in the TEMPO/tables/ directory. Copy this file to the directory where the model executable is located.'
      write(*,'(8A)') '   (2) Three tables:', trim(tempo_table_cfgs%qrqs_table_name), ', ', trim(tempo_table_cfgs%qrqg_table_name), ', and ', &
        trim(tempo_table_cfgs%freezewater_table_name), &
        ' can be build by compiling and running the executable "build_tables" in the main TEMPO directory. ', &
        'Then copy the file to the directory where the model executable is located.'
      write(*,'(A)') '   (3) Ask the developers for tables. They are willing to share.' 
      write(*,'(A)') ''      
      error stop '--- file "' // filename // '" needed for TEMPO microphysics was not found.'
    endif

    if (filesize /= table_size) then
      write(int_to_str1, '(I0)') filesize
      write(int_to_str2, '(I0)') table_size
      write(*,'(7A)') 'tempo_init() --- *** FATAL *** file "', filename, '" has a size of ', &
        trim(int_to_str1), ' bytes but the array allocated to hold the data expects a file size of ', &
          trim(int_to_str2), ' bytes.'
      write(*,'(A)') ''
      write(*,'(A)') 'How to fix issues with tables (datasets stored in files):'
      write(*,'(3A)') '   (1) The table ', trim(tempo_table_cfgs%ccn_table_name), &
        ' is located in the TEMPO/tables/ directory. Copy this file to the directory where the model executable is located.'
      write(*,'(8A)') '   (2) Three tables: ', trim(tempo_table_cfgs%qrqs_table_name), ', ', trim(tempo_table_cfgs%qrqg_table_name), ', and ', &
        trim(tempo_table_cfgs%freezewater_table_name), &
        ' can be build by compiling and running the executable "build_tables" in the main TEMPO directory. ', &
        'Then copy the file to the directory where the model executable is located.'
      write(*,'(A)') '   (3) Ask the developers for tables. They are willing to share.' 
      write(*,'(A)') ''
      error stop '--- size of file "' // filename // '" needed for TEMPO microphysics is inconsistent with expected size.'
    endif
  end subroutine check_before_table_read


  subroutine init_ml_data()
    !! initialize machine learning data for tempo microphysics
    type(ty_tempo_ml_data), dimension(1) :: tempo_ml_data

    ! cloud water
    tempo_ml_data(1)%input_size = nc_ml_input
    tempo_ml_data(1)%node_size = nc_ml_nodes
    tempo_ml_data(1)%output_size = nc_ml_output

    if (.not.allocated(tempo_ml_data(1)%transform_mean)) allocate(tempo_ml_data(1)%transform_mean(nc_ml_input))
    if (.not.allocated(tempo_ml_data(1)%transform_var)) allocate(tempo_ml_data(1)%transform_var(nc_ml_input))

    tempo_ml_data(1)%transform_mean = nc_ml_trans_mean
    tempo_ml_data(1)%transform_var = nc_ml_trans_var

    if (.not.allocated(tempo_ml_data(1)%weights00)) allocate(tempo_ml_data(1)%weights00(nc_ml_nodes,nc_ml_input))
    if (.not.allocated(tempo_ml_data(1)%weights01)) allocate(tempo_ml_data(1)%weights01(nc_ml_output,nc_ml_nodes))
    if (.not.allocated(tempo_ml_data(1)%bias00)) allocate(tempo_ml_data(1)%bias00(nc_ml_nodes))
    if (.not.allocated(tempo_ml_data(1)%bias01)) allocate(tempo_ml_data(1)%bias01(nc_ml_output))

    tempo_ml_data(1)%weights00 = reshape(nc_ml_w00, (/nc_ml_nodes, nc_ml_input/))
    tempo_ml_data(1)%weights01 = reshape(nc_ml_w01, (/nc_ml_output, nc_ml_nodes/))
    tempo_ml_data(1)%bias00 = nc_ml_b00
    tempo_ml_data(1)%bias01 = nc_ml_b01

    ! save neural network
    call save_or_read_ml_data(ml_data_in=tempo_ml_data)
  end subroutine init_ml_data

end module module_mp_tempo_driver
