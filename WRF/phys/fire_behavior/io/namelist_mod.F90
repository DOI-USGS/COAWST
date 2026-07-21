  module namelist_mod

#ifdef DM_PARALLEL
    use mpi_f08
#endif
      ! Get access to default options
    use interp_mod, only : HINTERP_BILINEAR, VINTERP_WINDS_FROM_10M_WINDS
    use fuel_mod, only : FUEL_ANDERSON
    use ros_mod, only : ROS_WRFFIRE
    use fmc_mod, only : FMC_WRFFIRE
    use emis_mod, only : EMIS_WRFFIRE
    use stderrout_mod, only : Stop_simulation, Print_message
    use mpi_mod, only : Convert_mpi_comm_to_f08

    implicit none

    private

    public :: namelist_t, FIRE_MAX_IGNITIONS_IN_NAMELIST

    integer, parameter :: FIRE_MAX_IGNITIONS_IN_NAMELIST = 5

    type :: namelist_t
        ! time block
      integer :: start_year = -1              ! start year of the simulation
      integer :: start_month = -1             ! Start month
      integer :: start_day = -1               ! Start day
      integer :: start_hour = -1              ! Start hour
      integer :: start_minute = -1            ! Start minute
      integer :: start_second = -1            ! Start sencond
      integer :: end_year = -1                ! End year of the simulation
      integer :: end_month = -1               ! End month
      integer :: end_day = -1                 ! End day
      integer :: end_hour = -1                ! End hour
      integer :: end_minute = -1              ! End minute
      integer :: end_second = -1              ! End second

      integer :: interval_output = -1         ! Frequency to save the output [s]
      real :: dt = 2.0                        ! Time step of the fire model [s]

      integer :: num_tiles = 1                ! Number of tiles in each patch
      integer :: tile_strategy = 0            ! Strategy for the tile decomposition: 0) ...

        ! fire block
      integer :: fire_print_msg = 0           ! "write fire statistics, 0 no writes, 1+ for more"  ""
      real :: fire_atm_feedback = 1.0         ! "the heat fluxes to the atmosphere are multiplied by this" "1"
      logical :: fire_is_real_perim = .false. ! .false. = point/line ignition, .true. = observed perimeter"

      integer :: fire_upwinding = 9           ! "Numerical method for normal spread: 1=standard, 2=godunov, 3=eno, 4=sethian, 5=2nd-order,
                                              ! 6=WENO3, 7=WENO5, 8=hybrid WENO3/ENO1, 9=hybrid WENO5/ENO1" "1"
      real :: fire_viscosity = 0.4            ! artificial viscosity in level set method farm from near-front region
      real :: fire_viscosity_bg = 0.4         ! artificial viscosity in level set method the near-front region. Should be lower/equal to fire_viscosity
      integer :: fire_viscosity_ngp = 2       ! number of grid points around lfn=0 where fire_viscosity_bg is used. Defines the near front region
      real :: fire_viscosity_band = 0.5       ! number of times the near front region, used to transition from fire_viscosity_bg to fire_viscosity

      logical :: fire_lsm_reinit = .true.     ! "flag to activate reinitialization of level set method"
      integer :: fire_lsm_reinit_iter = 1     ! "number of iterations for the reinitialization PDE"
      integer :: fire_upwinding_reinit = 4    ! "numerical scheme (space) for reinitialization PDE: 1=WENO3, 2=WENO5, 3=hybrid WENO3-ENO1, 4=hybrid WENO5-ENO1"
      integer :: fire_lsm_band_ngp = 4        ! "number of grid points around lfn=0 that WENO5/3 is used (ENO1 elsewhere),
                                              ! for fire_upwinding_reinit=4,5 and fire_upwinding=8,9 options"
      real :: reinit_pseudot_coef = 0.01      ! Coefficient for the pseudo time

      integer :: fast_dist_reinit_opt = 0     ! Fast distance reinitialization method (or eikonal solver): 0) None, 1) FSM
      integer :: fast_dist_reinit_freq = 600  ! Number of time steps to perform a reinit with fast distance reinit method

      real :: fire_wind_height = 6.096        ! "height of uah,vah wind in fire spread formula" "m"
      integer :: wind_vinterp_opt = VINTERP_WINDS_FROM_10M_WINDS ! "mid-flame height wind interpolation option: 0) Interp to specified height, 1) Use WAFs"
      integer :: hinterp_opt = HINTERP_BILINEAR ! "Horizontal interpolation from atm to fire (offline option): 1) nearest neighbour, 2) bi-linear"
      logical :: fire_lsm_zcoupling = .false. ! "flag to activate reference velocity at a different height from fire_wind_height"
      real :: fire_lsm_zcoupling_ref = 50.0   ! "reference height from wich u at fire_wind_hegiht is calculated using a logarithmic profile" "m"

      real :: frac_fburnt_to_smoke = 0.02     ! "parts per unit of burned fuel becoming smoke" "g_smoke/kg_air"
      real :: fuelmc_g = 0.08                 ! Fuel moisture content ground (Dead FMC)
      real :: fuelmc_g_live = 0.30            ! Fuel moisture content ground (Live FMC). 30% Completely cured, treat as dead fuel
      real :: fuelmc_c = 1.00                 ! Fuel moisture content canopy

      logical :: fmoist_run = .false.         ! run moisture model
      integer :: fmoist_freq = 0              ! frequency to run moisture model 0: use fmoist_dt, k>0: every k timesteps
      real :: fmoist_dt = 600.0               ! moisture model time step [s]

      integer :: ideal_opt = 0                ! 0) real world, 1) ideal
      integer :: devel_opt = 0                ! 0) Standard nml options, 1) reads options in the devel nml block

        ! Objects
      integer :: fuel_opt = FUEL_ANDERSON     !  1) Anderson 13 
      integer :: ros_opt = ROS_WRFFIRE        !  0) WRF-Fire ROS
      integer :: fmc_opt = FMC_WRFFIRE        ! -1) WRF-Fire FMC
      integer :: emis_opt = EMIS_WRFFIRE      !  0) WRF-Fire emiss, 1) PM2.5 as a function of FMC. Objects to be added

        ! Ignitions
      integer :: fire_num_ignitions = 0       ! "number of ignition lines"

      real :: fire_ignition_start_lon1 = 0.0  ! "long coord of start of ignition line" "deg"
      real :: fire_ignition_start_lat1 = 0.0  ! "lat coord of start of ignition line" "deg"
      real :: fire_ignition_end_lon1 = 0.0    ! "long coord of end of ignition line" "deg"
      real :: fire_ignition_end_lat1 = 0.0    ! "lat coord of end of ignition line" "deg"
      real :: fire_ignition_ros1 = 0.01       ! "rate of spread during ignition" "m/s"
      real :: fire_ignition_start_time1 = 0.0 ! "ignition line start time" "s"
      real :: fire_ignition_end_time1 = 0.0   ! "ignition line end time" "s"
      real :: fire_ignition_radius1 = 0.0     ! "ignite all within the radius" "m"

      real :: fire_ignition_start_lon2 = 0.0
      real :: fire_ignition_start_lat2 = 0.0
      real :: fire_ignition_end_lon2 = 0.0
      real :: fire_ignition_end_lat2 = 0.0
      real :: fire_ignition_ros2 = 0.01
      real :: fire_ignition_start_time2 = 0.0
      real :: fire_ignition_end_time2 = 0.0
      real :: fire_ignition_radius2 = 0.0

      real :: fire_ignition_start_lon3 = 0.0
      real :: fire_ignition_start_lat3 = 0.0
      real :: fire_ignition_end_lon3 = 0.0
      real :: fire_ignition_end_lat3 = 0.0
      real :: fire_ignition_ros3 = 0.01
      real :: fire_ignition_start_time3 = 0.0
      real :: fire_ignition_end_time3 = 0.0
      real :: fire_ignition_radius3 = 0.0

      real :: fire_ignition_start_lon4 = 0.0
      real :: fire_ignition_start_lat4 = 0.0
      real :: fire_ignition_end_lon4 = 0.0
      real :: fire_ignition_end_lat4 = 0.0
      real :: fire_ignition_ros4 = 0.01
      real :: fire_ignition_start_time4 = 0.0
      real :: fire_ignition_end_time4 = 0.0
      real :: fire_ignition_radius4 = 0.0

      real :: fire_ignition_start_lon5 = 0.0
      real :: fire_ignition_start_lat5 = 0.0
      real :: fire_ignition_end_lon5 = 0.0
      real :: fire_ignition_end_lat5 = 0.0
      real :: fire_ignition_ros5 = 0.01
      real :: fire_ignition_start_time5 = 0.0
      real :: fire_ignition_end_time5 = 0.0
      real :: fire_ignition_radius5 = 0.0

        ! Ideal block
      real :: dx = 100.0                      ! grid spacing in x direction [m]
      real :: dy = 100.0                      ! grid spacing in y direction [m]
      integer :: nx = 100                     ! number of grid points in x direction
      integer :: ny = 100                     ! number of grid points in y direction

      real :: zonal_wind = 5.0                ! zonal wind [m s-1]
      real :: meridional_wind = 0.0           ! meridional wind [m s-1]
      integer :: fuel_cat = 1                 ! fuel type
      real :: dz_dx = 0.0                     ! slope in x direction
      real :: dz_dy = 0.0                     ! slope in y direction
      real :: elevation = 0.0                 ! Elevation [m]

      real :: cen_lat = 40.3636               ! center latitude of the grid
      real :: cen_lon = -4.4035               ! center longitude of the grid
      real :: stand_lon = -4.4035             ! standard longitude
      real :: true_lat_1 = 40.363             ! true latitude 1
      real :: true_lat_2 = 40.363             ! true latitude 2

        ! Atm block
      integer :: interval_atm = -1            ! Time step [s] of the atm (or frequency to read atm data if offline)
      integer :: kde = 1                      ! Number of atm vertical levels

        ! Devel block
      integer :: check_isolated_neg_lfn = 0   ! 0) Nothing, 1) Check for isolated negative values of the level set function
      integer :: output_level = 0             ! 0) Standard output, >0) Specialized output
    contains
      procedure, public :: Broadcast_nml => Broadcast_nml
      procedure, public :: Check_nml => Check_nml
      procedure, public :: Initialization => Init_namelist
      procedure, public :: Init_devel_block => Init_devel_block
      procedure, public :: Init_fire_block => Init_fire_block
      procedure, public :: Init_ideal_block => Init_ideal_block
      procedure, public :: Init_time_block => Init_time_block
      procedure, public :: Init_atm_block => Init_atm_block
    end type namelist_t

  contains

    subroutine Broadcast_nml (this, mpi_comm_cfbm)

      implicit none

      class (namelist_t), intent (in out) :: this
      integer, intent (in) :: mpi_comm_cfbm
#ifdef DM_PARALLEL
      type(MPI_Comm) :: mpi_comm_cfbm_f08


      call Convert_mpi_comm_to_f08 (mpi_comm_cfbm, mpi_comm_cfbm_f08)

        ! Fire block
      call Broadcast_integer (this%fire_print_msg)
      call Broadcast_real (this%fire_atm_feedback)
      call Broadcast_integer (this%fire_upwinding)
      call Broadcast_real (this%fire_viscosity)
      call Broadcast_logical (this%fire_lsm_reinit)
      call Broadcast_integer (this%fire_lsm_reinit_iter)
      call Broadcast_integer (this%fire_upwinding_reinit)
      call Broadcast_integer (this%fire_lsm_band_ngp)
      call Broadcast_real (this%reinit_pseudot_coef)
      call Broadcast_integer (this%fast_dist_reinit_opt)
      call Broadcast_integer (this%fast_dist_reinit_freq)
      call Broadcast_logical (this%fire_lsm_zcoupling)
      call Broadcast_real (this%fire_lsm_zcoupling_ref)
      call Broadcast_real (this%fire_viscosity_bg)
      call Broadcast_real (this%fire_viscosity_band)
      call Broadcast_integer (this%fire_viscosity_ngp)
      call Broadcast_logical (this%fmoist_run)
      call Broadcast_integer (this%fmoist_freq)
      call Broadcast_real (this%fmoist_dt)
      call Broadcast_real (this%fire_wind_height)
      call Broadcast_logical (this%fire_is_real_perim)
      call Broadcast_real (this%frac_fburnt_to_smoke)
      call Broadcast_real (this%fuelmc_g)
      call Broadcast_real (this%fuelmc_g_live)
      call Broadcast_real (this%fuelmc_c)

      call Broadcast_integer (this%ideal_opt)
      call Broadcast_integer (this%devel_opt)
      call Broadcast_integer (this%fuel_opt)
      call Broadcast_integer (this%ros_opt)
      call Broadcast_integer (this%emis_opt)
      call Broadcast_integer (this%wind_vinterp_opt)
      call Broadcast_integer (this%hinterp_opt)

      call Broadcast_integer (this%fire_num_ignitions)

      call Broadcast_real (this%fire_ignition_start_lon1)
      call Broadcast_real (this%fire_ignition_start_lat1)
      call Broadcast_real (this%fire_ignition_end_lon1)
      call Broadcast_real (this%fire_ignition_end_lat1)
      call Broadcast_real (this%fire_ignition_ros1)
      call Broadcast_real (this%fire_ignition_start_time1)
      call Broadcast_real (this%fire_ignition_end_time1)
      call Broadcast_real (this%fire_ignition_radius1)

      call Broadcast_real (this%fire_ignition_start_lon2)
      call Broadcast_real (this%fire_ignition_start_lat2)
      call Broadcast_real (this%fire_ignition_end_lon2)
      call Broadcast_real (this%fire_ignition_end_lat2)
      call Broadcast_real (this%fire_ignition_ros2)
      call Broadcast_real (this%fire_ignition_start_time2)
      call Broadcast_real (this%fire_ignition_end_time2)
      call Broadcast_real (this%fire_ignition_radius2)

      call Broadcast_real (this%fire_ignition_start_lon3)
      call Broadcast_real (this%fire_ignition_start_lat3)
      call Broadcast_real (this%fire_ignition_end_lon3)
      call Broadcast_real (this%fire_ignition_end_lat3)
      call Broadcast_real (this%fire_ignition_ros3)
      call Broadcast_real (this%fire_ignition_start_time3)
      call Broadcast_real (this%fire_ignition_end_time3)
      call Broadcast_real (this%fire_ignition_radius3)

      call Broadcast_real (this%fire_ignition_start_lon4)
      call Broadcast_real (this%fire_ignition_start_lat4)
      call Broadcast_real (this%fire_ignition_end_lon4)
      call Broadcast_real (this%fire_ignition_end_lat4)
      call Broadcast_real (this%fire_ignition_ros4)
      call Broadcast_real (this%fire_ignition_start_time4)
      call Broadcast_real (this%fire_ignition_end_time4)
      call Broadcast_real (this%fire_ignition_radius4)

      call Broadcast_real (this%fire_ignition_start_lon5)
      call Broadcast_real (this%fire_ignition_start_lat5)
      call Broadcast_real (this%fire_ignition_end_lon5)
      call Broadcast_real (this%fire_ignition_end_lat5)
      call Broadcast_real (this%fire_ignition_ros5)
      call Broadcast_real (this%fire_ignition_start_time5)
      call Broadcast_real (this%fire_ignition_end_time5)
      call Broadcast_real (this%fire_ignition_radius5)

        ! Ideal block
      call Broadcast_real (this%dx)
      call Broadcast_real (this%dy)

      call Broadcast_integer (this%nx)
      call Broadcast_integer (this%ny)

      call Broadcast_real (this%zonal_wind)
      call Broadcast_real (this%meridional_wind)
      call Broadcast_integer (this%fuel_cat)
      call Broadcast_real (this%dz_dx)
      call Broadcast_real (this%dz_dy)
      call Broadcast_real (this%elevation)

      call Broadcast_real (this%cen_lat)
      call Broadcast_real (this%cen_lon)
      call Broadcast_real (this%stand_lon)
      call Broadcast_real (this%true_lat_1)
      call Broadcast_real (this%true_lat_2)

        ! time block
      call Broadcast_integer (this%start_year)
      call Broadcast_integer (this%start_month)
      call Broadcast_integer (this%start_day)
      call Broadcast_integer (this%start_hour)
      call Broadcast_integer (this%start_minute)
      call Broadcast_integer (this%start_second)
      call Broadcast_integer (this%end_year)
      call Broadcast_integer (this%end_month)
      call Broadcast_integer (this%end_day)
      call Broadcast_integer (this%end_hour)
      call Broadcast_integer (this%end_minute)
      call Broadcast_integer (this%end_second)
      call Broadcast_real (this%dt)
      call Broadcast_integer (this%interval_output)

      call Broadcast_integer (this%num_tiles)
      call Broadcast_integer (this%tile_strategy)

        ! Atm block
      call Broadcast_integer (this%interval_atm)
      call Broadcast_integer (this%kde)

        ! Devel block
      call Broadcast_integer (this%check_isolated_neg_lfn)
      call Broadcast_integer (this%output_level)
    contains

      subroutine Broadcast_integer (val)

        implicit none

        integer, intent(in out) :: val

        integer :: ierr
        integer :: sendbuf(1)

        sendbuf(1) = val
        call Mpi_bcast(sendbuf, 1, MPI_INTEGER, 0, mpi_comm_cfbm_f08, ierr)
        val = sendbuf(1)

        if (ierr /= MPI_SUCCESS) &
            call Stop_simulation ('Error broadcasting integer value')

      end subroutine Broadcast_integer

      subroutine Broadcast_logical (val)

        implicit none

        logical, intent(in out) :: val

        integer :: ierr
        logical :: sendbuf(1)

        sendbuf(1) = val
        call Mpi_bcast (sendbuf, 1, MPI_LOGICAL, 0, mpi_comm_cfbm_f08, ierr)
        val = sendbuf(1)

        if (ierr /= MPI_SUCCESS) &
            call Stop_simulation ('Error broadcasting logical value')

      end subroutine Broadcast_logical

      subroutine Broadcast_real (val)

        implicit none

        real, intent(in out) :: val

        integer :: ierr
        real    :: sendbuf(1)

        sendbuf(1) = val
        call Mpi_bcast (sendbuf, 1, MPI_REAL, 0, mpi_comm_cfbm_f08, ierr)
        val = sendbuf(1)

        if (ierr /= MPI_SUCCESS) &
            call Stop_simulation ('Error broadcasting real value')

      end subroutine Broadcast_real

#endif

    end subroutine Broadcast_nml

    subroutine Check_nml (this)

      implicit none

      class (namelist_t), intent (in out) :: this

      if (this%ideal_opt /= 0 .and. this%fmoist_run) &
          call Stop_simulation ('ideal runs do not support a FMC model')

    end subroutine Check_nml

    subroutine Init_atm_block (this, file_name)

      implicit none

      class (namelist_t), intent (in out) :: this
      character (len = *), intent (in) :: file_name

      integer :: kde, interval_atm
      integer :: unit_nml, io_stat
      character (len = :), allocatable :: msg

      namelist /atm/ kde, interval_atm


      interval_atm = this%interval_atm
      kde = this%kde

      open (newunit = unit_nml, file = trim (file_name), action = 'read', iostat = io_stat)
      if (io_stat /= 0) then
        msg = 'Problems opening namelist file ' // trim (file_name)
        call Stop_simulation (msg)
      end if

      read (unit_nml, nml = atm, iostat = io_stat)
      if (io_stat /= 0) call Stop_simulation ('Problems reading namelist atm block')
      close (unit_nml)

      this%interval_atm = interval_atm
      this%kde = kde

    end subroutine Init_atm_block

    subroutine Init_devel_block (this, file_name)

      implicit none

      class (namelist_t), intent (in out) :: this
      character (len = *), intent (in) :: file_name

      integer :: check_isolated_neg_lfn, output_level
      integer :: unit_nml, io_stat
      character (len = :), allocatable :: msg

      namelist /devel/ check_isolated_neg_lfn, output_level


      check_isolated_neg_lfn = this%check_isolated_neg_lfn
      output_level = this%output_level

      open (newunit = unit_nml, file = trim (file_name), action = 'read', iostat = io_stat)
      if (io_stat /= 0) then
        msg = 'Problems opening namelist file ' // trim (file_name)
        call Stop_simulation (msg)
      end if

      read (unit_nml, nml = devel, iostat = io_stat)
      if (io_stat /= 0) call Stop_simulation ('Problems reading namelist devel block')
      close (unit_nml)

      this%check_isolated_neg_lfn = check_isolated_neg_lfn
      this%output_level = output_level

    end subroutine Init_devel_block

    subroutine Init_fire_block (this, file_name)

      implicit none

      class (namelist_t), intent (in out) :: this
      character (len = *), intent (in) :: file_name

      integer :: fire_print_msg, fire_upwinding, fire_lsm_reinit_iter, fire_upwinding_reinit, fire_lsm_band_ngp, &
          fast_dist_reinit_opt, fast_dist_reinit_freq, fire_viscosity_ngp, wind_vinterp_opt, hinterp_opt, ideal_opt, devel_opt, &
          fuel_opt, ros_opt, fmc_opt, emis_opt, fmoist_freq
      real :: fire_atm_feedback, fire_viscosity, fire_lsm_zcoupling_ref, fire_viscosity_bg, fire_viscosity_band, &
          fmoist_dt, fire_wind_height, frac_fburnt_to_smoke, fuelmc_g, fuelmc_g_live, fuelmc_c, reinit_pseudot_coef
      logical :: fire_lsm_reinit, fire_lsm_zcoupling, fmoist_run, fire_is_real_perim

        ! ignitions
      integer :: fire_num_ignitions
      real :: fire_ignition_start_lon1, fire_ignition_start_lat1, fire_ignition_end_lon1, fire_ignition_end_lat1, &
          fire_ignition_ros1, fire_ignition_start_time1, fire_ignition_end_time1, fire_ignition_radius1
      real :: fire_ignition_start_lon2, fire_ignition_start_lat2, fire_ignition_end_lon2, fire_ignition_end_lat2, &
          fire_ignition_ros2, fire_ignition_start_time2, fire_ignition_end_time2, fire_ignition_radius2
      real :: fire_ignition_start_lon3, fire_ignition_start_lat3, fire_ignition_end_lon3, fire_ignition_end_lat3, &
          fire_ignition_ros3, fire_ignition_start_time3, fire_ignition_end_time3, fire_ignition_radius3
      real :: fire_ignition_start_lon4, fire_ignition_start_lat4, fire_ignition_end_lon4, fire_ignition_end_lat4, &
          fire_ignition_ros4, fire_ignition_start_time4, fire_ignition_end_time4, fire_ignition_radius4
      real :: fire_ignition_start_lon5, fire_ignition_start_lat5, fire_ignition_end_lon5, fire_ignition_end_lat5, &
          fire_ignition_ros5, fire_ignition_start_time5, fire_ignition_end_time5, fire_ignition_radius5

      namelist /fire/  fire_print_msg, fire_atm_feedback, fire_upwinding, fire_viscosity, fire_lsm_reinit, &
          fast_dist_reinit_opt, fast_dist_reinit_freq, fire_lsm_reinit_iter, fire_upwinding_reinit, &
          fire_lsm_band_ngp, fire_lsm_zcoupling, fire_lsm_zcoupling_ref, fire_viscosity_bg, fire_viscosity_band, &
          fire_viscosity_ngp, fmoist_run, fmoist_freq, fmoist_dt, fire_wind_height, fire_is_real_perim, &
          frac_fburnt_to_smoke, fuelmc_g, fuelmc_g_live, fuelmc_c, ideal_opt, devel_opt, fuel_opt, ros_opt, fmc_opt, emis_opt, &
          wind_vinterp_opt, hinterp_opt, reinit_pseudot_coef, &
            ! Ignitions
          fire_num_ignitions, &
            ! Ignition 1
          fire_ignition_start_lon1, fire_ignition_start_lat1, fire_ignition_end_lon1, fire_ignition_end_lat1, &
          fire_ignition_ros1, fire_ignition_start_time1, fire_ignition_end_time1, fire_ignition_radius1, &
            ! Ignition 2
          fire_ignition_start_lon2, fire_ignition_start_lat2, fire_ignition_end_lon2, fire_ignition_end_lat2, &
          fire_ignition_ros2, fire_ignition_start_time2, fire_ignition_end_time2, fire_ignition_radius2, &
            ! Ignition 3
          fire_ignition_start_lon3, fire_ignition_start_lat3, fire_ignition_end_lon3, fire_ignition_end_lat3, &
          fire_ignition_ros3, fire_ignition_start_time3, fire_ignition_end_time3, fire_ignition_radius3, &
            ! Ignition 4
          fire_ignition_start_lon4, fire_ignition_start_lat4, fire_ignition_end_lon4, fire_ignition_end_lat4, &
          fire_ignition_ros4, fire_ignition_start_time4, fire_ignition_end_time4, fire_ignition_radius4, &
            ! Ignition 5
          fire_ignition_start_lon5, fire_ignition_start_lat5, fire_ignition_end_lon5, fire_ignition_end_lat5, &
          fire_ignition_ros5, fire_ignition_start_time5, fire_ignition_end_time5, fire_ignition_radius5

      integer :: unit_nml, io_stat
      character (len = :), allocatable :: msg


        ! Set default values
      fire_print_msg = this%fire_print_msg
      fire_atm_feedback = this%fire_atm_feedback
      fire_upwinding = this%fire_upwinding
      fire_viscosity = this%fire_viscosity
      fire_lsm_reinit = this%fire_lsm_reinit
      fire_lsm_reinit_iter = this%fire_lsm_reinit_iter
      fire_upwinding_reinit = this%fire_upwinding_reinit
      fire_lsm_band_ngp = this%fire_lsm_band_ngp
      reinit_pseudot_coef = this%reinit_pseudot_coef
      fast_dist_reinit_opt = this%fast_dist_reinit_opt
      fast_dist_reinit_freq = this%fast_dist_reinit_freq
      fire_lsm_zcoupling = this%fire_lsm_zcoupling
      fire_lsm_zcoupling_ref = this%fire_lsm_zcoupling_ref
      fire_viscosity_bg = this%fire_viscosity_bg
      fire_viscosity_band = this%fire_viscosity_band
      fire_viscosity_ngp = this%fire_viscosity_ngp
      fmoist_run = this%fmoist_run
      fmoist_freq = this%fmoist_freq
      fmoist_dt = this%fmoist_dt
      fire_wind_height = this%fire_wind_height
      fire_is_real_perim = this%fire_is_real_perim
      frac_fburnt_to_smoke = this%frac_fburnt_to_smoke
      fuelmc_g = this%fuelmc_g
      fuelmc_g_live = this%fuelmc_g_live
      fuelmc_c = this%fuelmc_c

      ideal_opt = this%ideal_opt
      devel_opt = this%devel_opt

      fuel_opt = this%fuel_opt
      ros_opt = this%ros_opt
      fmc_opt = this%fmc_opt
      emis_opt = this%emis_opt
      wind_vinterp_opt = this%wind_vinterp_opt
      hinterp_opt = this%hinterp_opt

      fire_num_ignitions = this%fire_num_ignitions

      fire_ignition_start_lon1 = this%fire_ignition_start_lon1
      fire_ignition_start_lat1 = this%fire_ignition_start_lat1
      fire_ignition_end_lon1 = this%fire_ignition_end_lon1
      fire_ignition_end_lat1 = this%fire_ignition_end_lat1
      fire_ignition_ros1 = this%fire_ignition_ros1
      fire_ignition_start_time1 = this%fire_ignition_start_time1
      fire_ignition_end_time1 = this%fire_ignition_end_time1
      fire_ignition_radius1 = this%fire_ignition_radius1

      fire_ignition_start_lon2 = this%fire_ignition_start_lon2
      fire_ignition_start_lat2 = this%fire_ignition_start_lat2
      fire_ignition_end_lon2 = this%fire_ignition_end_lon2
      fire_ignition_end_lat2 = this%fire_ignition_end_lat2
      fire_ignition_ros2 = this%fire_ignition_ros2
      fire_ignition_start_time2 = this%fire_ignition_start_time2
      fire_ignition_end_time2 = this%fire_ignition_end_time2
      fire_ignition_radius2 = this%fire_ignition_radius2

      fire_ignition_start_lon3 = this%fire_ignition_start_lon3
      fire_ignition_start_lat3 = this%fire_ignition_start_lat3
      fire_ignition_end_lon3 = this%fire_ignition_end_lon3
      fire_ignition_end_lat3 = this%fire_ignition_end_lat3
      fire_ignition_ros3 = this%fire_ignition_ros3
      fire_ignition_start_time3 = this%fire_ignition_start_time3
      fire_ignition_end_time3 = this%fire_ignition_end_time3
      fire_ignition_radius3 = this%fire_ignition_radius3

      fire_ignition_start_lon4 = this%fire_ignition_start_lon4
      fire_ignition_start_lat4 = this%fire_ignition_start_lat4
      fire_ignition_end_lon4 = this%fire_ignition_end_lon4
      fire_ignition_end_lat4 = this%fire_ignition_end_lat4
      fire_ignition_ros4 = this%fire_ignition_ros4
      fire_ignition_start_time4 = this%fire_ignition_start_time4
      fire_ignition_end_time4 = this%fire_ignition_end_time4
      fire_ignition_radius4 = this%fire_ignition_radius4

      fire_ignition_start_lon5 = this%fire_ignition_start_lon5
      fire_ignition_start_lat5 = this%fire_ignition_start_lat5
      fire_ignition_end_lon5 = this%fire_ignition_end_lon5
      fire_ignition_end_lat5 = this%fire_ignition_end_lat5
      fire_ignition_ros5 = this%fire_ignition_ros5
      fire_ignition_start_time5 = this%fire_ignition_start_time5
      fire_ignition_end_time5 = this%fire_ignition_end_time5
      fire_ignition_radius5 = this%fire_ignition_radius5

        ! Read namelist
      open (newunit = unit_nml, file = trim (file_name), action = 'read', iostat = io_stat)
      if (io_stat /= 0) then
        msg = 'Problems opening namelist file ' // trim (file_name)
        call Stop_simulation (msg)
      end if

      read (unit_nml, nml = fire, iostat = io_stat)
      if (io_stat /= 0) call Stop_simulation ('Problems reading namelist fire block')

      close (unit_nml)

        ! Assign namelist values
      this%fire_print_msg = fire_print_msg
      this%fire_atm_feedback = fire_atm_feedback
      this%fire_upwinding = fire_upwinding
      this%fire_viscosity = fire_viscosity
      this%fire_lsm_reinit = fire_lsm_reinit
      this%fire_lsm_reinit_iter = fire_lsm_reinit_iter
      this%fire_upwinding_reinit = fire_upwinding_reinit
      this%fire_lsm_band_ngp = fire_lsm_band_ngp
      this%reinit_pseudot_coef = reinit_pseudot_coef
      this%fast_dist_reinit_opt = fast_dist_reinit_opt
      this%fast_dist_reinit_freq = fast_dist_reinit_freq
      this%fire_lsm_zcoupling = fire_lsm_zcoupling
      this%fire_lsm_zcoupling_ref = fire_lsm_zcoupling_ref
      this%fire_viscosity_bg = fire_viscosity_bg
      this%fire_viscosity_band = fire_viscosity_band
      this%fire_viscosity_ngp = fire_viscosity_ngp
      this%fmoist_run = fmoist_run
      this%fmoist_freq = fmoist_freq
      this%fmoist_dt = fmoist_dt
      this%fire_wind_height = fire_wind_height
      this%fire_is_real_perim = fire_is_real_perim
      this%frac_fburnt_to_smoke = frac_fburnt_to_smoke
      this%fuelmc_g = fuelmc_g
      this%fuelmc_g_live = fuelmc_g_live
      this%fuelmc_c = fuelmc_c

      this%ideal_opt = ideal_opt
      this%devel_opt = devel_opt

      this%fuel_opt = fuel_opt
      this%ros_opt = ros_opt
      this%fmc_opt = fmc_opt
      this%emis_opt = emis_opt
      this%wind_vinterp_opt = wind_vinterp_opt
      this%hinterp_opt = hinterp_opt

      this%fire_num_ignitions = fire_num_ignitions

      this%fire_ignition_start_lon1 = fire_ignition_start_lon1
      this%fire_ignition_start_lat1 = fire_ignition_start_lat1
      this%fire_ignition_end_lon1 = fire_ignition_end_lon1
      this%fire_ignition_end_lat1 = fire_ignition_end_lat1
      this%fire_ignition_ros1 = fire_ignition_ros1
      this%fire_ignition_start_time1 = fire_ignition_start_time1
      this%fire_ignition_end_time1 = fire_ignition_end_time1
      this%fire_ignition_radius1 = fire_ignition_radius1

      this%fire_ignition_start_lon2 = fire_ignition_start_lon2
      this%fire_ignition_start_lat2 = fire_ignition_start_lat2
      this%fire_ignition_end_lon2 = fire_ignition_end_lon2
      this%fire_ignition_end_lat2 = fire_ignition_end_lat2
      this%fire_ignition_ros2 = fire_ignition_ros2
      this%fire_ignition_start_time2 = fire_ignition_start_time2
      this%fire_ignition_end_time2 = fire_ignition_end_time2
      this%fire_ignition_radius2 = fire_ignition_radius2

      this%fire_ignition_start_lon3 = fire_ignition_start_lon3
      this%fire_ignition_start_lat3 = fire_ignition_start_lat3
      this%fire_ignition_end_lon3 = fire_ignition_end_lon3
      this%fire_ignition_end_lat3 = fire_ignition_end_lat3
      this%fire_ignition_ros3 = fire_ignition_ros3
      this%fire_ignition_start_time3 = fire_ignition_start_time3
      this%fire_ignition_end_time3 = fire_ignition_end_time3
      this%fire_ignition_radius3 = fire_ignition_radius3

      this%fire_ignition_start_lon4 = fire_ignition_start_lon4
      this%fire_ignition_start_lat4 = fire_ignition_start_lat4
      this%fire_ignition_end_lon4 = fire_ignition_end_lon4
      this%fire_ignition_end_lat4 = fire_ignition_end_lat4
      this%fire_ignition_ros4 = fire_ignition_ros4
      this%fire_ignition_start_time4 = fire_ignition_start_time4
      this%fire_ignition_end_time4 = fire_ignition_end_time4
      this%fire_ignition_radius4 = fire_ignition_radius4

      this%fire_ignition_start_lon5 = fire_ignition_start_lon5
      this%fire_ignition_start_lat5 = fire_ignition_start_lat5
      this%fire_ignition_end_lon5 = fire_ignition_end_lon5
      this%fire_ignition_end_lat5 = fire_ignition_end_lat5
      this%fire_ignition_ros5 = fire_ignition_ros5
      this%fire_ignition_start_time5 = fire_ignition_start_time5
      this%fire_ignition_end_time5 = fire_ignition_end_time5
      this%fire_ignition_radius5 = fire_ignition_radius5

    end subroutine Init_fire_block

    subroutine Init_ideal_block (this, file_name)

      implicit none

      class (namelist_t), intent (in out) :: this
      character (len = *), intent (in) :: file_name

      real :: dx, dy, zonal_wind, meridional_wind, cen_lat, cen_lon, stand_lon, true_lat_1, true_lat_2, &
          dz_dx, dz_dy, elevation
      integer :: nx, ny, fuel_cat

      character (len = :), allocatable :: msg
      integer :: unit_nml, io_stat

      namelist /ideal/ dx, dy, nx, ny, zonal_wind, meridional_wind, fuel_cat, &
          dz_dx, dz_dy, elevation, cen_lat, cen_lon, stand_lon, true_lat_1, true_lat_2


        ! Set default values
      dx = this%dx
      dy = this%dy
      nx = this%nx
      ny = this%ny

      zonal_wind = this%zonal_wind
      meridional_wind = this%meridional_wind
      fuel_cat = this%fuel_cat
      dz_dx = this%dz_dx
      dz_dy = this%dz_dy
      elevation = this%elevation

      cen_lat = this%cen_lat
      cen_lon = this%cen_lon
      stand_lon = this%stand_lon
      true_lat_1 = this%true_lat_1
      true_lat_2 = this%true_lat_2

        ! Read namelist
      open (newunit = unit_nml, file = trim (file_name), action = 'read', iostat = io_stat)
      if (io_stat /= 0) then
        msg = 'Problems opening namelist file ' // trim (file_name)
        call Stop_simulation (msg)
      end if

      read (unit_nml, nml = ideal, iostat = io_stat)
      if (io_stat /= 0) call Stop_simulation ('Problems reading namelist ideal block')

      close (unit_nml)

        ! Assign namelist values
      this%dx = dx
      this%dy = dy
      this%nx = nx
      this%ny = ny

      this%zonal_wind = zonal_wind
      this%meridional_wind = meridional_wind
      this%fuel_cat = fuel_cat
      this%dz_dx = dz_dx
      this%dz_dy = dz_dy
      this%elevation = elevation

      this%cen_lat = cen_lat
      this%cen_lon = cen_lon
      this%stand_lon = stand_lon
      this%true_lat_1 = true_lat_1
      this%true_lat_2 = true_lat_2

    end subroutine Init_ideal_block

    subroutine Init_time_block (this, file_name)

      implicit none

      class (namelist_t), intent (in out) :: this
      character (len = *), intent (in) :: file_name

      integer :: start_year, start_month, start_day, start_hour, start_minute, start_second, &
          end_year, end_month, end_day, end_hour, end_minute, end_second, interval_output, &
          num_tiles, tile_strategy
      real :: dt

      character (len = :), allocatable :: msg
      integer :: unit_nml, io_stat

      namelist /time/ start_year, start_month, start_day, start_hour, start_minute, start_second, &
          end_year, end_month, end_day, end_hour, end_minute, end_second, dt, interval_output, &
          num_tiles


        ! Set default values
      start_year = this%start_year
      start_month = this%start_month
      start_day = this%start_day
      start_hour = this%start_hour
      start_minute = this%start_minute
      start_second = this%start_second
      end_year = this%end_year
      end_month = this%end_month
      end_day = this%end_day
      end_hour = this%end_hour
      end_minute = this%end_minute
      end_second = this%end_second
      dt = this%dt
      interval_output = this%interval_output

      num_tiles = this%num_tiles
      tile_strategy = this%tile_strategy

        ! Read namelist
      open (newunit = unit_nml, file = trim (file_name), action = 'read', iostat = io_stat)
      if (io_stat /= 0) then
        msg = 'Problems opening namelist file ' // trim (file_name)
        call Stop_simulation (msg)
      end if

      read (unit_nml, nml = time, iostat = io_stat)
      if (io_stat /= 0) call Stop_simulation ('Problems reading namelist time block')

      close (unit_nml)

        ! Set namelist values
      this%start_year = start_year
      this%start_month = start_month
      this%start_day = start_day
      this%start_hour = start_hour
      this%start_minute = start_minute
      this%start_second = start_second
      this%end_year = end_year
      this%end_month = end_month
      this%end_day = end_day
      this%end_hour = end_hour
      this%end_minute = end_minute
      this%end_second = end_second
      this%dt = dt
      this%interval_output = interval_output

      this%num_tiles = num_tiles
      this%tile_strategy = tile_strategy

    end subroutine Init_time_block

    subroutine Init_namelist (this, file_name)

      implicit none

      class (namelist_t), intent (out) :: this
      character (len = *), intent (in) :: file_name

      logical, parameter :: DEBUG_LOCAL = .false.


      if (DEBUG_LOCAL) call Print_message ('  Entering subroutine Read_namelist')

      call this%Init_time_block (file_name = trim (file_name))
      call this%Init_fire_block (file_name = trim (file_name))
      call this%Init_atm_block (file_name = trim (file_name))
      if (this%ideal_opt > 0) call this%Init_ideal_block (file_name = trim (file_name))
      if (this%devel_opt > 0) call this%Init_devel_block (file_name = trim (file_name))

      call this%Check_nml ()

      if (DEBUG_LOCAL) call Print_message ('  Leaving subroutine Read_namelist')

    end subroutine Init_namelist

  end module namelist_mod
