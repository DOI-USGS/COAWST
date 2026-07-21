  module wrfdata_mod

    use datetime_mod, only : datetime_t
    use namelist_mod, only : namelist_t
    use netcdf_mod, only : Get_netcdf_var, Get_netcdf_att, Get_netcdf_dim, Is_netcdf_file_present
    use proj_lc_mod, only : proj_lc_t
    use stderrout_mod, only : Print_message, Stop_simulation
    use interp_mod, only : VINTERP_WINDS_FROM_3D_WINDS, VINTERP_WINDS_FROM_10M_WINDS
    use coupling_mod, only : Interp_horizontal, Calc_fire_wind

    implicit none

    private

    public :: wrfdata_t, G, RERADIUS

    real, parameter :: G = 9.81                   ! acceleration due to gravity [m s-2]
    real, parameter :: RERADIUS = 1.0 / 6370.0e03 ! reciprocal of earth radius (m^-1)

    type :: wrfdata_t
      character (len = 300) :: file_name
      real, dimension(:, :, :), allocatable :: u3d, v3d, phl
      real, dimension(:, :), allocatable :: lats, lons, lats_c, lons_c, t2, q2, z0, psfc, rain, ua, va, u10, v10
      integer :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte
      real :: cen_lat, cen_lon, dx, dy, truelat1, truelat2, stand_lon
    contains
      procedure, public :: Destroy_phl => Destroy_geopotential_levels
      procedure, public :: Destroy_psfc => Destroy_surface_pressure
      procedure, public :: Destroy_rain => Destroy_rain
      procedure, public :: Destroy_q2 => Destroy_specific_humidity_2m
      procedure, public :: Destroy_t2 => Destroy_temperature_2m
      procedure, public :: Destroy_u10 => Destroy_zonal_wind10
      procedure, public :: Destroy_v10 => Destroy_meridional_wind10
      procedure, public :: Destroy_u3d => Destroy_zonal_wind
      procedure, public :: Destroy_v3d => Destroy_meridional_wind
      procedure, public :: Destroy_z0 => Destroy_z0
      procedure, public :: Get_datetime_index => Get_datetime_index
      procedure, public :: Get_latlons => Get_latlons
      procedure, public :: Get_latcloncs => Get_latcloncs
      procedure, public :: Get_phl => Get_geopotential_levels
      procedure, public :: Get_projection => Get_projection
      procedure, public :: Get_rain => Get_rain
      procedure, public :: Get_psfc => Get_surface_pressure
      procedure, public :: Get_q2 => Get_specific_humidity_2m
      procedure, public :: Get_t2 => Get_temperature_2m
      procedure, public :: Get_u10 => Get_zonal_wind10
      procedure, public :: Get_u3d => Get_zonal_wind_3d
      procedure, public :: Get_u3d_stag => Get_zonal_wind_stag_3d
      procedure, public :: Get_v10 => Get_meridional_wind10
      procedure, public :: Get_v3d => Get_meridional_wind_3d
      procedure, public :: Get_v3d_stag => Get_meridional_wind_stag_3d
      procedure, public :: Get_z0 => Get_z0
      procedure, public :: Interp_var2grid => Interp_var2grid
      procedure, public :: Print_domain => Print_domain
      procedure, public :: Update_atm_state => Update_atm_state
    end type wrfdata_t

    interface wrfdata_t
      module procedure Wrfdata_t_const
    end interface wrfdata_t

  contains

    subroutine Destroy_geopotential_levels (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%phl)) deallocate (this%phl)

    end subroutine Destroy_geopotential_levels

    subroutine Destroy_meridional_wind (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%v3d)) deallocate (this%v3d)

    end subroutine Destroy_meridional_wind

    subroutine Destroy_meridional_wind10 (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%v10)) deallocate (this%v10)

    end subroutine Destroy_meridional_wind10

    subroutine Destroy_rain (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%rain)) deallocate (this%rain)

    end subroutine Destroy_rain

    subroutine Destroy_surface_pressure (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%psfc)) deallocate (this%psfc)

    end subroutine Destroy_surface_pressure

    subroutine Destroy_specific_humidity_2m (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%q2)) deallocate (this%q2)

    end subroutine Destroy_specific_humidity_2m

    subroutine Destroy_temperature_2m (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%t2)) deallocate (this%t2)

    end subroutine Destroy_temperature_2m

    subroutine Destroy_z0 (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%z0)) deallocate (this%z0)

    end subroutine Destroy_z0

    subroutine Destroy_zonal_wind (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%u3d)) deallocate (this%u3d)

    end subroutine Destroy_zonal_wind

    subroutine Destroy_zonal_wind10 (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      if (allocated(this%u10)) deallocate (this%u10)

    end subroutine Destroy_zonal_wind10

    function Get_datetime_index (this, datetime) result (return_value)

      use, intrinsic :: iso_fortran_env, only : ERROR_UNIT

      implicit none

      class (wrfdata_t), intent (in) :: this
      type (datetime_t), intent (in) :: datetime
      integer :: return_value

      character (len = :), dimension (:), allocatable :: times
      integer :: n


      call Get_netcdf_var (trim (this%file_name), 'Times', times)

      do n = 1, size (times)
        if (times(n) == datetime%datetime) then
          return_value = n
          return
        end if
      end do

      write (ERROR_UNIT, *) 'Datetime not found in WRF file:'
      call datetime%Print_datetime ()
      stop

    end function Get_datetime_index

    subroutine Get_latcloncs (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      type (proj_lc_t) :: proj
      logical, parameter :: OUTPUT_LATLON_CHECK = .false.
      integer :: nx, ny, i, j


      nx = size (this%lats, dim = 1) + 1
      ny = size (this%lats, dim = 2) + 1

      allocate (this%lats_c(nx, ny))
      allocate (this%lons_c(nx, ny))

      proj = this%Get_projection (stagger = .true.)

      do j = 1, ny
        do i = 1, nx
          call proj%Calc_latlon (i = real (i), j = real (j), lat = this%lats_c(i, j), lon = this%lons_c(i, j))
        end do
      end do

      if (OUTPUT_LATLON_CHECK) call Write_latlon_check ()

    contains

      subroutine Write_latlon_check ()

        implicit none

        integer :: i, j, unit1, unit2


        open (newunit = unit1, file = "latlons_wrf_mass.dat")
        do j = 1,  size (this%lats, dim = 2)
          do i = 1,  size (this%lats, dim = 1)
            write (unit1, *) i, j, this%lats(i, j), this%lons(i, j)
          end do
        end do
        close (unit1)

        open (newunit = unit2, file = "latlons_wrf_corners_estimated.dat")
        do j = 1,  size (this%lats, dim = 2) + 1
          do i = 1,  size (this%lats, dim = 1) + 1
            write (unit2, *) i, j, this%lats_c(i, j), this%lons_c(i, j)
          end do
        end do
        close (unit2)

      end subroutine Write_latlon_check

    end subroutine Get_latcloncs

    subroutine Get_latlons (this)

      implicit none

      class (wrfdata_t), intent (in out) :: this

      real, dimension(:, :, :), allocatable :: var3d


      call Get_netcdf_var (trim (this%file_name), 'XLAT', var3d)
      this%lats = var3d(:, :, 1)
      deallocate (var3d)

      call Get_netcdf_var (trim (this%file_name), 'XLONG', var3d)
      this%lons = var3d(:, :, 1)
      deallocate (var3d)

    end subroutine Get_latlons

    subroutine Get_geopotential_levels (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :, :), allocatable :: var4d, var4d2
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'PH', var4d)
      call Get_netcdf_var (trim (this%file_name), 'PHB', var4d2)

      this%phl = var4d(:, :, :, nt) + var4d2(:, :, :, nt)
      deallocate (var4d, var4d2)

    end subroutine Get_geopotential_levels

    subroutine Get_meridional_wind10 (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :), allocatable :: var3d
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'V10', var3d)
      this%v10 = var3d(:, :, nt)

    end subroutine Get_meridional_wind10

    subroutine Get_meridional_wind_3d (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :, :), allocatable :: var4d
      integer :: nt, nmass


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'V', var4d)
      nmass = size (var4d, dim = 2)
      this%v3d = 0.5 * (var4d(:, 1:nmass - 1, :, nt) + var4d(:, 2:nmass, :, nt))
      deallocate (var4d)

    end subroutine Get_meridional_wind_3d

    subroutine Get_meridional_wind_stag_3d (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :, :), allocatable :: var4d
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'V', var4d)
      this%v3d = var4d(:, :, :, nt)
      deallocate (var4d)

    end subroutine Get_meridional_wind_stag_3d

    function Get_projection (this, stagger) result (return_value)

      use, intrinsic :: iso_fortran_env, only : ERROR_UNIT

      implicit none

      class (wrfdata_t), intent (in) :: this
      logical, intent (in), optional :: stagger
      type (proj_lc_t) :: return_value

      integer :: nx, ny, offset


      offset = 0
      if (present (stagger)) then
        offset = 1
      end if

      if (allocated (this%lats)) then
        nx = size (this%lats, dim = 1) + offset
        ny = size (this%lats, dim = 2) + offset
      else
        write (ERROR_UNIT, *) 'lats array needs to be initialized to get the WRF projection'
      end if

      return_value = proj_lc_t (cen_lat = this%cen_lat , cen_lon = this%cen_lon, &
          dx = this%dx, dy = this%dy, standard_lon = this%stand_lon , true_lat_1 = this%truelat1 , true_lat_2 = this%truelat2, &
          nx = nx, ny = ny)

    end function Get_projection

    subroutine Get_specific_humidity_2m (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :), allocatable :: var3d
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'Q2', var3d)
      this%q2 = var3d(:, :, nt)

    end subroutine Get_specific_humidity_2m

    subroutine Get_rain (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :), allocatable :: var3d, var3d2
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'RAINC', var3d)
      call Get_netcdf_var (trim (this%file_name), 'RAINNC', var3d2)
      this%rain = var3d(:, :, nt) + var3d2(:, :, nt)
      deallocate (var3d, var3d2)

    end subroutine Get_rain

    subroutine Get_surface_pressure (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :), allocatable :: var3d
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'PSFC', var3d)
      this%psfc = var3d(:, :, nt)

    end subroutine Get_surface_pressure

    subroutine Get_temperature_2m (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :), allocatable :: var3d
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'T2', var3d)
      this%t2 = var3d(:, :, nt)

    end subroutine Get_temperature_2m

    subroutine Get_z0 (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :), allocatable :: var3d
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'ZNT', var3d)
      this%z0 = var3d(:, :, nt)

    end subroutine Get_z0

    subroutine Get_zonal_wind10 (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :), allocatable :: var3d
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'U10', var3d)
      this%u10 = var3d(:, :, nt)

    end subroutine Get_zonal_wind10

    subroutine Get_zonal_wind_3d (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :, :), allocatable :: var4d
      integer :: nt, nmass


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'U', var4d)
      nmass = size (var4d, dim = 1)
      this%u3d = 0.5 * (var4d(1:nmass - 1, :, :, nt) + var4d(2:nmass, :, :, nt))
      deallocate (var4d)

    end subroutine Get_zonal_wind_3d

    subroutine Get_zonal_wind_stag_3d (this, datetime)

      implicit none

      class (wrfdata_t), intent (in out) :: this
      type (datetime_t), intent (in) :: datetime

      real, dimension(:, :, :, :), allocatable :: var4d
      integer :: nt


      nt = this%Get_datetime_index (datetime)
      call Get_netcdf_var (trim (this%file_name), 'U', var4d)
      this%u3d = var4d(:, :, :, nt)
      deallocate (var4d)

    end subroutine Get_zonal_wind_stag_3d

    subroutine Interp_var2grid (this, lats_out, lons_out, ifms, ifme, jfms, jfme, &
        num_tiles, i_start, i_end, j_start, j_end, var_name, hinterp_opt, data_out)

      implicit none

      class (wrfdata_t), intent(in) :: this
      integer, intent (in) :: ifms, ifme, jfms, jfme, num_tiles
      integer, dimension (num_tiles), intent (in) :: i_start, i_end, j_start, j_end
      real, dimension(ifms:ifme, jfms:jfme), intent(in) :: lats_out, lons_out
      character (len = *), intent (in) :: var_name
      integer, intent (in) :: hinterp_opt
      real, dimension(ifms:ifme, jfms:jfme), intent(in out) :: data_out

      real, dimension(:, :), allocatable :: var_wrf
      type (proj_lc_t) :: proj
      integer :: ims, ime, jms, jme


        ! Init
      proj = this%Get_projection ()

        ! Get WRF data
      select case (var_name)
        case ('t2')
          var_wrf = this%t2

        case ('q2')
          var_wrf = this%q2

        case ('psfc')
          var_wrf = this%psfc

        case ('rain')
          var_wrf = this%rain

        case ('ua')
          var_wrf = this%ua

        case ('va')
          var_wrf = this%va

        case default
          call Stop_simulation ('Unknown variable name to interpolate')

      end select

        ! Interpolate
      ims = 1
      ime = size (var_wrf, dim = 1)
      jms = 1
      jme = size (var_wrf, dim = 2)
      call Interp_horizontal (var_wrf, proj, ims, ime, jms, jme, ifms, ifme, jfms, jfme, &
          num_tiles, i_start, i_end, j_start, j_end, hinterp_opt, lats_out, lons_out, data_out)

    end subroutine Interp_var2grid

    subroutine Print_domain (this)

      use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

      implicit none

      class (wrfdata_t), intent(in out) :: this


      write (OUTPUT_UNIT, *) ''
      write (OUTPUT_UNIT, *) 'ids = ', this%ids, 'ide = ', this%ide
      write (OUTPUT_UNIT, *) 'jds = ', this%jds, 'jde = ', this%jde
      write (OUTPUT_UNIT, *) 'kds = ', this%kds, 'kde = ', this%kde

      write (OUTPUT_UNIT, *) 'ims = ', this%ims, 'ime = ', this%ime
      write (OUTPUT_UNIT, *) 'jms = ', this%jms, 'jme = ', this%jme
      write (OUTPUT_UNIT, *) 'kms = ', this%kms, 'kme = ', this%kme

      write (OUTPUT_UNIT, *) 'its = ', this%its, 'ite = ', this%ite
      write (OUTPUT_UNIT, *) 'jts = ', this%jts, 'jte = ', this%jte
      write (OUTPUT_UNIT, *) 'kts = ', this%kts, 'kte = ', this%kte

    end subroutine Print_domain

    subroutine Update_atm_state (this, datetime_now, config_flags)

      implicit none

      class (wrfdata_t), intent(in out) :: this
      type (datetime_t), intent (in) :: datetime_now
      type (namelist_t), intent (in) :: config_flags

      integer :: iims, iime, jims, jime, kims, kime, ioms, iome, joms, jome, iops, iope, jops, jope


      call this%Get_t2 (datetime_now)
      call this%Get_q2 (datetime_now)
      call this%Get_psfc (datetime_now)
      call this%Get_rain (datetime_now)

      select case (config_flags%wind_vinterp_opt)
        case (VINTERP_WINDS_FROM_3D_WINDS)
          call this%Get_z0 (datetime_now)
          call this%Get_u3d (datetime_now)
          call this%Get_v3d (datetime_now)
          call this%Get_phl (datetime_now)

            ! Set input (i) and output (o) indices
          iims = this%ids
          iime = this%ide - 1
          jims = this%ids
          jime = this%ide - 1
          kims = this%kds
          kime = this%kde - 1

          ioms = this%ids
          iome = this%ide - 1
          joms = this%ids
          jome = this%ide - 1

          iops = this%ids
          iope = this%ide - 1
          jops = this%ids
          jope = this%ide - 1
                                                   ! For compatibility with nuopc couplings
                                                   ! pass z_at_w with vertical dim kde - 1 instead of kde
          call Calc_fire_wind (this%u3d, this%v3d, this%phl(iims:iime, jims:jime, kims:kime) / G, this%z0, &
              iims, iime, jims, jime, kims, kime, config_flags%fire_lsm_zcoupling, config_flags%fire_lsm_zcoupling_ref, &
              config_flags%fire_wind_height, ioms, iome, joms, jome, iops, &
              iope, jops, jope, this%ua, this%va)

          call this%Destroy_u3d ()
          call this%Destroy_v3d ()
          call this%Destroy_phl ()
          call this%Destroy_z0 ()

        case (VINTERP_WINDS_FROM_10M_WINDS)
          call this%Get_u10 (datetime_now)
          call this%Get_v10 (datetime_now)

          this%ua = this%u10
          this%va = this%v10

        case default
          call Stop_simulation ('Error: wrong wind_vinterp_opt')

      end select

    end subroutine Update_atm_state

    function Wrfdata_t_const (file_name, config_flags) result (return_value)

      use, intrinsic :: iso_fortran_env, only : REAL32, INT32

      implicit none

      character (len = *), intent (in) :: file_name
      type (namelist_t), intent (in) :: config_flags
      type (wrfdata_t) :: return_value

      logical, parameter :: DEBUG_LOCAL = .false.
      real, parameter :: DEFAULT_Z0 = 0.1, DEFAULT_ZSF = 0.0, DEFAULT_DZDXF = 0.0, DEFAULT_DZDYF = 0.0, &
          DEFAULT_T2 = 123.4, DEFAULT_Q2 = 0.0, DEFAULT_PSFC = 0.0, DEFAULT_RAIN = 0.0

      real (kind = REAL32) :: att_real32
      integer (kind = INT32) :: att_int32


      if (DEBUG_LOCAL) Call Print_message ('Entering wrfdata_t constructor')

      return_value%file_name = trim (file_name)
      call Is_netcdf_file_present (trim (file_name))

        ! Init projection
      call Get_netcdf_att (trim (return_value%file_name), 'global', 'CEN_LAT', att_real32)
      return_value%cen_lat = att_real32

      call Get_netcdf_att (trim (return_value%file_name), 'global', 'CEN_LON', att_real32)
      return_value%cen_lon = att_real32

      call Get_netcdf_att (trim (return_value%file_name), 'global', 'TRUELAT1', att_real32)
      return_value%truelat1 = att_real32

      call Get_netcdf_att (trim (return_value%file_name), 'global', 'TRUELAT2', att_real32)
      return_value%truelat2 = att_real32

      call Get_netcdf_att (trim (return_value%file_name), 'global', 'STAND_LON', att_real32)
      return_value%stand_lon = att_real32

      call Get_netcdf_att (trim (return_value%file_name), 'global', 'DX', att_real32)
      return_value%dx = att_real32

      call Get_netcdf_att (trim (return_value%file_name), 'global', 'DY', att_real32)
      return_value%dy = att_real32

        ! latlon at mass points
      call return_value%Get_latlons ()

        ! latlon at corners
      call return_value%Get_latcloncs ()

        ! Init domain dimensions
      return_value%ids = 1
      call Get_netcdf_dim (trim (file_name), 'west_east_stag', att_int32)
      return_value%ide = att_int32
      return_value%jds = 1
      call Get_netcdf_dim (trim (file_name), 'south_north_stag', att_int32)
      return_value%jde = att_int32
      return_value%kds = 1
      call Get_netcdf_dim (trim (file_name), 'bottom_top_stag', att_int32)
      return_value%kde = att_int32

        ! Init rest of dimensions
      return_value%ims = return_value%ids
      return_value%ime = return_value%ide
      return_value%kms = return_value%kds
      return_value%kme = return_value%kde
      return_value%jms = return_value%jds
      return_value%jme = return_value%jde

      return_value%its = return_value%ids
      return_value%ite = return_value%ide
      return_value%kts = return_value%kds
      return_value%kte = return_value%kde
      return_value%jts = return_value%jds
      return_value%jte = return_value%jde

      if (DEBUG_LOCAL) call return_value%Print_domain()

        ! Init some vars to default values
      allocate (return_value%z0(return_value%ids:return_value%ide - 1, return_value%jds:return_value%jde - 1))
      return_value%z0 = DEFAULT_Z0

      allocate (return_value%rain(return_value%ids:return_value%ide - 1, return_value%jds:return_value%jde - 1))
      return_value%rain = DEFAULT_RAIN

      allocate (return_value%t2(return_value%ids:return_value%ide - 1, return_value%jds:return_value%jde - 1))
      return_value%t2 = DEFAULT_T2

      allocate (return_value%q2(return_value%ids:return_value%ide - 1, return_value%jds:return_value%jde - 1))
      return_value%q2 = DEFAULT_Q2

      allocate (return_value%psfc(return_value%ids:return_value%ide - 1, return_value%jds:return_value%jde - 1))
      return_value%psfc = DEFAULT_PSFC

      allocate (return_value%ua(return_value%ids:return_value%ide - 1, return_value%jds:return_value%jde - 1))
      return_value%ua = 0.0

      allocate (return_value%va(return_value%ids:return_value%ide - 1, return_value%jds:return_value%jde - 1))
      return_value%va = 0.0

      if (DEBUG_LOCAL) Call Print_message ('Leaving wrfdata_t constructor')

    end function Wrfdata_t_const

  end module wrfdata_mod
