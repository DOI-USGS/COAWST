  module geogrid_mod

    use netcdf_mod, only : Get_netcdf_var, Get_netcdf_att, Get_netcdf_dim, Is_netcdf_file_present, Is_netcdf_var_present
    use proj_lc_mod, only : proj_lc_t
    use stderrout_mod, only: Print_message


    implicit none

    private

    public :: geogrid_t

    character (len = *), parameter :: VAR_NAME_LFNINI = 'lfn_init'

    type :: geogrid_t
      real, dimension (:, :), allocatable :: elevations, dz_dxs, dz_dys, fuel_cats, xlat, xlong, xlat_c, xlong_c, lfn_init
      real :: dx = 0.0, dy = 0.0, cen_lat = 0.0, cen_lon = 0.0, true_lat_1 = 0.0, true_lat_2 = 0.0, stand_lon = 0.0
      integer :: ids = 1, jds = 1, ide, jde = 0, ifds = 1, jfds = 1, ifde = 0, jfde = 0, sr_x = 0, sr_y = 0, map_proj = 0
    contains
      procedure, public :: Print => Print_geogrid
      procedure, public :: Get_atm_proj => Get_atm_proj
    end type geogrid_t

    interface geogrid_t
      module procedure Geogrid_t_const
    end interface geogrid_t

  contains

    function Geogrid_t_const (file_name) result (return_value)

      use, intrinsic :: iso_fortran_env, only : REAL32, INT32

      implicit none

      character (len = *), intent (in) :: file_name
      type (geogrid_t) :: return_value


      real (kind = REAL32), dimension(:, :, :), allocatable :: var_real32
      real (kind = REAL32) :: att_real32
      integer (kind = INT32) :: att_int32


      call Is_netcdf_file_present (file_name = trim (file_name))

      call Get_netcdf_var (trim (file_name), 'ZSF', var_real32)
      return_value%elevations = var_real32(:, :, 1)
      deallocate (var_real32)

      call Get_netcdf_var (trim (file_name), 'DZDXF', var_real32)
      return_value%dz_dxs = var_real32(:, :, 1)
      deallocate (var_real32)

      call Get_netcdf_var (trim (file_name), 'DZDYF', var_real32)
      return_value%dz_dys = var_real32(:, :, 1)
      deallocate (var_real32)

      call Get_netcdf_var (trim (file_name), 'NFUEL_CAT', var_real32)
      return_value%fuel_cats = var_real32(:, :, 1)
      deallocate (var_real32)

      call Get_netcdf_dim (trim (file_name), 'south_north_stag', att_int32)
      return_value%jde = att_int32

      call Get_netcdf_dim (trim (file_name), 'west_east_stag', att_int32)
      return_value%ide = att_int32

      call Get_netcdf_dim (trim (file_name), 'south_north_subgrid', att_int32)
      return_value%jfde = att_int32

      call Get_netcdf_dim (trim (file_name), 'west_east_subgrid', att_int32)
      return_value%ifde = att_int32

      call Get_netcdf_att (trim (file_name), 'global', 'DX', att_real32)
      return_value%dx = att_real32

      call Get_netcdf_att (trim (file_name), 'global', 'DY', att_real32)
      return_value%dy = att_real32

      call Get_netcdf_att (trim (file_name), 'global', 'CEN_LAT', att_real32)
      return_value%cen_lat = att_real32

      call Get_netcdf_att (trim (file_name), 'global', 'CEN_LON', att_real32)
      return_value%cen_lon = att_real32

      call Get_netcdf_att (trim (file_name), 'global', 'MAP_PROJ', att_int32)
      return_value%map_proj = att_int32

      call Get_netcdf_att (trim (file_name), 'global', 'TRUELAT1', att_real32)
      return_value%true_lat_1 = att_real32

      call Get_netcdf_att (trim (file_name), 'global', 'TRUELAT2', att_real32)
      return_value%true_lat_2 = att_real32

      call Get_netcdf_att (trim (file_name), 'global', 'STAND_LON', att_real32)
      return_value%stand_lon = att_real32

      call Get_netcdf_var (trim (file_name), 'XLAT_M', var_real32)
      return_value%xlat = var_real32(:, :, 1)
      deallocate (var_real32)

      call Get_netcdf_var (trim (file_name), 'XLONG_M', var_real32)
      return_value%xlong = var_real32(:, :, 1)
      deallocate (var_real32)

      call Get_netcdf_var (trim (file_name), 'XLAT_C', var_real32)
      return_value%xlat_c = var_real32(:, :, 1)
      deallocate (var_real32)

      call Get_netcdf_var (trim (file_name), 'XLONG_C', var_real32)
      return_value%xlong_c = var_real32(:, :, 1)
      deallocate (var_real32)

      call Get_netcdf_att (trim (file_name), 'global', 'sr_x', att_int32)
      return_value%sr_x = att_int32

      call Get_netcdf_att (trim (file_name), 'global', 'sr_y', att_int32)
      return_value%sr_y = att_int32

      if (Is_netcdf_var_present (trim (file_name), VAR_NAME_LFNINI)) &
          call Get_netcdf_var (trim (file_name), VAR_NAME_LFNINI, return_value%lfn_init)

    end function Geogrid_t_const

    function Get_atm_proj (this) result (return_value)

      use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

      implicit none

      class (geogrid_t), intent (in) :: this
      type (proj_lc_t) :: return_value

      logical, parameter :: DEBUG_LOCAL = .false.


      if (DEBUG_LOCAL) call Print_message ('Entering Get_atm_proj')

      if (DEBUG_LOCAL) then
        write (OUTPUT_UNIT, *) 'cen_lat = ', this%cen_lat
        write (OUTPUT_UNIT, *) 'cen_lon = ', this%cen_lon
        write (OUTPUT_UNIT, *) 'dx = ', this%dx
        write (OUTPUT_UNIT, *) 'dy = ', this%dy
        write (OUTPUT_UNIT, *) 'stand_lon = ', this%stand_lon
        write (OUTPUT_UNIT, *) 'true_lat_1 = ', this%true_lat_1
        write (OUTPUT_UNIT, *) 'true_lat_2 = ', this%true_lat_2
        write (OUTPUT_UNIT, *) 'nx = ', this%ide - 1
        write (OUTPUT_UNIT, *) 'ny = ', this%jde - 1
      end if

      return_value = proj_lc_t (cen_lat = this%cen_lat , cen_lon = this%cen_lon, &
          dx = this%dx, dy = this%dy, standard_lon = this%stand_lon, true_lat_1 = this%true_lat_1, &
          true_lat_2 = this%true_lat_2, nx = this%ide - 1, ny = this%jde - 1)

      if (DEBUG_LOCAL) call Print_message ('Leaving Get_atm_proj')

    end function Get_atm_proj

    subroutine Print_geogrid (this)

      use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

      implicit none

      class (geogrid_t), intent (in) :: this


      write (OUTPUT_UNIT, *) 'Contents of geogrid_t object:'
      if (allocated (this%elevations)) &
          write (OUTPUT_UNIT, *) 'Shape (elevations) = ', shape (this%elevations)
      if (allocated (this%dz_dxs)) &
          write (OUTPUT_UNIT, *) 'Shape (dz_dxs) = ', shape (this%dz_dxs)
      if (allocated (this%dz_dys)) &
          write (OUTPUT_UNIT, *) 'Shape (dz_dys) = ', shape (this%dz_dys)
      if (allocated (this%fuel_cats)) &
          write (OUTPUT_UNIT, *) 'Shape (fuel_cats) = ', shape (this%fuel_cats)
      if (allocated (this%xlat)) &
          write (OUTPUT_UNIT, *) 'Shape (xlat) = ', shape (this%xlat)
      if (allocated (this%xlong)) &
          write (OUTPUT_UNIT, *) 'Shape (xlong) = ', shape (this%xlong)

      write (OUTPUT_UNIT, *) 'ids = ', this%ids
      write (OUTPUT_UNIT, *) 'ide = ', this%ide
      write (OUTPUT_UNIT, *) 'jds = ', this%jds
      write (OUTPUT_UNIT, *) 'jde = ', this%jde

      write (OUTPUT_UNIT, *) 'dx = ', this%dx
      write (OUTPUT_UNIT, *) 'dy = ', this%dy
      write (OUTPUT_UNIT, *) 'sr_x = ', this%sr_x
      write (OUTPUT_UNIT, *) 'sr_y = ', this%sr_y

      write (OUTPUT_UNIT, *) 'map_proj = ', this%map_proj
      write (OUTPUT_UNIT, *) 'cen_lat = ', this%cen_lat
      write (OUTPUT_UNIT, *) 'cen_lon = ', this%cen_lon
      write (OUTPUT_UNIT, *) 'truelat1 = ', this%true_lat_1
      write (OUTPUT_UNIT, *) 'truelat2 = ', this%true_lat_2
      write (OUTPUT_UNIT, *) 'stand_lon = ', this%stand_lon
      write (OUTPUT_UNIT, *)

    end subroutine Print_geogrid

  end module geogrid_mod
