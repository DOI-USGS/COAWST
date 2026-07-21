  module proj_lc_mod

    implicit none

    private

    real, parameter :: EARTH_RADIUS_M = 6370000.0 ! m
    real, parameter :: DEG_TO_RAD = 0.01745329251994329444
    real, parameter :: RAD_TO_DEG = 57.29577951308232522583

    public ::  proj_lc_t

    type :: proj_lc_t
      real :: known_lat, known_lon, dx, dy, pole_i, pole_j, known_i, known_j, hemi, &
          true_lat_1, true_lat_2,  standard_lon, cone_factor, r_earth
      integer :: nx, ny
    contains
      procedure, public :: Calc_latlon => Calc_lc_latlon_at_ij
      procedure, public :: Calc_ij => Calc_lc_ij_from_latlon
    end type proj_lc_t

    interface proj_lc_t
      procedure Proj_lc_t_const
    end interface proj_lc_t

  contains

    subroutine Calc_lc_cone (truelat1, truelat2, cone)

      implicit none

      real, intent(in)  :: truelat1, truelat2
      real, intent(out) :: cone


      if (abs (truelat1 - truelat2) > 0.1) then
         cone = alog10 (cos (truelat1 * DEG_TO_RAD)) - &
                alog10 (cos (truelat2 * DEG_TO_RAD))
         cone = cone / (alog10 (tan ((45.0 - abs (truelat1) / 2.0) * DEG_TO_RAD)) - &
                alog10 (tan ((45.0 - abs (truelat2) / 2.0) * DEG_TO_RAD)))
      else
         cone = sin (abs (truelat1) * DEG_TO_RAD )
      end if

    end subroutine Calc_lc_cone

    pure subroutine Calc_lc_ij_from_latlon (this, lat, lon, i, j)

      implicit none

      class (proj_lc_t), intent (in) :: this
      real, intent (in) :: lat, lon
      real, intent (out) :: i, j

      real :: arg, deltalon, tl1r, rm, ctl1r, rebydx


      deltalon = lon - this%standard_lon
      if (deltalon > 180.0) deltalon = deltalon - 360.0
      if (deltalon < -180.0) deltalon = deltalon + 360.0

      tl1r = this%true_lat_1 * DEG_TO_RAD
      ctl1r = cos (tl1r)

      rebydx = this%r_earth / this%dx

      rm = rebydx * ctl1r / this%cone_factor * &
          (tan ((90.0 * this%hemi - lat) * DEG_TO_RAD / 2.0) / &
          tan ((90.0 * this%hemi - this%true_lat_1) * DEG_TO_RAD / 2.0)) ** this%cone_factor

      arg = this%cone_factor * (deltalon * DEG_TO_RAD)

      i = this%pole_i + this%hemi * rm * sin (arg)
      i = this%hemi * i

      j = this%pole_j - rm * cos (arg)
      j = this%hemi * j

    end subroutine Calc_lc_ij_from_latlon

    subroutine Calc_lc_latlon_at_ij (this, i, j, lat, lon)

      implicit none

      class (proj_lc_t), intent (in) :: this
      real, intent(in) :: i, j
      real, intent(out) :: lat, lon

      real :: inew, jnew, r, chi, chi1, chi2, r2, xx, yy, rebydx


      chi1 = (90.0 - this%hemi * this%true_lat_1) * DEG_TO_RAD
      chi2 = (90.0 - this%hemi * this%true_lat_2) * DEG_TO_RAD

        ! Flip indeces in southern hemisphere
      inew = this%hemi * i
      jnew = this%hemi * j

        ! Compute radius**2 to i/j location
      rebydx = this%r_earth / this%dx
      xx = inew - this%pole_i
      yy = this%pole_j - jnew
      r2 = (xx * xx + yy * yy)
      r = sqrt (r2) / rebydx

        ! Convert to lat/lon
      if (r2 == 0.0) then
         lat = this%hemi * 90.0
         lon = this%standard_lon
      else
        lon = this%standard_lon + RAD_TO_DEG * ATAN2 (this%hemi * xx, yy) / this%cone_factor
        lon = amod (lon + 360.0, 360.0)
        if (chi1 == chi2) THEN
           chi = 2.0 * atan ((r / tan (chi1)) ** (1.0 / this%cone_factor) * tan (chi1 * 0.5))
        else
            chi = 2.0 * atan ((r * this%cone_factor / sin (chi1)) ** (1.0 / this%cone_factor) * &
                tan (chi1 * 0.5))
        end if
        lat = (90.0 - chi * RAD_TO_DEG) * this%hemi
      end if

      if (lon >  180.0) lon = lon - 360.0
      if (lon < -180.0) lon = lon + 360.0

    end subroutine Calc_lc_latlon_at_ij

    function Proj_lc_t_const (cen_lat, cen_lon, dx, dy, standard_lon, true_lat_1, &
        true_lat_2, nx, ny) result (return_value)

      implicit none

      real, intent (in) :: cen_lat, cen_lon, dx, dy, standard_lon, true_lat_1, true_lat_2
      integer, intent (in) :: nx, ny

      type (proj_lc_t) :: return_value

      real :: deltalon, tl1r, ctl1r, rsw, rebydx, arg


      return_value%r_earth = EARTH_RADIUS_M

      return_value%dx = dx
      return_value%dy = dy
      return_value%nx = nx
      return_value%ny = ny

      return_value%known_lat = cen_lat
      return_value%known_lon = cen_lon
      return_value%known_i = (nx + 1) / 2.0
      return_value%known_j = (ny + 1) / 2.0

      return_value%true_lat_1 = true_lat_1
      return_value%true_lat_2 = true_lat_2
      return_value%standard_lon = standard_lon

      if (true_lat_1 < 0.0) then
        return_value%hemi = -1.0
      else
        return_value%hemi = 1.0
      end if

      call Calc_lc_cone (return_value%true_lat_1, return_value%true_lat_2, return_value%cone_factor)

      deltalon = return_value%known_lon - return_value%standard_lon
      if (deltalon > + 180.0) deltalon = deltalon - 360.0
      if (deltalon < - 180.0) deltalon = deltalon + 360.0

      tl1r = return_value%true_lat_1 * DEG_TO_RAD
      ctl1r = cos (tl1r)

      rebydx = return_value%r_earth / return_value%dx
      rsw = rebydx * ctl1r / return_value%cone_factor * &
          (tan ((90.0 * return_value%hemi - return_value%known_lat) * DEG_TO_RAD / 2.0) / &
          tan ((90.0 * return_value%hemi - return_value%true_lat_1) * DEG_TO_RAD / 2.0)) ** return_value%cone_factor

      arg = return_value%cone_factor * (deltalon * DEG_TO_RAD)
      return_value%pole_i = return_value%hemi * return_value%known_i - return_value%hemi * rsw * sin (arg)
      return_value%pole_j = return_value%hemi * return_value%known_j + rsw * cos (arg)

    end function Proj_lc_t_const

  end module proj_lc_mod
