  module ignition_line_mod

    use namelist_mod, only: namelist_t, FIRE_MAX_IGNITIONS_IN_NAMELIST
    use stderrout_mod, only: Stop_simulation

    implicit none

    private

    public :: ignition_line_t

    integer :: fire_num_ignitions

    type :: ignition_line_t
      real, dimension(:), allocatable :: start_x, start_y, end_x, end_y, start_time, end_time, radius, ros
    contains
      procedure, public :: Init => Initialization
      procedure, public :: Ignite_fire => Ignite_fire
    end type ignition_line_t

  contains

    subroutine Ignite_fire (this, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, &
        n_line, start_ts, end_ts, coord_xf, coord_yf, unit_xf, unit_yf, lfn, tign, ignited)

    ! Ignite a circular/line fire
    ! ignite fire in the region within radius r from the line (sx, sy) to (ex,e y).

      implicit none

      class (ignition_line_t), intent (in) :: this
      integer, intent (in) :: ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, n_line
      real, intent (in) :: unit_xf, unit_yf
        ! the time step start and end
      real, intent (in) :: start_ts, end_ts
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: coord_xf, coord_yf
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: lfn, tign
      integer, intent(out) :: ignited

      integer :: i, j
      real :: lfn_new, time_ign, ax, ay, rele, d, sx, sy, ex, ey, st, et, cx2, cy2, dmax, dmin, &
           end_x, end_y, radius, start_time, end_time, ros, tos


      sx = this%start_x(n_line)
      sy = this%start_y(n_line)
      end_x = this%end_x(n_line)
      end_y = this%end_y(n_line)
      start_time = this%start_time(n_line)
      end_time = this%end_time(n_line)
      radius = this%radius(n_line)
      ros = this%ros(n_line)

      tos = radius / ros          ! time of spread to the given radius
      st = start_time             ! the start time of ignition considered in this time step
      et = min (end_ts, end_time) ! the end time of the ignition segment in this time step

        ! (start_ts, end_ts) must be subset (start_time, end_time + tos)
      if (start_ts > et + tos .or. end_ts < st) return

      if (start_time < end_time) then  ! we really want to test start_time .ne. end_time, but avoiding test of floats on equality
        rele =  (et - start_time) / (end_time - start_time)    ! relative position of et in the segment (start,end)
        ex = sx + rele * (end_x - sx)
        ey = sy + rele * (end_y - sy)
      else
        ex = end_x
        ey = end_y
      end if

      cx2 = unit_xf * unit_xf
      cy2 = unit_yf * unit_yf

      ignited = 0
      dmax = 0.0
      dmin = huge (dmax)
      do j = jfts, jfte
        do i = ifts, ifte
          ax = coord_xf(i, j)
          ay = coord_yf(i, j)
            ! get d= distance from the nearest point on the ignition segment
            ! and time_ign = the ignition time there
          call Nearest (d, time_ign, ax, ay, sx, sy, st, ex, ey, et, cx2, cy2)
          dmax = max (d, dmax)
          dmin = min (d, dmin)
            ! lft at end_ts
          lfn_new = d - min (radius, ros * (end_ts - time_ign))

          if (.not. lfn_new > 0.0) ignited = ignited + 1

          if (lfn(i, j) > 0.0 .and. .not. lfn_new > 0.0) then
              ! newly ignited now
            tign(i, j) = time_ign + d / ros
            tign(i, j) = min (max (tign(i, j), start_ts), end_ts)
          end if
          lfn(i, j) = min (lfn(i, j), lfn_new)
        end do
      end do

      return

    contains

      subroutine Nearest (d, t, ax, ay, sx, sy, st, ex, ey, et, cx2, cy2)

      ! d the distance of a and the nearest point z on the segment [x,y]
      ! t linear interpolation from the values st, et to the point z
      !
      ! Compute d as the distance (ax, ay) from the midpoint (mx, my) of the
      ! line from (sx, xy) to (ex, ey) minus a correction: |a-c|^2 = |a-m|^2 - |m-c|^2
      ! when |m-c| >= |s-e|/2 the nearest point is one of the endpoints
      ! the computation work also for the case when s=e exactly or approximately
      !
      !           a
      !          /| \
      !     s---m-c--e
      !
      ! |m-c| = |a-m| cos (a-m,e-s)
      !       = |a-m| (a-m).(e-s))/(|a-m|*|e-s|)

        implicit none

        real, intent (in) :: ax, ay, sx, sy, st, ex, ey, et, cx2, cy2
        real, intent (out) :: d, t

        real :: mx, my, dam2, dames, am_es, cos2, dmc2, mcrel, mid_t, dif_t, des2, cx, cy


          ! midpoint m = (mx,my)
        mx = (sx + ex) * 0.5
        my = (sy + ey) * 0.5
          ! |a-m|^2
        dam2 = (ax - mx) * (ax - mx) * cx2 + (ay - my) * (ay - my) * cy2
          ! des2 = |e-s|^2
        des2 = (ex - sx) * (ex - sx) * cx2 + (ey - sy) * (ey - sy) * cy2
        dames = dam2 * des2
          ! am_es = (a-m).(e-s)
        am_es = (ax - mx) * (ex - sx) * cx2 + (ay - my) * (ey - sy) * cy2
        if (dames > 0.0) then
            ! cos2 = cos^2 (a-m,e-s)
          cos2 = (am_es * am_es) / dames
        else ! point a already is the midpoint
          cos2 = 0.0
        end if
          ! dmc2 = |m-c|^2
        dmc2 = dam2 * cos2
        if (4.0 * dmc2 < des2) then
            ! if |m-c|<=|e-s|/2
            ! d = sqrt(max(dam2 - dmc2,0.))     ! d=|a-m|^2 - |m-c|^2, guard rounding
            ! relative distance of c from m
          mcrel = sign (sqrt (4.0 * dmc2 / des2), am_es)
        else if (am_es > 0) then
            ! if cos > 0, closest is e
          mcrel = 1.0
        else
           ! closest is s
          mcrel = -1.0
        end if
          ! interpolate to c by going from m
        cx = (ex + sx) * 0.5 + mcrel * (ex - sx) * 0.5
        cy = (ey + sy) * 0.5 + mcrel * (ey - sy) * 0.5

          ! output:
          ! 1) |a-c|^2
        d = sqrt ((ax - cx) * (ax - cx) * cx2 + (ay - cy) * (ay - cy) * cy2)
          ! 2) interpolate to c by going from m
        t = (et + st) * 0.5 + mcrel * (et - st) * 0.5

      end subroutine Nearest

    end subroutine Ignite_fire

    subroutine Initialization (this, config_flags)

      implicit none

      class (ignition_line_t), intent (out) :: this
      type (namelist_t), intent (in) :: config_flags

      integer :: i, n_ignitions


      n_ignitions = config_flags%fire_num_ignitions
      if (n_ignitions <= 0) call Stop_simulation ('Not enough ignitions set')
      if (n_ignitions > FIRE_MAX_IGNITIONS_IN_NAMELIST) call Stop_simulation ('FIRE_MAX_IGNITIONS_IN_NAMELIST too small')

      allocate (this%start_x(n_ignitions))
      allocate (this%start_y(n_ignitions))
      allocate (this%end_x(n_ignitions))
      allocate (this%end_y(n_ignitions))
      allocate (this%start_time(n_ignitions))
      allocate (this%end_time(n_ignitions))
      allocate (this%radius(n_ignitions))
      allocate (this%ros(n_ignitions))

      if (n_ignitions >= 1) then
        this%start_x(1) = config_flags%fire_ignition_start_lon1
        this%start_y(1) = config_flags%fire_ignition_start_lat1
        this%end_x(1) = config_flags%fire_ignition_end_lon1
        this%end_y(1) = config_flags%fire_ignition_end_lat1
        this%ros(1) = config_flags%fire_ignition_ros1
        this%radius(1) = config_flags%fire_ignition_radius1
        this%start_time(1) = config_flags%fire_ignition_start_time1
        this%end_time(1) = config_flags%fire_ignition_end_time1
      end if

      if (n_ignitions >= 2) then
        this%start_x(2) = config_flags%fire_ignition_start_lon2
        this%start_y(2) = config_flags%fire_ignition_start_lat2
        this%end_x(2) = config_flags%fire_ignition_end_lon2
        this%end_y(2) = config_flags%fire_ignition_end_lat2
        this%ros(2) = config_flags%fire_ignition_ros2
        this%radius(2) = config_flags%fire_ignition_radius2
        this%start_time(2) = config_flags%fire_ignition_start_time2
        this%end_time(2) = config_flags%fire_ignition_end_time2
      end if

      if (n_ignitions >= 3) then
        this%start_x(3) = config_flags%fire_ignition_start_lon3
        this%start_y(3) = config_flags%fire_ignition_start_lat3
        this%end_x(3) = config_flags%fire_ignition_end_lon3
        this%end_y(3) = config_flags%fire_ignition_end_lat3
        this%ros(3) = config_flags%fire_ignition_ros3
        this%radius(3) = config_flags%fire_ignition_radius3
        this%start_time(3) = config_flags%fire_ignition_start_time3
        this%end_time(3) = config_flags%fire_ignition_end_time3
      end if

      if (n_ignitions >= 4) then
        this%start_x(4) = config_flags%fire_ignition_start_lon4
        this%start_y(4) = config_flags%fire_ignition_start_lat4
        this%end_x(4) = config_flags%fire_ignition_end_lon4
        this%end_y(4) = config_flags%fire_ignition_end_lat4
        this%ros(4) = config_flags%fire_ignition_ros4
        this%radius(4) = config_flags%fire_ignition_radius4
        this%start_time(4) = config_flags%fire_ignition_start_time4
        this%end_time(4) = config_flags%fire_ignition_end_time4
      end if

      if (n_ignitions >= 5) then
        this%start_x(5) = config_flags%fire_ignition_start_lon5
        this%start_y(5) = config_flags%fire_ignition_start_lat5
        this%end_x(5) = config_flags%fire_ignition_end_lon5
        this%end_y(5) = config_flags%fire_ignition_end_lat5
        this%ros(5) = config_flags%fire_ignition_ros5
        this%radius(5) = config_flags%fire_ignition_radius5
        this%start_time(5) = config_flags%fire_ignition_start_time5
        this%end_time(5) = config_flags%fire_ignition_end_time5
      end if

      do i = 1, n_ignitions
        if (this%radius(i) <= 0.0) call Stop_simulation ('Radius ignition line must be > 0')
          ! Expand ignition data given as zero
        if (this%end_x(i) == 0.0) this%end_x(i) = this%start_x(i)
        if (this%end_y(i) == 0.0) this%end_y(i) = this%start_y(i)
        if (this%end_time(i) == 0.0) this%end_time(i) = this%start_time(i)
      end do

    end subroutine Initialization

  end module ignition_line_mod
