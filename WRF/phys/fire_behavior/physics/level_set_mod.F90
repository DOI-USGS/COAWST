  module level_set_mod

  ! References:
  !
  ! Based on S. Osher and R. Fedkiw, Level set methods and dynamic implicit surfaces,
  ! Springer, 2003, Sec. 6.4, as implemented in toolboxLS for Matlab by
  ! I. Mitchell, A toolbox of Level Set Methods (Version 1.1), TR-2007-11,
  ! Dept. Computer Science, University of British Columbia, 2007
  ! http://www.cs.ubc.ca/\~mitchell/Toolbo\LS
  !
  ! D. Munoz-Esparza, B. Kosovic, P. Jimenez, J. Coen: "An accurate
  ! fire-spread algorithm in the Weather Research and Forecasting model using the
  ! level-set method", Journal of Advances in Modeling Earth Systems, 2018
  ! https://doi.org/10.1002/2017MS001108

    use ros_wrffire_mod, only: ros_wrffire_t
    use stderrout_mod, only: Stop_simulation, Print_message
    use state_mod, only: state_fire_t, N_POINTS_IN_HALO
    use ignition_line_mod, only : ignition_line_t
    use ros_mod, only : ros_t
    use constants_mod, only : PI

#ifdef DM_PARALLEL
    use mpi_mod, only : Do_halo_exchange, Sum_across_mpi_tasks, Max_across_mpi_tasks, Min_across_mpi_tasks
#endif

    implicit none

    private

    public :: Calc_fuel_left, Update_ignition_times, Reinit_level_set, Prop_level_set, Extrapol_var_at_bdys, Stop_if_close_to_bdy, &
        Copy_lfnout_to_lfn, Reinit_level_set_fast_dist, Check_isolated_negative_lfn

    integer, parameter :: BDY_ENO1 = 10, FAST_DIST_REINIT_FSM = 1
    integer, parameter :: FAR = 0, TRIAL = 1, KNOWN = 2

  contains

    subroutine Calc_fuel_left (ims, ime, jms, jme, its, ite, jts, jte, ifs, ife, jfs, jfe, &
        lfn, tign, fuel_time, time_now, fuel_frac, fire_area, fuel_frac_burnt_dt)

      implicit none

      !*** purpose: determine fraction of fuel remaining
      !*** NOTE: because variables are cell centered, need halo/sync width 1 before

      integer, intent (in) :: its, ite, jts, jte, ims, ime, jms, jme, ifs, ife, jfs, jfe
      real, intent (in), dimension (ims:ime, jms:jme) :: lfn,tign, fuel_time
      real, intent (in) :: time_now
      real, intent (in out), dimension (ims:ime, jms:jme) :: fuel_frac
      real, intent (out), dimension (ims:ime, jms:jme) :: fire_area, fuel_frac_burnt_dt

      real, dimension (ifs:ife, jfs:jfe) :: fuel_frac_end
      integer :: i, j, ir, jr, icl, jcl, isubcl, jsubcl, i2, j2, ii, jj
      real :: fmax, frat, helpsum1, helpsum2, fuel_left_ff, fire_area_ff, rx, ry, tignf(2,2)
         ! help variables instead of arrays fuel_left_f and fire_area_f 
      real :: lffij, lffi1j, lffij1, lffi1j1, tifij, tifi1j, tifij1, tifi1j1, tx, ty, txx, tyy
         ! variables for calculation instead of lff(i,j) and tif(i,j)is lffij,tifij etc..
      character (len = 128) :: msg
      integer :: m, omp_get_thread_num
      integer, parameter :: fuel_left_irl = 2       ! "submesh to compute fuel lwft, even, at least 2" ""
      integer, parameter :: fuel_left_jrl = 2       ! "submesh to compute fuel lwft, even, at least 2" ""

        ! refinement
      ir = fuel_left_irl
      jr = fuel_left_jrl

      if ((ir /= 2) .or. (jr /= 2)) call Stop_simulation ('fuel_left: ir.ne.2 or jr.ne.2 ')

      rx = 1.0 / ir 
      ry = 1.0 / jr

        ! example for ir=2:
        !
        !                     |      coarse cell      |
        !      its-1                     its                                   ite             ite+1
        ! -------X------------|-----.-----X-----.-----|--........----|----------X----------|---------X
        !           fine node 1           2           3                   2*(ite-its+1) 
        !                fine cell  1           2                             cell 2*(ite-its+1)

        !  Loop over cells in Tile 
        !  Changes made by Volodymyr Kondratenko 09/24/2009
      do icl = its, ite
        do jcl = jts, jte
          helpsum1 = 0.0
          helpsum2 = 0.0
            ! Loop over subcells in cell #(icl,jcl)
          do isubcl = 1, ir
            do jsubcl = 1, jr 
              i = (icl - its) * ir + isubcl
              j = (jcl - jts) * jr + jsubcl
                ! Direct calculation tif and lff, avoiding arrays, just for case ir=jr=2
              if ((isubcl == 1) .and. (jsubcl == 1)) then
                 i2 = icl - 1
                 j2 = jcl - 1
                 ty = 0.5
                 tx = 0.5
                 tyy = 1.0
                 txx = 1.0
              else if ((isubcl == 2).and. (jsubcl == 1)) then
                i2 = icl
                j2 = jcl - 1
                ty = 0.5
                tx = 0
                tyy = 1.0
                txx = 0.5
              else if ((isubcl == 1) .and. (jsubcl == 2)) then
                i2 = icl - 1
                j2 = jcl
                tx = 0.5
                ty = 0
                txx = 1.0
                tyy = 0.5
              else if ((isubcl == 2) .and. (jsubcl == 2)) then
                i2 = icl
                j2 = jcl
                tx = 0
                ty = 0
                txx = 0.5
                tyy = 0.5
              else
                call Stop_simulation ('fuel_left: isubcl,jsubcl should be only 1 or 2')
              endif 

              ! calculation lff and tif in 4 endpoints of subcell
              lffij=                                &    
                  (1 - tx) * (1 - ty) * lfn(i2, j2) &
                  + (1 - tx) * ty * lfn(i2, j2 + 1) &
                  + tx * (1 - ty) * lfn(i2 + 1, j2) &
                  + tx * ty * lfn(i2 + 1, j2 + 1)
              lffi1j=                                &
                  (1 - txx) * (1 - ty) * lfn(i2, j2) &
                  + (1 - txx) * ty * lfn(i2 ,j2 + 1) &
                  + (txx) * (1 - ty) * lfn(i2 + 1, j2) &
                  + (txx) * ty * lfn(i2 + 1, j2 + 1)
              lffij1=                            &
                  (1 - tx) * (1 - tyy) * lfn(i2, j2) &
                  + (1 - tx) * tyy * lfn(i2, j2 + 1) &
                  + tx * (1 - tyy) * lfn(i2 + 1, j2) &
                  + tx * tyy * lfn(i2 + 1, j2 + 1)
              lffi1j1 =                               &
                  (1 - txx) * (1 - tyy) * lfn(i2, j2) &
                  + (1 - txx) * tyy * lfn(i2, j2 + 1) &        
                  + txx * (1 - tyy) * lfn(i2 + 1, j2) &
                  + txx * tyy * lfn(i2 + 1, j2 + 1)
 
               ! get ready to fix up tign values to be interpolated
              do ii = 1, 2
                do jj = 1, 2
                  tignf(ii, jj) = tign(i2 + ii - 1, j2 + jj - 1)
                end do
              end do
              tifij =                                 &
                  (1 - tx) * (1 - ty) * tignf(1, 1) &
                  + (1 - tx) * ty * tignf(1, 1 + 1) &
                  + tx * (1 - ty) * tignf(1 + 1, 1) &
                  + tx * ty * tignf(1 + 1, 1 + 1)
              tifi1j =                               &
                  (1 - txx) * (1 - ty) * tignf(1, 1) &
                  + (1 - txx) * ty * tignf(1, 1 + 1) &
                  + txx * (1 - ty) * tignf(1 + 1, 1) &
                  + txx * ty * tignf(1 + 1, 1 + 1)
              tifij1 =                               &
                  (1 - tx) * (1 - tyy) * tignf(1, 1) &
                  + (1 - tx) * tyy * tignf(1, 1 + 1) &
                  + tx * (1 - tyy) * tignf(1+1, 1) &
                  + tx * tyy * tignf(1 + 1, 1 + 1)
             tifi1j1 =                               &
                  (1 - txx) * (1 - tyy) * tignf(1, 1) &
                  + (1 - txx) * tyy * tignf(1, 1 + 1) &
                  + txx * (1 - tyy) * tignf(1 + 1, 1) &
                  + txx * tyy * tignf(1 + 1, 1 + 1)
 
             call Calc_fuel_left_at_grid_point (lffij, lffij1, lffi1j, lffi1j1, &
                 tifij, tifij1, tifi1j, tifi1j1, time_now, fuel_time(icl,jcl), &
                 fuel_left_ff, fire_area_ff)

                ! consistency check
              if (fire_area_ff < -1e-6 .or.  &
                  (fire_area_ff == 0.0 .and. fuel_left_ff < 1.0 - 1e-6)) then
                !$OMP CRITICAL(FIRE_CORE_CRIT)
                write (msg, '(a, 2i6, 2(a, f11.8))') 'fuel_left: at node', i, j, &
                    ' of refined mesh fuel burnt', 1 - fuel_left_ff, ' fire area', fire_area_ff
                !$OMP END CRITICAL(FIRE_CORE_CRIT)
                call Stop_simulation (msg)
              endif

              helpsum1 = helpsum1 + fuel_left_ff
              helpsum2 = helpsum2 + fire_area_ff
            end do
          end do
          fuel_frac_end(icl, jcl) = helpsum1 / (ir * jr)
          fire_area(icl, jcl) = helpsum2 / (ir * jr)
          fuel_frac_burnt_dt(icl, jcl) = fuel_frac(icl, jcl) - fuel_frac_end(icl, jcl)
          fuel_frac(icl, jcl) = fuel_frac_end(icl, jcl)
        end do 
      end do

        ! consistency check after sum
      fmax = 0.0
      do j = jts, jte
          do i = its, ite        
             if (fire_area(i, j) == 0.0) then
               if (fuel_frac_end(i, j) < 1.-1e-6) then
                 !$OMP CRITICAL(FIRE_CORE_CRIT)
                 write (msg, '(a, 2i6, 2(a, f11.8))') 'fuel_left: at node', i, j, &
                     ' fuel burnt', 1 - fuel_frac_end(i, j), ' but fire area', fire_area(i, j)
                 !$OMP END CRITICAL(FIRE_CORE_CRIT)
                     call Stop_simulation (msg)
               end if
             else
               frat = (1 - fuel_frac_end(i, j)) / fire_area(i, j)
               fmax = max (fmax, frat)
             end if
          end do
      end do

      return
    end subroutine Calc_fuel_left

    subroutine Calc_fuel_left_at_grid_point (lfn00, lfn01, lfn10, &
        lfn11, tign00, tign01, tign10, tign11, time_now, fuel_time_cell, &
        fuel_frac_left, fire_frac_area)

      ! Example:
      !
      !        lfn<0         lfn>0
      !      (0,1) -----O--(1,1)      O = points on the fireline
      !            |      \ |         A = the burning area for computing
      !            |       \|             fuel_frac_left(i, j)
      !            |   A    O 
      !            |        |
      !            |        |
      !       (0,0)---------(1,0)
      !       lfn<0          lfn<0
      !

      implicit none

      real, intent (in) :: lfn00, lfn01, lfn10, lfn11, tign00, tign01, tign10, tign11, &
          time_now, fuel_time_cell
      real, intent (out) :: fuel_frac_left, fire_frac_area

      real :: ps, aps, ta, t00, t01, t10, t11


        ! minus time since ignition, 0 if no ignition yet
        ! it is possible to have 0 in fire region when ignition time falls in 
        ! inside the time step because lfn is updated at the beginning of the time step
      t00 = tign00 - time_now
      if (lfn00 > 0.0 .or. t00 > 0.0) t00 = 0.0
      t01 = tign01 - time_now
      if (lfn01 > 0.0 .or. t01 > 0.0) t01 = 0.0
      t10 = tign10 - time_now
      if (lfn10 > 0.0 .or. t10 > 0.0) t10 = 0.0
      t11 = tign11 - time_now
      if (lfn11 > 0.0 .or. t11 > 0.0) t11 = 0.0
        ! average negative time since ignition
      ta = 0.25 * (t00 + t01 + t10 + t11)

        ! Approx burning area
      ps = lfn00 + lfn01 + lfn10 + lfn11
      aps = abs (lfn00) + abs (lfn01) + abs (lfn10) + abs (lfn11)
      aps = max (aps, tiny (aps))
      fire_frac_area = (-ps / aps + 1.0) / 2.0
      fire_frac_area = max (fire_frac_area, 0.0)
      fire_frac_area = min (fire_frac_area, 1.0)
    
        ! Calc remaining fuel fraction
      fuel_frac_left = 1.0
      if (fire_frac_area > 0.0) fuel_frac_left = fire_frac_area * exp (ta / fuel_time_cell) + (1.0 - fire_frac_area)
      if (fuel_frac_left > 1.0) call Stop_simulation ('Remaining fuel fraction > 1')

    end subroutine Calc_fuel_left_at_grid_point

    pure subroutine Copy_lfnout_to_lfn (ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, lfn_out, lfn)

      implicit none

      integer, intent (in) :: ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: lfn_out
      real, dimension (ifms:ifme, jfms:jfme), intent (out) :: lfn

      integer :: i, j


      do j = jfts, jfte
        do i = ifts, ifte
          lfn(i, j) = lfn_out(i, j)
        end do
      end do

    end subroutine Copy_lfnout_to_lfn

    subroutine Check_isolated_negative_lfn (grid, threshold, min_size, max_size)

      implicit none

      type (state_fire_t), intent (in out) :: grid
      real, intent (in), optional :: threshold
      integer, intent (in), optional :: min_size, max_size

      character (len = 256) :: msg
      integer :: cluster_i, cluster_j, cluster_size, flagged_size
      integer :: head, tail, i, ic, inb, j, jc, jnb, max_queue, min_size_value, max_size_value
      integer :: nx_loc, ny_loc
      logical :: found_cluster, touches_boundary
      logical, dimension(:, :), allocatable :: visited
      integer, dimension(:), allocatable :: queue_i, queue_j
      real :: detection_flag, detection_sum, global_location_i, global_location_j, global_size_value
      real :: location_i, location_j, size_value, threshold_value


      threshold_value = 0.0
      if (present (threshold)) threshold_value = threshold

      min_size_value = 1
      if (present (min_size)) min_size_value = min_size

      max_size_value = 6
      if (present (max_size)) max_size_value = max_size

      if (max_size_value < min_size_value) return

      nx_loc = grid%ifpe - grid%ifps + 1
      ny_loc = grid%jfpe - grid%jfps + 1
      if (nx_loc <= 0 .or. ny_loc <= 0) return
      max_queue = nx_loc * ny_loc

      allocate (visited(grid%ifps:grid%ifpe, grid%jfps:grid%jfpe))
      allocate (queue_i(max_queue))
      allocate (queue_j(max_queue))

      visited = .false.
      found_cluster = .false.
      flagged_size = 0
      cluster_i = grid%ifps
      cluster_j = grid%jfps

      Cluster_search: do j = grid%jfps, grid%jfpe
        do i = grid%ifps, grid%ifpe
          if (visited(i, j)) cycle
          visited(i, j) = .true.
          if (grid%lfn(i, j) >= threshold_value) cycle

          head = 1
          tail = 1
          queue_i(1) = i
          queue_j(1) = j
          cluster_size = 0
          touches_boundary = .false.

          do while (head <= tail)
            ic = queue_i(head)
            jc = queue_j(head)
            head = head + 1

            cluster_size = cluster_size + 1
            if (ic == grid%ifps .or. ic == grid%ifpe .or. jc == grid%jfps .or. jc == grid%jfpe) touches_boundary = .true.

            do jnb = jc - 1, jc + 1
              do inb = ic - 1, ic + 1
                if (inb == ic .and. jnb == jc) cycle
                if (inb < grid%ifps .or. inb > grid%ifpe) cycle
                if (jnb < grid%jfps .or. jnb > grid%jfpe) cycle
                if (visited(inb, jnb)) cycle
                visited(inb, jnb) = .true.
                if (grid%lfn(inb, jnb) >= threshold_value) cycle

                tail = tail + 1
                queue_i(tail) = inb
                queue_j(tail) = jnb
              end do
            end do
          end do

          if (cluster_size >= min_size_value .and. cluster_size <= max_size_value .and. .not. touches_boundary) then
            found_cluster = .true.
            flagged_size = cluster_size
            cluster_i = i
            cluster_j = j
            exit Cluster_search
          end if
        end do
      end do Cluster_search

      deallocate (visited)
      deallocate (queue_i)
      deallocate (queue_j)

      detection_flag = 0.0
      location_i = -1.0
      location_j = -1.0
      size_value = -1.0
      if (found_cluster) then
        detection_flag = 1.0
        location_i = real (cluster_i)
        location_j = real (cluster_j)
        size_value = real (flagged_size)
      end if

#ifdef DM_PARALLEL
      call Sum_across_mpi_tasks (detection_flag, grid%cart_comm, detection_sum)
      call Max_across_mpi_tasks (location_i, grid%cart_comm, global_location_i)
      call Max_across_mpi_tasks (location_j, grid%cart_comm, global_location_j)
      call Max_across_mpi_tasks (size_value, grid%cart_comm, global_size_value)
#else
      detection_sum = detection_flag
      global_location_i = location_i
      global_location_j = location_j
      global_size_value = size_value
#endif

      if (detection_sum > 0.0) then
        grid%datetime_now = grid%datetime_start
        call grid%datetime_now%Add_seconds (grid%itimestep * grid%dt)
        call grid%Save_state ()

        write (msg, '(a, i6, a, i6, a, i6, a)') 'Isolated negative LFN cluster (size', int (global_size_value), ') near (i=', &
            int (global_location_i), ', j=', int (global_location_j), '). Simulation stopped after saving state.'
        call Stop_simulation (trim (msg))
      end if

    end subroutine Check_isolated_negative_lfn

    subroutine Prop_level_set (ifds, ifde, jfds, jfde, ifms, ifme, jfms, jfme, &
        num_tiles, i_start, i_end, j_start, j_end, ts, dt, dx, dy, fire_upwinding, fire_viscosity, &
        fire_viscosity_bg, fire_viscosity_band, fire_viscosity_ngp, fire_lsm_band_ngp, &
        tbound, lfn_in, lfn_0, lfn_1, lfn_2, lfn_out, tign, ros, uf, vf, dzdxf, dzdyf, ros_model, cart_comm, &
        ifps, ifpe, jfps, jfpe, grad_norm_ls, grad_norm_residual_sq_sum, grad_norm_residual_sq_sum_band, &
        grad_norm_residual_rms_band)

      ! Purpose: Advance the level set function from time ts to time ts + dt

      implicit none
      
      integer, intent(in) :: ifms, ifme, jfms, jfme, ifds, ifde, jfds, jfde, num_tiles, &
          fire_upwinding, fire_viscosity_ngp, fire_lsm_band_ngp, cart_comm, ifps, ifpe, jfps, jfpe
      integer, dimension (num_tiles), intent (in) :: i_start, i_end, j_start, j_end
      real, intent(in) :: fire_viscosity, fire_viscosity_bg, fire_viscosity_band
      real, dimension(ifms:ifme, jfms:jfme), intent (in) :: uf, vf, dzdxf, dzdyf
      real, dimension(ifms:ifme, jfms:jfme), intent (in out) :: lfn_in, tign, lfn_1, lfn_2, lfn_0
      real, dimension(ifms:ifme, jfms:jfme), intent (out) :: lfn_out, ros, grad_norm_ls
      real, intent (in) :: dx, dy, ts, dt
      real, intent (out) :: tbound, grad_norm_residual_sq_sum, grad_norm_residual_sq_sum_band, grad_norm_residual_rms_band
      class (ros_t), intent (in) :: ros_model

        ! to store tendency (rhs of the level set pde)
      real, dimension(ifms:ifme, jfms:jfme) :: tend
      real :: tbound2, tbound3, tbound_thread, tbound_min, grad_norm_residual_sq_sum_local, grad_norm_residual_sq_sum_band_local, &
          np_band, gmin, gsum
      integer :: i, j, ij, ifts, ifte, jfts, jfte, np_band_local
      character (len = 128) :: msg
      logical, parameter :: DEBUG_LOCAL = .false., PRINT_ERRORS = .false.


      if (DEBUG_LOCAL) call Print_message ('Entering sub Prop_level_set...')

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)
        do j = jfts, jfte
          do i = ifts, ifte
            lfn_0(i, j) = lfn_in(i, j)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
      call Do_halo_exchange (lfn_0, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

        ! Runge-Kutta step 1
      if (DEBUG_LOCAL) call Print_message ('call Calc_tend_ls 1...')

      tbound_min = huge (tbound_min)
      grad_norm_residual_sq_sum = 0.0
      grad_norm_residual_sq_sum_band = 0.0
      np_band = 0.0
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte, tbound_thread, grad_norm_residual_sq_sum_local, &
      !$OMP grad_norm_residual_sq_sum_band_local, np_band_local) &
      !$OMP REDUCTION(MIN: tbound_min) &
      !$OMP REDUCTION(+: grad_norm_residual_sq_sum, grad_norm_residual_sq_sum_band, np_band)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        call Calc_tend_ls (ifds, ifde, jfds, jfde, ifts, ifte, jfts, jfte, &
            ifms, ifme, jfms, jfme, ts, dt, dx, dy, fire_upwinding, &
            fire_viscosity, fire_viscosity_bg, fire_viscosity_band, &
            fire_viscosity_ngp, fire_lsm_band_ngp, lfn_0, tbound_thread, tend, ros, uf, vf, dzdxf, dzdyf, &
            ros_model, grad_norm_ls, grad_norm_residual_sq_sum_local, grad_norm_residual_sq_sum_band_local, np_band_local)

        tbound_min = min(tbound_min, tbound_thread)
        grad_norm_residual_sq_sum = grad_norm_residual_sq_sum + grad_norm_residual_sq_sum_local
        grad_norm_residual_sq_sum_band = grad_norm_residual_sq_sum_band + grad_norm_residual_sq_sum_band_local
        np_band = np_band + real (np_band_local)
      end do
      !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
      call Sum_across_mpi_tasks (grad_norm_residual_sq_sum, cart_comm, gsum)
      grad_norm_residual_sq_sum = gsum

      call Sum_across_mpi_tasks (grad_norm_residual_sq_sum_band, cart_comm, gsum)
      grad_norm_residual_sq_sum_band = gsum

      call Sum_across_mpi_tasks (np_band, cart_comm, gsum)
      np_band = gsum

      call Min_across_mpi_tasks (tbound_min, cart_comm, gmin)
      tbound_min = gmin
#endif

      tbound = tbound_min
      grad_norm_residual_sq_sum = sqrt (grad_norm_residual_sq_sum * dx * dy)
      grad_norm_residual_sq_sum_band = sqrt (grad_norm_residual_sq_sum_band * dx * dy)
      if (np_band > 0.0) then
        grad_norm_residual_rms_band = sqrt (grad_norm_residual_sq_sum_band / np_band)
      else
        grad_norm_residual_rms_band = 0.0
      end if

      if (PRINT_ERRORS) then
        write (msg, '(a, 3f12.6)') 'grad_norm_residual_sq_sum rk1, band, rms_band = ', grad_norm_residual_sq_sum, &
            grad_norm_residual_sq_sum_band, grad_norm_residual_rms_band
        call Print_message (msg)
      end if

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
            lfn_1(i, j) = lfn_0(i, j) + (dt / 3.0) * tend(i, j)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
      call Do_halo_exchange (lfn_1, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif


        ! Runge-Kutta step 2
      if (DEBUG_LOCAL) call Print_message ('call Calc_tend_ls 2...')

      tbound_min = huge (tbound_min)
      grad_norm_residual_sq_sum = 0.0
      grad_norm_residual_sq_sum_band = 0.0
      np_band = 0.0
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte, tbound_thread, grad_norm_residual_sq_sum_local, &
      !$OMP  grad_norm_residual_sq_sum_band_local, np_band_local) &
      !$OMP REDUCTION(MIN: tbound_min) &
      !$OMP REDUCTION(+: grad_norm_residual_sq_sum, grad_norm_residual_sq_sum_band, np_band)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        call Calc_tend_ls (ifds, ifde, jfds, jfde, ifts, ifte, jfts, jfte, &
            ifms,ifme,jfms,jfme, ts + dt, dt, dx, dy, fire_upwinding, &
            fire_viscosity, fire_viscosity_bg, fire_viscosity_band, &
            fire_viscosity_ngp, fire_lsm_band_ngp, lfn_1, tbound_thread, tend, ros, uf, vf, dzdxf, dzdyf, &
            ros_model, grad_norm_ls, grad_norm_residual_sq_sum_local, grad_norm_residual_sq_sum_band_local, np_band_local)

        tbound_min = min(tbound_min, tbound_thread)
        grad_norm_residual_sq_sum = grad_norm_residual_sq_sum + grad_norm_residual_sq_sum_local
        grad_norm_residual_sq_sum_band = grad_norm_residual_sq_sum_band + grad_norm_residual_sq_sum_band_local
        np_band = np_band + real (np_band_local)
      end do
      !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
      call Sum_across_mpi_tasks (grad_norm_residual_sq_sum, cart_comm, gsum)
      grad_norm_residual_sq_sum = gsum

      call Sum_across_mpi_tasks (grad_norm_residual_sq_sum_band, cart_comm, gsum)
      grad_norm_residual_sq_sum_band = gsum

      call Sum_across_mpi_tasks (np_band, cart_comm, gsum)
      np_band = gsum

      call Min_across_mpi_tasks (tbound_min, cart_comm, gmin)
      tbound_min = gmin
#endif

      tbound2 = tbound_min
      grad_norm_residual_sq_sum = sqrt (grad_norm_residual_sq_sum * dx * dy)
      grad_norm_residual_sq_sum_band = sqrt (grad_norm_residual_sq_sum_band * dx * dy)
      if (np_band > 0.0) then
        grad_norm_residual_rms_band = sqrt (grad_norm_residual_sq_sum_band / np_band)
      else
        grad_norm_residual_rms_band = 0.0
      end if

      if (PRINT_ERRORS) then
        write (msg, '(a, 3f12.6)') 'grad_norm_residual_sq_sum rk2, band, rms_band = ', grad_norm_residual_sq_sum, &
            grad_norm_residual_sq_sum_band, grad_norm_residual_rms_band
        call Print_message (msg)
      end if

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
            lfn_2(i, j) = lfn_0(i, j) + (dt / 2.0) * tend(i, j)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
      call Do_halo_exchange (lfn_2, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

        ! Runge-Kutta step 3
      if (DEBUG_LOCAL) call Print_message ('call Calc_tend_ls 3...')

      tbound_min = huge (tbound_min)
      grad_norm_residual_sq_sum = 0.0
      grad_norm_residual_sq_sum_band = 0.0
      np_band = 0.0
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte, tbound_thread, grad_norm_residual_sq_sum_local, &
      !$OMP grad_norm_residual_sq_sum_band_local, np_band_local) &
      !$OMP REDUCTION(MIN: tbound_min) &
      !$OMP REDUCTION(+: grad_norm_residual_sq_sum, grad_norm_residual_sq_sum_band, np_band)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        call Calc_tend_ls (ifds,ifde,jfds,jfde, ifts, ifte, jfts, jfte, &
            ifms, ifme, jfms, jfme, ts + dt, dt, dx, dy, fire_upwinding, &
            fire_viscosity, fire_viscosity_bg, fire_viscosity_band, &
            fire_viscosity_ngp, fire_lsm_band_ngp, lfn_2, tbound_thread, tend, ros, uf, vf, dzdxf, dzdyf, &
            ros_model, grad_norm_ls, grad_norm_residual_sq_sum_local, grad_norm_residual_sq_sum_band_local, np_band_local)

        tbound_min = min(tbound_min, tbound_thread)
        grad_norm_residual_sq_sum = grad_norm_residual_sq_sum + grad_norm_residual_sq_sum_local
        grad_norm_residual_sq_sum_band = grad_norm_residual_sq_sum_band + grad_norm_residual_sq_sum_band_local
        np_band = np_band + real (np_band_local)
      end do
      !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
      call Sum_across_mpi_tasks (grad_norm_residual_sq_sum, cart_comm, gsum)
      grad_norm_residual_sq_sum = gsum

      call Sum_across_mpi_tasks (grad_norm_residual_sq_sum_band, cart_comm, gsum)
      grad_norm_residual_sq_sum_band = gsum

      call Sum_across_mpi_tasks (np_band, cart_comm, gsum)
      np_band = gsum

      call Min_across_mpi_tasks (tbound_min, cart_comm, gmin)
      tbound_min = gmin
#endif

      tbound3 = tbound_min
      grad_norm_residual_sq_sum = sqrt (grad_norm_residual_sq_sum * dx * dy)
      grad_norm_residual_sq_sum_band = sqrt (grad_norm_residual_sq_sum_band * dx * dy)
      if (np_band > 0.0) then
        grad_norm_residual_rms_band = sqrt (grad_norm_residual_sq_sum_band / np_band)
      else
        grad_norm_residual_rms_band = 0.0
      end if
      if (PRINT_ERRORS) then
        write (msg, '(a, 3f12.6)') 'grad_norm_residual_sq_sum rk3, band, rms_band = ', grad_norm_residual_sq_sum, &
            grad_norm_residual_sq_sum_band, grad_norm_residual_rms_band
        call Print_message (msg)
      end if

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)
        do j = jfts, jfte
          do i = ifts, ifte
            lfn_out(i, j) = lfn_0(i, j) + dt * tend(i, j)
!            lfn_2(i,j) = lfn_out(i,j) ! lfn_2=lfn_out (needed for reinitialization purposes)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

!#ifdef DM_PARALLEL
!      include "HALO_FIRE_LFN_2.inc"
!#endif

        ! CFL check, tbound is the max allowed time step
      tbound = min(min(tbound, tbound2), tbound3)

      if (dt > tbound) then
        !$omp critical
        write (msg, '(a, f10.2, a, f10.2)') 'CFL violation: time step ', dt, ' > bound =', tbound
        call Stop_simulation (msg)
        !$omp end critical
      end if
    
      if (DEBUG_LOCAL) call Print_message ('leaving sub Prop_level_set...')

    end subroutine Prop_level_set

    subroutine Reinit_level_set (num_tiles, i_start, i_end, j_start, j_end, ifms, ifme, jfms, jfme, &
        ifds, ifde, jfds, jfde, ts, dt, dx, dy, fire_upwinding_reinit, &
        fire_lsm_reinit_iter, fire_lsm_band_ngp, lfn_in, lfn_2, lfn_s0, &
        lfn_s1, lfn_s2, lfn_s3, lfn_out, tign, cart_comm, &
        ifps, ifpe, jfps, jfpe, reinit_pseudot_coef, grad_norm_reinit)

    ! Purpose: Level-set function reinitialization
    !
    ! Referencess:
    ! Sussman, Smereka, Osher. Journal of Computational Physics 114, 146-159 (1994)
    !
    ! D. Munoz-Esparza, B. Kosovic, P. Jimenez, J. Coen: "An accurate
    ! fire-spread algorithm in the Weather Research and Forecasting model using the
    ! level-set method", Journal of Advances in Modeling Earth Systems, 2018
      

      implicit none

      integer, intent (in) :: num_tiles, cart_comm, ifps, ifpe, jfps, jfpe
      integer, dimension (num_tiles), intent (in) :: i_start, i_end, j_start, j_end
      integer, intent (in) :: ifms, ifme, jfms, jfme
      integer, intent (in) :: ifds, ifde, jfds, jfde
      integer, intent (in) :: fire_upwinding_reinit, fire_lsm_reinit_iter, fire_lsm_band_ngp
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: lfn_in, tign
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: lfn_2, lfn_s0, lfn_s1, lfn_s2, lfn_s3
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: lfn_out
      real, dimension (ifms:ifme, jfms:jfme), intent (out) :: grad_norm_reinit
      real, intent (in) :: reinit_pseudot_coef, dx, dy, ts, dt

      real :: dt_s, threshold_hlu
      integer :: nts, i, j, ij, ifts, ifte, jfts, jfte


      threshold_hlu = fire_lsm_band_ngp * dx

        ! Define S0 based on current lfn values
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
            lfn_s0(i, j) = lfn_out(i,j) / sqrt (lfn_out(i, j) ** 2.0 + dx ** 2.0)
            lfn_s3(i, j) = lfn_out(i,j)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
      call Do_halo_exchange (lfn_s3, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif
 
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        call Extrapol_var_at_bdys (ifms, ifme, jfms, jfme, ifds, ifde, &
            jfds, jfde, ifts, ifte, jfts, jfte, lfn_s3)
      end do
      !$OMP END PARALLEL DO

      dt_s = reinit_pseudot_coef * dx

        ! iterate to solve to steady state reinit PDE
        ! 1 iter each time step is enoguh
      Loop_iter: do nts = 1, fire_lsm_reinit_iter
          ! Runge-Kutta step 1
        !$OMP PARALLEL DO   &
        !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          call Advance_ls_reinit (ifms, ifme, jfms, jfme, ifds, ifde, jfds, jfde, &
              ifts, ifte, jfts, jfte, dx, dy, dt_s, threshold_hlu, &
              lfn_s0, lfn_s3, lfn_s3, lfn_s1, 1.0 / 3.0, & ! sign funcition, initial ls, current stage ls, next stage advanced ls, RK coefficient
              fire_upwinding_reinit, grad_norm_reinit)
        end do
        !$OMP END PARALLEL DO
 
#ifdef DM_PARALLEL
      call Do_halo_exchange (lfn_s1, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

        !$OMP PARALLEL DO   &
        !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          call Extrapol_var_at_bdys (ifms, ifme, jfms, jfme, ifds, ifde, &
              jfds, jfde, ifts, ifte, jfts, jfte, lfn_s1)
        end do
        !$OMP END PARALLEL DO

          ! Runge-Kutta step 2
        !$OMP PARALLEL DO   &
        !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          call Advance_ls_reinit (ifms, ifme, jfms, jfme, ifds, ifde, jfds, jfde, &
              ifts, ifte, jfts, jfte, dx, dy, dt_s, threshold_hlu, &
              lfn_s0, lfn_s3, lfn_s1, lfn_s2, 1.0 / 2.0, &
              fire_upwinding_reinit, grad_norm_reinit)
        end do
        !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
      call Do_halo_exchange (lfn_s2, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

        !$OMP PARALLEL DO   &
        !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          call Extrapol_var_at_bdys (ifms, ifme, jfms, jfme, ifds, ifde, &
              jfds, jfde, ifts, ifte, jfts, jfte, lfn_s2)
        end do
        !$OMP END PARALLEL DO

          ! Runge-Kutta step 3
        !$OMP PARALLEL DO   &
        !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          call Advance_ls_reinit (ifms, ifme, jfms, jfme, ifds, ifde, jfds, jfde, &
              ifts, ifte, jfts, jfte, dx, dy, dt_s, threshold_hlu, &
              lfn_s0, lfn_s3, lfn_s2, lfn_s3, 1.0, &
              fire_upwinding_reinit, grad_norm_reinit)
        end do
        !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
      call Do_halo_exchange (lfn_s3, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

        !$OMP PARALLEL DO   &
        !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          call Extrapol_var_at_bdys (ifms, ifme, jfms, jfme, ifds, ifde, &
              jfds,jfde,  ifts, ifte, jfts, jfte, lfn_s3)
        end do
        !$OMP END PARALLEL DO
      end do Loop_iter

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
              ! assing to lfn_out the reinitialized level-set function
            lfn_out(i, j) = lfn_s3(i, j)
              ! fire area can only increase
            lfn_out(i, j) = min (lfn_out(i, j), lfn_in(i, j))
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine Reinit_level_set

    subroutine Advance_ls_reinit (ifms, ifme, jfms, jfme, ifds, ifde, jfds, jfde, &
        ifts, ifte, jfts, jfte, dx, dy, dt_s, threshold_hlu, lfn_s0, &
        lfn_ini, lfn_curr, lfn_fin, rk_coeff, fire_upwinding_reinit, grad_norm_reinit)

      ! Calculates right-hand-side forcing and advances a RK-stage the level-set reinitialization PDE

      implicit none

      integer, intent (in) :: ifms, ifme, jfms, jfme, ifts, ifte, jfts, &
          jfte, ifds, ifde, jfds, jfde
      integer, intent (in) :: fire_upwinding_reinit
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: lfn_s0, lfn_ini, lfn_curr
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: lfn_fin
      real, dimension (ifms:ifme, jfms:jfme), intent (out) :: grad_norm_reinit
      real, intent (in) :: dx, dy, dt_s, threshold_hlu, rk_coeff

      integer :: i, j
      real :: diffLx, diffLy, diffRx, diffRy, diff2x, diff2y, grad, tend_r


      do j = jfts, jfte 
        do i = ifts, ifte 
          if (i < ifds + BDY_ENO1 .or. i > ifde - BDY_ENO1 .or. &
              j < jfds + BDY_ENO1 .or. j > jfde - BDY_ENO1) then
            diffLx = (lfn_curr(i, j) - lfn_curr(i - 1, j)) / dx
            diffLy = (lfn_curr(i, j) - lfn_curr(i, j - 1)) / dy
            diffRx = (lfn_curr(i + 1, j) - lfn_curr(i, j)) / dx
            diffRy = (lfn_curr(i, j + 1) - lfn_curr(i, j)) / dy
            diff2x = Select_eno (diffLx, diffRx)
            diff2y = Select_eno (diffLy, diffRy)
          else
            select case (fire_upwinding_reinit)
              case (1)
                diff2x = Select_4th (dx, lfn_curr(i, j), lfn_curr(i - 1, j), &
                    lfn_curr(i - 2, j), lfn_curr(i + 1, j), lfn_curr(i + 2, j))
                diff2y = Select_4th (dy, lfn_curr(i, j), lfn_curr(i, j - 1), &
                    lfn_curr(i, j - 2), lfn_curr(i, j + 1), lfn_curr(i, j + 2))
                diff2x = Select_weno3 (dx, lfn_curr(i, j), lfn_curr(i - 1, j), &
                    lfn_curr(i - 2, j), lfn_curr(i + 1, j), lfn_curr(i + 2, j), &
                    lfn_s0(i, j) * diff2x)
                diff2y = Select_weno3 (dy, lfn_curr(i, j), lfn_curr(i, j - 1), &
                    lfn_curr(i, j - 2), lfn_curr(i, j + 1), lfn_curr(i, j + 2), &
                    lfn_s0(i, j) * diff2y)

              case (2)
                diff2x = Select_4th (dx, lfn_curr(i, j), lfn_curr(i - 1, j), & 
                    lfn_curr(i - 2, j), lfn_curr(i + 1, j), lfn_curr(i + 2, j))
                diff2y = Select_4th (dy, lfn_curr(i, j), lfn_curr(i, j - 1), &
                    lfn_curr(i, j - 2), lfn_curr(i, j + 1), lfn_curr(i, j + 2))
                diff2x = Select_weno5 (dx, lfn_curr(i, j), lfn_curr(i - 1, j), &
                    lfn_curr(i - 2, j), lfn_curr(i - 3, j), lfn_curr(i + 1, j), &
                    lfn_curr(i + 2, j), lfn_curr(i + 3, j), lfn_s0(i, j) * diff2x)
                diff2y = Select_weno5 (dy, lfn_curr(i, j), lfn_curr(i, j - 1), &
                    lfn_curr(i, j - 2), lfn_curr(i, j - 3), lfn_curr(i, j + 1), &
                    lfn_curr(i, j + 2), lfn_curr(i, j + 3), lfn_s0(i, j) * diff2y)

              case (3)
                if (lfn_curr(i, j) < threshold_hlu) then
                  diff2x = Select_4th (dx, lfn_curr(i, j), lfn_curr(i - 1, j), &
                      lfn_curr(i - 2, j), lfn_curr(i + 1, j), lfn_curr(i + 2, j))
                  diff2y = Select_4th (dy, lfn_curr(i, j), lfn_curr(i, j - 1), & 
                      lfn_curr(i, j - 2), lfn_curr(i, j + 1), lfn_curr(i, j + 2))
                  diff2x = Select_weno3 (dx, lfn_curr(i, j), lfn_curr(i - 1, j), &
                      lfn_curr(i - 2, j), lfn_curr(i + 1, j), lfn_curr(i + 2, j), &
                      lfn_s0(i, j) * diff2x)
                  diff2y = Select_weno3 (dy, lfn_curr(i, j), lfn_curr(i, j - 1), &
                      lfn_curr(i, j - 2), lfn_curr(i, j + 1), lfn_curr(i, j + 2), &
                      lfn_s0(i, j) * diff2y)
                else
                  diffLx = (lfn_curr(i, j) - lfn_curr(i - 1, j)) / dx
                  diffLy = (lfn_curr(i, j) - lfn_curr(i, j - 1)) / dy
                  diffRx = (lfn_curr(i + 1, j) - lfn_curr(i, j)) / dx
                  diffRy = (lfn_curr(i, j + 1) - lfn_curr(i, j)) / dy
                  diff2x = Select_eno (diffLx, diffRx)
                  diff2y = Select_eno (diffLy, diffRy)
                endif

              case(4)
                if (lfn_curr(i, j) < threshold_hlu) then
                  diff2x = Select_4th (dx, lfn_curr(i, j), lfn_curr(i - 1, j), &
                      lfn_curr(i - 2, j), lfn_curr(i + 1, j), lfn_curr(i + 2, j))
                  diff2y = Select_4th (dy,lfn_curr(i, j), lfn_curr(i, j - 1), &
                      lfn_curr(i, j - 2), lfn_curr(i, j + 1), lfn_curr(i, j + 2))
                  diff2x = Select_weno5 (dx, lfn_curr(i, j), lfn_curr(i - 1, j), &
                      lfn_curr(i - 2, j), lfn_curr(i - 3, j), lfn_curr(i + 1, j), &
                      lfn_curr(i + 2, j), lfn_curr(i + 3, j), lfn_s0(i, j) * diff2x)
                  diff2y = Select_weno5 (dy, lfn_curr(i,j), lfn_curr(i, j - 1), &
                      lfn_curr(i, j - 2), lfn_curr(i, j - 3), lfn_curr(i, j + 1), &
                      lfn_curr(i, j + 2), lfn_curr(i, j + 3), lfn_s0(i, j) * diff2y)
                else
                  diffLx = (lfn_curr(i, j) - lfn_curr(i - 1, j)) / dx
                  diffLy = (lfn_curr(i, j) - lfn_curr(i, j - 1)) / dy
                  diffRx = (lfn_curr(i + 1, j) - lfn_curr(i, j)) / dx
                  diffRy = (lfn_curr(i, j + 1) - lfn_curr(i, j)) / dy
                  diff2x = Select_eno (diffLx, diffRx)
                  diff2y = Select_eno (diffLy, diffRy)
                endif

              case default
                if (lfn_curr(i,j) < threshold_hlu) then
                  diff2x = Select_4th (dx, lfn_curr(i, j), lfn_curr(i - 1, j), &
                      lfn_curr(i - 2, j), lfn_curr(i + 1, j), lfn_curr(i + 2, j))
                  diff2y = Select_4th (dy, lfn_curr(i, j), lfn_curr(i, j - 1), &
                      lfn_curr(i, j - 2), lfn_curr(i, j + 1), lfn_curr(i, j + 2))
                  diff2x = Select_weno5 (dx, lfn_curr(i, j), lfn_curr(i - 1, j), &
                      lfn_curr(i - 2, j), lfn_curr(i - 3, j), lfn_curr(i + 1, j), &
                      lfn_curr(i + 2, j), lfn_curr(i + 3, j), lfn_s0(i, j) * diff2x)
                  diff2y = Select_weno5 (dy, lfn_curr(i, j), lfn_curr(i, j - 1), &
                      lfn_curr(i, j - 2), lfn_curr(i, j - 3), lfn_curr(i, j + 1), &
                      lfn_curr(i, j + 2), lfn_curr(i, j + 3), lfn_s0(i, j) * diff2y)
                else
                  diffLx = (lfn_curr(i, j) - lfn_curr(i - 1, j)) / dx
                  diffLy = (lfn_curr(i, j) - lfn_curr(i, j - 1)) / dy
                  diffRx = (lfn_curr(i + 1, j) - lfn_curr(i, j)) / dx
                  diffRy = (lfn_curr(i, j + 1) - lfn_curr(i, j)) / dy
                  diff2x = Select_eno (diffLx, diffRx)
                  diff2y = Select_eno (diffLy, diffRy)
                end if

            end select
          end if
            grad = sqrt (diff2x * diff2x + diff2y * diff2y)
            grad_norm_reinit(i, j) = grad
            tend_r = lfn_s0(i, j) * (1.0 - grad)
            lfn_fin(i, j) = lfn_ini(i, j) + (dt_s * rk_coeff) * tend_r
        end do
      end do

    end subroutine Advance_ls_reinit

    subroutine Update_ignition_times (ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, &
        ifds, jfds, ifde, jfde, ts, dt, lfn_in, lfn_out, tign)

      ! compute ignition time for a node that was not burning at start but it is burning at end

      implicit none

      integer, intent (in) :: ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, &
          ifds, jfds, ifde, jfde
      real, intent(in) :: ts, dt
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: lfn_in, lfn_out
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: tign
      real :: time_now
      integer :: i, j


      time_now = ts + dt
      time_now = time_now + abs (time_now) * epsilon (time_now) * 2.0
      do j = jfts, jfte
        do i = ifts, ifte
            ! interpolate the cross-over time
          if (lfn_out(i, j) <= 0.0 .and. lfn_in(i, j) > 0.0) then
            tign(i, j) = ts + dt * lfn_in(i, j) / (lfn_in(i, j) - lfn_out(i, j))
          end if
          if (lfn_out(i,j) > 0.0) tign(i, j) = time_now
        end do
      end do

    end subroutine Update_ignition_times

    subroutine Calc_tend_ls (ids, ide, jds, jde, its, ite, jts, jte, ifms, ifme, jfms, jfme, &
        t, dt, dx, dy, fire_upwinding, fire_viscosity, fire_viscosity_bg, &
        fire_viscosity_band, fire_viscosity_ngp, fire_lsm_band_ngp, lfn, tbound, tend, ros, uf, vf, dzdxf, dzdyf, &
        ros_model, grad_norm_ls, grad_norm_residual_sq_sum_local, grad_norm_residual_sq_sum_band_local, np_band_local)

      ! compute the right hand side of the level set equation

      implicit none

      integer, intent (in) :: ifms, ifme, jfms, jfme, its, ite, jts, jte, ids, ide, jds, jde, &
          fire_upwinding,fire_viscosity_ngp, fire_lsm_band_ngp
      integer, intent (out) :: np_band_local
      real, intent (in) :: fire_viscosity, fire_viscosity_bg, fire_viscosity_band, t, dt, dx, dy
      real, dimension(ifms:ifme, jfms:jfme), intent (in) :: uf, vf, dzdxf, dzdyf
      real, dimension(ifms:ifme, jfms:jfme), intent (in out) :: lfn
      real, dimension(ifms:ifme, jfms:jfme), intent (out) :: tend, ros, grad_norm_ls
      real, intent (out) :: tbound, grad_norm_residual_sq_sum_local, grad_norm_residual_sq_sum_band_local
      class (ros_t), intent (in) :: ros_model

      real, parameter :: EPS = epsilon (0.0), TOL = 100.0 * EPS
      real :: difflx, diffly, diffrx, diffry, diffcx, diffcy, &
         diff2x, diff2y, grad, mask, diff2x_eno, diff2y_eno, &
         diff2x_weno, diff2y_weno, transition_width, band_width, &
         scale, nvx, nvy, a_valor, signo_x, signo_y, threshold_hll, &
         threshold_hlu, threshold_av, fire_viscosity_var
      integer :: i, j
      character (len = :), allocatable :: msg
      logical, parameter :: DEBUG_LOCAL = .false.


      if (DEBUG_LOCAL) call Print_message ('Entering subroutine Calc_tend_ls')

      threshold_hll = -fire_lsm_band_ngp * dx
      threshold_hlu = fire_lsm_band_ngp * dx
      threshold_av = fire_viscosity_ngp * dx

      if (DEBUG_LOCAL) call Print_message ('calling Extrapol_var_at_bdys')
      call Extrapol_var_at_bdys (ifms, ifme, jfms, jfme, ids, ide, jds, jde, &
          its, ite, jts, jte, lfn)

      if (DEBUG_LOCAL) call Print_message ('starting ij loops')
      tbound = 0.0
      grad_norm_residual_sq_sum_local = 0.0
      grad_norm_residual_sq_sum_band_local = 0.0
      np_band_local = 0
      do j = jts, jte
        do i = its, ite
            ! one sided differences
          diffrx = (lfn(i + 1, j) - lfn(i, j)) / dx
          difflx = (lfn(i, j) - lfn(i - 1, j)) / dx
          diffry = (lfn(i, j + 1) - lfn(i, j)) / dy
          diffly = (lfn(i, j) - lfn(i, j - 1)) / dy

            ! central differences
          diffcx = 0.5 * (difflx + diffrx)
          diffcy = 0.5 * (diffly + diffry)

            ! use eno1 near domain boundaries
          if (i < ids + BDY_ENO1 .or. i > ide - BDY_ENO1 .or. &
              j < jds + BDY_ENO1 .or. j > jde - BDY_ENO1) then 
            diff2x = Select_eno (difflx, diffrx)
            diff2y = Select_eno (diffly, diffry)
            grad = sqrt (diff2x * diff2x + diff2y * diff2y)

          else
            select case (fire_upwinding)
              case (0)
                  ! none
                diff2x = diffcx
                diff2y = diffcy
                grad = sqrt (diff2x * diff2x + diff2y * diff2y)

              case (1)
                  ! First-order sign-consistent upwind scheme (component wise). Standard
                diff2x = Select_upwind (difflx, diffrx)
                diff2y = Select_upwind (diffly, diffry)
                grad = sqrt (diff2x * diff2x + diff2y * diff2y)

              case (2)
                  ! Component-wise Godunov upwind (1st order) derivative selection
                diff2x = Select_godunov (difflx, diffrx)
                diff2y = Select_godunov (diffly, diffry)
                grad = sqrt (diff2x * diff2x + diff2y * diff2y)

              case (3)
                  ! First order ENO (minmod limited) componenet-wise gradient
                diff2x = Select_eno (difflx, diffrx)
                diff2y = Select_eno (diffly, diffry)
                grad = sqrt (diff2x * diff2x + diff2y * diff2y)

              case (4)
                  ! Classical Godunov discretization of the Hamiltonian
                grad = sqrt (max (difflx, 0.0) ** 2 + min (diffrx, 0.0) ** 2 + &
                    max (diffly, 0.0) ** 2 + min (diffry, 0.0) ** 2)
                diff2x = max (difflx, 0.0) - min (diffrx, 0.0)
                diff2y = max (diffly, 0.0) - min (diffry, 0.0)

              case (5)
                  ! 2nd order central differences
                diff2x = Select_2nd (dx, lfn(i, j), lfn(i - 1, j), lfn(i + 1, j))
                diff2y = Select_2nd (dy, lfn(i, j), lfn(i, j - 1), lfn(i, j + 1))
                grad = sqrt (diff2x * diff2x + diff2y * diff2y)

              case (6)
                  ! 3rd order WENO
                a_valor = Select_4th (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j)) * uf(i, j) + &
                    Select_4th (dy, lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), lfn(i, j + 1), lfn(i, j + 2)) * vf(i, j)
                signo_x = a_valor * Select_4th (dx, lfn(i, j), lfn(i - 1, j), &
                    lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j))
                signo_y = a_valor * Select_4th (dy, lfn(i, j), lfn(i, j - 1), &
                    lfn(i, j - 2), lfn(i, j + 1), lfn(i, j + 2))
                diff2x = Select_weno3 (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), &
                    lfn(i + 1, j), lfn(i + 2, j), signo_x)
                diff2y = Select_weno3 (dy, lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), & 
                    lfn(i, j + 1), lfn(i, j + 2), signo_y)
                grad = sqrt (diff2x * diff2x + diff2y * diff2y)

              case (7)
                  ! 5th order WENO
                a_valor = Select_4th (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j)) * uf(i, j)+ &
                    Select_4th (dy, lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), lfn(i, j + 1), lfn(i, j + 2)) * vf(i, j)
                signo_x = a_valor * Select_4th (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j))
                signo_y = a_valor * Select_4th (dy,lfn(i,j),lfn(i,j-1),lfn(i,j-2),lfn(i,j+1),lfn(i,j+2))
                diff2x = Select_weno5 (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), &
                    lfn(i - 3, j), lfn(i + 1, j), lfn(i + 2, j), lfn(i + 3, j), signo_x)
                diff2y = Select_weno5 (dy, lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), &
                    lfn(i ,j - 3), lfn(i, j + 1), lfn(i, j + 2), lfn(i, j + 3), signo_y)
                grad = sqrt (diff2x * diff2x + diff2y * diff2y)

              case (8)
                  ! WENO3/ENO1
                if (abs (lfn(i, j)) < threshold_hlu) then
                  a_valor = Select_4th (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j)) * uf(i, j) + &
                      Select_4th (dy, lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), lfn(i, j + 1), lfn(i, j + 2)) * vf(i, j)
                  signo_x = a_valor * Select_4th (dx, lfn(i, j), lfn(i - 1, j), & 
                      lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j))
                  signo_y = a_valor * Select_4th (dy, lfn(i, j), lfn(i, j - 1), &
                      lfn(i, j - 2), lfn(i, j + 1), lfn(i, j + 2))
                  diff2x = Select_weno3 (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), &
                      lfn(i + 1, j), lfn(i + 2, j), signo_x)
                  diff2y = Select_weno3 (dy, lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), &
                      lfn(i, j + 1), lfn(i, j + 2), signo_y)
                  grad = sqrt (diff2x * diff2x + diff2y * diff2y)
                else
                  diff2x = Select_eno (difflx, diffrx)
                  diff2y = Select_eno (diffly, diffry)
                  grad = sqrt (diff2x * diff2x + diff2y * diff2y)
                end if

              case (9)
                  ! WENO5/ENO1
                if (abs (lfn(i, j)) < threshold_hlu) then
                  a_valor = Select_4th (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j)) * uf(i, j) + &
                      Select_4th (dy,lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), lfn(i, j + 1), lfn(i, j + 2)) * vf(i, j)
                  signo_x = a_valor * Select_4th (dx, lfn(i, j), lfn(i - 1, j), &
                      lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j))
                  signo_y = a_valor * Select_4th (dy, lfn(i, j), lfn(i, j - 1), &
                      lfn(i, j - 2), lfn(i, j + 1), lfn(i, j + 2))
                  diff2x = Select_weno5 (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), &
                      lfn(i - 3, j), lfn(i + 1, j), lfn(i + 2, j), lfn(i + 3, j), signo_x)
                  diff2y = Select_weno5 (dy, lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), &
                      lfn(i, j - 3), lfn(i, j + 1), lfn(i, j + 2), lfn(i, j + 3), signo_y)
                  grad = sqrt (diff2x * diff2x + diff2y * diff2y)
                else
                  diff2x = Select_eno (difflx, diffrx)
                  diff2y = Select_eno (diffly, diffry)
                  grad = sqrt (diff2x * diff2x + diff2y * diff2y)
                end if

              case (10)
                  ! WENO5/ENO1 + blending zone
                band_width = 2 * threshold_hlu
                transition_width = threshold_hlu
                if (abs(lfn(i, j)) <= band_width) then
                  if (abs(lfn(i, j)) <= band_width - transition_width) then
                    mask = 1.0
                  else
                    mask = 0.5 * (1.0 + cos (PI * (abs (lfn(i, j)) -(band_width - transition_width)) / transition_width))
                  end if
                else
                  mask = 0.0
                end if

                diff2x_eno = Select_eno (difflx, diffrx)
                diff2y_eno = Select_eno (diffly, diffry)

                if (mask > 0.0) then
                  a_valor = Select_4th (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j)) * uf(i, j) + &
                      Select_4th (dy,lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), lfn(i, j + 1), lfn(i, j + 2)) * vf(i, j)
                  signo_x = a_valor * Select_4th (dx, lfn(i, j), lfn(i - 1, j), &
                      lfn(i - 2, j), lfn(i + 1, j), lfn(i + 2, j))
                  signo_y = a_valor * Select_4th (dy, lfn(i, j), lfn(i, j - 1), &
                      lfn(i, j - 2), lfn(i, j + 1), lfn(i, j + 2))
                  diff2x_weno = Select_weno5 (dx, lfn(i, j), lfn(i - 1, j), lfn(i - 2, j), &
                      lfn(i - 3, j), lfn(i + 1, j), lfn(i + 2, j), lfn(i + 3, j), signo_x)
                  diff2y_weno = Select_weno5 (dy, lfn(i, j), lfn(i, j - 1), lfn(i, j - 2), &
                      lfn(i, j - 3), lfn(i, j + 1), lfn(i, j + 2), lfn(i, j + 3), signo_y)
                  diff2x = mask * diff2x_weno + (1.0 - mask) * diff2x_eno
                  diff2y = mask * diff2y_weno + (1.0 - mask) * diff2y_eno
                else
                  diff2x = diff2x_eno
                  diff2y = diff2y_eno
                end if

                grad = sqrt (diff2x * diff2x + diff2y * diff2y)

              case default
                !$omp critical
                write (msg, '(a, i2)') 'Unknown upwinding option in level set : ', fire_upwinding
                call Stop_simulation (msg)
                !$omp end critical

            end select
          end if

          grad_norm_ls(i, j) = grad
          grad_norm_residual_sq_sum_local = grad_norm_residual_sq_sum_local + (grad - 1.0) ** 2
          if (abs(lfn(i, j)) < threshold_hlu) then
            grad_norm_residual_sq_sum_band_local = grad_norm_residual_sq_sum_band_local + (grad - 1.0) ** 2
            np_band_local = np_band_local + 1
          end if

            ! Calc normal
          scale = sqrt (grad ** 2.0 + EPS)
          nvx = diff2x / scale
          nvy = diff2y / scale

            ! Get rate of spread from wind speed and slope
          ros(i, j) = ros_model%Calc_ros (ifms, ifme, jfms, jfme, i, j, &
              nvx, nvy, uf(i, j), vf(i, j), dzdxf(i, j), dzdyf(i, j))

            ! CFL condition
          if (grad > 0.0) tbound = max (tbound, ros(i, j) * (abs (diff2x) / dx + &
              abs (diff2y) / dy) / grad)

            ! Tendency level set function
          tend(i, j) = -ros(i, j) * grad

            ! Add to tend effect Artificial viscosity
          if (abs (lfn(i,j)) < threshold_av .and. (i > ids + BDY_ENO1 .and. i < ide - BDY_ENO1) .and. &
              (j > jds + BDY_ENO1 .and. j < jde - BDY_ENO1)) then 
            fire_viscosity_var = fire_viscosity_bg
          else if (abs (lfn(i,j)) >= threshold_av .and. abs (lfn(i,j)) < threshold_av * (1.0 + fIre_viscosity_band) .and. &
              (i > ids + BDY_ENO1 .and. i < ide - BDY_ENO1) .and. (j > jds + BDY_ENO1 .and. j < jde - BDY_ENO1)) then
            fire_viscosity_var = min (fire_viscosity_bg + (fire_viscosity - fire_viscosity_bg) * &
                (abs (lfn(i, j)) - threshold_av) / (fire_viscosity_band * threshold_av), fire_viscosity)
          else
            fire_viscosity_var = fire_viscosity
          end if

          tend(i, j) = tend(i, j) + fire_viscosity_var * abs (ros(i, j)) * ((diffrx - difflx) + (diffry - diffly))
        end do
      end do

        ! final CFL bound
      tbound = 1.0 / (tbound + TOL)

      if (DEBUG_LOCAL) call Print_message ('Leaving subroutine Calc_tend_ls')

    end subroutine Calc_tend_ls

    subroutine Reinit_level_set_fast_dist (lfn_s0, lfn_out, i_start, i_end, j_start, j_end, ifms, ifme, jfms, &
            jfme, num_tiles, fast_dist_reinit_opt, dx, dy, ifps, ifpe, jfps, jfpe, ifds, ifde, jfds, jfde, cart_comm)

      implicit none

      integer, intent (in) :: ifms, ifme, jfms, jfme, num_tiles, fast_dist_reinit_opt, &
          ifps, ifpe, jfps, jfpe, ifds, ifde, jfds, jfde, cart_comm
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: lfn_s0, lfn_out
      integer, dimension(num_tiles), intent(in) :: i_start, i_end, j_start, j_end
      real, intent (in) :: dx, dy

      real, parameter :: INF = Huge (1.0)
      integer :: i, j, ij, ifts, ifte, jfts, jfte
      integer, dimension (ifms:ifme, jfms:jfme) :: status

      logical, parameter :: DEBUG_LOCAL = .false.


      if (DEBUG_LOCAL) call Print_message ('Entering Reinit_level_set_fast_dist')

        ! Init arrays
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
            lfn_s0(i, j) = INF
            status(i, j) = FAR
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      do j = jfms, jfme
        do i = ifms, ifme
          if (i < ifds .or. i > ifde .or. j < jfds .or. j > jfde) then
            lfn_s0(i, j) = INF
            status(i, j) = FAR
          end if
        end do
      end do

#ifdef DM_PARALLEL
      call Do_halo_exchange (lfn_s0, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

      ! Update halos and boundaries for lfn_out
#ifdef DM_PARALLEL
      call Do_halo_exchange (lfn_out, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        call Extrapol_var_at_bdys (ifms, ifme, jfms, jfme, ifds, ifde, &
            jfds, jfde, ifts, ifte, jfts, jfte, lfn_out)
      end do
      !$OMP END PARALLEL DO

        ! Set level set function near the interface
      if (DEBUG_LOCAL) call Print_message ('calling Init level set function at the interface...')
      call Init_level_set_at_interface (lfn_s0, status, lfn_out, i_start, i_end, &
          j_start, j_end, num_tiles, ifms, ifme, jfms, jfme, dx, dy)

        ! Set level set function in the rest of the grid
      select case (fast_dist_reinit_opt)
        case (FAST_DIST_REINIT_FSM)
          if (DEBUG_LOCAL) call Print_message ('calling FSM method...')
          call Reinit_ls_fsm (lfn_s0, status, i_start, i_end, j_start, j_end, num_tiles, dx, dy, &
              ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, cart_comm)

          !$OMP PARALLEL DO   &
          !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
          do ij = 1, num_tiles
            ifts = i_start(ij)
            ifte = i_end(ij)
            jfts = j_start(ij)
            jfte = j_end(ij)

            do j = jfts, jfte
              do i = ifts, ifte
                lfn_out(i, j) = sign (1.0, lfn_out(i, j)) * lfn_s0(i, j)
              end do
            end do
          end do
          !$OMP END PARALLEL DO

        case default
          call Stop_simulation ('Fast distatnce reinitialization option not valid.')

      end select

      if (DEBUG_LOCAL) call Print_message ('Leaving Reinit_level_set_fast_dist')

    end subroutine Reinit_level_set_fast_dist

    subroutine Init_level_set_at_interface (lfn_s0, status, lfn_out, i_start, i_end, &
        j_start, j_end, num_tiles, ifms, ifme, jfms, jfme, dx, dy)

      implicit none

      integer, intent(in) :: ifms, ifme, jfms, jfme, num_tiles
      integer, dimension(num_tiles), intent(in) :: i_start, i_end, j_start, j_end
      real, dimension(ifms:ifme, jfms:jfme), intent(in) :: lfn_out
      real, dimension(ifms:ifme, jfms:jfme), intent(in out) :: lfn_s0
      integer, dimension(ifms:ifme, jfms:jfme), intent(in out) :: status
      real, intent (in) :: dx, dy

      integer, parameter :: INIT_DEFAULT=0, INIT_ZERO = 1, INIT_SUBGRID = 2, INIT_QUADRATIC = 3, &
          INIT_GRADIENT = 4
      integer :: init_method


      init_method = INIT_GRADIENT

      select case (init_method)
        case (INIT_DEFAULT)
          call Set_ls_interface_default (lfn_s0, status, lfn_out, i_start, i_end, &
              j_start, j_end, num_tiles, ifms, ifme, jfms, jfme)

        case (INIT_ZERO)
          call Set_ls_interface_zero (lfn_s0, status, lfn_out, i_start, i_end, &
              j_start, j_end, num_tiles, ifms, ifme, jfms, jfme)

        case (INIT_SUBGRID)
          call Set_ls_interface_subgrid_inter (lfn_s0, status, lfn_out, i_start, i_end, &
              j_start, j_end, num_tiles, ifms, ifme, jfms, jfme, dx, dy)

        case (INIT_QUADRATIC)
          call Set_ls_interface_quadratic_inter (lfn_s0, status, lfn_out, i_start, i_end, &
              j_start, j_end, num_tiles, ifms, ifme, jfms, jfme, dx, dy)

        case (INIT_GRADIENT)
          call Set_ls_interface_gradient_inter (lfn_s0, status, lfn_out, i_start, i_end, &
              j_start, j_end, num_tiles, ifms, ifme, jfms, jfme, dx, dy)

        case default
          call Stop_simulation ('Error: phi_init method not available')

      end select

    end subroutine Init_level_set_at_interface

    subroutine Set_ls_interface_default (lfn_s0, status, lfn_out, i_start, i_end, &
        j_start, j_end, num_tiles, ifms, ifme, jfms, jfme)

      implicit none

      integer, intent(in) :: ifms, ifme, jfms, jfme, num_tiles
      integer, dimension(num_tiles), intent(in) :: i_start, i_end, j_start, j_end
      real, dimension(ifms:ifme, jfms:jfme), intent(in) :: lfn_out
      real, dimension(ifms:ifme, jfms:jfme), intent(in out) :: lfn_s0
      integer, dimension(ifms:ifme, jfms:jfme), intent(in out) :: status

      integer :: ifts, ifte, jfts, jfte, i, j, ij


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
            if (lfn_out(i, j) * lfn_out(i + 1, j) < 0.0 .or. &
                lfn_out(i, j) * lfn_out(i - 1, j) < 0.0 .or. &
                lfn_out(i, j) * lfn_out(i, j + 1) < 0.0 .or. &
                lfn_out(i, j) * lfn_out(i, j - 1) < 0.0) then
              lfn_s0(i, j) = lfn_out(i, j)
              status(i, j) = KNOWN
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine Set_ls_interface_default

    subroutine Set_ls_interface_zero (lfn_s0, status, lfn_out, i_start, i_end, &
        j_start, j_end, num_tiles, ifms, ifme, jfms, jfme)

      implicit none

      integer, intent(in) :: ifms, ifme, jfms, jfme, num_tiles
      integer, dimension(num_tiles), intent(in) :: i_start, i_end, j_start, j_end
      real, dimension(ifms:ifme, jfms:jfme), intent(in) :: lfn_out
      real, dimension(ifms:ifme, jfms:jfme), intent(in out) :: lfn_s0
      integer, dimension(ifms:ifme, jfms:jfme), intent(in out) :: status

      integer :: ifts, ifte, jfts, jfte, i, j, ij


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
            if (lfn_out(i, j) * lfn_out(i + 1, j) < 0.0 .or. &
                lfn_out(i, j) * lfn_out(i - 1, j) < 0.0 .or. &
                lfn_out(i, j) * lfn_out(i, j + 1) < 0.0 .or. &
                lfn_out(i, j) * lfn_out(i, j - 1) < 0.0) then
              lfn_s0(i, j) = 0.0
              status(i, j) = KNOWN
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine Set_ls_interface_zero

    subroutine Set_ls_interface_subgrid_inter (lfn_s0, status, lfn_out, i_start, i_end, &
        j_start, j_end, num_tiles, ifms, ifme, jfms, jfme, dx, dy)

      implicit none

      integer, intent(in) :: ifms, ifme, jfms, jfme, num_tiles
      integer, dimension(num_tiles), intent(in) :: i_start, i_end, j_start, j_end
      real, dimension(ifms:ifme, jfms:jfme), intent(in) :: lfn_out
      real, dimension(ifms:ifme, jfms:jfme), intent(in out) :: lfn_s0
      integer, dimension(ifms:ifme, jfms:jfme), intent(in out) :: status
      real, intent (in) :: dx, dy

      integer :: ifts, ifte, jfts, jfte, i, j, ij


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
            if (lfn_out(i, j) * lfn_out(i + 1, j) < 0.0) then
              lfn_s0(i, j) = abs (lfn_out(i, j)) * dx / (abs (lfn_out(i, j)) + abs (lfn_out(i + 1, j)))
              status(i, j) = KNOWN

            else if (lfn_out(i, j) * lfn_out(i - 1, j) < 0.0) then
              lfn_s0(i, j) = abs (lfn_out(i, j)) * dx / (abs (lfn_out(i, j)) + abs (lfn_out(i - 1, j)))
              status(i, j) = KNOWN

            else if (lfn_out(i, j) * lfn_out(i, j + 1) < 0.0) then
              lfn_s0(i, j) = abs (lfn_out(i, j)) * dy / (abs (lfn_out(i, j)) + abs (lfn_out(i, j + 1)))
              status(i, j) = KNOWN

            else if (lfn_out(i, j) * lfn_out(i, j - 1) < 0.0) then
              lfn_s0(i, j) = abs (lfn_out(i, j)) * dy / (abs (lfn_out(i, j)) + abs (lfn_out(i, j - 1)))
              status(i, j) = KNOWN
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine Set_ls_interface_subgrid_inter

    subroutine Set_ls_interface_quadratic_inter (lfn_s0, status, lfn_out, i_start, i_end, &
        j_start, j_end, num_tiles, ifms, ifme, jfms, jfme, dx, dy)

      ! Use 3 points for a quadratic interpolation: lfn(x) = a x^2 + b x + c

      implicit none

      integer, intent(in) :: ifms, ifme, jfms, jfme, num_tiles
      integer, dimension(num_tiles), intent(in) :: i_start, i_end, j_start, j_end
      real, dimension(ifms:ifme, jfms:jfme), intent(in) :: lfn_out
      real, dimension(ifms:ifme, jfms:jfme), intent(in out) :: lfn_s0
      integer, dimension(ifms:ifme, jfms:jfme), intent(in out) :: status
      real, intent (in) :: dx, dy

      integer :: ifts, ifte, jfts, jfte, i, j, ij
      real :: denom, ratio, a, b, c, x0, x1, x2, y0, y1, y2, root
      real, parameter :: TOL = 1.0e-12


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
              ! x-direction
            if (lfn_out(i, j) * lfn_out (i + 1, j) < 0.0) then
              x0 = -dx
              x1 = 0.0
              x2 = dx
              y0 = lfn_out(i - 1, j)
              y1 = lfn_out(i, j)
              y2 = lfn_out(i + 1, j)

              denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
              if (abs (denom) > TOL) then
                a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
                b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom
                c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / denom

                  ! Solve lfn(x) = 0 for x between [0, dx]
                root = (-b + sign (1.0, lfn_out(i + 1, j)) * sqrt (b * b - 4 * a * c)) / (2 * a)
                root = max (min (root, dx), 0.0)
                lfn_s0(i, j) = abs (root)
                status(i, j) = KNOWN
              end if

              ! y-direction
            else if (lfn_out(i, j) * lfn_out(i, j + 1) < 0.0) then
              y0 = -dy
              y1 = 0.0
              y2 = dy
              x0 = lfn_out(i, j - 1)
              x1 = lfn_out(i, j)
              x2 = lfn_out(i, j + 1)

              denom = (y0 - y1) * (y0 - y2) * (y1 - y2)
              if (abs(denom) > TOL) then
                a = (y2 * (x1 - x0) + y1 * (x0 - x2) + y0 * (x2 - x1)) / denom
                b = (y2 * y2 * (x0 - x1) + y1 * y1* (x2 - x0) + y0 * y0 * (x1 - x2)) / denom
                c = (y1 * y2 * (y1 - y2) * x0 + y2 * y0 * (y2 - y0 )* x1 + y0 * y1 * (y0 - y1) * x2) / denom

                root = (-b + sign(1.0, lfn_out(i, j + 1)) * sqrt (b * b - 4 * a * c)) / (2 * a)
                root = max (min (root, dy), 0.0)
                lfn_s0(i, j) = abs (root)
                status(i, j) = KNOWN
              end if
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine Set_ls_interface_quadratic_inter

    subroutine Set_ls_interface_gradient_inter (lfn_s0, status, lfn_out, i_start, i_end, &
        j_start, j_end, num_tiles, ifms, ifme, jfms, jfme, dx, dy)

      implicit none

      integer, intent(in) :: ifms, ifme, jfms, jfme, num_tiles
      integer, dimension(num_tiles), intent(in) :: i_start, i_end, j_start, j_end
      real, dimension(ifms:ifme, jfms:jfme), intent(in) :: lfn_out
      real, dimension(ifms:ifme, jfms:jfme), intent(in out) :: lfn_s0
      integer, dimension(ifms:ifme, jfms:jfme), intent(in out) :: status
      real, intent (in) :: dx, dy

      integer :: ifts, ifte, jfts, jfte, i, j, ij
      real dphidx, dphidy, gradphi


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij, i, j, ifts, ifte, jfts, jfte)
      do ij = 1, num_tiles
        ifts = i_start(ij)
        ifte = i_end(ij)
        jfts = j_start(ij)
        jfte = j_end(ij)

        do j = jfts, jfte
          do i = ifts, ifte
            if (lfn_out(i, j) == 0.0) then
              lfn_s0(i, j) = 0.0
              status(i, j) = KNOWN
            else if (lfn_out(i, j) * lfn_out(i + 1, j) < 0.0 .or. &
                lfn_out(i,j) * lfn_out(i - 1, j) < 0.0 .or. &
                lfn_out(i,j) * lfn_out(i, j + 1) < 0.0 .or. &
                lfn_out(i,j) * lfn_out(i, j - 1) < 0.0) then

               ! Calc grads
             dphidx = (lfn_out(i + 1, j) - lfn_out(i - 1, j)) / (2.0 * dx)
             dphidy = (lfn_out(i, j + 1) - lfn_out(i, j - 1)) / (2.0 * dy)
             gradphi = sqrt (dphidx ** 2 + dphidy ** 2 + 1.0e-12)

              lfn_s0(i, j) = abs (lfn_out(i, j)) / gradphi
              status(i, j) = KNOWN
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine Set_ls_interface_gradient_inter

    subroutine Solve_eikonal_eq_homo_dxys (i, j, lfn, dx, dy, ifms, ifme, jfms, jfme, val)

      implicit none

      integer, intent(in) :: i, j, ifms, ifme, jfms, jfme
      real, intent(in) :: dx, dy
      real, dimension(ifms:ifme, jfms:jfme), intent(in) :: lfn
      real, intent(out) :: val

      real :: ax, ay, tmp, h


      if (dx /= dy) call Stop_simulation ('This eikonal sover only works for homogenous grids (dx = dy)')

      ax = min (lfn(i - 1, j), lfn(i + 1, j))
      ay = min (lfn(i, j - 1), lfn(i, j + 1))

      tmp = abs(ax - ay)
      h = dx
      if (tmp >= h) then
        val = min(ax, ay) + h
      else
        val = (ax + ay + sqrt (2.0 * h ** 2 - tmp ** 2)) / 2.0
      end if

    end subroutine Solve_eikonal_eq_homo_dxys

    subroutine Reinit_ls_fsm (lfn, status, i_start, i_end, j_start, j_end, num_tiles, dx, dy, &
        ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, cart_comm)

    ! Purpose: Fast Sweeping Method. It is thread safe and there are no race conditions. 
    !          It will not give bit for bit results but will converge to the right solution.

#ifdef DM_PARALLEL
      use mpi
#endif

      implicit none

      integer, intent (in) :: num_tiles, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, cart_comm
      integer, dimension(ifms:ifme, jfms:jfme), intent(in) :: status
      real, dimension(ifms:ifme, jfms:jfme), intent(in out) :: lfn
      integer, dimension (num_tiles), intent (in) :: i_start, i_end, j_start, j_end
      real, intent(in) :: dx, dy
     
      real, parameter ::  TOL = 1.0e-6
      integer, parameter :: MAX_ITER = 20

      integer :: iter, i, j, ij, ifts, ifte, jfts, jfte, ierr
      real :: old, new_val, diff, max_diff, global_max_diff

      logical, parameter :: DEBUG_LOCAL = .false.

      do iter = 1, MAX_ITER
        max_diff = 0.0

          ! Sweep: i=1->nx, j=1->ny
        !$omp parallel do default(shared) private(i, j, old, new_val, diff, ifts, ifte, jfts, jfte, ij) &
        !$omp reduction(max:max_diff)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          do i = ifts, ifte
            do j = jfts, jfte
              if (status(i, j) /= KNOWN) then
                old = lfn(i,j)
                call Solve_eikonal_eq_homo_dxys (i, j, lfn, dx, dy, ifms, ifme, jfms, jfme, new_val)
                lfn(i,j) = min(old,new_val)
                diff = abs(lfn(i,j) - old)
                if (diff > max_diff) max_diff = diff
              end if
            end do
          end do
        end do
        !$omp end parallel do

#ifdef DM_PARALLEL
        call Do_halo_exchange (lfn, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

          ! Sweep: i=nx->1, j=1->ny
        !$omp parallel do default(shared) private(i, j, old, new_val, diff, ifts, ifte, jfts, jfte, ij) &
        !$omp reduction(max:max_diff)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          do i = ifte, ifts, -1
            do j = jfts, jfte
              if (status(i, j) /= KNOWN) then
                old = lfn(i,j)
                call Solve_eikonal_eq_homo_dxys (i, j, lfn, dx, dy, ifms, ifme, jfms, jfme, new_val)
                lfn(i,j) = min(old,new_val)
                diff = abs(lfn(i,j) - old)
                if (diff > max_diff) max_diff = diff
              end if
            end do
          end do
        end do
        !$omp end parallel do

#ifdef DM_PARALLEL
        call Do_halo_exchange (lfn, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

          ! Sweep: i=1->nx, j=ny->1
        !$omp parallel do default(shared) private(i, j, old, new_val, diff, ifts, ifte, jfts, jfte, ij) &
        !$omp reduction(max:max_diff)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          do i = ifts, ifte
            do j = jfte, jfts, -1
              if (status(i,j) /= KNOWN) then
                old = lfn(i,j)
                call Solve_eikonal_eq_homo_dxys (i, j, lfn, dx, dy, ifms, ifme, jfms, jfme, new_val)
                lfn(i,j) = min(old,new_val)
                diff = abs(lfn(i,j) - old)
                if (diff > max_diff) max_diff = diff
              end if
            end do
          end do
        end do
        !$omp end parallel do

#ifdef DM_PARALLEL
        call Do_halo_exchange (lfn, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

          ! Sweep: i=nx->1, j=ny->1
        !$omp parallel do default(shared) private(i, j, old, new_val, diff, ifts, ifte, jfts, jfte, ij) &
        !$omp reduction(max:max_diff)
        do ij = 1, num_tiles
          ifts = i_start(ij)
          ifte = i_end(ij)
          jfts = j_start(ij)
          jfte = j_end(ij)

          do i = ifte, ifts, -1
            do j = jfte, jfts, -1
              if (status(i,j) /= KNOWN) then
                old = lfn(i,j)
                call Solve_eikonal_eq_homo_dxys (i, j, lfn, dx, dy, ifms, ifme, jfms, jfme, new_val)
                lfn(i,j) = min(old,new_val)
                diff = abs(lfn(i,j) - old)
                if (diff > max_diff) max_diff = diff
              end if
            end do
          end do
        end do
        !$omp end parallel do

#ifdef DM_PARALLEL
        call Do_halo_exchange (lfn, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, N_POINTS_IN_HALO, cart_comm)
#endif

#ifdef DM_PARALLEL
        call MPI_Allreduce(max_diff, global_max_diff, 1, MPI_REAL, MPI_MAX, cart_comm, ierr)
#else
        global_max_diff = max_diff
#endif

        if (global_max_diff < TOL) exit
      end do

!      if (DEBUG_LOCAL) print *, 'FSM number iterations = ', iter

    end subroutine Reinit_ls_fsm

    pure function Select_upwind (diff_lx, diff_rx) result (return_value)

      ! upwind differences, L or R if both same sign, otherwise zero

      implicit none

      real, intent (in) :: diff_lx, diff_rx
      real :: return_value

      real :: diff2x


      diff2x = 0.0
      if (diff_lx > 0.0 .and. diff_rx > 0.0) diff2x = diff_lx
      if (diff_lx < 0.0 .and. diff_rx < 0.0) diff2x = diff_rx

      return_value = diff2x

    end function Select_upwind

    pure function Select_godunov (diff_lx, diff_rx) result (return_value)

      ! Godunov scheme: upwind differences, L or R or none    
      ! always test on > or < never = , much faster because of IEEE
      ! central diff >= 0 => take left diff if >0, ortherwise 0
      ! central diff <= 0 => take right diff if <0, ortherwise 0

      implicit none

      real, intent (in) :: diff_lx, diff_rx
      real :: return_value

      real :: diff2x, diff_cx


      diff2x = 0.0
      diff_cx = diff_rx + diff_lx
      if (diff_lx > 0.0 .and. .not. diff_cx < 0.0) diff2x = diff_lx
      if (diff_rx < 0.0 .and.       diff_cx < 0.0) diff2x = diff_rx

      return_value = diff2x

    end function Select_godunov

    pure function Select_eno (diff_lx, diff_rx) result (return_value)

      ! 1st order ENO scheme: choose the smaller magnitude derivative
      ! Both positive: pick smaller
      ! Both negative: pick larger (less negative)
      ! Mixed signs: zero (derivative crosses zero)

      implicit none

      real, intent(in) :: diff_lx, diff_rx
      real :: return_value


      if (diff_lx * diff_rx > 0.0) then
        return_value = sign(1.0, diff_lx) * min(abs(diff_lx), abs(diff_rx))
      else
        return_value = 0.0
      end if

    end function Select_eno

    pure function Select_2nd (dx, lfn_i, lfn_im1, lfn_ip1) result (return_value)

      ! 2nd-order advection scheme in the x,y-direction (DME)

      implicit none

      real, intent(in):: lfn_i, lfn_im1, lfn_ip1
      real, intent(in):: dx
      real :: return_value

      real  :: diff2x_p, diff2x_m


      diff2x_p = 0.0
      diff2x_m = 0.0
      diff2x_p = (lfn_ip1 + lfn_i) / (2.0 * dx)
      diff2x_m = (lfn_i + lfn_im1) / (2.0 * dx)
      return_value = diff2x_p - diff2x_m

    end function Select_2nd

    pure function Select_4th (dx, lfn_i, lfn_im1, lfn_im2, lfn_ip1, lfn_ip2) &
        result (return_value)

      ! 4th-order advection scheme in the x,y-direction (DME)

      implicit none

      real, intent(in) :: dx, lfn_i, lfn_im1, lfn_im2, lfn_ip1, lfn_ip2
      real :: return_value

      real :: diff2x_p, diff2x_m


      diff2x_p = 0.0
      diff2x_m = 0.0
      diff2x_p = (7.0 * lfn_ip1 + 7.0 * lfn_i - lfn_ip2 - lfn_im1) / (12.0 * dx)
      diff2x_m = (7.0 * lfn_i + 7.0 * lfn_im1 - lfn_ip1 - lfn_im2) / (12.0 * dx)
      return_value = diff2x_p - diff2x_m

    end function Select_4th

    pure function Select_weno3 (dx, lfn_it, lfn_im1t, lfn_im2t, lfn_ip1t, &
        lfn_ip2t, uf) result (return_value)

      ! 3rd-order advection WENO scheme in the x,y-direction (DME)

      implicit none

      real, intent(in) :: dx, lfn_it, lfn_im1t, lfn_im2t, lfn_ip1t, lfn_ip2t, uf
      real :: return_value

      real, parameter :: GAMMA1 = 1.0 / 3.0, GAMMA2 = 2.0 / 3.0, TOL = 1e-6
      real :: lfn_i, lfn_im1, lfn_im2, lfn_ip1, lfn_ip2, flux_p,flux_m, &
          w1, w2, w1t, w2t, beta1, beta2, fh_1, fh_2


      if (uf >= 0.0) then
        lfn_i = lfn_it
        lfn_im1 = lfn_im1t
        lfn_im2 = lfn_im2t
        lfn_ip1 = lfn_ip1t
      else
        lfn_i = lfn_it
        lfn_im1 = lfn_ip1t
        lfn_im2 = lfn_ip2t
        lfn_ip1 = lfn_im1t
      end if

        ! numerical flux at i,j+1/2 face
      fh_1 = -0.5 * lfn_im1 + 1.5 * lfn_i
      fh_2 = 0.5 * lfn_i + 0.5 * lfn_ip1
      beta1 = (lfn_i - lfn_im1) ** 2
      beta2 = (lfn_ip1 - lfn_i) ** 2
      w1t = GAMMA1 / (beta1 + TOL) ** 2
      w2t = GAMMA2 / (beta2 + TOL) ** 2
      w1 = w1t / (w1t + w2t)
      w2 = w2t / (w1t + w2t)
      flux_p = w1 * fh_1 + w2 * fh_2
        ! numerical flux at i,j-1/2 face
      fh_1 = -0.5 * lfn_im2 + 1.5 * lfn_im1
      fh_2 = 0.5 * lfn_im1 + 0.5 * lfn_i
      beta1 = (lfn_im1 - lfn_im2) ** 2
      beta2 = (lfn_i - lfn_im1) ** 2
      w1t = GAMMA1 / (beta1 + TOL) ** 2
      w2t = GAMMA2 / (beta2 + TOL) ** 2
      w1 = w1t / (w1t + w2t)
      w2 = w2t / (w1t + w2t)
      flux_m = w1 * fh_1 + w2 * fh_2

      if (uf >= 0.0 ) then
        return_value = (flux_p - flux_m) / dx
      else
        return_value = (flux_m - flux_p) / dx
      end if

    end function Select_weno3

    pure function Select_weno5 (dx, lfn_it, lfn_im1t, lfn_im2t, lfn_im3t, lfn_ip1t, &
        lfn_ip2t, lfn_ip3t, uf) result (return_value)

      ! 5th-order advection WENO scheme in the x,y-direction (DME)

      implicit none

      real, intent(in):: dx, lfn_it, lfn_im1t, lfn_im2t, lfn_im3t, lfn_ip1t, lfn_ip2t, lfn_ip3t, uf
      real :: return_value

      real, parameter:: GAMMA1 = 1.0 / 10.0, GAMMA2 = 3.0 / 5.0, GAMMA3 = 3.0 / 10.0, TOL = 1e-6
      real :: lfn_i, lfn_im1, lfn_im2, lfn_im3, lfn_ip1, lfn_ip2, lfn_ip3,  flux_p, flux_m , &
          w1, w2, w3, w1t, w2t, w3t, beta1, beta2, beta3, fh_1, fh_2, fh_3


      if (uf >= 0.0) then
        lfn_i = lfn_it
        lfn_im1 = lfn_im1t
        lfn_im2 = lfn_im2t
        lfn_im3 = lfn_im3t
        lfn_ip1 = lfn_ip1t
        lfn_ip2 = lfn_ip2t
      else
        lfn_i = lfn_it
        lfn_im1 = lfn_ip1t
        lfn_im2 = lfn_ip2t
        lfn_im3 = lfn_ip3t
        lfn_ip1 = lfn_im1t
        lfn_ip2 = lfn_im2t
      end if

        ! numerical flux at i,j+1/2 face
      fh_1 = (2.0 * lfn_im2 - 7.0 * lfn_im1 + 11.0 * lfn_i) / 6.0
      fh_2 = (-1.0 * lfn_im1 + 5.0 * lfn_i + 2.0 * lfn_ip1) / 6.0
      fh_3 = (2.0 * lfn_i + 5.0 * lfn_ip1 - 1.0 * lfn_ip2) / 6.0
      beta1 = (13.0 / 12.0) * (lfn_im2 - 2.0 * lfn_im1 + lfn_i) ** 2 + &
          0.25 * (lfn_im2 - 4.0 * lfn_im1 + 3.0 * lfn_i) ** 2
      beta2 = (13.0 / 12.0) * (lfn_im1 - 2.0 * lfn_i + lfn_ip1) ** 2 + &
          0.25 * (lfn_im1 - lfn_ip1) ** 2
      beta3 = (13.0 / 12.0) * (lfn_i - 2.0 * lfn_ip1 + lfn_ip2) ** 2 + &
          0.25 * (3.0 * lfn_i - 4.0 * lfn_ip1 + lfn_ip2) ** 2
      w1t = GAMMA1 / (beta1 + TOL) ** 2
      w2t = GAMMA2 / (beta2 + TOL) ** 2
      w3t = GAMMA3 / (beta3 + TOL) ** 2
      w1 = w1t / (w1t + w2t + w3t)
      w2 = w2t / (w1t + w2t + w3t)
      w3 = w3t / (w1t + w2t + w3t)
      flux_p = w1 * fh_1 + w2 * fh_2 + w3 * fh_3

        ! numerical flux at i,j-1/2 face
      fh_1 = (2.0 * lfn_im3 - 7.0 * lfn_im2 + 11.0 * lfn_im1) / 6.0
      fh_2 = (-1.0 * lfn_im2 + 5.0 * lfn_im1 + 2.0 * lfn_i) / 6.0
      fh_3 = (2.0 * lfn_im1 + 5.0 * lfn_i - 1.0 * lfn_ip1) / 6.0
      beta1 = (13.0 / 12.0) * (lfn_im3 - 2.0 * lfn_im2 + lfn_im1) ** 2 + &
          0.25 * (lfn_im3 - 4.0 * lfn_im2 + 3.0 * lfn_im1) ** 2
      beta2 = (13.0 / 12.0) * (lfn_im2 - 2.0 * lfn_im1 + lfn_i) ** 2 + &
          0.25 * (lfn_im2 - lfn_i) ** 2
      beta3 = (13.0 / 12.0) * (lfn_im1 - 2.0 * lfn_i + lfn_ip1) ** 2 + &
          0.25 * (3.0 * lfn_im1 - 4.0 * lfn_i + lfn_ip1) ** 2
      w1t = GAMMA1 / (beta1 + TOL) ** 2
      w2t = GAMMA2 / (beta2 + TOL) ** 2
      w3t = GAMMA3 / (beta3 + TOL) ** 2
      w1 = w1t / (w1t + w2t + w3t)
      w2 = w2t / (w1t + w2t + w3t)
      w3 = w3t / (w1t + w2t + w3t)
      flux_m = w1 * fh_1 + w2 * fh_2 + w3 * fh_3

      if (uf >= 0.0) then
        return_value = (flux_p - flux_m) / dx
      else
        return_value = (flux_m - flux_p) / dx
      end if

    end function Select_weno5

    subroutine Stop_if_close_to_bdy (ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, &
        ifds, jfds, ifde, jfde, lfn_out)

      implicit none

      integer, intent (in) :: ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme, &
          ifds, jfds, ifde, jfde
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: lfn_out

      integer, parameter :: BOUNDARY_GUARD = 8
      integer :: i, j


      do j = jfts, jfte
        if (j <= BOUNDARY_GUARD .or. j > (jfde - BOUNDARY_GUARD)) then
          do i = ifts, ifte
            if (lfn_out(i, j) < 0.0) call Stop_simulation ('Fire too close to domain boundary.')
          end do
        end if
      end do

      do i = ifts, ifte
        if (i <= BOUNDARY_GUARD .or. i > (ifde - BOUNDARY_GUARD)) then
          do j = jfts, jfte
            if (lfn_out(i, j) < 0.0) call Stop_simulation ('Fire too close to domain boundary.')
          end do
        end if
      end do

    end subroutine Stop_if_close_to_bdy

    subroutine Extrapol_var_at_bdys (ifms, ifme, jfms, jfme, ifds, ifde, jfds, jfde, &
        ifts, ifte, jfts, jfte, var)

      ! extend 2D array beyond domain boundaries

      implicit none

      integer, intent (in) :: ifms, ifme, jfms, jfme, ifds, ifde, jfds, jfde, ifts, ifte, jfts, jfte
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: var

      integer :: i, j
      integer :: ifts1, ifte1, jfts1, jfte1

        ! only 1 is needed since ENO1 is 
        ! used near domain boundaries
      integer, parameter :: HALO = 1


        ! Go HALO width beyond if not a tile in a domain corner
        ! assume we have halo need to compute the value we do not have
      ifts1 = ifts
      jfts1 = jfts
      ifte1 = ifte
      jfte1 = jfte
      if (.not. ifts == ifds) ifts1 = ifts - HALO
      if (.not. jfts == jfds) jfts1 = jfts - HALO
      if (.not. ifte == ifde) ifte1 = ifte + HALO
      if (.not. jfte == jfde) jfte1 = jfte + HALO

      if (ifts == ifds)then
        do j = jfts1, jfte1
          var(ifds - 1, j) = Extrapol (var(ifds, j), var(ifds + 1, j))
        end do
      end if

      if (ifte == ifde) then
        do j = jfts1, jfte1
          var(ifde + 1, j) = Extrapol (var(ifde, j), var(ifde - 1, j))
        end do
      end if

      if (jfts == jfds) then
        do i = ifts1, ifte1
          var(i, jfds - 1) = Extrapol (var(i, jfds), var(i, jfds + 1))
        end do
      end if

      if (jfte == jfde) then
        do i = ifts1, ifte1
          var(i, jfde + 1) = Extrapol (var(i, jfde), var(i, jfde - 1))
        end do
      end if

        ! corners
      if (ifts == ifds .and. jfts == jfds) var(ifds - 1, jfds - 1) = Extrapol (var(ifds, jfds), var(ifds + 1, jfds + 1))
      if (ifts == ifds .and. jfte == jfde) var(ifds - 1, jfde + 1) = Extrapol (var(ifds, jfde), var(ifds + 1, jfde - 1))
      if (ifte == ifde .and. jfts == jfds) var(ifde + 1, jfds - 1) = Extrapol (var(ifde, jfds), var(ifde - 1, jfds + 1))
      if (ifte == ifde .and. jfte == jfde) var(ifde + 1, jfde + 1) = Extrapol (var(ifde, jfde), var(ifde - 1, jfde - 1))

    contains

      pure function Extrapol (a, b) result (return_value)

        implicit none

        real, intent (in) :: a, b

        real :: return_value
        real, parameter :: BIAS = 1.0

          ! extrapolation, max quarded
        return_value = (1.0 - BIAS) * (2.0 * a - b) + BIAS * max (2.0 * a - b, a, b)

      end function Extrapol

    end subroutine Extrapol_var_at_bdys

  end module level_set_mod

