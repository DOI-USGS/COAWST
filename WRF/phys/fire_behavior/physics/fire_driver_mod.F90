  module fire_driver_mod

    use fire_model_mod, only : Advance_fire_model
    use level_set_mod, only : Extrapol_var_at_bdys
    use state_mod, only : state_fire_t
    use namelist_mod, only : namelist_t
    use stderrout_mod, only : Print_message, Stop_simulation

    use fuel_mod, only : FUEL_ANDERSON
    use fuel_anderson_mod, only : fuel_anderson_t

    use ros_mod, only : ROS_WRFFIRE
    use ros_wrffire_mod, only : ros_wrffire_t

    use fmc_mod, only : FMC_WRFFIRE
    use fmc_wrffire_mod, only : fmc_wrffire_t

#ifdef DM_PARALLEL
    use mpi_mod, only : Max_across_mpi_tasks, Sum_across_mpi_tasks
#endif

    implicit none

    private

    public :: Advance_fire_components, Init_fire_components

    integer, parameter:: REAL_SUM = 10, REAL_MAX = 20, RNRM_SUM = 30, RNRM_MAX = 40

  contains

    subroutine Init_fire_components (grid, config_flags)

      implicit none

      type (state_fire_t) :: grid
      type (namelist_t), intent(in) :: config_flags

      integer :: ij


        ! Ignitions lines
      call grid%Init_ignition_lines (config_flags)

        ! Fuel model
      select case (config_flags%fuel_opt)
        case (FUEL_ANDERSON)
          allocate (fuel_anderson_t::grid%fuels)

        case default
          call Stop_simulation ('The selected fuel_opt does not exist')
      end select
      call grid%fuels%Initialization (config_flags%fuelmc_c)
      call grid%Init_fuel_vars ()

        ! FMC model
      select case (config_flags%fmc_opt)
        case (FMC_WRFFIRE)
          allocate (fmc_wrffire_t::grid%fmc_param)

        case default
          call Stop_simulation ('The selected fmc_param does not exist')
      end select
      if (config_flags%fmoist_run) call grid%fmc_param%Init (grid%fuels, config_flags%fuelmc_g, config_flags%fuelmc_g_live, grid%ifms, &
          grid%ifme, grid%jfms, grid%jfme, grid%itimestep, grid%dt)

        ! Rate of spread parameterization
      select case (config_flags%ros_opt)
        case (ROS_WRFFIRE)
          allocate (ros_wrffire_t::grid%ros_param)

        case default
          call Stop_simulation ('The selected ros_opt does not exist')
      end select
      call grid%ros_param%Init (grid%ifms, grid%ifme, grid%jfms, grid%jfme)

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE (ij)
      do ij = 1, grid%num_tiles
        call Extrapol_var_at_bdys (grid%ifms, grid%ifme, grid%jfms, grid%jfme, grid%ifds, grid%ifde, &
            grid%jfds, grid%jfde, grid%i_start(ij), grid%i_end(ij), grid%j_start(ij), grid%j_end(ij), &
            grid%lfn)

        call Extrapol_var_at_bdys (grid%ifms, grid%ifme, grid%jfms, grid%jfme, grid%ifds, grid%ifde, &
            grid%jfds, grid%jfde, grid%i_start(ij), grid%i_end(ij), grid%j_start(ij), grid%j_end(ij), &
            grid%tign_g)

        call grid%ros_param%Set_params (grid%ifms, grid%ifme, grid%jfms, grid%jfme, grid%i_start(ij), grid%i_end(ij), &
            grid%j_start(ij), grid%j_end(ij), grid%fuels, grid%nfuel_cat, grid%fmc_g)
      end do
      !$OMP END PARALLEL DO

    end subroutine Init_fire_components

    subroutine Advance_fire_components (grid, config_flags)

      implicit none

      type (state_fire_t), intent (in out) :: grid
      type (namelist_t), intent (in) :: config_flags

      integer, parameter :: PRINT_LEVEL = 1
      integer :: ij
      logical, parameter :: DEBUG_LOCAL = .false.


      if (DEBUG_LOCAL) call Print_message ('Entering Advance_fire_components...') 

      if (config_flags%fmoist_run) call grid%fmc_param%Advance_fmc_model (config_flags%fmoist_freq, config_flags%fmoist_dt, &
          grid%itimestep, grid%dt, grid%ifms, grid%ifme, grid%jfms, grid%jfme, &
          grid%i_start, grid%i_end, grid%j_start, &
          grid%j_end, grid%num_tiles, grid%fire_rain, grid%fire_t2, grid%fire_q2, grid%fire_psfc, &
          grid%fire_rain_old, grid%fire_t2_old, grid%fire_q2_old, grid%fire_psfc_old, grid%fire_rh_fire, config_flags%fuelmc_g, &
          grid%fmc_g, grid%nfuel_cat, grid%fuels, grid%ros_param)

      call Advance_fire_model (config_flags, grid)

      if (config_flags%fire_print_msg >= PRINT_LEVEL) call Print_summary (config_flags, grid)

      if (DEBUG_LOCAL) call Print_message ('Leaving Advance_fire_components...') 

    end subroutine Advance_fire_components

    function Calc_domain_stats (fun, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, a, b, cart_comm) result (return_value)

      implicit none

      integer, intent (in) :: fun, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, cart_comm
      real, dimension(ifms:ifme, jfms:jfme), intent (in) :: a, b
      real :: return_value

      real :: lmax, lsum, gmax, gsum
      integer :: i, j


      if (fun == REAL_SUM) then
        lsum = 0.0
        do j = jfps, jfpe
          do i = ifps, ifpe
            lsum = lsum + a(i, j)
          end do
        end do
#ifdef DM_PARALLEL
        call Sum_across_mpi_tasks (lsum, cart_comm, gsum)
        return_value = gsum
#else
        return_value = lsum
#endif

      else if (fun == RNRM_SUM) then
        lsum = 0.0
        do j = jfps, jfpe
          do i = ifps, ifpe
            lsum = lsum + sqrt (a(i, j) * a(i, j) + b(i, j) * b(i, j))
          end do
        end do
#ifdef DM_PARALLEL
        call Sum_across_mpi_tasks (lsum, cart_comm, gsum)
        return_value = gsum
#else
        return_value = lsum
#endif

      else if (fun == REAL_MAX) then
        lmax = - huge (lmax)
        do j = jfps, jfpe
          do i = ifps, ifpe
            lmax = max (lmax, a(i, j))
          end do
        end do
#ifdef DM_PARALLEL
        call Max_across_mpi_tasks (lmax, cart_comm, gmax)
        return_value = gmax
#else
        return_value = lmax
#endif

      else if (fun == RNRM_MAX) then
        lmax = 0.0
        do j = jfps, jfpe
          do i = ifps, ifpe
            lmax = max (lmax, sqrt (a(i, j) * a(i, j) + b(i, j) * b(i, j)))
          end do
        end do
#ifdef DM_PARALLEL
        call Max_across_mpi_tasks (lmax, cart_comm, gmax)
        return_value = gmax
#else
        return_value = lmax
#endif

      else
        call Stop_simulation ('Value not supported for printing summary')

      end if

    end function Calc_domain_stats

    subroutine Print_summary (config_flags, grid)

      implicit none

      type (namelist_t), intent (in) :: config_flags
      type (state_fire_t), intent (in) :: grid

      real :: tfa, thf, mhf, tqf, mqf, aw, mw
      real :: time_start
      integer :: stat_lev = 1
      integer :: ifds, ifde, jfds, jfde, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe
      character (len = 256) :: msg


      ifds = grid%ifds
      ifde = grid%ifde
      jfds = grid%jfds
      jfde = grid%jfde

      ifms = grid%ifms
      ifme = grid%ifme
      jfms = grid%jfms
      jfme = grid%jfme

      ifps = grid%ifps
      ifpe = grid%ifpe
      jfps = grid%jfps
      jfpe = grid%jfpe

      time_start = grid%itimestep * grid%dt

      aw = Calc_domain_stats (RNRM_SUM, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, grid%uf,grid%vf, grid%cart_comm) / &
          ((ifde - ifds + 1) * (jfde - jfds + 1))
      mw = Calc_domain_stats (RNRM_MAX, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, grid%uf, grid%vf, grid%cart_comm)

      tfa = Calc_domain_stats (REAL_SUM, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, grid%fire_area, &
         grid%fire_area, grid%cart_comm) * grid%dx * grid%dy

      thf = Calc_domain_stats (REAL_SUM, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, grid%fgrnhfx, &
          grid%fgrnhfx, grid%cart_comm) * grid%dx * grid%dy

      mhf = Calc_domain_stats (REAL_MAX, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, grid%fgrnhfx, &
          grid%fgrnhfx, grid%cart_comm)

      tqf = Calc_domain_stats (REAL_SUM, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, grid%fgrnqfx, &
          grid%fgrnqfx, grid%cart_comm) * grid%dx * grid%dy

      mqf = Calc_domain_stats (REAL_MAX, ifms, ifme, jfms, jfme, ifps, ifpe, jfps, jfpe, grid%fgrnqfx, &
          grid%fgrnqfx, grid%cart_comm)

 91   format('Time ',f11.3,' s ',a,e12.3,1x,a)

      write (msg, 91) time_start, 'Average wind        ', aw, 'm/s'
      call Print_message (msg)

      write (msg, 91) time_start, 'Maximum wind        ', mw, 'm/s'
      call Print_message (msg)

      write (msg, 91) time_start, 'Fire area           ', tfa, 'm^2'
      call Print_message (msg)

      write (msg, 91) time_start, 'Heat output         ', thf, 'W'
      call Print_message (msg)

      write (msg, 91) time_start, 'Max heat flux       ', mhf, 'W/m^2'
      call Print_message (msg)

      write (msg, 91) time_start, 'Latent heat output  ', tqf, 'W'
      call Print_message (msg)

      write (msg, 91) time_start, 'Max latent heat flux', mqf, 'W/m^2'
      call Print_message (msg)

    end subroutine Print_summary

  end module fire_driver_mod
