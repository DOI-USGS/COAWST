  module advance_mod

    use state_mod, only : state_fire_t
    use namelist_mod, only : namelist_t
    use fire_driver_mod, only : Advance_fire_components
    use stderrout_mod, only: Print_message

    implicit none

    private

    public :: Advance_state

  contains

    subroutine Advance_state (grid, config_flags)

      implicit none

      type (state_fire_t), intent (in out) :: grid
      type (namelist_t), intent (in) :: config_flags

      logical, parameter :: DEBUG_LOCAL = .false.


      if (DEBUG_LOCAL) call Print_message ('  Entering Advance_state... ')

      grid%itimestep = grid%itimestep + 1

      call Advance_fire_components (grid, config_flags)

      grid%datetime_now = grid%datetime_start
      call grid%datetime_now%Add_seconds (grid%itimestep * grid%dt)
      call grid%datetime_now%Print_datetime ()

      if (DEBUG_LOCAL) call Print_message ('  Leaving Advance_state... ')

    end subroutine Advance_state

  end module advance_mod
