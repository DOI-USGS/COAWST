  module fmc_mod

    use fuel_mod, only : fuel_t
    use ros_mod, only : ros_t

    implicit none

    private
    
    public :: fmc_t, FMC_WRFFIRE

    integer, parameter :: FMC_WRFFIRE = -1

    type, abstract :: fmc_t
    contains
      procedure (Advance_fmc_model), deferred :: Advance_fmc_model
      procedure (Init), deferred :: Init
    end type fmc_t

    abstract interface
      subroutine Init (this, fuels, fuelmc_g, fuelmc_g_live, ifms, ifme, jfms, jfme, itimestep, dt)
        import :: fmc_t, fuel_t
        class (fmc_t), intent (in out) :: this
        class (fuel_t), intent (in) :: fuels
        real, intent (in) :: fuelmc_g, fuelmc_g_live, dt
        integer, intent (in) :: ifms, ifme, jfms, jfme, itimestep
      end subroutine Init

      subroutine Advance_fmc_model (this, fmoist_freq, fmoist_dt, itimestep, dt, ifms, ifme, jfms, jfme, i_start, i_end, &
            j_start, j_end, num_tiles, fire_rain, fire_t2, fire_q2, fire_psfc, fire_rain_old, fire_t2_old, fire_q2_old, &
            fire_psfc_old, fire_rh_fire, fuelmc_g, fmc_g, nfuel_cat, fuels, ros_param)
        import :: fmc_t, fuel_t, ros_t
        class (fmc_t), intent (in out) :: this
        class (fuel_t), intent (in) :: fuels
        class (ros_t), intent (in out) :: ros_param
        integer, intent (in) :: fmoist_freq, itimestep, ifms, ifme, jfms, jfme, num_tiles
        integer, dimension(num_tiles) :: i_start, i_end, j_start, j_end
        real, intent (in) ::  fmoist_dt, dt, fuelmc_g
        real, dimension (ifms:ifme, jfms:jfme), intent (in) :: nfuel_cat
        real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: fire_rain, fire_t2, fire_q2, fire_psfc, fire_rain_old, &
            fire_t2_old, fire_q2_old, fire_psfc_old, fire_rh_fire, fmc_g
      end subroutine Advance_fmc_model
    end interface

  end module fmc_mod
