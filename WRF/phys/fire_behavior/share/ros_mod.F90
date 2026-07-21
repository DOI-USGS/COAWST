  module ros_mod

    use fuel_mod, only : fuel_t

    implicit none

    private
    
    public :: ros_t, ROS_WRFFIRE

    integer, parameter :: ROS_WRFFIRE = 0

    type, abstract :: ros_t
      real, dimension(:, :), allocatable :: iboros
    contains
      procedure (Calc_ros), deferred :: Calc_ros
      procedure (Init), deferred :: Init
      procedure (Set_params), deferred :: Set_params
    end type ros_t

    abstract interface
      subroutine Init (this, ifms, ifme, jfms, jfme)
        import :: ros_t
        class (ros_t), intent (in out) :: this
        integer, intent (in) :: ifms, ifme, jfms, jfme
      end subroutine Init

      subroutine Set_params (this, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, fuels, nfuel_cat, fmc_g)
        import :: ros_t, fuel_t
        class (ros_t), intent (in out) :: this
        integer, intent (in) :: ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte
        class (fuel_t), intent (in) :: fuels
        real, dimension(ifms:ifme, jfms:jfme), intent (in) :: nfuel_cat, fmc_g
      end subroutine Set_params

      pure function Calc_ros (this, ifms, ifme, jfms, jfme, i, j, nvx, nvy, uf, vf, dzdxf, dzdyf) result (return_value)
        import :: ros_t, fuel_t
        class (ros_t), intent (in) :: this
        real, intent (in) :: nvx, nvy, uf, vf, dzdxf, dzdyf
        integer, intent (in) :: ifms, ifme, jfms, jfme, i, j
        real :: return_value
      end function Calc_ros

    end interface

  end module ros_mod
