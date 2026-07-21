  module fuel_anderson_mod

    ! =============================================================================
    ! Anderson 13 surface fire fuel models, along with some
    !          estimated canopy properties (for crown fire).
    ! =============================================================================
    !  --- Grass-dominated fuel models
    !  FUEL MODEL 1: Short grass (1 ft)
    !  FUEL MODEL 2: Timber (grass and understory)
    !  FUEL MODEL 3: Tall grass (2.5 ft)
    !  --- Shrub-dominated fuel models
    !  FUEL MODEL 4: Chaparral (6 ft)
    !  FUEL MODEL 5: Brush (2 ft)
    !  FUEL MODEL 6: Dormant brush, hardwood slash
    !  FUEL MODEL 7: Southern rough
    !  --- Forest litter-dominated fuel models
    !  FUEL MODEL 8: Closed timber litter
    !  FUEL MODEL 9: Hardwood litter
    !  FUEL MODEL 10: Timber (litter + understory)
    !  --- Logging debris-dominated fuel models
    !  FUEL MODEL 11: Light logging slash
    !  FUEL MODEL 12: Medium logging slash
    !  FUEL MODEL 13: Heavy logging slash
    !  --- Fuel-free
    !  FUEL MODEL 14: no fuel
    ! =============================================================================

    use fuel_mod, only : fuel_t

    implicit none

    private

    public :: fuel_anderson_t

    integer, parameter :: N_FUEL_CAT_ANDERSON = 13, NO_FUEL_CAT_ANDERSON = 14

    type, extends (fuel_t) :: fuel_anderson_t
      character (len = 80), dimension(N_FUEL_CAT_ANDERSON + 1) :: fuel_name
        ! Canopy
        ! fct: burn out time for canopy fuel, after dry [s]
      real, dimension(N_FUEL_CAT_ANDERSON + 1) :: fct = [ 60., 60., 60., 60., 60., 60., 60., &
                                                          60., 120., 180., 180., 180., 180. , 60. ]
        ! fci_d: initial dry mass of canopy fuel
        ! ----- 1.12083 is 5 tons/acre.  5-50 tons/acre orig., 100-300 after blowdown
      real, dimension(N_FUEL_CAT_ANDERSON + 1) :: fci_d = [ 0., 0., 0., 1.123, 0., 0., 0., &
                                                            1.121, 1.121, 1.121, 1.121, 1.121, 1.121, 0. ]
      real, dimension(N_FUEL_CAT_ANDERSON + 1) :: windrf = [ 0.36, 0.36, 0.44,  0.55,  0.42,  0.44,  0.44,  &
                                                             0.36, 0.36, 0.36,  0.36,  0.43,  0.46,  1.0e-7 ]
        ! fci: initial total mass of canopy fuel
        ! fcbr: fuel canopy burn rate (kg/m2/s)
      real, dimension(N_FUEL_CAT_ANDERSON + 1) :: fcbr, fci
    contains
       procedure, public :: Initialization => Init_anderson_fuel_model
    end type fuel_anderson_t

  contains

   subroutine Init_anderson_fuel_model (this, fuelmc_c)

      implicit none

      class (fuel_anderson_t), intent(in out) :: this
      real, intent (in) :: fuelmc_c

      integer :: i


      this%n_fuel_cat = N_FUEL_CAT_ANDERSON
      this%no_fuel_cat = NO_FUEL_CAT_ANDERSON

      this%fgi = [ 0.166, 0.896, 0.674, 3.591, 0.784, 1.344, 1.091, 1.120, 0.780, 2.692, 2.582, 7.749, 13.024, 1.e-7 ]
      this%fueldepthm = [ 0.305,  0.305,  0.762, 1.829, 0.61, 0.762, 0.762, 0.0610, 0.0610, 0.305, 0.305, 0.701, 0.914, 0.305 ]
      this%weight = [ 7.0,  7.0,  7.0, 180.0, 100.0, 100.0, 100.0, 900.0, 900.0, 900.0, 900.0, 900.0, 900.0, 7.0 ]
      this%ichap = [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
      allocate (this%fueldens(N_FUEL_CAT_ANDERSON + 1))
      this%fueldens = 32.0
      allocate (this%st(N_FUEL_CAT_ANDERSON + 1))
      this%st = 0.0555
      allocate (this%se(N_FUEL_CAT_ANDERSON + 1))
      this%se = 0.010
      this%savr = [ 3500.0, 2784.0, 1500.0, 1739.0, 1683.0, 1564.0, 1562.0, 1889.0, 2484.0, 1764.0, 1182.0, 1145.0, 1159.0, 3500.0 ]
      this%fuelmce = [ 0.12, 0.15, 0.25, 0.20, 0.20, 0.25, 0.40, 0.30, 0.25, 0.25, 0.15, 0.20, 0.25, 0.12 ]

        ! following Albini 1976 as reprinted in Anderson 1982 Table 1 (for proportions only)
        ! to convert ton/acre to 1 km/m2 multiply by 4.4609
      this%fgi_1h =  [ 0.74, 2.00, 3.01, 5.01, 1.00, 1.50, 1.13, 1.50, 2.92, 3.01, 1.50, 4.01, 7.01, 0.0 ]
      this%fgi_10h = [ 0.00, 1.00, 0.00, 4.01, 0.50, 2.50, 1.87, 1.00, 0.41, 2.00, 4.51, 14.03, 23.04, 0.0 ]
      this%fgi_100h = [ 0.00, 0.50, 0.00, 2.00, 0.00, 2.00, 1.50, 2.50, 0.15, 5.01, 5.51, 16.53, 28.05, 0.0 ]
      allocate (this%fgi_1000h(N_FUEL_CAT_ANDERSON + 1))
      this%fgi_1000h = 0.0
      this%fgi_live = [ 0.00, 0.50, 0.000, 5.01, 2.00, 0.00, 0.37, 0.00, 0.00, 2.00, 0.00, 0.0, 0.00, 0.0 ]

      this%fuel_name(1)  = '1: Short grass (1 ft)'
      this%fuel_name(2)  = '2: Timber (grass and understory)'
      this%fuel_name(3)  = '3: Tall grass (2.5 ft)'
      this%fuel_name(4)  = '4: Chaparral (6 ft)'
      this%fuel_name(5)  = '5: Brush (2 ft) '
      this%fuel_name(6)  = '6: Dormant brush, hardwood slash'
      this%fuel_name(7)  = '7: Southern rough'
      this%fuel_name(8)  = '8: Closed timber litter'
      this%fuel_name(9)  = '9: Hardwood litter'
      this%fuel_name(10) = '10: Timber (litter + understory)'
      this%fuel_name(11) = '11: Light logging slash'
      this%fuel_name(12) = '12: Medium logging slash'
      this%fuel_name(13) = '13: Heavy logging slash'
      this%fuel_name(14) = '14: no fuel'

      this%fci = 0.0
      this%fcbr = 0.0
      do i = 1, N_FUEL_CAT_ANDERSON
        this%fci(i) = (1.0 + fuelmc_c) * this%fci_d(i)
        if(this%fct(i) /= 0.0) then
          this%fcbr(i) = this%fci_d(i) / this%fct(i)
        else
          this%fcbr(i) = 0.0
        end if
      end do

      ! Initialize wind adjustment factor
      call this%Calc_wind_adjustment_factor()

    end subroutine Init_anderson_fuel_model

  end module fuel_anderson_mod

