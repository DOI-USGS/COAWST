  module ros_wrffire_mod

    use constants_mod, only : CMBCNST, CONVERT_J_PER_KG_TO_BTU_PER_POUND
    use fuel_mod, only : fuel_t
    use namelist_mod, only : namelist_t
    use ros_mod, only : ros_t
    use state_mod, only : state_fire_t

    implicit none

    private

    public :: ros_wrffire_t

    logical, parameter :: FIRE_GROWS_ONLY = .true.
    integer, parameter :: SLOPE_FACTOR = 1.0

      ! fuelheat: fuel particle low heat content [btu/lb]
    real, parameter :: FUELHEAT = CMBCNST * CONVERT_J_PER_KG_TO_BTU_PER_POUND
    integer, parameter :: FIRE_ADVECTION = 1 ! "0 = fire spread computed from normal wind speed/slope, 1 = fireline particle speed projected on normal" "0"

    type, extends(ros_t) :: ros_wrffire_t
      real, dimension(:, :), allocatable :: bbb, ischap, betafl, phiwc, r_0
    contains
      procedure, public :: Calc_ros => Calc_ros_wrffire
      procedure, public :: Init => Init_ros_wrffire
      procedure, public :: Set_params => Set_ros_parameters_wrffire
    end type ros_wrffire_t

  contains

    pure function Calc_ros_wrffire (this, ifms, ifme, jfms, jfme, i, j, nvx, nvy, uf, vf, dzdxf, dzdyf) result (return_value)

      implicit none

      ! m/s = (ft/min) * 0.3048 / 60.0 = (ft/min) * .00508
      ! ft/min = m/s * 2.2369 * 88.0 = m/s *  196.850

      class (ros_wrffire_t), intent (in) :: this
      integer, intent (in) :: ifms, ifme, jfms, jfme, i, j
      real, intent (in) :: nvx, nvy, uf, vf, dzdxf, dzdyf
      real :: return_value

      real :: speed, tanphi ! windspeed and slope in the direction normal to the fireline
      real :: umid, phis, phiw, spdms, umidm, excess, ros_back, cor_wind, cor_slope, ros_base, ros_wind, ros_slope
      real, parameter :: ROS_MAX = 6.0


      if (FIRE_ADVECTION /= 0) then
          ! wind speed is total speed 
        speed = sqrt (uf * uf + vf * vf) + tiny (speed)
          ! slope is total slope
        tanphi = sqrt (dzdxf * dzdxf + dzdyf * dzdyf) + tiny (tanphi)
          ! cos of wind and spread, if >0
        cor_wind =  max (0.0, (uf * nvx + vf * nvy) / speed)
          ! cos of slope and spread, if >0
        cor_slope = max (0.0, (dzdxf * nvx + dzdyf * nvy) / tanphi)
      else
          ! wind speed in spread direction
        speed = uf * nvx + vf * nvy
          ! slope in spread direction
        tanphi = dzdxf * nvx + dzdyf * nvy
        cor_wind = 1.0
        cor_slope = 1.0
      end if

      if (.not. this%ischap(i, j) > 0.0) then
          ! Rothermel
        spdms = max (speed, 0.0)
        umidm = min (spdms, 30.0)
        umid = umidm * 196.850 ! m/s to ft/min
        phiw = umid ** this%bbb(i, j) * this%phiwc(i, j)
        phis = 0.0
        if (tanphi > 0.0) phis = 5.275 * (this%betafl(i, j)) ** (-0.3) * tanphi ** 2
        ros_base = this%r_0(i, j) * 0.00508 ! ft/min to m/s
        ros_wind = ros_base * phiw
        ros_slope = ros_base * phis
      else
        spdms = max (speed, 0.0)
        ros_back = 0.03333    ! chaparral backing fire spread rate 0.033 m/s   ! param!
          ! spread rate, m/s
        ros_wind = 1.2974 * spdms ** 1.41
        ros_wind = max (ros_wind, ros_back)
        ros_slope = 0.0
        ros_base = 0.0
      end if

      ros_wind = ros_wind * cor_wind
      ros_slope = ros_slope * cor_slope

      return_value = min (ros_base + ros_wind + SLOPE_FACTOR * ros_slope, ROS_MAX)
      if (FIRE_GROWS_ONLY) return_value = max (return_value, 0.0)

    end function Calc_ros_wrffire

    subroutine Init_ros_wrffire (this, ifms, ifme, jfms, jfme)

      implicit none

      class (ros_wrffire_t), intent (in out) :: this
      integer, intent (in) :: ifms, ifme, jfms, jfme


      allocate (this%iboros(ifms:ifme, jfms:jfme))
 
      allocate (this%ischap(ifms:ifme, jfms:jfme))
      allocate (this%betafl(ifms:ifme, jfms:jfme))
      allocate (this%bbb(ifms:ifme, jfms:jfme))
      allocate (this%phiwc(ifms:ifme, jfms:jfme))
      allocate (this%r_0(ifms:ifme, jfms:jfme))

    end subroutine Init_ros_wrffire

    subroutine Set_ros_parameters_wrffire (this, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, &
        fuels, nfuel_cat, fmc_g)

      implicit none

      class (ros_wrffire_t), intent (in out) :: this
      integer, intent(in) :: ifts, ifte, jfts, jfte, ifms, ifme, jfms, jfme
      class (fuel_t), intent (in) :: fuels
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: nfuel_cat, fmc_g


      real ::  fuelload, fueldepth, rtemp1, rtemp2, qig, epsilon, rhob, wn, betaop, e, c, &
          xifr, etas, etam, a, gammax, gamma, ratio, ir, fuelloadm, bmst
      integer:: i, j, k, kk
      character (len = 128) :: msg


      Loop_j: do j = jfts, jfte
        Loop_i: do i = ifts, ifte
          k = int (nfuel_cat(i, j))
          if(k == fuels%no_fuel_cat) then
            this%ischap(i, j) = 0.0
              ! set to 1.0 to prevent grid%betafl(i,j)**(-0.3) to be Inf in fire_ros
            this%betafl(i, j) = 1.0
            this%bbb(i, j) = 1.0
            this%phiwc(i, j) = 0.0
            this%r_0(i, j) = 0.0
            this%iboros(i, j) = 0.0
          else
            this%ischap(i, j) = fuels%ichap(k)
              ! Settings of fire spread parameters from Rothermel
              ! No need to recalculate if FMC does not change
            bmst = fmc_g(i, j) / (1.0 + fmc_g(i, j))
              !  fuelload without moisture
            fuelloadm = (1.0 - bmst) * fuels%fgi(k)
            fuelload = fuelloadm * (0.3048) ** 2 * 2.205 ! to lb/ft^2
            fueldepth = fuels%fueldepthm(k) / 0.3048 ! to ft
              ! packing ratio
            this%betafl(i, j) = fuelload / (fueldepth * fuels%fueldens(k))
              ! optimum packing ratio
            betaop = 3.348 * fuels%savr(k) ** (-0.8189)
              ! heat of preignition, btu/lb
            qig = 250.0 + 1116.0 * fmc_g(i, j)
              ! effective heating number
            epsilon = exp (-138.0 / fuels%savr(k))
              ! ovendry bulk density, lb/ft^3
            rhob = fuelload/fueldepth

              ! const in wind coef
            c = 7.47 * exp (-0.133 * fuels%savr(k) ** 0.55)
            this%bbb(i,j) = 0.02526 * fuels%savr(k) ** 0.54
            e = 0.715 * exp (-3.59e-4 * fuels%savr(k))
            this%phiwc(i,j) = c * (this%betafl(i, j) / betaop) ** (-e)

            rtemp2 = fuels%savr(k) ** 1.5
              ! maximum rxn vel, 1/min
            gammax = rtemp2 / (495.0 + 0.0594 * rtemp2)
              ! coef for optimum rxn vel
            a = 1.0 / (4.774 * fuels%savr(k) ** 0.1 - 7.27)
            ratio = this%betafl(i,j)/betaop
              !optimum rxn vel, 1/min
            gamma = gammax * (ratio ** a) * exp(a * (1.0 - ratio))

             ! net fuel loading, lb/ft^2
            wn = fuelload/(1 + fuels%st(k))
            rtemp1 = fmc_g(i, j) / fuels%fuelmce(k)
              ! moist damp coef
            etam = 1.0 - 2.59 * rtemp1 + 5.11 * rtemp1 ** 2 - 3.52 * rtemp1 ** 3
              ! mineral damping coef
            etas = 0.174 * fuels%se(k) ** (-0.19)
              !rxn intensity,btu/ft^2 min
            ir = gamma * wn * FUELHEAT * etam * etas
            ! irm = ir * 1055./( 0.3048**2 * 60.) * 1.e-6     !for mw/m^2
            this%iboros(i, j) = ir * 1055.0 / ( 0.3048 ** 2 * 60.0) * 1.e-3 * (60.0 * 384.0 / fuels%savr(k)) ! I_R x t_r (kJ m^-2)
              ! propagating flux ratio
            xifr = exp((0.792 + 0.681 * fuels%savr(k) ** 0.5) &
                * (this%betafl(i, j) + 0.1)) / (192.0 + 0.2595 * fuels%savr(k))

              ! r_0 is the spread rate for a fire on flat ground with no wind.
              ! default spread rate in ft/min
            this%r_0(i, j) = ir * xifr / (rhob * epsilon * qig)
          end if
        end do Loop_i
      end do Loop_j

    end subroutine Set_ros_parameters_wrffire

  end module ros_wrffire_mod
