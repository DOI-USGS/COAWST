  module fmc_wrffire_mod

    ! Reference: J. Mandel, S. Amram, J.D. Beezley, G. Kelman, A.K. Kochanski, V.Y.
    !            Kondratenko, B.H. Lynn, B. Regev, M. Vejmelka, Recent advances and
    !            applications of WRF-SFIRE. Natural Hazards and Earth System Science, 14,
    !            2829-2845, 2014, doi:10.5194/nhess-14-2829-2014
    !
    ! live fuel moisture guidance following https://www.nwcg.gov/publications/pms437/fuel-moisture/live-fuel-moisture-content
    ! from Rothermel, 1983 Table II-2 p.13 https://www.fs.usda.gov/treesearch/pubs/24635
    ! 300% Fresh foliage, annuals developing early in the growing cycle
    ! 200% Maturing foliage, still developing, with full turgor
    ! 100% Mature foliage, new growth complete and comparable to older perennial foliage
    ! 50%  Entering dormancy, coloration starting, some leaves may have dropped from stem
    ! 30%  Completely cured, treat as dead fuel

    use stderrout_mod, only : Stop_simulation
    use ros_mod, only : ros_t
    use state_mod, only: state_fire_t
    use namelist_mod, only : namelist_t
    use fmc_mod, only : fmc_t
    use fuel_mod, only : fuel_t

    implicit none

    private

    public :: fmc_wrffire_t

    integer, parameter :: MOISTURE_CLASSES = 5, NUM_FMEP = 2
    real, parameter :: FMEP_DECAY_TLAG = 999999 ! time constant of assimilated adjustments of equilibria decay

    type, extends(fmc_t) :: fmc_wrffire_t
      real :: fmoist_lasttime, fmoist_nexttime, dt_moisture
      real, dimension(MOISTURE_CLASSES) :: rec_drying_lag_sec, rec_wetting_lag_sec, fmc_gc_initial_value
      logical :: run_advance_moisture
      real, dimension(:, :, :), allocatable :: fmc_gc ! fuel moisture contents by class
      real, dimension(:, :, :), allocatable :: fmep  ! fuel moisture extended model parameters
      real, dimension (:, :), allocatable :: fmc_gw ! fuel moisture class weights
      real, dimension(:, :, :), allocatable :: fmc_equi ! fuel moisture contents by class equilibrium (diagnostics only)
      real, dimension(:, :, :), allocatable :: fmc_lag ! fuel moisture contents by class time lag (diagnostics only) [h]
    contains
      procedure, public :: Advance_fmc_model => Advance_fmc_model
      procedure, public :: Advance_moisture_classes => Advance_moisture_classes
      procedure, public :: Average_moisture_classes => Average_moisture_classes
      procedure, public :: Init => Init_fmc_wrffire
    end type fmc_wrffire_t

!    character (len = 80), dimension (MOISTURE_CLASSES) :: moisture_class_name
!    data moisture_class_name /'1-h', '10-h', '100-h', '1000-h', 'Live'/
                                                          ! time lag [h] approaching equilibrium moisture
    real, dimension(MOISTURE_CLASSES), parameter :: drying_lag = [ 1.0, 10.0, 100.0, 1000.0, 1e9 ], &
                                                      ! time lag [h] for approaching saturation in rain
                                                    wetting_lag = [ 1.4, 14.0, 140.0, 1400.0, 1e9 ], &
                                                      ! saturation moisture contents (1) in rain
                                                    saturation_moisture = [ 2.5, 2.5, 2.5 ,2.5, 2.5 ], &
                                                      ! stronger rain matters only in duration [mm/h]
                                                    saturation_rain = [ 8.0, 8.0, 8.0, 8.0, 8.0 ], &
                                                      ! rain intensity this small is same as nothing
                                                    rain_threshold = [ 0.05, 0.05, 0.05, 0.05, 0.05 ]
    integer, dimension(MOISTURE_CLASSES), parameter :: drying_model = [ 1, 1, 1, 1, 1 ], &
                                                       wetting_model = [ 1, 1, 1, 1, 1 ], &
                                                         ! 0=input, 1=from fuelmc_g 2=from equilibrium, 3=from fmc_1h,...,fmc_live
                                                       fmc_gc_initialization = [ 2, 2, 2, 2, 3 ]

  contains

    subroutine Init_fmc_wrffire (this, fuels, fuelmc_g, fuelmc_g_live, ifms, ifme, jfms, jfme, itimestep, dt)

      implicit none

      class (fmc_wrffire_t), intent (in out) :: this
      class (fuel_t), intent (in) :: fuels
      integer, intent (in) :: ifms, ifme, jfms, jfme, itimestep
      real, intent (in) :: fuelmc_g, fuelmc_g_live, dt

      real, dimension (:), allocatable :: fgi_t
      real, dimension (:, :), allocatable :: fgi_c
      integer :: i, k, mfuelcats
      character (len = 128) :: msg


      allocate (this%fmc_gc(ifms:ifme, MOISTURE_CLASSES, jfms:jfme))
      allocate (this%fmc_equi(ifms:ifme, MOISTURE_CLASSES, jfms:jfme))
      allocate (this%fmc_lag(ifms:ifme, MOISTURE_CLASSES, jfms:jfme))
      allocate (this%fmep(ifms:ifme, NUM_FMEP, jfms:jfme))
      this%fmep = 0.0

      mfuelcats = fuels%n_fuel_cat + 1
      allocate (this%fmc_gw(mfuelcats, MOISTURE_CLASSES))
      allocate (fgi_t(mfuelcats))
      allocate (fgi_c(mfuelcats, MOISTURE_CLASSES))

      this%fmoist_lasttime = itimestep * dt
      this%fmoist_nexttime = this%fmoist_lasttime

      fgi_c = 0.0
      fgi_c(1:mfuelcats, 1) = fuels%fgi_1h
      fgi_c(1:mfuelcats, 2) = fuels%fgi_10h
      fgi_c(1:mfuelcats, 3) = fuels%fgi_100h
      fgi_c(1:mfuelcats, 4) = fuels%fgi_1000h
      fgi_c(1:mfuelcats, 5) = fuels%fgi_live

      this%fmc_gc_initial_value(1) = fuelmc_g
      this%fmc_gc_initial_value(2) = fuelmc_g
      this%fmc_gc_initial_value(3) = fuelmc_g
      this%fmc_gc_initial_value(4) = fuelmc_g
      this%fmc_gc_initial_value(5) = fuelmc_g_live

        ! Calc averaging weights
      do i = 1, mfuelcats
        fgi_t(i) = 0.0
        do k = 1, MOISTURE_CLASSES
          if (fgi_c(i, k) >= 0.0) then
            fgi_t(i) = fgi_t(i) + fgi_c(i, k)
          else
            write (msg, *) 'fuel load in category', i, ' fuel class ', k, ' is ', fgi_c(i, k), ',must be nonegative.'
            call Stop_simulation (msg)
          end if
        end do

        do k = 1, MOISTURE_CLASSES
          if (fgi_t(i) > 0.0) then
            this%fmc_gw(i, k) = fgi_c(i,k) / fgi_t(i)
          else
            this%fmc_gw(i, k) = 0.0
          end if
        end do
      end do

        ! moisture model derived scalars
      do i = 1, MOISTURE_CLASSES
        this%rec_drying_lag_sec(i) = 1.0 / (3600.0 * drying_lag(i))
        this%rec_wetting_lag_sec(i) = 1.0 / (3600.0 * wetting_lag(i))
      end do

    end subroutine Init_fmc_wrffire

    subroutine Advance_fmc_model (this, fmoist_freq, fmoist_dt, itimestep, dt, ifms, ifme, jfms, jfme, i_start, i_end, j_start, &
            j_end, num_tiles, fire_rain, fire_t2, fire_q2, fire_psfc, fire_rain_old, fire_t2_old, fire_q2_old, &
            fire_psfc_old, fire_rh_fire, fuelmc_g, fmc_g, nfuel_cat, fuels, ros_param)

      implicit none

      class (fmc_wrffire_t), intent (in out) :: this
      class (ros_t), intent (in out) :: ros_param
      class (fuel_t), intent (in) :: fuels
      integer, intent (in) :: fmoist_freq, itimestep, ifms, ifme, jfms, jfme, num_tiles
      integer, dimension(num_tiles) :: i_start, i_end, j_start, j_end
      real, intent (in) ::  fmoist_dt, dt, fuelmc_g
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: nfuel_cat
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: fire_rain, fire_t2, fire_q2, fire_psfc, fire_rain_old, fire_t2_old, &
          fire_q2_old, fire_psfc_old, fire_rh_fire, fmc_g

      integer :: ij
      real :: time_start, moisture_time


      time_start = itimestep * dt
      moisture_time = time_start

      this%run_advance_moisture = .false.
      if (fmoist_freq > 0) then
        if (mod (itimestep, fmoist_freq) == 0) this%run_advance_moisture = .true.
      else
        if (.not. moisture_time < this%fmoist_nexttime) this%run_advance_moisture = .true.
      end if

      if (this%run_advance_moisture) then
        this%dt_moisture = moisture_time - this%fmoist_lasttime
        this%fmoist_lasttime = moisture_time
        if (fmoist_freq > 0)then
            continue
          else
            this%fmoist_nexttime = moisture_time + fmoist_dt
        end if
      end if

      if (this%run_advance_moisture) then
        !$OMP PARALLEL DO   &
        !$OMP PRIVATE (ij)
        do ij = 1, num_tiles
          call this%Advance_moisture_classes (itimestep == 1, ifms, ifme, jfms, jfme, i_start(ij), i_end(ij), j_start(ij), j_end(ij), &
              fire_rain, fire_t2, fire_q2, fire_psfc, fire_rain_old, fire_t2_old, fire_q2_old, fire_psfc_old, fire_rh_fire, fuelmc_g)
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO   &
        !$OMP PRIVATE (ij)
        do ij = 1, num_tiles
          call this%Average_moisture_classes (ifms, ifme, jfms, jfme, i_start(ij), i_end(ij), j_start(ij), j_end(ij), nfuel_cat, fmc_g)
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO   &
        !$OMP PRIVATE (ij)
        do ij = 1, num_tiles
          call ros_param%Set_params (ifms, ifme, jfms, jfme, i_start(ij), i_end(ij), j_start(ij), j_end(ij), &
              fuels, nfuel_cat, fmc_g)
        end do
        !$OMP END PARALLEL DO
      end if

    end subroutine Advance_fmc_model

    subroutine Average_moisture_classes (this, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, nfuel_cat, fmc_g)

      implicit none

      class (fmc_wrffire_t), intent (in) :: this
      integer, intent (in) :: ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte
      real, dimension(ifms:ifme, jfms:jfme), intent (in) :: nfuel_cat
      real, dimension(ifms:ifme, jfms:jfme), intent (out) :: fmc_g

      integer :: i, j, k, n


      do j = jfts, jfte
        do i = ifts, ifte
          fmc_g(i, j) = 0.0
        end do
      end do

      do k = 1, MOISTURE_CLASSES
        do j = jfts, jfte
          do i = ifts, ifte
            n = nfuel_cat(i, j)
            fmc_g(i, j) = fmc_g(i, j) + this%fmc_gw(n, k) * this%fmc_gc(i, k, j)
          end do
        end do
      end do

    end subroutine Average_moisture_classes

    subroutine Advance_moisture_classes (this, initialize, ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte, &
        rain, t2, q2, psfc, rain_old, t2_old, q2_old, psfc_old, rh_fire, fuelmc_g)

      implicit none

      class (fmc_wrffire_t), intent (in out) :: this
      logical, intent (in) :: initialize
      integer, intent (in) :: ifms, ifme, jfms, jfme, ifts, ifte, jfts, jfte
      real, intent (in) :: fuelmc_g
      real, dimension (ifms:ifme, jfms:jfme), intent (in) :: t2, q2, psfc, rain
      real, dimension (ifms:ifme, jfms:jfme), intent (in out) :: t2_old, q2_old, psfc_old, rain_old
      real, intent (out), dimension (ifms:ifme, jfms:jfme) :: rh_fire 

      integer :: i, j, k
      real :: rain_int, T, P, Q, QRS, ES, RH, tend, EMC_d, EMC_w, EMC, R, rain_diff, fmc, rlag, equi, &
          d, w, rhmax, rhmin, change, rainmax,rainmin, fmc_old, H, deltaS, deltaE
      real, parameter :: TOL = 1e-2 ! relative change larger than that will switch to exponential ode solver 
      logical, parameter :: CHECK_RH = .false.
      real :: epsilon, Pws, Pw


      if (initialize) call Copy2old ()

      rhmax = -huge (rhmax)
      rhmin = huge (rhmin)
      rainmax = -huge (rainmax)
      rainmin = huge (rainmin)
      Loop_j: do j = jfts, jfte
        Loop_fmc_classes: do k = 1, MOISTURE_CLASSES
          Loop_i: do i = ifts, ifte
              ! old fuel moisture contents
              ! compute the rain intensity from the difference of accumulated rain
            rain_diff = rain(i, j) - rain_old(i, j)
            if (this%dt_moisture > 0.0)then
              rain_int = 3600.0 * rain_diff / this%dt_moisture
            else
              rain_int = 0.0
            endif
            rainmax = max (rainmax,rain_int)
            rainmin = min (rainmin,rain_int)
            r = rain_int - rain_threshold(k)

              ! average the inputs for second order accuracy
            t = 0.5 * (t2_old(i,j) + t2(i,j))
            p = 0.5 * (psfc_old(i,j) + psfc(i,j))
            q = 0.5 * (q2_old(i,j) + q2(i,j))

              ! compute the relative humidity
              ! ES=610.78*exp(17.269*(T-273.161)/(T-35.861))
              ! QRS=0.622*ES/(P-0.378*ES)
              ! RH = Q/QRS
              ! function rh_from_q from Adam Kochanski following Murphy and Koop, Q.J.R. Meteorol. Soc (2005) 131 1539-1565 eq. (10)
            epsilon = 0.622 ! Molecular weight of water (18.02 g/mol) to molecular weight of dry air (28.97 g/mol)
              ! vapor pressure [Pa]
            pw = q * p / (epsilon + (1.0 - epsilon) * q)
              ! saturation vapor pressure [Pa]
            pws = exp (54.842763 - 6763.22 / t - 4.210 * log (t) + 0.000367 * t + &
                tanh (0.0415 * (t - 218.8)) * (53.878 - 1331.22 / t - 9.44523 * log (t) + 0.014025 * t))
            rh = pw / pws
            rh_fire(i, j) = rh
            rhmax = max (rh, rhmax)         
            rhmin = min (rh, rhmin)         

            deltae = this%fmep(i, 1, j)
            deltas = this%fmep(i, 2, j)

            if (.not. CHECK_RH) then
              rh = min (rh, 1.0)
            else
              if (rh < 0.0 .or. rh > 1.0 .or. rh .ne. rh ) then
                call Stop_simulation ('Relative humidity must be between 0 and 1, saturated water contents must be >0')
              end if
            end if 

            if (r > 0.0) then
              select case (wetting_model(k))
                  ! saturation_moisture=2.5 wetting_lag=14h saturation_rain=8 mm/h calibrated to VanWagner&Pickett 1985 per 24 hours
                case (1)
                  emc_w = saturation_moisture(k) + deltas
                  emc_d = saturation_moisture(k) + deltas
                  rlag = this%rec_wetting_lag_sec(k) * (1.0 - exp(-r / saturation_rain(k)))
              end select
            else
                ! not raining
              select case (drying_model(k))
                  ! Van Wagner and Pickett (1972) per Viney (1991) eq (7) (8)
                case (1)
                  h = rh * 100.0
                    ! equilibrium moisture for drying
                  d = 0.942 * h ** 0.679 + 0.000499 * exp (0.1 * h) + 0.18 * (21.1 + 273.15 - t) * (1.0 - exp(-0.115 * h))
                    ! equilibrium moisture for adsorbtion
                  w = 0.618 * H ** 0.753 + 0.000454 * exp (0.1 * h) + 0.18 * (21.1 + 273.15 - t) * (1.0 - exp(-0.115 * h))
                  if (d .ne. d .or. w .ne. w) call Stop_simulation ('equilibrium moisture calculation failed, result is NaN')
                  d = d * 0.01
                  w = w * 0.01
                  emc_d = max (max (d, w) + deltae, 0.0)
                  emc_w = max (min (d, w) + deltae, 0.0)
                  rlag = this%rec_drying_lag_sec(k)
              end select
            end if
              !*** MODELS THAT ARE NOT OF THE EXPONENTIAL TIME LAG KIND 
              ! ARE RESPONSIBLE FOR THEIR OWN LOGIC, THESE MODELS
              ! SHOULD COMPUTE fmc_gc(i,k,j) DIRECTLY AND SET TLAG < 0
            If_rlag:  if (rlag > 0.0) then
                ! take old from before, no initialization
              if (.not. initialize .or. fmc_gc_initialization(k) == 0) then
                fmc_old = this%fmc_gc(i, k, j)
              else if (fmc_gc_initialization(k) == 1) then
                   ! from scalar fuelmc_g
                fmc_old = fuelmc_g
              else if (fmc_gc_initialization(k) == 2)then
                  ! from computed equilibrium
                fmc_old = 0.5 * (emc_d + emc_w)
              else if (fmc_gc_initialization(k) == 3)then
                  ! from scalar parameter
                fmc_old = this%fmc_gc_initial_value(k)
              else
                call Stop_simulation ('bad value of fmc_gc_initialization(k), must be between 0 and 2')
              end if
                 ! take lower or upper equilibrium value 
              equi = max (min (fmc_old, emc_d), emc_w)

              change = this%dt_moisture * rlag 

              if (change < TOL)then
                   ! 2nd order Taylor (midpoint method)
                 fmc = fmc_old + (equi - fmc_old) * change * (1.0 - 0.5 * change)
              else
                  ! exponential method
                fmc = fmc_old + (equi - fmc_old) * (1 - exp (-change))
              end if
              this%fmc_gc(i, k, j) = fmc

                ! diagnostics out
              this%fmc_equi(i, k, j) = equi
              this%fmc_lag(i, k, j) = 1.0 / (3600.0 * rlag)
               
            end if If_rlag
          end do Loop_i
        end do Loop_fmc_classes
      end do Loop_j

        ! assimilated differences decay
      do j = jfts, jfte
        do k = 1, NUM_FMEP
          do i = ifts, ifte
            change = this%dt_moisture / (FMEP_DECAY_TLAG * 3600.0)
            if (change < TOL) then
              this%fmep(i, k, j) = this%fmep(i, k, j) * (1.0 - change * (1.0 - 0.5 * change))
            else
              this%fmep(i, k, j) = this%fmep(i, k, j) * exp (-change)
            end if
          end do
        end do
      end do

      call Copy2old ()

      return

    contains

      subroutine Copy2old ()

        implicit none

        integer :: i, j


        do j = jfts, jfte
          do i = ifts ,ifte
            rain_old(i, j) = rain(i, j)
            t2_old(i, j) = t2(i, j)
            q2_old(i, j) = q2(i, j)
            psfc_old(i, j) = psfc(i, j)
          end do
        end do

      end subroutine Copy2old

    end subroutine Advance_moisture_classes

  end module fmc_wrffire_mod
