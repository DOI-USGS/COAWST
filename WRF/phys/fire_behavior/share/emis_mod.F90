  module emis_mod

    implicit none

    private

    public :: EMIS_WRFFIRE, EMIS_FMC_PM2P5, Calc_smoke_aod

    integer, parameter :: EMIS_WRFFIRE = 0, EMIS_FMC_PM2P5 = 1

  contains

    pure function Calc_rh (p, t, qv) result (rh)

      implicit none

      real, intent(in) :: p, t, qv
      real :: rh

      real, parameter :: PQ0 = 379.90516, A2 = 17.2693882, A3= 273.16, A4 = 35.86, RHMIN = 1.0, RHMAX = 100.0
      real :: q, qs
      integer :: i, j, k


      q = qv / (1.0 + qv)
      qs = PQ0 / p * exp (A2 * (t - A3) / (t - A4))
      rh = 100.0 * q / qs
      rh = max (min (rh, RHMAX), RHMIN)

    end function Calc_rh

    subroutine Calc_smoke_aod (dz8w, p_phy, t_phy, qv, rho, smoke_tracer, aod5502d_smoke, &
        ids, ide, kds, kde, jds, jde,          &
        ims, ime, kms, kme, jms, jme,          &
        its, ite, kts, kte, jts, jte)

      implicit none

      integer, intent (in) :: ids, ide, kds, kde, jds, jde, &
                              ims, ime, kms, kme, jms, jme, &
                              its, ite, kts, kte, jts, jte
      real, dimension(ims:ime, kms:kme, jms:jme), intent (in) :: dz8w, p_phy, t_phy, qv, rho, smoke_tracer
      real, dimension(ims:ime, jms:jme), intent (out) ::  aod5502d_smoke

                                   ! 4.0 abs + 0.5 scattering [m2 g-1]
      real, parameter :: MASS_EXT_COEF = 4.5, RH_CRIT = 0.3, RH_MAX = 0.95, CONVERT_PERCENT_TO_UNITLESS = 0.01
      real :: rh, augm_ext_coef
      integer :: i, k, j
      logical, parameter :: DEBUG_LOCAL = .false.


      Loop_j_aod : do j = jts, min (jte, jde - 1)
        Loop_i_aod : do i = its, min (ite, ide - 1)
          aod5502d_smoke(i, j) = 0.0
          Loop_k_aod : do k = kts, min (kte, kde - 1)
            rh = CONVERT_PERCENT_TO_UNITLESS * Calc_rh (p_phy(i, k, j), t_phy(i, k, j), qv(i, k, j))
            rh = min (rh, RH_MAX)
            if (rh > RH_CRIT) then
              augm_ext_coef = MASS_EXT_COEF * ((1.0 - RH_CRIT) / (1.0 - rh)) ** 0.18
            else
              augm_ext_coef = MASS_EXT_COEF
            end if                                      !  [m2 g-1]        [g smoke kg-1 air]     [kg air m-3]    [m]
            aod5502d_smoke(i, j) = aod5502d_smoke(i, j) + augm_ext_coef * smoke_tracer(i, k, j) * rho(i, k, j) * dz8w(i, k, j)
            if (DEBUG_LOCAL) call Print_profile (i, j)
          end do Loop_k_aod
        end do Loop_i_aod
      end do Loop_j_aod

    contains

      subroutine Print_profile (i, j)

        use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

        implicit none

        integer, intent (in) :: i, j
        integer, parameter :: I_TO_PRINT = 31, J_TO_PRINT = 34


        if (i == I_TO_PRINT .and. j == J_TO_PRINT) then
          write (OUTPUT_UNIT, *) k, dz8w(i, k, j), rh, augm_ext_coef, rho(i, k, j), smoke_tracer(i, k, j), aod5502d_smoke(i, j)
        end if

      end subroutine Print_profile

    end subroutine Calc_smoke_aod

  end module emis_mod
