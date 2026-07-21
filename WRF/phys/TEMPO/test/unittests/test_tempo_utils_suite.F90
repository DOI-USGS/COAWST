module test_tempo_utils_suite
  !! unit tests for module_mp_tempo_utils
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use module_mp_tempo_params, only : wp, sp, dp, t_efrw, t_efsw, &
    initialize_array_efrw, initialize_array_efsw, initialize_bins_for_tables
  use module_mp_tempo_utils, only : calc_gamma_p, compute_efrw, compute_efsw
  implicit none
  private

  public :: collect_tempo_utils_suite

  contains

  ! collect all exported unit tests
  subroutine collect_tempo_utils_suite(testsuite)
    ! collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("checking that procedure calc_gamma_p(8.,8.) = 0.547", &
        CalcGammaP_a8x8), &
      new_unittest("checking initialization of rain-cloud water collection efficiency data after allocation resulting in all(t_efrw == 0.) being false", &
        TableEfrw_AllZeroFalse), &
      new_unittest("checking initialization of snow-cloud water collection efficiency data after allocation resulting in all(t_efsw == 0.) being false", &
        TableEfsw_AllZeroFalse) &
      ]

  end subroutine collect_tempo_utils_suite

  subroutine CalcGammaP_a8x8(error)
    !! test calc_gamma_p with a=8 and x=8
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: value_from_wolfram = 0.54703919051300551 

    call check(error, abs(calc_gamma_p(8.,8.) - value_from_wolfram) < 1.e-6_wp)
    if (allocated(error)) return
  end subroutine CalcGammaP_a8x8


  subroutine TableEfrw_allZeroFalse(error)
    !! test for proper initialization of t_efrw
    type(error_type), allocatable, intent(out) :: error
    call initialize_bins_for_tables() 
    call initialize_array_efrw()
    write(*,*) shape(t_efrw)
    call compute_efrw()
    call check(error, all(t_efrw==0.0_dp), .false.)
    if (allocated(error)) return
  end subroutine TableEfrw_AllZeroFalse


  subroutine TableEfsw_allZeroFalse(error)
    !! test for proper initialization of t_efsw
    type(error_type), allocatable, intent(out) :: error
    call initialize_bins_for_tables() 
    call initialize_array_efsw()
    call compute_efsw()
    call check(error, all(t_efsw==0.0_dp), .false.)
    if (allocated(error)) return
  end subroutine TableEfsw_AllZeroFalse

end module test_tempo_utils_suite