module test_tempo_init_suite
  !! unit tests for module_mp_tempo_init
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use module_mp_tempo_params, only : wp, sp, dp, &
    initialize_graupel_vars, dim_nrhg, nrhg, nrhg1, rho_g, am_g, pi, &
    initialize_parameters, ccg, &
    initialize_bins_for_tables, dc, nbc
  implicit none
  private

  public :: collect_tempo_init_suite

  contains

  ! collect all exported unit tests
  subroutine collect_tempo_init_suite(testsuite)
    ! collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("checking thatdim_nrhg=9 when hail-aware = true", &
        InitializeGraupelVars_dimNrhg9_when_HailAwareTrue), &
      new_unittest("checking that dim_nrhg=1 with hail-aware = false", &
        InitializeGraupelVars_dimNrhg1_when_HailAwareFalse), &
      new_unittest("checking that ccg(2,4) = 7! = 5040", &
        InitializeParameters_ccgValue5040), &
      new_unittest("checking the value of the last dc bin, dc(nbc) = 100 microns", &
        InitializeBinsForTables_lastDcBin100microns) &
      ]
  end subroutine collect_tempo_init_suite

    
  subroutine InitializeGraupelVars_dimNrhg9_when_HailAwareTrue(error)
    !! test hail-aware initialization
    type(error_type), allocatable, intent(out) :: error

    call initialize_graupel_vars(hail_flag=.true.)
    call check(error, dim_nrhg, nrhg)
    if (allocated(error)) return
  end subroutine InitializeGraupelVars_dimNrhg9_when_HailAwareTrue


  subroutine InitializeGraupelVars_dimNrhg1_when_HailAwareFalse(error)
    !! test non-hail-aware initialization
    type(error_type), allocatable, intent(out) :: error

    call initialize_graupel_vars(hail_flag=.false.)
    call check(error, dim_nrhg, nrhg1)
    if (allocated(error)) return
  end subroutine InitializeGraupelVars_dimNrhg1_when_HailAwareFalse


  subroutine InitializeParameters_ccgValue5040(error)
    !! test value of ccg(2,4)
    type(error_type), allocatable, intent(out) :: error

    call initialize_parameters()
    call check(error, ccg(2,4), 5040._wp)
    if (allocated(error)) return
  end subroutine InitializeParameters_ccgValue5040


  subroutine InitializeBinsForTables_lastDcBin100microns(error)
    ! test value of last dc bin, dc(nbc)
    type(error_type), allocatable, intent(out) :: error

    call initialize_bins_for_tables()
    call check(error, dc(nbc), 100.e-6_wp)
    if (allocated(error)) return
  end subroutine InitializeBinsForTables_lastDcBin100microns

end module test_tempo_init_suite