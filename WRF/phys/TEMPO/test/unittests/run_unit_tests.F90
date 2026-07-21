program run_unit_tests
  !! runs tempo unit tests
  use iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  use test_tempo_utils_suite, only : collect_tempo_utils_suite
  use test_tempo_init_suite, only : collect_tempo_init_suite
  use test_tempo_main_suite, only : collect_tempo_main_suite
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
    new_testsuite("tempo_utils_suite", collect_tempo_utils_suite), &
    new_testsuite("tempo_init_suite", collect_tempo_init_suite), &
    new_testsuite("tempo_main_suite", collect_tempo_main_suite) &
    ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program run_unit_tests