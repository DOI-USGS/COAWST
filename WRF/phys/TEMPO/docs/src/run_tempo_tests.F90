program run_tempo_tests
  !! runs tempo tests
  use tests, only : test_tempo_init, test_graupel_sedimentation, &
    test_snow_sedimentation, test_cloud_number_aerosolaware, &
    test_cloud_number_non_aerosolaware, test_cloud_number_ml, &
    test_ml_cloud_effective_radius
  
  implicit none

  real, dimension(7) :: sedi_tests = &
    [1., 10., 20., 60., 120., 300., 600.]
  integer :: t

  ! tempo init
  call test_tempo_init()
  
  ! ml cloud effective radius
  call test_ml_cloud_effective_radius(dt=20.)
  
  ! test cloud number concentration
  call test_cloud_number_aerosolaware(dt=20.)
  call test_cloud_number_non_aerosolaware(dt=20.)
  call test_cloud_number_ml(dt=20.)

  ! graupel sedimentation
  do t = 1, size(sedi_tests)
    call test_graupel_sedimentation(dt=sedi_tests(t), semi_sedi=.false.)
  enddo
  do t = 1, size(sedi_tests)
    call test_graupel_sedimentation(dt=sedi_tests(t), semi_sedi=.true.)
  enddo

  ! snow sedimentation
  do t = 1, size(sedi_tests)
    call test_snow_sedimentation(dt=sedi_tests(t))
  enddo

end program run_tempo_tests
