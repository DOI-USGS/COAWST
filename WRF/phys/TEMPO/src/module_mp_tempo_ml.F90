module module_mp_tempo_ml
  !! neural network module designed for cloud droplet number concentration prediction 
  !!
  !! inspired by work from [David John Gagne (NCAR)](https://github.com/NCAR/mlmicrophysics)

  use module_mp_tempo_params, only : wp, sp, dp
  
  implicit none
  private
  
  public :: ty_tempo_ml_data, save_or_read_ml_data, tempo_ml_predict_cloud_number, &
    nc_ml_nodes, nc_ml_input, nc_ml_output, nc_ml_trans_mean, nc_ml_trans_var, &
    nc_ml_w00, nc_ml_w01, nc_ml_b00, nc_ml_b01
   
  type ty_tempo_ml_data
    integer :: input_size
    integer :: output_size
    integer :: node_size
    real(wp), dimension(:), allocatable :: transform_mean
    real(wp), dimension(:), allocatable :: transform_var
    real(wp), dimension(:,:), allocatable :: weights00
    real(wp), dimension(:), allocatable :: bias00
    real(wp), dimension(:,:), allocatable :: weights01
    real(wp), dimension(:), allocatable :: bias01
  end type ty_tempo_ml_data

  ! input / output dimensions
  integer, parameter :: nc_ml_input = 7
  integer, parameter :: nc_ml_nodes = 24
  integer, parameter :: nc_ml_output = 1

  real(wp), dimension(nc_ml_input), parameter :: &
    nc_ml_trans_mean = [0.000191556196486247_wp, 3.58145042772654e-05_wp, &
    3.12611085359273e-07_wp, 5.74303078738579e-05_wp, 84191.2092225319_wp, &
    279.070551773565_wp, 0.123679354084004_wp]
  real(wp), dimension(nc_ml_input), parameter :: &
    nc_ml_trans_var = [5.78143564171777e-08_wp, 3.22834309750552e-08_wp, &
    6.45745893307455e-11_wp, 4.16625579383794e-08_wp, 215694631.771185_wp, &
    94.6576255386858_wp, 0.384841247662964_wp]

  real(wp), dimension(nc_ml_input * nc_ml_nodes), parameter :: &
    nc_ml_w00 = [-2.006957_wp, -0.2812008_wp, -0.339073_wp, 1.596426_wp, 2.395225_wp, 1.76315_wp, &
    -0.0626798_wp, 1.267002_wp, -0.02234177_wp, -6.522605e-33_wp, 0.4792154_wp, -0.1253034_wp, &
    -3.217191_wp, -3.092887_wp, -0.0863651_wp, 1.071625_wp, 0.09741028_wp, 0.2255831_wp, &
    -0.6929023_wp, -0.02693799_wp, -3.432344e-33_wp, -0.8791879_wp, -0.9359049_wp, 1.083484_wp, &
    -0.07909214_wp, -0.0122418_wp, -0.02815927_wp, 0.1676407_wp, 0.08252326_wp, 0.6697816_wp, &
    -0.4019359_wp, 0.4687141_wp, 0.001813132_wp, -9.792186e-33_wp, 0.0409322_wp, 0.0113192_wp, &
    -0.01354596_wp, 0.00307771_wp, -0.4635534_wp, 0.03835761_wp, -0.1015553_wp, 0.7316446_wp, &
    -0.05791711_wp, -0.0002690362_wp, -7.920147e-33_wp, -0.1216918_wp, -0.3190572_wp, 0.09809405_wp, &
    -0.16476_wp, -0.03387314_wp, 0.005422261_wp, 0.04043967_wp, 0.03901243_wp, 0.07444729_wp, &
    0.01954299_wp, 0.06918761_wp, 0.04823543_wp, -8.637957e-33_wp, 0.06371575_wp, -0.09250915_wp, &
    -1.109653_wp, -1.373999_wp, -0.2412623_wp, -0.04482195_wp, 0.1584691_wp, 0.06353725_wp, &
    0.0006248798_wp, 0.04593191_wp, -8.878673e-33_wp, -0.4988684_wp, 0.01110262_wp, 0.04623203_wp, &
    0.006581791_wp, 0.03536217_wp, -0.1890567_wp, -0.08839592_wp, 0.1327181_wp, 0.03478973_wp, &
    -0.1565902_wp, -0.100401_wp, -0.1179777_wp, -8.879818e-33_wp, -0.1383738_wp, 0.02847495_wp, &
    -0.005902881_wp, 0.005615512_wp, -0.6308192_wp, -0.02431803_wp, -0.141971_wp, -0.3490018_wp, &
    -0.9850957_wp, -0.1449479_wp, -8.059166e-33_wp, -0.1186465_wp, -1.165381_wp, 0.069015_wp, &
    0.003388841_wp, -0.04041302_wp, 0.1638467_wp, 0.1147008_wp, -0.04833491_wp, -0.07755993_wp, &
    -0.5137688_wp, 0.04546477_wp, 0.04101883_wp, 6.752353e-33_wp, 0.3541977_wp, 0.04880851_wp, &
    -0.00102834_wp, -0.01280629_wp, -0.1116254_wp, -0.02204754_wp, 0.07100908_wp, 0.2354002_wp, &
    0.07129629_wp, 0.2489657_wp, 8.080785e-33_wp, -0.03449865_wp, 0.06037927_wp, -0.02023619_wp, &
    0.7779589_wp, 0.04680278_wp, 0.7492616_wp, 0.6545208_wp, -1.09497_wp, -1.176524_wp, &
    -0.451585_wp, 0.881124_wp, -0.4551499_wp, -8.708624e-33_wp, 1.006558_wp, -0.04979523_wp, &
    -0.0006915367_wp, -0.002993054_wp, 0.01654614_wp, -0.07141764_wp, -0.2216591_wp, 0.8637336_wp, &
    0.8358089_wp, -0.7576646_wp, 8.186339e-33_wp, 0.1161914_wp, 0.7121871_wp, -1.146734_wp, &
    0.03925339_wp, 0.9975697_wp, -0.06953461_wp, -0.07598846_wp, -0.06418022_wp, -0.01897495_wp, &
    -0.03612464_wp, -0.06703389_wp, -0.103049_wp, 9.364858e-33_wp, -0.03949085_wp, 1.005938_wp, &
    0.0001625419_wp, -0.01027657_wp, 0.03823901_wp, 0.3197246_wp, -0.1160718_wp, 0.04449431_wp, &
    0.009831783_wp, 0.6551342_wp, -8.470687e-33_wp, 0.1374123_wp, -0.01782138_wp, 0.01719949_wp]
  real(wp), dimension(nc_ml_nodes), parameter :: &
    nc_ml_w01 = [-4.045869_wp, 0.558111_wp, -1.567351_wp, -1.64972_wp, 2.748608_wp, &
    -1.909901_wp, 0.3955558_wp, 1.507247_wp, 0.3599722_wp, -0.0001173223_wp, -0.9398569_wp, &
    -0.6867028_wp, -61.72853_wp, -63.13766_wp, 1.165811_wp, -0.6848684_wp, 0.1931683_wp, &
    1.1208_wp, 2.63087_wp, 0.740169_wp, -21.62499_wp, 1.545568_wp, 3.575141_wp, -1.299604_wp]
  real(wp), dimension(nc_ml_nodes), parameter :: &
    nc_ml_b00 = [-0.9842531_wp, 0.3064759_wp, -0.4500185_wp, 1.28336_wp, 1.384105_wp, 0.528031_wp, &
    0.8453538_wp, 1.579872_wp, 2.245679_wp, -0.008679952_wp, 0.4549862_wp, -0.136581_wp, &
    -2.576741_wp, -2.483647_wp, -0.2089484_wp, 0.7607977_wp, 1.847745_wp, 0.7316047_wp, &
    -0.287945_wp, 2.227298_wp, -2.314714_wp, -0.2561245_wp, -0.6993448_wp, -0.1359731_wp]
  real(wp), dimension(nc_ml_output), parameter :: &
    nc_ml_b01 = [1.572826_wp]

  contains

  subroutine save_or_read_ml_data(ml_data_in, ml_data_out)
    !! initializes and saves or returns neural network information

    logical, save :: not_initialized = .true.
    type(ty_tempo_ml_data), dimension(1), intent(in), optional :: ml_data_in
    type(ty_tempo_ml_data), dimension(1), intent(out), optional :: ml_data_out
    type(ty_tempo_ml_data), dimension(1), save :: tempo_ml_data_save

    if (not_initialized) then
      if (present(ml_data_in)) tempo_ml_data_save = ml_data_in
      not_initialized = .false.
    endif

    if (present(ml_data_out)) then
      ml_data_out = tempo_ml_data_save
    endif
  end subroutine save_or_read_ml_data


  subroutine tempo_ml_predict_cloud_number(qc, qr, qi, qs, pres, temp, w, &
    predicted_number)
    !! predicts number concentration

    real(wp), dimension(:), intent(in) :: qc, qr, qi, qs, pres, temp, w
    real(wp), dimension(:), intent(inout) :: predicted_number

    type(ty_tempo_ml_data), dimension(1) :: get_ml_data
    type(ty_tempo_ml_data) :: ml_data
    integer, parameter :: input_rows = 1

    real(wp) :: input(nc_ml_input, size(qc))
    real(wp) :: input_transformed(nc_ml_input, size(qc))
    real(wp) :: output00(nc_ml_nodes, size(qc))
    real(wp) :: output00_activ(nc_ml_nodes, size(qc))
    real(wp) :: reshaped_bias00(nc_ml_nodes, size(qc))
    real(wp) :: output01(nc_ml_output, size(qc))
    real(wp) :: output01_activ(nc_ml_output, size(qc))
    real(wp) :: reshaped_bias01(nc_ml_output, size(qc))

    real(wp), parameter :: logMin = -6.0_wp       ! r2
    real(wp), parameter :: logMax = 9.3010299957_wp ! 2000 cm^-3
    real(wp) :: predicted_exp, bias_corr
    integer :: k, nz

    ! get neural network data
    call save_or_read_ml_data(ml_data_out=get_ml_data)
    ml_data = get_ml_data(1)
    nz = size(qc)

    ! collect input data
    input(1,:) = qc
    input(2,:) = qr
    input(3,:) = qi
    input(4,:) = qs
    input(5,:) = pres
    input(6,:) = temp
    input(7,:) = w

    ! transform input data
    call standard_scaler_transform(mean=ml_data%transform_mean, var=ml_data%transform_var, &
         raw_data=input, transformed_data=input_transformed)

    do k = 1, nz
      reshaped_bias00(:,k) = ml_data%bias00
      reshaped_bias01(1,k) = ml_data%bias01(1)
    enddo

    ! reconstruct neural network
    ! first layer
    output00 = matmul(ml_data%weights00, input_transformed) + reshaped_bias00
    call relu_activation(input=output00, output=output00_activ)

    ! second layer
    output01 = matmul(ml_data%weights01, output00_activ) + reshaped_bias01
    call relu_activation(input=output01, output=output01_activ)

    ! prediction
    do k = 1, nz
      predicted_exp = min(logMax, max(logMin, output01_activ(1,k)))

      ! bias correction
      bias_corr = 1.0_wp
      if ((predicted_exp >= 0._wp) .and. (predicted_exp < 3._wp)) then
        bias_corr = -0.2704_wp*predicted_exp**5 + 1.838_wp*predicted_exp**4 - &
          5.127_wp*predicted_exp**3 + 8.547_wp*predicted_exp**2 - &
          8.439_wp*predicted_exp + 4.297_wp
      endif
      predicted_number(k) = bias_corr * (10._wp**predicted_exp)
    enddo
  end subroutine tempo_ml_predict_cloud_number


  subroutine standard_scaler_transform(mean, var, raw_data, transformed_data)
    !! standard scaler transformer
  
    real(wp), dimension(:,:), intent(in) :: raw_data
    real(wp), dimension(:), intent(in) :: mean, var
    real(wp), dimension(:,:), intent(out) :: transformed_data
    integer :: i
    
    do i = 1, size(raw_data, 1)
      transformed_data(i,:) = (raw_data(i,:) - mean(i)) / sqrt(var(i))
    end do
  end subroutine standard_scaler_transform


  subroutine relu_activation(input, output)
   !! relu activation function

    real(wp), dimension(:,:), intent(in) :: input
    real(wp), dimension(:,:), intent(out) :: output
    integer :: i, j
    
    do i = 1, size(input, 1)
      do j = 1, size(input, 2)
        output(i, j) = max(input(i,j), 0._wp)
      end do
    end do
  end subroutine relu_activation
  
end module module_mp_tempo_ml
