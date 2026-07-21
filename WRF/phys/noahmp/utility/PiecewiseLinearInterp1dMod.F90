module PiecewiseLinearInterp1dMod

!!! Piecewise linear interpolation method for 1-dimensional data

  use Machine

  implicit none

contains

  subroutine PiecewiseLinearInterp1d(ND, XD, YD, XI, YI)

! ------------------------ Code history --------------------------------------------------
! Piecewise linear interpolation method for 1-dimensional data
! Original author: John Burkardt, Florida State University, 09/22/2012
! Added and modified by Cenlin He (NCAR) in CTSM, 01/27/2022
! Added in Noah-MP by T.-S. Lin (NCAR), 2024
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    integer               , intent(in)                  :: ND    ! number of data points of (XD)
    real(kind=kind_noahmp), dimension(1:ND), intent(in) :: XD    ! x-value of data points
    real(kind=kind_noahmp), dimension(1:ND), intent(in) :: YD    ! y-value of data points
    real(kind=kind_noahmp), intent(in)                  :: XI    ! x-value for to-be-interpolated point
    real(kind=kind_noahmp), intent(out)                 :: YI    ! the interpolated value at xi

! local variables
    integer                                             :: K     ! loop index
    real(kind=kind_noahmp)                              :: T
! ------------------------------------------------------------

     YI = 0.0

     ! if only one data point
     if ( ND == 1 ) then
        YI = YD(1)
        return
     endif

     ! if multiple data points
     if ( XI < XD(1) ) then ! extrapolate
        T  = ( XI - XD(1) ) / ( XD(2) - XD(1) )
        YI = (1.0 - T) * YD(1) + T * YD(2)
     elseif ( XI > XD(ND) ) then ! extrapolate
        T  = ( XI - XD(ND-1) ) / ( XD(ND) - XD(ND-1) )
        YI = (1.0 - T) * YD(ND-1) + T * YD(ND)
     else  ! piecsewise interpolate
        do K = 2, ND
           if ( (XD(K-1) <= XI) .and. (XI <= XD(K)) ) then
              T  = ( XI - XD(K-1) ) / ( XD(K) - XD(K-1) )
              YI = (1.0 - T) * YD(K-1) + T * YD(K)
              exit
           endif
        enddo
     endif


  end subroutine PiecewiseLinearInterp1d

end module PiecewiseLinearInterp1dMod
