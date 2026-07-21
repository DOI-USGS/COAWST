module SnowFreshRadiusMod

!!! Compute fresh fallen snow grain size for SNICAR albedo calculation
!!! Returns fresh snow grain radius, which is linearly dependent on temperature.
!!! This is implemented to remedy an outstanding bias that SNICAR has in initial
!!! grain size. See e.g. Sandells et al, 2017 for a discussion (10.5194/tc-11-229-2017).
!!! Yang et al. (2017), 10.1016/j.jqsrt.2016.03.033
!!! discusses grain size observations, which suggest a temperature dependence. 

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowFreshRadius(noahmp)

! ------------------------ Code history -----------------------------------
! Implementation: T.-S. Lin, C. He, et al. (2025, JHM)
! Adapted from CTSM function: FreshSnowRadius
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    real(kind=kind_noahmp)           :: Tmin   ! start of linear ramp
    real(kind=kind_noahmp)           :: Tmax   ! end of linear ramp

! --------------------------------------------------------------------
    associate(                                                                   &
              TemperatureAirRefHeight => noahmp%forcing%TemperatureAirRefHeight ,& ! in,  air temperature [K] at reference height
              SnowRadiusMin           => noahmp%water%param%SnowRadiusMin       ,& ! in,  minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
              FreshSnowRadiusMax      => noahmp%water%param%FreshSnowRadiusMax  ,& ! in,  maximum warm fresh snow effective radius [microns]
              SnowRadiusFresh         => noahmp%water%state%SnowRadiusFresh      & ! out, fresh snow radius [microns]
             )
! ----------------------------------------------------------------------

    Tmin = ConstFreezePoint - 30.0
    Tmax = ConstFreezePoint

    if ( FreshSnowRadiusMax <= SnowRadiusMin )then
        SnowRadiusFresh = SnowRadiusMin
    else
        if (TemperatureAirRefHeight < Tmin) then
            SnowRadiusFresh = SnowRadiusMin
        else if (TemperatureAirRefHeight > Tmax) then
            SnowRadiusFresh = FreshSnowRadiusMax
        else
            SnowRadiusFresh = (Tmax - TemperatureAirRefHeight) / (Tmax - Tmin) * SnowRadiusMin + &
                              (TemperatureAirRefHeight - Tmin) / (Tmax - Tmin) * FreshSnowRadiusMax
        end if
    end if

    end associate

  end subroutine SnowFreshRadius

end module SnowFreshRadiusMod
