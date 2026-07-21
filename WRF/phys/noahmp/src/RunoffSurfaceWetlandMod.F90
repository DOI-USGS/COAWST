module RunoffSurfaceWetlandMod

!!! Calculate surface runoff based on first-layer soil moisture with wetland scheme (Zhang et al., 2022)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine RunoffSurfaceWetland(noahmp)

! ------------------------ Code history --------------------------------------------------
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Implemented by Z. Zhang (Zhang et al. 2022)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! --------------------------------------------------------------------
    associate(                                                           &
              SoilSfcInflowMean => noahmp%water%flux%SoilSfcInflowMean  ,& ! in,  mean water input on soil surface [m/s]
              SoilSfcSatFracMax => noahmp%water%param%SoilSfcSatFracMax ,& ! in,  maximum surface saturated fraction (global mean)
              SoilImpervFrac    => noahmp%water%state%SoilImpervFrac    ,& ! in,  impervious fraction due to frozen soil
              SoilLiqWater      => noahmp%water%state%SoilLiqWater      ,& ! in,  soil water content [m3/m3] 
              SoilMoistureSat   => noahmp%water%param%SoilMoistureSat   ,& ! in,  saturated value of soil moisture [m3/m3]
              SoilSaturateFrac  => noahmp%water%state%SoilSaturateFrac  ,& ! out, fractional saturated area for soil moisture
              RunoffSurface     => noahmp%water%flux%RunoffSurface      ,& ! out, surface runoff [m/s]
              InfilRateSfc      => noahmp%water%flux%InfilRateSfc        & ! out, infiltration rate at surface [m/s]
             )
! ----------------------------------------------------------------------

    ! compute saturated area fraction, based on top-layer soil moisture
    SoilSaturateFrac = SoilSfcSatFracMax * (SoilLiqWater(1)/SoilMoistureSat(1))

    ! compute surface runoff and infiltration  m/s
    if ( SoilSfcInflowMean > 0.0 ) then
       RunoffSurface = SoilSfcInflowMean * ((1.0-SoilImpervFrac(1)) * SoilSaturateFrac + SoilImpervFrac(1))
       InfilRateSfc  = SoilSfcInflowMean - RunoffSurface 
    endif

    end associate

  end subroutine RunoffSurfaceWetland

end module RunoffSurfaceWetlandMod
