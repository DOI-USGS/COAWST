module WetlandWaterZhang22Mod

!!! Calculate wetland water processes (Zhang et al., 2022), evaporation is estimated based on Priestley Taylor (P-T) Method

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine WetlandWaterZhang22(noahmp,TimeStep)

! ------------------------ Code history --------------------------------------------------
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Implementation by Z. Zhang (Zhang et al. 2022)
! ----------------------------------------------------------------------------------------

    implicit none

    type(noahmp_type),      intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: TimeStep       ! timestep (may not be the same as main model timestep)

! local variables
    real(kind=kind_noahmp)                :: LatHeatSpec    ! latent heat vap./sublimation [J/kg]
    real(kind=kind_noahmp)                :: VapPresSlope   ! slope of saturation vapor pressure curve
    real(kind=kind_noahmp)                :: LatHeatPot     ! evaporation from P-T method, energy flux [W/m2]
    real(kind=kind_noahmp)                :: PsychroConst   ! gamma psychrometic constant
    real(kind=kind_noahmp)                :: EvapLatentHeat ! evaporation heat from surface [W/m2]
    real(kind=kind_noahmp)                :: EvapWaterFlux  ! evaporation water from surface [mm/s]

! --------------------------------------------------------------------
    associate(                                                                &
              TemperatureSfc      => noahmp%energy%state%TemperatureSfc      ,& ! in,    surface air temperature [K]
              RadSwAbsSfc         => noahmp%energy%flux%RadSwAbsSfc          ,& ! in,    total absorbed solar radiation [W/m2]
              RadSwReflSfc        => noahmp%energy%flux%RadSwReflSfc         ,& ! in,    total reflected solar radiation [W/m2]
              RadLwNetSfc         => noahmp%energy%flux%RadLwNetSfc          ,& ! in,    total net longwave rad [W/m2] (+ to atm)
              HeatGroundTot       => noahmp%energy%flux%HeatGroundTot        ,& ! in,    total ground heat flux [W/m2] (+ to soil/snow)
              HeatSensibleSfc     => noahmp%energy%flux%HeatSensibleSfc      ,& ! in,    total sensible heat [W/m2] (+ to atm)
              SoilSaturateFrac    => noahmp%water%state%SoilSaturateFrac     ,& ! in,    fractional saturated area for soil moisture
              WetlandCapMax       => noahmp%water%param%WetlandCapMax        ,& ! in,    maximum wetland capacity [m]
              WaterStorageWetland => noahmp%water%state%WaterStorageWetland  ,& ! inout, wetland water storage [mm] 
              RunoffSurface       => noahmp%water%flux%RunoffSurface         ,& ! inout, surface runoff [mm] per soil timestep
              EvapGroundNet       => noahmp%water%flux%EvapGroundNet         ,& ! inout, accumulated net ground evaporation per soil timestep [mm]
              HeatLatentGrd       => noahmp%energy%flux%HeatLatentGrd         & ! inout, ground evaporation heat flux [W/m2] (+ to atm)
             )
! ----------------------------------------------------------------------

    ! set initial value
    EvapWaterFlux  = 0.0
    EvapLatentHeat = 0.0

    ! set psychrometric constant and VapPresSlope
    if ( TemperatureSfc > ConstFreezePoint ) then
       LatHeatSpec = ConstLatHeatEvap
    else
       LatHeatSpec = ConstLatHeatSublim
    endif

    ! determine psychrometic constant
    if ( TemperatureSfc < (273.15+26.85) ) then
       PsychroConst = 0.00040
    else 
       PsychroConst = 0.00041
    endif

    ! calculate slope VapPresSlope based on surface temperature
    if ( TemperatureSfc < (273.15+6.85) ) then
       VapPresSlope = 0.00022
    elseif ( TemperatureSfc < (273.15+16.85) ) then
       VapPresSlope = 0.00042
    elseif ( TemperatureSfc < (273.15+26.85) ) then
       VapPresSlope = 0.00078 
    else 
       VapPresSlope = 0.00132
    endif 

    ! compute wetland water balance
    LatHeatPot = 1.26 * VapPresSlope * (RadSwAbsSfc - RadLwNetSfc - HeatGroundTot) / &
                 (VapPresSlope + PsychroConst)                                           ! LatHeatPot Potential latent heat W/m2 P-T method
    WaterStorageWetland = WaterStorageWetland + RunoffSurface     
    if ( WaterStorageWetland > (LatHeatPot*TimeStep/LatHeatSpec*SoilSaturateFrac) ) then ! if current wetland storage is larger than PET rate 
       EvapWaterFlux = max(LatHeatPot / LatHeatSpec * SoilSaturateFrac, 0.0)             ! evaporation rate mm/s
    else                                                                                 ! if current wetland storage is less than PET rate
       EvapWaterFlux = max(WaterStorageWetland/TimeStep, 0.0)                            ! use all of it
    endif 
    WaterStorageWetland = max(WaterStorageWetland-EvapWaterFlux*TimeStep,0.0)            ! adjust surface wetland storage

    ! adjuest energy and water balance 
    EvapLatentHeat      = EvapWaterFlux * LatHeatSpec                                    ! convert evaporation to latent heat flux 
    HeatSensibleSfc     = HeatSensibleSfc - EvapLatentHeat                               ! reduce sensible heat flux
    HeatLatentGrd       = HeatLatentGrd + EvapLatentHeat                                 ! increase direct evaporation
    RunoffSurface       = max(WaterStorageWetland - WetlandCapMax*1000.0, 0.0)           ! excessive storage becomes runoff
    EvapGroundNet       = EvapGroundNet + EvapWaterFlux                                  ! increase direct evaporation
    WaterStorageWetland = min(WetlandCapMax*1000.0, WaterStorageWetland)

    end associate

  end subroutine WetlandWaterZhang22

end module WetlandWaterZhang22Mod
