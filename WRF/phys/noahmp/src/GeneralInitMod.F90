module GeneralInitMod

!!! General initialization for variables

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
 
  implicit none

contains

  subroutine GeneralInit(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in NOAHMP_SFLX)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd   ! loop index

! --------------------------------------------------------------------
    associate(                                                                       &
              LandUseDataName        => noahmp%config%domain%LandUseDataName        ,& ! in,  landuse data name (USGS or MODIS_IGBP)
              VegType                => noahmp%config%domain%VegType                ,& ! in,  vegetation type
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,  number of soil layers
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,  depth [m] of layer-bottom from soil surface
              NumSoilLayerRoot       => noahmp%water%param%NumSoilLayerRoot         ,& ! in,  number of soil layers with root present
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg        ,& ! in,  actual number of snow layers (negative)
              DepthSnowSoilLayer     => noahmp%config%domain%DepthSnowSoilLayer     ,& ! in,  depth of snow/soil layer-bottom [m]
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow     ,& ! in,  snow and soil layer temperature [K]
              FlagCropland           => noahmp%config%domain%FlagCropland           ,& ! out, flag to identify croplands
              FlagWetland            => noahmp%config%domain%FlagWetland            ,& ! out, flag to identify wetlands
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! out, thickness of snow/soil layers [m]
              TemperatureRootZone    => noahmp%energy%state%TemperatureRootZone      & ! out, root-zone averaged temperature [K]
             )
! ----------------------------------------------------------------------

    ! initialize snow/soil layer thickness
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       if ( LoopInd == NumSnowLayerNeg+1 ) then
          ThicknessSnowSoilLayer(LoopInd) = - DepthSnowSoilLayer(LoopInd)
       else
          ThicknessSnowSoilLayer(LoopInd) = DepthSnowSoilLayer(LoopInd-1) - DepthSnowSoilLayer(LoopInd)
       endif
    enddo

    ! initialize root-zone soil temperature
    TemperatureRootZone = 0.0
    do LoopInd = 1, NumSoilLayerRoot
       TemperatureRootZone = TemperatureRootZone + &
                             TemperatureSoilSnow(LoopInd) * ThicknessSnowSoilLayer(LoopInd) / (-DepthSoilLayer(NumSoilLayerRoot))
    enddo

    ! initialize special land type flags
    FlagCropland = .false.
    FlagWetland  = .false.
    if ( trim(LandUseDataName) == "USGS" ) then
       if ( (VegType >= 3 ) .and. (VegType <= 6 ) ) FlagCropland = .true.
       if ( (VegType >= 17) .and. (VegType <= 18) ) FlagWetland  = .true.
    elseif ( trim(LandUseDataName) == "MODIFIED_IGBP_MODIS_NOAH") then
       if ( (VegType == 12) .or. (VegType == 14) )  FlagCropland = .true.
       if ( (VegType == 11) )                       FlagWetland  = .true.
    endif

    end associate

  end subroutine GeneralInit

end module GeneralInitMod
