module WaterVarInitMod

!!! Initialize column (1-D) Noah-MP water variables
!!! Water variables should be first defined in WaterVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpVarType

  implicit none

contains

!=== initialize with default values
  subroutine WaterVarInitDefault(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    associate(                                                                         &
              NumSnowLayerMax         => noahmp%config%domain%NumSnowLayerMax         ,&
              NumSoilLayer            => noahmp%config%domain%NumSoilLayer            ,&
              NumDensitySnwAgeSnicar  => noahmp%config%domain%NumDensitySnwAgeSnicar  ,&
              NumTempGradSnwAgeSnicar => noahmp%config%domain%NumTempGradSnwAgeSnicar ,&
              NumTempSnwAgeSnicar     => noahmp%config%domain%NumTempSnwAgeSnicar      &
             )

    ! water state variables
    noahmp%water%state%IrrigationCntSprinkler      = undefined_int
    noahmp%water%state%IrrigationCntMicro          = undefined_int
    noahmp%water%state%IrrigationCntFlood          = undefined_int
    noahmp%water%state%IrrigationFracFlood         = undefined_real
    noahmp%water%state%IrrigationAmtFlood          = undefined_real
    noahmp%water%state%IrrigationFracMicro         = undefined_real
    noahmp%water%state%IrrigationAmtMicro          = undefined_real
    noahmp%water%state%IrrigationFracSprinkler     = undefined_real
    noahmp%water%state%IrrigationAmtSprinkler      = undefined_real
    noahmp%water%state%IrrigationFracGrid          = undefined_real
    noahmp%water%state%CanopyLiqWater              = undefined_real
    noahmp%water%state%CanopyIce                   = undefined_real
    noahmp%water%state%CanopyTotalWater            = undefined_real
    noahmp%water%state%CanopyWetFrac               = undefined_real
    noahmp%water%state%CanopyIceMax                = undefined_real
    noahmp%water%state%CanopyLiqWaterMax           = undefined_real
    noahmp%water%state%SnowfallDensity             = undefined_real
    noahmp%water%state%SnowDepth                   = undefined_real
    noahmp%water%state%SnowWaterEquiv              = undefined_real
    noahmp%water%state%SnowWaterEquivPrev          = undefined_real
    noahmp%water%state%SnowCoverFrac               = undefined_real
    noahmp%water%state%PondSfcThinSnwMelt          = undefined_real
    noahmp%water%state%PondSfcThinSnwComb          = undefined_real
    noahmp%water%state%PondSfcThinSnwTrans         = undefined_real
    noahmp%water%state%SoilIceMax                  = undefined_real
    noahmp%water%state%SoilLiqWaterMin             = undefined_real
    noahmp%water%state%SoilSaturateFrac            = undefined_real
    noahmp%water%state%SoilImpervFracMax           = undefined_real
    noahmp%water%state%SoilMoistureToWT            = undefined_real
    noahmp%water%state%SoilTranspFacAcc            = undefined_real
    noahmp%water%state%SoilWaterRootZone           = undefined_real
    noahmp%water%state%SoilWaterStress             = undefined_real
    noahmp%water%state%SoilSaturationExcess        = undefined_real
    noahmp%water%state%RechargeGwDeepWT            = undefined_real
    noahmp%water%state%RechargeGwShallowWT         = undefined_real
    noahmp%water%state%WaterTableHydro             = undefined_real
    noahmp%water%state%WaterTableDepth             = undefined_real
    noahmp%water%state%WaterStorageAquifer         = undefined_real
    noahmp%water%state%WaterStorageSoilAqf         = undefined_real
    noahmp%water%state%WaterStorageLake            = undefined_real
    noahmp%water%state%WaterStorageTotBeg          = undefined_real
    noahmp%water%state%WaterBalanceError           = undefined_real
    noahmp%water%state%WaterStorageTotEnd          = undefined_real
    noahmp%water%state%WaterHeadSfc                = undefined_real
    noahmp%water%state%PrecipAreaFrac              = undefined_real
    noahmp%water%state%TileDrainFrac               = undefined_real
    noahmp%water%state%FrozenPrecipFrac            = undefined_real
    noahmp%water%state%WaterStorageWetland         = undefined_real

    if ( .not. allocated(noahmp%water%state%IndexPhaseChange) )     &
       allocate( noahmp%water%state%IndexPhaseChange(-NumSnowLayerMax+1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilSupercoolWater) )   &
       allocate( noahmp%water%state%SoilSupercoolWater(-NumSnowLayerMax+1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SnowIce) )              &
       allocate( noahmp%water%state%SnowIce(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%state%SnowLiqWater) )         &
       allocate( noahmp%water%state%SnowLiqWater(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%state%SnowIceVol) )           &
       allocate( noahmp%water%state%SnowIceVol(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%state%SnowLiqWaterVol) )      &
       allocate( noahmp%water%state%SnowLiqWaterVol(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%state%SnowIceFracPrev) )      &
       allocate( noahmp%water%state%SnowIceFracPrev(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%state%SnowIceFrac) )          &
       allocate( noahmp%water%state%SnowIceFrac(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%state%SnowEffPorosity) )      &
       allocate( noahmp%water%state%SnowEffPorosity(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%state%SoilLiqWater) )         &
       allocate( noahmp%water%state%SoilLiqWater(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilIce) )              &
       allocate( noahmp%water%state%SoilIce(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilMoisture) )         &
       allocate( noahmp%water%state%SoilMoisture(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilImpervFrac) )       &
       allocate( noahmp%water%state%SoilImpervFrac(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilWatConductivity) )  &
       allocate( noahmp%water%state%SoilWatConductivity(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilWatDiffusivity) )   &
       allocate( noahmp%water%state%SoilWatDiffusivity(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilEffPorosity) )      &
       allocate( noahmp%water%state%SoilEffPorosity(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilIceFrac) )          &
       allocate( noahmp%water%state%SoilIceFrac(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilMoistureEqui) )     &
       allocate( noahmp%water%state%SoilMoistureEqui(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilTranspFac) )        &
       allocate( noahmp%water%state%SoilTranspFac(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%state%SoilMatPotential) )     &
       allocate( noahmp%water%state%SoilMatPotential(1:NumSoilLayer) )

    noahmp%water%state%IndexPhaseChange   (:)      = undefined_int
    noahmp%water%state%SoilSupercoolWater (:)      = undefined_real
    noahmp%water%state%SnowIce            (:)      = undefined_real
    noahmp%water%state%SnowLiqWater       (:)      = undefined_real
    noahmp%water%state%SnowIceVol         (:)      = undefined_real
    noahmp%water%state%SnowLiqWaterVol    (:)      = undefined_real
    noahmp%water%state%SnowIceFracPrev    (:)      = undefined_real
    noahmp%water%state%SnowIceFrac        (:)      = undefined_real
    noahmp%water%state%SoilIceFrac        (:)      = undefined_real
    noahmp%water%state%SnowEffPorosity    (:)      = undefined_real
    noahmp%water%state%SoilLiqWater       (:)      = undefined_real
    noahmp%water%state%SoilIce            (:)      = undefined_real
    noahmp%water%state%SoilMoisture       (:)      = undefined_real
    noahmp%water%state%SoilImpervFrac     (:)      = undefined_real
    noahmp%water%state%SoilWatConductivity(:)      = undefined_real
    noahmp%water%state%SoilWatDiffusivity (:)      = undefined_real
    noahmp%water%state%SoilEffPorosity    (:)      = undefined_real
    noahmp%water%state%SoilMoistureEqui   (:)      = undefined_real
    noahmp%water%state%SoilTranspFac      (:)      = undefined_real
    noahmp%water%state%SoilMatPotential   (:)      = undefined_real

    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 )then
       if ( .not. allocated(noahmp%water%state%SnowRadius) )             &
       allocate( noahmp%water%state%SnowRadius(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassBChydropho) )         &
       allocate( noahmp%water%state%MassBChydropho(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassBChydrophi) )         &
       allocate( noahmp%water%state%MassBChydrophi(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassOChydropho) )         &
       allocate( noahmp%water%state%MassOChydropho(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassOChydrophi) )         &
       allocate( noahmp%water%state%MassOChydrophi(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassDust1) )              &
       allocate( noahmp%water%state%MassDust1(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassDust2) )              &
       allocate( noahmp%water%state%MassDust2(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassDust3) )              &
       allocate( noahmp%water%state%MassDust3(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassDust4) )              &
       allocate( noahmp%water%state%MassDust4(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassDust5) )              &
       allocate( noahmp%water%state%MassDust5(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassConcBChydropho) )     &
       allocate( noahmp%water%state%MassConcBChydropho(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassConcBChydrophi) )     &
       allocate( noahmp%water%state%MassConcBChydrophi(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassConcOChydropho) )     &
       allocate( noahmp%water%state%MassConcOChydropho(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassConcOChydrophi) )     &
       allocate( noahmp%water%state%MassConcOChydrophi(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassConcDust1) )          &
       allocate( noahmp%water%state%MassConcDust1(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassConcDust2) )          &
       allocate( noahmp%water%state%MassConcDust2(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassConcDust3) )          &
       allocate( noahmp%water%state%MassConcDust3(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassConcDust4) )          &
       allocate( noahmp%water%state%MassConcDust4(-NumSnowLayerMax+1:0) )
       if ( .not. allocated(noahmp%water%state%MassConcDust5) )          &
       allocate( noahmp%water%state%MassConcDust5(-NumSnowLayerMax+1:0) )

       noahmp%water%state%SnowRadiusFresh          = undefined_real
       noahmp%water%state%SnowRadius        (:)    = undefined_real
       noahmp%water%state%MassBChydropho    (:)    = undefined_real
       noahmp%water%state%MassBChydrophi    (:)    = undefined_real
       noahmp%water%state%MassOChydropho    (:)    = undefined_real
       noahmp%water%state%MassOChydrophi    (:)    = undefined_real
       noahmp%water%state%MassDust1         (:)    = undefined_real
       noahmp%water%state%MassDust2         (:)    = undefined_real
       noahmp%water%state%MassDust3         (:)    = undefined_real
       noahmp%water%state%MassDust4         (:)    = undefined_real
       noahmp%water%state%MassDust5         (:)    = undefined_real
       noahmp%water%state%MassConcBChydropho(:)    = undefined_real
       noahmp%water%state%MassConcBChydrophi(:)    = undefined_real
       noahmp%water%state%MassConcOChydropho(:)    = undefined_real
       noahmp%water%state%MassConcOChydrophi(:)    = undefined_real
       noahmp%water%state%MassConcDust1     (:)    = undefined_real
       noahmp%water%state%MassConcDust2     (:)    = undefined_real
       noahmp%water%state%MassConcDust3     (:)    = undefined_real
       noahmp%water%state%MassConcDust4     (:)    = undefined_real
       noahmp%water%state%MassConcDust5     (:)    = undefined_real
    endif

    ! water flux variables
    noahmp%water%flux%PrecipTotRefHeight           = undefined_real
    noahmp%water%flux%RainfallRefHeight            = undefined_real
    noahmp%water%flux%SnowfallRefHeight            = undefined_real
    noahmp%water%flux%PrecipConvTotRefHeight       = undefined_real
    noahmp%water%flux%PrecipLargeSclRefHeight      = undefined_real
    noahmp%water%flux%EvapCanopyNet                = undefined_real
    noahmp%water%flux%Transpiration                = undefined_real
    noahmp%water%flux%EvapCanopyLiq                = undefined_real
    noahmp%water%flux%DewCanopyLiq                 = undefined_real
    noahmp%water%flux%FrostCanopyIce               = undefined_real
    noahmp%water%flux%SublimCanopyIce              = undefined_real
    noahmp%water%flux%MeltCanopyIce                = undefined_real
    noahmp%water%flux%FreezeCanopyLiq              = undefined_real
    noahmp%water%flux%SnowfallGround               = undefined_real
    noahmp%water%flux%SnowDepthIncr                = undefined_real
    noahmp%water%flux%FrostSnowSfcIce              = undefined_real
    noahmp%water%flux%SublimSnowSfcIce             = undefined_real
    noahmp%water%flux%RainfallGround               = undefined_real
    noahmp%water%flux%SnowBotOutflow               = undefined_real
    noahmp%water%flux%GlacierExcessFlow            = undefined_real
    noahmp%water%flux%SoilSfcInflow                = undefined_real
    noahmp%water%flux%RunoffSurface                = undefined_real
    noahmp%water%flux%RunoffSubsurface             = undefined_real
    noahmp%water%flux%InfilRateSfc                 = undefined_real
    noahmp%water%flux%EvapSoilSfcLiq               = undefined_real
    noahmp%water%flux%DrainSoilBot                 = undefined_real
    noahmp%water%flux%RechargeGw                   = undefined_real
    noahmp%water%flux%DischargeGw                  = undefined_real
    noahmp%water%flux%VaporizeGrd                  = undefined_real
    noahmp%water%flux%CondenseVapGrd               = undefined_real
    noahmp%water%flux%DewSoilSfcLiq                = undefined_real
    noahmp%water%flux%InterceptCanopyRain          = undefined_real
    noahmp%water%flux%DripCanopyRain               = undefined_real
    noahmp%water%flux%ThroughfallRain              = undefined_real
    noahmp%water%flux%InterceptCanopySnow          = undefined_real
    noahmp%water%flux%DripCanopySnow               = undefined_real
    noahmp%water%flux%ThroughfallSnow              = undefined_real
    noahmp%water%flux%EvapGroundNet                = undefined_real
    noahmp%water%flux%MeltGroundSnow               = undefined_real
    noahmp%water%flux%WaterToAtmosTotal            = undefined_real
    noahmp%water%flux%EvapSoilSfcLiqAcc            = undefined_real
    noahmp%water%flux%SoilSfcInflowAcc             = undefined_real
    noahmp%water%flux%GlacierExcessFlowAcc         = undefined_real
    noahmp%water%flux%SfcWaterTotChgAcc            = undefined_real
    noahmp%water%flux%PrecipTotAcc                 = undefined_real
    noahmp%water%flux%EvapCanopyNetAcc             = undefined_real
    noahmp%water%flux%TranspirationAcc             = undefined_real
    noahmp%water%flux%EvapGroundNetAcc             = undefined_real
    noahmp%water%flux%EvapSoilSfcLiqMean           = undefined_real
    noahmp%water%flux%SoilSfcInflowMean            = undefined_real
    noahmp%water%flux%IrrigationRateFlood          = 0.0
    noahmp%water%flux%IrrigationRateMicro          = 0.0
    noahmp%water%flux%IrrigationRateSprinkler      = 0.0
    noahmp%water%flux%IrriEvapLossSprinkler        = 0.0
    noahmp%water%flux%EvapIrriSprinkler            = 0.0
    noahmp%water%flux%TileDrain                    = 0.0

    if ( .not. allocated(noahmp%water%flux%CompactionSnowAging) )   &
       allocate( noahmp%water%flux%CompactionSnowAging(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%flux%CompactionSnowBurden) )  &
       allocate( noahmp%water%flux%CompactionSnowBurden(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%flux%CompactionSnowMelt) )    &
       allocate( noahmp%water%flux%CompactionSnowMelt(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%flux%CompactionSnowTot) )     &
       allocate( noahmp%water%flux%CompactionSnowTot(-NumSnowLayerMax+1:0) )
    if ( .not. allocated(noahmp%water%flux%TranspWatLossSoil) )     &
       allocate( noahmp%water%flux%TranspWatLossSoil(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%flux%TranspWatLossSoilAcc) )  &
       allocate( noahmp%water%flux%TranspWatLossSoilAcc(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%flux%TranspWatLossSoilMean) ) &
       allocate( noahmp%water%flux%TranspWatLossSoilMean(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%flux%OutflowSnowLayer) )      &
       allocate( noahmp%water%flux%OutflowSnowLayer(-NumSnowLayerMax+1:0) )

    noahmp%water%flux%CompactionSnowAging  (:)     = undefined_real
    noahmp%water%flux%CompactionSnowBurden (:)     = undefined_real
    noahmp%water%flux%CompactionSnowMelt   (:)     = undefined_real
    noahmp%water%flux%CompactionSnowTot    (:)     = undefined_real
    noahmp%water%flux%TranspWatLossSoil    (:)     = undefined_real
    noahmp%water%flux%TranspWatLossSoilAcc (:)     = undefined_real
    noahmp%water%flux%TranspWatLossSoilMean(:)     = undefined_real
    noahmp%water%flux%OutflowSnowLayer     (:)     = undefined_real

    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 )then
       if ( .not. allocated(noahmp%water%flux%SnowFreezeRate) )     &
       allocate( noahmp%water%flux%SnowFreezeRate(-NumSnowLayerMax+1:0) )

       noahmp%water%flux%SnowFreezeRate(:)         = undefined_real
    endif

    ! water parameter variables
    noahmp%water%param%DrainSoilLayerInd           = undefined_int
    noahmp%water%param%TileDrainTubeDepth          = undefined_int
    noahmp%water%param%NumSoilLayerRoot            = undefined_int
    noahmp%water%param%IrriStopDayBfHarvest        = undefined_int
    noahmp%water%param%CanopyLiqHoldCap            = undefined_real
    noahmp%water%param%SnowCompactBurdenFac        = undefined_real
    noahmp%water%param%SnowCompactAgingFac1        = undefined_real
    noahmp%water%param%SnowCompactAgingFac2        = undefined_real
    noahmp%water%param%SnowCompactAgingFac3        = undefined_real
    noahmp%water%param%SnowCompactAgingMax         = undefined_real
    noahmp%water%param%SnowViscosityCoeff          = undefined_real
    noahmp%water%param%SnowCompactmAR24            = undefined_real
    noahmp%water%param%SnowCompactbAR24            = undefined_real
    noahmp%water%param%SnowCompactP1AR24           = undefined_real
    noahmp%water%param%SnowCompactP2AR24           = undefined_real
    noahmp%water%param%SnowCompactP3AR24           = undefined_real
    noahmp%water%param%SnowCoverM1AR25             = undefined_real
    noahmp%water%param%SnowCoverM2AR25             = undefined_real
    noahmp%water%param%SnowCoverFac1AR25           = undefined_real
    noahmp%water%param%SnowCoverFac2AR25           = undefined_real
    noahmp%water%param%BurdenFacUpAR24             = undefined_real
    noahmp%water%param%SnowLiqFracMax              = undefined_real
    noahmp%water%param%SnowLiqHoldCap              = undefined_real
    noahmp%water%param%SnowLiqReleaseFac           = undefined_real
    noahmp%water%param%IrriFloodRateFac            = undefined_real
    noahmp%water%param%IrriMicroRate               = undefined_real
    noahmp%water%param%SoilInfilMaxCoeff           = undefined_real
    noahmp%water%param%SoilImpervFracCoeff         = undefined_real
    noahmp%water%param%InfilFacVic                 = undefined_real
    noahmp%water%param%TensionWatDistrInfl         = undefined_real
    noahmp%water%param%TensionWatDistrShp          = undefined_real
    noahmp%water%param%FreeWatDistrShp             = undefined_real
    noahmp%water%param%InfilHeteroDynVic           = undefined_real
    noahmp%water%param%InfilCapillaryDynVic        = undefined_real
    noahmp%water%param%InfilFacDynVic              = undefined_real
    noahmp%water%param%SoilDrainSlope              = undefined_real
    noahmp%water%param%TileDrainCoeffSp            = undefined_real
    noahmp%water%param%DrainFacSoilWat             = undefined_real
    noahmp%water%param%TileDrainCoeff              = undefined_real
    noahmp%water%param%DrainDepthToImperv          = undefined_real
    noahmp%water%param%LateralWatCondFac           = undefined_real
    noahmp%water%param%TileDrainDepth              = undefined_real
    noahmp%water%param%DrainTubeDist               = undefined_real
    noahmp%water%param%DrainTubeRadius             = undefined_real
    noahmp%water%param%DrainWatDepToImperv         = undefined_real
    noahmp%water%param%RunoffDecayFac              = undefined_real
    noahmp%water%param%BaseflowCoeff               = undefined_real
    noahmp%water%param%GridTopoIndex               = undefined_real
    noahmp%water%param%SoilSfcSatFracMax           = undefined_real
    noahmp%water%param%SpecYieldGw                 = undefined_real
    noahmp%water%param%MicroPoreContent            = undefined_real
    noahmp%water%param%WaterStorageLakeMax         = undefined_real
    noahmp%water%param%SnoWatEqvMaxGlacier         = undefined_real
    noahmp%water%param%SoilConductivityRef         = undefined_real
    noahmp%water%param%SoilInfilFacRef             = undefined_real
    noahmp%water%param%GroundFrzCoeff              = undefined_real
    noahmp%water%param%IrriTriggerLaiMin           = undefined_real
    noahmp%water%param%SoilWatDeficitAllow         = undefined_real
    noahmp%water%param%IrriFloodLossFrac           = undefined_real
    noahmp%water%param%IrriSprinklerRate           = undefined_real
    noahmp%water%param%IrriFracThreshold           = undefined_real
    noahmp%water%param%IrriStopPrecipThr           = undefined_real
    noahmp%water%param%SnowfallDensityMax          = undefined_real
    noahmp%water%param%SnowMassFullCoverOld        = undefined_real
    noahmp%water%param%SoilMatPotentialWilt        = undefined_real
    noahmp%water%param%SnowMeltFac                 = undefined_real
    noahmp%water%param%SnowCoverFac                = undefined_real
    noahmp%water%param%WetlandCapMax               = undefined_real

    if ( .not. allocated(noahmp%water%param%SoilMoistureSat) )        &
       allocate( noahmp%water%param%SoilMoistureSat(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%param%SoilMoistureWilt) )       &
       allocate( noahmp%water%param%SoilMoistureWilt(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%param%SoilMoistureFieldCap) )   &
       allocate( noahmp%water%param%SoilMoistureFieldCap(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%param%SoilMoistureDry) )        &
       allocate( noahmp%water%param%SoilMoistureDry(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%param%SoilWatDiffusivitySat) )  &
       allocate( noahmp%water%param%SoilWatDiffusivitySat(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%param%SoilWatConductivitySat) ) &
       allocate( noahmp%water%param%SoilWatConductivitySat(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%param%SoilExpCoeffB) )          &
       allocate( noahmp%water%param%SoilExpCoeffB(1:NumSoilLayer) )
    if ( .not. allocated(noahmp%water%param%SoilMatPotentialSat) )    &
       allocate( noahmp%water%param%SoilMatPotentialSat(1:NumSoilLayer) )

    noahmp%water%param%SoilMoistureSat       (:)   = undefined_real
    noahmp%water%param%SoilMoistureWilt      (:)   = undefined_real
    noahmp%water%param%SoilMoistureFieldCap  (:)   = undefined_real
    noahmp%water%param%SoilMoistureDry       (:)   = undefined_real
    noahmp%water%param%SoilWatDiffusivitySat (:)   = undefined_real
    noahmp%water%param%SoilWatConductivitySat(:)   = undefined_real
    noahmp%water%param%SoilExpCoeffB         (:)   = undefined_real
    noahmp%water%param%SoilMatPotentialSat   (:)   = undefined_real

    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 )then
       if ( .not. allocated(noahmp%water%param%snowage_tau) )      &
       allocate( noahmp%water%param%snowage_tau(NumDensitySnwAgeSnicar,NumTempGradSnwAgeSnicar,NumTempSnwAgeSnicar) )
       if ( .not. allocated(noahmp%water%param%snowage_kappa) )    &
       allocate( noahmp%water%param%snowage_kappa(NumDensitySnwAgeSnicar,NumTempGradSnwAgeSnicar,NumTempSnwAgeSnicar) )
       if ( .not. allocated(noahmp%water%param%snowage_drdt0) )    &
       allocate( noahmp%water%param%snowage_drdt0(NumDensitySnwAgeSnicar,NumTempGradSnwAgeSnicar,NumTempSnwAgeSnicar) )

       noahmp%water%param%snowage_tau  (:,:,:)     = undefined_real
       noahmp%water%param%snowage_kappa(:,:,:)     = undefined_real
       noahmp%water%param%snowage_drdt0(:,:,:)     = undefined_real
       noahmp%water%param%SnowRadiusMin            = undefined_real
       noahmp%water%param%FreshSnowRadiusMax       = undefined_real
       noahmp%water%param%SnowRadiusRefrz          = undefined_real
       noahmp%water%param%ScavEffMeltScale         = undefined_real
       noahmp%water%param%ScavEffMeltBCphi         = undefined_real
       noahmp%water%param%ScavEffMeltBCpho         = undefined_real
       noahmp%water%param%ScavEffMeltOCphi         = undefined_real
       noahmp%water%param%ScavEffMeltOCpho         = undefined_real
       noahmp%water%param%ScavEffMeltDust1         = undefined_real
       noahmp%water%param%ScavEffMeltDust2         = undefined_real
       noahmp%water%param%ScavEffMeltDust3         = undefined_real
       noahmp%water%param%ScavEffMeltDust4         = undefined_real
       noahmp%water%param%ScavEffMeltDust5         = undefined_real
       noahmp%water%param%SnowRadiusMax            = undefined_real
       noahmp%water%param%SnowWetAgeC1Brun89       = undefined_real
       noahmp%water%param%SnowWetAgeC2Brun89       = undefined_real
       noahmp%water%param%SnowAgeScaleFac          = undefined_real
    endif

    end associate

  end subroutine WaterVarInitDefault

end module WaterVarInitMod
