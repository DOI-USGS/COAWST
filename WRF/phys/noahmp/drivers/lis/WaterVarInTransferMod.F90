module WaterVarInTransferMod

!!! Transfer input 2-D NoahmpIO Water variables to 1-D column variable
!!! 1-D variables should be first defined in /src/WaterVarType.F90
!!! 2-D variables should be first defined in NoahmpIOVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType
  use NoahmpVarType
  use PedoTransferSR2006Mod
  use LisNoahmpParamType

  implicit none

contains

!=== initialize with input data or table values

  subroutine WaterVarInTransfer(noahmp, NoahmpIO, LISparam)

    implicit none

    type(noahmp_type),         intent(inout) :: noahmp
    type(NoahmpIO_type),       intent(inout) :: NoahmpIO
    type(LisNoahmpParam_type), intent(in)    :: LISparam   ! lis/noahmp parameter

    ! local variables 
    integer                            :: IndexSoilLayer
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilSand
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilClay
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilOrg

! -------------------------------------------------------------------------
    associate(                                                         &
              I               => noahmp%config%domain%GridIndexI      ,&
              J               => noahmp%config%domain%GridIndexJ      ,&
              NumSnowLayerMax => noahmp%config%domain%NumSnowLayerMax ,&
              NumSoilLayer    => noahmp%config%domain%NumSoilLayer    ,&
              VegType         => noahmp%config%domain%VegType         ,&
              SoilType        => noahmp%config%domain%SoilType        ,&
              FlagUrban       => noahmp%config%domain%FlagUrban       ,&
              RunoffSlopeType => noahmp%config%domain%RunoffSlopeType ,&
              NumSnowLayerNeg => noahmp%config%domain%NumSnowLayerNeg  &
             )
! -------------------------------------------------------------------------

    ! water state variables
    noahmp%water%state%CanopyLiqWater                     = NoahmpIO%CANLIQXY   (I,J)
    noahmp%water%state%CanopyIce                          = NoahmpIO%CANICEXY   (I,J)
    noahmp%water%state%CanopyWetFrac                      = NoahmpIO%FWETXY     (I,J)
    noahmp%water%state%SnowWaterEquiv                     = NoahmpIO%SNOW       (I,J)
    noahmp%water%state%SnowWaterEquivPrev                 = NoahmpIO%SNEQVOXY   (I,J) 
    noahmp%water%state%SnowDepth                          = NoahmpIO%SNOWH      (I,J)
    noahmp%water%state%IrrigationFracFlood                = NoahmpIO%FIFRACT    (I,J)
    noahmp%water%state%IrrigationAmtFlood                 = NoahmpIO%IRWATFI    (I,J)
    noahmp%water%state%IrrigationFracMicro                = NoahmpIO%MIFRACT    (I,J)
    noahmp%water%state%IrrigationAmtMicro                 = NoahmpIO%IRWATMI    (I,J) 
    noahmp%water%state%IrrigationFracSprinkler            = NoahmpIO%SIFRACT    (I,J)
    noahmp%water%state%IrrigationAmtSprinkler             = NoahmpIO%IRWATSI    (I,J)  
    noahmp%water%state%WaterTableDepth                    = NoahmpIO%ZWTXY      (I,J) 
    noahmp%water%state%SoilMoistureToWT                   = NoahmpIO%SMCWTDXY   (I,J)
    noahmp%water%state%TileDrainFrac                      = NoahmpIO%TD_FRACTION(I,J)
    noahmp%water%state%WaterStorageAquifer                = NoahmpIO%WAXY       (I,J)
    noahmp%water%state%WaterStorageSoilAqf                = NoahmpIO%WTXY       (I,J)
    noahmp%water%state%WaterStorageLake                   = NoahmpIO%WSLAKEXY   (I,J)
    noahmp%water%state%IrrigationFracGrid                 = NoahmpIO%IRFRACT    (I,J)
    noahmp%water%state%IrrigationCntSprinkler             = NoahmpIO%IRNUMSI    (I,J)     
    noahmp%water%state%IrrigationCntMicro                 = NoahmpIO%IRNUMMI    (I,J)
    noahmp%water%state%IrrigationCntFlood                 = NoahmpIO%IRNUMFI    (I,J)
    noahmp%water%state%SnowIce     (-NumSnowLayerMax+1:0) = NoahmpIO%SNICEXY    (I,-NumSnowLayerMax+1:0,J)
    noahmp%water%state%SnowLiqWater(-NumSnowLayerMax+1:0) = NoahmpIO%SNLIQXY    (I,-NumSnowLayerMax+1:0,J)
    noahmp%water%state%SoilLiqWater      (1:NumSoilLayer) = NoahmpIO%SH2O       (I,1:NumSoilLayer,J)
    noahmp%water%state%SoilMoisture      (1:NumSoilLayer) = NoahmpIO%SMOIS      (I,1:NumSoilLayer,J)    
    noahmp%water%state%SoilMoistureEqui  (1:NumSoilLayer) = NoahmpIO%SMOISEQ    (I,1:NumSoilLayer,J)
    noahmp%water%state%RechargeGwDeepWT                   = 0.0
    noahmp%water%state%RechargeGwShallowWT                = 0.0
    if ( noahmp%config%nmlist%OptWetlandModel > 0 ) then
       noahmp%water%state%SoilSaturateFrac                = NoahmpIO%FSATXY     (I,J)
       noahmp%water%state%WaterStorageWetland             = NoahmpIO%WSURFXY    (I,J)
    endif
#ifdef WRF_HYDRO
    noahmp%water%state%WaterTableHydro                    = NoahmpIO%ZWATBLE2D  (I,J)
    noahmp%water%state%WaterHeadSfc                       = NoahmpIO%sfcheadrt  (I,J)
#endif
    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       noahmp%water%state%SnowRadius  (-NumSnowLayerMax+1:0)       = NoahmpIO%SNRDSXY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassBChydrophi(-NumSnowLayerMax+1:0)     = NoahmpIO%BCPHIXY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassBChydropho(-NumSnowLayerMax+1:0)     = NoahmpIO%BCPHOXY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassOChydrophi(-NumSnowLayerMax+1:0)     = NoahmpIO%OCPHIXY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassOChydropho(-NumSnowLayerMax+1:0)     = NoahmpIO%OCPHOXY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassDust1(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST1XY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassDust2(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST2XY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassDust3(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST3XY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassDust4(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST4XY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassDust5(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST5XY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassConcBChydrophi(-NumSnowLayerMax+1:0) = NoahmpIO%MassConcBCPHIXY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassConcBChydropho(-NumSnowLayerMax+1:0) = NoahmpIO%MassConcBCPHOXY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassConcOChydrophi(-NumSnowLayerMax+1:0) = NoahmpIO%MassConcOCPHIXY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassConcOChydropho(-NumSnowLayerMax+1:0) = NoahmpIO%MassConcOCPHOXY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassConcDust1(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST1XY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassConcDust2(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST2XY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassConcDust3(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST3XY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassConcDust4(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST4XY (I,-NumSnowLayerMax+1:0,J)
       noahmp%water%state%MassConcDust5(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST5XY (I,-NumSnowLayerMax+1:0,J)
    endif

    ! water flux variables
    noahmp%water%flux%EvapSoilSfcLiqAcc                   = NoahmpIO%ACC_QSEVAXY (I,J)
    noahmp%water%flux%SoilSfcInflowAcc                    = NoahmpIO%ACC_QINSURXY(I,J)
    noahmp%water%flux%SfcWaterTotChgAcc                   = NoahmpIO%ACC_DWATERXY(I,J)
    noahmp%water%flux%PrecipTotAcc                        = NoahmpIO%ACC_PRCPXY  (I,J)
    noahmp%water%flux%EvapCanopyNetAcc                    = NoahmpIO%ACC_ECANXY  (I,J)
    noahmp%water%flux%TranspirationAcc                    = NoahmpIO%ACC_ETRANXY (I,J)
    noahmp%water%flux%EvapGroundNetAcc                    = NoahmpIO%ACC_EDIRXY  (I,J)
    noahmp%water%flux%TranspWatLossSoilAcc(1:NumSoilLayer)= NoahmpIO%ACC_ETRANIXY(I,1:NumSoilLayer,J)
    noahmp%water%flux%GlacierExcessFlowAcc                = NoahmpIO%ACC_GLAFLWXY(I,J)
    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       noahmp%water%flux%SnowFreezeRate(-NumSnowLayerMax+1:0) = NoahmpIO%SNFRXY(I,-NumSnowLayerMax+1:0,J)
    endif

    ! water parameter variables
    noahmp%water%param%DrainSoilLayerInd                  = LISparam%DRAIN_LAYER_OPT
    noahmp%water%param%CanopyLiqHoldCap                   = LISparam%CH2OP
    noahmp%water%param%SnowCompactBurdenFac               = LISparam%C2_SNOWCOMPACT
    noahmp%water%param%SnowCompactAgingFac1               = LISparam%C3_SNOWCOMPACT
    noahmp%water%param%SnowCompactAgingFac2               = LISparam%C4_SNOWCOMPACT
    noahmp%water%param%SnowCompactAgingFac3               = LISparam%C5_SNOWCOMPACT
    noahmp%water%param%SnowCompactAgingMax                = LISparam%DM_SNOWCOMPACT
    noahmp%water%param%SnowViscosityCoeff                 = LISparam%ETA0_SNOWCOMPACT
    noahmp%water%param%SnowCompactmAR24                   = LISparam%SNOWCOMPACTm_AR24
    noahmp%water%param%SnowCompactbAR24                   = LISparam%SNOWCOMPACTb_AR24
    noahmp%water%param%SnowCompactP1AR24                  = LISparam%SNOWCOMPACT_P1_AR24
    noahmp%water%param%SnowCompactP2AR24                  = LISparam%SNOWCOMPACT_P2_AR24
    noahmp%water%param%SnowCompactP3AR24                  = LISparam%SNOWCOMPACT_P3_AR24
    noahmp%water%param%BurdenFacUpAR24                    = LISparam%SNOWCOMPACT_Up_AR24
    noahmp%water%param%SnowCoverM1AR25                    = LISparam%SCFm1_AR25
    noahmp%water%param%SnowCoverM2AR25                    = LISparam%SCFm2_AR25
    noahmp%water%param%SnowCoverFac1AR25                  = LISparam%SCfac1_AR25
    noahmp%water%param%SnowCoverFac2AR25                  = LISparam%SCfac2_AR25
    noahmp%water%param%SnowLiqFracMax                     = LISparam%SNLIQMAXFRAC
    noahmp%water%param%SnowLiqHoldCap                     = LISparam%SSI
    noahmp%water%param%SnowLiqReleaseFac                  = LISparam%SNOW_RET_FAC
    noahmp%water%param%IrriFloodRateFac                   = LISparam%FIRTFAC
    noahmp%water%param%IrriMicroRate                      = LISparam%MICIR_RATE
    noahmp%water%param%SoilConductivityRef                = LISparam%REFDK
    noahmp%water%param%SoilInfilFacRef                    = LISparam%REFKDT
    noahmp%water%param%GroundFrzCoeff                     = LISparam%FRZK
    noahmp%water%param%GridTopoIndex                      = LISparam%TIMEAN
    noahmp%water%param%SoilSfcSatFracMax                  = LISparam%FSATMX
    noahmp%water%param%SpecYieldGw                        = LISparam%ROUS
    noahmp%water%param%MicroPoreContent                   = LISparam%CMIC
    noahmp%water%param%WaterStorageLakeMax                = LISparam%WSLMAX
    noahmp%water%param%SnoWatEqvMaxGlacier                = LISparam%SWEMAXGLA
    noahmp%water%param%IrriStopDayBfHarvest               = LISparam%IRR_HAR
    noahmp%water%param%IrriTriggerLaiMin                  = LISparam%IRR_LAI
    noahmp%water%param%SoilWatDeficitAllow                = LISparam%IRR_MAD
    noahmp%water%param%IrriFloodLossFrac                  = LISparam%FILOSS
    noahmp%water%param%IrriSprinklerRate                  = LISparam%SPRIR_RATE
    noahmp%water%param%IrriFracThreshold                  = LISparam%IRR_FRAC
    noahmp%water%param%IrriStopPrecipThr                  = LISparam%IR_RAIN
    noahmp%water%param%SnowfallDensityMax                 = LISparam%SNOWDEN_MAX
    noahmp%water%param%SnowMassFullCoverOld               = LISparam%SWEMX
    noahmp%water%param%SoilMatPotentialWilt               = LISparam%PSIWLT
    noahmp%water%param%SnowMeltFac                        = LISparam%MFSNO
    noahmp%water%param%SnowCoverFac                       = LISparam%SCFFAC
    noahmp%water%param%InfilFacVic                        = LISparam%BVIC
    noahmp%water%param%TensionWatDistrInfl                = LISparam%AXAJ
    noahmp%water%param%TensionWatDistrShp                 = LISparam%BXAJ
    noahmp%water%param%FreeWatDistrShp                    = LISparam%XXAJ
    noahmp%water%param%InfilHeteroDynVic                  = LISparam%BBVIC
    noahmp%water%param%InfilCapillaryDynVic               = LISparam%GDVIC
    noahmp%water%param%InfilFacDynVic                     = LISparam%BDVIC
    noahmp%water%param%TileDrainCoeffSp                   = LISparam%TD_DC
    noahmp%water%param%TileDrainTubeDepth                 = LISparam%TD_DEPTH
    noahmp%water%param%DrainFacSoilWat                    = LISparam%TDSMC_FAC
    noahmp%water%param%TileDrainCoeff                     = LISparam%TD_DCOEF
    noahmp%water%param%DrainDepthToImperv                 = LISparam%TD_ADEPTH
    noahmp%water%param%LateralWatCondFac                  = LISparam%KLAT_FAC
    noahmp%water%param%TileDrainDepth                     = LISparam%TD_DDRAIN
    noahmp%water%param%DrainTubeDist                      = LISparam%TD_SPAC
    noahmp%water%param%DrainTubeRadius                    = LISparam%TD_RADI
    noahmp%water%param%DrainWatDepToImperv                = LISparam%TD_D
    noahmp%water%param%NumSoilLayerRoot                   = LISparam%NROOT
    noahmp%water%param%SoilDrainSlope                     = LISparam%SLOPE
    noahmp%water%param%WetlandCapMax                      = LISparam%WCAP

    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 )then
       noahmp%water%param%snowage_tau                     = NoahmpIO%snowage_tau
       noahmp%water%param%snowage_kappa                   = NoahmpIO%snowage_kappa
       noahmp%water%param%snowage_drdt0                   = NoahmpIO%snowage_drdt0
       noahmp%water%param%SnowRadiusMin                   = LISparam%SnowRadiusMin
       noahmp%water%param%FreshSnowRadiusMax              = LISparam%FreshSnowRadiusMax
       noahmp%water%param%SnowRadiusRefrz                 = LISparam%SnowRadiusRefrz
       noahmp%water%param%ScavEffMeltScale                = LISparam%ScavEffMeltScale
       noahmp%water%param%ScavEffMeltBCphi                = LISparam%ScavEffMeltBCphi
       noahmp%water%param%ScavEffMeltBCpho                = LISparam%ScavEffMeltBCpho
       noahmp%water%param%ScavEffMeltOCphi                = LISparam%ScavEffMeltOCphi
       noahmp%water%param%ScavEffMeltOCpho                = LISparam%ScavEffMeltOCpho
       noahmp%water%param%ScavEffMeltDust1                = LISparam%ScavEffMeltDust1
       noahmp%water%param%ScavEffMeltDust2                = LISparam%ScavEffMeltDust2
       noahmp%water%param%ScavEffMeltDust3                = LISparam%ScavEffMeltDust3
       noahmp%water%param%ScavEffMeltDust4                = LISparam%ScavEffMeltDust4
       noahmp%water%param%ScavEffMeltDust5                = LISparam%ScavEffMeltDust5
       noahmp%water%param%SnowRadiusMax                   = LISparam%SnowRadiusMax
       noahmp%water%param%SnowWetAgeC1Brun89              = LISparam%SnowWetAgeC1Brun89
       noahmp%water%param%SnowWetAgeC2Brun89              = LISparam%SnowWetAgeC2Brun89
       noahmp%water%param%SnowAgeScaleFac                 = LISparam%SnowAgeScaleFac
    endif

    ! soil properties
    do IndexSoilLayer = 1, size(SoilType)
       noahmp%water%param%SoilMoistureSat       (IndexSoilLayer) = LISparam%SMCMAX(IndexSoilLayer)
       noahmp%water%param%SoilMoistureWilt      (IndexSoilLayer) = LISparam%SMCWLT(IndexSoilLayer)
       noahmp%water%param%SoilMoistureFieldCap  (IndexSoilLayer) = LISparam%SMCREF(IndexSoilLayer)
       noahmp%water%param%SoilMoistureDry       (IndexSoilLayer) = LISparam%SMCDRY(IndexSoilLayer)
       noahmp%water%param%SoilWatDiffusivitySat (IndexSoilLayer) = LISparam%DWSAT(IndexSoilLayer)
       noahmp%water%param%SoilWatConductivitySat(IndexSoilLayer) = LISparam%DKSAT(IndexSoilLayer)
       noahmp%water%param%SoilExpCoeffB         (IndexSoilLayer) = LISparam%BEXP(IndexSoilLayer)
       noahmp%water%param%SoilMatPotentialSat   (IndexSoilLayer) = LISparam%PSISAT(IndexSoilLayer)
    enddo
   
    ! spatial varying soil texture and properties directly from input
    if ( noahmp%config%nmlist%OptSoilProperty == 4 ) then
       ! 3D soil properties
       noahmp%water%param%SoilExpCoeffB          = NoahmpIO%BEXP_3D  (I,1:NumSoilLayer,J) ! C-H B exponent
       noahmp%water%param%SoilMoistureDry        = NoahmpIO%SMCDRY_3D(I,1:NumSoilLayer,J) ! Soil Moisture Limit: Dry
       noahmp%water%param%SoilMoistureWilt       = NoahmpIO%SMCWLT_3D(I,1:NumSoilLayer,J) ! Soil Moisture Limit: Wilt
       noahmp%water%param%SoilMoistureFieldCap   = NoahmpIO%SMCREF_3D(I,1:NumSoilLayer,J) ! Soil Moisture Limit: Reference
       noahmp%water%param%SoilMoistureSat        = NoahmpIO%SMCMAX_3D(I,1:NumSoilLayer,J) ! Soil Moisture Limit: Max
       noahmp%water%param%SoilWatConductivitySat = NoahmpIO%DKSAT_3D (I,1:NumSoilLayer,J) ! Saturated Soil Conductivity
       noahmp%water%param%SoilWatDiffusivitySat  = NoahmpIO%DWSAT_3D (I,1:NumSoilLayer,J) ! Saturated Soil Diffusivity
       noahmp%water%param%SoilMatPotentialSat    = NoahmpIO%PSISAT_3D(I,1:NumSoilLayer,J) ! Saturated Matric Potential
       noahmp%water%param%SoilConductivityRef    = NoahmpIO%REFDK_2D (I,J)                ! Reference Soil Conductivity
       noahmp%water%param%SoilInfilFacRef        = NoahmpIO%REFKDT_2D(I,J)                ! Soil Infiltration Parameter
       ! 2D additional runoff6~8 parameters
       noahmp%water%param%InfilFacVic            = NoahmpIO%BVIC_2D (I,J)                 ! VIC model infiltration parameter
       noahmp%water%param%TensionWatDistrInfl    = NoahmpIO%AXAJ_2D (I,J)                 ! Xinanjiang: Tension water distribution inflection parameter
       noahmp%water%param%TensionWatDistrShp     = NoahmpIO%BXAJ_2D (I,J)                 ! Xinanjiang: Tension water distribution shape parameter
       noahmp%water%param%FreeWatDistrShp        = NoahmpIO%XXAJ_2D (I,J)                 ! Xinanjiang: Free water distribution shape parameter
       noahmp%water%param%InfilFacDynVic         = NoahmpIO%BDVIC_2D(I,J)                 ! VIC model infiltration parameter
       noahmp%water%param%InfilCapillaryDynVic   = NoahmpIO%GDVIC_2D(I,J)                 ! Mean Capillary Drive for infiltration models
       noahmp%water%param%InfilHeteroDynVic      = NoahmpIO%BBVIC_2D(I,J)                 ! DVIC heterogeniety parameter for infiltraton
       ! 2D irrigation params
       noahmp%water%param%IrriFracThreshold      = NoahmpIO%IRR_FRAC_2D  (I,J)            ! irrigation Fraction
       noahmp%water%param%IrriStopDayBfHarvest   = NoahmpIO%IRR_HAR_2D   (I,J)            ! number of days before harvest date to stop irrigation 
       noahmp%water%param%IrriTriggerLaiMin      = NoahmpIO%IRR_LAI_2D   (I,J)            ! Minimum lai to trigger irrigation
       noahmp%water%param%SoilWatDeficitAllow    = NoahmpIO%IRR_MAD_2D   (I,J)            ! management allowable deficit (0-1)
       noahmp%water%param%IrriFloodLossFrac      = NoahmpIO%FILOSS_2D    (I,J)            ! fraction of flood irrigation loss (0-1) 
       noahmp%water%param%IrriSprinklerRate      = NoahmpIO%SPRIR_RATE_2D(I,J)            ! mm/h, sprinkler irrigation rate
       noahmp%water%param%IrriMicroRate          = NoahmpIO%MICIR_RATE_2D(I,J)            ! mm/h, micro irrigation rate
       noahmp%water%param%IrriFloodRateFac       = NoahmpIO%FIRTFAC_2D   (I,J)            ! flood application rate factor
       noahmp%water%param%IrriStopPrecipThr      = NoahmpIO%IR_RAIN_2D   (I,J)            ! maximum precipitation to stop irrigation trigger
       ! 2D tile drainage parameters
       noahmp%water%param%LateralWatCondFac      = NoahmpIO%KLAT_FAC (I,J)                ! factor multiplier to hydraulic conductivity
       noahmp%water%param%DrainFacSoilWat        = NoahmpIO%TDSMC_FAC(I,J)                ! factor multiplier to field capacity
       noahmp%water%param%TileDrainCoeffSp       = NoahmpIO%TD_DC    (I,J)                ! drainage coefficient for simple
       noahmp%water%param%TileDrainCoeff         = NoahmpIO%TD_DCOEF (I,J)                ! drainge coefficient for Hooghoudt 
       noahmp%water%param%TileDrainDepth         = NoahmpIO%TD_DDRAIN(I,J)                ! depth of drain
       noahmp%water%param%DrainTubeRadius        = NoahmpIO%TD_RADI  (I,J)                ! tile tube radius
       noahmp%water%param%DrainTubeDist          = NoahmpIO%TD_SPAC  (I,J)                ! tile spacing
    endif

    ! spatial varying wetland parameters from input
    if ( noahmp%config%nmlist%OptWetlandModel == 2 ) then
       noahmp%water%param%SoilSfcSatFracMax      = NoahmpIO%FSATMX(I,J)
       noahmp%water%param%WetlandCapMax          = NoahmpIO%WCAP(I,J)
    endif

    ! derived water parameters
    noahmp%water%param%SoilInfilMaxCoeff  = noahmp%water%param%SoilInfilFacRef *           &
                                            noahmp%water%param%SoilWatConductivitySat(1) / &
                                            noahmp%water%param%SoilConductivityRef
    if ( FlagUrban .eqv. .true. ) then
       noahmp%water%param%SoilMoistureSat      = 0.45
       noahmp%water%param%SoilMoistureFieldCap = 0.42
       noahmp%water%param%SoilMoistureWilt     = 0.40
       noahmp%water%param%SoilMoistureDry      = 0.40
    endif

    if ( SoilType(1) /= 14 ) then
       noahmp%water%param%SoilImpervFracCoeff = noahmp%water%param%GroundFrzCoeff *       &
                                                ((noahmp%water%param%SoilMoistureSat(1) / &
                                                 noahmp%water%param%SoilMoistureFieldCap(1)) * (0.412/0.468))
    endif

    noahmp%water%state%SnowIceFracPrev = 0.0
    noahmp%water%state%SnowIceFracPrev(NumSnowLayerNeg+1:0) = NoahmpIO%SNICEXY(I,NumSnowLayerNeg+1:0,J) /  & 
                                                              (NoahmpIO%SNICEXY(I,NumSnowLayerNeg+1:0,J) + &
                                                               NoahmpIO%SNLIQXY(I,NumSnowLayerNeg+1:0,J))

    if ( (noahmp%config%nmlist%OptSoilProperty == 3) .and. (.not. noahmp%config%domain%FlagUrban) ) then
       if (.not. allocated(SoilSand)) allocate( SoilSand(1:NumSoilLayer) )
       if (.not. allocated(SoilClay)) allocate( SoilClay(1:NumSoilLayer) )
       if (.not. allocated(SoilOrg) ) allocate( SoilOrg (1:NumSoilLayer) )
       SoilSand = 0.01 * NoahmpIO%soilcomp(I,1:NumSoilLayer,J)
       SoilClay = 0.01 * NoahmpIO%soilcomp(I,(NumSoilLayer+1):(NumSoilLayer*2),J)
       SoilOrg  = 0.0
       if (noahmp%config%nmlist%OptPedotransfer == 1) &
          call PedoTransferSR2006(NoahmpIO,noahmp,SoilSand,SoilClay,SoilOrg)
       deallocate(SoilSand)
       deallocate(SoilClay)
       deallocate(SoilOrg )
    endif

    end associate

  end subroutine WaterVarInTransfer

end module WaterVarInTransferMod
