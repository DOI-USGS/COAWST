module ConfigVarInitMod

!!! Initialize column (1-D) Noah-MP configuration variables
!!! Configuration variables should be first defined in ConfigVarType.F90

! ------------------------ Code history ------------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! --------------------------------------------------------------------------

  use Machine
  use NoahmpVarType

  implicit none

contains

!=== initialize with default values
  subroutine ConfigVarInitDefault(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    ! config namelist variable
    noahmp%config%nmlist%OptDynamicVeg               = undefined_int
    noahmp%config%nmlist%OptRainSnowPartition        = undefined_int
    noahmp%config%nmlist%OptSoilWaterTranspiration   = undefined_int
    noahmp%config%nmlist%OptGroundResistanceEvap     = undefined_int
    noahmp%config%nmlist%OptSurfaceDrag              = undefined_int
    noahmp%config%nmlist%OptStomataResistance        = undefined_int
    noahmp%config%nmlist%OptSnowAlbedo               = undefined_int
    noahmp%config%nmlist%OptSnowCoverGround          = undefined_int
    noahmp%config%nmlist%OptCanopyRadiationTransfer  = undefined_int
    noahmp%config%nmlist%OptSnowSoilTempTime         = undefined_int
    noahmp%config%nmlist%OptSnowThermConduct         = undefined_int
    noahmp%config%nmlist%OptSoilTemperatureBottom    = undefined_int
    noahmp%config%nmlist%OptSoilSupercoolWater       = undefined_int
    noahmp%config%nmlist%OptRunoffSurface            = undefined_int
    noahmp%config%nmlist%OptRunoffSubsurface         = undefined_int
    noahmp%config%nmlist%OptSoilPermeabilityFrozen   = undefined_int
    noahmp%config%nmlist%OptDynVicInfiltration       = undefined_int
    noahmp%config%nmlist%OptTileDrainage             = undefined_int
    noahmp%config%nmlist%OptIrrigation               = undefined_int
    noahmp%config%nmlist%OptIrrigationMethod         = undefined_int
    noahmp%config%nmlist%OptCropModel                = undefined_int
    noahmp%config%nmlist%OptSoilProperty             = undefined_int
    noahmp%config%nmlist%OptPedotransfer             = undefined_int
    noahmp%config%nmlist%OptGlacierTreatment         = undefined_int
    noahmp%config%nmlist%OptSnowCompaction           = undefined_int 
    noahmp%config%nmlist%OptWetlandModel             = undefined_int
    noahmp%config%nmlist%OptSnicarSnowShape          = undefined_int
    noahmp%config%nmlist%OptSnicarRTSolver           = undefined_int
    noahmp%config%nmlist%OptSnicarBandNum            = undefined_int
    noahmp%config%nmlist%OptSnicarSolarSpec          = undefined_int
    noahmp%config%nmlist%OptSnicarSnwOptic           = undefined_int
    noahmp%config%nmlist%OptSnicarDustOptic          = undefined_int
    noahmp%config%nmlist%FlagSnicarSnowBCIntmix      = .false.
    noahmp%config%nmlist%FlagSnicarSnowDustIntmix    = .false.
    noahmp%config%nmlist%FlagSnicarUseAerosol        = .false.
    noahmp%config%nmlist%FlagSnicarUseOC             = .false.
    noahmp%config%nmlist%FlagSnicarAerosolReadTable  = .false.

    ! config domain variable
    noahmp%config%domain%LandUseDataName             = "MODIFIED_IGBP_MODIS_NOAH"
    noahmp%config%domain%FlagUrban                   = .false.
    noahmp%config%domain%FlagCropland                = .false.
    noahmp%config%domain%FlagWetland                 = .false.
    noahmp%config%domain%FlagDynamicCrop             = .false.
    noahmp%config%domain%FlagDynamicVeg              = .false.
    noahmp%config%domain%FlagSoilProcess             = .false.
    noahmp%config%domain%NumSoilTimeStep             = undefined_int
    noahmp%config%domain%NumSnowLayerMax             = undefined_int
    noahmp%config%domain%NumSnowLayerNeg             = undefined_int
    noahmp%config%domain%NumSoilLayer                = undefined_int
    noahmp%config%domain%GridIndexI                  = undefined_int
    noahmp%config%domain%GridIndexJ                  = undefined_int
    noahmp%config%domain%VegType                     = undefined_int
    noahmp%config%domain%CropType                    = undefined_int
    noahmp%config%domain%SurfaceType                 = undefined_int
    noahmp%config%domain%NumSwRadBand                = undefined_int
    noahmp%config%domain%SoilColor                   = undefined_int
    noahmp%config%domain%IndicatorIceSfc             = undefined_int
    noahmp%config%domain%NumCropGrowStage            = undefined_int
    noahmp%config%domain%IndexWaterPoint             = undefined_int
    noahmp%config%domain%IndexBarrenPoint            = undefined_int
    noahmp%config%domain%IndexIcePoint               = undefined_int
    noahmp%config%domain%IndexCropPoint              = undefined_int
    noahmp%config%domain%IndexEBLForest              = undefined_int
    noahmp%config%domain%NumDayInYear                = undefined_int
    noahmp%config%domain%RunoffSlopeType             = undefined_int
    noahmp%config%domain%MainTimeStep                = undefined_real
    noahmp%config%domain%SoilTimeStep                = undefined_real
    noahmp%config%domain%GridSize                    = undefined_real
    noahmp%config%domain%DayJulianInYear             = undefined_real
    noahmp%config%domain%CosSolarZenithAngle         = undefined_real
    noahmp%config%domain%RefHeightAboveSfc           = undefined_real
    noahmp%config%domain%ThicknessAtmosBotLayer      = undefined_real
    noahmp%config%domain%Latitude                    = undefined_real
    noahmp%config%domain%DepthSoilTempBottom         = undefined_real
    noahmp%config%domain%NumTempSnwAgeSnicar         = undefined_int
    noahmp%config%domain%NumTempGradSnwAgeSnicar     = undefined_int
    noahmp%config%domain%NumDensitySnwAgeSnicar      = undefined_int
    noahmp%config%domain%NumSnicarRadBand            = undefined_int
    noahmp%config%domain%NumRadiusSnwMieSnicar       = undefined_int

  end subroutine ConfigVarInitDefault

end module ConfigVarInitMod
