module ConfigVarType

!!! Define column (1-D) Noah-MP configuration variables
!!! Configuration variable initialization is done in ConfigVarInitMod.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine

  implicit none
  save
  private

!=== define "namelist" sub-type of config (config%nmlist%variable)
  type :: namelist_type

    integer :: OptDynamicVeg               ! options for dynamic vegetation
                                              ! 1 -> off (use table LeafAreaIndex; use VegFrac = VegFracGreen from input)
                                              ! 2 -> on  (together with OptStomataResistance = 1)
                                              ! 3 -> off (use table LeafAreaIndex; calculate VegFrac)
                                              ! 4 -> off (use table LeafAreaIndex; use maximum vegetation fraction)
                                              ! 5 -> on  (use maximum vegetation fraction)
                                              ! 6 -> on  (use VegFrac = VegFracGreen from input)
                                              ! 7 -> off (use input LeafAreaIndex; use VegFrac = VegFracGreen from input)
                                              ! 8 -> off (use input LeafAreaIndex; calculate VegFrac)
                                              ! 9 -> off (use input LeafAreaIndex; use maximum vegetation fraction)
    integer :: OptRainSnowPartition        ! options for partitioning  precipitation into rainfall & snowfall
                                              ! 1 -> Jordan (1991) scheme
                                              ! 2 -> BATS: when TemperatureAirRefHeight < freezing point+2.2 
                                              ! 3 -> TemperatureAirRefHeight < freezing point
                                              ! 4 -> Use WRF microphysics output
                                              ! 5 -> Use wetbulb temperature (Wang et al., 2019)
    integer :: OptSoilWaterTranspiration   ! options for soil moisture factor for stomatal resistance & evapotranspiration
                                              ! 1 -> Noah (soil moisture)
                                              ! 2 -> CLM  (matric potential)
                                              ! 3 -> SSiB (matric potential)
    integer :: OptGroundResistanceEvap     ! options for ground resistent to evaporation/sublimation
                                              ! 1 -> Sakaguchi and Zeng, 2009
                                              ! 2 -> Sellers (1992)
                                              ! 3 -> adjusted Sellers to decrease ResistanceGrdEvap for wet soil
                                              ! 4 -> option 1 for non-snow; rsurf = rsurf_snow for snow (set in table)
    integer :: OptSurfaceDrag              ! options for surface layer drag/exchange coefficient
                                              ! 1 -> Monin-Obukhov (M-O) Similarity Theory (MOST)
                                              ! 2 -> original Noah (Chen et al. 1997)
    integer :: OptStomataResistance        ! options for canopy stomatal resistance
                                              ! 1 -> Ball-Berry scheme
                                              ! 2 -> Jarvis scheme
    integer :: OptSnowAlbedo               ! options for ground snow surface albedo
                                              ! 1 -> BATS snow albedo scheme
                                              ! 2 -> CLASS snow albedo scheme
                                              ! 3 -> SNICAR snow albedo scheme (Lin et al., 2025 JHM)
    integer :: OptCanopyRadiationTransfer  ! options for canopy radiation transfer
                                              ! 1 -> modified two-stream (gap=F(solar angle,3D structure, etc)<1-VegFrac)
                                              ! 2 -> two-stream applied to grid-cell (gap = 0)
                                              ! 3 -> two-stream applied to vegetated fraction (gap=1-VegFrac)
    integer :: OptSnowSoilTempTime         ! options for snow/soil temperature time scheme (only layer 1)
                                              ! 1 -> semi-implicit; flux top boundary condition
                                              ! 2 -> full implicit (original Noah); temperature top boundary condition
                                              ! 3 -> same as 1, but snow cover for skin temperature calculation (generally improves snow)
    integer :: OptSnowThermConduct         ! options for snow thermal conductivity
                                              ! 1 -> Stieglitz(yen,1965) scheme
                                              ! 2 -> Anderson, 1976 scheme
                                              ! 3 -> constant
                                              ! 4 -> Verseghy (1991) scheme
                                              ! 5 -> Douvill(Yen, 1981) scheme
    integer :: OptSoilTemperatureBottom    ! options for lower boundary condition of soil temperature
                                              ! 1 -> zero heat flux from bottom (DepthSoilTempBottom & TemperatureSoilBottom not used)
                                              ! 2 -> TemperatureSoilBottom at DepthSoilTempBottom (8m) read from a file (original Noah)
    integer :: OptSoilSupercoolWater       ! options for soil supercooled liquid water
                                              ! 1 -> no iteration (Niu and Yang, 2006 JHM)
                                              ! 2 -> Koren's iteration (Koren et al., 1999 JGR)
    integer :: OptRunoffSurface            ! options for surface runoff
                                              ! 1 -> TOPMODEL with groundwater
                                              ! 2 -> TOPMODEL with an equilibrium water table
                                              ! 3 -> original surface and subsurface runoff (free drainage)
                                              ! 4 -> BATS surface and subsurface runoff (free drainage)
                                              ! 5 -> Miguez-Macho&Fan groundwater scheme
                                              ! 6 -> Variable Infiltration Capacity Model surface runoff scheme
                                              ! 7 -> Xinanjiang Infiltration and surface runoff scheme 
                                              ! 8 -> Dynamic VIC surface runoff scheme
    integer :: OptRunoffSubsurface         ! options for drainage & subsurface runoff 
                                              ! 1~8: similar to runoff option, separated from original NoahMP runoff option
                                              ! currently tested & recommended the same option# as surface runoff
    integer :: OptSoilPermeabilityFrozen   ! options for frozen soil permeability
                                              ! 1 -> linear effects, more permeable
                                              ! 2 -> nonlinear effects, less permeable
    integer :: OptDynVicInfiltration       ! options for infiltration in dynamic VIC runoff scheme
                                              ! 1 -> Philip scheme
                                              ! 2 -> Green-Ampt scheme 
                                              ! 3 -> Smith-Parlange scheme    
    integer :: OptTileDrainage             ! options for tile drainage 
                                              ! currently only tested & calibrated to work with runoff option=3
                                              ! 0 -> No tile drainage
                                              ! 1 -> on (simple scheme)
                                              ! 2 -> on (Hooghoudt's scheme)
    integer :: OptIrrigation               ! options for irrigation
                                              ! 0 -> No irrigation
                                              ! 1 -> Irrigation ON
                                              ! 2 -> irrigation trigger based on crop season Planting and harvesting dates
                                              ! 3 -> irrigation trigger based on LeafAreaIndex threshold
    integer :: OptIrrigationMethod         ! options for irrigation method
                                              ! only works when OptIrrigation > 0
                                              ! 0 -> method based on geo_em fractions
                                              ! 1 -> sprinkler method
                                              ! 2 -> micro/drip irrigation
                                              ! 3 -> surface flooding
    integer :: OptCropModel                ! options for crop model
                                              ! 0 -> No crop model
                                              ! 1 -> Liu, et al. 2016 crop scheme
    integer :: OptSoilProperty             ! options for defining soil properties
                                              ! 1 -> use input dominant soil texture
                                              ! 2 -> use input soil texture that varies with depth
                                              ! 3 -> use soil composition (sand, clay, orgm) and pedotransfer function
                                              ! 4 -> use input soil properties
    integer :: OptPedotransfer             ! options for pedotransfer functions 
                                              ! only works when OptSoilProperty = 3
                                              ! 1 -> Saxton and Rawls (2006) scheme
    integer :: OptGlacierTreatment         ! options for glacier treatment
                                              ! 1 -> include phase change of ice
                                              ! 2 -> ice treatment more like original Noah
    integer :: OptSnowCompaction           ! options for ground snow compaction
                                              ! 1 -> original scheme from Anderson (1976)
                                              ! 2 -> new scheme from Abolafia-Rosenzweig et al. (2024)
    integer :: OptSnowCoverGround          ! options for ground snow cover fraction
                                              ! 1 -> original scheme from Niu and Yang (07) using veg-class based MPTABLE parameters
                                              ! 2 -> enhanced scheme from Abolafia-Rosenzweig et al. (2025) adding scale-dependency to ground SCF parameters
    integer :: OptWetlandModel             ! option for wetland model
                                              ! 0 -> No Wetland model
                                              ! 1 -> Single-point/uniform parameter (Zhang, et al. 2022 WRR)
                                              ! 2 -> 2-D regional parameter input (Zhang, et al. 2022 WRR)
    integer :: OptSnicarSnowShape          ! options for snow grain shape in SNICAR (He et al. 2017 JC)
                                              ! 1 -> sphere
                                              ! 2 -> spheroid
                                              ! 3 -> hexagonal plate
                                              ! 4 -> Koch snowflake
    integer :: OptSnicarRTSolver           ! option for two different SNICAR radiative transfer solver
                                              ! 1 -> Toon et a 1989 2-stream (Flanner et al. 2007)
                                              ! 2 -> Adding-doubling 2-stream (Dang et al.2019)
    integer :: OptSnicarBandNum            ! option for SNICAR number of solar bands in RT solver
                                              ! 1 -> 5 bands
                                              ! 2 -> 480 bands (10-nm spectral resolution)   
    integer :: OptSnicarSolarSpec          ! option for SNICAR downward solar spectrum 
                                              ! 1 -> mid-latitude winter
                                              ! 2 -> mid-latitude summer
                                              ! 3 -> sub-Arctic winter
                                              ! 4 -> sub-Arctic summer
                                              ! 5 -> Summit,Greenland,summer
                                              ! 6 -> High Mountain summer
    integer :: OptSnicarSnwOptic           ! option for snow optics using different refractive index databases in SNICAR
                                              ! 1 -> Warren (1984)
                                              ! 2 -> Warren and Brandt (2008)
                                              ! 3 -> Picard et al (2016)
    integer :: OptSnicarDustOptic          ! option for dust optics for SNICAR snow albedo calculation
                                              ! 1 -> Saharan dust (Balkanski et al., 2007, central hematite)
                                              ! 2 -> San Juan Mountains dust, CO (Skiles et al, 2017)
                                              ! 3 -> Greenland dust (Polashenski et al., 2015, central absorptivity)
    logical :: FlagSnicarSnowBCIntmix      ! flag to determine SNICAR BC-snow mixing state
                                              ! .false. -> external mixing for all BC
                                              ! .true.  -> internal mixing for hydrophilic BC
    logical :: FlagSnicarSnowDustIntmix    ! flag to determine SNICAR dust-snow mixing state
                                              ! .false. -> external mixing for all dust
                                              ! .true.  -> internal mixing for all dust
    logical :: FlagSnicarUseAerosol        ! option to turn on/off aerosol deposition flux effect in snow in SNICAR
                                              ! .false. -> without aerosol deposition flux effect
                                              ! .true.  -> with aerosol deposition flux effect
    logical :: FlagSnicarUseOC             ! option to activate OC in snow in SNICAR
                                              ! .false. -> without organic carbon in snow
                                              ! .true.  -> with organic carbon in snow
    logical :: FlagSnicarAerosolReadTable  ! option to read aerosol deposition fluxes from table or not
                                              ! .false. -> data read from NetCDF forcing file
                                              ! .true.  -> data read from table
   end type namelist_type


!=== define "domain" sub-type of config (config%domain%variable)
  type :: domain_type

    character(len=256)     :: LandUseDataName             ! landuse dataset name (USGS or MODIFIED_IGBP_MODIS_NOAH)
    logical                :: FlagUrban                   ! flag for urban grid
    logical                :: FlagCropland                ! flag to identify croplands
    logical                :: FlagWetland                 ! flag to identify wetlands
    logical                :: FlagDynamicCrop             ! flag to activate dynamic crop model
    logical                :: FlagDynamicVeg              ! flag to activate dynamic vegetation scheme
    logical                :: FlagSoilProcess             ! flag to determine if calculating soil processes
    integer                :: GridIndexI                  ! model grid index in x-direction
    integer                :: GridIndexJ                  ! model grid index in y-direction
    integer                :: VegType                     ! vegetation type
    integer                :: CropType                    ! crop type
    integer                :: NumSoilLayer                ! number of soil layers
    integer                :: NumSnowLayerMax             ! maximum number of snow layers
    integer                :: NumSnowLayerNeg             ! actual number of snow layers (negative)
    integer                :: SurfaceType                 ! surface type (1=soil; 2=lake)
    integer                :: NumSwRadBand                ! number of shortwave radiation bands
    integer                :: SoilColor                   ! soil color type for albedo
    integer                :: IndicatorIceSfc             ! indicator for ice surface/point (1=sea ice, 0=non-ice, -1=land ice)
    integer                :: IndexWaterPoint             ! land type index for water point
    integer                :: IndexBarrenPoint            ! land type index for barren land point
    integer                :: IndexIcePoint               ! land type index for  ice point
    integer                :: IndexCropPoint              ! land type index for cropland point
    integer                :: IndexEBLForest              ! land type index for evergreen broadleaf (EBL) Forest
    integer                :: NumCropGrowStage            ! number of crop growth stages
    integer                :: NumDayInYear                ! Number of days in the particular year
    integer                :: RunoffSlopeType             ! underground runoff slope term type
    integer                :: NumSoilTimeStep             ! number of timesteps to calculate soil processes
    integer                :: NumTempSnwAgeSnicar         ! maxiumum temperature index used in aging lookup table [idx]
    integer                :: NumTempGradSnwAgeSnicar     ! maxiumum temperature gradient index used in aging lookup table [idx]
    integer                :: NumDensitySnwAgeSnicar      ! maxiumum snow density index used in aging lookup table [idx]
    integer                :: NumSnicarRadBand            ! wavelength bands used in SNICAR snow albedo calculation
    integer                :: NumRadiusSnwMieSnicar       ! number of effective radius indices used in Mie lookup table [idx] 
    real(kind=kind_noahmp) :: MainTimeStep                ! noahmp main timestep [sec]
    real(kind=kind_noahmp) :: SoilTimeStep                ! soil timestep [sec]
    real(kind=kind_noahmp) :: GridSize                    ! noahmp model grid spacing [m]
    real(kind=kind_noahmp) :: DayJulianInYear             ! julian day of the year
    real(kind=kind_noahmp) :: CosSolarZenithAngle         ! cosine solar zenith angle
    real(kind=kind_noahmp) :: RefHeightAboveSfc           ! reference height [m] above surface zero plane (including vegetation)
    real(kind=kind_noahmp) :: ThicknessAtmosBotLayer      ! thickness of atmospheric bottom layers [m]
    real(kind=kind_noahmp) :: Latitude                    ! latitude [degree]
    real(kind=kind_noahmp) :: DepthSoilTempBottom         ! depth [m, negative] from soil surface for lower boundary soil temperature forcing

    integer               , allocatable, dimension(:) :: SoilType                  ! soil type for each soil layer
    real(kind=kind_noahmp), allocatable, dimension(:) :: DepthSoilLayer            ! depth [m] of layer-bottom from soil surface
    real(kind=kind_noahmp), allocatable, dimension(:) :: ThicknessSnowSoilLayer    ! snow and soil layer thickness [m]
    real(kind=kind_noahmp), allocatable, dimension(:) :: DepthSnowSoilLayer        ! snow and soil layer-bottom depth [m]
    real(kind=kind_noahmp), allocatable, dimension(:) :: ThicknessSoilLayer        ! soil layer thickness [m]

  end type domain_type


!=== define config type that includes namelist & domain subtypes
  type, public :: config_type

    type(namelist_type) :: nmlist
    type(domain_type)   :: domain

  end type config_type

end module ConfigVarType
