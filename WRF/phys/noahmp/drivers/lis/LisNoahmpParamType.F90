module LisNoahmpParamType

!!! Define LIS-Noah-MP table parameter variables

! ------------------------ Code history -----------------------------------
! Refactered code: C. He, P. Valayamkunnath & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine

  implicit none
  save
  private

  integer, private, parameter :: MBAND  = 2
  integer, private, parameter :: NSOIL  = 4
  integer, private, parameter :: NSTAGE = 8

  type, public :: LisNoahmpParam_type

!----------------------------------------------------------------
! Noahmp Parameters Table 
!----------------------------------------------------------------

    ! vegetation parameters
    logical                     :: URBAN_FLAG          ! urban flag
    integer                     :: ISWATER             ! water flag
    integer                     :: ISBARREN            ! barren ground flag
    integer                     :: ISICE               ! ice flag
    integer                     :: ISCROP              ! cropland flag
    integer                     :: EBLFOREST           ! evergreen broadleaf forest flag
    real(kind=kind_noahmp)      :: CH2OP               ! maximum intercepted h2o per unit lai+sai (mm)
    real(kind=kind_noahmp)      :: DLEAF               ! characteristic leaf dimension (m)
    real(kind=kind_noahmp)      :: Z0MVT               ! momentum roughness length (m)
    real(kind=kind_noahmp)      :: HVT                 ! top of canopy (m)
    real(kind=kind_noahmp)      :: HVB                 ! bottom of canopy (m)
    real(kind=kind_noahmp)      :: DEN                 ! tree density (no. of trunks per m2)
    real(kind=kind_noahmp)      :: RC                  ! tree crown radius (m)
    real(kind=kind_noahmp)      :: MFSNO               ! snowmelt curve parameter
    real(kind=kind_noahmp)      :: SCFFAC              ! snow cover factor (m) (replace original hard-coded 2.5*z0 in SCF formulation)
    real(kind=kind_noahmp)      :: CBIOM               ! canopy biomass heat capacity parameter (m) 
    real(kind=kind_noahmp)      :: SAIM(12)            ! monthly stem area index, one-sided
    real(kind=kind_noahmp)      :: LAIM(12)            ! monthly leaf area index, one-sided
    real(kind=kind_noahmp)      :: SLA                 ! single-side leaf area per Kg [m2/kg]
    real(kind=kind_noahmp)      :: DILEFC              ! coeficient for leaf stress death [1/s]
    real(kind=kind_noahmp)      :: DILEFW              ! coeficient for leaf stress death [1/s]
    real(kind=kind_noahmp)      :: FRAGR               ! fraction of growth respiration  !original was 0.3 
    real(kind=kind_noahmp)      :: LTOVRC              ! leaf turnover [1/s]
    real(kind=kind_noahmp)      :: C3PSN               ! photosynthetic pathway: 0. = c4, 1. = c3
    real(kind=kind_noahmp)      :: KC25                ! co2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp)      :: AKC                 ! q10 for kc25
    real(kind=kind_noahmp)      :: KO25                ! o2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp)      :: AKO                 ! q10 for ko25
    real(kind=kind_noahmp)      :: VCMX25              ! maximum rate of carboxylation at 25c (umol co2/m2/s)
    real(kind=kind_noahmp)      :: AVCMX               ! q10 for vcmx25
    real(kind=kind_noahmp)      :: BP                  ! minimum leaf conductance (umol/m2/s)
    real(kind=kind_noahmp)      :: MP                  ! slope of conductance-to-photosynthesis relationship
    real(kind=kind_noahmp)      :: QE25                ! quantum efficiency at 25c (umol co2 / umol photon)
    real(kind=kind_noahmp)      :: AQE                 ! q10 for qe25
    real(kind=kind_noahmp)      :: RMF25               ! leaf maintenance respiration at 25c (umol co2/m2/s)
    real(kind=kind_noahmp)      :: RMS25               ! stem maintenance respiration at 25c (umol co2/kg bio/s)
    real(kind=kind_noahmp)      :: RMR25               ! root maintenance respiration at 25c (umol co2/kg bio/s)
    real(kind=kind_noahmp)      :: ARM                 ! q10 for maintenance respiration
    real(kind=kind_noahmp)      :: FOLNMX              ! foliage nitrogen concentration when f(n)=1 (%)
    real(kind=kind_noahmp)      :: TMIN                ! minimum temperature for photosynthesis (k)
    real(kind=kind_noahmp)      :: XL                  ! leaf/stem orientation index
    real(kind=kind_noahmp)      :: RHOL(MBAND)         ! leaf reflectance: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: RHOS(MBAND)         ! stem reflectance: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: TAUL(MBAND)         ! leaf transmittance: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: TAUS(MBAND)         ! stem transmittance: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: MRP                 ! microbial respiration parameter (umol co2 /kg c/ s)
    real(kind=kind_noahmp)      :: CWPVT               ! empirical canopy wind parameter
    real(kind=kind_noahmp)      :: WRRAT               ! wood to non-wood ratio
    real(kind=kind_noahmp)      :: WDPOOL              ! wood pool (switch 1 or 0) depending on woody or not [-]
    real(kind=kind_noahmp)      :: TDLEF               ! characteristic T for leaf freezing [K]
    real(kind=kind_noahmp)      :: NROOT               ! number of soil layers with root present
    real(kind=kind_noahmp)      :: RGL                 ! Parameter used in radiation stress function
    real(kind=kind_noahmp)      :: RSMIN               ! Minimum stomatal resistance [s m-1]
    real(kind=kind_noahmp)      :: HS                  ! Parameter used in vapor pressure deficit function
    real(kind=kind_noahmp)      :: TOPT                ! Optimum transpiration air temperature [K]
    real(kind=kind_noahmp)      :: RSMAX               ! Maximal stomatal resistance [s m-1]
    real(kind=kind_noahmp)      :: RTOVRC              ! root turnover coefficient [1/s]
    real(kind=kind_noahmp)      :: RSWOODC             ! wood respiration coeficient [1/s]
    real(kind=kind_noahmp)      :: BF                  ! parameter for present wood allocation [-]
    real(kind=kind_noahmp)      :: WSTRC               ! water stress coeficient [-]
    real(kind=kind_noahmp)      :: LAIMIN              ! minimum leaf area index [m2/m2]
    real(kind=kind_noahmp)      :: XSAMIN              ! minimum stem area index [m2/m2]

    ! radiation parameters
    real(kind=kind_noahmp)      :: ALBSAT(MBAND)       ! saturated soil albedos: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: ALBDRY(MBAND)       ! dry soil albedos: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: ALBICE(MBAND)       ! albedo land ice: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: ALBLAK(MBAND)       ! albedo frozen lakes: 1=vis, 2=nir
    real(kind=kind_noahmp)      :: OMEGAS(MBAND)       ! two-stream parameter omega for snow
    real(kind=kind_noahmp)      :: BETADS              ! two-stream parameter betad for snow
    real(kind=kind_noahmp)      :: BETAIS              ! two-stream parameter betad for snow
    real(kind=kind_noahmp)      :: EG(2)               ! emissivity soil surface
    real(kind=kind_noahmp)      :: EICE                ! ice surface emissivity

    ! global parameters
    real(kind=kind_noahmp)      :: CO2                 ! co2 partial pressure
    real(kind=kind_noahmp)      :: O2                  ! o2 partial pressure
    real(kind=kind_noahmp)      :: TIMEAN              ! gridcell mean topgraphic index (global mean)
    real(kind=kind_noahmp)      :: FSATMX              ! maximum surface saturated fraction (global mean)
    real(kind=kind_noahmp)      :: Z0SNO               ! snow surface roughness length (m) (0.002)
    real(kind=kind_noahmp)      :: SSI                 ! liquid water holding capacity for snowpack (m3/m3) (0.03)
    real(kind=kind_noahmp)      :: SNOW_RET_FAC        ! snowpack water release timescale factor (1/s)
    real(kind=kind_noahmp)      :: SNOW_EMIS           ! snow emissivity
    real(kind=kind_noahmp)      :: SWEMX               ! new snow mass to fully cover old snow (mm)
    real(kind=kind_noahmp)      :: RSURF_SNOW          ! surface resistance for snow(s/m)
    real(kind=kind_noahmp)      :: TAU0                ! tau0 from Yang97 eqn. 10a
    real(kind=kind_noahmp)      :: GRAIN_GROWTH        ! growth from vapor diffusion Yang97 eqn. 10b
    real(kind=kind_noahmp)      :: EXTRA_GROWTH        ! extra growth near freezing Yang97 eqn. 10c
    real(kind=kind_noahmp)      :: DIRT_SOOT           ! dirt and soot term Yang97 eqn. 10d
    real(kind=kind_noahmp)      :: BATS_COSZ           ! zenith angle snow albedo adjustment; b in Yang97 eqn. 15
    real(kind=kind_noahmp)      :: BATS_VIS_NEW        ! new snow visible albedo
    real(kind=kind_noahmp)      :: BATS_NIR_NEW        ! new snow NIR albedo
    real(kind=kind_noahmp)      :: BATS_VIS_AGE        ! age factor for diffuse visible snow albedo Yang97 eqn. 17
    real(kind=kind_noahmp)      :: BATS_NIR_AGE        ! age factor for diffuse NIR snow albedo Yang97 eqn. 18
    real(kind=kind_noahmp)      :: BATS_VIS_DIR        ! cosz factor for direct visible snow albedo Yang97 eqn. 15
    real(kind=kind_noahmp)      :: BATS_NIR_DIR        ! cosz factor for direct NIR snow albedo Yang97 eqn. 16
    real(kind=kind_noahmp)      :: RSURF_EXP           ! exponent in the shape parameter for soil resistance option 1
    real(kind=kind_noahmp)      :: C2_SNOWCOMPACT      ! overburden snow compaction parameter (m3/kg)
    real(kind=kind_noahmp)      :: C3_SNOWCOMPACT      ! snow desctructive metamorphism compaction parameter1 [1/s]
    real(kind=kind_noahmp)      :: C4_SNOWCOMPACT      ! snow desctructive metamorphism compaction parameter2 [1/k]
    real(kind=kind_noahmp)      :: C5_SNOWCOMPACT      ! snow desctructive metamorphism compaction parameter3
    real(kind=kind_noahmp)      :: DM_SNOWCOMPACT      ! upper Limit on destructive metamorphism compaction [kg/m3]
    real(kind=kind_noahmp)      :: ETA0_SNOWCOMPACT    ! snow viscosity coefficient [kg-s/m2]
    real(kind=kind_noahmp)      :: SNOWCOMPACTm_AR24   ! snow compaction m parameter for linear sfc temp fitting from AR24
    real(kind=kind_noahmp)      :: SNOWCOMPACTb_AR24   ! snow compaction b parameter for linear sfc temp fitting from AR24
    real(kind=kind_noahmp)      :: SNOWCOMPACT_P1_AR24 ! lower constraint for SnowCompactBurdenFac for high pressure bin from AR24
    real(kind=kind_noahmp)      :: SNOWCOMPACT_P2_AR24 ! lower constraint for SnowCompactBurdenFac for mid pressure bin from AR24
    real(kind=kind_noahmp)      :: SNOWCOMPACT_P3_AR24 ! lower constraint for SnowCompactBurdenFac for low pressure bin from AR24
    real(kind=kind_noahmp)      :: SNOWCOMPACT_Up_AR24 ! upper constraint on SnowCompactBurdenFac from AR24
    real(kind=kind_noahmp)      :: SCFm1_AR25          ! m1 parameter for ground SCF from AR2025
    real(kind=kind_noahmp)      :: SCFm2_AR25          ! m2 parameter for ground SCF from AR2025
    real(kind=kind_noahmp)      :: SCfac1_AR25         ! SCfac1 parameter for ground SCF from AR2025
    real(kind=kind_noahmp)      :: SCfac2_AR25         ! SCfac2 parameter for ground SCF from AR2025
    real(kind=kind_noahmp)      :: SNLIQMAXFRAC        ! maximum liquid water fraction in snow
    real(kind=kind_noahmp)      :: SWEMAXGLA           ! Maximum SWE allowed at glaciers (mm)
    real(kind=kind_noahmp)      :: WSLMAX              ! maximum lake water storage (mm)
    real(kind=kind_noahmp)      :: ROUS                ! specific yield [-] for Niu et al. 2007 groundwater scheme
    real(kind=kind_noahmp)      :: CMIC                ! microprore content (0.0-1.0), 0.0: close to free drainage
    real(kind=kind_noahmp)      :: SNOWDEN_MAX         ! maximum fresh snowfall density (kg/m3)
    real(kind=kind_noahmp)      :: CLASS_ALB_REF       ! reference snow albedo in CLASS scheme
    real(kind=kind_noahmp)      :: CLASS_SNO_AGE       ! snow aging e-folding time (s) in CLASS albedo scheme
    real(kind=kind_noahmp)      :: CLASS_ALB_NEW       ! fresh snow albedo in CLASS scheme
    real(kind=kind_noahmp)      :: PSIWLT              ! soil metric potential for wilting point (m)
    real(kind=kind_noahmp)      :: Z0SOIL              ! Bare-soil roughness length (m) (i.e., under the canopy)
    real(kind=kind_noahmp)      :: Z0LAKE              ! Lake surface roughness length (m)

    ! irrigation parameters
    integer                     :: IRR_HAR             ! number of days before harvest date to stop irrigation 
    real(kind=kind_noahmp)      :: IRR_FRAC            ! irrigation Fraction
    real(kind=kind_noahmp)      :: IRR_LAI             ! Minimum lai to trigger irrigation
    real(kind=kind_noahmp)      :: IRR_MAD             ! management allowable deficit (0-1)
    real(kind=kind_noahmp)      :: FILOSS              ! factor of flood irrigation loss
    real(kind=kind_noahmp)      :: SPRIR_RATE          ! mm/h, sprinkler irrigation rate
    real(kind=kind_noahmp)      :: MICIR_RATE          ! mm/h, micro irrigation rate
    real(kind=kind_noahmp)      :: FIRTFAC             ! flood application rate factor
    real(kind=kind_noahmp)      :: IR_RAIN             ! maximum precipitation to stop irrigation trigger

    ! tile drainage parameters
    integer                     :: DRAIN_LAYER_OPT     ! tile drainage layer
    integer                     :: TD_DEPTH            ! tile drainage depth (layer number) from soil surface
    real(kind=kind_noahmp)      :: TDSMC_FAC           ! tile drainage soil moisture factor
    real(kind=kind_noahmp)      :: TD_DC               ! tile drainage coefficient [mm/d]
    real(kind=kind_noahmp)      :: TD_DCOEF            ! tile drainage coefficient [mm/d]
    real(kind=kind_noahmp)      :: TD_D                ! depth to impervious layer from drain water level [m]
    real(kind=kind_noahmp)      :: TD_ADEPTH           ! actual depth of impervious layer from land surface [m]
    real(kind=kind_noahmp)      :: TD_RADI             ! effective radius of drain tubes [m]
    real(kind=kind_noahmp)      :: TD_SPAC             ! distance between two drain tubes or tiles [m]
    real(kind=kind_noahmp)      :: TD_DDRAIN           ! tile drainage depth [m]
    real(kind=kind_noahmp)      :: KLAT_FAC            ! hydraulic conductivity mutiplification factor

    ! crop parameters
    integer                     :: PLTDAY              ! Planting date
    integer                     :: HSDAY               ! Harvest date
    real(kind=kind_noahmp)      :: PLANTPOP            ! Plant density [per ha] - used?
    real(kind=kind_noahmp)      :: IRRI                ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
    real(kind=kind_noahmp)      :: GDDTBASE            ! Base temperature for GDD accumulation [C]
    real(kind=kind_noahmp)      :: GDDTCUT             ! Upper temperature for GDD accumulation [C]
    real(kind=kind_noahmp)      :: GDDS1               ! GDD from seeding to emergence
    real(kind=kind_noahmp)      :: GDDS2               ! GDD from seeding to initial vegetative 
    real(kind=kind_noahmp)      :: GDDS3               ! GDD from seeding to post vegetative 
    real(kind=kind_noahmp)      :: GDDS4               ! GDD from seeding to intial reproductive
    real(kind=kind_noahmp)      :: GDDS5               ! GDD from seeding to pysical maturity 
    real(kind=kind_noahmp)      :: AREF                ! reference maximum CO2 assimulation rate 
    real(kind=kind_noahmp)      :: PSNRF               ! CO2 assimulation reduction factor(0-1) (caused by non-modeled part, pest,weeds)
    real(kind=kind_noahmp)      :: I2PAR               ! Fraction of incoming solar radiation to photosynthetically active radiation
    real(kind=kind_noahmp)      :: TASSIM0             ! Minimum temperature for CO2 assimulation [C]
    real(kind=kind_noahmp)      :: TASSIM1             ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
    real(kind=kind_noahmp)      :: TASSIM2             ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
    real(kind=kind_noahmp)      :: K                   ! light extinction coefficient
    real(kind=kind_noahmp)      :: EPSI                ! initial light use efficiency
    real(kind=kind_noahmp)      :: Q10MR               ! q10 for maintainance respiration
    real(kind=kind_noahmp)      :: LEFREEZ             ! characteristic T for leaf freezing [K]
    real(kind=kind_noahmp)      :: DILE_FC(NSTAGE)     ! coeficient for temperature leaf stress death [1/s]
    real(kind=kind_noahmp)      :: DILE_FW(NSTAGE)     ! coeficient for water leaf stress death [1/s]
    real(kind=kind_noahmp)      :: FRA_GR              ! fraction of growth respiration
    real(kind=kind_noahmp)      :: LF_OVRC(NSTAGE)     ! fraction of leaf turnover  [1/s]
    real(kind=kind_noahmp)      :: ST_OVRC(NSTAGE)     ! fraction of stem turnover  [1/s]
    real(kind=kind_noahmp)      :: RT_OVRC(NSTAGE)     ! fraction of root tunrover  [1/s]
    real(kind=kind_noahmp)      :: LFMR25              ! leaf maintenance respiration at 25C [umol CO2/m2/s]
    real(kind=kind_noahmp)      :: STMR25              ! stem maintenance respiration at 25C [umol CO2/kg bio/s]
    real(kind=kind_noahmp)      :: RTMR25              ! root maintenance respiration at 25C [umol CO2/kg bio/s]
    real(kind=kind_noahmp)      :: GRAINMR25           ! grain maintenance respiration at 25C [umol CO2/kg bio/s]
    real(kind=kind_noahmp)      :: LFPT(NSTAGE)        ! fraction of carbohydrate flux to leaf
    real(kind=kind_noahmp)      :: STPT(NSTAGE)        ! fraction of carbohydrate flux to stem
    real(kind=kind_noahmp)      :: RTPT(NSTAGE)        ! fraction of carbohydrate flux to root
    real(kind=kind_noahmp)      :: GRAINPT(NSTAGE)     ! fraction of carbohydrate flux to grain
    real(kind=kind_noahmp)      :: LFCT(NSTAGE)        ! fraction of carbohydrate translocation from leaf to grain 
    real(kind=kind_noahmp)      :: STCT(NSTAGE)        ! fraction of carbohydrate translocation from stem to grain
    real(kind=kind_noahmp)      :: RTCT(NSTAGE)        ! fraction of carbohydrate translocation from root to grain
    real(kind=kind_noahmp)      :: BIO2LAI             ! leaf area per living leaf biomass [m2/kg]

    ! wetland parameters
    real(kind=kind_noahmp)      :: WCAP                ! maximum wetland water holding capacity [m] (tunable) for opt_wetland=1

    ! soil parameters
    real(kind=kind_noahmp)      :: BEXP(NSOIL)         ! soil B parameter
    real(kind=kind_noahmp)      :: SMCDRY(NSOIL)       ! dry soil moisture threshold
    real(kind=kind_noahmp)      :: SMCMAX(NSOIL)       ! porosity, saturated value of soil moisture (volumetric)
    real(kind=kind_noahmp)      :: SMCREF(NSOIL)       ! reference soil moisture (field capacity) (volumetric)
    real(kind=kind_noahmp)      :: PSISAT(NSOIL)       ! saturated soil matric potential
    real(kind=kind_noahmp)      :: DKSAT(NSOIL)        ! saturated soil hydraulic conductivity
    real(kind=kind_noahmp)      :: DWSAT(NSOIL)        ! saturated soil hydraulic diffusivity
    real(kind=kind_noahmp)      :: SMCWLT(NSOIL)       ! wilting point soil moisture (volumetric)
    real(kind=kind_noahmp)      :: QUARTZ(NSOIL)       ! soil quartz content
    real(kind=kind_noahmp)      :: BVIC                ! VIC model infiltration parameter (-) for opt_run=6
    real(kind=kind_noahmp)      :: AXAJ                ! Xinanjiang: Tension water distribution inflection parameter [-] for opt_run=7
    real(kind=kind_noahmp)      :: BXAJ                ! Xinanjiang: Tension water distribution shape parameter [-] for opt_run=7
    real(kind=kind_noahmp)      :: XXAJ                ! Xinanjiang: Free water distribution shape parameter [-] for opt_run=7
    real(kind=kind_noahmp)      :: BDVIC               ! VIC model infiltration parameter (-)
    real(kind=kind_noahmp)      :: GDVIC               ! mean capilary drive (m)
    real(kind=kind_noahmp)      :: BBVIC               ! heterogeniety parameter for DVIC infiltration [-]

    ! general parameters
    real(kind=kind_noahmp)      :: SLOPE               ! slope factor for soil drainage
    real(kind=kind_noahmp)      :: CSOIL               ! Soil heat capacity [J m-3 K-1]
    real(kind=kind_noahmp)      :: REFDK               ! Parameter in the surface runoff parameterization
    real(kind=kind_noahmp)      :: REFKDT              ! Parameter in the surface runoff parameterization
    real(kind=kind_noahmp)      :: KDT                 ! used in compute maximum infiltration rate (in INFIL)
    real(kind=kind_noahmp)      :: FRZX                ! used in compute maximum infiltration rate (in INFIL)
    real(kind=kind_noahmp)      :: FRZK                ! Frozen ground parameter
    real(kind=kind_noahmp)      :: ZBOT                ! Depth [m] of lower boundary soil temperature
    real(kind=kind_noahmp)      :: CZIL                ! Parameter used in the calculation of the roughness length for heat
    real(kind=kind_noahmp)      :: mxsnalb             ! LIS specific: max snow albedo
    real(kind=kind_noahmp)      :: mnsnalb             ! LIS specific: min snow albedo
    real(kind=kind_noahmp)      :: sndecayexp          ! LIS specific: snow age exponential decay
    real(kind=kind_noahmp)      :: t_ulimit            ! LIS specific: 
    real(kind=kind_noahmp)      :: t_mlimit            ! LIS specific: 
    real(kind=kind_noahmp)      :: t_llimit            ! LIS specific: 
    real(kind=kind_noahmp)      :: snowf_scalef        ! LIS specific: snow cover scaling factor

    ! SNICAR parameter
    real(kind=kind_noahmp)      :: DepBChydropho       ! hydrophobic Black Carbon deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)      :: DepBChydrophi       ! hydrophillic Black Carbon deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)      :: DepOChydropho       ! hydrophobic Organic Carbon deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)      :: DepOChydrophi       ! hydrophillic Organic Carbon deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)      :: DepDust1            ! dust species 1 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)      :: DepDust2            ! dust species 2 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)      :: DepDust3            ! dust species 3 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)      :: DepDust4            ! dust species 4 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)      :: DepDust5            ! dust species 5 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)      :: SnowRadiusMin       ! minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
    real(kind=kind_noahmp)      :: FreshSnowRadiusMax  ! maximum warm fresh snow effective radius [microns]
    real(kind=kind_noahmp)      :: SnowRadiusRefrz     ! effective radius of re-frozen snow [microns]
    real(kind=kind_noahmp)      :: ScavEffMeltScale    ! Scaling factor modifying scavenging factors for aerosol in meltwater (-)
    real(kind=kind_noahmp)      :: ScavEffMeltBCphi    ! scavenging factor for hydrophillic BC inclusion in meltwater [frc]
    real(kind=kind_noahmp)      :: ScavEffMeltBCpho    ! scavenging factor for hydrophobic BC inclusion in meltwater  [frc]
    real(kind=kind_noahmp)      :: ScavEffMeltOCphi    ! scavenging factor for hydrophillic OC inclusion in meltwater [frc]
    real(kind=kind_noahmp)      :: ScavEffMeltOCpho    ! scavenging factor for hydrophobic OC inclusion in meltwater  [frc]
    real(kind=kind_noahmp)      :: ScavEffMeltDust1    ! scavenging factor for dust species 1 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)      :: ScavEffMeltDust2    ! scavenging factor for dust species 2 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)      :: ScavEffMeltDust3    ! scavenging factor for dust species 3 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)      :: ScavEffMeltDust4    ! scavenging factor for dust species 4 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)      :: ScavEffMeltDust5    ! scavenging factor for dust species 5 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)      :: SnowRadiusMax       ! maximum allowed snow effective radius [microns]
    real(kind=kind_noahmp)      :: SnowWetAgeC1Brun89  ! constant for liquid water grain growth [m3 s-1], from Brun89
    real(kind=kind_noahmp)      :: SnowWetAgeC2Brun89  ! Constant for liquid water grain growth [m3 s-1], from Brun89: corrected for LWC 
    real(kind=kind_noahmp)      :: SnowAgeScaleFac     ! Arbitrary scaling factor applied to snow aging rate (-)

  end type LisNoahmpParam_type

end module LisNoahmpParamType
