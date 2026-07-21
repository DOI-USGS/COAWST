module EnergyVarInTransferMod

!!! Transfer input 2-D NoahmpIO Energy variables to 1-D column variable
!!! 1-D variables should be first defined in /src/EnergyVarType.F90
!!! 2-D variables should be first defined in NoahmpIOVarType.F90

! ------------------------ Code history -----------------------------------
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

  use Machine
  use NoahmpIOVarType
  use NoahmpVarType
  use LisNoahmpParamType

  implicit none

contains

!=== initialize with input data or table values

  subroutine EnergyVarInTransfer(noahmp, NoahmpIO, LISparam)

    implicit none

    type(NoahmpIO_type),       intent(inout) :: NoahmpIO
    type(noahmp_type),         intent(inout) :: noahmp
    type(LisNoahmpParam_type), intent(in)    :: LISparam   ! lis/noahmp parameter

    ! local loop index
    integer                                  :: SoilLayerIndex

! -------------------------------------------------------------------------
    associate(                                                                    &
              I                     => noahmp%config%domain%GridIndexI           ,&
              J                     => noahmp%config%domain%GridIndexJ           ,&
              VegType               => noahmp%config%domain%VegType              ,&
              SoilType              => noahmp%config%domain%SoilType             ,&
              CropType              => noahmp%config%domain%CropType             ,&
              SoilColor             => noahmp%config%domain%SoilColor            ,&
              FlagUrban             => noahmp%config%domain%FlagUrban            ,&
              NumSnowLayerMax       => noahmp%config%domain%NumSnowLayerMax      ,&
              NumSoilLayer          => noahmp%config%domain%NumSoilLayer         ,&
              NumSwRadBand          => noahmp%config%domain%NumSwRadBand         ,&
              NumSnicarRadBand      => noahmp%config%domain%NumSnicarRadBand     ,&
              NumRadiusSnwMieSnicar => noahmp%config%domain%NumRadiusSnwMieSnicar &
              )
! -------------------------------------------------------------------------

    ! energy state variables
    noahmp%energy%state%LeafAreaIndex                             = NoahmpIO%LAI     (I,J)
    noahmp%energy%state%StemAreaIndex                             = NoahmpIO%XSAIXY  (I,J)
    noahmp%energy%state%SpecHumiditySfcMean                       = NoahmpIO%QSFC    (I,J)
    noahmp%energy%state%TemperatureGrd                            = NoahmpIO%TGXY    (I,J)
    noahmp%energy%state%TemperatureCanopy                         = NoahmpIO%TVXY    (I,J)
    noahmp%energy%state%SnowAgeNondim                             = NoahmpIO%TAUSSXY (I,J)
    noahmp%energy%state%AlbedoSnowPrev                            = NoahmpIO%ALBOLDXY(I,J)
    noahmp%energy%state%PressureVaporCanAir                       = NoahmpIO%EAHXY   (I,J)
    noahmp%energy%state%TemperatureCanopyAir                      = NoahmpIO%TAHXY   (I,J)
    noahmp%energy%state%ExchCoeffShSfc                            = NoahmpIO%CHXY    (I,J) 
    noahmp%energy%state%ExchCoeffMomSfc                           = NoahmpIO%CMXY    (I,J)
    noahmp%energy%state%TemperatureSoilSnow(-NumSnowLayerMax+1:0) = NoahmpIO%TSNOXY  (I,-NumSnowLayerMax+1:0,J)
    noahmp%energy%state%TemperatureSoilSnow(1:NumSoilLayer)       = NoahmpIO%TSLB    (I,1:NumSoilLayer,J)
    noahmp%energy%state%PressureAtmosCO2                          = LISparam%CO2 * noahmp%forcing%PressureAirRefHeight
    noahmp%energy%state%PressureAtmosO2                           = LISparam%O2  * noahmp%forcing%PressureAirRefHeight
    noahmp%energy%state%AlbedoSoilDir(1:NumSwRadBand)             = NoahmpIO%ALBSOILDIRXY(I,1:NumSwRadBand,J)
    noahmp%energy%state%AlbedoSoilDif(1:NumSwRadBand)             = NoahmpIO%ALBSOILDIFXY(I,1:NumSwRadBand,J)
    ! vegetation treatment for USGS land types (playa, lava, sand to bare)
    if ( (VegType == 25) .or. (VegType == 26) .or. (VegType == 27) ) then
       noahmp%energy%state%VegFrac       = 0.0
       noahmp%energy%state%LeafAreaIndex = 0.0
    endif

    ! energy flux variables
    noahmp%energy%flux%HeatGroundTotAcc                           = NoahmpIO%ACC_SSOILXY(I,J)

    ! energy parameter variables
    noahmp%energy%param%SoilHeatCapacity                          = LISparam%CSOIL
    noahmp%energy%param%SnowAgeFacBats                            = LISparam%TAU0
    noahmp%energy%param%SnowGrowVapFacBats                        = LISparam%GRAIN_GROWTH
    noahmp%energy%param%SnowSootFacBats                           = LISparam%DIRT_SOOT
    noahmp%energy%param%SnowGrowFrzFacBats                        = LISparam%EXTRA_GROWTH
    noahmp%energy%param%SolarZenithAdjBats                        = LISparam%BATS_COSZ
    noahmp%energy%param%FreshSnoAlbVisBats                        = LISparam%BATS_VIS_NEW
    noahmp%energy%param%FreshSnoAlbNirBats                        = LISparam%BATS_NIR_NEW
    noahmp%energy%param%SnoAgeFacDifVisBats                       = LISparam%BATS_VIS_AGE
    noahmp%energy%param%SnoAgeFacDifNirBats                       = LISparam%BATS_NIR_AGE
    noahmp%energy%param%SzaFacDirVisBats                          = LISparam%BATS_VIS_DIR
    noahmp%energy%param%SzaFacDirNirBats                          = LISparam%BATS_NIR_DIR
    noahmp%energy%param%SnowAlbRefClass                           = LISparam%CLASS_ALB_REF
    noahmp%energy%param%SnowAgeFacClass                           = LISparam%CLASS_SNO_AGE
    noahmp%energy%param%SnowAlbFreshClass                         = LISparam%CLASS_ALB_NEW
    noahmp%energy%param%UpscatterCoeffSnowDir                     = LISparam%BETADS
    noahmp%energy%param%UpscatterCoeffSnowDif                     = LISparam%BETAIS
    noahmp%energy%param%ZilitinkevichCoeff                        = LISparam%CZIL
    noahmp%energy%param%EmissivitySnow                            = LISparam%SNOW_EMIS
    noahmp%energy%param%EmissivitySoilLake                        = LISparam%EG
    noahmp%energy%param%AlbedoLandIce                             = LISparam%ALBICE
    noahmp%energy%param%RoughLenMomSnow                           = LISparam%Z0SNO
    noahmp%energy%param%RoughLenMomSoil                           = LISparam%Z0SOIL
    noahmp%energy%param%RoughLenMomLake                           = LISparam%Z0LAKE
    noahmp%energy%param%EmissivityIceSfc                          = LISparam%EICE
    noahmp%energy%param%ResistanceSoilExp                         = LISparam%RSURF_EXP
    noahmp%energy%param%ResistanceSnowSfc                         = LISparam%RSURF_SNOW
    noahmp%energy%param%VegFracAnnMax                             = NoahmpIO%GVFMAX(I,J) / 100.0
    noahmp%energy%param%VegFracGreen                              = NoahmpIO%VEGFRA(I,J) / 100.0
    noahmp%energy%param%TreeCrownRadius                           = LISparam%RC
    noahmp%energy%param%HeightCanopyTop                           = LISparam%HVT
    noahmp%energy%param%HeightCanopyBot                           = LISparam%HVB
    noahmp%energy%param%RoughLenMomVeg                            = LISparam%Z0MVT
    noahmp%energy%param%CanopyWindExtFac                          = LISparam%CWPVT
    noahmp%energy%param%TreeDensity                               = LISparam%DEN
    noahmp%energy%param%CanopyOrientIndex                         = LISparam%XL
    noahmp%energy%param%ConductanceLeafMin                        = LISparam%BP
    noahmp%energy%param%Co2MmConst25C                             = LISparam%KC25
    noahmp%energy%param%O2MmConst25C                              = LISparam%KO25
    noahmp%energy%param%Co2MmConstQ10                             = LISparam%AKC
    noahmp%energy%param%O2MmConstQ10                              = LISparam%AKO
    noahmp%energy%param%RadiationStressFac                        = LISparam%RGL
    noahmp%energy%param%ResistanceStomataMin                      = LISparam%RSMIN
    noahmp%energy%param%ResistanceStomataMax                      = LISparam%RSMAX
    noahmp%energy%param%AirTempOptimTransp                        = LISparam%TOPT
    noahmp%energy%param%VaporPresDeficitFac                       = LISparam%HS
    noahmp%energy%param%LeafDimLength                             = LISparam%DLEAF
    noahmp%energy%param%HeatCapacCanFac                           = LISparam%CBIOM
    noahmp%energy%param%LeafAreaIndexMon                          = LISparam%LAIM
    noahmp%energy%param%StemAreaIndexMon                          = LISparam%SAIM
    noahmp%energy%param%ReflectanceLeaf                           = LISparam%RHOL
    noahmp%energy%param%ReflectanceStem                           = LISparam%RHOS
    noahmp%energy%param%TransmittanceLeaf                         = LISparam%TAUL
    noahmp%energy%param%TransmittanceStem                         = LISparam%TAUS
    noahmp%energy%param%AlbedoSoilSat                             = LISparam%ALBSAT
    noahmp%energy%param%AlbedoSoilDry                             = LISparam%ALBDRY
    noahmp%energy%param%AlbedoLakeFrz                             = LISparam%ALBLAK
    noahmp%energy%param%ScatterCoeffSnow                          = LISparam%OMEGAS

    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then ! SNICAR variables
       noahmp%energy%param%RadSwWgtDif        (1:NumSnicarRadBand) = NoahmpIO%flx_wgt_dif(1:NumSnicarRadBand)
       noahmp%energy%param%RadSwWgtDir        (1:NumSnicarRadBand) = NoahmpIO%flx_wgt_dir(1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbBCphi         (1:NumSnicarRadBand) = NoahmpIO%ss_alb_bc1       (1:NumSnicarRadBand) 
       noahmp%energy%param%AsyPrmBCphi        (1:NumSnicarRadBand) = NoahmpIO%asm_prm_bc1      (1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassBCphi    (1:NumSnicarRadBand) = NoahmpIO%ext_cff_mss_bc1  (1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbBCpho         (1:NumSnicarRadBand) = NoahmpIO%ss_alb_bc2       (1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmBCpho        (1:NumSnicarRadBand) = NoahmpIO%asm_prm_bc2      (1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassBCpho    (1:NumSnicarRadBand) = NoahmpIO%ext_cff_mss_bc2  (1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbOCphi         (1:NumSnicarRadBand) = NoahmpIO%ss_alb_oc1       (1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmOCphi        (1:NumSnicarRadBand) = NoahmpIO%asm_prm_oc1      (1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassOCphi    (1:NumSnicarRadBand) = NoahmpIO%ext_cff_mss_oc1  (1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbOCpho         (1:NumSnicarRadBand) = NoahmpIO%ss_alb_oc2       (1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmOCpho        (1:NumSnicarRadBand) = NoahmpIO%asm_prm_oc2      (1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassOCpho    (1:NumSnicarRadBand) = NoahmpIO%ext_cff_mss_oc2  (1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbDustB1        (1:NumSnicarRadBand) = NoahmpIO%ss_alb_dst1      (1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmDustB1       (1:NumSnicarRadBand) = NoahmpIO%asm_prm_dst1     (1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassDustB1   (1:NumSnicarRadBand) = NoahmpIO%ext_cff_mss_dst1 (1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbDustB2        (1:NumSnicarRadBand) = NoahmpIO%ss_alb_dst2      (1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmDustB2       (1:NumSnicarRadBand) = NoahmpIO%asm_prm_dst2     (1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassDustB2   (1:NumSnicarRadBand) = NoahmpIO%ext_cff_mss_dst2 (1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbDustB3        (1:NumSnicarRadBand) = NoahmpIO%ss_alb_dst3      (1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmDustB3       (1:NumSnicarRadBand) = NoahmpIO%asm_prm_dst3     (1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassDustB3   (1:NumSnicarRadBand) = NoahmpIO%ext_cff_mss_dst3 (1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbDustB4        (1:NumSnicarRadBand) = NoahmpIO%ss_alb_dst4      (1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmDustB4       (1:NumSnicarRadBand) = NoahmpIO%asm_prm_dst4     (1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassDustB4   (1:NumSnicarRadBand) = NoahmpIO%ext_cff_mss_dst4 (1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbDustB5        (1:NumSnicarRadBand) = NoahmpIO%ss_alb_dst5      (1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmDustB5       (1:NumSnicarRadBand) = NoahmpIO%asm_prm_dst5     (1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassDustB5   (1:NumSnicarRadBand) = NoahmpIO%ext_cff_mss_dst5 (1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbSnwRadDir     (1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand) = &
                  NoahmpIO%ss_alb_snw_drc     (1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmSnwRadDir    (1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand) = &
                  NoahmpIO%asm_prm_snw_drc    (1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassSnwRadDir(1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand) = &
                  NoahmpIO%ext_cff_mss_snw_drc(1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand)
       noahmp%energy%param%SsAlbSnwRadDif     (1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand) = &
                  NoahmpIO%ss_alb_snw_dfs     (1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand)
       noahmp%energy%param%AsyPrmSnwRadDif    (1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand) = &
                  NoahmpIO%asm_prm_snw_dfs    (1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand)
       noahmp%energy%param%ExtCffMassSnwRadDif(1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand) = &
                  NoahmpIO%ext_cff_mss_snw_dfs(1:NumRadiusSnwMieSnicar,1:NumSnicarRadBand)
    endif

    do SoilLayerIndex = 1, size(SoilType)
       noahmp%energy%param%SoilQuartzFrac(SoilLayerIndex) = LISparam%QUARTZ(SoilLayerIndex)
    enddo

    ! spatial varying soil input
    if ( noahmp%config%nmlist%OptSoilProperty == 4 ) then
       noahmp%energy%param%SoilQuartzFrac(1:NumSoilLayer) = NoahmpIO%QUARTZ_3D(I,1:NumSoilLayer,J)
    endif

    if ( FlagUrban .eqv. .true. ) noahmp%energy%param%SoilHeatCapacity = 3.0e6

    if ( CropType > 0 ) then
       noahmp%energy%param%ConductanceLeafMin             = LISparam%BP
       noahmp%energy%param%Co2MmConst25C                  = LISparam%KC25
       noahmp%energy%param%O2MmConst25C                   = LISparam%KO25
       noahmp%energy%param%Co2MmConstQ10                  = LISparam%AKC
       noahmp%energy%param%O2MmConstQ10                   = LISparam%AKO
    endif

    end associate

  end subroutine EnergyVarInTransfer

end module EnergyVarInTransferMod
