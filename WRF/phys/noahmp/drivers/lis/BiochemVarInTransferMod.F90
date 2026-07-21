module BiochemVarInTransferMod

!!! Transfer input 2-D NoahmpIO Biochemistry variables to 1-D column variable
!!! 1-D variables should be first defined in /src/BiochemVarType.F90
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

  subroutine BiochemVarInTransfer(noahmp, NoahmpIO, LISparam)

    implicit none

    type(noahmp_type),         intent(inout) :: noahmp
    type(NoahmpIO_type),       intent(inout) :: NoahmpIO
    type(LisNoahmpParam_type), intent(in)    :: LISparam

! -------------------------------------------------------------------------
    associate(                                                   &
              I            => noahmp%config%domain%GridIndexI   ,&
              J            => noahmp%config%domain%GridIndexJ   ,&
              VegType      => noahmp%config%domain%VegType      ,&
              CropType     => noahmp%config%domain%CropType     ,&
              OptCropModel => noahmp%config%nmlist%OptCropModel  &
             )
! -------------------------------------------------------------------------

    ! biochem state variables
    noahmp%biochem%state%PlantGrowStage             = NoahmpIO%PGSXY   (I,J)   
    noahmp%biochem%state%LeafMass                   = NoahmpIO%LFMASSXY(I,J)
    noahmp%biochem%state%RootMass                   = NoahmpIO%RTMASSXY(I,J)
    noahmp%biochem%state%StemMass                   = NoahmpIO%STMASSXY(I,J) 
    noahmp%biochem%state%WoodMass                   = NoahmpIO%WOODXY  (I,J) 
    noahmp%biochem%state%CarbonMassDeepSoil         = NoahmpIO%STBLCPXY(I,J) 
    noahmp%biochem%state%CarbonMassShallowSoil      = NoahmpIO%FASTCPXY(I,J)
    noahmp%biochem%state%GrainMass                  = NoahmpIO%GRAINXY (I,J)  
    noahmp%biochem%state%GrowDegreeDay              = NoahmpIO%GDDXY   (I,J)  
    noahmp%biochem%state%NitrogenConcFoliage        = 1.0  ! for now, set to nitrogen saturation

    ! biochem parameter variables
    noahmp%biochem%param%NitrogenConcFoliageMax     = LISparam%FOLNMX
    noahmp%biochem%param%QuantumEfficiency25C       = LISparam%QE25
    noahmp%biochem%param%CarboxylRateMax25C         = LISparam%VCMX25
    noahmp%biochem%param%CarboxylRateMaxQ10         = LISparam%AVCMX
    noahmp%biochem%param%PhotosynPathC3             = LISparam%C3PSN
    noahmp%biochem%param%SlopeConductToPhotosyn     = LISparam%MP
    noahmp%biochem%param%RespMaintQ10               = LISparam%ARM
    noahmp%biochem%param%RespMaintLeaf25C           = LISparam%RMF25
    noahmp%biochem%param%RespMaintStem25C           = LISparam%RMS25
    noahmp%biochem%param%RespMaintRoot25C           = LISparam%RMR25
    noahmp%biochem%param%WoodToRootRatio            = LISparam%WRRAT
    noahmp%biochem%param%WoodPoolIndex              = LISparam%WDPOOL
    noahmp%biochem%param%TurnoverCoeffLeafVeg       = LISparam%LTOVRC
    noahmp%biochem%param%TemperaureLeafFreeze       = LISparam%TDLEF
    noahmp%biochem%param%LeafDeathWaterCoeffVeg     = LISparam%DILEFW
    noahmp%biochem%param%LeafDeathTempCoeffVeg      = LISparam%DILEFC
    noahmp%biochem%param%GrowthRespFrac             = LISparam%FRAGR
    noahmp%biochem%param%MicroRespCoeff             = LISparam%MRP
    noahmp%biochem%param%TemperatureMinPhotosyn     = LISparam%TMIN
    noahmp%biochem%param%LeafAreaPerMass1side       = LISparam%SLA
    noahmp%biochem%param%StemAreaIndexMin           = LISparam%XSAMIN
    noahmp%biochem%param%WoodAllocFac               = LISparam%BF
    noahmp%biochem%param%WaterStressCoeff           = LISparam%WSTRC
    noahmp%biochem%param%LeafAreaIndexMin           = LISparam%LAIMIN
    noahmp%biochem%param%TurnoverCoeffRootVeg       = LISparam%RTOVRC
    noahmp%biochem%param%WoodRespCoeff              = LISparam%RSWOODC
    ! crop model specific parameters
    if ( (OptCropModel > 0) .and. (CropType > 0) ) then
       noahmp%biochem%param%DatePlanting            = LISparam%PLTDAY
       noahmp%biochem%param%DateHarvest             = LISparam%HSDAY
       noahmp%biochem%param%NitrogenConcFoliageMax  = LISparam%FOLNMX
       noahmp%biochem%param%QuantumEfficiency25C    = LISparam%QE25
       noahmp%biochem%param%CarboxylRateMax25C      = LISparam%VCMX25
       noahmp%biochem%param%CarboxylRateMaxQ10      = LISparam%AVCMX
       noahmp%biochem%param%PhotosynPathC3          = LISparam%C3PSN
       noahmp%biochem%param%SlopeConductToPhotosyn  = LISparam%MP
       noahmp%biochem%param%RespMaintQ10            = LISparam%Q10MR
       noahmp%biochem%param%RespMaintLeaf25C        = LISparam%LFMR25
       noahmp%biochem%param%RespMaintStem25C        = LISparam%STMR25
       noahmp%biochem%param%RespMaintRoot25C        = LISparam%RTMR25
       noahmp%biochem%param%GrowthRespFrac          = LISparam%FRA_GR
       noahmp%biochem%param%TemperaureLeafFreeze    = LISparam%LEFREEZ
       noahmp%biochem%param%LeafAreaPerBiomass      = LISparam%BIO2LAI
       noahmp%biochem%param%TempBaseGrowDegDay      = LISparam%GDDTBASE
       noahmp%biochem%param%TempMaxGrowDegDay       = LISparam%GDDTCUT
       noahmp%biochem%param%GrowDegDayEmerg         = LISparam%GDDS1
       noahmp%biochem%param%GrowDegDayInitVeg       = LISparam%GDDS2
       noahmp%biochem%param%GrowDegDayPostVeg       = LISparam%GDDS3
       noahmp%biochem%param%GrowDegDayInitReprod    = LISparam%GDDS4
       noahmp%biochem%param%GrowDegDayMature        = LISparam%GDDS5
       noahmp%biochem%param%PhotosynRadFrac         = LISparam%I2PAR
       noahmp%biochem%param%TempMinCarbonAssim      = LISparam%TASSIM0
       noahmp%biochem%param%TempMaxCarbonAssim      = LISparam%TASSIM1
       noahmp%biochem%param%TempMaxCarbonAssimMax   = LISparam%TASSIM2
       noahmp%biochem%param%CarbonAssimRefMax       = LISparam%AREF
       noahmp%biochem%param%LightExtCoeff           = LISparam%K
       noahmp%biochem%param%LightUseEfficiency      = LISparam%EPSI
       noahmp%biochem%param%CarbonAssimReducFac     = LISparam%PSNRF
       noahmp%biochem%param%RespMaintGrain25C       = LISparam%GRAINMR25
       noahmp%biochem%param%LeafDeathTempCoeffCrop  = LISparam%DILE_FC
       noahmp%biochem%param%LeafDeathWaterCoeffCrop = LISparam%DILE_FW
       noahmp%biochem%param%CarbohydrLeafToGrain    = LISparam%LFCT
       noahmp%biochem%param%CarbohydrStemToGrain    = LISparam%STCT
       noahmp%biochem%param%CarbohydrRootToGrain    = LISparam%RTCT
       noahmp%biochem%param%CarbohydrFracToLeaf     = LISparam%LFPT
       noahmp%biochem%param%CarbohydrFracToStem     = LISparam%STPT
       noahmp%biochem%param%CarbohydrFracToRoot     = LISparam%RTPT
       noahmp%biochem%param%CarbohydrFracToGrain    = LISparam%GRAINPT
       noahmp%biochem%param%TurnoverCoeffLeafCrop   = LISparam%LF_OVRC
       noahmp%biochem%param%TurnoverCoeffStemCrop   = LISparam%ST_OVRC
       noahmp%biochem%param%TurnoverCoeffRootCrop   = LISparam%RT_OVRC

       if ( OptCropModel == 1 ) then
          if ( (NoahmpIO%PLANTING(I,J)>0) .and. (NoahmpIO%PLANTING(I,J)<367) ) then
             noahmp%biochem%param%DatePlanting      = NoahmpIO%PLANTING(I,J)
          endif ! 2D input map exist
          if ( (NoahmpIO%HARVEST(I,J)>0) .and. (NoahmpIO%HARVEST(I,J)<367) ) then
             noahmp%biochem%param%DateHarvest       = NoahmpIO%HARVEST(I,J)
          endif ! 2D input map exist
          if ( (NoahmpIO%SEASON_GDD(I,J)>0.0) .and. (NoahmpIO%SEASON_GDD(I,J)<10000.0) ) then
             noahmp%biochem%param%GrowDegDayEmerg   = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayEmerg
             noahmp%biochem%param%GrowDegDayInitVeg = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayInitVeg
             noahmp%biochem%param%GrowDegDayPostVeg = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayPostVeg
             noahmp%biochem%param%GrowDegDayInitReprod = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                         noahmp%biochem%param%GrowDegDayInitReprod
             noahmp%biochem%param%GrowDegDayMature  = NoahmpIO%SEASON_GDD(I,J) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayMature
          endif ! 2D input map exist
       endif ! OptCropModel == 1
    endif ! activate crop parameters

    if ( noahmp%config%nmlist%OptIrrigation == 2 ) then
       if ( (NoahmpIO%PLANTING(I,J)>0) .and. (NoahmpIO%PLANTING(I,J)<367) ) then
          noahmp%biochem%param%DatePlanting = NoahmpIO%PLANTING(I,J)
       endif ! 2D input map exist
       if ( (NoahmpIO%HARVEST(I,J)>0) .and. (NoahmpIO%HARVEST(I,J)<367) ) then
          noahmp%biochem%param%DateHarvest  = NoahmpIO%HARVEST (I,J)
       endif ! 2D input map exist
    endif
    
    end associate

  end subroutine BiochemVarInTransfer

end module BiochemVarInTransferMod
