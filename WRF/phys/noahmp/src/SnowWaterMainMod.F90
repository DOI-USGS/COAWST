module SnowWaterMainMod

!!! Main snow water module including all snowpack processes
!!! Snowfall -> Snowpack compaction -> Snow layer combination -> Snow layer division -> Snow Hydrology

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowfallBelowCanopyMod,    only : SnowfallAfterCanopyIntercept
  use SnowpackCompactionMod,     only : SnowpackCompaction
  use SnowpackCompactionAR24Mod, only : SnowpackCompactionAR24 
  use SnowLayerCombineMod,       only : SnowLayerCombine
  use SnowLayerDivideMod,        only : SnowLayerDivide
  use SnowpackHydrologyMod,      only : SnowpackHydrology
  use SnowAerosolSnicarMod,      only : SnowAerosolSnicar

  implicit none

contains

  subroutine SnowWaterMain(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: SNOWWATER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd                 ! do loop/array indices
    real(kind=kind_noahmp)           :: SnowDensBulk            ! bulk density of snow [kg/m3]
    real(kind=kind_noahmp)           :: GlacierExcessRemainFrac ! fraction of mass remaining after glacier excess flow

! --------------------------------------------------------------------
    associate(                                                                       &
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax        ,& ! in,    maximum number of snow layers
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              MainTimeStep           => noahmp%config%domain%MainTimeStep           ,& ! in,    noahmp main time step [s]
              DepthSoilLayer         => noahmp%config%domain%DepthSoilLayer         ,& ! in,    depth [m] of layer-bottom from soil surface
              OptSnowAlbedo          => noahmp%config%nmlist%OptSnowAlbedo          ,& ! in,    options for ground snow surface albedo
              OptSnowCompaction      => noahmp%config%nmlist%OptSnowCompaction      ,& ! in,    options for ground snowpack compaction
              SnoWatEqvMaxGlacier    => noahmp%water%param%SnoWatEqvMaxGlacier      ,& ! in,    Maximum SWE allowed at glaciers [mm]
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
              DepthSnowSoilLayer     => noahmp%config%domain%DepthSnowSoilLayer     ,& ! inout, depth of snow/soil layer-bottom [m]
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg        ,& ! inout, actual number of snow layers (negative)
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow     ,& ! inout, snow and soil layer temperature [K]
              SnowDepth              => noahmp%water%state%SnowDepth                ,& ! inout, snow depth [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv           ,& ! inout, snow water equivalent [mm]
              SnowIce                => noahmp%water%state%SnowIce                  ,& ! inout, snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater             ,& ! inout, snow layer liquid water [mm]
              MassBChydropho         => noahmp%water%state%MassBChydropho           ,& ! inout, mass of hydrophobic Black Carbon in snow [kg m-2]
              MassBChydrophi         => noahmp%water%state%MassBChydrophi           ,& ! inout, mass of hydrophillic Black Carbon in snow [kg m-2]
              MassOChydropho         => noahmp%water%state%MassOChydropho           ,& ! inout, mass of hydrophobic Organic Carbon in snow [kg m-2]
              MassOChydrophi         => noahmp%water%state%MassOChydrophi           ,& ! inout, mass of hydrophillic Organic Carbon in snow [kg m-2]
              MassDust1              => noahmp%water%state%MassDust1                ,& ! inout, mass of dust species 1 in snow [kg m-2]
              MassDust2              => noahmp%water%state%MassDust2                ,& ! inout, mass of dust species 2 in snow [kg m-2]
              MassDust3              => noahmp%water%state%MassDust3                ,& ! inout, mass of dust species 3 in snow [kg m-2]
              MassDust4              => noahmp%water%state%MassDust4                ,& ! inout, mass of dust species 4 in snow [kg m-2]
              MassDust5              => noahmp%water%state%MassDust5                ,& ! inout, mass of dust species 5 in snow [kg m-2]
              GlacierExcessFlow      => noahmp%water%flux%GlacierExcessFlow         ,& ! out,   glacier snow excess flow [mm/s]
              PondSfcThinSnwComb     => noahmp%water%state%PondSfcThinSnwComb       ,& ! out,   surface ponding [mm] from liquid in thin snow layer combination
              PondSfcThinSnwTrans    => noahmp%water%state%PondSfcThinSnwTrans       & ! out,   surface ponding [mm] from thin snow when changing from multilayer to no layer
             )
! ----------------------------------------------------------------------

    ! initialize out-only variables
    GlacierExcessFlow       = 0.0
    PondSfcThinSnwComb      = 0.0
    PondSfcThinSnwTrans     = 0.0
    GlacierExcessRemainFrac = 1.0

    ! snowfall after canopy interception
    call SnowfallAfterCanopyIntercept(noahmp)

    ! do following snow layer compaction, combination, and division only for multi-layer snowpack

    ! snowpack compaction (option: 1->original,Anderson1976; 2->new,Abolafia-Rosenzweig2024)
    if ( NumSnowLayerNeg < 0 .and. OptSnowCompaction == 1) call SnowpackCompaction(noahmp)
    if ( NumSnowLayerNeg < 0 .and. OptSnowCompaction == 2) call SnowpackCompactionAR24(noahmp)

    ! snow layer combination
    if ( NumSnowLayerNeg < 0 ) call SnowLayerCombine(noahmp)

    ! snow layer division
    if ( NumSnowLayerNeg < 0 ) call SnowLayerDivide(noahmp)

    ! snow hydrology for all snow cases
    call SnowpackHydrology(noahmp)

    ! set empty snow layer properties to zero
    do LoopInd = -NumSnowLayerMax+1, NumSnowLayerNeg
       SnowIce(LoopInd)                = 0.0
       SnowLiqWater(LoopInd)           = 0.0
       TemperatureSoilSnow(LoopInd)    = 0.0
       ThicknessSnowSoilLayer(LoopInd) = 0.0
       DepthSnowSoilLayer(LoopInd)     = 0.0

       if ( (OptSnowAlbedo == 3) .and. (NumSnowLayerNeg < 0) ) then
          MassBChydropho(LoopInd)      = 0.0
          MassBChydrophi(LoopInd)      = 0.0
          MassOChydropho(LoopInd)      = 0.0
          MassOChydrophi(LoopInd)      = 0.0
          MassDust1(LoopInd)           = 0.0
          MassDust2(LoopInd)           = 0.0
          MassDust3(LoopInd)           = 0.0
          MassDust4(LoopInd)           = 0.0
          MassDust5(LoopInd)           = 0.0
       endif
    enddo

    ! to obtain equilibrium state of snow in glacier region
    if ( SnowWaterEquiv > SnoWatEqvMaxGlacier ) then
       SnowDensBulk              = SnowIce(0) / ThicknessSnowSoilLayer(0)
       GlacierExcessFlow         = SnowWaterEquiv - SnoWatEqvMaxGlacier
       SnowIce(0)                = SnowIce(0)  - GlacierExcessFlow
       ThicknessSnowSoilLayer(0) = ThicknessSnowSoilLayer(0) - GlacierExcessFlow / SnowDensBulk
       GlacierExcessFlow         = GlacierExcessFlow / MainTimeStep

       if ( OptSnowAlbedo == 3 ) then
          GlacierExcessRemainFrac = SnowIce(0) / (SnowIce(0) + GlacierExcessFlow)
          MassBChydropho(0)       = MassBChydropho(0) * GlacierExcessRemainFrac
          MassBChydrophi(0)       = MassBChydrophi(0) * GlacierExcessRemainFrac
          MassOChydropho(0)       = MassOChydropho(0) * GlacierExcessRemainFrac
          MassOChydrophi(0)       = MassOChydrophi(0) * GlacierExcessRemainFrac
          MassDust1(0)            = MassDust1(0) * GlacierExcessRemainFrac
          MassDust2(0)            = MassDust2(0) * GlacierExcessRemainFrac
          MassDust3(0)            = MassDust3(0) * GlacierExcessRemainFrac
          MassDust4(0)            = MassDust4(0) * GlacierExcessRemainFrac
          MassDust5(0)            = MassDust5(0) * GlacierExcessRemainFrac
       endif
    endif

    ! SNICAR
    if ( OptSnowAlbedo == 3 ) call SnowAerosolSnicar(noahmp)

    ! sum up snow mass for layered snow
    if ( NumSnowLayerNeg < 0 ) then  ! MB: only do for multi-layer
       SnowWaterEquiv = 0.0
       do LoopInd = NumSnowLayerNeg+1, 0
          SnowWaterEquiv = SnowWaterEquiv + SnowIce(LoopInd) + SnowLiqWater(LoopInd)
       enddo
    endif

    ! Reset DepthSnowSoilLayer and ThicknessSnowSoilLayer
    do LoopInd = NumSnowLayerNeg+1, 0
       ThicknessSnowSoilLayer(LoopInd) = -ThicknessSnowSoilLayer(LoopInd)
    enddo

    ThicknessSnowSoilLayer(1) = DepthSoilLayer(1)
    do LoopInd = 2, NumSoilLayer
       ThicknessSnowSoilLayer(LoopInd) = DepthSoilLayer(LoopInd) - DepthSoilLayer(LoopInd-1)
    enddo

    DepthSnowSoilLayer(NumSnowLayerNeg+1) = ThicknessSnowSoilLayer(NumSnowLayerNeg+1)
    do LoopInd = NumSnowLayerNeg+2, NumSoilLayer
       DepthSnowSoilLayer(LoopInd) = DepthSnowSoilLayer(LoopInd-1) + ThicknessSnowSoilLayer(LoopInd)
    enddo

    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       ThicknessSnowSoilLayer(LoopInd) = -ThicknessSnowSoilLayer(LoopInd)
    enddo

    ! Update SnowDepth for multi-layer snow
    if ( NumSnowLayerNeg < 0 ) then
       SnowDepth = 0.0
       do LoopInd = NumSnowLayerNeg+1, 0
          SnowDepth = SnowDepth + ThicknessSnowSoilLayer(LoopInd)
       enddo
    endif

    ! update snow quantity
    if ( (SnowDepth <= 1.0e-6) .or. (SnowWaterEquiv <= 1.0e-6) ) then
       SnowDepth      = 0.0
       SnowWaterEquiv = 0.0
    endif

    end associate

  end subroutine SnowWaterMain

end module SnowWaterMainMod
