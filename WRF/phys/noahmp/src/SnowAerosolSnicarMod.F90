module SnowAerosolSnicarMod

!!! compute aerosol content in snow and its evolution prepared for SNICAR snow albedo calculation

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowAerosolSnicar(noahmp)

! --------------------------------- Code history -----------------------------------
! Implementation: T.-S. Lin, C. He, et al. (2025, JHM)
! Adapted from CTSM, AerosolFluxes, AerosolMasses, CalcAndApplyAerosolFluxes modules
! ----------------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd                       ! do loop/array indices
    real(kind=kind_noahmp)           :: SnowMass                      ! liquid+ice snow mass in a layer [kg/m2]
    real(kind=kind_noahmp)           :: FluxInBChydrophi              ! flux of hydrophilic BC into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutBChydrophi             ! flux of hydrophilic BC out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInBChydropho              ! flux of hydrophobic BC into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutBChydropho             ! flux of hydrophobic BC out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInOChydrophi              ! flux of hydrophilic OC into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutOChydrophi             ! flux of hydrophilic OC out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInOChydropho              ! flux of hydrophobic OC into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutOChydropho             ! flux of hydrophobic OC out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust1                   ! flux of dust species 1 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust1                  ! flux of dust species 1 out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust2                   ! flux of dust species 2 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust2                  ! flux of dust species 2 out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust3                   ! flux of dust species 3 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust3                  ! flux of dust species 3 out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust4                   ! flux of dust species 4 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust4                  ! flux of dust species 4 out of layer [kg/s]
    real(kind=kind_noahmp)           :: FluxInDust5                   ! flux of dust species 5 into   layer [kg/s]
    real(kind=kind_noahmp)           :: FluxOutDust5                  ! flux of dust species 5 out of layer [kg/s]

! --------------------------------------------------------------------
    associate(                                                                          &
              MainTimeStep           => noahmp%config%domain%MainTimeStep              ,& ! in,    noahmp main time step [s]
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax           ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg           ,& ! in,    actual number of snow layers (negative)
              DepBChydropho          => noahmp%forcing%DepBChydropho                   ,& ! in,    hydrophobic Black Carbon deposition [kg m-2 s-1] 
              DepBChydrophi          => noahmp%forcing%DepBChydrophi                   ,& ! in,    hydrophillic Black Carbon deposition [kg m-2 s-1]
              DepOChydropho          => noahmp%forcing%DepOChydropho                   ,& ! in,    hydrophobic Organic Carbon deposition [kg m-2 s-1]
              DepOChydrophi          => noahmp%forcing%DepOChydrophi                   ,& ! in,    hydrophillic Organic Carbon deposition [kg m-2 s-1]
              DepDust1               => noahmp%forcing%DepDust1                        ,& ! in,    dust species 1 deposition [kg m-2 s-1]
              DepDust2               => noahmp%forcing%DepDust2                        ,& ! in,    dust species 2 deposition [kg m-2 s-1]
              DepDust3               => noahmp%forcing%DepDust3                        ,& ! in,    dust species 3 deposition [kg m-2 s-1]
              DepDust4               => noahmp%forcing%DepDust4                        ,& ! in,    dust species 4 deposition [kg m-2 s-1]
              DepDust5               => noahmp%forcing%DepDust5                        ,& ! in,    dust species 5 deposition [kg m-2 s-1]
              ScavEffMeltScale       => noahmp%water%param%ScavEffMeltScale            ,& ! in,    Scaling factor modifying scavenging factors for aerosol in meltwater (-)
              ScavEffMeltBCphi       => noahmp%water%param%ScavEffMeltBCphi            ,& ! in,    scavenging factor for hydrophillic BC inclusion in meltwater [frc]
              ScavEffMeltBCpho       => noahmp%water%param%ScavEffMeltBCpho            ,& ! in,    scavenging factor for hydrophobic BC inclusion in meltwater  [frc]
              ScavEffMeltOCphi       => noahmp%water%param%ScavEffMeltOCphi            ,& ! in,    scavenging factor for hydrophillic OC inclusion in meltwater [frc]
              ScavEffMeltOCpho       => noahmp%water%param%ScavEffMeltOCpho            ,& ! in,    scavenging factor for hydrophobic OC inclusion in meltwater  [frc]
              ScavEffMeltDust1       => noahmp%water%param%ScavEffMeltDust1            ,& ! in,    scavenging factor for dust species 1 inclusion in meltwater  [frc]
              ScavEffMeltDust2       => noahmp%water%param%ScavEffMeltDust2            ,& ! in,    scavenging factor for dust species 2 inclusion in meltwater  [frc]
              ScavEffMeltDust3       => noahmp%water%param%ScavEffMeltDust3            ,& ! in,    scavenging factor for dust species 3 inclusion in meltwater  [frc]
              ScavEffMeltDust4       => noahmp%water%param%ScavEffMeltDust4            ,& ! in,    scavenging factor for dust species 4 inclusion in meltwater  [frc]
              ScavEffMeltDust5       => noahmp%water%param%ScavEffMeltDust5            ,& ! in,    scavenging factor for dust species 5 inclusion in meltwater  [frc]
              SnowIce                => noahmp%water%state%SnowIce                     ,& ! in,    snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater                ,& ! in,    snow layer liquid water [mm]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv              ,& ! in,    snow water equivalent [mm]
              OutflowSnowLayer       => noahmp%water%flux%OutflowSnowLayer             ,& ! in,    water flow out of each snow layer [mm/s]
              MassBChydropho         => noahmp%water%state%MassBChydropho              ,& ! inout, mass of hydrophobic Black Carbon in snow [kg m-2]
              MassBChydrophi         => noahmp%water%state%MassBChydrophi              ,& ! inout, mass of hydrophillic Black Carbon in snow [kg m-2]
              MassOChydropho         => noahmp%water%state%MassOChydropho              ,& ! inout, mass of hydrophobic Organic Carbon in snow [kg m-2]
              MassOChydrophi         => noahmp%water%state%MassOChydrophi              ,& ! inout, mass of hydrophillic Organic Carbon in snow [kg m-2]
              MassDust1              => noahmp%water%state%MassDust1                   ,& ! inout, mass of dust species 1 in snow [kg m-2]
              MassDust2              => noahmp%water%state%MassDust2                   ,& ! inout, mass of dust species 2 in snow [kg m-2]
              MassDust3              => noahmp%water%state%MassDust3                   ,& ! inout, mass of dust species 3 in snow [kg m-2]
              MassDust4              => noahmp%water%state%MassDust4                   ,& ! inout, mass of dust species 4 in snow [kg m-2]
              MassDust5              => noahmp%water%state%MassDust5                   ,& ! inout, mass of dust species 5 in snow [kg m-2]
              MassConcBChydropho     => noahmp%water%state%MassConcBChydropho          ,& ! inout, mass concentration of hydrophobic Black Carbon in snow [kg/kg]
              MassConcBChydrophi     => noahmp%water%state%MassConcBChydrophi          ,& ! inout, mass concentration of hydrophillic Black Carbon in snow [kg/kg]
              MassConcOChydropho     => noahmp%water%state%MassConcOChydropho          ,& ! inout, mass concentration of hydrophobic Organic Carbon in snow [kg/kg]
              MassConcOChydrophi     => noahmp%water%state%MassConcOChydrophi          ,& ! inout, mass concentration of hydrophillic Organic Carbon in snow [kg/kg]
              MassConcDust1          => noahmp%water%state%MassConcDust1               ,& ! inout, mass concentration of dust species 1 in snow [kg/kg]
              MassConcDust2          => noahmp%water%state%MassConcDust2               ,& ! inout, mass concentration of dust species 2 in snow [kg/kg]
              MassConcDust3          => noahmp%water%state%MassConcDust3               ,& ! inout, mass concentration of dust species 3 in snow [kg/kg]
              MassConcDust4          => noahmp%water%state%MassConcDust4               ,& ! inout, mass concentration of dust species 4 in snow [kg/kg]
              MassConcDust5          => noahmp%water%state%MassConcDust5                & ! inout, mass concentration of dust species 5 in snow [kg/kg]
             )
! ----------------------------------------------------------------------

    ! initialize
    FluxInBChydropho  = 0.0
    FluxInBChydrophi  = 0.0
    FluxInOChydropho  = 0.0
    FluxInOChydrophi  = 0.0
    FluxInDust1       = 0.0
    FluxInDust2       = 0.0
    FluxInDust3       = 0.0
    FluxInDust4       = 0.0
    FluxInDust5       = 0.0
    FluxOutBChydropho = 0.0
    FluxOutBChydrophi = 0.0
    FluxOutOChydropho = 0.0
    FluxOutOChydrophi = 0.0
    FluxOutDust1      = 0.0
    FluxOutDust2      = 0.0
    FluxOutDust3      = 0.0
    FluxOutDust4      = 0.0
    FluxOutDust5      = 0.0

    ! compute aerosol mass in snow for each layer from interlayer flux
    do LoopInd = -NumSnowLayerMax+1, 0

       SnowMass = SnowLiqWater(LoopInd) + SnowIce(LoopInd)

       if (LoopInd >= NumSnowLayerNeg+1 .and. SnowMass > 0.0) then

          MassBChydropho(LoopInd) =  MassBChydropho (LoopInd) + FluxInBChydropho * MainTimeStep
          MassBChydrophi(LoopInd) =  MassBChydrophi (LoopInd) + FluxInBChydrophi * MainTimeStep
          MassOChydropho(LoopInd) =  MassOChydropho (LoopInd) + FluxInOChydropho * MainTimeStep
          MassOChydrophi(LoopInd) =  MassOChydrophi (LoopInd) + FluxInOChydrophi * MainTimeStep
          MassDust1(LoopInd)      =  MassDust1 (LoopInd)      + FluxInDust1 * MainTimeStep
          MassDust2(LoopInd)      =  MassDust2 (LoopInd)      + FluxInDust2 * MainTimeStep
          MassDust3(LoopInd)      =  MassDust3 (LoopInd)      + FluxInDust3 * MainTimeStep
          MassDust4(LoopInd)      =  MassDust4 (LoopInd)      + FluxInDust4 * MainTimeStep
          MassDust5(LoopInd)      =  MassDust5 (LoopInd)      + FluxInDust5 * MainTimeStep

          !BCPHO
          FluxOutBChydropho = OutflowSnowLayer(LoopInd) * ScavEffMeltScale * &
                              ScavEffMeltBCpho * (MassBChydropho(LoopInd) / SnowMass)
          if (FluxOutBChydropho * MainTimeStep > MassBChydropho(LoopInd)) then
             FluxOutBChydropho = MassBChydropho(LoopInd) / MainTimeStep
             MassBChydropho(LoopInd) = 0.0
          else
             MassBChydropho(LoopInd) = MassBChydropho(LoopInd) - FluxOutBChydropho * MainTimeStep
          end if
          FluxInBChydropho = FluxOutBChydropho

          !BCPHI
          FluxOutBChydrophi = OutflowSnowLayer(LoopInd) * ScavEffMeltScale * &
                              ScavEffMeltBCphi * (MassBChydrophi(LoopInd) / SnowMass)
          if (FluxOutBChydrophi * MainTimeStep > MassBChydrophi(LoopInd)) then
             FluxOutBChydrophi = MassBChydrophi(LoopInd) / MainTimeStep
             MassBChydrophi(LoopInd) = 0.0
          else
             MassBChydrophi(LoopInd) = MassBChydrophi(LoopInd) - FluxOutBChydrophi * MainTimeStep
          end if
          FluxInBChydrophi = FluxOutBChydrophi

          !OCPHO
          FluxOutOChydropho = OutflowSnowLayer(LoopInd) * ScavEffMeltScale * &
                              ScavEffMeltOCpho * (MassOChydropho(LoopInd) / SnowMass)
          if (FluxOutOChydropho * MainTimeStep > MassOChydropho(LoopInd)) then
             FluxOutOChydropho = MassOChydropho(LoopInd) / MainTimeStep
             MassOChydropho(LoopInd) = 0.0
          else
             MassOChydropho(LoopInd) = MassOChydropho(LoopInd) - FluxOutOChydropho * MainTimeStep
          end if
          FluxInOChydropho = FluxOutOChydropho

          !OCPHI
          FluxOutOChydrophi = OutflowSnowLayer(LoopInd) * ScavEffMeltScale * &
                              ScavEffMeltOCphi * (MassOChydrophi(LoopInd) / SnowMass)
          if (FluxOutOChydrophi * MainTimeStep > MassOChydrophi(LoopInd)) then
             FluxOutOChydrophi = MassOChydrophi(LoopInd) / MainTimeStep
             MassOChydrophi(LoopInd) = 0.0
          else
             MassOChydrophi(LoopInd) = MassOChydrophi(LoopInd) - FluxOutOChydrophi * MainTimeStep
          end if
          FluxInOChydrophi = FluxOutOChydrophi

          !Dust 1
          FluxOutDust1 = OutflowSnowLayer(LoopInd) * ScavEffMeltScale * &
                         ScavEffMeltDust1 * (MassDust1(LoopInd) / SnowMass)
          if (FluxOutDust1 * MainTimeStep > MassDust1(LoopInd)) then
             FluxOutDust1 = MassDust1(LoopInd) / MainTimeStep
             MassDust1(LoopInd) = 0.0
          else
             MassDust1(LoopInd) = MassDust1(LoopInd) - FluxOutDust1 * MainTimeStep
          end if
          FluxInDust1 = FluxOutDust1

          !Dust 2
          FluxOutDust2 = OutflowSnowLayer(LoopInd) * ScavEffMeltScale * &
                         ScavEffMeltDust2 * (MassDust2(LoopInd) / SnowMass)
          if (FluxOutDust2 * MainTimeStep > MassDust2(LoopInd)) then
             FluxOutDust2 = MassDust2(LoopInd) / MainTimeStep
             MassDust2(LoopInd) = 0.0
          else
             MassDust2(LoopInd) = MassDust2(LoopInd) - FluxOutDust2 * MainTimeStep
          end if
          FluxInDust2 = FluxOutDust2

          !Dust 3
          FluxOutDust3 = OutflowSnowLayer(LoopInd) * ScavEffMeltScale * &
                         ScavEffMeltDust3 * (MassDust3(LoopInd) / SnowMass)
          if (FluxOutDust3 * MainTimeStep > MassDust3(LoopInd)) then
             FluxOutDust3 = MassDust3(LoopInd) / MainTimeStep
             MassDust3(LoopInd) = 0.0
          else
             MassDust3(LoopInd) = MassDust3(LoopInd) - FluxOutDust3 * MainTimeStep
          end if
          FluxInDust3 = FluxOutDust3

          !Dust 4
          FluxOutDust4 = OutflowSnowLayer(LoopInd) * ScavEffMeltScale * &
                         ScavEffMeltDust4 * (MassDust4(LoopInd) / SnowMass)
          if (FluxOutDust4 * MainTimeStep > MassDust4(LoopInd)) then
             FluxOutDust4 = MassDust4(LoopInd) / MainTimeStep
             MassDust4(LoopInd) = 0.0
          else
             MassDust4(LoopInd) = MassDust4(LoopInd) - FluxOutDust4 * MainTimeStep
          end if
          FluxInDust4 = FluxOutDust4

          !Dust 5
          FluxOutDust5 = OutflowSnowLayer(LoopInd) * ScavEffMeltScale * &
                         ScavEffMeltDust5 * (MassDust5(LoopInd) / SnowMass)
          if (FluxOutDust5 * MainTimeStep > MassDust5(LoopInd)) then
             FluxOutDust5 = MassDust5(LoopInd) / MainTimeStep
             MassDust5(LoopInd) = 0.0
          else
             MassDust5(LoopInd) = MassDust5(LoopInd) - FluxOutDust5 * MainTimeStep
          end if
          FluxInDust5 = FluxOutDust5

       else ! SnowMass <=0 or non-existence snow layer

          MassBChydropho(LoopInd) = 0.0
          MassBChydrophi(LoopInd) = 0.0
          MassOChydropho(LoopInd) = 0.0
          MassOChydrophi(LoopInd) = 0.0
          MassDust1(LoopInd)      = 0.0
          MassDust2(LoopInd)      = 0.0
          MassDust3(LoopInd)      = 0.0
          MassDust4(LoopInd)      = 0.0
          MassDust5(LoopInd)      = 0.0
          FluxInBChydropho        = 0.0
          FluxInBChydrophi        = 0.0
          FluxInOChydropho        = 0.0
          FluxInOChydrophi        = 0.0
          FluxInDust1             = 0.0
          FluxInDust2             = 0.0
          FluxInDust3             = 0.0
          FluxInDust4             = 0.0
          FluxInDust5             = 0.0
          FluxOutBChydropho       = 0.0
          FluxOutBChydrophi       = 0.0
          FluxOutOChydropho       = 0.0
          FluxOutOChydrophi       = 0.0
          FluxOutDust1            = 0.0
          FluxOutDust2            = 0.0
          FluxOutDust3            = 0.0
          FluxOutDust4            = 0.0
          FluxOutDust5            = 0.0

       endif
    enddo

    ! update aerosol mass for the top snow layer from atmos deposition flux
    if (NumSnowLayerNeg < 0 ) then
       MassBChydropho(NumSnowLayerNeg+1) =  MassBChydropho (NumSnowLayerNeg+1) +  DepBChydropho * MainTimeStep
       MassBChydrophi(NumSnowLayerNeg+1) =  MassBChydrophi (NumSnowLayerNeg+1) +  DepBChydrophi * MainTimeStep
       MassOChydropho(NumSnowLayerNeg+1) =  MassOChydropho (NumSnowLayerNeg+1) +  DepOChydropho * MainTimeStep
       MassOChydrophi(NumSnowLayerNeg+1) =  MassOChydrophi (NumSnowLayerNeg+1) +  DepOChydrophi * MainTimeStep
       MassDust1(NumSnowLayerNeg+1)      =  MassDust1 (NumSnowLayerNeg+1)      +  DepDust1 * MainTimeStep
       MassDust2(NumSnowLayerNeg+1)      =  MassDust2 (NumSnowLayerNeg+1)      +  DepDust2 * MainTimeStep
       MassDust3(NumSnowLayerNeg+1)      =  MassDust3 (NumSnowLayerNeg+1)      +  DepDust3 * MainTimeStep
       MassDust4(NumSnowLayerNeg+1)      =  MassDust4 (NumSnowLayerNeg+1)      +  DepDust4 * MainTimeStep
       MassDust5(NumSnowLayerNeg+1)      =  MassDust5 (NumSnowLayerNeg+1)      +  DepDust5 * MainTimeStep
    endif

    ! update aerosol mass concentration in snow for each layer
    do LoopInd = -NumSnowLayerMax+1, 0

       SnowMass = SnowLiqWater(LoopInd) + SnowIce(LoopInd) 

       if (LoopInd >= NumSnowLayerNeg+1 .and. SnowMass > 0.0) then
          MassConcBChydropho(LoopInd) =  MassBChydropho(LoopInd) / SnowMass
          MassConcBChydrophi(LoopInd) =  MassBChydrophi(LoopInd) / SnowMass
          MassConcOChydropho(LoopInd) =  MassOChydropho(LoopInd) / SnowMass
          MassConcOChydrophi(LoopInd) =  MassOChydrophi(LoopInd) / SnowMass
          MassConcDust1(LoopInd)      =  MassDust1(LoopInd) / SnowMass
          MassConcDust2(LoopInd)      =  MassDust2(LoopInd) / SnowMass
          MassConcDust3(LoopInd)      =  MassDust3(LoopInd) / SnowMass
          MassConcDust4(LoopInd)      =  MassDust4(LoopInd) / SnowMass
          MassConcDust5(LoopInd)      =  MassDust5(LoopInd) / SnowMass
       else
          MassConcBChydropho(LoopInd) =  0.0
          MassConcBChydrophi(LoopInd) =  0.0
          MassConcOChydropho(LoopInd) =  0.0
          MassConcOChydrophi(LoopInd) =  0.0
          MassConcDust1(LoopInd)      =  0.0
          MassConcDust2(LoopInd)      =  0.0
          MassConcDust3(LoopInd)      =  0.0
          MassConcDust4(LoopInd)      =  0.0
          MassConcDust5(LoopInd)      =  0.0
          MassBChydropho(LoopInd)     =  0.0
          MassBChydrophi(LoopInd)     =  0.0
          MassOChydropho(LoopInd)     =  0.0
          MassOChydrophi(LoopInd)     =  0.0
          MassDust1(LoopInd)          =  0.0
          MassDust2(LoopInd)          =  0.0
          MassDust3(LoopInd)          =  0.0
          MassDust4(LoopInd)          =  0.0
          MassDust5(LoopInd)          =  0.0
       endif

    enddo

    ! special treatment for very shallow snowpack (NumSnowLayerNeg = 0 and SnowMass > 0.0)
    if ( NumSnowLayerNeg == 0 ) then

       SnowMass = SnowWaterEquiv

       if ( SnowMass > 0.1 ) then ! set minimum threshold (0.1 mm SWE) for computing aerosol-snow albedo
          MassBChydropho(0)     =  MassBChydropho (0) +  DepBChydropho * MainTimeStep
          MassBChydrophi(0)     =  MassBChydrophi (0) +  DepBChydrophi * MainTimeStep
          MassOChydropho(0)     =  MassOChydropho (0) +  DepOChydropho * MainTimeStep
          MassOChydrophi(0)     =  MassOChydrophi (0) +  DepOChydrophi * MainTimeStep
          MassDust1(0)          =  MassDust1 (0)      +  DepDust1 * MainTimeStep
          MassDust2(0)          =  MassDust2 (0)      +  DepDust2 * MainTimeStep
          MassDust3(0)          =  MassDust3 (0)      +  DepDust3 * MainTimeStep
          MassDust4(0)          =  MassDust4 (0)      +  DepDust4 * MainTimeStep
          MassDust5(0)          =  MassDust5 (0)      +  DepDust5 * MainTimeStep
          MassConcBChydropho(0) =  MassBChydropho(0) / SnowMass
          MassConcBChydrophi(0) =  MassBChydrophi(0) / SnowMass
          MassConcOChydropho(0) =  MassOChydropho(0) / SnowMass
          MassConcOChydrophi(0) =  MassOChydrophi(0) / SnowMass
          MassConcDust1(0)      =  MassDust1(0) / SnowMass
          MassConcDust2(0)      =  MassDust2(0) / SnowMass
          MassConcDust3(0)      =  MassDust3(0) / SnowMass
          MassConcDust4(0)      =  MassDust4(0) / SnowMass
          MassConcDust5(0)      =  MassDust5(0) / SnowMass
       else
          MassBChydropho(0)     =  0.0
          MassBChydrophi(0)     =  0.0
          MassOChydropho(0)     =  0.0
          MassOChydrophi(0)     =  0.0
          MassDust1(0)          =  0.0
          MassDust2(0)          =  0.0
          MassDust3(0)          =  0.0
          MassDust4(0)          =  0.0
          MassDust5(0)          =  0.0
          MassConcBChydropho(0) =  0.0
          MassConcBChydrophi(0) =  0.0
          MassConcOChydropho(0) =  0.0
          MassConcOChydrophi(0) =  0.0
          MassConcDust1(0)      =  0.0
          MassConcDust2(0)      =  0.0
          MassConcDust3(0)      =  0.0
          MassConcDust4(0)      =  0.0
          MassConcDust5(0)      =  0.0
       endif

    endif

    end associate

  end subroutine SnowAerosolSnicar

end module SnowAerosolSnicarMod
