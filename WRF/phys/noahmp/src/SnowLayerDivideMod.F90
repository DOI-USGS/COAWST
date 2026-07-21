module SnowLayerDivideMod

!!! Snowpack layer division process
!!! Update snow ice, snow water, snow thickness, snow temperature

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowLayerWaterComboMod, only: SnowLayerWaterCombo

  implicit none

contains

  subroutine SnowLayerDivide(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: DIVIDE
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd                              ! snow layer loop index
    integer                          :: NumSnowLayerTmp                      ! number of snow layer top to bottom
    real(kind=kind_noahmp)           :: SnowThickCombTmp                     ! thickness of the combined [m]
    real(kind=kind_noahmp)           :: SnowIceExtra                         ! extra snow ice to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: SnowLiqExtra                         ! extra snow liquid water to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: SnowFracExtra                        ! fraction of extra snow to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: SnowTempGrad                         ! temperature gradient between two snow layers
    real(kind=kind_noahmp)           :: MassBChydrophoExtra                  ! extra mass of hydrophobic BC in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassBChydrophiExtra                  ! extra mass of hydrophillic BC in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassOChydrophoExtra                  ! extra mass of hydrophobic OC in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassOChydrophiExtra                  ! extra mass of hydrophillic OC in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust1Extra                       ! extra mass of dust species 1 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust2Extra                       ! extra mass of dust species 2 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust3Extra                       ! extra mass of dust species 3 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust4Extra                       ! extra mass of dust species 4 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp)           :: MassDust5Extra                       ! extra mass of dust species 5 in snow [kg m-2] to be divided compared to allowed layer thickness
    real(kind=kind_noahmp), allocatable, dimension(:) :: SnowThickTmp        ! snow layer thickness [m]
    real(kind=kind_noahmp), allocatable, dimension(:) :: SnowIceTmp          ! partial volume of ice [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:) :: SnowLiqTmp          ! partial volume of liquid water [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:) :: TemperatureSnowTmp  ! node temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:) :: MassBChydrophoTmp   ! mass of hydrophobic Black Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:) :: MassBChydrophiTmp   ! mass of hydrophillic Black Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:) :: MassOChydrophoTmp   ! mass of hydrophobic Organic Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:) :: MassOChydrophiTmp   ! mass of hydrophillic Organic Carbon in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:) :: MassDust1Tmp        ! mass of dust species 1 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:) :: MassDust2Tmp        ! mass of dust species 2 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:) :: MassDust3Tmp        ! mass of dust species 3 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:) :: MassDust4Tmp        ! mass of dust species 4 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:) :: MassDust5Tmp        ! mass of dust species 5 in snow [kg m-2]
    real(kind=kind_noahmp), allocatable, dimension(:) :: SnowRadiusTmp       ! effective grain radius [microns, m-6]

! --------------------------------------------------------------------
    associate(                                                                       &
              OptSnowAlbedo          => noahmp%config%nmlist%OptSnowAlbedo          ,& ! in,    options for ground snow surface albedo
              NumSnowLayerMax        => noahmp%config%domain%NumSnowLayerMax        ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg        ,& ! inout, actual number of snow layers (negative)
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow     ,& ! inout, snow and soil layer temperature [K]
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
              SnowRadius             => noahmp%water%state%SnowRadius                & ! inout, effective grain radius [microns, m-6]
             )
! ----------------------------------------------------------------------

    ! initialization
    if (.not. allocated(SnowIceTmp)        ) allocate(SnowIceTmp        (1:NumSnowLayerMax))
    if (.not. allocated(SnowLiqTmp)        ) allocate(SnowLiqTmp        (1:NumSnowLayerMax))
    if (.not. allocated(TemperatureSnowTmp)) allocate(TemperatureSnowTmp(1:NumSnowLayerMax))
    if (.not. allocated(SnowThickTmp)      ) allocate(SnowThickTmp      (1:NumSnowLayerMax))

    if ( OptSnowAlbedo == 3 ) then
       if (.not. allocated(MassBChydrophoTmp)) allocate(MassBChydrophoTmp (1:NumSnowLayerMax))
       if (.not. allocated(MassBChydrophiTmp)) allocate(MassBChydrophiTmp (1:NumSnowLayerMax))
       if (.not. allocated(MassOChydrophoTmp)) allocate(MassOChydrophoTmp (1:NumSnowLayerMax))
       if (.not. allocated(MassOChydrophiTmp)) allocate(MassOChydrophiTmp (1:NumSnowLayerMax))
       if (.not. allocated(MassDust1Tmp)     ) allocate(MassDust1Tmp      (1:NumSnowLayerMax))
       if (.not. allocated(MassDust2Tmp)     ) allocate(MassDust2Tmp      (1:NumSnowLayerMax))
       if (.not. allocated(MassDust3Tmp)     ) allocate(MassDust3Tmp      (1:NumSnowLayerMax))
       if (.not. allocated(MassDust4Tmp)     ) allocate(MassDust4Tmp      (1:NumSnowLayerMax))
       if (.not. allocated(MassDust5Tmp)     ) allocate(MassDust5Tmp      (1:NumSnowLayerMax))
       if (.not. allocated(SnowRadiusTmp)    ) allocate(SnowRadiusTmp     (1:NumSnowLayerMax))
    endif

    SnowIceTmp        (:) = 0.0
    SnowLiqTmp        (:) = 0.0
    TemperatureSnowTmp(:) = 0.0
    SnowThickTmp      (:) = 0.0

    if ( OptSnowAlbedo == 3 ) then
       MassBChydrophoTmp(:) = 0.0
       MassBChydrophiTmp(:) = 0.0
       MassOChydrophoTmp(:) = 0.0
       MassOChydrophiTmp(:) = 0.0
       MassDust1Tmp     (:) = 0.0
       MassDust2Tmp     (:) = 0.0
       MassDust3Tmp     (:) = 0.0
       MassDust4Tmp     (:) = 0.0
       MassDust5Tmp     (:) = 0.0
       SnowRadiusTmp    (:) = 0.0
    endif

    do LoopInd = 1, NumSnowLayerMax
       if ( LoopInd <= abs(NumSnowLayerNeg) ) then
          SnowThickTmp(LoopInd)       = ThicknessSnowSoilLayer(LoopInd+NumSnowLayerNeg)
          SnowIceTmp(LoopInd)         = SnowIce(LoopInd+NumSnowLayerNeg)
          SnowLiqTmp(LoopInd)         = SnowLiqWater(LoopInd+NumSnowLayerNeg)
          TemperatureSnowTmp(LoopInd) = TemperatureSoilSnow(LoopInd+NumSnowLayerNeg)

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(LoopInd) = MassBChydropho(LoopInd+NumSnowLayerNeg)
             MassBChydrophiTmp(LoopInd) = MassBChydrophi(LoopInd+NumSnowLayerNeg)
             MassOChydrophoTmp(LoopInd) = MassOChydropho(LoopInd+NumSnowLayerNeg)
             MassOChydrophiTmp(LoopInd) = MassOChydrophi(LoopInd+NumSnowLayerNeg)
             MassDust1Tmp(LoopInd)      = MassDust1(LoopInd+NumSnowLayerNeg)
             MassDust2Tmp(LoopInd)      = MassDust2(LoopInd+NumSnowLayerNeg)
             MassDust3Tmp(LoopInd)      = MassDust3(LoopInd+NumSnowLayerNeg)
             MassDust4Tmp(LoopInd)      = MassDust4(LoopInd+NumSnowLayerNeg)
             MassDust5Tmp(LoopInd)      = MassDust5(LoopInd+NumSnowLayerNeg)
             SnowRadiusTmp(LoopInd)     = SnowRadius(LoopInd+NumSnowLayerNeg)
          endif
       endif
    enddo

    ! start snow layer division
    NumSnowLayerTmp = abs(NumSnowLayerNeg)

    if ( NumSnowLayerTmp == 1 ) then
       ! Specify a new snow layer
       if ( SnowThickTmp(1) > 0.05 ) then
          NumSnowLayerTmp       = 2
          SnowThickTmp(1)       = SnowThickTmp(1)/2.0
          SnowIceTmp(1)         = SnowIceTmp(1)/2.0
          SnowLiqTmp(1)         = SnowLiqTmp(1)/2.0
          SnowThickTmp(2)       = SnowThickTmp(1)
          SnowIceTmp(2)         = SnowIceTmp(1)
          SnowLiqTmp(2)         = SnowLiqTmp(1)
          TemperatureSnowTmp(2) = TemperatureSnowTmp(1)
         
          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(1) = MassBChydrophoTmp(1)/2.0
             MassBChydrophoTmp(2) = MassBChydrophoTmp(1)
             MassBChydrophiTmp(1) = MassBChydrophiTmp(1)/2.0
             MassBChydrophiTmp(2) = MassBChydrophiTmp(1)
             MassOChydrophoTmp(1) = MassOChydrophoTmp(1)/2.0
             MassOChydrophoTmp(2) = MassOChydrophoTmp(1)
             MassOChydrophiTmp(1) = MassOChydrophiTmp(1)/2.0
             MassOChydrophiTmp(2) = MassOChydrophiTmp(1)
             MassDust1Tmp(1)      = MassDust1Tmp(1)/2.0
             MassDust1Tmp(2)      = MassDust1Tmp(1)
             MassDust2Tmp(1)      = MassDust2Tmp(1)/2.0
             MassDust2Tmp(2)      = MassDust2Tmp(1)
             MassDust3Tmp(1)      = MassDust3Tmp(1)/2.0
             MassDust3Tmp(2)      = MassDust3Tmp(1)
             MassDust4Tmp(1)      = MassDust4Tmp(1)/2.0
             MassDust4Tmp(2)      = MassDust4Tmp(1)
             MassDust5Tmp(1)      = MassDust5Tmp(1)/2.0
             MassDust5Tmp(2)      = MassDust5Tmp(1)
             SnowRadiusTmp(2)     = SnowRadiusTmp(1)
          endif
       endif
    endif

    if ( NumSnowLayerTmp > 1 ) then
       if ( SnowThickTmp(1) > 0.05 ) then     ! maximum allowed thickness (5cm) for top snow layer
          SnowThickCombTmp     = SnowThickTmp(1) - 0.05
          SnowFracExtra        = SnowThickCombTmp / SnowThickTmp(1)
          SnowIceExtra         = SnowFracExtra * SnowIceTmp(1)
          SnowLiqExtra         = SnowFracExtra * SnowLiqTmp(1)

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoExtra = SnowFracExtra * MassBChydrophoTmp(1)
             MassBChydrophiExtra = SnowFracExtra * MassBChydrophiTmp(1)
             MassOChydrophoExtra = SnowFracExtra * MassOChydrophoTmp(1)
             MassOChydrophiExtra = SnowFracExtra * MassOChydrophiTmp(1)
             MassDust1Extra      = SnowFracExtra * MassDust1Tmp(1)
             MassDust2Extra      = SnowFracExtra * MassDust2Tmp(1)
             MassDust3Extra      = SnowFracExtra * MassDust3Tmp(1)
             MassDust4Extra      = SnowFracExtra * MassDust4Tmp(1)
             MassDust5Extra      = SnowFracExtra * MassDust5Tmp(1)
          endif

          SnowFracExtra        = 0.05 / SnowThickTmp(1)
          SnowIceTmp(1)        = SnowFracExtra * SnowIceTmp(1)
          SnowLiqTmp(1)        = SnowFracExtra * SnowLiqTmp(1)
          SnowThickTmp(1)      = 0.05

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(1) = SnowFracExtra * MassBChydrophoTmp(1)
             MassBChydrophiTmp(1) = SnowFracExtra * MassBChydrophiTmp(1)
             MassOChydrophoTmp(1) = SnowFracExtra * MassOChydrophoTmp(1)
             MassOChydrophiTmp(1) = SnowFracExtra * MassOChydrophiTmp(1)
             MassDust1Tmp(1)      = SnowFracExtra * MassDust1Tmp(1)
             MassDust2Tmp(1)      = SnowFracExtra * MassDust2Tmp(1)
             MassDust3Tmp(1)      = SnowFracExtra * MassDust3Tmp(1)
             MassDust4Tmp(1)      = SnowFracExtra * MassDust4Tmp(1)
             MassDust5Tmp(1)      = SnowFracExtra * MassDust5Tmp(1)

             MassBChydrophoTmp(2) = MassBChydrophoTmp(2) + MassBChydrophoExtra
             MassBChydrophiTmp(2) = MassBChydrophiTmp(2) + MassBChydrophiExtra
             MassOChydrophoTmp(2) = MassOChydrophoTmp(2) + MassOChydrophoExtra
             MassOChydrophiTmp(2) = MassOChydrophiTmp(2) + MassOChydrophiExtra
             MassDust1Tmp(2)      = MassDust1Tmp(2) + MassDust1Extra 
             MassDust2Tmp(2)      = MassDust2Tmp(2) + MassDust2Extra
             MassDust3Tmp(2)      = MassDust3Tmp(2) + MassDust3Extra
             MassDust4Tmp(2)      = MassDust4Tmp(2) + MassDust4Extra
             MassDust5Tmp(2)      = MassDust5Tmp(2) + MassDust5Extra
             SnowRadiusTmp(2)     = (SnowRadiusTmp(2)*(SnowLiqTmp(2)+SnowIceTmp(2))+SnowRadiusTmp(1)*(SnowLiqExtra+SnowIceExtra)) / &
                                    (SnowLiqTmp(2) + SnowIceTmp(2) + SnowLiqExtra + SnowIceExtra) 
          endif

          ! update combined snow water & temperature
          call SnowLayerWaterCombo(SnowThickTmp(2), SnowLiqTmp(2), SnowIceTmp(2), TemperatureSnowTmp(2), &
                                   SnowThickCombTmp, SnowLiqExtra, SnowIceExtra, TemperatureSnowTmp(1))

          ! subdivide a new layer, maximum allowed thickness (20cm) for second snow layer
          if ( (NumSnowLayerTmp <= 2) .and. (SnowThickTmp(2) > 0.20) ) then  ! MB: change limit
         !if ( (NumSnowLayerTmp <= 2) .and. (SnowThickTmp(2) > 0.10) ) then
             NumSnowLayerTmp       = 3
             SnowTempGrad          = (TemperatureSnowTmp(1) - TemperatureSnowTmp(2)) / &
                                     ((SnowThickTmp(1)+SnowThickTmp(2)) / 2.0)
             SnowThickTmp(2)       = SnowThickTmp(2) / 2.0
             SnowIceTmp(2)         = SnowIceTmp(2) / 2.0
             SnowLiqTmp(2)         = SnowLiqTmp(2) / 2.0
             SnowThickTmp(3)       = SnowThickTmp(2)
             SnowIceTmp(3)         = SnowIceTmp(2)
             SnowLiqTmp(3)         = SnowLiqTmp(2)
             TemperatureSnowTmp(3) = TemperatureSnowTmp(2) - SnowTempGrad * SnowThickTmp(2) / 2.0
             if ( TemperatureSnowTmp(3) >= ConstFreezePoint ) then
                TemperatureSnowTmp(3) = TemperatureSnowTmp(2)
             else
                TemperatureSnowTmp(2) = TemperatureSnowTmp(2) + SnowTempGrad * SnowThickTmp(2) / 2.0
             endif

             if ( OptSnowAlbedo == 3 ) then
                MassBChydrophoTmp(2) = MassBChydrophoTmp(2) / 2.0
                MassBChydrophoTmp(3) = MassBChydrophoTmp(2)
                MassBChydrophiTmp(2) = MassBChydrophiTmp(2) / 2.0
                MassBChydrophiTmp(3) = MassBChydrophiTmp(2)
                MassOChydrophoTmp(2) = MassOChydrophoTmp(2) / 2.0
                MassOChydrophoTmp(3) = MassOChydrophoTmp(2)
                MassOChydrophiTmp(2) = MassOChydrophiTmp(2) / 2.0
                MassOChydrophiTmp(3) = MassOChydrophiTmp(2)
                MassDust1Tmp(2)      = MassDust1Tmp(2) / 2.0
                MassDust1Tmp(3)      = MassDust1Tmp(2)
                MassDust2Tmp(2)      = MassDust2Tmp(2) / 2.0
                MassDust2Tmp(3)      = MassDust2Tmp(2)
                MassDust3Tmp(2)      = MassDust3Tmp(2) / 2.0
                MassDust3Tmp(3)      = MassDust3Tmp(2)
                MassDust4Tmp(2)      = MassDust4Tmp(2) / 2.0
                MassDust4Tmp(3)      = MassDust4Tmp(2)
                MassDust5Tmp(2)      = MassDust5Tmp(2) / 2.0
                MassDust5Tmp(3)      = MassDust5Tmp(2)
                SnowRadiusTmp(3)     = SnowRadiusTmp(2)
             endif

          endif
       endif ! if(SnowThickTmp(1) > 0.05)
    endif  ! if (NumSnowLayerTmp > 1)

    if ( NumSnowLayerTmp > 2 ) then
       if ( SnowThickTmp(2) > 0.2 ) then
          SnowThickCombTmp = SnowThickTmp(2) - 0.2
          SnowFracExtra    = SnowThickCombTmp / SnowThickTmp(2)
          SnowIceExtra     = SnowFracExtra * SnowIceTmp(2)
          SnowLiqExtra     = SnowFracExtra * SnowLiqTmp(2)

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoExtra = SnowFracExtra * MassBChydrophoTmp(2)
             MassBChydrophiExtra = SnowFracExtra * MassBChydrophiTmp(2)
             MassOChydrophoExtra = SnowFracExtra * MassOChydrophoTmp(2)
             MassOChydrophiExtra = SnowFracExtra * MassOChydrophiTmp(2)
             MassDust1Extra      = SnowFracExtra * MassDust1Tmp(2)
             MassDust2Extra      = SnowFracExtra * MassDust2Tmp(2)
             MassDust3Extra      = SnowFracExtra * MassDust3Tmp(2)
             MassDust4Extra      = SnowFracExtra * MassDust4Tmp(2)
             MassDust5Extra      = SnowFracExtra * MassDust5Tmp(2)
          endif

          SnowFracExtra    = 0.2 / SnowThickTmp(2)
          SnowIceTmp(2)    = SnowFracExtra * SnowIceTmp(2)
          SnowLiqTmp(2)    = SnowFracExtra * SnowLiqTmp(2)
          SnowThickTmp(2)  = 0.2

          if ( OptSnowAlbedo == 3 ) then
             MassBChydrophoTmp(2) = SnowFracExtra * MassBChydrophoTmp(2)
             MassBChydrophiTmp(2) = SnowFracExtra * MassBChydrophiTmp(2)
             MassOChydrophoTmp(2) = SnowFracExtra * MassOChydrophoTmp(2)
             MassOChydrophiTmp(2) = SnowFracExtra * MassOChydrophiTmp(2)
             MassDust1Tmp(2)      = SnowFracExtra * MassDust1Tmp(2)
             MassDust2Tmp(2)      = SnowFracExtra * MassDust2Tmp(2)
             MassDust3Tmp(2)      = SnowFracExtra * MassDust3Tmp(2)
             MassDust4Tmp(2)      = SnowFracExtra * MassDust4Tmp(2)
             MassDust5Tmp(2)      = SnowFracExtra * MassDust5Tmp(2)

             MassBChydrophoTmp(3) = MassBChydrophoTmp(3) + MassBChydrophoExtra
             MassBChydrophiTmp(3) = MassBChydrophiTmp(3) + MassBChydrophiExtra
             MassOChydrophoTmp(3) = MassOChydrophoTmp(3) + MassOChydrophoExtra
             MassOChydrophiTmp(3) = MassOChydrophiTmp(3) + MassOChydrophiExtra
             MassDust1Tmp(3)      = MassDust1Tmp(3) + MassDust1Extra
             MassDust2Tmp(3)      = MassDust2Tmp(3) + MassDust2Extra
             MassDust3Tmp(3)      = MassDust3Tmp(3) + MassDust3Extra
             MassDust4Tmp(3)      = MassDust4Tmp(3) + MassDust4Extra
             MassDust5Tmp(3)      = MassDust5Tmp(3) + MassDust5Extra
             SnowRadiusTmp(3)     = (SnowRadiusTmp(3)*(SnowLiqTmp(3)+SnowIceTmp(3))+SnowRadiusTmp(2)*(SnowLiqExtra+SnowIceExtra)) / &
                                    (SnowLiqTmp(3) + SnowIceTmp(3) + SnowLiqExtra + SnowIceExtra)
          endif

          ! update combined snow water & temperature
          call SnowLayerWaterCombo(SnowThickTmp(3), SnowLiqTmp(3), SnowIceTmp(3), TemperatureSnowTmp(3), &
                                   SnowThickCombTmp, SnowLiqExtra, SnowIceExtra, TemperatureSnowTmp(2))
       endif
    endif

    NumSnowLayerNeg = -NumSnowLayerTmp

    do LoopInd = NumSnowLayerNeg+1, 0
       ThicknessSnowSoilLayer(LoopInd) = SnowThickTmp(LoopInd-NumSnowLayerNeg)
       SnowIce(LoopInd)                = SnowIceTmp(LoopInd-NumSnowLayerNeg)
       SnowLiqWater(LoopInd)           = SnowLiqTmp(LoopInd-NumSnowLayerNeg)
       TemperatureSoilSnow(LoopInd)    = TemperatureSnowTmp(LoopInd-NumSnowLayerNeg)

       if ( OptSnowAlbedo == 3 ) then
          MassBChydropho(LoopInd)      = MassBChydrophoTmp(LoopInd-NumSnowLayerNeg)
          MassBChydrophi(LoopInd)      = MassBChydrophiTmp(LoopInd-NumSnowLayerNeg)
          MassOChydropho(LoopInd)      = MassOChydrophoTmp(LoopInd-NumSnowLayerNeg)
          MassOChydrophi(LoopInd)      = MassOChydrophiTmp(LoopInd-NumSnowLayerNeg)
          MassDust1(LoopInd)           = MassDust1Tmp(LoopInd-NumSnowLayerNeg)
          MassDust2(LoopInd)           = MassDust2Tmp(LoopInd-NumSnowLayerNeg)
          MassDust3(LoopInd)           = MassDust3Tmp(LoopInd-NumSnowLayerNeg)
          MassDust4(LoopInd)           = MassDust4Tmp(LoopInd-NumSnowLayerNeg)
          MassDust5(LoopInd)           = MassDust5Tmp(LoopInd-NumSnowLayerNeg)
          SnowRadius(LoopInd)          = SnowRadiusTmp(LoopInd-NumSnowLayerNeg)
       endif
    enddo

    ! deallocate local arrays to avoid memory leaks
    deallocate(SnowIceTmp        )
    deallocate(SnowLiqTmp        )
    deallocate(TemperatureSnowTmp)
    deallocate(SnowThickTmp      )

    if ( OptSnowAlbedo == 3 ) then
       deallocate(MassBChydrophoTmp)
       deallocate(MassBChydrophiTmp)
       deallocate(MassOChydrophoTmp)
       deallocate(MassOChydrophiTmp)
       deallocate(MassDust1Tmp     )
       deallocate(MassDust2Tmp     )
       deallocate(MassDust3Tmp     )
       deallocate(MassDust4Tmp     )
       deallocate(MassDust5Tmp     )
       deallocate(SnowRadiusTmp    )
    endif

    end associate

  end subroutine SnowLayerDivide

end module SnowLayerDivideMod
