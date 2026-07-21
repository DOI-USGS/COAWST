module SnowLayerCombineMod

!!! Snowpack layer combination process
!!! Update snow ice, snow water, snow thickness, snow temperature

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowLayerWaterComboMod, only: SnowLayerWaterCombo

  implicit none

contains

  subroutine SnowLayerCombine(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: COMBINE
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: I,J,K,L           ! node indices
    integer                          :: NumSnowLayerOld   ! number of snow layer
    integer                          :: IndLayer          ! node index
    integer                          :: IndNeighbor       ! adjacent node selected for combination
    real(kind=kind_noahmp)           :: SnowIceTmp        ! total ice mass in snow
    real(kind=kind_noahmp)           :: SnowLiqTmp        ! total liquid water in snow
    real(kind=kind_noahmp)           :: SnowThickMin(3)   ! minimum thickness of each snow layer
    data SnowThickMin /0.025, 0.025, 0.1/                 ! MB: change limit
    !data SnowThickMin /0.045, 0.05, 0.2/

! --------------------------------------------------------------------
    associate(                                                                       &
              OptSnowAlbedo          => noahmp%config%nmlist%OptSnowAlbedo          ,& ! in,    options for ground snow surface albedo
              NumSnowLayerNeg        => noahmp%config%domain%NumSnowLayerNeg        ,& ! inout, actual number of snow layers (negative)
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! inout, thickness of snow/soil layers [m]
              TemperatureSoilSnow    => noahmp%energy%state%TemperatureSoilSnow     ,& ! inout, snow and soil layer temperature [K]
              SnowDepth              => noahmp%water%state%SnowDepth                ,& ! inout, snow depth [m]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv           ,& ! inout, snow water equivalent [mm]
              SnowIce                => noahmp%water%state%SnowIce                  ,& ! inout, snow layer ice [mm]
              SnowLiqWater           => noahmp%water%state%SnowLiqWater             ,& ! inout, snow layer liquid water [mm]
              SoilLiqWater           => noahmp%water%state%SoilLiqWater             ,& ! inout, soil liquid moisture [m3/m3]
              SoilIce                => noahmp%water%state%SoilIce                  ,& ! inout, soil ice moisture [m3/m3]
              MassBChydropho         => noahmp%water%state%MassBChydropho           ,& ! inout, mass of hydrophobic Black Carbon in snow [kg m-2]
              MassBChydrophi         => noahmp%water%state%MassBChydrophi           ,& ! inout, mass of hydrophillic Black Carbon in snow [kg m-2]
              MassOChydropho         => noahmp%water%state%MassOChydropho           ,& ! inout, mass of hydrophobic Organic Carbon in snow [kg m-2]
              MassOChydrophi         => noahmp%water%state%MassOChydrophi           ,& ! inout, mass of hydrophillic Organic Carbon in snow [kg m-2]
              MassDust1              => noahmp%water%state%MassDust1                ,& ! inout, mass of dust species 1 in snow [kg m-2]
              MassDust2              => noahmp%water%state%MassDust2                ,& ! inout, mass of dust species 2 in snow [kg m-2]
              MassDust3              => noahmp%water%state%MassDust3                ,& ! inout, mass of dust species 3 in snow [kg m-2]
              MassDust4              => noahmp%water%state%MassDust4                ,& ! inout, mass of dust species 4 in snow [kg m-2]
              MassDust5              => noahmp%water%state%MassDust5                ,& ! inout, mass of dust species 5 in snow [kg m-2]
              SnowRadius             => noahmp%water%state%SnowRadius               ,& ! inout, effective grain radius [microns, m-6]
              PondSfcThinSnwComb     => noahmp%water%state%PondSfcThinSnwComb       ,& ! out,   surface ponding [mm] from liquid in thin snow layer combination
              PondSfcThinSnwTrans    => noahmp%water%state%PondSfcThinSnwTrans       & ! out,   surface ponding [mm] from thin snow when changing from multilayer to no layer
             )
! ----------------------------------------------------------------------

    ! check and combine small ice content layer
    NumSnowLayerOld = NumSnowLayerNeg

    do J = NumSnowLayerOld+1,0
       if ( SnowIce(J) <= 0.1 ) then
          if ( J /= 0 ) then
             SnowLiqWater(J+1)           = SnowLiqWater(J+1) + SnowLiqWater(J)
             SnowIce(J+1)                = SnowIce(J+1) + SnowIce(J)
             ThicknessSnowSoilLayer(J+1) = ThicknessSnowSoilLayer(J+1) + ThicknessSnowSoilLayer(J)

             if ( OptSnowAlbedo == 3 ) then
                MassBChydropho(J+1)      = MassBChydropho(J+1) +  MassBChydropho(J)
                MassBChydrophi(J+1)      = MassBChydrophi(J+1) +  MassBChydrophi(J)
                MassOChydropho(J+1)      = MassOChydropho(J+1) +  MassOChydropho(J)
                MassOChydrophi(J+1)      = MassOChydrophi(J+1) +  MassOChydrophi(J)
                MassDust1(J+1)           = MassDust1(J+1) +  MassDust1(J)
                MassDust2(J+1)           = MassDust2(J+1) +  MassDust2(J)
                MassDust3(J+1)           = MassDust3(J+1) +  MassDust3(J)
                MassDust4(J+1)           = MassDust4(J+1) +  MassDust4(J)
                MassDust5(J+1)           = MassDust5(J+1) +  MassDust5(J)
             endif

          else
             if ( NumSnowLayerNeg < -1 ) then    ! MB/KM: change to NumSnowLayerNeg
                SnowLiqWater(J-1)           = SnowLiqWater(J-1) + SnowLiqWater(J)
                SnowIce(J-1)                = SnowIce(J-1) + SnowIce(J)
                ThicknessSnowSoilLayer(J-1) = ThicknessSnowSoilLayer(J-1) + ThicknessSnowSoilLayer(J)

                if ( OptSnowAlbedo == 3 ) then
                   MassBChydropho(J-1)      = MassBChydropho(J-1) +  MassBChydropho(J)
                   MassBChydrophi(J-1)      = MassBChydrophi(J-1) +  MassBChydrophi(J)
                   MassOChydropho(J-1)      = MassOChydropho(J-1) +  MassOChydropho(J)
                   MassOChydrophi(J-1)      = MassOChydrophi(J-1) +  MassOChydrophi(J)
                   MassDust1(J-1)           = MassDust1(J-1) +  MassDust1(J)
                   MassDust2(J-1)           = MassDust2(J-1) +  MassDust2(J)
                   MassDust3(J-1)           = MassDust3(J-1) +  MassDust3(J)
                   MassDust4(J-1)           = MassDust4(J-1) +  MassDust4(J)
                   MassDust5(J-1)           = MassDust5(J-1) +  MassDust5(J)
                endif

             else
                if ( SnowIce(J) >= 0.0 ) then
                   PondSfcThinSnwComb = SnowLiqWater(J)                ! NumSnowLayerNeg WILL GET SET TO ZERO BELOW; PondSfcThinSnwComb WILL GET 
                   SnowWaterEquiv     = SnowIce(J)                     ! ADDED TO PONDING FROM PHASECHANGE PONDING SHOULD BE
                   SnowDepth          = ThicknessSnowSoilLayer(J)      ! ZERO HERE BECAUSE IT WAS CALCULATED FOR THIN SNOW
                else  ! SnowIce OVER-SUBLIMATED EARLIER
                   PondSfcThinSnwComb = SnowLiqWater(J) + SnowIce(J)
                   if ( PondSfcThinSnwComb < 0.0 ) then                ! IF SnowIce AND SnowLiqWater SUBLIMATES REMOVE FROM SOIL
                      SoilIce(1) = SoilIce(1) + PondSfcThinSnwComb/(ThicknessSnowSoilLayer(1)*1000.0) ! negative SoilIce from oversublimation is adjusted below
                      PondSfcThinSnwComb = 0.0
                   endif
                   SnowWaterEquiv = 0.0
                   SnowDepth      = 0.0
                endif ! if(SnowIce(J) >= 0.0)
                SnowLiqWater(J)   = 0.0
                SnowIce(J)        = 0.0
                ThicknessSnowSoilLayer(J) = 0.0

                ! SNICAR, aerosol flux may infiltrate into top soil like PondSfcThinSnwComb, it
                ! would be more thorough to do so later  
                if ( OptSnowAlbedo == 3 ) then
                   MassBChydropho(J) = 0.0 
                   MassBChydrophi(J) = 0.0       
                   MassOChydropho(J) = 0.0    
                   MassOChydrophi(J) = 0.0    
                   MassDust1(J)      = 0.0        
                   MassDust2(J)      = 0.0        
                   MassDust3(J)      = 0.0     
                   MassDust4(J)      = 0.0  
                   MassDust5(J)      = 0.0
                endif

             endif ! if(NumSnowLayerNeg < -1)
          endif ! if(J /= 0)

          ! shift all elements above this down by one.
          if ( (J > NumSnowLayerNeg+1) .and. (NumSnowLayerNeg < -1) ) then
             do I = J, NumSnowLayerNeg+2, -1
                TemperatureSoilSnow(I)    = TemperatureSoilSnow(I-1)
                SnowLiqWater(I)           = SnowLiqWater(I-1)
                SnowIce(I)                = SnowIce(I-1)
                ThicknessSnowSoilLayer(I) = ThicknessSnowSoilLayer(I-1)

                if ( OptSnowAlbedo == 3 ) then
                   MassBChydropho(I)      = MassBChydropho(I-1)
                   MassBChydrophi(I)      = MassBChydrophi(I-1)
                   MassOChydropho(I)      = MassOChydropho(I-1)
                   MassOChydrophi(I)      = MassOChydrophi(I-1)
                   MassDust1(I)           = MassDust1(I-1)
                   MassDust2(I)           = MassDust2(I-1)
                   MassDust3(I)           = MassDust3(I-1)
                   MassDust4(I)           = MassDust4(I-1)
                   MassDust5(I)           = MassDust5(I-1)
                   SnowRadius(I)          = SnowRadius(I-1)
                endif

             enddo
          endif
          NumSnowLayerNeg = NumSnowLayerNeg + 1

       endif ! if(SnowIce(J) <= 0.1)
    enddo ! do J

    ! to conserve water in case of too large surface sublimation
    if ( SoilIce(1) < 0.0) then
       SoilLiqWater(1) = SoilLiqWater(1) + SoilIce(1)
       SoilIce(1)      = 0.0
    endif

    if ( NumSnowLayerNeg ==0 ) return   ! MB: get out if no longer multi-layer

    SnowWaterEquiv = 0.0
    SnowDepth      = 0.0
    SnowIceTmp     = 0.0
    SnowLiqTmp     = 0.0

    do J = NumSnowLayerNeg+1, 0
       SnowWaterEquiv = SnowWaterEquiv + SnowIce(J) + SnowLiqWater(J)
       SnowDepth      = SnowDepth + ThicknessSnowSoilLayer(J)
       SnowIceTmp     = SnowIceTmp + SnowIce(J)
       SnowLiqTmp     = SnowLiqTmp + SnowLiqWater(J)
    enddo

    ! check the snow depth - all snow gone, the liquid water assumes ponding on soil surface.
    ! if ( (SnowDepth < 0.05) .and. (NumSnowLayerNeg < 0) ) then
    if ( (SnowDepth < 0.025) .and. (NumSnowLayerNeg < 0) ) then ! MB: change limit
       NumSnowLayerNeg     = 0
       SnowWaterEquiv      = SnowIceTmp
       PondSfcThinSnwTrans = SnowLiqTmp                ! LIMIT OF NumSnowLayerNeg < 0 MEANS INPUT PONDING
       if ( SnowWaterEquiv <= 0.0 ) SnowDepth = 0.0    ! SHOULD BE ZERO; SEE ABOVE
    endif

    ! check the snow depth - snow layers combined
    if ( NumSnowLayerNeg < -1 ) then
       NumSnowLayerOld = NumSnowLayerNeg
       IndLayer        = 1
       do I = NumSnowLayerOld+1, 0
          if ( ThicknessSnowSoilLayer(I) < SnowThickMin(IndLayer) ) then
             if ( I == NumSnowLayerNeg+1 ) then
                IndNeighbor = I + 1
             else if ( I == 0 ) then
                IndNeighbor = I - 1
             else
                IndNeighbor = I + 1
                if ( (ThicknessSnowSoilLayer(I-1)+ThicknessSnowSoilLayer(I)) < &
                     (ThicknessSnowSoilLayer(I+1)+ThicknessSnowSoilLayer(I)) ) IndNeighbor = I-1
             endif
             ! Node l and j are combined and stored as node j.
             if ( IndNeighbor > I ) then
                J = IndNeighbor
                L = I
             else
                J = I
                L = IndNeighbor
             endif

             if ( OptSnowAlbedo == 3 ) then
                MassBChydropho(J) = MassBChydropho(J) +  MassBChydropho(L)
                MassBChydrophi(J) = MassBChydrophi(J) +  MassBChydrophi(L)
                MassOChydropho(J) = MassOChydropho(J) +  MassOChydropho(L)
                MassOChydrophi(J) = MassOChydrophi(J) +  MassOChydrophi(L)
                MassDust1(J)      = MassDust1(J) +  MassDust1(L)
                MassDust2(J)      = MassDust2(J) +  MassDust2(L)
                MassDust3(J)      = MassDust3(J) +  MassDust3(L)
                MassDust4(J)      = MassDust4(J) +  MassDust4(L)
                MassDust5(J)      = MassDust5(J) +  MassDust5(L)
                SnowRadius(J)     = (SnowRadius(J)*(SnowLiqWater(J)+SnowIce(J)) + SnowRadius(L)*(SnowLiqWater(L)+SnowIce(L))) / &
                                    (SnowLiqWater(J) + SnowIce(J) + SnowLiqWater(L) + SnowIce(L))
             endif

             ! update combined snow water & temperature
             call SnowLayerWaterCombo(ThicknessSnowSoilLayer(J), SnowLiqWater(J), SnowIce(J), TemperatureSoilSnow(J), &
                                      ThicknessSnowSoilLayer(L), SnowLiqWater(L), SnowIce(L), TemperatureSoilSnow(L) )


             ! Now shift all elements above this down one.
             if ( (J-1) > (NumSnowLayerNeg+1) ) then
                do K = J-1, NumSnowLayerNeg+2, -1
                   TemperatureSoilSnow(K)    = TemperatureSoilSnow(K-1)
                   SnowIce(K)                = SnowIce(K-1)
                   SnowLiqWater(K)           = SnowLiqWater(K-1)
                   ThicknessSnowSoilLayer(K) = ThicknessSnowSoilLayer(K-1)

                   if ( OptSnowAlbedo == 3 ) then
                      MassBChydropho(K)      = MassBChydropho(K-1) 
                      MassBChydrophi(K)      = MassBChydrophi(K-1) 
                      MassOChydropho(K)      = MassOChydropho(K-1) 
                      MassOChydrophi(K)      = MassOChydrophi(K-1) 
                      MassDust1(K)           = MassDust1(K-1) 
                      MassDust2(K)           = MassDust2(K-1) 
                      MassDust3(K)           = MassDust3(K-1) 
                      MassDust4(K)           = MassDust4(K-1) 
                      MassDust5(K)           = MassDust5(K-1)
                      SnowRadius(K)          = SnowRadius(K-1)
                   endif

                enddo
             endif
             ! Decrease the number of snow layers
             NumSnowLayerNeg = NumSnowLayerNeg + 1
             if ( NumSnowLayerNeg >= -1 ) Exit
          else 
             ! The layer thickness is greater than the prescribed minimum value
             IndLayer = IndLayer + 1
          endif
       enddo
    endif

    end associate

  end subroutine SnowLayerCombine

end module SnowLayerCombineMod
