module SnowAgingSnicarMod

!!! Compute snow effective grain size (radius) based on SNICAR scheme (Flanner et al. (2021) GMD)
!!! Description: SNICAR snow aging process, contributions to grain size evolution are from:
!!! 1. vapor redistribution (dry snow) 
!!! 2. liquid water redistribution (wet snow)
!!! 3. re-freezing of liquid water
!!! Vapor redistribution: Method is to retrieve 3 best-bit parameters that
!!! depend on snow temperature, temperature gradient, and density,
!!! that are derived from the microphysical model described in: 
!!! Flanner and Zender (2006), Linking snowpack microphysics and albedo
!!! evolution, J. Geophys. Res., 111, D12208, doi:10.1029/2005JD006834. 
!!! The parametric equation has the form: 
!!! dr/dt = drdt_0*(tau/(dr_fresh+tau))^(1/kappa), where: r is the effective radius,
!!! tau and kappa are best-fit parameters,
!!! drdt_0 is the initial rate of change of effective radius, and
!!! dr_fresh is the difference between the current and fresh snow states (r_current - r_fresh).
!!! Liquid water redistribution: Apply the grain growth function from:
!!! Brun, E. (1989), Investigation of wet-snow metamorphism in respect of 
!!! liquid-water content, Annals of Glaciology, 13, 22-26.
!!! There are two parameters that describe the grain growth rate as 
!!! a function of snow liquid water content (LWC). The "LWC=0" parameter
!!! is zeroed here because we are accounting for dry snowing with a different representation
!!! Re-freezing of liquid water: Assume that re-frozen liquid water clumps
!!! into an arbitrarily large effective grain size (SnowRadiusRefrz). 
!!! The phenomenon is observed (Grenfell), but so far not well quantified.

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowAgingSnicar(noahmp)

! ------------------------ Code history -----------------------------------
! Implementation: T.-S. Lin, C. He, et al. (2025, JHM)
! Adapted from Flanner and Zender (2006) in CTSM, SnowAge_grain module
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: SnowLayerTop                         ! top snow layer index [idx]
    integer                          :: SnowLayerBottom                      ! bottom snow layer index [idx]
    integer                          :: TemperatureInd                       ! snow aging lookup table temperature index [idx]
    integer                          :: TemperatureGradientInd               ! snow aging lookup table temperature gradient index [idx]
    integer                          :: SnowDensityInd                       ! snow aging lookup table snow density index [idx]
    integer                          :: LoopInd                              ! do loop/array indices
    integer, parameter               :: IndTempSnwAgeMin     = 1             ! minimum temperature index used in aging lookup table [idx]
    integer, parameter               :: IndTempGradSnwAgeMin = 1             ! minimum temperature gradient index used in aging lookup table [idx]
    integer, parameter               :: IndDensitySnwAgeMin  = 1             ! minimum snow density index used in aging lookup table [idx]
    real(kind=kind_noahmp)           :: bst_tau                              ! best fit snow aging parameter retrieved from lookup table [hour]
    real(kind=kind_noahmp)           :: bst_kappa                            ! best fit snow aging parameter retrieved from lookup table [unitless]
    real(kind=kind_noahmp)           :: bst_drdt0                            ! best fit snow aging parameter retrieved from lookup table [um hr-1]
    real(kind=kind_noahmp)           :: SnowMassLayer                        ! liquid + solid H2O in snow layer [kg m-2]
    real(kind=kind_noahmp)           :: TemperatureSnowLayerTop              ! temperature at upper layer boundary [K]
    real(kind=kind_noahmp)           :: TemperatureSnowLayerBottom           ! temperature at lower layer boundary [K]
    real(kind=kind_noahmp)           :: SnowDensity                          ! snow density [kg m-3]
    real(kind=kind_noahmp)           :: SnowRadiusChgTot                     ! incremental change in snow effective radius [um]
    real(kind=kind_noahmp)           :: SnowRadiusChgWet                     ! incremental change in snow effective radius from wet growth [um]
    real(kind=kind_noahmp)           :: SnowRadiusChgFresh                   ! difference between fresh snow r_e and current r_e [um]
    real(kind=kind_noahmp)           :: NewSnow                              ! fresh snowfall [kg m-2]
    real(kind=kind_noahmp)           :: RefrzSnow                            ! re-frozen snow [kg m-2]
    real(kind=kind_noahmp)           :: FracRefrz                            ! fraction of layer mass that is re-frozen snow [frc]
    real(kind=kind_noahmp)           :: FracNewSnow                          ! fraction of layer mass that is new snow [frc]
    real(kind=kind_noahmp)           :: FracOldSnow                          ! fraction of layer mass that is old snow [frc]
    real(kind=kind_noahmp)           :: FracLiqWater                         ! fraction of layer mass that is liquid water[frc]    
    real(kind=kind_noahmp), allocatable, dimension(:) :: TemperatureGradient ! snow temperature gradient (lyr) [K m-1]

! --------------------------------------------------------------------
    associate(                                                                         &
              MainTimeStep            => noahmp%config%domain%MainTimeStep            ,& ! in,  noahmp main time step [s]
              NumSnowLayerMax         => noahmp%config%domain%NumSnowLayerMax         ,& ! in,  maximum number of snow layers
              NumSnowLayerNeg         => noahmp%config%domain%NumSnowLayerNeg         ,& ! in,  actual number of snow layers (negative)
              ThicknessSnowSoilLayer  => noahmp%config%domain%ThicknessSnowSoilLayer  ,& ! in,  thickness of snow/soil layers [m]
              NumTempSnwAgeSnicar     => noahmp%config%domain%NumTempSnwAgeSnicar     ,& ! in,  maxiumum temperature index used in aging lookup table [idx]
              NumTempGradSnwAgeSnicar => noahmp%config%domain%NumTempGradSnwAgeSnicar ,& ! in,  maxiumum temperature gradient index used in aging lookup table [idx]   
              NumDensitySnwAgeSnicar  => noahmp%config%domain%NumDensitySnwAgeSnicar  ,& ! in,  maxiumum snow density index used in aging lookup table [idx]
              TemperatureSoilSnow     => noahmp%energy%state%TemperatureSoilSnow      ,& ! in,  snow and soil layer temperature [K] 
              SnowIce                 => noahmp%water%state%SnowIce                   ,& ! in,  snow layer ice [mm]
              SnowLiqWater            => noahmp%water%state%SnowLiqWater              ,& ! in,  snow layer liquid water [mm]
              SnowRadiusFresh         => noahmp%water%state%SnowRadiusFresh           ,& ! in,  fresh snow radius [microns]
              SnowWaterEquiv          => noahmp%water%state%SnowWaterEquiv            ,& ! in,  snow water equivalent [mm]
              SnowDepth               => noahmp%water%state%SnowDepth                 ,& ! in,  snow depth [m]
              SnowfallGround          => noahmp%water%flux%SnowfallGround             ,& ! in,  snowfall at ground surface [mm/s]
              SnowFreezeRate          => noahmp%water%flux%SnowFreezeRate             ,& ! in,  rate of snow freezing [mm/s]
              SnowRadiusMin           => noahmp%water%param%SnowRadiusMin             ,& ! in,  minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
              SnowRadiusMax           => noahmp%water%param%SnowRadiusMax             ,& ! in,  maximum allowed snow effective radius [microns]
              SnowWetAgeC1Brun89      => noahmp%water%param%SnowWetAgeC1Brun89        ,& ! in,  constant for liquid water grain growth [m3 s-1], from Brun89
              SnowWetAgeC2Brun89      => noahmp%water%param%SnowWetAgeC2Brun89        ,& ! in,  constant for liquid water grain growth [m3 s-1], from Brun89 corrected for LWC
              SnowAgeScaleFac         => noahmp%water%param%SnowAgeScaleFac           ,& ! in,  arbitrary tuning/scaling factor applied to snow aging rate (-)
              SnowRadiusRefrz         => noahmp%water%param%SnowRadiusRefrz           ,& ! in,  effective radius of re-frozen snow [microns]
              snowage_tau             => noahmp%water%param%snowage_tau               ,& ! in,  snowage tau from table [hours]
              snowage_kappa           => noahmp%water%param%snowage_kappa             ,& ! in,  snowage kappa from table [unitless]
              snowage_drdt0           => noahmp%water%param%snowage_drdt0             ,& ! in,  snowage dr/dt_0 from table [m2 kg-1 hr-1]
              SnowRadius              => noahmp%water%state%SnowRadius                 & ! out, effective grain radius [microns, m-6]
             )
! ----------------------------------------------------------------------

    ! initialize
    if (.not. allocated(TemperatureGradient)) allocate(TemperatureGradient(-NumSnowLayerMax:0))
    SnowLayerBottom = 0
    SnowLayerTop = NumSnowLayerNeg + 1

    ! loop over snow layers
    do LoopInd = SnowLayerTop, SnowLayerBottom, 1
     
      !**********  1. DRY SNOW AGING  ***********
      SnowMassLayer = SnowLiqWater(LoopInd) + SnowIce(LoopInd)
 
      ! temperature gradient
      if (LoopInd == SnowLayerTop) then ! top layer
         TemperatureSnowLayerTop    = TemperatureSoilSnow(SnowLayerTop)
         TemperatureSnowLayerBottom = (TemperatureSoilSnow(LoopInd+1) * ThicknessSnowSoilLayer(LoopInd) &
                                      + TemperatureSoilSnow(LoopInd) * ThicknessSnowSoilLayer(LoopInd+1)) &
                                      / (ThicknessSnowSoilLayer(LoopInd) + ThicknessSnowSoilLayer(LoopInd+1))
      else
         TemperatureSnowLayerTop    = (TemperatureSoilSnow(LoopInd-1) * ThicknessSnowSoilLayer(LoopInd) &
                                      + TemperatureSoilSnow(LoopInd) * ThicknessSnowSoilLayer(LoopInd-1)) &
                                      / (ThicknessSnowSoilLayer(LoopInd) + ThicknessSnowSoilLayer(LoopInd-1))
         TemperatureSnowLayerBottom = (TemperatureSoilSnow(LoopInd+1) * ThicknessSnowSoilLayer(LoopInd) &
                                      + TemperatureSoilSnow(LoopInd) * ThicknessSnowSoilLayer(LoopInd+1)) &
                                      / (ThicknessSnowSoilLayer(LoopInd) + ThicknessSnowSoilLayer(LoopInd+1))
      endif
      TemperatureGradient(LoopInd) = abs((TemperatureSnowLayerTop - TemperatureSnowLayerBottom) / &
                                         ThicknessSnowSoilLayer(LoopInd))

      ! snow density
      SnowDensity = SnowMassLayer / ThicknessSnowSoilLayer(LoopInd)

      ! make sure snow density doesn't drop below 50 (see SnowDensityInd below)
      SnowDensity = max(50.0, SnowDensity)

      ! best-fit table indices
      TemperatureInd         = nint((TemperatureSoilSnow(LoopInd)-223.15) / 5) + 1
      TemperatureGradientInd = nint(TemperatureGradient(LoopInd) / 10) + 1
      SnowDensityInd         = nint((SnowDensity-50.0) / 50.0) + 1

      ! boundary check:
      if (TemperatureInd < IndTempSnwAgeMin) then
         TemperatureInd = IndTempSnwAgeMin
      endif
      if (TemperatureInd > NumTempSnwAgeSnicar) then
         TemperatureInd = NumTempSnwAgeSnicar
      endif
      if (TemperatureGradientInd < IndTempGradSnwAgeMin) then
         TemperatureGradientInd = IndTempGradSnwAgeMin
      endif
      if (TemperatureGradientInd > NumTempGradSnwAgeSnicar) then
         TemperatureGradientInd = NumTempGradSnwAgeSnicar
      endif
      if (SnowDensityInd < IndDensitySnwAgeMin) then
         SnowDensityInd = IndDensitySnwAgeMin
      endif
      if (SnowDensityInd > NumDensitySnwAgeSnicar) then
         SnowDensityInd = NumDensitySnwAgeSnicar
      endif

      ! best-fit parameters
      bst_tau   = snowage_tau(SnowDensityInd,TemperatureGradientInd,TemperatureInd)
      bst_kappa = snowage_kappa(SnowDensityInd,TemperatureGradientInd,TemperatureInd)
      bst_drdt0 = snowage_drdt0(SnowDensityInd,TemperatureGradientInd,TemperatureInd)

      ! extra boundary check, to prevent when using old restart file with lower SnowRadiusMin than current run
      if (SnowRadius(LoopInd) < SnowRadiusMin) then
         SnowRadius(LoopInd) = SnowRadiusMin
      endif
      if (SnowRadius(LoopInd) < SnowRadiusFresh) then
         SnowRadius(LoopInd) = SnowRadiusFresh
      endif

      ! change in snow effective radius, using best-fit parameters
      SnowRadiusChgFresh = SnowRadius(LoopInd) - SnowRadiusMin
      SnowRadiusChgTot   = (bst_drdt0 * (bst_tau/(SnowRadiusChgFresh+bst_tau))**(1.0/bst_kappa)) * &
                           (MainTimeStep/3600.0)


      !**********  2. WET SNOW AGING  ***********
      ! We are assuming wet and dry evolution occur simultaneously, and 
      ! the contributions from both can be summed. 
      ! This is justified by setting the linear offset constant SnowWetAgeC1Brun89 to zero [Brun, 1989]

      ! liquid water faction
      FracLiqWater = min(0.1, (SnowLiqWater(LoopInd)/SnowMassLayer))
      SnowRadiusChgWet = 1.0e18 * ( MainTimeStep*(SnowWetAgeC1Brun89 + SnowWetAgeC2Brun89*(FracLiqWater**(3))) / &
                         (4.0 * ConstPI * SnowRadius(LoopInd)**(2)) )
      SnowRadiusChgTot = SnowRadiusChgTot + SnowRadiusChgWet


      !**********  3. SNOWAGE SCALING (TUNING OPTION)              ***********
      ! Multiply rate of change of effective radius by some constant, SnowAgeScaleFac

      SnowRadiusChgTot = SnowRadiusChgTot * SnowAgeScaleFac


      !**********  4. INCREMENT EFFECTIVE RADIUS, ACCOUNTING FOR:  ***********
      ! DRY AGING, WET AGING, FRESH SNOW, RE-FREEZING

      ! new snowfall [kg/m2]
      NewSnow = max(0.0, (SnowfallGround*MainTimeStep))

      ! snow that has re-frozen [kg/m2]
      RefrzSnow = max(0.0, (SnowFreezeRate(LoopInd)*MainTimeStep)) 

      ! fraction of layer mass that is re-frozen
      FracRefrz = RefrzSnow / SnowMassLayer

      ! fraction of layer mass that is new snow
      if (LoopInd == SnowLayerTop) then
         FracNewSnow = NewSnow / SnowMassLayer
      else
         FracNewSnow = 0.0
      endif

      if ((FracRefrz + FracNewSnow) > 1.0) then
         FracRefrz = FracRefrz / (FracRefrz + FracNewSnow)
         FracNewSnow = 1.0 - FracRefrz
         FracOldSnow = 0.0
      else
         FracOldSnow = 1.0 - FracRefrz - FracNewSnow
      endif

      ! mass-weighted mean of fresh snow, old snow, and re-frozen snow effective radius
      SnowRadius(LoopInd) = (SnowRadius(LoopInd) + SnowRadiusChgTot)*FracOldSnow + &
                            SnowRadiusFresh*FracNewSnow + SnowRadiusRefrz*FracRefrz


      !**********  5. CHECK BOUNDARIES   ***********
      ! boundary check

      if (SnowRadius(LoopInd) < SnowRadiusMin) then
         SnowRadius(LoopInd) = SnowRadiusMin
      endif

      if (SnowRadius(LoopInd) > SnowRadiusMax) then
         SnowRadius(LoopInd) = SnowRadiusMax
      end if

    enddo ! layer loop

    ! sanity check for snow layer
    if (-NumSnowLayerMax /= NumSnowLayerNeg) then
       SnowRadius(-NumSnowLayerMax:NumSnowLayerNeg) = 0.0
    endif

    if (NumSnowLayerNeg == 0) then
       SnowRadius(:) = 0.0
    endif

    ! special case: snow on ground, but not enough to have defined a snow layer:
    ! set SnowRadius to fresh snow grain size
    if (NumSnowLayerNeg == 0 .and. &
        ((SnowfallGround > 0.0) .or. (SnowWaterEquiv > 0.0) .or. (SnowDepth > 0.0))) then
       SnowRadius(SnowLayerBottom) = SnowRadiusFresh  !SnowRadiusMin
    endif

    deallocate(TemperatureGradient)

    end associate

  end subroutine SnowAgingSnicar

end module SnowAgingSnicarMod
