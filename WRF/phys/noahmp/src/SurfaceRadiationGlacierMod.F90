module SurfaceRadiationGlacierMod

!!! Compute glacier surface radiative fluxes (absorption and reflection)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SurfaceRadiationGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: RADIATION_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                                             :: IndBand                 ! waveband indices (1=vis, 2=nir)
    integer                                             :: IndLoop                 ! snow and soil layer loop
    real(kind=kind_noahmp)                              :: RadSwAbsGrdTmp          ! ground absorbed solar radiation [W/m2]
    real(kind=kind_noahmp)                              :: RadSwReflGrdTmp         ! ground reflected solar radiation [W/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FracRadSwAbsSnowDirMean ! direct solar flux factor absorbed by snow [frc] scaling
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: FracRadSwAbsSnowDifMean ! diffuse solar flux factor absorbed by snow [frc] scaling

! -----------------------------------------------------------------
    associate(                                                                  &
              NumSwRadBand          => noahmp%config%domain%NumSwRadBand       ,& ! in,  number of solar radiation wave bands
              NumSnowLayerMax       => noahmp%config%domain%NumSnowLayerMax    ,& ! in,  maximum number of snow layers
              NumSnowLayerNeg       => noahmp%config%domain%NumSnowLayerNeg    ,& ! in,  actual number of snow layers (negative)
              OptSnowAlbedo         => noahmp%config%nmlist%OptSnowAlbedo      ,& ! in,  options for ground snow surface albedo
              RadSwDownDir          => noahmp%energy%flux%RadSwDownDir         ,& ! in,  incoming direct solar radiation [W/m2]
              RadSwDownDif          => noahmp%energy%flux%RadSwDownDif         ,& ! in,  incoming diffuse solar radiation [W/m2]
              FracRadSwAbsSnowDir   => noahmp%energy%flux%FracRadSwAbsSnowDir  ,& ! in,  direct solar flux factor absorbed by snow [frc]
              FracRadSwAbsSnowDif   => noahmp%energy%flux%FracRadSwAbsSnowDif  ,& ! in,  diffuse solar flux factor absorbed by snow [frc]
              AlbedoGrdDir          => noahmp%energy%state%AlbedoGrdDir        ,& ! in,  ground albedo (direct beam: vis, nir)
              AlbedoGrdDif          => noahmp%energy%state%AlbedoGrdDif        ,& ! in,  ground albedo (diffuse: vis, nir)
              AlbedoSnowDir         => noahmp%energy%state%AlbedoSnowDir       ,& ! in,  snow albedo for direct(1=vis, 2=nir)
              AlbedoSnowDif         => noahmp%energy%state%AlbedoSnowDif       ,& ! in,  snow albedo for diffuse(1=vis, 2=nir)
              SnowCoverFrac         => noahmp%water%state%SnowCoverFrac        ,& ! in,  snow cover fraction
              AlbedoLandIce         => noahmp%energy%param%AlbedoLandIce       ,& ! in,  albedo land ice: 1=vis, 2=nir
              RadSwAbsGrd           => noahmp%energy%flux%RadSwAbsGrd          ,& ! out, solar radiation absorbed by ground [W/m2]
              RadSwAbsSfc           => noahmp%energy%flux%RadSwAbsSfc          ,& ! out, total absorbed solar radiation [W/m2]
              RadSwReflSfc          => noahmp%energy%flux%RadSwReflSfc         ,& ! out, total reflected solar radiation [W/m2]
              RadSwAbsSnowSoilLayer => noahmp%energy%flux%RadSwAbsSnowSoilLayer & ! out, total absorbed solar radiation by snow for each layer [W/m2]
             )
! ----------------------------------------------------------------------

    ! initialization
    RadSwAbsGrd  = 0.0
    RadSwAbsSfc  = 0.0
    RadSwReflSfc = 0.0
    if ( OptSnowAlbedo == 3 ) then
       if (.not. allocated(FracRadSwAbsSnowDirMean)) allocate(FracRadSwAbsSnowDirMean(-NumSnowLayerMax+1:1,1:NumSwRadBand))
       if (.not. allocated(FracRadSwAbsSnowDifMean)) allocate(FracRadSwAbsSnowDifMean(-NumSnowLayerMax+1:1,1:NumSwRadBand))
       FracRadSwAbsSnowDirMean(:,:) = 0.0
       FracRadSwAbsSnowDifMean(:,:) = 0.0
       RadSwAbsSnowSoilLayer(:)     = 0.0
    endif

    ! aggretate radiative flux
    do IndBand = 1, NumSwRadBand
       ! solar radiation absorbed by glacier surface
       RadSwAbsGrdTmp  = RadSwDownDir(IndBand) * (1.0 - AlbedoGrdDir(IndBand)) + &
                         RadSwDownDif(IndBand) * (1.0 - AlbedoGrdDif(IndBand))
       RadSwAbsGrd     = RadSwAbsGrd + RadSwAbsGrdTmp
       RadSwAbsSfc     = RadSwAbsSfc + RadSwAbsGrdTmp
      
       ! solar radiation reflected by glacier surface
       RadSwReflGrdTmp = RadSwDownDir(IndBand) * AlbedoGrdDir(IndBand) + &
                         RadSwDownDif(IndBand) * AlbedoGrdDif(IndBand)
       RadSwReflSfc    = RadSwReflSfc + RadSwReflGrdTmp

       if ( OptSnowAlbedo == 3 ) then
          do IndLoop = -NumSnowLayerMax+1, 1, 1
             FracRadSwAbsSnowDirMean(IndLoop,IndBand) = FracRadSwAbsSnowDir(IndLoop,IndBand) * SnowCoverFrac +    &
                                                        ((1.0 - SnowCoverFrac) * (1.0 - AlbedoLandIce(IndBand)) * &
                                                        (FracRadSwAbsSnowDir(IndLoop,IndBand)/(1.0-AlbedoSnowDir(IndBand))))
             FracRadSwAbsSnowDifMean(IndLoop,IndBand) = FracRadSwAbsSnowDif(IndLoop,IndBand) * SnowCoverFrac +    &
                                                        ((1.0 - SnowCoverFrac) * (1.0 - AlbedoLandIce(IndBand)) * &
                                                        (FracRadSwAbsSnowDif(IndLoop,IndBand)/(1.0-AlbedoSnowDif(IndBand))))
             RadSwAbsSnowSoilLayer(IndLoop) = RadSwAbsSnowSoilLayer(IndLoop) + &
                                              RadSwDownDir(IndBand) * FracRadSwAbsSnowDirMean(IndLoop,IndBand) +  &
                                              RadSwDownDif(IndBand) * FracRadSwAbsSnowDifMean(IndLoop,IndBand)
          enddo
       endif
    enddo

    if (OptSnowAlbedo == 3 .and. NumSnowLayerNeg == 0) then
       RadSwAbsSnowSoilLayer(:) = 0.0
       RadSwAbsSnowSoilLayer(1) = RadSwAbsGrd
    endif

    if ( OptSnowAlbedo == 3 ) then
       deallocate(FracRadSwAbsSnowDirMean)
       deallocate(FracRadSwAbsSnowDifMean)
    endif

    end associate

  end subroutine SurfaceRadiationGlacier

end module SurfaceRadiationGlacierMod
