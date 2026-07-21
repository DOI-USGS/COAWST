module SnowAlbedoSnicarMod

!!! Compute snow albedo based on SNICAR scheme (Flanner et al. (2021) GMD)

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use SnowRadiationSnicarMod, only : SnowRadiationSnicar

  implicit none

contains

  subroutine SnowAlbedoSnicar(noahmp)

! ------------------------ Code history -----------------------------------
! Implementation: T.-S. Lin, C. He, et al. (2025, JHM)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: FlagSwRadType  ! flag: 1 for direct-beam incident flux, 2 for diffuse incident flux

! --------------------------------------------------------------------
    associate(                                                          &
              NumSwRadBand  => noahmp%config%domain%NumSwRadBand ,& ! in,  number of solar radiation wave bands
              AlbedoSnowDir => noahmp%energy%state%AlbedoSnowDir ,& ! out, snow albedo for direct (1=vis, 2=nir)
              AlbedoSnowDif => noahmp%energy%state%AlbedoSnowDif  & ! out, snow albedo for diffuse (1=vis, 2=nir)
             )
! ----------------------------------------------------------------------

    ! initialization
    AlbedoSnowDir(1:NumSwRadBand) = 0.0
    AlbedoSnowDif(1:NumSwRadBand) = 0.0

    FlagSwRadType = 1 ! Direct
    call SnowRadiationSnicar(noahmp,FlagSwRadType) 

    FlagSwRadType = 2 ! Diffuse
    call SnowRadiationSnicar(noahmp,FlagSwRadType)
    
    end associate

  end subroutine SnowAlbedoSnicar

end module SnowAlbedoSnicarMod
