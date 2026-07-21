module SnowCoverGroundAR25Mod

!!! Compute ground snow cover fraction based on Abolafia-Rosenzweig et al. (2025, JAMES) scheme

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowCoverGroundAR25(noahmp)

! ------------------------ Code history -----------------------------------
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Code updated: R. Abolafia-Rosenzweig and C. He (Abolafia-Rosenzweig et al., 2025)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    real(kind=kind_noahmp)           :: SnowDensBulk   ! bulk density of snow [Kg/m3]
    real(kind=kind_noahmp)           :: MeltFac        ! melting factor for snow cover frac
    real(kind=kind_noahmp)           :: SnowMeltFac    ! snowmelt m parameter (scale-dependent)
    real(kind=kind_noahmp)           :: SnowCoverFac   ! snow cover factor [scfac] (scale-dependent)
! --------------------------------------------------------------------
    associate(                                                           &
              SnowDepth         => noahmp%water%state%SnowDepth         ,& ! in,  snow depth [m]
              SnowWaterEquiv    => noahmp%water%state%SnowWaterEquiv    ,& ! in,  snow water equivalent [mm]
              GridSize          => noahmp%config%domain%GridSize        ,& ! in,  noahmp model grid spacing [m]
              SnowCoverM1AR25   => noahmp%water%param%SnowCoverM1AR25   ,& ! in,  SCFm1 parameter from AR2025
              SnowCoverM2AR25   => noahmp%water%param%SnowCoverM2AR25   ,& ! in,  SCFm2 parameter from AR2025
              SnowCoverFac1AR25 => noahmp%water%param%SnowCoverFac1AR25 ,& ! in,  SCfac1 parameter from AR2025
              SnowCoverFac2AR25 => noahmp%water%param%SnowCoverFac2AR25 ,& ! in,  SCfac2 parameter from AR2025
              SnowCoverFrac     => noahmp%water%state%SnowCoverFrac      & ! out, snow cover fraction
             )
! ----------------------------------------------------------------------

    SnowCoverFrac = 0.0
    if ( SnowDepth > 0.0 ) then
         ! calculate SCF parameters as a function of grid size (limit gridsize to 500m~36km due to parameterization limitation)
         SnowMeltFac  = SnowCoverM1AR25 + tanh(SnowCoverM2AR25 * min((max(GridSize,500.0)/1000),36.0))
         SnowCoverFac = SnowCoverFac1AR25 * sinh(SnowCoverFac2AR25*min((max(GridSize,500.0)/1000),36.0)) + SnowCoverFac2AR25
        
         ! using scale-dependent parameters, employ the Niu-Yang 07 SCF soluiton
         SnowDensBulk  = SnowWaterEquiv / SnowDepth
         MeltFac       = (SnowDensBulk / 100.0)**SnowMeltFac
         SnowCoverFrac = tanh( SnowDepth /(SnowCoverFac * MeltFac))
    endif

    end associate

  end subroutine SnowCoverGroundAR25

end module SnowCoverGroundAR25Mod
