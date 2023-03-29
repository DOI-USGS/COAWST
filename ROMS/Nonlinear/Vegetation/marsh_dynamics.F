#include "cppdefs.h"
 
       MODULE marsh_dynamics_mod
#if defined VEGETATION && defined MARSH_DYNAMICS
!
!svn $Id: marsh_dynamics.F 429 2015-04-20 17:30:26Z arango $
!======================================================================!
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==============================================Tarandeep S. Kalra======!
!================================================Neil K. Ganju  =======!
!================================================John C. Warner========!
!==============================================Julia M. Moriarty=======!
!===========================================Alfredo Aretxabaleta=======!
!                                                                      ! 
!  This routine is the main driver for marsh dynamics and accounts     !
!  for two processes through various sub-routines:                     !
!======================================================================!
!   1. Calculate lateral wave thrust and sediment release from         !
!      edge erosion                                                    !
!   2. Calculate vertical biomass production (autochthonous supply)    !
!      (organic sediment)                                              !
!======================================================================!
!                                                                      !
!  This is done through the following modules:                         !
!                                                                      !
!  1. Compute wave thrust on marshes                                   !  
!     --> marsh_wave_thrust.F                                          !
!  2. Obtain sediment release (i.e. lateral erosion) out of marsh cells!
!     --> marsh_sed_erosion.F                                          !
!  3. Compute the mean tidal range and marsh high water                !
!     --> marsh_tidal_range.F                                          !
!  4. Compute marsh biomass production and vertical accretion          !
!     --> marsh_vert_growth.F                                          !
!                                                                      !
!  References:                                                         !   
!  Kalra, T.S., Ganju, N.K., Aretxabaleta, A., Moriarty, J.M.,         !
!  Aretxabaleta, A., Warner, J.C., Carr, J., Leonardi., N.,            !
!  Modeling Marsh Dynamics Using a 3-D Coupled Wave-Flow-Sediment      !
!  Model (In process).                                                 !
!                                                                      !
!======================================================================! 
!                                                                      !
      implicit none

      PRIVATE
      PUBLIC  :: marsh_dynamics

      CONTAINS
!
!***********************************************************************
!
      SUBROUTINE marsh_dynamics (ng, tile)
!
!***********************************************************************
!
      USE mod_stepping 
      USE mod_scalars
# if defined MARSH_WAVE_THRUST
      USE marsh_wave_thrust_mod, ONLY : marsh_wave_thrust
# endif 
# if defined MARSH_SED_EROSION
      USE marsh_sed_erosion_mod, ONLY : marsh_sed_erosion 
# endif 
# if defined MARSH_TIDAL_RANGE_CALC
      USE marsh_tidal_range_calc_mod, ONLY : marsh_tidal_range_calc
# endif 
# if defined MARSH_VERT_GROWTH 
      USE marsh_vert_growth_mod, ONLY : marsh_vert_growth
# endif 
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!

# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 16)
# endif
!
!***********************************************************************
! Wave induced thrust calculation on marsh boundary.
!***********************************************************************
!
# if defined MARSH_WAVE_THRUST
      CALL marsh_wave_thrust(ng, tile)
# endif 
!
# if defined MARSH_SED_EROSION
!
!***********************************************************************
! Marsh sediment erosion based on bathy change in marsh cells.
!***********************************************************************
! 
      CALL marsh_sed_erosion(ng, tile) 
# endif 
!
# if defined MARSH_TIDAL_RANGE_CALC
!
!***********************************************************************
! Compute tidal range and mean high water over a given frequency. 
!***********************************************************************
!
      CALL marsh_tidal_range_calc(ng, tile) 
#  endif 
!

#  if defined MARSH_VERT_GROWTH 
!
!***********************************************************************
! Compute vertical biomass production of marsh.  
!***********************************************************************
!
      CALL marsh_vert_growth(ng, tile)
# endif 
!
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)!
# endif
      RETURN
      END SUBROUTINE marsh_dynamics
#endif
      END MODULE marsh_dynamics_mod
