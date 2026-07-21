#include "cppdefs.h"
       MODULE marsh_dynamics_mod

#if defined VEGETATION && defined MARSH_DYNAMICS
!
!git $Id$
!================================================ Tarandeep S. Kalra ===
!  Copyright (c) 2002-2026 The ROMS Group              Neil K. Ganju   !
!    Licensed under a MIT/X style license             John C. Warner   !
!    See License_ROMS.md                           Julia M. Moriarty   !
!============================================== Alfredo Aretxabaleta ===
!                                                                      !
!  This routine is the main driver for marsh dynamics and accounts     !
!  for two processes through various sub-routines:                     !
!                                                                      !
!   1. Calculate lateral wave thrust and sediment release from         !
!      edge erosion                                                    !
!                                                                      !
!   2. Calculate vertical biomass production (autochthonous supply)    !
!      (organic sediment)                                              !
!                                                                      !
!  This is done through the following modules:                         !
!                                                                      !
!    marsh_wave_thrust.F      Computes wave thrust on marshes          !
!                                                                      !
!    marsh_sed_erosion.F      Computes sediment release (i.e. lateral  !
!                             erosion) out of marsh cells              !
!                                                                      !
!    marsh_tidal_range.F      Computes the mean tidal range and marsh  !
!                             high water                               !
!                                                                      !
!    marsh_vert_growth.F      Computes marsh biomass production and    !
!                             vertical accretion                       !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Beudin, A.,  Kalra, T.S., Ganju, N.K., Warner, J.C., 2017:          !
!    Development of a coupled wave-flow-vegetation interaction         !
!    model, Computers & Geosciences, Vol 100, 76-86,                   !
!    doi:10.1016/j.cageo.2016.12.010.                                  !
!                                                                      !
!  Kalra, T.S., Ganju, N.K., Aretxabaleta, A.L., Carr, J.A., Zafer,    !
!    D., Moriarty, J.M., 2021: Modeling Marsh Dynamics Using a 3-D     !
!    Coupled Wave-Flow-Sediment Model, Front. Mar. Sci., Vol 8,        !
!    doi:10.3389/fmars.2021.740921.                                    !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: marsh_dynamics
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE marsh_dynamics (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_stepping
      USE mod_scalars
!
# if defined MARSH_WAVE_THRUST
      USE marsh_wave_thrust_mod,  ONLY : marsh_wave_thrust
# endif
# if defined MARSH_SED_EROSION
      USE marsh_sed_erosion_mod,  ONLY : marsh_sed_erosion
# endif
# if defined MARSH_TIDAL_RANGE
      USE marsh_tidal_range_mod,  ONLY : marsh_tidal_range
# endif
# if defined MARSH_VERT_GROWTH
      USE marsh_vert_growth_mod,  ONLY : marsh_vert_growth
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!-----------------------------------------------------------------------
!  Wave induced thrust calculation on marsh boundary.
!-----------------------------------------------------------------------
!
# if defined MARSH_WAVE_THRUST
      CALL marsh_wave_thrust (ng, tile)
# endif

# if defined MARSH_SED_EROSION
!
!-----------------------------------------------------------------------
!  Marsh sediment erosion based on bathy change in marsh cells.
!-----------------------------------------------------------------------
!
      CALL marsh_sed_erosion (ng, tile)
# endif

# if defined MARSH_TIDAL_RANGE
!
!-----------------------------------------------------------------------
! Compute tidal range and mean high water over a given frequency.
!-----------------------------------------------------------------------
!
      CALL marsh_tidal_range (ng, tile)
#  endif

#  if defined MARSH_VERT_GROWTH
!
!-----------------------------------------------------------------------
! Compute vertical biomass production of marsh.
!-----------------------------------------------------------------------
!
      CALL marsh_vert_growth (ng, tile)
# endif
!
      RETURN
      END SUBROUTINE marsh_dynamics
#endif
      END MODULE marsh_dynamics_mod
