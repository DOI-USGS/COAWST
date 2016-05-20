#include "cppdefs.h"

      MODULE sediment_mod

#if defined NONLINEAR && defined SEDIMENT
!
!svn $Id: sediment.F 429 2009-12-20 17:30:26Z arango $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine it is the main driver for the sediment-transport       !
!  model. Currently, it includes calls to the following routines:      !
!                                                                      !
!  * Vertical settling of sediment in the water column.                !
!  * Erosive and depositional flux interactions of sediment            !
!    between water column and the bed.                                 !
!  * Transport of multiple grain sizes.                                !
!  * Bed layer stratigraphy.                                           !
!  * Bed morphology.                                                   !
!  * Bedload based on Meyer Peter Mueller.                             !
!  * Bedload based on Soulsby combined waves + currents                !
!    (p166 Soulsby 1997)                                               !
!  * Bedload slope term options: Nemeth et al, 2006, Coastal           !
!    Engineering, v 53, p 265-275; Lesser et al, 2004, Coastal         !
!    Engineering, v 51, p 883-915.                                     !
!                                                                      !
!  * Seawater/sediment vertical level distribution:                    !
!                                                                      !
!         W-level  RHO-level                                           !
!                                                                      !
!            N     _________                                           !
!                 |         |                                          !
!                 |    N    |                                          !
!          N-1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |   N-1   |  A                                       !
!            2    |_________|  W                                       !
!                 |         |  A                                       !
!                 |    2    |  T                                       !
!            1    |_________|  E                                       !
!                 |         |  R                                       !
!                 |    1    |                                          !
!            0    |_________|_____ bathymetry                          !
!                 |/////////|                                          !
!                 |    1    |                                          !
!            1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |    2    |  D                                       !
!            2    |_________|  I                                       !
!                 |         |  M                                       !
!                 |  Nbed-1 |  E                                       !
!        Nbed-1   |_________|  N                                       !
!                 |         |  T                                       !
!                 |  Nbed   |                                          !
!         Nbed    |_________|                                          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sediment

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sediment (ng, tile)
!***********************************************************************
!
# if defined COHESIVE_BED  || defined MIXED_BED
      USE sed_bed_cohesive_mod, ONLY : sed_bed_cohesive
# elif defined NONCOHESIVE_BED2
      USE sed_bed_mod2, ONLY : sed_bed2
# else
      USE sed_bed_mod, ONLY : sed_bed
# endif
# ifdef BEDLOAD
      USE sed_bedload_mod, ONLY : sed_bedload
# endif
# if defined SED_BIODIFF
    USE sed_biodiff_mod, ONLY : sed_biodiff
# endif
# ifdef SUSPLOAD
      USE sed_fluxes_mod, ONLY : sed_fluxes
      USE sed_settling_mod, ONLY : sed_settling
#  if defined SED_FLOCS
    USE sed_flocs_mod, ONLY : sed_flocmod
#  endif
# endif
      USE sed_surface_mod, ONLY : sed_surface
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile

# ifdef BEDLOAD
!
!-----------------------------------------------------------------------
!  Compute sediment bedload transport.
!-----------------------------------------------------------------------
!
      CALL sed_bedload (ng, tile)
# endif

# ifdef SUSPLOAD
#  ifdef SED_FLOCS
!
!-----------------------------------------------------------------------
!  Compute sediment flocculation
!-----------------------------------------------------------------------
!
      CALL sed_flocmod (ng, tile)
#  endif
!
!-----------------------------------------------------------------------
!  Compute sediment vertical settling.
!-----------------------------------------------------------------------
!
      CALL sed_settling (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute bed-water column exchanges: erosion and deposition.
!-----------------------------------------------------------------------
!
      CALL sed_fluxes (ng, tile)
!
# endif
!
!-----------------------------------------------------------------------
!  Compute sediment bed stratigraphy.
!-----------------------------------------------------------------------
!
# if defined COHESIVE_BED || defined MIXED_BED
      CALL sed_bed_cohesive (ng, tile)
# elif defined NONCOHESIVE_BED2
      CALL sed_bed2 (ng, tile)
# else
      CALL sed_bed (ng, tile)
# endif
!
!-----------------------------------------------------------------------
!  Compute sediment bed biodiffusivity.
!-----------------------------------------------------------------------
!
# if defined SED_BIODIFF
      CALL sed_biodiff (ng, tile)
# endif
!
!-----------------------------------------------------------------------
!  Compute sediment surface layer properties.
!-----------------------------------------------------------------------
!
      CALL sed_surface (ng, tile)

      RETURN
      END SUBROUTINE sediment
#endif
      END MODULE sediment_mod
