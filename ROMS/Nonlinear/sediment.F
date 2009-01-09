#include "cppdefs.h"
      MODULE sediment_mod
#if defined NONLINEAR && defined SEDIMENT
!
!svn $Id: sediment.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This  routine computes the sediment sources and sinks and adds      !
!  then the global sediment tracer fields. Currently, it includes      !
!  calls to the following routines:                                    !
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
# ifdef BEDLOAD
      USE bedload_mod, ONLY : bedload
# endif
# ifdef SUSPLOAD
      USE settling_mod, ONLY : settling
      USE dep_ero_mod, ONLY : dep_ero
# endif
# if defined BED_COHESIVE
      USE bed_cohesive_mod, ONLY : bed_cohesive
# else
      USE bed_mod, ONLY : bed
# endif
# if defined SED_BIODIFF
      USE bedbiodiff_mod, ONLY : bedbiodiff
# endif
      USE bottom_mod, ONLY : bottom
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
!-----------------------------------------------------------------------
!  Compute bedload transport.
!-----------------------------------------------------------------------
!
# ifdef BEDLOAD
      CALL bedload (ng, tile)
# endif
!
!-----------------------------------------------------------------------
!  Compute vertical settling.
!-----------------------------------------------------------------------
!
# ifdef SUSPLOAD
      CALL settling (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute bed-water column exchanges.
!-----------------------------------------------------------------------
!
      CALL dep_ero (ng, tile)
# endif
!
!-----------------------------------------------------------------------
!  Compute bed stratigraphy.
!-----------------------------------------------------------------------
!
# if defined BED_COHESIVE
     CALL bed_cohesive (ng, tile)
# else
     CALL bed (ng, tile)
# endif
!
!-----------------------------------------------------------------------
!  Compute bed biodiffusivity.
!-----------------------------------------------------------------------
!
# if defined SED_BIODIFF
      CALL bedbiodiff (ng, tile)
# endif
!
!-----------------------------------------------------------------------
!  Compute bottom surface properties.
!-----------------------------------------------------------------------
!
      CALL bottom (ng, tile)
!
      RETURN
      END SUBROUTINE sediment
#endif
      END MODULE sediment_mod
