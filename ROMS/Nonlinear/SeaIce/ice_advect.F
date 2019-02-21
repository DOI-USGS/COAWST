#include "cppdefs.h"
      MODULE ice_advect_mod
#if defined ICE_MODEL && defined ICE_ADVECT
!
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the advection of the ice tracer fields.       !
!  Currently, the only option is to use the MPDATA advection scheme.   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC ice_advect

      CONTAINS

# if defined ICE_SMOLAR
#  include "ice_smolar.h"
# endif

#endif

      END MODULE ice_advect_mod
