#include "cppdefs.h"
      MODULE rp_prsgrd_mod
#if defined TL_IOMS && defined SOLVE3D
!
!svn $Id: rp_prsgrd.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the representers tangent linear baroclinic    !
!  hydrostatic pressure gradient term.                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: rp_prsgrd

      CONTAINS

# if defined PJ_GRADPQ4
!!#  include "rp_prsgrd44.h"
# elif defined PJ_GRADPQ2
!!#  include "rp_prsgrd42.h"
# elif defined PJ_GRADP
#  include "rp_prsgrd40.h"
# elif defined DJ_GRADPS
#  include "rp_prsgrd32.h"
# else
#  include "rp_prsgrd31.h"
# endif

#endif

      END MODULE rp_prsgrd_mod
