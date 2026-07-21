#include "cppdefs.h"
#if defined TL_IOMS && defined SOLVE3D
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!! This routine computes the representers tangent linear baroclinic    !
!! hydrostatic pressure gradient term.                                 !
!!                                                                     !
!!======================================================================
!!
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
#else
      MODULE rp_prsgrd_mod
      END MODULE rp_prsgrd_mod
#endif
