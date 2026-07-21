#include "cppdefs.h"
#if defined TANGENT && defined SOLVE3D
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!! This routine computes the tangent linear baroclinic hydrostatic     !
!! pressure gradient term.                                             !
!!                                                                     !
!!======================================================================
!!
# if defined PJ_GRADPQ4
!!#  include "tl_prsgrd44.h"
# elif defined PJ_GRADPQ2
!!#  include "tl_prsgrd42.h"
# elif defined PJ_GRADP
#  include "tl_prsgrd40.h"
# elif defined DJ_GRADPS
#  include "tl_prsgrd32.h"
# else
#  include "tl_prsgrd31.h"
# endif
#else
      MODULE tl_prsgrd_mod
      END MODULE tl_prsgrd_mod
#endif
