#include "cppdefs.h"
#if defined NONLINEAR && defined SOLVE3D
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!! This routine computes the baroclinic hydrostatic pressure gradient  !
!! term.                                                               !
!!                                                                     !
!!======================================================================
!!
# if defined PJ_GRADPQ4
#  include "prsgrd44.h"
# elif defined PJ_GRADPQ2
#  include "prsgrd42.h"
# elif defined PJ_GRADP
#  include "prsgrd40.h"
# elif defined DJ_GRADPS
#  include "prsgrd32.h"
# else
#  include "prsgrd31.h"
# endif
#else
      MODULE prsgrd_mod
      END MODULE prsgrd_mod
#endif
