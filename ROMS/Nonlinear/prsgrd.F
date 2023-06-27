#include "cppdefs.h"
#if defined NONLINEAR && defined SOLVE3D
!!
!!git $Id$
!!svn $Id: prsgrd.F 1151 2023-02-09 03:08:53Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
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
