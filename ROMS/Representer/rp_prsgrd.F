#include "cppdefs.h"
#if defined TL_IOMS && defined SOLVE3D
!!
!!git $Id$
!!svn $Id: rp_prsgrd.F 1151 2023-02-09 03:08:53Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
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
