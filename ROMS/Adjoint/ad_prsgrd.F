#include "cppdefs.h"
#if defined ADJOINT && defined SOLVE3D
!!
!!git $Id$
!!svn $Id: ad_prsgrd.F 1151 2023-02-09 03:08:53Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!!                                                                     !
!! This routine computes the adjoint baroclinic hydrostatic pressure   !
!! gradient term.                                                      !
!!                                                                     !
!!======================================================================
!!
# if defined PJ_GRADPQ4
!!#  include "ad_prsgrd44.h"
# elif defined PJ_GRADPQ2
!!#  include "ad_prsgrd42.h"
# elif defined PJ_GRADP
#  include "ad_prsgrd40.h"
# elif defined DJ_GRADPS
#  include "ad_prsgrd32.h"
# else
#  include "ad_prsgrd31.h"
# endif
#else
      MODULE ad_prsgrd_mod
      END MODULE ad_prsgrd_mod
#endif
