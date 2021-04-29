#include "cppdefs.h"
      MODULE prsgrd_mod
#if defined NONLINEAR && defined SOLVE3D
!
!svn $Id: prsgrd.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the baroclinic hydrostatic pressure gradient  !
!  term.                                                               !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: prsgrd
!
      CONTAINS
!
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

#endif

      END MODULE prsgrd_mod
