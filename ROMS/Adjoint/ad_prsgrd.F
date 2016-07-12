#include "cppdefs.h"
      MODULE ad_prsgrd_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_prsgrd.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the adjoint baroclinic hydrostatic pressure   !
!  gradient term.                                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_prsgrd

      CONTAINS

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

#endif

      END MODULE ad_prsgrd_mod
