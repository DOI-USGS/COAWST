#include "cppdefs.h"
      MODULE tl_prsgrd_mod
#if defined TANGENT && defined SOLVE3D
!
!svn $Id: tl_prsgrd.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the tangent linear baroclinic hydrostatic     !
!  pressure gradient term.                                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: tl_prsgrd

      CONTAINS

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

#endif

      END MODULE tl_prsgrd_mod
