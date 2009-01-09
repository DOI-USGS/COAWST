#include "cppdefs.h"
      MODULE biology_mod
#if defined NONLINEAR && defined BIOLOGY
!
!svn $Id: biology.F 1311 2008-01-10 04:13:52Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the source and sink terms for selected        !
!   biology model.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: biology

      CONTAINS

# if defined BIO_FASHAM
#  include "fasham.h"
# elif defined NEMURO
#  include "nemuro.h"
# elif defined NPZD_FRANKS
#  include "npzd_Franks.h"
# elif defined NPZD_POWELL
#  include "npzd_Powell.h"
# elif defined ECOSIM
#  include "ecosim.h"
# endif

#endif
      END MODULE biology_mod
