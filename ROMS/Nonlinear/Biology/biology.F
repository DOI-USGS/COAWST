#include "cppdefs.h"
      MODULE biology_mod
#if defined NONLINEAR && defined BIOLOGY
!
!svn $Id: biology.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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

# if defined BIO_FENNEL
#  include <fennel.h>
# elif defined ECOSIM
#  include <ecosim.h>
# elif defined NEMURO
#  include <nemuro.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks.h>
# elif defined NPZD_IRON
#  include <npzd_iron.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell.h>
# endif

#endif
      END MODULE biology_mod
