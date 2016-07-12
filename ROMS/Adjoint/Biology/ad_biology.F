#include "cppdefs.h"
      MODULE ad_biology_mod
#if defined ADJOINT && defined BIOLOGY
!
!svn $Id: ad_biology.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the adjoint source and sink terms for the     !
!  selected biology model.                                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_biology

      CONTAINS

# if defined BIO_FENNEL
!!#  include <ad_fennel.h>
# elif defined NEMURO
!!#  include <ad_nemuro.h>
# elif defined NPZD_FRANKS
#  include <ad_npzd_Franks.h>
# elif defined NPZD_IRON
#  include <ad_npzd_iron.h>
# elif defined NPZD_POWELL
#  include <ad_npzd_Powell.h>
# endif

#endif
      END MODULE ad_biology_mod
