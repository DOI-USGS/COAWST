#include "cppdefs.h"
      MODULE rp_biology_mod
#if defined TL_IOMS && defined BIOLOGY
!
!svn $Id: rp_biology.F 503 2008-01-10 00:11:51Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the representers tangent linear source and    !
!  sink terms for the selected biology model.                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: rp_biology

      CONTAINS

# if defined BIO_FASHAM
!!#  include "rp_fasham.h"
# elif defined NEMURO
!!#  include "rp_nemuro.h"
# elif defined NPZD_FRANKS
#  include "rp_npzd_Franks.h"
# elif defined NPZD_POWELL
#  include "rp_npzd_Powell.h"
# endif

#endif
      END MODULE rp_biology_mod
