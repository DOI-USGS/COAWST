#include "cppdefs.h"
      MODULE ad_biology_mod
#if defined ADJOINT && defined BIOLOGY
!
!svn $Id: ad_biology.F 503 2008-01-10 00:11:51Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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

# if defined BIO_FASHAM
!!#  include "ad_fasham.h"
# elif defined NEMURO
!!#  include "ad_nemuro.h"
# elif defined NPZD_FRANKS
#  include "ad_npzd_Franks.h"
# elif defined NPZD_POWELL
#  include "ad_npzd_Powell.h"
# endif

#endif
      END MODULE ad_biology_mod
