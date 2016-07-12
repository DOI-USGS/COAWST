#include "cppdefs.h"

      MODULE rp_t3dmix_mod

#ifdef TL_IOMS
# if !defined TS_FIXED && defined SOLVE3D &&  \
     (defined TS_DIF2  || defined TS_DIF4)

!
!svn $Id: rp_t3dmix.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes representers tangent linear horizontal mixing !
!  of tracers.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
#  ifdef TS_DIF2
      PUBLIC rp_t3dmix2
#  endif
#  ifdef TS_DIF4
      PUBLIC rp_t3dmix4
#  endif

      CONTAINS

#  ifdef TS_DIF2
#   if defined MIX_S_TS
#    include "rp_t3dmix2_s.h"
#   elif defined MIX_GEO_TS
#    include "rp_t3dmix2_geo.h"
#   elif defined MIX_ISO_TS
#    include "rp_t3dmix2_iso.h"
#   endif
#  endif

#  ifdef TS_DIF4
#   if defined MIX_S_TS
#    include "rp_t3dmix4_s.h"
#   elif defined MIX_GEO_TS
#    include "rp_t3dmix4_geo.h"
#   elif defined MIX_ISO_TS
#    include "rp_t3dmix4_iso.h"
#   endif
#  endif

# endif
#endif

      END MODULE rp_t3dmix_mod
