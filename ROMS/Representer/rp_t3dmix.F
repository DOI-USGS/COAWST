#include "cppdefs.h"

#ifdef TL_IOMS
# if !defined TS_FIXED && defined SOLVE3D &&  \
     (defined TS_DIF2  || defined TS_DIF4)
!!
!!git $Id$
!!svn $Id: rp_t3dmix.F 1151 2023-02-09 03:08:53Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!!                                                                     !
!! This routine computes representers tangent linear horizontal mixing !
!! of tracers.                                                         !
!!                                                                     !
!!======================================================================
!!
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
#else
      MODULE rp_t3dmix2_mod
      END MODULE rp_t3dmix2_mod

      MODULE rp_t3dmix4_mod
      END MODULE rp_t3dmix4_mod
#endif
