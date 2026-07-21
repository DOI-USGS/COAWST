#include "cppdefs.h"

#ifdef ADJOINT
# if !defined TS_FIXED && defined SOLVE3D &&  \
     (defined TS_DIF2  || defined TS_DIF4)
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!! This routine computes adjoint horizontal mixing of tracers.         !
!!                                                                     !
!!======================================================================
!!
#  ifdef TS_DIF4
#   if defined MIX_S_TS
#    include "ad_t3dmix4_s.h"
#   elif defined MIX_GEO_TS
#    include "ad_t3dmix4_geo.h"
#   elif defined MIX_ISO_TS
#    include "ad_t3dmix4_iso.h"
#   endif
#  endif

#  ifdef TS_DIF2
#   if defined MIX_S_TS
#    include "ad_t3dmix2_s.h"
#   elif defined MIX_GEO_TS
#    include "ad_t3dmix2_geo.h"
#   elif defined MIX_ISO_TS
#    include "ad_t3dmix2_iso.h"
#   endif
#  endif
# endif
#else
      MODULE ad_t3dmix2_mod
      END MODULE ad_t3dmix2_mod

      MODULE ad_t3dmix4_mod
      END MODULE ad_t3dmix4_mod
#endif
