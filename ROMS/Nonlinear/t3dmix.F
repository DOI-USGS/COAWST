#include "cppdefs.h"

#if !defined TS_FIXED && defined SOLVE3D && defined NONLINEAR && \
    (defined TS_DIF2  || defined TS_DIF4)

!!
!!git $Id$
!!svn $Id: t3dmix.F 1151 2023-02-09 03:08:53Z arango $
!!======================================================================
!! Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                           Hernan G. Arango   !
!!========================================= Alexander F. Shchepetkin ===
!!                                                                     !
!! This routine computes horizontal mixing of tracers.                 !
!!                                                                     !
!!======================================================================
!!
# ifdef TS_DIF2
#  if defined MIX_S_TS
#   include "t3dmix2_s.h"
#  elif defined MIX_GEO_TS
#   include "t3dmix2_geo.h"
#  elif defined MIX_ISO_TS
#   include "t3dmix2_iso.h"
#  else
      T3DMIX: must define one of MIX_S_TS, MIX_GEO_TS, MIX_ISO_TS
#  endif
# endif

# ifdef TS_DIF4
#  if defined MIX_S_TS
#   include "t3dmix4_s.h"
#  elif defined MIX_GEO_TS
#   include "t3dmix4_geo.h"
#  elif defined MIX_ISO_TS
#   include "t3dmix4_iso.h"
#  else
      T3DMIX: must define one of MIX_S_TS, MIX_GEO_TS, MIX_ISO_TS
#  endif
# endif
#else
      MODULE t3dmix2_mod
      END MODULE t3dmix2_mod

      MODULE t3dmix4_mod
      END MODULE t3dmix4_mod
#endif
