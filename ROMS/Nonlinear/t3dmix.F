#include "cppdefs.h"

      MODULE t3dmix_mod

#if !defined TS_FIXED && defined SOLVE3D && defined NONLINEAR && \
    (defined TS_DIF2  || defined TS_DIF4)

!
!svn $Id: t3dmix.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes horizontal mixing of tracers.                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
# ifdef TS_DIF2
      PUBLIC t3dmix2
# endif
# ifdef TS_DIF4
      PUBLIC t3dmix4
# endif

      CONTAINS

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

#endif

      END MODULE t3dmix_mod
