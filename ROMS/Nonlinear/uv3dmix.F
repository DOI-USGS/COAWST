#include "cppdefs.h"

#if defined SOLVE3D && defined NONLINEAR && \
   (defined UV_VIS2 || defined UV_VIS4)

!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!! This routine computes horizontal viscosity of momentum.             !
!!                                                                     !
!!======================================================================
!!
# ifdef UV_VIS2
#  if defined MIX_S_UV
#   include "uv3dmix2_s.h"
#  elif defined MIX_GEO_UV
#   include "uv3dmix2_geo.h"
#  else
      UV3DMIX: must define one of MIX_S_UV, MIX_GEO_UV
#  endif
# endif

# ifdef UV_VIS4
#  if defined MIX_S_UV
#   include "uv3dmix4_s.h"
#  elif defined MIX_GEO_UV
#   include "uv3dmix4_geo.h"
#  else
      UV3DMIX: must define one of MIX_S_UV, MIX_GEO_UV
#  endif
# endif
#else
      MODULE uv3dmix2_mod
      END MODULE uv3dmix2_mod

      MODULE uv3dmix4_mod
      END MODULE uv3dmix4_mod
#endif
