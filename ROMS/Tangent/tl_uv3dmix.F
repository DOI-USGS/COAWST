#include "cppdefs.h"

#if defined TANGENT && defined SOLVE3D && \
   (defined UV_VIS2 || defined UV_VIS4)
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!! This routine computes tangent linear horizontal viscosity of        !
!! momentum.                                                           !
!!                                                                     !
!!======================================================================
!!
# ifdef UV_VIS2
#  if defined MIX_S_UV
#   include "tl_uv3dmix2_s.h"
#  elif defined MIX_GEO_UV
#   include "tl_uv3dmix2_geo.h"
#  endif
# endif

# ifdef UV_VIS4
#  if defined MIX_S_UV
#   include "tl_uv3dmix4_s.h"
#  elif defined MIX_GEO_UV
#   include "tl_uv3dmix4_geo.h"
#  endif
# endif
#else
      MODULE tl_uv3dmix2_mod
      END MODULE tl_uv3dmix2_mod

      MODULE tl_uv3dmix4_mod
      END MODULE tl_uv3dmix4_mod
#endif
