#include "cppdefs.h"

#if defined SOLVE3D && defined NONLINEAR && \
   (defined UV_VIS2 || defined UV_VIS4)

!!
!!git $Id$
!!svn $Id: uv3dmix.F 1151 2023-02-09 03:08:53Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
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
