#include "cppdefs.h"

      MODULE uv3dmix_mod

#if defined SOLVE3D && defined NONLINEAR && \
   (defined UV_VIS2 || defined UV_VIS4)

!
!svn $Id: uv3dmix.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes horizontal viscosity of momentum.             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
# ifdef UV_VIS2
      PUBLIC uv3dmix2
# endif
# ifdef UV_VIS4
      PUBLIC uv3dmix4
# endif

      CONTAINS

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

#endif

      END MODULE uv3dmix_mod
