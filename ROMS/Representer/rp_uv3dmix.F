#include "cppdefs.h"
      MODULE rp_uv3dmix_mod

#if defined TL_IOMS && defined SOLVE3D && \
   (defined UV_VIS2 || defined UV_VIS4)
!
!svn $Id: rp_uv3dmix.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine computes representers tangent linear horizontal     !
!  viscosity of momentum.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
# ifdef UV_VIS2
      PUBLIC rp_uv3dmix2
# endif
# ifdef UV_VIS4
      PUBLIC rp_uv3dmix4
# endif

      CONTAINS

# ifdef UV_VIS2
#  if defined MIX_S_UV
#   include "rp_uv3dmix2_s.h"
#  elif defined MIX_GEO_UV
#   include "rp_uv3dmix2_geo.h"
#  endif
# endif

# ifdef UV_VIS4
#  if defined MIX_S_UV
#   include "rp_uv3dmix4_s.h"
#  elif defined MIX_GEO_UV
#   include "rp_uv3dmix4_geo.h"
#  endif
# endif

#endif

      END MODULE rp_uv3dmix_mod
