#include "cppdefs.h"
#if defined NONLINEAR && defined BBL_MODEL
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!! This routine computes bottom momentum stress via a bottom boundary  !
!! layer formulation.                                                  !
!!                                                                     !
!!======================================================================
!!
# if defined SSW_BBL
#  include "ssw_bbl.h"
# elif defined MB_BBL
#  include "mb_bbl.h"
# elif defined SG_BBL
#  include "sg_bbl.h"
# else
      MODULE bbl_mod
      END MODULE bbl_mod
# endif
#else
      MODULE bbl_mod
      END MODULE bbl_mod
#endif
