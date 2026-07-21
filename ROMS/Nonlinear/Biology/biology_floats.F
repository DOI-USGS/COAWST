#include "cppdefs.h"
#if defined NONLINEAR && defined FLOATS && defined FLOAT_BIOLOGY
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  These routine computes the behavior for biological floats.          !
!                                                                      !
!=======================================================================
!
# ifdef FLOAT_OYSTER
#  include <oyster_floats.h>
# endif
#else
      MODULE biology_floats_mod
      END MODULE biology_floats_mod
#endif
