#include "cppdefs.h"
/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
** This module computes the advection of the ice tracer fields.       **
**                                                                    **
************************************************************************
*/
#if defined ICE_MODEL && defined ICE_ADVECT
# if defined ICE_SMOLAR
#  include <ice_smolar.h>
# endif
#else
      MODULE ice_advect_mod
      END MODULE ice_advect_mod
#endif
