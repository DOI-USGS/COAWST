#include "cppdefs.h"
/*
** svn $Id: coupler.F 889 2018-02-10 03:32:52Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2019 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Module include several routines used for the couplig of ROMS/TOMS **
**  to other models using the ESMF/NUOPC or MCT libraries.            **
**                                                                    **
************************************************************************
*/

#if defined MODEL_COUPLING && ! defined COAWST_MODEL
# if defined MCT_LIB
#  include "mct_coupler.h"
# elif defined ESMF_LIB
#  include "esmf_coupler.h"
# endif
#else
      MODULE coupler_mod
      END MODULE coupler_mod
#endif
