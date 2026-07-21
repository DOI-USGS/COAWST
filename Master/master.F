#include "cppdefs.h"
/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Master program to run ROMS as single ocean model or coupled       **
**  to other models using the MCT or ESMF libraries.                  **
**                                                                    **
************************************************************************
*/

#if defined COAWST_COUPLING && !defined CMEPS
# if defined MCT_LIB
#  include "mct_driver.h"
# elif defined ESMF_LIB
#  include "esmf_driver.h"
# else
      PROGRAM master         !  ROMS is coupled with other external
      END PROGRAM master     !  software and drivers not written by us
# endif
#elif defined ROMS_MODEL
# include "roms.h"
#elif (defined SWAN_MODEL || defined WW3_MODEL)
# include "waves.h"
#endif
