#include "cppdefs.h"
/*
** svn $Id: master.F 1054 2021-03-06 19:47:12Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2021 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Master program to run ROMS/TOMS as single ocean model or coupled  **
**  to other models using the MCT or ESMF libraries.                  **
**                                                                    **
************************************************************************
*/

#if defined COAWST_COUPLING
# if defined MCT_LIB
#  include "mct_driver.h"
# elif defined ESMF_LIB
#  include "esmf_driver.h"
# else
      PROGRAM master         !  ROMS is coupled with other external
      END PROGRAM master     !  software and drivers not written by us
# endif
#elif defined ROMS_MODEL
# include "ocean.h"
#elif (defined SWAN_MODEL || defined WW3_MODEL)
# include "waves.h"
#endif
