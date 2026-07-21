#include "cppdefs.h"
/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
** This module declares Lagrangian drifters biological behavior       **
** parameteres. They are usually read from the appropriate input      **
** script.                                                            **
**                                                                    **
** The current design allows the user to have a lot of flexibility    **
** for customizing or adding any particular behavior model.           **
**                                                                    **
** Note that all the *.h files are located in ROMS/Nonlinear/Biology  **
** and included within  <...>  to allow the user to customize any of  **
** them in the project directory  while keeping the distributed code  **
** intact (check the build script for details).                       **
**                                                                    **
************************************************************************
*/

#if defined FLOATS && defined FLOAT_BIOLOGY
# if defined FLOAT_OYSTER
#  include <oyster_floats_mod.h>
#  else
      MODULE mod_behavior
      END MODULE mod_behavior
# endif
#else
      MODULE mod_behavior
      END MODULE mod_behavior
#endif
