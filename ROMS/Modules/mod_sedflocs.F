#include "cppdefs.h"
      MODULE mod_sedflocs

/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
** This module declares  Sediment Floc Model Kernel Structure         **
** containing several variables describing flocculation interaction   **
** properties.                                                        **
**                                                                    **
** Note that the  *.h file is located in ROMS/Nonlinear/Sediment and  **
** it is included within <...> to allow the user customize it in the  **
** project directory, while keeping the distributed file intact       **
** (check build scripts for details).                                 **
**                                                                    **
************************************************************************
*/

#if defined SEDIMENT
# include <sedflocs_mod.h>
#endif

      END MODULE mod_sedflocs
