#include "cppdefs.h"

/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
** This module declares  Sediment Model Kernel Structure  containing  **
** several variables describing bed layers, stratigraphy, and bottom  **
** properties.                                                        **
**                                                                    **
** Note that the  *.h file is located in ROMS/Nonlinear/Sediment and  **
** it is included within <...> to allow the user customize it in the  **
** project directory, while keeping the distributed file intact       **
** (check build scripts for details).                                 **
**                                                                    **
************************************************************************
*/

#if defined SEDIMENT || defined BBL_MODEL
# include <sedbed_mod.h>
#else
      MODULE mod_sedbed
      END MODULE mod_sedbed
#endif
