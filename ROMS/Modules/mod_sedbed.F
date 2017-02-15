#include "cppdefs.h"
      MODULE mod_sedbed

/*
** svn $Id: mod_sedbed.F 830 2017-01-24 21:21:11Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2017 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
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
#endif

      END MODULE mod_sedbed
