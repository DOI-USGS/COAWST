#include "cppdefs.h"
      MODULE mod_vegarr

/*
** svn $Id: mod_vegarr.F 429 2015-06-10 17:30:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This module declares  Vegetation Model Kernel Structure  containing**
** several variables describing plants, momentum and turbulence terms **
** along wave thrust marsh output arrays.                             **
**                                                                    **
** Note that the  *.h file is located in ROMS/Nonlinear/Vegetation    **
** it is included within <...> to allow the user customize it in the  **
** project directory, while keeping the distributed file intact       **
** (check build scripts for details).                                 **
**                                                                    **
************************************************************************
*/

#if defined VEGETATION
# include <vegarr_mod.h>
#endif

      END MODULE mod_vegarr
