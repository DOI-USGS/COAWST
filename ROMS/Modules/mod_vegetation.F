#include "cppdefs.h"
      MODULE mod_vegetation 

/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
** This module declares submerged aquatic vegetation model state      **
** variables and internal parameters. Some of these parameters are    **
** usually read from the appropriate input script.                    **
**                                                                    **
** The current design allows the user to have a lot of latitude for   **
** customizing or expanding the vegetation model.                     **
**                                                                    **
** The vegetatation model is composed of several files:               **
**                                                                    **
**   *  Vegetation modifies rhs3d terms:                              **
**                                                                    **
**        vegetation_drag.F                                           **
**                                                                    **
**   *  Vegetation modifies turbulence terms:                         **
**                                                                    **
**        vegetation_turb.F                                           **
**                                                                    **
**   *  Vegetation modifies streaming terms:                          **
**                                                                    **
**        vegetation_stream.F                                         **
**                                                                    **
**   *  Vegetation biomass calculation:                               **
**                                                                    **
**        vegetation_biomass.F                                        **
**                                                                    **
**   *  Vegetation output variables into NetCDF files:                **
**                                                                    **
**        vegetation_output.F                                         **
**                                                                    **
**   *  Marsh dynamics (main calling routine) calculation:            **
**                                                                    **
**        marsh_dynamics.F                                            **
**                                                                    **
**   *  Marsh wave thrust calculation:                                **
**                                                                    **
**        marsh_wave_thrust.F                                         **
**                                                                    **
**   *  Marsh sediment erosion calculations:                          **
**                                                                    **
**        marsh_sed_erosion.F                                         **
**                                                                    **
**   *  Marsh growth mean high water/tidal calculations:              **
**                                                                    **
**        marsh_tidal_range.F                                         **
**                                                                    **
**   *  Marsh vertical growth calculations:                           **
**                                                                    **
**        marsh_vert_growth.F                                         **
**                                                                    **
**   *  Internal model parameters declarations:                       **
**                                                                    **
**        vegetation_mod.h                                            **
**                                                                    **
**   *  Model parameters standard input script:                       **
**                                                                    **
**        vegetation.in                                               **
**                                                                    **
**   *  Code to read input model parameters:                          **
**                                                                    **
**        vegetation_inp.h                                            **
**                                                                    **
**   *  Code to assign indices to model variables during the          **
**      reading of metadata information from "varinfo.yaml":          **
**                                                                    **
**        vegetation_var.h                                            **
**                                                                    **
** Note that all the files are located in ROMS/Nonlinear/Vegetation   **
** and the *.h files are included within  <...> to allow the user     **
** customize any of them in the project directory, while  keeping     **
** the distributed code intact (check build scripts for details).     **
**                                                                    **
************************************************************************
*/

#if defined VEGETATION 
# include <vegetation_mod.h>
#endif

      END MODULE mod_vegetation
