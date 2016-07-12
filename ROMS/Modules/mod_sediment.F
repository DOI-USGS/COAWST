#include "cppdefs.h"
      MODULE mod_sediment

/*
** svn $Id: mod_sediment.F 795 2016-05-11 01:42:43Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This module declares  sediment model internal parameters.  Some    **
** of these parameters are usually read from the appropriate input    **
** script.                                                            **
**                                                                    **
** The current design allows the user to have a lot of latitude for   **
** customizing or expanding the sediment model.                       **
**                                                                    **
** The sediment model is composed of several files:                   **
**                                                                    **
**   *  Model source and sink discretized equations:                  **
**                                                                    **
**        sediment.F           sediment model driver                  **
**        sed_bed.F            bed layer stratigraphy                 **
**        sed_bedload.F        bedload transport                      **
**        sed_fluxes.F         deposition, resuspension, erosion      **
**        sed_settling.F       vertical settling                      **
**        sed_surface.F        sediment-water interface properties    **
**                                                                    **
**   *  Internal model parameters declarations:                       **
**                                                                    **
**        sediment_mod.h                                              **
**                                                                    **
**   *  Model parameters standard input script:                       **
**                                                                    **
**        sediment.in                                                 **
**                                                                    **
**   *  Code to read input model parameters:                          **
**                                                                    **
**        sediment_inp.h                                              **
**                                                                    **
**   *  Code to assign indices to model variables during the          **
**      reading of metadata information from "varinfo.dat":           **
**                                                                    **
**        sediment_var.h                                              **
**                                                                    **
**   *  Code to define input model parameters in all output           **
**      NetCDF files:                                                 **
**                                                                    **
**        sediment_def.h                                              **
**                                                                    **
**   *  Code to write out input model parameters in all output        **
**      NetCDF files:                                                 **
**                                                                    **
**        sediment_wrt.h                                              **
**                                                                    **
** Note that all the files are located in ROMS/Nonlinear/Sediment     **
** and the *.h files are included within  <...> to allow the user     **
** customize any of them in the project directory, while  keeping     **
** the distributed code intact (check build scripts for details).     **
**                                                                    **
************************************************************************
*/

#if defined SEDIMENT || defined BBL_MODEL
# include <sediment_mod.h>
#endif

      END MODULE mod_sediment
