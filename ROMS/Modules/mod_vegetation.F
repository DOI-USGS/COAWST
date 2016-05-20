#include "cppdefs.h"
      MODULE mod_vegetation 

/*
** svn $Id: mod_vegetation.F 429 2015-05-20 10:10:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This module declares vegetation model internal parameters.  Some   **
** of these parameters are usually read from the appropriate input    **
** script.                                                            **
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
**        vegetation_turb_cal.F                                       **
**                                                                    **
**   *  Vegetation modifies streaming terms:                          **
**                                                                    **
**        vegetation_stream.F                                         **
**                                                                    **
**   *  Wave Thrust on Marsh calculation:                             **
**                                                                    **
**        marsh_wave_thrust.F                                         ** 
**                                                                    **
**   *  Vegetation biomass calculation:                               **
**                                                                    **
**        vegetation_biomass.F                                        **
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
**      reading of metadata information from "varinfo.dat":           **
**                                                                    **
**        vegetation_var.h                                            **
**                                                                    **
**   *  Code to define input model parameters in all output           **
**      NetCDF files:                                                 **
**                                                                    **
**        vegetation_def_his.h                                        **
**                                                                    **
**   *  Code to write out input model parameters in all output        **
**      NetCDF files:                                                 **
**                                                                    **
**        vegetation_wrt_his.h                                        **
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
