#include "cppdefs.h"
      MODULE mod_biology
/*
** svn $Id: mod_biology.F 795 2016-05-11 01:42:43Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This module declares ecosystem model internal parameteres. They    **
** are usually read from the appropriate input script.                **
**                                                                    **
** The current design allows the user to have a lot of latitude for   **
** customizing or adding any particular ecosystem model.              **
**                                                                    **
** Each ecosystem model is composed of seven files:                   **
**                                                                    **
**   (1) Model source and sink discretized equations:                 **
**                                                                    **
**         ecosim.h           ECOSIM                                  **
**         fennel.h           BIO_FENNEL                              **
**         nemuro.h           NEMURO                                  **
**         npzd_Franks.h      NPZD_FRANKS                             **
**         npzd_iron.h        NPZD_IRON                               **
**         npzd_Powell.h      NPZD_POWELL                             **
**         umaine.h           BIO_UMAINE                              **
**         estuarybgc.h       ESTUARINE BIOGEOCHEMISTRY               **
**         red_tide.h         RED_TIDE                                **
**                                                                    **
**   (2) Internal model parameters declaration:                       **
**                                                                    **
**         ecosim_mod.h                                               **
**         fennel_mod.h                                               **
**         nemuro_mod.h                                               **
**         npzd_Franks_mod.h                                          **
**         npzd_iron_mod.h                                            **
**         npzd_Powell_mod.h                                          **
**         umaine_mod.h                                               **
**         estuarybgc_mod.h                                           **
**         red_tide_mod.h                                             **
**                                                                    **
**   (3) Model parameters standard input script:                      **
**                                                                    **
**         ecosim.in                                                  **
**         bio_Fennel.in                                              **
**         nemuro.in                                                  **
**         npzd_Franks.in                                             **
**         npzd_iron.in                                               **
**         npzd_Powell.in                                             **
**         red_tide.in                                                **
**                                                                    **
**   (4) Code to read input model parameters:                         **
**                                                                    **
**         ecosim_inp.h                                               **
**         fennel_inp.h                                               **
**         nemuro_inp.h                                               **
**         npzd_Franks_inp.h                                          **
**         npzd_iron_inp.h                                            **
**         npzd_Powell_inp.h                                          **
**         umaine_inp.h                                               **
**         estuarybgc_inp.h                                           **
**         red_tide_inp.h                                             **
**                                                                    **
**   (5) Code to assign indices to model variables during the         **
**       reading of metadata information from "varinfo.dat":          **
**                                                                    **
**         ecosim_var.h                                               **
**         fennel_var.h                                               **
**         nemuro_var.h                                               **
**         npzd_Franks_var.h                                          **
**         npzd_iron_var.h                                            **
**         npzd_Powell_var.h                                          **
**         umaine_var.h                                               **
**         estuarybgc_var.h                                           **
**         red_tide_var.h                                             **
**                                                                    **
**   (6) Code to define input model parameters in all output          **
**       NetCDF files:                                                **
**                                                                    **
**         ecosim_def.h                                               **
**         fennel_def.h                                               **
**         nemuro_def.h                                               **
**         npzd_Franks_def.h                                          **
**         npzd_iron_def.h                                            **
**         npzd_Powell_def.h                                          **
**         umaine_def.h                                               **
**         estuarybgc_def.h                                           **
**         red_tide_def.h                                             **
**                                                                    **
**   (7) Code to write out input model parameters in all output       **
**       NetCDF files:                                                **
**                                                                    **
**         ecosim_wrt.h                                               **
**         fennel_wrt.h                                               **
**         nemuro_wrt.h                                               **
**         npzd_Franks_wrt.h                                          **
**         npzd_iron_wrt.h                                            **
**         npzd_Powell_wrt.h                                          **
**         umaine_wrt.h                                               **
**         estuarybgc_wrt.h                                           **
**         red_tide_wrt.h                                             **
**                                                                    **
**                                                                    **
** Note that all the *.h files are located in ROMS/Nonlinear/Biology  **
** and included within  <...>  to allow the user to customize any of  **
** them in the project directory  while keeping the distributed code  **
** intact (check the build script for details).                       **
**                                                                    **
************************************************************************
*/

#ifdef BIOLOGY
# if defined BIO_FENNEL
#  include <fennel_mod.h>
# elif defined ESTUARYBGC
#  include <estuarybgc_mod.h>
# elif defined ECOSIM
#  include <ecosim_mod.h>
# elif defined NEMURO
#  include <nemuro_mod.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks_mod.h>
# elif defined NPZD_IRON
#  include <npzd_iron_mod.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell_mod.h>
# elif defined BIO_UMAINE
#  include <umaine_mod.h>
# elif defined RED_TIDE
#  include <red_tide_mod.h>
# endif
#endif

      END MODULE mod_biology
