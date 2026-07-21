#include "cppdefs.h"
#ifdef BIOLOGY
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!!  This routine reads and reports selected biological model input     !
!!  parameters.                                                        !
!!                                                                     !
!!======================================================================
!!
# if defined BIO_FENNEL
#  include <fennel_inp.h>
# elif defined ECB
#  include <ecb_inp.h>
# elif defined ESTUARYBGC 
#  include <estuarybgc_inp.h>
# elif defined ECOSIM
#  include <ecosim_inp.h>
# elif defined HYPOXIA_SRM
#  include <hypoxia_srm_inp.h>
# elif defined BIO_GOANPZ
#  include <goanpz_inp.h>
# elif defined BIO_UMAINE
#  include <umaine_inp.h>
# elif defined NEMURO
#  include <nemuro_inp.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks_inp.h>
# elif defined NPZD_IRON
#  include <npzd_iron_inp.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell_inp.h>
# elif defined RED_TIDE
#  include <red_tide_inp.h>
# endif
#else
      SUBROUTINE read_BioPar
      END SUBROUTINE read_BioPar
#endif
