#include "cppdefs.h"
#ifdef SEAICE
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!!  This routine reads and reports seaice model input parameters.      !
!!                                                                     !
!!======================================================================
!!
# if defined ICE_MODEL
#  include <ice_inp.h>
# endif
#else
      SUBROUTINE read_IcePar
      END SUBROUTINE read_IcePar
#endif
