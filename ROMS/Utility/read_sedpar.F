#include "cppdefs.h"
#ifdef SEDIMENT
!!
!!svn $Id: read_sedpar.F 1054 2021-03-06 19:47:12Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!!                                                                     !
!!  This routine reads and reports sediment model input parameters.    !
!!                                                                     !
!!======================================================================
!!
# include <sediment_inp.h>
#else
      SUBROUTINE read_SedPar
      END SUBROUTINE read_SedPar
#endif
