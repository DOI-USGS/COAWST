#include "cppdefs.h"
#ifdef SEDIMENT
!!
!!svn $Id: read_sedpar.F 830 2017-01-24 21:21:11Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
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
