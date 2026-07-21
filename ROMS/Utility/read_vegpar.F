#include "cppdefs.h"
#ifdef VEGETATION
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2026 The ROMS Group                              !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!!  This routine reads and reports vegetation model input parameters.  !
!!                                                                     !
!!======================================================================
!!
# include <vegetation_inp.h>
#else
      SUBROUTINE read_vegpar
      END SUBROUTINE read_vegpar
#endif
