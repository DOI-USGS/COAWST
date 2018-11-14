#include "cppdefs.h"
#ifdef VEGETATION
!!
!!svn $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
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
