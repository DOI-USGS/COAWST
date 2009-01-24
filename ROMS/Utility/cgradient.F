#include "cppdefs.h"
#define ORTHOGONALIZATION
#ifdef ORTHOGONALIZATION
# undef TEST_ORTHOGONALIZATION
# undef NORMALIZATION
#endif

      MODULE cgradient_mod
!!
!!svn $Id: cgradient.F 503 2008-01-10 00:11:51Z arango $
!!================================================== Hernan G. Arango ==
!! Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!!                                                                     !
!! This module minimizes a quadratic cost function using the conjugate !
!! gradient algorithm proposed by Mike Fisher (ECMWF).                 !
!!                                                                     !
!! Reference:                                                          !
!!                                                                     !
!!   Fisher, M., 1997: Efficient Minimization of Quadratic Penalty     !
!!     funtions, unpublish manuscript, 1-14.                           !
!!                                                                     !
!!======================================================================
!!

#if defined IS4DVAR
# if defined LANCZOS
#  include "cgradient_lanczos.h"
# else
#  include "cgradient.h"
# endif
#endif

      END MODULE cgradient_mod
