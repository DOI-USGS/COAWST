#include "cppdefs.h"
      MODULE step2d_mod
#ifdef NONLINEAR
/*
** svn $Id: step2d.F 1311 2008-01-10 04:13:52Z jcwarner $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2008 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This routine performs a fast (predictor or corrector) time-step    **
** for the free-surface and 2D momentum nonlinear equations.          **
**                                                                    **
************************************************************************
*/
      implicit none

      PRIVATE
      PUBLIC  :: step2d

      CONTAINS

# include "step2d_LF_AM3.h"

#endif
      END MODULE step2d_mod
