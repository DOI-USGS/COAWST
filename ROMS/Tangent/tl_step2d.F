#include "cppdefs.h"
      MODULE tl_step2d_mod
#ifdef TANGENT
/*
** svn $Id: tl_step2d.F 503 2008-01-10 00:11:51Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2008 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************ 
**                                                                    **
** This routine performs a fast (predictor or corrector) time-step    **
** for the free-surface and 2D momentum tangent linear equations.     **
**                                                                    **
************************************************************************
*/
      implicit none

      PRIVATE
      PUBLIC  :: tl_step2d

      CONTAINS

# include "tl_step2d_LF_AM3.h"

#endif
      END MODULE tl_step2d_mod
