#include "cppdefs.h"
#ifdef ADJOINT
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group            Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This subroutine performs a fast (predictor or corrector) time-step  !
!  for the free-surface and 2D momentum adjoint equations.             !
# ifdef SOLVE3D
!  It also calculates the time filtering variables over all fast-time  !
!  steps to damp high frequency signals in 3D applications.            !
# endif
!                                                                      !
!=======================================================================
!
# include "ad_step2d_LF_AM3.h"
#else
      MODULE ad_step2d_mod
      END MODULE ad_step2d_mod
#endif
