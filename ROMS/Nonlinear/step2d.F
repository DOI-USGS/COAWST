#include "cppdefs.h"
#ifdef NONLINEAR
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.md                    Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This subroutine performs a fast (predictor or corrector) time-step  !
!  for the free-surface  and 2D momentum nonlinear equations.
# ifdef SOLVE3D
!  It also calculates the time filtering variables over all fast-time  !
!  steps  to damp high frequency signals in 3D applications.           !
# endif
!                                                                      !
!=======================================================================
!
# include "step2d_LF_AM3.h"
#else
      MODULE step2d_mod
      END MODULE step2d_mod
#endif
