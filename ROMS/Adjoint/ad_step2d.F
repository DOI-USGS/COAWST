#include "cppdefs.h"
      MODULE ad_step2d_mod
#ifdef ADJOINT
!
!svn $Id: ad_step2d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
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
      implicit none
!
      PRIVATE
      PUBLIC  :: ad_step2d

      CONTAINS

# include "ad_step2d_LF_AM3.h"

#endif
      END MODULE ad_step2d_mod
