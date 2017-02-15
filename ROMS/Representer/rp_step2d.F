#include "cppdefs.h"
      MODULE rp_step2d_mod
#ifdef TL_IOMS
!
!svn $Id: rp_step2d.F 830 2017-01-24 21:21:11Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine performs a fast (predictor or corrector) time-step  !
!  for the free-surface and 2D momentum  representers  tangent linear  !
!  equations.                                                          !
# ifdef SOLVE3D
!  It also calculates the time filtering variables over all fast-time  !
!  steps  to damp high frequency signals in 3D applications.           !
# endif
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rp_step2d

      CONTAINS

# include "rp_step2d_LF_AM3.h"

#endif
      END MODULE rp_step2d_mod
