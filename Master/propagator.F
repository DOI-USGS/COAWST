#include "cppdefs.h"
      MODULE propagator_mod
#ifdef PROPAGATOR
!
!svn $Id: propagator.F 503 2008-01-10 00:11:51Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!   These drivers are part of the Generalized Stability Theory (GST)   !
!   package to study the dynamics, sensitivity, and stability of the   !
!   ocean circulation to naturally ocurring perturbations, errors or   !
!   uncetainties in the forecasting system, and adaptive sampling.     !
!                                                                      !
!!  Available Drivers:                                                 !
!!                                                                     !
!!  propagator_adte      Adjoint Finite time Eigenmodes                !
!!  propagator_fte       Finite Time Eigenmodes                        !
!!  propagator_fsv       Forcing singular vectors                      !
!!  propagator_op        Optimal Perturbations, singular vectors       !
!!  propagator_so        Stochactic optimals                           !
!!  propagator_so_semi   Stochactic optimals, semi-norm estimation     !
!!  propagator_so_trace  Stochactic optimals, randomized trace         !
!!                                                                     !
!   Reference:                                                         !
!                                                                      !
!     Moore, A.M. et al., 2004: A comprehensive ocean prediction and   !
!       analysis system based on the tangent linear and adjoint of a   !
!       regional ocean model, Ocean Modelling, 7, 227-258.             !
!                                                                      !
!=======================================================================
!
      implicit none

      PUBLIC

      CONTAINS

# if defined AFT_EIGENMODES
#  include "propagator_afte.h"
# elif defined FT_EIGENMODES
#  include "propagator_fte.h"
# elif defined OPT_PERTURBATION
#  include "propagator_op.h"
# elif defined FORCING_SV
#  include "propagator_fsv.h"
# elif defined STOCHASTIC_OPT
!!#  include "propagator_so.h"
# elif defined SO_SEMI
#  include "propagator_so_semi.h"
# elif defined SO_TRACE
!!#  include "propagator_so_trace.h"
# endif
#endif

      END MODULE propagator_mod
