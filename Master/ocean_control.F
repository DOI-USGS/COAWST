#include "cppdefs.h"
!!
!!svn $Id: ocean_control.F 645 2008-07-16 19:26:48Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!====================================================================== 
!!                                                                     !
!!  Regional Ocean Model System (ROMS) Drivers:                        !
!!                                                                     !
!!  convolution.h     ROMS/TOMS convolve an adjoint solution.          !
!!  correlation.h     ROMS/TOMS background-error correlation.          !
!!  nl_ocean.h        ROMS/TOMS nonlinear model.                       !
!!  tl_ocean.h        ROMS/TOMS tangent linear model.                  !
!!  rp_ocean.h        ROMS/TOMS representer tangent linear model.      !
!!  ad_ocean.h        ROMS/TOMS adjoint model.                         !
!!  adsen_ocean.h     ROMS/TOMS adjoint sensitivity analysis driver.   !
!!  afte_ocean.h      ROMS/TOMS adjoint finite time eigenmodes driver. !
!!  fte_ocean.h       ROMS/TOMS finite time eigenmodes driver.         !
!!  fsv_ocean.h       ROMS/TOMS forcing singular vectors driver.       !
!!  grad_ocean.h      ROMS/TOMS tangent linear and adjoint models      !
!!                      gradient test.                                 !
!!  tlcheck_ocean.h   ROMS/TOMS tangent linear model linearization     !
!!                      test.                                          !
!!  pert_ocean.h      ROMS/TOMS tangent linear and adjoint models      !
!!                      sanity test.                                   !
!!  picard_ocean.h    ROMS/TOMS Picard test for representers tangent   !
!!                      linear model.                                  ! 
!!  obs_sen_ocean.h   ROMS/TOMS 4DVAR observations sensitivity.        !
!!  op_ocean.h        ROMS/TOMS optimal perturbations.                 !
!!  optobs_ocean.h    ROMS/TOMS optimal observations.                  !
!!  symmetry.h        ROMS/TOMS representer matrix symmetry driver.    !
!!  s4dvar_ocean.h    ROMS/TOMS strong constraint 4DVAR assimilation.  !
!!  so_semi_ocean.h   ROMS/TOMS stochastic optimals, semi-norm.        !
!!  is4dvar_*.h       ROMS/TOMS Strong constraint, incremental 4DVAR   !
!!                      assimilation with various conjugate gradient   !
!!                      algorithms (cgradient.h, cgradient_lanczos.h,  !
!!                      descent.F).                                    !
!!  w4dpsas_ocean.h   ROMS/TOMS weak constraint 4D-PSAS assimilation.  !
!!  w4dvar_ocean.h    ROMS/TOMS weak constraint 4DVAR assimilation.    !
!!                                                                     !
!!======================================================================
!!
#if defined AD_SENSITIVITY
# include "adsen_ocean.h"
#elif defined AFT_EIGENMODES
# include "afte_ocean.h"
#elif defined CONVOLUTION
# include "convolution.h"
#elif defined CORRELATION
# include "correlation.h"
#elif defined FT_EIGENMODES
# include "fte_ocean.h"
#elif defined FORCING_SV
# include "fsv_ocean.h"
#elif defined GRADIENT_CHECK
# include "grad_ocean.h"
#elif defined OBS_SENSITIVITY
# include "obs_sen_ocean.h"
#elif defined OPT_PERTURBATION
# include "op_ocean.h"
#elif defined OPT_OBSERVATIONS
# include "optobs_ocean.h"
#elif defined SO_SEMI
# include "so_semi_ocean.h"
#elif defined TLM_CHECK
# include "tlcheck_ocean.h"
#elif defined INNER_PRODUCT || defined SANITY_CHECK
# include "pert_ocean.h"
#elif defined PICARD_TEST
# include "picard_ocean.h"
#elif defined R_SYMMETRY
# include "symmetry.h"
#elif defined S4DVAR
# include "s4dvar_ocean.h"
#elif defined IS4DVAR
# if defined LANCZOS
#  include "is4dvar_lanczos_ocean.h"
# else
#  include "is4dvar_ocean.h"
# endif
#elif defined IS4DVAR_OLD
# include "is4dvar_ocean_old.h"
#elif defined W4DPSAS
# include "w4dpsas_ocean.h"
#elif defined W4DVAR
# include "w4dvar_ocean.h"
#else
# if defined TLM_DRIVER
#  include "tl_ocean.h"
# elif defined RPM_DRIVER
#  include "rp_ocean.h"
# elif defined ADM_DRIVER
#  include "ad_ocean.h"
# else
#  include "nl_ocean.h"
# endif
#endif
