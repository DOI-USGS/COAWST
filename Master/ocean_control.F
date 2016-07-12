#include "cppdefs.h"
!!
!!svn $Id: ocean_control.F 795 2016-05-11 01:42:43Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!!                                                                     !
!!  Regional Ocean Model System (ROMS) Drivers:                        !
!!                                                                     !
!!  ad_ocean.h             Adjoint model driver                        !
!!  adsen_ocean.h          Adjoint sensitivity analysis driver         !
!!  afte_ocean.h           Adjoint finite time eigenmodes driver       !
!!  array_modes_w4dvar.h   W4DVAR array modes driver                   !
!!  correlation.h          Error covariance correlation driver         !
!!  fte_ocean.h            Finite time eigenmodes driver               !
!!  fsv_ocean.h            Forcing singular vectors driver             !
!!  hessian_fsv_ocean.h    Hessian forcing singular vectors driver     !
!!  hessian_op_ocean.h     Hessian singular vectors driver             !
!!  hessian_so_ocean.h     Hessian stochastic optimals driver          !
!!  is4dvar_ocean.h        Strong constraint, incremental 4DVar data   !
!!                           assimilation driver                       !
!!  nl_ocean.h             Nonlinear model driver (default)            !
!!  op_ocean.h             Optimal perturbations driver                !
!!  optobs_ocean.h         Optimal observations driver                 !
!!  obs_sen_is4dvar.h      Observations sensitivity driver to the      !
!!                           IS4DVAR data assimilation system          !
!!  obs_sen_w4dpsas.h      Observations sensitivity driver to the      !
!!                           W4DPSAS data assimilation system          !
!!  obs_sen_w4dvar.h       Observations sensitivity driver to the      !
!!                           W4DVAR  data assimilation system          !
!!  rp_ocean.h             Representer tangent linear model driver     !
!!  so_ocean.h             Stochastic optimals driver                  !
!!  so_semi_ocean.h        Stochastic optimals, semi-norm driver       !
!!  symmetry.h             Representer matrix, symmetry driver         !
!!  pert_ocean.h           Tangent linear and adjoint models sanity    !
!!                           test driver                               !
!!  picard_ocean.h         Picard test for representers tangent linear !
!!                           model driver                              !
!!  tlcheck_ocean.h        Tangent linear model linearization test     !
!!                           driver                                    !
!!  tl_ocean.h             Tangent linear model driver                 !
!!  tl_w4dpsas_ocean.h     Tangent linear driver to the W4DPSAS        !
!!                           data assimilation system                  !
!!  tl_w4dvar_ocean.h      Tangent linear driver to the W4DVAR         !
!!                           data assimilation system                  !
!!  w4dpsas_ocean.h        Weak constraint 4D-PSAS assimilation driver !
!!  w4dvar_ocean.h         Weak constraint 4DVAR assimilation,         !
!!                           indirect representer method               !
!!                                                                     !
!!======================================================================
!!
#ifdef ROMS_MODEL
# if defined AD_SENSITIVITY
#  include "adsen_ocean.h"
# elif defined AFT_EIGENMODES
#  include "afte_ocean.h"
# elif defined CORRELATION
#  include "correlation.h"
# elif defined FT_EIGENMODES
#  include "fte_ocean.h"
# elif defined FORCING_SV
#  include "fsv_ocean.h"
# elif defined IS4DVAR_SENSITIVITY
#  include "obs_sen_is4dvar.h"
# elif defined OPT_PERTURBATION
#  include "op_ocean.h"
# elif defined HESSIAN_SV
#  include "hessian_op_ocean.h"
# elif defined HESSIAN_FSV
#  include "hessian_fsv_ocean.h"
# elif defined OPT_OBSERVATIONS
#  include "optobs_ocean.h"
# elif defined STOCHASTIC_OPT
#  ifdef HESSIAN_SO
#    include "hessian_so_ocean.h"
#  else
#    include "so_ocean.h"
#  endif
# elif defined SO_SEMI
#  include "so_semi_ocean.h"
# elif defined TLM_CHECK
#  include "tlcheck_ocean.h"
# elif defined INNER_PRODUCT || defined SANITY_CHECK
#  include "pert_ocean.h"
# elif defined PICARD_TEST
#  include "picard_ocean.h"
# elif defined R_SYMMETRY
#  include "symmetry.h"
# elif defined IS4DVAR
#  include "is4dvar_ocean.h"
# elif defined W4DPSAS
#  include "w4dpsas_ocean.h"
# elif defined W4DVAR
#  include "w4dvar_ocean.h"
# elif defined W4DPSAS_SENSITIVITY
#  include "obs_sen_w4dpsas.h"
# elif defined TL_W4DPSAS
#  include "tl_w4dpsas_ocean.h"
# elif defined W4DVAR_SENSITIVITY
#  include "obs_sen_w4dvar.h"
# elif defined TL_W4DVAR
#  include "tl_w4dvar_ocean.h"
# elif defined ARRAY_MODES || defined CLIPPING
#  include "array_modes_w4dvar.h"
# else
#  if defined TLM_DRIVER
#   include "tl_ocean.h"
#  elif defined RPM_DRIVER
#   include "rp_ocean.h"
#  elif defined ADM_DRIVER
#   include "ad_ocean.h"
#  else
#   include "nl_ocean.h"
#  endif
# endif
#else
      MODULE ocean_control
      END MODULE ocean_control
#endif
