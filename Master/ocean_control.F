#include "cppdefs.h"
!!
!!svn $Id: ocean_control.F 1054 2021-03-06 19:47:12Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2021 The ROMS/TOMS Group       Andrew M. Moore   !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!!                                                                     !
!! Regional Ocean Model System (ROMS) Drivers:                         !
!!                                                                     !
!! ad_ocean.h                   Adjoint model driver                   !
!!                                                                     !
!! adsen_ocean.h                Adjoint sensitivity analysis driver    !
!!                                                                     !
!! afte_ocean.h                 Adjoint finite time eigenmodes driver  !
!!                                                                     !
!! array_modes.h                Dual 4D-Var array modes driver         !
!!                                                                     !
!! correlation.h                Error covariance correlation driver    !
!!                                                                     !
!! fte_ocean.h                  Finite time eigenmodes driver          !
!!                                                                     !
!! fsv_ocean.h                  Forcing singular vectors driver        !
!!                                                                     !
!! hessian_fsv_ocean.h          Hessian forcing singular vectors       !
!!                                driver                               !
!!                                                                     !
!! hessian_op_ocean.h           Hessian singular vectors driver        !
!!                                                                     !
!! hessian_so_ocean.h           Hessian stochastic optimals driver     !
!!                                                                     !
!! i4dvar_ocean.h               Strong constraint, incremental 4D-Var  !
!!                                data assimilation driver (I4D-Var)   !
!!                                                                     !
!! nl_ocean.h                   Nonlinear model driver (default)       !
!!                                                                     !
!! op_ocean.h                   Optimal perturbations driver           !
!!                                                                     !
!! optobs_ocean.h               Optimal observations driver            !
!!                                                                     !
!! obs_sen_i4dvar_analysis.h    Observations sensitivity driver to the !
!!                                I4D-Var data assimilation system     !
!!                                                                     !
!! obs_sen_rbl4dvar_analysis.h  Observations sensitivity driver to the !
!!                                RBL4D-Var data assimilation system   !
!!                                                                     !
!! obs_sen_rbl4dvar_forecast.h  Observations analysis and forecast     !
!!                                sensitivity driver for RBL4D-Var     !
!!                                                                     !
!! obs_sen_r4dvar_analysis.h    Observations sensitivity driver to the !
!!                                R4D-Var  data assimilation system    !
!!                                                                     !
!! pert_ocean.h                 Tangent linear and adjoint models      !
!!                                sanity test driver                   !
!!                                                                     !
!! picard_ocean.h               Picard test for representers tangent   !
!!                                linear model driver                  !
!!                                                                     !
!! rp_ocean.h                   Representer tangent linear model       !
!!                                driver                               !
!!                                                                     !
!! rbl4dvar_ocean.h             Strong/Weak constraint 4D-Var data     !
!!                                assimilation driver: Restricted,     !
!!                                B-preconditioned Lanczos method      !
!!                                (RBL4D-Var)                          !
!!                                                                     !
!! r4dvar_ocean.h               Strong/Weak constraint 4D-Var data     !
!!                                assimilation, indirect representer   !
!!                                method (R4D-Var)                     !
!!                                                                     !
!! so_ocean.h                   Stochastic optimals driver             !
!!                                                                     !
!! so_semi_ocean.h              Stochastic optimals, semi-norm driver  !
!!                                                                     !
!! sp4dvar_ocean.h              Saddle Point 4D-Var data assimilation  !
!!                                driver (SP4D-Var)                    !
!!                                                                     !
!! split_i4dvar_ocean.h         Multiple executables Strong constraint,!
!!                                incremental 4D-Var data assimilation !
!!                                driver (I4D-Var)                     !
!!                                                                     !
!! split_rbl4dvar_ocean.h       Multiple executables Strong/Weak       !
!!                                constraint 4D-Var data assimilation  !
!!                                driver: Restricted, B-preconditioned !
!!                                Lanczos method (RBL4D-Var)           !
!!                                                                     !
!! split_r4dvar_ocean.h         Multiple executables Strong/Weak       !
!!                                constraint 4D-Var data assimilation, !
!!                                indirect representer method (R4D-Var)!
!!                                                                     !
!! split_sp4dvar_ocean.h        Multiple executables Saddle Point      !
!!                                4D-Var data assimilation  driver     !
!!                                (SP4D-Var)                           !
!!                                                                     !
!! symmetry.h                   Representer matrix, symmetry driver    !
!!                                                                     !
!! tlcheck_ocean.h              Tangent linear model linearization     !
!!                                test driver                          !
!!                                                                     !
!! tl_ocean.h                   Tangent linear model driver            !
!!                                                                     !
!! tl_rbl4dvar_ocean.           Tangent linear driver to the RBL4D-Var !
!!                                data assimilation system             !
!!                                                                     !
!! tl_r4dvar_ocean.h            Tangent linear driver to the R4D-Var   !
!!                                data assimilation system             !
!!                                                                     !
!!======================================================================
!!
#ifdef ROMS_MODEL
# if defined AD_SENSITIVITY
#  include "adsen_ocean.h"

# elif defined AFT_EIGENMODES
#  include "afte_ocean.h"
# elif defined FT_EIGENMODES
#  include "fte_ocean.h"
# elif defined FORCING_SV
#  include "fsv_ocean.h"
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

# elif defined CORRELATION
#  include "correlation.h"
# elif defined INNER_PRODUCT || defined SANITY_CHECK
#  include "pert_ocean.h"
# elif defined PICARD_TEST
#  include "picard_ocean.h"
# elif defined R_SYMMETRY
#  include "symmetry.h"
# elif defined TLM_CHECK
#  include "tlcheck_ocean.h"

# elif defined ARRAY_MODES || defined CLIPPING
#  include "array_modes.h"
# elif defined I4DVAR && !defined SPLIT_I4DVAR
#   include "i4dvar_ocean.h"
# elif defined I4DVAR_ANA_SENSITIVITY
#  include "obs_sen_i4dvar_analysis.h"
# elif defined RBL4DVAR && !defined SPLIT_RBL4DVAR
#  include "rbl4dvar_ocean.h"
# elif defined RBL4DVAR_ANA_SENSITIVITY
#  include "obs_sen_rbl4dvar_analysis.h"
# elif defined RBL4DVAR_FCT_SENSITIVITY
#  include "obs_sen_rbl4dvar_forecast.h"
# elif defined R4DVAR && !defined SPLIT_R4DVAR
#  include "r4dvar_ocean.h"
# elif defined R4DVAR_ANA_SENSITIVITY
#  include "obs_sen_r4dvar_analysis.h"
# elif defined SPLIT_I4DVAR
#  include "split_i4dvar_ocean.h"
# elif defined SPLIT_RBL4DVAR
#  include "split_rbl4dvar_ocean.h"
# elif defined SPLIT_R4DVAR
#  include "split_r4dvar_ocean.h"
# elif defined SPLIT_SP4DVAR
#  include "split_sp4dvar_ocean.h"
# elif defined SP4DVAR && !defined SPLIT_SP4DVAR
#  include "sp4dvar_ocean.h"
# elif defined TL_RBL4DVAR
#  include "tl_rbl4dvar_ocean.h"
# elif defined TL_R4DVAR
#  include "tl_r4dvar_ocean.h"
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
