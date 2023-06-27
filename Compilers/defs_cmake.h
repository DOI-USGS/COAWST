/*
** Include file "defs_cmake.h"
**
** git $Id$
** svn $Id: defs_cmake.h 1151 2023-02-09 03:08:53Z arango $
********************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2023 The ROMS/TOMS Group              David Robertson  **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**                                                                           **
** This file is used by CMake to determine whether additional models         **
** (adjoint, tangent, an representer) should be enabled. It also checks      **
** if the ARPACK library is needed.                                          **
**                                                                           **
*******************************************************************************
*/

#if defined ROMS_HEADER
# include ROMS_HEADER
#else
   CPPDEFS - Choose an appropriate ROMS application.
#endif

/*
** Check for 4D-Var deprecated CPP options.
*/

#ifdef IS4DVAR
# ifndef I4DVAR
#  define I4DVAR
# endif
#endif

#ifdef IS4DVAR_SENSITIVITY
# ifndef I4DVAR_ANA_SENSITIVITY
#  define I4DVAR_ANA_SENSITIVITY
# endif
#endif

#ifdef W4DPSAS
# ifndef RBL4DVAR
#  define RBL4DVAR
# endif
#endif

#ifdef W4DPSAS_SENSITIVITY
# ifndef RBL4DVAR_ANA_SENSITIVITY
#  define RBL4DVAR_ANA_SENSITIVITY
# endif
#endif

#ifdef W4DPSAS_FCT_SENSITIVITY
# ifndef RBL4DVAR_FCT_SENSITIVITY
#  define RBL4DVAR_FCT_SENSITIVITY
# endif
#endif

#ifdef W4DVAR
# ifndef R4DVAR
#  define R4DVAR
# endif
#endif

#ifdef W4DVAR_SENSITIVITY
# ifndef R4DVAR_ANA_SENSITIVITY
#  define R4DVAR_ANA_SENSITIVITY
# endif
#endif

/*
** Set multiple excutables split 4D-Var.
*/

#if defined SPLIT_I4DVAR   || \
    defined SPLIT_RBL4DVAR || \
    defined SPLIT_R4DVAR   || \
    defined SPLIT_SP4DVAR
# define SPLIT_4DVAR
#endif

/*
** If split 4D-Var, activate the unsplit option since both share
** identical configuration to avoid too many directives changes.
*/

#if !defined I4DVAR && defined SPLIT_I4DVAR
# define I4DVAR
#endif

#if !defined RBL4DVAR && defined SPLIT_RBL4DVAR
# define RBL4DVAR
#endif

#if !defined R4DVAR && defined SPLIT_R4DVAR
# define R4DVAR
#endif

#if !defined SP4DVAR && defined SPLIT_SP4DVAR
# define SP4DVAR
#endif


/*
** Set 4D-Var sensitivity switch.
*/

#if defined RBL4DVAR_ANA_SENSITIVITY || \
    defined RBL4DVAR_FCT_SENSITIVITY || \
    defined R4DVAR_ANA_SENSITIVITY
# define SENSITIVITY_4DVAR
#endif

/*
** Set perturbation tangent linear, fine amplitude tangent linear,
** and adjoint model switches.
*/

#if defined ARRAY_MODES            || \
    defined CLIPPING               || \
    defined CORRELATION            || \
    defined FORCING_SV             || \
    defined FT_EIGENMODES          || \
    defined HESSIAN_FSV            || \
    defined HESSIAN_SO             || \
    defined HESSIAN_SV             || \
    defined INNER_PRODUCT          || \
    defined I4DVAR                 || \
    defined I4DVAR_ANA_SENSITIVITY || \
    defined JEDI                   || \
    defined OPT_PERTURBATION       || \
    defined OPT_OBSERVATIONS       || \
    defined PICARD_TEST            || \
    defined RBL4DVAR               || \
    defined RPM_DRIVER             || \
    defined R4DVAR                 || \
    defined R_SYMMETRY             || \
    defined SANITY_CHECK           || \
    defined SENSITIVITY_4DVAR      || \
    defined SPLIT_I4DVAR           || \
    defined SPLIT_RBL4DVAR         || \
    defined SPLIT_R4DVAR           || \
    defined SPLIT_SP4DVAR          || \
    defined SP4DVAR                || \
    defined STOCHASTIC_OPT         || \
    defined TLM_CHECK              || \
    defined TLM_DRIVER             || \
    defined TL_RBL4DVAR            || \
    defined TL_R4DVAR
TANGENT
#define FOUND
#endif

#if defined AD_SENSITIVITY         || \
    defined ADM_DRIVER             || \
    defined AFT_EIGENMODES         || \
    defined ARRAY_MODES            || \
    defined CLIPPING               || \
    defined CORRELATION            || \
    defined FORCING_SV             || \
    defined HESSIAN_SO             || \
    defined HESSIAN_FSV            || \
    defined HESSIAN_SV             || \
    defined INNER_PRODUCT          || \
    defined I4DVAR                 || \
    defined I4DVAR_ANA_SENSITIVITY || \
    defined JEDI                   || \
    defined OPT_PERTURBATION       || \
    defined OPT_OBSERVATIONS       || \
    defined RBL4DVAR               || \
    defined R4DVAR                 || \
    defined R_SYMMETRY             || \
    defined SANITY_CHECK           || \
    defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI                || \
    defined SPLIT_I4DVAR           || \
    defined SPLIT_RBL4DVAR         || \
    defined SPLIT_R4DVAR           || \
    defined SPLIT_SP4DVAR          || \
    defined SP4DVAR                || \
    defined STOCHASTIC_OPT         || \
    defined TLM_CHECK              || \
    defined TL_RBL4DVAR            || \
    defined TL_R4DVAR
ADJOINT
#define FOUND
#endif

#if defined ARRAY_MODES            || \
    defined CLIPPING               || \
    defined PICARD_TEST            || \
    defined RPM_DRIVER             || \
    defined TL_R4DVAR              || \
    defined R4DVAR                 || \
    defined R4DVAR_ANA_SENSITIVITY
REPRESENTER
#define FOUND
#endif

#if defined ARRAY_MODES              || \
    defined CLIPPING                 || \
    defined I4DVAR                   || \
    defined PROPAGATOR               || \
    defined RBL4DVAR                 || \
    defined RBL4DVAR_ANA_SENSITIVITY || \
    defined R4DVAR                   || \
    defined R4DVAR_ANA_SENSITIVITY   || \
    defined SP4DVAR                  || \
    defined TL_RBL4DVAR              || \
    defined TL_R4DVAR
ARPACK
#define FOUND
#endif

/*
** Determine if using the NCAR parallel-IO (SCORPIO) library is needed.
*/

#ifdef PIO_LIB
SCORPIO
#define FOUND
#endif

#ifndef FOUND
NONE
#endif

