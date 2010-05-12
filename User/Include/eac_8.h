/*
** svn $Id: eac_8.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for East Australia Current, 1/8 degree resolution
**
** Application flag:   EAC_8
** Input script:       ocean_eac_8.in
*/


#undef  AFT_EIGENMODES          /* Adjoint Finite Time Eigenmodes */
#undef  CORRELATION             /* Background-error Correlation Check */
#undef  FORCING_SV              /* Forcing Singular Vectors */
#undef  FT_EIGENMODES           /* Finite Time Eigenmodes */
#undef  IS4DVAR_OLD             /* Old Incremental, strong constraint 4DVAR */
#define IS4DVAR                 /* Incremental, strong constraint 4DVAR */
#undef  NLM_DRIVER              /* Nonlinear Basic State trajectory */
#undef  OPT_PERTURBATION        /* Optimal perturbations */
#undef  PICARD_TEST             /* Picard Iterations Test */
#undef  R_SYMMETRY              /* Representer Matrix Symmetry Test */
#undef  SANITY_CHECK            /* Sanity Check */
#undef  SO_SEMI                 /* Stochastic Optimals: Semi-norm */
#undef  TLM_CHECK               /* Tangent Linear Model Check */
#undef  W4DPSAS                 /* Weak constraint 4D-PSAS */
#undef  W4DVAR                  /* Weak constraint 4DVAR */

#define UV_ADV
#define UV_COR
#ifdef NLM_DRIVER
# define UV_QDRAG
#else
# define UV_LDRAG
#endif
#define DJ_GRADPS
#define UV_VIS2
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  TS_DIF2
#undef  MIX_GEO_TS
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPLINES
#define MASKING
#ifdef NLM_DRIVER
# define SOLAR_SOURCE
#endif

#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3CLAMPED
#define EAST_TCLAMPED
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3CLAMPED
#define WEST_TCLAMPED
#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_M3CLAMPED
#define NORTH_TCLAMPED
#define SOUTH_FSCHAPMAN
#define SOUTH_M2FLATHER
#define SOUTH_M3CLAMPED
#define SOUTH_TCLAMPED
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX

#undef  BULK_FLUXES
#ifdef BULK_FLUXES
# undef LONGWAVE
# ifdef LONGWAVE
#  define ANA_CLOUD
# endif
# undef  SRELAXATION
# define ANA_RAIN
# define ANA_SSFLUX
#endif

#ifdef NLM_DRIVER
# define LMD_MIXING
#endif
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif

#ifndef NLM_DRIVER
# undef  FULL_GRID
# define VCONVOLUTION
# define IMPLICIT_VCONV
#endif
#define FORWARD_WRITE
#define FORWARD_MIXING
#ifndef NLM_DRIVER
# define FORWARD_READ
#endif
#define OUT_DOUBLE
