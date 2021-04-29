/*
** svn $Id: scb.h 1054 2021-03-06 19:47:12Z arango $
*******************************************************************************
** Copyright (c) 2002-2021 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Southern California Bight.
**
** Application flag:   SCB
** Input script:       ocean_scb.in
*/


#define ROMS_MODEL

#undef  AFT_EIGENMODES          /* Adjoint Finite Time Eigenmodes */
#undef  CORRELATION             /* Background-error Correlation Check */
#undef  FORCING_SV              /* Forcing Singular Vectors */
#undef  FT_EIGENMODES           /* Finite Time Eigenmodes */
#undef  I4DVAR                  /* Incremental, strong constraint 4D-Var */
#define NLM_DRIVER              /* Nonlinear Basic State trajectory */
#undef  OPT_PERTURBATION        /* Optimal perturbations */
#undef  PICARD_TEST             /* Picard Iterations Test */
#undef  RBL4DVAR                /* Strong/Weak constraint RBL4D-Var */
#undef  R4DVAR                  /* Strong/Weak constraint R4D-Var */
#undef  R_SYMMETRY              /* Representer Matrix Symmetry Test */
#undef  SANITY_CHECK            /* Sanity Check */
#undef  SO_SEMI                 /* Stochastic Optimals: Semi-norm */
#undef  TLM_CHECK               /* Tangent Linear Model Check */

/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state tracjectory.
**-----------------------------------------------------------------------------
*/

#if defined NLM_DRIVER
# define UV_ADV
# define UV_COR
# define UV_LDRAG
# define UV_VIS2
# define MIX_S_UV
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define DJ_GRADPS
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# define TS_DIF2
# define MIX_S_TS
# undef  MIX_GEO_TS
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# define MASKING
# define SOLVE3D
# define AVERAGES
# define AVERAGES_QUADRATIC
# define DIAGNOSTICS_UV
# define DIAGNOSTICS_TS

# define LMD_MIXING
# ifdef LMD_MIXING
#  define SOLAR_SOURCE
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_NONLOCAL
#  define LMD_SKPP
#  define RI_SPLINES
# endif

# define ANA_BSFLUX
# define ANA_BTFLUX

# define FORWARD_MIXING
# define FORWARD_READ
# define FORWARD_WRITE
# define OUT_DOUBLE

#else

/*
**-----------------------------------------------------------------------------
**  Adjoint-based drivers.
**-----------------------------------------------------------------------------
*/

# define UV_ADV
# define UV_COR
# define UV_LDRAG
# define UV_VIS2
# define MIX_S_UV
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define DJ_GRADPS
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# define TS_DIF2
# define MIX_S_TS
# undef  MIX_GEO_TS
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# define MASKING
# define SOLVE3D

# define LMD_MIXING
# ifdef LMD_MIXING
#  define SOLAR_SOURCE
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_NONLOCAL
#  define LMD_SKPP
#  define RI_SPLINES
# endif

# define ANA_BSFLUX
# define ANA_BTFLUX

# if defined RBL4DVAR || defined R4DVAR
#  define CONVOLVE
# endif

# define VCONVOLUTION
# define IMPLICIT_VCONV
# undef  FULL_GRID

# define FORWARD_MIXING
# define FORWARD_READ
# define FORWARD_WRITE
# define OUT_DOUBLE
#endif
