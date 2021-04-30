/*
** svn $Id: sw06_coarse.h 1054 2021-03-06 19:47:12Z arango $
*******************************************************************************
** Copyright (c) 2002-2021 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Shallow Water 2006 Experiment, coarse grid resolution.
**
** Application flag:   SW06_COARSE
** Input script:       ocean_sw05_coarse.in
**                     s4dvar.in
*/

#define ROMS_MODEL

#undef  AFT_EIGENMODES          /* Adjoint Finite Time Eigenmodes */
#undef  CORRELATION             /* Background-error Correlation Check */
#undef  FORCING_SV              /* Forcing Singular Vectors */
#undef  FT_EIGENMODES           /* Finite Time Eigenmodes */
#define I4DVAR                  /* Incremental, strong constraint 4D-Var */
#undef  NLM_DRIVER              /* Nonlinear Basic State trajectory */
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
**  Nonlinear model or basic state tracjectory.
**-----------------------------------------------------------------------------
*/

#if defined NLM_DRIVER
# define UV_ADV
# define UV_COR
# define UV_QDRAG
# define UV_VIS2
# define MIX_S_UV
# define TS_DIF2
# define MIX_GEO_TS
# define DJ_GRADPS
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define TS_U3HADVECTION
# define SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# undef  AVERAGES
# define MASKING
# define SOLAR_SOURCE

# define GLS_MIXING
# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  undef  CANUTO_A
#  define N2S2_HORAVG
#  define RI_SPLINES
# endif

# define BULK_FLUXES
# ifdef BULK_FLUXES
#  define LONGWAVE_OUT
#  define ANA_CLOUD
#  define ANA_RAIN
# endif

# undef  RAMP_TIDES
# define SSH_TIDES
# undef  UV_TIDES
# ifdef SSH_TIDES
#  define FSOBC_REDUCED
#  define ADD_FSOBC
# endif
# if defined UV_TIDES || defined SSH_TIDES
#  define ADD_M2OBC
# endif

# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX

# define OUT_DOUBLE

#else

/*
**-----------------------------------------------------------------------------
**  Adjoint-based drivers.
**-----------------------------------------------------------------------------
*/

# define UV_ADV
# define UV_COR
# define UV_QDRAG
# define UV_VIS2
# define MIX_S_UV
# define TS_DIF2
# define MIX_GEO_TS
# define DJ_GRADPS
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define TS_U3HADVECTION
# define SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# undef  AVERAGES
# define MASKING
# define SOLAR_SOURCE

# define GLS_MIXING
# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  undef  CANUTO_A
#  define N2S2_HORAVG
#  define RI_SPLINES
# endif

# undef  RAMP_TIDES
# define SSH_TIDES
# undef  UV_TIDES
# ifdef SSH_TIDES
#  define FSOBC_REDUCED
#  define ADD_FSOBC
# if defined UV_TIDES || defined SSH_TIDES
#  define ADD_M2OBC
# endif

# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX

# define VCONVOLUTION
# define IMPLICIT_VCONV
# undef  FULL_GRID

# define FORWARD_READ
# define FORWARD_WRITE
# define FORWARD_MIXING
# define OUT_DOUBLE

#endif
