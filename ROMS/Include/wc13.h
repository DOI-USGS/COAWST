/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for the California Current System, 1/3 degree resolution.
**
** Application flag:   WC13
** Input script:       ocean_wc13.in
**                     s4dvar.in
**
** Available Drivers options: choose only one and activate it in the
**                            build.sh script (MY_CPP_FLAGS definition)
**
**  AD_SENSITIVITY            Adjoint Sensitivity Driver
**  AFT_EIGENMODES            Adjoint Finite Time Eigenmodes
**  ARRAY_MODES               Stabilized representer matrix array modes
**  CLIPPING                  Stabilized representer matrix clipped analysis
**  CORRELATION               Background-error Correlation Check
**  GRADIENT_CHECK            TLM/ADM Gradient Check
**  FORCING_SV                Forcing Singular Vectors
**  FT_EIGENMODES             Finite Time Eigenmodes
**  IS4DVAR                   Incremental, strong constraint 4DVAR
**  NLM_DRIVER                Nonlinear Basic State trajectory
**  OPT_PERTURBATION          Optimal perturbations
**  PICARD_TEST               Picard Iterations Test
**  R_SYMMETRY                Representer Matrix Symmetry Test
**  SANITY_CHECK              Sanity Check
**  SO_SEMI                   Stochastic Optimals: Semi-norm
**  TLM_CHECK                 Tangent Linear Model Check
**  W4DPSAS                   Weak constraint 4D-PSAS
**  W4DVAR                    Weak constraint 4DVAR
**  VERIFICATION              NL Observation Verification Driver
**  NORMALIZATION             Background error Covariance Normalization
*/

/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state settings.
**-----------------------------------------------------------------------------
*/

#ifdef VERIFICATION
# define FULL_GRID
#endif

#define ANA_BSFLUX
#define ANA_BTFLUX

#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define MIX_GEO_TS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_DIF2
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define PROFILE
#define SPHERICAL
#define MASKING

#ifdef NLM_DRIVER
# define AVERAGES           /* define if writing out time-averaged data */
#endif

/*
**  Vertical Mixing parameterization
*/

#define GLS_MIXING
#ifdef GLS_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
# define RI_SPLINES
#endif

/*
** Surface atmospheric fluxes. Note, that we must define DIURNAL_SRFLUX
** when using daily averaged fields.
*/

#define BULK_FLUXES      /* turn ON or OFF bulk fluxes computation */

#define DIURNAL_SRFLUX   /* impose shortwave radiation local diurnal cycle */
#define SOLAR_SOURCE     /* define solar radiation source term */
#define LONGWAVE_OUT     /* Compute net longwave radiation internally */
#define EMINUSP          /* turn ON internal calculation of E-P */

/*
**-----------------------------------------------------------------------------
**  Variational Data Assimilation.
**-----------------------------------------------------------------------------
*/

/*
**  Options to compute error covariance normalization coefficients.
*/

#ifdef NORMALIZATION
# define ADJUST_BOUNDARY
# define ADJUST_WSTRESS
# define ADJUST_STFLUX
# define CORRELATION
# define VCONVOLUTION
# define IMPLICIT_VCONV
# define FULL_GRID
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

/*
**  Options for adjoint-based algorithms sanity checks.
*/

#ifdef SANITY_CHECK
# define FULL_GRID
# define FORWARD_READ
# define FORWARD_WRITE
# define FORWARD_MIXING
# define OUT_DOUBLE
# define ANA_PERTURB
# define ANA_INITIAL
#endif


/*
**  Common options to all 4DVAR algorithms.
*/

#if defined ARRAY_MODES || defined CLIPPING            || \
    defined IS4DVAR     || defined IS4DVAR_SENSITIVITY || \
    defined W4DPSAS     || defined W4DPSAS_SENSITIVITY || \
    defined W4DVAR      || defined W4DVAR_SENSITIVITY
# define ADJUST_BOUNDARY
# define ADJUST_WSTRESS
# define ADJUST_STFLUX
# define NL_BULK_FLUXES
# define VCONVOLUTION
# define IMPLICIT_VCONV
# define BALANCE_OPERATOR
# ifdef BALANCE_OPERATOR
#  define ZETA_ELLIPTIC
# endif
# define FORWARD_WRITE
# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

/*
**  Special options for each 4DVAR algorithm.
*/

#if defined ARRAY_MODES || \
    defined W4DVAR      || defined W4DVAR_SENSITIVITY
# define RPM_RELAXATION
#endif
