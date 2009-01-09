/*
** svn $Id: double_gyre.h 503 2008-01-10 00:11:51Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for 4DVar Data Assimilation Toy.
**
** Application flag:   DOUBLE_GYRE
** Input script:       ocean_double_gyre.in
*/

#define SOLVE3D

#undef  AD_SENSITIVITY          /* Adjoint Sensitivity */
#undef  AFT_EIGENMODES          /* Adjoint Finite Time Eigenmodes */
#undef  CORRELATION             /* Background-error Correlation Check */
#undef  GRADIENT_CHECK          /* TLM/ADM Gradient Check */
#undef  FORCING_SV              /* Forcing Singular Vectors */
#undef  FT_EIGENMODES           /* Finite Time Eigenmodes */
#undef  IS4DVAR_OLD             /* Old Incremental, strong constraint 4DVAR */
#define IS4DVAR                 /* Incremental, strong constraint 4DVAR */
#undef  NLM_DRIVER              /* Nonlinear Basic State trajectory */
#undef  OPT_PERTURBATION        /* Optimal perturbations */
#undef  PICARD_TEST             /* Picard Iterations Test */
#undef  R_SYMMETRY              /* Representer Matrix Symmetry Test */
#undef  S4DVAR                  /* Strong constraint 4DVAR */
#undef  SANITY_CHECK            /* Sanity Check */
#undef  SO_SEMI                 /* Stochastic Optimals: Semi-norm */
#undef  TLM_CHECK               /* Tangent Linear Model Check */
#undef  W4DPSAS                 /* Weak constraint 4D-PSAS */
#undef  W4DVAR                  /* Weak constraint 4DVAR */

/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state tracjectory.
**-----------------------------------------------------------------------------
*/

#if defined NLM_DRIVER
# if defined SOLVE3D                   /* 3D Application */
#  undef  UV_C2ADVECTION
#  define UV_ADV
#  define UV_COR
#  define UV_LDRAG
#  define UV_VIS2
#  define MIX_S_UV
#  undef  MIX_GEO_UV
#  undef  WJ_GRADP
#  define DJ_GRADPS
#  undef  TS_C2HADVECTION
#  undef  TS_C2VADVECTION
#  undef  TS_A4HADVECTION
#  undef  TS_A4VADVECTION
#  define TS_U3HADVECTION
#  define TS_C4VADVECTION
#  define TS_DIF2
#  undef  TS_DIF4
#  define MIX_S_TS
#  undef  MIX_GEO_TS
#  undef  MIX_ISO_TS
#  define NONLIN_EOS
#  define SALINITY
#  define AVERAGES
#  define SPLINES
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define TCLIMATOLOGY
#  define TCLM_NUDGING
#  define ANA_GRID
#  undef  MY25_MIXING
#  if defined MY25_MIXING
#   undef KANTHA_CLAYSON
#   undef CANUTO_A
#   undef N2S2_HORAVG
#  endif
#  undef  ANA_INITIAL
#  define ANA_TCLIMA
#  define ANA_SMFLUX
#  define ANA_STFLUX
#  define ANA_SSFLUX
#  define ANA_BSFLUX
#  define ANA_BTFLUX
#  undef  ANA_VMIX
#  undef  VERIFICATION
#  define FORWARD_RHS
#  define FORWARD_WRITE
#  define OUT_DOUBLE
# else                                 /* 2D Application */
#  undef  UV_C2ADVECTION
#  define UV_ADV
#  define UV_COR
#  define UV_LDRAG
#  define UV_VIS2
#  define AVERAGES
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define ANA_GRID
#  undef  ANA_INITIAL
#  define ANA_SMFLUX
#  undef  FORWARD_RHS
#  undef  FORWARD_WRITE
#  undef  OUT_DOUBLE
# endif
#endif

/*
**-----------------------------------------------------------------------------
**  Picard iteration test.
**-----------------------------------------------------------------------------
*/

#if defined PICARD_TEST
# if defined SOLVE3D                   /* 3D Application */
#  undef  UV_C2ADVECTION
#  define UV_ADV
#  define UV_COR
#  define UV_LDRAG
#  define UV_VIS2
#  define MIX_S_UV
#  undef  MIX_GEO_UV
#  undef  WJ_GRADP
#  define DJ_GRADPS
#  undef  TS_C2HADVECTION
#  undef  TS_C2VADVECTION
#  undef  TS_A4HADVECTION
#  undef  TS_A4VADVECTION
#  define TS_U3HADVECTION
#  define TS_C4VADVECTION
#  define TS_DIF2
#  undef  TS_DIF4
#  undef  MIX_S_TS
#  undef  MIX_GEO_TS
#  define MIX_ISO_TS
#  define NONLIN_EOS
#  define SPLINES
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  undef  TCLIMATOLOGY
#  undef  TCLM_NUDGING
#  define ANA_GRID
#  undef  ANA_INITIAL
#  undef  ANA_TCLIMA
#  define ANA_SMFLUX
#  define ANA_STFLUX
#  define ANA_BTFLUX
#  undef  ANA_VMIX
#  undef  IMPULSE
#  undef  FORWARD_RHS
#  define FORWARD_READ
#  define FORWARD_WRITE
#  undef  RST_SINGLE
#  define OUT_DOUBLE
# else                                 /* 2D Application */
#  undef  UV_C2ADVECTION
#  define UV_ADV
#  define UV_COR
#  define UV_LDRAG
#  define UV_VIS2
#  undef  UV_VIS4
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define ANA_GRID
#  undef  ANA_INITIAL
#  define ANA_SMFLUX
#  undef  FORWARD_RHS
#  define FORWARD_READ
#  define FORWARD_WRITE
#  define OUT_DOUBLE
# endif
#endif

/*
**-----------------------------------------------------------------------------
**  Generalized Stability Theory analysis.
**-----------------------------------------------------------------------------
*/

#if defined AFT_EIGENMODES || defined FT_EIGENMODES    || \
    defined FORCING_SV     || defined OPT_PERTURBATION || \
    defined SO_SEMI
# if defined SOLVE3D                   /* 3D Application */
#  define UV_ADV
#  define UV_LDRAG
#  define UV_VIS2
#  define MIX_S_UV
#  define UV_COR
#  define TS_U3HADVECTION
#  define TS_DIF2
#  define MIX_S_TS
#  define DJ_GRADPS
#  define NONLIN_EOS
#  define ANA_GRID
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define TCLIMATOLOGY
#  define TCLM_NUDGING
#  define ANA_TCLIMA
#  define ANA_GRID
#  define ANA_INITIAL
#  define ANA_SMFLUX
#  define ANA_STFLUX
#  define ANA_BTFLUX
#  define FORWARD_READ
# else                                 /* 2D Application */
#  define UV_ADV
#  define UV_VIS2
#  define UV_COR
#  define UV_LDRAG
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define ANA_GRID
#  define ANA_INITIAL
#  define ANA_SMFLUX
#  define FORWARD_READ
# endif
#endif

/*
**-----------------------------------------------------------------------------
**  Variational Data Assimilation.
**-----------------------------------------------------------------------------
*/

#if defined CORRELATION || defined GRADIENT_CHECK || \
    defined IS4DVAR     || defined IS4DVAR_OLD    || \
    defined R_SYMMETRY  || defined TLM_CHECK      || \
    defined W4DPSAS     || defined W4DVAR
# if defined SOLVE3D                   /* 3D Application */
#  undef  UV_C2ADVECTION
#  undef  UV_C4ADVECTION
#  undef  UV_SADVECTION
#  define UV_ADV
#  define UV_COR
#  define UV_LDRAG
#  define UV_VIS2
#  define MIX_S_UV
#  undef  MIX_GEO_UV
#  undef  WJ_GRADP
#  define DJ_GRADPS
#  undef  TS_FIXED
#  undef  TS_C2HADVECTION
#  undef  TS_C2VADVECTION
#  define TS_U3HADVECTION
#  define TS_C4VADVECTION
#  undef  TS_A4HADVECTION
#  undef  TS_A4VADVECTION
#  define TS_DIF2
#  define MIX_S_TS
#  undef  MIX_GEO_TS
#  undef  NONLIN_EOS
#  define SPLINES
#  define SALINITY
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define TCLIMATOLOGY
#  define TCLM_NUDGING
#  undef  ANA_INITIAL
#  define ANA_GRID
#  define ANA_TCLIMA
#  define ANA_SMFLUX
#  define ANA_STFLUX
#  define ANA_SSFLUX
#  define ANA_BTFLUX
#  define ANA_BSFLUX
#  define FORWARD_WRITE
#  define FORWARD_READ
#  undef  FULL_GRID
#  define OUT_DOUBLE
#  if defined W4DPSAS || defined W4DVAR
#    define CONVOLVE
#  endif
#  define VCONVOLUTION
#  define IMPLICIT_VCONV
#  ifdef S4DVAR
#   define ENERGY1_NORM
#   undef  ENERGY2_NORM
#   undef  ENERGY3_NORM
#   define N2NORM_PROFILE
#  endif
#  ifdef W4DVAR
#   define RPM_RELAXATION
#   undef  CONVOLVE
#  endif
#  ifdef IS4DVAR_OLD
#   undef  MULTIPLE_TLM
#  endif
# else                                 /* 2D Application */
#  undef  UV_C2ADVECTION
#  define UV_ADV
#  define UV_LDRAG
#  define UV_VIS2
#  undef  UV_VIS4
#  define UV_COR
#  undef  AVERAGES
#  define ANA_GRID
#  undef  ANA_INITIAL
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define ANA_SMFLUX
#  undef  ENERGY1_NORM
#  undef  ENERGY2_NORM
#  undef  ENERGY3_NORM
#  define FORWARD_WRITE
#  define FORWARD_READ
#  define OUT_DOUBLE
# endif
#endif

/*
**-----------------------------------------------------------------------------
**  Sanity check.
**-----------------------------------------------------------------------------
*/

#if defined SANITY_CHECK
# if defined SOLVE3D                   /* 3D Application */
#  undef  UV_C2ADVECTION
#  define UV_C4ADVECTION
#  undef  UV_SADVECTION
#  define UV_ADV
#  define UV_COR
#  define UV_LDRAG
#  undef  UV_QDRAG
#  undef  UV_VIS2
#  define UV_VIS4
#  undef  MIX_S_UV
#  define MIX_GEO_UV
#  define WJ_GRADP
#  undef  DJ_GRADPS
#  define TS_C2HADVECTION
#  define TS_C2VADVECTION
#  undef  TS_C4HADVECTION
#  undef  TS_C4VADVECTION
#  undef  TS_U3HADVECTION
#  undef  TS_C4VADVECTION
#  undef  TS_A4HADVECTION
#  undef  TS_A4VADVECTION
#  undef  TS_SVADVECTION
#  define TS_DIF2
#  undef  TS_DIF4
#  define MIX_S_TS
#  undef  MIX_GEO_TS
#  undef  MIX_ISO_TS
#  define NONLIN_EOS
#  undef  SPLINES
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define TCLIMATOLOGY
#  define TCLM_NUDGING
#  define ANA_INITIAL
#  define ANA_GRID
#  define ANA_TCLIMA
#  define ANA_SMFLUX
#  define ANA_STFLUX
#  define ANA_BTFLUX
#  define ANA_PERTURB
#  define SANITY_CHECK
#  define FORWARD_READ
#  undef  OUT_DOUBLE
# else                                 /* 2D Application */
#  define UV_ADV
#  define UV_VIS2
#  undef  UV_VIS4
#  define UV_COR
#  define UV_LDRAG
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define ANA_GRID
#  define ANA_INITIAL
#  define ANA_SMFLUX
#  define ANA_PERTURB
#  define FORWARD_READ
# endif
#endif

/*
**-----------------------------------------------------------------------------
**  3D Double-Gyre Adjoint sensitivity.
**-----------------------------------------------------------------------------
*/

#if defined AD_SENSITIVITY
# if defined SOLVE3D                   /* 3D Application */
#  undef  UV_C2ADVECTION
#  undef  UV_C4ADVECTION
#  undef  UV_SADVECTION
#  define UV_ADV
#  define UV_COR
#  define UV_LDRAG
#  define UV_VIS2
#  define MIX_S_UV
#  undef  MIX_GEO_UV
#  undef  WJ_GRADP
#  undef  DJ_GRADPS
#  define TS_U3HADVECTION
#  define TS_C4VADVECTION
#  undef  TS_A4HADVECTION
#  undef  TS_A4VADVECTION
#  define TS_DIF2
#  define MIX_S_TS
#  undef  MIX_GEO_TS
#  define NONLIN_EOS
#  undef  SPLINES
#  define AVERAGES
#  define AVERAGES_FLUXES
#  define AVERAGES_AKV
#  define AVERAGES_AKT
#  define AVERAGES_QUADRATIC
#  define WESTERN_WALL
#  define EASTERN_WALL
#  define NORTHERN_WALL
#  define SOUTHERN_WALL
#  define TCLIMATOLOGY
#  define TCLM_NUDGING
#  undef  ANA_INITIAL
#  define ANA_GRID
#  define ANA_TCLIMA
#  define ANA_SMFLUX
#  define ANA_STFLUX
#  define ANA_BTFLUX
#  define ANA_SCOPE
#  define AD_SENSITIVITY
#  define FORWARD_READ
#  define OUT_DOUBLE
# endif
#endif
