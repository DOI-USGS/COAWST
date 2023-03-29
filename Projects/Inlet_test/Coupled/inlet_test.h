/*
** svn $Id: inlet_test.h 838 2008-11-17 04:22:18Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Inlet Test Case, waves-ocean (SWAN/ROMS) two-way coupling.
**
** Application flag:   INLET_TEST
** Input script:       ocean_inlet_test.in
**                     coupling_inlet_test.in
**                     sediment_inlet_test.in
*/

#define ROMS_MODEL
#define SWAN_MODEL
#define MCT_LIB

#define UV_VIS2
#define MIX_S_UV
#define MASKING
#define UV_ADV
#undef  UV_COR
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define OMEGA_IMPLICIT

#define SOLVE3D
#undef  WEC_MELLOR
#define WEC_VF
#define WDISS_WAVEMOD
#define UV_KIRBY
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_FSOBC
#define ANA_M2OBC
#undef  SALINITY
#undef  OUT_DOUBLE

/* define only one of the following 5 */
#undef  UV_LOGDRAG
#undef  UV_QDRAG
#undef  MB_BBL
#undef  SG_BBL
#define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# define SSW_LOGINT
/* define one of these 2 */
# define SSW_LOGINT_WBL
# undef  SSW_LOGINT_DIRECT
#endif

#ifdef SOLVE3D
# define GLS_MIXING
# ifdef GLS_MIXING
#  define RI_SPLINES
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
# define SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD_SOULSBY
#  undef  BEDLOAD_MPM
#  define BEDLOAD_VANDERA
#  ifdef BEDLOAD_VANDERA
/* select any or all of these 3 */
#   define BEDLOAD_VANDERA_ASYM_LIMITS
#   define BEDLOAD_VANDERA_SURFACE_WAVE
#   define BEDLOAD_VANDERA_WAVE_AVGD_STRESS
/* define one of these 2 */
#   define BEDLOAD_VANDERA_MADSEN_UDELTA
#   undef  BEDLOAD_VANDERA_DIRECT_UDELTA
#  endif
#  define SED_MORPH
#  undef  SED_SLUMP
#  undef  SLOPE_KIRWAN
#  undef  SLOPE_NEMETH
#  undef  SLOPE_LESSER
# endif
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define ANA_SEDIMENT
# endif
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
#endif

#define AVERAGES
#define STATIONS
#define DIAGNOSTICS_UV
