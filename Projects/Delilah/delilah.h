/*
** svn $Id: inlet_test.h 838 2008-11-17 04:22:18Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Inlet Test Case, waves-ocean (INWAVE/ROMS) coupling.
**
** Application flag:   INLET_TEST
** Input script:       ocean_inlet_test.in
**                     sediment_inlet_test.in
*/

#define ROMS_MODEL
#define INWAVE_MODEL

#ifdef  INWAVE_MODEL
# define INWAVE_SWAN_COUPLING
# define ACX_ADVECTION
# define ACY_ADVECTION
# define ACT_ADVECTION
# undef  DOPPLER
# undef  WDISS_GAMMA
# define WDISS_ROELVINK
# define VARY_ACBC
# define ROLLER_RENIERS
# define WDISS_FRIC
# define RAMP_INWAVE
#endif

#define WEC_VF
#define WET_DRY
#define CURVGRID

#define UV_VIS2
#define MIX_S_UV
#define MASKING
#define UV_ADV
#undef  UV_COR
#define DJ_GRADPS
#define SOLVE3D
#undef  SPLINES_VVISC
#undef  SPLINES_VDIFF
#undef  SALINITY

#define ANA_INITIAL
#undef  ANA_SMFLUX
#undef  ANA_FSOBC
#define ANA_M2OBC

/* define only one of the following 5 */
#undef  UV_LOGDRAG
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
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define RI_SPLINES
#  define TKE_WAVEDISS
#  define ZOS_HSIG
# endif
# undef SEDIMENT
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
#undef DIAGNOSTICS_UV
#define STATIONS
