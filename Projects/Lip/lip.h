/*
** svn $Id: shoreface.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Shore Face Planar Beach Test Case.
**
** Application flag:   SHOREFACE
** Input scripts:      ocean_shoreface.h
**                     sediment_shoreface.h
*/

#define ROMS_MODEL
#define INWAVE_MODEL

#ifdef  INWAVE_MODEL
# undef  INWAVE_SWAN_COUPLING
# define ACX_ADVECTION
# define ACY_ADVECTION
# undef  ACT_ADVECTION
# undef  THETA_AC_PERIODIC
# define DOPPLER
# undef  WDISS_GAMMA
# define WDISS_ROELVINK
/* # define WDISS_INWAVE */
# define RAMP_INWAVE
# define WEC_VF
# define UV_KIRBY
# define ROLLER_RENIERS /*  this is a problem */
# undef  ROLLER_SVENDSEN
#endif

#define UV_VIS2
#define MIX_S_UV
#undef  DIAGNOSTICS_UV
#undef  AVERAGES
#define WET_DRY
#undef  OUT_DOUBLE
#define UV_ADV
#define DJ_GRADPS
#undef  SALINITY
#define SOLVE3D
#undef  SPLINES_VDIFF
#undef  SPLINES_VVISC
#undef  TS_FIXED

#define MASKING
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX

#ifdef SOLVE3D
# define SSW_BBL
# ifdef SSW_BBL
#  define SSW_CALC_ZNOT
#  define SSW_LOGINT
/* define one of these 2 */
#  define SSW_LOGINT_WBL
#  undef  SSW_LOGINT_DIRECT
# else
#  define UV_LOGDRAG
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
# undef  ANA_VMIX
# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define RI_SPLINES
#  undef CRAIG_BANNER
#  undef CHARNOK
#  undef ZOS_HSIG
#  undef TKE_WAVEDISS
# endif
#endif
