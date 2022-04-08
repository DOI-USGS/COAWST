/*
** svn $Id: shoreface.h 139 2008-01-10 00:17:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Kumar, N., Voulgaris, G., Warner, J.C., and M., Olabarrieta (2012).
** Implementation of a vortex force formalism in a coupled modeling system
** for inner-shelf and surf-zone applications. Ocean Modelling, 47, 65-95.
**
** Options for Rip_current.
** Application flag:   RIP_CURRENT
** Input scripts:      rip_current.h
*/

#define ROMS_MODEL
#define SWAN_MODEL
#define MCT_LIB

#define WEC_VF
#define WDISS_WAVEMOD
#define ROLLER_RENIERS
#define UV_KIRBY
#define WET_DRY
#define MASKING

#define OUT_DOUBLE
#define UV_ADV
#define DJ_GRADPS
#undef  SALINITY
#define UV_VIS2
#define MIX_S_UV

#undef  DIAGNOSTICS_UV
#undef  AVERAGES
#undef  AVERAGES_WEC

#undef  ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX

#define SOLVE3D
#ifdef SOLVE3D
# define SPLINES_VVISC
# define SPLINES_VDIFF
# define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# define SSW_LOGINT
/* define one of these 2 */
# define SSW_LOGINT_WBL
# undef  SSW_LOGINT_DIRECT
#endif

# undef  SEDIMENT
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
# else
#  define TS_FIXED
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

# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define RI_SPLINES
#  undef  ZOS_HSIG
#  undef  TKE_WAVEDISS
# endif

#else
# define UV_QDRAG
#endif
