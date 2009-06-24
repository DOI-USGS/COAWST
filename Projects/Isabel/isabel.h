/*
** svn $Id: isabel.h 996 2007-08-10 00:57:22Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Isabel Application.
**
*/

!# define INLINE_2DIO
!# define PERFECT_RESTART
!# define UV_SADVECTION
!# define TS_FIXED
!# define AVERAGES


#define SWAN_COUPLING
#define MCT_INTERP_WV2AT
#define MCT_INTERP_OC2AT
#undef MCT_INTERP_OC2WV
#define WRF_COUPLING
#if defined SWAN_COUPLING || defined WRF_COUPLING
# define MCT_LIB
#endif
#define WET_DRY
#define ATM_PRESS
#undef  NEARSHORE_MELLOR

#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#undef  UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
/*# define TS_U3HADVECTION
  # define TS_SVADVECTION */
#define TS_MPDATA
#define TS_DIF2
#define MIX_GEO_TS

#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPLINES            
#define MASKING

/* define River flows */
!#define  UV_PSOURCE
!#define  TS_PSOURCE

/* surface forcing */
#define BULK_FLUXES
#ifdef BULK_FLUXES
#undef LONGWAVE_OUT
!# define LONGWAVE
#  define ANA_SSFLUX
#  define COARE_TAYLOR_YELLAND
#  define EMINUSP
#  define SOLAR_SOURCE
#else
#  define ANA_SSFLUX
#  define ANA_SMFLUX
#  define ANA_STFLUX
#endif

#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX

/* Select turbulence mixing scheme */

#define GLS_MIXING
#ifdef GLS_MIXING
# define AKLIMIT
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define CRAIG_BANNER
# define CHARNOK
# undef  ZOS_HSIG
# undef  TKE_WAVEDISS
#endif

/*define boundary conditon */
#define NORTHERN_WALL
#define WESTERN_WALL

#define SSH_TIDES
#define ADD_FSOBC
#define UV_TIDES
#define ADD_M2OBC
#define RAMP_TIDES

#define RADIATION_2D

#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3RADIATION
#define EAST_TRADIATION
#define EAST_TNUDGING
#define EAST_M3NUDGING

#define SOUTH_FSCHAPMAN
#define SOUTH_M2FLATHER
#define SOUTH_M3RADIATION
#define SOUTH_TRADIATION
#define SOUTH_TNUDGING
#define SOUTH_M3NUDGING

/* define water mass relaxation */
#define TCLIMATOLOGY
#define TCLM_NUDGING

#define ANA_FSOBC
#define ANA_M2OBC

#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
# define BEDLOAD_SOULSBY
# undef  ANA_SEDIMENT
# undef  UV_LOGDRAG
#endif
#define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# undef  SSW_CALC_UB
# undef  SSW_LOGINT
# undef  SSW_Z0RIP
#endif


