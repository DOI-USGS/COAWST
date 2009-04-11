/*
** svn $Id: south_car.h 996 2007-08-10 00:57:22Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for SOUTH_CAR Application.
**
*/
#define NEARSHORE_MELLOR
#define WET_DRY
#define SWAN_COUPLING
#define MCT_LIB
#define MCT_INTERP_OC2WV
#undef  SVENDSEN_ROLLER
#undef  DIAGNOSTICS_UV
#undef  UV_VIS2
#undef  MIX_S_UV
#define UV_ADV
#undef  UV_COR
#define SALINITY
#undef  CURVGRID
#define SOLVE3D
#define MASKING
#define DJ_GRADPS
#undef  SPLINES
#define TS_U3HADVECTION
#undef  TS_MPDATA
#define ANA_INITIAL
#define ANA_SMFLUX
#ifdef SOLVE3D
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif
#define WESTERN_WALL
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3GRADIENT
#define EAST_TGRADIENT

#define RADIATION_2D
/* #define FSOBC_REDUCED*/
#define NORTH_FSRADIATION
#define NORTH_M2RADIATION
#define NORTH_M3RADIATION
#define NORTH_TGRADIENT
#define SOUTH_FSRADIATION
#define SOUTH_M2RADIATION
#define SOUTH_M3RADIATION
#define SOUTH_TGRADIENT

/*#define NORTH_FSGRADIENT
  #define NORTH_M2GRADIENT
  #define NORTH_M3GRADIENT
  #define NORTH_TGRADIENT
  #define SOUTH_FSGRADIENT
  #define SOUTH_M2GRADIENT
  #define SOUTH_M3GRADIENT
  #define SOUTH_TGRADIENT*/

/* #define NORTH_FSGRADIENT
   #define NORTH_M2REDUCED
   #define NORTH_M3GRADIENT
   #define NORTH_TGRADIENT
   #define SOUTH_FSGRADIENT
   #define SOUTH_M2REDUCED
   #define SOUTH_M3GRADIENT
   #define SOUTH_TGRADIENT */

#define ANA_FSOBC
#define ANA_M2OBC
#undef  ANA_TOBC
#define GLS_MIXING
#if defined GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# undef  CRAIG_BANNER
# undef  CHARNOK
# undef  ZOS_HSIG
# undef  TKE_WAVEDISS
#endif
#define ANA_SEDIMENT
#undef SEDIMENT
#ifdef SEDIMENT
# undef SUSPLOAD
# undef BEDLOAD_SOULSBY
#endif
#undef UV_LOGDRAG
#define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# undef  SSW_CALC_UB
# undef  SSW_LOGINT
# undef  SSW_Z0RIP
#endif
#undef  ANA_WWAVE

