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
#undef  NEARSHORE_MELLOR
#undef  WET_DRY
#define UV_VIS2
#define MIX_S_UV
#define UV_ADV
#define UV_COR
#define SALINITY
#define CURVGRID
#define SOLVE3D
#define MASKING
#define DJ_GRADPS
#define SPLINES
#undef  TS_U3HADVECTION
#define TS_MPDATA
#undef  ANA_INITIAL
#undef  ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
#define ATM_PRESS
#define WESTERN_WALL
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3GRADIENT
#define EAST_TGRADIENT
#define NORTH_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_M3GRADIENT
#define NORTH_TGRADIENT
/*#undef  SOUTH_FSRADIATION
  #undef  SOUTH_M2RADIATION
  #define FSOBC_REDUCED
  #define SOUTH_FSGRADIENT
  #define SOUTH_M2GRADIENT*/
#define SOUTH_FSCHAPMAN
#define SOUTH_M2FLATHER
#define SOUTH_M3GRADIENT
#define SOUTH_TGRADIENT
#define SSH_TIDES
#define UV_TIDES
#define RAMP_TIDES
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
#define SEDIMENT
#define SUSPLOAD
#undef  BEDLOAD_SOULSBY
#undef  ANA_SEDIMENT
#undef  UV_LOGDRAG
#undef  MB_BBL
#ifdef  MB_BBL
# define MB_CALC_ZNOT
# undef  MB_CALC_UB
# define MB_LOGINT
# undef  MB_Z0RIP
#endif
#define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# undef  SSW_CALC_UB
# undef  SSW_LOGINT
# undef  SSW_Z0RIP
#endif
#undef  ANA_WWAVE

