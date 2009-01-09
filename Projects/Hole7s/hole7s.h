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
#define WET_DRY
#undef  SWAN_COUPLING
#undef  SVENDSEN_ROLLER

#undef  DIAGNOSTICS_UV
#define UV_VIS2
#define MIX_S_UV
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

#undef  ANA_INITIAL

#undef  ANA_SMFLUX
#define ANA_FSOBC
#define ANA_M2OBC
#undef  ANA_TOBC
#ifdef SOLVE3D
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif

/* BOUNDARY CONDITIONS*/
#define WESTERN_WALL
#define EAST_FSCHAPMAN
#define EAST_M2REDUCED
#define FSOBC_REDUCED
#define EAST_M3GRADIENT
#define EAST_TGRADIENT

#define RADIATION_2D
#define NORTH_FSRADIATION
#define NORTH_M2RADIATION
#define NORTH_M3RADIATION
#define NORTH_TGRADIENT
#define SOUTH_FSRADIATION
#define SOUTH_M2RADIATION
#define SOUTH_M3RADIATION
#define SOUTH_TGRADIENT


#define GLS_MIXING
#if defined GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# undef  CRAIG_BANNER
# undef  CHARNOK
# undef  ZOS_HSIG
# undef  TKE_WAVEDISS
#endif

#undef SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
# define BEDLOAD_SOULSBY
# define ANA_SEDIMENT
#endif

#define UV_LOGDRAG
#undef  SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# undef  SSW_CALC_UB
# undef  SSW_LOGINT
# undef  SSW_Z0RIP
#endif

#undef  ANA_WWAVE

