/*
** svn $Id: rip_current.h 1009 2007-08-19 16:51:29Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for RIP_CURRENT Application.
**
** Rip Current - test for rip current.
**
*/

#define UV_VIS2
#define MIX_S_UV
#define DIAGNOSTICS_UV
#define AVERAGES
#define AVERAGES_NEARSHORE
#define NEARSHORE_MELLOR
#undef  SVENDSEN_ROLLER
#undef  MONO_ROLLER
#define REFDIF_COUPLING
#define UV_ADV
#define TS_U3HADVECTION
#define DJ_GRADPS
#undef  SALINITY
#define SOLVE3D
#define SPLINES
#define EASTERN_WALL
#define WESTERN_WALL
#define NORTHERN_WALL
#define SOUTHERN_WALL

#undef UV_QDRAG
#define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# define SSW_CALC_UB
# define ANA_SEDIMENT
# undef  SSW_LOGINT
#endif

#ifdef SOLVE3D
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
#endif
#define ANA_INITIAL
#define ANA_SMFLUX
#define GLS_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define ZOS_HSIG
# define TKE_WAVEDISS
#endif

