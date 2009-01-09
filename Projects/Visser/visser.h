/*
** svn $Id: visser.h 999 2007-08-13 14:17:47Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for VISSER Application.
**
** VISSER - test for simple beach
**
*/

#define DIAGNOSTICS_UV
#define UV_VIS2
#define MIX_GEO_UV
#define SPONGE
#define NEARSHORE_MELLOR
#define SVENDSEN_ROLLER
#define MONO_ROLLER
#define UV_ADV
#define TS_U3HADVECTION
#define DJ_GRADPS
#define SOLVE3D
#define SPLINES
#define EASTERN_WALL
#define WESTERN_WALL
/*#define WEST_M2FLATHER*/
#define NS_PERIODIC
#undef  UV_QDRAG
#undef  UV_LOGDRAG
#define SSW_BBL
#ifdef  SSW_BBL
# define SSW_CALC_ZNOT
# define SSW_CALC_UB
# define ANA_SEDIMENT
#endif
#ifdef SOLVE3D
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
#  define ZOS_HSIG
#  define TKE_WAVEDISS
# endif
#endif
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX

