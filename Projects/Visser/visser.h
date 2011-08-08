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

#define ROMS_MODEL
#define DIAGNOSTICS_UV

/* WEC forcings */
#define WEC_VF
#define ROLLER_SVENDSEN
#define ROLLER_MONO

#define UV_ADV
#define TS_U3HADVECTION
#define DJ_GRADPS
#define UV_QDRAG

#define EASTERN_WALL
#define WESTERN_WALL
#define NS_PERIODIC

#define SOLVE3D
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

