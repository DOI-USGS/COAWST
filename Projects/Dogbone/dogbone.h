/*
** svn $Id: dogbone.h 8 2007-02-06 19:00:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for DOGBONE Application.
**
** DOGBONE - composite grid test
**           Just change Ngrids in the makefile, do not make any changes here.
*/

#define OUT_DOUBLE
#undef  UV_ADV
#undef  UV_VIS4
#undef  MIX_S_UV
#define UV_QDRAG
#define COMPOSED_GRID
#undef  REFINED_GRID
#define ANA_SMFLUX
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EASTERN_WALL
#define WESTERN_WALL
#define MASKING

#define SOLVE3D
#ifdef SOLVE3D
# undef  TS_FIXED
# define SALINITY
# define SPLINES
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_BPFLUX
# define TS_U3HADVECTION
# undef  TS_MPDATA
# undef  ANA_VMIX
# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
# undef SEDIMENT
# undef SUSPLOAD
#endif


