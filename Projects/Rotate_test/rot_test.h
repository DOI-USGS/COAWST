/*
** svn $Id: upwelling.h 25 2007-04-09 23:43:58Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for ROT_TEST Test.
**
** Application flag:   ROT_TEST
*/

#define ROMS_MODEL
#define SWAN_MODEL
#define MCT_LIB

/* Physics + numerics */
#define UV_ADV
#undef  UV_COR
#define UV_VIS2
#define MIX_S_UV
#define TS_U3HADVECTION
#define UV_QDRAG
#define WEC_VF
#define WDISS_WAVEMOD
#define UV_KIRBY

#define CURVGRID
#define DJ_GRADPS

#define SALINITY
#define SOLVE3D
#define SPLINES

/* Grid and Initial */
#define MASKING

/* Forcing */
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX

/* Turbulence closure */
#define GLS_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif

/* Boundary condition */

#define NORTHERN_WALL
#define SOUTHERN_WALL
#undef WESTERN_WALL
#undef EASTERN_WALL

#define WEST_FSCHAPMAN
#define WEST_M2GRADIENT
#define WEST_M2SGRADIENT
#define WEST_M3SGRADIENT
#define WEST_M3GRADIENT
#define WEST_TGRADIENT

#define FSOBC_REDUCED
#define EAST_FSCHAPMAN
#define EAST_M2CLAMPED
#define EAST_M2SGRADIENT
#define EAST_M3SGRADIENT
#define EAST_M3GRADIENT
#define EAST_TGRADIENT

#define ANA_FSOBC
#define ANA_M2OBC
