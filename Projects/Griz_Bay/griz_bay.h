/*
** svn $Id: grizzly_bay.h $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Grizzly Bay wet_dry test case.
**
** Application flag:   GRIZZLY_BAY
*/


#define UV_COR
#define UV_ADV
#define MASKING
#define CURVGRID
#define SOLVE3D
#define UV_QDRAG
#define WET_DRY

#undef  SALINITY
#define SPLINES
#define TS_FIXED
#define ANA_INITIAL
#define ANA_SMFLUX
#ifdef SOLVE3D
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif

#define NORTHERN_WALL
#define SOUTHERN_WALL
#define FSOBC_REDUCED
#define ANA_M2OBC
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3GRADIENT
#define EAST_TGRADIENT
#define WEST_FSCHAPMAN
#define WEST_M2REDUCED
#define WEST_M3GRADIENT
#define WEST_TGRADIENT

#define GLS_MIXING
#if defined GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif

#define DIAGNOSTICS_UV
