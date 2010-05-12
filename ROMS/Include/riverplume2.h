/*
** svn $Id: riverplume2.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for River Plume example by Hyatt and Signell described at
** http://smig.usgs.gov/SMIG/features_0300/plumes_inline.html
**
** Application flag:   RIVERPLUME2
** Input script:       ocean_riverplume2.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_PSOURCE
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define TS_PSOURCE
#define NONLIN_EOS
#define SALINITY
#define MASKING
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS

#define SOUTH_FSCHAPMAN
#define SOUTH_M2GRADIENT
#define SOUTH_M3GRADIENT
#define SOUTH_TGRADIENT
#define NORTH_FSCHAPMAN
#define NORTH_M2GRADIENT
#define NORTH_M3GRADIENT
#define NORTH_TGRADIENT
#define WESTERN_WALL
#define EASTERN_WALL

#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
#endif

#define ANA_GRID
#define ANA_INITIAL
#define ANA_PSOURCE
#define ANA_SMFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

