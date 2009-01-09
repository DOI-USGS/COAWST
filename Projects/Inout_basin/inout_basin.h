/*
** svn $Id: inout_basin.h 682 2007-04-24 18:32:25Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for INFLOW OUTLFLOW BASIN Application.
**
**  INOUT BASIN
*/

#define ANA_GRID
#define MASKING
#define ANA_MASK
#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define SALINITY
#define SOLVE3D
#define SPLINES
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
#undef  TS_U3HADVECTION
#define TS_MPDATA
/*#define NORTHERN_WALL
  #define SOUTHERN_WALL
  #define EASTERN_WALL*/
#define EAST_FSRADIATION
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define NORTH_FSRADIATION
#define NORTH_M2RADIATION
#define NORTH_M3RADIATION
#define SOUTH_FSRADIATION
#define SOUTH_M2RADIATION
#define SOUTH_M3RADIATION
#define WESTERN_WALL
#define UV_PSOURCE
#define TS_PSOURCE
#define ANA_PSOURCE
#define ANA_VMIX
#undef  GLS_MIXING
#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif
#define UV_QDRAG
