/*
** svn $Id: sed_test1.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Suspended Sediment Test in a Channel.
**
** Application flag:   SED_TEST1
** Input scripts:      ocean_sed_test1.in
**                     sediment_sed_test1.in
*/

#define UV_ADV
#define UV_PSOURCE
#define UV_LOGDRAG
#define UV_VIS4
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define TS_DIF4
#define MIX_S_TS
#define SALINITY
#define SOLVE3D
#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
#endif
#define  SPLINES
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define WEST_FSRADIATION
#define WEST_M2RADIATION
#define WEST_M3RADIATION
#define WEST_TGRADIENT
#define EAST_FSCLAMPED
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define EAST_TCLAMPED
#define MY25_MIXING
#ifdef MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif
#define ANA_BPFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SEDIMENT
#define ANA_SMFLUX
#define ANA_SPFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_PSOURCE
#define ANA_TOBC
#define ANA_FSOBC
