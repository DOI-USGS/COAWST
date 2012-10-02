/*
** svn $Id: estuary_test.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Estuary with Sediment Transport Test.
**
** Application flag:   ESTUARY_TEST
** Input script:       ocean_estuary_test.in
**                     sediment_estuary_test.in
*/

#define ROMS_MODEL
#define UV_ADV
#define UV_LOGDRAG
#define TS_U3HADVECTION
#define SALINITY
#define SOLVE3D
#define SPLINES
#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
#endif
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKS
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EAST_FSGRADIENT
#define EAST_M2CLAMPED
#define EAST_M3GRADIENT
#define EAST_TCLAMPED
#define WEST_FSCHAPMAN
#define WEST_M2CLAMPED
#undef  WEST_M2REDUCED
#undef  FSOBC_REDUCED
#define WEST_M3GRADIENT
#define WEST_TRADIATION
#define WEST_TNUDGING
#define GLS_MIXING
#undef  MY25_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# undef  CANUTO_A
# define N2S2_HORAVG
#endif
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SEDIMENT
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_TOBC
