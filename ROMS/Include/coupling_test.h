/*
** svn $Id: coupling_test.h 503 2008-01-10 00:11:51Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Atmosphere-Ocean Two-Way Coupling Test (WRF/ROMS).
**
** Application flag:   COUPLING_TEST
** Input script:       ocean_coupling_test.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#undef  UV_VIS2
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define DJ_GRADPS
#undef  TS_DIF2
#undef  MIX_GEO_TS
#define SALINITY
#define EW_PERIODIC
#define NS_PERIODIC
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define AIR_OCEAN
#define ANA_GRID
#define ANA_INITIAL
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_VMIX
