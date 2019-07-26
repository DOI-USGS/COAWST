/*
** svn $Id: coupling_test.h 889 2018-02-10 03:32:52Z arango $
*******************************************************************************
** Copyright (c) 2002-2019 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Atmosphere-Ocean Two-Way Coupling Test (WRF/ROMS).
**
** Application flag:   COUPLING_TEST
** Input script:       roms_coupling_test.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#undef  UV_VIS2
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define DJ_GRADPS
#undef  TS_DIF2
#undef  MIX_GEO_TS
#define SALINITY
#define SOLVE3D
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_VMIX
