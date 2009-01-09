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
#define SOLVE3D
#define SPLINES
#define UV_COR
#define UV_LDRAG
#define UV_ADV
#define UV_VIS2
#define UV_SADVECTION
#define TS_FIXED
#define MIX_GEO_UV
#define ANA_VMIX
#define SALINITY
#define ANA_GRID
#define ANA_INITIAL
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
#define ANA_SMFLUX
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EW_PERIODIC

#undef MASKING
#undef TS_U3HADVECTION


