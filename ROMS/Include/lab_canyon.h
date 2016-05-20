/*
** svn $Id: lab_canyon.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Lab Canyon (Polar Coordinates).
**
** Application flag:   LAB_CANYON
** Input script:       ocean_lab_canyon.in
*/

#define ROMS_MODEL
#define UV_COR
#define UV_ADV
#define UV_LDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define SOLVE3D
#define CURVGRID
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
