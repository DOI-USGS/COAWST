/*
** svn $Id: overflow.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Gravitational Overflow Test.
**
** Application flag:   OVERFLOW
** Input script:       ocean_overflow.in
*/

#define ROMS_MODEL
#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define TS_DIF2
#define MIX_ISO_TS
#define SOLVE3D
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
