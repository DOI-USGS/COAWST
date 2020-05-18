/*
** svn $Id: kelvin.h 995 2020-01-10 04:01:28Z arango $
*******************************************************************************
** Copyright (c) 2002-2020 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Kelvin Wave Test.
**
** Application flag:   KELVIN
** Input script:       roms_kelvin.in
*/

#define ROMS_MODEL
#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_DIF2
#define MIX_S_TS
#define SOLVE3D
#define RADIATION_2D
#define ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SRFLUX
#define ANA_BTFLUX
