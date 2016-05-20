/*
** svn $Id: windbasin.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Wind-Driven Constant Coriolis Basin Test.
**
** Application flag:   WINDBASIN
** Input script:       ocean_windbasin.in
*/

#define ROMS_MODEL
#undef UV_ADV
#define UV_COR
#define UV_QDRAG
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define SOLVE3D
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX

