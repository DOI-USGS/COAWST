/*
** svn $Id: soliton.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Equatorial Soliton Test.
**
** Application flag:   SOLITON
** Input script:       ocean_soliton.in
*/

#define ROMS_MODEL
#define UV_ADV
#define UV_C4ADVECTION
#define UV_VIS2
#define UV_COR
#define UV_QDRAG
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
