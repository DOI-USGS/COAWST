/*
** svn $Id: soliton.h 995 2020-01-10 04:01:28Z arango $
*******************************************************************************
** Copyright (c) 2002-2020 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Equatorial Soliton Test.
**
** Application flag:   SOLITON
** Input script:       roms_soliton.in
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
