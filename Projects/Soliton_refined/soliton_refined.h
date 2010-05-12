/*
** svn $Id: soliton.h 503 2008-01-10 00:11:51Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
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
#define REFINED_GRID

#define UV_ADV
#define UV_C4ADVECTION
#define UV_VIS2
#define UV_COR
#define UV_QDRAG
#define ANA_INITIAL
#define EW_PERIODIC
#define EW_PERIODIC_REFINED
#define ANA_SMFLUX

