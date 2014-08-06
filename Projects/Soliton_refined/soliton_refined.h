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
#define NESTING
#define OUT_DOUBLE

#define UV_ADV
#define UV_COR
#undef  UV_C4ADVECTION
#define UV_C2ADVECTION
#define UV_VIS2

#define UV_QDRAG
#define ANA_INITIAL
#define EW_PERIODIC
#define ANA_SMFLUX

