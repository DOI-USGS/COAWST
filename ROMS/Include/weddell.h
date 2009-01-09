/*
** svn $Id: weddell.h 503 2008-01-10 00:11:51Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Idealized Weddell Sea Application: Tides and Ice Shelf Test.
**
** Application flag:   WEDDELL
** Input script:       ocean_wedell.in
*/

#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#undef  UV_VIS4
#undef  MIX_S_UV
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#undef  TS_DIF4
#undef  MIX_GEO_TS
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPLINES
#define ICESHELF
#define AVERAGES
#define NS_PERIODIC
#define EAST_VOLCONS
#define WEST_VOLCONS
#define RADIATION_2D
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3RADIATION
#define EAST_TRADIATION
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3RADIATION
#define WEST_TRADIATION
#define ZCLIMATOLOGY
#define M2CLIMATOLOGY
#define ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_SRFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
