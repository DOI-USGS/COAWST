/*
** svn $Id: weddell.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Idealized Weddell Sea Application: Tides and Ice Shelf Test.
**
** Application flag:   WEDDELL
** Input script:       ocean_wedell.in
*/

#define ROMS_MODEL
#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#undef  UV_VIS4
#undef  MIX_S_UV
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#undef  TS_DIF4
#undef  MIX_GEO_TS
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define ICESHELF
#define AVERAGES
#define RADIATION_2D
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
