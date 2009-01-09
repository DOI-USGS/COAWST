/*
** svn $Id: canyon.h 503 2008-01-10 00:11:51Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for an Idealized Canyon.
**
** Application flag:   CANYON
** Input script:       ocean_canyon2d.in, ocean_canyon3d.in
*/

#ifndef SOLVE3D                   /* 2D set-up */
# define UV_ADV
# define UV_QDRAG
# define UV_VIS2
# define UV_COR
# define EW_PERIODIC
# define BODYFORCE
# define ANA_DIAG
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
#else                             /* 3D set-up */
# define UV_ADV
# define UV_COR
# define UV_QDRAG
# define UV_VIS2
# define MIX_S_UV
# define DJ_GRADPS
# define TS_A4HADVECTION
# define TS_A4VADVECTION
# define TS_DIF2
# define MIX_GEO_TS
# define SPLINES
# define EW_PERIODIC
# define ANA_DIAG
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_VMIX
#endif
