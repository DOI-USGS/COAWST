/*
** svn $Id: shoreface.h 139 2008-01-10 00:17:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Shore Face Planar Beach Test Case.
**
** Application flag:   WETDRY_SLOPE_CHAN
*/

#define ROMS_MODEL

#define WET_DRY
#undef  UV_ADV
#define TS_FIXED
#define ANA_SMFLUX
#define UV_LDRAG

#undef SOLVE3D
#ifdef SOLVE3D
# define DJ_GRADPS
# undef  SALINITY
# define EAST_M3GRADIENT
# define SPLINES
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
# undef  ANA_VMIX
# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
#endif
