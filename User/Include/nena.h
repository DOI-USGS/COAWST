/*
** svn $Id: nena.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for North East North America Application.
**
** Application flag:   NENA
** Input script:       ocean_nena.in
*/

#define UV_ADV
#define UV_SADVECTION
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define UV_PSOURCE
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define TS_PSOURCE
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPLINES
#define MASKING
#define AVERAGES
#define SRELAXATION
#define QCORRECTION
#define SOLAR_SOURCE
#define DIURNAL_SRFLUX
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif
#undef  GLS_MIXING
#define BIO_FENNEL
#ifdef BIO_FENNEL
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif
#undef  M2CLIMATOLOGY
#undef  M3CLIMATOLOGY
#undef  TCLIMATOLOGY
#undef  ZCLIMATOLOGY
#undef  EAST_VOLCONS
#undef  WEST_VOLCONS
#undef  SOUTH_VOLCONS
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3CLAMPED
#define EAST_TCLAMPED
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3CLAMPED
#define WEST_TCLAMPED
#define NORTHERN_WALL
#define SOUTH_FSCHAPMAN
#define SOUTH_M2FLATHER
#define SOUTH_M3CLAMPED
#define SOUTH_TCLAMPED
#define ANA_BSFLUX
#define ANA_BTFLUX
