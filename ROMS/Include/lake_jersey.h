/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Lake Jersey forced with wind.
**
** Application flag: LAKE_JERSEY
**
** Input scripts:    ocean_lake_jersey.in
**                   sediment_lake_jersey.in
*/

#define ROMS_MODEL
#define UV_ADV
#define UV_COR
#define DJ_GRADPS
#define NONLIN_EOS
#define SALINITY
#define SOLVE3D
#define MASKING

#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BPFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_SRFLUX

/* Bottom boundary layer options */

#defined SSW_BBL
#ifdef SG_BBL
# define SG_CALC_ZNOT
# undef  SG_LOGINT
#endif
#ifdef MB_BBL
# define MB_CALC_ZNOT
# undef  MB_Z0BIO
# undef  MB_Z0BL
# undef  MB_Z0RIP
#endif
#ifdef SSW_BBL
# define SSW_CALC_UB
# define SSW_CALC_ZNOT
# undef  SSW_LOGINT
#endif

#if defined MB_BBL || defined SG_BBL || defined SSW_BBL
# define ANA_WWAVE
#endif

#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
# define BEDLOAD_SOULSBY
# undef  SED_MORPH
#endif
#if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
# define ANA_SEDIMENT
#endif

