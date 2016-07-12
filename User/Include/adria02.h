/*
** svn $Id: adria02.h 795 2016-05-11 01:42:43Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Adriatic Sea Application.
**
** Application flag:   ADRIA02
** Input script:       ocean_adria02.in
**                     floats_adria02.in
**                     sediment_adria02.in
**                     stations_adria02.in
*/

#define UV_ADV
#define UV_COR
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_MPDATA
#define TS_DIF2
#define MIX_GEO_TS
#define NONLIN_EOS
#define SALINITY
#define MASKING
#define SOLVE3D
#define STATIONS
#define CURVGRID
#define FLOATS
#define AVERAGES

#undef NOSEDBBL
#ifdef NOSEDBBL
# undef SEDIMENT
# undef SUSPLOAD
# define ANA_SEDIMENT
# undef  ANA_WWAVE
#else
# define SEDIMENT
# define SUSPLOAD
# undef  ANA_SEDIMENT
# undef  ANA_WWAVE
#endif

#undef  UV_LOGDRAG
#undef  MB_BBL
#undef  SG_BBL
#define SSW_BBL

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

#undef MY25_MIXING
#define GLS_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define CRAIG_BANNER
# define CHARNOK
# define RI_SPLINES
#endif

#undef ANA_SRFLUX
#undef ALBEDO_CLOUD
#define DIURNAL_SRFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BPFLUX
#define ANA_BTFLUX
#define ANA_SPFLUX

#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE
# undef SOLAR_SOURCE
# define ANA_RAIN
# undef COOL_SKIN
#endif

#define RADIATION_2D

#define RAMP_TIDES
#define SSH_TIDES
#ifdef SSH_TIDES
# define ANA_FSOBC
#endif
#define UV_TIDES
#ifdef UV_TIDES
# define ANA_M2OBC
#endif
