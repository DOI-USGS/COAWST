/*
** svn $Id: nj_bight.h 795 2016-05-11 01:42:43Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for New York/New Jersey Bight.
**
** Application flag:   NJ_BIGHT
** Input script:       ocean_nj_bight.in
*/

#define UV_ADV
#define UV_SADVECTION
#define UV_QDRAG
#define UV_COR
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  TS_SVADVECTION
#undef  TS_A4HADVECTION
#undef  TS_A4VADVECTION
#undef  TS_DIF2
#undef  MIX_GEO_TS
#define SOLAR_SOURCE

#define NONLIN_EOS
#define SALINITY
#define CURVGRID
#define SOLVE3D
#define MASKING
#define AVERAGES
#define STATIONS
#undef  FLOATS

#define RADIATION_2D

#define GLS_MIXING
#undef  MY25_MIXING
#undef  LMD_MIXING

#ifdef GLS_MIXING
# define N2S2_HORAVG
# define RI_SPLINES
# undef  KANTHA_CLAYSON
#endif

#ifdef MY25_MIXING
# define N2S2_HORAVG
# define RI_SPLINES
# define KANTHA_CLAYSON
#endif

#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef  LMD_BKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif

#define BULK_FLUXES
#ifdef BULK_FLUXES
# define ANA_RAIN
# define LONGWAVE
# ifdef LONGWAVE
#  define ANA_CLOUD
# endif
#endif

#undef SG_BBL
#ifdef SG_BBL
# define SG_CALC_UB
# define SG_CALC_ZNOT
# define ANA_SEDIMENT
# define ANA_WWAVE
#endif

#define RAMP_TIDES
#define SSH_TIDES
#ifdef SSH_TIDES
# define ANA_FSOBC
#endif

#define UV_TIDES
#ifdef UV_TIDES
# define ANA_M2OBC
#endif

#undef  BIO_FENNEL
#undef  ECOSIM

#ifdef BIO_FENNEL
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
#endif

#if defined BIO_FENNEL || defined ECOSIM
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif

#define ANA_SRFLUX
#define ALBEDO_CLOUD
#define ANA_SMFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#define INLINE_2DIO
