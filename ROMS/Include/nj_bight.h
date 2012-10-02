/*
** svn $Id: nj_bight.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for 4DVar Data Assimilation Toy.
**
** Application flag:   NJ_BIGHT
** Input script:       ocean_nj_bight.in
*/

#define ROMS_MODEL
#define UV_ADV
#define UV_SADVECTION
#define UV_QDRAG
#define UV_COR
#undef  UV_PSOURCE
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  TS_SVADVECTION
#undef  TS_A4HADVECTION
#undef  TS_A4VADVECTION
#undef  TS_DIF2
#undef  MIX_GEO_TS
#undef  TS_PSOURCE
#define SOLAR_SOURCE

#define NONLIN_EOS
#define SALINITY
#define CURVGRID
#define SOLVE3D
#define MASKING
#define SPLINES
#define AVERAGES
#ifdef AVERAGES
# define AVERAGES_AKV
# define AVERAGES_AKT
# define AVERAGES_AKS
# define AVERAGES_FLUXES
#endif
#define STATIONS
#undef  FLOATS

#undef  ASSIMILATION_UVsur
#undef  ASSIMILATION_T
#undef  NUDGING_T
#undef  NUDGING_UVsur

#define WESTERN_WALL
#define NORTHERN_WALL
#define RADIATION_2D
#define EAST_M3RADIATION
#define EAST_TRADIATION
#define SOUTH_M3RADIATION
#define SOUTH_TRADIATION

#define GLS_MIXING
#undef  MY25_MIXING
#undef  LMD_MIXING

#ifdef GLS_MIXING
# define N2S2_HORAVG
# undef  KANTHA_CLAYSON
#endif

#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif

#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef  LMD_BKPP
# define LMD_NONLOCAL
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
# define EAST_FSCHAPMAN
# define SOUTH_FSCHAPMAN
#else
# define EAST_FSGRADIENT
# define SOUTH_M2RADIATION
#endif

#define UV_TIDES
#ifdef UV_TIDES
# define ANA_M2OBC
# define EAST_M2FLATHER
# define SOUTH_M2FLATHER
#else
# define EAST_M2RADIATION
# define SOUTH_FSGRADIENT
#endif
#if defined SSH_TIDES || defined UV_TIDES
# undef  EAST_VOLCONS
# undef  SOUTH_VOLCONS
#else
# define EAST_VOLCONS
# define SOUTH_VOLCONS
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
#define ALBEDO
#define ANA_SMFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#define INLINE_2DIO
