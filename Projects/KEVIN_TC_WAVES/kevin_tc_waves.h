/*
** svn $Id: upwelling.h 25 2007-04-09 23:43:58Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for KEVIN_TC Test.
**
** Application flag:   kevin_tc
** Input script:       ocean_upwelling.in
*/

#define WRF_COUPLING
#define SWAN_COUPLING
#undef  NEARSHORE_MELLOR
#undef  SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# define ANA_SEDIMENT
#else
# define UV_LDRAG
#endif
#define UV_ADV
#define UV_COR
#define UV_VIS2
#undef  MIX_GEO_UV
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  TS_MPDATA
#undef  WJ_GRADP
#define DJ_GRADPS
#define TS_DIF2
#undef  TS_DIF4
#undef  MIX_GEO_TS
#define MIX_S_TS

#define SALINITY
#define SOLVE3D
#define SPLINES
#define AVERAGES
#undef DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV

#define EW_PERIODIC

#define ANA_GRID
#define ANA_INITIAL
#if !defined SWAN_COUPLING && defined SSW_BBL
# define ANA_WWAVE
#endif
#ifdef WRF_COUPLING
# undef ANA_SMFLUX
# undef ANA_STFLUX
# define BULK_FLUXES
# define LONGWAVE_OUT
#else
# define ANA_SMFLUX
# define ANA_STFLUX
#endif
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_SSFLUX

#define ANA_VMIX
#undef  GLS_MIXING
#undef  MY25_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# undef  CANUTO_A
# define N2S2_HORAVG
#endif

#undef  BIO_FASHAM
#undef  NPZD_POWELL
#undef  ECOSIM

#if defined BIO_FASHAM || defined ECOSIM || defined NPZD_POWELL
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_SRFLUX
#endif

#ifdef BIO_FASHAM
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
#endif

#undef  PERFECT_RESTART
#ifdef PERFECT_RESTART
# undef  AVERAGES
# undef  DIAGNOSTICS_BIO
# undef  DIAGNOSTICS_TS
# undef  DIAGNOSTICS_UV
# define OUT_DOUBLE
#endif
