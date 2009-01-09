/*
** svn $Id: upwelling.h 25 2007-04-09 23:43:58Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for JOE_TC Test.
**
** Application flag:   JOE_tc
*/

#define WRF_COUPLING
/*#define DIAGNOSTICS_TS
  #define DIAGNOSTICS_UV*/
#define SWAN_COUPLING

#define NEARSHORE_MELLOR
#define SSW_BBL

/* Physics + numerics */
#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  TS_MPDATA
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# define ANA_SEDIMENT
#else
# define UV_QDRAG
#endif
#if !defined SWAN_COUPLING && defined SSW_BBL
# define ANA_WWAVE
#endif

#define DJ_GRADPS
#define TS_DIF2
#define MIX_GEO_TS

#define SALINITY
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define NONLIN_EOS

/* Grid and Initial */
#define ANA_GRID
#define ANA_INITIAL

/* Forcing */
#ifdef WRF_COUPLING
# define BULK_FLUXES
# define ANA_SSFLUX
# define LONGWAVE_OUT
#else
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
#endif

#define ANA_BTFLUX
#define ANA_BSFLUX

/* Turbulence closure */
#undef  GLS_MIXING
#define MY25_MIXING

#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif

/* Boundary condition */
#define WESTERN_WALL
#define RADIATION_2D

#define EAST_FSRADIATION
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define EAST_TRADIATION

#define NORTH_FSRADIATION
#define NORTH_M2RADIATION
#define NORTH_M3RADIATION
#define NORTH_TRADIATION

#define SOUTH_FSRADIATION
#define SOUTH_M2RADIATION
#define SOUTH_M3RADIATION
#define SOUTH_TRADIATION

/* Biological module */
#undef  NPZD_POWELL

#if defined NPZD_POWELL
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_SRFLUX
#endif
