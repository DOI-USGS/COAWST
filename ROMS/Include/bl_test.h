/*
** svn $Id: bl_test.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Boundary Layers Test.
**
** Application flag:   BL_TEST
** Input script:       ocean_bl_test.in
*/

#define UV_ADV
#define UV_SADVECTION
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define SPLINES
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS
#define STATIONS
#define SOLVE3D
#define WESTERN_WALL
#define NS_PERIODIC
#define RADIATION_2D
#define EAST_FSGRADIENT
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define EAST_KRADIATION
#define EAST_TRADIATION
#define EAST_VOLCONS
#undef  MY25_MIXING
#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
# define LMD_DDMIX
#endif
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE
# define ANA_CLOUD
# define ANA_HUMIDITY
# define ANA_PAIR
# define ANA_TAIR
# define ANA_RAIN
# define ANA_WINDS
#else
# define ANA_SMFLUX
# define ANA_STFLUX
#endif
#define SG_BBL
#ifdef SG_BBL
# define SG_CALC_UB
# define SG_CALC_ZNOT
# define ANA_SEDIMENT
# define ANA_WWAVE
#else
# define UV_QDRAG
#endif
#define ANA_GRID
#define ANA_INITIAL
#define ALBEDO
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
