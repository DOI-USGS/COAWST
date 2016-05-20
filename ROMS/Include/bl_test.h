/*
** svn $Id: bl_test.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Boundary Layers Test.
**
** Application flag:   BL_TEST
** Input scripts:      ocean_bl_test.in
**                     stations_bl_test.in
*/

#define ROMS_MODEL
#define UV_ADV
#define UV_SADVECTION
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define AVERAGES
#define STATIONS
#define SOLVE3D

#define RADIATION_2D

#undef  MY25_MIXING
#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
# define RI_SPLINES
#endif

#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
# define LMD_DDMIX
# define RI_SPLINES
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

#undef SG_BBL
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
#define ALBEDO_CLOUD
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
