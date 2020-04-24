/*
** svn $Id: benchmark.h 1001 2020-01-10 22:41:16Z arango $
*******************************************************************************
** Copyright (c) 2002-2020 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Benchmark Tests.  There are several grids configurations to run,
** choose the appropriate standard input script: small ("roms_benchmark1.in"),
** medium ("roms_benchmark2.in"), and large (roms_benchmark3.in).
**
** Application flag:   BENCHMARK
** Input script:       benchmark1.in, benchmark2.in, benchmark3.in
*/

#define ROMS_MODEL
#undef PARALLEL_IO
#define NETCDF4
#define AVERAGES

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_DIF2
#define MIX_GEO_TS
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define CURVGRID
#define SOLVE3D
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define ANA_WINDS
# define ANA_TAIR
# define ANA_PAIR
# define ANA_HUMIDITY
# define ANA_RAIN
# define LONGWAVE
# define ANA_CLOUD
#endif
#define SPHERICAL
#define ANA_GRID
#define ANA_INITIAL
#define ALBEDO_CLOUD
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
