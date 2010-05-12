/*
** svn $Id: cblast.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Coupled Boundary Layers and Air-Sea Transfer Application.
**
** Application flag:   CBLAST
** Input script:       ocean_cblast.in
*/

/* Basic physics options */
#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  MIX_S_UV
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS

/* Basic numerics options */
#define UV_SADVECTION
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define DJ_GRADPS
#define SPLINES
#define CURVGRID
#define MASKING

/* Outputs */
#define AVERAGES
#define AVERAGES_QUADRATIC
#define AVERAGES_FLUXES
#define DIAGNOSTICS_UV
#define DIAGNOSTICS_TS
#define STATIONS
#undef  FLOATS

/* Surface and bottom boundary conditions */
#define BULK_FLUXES
#define SOLAR_SOURCE
#define LONGWAVE_OUT /* input is lwrad downward - model computes upward */
#define ANA_RAIN
#define UV_QDRAG
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

/* Vertical subgridscale turbuelnce closure */
#undef  LMD_MIXING
#define MY25_MIXING
#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif
#ifdef  LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define ANA_CLOUD
#endif

/* Open boundary conditions */
#define RADIATION_2D
#define RAMP_TIDES
#define SSH_TIDES
# define ADD_FSOBC
# define EAST_FSCHAPMAN
# define WEST_FSCHAPMAN
# define SOUTH_FSCHAPMAN
# define NORTH_FSCHAPMAN
#define UV_TIDES
# define ADD_M2OBC
# define EAST_M2FLATHER
# define WEST_M2FLATHER
# define SOUTH_M2FLATHER
# define NORTH_M2FLATHER
#define EAST_M3RADIATION
#define EAST_M3NUDGING
#define EAST_TRADIATION
#define EAST_TNUDGING
#define WEST_M3RADIATION
#define WEST_M3NUDGING
#define WEST_TRADIATION
#define WEST_TNUDGING
#define SOUTH_M3RADIATION
#define SOUTH_M3NUDGING
#define SOUTH_TRADIATION
#define SOUTH_TNUDGING
#define NORTH_M3RADIATION
#define NORTH_M3NUDGING
#define NORTH_TRADIATION
#define NORTH_TNUDGING
