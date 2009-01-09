/*
** svn $Id: mixed_layer.h 513 2007-04-01 18:34:18Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Mixed Layer Deepening on a rectangular periodic grid with
** thermal stratification and wind stress.
**
** Application flag:   Drifter
** Input script:       ocean_mixed_layer.in
*/

#define DIAGNOSTICS_UV
#define OUT_DOUBLE
#define ANA_GRID
#define SALINITY
#define SOLVE3D
#undef  SPLINES
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define TS_U3HADVECTION
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EW_PERIODIC
#define UV_ADV

#undef  MY25_MIXING
#define GLS_MIXING

#if defined GLS_MIXING
# define KANTHA_CLAYSON
# undef  CANUTO_A
# define N2S2_HORAVG
# define CRAIG_BANNER
# define CHARNOK
# undef  ZOS_HSIG
# undef  TKE_WAVEDISS
# undef  ANA_WWAVE
#endif

#define UV_LOGDRAG
