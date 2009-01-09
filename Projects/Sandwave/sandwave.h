/*
** svn $Id: sandwave.h 683 2007-04-24 18:38:59Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Sandwave Test Case in a Channel.
**
** Application flag:   SANDWAVE
**                     ocean_sandwave.in
**                     sediment_sandwave.in
*/

#define UV_ADV
#define BODYFORCE
#undef  UV_PSOURCE
#define UV_LOGDRAG
#define UV_VIS4
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  TS_MPDATA
#define TS_DIF4
#define MIX_S_TS
#undef  SALINITY
#undef  NONLIN_EOS
#define SOLVE3D
#undef  SPLINES
#define AVERAGES
#ifdef AVERAGES
# define AVERAGES_AKV
# define AVERAGES_AKT
# define AVERAGES_AKS
# define AVERAGES_BEDLOAD   
#endif
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EW_PERIODIC

#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
# define BEDLOAD_MPM
# undef  BEDLOAD_SOULSBY
# define SED_MORPH
#endif

#undef  MY25_MIXING
#define GLS_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif

#undef  LMD_MIXING
#ifdef  LMD_MIXING
# define LMD_RIMIX
# define LMD_SKPP
# define LMD_BKPP
#endif

#define ANA_BPFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SEDIMENT
#define ANA_SMFLUX
#define ANA_SPFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#undef  ANA_PSOURCE
#define ANA_TOBC
#define ANA_FSOBC
#undef  ANA_VMIX
