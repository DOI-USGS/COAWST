/*
** svn $Id: bevano.h 729 2007-04-27 13:38:15Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Bevano Application.
**
** Test case for closed basin with wind, nearshore processes and sediments.
**
*/

#define WET_DRY
#define SWAN_COUPLING
#define NEARSHORE_MELLOR

#define TS_DIF2
#define MIX_GEO_TS
#define UV_VIS2
#define MIX_S_UV
#define UV_ADV
#define UV_COR
#define DJ_GRADPS

#define TS_MPDATA
#undef  TS_U3HADVECTION
#undef  TS_C4VADVECTION

#undef  OUT_DOUBLE

#define SALINITY
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define AVERAGES_NEARSHORE
#if defined AVERAGES && defined BEDLOAD
# define AVERAGES_BEDLOAD
#endif
#undef FLOATS

#define WESTERN_WALL
#define NORTH_FSRADIATION
#define NORTH_M2RADIATION
#define NORTH_M3RADIATION
#define SOUTH_FSRADIATION
#define SOUTH_M2RADIATION
#define SOUTH_M3RADIATION
#define EAST_FSCHAPMAN
#define EAST_M2REDUCED
#define FSOBC_REDUCED
#define EAST_M3RADIATION
#define ANA_FSOBC

#undef Q_PSOURCE

#undef MB_BBL
#ifdef MB_BBL 
# define MB_CALC_ZNOT
# undef MB_LOGINT
#endif
#define SSW_BBL
#ifdef SSW_BBL 
# define SSW_CALC_ZNOT
# undef SSW_LOGINT
#endif
#undef UV_QDRAG

#define SEDIMENT
#ifdef SEDIMENT
# undef ANA_SEDIMENT
# undef  ANA_WWAVE
# define SED_MORPH
# define SUSPLOAD
# undef  BEDLOAD_MPM
# define BEDLOAD_SOULSBY
#endif

#define ANA_PSOURCE
#define UV_PSOURCE
#define TS_PSOURCE
#define RIVER_SEDIMENT

#define MASKING
#define CURVGRID
#undef ANA_GRID

#undef ANA_SMFLUX
#undef ANA_INITIAL

#ifdef SOLVE3D
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX

/* define one vertical mixing scheme here*/
# undef ANA_VMIX
# undef MY25_MIXING
# define GLS_MIXING
# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define  N2S2_HORAVG
#  undef CRAIG_BANNER
#  undef CHARNOK
#  undef ZOS_HSIG
#  undef TKE_WAVEDISS
# endif
#endif


