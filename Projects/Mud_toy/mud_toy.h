/*
** svn $Id: mud_toy.h 1162 2007-10-22 02:18:14Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for MUD_TOY Application.
**
** One-dimensional (vertical) Sediment Test.
*/

#undef  BODYFORCE
#undef  LOG_PROFILE
#define DJ_GRADPS
#undef  TS_MPDATA
#define TS_U3HADVECTION
#define SALINITY
#define SPLINES
#define OUT_DOUBLE
/*#undef  EW_PERIODIC
  #define NS_PERIODIC */
#define EASTERN_WALL
#define WESTERN_WALL
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define MASKING
#undef SOLVE3D
#ifdef SOLVE3D
# define ANA_SEDIMENT
# define ANA_BPFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_STFLUX
# define TS_FIXED
# undef  SG_BBL
# ifdef SG_BBL
#  undef  SG_CALC_ZNOT
#  undef  SG_LOGINT
# endif
# undef  MB_BBL
# ifdef MB_BBL
#  undef  MB_CALC_ZNOT
#  undef  MB_Z0BIO
#  undef  MB_Z0BL
#  undef  MB_Z0RIP
# endif
# undef SSW_BBL
# ifdef SSW_BBL
#  define SSW_CALC_ZNOT
#  undef  SSW_LOGINT
# endif
# undef ANA_VMIX
# define GLS_MIXING
# ifdef GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  undef  CRAIG_BANNER
#  undef  CHARNOK
#  undef  ZOS_HSIG
#  undef  TKE_WAVEDISS
# endif
# define SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  undef BEDLOAD
#  undef BEDLOAD_SOULSBY
#  undef BEDLOAD_MPM
#  undef SED_DENS
#  undef SED_MORPH
#  undef COHESIVE_BED
#  undef SED_BIODIFF
# endif
#endif

#undef  ANA_WWAVE
#undef  UV_LOGDRAG
#undef  UV_LDRAG
#define UV_QDRAG

#undef  Q_PSOURCE
#define UV_PSOURCE
#define TS_PSOURCE
#define ANA_PSOURCE
