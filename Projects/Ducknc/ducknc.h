/*
** svn $Id: shoreface.h 139 2008-01-10 00:17:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Hole7s.
**
** Application flag:   HOLE7S
** Input scripts:      hole7s.h
**                     sediment_shoreface.h
*/

#define UV_VIS2
#define MIX_S_UV
#undef DIAGNOSTICS_UV
#undef  AVERAGES
#undef  AVERAGES_NEARSHORE
#undef  AVERAGES_BEDLOAD
#define WET_DRY
#define NEARSHORE_MELLOR
#define SWAN_COUPLING
#ifdef  SWAN_COUPLING
# define MCT_LIB
# define UV_KIRBY
#endif

#define OUT_DOUBLE
#define UV_ADV
#define TS_MPDATA
#define DJ_GRADPS
#undef  SALINITY
#define SOLVE3D
#define SPLINES

/* Boundary Conditions */
/*#define WAVE_BOUNDARY
  #define NS_PERIODIC*/

#define NORTH_FSGRADIENT
#define NORTH_M2GRADIENT
#define NORTH_M3GRADIENT
#define SOUTH_FSGRADIENT
#define SOUTH_M2GRADIENT
#define SOUTH_M3GRADIENT

#define EAST_FSGRADIENT
#define EAST_M2REDUCED
#define FSOBC_REDUCED
/*#define EAST_M2GRADIENT*/
#define EAST_M3GRADIENT
#define EAST_TGRADIENT
#define WESTERN_WALL

/*Define Roller Model*/
#define SVENDSEN_ROLLER

#define MASKING

#undef  ANA_GRID
#undef  ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#undef  ANA_SMFLUX
#undef  UV_LOGDRAG
#undef UV_QDRAG

#ifdef SOLVE3D
#define SSW_BBL
# ifdef SSW_BBL
# define SSW_CALC_ZNOT
# undef  SSW_ZORIP 
# undef  SSW_LOGINT
# endif

# undef  SEDIMENT
# ifdef SEDIMENT
#  undef  SED_MORPH
#  define SUSPLOAD
#  define BEDLOAD_MPM
#  undef  BEDLOAD_SOULSBY
# endif
/*
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define ANA_SEDIMENT
# endif
*/
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX


# undef  ANA_VMIX
# define GLS_MIXING
# if defined GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# undef  CRAIG_BANNER
# undef  CHARNOK
# define  ZOS_HSIG
# define  TKE_WAVEDISS
# endif

#endif
