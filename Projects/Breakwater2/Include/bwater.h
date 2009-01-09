/*
** svn $Id: adriatic.h 8 2007-02-06 19:00:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for BREAKWATER Application.
**
**  Waves-ocean (SWAN/ROMS) two-way coupling test.
*/

# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX

#undef  ANA_GRID
#define UV_VIS2
#define MIX_S_UV
#define MASKING
#define UV_ADV
#define UV_COR
#define TS_MPDATA
/*
#define TS_U3HADVECTION
*/
#define DJ_GRADPS

/*
#define RADIATION_2D
#define EAST_FSRADIATION
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define EAST_KRADIATION
#define EAST_TRADIATION
#define WESTERN_WALL
*/


#define EASTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define WESTERN_WALL
#define SOLVE3D
#define SPLINES
#define SWAN_COUPLING
#define NEARSHORE_MELLOR
#undef WET_DRY

#undef UV_QDRAG
#define SSW_BBL
#define SSW_CALC_ZNOT

/* define only one of the following 5 */
#undef MB_BBL
#ifdef MB_BBL
# define MB_CALC_ZNOT
#endif

#ifdef SOLVE3D
# define GLS_MIXING
# ifdef GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif

# define SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  define  BEDLOAD_SOULSBY
#  undef  BEDLOAD_MPM
#  undef  SED_DENS
#  define SED_MORPH
# endif
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define ANA_SEDIMENT
# endif
#endif

