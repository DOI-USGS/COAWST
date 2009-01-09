/*
** svn $Id: breakwater.h 682 2007-04-24 18:32:25Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for INLET_TEST Application.
**
**  Waves-ocean (SWAN/ROMS) two-way coupling test.
*/

#define UV_VIS2
#define MIX_S_UV
#undef  ANA_GRID
#define MASKING
#define UV_ADV
#undef  UV_COR
#define TS_U3HADVECTION
#define DJ_GRADPS
#define SOUTHERN_WALL
#define EASTERN_WALL
#define WESTERN_WALL
#define NORTH_FSCHAPMAN
#define NORTH_M2REDUCED
#define NORTH_M3GRADIENT
#define ANA_FSOBC
#define ANA_M2OBC
#define FSOBC_REDUCED
#define SOLVE3D
#define SPLINES
#define SWAN_COUPLING
#define NEARSHORE_MELLOR
#define WET_DRY
#define ANA_INITIAL
#define ANA_SMFLUX

/* define only one of the following 5 */
#define MB_BBL
#ifdef MB_BBL
# define MB_CALC_ZNOT
#endif

#ifdef SOLVE3D
# define GLS_MIXING
# ifdef GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
# undef SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD_SOULSBY
#  undef  BEDLOAD_MPM
#  undef  SED_DENS
#  define SED_MORPH
# endif
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define ANA_SEDIMENT
# endif
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
#endif

