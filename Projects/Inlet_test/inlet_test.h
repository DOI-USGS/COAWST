/*
** svn $Id: inlet_test.h 838 2008-11-17 04:22:18Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Inlet Test Case, waves-ocean (SWAN/ROMS) two-way coupling.
**
** Application flag:   INLET_TEST
** Input script:       ocean_inlet_test.in
**                     coupling_inlet_test.in
**                     sediment_inlet_test.in
*/

#define UV_VIS2
#define MIX_S_UV
#undef  ANA_GRID
#define MASKING
#define UV_ADV
#undef  UV_COR
#define TS_MPDATA
#define DJ_GRADPS
#define SOUTHERN_WALL
#define FSOBC_REDUCED
#define NORTH_FSGRADIENT
#define NORTH_M2REDUCED
#define NORTH_M3GRADIENT
#define WEST_FSGRADIENT
#define WEST_M2GRADIENT
#define WEST_M3GRADIENT
#define EAST_FSGRADIENT
#define EAST_M2GRADIENT
#define EAST_M3GRADIENT
#define SOLVE3D
#define SPLINES
#define SWAN_COUPLING
#define NEARSHORE_MELLOR
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_FSOBC
#define ANA_M2OBC

#ifdef SWAN_COUPLING
# define MCT_LIB
# undef MCT_INTERP_OC2WV
#endif
#undef REFINED_GRID

/* define only one of the following 5 */
#undef  UV_LOGDRAG
#undef  UV_QDRAG
#undef  MB_BBL
#undef  SG_BBL
#define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
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
#  undef  BEDLOAD_SOULSBY
#  undef  BEDLOAD_MPM
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

