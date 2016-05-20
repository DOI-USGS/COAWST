/*
** svn $Id: inlet_test.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
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

#define ROMS_MODEL
#define UV_VIS2
#define MIX_S_UV
#undef  ANA_GRID
#define MASKING
#define UV_ADV
#undef  UV_COR
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_MPDATA
#define DJ_GRADPS
#define FSOBC_REDUCED
#define SOLVE3D
#define SWAN_COUPLING
#undef  WEC_MELLOR
#define WEC_VF
#define WDISS_WAVEMOD
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_FSOBC
#define ANA_M2OBC

#ifdef SWAN_COUPLING
# define MCT_LIB
#endif

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
#  define RI_SPLINES
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

