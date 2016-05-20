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

#define ROMS_MODEL
#define INWAVE_MODEL
#undef  INWAVE_SWAN_COUPLING

#define ACX_ADVECTION
#define ACY_ADVECTION
#define ACT_ADVECTION
#undef DOPPLER

#define WEC_VF
#define ENERGY_DISSIPATION
#define WDISS_INWAVE
#undef ROELVINK

#define UV_VIS2
#define MIX_S_UV
#define MASKING
#define UV_ADV
#undef  UV_COR
#define TS_MPDATA
#define DJ_GRADPS
#define SOLVE3D
#define SPLINES_VVISC
#define SPLINES_VDIFF

#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_FSOBC
#define ANA_M2OBC

/* define only one of the following 5 */
#define UV_LOGDRAG

#ifdef SOLVE3D
# define GLS_MIXING
# ifdef GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define RI_SPLINES 
# endif
# undef SEDIMENT
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

#undef AVERAGES
#undef STATIONS
#undef DIAGNOSTICS_UV
