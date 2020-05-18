/*
** svn $Id: joe_tc.h 25 2007-04-09 23:43:58Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for JOE_TC Test.
**
** Application flag:   JOE_tc
*/

# undef ExpA1
# undef ExpA
# undef ExpB
# undef ExpC
# undef ExpD
# undef ExpE
# define ExpF

#ifdef ExpA1            /* WRF<->SWAN */
# define WRF_MODEL
# define SWAN_MODEL
#endif

#ifdef ExpA            /* WRF->ROMS */
# define AKLIMIT
# define WRF_MODEL
# define SST_CONST
# define ROMS_MODEL
#endif

#ifdef ExpB  /* WRF<->ROMS */
# define WRF_MODEL
# define ROMS_MODEL
#endif

#ifdef ExpC /* WRF<->ROMS<- SWAN :enhanced surface stress, no currents from ROMS to SWAN*/
# define WRF_MODEL
# define ROMS_MODEL
# define SWAN_MODEL
# define UV_CONST
# define COARE_TAYLOR_YELLAND
#endif

#ifdef ExpD /* WRF<->ROMS<-> SWAN :enhanced surface stress, currents from ROMS to SWAN*/
# define WRF_MODEL
# define ROMS_MODEL
# define SWAN_MODEL
# define COARE_TAYLOR_YELLAND
#endif

#ifdef ExpE  /*Same as ExpD, but will SWAN BBL dynamics */
# define WRF_MODEL
# define ROMS_MODEL
# define SWAN_MODEL
# define COARE_TAYLOR_YELLAND
# define SSW_BBL
#endif

#ifdef ExpF /*Same as ExpE, but will both SWAN BBL dynamics and SWAN radiation stress*/
# define SWAN_MODEL
# define WRF_MODEL
# define ROMS_MODEL
# define COARE_TAYLOR_YELLAND
# define WEC_VF
# define SSW_BBL
#endif

#define MCT_LIB
#if defined WRF_MODEL && defined ROMS_MODEL
# define MCT_INTERP_OC2AT
#endif
#if defined WRF_MODEL && defined SWAN_MODEL
# define MCT_INTERP_WV2AT
#endif

#if defined WRF_MODEL && defined SWAN_MODEL
# define DRAGLIM_DAVIS
#endif

#ifdef ROMS_MODEL
/* Physics + numerics */
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define MIX_S_UV

# ifdef SSW_BBL
#  define SSW_CALC_ZNOT
# /*define ANA_SEDIMENT*/
# else
#  define UV_LOGDRAG
# endif
# if !defined SWAN_MODEL && defined SSW_BBL
#  define ANA_WWAVE
# endif

# define DJ_GRADPS
# define TS_DIF2
# define MIX_GEO_TS

# define SALINITY
# define SOLVE3D
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define AVERAGES
# define NONLIN_EOS

/* Grid and Initial */
# define MASKING

/* Forcing */
# ifdef WRF_MODEL
#  undef BULK_FLUXES
#  define ATM2OCN_FLUXES
#  undef  ANA_SSFLUX
#  undef  LONGWAVE_OUT
# else
#  define ANA_SMFLUX
#  define ANA_STFLUX
#  define ANA_SSFLUX
# endif
# define ATM_PRESS
# define ANA_BTFLUX
# define ANA_BSFLUX
# define EMINUSP

/* Turbulence closure */
# define GLS_MIXING
# undef  MY25_MIXING

# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define RI_SPLINES
# endif

/* Output */
# define DIAGNOSTICS_UV
# define DIAGNOSTICS_TS

#endif
