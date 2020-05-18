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
# undef ExpF
# define ExpG

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
            /* WRF<->ROMS<-> SWAN; WRF<->SWAN: WRF wind to SWAN, SWAN roughness coef to WRF PBL scheme */
#ifdef ExpE /* Joe finds out WRF PBL condition, if it accepts wave info*/
# define WRF_MODEL
# define ROMS_MODEL
# define SWAN_MODEL
# define COARE_TAYLOR_YELLAND
#endif

#ifdef ExpF  /*Same as ExpD, but will SWAN BBL dynamics */
# define WRF_MODEL
# define ROMS_MODEL
# define SWAN_MODEL
# define COARE_TAYLOR_YELLAND
# undef WEC_VF
# define SSW_BBL
#endif

#ifdef ExpG /*Same as ExpD, but will both SWAN BBL dynamics and SWAN radiation stress*/
# define SWAN_MODEL
# define WRF_MODEL
# define ROMS_MODEL
# define COARE_TAYLOR_YELLAND
# undef  WEC_VF
# define SSW_BBL
#endif

#define MCT_LIB

#if defined WRF_MODEL || defined SWAN_MODEL
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
#  /*define ANA_SEDIMENT*/
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
# define SPLINES_VVISC
# define SPLINES_VDIFF
# define AVERAGES
# define NONLIN_EOS

/* Grid and Initial */
# define MASKING
/*#define ANA_GRID
  #define ANA_INITIAL
  #define ANA_MASK*/
# undef NO_MASK_TEMP        /* JBZ 28 Jan 09, undefined */

/* Forcing */
# ifdef WRF_MODEL
#  undef  BULK_FLUXES
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
# define RI_SPLINES
# undef  MY25_MIXING

# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif

/* Output */
# define DIAGNOSTICS_UV
# define DIAGNOSTICS_TS

/* Biological module */
# undef  NPZD_POWELL

# if defined NPZD_POWELL
#  define ANA_BIOLOGY
#  define ANA_SPFLUX
#  define ANA_BPFLUX
#  define ANA_SRFLUX
# endif
#endif
