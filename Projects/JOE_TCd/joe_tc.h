/*
** svn $Id: upwelling.h 25 2007-04-09 23:43:58Z arango $
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

# undef ExpA
# undef ExpB
# undef ExpC
# undef ExpD
# undef ExpE
# undef ExpF
# undef ExpG
# define ExpH
#define ROMS_MODEL

#ifdef ExpA            /*WRF->ROMS */
# define AKLIMIT
# define WRF_MODEL
# define SST_CONST
#endif

#ifdef ExpB  /*WRF<->ROMS */
# define WRF_MODEL
#endif

#ifdef ExpC /* WRF<->ROMS<- SWAN :enhanced surface stress, no currents from ROMS to SWAN*/
# define WRF_MODEL
# define SWAN_MODEL
# define UV_CONST
# define COARE_TAYLOR_YELLAND
#endif

#ifdef ExpD /* WRF<->ROMS<-> SWAN :enhanced surface stress, currents from ROMS to SWAN*/
# define WRF_MODEL
# define SWAN_MODEL
# define COARE_TAYLOR_YELLAND
#endif
            /* WRF<->ROMS<-> SWAN; WRF<->SWAN: WRF wind to SWAN, SWAN roughness coef to WRF PBL scheme */
#ifdef ExpE /* Joe finds out WRF PBL condition, if it accepts wave info*/
# define WRF_MODEL
# define SWAN_MODEL
# define COARE_TAYLOR_YELLAND
#endif

#ifdef ExpF  /*Same as ExpD, but will SWAN BBL dynamics */
# define WRF_MODEL
# define SWAN_MODEL
# define COARE_TAYLOR_YELLAND
# define SSW_BBL
#endif

#ifdef ExpG /*Same as ExpD, but will both SWAN BBL dynamics and SWAN radiation stress*/
# define WRF_MODEL
# define SWAN_MODEL
# define COARE_TAYLOR_YELLAND
# define WEC_VF
# define SSW_BBL
#endif

#ifdef ExpH /*Same as ExpG, but roms+swan on a different grid than wrf*/
# define WRF_MODEL
# define SWAN_MODEL
# define COARE_TAYLOR_YELLAND
# define WEC_VF
# define SSW_BBL
# define MCT_INTERP_WV2AT
# define MCT_INTERP_OC2AT
#endif


/* Physics + numerics */
#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  TS_MPDATA

#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# /*define ANA_SEDIMENT*/
#else
# define UV_LOGDRAG
#endif
#if !defined SWAN_MODEL && defined SSW_BBL
# define ANA_WWAVE
#endif

#define DJ_GRADPS
#define TS_DIF2
#define MIX_GEO_TS

#define SALINITY
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define NONLIN_EOS

/* Grid and Initial */
#define MASKING
/*#define ANA_GRID
  #define ANA_INITIAL
  #define ANA_MASK*/
#undef NO_MASK_TEMP        /* JBZ 28 Jan 09, undefined */

/* Forcing */
#if defined WRF_MODEL || defined SWAN_MODEL
# define MCT_LIB
#endif
#ifdef WRF_MODEL
# undef BULK_FLUXES
# define ATM2OCN_FLUXES
# define ANA_SSFLUX
# undef LONGWAVE_OUT
#else
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
#endif
# define ATM_PRESS
#define ANA_BTFLUX
#define ANA_BSFLUX

/* Turbulence closure */
#define GLS_MIXING
#undef  MY25_MIXING

#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#endif

/* Boundary condition */
#define WESTERN_WALL
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EASTERN_WALL

/* Output */
#define DIAGNOSTICS_UV
#define DIAGNOSTICS_TS

/* Biological module */
#undef  NPZD_POWELL

#if defined NPZD_POWELL
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_SRFLUX
#endif
