/*
** svn $Id: vegetation_growth_test.h 838 2019-11-17 04:22:18Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2019 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Vegetation Growth Test Case, ROMS only.
**
** Application flag:   VEG_GROWTH_TEST
** Input script:       ocean_veg_test.in
**                     vegetation.in
**                     sediment_veg_test.in
*/
#define ROMS_MODEL

#define ESTUARYBGC
#define BIO_SEDIMENT
#define SPECTRAL_LIGHT
#define DENITRIFICATION
#define ALGAL_RESP
#define SAV_BIOMASS
#define OXYGEN
#define CARBON
#define CDOM_DEFAULT
#define MOD_SWR_SPECT
#define ANA_BIOLOGY
#define ANA_SPECIR

#define NO_LBC_ATT

#define MASKING
#define WET_DRY
#define OUT_DOUBLE
#define UV_ADV

#define UV_VIS2
#define MIX_S_UV
#define TS_DIF2
#define MIX_S_TS

#define DJ_GRADPS
#define SALINITY
#define SOLVE3D
#define SPLINES_VVISC
#define SPLINES_VDIFF
#define CURVGRID

#define UV_LOGDRAG

#define SOLAR_SOURCE
#define LONGWAVE
#define CLOUDS

#define BULK_FLUXES
#define ANA_LRFLUX

#define ANA_TOBC
#define ANA_HUMIDITY
#define ANA_RAIN
#define ANA_CLOUD
#define ANA_PAIR
#define ANA_TAIR
#define ANA_WINDS

# define ANA_FSOBC
# define ANA_M2OBC
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX

#define VEGETATION
# ifdef VEGETATION
#  define VEG_DRAG
#  ifdef VEG_DRAG
#   define VEG_TURB
#  endif
#  define VEG_BIOMASS
#  endif

#ifdef SOLVE3D
# define GLS_MIXING
# ifdef GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# ifdef VEG_TURB
#  undef N2S2_HORAVG
# endif
#  define RI_SPLINES
# endif

# define  SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD_SOULSBY
#  undef  BEDLOAD_MPM
#  undef  SED_MORPH
# endif
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define ANA_SEDIMENT
# endif
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
#endif

#undef  DIAGNOSTICS_UV
