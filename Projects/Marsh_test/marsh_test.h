/*
** svn $Id: marsh_test.h 838 2019-11-17 04:22:18Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2019 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Marsh Dynamics Test Case, ROMS only.
**
** Application flag:   MARSH_TEST
** Input script:       ocean_marsh_test.in
**                     vegetation.in
**                     sediment_marsh_test.in
*/

#define ROMS_MODEL

#define UV_VIS2
#define MIX_S_UV
#define MASKING
#define WET_DRY
#define UV_ADV
#undef  UV_COR
#define DJ_GRADPS
#define SOLVE3D
#define SPLINES_VVISC
#define SPLINES_VDIFF
#undef  SALINITY

#define ANA_SMFLUX
#define ANA_FSOBC
#define ANA_M2OBC

#define WAVES_HEIGHT
#define WAVES_LENGTH

#define ANA_WWAVE

#define VEGETATION
# define VEG_DRAG

# define MARSH_DYNAMICS

# define MARSH_WAVE_THRUST
# define MARSH_SED_EROSION

# define MARSH_TIDAL_RANGE_CALC
# define MARSH_VERT_GROWTH
/** If want internal calculation**/
/** Choose one of the two formulation **/
#  define MARSH_KIRWAN_FORMULATION
#  define MARSH_TIDAL_RANGE_INTERNAL 
#  undef MARSH_MCKEE_FORMULATION

#  define MARSH_BIOMASS_VEG 

#define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
#endif

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
#endif

#  define ANA_SEDIMENT
#define SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  define SED_MORPH
#  define BEDLOAD
#  define  BEDLOAD_SOULSBY
# endif

/*
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define SED_MORPH*/

# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX

