/*
** svn $Id: ducknc.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for DUCK, NC Beach Test Case.
**
** Application flag:   DUCKNC
** Input scripts:      ocean_ducknc.h
**                     sediment_shoreface.h
*/

#define ROMS_MODEL
#define WEC_VF
#define ROLLER_RENIERS
#undef  ROLLER_SVENDSEN
#undef  WEC_MELLOR

#define UV_VIS2
#define MIX_S_UV
#define DIAGNOSTICS_UV
#define AVERAGES
#define AVERAGES_WEC
#define WET_DRY
#define UV_ADV
#define UV_C2ADVECTION
#undef  TS_MPDATA
#define TS_U3HADVECTION
#define DJ_GRADPS
#undef  SALINITY
#define SOLVE3D
#define SPLINES_VVISC
#define SPLINES_VDIFF
#define MASKING

#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#undef  UV_LOGDRAG
#ifdef SOLVE3D
# define SSW_BBL
# ifdef SSW_BBL
#  define SSW_CALC_ZNOT
#  undef  SSW_LOGINT
#  undef  SSW_ZORIP
# endif

# undef SEDIMENT
# ifdef SEDIMENT
#  undef  SED_MORPH
#  define SUSPLOAD
#  define BEDLOAD_MPM
#  undef  BEDLOAD_SOULSBY
#  define AVERAGES_BEDLOAD
# endif

/*# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define ANA_SEDIMENT
# endif
*/

# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
# undef  ANA_VMIX

# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define RI_SPLINES 
#  undef CRAIG_BANNER
#  undef CHARNOK
#  undef ZOS_HSIG
#  undef TKE_WAVEDISS
# endif

# undef LMD_MIXING
# ifdef LMD_MIXING
#  undef LMD_RIMIX
#  undef LMD_CONVEC
#  undef LMD_DDMIX
#  define LMD_SKPP
#  define LMD_BKPP
#  undef LMD_NONLOCAL
# endif

#endif
