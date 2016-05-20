/*
** svn $Id: sed_toy.h 2232 2012-01-03 18:55:20Z arango $
*******************************************************************************
** Copyright (c) 2002-2012 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for One-Dimensional (vertical) Sediment Toy.
**
** Application flag:   SEDBED_TOY
** Input scripts:      ocean_sed_toy.in
**                     sediment_sed_toy.in
*/

#define ROMS_MODEL

#undef  BODYFORCE
#undef  LOG_PROFILE
#define DJ_GRADPS
#undef  TS_U3HADVECTION
#undef  TS_C2VADVECTION
#define TS_MPDATA
#undef  SALINITY
#define SPLINES_VVISC
#define SPLINES_VDIFF
#define OUT_DOUBLE

#define ANA_GRID
#undef  ANA_INITIAL
#undef  ANA_SMFLUX
#define SOLVE3D
#ifdef SOLVE3D
# undef ANA_SEDIMENT
# define ANA_BPFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_STFLUX
#endif
#undef  ANA_VMIX
#undef  ANA_WWAVE

/* select one of six bottom stress methods */
#define UV_LOGDRAG
#undef  UV_LDRAG
#undef  UV_QDRAG
#undef  SG_BBL
#undef  MB_BBL
#undef SSW_BBL

#ifdef SG_BBL
# undef  SG_CALC_ZNOT
# undef  SG_LOGINT
#endif
#ifdef MB_BBL
# undef  MB_CALC_ZNOT
# undef  MB_Z0BIO
# undef  MB_Z0BL
# undef  MB_Z0RIP
#endif
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# undef  SSW_LOGINT
#endif

/* turb closure */
#define GLS_MIXING
#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define RI_SPLINES
# undef  CRAIG_BANNER
# undef  CHARNOK
# undef  ZOS_HSIG
# undef  TKE_WAVEDISS
#endif

/* sediment choices */
#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
# undef  BEDLOAD_SOULSBY
# undef  BEDLOAD_MPM
# define SED_DENS
# define SED_MORPH
# undef  SED_BIODIFF
# undef  NONCOHESIVE_BED1
# undef  NONCOHESIVE_BED2
# define COHESIVE_BED
# undef  MIXED_BED
#endif

