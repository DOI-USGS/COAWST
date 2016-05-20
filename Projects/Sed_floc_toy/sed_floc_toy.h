/*
** svn $Id: sed_toy.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for One-Dimensional (vertical) Sediment Toy.
**
** Application flag:   SED_TOY
** Input scripts:      ocean_sed_toy.in
**                     sediment_floc_toy.in
*/

#define ROMS_MODEL
#define BODYFORCE
#undef  LOG_PROFILE
#define DJ_GRADPS
#undef  TS_U3HADVECTION
#undef  TS_C2VADVECTION
#define TS_MPDATA
#undef  SALINITY
#undef  SPLINES
#define OUT_DOUBLE
#define ANA_GRID
#undef  ANA_INITIAL
#define ANA_SMFLUX
#define ATM_PRESS
#undef  ANA_PAIR
#define SOLVE3D
#ifdef SOLVE3D
# undef  ANA_SEDIMENT
# define ANA_BPFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_STFLUX
#endif
#undef   ANA_VMIX
#undef   ANA_WWAVE

#undef  UV_LOGDRAG
#undef  UV_LDRAG
#undef  UV_QDRAG


#undef  SG_BBL
#ifdef SG_BBL
# undef  SG_CALC_ZNOT
# undef  SG_LOGINT
#endif

#undef  MB_BBL
#ifdef MB_BBL
# undef  MB_CALC_ZNOT
# undef  MB_Z0BIO
# undef  MB_Z0BL
# undef  MB_Z0RIP
#endif

#define SSW_BBL
#ifdef SSW_BBL
# undef  SSW_CALC_UB
# define SSW_CALC_ZNOT
# define SSW_LOGINT
#endif

#define GLS_MIXING
#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# undef  CRAIG_BANNER
# undef  CHARNOK
# undef  ZOS_HSIG
# undef  TKE_WAVEDISS
#endif

#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
# undef  BEDLOAD_SOULSBY
# undef  BEDLOAD_MPM
# define SED_DENS
# undef  COHESIVE_BED
# undef  BF_TCR
# undef  LINEAR_TCR
# undef  POWERLAW_TCR
# define MIXED_BED
# undef  SED_MORPH
# define SED_FLOCS
# undef  FLOC_TURB_DISS
# define FLOC_BBL_DISS
# define SED_DEFLOC
# define SED_TAU_CD_CONST
# undef  SED_TAU_CD_LIN
# undef  SED_BIODIFF
#endif
