/*
** svn $Id: sed_toy.h 429 2009-12-20 17:30:26Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for two-dimensional Lentz test.
**
** Application flag:   LENTZ_TEST
** Input scripts:      ocean_sed_toy.in
**                     sediment_sed_toy.in
*/


#define ROMS_MODEL
#undef  WEC_MELLOR
#define WEC_VF
#undef  WDISS_CHURTHOR
#undef  BOTTOM_STREAMING
#undef  SURFACE_STREAMING

#define BODYFORCE
#define UV_ADV
#define UV_COR

#define TS_FIXED
#define DJ_GRADPS
#define SPLINES
#define OUT_DOUBLE
#define NS_PERIODIC
#define EW_PERIODIC

#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_WWAVE
#define SOLVE3D
#ifdef SOLVE3D
# define ANA_BPFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_STFLUX
#endif

#define UV_LOGDRAG

#define GLS_MIXING
#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# undef  CRAIG_BANNER
# undef  CHARNOK
# undef  ZOS_HSIG
# undef  TKE_WAVEDISS
#endif
