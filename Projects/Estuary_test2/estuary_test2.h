/*
** svn $Id: estuary_test.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Estuary test.
**
** Application flag:   ESTUARY_TEST2
** Input script:       ocean_estuary_test2.in
*/

#define ROMS_MODEL

#define UV_LOGDRAG
#define TS_U3HADVECTION
#define SALINITY
#define SOLVE3D
#define SPLINES_VVISC
#define SPLINES_VDIFF
#define DJ_GRADPS
#define CURVGRID
#define UV_ADV
#undef  UV_COR
#define UV_VIS2
#define MIX_S_UV

#define GLS_MIXING
#if defined GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
#define RI_SPLINES
#endif

#define MASKING
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_TOBC
