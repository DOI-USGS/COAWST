/*
** svn $Id: test_chan.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Sediment Test Channel Case.
**
** Application flag:   REFINED_CHAN
** Input scripts:      ocean_refined_chan.in
**                     sediment_refined_chan.in
*/

#define ROMS_MODEL
#define REFINED_GRID

#define WRITE_GRID
#define OUT_DOUBLE

#define  UV_ADV
#undef  SALINITY
#define SOLVE3D
#define MASKING
#define DJ_GRADPS

#define ANA_INITIAL
#define ANA_SMFLUX
#ifdef SOLVE3D
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_SSFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_TOBC
#endif
#define GLS_MIXING
#define KANTHA_CLAYSON
#undef  N2S2_HORAVG   /*keep this off for uniformity in this test case*/
#define UV_LOGDRAG
#define TS_U3HADVECTION

#define EAST_FSCHAPMAN
#define EAST_M2REDUCED
#define FSOBC_REDUCED
#define WEST_FSCHAPMAN
#define WEST_M2GRADIENT
#ifdef SOLVE3D
# define EAST_M3GRADIENT
# define WEST_M3GRADIENT
# define WEST_TCLAMPED
#endif
!#define ANA_M2OBC
#define ANA_FSOBC

#undef  SEDIMENT
#ifdef SEDIMENT
# define ANA_SEDIMENT
# define SUSPLOAD
#endif

