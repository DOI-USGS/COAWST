/*
** svn $Id: riverplume1.h 1001 2020-01-10 22:41:16Z arango $
*******************************************************************************
** Copyright (c) 2002-2020 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for River Plume Test (original version).
**
** Application flag:   RIVERPLUME1
** Input script:       roms_riverplume1.in
*/

#define ROMS_MODEL
#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_DIF2
#define MIX_GEO_TS
#define NONLIN_EOS
#define SALINITY
#define MASKING
#define SOLVE3D
#define AVERAGES

#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif

#define ANA_GRID
#define ANA_INITIAL
#define ANA_PSOURCE
#define ANA_SMFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
