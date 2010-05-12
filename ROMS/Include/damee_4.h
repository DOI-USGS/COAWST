/*
** svn $Id: damee_4.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for North Atlantic DAMEE Application, 3/4 degree resolution
**
** Application flag:   DAMEE_4
** Input script:       ocean_damee_4.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define NONLIN_EOS
#define SALINITY
#define SOLVE3D
#define MASKING
#define SPLINES
#define QCORRECTION
#define SRELAXATION
#define CURVGRID
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif
#define TCLIMATOLOGY
#define TCLM_NUDGING
#define ANA_BSFLUX
#define ANA_BTFLUX
