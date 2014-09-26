
/*
** svn $Id: upwelling.h 25 2007-04-09 23:43:58Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for TRENCH
**
** Application flag:   TRENCH
**
**  Trench migration suspended sediment test.
*/
#define ROMS_MODEL

#undef  LOG_PROFILE
#define UV_ADV
#define UV_LOGDRAG
#define TS_U3HADVECTION
#undef  SALINITY
#define SOLVE3D
#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
# define BEDLOAD_MPM
# define SED_MORPH
#endif
#undef SPLINES
#undef  ANA_VMIX
#define GLS_MIXING
#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# define  N2S2_HORAVG
#endif
#define ANA_BPFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_SMFLUX
#define ANA_SPFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_TOBC
#define ANA_M2OBC

