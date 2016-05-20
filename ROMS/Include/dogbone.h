/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for DOGBONE Application.
**
** DOGBONE - composite grid test
**           Just change Ngrids in the makefile, do not make any changes here.
*/

#define ROMS_MODEL
#define NESTING
#define SOLVE3D

#define UV_ADV
#define UV_QDRAG

#define ANA_SMFLUX
#define MASKING

#ifdef SOLVE3D
# define DJ_GRADPS
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# define SALINITY
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define RI_SPLINES
# endif
#endif

#define OUT_DOUBLE

