/*
** svn $Id: test_chan.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Sediment Test Channel Case.
**
** Application flag:   TEST_CHAN
** Input scripts:      ocean_test_chan.in
**                     sediment_test_chan.in
*/

#define ROMS_MODEL
#define WRITE_GRID
#define OUT_DOUBLE
#define ANA_GRID
#define UV_ADV
#undef  SALINITY
#define SOLVE3D
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
#define TS_U3HADVECTION
#define ANA_FSOBC
#define ANA_M2OBC
#define SEDIMENT
#define ANA_SEDIMENT
#define SUSPLOAD
#undef  ANA_VMIX
#define GLS_MIXING
#define KANTHA_CLAYSON
#define N2S2_HORAVG
#define RI_SPLINES
#define UV_LOGDRAG
#undef  PERFECT_RESTART
