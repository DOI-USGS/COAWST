/*
** svn $Id: a4dvar_toy.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for 4DVar Data Assimilation Toy
**
** Application flag:   A4DVAR_TOY
** Input script:       ocean_a4dvar_toy.in
**                     s4dvar.in
*/

#undef  NLM_DRIVER
#undef  IS4DVAR
#define W4DPSAS
#undef  W4DVAR
#undef  CORRELATION
#undef  CONVOLUTION
#undef  R_SYMMETRY
#undef  SANITY_CHECK

#define UV_ADV
#define UV_COR
#define UV_LDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_C4HADVECTION
#define TS_C4VADVECTION
#define TS_DIF2
#define MIX_S_TS
#define SALINITY
#define SOLVE3D
#define SPLINES
#define NS_PERIODIC
#define EASTERN_WALL
#define WESTERN_WALL
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_BMFLUX

#ifdef W4DPSAS
# define CONVOLVE
#endif
#ifdef W4DVAR
# define RPM_RELAXATION
# define CONVOLVE
#endif
#ifndef NLM_DRIVER
# define VCONVOLUTION
# define IMPLICIT_VCONV
#endif
#ifndef NLM_DRIVER
# define FORWARD_READ
# define FORWARD_WRITE
#endif
#ifdef NLM_DRIVER
# define VERIFICATION
#endif
#define OUT_DOUBLE
