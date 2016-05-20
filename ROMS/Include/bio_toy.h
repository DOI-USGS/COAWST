/*
** svn $Id: bio_toy.h 440 2010-01-25 06:36:07Z arango $
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options one-dimensional (vertical) Biology Toy.
**
** Application flag:   BIO_TOY
** Input script:       ocean_bio_toy.in
**                     bioFennel.in, ecosim.in, npzd_Franks.in, npzd_Powell.in
*/

#define ROMS_MODEL
#define UV_ADV
#define UV_SADVECTION
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define AVERAGES
#define SOLVE3D

#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif

#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE
# define ANA_RAIN
#else
# define ANA_SMFLUX
# define ANA_STFLUX
#endif

#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

/*
**  Biological model options.
*/

#define BIO_FENNEL
#undef  ECOSIM
#undef  NEMURO
#undef  NPZD_FRANKS
#undef  NPZD_IRON
#undef  NPZD_POWELL

#ifdef BIO_FENNEL
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
#endif

#if defined ECOSIM || defined BIO_FENNEL
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_CLOUD
#endif

#if defined NEMURO
# define HOLLING_GRAZING
# undef  IVLEV_EXPLICIT
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif

#if defined NPZD_FRANKS || defined NPZD_POWELL
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif

#if defined NPZD_IRON
# define ANA_SPFLUX
# define ANA_BPFLUX
# undef  IRON_LIMIT
# undef  IRON_RELAX
#endif

#if defined BULK_FLUXES || defined ECOSIM
# define ANA_CLOUD
# define PAPA_CLM
#endif
