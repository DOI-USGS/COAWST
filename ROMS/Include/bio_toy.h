/*
** svn $Id: bio_toy.h 503 2008-01-10 00:11:51Z arango $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options one-dimensional (vertical) Biology Toy.
**
** Application flag:   BIO_TOY
** Input script:       ocean_bio_toy.in
**                     bioFasham.in, ecosim.in, npzd_Franks.in, npzd_Powell.in
*/

#define UV_ADV
#define UV_SADVECTION
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define SPLINES
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS
#define SOLVE3D
#define EW_PERIODIC
#define NS_PERIODIC
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
#endif
#undef BIO_FASHAM
#ifdef BIO_FASHAM
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
#endif
#define ECOSIM
#if defined ECOSIM || defined BIO_FASHAM
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_CLOUD
#endif
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE
# define ANA_RAIN
#else
# define ANA_SMFLUX
# define ANA_STFLUX
#endif
#if defined BULK_FLUXES || defined ECOSIM
# define ANA_CLOUD
# define PAPA_CLM
#endif
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
