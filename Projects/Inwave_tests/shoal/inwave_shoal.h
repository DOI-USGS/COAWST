/*
** svn $Id: inwave_shoal.h 503 2008-01-10 00:11:51Z warner $
*******************************************************************************
**
** Options for Inwave shoal test.
**
*/

#define INWAVE_COUPLING
#ifdef INWAVE_COUPLING
# define DOPPLER
# define ACX_ADVECTION
# define ACY_ADVECTION
# define ACT_ADVECTION

# undef WEST_CX_WALL
# undef WEST_CX_GRADIENT
# undef WEST_CY_WALL
# undef WEST_CY_GRADIENT
# undef WEST_CT_WALL
# undef WEST_CT_GRADIENT
# undef EW_C_PERIODIC
# undef NS_C_PERIODIC

# undef SOUTH_AC_CLAMPED
# undef SOUTH_AC_GRADIENT
# undef SOUTH_AC_RADIATION
# undef SOUTH_AC_WALL

# undef NORTH_AC_CLAMPED
# undef NORTH_AC_GRADIENT
# undef NORTH_AC_RADIATION
# undef NORTH_AC_WALL

# undef EAST_AC_CLAMPED
# undef EAST_AC_GRADIENT
# undef EAST_AC_RADIATION
# undef EAST_AC_WALL

# undef WEST_AC_CLAMPED
# undef WEST_AC_GRADIENT
# undef WEST_AC_RADIATION
# undef WEST_AC_WALL

# undef EW_AC_PERIODIC
# undef NS_AC_PERIODIC
#endif

#define UV_VIS2
#define MIX_S_UV
#define WET_DRY
#undef  NEARSHORE_MELLOR
#define UV_ADV
#define TS_MPDATA
#define DJ_GRADPS
#undef  SALINITY
#define SOLVE3D
#define SPLINES
#define NS_PERIODIC
#define EASTERN_WALL
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3GRADIENT

#undef MASKING
#ifdef MASKING
# define ANA_MASK
#endif
#undef ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define UV_QDRAG

#ifdef SOLVE3D
# undef  SSW_BBL
# ifdef SSW_BBL
#  define SSW_CALC_ZNOT
#  undef  SSW_LOGINT
# endif

# undef SEDIMENT
# ifdef SEDIMENT
#  undef  SED_MORPH
#  define SUSPLOAD
#  define BEDLOAD_MPM
#  undef  BEDLOAD_SOULSBY
# endif
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define ANA_SEDIMENT
# endif

# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
# undef  ANA_VMIX

# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  undef CRAIG_BANNER
#  undef CHARNOK
#  undef ZOS_HSIG
#  undef TKE_WAVEDISS
# endif

#endif
