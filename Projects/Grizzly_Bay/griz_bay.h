#define UV_COR
#define UV_ADV
#define MASKING
#define CURVGRID
#define SOLVE3D
#define UV_LOGDRAG
#define WET_DRY

#undef  SALINITY
#undef  NONLIN_EOS
#define SPLINES
#define TS_FIXED
#define ANA_INITIAL
#define ANA_SMFLUX
#ifdef SOLVE3D
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif

#define FSOBC_REDUCED
#define ANA_FSOBC
#define ANA_M2OBC
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3GRADIENT
#define EAST_TGRADIENT
#define WEST_FSCHAPMAN
#define WEST_M2REDUCED
#define WEST_M3GRADIENT
#define WEST_TGRADIENT
