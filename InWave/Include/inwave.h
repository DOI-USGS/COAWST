#if defined INWAVE_MODEL

# if defined WEST_AC_WALL
#  define WEST_CX_WALL
# endif
# if defined EAST_AC_WALL
#  define EAST_CX_WALL
# endif
# if defined NORTH_AC_WALL
#  define NORTH_CY_WALL
# endif
# if defined SOUTH_AC_WALL
#  define SOUTH_CY_WALL
# endif

# if defined NS_AC_PERIODIC
#  define NS_C_PERIODIC
# endif
# if defined EW_AC_PERIODIC
#  define EW_C_PERIODIC
# endif

# if !defined WEST_CX_WALL && !defined EW_C_PERIODIC
# define WEST_CX_GRADIENT
# endif
# if !defined WEST_CY_WALL && !defined EW_C_PERIODIC
# define WEST_CY_GRADIENT
# endif
# if !defined WEST_CT_WALL && !defined EW_C_PERIODIC
# define WEST_CT_GRADIENT
# endif

# if !defined EAST_CX_WALL && !defined EW_C_PERIODIC
# define EAST_CX_GRADIENT
# endif
# if !defined EAST_CY_WALL && !defined EW_C_PERIODIC
# define EAST_CY_GRADIENT
# endif
# if !defined EAST_CT_WALL && !defined EW_C_PERIODIC
# define EAST_CT_GRADIENT
# endif

# if !defined NORTH_CX_WALL && !defined NS_C_PERIODIC
# define NORTH_CX_GRADIENT
# endif
# if !defined NORTH_CY_WALL && !defined NS_C_PERIODIC
# define NORTH_CY_GRADIENT
# endif
# if !defined NORTH_CT_WALL && !defined NS_C_PERIODIC
# define NORTH_CT_GRADIENT
# endif

# if !defined SOUTH_CX_WALL && !defined NS_C_PERIODIC
# define SOUTH_CX_GRADIENT
# endif
# if !defined SOUTH_CY_WALL && !defined NS_C_PERIODIC
# define SOUTH_CY_GRADIENT
# endif
# if !defined SOUTH_CT_WALL && !defined NS_C_PERIODIC
# define SOUTH_CT_GRADIENT
# endif

# if defined WEST_AC_CLAMPED || defined EAST_AC_CLAMPED || defined NORTH_AC_CLAMPED || defined SOUTH_AC_CLAMPED
#  define AC_OBC
# endif

# if defined WEST_AC_CLAMPED
#  define WEST_AC_OBC
# endif
# if defined EAST_AC_CLAMPED
#  define EAST_AC_OBC
# endif
# if defined NORTH_AC_CLAMPED
#  define NORTH_AC_OBC
# endif
# if defined SOUTH_AC_CLAMPED
#  define SOUTH_AC_OBC
# endif
# if defined WEST_AC_WALL
#  define WEST_CX_WALL
# endif
# if defined EAST_AC_WALL
#  define EAST_CX_WALL
# endif
# if defined NORTH_AC_WALL
#  define NORTH_CY_WALL
# endif
# if defined SOUTH_AC_WALL
#  define SOUTH_CY_WALL
# endif

#endif
