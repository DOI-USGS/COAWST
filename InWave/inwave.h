/*****
Global params file for InWave
***/

/*   Master list of all InWave cpp options
*#ifdef INWAVE_COUPLING
*# define DOPPLER
*# define ACX_ADVECTION
*# define ACY_ADVECTION
*# define ACT_ADVECTION

*# undef WEST_CX_GRADIENT
*# undef WEST_CY_GRADIENT
*# undef WEST_CT_GRADIENT

*# undef EW_C_PERIODIC
*# undef NS_C_PERIODIC

*# undef SOUTH_AC_CLAMPED
*# undef SOUTH_AC_GRADIENT
*# undef SOUTH_AC_RADIATION
*# undef SOUTH_AC_WALL

*# undef NORTH_AC_CLAMPED
*# undef NORTH_AC_GRADIENT
*# undef NORTH_AC_RADIATION
*# undef NORTH_AC_WALL

*# undef EAST_AC_CLAMPED
*# undef EAST_AC_GRADIENT
*# undef EAST_AC_RADIATION
*# undef EAST_AC_WALL

*# undef WEST_AC_CLAMPED
*# undef WEST_AC_GRADIENT
*# undef WEST_AC_RADIATION
*# undef WEST_AC_WALL

*# undef EW_AC_PERIODIC
*# undef NS_AC_PERIODIC
*/


#if defined WEST_AC_CLAMPED
# define WEST_AC_OBC
#endif
#if defined EAST_AC_CLAMPED
# define EAST_AC_OBC
#endif
#if defined NORTH_AC_CLAMPED
# define NORTH_AC_OBC
#endif
#if defined SOUTH_AC_CLAMPED
# define SOUTH_AC_OBC
#endif
#if defined NORTH_AC_WALL
# define NORTH_CY_WALL
#endif

*# undef WEST_CT_WALL
*# undef WEST_CY_WALL

