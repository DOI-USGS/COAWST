/*****
Global params file for InWave
***/

/*   Master list of all InWave cpp options

/********************************************************************************************************/
/***********************************   USER DEFINED OPTIONS   *******************************************/
/********************************************************************************************************/

*# define INWAVE_MODEL       use to turn ON or OFF InWave model
*# define INWAVE_SWAN_COUPLING       use to turn ON or OFF InWave SWAN coupling
*# define DOPPLER               use to turn ON or OFF the effect of currents on the dispersion relation
*# define ACX_ADVECTION         use to turn ON or OFF advection of Ac in the xi direction
*# define ACY_ADVECTION         use to turn ON or OFF advection of Ac in the etai direction
*# define ACT_ADVECTION         use to turn ON or OFF advection of Ac in the directional direction
*# define ENERGY_DISSIPATION    use to turn ON or OFF energy dissipation
*# define ROELVINK              use to turn ON or OFF Roelvink energy dissipation


/*****************************   BOUNDARY CONDITION OPTIONS   *******************************************/


/******** AC IN X-Y SPACE ********/

*# undef WEST_AC_CLAMPED        Western edge, clamped condition for Ac
*# undef WEST_AC_GRADIENT       Western edge, gradient condition for Ac
*# undef WEST_AC_RADIATION      Western edge, radiation condition for Ac
*# undef WEST_AC_WALL           Western edge, wall condition for Ac

*# undef EAST_AC_CLAMPED        Easthern edge, clamped condition for Ac
*# undef EAST_AC_GRADIENT       Easthern edge, gradient condition for Ac
*# undef EAST_AC_RADIATION      Easthern edge, radiation condition for Ac
*# undef EAST_AC_WALL           Easthern edge, wall condition for Ac

*# undef NORTH_AC_CLAMPED        Northern edge, clamped condition for Ac
*# undef NORTH_AC_GRADIENT       Northern edge, gradient condition for Ac
*# undef NORTH_AC_RADIATION      Northern edge, radiation condition for Ac
*# undef NORTH_AC_WALL           Northern edge, wall condition for Ac

*# undef SOUTH_AC_CLAMPED        Southern edge, clamped condition for Ac
*# undef SOUTH_AC_GRADIENT       Southern edge, gradient condition for Ac
*# undef SOUTH_AC_RADIATION      Southern edge, radiation condition for Ac
*# undef SOUTH_AC_WALL           Southern edge, wall condition for Ac

*# undef EW_AC_PERIODIC
*# undef NS_AC_PERIODIC

/******** AC IN THETA SPACE ********/

*# undef THETA_AC_GRADIENT       Start edge, gradient condition for the Ac directional change
*# undef THETA_AC_WALL           Start edge, wall condition for the Ac directional change
*# undef THETA_AC_PERIODIC


/********************************************************************************************************/
/*****************************   DO NOT MODIFY: NO USER DEFINED    **************************************/
/********************************************************************************************************/


/*****************************   LATERAL BOUNDARY CONDITION   *******************************************/

*# undef EW_C_PERIODIC           East-West periodic boundaries for the wave group celerities
*# undef NS_C_PERIODIC           North-South periodic boundaries for the wave group celerities

*# undef WEST_CX_GRADIENT        Western edge, gradient condition for the xi component of the group velocity
*# undef WEST_CY_GRADIENT        Western edge, gradient condition for the etai component of the group velocity
*# undef WEST_CT_GRADIENT        Western edge, gradient condition for the directional component of the group velocity

*# undef EAST_CX_GRADIENT        Eastern edge, gradient condition for the xi component of the group velocity
*# undef EAST_CY_GRADIENT        Eastern edge, gradient condition for the etai component of the group velocity
*# undef EAST_CT_GRADIENT        Eastern edge, gradient condition for the directional component of the group velocity

*# undef NORTH_CX_GRADIENT       Northern edge, gradient condition for the xi component of the group velocity
*# undef NORTH_CY_GRADIENT       Northern edge, gradient condition for the etai component of the group velocity
*# undef NORTH_CT_GRADIENT       Northern edge, gradient condition for the directional component of the group velocity

*# undef SOUTH_CX_GRADIENT        Southern edge, gradient condition for the xi component of the group velocity
*# undef SOUTH_CY_GRADIENT       Southern edge, gradient condition for the etai component of the group velocity
*# undef SOUTH_CT_GRADIENT       Southern edge, gradient condition for the directional component of the group velocity

*# undef WEST_CX_WALL            Western edge, wall condition for the xi component of the group velocity
*# undef EAST_CX_WALL            Eastern edge, wall condition for the xi component of the group velocity
*# undef NORTH_CY_WALL            Northern edge, wall condition for the xi component of the group velocity
*# undef SOUTH_CY_WALL            Southern edge, wall condition for the xi component of the group velocity

*# undef WEST_AC_OBC       Western edge, open boundary condition
*# undef EAST_AC_OBC       Eastern edge, open boundary condition
*# undef NORTH_AC_OBC       Northern edge, open boundary condition
*# undef SOUTH_AC_OBC       Southern edge, open boundary condition

*/

/********************************************************************************************************/
/**********************************   DO NOT MODIFY: SET DEFAULT   **************************************/
/********************************************************************************************************/

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
#  define WEST_CX_GRADIENT        Western edge, gradient condition for the xi component of the group velocity
# endif
# if !defined WEST_CY_WALL && !defined EW_C_PERIODIC
#  define WEST_CY_GRADIENT        Western edge, gradient condition for the etai component of the group velocity
# endif
# if !defined WEST_CT_WALL && !defined EW_C_PERIODIC
#  define WEST_CT_GRADIENT        Western edge, gradient condition for the directional component of the group velocity
# endif

# if !defined EAST_CX_WALL && !defined EW_C_PERIODIC
#  define EAST_CX_GRADIENT        Eastern edge, gradient condition for the xi component of the group velocity
# endif
# if !defined EAST_CY_WALL && !defined EW_C_PERIODIC
#  define EAST_CY_GRADIENT        Eastern edge, gradient condition for the etai component of the group velocity
# endif
# if !defined EAST_CT_WALL && !defined EW_C_PERIODIC
#  define EAST_CT_GRADIENT        Eastern edge, gradient condition for the directional component of the group velocity
# endif

# if !defined NORTH_CX_WALL && !defined NS_C_PERIODIC
#  define NORTH_CX_GRADIENT       Northern edge, gradient condition for the xi component of the group velocity
# endif
# if !defined NORTH_CY_WALL && !defined NS_C_PERIODIC
#  define NORTH_CY_GRADIENT       Northern edge, gradient condition for the etai component of the group velocity
# endif
# if !defined NORTH_CT_WALL && !defined NS_C_PERIODIC
#  define NORTH_CT_GRADIENT       Northern edge, gradient condition for the directional component of the group velocity
# endif

# if !defined SOUTH_CX_WALL && !defined NS_C_PERIODIC
#  define SOUTH_CX_GRADIENT        Southern edge, gradient condition for the xi component of the group velocity
# endif
# if !defined SOUTH_CY_WALL && !defined NS_C_PERIODIC
#  define SOUTH_CY_GRADIENT       Southern edge, gradient condition for the etai component of the group velocity
# endif
# if !defined SOUTH_CT_WALL && !defined NS_C_PERIODIC
#  define SOUTH_CT_GRADIENT       Southern edge, gradient condition for the directional component of the group velocity
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
