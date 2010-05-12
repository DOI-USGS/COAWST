#include "cppdefs.h"
      MODULE mod_boundary

#ifdef OBC
!
!svn $Id: mod_boundary.F 793 2008-10-15 13:57:18Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Open boundary conditions arrays:                                    !
!                                                                      !
# ifdef WEST_FSOBC
!  zeta_west      Free-surface (m) western boundary conditions.        !
#  if defined CELERITY_READ || defined CELERITY_WRITE
!  zeta_west_C2   Free-surface western boundary, quadratic radition    !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_west_Cx   Free-surface western boundary, celerity in the       !
!                   XI-direction (m).                                  !
!  zeta_west_Ce   Free-surface western boundary, celerity in the       !
!                   ETA-direction (m).                                 !
#  endif
#  ifndef ANA_FSOBC
!  zetaG_west     Latest two-time snapshots of input free-surface (m)  !
!                   western boundary data.                             !
#  endif
# endif
# ifdef EAST_FSOBC
!  zeta_east      Free-surface (m) eastern boundary conditions.        !
#  if defined CELERITY_READ || defined CELERITY_WRITE
!  zeta_east_C2   Free-surface eastern boundary, quadratic radition    !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_east_Cx   Free-surface eastern boundary, celerity in the       !
!                   XI-direction (m).                                  !
!  zeta_east_Ce   Free-surface eastern boundary, celerity in the       !
!                   ETA-direction (m).                                 !
#  endif
#  ifndef ANA_FSOBC
!  zetaG_east     Latest two-time snapshots of input free-surface (m)  !
!                   eastern boundary data.                             !
#  endif
# endif
# ifdef SOUTH_FSOBC
!  zeta_south     Free-surface (m) southern boundary conditions.       !
#  if defined CELERITY_READ || defined CELERITY_WRITE
!  zeta_south_C2  Free-surface southern boundary, quadratic radition   !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_south_Cx  Free-surface southern boundary, celerity in the      !
!                   XI-direction (m).                                  !
!  zeta_south_Ce  Free-surface northern boundary, celerity in the      !
!                   ETA-direction (m).                                 !
#  endif
#  ifndef ANA_FSOBC
!  zetaG_south    Latest two-time snapshots of input free-surface (m)  !
!                   southern boundary data.                            !
#  endif
# endif
# ifdef NORTH_FSOBC
!  zeta_north     Free-surface (m) northern boundary conditions.       !
#  if defined CELERITY_READ || defined CELERITY_WRITE
!  zeta_north_C2  Free-surface northern boundary, quadratic radition   !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_north_Cx  Free-surface northern boundary, celerity in the      !
!                   XI-direction (m).                                  !
!  zeta_north_Ce  Free-surface northern boundary, celerity in the      !
!                   ETA-direction (m).                                 !
#  endif
#  ifndef ANA_FSOBC
!  zetaG_north    Latest two-time snapshots of input free-surface (m)  !
!                   northern boundary data.                            !
#  endif
# endif
# ifdef WEST_M2OBC
!  ubar_west      2D u-momentum (m/s) western boundary conditions.     !
!  vbar_west      2D v-momentum (m/s) western boundary conditions.     !
#  if defined CELERITY_READ || defined CELERITY_WRITE
!  ubar_west_C2   2D U-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_west_C2   2D V-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_west_Ce   2D U-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!  vbar_west_Ce   2D V-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!  ubar_west_Cx   2D U-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  vbar_west_Cx   2D V-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
#  endif
#  ifndef ANA_M2OBC
!  ubarG_west     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) western boundary data.                       !
!  vbarG_west     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) western boundary data.                       !
#  endif
# endif
# ifdef EAST_M2OBC
!  ubar_east      2D u-momentum (m/s) eastern boundary conditions.     !
!  vbar_east      2D v-momentum (m/s) eastern boundary conditions.     !
#  if defined CELERITY_READ || defined CELERITY_WRITE
!  ubar_east_C2   2D U-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_east_C2   2D V-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_east_Ce   2D U-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!  vbar_east_Ce   2D V-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!  ubar_east_Cx   2D U-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  vbar_east_Cx   2D V-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
#  endif
#  ifndef ANA_M2OBC
!  ubarG_east     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  vbarG_east     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) eastern boundary data.                       !
#  endif
# endif
# ifdef SOUTH_M2OBC
!  ubar_south     2D u-momentum (m/s) southern boundary conditions.    !
!  vbar_south     2D v-momentum (m/s) southern boundary conditions.    !
#  if defined CELERITY_READ || defined CELERITY_WRITE
!  ubar_south_C2  2D U-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_south_C2  2D V-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_south_Ce  2D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!  vbar_south_Ce  2D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!  ubar_south_Cx  2D U-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  vbar_south_Cx  2D V-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
#  endif
#  ifndef ANA_M2OBC
!  ubarG_south    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  vbarG_south    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) southern boundary data.                      !
#  endif
# endif
# ifdef NORTH_M2OBC
!  ubar_north     2D u-momentum (m/s) northern boundary conditions.    !
!  vbar_north     2D v-momentum (m/s) northern boundary conditions.    !
#  if defined CELERITY_READ || defined CELERITY_WRITE
!  ubar_north_C2  2D U-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_north_C2  2D V-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_north_Ce  2D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!  vbar_north_Ce  2D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!  ubar_north_Cx  2D U-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  vbar_north_Cx  2D V-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
#  endif
#  ifndef ANA_M2OBC
!  ubarG_north    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  vbarG_north    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) northern boundary data.                      !
#  endif
# endif
# ifdef SOLVE3D
#  ifdef WEST_M3OBC
!  u_west         3D u-momentum (m/s) western boundary conditions.     !
!  v_west         3D v-momentum (m/s) western boundary conditions.     !
#   if defined CELERITY_READ || defined CELERITY_WRITE
!  u_west_C2      3D U-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_west_C2      3D V-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_west_Ce      3D U-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!  v_west_Ce      3D V-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!  u_west_Cx      3D U-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  v_west_Cx      3D V-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
#   endif
#   ifndef ANA_M3OBC
!  uG_west        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) western boundary data.                       !
!  vG_west        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) western boundary data.                       !
#   endif
#  endif
#  ifdef EAST_M3OBC
!  u_east         3D u-momentum (m/s) eastern boundary conditions.     !
!  v_east         3D v-momentum (m/s) eastern boundary conditions.     !
#   if defined CELERITY_READ || defined CELERITY_WRITE
!  u_east_C2      3D U-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_east_C2      3D V-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_east_Ce      3D U-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!  v_east_Ce      3D V-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!  u_east_Cx      3D U-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  v_east_Cx      3D V-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
#   endif
#   ifndef ANA_M3OBC
!  uG_east        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  vG_east        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) eastern boundary data.                       !
#   endif
#  endif
#  ifdef SOUTH_M3OBC
!  u_south        3D u-momentum (m/s) southern boundary conditions.    !
!  v_south        3D v-momentum (m/s) southern boundary conditions.    !
#   if defined CELERITY_READ || defined CELERITY_WRITE
!  u_south_C2     3D U-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_south_C2     3D V-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_south_Ce     3D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!  v_south_Ce     3D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!  u_south_Cx     3D U-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  v_south_Cx     3D V-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
#   endif
#   ifndef ANA_M3OBC
!  uG_south       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  vG_south       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) southern boundary data.                      !
#   endif
#  endif
#  ifdef NORTH_M3OBC
!  u_north        3D u-momentum (m/s) northern boundary conditions.    !
!  v_north        3D v-momentum (m/s) northern boundary conditions.    !
#   if defined CELERITY_READ || defined CELERITY_WRITE
!  u_north_C2     3D U-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_north_C2     3D V-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_north_Ce     3D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!  v_north_Ce     3D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!  u_north_Cx     3D U-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  v_north_Cx     3D V-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
#   endif
#   ifndef ANA_M3OBC
!  uG_north       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  vG_north       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) northern boundary data.                      !
#   endif
#  endif
#  ifdef WEST_TOBC
!  t_west         Tracer (T units) western boundary conditions.        !
#   if defined CELERITY_READ || defined CELERITY_WRITE
!  t_west_C2      Tracers western boundary, quadratic radition term    !
!                   dTdx*dTdx+dTde*dTde (Tunits squared).              !
!  t_west_Cx      Tracers western boundary, celerity in the            !
!                   XI-direction (Tunits).                             !
!  t_west_Ce      Tracers western boundary, celerity in the            !
!                   ETA-direction (Tunits).                            !
#   endif
#   ifndef ANA_TOBC
!  tG_west        Latest two-time snapshots of input tracer (Tunits)   !
!                   western boundary data.                             !
#   endif
#  endif
#  ifdef EAST_TOBC
!  t_east         Tracer (T units) eastern boundary conditions.        !
#   if defined CELERITY_READ || defined CELERITY_WRITE
!  t_east_C2      Tracers eastern boundary, quadratic radition term    !
!                   term dTdx*dTdx+dTde*dTde (Tunits squared).         !
!  t_east_Cx      Tracers eastern boundary, celerity in the            !
!                   XI-direction (Tunits).                             !
!  t_east_Ce      Tracers eastern boundary, celerity in the            !
!                   ETA-direction (Tunits).                            !
#   endif
#   ifndef ANA_TOBC
!  tG_east        Latest two-time snapshots of input tracer (Tunits)   !
!                   eastern boundary data.                             !
#   endif
#  endif
#  ifdef SOUTH_TOBC
!  t_south        Tracer (T units) southern boundary conditions.       !
#   if defined CELERITY_READ || defined CELERITY_WRITE
!  t_south_C2     Tracers southern boundary, quadratic radition term   !
!                   term dTdx*dTdx+dTde*dTde (Tunits squared).         !
!  t_south_Cx     Tracers southern boundary, celerity in the           !
!                   XI-direction (Tunits).                             !
!  t_south_Ce     Tracers northern boundary, celerity in the           !
!                   ETA-direction (Tunits).                            !
#   endif
#   ifndef ANA_TOBC
!  tG_south       Latest two-time snapshots of input tracer (Tunits)   !
!                   southern boundary data.                            !
#   endif
#  endif
#  ifdef NORTH_TOBC
!  t_north        Tracer (T units) northern boundary conditions.       !
#   if defined CELERITY_READ || defined CELERITY_WRITE
!  t_north_C2     Tracers northern boundary, quadratic radition term   !
!                   dTdx*dTdx+dTde*dTde (Tunits squared).              !
!  t_north_Cx     Tracers northern boundary, celerity in the           !
!                   XI-direction (Tunits).                             !
!  t_north_Ce     Tracers northern boundary, celerity in the           !
!                   ETA-direction (Tunits).                            !
#   endif
#   ifndef ANA_TOBC
!  tG_north       Latest two-time snapshots of input tracer (Tunits)   !
!                   northern boundary data.                            !
#   endif
#  endif
# endif
# ifdef ADJUST_BOUNDARY
!                                                                      !
!  Boundary arrays for 4DVar data assimilation adjustments:            !
!                                                                      !
#  ifdef SOLVE3D
!  t_obc          Tracers open boundaries values used during 4DVar     !
!                   adjustments, nonlinear model.                      !
!  u_obc          3D U-momentum open boundaries values used during     !
!                   4DVar adjustments, nonlinear model.                !
!  v_obc          3D V-momentum open boundaries values used during     !
!                   4DVar adjustments, nonlinear model.                !
#  endif
!  ubar_obc       2D U-momentum open boundaries values used during     !
!                   4DVar adjustments, nonlinear model.                !
!  vbar_obc       2D V-momentum open boundaries values used during     !
!                   4DVar adjustments, nonlinear model.                !
!  zeta_obc       Free-surface open boundaries values used during      !
!                   4DVar adjustments, nonlinear model.                !
!                                                                      !
#  ifdef SOLVE3D
!  ad_t_obc      Tracers open boundaries values used during 4DVar      !
!                   adjustments, adjoint model.                        !
!  ad_u_obc      3D U-momentum open boundaries values used during      !
!                   4DVar adjustments, adjoint model.                  !
!  ad_v_obc      3D V-momentum open boundaries values used during      !
!                   4DVar adjustments, adjoint model.                  !
#  else
!  ad_ubar_obc   2D U-momentum open boundaries values used during      !
!                   4DVar adjustments, adjoint model.                  !
!  ad_vbar_obc   2D V-momentum open boundaries values used during      !
!                   4DVar adjustments, adjoint model.                  !
#  endif
!  ad_zeta_obc   Free-surface open boundaries values used during       !
!                   4DVar adjustments, adjoint model.                  !
!                                                                      !
#  ifdef SOLVE3D
!  tl_t_obc      Tracers open boundaries values used during 4DVar      !
!                   adjustments, tangent linear model.                 !
!  tl_u_obc      3D U-momentum open boundaries values used during      !
!                   4DVar adjustments, tangent linear model.           !
!  tl_v_obc      3D V-momentum open boundaries values used during      !
!                   4DVar adjustments, tangent linear model.           !
#  else
!  tl_ubar_obc   2D U-momentum open boundaries values used during      !
!                   4DVar adjustments, tangent linear model.           !
!  tl_vbar_obc   2D V-momentum open boundaries values used during      !
!                   4DVar adjustments, tangent linear model.           !
#  endif
!  tl_zeta_obc   Free-surface open boundaries values used during       !
!                   4DVar adjustments, tangent linear model.           !
# endif
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_BOUNDARY
!
!  Nonlinear model state.
!
# ifdef WEST_FSOBC
          real(r8), pointer :: zeta_west(:)
#  if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_west_C2(:)
          real(r8), pointer :: zeta_west_Ce(:)
          real(r8), pointer :: zeta_west_Cx(:)
#  endif
#  ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_west(:,:)
#  endif
# endif

# ifdef EAST_FSOBC
          real(r8), pointer :: zeta_east(:)
#  if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_east_C2(:)
          real(r8), pointer :: zeta_east_Ce(:)
          real(r8), pointer :: zeta_east_Cx(:)
#  endif
#  ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_east(:,:)
#  endif
# endif

# ifdef SOUTH_FSOBC
          real(r8), pointer :: zeta_south(:)
#  if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_south_C2(:)
          real(r8), pointer :: zeta_south_Ce(:)
          real(r8), pointer :: zeta_south_Cx(:)
#  endif
#  ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_south(:,:)
#  endif
# endif

# ifdef NORTH_FSOBC
          real(r8), pointer :: zeta_north(:)
#  if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_north_C2(:)
          real(r8), pointer :: zeta_north_Ce(:)
          real(r8), pointer :: zeta_north_Cx(:)
#  endif
#  ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_north(:,:)
#  endif
# endif

# ifdef WEST_M2OBC
          real(r8), pointer :: ubar_west(:)
          real(r8), pointer :: vbar_west(:)
#  if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_west_C2(:)
          real(r8), pointer :: vbar_west_C2(:)
          real(r8), pointer :: ubar_west_Ce(:)
          real(r8), pointer :: vbar_west_Ce(:)
          real(r8), pointer :: ubar_west_Cx(:)
          real(r8), pointer :: vbar_west_Cx(:)
#  endif
#  ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_west(:,:)
          real(r8), pointer :: vbarG_west(:,:)
#  endif
# endif

# ifdef EAST_M2OBC
          real(r8), pointer :: ubar_east(:)
          real(r8), pointer :: vbar_east(:)
#  if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_east_C2(:)
          real(r8), pointer :: vbar_east_C2(:)
          real(r8), pointer :: ubar_east_Ce(:)
          real(r8), pointer :: vbar_east_Ce(:)
          real(r8), pointer :: ubar_east_Cx(:)
          real(r8), pointer :: vbar_east_Cx(:)
#  endif
#  ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_east(:,:)
          real(r8), pointer :: vbarG_east(:,:)
#  endif
# endif

# ifdef SOUTH_M2OBC
          real(r8), pointer :: ubar_south(:)
          real(r8), pointer :: vbar_south(:)
#  if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_south_C2(:)
          real(r8), pointer :: vbar_south_C2(:)
          real(r8), pointer :: ubar_south_Ce(:)
          real(r8), pointer :: vbar_south_Ce(:)
          real(r8), pointer :: ubar_south_Cx(:)
          real(r8), pointer :: vbar_south_Cx(:)
#  endif
#  ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_south(:,:)
          real(r8), pointer :: vbarG_south(:,:)
#  endif
# endif

# ifdef NORTH_M2OBC
          real(r8), pointer :: ubar_north(:)
          real(r8), pointer :: vbar_north(:)
#  if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_north_C2(:)
          real(r8), pointer :: vbar_north_C2(:)
          real(r8), pointer :: ubar_north_Ce(:)
          real(r8), pointer :: vbar_north_Ce(:)
          real(r8), pointer :: ubar_north_Cx(:)
          real(r8), pointer :: vbar_north_Cx(:)
#  endif
#  ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_north(:,:)
          real(r8), pointer :: vbarG_north(:,:)
#  endif
# endif

# ifdef SOLVE3D
#  ifdef WEST_M3OBC
          real(r8), pointer :: u_west(:,:)
          real(r8), pointer :: v_west(:,:)
#   if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_west_C2(:,:)
          real(r8), pointer :: v_west_C2(:,:)
          real(r8), pointer :: u_west_Ce(:,:)
          real(r8), pointer :: v_west_Ce(:,:)
          real(r8), pointer :: u_west_Cx(:,:)
          real(r8), pointer :: v_west_Cx(:,:)
#   endif
#   ifndef ANA_M3OBC
          real(r8), pointer :: uG_west(:,:,:)
          real(r8), pointer :: vG_west(:,:,:)
#   endif
#  endif

#  ifdef EAST_M3OBC
          real(r8), pointer :: u_east(:,:)
          real(r8), pointer :: v_east(:,:)
#   if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_east_C2(:,:)
          real(r8), pointer :: v_east_C2(:,:)
          real(r8), pointer :: u_east_Ce(:,:)
          real(r8), pointer :: v_east_Ce(:,:)
          real(r8), pointer :: u_east_Cx(:,:)
          real(r8), pointer :: v_east_Cx(:,:)
#   endif
#   ifndef ANA_M3OBC
          real(r8), pointer :: uG_east(:,:,:)
          real(r8), pointer :: vG_east(:,:,:)
#   endif
#  endif

#  ifdef SOUTH_M3OBC
          real(r8), pointer :: u_south(:,:)
          real(r8), pointer :: v_south(:,:)
#   if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_south_C2(:,:)
          real(r8), pointer :: v_south_C2(:,:)
          real(r8), pointer :: u_south_Ce(:,:)
          real(r8), pointer :: v_south_Ce(:,:)
          real(r8), pointer :: u_south_Cx(:,:)
          real(r8), pointer :: v_south_Cx(:,:)
#   endif
#   ifndef ANA_M3OBC
          real(r8), pointer :: uG_south(:,:,:)
          real(r8), pointer :: vG_south(:,:,:)
#   endif
#  endif

#  ifdef NORTH_M3OBC
          real(r8), pointer :: u_north(:,:)
          real(r8), pointer :: v_north(:,:)
#   if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_north_C2(:,:)
          real(r8), pointer :: v_north_C2(:,:)
          real(r8), pointer :: u_north_Ce(:,:)
          real(r8), pointer :: v_north_Ce(:,:)
          real(r8), pointer :: u_north_Cx(:,:)
          real(r8), pointer :: v_north_Cx(:,:)
#   endif
#   ifndef ANA_M3OBC
          real(r8), pointer :: uG_north(:,:,:)
          real(r8), pointer :: vG_north(:,:,:)
#   endif
#  endif

#  ifdef WEST_TOBC
          real(r8), pointer :: t_west(:,:,:)
#   if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_west_C2(:,:,:)
          real(r8), pointer :: t_west_Ce(:,:,:)
          real(r8), pointer :: t_west_Cx(:,:,:)
#   endif
#   ifndef ANA_TOBC
          real(r8), pointer :: tG_west(:,:,:,:)
#   endif
#  endif

#  ifdef EAST_TOBC
          real(r8), pointer :: t_east(:,:,:)
#   if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_east_C2(:,:,:)
          real(r8), pointer :: t_east_Ce(:,:,:)
          real(r8), pointer :: t_east_Cx(:,:,:)
#   endif
#   ifndef ANA_TOBC
          real(r8), pointer :: tG_east(:,:,:,:)
#   endif
#  endif

#  ifdef SOUTH_TOBC
          real(r8), pointer :: t_south(:,:,:)
#   if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_south_C2(:,:,:)
          real(r8), pointer :: t_south_Ce(:,:,:)
          real(r8), pointer :: t_south_Cx(:,:,:)
#   endif
#   ifndef ANA_TOBC
          real(r8), pointer :: tG_south(:,:,:,:)
#   endif
#  endif

#  ifdef NORTH_TOBC
          real(r8), pointer :: t_north(:,:,:)
#   if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_north_C2(:,:,:)
          real(r8), pointer :: t_north_Ce(:,:,:)
          real(r8), pointer :: t_north_Cx(:,:,:)
#   endif
#   ifndef ANA_TOBC
          real(r8), pointer :: tG_north(:,:,:,:)
#   endif
#  endif
# endif

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
#  ifdef WEST_FSOBC
          real(r8), pointer :: tl_zeta_west(:)
#  endif
#  ifdef EAST_FSOBC
          real(r8), pointer :: tl_zeta_east(:)
#  endif
#  ifdef SOUTH_FSOBC
          real(r8), pointer :: tl_zeta_south(:)
#  endif
#  ifdef NORTH_FSOBC
          real(r8), pointer :: tl_zeta_north(:)
#  endif

#  ifdef WEST_M2OBC
          real(r8), pointer :: tl_ubar_west(:)
          real(r8), pointer :: tl_vbar_west(:)
#  endif
#  ifdef EAST_M2OBC
          real(r8), pointer :: tl_ubar_east(:)
          real(r8), pointer :: tl_vbar_east(:)
#  endif
#  ifdef SOUTH_M2OBC
          real(r8), pointer :: tl_ubar_south(:)
          real(r8), pointer :: tl_vbar_south(:)
#  endif
#  ifdef NORTH_M2OBC
          real(r8), pointer :: tl_ubar_north(:)
          real(r8), pointer :: tl_vbar_north(:)
#  endif

#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
          real(r8), pointer :: tl_u_west(:,:)
          real(r8), pointer :: tl_v_west(:,:)
#   endif
#   ifdef EAST_M3OBC
          real(r8), pointer :: tl_u_east(:,:)
          real(r8), pointer :: tl_v_east(:,:)
#   endif
#   ifdef SOUTH_M3OBC
          real(r8), pointer :: tl_u_south(:,:)
          real(r8), pointer :: tl_v_south(:,:)
#   endif
#   ifdef NORTH_M3OBC
          real(r8), pointer :: tl_u_north(:,:)
          real(r8), pointer :: tl_v_north(:,:)
#   endif

#   ifdef WEST_TOBC
          real(r8), pointer :: tl_t_west(:,:,:)
#   endif
#   ifdef EAST_TOBC
          real(r8), pointer :: tl_t_east(:,:,:)
#   endif
#   ifdef SOUTH_TOBC
          real(r8), pointer :: tl_t_south(:,:,:)
#   endif
#   ifdef NORTH_TOBC
          real(r8), pointer :: tl_t_north(:,:,:)
#   endif
#  endif
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
#  ifdef WEST_FSOBC
          real(r8), pointer :: ad_zeta_west(:)
#  endif
#  ifdef EAST_FSOBC
          real(r8), pointer :: ad_zeta_east(:)
#  endif
#  ifdef SOUTH_FSOBC
          real(r8), pointer :: ad_zeta_south(:)
#  endif
#  ifdef NORTH_FSOBC
          real(r8), pointer :: ad_zeta_north(:)
#  endif

#  ifdef WEST_M2OBC
          real(r8), pointer :: ad_ubar_west(:)
          real(r8), pointer :: ad_vbar_west(:)
#  endif
#  ifdef EAST_M2OBC
          real(r8), pointer :: ad_ubar_east(:)
          real(r8), pointer :: ad_vbar_east(:)
#  endif
#  ifdef SOUTH_M2OBC
          real(r8), pointer :: ad_ubar_south(:)
          real(r8), pointer :: ad_vbar_south(:)
#  endif
#  ifdef NORTH_M2OBC
          real(r8), pointer :: ad_ubar_north(:)
          real(r8), pointer :: ad_vbar_north(:)
#  endif

#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
          real(r8), pointer :: ad_u_west(:,:)
          real(r8), pointer :: ad_v_west(:,:)
#   endif
#   ifdef EAST_M3OBC
          real(r8), pointer :: ad_u_east(:,:)
          real(r8), pointer :: ad_v_east(:,:)
#   endif
#   ifdef SOUTH_M3OBC
          real(r8), pointer :: ad_u_south(:,:)
          real(r8), pointer :: ad_v_south(:,:)
#   endif
#   ifdef NORTH_M3OBC
          real(r8), pointer :: ad_u_north(:,:)
          real(r8), pointer :: ad_v_north(:,:)
#   endif

#   ifdef WEST_TOBC
          real(r8), pointer :: ad_t_west(:,:,:)
#   endif
#   ifdef EAST_TOBC
          real(r8), pointer :: ad_t_east(:,:,:)
#   endif
#   ifdef SOUTH_TOBC
          real(r8), pointer :: ad_t_south(:,:,:)
#   endif
#   ifdef NORTH_TOBC
          real(r8), pointer :: ad_t_north(:,:,:)
#   endif
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!  Open boundaries arrays used in 4DVar adjustments.
!
#  ifdef SOLVE3D
          real(r8), pointer :: t_obc(:,:,:,:,:,:)
          real(r8), pointer :: b_t_obc(:,:,:,:)
          real(r8), pointer :: d_t_obc(:,:,:,:,:)
          real(r8), pointer :: e_t_obc(:,:,:,:)
          real(r8), pointer :: ad_t_obc(:,:,:,:,:,:)
          real(r8), pointer :: tl_t_obc(:,:,:,:,:,:)

          real(r8), pointer :: u_obc(:,:,:,:,:)
          real(r8), pointer :: b_u_obc(:,:,:)
          real(r8), pointer :: d_u_obc(:,:,:,:)
          real(r8), pointer :: e_u_obc(:,:,:)
          real(r8), pointer :: ad_u_obc(:,:,:,:,:)
          real(r8), pointer :: tl_u_obc(:,:,:,:,:)

          real(r8), pointer :: v_obc(:,:,:,:,:)
          real(r8), pointer :: b_v_obc(:,:,:)
          real(r8), pointer :: d_v_obc(:,:,:,:)
          real(r8), pointer :: e_v_obc(:,:,:)
          real(r8), pointer :: ad_v_obc(:,:,:,:,:)
          real(r8), pointer :: tl_v_obc(:,:,:,:,:)
#  endif
          real(r8), pointer :: ubar_obc(:,:,:,:)
          real(r8), pointer :: b_ubar_obc(:,:)
          real(r8), pointer :: d_ubar_obc(:,:,:)
          real(r8), pointer :: e_ubar_obc(:,:)
          real(r8), pointer :: ad_ubar_obc(:,:,:,:)
          real(r8), pointer :: tl_ubar_obc(:,:,:,:)

          real(r8), pointer :: vbar_obc(:,:,:,:)
          real(r8), pointer :: b_vbar_obc(:,:)
          real(r8), pointer :: d_vbar_obc(:,:,:)
          real(r8), pointer :: e_vbar_obc(:,:)
          real(r8), pointer :: ad_vbar_obc(:,:,:,:)
          real(r8), pointer :: tl_vbar_obc(:,:,:,:)

          real(r8), pointer :: zeta_obc(:,:,:,:)
          real(r8), pointer :: b_zeta_obc(:,:)
          real(r8), pointer :: d_zeta_obc(:,:,:)
          real(r8), pointer :: e_zeta_obc(:,:)
          real(r8), pointer :: ad_zeta_obc(:,:,:,:)
          real(r8), pointer :: tl_zeta_obc(:,:,:,:)
# endif

        END TYPE T_BOUNDARY

        TYPE (T_BOUNDARY), allocatable ::BOUNDARY(:)

      CONTAINS

      SUBROUTINE allocate_boundary (ng)
!
!=======================================================================
!                                                                      !
!  This routine initializes all variables in the module for all nested !
!  grids.  Currently, there is not parallel tiling in boundary arrays. !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_scalars
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

# ifdef ADJUST_BOUNDARY
!
!  Local variable declarations.
!
      integer :: LBij, UBij
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( BOUNDARY(Ngrids) )
!
!  Nonlinear model state.
!
# ifdef WEST_FSOBC
      allocate ( BOUNDARY(ng) % zeta_west(ETA_DIM) )
#  if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % zeta_west_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_west_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_west_Cx(ETA_DIM) )
#  endif
#  ifndef ANA_FSOBC
      allocate ( BOUNDARY(ng) % zetaG_west(ETA_DIM,2) )
#  endif
# endif

# ifdef EAST_FSOBC
      allocate ( BOUNDARY(ng) % zeta_east(ETA_DIM) )
#  if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % zeta_east_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_east_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_east_Cx(ETA_DIM) )
#  endif
#  ifndef ANA_FSOBC
      allocate ( BOUNDARY(ng) % zetaG_east(ETA_DIM,2) )
#  endif
# endif

# ifdef SOUTH_FSOBC
      allocate ( BOUNDARY(ng) % zeta_south(XI_DIM) )
#  if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % zeta_south_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_south_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_south_Cx(XI_DIM) )
#  endif
#  ifndef ANA_FSOBC
      allocate ( BOUNDARY(ng) % zetaG_south(XI_DIM,2) )
#  endif
# endif

# ifdef NORTH_FSOBC
      allocate ( BOUNDARY(ng) % zeta_north(XI_DIM) )
#  if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % zeta_north_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_north_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_north_Cx(XI_DIM) )
#  endif
#  ifndef ANA_FSOBC
      allocate ( BOUNDARY(ng) % zetaG_north(XI_DIM,2) )
#  endif
# endif

# ifdef WEST_M2OBC
      allocate ( BOUNDARY(ng) % ubar_west(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_west(ETA_DIM) )
#  if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % ubar_west_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_west_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_west_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_west_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_west_Cx(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_west_Cx(ETA_DIM) )
#  endif
#  ifndef ANA_M2OBC
      allocate ( BOUNDARY(ng) % ubarG_west(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_west(ETA_DIM,2) )
#  endif
# endif

# ifdef EAST_M2OBC
      allocate ( BOUNDARY(ng) % ubar_east(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east(ETA_DIM) )
#  if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % ubar_east_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_east_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_east_Cx(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east_Cx(ETA_DIM) )
#  endif
#  ifndef ANA_M2OBC
      allocate ( BOUNDARY(ng) % ubarG_east(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_east(ETA_DIM,2) )
#  endif
# endif

# ifdef SOUTH_M2OBC
      allocate ( BOUNDARY(ng) % ubar_south(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south(XI_DIM) )
#  if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % ubar_south_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_south_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_south_Cx(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south_Cx(XI_DIM) )
#  endif
#  ifndef ANA_M2OBC
      allocate ( BOUNDARY(ng) % ubarG_south(XI_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_south(XI_DIM,2) )
#  endif
# endif

# ifdef NORTH_M2OBC
      allocate ( BOUNDARY(ng) % ubar_north(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north(XI_DIM) )
#  if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % ubar_north_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_north_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_north_Cx(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north_Cx(XI_DIM) )
#  endif
#  ifndef ANA_M2OBC
      allocate ( BOUNDARY(ng) % ubarG_north(XI_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_north(XI_DIM,2) )
#  endif
# endif

# ifdef SOLVE3D
#  ifdef WEST_M3OBC
      allocate ( BOUNDARY(ng) % u_west(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_west(ETA_DIM,N(ng)) )
#   if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % u_west_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_west_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_west_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_west_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_west_Cx(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_west_Cx(ETA_DIM,N(ng)) )
#   endif
#   ifndef ANA_M3OBC
      allocate ( BOUNDARY(ng) % uG_west(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_west(ETA_DIM,N(ng),2) )
#   endif
#  endif

#  ifdef EAST_M3OBC
      allocate ( BOUNDARY(ng) % u_east(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east(ETA_DIM,N(ng)) )
#   if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % u_east_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_east_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_east_Cx(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east_Cx(ETA_DIM,N(ng)) )
#   endif
#   ifndef ANA_M3OBC
      allocate ( BOUNDARY(ng) % uG_east(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_east(ETA_DIM,N(ng),2) )
#   endif
#  endif

#  ifdef SOUTH_M3OBC
      allocate ( BOUNDARY(ng) % u_south(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south(XI_DIM,N(ng)) )
#   if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % u_south_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_south_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_south_Cx(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south_Cx(XI_DIM,N(ng)) )
#   endif
#   ifndef ANA_M3OBC
      allocate ( BOUNDARY(ng) % uG_south(XI_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_south(XI_DIM,N(ng),2) )
#   endif
#  endif

#  ifdef NORTH_M3OBC
      allocate ( BOUNDARY(ng) % u_north(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north(XI_DIM,N(ng)) )
#   if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % u_north_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_north_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_north_Cx(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north_Cx(XI_DIM,N(ng)) )
#   endif
#   ifndef ANA_M3OBC
      allocate ( BOUNDARY(ng) % uG_north(XI_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_north(XI_DIM,N(ng),2) )
#   endif
#  endif

#  ifdef WEST_TOBC
      allocate ( BOUNDARY(ng) % t_west(ETA_DIM,N(ng),NT(ng)) )
#   if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % t_west_C2(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_west_Ce(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_west_Cx(ETA_DIM,N(ng),NT(ng)) )
#   endif
#   ifndef ANA_TOBC
      allocate ( BOUNDARY(ng) % tG_west(ETA_DIM,N(ng),2,NT(ng)) )
#   endif
#  endif

#  ifdef EAST_TOBC
      allocate ( BOUNDARY(ng) % t_east(ETA_DIM,N(ng),NT(ng)) )
#   if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % t_east_C2(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_east_Ce(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_east_Cx(ETA_DIM,N(ng),NT(ng)) )
#   endif
#   ifndef ANA_TOBC
      allocate ( BOUNDARY(ng) % tG_east(ETA_DIM,N(ng),2,NT(ng)) )
#   endif
#  endif

#  ifdef SOUTH_TOBC
      allocate ( BOUNDARY(ng) % t_south(XI_DIM,N(ng),NT(ng)) )
#   if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % t_south_C2(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_south_Ce(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_south_Cx(XI_DIM,N(ng),NT(ng)) )
#   endif
#   ifndef ANA_TOBC
      allocate ( BOUNDARY(ng) % tG_south(XI_DIM,N(ng),2,NT(ng)) )
#   endif
#  endif

#  ifdef NORTH_TOBC
      allocate ( BOUNDARY(ng) % t_north(XI_DIM,N(ng),NT(ng)) )
#   if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % t_north_C2(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_north_Ce(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_north_Cx(XI_DIM,N(ng),NT(ng)) )
#   endif
#   ifndef ANA_TOBC
      allocate ( BOUNDARY(ng) % tG_north(XI_DIM,N(ng),2,NT(ng)) )
#   endif
#  endif
# endif

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
#  ifdef WEST_FSOBC
      allocate ( BOUNDARY(ng) % tl_zeta_west(ETA_DIM) )
#  endif
#  ifdef EAST_FSOBC
      allocate ( BOUNDARY(ng) % tl_zeta_east(ETA_DIM) )
#  endif
#  ifdef SOUTH_FSOBC
      allocate ( BOUNDARY(ng) % tl_zeta_south(XI_DIM) )
#  endif
#  ifdef NORTH_FSOBC
      allocate ( BOUNDARY(ng) % tl_zeta_north(XI_DIM) )
#  endif

#  ifdef WEST_M2OBC
      allocate ( BOUNDARY(ng) % tl_ubar_west(ETA_DIM) )
      allocate ( BOUNDARY(ng) % tl_vbar_west(ETA_DIM) )
#  endif
#  ifdef EAST_M2OBC
      allocate ( BOUNDARY(ng) % tl_ubar_east(ETA_DIM) )
      allocate ( BOUNDARY(ng) % tl_vbar_east(ETA_DIM) )
#  endif
#  ifdef SOUTH_M2OBC
      allocate ( BOUNDARY(ng) % tl_ubar_south(XI_DIM) )
      allocate ( BOUNDARY(ng) % tl_vbar_south(XI_DIM) )
#  endif
#  ifdef NORTH_M2OBC
      allocate ( BOUNDARY(ng) % tl_ubar_north(XI_DIM) )
      allocate ( BOUNDARY(ng) % tl_vbar_north(XI_DIM) )
#  endif

#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      allocate ( BOUNDARY(ng) % tl_u_west(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % tl_v_west(ETA_DIM,N(ng)) )
#   endif
#   ifdef EAST_M3OBC
      allocate ( BOUNDARY(ng) % tl_u_east(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % tl_v_east(ETA_DIM,N(ng)) )
#   endif
#   ifdef SOUTH_M3OBC
      allocate ( BOUNDARY(ng) % tl_u_south(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % tl_v_south(XI_DIM,N(ng)) )
#   endif
#   ifdef NORTH_M3OBC
      allocate ( BOUNDARY(ng) % tl_u_north(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % tl_v_north(XI_DIM,N(ng)) )
#   endif

#   ifdef WEST_TOBC
      allocate ( BOUNDARY(ng) % tl_t_west(ETA_DIM,N(ng),NT(ng)) )
#   endif
#   ifdef EAST_TOBC
      allocate ( BOUNDARY(ng) % tl_t_east(ETA_DIM,N(ng),NT(ng)) )
#   endif
#   ifdef SOUTH_TOBC
      allocate ( BOUNDARY(ng) % tl_t_south(XI_DIM,N(ng),NT(ng)) )
#   endif
#   ifdef NORTH_TOBC
      allocate ( BOUNDARY(ng) % tl_t_north(XI_DIM,N(ng),NT(ng)) )
#   endif
#  endif
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
#  ifdef WEST_FSOBC
      allocate ( BOUNDARY(ng) % ad_zeta_west(ETA_DIM) )
#  endif
#  ifdef EAST_FSOBC
      allocate ( BOUNDARY(ng) % ad_zeta_east(ETA_DIM) )
#  endif
#  ifdef SOUTH_FSOBC
      allocate ( BOUNDARY(ng) % ad_zeta_south(XI_DIM) )
#  endif
#  ifdef NORTH_FSOBC
      allocate ( BOUNDARY(ng) % ad_zeta_north(XI_DIM) )
#  endif

#  ifdef WEST_M2OBC
      allocate ( BOUNDARY(ng) % ad_ubar_west(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ad_vbar_west(ETA_DIM) )
#  endif
#  ifdef EAST_M2OBC
      allocate ( BOUNDARY(ng) % ad_ubar_east(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ad_vbar_east(ETA_DIM) )
#  endif
#  ifdef SOUTH_M2OBC
      allocate ( BOUNDARY(ng) % ad_ubar_south(XI_DIM) )
      allocate ( BOUNDARY(ng) % ad_vbar_south(XI_DIM) )
#  endif
#  ifdef NORTH_M2OBC
      allocate ( BOUNDARY(ng) % ad_ubar_north(XI_DIM) )
      allocate ( BOUNDARY(ng) % ad_vbar_north(XI_DIM) )
#  endif

#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
      allocate ( BOUNDARY(ng) % ad_u_west(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % ad_v_west(ETA_DIM,N(ng)) )
#   endif
#   ifdef EAST_M3OBC
      allocate ( BOUNDARY(ng) % ad_u_east(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % ad_v_east(ETA_DIM,N(ng)) )
#   endif
#   ifdef SOUTH_M3OBC
      allocate ( BOUNDARY(ng) % ad_u_south(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % ad_v_south(XI_DIM,N(ng)) )
#   endif
#   ifdef NORTH_M3OBC
      allocate ( BOUNDARY(ng) % ad_u_north(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % ad_v_north(XI_DIM,N(ng)) )
#   endif

#   ifdef WEST_TOBC
      allocate ( BOUNDARY(ng) % ad_t_west(ETA_DIM,N(ng),NT(ng)) )
#   endif
#   ifdef EAST_TOBC
      allocate ( BOUNDARY(ng) % ad_t_east(ETA_DIM,N(ng),NT(ng)) )
#   endif
#   ifdef SOUTH_TOBC
      allocate ( BOUNDARY(ng) % ad_t_south(XI_DIM,N(ng),NT(ng)) )
#   endif
#   ifdef NORTH_TOBC
      allocate ( BOUNDARY(ng) % ad_t_north(XI_DIM,N(ng),NT(ng)) )
#   endif
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!  Open boundaries arrays used in 4DVar adjustments.
!
#  ifdef SOLVE3D
      allocate ( BOUNDARY(ng) % t_obc(LBij:UBij,N(ng),4,                &
     &                                Nbrec(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % b_t_obc(LBij:UBij,N(ng),4,NT(ng)) )
      allocate ( BOUNDARY(ng) % d_t_obc(LBij:UBij,N(ng),4,              &
     &                                  Nbrec(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % e_t_obc(LBij:UBij,N(ng),4,NT(ng)) )
      allocate ( BOUNDARY(ng) % ad_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tl_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng)) )

      allocate ( BOUNDARY(ng) % u_obc(LBij:UBij,N(ng),                  &
     &                                4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % b_u_obc(LBij:UBij,N(ng),4) )
      allocate ( BOUNDARY(ng) % d_u_obc(LBij:UBij,N(ng),                &
     &                                   4,Nbrec(ng)) )
      allocate ( BOUNDARY(ng) % e_u_obc(LBij:UBij,N(ng),4) )
      allocate ( BOUNDARY(ng) % ad_u_obc(LBij:UBij,N(ng),               &
     &                                   4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % tl_u_obc(LBij:UBij,N(ng),               &
     &                                   4,Nbrec(ng),2) )

      allocate ( BOUNDARY(ng) % v_obc(LBij:UBij,N(ng),                  &
     &                                4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % b_v_obc(LBij:UBij,N(ng),4) )
      allocate ( BOUNDARY(ng) % d_v_obc(LBij:UBij,N(ng),                &
     &                                  4,Nbrec(ng)) )
      allocate ( BOUNDARY(ng) % e_v_obc(LBij:UBij,N(ng),4) )
      allocate ( BOUNDARY(ng) % ad_v_obc(LBij:UBij,N(ng),               &
     &                                   4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % tl_v_obc(LBij:UBij,N(ng),               &
     &                                   4,Nbrec(ng),2) )
#  endif
      allocate ( BOUNDARY(ng) % ubar_obc(LBij:UBij,4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % b_ubar_obc(LBij:UBij,4) )
      allocate ( BOUNDARY(ng) % d_ubar_obc(LBij:UBij,4,Nbrec(ng)) )
      allocate ( BOUNDARY(ng) % e_ubar_obc(LBij:UBij,4) )
      allocate ( BOUNDARY(ng) % ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2) )

      allocate ( BOUNDARY(ng) % vbar_obc(LBij:UBij,4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % b_vbar_obc(LBij:UBij,4) )
      allocate ( BOUNDARY(ng) % d_vbar_obc(LBij:UBij,4,Nbrec(ng)) )
      allocate ( BOUNDARY(ng) % e_vbar_obc(LBij:UBij,4) )
      allocate ( BOUNDARY(ng) % ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2) )

      allocate ( BOUNDARY(ng) % zeta_obc(LBij:UBij,4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % b_zeta_obc(LBij:UBij,4) )
      allocate ( BOUNDARY(ng) % d_zeta_obc(LBij:UBij,4,Nbrec(ng)) )
      allocate ( BOUNDARY(ng) % e_zeta_obc(LBij:UBij,4) )
      allocate ( BOUNDARY(ng) % ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2) )
      allocate ( BOUNDARY(ng) % tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2) )
# endif

      RETURN
      END SUBROUTINE allocate_boundary

      SUBROUTINE initialize_boundary (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      real(r8), parameter :: IniVal = 0.0_r8

# ifndef DISTRIBUTE
#  include "set_bounds.h"
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN

# ifdef WEST_FSOBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % zeta_west = IniVal
#  if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_west_C2 = IniVal
          BOUNDARY(ng) % zeta_west_Ce = IniVal
          BOUNDARY(ng) % zeta_west_Cx = IniVal
#  endif
#  ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_west = IniVal
#  endif
        END IF
# endif

# ifdef EAST_FSOBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % zeta_east = IniVal
#  if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_east_C2 = IniVal
          BOUNDARY(ng) % zeta_east_Ce = IniVal
          BOUNDARY(ng) % zeta_east_Cx = IniVal
#  endif
#  ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_east = IniVal
#  endif
        END IF
# endif

# ifdef SOUTH_FSOBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % zeta_south = IniVal
#  if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_south_C2 = IniVal
          BOUNDARY(ng) % zeta_south_Ce = IniVal
          BOUNDARY(ng) % zeta_south_Cx = IniVal
#  endif
#  ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_south = IniVal
#  endif
        END IF
# endif

# ifdef NORTH_FSOBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % zeta_north = IniVal
#  if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_north_C2 = IniVal
          BOUNDARY(ng) % zeta_north_Ce = IniVal
          BOUNDARY(ng) % zeta_north_Cx = IniVal
#  endif
#  ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_north = IniVal
#  endif
        END IF
# endif

# ifdef WEST_M2OBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % ubar_west = IniVal
          BOUNDARY(ng) % vbar_west = IniVal
#  if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % ubar_west_C2 = IniVal
          BOUNDARY(ng) % vbar_west_C2 = IniVal
          BOUNDARY(ng) % ubar_west_Ce = IniVal
          BOUNDARY(ng) % vbar_west_Ce = IniVal
          BOUNDARY(ng) % ubar_west_Cx = IniVal
          BOUNDARY(ng) % vbar_west_Cx = IniVal
#  endif
#  ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_west = IniVal
          BOUNDARY(ng) % vbarG_west = IniVal
#  endif
        END IF
# endif

# ifdef EAST_M2OBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % ubar_east = IniVal
          BOUNDARY(ng) % vbar_east = IniVal
#  if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % ubar_east_C2 = IniVal
          BOUNDARY(ng) % vbar_east_C2 = IniVal
          BOUNDARY(ng) % ubar_east_Ce = IniVal
          BOUNDARY(ng) % vbar_east_Ce = IniVal
          BOUNDARY(ng) % ubar_east_Cx = IniVal
          BOUNDARY(ng) % vbar_east_Cx = IniVal
#  endif
#  ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_east = IniVal
          BOUNDARY(ng) % vbarG_east = IniVal
#  endif
        END IF
# endif

# ifdef SOUTH_M2OBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % ubar_south = IniVal
          BOUNDARY(ng) % vbar_south = IniVal
#  if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % ubar_south_C2 = IniVal
          BOUNDARY(ng) % vbar_south_C2 = IniVal
          BOUNDARY(ng) % ubar_south_Ce = IniVal
          BOUNDARY(ng) % vbar_south_Ce = IniVal
          BOUNDARY(ng) % ubar_south_Cx = IniVal
          BOUNDARY(ng) % vbar_south_Cx = IniVal
#  endif
#  ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_south = IniVal
          BOUNDARY(ng) % vbarG_south = IniVal
#  endif
        END IF
# endif

# ifdef NORTH_M2OBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % ubar_north = IniVal
          BOUNDARY(ng) % vbar_north = IniVal
#  if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % ubar_north_C2 = IniVal
          BOUNDARY(ng) % vbar_north_C2 = IniVal
          BOUNDARY(ng) % ubar_north_Ce = IniVal
          BOUNDARY(ng) % vbar_north_Ce = IniVal
          BOUNDARY(ng) % ubar_north_Cx = IniVal
          BOUNDARY(ng) % vbar_north_Cx = IniVal
#  endif
#  ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_north = IniVal
          BOUNDARY(ng) % vbarG_north = IniVal
#  endif
        END IF
# endif

# ifdef SOLVE3D
#  ifdef WEST_M3OBC
!
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % u_west = IniVal
          BOUNDARY(ng) % v_west = IniVal
#   if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % u_west_C2 = IniVal
          BOUNDARY(ng) % v_west_C2 = IniVal
          BOUNDARY(ng) % u_west_Ce = IniVal
          BOUNDARY(ng) % v_west_Ce = IniVal
          BOUNDARY(ng) % u_west_Cx = IniVal
          BOUNDARY(ng) % v_west_Cx = IniVal
#   endif
#   ifndef ANA_M3OBC
          BOUNDARY(ng) % uG_west = IniVal
          BOUNDARY(ng) % vG_west = IniVal
#   endif
        END IF
#  endif

#  ifdef EAST_M3OBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % u_east = IniVal
          BOUNDARY(ng) % v_east = IniVal
#   if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % u_east_C2 = IniVal
          BOUNDARY(ng) % v_east_C2 = IniVal
          BOUNDARY(ng) % u_east_Ce = IniVal
          BOUNDARY(ng) % v_east_Ce = IniVal
          BOUNDARY(ng) % u_east_Cx = IniVal
          BOUNDARY(ng) % v_east_Cx = IniVal
#   endif
#   ifndef ANA_M3OBC
          BOUNDARY(ng) % uG_east = IniVal
          BOUNDARY(ng) % vG_east = IniVal
#   endif
      END IF
#  endif

#  ifdef SOUTH_M3OBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % u_south = IniVal
          BOUNDARY(ng) % v_south = IniVal
#   if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % u_south_C2 = IniVal
          BOUNDARY(ng) % v_south_C2 = IniVal
          BOUNDARY(ng) % u_south_Ce = IniVal
          BOUNDARY(ng) % v_south_Ce = IniVal
          BOUNDARY(ng) % u_south_Cx = IniVal
          BOUNDARY(ng) % v_south_Cx = IniVal
#   endif
#   ifndef ANA_M3OBC
          BOUNDARY(ng) % uG_south = IniVal
          BOUNDARY(ng) % vG_south = IniVal
#   endif
        END IF
#  endif

#  ifdef NORTH_M3OBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % u_north = IniVal
          BOUNDARY(ng) % v_north = IniVal
#   if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % u_north_C2 = IniVal
          BOUNDARY(ng) % v_north_C2 = IniVal
          BOUNDARY(ng) % u_north_Ce = IniVal
          BOUNDARY(ng) % v_north_Ce = IniVal
          BOUNDARY(ng) % u_north_Cx = IniVal
          BOUNDARY(ng) % v_north_Cx = IniVal
#   endif
#   ifndef ANA_M3OBC
          BOUNDARY(ng) % uG_north = IniVal
          BOUNDARY(ng) % vG_north = IniVal
#   endif
        END IF
#  endif

#  ifdef WEST_TOBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % t_west = IniVal
#   if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % t_west_C2 = IniVal
          BOUNDARY(ng) % t_west_Ce = IniVal
          BOUNDARY(ng) % t_west_Cx = IniVal
#   endif
#   ifndef ANA_TOBC
          BOUNDARY(ng) % tG_west = IniVal
#   endif
        END IF
#  endif

#  ifdef EAST_TOBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % t_east = IniVal
#   if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % t_east_C2 = IniVal
          BOUNDARY(ng) % t_east_Ce = IniVal
          BOUNDARY(ng) % t_east_Cx = IniVal
#   endif
#   ifndef ANA_TOBC
          BOUNDARY(ng) % tG_east = IniVal
#   endif
        END IF
#  endif

#  ifdef SOUTH_TOBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % t_south = IniVal
#   if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % t_south_C2 = IniVal
          BOUNDARY(ng) % t_south_Ce = IniVal
          BOUNDARY(ng) % t_south_Cx = IniVal
#   endif
#   ifndef ANA_TOBC
          BOUNDARY(ng) % tG_south = IniVal
#   endif
        END IF
#  endif

#  ifdef NORTH_TOBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % t_north = IniVal
#   if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % t_north_C2 = IniVal
          BOUNDARY(ng) % t_north_Ce = IniVal
          BOUNDARY(ng) % t_north_Cx = IniVal
#   endif
#   ifndef ANA_TOBC
          BOUNDARY(ng) % tG_north = IniVal
#   endif
        END IF
#  endif
# endif

      END IF

# if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN
#  ifdef WEST_FSOBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % tl_zeta_west = IniVal
        END IF
#  endif
#  ifdef EAST_FSOBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % tl_zeta_east = IniVal
        END IF
#  endif
#  ifdef SOUTH_FSOBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % tl_zeta_south = IniVal
        END IF
#  endif
#  ifdef NORTH_FSOBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % tl_zeta_north = IniVal
        END IF
#  endif
#  ifdef WEST_M2OBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % tl_ubar_west = IniVal
          BOUNDARY(ng) % tl_vbar_west = IniVal
        END IF
#  endif
#  ifdef EAST_M2OBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % tl_ubar_east = IniVal
          BOUNDARY(ng) % tl_vbar_east = IniVal
        END IF
#  endif
#  ifdef SOUTH_M2OBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % tl_ubar_south = IniVal
          BOUNDARY(ng) % tl_vbar_south = IniVal
        END IF
#  endif
#  ifdef NORTH_M2OBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % tl_ubar_north = IniVal
          BOUNDARY(ng) % tl_vbar_north = IniVal
        END IF
#  endif

#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % tl_u_west = IniVal
          BOUNDARY(ng) % tl_v_west = IniVal
        END IF
#   endif
#   ifdef EAST_M3OBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % tl_u_east = IniVal
          BOUNDARY(ng) % tl_v_east = IniVal
      END IF
#   endif
#   ifdef SOUTH_M3OBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % tl_u_south = IniVal
          BOUNDARY(ng) % tl_v_south = IniVal
        END IF
#   endif
#   ifdef NORTH_M3OBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % tl_u_north = IniVal
          BOUNDARY(ng) % tl_v_north = IniVal
        END IF
#   endif

#   ifdef WEST_TOBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % tl_t_west = IniVal
        END IF
#   endif
#   ifdef EAST_TOBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % tl_t_east = IniVal
        END IF
#   endif
#   ifdef SOUTH_TOBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % tl_t_south = IniVal
        END IF
#   endif
#   ifdef NORTH_TOBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % tl_t_north = IniVal
        END IF
#   endif
#  endif
      END IF
# endif

# ifdef ADJOINT
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN
#  ifdef WEST_FSOBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % ad_zeta_west = IniVal
        END IF
#  endif
#  ifdef EAST_FSOBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % ad_zeta_east = IniVal
        END IF
#  endif
#  ifdef SOUTH_FSOBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % ad_zeta_south = IniVal
        END IF
#  endif
#  ifdef NORTH_FSOBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % ad_zeta_north = IniVal
        END IF
#  endif
#  ifdef WEST_M2OBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % ad_ubar_west = IniVal
          BOUNDARY(ng) % ad_vbar_west = IniVal
        END IF
#  endif
#  ifdef EAST_M2OBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % ad_ubar_east = IniVal
          BOUNDARY(ng) % ad_vbar_east = IniVal
        END IF
#  endif
#  ifdef SOUTH_M2OBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % ad_ubar_south = IniVal
          BOUNDARY(ng) % ad_vbar_south = IniVal
        END IF
#  endif
#  ifdef NORTH_M2OBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % ad_ubar_north = IniVal
          BOUNDARY(ng) % ad_vbar_north = IniVal
        END IF
#  endif

#  ifdef SOLVE3D
#   ifdef WEST_M3OBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % ad_u_west = IniVal
          BOUNDARY(ng) % ad_v_west = IniVal
        END IF
#   endif
#   ifdef EAST_M3OBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % ad_u_east = IniVal
          BOUNDARY(ng) % ad_v_east = IniVal
      END IF
#   endif
#   ifdef SOUTH_M3OBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % ad_u_south = IniVal
          BOUNDARY(ng) % ad_v_south = IniVal
        END IF
#   endif
#   ifdef NORTH_M3OBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % ad_u_north = IniVal
          BOUNDARY(ng) % ad_v_north = IniVal
        END IF
#   endif

#   ifdef WEST_TOBC
        IF (NORTH_WEST_TEST) THEN
          BOUNDARY(ng) % ad_t_west = IniVal
        END IF
#   endif
#   ifdef EAST_TOBC
        IF (SOUTH_EAST_TEST) THEN
          BOUNDARY(ng) % ad_t_east = IniVal
        END IF
#   endif
#   ifdef SOUTH_TOBC
        IF (SOUTH_WEST_TEST) THEN
          BOUNDARY(ng) % ad_t_south = IniVal
        END IF
#   endif
#   ifdef NORTH_TOBC
        IF (NORTH_EAST_TEST) THEN
          BOUNDARY(ng) % ad_t_north = IniVal
        END IF
#   endif
#  endif
      END IF
# endif

# ifdef ADJUST_BOUNDARY
!
!  Open boundaries arrays used in 4DVar adjustments.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        IF (NORTH_EAST_TEST) THEN
#  ifdef SOLVE3D
          BOUNDARY(ng) % b_t_obc = IniVal
          BOUNDARY(ng) % b_u_obc = IniVal
          BOUNDARY(ng) % b_v_obc = IniVal

          BOUNDARY(ng) % d_t_obc = IniVal
          BOUNDARY(ng) % d_u_obc = IniVal
          BOUNDARY(ng) % d_v_obc = IniVal

          BOUNDARY(ng) % e_t_obc = IniVal
          BOUNDARY(ng) % e_u_obc = IniVal
          BOUNDARY(ng) % e_v_obc = IniVal

          BOUNDARY(ng) % t_obc = Inival
          BOUNDARY(ng) % u_obc = IniVal
          BOUNDARY(ng) % v_obc = IniVal
#  endif
          BOUNDARY(ng) % b_ubar_obc = IniVal
          BOUNDARY(ng) % b_vbar_obc = IniVal
          BOUNDARY(ng) % b_zeta_obc = IniVal

          BOUNDARY(ng) % d_ubar_obc = IniVal
          BOUNDARY(ng) % d_vbar_obc = IniVal
          BOUNDARY(ng) % d_zeta_obc = IniVal

          BOUNDARY(ng) % e_ubar_obc = IniVal
          BOUNDARY(ng) % e_vbar_obc = IniVal
          BOUNDARY(ng) % e_zeta_obc = IniVal

          BOUNDARY(ng) % ubar_obc = IniVal
          BOUNDARY(ng) % vbar_obc = IniVal
          BOUNDARY(ng) % zeta_obc = IniVal
        END IF
      END IF

      IF ((model.eq.0).or.(model.eq.iADM)) THEN
        IF (NORTH_EAST_TEST) THEN
#  ifdef SOLVE3D
          BOUNDARY(ng) % ad_t_obc = IniVal
          BOUNDARY(ng) % ad_u_obc = IniVal
          BOUNDARY(ng) % ad_v_obc = IniVal
#  endif
          BOUNDARY(ng) % ad_ubar_obc = IniVal
          BOUNDARY(ng) % ad_vbar_obc = IniVal
          BOUNDARY(ng) % ad_zeta_obc = IniVal
        END IF
      END IF

      IF ((model.eq.0).or.(model.eq.iTLM)) THEN
        IF (NORTH_EAST_TEST) THEN
#  ifdef SOLVE3D
          BOUNDARY(ng) % tl_t_obc = IniVal
          BOUNDARY(ng) % tl_u_obc = IniVal
          BOUNDARY(ng) % tl_v_obc = IniVal
#  endif
          BOUNDARY(ng) % tl_ubar_obc = IniVal
          BOUNDARY(ng) % tl_vbar_obc = IniVal
          BOUNDARY(ng) % tl_zeta_obc = IniVal
        END IF
      END IF
# endif

      RETURN
      END SUBROUTINE initialize_boundary
#endif
      END MODULE mod_boundary
