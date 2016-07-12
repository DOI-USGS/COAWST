#include "cppdefs.h"
      MODULE mod_boundary
!
!svn $Id: mod_boundary.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Open boundary conditions arrays:                                    !
!                                                                      !
!  zeta_west      Free-surface (m) western boundary conditions.        !
#if defined CELERITY_READ || defined CELERITY_WRITE
!  zeta_west_C2   Free-surface western boundary, quadratic radition    !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_west_Cx   Free-surface western boundary, celerity in the       !
!                   XI-direction (m).                                  !
!  zeta_west_Ce   Free-surface western boundary, celerity in the       !
!                   ETA-direction (m).                                 !
#endif
#ifndef ANA_FSOBC
!  zetaG_west     Latest two-time snapshots of input free-surface (m)  !
!                   western boundary data.                             !
#endif
!  zeta_east      Free-surface (m) eastern boundary conditions.        !
#if defined CELERITY_READ || defined CELERITY_WRITE
!  zeta_east_C2   Free-surface eastern boundary, quadratic radition    !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_east_Cx   Free-surface eastern boundary, celerity in the       !
!                   XI-direction (m).                                  !
!  zeta_east_Ce   Free-surface eastern boundary, celerity in the       !
!                   ETA-direction (m).                                 !
#endif
#ifndef ANA_FSOBC
!  zetaG_east     Latest two-time snapshots of input free-surface (m)  !
!                   eastern boundary data.                             !
#endif
!  zeta_south     Free-surface (m) southern boundary conditions.       !
#if defined CELERITY_READ || defined CELERITY_WRITE
!  zeta_south_C2  Free-surface southern boundary, quadratic radition   !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_south_Cx  Free-surface southern boundary, celerity in the      !
!                   XI-direction (m).                                  !
!  zeta_south_Ce  Free-surface northern boundary, celerity in the      !
!                   ETA-direction (m).                                 !
#endif
#ifndef ANA_FSOBC
!  zetaG_south    Latest two-time snapshots of input free-surface (m)  !
!                   southern boundary data.                            !
#endif
!  zeta_north     Free-surface (m) northern boundary conditions.       !
#if defined CELERITY_READ || defined CELERITY_WRITE
!  zeta_north_C2  Free-surface northern boundary, quadratic radition   !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_north_Cx  Free-surface northern boundary, celerity in the      !
!                   XI-direction (m).                                  !
!  zeta_north_Ce  Free-surface northern boundary, celerity in the      !
!                   ETA-direction (m).                                 !
#endif
#ifndef ANA_FSOBC
!  zetaG_north    Latest two-time snapshots of input free-surface (m)  !
!                   northern boundary data.                            !
#endif
!  ubar_west      2D u-momentum (m/s) western boundary conditions.     !
!  vbar_west      2D v-momentum (m/s) western boundary conditions.     !
#if defined CELERITY_READ || defined CELERITY_WRITE
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
#endif
#ifndef ANA_M2OBC
!  ubarG_west     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) western boundary data.                       !
!  vbarG_west     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) western boundary data.                       !
#endif
!  ubar_east      2D u-momentum (m/s) eastern boundary conditions.     !
!  vbar_east      2D v-momentum (m/s) eastern boundary conditions.     !
#if defined CELERITY_READ || defined CELERITY_WRITE
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
#endif
#ifndef ANA_M2OBC
!  ubarG_east     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  vbarG_east     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) eastern boundary data.                       !
#endif
!  ubar_south     2D u-momentum (m/s) southern boundary conditions.    !
!  vbar_south     2D v-momentum (m/s) southern boundary conditions.    !
#if defined CELERITY_READ || defined CELERITY_WRITE
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
#endif
#ifndef ANA_M2OBC
!  ubarG_south    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  vbarG_south    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) southern boundary data.                      !
#endif
!  ubar_north     2D u-momentum (m/s) northern boundary conditions.    !
!  vbar_north     2D v-momentum (m/s) northern boundary conditions.    !
#if defined CELERITY_READ || defined CELERITY_WRITE
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
#endif
#ifndef ANA_M2OBC
!  ubarG_north    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  vbarG_north    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) northern boundary data.                      !
#endif
#ifdef WEC
!  ubarstokes_east  2D u-momentum (m/s) eastern boundary conditions.   !
!  vbarstokes_east  2D v-momentum (m/s) eastern boundary conditions.   !
!  ubarstokes_west  2D u-momentum (m/s) western boundary conditions.   !
!  vbarstokes_west  2D v-momentum (m/s) western boundary conditions.   !
!  ubarstokes_north 2D u-momentum (m/s) northern boundary conditions.  !
!  vbarstokes_north 2D v-momentum (m/s) northern boundary conditions.  !
!  ubarstokes_south 2D u-momentum (m/s) southern boundary conditions.  !
!  vbarstokes_south 2D v-momentum (m/s) southern boundary conditions.  !
#endif
#ifdef SOLVE3D
!  u_west         3D u-momentum (m/s) western boundary conditions.     !
!  v_west         3D v-momentum (m/s) western boundary conditions.     !
# if defined CELERITY_READ || defined CELERITY_WRITE
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
# endif
# ifndef ANA_M3OBC
!  uG_west        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) western boundary data.                       !
!  vG_west        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) western boundary data.                       !
# endif
!  u_east         3D u-momentum (m/s) eastern boundary conditions.     !
!  v_east         3D v-momentum (m/s) eastern boundary conditions.     !
# if defined CELERITY_READ || defined CELERITY_WRITE
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
# endif
# ifndef ANA_M3OBC
!  uG_east        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  vG_east        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) eastern boundary data.                       !
# endif
!  u_south        3D u-momentum (m/s) southern boundary conditions.    !
!  v_south        3D v-momentum (m/s) southern boundary conditions.    !
# if defined CELERITY_READ || defined CELERITY_WRITE
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
# endif
# ifndef ANA_M3OBC
!  uG_south       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  vG_south       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) southern boundary data.                      !
# endif
!  u_north        3D u-momentum (m/s) northern boundary conditions.    !
!  v_north        3D v-momentum (m/s) northern boundary conditions.    !
# if defined CELERITY_READ || defined CELERITY_WRITE
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
# endif
# ifndef ANA_M3OBC
!  uG_north       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  vG_north       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) northern boundary data.                      !
# endif
# ifdef WEC
!  ustokes_east     3D u-momentum (m/s) eastern boundary conditions.   !
!  vstokes_east     3D v-momentum (m/s) eastern boundary conditions.   !
!  ustokes_west     3D u-momentum (m/s) western boundary conditions.   !
!  vstokes_west     3D v-momentum (m/s) western boundary conditions.   !
!  ustokes_north    3D u-momentum (m/s) northern boundary conditions.  !
!  vstokes_north    3D v-momentum (m/s) northern boundary conditions.  !
!  ustokes_south    3D u-momentum (m/s) southern boundary conditions.  !
!  vstokes_south    3D v-momentum (m/s) southern boundary conditions.  !
# endif
!  t_west         Tracer (T units) western boundary conditions.        !
# if defined CELERITY_READ || defined CELERITY_WRITE
!  t_west_C2      Tracers western boundary, quadratic radition term    !
!                   dTdx*dTdx+dTde*dTde (Tunits squared).              !
!  t_west_Cx      Tracers western boundary, celerity in the            !
!                   XI-direction (Tunits).                             !
!  t_west_Ce      Tracers western boundary, celerity in the            !
!                   ETA-direction (Tunits).                            !
# endif
# ifndef ANA_TOBC
!  tG_west        Latest two-time snapshots of input tracer (Tunits)   !
!                   western boundary data.                             !
# endif
!  t_east         Tracer (T units) eastern boundary conditions.        !
# if defined CELERITY_READ || defined CELERITY_WRITE
!  t_east_C2      Tracers eastern boundary, quadratic radition term    !
!                   term dTdx*dTdx+dTde*dTde (Tunits squared).         !
!  t_east_Cx      Tracers eastern boundary, celerity in the            !
!                   XI-direction (Tunits).                             !
!  t_east_Ce      Tracers eastern boundary, celerity in the            !
!                   ETA-direction (Tunits).                            !
# endif
# ifndef ANA_TOBC
!  tG_east        Latest two-time snapshots of input tracer (Tunits)   !
!                   eastern boundary data.                             !
# endif
!  t_south        Tracer (T units) southern boundary conditions.       !
# if defined CELERITY_READ || defined CELERITY_WRITE
!  t_south_C2     Tracers southern boundary, quadratic radition term   !
!                   term dTdx*dTdx+dTde*dTde (Tunits squared).         !
!  t_south_Cx     Tracers southern boundary, celerity in the           !
!                   XI-direction (Tunits).                             !
!  t_south_Ce     Tracers northern boundary, celerity in the           !
!                   ETA-direction (Tunits).                            !
# endif
# ifndef ANA_TOBC
!  tG_south       Latest two-time snapshots of input tracer (Tunits)   !
!                   southern boundary data.                            !
# endif
!  t_north        Tracer (T units) northern boundary conditions.       !
# if defined CELERITY_READ || defined CELERITY_WRITE
!  t_north_C2     Tracers northern boundary, quadratic radition term   !
!                   dTdx*dTdx+dTde*dTde (Tunits squared).              !
!  t_north_Cx     Tracers northern boundary, celerity in the           !
!                   XI-direction (Tunits).                             !
!  t_north_Ce     Tracers northern boundary, celerity in the           !
!                   ETA-direction (Tunits).                            !
# endif
# ifndef ANA_TOBC
!  tG_north       Latest two-time snapshots of input tracer (Tunits)   !
!                   northern boundary data.                            !
# endif
#endif
#ifdef ADJUST_BOUNDARY
!                                                                      !
!  Boundary arrays for 4DVar data assimilation adjustments:            !
!                                                                      !
# ifdef SOLVE3D
!  t_obc          Tracers open boundaries values used during 4DVar     !
!                   adjustments, nonlinear model.                      !
!  u_obc          3D U-momentum open boundaries values used during     !
!                   4DVar adjustments, nonlinear model.                !
!  v_obc          3D V-momentum open boundaries values used during     !
!                   4DVar adjustments, nonlinear model.                !
# endif
!  ubar_obc       2D U-momentum open boundaries values used during     !
!                   4DVar adjustments, nonlinear model.                !
!  vbar_obc       2D V-momentum open boundaries values used during     !
!                   4DVar adjustments, nonlinear model.                !
!  zeta_obc       Free-surface open boundaries values used during      !
!                   4DVar adjustments, nonlinear model.                !
!                                                                      !
# ifdef SOLVE3D
!  ad_t_obc      Tracers open boundaries values used during 4DVar      !
!                   adjustments, adjoint model.                        !
!  ad_u_obc      3D U-momentum open boundaries values used during      !
!                   4DVar adjustments, adjoint model.                  !
!  ad_v_obc      3D V-momentum open boundaries values used during      !
!                   4DVar adjustments, adjoint model.                  !
# else
!  ad_ubar_obc   2D U-momentum open boundaries values used during      !
!                   4DVar adjustments, adjoint model.                  !
!  ad_vbar_obc   2D V-momentum open boundaries values used during      !
!                   4DVar adjustments, adjoint model.                  !
# endif
!  ad_zeta_obc   Free-surface open boundaries values used during       !
!                   4DVar adjustments, adjoint model.                  !
!                                                                      !
# ifdef SOLVE3D
!  tl_t_obc      Tracers open boundaries values used during 4DVar      !
!                   adjustments, tangent linear model.                 !
!  tl_u_obc      3D U-momentum open boundaries values used during      !
!                   4DVar adjustments, tangent linear model.           !
!  tl_v_obc      3D V-momentum open boundaries values used during      !
!                   4DVar adjustments, tangent linear model.           !
# else
!  tl_ubar_obc   2D U-momentum open boundaries values used during      !
!                   4DVar adjustments, tangent linear model.           !
!  tl_vbar_obc   2D V-momentum open boundaries values used during      !
!                   4DVar adjustments, tangent linear model.           !
# endif
!  tl_zeta_obc   Free-surface open boundaries values used during       !
!                   4DVar adjustments, tangent linear model.           !
#endif
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none
!
!-----------------------------------------------------------------------
!  Lateral boundary condition apply switches.
!-----------------------------------------------------------------------
!
!  The following switches are used to control which grid points are
!  processed by the lateral boundary conditions. These switches are
!  set to TRUE by default.  However in composite grids, the points
!  processed by nesting are set to FALSE to allow mixed boundary
!  conditions along the grid edges.
!
        TYPE T_APPLY
          logical, pointer :: west(:)
          logical, pointer :: east(:)
          logical, pointer :: south(:)
          logical, pointer :: north(:)
        END TYPE

        TYPE (T_APPLY), allocatable :: LBC_apply(:)
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions structure.
!-----------------------------------------------------------------------
!
        TYPE T_BOUNDARY
!
!  Nonlinear model state.
!
          real(r8), pointer :: zeta_west(:)
#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_west_C2(:)
          real(r8), pointer :: zeta_west_Ce(:)
          real(r8), pointer :: zeta_west_Cx(:)
#endif
#ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_west(:,:)
#endif

          real(r8), pointer :: zeta_east(:)
#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_east_C2(:)
          real(r8), pointer :: zeta_east_Ce(:)
          real(r8), pointer :: zeta_east_Cx(:)
#endif
#ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_east(:,:)
#endif

          real(r8), pointer :: zeta_south(:)
#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_south_C2(:)
          real(r8), pointer :: zeta_south_Ce(:)
          real(r8), pointer :: zeta_south_Cx(:)
#endif
#ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_south(:,:)
#endif

          real(r8), pointer :: zeta_north(:)
#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_north_C2(:)
          real(r8), pointer :: zeta_north_Ce(:)
          real(r8), pointer :: zeta_north_Cx(:)
#endif
#ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_north(:,:)
#endif

          real(r8), pointer :: ubar_west(:)
          real(r8), pointer :: vbar_west(:)
#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_west_C2(:)
          real(r8), pointer :: vbar_west_C2(:)
          real(r8), pointer :: ubar_west_Ce(:)
          real(r8), pointer :: vbar_west_Ce(:)
          real(r8), pointer :: ubar_west_Cx(:)
          real(r8), pointer :: vbar_west_Cx(:)
#endif
#ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_west(:,:)
          real(r8), pointer :: vbarG_west(:,:)
#endif

          real(r8), pointer :: ubar_east(:)
          real(r8), pointer :: vbar_east(:)
#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_east_C2(:)
          real(r8), pointer :: vbar_east_C2(:)
          real(r8), pointer :: ubar_east_Ce(:)
          real(r8), pointer :: vbar_east_Ce(:)
          real(r8), pointer :: ubar_east_Cx(:)
          real(r8), pointer :: vbar_east_Cx(:)
#endif
#ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_east(:,:)
          real(r8), pointer :: vbarG_east(:,:)
#endif

          real(r8), pointer :: ubar_south(:)
          real(r8), pointer :: vbar_south(:)
#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_south_C2(:)
          real(r8), pointer :: vbar_south_C2(:)
          real(r8), pointer :: ubar_south_Ce(:)
          real(r8), pointer :: vbar_south_Ce(:)
          real(r8), pointer :: ubar_south_Cx(:)
          real(r8), pointer :: vbar_south_Cx(:)
#endif
#ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_south(:,:)
          real(r8), pointer :: vbarG_south(:,:)
#endif

          real(r8), pointer :: ubar_north(:)
          real(r8), pointer :: vbar_north(:)
#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_north_C2(:)
          real(r8), pointer :: vbar_north_C2(:)
          real(r8), pointer :: ubar_north_Ce(:)
          real(r8), pointer :: vbar_north_Ce(:)
          real(r8), pointer :: ubar_north_Cx(:)
          real(r8), pointer :: vbar_north_Cx(:)
#endif
#ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_north(:,:)
          real(r8), pointer :: vbarG_north(:,:)
#endif
#ifdef WEC
          real(r8), pointer :: ubarstokes_east(:)
          real(r8), pointer :: vbarstokes_east(:)
          real(r8), pointer :: ubarstokes_west(:)
          real(r8), pointer :: vbarstokes_west(:)
          real(r8), pointer :: ubarstokes_north(:)
          real(r8), pointer :: vbarstokes_north(:)
          real(r8), pointer :: ubarstokes_south(:)
          real(r8), pointer :: vbarstokes_south(:)
#endif

#ifdef SOLVE3D
          real(r8), pointer :: u_west(:,:)
          real(r8), pointer :: v_west(:,:)
# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_west_C2(:,:)
          real(r8), pointer :: v_west_C2(:,:)
          real(r8), pointer :: u_west_Ce(:,:)
          real(r8), pointer :: v_west_Ce(:,:)
          real(r8), pointer :: u_west_Cx(:,:)
          real(r8), pointer :: v_west_Cx(:,:)
# endif
# ifndef ANA_M3OBC
          real(r8), pointer :: uG_west(:,:,:)
          real(r8), pointer :: vG_west(:,:,:)
# endif

          real(r8), pointer :: u_east(:,:)
          real(r8), pointer :: v_east(:,:)
# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_east_C2(:,:)
          real(r8), pointer :: v_east_C2(:,:)
          real(r8), pointer :: u_east_Ce(:,:)
          real(r8), pointer :: v_east_Ce(:,:)
          real(r8), pointer :: u_east_Cx(:,:)
          real(r8), pointer :: v_east_Cx(:,:)
# endif
# ifndef ANA_M3OBC
          real(r8), pointer :: uG_east(:,:,:)
          real(r8), pointer :: vG_east(:,:,:)
# endif

          real(r8), pointer :: u_south(:,:)
          real(r8), pointer :: v_south(:,:)
# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_south_C2(:,:)
          real(r8), pointer :: v_south_C2(:,:)
          real(r8), pointer :: u_south_Ce(:,:)
          real(r8), pointer :: v_south_Ce(:,:)
          real(r8), pointer :: u_south_Cx(:,:)
          real(r8), pointer :: v_south_Cx(:,:)
# endif
# ifndef ANA_M3OBC
          real(r8), pointer :: uG_south(:,:,:)
          real(r8), pointer :: vG_south(:,:,:)
# endif

          real(r8), pointer :: u_north(:,:)
          real(r8), pointer :: v_north(:,:)
# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_north_C2(:,:)
          real(r8), pointer :: v_north_C2(:,:)
          real(r8), pointer :: u_north_Ce(:,:)
          real(r8), pointer :: v_north_Ce(:,:)
          real(r8), pointer :: u_north_Cx(:,:)
          real(r8), pointer :: v_north_Cx(:,:)
# endif
# ifndef ANA_M3OBC
          real(r8), pointer :: uG_north(:,:,:)
          real(r8), pointer :: vG_north(:,:,:)
# endif
# ifdef WEC
          real(r8), pointer :: ustokes_east(:,:)
          real(r8), pointer :: vstokes_east(:,:)
          real(r8), pointer :: ustokes_west(:,:)
          real(r8), pointer :: vstokes_west(:,:)
          real(r8), pointer :: ustokes_north(:,:)
          real(r8), pointer :: vstokes_north(:,:)
          real(r8), pointer :: ustokes_south(:,:)
          real(r8), pointer :: vstokes_south(:,:)
# endif

          real(r8), pointer :: t_west(:,:,:)
# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_west_C2(:,:,:)
          real(r8), pointer :: t_west_Ce(:,:,:)
          real(r8), pointer :: t_west_Cx(:,:,:)
# endif
# ifndef ANA_TOBC
          real(r8), pointer :: tG_west(:,:,:,:)
# endif

          real(r8), pointer :: t_east(:,:,:)
# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_east_C2(:,:,:)
          real(r8), pointer :: t_east_Ce(:,:,:)
          real(r8), pointer :: t_east_Cx(:,:,:)
# endif
# ifndef ANA_TOBC
          real(r8), pointer :: tG_east(:,:,:,:)
# endif

          real(r8), pointer :: t_south(:,:,:)
# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_south_C2(:,:,:)
          real(r8), pointer :: t_south_Ce(:,:,:)
          real(r8), pointer :: t_south_Cx(:,:,:)
# endif
# ifndef ANA_TOBC
          real(r8), pointer :: tG_south(:,:,:,:)
# endif

          real(r8), pointer :: t_north(:,:,:)
# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_north_C2(:,:,:)
          real(r8), pointer :: t_north_Ce(:,:,:)
          real(r8), pointer :: t_north_Cx(:,:,:)
# endif
# ifndef ANA_TOBC
          real(r8), pointer :: tG_north(:,:,:,:)
# endif
#endif

#ifdef ICE_MODEL
          real(r8), pointer :: ai_west(:)
          real(r8), pointer :: ai_east(:)
          real(r8), pointer :: ai_south(:)
          real(r8), pointer :: ai_north(:)

# ifndef ANA_AIOBC
          real(r8), pointer :: aiG_west(:,:)
          real(r8), pointer :: aiG_east(:,:)
          real(r8), pointer :: aiG_south(:,:)
          real(r8), pointer :: aiG_north(:,:)
# endif

          real(r8), pointer :: hi_west(:)
          real(r8), pointer :: hi_east(:)
          real(r8), pointer :: hi_south(:)
          real(r8), pointer :: hi_north(:)

# ifndef ANA_HIOBC
          real(r8), pointer :: hiG_west(:,:)
          real(r8), pointer :: hiG_east(:,:)
          real(r8), pointer :: hiG_south(:,:)
          real(r8), pointer :: hiG_north(:,:)
# endif

          real(r8), pointer :: hsn_west(:)
          real(r8), pointer :: hsn_east(:)
          real(r8), pointer :: hsn_south(:)
          real(r8), pointer :: hsn_north(:)

# ifndef ANA_HSNOBC
          real(r8), pointer :: hsnG_west(:,:)
          real(r8), pointer :: hsnG_east(:,:)
          real(r8), pointer :: hsnG_south(:,:)
          real(r8), pointer :: hsnG_north(:,:)
# endif

          real(r8), pointer :: ti_west(:)
          real(r8), pointer :: ti_east(:)
          real(r8), pointer :: ti_south(:)
          real(r8), pointer :: ti_north(:)

# ifndef ANA_TIOBC
          real(r8), pointer :: tiG_west(:,:)
          real(r8), pointer :: tiG_east(:,:)
          real(r8), pointer :: tiG_south(:,:)
          real(r8), pointer :: tiG_north(:,:)
# endif

# ifdef MELT_PONDS
          real(r8), pointer :: apond_west(:)
          real(r8), pointer :: apond_east(:)
          real(r8), pointer :: apond_south(:)
          real(r8), pointer :: apond_north(:)
          real(r8), pointer :: hpond_west(:)
          real(r8), pointer :: hpond_east(:)
          real(r8), pointer :: hpond_south(:)
          real(r8), pointer :: hpond_north(:)

#  ifndef ANA_SFWATOBC
          real(r8), pointer :: apondG_west(:,:)
          real(r8), pointer :: apondG_east(:,:)
          real(r8), pointer :: apondG_south(:,:)
          real(r8), pointer :: apondG_north(:,:)
          real(r8), pointer :: hpondG_west(:,:)
          real(r8), pointer :: hpondG_east(:,:)
          real(r8), pointer :: hpondG_south(:,:)
          real(r8), pointer :: hpondG_north(:,:)
#  endif
# endif

          real(r8), pointer :: sig11_west(:)
          real(r8), pointer :: sig11_east(:)
          real(r8), pointer :: sig11_south(:)
          real(r8), pointer :: sig11_north(:)

# ifndef ANA_SIG11OBC
          real(r8), pointer :: sig11G_west(:,:)
          real(r8), pointer :: sig11G_east(:,:)
          real(r8), pointer :: sig11G_south(:,:)
          real(r8), pointer :: sig11G_north(:,:)
# endif

          real(r8), pointer :: sig22_west(:)
          real(r8), pointer :: sig22_east(:)
          real(r8), pointer :: sig22_south(:)
          real(r8), pointer :: sig22_north(:)

# ifndef ANA_SIG22OBC
          real(r8), pointer :: sig22G_west(:,:)
          real(r8), pointer :: sig22G_east(:,:)
          real(r8), pointer :: sig22G_south(:,:)
          real(r8), pointer :: sig22G_north(:,:)
# endif

          real(r8), pointer :: sig12_west(:)
          real(r8), pointer :: sig12_east(:)
          real(r8), pointer :: sig12_south(:)
          real(r8), pointer :: sig12_north(:)

# ifndef ANA_SIG12OBC
          real(r8), pointer :: sig12G_west(:,:)
          real(r8), pointer :: sig12G_east(:,:)
          real(r8), pointer :: sig12G_south(:,:)
          real(r8), pointer :: sig12G_north(:,:)
# endif

          real(r8), pointer :: ui_west(:)
          real(r8), pointer :: ui_east(:)
          real(r8), pointer :: ui_south(:)
          real(r8), pointer :: ui_north(:)

# ifndef ANA_MIOBC
          real(r8), pointer :: uiG_west(:,:)
          real(r8), pointer :: uiG_east(:,:)
          real(r8), pointer :: uiG_south(:,:)
          real(r8), pointer :: uiG_north(:,:)
# endif

          real(r8), pointer :: vi_west(:)
          real(r8), pointer :: vi_east(:)
          real(r8), pointer :: vi_south(:)
          real(r8), pointer :: vi_north(:)

# ifndef ANA_MIOBC
          real(r8), pointer :: viG_west(:,:)
          real(r8), pointer :: viG_east(:,:)
          real(r8), pointer :: viG_south(:,:)
          real(r8), pointer :: viG_north(:,:)
# endif

#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
          real(r8), pointer :: tl_zeta_west(:)
          real(r8), pointer :: tl_zeta_east(:)
          real(r8), pointer :: tl_zeta_south(:)
          real(r8), pointer :: tl_zeta_north(:)

          real(r8), pointer :: tl_ubar_west(:)
          real(r8), pointer :: tl_vbar_west(:)
          real(r8), pointer :: tl_ubar_east(:)
          real(r8), pointer :: tl_vbar_east(:)
          real(r8), pointer :: tl_ubar_south(:)
          real(r8), pointer :: tl_vbar_south(:)
          real(r8), pointer :: tl_ubar_north(:)
          real(r8), pointer :: tl_vbar_north(:)

# ifdef SOLVE3D
          real(r8), pointer :: tl_u_west(:,:)
          real(r8), pointer :: tl_v_west(:,:)
          real(r8), pointer :: tl_u_east(:,:)
          real(r8), pointer :: tl_v_east(:,:)
          real(r8), pointer :: tl_u_south(:,:)
          real(r8), pointer :: tl_v_south(:,:)
          real(r8), pointer :: tl_u_north(:,:)
          real(r8), pointer :: tl_v_north(:,:)

          real(r8), pointer :: tl_t_west(:,:,:)
          real(r8), pointer :: tl_t_east(:,:,:)
          real(r8), pointer :: tl_t_south(:,:,:)
          real(r8), pointer :: tl_t_north(:,:,:)
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
          real(r8), pointer :: ad_zeta_west(:)
          real(r8), pointer :: ad_zeta_east(:)
          real(r8), pointer :: ad_zeta_south(:)
          real(r8), pointer :: ad_zeta_north(:)

          real(r8), pointer :: ad_ubar_west(:)
          real(r8), pointer :: ad_vbar_west(:)
          real(r8), pointer :: ad_ubar_east(:)
          real(r8), pointer :: ad_vbar_east(:)
          real(r8), pointer :: ad_ubar_south(:)
          real(r8), pointer :: ad_vbar_south(:)
          real(r8), pointer :: ad_ubar_north(:)
          real(r8), pointer :: ad_vbar_north(:)

# ifdef SOLVE3D
          real(r8), pointer :: ad_u_west(:,:)
          real(r8), pointer :: ad_v_west(:,:)
          real(r8), pointer :: ad_u_east(:,:)
          real(r8), pointer :: ad_v_east(:,:)
          real(r8), pointer :: ad_u_south(:,:)
          real(r8), pointer :: ad_v_south(:,:)
          real(r8), pointer :: ad_u_north(:,:)
          real(r8), pointer :: ad_v_north(:,:)

          real(r8), pointer :: ad_t_west(:,:,:)
          real(r8), pointer :: ad_t_east(:,:,:)
          real(r8), pointer :: ad_t_south(:,:,:)
          real(r8), pointer :: ad_t_north(:,:,:)
# endif
#endif

#ifdef ADJUST_BOUNDARY
!
!  Open boundaries arrays used in 4DVar adjustments.
!
# ifdef SOLVE3D
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
# endif
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
#endif

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
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
#ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
#endif
      integer :: my_tile
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  See dimension ranges. Notice that the boundary arrays are dimensioned
!  with the global dimensions of grid. That is, no tiling ranges in
!  distributed-memory. This is done to facilitate processing.
!
      my_tile=-1                           ! for global values
      LBi=BOUNDS(ng)%LBi(my_tile)
      UBi=BOUNDS(ng)%UBi(my_tile)
      LBj=BOUNDS(ng)%LBj(my_tile)
      UBj=BOUNDS(ng)%UBj(my_tile)
#ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
#endif
!
!  Allocate structures.
!
      IF (ng.eq.1) THEN
        allocate ( LBC_apply(Ngrids) )
        allocate ( BOUNDARY(Ngrids) )
      END IF
!
!  Lateral boundary conditions apply switches.  These switches need to
!  be initilized to TRUE here because 'initialize_boundary' is called
!  several times in adjoint-based application to clear state arrays.
!  These switches are part of the application grid and will be set to
!  FALSE elsewhere, if the boundary point is assigned by a nested grid.
!
      allocate ( LBC_apply(ng) % west(LBj:UBj) )
      LBC_apply(ng) % west = .TRUE.

      allocate ( LBC_apply(ng) % east(LBj:UBj) )
      LBC_apply(ng) % east = .TRUE.

      allocate ( LBC_apply(ng) % south(LBi:UBi) )
      LBC_apply(ng) % south = .TRUE.

      allocate ( LBC_apply(ng) % north(LBi:UBi) )
      LBC_apply(ng) % north = .TRUE.
!
!  Nonlinear model state.
!
      IF (LBC(iwest,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_west(LBj:UBj) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % zeta_west_C2(LBj:UBj) )
        allocate ( BOUNDARY(ng) % zeta_west_Ce(LBj:UBj) )
        allocate ( BOUNDARY(ng) % zeta_west_Cx(LBj:UBj) )
#endif
#ifndef ANA_FSOBC
        allocate ( BOUNDARY(ng) % zetaG_west(LBj:UBj,2) )
#endif
      END IF

      IF (LBC(ieast,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_east(LBj:UBj) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % zeta_east_C2(LBj:UBj) )
        allocate ( BOUNDARY(ng) % zeta_east_Ce(LBj:UBj) )
        allocate ( BOUNDARY(ng) % zeta_east_Cx(LBj:UBj) )
#endif
#ifndef ANA_FSOBC
        allocate ( BOUNDARY(ng) % zetaG_east(LBj:UBj,2) )
#endif
      END IF

      IF (LBC(isouth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_south(LBi:UBi) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % zeta_south_C2(LBi:UBi) )
        allocate ( BOUNDARY(ng) % zeta_south_Ce(LBi:UBi) )
        allocate ( BOUNDARY(ng) % zeta_south_Cx(LBi:UBi) )
#endif
#ifndef ANA_FSOBC
        allocate ( BOUNDARY(ng) % zetaG_south(LBi:UBi,2) )
#endif
      END IF

      IF (LBC(inorth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_north(LBi:UBi) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % zeta_north_C2(LBi:UBi) )
        allocate ( BOUNDARY(ng) % zeta_north_Ce(LBi:UBi) )
        allocate ( BOUNDARY(ng) % zeta_north_Cx(LBi:UBi) )
#endif
#ifndef ANA_FSOBC
        allocate ( BOUNDARY(ng) % zetaG_north(LBi:UBi,2) )
#endif
      END IF
!
      IF (LBC(iwest,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_west(LBj:UBj) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % ubar_west_C2(LBj:UBj) )
        allocate ( BOUNDARY(ng) % ubar_west_Ce(LBj:UBj) )
        allocate ( BOUNDARY(ng) % ubar_west_Cx(LBj:UBj) )
#endif
#ifndef ANA_M2OBC
        allocate ( BOUNDARY(ng) % ubarG_west(LBj:UBj,2) )
#endif
      END IF

      IF (LBC(ieast,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_east(LBj:UBj) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % ubar_east_C2(LBj:UBj) )
        allocate ( BOUNDARY(ng) % ubar_east_Ce(LBj:UBj) )
        allocate ( BOUNDARY(ng) % ubar_east_Cx(LBj:UBj) )
#endif
#ifndef ANA_M2OBC
        allocate ( BOUNDARY(ng) % ubarG_east(LBj:UBj,2) )
#endif
      END IF

      IF (LBC(isouth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_south(LBi:UBi) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % ubar_south_C2(LBi:UBi) )
        allocate ( BOUNDARY(ng) % ubar_south_Ce(LBi:UBi) )
        allocate ( BOUNDARY(ng) % ubar_south_Cx(LBi:UBi) )
#endif
#ifndef ANA_M2OBC
        allocate ( BOUNDARY(ng) % ubarG_south(LBi:UBi,2) )
#endif
      END IF

      IF (LBC(inorth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_north(LBi:UBi) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % ubar_north_C2(LBi:UBi) )
        allocate ( BOUNDARY(ng) % ubar_north_Ce(LBi:UBi) )
        allocate ( BOUNDARY(ng) % ubar_north_Cx(LBi:UBi) )
#endif
#ifndef ANA_M2OBC
        allocate ( BOUNDARY(ng) % ubarG_north(LBi:UBi,2) )
#endif
      END IF
!
      IF (LBC(iwest,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_west(LBj:UBj) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % vbar_west_C2(LBj:UBj) )
        allocate ( BOUNDARY(ng) % vbar_west_Ce(LBj:UBj) )
        allocate ( BOUNDARY(ng) % vbar_west_Cx(LBj:UBj) )
#endif
#ifndef ANA_M2OBC
        allocate ( BOUNDARY(ng) % vbarG_west(LBj:UBj,2) )
#endif
      END IF

      IF (LBC(ieast,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_east(LBj:UBj) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % vbar_east_C2(LBj:UBj) )
        allocate ( BOUNDARY(ng) % vbar_east_Ce(LBj:UBj) )
        allocate ( BOUNDARY(ng) % vbar_east_Cx(LBj:UBj) )
#endif
#ifndef ANA_M2OBC
        allocate ( BOUNDARY(ng) % vbarG_east(LBj:UBj,2) )
#endif
      END IF

      IF (LBC(isouth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_south(LBi:UBi) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % vbar_south_C2(LBi:UBi) )
        allocate ( BOUNDARY(ng) % vbar_south_Ce(LBi:UBi) )
        allocate ( BOUNDARY(ng) % vbar_south_Cx(LBi:UBi) )
#endif
#ifndef ANA_M2OBC
        allocate ( BOUNDARY(ng) % vbarG_south(LBi:UBi,2) )
#endif
      END IF

      IF (LBC(inorth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_north(LBi:UBi) )
#if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % vbar_north_C2(LBi:UBi) )
        allocate ( BOUNDARY(ng) % vbar_north_Ce(LBi:UBi) )
        allocate ( BOUNDARY(ng) % vbar_north_Cx(LBi:UBi) )
#endif
#ifndef ANA_M2OBC
        allocate ( BOUNDARY(ng) % vbarG_north(LBi:UBi,2) )
#endif
      END IF

#ifdef SOLVE3D
!
      IF (LBC(iwest,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % u_west(LBj:UBj,N(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % u_west_C2(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % u_west_Ce(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % u_west_Cx(LBj:UBj,N(ng)) )
# endif
# ifndef ANA_M3OBC
        allocate ( BOUNDARY(ng) % uG_west(LBj:UBj,N(ng),2) )
# endif
      END IF

      IF (LBC(ieast,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % u_east(LBj:UBj,N(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % u_east_C2(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % u_east_Ce(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % u_east_Cx(LBj:UBj,N(ng)) )
# endif
# ifndef ANA_M3OBC
        allocate ( BOUNDARY(ng) % uG_east(LBj:UBj,N(ng),2) )
# endif
      END IF

      IF (LBC(isouth,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % u_south(LBi:UBi,N(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % u_south_C2(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % u_south_Ce(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % u_south_Cx(LBi:UBi,N(ng)) )
# endif
# ifndef ANA_M3OBC
        allocate ( BOUNDARY(ng) % uG_south(LBi:UBi,N(ng),2) )
# endif
      END IF

      IF (LBC(inorth,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % u_north(LBi:UBi,N(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % u_north_C2(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % u_north_Ce(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % u_north_Cx(LBi:UBi,N(ng)) )
# endif
# ifndef ANA_M3OBC
        allocate ( BOUNDARY(ng) % uG_north(LBi:UBi,N(ng),2) )
# endif
      END IF
!
      IF (LBC(iwest,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % v_west(LBj:UBj,N(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % v_west_C2(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % v_west_Ce(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % v_west_Cx(LBj:UBj,N(ng)) )
# endif
# ifndef ANA_M3OBC
        allocate ( BOUNDARY(ng) % vG_west(LBj:UBj,N(ng),2) )
# endif
      END IF

      IF (LBC(ieast,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % v_east(LBj:UBj,N(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % v_east_C2(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % v_east_Ce(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % v_east_Cx(LBj:UBj,N(ng)) )
# endif
# ifndef ANA_M3OBC
        allocate ( BOUNDARY(ng) % vG_east(LBj:UBj,N(ng),2) )
# endif
      END IF

      IF (LBC(isouth,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % v_south(LBi:UBi,N(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % v_south_C2(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % v_south_Ce(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % v_south_Cx(LBi:UBi,N(ng)) )
# endif
# ifndef ANA_M3OBC
        allocate ( BOUNDARY(ng) % vG_south(LBi:UBi,N(ng),2) )
# endif
      END IF

      IF (LBC(inorth,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % v_north(LBi:UBi,N(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % u_north_C2(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % u_north_Ce(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % u_north_Cx(LBi:UBi,N(ng)) )
# endif
# ifndef ANA_M3OBC
        allocate ( BOUNDARY(ng) % vG_north(LBi:UBi,N(ng),2) )
# endif
      END IF
!

      IF (ANY(LBC(iwest,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % t_west(LBj:UBj,N(ng),NT(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % t_west_C2(LBj:UBj,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % t_west_Ce(LBj:UBj,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % t_west_Cx(LBj:UBj,N(ng),NT(ng)) )
# endif
# ifndef ANA_TOBC
        allocate ( BOUNDARY(ng) % tG_west(LBj:UBj,N(ng),2,NT(ng)) )
# endif
      END IF

      IF (ANY(LBC(ieast,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % t_east(LBj:UBj,N(ng),NT(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % t_east_C2(LBj:UBj,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % t_east_Ce(LBj:UBj,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % t_east_Cx(LBj:UBj,N(ng),NT(ng)) )
# endif
# ifndef ANA_TOBC
        allocate ( BOUNDARY(ng) % tG_east(LBj:UBj,N(ng),2,NT(ng)) )
# endif
      END IF

      IF (ANY(LBC(isouth,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % t_south(LBi:UBi,N(ng),NT(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % t_south_C2(LBi:UBi,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % t_south_Ce(LBi:UBi,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % t_south_Cx(LBi:UBi,N(ng),NT(ng)) )
# endif
# ifndef ANA_TOBC
        allocate ( BOUNDARY(ng) % tG_south(LBi:UBi,N(ng),2,NT(ng)) )
# endif
      END IF

      IF (ANY(LBC(inorth,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % t_north(LBi:UBi,N(ng),NT(ng)) )
# if defined CELERITY_READ || defined CELERITY_WRITE
        allocate ( BOUNDARY(ng) % t_north_C2(LBi:UBi,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % t_north_Ce(LBi:UBi,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % t_north_Cx(LBi:UBi,N(ng),NT(ng)) )
# endif
# ifndef ANA_TOBC
        allocate ( BOUNDARY(ng) % tG_north(LBi:UBi,N(ng),2,NT(ng)) )
# endif
      END IF
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF (tl_LBC(iwest,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_zeta_west(LBj:UBj) )
      END IF

      IF (tl_LBC(ieast,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_zeta_east(LBj:UBj) )
      END IF

      IF (tl_LBC(isouth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_zeta_south(LBi:UBi) )
      END IF

      IF (tl_LBC(inorth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_zeta_north(LBi:UBi) )
      END IF
!
      IF (tl_LBC(iwest,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_ubar_west(LBj:UBj) )
      END IF

      IF (tl_LBC(ieast,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_ubar_east(LBj:UBj) )
      END IF

      IF (tl_LBC(isouth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_ubar_south(LBi:UBi) )
      END IF

      IF (tl_LBC(inorth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_ubar_north(LBi:UBi) )
      END IF
!
      IF (tl_LBC(iwest,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_vbar_west(LBj:UBj) )
      END IF

      IF (tl_LBC(ieast,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_vbar_east(LBj:UBj) )
      END IF

      IF (tl_LBC(isouth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_vbar_south(LBi:UBi) )
      END IF

      IF (tl_LBC(inorth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_vbar_north(LBi:UBi) )
      END IF

# ifdef SOLVE3D
!
      IF (tl_LBC(iwest,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_u_west(LBj:UBj,N(ng)) )
      END IF

      IF (tl_LBC(ieast,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_u_east(LBj:UBj,N(ng)) )
      END IF

      IF (tl_LBC(isouth,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_u_south(LBi:UBi,N(ng)) )
      END IF

      IF (tl_LBC(inorth,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_u_north(LBi:UBi,N(ng)) )
      END IF
!
      IF (tl_LBC(iwest,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_v_west(LBj:UBj,N(ng)) )
      END IF

      IF (tl_LBC(ieast,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_v_east(LBj:UBj,N(ng)) )
      END IF

      IF (tl_LBC(isouth,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_v_south(LBi:UBi,N(ng)) )
      END IF

      IF (tl_LBC(inorth,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tl_v_north(LBi:UBi,N(ng)) )
      END IF
!
      IF (ANY(tl_LBC(iwest,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % tl_t_west(LBj:UBj,N(ng),NT(ng)) )
      END IF

      IF (ANY(tl_LBC(ieast,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % tl_t_east(LBj:UBj,N(ng),NT(ng)) )
      END IF

      IF (ANY(tl_LBC(isouth,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % tl_t_south(LBi:UBi,N(ng),NT(ng)) )
      END IF

      IF (ANY(tl_LBC(inorth,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % tl_t_north(LBi:UBi,N(ng),NT(ng)) )
      END IF
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      IF (ad_LBC(iwest,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_zeta_west(LBj:UBj) )
      END IF

      IF (ad_LBC(ieast,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_zeta_east(LBj:UBj) )
      END IF

      IF (ad_LBC(isouth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_zeta_south(LBi:UBi) )
      END IF

      IF (ad_LBC(inorth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_zeta_north(LBi:UBi) )
      END IF
!
      IF (ad_LBC(iwest,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_ubar_west(LBj:UBj) )
      END IF

      IF (ad_LBC(ieast,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_ubar_east(LBj:UBj) )
      END IF

      IF (ad_LBC(isouth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_ubar_south(LBi:UBi) )
      END IF

      IF (ad_LBC(inorth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_ubar_north(LBi:UBi) )
      END IF
!
      IF (ad_LBC(iwest,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_vbar_west(LBj:UBj) )
      END IF

      IF (ad_LBC(ieast,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_vbar_east(LBj:UBj) )
      END IF

      IF (ad_LBC(isouth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_vbar_south(LBi:UBi) )
      END IF

      IF (ad_LBC(inorth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_vbar_north(LBi:UBi) )
      END IF

# ifdef SOLVE3D
!
      IF (ad_LBC(iwest,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_u_west(LBj:UBj,N(ng)) )
      END IF

      IF (ad_LBC(ieast,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_u_east(LBj:UBj,N(ng)) )
      END IF

      IF (ad_LBC(isouth,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_u_south(LBi:UBi,N(ng)) )
      END IF

      IF (ad_LBC(inorth,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_u_north(LBi:UBi,N(ng)) )
      END IF
!
      IF (ad_LBC(iwest,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_v_west(LBj:UBj,N(ng)) )
      END IF

      IF (ad_LBC(ieast,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_v_east(LBj:UBj,N(ng)) )
      END IF

      IF (ad_LBC(isouth,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_v_south(LBi:UBi,N(ng)) )
      END IF

      IF (ad_LBC(inorth,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ad_v_north(LBi:UBi,N(ng)) )
      END IF
!
      IF (ANY(ad_LBC(iwest,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % ad_t_west(LBj:UBj,N(ng),NT(ng)) )
      END IF

      IF (ANY(ad_LBC(ieast,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % ad_t_east(LBj:UBj,N(ng),NT(ng)) )
      END IF

      IF (ANY(ad_LBC(isouth,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % ad_t_south(LBi:UBi,N(ng),NT(ng)) )
      END IF

      IF (ANY(ad_LBC(inorth,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % ad_t_north(LBi:UBi,N(ng),NT(ng)) )
      END IF
# endif
#endif

#ifdef ADJUST_BOUNDARY
!
!  Open boundaries arrays used in 4DVar adjustments.
!
# ifdef SOLVE3D
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
# endif
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
#endif

#ifdef ICE_MODEL
      allocate ( BOUNDARY(ng) % ai_west(LBj:UBj) )
      allocate ( BOUNDARY(ng) % ai_east(LBj:UBj) )
      allocate ( BOUNDARY(ng) % ai_south(LBi:UBi) )
      allocate ( BOUNDARY(ng) % ai_north(LBi:UBi) )

# ifndef ANA_AIOBC
      IF (LBC(iwest,isAice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % aiG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isAice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % aiG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isAice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % aiG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isAice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % aiG_north(LBi:UBi,2) )
      END IF
# endif

      allocate ( BOUNDARY(ng) % hi_west(LBj:UBj) )
      allocate ( BOUNDARY(ng) % hi_east(LBj:UBj) )
      allocate ( BOUNDARY(ng) % hi_south(LBi:UBi) )
      allocate ( BOUNDARY(ng) % hi_north(LBi:UBi) )

# ifndef ANA_HIOBC
      IF (LBC(iwest,isHice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hiG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isHice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hiG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isHice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hiG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isHice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hiG_north(LBi:UBi,2) )
      END IF
# endif

      allocate ( BOUNDARY(ng) % hsn_west(LBj:UBj) )
      allocate ( BOUNDARY(ng) % hsn_east(LBj:UBj) )
      allocate ( BOUNDARY(ng) % hsn_south(LBi:UBi) )
      allocate ( BOUNDARY(ng) % hsn_north(LBi:UBi) )

# ifndef ANA_HSNOBC
      IF (LBC(iwest,isHsno,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hsnG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isHsno,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hsnG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isHsno,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hsnG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isHsno,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hsnG_north(LBi:UBi,2) )
      END IF
# endif

      IF (LBC(iwest,isTice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ti_west(LBj:UBj) )
      END IF
      IF (LBC(ieast,isTice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ti_east(LBj:UBj) )
      END IF
      IF (LBC(isouth,isTice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ti_south(LBi:UBi) )
      END IF
      IF (LBC(inorth,isTice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ti_north(LBi:UBi) )
      END IF

# ifndef ANA_TIOBC
      IF (LBC(iwest,isTice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tiG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isTice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tiG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isTice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tiG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isTice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % tiG_north(LBi:UBi,2) )
      END IF
# endif

# ifdef MELT_PONDS
      allocate ( BOUNDARY(ng) % apond_west(LBj:UBj) )
      allocate ( BOUNDARY(ng) % apond_east(LBj:UBj) )
      allocate ( BOUNDARY(ng) % apond_south(LBi:UBi) )
      allocate ( BOUNDARY(ng) % apond_north(LBi:UBi) )
      allocate ( BOUNDARY(ng) % hpond_west(LBj:UBj) )
      allocate ( BOUNDARY(ng) % hpond_east(LBj:UBj) )
      allocate ( BOUNDARY(ng) % hpond_south(LBi:UBi) )
      allocate ( BOUNDARY(ng) % hpond_north(LBi:UBi) )

#  ifndef ANA_SFWATOBC
      IF (LBC(iwest,isApond,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % apondG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isApond,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % apondG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isApond,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % apondG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isApond,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % apondG_north(LBi:UBi,2) )
      END IF
      IF (LBC(iwest,isHpond,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hpondG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isHpond,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hpondG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isHpond,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hpondG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isHpond,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % hpondG_north(LBi:UBi,2) )
      END IF
#  endif
# endif

      allocate ( BOUNDARY(ng) % sig11_west(LBj:UBj) )
      allocate ( BOUNDARY(ng) % sig11_east(LBj:UBj) )
      allocate ( BOUNDARY(ng) % sig11_south(LBi:UBi) )
      allocate ( BOUNDARY(ng) % sig11_north(LBi:UBi) )

# ifndef ANA_SIG11OBC
      IF (LBC(iwest,isSig11,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig11G_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isSig11,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig11G_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isSig11,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig11G_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isSig11,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig11G_north(LBi:UBi,2) )
      END IF
# endif

      allocate ( BOUNDARY(ng) % sig22_west(LBj:UBj) )
      allocate ( BOUNDARY(ng) % sig22_east(LBj:UBj) )
      allocate ( BOUNDARY(ng) % sig22_south(LBi:UBi) )
      allocate ( BOUNDARY(ng) % sig22_north(LBi:UBi) )

# ifndef ANA_SIG22OBC
      IF (LBC(iwest,isSig22,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig22G_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isSig22,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig22G_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isSig22,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig22G_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isSig22,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig22G_north(LBi:UBi,2) )
      END IF
# endif

      allocate ( BOUNDARY(ng) % sig12_west(LBj:UBj) )
      allocate ( BOUNDARY(ng) % sig12_east(LBj:UBj) )
      allocate ( BOUNDARY(ng) % sig12_south(LBi:UBi) )
      allocate ( BOUNDARY(ng) % sig12_north(LBi:UBi) )

# ifndef ANA_SIG12OBC
      IF (LBC(iwest,isSig12,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig12G_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isSig12,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig12G_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isSig12,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig12G_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isSig12,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % sig12G_north(LBi:UBi,2) )
      END IF
# endif

      IF (LBC(iwest,isUice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ui_west(LBj:UBj) )
      END IF
      IF (LBC(ieast,isUice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ui_east(LBj:UBj) )
      END IF
      IF (LBC(isouth,isUice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ui_south(LBi:UBi) )
      END IF
      IF (LBC(inorth,isUice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ui_north(LBi:UBi) )
      END IF

# ifndef ANA_MIOBC
      IF (LBC(iwest,isUice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % uiG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isUice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % uiG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isUice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % uiG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isUice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % uiG_north(LBi:UBi,2) )
      END IF
# endif

      IF (LBC(iwest,isVice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vi_west(LBj:UBj) )
      END IF
      IF (LBC(ieast,isVice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vi_east(LBj:UBj) )
      END IF
      IF (LBC(isouth,isVice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vi_south(LBi:UBi) )
      END IF
      IF (LBC(inorth,isVice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vi_north(LBi:UBi) )
      END IF

# ifndef ANA_MIOBC
      IF (LBC(iwest,isVice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % viG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isVice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % viG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isVice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % viG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isVice,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % viG_north(LBi:UBi,2) )
      END IF
# endif

#endif

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
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      real(r8), parameter :: IniVal = 0.0_r8

#ifndef DISTRIBUTE
# include "set_bounds.h"
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN

        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_west = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_west_C2 = IniVal
          BOUNDARY(ng) % zeta_west_Ce = IniVal
          BOUNDARY(ng) % zeta_west_Cx = IniVal
#endif
#ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_west = IniVal
#endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_east = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_east_C2 = IniVal
          BOUNDARY(ng) % zeta_east_Ce = IniVal
          BOUNDARY(ng) % zeta_east_Cx = IniVal
#endif
#ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_east = IniVal
#endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_south = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_south_C2 = IniVal
          BOUNDARY(ng) % zeta_south_Ce = IniVal
          BOUNDARY(ng) % zeta_south_Cx = IniVal
#endif
#ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_south = IniVal
#endif
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_north = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_north_C2 = IniVal
          BOUNDARY(ng) % zeta_north_Ce = IniVal
          BOUNDARY(ng) % zeta_north_Cx = IniVal
#endif
#ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_north = IniVal
#endif
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_west = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % ubar_west_C2 = IniVal
          BOUNDARY(ng) % ubar_west_Ce = IniVal
          BOUNDARY(ng) % ubar_west_Cx = IniVal
#endif
#ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_west = IniVal
#endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_east = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % ubar_east_C2 = IniVal
          BOUNDARY(ng) % ubar_east_Ce = IniVal
          BOUNDARY(ng) % ubar_east_Cx = IniVal
#endif
#ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_east = IniVal
#endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_south = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % ubar_south_C2 = IniVal
          BOUNDARY(ng) % ubar_south_Ce = IniVal
          BOUNDARY(ng) % ubar_south_Cx = IniVal
#endif
#ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_south = IniVal
#endif
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_north = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % ubar_north_C2 = IniVal
          BOUNDARY(ng) % ubar_north_Ce = IniVal
          BOUNDARY(ng) % ubar_north_Cx = IniVal
#endif
#ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_north = IniVal
          BOUNDARY(ng) % vbarG_north = IniVal
#endif
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_west = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % vbar_west_C2 = IniVal
          BOUNDARY(ng) % vbar_west_Ce = IniVal
          BOUNDARY(ng) % vbar_west_Cx = IniVal
#endif
#ifndef ANA_M2OBC
          BOUNDARY(ng) % vbarG_west = IniVal
#endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_east = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % vbar_east_C2 = IniVal
          BOUNDARY(ng) % vbar_east_Ce = IniVal
          BOUNDARY(ng) % vbar_east_Cx = IniVal
#endif
#ifndef ANA_M2OBC
          BOUNDARY(ng) % vbarG_east = IniVal
#endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_south = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % vbar_south_C2 = IniVal
          BOUNDARY(ng) % vbar_south_Ce = IniVal
          BOUNDARY(ng) % vbar_south_Cx = IniVal
#endif
#ifndef ANA_M2OBC
          BOUNDARY(ng) % vbarG_south = IniVal
#endif
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_north = IniVal
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % vbar_north_C2 = IniVal
          BOUNDARY(ng) % vbar_north_Ce = IniVal
          BOUNDARY(ng) % vbar_north_Cx = IniVal
#endif
#ifndef ANA_M2OBC
          BOUNDARY(ng) % vbarG_north = IniVal
#endif
        END IF

#ifdef SOLVE3D
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % u_west = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % u_west_C2 = IniVal
          BOUNDARY(ng) % u_west_Ce = IniVal
          BOUNDARY(ng) % u_west_Cx = IniVal
# endif
# ifndef ANA_M3OBC
          BOUNDARY(ng) % uG_west = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % u_east = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % u_east_C2 = IniVal
          BOUNDARY(ng) % u_east_Ce = IniVal
          BOUNDARY(ng) % u_east_Cx = IniVal
# endif
# ifndef ANA_M3OBC
          BOUNDARY(ng) % uG_east = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % u_south = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % u_south_C2 = IniVal
          BOUNDARY(ng) % u_south_Ce = IniVal
          BOUNDARY(ng) % u_south_Cx = IniVal
# endif
# ifndef ANA_M3OBC
          BOUNDARY(ng) % uG_south = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % u_north = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % u_north_C2 = IniVal
          BOUNDARY(ng) % u_north_Ce = IniVal
          BOUNDARY(ng) % u_north_Cx = IniVal
# endif
# ifndef ANA_M3OBC
          BOUNDARY(ng) % uG_north = IniVal
# endif
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % v_west = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % v_west_C2 = IniVal
          BOUNDARY(ng) % v_west_Ce = IniVal
          BOUNDARY(ng) % v_west_Cx = IniVal
# endif
# ifndef ANA_M3OBC
          BOUNDARY(ng) % vG_west = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % v_east = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % v_east_C2 = IniVal
          BOUNDARY(ng) % v_east_Ce = IniVal
          BOUNDARY(ng) % v_east_Cx = IniVal
# endif
# ifndef ANA_M3OBC
          BOUNDARY(ng) % vG_east = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % v_south = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % v_south_C2 = IniVal
          BOUNDARY(ng) % v_south_Ce = IniVal
          BOUNDARY(ng) % v_south_Cx = IniVal
# endif
# ifndef ANA_M3OBC
          BOUNDARY(ng) % vG_south = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % v_north = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % v_north_C2 = IniVal
          BOUNDARY(ng) % v_north_Ce = IniVal
          BOUNDARY(ng) % v_north_Cx = IniVal
# endif
# ifndef ANA_M3OBC
          BOUNDARY(ng) % vG_north = IniVal
# endif
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      ANY(LBC(iwest,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % t_west = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % t_west_C2 = IniVal
          BOUNDARY(ng) % t_west_Ce = IniVal
          BOUNDARY(ng) % t_west_Cx = IniVal
# endif
# ifndef ANA_TOBC
          BOUNDARY(ng) % tG_west = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      ANY(LBC(ieast,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % t_east = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % t_east_C2 = IniVal
          BOUNDARY(ng) % t_east_Ce = IniVal
          BOUNDARY(ng) % t_east_Cx = IniVal
# endif
# ifndef ANA_TOBC
          BOUNDARY(ng) % tG_east = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      ANY(LBC(isouth,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % t_south = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % t_south_C2 = IniVal
          BOUNDARY(ng) % t_south_Ce = IniVal
          BOUNDARY(ng) % t_south_Cx = IniVal
# endif
# ifndef ANA_TOBC
          BOUNDARY(ng) % tG_south = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      ANY(LBC(inorth,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % t_north = IniVal
# if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % t_north_C2 = IniVal
          BOUNDARY(ng) % t_north_Ce = IniVal
          BOUNDARY(ng) % t_north_Cx = IniVal
# endif
# ifndef ANA_TOBC
          BOUNDARY(ng) % tG_north = IniVal
# endif
        END IF
#endif

# ifdef ICE_MODEL
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isAice,ng)%acquire) THEN
          BOUNDARY(ng) % ai_west = IniVal
        END IF
#  ifndef ANA_AIOBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isAice,ng)%acquire) THEN
          BOUNDARY(ng) % aiG_west = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isHice,ng)%acquire) THEN
          BOUNDARY(ng) % hi_west = IniVal
        END IF
#  ifndef ANA_HIOBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isHice,ng)%acquire) THEN
          BOUNDARY(ng) % hiG_west = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isHsno,ng)%acquire) THEN
          BOUNDARY(ng) % hsn_west = IniVal
        END IF
#  ifndef ANA_HSNOBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isHsno,ng)%acquire) THEN
          BOUNDARY(ng) % hsnG_west = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isTice,ng)%acquire) THEN
          BOUNDARY(ng) % ti_west = IniVal
        END IF
#  ifndef ANA_TIOBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isTice,ng)%acquire) THEN
          BOUNDARY(ng) % tiG_west = IniVal
        END IF
#  endif
#  ifdef MELT_PONDS
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isApond,ng)%acquire) THEN
          BOUNDARY(ng) % apond_west = IniVal
        END IF
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isHpond,ng)%acquire) THEN
          BOUNDARY(ng) % hpond_west = IniVal
        END IF
#   ifndef ANA_SFWATOBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isApond,ng)%acquire) THEN
          BOUNDARY(ng) % apondG_west = IniVal
        END IF
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isHpond,ng)%acquire) THEN
          BOUNDARY(ng) % hpondG_west = IniVal
        END IF
#   endif
#  endif
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isSig11,ng)%acquire) THEN
          BOUNDARY(ng) % sig11_west = IniVal
        END IF
#  ifndef ANA_SIG11OBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isSig11,ng)%acquire) THEN
          BOUNDARY(ng) % sig11G_west = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isSig22,ng)%acquire) THEN
          BOUNDARY(ng) % sig22_west = IniVal
        END IF
#  ifndef ANA_SIG22OBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isSig22,ng)%acquire) THEN
          BOUNDARY(ng) % sig22G_west = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isSig12,ng)%acquire) THEN
          BOUNDARY(ng) % sig12_west = IniVal
        END IF
#  ifndef ANA_SIG12OBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isSig12,ng)%acquire) THEN
          BOUNDARY(ng) % sig12G_west = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isUice,ng)%acquire) THEN
          BOUNDARY(ng) % ui_west = IniVal
        END IF
#  ifndef ANA_MIOBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isUice,ng)%acquire) THEN
          BOUNDARY(ng) % uiG_west = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isVice,ng)%acquire) THEN
          BOUNDARY(ng) % vi_west = IniVal
        END IF
#  ifndef ANA_MIOBC
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isVice,ng)%acquire) THEN
          BOUNDARY(ng) % viG_west = IniVal
        END IF
#  endif

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isAice,ng)%acquire) THEN
          BOUNDARY(ng) % ai_east = IniVal
        END IF
#  ifndef ANA_AIOBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isAice,ng)%acquire) THEN
          BOUNDARY(ng) % aiG_east = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isHice,ng)%acquire) THEN
          BOUNDARY(ng) % hi_east = IniVal
        END IF
#  ifndef ANA_HIOBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isHice,ng)%acquire) THEN
          BOUNDARY(ng) % hiG_east = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isHsno,ng)%acquire) THEN
          BOUNDARY(ng) % hsn_east = IniVal
        END IF
#  ifndef ANA_HSNOBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isHsno,ng)%acquire) THEN
          BOUNDARY(ng) % hsnG_east = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isTice,ng)%acquire) THEN
          BOUNDARY(ng) % ti_east = IniVal
        END IF
#  ifndef ANA_TIOBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isTice,ng)%acquire) THEN
          BOUNDARY(ng) % tiG_east = IniVal
        END IF
#  endif
#  ifdef MELT_PONDS
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isApond,ng)%acquire) THEN
          BOUNDARY(ng) % apond_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isHpond,ng)%acquire) THEN
          BOUNDARY(ng) % hpond_east = IniVal
        END IF
#   ifndef ANA_SFWATOBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isApond,ng)%acquire) THEN
          BOUNDARY(ng) % apondG_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isHpond,ng)%acquire) THEN
          BOUNDARY(ng) % hpondG_east = IniVal
        END IF
#   endif
#  endif
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isSig11,ng)%acquire) THEN
          BOUNDARY(ng) % sig11_east = IniVal
        END IF
#  ifndef ANA_SIG11OBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isSig11,ng)%acquire) THEN
          BOUNDARY(ng) % sig11G_east = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isSig22,ng)%acquire) THEN
          BOUNDARY(ng) % sig22_east = IniVal
        END IF
#  ifndef ANA_SIG22OBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isSig22,ng)%acquire) THEN
          BOUNDARY(ng) % sig22G_east = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isSig12,ng)%acquire) THEN
          BOUNDARY(ng) % sig12_east = IniVal
        END IF
#  ifndef ANA_SIG12OBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isSig12,ng)%acquire) THEN
          BOUNDARY(ng) % sig12G_east = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isUice,ng)%acquire) THEN
          BOUNDARY(ng) % ui_east = IniVal
        END IF
#  ifndef ANA_MIOBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isUice,ng)%acquire) THEN
          BOUNDARY(ng) % uiG_east = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isVice,ng)%acquire) THEN
          BOUNDARY(ng) % vi_east = IniVal
        END IF
#  ifndef ANA_MIOBC
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isVice,ng)%acquire) THEN
          BOUNDARY(ng) % viG_east = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isAice,ng)%acquire) THEN
          BOUNDARY(ng) % ai_south = IniVal
        END IF
#  ifndef ANA_AIOBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isAice,ng)%acquire) THEN
          BOUNDARY(ng) % aiG_south = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isHice,ng)%acquire) THEN
          BOUNDARY(ng) % hi_south = IniVal
        END IF
#  ifndef ANA_HIOBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isHice,ng)%acquire) THEN
          BOUNDARY(ng) % hiG_south = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isHsno,ng)%acquire) THEN
          BOUNDARY(ng) % hsn_south = IniVal
        END IF
#  ifndef ANA_HSNOBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isHsno,ng)%acquire) THEN
          BOUNDARY(ng) % hsnG_south = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isTice,ng)%acquire) THEN
          BOUNDARY(ng) % ti_south = IniVal
        END IF
#  ifndef ANA_TIOBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isTice,ng)%acquire) THEN
          BOUNDARY(ng) % tiG_south = IniVal
        END IF
#  endif
#  ifdef MELT_PONDS
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isApond,ng)%acquire) THEN
          BOUNDARY(ng) % apond_south = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isHpond,ng)%acquire) THEN
          BOUNDARY(ng) % hpond_south = IniVal
        END IF
#   ifndef ANA_SFWATOBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isApond,ng)%acquire) THEN
          BOUNDARY(ng) % apondG_south = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isHpond,ng)%acquire) THEN
          BOUNDARY(ng) % hpondG_south = IniVal
        END IF
#   endif
#  endif
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isSig11,ng)%acquire) THEN
          BOUNDARY(ng) % sig11_south = IniVal
        END IF
#  ifndef ANA_SIG11OBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isSig11,ng)%acquire) THEN
          BOUNDARY(ng) % sig11G_south = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isSig22,ng)%acquire) THEN
          BOUNDARY(ng) % sig22_south = IniVal
        END IF
#  ifndef ANA_SIG22OBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isSig22,ng)%acquire) THEN
          BOUNDARY(ng) % sig22G_south = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isSig12,ng)%acquire) THEN
          BOUNDARY(ng) % sig12_south = IniVal
        END IF
#  ifndef ANA_SIG12OBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isSig12,ng)%acquire) THEN
          BOUNDARY(ng) % sig12G_south = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isUice,ng)%acquire) THEN
          BOUNDARY(ng) % ui_south = IniVal
        END IF
#  ifndef ANA_MIOBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isUice,ng)%acquire) THEN
          BOUNDARY(ng) % uiG_south = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isVice,ng)%acquire) THEN
          BOUNDARY(ng) % vi_south = IniVal
        END IF
#  ifndef ANA_MIOBC
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isVice,ng)%acquire) THEN
          BOUNDARY(ng) % viG_south = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isAice,ng)%acquire) THEN
          BOUNDARY(ng) % ai_north = IniVal
        END IF
#  ifndef ANA_AIOBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isAice,ng)%acquire) THEN
          BOUNDARY(ng) % aiG_north = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isHice,ng)%acquire) THEN
          BOUNDARY(ng) % hi_north = IniVal
        END IF
#  ifndef ANA_HIOBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isHice,ng)%acquire) THEN
          BOUNDARY(ng) % hiG_north = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isHsno,ng)%acquire) THEN
          BOUNDARY(ng) % hsn_north = IniVal
        END IF
#  ifndef ANA_HSNOBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isHsno,ng)%acquire) THEN
          BOUNDARY(ng) % hsnG_north = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isTice,ng)%acquire) THEN
          BOUNDARY(ng) % ti_north = IniVal
        END IF
#  ifndef ANA_TIOBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isTice,ng)%acquire) THEN
          BOUNDARY(ng) % tiG_north = IniVal
        END IF
#  endif
#  ifdef MELT_PONDS
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isApond,ng)%acquire) THEN
          BOUNDARY(ng) % apond_north = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isHpond,ng)%acquire) THEN
          BOUNDARY(ng) % hpond_north = IniVal
        END IF
#   ifndef ANA_SFWATOBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isApond,ng)%acquire) THEN
          BOUNDARY(ng) % apondG_north = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isHpond,ng)%acquire) THEN
          BOUNDARY(ng) % hpondG_north = IniVal
        END IF
#   endif
#  endif
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isSig11,ng)%acquire) THEN
          BOUNDARY(ng) % sig11_north = IniVal
        END IF
#  ifndef ANA_SIG11OBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isSig11,ng)%acquire) THEN
          BOUNDARY(ng) % sig11G_north = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isSig22,ng)%acquire) THEN
          BOUNDARY(ng) % sig22_north = IniVal
        END IF
#  ifndef ANA_SIG22OBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isSig22,ng)%acquire) THEN
          BOUNDARY(ng) % sig22G_north = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isSig12,ng)%acquire) THEN
          BOUNDARY(ng) % sig12_north = IniVal
        END IF
#  ifndef ANA_SIG12OBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isSig12,ng)%acquire) THEN
          BOUNDARY(ng) % sig12G_north = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isUice,ng)%acquire) THEN
          BOUNDARY(ng) % ui_north = IniVal
        END IF
#  ifndef ANA_MIOBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isUice,ng)%acquire) THEN
          BOUNDARY(ng) % uiG_north = IniVal
        END IF
#  endif
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isVice,ng)%acquire) THEN
          BOUNDARY(ng) % vi_north = IniVal
        END IF
#  ifndef ANA_MIOBC
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isVice,ng)%acquire) THEN
          BOUNDARY(ng) % viG_north = IniVal
        END IF
#  endif
# endif
      END IF

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      tl_LBC(iwest,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % tl_zeta_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      tl_LBC(ieast,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % tl_zeta_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      tl_LBC(isouth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % tl_zeta_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      tl_LBC(inorth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % tl_zeta_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      tl_LBC(iwest,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % tl_ubar_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      tl_LBC(ieast,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % tl_ubar_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      tl_LBC(isouth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % tl_ubar_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      tl_LBC(inorth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % tl_ubar_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      tl_LBC(iwest,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % tl_vbar_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      tl_LBC(ieast,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % tl_vbar_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      tl_LBC(isouth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % tl_vbar_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      tl_LBC(inorth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % tl_vbar_north = IniVal
        END IF

# ifdef SOLVE3D
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      tl_LBC(iwest,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % tl_u_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      tl_LBC(ieast,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % tl_u_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      tl_LBC(isouth,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % tl_u_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      tl_LBC(inorth,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % tl_u_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      tl_LBC(iwest,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % tl_v_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      tl_LBC(ieast,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % tl_v_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      tl_LBC(isouth,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % tl_v_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      tl_LBC(inorth,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % tl_v_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      ANY(tl_LBC(iwest,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % tl_t_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      ANY(tl_LBC(ieast,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % tl_t_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      ANY(tl_LBC(isouth,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % tl_t_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      ANY(tl_LBC(inorth,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % tl_t_north = IniVal
        END IF
# endif
      END IF
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      ad_LBC(iwest,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % ad_zeta_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      ad_LBC(ieast,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % ad_zeta_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      ad_LBC(isouth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % ad_zeta_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      ad_LBC(inorth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % ad_zeta_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      ad_LBC(iwest,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ad_ubar_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      ad_LBC(ieast,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ad_ubar_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      ad_LBC(isouth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ad_ubar_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      ad_LBC(inorth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ad_ubar_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      ad_LBC(iwest,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % ad_vbar_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      ad_LBC(ieast,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % ad_vbar_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      ad_LBC(isouth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % ad_vbar_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      ad_LBC(inorth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % ad_vbar_north = IniVal
        END IF

# ifdef SOLVE3D
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      ad_LBC(iwest,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % ad_u_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      ad_LBC(ieast,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % ad_u_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      ad_LBC(isouth,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % ad_u_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      ad_LBC(inorth,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % ad_u_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      ad_LBC(iwest,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % ad_v_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      ad_LBC(ieast,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % ad_v_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      ad_LBC(isouth,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % ad_v_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      ad_LBC(inorth,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % ad_v_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      ANY(ad_LBC(iwest,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % ad_t_west = IniVal
        END IF

        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      ANY(ad_LBC(ieast,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % ad_t_east = IniVal
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      ANY(ad_LBC(isouth,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % ad_t_south = IniVal
        END IF

        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      ANY(ad_LBC(inorth,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % ad_t_north = IniVal
        END IF
# endif
      END IF
#endif

#ifdef ADJUST_BOUNDARY
!
!  Open boundaries arrays used in 4DVar adjustments.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
# ifdef SOLVE3D
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
# endif
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
        IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
# ifdef SOLVE3D
          BOUNDARY(ng) % ad_t_obc = IniVal
          BOUNDARY(ng) % ad_u_obc = IniVal
          BOUNDARY(ng) % ad_v_obc = IniVal
# endif
          BOUNDARY(ng) % ad_ubar_obc = IniVal
          BOUNDARY(ng) % ad_vbar_obc = IniVal
          BOUNDARY(ng) % ad_zeta_obc = IniVal
        END IF
      END IF

      IF ((model.eq.0).or.(model.eq.iTLM)) THEN
        IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
# ifdef SOLVE3D
          BOUNDARY(ng) % tl_t_obc = IniVal
          BOUNDARY(ng) % tl_u_obc = IniVal
          BOUNDARY(ng) % tl_v_obc = IniVal
# endif
          BOUNDARY(ng) % tl_ubar_obc = IniVal
          BOUNDARY(ng) % tl_vbar_obc = IniVal
          BOUNDARY(ng) % tl_zeta_obc = IniVal
        END IF
      END IF
#endif

      RETURN
      END SUBROUTINE initialize_boundary

      END MODULE mod_boundary
