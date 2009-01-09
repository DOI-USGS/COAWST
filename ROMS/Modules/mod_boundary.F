#include "cppdefs.h"
      MODULE mod_boundary
!
!svn $Id: mod_boundary.F 793 2008-10-15 13:57:18Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Open boundary conditions arrays:                                    !
!                                                                      !
!  t_east         Tracer (T units) eastern boundary conditions.        !
!  tG_east        Latest two-time snapshots of input tracer (Tunits)   !
!                   eastern boundary data.                             !
!  t_north        Tracer (T units) northern boundary conditions.       !
!  tG_north       Latest two-time snapshots of input tracer (Tunits)   !
!                   northern boundary data.                            !
!  t_south        Tracer (T units) southern boundary conditions.       !
!  tG_south       Latest two-time snapshots of input tracer (Tunits)   !
!                   southern boundary data.                            !
!  t_west         Tracer (T units) western boundary conditions.        !
!  tG_west        Latest two-time snapshots of input tracer (Tunits)   !
!                   western boundary data.                             !
!                                                                      !
!  ubar_east      2D u-momentum (m/s) eastern boundary conditions.     !
!  ubarG_east     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  ubar_north     2D u-momentum (m/s) northern boundary conditions.    !
!  ubarG_north    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  ubar_south     2D u-momentum (m/s) southern boundary conditions.    !
!  ubarG_south    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  ubar_west      2D u-momentum (m/s) western boundary conditions.     !
!  ubarG_west     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) western boundary data.                       !
!  vbar_east      2D v-momentum (m/s) eastern boundary conditions.     !
!  vbarG_east     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) eastern boundary data.                       !
!  vbar_north     2D v-momentum (m/s) northern boundary conditions.    !
!  vbarG_north    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) northern boundary data.                      !
!  vbar_south     2D v-momentum (m/s) southern boundary conditions.    !
!  vbarG_south    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) southern boundary data.                      !
!  vbar_west      2D v-momentum (m/s) western boundary conditions.     !
!  vbarG_west     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) western boundary data.                       !
!                                                                      !
!  u_east         3D u-momentum (m/s) eastern boundary conditions.     !
!  uG_east        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  u_north        3D u-momentum (m/s) northern boundary conditions.    !
!  uG_north       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  u_south        3D u-momentum (m/s) southern boundary conditions.    !
!  uG_south       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  u_west         3D u-momentum (m/s) western boundary conditions.     !
!  uG_west        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) western boundary data.                       !
!  v_east         3D v-momentum (m/s) eastern boundary conditions.     !
!  vG_east        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) eastern boundary data.                       !
!  v_north        3D v-momentum (m/s) northern boundary conditions.    !
!  vG_north       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) northern boundary data.                      !
!  v_south        3D v-momentum (m/s) southern boundary conditions.    !
!  vG_south       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) southern boundary data.                      !
!  v_west         3D v-momentum (m/s) western boundary conditions.     !
!  vG_west        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) western boundary data.                       !
!                                                                      !
!  zeta_east      Free-surface (m) eastern boundary conditions.        !
!  zetaG_east     Latest two-time snapshots of input free-surface (m)  !
!                   eastern boundary data.                             !
!  zeta_north     Free-surface (m) northern boundary conditions.       !
!  zetaG_north    Latest two-time snapshots of input free-surface (m)  !
!                   northern boundary data.                            !
!  zeta_south     Free-surface (m) southern boundary conditions.       !
!  zetaG_south    Latest two-time snapshots of input free-surface (m)  !
!                   southern boundary data.                            !
!  zeta_west      Free-surface (m) western boundary conditions.        !
!  zetaG_west     Latest two-time snapshots of input free-surface (m)  !
!                   western boundary data.                             !
!                                                                      !
!  Boundary arrays for radiation celerities needed for tangent linear  !
!  and adjoint models applications:                                    !
!                                                                      !
!  t_east_C2      Tracers eastern boundary, quadratic radition term    !
!                   term dTdx*dTdx+dTde*dTde (Tunits squared).         !
!  t_east_Cx      Tracers eastern boundary, celerity in the            !
!                   XI-direction (Tunits).                             !
!  t_east_Ce      Tracers eastern boundary, celerity in the            !
!                   ETA-direction (Tunits).                            !
!                                                                      !
!  t_north_C2     Tracers northern boundary, quadratic radition term   !
!                   dTdx*dTdx+dTde*dTde (Tunits squared).              !
!  t_north_Cx     Tracers northern boundary, celerity in the           !
!                   XI-direction (Tunits).                             !
!  t_north_Ce     Tracers northern boundary, celerity in the           !
!                   ETA-direction (Tunits).                            !
!                                                                      !
!  t_south_C2     Tracers southern boundary, quadratic radition term   !
!                   term dTdx*dTdx+dTde*dTde (Tunits squared).         !
!  t_south_Cx     Tracers southern boundary, celerity in the           !
!                   XI-direction (Tunits).                             !
!  t_south_Ce     Tracers northern boundary, celerity in the           !
!                   ETA-direction (Tunits).                            !
!                                                                      !
!  t_west_C2      Tracers western boundary, quadratic radition term    !
!                   dTdx*dTdx+dTde*dTde (Tunits squared).              !
!  t_west_Cx      Tracers western boundary, celerity in the            !
!                   XI-direction (Tunits).                             !
!  t_west_Ce      Tracers western boundary, celerity in the            !
!                   ETA-direction (Tunits).                            !
!                                                                      !
!  ubar_east_C2   2D U-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_east_Cx   2D U-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  ubar_east_Ce   2D U-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  ubar_north_C2  2D U-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_north_Cx  2D U-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  ubar_north_Ce  2D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  ubar_south_C2  2D U-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_south_Cx  2D U-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  ubar_south_Ce  2D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  ubar_west_C2   2D U-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  ubar_west_Cx   2D U-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  ubar_west_Ce   2D U-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  vbar_east_C2   2D V-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_east_Cx   2D V-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  vbar_east_Ce   2D V-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  vbar_north_C2  2D V-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_north_Cx  2D V-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  vbar_north_Ce  2D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  vbar_south_C2  2D V-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_south_Cx  2D V-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  vbar_south_Ce  2D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  vbar_west_C2   2D V-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  vbar_west_Cx   2D V-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  vbar_west_Ce   2D V-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!
!  u_east_C2      3D U-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_east_Cx      3D U-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  u_east_Ce      3D U-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  u_north_C2     3D U-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_north_Cx     3D U-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  u_north_Ce     3D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  u_south_C2     3D U-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_south_Cx     3D U-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  u_south_Ce     3D U-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  u_west_C2      3D U-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  u_west_Cx      3D U-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  u_west_Ce      3D U-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  v_east_C2      3D V-momentum eastern boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_east_Cx      3D V-momentum eastern boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  v_east_Ce      3D V-momentum eastern boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  v_north_C2     3D V-momentum northern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_north_Cx     3D V-momentum northern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  v_north_Ce     3D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  v_south_C2     3D V-momentum southern boundary, quadratic radition  !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_south_Cx     3D V-momentum southern boundary, celerity in the     !
!                   XI-direction (m/s).                                !
!  v_south_Ce     3D V-momentum northern boundary, celerity in the     !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  v_west_C2      3D V-momentum western boundary, quadratic radition   !
!                   term dUdx*dUdx+dUde*dUde (m2/s2).                  !
!  v_west_Cx      3D V-momentum western boundary, celerity in the      !
!                   XI-direction (m/s).                                !
!  v_west_Ce      3D V-momentum western boundary, celerity in the      !
!                   ETA-direction (m/s).                               !
!                                                                      !
!  zeta_east_C2   Free-surface eastern boundary, quadratic radition    !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_east_Cx   Free-surface eastern boundary, celerity in the       !
!                   XI-direction (m).                                  !
!  zeta_east_Ce   Free-surface eastern boundary, celerity in the       !
!                   ETA-direction (m).                                 !
!                                                                      !
!  zeta_north_C2  Free-surface northern boundary, quadratic radition   !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_north_Cx  Free-surface northern boundary, celerity in the      !
!                   XI-direction (m).                                  !
!  zeta_north_Ce  Free-surface northern boundary, celerity in the      !
!                   ETA-direction (m).                                 !
!                                                                      !
!  zeta_south_C2  Free-surface southern boundary, quadratic radition   !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_south_Cx  Free-surface southern boundary, celerity in the      !
!                   XI-direction (m).                                  !
!  zeta_south_Ce  Free-surface northern boundary, celerity in the      !
!                   ETA-direction (m).                                 !
!                                                                      !
!  zeta_west_C2   Free-surface western boundary, quadratic radition    !
!                   term dZdx*dZdx+dZde*dZde (m2).                     !
!  zeta_west_Cx   Free-surface western boundary, celerity in the       !
!                   XI-direction (m).                                  !
!  zeta_west_Ce   Free-surface western boundary, celerity in the       !
!                   ETA-direction (m).                                 !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_BOUNDARY

          real(r8), pointer :: zeta_west(:)
          real(r8), pointer :: zeta_east(:)
          real(r8), pointer :: zeta_south(:)
          real(r8), pointer :: zeta_north(:)

#ifndef ANA_FSOBC
          real(r8), pointer :: zetaG_west(:,:)
          real(r8), pointer :: zetaG_east(:,:)
          real(r8), pointer :: zetaG_south(:,:)
          real(r8), pointer :: zetaG_north(:,:)
#endif

#ifdef REFINED_GRID
          real(r8), pointer :: zetaR_west(:,:,:)
          real(r8), pointer :: zetaR_east(:,:,:)
          real(r8), pointer :: zetaR_south(:,:,:)
          real(r8), pointer :: zetaR_north(:,:,:)
!# if !defined SOLVE3D
          real(r8), pointer :: hR_west(:,:,:)
          real(r8), pointer :: hR_east(:,:,:)
          real(r8), pointer :: hR_south(:,:,:)
          real(r8), pointer :: hR_north(:,:,:)
!# endif
#endif

#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: zeta_west_C2(:)
          real(r8), pointer :: zeta_west_Ce(:)
          real(r8), pointer :: zeta_west_Cx(:)

          real(r8), pointer :: zeta_east_C2(:)
          real(r8), pointer :: zeta_east_Ce(:)
          real(r8), pointer :: zeta_east_Cx(:)

          real(r8), pointer :: zeta_south_C2(:)
          real(r8), pointer :: zeta_south_Ce(:)
          real(r8), pointer :: zeta_south_Cx(:)

          real(r8), pointer :: zeta_north_C2(:)
          real(r8), pointer :: zeta_north_Ce(:)
          real(r8), pointer :: zeta_north_Cx(:)
#endif

          real(r8), pointer :: ubar_west(:)
          real(r8), pointer :: vbar_west(:)
          real(r8), pointer :: ubar_east(:)
          real(r8), pointer :: vbar_east(:)
          real(r8), pointer :: ubar_south(:)
          real(r8), pointer :: vbar_south(:)
          real(r8), pointer :: ubar_north(:)
          real(r8), pointer :: vbar_north(:)

#ifndef ANA_M2OBC
          real(r8), pointer :: ubarG_west(:,:)
          real(r8), pointer :: vbarG_west(:,:)
          real(r8), pointer :: ubarG_east(:,:)
          real(r8), pointer :: vbarG_east(:,:)
          real(r8), pointer :: ubarG_south(:,:)
          real(r8), pointer :: vbarG_south(:,:)
          real(r8), pointer :: ubarG_north(:,:)
          real(r8), pointer :: vbarG_north(:,:)
#endif

#ifdef REFINED_GRID
          real(r8), pointer :: ubarR_west(:,:,:)
          real(r8), pointer :: vbarR_west(:,:,:)
          real(r8), pointer :: ubarR_east(:,:,:)
          real(r8), pointer :: vbarR_east(:,:,:)
          real(r8), pointer :: ubarR_south(:,:,:)
          real(r8), pointer :: vbarR_south(:,:,:)
          real(r8), pointer :: ubarR_north(:,:,:)
          real(r8), pointer :: vbarR_north(:,:,:)
          real(r8), pointer :: DUonR_west(:,:,:)
          real(r8), pointer :: DUonR_east(:,:,:)
          real(r8), pointer :: DUonR_south(:,:,:)
          real(r8), pointer :: DUonR_north(:,:,:)
          real(r8), pointer :: DVomR_west(:,:,:)
          real(r8), pointer :: DVomR_east(:,:,:)
          real(r8), pointer :: DVomR_north(:,:,:)
          real(r8), pointer :: DVomR_south(:,:,:)
#endif

#if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: ubar_west_C2(:)
          real(r8), pointer :: ubar_west_Ce(:)
          real(r8), pointer :: ubar_west_Cx(:)

          real(r8), pointer :: ubar_east_C2(:)
          real(r8), pointer :: ubar_east_Ce(:)
          real(r8), pointer :: ubar_east_Cx(:)

          real(r8), pointer :: ubar_south_C2(:)
          real(r8), pointer :: ubar_south_Ce(:)
          real(r8), pointer :: ubar_south_Cx(:)

          real(r8), pointer :: ubar_north_C2(:)
          real(r8), pointer :: ubar_north_Ce(:)
          real(r8), pointer :: ubar_north_Cx(:)

          real(r8), pointer :: vbar_west_C2(:)
          real(r8), pointer :: vbar_west_Ce(:)
          real(r8), pointer :: vbar_west_Cx(:)

          real(r8), pointer :: vbar_east_C2(:)
          real(r8), pointer :: vbar_east_Ce(:)
          real(r8), pointer :: vbar_east_Cx(:)

          real(r8), pointer :: vbar_south_C2(:)
          real(r8), pointer :: vbar_south_Ce(:)
          real(r8), pointer :: vbar_south_Cx(:)

          real(r8), pointer :: vbar_north_C2(:)
          real(r8), pointer :: vbar_north_Ce(:)
          real(r8), pointer :: vbar_north_Cx(:)
#endif

#ifdef SOLVE3D
          real(r8), pointer :: u_west(:,:)
          real(r8), pointer :: v_west(:,:)
          real(r8), pointer :: u_east(:,:)
          real(r8), pointer :: v_east(:,:)

# ifndef ANA_M3OBC
          real(r8), pointer :: uG_west(:,:,:)
          real(r8), pointer :: vG_west(:,:,:)
          real(r8), pointer :: uG_east(:,:,:)
          real(r8), pointer :: vG_east(:,:,:)
# endif

          real(r8), pointer :: u_south(:,:)
          real(r8), pointer :: v_south(:,:)
          real(r8), pointer :: u_north(:,:)
          real(r8), pointer :: v_north(:,:)

# ifndef ANA_M3OBC
          real(r8), pointer :: uG_south(:,:,:)
          real(r8), pointer :: vG_south(:,:,:)
          real(r8), pointer :: uG_north(:,:,:)
          real(r8), pointer :: vG_north(:,:,:)
# endif

#ifdef REFINED_GRID
          real(r8), pointer :: uR_west(:,:,:,:)
          real(r8), pointer :: vR_west(:,:,:,:)
          real(r8), pointer :: uR_east(:,:,:,:)
          real(r8), pointer :: vR_east(:,:,:,:)
          real(r8), pointer :: uR_south(:,:,:,:)
          real(r8), pointer :: vR_south(:,:,:,:)
          real(r8), pointer :: uR_north(:,:,:,:)
          real(r8), pointer :: vR_north(:,:,:,:)
#endif

# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: u_west_C2(:,:)
          real(r8), pointer :: u_west_Ce(:,:)
          real(r8), pointer :: u_west_Cx(:,:)

          real(r8), pointer :: u_east_C2(:,:)
          real(r8), pointer :: u_east_Ce(:,:)
          real(r8), pointer :: u_east_Cx(:,:)

          real(r8), pointer :: u_south_C2(:,:)
          real(r8), pointer :: u_south_Ce(:,:)
          real(r8), pointer :: u_south_Cx(:,:)

          real(r8), pointer :: u_north_C2(:,:)
          real(r8), pointer :: u_north_Ce(:,:)
          real(r8), pointer :: u_north_Cx(:,:)

          real(r8), pointer :: v_west_C2(:,:)
          real(r8), pointer :: v_west_Ce(:,:)
          real(r8), pointer :: v_west_Cx(:,:)

          real(r8), pointer :: v_east_C2(:,:)
          real(r8), pointer :: v_east_Ce(:,:)
          real(r8), pointer :: v_east_Cx(:,:)

          real(r8), pointer :: v_south_C2(:,:)
          real(r8), pointer :: v_south_Ce(:,:)
          real(r8), pointer :: v_south_Cx(:,:)

          real(r8), pointer :: v_north_C2(:,:)
          real(r8), pointer :: v_north_Ce(:,:)
          real(r8), pointer :: v_north_Cx(:,:)
# endif

          real(r8), pointer :: t_west(:,:,:)
          real(r8), pointer :: t_east(:,:,:)
          real(r8), pointer :: t_south(:,:,:)
          real(r8), pointer :: t_north(:,:,:)

# ifndef ANA_TOBC
          real(r8), pointer :: tG_west(:,:,:,:)
          real(r8), pointer :: tG_east(:,:,:,:)
          real(r8), pointer :: tG_south(:,:,:,:)
          real(r8), pointer :: tG_north(:,:,:,:)
# endif

#ifdef REFINED_GRID
          real(r8), pointer :: tR_west(:,:,:,:,:)
          real(r8), pointer :: tR_east(:,:,:,:,:)
          real(r8), pointer :: tR_south(:,:,:,:,:)
          real(r8), pointer :: tR_north(:,:,:,:,:)
#endif

# if defined CELERITY_READ || defined CELERITY_WRITE
          real(r8), pointer :: t_west_C2(:,:,:)
          real(r8), pointer :: t_west_Ce(:,:,:)
          real(r8), pointer :: t_west_Cx(:,:,:)

          real(r8), pointer :: t_east_C2(:,:,:)
          real(r8), pointer :: t_east_Ce(:,:,:)
          real(r8), pointer :: t_east_Cx(:,:,:)

          real(r8), pointer :: t_south_C2(:,:,:)
          real(r8), pointer :: t_south_Ce(:,:,:)
          real(r8), pointer :: t_south_Cx(:,:,:)

          real(r8), pointer :: t_north_C2(:,:,:)
          real(r8), pointer :: t_north_Ce(:,:,:)
          real(r8), pointer :: t_north_Cx(:,:,:)
# endif

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
!
!  Local variable declarations.
!
      integer, intent(in) :: ng
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( BOUNDARY(Ngrids) )
!
      allocate ( BOUNDARY(ng) % zeta_west(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_east(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_south(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_north(XI_DIM) )

#ifndef ANA_FSOBC
      allocate ( BOUNDARY(ng) % zetaG_west(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % zetaG_east(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % zetaG_south(XI_DIM,2) )
      allocate ( BOUNDARY(ng) % zetaG_north(XI_DIM,2) )
#endif

#ifdef REFINED_GRID
      allocate ( BOUNDARY(ng) % zetaR_west(5,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % zetaR_east(4,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % zetaR_south(XI_DIM,5,2) )
      allocate ( BOUNDARY(ng) % zetaR_north(XI_DIM,4,2) )
!# if !defined SOLVE3D
      allocate ( BOUNDARY(ng) % hR_west(5,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % hR_east(4,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % hR_south(XI_DIM,5,2) )
      allocate ( BOUNDARY(ng) % hR_north(XI_DIM,4,2) )
!# endif
#endif

#if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % zeta_west_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_west_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_west_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % zeta_east_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_east_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % zeta_east_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % zeta_south_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_south_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_south_Cx(XI_DIM) )

      allocate ( BOUNDARY(ng) % zeta_north_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_north_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % zeta_north_Cx(XI_DIM) )
#endif

      allocate ( BOUNDARY(ng) % ubar_west(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_east(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_south(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_north(XI_DIM) )

      allocate ( BOUNDARY(ng) % vbar_west(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north(XI_DIM) )

#ifndef ANA_M2OBC
      allocate ( BOUNDARY(ng) % ubarG_west(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % ubarG_east(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % ubarG_south(XI_DIM,2) )
      allocate ( BOUNDARY(ng) % ubarG_north(XI_DIM,2) )

      allocate ( BOUNDARY(ng) % vbarG_west(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_east(ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_south(XI_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarG_north(XI_DIM,2) )
#endif

#ifdef REFINED_GRID
      allocate ( BOUNDARY(ng) % ubarR_west(4,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % ubarR_east(3,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % ubarR_south(XI_DIM,4,2) )
      allocate ( BOUNDARY(ng) % ubarR_north(XI_DIM,3,2) )

      allocate ( BOUNDARY(ng) % vbarR_west(4,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarR_east(3,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % vbarR_south(XI_DIM,4,2) )
      allocate ( BOUNDARY(ng) % vbarR_north(XI_DIM,3,2) )

      allocate ( BOUNDARY(ng) % DUonR_west(4,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % DUonR_east(3,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % DUonR_south(XI_DIM,4,2) )
      allocate ( BOUNDARY(ng) % DUonR_north(XI_DIM,3,2) )

      allocate ( BOUNDARY(ng) % DVomR_west(4,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % DVomR_east(3,ETA_DIM,2) )
      allocate ( BOUNDARY(ng) % DVomR_south(XI_DIM,4,2) )
      allocate ( BOUNDARY(ng) % DVomR_north(XI_DIM,3,2) )
#endif

#if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % ubar_west_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_west_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_west_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % ubar_east_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_east_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % ubar_east_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % ubar_south_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_south_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_south_Cx(XI_DIM) )

      allocate ( BOUNDARY(ng) % ubar_north_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_north_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % ubar_north_Cx(XI_DIM) )

      allocate ( BOUNDARY(ng) % vbar_west_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_west_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_west_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % vbar_east_C2(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east_Ce(ETA_DIM) )
      allocate ( BOUNDARY(ng) % vbar_east_Cx(ETA_DIM) )

      allocate ( BOUNDARY(ng) % vbar_south_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_south_Cx(XI_DIM) )

      allocate ( BOUNDARY(ng) % vbar_north_C2(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north_Ce(XI_DIM) )
      allocate ( BOUNDARY(ng) % vbar_north_Cx(XI_DIM) )
#endif

#ifdef SOLVE3D

      allocate ( BOUNDARY(ng) % u_west(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_east(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_south(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_north(XI_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_west(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north(XI_DIM,N(ng)) )

# ifndef ANA_M3OBC
      allocate ( BOUNDARY(ng) % uG_west(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % uG_east(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % uG_south(XI_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % uG_north(XI_DIM,N(ng),2) )

      allocate ( BOUNDARY(ng) % vG_west(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_east(ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_south(XI_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vG_north(XI_DIM,N(ng),2) )
# endif

#ifdef REFINED_GRID
      allocate ( BOUNDARY(ng) % uR_west(4,ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % uR_east(3,ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % uR_south(XI_DIM,4,N(ng),2) )
      allocate ( BOUNDARY(ng) % uR_north(XI_DIM,3,N(ng),2) )

      allocate ( BOUNDARY(ng) % vR_west(4,ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vR_east(3,ETA_DIM,N(ng),2) )
      allocate ( BOUNDARY(ng) % vR_south(XI_DIM,4,N(ng),2) )
      allocate ( BOUNDARY(ng) % vR_north(XI_DIM,3,N(ng),2) )
#endif

# if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % u_west_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_west_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_west_Cx(ETA_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % u_east_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_east_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_east_Cx(ETA_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % u_south_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_south_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_south_Cx(XI_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % u_north_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_north_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % u_north_Cx(XI_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_west_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_west_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_west_Cx(ETA_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_east_C2(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east_Ce(ETA_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_east_Cx(ETA_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_south_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_south_Cx(XI_DIM,N(ng)) )

      allocate ( BOUNDARY(ng) % v_north_C2(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north_Ce(XI_DIM,N(ng)) )
      allocate ( BOUNDARY(ng) % v_north_Cx(XI_DIM,N(ng)) )
# endif

      allocate ( BOUNDARY(ng) % t_west(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_east(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_south(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_north(XI_DIM,N(ng),NT(ng)) )

# ifndef ANA_TOBC
      allocate ( BOUNDARY(ng) % tG_west(ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tG_east(ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tG_south(XI_DIM,N(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tG_north(XI_DIM,N(ng),2,NT(ng)) )
# endif

#ifdef REFINED_GRID
      allocate ( BOUNDARY(ng) % tR_west(4,ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tR_east(3,ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tR_south(XI_DIM,4,N(ng),2,NT(ng)) )
      allocate ( BOUNDARY(ng) % tR_north(XI_DIM,3,N(ng),2,NT(ng)) )
#endif

# if defined CELERITY_READ || defined CELERITY_WRITE
      allocate ( BOUNDARY(ng) % t_west_C2(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_west_Ce(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_west_Cx(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_east_C2(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_east_Ce(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_east_Cx(ETA_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_south_C2(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_south_Ce(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_south_Cx(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_north_C2(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_north_Ce(XI_DIM,N(ng),NT(ng)) )
      allocate ( BOUNDARY(ng) % t_north_Cx(XI_DIM,N(ng),NT(ng)) )
# endif
#endif

      RETURN
      END SUBROUTINE allocate_boundary

      SUBROUTINE initialize_boundary (ng, tile)
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
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: i, j
#ifdef SOLVE3D
      integer :: itrc, k
#endif
      real(r8), parameter :: IniVal = 0.0_r8

#ifndef DISTRIBUTE
# include "tile.h"
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Western boundary variables.
!
      IF (NORTH_WEST_TEST) THEN
        DO j=LOWER_BOUND_J,UPPER_BOUND_J
          BOUNDARY(ng) % zeta_west(j) = IniVal
#ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_west(j,1) = IniVal
          BOUNDARY(ng) % zetaG_west(j,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % zetaR_west(1,j,1) = IniVal
          BOUNDARY(ng) % zetaR_west(2,j,1) = IniVal
          BOUNDARY(ng) % zetaR_west(3,j,1) = IniVal
          BOUNDARY(ng) % zetaR_west(4,j,1) = IniVal
          BOUNDARY(ng) % zetaR_west(5,j,1) = IniVal
          BOUNDARY(ng) % zetaR_west(1,j,2) = IniVal
          BOUNDARY(ng) % zetaR_west(2,j,2) = IniVal
          BOUNDARY(ng) % zetaR_west(3,j,2) = IniVal
          BOUNDARY(ng) % zetaR_west(4,j,2) = IniVal
          BOUNDARY(ng) % zetaR_west(5,j,2) = IniVal
!# if !defined SOLVE3D
          BOUNDARY(ng) % hR_west(1,j,1) = IniVal
          BOUNDARY(ng) % hR_west(2,j,1) = IniVal
          BOUNDARY(ng) % hR_west(3,j,1) = IniVal
          BOUNDARY(ng) % hR_west(4,j,1) = IniVal
          BOUNDARY(ng) % hR_west(5,j,1) = IniVal
          BOUNDARY(ng) % hR_west(1,j,2) = IniVal
          BOUNDARY(ng) % hR_west(2,j,2) = IniVal
          BOUNDARY(ng) % hR_west(3,j,2) = IniVal
          BOUNDARY(ng) % hR_west(4,j,2) = IniVal
          BOUNDARY(ng) % hR_west(5,j,2) = IniVal
!# endif
#endif
          BOUNDARY(ng) % ubar_west(j) = IniVal
#ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_west(j,1) = IniVal
          BOUNDARY(ng) % ubarG_west(j,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % ubarR_west(1,j,1) = IniVal
          BOUNDARY(ng) % ubarR_west(2,j,1) = IniVal
          BOUNDARY(ng) % ubarR_west(3,j,1) = IniVal
          BOUNDARY(ng) % ubarR_west(4,j,1) = IniVal
          BOUNDARY(ng) % ubarR_west(1,j,2) = IniVal
          BOUNDARY(ng) % ubarR_west(2,j,2) = IniVal
          BOUNDARY(ng) % ubarR_west(3,j,2) = IniVal
          BOUNDARY(ng) % ubarR_west(4,j,2) = IniVal

          BOUNDARY(ng) % DUonR_west(1,j,1) = IniVal
          BOUNDARY(ng) % DUonR_west(2,j,1) = IniVal
          BOUNDARY(ng) % DUonR_west(3,j,1) = IniVal
          BOUNDARY(ng) % DUonR_west(4,j,1) = IniVal
          BOUNDARY(ng) % DUonR_west(1,j,2) = IniVal
          BOUNDARY(ng) % DUonR_west(2,j,2) = IniVal
          BOUNDARY(ng) % DUonR_west(3,j,2) = IniVal
          BOUNDARY(ng) % DUonR_west(4,j,2) = IniVal
#endif
          BOUNDARY(ng) % vbar_west(j) = IniVal
#ifndef ANA_M2OBC
          BOUNDARY(ng) % vbarG_west(j,1) = IniVal
          BOUNDARY(ng) % vbarG_west(j,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % vbarR_west(1,j,1) = IniVal
          BOUNDARY(ng) % vbarR_west(2,j,1) = IniVal
          BOUNDARY(ng) % vbarR_west(3,j,1) = IniVal
          BOUNDARY(ng) % vbarR_west(4,j,1) = IniVal
          BOUNDARY(ng) % vbarR_west(1,j,2) = IniVal
          BOUNDARY(ng) % vbarR_west(2,j,2) = IniVal
          BOUNDARY(ng) % vbarR_west(3,j,2) = IniVal
          BOUNDARY(ng) % vbarR_west(4,j,2) = IniVal

          BOUNDARY(ng) % DVomR_west(1,j,1) = IniVal
          BOUNDARY(ng) % DVomR_west(2,j,1) = IniVal
          BOUNDARY(ng) % DVomR_west(3,j,1) = IniVal
          BOUNDARY(ng) % DVomR_west(4,j,1) = IniVal
          BOUNDARY(ng) % DVomR_west(1,j,2) = IniVal
          BOUNDARY(ng) % DVomR_west(2,j,2) = IniVal
          BOUNDARY(ng) % DVomR_west(3,j,2) = IniVal
          BOUNDARY(ng) % DVomR_west(4,j,2) = IniVal
#endif
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_west_C2(j) = IniVal
          BOUNDARY(ng) % zeta_west_Ce(j) = IniVal
          BOUNDARY(ng) % zeta_west_Cx(j) = IniVal
          BOUNDARY(ng) % ubar_west_C2(j) = IniVal
          BOUNDARY(ng) % ubar_west_Ce(j) = IniVal
          BOUNDARY(ng) % ubar_west_Cx(j) = IniVal
          BOUNDARY(ng) % vbar_west_C2(j) = IniVal
          BOUNDARY(ng) % vbar_west_Ce(j) = IniVal
          BOUNDARY(ng) % vbar_west_Cx(j) = IniVal
#endif
        END DO
#ifdef SOLVE3D
        DO k=1,N(ng)
          DO j=LOWER_BOUND_J,UPPER_BOUND_J
            BOUNDARY(ng) % u_west(j,k) = IniVal
# ifndef ANA_M3OBC
            BOUNDARY(ng) % uG_west(j,k,1) = IniVal
            BOUNDARY(ng) % uG_west(j,k,2) = IniVal
# endif
#ifdef REFINED_GRID
            BOUNDARY(ng) % uR_west(1,j,k,1) = IniVal
            BOUNDARY(ng) % uR_west(2,j,k,1) = IniVal
            BOUNDARY(ng) % uR_west(3,j,k,1) = IniVal
            BOUNDARY(ng) % uR_west(4,j,k,1) = IniVal
            BOUNDARY(ng) % uR_west(1,j,k,2) = IniVal
            BOUNDARY(ng) % uR_west(2,j,k,2) = IniVal
            BOUNDARY(ng) % uR_west(3,j,k,2) = IniVal
            BOUNDARY(ng) % uR_west(4,j,k,2) = IniVal
#endif
            BOUNDARY(ng) % v_west(j,k) = IniVal
# ifndef ANA_M3OBC
            BOUNDARY(ng) % vG_west(j,k,1) = IniVal
            BOUNDARY(ng) % vG_west(j,k,2) = IniVal
# endif
#ifdef REFINED_GRID
            BOUNDARY(ng) % vR_west(1,j,k,1) = IniVal
            BOUNDARY(ng) % vR_west(2,j,k,1) = IniVal
            BOUNDARY(ng) % vR_west(3,j,k,1) = IniVal
            BOUNDARY(ng) % vR_west(4,j,k,1) = IniVal
            BOUNDARY(ng) % vR_west(1,j,k,2) = IniVal
            BOUNDARY(ng) % vR_west(2,j,k,2) = IniVal
            BOUNDARY(ng) % vR_west(3,j,k,2) = IniVal
            BOUNDARY(ng) % vR_west(4,j,k,2) = IniVal
#endif
# if defined CELERITY_READ || defined CELERITY_WRITE
            BOUNDARY(ng) % u_west_C2(j,k) = IniVal
            BOUNDARY(ng) % u_west_Ce(j,k) = IniVal
            BOUNDARY(ng) % u_west_Cx(j,k) = IniVal
            BOUNDARY(ng) % v_west_C2(j,k) = IniVal
            BOUNDARY(ng) % v_west_Ce(j,k) = IniVal
            BOUNDARY(ng) % v_west_Cx(j,k) = IniVal
# endif
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=LOWER_BOUND_J,UPPER_BOUND_J
              BOUNDARY(ng) % t_west(j,k,itrc) = IniVal
# ifndef ANA_TOBC
              BOUNDARY(ng) % tG_west(j,k,1,itrc) = IniVal
              BOUNDARY(ng) % tG_west(j,k,2,itrc) = IniVal
# endif
#ifdef REFINED_GRID
              BOUNDARY(ng) % tR_west(1,j,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_west(2,j,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_west(3,j,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_west(4,j,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_west(1,j,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_west(2,j,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_west(3,j,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_west(4,j,k,2,itrc) = IniVal
#endif
# if defined CELERITY_READ || defined CELERITY_WRITE
              BOUNDARY(ng) % t_west_C2(j,k,itrc) = IniVal
              BOUNDARY(ng) % t_west_Ce(j,k,itrc) = IniVal
              BOUNDARY(ng) % t_west_Cx(j,k,itrc) = IniVal
# endif
            END DO
          END DO
        END DO
#endif
      END IF
!
!  Eastern boundary variables.
!
      IF (SOUTH_EAST_TEST) THEN
        DO j=LOWER_BOUND_J,UPPER_BOUND_J
          BOUNDARY(ng) % zeta_east(j) = IniVal
#ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_east(j,1) = IniVal
          BOUNDARY(ng) % zetaG_east(j,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % zetaR_east(1,j,1) = IniVal
          BOUNDARY(ng) % zetaR_east(2,j,1) = IniVal
          BOUNDARY(ng) % zetaR_east(3,j,1) = IniVal
          BOUNDARY(ng) % zetaR_east(4,j,1) = IniVal
          BOUNDARY(ng) % zetaR_east(1,j,2) = IniVal
          BOUNDARY(ng) % zetaR_east(2,j,2) = IniVal
          BOUNDARY(ng) % zetaR_east(3,j,2) = IniVal
          BOUNDARY(ng) % zetaR_east(4,j,2) = IniVal
!# if !defined SOLVE3D
          BOUNDARY(ng) % hR_east(1,j,1) = IniVal
          BOUNDARY(ng) % hR_east(2,j,1) = IniVal
          BOUNDARY(ng) % hR_east(3,j,1) = IniVal
          BOUNDARY(ng) % hR_east(4,j,1) = IniVal
          BOUNDARY(ng) % hR_east(1,j,2) = IniVal
          BOUNDARY(ng) % hR_east(2,j,2) = IniVal
          BOUNDARY(ng) % hR_east(3,j,2) = IniVal
          BOUNDARY(ng) % hR_east(4,j,2) = IniVal
!# endif
#endif
          BOUNDARY(ng) % ubar_east(j) = IniVal
#ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_east(j,1) = IniVal
          BOUNDARY(ng) % ubarG_east(j,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % ubarR_east(1,j,1) = IniVal
          BOUNDARY(ng) % ubarR_east(2,j,1) = IniVal
          BOUNDARY(ng) % ubarR_east(3,j,1) = IniVal
          BOUNDARY(ng) % ubarR_east(1,j,2) = IniVal
          BOUNDARY(ng) % ubarR_east(2,j,2) = IniVal
          BOUNDARY(ng) % ubarR_east(3,j,2) = IniVal

          BOUNDARY(ng) % DUonR_east(1,j,1) = IniVal
          BOUNDARY(ng) % DUonR_east(2,j,1) = IniVal
          BOUNDARY(ng) % DUonR_east(3,j,1) = IniVal
          BOUNDARY(ng) % DUonR_east(1,j,2) = IniVal
          BOUNDARY(ng) % DUonR_east(2,j,2) = IniVal
          BOUNDARY(ng) % DUonR_east(3,j,2) = IniVal
#endif
          BOUNDARY(ng) % vbar_east(j) = IniVal
#ifndef ANA_M2OBC
          BOUNDARY(ng) % vbarG_east(j,1) = IniVal
          BOUNDARY(ng) % vbarG_east(j,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % vbarR_east(1,j,1) = IniVal
          BOUNDARY(ng) % vbarR_east(2,j,1) = IniVal
          BOUNDARY(ng) % vbarR_east(3,j,1) = IniVal
          BOUNDARY(ng) % vbarR_east(1,j,2) = IniVal
          BOUNDARY(ng) % vbarR_east(2,j,2) = IniVal
          BOUNDARY(ng) % vbarR_east(3,j,2) = IniVal

          BOUNDARY(ng) % DVomR_east(1,j,1) = IniVal
          BOUNDARY(ng) % DVomR_east(2,j,1) = IniVal
          BOUNDARY(ng) % DVomR_east(3,j,1) = IniVal
          BOUNDARY(ng) % DVomR_east(1,j,2) = IniVal
          BOUNDARY(ng) % DVomR_east(2,j,2) = IniVal
          BOUNDARY(ng) % DVomR_east(3,j,2) = IniVal
#endif
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_east_C2(j) = IniVal
          BOUNDARY(ng) % zeta_east_Ce(j) = IniVal
          BOUNDARY(ng) % zeta_east_Cx(j) = IniVal
          BOUNDARY(ng) % ubar_east_C2(j) = IniVal
          BOUNDARY(ng) % ubar_east_Ce(j) = IniVal
          BOUNDARY(ng) % ubar_east_Cx(j) = IniVal
          BOUNDARY(ng) % vbar_east_C2(j) = IniVal
          BOUNDARY(ng) % vbar_east_Ce(j) = IniVal
          BOUNDARY(ng) % vbar_east_Cx(j) = IniVal
#endif
        END DO
#ifdef SOLVE3D
        DO k=1,N(ng)
          DO j=LOWER_BOUND_J,UPPER_BOUND_J
            BOUNDARY(ng) % u_east(j,k) = IniVal
# ifndef ANA_M3OBC
            BOUNDARY(ng) % uG_east(j,k,1) = IniVal
            BOUNDARY(ng) % uG_east(j,k,2) = IniVal
# endif
#ifdef REFINED_GRID
            BOUNDARY(ng) % uR_east(1,j,k,1) = IniVal
            BOUNDARY(ng) % uR_east(2,j,k,1) = IniVal
            BOUNDARY(ng) % uR_east(3,j,k,1) = IniVal
            BOUNDARY(ng) % uR_east(1,j,k,2) = IniVal
            BOUNDARY(ng) % uR_east(2,j,k,2) = IniVal
            BOUNDARY(ng) % uR_east(3,j,k,2) = IniVal
#endif
            BOUNDARY(ng) % v_east(j,k) = IniVal
# ifndef ANA_M3OBC
            BOUNDARY(ng) % vG_east(j,k,1) = IniVal
            BOUNDARY(ng) % vG_east(j,k,2) = IniVal
# endif
#ifdef REFINED_GRID
            BOUNDARY(ng) % vR_east(1,j,k,1) = IniVal
            BOUNDARY(ng) % vR_east(2,j,k,1) = IniVal
            BOUNDARY(ng) % vR_east(3,j,k,1) = IniVal
            BOUNDARY(ng) % vR_east(1,j,k,2) = IniVal
            BOUNDARY(ng) % vR_east(2,j,k,2) = IniVal
            BOUNDARY(ng) % vR_east(3,j,k,2) = IniVal
#endif
# if defined CELERITY_READ || defined CELERITY_WRITE
            BOUNDARY(ng) % u_east_C2(j,k) = IniVal
            BOUNDARY(ng) % u_east_Ce(j,k) = IniVal
            BOUNDARY(ng) % u_east_Cx(j,k) = IniVal
            BOUNDARY(ng) % v_east_C2(j,k) = IniVal
            BOUNDARY(ng) % v_east_Ce(j,k) = IniVal
            BOUNDARY(ng) % v_east_Cx(j,k) = IniVal
# endif
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=LOWER_BOUND_J,UPPER_BOUND_J
              BOUNDARY(ng) % t_east(j,k,itrc) = IniVal
# ifndef ANA_TOBC
              BOUNDARY(ng) % tG_east(j,k,1,itrc) = IniVal
              BOUNDARY(ng) % tG_east(j,k,2,itrc) = IniVal
# endif
#ifdef REFINED_GRID
              BOUNDARY(ng) % tR_east(1,j,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_east(2,j,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_east(3,j,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_east(1,j,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_east(2,j,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_east(3,j,k,2,itrc) = IniVal
#endif
# if defined CELERITY_READ || defined CELERITY_WRITE
              BOUNDARY(ng) % t_east_C2(j,k,itrc) = IniVal
              BOUNDARY(ng) % t_east_Ce(j,k,itrc) = IniVal
              BOUNDARY(ng) % t_east_Cx(j,k,itrc) = IniVal
# endif
            END DO
          END DO
        END DO
#endif
      END IF
!
!  Southern boundary variables.
!
      IF (SOUTH_WEST_TEST) THEN
        DO i=LOWER_BOUND_I,UPPER_BOUND_I
          BOUNDARY(ng) % zeta_south(i) = IniVal
#ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_south(i,1) = IniVal
          BOUNDARY(ng) % zetaG_south(i,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % zetaR_south(i,1,1) = IniVal
          BOUNDARY(ng) % zetaR_south(i,2,1) = IniVal
          BOUNDARY(ng) % zetaR_south(i,3,1) = IniVal
          BOUNDARY(ng) % zetaR_south(i,4,1) = IniVal
          BOUNDARY(ng) % zetaR_south(i,5,1) = IniVal
          BOUNDARY(ng) % zetaR_south(i,1,2) = IniVal
          BOUNDARY(ng) % zetaR_south(i,2,2) = IniVal
          BOUNDARY(ng) % zetaR_south(i,3,2) = IniVal
          BOUNDARY(ng) % zetaR_south(i,4,2) = IniVal
          BOUNDARY(ng) % zetaR_south(i,5,2) = IniVal
!# if !defined SOLVE3D
          BOUNDARY(ng) % hR_south(i,1,1) = IniVal
          BOUNDARY(ng) % hR_south(i,2,1) = IniVal
          BOUNDARY(ng) % hR_south(i,3,1) = IniVal
          BOUNDARY(ng) % hR_south(i,4,1) = IniVal
          BOUNDARY(ng) % hR_south(i,5,1) = IniVal
          BOUNDARY(ng) % hR_south(i,1,2) = IniVal
          BOUNDARY(ng) % hR_south(i,2,2) = IniVal
          BOUNDARY(ng) % hR_south(i,3,2) = IniVal
          BOUNDARY(ng) % hR_south(i,4,2) = IniVal
          BOUNDARY(ng) % hR_south(i,5,2) = IniVal
!# endif
#endif
          BOUNDARY(ng) % ubar_south(i) = IniVal
#ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_south(i,1) = IniVal
          BOUNDARY(ng) % ubarG_south(i,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % ubarR_south(i,1,1) = IniVal
          BOUNDARY(ng) % ubarR_south(i,2,1) = IniVal
          BOUNDARY(ng) % ubarR_south(i,3,1) = IniVal
          BOUNDARY(ng) % ubarR_south(i,4,1) = IniVal
          BOUNDARY(ng) % ubarR_south(i,1,2) = IniVal
          BOUNDARY(ng) % ubarR_south(i,2,2) = IniVal
          BOUNDARY(ng) % ubarR_south(i,3,2) = IniVal
          BOUNDARY(ng) % ubarR_south(i,4,2) = IniVal

          BOUNDARY(ng) % DUonR_south(i,1,1) = IniVal
          BOUNDARY(ng) % DUonR_south(i,2,1) = IniVal
          BOUNDARY(ng) % DUonR_south(i,3,1) = IniVal
          BOUNDARY(ng) % DUonR_south(i,4,1) = IniVal
          BOUNDARY(ng) % DUonR_south(i,1,2) = IniVal
          BOUNDARY(ng) % DUonR_south(i,2,2) = IniVal
          BOUNDARY(ng) % DUonR_south(i,3,2) = IniVal
          BOUNDARY(ng) % DUonR_south(i,4,2) = IniVal
#endif
          BOUNDARY(ng) % vbar_south(i) = IniVal
#ifndef ANA_M2OBC
          BOUNDARY(ng) % vbarG_south(i,1) = IniVal
          BOUNDARY(ng) % vbarG_south(i,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % vbarR_south(i,1,1) = IniVal
          BOUNDARY(ng) % vbarR_south(i,2,1) = IniVal
          BOUNDARY(ng) % vbarR_south(i,3,1) = IniVal
          BOUNDARY(ng) % vbarR_south(i,4,1) = IniVal
          BOUNDARY(ng) % vbarR_south(i,1,2) = IniVal
          BOUNDARY(ng) % vbarR_south(i,2,2) = IniVal
          BOUNDARY(ng) % vbarR_south(i,3,2) = IniVal
          BOUNDARY(ng) % vbarR_south(i,4,2) = IniVal

          BOUNDARY(ng) % DVomR_south(i,1,1) = IniVal
          BOUNDARY(ng) % DVomR_south(i,2,1) = IniVal
          BOUNDARY(ng) % DVomR_south(i,3,1) = IniVal
          BOUNDARY(ng) % DVomR_south(i,4,1) = IniVal
          BOUNDARY(ng) % DVomR_south(i,1,2) = IniVal
          BOUNDARY(ng) % DVomR_south(i,2,2) = IniVal
          BOUNDARY(ng) % DVomR_south(i,3,2) = IniVal
          BOUNDARY(ng) % DVomR_south(i,4,2) = IniVal
#endif
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_south_C2(i) = IniVal
          BOUNDARY(ng) % zeta_south_Ce(i) = IniVal
          BOUNDARY(ng) % zeta_south_Cx(i) = IniVal
          BOUNDARY(ng) % ubar_south_C2(i) = IniVal
          BOUNDARY(ng) % ubar_south_Ce(i) = IniVal
          BOUNDARY(ng) % ubar_south_Cx(i) = IniVal
          BOUNDARY(ng) % vbar_south_C2(i) = IniVal
          BOUNDARY(ng) % vbar_south_Ce(i) = IniVal
          BOUNDARY(ng) % vbar_south_Cx(i) = IniVal
#endif
        END DO
#ifdef SOLVE3D
        DO k=1,N(ng)
          DO i=LOWER_BOUND_I,UPPER_BOUND_I
            BOUNDARY(ng) % u_south(i,k) = IniVal
# ifndef ANA_M3OBC
            BOUNDARY(ng) % uG_south(i,k,1) = IniVal
            BOUNDARY(ng) % uG_south(i,k,2) = IniVal
# endif
#ifdef REFINED_GRID
            BOUNDARY(ng) % uR_south(i,1,k,1) = IniVal
            BOUNDARY(ng) % uR_south(i,2,k,1) = IniVal
            BOUNDARY(ng) % uR_south(i,3,k,1) = IniVal
            BOUNDARY(ng) % uR_south(i,4,k,1) = IniVal
            BOUNDARY(ng) % uR_south(i,1,k,2) = IniVal
            BOUNDARY(ng) % uR_south(i,2,k,2) = IniVal
            BOUNDARY(ng) % uR_south(i,3,k,2) = IniVal
            BOUNDARY(ng) % uR_south(i,4,k,2) = IniVal
#endif
            BOUNDARY(ng) % v_south(i,k) = IniVal
# ifndef ANA_M3OBC
            BOUNDARY(ng) % vG_south(i,k,1) = IniVal
            BOUNDARY(ng) % vG_south(i,k,2) = IniVal
# endif
#ifdef REFINED_GRID
            BOUNDARY(ng) % vR_south(i,1,k,1) = IniVal
            BOUNDARY(ng) % vR_south(i,2,k,1) = IniVal
            BOUNDARY(ng) % vR_south(i,3,k,1) = IniVal
            BOUNDARY(ng) % vR_south(i,4,k,1) = IniVal
            BOUNDARY(ng) % vR_south(i,1,k,2) = IniVal
            BOUNDARY(ng) % vR_south(i,2,k,2) = IniVal
            BOUNDARY(ng) % vR_south(i,3,k,2) = IniVal
            BOUNDARY(ng) % vR_south(i,4,k,2) = IniVal
#endif
# if defined CELERITY_READ || defined CELERITY_WRITE
            BOUNDARY(ng) % u_south_C2(i,k) = IniVal
            BOUNDARY(ng) % u_south_Ce(i,k) = IniVal
            BOUNDARY(ng) % u_south_Cx(i,k) = IniVal
            BOUNDARY(ng) % v_south_C2(i,k) = IniVal
            BOUNDARY(ng) % v_south_Ce(i,k) = IniVal
            BOUNDARY(ng) % v_south_Cx(i,k) = IniVal
# endif
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=LOWER_BOUND_I,UPPER_BOUND_I
              BOUNDARY(ng) % t_south(i,k,itrc) = IniVal
# ifndef ANA_TOBC
              BOUNDARY(ng) % tG_south(i,k,1,itrc) = IniVal
              BOUNDARY(ng) % tG_south(i,k,2,itrc) = IniVal
# endif
#ifdef REFINED_GRID
              BOUNDARY(ng) % tR_south(i,1,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_south(i,2,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_south(i,3,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_south(i,4,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_south(i,1,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_south(i,2,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_south(i,3,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_south(i,4,k,2,itrc) = IniVal
#endif
# if defined CELERITY_READ || defined CELERITY_WRITE
              BOUNDARY(ng) % t_south_C2(i,k,itrc) = IniVal
              BOUNDARY(ng) % t_south_Ce(i,k,itrc) = IniVal
              BOUNDARY(ng) % t_south_Cx(i,k,itrc) = IniVal
# endif
            END DO
          END DO
        END DO
#endif
      END IF
!
!  Northern boundary variables.
!
      IF (NORTH_EAST_TEST) THEN
        DO i=LOWER_BOUND_I,UPPER_BOUND_I
          BOUNDARY(ng) % zeta_north(i) = IniVal
#ifndef ANA_FSOBC
          BOUNDARY(ng) % zetaG_north(i,1) = IniVal
          BOUNDARY(ng) % zetaG_north(i,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % zetaR_north(i,1,1) = IniVal
          BOUNDARY(ng) % zetaR_north(i,2,1) = IniVal
          BOUNDARY(ng) % zetaR_north(i,3,1) = IniVal
          BOUNDARY(ng) % zetaR_north(i,4,1) = IniVal
          BOUNDARY(ng) % zetaR_north(i,1,2) = IniVal
          BOUNDARY(ng) % zetaR_north(i,2,2) = IniVal
          BOUNDARY(ng) % zetaR_north(i,3,2) = IniVal
          BOUNDARY(ng) % zetaR_north(i,4,2) = IniVal
!# if !defined SOLVE3D
          BOUNDARY(ng) % hR_north(i,1,1) = IniVal
          BOUNDARY(ng) % hR_north(i,2,1) = IniVal
          BOUNDARY(ng) % hR_north(i,3,1) = IniVal
          BOUNDARY(ng) % hR_north(i,4,1) = IniVal
          BOUNDARY(ng) % hR_north(i,1,2) = IniVal
          BOUNDARY(ng) % hR_north(i,2,2) = IniVal
          BOUNDARY(ng) % hR_north(i,3,2) = IniVal
          BOUNDARY(ng) % hR_north(i,4,2) = IniVal
!# endif
#endif
          BOUNDARY(ng) % ubar_north(i) = IniVal
#ifndef ANA_M2OBC
          BOUNDARY(ng) % ubarG_north(i,1) = IniVal
          BOUNDARY(ng) % ubarG_north(i,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % ubarR_north(i,1,1) = IniVal
          BOUNDARY(ng) % ubarR_north(i,2,1) = IniVal
          BOUNDARY(ng) % ubarR_north(i,3,1) = IniVal
          BOUNDARY(ng) % ubarR_north(i,1,2) = IniVal
          BOUNDARY(ng) % ubarR_north(i,2,2) = IniVal
          BOUNDARY(ng) % ubarR_north(i,3,2) = IniVal

          BOUNDARY(ng) % DUonR_north(i,1,1) = IniVal
          BOUNDARY(ng) % DUonR_north(i,2,1) = IniVal
          BOUNDARY(ng) % DUonR_north(i,3,1) = IniVal
          BOUNDARY(ng) % DUonR_north(i,1,2) = IniVal
          BOUNDARY(ng) % DUonR_north(i,2,2) = IniVal
          BOUNDARY(ng) % DUonR_north(i,3,2) = IniVal
#endif
          BOUNDARY(ng) % vbar_north(i) = IniVal
#ifndef ANA_M2OBC
          BOUNDARY(ng) % vbarG_north(i,1) = IniVal
          BOUNDARY(ng) % vbarG_north(i,2) = IniVal
#endif
#ifdef REFINED_GRID
          BOUNDARY(ng) % vbarR_north(i,1,1) = IniVal
          BOUNDARY(ng) % vbarR_north(i,2,1) = IniVal
          BOUNDARY(ng) % vbarR_north(i,3,1) = IniVal
          BOUNDARY(ng) % vbarR_north(i,1,2) = IniVal
          BOUNDARY(ng) % vbarR_north(i,2,2) = IniVal
          BOUNDARY(ng) % vbarR_north(i,3,2) = IniVal

          BOUNDARY(ng) % DVomR_north(i,1,1) = IniVal
          BOUNDARY(ng) % DVomR_north(i,2,1) = IniVal
          BOUNDARY(ng) % DVomR_north(i,3,1) = IniVal
          BOUNDARY(ng) % DVomR_north(i,1,2) = IniVal
          BOUNDARY(ng) % DVomR_north(i,2,2) = IniVal
          BOUNDARY(ng) % DVomR_north(i,3,2) = IniVal
#endif
#if defined CELERITY_READ || defined CELERITY_WRITE
          BOUNDARY(ng) % zeta_north_C2(i) = IniVal
          BOUNDARY(ng) % zeta_north_Ce(i) = IniVal
          BOUNDARY(ng) % zeta_north_Cx(i) = IniVal
          BOUNDARY(ng) % ubar_north_C2(i) = IniVal
          BOUNDARY(ng) % ubar_north_Ce(i) = IniVal
          BOUNDARY(ng) % ubar_north_Cx(i) = IniVal
          BOUNDARY(ng) % vbar_north_C2(i) = IniVal
          BOUNDARY(ng) % vbar_north_Ce(i) = IniVal
          BOUNDARY(ng) % vbar_north_Cx(i) = IniVal
#endif
        END DO
#ifdef SOLVE3D
        DO k=1,N(ng)
          DO i=LOWER_BOUND_I,UPPER_BOUND_I
            BOUNDARY(ng) % u_north(i,k) = IniVal
# ifndef ANA_M3OBC
            BOUNDARY(ng) % uG_north(i,k,1) = IniVal
            BOUNDARY(ng) % uG_north(i,k,2) = IniVal
# endif
#ifdef REFINED_GRID
            BOUNDARY(ng) % uR_north(i,1,k,1) = IniVal
            BOUNDARY(ng) % uR_north(i,2,k,1) = IniVal
            BOUNDARY(ng) % uR_north(i,3,k,1) = IniVal
            BOUNDARY(ng) % uR_north(i,1,k,2) = IniVal
            BOUNDARY(ng) % uR_north(i,2,k,2) = IniVal
            BOUNDARY(ng) % uR_north(i,3,k,2) = IniVal
#endif
            BOUNDARY(ng) % v_north(i,k) = IniVal
# ifndef ANA_M3OBC
            BOUNDARY(ng) % vG_north(i,k,1) = IniVal
            BOUNDARY(ng) % vG_north(i,k,2) = IniVal
# endif
#ifdef REFINED_GRID
            BOUNDARY(ng) % vR_north(i,1,k,1) = IniVal
            BOUNDARY(ng) % vR_north(i,2,k,1) = IniVal
            BOUNDARY(ng) % vR_north(i,3,k,1) = IniVal
            BOUNDARY(ng) % vR_north(i,1,k,2) = IniVal
            BOUNDARY(ng) % vR_north(i,2,k,2) = IniVal
            BOUNDARY(ng) % vR_north(i,3,k,2) = IniVal
#endif
# if defined CELERITY_READ || defined CELERITY_WRITE
            BOUNDARY(ng) % u_north_C2(i,k) = IniVal
            BOUNDARY(ng) % u_north_Ce(i,k) = IniVal
            BOUNDARY(ng) % u_north_Cx(i,k) = IniVal
            BOUNDARY(ng) % v_north_C2(i,k) = IniVal
            BOUNDARY(ng) % v_north_Ce(i,k) = IniVal
            BOUNDARY(ng) % v_north_Cx(i,k) = IniVal
# endif
          END DO
        END DO
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=LOWER_BOUND_I,UPPER_BOUND_I
              BOUNDARY(ng) % t_north(i,k,itrc) = IniVal
# ifndef ANA_TOBC
              BOUNDARY(ng) % tG_north(i,k,1,itrc) = IniVal
              BOUNDARY(ng) % tG_north(i,k,2,itrc) = IniVal
# endif
#ifdef REFINED_GRID
              BOUNDARY(ng) % tR_north(i,1,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_north(i,2,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_north(i,3,k,1,itrc) = IniVal
              BOUNDARY(ng) % tR_north(i,1,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_north(i,2,k,2,itrc) = IniVal
              BOUNDARY(ng) % tR_north(i,3,k,2,itrc) = IniVal
#endif
# if defined CELERITY_READ || defined CELERITY_WRITE
              BOUNDARY(ng) % t_north_C2(i,k,itrc) = IniVal
              BOUNDARY(ng) % t_north_Ce(i,k,itrc) = IniVal
              BOUNDARY(ng) % t_north_Cx(i,k,itrc) = IniVal
# endif
            END DO
          END DO
        END DO
#endif
      END IF
      RETURN
      END SUBROUTINE initialize_boundary

      END MODULE mod_boundary
