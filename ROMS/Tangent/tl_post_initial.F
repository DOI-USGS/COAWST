#include "cppdefs.h"
      MODULE tl_post_initial_mod
#ifdef TANGENT
!
!git $Id$
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  On the first timestep, it computes the initial depths and level     !
!  thicknesses from the initial free-surface field. Additionally, it   !
!  initializes the tangent linear state variables for all time levels  !
!  and applies lateral boundary conditions.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef NESTING
      USE mod_scalars
# endif
!
      USE ini_fields_mod,    ONLY : set_zeta_timeavg
# ifdef SOLVE3D
      USE set_depth_mod,     ONLY : set_depth
# endif
      USE tl_ini_fields_mod, ONLY : tl_ini_fields, tl_ini_zeta
# ifdef SOLVE3D
      USE tl_set_depth_mod,  ONLY : tl_set_depth
# endif
!
      implicit none
!
      PUBLIC  :: tl_post_initial
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_post_initial (ng, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: tile
!
!-----------------------------------------------------------------------
!  Initialize free-surface and compute initial level thicknesses and
!  depths. Here we need nonlinear 'Zt_avg1' and 'Hz' to compute initial
!  'tl_ubar' and 'tl_vbar'.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_zeta_timeavg (ng, tile, iNLM)
        CALL set_depth (ng, tile, iNLM)
        CALL tl_ini_zeta (ng, tile, model)
        CALL tl_set_depth (ng, tile, model)
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Initialize other state variables.
!-----------------------------------------------------------------------
!
      DO tile=last_tile(ng),first_tile(ng),-1
        CALL tl_ini_fields (ng, tile, model)
      END DO
!$OMP BARRIER

# ifdef NESTING
!
!-----------------------------------------------------------------------
!  Extract donor grid initial data at contact points and store it in
!  REFINED structure so it can be used for the space-time interpolation.
!-----------------------------------------------------------------------
!
      IF (RefinedGrid(ng)) THEN
        CALL tl_nesting (ng, model, ngetD)
      END IF
# endif
!
      RETURN
      END SUBROUTINE tl_post_initial
#endif
      END MODULE tl_post_initial_mod
