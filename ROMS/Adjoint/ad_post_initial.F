#include "cppdefs.h"
      MODULE ad_post_initial_mod
#ifdef ADJOINT
!
!git $Id$
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  On the last timestep, it computes the adjoint of initial depths and !
!  level thicknesses from the initial free-surface field. Additionally,!
!  it initializes the adjoint state variables for all time levels and  !
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
      USE ad_ini_fields_mod, ONLY : ad_ini_fields, ad_ini_zeta
# ifdef SOLVE3D
      USE ad_set_depth_mod,  ONLY : ad_set_depth
# endif
!
      implicit none
!
      PUBLIC  :: ad_post_initial
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_post_initial (ng, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: tile

# ifdef NESTING
!
!-----------------------------------------------------------------------
!  Extract donor grid initial data at contact points and store it in
!  REFINED structure so it can be used for the space-time interpolation.
!-----------------------------------------------------------------------
!
      IF (RefinedGrid(ng)) THEN
        CALL ad_nesting (ng, model, ngetD)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Initialize other state variables.
!-----------------------------------------------------------------------
!
      DO tile=last_tile(ng),first_tile(ng),-1
        CALL ad_ini_fields (ng, tile, model)
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Initialize free-surface and compute initial level thicknesses and
!  depths. Here we need nonlinear 'Zt_avg1' and 'Hz' to compute initial
!  'ad_ubar' and 'ad_vbar'.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL set_zeta_timeavg (ng, tile, iNLM)
        CALL set_depth (ng, tile, iNLM)
        CALL ad_set_depth (ng, tile, model)
        CALL ad_ini_zeta (ng, tile, model)
      END DO
!$OMP BARRIER
!
      RETURN
      END SUBROUTINE ad_post_initial
#endif
      END MODULE ad_post_initial_mod
