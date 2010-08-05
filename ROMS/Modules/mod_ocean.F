#include "cppdefs.h"
      MODULE mod_ocean
!
!svn $Id: mod_ocean.F 787 2008-10-14 00:37:50Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  2D Primitive Variables.                                             !
!                                                                      !
!  rubar        Right-hand-side of 2D U-momentum equation (m4/s2).     !
!  rvbar        Right-hand-side of 2D V-momentum equation (m4/s2).     !
!  rzeta        Right-hand-side of free surface equation (m3/s).       !
!  ubar         Vertically integrated U-momentum component (m/s).      !
!  vbar         Vertically integrated V-momentum component (m/s).      !
!  zeta         Free surface (m).                                      !
!                                                                      !
!  3D Primitive Variables.                                             !
!                                                                      !
!  pden         Potential Density anomaly (kg/m3).                     !
!  rho          Density anomaly (kg/m3).                               !
!  ru           Right-hand-side of 3D U-momentum equation (m4/s2).     !
!  rv           Right hand side of 3D V-momentum equation (m4/s2).     !
!  t            Tracer type variables (active and passive).            !
!  u            3D U-momentum component (m/s).                         !
!  v            3D V-momentum component (m/s).                         !
!  W            S-coordinate (omega*Hz/mn) vertical velocity (m3/s).   !
!                                                                      !
!  Biology Variables.                                                  !
!                                                                      !
!  pH           Surface concentration of hydrogen ions.                !
!                                                                      !
#ifdef NEARSHORE
!  Nearshore radiation stresses:                                       !
!                                                                      !
!  rulag2d      2D U-Stokes tendency term (m4/s2).                     !
!  rvlag2d      2D V-Stokes tendency term (m4/s2).                     !
!  ubar_stokes  2D U-Stokes drift velocity (m/s).                      !
!  vbar_stokes  2D V-Stokes drift velocity (m/s).                      !
!  rulag3d      3D U-Stokes tendency term (m4/s2).                     !
!  rvlag3d      3D V-Stokes tendency term (m4/s2).                     !
!  u_stokes     3D U-Stokes drift velocity (m/s).                      !
!  v_stokes     3D V-Stokes drift velocity (m/s).                      !
!  w_stokes     3D W-Stokes drift velocity (m/s).                      !
!                                                                      !
#endif
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_OCEAN
!
!  Nonlinear model state.
!
          real(r8), pointer :: rubar(:,:,:)
          real(r8), pointer :: rvbar(:,:,:)
          real(r8), pointer :: rzeta(:,:,:)
          real(r8), pointer :: ubar(:,:,:)
          real(r8), pointer :: vbar(:,:,:)
          real(r8), pointer :: zeta(:,:,:)
#if defined NEARSHORE_WEC
          real(r8), pointer :: zetat(:,:)
          real(r8), pointer :: zetaw(:,:)
#endif
#if defined NEARSHORE
          real(r8), pointer :: rulag2d(:,:)
          real(r8), pointer :: rvlag2d(:,:)
          real(r8), pointer :: ubar_stokes(:,:)
          real(r8), pointer :: vbar_stokes(:,:)
#endif
#ifdef SOLVE3D
          real(r8), pointer :: pden(:,:,:)
          real(r8), pointer :: rho(:,:,:)
          real(r8), pointer :: ru(:,:,:,:)
          real(r8), pointer :: rv(:,:,:,:)
          real(r8), pointer :: t(:,:,:,:,:)
          real(r8), pointer :: u(:,:,:,:)
          real(r8), pointer :: v(:,:,:,:)
          real(r8), pointer :: W(:,:,:)
          real(r8), pointer :: wvel(:,:,:)
# if defined NEARSHORE
          real(r8), pointer :: rulag3d(:,:,:)
          real(r8), pointer :: rvlag3d(:,:,:)
          real(r8), pointer :: u_stokes(:,:,:)
          real(r8), pointer :: v_stokes(:,:,:)
          real(r8), pointer :: w_stokes(:,:,:)
# endif
# ifdef UV_KIRBY
          real(r8), pointer :: uwave(:,:)
          real(r8), pointer :: vwave(:,:)
# endif
# if defined BIO_FENNEL && defined CARBON
          real(r8), pointer :: pH(:,:)
# endif
#endif

#if defined REFINED_GRID
          real(r8), pointer :: hindpi_u(:,:)
          real(r8), pointer :: hindpj_u(:,:)
          real(r8), pointer :: hindpi_v(:,:)
          real(r8), pointer :: hindpj_v(:,:)
          real(r8), pointer :: hindpi_r(:,:)
          real(r8), pointer :: hindpj_r(:,:)
          integer, pointer :: hindpr(:)
          integer, pointer :: hindpu(:)
          integer, pointer :: hindpv(:)
          real(r8), pointer :: hindci_u(:,:)
          real(r8), pointer :: hindcj_u(:,:)
          real(r8), pointer :: hindci_v(:,:)
          real(r8), pointer :: hindcj_v(:,:)
          real(r8), pointer :: hindci_r(:,:)
          real(r8), pointer :: hindcj_r(:,:)
          integer, pointer :: hindcr(:)
          integer, pointer :: hindcu(:)
          integer, pointer :: hindcv(:)
#endif
#if defined COMPOSED_GRID
          real(r8), pointer :: hind2di_u(:,:)
          real(r8), pointer :: hind2dj_u(:,:)
          real(r8), pointer :: hind2di_v(:,:)
          real(r8), pointer :: hind2dj_v(:,:)
          real(r8), pointer :: hind2di_r(:,:)
          real(r8), pointer :: hind2dj_r(:,:)
          integer, pointer :: hindr(:)
          integer, pointer :: hindu(:)
          integer, pointer :: hindv(:)
          logical, pointer :: hind_init(:)
          real(r8), pointer :: hind3di_u(:,:)
          real(r8), pointer :: hind3dj_u(:,:)
          real(r8), pointer :: hind3di_v(:,:)
          real(r8), pointer :: hind3dj_v(:,:)
          real(r8), pointer :: hind3di_r(:,:)
          real(r8), pointer :: hind3dj_r(:,:)
          real(r8), pointer :: hind4di_u(:,:)
          real(r8), pointer :: hind4dj_u(:,:)
          real(r8), pointer :: hind4di_v(:,:)
          real(r8), pointer :: hind4dj_v(:,:)
          real(r8), pointer :: hind4di_r(:,:)
          real(r8), pointer :: hind4dj_r(:,:)
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
          real(r8), pointer :: tl_rubar(:,:,:)
          real(r8), pointer :: tl_rvbar(:,:,:)
          real(r8), pointer :: tl_rzeta(:,:,:)
          real(r8), pointer :: tl_ubar(:,:,:)
          real(r8), pointer :: tl_vbar(:,:,:)
          real(r8), pointer :: tl_zeta(:,:,:)
# if defined NEARSHORE_MELLOR
          real(r8), pointer :: tl_rulag2d(:,:)
          real(r8), pointer :: tl_rvlag2d(:,:)
          real(r8), pointer :: tl_ubar_stokes(:,:)
          real(r8), pointer :: tl_vbar_stokes(:,:)
# endif
# ifdef SOLVE3D
          real(r8), pointer :: tl_pden(:,:,:)
          real(r8), pointer :: tl_rho(:,:,:)
          real(r8), pointer :: tl_ru(:,:,:,:)
          real(r8), pointer :: tl_rv(:,:,:,:)
          real(r8), pointer :: tl_t(:,:,:,:,:)
          real(r8), pointer :: tl_u(:,:,:,:)
          real(r8), pointer :: tl_v(:,:,:,:)
          real(r8), pointer :: tl_W(:,:,:)
#  if defined NEARSHORE_MELLOR
          real(r8), pointer :: tl_rulag3d(:,:,:)
          real(r8), pointer :: tl_rvlag3d(:,:,:)
          real(r8), pointer :: tl_u_stokes(:,:,:)
          real(r8), pointer :: tl_v_stokes(:,:,:)
#  endif
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
          real(r8), pointer :: ad_rubar(:,:,:)
          real(r8), pointer :: ad_rvbar(:,:,:)
          real(r8), pointer :: ad_rzeta(:,:,:)
          real(r8), pointer :: ad_ubar(:,:,:)
          real(r8), pointer :: ad_vbar(:,:,:)
          real(r8), pointer :: ad_zeta(:,:,:)
          real(r8), pointer :: ad_ubar_sol(:,:)
          real(r8), pointer :: ad_vbar_sol(:,:)
          real(r8), pointer :: ad_zeta_sol(:,:)
# if defined NEARSHORE_MELLOR
          real(r8), pointer :: ad_rulag2d(:,:)
          real(r8), pointer :: ad_rvlag2d(:,:)
          real(r8), pointer :: ad_ubar_stokes(:,:)
          real(r8), pointer :: ad_vbar_stokes(:,:)
# endif
# ifdef SOLVE3D
          real(r8), pointer :: ad_pden(:,:,:)
          real(r8), pointer :: ad_rho(:,:,:)
          real(r8), pointer :: ad_ru(:,:,:,:)
          real(r8), pointer :: ad_rv(:,:,:,:)
          real(r8), pointer :: ad_t(:,:,:,:,:)
          real(r8), pointer :: ad_u(:,:,:,:)
          real(r8), pointer :: ad_v(:,:,:,:)
          real(r8), pointer :: ad_W(:,:,:)
#  if defined NEARSHORE_MELLOR
          real(r8), pointer :: ad_rulag3d(:,:,:)
          real(r8), pointer :: ad_rvlag3d(:,:,:)
          real(r8), pointer :: ad_u_stokes(:,:,:)
          real(r8), pointer :: ad_v_stokes(:,:,:)
#  endif
# endif
#endif

#if defined FOUR_DVAR || defined IMPULSE
!
!  Working arrays to store adjoint impulse forcing, error covariance,
!  standard deviations, or descent conjugate vectors (directions).
!
          real(r8), pointer :: b_ubar(:,:,:)
          real(r8), pointer :: b_vbar(:,:,:)
          real(r8), pointer :: b_zeta(:,:,:)
# ifdef SOLVE3D
          real(r8), pointer :: b_t(:,:,:,:,:)
          real(r8), pointer :: b_u(:,:,:,:)
          real(r8), pointer :: b_v(:,:,:,:)
# endif
# ifdef FOUR_DVAR
          real(r8), pointer :: d_ubar(:,:)
          real(r8), pointer :: d_vbar(:,:)
          real(r8), pointer :: d_zeta(:,:)
#  ifdef SOLVE3D
          real(r8), pointer :: d_t(:,:,:,:)
          real(r8), pointer :: d_u(:,:,:)
          real(r8), pointer :: d_v(:,:,:)
#  endif
          real(r8), pointer :: e_ubar(:,:,:)
          real(r8), pointer :: e_vbar(:,:,:)
          real(r8), pointer :: e_zeta(:,:,:)
#  ifdef SOLVE3D
          real(r8), pointer :: e_t(:,:,:,:,:)
          real(r8), pointer :: e_u(:,:,:,:)
          real(r8), pointer :: e_v(:,:,:,:)
#  endif
# endif
# ifdef WEAK_CONSTRAINT
          real(r8), pointer :: f_ubar(:,:)
          real(r8), pointer :: f_vbar(:,:)
          real(r8), pointer :: f_zeta(:,:)
#  ifdef SOLVE3D
          real(r8), pointer :: f_t(:,:,:,:)
          real(r8), pointer :: f_u(:,:,:)
          real(r8), pointer :: f_v(:,:,:)
#  endif
# endif
#endif

#if defined FORWARD_READ && \
   (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
# ifdef FORWARD_RHS
          real(r8), pointer :: rubarG(:,:,:)
          real(r8), pointer :: rvbarG(:,:,:)
          real(r8), pointer :: rzetaG(:,:,:)
# endif
          real(r8), pointer :: ubarG(:,:,:)
          real(r8), pointer :: vbarG(:,:,:)
          real(r8), pointer :: zetaG(:,:,:)
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
          real(r8), pointer :: ruG(:,:,:,:)
          real(r8), pointer :: rvG(:,:,:,:)
#  endif
          real(r8), pointer :: tG(:,:,:,:,:)
          real(r8), pointer :: uG(:,:,:,:)
          real(r8), pointer :: vG(:,:,:,:)
# endif
# ifdef WEAK_CONSTRAINT
          real(r8), pointer :: f_zetaG(:,:,:)
#  ifdef SOLVE3D
          real(r8), pointer :: f_tG(:,:,:,:,:)
          real(r8), pointer :: f_uG(:,:,:,:)
          real(r8), pointer :: f_vG(:,:,:,:)
#  else
          real(r8), pointer :: f_ubarG(:,:,:)
          real(r8), pointer :: f_vbarG(:,:,:)
#  endif
# endif
#endif

        END TYPE T_OCEAN

        TYPE (T_OCEAN), allocatable :: OCEAN(:)

      CONTAINS

      SUBROUTINE allocate_ocean (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( OCEAN(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( OCEAN(ng) % rubar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rvbar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rzeta(LBi:UBi,LBj:UBj,2) )

      allocate ( OCEAN(ng) % ubar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % vbar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % zeta(LBi:UBi,LBj:UBj,3) )
#if defined NEARSHORE_WEC
      allocate ( OCEAN(ng) % zetat(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % zetaw(LBi:UBi,LBj:UBj) )
#endif
#if defined NEARSHORE
      allocate ( OCEAN(ng) % rulag2d(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % rvlag2d(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % ubar_stokes(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % vbar_stokes(LBi:UBi,LBj:UBj) )
#endif

#ifdef SOLVE3D
      allocate ( OCEAN(ng) % pden(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % rho(LBi:UBi,LBj:UBj,N(ng)) )

      allocate ( OCEAN(ng) % ru(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % rv(LBi:UBi,LBj:UBj,0:N(ng),2) )

      allocate ( OCEAN(ng) % t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng)) )
      allocate ( OCEAN(ng) % u(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % v(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % W(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( OCEAN(ng) % wvel(LBi:UBi,LBj:UBj,0:N(ng)) )

# if defined NEARSHORE
      allocate ( OCEAN(ng) % rulag3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % rvlag3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % u_stokes(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % v_stokes(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % w_stokes(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif
# ifdef UV_KIRBY
      allocate ( OCEAN(ng) % uwave(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % vwave(LBi:UBi,LBj:UBj) )
# endif

# if defined BIO_FENNEL && defined CARBON
      allocate ( OCEAN(ng) % pH(LBi:UBi,LBj:UBj) )
# endif
#endif


#if defined REFINED_GRID
          allocate ( OCEAN(ng) % hindpi_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindpj_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindpi_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindpj_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindpi_r(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindpj_r(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindpr(1:4) )
          allocate ( OCEAN(ng) % hindpu(1:4) )
          allocate ( OCEAN(ng) % hindpv(1:4) )
          allocate ( OCEAN(ng) % hindci_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindcj_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindci_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindcj_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindci_r(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindcj_r(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindcr(1:4) )
          allocate ( OCEAN(ng) % hindcu(1:4) )
          allocate ( OCEAN(ng) % hindcv(1:4) )
#endif
#ifdef COMPOSED_GRID
          allocate ( OCEAN(ng) % hind2di_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind2dj_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind2di_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind2dj_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind2di_r(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind2dj_r(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hindr(1:4) )
          allocate ( OCEAN(ng) % hindu(1:4) )
          allocate ( OCEAN(ng) % hindv(1:4) )
          allocate ( OCEAN(ng) % hind_init(1:12) )
          allocate ( OCEAN(ng) % hind3di_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind3dj_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind3di_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind3dj_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind3di_r(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind3dj_r(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind4di_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind4dj_u(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind4di_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind4dj_v(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind4di_r(LBi:UBi,LBj:UBj) )
          allocate ( OCEAN(ng) % hind4dj_r(LBi:UBi,LBj:UBj) )
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      allocate ( OCEAN(ng) % tl_rubar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % tl_rvbar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % tl_rzeta(LBi:UBi,LBj:UBj,2) )

      allocate ( OCEAN(ng) % tl_ubar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % tl_vbar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % tl_zeta(LBi:UBi,LBj:UBj,3) )

# if defined NEARSHORE_MELLOR
      allocate ( OCEAN(ng) % tl_rulag2d(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % tl_rvlag2d(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % tl_ubar_stokes(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % tl_vbar_stokes(LBi:UBi,LBj:UBj) )
# endif

# ifdef SOLVE3D
      allocate ( OCEAN(ng) % tl_pden(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % tl_rho(LBi:UBi,LBj:UBj,N(ng)) )

      allocate ( OCEAN(ng) % tl_ru(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % tl_rv(LBi:UBi,LBj:UBj,0:N(ng),2) )

      allocate ( OCEAN(ng) % tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng)) )
      allocate ( OCEAN(ng) % tl_u(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % tl_v(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % tl_W(LBi:UBi,LBj:UBj,0:N(ng)) )

#  if defined NEARSHORE_MELLOR
      allocate ( OCEAN(ng) % tl_rulag3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % tl_rvlag3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % tl_u_stokes(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % tl_v_stokes(LBi:UBi,LBj:UBj,N(ng)) )
#  endif
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      allocate ( OCEAN(ng) % ad_rubar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % ad_rvbar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % ad_rzeta(LBi:UBi,LBj:UBj,2) )

      allocate ( OCEAN(ng) % ad_ubar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % ad_vbar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % ad_zeta(LBi:UBi,LBj:UBj,3) )

      allocate ( OCEAN(ng) % ad_ubar_sol(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % ad_vbar_sol(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % ad_zeta_sol(LBi:UBi,LBj:UBj) )

# if defined NEARSHORE_MELLOR
      allocate ( OCEAN(ng) % ad_rulag2d(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % ad_rvlag2d(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % ad_ubar_stokes(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % ad_vbar_stokes(LBi:UBi,LBj:UBj) )
# endif

# ifdef SOLVE3D
      allocate ( OCEAN(ng) % ad_pden(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % ad_rho(LBi:UBi,LBj:UBj,N(ng)) )

      allocate ( OCEAN(ng) % ad_ru(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % ad_rv(LBi:UBi,LBj:UBj,0:N(ng),2) )

      allocate ( OCEAN(ng) % ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng)) )
      allocate ( OCEAN(ng) % ad_u(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % ad_v(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % ad_W(LBi:UBi,LBj:UBj,0:N(ng)) )

#  if defined NEARSHORE_MELLOR
      allocate ( OCEAN(ng) % ad_rulag3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % ad_rvlag3d(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % ad_u_stokes(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % ad_v_stokes(LBi:UBi,LBj:UBj,N(ng)) )
#  endif
# endif
#endif

#if defined FOUR_DVAR || defined IMPULSE
!
!  Working arrays to store adjoint impulse forcing, background error
!  covariance, background-error standard deviations, or descent
!  conjugate vectors (directions).
!
      allocate ( OCEAN(ng) % b_ubar(LBi:UBi,LBj:UBj,NSA) )
      allocate ( OCEAN(ng) % b_vbar(LBi:UBi,LBj:UBj,NSA) )
      allocate ( OCEAN(ng) % b_zeta(LBi:UBi,LBj:UBj,NSA) )

# ifdef SOLVE3D
      allocate ( OCEAN(ng) % b_t(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng)) )
      allocate ( OCEAN(ng) % b_u(LBi:UBi,LBj:UBj,N(ng),NSA) )
      allocate ( OCEAN(ng) % b_v(LBi:UBi,LBj:UBj,N(ng),NSA) )
# endif

# ifdef FOUR_DVAR
      allocate ( OCEAN(ng) % d_ubar(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % d_vbar(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % d_zeta(LBi:UBi,LBj:UBj) )

#  ifdef SOLVE3D
      allocate ( OCEAN(ng) % d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      allocate ( OCEAN(ng) % d_u(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % d_v(LBi:UBi,LBj:UBj,N(ng)) )
#  endif

      allocate ( OCEAN(ng) % e_ubar(LBi:UBi,LBj:UBj,NSA) )
      allocate ( OCEAN(ng) % e_vbar(LBi:UBi,LBj:UBj,NSA) )
      allocate ( OCEAN(ng) % e_zeta(LBi:UBi,LBj:UBj,NSA) )

#  ifdef SOLVE3D
      allocate ( OCEAN(ng) % e_t(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng)) )
      allocate ( OCEAN(ng) % e_u(LBi:UBi,LBj:UBj,N(ng),NSA) )
      allocate ( OCEAN(ng) % e_v(LBi:UBi,LBj:UBj,N(ng),NSA) )
#  endif

#  ifdef WEAK_CONSTRAINT
      allocate ( OCEAN(ng) % f_ubar(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % f_vbar(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % f_zeta(LBi:UBi,LBj:UBj) )

#   ifdef SOLVE3D
      allocate ( OCEAN(ng) % f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      allocate ( OCEAN(ng) % f_u(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % f_v(LBi:UBi,LBj:UBj,N(ng)) )
#   endif
#  endif
# endif
#endif

#if defined FORWARD_READ && \
   (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
# ifdef FORWARD_RHS
      allocate ( OCEAN(ng) % rubarG(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rvbarG(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rzetaG(LBi:UBi,LBj:UBj,2) )
# endif
      allocate ( OCEAN(ng) % ubarG(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % vbarG(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % zetaG(LBi:UBi,LBj:UBj,2) )

# ifdef SOLVE3D
#  ifdef FORWARD_RHS
      allocate ( OCEAN(ng) % ruG(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % rvG(LBi:UBi,LBj:UBj,0:N(ng),2) )
#  endif
      allocate ( OCEAN(ng) % tG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)) )
      allocate ( OCEAN(ng) % uG(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % vG(LBi:UBi,LBj:UBj,N(ng),2) )
# endif
# ifdef WEAK_CONSTRAINT
      allocate ( OCEAN(ng) % f_zetaG(LBi:UBi,LBj:UBj,2) )
#  ifdef SOLVE3D
      allocate ( OCEAN(ng) % f_tG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)) )
      allocate ( OCEAN(ng) % f_uG(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % f_vG(LBi:UBi,LBj:UBj,N(ng),2) )
#  else
      allocate ( OCEAN(ng) % f_ubarG(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % f_vbarG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif
#endif

      RETURN
      END SUBROUTINE allocate_ocean

      SUBROUTINE initialize_ocean (ng, tile, model)
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
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, rec
#ifdef SOLVE3D
      integer :: itrc, k
#endif

      real(r8), parameter :: IniVal = 0.0_r8

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
#if defined REFINED_GRID
            OCEAN(ng) % hindpr(1) = IniVal
            OCEAN(ng) % hindpr(2) = IniVal
            OCEAN(ng) % hindpr(3) = IniVal
            OCEAN(ng) % hindpr(4) = IniVal
            OCEAN(ng) % hindpu(1) = IniVal
            OCEAN(ng) % hindpu(2) = IniVal
            OCEAN(ng) % hindpu(3) = IniVal
            OCEAN(ng) % hindpu(4) = IniVal
            OCEAN(ng) % hindpv(1) = IniVal
            OCEAN(ng) % hindpv(2) = IniVal
            OCEAN(ng) % hindpv(3) = IniVal
            OCEAN(ng) % hindpv(4) = IniVal
            OCEAN(ng) % hindcr(1) = IniVal
            OCEAN(ng) % hindcr(2) = IniVal
            OCEAN(ng) % hindcr(3) = IniVal
            OCEAN(ng) % hindcr(4) = IniVal
            OCEAN(ng) % hindcu(1) = IniVal
            OCEAN(ng) % hindcu(2) = IniVal
            OCEAN(ng) % hindcu(3) = IniVal
            OCEAN(ng) % hindcu(4) = IniVal
            OCEAN(ng) % hindcv(1) = IniVal
            OCEAN(ng) % hindcv(2) = IniVal
            OCEAN(ng) % hindcv(3) = IniVal
            OCEAN(ng) % hindcv(4) = IniVal
#endif
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            OCEAN(ng) % rubar(i,j,1) = IniVal
            OCEAN(ng) % rubar(i,j,2) = IniVal
            OCEAN(ng) % rvbar(i,j,1) = IniVal
            OCEAN(ng) % rvbar(i,j,2) = IniVal
            OCEAN(ng) % rzeta(i,j,1) = IniVal
            OCEAN(ng) % rzeta(i,j,2) = IniVal

            OCEAN(ng) % ubar(i,j,1) = IniVal
            OCEAN(ng) % ubar(i,j,2) = IniVal
            OCEAN(ng) % ubar(i,j,3) = IniVal
            OCEAN(ng) % vbar(i,j,1) = IniVal
            OCEAN(ng) % vbar(i,j,2) = IniVal
            OCEAN(ng) % vbar(i,j,3) = IniVal
            OCEAN(ng) % zeta(i,j,1) = IniVal
            OCEAN(ng) % zeta(i,j,2) = IniVal
            OCEAN(ng) % zeta(i,j,3) = IniVal
#if defined NEARSHORE_WEC
            OCEAN(ng) % zetat(i,j) = IniVal
            OCEAN(ng) % zetaw(i,j) = IniVal
#endif
#if defined NEARSHORE
            OCEAN(ng) % rulag2d(i,j) = IniVal
            OCEAN(ng) % rvlag2d(i,j) = IniVal
            OCEAN(ng) % ubar_stokes(i,j) = IniVal
            OCEAN(ng) % vbar_stokes(i,j) = IniVal
#endif
#if defined BIO_FENNEL && defined CARBON && defined SOLVE3D
            OCEAN(ng) % pH(i,j) = 8.0_r8
#endif
#ifdef UV_KIRBY
            OCEAN(ng) % uwave(i,j) = IniVal
            OCEAN(ng) % vwave(i,j) = IniVal
#endif
          END DO
#ifdef SOLVE3D
          DO k=1,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % pden(i,j,k) = IniVal
              OCEAN(ng) % rho(i,j,k) = IniVal

              OCEAN(ng) % u(i,j,k,1) = IniVal
              OCEAN(ng) % u(i,j,k,2) = IniVal
              OCEAN(ng) % v(i,j,k,1) = IniVal
              OCEAN(ng) % v(i,j,k,2) = IniVal
# if defined NEARSHORE
              OCEAN(ng) % rulag3d(i,j,k) = IniVal
              OCEAN(ng) % rvlag3d(i,j,k) = IniVal
              OCEAN(ng) % u_stokes(i,j,k) = IniVal
              OCEAN(ng) % v_stokes(i,j,k) = IniVal
# endif
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % ru(i,j,k,1) = IniVal
              OCEAN(ng) % ru(i,j,k,2) = IniVal
              OCEAN(ng) % rv(i,j,k,1) = IniVal
              OCEAN(ng) % rv(i,j,k,2) = IniVal

              OCEAN(ng) % W(i,j,k) = IniVal
              OCEAN(ng) % wvel(i,j,k) = IniVal
# if defined NEARSHORE
              OCEAN(ng) % w_stokes(i,j,k) = IniVal
# endif
            END DO
          END DO
          DO itrc=1,NT(ng)
            DO k=1,N(ng)
              DO i=Imin,Imax
                OCEAN(ng) % t(i,j,k,1,itrc) = IniVal
                OCEAN(ng) % t(i,j,k,2,itrc) = IniVal
                OCEAN(ng) % t(i,j,k,3,itrc) = IniVal
              END DO
            END DO
          END DO 
# if defined REFINED_GRID
          DO i=Imin,Imax
            OCEAN(ng) % hindpi_u(i,j) = IniVal
            OCEAN(ng) % hindpj_u(i,j) = IniVal
            OCEAN(ng) % hindpi_v(i,j) = IniVal
            OCEAN(ng) % hindpj_v(i,j) = IniVal
            OCEAN(ng) % hindpi_r(i,j) = IniVal
            OCEAN(ng) % hindpj_r(i,j) = IniVal
            OCEAN(ng) % hindci_u(i,j) = IniVal
            OCEAN(ng) % hindcj_u(i,j) = IniVal
            OCEAN(ng) % hindci_v(i,j) = IniVal
            OCEAN(ng) % hindcj_v(i,j) = IniVal
            OCEAN(ng) % hindci_r(i,j) = IniVal
            OCEAN(ng) % hindcj_r(i,j) = IniVal
          END DO
# endif
# ifdef COMPOSED_GRID
          OCEAN(ng) % hind_init(1:12)=.FALSE.
          DO i=Imin,Imax
            OCEAN(ng) % hind2di_u(i,j) = IniVal
            OCEAN(ng) % hind2dj_u(i,j) = IniVal
            OCEAN(ng) % hind2di_v(i,j) = IniVal
            OCEAN(ng) % hind2dj_v(i,j) = IniVal
            OCEAN(ng) % hind2di_r(i,j) = IniVal
            OCEAN(ng) % hind2dj_r(i,j) = IniVal
            OCEAN(ng) % hind3di_u(i,j) = IniVal
            OCEAN(ng) % hind3dj_u(i,j) = IniVal
            OCEAN(ng) % hind3di_v(i,j) = IniVal
            OCEAN(ng) % hind3dj_v(i,j) = IniVal
            OCEAN(ng) % hind3di_r(i,j) = IniVal
            OCEAN(ng) % hind3dj_r(i,j) = IniVal
            OCEAN(ng) % hind4di_u(i,j) = IniVal
            OCEAN(ng) % hind4dj_u(i,j) = IniVal
            OCEAN(ng) % hind4di_v(i,j) = IniVal
            OCEAN(ng) % hind4dj_v(i,j) = IniVal
            OCEAN(ng) % hind4di_r(i,j) = IniVal
            OCEAN(ng) % hind4dj_r(i,j) = IniVal
          END DO
# endif
#endif
        END DO
      END IF

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            OCEAN(ng) % tl_rubar(i,j,1) = IniVal
            OCEAN(ng) % tl_rubar(i,j,2) = IniVal
            OCEAN(ng) % tl_rvbar(i,j,1) = IniVal
            OCEAN(ng) % tl_rvbar(i,j,2) = IniVal
            OCEAN(ng) % tl_rzeta(i,j,1) = IniVal
            OCEAN(ng) % tl_rzeta(i,j,2) = IniVal

            OCEAN(ng) % tl_ubar(i,j,1) = IniVal
            OCEAN(ng) % tl_ubar(i,j,2) = IniVal
            OCEAN(ng) % tl_ubar(i,j,3) = IniVal
            OCEAN(ng) % tl_vbar(i,j,1) = IniVal
            OCEAN(ng) % tl_vbar(i,j,2) = IniVal
            OCEAN(ng) % tl_vbar(i,j,3) = IniVal
            OCEAN(ng) % tl_zeta(i,j,1) = IniVal
            OCEAN(ng) % tl_zeta(i,j,2) = IniVal
            OCEAN(ng) % tl_zeta(i,j,3) = IniVal
# if defined NEARSHORE_MELLOR
            OCEAN(ng) % tl_rulag2d(i,j) = IniVal
            OCEAN(ng) % tl_rvlag2d(i,j) = IniVal
            OCEAN(ng) % tl_ubar_stokes(i,j) = IniVal
            OCEAN(ng) % tl_vbar_stokes(i,j) = IniVal
# endif
          END DO
# ifdef SOLVE3D
          DO k=1,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % tl_pden(i,j,k) = IniVal
              OCEAN(ng) % tl_rho(i,j,k) = IniVal

              OCEAN(ng) % tl_u(i,j,k,1) = IniVal
              OCEAN(ng) % tl_u(i,j,k,2) = IniVal
              OCEAN(ng) % tl_v(i,j,k,1) = IniVal
              OCEAN(ng) % tl_v(i,j,k,2) = IniVal
#  if defined NEARSHORE_MELLOR
              OCEAN(ng) % tl_rulag3d(i,j,k) = IniVal
              OCEAN(ng) % tl_rvlag3d(i,j,k) = IniVal
              OCEAN(ng) % tl_u_stokes(i,j,k) = IniVal
              OCEAN(ng) % tl_v_stokes(i,j,k) = IniVal
#  endif
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % tl_ru(i,j,k,1) = IniVal
              OCEAN(ng) % tl_ru(i,j,k,2) = IniVal
              OCEAN(ng) % tl_rv(i,j,k,1) = IniVal
              OCEAN(ng) % tl_rv(i,j,k,2) = IniVal

              OCEAN(ng) % tl_W(i,j,k) = IniVal
            END DO
          END DO
          DO itrc=1,NT(ng)
            DO k=1,N(ng)
              DO i=Imin,Imax
                OCEAN(ng) % tl_t(i,j,k,1,itrc) = IniVal
                OCEAN(ng) % tl_t(i,j,k,2,itrc) = IniVal
                OCEAN(ng) % tl_t(i,j,k,3,itrc) = IniVal
              END DO
            END DO
          END DO
# endif
        END DO
      END IF
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            OCEAN(ng) % ad_rubar(i,j,1) = IniVal
            OCEAN(ng) % ad_rubar(i,j,2) = IniVal
            OCEAN(ng) % ad_rvbar(i,j,1) = IniVal
            OCEAN(ng) % ad_rvbar(i,j,2) = IniVal
            OCEAN(ng) % ad_rzeta(i,j,1) = IniVal
            OCEAN(ng) % ad_rzeta(i,j,2) = IniVal

            OCEAN(ng) % ad_ubar(i,j,1) = IniVal
            OCEAN(ng) % ad_ubar(i,j,2) = IniVal
            OCEAN(ng) % ad_ubar(i,j,3) = IniVal
            OCEAN(ng) % ad_vbar(i,j,1) = IniVal
            OCEAN(ng) % ad_vbar(i,j,2) = IniVal
            OCEAN(ng) % ad_vbar(i,j,3) = IniVal
            OCEAN(ng) % ad_zeta(i,j,1) = IniVal
            OCEAN(ng) % ad_zeta(i,j,2) = IniVal
            OCEAN(ng) % ad_zeta(i,j,3) = IniVal

            OCEAN(ng) % ad_ubar_sol(i,j) = IniVal
            OCEAN(ng) % ad_vbar_sol(i,j) = IniVal
            OCEAN(ng) % ad_zeta_sol(i,j) = IniVal
# if defined NEARSHORE_MELLOR
            OCEAN(ng) % ad_rulag2d(i,j) = IniVal
            OCEAN(ng) % ad_rvlag2d(i,j) = IniVal
            OCEAN(ng) % ad_ubar_stokes(i,j) = IniVal
            OCEAN(ng) % ad_vbar_stokes(i,j) = IniVal
# endif
          END DO
# ifdef SOLVE3D
          DO k=1,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % ad_pden(i,j,k) = IniVal
              OCEAN(ng) % ad_rho(i,j,k) = IniVal

              OCEAN(ng) % ad_u(i,j,k,1) = IniVal
              OCEAN(ng) % ad_u(i,j,k,2) = IniVal
              OCEAN(ng) % ad_v(i,j,k,1) = IniVal
              OCEAN(ng) % ad_v(i,j,k,2) = IniVal
#  if defined NEARSHORE_MELLOR
              OCEAN(ng) % ad_rulag3d(i,j,k) = IniVal
              OCEAN(ng) % ad_rvlag3d(i,j,k) = IniVal
              OCEAN(ng) % ad_u_stokes(i,j,k) = IniVal
              OCEAN(ng) % ad_v_stokes(i,j,k) = IniVal
#  endif
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % ad_ru(i,j,k,1) = IniVal
              OCEAN(ng) % ad_ru(i,j,k,2) = IniVal
              OCEAN(ng) % ad_rv(i,j,k,1) = IniVal
              OCEAN(ng) % ad_rv(i,j,k,2) = IniVal

              OCEAN(ng) % ad_W(i,j,k) = IniVal
            END DO
          END DO
          DO itrc=1,NT(ng)
            DO k=1,N(ng)
              DO i=Imin,Imax
                OCEAN(ng) % ad_t(i,j,k,1,itrc) = IniVal
                OCEAN(ng) % ad_t(i,j,k,2,itrc) = IniVal
                OCEAN(ng) % ad_t(i,j,k,3,itrc) = IniVal
              END DO
            END DO
          END DO
# endif
        END DO
      END IF
#endif

#if defined FOUR_DVAR || defined IMPULSE
!
!  Working arrays to store adjoint impulse forcing, background error
!  covariance, background-error standard deviations, or descent
!  conjugate vectors (directions).
!
      IF (model.eq.0) THEN
        DO j=Jmin,Jmax
          DO rec=1,NSA
            DO i=Imin,Imax
              OCEAN(ng) % b_ubar(i,j,rec) = IniVal
              OCEAN(ng) % b_vbar(i,j,rec) = IniVal
              OCEAN(ng) % b_zeta(i,j,rec) = IniVal
# ifdef FOUR_DVAR
              OCEAN(ng) % e_ubar(i,j,rec) = IniVal
              OCEAN(ng) % e_vbar(i,j,rec) = IniVal
              OCEAN(ng) % e_zeta(i,j,rec) = IniVal
# endif
            END DO
          END DO
# ifdef FOUR_DVAR
          DO i=Imin,Imax
            OCEAN(ng) % d_ubar(i,j) = IniVal
            OCEAN(ng) % d_vbar(i,j) = IniVal
            OCEAN(ng) % d_zeta(i,j) = IniVal

#  ifdef WEAK_CONSTRAINT
            OCEAN(ng) % f_ubar(i,j) = IniVal
            OCEAN(ng) % f_vbar(i,j) = IniVal
            OCEAN(ng) % f_zeta(i,j) = IniVal
#  endif
          END DO
# endif
# ifdef SOLVE3D
          DO rec=1,NSA
            DO k=1,N(ng)
              DO i=Imin,Imax
                OCEAN(ng) % b_u(i,j,k,rec) = IniVal
                OCEAN(ng) % b_v(i,j,k,rec) = IniVal
#  ifdef FOUR_DVAR
                OCEAN(ng) % e_u(i,j,k,rec) = IniVal
                OCEAN(ng) % e_v(i,j,k,rec) = IniVal
#  endif
              END DO
            END DO
          END DO
#  ifdef FOUR_DVAR
          DO k=1,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % d_u(i,j,k) = IniVal
              OCEAN(ng) % d_v(i,j,k) = IniVal

#   ifdef WEAK_CONSTRAINT
              OCEAN(ng) % f_u(i,j,k) = IniVal
              OCEAN(ng) % f_v(i,j,k) = IniVal
#   endif
            END DO
          END DO
#  endif
          DO itrc=1,NT(ng)
            DO rec=1,NSA
              DO k=1,N(ng)
                DO i=Imin,Imax
                  OCEAN(ng) % b_t(i,j,k,rec,itrc) = IniVal
#  ifdef FOUR_DVAR
                  OCEAN(ng) % e_t(i,j,k,rec,itrc) = IniVal
#  endif
                END DO
              END DO
            END DO
#  ifdef FOUR_DVAR
            DO k=1,N(ng)
              DO i=Imin,Imax
                OCEAN(ng) % d_t(i,j,k,itrc) = IniVal

#   ifdef WEAK_CONSTRAINT
                OCEAN(ng) % f_t(i,j,k,itrc) = IniVal
#   endif
              END DO
            END DO
#  endif
          END DO
# endif
        END DO
      END IF
#endif

#if defined FORWARD_READ && \
   (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
      IF (model.eq.0) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
# ifdef FORWARD_RHS
            OCEAN(ng) % rubarG(i,j,1) = IniVal
            OCEAN(ng) % rubarG(i,j,2) = IniVal
            OCEAN(ng) % rvbarG(i,j,1) = IniVal
            OCEAN(ng) % rvbarG(i,j,2) = IniVal
            OCEAN(ng) % rzetaG(i,j,1) = IniVal
            OCEAN(ng) % rzetaG(i,j,2) = IniVal
# endif
            OCEAN(ng) % ubarG(i,j,1) = IniVal
            OCEAN(ng) % ubarG(i,j,2) = IniVal
            OCEAN(ng) % vbarG(i,j,1) = IniVal
            OCEAN(ng) % vbarG(i,j,2) = IniVal
            OCEAN(ng) % zetaG(i,j,1) = IniVal
            OCEAN(ng) % zetaG(i,j,2) = IniVal
# ifdef WEAK_CONSTRAINT
            OCEAN(ng) % f_zetaG(i,j,1) = IniVal
            OCEAN(ng) % f_zetaG(i,j,2) = IniVal
#  ifndef SOLVE3D
            OCEAN(ng) % f_ubarG(i,j,1) = IniVal
            OCEAN(ng) % f_ubarG(i,j,2) = IniVal
            OCEAN(ng) % f_vbarG(i,j,1) = IniVal
            OCEAN(ng) % f_vbarG(i,j,2) = IniVal
#  endif
# endif
          END DO
# ifdef SOLVE3D
          DO k=1,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % uG(i,j,k,1) = IniVal
              OCEAN(ng) % uG(i,j,k,2) = IniVal
              OCEAN(ng) % vG(i,j,k,1) = IniVal
              OCEAN(ng) % vG(i,j,k,2) = IniVal
#  ifdef WEAK_CONSTRAINT
              OCEAN(ng) % f_uG(i,j,k,1) = IniVal
              OCEAN(ng) % f_uG(i,j,k,2) = IniVal
              OCEAN(ng) % f_vG(i,j,k,1) = IniVal
              OCEAN(ng) % f_vG(i,j,k,2) = IniVal
#  endif
            END DO
          END DO
#  ifdef FORWARD_RHS
          DO k=0,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % ruG(i,j,k,1) = IniVal
              OCEAN(ng) % ruG(i,j,k,2) = IniVal
              OCEAN(ng) % rvG(i,j,k,1) = IniVal
              OCEAN(ng) % rvG(i,j,k,2) = IniVal
            END DO
          END DO
#  endif
          DO itrc=1,NT(ng)
            DO k=1,N(ng)
              DO i=Imin,Imax
                OCEAN(ng) % tG(i,j,k,1,itrc) = IniVal
                OCEAN(ng) % tG(i,j,k,2,itrc) = IniVal
#  ifdef WEAK_CONSTRAINT
                OCEAN(ng) % f_tG(i,j,k,1,itrc) = IniVal
                OCEAN(ng) % f_tG(i,j,k,2,itrc) = IniVal
#  endif
              END DO
            END DO
          END DO
# endif
        END DO
      END IF
#endif

      RETURN
      END SUBROUTINE initialize_ocean

      END MODULE mod_ocean
