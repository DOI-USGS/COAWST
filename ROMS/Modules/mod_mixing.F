#include "cppdefs.h"
      MODULE mod_mixing
!
!svn $Id: mod_mixing.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Horizontal and vertical mixing coefficients:                        !
!                                                                      !
!  Akt          Vertical mixing coefficient (m2/s) for tracers.        !
!  Akv          Vertical mixing coefficient (m2/s) for momentum.       !
#ifdef DIFF_3DCOEF
!  Hdiffusion   Time invariant, horizontal diffusion (harmonic or      !
!                 (biharmonic) term at RHO-points.                     !
#endif
#ifdef VISC_3DCOEF
!  Hviscosity   Time invariant, horizontal viscosity (harmonic or      !
!                 (biharmonic) term at RHO-points.                     !
#endif
!  dAktdz       Vertical gradient in mixing coefficient (m/s) for      !
!                 tracer 1, used in float random walk calculations     !
!  diff2        Horizontal, time invariant harmonic coefficient        !
!                 (m2/s) for tracers.                                  !
!  diff4        Horizontal, time invariant biharmonic coefficient      !
#ifdef DIFF_3DCOEF
!                 (m4/s) for tracers.                                  !
#else
!                 SQRT(m4/s) for tracers.                              !
#endif
#ifdef DIFF_3DCOEF
!  diff3d_r     Horizontal, time-dependent 3D diffusion coefficient at !
!                 RHO-points.                                          !
# ifdef TS_U3ADV_SPLIT
!  diff3d_u     Horizontal, time-dependent 3D diffusion coefficient at !
!                 U-points.                                            !
!  diff3d_v     Horizontal, time-dependent 3D diffusiob coefficient at !
!                 V-points.                                            !
# endif
#endif
!  visc2_r      Horizontal, time invariant harmonic viscosity          !
!                 coefficient (m2/s) at RHO-points.                    !
!  visc2_p      Horizontal, time invariant harmonic viscosity          !
!                 coefficient (m2/s) at PSI-points.                    !
!  visc4_r      Horizontal, time invariant harmonic viscosity          !
#ifdef DIFF_3DCOEF
!                 coefficient (m4/s) at RHO-points.                    !
#else
!                 coefficient SQRT(m4/s) at RHO-points.                !
#endif
!  visc4_p      Horizontal, time invariant harmonic viscosity          !
#ifdef DIFF_3DCOEF
!                 coefficient (m4/s) at RHO-points.                    !
#else
!                 coefficient SQRT(m4/s) at RHO-points.                !
#endif
#ifdef DIFF_3DCOEF
!  visc3d_r     Horizontal, time-dependent 3D viscosity coefficient    !
!                 at RHO-points.                                       !
# ifdef UV_U3ADV_SPLIT
!  Uvis3d_r     Horizontal, time-dependent 3D U-viscosity coefficient  !
!                 at RHO-points.                                       !
!  Vvis3d_r     Horizontal, time-dependent 3D V-viscosity coefficient  !
!                 at RHO-points.                                       !
# endif
#endif
!                                                                      !
#ifdef FOUR_DVAR
!  Convolutions diffusion coefficients:                                !
!                                                                      !
!  Kh           Convolution horizontal diffusion coefficient (m2/s).   !
!  Kv           Convolution vertical diffusion coefficient (m2/s).     !
!
#endif
!  Variables associated with the equation of state:                    !
!                                                                      !
!  alpha        Surface thermal expansion coefficient (1/Celsius).     !
!  beta         Surface saline contraction coefficient (1/PSU).        !
!  bvf          Brunt-Vaisala frequency squared (1/s2).                !
!  neutral      Coefficient to convert "in situ" density to neutral    !
!                 surface.                                             !
!                                                                      !
!  tke          Turbulent energy squared (m2/s2) at horizontal         !
!                 at W-points.                                         !
!  gls          Turbulent energy squared times turbulent length        !
!                 scale (m3/s2) at W-points.                           !
!                                                                      !
!  Large/McWilliams/Doney interior vertical mixing variables:          !
!                                                                      !
!  alfaobeta    Ratio of thermal expansion and saline contraction      !
!                 coefficients (Celsius/PSU) used in double            !
!                 diffusion.                                           !
!                                                                      !
!  Water clarity parameters:                                           !
!                                                                      !
!  Jwtype       Water clarity (Jerlov water type classification).      !
!                                                                      !
!  Large/McWilliams/Doney oceanic boundary layer variables:            !
!                                                                      !
!  ghats        Boundary layer nonlocal transport (T units/m).         !
!  hbbl         Depth of bottom oceanic boundary layer (m).            !
!  hsbl         Depth of surface oceanic boundary layer (m).           !
!  kbbl         Index of grid level above bottom  boundary layer.      !
!  ksbl         Index of grid level below surface boundary layer.      !
!                                                                      !
#ifdef NEARSHORE_MELLOR
!  Nearshore radiation stresses.                                       !
!                                                                      !
!  Sxx_bar      2D nearshore radiation stress, xx-component.           !
!  Sxy_bar      2D nearshore radiation stress, xy-component.           !
!  Syy_bar      2D nearshore radiation stress, yy-component.           !
!  Sxx          3D Nearshore horizontal radiation stress, xx-component.!
!  Sxy          3D Nearshore horizontal radiation stress, xx-component.!
!  Syy          3D Nearshore horizontal radiation stress, xx-component.!
!  Szx          3D Nearshore vertical radiation stress, xx-component.  !
!  Szy          3D Nearshore vertical radiation stress, xx-component.  !
!  rustr2d      2D radiation stress tensor in Xi-direction.            !
!  rvstr2d      2D radiation stress tensor in ETA-direction.           !
!  rustr3d      3D radiation stress tensor in Xi-direction.            !
!  rvstr3d      3D radiation stress tensor in ETA-direction.           !
!                                                                      !
#endif
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_MIXING
!
!  Nonlinear model state.
!
#ifdef SOLVE3D
# if defined LMD_SKPP || defined SOLAR_SOURCE
          integer,  pointer :: Jwtype(:,:)
# endif
# if defined LMD_SKPP || defined LMD_BKPP
          integer,  pointer :: ksbl(:,:)
#  ifdef LMD_BKPP
          integer,  pointer :: kbbl(:,:)
#  endif
# endif
#endif
#if defined DIFF_3DCOEF && defined SOLVE3D
          real(r8), pointer :: Hdiffusion(:,:)
#endif
#if defined VISC_3DCOEF && defined SOLVE3D
          real(r8), pointer :: Hviscosity(:,:)
#endif
#if defined UV_VIS2 || !defined SOLVE3D
          real(r8), pointer :: visc2_p(:,:)
          real(r8), pointer :: visc2_r(:,:)
#endif
#ifdef UV_VIS4
          real(r8), pointer :: visc4_p(:,:)
          real(r8), pointer :: visc4_r(:,:)
#endif
# ifdef VISC_3DCOEF
#  ifdef UV_U3ADV_SPLIT
          real(r8), pointer :: Uvis3d_r(:,:,:)
          real(r8), pointer :: Vvis3d_r(:,:,:)
#  else
          real(r8), pointer :: visc3d_r(:,:,:)
#  endif
# endif
#ifdef NEARSHORE_MELLOR
          real(r8), pointer :: Sxx_bar(:,:)
          real(r8), pointer :: Sxy_bar(:,:)
          real(r8), pointer :: Syy_bar(:,:)
          real(r8), pointer :: rustr2d(:,:)
          real(r8), pointer :: rvstr2d(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: Sxx(:,:,:)
          real(r8), pointer :: Sxy(:,:,:)
          real(r8), pointer :: Syy(:,:,:)
          real(r8), pointer :: Szx(:,:,:)
          real(r8), pointer :: Szy(:,:,:)
          real(r8), pointer :: rustr3d(:,:,:)
          real(r8), pointer :: rvstr3d(:,:,:)
# endif
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
          real(r8), pointer :: diff2(:,:,:)
# endif
# ifdef TS_DIF4
          real(r8), pointer :: diff4(:,:,:)
# endif
# ifdef DIFF_3DCOEF
#  ifdef TS_U3ADV_SPLIT
          real(r8), pointer :: diff3d_u(:,:,:)
          real(r8), pointer :: diff3d_v(:,:,:)
#  else
          real(r8), pointer :: diff3d_r(:,:,:)
#  endif
# endif
          real(r8), pointer :: Akv(:,:,:)
          real(r8), pointer :: Akt(:,:,:,:)
# ifdef FLOAT_VWALK
          real(r8), pointer :: dAktdz(:,:,:)
# endif
# if defined LMD_SKPP    || defined LMD_BKPP         || \
     defined BULK_FLUXES || defined BALANCE_OPERATOR
          real(r8), pointer :: alpha(:,:)
          real(r8), pointer :: beta(:,:)
# endif
# ifdef BV_FREQUENCY
          real(r8), pointer :: bvf(:,:,:)
# endif
# if defined MIX_ISO_TS || defined GENT_McWILLIAMS
          real(r8), pointer :: neutral(:,:,:)
# endif
# if defined GLS_MIXING || defined MY25_MIXING
          real(r8), pointer :: tke(:,:,:,:)
          real(r8), pointer :: gls(:,:,:,:)
          real(r8), pointer :: Lscale(:,:,:)
          real(r8), pointer :: Akk(:,:,:)
#  ifdef GLS_MIXING
          real(r8), pointer :: Akp(:,:,:)
#  endif
# endif
# if defined LMD_MIXING && defined LMD_DDMIX
          real(r8), pointer :: alfaobeta(:,:,:)
# endif
# if defined LMD_SKPP || defined LMD_BKPP
          real(r8), pointer :: hsbl(:,:)
#  ifdef LMD_BKPP
          real(r8), pointer :: hbbl(:,:)
#  endif
#  ifdef LMD_NONLOCAL
          real(r8), pointer :: ghats(:,:,:,:)
#  endif
# endif
#endif

#ifdef FOUR_DVAR
!
!  Spatial convolution diffusion coefficients.
!
          real(r8), pointer :: Kh(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: Kv(:,:,:)
# endif
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
# ifdef SOLVE3D
#  ifdef DIFF_3DCOEF
#   ifdef TS_U3ADV_SPLIT
          real(r8), pointer :: tl_diff3d_u(:,:,:)
          real(r8), pointer :: tl_diff3d_v(:,:,:)
#   else
          real(r8), pointer :: tl_diff3d_r(:,:,:)
#   endif
#  endif
#  ifdef VISC_3DCOEF
#   ifdef UV_U3ADV_SPLIT
          real(r8), pointer :: tl_Uvis3d_r(:,:,:)
          real(r8), pointer :: tl_Vvis3d_r(:,:,:)
#   else
          real(r8), pointer :: tl_visc3d_r(:,:,:)
#   endif
#  endif
          real(r8), pointer :: tl_Akv(:,:,:)
          real(r8), pointer :: tl_Akt(:,:,:,:)
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
          real(r8), pointer :: tl_alpha(:,:)
          real(r8), pointer :: tl_beta(:,:)
#  endif
#  ifdef BV_FREQUENCY
          real(r8), pointer :: tl_bvf(:,:,:)
#  endif
#  if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
          real(r8), pointer :: tl_tke(:,:,:,:)
          real(r8), pointer :: tl_gls(:,:,:,:)
          real(r8), pointer :: tl_Lscale(:,:,:)
          real(r8), pointer :: tl_Akk(:,:,:)
#   ifdef GLS_MIXING_NOT_YET
          real(r8), pointer :: tl_Akp(:,:,:)
#   endif
#  endif
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
# ifdef SOLVE3D
#  ifdef DIFF_3DCOEF
#   ifdef TS_U3ADV_SPLIT
          real(r8), pointer :: ad_diff3d_u(:,:,:)
          real(r8), pointer :: ad_diff3d_v(:,:,:)
#   else
          real(r8), pointer :: ad_diff3d_r(:,:,:)
#   endif
#  endif
#  ifdef VISC_3DCOEF
#   ifdef UV_U3ADV_SPLIT
          real(r8), pointer :: ad_Uvis3d_r(:,:,:)
          real(r8), pointer :: ad_Vvis3d_r(:,:,:)
#   else
          real(r8), pointer :: ad_visc3d_r(:,:,:)
#   endif
#  endif
          real(r8), pointer :: ad_Akv(:,:,:)
          real(r8), pointer :: ad_Akt(:,:,:,:)
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
          real(r8), pointer :: ad_alpha(:,:)
          real(r8), pointer :: ad_beta(:,:)
#  endif
#  ifdef BV_FREQUENCY
          real(r8), pointer :: ad_bvf(:,:,:)
#  endif
#  if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
          real(r8), pointer :: ad_tke(:,:,:,:)
          real(r8), pointer :: ad_gls(:,:,:,:)
          real(r8), pointer :: ad_Lscale(:,:,:)
          real(r8), pointer :: ad_Akk(:,:,:)
#   ifdef GLS_MIXING_NOT_YET
          real(r8), pointer :: ad_Akp(:,:,:)
#   endif
#  endif
# endif
#endif

#if defined FORWARD_READ && \
   (defined TANGENT || defined TL_IOMS || defined ADJOINT)
# ifdef FORWARD_MIXING
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
          real(r8), pointer :: AkvG(:,:,:,:)
          real(r8), pointer :: AktG(:,:,:,:,:)
# endif
# if defined LMD_MIXING_NOT_YET
          real(r8), pointer :: hsblG(:,:,:)
# endif
# if defined LMD_BKPP_NOT_YET
          real(r8), pointer :: hbblG(:,:,:)
# endif
# if defined LMD_NONLOCAL_NOT_YET
          real(r8), pointer :: ghatsG(:,:,:,:,:)
# endif
# if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
          real(r8), pointer :: tkeG(:,:,:,:)
          real(r8), pointer :: glsG(:,:,:,:)
          real(r8), pointer :: LscaleG(:,:,:,:)
          real(r8), pointer :: AkkG(:,:,:,:)
#  ifdef GLS_MIXING_NOT_YET
          real(r8), pointer :: AkpG(:,:,:,:)
#  endif
# endif
#endif

        END TYPE T_MIXING

        TYPE (T_MIXING), allocatable :: MIXING(:)

      CONTAINS
        
      SUBROUTINE allocate_mixing (ng, LBi, UBi, LBj, UBj)
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
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( MIXING(Ngrids) )
! 
!  Nonlinear model state.
!
#if defined UV_VIS2 || !defined SOLVE3D
      allocate ( MIXING(ng) % visc2_p(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % visc2_r(LBi:UBi,LBj:UBj) )
#endif

#ifdef UV_VIS4
      allocate ( MIXING(ng) % visc4_p(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % visc4_r(LBi:UBi,LBj:UBj) )
#endif

# ifdef VISC_3DCOEF
      allocate ( MIXING(ng) % Hviscosity(LBi:UBi,LBj:UBj) )
#  ifdef UV_U3ADV_SPLIT
      allocate ( MIXING(ng) % Uvis3d_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( MIXING(ng) % Vvis3d_r(LBi:UBi,LBj:UBj,N(ng)) )
#  else
      allocate ( MIXING(ng) % visc3d_r(LBi:UBi,LBj:UBj,N(ng)) )
#  endif
# endif

#ifdef NEARSHORE_MELLOR
      allocate ( MIXING(ng) % Sxx_bar(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % Sxy_bar(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % Syy_bar(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % rustr2d(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % rvstr2d(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( MIXING(ng) % Sxx(LBi:UBi,LBj:UBj,1:N(ng)) )
      allocate ( MIXING(ng) % Sxy(LBi:UBi,LBj:UBj,1:N(ng)) )
      allocate ( MIXING(ng) % Syy(LBi:UBi,LBj:UBj,1:N(ng)) )
      allocate ( MIXING(ng) % Szx(LBi:UBi,LBj:UBj,1:N(ng)) )
      allocate ( MIXING(ng) % Szy(LBi:UBi,LBj:UBj,1:N(ng)) )
      allocate ( MIXING(ng) % rustr3d(LBi:UBi,LBj:UBj,1:N(ng)) )
      allocate ( MIXING(ng) % rvstr3d(LBi:UBi,LBj:UBj,1:N(ng)) )
# endif
#endif

#ifdef SOLVE3D
# ifdef TS_DIF2
      allocate ( MIXING(ng) % diff2(LBi:UBi,LBj:UBj,NT(ng)) )
# endif

# ifdef TS_DIF4
      allocate ( MIXING(ng) % diff4(LBi:UBi,LBj:UBj,NT(ng)) )
# endif

# ifdef DIFF_3DCOEF
      allocate ( MIXING(ng) % Hdiffusion(LBi:UBi,LBj:UBj) )
#  ifdef TS_U3ADV_SPLIT
      allocate ( MIXING(ng) % diff3d_u(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( MIXING(ng) % diff3d_v(LBi:UBi,LBj:UBj,N(ng)) )
#  else
      allocate ( MIXING(ng) % diff3d_r(LBi:UBi,LBj:UBj,N(ng)) )
#  endif
# endif

      allocate ( MIXING(ng) % Akv(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT) )

# ifdef FLOAT_VWALK
      allocate ( MIXING(ng) % dAktdz(LBi:UBi,LBj:UBj,N(ng)) )
# endif

# if defined LMD_SKPP    || defined LMD_BKPP         || \
     defined BULK_FLUXES || defined BALANCE_OPERATOR
      allocate ( MIXING(ng) % alpha(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % beta(LBi:UBi,LBj:UBj) )
# endif

# ifdef BV_FREQUENCY
      allocate ( MIXING(ng) % bvf(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif

# if defined MIX_ISO_TS || defined GENT_McWILLIAMS
      allocate ( MIXING(ng) % neutral(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif

# if defined GLS_MIXING || defined MY25_MIXING
      allocate ( MIXING(ng) % tke(LBi:UBi,LBj:UBj,0:N(ng),3) )
      allocate ( MIXING(ng) % gls(LBi:UBi,LBj:UBj,0:N(ng),3) )
      allocate ( MIXING(ng) % Lscale(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % Akk(LBi:UBi,LBj:UBj,0:N(ng)) )

#  ifdef GLS_MIXING
      allocate ( MIXING(ng) % Akp(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
# endif

# if defined LMD_MIXING && defined LMD_DDMIX
      allocate ( MIXING(ng) % alfaobeta(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif

# if defined LMD_SKPP || defined SOLAR_SOURCE
      allocate ( MIXING(ng) % Jwtype(LBi:UBi,LBj:UBj) )
# endif

# if defined LMD_SKPP || defined LMD_BKPP
      allocate ( MIXING(ng) % ksbl(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % hsbl(LBi:UBi,LBj:UBj) )

#  ifdef LMD_BKPP
      allocate ( MIXING(ng) % kbbl(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % hbbl(LBi:UBi,LBj:UBj) )
#  endif

#  ifdef LMD_NONLOCAL
      allocate ( MIXING(ng) % ghats(LBi:UBi,LBj:UBj,0:N(ng),NAT) )
#  endif

# endif
#endif

#ifdef FOUR_DVAR
!
!  Spatial convolution diffusion coefficients.
!
      allocate ( MIXING(ng) % Kh(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( MIXING(ng) % Kv(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
# ifdef SOLVE3D
#  ifdef DIFF_3DCOEF
#   ifdef TS_U3ADV_SPLIT
      allocate ( MIXING(ng) % tl_diff3d_u(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( MIXING(ng) % tl_diff3d_v(LBi:UBi,LBj:UBj,N(ng)) )
#   else
      allocate ( MIXING(ng) % tl_diff3d_r(LBi:UBi,LBj:UBj,N(ng)) )
#   endif
#  endif

#  ifdef VISC_3DCOEF
#   ifdef UV_U3ADV_SPLIT
      allocate ( MIXING(ng) % tl_Uvis3d_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( MIXING(ng) % tl_Vvis3d_r(LBi:UBi,LBj:UBj,N(ng)) )
#   else
      allocate ( MIXING(ng) % tl_visc3d_r(LBi:UBi,LBj:UBj,N(ng)) )
#   endif
#  endif

      allocate ( MIXING(ng) % tl_Akv(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % tl_Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT) )

#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
      allocate ( MIXING(ng) % tl_alpha(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % tl_beta(LBi:UBi,LBj:UBj) )
#  endif

#  ifdef BV_FREQUENCY
      allocate ( MIXING(ng) % tl_bvf(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif

#  if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
      allocate ( MIXING(ng) % tl_tke(LBi:UBi,LBj:UBj,0:N(ng),3) )
      allocate ( MIXING(ng) % tl_gls(LBi:UBi,LBj:UBj,0:N(ng),3) )
      allocate ( MIXING(ng) % tl_Lscale(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % tl_Akk(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef GLS_MIXING_NOT_YET
      allocate ( MIXING(ng) % tl_Akp(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
# ifdef SOLVE3D
#  ifdef DIFF_3DCOEF
#   ifdef TS_U3ADV_SPLIT
      allocate ( MIXING(ng) % ad_diff3d_u(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( MIXING(ng) % ad_diff3d_v(LBi:UBi,LBj:UBj,N(ng)) )
#   else
      allocate ( MIXING(ng) % ad_diff3d_r(LBi:UBi,LBj:UBj,N(ng)) )
#   endif
#  endif

#  ifdef VISC_3DCOEF
#   ifdef UV_U3ADV_SPLIT
      allocate ( MIXING(ng) % ad_Uvis3d_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( MIXING(ng) % ad_Vvis3d_r(LBi:UBi,LBj:UBj,N(ng)) )
#   else
      allocate ( MIXING(ng) % ad_visc3d_r(LBi:UBi,LBj:UBj,N(ng)) )
#   endif
#  endif

      allocate ( MIXING(ng) % ad_Akv(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % ad_Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT) )

#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
      allocate ( MIXING(ng) % ad_alpha(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % ad_beta(LBi:UBi,LBj:UBj) )
#  endif

#  ifdef BV_FREQUENCY
      allocate ( MIXING(ng) % ad_bvf(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif

#  if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
      allocate ( MIXING(ng) % ad_tke(LBi:UBi,LBj:UBj,0:N(ng),3) )
      allocate ( MIXING(ng) % ad_gls(LBi:UBi,LBj:UBj,0:N(ng),3) )
      allocate ( MIXING(ng) % ad_Lscale(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % ad_Akk(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
#  ifdef GLS_MIXING_NOT_YET
      allocate ( MIXING(ng) % ad_Akp(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
# endif
#endif

#if defined FORWARD_READ && \
   (defined TANGENT || defined TL_IOMS || defined ADJOINT)
# ifdef FORWARD_MIXING
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
      allocate ( MIXING(ng) % AkvG(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( MIXING(ng) % AktG(LBi:UBi,LBj:UBj,0:N(ng),2,NAT) )
# endif

# if defined LMD_MIXING_NOT_YET
      allocate ( MIXING(ng) % hsblG(LBi:UBi,LBj:UBj,2) )
# endif

# if defined LMD_BKPP_NOT_YET
      allocate ( MIXING(ng) % hbblG(LBi:UBi,LBj:UBj,2) )
# endif

# if defined LMD_NONLOCAL_NOT_YET
      allocate ( MIXING(ng) % ghatsG(LBi:UBi,LBj:UBj,0:N(ng),2,NAT) )
# endif

# if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
      allocate ( MIXING(ng) % tkeG(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( MIXING(ng) % glsG(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( MIXING(ng) % LscaleG(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( MIXING(ng) % AkkG(LBi:UBi,LBj:UBj,0:N(ng),2) )
#  ifdef GLS_MIXING_NOT_YET
      allocate ( MIXING(ng) % AkpG(LBi:UBi,LBj:UBj,0:N(ng),2) )
#  endif
# endif
#endif

      RETURN
      END SUBROUTINE allocate_mixing

      SUBROUTINE initialize_mixing (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes all variables in module      !
!  "mod_mixing" for all nested grids.                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
#ifdef SOLVE3D
      integer :: itrc, k
#endif

      real(r8), parameter :: IniVal = 0.0_r8
      real(r8) :: cff1, cff2, cff3, cff4
!
#include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=LBi
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=UBi
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=LBj
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=UBj
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
#if defined UV_VIS2 || !defined SOLVE3D
          DO i=Imin,Imax
            MIXING(ng) % visc2_p(i,j) = visc2(ng)
            MIXING(ng) % visc2_r(i,j) = visc2(ng)
          END DO
#endif
#ifdef UV_VIS4
          DO i=Imin,Imax
            MIXING(ng) % visc4_p(i,j) = visc4(ng)
            MIXING(ng) % visc4_r(i,j) = visc4(ng)
          END DO
#endif
#ifdef SOLVE3D
# ifdef VISC_3DCOEF
          DO k=1,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % Hviscosity(i,j) = IniVal
#  ifdef UV_U3ADV_SPLIT
              MIXING(ng) % Uvis3d_r(i,j,k) = IniVal
              MIXING(ng) % Vvis3d_r(i,j,k) = IniVal
#  else
              MIXING(ng) % visc3d_r(i,j,k) = IniVal
#  endif
            END DO
          END DO
# endif
#endif
#ifdef NEARSHORE_MELLOR
          DO i=Imin,Imax
            MIXING(ng) % Sxx_bar(i,j) = IniVal
            MIXING(ng) % Sxy_bar(i,j) = IniVal
            MIXING(ng) % Syy_bar(i,j) = IniVal
            MIXING(ng) % rustr2d(i,j) = IniVal
            MIXING(ng) % rvstr2d(i,j) = IniVal
# ifdef SOLVE3D
            DO k=1,N(ng)
              MIXING(ng) % Sxx(i,j,k) = IniVal
              MIXING(ng) % Sxy(i,j,k) = IniVal
              MIXING(ng) % Syy(i,j,k) = IniVal
              MIXING(ng) % Szx(i,j,k) = IniVal
              MIXING(ng) % Szy(i,j,k) = IniVal
              MIXING(ng) % rustr3d(i,j,k) = IniVal
              MIXING(ng) % rvstr3d(i,j,k) = IniVal
            END DO
#endif
          END DO
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
          DO itrc=1,NT(ng)
            DO i=Imin,Imax
              MIXING(ng) % diff2(i,j,itrc) = tnu2(itrc,ng)
            END DO
          END DO
# endif
# ifdef TS_DIF4
          DO itrc=1,NT(ng)
            DO i=Imin,Imax
              MIXING(ng) % diff4(i,j,itrc) = tnu4(itrc,ng)
            END DO
          END DO
# endif
# ifdef DIFF_3DCOEF
          DO k=1,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % Hdiffusion(i,j) = IniVal
#  ifdef TS_U3ADV_SPLIT
              MIXING(ng) % diff3d_u(i,j,k) = IniVal
              MIXING(ng) % diff3d_v(i,j,k) = IniVal
#  else
              MIXING(ng) % diff3d_r(i,j,k) = IniVal
#  endif
            END DO
          END DO
# endif
          DO i=Imin,Imax
            MIXING(ng) % Akv(i,j,0) = IniVal
            MIXING(ng) % Akv(i,j,N(ng)) = IniVal
          END DO
          DO k=1,N(ng)-1
            DO i=Imin,Imax
              MIXING(ng) % Akv(i,j,k) = Akv_bak(ng)
            END DO
          END DO
          DO itrc=1,NAT
            DO i=Imin,Imax
              MIXING(ng) % Akt(i,j,0,itrc) = IniVal
              MIXING(ng) % Akt(i,j,N(ng),itrc) = IniVal
            END DO
            DO k=1,N(ng)-1
              DO i=Imin,Imax
                MIXING(ng) % Akt(i,j,k,itrc) = Akt_bak(itrc,ng)
              END DO
            END DO
          END DO
# ifdef FLOAT_VWALK
          DO k=1,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % dAktdz(i,j,k) = IniVal
            END DO
          END DO
# endif
# if defined LMD_SKPP    || defined LMD_BKPP         || \
     defined BULK_FLUXES || defined BALANCE_OPERATOR
          DO i=Imin,Imax
            MIXING(ng) % alpha(i,j) = IniVal
            MIXING(ng) % beta(i,j) = IniVal
          END DO
# endif
# ifdef BV_FREQUENCY
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % bvf(i,j,k) = IniVal
            END DO
          END DO
# endif
# if defined MIX_ISO_TS || defined GENT_McWILLIAMS
          DO k=1,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % neutral(i,j,k) = IniVal
            END DO
          END DO
# endif
# if defined GLS_MIXING || defined MY25_MIXING
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % tke(i,j,k,1) = gls_Kmin(ng)
              MIXING(ng) % tke(i,j,k,2) = gls_Kmin(ng)
              MIXING(ng) % tke(i,j,k,3) = gls_Kmin(ng)
              MIXING(ng) % gls(i,j,k,1) = gls_Pmin(ng)
              MIXING(ng) % gls(i,j,k,2) = gls_Pmin(ng)
              MIXING(ng) % gls(i,j,k,3) = gls_Pmin(ng)
              MIXING(ng) % Lscale(i,j,k) = IniVal
            END DO
          END DO
          DO i=Imin,Imax
            MIXING(ng) % Akk(i,j,0) = IniVal
            MIXING(ng) % Akk(i,j,N(ng)) = IniVal
#  ifdef GLS_MIXING
            MIXING(ng) % Akp(i,j,0) = IniVal
            MIXING(ng) % Akp(i,j,N(ng)) = IniVal
#  endif
          END DO
          DO k=1,N(ng)-1
            DO i=Imin,Imax
              MIXING(ng) % Akk(i,j,k) = Akk_bak(ng)
#  ifdef GLS_MIXING
              MIXING(ng) % Akp(i,j,k) = Akp_bak(ng)
#  endif
            END DO
          END DO
# endif
# if defined LMD_MIXING && defined LMD_DDMIX
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % alfaobeta(i,j,k) = IniVal
            END DO
          END DO
# endif
# if defined LMD_SKPP || defined SOLAR_SOURCE
          DO i=Imin,Imax
            MIXING(ng) % Jwtype(i,j) = lmd_Jwt(ng)
          END DO
# endif
# if defined LMD_SKPP || defined LMD_BKPP
          DO i=Imin,Imax
            MIXING(ng) % ksbl(i,j) = 0
            MIXING(ng) % hsbl(i,j) = IniVal
          END DO
#  ifdef LMD_BKPP
          DO i=Imin,Imax
            MIXING(ng) % kbbl(i,j) = 0
            MIXING(ng) % hbbl(i,j) = IniVal
          END DO
#  endif
#  ifdef LMD_NONLOCAL
          DO itrc=1,NAT
            DO k=0,N(ng)
              DO i=Imin,Imax
                MIXING(ng) % ghats(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
#  endif
# endif
#endif
        END DO

#ifdef FOUR_DVAR
!
!  Spatial convolution diffusion coefficients.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
             MIXING(ng) % Kh(i,j) = 1.0_r8
          END DO
# ifdef SOLVE3D
          DO k=0,N(ng)           
            DO i=Imin,Imax
              MIXING(ng) % Kv(i,j,k) = 1.0_r8
            END DO
          END DO
# endif
        END DO
#endif
      END IF

#if defined TANGENT || defined TL_IOMS
# ifdef SOLVE3D
!
!  Tangent linear model state.
!
#  if defined GLS_MIXING || defined MY25_MIXING
      IF (model.eq.iRPM) THEN
        cff1=gls_Kmin(ng)
        cff2=gls_Pmin(ng)
        cff3=Akk_bak(ng)
        cff4=Akp_bak(ng)
      ELSE
        cff1=IniVal
        cff2=IniVal
        cff3=IniVal
        cff4=IniVal
      END IF
#  endif
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN
        DO j=Jmin,Jmax
#  ifdef DIFF_3DCOEF
          DO k=1,N(ng)
            DO i=Imin,Imax
#   ifdef TS_U3ADV_SPLIT
              MIXING(ng) % tl_diff3d_u(i,j,k) = IniVal
              MIXING(ng) % tl_diff3d_v(i,j,k) = IniVal
#   else
              MIXING(ng) % tl_diff3d_r(i,j,k) = IniVal
#   endif
            END DO
          END DO
#  endif
#  ifdef VISC_3DCOEF
          DO k=1,N(ng)
            DO i=Imin,Imax
#   ifdef UV_U3ADV_SPLIT
              MIXING(ng) % tl_Uvis3d_r(i,j,k) = IniVal
              MIXING(ng) % tl_Vvis3d_r(i,j,k) = IniVal
#   else
              MIXING(ng) % tl_visc3d_r(i,j,k) = IniVal
#   endif
            END DO
          END DO
#  endif
          DO k=0,N(ng)           
            DO i=Imin,Imax
              MIXING(ng) % tl_Akv(i,j,k) = IniVal
            END DO
          END DO
          DO itrc=1,NAT
            DO k=0,N(ng)
              DO i=Imin,Imax
                MIXING(ng) % tl_Akt(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
          DO i=Imin,Imax
            MIXING(ng) % tl_alpha(i,j) = IniVal
            MIXING(ng) % tl_beta(i,j) = IniVal
          END DO
#  endif
#  ifdef BV_FREQUENCY
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % tl_bvf(i,j,k) = IniVal
            END DO
          END DO
#  endif
#  if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % tl_tke(i,j,k,1) = cff1
              MIXING(ng) % tl_tke(i,j,k,2) = cff1
              MIXING(ng) % tl_tke(i,j,k,3) = cff1
              MIXING(ng) % tl_gls(i,j,k,1) = cff2
              MIXING(ng) % tl_gls(i,j,k,2) = cff2
              MIXING(ng) % tl_gls(i,j,k,3) = cff2
              MIXING(ng) % tl_Lscale(i,j,k) = IniVal
            END DO
          END DO
          DO i=Imin,Imax
            MIXING(ng) % tl_Akk(i,j,0) = IniVal
            MIXING(ng) % tl_Akk(i,j,N(ng)) = IniVal
#   ifdef GLS_MIXING_NOT_YET
            MIXING(ng) % tl_Akp(i,j,0) = IniVal
            MIXING(ng) % tl_Akp(i,j,N(ng)) = IniVal
#   endif
          END DO
          DO k=1,N(ng)-1
            DO i=Imin,Imax
              MIXING(ng) % tl_Akk(i,j,k) = cff3
#   ifdef GLS_MIXING
              MIXING(ng) % tl_Akp(i,j,k) = cff4
#   endif
            END DO
          END DO
#  endif
        END DO
      END IF
# endif
#endif

#ifdef ADJOINT
# ifdef SOLVE3D
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN
        DO j=Jmin,Jmax
#  ifdef DIFF_3DCOEF
          DO k=1,N(ng)
            DO i=Imin,Imax
#   ifdef TS_U3ADV_SPLIT
              MIXING(ng) % ad_diff3d_u(i,j,k) = IniVal
              MIXING(ng) % ad_diff3d_v(i,j,k) = IniVal
#   else
              MIXING(ng) % ad_diff3d_r(i,j,k) = IniVal
#   endif
            END DO
          END DO
#  endif
#  ifdef VISC_3DCOEF
          DO k=1,N(ng)
            DO i=Imin,Imax
#   ifdef UV_U3ADV_SPLIT
              MIXING(ng) % ad_Uvis3d_r(i,j,k) = IniVal
              MIXING(ng) % ad_Vvis3d_r(i,j,k) = IniVal
#   else
              MIXING(ng) % ad_visc3d_r(i,j,k) = IniVal
#   endif
            END DO
          END DO
#  endif
          DO k=0,N(ng)           
            DO i=Imin,Imax
              MIXING(ng) % ad_Akv(i,j,k) = IniVal
            END DO
          END DO
          DO itrc=1,NAT
            DO k=0,N(ng)
              DO i=Imin,Imax
                MIXING(ng) % ad_Akt(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
#  if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
          DO i=Imin,Imax
            MIXING(ng) % ad_alpha(i,j) = IniVal
            MIXING(ng) % ad_beta(i,j) = IniVal
          END DO
#  endif
#  ifdef BV_FREQUENCY
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % ad_bvf(i,j,k) = IniVal
            END DO
          END DO
#  endif
#  if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % ad_tke(i,j,k,1) = IniVal
              MIXING(ng) % ad_tke(i,j,k,2) = IniVal
              MIXING(ng) % ad_tke(i,j,k,3) = IniVal
              MIXING(ng) % ad_gls(i,j,k,1) = IniVal
              MIXING(ng) % ad_gls(i,j,k,2) = IniVal
              MIXING(ng) % ad_gls(i,j,k,3) = IniVal
              MIXING(ng) % ad_Lscale(i,j,k) = IniVal
              MIXING(ng) % ad_Akk(i,j,k) = IniVal
#   ifdef GLS_MIXING
              MIXING(ng) % ad_Akp(i,j,k) = IniVal
#   endif
            END DO
          END DO
#  endif
        END DO
      END IF
# endif
#endif

#if defined FORWARD_READ && defined FOWARD_MIXING && defined SOLVE3D && \
   (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
      IF (model.eq.0) THEN
        DO j=Jmin,Jmax
# ifdef FORWARD_MIXING
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % AkvG(i,j,k,1) = IniVal
              MIXING(ng) % AkvG(i,j,k,2) = IniVal
            END DO
          END DO
          DO itrc=1,NAT
            DO k=0,N(ng)
              DO i=Imin,Imax
                MIXING(ng) % AktG(i,j,k,1,itrc) = IniVal
                MIXING(ng) % AktG(i,j,k,2,itrc) = IniVal
              END DO
            END DO
          END DO
# endif
# if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
          DO k=0,N(ng)
            DO i=Imin,Imax
              MIXING(ng) % tkeG(i,j,k,1) = IniVal
              MIXING(ng) % tkeG(i,j,k,2) = IniVal
              MIXING(ng) % glsG(i,j,k,1) = IniVal
              MIXING(ng) % glsG(i,j,k,2) = IniVal
              MIXING(ng) % LscaleG(i,j,k,1) = IniVal
              MIXING(ng) % LscaleG(i,j,k,2) = IniVal
              MIXING(ng) % AkkG(i,j,k,1) = IniVal
              MIXING(ng) % AkkG(i,j,k,2) = IniVal
#  ifdef GLS_MIXING_NOT_YET
              MIXING(ng) % AkpG(i,j,k,1) = IniVal
              MIXING(ng) % AkpG(i,j,k,2) = IniVal
#  endif
            END DO
          END DO
# endif
# if defined LMD_MIXING_NOT_YET
          DO i=Imin,Imax
            MIXING(ng) % hsblG(i,j,1) = IniVal
            MIXING(ng) % hsblG(i,j,2) = IniVal
          END DO
# endif
# if defined LMD_BKPP_NOT_YET
          DO i=Imin,Imax
            MIXING(ng) % hbblG(i,j,1) = IniVal
            MIXING(ng) % hbblG(i,j,2) = IniVal
          END DO
# endif
# if defined LMD_NONLOCAL_NOT_YET
          DO itrc=1,NAT
            DO k=0,N(ng)
              DO i=Imin,Imax
                MIXING(ng) % ghatsG(i,j,0:N(ng),1,itrc) = IniVal
                MIXING(ng) % ghatsG(i,j,0:N(ng),2,itrc) = IniVal
              END DO
            END DO
          END DO
# endif
        END DO
      END IF
#endif

      RETURN
      END SUBROUTINE initialize_mixing

      END MODULE mod_mixing
