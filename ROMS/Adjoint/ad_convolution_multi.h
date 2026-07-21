#include "cppdefs.h"
      MODULE ad_convolution_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2024 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Multi-scale, Background-Error Covariace modeling:                   !
!                                                                      !
!  This module simulates the propagation effects of the background-    !
!  error covariance matrix (B) for each variable in the data           !
!  assimilation control vector by convolving generalized pseudo-       !
!  diffusion tangent-linear operators. In data assimilation, the       !
!  observational information is spatially distributed according to     !
!  the specified correlation (C) functions defined by the application  !
!  error hypothesis.                                                   !
!                                                                      !
!  The background error covariance is defined as:                      !
!                                                                      !
!    B = S C S                                                         !
!                                                                      !
!    C = C^(1/2) C^(T/2)                                               !
!                                                                      !
!    C^(1/2) = G L^(1/2) W^(-1/2)                               TLM    !
!    C^(T/2) = W^(-1/2) L^(T/2) G                               ADM    !
!                                                                      !
!  where                                                               !
!                                                                      !
!    B : background-error covariance matrix                            !
!    S : diagonal matrix of background-error standard deviations       !
!    C : symmetric matrix of background-error correlations             !
!    G : normalization coefficients matrix used to ensure that the     !
!          diagonal variances of C are equal to unity.                 !
!    L : tangent linear and adjoint diffusion operators                !
!    W : diagonal matrix of local area or volume metrics used to       !
!          convert L into a symmetric matrix: LW^(-1).                 !
!                                                                      !
!  Here, T/2 denote the transpose (adjoint) of the squared-root        !
!  operator.                                                           !
!                                                                      !
!  Spreading is modeled using pseudo-diffusion operators in spatial    !
!  correlation space, where diffusion coefficients (K) are proportional!
!  to the square of the correlation length scale (Daley, 1992).        !
!  Correlation scales may be constant or spatially varying for each    !
!  variable in the control vector.                                     !
!                                                                      !
!  In the multiscale approach, B is expressed as a linear combination  !
!  of distinct spatial scales, ranging from large to small. This       !
!  formulation enables the representation of both broad structures     !
!  and fine-scale features, while reducing scale aliasing in the data  !
!  assimilation cost function (Weaver et al., 2016).                   !
!                                                                      !
!  The multiscale approach uses an implicit horizontal pseudo-         !
!  diffusion operator implemented via Conjugate Gradient (CG) and      !
!  Chebyshev Iterations (CI).                                          !
!                                                                      !
!  Error correlations are considered separable in the horizontal and   !
!  vertical directions. The vertical diffusion operator is also        !
!  implicit.                                                           !
!                                                                      !
!  The multiscale concept, expressed as B = SUM(Wi*Bi), represents a   !
!  weighted sum of B values corresponding to various scales. Typically,!
!  two to few different scales are combined in practical applications. !
!  The weight coefficients Wi are constrained to sum to unity. The     !
!  resulting correlation functions belong to the Matern-class family,  !
!  which accommodates complex functional shapes.                       !
!                                                                      !
!  Horizontally spatially varying correlation length scales for the    !
!  K-diffusion coefficient are permitted. A horizontal map of isotropic!
!  or anisotropic correlations can be specified and is read from an    !
!  input NetCDF file. For computational efficiency, these correlations !
!  require an implicit diffusion CG/CI solver.                         !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Weaver, A. and P. Courtier, 2001: Correlation modeling on the     !
!      sphere using a generalized diffusion equation, Q.J.R. Meteorol. !
!      Soc, 127, 1815-1846, doi:10.1002/qj.49712757518.                !
!                                                                      !
!    Weaver, A.T., J. Tshimanga, and A. Piacentini, 2016: Correlation  !
!      operators based on an implicitly formulated diffusion equation  !
!      solved with the Chebyshev iteration, Q.J.R. Meteorol. Soc.,     !
!      142, 455-471, doi:10.1002/qj.2664.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
#ifdef ADJUST_BOUNDARY
      USE mod_boundary
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
#endif
      USE mod_fourdvar
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
#endif
!
      USE ad_bc_2d_mod,        ONLY : ad_dabc_r2d_tile,                 &
     &                                ad_dabc_u2d_tile,                 &
     &                                ad_dabc_v2d_tile
#ifdef SOLVE3D
      USE ad_bc_3d_mod,        ONLY : ad_dabc_r3d_tile,                 &
     &                                ad_dabc_u3d_tile,                 &
     &                                ad_dabc_v3d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod,     ONLY : ad_mp_exchange2d,                 &
     &                                ad_mp_exchange3d,                 &
     &                                ad_mp_exchange4d
# ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,     ONLY : ad_mp_exchange2d_bry,             &
     &                                ad_mp_exchange3d_bry
# endif
#endif
#ifdef SOLVE3D
      USE set_depth_mod,       ONLY : set_depth_tile
#endif
!
      USE roms_multiscale_mod, ONLY : multiscale         ! CLASS object
!
      implicit none
!
      PUBLIC  :: ad_convolution
      PRIVATE :: ad_convolution_tile
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_convolution (ng, tile, ns, Linp, Lweak, ifac)
!***********************************************************************
!
      USE mod_stepping,        ONLY : nnew, nstp
      USE roms_multiscale_mod, ONLY : MSB
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak
!
      integer, intent(in) :: ng, tile, Linp, ifac, ns
!
!  Local variable declarations.
!  (Note: Error covariance normalization coefficients are in b_arrays)
!
#include "tile.h"
!
      CALL ad_convolution_tile (MSB(ng), ng, tile, iADM,                &
     &                          ns, Lweak, ifac,                        &
     &                          nstp(ng), nnew(ng), Linp,               &
     &                          LBi, UBi, LBj, UBj,                     &
#ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
#endif
     &                          IminS, ImaxS, JminS, JmaxS,             &
#ifdef ADJUST_BOUNDARY
     &                          BOUNDARY(ng) % b_ubar_obc,              &
     &                          BOUNDARY(ng) % b_vbar_obc,              &
     &                          BOUNDARY(ng) % b_zeta_obc,              &
# ifdef SOLVE3D
     &                          BOUNDARY(ng) % b_t_obc,                 &
     &                          BOUNDARY(ng) % b_u_obc,                 &
     &                          BOUNDARY(ng) % b_v_obc,                 &
# endif
#endif
#ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % b_sustr,                   &
     &                          FORCES(ng) % b_svstr,                   &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % b_stflx,                   &
#endif
     &                          OCEAN(ng) % b_zeta,                     &
     &                          OCEAN(ng) % b_ubar,                     &
     &                          OCEAN(ng) % b_vbar,                     &
#ifdef SOLVE3D
     &                          OCEAN(ng) % b_t,                        &
     &                          OCEAN(ng) % b_u,                        &
     &                          OCEAN(ng) % b_v,                        &
#endif
#ifdef ADJUST_BOUNDARY
     &                          BOUNDARY(ng) % ad_ubar_obc,             &
     &                          BOUNDARY(ng) % ad_vbar_obc,             &
     &                          BOUNDARY(ng) % ad_zeta_obc,             &
# ifdef SOLVE3D
     &                          BOUNDARY(ng) % ad_t_obc,                &
     &                          BOUNDARY(ng) % ad_u_obc,                &
     &                          BOUNDARY(ng) % ad_v_obc,                &
# endif
#endif
#ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % ad_ustr,                   &
     &                          FORCES(ng) % ad_vstr,                   &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % ad_tflux,                  &
#endif
#ifdef SOLVE3D
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
#endif
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
     &                          OCEAN(ng) % ad_zeta)
!
      RETURN
      END SUBROUTINE ad_convolution
!
!***********************************************************************
      SUBROUTINE ad_convolution_tile (self, ng, tile, model,            &
     &                                ns, Lweak, ifac,                  &
     &                                nstp, nnew, Linp,                 &
     &                                LBi, UBi, LBj, UBj,               &
#ifdef ADJUST_BOUNDARY
     &                                LBij, UBij,                       &
#endif
     &                                IminS, ImaxS, JminS, JmaxS,       &
#ifdef ADJUST_BOUNDARY
     &                                HnormRobc, HnormUobc, HnormVobc,  &
# ifdef SOLVE3D
     &                                VnormRobc, VnormUobc, VnormVobc,  &
# endif
#endif
#ifdef ADJUST_WSTRESS
     &                                HnormSUS, HnormSVS,               &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                                HnormSTF,                         &
#endif
     &                                HnormR, HnormU, HnormV,           &
#ifdef SOLVE3D
     &                                VnormR, VnormU, VnormV,           &
#endif
#ifdef ADJUST_BOUNDARY
     &                                ad_ubar_obc, ad_vbar_obc,         &
     &                                ad_zeta_obc,                      &
# ifdef SOLVE3D
     &                                ad_t_obc, ad_u_obc, ad_v_obc,     &
# endif
#endif
#ifdef ADJUST_WSTRESS
     &                                ad_ustr, ad_vstr,                 &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                                ad_tflux,                         &
#endif
#ifdef SOLVE3D
     &                                ad_t, ad_u, ad_v,                 &
#endif
     &                                ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
!
!  Imported variable declarations.
!
      CLASS (multiscale), intent(inout) :: self
!
      logical, intent(in) :: Lweak
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: ns, ifac
      integer, intent(in) :: nstp, nnew, Linp
      integer, intent(in) :: LBi, UBi, LBj, UBj
#ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#endif
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
# ifdef ADJUST_BOUNDARY
      real(r8), intent (in) :: HnormRobc(LBij:,:)
      real(r8), intent (in) :: HnormUobc(LBij:,:)
      real(r8), intent (in) :: HnormVobc(LBij:,:)
#  ifdef SOLVE3D
      real(r8), intent (in) :: VnormRobc(LBij:,:,:,:)
      real(r8), intent (in) :: VnormUobc(LBij:,:,:)
      real(r8), intent (in) :: VnormVobc(LBij:,:,:)
#  endif
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(in) :: HnormSVS(LBi:,LBj:)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: HnormSTF(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
      real(r8), intent(in) :: VnormR(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: VnormU(LBi:,LBj:,:,:)
      real(r8), intent(in) :: VnormV(LBi:,LBj:,:,:)
# endif
      real(r8), intent(in) :: HnormR(LBi:,LBj:,:)
      real(r8), intent(in) :: HnormU(LBi:,LBj:,:)
      real(r8), intent(in) :: HnormV(LBi:,LBj:,:)
# ifdef ADJUST_BOUNDARY
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#  endif
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
# endif

#else

# ifdef ADJUST_BOUNDARY
      real(r8), intent (in) :: HnormRobc(LBij:UBij,4)
      real(r8), intent (in) :: HnormUobc(LBij:UBij,4)
      real(r8), intent (in) :: HnormVobc(LBij:UBij,4)
#  ifdef SOLVE3D
      real(r8), intent (in) :: VnormRobc(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent (in) :: VnormUobc(LBij:UBij,N(ng),4)
      real(r8), intent (in) :: VnormVobc(LBij:UBij,N(ng),4)
#  endif
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: HnormSVS(LBi:UBi,LBj:UBj)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
# endif
      real(r8), intent(in) :: HnormR(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: HnormU(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: HnormV(LBi:UBi,LBj:UBj,NSA)
# ifdef SOLVE3D
      real(r8), intent(in) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(in) :: VnormU(LBi:UBi,LBj:UBj,NSA,N(ng))
      real(r8), intent(in) :: VnormV(LBi:UBi,LBj:UBj,NSA,N(ng))
# endif
# ifdef ADJUST_BOUNDARY
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
# endif
# ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,:)
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
#endif
!
!  Local variable declarations.
!
#ifdef ADJUST_BOUNDARY
      logical, dimension(4) :: Lconvolve
!
#endif
      integer :: i, ib, ir, j, k, rec
#ifdef SOLVE3D
      integer :: ifield, itrc
#endif
!
      real(r8) :: cff, fac
#ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj) :: work
#endif
!
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine error covariance normalization factors to use.
!-----------------------------------------------------------------------
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint, Model error
      ELSE
        rec=1                        ! strong constraint, Prior error
      END IF

#ifdef ADJUST_BOUNDARY
!
!  Set switch to convolve boundary segments by the appropriate
!  tiles.
!
      Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time invariant depths (use zero free-surface).
!-----------------------------------------------------------------------
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          work(i,j)=0.0_r8                  ! free surface
        END DO
      END DO
!
      CALL set_depth_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     GRID(ng) % h,                                &
# ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     SEDBED(ng) % bed_thick,                      &
# endif
     &                     work,                                        &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Adjoint of multiply data assimilation state vector by its
!  corresponding normalization factor.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Adjoint free-surface.
!
#ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,Linp))
!
#endif
      CALL ad_dabc_r2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_zeta(:,:,Linp))
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)*HnormR(i,j,rec)
#ifdef MASKING
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)*GRID(ng)%rmask(i,j)
#endif
        END DO
      END DO

#ifndef SOLVE3D
!
!  Adjoint 2D momentum.
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ubar(:,:,Linp),                         &
     &                       ad_vbar(:,:,Linp))
!
# endif
      CALL ad_dabc_u2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_ubar(:,:,Linp))
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)*HnormU(i,j,rec)
# ifdef MASKING
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)*GRID(ng)%umask(i,j)
# endif
        END DO
      END DO
!
      CALL ad_dabc_v2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_vbar(:,:,Linp))
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          ad_vbar(i,j,Linp)=ad_vbar(i,j,Linp)*HnormV(i,j,rec)
# ifdef MASKING
          ad_vbar(i,j,Linp)=ad_vbar(i,j,Linp)*GRID(ng)%vmask(i,j)
# endif
        END DO
      END DO
#endif

#ifdef SOLVE3D
!
!  Adjoint 3D momentum.
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_u(:,:,:,Linp),                          &
     &                       ad_v(:,:,:,Linp))
!
# endif
!
      CALL ad_dabc_u3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       ad_u(:,:,:,Linp))
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)*VnormU(i,j,k,rec)
# ifdef MASKING
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)*GRID(ng)%umask(i,j)
# endif
          END DO
        END DO
      END DO
!
      CALL ad_dabc_v3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       ad_v(:,:,:,Linp))
!
      DO k=1,N(ng)
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)*VnormV(i,j,k,rec)
# ifdef MASKING
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)*GRID(ng)%vmask(i,j)
# endif
          END DO
        END DO
      END DO
!
!  Adjoint tracer variables.
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_t(:,:,:,Linp,:))
!
# endif
      DO itrc=1,NT(ng)
        CALL ad_dabc_r3d_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         ad_t(:,:,:,Linp,itrc))
!
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              ad_t(i,j,k,Linp,itrc)=ad_t(i,j,k,Linp,itrc)*              &
     &                              VnormR(i,j,k,rec,itrc)
# ifdef MASKING
              ad_t(i,j,k,Linp,itrc)=ad_t(i,j,k,Linp,itrc)*              &
     &                              GRID(ng)%rmask(i,j)
# endif
            END DO
          END DO
        END DO
      END DO
#endif

#ifdef ADJUST_BOUNDARY
!
!  Adjoint free-surface open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
# ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_zeta_obc(:,ib,ir,Linp))
!
# endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO j=JstrT,JendT
                    ad_zeta_obc(j,ib,ir,Linp)=HnormRobc(j,ib)*          &
     &                                        ad_zeta_obc(j,ib,ir,Linp)
                  END DO
                CASE (isouth, inorth)
                  DO i=IstrT,IendT
                    ad_zeta_obc(i,ib,ir,Linp)=HnormRobc(i,ib)*          &
     &                                        ad_zeta_obc(i,ib,ir,Linp)
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
!
!  Adjoint 2D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
# ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_ubar_obc(:,ib,ir,Linp))
!
# endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO j=JstrT,JendT
                    ad_ubar_obc(j,ib,ir,Linp)=HnormUobc(j,ib)*          &
     &                                        ad_ubar_obc(j,ib,ir,Linp)
                  END DO
                CASE (isouth, inorth)
                  DO i=IstrP,IendT
                    ad_ubar_obc(i,ib,ir,Linp)=HnormUobc(i,ib)*          &
     &                                        ad_ubar_obc(i,ib,ir,Linp)
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
!
!  Adjoint 2D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
# ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_vbar_obc(:,ib,ir,Linp))
!
# endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO j=JstrP,JendT
                    ad_vbar_obc(j,ib,ir,Linp)=HnormVobc(j,ib)*          &
     &                                        ad_vbar_obc(j,ib,ir,Linp)
                  END DO
                CASE (isouth, inorth)
                  DO i=IstrT,IendT
                    ad_vbar_obc(i,ib,ir,Linp)=HnormVobc(i,ib)*          &
     &                                        ad_vbar_obc(i,ib,ir,Linp)
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO

# ifdef SOLVE3D
!
!  Adjoint 3D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij, 1, N(ng),            &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_u_obc(:,:,ib,ir,Linp))
!
#  endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO k=1,N(ng)
                    DO j=JstrT,JendT
                      ad_u_obc(j,k,ib,ir,Linp)=VnormUobc(j,k,ib)*       &
     &                                         ad_u_obc(j,k,ib,ir,Linp)
                    END DO
                  END DO
                CASE (isouth, inorth)
                  DO k=1,N(ng)
                    DO i=IstrP,IendT
                      ad_u_obc(i,k,ib,ir,Linp)=VnormUobc(i,k,ib)*       &
     &                                         ad_u_obc(i,k,ib,ir,Linp)
                    END DO
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
!
!  Adjoint 3D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij, 1, N(ng),            &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_v_obc(:,:,ib,ir,Linp))
!
#  endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO k=1,N(ng)
                    DO j=JstrP,JendT
                      ad_v_obc(j,k,ib,ir,Linp)=VnormVobc(j,k,ib)*       &
     &                                         ad_v_obc(j,k,ib,ir,Linp)
                    END DO
                  END DO
                CASE (isouth, inorth)
                  DO k=1,N(ng)
                    DO i=IstrT,IendT
                      ad_v_obc(i,k,ib,ir,Linp)=VnormVobc(i,k,ib)*       &
     &                                         ad_v_obc(i,k,ib,ir,Linp)
                    END DO
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
!
!  Adjoint tracer variables open boundaries.
!
      DO itrc=1,NT(ng)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,isTvar(itrc),ng)) THEN
#  ifdef DISTRIBUTE
              CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ib,        &
     &                                   LBij, UBij, 1, N(ng),          &
     &                                   NghostPoints,                  &
     &                                   EWperiodic(ng), NSperiodic(ng),&
     &                                   ad_t_obc(:,:,ib,ir,Linp,itrc))
!
#  endif
              IF (Lconvolve(ib)) THEN
                SELECT CASE (ib)
                  CASE (iwest, ieast)
                    DO k=1,N(ng)
                      DO j=JstrT,JendT
                        ad_t_obc(j,k,ib,ir,Linp,itrc)=                  &
     &                                  VnormRobc(j,k,ib,itrc)*         &
     &                                  ad_t_obc(j,k,ib,ir,Linp,itrc)
                      END DO
                    END DO
                  CASE (isouth, inorth)
                    DO k=1,N(ng)
                      DO i=IstrT,IendT
                        ad_t_obc(i,k,ib,ir,Linp,itrc)=                  &
     &                                  VnormRobc(i,k,ib,itrc)*         &
     &                                  ad_t_obc(i,k,ib,ir,Linp,itrc)
                      END DO
                    END DO
                END SELECT
              END IF
            END IF
          END DO
        END DO
      END DO
# endif
#endif

#ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress.
!
      IF (.not.Lweak) THEN
# ifdef DISTRIBUTE
        CALL ad_mp_exchange3d (ng, tile, model, 2,                      &
     &                         LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_ustr(:,:,:,Linp),                     &
     &                         ad_vstr(:,:,:,Linp))
# endif
        DO ir=1,Nfrec(ng)
          DO j=JstrT,JendT
            DO i=IstrP,IendT
              ad_ustr(i,j,ir,Linp)=HnormSUS(i,j)*ad_ustr(i,j,ir,Linp)
            END DO
          END DO
!
          DO j=JstrP,JendT
            DO i=IstrT,IendT
              ad_vstr(i,j,ir,Linp)=HnormSVS(i,j)*ad_vstr(i,j,ir,Linp)
            END DO
          END DO
        END DO
      END IF
#endif

#ifdef ADJUST_STFLUX
!
!  Adjoint surface tracers flux.
!
      IF (.not.Lweak) THEN
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
# ifdef DISTRIBUTE
            CALL ad_mp_exchange3d (ng, tile, model, 1,                  &
     &                             LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                             NghostPoints,                        &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_tflux(:,:,:,Linp,itrc))
!
# endif
            DO ir=1,Nfrec(ng)
              DO j=JstrT,JendT
                DO i=IstrT,IendT
                  ad_tflux(i,j,ir,Linp,itrc)=HnormSTF(i,j,itrc)*        &
     &                                       ad_tflux(i,j,ir,Linp,itrc)
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Prior/model multiscale error convariance: Convolve adjoint state
!  vector with a generalized implicit diffusion operators to filter
!  solution with specified horizontal scales. Convert from minimization
!  space (v-space) to model space.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Adjoint free-surface:  CG/CI horizontal convolution.
!
      CALL self%ad_CI_2d (ng, tile, model, isFsur, r2dvar,              &
     &                    ns, NiterCI(ns,ng), ifac, Lweak,              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ad_zeta(:,:,Linp))

#ifndef SOLVE3D
!
!  Adjoint 2D U-momentum:  CG/CI horizontal convolution.
!
      CALL self%ad_CI_2d (ng, tile, model, isUbar, u2dvar,              &
     &                    ns, NiterCI(ns,ng), ifac, Lweak,              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ad_ubar(:,:,Linp))
!
!  Adjoint 2D V-momentum:  CG/CI horizontal convolution.
!
      CALL self%ad_CI_2d (ng, tile, model, isVbar, v2dvar,              &
     &                    ns, NiterCI(ns,ng), ifac, Lweak,              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ad_vbar(:,:,Linp))
#endif

#ifdef SOLVE3D
!
!  Adjoint 3D U-momentum:  Implicit vertical diffusion and
!                          CG/CI horizontal convolution.
!
      CALL self%ad_Vdiff_u3d (ng, tile, model, isUvel,                  &
     &                        NVsteps(rec,isUvel)/ifac,                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        DTsizeV(rec,isUvel),                      &
     &                        MIXING(ng) % Kv,                          &
     &                        ad_u(:,:,:,Linp))
!
      CALL self%ad_CI_3d (ng, tile, model, isUvel, u3dvar,              &
     &                    ns, NiterCI(ns,ng), ifac, Lweak,              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ad_u(:,:,:,Linp))
!
!  Adjoint 3D V-momentum:  Implicit vertical diffusion and
!                          CG/CI horizontal convolution.
!
      CALL self%ad_Vdiff_v3d (ng, tile, model, isVvel,                  &
     &                        NVsteps(rec,isUvel)/ifac,                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        DTsizeV(rec,isUvel),                      &
     &                        MIXING(ng) % Kv,                          &
     &                        ad_v(:,:,:,Linp))
!
      CALL self%ad_CI_3d (ng, tile, model, isVvel, v3dvar,              &
     &                    ns, NiterCI(ns,ng), ifac, Lweak,              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ad_v(:,:,:,Linp))
!
!  Adjoint tracer variables:  Implicit vertical diffusion and
!                             CG/CI horizontal convolution.
!

      TRACER_LOOP : DO itrc=1,NT(ng)
        ifield=isTvar(itrc)
        CALL self%ad_Vdiff_r3d (ng, tile, model, ifield,                &
     &                          NVsteps(rec,ifield)/ifac,               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          DTsizeV(rec,ifield),                    &
     &                          MIXING(ng) % Kv,                        &
     &                          ad_t(:,:,:,Linp,itrc))
!
        CALL self%ad_CI_3d (ng, tile, model, ifield, r3dvar,            &
     &                      ns, NiterCI(ns,ng), ifac, Lweak,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      ad_t(:,:,:,Linp,itrc))
      END DO TRACER_LOOP
#endif

#ifdef ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Open boundaries error convariance: Convolve adjoint boundary edges
!  with a generalized diffusion equation to filter solution with
!  specified horizontal scales. Convert from minimization space
!  (v-space) to model space.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Adjoint free-surface open boundaries: CG/CI horizontal convolution.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
            CALL self%ad_CI_b1d (ng, tile, model, isFsur, ib,           &
     &                           r2dvar, ns, NiterCI(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ad_zeta_obc(:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Adjoint 2D U-momentum open boundaries: CG/CI horizontal convolution.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
            CALL self%ad_CI_b1d (ng, tile, model, isUbar, ib,           &
     &                           u2dvar, ns, NiterCI(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ad_ubar_obc(:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Adjoint 2D V-momentum open boundaries: CG/CI horizontal convolution.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
            CALL self%ad_CI_b1d (ng, tile, model, isVbar, ib,           &
     &                           v2dvar, ns, NiterCI(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ad_vbar_obc(:,ib,ir,Linp))
          END IF
        END DO
      END DO

# ifdef SOLVE3D
!
!  Adjoint 3D U-momentum open boundaries: Implicit vertical diffusion
!  and CG/CI horizontal convolution.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
            CALL self%ad_bry_Vdiff (ng, tile, model, isUvel, ib,        &
     &                              u3dvar, NVstepsB(ib,isUvel)/ifac,   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              LBij, UBij,                         &
     &                              DTsizeVB(ib,isUvel),                &
     &                              MIXING(ng) % Kv,                    &
     &                              ad_u_obc(:,:,ib,ir,Linp))
!
            CALL self%ad_CI_b2d (ng, tile, model, isUvel, ib,           &
     &                           u3dvar, ns, NiterCI(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ad_u_obc(:,:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Adjoint 3D V-momentum open boundaries: Implicit vertical diffusion
!  and CG/CI horizontal convolution.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
            CALL self%ad_bry_Vdiff (ng, tile, model, isVvel, ib,        &
     &                              v3dvar, NVstepsB(ib,isUvel)/ifac,   &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              LBij, UBij,                         &
     &                              DTsizeVB(ib,isVvel),                &
     &                              MIXING(ng) % Kv,                    &
     &                              ad_v_obc(:,:,ib,ir,Linp))
!
            CALL self%ad_CI_b2d (ng, tile, model, isVvel, ib,           &
     &                           v3dvar, ns, NiterCI(ns,ng), ifac,      &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ad_v_obc(:,:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Adjoint tracer varibles open boundaries: Implicit vertical diffusion
!  and CG/CI horizontal convolution.
!
      BRY_TRACER_LOOP : DO itrc=1,NT(ng)
        ifield=isTvar(itrc)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,ifield,ng)) THEN
              CALL self%ad_bry_Vdiff (ng, tile, model, ifield, ib,      &
     &                                r3dvar, NVstepsB(ib,ifield)/ifac, &
     &                                LBi, UBi, LBj, UBj,               &
     &                                LBij, UBij,                       &
     &                                DTsizeVB(ib,ifield),              &
     &                                MIXING(ng) % Kv,                  &
     &                                ad_t_obc(:,:,ib,ir,Linp,itrc))
!
              CALL self%ad_CI_b2d (ng, tile, model, ifield, ib,         &
     &                             r3dvar, ns, NiterCI(ns,ng), ifac,    &
     &                             LBij, UBij,                          &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             ad_t_obc(:,:,ib,ir,Linp,itrc))
!
            END IF
          END DO
        END DO
      END DO BRY_TRACER_LOOP
# endif
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Surface forcing error convariance: Convolve adjoint state vector
!  with a generalized diffusion equation to filter solution with
!  specified horizontal scales. Convert from minimization space
!  (v-space) to model space.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress: CG/CI horizontal convolution.
!
      IF (.not.Lweak) THEN
        DO ir=1,Nfrec(ng)
          CALL self%ad_CI_2d (ng, tile, model, isUstr, u2dvar,          &
     &                        ns, NiterCI(ns,ng), ifac, Lweak,          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        ad_ustr(:,:,ir,Linp))
!
          CALL self%ad_CI_2d (ng, tile, model, isVstr, v2dvar,          &
     &                        ns, NiterCI(ns,ng), ifac, Lweak,          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        ad_vstr(:,:,ir,Linp))
        END DO
      END IF
# endif

# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Adjoint surface tracers flux: CG/CI horizontal convolution.
!
      IF (.not.Lweak) THEN
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            ifield=isTsur(itrc)
            DO ir=1,Nfrec(ng)
              CALL self%ad_CI_2d (ng, tile, model, ifield, r2dvar,      &
     &                            ns, NiterCI(ns,ng), ifac, Lweak,      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            ad_tflux(:,:,ir,Linp,itrc))
            END DO
          END IF
        END DO
      END IF
# endif
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Multiply adjoint state by the inverse squared root of its associated
!  area (2D) or volume (3D).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Adjoint free-surface.
!
#ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,Linp))
!
#endif
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)/                          &
     &                      SQRT(GRID(ng)%om_r(i,j)*                    &
     &                           GRID(ng)%on_r(i,j))
#ifdef MASKING
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)*GRID(ng)%rmask(i,j)
#endif
        END DO
      END DO

#ifndef SOLV3D
!
!  Adjoint 2D momentum.
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ubar(:,:,Linp),                         &
     &                       ad_vbar(:,:,Linp))
!
# endif
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)/                          &
     &                      SQRT(GRID(ng)%om_u(i,j)*                    &
     &                           GRID(ng)%on_u(i,j))
# ifdef MASKING
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)*GRID(ng)%umask(i,j)
# endif
        END DO
      END DO
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          ad_vbar(i,j,Linp)=ad_vbar(i,j,Linp)/                          &
     &                      SQRT(GRID(ng)%om_v(i,j)*                    &
     &                           GRID(ng)%on_v(i,j))
# ifdef MASKING
          ad_vbar(i,j,Linp)=GRID(ng)%vmask(i,j)*ad_vbar(i,j,Linp)
# endif
        END DO
      END DO
#endif

#ifdef SOLVE3D
!
!  Adjoint 3D momentum.
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_u(:,:,:,Linp),                          &
     &                       ad_v(:,:,:,Linp))
!
# endif
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          cff=GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*0.5_r8
          DO k=1,N(ng)
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)/                          &
     &                       SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+            &
     &                                 GRID(ng)%Hz(i  ,j,k)))
# ifdef MASKING
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)*GRID(ng)%umask(i,j)
# endif
          END DO
        END DO
      END DO
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          cff=GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*0.5_r8
          DO k=1,N(ng)
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)/                          &
     &                       SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+            &
     &                                 GRID(ng)%Hz(i,j  ,k)))
# ifdef MASKING
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)*GRID(ng)%vmask(i,j)
# endif
          END DO
        END DO
      END DO
!
!  Adjoint tracer variables.
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_t(:,:,:,Linp,:))
!
# endif
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          cff=GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j)
          DO k=1,N(ng)
            fac=1.0_r8/SQRT(cff*GRID(ng)%Hz(i,j,k))
            DO itrc=1,NT(ng)
              ad_t(i,j,k,Linp,itrc)=fac*ad_t(i,j,k,Linp,itrc)
# ifdef MASKING
              ad_t(i,j,k,Linp,itrc)=ad_t(i,j,k,Linp,itrc)*              &
     &                              GRID(ng)%rmask(i,j)
# endif
            END DO
          END DO
        END DO
      END DO
#endif

#ifdef ADJUST_BOUNDARY
!
!  Adjoint free-surface open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
# ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_zeta_obc(:,ib,ir,Linp))
!
# endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,r2dvar)
                  DO j=JstrT,JendT
                    ad_zeta_obc(j,ib,ir,Linp)=                          &
     &                                  ad_zeta_obc(j,ib,ir,Linp)/      &
     &                                  SQRT(GRID(ng)%on_r(i,j))
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,r2dvar)
                  DO i=IstrT,IendT
                    ad_zeta_obc(i,ib,ir,Linp)=                          &
     &                                  ad_zeta_obc(i,ib,ir,Linp)/      &
     &                                  SQRT(GRID(ng)%om_r(i,j))

                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
!
!  Adjoint 2D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
# ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_ubar_obc(:,ib,ir,Linp))
!
# endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,u2dvar)
                  DO j=JstrT,JendT
                    ad_ubar_obc(j,ib,ir,Linp)=                          &
     &                                  ad_ubar_obc(j,ib,ir,Linp)/      &
     &                                  SQRT(GRID(ng)%on_u(i,j))
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,u2dvar)
                  DO i=IstrP,IendT
                    ad_ubar_obc(i,ib,ir,Linp)=                          &
     &                                  ad_ubar_obc(i,ib,ir,Linp)/      &
     &                                  SQRT(GRID(ng)%om_u(i,j))
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
!
!  Adjoint 2D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
# ifdef DISTRIBUTE
            CALL ad_mp_exchange2d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij,                      &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_vbar_obc(:,ib,ir,Linp))
# endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,v2dvar)
                  DO j=JstrP,JendT
                    ad_vbar_obc(j,ib,ir,Linp)=                          &
     &                                  ad_vbar_obc(j,ib,ir,Linp)/      &
     &                                  SQRT(GRID(ng)%on_v(i,j))
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,v2dvar)
                  DO i=IstrT,IendT
                    ad_vbar_obc(i,ib,ir,Linp)=                          &
     &                                  ad_vbar_obc(i,ib,ir,Linp)/      &
     &                                  SQRT(GRID(ng)%om_v(i,j))
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO

# ifdef SOLVE3D
!
!  Adjoint 3D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij, 1, N(ng),            &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_u_obc(:,:,ib,ir,Linp))
!
#  endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,u2dvar)
                  DO j=JstrT,JendT
                    cff=GRID(ng)%on_u(i,j)*0.5_r8
                    DO k=1,N(ng)
                      ad_u_obc(j,k,ib,ir,Linp)=                         &
     &                                ad_u_obc(j,k,ib,ir,Linp)/         &
     &                                SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+   &
     &                                          GRID(ng)%Hz(i  ,j,k)))
                    END DO
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,u2dvar)
                  DO i=IstrP,IendT
                    cff=GRID(ng)%om_u(i,j)*0.5_r8
                    DO k=1,N(ng)
                      ad_u_obc(i,k,ib,ir,Linp)=                         &
     &                                ad_u_obc(i,k,ib,ir,Linp)/         &
     &                                SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+   &
     &                                          GRID(ng)%Hz(i  ,j,k)))
                    END DO
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
!
!  Adjoint 3D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
#  ifdef DISTRIBUTE
            CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ib,          &
     &                                 LBij, UBij, 1, N(ng),            &
     &                                 NghostPoints,                    &
     &                                 EWperiodic(ng), NSperiodic(ng),  &
     &                                 ad_v_obc(:,:,ib,ir,Linp))
!
#  endif
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,v2dvar)
                  DO j=JstrP,JendT
                    cff=GRID(ng)%on_v(i,j)*0.5_r8
                    DO k=1,N(ng)
                      ad_v_obc(j,k,ib,ir,Linp)=                         &
     &                                ad_v_obc(j,k,ib,ir,Linp)/         &
     &                                SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+   &
     &                                          GRID(ng)%Hz(i,j  ,k)))
                    END DO
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,v2dvar)
                  DO i=IstrT,IendT
                    cff=GRID(ng)%om_v(i,j)*0.5_r8
                    DO k=1,N(ng)
                      ad_v_obc(i,k,ib,ir,Linp)=                         &
     &                                ad_v_obc(i,k,ib,ir,Linp)/         &
     &                                SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+   &
     &                                          GRID(ng)%Hz(i,j  ,k)))
                    END DO
                  END DO
              END SELECT
            END IF
          END IF
        END DO
      END DO
!
!  Adjoint tracer variables open boundaries.
!
      DO itrc=1,NT(ng)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,isTvar(itrc),ng)) THEN
#  ifdef DISTRIBUTE
              CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ib,        &
     &                                   LBij, UBij, 1, N(ng),          &
     &                                   NghostPoints,                  &
     &                                   EWperiodic(ng), NSperiodic(ng),&
     &                                   ad_t_obc(:,:,ib,ir,Linp,itrc))
#  endif
              IF (Lconvolve(ib)) THEN
                SELECT CASE (ib)
                  CASE (iwest, ieast)
                    i=BOUNDS(ng)%edge(ib,r2dvar)
                    DO j=JstrT,JendT
                      cff=GRID(ng)%on_r(i,j)
                      DO k=1,N(ng)
                        ad_t_obc(j,k,ib,ir,Linp,itrc)=                  &
     &                                  ad_t_obc(j,k,ib,ir,Linp,itrc)/  &
     &                                  SQRT(cff*GRID(ng)%Hz(i,j,k))
                      END DO
                    END DO
                  CASE (isouth, inorth)
                    j=BOUNDS(ng)%edge(ib,r2dvar)
                    DO i=IstrT,IendT
                      cff=GRID(ng)%om_r(i,j)
                      DO k=1,N(ng)
                        ad_t_obc(i,k,ib,ir,Linp,itrc)=                  &
     &                                  ad_t_obc(i,k,ib,ir,Linp,itrc)/  &
     &                                  SQRT(cff*GRID(ng)%Hz(i,j,k))
                      END DO
                    END DO
                END SELECT
              END IF
            END IF
          END DO
        END DO
      END DO
# endif
#endif

#ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress.
!
      IF (.not.Lweak) THEN
# ifdef DISTRIBUTE
        CALL ad_mp_exchange3d (ng, tile, model, 2,                      &
     &                         LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_ustr(:,:,:,Linp),                     &
     &                         ad_vstr(:,:,:,Linp))
!
# endif
        DO ir=1,Nfrec(ng)
          DO j=JstrT,JendT
            DO i=IstrP,IendT
              ad_ustr(i,j,ir,Linp)=ad_ustr(i,j,ir,Linp)/                &
     &                             SQRT(GRID(ng)%om_u(i,j)*             &
     &                                  GRID(ng)%on_u(i,j))
            END DO
          END DO
!
          DO j=JstrP,JendT
            DO i=IstrT,IendT
              ad_vstr(i,j,ir,Linp)=ad_vstr(i,j,ir,Linp)/                &
     &                             SQRT(GRID(ng)%om_v(i,j)*             &
     &                                  GRID(ng)%on_v(i,j))
            END DO
          END DO
        END DO
      END IF
#endif

#if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Adjoint surface tracers flux.
!
      IF (.not.Lweak) THEN
# ifdef DISTRIBUTE
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            CALL ad_mp_exchange3d (ng, tile, model, 1,                  &
     &                             LBi, UBi, LBj, UBj, 1, Nfrec(ng),    &
     &                             NghostPoints,                        &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_tflux(:,:,:,Linp,itrc))
          END IF
        END DO
!
# endif
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            fac=1.0_r8/SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            DO itrc=1,NT(ng)
              IF (Lstflux(itrc,ng)) THEN
                DO ir=1,Nfrec(ng)
                  ad_tflux(i,j,ir,Linp,itrc)=fac*                       &
     &                                       ad_tflux(i,j,ir,Linp,itrc)
                END DO
              END IF
            END DO
          END DO
        END DO
      END IF
#endif
!
      RETURN
      END SUBROUTINE ad_convolution_tile
!
      END MODULE ad_convolution_mod
