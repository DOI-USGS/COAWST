#include "cppdefs.h"
      MODULE tl_convolution_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group            Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Mono-scale, Background-Error Covariace modeling:                    !
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
!                                                                      !
!  Spreading is modeled using pseudo-diffusion operators in spatial    !
!  correlation space, where diffusion coefficients (K) are proportional!
!  to the square of the correlation length scale (Daley, 1992).        !
!  Correlation scales may be constant or spatially varying for each    !
!  variable in the control vector.                                     !
!                                                                      !
!  The monoscale approach uses a single spatial scale to represent B.  !
!  Error correlations are considered separable in the horizontal and   !
!  vertical directions. The horizontal pseudo-diffusion operator is    !
!  explicit while the vertical operator is implicit.                   !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Weaver, A. and P. Courtier, 2001: Correlation modeling on the     !
!      sphere using a generalized diffusion equation, Q.J.R. Meteorol. !
!      Soc, 127, 1815-1846, doi:10.1002/qj.49712757518.                !
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
#ifdef DISTRIBUTE
      USE mp_exchange_mod,   ONLY : mp_exchange2d,                      &
     &                              mp_exchange3d,                      &
     &                              mp_exchange4d
# ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,   ONLY : mp_exchange2d_bry,                  &
     &                              mp_exchange3d_bry
# endif
#endif
#ifdef SOLVE3D
      USE set_depth_mod,     ONLY : set_depth_tile
#endif
      USE tl_conv_2d_mod,    ONLY : tl_conv_r2d_tile,                   &
     &                              tl_conv_u2d_tile,                   &
     &                              tl_conv_v2d_tile
#ifdef SOLVE3D
      USE tl_conv_3d_mod,    ONLY : tl_conv_r3d_tile,                   &
     &                              tl_conv_u3d_tile,                   &
     &                              tl_conv_v3d_tile
#endif
#ifdef ADJUST_BOUNDARY
      USE tl_conv_bry2d_mod, ONLY : tl_conv_r2d_bry_tile,               &
     &                              tl_conv_u2d_bry_tile,               &
     &                              tl_conv_v2d_bry_tile
# ifdef SOLVE3D
      USE tl_conv_bry3d_mod, ONLY : tl_conv_r3d_bry_tile,               &
     &                              tl_conv_u3d_bry_tile,               &
     &                              tl_conv_v3d_bry_tile
# endif
#endif
!
      implicit none
!
      PUBLIC  :: tl_convolution
      PRIVATE :: tl_convolution_tile
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_convolution (ng, tile, Linp, Lweak, ifac)
!***********************************************************************
!
      USE mod_stepping, ONLY : nnew, nstp
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak

      integer, intent(in) :: ng, tile, Linp, ifac
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL tl_convolution_tile (ng, tile, iTLM,                         &
     &                          Lweak, ifac,                            &
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
     &                          BOUNDARY(ng) % tl_ubar_obc,             &
     &                          BOUNDARY(ng) % tl_vbar_obc,             &
     &                          BOUNDARY(ng) % tl_zeta_obc,             &
# ifdef SOLVE3D
     &                          BOUNDARY(ng) % tl_t_obc,                &
     &                          BOUNDARY(ng) % tl_u_obc,                &
     &                          BOUNDARY(ng) % tl_v_obc,                &
# endif
#endif
#ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % tl_ustr,                   &
     &                          FORCES(ng) % tl_vstr,                   &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % tl_tflux,                  &
#endif
#ifdef SOLVE3D
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
#endif
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
     &                          OCEAN(ng) % tl_zeta)

      RETURN
      END SUBROUTINE tl_convolution
!
!***********************************************************************
      SUBROUTINE tl_convolution_tile (ng, tile, model,                  &
     &                                Lweak, ifac,                      &
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
     &                                tl_ubar_obc, tl_vbar_obc,         &
     &                                tl_zeta_obc,                      &
# ifdef SOLVE3D
     &                                tl_t_obc, tl_u_obc, tl_v_obc,     &
# endif
#endif
#ifdef ADJUST_WSTRESS
     &                                tl_ustr, tl_vstr,                 &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                                tl_tflux,                         &
#endif
#ifdef SOLVE3D
     &                                tl_t, tl_u, tl_v,                 &
#endif
     &                                tl_ubar, tl_vbar,                 &
     &                                tl_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: nstp, nnew, Linp, ifac
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
      real(r8), intent(in) :: HnormR(LBi:,LBj:,:)
      real(r8), intent(in) :: HnormU(LBi:,LBj:,:)
      real(r8), intent(in) :: HnormV(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(in) :: VnormR(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: VnormU(LBi:,LBj:,:,:)
      real(r8), intent(in) :: VnormV(LBi:,LBj:,:,:)
# endif
# ifdef ADJUST_BOUNDARY
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#  endif
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
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
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
# endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,:)
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,:)
# ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
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
!  Multiply tangent linear state by the inverse squared root of its
!  associated area (2D) or volume (3D).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Tangent linear free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)/                          &
     &                      SQRT(GRID(ng)%om_r(i,j)*                    &
     &                           GRID(ng)%on_r(i,j))
        END DO
      END DO
!
!  Tangent linear 2D momentum.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          tl_ubar(i,j,Linp)=tl_ubar(i,j,Linp)/                          &
     &                      SQRT(GRID(ng)%om_u(i,j)*                    &
     &                           GRID(ng)%on_u(i,j))
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          tl_vbar(i,j,Linp)=tl_vbar(i,j,Linp)/                          &
     &                      SQRT(GRID(ng)%om_v(i,j)*                    &
     &                           GRID(ng)%on_v(i,j))
        END DO
      END DO
#ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_zeta(:,:,Linp),                            &
     &                    tl_ubar(:,:,Linp),                            &
     &                    tl_vbar(:,:,Linp))
#endif

#ifdef SOLVE3D
!
!  Tangent linear 3D momentum.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          cff=GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*0.5_r8
          DO k=1,N(ng)
            tl_u(i,j,k,Linp)=tl_u(i,j,k,Linp)/                          &
     &                       SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+            &
     &                                 GRID(ng)%Hz(i  ,j,k)))
          END DO
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          cff=GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*0.5_r8
          DO k=1,N(ng)
            tl_v(i,j,k,Linp)=tl_v(i,j,k,Linp)/                          &
     &                       SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+            &
     &                                 GRID(ng)%Hz(i,j  ,k)))
          END DO
        END DO
      END DO
# ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_u(:,:,:,Linp),                             &
     &                    tl_v(:,:,:,Linp))
# endif
!
!  Tangent linear tracers.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          cff=GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j)
          DO k=1,N(ng)
            fac=1.0_r8/SQRT(cff*GRID(ng)%Hz(i,j,k))
            DO itrc=1,NT(ng)
              tl_t(i,j,k,Linp,itrc)=fac*tl_t(i,j,k,Linp,itrc)
            END DO
          END DO
        END DO
      END DO
# ifdef DISTRIBUTE
!
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_t(:,:,:,Linp,:))
# endif
#endif

#ifdef ADJUST_BOUNDARY
!
!  Tangent linear free-surface open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,r2dvar)
                  DO j=JstrT,JendT
                    tl_zeta_obc(j,ib,ir,Linp)=tl_zeta_obc(j,ib,ir,Linp)/&
     &                                        SQRT(GRID(ng)%on_r(i,j))
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,r2dvar)
                  DO i=IstrT,IendT
                    tl_zeta_obc(i,ib,ir,Linp)=tl_zeta_obc(i,ib,ir,Linp)/&
     &                                        SQRT(GRID(ng)%om_r(i,j))
                  END DO
              END SELECT
            END IF
# ifdef DISTRIBUTE
!
            CALL mp_exchange2d_bry (ng, tile, model, 1, ib,             &
     &                              LBij, UBij,                         &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_zeta_obc(:,ib,ir,Linp))
# endif
          END IF
        END DO
      END DO
!
!  Tangent linear 2D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,u2dvar)
                  DO j=JstrT,JendT
                    tl_ubar_obc(j,ib,ir,Linp)=tl_ubar_obc(j,ib,ir,Linp)/&
     &                                        SQRT(GRID(ng)%on_u(i,j))
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,u2dvar)
                  DO i=IstrP,IendT
                    tl_ubar_obc(i,ib,ir,Linp)=tl_ubar_obc(i,ib,ir,Linp)/&
     &                                        SQRT(GRID(ng)%om_u(i,j))
                  END DO
              END SELECT
            END IF
# ifdef DISTRIBUTE
!
            CALL mp_exchange2d_bry (ng, tile, model, 1, ib,             &
     &                              LBij, UBij,                         &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_ubar_obc(:,ib,ir,Linp))
# endif
          END IF
        END DO
      END DO
!
!  Tangent linear 2D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,v2dvar)
                  DO j=JstrP,JendT
                    tl_vbar_obc(j,ib,ir,Linp)=tl_vbar_obc(j,ib,ir,Linp)/&
     &                                        SQRT(GRID(ng)%on_v(i,j))
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,v2dvar)
                  DO i=IstrT,IendT
                    tl_vbar_obc(i,ib,ir,Linp)=tl_vbar_obc(i,ib,ir,Linp)/&
     &                                        SQRT(GRID(ng)%om_v(i,j))
                  END DO
              END SELECT
            END IF
# ifdef DISTRIBUTE
!
            CALL mp_exchange2d_bry (ng, tile, model, 1, ib,             &
     &                              LBij, UBij,                         &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_vbar_obc(:,ib,ir,Linp))
# endif
          END IF
        END DO
      END DO

# ifdef SOLVE3D
!
!  Tangent linear 3D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,u2dvar)
                  DO j=JstrT,JendT
                    cff=GRID(ng)%on_u(i,j)*0.5_r8
                    DO k=1,N(ng)
                      tl_u_obc(j,k,ib,ir,Linp)=                         &
     &                                tl_u_obc(j,k,ib,ir,Linp)/         &
     &                                SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+   &
     &                                          GRID(ng)%Hz(i  ,j,k)))
                    END DO
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,u2dvar)
                  DO i=IstrP,IendT
                    cff=GRID(ng)%om_u(i,j)*0.5_r8
                    DO k=1,N(ng)
                      tl_u_obc(i,k,ib,ir,Linp)=                         &
     &                                tl_u_obc(i,k,ib,ir,Linp)/         &
     &                                SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+   &
     &                                          GRID(ng)%Hz(i  ,j,k)))
                    END DO
                  END DO
              END SELECT
            END IF
#  ifdef DISTRIBUTE
!
            CALL mp_exchange3d_bry (ng, tile, model, 1, ib,             &
     &                              LBij, UBij, 1, N(ng),               &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_u_obc(:,:,ib,ir,Linp))
#  endif
          END IF
        END DO
      END DO
!
!  Tangent linear 3D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ib,v2dvar)
                  DO j=JstrP,JendT
                    cff=GRID(ng)%on_v(i,j)*0.5_r8
                    DO k=1,N(ng)
                      tl_v_obc(j,k,ib,ir,Linp)=                         &
     &                                tl_v_obc(j,k,ib,ir,Linp)/         &
     &                                SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+   &
     &                                          GRID(ng)%Hz(i,j  ,k)))
                    END DO
                  END DO
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ib,v2dvar)
                  DO i=IstrT,IendT
                    cff=GRID(ng)%om_v(i,j)*0.5_r8
                    DO k=1,N(ng)
                      tl_v_obc(i,k,ib,ir,Linp)=                         &
     &                                tl_v_obc(i,k,ib,ir,Linp)/         &
     &                                SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+   &
     &                                          GRID(ng)%Hz(i,j  ,k)))
                    END DO
                  END DO
              END SELECT
            END IF
#  ifdef DISTRIBUTE
!
            CALL mp_exchange3d_bry (ng, tile, model, 1, ib,             &
     &                              LBij, UBij, 1, N(ng),               &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_v_obc(:,:,ib,ir,Linp))
#  endif
          END IF
        END DO
      END DO
!
!  Tangent linear tracers open boundaries.
!
      BRY_TRACER_LOOP : DO itrc=1,NT(ng)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,isTvar(itrc),ng)) THEN
              IF (Lconvolve(ib)) THEN
                SELECT CASE (ib)
                  CASE (iwest, ieast)
                    i=BOUNDS(ng)%edge(ib,r2dvar)
                    DO j=JstrT,JendT
                      cff=GRID(ng)%on_r(i,j)
                      DO k=1,N(ng)
                        tl_t_obc(j,k,ib,ir,Linp,itrc)=                  &
     &                                   tl_t_obc(j,k,ib,ir,Linp,itrc)/ &
     &                                   SQRT(cff*GRID(ng)%Hz(i,j,k))
                      END DO
                    END DO
                  CASE (isouth, inorth)
                    j=BOUNDS(ng)%edge(ib,r2dvar)
                    DO i=IstrT,IendT
                      cff=GRID(ng)%om_r(i,j)
                      DO k=1,N(ng)
                        tl_t_obc(i,k,ib,ir,Linp,itrc)=                  &
     &                                   tl_t_obc(i,k,ib,ir,Linp,itrc)/ &
     &                                   SQRT(cff*GRID(ng)%Hz(i,j,k))
                      END DO
                    END DO
                END SELECT
              END IF
#  ifdef DISTRIBUTE
!
              CALL mp_exchange3d_bry (ng, tile, model, 1, ib,           &
     &                                LBij, UBij, 1, N(ng),             &
     &                                NghostPoints,                     &
     &                                EWperiodic(ng), NSperiodic(ng),   &
     &                                tl_t_obc(:,:,ib,ir,Linp,itrc))
#  endif
            END IF
          END DO
        END DO
      END DO BRY_TRACER_LOOP
# endif
#endif

#ifdef ADJUST_WSTRESS
!
!  Tangent linear surface momentum stress.
!
      IF (.not.Lweak) THEN
        DO ir=1,Nfrec(ng)
          DO j=JstrT,JendT
            DO i=IstrP,IendT
              tl_ustr(i,j,ir,Linp)=tl_ustr(i,j,ir,Linp)/                &
     &                             SQRT(GRID(ng)%om_u(i,j)*             &
     &                                  GRID(ng)%on_u(i,j))
            END DO
          END DO
          DO j=JstrP,JendT
            DO i=IstrT,IendT
              tl_vstr(i,j,ir,Linp)=tl_vstr(i,j,ir,Linp)/                &
     &                             SQRT(GRID(ng)%om_v(i,j)*             &
     &                                  GRID(ng)%on_v(i,j))
            END DO
          END DO
        END DO
# ifdef DISTRIBUTE
!
        CALL mp_exchange3d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_ustr(:,:,:,Linp),                        &
     &                      tl_vstr(:,:,:,Linp))
# endif
      END IF
#endif

#if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Tangent linear surface tracers flux.
!
      IF (.not.Lweak) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            fac=1.0_r8/SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            DO itrc=1,NT(ng)
              IF (Lstflux(itrc,ng)) THEN
                DO ir=1,Nfrec(ng)
                  tl_tflux(i,j,ir,Linp,itrc)=fac*                       &
     &                                       tl_tflux(i,j,ir,Linp,itrc)
                END DO
              END IF
            END DO
          END DO
        END DO
# ifdef DISTRIBUTE
!
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            CALL mp_exchange3d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj, 1, Nfrec(ng),       &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          tl_tflux(:,:,:,Linp,itrc))
          END IF
        END DO
# endif
      END IF
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Initial conditions and model error convariance: Convolve tangent
!  linear state vector with a generalized diffusion equation to filter
!  solution with specified horizontal scales. Convert from minimization
!  space (v-space) to model space.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Tangent linear free-surface.
!
      CALL tl_conv_r2d_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isFsur)/ifac,                  &
     &                       DTsizeH(rec,isFsur),                       &
     &                       MIXING(ng) % Kh,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % pmon_u,                         &
     &                       GRID(ng) % pnom_v,                         &
#ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
#endif
     &                       tl_zeta(:,:,Linp))
!
!  Tangent linear 2D momentum.
!
      CALL tl_conv_u2d_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isUbar)/ifac,                  &
     &                       DTsizeH(rec,isUbar),                       &
     &                       MIXING(ng) % Kh,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % pmon_r,                         &
     &                       GRID(ng) % pnom_p,                         &
#ifdef MASKING
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % pmask,                          &
#endif
     &                       tl_ubar(:,:,Linp))
!
      CALL tl_conv_v2d_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isVbar)/ifac,                  &
     &                       DTsizeH(rec,isVbar),                       &
     &                       MIXING(ng) % Kh,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % pmon_p,                         &
     &                       GRID(ng) % pnom_r,                         &
#ifdef MASKING
     &                       GRID(ng) % vmask,                          &
     &                       GRID(ng) % pmask,                          &
#endif
     &                       tl_vbar(:,:,Linp))
#ifdef SOLVE3D
!
!  Tangent linear 3D momentum.
!
      CALL tl_conv_u3d_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isUvel)/ifac,                  &
     &                       NVsteps(rec,isUvel)/ifac,                  &
     &                       DTsizeH(rec,isUvel),                       &
     &                       DTsizeV(rec,isUvel),                       &
     &                       MIXING(ng) % Kh,                           &
     &                       MIXING(ng) % Kv,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
# ifdef GEOPOTENTIAL_HCONV
     &                       GRID(ng) % on_r,                           &
     &                       GRID(ng) % om_p,                           &
# else
     &                       GRID(ng) % pmon_r,                         &
     &                       GRID(ng) % pnom_p,                         &
# endif
# ifdef MASKING
#  ifdef GEOPOTENTIAL_HCONV
     &                       GRID(ng) % pmask,                          &
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
#  else
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % pmask,                          &
#  endif
# endif
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % z_r,                            &
     &                       tl_u(:,:,:,Linp))
!
      CALL tl_conv_v3d_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       NghostPoints,                              &
     &                       NHsteps(rec,isVvel)/ifac,                  &
     &                       NVsteps(rec,isVvel)/ifac,                  &
     &                       DTsizeH(rec,isVvel),                       &
     &                       DTsizeV(rec,isVvel),                       &
     &                       MIXING(ng) % Kh,                           &
     &                       MIXING(ng) % Kv,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
# ifdef GEOPOTENTIAL_HCONV
     &                       GRID(ng) % on_p,                           &
     &                       GRID(ng) % om_r,                           &
# else
     &                       GRID(ng) % pmon_p,                         &
     &                       GRID(ng) % pnom_r,                         &
# endif
# ifdef MASKING
#  ifdef GEOPOTENTIAL_HCONV
     &                       GRID(ng) % pmask,                          &
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
#  else
     &                       GRID(ng) % vmask,                          &
     &                       GRID(ng) % pmask,                          &
#  endif
# endif
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % z_r,                            &
     &                       tl_v(:,:,:,Linp))
!
!  Tangent linear tracers.
!
      DO itrc=1,NT(ng)
        ifield=isTvar(itrc)
        CALL tl_conv_r3d_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         NghostPoints,                            &
     &                         NHsteps(rec,ifield)/ifac,                &
     &                         NVsteps(rec,ifield)/ifac,                &
     &                         DTsizeH(rec,ifield),                     &
     &                         DTsizeV(rec,ifield),                     &
     &                         MIXING(ng) % Kh,                         &
     &                         MIXING(ng) % Kv,                         &
     &                         GRID(ng) % pm,                           &
     &                         GRID(ng) % pn,                           &
# ifdef GEOPOTENTIAL_HCONV
     &                         GRID(ng) % on_u,                         &
     &                         GRID(ng) % om_v,                         &
# else
     &                         GRID(ng) % pmon_u,                       &
     &                         GRID(ng) % pnom_v,                       &
# endif
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
# endif
     &                         GRID(ng) % Hz,                           &
     &                         GRID(ng) % z_r,                          &
     &                         tl_t(:,:,:,Linp,itrc))
      END DO
#endif

#ifdef ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Open boundaries error convariance: Convolve tangent linear boundary
!  edges with a generalized diffusion equation to filter solution with
!  specified horizontal scales. Convert from minimization space
!  (v-space) to model space.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Tangent linear free-surface open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
            CALL tl_conv_r2d_bry_tile (ng, tile, model, ib,             &
     &                                 BOUNDS(ng)%edge(:,r2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isFsur)/ifac,        &
     &                                 DTsizeHB(ib,isFsur),             &
     &                                 MIXING(ng) % Kh,                 &
     &                                 GRID(ng) % pm,                   &
     &                                 GRID(ng) % pn,                   &
     &                                 GRID(ng) % pmon_u,               &
     &                                 GRID(ng) % pnom_v,               &
# ifdef MASKING
     &                                 GRID(ng) % rmask,                &
     &                                 GRID(ng) % umask,                &
     &                                 GRID(ng) % vmask,                &
# endif
     &                                 tl_zeta_obc(:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Tangent linear 2D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
            CALL tl_conv_u2d_bry_tile (ng, tile, model, ib,             &
     &                                 BOUNDS(ng)%edge(:,u2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isUbar)/ifac,        &
     &                                 DTsizeHB(ib,isUbar),             &
     &                                 MIXING(ng) % Kh,                 &
     &                                 GRID(ng) % pm,                   &
     &                                 GRID(ng) % pn,                   &
     &                                 GRID(ng) % pmon_r,               &
     &                                 GRID(ng) % pnom_p,               &
# ifdef MASKING
     &                                 GRID(ng) % umask,                &
     &                                 GRID(ng) % pmask,                &
# endif
     &                                 tl_ubar_obc(:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Tangent linear 2D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
            CALL tl_conv_v2d_bry_tile (ng, tile, model, ib,             &
     &                                 BOUNDS(ng)%edge(:,v2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isVbar)/ifac,        &
     &                                 DTsizeHB(ib,isVbar),             &
     &                                 MIXING(ng) % Kh,                 &
     &                                 GRID(ng) % pm,                   &
     &                                 GRID(ng) % pn,                   &
     &                                 GRID(ng) % pmon_p,               &
     &                                 GRID(ng) % pnom_r,               &
# ifdef MASKING
     &                                 GRID(ng) % vmask,                &
     &                                 GRID(ng) % pmask,                &
# endif
     &                                 tl_vbar_obc(:,ib,ir,Linp))
          END IF
        END DO
      END DO

# ifdef SOLVE3D
!
!  Tangent linear 3D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
            CALL tl_conv_u3d_bry_tile (ng, tile, model, ib,             &
     &                                 BOUNDS(ng)%edge(:,u2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isUvel)/ifac,        &
     &                                 NVstepsB(ib,isUvel)/ifac,        &
     &                                 DTsizeHB(ib,isUvel),             &
     &                                 DTsizeVB(ib,isUvel),             &
     &                                 MIXING(ng) % Kh,                 &
     &                                 MIXING(ng) % Kv,                 &
     &                                 GRID(ng) % pm,                   &
     &                                 GRID(ng) % pn,                   &
     &                                 GRID(ng) % pmon_r,               &
     &                                 GRID(ng) % pnom_p,               &
#  ifdef MASKING
     &                                 GRID(ng) % umask,                &
     &                                 GRID(ng) % pmask,                &
#  endif
     &                                 GRID(ng) % Hz,                   &
     &                                 GRID(ng) % z_r,                  &
     &                                 tl_u_obc(:,:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Tangent linear 3D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
            CALL tl_conv_v3d_bry_tile (ng, tile, model, ib,             &
     &                                 BOUNDS(ng)%edge(:,v2dvar),       &
     &                                 LBij, UBij,                      &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 NghostPoints,                    &
     &                                 NHstepsB(ib,isVvel)/ifac,        &
     &                                 NVstepsB(ib,isVvel)/ifac,        &
     &                                 DTsizeHB(ib,isVvel),             &
     &                                 DTsizeVB(ib,isVvel),             &
     &                                 MIXING(ng) % Kh,                 &
     &                                 MIXING(ng) % Kv,                 &
     &                                 GRID(ng) % pm,                   &
     &                                 GRID(ng) % pn,                   &
     &                                 GRID(ng) % pmon_p,               &
     &                                 GRID(ng) % pnom_r,               &
#  ifdef MASKING
     &                                 GRID(ng) % vmask,                &
     &                                 GRID(ng) % pmask,                &
#  endif
     &                                 GRID(ng) % Hz,                   &
     &                                 GRID(ng) % z_r,                  &
     &                                 tl_v_obc(:,:,ib,ir,Linp))
          END IF
        END DO
      END DO
!
!  Tangent linear tracers open boundaries.
!
      DO itrc=1,NT(ng)
        ifield=isTvar(itrc)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,ifield,ng)) THEN
              CALL tl_conv_r3d_bry_tile (ng, tile, model, ib,           &
     &                                   BOUNDS(ng)%edge(:,r2dvar),     &
     &                                   LBij, UBij,                    &
     &                                   LBi, UBi, LBj, UBj, 1, N(ng),  &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   NghostPoints,                  &
     &                                   NHstepsB(ib,ifield)/ifac,      &
     &                                   NVstepsB(ib,ifield)/ifac,      &
     &                                   DTsizeHB(ib,ifield),           &
     &                                   DTsizeVB(ib,ifield),           &
     &                                   MIXING(ng) % Kh,               &
     &                                   MIXING(ng) % Kv,               &
     &                                   GRID(ng) % pm,                 &
     &                                   GRID(ng) % pn,                 &
     &                                   GRID(ng) % pmon_u,             &
     &                                   GRID(ng) % pnom_v,             &
#  ifdef MASKING
     &                                   GRID(ng) % rmask,              &
     &                                   GRID(ng) % umask,              &
     &                                   GRID(ng) % vmask,              &
#  endif
     &                                   GRID(ng) % Hz,                 &
     &                                   GRID(ng) % z_r,                &
     &                                   tl_t_obc(:,:,ib,ir,Linp,itrc))
            END IF
          END DO
        END DO
      END DO
# endif
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Surface forcing error convariance: Convolve tangent linear state
!  vector with a generalized diffusion equation to filter solution with
!  specified horizontal scales. Convert from minimization space
!  (v-space) to model space.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ifdef ADJUST_WSTRESS
!
!  Tangent linear surface momentum stress.
!
      IF (.not.Lweak) THEN
        DO k=1,Nfrec(ng)
          CALL tl_conv_u2d_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           NghostPoints,                          &
     &                           NHsteps(rec,isUstr)/ifac,              &
     &                           DTsizeH(rec,isUstr),                   &
     &                           MIXING(ng) % Kh,                       &
     &                           GRID(ng) % pm,                         &
     &                           GRID(ng) % pn,                         &
     &                           GRID(ng) % pmon_r,                     &
     &                           GRID(ng) % pnom_p,                     &
#  ifdef MASKING
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % pmask,                      &
#  endif
     &                           tl_ustr(:,:,k,Linp))
!
          CALL tl_conv_v2d_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           NghostPoints,                          &
     &                           NHsteps(rec,isVstr)/ifac,              &
     &                           DTsizeH(rec,isVstr),                   &
     &                           MIXING(ng) % Kh,                       &
     &                           GRID(ng) % pm,                         &
     &                           GRID(ng) % pn,                         &
     &                           GRID(ng) % pmon_p,                     &
     &                           GRID(ng) % pnom_r,                     &
#  ifdef MASKING
     &                           GRID(ng) % vmask,                      &
     &                           GRID(ng) % pmask,                      &
#  endif
     &                           tl_vstr(:,:,k,Linp))
        END DO
      END IF
# endif

# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Tangent linear surface tracers flux.
!
      IF (.not.Lweak) THEN
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            ifield=isTsur(itrc)
            DO k=1,Nfrec(ng)
              CALL tl_conv_r2d_tile (ng, tile, model,                   &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               NghostPoints,                      &
     &                               NHsteps(rec,ifield)/ifac,          &
     &                               DTsizeH(rec,ifield),               &
     &                               MIXING(ng) % Kh,                   &
     &                               GRID(ng) % pm,                     &
     &                               GRID(ng) % pn,                     &
     &                               GRID(ng) % pmon_u,                 &
     &                               GRID(ng) % pnom_v,                 &
#  ifdef MASKING
     &                               GRID(ng) % rmask,                  &
     &                               GRID(ng) % umask,                  &
     &                               GRID(ng) % vmask,                  &
#  endif
     &                               tl_tflux(:,:,k,Linp,itrc))
            END DO
          END IF
        END DO
      END IF
# endif
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Multiply convolved tangent linear state by its corresponding
!  normalization factor.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Tangent linear free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          tl_zeta(i,j,Linp)=HnormR(i,j,rec)*tl_zeta(i,j,Linp)
        END DO
      END DO
!
!  Tangent linear 2D momentum.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          tl_ubar(i,j,Linp)=HnormU(i,j,rec)*tl_ubar(i,j,Linp)
        END DO
      END DO
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          tl_vbar(i,j,Linp)=HnormV(i,j,rec)*tl_vbar(i,j,Linp)
        END DO
      END DO
#ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_zeta(:,:,Linp),                            &
     &                    tl_ubar(:,:,Linp),                            &
     &                    tl_vbar(:,:,Linp))
#endif

#ifdef SOLVE3D
!
!  Tangent linear 3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            tl_u(i,j,k,Linp)=VnormU(i,j,k,rec)*tl_u(i,j,k,Linp)
          END DO
        END DO
!
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            tl_v(i,j,k,Linp)=VnormV(i,j,k,rec)*tl_v(i,j,k,Linp)
          END DO
        END DO
      END DO
# ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_u(:,:,:,Linp),                             &
     &                    tl_v(:,:,:,Linp))
# endif
!
!  Tangent linear tracer variables.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              tl_t(i,j,k,Linp,itrc)=VnormR(i,j,k,rec,itrc)*             &
     &                              tl_t(i,j,k,Linp,itrc)
            END DO
          END DO
        END DO
      END DO
# ifdef DISTRIBUTE
!
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_t(:,:,:,Linp,:))
# endif
#endif

#ifdef ADJUST_BOUNDARY
!
!  Tangent linear free-surface open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isFsur,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO j=JstrT,JendT
                    tl_zeta_obc(j,ib,ir,Linp)=HnormRobc(j,ib)*          &
     &                                        tl_zeta_obc(j,ib,ir,Linp)
                  END DO
                CASE (isouth, inorth)
                  DO i=IstrT,IendT
                    tl_zeta_obc(i,ib,ir,Linp)=HnormRobc(i,ib)*          &
     &                                        tl_zeta_obc(i,ib,ir,Linp)
                  END DO
              END SELECT
            END IF
# ifdef DISTRIBUTE
!
            CALL mp_exchange2d_bry (ng, tile, model, 1, ib,              &
     &                              LBij, UBij,                         &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_zeta_obc(:,ib,ir,Linp))
# endif
          END IF
        END DO
      END DO
!
!  Tangent linear 2D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUbar,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO j=JstrT,JendT
                    tl_ubar_obc(j,ib,ir,Linp)=HnormUobc(j,ib)*          &
     &                                        tl_ubar_obc(j,ib,ir,Linp)
                  END DO
                CASE (isouth, inorth)
                  DO i=IstrP,IendT
                    tl_ubar_obc(i,ib,ir,Linp)=HnormUobc(i,ib)*          &
     &                                        tl_ubar_obc(i,ib,ir,Linp)
                  END DO
              END SELECT
            END IF
# ifdef DISTRIBUTE
!
            CALL mp_exchange2d_bry (ng, tile, model, 1, ib,             &
     &                              LBij, UBij,                         &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_ubar_obc(:,ib,ir,Linp))
# endif
          END IF
        END DO
      END DO
!
!  Tangent linear 2D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVbar,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO j=JstrP,JendT
                    tl_vbar_obc(j,ib,ir,Linp)=HnormVobc(j,ib)*          &
     &                                        tl_vbar_obc(j,ib,ir,Linp)
                  END DO
                CASE (isouth, inorth)
                  DO i=IstrT,IendT
                    tl_vbar_obc(i,ib,ir,Linp)=HnormVobc(i,ib)*          &
     &                                        tl_vbar_obc(i,ib,ir,Linp)
                  END DO
              END SELECT
            END IF
# ifdef DISTRIBUTE
!
            CALL mp_exchange2d_bry (ng, tile, model, 1, ib,             &
     &                              LBij, UBij,                         &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_vbar_obc(:,ib,ir,Linp))
# endif
          END IF
        END DO
      END DO

# ifdef SOLVE3D
!
!  Tangent linear 3D U-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isUvel,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO k=1,N(ng)
                    DO j=JstrT,JendT
                      tl_u_obc(j,k,ib,ir,Linp)=VnormUobc(j,k,ib)*       &
     &                                         tl_u_obc(j,k,ib,ir,Linp)
                    END DO
                  END DO
                CASE (isouth, inorth)
                  DO k=1,N(ng)
                    DO i=IstrP,IendT
                      tl_u_obc(i,k,ib,ir,Linp)=VnormUobc(i,k,ib)*       &
     &                                         tl_u_obc(i,k,ib,ir,Linp)
                    END DO
                  END DO
              END SELECT
            END IF
#  ifdef DISTRIBUTE
!
            CALL mp_exchange3d_bry (ng, tile, model, 1, ib,             &
     &                              LBij, UBij, 1, N(ng),               &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_u_obc(:,:,ib,ir,Linp))
#  endif
          END IF
        END DO
      END DO
!
!  Tangent linear 3D V-momentum open boundaries.
!
      DO ir=1,Nbrec(ng)
        DO ib=1,4
          IF (.not.Lweak.and.Lobc(ib,isVvel,ng)) THEN
            IF (Lconvolve(ib)) THEN
              SELECT CASE (ib)
                CASE (iwest, ieast)
                  DO k=1,N(ng)
                    DO j=JstrP,JendT
                      tl_v_obc(j,k,ib,ir,Linp)=VnormVobc(j,k,ib)*       &
     &                                         tl_v_obc(j,k,ib,ir,Linp)
                    END DO
                  END DO
                CASE (isouth, inorth)
                  DO k=1,N(ng)
                    DO i=IstrT,IendT
                      tl_v_obc(i,k,ib,ir,Linp)=VnormVobc(i,k,ib)*       &
     &                                         tl_v_obc(i,k,ib,ir,Linp)
                    END DO
                  END DO
              END SELECT
            END IF
#  ifdef DISTRIBUTE
!
            CALL mp_exchange3d_bry (ng, tile, model, 1, ib,             &
     &                              LBij, UBij, 1, N(ng),               &
     &                              NghostPoints,                       &
     &                              EWperiodic(ng), NSperiodic(ng),     &
     &                              tl_v_obc(:,:,ib,ir,Linp))
#  endif
          END IF
        END DO
      END DO
!
!  Tangent linear tracer variables open boundaries.
!
      DO itrc=1,NT(ng)
        DO ir=1,Nbrec(ng)
          DO ib=1,4
            IF (.not.Lweak.and.Lobc(ib,isTvar(itrc),ng)) THEN
              IF (Lconvolve(ib)) THEN
                SELECT CASE (ib)
                  CASE (iwest, ieast)
                    DO k=1,N(ng)
                      DO j=JstrT,JendT
                        tl_t_obc(j,k,ib,ir,Linp,itrc)=                  &
     &                                     VnormRobc(j,k,ib,itrc)*      &
     &                                     tl_t_obc(j,k,ib,ir,Linp,itrc)
                      END DO
                    END DO
                  CASE (isouth, inorth)
                    DO k=1,N(ng)
                      DO i=IstrT,IendT
                        tl_t_obc(i,k,ib,ir,Linp,itrc)=                  &
     &                                     VnormRobc(i,k,ib,itrc)*      &
     &                                     tl_t_obc(i,k,ib,ir,Linp,itrc)
                      END DO
                    END DO
                END SELECT
              END IF
#  ifdef DISTRIBUTE
!
              CALL mp_exchange3d_bry (ng, tile, model, 1, ib,           &
     &                                LBij, UBij, 1, N(ng),             &
     &                                NghostPoints,                     &
     &                                EWperiodic(ng), NSperiodic(ng),   &
     &                                tl_t_obc(:,:,ib,ir,Linp,itrc))
#  endif
            END IF
          END DO
        END DO
      END DO
# endif
#endif

#ifdef ADJUST_WSTRESS
!
!  Tangent linear surface momentum stress.
!
      IF (.not.Lweak) THEN
        DO k=1,Nfrec(ng)
          DO j=JstrT,JendT
            DO i=IstrP,IendT
              tl_ustr(i,j,k,Linp)=HnormSUS(i,j)*tl_ustr(i,j,k,Linp)
            END DO
          END DO
          DO j=JstrP,JendT
            DO i=IstrT,IendT
              tl_vstr(i,j,k,Linp)=HnormSVS(i,j)*tl_vstr(i,j,k,Linp)
            END DO
          END DO
        END DO
# ifdef DISTRIBUTE
!
        CALL mp_exchange3d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_ustr(:,:,:,Linp),                        &
     &                      tl_vstr(:,:,:,Linp))
# endif
      END IF
#endif

#if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Tangent linear surface tracers flux.
!
      IF (.not.Lweak) THEN
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            DO k=1,Nfrec(ng)
              DO j=JstrT,JendT
                DO i=IstrT,IendT
                  tl_tflux(i,j,k,Linp,itrc)=HnormSTF(i,j,itrc)*         &
     &                                     tl_tflux(i,j,k,Linp,itrc)
                END DO
              END DO
            END DO
# ifdef DISTRIBUTE
!
            CALL mp_exchange3d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj, 1, Nfrec(ng),       &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          tl_tflux(:,:,:,Linp,itrc))
# endif
          END IF
        END DO
      END IF
#endif
!
      RETURN
      END SUBROUTINE tl_convolution_tile
!
      END MODULE tl_convolution_mod
