#include "cppdefs.h"
      MODULE ad_variability_mod

#if defined ADJOINT && defined FOUR_DVAR
!
!svn $Id: ad_variability.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine converts adjoint error correlations, C, to background  !
!  error covariances,  B,  by multiplying by the standard deviations,  !
!  S.                                                                  !
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
!    G : normalization coefficients matrix (convert to correlations)   !
!    L : tangent linear and adjoint diffusion operators                !
!    W : diagonal matrix of local area or volume metrics               !
!                                                                      !
!  Here, T/2 denote the transpose of a squared-root factor.            !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Weaver, A. and P. Courtier, 2001: Correlation modeling on the     !
!      sphere using a generalized diffusion equation, Q.J.R. Meteo.    !
!      Soc, 127, 1815-1846.                                            !
!                                                                      !
!======================================================================!
!
      USE mod_kinds

      implicit none

      PRIVATE
      PUBLIC :: ad_variability

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_variability (ng, tile, Linp, Lweak)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_ocean
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak

      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_variability_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Linp, Lweak,                            &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          BOUNDARY(ng) % e_t_obc,                 &
     &                          BOUNDARY(ng) % e_u_obc,                 &
     &                          BOUNDARY(ng) % e_v_obc,                 &
#  endif
     &                          BOUNDARY(ng) % e_ubar_obc,              &
     &                          BOUNDARY(ng) % e_vbar_obc,              &
     &                          BOUNDARY(ng) % e_zeta_obc,              &
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % e_sustr,                   &
     &                          FORCES(ng) % e_svstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % e_stflx,                   &
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % e_t,                        &
     &                          OCEAN(ng) % e_u,                        &
     &                          OCEAN(ng) % e_v,                        &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                          OCEAN(ng) % e_ubar,                     &
     &                          OCEAN(ng) % e_vbar,                     &
#  endif
# else
     &                          OCEAN(ng) % e_ubar,                     &
     &                          OCEAN(ng) % e_vbar,                     &
# endif
     &                          OCEAN(ng) % e_zeta,                     &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          BOUNDARY(ng) % ad_t_obc,                &
     &                          BOUNDARY(ng) % ad_u_obc,                &
     &                          BOUNDARY(ng) % ad_v_obc,                &
#  endif
     &                          BOUNDARY(ng) % ad_ubar_obc,             &
     &                          BOUNDARY(ng) % ad_vbar_obc,             &
     &                          BOUNDARY(ng) % ad_zeta_obc,             &
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % ad_ustr,                   &
     &                          FORCES(ng) % ad_vstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % ad_tflux,                  &
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
#  endif
# else
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
# endif
     &                          OCEAN(ng) % ad_zeta)

      RETURN
      END SUBROUTINE ad_variability
!
!***********************************************************************
      SUBROUTINE ad_variability_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, LBij, UBij,   &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                Linp, Lweak,                      &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                                t_obc_std, u_obc_std, v_obc_std,  &
#  endif
     &                                ubar_obc_std, vbar_obc_std,       &
     &                                zeta_obc_std,                     &
# endif
# ifdef ADJUST_WSTRESS
     &                                sustr_std, svstr_std,             &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                stflx_std,                        &
# endif
# ifdef SOLVE3D
     &                                t_std, u_std, v_std,              &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                                ubar_std, vbar_std,               &
#  endif
# else
     &                                ubar_std, vbar_std,               &
# endif
     &                                zeta_std,                         &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                                ad_t_obc, ad_u_obc, ad_v_obc,     &
#  endif
     &                                ad_ubar_obc, ad_vbar_obc,         &
     &                                ad_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                                ad_ustr, ad_vstr,                 &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                ad_tflux,                         &
# endif
# ifdef SOLVE3D
     &                                ad_t, ad_u, ad_v,                 &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                                ad_ubar, ad_vbar,                 &
#  endif
# else
     &                                ad_ubar, ad_vbar,                 &
# endif
     &                                ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lweak

      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: t_obc_std(LBij:,:,:,:)
      real(r8), intent(in) :: u_obc_std(LBij:,:,:)
      real(r8), intent(in) :: v_obc_std(LBij:,:,:)
#   endif
      real(r8), intent(in) :: ubar_obc_std(LBij:,:)
      real(r8), intent(in) :: vbar_obc_std(LBij:,:)
      real(r8), intent(in) :: zeta_obc_std(LBij:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx_std(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u_std(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v_std(LBi:,LBj:,:,:)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: ubar_std(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar_std(LBi:,LBj:,:)
#   endif
#  else
      real(r8), intent(in) :: ubar_std(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar_std(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: zeta_std(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)

# else

#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(in) :: t_obc_std(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(in) :: u_obc_std(LBij:UBij,N(ng),4)
      real(r8), intent(in) :: v_obc_std(LBij:UBij,N(ng),4)
#   endif
      real(r8), intent(in) :: ubar_obc_std(LBij:UBij,4)
      real(r8), intent(in) :: vbar_obc_std(LBij:UBij,4)
      real(r8), intent(in) :: zeta_obc_std(LBij:UBij,4)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr_std(LBi:,LBj:)
      real(r8), intent(in) :: svstr_std(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx_std(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: t_std(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(in) :: u_std(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(in) :: v_std(LBi:UBi,LBj:UBj,N(ng),NSA)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(in) :: ubar_std(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: vbar_std(LBi:UBi,LBj:UBj,NSA)
#   endif
#  else
      real(r8), intent(in) :: ubar_std(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: vbar_std(LBi:UBi,LBj:UBj,NSA)
#  endif
      real(r8), intent(in) :: zeta_std(LBi:UBi,LBj:UBj,NSA)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, ib, ir, j, rec
# ifdef SOLVE3D
      integer :: it, k
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Determine error covariance standard deviation factors to use.
!-----------------------------------------------------------------------
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint
      ELSE
        rec=1                        ! strong constraint
      END IF
!
!-----------------------------------------------------------------------
!  Initial conditions and model error covariance: Multiply adjoint state
!  by its corresponding  standard deviation.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
# ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,Linp))
# endif
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)*zeta_std(i,j,rec)
        END DO
      END DO
# if !defined SOLVE3D || (defined WEAK_CONSTRAINT && defined TIME_CONV)
!
!  2D U-momentum.
!
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ubar(:,:,Linp),                         &
     &                       ad_vbar(:,:,Linp))
#  endif
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar(i,j,Linp)=ad_ubar(i,j,Linp)*ubar_std(i,j,rec)
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          ad_vbar(i,j,Linp)=ad_vbar(i,j,Linp)*vbar_std(i,j,rec)
        END DO
      END DO
# endif
# ifdef SOLVE3D
!
!  3D U-momentum.
!
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_u(:,:,:,Linp),                          &
     &                       ad_v(:,:,:,Linp))
#  endif
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)*u_std(i,j,k,rec)
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)*v_std(i,j,k,rec)
          END DO
        END DO
      END DO
!
!  Tracers.
!
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_t(:,:,:,Linp,:))
#  endif
      DO it=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              ad_t(i,j,k,Linp,it)=ad_t(i,j,k,Linp,it)*                  &
     &                            t_std(i,j,k,rec,it)
            END DO
          END DO
        END DO
      END DO
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Boundary conditions error covariance: Multiply adjoint state by its
!  corresponding  standard deviation.
!-----------------------------------------------------------------------
!
!  Free-surface open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isFsur,ng)))) THEN
        DO ir=1,Nbrec(ng)
#  ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isFsur,ng)) THEN
              CALL ad_mp_exchange2d_bry (ng, tile, iADM, 1, ib,         &
     &                                   LBij, UBij,                    &
     &                                   NghostPoints,                  &
     &                                   EWperiodic(ng), NSperiodic(ng),&
     &                                   ad_zeta_obc(:,ib,ir,Linp))
            END IF
          END DO
#  endif
          IF ((Lobc(iwest,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              ad_zeta_obc(j,ib,ir,Linp)=ad_zeta_obc(j,ib,ir,Linp)*      &
     &                                  zeta_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(ieast,isFsur,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              ad_zeta_obc(j,ib,ir,Linp)=ad_zeta_obc(j,ib,ir,Linp)*      &
     &                                  zeta_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(isouth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              ad_zeta_obc(i,ib,ir,Linp)=ad_zeta_obc(i,ib,ir,Linp)*      &
     &                                  zeta_obc_std(i,ib)
            END DO
          END IF
          IF ((Lobc(inorth,isFsur,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              ad_zeta_obc(i,ib,ir,Linp)=ad_zeta_obc(i,ib,ir,Linp)*      &
     &                                  zeta_obc_std(i,ib)
            END DO
          END IF
        END DO
      END IF
!
!  2D U-momentum open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isUbar,ng)))) THEN
        DO ir=1,Nbrec(ng)
#  ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isUbar,ng)) THEN
              CALL ad_mp_exchange2d_bry (ng, tile, iADM, 1, ib,         &
     &                                   LBij, UBij,                    &
     &                                   NghostPoints,                  &
     &                                   EWperiodic(ng), NSperiodic(ng),&
     &                                   ad_ubar_obc(:,ib,ir,Linp))
            END IF
          END DO
#  endif
          IF ((Lobc(iwest,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=Jstr,Jend
              ad_ubar_obc(j,ib,ir,Linp)=ad_ubar_obc(j,ib,ir,Linp)*      &
     &                                  ubar_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(ieast,isUbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=Jstr,Jend
              ad_ubar_obc(j,ib,ir,Linp)=ad_ubar_obc(j,ib,ir,Linp)*      &
     &                                  ubar_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(isouth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=IstrU,Iend
              ad_ubar_obc(i,ib,ir,Linp)=ad_ubar_obc(i,ib,ir,Linp)*      &
     &                                  ubar_obc_std(i,ib)
            END DO
          END IF
          IF ((Lobc(inorth,isUbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=IstrU,Iend
              ad_ubar_obc(i,ib,ir,Linp)=ad_ubar_obc(i,ib,ir,Linp)*      &
     &                                  ubar_obc_std(i,ib)
            END DO
          END IF
        END DO
      END IF
!
!  2D V-momentum open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isVbar,ng)))) THEN
        DO ir=1,Nbrec(ng)
#  ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isVbar,ng)) THEN
              CALL ad_mp_exchange2d_bry (ng, tile, iADM, 1, ib,         &
     &                                   LBij, UBij,                    &
     &                                   NghostPoints,                  &
     &                                   EWperiodic(ng), NSperiodic(ng),&
     &                                   ad_vbar_obc(:,ib,ir,Linp))
            END IF
          END DO
#  endif
          IF ((Lobc(iwest,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO j=JstrV,Jend
              ad_vbar_obc(j,ib,ir,Linp)=ad_vbar_obc(j,ib,ir,Linp)*      &
     &                                  vbar_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(ieast,isVbar,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO j=JstrV,Jend
              ad_vbar_obc(j,ib,ir,Linp)=ad_vbar_obc(j,ib,ir,Linp)*      &
     &                                  vbar_obc_std(j,ib)
            END DO
          END IF
          IF ((Lobc(isouth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO i=Istr,Iend
              ad_vbar_obc(i,ib,ir,Linp)=ad_vbar_obc(i,ib,ir,Linp)*      &
     &                                  vbar_obc_std(i,ib)
            END DO
          END IF
          IF ((Lobc(inorth,isVbar,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO i=Istr,Iend
              ad_vbar_obc(i,ib,ir,Linp)=ad_vbar_obc(i,ib,ir,Linp)*      &
     &                                  vbar_obc_std(i,ib)
            END DO
          END IF
        END DO
      END IF

#  ifdef SOLVE3D
!
!  3D U-momentum open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isUvel,ng)))) THEN
        DO ir=1,Nbrec(ng)
#   ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isUvel,ng)) THEN
              CALL ad_mp_exchange3d_bry (ng, tile, iADM, 1, ib,         &
     &                                   LBij, UBij, 1, N(ng),          &
     &                                   NghostPoints,                  &
     &                                   EWperiodic(ng), NSperiodic(ng),&
     &                                   ad_u_obc(:,:,ib,ir,Linp))
            END IF
          END DO
#   endif
          IF ((Lobc(iwest,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=Jstr,Jend
                ad_u_obc(j,k,ib,ir,Linp)=ad_u_obc(j,k,ib,ir,Linp)*      &
     &                                   u_obc_std(j,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isUvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=Jstr,Jend
                ad_u_obc(j,k,ib,ir,Linp)=ad_u_obc(j,k,ib,ir,Linp)*      &
     &                                   u_obc_std(j,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                ad_u_obc(i,k,ib,ir,Linp)=ad_u_obc(i,k,ib,ir,Linp)*      &
     &                                   u_obc_std(i,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isUvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=IstrU,Iend
                ad_u_obc(i,k,ib,ir,Linp)=ad_u_obc(i,k,ib,ir,Linp)*      &
     &                                   u_obc_std(i,k,ib)
              END DO
            END DO
          END IF
        END DO
      END IF
!
!  3D V-momentum open boundaries.
!
      IF (.not.Lweak.and.(ANY(Lobc(:,isVvel,ng)))) THEN
        DO ir=1,Nbrec(ng)
#   ifdef DISTRIBUTE
          DO ib=1,4
            IF (Lobc(ib,isVvel,ng)) THEN
              CALL ad_mp_exchange3d_bry (ng, tile, iADM, 1, ib,         &
     &                                   LBij, UBij, 1, N(ng),          &
     &                                   NghostPoints,                  &
     &                                   EWperiodic(ng), NSperiodic(ng),&
     &                                   ad_v_obc(:,:,ib,ir,Linp))
            END IF
          END DO
#   endif
          IF ((Lobc(iwest,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Western_Edge(tile)) THEN
            ib=iwest
            DO k=1,N(ng)
              DO j=JstrV,Jend
                ad_v_obc(j,k,ib,ir,Linp)=ad_v_obc(j,k,ib,ir,Linp)*      &
     &                                   v_obc_std(j,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(ieast,isVvel,ng)).and.                              &
     &        DOMAIN(ng)%Eastern_Edge(tile)) THEN
            ib=ieast
            DO k=1,N(ng)
              DO j=JstrV,Jend
                ad_v_obc(j,k,ib,ir,Linp)=ad_v_obc(j,k,ib,ir,Linp)*      &
     &                                   v_obc_std(j,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(isouth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Southern_Edge(tile)) THEN
            ib=isouth
            DO k=1,N(ng)
              DO i=Istr,Iend
                ad_v_obc(i,k,ib,ir,Linp)=ad_v_obc(i,k,ib,ir,Linp)*      &
     &                                   v_obc_std(i,k,ib)
              END DO
            END DO
          END IF
          IF ((Lobc(inorth,isVvel,ng)).and.                             &
     &        DOMAIN(ng)%Northern_Edge(tile)) THEN
            ib=inorth
            DO k=1,N(ng)
              DO i=Istr,Iend
                ad_v_obc(i,k,ib,ir,Linp)=ad_v_obc(i,k,ib,ir,Linp)*      &
     &                                   v_obc_std(i,k,ib)
              END DO
            END DO
          END IF
        END DO
      END IF
!
!  Tracers open boundaries.
!
      DO it=1,NT(ng)
        IF (.not.Lweak.and.(ANY(Lobc(:,isTvar(it),ng)))) THEN
          DO ir=1,Nbrec(ng)
#   ifdef DISTRIBUTE
            DO ib=1,4
              IF (Lobc(ib,isTvar(it),ng)) THEN
                CALL ad_mp_exchange3d_bry (ng, tile, iADM, 1, ib,       &
     &                                     LBij, UBij, 1, N(ng),        &
     &                                     NghostPoints,                &
     &                                     EWperiodic(ng),              &
     &                                     NSperiodic(ng),              &
     &                                     ad_t_obc(:,:,ib,ir,Linp,it))
              END IF
            END DO
#   endif
            IF ((Lobc(iwest,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Western_Edge(tile)) THEN
              ib=iwest
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  ad_t_obc(j,k,ib,ir,Linp,it)=                          &
     &                                    ad_t_obc(j,k,ib,ir,Linp,it)*  &
     &                                    t_obc_std(j,k,ib,it)
                END DO
              END DO
            END IF
            IF ((Lobc(ieast,isTvar(it),ng)).and.                        &
     &          DOMAIN(ng)%Eastern_Edge(tile)) THEN
              ib=ieast
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  ad_t_obc(j,k,ib,ir,Linp,it)=                          &
     &                                    ad_t_obc(j,k,ib,ir,Linp,it)*  &
     &                                    t_obc_std(j,k,ib,it)
                END DO
              END DO
            END IF
            IF ((Lobc(isouth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Southern_Edge(tile)) THEN
              ib=isouth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  ad_t_obc(i,k,ib,ir,Linp,it)=                          &
     &                                    ad_t_obc(i,k,ib,ir,Linp,it)*  &
     &                                    t_obc_std(i,k,ib,it)
                END DO
              END DO
            END IF
            IF ((Lobc(inorth,isTvar(it),ng)).and.                       &
     &          DOMAIN(ng)%Northern_Edge(tile)) THEN
              ib=inorth
              DO k=1,N(ng)
                DO i=Istr,Iend
                  ad_t_obc(i,k,ib,ir,Linp,it)=                          &
     &                                    ad_t_obc(i,k,ib,ir,Linp,it)*  &
     &                                    t_obc_std(i,k,ib,it)
                END DO
              END DO
            END IF
          END DO
        END IF
      END DO
#  endif
# endif

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!-----------------------------------------------------------------------
!  Surface forcing error covariance: Multiply adjoint state by its
!  corresponding  standard deviation.
!-----------------------------------------------------------------------

#  ifdef ADJUST_WSTRESS
!
!  Adjoint surface momentum stress.
!
      IF (.not.Lweak) THEN
#   ifdef DISTRIBUTE
        CALL ad_mp_exchange3d (ng, tile, iADM, 2,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng),          &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ustr(:,:,:,Linp),                       &
     &                       ad_vstr(:,:,:,Linp))
#   endif
        DO ir=1,Nfrec(ng)
          DO j=JstrT,JendT
            DO i=IstrP,IendT
              ad_ustr(i,j,ir,Linp)=ad_ustr(i,j,ir,Linp)*sustr_std(i,j)
            END DO
          END DO
          DO j=JstrP,JendT
            DO i=IstrT,IendT
              ad_vstr(i,j,ir,Linp)=ad_vstr(i,j,ir,Linp)*svstr_std(i,j)
            END DO
          END DO
        END DO
      END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Surface tracer flux.
!
      IF (.not.Lweak) THEN
        DO it=1,NT(ng)
          IF (Lstflux(it,ng)) THEN
#   ifdef DISTRIBUTE
            CALL ad_mp_exchange3d (ng, tile, iADM, 1,                   &
     &                           LBi, UBi, LBj, UBj, 1,Nfrec(ng),       &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_tflux(:,:,:,Linp,it))
#   endif
            DO ir=1,Nfrec(ng)
              DO j=JstrT,JendT
                DO i=IstrT,IendT
                  ad_tflux(i,j,ir,Linp,it)=ad_tflux(i,j,ir,Linp,it)*    &
     &                                   stflx_std(i,j,it)
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
# endif

      RETURN
      END SUBROUTINE ad_variability_tile
#endif
      END MODULE ad_variability_mod
