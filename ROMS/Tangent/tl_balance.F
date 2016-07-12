#include "cppdefs.h"
      MODULE tl_balance_mod

#if defined TANGENT && defined BALANCE_OPERATOR
!
!svn $Id: tl_balance.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines impose a multivariate balance operator to constraint !
!  the  background/model  error  covariance matrix,  B,  such that the !
!  unobserved variables information is extracted from  observed  data. !
!  It follows the approach proposed by Weaver et al. (2006). The state !
!  vector is split between balanced and unbalanced components,  except !
!  for temperature,  which is used to  establish the  balanced part of !
!  other state variables.                                              !
!                                                                      !
!  The background/model error covariance is represented as:            !
!                                                                      !
!     B = K Bu K^(T)                                                   !
!                                                                      !
!  where                                                               !
!                                                                      !
!     B : background/model error covariance matrix.                    !
!     Bu: unbalanced background/model error covariance matrix modeled  !
!         with the generalized diffusion operator.                     !
!     K : balance matrix operator.                                     !
!                                                                      !
!  Here, T denotes the transpose.                                      !
!                                                                      !
!  The multivariate formulation is obtained by  establishing  balance  !
!  relationships with the other state variables  using  T-S empirical  !
!  formulas, hydrostatic balance, and geostrophic balance.             !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Weaver, A.T., C. Deltel, E. Machu, S. Ricci, and N. Daget, 2006:  !
!      A multivariate balance operator for variational data assimila-  !
!      tion, Q. J. R. Meteorol. Soc., 131, 3605-3625.                  !
!                                                                      !
!      (See also, ECMWR Technical Memorandum # 491, April 2006)        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      PRIVATE
      PUBLIC :: tl_balance
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_balance (ng, tile, Lbck, Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
      USE mod_mixing
# endif
# ifdef ZETA_ELLIPTIC
      USE mod_fourdvar
# endif
      USE mod_ocean
      USE mod_stepping
# ifdef SOLVE3D
!
      USE rho_eos_mod
      USE set_depth_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lbck, Linp
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef SOLVE3D
!
!  Compute background state thickness, depth arrays, thermal expansion,
!  and saline contraction coefficients.
!
      COUPLING(ng) % Zt_avg1 = 0.0_r8

      CALL set_depth (ng, tile)
      nrhs(ng)=Lbck
      CALL rho_eos (ng, tile)
!
# endif
      CALL tl_balance_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Lbck, Linp,                                 &
     &                      GRID(ng) % f,                               &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef ZETA_ELLIPTIC
     &                      GRID(ng) % pmon_u,                          &
     &                      GRID(ng) % pnom_v,                          &
     &                      GRID(ng) % h,                               &
# endif
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef SOLVE3D
     &                      MIXING(ng) % alpha,                         &
     &                      MIXING(ng) % beta,                          &
     &                      OCEAN(ng) % t,                              &
# endif
# ifdef ZETA_ELLIPTIC
     &                      FOURDVAR(ng) % tl_zeta_ref,                 &
     &                      FOURDVAR(ng) % tl_rhs_r2d,                  &
     &                      FOURDVAR(ng) % pc_r2d,                      &
     &                      FOURDVAR(ng) % r_r2d,                       &
     &                      FOURDVAR(ng) % br_r2d,                      &
     &                      FOURDVAR(ng) % p_r2d,                       &
     &                      FOURDVAR(ng) % bp_r2d,                      &
     &                      FOURDVAR(ng) % zdf1,                        &
     &                      FOURDVAR(ng) % zdf2,                        &
     &                      FOURDVAR(ng) % zdf3,                        &
     &                      FOURDVAR(ng) % bc_ak,                       &
     &                      FOURDVAR(ng) % bc_bk,                       &
# endif
# ifdef SOLVE3D
     &                      OCEAN(ng) % tl_rho,                         &
     &                      OCEAN(ng) % tl_t,                           &
     &                      OCEAN(ng) % tl_u,                           &
     &                      OCEAN(ng) % tl_v,                           &
# endif
     &                      OCEAN(ng) % tl_zeta)

      RETURN
      END SUBROUTINE tl_balance
!
!***********************************************************************
      SUBROUTINE tl_balance_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, LBij, UBij,       &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Lbck, Linp,                           &
     &                            f, pm, pn,                            &
# ifdef ZETA_ELLIPTIC
     &                            pmon_u, pnom_v, h,                    &
# endif
# ifdef SOLVE3D
     &                            Hz, z_r, z_w,                         &
# endif
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# ifdef SOLVE3D
     &                            alpha, beta, t,                       &
# endif
# ifdef ZETA_ELLIPTIC
     &                            tl_zeta_ref, tl_rhs_r2d,              &
     &                            pc_r2d, r_r2d, br_r2d,                &
     &                            p_r2d, bp_r2d,                        &
     &                            zdf1, zdf2, zdf3, bc_ak, bc_bk,       &
# endif
# ifdef SOLVE3D
     &                            tl_rho, tl_t, tl_u, tl_v,             &
# endif
     &                            tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE bc_2d_mod
      USE exchange_2d_mod
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
# ifdef ZETA_ELLIPTIC
      USE zeta_balance_mod, ONLY: r2d_bc, u2d_bc, v2d_bc
      USE zeta_balance_mod, ONLY: tl_biconj_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck, Linp
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  ifdef ZETA_ELLIPTIC
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: alpha(LBi:,LBj:)
      real(r8), intent(in) :: beta(LBi:,LBj:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(out) :: tl_rho(LBi:,LBj:,:)
#  endif
#  ifdef ZETA_ELLIPTIC
      real(r8), intent(in) :: bc_ak(:)
      real(r8), intent(in) :: bc_bk(:)
      real(r8), intent(in) :: zdf1(:)
      real(r8), intent(in) :: zdf2(:)
      real(r8), intent(in) :: zdf3(:)
      real(r8), intent(inout) :: pc_r2d(LBi:,LBj:)
      real(r8), intent(inout) :: r_r2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: br_r2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: p_r2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: bp_r2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_rhs_r2d(LBi:,LBj:)
      real(r8), intent(inout) :: tl_zeta_ref(LBi:,LBj:)
#  endif

# else

      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  ifdef ZETA_ELLIPTIC
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: beta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,2,N(ng))
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,2,N(ng))
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(out) :: tl_rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
#  ifdef ZETA_ELLIPTIC
      real(r8), intent(in) :: bc_ak(Nbico(ng))
      real(r8), intent(in) :: bc_bk(Nbico(ng))
      real(r8), intent(in) :: zdf1(Nbico(ng))
      real(r8), intent(in) :: zdf2(Nbico(ng))
      real(r8), intent(in) :: zdf3(Nbico(ng))
      real(r8), intent(inout) :: pc_r2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: r_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: br_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: p_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: bp_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: tl_rhs_r2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_zeta_ref(LBi:UBi,LBj:UBj)
#  endif

# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, order
      integer :: Norder = 2                 ! Shapiro filter order

      real(r8) :: fac, fac1, fac2, fac3, gamma
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: tl_cff, tl_cff1, tl_cff2
      real(r8) :: dzdT, zphi, zphi1, zwbot, zwtop

      real(r8), dimension(20) ::  filter_coef =                         &
     &   (/ 2.500000E-1_r8,    6.250000E-2_r8,     1.562500E-2_r8,      &
     &      3.906250E-3_r8,    9.765625E-4_r8,     2.44140625E-4_r8,    &
     &      6.103515625E-5_r8, 1.5258789063E-5_r8, 3.814697E-6_r8,      &
     &      9.536743E-7_r8,    2.384186E-7_r8,     5.960464E-8_r8,      &
     &      1.490116E-8_r8,    3.725290E-9_r8,     9.313226E-10_r8,     &
     &      2.328306E-10_r8,   5.820766E-11_r8,    1.455192E-11_r8,     &
     &      3.637979E-12_r8,   9.094947E-13_r8 /)

      real(r8), dimension(N(ng)) :: dSdT, dSdT_filter

      real(r8), dimension(IminS:ImaxS) :: tl_phie
      real(r8), dimension(IminS:ImaxS) :: tl_phix

# ifdef SALINITY
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dTdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dSdz
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: tl_gradP

# ifdef ZETA_ELLIPTIC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_phi

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_gradPx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_gradPy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_r2d_rhs
# endif

# include "set_bounds.h"

# ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Compute balance salinity contribution.
!-----------------------------------------------------------------------
!
      IF (balance(isTvar(isalt))) THEN
!
!  Compute temperature (dTdz) and salinity (dSdz) shears.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            FC(i,0)=0.0_r8
            dTdz(i,j,0)=0.0_r8
            dSdz(i,j,0)=0.0_r8
          END DO
          DO k=1,N(ng)-1
            DO i=IstrR,IendR
              cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                           &
     &            Hz(i,j,k)*(2.0_r8-FC(i,k-1)))
              FC(i,k)=cff*Hz(i,j,k+1)
              dTdz(i,j,k)=cff*(6.0_r8*(t(i,j,k+1,Lbck,itemp)-           &
     &                                 t(i,j,k  ,Lbck,itemp))-          &
     &                         Hz(i,j,k)*dTdz(i,j,k-1))
              dSdz(i,j,k)=cff*(6.0_r8*(t(i,j,k+1,Lbck,isalt)-           &
     &                                 t(i,j,k  ,Lbck,isalt))-          &
     &                         Hz(i,j,k)*dSdz(i,j,k-1))
            END DO
          END DO
          DO i=IstrR,IendR
            dTdz(i,j,N(ng))=0.0_r8
            dSdz(i,j,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=IstrR,IendR
              dTdz(i,j,k)=dTdz(i,j,k)-FC(i,k)*dTdz(i,j,k+1)
              dSdz(i,j,k)=dSdz(i,j,k)-FC(i,k)*dSdz(i,j,k+1)
            END DO
          END DO
        END DO
!
!  Add balanced salinity (deltaS_b) contribution to unbalanced salinity
!  increment. The unbalanced salinity increment is related related to
!  temperature increment:
!
!       deltaS_b = cff * dS/dz * dz/dT * deltaT
!
!  Here, cff is a coefficient that depends on the mixed layer depth:
!
!       cff = 1.0 - EXP (z_r / ml_depth)
!
!  the coefficient is smoothly reduced to zero at the surface and below
!  the mixed layer.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            DO k=1,N(ng)
              cff=0.5_r8*(dTdz(i,j,k-1)+dTdz(i,j,k))
              IF (ABS(cff).lt.dTdz_min(ng)) THEN
                dzdT=0.0_r8
              ELSE
                dzdT=1.0_r8/cff
              END IF
              dSdT(k)=(0.5_r8*(dSdz(i,j,k-1)+                           &
     &                         dSdz(i,j,k  )))*dzdT
            END DO
!
!  Shapiro filter.
!
            DO order=1,Norder/2
              IF (order.ne.Norder/2) THEN
                dSdT_filter(1)=2.0_r8*(dSdT(1)-dSdT(2))
                dSdT_filter(N(ng))=2.0_r8*(dSdT(N(ng))-dSdT(N(ng)-1))
              ELSE
                dSdT_filter(1)=0.0_r8
                dSdT_filter(N(ng))=0.0_r8
              END IF
              DO k=2,N(ng)-1
                dSdT_filter(k)=2.0_r8*dSdT(k)-dSdT(k-1)-dSdT(k+1)
              END DO
              DO k=1,N(ng)
                dSdT(k)=dSdT(k)-filter_coef(Norder/2)*dSdT_filter(k)
              END DO
            END DO

            DO k=1,N(ng)
              cff=(1.0_r8-EXP(z_r(i,j,k)/ml_depth(ng)))*dSdT(k)
              tl_t(i,j,k,Linp,isalt)=tl_t(i,j,k,Linp,isalt)+            &
     &                               cff*tl_t(i,j,k,Linp,itemp)
#  ifdef MASKING
              tl_t(i,j,k,Linp,isalt)=tl_t(i,j,k,Linp,isalt)*rmask(i,j)
#  endif
            END DO
          END DO
        END DO

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            tl_t(:,:,:,Linp,isalt))
        END IF
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_t(:,:,:,Linp,isalt))
#  endif
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Compute balanced density anomaly increment using linearized equation
!  of state.  The thermal expansion and saline contraction coefficients
!  are computed from the background state.
!-----------------------------------------------------------------------
!
      IF (balance(isFsur).or.balance(isVvel)) THEN
        DO j=JstrR,JendR
          DO k=1,N(ng)
            DO i=IstrR,IendR
              tl_rho(i,j,k)=-rho0*alpha(i,j)*tl_t(i,j,k,Linp,itemp)
# ifdef SALINITY
              tl_rho(i,j,k)=tl_rho(i,j,k)+                              &
     &                      rho0*beta(i,j)*tl_t(i,j,k,Linp,isalt)
# endif
# ifdef MASKING
              tl_rho(i,j,k)=tl_rho(i,j,k)*rmask(i,j)
# endif
            END DO
          END DO
        END DO

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            tl_rho)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_rho)
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Add balanced velocity contributions to unbalanced velocity
!  increments. Use linear pressure gradient formulation based
!  to that found in routine "prsgrd31.h".
!-----------------------------------------------------------------------
!
!  NOTE: fac2=0 because the balanced component should consist of the
!  baroclinic pressure gradient only.
!
      fac1=0.5_r8*g/rho0
!!    fac2=g
      fac2=0.0_r8
      fac3=0.25_r8*g/rho0

# ifdef ZETA_ELLIPTIC
!
!  Initialize vertical intergal local variables.
!
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          tl_gradPx(i,j)=0.0_r8
          tl_gradPy(i,j)=0.0_r8
        END DO
      END DO
# endif
!
!  Compute balanced, surface U-momentum from baroclinic and barotropic
!  surface pressure gradient.
!
      IF (balance(isFsur).or.balance(isUvel)) THEN
        DO j=Jstr,Jend+1
          DO i=Istr-1,Iend
            cff1=z_w(i,j  ,N(ng))-z_r(i,j  ,N(ng))+                     &
     &           z_w(i,j-1,N(ng))-z_r(i,j-1,N(ng))
            tl_phie(i)=fac1*(tl_rho(i,j  ,N(ng))-                       &
     &                       tl_rho(i,j-1,N(ng)))*cff1+                 &
     &                 fac2*(tl_zeta(i,j  ,Linp)-                       &
     &                       tl_zeta(i,j-1,Linp))
            tl_gradP(i,j,N(ng))=0.5_r8*tl_phie(i)*                      &
     &                          (pn(i,j-1)+pn(i,j))/(f(i,j-1)+f(i,j))
# ifdef ZETA_ELLIPTIC
            tl_phi(i,N(ng))=tl_phie(i)
# endif
          END DO
!
!  Compute balance, interior U-momentum from baroclinic pressure
!  gradient (differentiate and then vertically integrate).
!
          DO k=N(ng)-1,1,-1
            DO i=Istr-1,Iend
              cff1=1.0_r8/((z_r(i,j  ,k+1)-z_r(i,j  ,k))*               &
     &                     (z_r(i,j-1,k+1)-z_r(i,j-1,k)))
              cff2=z_r(i,j  ,k  )-z_r(i,j-1,k  )+                       &
     &             z_r(i,j  ,k+1)-z_r(i,j-1,k+1)
              cff3=z_r(i,j  ,k+1)-z_r(i,j  ,k  )-                       &
     &             z_r(i,j-1,k+1)+z_r(i,j-1,k  )
              gamma=0.125_r8*cff1*cff2*cff3
!
              tl_cff1=(1.0_r8+gamma)*(tl_rho(i,j  ,k+1)-                &
     &                                tl_rho(i,j-1,k+1))+               &
     &                (1.0_r8-gamma)*(tl_rho(i,j  ,k  )-                &
     &                                tl_rho(i,j-1,k  ))
              tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-                &
     &                tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
              cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i,j-1,k+1))+        &
     &             (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i,j-1,k  ))
              tl_phie(i)=tl_phie(i)+                                    &
     &                   fac3*(tl_cff1*cff3-tl_cff2*cff4)
              tl_gradP(i,j,k)=0.5_r8*tl_phie(i)*                        &
     &                        (pn(i,j-1)+pn(i,j))/(f(i,j-1)+f(i,j))
# ifdef ZETA_ELLIPTIC
              tl_phi(i,k)=tl_phie(i)
# endif
            END DO
          END DO

# ifdef ZETA_ELLIPTIC
!
!  Compute right-hand-side term used in the elliptic equation for
!  unbalanced free-surface.  Integrate from bottom to surface.
!
          IF (balance(isFsur)) THEN
            DO k=1,N(ng)
              DO i=Istr-1,Iend
                tl_cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))*tl_phi(i,k)
#  ifdef MASKING
                tl_cff=tl_cff*vmask(i,j)
#  endif
                tl_gradPy(i,j)=tl_gradPy(i,j)+tl_cff
              END DO
            END DO
          END IF
# endif
        END DO
      END IF
!
!  Compute 3D U-velocity error covariance constraint.
!
      IF (balance(isUvel)) THEN
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              tl_u(i,j,k,Linp)=tl_u(i,j,k,Linp)-                        &
     &                         0.25_r8*(tl_gradP(i-1,j  ,k)+            &
     &                                  tl_gradP(i  ,j  ,k)+            &
     &                                  tl_gradP(i-1,j+1,k)+            &
     &                                  tl_gradP(i  ,j+1,k))
# ifdef MASKING
              tl_u(i,j,k,Linp)=tl_u(i,j,k,Linp)*umask(i,j)
# endif
            END DO
          END DO
        END DO
      END IF
!
!  Compute balanced, surface V-momentum from baroclinic and barotropic
!  surface pressure gradient.
!
      IF (balance(isFsur).or.balance(isVvel)) THEN
        DO j=Jstr-1,Jend
          DO i=Istr,Iend+1
            cff1=z_w(i  ,j,N(ng))-z_r(i  ,j,N(ng))+                     &
     &           z_w(i-1,j,N(ng))-z_r(i-1,j,N(ng))
            tl_phix(i)=fac1*(tl_rho(i  ,j,N(ng))-                       &
     &                       tl_rho(i-1,j,N(ng)))*cff1+                 &
     &                 fac2*(tl_zeta(i  ,j,Linp)-                       &
     &                       tl_zeta(i-1,j,Linp))
            tl_gradP(i,j,N(ng))=0.5_r8*tl_phix(i)*                      &
     &                          (pm(i-1,j)+pm(i,j))/(f(i-1,j)+f(i,j))
# ifdef ZETA_ELLIPTIC
            tl_phi(i,N(ng))=tl_phix(i)
# endif
          END DO
!
!  Compute balance, interior V-momentum from baroclinic pressure
!  gradient (differentiate and then vertically integrate).
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend+1
              cff1=1.0_r8/((z_r(i  ,j,k+1)-z_r(i  ,j,k))*               &
     &                     (z_r(i-1,j,k+1)-z_r(i-1,j,k)))
              cff2=z_r(i  ,j,k  )-z_r(i-1,j,k  )+                       &
     &             z_r(i  ,j,k+1)-z_r(i-1,j,k+1)
              cff3=z_r(i  ,j,k+1)-z_r(i  ,j,k  )-                       &
     &             z_r(i-1,j,k+1)+z_r(i-1,j,k  )
              gamma=0.125_r8*cff1*cff2*cff3
!
              tl_cff1=(1.0_r8+gamma)*(tl_rho(i  ,j,k+1)-                &
     &                                tl_rho(i-1,j,k+1))+               &
     &                (1.0_r8-gamma)*(tl_rho(i  ,j,k  )-                &
     &                                tl_rho(i-1,j,k  ))
              tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-                &
     &                tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
              cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                         &
     &             z_r(i,j,k  )-z_r(i-1,j,k  )
              cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i-1,j,k+1))+        &
     &             (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i-1,j,k  ))
              tl_phix(i)=tl_phix(i)+                                    &
     &                   fac3*(tl_cff1*cff3-tl_cff2*cff4)
              tl_gradP(i,j,k)=0.5_r8*tl_phix(i)*                        &
     &                        (pm(i-1,j)+pm(i,j))/(f(i-1,j)+f(i,j))
# ifdef ZETA_ELLIPTIC
              tl_phi(i,k)=tl_phix(i)
# endif
            END DO
          END DO

# ifdef ZETA_ELLIPTIC
!
!  Compute right-hand-side term used in the elliptic equation for
!  unbalanced free-surface.  Integrate from bottom to surface.
!
          IF (balance(isFsur)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend+1
                tl_cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))*tl_phi(i,k)
#  ifdef MASKING
                tl_cff=tl_cff*umask(i,j)
#  endif
                tl_gradPx(i,j)=tl_gradPx(i,j)+tl_cff
              END DO
            END DO
          END IF
# endif
        END DO
      END IF
!
!  Compute 3D V-velocity error covariance constraint.
!
      IF (balance(isVvel)) THEN
        DO k=1,N(ng)
          DO j=JstrV,Jend
            DO i=Istr,Iend
              tl_v(i,j,k,Linp)=tl_v(i,j,k,Linp)+                        &
     &                         0.25_r8*(tl_gradP(i  ,j-1,k)+            &
     &                                  tl_gradP(i+1,j-1,k)+            &
     &                                  tl_gradP(i  ,j  ,k)+            &
     &                                  tl_gradP(i+1,j  ,k))
# ifdef MASKING
              tl_v(i,j,k,Linp)=tl_v(i,j,k,Linp)*vmask(i,j)
# endif
            END DO
          END DO
        END DO

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_u3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            tl_u(:,:,:,Linp))
          CALL exchange_v3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            tl_v(:,:,:,Linp))
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iTLM, 2,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_u(:,:,:,Linp), tl_v(:,:,:,Linp))
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Add balanced free-surface contribution to unbalanced free-surface
!  increment.
!-----------------------------------------------------------------------
!
      IF (balance(isFsur)) THEN

# ifdef ZETA_ELLIPTIC
!
!  Solve elliptic equation (Fukumori et al., 1998).
!
        CALL u2d_bc (ng, tile,                                          &
     &               IminS, ImaxS, JminS, JmaxS,                        &
     &               tl_gradPx)
        CALL v2d_bc (ng, tile,                                          &
     &               IminS, ImaxS, JminS, JmaxS,                        &
     &               tl_gradPy)
!
!  Compute the RHS term for the biconjugate gradient solver.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tl_rhs_r2d(i,j)=-pm(i,j)*pn(i,j)*                           &
     &                      (pmon_u(i+1,j)*tl_gradPx(i+1,j)-            &
     &                       pmon_u(i  ,j)*tl_gradPx(i  ,j)+            &
     &                       pnom_v(i,j+1)*tl_gradPy(i,j+1)-            &
     &                       pnom_v(i,j  )*tl_gradPy(i,j  ))
#  ifdef MASKING
            tl_rhs_r2d(i,j)=tl_rhs_r2d(i,j)*rmask(i,j)
#  endif
          END DO
        END DO

        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               tl_rhs_r2d)
#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_rhs_r2d)
#  endif
!
!  Choose the starting value of zeta_ref.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tl_zeta_ref(i,j)=0.0_r8
          END DO
        END DO
!
!  Call the tangent linear biconjugate gradient solver.
!
        CALL tl_biconj_tile (ng, tile, iTLM,                            &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Lbck,                                      &
     &                       h, pmon_u, pnom_v, pm, pn,                 &
#  ifdef MASKING
     &                       umask, vmask, rmask,                       &
#  endif
     &                       bc_ak, bc_bk, zdf1, zdf2, zdf3,            &
     &                       pc_r2d, r_r2d, br_r2d, p_r2d, bp_r2d,      &
     &                       tl_zeta_ref, tl_rhs_r2d)
!
!  Add balanced free-surface contribution to unbalanced free-surface
!  increment.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)+tl_zeta_ref(i,j)
          END DO
        END DO

# else
!
!  Integrate hydrostatic equation from bottom to surface.
!
        IF (LNM_flag.eq.0) THEN
          cff1=1.0_r8/rho0
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                tl_cff=-cff1*tl_rho(i,j,k)*Hz(i,j,k)
#  ifdef MASKING
                tl_cff=tl_cff*rmask(i,j)
#  endif
                tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)+tl_cff
              END DO
            END DO
          END DO
!
!  Integrate from level of no motion (LNM_depth) to surface or
!  integrate from local bottom if shallower than LNM_depth.
!  Notice that the level of motion depth is tested against the
!  bottom of the grid cell (at W-points) and bracketed between
!  W-points during the interpolation. Also positive depths are
!  used for clarity.
!
        ELSE IF (LNM_flag.eq.1) THEN
          cff1=1.0_r8/rho0
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DO k=1,N(ng)
                zwtop=ABS(z_w(i,j,k  ))
                zwbot=ABS(z_w(i,j,k-1))
                IF (zwbot.lt.LNM_depth(ng)) THEN
                  tl_cff=-cff1*tl_rho(i,j,k)*Hz(i,j,k)
#  ifdef MASKING
                  tl_cff=tl_cff*rmask(i,j)
#  endif
                  tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)+tl_cff
                ELSE IF ((Zwbot.gt.LNM_depth(ng)).and.                  &
     &                   (LNM_depth(ng).ge.Zwtop)) THEN    ! interpolate
                  zphi=ABS(z_r(i,j,k))
                  IF (zphi.ge.LNM_depth(ng)) THEN    ! above cell center
                    IF (k.eq.N(ng)) THEN
                      tl_cff=-cff1*tl_rho(i,j,k)*Hz(i,j,k)
                    ELSE
                      zphi1=ABS(z_r(i,j,k+1))
                      fac=(LNM_depth(ng)-zphi1)/(zphi-zphi1)
                      tl_cff=-cff1*                                     &
     &                       (tl_rho(i,j,k+1)+                          &
     &                        fac*(tl_rho(i,j,k)-tl_rho(i,j,k+1)))*     &
     &                       (LNM_depth(ng)-zwtop)
                    END IF
                  ELSE                               ! below cell center
                    IF (k.eq.1) THEN
                      tl_cff=-cff1*tl_rho(i,j,k)*Hz(i,j,k)
                    ELSE
                      zphi1=ABS(z_r(i,j,k-1))
                      fac=(LNM_depth(ng)-zphi)/(zphi1-zphi)
                      tl_cff=-cff1*                                     &
     &                       (tl_rho(i,j,k)+                            &
     &                        fac*(tl_rho(i,j,k-1)-tl_rho(i,j,k)))*     &
     &                       (zwbot-LNM_depth(ng))
                    END IF
                  END IF
#  ifdef MASKING
                  tl_cff=tl_cff*rmask(i,j)
#  endif
                  tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)+tl_cff
                END IF
              END DO
            END DO
          END DO
        END IF
# endif

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_zeta(:,:,Linp))
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_zeta(:,:,Linp))
# endif
      END IF

      RETURN
      END SUBROUTINE tl_balance_tile

#endif
      END MODULE tl_balance_mod
