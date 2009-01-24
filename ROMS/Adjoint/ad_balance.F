#include "cppdefs.h"
      MODULE ad_balance_mod

#if defined TANGENT && defined BALANCE_OPERATOR
!
!svn $Id: ad_balance.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
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
!      tion, Q. J. R. Meteorol. Soc., submitted.                       !
!      (See also, ECMWR Technical Memorandum # 491, April 2006)        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PRIVATE
      PUBLIC :: ad_balance

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_balance (ng, tile, Lbck, Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
      USE mod_mixing
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
      COUPLING(ng) % Zt_avg1 = OCEAN(ng) % zeta(:,:,Lbck)

      CALL set_depth (ng, tile)
      nrhs(ng)=Lbck
      CALL rho_eos (ng, tile)     
!
# endif
      CALL ad_balance_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Lbck, Linp,                                 &
     &                      GRID(ng) % f,                               &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
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
# ifdef SOLVE3D
     &                      OCEAN(ng) % ad_rho,                         &
     &                      OCEAN(ng) % ad_t,                           &
     &                      OCEAN(ng) % ad_u,                           &
     &                      OCEAN(ng) % ad_v,                           &
# endif
     &                      OCEAN(ng) % ad_zeta)

      RETURN
      END SUBROUTINE ad_balance
!
!***********************************************************************
      SUBROUTINE ad_balance_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Lbck, Linp,                           &
     &                            f, pm, pn,                            &
# ifdef SOLVE3D
     &                            Hz, z_r, z_w,                         &
# endif
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# ifdef SOLVE3D
     &                            alpha, beta, t,                       &
# endif
# ifdef SOLVE3D
     &                            ad_rho, ad_t, ad_u, ad_v,             &
# endif
     &                            ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
      USE ad_exchange_2d_mod
      USE ad_exchange_3d_mod
#  endif
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d, ad_mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck, Linp
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
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
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(out) :: ad_rho(LBi:,LBj:,:)
#  endif

# else

      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
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
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,2,N(ng))
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,2,N(ng))
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(out) :: ad_rho(LBi:UBi,LBj:UBj,N(ng))
#  endif

# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j, k, order

      integer :: Norder = 2                 ! Shapiro filter order

      real(r8) :: fac, fac1, fac2, fac3, gamma
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: ad_cff1, ad_cff2, adfac, adfac1, adfac2
      real(r8) :: dzdT, ml_depth

      real(r8) :: dTdz_min = 0.001_r8
      real(r8) :: ml_depth = 100.0_r8       ! hack value for now

      real(r8), dimension(20) ::  filter_coef =                         &
     &   (/ 2.500000E-1_r8,    6.250000E-2_r8,     1.562500E-2_r8,      &
     &      3.906250E-3_r8,    9.765625E-4_r8,     2.44140625E-4_r8,    &
     &      6.103515625E-5_r8, 1.5258789063E-5_r8, 3.814697E-6_r8,      &
     &      9.536743E-7_r8,    2.384186E-7_r8,     5.960464E-8_r8,      &
     &      1.490116E-8_r8,    3.725290E-9_r8,     9.313226E-10_r8,     &
     &      2.328306E-10_r8,   5.820766E-11_r8,    1.455192E-11_r8,     &
     &      3.637979E-12_r8,   9.094947E-13_r8 /)

      real(r8), dimension(N(ng)) :: dSdT, dSdT_filter

      real(r8), dimension(IminS:ImaxS) :: ad_phie
      real(r8), dimension(IminS:ImaxS) :: ad_phix

# ifdef SALINITY
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dTdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: dSdz
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: ad_gradP
      
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      DO i=IminS,ImaxS
        ad_phie(i)=0.0_r8
        ad_phix(i)=0.0_r8
      END DO
      DO k=1,N(ng)
        DO j=LBj,UBj
          DO i=LBi,UBi
            ad_rho(i,j,k)=0.0_r8
          END DO
        END DO
        DO j=JminS,JmaxS
          DO i=IminS,ImaxS
            ad_gradP(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Add balanced free-surface contribution to unbalanced free-surface
!  increment.  For simpliciyt, integrate hydrostatic equation from
!  surface to bottom. A more rigorous approach involves the solution
!  of an elliptic equation (Fukumori et al., 1998).
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_zeta(:,:,Linp))
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_zeta(:,:,Linp))
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_u2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_zeta(:,:,Linp))
!>
      CALL ad_exchange_u2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_zeta(:,:,Linp))
# endif

      cff1=1.0_r8/rho0
      DO j=Jstr,Jend
# ifdef MASKING
        DO i=Istr,Iend
!>        tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)*rmask(i,j)
!>
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)*rmask(i,j)
        END DO
# endif
        DO k=1,N(ng)-1
          DO i=Istr,Iend
!>          tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)+tl_cff2
!>
            ad_cff2=ad_zeta(i,j,Linp)
!>          tl_cff2=-cff1*tl_rho(i,j,k)*Hz(i,j,k)
!>
            ad_rho(i,j,k)=ad_rho(i,j,k)-cff1*Hz(i,j,k)*ad_cff2
            ad_cff2=0.0_r8
          END DO
        END DO
        DO i=Istr,Iend
!>        tl_zeta(i,j,Linp)=tl_zeta(i,j,Linp)+tl_cff2
!>
          ad_cff2=ad_zeta(i,j,Linp)
!>        tl_cff2=-cff1*tl_rho(i,j,N(ng))*Hz(i,j,N(ng))
!>
          ad_rho(i,j,N(ng))=ad_rho(i,j,N(ng))-cff1*Hz(i,j,N(ng))*ad_cff2
          ad_cff2=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Add balanced velocity contributions to unbalances velocity
!  increments. Use linear pressure gradient formulation based
!  to that found in routine "prsgrd31.h".
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_u(:,:,:,Linp), tl_v(:,:,:,Linp))
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_u(:,:,:,Linp), ad_v(:,:,:,Linp))
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_v3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_v(:,:,:,Linp))
!>
      CALL ad_exchange_v3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_v(:,:,:,Linp))
!>    CALL exchange_u3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_u(:,:,:,Linp))
!>
      CALL ad_exchange_u3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_u(:,:,:,Linp))
# endif
!
!  Compute balanced, interior V-momemtum from baroclinic pressure
!  gradient (differentiate and then vertically integrate).
!
      fac1=0.5_r8*g/rho0
      fac2=g
      fac3=0.25_r8*g/rho0
!
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
# ifdef MASKING
!>          tl_v(i,j,k,Linp)=tl_v(i,j,k,Linp)*vmask(i,j)
!>
            ad_v(i,j,k,Linp)=ad_v(i,j,k,Linp)*vmask(i,j)
# endif
!>          tl_v(i,j,k,Linp)=tl_v(i,j,k,Linp)+                          &
!>   &                       0.25_r8*(tl_gradP(i  ,j-1,k)+              &
!>   &                                tl_gradP(i+1,j-1,k)+              &
!>   &                                tl_gradP(i  ,j  ,k)+              &
!>   &                                tl_gradP(i+1,j  ,k))
!>
            adfac=0.25_r8*ad_v(i,j,k,Linp)
            ad_gradP(i  ,j-1,k)=ad_gradP(i  ,j-1,k)+adfac
            ad_gradP(i+1,j-1,k)=ad_gradP(i+1,j-1,k)+adfac
            ad_gradP(i  ,j  ,k)=ad_gradP(i  ,j  ,k)+adfac
            ad_gradP(i+1,j  ,k)=ad_gradP(i+1,j  ,k)+adfac
          END DO
        END DO
      END DO

      DO j=Jstr-1,Jend
!>      DO k=1,N(ng)-1
!>
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend+1
            cff1=1.0_r8/((z_r(i  ,j,k+1)-z_r(i  ,j,k))*                 &
     &                   (z_r(i-1,j,k+1)-z_r(i-1,j,k)))
            cff2=z_r(i  ,j,k  )-z_r(i-1,j,k  )+                         &
     &           z_r(i  ,j,k+1)-z_r(i-1,j,k+1)
            cff3=z_r(i  ,j,k+1)-z_r(i  ,j,k  )-                         &
     &           z_r(i-1,j,k+1)+z_r(i-1,j,k  )
            gamma=0.125_r8*cff1*cff2*cff3
!
            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i-1,j,k+1))+          &
     &           (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i-1,j,k  ))
!>          tl_gradP(i,j,k)=tl_phix(i)*(pm(i-1,j)+pm(i,j))/             &
!>   &                      (f(i-1,j)+f(i,j))
!>
            ad_phix(i)=ad_phix(i)+                                      &
     &                 ad_gradP(i,j,k)*(pm(i-1,j)+pm(i,j))/             &
     &                 (f(i-1,j)+f(i,j))
            ad_gradP(i,j,k)=0.0_r8
!>          tl_phix(i)=tl_phix(i)+                                      &
!>   &                 fac3*(tl_cff1*cff3-tl_cff2*cff4)
!>
            ad_cff1=fac3*cff3*ad_phix(i)
            ad_cff2=-fac3*cff4*ad_phix(i)
!>          tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-                  &
!>   &              tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
!>          tl_cff1=(1.0_r8+gamma)*(tl_rho(i,j,k+1)-tl_rho(i-1,j,k+1))+ &
!>   &              (1.0_r8-gamma)*(tl_rho(i,j,k  )-tl_rho(i-1,j,k  ))
!>
            adfac1=(1.0_r8+gamma)*ad_cff1
            adfac2=(1.0_r8-gamma)*ad_cff1
            ad_rho(i-1,j,k  )=ad_rho(i-1,j,k  )-adfac2-ad_cff2
            ad_rho(i  ,j,k  )=ad_rho(i  ,j,k  )+adfac2-ad_cff2
            ad_rho(i-1,j,k+1)=ad_rho(i-1,j,k+1)-adfac1+ad_cff2
            ad_rho(i  ,j,k+1)=ad_rho(i  ,j,k+1)+adfac1+ad_cff2
            ad_cff1=0.0_r8
            ad_cff2=0.0_r8
          END DO
        END DO
!
!  Compute balanced, surface V-momentum from baroclinic and barotropic
!  surface pressure gradient.
!
        DO i=Istr,Iend+1
          cff1=z_w(i  ,j,N(ng))-z_r(i  ,j,N(ng))+                        &
     &         z_w(i-1,j,N(ng))-z_r(i-1,j,N(ng))
!>        tl_gradP(i,j,N(ng))=tl_phix(i)*(pm(i-1,j)+pm(i,j))/            &
!>   &                        (f(i-1,j)+f(i,j))
!>
          ad_phix(i)=ad_phix(i)+                                         &
     &               ad_gradP(i,j,N(ng))*(pm(i-1,j)+pm(i,j))/            &
     &               (f(i-1,j)+f(i,j))
          ad_gradP(i,j,N(ng))=0.0_r8
!>        tl_phix(i)=fac1*(tl_rho(i,j,N(ng))-tl_rho(i-1,j,N(ng)))*cff1+  &
!>   &               fac2*(tl_zeta(i,j,Linp)-tl_zeta(i-1,j,Linp))
!>
          adfac1=fac1*cff1*ad_phix(i)
          adfac2=fac2*ad_phix(i)
          ad_rho(i,j,N(ng))=ad_rho(i,j,N(ng))+adfac1
          ad_rho(i-1,j,N(ng))=ad_rho(i-1,j,N(ng))-adfac1
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)+adfac2
          ad_zeta(i-1,j,Linp)=ad_zeta(i-1,j,Linp)-adfac2
          ad_phix(i)=0.0_r8
        END DO
      END DO
!
!  Compute balanced, interior U-momemtum from baroclinic pressure
!  gradient (differentiate and then vertically integrate).
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
# ifdef MASKING
!>          tl_u(i,j,k,Linp)=tl_u(i,j,k,Linp)*umask(i,j)
!>
            ad_u(i,j,k,Linp)=ad_u(i,j,k,Linp)*umask(i,j)
# endif
!>          tl_u(i,j,k,Linp)=tl_u(i,j,k,Linp)-                          &
!>   &                       0.25_r8*(tl_gradP(i-1,j  ,k)+              &
!>   &                                tl_gradP(i  ,j  ,k)+              &
!>   &                                tl_gradP(i-1,j+1,k)+              &
!>   &                                tl_gradP(i  ,j+1,k))
!>
            adfac=0.25_r8*ad_u(i,j,k,Linp)
            ad_gradP(i-1,j  ,k)=ad_gradP(i-1,j  ,k)-adfac
            ad_gradP(i  ,j  ,k)=ad_gradP(i  ,j  ,k)-adfac
            ad_gradP(i-1,j+1,k)=ad_gradP(i-1,j+1,k)-adfac
            ad_gradP(i  ,j+1,k)=ad_gradP(i  ,j+1,k)-adfac
          END DO
        END DO
      END DO

      DO j=Jstr,Jend+1
!>      DO k=1,N(ng)-1
!>
        DO k=N(ng)-1,1,-1
          DO i=Istr-1,Iend
            cff1=1.0_r8/((z_r(i,j  ,k+1)-z_r(i,j  ,k))*                 &
     &                   (z_r(i,j-1,k+1)-z_r(i,j-1,k)))
            cff2=z_r(i,j  ,k  )-z_r(i,j-1,k  )+                         &
     &           z_r(i,j  ,k+1)-z_r(i,j-1,k+1)
            cff3=z_r(i,j  ,k+1)-z_r(i,j  ,k  )-                         &
     &           z_r(i,j-1,k+1)+z_r(i,j-1,k  )
            gamma=0.125_r8*cff1*cff2*cff3
!
            cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i,j-1,k  )
            cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i,j-1,k+1))+          &
     &           (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i,j-1,k  ))
!>          tl_gradP(i,j,k)=tl_phie(i)*(pn(i,j-1)+pn(i,j))/             &
!>   &                      (f(i,j-1)+f(i,j))
!>
            ad_phie(i)=ad_phie(i)+                                      &
     &                 ad_gradP(i,j,k)*(pn(i,j-1)+pn(i,j))/             &
     &                 (f(i,j-1)+f(i,j))
            ad_gradP(i,j,k)=0.0_r8            
!>          tl_phie(i)=tl_phie(i)+                                      &
!>   &                 fac3*(tl_cff1*cff3-tl_cff2*cff4)
!>
            ad_cff1=fac3*cff3*ad_phie(i)
            ad_cff2=-fac3*cff4*ad_phie(i)
!>          tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-                  &
!>   &              tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
!>          tl_cff1=(1.0_r8+gamma)*(tl_rho(i,j,k+1)-tl_rho(i,j-1,k+1))+ &
!>   &              (1.0_r8-gamma)*(tl_rho(i,j,k  )-tl_rho(i,j-1,k  ))
!>
              adfac1=(1.0_r8+gamma)*ad_cff1
              adfac2=(1.0_r8-gamma)*ad_cff1
              ad_rho(i,j-1,k  )=ad_rho(i,j-1,k  )-adfac2-ad_cff2
              ad_rho(i,j  ,k  )=ad_rho(i,j  ,k  )+adfac2-ad_cff2
              ad_rho(i,j-1,k+1)=ad_rho(i,j-1,k+1)-adfac1+ad_cff2
              ad_rho(i,j  ,k+1)=ad_rho(i,j  ,k+1)+adfac1+ad_cff2
              ad_cff1=0.0_r8
              ad_cff2=0.0_r8
          END DO
        END DO
!
!  Compute balanced, surface U-momentum from baroclinic and barotropic
!  surface pressure gradient.
!
        DO i=Istr-1,Iend
          cff1=z_w(i,j  ,N(ng))-z_r(i,j  ,N(ng))+                       &
     &         z_w(i,j-1,N(ng))-z_r(i,j-1,N(ng))
!>        tl_gradP(i,j,N(ng))=tl_phie(i)*(pn(i,j-1)+pn(i,j))/           &
!>   &                        (f(i,j-1)+f(i,j))
!>
          ad_phie(i)=ad_phie(i)+                                        &
     &               ad_gradP(i,j,N(ng))*(pn(i,j-1)+pn(i,j))/           &
     &               (f(i,j-1)+f(i,j))
          ad_gradP(i,j,N(ng))=0.0_r8
!>        tl_phie(i)=fac1*(tl_rho(i,j,N(ng))-tl_rho(i,j-1,N(ng)))*cff1+ &
!>   &               fac2*(tl_zeta(i,j,Linp)-tl_zeta(i,j-1,Linp))
!>
          adfac1=fac1*cff1*ad_phie(i)
          adfac2=fac2*ad_phie(i)
          ad_rho(i,j,N(ng))=ad_rho(i,j,N(ng))+adfac1
          ad_rho(i,j-1,N(ng))=ad_rho(i,j-1,N(ng))-adfac1
          ad_zeta(i,j,Linp)=ad_zeta(i,j,Linp)+adfac2
          ad_zeta(i,j-1,Linp)=ad_zeta(i,j-1,Linp)-adfac2
          ad_phie(i)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute balanced density anomaly increment using linearized equation
!  of state.  The thermal expansion and saline contraction coefficients
!  are computed from the background state.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_rho)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_rho)
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_r3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_rho)
!>
      CALL ad_exchange_r3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_rho)
# endif

      DO j=JstrR,JendR
        DO k=1,N(ng)
          DO i=IstrR,IendR
# ifdef MASKING
!>          tl_rho(i,j,k)=tl_rho(i,j,k)*rmask(i,j)
!>
            ad_rho(i,j,k)=ad_rho(i,j,k)*rmask(i,j)
# endif
# ifdef SALINITY
!>          tl_rho(i,j,k)=tl_rho(i,j,k)+                                &
!>   &                    rho0*beta(i,j)*tl_t(i,j,k,Linp,isalt)
!>
            ad_t(i,j,k,Linp,isalt)=ad_t(i,j,k,Linp,isalt)+              & 
     &                             rho0*beta(i,j)*ad_rho(i,j,k)
# endif
!>          tl_rho(i,j,k)=-rho0*alpha(i,j)*tl_t(i,j,k,Linp,itemp)
!>
            ad_t(i,j,k,Linp,itemp)=ad_t(i,j,k,Linp,itemp)-              &
     &                             rho0*alpha(i,j)*ad_rho(i,j,k)
            ad_rho(i,j,k)=0.0
          END DO
        END DO
      END DO

# ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Compute balance salinity contribution.
!-----------------------------------------------------------------------
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
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                             &
     &                  Hz(i,j,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*Hz(i,j,k+1)
            dTdz(i,j,k)=cff*(6.0_r8*(t(i,j,k+1,Lbck,itemp)-             &
     &                               t(i,j,k  ,Lbck,itemp))-            &
     &                       Hz(i,j,k)*dTdz(i,j,k-1))
            dSdz(i,j,k)=cff*(6.0_r8*(t(i,j,k+1,Lbck,isalt)-             &
     &                               t(i,j,k  ,Lbck,isalt))-            &
     &                       Hz(i,j,k)*dSdz(i,j,k-1))
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
!!
!! Compute depth of the mixed layer from temperature shear.
!!
!!      DO i=IstrR,IendR
!!        dTdz_min=dTdz(i,j,N(ng))
!!        DO k=N(ng)-1,1,-1
!!          IF (ABS(dTdz(i,j,)).lt.dTdz_min) THEN
!!            dTdz_min=dTdz(i,j,k)
!!            ml_depth=ABS(z_w(i,j,k))
!!            EXIT
!!          END IF
!!        END DO
!!      END DO
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
#  ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_t(:,:,:,Linp,isalt))
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_t(:,:,:,Linp,isalt))
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_r3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_t(:,:,:,Linp,isalt)
!>
      CALL ad_exchange_r3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_t(:,:,:,Linp,isalt)
#  endif
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          DO k=1,N(ng)
            IF (ABS(dTdz(i,j,k)).lt.dTdz_min) THEN
              dzdT=0.0_r8
            ELSE
              dzdT=1.0_r8/(0.5_r8*(dTdz(i,j,k-1)+                       &
     &                             dTdz(i,j,k  )))
            END IF
            dSdT(k)=(0.5_r8*(dSdz(i,j,k-1)+                             &
     &                       dSdz(i,j,k  )))*dzdT
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
            cff=(1.0_r8-EXP(z_r(i,j,k)/ml_depth))*dSdT(k)
#  ifdef MASKING
!>          tl_t(i,j,k,Linp,isalt)=tl_t(i,j,k,Linp,isalt)*rmask(i,j)
!>
            ad_t(i,j,k,Linp,isalt)=ad_t(i,j,k,Linp,isalt)*rmask(i,j)
#  endif
!>          tl_t(i,j,k,Linp,isalt)=tl_t(i,j,k,Linp,isalt)+              &
!>   &                             cff*tl_t(i,j,k,Linp,itemp)
!>
            ad_t(i,j,k,Linp,itemp)=ad_t(i,j,k,Linp,itemp)+              &
     &                             cff*ad_t(i,j,k,Linp,isalt)
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE ad_balance_tile

#endif
      END MODULE ad_balance_mod

