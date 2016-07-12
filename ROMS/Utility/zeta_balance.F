#include "cppdefs.h"
      MODULE zeta_balance_mod

#if defined BALANCE_OPERATOR  && defined ZETA_ELLIPTIC
!svn $Id: zeta_balance.F 795 2016-05-11 01:42:43Z arango $
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains routines to solve the elliptic equation for    !
!  free surface (SSH) that is used as part of the  balance operator    !
!  during 4D-Var.                                                      !
!                                                                      !
!  The balanced (baroclinic) SSH 4DVar increment is approximated as    !
!  the sum of two terms:  a baroclinic term that depends on density    !
!  and a  barotropic  term  that  depends  on  the depth-integrated    !
!  transport.                                                          !
!                                                                      !
!   div (h grad(zeta)) = - div (int{int{grad(rho)z'/rho0) dz'} dz})    !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fukumori, I., R. Raghunath and L. Fu, 1998:   Nature of global    !
!      large-scale sea level variability in relation to atmospheric    !
!      forcing: a modeling study, J. Geophys. Res., 103, 5493-5512.    !
!                                                                      !
!    Weaver, A.T., C. Deltel, E. Machu, S. Ricci, and N. Daget, 2005:  !
!      A multivariate balance operator for variational data assimila-  !
!      tion, Q. J. R. Meteorol. Soc., 131, 3605-3625.                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE

      PUBLIC :: balance_ref
      PUBLIC :: biconj
      PUBLIC :: ad_biconj_tile
      PUBLIC :: tl_biconj_tile
      PUBLIC :: r2d_bc
      PUBLIC :: ad_r2d_bc
      PUBLIC :: u2d_bc
      PUBLIC :: v2d_bc

      CONTAINS
!
!***********************************************************************
      SUBROUTINE balance_ref (ng, tile, Lbck)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
      USE mod_mixing
# endif
      USE mod_fourdvar
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
      integer, intent(in) :: ng, tile, Lbck
!
!  Local variable declarations.
!
# include "tile.h"

# ifdef SOLVE3D
!
!  Compute background state thickness, depth arrays, and density fields.
!  Use zero value free-surface.
!
      COUPLING(ng) % Zt_avg1 = 0.0_r8

      CALL set_depth (ng, tile)
      nrhs(ng)=Lbck
      CALL rho_eos (ng, tile)
!
# endif
      CALL balance_ref_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Lbck,                                      &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % pmon_r,                         &
     &                       GRID(ng) % pnom_r,                         &
     &                       GRID(ng) % pmon_u,                         &
     &                       GRID(ng) % pnom_v,                         &
     &                       GRID(ng) % h,                              &
# ifdef SOLVE3D
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % z_w,                            &
# endif
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef SOLVE3D
     &                       MIXING(ng) % alpha,                        &
     &                       MIXING(ng) % beta,                         &
     &                       OCEAN(ng) % rho,                           &
# endif
     &                       OCEAN(ng) %  zeta,                         &
     &                       FOURDVAR(ng) % rhs_r2d)

      RETURN
      END SUBROUTINE balance_ref
!
!***********************************************************************
      SUBROUTINE balance_ref_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Lbck,                                &
     &                             pm, pn, pmon_r, pnom_r,              &
     &                             pmon_u, pnom_v, h,                   &
# ifdef SOLVE3D
     &                             Hz, z_r, z_w,                        &
# endif
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef SOLVE3D
     &                             alpha, beta, rho,                    &
# endif
     &                             zeta,                                &
     &                             rhs_r2d)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_grid
!
      USE exchange_2d_mod
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
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
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: rhs_r2d(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
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
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: rhs_r2d(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: fac, fac1, fac2, fac3, gamma
      real(r8) :: cff, cff1, cff2, cff3, cff4

      real(r8), dimension(IminS:ImaxS) :: phie
      real(r8), dimension(IminS:ImaxS) :: phix

      real(r8), dimension(IminS:ImaxS,1:N(ng)) :: phi

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradPx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradPy

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute the RHS side term of the elliptic equation for the
!  biconjugate gradient solver.
!-----------------------------------------------------------------------
!
!  NOTE: fac2=0 because the balanced component should consist of the
!  baroclinic pressure gradient only.
!
      fac1=0.5_r8*g/rho0
!!    fac2=1000.0_r8*g/rho0
      fac2=0.0_r8
      fac3=0.25_r8*g/rho0

      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          gradPx(i,j)=0.0_r8
          gradPy(i,j)=0.0_r8
        END DO
      END DO
!
!  Compute balanced, surface U-momentum from baroclinic and barotropic
!  surface pressure gradient.
!
      DO j=Jstr,Jend+1
        DO i=Istr-1,Iend
          cff1=z_w(i,j  ,N(ng))-z_r(i,j  ,N(ng))+                       &
     &         z_w(i,j-1,N(ng))-z_r(i,j-1,N(ng))
          phie(i)=fac1*(rho(i,j,N(ng))-rho(i,j-1,N(ng)))*cff1+          &
     &            fac2*(zeta(i,j,Lbck)-zeta(i,j-1,Lbck))
          phi(i,N(ng))=phie(i)
        END DO
!
!  Compute balance, interior U-momentum from baroclinic pressure
!  gradient (differentiate and then vertically integrate).
!
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
            cff1=(1.0_r8+gamma)*(rho(i,j,k+1)-rho(i,j-1,k+1))+          &
     &              (1.0_r8-gamma)*(rho(i,j,k  )-rho(i,j-1,k  ))
            cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-                           &
     &           rho(i,j,k  )-rho(i,j-1,k  )
            cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i,j-1,k  )
            cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i,j-1,k+1))+          &
     &           (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i,j-1,k  ))
            phie(i)=phie(i)+fac3*(cff1*cff3-cff2*cff4)
            phi(i,k)=phie(i)
          END DO
        END DO
!
!  Integrate from bottom to surface.
!
        DO k=1,N(ng)
          DO i=Istr-1,Iend
            cff=0.5_r8*(Hz(i,j-1,k)+Hz(i,j,k))*phi(i,k)
# ifdef MASKING
            cff=cff*vmask(i,j)
# endif
            gradPy(i,j)=gradPy(i,j)+cff
          END DO
        END DO
      END DO
!
!  Compute balanced, surface V-momentum from baroclinic and barotropic
!  surface pressure gradient.
!
      DO j=Jstr-1,Jend
        DO i=Istr,Iend+1
          cff1=z_w(i  ,j,N(ng))-z_r(i  ,j,N(ng))+                       &
     &         z_w(i-1,j,N(ng))-z_r(i-1,j,N(ng))
          phix(i)=fac1*(rho(i,j,N(ng))-rho(i-1,j,N(ng)))*cff1+          &
     &            fac2*(zeta(i,j,Lbck)-zeta(i-1,j,Lbck))
          phi(i,N(ng))=phix(i)
        END DO
!
!  Compute balance, interior V-momentum from baroclinic pressure
!  gradient (differentiate and then vertically integrate).
!
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
            cff1=(1.0_r8+gamma)*(rho(i,j,k+1)-rho(i-1,j,k+1))+          &
     &           (1.0_r8-gamma)*(rho(i,j,k  )-rho(i-1,j,k  ))
            cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i-1,j,k+1))+          &
     &           (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i-1,j,k  ))
            phix(i)=phix(i)+fac3*(cff1*cff3-cff2*cff4)
            phi(i,k)=phix(i)
          END DO
        END DO
!
!  Integrate from bottom to surface.
!
        DO k=1,N(ng)
          DO i=Istr,Iend+1
            cff=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k))*phi(i,k)
#  ifdef MASKING
            cff=cff*umask(i,j)
#  endif
            gradPx(i,j)=gradPx(i,j)+cff
          END DO
        END DO
      END DO
!
!  Impose zero gradient conditions at open boundaries.
!
      CALL u2d_bc (ng, tile,                                            &
     &             IminS, ImaxS, JminS, JmaxS,                          &
     &             gradPx)
      CALL v2d_bc (ng, tile,                                            &
     &             IminS, ImaxS, JminS, JmaxS,                          &
     &             gradPy)
!
!  Compute the RHS term, -div(phi_bar), for the biconjugate gradient
!  solver.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          rhs_r2d(i,j)=-pm(i,j)*pn(i,j)*                                &
     &                 (pmon_u(i+1,j)*gradPx(i+1,j)-                    &
     &                  pmon_u(i  ,j)*gradPx(i  ,j)+                    &
     &                  pnom_v(i,j+1)*gradPy(i,j+1)-                    &
     &                  pnom_v(i,j  )*gradPy(i,j  ))
# ifdef MASKING
          rhs_r2d(i,j)=rhs_r2d(i,j)*rmask(i,j)
# endif
        END DO
      END DO

      CALL r2d_bc (ng, tile,                                            &
     &             LBi, UBi, LBj, UBj,                                  &
     &             rhs_r2d)
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rhs_r2d)
# endif

      RETURN
      END SUBROUTINE balance_ref_tile

!
!***********************************************************************
      SUBROUTINE biconj (ng, tile, model, Lbck)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
      USE mod_mixing
# endif
      USE mod_fourdvar
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
      integer, intent(in) :: ng, tile, Lbck, model
!
!  Local variable declarations.
!
# include "tile.h"
!
!   Call the biconjugate gradient solver to compute "zeta_ref".
!
      CALL biconj_tile (ng, tile, model,                                &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  Lbck,                                           &
     &                  GRID(ng)% h,                                    &
     &                  GRID(ng) % pmon_u,                              &
     &                  GRID(ng) % pnom_v,                              &
     &                  GRID(ng) % pm,                                  &
     &                  GRID(ng) % pn,                                  &
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % vmask,                               &
     &                  GRID(ng) % rmask,                               &
# endif
     &                  FOURDVAR(ng) % bc_ak,                           &
     &                  FOURDVAR(ng) % bc_bk,                           &
     &                  FOURDVAR(ng) % zdf1,                            &
     &                  FOURDVAR(ng) % zdf2,                            &
     &                  FOURDVAR(ng) % zdf3,                            &
     &                  FOURDVAR(ng) % pc_r2d,                          &
     &                  FOURDVAR(ng) % r_r2d,                           &
     &                  FOURDVAR(ng) % br_r2d,                          &
     &                  FOURDVAR(ng) % p_r2d,                           &
     &                  FOURDVAR(ng) % bp_r2d,                          &
     &                  OCEAN(ng) % zeta,                               &
     &                  FOURDVAR(ng) % zeta_ref,                        &
     &                  FOURDVAR(ng) % rhs_r2d)

      RETURN
      END SUBROUTINE biconj
!
!***********************************************************************
      SUBROUTINE biconj_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        Lbck,                                     &
     &                        h, pmon_u, pnom_v, pm, pn,                &
# ifdef MASKING
     &                        umask, vmask, rmask,                      &
# endif
     &                        bc_ak, bc_bk, zdf1, zdf2, zdf3,           &
     &                        pc_r2d, r_r2d, br_r2d, p_r2d, bp_r2d,     &
     &                        zeta, r2d_ref, rhs_r2d)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_grid
      USE mod_iounits
      USE mod_parallel
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: bc_ak(:)
      real(r8), intent(inout) :: bc_bk(:)
      real(r8), intent(inout) :: zdf1(:)
      real(r8), intent(inout) :: zdf2(:)
      real(r8), intent(inout) :: zdf3(:)
      real(r8), intent(inout) :: pc_r2d(LBi:,LBj:)
      real(r8), intent(inout) :: r_r2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: br_r2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: p_r2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: bp_r2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rhs_r2d(LBi:,LBj:)

      real(r8), intent(inout) :: r2d_ref(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: bc_ak(Nbico(ng))
      real(r8), intent(inout) :: bc_bk(Nbico(ng))
      real(r8), intent(inout) :: zdf1(Nbico(ng))
      real(r8), intent(inout) :: zdf2(Nbico(ng))
      real(r8), intent(inout) :: zdf3(Nbico(ng))
      real(r8), intent(inout) :: pc_r2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: r_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: br_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: p_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: bp_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: rhs_r2d(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: r2d_ref(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      logical :: Ltrans

      integer :: i, j, it

      real(r8) :: error, zdf4, zdf5

      real(r8), dimension(LBi:UBi,LBj:UBj) :: bv_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: v_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: z1_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: z2_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: z3_r2d

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Biconjugate gradient algorithm.
!-----------------------------------------------------------------------
!
!  Initialize local arrays.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          v_r2d(i,j)=0.0_r8
          bv_r2d(i,j)=0.0_r8
          z1_r2d(i,j)=0.0_r8
          z2_r2d(i,j)=0.0_r8
          z3_r2d(i,j)=0.0_r8
        END DO
      END DO
!
!  Choose the starting value of "zeta_ref". Use background free-surface.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          r2d_ref(i,j)=zeta(i,j,Lbck)
# ifdef MASKING
          r2d_ref(i,j)=r2d_ref(i,j)*rmask(i,j)
# endif
        END DO
      END DO

      CALL r2d_bc (ng, tile,                                            &
     &             LBi, UBi, LBj, UBj,                                  &
     &             r2d_ref)
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    r2d_ref)
# endif
!
!  Compute starting value of divergence operator:
!     z1_r2d = div[h grad(r2d_ref)].
!
      Ltrans=.FALSE.
      CALL r2d_oper (ng, tile,                                          &
     &               Ltrans, Lbck,                                      &
     &               LBi, UBi, LBj, UBj,                                &
     &               IminS, ImaxS, JminS, JmaxS,                        &
# ifdef MASKING
     &               umask, vmask, rmask,                               &
# endif
     &               h,                                                 &
     &               pmon_u, pnom_v, pm, pn,                            &
     &               pc_r2d, r2d_ref, z1_r2d)
!
!  Set the initial values for residual vectors "r" and "br". Then, use
!  recurrence relationship to compute direction vectors "p" and "bp".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          r_r2d(i,j,1)=rhs_r2d(i,j)-z1_r2d(i,j)
          br_r2d(i,j,1)=r_r2d(i,j,1)
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
          p_r2d(i,j,1)=r_r2d(i,j,1)/pc_r2d(i,j)
          bp_r2d(i,j,1)=br_r2d(i,j,1)/pc_r2d(i,j)
        END DO
      END DO
!
!=======================================================================
!  Iterate.
!=======================================================================
!
      ITERATE : DO it=1,Nbico(ng)-1

        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=p_r2d(i,j,it)
          END DO
        END DO

        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               z1_r2d)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      z1_r2d)
# endif
!
!  Compute divergence operator: v_r2d = div[h grad(z1_r2d)].
!
        Ltrans=.FALSE.
        CALL r2d_oper (ng, tile,                                        &
     &                 Ltrans, Lbck,                                    &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
# ifdef MASKING
     &                 umask, vmask, rmask,                             &
# endif
     &                 h,                                               &
     &                 pmon_u, pnom_v, pm, pn,                          &
     &                 pc_r2d, z1_r2d, v_r2d)
!
!  Compute dot products and "bc_ak" coefficient.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=r_r2d(i,j,it)/pc_r2d(i,j)
            z2_r2d(i,j)=br_r2d(i,j,it)
            z3_r2d(i,j)=bp_r2d(i,j,it)
          END DO
        END DO

        CALL r2d_dotp (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 zdf1(it),                                        &
# ifdef MASKING
     &                 rmask,                                           &
# endif
     &                 z2_r2d, z1_r2d)

        CALL r2d_dotp (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 zdf2(it),                                        &
# ifdef MASKING
     &                 rmask,                                           &
# endif
     &                 z3_r2d, v_r2d)

        bc_ak(it)=zdf1(it)/zdf2(it)
!
!  Solve for new iterate of "r2d_ref".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            r2d_ref(i,j)=r2d_ref(i,j)+bc_ak(it)*p_r2d(i,j,it)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Test for convergence: use "bv_r2d" as temporary storage.
!-----------------------------------------------------------------------
!
        IF (it.eq.Nbico(ng)-1) THEN

          CALL r2d_bc (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 r2d_ref)
# ifdef DISTRIBUTE
          CALL mp_exchange2d (ng, tile, model, 1,                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        r2d_ref)
# endif
!
!  Compute divergence operator: bv_r2d = div[h grad(r2d_ref)].
!
          Ltrans=.FALSE.
          CALL r2d_oper (ng, tile,                                      &
     &                   Ltrans, Lbck,                                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
# ifdef MASKING
     &                   umask, vmask, rmask,                           &
# endif
     &                   h,                                             &
     &                   pmon_u, pnom_v, pm, pn,                        &
     &                   pc_r2d, r2d_ref, bv_r2d)
!
!  Compute dot products and report convergence value.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              bv_r2d(i,j)=bv_r2d(i,j)-rhs_r2d(i,j)
            END DO
          END DO

          CALL r2d_dotp (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   zdf4,                                          &
# ifdef MASKING
     &                   rmask,                                         &
# endif
     &                   bv_r2d, bv_r2d)

          CALL r2d_dotp (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   zdf5,                                          &
# ifdef MASKING
     &                   rmask,                                         &
# endif
     &                   rhs_r2d, rhs_r2d)

          IF (Master) THEN
            error=SQRT(zdf4/zdf5)
            WRITE (stdout,10) error
 10         FORMAT (/,'  BICONJ - balance operator, error in ',         &
     &              'reference free-surface = ',1p,e14.7)
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Compute new (it+1) residual and direction vectors.
!-----------------------------------------------------------------------
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=bp_r2d(i,j,it)
          END DO
        END DO

        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               z1_r2d)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      z1_r2d)
# endif
!
!  Compute divergence operator, tranpose: bv_r2d = div[h grad(z1_r2d)].
!  Need to call "ad_r2d_bc" here since Ltrans is TRUE.
# ifdef DISTRIBUTE
!  Need to call "ad_mp_exchange" here since Ltrans is TRUE.
# endif
!
        Ltrans=.TRUE.
        CALL r2d_oper (ng, tile,                                        &
     &                 Ltrans, Lbck,                                    &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
# ifdef MASKING
     &                 umask, vmask, rmask,                             &
# endif
     &                 h,                                               &
     &                 pmon_u, pnom_v, pm, pn,                          &
     &                 pc_r2d, z1_r2d, bv_r2d)

# ifdef DISTRIBUTE
        CALL ad_mp_exchange2d (ng, tile, model, 1,                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         bv_r2d)
# endif
        CALL ad_r2d_bc (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bv_r2d)
!
!  Compute new residual vectors "r" and "br".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            r_r2d(i,j,it+1)=r_r2d(i,j,it)-bc_ak(it)*v_r2d(i,j)
            br_r2d(i,j,it+1)=br_r2d(i,j,it)-bc_ak(it)*bv_r2d(i,j)
          END DO
        END DO
!
!  Compute dot product and "bc_ak" coefficient.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=r_r2d(i,j,it+1)/pc_r2d(i,j)
            z2_r2d(i,j)=br_r2d(i,j,it+1)
          END DO
        END DO

        CALL r2d_dotp (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 zdf3(it),                                        &
# ifdef MASKING
     &                 rmask,                                           &
# endif
     &                 z1_r2d, z2_r2d)

        bc_bk(it)=zdf3(it)/zdf1(it)
!
!  Use recurrence relationship to compute new direction vectors
!  "p" and "bp".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            p_r2d(i,j,it+1)=r_r2d(i,j,it+1)/pc_r2d(i,j)+                &
     &                      bc_bk(it)*p_r2d(i,j,it)
            bp_r2d(i,j,it+1)=br_r2d(i,j,it+1)/pc_r2d(i,j)+              &
     &                       bc_bk(it)*bp_r2d(i,j,it)
          END DO
        END DO

      END DO ITERATE

      RETURN
      END SUBROUTINE biconj_tile

!
!***********************************************************************
      SUBROUTINE r2d_oper (ng, tile,                                    &
     &                     Ltrans, Lbck,                                &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
# ifdef MASKING
     &                     umask, vmask, rmask,                         &
# endif
     &                     h,                                           &
     &                     pmon_u, pnom_v, pm, pn,                      &
     &                     pc_r2d,                                      &
     &                     r2d_in, r2d_out)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lbck
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      logical, intent(in) :: Ltrans
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(inout) :: pc_r2d(LBi:,LBj:)
      real(r8), intent(inout) :: r2d_in(LBi:,LBj:)

      real(r8), intent(out) :: r2d_out(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: r2d_in(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: r2d_out(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff, fac

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute divergence XI- and ETA-components: transposed operator.
!                                             (Ltrans=.TRUE.)
!-----------------------------------------------------------------------
!
      cff=0.5_r8*g

      IF (Ltrans) THEN
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
            FX(i,j)=0.0_r8
            FE(i,j)=0.0_r8
            r2d_out(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          r2d_out(i,j)=pm(i,j)*pn(i,j)*                               &
!>   &                   (FX(i+1,j)-FX(i,j)+                            &
!>   &                    FE(i,j+1)-FE(i,j))
!>
            fac=pm(i,j)*pn(i,j)*r2d_in(i,j)
# ifdef MASKING
            fac=fac*rmask(i,j)
# endif
            FX(i  ,j  )=FX(i  ,j  )-fac
            FX(i+1,j  )=FX(i+1,j  )+fac
            FE(i  ,j  )=FE(i  ,j  )-fac
            FE(i  ,j+1)=FE(i  ,j+1)+fac
          END DO
        END DO

        DO j=Jstr,Jend+1
          DO i=Istr,Iend
!>          FE(i,j)=cff*pnom_v(i,j)*(h(i,j)+h(i,j-1))*                  &
!>   &              (r2d_in(i,j)-r2d_in(i,j-1))
!>
            fac=cff*pnom_v(i,j)*(h(i,j)+h(i,j-1))*FE(i,j)
# ifdef MASKING
            fac=fac*vmask(i,j)
# endif
            r2d_out(i,j-1)=r2d_out(i,j-1)-fac
            r2d_out(i,j  )=r2d_out(i,j  )+fac
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
!>          FX(i,j)=cff*pmon_u(i,j)*(h(i,j)+h(i-1,j))*                  &
!>  &               (r2d_in(i,j)-r2d_in(i-1,j))
!>
            fac=cff*pmon_u(i,j)*(h(i,j)+h(i-1,j))*FX(i,j)
# ifdef MASKING
            fac=fac*umask(i,j)
# endif
            r2d_out(i-1,j)=r2d_out(i-1,j)-fac
            r2d_out(i  ,j)=r2d_out(i  ,j)+fac
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute divergence XI- and ETA-components: regular operator
!                                             (Ltrans=.FALSE.).
!-----------------------------------------------------------------------
!
      ELSE

        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            FX(i,j)=cff*pmon_u(i,j)*(h(i,j)+h(i-1,j))*                  &
     &              (r2d_in(i,j)-r2d_in(i-1,j))
# ifdef MASKING
            FX(i,j)=FX(i,j)*umask(i,j)
# endif
          END DO
        END DO

        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            FE(i,j)=cff*pnom_v(i,j)*(h(i,j)+h(i,j-1))*                  &
     &              (r2d_in(i,j)-r2d_in(i,j-1))
# ifdef MASKING
            FE(i,j)=FE(i,j)*vmask(i,j)
# endif
          END DO
        END DO

        DO j=Jstr,Jend
          DO i=Istr,Iend
             r2d_out(i,j)=pm(i,j)*pn(i,j)*                              &
     &                    (FX(i+1,j)-FX(i,j)+                           &
     &                     FE(i,j+1)-FE(i,j))
# ifdef MASKING
             r2d_out(i,j)=r2d_out(i,j)*rmask(i,j)
# endif
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Compute scale coefficient.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          pc_r2d(i,j)=-pm(i,j)*pn(i,j)*                                 &
     &                cff*(pnom_v(i  ,j+1)*(h(i  ,j+1)+h(i  ,j  ))+     &
     &                     pnom_v(i  ,j  )*(h(i  ,j  )+h(i  ,j-1))+     &
     &                     pmon_u(i+1,j  )*(h(i+1,j  )+h(i  ,j  ))+     &
     &                     pmon_u(i  ,j  )*(h(i  ,j  )+h(i-1,j  )))
        END DO
      END DO

      RETURN
      END SUBROUTINE r2d_oper

!
!***********************************************************************
      SUBROUTINE r2d_bc (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
# ifdef R2D_GRADIENT
            A(Iend+1,j)=A(Iend,j)
# else
            A(Iend+1,j)=0.0_r8
# endif
          END DO
        END IF
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
# ifdef R2D_GRADIENT
            A(Istr-1,j)=A(Istr,j)
# else
            A(Istr-1,j)=0.0_r8
# endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
# ifdef R2D_GRADIENT
            A(i,Jend+1)=A(i,Jend)
# else
            A(i,Jend+1)=0.0_r8
# endif
          END DO
        END IF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
# ifdef R2D_GRADIENT
            A(i,Jstr-1)=A(i,Jstr)
# else
            A(i,Jstr-1)=0.0_r8
# endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Istr-1,Jstr-1)=0.5_r8*(A(Istr  ,Jstr-1)+                    &
     &                             A(Istr-1,Jstr  ))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                    &
     &                             A(Iend+1,Jstr  ))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                    &
     &                             A(Istr  ,Jend+1))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                    &
     &                             A(Iend  ,Jend+1))
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE r2d_bc

!
!***********************************************************************
      SUBROUTINE u2d_bc (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_u2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed or gradient
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            A(Iend+1,j)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            A(Istr,j)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrU,Iend
            A(i,Jend+1)=0.0_r8
          END DO
        END IF

        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            A(i,Jstr-1)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Istr  ,Jstr-1)=0.5_r8*(A(Istr+1,Jstr-1)+                    &
     &                             A(Istr  ,Jstr  ))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(Iend+1,Jstr-1)=0.5_r8*(A(Iend  ,Jstr-1)+                    &
     &                             A(Iend+1,Jstr  ))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Istr  ,Jend+1)=0.5_r8*(A(Istr  ,Jend  )+                    &
     &                             A(Istr+1,Jend+1))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                    &
     &                             A(Iend  ,Jend+1))
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE u2d_bc

!
!***********************************************************************
      SUBROUTINE v2d_bc (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   A)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_v2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: A(LBi:,LBj:)
#else
      real(r8), intent(inout) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  East-West boundary conditions: Closed (free-slip/no-slip) or
!  gradient.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            A(Iend+1,j)=0.0_r8
          END DO
        END IF

        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrV,Jend
            A(Istr-1,j)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South boundary conditions: Closed or Gradient.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            A(i,Jend+1)=0.0_r8
          END DO
        END IF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            A(i,Jstr)=0.0_r8
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          A(Istr-1,Jstr  )=0.5_r8*(A(Istr  ,Jstr  )+                    &
     &                             A(Istr-1,Jstr+1))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          A(Iend+1,Jstr  )=0.5_r8*(A(Iend  ,Jstr  )+                    &
     &                             A(Iend+1,Jstr+1))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          A(Istr-1,Jend+1)=0.5_r8*(A(Istr-1,Jend  )+                    &
     &                             A(Istr  ,Jend+1))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          A(Iend+1,Jend+1)=0.5_r8*(A(Iend+1,Jend  )+                    &
     &                             A(Iend  ,Jend+1))
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          A)
      END IF

      RETURN
      END SUBROUTINE v2d_bc

!
!***********************************************************************
      SUBROUTINE r2d_dotp (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     DotProd,                                     &
# ifdef MASKING
     &                     rmask,                                       &
# endif
     &                     s1_zeta, s2_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: s1_zeta(LBi:,LBj:)
      real(r8), intent(in) :: s2_zeta(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: s1_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s2_zeta(LBi:UBi,LBj:UBj)
# endif
!
      real(r8), intent(out) :: DotProd
!
!  Local variable declarations.
!
      integer :: NSUB, i, j

      real(r8) :: cff
      real(r8) :: my_DotProd
# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute dot product between "s1_zeta" and "s2_zeta".
!-----------------------------------------------------------------------
!
      my_DotProd=0.0_r8

      DO j=Jstr,Jend
        DO i=Istr,Iend
          cff=s1_zeta(i,j)*s2_zeta(i,j)
# ifdef MASKING
          cff=cff*rmask(i,j)
# endif
          my_DotProd=my_DotProd+cff
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Perform parallel global reduction operations.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
# else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
# endif
!$OMP CRITICAL (DOT_PROD)
      IF (tile_count.eq.0) THEN
        DotProd=0.0_r8
      END IF
      DotProd=DotProd+my_DotProd
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, DotProd, op_handle)
# endif
      END IF
!$OMP END CRITICAL (DOT_PROD)

      RETURN
      END SUBROUTINE r2d_dotp

!
!***********************************************************************
      SUBROUTINE tl_biconj_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Lbck,                                  &
     &                           h, pmon_u, pnom_v, pm, pn,             &
# ifdef MASKING
     &                           umask, vmask, rmask,                   &
# endif
     &                           bc_ak, bc_bk, zdf1, zdf2, zdf3,        &
     &                           pc_r2d, r_r2d, br_r2d, p_r2d, bp_r2d,  &
     &                           tl_r2d_ref, tl_rhs_r2d)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_grid
      USE mod_parallel
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
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
      real(r8), intent(inout) :: tl_r2d_ref(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: bc_ak(Nbico(ng))
      real(r8), intent(inout) :: bc_bk(Nbico(ng))
      real(r8), intent(inout) :: zdf1(Nbico(ng))
      real(r8), intent(inout) :: zdf2(Nbico(ng))
      real(r8), intent(inout) :: zdf3(Nbico(ng))
      real(r8), intent(inout) :: pc_r2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: r_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: br_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: p_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: tl_bp_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: tl_rhs_r2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_r2d_ref(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      logical :: Ltrans

      integer :: i, j, it

      real(r8) :: tl_zdf1, tl_zdf2, tl_zdf3, tl_bc_ak, tl_bc_bk

      real(r8), dimension(LBi:UBi,LBj:UBj) :: bv_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: v_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: z1_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: z2_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: z3_r2d

      real(r8), dimension(LBi:UBi,LBj:UBj) :: tl_bp_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: tl_br_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: tl_bv_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: tl_p_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: tl_r_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: tl_v_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: tl_z1_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: tl_z2_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: tl_z3_r2d

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Tangent linear of biconjugate gradient algorithm.
!-----------------------------------------------------------------------
!
!  Initialize local arrays.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          bv_r2d(i,j)=0.0_r8
          v_r2d(i,j)=0.0_r8
          z1_r2d(i,j)=0.0_r8
          z2_r2d(i,j)=0.0_r8
          z3_r2d(i,j)=0.0_r8

          tl_bp_r2d(i,j)=0.0_r8
          tl_br_r2d(i,j)=0.0_r8
          tl_bv_r2d(i,j)=0.0_r8
          tl_p_r2d(i,j)=0.0_r8
          tl_r_r2d(i,j)=0.0_r8
          tl_v_r2d(i,j)=0.0_r8
          tl_z1_r2d(i,j)=0.0_r8
          tl_z2_r2d(i,j)=0.0_r8
          tl_z3_r2d(i,j)=0.0_r8
        END DO
      END DO
!
!  Compute starting value of divergence TLM operator:
!    tl_z1_r2d = div[h grad(tl_r2d_ref)].
!
      CALL r2d_bc (ng, tile,                                            &
     &             LBi, UBi, LBj, UBj,                                  &
     &             tl_r2d_ref)
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_r2d_ref)
# endif

      CALL tl_r2d_oper (ng, tile,                                       &
     &                  Lbck,                                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
# ifdef MASKING
     &                  umask, vmask, rmask,                            &
# endif
     &                  h,                                              &
     &                  pmon_u, pnom_v, pm, pn,                         &
     &                  pc_r2d, tl_r2d_ref, tl_z1_r2d)
!
!  Set the initial values for residual vectors "tl_r" and "tl_br". Then,
!  use recurrence relationship to compute direction vectors "tl_p" and
!  "tl_bp".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          tl_r_r2d(i,j)=tl_rhs_r2d(i,j)-tl_z1_r2d(i,j)
          tl_br_r2d(i,j)=tl_r_r2d(i,j)
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
          tl_p_r2d(i,j)=tl_r_r2d(i,j)/pc_r2d(i,j)
          tl_bp_r2d(i,j)=tl_br_r2d(i,j)/pc_r2d(i,j)
        END DO
      END DO
!
!=======================================================================
!  Iterate.
!=======================================================================
!
      ITERATE : DO it=1,Nbico(ng)-1

        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=p_r2d(i,j,it)
            tl_z1_r2d(i,j)=tl_p_r2d(i,j)
          END DO
        END DO

        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               z1_r2d)
        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               tl_z1_r2d)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      z1_r2d,                                     &
     &                      tl_z1_r2d)
# endif
!
!  Compute divergence NLM operator: v_r2d = div[h grad(z1_r2d)].
!
        Ltrans=.FALSE.
        CALL r2d_oper (ng, tile,                                        &
     &                 Ltrans, Lbck,                                    &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
# ifdef MASKING
     &                 umask, vmask, rmask,                             &
# endif
     &                 h,                                               &
     &                 pmon_u, pnom_v, pm, pn,                          &
     &                 pc_r2d, z1_r2d, v_r2d)
!
!  Compute divergence TLM operator: tl_v_r2d = div[h grad(tl_z1_r2d)].
!
        CALL tl_r2d_oper (ng, tile,                                     &
     &                    Lbck,                                         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
# ifdef MASKING
     &                    umask, vmask, rmask,                          &
# endif
     &                    h,                                            &
     &                    pmon_u, pnom_v, pm, pn,                       &
     &                    pc_r2d, tl_z1_r2d, tl_v_r2d)
!
!  Compute dot products and "tl_bc_ak" coefficient.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=r_r2d(i,j,it)/pc_r2d(i,j)
            z2_r2d(i,j)=br_r2d(i,j,it)
            z3_r2d(i,j)=bp_r2d(i,j,it)
            tl_z1_r2d(i,j)=tl_r_r2d(i,j)/pc_r2d(i,j)
            tl_z2_r2d(i,j)=tl_br_r2d(i,j)
            tl_z3_r2d(i,j)=tl_bp_r2d(i,j)
          END DO
        END DO

        CALL tl_r2d_dotp (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_zdf1,                                      &
# ifdef MASKING
     &                    rmask,                                        &
# endif
     &                    z2_r2d, z1_r2d, tl_z2_r2d, tl_z1_r2d)

        CALL tl_r2d_dotp (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_zdf2,                                      &
# ifdef MASKING
     &                    rmask,                                        &
# endif
     &                    z3_r2d, v_r2d, tl_z3_r2d, tl_v_r2d)

        tl_bc_ak=tl_zdf1/zdf2(it)-tl_zdf2*bc_ak(it)/zdf2(it)
!
!  Solve for new iterate of "tl_r2d_ref".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tl_r2d_ref(i,j)=tl_r2d_ref(i,j)+                            &
     &                      tl_bc_ak*p_r2d(i,j,it)+                     &
     &                      bc_ak(it)*tl_p_r2d(i,j)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute new (it+1) residual and direction vectors.
!-----------------------------------------------------------------------
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=bp_r2d(i,j,it)
            tl_z1_r2d(i,j)=tl_bp_r2d(i,j)
            tl_bv_r2d(i,j)=0.0_r8
          END DO
        END DO

        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               z1_r2d)
        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               tl_z1_r2d)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      z1_r2d)
# endif
!
!  Compute divergence operator, tranpose: bv_r2d = div[h grad(z1_r2d)]
!  for basic state and TLM (need to use "ad_r2d_oper").
!  Need to call "ad_r2d_bc" here since Ltrans is TRUE.
# ifdef DISTRIBUTE
!  Need to call "ad_mp_exchange" here since Ltrans is TRUE.
# endif
!
        Ltrans=.TRUE.
        CALL r2d_oper (ng, tile,                                        &
     &                 Ltrans, Lbck,                                    &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
# ifdef MASKING
     &                 umask, vmask, rmask,                             &
# endif
     &                 h,                                               &
     &                 pmon_u, pnom_v, pm, pn,                          &
     &                 pc_r2d, z1_r2d, bv_r2d)

        CALL ad_r2d_oper (ng, tile,                                     &
     &                    Lbck,                                         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
# ifdef MASKING
     &                    umask, vmask, rmask,                          &
# endif
     &                    h,                                            &
     &                    pmon_u, pnom_v, pm, pn,                       &
     &                    pc_r2d, tl_bv_r2d, tl_z1_r2d)

# ifdef DISTRIBUTE
        CALL ad_mp_exchange2d (ng, tile, model, 2,                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         bv_r2d,                                  &
     &                         tl_bv_r2d)
# endif
        CALL ad_r2d_bc (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bv_r2d)
        CALL ad_r2d_bc (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tl_bv_r2d)
!
!  Compute new residual vectors "tl_r" and "tl_br".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tl_r_r2d(i,j)=tl_r_r2d(i,j)-                                &
     &                    tl_bc_ak*v_r2d(i,j)-bc_ak(it)*tl_v_r2d(i,j)
            tl_br_r2d(i,j)=tl_br_r2d(i,j)-                              &
     &                     tl_bc_ak*bv_r2d(i,j)-bc_ak(it)*tl_bv_r2d(i,j)
          END DO
        END DO
!
!  Compute dot product and "tl_bc_ak" coefficient.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=r_r2d(i,j,it+1)/pc_r2d(i,j)
            z2_r2d(i,j)=br_r2d(i,j,it+1)
            tl_z1_r2d(i,j)=tl_r_r2d(i,j)/pc_r2d(i,j)
            tl_z2_r2d(i,j)=tl_br_r2d(i,j)
          END DO
        END DO

        CALL tl_r2d_dotp (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_zdf3,                                      &
# ifdef MASKING
     &                    rmask,                                        &
# endif
     &                    z1_r2d, z2_r2d, tl_z1_r2d, tl_z2_r2d)

        tl_bc_bk=tl_zdf3/zdf1(it)-tl_zdf1*bc_bk(it)/zdf1(it)
!
!  Use recurrence relationship to compute new direction vectors
!  "tl_p" and "tl_bp".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tl_p_r2d(i,j)=tl_r_r2d(i,j)/pc_r2d(i,j)+                    &
     &                    tl_bc_bk*p_r2d(i,j,it)+                       &
     &                    bc_bk(it)*tl_p_r2d(i,j)
            tl_bp_r2d(i,j)=tl_br_r2d(i,j)/pc_r2d(i,j)+                  &
     &                     tl_bc_bk*bp_r2d(i,j,it)+                     &
     &                     bc_bk(it)*tl_bp_r2d(i,j)
          END DO
        END DO

      END DO ITERATE

      RETURN
      END SUBROUTINE tl_biconj_tile

!
!***********************************************************************
      SUBROUTINE tl_r2d_oper (ng, tile,                                 &
     &                        Lbck,                                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        umask, vmask, rmask,                      &
# endif
     &                        h,                                        &
     &                        pmon_u, pnom_v, pm, pn,                   &
     &                        pc_r2d,                                   &
     &                        tl_r2d_in, tl_r2d_out)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lbck
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(inout) :: pc_r2d(LBi:,LBj:)
      real(r8), intent(inout) :: tl_r2d_in(LBi:,LBj:)

      real(r8), intent(out) :: tl_r2d_out(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_r2d_in(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: tl_r2d_out(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff, tl_fac

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Compute tl_r2d_out = div ( h grad(tl_r2d_in) ).
!----------------------------------------------------------------------
!
      cff=0.5_r8*g

      DO j=Jstr,Jend
        DO i=Istr,Iend+1
          tl_FX(i,j)=cff*pmon_u(i,j)*(h(i,j)+h(i-1,j))*                 &
     &               (tl_r2d_in(i,j)-tl_r2d_in(i-1,j))
# ifdef MASKING
          tl_FX(i,j)=tl_FX(i,j)*umask(i,j)
# endif
        END DO
      END DO

      DO j=Jstr,Jend+1
        DO i=Istr,Iend
          tl_FE(i,j)=cff*pnom_v(i,j)*(h(i,j)+h(i,j-1))*                 &
     &               (tl_r2d_in(i,j)-tl_r2d_in(i,j-1))
# ifdef MASKING
          tl_FE(i,j)=tl_FE(i,j)*vmask(i,j)
# endif
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
          tl_r2d_out(i,j)=pm(i,j)*pn(i,j)*                              &
     &                    (tl_FX(i+1,j)-tl_FX(i,j)+                     &
     &                     tl_FE(i,j+1)-tl_FE(i,j))
# ifdef MASKING
          tl_r2d_out(i,j)=tl_r2d_out(i,j)*rmask(i,j)
# endif
        END DO
      END DO

      RETURN
      END SUBROUTINE tl_r2d_oper

!
!***********************************************************************
      SUBROUTINE tl_r2d_dotp (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_DotProd,                               &
# ifdef MASKING
     &                        rmask,                                    &
# endif
     &                        s1_zeta, s2_zeta,                         &
     &                        tl_s1_zeta, tl_s2_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: s1_zeta(LBi:,LBj:)
      real(r8), intent(in) :: s2_zeta(LBi:,LBj:)
      real(r8), intent(in) :: tl_s1_zeta(LBi:,LBj:)
      real(r8), intent(in) :: tl_s2_zeta(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: s1_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s2_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_s1_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_s2_zeta(LBi:UBi,LBj:UBj)
# endif
!
      real(r8), intent(out) :: tl_DotProd
!
!  Local variable declarations.
!
      integer :: NSUB, i, j

      real(r8) :: tl_cff
      real(r8) :: tl_my_DotProd
# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute dot product between "tl_s1_zeta" and "tl_s2_zeta".
!-----------------------------------------------------------------------
!
      tl_my_DotProd=0.0_r8

      DO j=Jstr,Jend
        DO i=Istr,Iend
          tl_cff=tl_s1_zeta(i,j)*s2_zeta(i,j)+                          &
     &           s1_zeta(i,j)*tl_s2_zeta(i,j)
# ifdef MASKING
          tl_cff=tl_cff*rmask(i,j)
# endif
          tl_my_DotProd=tl_my_DotProd+tl_cff
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Perform parallel global reduction operations.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
# else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
# endif
!$OMP CRITICAL (TL_DOT_PROD)
      IF (tile_count.eq.0) THEN
        tl_DotProd=0.0_r8
      END IF
      tl_DotProd=tl_DotProd+tl_my_DotProd
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, tl_DotProd, op_handle)
# endif
      END IF
!$OMP END CRITICAL (TL_DOT_PROD)

      RETURN
      END SUBROUTINE tl_r2d_dotp

!
!***********************************************************************
      SUBROUTINE ad_biconj_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Lbck,                                  &
     &                           h, pmon_u, pnom_v, pm, pn,             &
# ifdef MASKING
     &                           umask, vmask, rmask,                   &
# endif
     &                           bc_ak, bc_bk, zdf1, zdf2, zdf3,        &
     &                           pc_r2d, r_r2d, br_r2d, p_r2d, bp_r2d,  &
     &                           ad_r2d_ref, ad_rhs_r2d)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_grid
      USE mod_parallel
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
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
      real(r8), intent(inout) :: ad_rhs_r2d(LBi:,LBj:)

      real(r8), intent(inout) :: ad_r2d_ref(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bc_ak(Nbico(ng))
      real(r8), intent(in) :: bc_bk(Nbico(ng))
      real(r8), intent(in) :: zdf1(Nbico(ng))
      real(r8), intent(in) :: zdf2(Nbico(ng))
      real(r8), intent(in) :: zdf3(Nbico(ng))
      real(r8), intent(inout) :: pc_r2d(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: r_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: br_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: p_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: ad_bp_r2d(LBi:UBi,LBj:UBj,Nbico(ng))
      real(r8), intent(inout) :: ad_rhs_r2d(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: ad_r2d_ref(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif
      logical :: Ltrans

      integer :: i, j, it, NSUB

      real(r8) :: ad_zdf1, ad_zdf2, ad_zdf3, ad_bc_ak, ad_bc_bk
      real(r8) :: my_ad_bc_ak, my_ad_bc_bk

      real(r8), dimension(LBi:UBi,LBj:UBj) :: bv_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: v_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: z1_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: z2_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: z3_r2d

      real(r8), dimension(LBi:UBi,LBj:UBj) :: ad_bp_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: ad_br_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: ad_bv_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: ad_p_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: ad_r_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: ad_v_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: ad_z1_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: ad_z2_r2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: ad_z3_r2d

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjoint of biconjugate gradient algorithm.
!-----------------------------------------------------------------------
!
!  Initialize local variables.
!
      DO j=LBj,UBj
        DO i=LBi,UBi
          bv_r2d(i,j)=0.0_r8
          v_r2d(i,j)=0.0_r8
          z1_r2d(i,j)=0.0_r8
          z2_r2d(i,j)=0.0_r8
          z3_r2d(i,j)=0.0_r8

          ad_bp_r2d(i,j)=0.0_r8
          ad_br_r2d(i,j)=0.0_r8
          ad_bv_r2d(i,j)=0.0_r8
          ad_p_r2d(i,j)=0.0_r8
          ad_r_r2d(i,j)=0.0_r8
          ad_v_r2d(i,j)=0.0_r8
          ad_z1_r2d(i,j)=0.0_r8
          ad_z2_r2d(i,j)=0.0_r8
          ad_z3_r2d(i,j)=0.0_r8
        END DO
      END DO
      ad_bc_ak=0.0_r8
      ad_bc_bk=0.0_r8
      ad_zdf1=0.0_r8
      ad_zdf2=0.0_r8
      ad_zdf3=0.0_r8
      my_ad_bc_ak=0.0_r8
      my_ad_bc_bk=0.0_r8
!
!=======================================================================
!  Iterate in reverse order.
!=======================================================================
!
      ITERATE : DO it=Nbico(ng)-1,1,-1
!
!  Compute "v_r2d" and "bv_r2d"
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=p_r2d(i,j,it)
            z2_r2d(i,j)=bp_r2d(i,j,it)
          END DO
        END DO

        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               z1_r2d)
        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               z2_r2d)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      z1_r2d,                                     &
     &                      z2_r2d)
# endif
!
!  Compute divergence operator: v_r2d = div[h grad(z1_r2d)].
!
        Ltrans=.FALSE.
        CALL r2d_oper (ng, tile,                                        &
     &                 Ltrans, Lbck,                                    &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
# ifdef MASKING
     &                 umask, vmask, rmask,                             &
# endif
     &                 h,                                               &
     &                 pmon_u, pnom_v, pm, pn,                          &
     &                 pc_r2d, z1_r2d, v_r2d)
!
!  Compute divergence operator, tranpose: bv_r2d = div[h grad(z2_r2d)].
!
        Ltrans=.TRUE.
        CALL r2d_oper (ng, tile,                                        &
     &                 Ltrans, Lbck,                                    &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
# ifdef MASKING
     &                 umask, vmask, rmask,                             &
# endif
     &                 h,                                               &
     &                 pmon_u, pnom_v, pm, pn,                          &
     &                 pc_r2d, z2_r2d, bv_r2d)

# ifdef DISTRIBUTE
        CALL ad_mp_exchange2d (ng, tile, model, 1,                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         bv_r2d)
# endif
        CALL ad_r2d_bc (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bv_r2d)
!
!  Adjoint of use recurrence relationship to compute new direction
!  vectors "tl_p" and "tl_bp".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_bp_r2d(i,j)=tl_br_r2d(i,j)/pc_r2d(i,j)+                  &
!>   &                     tl_bc_bk*bp_r2d(i,j,it)+                     &
!>   &                     bc_bk(it)*tl_bp_r2d(i,j)
!>
            my_ad_bc_bk=my_ad_bc_bk+bp_r2d(i,j,it)*ad_bp_r2d(i,j)
            ad_br_r2d(i,j)=ad_br_r2d(i,j)+ad_bp_r2d(i,j)/pc_r2d(i,j)
            ad_bp_r2d(i,j)=bc_bk(it)*ad_bp_r2d(i,j)

!>          tl_p_r2d(i,j)=tl_r_r2d(i,j)/pc_r2d(i,j)+                    &
!>   &                    tl_bc_bk*p_r2d(i,j,it)+                       &
!>   &                    bc_bk(it)*tl_p_r2d(i,j)
!>
            my_ad_bc_bk=my_ad_bc_bk+p_r2d(i,j,it)*ad_p_r2d(i,j)
            ad_r_r2d(i,j)=ad_r_r2d(i,j)+ad_p_r2d(i,j)/pc_r2d(i,j)
            ad_p_r2d(i,j)=bc_bk(it)*ad_p_r2d(i,j)
          END DO
        END DO
!
!  Perform parallel global reduction operation on "my_ad_bc_bk".
!
# ifdef DISTRIBUTE
        NSUB=1                              ! distributed-memory
# else
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                            ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)        ! tiled application
        END IF
# endif
!$OMP CRITICAL (AD_DOT_PROD1)
        IF (tile_count.eq.0) THEN
          ad_bc_bk=0.0_r8
        END IF
        ad_bc_bk=ad_bc_bk+my_ad_bc_bk
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
# ifdef DISTRIBUTE
          op_handle='SUM'
          CALL mp_reduce (ng, model, 1, ad_bc_bk, op_handle)
# endif
        END IF
        my_ad_bc_bk=0.0_r8
!$OMP END CRITICAL (AD_DOT_PROD1)
!
!  Adjoint of compute dot product and "tl_bc_ak" coefficient.
!
!>      tl_bc_bk=tl_zdf3/zdf1(it)-                                      &
!>   &           tl_zdf1*bc_bk(it)/zdf1(it)
!>
        ad_zdf1=ad_zdf1-ad_bc_bk*bc_bk(it)/zdf1(it)
        ad_zdf3=ad_zdf3+ad_bc_bk/zdf1(it)
        ad_bc_bk=0.0

        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=r_r2d(i,j,it+1)/pc_r2d(i,j)
            z2_r2d(i,j)=br_r2d(i,j,it+1)
          END DO
        END DO

        CALL ad_r2d_dotp (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    ad_zdf3,                                      &
# ifdef MASKING
     &                    rmask,                                        &
# endif
     &                    z1_r2d, z2_r2d, ad_z1_r2d, ad_z2_r2d)

        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_z2_r2d(i,j)=tl_br_r2d(i,j)
!>
            ad_br_r2d(i,j)=ad_br_r2d(i,j)+ad_z2_r2d(i,j)
            ad_z2_r2d(i,j)=0.0_r8
!>          tl_z1_r2d(i,j)=tl_r_r2d(i,j)/pc_r2d(i,j)
!>
            ad_r_r2d(i,j)=ad_r_r2d(i,j)+ad_z1_r2d(i,j)/pc_r2d(i,j)
            ad_z1_r2d(i,j)=0.0_r8
          END DO
        END DO
!
!  Adjoint of compute new residual vectors "tl_r" and "tl_br".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_br_r2d(i,j)=tl_br_r2d(i,j)-tl_bc_ak*bv_r2d(i,j)-         &
!>   &                     bc_ak(it)*tl_bv_r2d(i,j)
!>
            ad_bv_r2d(i,j)=ad_bv_r2d(i,j)-bc_ak(it)*ad_br_r2d(i,j)
            my_ad_bc_ak=my_ad_bc_ak-bv_r2d(i,j)*ad_br_r2d(i,j)
!>          tl_r_r2d(i,j)=tl_r_r2d(i,j)-tl_bc_ak*v_r2d(i,j)-            &
!>   &                    bc_ak(it)*tl_v_r2d(i,j)
!>
            ad_v_r2d(i,j)=ad_v_r2d(i,j)-bc_ak(it)*ad_r_r2d(i,j)
            my_ad_bc_ak=my_ad_bc_ak-v_r2d(i,j)*ad_r_r2d(i,j)
          END DO
        END DO
!
!  Adjoint of compute divergence operator, tranpose:
!  tl_bv_r2d = div[h grad(tl_z1_r2d)]. Need to use "tl_r2d_oper" here.
!
        CALL r2d_bc (ng, tile,                                          &
     &               LBi, UBi, LBj, UBj,                                &
     &               ad_bv_r2d)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      ad_bv_r2d)
# endif
        CALL tl_r2d_oper (ng, tile,                                     &
     &                    Lbck,                                         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
# ifdef MASKING
     &                    umask, vmask, rmask,                          &
# endif
     &                    h,                                            &
     &                    pmon_u, pnom_v, pm, pn,                       &
     &                    pc_r2d, ad_bv_r2d, ad_z1_r2d)
!
!  Adjoint of loading "tl_z1_r2d".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_z1_r2d(i,j)=tl_bp_r2d(i,j)
!>
            ad_bp_r2d(i,j)=ad_bp_r2d(i,j)+ad_z1_r2d(i,j)
            ad_z1_r2d(i,j)=0.0_r8
            ad_bv_r2d(i,j)=0.0_r8     ! cleared because it was used in
                                      ! "tl_r2d_oper"
          END DO
        END DO
!
!  Adjoint of solve for new iterate "tl_r2d_ref".
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_r2d_ref(i,j)=tl_r2d_ref(i,j)+tl_bc_ak*p_r2d(i,j,it)+     &
!>   &                      bc_ak(it)*tl_p_r2d(i,j)
!>
            ad_p_r2d(i,j)=ad_p_r2d(i,j)+bc_ak(it)*ad_r2d_ref(i,j)
            my_ad_bc_ak=my_ad_bc_ak+p_r2d(i,j,it)*ad_r2d_ref(i,j)
          END DO
        END DO
!
!  Perform parallel global reduction operation on "my_ad_bc_ak".
!
# ifdef DISTRIBUTE
        NSUB=1                              ! distributed-memory
# else
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                            ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)        ! tiled application
        END IF
# endif
!$OMP CRITICAL (AD_DOT_PROD2)
        IF (tile_count.eq.0) THEN
          ad_bc_ak=0.0_r8
        END IF
        ad_bc_ak=ad_bc_ak+my_ad_bc_ak
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
# ifdef DISTRIBUTE
          op_handle='SUM'
          CALL mp_reduce (ng, model, 1, ad_bc_ak, op_handle)
# endif
        END IF
        my_ad_bc_ak=0.0_r8
!$OMP END CRITICAL (AD_DOT_PROD2)
!
!  Adjoint of compute dot products and "tl_bc_ak" coefficient.
!
!>      tl_bc_ak=tl_zdf1/zdf2(it)-                                      &
!>   &           tl_zdf2*bc_ak(it)/zdf2(it)
!>
        ad_zdf2=ad_zdf2-ad_bc_ak*bc_ak(it)/zdf2(it)
        ad_zdf1=ad_zdf1+ad_bc_ak/zdf2(it)
        ad_bc_ak=0.0

        DO j=Jstr,Jend
          DO i=Istr,Iend
            z1_r2d(i,j)=r_r2d(i,j,it)/pc_r2d(i,j)
            z2_r2d(i,j)=br_r2d(i,j,it)
            z3_r2d(i,j)=bp_r2d(i,j,it)
          END DO
        END DO

        CALL ad_r2d_dotp (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    ad_zdf2,                                      &
# ifdef MASKING
     &                    rmask,                                        &
# endif
     &                    z3_r2d, v_r2d, ad_z3_r2d, ad_v_r2d)

        CALL ad_r2d_dotp (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    ad_zdf1,                                      &
# ifdef MASKING
     &                    rmask,                                        &
# endif
     &                    z2_r2d, z1_r2d, ad_z2_r2d, ad_z1_r2d)

        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_z1_r2d(i,j)=tl_r_r2d(i,j)/pc_r2d(i,j)
!>
            ad_r_r2d(i,j)=ad_r_r2d(i,j)+ad_z1_r2d(i,j)/pc_r2d(i,j)
            ad_z1_r2d(i,j)=0.0_r8
!>          tl_z2_r2d(i,j)=tl_br_r2d(i,j)
!>
            ad_br_r2d(i,j)=ad_br_r2d(i,j)+ad_z2_r2d(i,j)
            ad_z2_r2d(i,j)=0.0_r8
!>          tl_z3_r2d(i,j)=tl_bp_r2d(i,j)
!>
            ad_bp_r2d(i,j)=ad_bp_r2d(i,j)+ad_z3_r2d(i,j)
            ad_z3_r2d(i,j)=0.0_r8
          END DO
        END DO
!
!  Adjoint of compute divergence operator:
!    tl_v_r2d = div[h grad(tl_z1_r2d)].
!
        CALL ad_r2d_oper (ng, tile,                                     &
     &                    Lbck,                                         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
# ifdef MASKING
     &                    umask, vmask, rmask,                          &
# endif
     &                    h,                                            &
     &                    pmon_u, pnom_v, pm, pn,                       &
     &                    pc_r2d, ad_z1_r2d, ad_v_r2d)

# ifdef DISTRIBUTE
        CALL ad_mp_exchange2d (ng, tile, model, 1,                      &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         NghostPoints,                            &
     &                         EWperiodic(ng), NSperiodic(ng),          &
     &                         ad_z1_r2d)
# endif
        CALL ad_r2d_bc (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  ad_z1_r2d)

        DO j=Jstr,Jend
          DO i=Istr,Iend
!>          tl_z1_r2d(i,j)=tl_p_r2d(i,j)
!>
            ad_p_r2d(i,j)=ad_p_r2d(i,j)+ad_z1_r2d(i,j)
            ad_z1_r2d(i,j)=0.0_r8
          END DO
        END DO

      END DO ITERATE
!
!  Adjoint of set the initial values for residual vectors "tl_r" and
!  "tl_br". Then, use recurrence relationship to compute direction
!  vectors "tl_p" and "tl_bp".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!>        tl_p_r2d(i,j)=tl_r_r2d(i,j)/pc_r2d(i,j)
!>
          ad_r_r2d(i,j)=ad_r_r2d(i,j)+ad_p_r2d(i,j)/pc_r2d(i,j)
          ad_p_r2d(i,j)=0.0
!>        tl_bp_r2d(i,j)=tl_br_r2d(i,j)/pc_r2d(i,j)
!>
          ad_br_r2d(i,j)=ad_br_r2d(i,j)+ad_bp_r2d(i,j)/pc_r2d(i,j)
          ad_bp_r2d(i,j)=0.0
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
!>        tl_br_r2d(i,j)=tl_r_r2d(i,j)
!>
          ad_r_r2d(i,j)=ad_r_r2d(i,j)+ad_br_r2d(i,j)
          ad_br_r2d(i,j)=0.0
!>        tl_r_r2d(i,j)=tl_rhs_r2d(i,j)-tl_z1_r2d(i,j)
!>
          ad_rhs_r2d(i,j)=ad_rhs_r2d(i,j)+ad_r_r2d(i,j)
          ad_z1_r2d(i,j)=ad_z1_r2d(i,j)-ad_r_r2d(i,j)
          ad_r_r2d(i,j)=0.0_r8
        END DO
      END DO
!
!  Adjoint of compute starting value of divergence TLM operator:
!    tl_z1_r2d = div[h grad(tl_r2d_ref)].
!
      CALL ad_r2d_oper (ng, tile,                                       &
     &                  Lbck,                                           &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
# ifdef MASKING
     &                  umask, vmask, rmask,                            &
# endif
     &                  h,                                              &
     &                  pmon_u, pnom_v, pm, pn,                         &
     &                  pc_r2d, ad_r2d_ref, ad_z1_r2d)

# ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_r2d_ref)
# endif
      CALL ad_r2d_bc (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                ad_r2d_ref)

      RETURN
      END SUBROUTINE ad_biconj_tile

!
!***********************************************************************
      SUBROUTINE ad_r2d_oper (ng, tile,                                 &
     &                        Lbck,                                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        umask, vmask, rmask,                      &
# endif
     &                        h,                                        &
     &                        pmon_u, pnom_v, pm, pn,                   &
     &                        pc_r2d,                                   &
     &                        ad_r2d_in, ad_r2d_out)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lbck
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pmon_u(LBi:,LBj:)
      real(r8), intent(in) :: pnom_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(inout) :: pc_r2d(LBi:,LBj:)
      real(r8), intent(inout) :: ad_r2d_in(LBi:,LBj:)

      real(r8), intent(inout) :: ad_r2d_out(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_r2d_in(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: ad_r2d_out(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: adfac, cff

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX

# include "set_bounds.h"
!
!----------------------------------------------------------------------
!  Adjoint of compute tl_r2d_out = div( h grad(ad_r2d_in) ).
!----------------------------------------------------------------------
!
!  Initialize local arrays.
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
          ad_FX(i,j)=0.0_r8
          ad_FE(i,j)=0.0_r8
        END DO
      END DO
!
!  Compute adjoint of divergence XI- and ETA-components.
!
      cff=0.5_r8*g

      DO j=Jstr,Jend
        DO i=Istr,Iend
# ifdef MASKING
!>        tl_r2d_out(i,j)=tl_r2d_out(i,j)*rmask(i,j)
!>
          ad_r2d_out(i,j)=ad_r2d_out(i,j)*rmask(i,j)
# endif
!>        tl_r2d_out(i,j)=pm(i,j)*pn(i,j)*                              &
!>   &                    (tl_FX(i+1,j)-tl_FX(i,j)+                     &
!>   &                     tl_FE(i,j+1)-tl_FE(i,j))
!>
          adfac=pm(i,j)*pn(i,j)*ad_r2d_out(i,j)
          ad_FX(i  ,j)=ad_FX(i  ,j)-adfac
          ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
          ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
          ad_FE(i,j  )=ad_FE(i,j  )-adfac
          ad_r2d_out(i,j)=0.0_r8
        END DO
      END DO

      DO j=Jstr,Jend+1
        DO i=Istr,Iend
# ifdef MASKING
!>        tl_FE(i,j)=tl_FE(i,j)*vmask(i,j)
!>
          ad_FE(i,j)=ad_FE(i,j)*vmask(i,j)
# endif
!>        tl_FE(i,j)=cff*pnom_v(i,j)*(h(i,j)+h(i,j-1))*                 &
!>   &               (tl_r2d_in(i,j)-tl_r2d_in(i,j-1))
!>
          adfac=cff*pnom_v(i,j)*(h(i,j)+h(i,j-1))*ad_FE(i,j)
          ad_r2d_in(i,j-1)=ad_r2d_in(i,j-1)-adfac
          ad_r2d_in(i,j  )=ad_r2d_in(i,j  )+adfac
          ad_FE(i,j)=0.0_r8
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend+1
# ifdef MASKING
!>        tl_FX(i,j)=tl_FX(i,j)*umask(i,j)
!>
          ad_FX(i,j)=ad_FX(i,j)*umask(i,j)
# endif
!>        tl_FX(i,j)=cff*pmon_u(i,j)*(h(i,j)+h(i-1,j))*                 &
!>   &               (tl_r2d_in(i,j)-tl_r2d_in(i-1,j))
!>
          adfac=cff*pmon_u(i,j)*(h(i,j)+h(i-1,j))*ad_FX(i,j)
          ad_r2d_in(i-1,j)=ad_r2d_in(i-1,j)-adfac
          ad_r2d_in(i  ,j)=ad_r2d_in(i  ,j)+adfac
          ad_FX(i,j)=0.0_r8
        END DO
      END DO

      RETURN
      END SUBROUTINE ad_r2d_oper

!
!***********************************************************************
      SUBROUTINE ad_r2d_bc (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      ad_A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_r2d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ad_A(LBi:,LBj:)
# else
      real(r8), intent(inout) :: ad_A(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: adfac

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set adjoint periodic boundary conditons.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_A)
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
!>        tl_A(Iend+1,Jend+1)=0.5_r8*(tl_A(Iend+1,Jend  )+              &
!>   &                                tl_A(Iend  ,Jend+1))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jend+1)
          ad_A(Iend+1,Jend  )=ad_A(Iend+1,Jend  )+adfac
          ad_A(Iend  ,Jend+1)=ad_A(Iend  ,Jend+1)+adfac
          ad_A(Iend+1,Jend+1)=0.0_r8
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
!>        tl_A(Istr-1,Jend+1)=0.5_r8*(tl_A(Istr-1,Jend  )+              &
!>   &                                tl_A(Istr  ,Jend+1))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jend+1)
          ad_A(Istr-1,Jend  )=ad_A(Istr-1,Jend  )+adfac
          ad_A(Istr  ,Jend+1)=ad_A(Istr  ,Jend+1)+adfac
          ad_A(Istr-1,Jend+1)=0.0_r8
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
!>        tl_A(Iend+1,Jstr-1)=0.5_r8*(tl_A(Iend  ,Jstr-1)+              &
!>   &                                tl_A(Iend+1,Jstr  ))
!>
          adfac=0.5_r8*ad_A(Iend+1,Jstr-1)
          ad_A(Iend  ,Jstr-1)=ad_A(Iend  ,Jstr-1)+adfac
          ad_A(Iend+1,Jstr  )=ad_A(Iend+1,Jstr  )+adfac
          ad_A(Iend+1,Jstr-1)=0.0_r8
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
!>        tl_A(Istr-1,Jstr-1)=0.5_r8*(tl_A(Istr  ,Jstr-1)+              &
!>   &                                tl_A(Istr-1,Jstr  ))
!>
          adfac=0.5_r8*ad_A(Istr-1,Jstr-1)
          ad_A(Istr  ,Jstr-1)=ad_A(Istr  ,Jstr-1)+adfac
          ad_A(Istr-1,Jstr  )=ad_A(Istr-1,Jstr  )+adfac
          ad_A(Istr-1,Jstr-1)=0.0_r8
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint North-South gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
# ifdef R2D_GRADIENT
!>          tl_A(i,Jstr-1)=tl_A(i,Jstr)
!>
            ad_A(i,Jstr  )=ad_A(i,Jstr)+ad_A(i,Jstr-1)
            ad_A(i,Jstr-1)=0.0_r8
# else
!>          tl_A(i,Jstr-1)=0.0_r8
!>
            ad_A(i,Jstr-1)=0.0_r8
# endif
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
# ifdef R2D_GRADIENT
!>          tl_A(i,Jend+1)=tl_A(i,Jend)
!>
            ad_A(i,Jend  )=ad_A(i,Jend)+ad_A(i,Jend+1)
            ad_A(i,Jend+1)=0.0_r8
# else
!>          tl_A(i,Jend+1)=0.0_r8
!>
            ad_A(i,Jend+1)=0.0_r8
# endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Adjoint East-West gradient boundary conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
# ifdef R2D_GRADIENT
!>          tl_A(Istr-1,j)=tl_A(Istr,j)
!>
            ad_A(Istr  ,j)=ad_A(Istr,j)+ad_A(Istr-1,j)
            ad_A(Istr-1,j)=0.0_r8
# else
!>          tl_A(Istr-1,j)=0.0_r8
!>
            ad_A(Istr-1,j)=0.0_r8
# endif
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
# ifdef R2D_GRADIENT
!>          tl_A(Iend+1,j)=tl_A(Iend,j)
!>
            ad_A(Iend  ,j)=ad_A(Iend,j)+ad_A(Iend+1,j)
            ad_A(Iend+1,j)=0.0_r8
# else
!>          tl_A(Iend+1,j)=0.0_r8
!>
            ad_A(Iend+1,j)=0.0_r8
# endif
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_r2d_bc

!
!***********************************************************************
      SUBROUTINE ad_r2d_dotp (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ad_DotProd,                               &
# ifdef MASKING
     &                        rmask,                                    &
# endif
     &                        s1_zeta, s2_zeta,                         &
     &                        ad_s1_zeta, ad_s2_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: s1_zeta(LBi:,LBj:)
      real(r8), intent(in) :: s2_zeta(LBi:,LBj:)
      real(r8), intent(inout) :: ad_s1_zeta(LBi:,LBj:)
      real(r8), intent(inout) :: ad_s2_zeta(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: s1_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: s2_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_s1_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_s2_zeta(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: ad_DotProd
!
!  Local variable declarations.
!
      integer :: NSUB, i, j

      real(r8) :: ad_cff
      real(r8) :: ad_my_DotProd
# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute adjoint dot product between ad_s1_zeta and ad_s2_zeta.
!-----------------------------------------------------------------------
!
      ad_my_DotProd=ad_DotProd
      ad_DotProd=0.0_r8
      ad_cff=0.0_r8
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!>        tl_my_DotProd=tl_my_DotProd+tl_cff
!>
          ad_cff=ad_cff+ad_my_DotProd
# ifdef MASKING
!>        tl_cff=tl_cff*rmask(i,j)
!>
          ad_cff=ad_cff*rmask(i,j)
# endif
!>        tl_cff=tl_s1_zeta(i,j)*s2_zeta(i,j)+                          &
!>   &           s1_zeta(i,j)*tl_s2_zeta(i,j)
!>
          ad_s1_zeta(i,j)=ad_s1_zeta(i,j)+ad_cff*s2_zeta(i,j)
          ad_s2_zeta(i,j)=ad_s2_zeta(i,j)+ad_cff*s1_zeta(i,j)
          ad_cff=0.0_r8
        END DO
      END DO
!>    tl_my_DotProd=0.0_r8
!>
      ad_my_DotProd=0.0_r8

      RETURN
      END SUBROUTINE ad_r2d_dotp
#endif
      END MODULE zeta_balance_mod
