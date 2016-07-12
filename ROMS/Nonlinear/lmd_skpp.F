#include "cppdefs.h"
#undef LMD_BOUND
#undef QUADRATIC
#define SASHA
      MODULE lmd_skpp_mod
#if defined NONLINEAR && defined LMD_SKPP && defined SOLVE3D
!
!svn $Id: lmd_skpp.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine determines the depth of surface  oceanic boundary      !
!  layer, hsbl, as the shallowest depth where the bulk Richardson      !
!  number is equal to the critical value, Ric.                         !
!                                                                      !
!  Then,  it computes the vertical mixing coefficients  within the     !
!  boundary layer. They depend on surface forcing and the magnitude    !
!  and gradient of interior mixing below  the boundary layer.  The     !
!  ocean interior is allowed to force the boundary layer through a     !
!  dependence of the nondimensional vertical shape function G(sigma)   !
!  and its vertical derivative at  sigma=1  on the interior  mixing    !
!  coefficients, and it vertical derivative at d=hsbl. The boundary    !
!  layer mixing coefficients are computed by matching these values.    !
!                                                                      !
! Reference:                                                           !
!                                                                      !
!  Large, W.G., J.C. McWilliams, and S.C. Doney, 1994: A Review        !
!    and model with a nonlocal boundary layer parameterization,        !
!    Reviews of Geophysics, 32,363-403.                                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: lmd_skpp

      CONTAINS
!
!***********************************************************************
      SUBROUTINE lmd_skpp (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
      USE mod_clima
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL lmd_skpp_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng),                                     &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    GRID(ng) % f,                                 &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_r,                               &
     &                    GRID(ng) % z_w,                               &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v,                                &
     &                    OCEAN(ng) % pden,                             &
     &                    FORCES(ng) % srflx,                           &
     &                    FORCES(ng) % stflx,                           &
     &                    FORCES(ng) % bustr,                           &
     &                    FORCES(ng) % bvstr,                           &
     &                    FORCES(ng) % sustr,                           &
     &                    FORCES(ng) % svstr,                           &
     &                    MIXING(ng) % alpha,                           &
# ifdef SALINITY
     &                    MIXING(ng) % beta,                            &
# endif
     &                    MIXING(ng) % bvf,                             &
# ifdef LMD_NONLOCAL
     &                    MIXING(ng) % ghats,                           &
# endif
     &                    MIXING(ng) % Akt,                             &
     &                    MIXING(ng) % Akv,                             &
#ifdef AKT_3D
     &                    CLIMA(ng) % tclm,                             &
#endif
     &                    MIXING(ng) % hsbl,                            &
     &                    MIXING(ng) % ksbl)
      RETURN
      END SUBROUTINE lmd_skpp
!
!***********************************************************************
      SUBROUTINE lmd_skpp_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp,                                   &
# ifdef MASKING
     &                          rmask,                                  &
# endif
     &                          f, Hz, z_r, z_w,                        &
     &                          u, v, pden,                             &
     &                          srflx, stflx,                           &
     &                          bustr, bvstr, sustr, svstr,             &
     &                          alpha,                                  &
# ifdef SALINITY
     &                          beta,                                   &
# endif
     &                          bvf,                                    &
# ifdef LMD_NONLOCAL
     &                          ghats,                                  &
# endif
     &                          Akt, Akv,                               &
#ifdef AKT_3D
     &                          tclm,                                   &
#endif
     &                          hsbl, ksbl)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
# ifdef LMD_SHAPIRO
      USE shapiro_mod
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: pden(LBi:,LBj:,:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(in) :: alpha(LBi:,LBj:)
#  ifdef SALINITY
      real(r8), intent(in) :: beta(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
#ifdef AKT_3D
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: hsbl(LBi:,LBj:)
      integer,  intent(out) :: ksbl(LBi:,LBj:)
#  ifdef LMD_NONLOCAL
      real(r8), intent(out) :: ghats(LBi:,LBj:,0:,:)
#  endif
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: pden(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: alpha(LBi:UBi,LBj:UBj)
#  ifdef SALINITY
      real(r8), intent(in) :: beta(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
#ifdef AKT_3D
      real(r8), intent(in) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng)+2)
# endif
      real(r8), intent(inout) :: hsbl(LBi:UBi,LBj:UBj)
      integer,  intent(out) :: ksbl(LBi:UBi,LBj:UBj)
#  ifdef LMD_NONLOCAL
      real(r8), intent(out) :: ghats(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8
      real(r8), parameter :: small = 1.0E-20_r8

      real(r8) :: Gm, Gt, Gs, K_bl, Ribot, Ritop, Rk
      real(r8) :: Uk, Ustarb, Ustar3, Vk, Vtc
      real(r8) :: a1, a2, a3, cff, cff1, cff2,cff_up, cff_dn
      real(r8) :: depth, dK_bl, hekman, hmonob, sigma, zbl
      real(r8) :: zetahat, zetapar
# ifdef QUADRATIC
      real(r8) :: slope_up, slope_down
      real(r8) :: a_co, b_co, c_co, z_up, sqrt_arg
# endif

      real(r8), dimension (IminS:ImaxS) :: Rref
      real(r8), dimension (IminS:ImaxS) :: Uref
      real(r8), dimension (IminS:ImaxS) :: Vref

      real(r8), dimension (IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: Bflux

      real(r8), dimension (IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension (IminS:ImaxS,0:N(ng)) :: dR
      real(r8), dimension (IminS:ImaxS,0:N(ng)) :: dU
      real(r8), dimension (IminS:ImaxS,0:N(ng)) :: dV

      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Bo
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Bosol
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Bfsfc
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Gm1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Gt1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Gs1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Ustar
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: dGm1dS
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: dGt1dS
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: dGs1dS
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: f1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: sl_dpth
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: swdk
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: wm
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: ws
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: zgrid

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize relevant parameters.
!-----------------------------------------------------------------------
!
      Vtc=lmd_Cv*SQRT(-lmd_betaT)/(SQRT(lmd_cs*lmd_epsilon)*            &
     &                             lmd_Ric*vonKar*vonKar)
!
!-----------------------------------------------------------------------
!  Get approximation of surface layer depth using "lmd_eps" and
!  boundary layer depth from previous time step.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          sl_dpth(i,j)=lmd_epsilon*(z_w(i,j,N(ng))-hsbl(i,j))
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute turbulent friction velocity (m/s) "Ustar" from wind stress
!  at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Ustar(i,j)=SQRT(SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+   &
     &                         (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2))
# ifdef MASKING
          Ustar(i,j)=Ustar(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute surface turbulent buoyancy forcing "Bo" (m2/s3). Remove
!  incoming solar shortwave radiation because this contribution is
!  included in "Bosol".  Compute surface radiative buoyancy forcing
!  "Bosol" (m2/s3).
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
# ifdef SALINITY
          Bo(i,j)=g*(alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))-          &
     &               beta (i,j)*stflx(i,j,isalt))
# else
          Bo(i,j)=g*alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))
# endif
          Bosol(i,j)=g*alpha(i,j)*srflx(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux (m2/s3) at W-points.  Notice that the
!  radiative bouyancy flux is distributed vertically using decay
!  function, swdk. Begin computation of nonlocal transport, "ghats".
!-----------------------------------------------------------------------
!
      DO k=0,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            zgrid(i,j)=z_w(i,j,N(ng))-z_w(i,j,k)
          END DO
        END DO
        CALL lmd_swfrac_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        -1.0_r8, zgrid, swdk)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Bflux(i,j,k)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
# ifdef MASKING
            Bflux(i,j,k)=Bflux(i,j,k)*rmask(i,j)
# endif
# ifdef LMD_NONLOCAL
            cff=1.0_r8-(0.5_r8+SIGN(0.5_r8,Bflux(i,j,k)))
            ghats(i,j,k,itemp)=-cff*(stflx(i,j,itemp)-srflx(i,j)+       &
     &                               srflx(i,j)*(1.0_r8-swdk(i,j)))
#  ifdef SALINITY
            ghats(i,j,k,isalt)=cff*stflx(i,j,isalt)
#  endif
# endif
          END DO
        END DO
      END DO
!
!=======================================================================
!  Compute bulk Richardson number "Rib" and then find depth of the
!  oceanic surface boundary layer "hsbl", such that Rib(hsbl)=Ric.
!=======================================================================
!
      DO j=Jstr,Jend
# ifdef RI_SPLINES
!
! Construct parabolic splines for vertical derivatives of potential
! density and velocity components at W-points.  FC is a scratch array.
!
        DO i=Istr,Iend
          FC(i,0)=0.0_r8
          dR(i,0)=0.0_r8
          dU(i,0)=0.0_r8
          dV(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                             &
     &                  Hz(i,j,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*Hz(i,j,k+1)
            dR(i,k)=cff*(6.0_r8*(pden(i,j,k+1)-pden(i,j,k))-            &
     &                   Hz(i,j,k)*dR(i,k-1))
            dU(i,k)=cff*(3.0_r8*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp)+     &
     &                           u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))-    &
     &                   Hz(i,j,k)*dU(i,k-1))
            dV(i,k)=cff*(3.0_r8*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+     &
     &                           v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp))-    &
     &                   Hz(i,j,k)*dV(i,k-1))
          END DO
        END DO
        DO i=Istr,Iend
          dR(i,N(ng))=0.0_r8
          dU(i,N(ng))=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            dR(i,k)=dR(i,k)-FC(i,k)*dR(i,k+1)
            dU(i,k)=dU(i,k)-FC(i,k)*dU(i,k+1)
            dV(i,k)=dV(i,k)-FC(i,k)*dV(i,k+1)
          END DO
        END DO
# else
!
! Compute vertical derivatives of potential density and velocity
! components at W-points.
!
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            dR(i,k)=cff*(pden(i,j,k+1)-pden(i,j,k))
            cff=0.5_r8*cff
            dU(i,k)=cff*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp)+             &
     &                   u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))
            dV(i,k)=cff*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+             &
     &                   v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp))
          END DO
        END DO
        DO i=Istr,Iend
          dR(i,0)=0.0_r8
          dR(i,N(ng))=0.0_r8
          dU(i,0)=0.0_r8
          dU(i,N(ng))=0.0_r8
          dV(i,0)=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute bulk Richardson number "Rib" and then find depth of oceanic
!  surface boundary layer "hsbl".
!
!                  [Br - B(d)] * d
!     Rib(d) = ----------------------- ;     Rib(hsbl)=Ric     (1)
!              |Vr - V(d)|^2 + Vt(d)^2
!
!  where "Br" and "Vr" are the surface reference buoyancy and velocity
!  while "B(d)" and "V(d)" are the bouyancy and velocity at depth "d".
!
!  In the code below, the criterion "Rib(hsbl)=Ric" is reformulated
!  as follows:
!
!     Rib(d)       Ritop(d)
!     ------ = --------------- = 1                             (2)
!      Ric      Ric * Ribot(d)
!
!  where "Ritop" and "Ribot" are numerator and denominator in Eq. (1).
!  In its turn, Eq. (2) is rewritten in the following form:
!
!     FC(d) = Ritop(d) - Ric * Ribot(d) = 0                    (3)
!
!  That is, the planetary boundary layer extends to the depth where
!  the critical function "FC(d)" changes its sign.
!-----------------------------------------------------------------------
!
!  Compute potential density and velocity components surface reference
!  values.
!
        cff1=1.0_r8/3.0_r8
        cff2=1.0_r8/6.0_r8
        DO i=Istr,Iend
          Rref(i)=pden(i,j,N(ng))+                                      &
     &            Hz(i,j,N(ng))*(cff1*dR(i,N(ng))+cff2*dR(i,N(ng)-1))
          Uref(i)=0.5_r8*(u(i,j,N(ng),nstp)+u(i+1,j,N(ng),nstp))+       &
     &            Hz(i,j,N(ng))*(cff1*dU(i,N(ng))+cff2*dU(i,N(ng)-1))
          Vref(i)=0.5_r8*(v(i,j,N(ng),nstp)+v(i,j+1,N(ng),nstp))+       &
     &            Hz(i,j,N(ng))*(cff1*dV(i,N(ng))+cff2*dV(i,N(ng)-1))
        END DO
!
!  Compute turbulent velocity scales for momentum (wm) and tracers (ws).
!  Then, compute critical function (FC) for bulk Richardson number.
!
        DO i=Istr,Iend
          FC(i,N(ng))=0.0_r8
          DO k=N(ng),1,-1
            depth=z_w(i,j,N(ng))-z_w(i,j,k-1)
            IF (Bflux(i,j,k-1).lt.0.0_r8) THEN
              sigma=MIN(sl_dpth(i,j),depth)
            ELSE
              sigma=depth
            END IF
            Ustar3=Ustar(i,j)*Ustar(i,j)*Ustar(i,j)
            zetahat=vonKar*sigma*Bflux(i,j,k-1)
            zetapar=zetahat/(Ustar3+small)
            IF (zetahat.ge.0.0_r8) THEN                         ! stable
              wm(i,j)=vonKar*Ustar(i,j)/(1.0_r8+5.0_r8*zetapar)
              ws(i,j)=wm(i,j)
            ELSE                                              ! unstable
              IF (zetapar.gt.lmd_zetam) THEN
                wm(i,j)=vonKar*Ustar(i,j)*                              &
     &                  (1.0_r8-16.0_r8*zetapar)**0.25_r8
              ELSE
                wm(i,j)=vonKar*(lmd_am*Ustar3-lmd_cm*zetahat)**r3
              END IF
              IF (zetapar.gt.lmd_zetas) THEN
                ws(i,j)=vonKar*Ustar(i,j)*                              &
     &                  (1.0_r8-16.0_r8*zetapar)**0.5_r8
              ELSE
                ws(i,j)=vonKar*(lmd_as*Ustar3-lmd_cs*zetahat)**r3
              END IF
            END IF
!
            Rk=pden(i,j,k)-                                             &
     &         Hz(i,j,k)*(cff1*dR(i,k-1)+cff2*dR(i,k))
            Uk=0.5_r8*(u(i,j,k,nstp)+u(i+1,j,k,nstp))-                  &
     &         Hz(i,j,k)*(cff1*dU(i,k-1)+cff2*dU(i,k))
            Vk=0.5_r8*(v(i,j,k,nstp)+v(i,j+1,k,nstp))-                  &
     &         Hz(i,j,k)*(cff1*dV(i,k-1)+cff2*dV(i,k))
!
            Ritop=-gorho0*(Rref(i)-Rk)*depth
            Ribot=(Uref(i)-Uk)**2+(Vref(i)-Vk)**2+                      &
     &            Vtc*depth*ws(i,j)*SQRT(ABS(bvf(i,j,k-1)))
# ifdef SASHA
            FC(i,k-1)=Ritop-lmd_Ric*Ribot
# else
            FC(i,k-1)=Ritop/(Ribot+eps)
# endif
          END DO
        END DO
!
! Linearly interpolate to find "hsbl" where Rib/Ric=1.
!
        DO i=Istr,Iend
          ksbl(i,j)=1
          hsbl(i,j)=z_w(i,j,1)
        END DO
# ifdef SASHA
        DO k=N(ng),2,-1
          DO i=Istr,Iend
            IF ((ksbl(i,j).eq.1).and.(FC(i,k-1).gt.0.0_r8)) THEN
              hsbl(i,j)=(z_w(i,j,k)*FC(i,k-1)-z_w(i,j,k-1)*FC(i,k))/    &
     &                  (FC(i,k-1)-FC(i,k))
              ksbl(i,j)=k
            END IF
          END DO
        END DO
# else
#  ifdef QUADRATIC
#ifdef KATE
!    Quadratic interpolation to find hsbl, using a scheme from Bill
!    Large et al. (personal communication)
        DO k=N(ng),2,-1
          DO i=Istr,Iend
!    CAE version.  Assumes that FC is calculated on w-points.
!    Solve FC = a*z^2 + b*z + c,
!    If F' is d(FC)/dz, either up or down a half grid-point,
!    a=(F'(k+.5)-F'(k-.5))/(z(k+1)-z(k-1)
!    b=F'(k+.5)-2*a*z(k+.5)
!    c=F(k)-a*z(k)*z(k)-bz(k)
!    or 0 = a*z^2 + b*z + (c-FC).
!    Let a_co=a, b_co=b, and c_co=c-Ric.
            IF (k .eq. N(ng)) THEN
              slope_up = 0.0_r8
            ELSE
              slope_up = (FC(i,k+1) - FC(i,k)) /                      &
     &                   (z_w(i,j,k+1) - z_w(i,j,k))
            END IF
            slope_down = (FC(i,k) - FC(i,k-1)) /                        &
     &                 (z_w(i,j,k) - z_w(i,j,k-1))
            a_co = (slope_up-slope_down)/(z_w(i,j,k+1)-z_w(i,j,k-1))
            b_co = slope_up - a_co*(z_w(i,j,k+1)+z_w(i,j,k))
            c_co = FC(i,k) - a_co*z_w(i,j,k)*z_w(i,j,k) -               &
     &               b_co*z_w(i,j,k)-lmd_Ric
            sqrt_arg = b_co*b_co - 4.0_r8*a_co*c_co
            IF ( ( (abs(b_co) .gt. eps) .and.                           &
     &           (abs(a_co/b_co) .le. eps) ) .or.                       &
     &           (sqrt_arg .le.  0.0_r8 ) ) THEN
!   If problem using quadratic, then use linear approx.
              hsbl(i,j)=((FC(i,k-1)-lmd_Ric)*z_w(i,j,k  )+              &
     &                    (lmd_Ric-FC(i,k  ))*z_w(i,j,k-1))/            &
     &                    (FC(i,k-1)-FC(i,k))
            ELSE
              hsbl(i,j) = (-b_co - sqrt(sqrt_arg)) / (2.0_r8*a_co)
              IF (hsbl(i,j) .gt. z_w(i,j,k)) THEN
                write(6,*) 'Bad calculation of hsbl. Exiting', i,j
                stop
              END IF
            END IF
            ksbl(i,j)=k
          END DO
        END DO
#endif
!    Quadratic interpolation to find hsbl, using a scheme from Bill
!    Large et al. (personal communication)
        DO k=N(ng),2,-1
          DO i=Istr,Iend
!              z_up = z_w(i,j,k)
!              IF (k .eq. N(ng)) THEN
!                slope_up = 0.0_r8
!              ELSE
!                slope_up = (FC(i,k+1) - FC(i,k)) /                     &
!     &                   (z_up - z_w(i,j,k+1))
!              END IF
!              a_co = (FC(i,k-1) - FC(i,k) +                            &
!     &              slope_up*(z_w(i,j,k-1) - z_up)) /                  &
!     &              (z_up - z_w(i,j,k-1))**2
!              b_co = slope_up + 2.0_r8*a_co*z_up
!              c_co = FC(i,k) + z_up*                                   &
!     &               (a_co*z_up + slope_up) - lmd_Ric
!              sqrt_arg = b_co**2 - 4.0_r8*a_co*c_co
!              if ( ( (abs(b_co) .gt. eps) .and.                        &
!     &           (abs(a_co/b_co) .le. eps) ) .or.                      &
!     &           (sqrt_arg .le.  0.0_r8 ) ) then
!                hsbl(i,j) = - z_up + (z_up - z_w(i,j,k-1))*            &
!     &                 (lmd_Ric - FC(i,k)) / (FC(i,k-1) - Fc(i,k))
!              else
!                hsbl(i,j) = (-b_co + sqrt(sqrt_arg)) / (2.0_r8*a_co)
!              endif
!              ksbl(i,j)=k
!
!             CAE version.  Assumes that FC is calculated on w-points
!             It's possible my version is the same as that above.  Not sure.
!             Solve FC = a*z^2 + b*z + c,
!             If F' is d(FC)/dz, either up or down a half grid-point,
!             a=(F'(k+.5)-F'(k-.5))/(z(k+1)-z(k-1)
!             b=F'(k+.5)-2*a*z(k+.5)
!             c=F(k)-a*z(k)*z(k)-bz(k)
!             or 0 = a*z^2 + b*z + (c-FC).
!             Let a_co=a, b_co=b, and c_co=c-Ric.
            IF (k .eq. N(ng)) THEN
              slope_up = 0.0_r8
            ELSE
              slope_up = (FC(i,k+1) - FC(i,k)) /                        &
     &                   (z_w(i,j,k+1) - z_w(i,j,k))
            END IF
            slope_down = (FC(i,k) - FC(i,k-1)) /                        &
     &                 (z_w(i,j,k) - z_w(i,j,k-1))
            a_co = (slope_up-slope_down)/(z_w(i,j,k+1)-z_w(i,j,k-1))
            b_co = slope_up - a_co*(z_w(i,j,k+1)+z_w(i,j,k))
            c_co = FC(i,k) - a_co*z_w(i,j,k)*z_w(i,j,k) -               &
     &               b_co*z_w(i,j,k)-lmd_Ric
            sqrt_arg = b_co*b_co - 4.0_r8*a_co*c_co
            if ( ( (abs(b_co) .gt. eps) .and.                           &
     &           (abs(a_co/b_co) .le. eps) ) .or.                       &
     &           (sqrt_arg .le.  0.0_r8 ) ) then
!               CAE: if problem using quadratic, then use linear approx.
              hsbl(i,j)=((FC(i,k-1)-lmd_Ric)*z_w(i,j,k  )+              &
     &                    (lmd_Ric-FC(i,k  ))*z_w(i,j,k-1))/            &
     &                    (FC(i,k-1)-FC(i,k))
            else
              hsbl(i,j) = (-b_co - sqrt(sqrt_arg)) / (2.0_r8*a_co)
              if ( (hsbl(i,j) .gt. z_w(i,j,k)) .or.                     &
     &               (hsbl(i,j) .gt. z_w(i,j,k)) ) then
                write(6,*) 'Bad calculation of hsbl Q. Exiting', i,j
! or bad test?
                stop
              endif
            endif
            ksbl(i,j)=k
          END DO
        END DO
#  else
        DO k=N(ng),2,-1
          DO i=Istr,Iend
            IF ((ksbl(i,j).eq.1).and.((FC(i,k  ).lt.lmd_Ric).and.       &
     &                                (FC(i,k-1).ge.lmd_Ric))) THEN
              hsbl(i,j)=((FC(i,k-1)-lmd_Ric)*z_w(i,j,k  )+              &
     &                  (lmd_Ric-FC(i,k  ))*z_w(i,j,k-1))/              &
     &                  (FC(i,k-1)-FC(i,k))
              ksbl(i,j)=k
            END IF
          END DO
        END DO
#  endif
# endif
      END DO
!
!  Compute total buoyancy flux at surface boundary layer depth,
!  "Bfsfc".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zgrid(i,j)=z_w(i,j,N(ng))-hsbl(i,j)
# ifdef MASKING
          zgrid(i,j)=zgrid(i,j)*rmask(i,j)
# endif
        END DO
      END DO
      CALL lmd_swfrac_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      -1.0_r8, zgrid, swdk)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Bfsfc(i,j)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
# ifdef MASKING
          Bfsfc(i,j)=Bfsfc(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!  Under neutral and stable conditions, the depth of the surface
!  boundary layer is required to be less than Ekman and Monin-Obukov
!  depths.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF ((Ustar(i,j).gt.0.0_r8).and.(Bfsfc(i,j).gt.0.0_r8)) THEN
            hekman=lmd_cekman*Ustar(i,j)/MAX(ABS(f(i,j)),eps)
            hmonob=lmd_cmonob*Ustar(i,j)*Ustar(i,j)*Ustar(i,j)/         &
     &             MAX(vonKar*Bfsfc(i,j),eps)
            hsbl(i,j)=(z_w(i,j,N(ng))-                                  &
     &                 MIN(hekman,hmonob,z_w(i,j,N(ng))-hsbl(i,j)))
          END IF
          hsbl(i,j)=MIN(hsbl(i,j),z_w(i,j,N(ng)))
          hsbl(i,j)=MAX(hsbl(i,j),z_w(i,j,0))
# ifdef MASKING
          hsbl(i,j)=hsbl(i,j)*rmask(i,j)
# endif
        END DO
      END DO
# ifdef LMD_SHAPIRO
!
!  Apply gradient or periodic boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  hsbl)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    hsbl)
#  endif
!
!  Apply Shapiro filter
!
      CALL shapiro2d_tile (ng, tile, iNLM,                              &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#  ifdef MASKING
     &                     rmask,                                       &
#  endif
     &                     hsbl)
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          hsbl(i,j)=MIN(hsbl(i,j),z_w(i,j,N(ng)))
          hsbl(i,j)=MAX(hsbl(i,j),z_w(i,j,0))
#  ifdef MASKING
          hsbl(i,j)=hsbl(i,j)*rmask(i,j)
#  endif
        END DO
      END DO
# endif
!
!  Apply gradient or periodic boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  hsbl)
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    hsbl)
# endif
!
!  Find new boundary layer index "ksbl".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ksbl(i,j)=1
          DO k=N(ng),2,-1
            IF ((ksbl(i,j).eq.1).and.(z_w(i,j,k-1).lt.hsbl(i,j))) THEN
              ksbl(i,j)=k
            END IF
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux at final surface boundary layer depth.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zgrid(i,j)=z_w(i,j,N(ng))-hsbl(i,j)
# ifdef MASKING
          zgrid(i,j)=zgrid(i,j)*rmask(i,j)
# endif
        END DO
      END DO
      CALL lmd_swfrac_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      -1.0_r8, zgrid, swdk)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Bfsfc(i,j)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
# ifdef MASKING
          Bfsfc(i,j)=Bfsfc(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!=======================================================================
!  Compute vertical mixing coefficients within the planetary boundary
!  layer.
!=======================================================================
!
!  Compute tubulent velocity scales (wm,ws) at "hsbl".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          sl_dpth(i,j)=lmd_epsilon*(z_w(i,j,N(ng))-hsbl(i,j))
          IF (Bfsfc(i,j).gt.0.0_r8) THEN
            cff=1.0_r8
          ELSE
            cff=lmd_epsilon
          END IF
          sigma=cff*(z_w(i,j,N(ng))-hsbl(i,j))
          Ustar3=Ustar(i,j)*Ustar(i,j)*Ustar(i,j)
          zetahat=vonKar*sigma*Bfsfc(i,j)
          zetapar=zetahat/(Ustar3+small)
          IF (zetahat.ge.0.0_r8) THEN                           ! stable
            wm(i,j)=vonKar*Ustar(i,j)/(1.0_r8+5.0_r8*zetapar)
            ws(i,j)=wm(i,j)
          ELSE                                                ! unstable
            IF (zetapar.gt.lmd_zetam) THEN
              wm(i,j)=vonKar*Ustar(i,j)*                                &
     &                (1.0_r8-16.0_r8*zetapar)**0.25_r8
            ELSE
              wm(i,j)=vonKar*(lmd_am*Ustar3-lmd_cm*zetahat)**r3
            END IF
            IF (zetapar.gt.lmd_zetas) THEN
              ws(i,j)=vonKar*Ustar(i,j)*                                &
     &                (1.0_r8-16.0_r8*zetapar)**0.5_r8
            ELSE
              ws(i,j)=vonKar*(lmd_as*Ustar3-lmd_cs*zetahat)**r3
            END IF
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute nondimensional shape function Gx(sigma) in terms of the
!  interior diffusivities at sigma=1 (Gm1, Gs1, Gt1) and its vertical
!  derivative evaluated "hsbl" via interpolation.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          f1(i,j)=5.0_r8*MAX(0.0_r8,Bfsfc(i,j))*vonKar/                 &
     &            (Ustar(i,j)*Ustar(i,j)*Ustar(i,j)*Ustar(i,j)+eps)
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zbl=z_w(i,j,N(ng))-hsbl(i,j)
          IF (hsbl(i,j).gt.z_w(i,j,1)) THEN
            k=ksbl(i,j)
            cff=1.0_r8/(z_w(i,j,k)-z_w(i,j,k-1))
            cff_dn=cff*(hsbl(i,j)-z_w(i,j,k-1))
            cff_up=cff*(z_w(i,j,k)-hsbl(i,j))
!
!  Compute nondimensional shape function for viscosity "Gm1" and its
!  vertical derivative "dGm1dS" evaluated at "hsbl".
!
            K_bl=cff_dn*Akv(i,j,k)+cff_up*Akv(i,j,k-1)
            dK_bl=cff*(Akv(i,j,k)-Akv(i,j,k-1))
            Gm1(i,j)=K_bl/(zbl*wm(i,j)+eps)
# ifdef MASKING
            Gm1(i,j)=Gm1(i,j)*rmask(i,j)
# endif
            dGm1dS(i,j)=MIN(0.0_r8,-dK_bl/(wm(i,j)+eps)-K_bl*f1(i,j))
!
!  Compute nondimensional shape function for diffusion of temperature
!  "Gt1" and its vertical derivative "dGt1dS" evaluated at "hsbl".
!
            K_bl=cff_dn*Akt(i,j,k,itemp)+cff_up*Akt(i,j,k-1,itemp)
            dK_bl=cff*(Akt(i,j,k,itemp)-Akt(i,j,k-1,itemp))
            Gt1(i,j)=K_bl/(zbl*ws(i,j)+eps)
# ifdef MASKING
            Gt1(i,j)=Gt1(i,j)*rmask(i,j)
# endif
            dGt1dS(i,j)=MIN(0.0_r8,-dK_bl/(ws(i,j)+eps)-K_bl*f1(i,j))
# ifdef SALINITY
!
!  Compute nondimensional shape function for diffusion of salinity
!  "Gs1" and its vertical derivative "dGs1dS" evaluated at "hsbl".
!
            K_bl=cff_dn*Akt(i,j,k,isalt)+cff_up*Akt(i,j,k-1,isalt)
            dK_bl=cff*(Akt(i,j,k,isalt)-Akt(i,j,k-1,isalt))
            Gs1(i,j)=K_bl/(zbl*ws(i,j)+eps)
#  ifdef MASKING
            Gs1(i,j)=Gs1(i,j)*rmask(i,j)
#  endif
            dGs1dS(i,j)=MIN(0.0_r8,-dK_bl/(ws(i,j)+eps)-K_bl*f1(i,j))
# endif
          ELSE
!
!  If the surface boundary layer extends to the bottom, assume that
!  the neutral boundary layer similarity theory holds at the bottom.
!
            ksbl(i,j)=0
!
!  Compute nondimensional bottom shape function for viscosity.
!
            Ustarb=SQRT(SQRT((0.5_r8*(bustr(i,j)+bustr(i+1,j)))**2+     &
     &                       (0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))**2))
# ifdef MASKING
            Ustarb=Ustarb*rmask(i,j)
# endif
            dK_bl=vonKar*Ustarb
            K_bl=dK_bl*(hsbl(i,j)-z_w(i,j,0))
            Gm1(i,j)=K_bl/(zbl*wm(i,j)+eps)
# ifdef MASKING
            Gm1(i,j)=Gm1(i,j)*rmask(i,j)
# endif
            dGm1dS(i,j)=MIN(0.0_r8,-dK_bl/(wm(i,j)+eps)-K_bl*f1(i,j))
!
!  Compute nondimensional bottom shape function for diffusion of
!  temperature.
!
            Gt1(i,j)=K_bl/(zbl*ws(i,j)+eps)
# ifdef MASKING
            Gt1(i,j)=Gt1(i,j)*rmask(i,j)
# endif
            dGt1dS(i,j)=MIN(0.0_r8,-dK_bl/(ws(i,j)+eps)-K_bl*f1(i,j))
!
!  Compute nondimensional bottom shape function for diffusion of
!  salinity.
!
# ifdef SALINITY
            Gs1(i,j)=Gt1(i,j)
            dGs1dS(i,j)=dGt1dS(i,j)
# endif
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute surface boundary layer mixing coefficients.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            zbl=z_w(i,j,N(ng))-hsbl(i,j)
            IF (k.gt.ksbl(i,j)) THEN
!
!  Compute turbulent velocity scales at vertical W-points.
!
              depth=z_w(i,j,N(ng))-z_w(i,j,k)
              IF (Bflux(i,j,k).lt.0.0_r8) THEN
                sigma=MIN(sl_dpth(i,j),depth)
              ELSE
                sigma=depth
              END IF
              Ustar3=Ustar(i,j)*Ustar(i,j)*Ustar(i,j)
              zetahat=vonKar*sigma*Bflux(i,j,k)
              zetapar=zetahat/(Ustar3+small)
              IF (zetahat.ge.0.0_r8) THEN                       ! stable
                wm(i,j)=vonKar*Ustar(i,j)/(1.0_r8+5.0_r8*zetapar)
                ws(i,j)=wm(i,j)
              ELSE                                            ! unstable
                IF (zetapar.gt.lmd_zetam) THEN
                  wm(i,j)=vonKar*Ustar(i,j)*                            &
     &                    (1.0_r8-16.0_r8*zetapar)**0.25_r8
                ELSE
                  wm(i,j)=vonKar*(lmd_am*Ustar3-lmd_cm*zetahat)**r3
                END IF
                IF (zetapar.gt.lmd_zetas) THEN
                  ws(i,j)=vonKar*Ustar(i,j)*                            &
     &                    (1.0_r8-16.0_r8*zetapar)**0.5_r8
                ELSE
                  ws(i,j)=vonKar*(lmd_as*Ustar3-lmd_cs*zetahat)**r3
                END IF
              END IF
!
!  Set polynomial coefficients for shape function.
!
              sigma=depth/(zbl+eps)
# ifdef MASKING
              sigma=sigma*rmask(i,j)
# endif
              a1=sigma-2.0_r8
              a2=3.0_r8-2.0_r8*sigma
              a3=sigma-1.0_r8
!
!  Compute nondimesional shape functions.
!
              Gm=a1+a2*Gm1(i,j)+a3*dGm1dS(i,j)
              Gt=a1+a2*Gt1(i,j)+a3*dGt1dS(i,j)
# ifdef SALINITY
              Gs=a1+a2*Gs1(i,j)+a3*dGs1dS(i,j)
# endif
!
!  Compute boundary layer mixing coefficients, combine them
!  with interior mixing coefficients.
!
# ifdef LMD_BOUND
              Akv(i,j,k)=MIN(lmd_nu0c,                                  &
     &                       depth*wm(i,j)*(1.0_r8+sigma*Gm))
              Akt(i,j,k,itemp)=MIN(lmd_nu0c,                            &
     &                             depth*ws(i,j)*(1.0_r8+sigma*Gt))
#  ifdef SALINITY
              AKt(i,j,k,isalt)=MIN(lmd_nu0c,                            &
     &                             depth*ws(i,j)*(1.0_r8+sigma*Gs))
#  endif
# else
              Akv(i,j,k)=depth*wm(i,j)*(1.0_r8+sigma*Gm)
              Akt(i,j,k,itemp)=depth*ws(i,j)*(1.0_r8+sigma*Gt)
#  ifdef SALINITY
              AKt(i,j,k,isalt)=depth*ws(i,j)*(1.0_r8+sigma*Gs)
#  endif
# endif
# ifdef LMD_NONLOCAL
!
!  Compute boundary layer nonlocal transport (m/s2).
!
              cff=lmd_Cg*(1.0_r8-(0.5_r8+SIGN(0.5_r8,Bflux(i,j,k))))/   &
     &            (zbl*ws(i,j)+eps)
              ghats(i,j,k,itemp)=cff*ghats(i,j,k,itemp)
#  ifdef SALINITY
              ghats(i,j,k,isalt)=cff*ghats(i,j,k,isalt)
#  endif
# endif
!
!  Set vertical mixing coefficients to fit neutral log layer
!  similarity theory.
!
            ELSE
# ifdef LMD_NONLOCAL
              ghats(i,j,k,itemp)=0.0_r8
#  ifdef SALINITY
              ghats(i,j,k,isalt)=0.0_r8
#  endif
# endif
            END IF
          END DO
        END DO
      END DO
#ifdef AKT_3D
      DO itrc=1,2
        DO j=Jstr-2,Jend+2
          DO k=1,N(ng)
            DO i=Istr-2,Iend+2
              Akt(i,j,k,itrc)=tclm(i,j,k,iAKt3)
            END DO
          END DO
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE lmd_skpp_tile
#endif
      END MODULE lmd_skpp_mod
