#include "cppdefs.h"
#define SASHA
      MODULE lmd_bkpp_mod
#if defined NONLINEAR && defined LMD_BKPP && defined SOLVE3D
!
!svn $Id: lmd_bkpp.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Scott M. Durski   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine determines the depth of  bottom  oceanic boundary      !
!  layer,  hbbl,  as the deepest depth  where the bulk Richardson      !
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
!  This routine was adapted from Bill Large 1995 code.                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: lmd_bkpp

      CONTAINS
!
!***********************************************************************
      SUBROUTINE lmd_bkpp (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL lmd_bkpp_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng),                                     &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    GRID(ng) % f,                                 &
     &                    GRID(ng) % h,                                 &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_r,                               &
     &                    GRID(ng) % z_w,                               &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v,                                &
     &                    OCEAN(ng) % pden,                             &
     &                    FORCES(ng) % srflx,                           &
     &                    FORCES(ng) % btflx,                           &
     &                    FORCES(ng) % bustr,                           &
     &                    FORCES(ng) % bvstr,                           &
     &                    MIXING(ng) % alpha,                           &
# ifdef SALINITY
     &                    MIXING(ng) % beta,                            &
# endif
     &                    MIXING(ng) % bvf,                             &
     &                    MIXING(ng) % ksbl,                            &
     &                    MIXING(ng) % Akt,                             &
     &                    MIXING(ng) % Akv,                             &
     &                    MIXING(ng) % kbbl,                            &
     &                    MIXING(ng) % hbbl)
      RETURN
      END SUBROUTINE lmd_bkpp
!
!***********************************************************************
      SUBROUTINE lmd_bkpp_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp,                                   &
# ifdef MASKING
     &                          rmask,                                  &
# endif
     &                          f, h, Hz, z_r, z_w,                     &
     &                          u, v, pden,                             &
     &                          srflx, btflx,  bustr, bvstr,            &
     &                          alpha,                                  &
# ifdef SALINITY
     &                          beta,                                   &
# endif
     &                          bvf,                                    &
     &                          ksbl, Akt, Akv, kbbl, hbbl)
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
      integer, Intent(in) :: nstp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: f(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: pden(LBi:,LBj:,:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: btflx(LBi:,LBj:,:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: alpha(LBi:,LBj:)
#  ifdef SALINITY
      real(r8), intent(in) :: beta(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)
      integer,  intent(in) :: ksbl(LBi:,LBj:)
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: hbbl(LBi:,LBj:)
      integer,  intent(out) :: kbbl(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: pden(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: alpha(LBi:UBi,LBj:UBj)
#  ifdef SALINITY
      real(r8), intent(in) :: beta(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
      integer,  intent(in) :: ksbl(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: hbbl(LBi:UBi,LBj:UBj)
      integer,  intent(out) :: kbbl(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8
      real(r8), parameter :: small = 1.0E-20_r8

      real(r8) :: Gm, Gt, Gs, K_bl, Ribot, Ritop, Rk
      real(r8) :: Uk, Ustar3, Vk, Vtc
      real(r8) :: a1, a2, a3, cff, cff1, cff2, cff_up, cff_dn
      real(r8) :: depth, dK_bl, hekman, sigma, zbl
      real(r8) :: zetahat, zetapar

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
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Bfbot
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Gm1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Gt1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Gs1
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: Ustar
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: bl_dpth
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: dGm1dS
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: dGt1dS
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: dGs1dS
      real(r8), dimension (IminS:ImaxS,JminS:JmaxS) :: f1
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
!  Get approximation of bottom layer depth using "lmd_eps" and
!  boundary layer depth from previous time step.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          bl_dpth(i,j)=lmd_epsilon*(hbbl(i,j)-z_w(i,j,0))
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute turbulent friction velocity (m/s) "Ustar" from bottom
!  stress at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Ustar(i,j)=SQRT(SQRT((0.5_r8*(bustr(i,j)+bustr(i+1,j)))**2+   &
     &                         (0.5_r8*(bvstr(i,j)+bvstr(i,j+1)))**2))
# ifdef MASKING
          Ustar(i,j)=Ustar(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute bottom turbulent buoyancy forcing "Bo" (m2/s3). Compute
!  surface radiative buoyancy forcing "Bosol" (m2/s3) that can be
!  used in shallow areas when it can penetrate all the way to the
!  bottom.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
# ifdef SALINITY
          Bo(i,j)=g*(alpha(i,j)*btflx(i,j,itemp)-                       &
     &               beta (i,j)*btflx(i,j,isalt))
# else
          Bo(i,j)=g*alpha(i,j)*btflx(i,j,itemp)
# endif
          Bosol(i,j)=g*alpha(i,j)*srflx(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux (m2/s3) at W-points.
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
          END DO
        END DO
      END DO
!
!=======================================================================
!  Compute bulk Richardson number "Rib" and then find depth of the
!  oceanic bottom boundary layer "hbbl", such that Rib(hbbl)=Ric.
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
!  bottom boundary layer "hbbl".
!
!                  [Br - B(d)] * d
!     Rib(d) = ----------------------- ;     Rib(hbbl)=Ric     (1)
!              |Vr - V(d)|^2 + Vt(d)^2
!
!  where "Br" and "Vr" are the bottom reference buoyancy and velocity
!  while "B(d)" and "V(d)" are the buoyancy and velocity at depth "d".
!
!  In the code below, the criterion "Rib(hbbl)=Ric" is reformulated
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
!  Compute potential density and velocity components bottom reference
!  values.
!
        cff1=1.0_r8/3.0_r8
        cff2=1.0_r8/6.0_r8
        DO i=Istr,Iend
          Rref(i)=pden(i,j,1)-                                          &
     &            Hz(i,j,1)*(cff1*dR(i,0)+cff2*dR(i,1))
          Uref(i)=0.5_r8*(u(i,j,1,nstp)+u(i+1,j,1,nstp))-               &
     &            Hz(i,j,1)*(cff1*dU(i,0)+cff2*dU(i,1))
          Vref(i)=0.5_r8*(v(i,j,1,nstp)+v(i,j+1,1,nstp))-               &
     &            Hz(i,j,1)*(cff1*dV(i,0)+cff2*dV(i,1))
        END DO
!
!  Compute turbulent velocity scales for momentum (wm) and tracers (ws).
!  Then, compute critical function (FC) for bulk Richardson number.
!
        DO i=Istr,Iend
          FC(i,0)=0.0_r8
          DO k=1,N(ng)
            depth=z_w(i,j,k)-z_w(i,j,0)
            IF (Bflux(i,j,k).lt.0.0_r8) THEN
              sigma=MIN(bl_dpth(i,j),depth)
            ELSE
              sigma=depth
            END IF
            Ustar3=Ustar(i,j)*Ustar(i,j)*Ustar(i,j)
            zetahat=vonKar*sigma*Bflux(i,j,k)
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
            Rk=pden(i,j,k)+                                             &
     &         Hz(i,j,k)*(cff1*dR(i,k)+cff2*dR(i,k-1))
            Uk=0.5_r8*(u(i,j,k,nstp)+u(i+1,j,k,nstp))+                  &
     &         Hz(i,j,k)*(cff1*dU(i,k)+cff2*dU(i,k-1))
            Vk=0.5_r8*(v(i,j,k,nstp)+v(i,j+1,k,nstp))+                  &
     &         Hz(i,j,k)*(cff1*dV(i,k)+cff2*dV(i,k-1))
!
            Ritop=-gorho0*(Rk-Rref(i))*depth
            Ribot=(Uk-Uref(i))**2+(Vk-Vref(i))**2+                      &
     &            Vtc*depth*ws(i,j)*SQRT(ABS(bvf(i,j,k)))
# ifdef SASHA
            FC(i,k)=Ritop-lmd_Ric*Ribot
# else
            FC(i,k)=Ritop/(Ribot+eps)
# endif
          END DO
        END DO
!
! Linearly interpolate to find "hbbl" where Rib/Ric=1.
!
        DO i=Istr,Iend
          kbbl(i,j)=N(ng)
          hbbl(i,j)=z_w(i,j,N(ng))
        END DO
# ifdef SASHA
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            IF ((kbbl(i,j).eq.N(ng)).and.(FC(i,k).gt.0.0_r8)) THEN
              hbbl(i,j)=(z_w(i,j,k)*FC(i,k-1)-z_w(i,j,k-1)*FC(i,k))/    &
     &                  (FC(i,k-1)-FC(i,k))
              kbbl(i,j)=k
            END IF
          END DO
        END DO
# else
        DO k=1,N(ng)
          DO i=Istr,Iend
            IF ((kbbl(i,j).eq.N(ng)).and.((FC(i,k-1).lt.lmd_Ric).and.   &
     &                                    (FC(i,k  ).ge.lmd_Ric))) THEN
              hbbl(i,j)=((lmd_Ric-FC(i,k-1))*z_w(i,j,k  )+              &
     &                   (FC(i,k  )-lmd_Ric)*z_w(i,j,k-1))/             &
     &                  (FC(i,k)-FC(i,k-1))
              kbbl(i,j)=k
            END IF
          END DO
        END DO
# endif
      END DO
!
!  Compute total buoyancy flux at bottom boundary layer depth,
!  "Bfbot".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zgrid(i,j)=z_w(i,j,N(ng))-hbbl(i,j)
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
          Bfbot(i,j)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
# ifdef MASKING
          Bfbot(i,j)=Bfbot(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!  Under neutral and stable conditions, the depth of the bottom
!  boundary layer is required to be less than Ekman and Monin-Obukov
!  depths.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Ustar(i,j).ge.0.0_r8) THEN
            hekman=lmd_cekman*Ustar(i,j)/MAX(ABS(f(i,j)),eps)-h(i,j)
            hbbl(i,j)= MIN(hekman,hbbl(i,j))
          END IF
          hbbl(i,j)=MIN(hbbl(i,j),z_w(i,j,N(ng)))
          hbbl(i,j)=MAX(hbbl(i,j),z_w(i,j,0))
# ifdef MASKING
          hbbl(i,j)=hbbl(i,j)*rmask(i,j)
# endif
        END DO
      END DO
# ifdef LMD_SHAPIRO
!
!  Apply gradient or periodic boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  hbbl)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    hbbl)
#  endif
!
!  Shapiro filter on boundary layer thickness
!
      CALL shapiro2d_tile (ng, tile, iNLM,                              &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#  ifdef MASKING
     &                     rmask,                                       &
#  endif
     &                     hbbl)
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          hbbl(i,j)=MIN(hbbl(i,j),z_w(i,j,N(ng)))
          hbbl(i,j)=MAX(hbbl(i,j),z_w(i,j,0))
#  ifdef MASKING
          hbbl(i,j)=hbbl(i,j)*rmask(i,j)
#  endif
        END DO
      END DO
# endif
!
!  Apply gradient or periodic boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  hbbl)
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    hbbl)
# endif
!
!  Find new boundary layer index "kbbl".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          kbbl(i,j)=N(ng)
          DO k=1,N(ng)
            IF ((kbbl(i,j).eq.N(ng)).and.(z_w(i,j,k).gt.hbbl(i,j))) THEN
              kbbl(i,j)=k
            END IF
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute total buoyancy flux at final bottom boundary layer depth.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zgrid(i,j)=z_w(i,j,N(ng))-hbbl(i,j)
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
          Bfbot(i,j)=(Bo(i,j)+Bosol(i,j)*(1.0_r8-swdk(i,j)))
# ifdef MASKING
          Bfbot(i,j)=Bfbot(i,j)*rmask(i,j)
# endif
        END DO
      END DO
!
!=======================================================================
!  Compute vertical mixing coefficients within the planetary boundary
!  layer.
!=======================================================================
!
!  Compute tubulent velocity scales (wm,ws) at "hbbl".
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          bl_dpth(i,j)=lmd_epsilon*(hbbl(i,j)-z_w(i,j,0))
          IF (Bfbot(i,j).gt.0.0_r8) THEN
            cff=1.0_r8
          ELSE
            cff=lmd_epsilon
          END IF
          sigma=cff*(hbbl(i,j)-z_w(i,j,0))
          Ustar3=Ustar(i,j)*Ustar(i,j)*Ustar(i,j)
          zetahat=vonKar*sigma*Bfbot(i,j)
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
!  derivative evaluated "hbbl" via interpolation.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          f1(i,j)=5.0_r8*MAX(0.0_r8,Bfbot(i,j))*vonKar/                 &
     &            (Ustar(i,j)*Ustar(i,j)*Ustar(i,j)*Ustar(i,j)+eps)
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          zbl=hbbl(i,j)-z_w(i,j,0)
          k=kbbl(i,j)
          cff=1.0_r8/(z_w(i,j,k)-z_w(i,j,k-1))
          cff_dn=cff*(hbbl(i,j)-z_w(i,j,k-1))
          cff_up=cff*(z_w(i,j,k)-hbbl(i,j))
!
!  Compute nondimensional shape function for viscosity "Gm1" and its
!  vertical derivative "dGm1dS" evaluated at "hbbl".
!
          K_bl=cff_dn*Akv(i,j,k)+cff_up*Akv(i,j,k-1)
          dK_bl=-cff*(Akv(i,j,k)-Akv(i,j,k-1))
          Gm1(i,j)=K_bl/(zbl*wm(i,j)+eps)
# ifdef MASKING
          Gm1(i,j)=Gm1(i,j)*rmask(i,j)
# endif
          dGm1dS(i,j)=MIN(0.0_r8, K_bl*f1(i,j)-dK_bl/(wm(i,j)+eps))
!
!  Compute nondimensional shape function for diffusion of temperature
!  "Gt1" and its vertical derivative "dGt1dS" evaluated at "hbbl".
!
          K_bl=cff_dn*Akt(i,j,k,itemp)+cff_up*Akt(i,j,k-1,itemp)
          dK_bl=-cff*(Akt(i,j,k,itemp)-Akt(i,j,k-1,itemp))
          Gt1(i,j)=K_bl/(zbl*ws(i,j)+eps)
# ifdef MASKING
          Gt1(i,j)=Gt1(i,j)*rmask(i,j)
# endif
          dGt1dS(i,j)=MIN(0.0_r8, K_bl*f1(i,j)-dK_bl/(ws(i,j)+eps))
# ifdef SALINITY
!
!  Compute nondimensional shape function for diffusion of salinity
!  "Gs1" and its vertical derivative "dGs1dS" evaluated at "hbbl".
!
          K_bl=cff_dn*Akt(i,j,k,isalt)+cff_up*Akt(i,j,k-1,isalt)
          dK_bl=-cff*(Akt(i,j,k,isalt)-Akt(i,j,k-1,isalt))
          Gs1(i,j)=K_bl/(zbl*ws(i,j)+eps)
#  ifdef MASKING
          Gs1(i,j)=Gs1(i,j)*rmask(i,j)
#  endif
          dGs1dS(i,j)=MIN(0.0_r8, K_bl*f1(i,j)-dK_bl/(ws(i,j)+eps))
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute bottom boundary layer mixing coefficients.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            IF (z_w(i,j,k).lt.hbbl(i,j)) THEN
!
!  Compute turbulent velocity scales at vertical W-points.
!
              depth=z_w(i,j,k)-z_w(i,j,0)
              IF (Bflux(i,j,k).lt.0.0_r8) THEN
                sigma=MIN(bl_dpth(i,j),depth)
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
              sigma=depth/(hbbl(i,j)-z_w(i,j,0)+eps)
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
!  with interior mixing coefficients. Take the maximum estimate
!  of vertical mixing only where the surface and bottom boundary
!  layers overlap. (Do not let the interior overwrite the boundary
!  layer estimate).
!
              IF (k.gt.ksbl(i,j)) THEN
                Akv(i,j,k)=MAX(Akv(i,j,k),                              &
     &                         depth*wm(i,j)*(1.0_r8+sigma*Gm))
                Akt(i,j,k,itemp)=MAX(Akt(i,j,k,itemp),                  &
     &                               depth*ws(i,j)*(1.0_r8+sigma*Gt))
# ifdef SALINITY
                Akt(i,j,k,isalt)=MAX(Akt(i,j,k,isalt),                  &
     &                               depth*ws(i,j)*(1.0_r8+sigma*Gs))
# endif
              ELSE
                Akv(i,j,k)=depth*wm(i,j)*(1.0_r8+sigma*Gm)
                Akt(i,j,k,itemp)=depth*ws(i,j)*(1.0_r8+sigma*Gt)
# ifdef SALINITY
                Akt(i,j,k,isalt)=depth*ws(i,j)*(1.0_r8+sigma*Gs)
# endif
              END IF
            END IF
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE lmd_bkpp_tile
#endif
      END MODULE lmd_bkpp_mod
