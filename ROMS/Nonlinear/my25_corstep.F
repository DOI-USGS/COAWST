#include "cppdefs.h"
      MODULE my25_corstep_mod
#if defined NONLINEAR && defined MY25_MIXING && defined SOLVE3D
!
!svn $Id: my25_corstep.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine perfoms the corrector step for turbulent kinetic       !
!  energy and length scale prognostic variables, tke and gls.          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Mellor, G.L. and T. Yamada, 1982:  Development of turbulence        !
!    closure model for geophysical fluid problems, Rev. Geophys.       !
!    Space Phys., 20, 851-875.                                         !
!                                                                      !
!  Galperin, B., L.H. Kantha, S. Hassid, and A.Rosati, 1988:  A        !
!    quasi-equilibrium  turbulent  energy model for geophysical        !
!    flows, J. Atmos. Sci., 45, 55-62.                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: my25_corstep

      CONTAINS
!
!***********************************************************************
      SUBROUTINE my25_corstep (ng, tile)
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
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 20)
# endif
      CALL my25_corstep_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
# endif
     &                        GRID(ng) % Huon,                          &
     &                        GRID(ng) % Hvom,                          &
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % pm,                            &
     &                        GRID(ng) % pn,                            &
     &                        GRID(ng) % z_r,                           &
     &                        GRID(ng) % z_w,                           &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
     &                        OCEAN(ng) % W,                            &
     &                        FORCES(ng) % bustr,                       &
     &                        FORCES(ng) % bvstr,                       &
     &                        FORCES(ng) % sustr,                       &
     &                        FORCES(ng) % svstr,                       &
     &                        MIXING(ng) % bvf,                         &
     &                        MIXING(ng) % Akt,                         &
     &                        MIXING(ng) % Akv,                         &
     &                        MIXING(ng) % Akk,                         &
     &                        MIXING(ng) % Lscale,                      &
     &                        MIXING(ng) % gls,                         &
     &                        MIXING(ng) % tke)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 20)
# endif
      RETURN
      END SUBROUTINE my25_corstep
!
!***********************************************************************
      SUBROUTINE my25_corstep_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
# ifdef MASKING
     &                              umask, vmask,                       &
# endif
     &                              Huon, Hvom, Hz, pm, pn, z_r, z_w,   &
     &                              u, v, W,                            &
     &                              bustr, bvstr, sustr, svstr,         &
     &                              bvf, Akt, Akv,                      &
     &                              Akk, Lscale, gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
      USE tkebc_mod, ONLY : tkebc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)

      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akk(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Lscale(LBi:,LBj:,0:)
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: Akk(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: Lscale(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: gls(LBi:UBi,LBj:UBj,0:N(ng),3)
      real(r8), intent(inout) :: tke(LBi:UBi,LBj:UBj,0:N(ng),3)
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k

      real(r8), parameter :: Gadv = 1.0_r8/3.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8) :: Gh, Ls_unlmt, Ls_lmt, Qprod, Qdiss, Sh, Sm, Wscale
      real(r8) :: cff, cff1, cff2, cff3, ql, strat2

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BCK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BCP
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCP
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dU
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dV

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shear2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: buoy2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curvK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curvP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradP

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute vertical velocity shear at W-points.
!-----------------------------------------------------------------------
!
# ifdef RI_SPLINES
      DO j=Jstrm1,Jendp1
        DO i=Istrm1,Iendp1
          CF(i,0)=0.0_r8
          dU(i,0)=0.0_r8
          dV(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=Istrm1,Iendp1
            cff=1.0_r8/(2.0_r8*Hz(i,j,k+1)+                             &
     &                  Hz(i,j,k)*(2.0_r8-CF(i,k-1)))
            CF(i,k)=cff*Hz(i,j,k+1)
            dU(i,k)=cff*(3.0_r8*(u(i  ,j,k+1,nstp)-u(i,  j,k,nstp)+     &
     &                           u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))-    &
     &                   Hz(i,j,k)*dU(i,k-1))
            dV(i,k)=cff*(3.0_r8*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+     &
     &                           v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp))-    &
     &                   Hz(i,j,k)*dV(i,k-1))
          END DO
        END DO
        DO i=Istrm1,Iendp1
          dU(i,N(ng))=0.0_r8
          dV(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=Istrm1,Iendp1
            dU(i,k)=dU(i,k)-CF(i,k)*dU(i,k+1)
            dV(i,k)=dV(i,k)-CF(i,k)*dV(i,k+1)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istrm1,Iendp1
            shear2(i,j,k)=dU(i,k)*dU(i,k)+dV(i,k)*dV(i,k)
          END DO
        END DO
      END DO
# else
      DO k=1,N(ng)-1
        DO j=Jstrm1,Jendp1
          DO i=Istrm1,Iendp1
            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            shear2(i,j,k)=(cff*(u(i  ,j,k+1,nstp)-u(i  ,j,k,nstp)+      &
     &                          u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp)))**2+ &
     &                    (cff*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+      &
     &                          v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))**2
          END DO
        END DO
      END DO
# endif
!
! Load Brunt-Vaisala frequency.
!
      DO k=1,N(ng)-1
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            buoy2(i,j,k)=bvf(i,j,k)
          END DO
        END DO
      END DO
# ifdef N2S2_HORAVG
!
!-----------------------------------------------------------------------
!  Smooth horizontally buoyancy and shear.  Use buoy2(:,:,0) and
!  shear2(:,:,0) as scratch utility array.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
            shear2(Istr-1,j,k)=shear2(Istr,j,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=MAX(1,Jstr-1),MIN(Jend+1,Mm(ng))
            shear2(Iend+1,j,k)=shear2(Iend,j,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            shear2(i,Jstr-1,k)=shear2(i,Jstr,k)
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=MAX(1,Istr-1),MIN(Iend+1,Lm(ng))
            shear2(i,Jend+1,k)=shear2(i,Jend,k)
          END DO
        END IF
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          shear2(Istr-1,Jstr-1,k)=shear2(Istr,Jstr,k)
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          shear2(Istr-1,Jend+1,k)=shear2(Istr,Jend,k)
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          shear2(Iend+1,Jstr-1,k)=shear2(Iend,Jstr,k)
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          shear2(Iend+1,Jend+1,k)=shear2(Iend,Jend,k)
        END IF
!
!  Average horizontally.
!
        DO j=Jstr-1,Jend
          DO i=Istr-1,Iend
            buoy2(i,j,0)=0.25_r8*(buoy2(i,j  ,k)+buoy2(i+1,j  ,k)+      &
     &                            buoy2(i,j+1,k)+buoy2(i+1,j+1,k))
            shear2(i,j,0)=0.25_r8*(shear2(i,j  ,k)+shear2(i+1,j  ,k)+   &
     &                             shear2(i,j+1,k)+shear2(i+1,j+1,k))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            buoy2(i,j,k)=0.25_r8*(buoy2(i,j  ,0)+buoy2(i-1,j  ,0)+      &
     &                            buoy2(i,j-1,0)+buoy2(i-1,j-1,0))
            shear2(i,j,k)=0.25_r8*(shear2(i,j  ,0)+shear2(i-1,j  ,0)+   &
     &                             shear2(i,j-1,0)+shear2(i-1,j-1,0))
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Time-step advective terms.
!-----------------------------------------------------------------------
!
!  At entry, it is assumed that the turbulent kinetic energy fields
!  "tke" and "gls", at time level "nnew", are set to its values at
!  time level "nstp" multiplied by the grid box thicknesses Hz
!  (from old time step and at W-points).
!
      DO k=1,N(ng)-1
# ifdef K_C2ADVECTION
!
!  Second-order, centered differences advection.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.25_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            FXK(i,j)=cff*(tke(i,j,k,3)+tke(i-1,j,k,3))
            FXP(i,j)=cff*(gls(i,j,k,3)+gls(i-1,j,k,3))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.25_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FEK(i,j)=cff*(tke(i,j,k,3)+tke(i,j-1,k,3))
            FEP(i,j)=cff*(gls(i,j,k,3)+gls(i,j-1,k,3))
          END DO
        END DO
# else
        DO j=Jstr,Jend
          DO i=Istrm1,Iendp2
            gradK(i,j)=(tke(i,j,k,3)-tke(i-1,j,k,3))
#  ifdef MASKING
            gradK(i,j)=gradK(i,j)*umask(i,j)
#  endif
            gradP(i,j)=(gls(i,j,k,3)-gls(i-1,j,k,3))
#  ifdef MASKING
            gradP(i,j)=gradP(i,j)*umask(i,j)
#  endif
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              gradK(Istr-1,j)=gradK(Istr,j)
              gradP(Istr-1,j)=gradP(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              gradK(Iend+2,j)=gradK(Iend+1,j)
              gradP(Iend+2,j)=gradP(Iend+1,j)
            END DO
          END IF
        END IF
#  ifdef K_C4ADVECTION
!
!  Fourth-order, centered differences advection.
!
        cff1=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            FXK(i,j)=cff*0.5_r8*(tke(i-1,j,k,3)+tke(i,j,k,3)-           &
     &                           cff1*(gradK(i+1,j)-gradK(i-1,j)))
            FXP(i,j)=cff*0.5_r8*(gls(i-1,j,k,3)+gls(i,j,k,3)-           &
     &                           cff1*(gradP(i+1,j)-gradP(i-1,j)))
          END DO
        END DO
#  else
!
!  Third-order, upstream bias advection with velocity dependent
!  hyperdiffusion.
!
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
            curvK(i,j)=gradK(i+1,j)-gradK(i,j)
            curvP(i,j)=gradP(i+1,j)-gradP(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            IF (cff.gt.0.0_r8) THEN
              cff1=curvK(i-1,j)
              cff2=curvP(i-1,j)
            ELSE
              cff1=curvK(i,j)
              cff2=curvP(i,j)
            END IF
            FXK(i,j)=cff*0.5_r8*(tke(i-1,j,k,3)+tke(i,j,k,3)-           &
     &                           Gadv*cff1)
            FXP(i,j)=cff*0.5_r8*(gls(i-1,j,k,3)+gls(i,j,k,3)-           &
     &                           Gadv*cff2)
          END DO
        END DO
#  endif
        DO j=Jstrm1,Jendp2
          DO i=Istr,Iend
            gradK(i,j)=(tke(i,j,k,3)-tke(i,j-1,k,3))
#  ifdef MASKING
            gradK(i,j)=gradK(i,j)*vmask(i,j)
#  endif
            gradP(i,j)=(gls(i,j,k,3)-gls(i,j-1,k,3))
#  ifdef MASKING
            gradP(i,j)=gradP(i,j)*vmask(i,j)
#  endif
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              gradK(i,Jstr-1)=gradK(i,Jstr)
              gradP(i,Jstr-1)=gradP(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              gradK(i,Jend+2)=gradK(i,Jend+1)
              gradP(i,Jend+2)=gradP(i,Jend+1)
            END DO
          END IF
        END IF
#  ifdef K_C4ADVECTION
        cff1=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.5_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FEK(i,j)=cff*0.5_r8*(tke(i,j-1,k,3)+tke(i,j,k,3)-           &
     &                           cff1*(gradK(i,j+1)-gradK(i,j-1)))
            FEP(i,j)=cff*0.5_r8*(gls(i,j-1,k,3)+gls(i,j,k,3)-           &
     &                           cff1*(gradP(i,j+1)-gradP(i,j-1)))
          END DO
        END DO
#  else
        DO j=Jstr-1,Jend+1
          DO i=Istr,Iend
            curvK(i,j)=gradK(i,j+1)-gradK(i,j)
            curvP(i,j)=gradP(i,j+1)-gradP(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.5_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            IF (cff.gt.0.0_r8) THEN
              cff1=curvK(i,j-1)
              cff2=curvP(i,j-1)
            ELSE
              cff1=curvK(i,j)
              cff2=curvP(i,j)
            END IF
            FEK(i,j)=cff*0.5_r8*(tke(i,j-1,k,3)+tke(i,j,k,3)-           &
     &                           Gadv*cff1)
            FEP(i,j)=cff*0.5_r8*(gls(i,j-1,k,3)+gls(i,j,k,3)-           &
     &                           Gadv*cff2)
          END DO
        END DO
#  endif
# endif
!
!  Time-step horizontal advection.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-                            &
     &                      cff*(FXK(i+1,j)-FXK(i,j)+                   &
     &                           FEK(i,j+1)-FEK(i,j))
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-                            &
     &                      cff*(FXP(i+1,j)-FXP(i,j)+                   &
     &                           FEP(i,j+1)-FEP(i,j))
          END DO
        END DO
      END DO
!
! Compute vertical advection.
!
      DO j=Jstr,Jend
# ifdef K_C2ADVECTION
        DO k=1,N(ng)
          DO i=Istr,Iend
            cff=0.25_r8*(W(i,j,k)+W(i,j,k-1))
            FCK(i,k)=cff*(tke(i,j,k,3)+tke(i,j,k-1,3))
            FCP(i,k)=cff*(gls(i,j,k,3)+gls(i,j,k-1,3))
          END DO
        END DO
# else
        cff1=7.0_r8/12.0_r8
        cff2=1.0_r8/12.0_r8
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            cff=0.5*(W(i,j,k)+W(i,j,k-1))
            FCK(i,k)=cff*(cff1*(tke(i,j,k-1,3)+                         &
     &                          tke(i,j,k  ,3))-                        &
     &                    cff2*(tke(i,j,k-2,3)+                         &
     &                          tke(i,j,k+1,3)))
            FCP(i,k)=cff*(cff1*(gls(i,j,k-1,3)+                         &
     &                          gls(i,j,k  ,3))-                        &
     &                    cff2*(gls(i,j,k-2,3)+                         &
     &                          gls(i,j,k+1,3)))
          END DO
        END DO
        cff1=1.0_r8/3.0_r8
        cff2=5.0_r8/6.0_r8
        cff3=1.0_r8/6.0_r8
         DO i=Istr,Iend
          cff=0.5_r8*(W(i,j,0)+W(i,j,1))
          FCK(i,1)=cff*(cff1*tke(i,j,0,3)+                              &
     &                  cff2*tke(i,j,1,3)-                              &
     &                  cff3*tke(i,j,2,3))
          FCP(i,1)=cff*(cff1*gls(i,j,0,3)+                              &
     &                  cff2*gls(i,j,1,3)-                              &
     &                  cff3*gls(i,j,2,3))
          cff=0.5_r8*(W(i,j,N(ng))+W(i,j,N(ng)-1))
          FCK(i,N(ng))=cff*(cff1*tke(i,j,N(ng)  ,3)+                    &
     &                      cff2*tke(i,j,N(ng)-1,3)-                    &
     &                      cff3*tke(i,j,N(ng)-2,3))
          FCP(i,N(ng))=cff*(cff1*gls(i,j,N(ng)  ,3)+                    &
     &                      cff2*gls(i,j,N(ng)-1,3)-                    &
     &                      cff3*gls(i,j,N(ng)-2,3))
        END DO
# endif
!
!  Time-step vertical advection term.
!
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-                            &
     &                      cff*(FCK(i,k+1)-FCK(i,k))
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-                            &
     &                      cff*(FCP(i,k+1)-FCP(i,k))
          END DO
        END DO
!
!----------------------------------------------------------------------
!  Compute vertical mixing, turbulent production and turbulent
!  dissipation terms.
!----------------------------------------------------------------------
!
!  Set term for vertical mixing of turbulent fields.
!
        cff=-0.5_r8*dt(ng)
        DO k=1,N(ng)
          DO i=Istr,Iend
            FCK(i,k)=cff*(Akk(i,j,k)+Akk(i,j,k-1))/Hz(i,j,k)
            CF(i,k)=0.0_r8
          END DO
        END DO
!
!  Compute production and dissipation terms.
!
        cff3=my_E2/(vonKar*vonKar)
        DO k=1,N(ng)-1
          DO i=Istr,Iend
!
!  Compute shear and bouyant production of turbulent energy (m3/s3)
!  at W-points (ignore small negative values of buoyancy).
!
            IF ((buoy2(i,j,k).gt.-5.0E-5_r8).and.                       &
     &          (buoy2(i,j,k).lt.0.0_r8)) THEN
              strat2=0.0_r8
            ELSE
              strat2=buoy2(i,j,k)
            END IF
            Qprod=shear2(i,j,k)*(Akv(i,j,k)-Akv_bak(ng))-               &
     &            strat2*(Akt(i,j,k,itemp)-Akt_bak(itemp,ng))
!
!  Recalculate old time-step unlimited length scale.
!
            Ls_unlmt=MAX(eps,                                           &
     &                   gls(i,j,k,nstp)/(MAX(tke(i,j,k,nstp),eps)))
!
!  Time-step production term.
!
            cff1=0.5_r8*(Hz(i,j,k)+Hz(i,j,k+1))
            tke(i,j,k,nnew)=tke(i,j,k,nnew)+                            &
     &                      dt(ng)*cff1*Qprod*2.0_r8
            gls(i,j,k,nnew)=gls(i,j,k,nnew)+                            &
     &                      dt(ng)*cff1*Qprod*my_E1*Ls_unlmt
!
!  Compute dissipation of turbulent energy (m3/s3).  Add in vertical
!  mixing term.
!
            Qdiss=dt(ng)*SQRT(tke(i,j,k,nstp))/(my_B1*Ls_unlmt)
            cff=Ls_unlmt*(1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,k))+           &
     &                    1.0_r8/(z_w(i,j,k)-z_w(i,j,0)))
            Wscale=1.0_r8+cff3*cff*cff
            BCK(i,k)=cff1*(1.0_r8+2.0_r8*Qdiss)-FCK(i,k)-FCK(i,k+1)
            BCP(i,k)=cff1*(1.0_r8+Wscale*Qdiss)-FCK(i,k)-FCK(i,k+1)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step dissipation and vertical diffusion terms implicitly.
!-----------------------------------------------------------------------
!
!  Set surface and bottom boundary conditions.
!
        DO i=Istr,Iend
          tke(i,j,N(ng),nnew)=my_B1p2o3*0.5_r8*                         &
     &                        SQRT((sustr(i,j)+sustr(i+1,j))**2+        &
     &                             (svstr(i,j)+svstr(i,j+1))**2)
          gls(i,j,N(ng),nnew)=0.0_r8
          tke(i,j,0,nnew)=my_B1p2o3*0.5_r8*                             &
     &                    SQRT((bustr(i,j)+bustr(i+1,j))**2+            &
     &                         (bvstr(i,j)+bvstr(i,j+1))**2)
          gls(i,j,0,nnew)=0.0_r8
        END DO
!
!  Solve tri-diagonal system for "tke".
!
        DO i=Istr,Iend
          cff=1.0_r8/BCK(i,N(ng)-1)
          CF(i,N(ng)-1)=cff*FCK(i,N(ng)-1)
          tke(i,j,N(ng)-1,nnew)=cff*(tke(i,j,N(ng)-1,nnew)-             &
     &                               FCK(i,N(ng))*tke(i,j,N(ng),nnew))
        END DO
        DO k=N(ng)-2,1,-1
          DO i=Istr,Iend
            cff=1.0_r8/(BCK(i,k)-CF(i,k+1)*FCK(i,k+1))
            CF(i,k)=cff*FCK(i,k)
            tke(i,j,k,nnew)=cff*(tke(i,j,k,nnew)-                       &
     &                           FCK(i,k+1)*tke(i,j,k+1,nnew))
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-CF(i,k)*tke(i,j,k-1,nnew)
          END DO
        END DO
!
!  Solve tri-diagonal system for "gls".
!
        DO i=Istr,Iend
          cff=1.0_r8/BCP(i,N(ng)-1)
          CF(i,N(ng)-1)=cff*FCK(i,N(ng)-1)
          gls(i,j,N(ng)-1,nnew)=cff*(gls(i,j,N(ng)-1,nnew)-             &
     &                               FCK(i,N(ng))*gls(i,j,N(ng),nnew))
        END DO
        DO k=N(ng)-2,1,-1
          DO i=Istr,Iend
            cff=1.0_r8/(BCP(i,k)-CF(i,k+1)*FCK(i,k+1))
            CF(i,k)=cff*FCK(i,k)
            gls(i,j,k,nnew)=cff*(gls(i,j,k,nnew)-                       &
     &                           FCK(i,k+1)*gls(i,j,k+1,nnew))
          END DO
        END DO
        DO k=1,N(ng)-1,+1
          DO i=Istr,Iend
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-CF(i,k)*gls(i,j,k-1,nnew)
          END DO
        END DO
!
!---------------------------------------------------------------------
!  Compute vertical mixing coefficients (m2/s).
!---------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=Istr,Iend
!
!  Compute turbulent length scale (m).  The length scale is only
!  limited in the K-related calculations and not in QL production,
!  dissipation, wall-proximity, etc.
!
            tke(i,j,k,nnew)=MAX(tke(i,j,k,nnew),my_qmin)
            gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),my_qmin)
            Ls_unlmt=gls(i,j,k,nnew)/tke(i,j,k,nnew)
            Ls_lmt=MIN(Ls_unlmt,                                        &
     &                 my_lmax*SQRT(tke(i,j,k,nnew)/                    &
     &                         (MAX(0.0_r8,buoy2(i,j,k))+eps)))
!
!  Compute Galperin et al. (1988) nondimensional stability function,
!  Gh.  Then, compute nondimensional stability functions for tracers
!  (Sh) and momentum (Sm).  The limit on length scale, sets the lower
!  limit on Gh.  Use Kantha and Clayton or Galperin et al. expression
!  for Sm.
!
            Gh=MIN(my_Gh0,-buoy2(i,j,k)*Ls_lmt*Ls_lmt/                  &
     &                    tke(i,j,k,nnew))
            cff=1.0_r8-my_Sh2*Gh
            Sh=my_Sh1/cff
# ifdef KANTHA_CLAYSON
            Sm=(my_B1pm1o3+Sh*Gh*my_Sm4)/(1.0_r8-my_Sm2*Gh)
# else
            Sm=(my_Sm3+Sh*Gh*my_Sm4)/(1.0_r8-my_Sm2*Gh)
# endif
!
!  Compute vertical mixing (m2/s) coefficients of momentum and
!  tracers.  Average ql over the two timesteps rather than using
!  the new Lscale and just averaging tke.
!
            ql=0.5_r8*(Ls_lmt*SQRT(tke(i,j,k,nnew))+                    &
     &                 Lscale(i,j,k)*SQRT(tke(i,j,k,nstp)))
            Akv(i,j,k)=Akv_bak(ng)+ql*Sm
            DO itrc=1,NAT
              Akt(i,j,k,itrc)=Akt_bak(itrc,ng)+ql*Sh
            END DO
!
!  Compute vertical mixing (m2/s) coefficient of turbulent kinetic
!  energy.  Use original formulation (Mellor and Yamada 1982;
!  Blumberg 1991; Kantha and Clayson 1994).
!
            Akk(i,j,k)=Akk_bak(ng)+ql*my_Sq
!
!  Save limited length scale.
!
            Lscale(i,j,k)=Ls_lmt
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL tkebc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nnew, nstp,                                      &
     &                 gls, tke)
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tke(:,:,:,nnew))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          gls(:,:,:,nnew))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akv)
        DO itrc=1,NAT
          CALL exchange_w3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            Akt(:,:,:,itrc))
        END DO
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tke(:,:,:,nnew),                              &
     &                    gls(:,:,:,nnew),                              &
     &                    Akv)
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng), 1, NAT,         &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akt)
# endif

      RETURN
      END SUBROUTINE my25_corstep_tile
#endif
      END MODULE my25_corstep_mod
