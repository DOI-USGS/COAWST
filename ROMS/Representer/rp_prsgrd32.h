      SUBROUTINE rp_prsgrd (ng, tile)
!
!svn $Id: rp_prsgrd32.h 995 2020-01-10 04:01:28Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2020 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine evaluates the representers tangent linear baroclinic,  !
!  hydrostatic pressure gradient term using a nonconservative Density  !
!  Jacobian scheme,  based on  cubic polynomial fits  for  "rho"  and  !
!  "z_r" as functions of nondimensional coordinates (XI,ETA,s),  that  !
!  is,  its respective  array  indices.  The  cubic  polynomials  are  !
!  monotonized  by using  harmonic mean instead of linear averages to  !
!  interpolate slopes. This scheme retains exact anti-symmetry:        !
!                                                                      !
!        J(rho,z_r)=-J(z_r,rho).                                       !
!                                                                      !
!  If parameter OneFifth  (below) is set to zero,  the scheme becomes  !
!  identical to standard Jacobian.                                     !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Shchepetkin A.F and J.C. McWilliams, 2003:  A method for          !
!      computing horizontal pressure gradient force in an ocean        !
!      model with non-aligned vertical coordinate, JGR, 108,           !
!      1-34.                                                           !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef DIAGNOSTICS
!!    USE mod_diags
#endif
#ifdef ATM_PRESS
      USE mod_forces
#endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iRPM, 23, __LINE__, __FILE__)
#endif
      CALL rp_prsgrd_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#endif
     &                     GRID(ng) % om_v,                             &
     &                     GRID(ng) % on_u,                             &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % tl_Hz,                            &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % tl_z_r,                           &
     &                     GRID(ng) % z_w,                              &
     &                     GRID(ng) % tl_z_w,                           &
     &                     OCEAN(ng) % rho,                             &
     &                     OCEAN(ng) % tl_rho,                          &
#ifdef ATM_PRESS
     &                     FORCES(ng) % Pair,                           &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                     DIAGS(ng) % DiaRU,                           &
!!   &                     DIAGS(ng) % DiaRV,                           &
#endif
     &                     OCEAN(ng) % tl_ru,                           &
     &                     OCEAN(ng) % tl_rv)
#ifdef PROFILE
      CALL wclock_off (ng, iRPM, 23, __LINE__, __FILE__)
#endif
      RETURN
      END SUBROUTINE rp_prsgrd
!
!***********************************************************************
      SUBROUTINE rp_prsgrd_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
#ifdef MASKING
     &                           umask, vmask,                          &
#endif
     &                           om_v, on_u,                            &
     &                           Hz, tl_Hz,                             &
     &                           z_r, tl_z_r,                           &
     &                           z_w, tl_z_w,                           &
     &                           rho, tl_rho,                           &
#ifdef ATM_PRESS
     &                           Pair,                                  &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                           DiaRU, DiaRV,                          &
#endif
     &                           tl_ru, tl_rv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)

      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_rho(LBi:,LBj:,:)
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:,LBj:)
# endif
# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
!!    real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: tl_ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tl_rv(LBi:,LBj:,0:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_rho(LBi:UBi,LBj:UBj,N(ng))
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
# endif
# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
!!    real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
# endif
      real(r8), intent(inout) :: tl_ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: tl_rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: OneFifth = 0.2_r8
      real(r8), parameter :: OneTwelfth = 1.0_r8/12.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8) :: GRho, GRho0, HalfGRho
      real(r8) :: cff, cff1, cff2
      real(r8) :: tl_cff, tl_cff1, tl_cff2
#ifdef ATM_PRESS
      real(r8) :: OneAtm, fac
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: P

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: tl_P

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR1
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dZ
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dZ1

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_dR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_dZ

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dRx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dZx

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_aux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_dRx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_dZx

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Preliminary step (same for XI- and ETA-components:
!-----------------------------------------------------------------------
!
      GRho=g/rho0
      GRho0=1000.0_r8*GRho
      HalfGRho=0.5_r8*GRho
#ifdef ATM_PRESS
      OneAtm=1013.25_r8                  ! 1 atm = 1013.25 mb
      fac=100.0_r8/rho0
#endif
!
      DO j=JstrV-1,Jend
        DO k=1,N(ng)-1
          DO i=IstrU-1,Iend
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
            tl_dR(i,k)=tl_rho(i,j,k+1)-tl_rho(i,j,k)
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
            tl_dZ(i,k)=tl_z_r(i,j,k+1)-tl_z_r(i,j,k)
          END DO
        END DO
        DO i=IstrU-1,Iend
          dR(i,N(ng))=dR(i,N(ng)-1)
          tl_dR(i,N(ng))=tl_dR(i,N(ng)-1)
          dZ(i,N(ng))=dZ(i,N(ng)-1)
          tl_dZ(i,N(ng))=tl_dZ(i,N(ng)-1)
          dR(i,0)=dR(i,1)
          tl_dR(i,0)=tl_dR(i,1)
          dZ(i,0)=dZ(i,1)
          tl_dZ(i,0)=tl_dZ(i,1)
        END DO
        DO k=N(ng),1,-1
          DO i=IstrU-1,Iend
            cff=2.0_r8*dR(i,k)*dR(i,k-1)
            tl_cff=2.0_r8*(tl_dR(i,k)*dR(i,k-1)+                        &
     &                     dR(i,k)*tl_dR(i,k-1))-                       &
#ifdef TL_IOMS
     &             cff
#endif
            dR1(i,k)=dR(i,k)
            IF (cff.gt.eps) THEN
              dR(i,k)=cff/(dR(i,k)+dR(i,k-1))
              tl_dR(i,k)=(tl_cff-dR(i,k)*(tl_dR(i,k)+tl_dR(i,k-1)))/    &
     &                   (dR1(i,k)+dR(i,k-1))+                          &
#ifdef TL_IOMS
     &                   dR(i,k)
#endif
            ELSE
              dR(i,k)=0.0_r8
              tl_dR(i,k)=0.0_r8
            END IF
            dZ1(i,k)=dZ(i,k)
            dZ(i,k)=2.0_r8*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1))
            tl_dZ(i,k)=(2.0_r8*(tl_dZ(i,k)*dZ(i,k-1)+                   &
     &                          dZ1(i,k)*tl_dZ(i,k-1))-                 &
     &                  dZ(i,k)*(tl_dZ(i,k)+tl_dZ(i,k-1)))/             &
     &                 (dZ1(i,k)+dZ(i,k-1))
          END DO
        END DO
        DO i=IstrU-1,Iend
          cff1=1.0_r8/(z_r(i,j,N(ng))-z_r(i,j,N(ng)-1))
          tl_cff1=-cff1*cff1*(tl_z_r(i,j,N(ng))-tl_z_r(i,j,N(ng)-1))+   &
#ifdef TL_IOMS
     &            2.0_r8*cff1
#endif
          cff2=0.5_r8*(rho(i,j,N(ng))-rho(i,j,N(ng)-1))*                &
     &         (z_w(i,j,N(ng))-z_r(i,j,N(ng)))*cff1
          tl_cff2=0.5_r8*((tl_rho(i,j,N(ng))-tl_rho(i,j,N(ng)-1))*      &
     &                    (z_w(i,j,N(ng))-z_r(i,j,N(ng)))*cff1+         &
     &                    (rho(i,j,N(ng))-rho(i,j,N(ng)-1))*            &
     &                    ((tl_z_w(i,j,N(ng))-tl_z_r(i,j,N(ng)))*cff1+  &
     &                     (z_w(i,j,N(ng))-z_r(i,j,N(ng)))*tl_cff1))-   &
#ifdef TL_IOMS
     &            2.0_r8*cff2
#endif
          P(i,j,N(ng))=GRho0*z_w(i,j,N(ng))+                            &
#ifdef ATM_PRESS
     &                 fac*(Pair(i,j)-OneAtm)+                          &
#endif
     &                 GRho*(rho(i,j,N(ng))+cff2)*                      &
     &                 (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
          tl_P(i,j,N(ng))=GRho0*tl_z_w(i,j,N(ng))+                      &
     &                    GRho*((tl_rho(i,j,N(ng))+tl_cff2)*            &
     &                          (z_w(i,j,N(ng))-z_r(i,j,N(ng)))+        &
     &                          (rho(i,j,N(ng))+cff2)*                  &
     &                          (tl_z_w(i,j,N(ng))-tl_z_r(i,j,N(ng))))- &
#ifdef TL_IOMS
     &                    GRho*(rho(i,j,N(ng))+cff2)*                   &
     &                    (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
#endif
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU-1,Iend
            cff=HalfGRho*((rho(i,j,k+1)+rho(i,j,k))*                    &
     &                    (z_r(i,j,k+1)-z_r(i,j,k))-                    &
     &                    OneFifth*                                     &
     &                    ((dR(i,k+1)-dR(i,k))*                         &
     &                     (z_r(i,j,k+1)-z_r(i,j,k)-                    &
     &                      OneTwelfth*                                 &
     &                      (dZ(i,k+1)+dZ(i,k)))-                       &
     &                     (dZ(i,k+1)-dZ(i,k))*                         &
     &                     (rho(i,j,k+1)-rho(i,j,k)-                    &
     &                      OneTwelfth*                                 &
     &                      (dR(i,k+1)+dR(i,k)))))
            tl_cff=HalfGRho*((tl_rho(i,j,k+1)+tl_rho(i,j,k))*           &
     &                       (z_r(i,j,k+1)-z_r(i,j,k))+                 &
     &                       (rho(i,j,k+1)+rho(i,j,k))*                 &
     &                       (tl_z_r(i,j,k+1)-tl_z_r(i,j,k))-           &
     &                       OneFifth*                                  &
     &                       ((tl_dR(i,k+1)-tl_dR(i,k))*                &
     &                        (z_r(i,j,k+1)-z_r(i,j,k)-                 &
     &                         OneTwelfth*                              &
     &                         (dZ(i,k+1)+dZ(i,k)))+                    &
     &                        (dR(i,k+1)-dR(i,k))*                      &
     &                        (tl_z_r(i,j,k+1)-tl_z_r(i,j,k)-           &
     &                         OneTwelfth*                              &
     &                         (tl_dZ(i,k+1)+tl_dZ(i,k)))-              &
     &                        (tl_dZ(i,k+1)-tl_dZ(i,k))*                &
     &                        (rho(i,j,k+1)-rho(i,j,k)-                 &
     &                         OneTwelfth*                              &
     &                         (dR(i,k+1)+dR(i,k)))-                    &
     &                        (dZ(i,k+1)-dZ(i,k))*                      &
     &                        (tl_rho(i,j,k+1)-tl_rho(i,j,k)-           &
     &                         OneTwelfth*                              &
     &                         (tl_dR(i,k+1)+tl_dR(i,k)))))-            &
#ifdef TL_IOMS
     &             cff
#endif
            P(i,j,k)=P(i,j,k+1)+cff
            tl_P(i,j,k)=tl_P(i,j,k+1)+tl_cff
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute XI-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=N(ng),1,-1
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend+1
            aux(i,j)=z_r(i,j,k)-z_r(i-1,j,k)
            tl_aux(i,j)=tl_z_r(i,j,k)-tl_z_r(i-1,j,k)
#ifdef MASKING
            aux(i,j)=aux(i,j)*umask(i,j)
            tl_aux(i,j)=tl_aux(i,j)*umask(i,j)
#endif
            FC(i,j)=rho(i,j,k)-rho(i-1,j,k)
            tl_FC(i,j)=tl_rho(i,j,k)-tl_rho(i-1,j,k)
#ifdef MASKING
            FC(i,j)=FC(i,j)*umask(i,j)
            tl_FC(i,j)=tl_FC(i,j)*umask(i,j)
#endif
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff=2.0_r8*aux(i,j)*aux(i+1,j)
            tl_cff=2.0_r8*(tl_aux(i,j)*aux(i+1,j)+                      &
     &                     aux(i,j)*tl_aux(i+1,j))-                     &
#ifdef TL_IOMS
     &             cff
#endif
            IF (cff.gt.eps) THEN
              cff1=1.0_r8/(aux(i,j)+aux(i+1,j))
              tl_cff1=-cff1*cff1*(tl_aux(i,j)+tl_aux(i+1,j))+           &
#ifdef TL_IOMS
     &                2.0_r8*cff1
#endif
              dZx(i,j)=cff*cff1
              tl_dZx(i,j)=tl_cff*cff1+cff*tl_cff1-                      &
#ifdef TL_IOMS
     &                    dZx(i,j)
#endif
            ELSE
              dZx(i,j)=0.0_r8
              tl_dZx(i,j)=0.0_r8
            END IF
            cff1=2.0_r8*FC(i,j)*FC(i+1,j)
            tl_cff1=2.0_r8*(tl_FC(i,j)*FC(i+1,j)+                       &
     &                      FC(i,j)*tl_FC(i+1,j))-                      &
#ifdef TL_IOMS
     &              cff1
#endif
            IF (cff1.gt.eps) THEN
              cff2=1.0_r8/(FC(i,j)+FC(i+1,j))
              tl_cff2=-cff2*cff2*(tl_FC(i,j)+tl_FC(i+1,j))+             &
#ifdef TL_IOMS
     &                2.0_r8*cff2
#endif
              dRx(i,j)=cff1*cff2
              tl_dRx(i,j)=tl_cff1*cff2+cff1*tl_cff2-                    &
#ifdef TL_IOMS
     &                    dRx(i,j)
#endif
            ELSE
              dRx(i,j)=0.0_r8
              tl_dRx(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!>          ru(i,j,k,nrhs)=on_u(i,j)*0.5_r8*                            &
!>   &                     (Hz(i,j,k)+Hz(i-1,j,k))*                     &
!>   &                     (P(i-1,j,k)-P(i,j,k)-                        &
!>   &                      HalfGRho*                                   &
!>   &                      ((rho(i,j,k)+rho(i-1,j,k))*                 &
!>   &                       (z_r(i,j,k)-z_r(i-1,j,k))-                 &
!>   &                        OneFifth*                                 &
!>   &                        ((dRx(i,j)-dRx(i-1,j))*                   &
!>   &                         (z_r(i,j,k)-z_r(i-1,j,k)-                &
!>   &                          OneTwelfth*                             &
!>   &                          (dZx(i,j)+dZx(i-1,j)))-                 &
!>   &                         (dZx(i,j)-dZx(i-1,j))*                   &
!>   &                         (rho(i,j,k)-rho(i-1,j,k)-                &
!>   &                          OneTwelfth*                             &
!>   &                          (dRx(i,j)+dRx(i-1,j))))))
!>
            tl_ru(i,j,k,nrhs)=on_u(i,j)*0.5_r8*                         &
     &                        ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*           &
     &                         (P(i-1,j,k)-P(i,j,k)-                    &
     &                          HalfGRho*                               &
     &                          ((rho(i,j,k)+rho(i-1,j,k))*             &
     &                           (z_r(i,j,k)-z_r(i-1,j,k))-             &
     &                            OneFifth*                             &
     &                            ((dRx(i,j)-dRx(i-1,j))*               &
     &                             (z_r(i,j,k)-z_r(i-1,j,k)-            &
     &                              OneTwelfth*                         &
     &                              (dZx(i,j)+dZx(i-1,j)))-             &
     &                             (dZx(i,j)-dZx(i-1,j))*               &
     &                             (rho(i,j,k)-rho(i-1,j,k)-            &
     &                              OneTwelfth*                         &
     &                              (dRx(i,j)+dRx(i-1,j))))))+          &
     &                         (Hz(i,j,k)+Hz(i-1,j,k))*                 &
     &                         (tl_P(i-1,j,k)-tl_P(i,j,k)-              &
     &                          HalfGRho*                               &
     &                          ((tl_rho(i,j,k)+tl_rho(i-1,j,k))*       &
     &                           (z_r(i,j,k)-z_r(i-1,j,k))+             &
     &                           (rho(i,j,k)+rho(i-1,j,k))*             &
     &                           (tl_z_r(i,j,k)-tl_z_r(i-1,j,k))-       &
     &                            OneFifth*                             &
     &                            ((tl_dRx(i,j)-tl_dRx(i-1,j))*         &
     &                             (z_r(i,j,k)-z_r(i-1,j,k)-            &
     &                              OneTwelfth*                         &
     &                              (dZx(i,j)+dZx(i-1,j)))+             &
     &                             (dRx(i,j)-dRx(i-1,j))*               &
     &                             (tl_z_r(i,j,k)-tl_z_r(i-1,j,k)-      &
     &                              OneTwelfth*                         &
     &                              (tl_dZx(i,j)+tl_dZx(i-1,j)))-       &
     &                             (tl_dZx(i,j)-tl_dZx(i-1,j))*         &
     &                             (rho(i,j,k)-rho(i-1,j,k)-            &
     &                              OneTwelfth*                         &
     &                              (dRx(i,j)+dRx(i-1,j)))-             &
     &                             (dZx(i,j)-dZx(i-1,j))*               &
     &                             (tl_rho(i,j,k)-tl_rho(i-1,j,k)-      &
     &                              OneTwelfth*                         &
     &                              (tl_dRx(i,j)+tl_dRx(i-1,j)))))))-   &
#ifdef TL_IOMS
     &                        on_u(i,j)*0.5_r8*                         &
     &                        (Hz(i,j,k)+Hz(i-1,j,k))*                  &
     &                        (P(i-1,j,k)-P(i,j,k)-                     &
     &                         2.0_r8*HalfGRho*                         &
     &                         ((rho(i,j,k)+rho(i-1,j,k))*              &
     &                          (z_r(i,j,k)-z_r(i-1,j,k))-              &
     &                           OneFifth*                              &
     &                           ((dRx(i,j)-dRx(i-1,j))*                &
     &                            (z_r(i,j,k)-z_r(i-1,j,k)-             &
     &                             OneTwelfth*                          &
     &                             (dZx(i,j)+dZx(i-1,j)))-              &
     &                            (dZx(i,j)-dZx(i-1,j))*                &
     &                            (rho(i,j,k)-rho(i-1,j,k)-             &
     &                             OneTwelfth*                          &
     &                             (dRx(i,j)+dRx(i-1,j))))))
#endif
#ifdef DIAGNOSTICS_UV
!!          DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
#endif
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  ETA-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=N(ng),1,-1
        DO j=JstrV-1,Jend+1
          DO i=Istr,Iend
            aux(i,j)=z_r(i,j,k)-z_r(i,j-1,k)
            tl_aux(i,j)=tl_z_r(i,j,k)-tl_z_r(i,j-1,k)
#ifdef MASKING
            aux(i,j)=aux(i,j)*vmask(i,j)
            tl_aux(i,j)=tl_aux(i,j)*vmask(i,j)
#endif
            FC(i,j)=rho(i,j,k)-rho(i,j-1,k)
            tl_FC(i,j)=tl_rho(i,j,k)-tl_rho(i,j-1,k)
#ifdef MASKING
            FC(i,j)=FC(i,j)*vmask(i,j)
            tl_FC(i,j)=tl_FC(i,j)*vmask(i,j)
#endif
          END DO
        END DO
!
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff=2.0_r8*aux(i,j)*aux(i,j+1)
            tl_cff=2.0_r8*(tl_aux(i,j)*aux(i,j+1)+                      &
     &                     aux(i,j)*tl_aux(i,j+1))-                     &
#ifdef TL_IOMS
     &             cff
#endif
            IF (cff.gt.eps) THEN
              cff1=1.0_r8/(aux(i,j)+aux(i,j+1))
              tl_cff1=-cff1*cff1*(tl_aux(i,j)+tl_aux(i,j+1))+           &
#ifdef TL_IOMS
     &                2.0_r8*cff1
#endif
              dZx(i,j)=cff*cff1
              tl_dZx(i,j)=tl_cff*cff1+cff*tl_cff1-                      &
#ifdef TL_IOMS
     &                    dZx(i,j)
#endif
            ELSE
              dZx(i,j)=0.0_r8
              tl_dZx(i,j)=0.0_r8
            END IF
            cff1=2.0_r8*FC(i,j)*FC(i,j+1)
            tl_cff1=2.0_r8*(tl_FC(i,j)*FC(i,j+1)+                       &
     &                      FC(i,j)*tl_FC(i,j+1))-                      &
#ifdef TL_IOMS
     &              cff1
#endif
            IF (cff1.gt.eps) THEN
              cff2=1.0_r8/(FC(i,j)+FC(i,j+1))
              tl_cff2=-cff2*cff2*(tl_FC(i,j)+tl_FC(i,j+1))+             &
#ifdef TL_IOMS
     &                2.0_r8*cff2
#endif
              dRx(i,j)=cff1*cff2
              tl_dRx(i,j)=tl_cff1*cff2+cff1*tl_cff2-                    &
#ifdef TL_IOMS
     &                    dRx(i,j)
#endif
            ELSE
              dRx(i,j)=0.0_r8
              tl_dRx(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
!>          rv(i,j,k,nrhs)=om_v(i,j)*0.5_r8*                            &
!>   &                     (Hz(i,j,k)+Hz(i,j-1,k))*                     &
!>   &                     (P(i,j-1,k)-P(i,j,k)-                        &
!>   &                      HalfGRho*                                   &
!>   &                      ((rho(i,j,k)+rho(i,j-1,k))*                 &
!>   &                       (z_r(i,j,k)-z_r(i,j-1,k))-                 &
!>   &                        OneFifth*                                 &
!>   &                        ((dRx(i,j)-dRx(i,j-1))*                   &
!>   &                         (z_r(i,j,k)-z_r(i,j-1,k)-                &
!>   &                          OneTwelfth*                             &
!>   &                          (dZx(i,j)+dZx(i,j-1)))-                 &
!>   &                         (dZx(i,j)-dZx(i,j-1))*                   &
!>   &                         (rho(i,j,k)-rho(i,j-1,k)-                &
!>   &                          OneTwelfth*                             &
!>   &                          (dRx(i,j)+dRx(i,j-1))))))
!>
            tl_rv(i,j,k,nrhs)=om_v(i,j)*0.5_r8*                         &
     &                        ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*           &
     &                         (P(i,j-1,k)-P(i,j,k)-                    &
     &                          HalfGRho*                               &
     &                          ((rho(i,j,k)+rho(i,j-1,k))*             &
     &                           (z_r(i,j,k)-z_r(i,j-1,k))-             &
     &                           OneFifth*                              &
     &                           ((dRx(i,j)-dRx(i,j-1))*                &
     &                            (z_r(i,j,k)-z_r(i,j-1,k)-             &
     &                             OneTwelfth*                          &
     &                             (dZx(i,j)+dZx(i,j-1)))-              &
     &                            (dZx(i,j)-dZx(i,j-1))*                &
     &                            (rho(i,j,k)-rho(i,j-1,k)-             &
     &                             OneTwelfth*                          &
     &                             (dRx(i,j)+dRx(i,j-1))))))+           &
     &                         (Hz(i,j,k)+Hz(i,j-1,k))*                 &
     &                         (tl_P(i,j-1,k)-tl_P(i,j,k)-              &
     &                          HalfGRho*                               &
     &                          ((tl_rho(i,j,k)+tl_rho(i,j-1,k))*       &
     &                           (z_r(i,j,k)-z_r(i,j-1,k))+             &
     &                           (rho(i,j,k)+rho(i,j-1,k))*             &
     &                           (tl_z_r(i,j,k)-tl_z_r(i,j-1,k))-       &
     &                           OneFifth*                              &
     &                           ((tl_dRx(i,j)-tl_dRx(i,j-1))*          &
     &                            (z_r(i,j,k)-z_r(i,j-1,k)-             &
     &                             OneTwelfth*                          &
     &                             (dZx(i,j)+dZx(i,j-1)))+              &
     &                            (dRx(i,j)-dRx(i,j-1))*                &
     &                            (tl_z_r(i,j,k)-tl_z_r(i,j-1,k)-       &
     &                             OneTwelfth*                          &
     &                             (tl_dZx(i,j)+tl_dZx(i,j-1)))-        &
     &                            (tl_dZx(i,j)-tl_dZx(i,j-1))*          &
     &                            (rho(i,j,k)-rho(i,j-1,k)-             &
     &                             OneTwelfth*                          &
     &                             (dRx(i,j)+dRx(i,j-1)))-              &
     &                            (dZx(i,j)-dZx(i,j-1))*                &
     &                            (tl_rho(i,j,k)-tl_rho(i,j-1,k)-       &
     &                             OneTwelfth*                          &
     &                             (tl_dRx(i,j)+tl_dRx(i,j-1)))))))-    &
#ifdef TL_IOMS
     &                        om_v(i,j)*0.5_r8*                         &
     &                        (Hz(i,j,k)+Hz(i,j-1,k))*                  &
     &                        (P(i,j-1,k)-P(i,j,k)-                     &
     &                         2.0_r8*HalfGRho*                         &
     &                         ((rho(i,j,k)+rho(i,j-1,k))*              &
     &                          (z_r(i,j,k)-z_r(i,j-1,k))-              &
     &                           OneFifth*                              &
     &                           ((dRx(i,j)-dRx(i,j-1))*                &
     &                            (z_r(i,j,k)-z_r(i,j-1,k)-             &
     &                             OneTwelfth*                          &
     &                             (dZx(i,j)+dZx(i,j-1)))-              &
     &                            (dZx(i,j)-dZx(i,j-1))*                &
     &                            (rho(i,j,k)-rho(i,j-1,k)-             &
     &                             OneTwelfth*                          &
     &                             (dRx(i,j)+dRx(i,j-1))))))
#endif
#ifdef DIAGNOSTICS_UV
!!          DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
#endif
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE rp_prsgrd_tile
