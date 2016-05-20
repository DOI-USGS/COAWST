      SUBROUTINE prsgrd (ng, tile)
!
!svn $Id: prsgrd32.h 732 2008-09-07 01:55:51Z jcwarner $
!***********************************************************************
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This subroutine evaluates the nonlinear  baroclinic,  hydrostatic   !
!  pressure gradient term using a  nonconservative  Density-Jacobian   !
!  scheme,  based on  cubic polynomial fits for  "rho" and  "z_r" as   !
!  functions of nondimensional coordinates (XI,ETA,s), that is,  its   !
!  respective array indices. The  cubic polynomials  are monotonized   !
!  by using  harmonic mean instead of linear averages to interpolate   !
!  slopes. This scheme retains exact anti-symmetry:                    !
!                                                                      !
!        J(rho,z_r)=-J(z_r,rho).                                       !
!                                                                      !
!  If parameter OneFifth (below) is set to zero,  the scheme becomes   !
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
      USE mod_diags
#endif
#ifdef ATM_PRESS
      USE mod_forces
#endif
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#ifdef POT_TIDES
      USE mod_tides
#endif
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
      CALL wclock_on (ng, iNLM, 23)
#endif
      CALL prsgrd_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  nrhs(ng),                                       &
#ifdef MASKING
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % vmask,                               &
#endif
#ifdef WET_DRY
     &                  GRID(ng) % umask_wet,                           &
     &                  GRID(ng) % vmask_wet,                           &
#endif
     &                  GRID(ng) % om_v,                                &
     &                  GRID(ng) % on_u,                                &
     &                  GRID(ng) % Hz,                                  &
     &                  GRID(ng) % z_r,                                 &
     &                  GRID(ng) % z_w,                                 &
#ifdef ICESHELF
     &                  GRID(ng) % zice,                                &
#endif
     &                  OCEAN(ng) % rho,                                &
#ifdef WEC_VF
     &                  OCEAN(ng) % zetat,                              &
#endif

#ifdef ATM_PRESS
     &                  FORCES(ng) % Pair,                              &
#endif
#ifdef POT_TIDES
     &                  TIDES(ng) % Ptide,                              &
#endif
#ifdef DIAGNOSTICS_UV
     &                  DIAGS(ng) % DiaRU,                              &
     &                  DIAGS(ng) % DiaRV,                              &
#endif
     &                  OCEAN(ng) % ru,                                 &
     &                  OCEAN(ng) % rv)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 23)
#endif
      RETURN
      END SUBROUTINE prsgrd
!
!***********************************************************************
      SUBROUTINE prsgrd_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs,                                     &
#ifdef MASKING
     &                        umask, vmask,                             &
#endif
#ifdef WET_DRY
     &                        umask_wet, vmask_wet,                     &
#endif
     &                        om_v, on_u,                               &
     &                        Hz, z_r, z_w,                             &
# ifdef ICESHELF
     &                        zice,                                     &
# endif
     &                        rho,                                      &
#ifdef WEC_VF
     &                        zetat,                                    &
#endif
#ifdef ATM_PRESS
     &                        Pair,                                     &
#endif
#ifdef POT_TIDES
     &                        Ptide,                                    &
#endif
#ifdef DIAGNOSTICS_UV
     &                        DiaRU, DiaRV,                             &
#endif
     &                        ru, rv)
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
# ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
# endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
# endif

      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#ifdef WEC_VF
      real(r8), intent(in) :: zetat(LBi:,LBj:)
#endif
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:,LBj:)
# endif
# ifdef POT_TIDES
      real(r8), intent(in) :: Ptide(LBi:,LBj:)
# endif
# ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
#ifdef WEC_VF
      real(r8), intent(in) :: zetat(LBi:UBi,LBj:UBj)
#endif
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
# endif
# ifdef POT_TIDES
      real(r8), intent(in) :: Ptide(LBi:UBi,LBj:UBj)
# endif
# ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
# endif
      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: OneFifth = 0.2_r8
      real(r8), parameter :: OneTwelfth = 1.0_r8/12.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8
#ifdef ICESHELF
      real(r8), parameter :: drhodz = 0.00478_r8
#endif

      real(r8) :: GRho, GRho0,  HalfGRho
      real(r8) :: cff, cff1, cff2
#ifdef ATM_PRESS
      real(r8) :: OneAtm, fac
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: P

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dZ

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dRx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dZx
!
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
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
          END DO
        END DO
        DO i=IstrU-1,Iend
          dR(i,N(ng))=dR(i,N(ng)-1)
          dZ(i,N(ng))=dZ(i,N(ng)-1)
          dR(i,0)=dR(i,1)
          dZ(i,0)=dZ(i,1)
        END DO
        DO k=N(ng),1,-1
          DO i=IstrU-1,Iend
            cff=2.0_r8*dR(i,k)*dR(i,k-1)
            IF (cff.gt.eps) THEN
              dR(i,k)=cff/(dR(i,k)+dR(i,k-1))
            ELSE
              dR(i,k)=0.0_r8
            END IF
            dZ(i,k)=2.0_r8*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1))
          END DO
        END DO
        DO i=IstrU-1,Iend
          cff1=1.0_r8/(z_r(i,j,N(ng))-z_r(i,j,N(ng)-1))
          cff2=0.5_r8*(rho(i,j,N(ng))-rho(i,j,N(ng)-1))*                &
     &         (z_w(i,j,N(ng))-z_r(i,j,N(ng)))*cff1
#ifdef ICESHELF
          P(i,j,N(ng))=GRho0*(z_w(i,j,N(ng))-zice(i,j))-                &
     &                 GRho*(rho(i,j,N(ng))+0.5_r8*drhodz*zice(i,j))*   &
     &                 zice(i,j)+                                       &
     &                 GRho*(rho(i,j,N(ng))+cff2)*                      &
     &                 (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
#else
          P(i,j,N(ng))=GRho0*z_w(i,j,N(ng))+                            &
#ifdef WEC_VF
     &                 zetat(i,j)+                                      &
#endif
#ifdef ATM_PRESS
     &                 fac*(Pair(i,j)-OneAtm)+                          &
#endif
     &                 GRho*(rho(i,j,N(ng))+cff2)*                      &
     &                 (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
#endif
#ifdef POT_TIDES
          P(i,j,N(ng)) = P(i,j,N(ng)) - g*Ptide(i,j)
#endif
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU-1,Iend
            P(i,j,k)=P(i,j,k+1)+                                        &
     &               HalfGRho*((rho(i,j,k+1)+rho(i,j,k))*               &
     &                         (z_r(i,j,k+1)-z_r(i,j,k))-               &
     &                         OneFifth*                                &
     &                         ((dR(i,k+1)-dR(i,k))*                    &
     &                          (z_r(i,j,k+1)-z_r(i,j,k)-               &
     &                           OneTwelfth*                            &
     &                           (dZ(i,k+1)+dZ(i,k)))-                  &
     &                          (dZ(i,k+1)-dZ(i,k))*                    &
     &                          (rho(i,j,k+1)-rho(i,j,k)-               &
     &                           OneTwelfth*                            &
     &                           (dR(i,k+1)+dR(i,k)))))
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
#ifdef MASKING
            aux(i,j)=aux(i,j)*umask(i,j)
#endif
            FC(i,j)=rho(i,j,k)-rho(i-1,j,k)
#ifdef MASKING
            FC(i,j)=FC(i,j)*umask(i,j)
#endif
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff=2.0_r8*aux(i,j)*aux(i+1,j)
            IF (cff.gt.eps) THEN
              cff1=1.0_r8/(aux(i,j)+aux(i+1,j))
              dZx(i,j)=cff*cff1
            ELSE
              dZx(i,j)=0.0_r8
            END IF
            cff1=2.0_r8*FC(i,j)*FC(i+1,j)
            IF (cff1.gt.eps) THEN
              cff2=1.0_r8/(FC(i,j)+FC(i+1,j))
              dRx(i,j)=cff1*cff2
            ELSE
              dRx(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            ru(i,j,k,nrhs)=on_u(i,j)*0.5_r8*                            &
     &                     (Hz(i,j,k)+Hz(i-1,j,k))*                     &
     &                     (P(i-1,j,k)-P(i,j,k)-                        &
     &                      HalfGRho*                                   &
     &                      ((rho(i,j,k)+rho(i-1,j,k))*                 &
     &                       (z_r(i,j,k)-z_r(i-1,j,k))-                 &
     &                        OneFifth*                                 &
     &                        ((dRx(i,j)-dRx(i-1,j))*                   &
     &                         (z_r(i,j,k)-z_r(i-1,j,k)-                &
     &                          OneTwelfth*                             &
     &                          (dZx(i,j)+dZx(i-1,j)))-                 &
     &                         (dZx(i,j)-dZx(i-1,j))*                   &
     &                         (rho(i,j,k)-rho(i-1,j,k)-                &
     &                          OneTwelfth*                             &
     &                          (dRx(i,j)+dRx(i-1,j))))))
#ifdef WET_DRY
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)*umask_wet(i,j)
#endif
#ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
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
#ifdef MASKING
            aux(i,j)=aux(i,j)*vmask(i,j)
#endif
            FC(i,j)=rho(i,j,k)-rho(i,j-1,k)
#ifdef MASKING
            FC(i,j)=FC(i,j)*vmask(i,j)
#endif
          END DO
        END DO
!
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff=2.0_r8*aux(i,j)*aux(i,j+1)
            IF (cff.gt.eps) THEN
              cff1=1.0_r8/(aux(i,j)+aux(i,j+1))
              dZx(i,j)=cff*cff1
            ELSE
              dZx(i,j)=0.0_r8
            END IF
            cff1=2.0_r8*FC(i,j)*FC(i,j+1)
            IF (cff1.gt.eps) THEN
              cff2=1.0_r8/(FC(i,j)+FC(i,j+1))
              dRx(i,j)=cff1*cff2
            ELSE
              dRx(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            rv(i,j,k,nrhs)=om_v(i,j)*0.5_r8*                            &
     &                     (Hz(i,j,k)+Hz(i,j-1,k))*                     &
     &                     (P(i,j-1,k)-P(i,j,k)-                        &
     &                      HalfGRho*                                   &
     &                      ((rho(i,j,k)+rho(i,j-1,k))*                 &
     &                       (z_r(i,j,k)-z_r(i,j-1,k))-                 &
     &                        OneFifth*                                 &
     &                        ((dRx(i,j)-dRx(i,j-1))*                   &
     &                         (z_r(i,j,k)-z_r(i,j-1,k)-                &
     &                          OneTwelfth*                             &
     &                          (dZx(i,j)+dZx(i,j-1)))-                 &
     &                         (dZx(i,j)-dZx(i,j-1))*                   &
     &                         (rho(i,j,k)-rho(i,j-1,k)-                &
     &                          OneTwelfth*                             &
     &                          (dRx(i,j)+dRx(i,j-1))))))
#ifdef WET_DRY
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)*vmask_wet(i,j)
#endif
#ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
#endif
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE prsgrd_tile
