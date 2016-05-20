# undef NEUMANN
      SUBROUTINE prsgrd (ng, tile)
!
!svn $Id: prsgrd44.h 732 2008-09-07 01:55:51Z jcwarner $
!***********************************************************************
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This subroutine evaluates the baroclinic,  hydrostatic pressure     !
!  gradient term using a  finite-volume  pressure Jacobian scheme.     !
!  The scheme is based on local, conservative, limited-oscillation     !
!  vertical quartic polynomial reconstruction of density field and     !
!  subsequent  projection fits of the derivatives of  density into     !
!  isosurface of vertical coordinate.  The monotonicity constraint     !
!  uses a  PPM-style  limitting  algorithm with a  power-law slope     !
!  reconciliation step.                                                !
!                                                                      !
!  isosurface of vertical coordinate.                                  !
!                                                                      !
!  The pressure gradient terms (m4/s2) are loaded into right-hand-     !
!  side arrays "ru" and "rv".                                          !
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
# ifdef DIAGNOSTICS
      USE mod_diags
# endif
# ifdef ATM_PRESS
      USE mod_forces
# endif
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 23)
# endif
      CALL prsgrd_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  nrhs(ng),                                       &
# ifdef WET_DRY
     &                  GRID(ng)%umask_wet,                             &
     &                  GRID(ng)%vmask_wet,                             &
# endif
     &                  GRID(ng) % Hz,                                  &
     &                  GRID(ng) % om_v,                                &
     &                  GRID(ng) % on_u,                                &
     &                  GRID(ng) % z_w,                                 &
     &                  OCEAN(ng) % rho,                                &
# ifdef WEC_VF
     &                  OCEAN(ng) % zetat,                              &
# endif
# ifdef ATM_PRESS
     &                  FORCES(ng) % Pair,                              &
# endif
# ifdef DIAGNOSTICS_UV
     &                  DIAGS(ng) % DiaRU,                              &
     &                  DIAGS(ng) % DiaRV,                              &
# endif
     &                  OCEAN(ng) % ru,                                 &
     &                  OCEAN(ng) % rv)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 23)
# endif
      RETURN
      END SUBROUTINE prsgrd
!
!***********************************************************************
      SUBROUTINE prsgrd_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs,                                     &
# ifdef WET_DRY
     &                        umask_wet, vmask_wet,                     &
# endif
     &                        Hz, om_v, on_u, z_w,                      &
     &                        rho,                                      &
# ifdef WEC_VF
     &                        zetat,                                    &
# endif
# ifdef ATM_PRESS
     &                        Pair,                                     &
# endif
# ifdef DIAGNOSTICS_UV
     &                        DiaRU, DiaRV,                             &
# endif
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
!
# ifdef ASSUMED_SHAPE
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#  ifdef WEC_VF
      real(r8), intent(in) :: zetat(LBi:,LBj:)
#  endif
#  ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:,LBj:)
#  endif
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
# else
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
#  ifdef WEC_VF
      real(r8), intent(in) :: zetat(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
#  endif
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif
      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: eps = 1.0E-8_r8

      real(r8) :: Ampl, Hdd, cff, cff1, cff2, cff3, cffL, cffR
      real(r8) :: deltaL, deltaR, dh, dP, limtr, rr
#ifdef ATM_PRESS
      real(r8) :: OneAtm, fac
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: P
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: d

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: FX

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: aL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: aR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: d1
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: r1

# include "set_bounds.h"
!
!---------------------------------------------------------------------
!  Finite-volume pressure gradient force algorithm.
!---------------------------------------------------------------------
!
# ifdef ATM_PRESS
      OneAtm=1013.25_r8                  ! 1 atm = 1013.25 mb
      fac=100.0_r8/g
# endif
      DO j=JstrV-1,Jend
        DO k=N(ng)-1,1,-1
          DO i=IstrU-1,Iend
            FC(i,k)=1.0_r8/(Hz(i,j,k+1)+Hz(i,j,k))
            r(i,j,k)=FC(i,k)*(rho(i,j,k+1)*Hz(i,j,k  )+                 &
     &                        rho(i,j,k  )*Hz(i,j,k+1))
            d(i,j,k)=FC(i,k)*(rho(i,j,k+1)-rho(i,j,k))
          END DO
        END DO
!
!  Parabolic WENO reconstruction of density field. Compute left and
!  right side limits aL and aR for the density assuming monotonized
!  quartic polynomial distributions within each grid box.  Also
!  compute dL and dR which are used as a measure of quadratic
!  variation during subsquent WENO reconciliation of side limits.
!
        DO k=2,N(ng)-1
          DO i=IstrU-1,Iend
            deltaR=Hz(i,j,k)*d(i,j,k  )
            deltaL=Hz(i,j,k)*d(i,j,k-1)
            IF ((deltaR*deltaL).lt.0.0_r8) THEN
              deltaR=0.0_r8
              deltaL=0.0_r8
            END IF
            cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
            cffR=cff*d(i,j,k  )
            cffL=cff*d(i,j,k-1)
            IF (ABS(deltaR).gt.ABS(cffL)) deltaR=cffL
            IF (ABS(deltaL).gt.ABS(cffR)) deltaL=cffR
            cff=(deltaR-deltaL)/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
            deltaR=deltaR-cff*Hz(i,j,k+1)
            deltaL=deltaL+cff*Hz(i,j,k-1)
            aR(i,k)=rho(i,j,k)+deltaR
            aL(i,k)=rho(i,j,k)-deltaL
            dR(i,k)=(2.0_r8*deltaR-deltaL)**2
            dL(i,k)=(2.0_r8*deltaL-deltaR)**2
          END DO
        END DO
!
        DO i=IstrU-1,Iend
          aL(i,N(ng))=aR(i,N(ng)-1)
          aR(i,N(ng))=2.0_r8*rho(i,j,N(ng))-aL(i,N(ng))
          dR(i,N(ng))=(2.0_r8*aR(i,N(ng))+aL(i,N(ng))-                  &
     &                 3.0_r8*rho(i,j,N(ng)))**2
          dL(i,N(ng))=(3.0_r8*rho(i,j,N(ng))-                           &
     &                 2.0_r8*aL(i,N(ng))-aR(i,N(ng)))**2
          aR(i,1)=aL(i,2)
          aL(i,1)=2.0_r8*rho(i,j,1)-aR(i,1)
          dR(i,1)=(2.0_r8*aR(i,1)+aL(i,1)-3.0_r8*rho(i,j,1))**2
          dL(i,1)=(3.0_r8*rho(i,j,1)-2.0_r8*aL(i,1)-aR(i,1))**2
        END DO
!
        DO k=1,N(ng)-1
          DO i=IstrU-1,Iend
            deltaL=MAX(dL(i,k  ),eps)
            deltaR=MAX(dR(i,k+1),eps)
            r1(i,k)=(deltaR*aR(i,k)+deltaL*aL(i,k+1))/(deltaR+deltaL)
          END DO
        END DO
!
        DO i=IstrU-1,Iend
# ifdef NEUMANN
          r1(i,N(ng))=1.5_r8*rho(i,j,N(ng))-0.5_r8*r1(i,N(ng)-1)
          r1(i,0)=1.5_r8*rho(i,j,1)-0.5_r8*r1(i,1)
# else
          r1(i,N(ng))=2.0_r8*rho(i,j,N(ng))-r1(i,N(ng)-1)
          r1(i,0)=2.0_r8*rho(i,j,1)-r1(i,1)
# endif
        END DO
!
!  Power-law reconciliation step.  It starts with the computation of
!  side limits dR and dL of the first derivative assuming parabolic
!  distributions within each grid box.  In this version of the code,
!  before doing so (see "else" branch of 3-way switch below), in the
!  situation when interfacial deviations deltaR and deltaL differ by
!  more than a factor of two (hence monotonic parabolic fit becomes
!  impossible), the parabolic assumption is switched to power-law
!  function,  such that its derivative is zero at one end and,
!  consequently, larger than that of (would be) limited parabolic
!  on the other end.  The basic parabolic version of the code is
!  commented out, but left here for reference.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
!!          cff=2.0_r8/Hz(i,j,k)
!!          dR(i,k)=cff*(2.0_r8*r1(i,k)+r1(i,k-1)-3.0_r8*rho(i,j,k))
!!          dL(i,k)=cff*(3.0_r8*rho(i,j,k)-2.0_r8*r1(i,k-1)-r1(i,k))
!!          cff=r(i,j,k)-r(i,j,k-1)
!!          if (cff*dR(i,k).lt.0.0_r8) dR(i,k)=0.0_r8
!!          if (cff*dL(i,k).lt.0.0_r8) dL(i,k)=0.0_r8
            deltaR=r1(i,k)-rho(i,j,k)
            deltaL=rho(i,j,k)-r1(i,k-1)
            cff=deltaR*deltaL
            IF (cff.gt.eps) THEN
              cff=(deltaR+deltaL)/cff
            ELSE
              cff=0.0_r8
            END IF
            cffL=cff*deltaL
            cffR=cff*deltaR
            IF (cffL.gt.3.0_r8) THEN
              cffL=cffL*deltaL
              cffR=0.0_r8
            ELSE IF (cffR.gt.3.0_r8) THEN
              cffL=0.0_r8
              cffR=cffR*deltaR
            ELSE
              cffL=4.0_r8*deltaL-2.0_r8*deltaR
              cffR=4.0_r8*deltaR-2.0_r8*deltaL
            END IF
            cff=1.0_r8/Hz(i,j,k)
            dR(i,k)=cff*cffR
            dL(i,k)=cff*cffL
          END DO
        END DO
!
!  Compute final value of derivative at each interface by reconciling
!  two side limits dR(k) and dL(k+1) coming from adjacent grid boxes.
!  The difference between these two also causes change of interfacial
!  value r(k) by Ampl. The commented code (left here for reference)
!  computes the exact value of Ampl assuming power law reconciliation
!  and solving associated quadratic equation. The code segment below
!  corresponds to Pade fit to exact solution, which avoids computation
!  of SQRT for the sake of computational efficiency.
!
        DO k=N(ng)-1,1,-1
          DO i=IstrU-1,Iend
            d(i,j,k)=FC(i,k)*(Hz(i,j,k+1)*dL(i,k+1)+Hz(i,j,k)*dR(i,k))
            cffR=8.0_r8*(dR(i,k  )+2.0_r8*dL(i,k  ))
            cffL=8.0_r8*(dL(i,k+1)+2.0_r8*dR(i,k+1))
            IF (ABS(d(i,j,k)).gt.ABS(cffR)) d(i,j,k)=cffR
            IF (ABS(d(i,j,k)).gt.ABS(cffL)) d(i,j,k)=cffL
            IF ((dL(i,k+1)-dR(i,k))*                                    &
     &          (rho(i,j,k+1)-rho(i,j,k)).gt.0.0_r8) THEN
              Hdd=Hz(i,j,k)*(d(i,j,k)-dR(i,k))
              rr=rho(i,j,k)-r1(i,k-1)
            ELSE
              Hdd=Hz(i,j,k+1)*(dL(i,k+1)-d(i,j,k))
              rr=r1(i,k+1)-rho(i,j,k+1)
            END IF
            rr=abs(rr)
!!          Ampl=0.4_r8*Hdd*rr
!!          Hdd=ABS(Hdd)
!!          cff=rr*(rr+0.16_r8*Hdd)
!!          if (cff.gt.eps) Ampl=Ampl/(rr+sqrt(cff))
            Ampl=0.2_r8*Hdd*rr
            Hdd=ABS(Hdd)
            cff=rr*rr+0.0763636363636363636_r8*Hdd*                     &
     &                (rr+0.004329004329004329_r8*Hdd)
            IF (cff.gt.eps) THEN
              Ampl=Ampl*(rr+0.0363636363636363636_r8*Hdd)/cff
            ELSE
              Ampl=0.0_r8
            END IF
            r(i,j,k)=r1(i,k)+Ampl
          END DO
        END DO
        DO i=IstrU-1,Iend
# ifdef NEUMANN
          r(i,j,0)=1.5_r8*rho(i,j,1)-0.5_r8*r(i,j,1)
          r(i,j,N(ng))=1.5_r8*rho(i,j,N(ng))-0.5_r8*r(i,j,N(ng)-1)
          d(i,j,0)=0.0_r8
          d(i,j,N(ng))=0.0_r8
# else
          r(i,j,0)=2.0_r8*rho(i,j,1)-r(i,j,1)
          r(i,j,N(ng))=2.0_r8*rho(i,j,N(ng))-r(i,j,N(ng)-1)
          d(i,j,0)=d(i,j,1)
          d(i,j,N(ng))=d(i,j,N(ng)-1)
# endif
        END DO
!
!  Compute pressure (P) and lateral pressure force (FX). Initialize
!  pressure at the free-surface as zero
!
        DO i=IstrU-1,Iend
          P(i,j,N(ng))=0.0_r8
#ifdef WEC_VF
          P(i,j,N(ng))=P(i,j,N(ng))+zetat(i,j)
#endif
#ifdef ATM_PRESS
          P(i,j,N(ng))=P(i,j,N(ng))+fac*(Pair(i,j)-OneAtm)
#endif
        END DO
        cff3=1.0_r8/12.0_r8
        DO k=N(ng),1,-1
          DO i=IstrU-1,Iend
            P(i,j,k-1)=P(i,j,k)+Hz(i,j,k)*rho(i,j,k)
            FX(i,j,k)=0.5_r8*Hz(i,j,k)*                                 &
     &                (P(i,j,k)+P(i,j,k-1)+                             &
     &                 0.2_r8*Hz(i,j,k)*                                &
     &                 (r(i,j,k)-r(i,j,k-1)-                            &
     &                  cff3*Hz(i,j,k)*(d(i,j,k)+d(i,j,k-1))))
          END DO
        END DO
!
!  Compute net pressure gradient forces in the XI- and ETA-directions.
!
        IF (j.ge.Jstr) THEN
          DO i=IstrU,Iend
            FC(i,N(ng))=0.0_r8
          END DO
          cff=0.5_r8*g
          cff1=g/rho0
          cff2=1.0_r8/6.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=N(ng),1,-1
            DO i=IstrU,Iend
              dh=z_w(i,j,k-1)-z_w(i-1,j,k-1)
              dP=P(i-1,j,k-1)-P(i,j,k-1)
              rr=0.5_r8*dh*(r(i,j,k-1)+r(i-1,j,k-1)-                    &
     &                      cff2*dh*(d(i,j,k-1)-d(i-1,j,k-1)))
              limtr=2.0_r8*dP*rr
              rr=rr*rr+dP*dP
              IF (limtr.gt.eps*rr) THEN
                limtr=limtr/rr
              ELSE
                limtr=0.0_r8
              END IF
              FC(i,k-1)=0.5_r8*dh*                                      &
     &                  (P(i,j,k-1)+P(i-1,j,k-1)+                       &
     &                   limtr*0.2_r8*dh*                               &
     &                   (r(i,j,k-1)-r(i-1,j,k-1)-                      &
     &                    cff3*dh*(d(i,j,k-1)+d(i-1,j,k-1))))
              ru(i,j,k,nrhs)=(cff*(Hz(i-1,j,k)+Hz(i,j,k))*              &
     &                            (z_w(i-1,j,N(ng))-z_w(i,j,N(ng)))+    &
     &                        cff1*(FX(i-1,j,k)-FX(i,j,k)+              &
     &                              FC(i,k)-FC(i,k-1)))*on_u(i,j)
# ifdef WET_DRY
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)*umask_wet(i,j)
# endif
# ifdef DIAGNOSTICS_UV
              DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
# endif
            END DO
          END DO
        END IF
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
          END DO
          cff=0.5_r8*g
          cff1=g/rho0
          cff2=1.0_r8/6.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=N(ng),1,-1
            DO i=Istr,Iend
              dh=z_w(i,j,k-1)-z_w(i,j-1,k-1)
              dP=P(i,j-1,k-1)-P(i,j,k-1)
              rr=0.5_r8*dh*(r(i,j,k-1)+r(i,j-1,k-1)-                    &
     &                      cff2*dh*(d(i,j,k-1)-d(i,j-1,k-1)))
              limtr=2.0_r8*dP*rr
              rr=rr*rr+dP*dP
              IF (limtr.gt.eps*rr) THEN
                limtr=limtr/rr
              ELSE
                limtr=0.0_r8
              END IF
              FC(i,k-1)=0.5_r8*dh*                                      &
     &                  (P(i,j,k-1)+P(i,j-1,k-1)+                       &
     &                   limtr*0.2_r8*dh*                               &
     &                   (r(i,j,k-1)-r(i,j-1,k-1)-                      &
     &                    cff3*dh*(d(i,j,k-1)+d(i,j-1,k-1))))
              rv(i,j,k,nrhs)=(cff*(Hz(i,j-1,k)+Hz(i,j,k))*              &
     &                            (z_w(i,j-1,N(ng))-z_w(i,j,N(ng)))+    &
     &                        cff1*(FX(i,j-1,k)-FX(i,j,k)+              &
     &                              FC(i,k)-FC(i,k-1)))*om_v(i,j)
# ifdef WET_DRY
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)*vmask_wet(i,j)
# endif
# ifdef DIAGNOSTICS_UV
              DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
# endif
           END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE prsgrd_tile
