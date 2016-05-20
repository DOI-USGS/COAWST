#undef NEUMANN
      SUBROUTINE prsgrd (ng, tile)
!
!svn $Id: prsgrd42.h 732 2008-09-07 01:55:51Z jcwarner $
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
!  uses a PPM-style limitting algorithm.                               !
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
# ifdef MASKING
     &                  GRID(ng) % umask,                               &
     &                  GRID(ng) % vmask,                               &
# endif
# ifdef WET_DRY
     &                  GRID(ng)%umask_wet,                             &
     &                  GRID(ng)%vmask_wet,                             &
# endif
     &                  GRID(ng) % Hz,                                  &
     &                  GRID(ng) % om_v,                                &
     &                  GRID(ng) % on_u,                                &
     &                  GRID(ng) % z_w,                                 &
     &                  OCEAN(ng) % rho,                                &
#ifdef WEC_VF
     &                  OCEAN(ng) % zetat,                              &
#endif
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
# ifdef MASKING
     &                        umask, vmask,                             &
# endif
# ifdef WET_DRY
     &                        umask_wet, vmask_wet,                     &
# endif
     &                        Hz, om_v, on_u, z_w,                      &
     &                        rho,                                      &
#ifdef WEC_VF
     &                        zetat,                                    &
#endif
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
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
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
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
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

      real(r8) :: cff, cff1, cff2, cffL, cffR
      real(r8) :: deltaL, deltaR, dh, dP, rr
#ifdef ATM_PRESS
      real(r8) :: OneAtm, fac
#endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: P
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: r

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: aL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: aR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR

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
      cff2=1.0_r8/6.0_r8
      DO j=JstrV-2,Jend+1
        DO k=N(ng)-1,1,-1
          DO i=IstrU-2,Iend+1
            FC(i,k)=(rho(i,j,k+1)-rho(i,j,k))/(Hz(i,j,k+1)+Hz(i,j,k))
          END DO
        END DO
!
!  Parabolic WENO reconstruction of density field. Compute left and
!  right side limits aL and aR for the density assuming monotonized
!  parabolic distributions within each grid box.  Also compute dL and
!  dR which are used as a measure of quadratic variation during
!  subsquent WENO reconciliation of side limits.
!
        DO k=2,N(ng)-1
          DO i=IstrU-2,Iend+1
            deltaR=Hz(i,j,k)*FC(i,k)
            deltaL=Hz(i,j,k)*FC(i,k-1)
            IF ((deltaR*deltaL).lt.0.0_r8) THEN
              deltaR=0.0_r8
              deltaL=0.0_r8
            END IF
            cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
            cffR=cff*FC(i,k)
            cffL=cff*FC(i,k-1)
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
        DO i=IstrU-2,Iend+1
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
          DO i=IstrU-2,Iend+1
             deltaL=MAX(dL(i,k  ),eps)
             deltaR=MAX(dR(i,k+1),eps)
             r(i,j,k)=(deltaR*aR(i,k)+deltaL*aL(i,k+1))/                &
     &                (deltaR+deltaL)
          END DO
        END DO
!
        DO i=IstrU-2,Iend+1
# ifdef NEUMANN
          r(i,j,N(ng))=1.5_r8*rho(i,j,N(ng))-0.5_r8*r(i,j,N(ng)-1)
          r(i,j,0)=1.5_r8*rho(i,j,1)-0.5_r8*r(i,j,1  )
# else
          r(i,j,N(ng))=2.0_r8*rho(i,j,N(ng))-r(i,j,N(ng)-1)
          r(i,j,0)=2.0_r8*rho(i,j,1)-r(i,j,1  )
# endif
        END DO
!
!  Compute pressure (P) and lateral pressure force (FX). Initialize
!  pressure at the free-surface as zero
!
        DO i=IstrU-2,Iend+1
          P(i,j,N(ng))=0.0_r8
#ifdef WEC_VF
          P(i,j,N(ng))=P(i,j,N(ng))+zetat(i,j)
#endif
#ifdef ATM_PRESS
          P(i,j,N(ng))=P(i,j,N(ng))+fac*(Pair(i,j)-OneAtm)
#endif
        END DO
        DO k=N(ng),1,-1
          DO i=IstrU-2,Iend+1
            P(i,j,k-1)=P(i,j,k)+Hz(i,j,k)*rho(i,j,k)
            deltaR=r(i,j,k)-rho(i,j,k)
            deltaL=rho(i,j,k)-r(i,j,k-1)
            IF ((deltaR*deltaL).lt.0.0_r8) THEN
              rr=0.0_r8
            ELSE IF (ABS(deltaR).gt.(2.0_r8*ABS(deltaL))) THEN
              rr=3.0_r8*deltaL
            ELSE IF (ABS(deltaL).gt.(2.0_r8*ABS(deltaR))) THEN
              rr=3.0_r8*deltaR
            ELSE
              rr=deltaR+deltaL
            END IF
            FX(i,j,k)=0.5_r8*Hz(i,j,k)*                                 &
     &                (P(i,j,k)+P(i,j,k-1)+cff2*rr*Hz(i,j,k))
          END DO
        END DO
!
!  Compute net pressure gradient forces in the XI-directions.
!  Set pressure at free-surface as zero.
!
        IF ((j.ge.Jstr).and.(j.le.Jend)) THEN
          DO i=IstrU-1,Iend+1
            FC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng),1,-1
            DO i=IstrU-1,Iend+1
              dP=P(i-1,j,k-1)-P(i,j,k-1)
              dh=z_w(i,j,k-1)-z_w(i-1,j,k-1)
              deltaR=dh*r(i,j,k-1)-dP
              deltaL=dP-dh*r(i-1,j,k-1)
              IF ((deltaR*deltaL).lt.0.0_r8) THEN
                rr=0.0_r8
              ELSE IF (ABS(deltaR).gt.(2.0_r8*ABS(deltaL))) THEN
                rr=3.0_r8*deltaL
              ELSE IF (ABS(deltaL).gt.(2.0_r8*ABS(deltaR))) THEN
                rr=3.0_r8*deltaR
              ELSE
                rr=deltaR+deltaL
              END IF
              FC(i,k-1)=0.5_r8*dh*(P(i,j,k-1)+P(i-1,j,k-1)+cff2*rr)
              ru(i,j,k,nrhs)=2.0_r8*(FX(i-1,j,k)-FX(i,j,k)+             &
     &                               FC(i,k)-FC(i,k-1))/                &
     &                       (Hz(i-1,j,k)+Hz(i,j,k))
# ifdef MASKING
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)*umask(i,j)
# endif
# ifdef WET_DRY
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)*umask_wet(i,j)
# endif
            END DO
          END DO
        END IF
!
!  Compute net pressure gradient forces in the ETA-directions.
!  Set pressure at free-surface as zero.
!
        IF (j.ge.JstrV-1) THEN
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng),1,-1
            DO i=Istr,Iend
              dP=P(i,j-1,k-1)-P(i,j,k-1)
              dh=z_w(i,j,k-1)-z_w(i,j-1,k-1)
              deltaR=dh*r(i,j,k-1)-dP
              deltaL=dP-dh*r(i,j-1,k-1)
              IF ((deltaR*deltaL).lt.0.0_r8) THEN
                rr=0.0_r8
              ELSE IF (ABS(deltaR).gt.(2.0_r8*ABS(deltaL))) THEN
                rr=3.0_r8*deltaL
              ELSE IF (ABS(deltaL).gt.(2.0_r8*ABS(deltaR))) THEN
                rr=3.0_r8*deltaR
              ELSE
                rr=deltaR+deltaL
              END IF
              FC(i,k-1)=0.5_r8*dh*(P(i,j,k-1)+P(i,j-1,k-1)+cff2*rr)
              rv(i,j,k,nrhs)=2.0_r8*(FX(i,j-1,k)-FX(i,j,k)+             &
     &                               FC(i,k)-FC(i,k-1))/                &
     &                       (Hz(i,j-1,k)+Hz(i,j,k))
# ifdef MASKING
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)*vmask(i,j)
# endif
# ifdef WET_DRY
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)*vmask_wet(i,j)
# endif
            END DO
          END DO
        END IF
      END DO
!
      rr=g/(24.0_r8*rho0)
      cff=0.5_r8*g
      cff1=0.5_r8*g/rho0
      DO j=Jstr,Jend
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            dh=rr*(z_w(i,j,k)-z_w(i-1,j,k))
            FC(i,k)=MAX(dh,0.0_r8)*                                     &
     &                 (ru(i,j,k+1,nrhs)+ru(i+1,j,k  ,nrhs)-            &
     &                  ru(i,j,k  ,nrhs)-ru(i-1,j,k+1,nrhs))+           &
     &              MIN(dh,0.0_r8)*                                     &
     &                 (ru(i,j,k  ,nrhs)+ru(i+1,j,k+1,nrhs)-            &
     &                  ru(i,j,k+1,nrhs)-ru(i-1,j,k  ,nrhs))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          dh=rr*(z_w(i,j,0)-z_w(i-1,j,0))
          FC(i,0)=MAX(dh,0.0_r8)*                                       &
     &               (ru(i  ,j,1,nrhs)-ru(i-1,j,1,nrhs))+               &
     &            MIN(dh,0.0_r8)*                                       &
     &               (ru(i+1,j,1,nrhs)-ru(i  ,j,1,nrhs))
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            ru(i,j,k,nrhs)=(cff*(z_w(i-1,j,N(ng))-z_w(i,j,N(ng)))+      &
     &                      cff1*ru(i,j,k,nrhs))*                       &
     &                     (Hz(i-1,j,k)+Hz(i,j,k))*on_u(i,j)+           &
     &                     (FC(i,k)-FC(i,k-1))*on_u(i,j)
# ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
# endif
          END DO
        END DO
      END DO
!
      DO j=JstrV,Jend
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            dh=rr*(z_w(i,j,k)-z_w(i,j-1,k))
            FX(i,j,k)=MAX(dh,0.0_r8)*                                   &
     &                   (rv(i,j,k+1,nrhs)+rv(i+1,j  ,k  ,nrhs)-        &
     &                    rv(i,j,k  ,nrhs)-rv(i  ,j-1,k+1,nrhs))+       &
     &                MIN(dh,0.0_r8)*                                   &
     &                   (rv(i,j,k  ,nrhs)+rv(i+1,j  ,k+1,nrhs)-        &
     &                    rv(i,j,k+1,nrhs)-rv(i  ,j-1,k  ,nrhs))
          END DO
        END DO
        DO i=Istr,Iend
          FX(i,j,N(ng))=0.0_r8
          dh=rr*(z_w(i,j,0)-z_w(i,j-1,0))
          FX(i,j,0)=MAX(dh,0.0_r8)*                                     &
     &                 (rv(i  ,j,1,nrhs)-rv(i,j-1,1,nrhs))+             &
     &              MIN(dh,0.0_r8)*                                     &
     &                 (rv(i+1,j,1,nrhs)-rv(i,j  ,1,nrhs))
        END DO
      END DO
      DO j=JstrV,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            rv(i,j,k,nrhs)=(cff*(z_w(i,j-1,N(ng))-z_w(i,j,N(ng)))+      &
     &                      cff1*rv(i,j,k,nrhs))*                       &
     &                     (Hz(i,j-1,k)+Hz(i,j,k))*om_v(i,j)+           &
     &                     (FX(i,j,k)-FX(i,j,k-1))*om_v(i,j)
# ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
# endif
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE prsgrd_tile
