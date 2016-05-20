       SUBROUTINE prsgrd (ng, tile)
!
!svn $Id: prsgrd40.h 732 2008-09-07 01:55:51Z jcwarner $
!***********************************************************************
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This subroutine evaluates the  nonlinear baroclinic,  hydrostatic   !
!  pressure gradient term using the finite-volume pressure  Jacobian   !
!  scheme of Lin (1997).                                               !
!                                                                      !
!  The pressure gradient terms (m4/s2) are loaded into right-hand-     !
!  side arrays "ru" and "rv".                                          !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Lin, Shian-Jiann, 1997:  A finite volume integration method       !
!      for computing pressure gradient force in general vertical       !
!      coordinates, Q. J. R. Meteorol. Soc., 123, 1749-1762.           !
!                                                                      !
!**********************************************************************=
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
#ifdef WET_DRY
     &                  GRID(ng)%umask_wet,                             &
     &                  GRID(ng)%vmask_wet,                             &
#endif
     &                  GRID(ng) % om_v,                                &
     &                  GRID(ng) % on_u,                                &
     &                  GRID(ng) % Hz,                                  &
     &                  GRID(ng) % z_w,                                 &
     &                  OCEAN(ng) % rho,                                &
#ifdef WEC_VF
     &                  OCEAN(ng) % zetat,                              &
#endif
#ifdef ATM_PRESS
     &                  FORCES(ng) % Pair,                              &
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
#ifdef WET_DRY
     &                        umask_wet, vmask_wet,                     &
#endif
     &                        om_v, on_u,                               &
     &                        Hz, z_w,                                  &
     &                        rho,                                      &
#ifdef WEC_VF
     &                        zetat,                                    &
#endif
#ifdef ATM_PRESS
     &                        Pair,                                     &
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
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
# ifdef WEC_VF
      real(r8), intent(in) :: zetat(LBi:,LBj:)
# endif
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:,LBj:)
# endif
# ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
#else
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
# ifdef WEC_VF
      real(r8), intent(in) :: zetat(LBi:UBi,LBj:UBj)
# endif
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
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

      real(r8) :: cff, cff1, dh
#ifdef WEC_VF
      real(r8) :: fac1
#endif
#ifdef ATM_PRESS
      real(r8) :: OneAtm, fac
#endif
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: FX

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: P

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Finite Volume pressure gradient algorithm (Lin, 1997).
!-----------------------------------------------------------------------
!
!  Compute pressure and its vertical integral.  Initialize pressure at
!  the free-surface as zero.
!
#ifdef WEC_VF
      fac1=rho0/g
#endif
#ifdef ATM_PRESS
      OneAtm=1013.25_r8                  ! 1 atm = 1013.25 mb
      fac=100.0_r8/g
#endif
      J_LOOP : DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          P(i,j,N(ng))=0.0_r8
#ifdef WEC_VF
          P(i,j,N(ng))=P(i,j,N(ng))+fac1*zetat(i,j)
#endif
#ifdef ATM_PRESS
          P(i,j,N(ng))=P(i,j,N(ng))+fac*(Pair(i,j)-OneAtm)
#endif
        END DO
        DO k=N(ng),1,-1
          DO i=IstrU-1,Iend
            P(i,j,k-1)=P(i,j,k)+                                        &
     &                 Hz(i,j,k)*rho(i,j,k)
            FX(i,j,k)=0.5_r8*Hz(i,j,k)*(P(i,j,k)+P(i,j,k-1))
          END DO
        END DO
!
!  Calculate pressure gradient in the XI-direction (m4/s2).
!
        IF (j.ge.Jstr) THEN
          DO i=IstrU,Iend
            FC(i,N(ng))=0.0_r8
          END DO
          cff=0.5_r8*g
          cff1=g/rho0
          DO k=N(ng),1,-1
            DO i=IstrU,Iend
              dh=z_w(i,j,k-1)-z_w(i-1,j,k-1)
              FC(i,k-1)=0.5_r8*dh*(P(i,j,k-1)+P(i-1,j,k-1))
              ru(i,j,k,nrhs)=(cff*(Hz(i-1,j,k)+                         &
     &                             Hz(i  ,j,k))*                        &
     &                            (z_w(i-1,j,N(ng))-                    &
     &                             z_w(i  ,j,N(ng)))+                   &
     &                        cff1*(FX(i-1,j,k)-                        &
     &                              FX(i  ,j,k)+                        &
     &                              FC(i,k  )-                          &
     &                              FC(i,k-1)))*on_u(i,j)
#ifdef WET_DRY
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)*umask_wet(i,j)
#endif
#ifdef DIAGNOSTICS_UV
              DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
#endif
            END DO
          END DO
        END IF
!
!  Calculate pressure gradient in the ETA-direction (m4/s2).
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
          END DO
          cff=0.5_r8*g
          cff1=g/rho0
          DO k=N(ng),1,-1
            DO i=Istr,Iend
              dh=z_w(i,j,k-1)-z_w(i,j-1,k-1)
              FC(i,k-1)=0.5_r8*dh*(P(i,j,k-1)+P(i,j-1,k-1))
              rv(i,j,k,nrhs)=(cff*(Hz(i,j-1,k)+                         &
     &                             Hz(i,j  ,k))*                        &
     &                            (z_w(i,j-1,N(ng))-                    &
     &                             z_w(i,j  ,N(ng)))+                   &
     &                        cff1*(FX(i,j-1,k)-                        &
     &                              FX(i,j  ,k)+                        &
     &                              FC(i,k  )-                          &
     &                              FC(i,k-1)))*om_v(i,j)
#ifdef WET_DRY
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)*vmask_wet(i,j)
#endif
#ifdef DIAGNOSTICS_UV
              DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
#endif
            END DO
          END DO
        END IF
      END DO J_LOOP
      RETURN
      END SUBROUTINE prsgrd_tile
