       SUBROUTINE ad_prsgrd (ng, tile)
!
!svn $Id: ad_prsgrd40.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This subroutine  evaluates the  adjoint  baroclinic,  hydrostatic   !
!  pressure gradient term using the finite-volume pressure  Jacobian   !
!  scheme of Lin (1997).                                               !
!                                                                      !
!  The pressure gradient terms (m4/s2) are loaded into right-hand-     !
!  side arrays "tl_ru" and "tl_rv".                                    !
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
      CALL wclock_on (ng, iADM, 23)
#endif
      CALL ad_prsgrd_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
     &                     GRID(ng) % om_v,                             &
     &                     GRID(ng) % on_u,                             &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % ad_Hz,                            &
     &                     GRID(ng) % z_w,                              &
     &                     GRID(ng) % ad_z_w,                           &
     &                     OCEAN(ng) % rho,                             &
     &                     OCEAN(ng) % ad_rho,                          &
#ifdef ATM_PRESS
     &                     FORCES(ng) % Pair,                           &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                     DIAGS(ng) % DiaRU,                           &
!!   &                     DIAGS(ng) % DiaRV,                           &
#endif
     &                     OCEAN(ng) % ad_ru,                           &
     &                     OCEAN(ng) % ad_rv)
#ifdef PROFILE
      CALL wclock_off (ng, iADM, 23)
#endif
      RETURN
      END SUBROUTINE ad_prsgrd
!
!***********************************************************************
      SUBROUTINE ad_prsgrd_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
     &                           om_v, on_u,                            &
     &                           Hz, ad_Hz,                             &
     &                           z_w, ad_z_w,                           &
     &                           rho, ad_rho,                           &
#ifdef ATM_PRESS
     &                           Pair,                                  &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                           DiaRU, DiaRV,                          &
#endif
     &                           ad_ru, ad_rv)
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
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:,LBj:)
# endif
# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
!!    real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: ad_rv(LBi:,LBj:,0:,:)
#else
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
# ifdef ATM_PRESS
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
# endif
# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
!!    real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
# endif
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: ad_rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: cff, cff1, dh
      real(r8) :: ad_dh
#ifdef ATM_PRESS
      real(r8) :: OneAtm, fac
#endif
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: ad_FX

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: P
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: ad_P

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_dh=0.0_r8
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          ad_FX=(i,j)=0.0_r8
        END DO
        DO k=0,N(ng)
          DO i=IminS,ImaxS
            ad_P(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
      DO k=0,N(ng)
        DO i=IminS,ImaxS
          ad_FC=(i,k)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Finite Volume pressure gradient algorithm (Lin, 1997).
!-----------------------------------------------------------------------
!
!  Compute BASIC STATE dynamic pressure, P.
!
#ifdef ATM_PRESS
      OneAtm=1013.25_r8                  ! 1 atm = 1013.25 mb
      fac=100.0_r8/g
#endif
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
#ifdef ATM_PRESS
          P(i,j,N(ng))=fac*(Pair(i,j)-OneAtm)
#else
          P(i,j,N(ng))=0.0_r8
#endif
        END DO
        DO k=N(ng),1,-1
          DO i=IstrU-1,Iend
            P(i,j,k-1)=P(i,j,k)+Hz(i,j,k)*rho(i,j,k)
          END DO
        END DO
      END DO
!
!  Calculate adjoint of pressure gradient in the ETA-direction (m4/s2).
!
      J_LOOP: DO j=Jend,JstrV-1,-1
        IF (j.ge.JstrV) THEN
          cff=0.5_r8*g
          cff1=g/rho0
          DO k=1,N(ng)
            DO i=Istr,Iend
              dh=z_w(i,j,k-1)-z_w(i,j-1,k-1)
#ifdef DIAGNOSTICS_UV
!!            DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
#endif
!>            tl_rv(i,j,k,nrhs)=(cff*((tl_Hz(i,j-1,k)+                  &
!>   &                                 tl_Hz(i,j  ,k))*                 &
!>   &                                (z_w(i,j-1,N(ng))-                &
!>   &                                 z_w(i,j  ,N(ng)))+               &
!>   &                                (Hz(i,j-1,k)+                     &
!>   &                                 Hz(i,j  ,k))*                    &
!>   &                                (tl_z_w(i,j-1,N(ng))-             &
!>   &                                 tl_z_w(i,j  ,N(ng))))+           &
!>   &                           cff1*(tl_FX(i,j-1,k)-                  &
!>   &                                 tl_FX(i,j  ,k)+                  &
!>   &                                 tl_FC(i,k  )-                    &
!>   &                                 tl_FC(i,k-1)))*om_v(i,j)
!>
              adfac=om_v(i,j)*ad_rv(i,j,k,nrhs)
              adfac1=adfac*cff
              adfac2=adfac1*(z_w(i,j-1,N(ng))-                          &
     &                       z_w(i,j  ,N(ng)))
              adfac3=adfac1*(Hz(i,j-1,k)+                               &
     &                       Hz(i,j  ,k))
              adfac4=adfac*cff1
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac2
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac2
              ad_z_w(i,j-1,N(ng))=ad_z_w(i,j-1,N(ng))+adfac3
              ad_z_w(i,j  ,N(ng))=ad_z_w(i,j  ,N(ng))-adfac3
              ad_FX(i,j-1,k)=ad_FX(i,j-1,k)+adfac4
              ad_FX(i,j  ,k)=ad_FX(i,j  ,k)-adfac4
              ad_FC(i,k-1)=ad_FC(i,k-1)-adfac4
              ad_FC(i,k  )=ad_FC(i,k  )+adfac4
              ad_rv(i,j,k,nrhs)=0.0_r8
!>            tl_FC(i,k-1)=0.5_r8*                                      &
!>   &                     (tl_dh*(P(i,j,k-1)+P(i,j-1,k-1))+            &
!>   &                      dh*(tl_P(i,j,k-1)+tl_P(i,j-1,k-1))
!>
              adfac=0.5_r8*ad_FC(i,k-1)
              adfac1=adfac*dh
              ad_dh=ad_dh+(P(i,j,k-1)+P(i,j-1,k-1))*adfac
              ad_P(i,j-1,k-1)=ad_P(i,j-1,k-1)+adfac1
              ad_P(i,j  ,k-1)=ad_P(i,j  ,k-1)+adfac1
              ad_FC(i,k-1)=0.0_r8
!>            tl_dh=tl_z_w(i,j,k-1)-tl_z_w(i,j-1,k-1)
!>
              ad_z_w(i,j-1,k-1)=ad_z_w(i,j-1,k-1)-ad_dh
              ad_z_w(i,j  ,k-1)=ad_z_w(i,j  ,k-1)+ad_dh
              ad_dh=0.0_r8
            END DO
          END DO
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
          END DO
        END IF
!
!  Calculate the adjoint of pressure gradient in the XI-direction (m4/s2).
!
        IF (j.ge.Jstr) THEN
          cff=0.5_r8*g
          cff1=g/rho0
          DO k=1,N(ng)
            DO i=Iend,IstrU,-1
              dh=z_w(i,j,k-1)-z_w(i-1,j,k-1)
#ifdef DIAGNOSTICS_UV
!!            DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
#endif
!>            tl_ru(i,j,k,nrhs)=(cff*((tl_Hz(i-1,j,k)+                  &
!>   &                                 tl_Hz(i  ,j,k))*                 &
!>   &                                (z_w(i-1,j,N(ng))-                &
!>   &                                 z_w(i  ,j,N(ng)))+               &
!>   &                                (Hz(i-1,j,k)+                     &
!>   &                                 Hz(i  ,j,k))*                    &
!>   &                                (tl_z_w(i-1,j,N(ng))-             &
!>   &                                 tl_z_w(i  ,j,N(ng))))+           &
!>   &                           cff1*(tl_FX(i-1,j,k)-                  &
!>   &                                 tl_FX(i  ,j,k)+                  &
!>   &                                 tl_FC(i,k  )-                    &
!>   &                                 tl_FC(i,k-1)))*on_u(i,j)
!>
              adfac=on_u(i,j)*tl_ru(i,j,k,nrhs)
              adfac1=adfac*cff
              adfac2=adfac1*(z_w(i-1,j,N(ng))-                          &
     &                       z_w(i  ,j,N(ng)))
              adfac3=adfac1*(Hz(i-1,j,k)+                               &
     &                       Hz(i  ,j,k))
              adfac4=adfac*cff1
              ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac2
              ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac2
              ad_z_w(i-1,j,N(ng))=ad_z_w(i-1,j,N(ng))+adfac3
              ad_z_w(i  ,j,N(ng))=ad_z_w(i  ,j,N(ng))-adfac3
              ad_FX(i-1,j,k)=ad_FX(i-1,j,k)+adfac4
              ad_FX(i  ,j,k)=ad_FX(i  ,j,k)-adfac4
              ad_FC(i,k  )=ad_FC(i,k  )+adfac4
              ad_FC(i,k-1)=ad_FC(i,k-1)-adfac4
              ad_ru(i,j,k,nrhs)=0.0_r8
!>            tl_FC(i,k-1)=0.5_r8*                                      &
!>   &                     (tl_dh*(P(i,j,k-1)+P(i-1,j,k-1))+            &
!>   &                      dh*(tl_P(i,j,k-1)+tl_P(i-1,j,k-1)))
!>
              adfac=0.5_r8*tl_FC(i,k-1)
              adfac1=adfac*dh
              ad_dh=ad_dh+(P(i,j,k-1)+P(i-1,j,k-1))*adfac
              ad_P(i-1,j,k-1)=ad_P(i-1,j,k-1)+adfac1
              ad_P(i  ,j,k-1)=ad_P(i  ,j,k-1)+adfac1
              ad_FC(i,k-1)=0.0_r8
!>            tl_dh=tl_z_w(i,j,k-1)-tl_z_w(i-1,j,k-1)
!>
              ad_z_w(i-1,j,k-1)=ad_z_w(i-1,j,k-1)-ad_dh
              ad_z_w(i  ,j,k-1)=ad_z_w(i  ,j,k-1)+ad_dh
              ad_dh=0.0_r8
            END DO
          END DO
          DO i=IstrU,Iend
!>          tl_FC(i,N(ng))=0.0_r8
            ad_FC(i,N(ng))=0.0_r8
          END DO
        END IF
!
!  Compute adjoint pressure and its vertical integral.
!
        DO k=1,N(ng)
          DO i=IstrU-1,Iend
!>          tl_FX(i,j,k)=0.5_r8*                                        &
!>   &                   (tl_Hz(i,j,k)*(P(i,j,k)+P(i,j,k-1))+           &
!>   &                    Hz(i,j,k)*(tl_P(i,j,k)+tl_P(i,j,k-1)))
!>
            adfac=0.5_r8*ad_FX(i,j,k)
            adfac1=adfac*Hz(i,j,k)
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+(P(i,j,k)+P(i,j,k-1))*adfac
            ad_P(i,j,k-1)=ad_P(i,j,k-1)+adfac1
            ad_P(i,j,k  )=ad_P(i,j,k  )+adfac1
            ad_FX(i,j,k)=0.0_r8
!>          tl_P(i,j,k-1)=tl_P(i,j,k)+                                  &
!>   &                    tl_Hz(i,j,k)*rho(i,j,k)+                      &
!>   &                    Hz(i,j,k)*tl_rho(i,j,k)
!>
            ad_P(i,j,k)=ad_P(i,j,k)+ad_P(i,j,k-1)
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+rho(i,j,k)*ad_P(i,j,k-1)
            ad_rho(i,j,k)=ad_rho(i,j,k)+Hz(i,j,k)*ad_P(i,j,k-1)
          END DO
        END DO
        DO i=IstrU-1,Iend
!>        tl_P(i,j,N(ng))=0.0_r8
!>
          ad_P(i,j,N(ng))=0.0_r8
        END DO
      END DO J_LOOP
      RETURN
      END SUBROUTINE ad_prsgrd_tile
