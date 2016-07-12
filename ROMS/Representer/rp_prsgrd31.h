      SUBROUTINE rp_prsgrd (ng, tile)
!
!svn $Id: rp_prsgrd31.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine evalutes the representers tangent linear baroclinic    !
!  hydrostatic  pressure  gradient  term  using  the  standard  and    !
!  weighted Jacobians scheme  of Song and Wright (1997).  Note that    !
!  horizontal  gradients  are  computed  before   of  the  vertical    !
!  integration.                                                        !
!                                                                      !
!  The pressure gradient terms (m4/s2) are loaded into  right-hand-    !
!  side arrays "tl_ru" and "tl_rv".                                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Song, Y.T., 1998:  A general pressure gradient formulation for    !
!      numerical ocean models. Part I: Scheme design and diagnostic    !
!      analysis, Monthly Weather Rev., 126, 3213-3230.                 !
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
      CALL wclock_on (ng, iRPM, 23)
#endif
      CALL rp_prsgrd_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
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
      CALL wclock_off (ng, iRPM, 23)
#endif
      RETURN
      END SUBROUTINE rp_prsgrd
!
!***********************************************************************
      SUBROUTINE rp_prsgrd_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
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

      real(r8) :: fac, fac1, fac2, fac3
      real(r8) :: cff1, cff2, cff3, cff4
      real(r8) :: tl_cff1, tl_cff2, tl_cff3, tl_cff4
#ifdef WJ_GRADP
      real(r8) :: gamma, tl_gamma
#endif

      real(r8), dimension(IminS:ImaxS) :: phie
      real(r8), dimension(IminS:ImaxS) :: phix

      real(r8), dimension(IminS:ImaxS) :: tl_phie
      real(r8), dimension(IminS:ImaxS) :: tl_phix

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Calculate pressure gradient in the XI-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute surface baroclinic pressure gradient.
!
#ifdef ATM_PRESS
      fac=100.0_r8/rho0
#endif
      fac1=0.5_r8*g/rho0
      fac2=1000.0_r8*g/rho0
      fac3=0.25_r8*g/rho0

      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=z_w(i  ,j,N(ng))-z_r(i  ,j,N(ng))+                       &
     &         z_w(i-1,j,N(ng))-z_r(i-1,j,N(ng))
          tl_cff1=tl_z_w(i  ,j,N(ng))-tl_z_r(i  ,j,N(ng))+              &
     &            tl_z_w(i-1,j,N(ng))-tl_z_r(i-1,j,N(ng))
          phix(i)=fac1*(rho(i,j,N(ng))-rho(i-1,j,N(ng)))*cff1
          tl_phix(i)=fac1*                                              &
     &               ((tl_rho(i,j,N(ng))-tl_rho(i-1,j,N(ng)))*cff1+     &
     &                (rho(i,j,N(ng))-rho(i-1,j,N(ng)))*tl_cff1)-       &
#ifdef TL_IOMS
     &               phix(i)
#endif
#ifdef ATM_PRESS
          phix(i)=phix(i)+fac*(Pair(i,j)-Pair(i-1,j))
#endif
#ifdef RHO_SURF
          phix(i)=phix(i)+                                              &
     &            (fac2+fac1*(rho(i,j,N(ng))+rho(i-1,j,N(ng))))*        &
     &            (z_w(i,j,N(ng))-z_w(i-1,j,N(ng)))
          tl_phix(i)=tl_phix(i)+                                        &
     &               (fac1*(tl_rho(i,j,N(ng))+tl_rho(i-1,j,N(ng))))*    &
     &               (z_w(i,j,N(ng))-z_w(i-1,j,N(ng)))+                 &
     &               (fac2+fac1*(rho(i,j,N(ng))+rho(i-1,j,N(ng))))*     &
     &               (tl_z_w(i,j,N(ng))-tl_z_w(i-1,j,N(ng)))-           &
# ifdef TL_IOMS
     &               fac1*(rho(i,j,N(ng))+rho(i-1,j,N(ng)))*            &
     &                    (z_w(i,j,N(ng))-z_w(i-1,j,N(ng)))
# endif
#endif
!>        ru(i,j,N(ng),nrhs)=-0.5_r8*(Hz(i,j,N(ng))+Hz(i-1,j,N(ng)))*   &
!>   &                       phix(i)*on_u(i,j)
!>
          tl_ru(i,j,N(ng),nrhs)=-0.5_r8*on_u(i,j)*                      &
     &                          ((tl_Hz(i  ,j,N(ng))+                   &
     &                            tl_Hz(i-1,j,N(ng)))*phix(i)+          &
     &                           (Hz(i  ,j,N(ng))+                      &
     &                            Hz(i-1,j,N(ng)))*tl_phix(i))+         &
#ifdef TL_IOMS
     &                          0.5_r8*on_u(i,j)*                       &
     &                          (Hz(i  ,j,N(ng))+                       &
     &                           Hz(i-1,j,N(ng)))*phix(i)
#endif
#ifdef DIAGNOSTICS_UV
!!        DiaRU(i,j,N(ng),nrhs,M3pgrd)=ru(i,j,N(ng),nrhs)
#endif
        END DO
!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate.
!
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
#ifdef WJ_GRADP
            cff1=1.0_r8/((z_r(i  ,j,k+1)-z_r(i  ,j,k))*                 &
     &                   (z_r(i-1,j,k+1)-z_r(i-1,j,k)))
            tl_cff1=-cff1*cff1*((tl_z_r(i  ,j,k+1)-tl_z_r(i  ,j,k))*    &
     &                          (z_r(i-1,j,k+1)-z_r(i-1,j,k))+          &
     &                          (z_r(i  ,j,k+1)-z_r(i  ,j,k))*          &
     &                          (tl_z_r(i-1,j,k+1)-tl_z_r(i-1,j,k)))+   &
# ifdef TL_IOMS
     &              3.0_r8*cff1
# endif
            cff2=z_r(i  ,j,k  )-z_r(i-1,j,k  )+                         &
     &           z_r(i  ,j,k+1)-z_r(i-1,j,k+1)
            tl_cff2=tl_z_r(i  ,j,k  )-tl_z_r(i-1,j,k  )+                &
     &              tl_z_r(i  ,j,k+1)-tl_z_r(i-1,j,k+1)
            cff3=z_r(i  ,j,k+1)-z_r(i  ,j,k  )-                         &
     &           z_r(i-1,j,k+1)+z_r(i-1,j,k  )
            tl_cff3=tl_z_r(i  ,j,k+1)-tl_z_r(i  ,j,k  )-                &
     &              tl_z_r(i-1,j,k+1)+tl_z_r(i-1,j,k  )
            gamma=0.125_r8*cff1*cff2*cff3
            tl_gamma=0.125_r8*(tl_cff1*cff2*cff3+                       &
     &                         cff1*(tl_cff2*cff3+                      &
     &                               cff2*tl_cff3))-                    &
# ifdef TL_IOMS
     &               2.0_r8*gamma
# endif

            cff1=(1.0_r8+gamma)*(rho(i,j,k+1)-rho(i-1,j,k+1))+          &
     &           (1.0_r8-gamma)*(rho(i,j,k  )-rho(i-1,j,k  ))
            tl_cff1=tl_gamma*(rho(i,j,k+1)-rho(i-1,j,k+1)-              &
     &                        rho(i,j,k  )+rho(i-1,j,k  ))+             &
     &              (1.0_r8+gamma)*(tl_rho(i  ,j,k+1)-                  &
     &                              tl_rho(i-1,j,k+1))+                 &
     &              (1.0_r8-gamma)*(tl_rho(i  ,j,k  )-                  &
     &                              tl_rho(i-1,j,k  ))-                 &
# ifdef TL_IOMS
     &              gamma*((rho(i,j,k+1)-rho(i-1,j,k+1))-               &
     &                     (rho(i,j,k  )-rho(i-1,j,k  )))
# endif
            cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-                  &
     &              tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i-1,j,k+1)-                  &
     &              tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )
            cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i-1,j,k+1))+          &
     &           (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i-1,j,k  ))
            tl_cff4=tl_gamma*(z_r(i,j,k+1)-z_r(i-1,j,k+1)-              &
     &                        z_r(i,j,k  )+z_r(i-1,j,k  ))+             &
     &              (1.0_r8+gamma)*(tl_z_r(i  ,j,k+1)-                  &
     &                              tl_z_r(i-1,j,k+1))+                 &
     &              (1.0_r8-gamma)*(tl_z_r(i  ,j,k  )-                  &
     &                              tl_z_r(i-1,j,k  ))-                 &
# ifdef TL_IOMS
     &              gamma*((z_r(i,j,k+1)-z_r(i-1,j,k+1))-               &
     &                     (z_r(i,j,k  )-z_r(i-1,j,k  )))
# endif
            phix(i)=phix(i)+                                            &
     &              fac3*(cff1*cff3-cff2*cff4)
            tl_phix(i)=tl_phix(i)+                                      &
     &                 fac3*(tl_cff1*cff3+                              &
     &                       cff1*tl_cff3-                              &
     &                       tl_cff2*cff4-                              &
     &                       cff2*tl_cff4)-                             &
# ifdef TL_IOMS
     &                 fac3*(cff1*cff3-                                 &
     &                       cff2*cff4)
# endif
#else
            cff1=rho(i,j,k+1)-rho(i-1,j,k+1)+                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            tl_cff1=tl_rho(i,j,k+1)-tl_rho(i-1,j,k+1)+                  &
     &              tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
            tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-                  &
     &              tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            cff4=z_r(i,j,k+1)-z_r(i-1,j,k+1)+                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i-1,j,k+1)-                  &
     &              tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )
            tl_cff4=tl_z_r(i,j,k+1)-tl_z_r(i-1,j,k+1)+                  &
     &              tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )
            phix(i)=phix(i)+                                            &
     &              fac3*(cff1*cff3-cff2*cff4)
            tl_phix(i)=tl_phix(i)+                                      &
     &                 fac3*(tl_cff1*cff3+                              &
     &                       cff1*tl_cff3-                              &
     &                       tl_cff2*cff4-                              &
     &                       cff2*tl_cff4)-                             &
# ifdef TL_IOMS
     &                 fac3*(cff1*cff3-                                 &
     &                       cff2*cff4)
# endif
#endif
!>          ru(i,j,k,nrhs)=-0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*             &
!>   &                     phix(i)*on_u(i,j)
!>
            tl_ru(i,j,k,nrhs)=-0.5_r8*on_u(i,j)*                        &
     &                        ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*           &
     &                         phix(i)+                                 &
     &                         (Hz(i,j,k)+Hz(i-1,j,k))*                 &
     &                         tl_phix(i))+                             &
#ifdef TL_IOMS
     &                        0.5_r8*on_u(i,j)*                         &
     &                        (Hz(i,j,k)+Hz(i-1,j,k))*phix(i)
#endif
#ifdef DIAGNOSTICS_UV
!!          DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
#endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Calculate pressure gradient in the ETA-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute surface baroclinic pressure gradient.
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            cff1=z_w(i,j  ,N(ng))-z_r(i,j  ,N(ng))+                     &
     &           z_w(i,j-1,N(ng))-z_r(i,j-1,N(ng))
            tl_cff1=tl_z_w(i,j  ,N(ng))-tl_z_r(i,j  ,N(ng))+            &
     &              tl_z_w(i,j-1,N(ng))-tl_z_r(i,j-1,N(ng))
            phie(i)=fac1*(rho(i,j,N(ng))-rho(i,j-1,N(ng)))*cff1
            tl_phie(i)=fac1*                                            &
     &                 ((tl_rho(i,j,N(ng))-tl_rho(i,j-1,N(ng)))*cff1+   &
     &                  (rho(i,j,N(ng))-rho(i,j-1,N(ng)))*tl_cff1)-     &
#ifdef TL_IOMS
     &                 phie(i)
#endif
#ifdef ATM_PRESS
            phie(i)=phie(i)+fac*(Pair(i,j)-Pair(i,j-1))
#endif
#ifdef RHO_SURF
            phie(i)=phie(i)+                                            &
     &              (fac2+fac1*(rho(i,j,N(ng))+rho(i,j-1,N(ng))))*      &
     &              (z_w(i,j,N(ng))-z_w(i,j-1,N(ng)))
            tl_phie(i)=tl_phie(i)+                                      &
     &                 (fac1*(tl_rho(i,j,N(ng))+tl_rho(i,j-1,N(ng))))*  &
     &                 (z_w(i,j,N(ng))-z_w(i,j-1,N(ng)))+               &
     &                 (fac2+fac1*(rho(i,j,N(ng))+rho(i,j-1,N(ng))))*   &
     &                 (tl_z_w(i,j,N(ng))-tl_z_w(i,j-1,N(ng)))-         &
# ifdef TL_IOMS
     &                 fac1*(rho(i,j,N(ng))+rho(i,j-1,N(ng)))*          &
     &                      (z_w(i,j,N(ng))-z_w(i,j-1,N(ng)))
# endif
#endif
!>          rv(i,j,N(ng),nrhs)=-0.5_r8*(Hz(i,j,N(ng))+Hz(i,j-1,N(ng)))* &
!>   &                         phie(i)*om_v(i,j)
!>
            tl_rv(i,j,N(ng),nrhs)=-0.5_r8*om_v(i,j)*                    &
     &                            ((tl_Hz(i,j  ,N(ng))+                 &
     &                              tl_Hz(i,j-1,N(ng)))*phie(i)+        &
     &                             (Hz(i,j  ,N(ng))+                    &
     &                              Hz(i,j-1,N(ng)))*tl_phie(i))+       &
#ifdef TL_IOMS
     &                            0.5_r8*om_v(i,j)*                     &
     &                            (Hz(i,j  ,N(ng))+                     &
     &                             Hz(i,j-1,N(ng)))*phie(i)
#endif
#ifdef DIAGNOSTICS_UV
!!          DiaRV(i,j,N(ng),nrhs,M3pgrd)=rv(i,j,N(ng),nrhs)
#endif
          END DO
!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate.
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
#ifdef WJ_GRADP
              cff1=1.0_r8/((z_r(i,j  ,k+1)-z_r(i,j  ,k))*               &
     &                     (z_r(i,j-1,k+1)-z_r(i,j-1,k)))
              tl_cff1=-cff1*cff1*((tl_z_r(i,j  ,k+1)-tl_z_r(i,j  ,k))*  &
     &                            (z_r(i,j-1,k+1)-z_r(i,j-1,k))+        &
     &                            (z_r(i,j  ,k+1)-z_r(i,j  ,k))*        &
     &                            (tl_z_r(i,j-1,k+1)-tl_z_r(i,j-1,k)))+ &
#  ifdef TL_IOMS
     &                3.0_r8*cff1
#  endif
              cff2=z_r(i,j  ,k  )-z_r(i,j-1,k  )+                       &
     &             z_r(i,j  ,k+1)-z_r(i,j-1,k+1)
              tl_cff2=tl_z_r(i,j  ,k  )-tl_z_r(i,j-1,k  )+              &
     &                tl_z_r(i,j  ,k+1)-tl_z_r(i,j-1,k+1)
              cff3=z_r(i,j  ,k+1)-z_r(i,j  ,k  )-                       &
     &             z_r(i,j-1,k+1)+z_r(i,j-1,k  )
              tl_cff3=tl_z_r(i,j  ,k+1)-tl_z_r(i,j  ,k  )-              &
     &                tl_z_r(i,j-1,k+1)+tl_z_r(i,j-1,k  )
              gamma=0.125_r8*cff1*cff2*cff3
              tl_gamma=0.125_r8*(tl_cff1*cff2*cff3+                     &
     &                           cff1*(tl_cff2*cff3+                    &
     &                                 cff2*tl_cff3))-                  &
# ifdef TL_IOMS
     &                 2.0_r8*gamma
# endif

              cff1=(1.0_r8+gamma)*(rho(i,j,k+1)-rho(i,j-1,k+1))+        &
     &             (1.0_r8-gamma)*(rho(i,j,k  )-rho(i,j-1,k  ))
              tl_cff1=tl_gamma*(rho(i,j,k+1)-rho(i,j-1,k+1)-            &
     &                          rho(i,j,k  )+rho(i,j-1,k  ))+           &
     &                (1.0_r8+gamma)*(tl_rho(i,j  ,k+1)-                &
     &                                tl_rho(i,j-1,k+1))+               &
     &                (1.0_r8-gamma)*(tl_rho(i,j  ,k  )-                &
     &                                tl_rho(i,j-1,k  ))-               &
# ifdef TL_IOMS
     &                gamma*((rho(i,j,k+1)-rho(i,j-1,k+1))-             &
     &                       (rho(i,j,k  )-rho(i,j-1,k  )))
# endif
              cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-                &
     &                tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
              cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i,j-1,k+1)-                &
     &                tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )
              cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i,j-1,k+1))+        &
     &             (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i,j-1,k  ))
              tl_cff4=tl_gamma*(z_r(i,j,k+1)-z_r(i,j-1,k+1)-            &
     &                          z_r(i,j,k  )+z_r(i,j-1,k  ))+           &
     &                (1.0_r8+gamma)*(tl_z_r(i,j  ,k+1)-                &
     &                                tl_z_r(i,j-1,k+1))+               &
     &                (1.0_r8-gamma)*(tl_z_r(i,j  ,k  )-                &
     &                                tl_z_r(i,j-1,k  ))-               &
# ifdef TL_IOMS
     &                gamma*((z_r(i,j,k+1)-z_r(i,j-1,k+1))-             &
     &                       (z_r(i,j,k  )-z_r(i,j-1,k  )))
# endif
              phie(i)=phie(i)+                                          &
     &                fac3*(cff1*cff3-cff2*cff4)
              tl_phie(i)=tl_phie(i)+                                    &
     &                   fac3*(tl_cff1*cff3+                            &
     &                         cff1*tl_cff3-                            &
     &                         tl_cff2*cff4-                            &
     &                         cff2*tl_cff4)-                           &
# ifdef TL_IOMS
     &                   fac3*(cff1*cff3-                               &
     &                         cff2*cff4)
# endif
#else
              cff1=rho(i,j,k+1)-rho(i,j-1,k+1)+                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              tl_cff1=tl_rho(i,j,k+1)-tl_rho(i,j-1,k+1)+                &
     &                tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
              tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-                &
     &                tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
              cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              cff4=z_r(i,j,k+1)-z_r(i,j-1,k+1)+                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i,j-1,k+1)-                &
     &                tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )
              tl_cff4=tl_z_r(i,j,k+1)-tl_z_r(i,j-1,k+1)+                &
     &                tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )
              phie(i)=phie(i)+                                          &
     &                fac3*(cff1*cff3-cff2*cff4)
              tl_phie(i)=tl_phie(i)+                                    &
     &                   fac3*(tl_cff1*cff3+                            &
     &                         cff1*tl_cff3-                            &
     &                         tl_cff2*cff4-                            &
     &                         cff2*tl_cff4)-                           &
# ifdef TL_IOMS
     &                   fac3*(cff1*cff3-                               &
     &                         cff2*cff4)
# endif
#endif
!>            rv(i,j,k,nrhs)=-0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*           &
!>   &                       phie(i)*om_v(i,j)
!>
              tl_rv(i,j,k,nrhs)=-0.5_r8*om_v(i,j)*                      &
     &                          ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*         &
     &                           phie(i)+                               &
     &                           (Hz(i,j,k)+Hz(i,j-1,k))*               &
     &                           tl_phie(i))+                           &
#ifdef TL_IOMS
     &                          0.5_r8*om_v(i,j)*                       &
     &                          (Hz(i,j,k)+Hz(i,j-1,k))*phie(i)
#endif
#ifdef DIAGNOSTICS_UV
!!            DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
#endif
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE rp_prsgrd_tile
