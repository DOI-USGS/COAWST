      SUBROUTINE ad_prsgrd (ng, tile)
!
!svn $Id: ad_prsgrd31.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine evalutes the adjoint baroclinic hydrostatic pressure   !
!  gradient term using the standard and weighted Jacobians scheme of   !
!  Song and Wright (1997).  Notice  that  horizontal  gradients  are   !
!  computed before of the vertical integration.                        !
!                                                                      !
!  The pressure gradient terms (m4/s2) are loaded into right-hand-     !
!  side arrays "ad_ru" and "ad_rv".                                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Song, Y.T. and D.G. Wright, 1997: A general pressure gradient     !
!          formutlation for numerical ocean models. Part I: Scheme     !
!          design and diagnostic analysis.  DRAFT.                     !
!                                                                      !
!  BASIC STATE variables needed: Hz, rho, z_r, z_w                     !
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
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % ad_z_r,                           &
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
     &                           z_r, ad_z_r,                           &
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
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
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
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: ad_rv(LBi:,LBj:,0:,:)
#else
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
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
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: ad_rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k, kk

      real(r8) :: fac, fac1, fac2, fac3
      real(r8) :: cff1, cff2, cff3, cff4
      real(r8) :: adfac, adfac1, adfac2
      real(r8) :: ad_cff1, ad_cff2, ad_cff3, ad_cff4
#ifdef WJ_GRADP
      real(r8) :: gamma, ad_gamma
#endif

      real(r8), dimension(IminS:ImaxS) :: phie
      real(r8), dimension(IminS:ImaxS) :: phix

      real(r8), dimension(IminS:ImaxS) :: ad_phie
      real(r8), dimension(IminS:ImaxS) :: ad_phix

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_cff3=0.0_r8
      ad_cff4=0.0_r8

#ifdef WJ_GRADP
      ad_gamma=0.0_r8
#endif
      DO i=IminS,ImaxS
        ad_phie(i)=0.0_r8
        ad_phix(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Calculate adjoint pressure gradient in the ETA-direction (m4/s2).
!-----------------------------------------------------------------------
!
#ifdef ATM_PRESS
      fac=100.0_r8/rho0
#endif
      fac1=0.5_r8*g/rho0
      fac2=1000.0_r8*g/rho0
      fac3=0.25_r8*g/rho0

      J_LOOP : DO j=Jstr,Jend
        IF (j.ge.JstrV) THEN
!
!  Compute appropriate BASIC STATE "phie".  Notice that a reverse
!  vertical integration of "phie" is carried out over kk-index.
!
!>        DO k=N(ng)-1,1,-1
!>
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=z_w(i,j  ,N(ng))-z_r(i,j  ,N(ng))+                   &
     &             z_w(i,j-1,N(ng))-z_r(i,j-1,N(ng))
              phie(i)=fac1*(rho(i,j,N(ng))-rho(i,j-1,N(ng)))*cff1
#ifdef ATM_PRESS
              phie(i)=phie(i)+fac*(Pair(i,j)-Pair(i,j-1))
#endif
#ifdef RHO_SURF
              phie(i)=phie(i)+                                          &
     &                (fac2+fac1*(rho(i,j,N(ng))+rho(i,j-1,N(ng))))*    &
     &                (z_w(i,j,N(ng))-z_w(i,j-1,N(ng)))
#endif
            END DO
            DO kk=N(ng)-1,k,-1
              DO i=Istr,Iend
#ifdef WJ_GRADP
                cff1=1.0_r8/((z_r(i,j  ,kk+1)-z_r(i,j  ,kk))*           &
     &                       (z_r(i,j-1,kk+1)-z_r(i,j-1,kk)))
                cff2=z_r(i,j  ,kk  )-z_r(i,j-1,kk  )+                   &
     &               z_r(i,j  ,kk+1)-z_r(i,j-1,kk+1)
                cff3=z_r(i,j  ,kk+1)-z_r(i,j  ,kk  )-                   &
     &               z_r(i,j-1,kk+1)+z_r(i,j-1,kk  )
                gamma=0.125_r8*cff1*cff2*cff3

                cff1=(1.0_r8+gamma)*(rho(i,j,kk+1)-rho(i,j-1,kk+1))+    &
     &               (1.0_r8-gamma)*(rho(i,j,kk  )-rho(i,j-1,kk  ))
                cff2=rho(i,j,kk+1)+rho(i,j-1,kk+1)-                     &
     &               rho(i,j,kk  )-rho(i,j-1,kk  )
                cff3=z_r(i,j,kk+1)+z_r(i,j-1,kk+1)-                     &
     &               z_r(i,j,kk  )-z_r(i,j-1,kk  )
                cff4=(1.0_r8+gamma)*(z_r(i,j,kk+1)-z_r(i,j-1,kk+1))+    &
     &               (1.0_r8-gamma)*(z_r(i,j,kk  )-z_r(i,j-1,kk  ))
                phie(i)=phie(i)+                                        &
     &                  fac3*(cff1*cff3-cff2*cff4)
#else
                cff1=rho(i,j,kk+1)-rho(i,j-1,kk+1)+                     &
     &               rho(i,j,kk  )-rho(i,j-1,kk  )
                cff2=rho(i,j,kk+1)+rho(i,j-1,kk+1)-                     &
     &               rho(i,j,kk  )-rho(i,j-1,kk  )
                cff3=z_r(i,j,kk+1)+z_r(i,j-1,kk+1)-                     &
     &               z_r(i,j,kk  )-z_r(i,j-1,kk  )
                cff4=z_r(i,j,kk+1)-z_r(i,j-1,kk+1)+                     &
     &               z_r(i,j,kk  )-z_r(i,j-1,kk  )
                phie(i)=phie(i)+                                        &
     &                  fac3*(cff1*cff3-cff2*cff4)
#endif
              END DO
            END DO
!
!  Compute interior adjoint baroclinic pressure gradient.  Differentiate
!  and then vertically integrate.
!
            DO i=Istr,Iend
#ifdef DIAGNOSTICS_UV
!!            DiaRV(i,j,k,nrhs,M3pgrd)=rv(i,j,k,nrhs)
#endif
!>            tl_rv(i,j,k,nrhs)=-0.5_r8*om_v(i,j)*                      &
!>   &                          ((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*         &
!>   &                           phie(i)+                               &
!>   &                           (Hz(i,j,k)+Hz(i,j-1,k))*               &
!>   &                           tl_phie(i))
!>
              adfac=-0.5_r8*om_v(i,j)*ad_rv(i,j,k,nrhs)
              adfac1=adfac*phie(i)
              ad_phie(i)=ad_phie(i)+                                    &
     &                   (Hz(i,j,k)+Hz(i,j-1,k))*adfac
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
              ad_rv(i,j,k,nrhs)=0.0_r8
#ifdef WJ_GRADP
              cff1=1.0_r8/((z_r(i,j  ,k+1)-z_r(i,j  ,k))*               &
     &                     (z_r(i,j-1,k+1)-z_r(i,j-1,k)))
              cff2=z_r(i,j  ,k  )-z_r(i,j-1,k  )+                       &
     &             z_r(i,j  ,k+1)-z_r(i,j-1,k+1)
              cff3=z_r(i,j  ,k+1)-z_r(i,j  ,k  )-                       &
     &             z_r(i,j-1,k+1)+z_r(i,j-1,k  )
              gamma=0.125_r8*cff1*cff2*cff3

              cff1=(1.0_r8+gamma)*(rho(i,j,k+1)-rho(i,j-1,k+1))+        &
     &             (1.0_r8-gamma)*(rho(i,j,k  )-rho(i,j-1,k  ))
              cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i,j-1,k+1))+        &
     &             (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i,j-1,k  ))
!>            tl_phie(i)=tl_phie(i)+                                    &
!>   &                   fac3*(tl_cff1*cff3+                            &
!>   &                         cff1*tl_cff3-                            &
!>   &                         tl_cff2*cff4-                            &
!>   &                         cff2*tl_cff4)
!>
              adfac=fac3*ad_phie(i)
              ad_cff1=ad_cff1+cff3*adfac
              ad_cff2=ad_cff2-cff4*adfac
              ad_cff3=ad_cff3+cff1*adfac
              ad_cff4=ad_cff4-cff2*adfac
!>            tl_cff4=tl_gamma*(z_r(i,j,k+1)-z_r(i,j-1,k+1)-            &
!>   &                          z_r(i,j,k  )+z_r(i,j-1,k  ))+           &
!>   &                (1.0_r8+gamma)*(tl_z_r(i,j  ,k+1)-                &
!>   &                                tl_z_r(i,j-1,k+1))+               &
!>   &                (1.0_r8-gamma)*(tl_z_r(i,j  ,k  )-                &
!>   &                                tl_z_r(i,j-1,k  ))
!>            tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i,j-1,k+1)-                &
!>   &                tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )
!>
              adfac1=(1.0_r8+gamma)*ad_cff4
              adfac2=(1.0_r8-gamma)*ad_cff4
              ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-adfac2-ad_cff3
              ad_z_r(i,j  ,k  )=ad_z_r(i,j  ,k  )+adfac2-ad_cff3
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)-adfac1+ad_cff3
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+adfac1+ad_cff3
              ad_gamma=ad_gamma+                                        &
     &                 (z_r(i,j,k+1)-z_r(i,j-1,k+1)-                    &
     &                  z_r(i,j,k  )+z_r(i,j-1,k  ))*ad_cff4
              ad_cff4=0.0_r8
              ad_cff3=0.0_r8
!>            tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-                &
!>   &                tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
!>            tl_cff1=tl_gamma*(rho(i,j,k+1)-rho(i,j-1,k+1)-            &
!>   &                          rho(i,j,k  )+rho(i,j-1,k  ))+           &
!>   &                (1.0_r8+gamma)*(tl_rho(i,j  ,k+1)-                &
!>   &                                tl_rho(i,j-1,k+1))+               &
!>   &                (1.0_r8-gamma)*(tl_rho(i,j  ,k  )-                &
!>   &                                tl_rho(i,j-1,k  ))
!>
              adfac1=(1.0_r8+gamma)*ad_cff1
              adfac2=(1.0_r8-gamma)*ad_cff1
              ad_rho(i,j-1,k  )=ad_rho(i,j-1,k  )-adfac2-ad_cff2
              ad_rho(i,j  ,k  )=ad_rho(i,j  ,k  )+adfac2-ad_cff2
              ad_rho(i,j-1,k+1)=ad_rho(i,j-1,k+1)-adfac1+ad_cff2
              ad_rho(i,j  ,k+1)=ad_rho(i,j  ,k+1)+adfac1+ad_cff2
              ad_gamma=ad_gamma+                                        &
     &                  (rho(i,j,k+1)-rho(i,j-1,k+1)-                   &
     &                   rho(i,j,k  )+rho(i,j-1,k  ))*ad_cff1
              ad_cff2=0.0_r8
              ad_cff1=0.0_r8
!
              cff1=1.0_r8/((z_r(i,j  ,k+1)-z_r(i,j  ,k))*               &
     &                     (z_r(i,j-1,k+1)-z_r(i,j-1,k)))
              cff2=z_r(i,j  ,k  )-z_r(i,j-1,k  )+                       &
     &             z_r(i,j  ,k+1)-z_r(i,j-1,k+1)
              cff3=z_r(i,j  ,k+1)-z_r(i,j  ,k  )-                       &
     &             z_r(i,j-1,k+1)+z_r(i,j-1,k  )

!>            tl_gamma=0.125_r8*(tl_cff1*cff2*cff3+                     &
!>   &                           cff1*(tl_cff2*cff3+                    &
!>   &                                 cff2*tl_cff3))
!>
              adfac=0.125_r8*ad_gamma
              adfac1=adfac*cff1
              ad_cff3=ad_cff3+cff2*adfac1
              ad_cff2=ad_cff2+cff3*adfac1
              ad_cff1=ad_cff1+cff2*cff3*adfac
              ad_gamma=0.0_r8
!>            tl_cff3=tl_z_r(i,j  ,k+1)-tl_z_r(i,j  ,k  )-              &
!>   &                tl_z_r(i,j-1,k+1)+tl_z_r(i,j-1,k  )
!>            tl_cff2=tl_z_r(i,j  ,k  )-tl_z_r(i,j-1,k  )+              &
!>   &                tl_z_r(i,j  ,k+1)-tl_z_r(i,j-1,k+1)
!>
              ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-ad_cff2+ad_cff3
              ad_z_r(i,j  ,k  )=ad_z_r(i,j  ,k  )+ad_cff2-ad_cff3
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)-ad_cff2-ad_cff3
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+ad_cff2+ad_cff3
              ad_cff3=0.0_r8
              ad_cff2=0.0_r8
!>            tl_cff1=-cff1*cff1*((tl_z_r(i,j  ,k+1)-tl_z_r(i,j  ,k))*  &
!>   &                            (z_r(i,j-1,k+1)-z_r(i,j-1,k))+        &
!>   &                            (z_r(i,j  ,k+1)-z_r(i,j  ,k))*        &
!>   &                            (tl_z_r(i,j-1,k+1)-tl_z_r(i,j-1,k)))
!>
              adfac=-cff1*cff1*ad_cff1
              adfac1=adfac*(z_r(i,j-1,k+1)-z_r(i,j-1,k))
              adfac2=adfac*(z_r(i,j  ,k+1)-z_r(i,j  ,k))
              ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-adfac2
              ad_z_r(i,j  ,k  )=ad_z_r(i,j  ,k  )-adfac1
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)+adfac2
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+adfac1
              ad_cff1=0.0_r8
#else
!
              cff1=rho(i,j,k+1)-rho(i,j-1,k+1)+                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              cff4=z_r(i,j,k+1)-z_r(i,j-1,k+1)+                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )

!>            tl_phie(i)=tl_phie(i)+                                    &
!>   &                   fac3*(tl_cff1*cff3+                            &
!>   &                         cff1*tl_cff3-                            &
!>   &                         tl_cff2*cff4-                            &
!>   &                         cff2*tl_cff4)
!>
              adfac=fac3*ad_phie(i)
              ad_cff1=ad_cff1+cff3*adfac
              ad_cff2=ad_cff2-cff4*adfac
              ad_cff3=ad_cff3+cff1*adfac
              ad_cff4=ad_cff4-cff2*adfac
!>            tl_cff4=tl_z_r(i,j,k+1)-tl_z_r(i,j-1,k+1)+                &
!>   &                tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )
!>            tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i,j-1,k+1)-                &
!>   &                tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )
!>
              ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-ad_cff3-ad_cff4
              ad_z_r(i,j  ,k  )=ad_z_r(i,j  ,k  )-ad_cff3+ad_cff4
              ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)+ad_cff3-ad_cff4
              ad_z_r(i,j  ,k+1)=ad_z_r(i,j  ,k+1)+ad_cff3+ad_cff4
              ad_cff4=0.0_r8
              ad_cff3=0.0_r8
!>            tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-                &
!>   &                tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
!>            tl_cff1=tl_rho(i,j,k+1)-tl_rho(i,j-1,k+1)+                &
!>   &                tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
!>
              ad_rho(i,j-1,k  )=ad_rho(i,j-1,k  )-ad_cff1-ad_cff2
              ad_rho(i,j  ,k  )=ad_rho(i,j  ,k  )+ad_cff1-ad_cff2
              ad_rho(i,j-1,k+1)=ad_rho(i,j-1,k+1)-ad_cff1+ad_cff2
              ad_rho(i,j  ,k+1)=ad_rho(i,j  ,k+1)+ad_cff1+ad_cff2
              ad_cff2=0.0_r8
              ad_cff1=0.0_r8
#endif
            END DO
          END DO
!
!  Compute surface adjoint baroclinic pressure gradient.
!
          DO i=Istr,Iend
            cff1=z_w(i,j  ,N(ng))-z_r(i,j  ,N(ng))+                     &
     &           z_w(i,j-1,N(ng))-z_r(i,j-1,N(ng))
            phie(i)=fac1*(rho(i,j,N(ng))-rho(i,j-1,N(ng)))*cff1
#ifdef ATM_PRESS
            phie(i)=phie(i)+fac*(Pair(i,j)-Pair(i,j-1))
#endif
#ifdef RHO_SURF
            phie(i)=phie(i)+                                            &
     &              (fac2+fac1*(rho(i,j,N(ng))+rho(i,j-1,N(ng))))*      &
     &              (z_w(i,j,N(ng))-z_w(i,j-1,N(ng)))
#endif
# ifdef DIAGNOSTICS_UV
!!          DiaRV(i,j,N(ng),nrhs,M3pgrd)=rv(i,j,N(ng),nrhs)
# endif
!>          tl_rv(i,j,N(ng),nrhs)=-0.5_r8*om_v(i,j)*                    &
!>   &                            ((tl_Hz(i,j  ,N(ng))+                 &
!>   &                              tl_Hz(i,j-1,N(ng)))*phie(i)+        &
!>   &                             (Hz(i,j  ,N(ng))+                    &
!>   &                              Hz(i,j-1,N(ng)))*tl_phie(i))
!>
            adfac=-0.5_r8*om_v(i,j)*ad_rv(i,j,N(ng),nrhs)
            adfac1=adfac*phie(i)
            ad_phie(i)=ad_phie(i)+(Hz(i,j  ,N(ng))+                     &
     &                             Hz(i,j-1,N(ng)))*adfac
            ad_Hz(i,j-1,N(ng))=ad_Hz(i,j-1,N(ng))+adfac1
            ad_Hz(i,j  ,N(ng))=ad_Hz(i,j  ,N(ng))+adfac1
            ad_rv(i,j,N(ng),nrhs)=0.0
#ifdef RHO_SURF
!>          tl_phie(i)=tl_phie(i)+                                      &
!>   &                 (fac1*(tl_rho(i,j,N(ng))+tl_rho(i,j-1,N(ng))))*  &
!>   &                 (z_w(i,j,N(ng))-z_w(i,j-1,N(ng)))+               &
!>   &                 (fac2+fac1*(rho(i,j,N(ng))+rho(i,j-1,N(ng))))*   &
!>   &                 (tl_z_w(i,j,N(ng))-tl_z_w(i,j-1,N(ng)))
!>
            adfac1=fac1*(z_w(i,j,N(ng))-z_w(i,j-1,N(ng)))*              &
     &             ad_phie(i)
            adfac2=(fac2+fac1*(rho(i,j,N(ng))+rho(i,j-1,N(ng))))*       &
     &             ad_phie(i)
            ad_rho(i,j-1,N(ng))=ad_rho(i,j-1,N(ng))+adfac1
            ad_rho(i,j  ,N(ng))=ad_rho(i,j  ,N(ng))+adfac1
            ad_z_w(i,j-1,N(ng))=ad_z_w(i,j-1,N(ng))-adfac2
            ad_z_w(i,j  ,N(ng))=ad_z_w(i,j  ,N(ng))+adfac2
#endif
!>          tl_phie(i)=fac1*                                            &
!>   &                 ((tl_rho(i,j,N(ng))-tl_rho(i,j-1,N(ng)))*cff1+   &
!>   &                  (rho(i,j,N(ng))-rho(i,j-1,N(ng)))*tl_cff1)
!>
            adfac=fac1*ad_phie(i)
            adfac1=adfac*cff1
            ad_rho(i,j-1,N(ng))=ad_rho(i,j-1,N(ng))-adfac1
            ad_rho(i,j  ,N(ng))=ad_rho(i,j  ,N(ng))+adfac1
            ad_cff1=ad_cff1+                                            &
     &              (rho(i,j,N(ng))-rho(i,j-1,N(ng)))*adfac
            ad_phie(i)=0.0_r8
!>          tl_cff1=tl_z_w(i,j  ,N(ng))-tl_z_r(i,j  ,N(ng))+            &
!>   &              tl_z_w(i,j-1,N(ng))-tl_z_r(i,j-1,N(ng))
!>
            ad_z_r(i,j-1,N(ng))=ad_z_r(i,j-1,N(ng))-ad_cff1
            ad_z_r(i,j  ,N(ng))=ad_z_r(i,j  ,N(ng))-ad_cff1
            ad_z_w(i,j-1,N(ng))=ad_z_w(i,j-1,N(ng))+ad_cff1
            ad_z_w(i,j  ,N(ng))=ad_z_w(i,j  ,N(ng))+ad_cff1
            ad_cff1=0.0_r8
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Calculate adjoint pressure gradient in the XI-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute appropriate BASIC STATE "phix".  Notice that a reverse
!  vertical integration of "phix" is carried out over kk-index.
!
!>      DO k=N(ng)-1,1,-1
!>
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff1=z_w(i  ,j,N(ng))-z_r(i  ,j,N(ng))+                     &
     &           z_w(i-1,j,N(ng))-z_r(i-1,j,N(ng))
            phix(i)=fac1*(rho(i,j,N(ng))-rho(i-1,j,N(ng)))*cff1
#ifdef ATM_PRESS
            phix(i)=phix(i)+fac*(Pair(i,j)-Pair(i-1,j))
#endif
#ifdef RHO_SURF
            phix(i)=phix(i)+                                            &
     &              (fac2+fac1*(rho(i,j,N(ng))+rho(i-1,j,N(ng))))*      &
     &              (z_w(i,j,N(ng))-z_w(i-1,j,N(ng)))
#endif
          END DO
          DO kk=N(ng)-1,k,-1
            DO i=IstrU,Iend
#ifdef WJ_GRADP
              cff1=1.0_r8/((z_r(i  ,j,kk+1)-z_r(i  ,j,kk))*             &
     &                     (z_r(i-1,j,kk+1)-z_r(i-1,j,kk)))
              cff2=z_r(i  ,j,kk  )-z_r(i-1,j,kk )+                      &
     &             z_r(i  ,j,kk+1)-z_r(i-1,j,k+1)
              cff3=z_r(i  ,j,kk+1)-z_r(i  ,j,kk )-                      &
     &             z_r(i-1,j,kk+1)+z_r(i-1,j,kk )
              gamma=0.125_r8*cff1*cff2*cff3

              cff1=(1.0_r8+gamma)*(rho(i,j,kk+1)-rho(i-1,j,kk+1))+      &
     &             (1.0_r8-gamma)*(rho(i,j,kk  )-rho(i-1,j,kk  ))
              cff2=rho(i,j,kk+1)+rho(i-1,j,kk+1)-                       &
     &             rho(i,j,kk  )-rho(i-1,j,kk  )
              cff3=z_r(i,j,kk+1)+z_r(i-1,j,kk+1)-                       &
     &             z_r(i,j,kk  )-z_r(i-1,j,kk  )
              cff4=(1.0_r8+gamma)*(z_r(i,j,kk+1)-z_r(i-1,j,kk+1))+      &
     &             (1.0_r8-gamma)*(z_r(i,j,kk  )-z_r(i-1,j,kk  ))
              phix(i)=phix(i)+                                          &
     &                fac3*(cff1*cff3-cff2*cff4)
#else
              cff1=rho(i,j,kk+1)-rho(i-1,j,kk+1)+                       &
     &             rho(i,j,kk  )-rho(i-1,j,kk  )
              cff2=rho(i,j,kk+1)+rho(i-1,j,kk+1)-                       &
     &             rho(i,j,kk  )-rho(i-1,j,kk  )
              cff3=z_r(i,j,kk+1)+z_r(i-1,j,kk+1)-                       &
     &             z_r(i,j,kk  )-z_r(i-1,j,kk  )
              cff4=z_r(i,j,kk+1)-z_r(i-1,j,kk+1)+                       &
     &             z_r(i,j,kk  )-z_r(i-1,j,kk  )
              phix(i)=phix(i)+                                          &
     &                fac3*(cff1*cff3-cff2*cff4)
#endif
            END DO
          END DO
!
!  Compute interior adjoint baroclinic pressure gradient.  Differentiate
!  and then vertically integrate.
!
          DO i=IstrU,Iend
# ifdef DIAGNOSTICS_UV
!!          DiaRU(i,j,k,nrhs,M3pgrd)=ru(i,j,k,nrhs)
# endif
!>          tl_ru(i,j,k,nrhs)=-0.5_r8*on_u(i,j)*                        &
!>   &                        ((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*           &
!>   &                         phix(i)+                                 &
!>   &                         (Hz(i,j,k)+Hz(i-1,j,k))*                 &
!>   &                         tl_phix(i))
!>
            adfac=-0.5_r8*on_u(i,j)*ad_ru(i,j,k,nrhs)
            adfac1=adfac*phix(i)
            ad_phix(i)=ad_phix(i)+                                      &
     &                (Hz(i,j,k)+Hz(i-1,j,k))*adfac
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
            ad_ru(i,j,k,nrhs)=0.0
#ifdef WJ_GRADP
            cff1=1.0_r8/((z_r(i  ,j,k+1)-z_r(i  ,j,k))*                 &
     &                   (z_r(i-1,j,k+1)-z_r(i-1,j,k)))
            cff2=z_r(i  ,j,k  )-z_r(i-1,j,k  )+                         &
     &           z_r(i  ,j,k+1)-z_r(i-1,j,k+1)
            cff3=z_r(i  ,j,k+1)-z_r(i  ,j,k  )-                         &
     &           z_r(i-1,j,k+1)+z_r(i-1,j,k  )
            gamma=0.125_r8*cff1*cff2*cff3

            cff1=(1.0_r8+gamma)*(rho(i,j,k+1)-rho(i-1,j,k+1))+          &
     &           (1.0_r8-gamma)*(rho(i,j,k  )-rho(i-1,j,k  ))
            cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            cff4=(1.0_r8+gamma)*(z_r(i,j,k+1)-z_r(i-1,j,k+1))+          &
     &           (1.0_r8-gamma)*(z_r(i,j,k  )-z_r(i-1,j,k  ))
!>          tl_phix(i)=tl_phix(i)+                                      &
!>   &                 fac3*(tl_cff1*cff3+                              &
!>   &                       cff1*tl_cff3-                              &
!>   &                       tl_cff2*cff4-                              &
!>   &                       cff2*tl_cff4)
!>
            adfac=fac3*ad_phix(i)
            ad_cff1=ad_cff1+cff3*adfac
            ad_cff2=ad_cff2-cff4*adfac
            ad_cff3=ad_cff3+cff1*adfac
            ad_cff4=ad_cff4-cff2*adfac
!>          tl_cff4=tl_gamma*(z_r(i,j,k+1)-z_r(i-1,j,k+1)-              &
!>   &                        z_r(i,j,k  )+z_r(i-1,j,k  ))+             &
!>   &              (1.0_r8+gamma)*(tl_z_r(i  ,j,k+1)-                  &
!>   &                              tl_z_r(i-1,j,k+1))+                 &
!>   &              (1.0_r8-gamma)*(tl_z_r(i  ,j,k  )-                  &
!>   &                              tl_z_r(i-1,j,k  ))
!>          tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i-1,j,k+1)-                  &
!>   &              tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )
!>
            adfac1=(1.0_r8+gamma)*ad_cff4
            adfac2=(1.0_r8-gamma)*ad_cff4
            ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-adfac2-ad_cff3
            ad_z_r(i  ,j,k  )=ad_z_r(i  ,j,k  )+adfac2-ad_cff3
            ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)-adfac1+ad_cff3
            ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+adfac1+ad_cff3
            ad_gamma=ad_gamma+                                          &
     &               (z_r(i,j,k+1)-z_r(i-1,j,k+1)-                      &
     &                z_r(i,j,k  )+z_r(i-1,j,k  ))*ad_cff4
            ad_cff4=0.0_r8
            ad_cff3=0.0_r8
!>   &      tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-                  &
!>   &              tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
!>          tl_cff1=tl_gamma*(rho(i,j,k+1)-rho(i-1,j,k+1)-              &
!>   &                        rho(i,j,k  )+rho(i-1,j,k  ))+             &
!>   &              (1.0_r8+gamma)*(tl_rho(i  ,j,k+1)-                  &
!>   &                              tl_rho(i-1,j,k+1))+                 &
!>   &              (1.0_r8-gamma)*(tl_rho(i  ,j,k  )-                  &
!>   &                              tl_rho(i-1,j,k  ))
!>
            adfac1=(1.0_r8+gamma)*ad_cff1
            adfac2=(1.0_r8-gamma)*ad_cff1
            ad_rho(i-1,j,k  )=ad_rho(i-1,j,k  )-adfac2-ad_cff2
            ad_rho(i  ,j,k  )=ad_rho(i  ,j,k  )+adfac2-ad_cff2
            ad_rho(i-1,j,k+1)=ad_rho(i-1,j,k+1)-adfac1+ad_cff2
            ad_rho(i  ,j,k+1)=ad_rho(i  ,j,k+1)+adfac1+ad_cff2
            ad_gamma=ad_gamma+                                          &
     &               (rho(i,j,k+1)-rho(i-1,j,k+1)-                      &
     &                rho(i,j,k  )+rho(i-1,j,k  ))*ad_cff1
            ad_cff2=0.0_r8
            ad_cff1=0.0_r8
!
            cff1=1.0_r8/((z_r(i  ,j,k+1)-z_r(i  ,j,k))*                 &
     &                   (z_r(i-1,j,k+1)-z_r(i-1,j,k)))
            cff2=z_r(i  ,j,k  )-z_r(i-1,j,k  )+                         &
     &           z_r(i  ,j,k+1)-z_r(i-1,j,k+1)
            cff3=z_r(i  ,j,k+1)-z_r(i  ,j,k  )-                         &
     &           z_r(i-1,j,k+1)+z_r(i-1,j,k  )

!>          tl_gamma=0.125_r8*(tl_cff1*cff2*cff3+                       &
!>   &                         cff1*(tl_cff2*cff3+                      &
!>   &                               cff2*tl_cff3))
!>
            adfac=0.125_r8*ad_gamma
            adfac1=adfac*cff1
            ad_cff3=ad_cff3+cff2*adfac1
            ad_cff2=ad_cff2+cff3*adfac1
            ad_cff1=ad_cff1+cff2*cff3*adfac
            ad_gamma=0.0_r8
!>   &      tl_cff3=tl_z_r(i  ,j,k+1)-tl_z_r(i  ,j,k  )-                &
!>   &              tl_z_r(i-1,j,k+1)+tl_z_r(i-1,j,k  )
!>          tl_cff2=tl_z_r(i  ,j,k  )-tl_z_r(i-1,j,k  )+                &
!>   &              tl_z_r(i  ,j,k+1)-tl_z_r(i-1,j,k+1)
!>
            ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-ad_cff2+ad_cff3
            ad_z_r(i  ,j,k  )=ad_z_r(i  ,j,k  )+ad_cff2-ad_cff3
            ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)-ad_cff2-ad_cff3
            ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+ad_cff2+ad_cff3
            ad_cff3=0.0
            ad_cff2=0.0
!>          tl_cff1=-cff1*cff1*((tl_z_r(i  ,j,k+1)-tl_z_r(i  ,j,k))*    &
!>   &                          (z_r(i-1,j,k+1)-z_r(i-1,j,k))+          &
!>   &                          (z_r(i  ,j,k+1)-z_r(i  ,j,k))*          &
!>   &                          (tl_z_r(i-1,j,k+1)-tl_z_r(i-1,j,k)))
!>
            adfac=-cff1*cff1*ad_cff1
            adfac1=adfac*(z_r(i-1,j,k+1)-z_r(i-1,j,k))
            adfac2=adfac*(z_r(i  ,j,k+1)-z_r(i  ,j,k))
            ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-adfac2
            ad_z_r(i  ,j,k  )=ad_z_r(i  ,j,k  )-adfac1
            ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)+adfac2
            ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+adfac1
            ad_cff1=0.0_r8
#else
            cff1=rho(i,j,k+1)-rho(i-1,j,k+1)+                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            cff4=z_r(i,j,k+1)-z_r(i-1,j,k+1)+                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
!>          tl_phix(i)=tl_phix(i)+                                      &
!>   &                 fac3*(tl_cff1*cff3+                              &
!>   &                       cff1*tl_cff3-                              &
!>   &                       tl_cff2*cff4-                              &
!>   &                       cff2*tl_cff4)
!>
            adfac=fac3*ad_phix(i)
            ad_cff1=ad_cff1+cff3*adfac
            ad_cff2=ad_cff2-cff4*adfac
            ad_cff3=ad_cff3+cff1*adfac
            ad_cff4=ad_cff4-cff2*adfac
!>          tl_cff4=tl_z_r(i,j,k+1)-tl_z_r(i-1,j,k+1)+                  &
!>   &              tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )
!>          tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i-1,j,k+1)-                  &
!>   &              tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )
!>
            ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-ad_cff3-ad_cff4
            ad_z_r(i  ,j,k  )=ad_z_r(i  ,j,k  )-ad_cff3+ad_cff4
            ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)+ad_cff3-ad_cff4
            ad_z_r(i  ,j,k+1)=ad_z_r(i  ,j,k+1)+ad_cff3+ad_cff4
            ad_cff4=0.0_r8
            ad_cff3=0.0_r8
!>          tl_cff1=tl_rho(i,j,k+1)-tl_rho(i-1,j,k+1)+                  &
!>   &              tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
!>          tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-                  &
!>   &              tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
!>
            ad_rho(i-1,j,k  )=ad_rho(i-1,j,k  )-ad_cff2-ad_cff1
            ad_rho(i  ,j,k  )=ad_rho(i  ,j,k  )-ad_cff2+ad_cff1
            ad_rho(i-1,j,k+1)=ad_rho(i-1,j,k+1)+ad_cff2-ad_cff1
            ad_rho(i  ,j,k+1)=ad_rho(i  ,j,k+1)+ad_cff2+ad_cff1
            ad_cff2=0.0_r8
            ad_cff1=0.0_r8
#endif
          END DO
        END DO
!
!  Compute surface adjoint baroclinic pressure gradient.
!
        DO i=IstrU,Iend
          cff1=z_w(i  ,j,N(ng))-z_r(i  ,j,N(ng))+                       &
     &         z_w(i-1,j,N(ng))-z_r(i-1,j,N(ng))
          phix(i)=fac1*(rho(i,j,N(ng))-rho(i-1,j,N(ng)))*cff1
#ifdef ATM_PRESS
          phix(i)=phix(i)+fac*(Pair(i,j)-Pair(i-1,j))
#endif
#ifdef RHO_SURF
          phix(i)=phix(i)+                                              &
     &            (fac2+fac1*(rho(i,j,N(ng))+rho(i-1,j,N(ng))))*        &
     &            (z_w(i,j,N(ng))-z_w(i-1,j,N(ng)))
#endif
#ifdef DIAGNOSTICS_UV
!!        DiaRU(i,j,N(ng),nrhs,M3pgrd)=ru(i,j,N(ng),nrhs)
#endif
!>        tl_ru(i,j,N(ng),nrhs)=-0.5_r8*on_u(i,j)*                      &
!>   &                          ((tl_Hz(i  ,j,N(ng))+                   &
!>   &                            tl_Hz(i-1,j,N(ng)))*phix(i)+          &
!>   &                           (Hz(i  ,j,N(ng))+                      &
!>   &                            Hz(i-1,j,N(ng)))*tl_phix(i))
!>
          adfac=-0.5_r8*on_u(i,j)*ad_ru(i,j,N(ng),nrhs)
          adfac1=adfac*phix(i)
          ad_phix(i)=ad_phix(i)+(Hz(i  ,j,N(ng))+                       &
     &                           Hz(i-1,j,N(ng)))*adfac
          ad_Hz(i-1,j,N(ng))=ad_Hz(i-1,j,N(ng))+adfac1
          ad_Hz(i  ,j,N(ng))=ad_Hz(i  ,j,N(ng))+adfac1
          ad_ru(i,j,N(ng),nrhs)=0.0_r8
#ifdef RHO_SURF
!>        tl_phix(i)=tl_phix(i)+                                        &
!>   &               (fac1*(tl_rho(i,j,N(ng))+tl_rho(i-1,j,N(ng))))*    &
!>   &               (z_w(i,j,N(ng))-z_w(i-1,j,N(ng)))+                 &
!>   &               (fac2+fac1*(rho(i,j,N(ng))+rho(i-1,j,N(ng))))*
!>   &               (tl_z_w(i,j,N(ng))-tl_z_w(i-1,j,N(ng)))
!>
          adfac1=fac1*(z_w(i,j,N(ng))-z_w(i-1,j,N(ng)))*                &
     &           ad_phix(i)
          adfac2=(fac2+fac1*(rho(i,j,N(ng))+rho(i-1,j,N(ng))))*         &
     &           ad_phix(i)
          ad_rho(i-1,j,N(ng))=ad_rho(i-1,j,N(ng))+adfac1
          ad_rho(i  ,j,N(ng))=ad_rho(i  ,j,N(ng))+adfac1
          ad_z_w(i-1,j,N(ng))=ad_z_w(i-1,j,N(ng))-adfac2
          ad_z_w(i  ,j,N(ng))=ad_z_w(i  ,j,N(ng))+adfac2
#endif
!>        tl_phix(i)=fac1*                                              &
!>   &               ((tl_rho(i,j,N(ng))-tl_rho(i-1,j,N(ng)))*cff1+     &
!>   &                (rho(i,j,N(ng))-rho(i-1,j,N(ng)))*tl_cff1)
!>
          adfac=fac1*ad_phix(i)
          adfac1=adfac*cff1
          ad_rho(i-1,j,N(ng))=ad_rho(i-1,j,N(ng))-adfac1
          ad_rho(i  ,j,N(ng))=ad_rho(i  ,j,N(ng))+adfac1
          ad_cff1=ad_cff1+                                              &
     &            (rho(i,j,N(ng))-rho(i-1,j,N(ng)))*adfac
          ad_phix(i)=0.0
!>        tl_cff1=tl_z_w(i  ,j,N(ng))-tl_z_r(i  ,j,N(ng))+              &
!>   &            tl_z_w(i-1,j,N(ng))-tl_z_r(i-1,j,N(ng))
!>
          ad_z_r(i-1,j,N(ng))=ad_z_r(i-1,j,N(ng))-ad_cff1
          ad_z_r(i  ,j,N(ng))=ad_z_r(i  ,j,N(ng))-ad_cff1
          ad_z_w(i-1,j,N(ng))=ad_z_w(i-1,j,N(ng))+ad_cff1
          ad_z_w(i  ,j,N(ng))=ad_z_w(i  ,j,N(ng))+ad_cff1
          ad_cff1=0.0_r8
        END DO
      END DO J_LOOP
      RETURN
      END SUBROUTINE ad_prsgrd_tile
