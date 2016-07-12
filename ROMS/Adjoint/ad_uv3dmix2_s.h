      SUBROUTINE ad_uv3dmix2 (ng, tile)
!
!svn $Id: ad_uv3dmix2_s.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine computes adjoint harmonic mixing of momentum,  along   !
!  constant S-surfaces, from the horizontal divergence of the stress   !
!  tensor.  A transverse isotropy is assumed so the stress tensor is   !
!  splitted into vertical and horizontal subtensors.                   !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!      Wajsowicz, R.C, 1993: A consistent formulation of the           !
!         anisotropic stress tensor for use in models of the           !
!         large-scale ocean circulation, JCP, 105, 333-338.            !
!                                                                      !
!      Sadourny, R. and K. Maynard, 1997: Formulations of              !
!         lateral diffusion in geophysical fluid dynamics              !
!         models, In Numerical Methods of Atmospheric and              !
!         Oceanic Modelling. Lin, Laprise, and Ritchie,                !
!         Eds., NRC Research Press, 547-556.                           !
!                                                                      !
!      Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic              !
!         friction with a Smagorinsky-like viscosity for               !
!         use in large-scale eddy-permitting ocean models,             !
!         Monthly Weather Rev., 128, 8, 2935-2946.                     !
!                                                                      !
!  BASIC STATE variables needed: visc2, u, v, Hz                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_coupling
!!#ifdef DIAGNOSTICS_UV
!!    USE mod_diags
!!#endif
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
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iADM, 30)
#endif
      CALL ad_uv3dmix2_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs(ng), nnew(ng),                        &
#ifdef MASKING
     &                       GRID(ng) % pmask,                          &
#endif
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % ad_Hz,                          &
     &                       GRID(ng) % om_p,                           &
     &                       GRID(ng) % om_r,                           &
     &                       GRID(ng) % on_p,                           &
     &                       GRID(ng) % on_r,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pmon_p,                         &
     &                       GRID(ng) % pmon_r,                         &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % pnom_p,                         &
     &                       GRID(ng) % pnom_r,                         &
     &                       MIXING(ng) % visc2_p,                      &
     &                       MIXING(ng) % visc2_r,                      &
!!#ifdef DIAGNOSTICS_UV
!!   &                       DIAGS(ng) % DiaRUfrc,                      &
!!   &                       DIAGS(ng) % DiaRVfrc,                      &
!!   &                       DIAGS(ng) % DiaU3wrk,                      &
!!   &                       DIAGS(ng) % DiaV3wrk,                      &
!!#endif
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       COUPLING(ng) % ad_rufrc,                   &
     &                       COUPLING(ng) % ad_rvfrc,                   &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v)
#ifdef PROFILE
      CALL wclock_off (ng, iADM, 30)
#endif
      RETURN
      END SUBROUTINE ad_uv3dmix2

!
!***********************************************************************
      SUBROUTINE ad_uv3dmix2_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs, nnew,                          &
#ifdef MASKING
     &                             pmask,                               &
#endif
     &                             Hz, ad_Hz,                           &
     &                             om_p, om_r, on_p, on_r,              &
     &                             pm, pmon_p, pmon_r,                  &
     &                             pn, pnom_p, pnom_r,                  &
     &                             visc2_p, visc2_r,                    &
!!#ifdef DIAGNOSTICS_UV
!!   &                             DiaRUfrc, DiaRVfrc,                  &
!!   &                             DiaU3wrk, DiaV3wrk,                  &
!!#endif
     &                             u, v,                                &
     &                             ad_rufrc, ad_rvfrc,                  &
     &                             ad_u, ad_v)
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
      integer, intent(in) :: nrhs, nnew

#ifdef ASSUMED_SHAPE
#ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
#endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pmon_p(LBi:,LBj:)
      real(r8), intent(in) :: pmon_r(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: pnom_p(LBi:,LBj:)
      real(r8), intent(in) :: pnom_r(LBi:,LBj:)
      real(r8), intent(in) :: visc2_p(LBi:,LBj:)
      real(r8), intent(in) :: visc2_r(LBi:,LBj:)

      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)

      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#else
#ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmon_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pnom_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc2_r(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)

      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: cff, ad_cff, ad_cff1, ad_cff2
      real(r8) :: adfac, adfac1, adfac2, adfac3, adfac4

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFx

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8

      ad_UFe(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_VFe(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_UFx(IminS:ImaxS,JminS:JmaxS)=0.0_r8
      ad_VFx(IminS:ImaxS,JminS:JmaxS)=0.0_r8
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic viscosity along constant S-surfaces.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
! Time-step harmonic, S-surfaces viscosity term.  Notice that momentum
! at this stage is HzU and HzV and has m2/s units.  Add contribution for
! barotropic forcing terms.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!>          tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+tl_cff2
!>
            ad_cff2=ad_cff2+ad_v(i,j,k,nnew)
!>          tl_rvfrc(i,j)=tl_rvfrc(i,j)+tl_cff1
!>
            ad_cff1=ad_cff1+ad_rvfrc(i,j)
!>          tl_cff2=dt(ng)*cff*tl_cff1
!>
            ad_cff1=ad_cff1+dt(ng)*cff*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=0.5_r8*((pn(i,j-1)+pn(i,j))*                        &
!>   &                      (tl_VFx(i+1,j)-tl_VFx(i,j  ))-              &
!>   &                      (pm(i,j-1)+pm(i,j))*                        &
!>   &                      (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
!>
            adfac=0.5_r8*ad_cff1
            adfac1=adfac*(pn(i,j-1)+pn(i,j))
            adfac2=adfac*(pm(i,j-1)+pm(i,j))
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-adfac1
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+adfac1
            ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac2
            ad_VFe(i,j  )=ad_VFe(i,j  )-adfac2
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
!>          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+tl_cff2
!>
            ad_cff2=ad_cff2+ad_u(i,j,k,nnew)
!>          tl_rufrc(i,j)=tl_rufrc(i,j)+tl_cff1
!>
            ad_cff1=ad_cff1+ad_rufrc(i,j)
!>          tl_cff2=dt(ng)*cff*tl_cff1
!>
            ad_cff1=ad_cff1+dt(ng)*cff*ad_cff2
            ad_cff2=0.0_r8
!>          cff1=0.5_r8*((pn(i-1,j)+pn(i,j))*                           &
!>   &                   (UFx(i,j  )-UFx(i-1,j))+                       &
!>   &                   (pm(i-1,j)+pm(i,j))*                           &
!>   &                   (UFe(i,j+1)-UFe(i  ,j)))
!>
!>
            adfac=0.5_r8*ad_cff1
            adfac1=adfac*(pn(i-1,j)+pn(i,j))
            adfac2=adfac*(pm(i-1,j)+pm(i,j))
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-adfac1
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac1
            ad_UFe(i,j  )=ad_UFe(i,j  )-adfac2
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+adfac2
            ad_cff1=0.0_r8
          END DO
        END DO
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
!>          tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
!>          tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
!>
            ad_cff=ad_cff+                                              &
     &             on_p(i,j)*on_p(i,j)*ad_VFx(i,j)+                     &
     &             om_p(i,j)*om_p(i,j)*ad_UFe(i,j)
            ad_VFx(i,j)=0.0_r8
            ad_UFe(i,j)=0.0_r8
#ifdef MASKING
!>          tl_cff=tl_cff*pmask(i,j)
!>
            ad_cff=ad_cff*pmask(i,j)
#endif
!>          tl_cff=visc2_p(i,j)*0.125_r8*                               &
!>   &             ((tl_Hz(i-1,j  ,k)+tl_Hz(i,j  ,k)+                   &
!>   &               tl_Hz(i-1,j-1,k)+tl_Hz(i,j-1,k))*                  &
!>   &              (pmon_p(i,j)*                                       &
!>   &               ((pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nrhs)-          &
!>   &                (pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nrhs))+         &
!>   &               pnom_p(i,j)*                                       &
!>   &               ((pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nrhs)-          &
!>   &                (pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nrhs)))+        &
!>   &              (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                         &
!>   &               Hz(i-1,j-1,k)+Hz(i,j-1,k))*                        &
!>   &              (pmon_p(i,j)*                                       &
!>   &               ((pn(i  ,j-1)+pn(i  ,j))*tl_v(i  ,j,k,nrhs)-       &
!>   &                (pn(i-1,j-1)+pn(i-1,j))*tl_v(i-1,j,k,nrhs))+      &
!>   &               pnom_p(i,j)*                                       &
!>   &               ((pm(i-1,j  )+pm(i,j  ))*tl_u(i,j  ,k,nrhs)-       &
!>   &                (pm(i-1,j-1)+pm(i,j-1))*tl_u(i,j-1,k,nrhs))))
!>

            adfac=visc2_p(i,j)*0.125_r8*ad_cff
            adfac1=adfac*                                               &
     &             (pmon_p(i,j)*                                        &
     &              ((pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nrhs)-           &
     &               (pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nrhs))+          &
     &              pnom_p(i,j)*                                        &
     &              ((pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nrhs)-           &
     &               (pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nrhs)))
            adfac2=adfac*(Hz(i-1,j  ,k)+Hz(i,j  ,k)+                    &
     &                    Hz(i-1,j-1,k)+Hz(i,j-1,k))
            adfac3=adfac2*pmon_p(i,j)
            adfac4=adfac2*pnom_p(i,j)
            ad_Hz(i-1,j-1,k)=ad_Hz(i-1,j-1,k)+adfac1
            ad_Hz(i  ,j-1,k)=ad_Hz(i  ,j-1,k)+adfac1
            ad_Hz(i-1,j  ,k)=ad_Hz(i-1,j  ,k)+adfac1
            ad_Hz(i  ,j  ,k)=ad_Hz(i  ,j  ,k)+adfac1
            ad_v(i-1,j,k,nrhs)=ad_v(i-1,j,k,nrhs)-                      &
     &                         (pn(i-1,j-1)+pn(i-1,j))*adfac3
            ad_v(i  ,j,k,nrhs)=ad_v(i  ,j,k,nrhs)+                      &
     &                         (pn(i  ,j-1)+pn(i  ,j))*adfac3
            ad_u(i,j-1,k,nrhs)=ad_u(i,j-1,k,nrhs)-                      &
     &                         (pm(i-1,j-1)+pm(i,j-1))*adfac4
            ad_u(i,j  ,k,nrhs)=ad_u(i,j  ,k,nrhs)+                      &
     &                         (pm(i-1,j  )+pm(i,j  ))*adfac4
            ad_cff=0.0_r8
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
!>          tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
!>          tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
!>
            ad_cff=ad_cff+                                              &
     &             om_r(i,j)*om_r(i,j)*ad_VFe(i,j)+                     &
     &             on_r(i,j)*on_r(i,j)*ad_UFx(i,j)
            ad_VFe(i,j)=0.0_r8
            ad_UFx(i,j)=0.0_r8
!>          tl_cff=visc2_r(i,j)*0.5_r8*                                 &
!>   &             (tl_Hz(i,j,k)*                                       &
!>   &              (pmon_r(i,j)*                                       &
!>   &               ((pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nrhs)-            &
!>   &                (pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nrhs))-           &
!>   &               pnom_r(i,j)*                                       &
!>   &               ((pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nrhs)-            &
!>   &                (pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nrhs)))+          &
!>   &              Hz(i,j,k)*                                          &
!>   &              (pmon_r(i,j)*                                       &
!>   &               ((pn(i  ,j)+pn(i+1,j))*tl_u(i+1,j,k,nrhs)-         &
!>   &                (pn(i-1,j)+pn(i  ,j))*tl_u(i  ,j,k,nrhs))-        &
!>   &               pnom_r(i,j)*                                       &
!>   &               ((pm(i,j  )+pm(i,j+1))*tl_v(i,j+1,k,nrhs)-         &
!>   &                (pm(i,j-1)+pm(i,j  ))*tl_v(i,j  ,k,nrhs))))
!>
            adfac=visc2_r(i,j)*0.5_r8*ad_cff
            adfac1=adfac*Hz(i,j,k)
            adfac2=adfac1*pmon_r(i,j)
            adfac3=adfac1*pnom_r(i,j)
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+                                  &
     &                   (pmon_r(i,j)*                                  &
     &                    ((pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nrhs)-       &
     &                     (pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nrhs))-      &
     &                    pnom_r(i,j)*                                  &
     &                   ((pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nrhs)-        &
     &                    (pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nrhs)))*      &
     &                   adfac
            ad_u(i  ,j,k,nrhs)=ad_u(i  ,j,k,nrhs)-                      &
     &                         (pn(i-1,j)+pn(i  ,j))*adfac2
            ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+                      &
     &                         (pn(i  ,j)+pn(i+1,j))*adfac2
            ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)+                      &
     &                         (pm(i,j-1)+pm(i,j  ))*adfac3
            ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)-                      &
     &                         (pm(i,j  )+pm(i,j+1))*adfac3
            ad_cff=0.0_r8
          END DO
        END DO
      END DO K_LOOP
      RETURN
      END SUBROUTINE ad_uv3dmix2_tile
