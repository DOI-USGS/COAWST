      SUBROUTINE rp_uv3dmix2 (ng, tile)
!
!svn $Id: rp_uv3dmix2_s.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine computes representers tangent linear harmonic mixing   !
!  of momentum,  along  constant  S-surfaces,  from  the  horizontal   !
!  divergence of the stress tensor. A transverse isotropy is assumed   !
!  so the stress tensor is splitted  into  vertical  and  horizontal   !
!  subtensors.                                                         !
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
!***********************************************************************
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
      CALL wclock_on (ng, iRPM, 30)
#endif
      CALL rp_uv3dmix2_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs(ng), nnew(ng),                        &
#ifdef MASKING
     &                       GRID(ng) % pmask,                          &
#endif
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % tl_Hz,                          &
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
     &                       COUPLING(ng) % tl_rufrc,                   &
     &                       COUPLING(ng) % tl_rvfrc,                   &
     &                       OCEAN(ng) % tl_u,                          &
     &                       OCEAN(ng) % tl_v)
#ifdef PROFILE
      CALL wclock_off (ng, iRPM, 30)
#endif
      RETURN
      END SUBROUTINE rp_uv3dmix2

!
!***********************************************************************
      SUBROUTINE rp_uv3dmix2_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nrhs, nnew,                          &
#ifdef MASKING
     &                             pmask,                               &
#endif
     &                             Hz, tl_Hz,                           &
     &                             om_p, om_r, on_p, on_r,              &
     &                             pm, pmon_p, pmon_r,                  &
     &                             pn, pnom_p, pnom_r,                  &
     &                             visc2_p, visc2_r,                    &
!!#ifdef DIAGNOSTICS_UV
!!   &                             DiaRUfrc, DiaRVfrc,                  &
!!   &                             DiaU3wrk, DiaV3wrk,                  &
!!#endif
     &                             u, v,                                &
     &                             tl_rufrc, tl_rvfrc,                  &
     &                             tl_u, tl_v)
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
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
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

      real(r8), intent(inout) :: tl_rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: tl_rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
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

      real(r8), intent(inout) :: tl_rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: cff, tl_cff, tl_cff1, tl_cff2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_VFx

#include "set_bounds.h"
!

!-----------------------------------------------------------------------
!  Compute horizontal harmonic viscosity along constant S-surfaces.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=visc2_r(i,j)*Hz(i,j,k)*0.5_r8*                          &
     &          (pmon_r(i,j)*                                           &
     &           ((pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nrhs)-                &
     &            (pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nrhs))-               &
     &           pnom_r(i,j)*                                           &
     &           ((pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nrhs)-                &
     &            (pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nrhs)))
            tl_cff=visc2_r(i,j)*0.5_r8*                                 &
     &             (tl_Hz(i,j,k)*                                       &
     &              (pmon_r(i,j)*                                       &
     &               ((pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nrhs)-            &
     &                (pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nrhs))-           &
     &               pnom_r(i,j)*                                       &
     &               ((pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nrhs)-            &
     &                (pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nrhs)))+          &
     &              Hz(i,j,k)*                                          &
     &              (pmon_r(i,j)*                                       &
     &               ((pn(i  ,j)+pn(i+1,j))*tl_u(i+1,j,k,nrhs)-         &
     &                (pn(i-1,j)+pn(i  ,j))*tl_u(i  ,j,k,nrhs))-        &
     &               pnom_r(i,j)*                                       &
     &               ((pm(i,j  )+pm(i,j+1))*tl_v(i,j+1,k,nrhs)-         &
     &                (pm(i,j-1)+pm(i,j  ))*tl_v(i,j  ,k,nrhs))))-      &
#ifdef TL_IOMS
     &             cff
#endif
!>          UFx(i,j)=on_r(i,j)*on_r(i,j)*cff
!>
            tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*tl_cff
!>          VFe(i,j)=om_r(i,j)*om_r(i,j)*cff
!>
            tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*tl_cff
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend+1
            cff=visc2_p(i,j)*0.125_r8*(Hz(i-1,j  ,k)+Hz(i,j  ,k)+       &
     &                                 Hz(i-1,j-1,k)+Hz(i,j-1,k))*      &
     &          (pmon_p(i,j)*                                           &
     &           ((pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nrhs)-              &
     &            (pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nrhs))+             &
     &           pnom_p(i,j)*                                           &
     &           ((pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nrhs)-              &
     &            (pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nrhs)))
            tl_cff=visc2_p(i,j)*0.125_r8*                               &
     &             ((tl_Hz(i-1,j  ,k)+tl_Hz(i,j  ,k)+                   &
     &               tl_Hz(i-1,j-1,k)+tl_Hz(i,j-1,k))*                  &
     &              (pmon_p(i,j)*                                       &
     &               ((pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nrhs)-          &
     &                (pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nrhs))+         &
     &               pnom_p(i,j)*                                       &
     &               ((pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nrhs)-          &
     &                (pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nrhs)))+        &
     &              (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                         &
     &               Hz(i-1,j-1,k)+Hz(i,j-1,k))*                        &
     &              (pmon_p(i,j)*                                       &
     &               ((pn(i  ,j-1)+pn(i  ,j))*tl_v(i  ,j,k,nrhs)-       &
     &                (pn(i-1,j-1)+pn(i-1,j))*tl_v(i-1,j,k,nrhs))+      &
     &               pnom_p(i,j)*                                       &
     &               ((pm(i-1,j  )+pm(i,j  ))*tl_u(i,j  ,k,nrhs)-       &
     &                (pm(i-1,j-1)+pm(i,j-1))*tl_u(i,j-1,k,nrhs))))-    &
#ifdef TL_IOMS
     &             cff
#endif
#ifdef MASKING
!>          cff=cff*pmask(i,j)
!>
            tl_cff=tl_cff*pmask(i,j)
#endif
!>          UFe(i,j)=om_p(i,j)*om_p(i,j)*cff
!>
            tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*tl_cff
!>          VFx(i,j)=on_p(i,j)*on_p(i,j)*cff
!>
            tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*tl_cff
          END DO
        END DO
!
! Time-step harmonic, S-surfaces viscosity term.  Notice that momentum
! at this stage is HzU and HzV and has m2/s units.  Add contribution for
! barotropic forcing terms.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
!>          cff1=0.5_r8*((pn(i-1,j)+pn(i,j))*                           &
!>   &                   (UFx(i,j  )-UFx(i-1,j))+                       &
!>   &                   (pm(i-1,j)+pm(i,j))*                           &
!>   &                   (UFe(i,j+1)-UFe(i  ,j)))
!>
            tl_cff1=0.5_r8*((pn(i-1,j)+pn(i,j))*                        &
     &                      (tl_UFx(i,j  )-tl_UFx(i-1,j))+              &
     &                      (pm(i-1,j)+pm(i,j))*                        &
     &                      (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
!>          cff2=dt(ng)*cff*cff1
!>
            tl_cff2=dt(ng)*cff*tl_cff1
!>          rufrc(i,j)=rufrc(i,j)+cff1
!>
            tl_rufrc(i,j)=tl_rufrc(i,j)+tl_cff1

!>          u(i,j,k,nnew)=u(i,j,k,nnew)+cff2
!!#ifdef DIAGNOSTICS_UV
!!          DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)+cff1
!!          DiaU3wrk(i,j,k,M3hvis)=cff2
!!#endif
!>
            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+tl_cff2
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!>          cff1=0.5_r8*((pn(i,j-1)+pn(i,j))*                           &
!>   &                   (VFx(i+1,j)-VFx(i,j  ))-                       &
!>   &                   (pm(i,j-1)+pm(i,j))*                           &
!>   &                   (VFe(i  ,j)-VFe(i,j-1)))
!>
            tl_cff1=0.5_r8*((pn(i,j-1)+pn(i,j))*                        &
     &                      (tl_VFx(i+1,j)-tl_VFx(i,j  ))-              &
     &                      (pm(i,j-1)+pm(i,j))*                        &
     &                      (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
!>          cff2=dt(ng)*cff*cff1
!>
            tl_cff2=dt(ng)*cff*tl_cff1
!>          rvfrc(i,j)=rvfrc(i,j)+cff1
!>
            tl_rvfrc(i,j)=tl_rvfrc(i,j)+tl_cff1
!>          v(i,j,k,nnew)=v(i,j,k,nnew)+cff2
!!#ifdef DIAGNOSTICS_UV
!!          DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)+cff1
!!          DiaV3wrk(i,j,k,M3hvis)=cff2
!!#endif
!>
            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+tl_cff2
          END DO
        END DO
      END DO K_LOOP
      RETURN
      END SUBROUTINE rp_uv3dmix2_tile
