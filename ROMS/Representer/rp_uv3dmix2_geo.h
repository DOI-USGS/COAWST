      SUBROUTINE rp_uv3dmix2 (ng, tile)
!
!svn $Id: rp_uv3dmix2_geo.h 795 2016-05-11 01:42:43Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine computes representers tangent linear harmonic mixing   !
!  of momentum,  rotated along  geopotentials,  from the  horizontal   !
!  divergence of the stress tensor. A transverse isotropy is assumed   !
!  so the  stress tensor  is  split  into  vertical  and  horizontal   !
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
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
#ifdef DIAGNOSTICS_UV
!!    USE mod_diags
#endif
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
      CALL wclock_on (ng, iRPM, 31)
#endif
      CALL rp_uv3dmix2_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs(ng), nnew(ng),                        &
#ifdef MASKING
     &                       GRID(ng) % pmask,                          &
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
#endif
     &                       GRID(ng) % om_p,                           &
     &                       GRID(ng) % om_r,                           &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % om_v,                           &
     &                       GRID(ng) % on_p,                           &
     &                       GRID(ng) % on_r,                           &
     &                       GRID(ng) % on_u,                           &
     &                       GRID(ng) % on_v,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % tl_Hz,                          &
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % tl_z_r,                         &
#ifdef VISC_3DCOEF
     &                       MIXING(ng) % visc3d_r,                     &
     &                       MIXING(ng) % tl_visc3d_r,                  &
#else
     &                       MIXING(ng) % visc2_p,                      &
     &                       MIXING(ng) % visc2_r,                      &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                       DIAGS(ng) % DiaRUfrc,                      &
!!   &                       DIAGS(ng) % DiaRVfrc,                      &
!!   &                       DIAGS(ng) % DiaU3wrk,                      &
!!   &                       DIAGS(ng) % DiaV3wrk,                      &
#endif
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % tl_u,                          &
     &                       OCEAN(ng) % tl_v,                          &
     &                       COUPLING(ng) % tl_rufrc,                   &
     &                       COUPLING(ng) % tl_rvfrc)
#ifdef PROFILE
      CALL wclock_off (ng, iRPM, 31)
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
     &                             pmask, rmask, umask, vmask,          &
#endif
     &                             om_p, om_r, om_u, om_v,              &
     &                             on_p, on_r, on_u, on_v,              &
     &                             pm, pn,                              &
     &                             Hz, tl_Hz,                           &
     &                             z_r, tl_z_r,                         &
#ifdef VISC_3DCOEF
     &                             visc3d_r, tl_visc3d_r,               &
#else
     &                             visc2_p, visc2_r,                    &
#endif
#ifdef DIAGNOSTICS_UV
!!   &                             DiaRUfrc, DiaRVfrc,                  &
!!   &                             DiaU3wrk, DiaV3wrk,                  &
#endif
     &                             u, v,                                &
     &                             tl_u, tl_v,                          &
     &                             tl_rufrc, tl_rvfrc)
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
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: om_p(LBi:,LBj:)
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_p(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
# ifdef VISC_3DCOEF
      real(r8), intent(in) :: visc3d_r(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: visc2_p(LBi:,LBj:)
      real(r8), intent(in) :: visc2_r(LBi:,LBj:)
# endif
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)

# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: tl_rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: tl_rvfrc(LBi:,LBj:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
# ifdef VISC_3DCOEF
      real(r8), intent(inout) :: tl_visc3d_r(LBi:,LBj:,:)
# endif

#else
# ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: om_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
# ifdef VISC_3DCOEF
      real(r8), intent(in) :: visc3d_r(LBi:UBi,LBj:UBj,N(ng))
# else
      real(r8), intent(in) :: visc2_p(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: visc2_r(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)

# ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
!!    real(r8), intent(inout) :: DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
# endif
      real(r8), intent(inout) :: tl_rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_rvfrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
# ifdef VISC_3DCOEF
      real(r8), intent(inout) :: tl_visc3d_r(LBi:UBi,LBj:UBj,N(ng))
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k, k1, k2

      real(r8) :: cff, fac1, fac2, pm_p, pn_p
      real(r8) :: cff1, cff2, cff3, cff4
      real(r8) :: cff5, cff6, cff7, cff8
      real(r8) :: dmUdz, dnUdz, dmVdz, dnVdz
#ifdef VISC_3DCOEF
      real(r8) :: visc_p
      real(r8) :: tl_fac1, tl_fac2, tl_visc_p
#endif
      real(r8) :: tl_cff
      real(r8) :: tl_cff1, tl_cff2, tl_cff3, tl_cff4
      real(r8) :: tl_cff5, tl_cff6, tl_cff7, tl_cff8
      real(r8) :: tl_dmUdz, tl_dnUdz, tl_dmVdz, tl_dnVdz

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_VFx

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dmUde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dmVde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dnUdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dnVdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dUdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dVdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx_r

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_UFse
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_UFsx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_VFse
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_VFsx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dmUde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dmVde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dnUdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dnVdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dUdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dVdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dZde_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dZde_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dZdx_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dZdx_r

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic viscosity along geopotential surfaces.
!-----------------------------------------------------------------------
!
!  Compute horizontal and vertical gradients.  Notice the recursive
!  blocking sequence.  The vertical placement of the gradients is:
!
!    dZdx_r, dZde_r, dnUdx, dmVde(:,:,k1) k      rho-points
!    dZdx_r, dZde_r, dnUdx, dmVde(:,:,k2) k+1    rho-points
!    dZdx_p, dZde_p, dnVdx, dmUde(:,:,k1) k      psi-points
!    dZdx_p, dZde_p, dnVdx, dmUde(:,:,k2) k+1    psi-points
!                UFse, UFsx, dUdz(:,:,k1) k-1/2  WU-points
!                UFse, UFsx, dUdz(:,:,k2) k+1/2  WU-points
!                VFse, VFsx, dVdz(:,:,k1) k-1/2  WV-points
!                VFse, VFsx, dVdz(:,:,k2) k+1/2  WV-points
!
      k2=1
      K_LOOP : DO k=0,N(ng)
        k1=k2
        k2=3-k1
        IF (k.lt.N(ng)) THEN
!
!  Compute slopes (nondimensional) at RHO- and PSI-points.
!
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              cff=0.5_r8*(pm(i-1,j)+pm(i,j))
#ifdef MASKING
              cff=cff*umask(i,j)
#endif
              UFx(i,j)=cff*(z_r(i  ,j,k+1)-                             &
     &                      z_r(i-1,j,k+1))
              tl_UFx(i,j)=cff*(tl_z_r(i  ,j,k+1)-                       &
     &                         tl_z_r(i-1,j,k+1))
            END DO
          END DO
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=0.5_r8*(pn(i,j-1)+pn(i,j))
#ifdef MASKING
              cff=cff*vmask(i,j)
#endif
              VFe(i,j)=cff*(z_r(i,j  ,k+1)-                             &
     &                      z_r(i,j-1,k+1))
              tl_VFe(i,j)=cff*(tl_z_r(i,j  ,k+1)-                       &
     &                         tl_z_r(i,j-1,k+1))
            END DO
          END DO
!
          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
              dZdx_p(i,j,k2)=0.5_r8*(UFx(i,j-1)+                        &
     &                               UFx(i,j  ))
              tl_dZdx_p(i,j,k2)=0.5_r8*(tl_UFx(i,j-1)+                  &
     &                                  tl_UFx(i,j  ))
              dZde_p(i,j,k2)=0.5_r8*(VFe(i-1,j)+                        &
     &                               VFe(i  ,j))
              tl_dZde_p(i,j,k2)=0.5_r8*(tl_VFe(i-1,j)+                  &
     &                                  tl_VFe(i  ,j))
            END DO
          END DO
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              dZdx_r(i,j,k2)=0.5_r8*(UFx(i  ,j)+                        &
     &                               UFx(i+1,j))
              tl_dZdx_r(i,j,k2)=0.5_r8*(tl_UFx(i  ,j)+                  &
     &                                  tl_UFx(i+1,j))
              dZde_r(i,j,k2)=0.5_r8*(VFe(i,j  )+                        &
     &                               VFe(i,j+1))
              tl_dZde_r(i,j,k2)=0.5_r8*(tl_VFe(i,j  )+                  &
     &                                  tl_VFe(i,j+1))
            END DO
          END DO
!
!  Compute momentum horizontal (1/m/s) and vertical (1/s) gradients.
!
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff=0.5_r8*pm(i,j)
#ifdef MASKING
              cff=cff*rmask(i,j)
#endif
              dnUdx(i,j,k2)=cff*((pn(i  ,j)+pn(i+1,j))*                 &
     &                           u(i+1,j,k+1,nrhs)-                     &
     &                           (pn(i-1,j)+pn(i  ,j))*                 &
     &                           u(i  ,j,k+1,nrhs))
              tl_dnUdx(i,j,k2)=cff*((pn(i  ,j)+pn(i+1,j))*              &
     &                              tl_u(i+1,j,k+1,nrhs)-               &
     &                              (pn(i-1,j)+pn(i  ,j))*              &
     &                              tl_u(i  ,j,k+1,nrhs))
            END DO
          END DO

          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
              cff=0.125_r8*(pn(i-1,j  )+pn(i,j  )+                      &
     &                      pn(i-1,j-1)+pn(i,j-1))
#ifdef MASKING
              cff=cff*pmask(i,j)
#endif
              dmUde(i,j,k2)=cff*((pm(i-1,j  )+pm(i,j  ))*               &
     &                           u(i,j  ,k+1,nrhs)-                     &
     &                           (pm(i-1,j-1)+pm(i,j-1))*               &
     &                           u(i,j-1,k+1,nrhs))
              tl_dmUde(i,j,k2)=cff*((pm(i-1,j  )+pm(i,j  ))*            &
     &                              tl_u(i,j  ,k+1,nrhs)-               &
     &                              (pm(i-1,j-1)+pm(i,j-1))*            &
     &                              tl_u(i,j-1,k+1,nrhs))
            END DO
          END DO

          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
              cff=0.125_r8*(pm(i-1,j  )+pm(i,j  )+                      &
     &                      pm(i-1,j-1)+pm(i,j-1))
#ifdef MASKING
              cff=cff*pmask(i,j)
#endif
              dnVdx(i,j,k2)=cff*((pn(i  ,j-1)+pn(i  ,j))*               &
     &                           v(i  ,j,k+1,nrhs)-                     &
     &                           (pn(i-1,j-1)+pn(i-1,j))*               &
     &                           v(i-1,j,k+1,nrhs))
              tl_dnVdx(i,j,k2)=cff*((pn(i  ,j-1)+pn(i  ,j))*            &
     &                              tl_v(i  ,j,k+1,nrhs)-               &
     &                              (pn(i-1,j-1)+pn(i-1,j))*            &
     &                              tl_v(i-1,j,k+1,nrhs))
            END DO
          END DO

          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff=0.5_r8*pn(i,j)
#ifdef MASKING
              cff=cff*rmask(i,j)
#endif
              dmVde(i,j,k2)=cff*((pm(i,j  )+pm(i,j+1))*                 &
     &                           v(i,j+1,k+1,nrhs)-                     &
     &                           (pm(i,j-1)+pm(i,j  ))*                 &
     &                           v(i,j  ,k+1,nrhs))
              tl_dmVde(i,j,k2)=cff*((pm(i,j  )+pm(i,j+1))*              &
     &                              tl_v(i,j+1,k+1,nrhs)-               &
     &                              (pm(i,j-1)+pm(i,j  ))*              &
     &                              tl_v(i,j  ,k+1,nrhs))
            END DO
          END DO
        END IF

        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              dUdz(i,j,k2)=0.0_r8
              tl_dUdz(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              dVdz(i,j,k2)=0.0_r8
              tl_dVdz(i,j,k2)=0.0_r8
            END DO
          END DO

          DO j=Jstr,Jend
            DO i=IstrU,Iend
!>            UFsx(i,j,k2)=0.0_r8
!>
              tl_UFsx(i,j,k2)=0.0_r8
!>            UFse(i,j,k2)=0.0_r8
!>
              tl_UFse(i,j,k2)=0.0_r8
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
!>            VFsx(i,j,k2)=0.0_r8
!>
              tl_VFsx(i,j,k2)=0.0_r8
!>            VFse(i,j,k2)=0.0_r8
!>
              tl_VFse(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i-1,j,k+1)-z_r(i-1,j,k)+          &
     &                            z_r(i  ,j,k+1)-z_r(i  ,j,k)))
              tl_cff=-cff*cff*(0.5_r8*(tl_z_r(i-1,j,k+1)-               &
     &                                 tl_z_r(i-1,j,k  )+               &
     &                                 tl_z_r(i  ,j,k+1)-               &
     &                                 tl_z_r(i  ,j,k  )))+             &
#ifdef TL_IOMS
     &               2.0_r8*cff
#endif
              dUdz(i,j,k2)=cff*(u(i,j,k+1,nrhs)-                        &
     &                          u(i,j,k  ,nrhs))
              tl_dUdz(i,j,k2)=tl_cff*(u(i,j,k+1,nrhs)-                  &
     &                                u(i,j,k  ,nrhs))+                 &
     &                        cff*(tl_u(i,j,k+1,nrhs)-                  &
     &                             tl_u(i,j,k  ,nrhs))-                 &
#ifdef TL_IOMS
     &                        dUdz(i,j,k2)
#endif
            END DO
          END DO

          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(0.5_r8*(z_r(i,j-1,k+1)-z_r(i,j-1,k)+          &
     &                            z_r(i,j  ,k+1)-z_r(i,j  ,k)))
              tl_cff=-cff*cff*(0.5_r8*(tl_z_r(i,j-1,k+1)-               &
     &                                 tl_z_r(i,j-1,k  )+               &
     &                                 tl_z_r(i,j  ,k+1)-               &
     &                                 tl_z_r(i,j  ,k  )))+             &
#ifdef TL_IOMS
     &               2.0_r8*cff
#endif
              dVdz(i,j,k2)=cff*(v(i,j,k+1,nrhs)-                        &
     &                          v(i,j,k  ,nrhs))
              tl_dVdz(i,j,k2)=tl_cff*(v(i,j,k+1,nrhs)-                  &
     &                                v(i,j,k  ,nrhs))+                 &
     &                        cff*(tl_v(i,j,k+1,nrhs)-                  &
     &                             tl_v(i,j,k  ,nrhs))-                 &
#ifdef TL_IOMS
     &                        dVdz(i,j,k2)
#endif
            END DO
          END DO
        END IF
!
!  Compute components of the rotated viscous flux (m5/s2) along
!  geopotential surfaces in the XI- and ETA-directions.
!
        IF (k.gt.0) THEN
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
              cff1=MIN(dZdx_r(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_r(i,j,k1),0.0_r8)
              cff3=MIN(dZde_r(i,j,k1),0.0_r8)
              cff4=MAX(dZde_r(i,j,k1),0.0_r8)
              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i,j,k1)))*            &
     &                tl_dZdx_r(i,j,k1)
              tl_cff2=(0.5_r8+SIGN(0.5_r8, dZdx_r(i,j,k1)))*            &
     &                tl_dZdx_r(i,j,k1)
              tl_cff3=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j,k1)))*            &
     &                tl_dZde_r(i,j,k1)
              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j,k1)))*            &
     &                tl_dZde_r(i,j,k1)
#ifdef VISC_3DCOEF
              cff=Hz(i,j,k)*                                            &
     &            (on_r(i,j)*(dnUdx(i,j,k1)-                            &
     &                        0.5_r8*pn(i,j)*                           &
     &                        (cff1*(dUdz(i  ,j,k1)+                    &
     &                               dUdz(i+1,j,k2))+                   &
     &                         cff2*(dUdz(i  ,j,k2)+                    &
     &                               dUdz(i+1,j,k1))))-                 &
     &             om_r(i,j)*(dmVde(i,j,k1)-                            &
     &                        0.5_r8*pm(i,j)*                           &
     &                        (cff3*(dVdz(i,j  ,k1)+                    &
     &                               dVdz(i,j+1,k2))+                   &
     &                         cff4*(dVdz(i,j  ,k2)+                    &
     &                               dVdz(i,j+1,k1)))))
#else
!>            cff=Hz(i,j,k)*                                            &
!>   &            (on_r(i,j)*(dnUdx(i,j,k1)-                            &
!>   &                        0.5_r8*pn(i,j)*                           &
!>   &                        (cff1*(dUdz(i  ,j,k1)+                    &
!>   &                               dUdz(i+1,j,k2))+                   &
!>   &                         cff2*(dUdz(i  ,j,k2)+                    &
!>   &                               dUdz(i+1,j,k1))))-                 &
!>   &             om_r(i,j)*(dmVde(i,j,k1)-                            &
!>   &                        0.5_r8*pm(i,j)*                           &
!>   &                        (cff3*(dVdz(i,j  ,k1)+                    &
!>   &                               dVdz(i,j+1,k2))+                   &
!>   &                         cff4*(dVdz(i,j  ,k2)+                    &
!>   &                               dVdz(i,j+1,k1)))))
!>
#endif
              tl_cff=tl_Hz(i,j,k)*                                      &
     &               (on_r(i,j)*(dnUdx(i,j,k1)-                         &
     &                           0.5_r8*pn(i,j)*                        &
     &                           (cff1*(dUdz(i  ,j,k1)+                 &
     &                                  dUdz(i+1,j,k2))+                &
     &                            cff2*(dUdz(i  ,j,k2)+                 &
     &                                  dUdz(i+1,j,k1))))-              &
     &                om_r(i,j)*(dmVde(i,j,k1)-                         &
     &                           0.5_r8*pm(i,j)*                        &
     &                           (cff3*(dVdz(i,j  ,k1)+                 &
     &                                  dVdz(i,j+1,k2))+                &
     &                            cff4*(dVdz(i,j  ,k2)+                 &
     &                                  dVdz(i,j+1,k1)))))+             &
     &                Hz(i,j,k)*                                        &
     &                (on_r(i,j)*(tl_dnUdx(i,j,k1)-                     &
     &                            0.5_r8*pn(i,j)*                       &
     &                            (tl_cff1*(dUdz(i  ,j,k1)+             &
     &                                      dUdz(i+1,j,k2))+            &
     &                             cff1*(tl_dUdz(i  ,j,k1)+             &
     &                                   tl_dUdz(i+1,j,k2))+            &
     &                             tl_cff2*(dUdz(i  ,j,k2)+             &
     &                                      dUdz(i+1,j,k1))+            &
     &                             cff2*(tl_dUdz(i  ,j,k2)+             &
     &                                   tl_dUdz(i+1,j,k1))))-          &
     &                 om_r(i,j)*(tl_dmVde(i,j,k1)-                     &
     &                            0.5_r8*pm(i,j)*                       &
     &                            (tl_cff3*(dVdz(i,j  ,k1)+             &
     &                                      dVdz(i,j+1,k2))+            &
     &                             cff3*(tl_dVdz(i,j  ,k1)+             &
     &                                   tl_dVdz(i,j+1,k2))+            &
     &                             tl_cff4*(dVdz(i,j  ,k2)+             &
     &                                      dVdz(i,j+1,k1))+            &
     &                             cff4*(tl_dVdz(i,j  ,k2)+             &
     &                                   tl_dVdz(i,j+1,k1)))))-         &
#ifdef TL_IOMS
     &                Hz(i,j,k)*                                        &
     &                (on_r(i,j)*(dnUdx(i,j,k1)-                        &
     &                            pn(i,j)*                              &
     &                            (cff1*(dUdz(i  ,j,k1)+                &
     &                                   dUdz(i+1,j,k2))+               &
     &                             cff2*(dUdz(i  ,j,k2)+                &
     &                                   dUdz(i+1,j,k1))))-             &
     &                 om_r(i,j)*(dmVde(i,j,k1)-                        &
     &                            pm(i,j)*                              &
     &                            (cff3*(dVdz(i,j  ,k1)+                &
     &                                   dVdz(i,j+1,k2))+               &
     &                             cff4*(dVdz(i,j  ,k2)+                &
     &                                   dVdz(i,j+1,k1)))))
#endif
#ifdef MASKING
# ifdef VISC_3DCOEF
              cff=cff*rmask(i,j)
# else
!>            cff=cff*rmask(i,j)
!>
# endif
              tl_cff=tl_cff*rmask(i,j)
#endif
#ifdef VISC_3DCOEF
              UFx(i,j)=on_r(i,j)*on_r(i,j)*visc3d_r(i,j,k)*cff
              tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*                          &
     &                    (tl_visc3d_r(i,j,k)*cff+                      &
     &                     visc3d_r(i,j,k)*tl_cff)-                     &
# ifdef TL_IOMS
     &                    UFx(i,j)
# endif
              VFe(i,j)=om_r(i,j)*om_r(i,j)*visc3d_r(i,j,k)*cff
              tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*                          &
     &                    (tl_visc3d_r(i,j,k)*cff+                      &
     &                     visc3d_r(i,j,k)*tl_cff)-                     &
# ifdef TL_IOMS
     &                    VFe(i,j)
# endif
#else
!>            UFx(i,j)=on_r(i,j)*on_r(i,j)*visc2_r(i,j)*cff
!>
              tl_UFx(i,j)=on_r(i,j)*on_r(i,j)*visc2_r(i,j)*tl_cff
!>            VFe(i,j)=om_r(i,j)*om_r(i,j)*visc2_r(i,j)*cff
!>
              tl_VFe(i,j)=om_r(i,j)*om_r(i,j)*visc2_r(i,j)*tl_cff
#endif
            END DO
          END DO

          DO j=Jstr,Jend+1
            DO i=Istr,Iend+1
              pm_p=0.25_r8*(pm(i-1,j-1)+pm(i-1,j)+                      &
     &                      pm(i  ,j-1)+pm(i  ,j))
              pn_p=0.25_r8*(pn(i-1,j-1)+pn(i-1,j)+                      &
     &                      pn(i  ,j-1)+pn(i  ,j))
              cff1=MIN(dZdx_p(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_p(i,j,k1),0.0_r8)
              cff3=MIN(dZde_p(i,j,k1),0.0_r8)
              cff4=MAX(dZde_p(i,j,k1),0.0_r8)
              tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i,j,k1)))*            &
     &                tl_dZdx_p(i,j,k1)
              tl_cff2=(0.5_r8+SIGN(0.5_r8, dZdx_p(i,j,k1)))*            &
     &                tl_dZdx_p(i,j,k1)
              tl_cff3=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j,k1)))*            &
     &                tl_dZde_p(i,j,k1)
              tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j,k1)))*            &
     &                tl_dZde_p(i,j,k1)
#ifdef VISC_3DCOEF
              cff=0.25_r8*                                              &
     &            (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                           &
     &             Hz(i-1,j-1,k)+Hz(i,j-1,k))*                          &
     &            (on_p(i,j)*(dnVdx(i,j,k1)-                            &
     &                        0.5_r8*pn_p*                              &
     &                        (cff1*(dVdz(i-1,j,k1)+                    &
     &                               dVdz(i  ,j,k2))+                   &
     &                         cff2*(dVdz(i-1,j,k2)+                    &
     &                               dVdz(i  ,j,k1))))+                 &
     &             om_p(i,j)*(dmUde(i,j,k1)-                            &
     &                        0.5_r8*pm_p*                              &
     &                        (cff3*(dUdz(i,j-1,k1)+                    &
     &                               dUdz(i,j  ,k2))+                   &
     &                         cff4*(dUdz(i,j-1,k2)+                    &
     &                               dUdz(i,j  ,k1)))))
#else
!>            cff=0.25_r8*                                              &
!>   &            (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                           &
!>   &             Hz(i-1,j-1,k)+Hz(i,j-1,k))*                          &
!>   &            (on_p(i,j)*(dnVdx(i,j,k1)-                            &
!>   &                        0.5_r8*pn_p*                              &
!>   &                        (cff1*(dVdz(i-1,j,k1)+                    &
!>   &                               dVdz(i  ,j,k2))+                   &
!>   &                         cff2*(dVdz(i-1,j,k2)+                    &
!>   &                               dVdz(i  ,j,k1))))+                 &
!>   &             om_p(i,j)*(dmUde(i,j,k1)-                            &
!>   &                        0.5_r8*pm_p*                              &
!>   &                        (cff3*(dUdz(i,j-1,k1)+                    &
!>   &                               dUdz(i,j  ,k2))+                   &
!>   &                         cff4*(dUdz(i,j-1,k2)+                    &
!>   &                               dUdz(i,j  ,k1)))))
!>
#endif
              tl_cff=0.25_r8*                                           &
     &               ((tl_Hz(i-1,j  ,k)+tl_Hz(i,j  ,k)+                 &
     &                 tl_Hz(i-1,j-1,k)+tl_Hz(i,j-1,k))*                &
     &                (on_p(i,j)*(dnVdx(i,j,k1)-                        &
     &                            0.5_r8*pn_p*                          &
     &                            (cff1*(dVdz(i-1,j,k1)+                &
     &                                   dVdz(i  ,j,k2))+               &
     &                             cff2*(dVdz(i-1,j,k2)+                &
     &                                   dVdz(i  ,j,k1))))+             &
     &                 om_p(i,j)*(dmUde(i,j,k1)-                        &
     &                            0.5_r8*pm_p*                          &
     &                            (cff3*(dUdz(i,j-1,k1)+                &
     &                                   dUdz(i,j  ,k2))+               &
     &                             cff4*(dUdz(i,j-1,k2)+                &
     &                                   dUdz(i,j  ,k1)))))+            &
     &                (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                       &
     &                 Hz(i-1,j-1,k)+Hz(i,j-1,k))*                      &
     &                (on_p(i,j)*(tl_dnVdx(i,j,k1)-                     &
     &                            0.5_r8*pn_p*                          &
     &                            (tl_cff1*(dVdz(i-1,j,k1)+             &
     &                                      dVdz(i  ,j,k2))+            &
     &                             cff1*(tl_dVdz(i-1,j,k1)+             &
     &                                   tl_dVdz(i  ,j,k2))+            &
     &                             tl_cff2*(dVdz(i-1,j,k2)+             &
     &                                      dVdz(i  ,j,k1))+            &
     &                             cff2*(tl_dVdz(i-1,j,k2)+             &
     &                                   tl_dVdz(i  ,j,k1))))+          &
     &                 om_p(i,j)*(tl_dmUde(i,j,k1)-                     &
     &                            0.5_r8*pm_p*                          &
     &                            (tl_cff3*(dUdz(i,j-1,k1)+             &
     &                                      dUdz(i,j  ,k2))+            &
     &                             cff3*(tl_dUdz(i,j-1,k1)+             &
     &                                   tl_dUdz(i,j  ,k2))+            &
     &                             tl_cff4*(dUdz(i,j-1,k2)+             &
     &                                      dUdz(i,j  ,k1))+            &
     &                             cff4*(tl_dUdz(i,j-1,k2)+             &
     &                                   tl_dUdz(i,j  ,k1))))))-        &
#ifdef TL_IOMS
     &               0.25_r8*                                           &
     &               (Hz(i-1,j  ,k)+Hz(i,j  ,k)+                        &
     &                Hz(i-1,j-1,k)+Hz(i,j-1,k))*                       &
     &               (on_p(i,j)*(dnVdx(i,j,k1)-                         &
     &                           pn_p*                                  &
     &                           (cff1*(dVdz(i-1,j,k1)+                 &
     &                                  dVdz(i  ,j,k2))+                &
     &                            cff2*(dVdz(i-1,j,k2)+                 &
     &                                  dVdz(i  ,j,k1))))+              &
     &                om_p(i,j)*(dmUde(i,j,k1)-                         &
     &                           pm_p*                                  &
     &                           (cff3*(dUdz(i,j-1,k1)+                 &
     &                                  dUdz(i,j  ,k2))+                &
     &                            cff4*(dUdz(i,j-1,k2)+                 &
     &                                  dUdz(i,j  ,k1)))))
#endif
#ifdef MASKING
# ifdef VISC_3DCOEF
              cff=cff*pmask(i,j)
# else
!>            cff=cff*pmask(i,j)
!>
# endif
              tl_cff=tl_cff*pmask(i,j)
#endif
#ifdef VISC_3DCOEF
              visc_p=0.25_r8*                                           &
     &               (visc3d_r(i-1,j-1,k)+visc3d_r(i-1,j,k)+            &
     &                visc3d_r(i  ,j-1,k)+visc3d_r(i  ,j,k))
              tl_visc_p=0.25_r8*                                        &
     &                  (tl_visc3d_r(i-1,j-1,k)+tl_visc3d_r(i-1,j,k)+   &
     &                   tl_visc3d_r(i  ,j-1,k)+tl_visc3d_r(i  ,j,k))
              UFe(i,j)=om_p(i,j)*om_p(i,j)*visc_p*cff
              tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*                          &
     &                    (tl_visc_p*cff+visc_p*tl_cff)-                &
# ifdef TL_IOMS
     &                    UFe(i,j)
# endif
              VFx(i,j)=on_p(i,j)*on_p(i,j)*visc_p*cff
              tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*                          &
     &                    (tl_visc_p*cff+visc_p*tl_cff)-                &
# ifdef TL_IOMS
     &                    VFx(i,j)
# endif
#else
!>            UFe(i,j)=om_p(i,j)*om_p(i,j)*visc2_p(i,j)*cff
!>
              tl_UFe(i,j)=om_p(i,j)*om_p(i,j)*visc2_p(i,j)*tl_cff
!>            VFx(i,j)=on_p(i,j)*on_p(i,j)*visc2_p(i,j)*cff
!>
              tl_VFx(i,j)=on_p(i,j)*on_p(i,j)*visc2_p(i,j)*tl_cff
#endif
            END DO
          END DO
!
!  Compute vertical flux (m2/s2) due to sloping terrain-following
!  surfaces.
!
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=IstrU,Iend
#ifdef VISC_3DCOEF
                cff=0.125_r8*                                           &
     &              (visc3d_r(i-1,j,k  )+visc3d_r(i,j,k  )+             &
     &               visc3d_r(i-1,j,k+1)+visc3d_r(i,j,k+1))
                tl_cff=0.125_r8*                                        &
     &                 (tl_visc3d_r(i-1,j,k  )+tl_visc3d_r(i,j,k  )+    &
     &                  tl_visc3d_r(i-1,j,k+1)+tl_visc3d_r(i,j,k+1))
                fac1=cff*on_u(i,j)
                tl_fac1=tl_cff*on_u(i,j)
                fac2=cff*om_u(i,j)
                tl_fac2=tl_cff*om_u(i,j)
#else
                cff=0.25_r8*(visc2_r(i-1,j)+visc2_r(i,j))
                fac1=cff*on_u(i,j)
                fac2=cff*om_u(i,j)
#endif
                cff=0.5_r8*(pn(i-1,j)+pn(i,j))
                dnUdz=cff*dUdz(i,j,k2)
                tl_dnUdz=cff*tl_dUdz(i,j,k2)
                dnVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))
                tl_dnVdz=cff*0.25_r8*(tl_dVdz(i-1,j+1,k2)+              &
     &                                tl_dVdz(i  ,j+1,k2)+              &
     &                                tl_dVdz(i-1,j  ,k2)+              &
     &                                tl_dVdz(i  ,j  ,k2))
                cff=0.5_r8*(pm(i-1,j)+pm(i,j))
                dmUdz=cff*dUdz(i,j,k2)
                tl_dmUdz=cff*tl_dUdz(i,j,k2)
                dmVdz=cff*0.25_r8*(dVdz(i-1,j+1,k2)+                    &
     &                             dVdz(i  ,j+1,k2)+                    &
     &                             dVdz(i-1,j  ,k2)+                    &
     &                             dVdz(i  ,j  ,k2))
                tl_dmVdz=cff*0.25_r8*(tl_dVdz(i-1,j+1,k2)+              &
     &                                tl_dVdz(i  ,j+1,k2)+              &
     &                                tl_dVdz(i-1,j  ,k2)+              &
     &                                tl_dVdz(i  ,j  ,k2))

                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*        &
     &                  tl_dZdx_r(i-1,j,k1)
                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*        &
     &                  tl_dZdx_r(i  ,j,k2)
                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*        &
     &                  tl_dZdx_r(i-1,j,k2)
                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*        &
     &                  tl_dZdx_r(i  ,j,k1)
!>              UFsx(i,j,k2)=fac1*                                      &
!>   &                       (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+        &
!>   &                        cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+        &
!>   &                        cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+        &
!>   &                        cff4*(cff4*dnUdz-dnUdx(i  ,j,k1)))
!>
                tl_UFsx(i,j,k2)=fac1*                                   &
     &                          (tl_cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+  &
     &                           tl_cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+  &
     &                           tl_cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+  &
     &                           tl_cff4*(cff4*dnUdz-dnUdx(i  ,j,k1))+  &
     &                           cff1*(tl_cff1*dnUdz+cff1*tl_dnUdz-     &
     &                                 tl_dnUdx(i-1,j,k1))+             &
     &                           cff2*(tl_cff2*dnUdz+cff2*tl_dnUdz-     &
     &                                 tl_dnUdx(i  ,j,k2))+             &
     &                           cff3*(tl_cff3*dnUdz+cff3*tl_dnUdz-     &
     &                                 tl_dnUdx(i-1,j,k2))+             &
     &                           cff4*(tl_cff4*dnUdz+cff4*tl_dnUdz-     &
     &                                 tl_dnUdx(i  ,j,k1)))-            &
#ifdef TL_IOMS
     &                          fac1*                                   &
     &                          (cff1*(2.0_r8*cff1*dnUdz-               &
     &                                 dnUdx(i-1,j,k1))+                &
     &                           cff2*(2.0_r8*cff2*dnUdz-               &
     &                                 dnUdx(i  ,j,k2))+                &
     &                           cff3*(2.0_r8*cff3*dnUdz-               &
     &                                 dnUdx(i-1,j,k2))+                &
     &                           cff4*(2.0_r8*cff4*dnUdz-               &
     &                                 dnUdx(i  ,j,k1)))
#endif
#ifdef VISC_3DCOEF
                tl_UFsx(i,j,k2)=tl_UFsx(i,j,k2)+                        &
     &                          tl_fac1*                                &
     &                          (cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))+     &
     &                           cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))+     &
     &                           cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))+     &
     &                           cff4*(cff4*dnUdz-dnUdx(i  ,j,k1)))
#endif

                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*        &
     &                  tl_dZde_p(i,j  ,k1)
                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*        &
     &                  tl_dZde_p(i,j+1,k2)
                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j  ,k2)))*        &
     &                  tl_dZde_p(i,j  ,k2)
                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j+1,k1)))*        &
     &                  tl_dZde_p(i,j+1,k1)
!>              UFse(i,j,k2)=fac2*                                      &
!>   &                       (cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+        &
!>   &                        cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+        &
!>   &                        cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+        &
!>   &                        cff4*(cff4*dmUdz-dmUde(i,j+1,k1)))
!>
                tl_UFse(i,j,k2)=fac2*                                   &
     &                          (tl_cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+  &
     &                           tl_cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+  &
     &                           tl_cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+  &
     &                           tl_cff4*(cff4*dmUdz-dmUde(i,j+1,k1))+  &
     &                           cff1*(tl_cff1*dmUdz+cff1*tl_dmUdz-     &
     &                                 tl_dmUde(i,j  ,k1))+             &
     &                           cff2*(tl_cff2*dmUdz+cff2*tl_dmUdz-     &
     &                                 tl_dmUde(i,j+1,k2))+             &
     &                           cff3*(tl_cff3*dmUdz+cff3*tl_dmUdz-     &
     &                                 tl_dmUde(i,j  ,k2))+             &
     &                           cff4*(tl_cff4*dmUdz+cff4*tl_dmUdz-     &
     &                                 tl_dmUde(i,j+1,k1)))-            &
#ifdef TL_IOMS
     &                          fac2*                                   &
     &                          (cff1*(2.0_r8*cff1*dmUdz-               &
     &                                 dmUde(i,j  ,k1))+                &
     &                           cff2*(2.0_r8*cff2*dmUdz-               &
     &                                 dmUde(i,j+1,k2))+                &
     &                           cff3*(2.0_r8*cff3*dmUdz-               &
     &                                 dmUde(i,j  ,k2))+                &
     &                           cff4*(2.0_r8*cff4*dmUdz-               &
     &                                 dmUde(i,j+1,k1)))
#endif
#ifdef VISC_3DCOEF
                tl_UFse(i,j,k2)=tl_UFse(i,j,k2)+                        &
     &                          tl_fac2*                                &
     &                          (cff1*(cff1*dmUdz-dmUde(i,j  ,k1))+     &
     &                           cff2*(cff2*dmUdz-dmUde(i,j+1,k2))+     &
     &                           cff3*(cff3*dmUdz-dmUde(i,j  ,k2))+     &
     &                           cff4*(cff4*dmUdz-dmUde(i,j+1,k1)))
#endif

                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                cff5=MIN(dZdx_p(i,j  ,k1),0.0_r8)
                cff6=MIN(dZdx_p(i,j+1,k2),0.0_r8)
                cff7=MAX(dZdx_p(i,j  ,k2),0.0_r8)
                cff8=MAX(dZdx_p(i,j+1,k1),0.0_r8)
                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j  ,k1)))*        &
     &                  tl_dZde_p(i,j  ,k1)
                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_p(i,j+1,k2)))*        &
     &                  tl_dZde_p(i,j+1,k2)
                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j  ,k2)))*        &
     &                  tl_dZde_p(i,j  ,k2)
                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_p(i,j+1,k1)))*        &
     &                  tl_dZde_p(i,j+1,k1)
                tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i,j  ,k1)))*        &
     &                  tl_dZdx_p(i,j  ,k1)
                tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i,j+1,k2)))*        &
     &                  tl_dZdx_p(i,j+1,k2)
                tl_cff7=(0.5_r8+SIGN(0.5_r8, dZdx_p(i,j  ,k2)))*        &
     &                  tl_dZdx_p(i,j  ,k2)
                tl_cff8=(0.5_r8+SIGN(0.5_r8, dZdx_p(i,j+1,k1)))*        &
     &                  tl_dZdx_p(i,j+1,k1)
!>              UFsx(i,j,k2)=UFsx(i,j,k2)+                              &
!>   &                       fac1*                                      &
!>   &                       (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+        &
!>   &                        cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+        &
!>   &                        cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+        &
!>   &                        cff4*(cff8*dnVdz-dnVdx(i,j+1,k1)))
!>
                tl_UFsx(i,j,k2)=tl_UFsx(i,j,k2)+                        &
     &                          fac1*                                   &
     &                          (tl_cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+  &
     &                           tl_cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+  &
     &                           tl_cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+  &
     &                           tl_cff4*(cff8*dnVdz-dnVdx(i,j+1,k1))+  &
     &                           cff1*(tl_cff5*dnVdz+cff5*tl_dnVdz-     &
     &                                 tl_dnVdx(i,j  ,k1))+             &
     &                           cff2*(tl_cff6*dnVdz+cff6*tl_dnVdz-     &
     &                                 tl_dnVdx(i,j+1,k2))+             &
     &                           cff3*(tl_cff7*dnVdz+cff7*tl_dnVdz-     &
     &                                 tl_dnVdx(i,j  ,k2))+             &
     &                           cff4*(tl_cff8*dnVdz+cff8*tl_dnVdz-     &
     &                                 tl_dnVdx(i,j+1,k1)))-            &
#ifdef TL_IOMS
     &                          fac1*                                   &
     &                          (cff1*(2.0_r8*cff5*dnVdz-               &
     &                                 dnVdx(i,j  ,k1))+                &
     &                           cff2*(2.0_r8*cff6*dnVdz-               &
     &                                 dnVdx(i,j+1,k2))+                &
     &                           cff3*(2.0_r8*cff7*dnVdz-               &
     &                                 dnVdx(i,j  ,k2))+                &
     &                           cff4*(2.0_r8*cff8*dnVdz-               &
     &                                 dnVdx(i,j+1,k1)))
#endif
#ifdef VISC_3DCOEF
                tl_UFsx(i,j,k2)=tl_UFsx(i,j,k2)+                        &
     &                          tl_fac1*                                &
     &                          (cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))+     &
     &                           cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))+     &
     &                           cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))+     &
     &                           cff4*(cff8*dnVdz-dnVdx(i,j+1,k1)))
#endif

                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                cff5=MIN(dZde_r(i-1,j,k1),0.0_r8)
                cff6=MIN(dZde_r(i  ,j,k2),0.0_r8)
                cff7=MAX(dZde_r(i-1,j,k2),0.0_r8)
                cff8=MAX(dZde_r(i  ,j,k1),0.0_r8)
                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i-1,j,k1)))*        &
     &                  tl_dZdx_r(i-1,j,k1)
                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i  ,j,k2)))*        &
     &                  tl_dZdx_r(i  ,j,k2)
                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_r(i-1,j,k2)))*        &
     &                  tl_dZdx_r(i-1,j,k2)
                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_r(i  ,j,k1)))*        &
     &                  tl_dZdx_r(i  ,j,k1)
                tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZde_r(i-1,j,k1)))*        &
     &                  tl_dZde_r(i-1,j,k1)
                tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZde_r(i  ,j,k2)))*        &
     &                  tl_dZde_r(i  ,j,k2)
                tl_cff7=(0.5_r8+SIGN(0.5_r8, dZde_r(i-1,j,k2)))*        &
     &                  tl_dZde_r(i-1,j,k2)
                tl_cff8=(0.5_r8+SIGN(0.5_r8, dZde_r(i  ,j,k1)))*        &
     &                  tl_dZde_r(i  ,j,k1)
!>              UFse(i,j,k2)=UFse(i,j,k2)-                              &
!>   &                       fac2*                                      &
!>   &                       (cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+        &
!>   &                        cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+        &
!>   &                        cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+        &
!>   &                        cff4*(cff8*dmVdz-dmVde(i  ,j,k1)))
!>
                tl_UFse(i,j,k2)=tl_UFse(i,j,k2)-                        &
     &                          fac2*                                   &
     &                          (tl_cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+  &
     &                           tl_cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+  &
     &                           tl_cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+  &
     &                           tl_cff4*(cff8*dmVdz-dmVde(i  ,j,k1))+  &
     &                           cff1*(tl_cff5*dmVdz+cff5*tl_dmVdz-     &
     &                                 tl_dmVde(i-1,j,k1))+             &
     &                           cff2*(tl_cff6*dmVdz+cff6*tl_dmVdz-     &
     &                                 tl_dmVde(i  ,j,k2))+             &
     &                           cff3*(tl_cff7*dmVdz+cff7*tl_dmVdz-     &
     &                                 tl_dmVde(i-1,j,k2))+             &
     &                           cff4*(tl_cff8*dmVdz+cff8*tl_dmVdz-     &
     &                                 tl_dmVde(i  ,j,k1)))+            &
#ifdef TL_IOMS
     &                          fac2*                                   &
     &                          (cff1*(2.0_r8*cff5*dmVdz-               &
     &                                 dmVde(i-1,j,k1))+                &
     &                           cff2*(2.0_r8*cff6*dmVdz-               &
     &                                 dmVde(i  ,j,k2))+                &
     &                           cff3*(2.0_r8*cff7*dmVdz-               &
     &                                 dmVde(i-1,j,k2))+                &
     &                           cff4*(2.0_r8*cff8*dmVdz-               &
     &                                 dmVde(i  ,j,k1)))
#endif
#ifdef VISC_3DCOEF
                tl_UFse(i,j,k2)=tl_UFse(i,j,k2)-                        &
     &                          tl_fac2*                                &
     &                          (cff1*(cff5*dmVdz-dmVde(i-1,j,k1))+     &
     &                           cff2*(cff6*dmVdz-dmVde(i  ,j,k2))+     &
     &                           cff3*(cff7*dmVdz-dmVde(i-1,j,k2))+     &
     &                           cff4*(cff8*dmVdz-dmVde(i  ,j,k1)))
#endif
              END DO
            END DO
!
            DO j=JstrV,Jend
              DO i=Istr,Iend
#ifdef VISC_3DCOEF
                cff=0.125_r8*                                           &
     &              (visc3d_r(i,j-1,k  )+visc3d_r(i,j,k  )+             &
     &               visc3d_r(i,j-1,k+1)+visc3d_r(i,j,k+1))
                tl_cff=0.125_r8*                                        &
     &                 (tl_visc3d_r(i,j-1,k  )+tl_visc3d_r(i,j,k  )+    &
     &                  tl_visc3d_r(i,j-1,k+1)+tl_visc3d_r(i,j,k+1))
                fac1=cff*on_v(i,j)
                tl_fac1=tl_cff*on_v(i,j)
                fac2=cff*om_v(i,j)
                tl_fac2=tl_cff*om_v(i,j)
#else
                cff=0.25_r8*(visc2_r(i,j-1)+visc2_r(i,j))
                fac1=cff*on_v(i,j)
                fac2=cff*om_v(i,j)
#endif
                cff=0.5_r8*(pn(i,j-1)+pn(i,j))
                dnUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                tl_dnUdz=cff*0.25_r8*(tl_dUdz(i  ,j  ,k2)+              &
     &                                tl_dUdz(i+1,j  ,k2)+              &
     &                                tl_dUdz(i  ,j-1,k2)+              &
     &                                tl_dUdz(i+1,j-1,k2))
                dnVdz=cff*dVdz(i,j,k2)
                tl_dnVdz=cff*tl_dVdz(i,j,k2)
                cff=0.5_r8*(pm(i,j-1)+pm(i,j))
                dmUdz=cff*0.25_r8*(dUdz(i  ,j  ,k2)+                    &
     &                             dUdz(i+1,j  ,k2)+                    &
     &                             dUdz(i  ,j-1,k2)+                    &
     &                             dUdz(i+1,j-1,k2))
                tl_dmUdz=cff*0.25_r8*(tl_dUdz(i  ,j  ,k2)+              &
     &                                tl_dUdz(i+1,j  ,k2)+              &
     &                                tl_dUdz(i  ,j-1,k2)+              &
     &                                tl_dUdz(i+1,j-1,k2))
                dmVdz=cff*dVdz(i,j,k2)
                tl_dmVdz=cff*tl_dVdz(i,j,k2)

                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*        &
     &                  tl_dZdx_p(i  ,j,k1)
                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*        &
     &                  tl_dZdx_p(i+1,j,k2)
                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*        &
     &                  tl_dZdx_p(i  ,j,k2)
                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*        &
     &                  tl_dZdx_p(i+1,j,k1)
!>              VFsx(i,j,k2)=fac1*                                      &
!>   &                       (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+        &
!>   &                        cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+        &
!>   &                        cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+        &
!>   &                        cff4*(cff4*dnVdz-dnVdx(i+1,j,k1)))
!>
                tl_VFsx(i,j,k2)=fac1*                                   &
     &                          (tl_cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+  &
     &                           tl_cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+  &
     &                           tl_cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+  &
     &                           tl_cff4*(cff4*dnVdz-dnVdx(i+1,j,k1))+  &
     &                           cff1*(tl_cff1*dnVdz+cff1*tl_dnVdz-     &
     &                                 tl_dnVdx(i  ,j,k1))+             &
     &                           cff2*(tl_cff2*dnVdz+cff2*tl_dnVdz-     &
     &                                 tl_dnVdx(i+1,j,k2))+             &
     &                           cff3*(tl_cff3*dnVdz+cff3*tl_dnVdz-     &
     &                                 tl_dnVdx(i  ,j,k2))+             &
     &                           cff4*(tl_cff4*dnVdz+cff4*tl_dnVdz-     &
     &                                 tl_dnVdx(i+1,j,k1)))-            &
#ifdef TL_IOMS
     &                          fac1*                                   &
     &                          (cff1*(2.0_r8*cff1*dnVdz-               &
     &                                 dnVdx(i  ,j,k1))+                &
     &                           cff2*(2.0_r8*cff2*dnVdz-               &
     &                                 dnVdx(i+1,j,k2))+                &
     &                           cff3*(2.0_r8*cff3*dnVdz-               &
     &                                 dnVdx(i  ,j,k2))+                &
     &                           cff4*(2.0_r8*cff4*dnVdz-               &
     &                                 dnVdx(i+1,j,k1)))
#endif
#ifdef VISC_3DCOEF
                tl_VFsx(i,j,k2)=tl_VFsx(i,j,k2)+                        &
     &                          tl_fac1*                                &
     &                          (cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))+     &
     &                           cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))+     &
     &                           cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))+     &
     &                           cff4*(cff4*dnVdz-dnVdx(i+1,j,k1)))
#endif

                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*        &
     &                  tl_dZde_r(i,j-1,k1)
                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*        &
     &                  tl_dZde_r(i,j  ,k2)
                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j-1,k2)))*        &
     &                  tl_dZde_r(i,j-1,k2)
                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j  ,k1)))*        &
     &                  tl_dZde_r(i,j  ,k1)
!>              VFse(i,j,k2)=fac2*                                      &
!>   &                       (cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+        &
!>   &                        cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+        &
!>   &                        cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+        &
!>   &                        cff4*(cff4*dmVdz-dmVde(i,j  ,k1)))
!>
                tl_VFse(i,j,k2)=fac2*                                   &
     &                          (tl_cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+  &
     &                           tl_cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+  &
     &                           tl_cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+  &
     &                           tl_cff4*(cff4*dmVdz-dmVde(i,j  ,k1))+  &
     &                           cff1*(tl_cff1*dmVdz+cff1*tl_dmVdz-     &
     &                                 tl_dmVde(i,j-1,k1))+             &
     &                           cff2*(tl_cff2*dmVdz+cff2*tl_dmVdz-     &
     &                                 tl_dmVde(i,j  ,k2))+             &
     &                           cff3*(tl_cff3*dmVdz+cff3*tl_dmVdz-     &
     &                                 tl_dmVde(i,j-1,k2))+             &
     &                           cff4*(tl_cff4*dmVdz+cff4*tl_dmVdz-     &
     &                                 tl_dmVde(i,j  ,k1)))-            &
#ifdef TL_IOMS
     &                          fac2*                                   &
     &                          (cff1*(2.0_r8*cff1*dmVdz-               &
     &                                 dmVde(i,j-1,k1))+                &
     &                           cff2*(2.0_r8*cff2*dmVdz-               &
     &                                 dmVde(i,j  ,k2))+                &
     &                           cff3*(2.0_r8*cff3*dmVdz-               &
     &                                 dmVde(i,j-1,k2))+                &
     &                           cff4*(2.0_r8*cff4*dmVdz-               &
     &                                 dmVde(i,j  ,k1)))
#endif
#ifdef VISC_3DCOEF
                tl_VFse(i,j,k2)=tl_VFse(i,j,k2)+                        &
     &                          tl_fac2*                                &
     &                          (cff1*(cff1*dmVdz-dmVde(i,j-1,k1))+     &
     &                           cff2*(cff2*dmVdz-dmVde(i,j  ,k2))+     &
     &                           cff3*(cff3*dmVdz-dmVde(i,j-1,k2))+     &
     &                           cff4*(cff4*dmVdz-dmVde(i,j  ,k1)))
#endif

                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                cff5=MIN(dZdx_r(i,j-1,k1),0.0_r8)
                cff6=MIN(dZdx_r(i,j  ,k2),0.0_r8)
                cff7=MAX(dZdx_r(i,j-1,k2),0.0_r8)
                cff8=MAX(dZdx_r(i,j  ,k1),0.0_r8)
                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j-1,k1)))*        &
     &                  tl_dZde_r(i,j-1,k1)
                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZde_r(i,j  ,k2)))*        &
     &                  tl_dZde_r(i,j  ,k2)
                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j-1,k2)))*        &
     &                  tl_dZde_r(i,j-1,k2)
                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZde_r(i,j  ,k1)))*        &
     &                  tl_dZde_r(i,j  ,k1)
                tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i,j-1,k1)))*        &
     &                  tl_dZdx_r(i,j-1,k1)
                tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZdx_r(i,j  ,k2)))*        &
     &                  tl_dZdx_r(i,j  ,k2)
                tl_cff7=(0.5_r8+SIGN(0.5_r8, dZdx_r(i,j-1,k2)))*        &
     &                  tl_dZdx_r(i,j-1,k2)
                tl_cff8=(0.5_r8+SIGN(0.5_r8, dZdx_r(i,j  ,k1)))*        &
     &                  tl_dZdx_r(i,j  ,k1)
!>              VFsx(i,j,k2)=VFsx(i,j,k2)-                              &
!>   &                       fac1*                                      &
!>   &                       (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+        &
!>   &                        cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+        &
!>   &                        cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+        &
!>   &                        cff4*(cff8*dnUdz-dnUdx(i,j  ,k1)))
!>
                tl_VFsx(i,j,k2)=tl_VFsx(i,j,k2)-                        &
     &                          fac1*                                   &
     &                          (tl_cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+  &
     &                           tl_cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+  &
     &                           tl_cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+  &
     &                           tl_cff4*(cff8*dnUdz-dnUdx(i,j  ,k1))+  &
     &                           cff1*(tl_cff5*dnUdz+cff5*tl_dnUdz-     &
     &                                 tl_dnUdx(i,j-1,k1))+             &
     &                           cff2*(tl_cff6*dnUdz+cff6*tl_dnUdz-     &
     &                                 tl_dnUdx(i,j  ,k2))+             &
     &                           cff3*(tl_cff7*dnUdz+cff7*tl_dnUdz-     &
     &                                 tl_dnUdx(i,j-1,k2))+             &
     &                           cff4*(tl_cff8*dnUdz+cff8*tl_dnUdz-     &
     &                                 tl_dnUdx(i,j  ,k1)))+            &
#ifdef TL_IOMS
     &                          fac1*                                   &
     &                          (cff1*(2.0_r8*cff5*dnUdz-               &
     &                                 dnUdx(i,j-1,k1))+                &
     &                           cff2*(2.0_r8*cff6*dnUdz-               &
     &                                 dnUdx(i,j  ,k2))+                &
     &                           cff3*(2.0_r8*cff7*dnUdz-               &
     &                                 dnUdx(i,j-1,k2))+                &
     &                           cff4*(2.0_r8*cff8*dnUdz-               &
     &                                 dnUdx(i,j  ,k1)))
#endif
#ifdef VISC_3DCOEF
                tl_VFsx(i,j,k2)=tl_VFsx(i,j,k2)-                        &
     &                          tl_fac1*                                &
     &                          (cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))+     &
     &                           cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))+     &
     &                           cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))+     &
     &                           cff4*(cff8*dnUdz-dnUdx(i,j  ,k1)))
#endif

                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                cff5=MIN(dZde_p(i  ,j,k1),0.0_r8)
                cff6=MIN(dZde_p(i+1,j,k2),0.0_r8)
                cff7=MAX(dZde_p(i  ,j,k2),0.0_r8)
                cff8=MAX(dZde_p(i+1,j,k1),0.0_r8)
                tl_cff1=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i  ,j,k1)))*        &
     &                  tl_dZdx_p(i  ,j,k1)
                tl_cff2=(0.5_r8+SIGN(0.5_r8,-dZdx_p(i+1,j,k2)))*        &
     &                  tl_dZdx_p(i+1,j,k2)
                tl_cff3=(0.5_r8+SIGN(0.5_r8, dZdx_p(i  ,j,k2)))*        &
     &                  tl_dZdx_p(i  ,j,k2)
                tl_cff4=(0.5_r8+SIGN(0.5_r8, dZdx_p(i+1,j,k1)))*        &
     &                  tl_dZdx_p(i+1,j,k1)
                tl_cff5=(0.5_r8+SIGN(0.5_r8,-dZde_p(i  ,j,k1)))*        &
     &                  tl_dZde_p(i  ,j,k1)
                tl_cff6=(0.5_r8+SIGN(0.5_r8,-dZde_p(i+1,j,k2)))*        &
     &                  tl_dZde_p(i+1,j,k2)
                tl_cff7=(0.5_r8+SIGN(0.5_r8, dZde_p(i  ,j,k2)))*        &
     &                  tl_dZde_p(i  ,j,k2)
                tl_cff8=(0.5_r8+SIGN(0.5_r8, dZde_p(i+1,j,k1)))*        &
     &                  tl_dZde_p(i+1,j,k1)
!>              VFse(i,j,k2)=VFse(i,j,k2)+                              &
!>   &                       fac2*                                      &
!>   &                       (cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+        &
!>   &                        cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+        &
!>   &                        cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+        &
!>   &                        cff4*(cff8*dmUdz-dmUde(i+1,j,k1)))
!>
                tl_VFse(i,j,k2)=tl_VFse(i,j,k2)+                        &
     &                          fac2*                                   &
     &                          (tl_cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+  &
     &                           tl_cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+  &
     &                           tl_cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+  &
     &                           tl_cff4*(cff8*dmUdz-dmUde(i+1,j,k1))+  &
     &                           cff1*(tl_cff5*dmUdz+cff5*tl_dmUdz-     &
     &                                 tl_dmUde(i  ,j,k1))+             &
     &                           cff2*(tl_cff6*dmUdz+cff6*tl_dmUdz-     &
     &                                 tl_dmUde(i+1,j,k2))+             &
     &                           cff3*(tl_cff7*dmUdz+cff7*tl_dmUdz-     &
     &                                 tl_dmUde(i  ,j,k2))+             &
     &                           cff4*(tl_cff8*dmUdz+cff8*tl_dmUdz-     &
     &                                 tl_dmUde(i+1,j,k1)))-            &
#ifdef TL_IOMS
     &                          fac2*                                   &
     &                          (cff1*(2.0_r8*cff5*dmUdz-               &
     &                                 dmUde(i  ,j,k1))+                &
     &                           cff2*(2.0_r8*cff6*dmUdz-               &
     &                                 dmUde(i+1,j,k2))+                &
     &                           cff3*(2.0_r8*cff7*dmUdz-               &
     &                                 dmUde(i  ,j,k2))+                &
     &                           cff4*(2.0_r8*cff8*dmUdz-               &
     &                                 dmUde(i+1,j,k1)))
#endif
#ifdef VISC_3DCOEF
                tl_VFse(i,j,k2)=tl_VFse(i,j,k2)+                        &
     &                          tl_fac2*                                &
     &                          (cff1*(cff5*dmUdz-dmUde(i  ,j,k1))+     &
     &                           cff2*(cff6*dmUdz-dmUde(i+1,j,k2))+     &
     &                           cff3*(cff7*dmUdz-dmUde(i  ,j,k2))+     &
     &                           cff4*(cff8*dmUdz-dmUde(i+1,j,k1)))
#endif
              END DO
            END DO
          END IF
!
!  Time-step harmonic, geopotential viscosity term.  Notice that
!  momentum at this stage is HzU and HzV and has m2/s units.  Add
!  contribution for barotropic forcing terms.
#ifdef DIAGNOSTICS_UV
!! The rotated vertical term cannot be split from the horizontal
!! terms because of the 2D/3D momentum coupling.
#endif
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=dt(ng)*0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
!>            cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*(UFx(i,j  )-UFx(i-1,j))
!>
              tl_cff1=0.5_r8*(pn(i-1,j)+pn(i,j))*                       &
     &                (tl_UFx(i,j  )-tl_UFx(i-1,j))
!>            cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i  ,j))
!>
              tl_cff2=0.5_r8*(pm(i-1,j)+pm(i,j))*                       &
     &                (tl_UFe(i,j+1)-tl_UFe(i  ,j))
!>            cff3=UFsx(i,j,k2)-UFsx(i,j,k1)
!>
              tl_cff3=tl_UFsx(i,j,k2)-tl_UFsx(i,j,k1)
!>            cff4=UFse(i,j,k2)-UFse(i,j,k1)
!>
              tl_cff4=tl_UFse(i,j,k2)-tl_UFse(i,j,k1)
!>            cff5=cff*(cff1+cff2)
!>
              tl_cff5=cff*(tl_cff1+tl_cff2)
!>            cff6=dt(ng)*(cff3+cff4)
!>
              tl_cff6=dt(ng)*(tl_cff3+tl_cff4)
!>            rufrc(i,j)=rufrc(i,j)+cff1+cff2+cff3+cff4
!>
              tl_rufrc(i,j)=tl_rufrc(i,j)+                              &
     &                      tl_cff1+tl_cff2+tl_cff3+tl_cff4
!>            u(i,j,k,nnew)=u(i,j,k,nnew)+cff5+cff6
!>
              tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+tl_cff5+tl_cff6
#ifdef DIAGNOSTICS_UV
!!            DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)+cff1+cff2+  &
!!   &                                                      cff3+cff4
!!            DiaRUfrc(i,j,3,M2xvis)=DiaRUfrc(i,j,3,M2xvis)+cff1+cff3
!!            DiaRUfrc(i,j,3,M2yvis)=DiaRUfrc(i,j,3,M2yvis)+cff2+cff4
!!            DiaU3wrk(i,j,k,M3hvis)=cff5+cff6
!!            DiaU3wrk(i,j,k,M3xvis)=cff*cff1+dt(ng)*cff3
!!            DiaU3wrk(i,j,k,M3yvis)=cff*cff2+dt(ng)*cff4
#endif
            END DO
          END DO

          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=dt(ng)*0.25_r8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
!>            cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j  ))
!>
              tl_cff1=0.5_r8*(pn(i,j-1)+pn(i,j))*                       &
     &                (tl_VFx(i+1,j)-tl_VFx(i,j  ))
!>            cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*(VFe(i  ,j)-VFe(i,j-1))
!>
              tl_cff2=0.5_r8*(pm(i,j-1)+pm(i,j))*                       &
     &                (tl_VFe(i  ,j)-tl_VFe(i,j-1))
!>            cff3=VFsx(i,j,k2)-VFsx(i,j,k1)
!>
              tl_cff3=tl_VFsx(i,j,k2)-tl_VFsx(i,j,k1)
!>            cff4=VFse(i,j,k2)-VFse(i,j,k1)
!>
              tl_cff4=tl_VFse(i,j,k2)-tl_VFse(i,j,k1)
!>            cff5=cff*(cff1-cff2)
!>
              tl_cff5=cff*(tl_cff1-tl_cff2)
!>            cff6=dt(ng)*(cff3+cff4)
!>
              tl_cff6=dt(ng)*(tl_cff3+tl_cff4)
!>            rvfrc(i,j)=rvfrc(i,j)+cff1-cff2+cff3+cff4
!>
              tl_rvfrc(i,j)=tl_rvfrc(i,j)+                              &
     &                      tl_cff1-tl_cff2+tl_cff3+tl_cff4
!>            v(i,j,k,nnew)=v(i,j,k,nnew)+cff5+cff6
!>
              tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+tl_cff5+tl_cff6
#ifdef DIAGNOSTICS_UV
!!            DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)+cff1-cff2+  &
!!   &                                                      cff3+cff4
!!            DiaRVfrc(i,j,3,M2xvis)=DiaRVfrc(i,j,3,M2xvis)+cff1+cff3
!!            DiaRVfrc(i,j,3,M2yvis)=DiaRVfrc(i,j,3,M2yvis)-cff2+cff4
!!            DiaV3wrk(i,j,k,M3hvis)=cff5+cff6
!!            DiaV3wrk(i,j,k,M3xvis)= cff*cff1+dt(ng)*cff3
!!            DiaV3wrk(i,j,k,M3yvis)=-cff*cff2+dt(ng)*cff4
#endif
            END DO
          END DO
        END IF
      END DO K_LOOP
      RETURN
      END SUBROUTINE rp_uv3dmix2_tile
